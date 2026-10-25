# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Source and wheel provenance producer/consumer regression tests."""
import hashlib
import base64
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import textwrap
import types
import unittest
from unittest import mock
import zipfile

SCRIPTS = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD = load("manifest_build", SCRIPTS / "build_python_wheels.py")
COLLECTOR = load("manifest_collector", SCRIPTS / "collect-release-wheels.py")
GATE = load("manifest_gate", SCRIPTS / "check_pypi_release.py")
WHEELS = load("manifest_wheel_fixtures", Path(__file__).with_name("test_collect_release_wheels.py"))
SHA = "1234567890abcdef" * 2 + "12345678"


def archive(directory, files):
    path = directory / "smsd-7.2.0.tar.gz"
    with tarfile.open(path, "w:gz") as output:
        for name, data in files.items():
            item = tarfile.TarInfo("smsd-7.2.0/" + name)
            item.size = len(data)
            output.addfile(item, io.BytesIO(data))
    return path


class ReleaseManifestsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.directory = self.root / "release"
        self.directory.mkdir()
        self.checkout = self.root / "checkout"
        self.sources = {**WHEELS.SOURCES,
            "pyproject.toml": b'[project]\nname = "smsd"\nversion = "7.2.0"\nrequires-python = ">=3.10"\n[tool.scikit-build]\nsdist.exclude = ["packaging/**"]\n',
            "pom.xml": b'<project xmlns="http://maven.apache.org/POM/4.0.0"><version>7.2.0</version></project>',
            "java/pom.xml": b'<project xmlns="http://maven.apache.org/POM/4.0.0"><version>7.2.0</version></project>',
            "cpp/CMakeLists.txt": b'project(smsd VERSION 7.2.0 LANGUAGES CXX)\n',
            "python/README.md": b'Fixture package description\n',
            "PKG-INFO": b'Name: smsd\nVersion: 7.2.0\nRequires-Python: >=3.10\n\nFixture package description\n',
        }
        for name in ("cpp/include/smsd/mcs.hpp", "cpp/bindings/pybind11/smsd_bindings.cpp",
                     "python/tests/test_smsd.py", "python/tests/test_aromaticity_parity.py",
                     "java/src/test/resources/com/bioinception/smsd/aromaticity_parity.json",
                     "scripts/check_python_wheel.py", "scripts/check_native_wheel_build.py",
                     "scripts/build_python_wheels.py", "scripts/collect-release-wheels.py", "scripts/cibuildwheel.toml",
                     "scripts/repair_windows_wheel.py", "licenses/msvc/README.md",
                     "licenses/msvc/LICENSE-2022.docx", "licenses/msvc/LICENSE-2026.docx",
                     "Dockerfile", "Dockerfile.cuda", ".dockerignore", "scripts/docker_install.py",
                     "scripts/docker_smoke.py", "scripts/prepare_docker_context.py",
                     "java/core/pom.xml", "java/cli/pom.xml",
                     "java/src/main/java/com/bioinception/smsd/core/MCS.java",
                     "java/src/main/java/com/bioinception/smsd/core/Substructure.java",
                     "java/cli/src/main/java/com/bioinception/smsd/cli/SMSDcli.java",
                     "licenses/cdk/LICENSE", "licenses/beam/LICENSE", "licenses/vecmath/LICENSE"):
            self.sources[name] = b'fixture\n'
        for name, data in self.sources.items():
            path = self.checkout / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        self.sdist = archive(self.directory, self.sources)
        self.source = {**BUILD.inspect_source(self.sdist), "schema_version": 1,
                       "source_commit": SHA, "release_ref": SHA, "working_tree_clean": True}
        self.manifest = self.directory / "SOURCE_MANIFEST.json"
        self.manifest.write_text(json.dumps(self.source))
        self.proofs = []
        self.cuda_sources = {}
        self.targets = {}
        self.attributes = {}
        self.info_attributes = self.root / "info-attributes"

    def git(self, command, **kwargs):
        arguments = command[1:]
        if "status" in arguments:
            return ""
        if "--git-path" in arguments:
            return str(self.info_attributes) + "\n"
        if "rev-parse" in arguments:
            return SHA + "\n"
        if "check-attr" in arguments:
            self.assertIn("--source=" + SHA, arguments)
            self.assertIn("core.attributesFile=" + os.devnull, arguments)
            self.assertEqual(kwargs["env"]["GIT_ATTR_NOSYSTEM"], "1")
            values = []
            for name in kwargs["input"].decode().split("\0"):
                if name:
                    for attribute in ("text", "eol", "filter"):
                        values.extend((name, attribute, self.attributes.get(name, {}).get(attribute, "unspecified")))
            return ("\0".join(values) + "\0").encode()
        if "ls-tree" in arguments:
            entries = []
            for name, data in {**self.sources, **self.cuda_sources}.items():
                if name == "PKG-INFO":
                    continue
                blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
                entries.append(f"100644 blob {blob}\t{name}".encode())
            return b"\0".join(entries) + b"\0"
        raise AssertionError(command)

    def wheel(self, family="macos", python_version="3.14", architecture=None, **proof_changes):
        architecture = architecture or BUILD.ARCHITECTURES[family][0]
        tag = {"macos": "macosx_arm64", "linux": "manylinux_" + architecture, "windows": "win_amd64"}[family]
        identifier = "cp" + python_version.replace(".", "") + "-" + tag
        directory = self.root / identifier
        directory.mkdir(exist_ok=True)
        wheel = WHEELS.make_wheel(directory, family=family, repaired=False, sources=self.sources,
                                 python_version=python_version, architecture=architecture,
                                 changes={WHEELS.PREFIX + "METADATA": self.sources["PKG-INFO"],
                                          "include/smsd/mcs.hpp": self.sources["cpp/include/smsd/mcs.hpp"]})
        proof = {**self.source, **BUILD.BUILD_IDENTIFIERS[identifier], "wheel": wheel.name, "wheel_sha256": BUILD.sha256(wheel),
                 "build_identifier": identifier, "build_system": "cibuildwheel", "cibuildwheel": "4.2.1",
                 "native_debug_tests": "passed", "installed_wheel_tests": "passed"}
        proof.update(proof_changes)
        path = directory / "wheel-build-provenance.json"
        path.write_text(json.dumps(proof))
        self.proofs.append(path)
        self.targets[identifier] = wheel.name
        return wheel, path

    def all_wheels(self):
        for target in BUILD.BUILD_IDENTIFIERS.values():
            self.wheel(target["platform"], target["python_version"], target["architecture"])

    def cuda_wheel(self, changes=None):
        root = "packaging/smsd-cuda/"
        files = {"pyproject.toml": b'[project]\nname="smsd-cuda"\nversion="7.2.0"\nrequires-python=">=3.10"\ndependencies=["smsd==7.2.0", "nvidia-cuda-runtime-cu12>=12.8,<13"]\n',
                 "README.md": b"CUDA package description\n", "LICENSE": b"CUDA licence\n", "NOTICE": b"CUDA notice\n",
                 "python/smsd_cuda/__init__.py": b'__version__ = "7.2.0"\n__smsd_version__ = "7.2.0"\n'}
        self.cuda_sources = {root + name: data for name, data in files.items()}
        for name, data in files.items():
            path = self.checkout / root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        prefix = "smsd_cuda-7.2.0.dist-info/"
        binary = bytearray(64); binary[:6] = b"\x7fELF\x02\x01"; binary[18:20] = b"\x3e\x00"
        entries = {"smsd_cuda/__init__.py": files["python/smsd_cuda/__init__.py"],
                   "smsd_cuda/libsmsd_cuda_screen.so": bytes(binary),
                   prefix + "METADATA": b"Metadata-Version: 2.4\nName: smsd-cuda\nVersion: 7.2.0\nRequires-Python: >=3.10\nRequires-Dist: smsd==7.2.0\nRequires-Dist: nvidia-cuda-runtime-cu12<13,>=12.8\n\n" + files["README.md"],
                   prefix + "WHEEL": b"Wheel-Version: 1.0\nRoot-Is-Purelib: false\nTag: py3-none-manylinux_2_34_x86_64\n",
                   prefix + "licenses/LICENSE": files["LICENSE"], prefix + "licenses/NOTICE": files["NOTICE"]}
        entries.update(changes or {})
        record = io.StringIO()
        writer = csv.writer(record, lineterminator="\n")
        for name, data in entries.items():
            if name != prefix + "RECORD":
                digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
                writer.writerow([name, "sha256=" + digest, str(len(data))])
        writer.writerow([prefix + "RECORD", "", ""])
        entries.setdefault(prefix + "RECORD", record.getvalue().encode())
        wheel = self.directory / GATE.cuda_wheel_name("7.2.0")
        with zipfile.ZipFile(wheel, "w") as output:
            for name, data in entries.items():
                output.writestr(name, data)
        proof = {"schema_version": 1, "source_commit": SHA, "sdist_sha256": self.source["sdist_sha256"],
                 "smsd_version": "7.2.0", "wheel": wheel.name, "wheel_sha256": BUILD.sha256(wheel),
                 "installed_wheel_tests": "passed", "native_gpu_tests": "passed"}
        (self.directory / "CUDA_WHEEL_MANIFEST.json").write_text(json.dumps(proof))
        return wheel, proof

    def collect(self, *extra):
        arguments = ["collect-release-wheels.py", "--release-dir", str(self.directory)]
        for identifier in self.targets:
            arguments += ["--wheel-dir", str(self.root / identifier)]
        arguments += list(extra)
        with mock.patch.object(COLLECTOR, "__file__", str(self.checkout / "scripts/collect-release-wheels.py")), mock.patch("sys.argv", arguments):
            COLLECTOR.main()

    def test_prepare_resolves_tag_and_emits_consumer_filename(self):
        output = self.root / "prepared"
        args = types.SimpleNamespace(output_dir=output, release_ref="v7.2.0", no_isolation=True)
        def build_source(command, **kwargs):
            self.assertTrue(kwargs["check"])
            archive(output, self.sources)
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git), mock.patch.object(BUILD.subprocess, "run", side_effect=build_source):
            BUILD.prepare(args)
        result = json.loads((output / "SOURCE_MANIFEST.json").read_text())
        self.assertEqual(result["release_ref"], SHA)
        self.assertEqual(result["source_commit"], SHA)
        self.assertTrue(result["working_tree_clean"])
        self.assertEqual(result["source_file_count"], len(self.sources))
        GATE.check_source(output, "7.2.0", result)

    def test_prepare_requires_clean_frozen_identity(self):
        for release_ref, dirty in ((None, False), ("v7.2.0", True)):
            with self.subTest(release_ref=release_ref, dirty=dirty):
                args = types.SimpleNamespace(output_dir=self.root / "dirty", release_ref=release_ref, no_isolation=True)
                def git(command, **kwargs):
                    return " M cpp/file.hpp\n" if dirty and "status" in command else self.git(command, **kwargs)
                with mock.patch.object(BUILD.subprocess, "check_output", side_effect=git), mock.patch.object(BUILD.subprocess, "run") as run:
                    with self.assertRaises(ValueError):
                        BUILD.prepare(args)
                    run.assert_not_called()

    def test_check_source_requires_the_canonical_archive_and_clean_frozen_checkout(self):
        args = types.SimpleNamespace(sdist=self.sdist, source_manifest=self.manifest, release_commit=SHA)
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git), mock.patch.object(BUILD.subprocess, "run") as run:
            BUILD.check_source(args)
        run.assert_not_called()
        for commit in ("HEAD", SHA[:12], "v7.2.0", "0" * 40):
            with self.subTest(commit=commit):
                args.release_commit = commit
                with mock.patch.object(BUILD.subprocess, "check_output") as git:
                    with self.assertRaises(ValueError):
                        BUILD.check_source(args)
                    git.assert_not_called()

    def test_check_source_rejects_checkout_changes_before_or_during_archive_check(self):
        args = types.SimpleNamespace(sdist=self.sdist, source_manifest=self.manifest, release_commit=SHA)
        for change, after in (("dirty", False), ("head", False), ("dirty", True), ("head", True)):
            with self.subTest(change=change, after=after):
                archive_checked = False
                def verify(*arguments):
                    nonlocal archive_checked
                    archive_checked = True
                def git(command, **kwargs):
                    changed = not after or archive_checked
                    if changed and change == "dirty" and "status" in command:
                        return " M cpp/file.hpp\n"
                    if changed and change == "head" and "rev-parse" in command:
                        return "0" * 40 + "\n"
                    return self.git(command, **kwargs)
                with mock.patch.object(BUILD.subprocess, "check_output", side_effect=git), mock.patch.object(BUILD, "verify_source_tree", side_effect=verify):
                    with self.assertRaisesRegex(ValueError, "clean checkout"):
                        BUILD.check_source(args)

    def test_check_source_rejects_archive_files_outside_the_frozen_tree(self):
        args = types.SimpleNamespace(sdist=self.sdist, source_manifest=self.manifest, release_commit=SHA)
        archive(self.directory, {**self.sources, "python/smsd/search.py": b"changed\n"})
        self.manifest.write_text(json.dumps({**self.source, **BUILD.inspect_source(self.sdist)}))
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
            with self.assertRaisesRegex(ValueError, "frozen"):
                BUILD.check_source(args)

    def test_prepare_rejects_sdist_file_not_in_frozen_tree(self):
        output = self.root / "different"
        args = types.SimpleNamespace(output_dir=output, release_ref="v7.2.0", no_isolation=True)
        def changed_source(command, **kwargs):
            archive(output, {**self.sources, "python/smsd/search.py": b'changed after freeze\n'})
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git), mock.patch.object(BUILD.subprocess, "run", side_effect=changed_source):
            with self.assertRaises(ValueError):
                BUILD.prepare(args)
        self.assertFalse((output / "SOURCE_MANIFEST.json").exists())

    def test_source_tree_accepts_committed_crlf_checkout(self):
        filename = "java/cli/src/scripts/smsd.bat"
        self.sources[".gitattributes"] = b"*.bat text eol=crlf\n"
        self.sources[filename] = b"@echo off\necho ready\n"
        self.attributes[filename] = {"text": "set", "eol": "crlf"}
        source = archive(self.directory, {**self.sources, filename: b"@echo off\r\necho ready\r\n"})
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
            BUILD.verify_source_tree(source, SHA)

    def test_source_tree_rejects_crlf_without_committed_text_policy(self):
        filename = "java/cli/src/scripts/smsd.bat"
        self.sources[filename] = b"@echo off\necho ready\n"
        for rule, attributes in ((b"", {}), (b"*.bat binary\n", {"text": "unset"}),
                                 (b"*.bat text eol=lf\n", {"text": "set", "eol": "lf"}),
                                 (b"*.bat text eol=crlf filter=custom\n", {"text": "set", "eol": "crlf", "filter": "custom"})):
            with self.subTest(rule=rule):
                self.sources[".gitattributes"] = rule
                self.attributes[filename] = attributes
                source = archive(self.directory, {**self.sources, filename: b"@echo off\r\necho ready\r\n"})
                with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
                    with self.assertRaises(ValueError):
                        BUILD.verify_source_tree(source, SHA)

    def test_source_tree_rejects_changed_or_mixed_crlf_checkout(self):
        filename = "java/cli/src/scripts/smsd.bat"
        self.sources[".gitattributes"] = b"*.bat text eol=crlf\n"
        self.sources[filename] = b"@echo off\necho ready\n"
        self.attributes[filename] = {"text": "set", "eol": "crlf"}
        for data in (b"@echo off\r\necho changed\r\n", b"@echo off\r\necho ready\n"):
            with self.subTest(data=data):
                source = archive(self.directory, {**self.sources, filename: data})
                with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
                    with self.assertRaises(ValueError):
                        BUILD.verify_source_tree(source, SHA)

    def test_source_tree_rejects_local_attribute_override_for_crlf(self):
        filename = "java/cli/src/scripts/smsd.bat"
        self.sources[".gitattributes"] = b"*.bat text eol=crlf\n"
        self.sources[filename] = b"@echo off\necho ready\n"
        self.attributes[filename] = {"text": "set", "eol": "crlf"}
        self.info_attributes.write_text("*.bat text eol=crlf\n")
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
            BUILD.verify_source_tree(archive(self.directory, self.sources), SHA)
            source = archive(self.directory, {**self.sources, filename: b"@echo off\r\necho ready\r\n"})
            with self.assertRaises(ValueError):
                BUILD.verify_source_tree(source, SHA)

    def test_source_distribution_cannot_drop_another_frozen_java_build_input(self):
        filename = "java/src/main/java/com/bioinception/smsd/core/Other.java"
        self.sources[filename] = b"class Other {}\n"
        output = self.root / "missing-java"
        args = types.SimpleNamespace(output_dir=output, release_ref="v7.2.0", no_isolation=True)
        def incomplete_source(command, **kwargs):
            archive(output, {name: value for name, value in self.sources.items() if name != filename})
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git), mock.patch.object(BUILD.subprocess, "run", side_effect=incomplete_source):
            with self.assertRaisesRegex(ValueError, "omits frozen Java/Docker inputs"):
                BUILD.prepare(args)
        self.assertFalse((output / "SOURCE_MANIFEST.json").exists())

    def test_collector_rejects_sdk_missing_the_checkouts_c_abi_header(self):
        del self.sources["cpp/include/smsd/screen_accelerator_abi.h"]
        self.sdist = archive(self.directory, self.sources)
        self.source.update(BUILD.inspect_source(self.sdist))
        self.manifest.write_text(json.dumps(self.source))
        self.wheel("macos")
        with self.assertRaisesRegex(ValueError, "Source distribution and release checkout differ"):
            self.collect("--allow-incomplete")

    def test_collector_rejects_missing_provenance(self):
        _, proof = self.wheel()
        proof.unlink()
        with self.assertRaises(ValueError):
            self.collect("--allow-incomplete")
        self.assertFalse((self.directory / "WHEEL_MANIFEST.json").exists())

    def test_collector_emits_all_platforms_from_valid_receipts(self):
        self.all_wheels()
        self.collect()
        result = json.loads((self.directory / "WHEEL_MANIFEST.json").read_text())
        self.assertTrue(result["complete"])
        self.assertEqual(len(result["wheels"]), 20)
        self.assertEqual({entry["build_identifier"] for entry in result["wheels"]}, set(self.targets))
        self.assertEqual({entry["python_version"] for entry in result["wheels"]}, {"3.10", "3.11", "3.12", "3.13", "3.14"})
        GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.targets, self.sources["PKG-INFO"])
        self.assertNotIn(str(self.root), json.dumps(result))

    def test_collector_rejects_mismatched_or_unpassed_receipts(self):
        wheel, proof = self.wheel()
        original = json.loads(proof.read_text())
        for field, value in (("wheel_sha256", "0" * 64), ("source_commit", "0" * 40),
                             ("sdist_sha256", "0" * 64), ("release_ref", "v7.2.0"),
                             ("native_debug_tests", "not run"), ("installed_wheel_tests", None),
                             ("build_identifier", "cp314-win_amd64"), ("python_version", "3.10"),
                             ("architecture", "aarch64"), ("platform", "linux")):
            with self.subTest(field=field):
                proof.write_text(json.dumps({**original, field: value}))
                with self.assertRaises(ValueError):
                    self.collect("--allow-incomplete")
        self.assertFalse((self.directory / wheel.name).exists())

    def test_incomplete_collection_never_invents_platform_test_results(self):
        self.wheel("macos")
        self.collect("--allow-incomplete")
        result = json.loads((self.directory / "WHEEL_MANIFEST.json").read_text())
        self.assertFalse(result["complete"])
        self.assertEqual(result["missing_targets"], sorted(set(BUILD.BUILD_IDENTIFIERS) - {"cp314-macosx_arm64"}))
        self.assertEqual([entry["platform"] for entry in result["wheels"]], ["macos"])

    def test_check_only_accepts_published_manifest_without_sidecars(self):
        self.all_wheels()
        self.collect()
        for proof in self.proofs:
            proof.unlink()
        before = {path.name: path.read_bytes() for path in self.directory.iterdir()}
        self.targets.clear()
        self.collect("--check-only")
        self.assertEqual(before, {path.name: path.read_bytes() for path in self.directory.iterdir()})

    def test_check_only_rejects_hand_assembled_passed_strings(self):
        self.all_wheels()
        self.collect()
        path = self.directory / "WHEEL_MANIFEST.json"
        manifest = json.loads(path.read_text())
        for entry in manifest["wheels"]:
            entry.pop("build_provenance", None)
        path.write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            self.collect("--check-only")
        with self.assertRaises((ValueError, SystemExit)):
            GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.targets, self.sources["PKG-INFO"])

    def test_publish_gate_rejects_status_only_test_evidence(self):
        entries = []
        for identifier, target in BUILD.BUILD_IDENTIFIERS.items():
            wheel, _ = self.wheel(target["platform"], target["python_version"], target["architecture"])
            (self.directory / wheel.name).write_bytes(wheel.read_bytes())
            entries.append({**target, "build_identifier": identifier, "wheel": wheel.name,
                            "wheel_sha256": BUILD.sha256(wheel), "native_debug_tests": "passed",
                            "installed_wheel_tests": "passed"})
        manifest = {"schema_version": 2, "complete": True, "missing_targets": [], "wheels": entries,
                    **{field: self.source[field] for field in ("version", "source_commit", "sdist", "sdist_sha256")}}
        (self.directory / "WHEEL_MANIFEST.json").write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.targets, self.sources["PKG-INFO"])

    def test_one_missing_interpreter_is_incomplete_even_with_every_platform(self):
        missing = "cp310-manylinux_aarch64"
        for identifier, target in BUILD.BUILD_IDENTIFIERS.items():
            if identifier != missing:
                self.wheel(target["platform"], target["python_version"], target["architecture"])
        with self.assertRaisesRegex(ValueError, "Missing release wheels"):
            self.collect()
        self.collect("--allow-incomplete")
        manifest = json.loads((self.directory / "WHEEL_MANIFEST.json").read_text())
        self.assertFalse(manifest["complete"])
        self.assertEqual(manifest["missing_targets"], [missing])

    def published_metadata(self):
        return {"tag_name": "v7.2.0", "draft": False, "prerelease": False, "published_at": "2026-10-23T00:00:00Z",
                "assets": [{"name": path.name, "state": "uploaded", "size": path.stat().st_size,
                            "digest": "sha256:" + BUILD.sha256(path)} for path in self.directory.iterdir()]}

    def test_publish_selects_twenty_cpu_wheels_without_uploading_the_separate_cuda_package(self):
        self.all_wheels()
        self.collect()
        metadata = self.published_metadata()
        cuda = "smsd_cuda-7.2.0-py3-none-manylinux_2_34_x86_64.whl"
        metadata["assets"].append({"name": cuda, "state": "uploaded", "size": 7, "digest": "sha256:" + "a" * 64})
        version, selected, targets = GATE.selected_assets(metadata, "v7.2.0")
        self.assertEqual(version, "7.2.0")
        self.assertEqual(set(targets), set(BUILD.BUILD_IDENTIFIERS))
        self.assertEqual(len(selected), 24)
        self.assertNotIn(cuda, selected)
        metadata["assets"] = [row for row in metadata["assets"] if row["name"] != targets["cp310-win_amd64"]]
        with self.assertRaisesRegex(ValueError, "Missing release wheels"):
            GATE.selected_assets(metadata, "v7.2.0")

    def test_cuda_publish_selects_only_the_declared_wheel_and_source_evidence(self):
        wheel, _ = self.cuda_wheel()
        (self.directory / "SHA256SUMS").write_text("fixture\n")
        metadata = self.published_metadata()
        _, selected, targets = GATE.selected_assets(metadata, "v7.2.0", "smsd-cuda")
        self.assertEqual(targets, {"cuda": wheel.name})
        self.assertEqual(set(selected), {wheel.name, self.sdist.name, "SOURCE_MANIFEST.json", "CUDA_WHEEL_MANIFEST.json", "SHA256SUMS"})
        with self.assertRaisesRegex(ValueError, "Missing release wheels"):
            GATE.selected_assets(metadata, "v7.2.0")
        for name in (wheel.name, "CUDA_WHEEL_MANIFEST.json", "SOURCE_MANIFEST.json"):
            with self.subTest(missing=name), self.assertRaises(ValueError):
                GATE.selected_assets({**metadata, "assets": [row for row in metadata["assets"] if row["name"] != name]}, "v7.2.0", "smsd-cuda")
        for extra in (wheel.name.replace("2_34", "2_35"), wheel.name.replace("x86_64", "aarch64"), "smsd_cuda-7.2.0.tar.gz"):
            with self.subTest(extra=extra), self.assertRaises(ValueError):
                GATE.selected_assets({**metadata, "assets": metadata["assets"] + [{"name": extra}]}, "v7.2.0", "smsd-cuda")

    def test_cuda_proof_requires_exact_source_wheel_and_successful_actual_checks(self):
        _, proof = self.cuda_wheel()
        path = self.directory / "CUDA_WHEEL_MANIFEST.json"
        GATE.check_cuda_wheel(self.directory, "7.2.0", self.source, self.cuda_sources)
        mutations = [{key: value for key, value in proof.items() if key != missing} for missing in proof]
        mutations += [{**proof, key: value} for key, value in (("schema_version", True), ("source_commit", "0" * 40),
                      ("sdist_sha256", "0" * 64), ("wheel_sha256", "0" * 64), ("smsd_version", "7.2.1"),
                      ("wheel", "another.whl"), ("installed_wheel_tests", False), ("native_gpu_tests", "failed"), ("unexpected", "value"))]
        for index, changed in enumerate(mutations):
            with self.subTest(index=index):
                path.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, "proof"):
                    GATE.check_cuda_wheel(self.directory, "7.2.0", self.source, self.cuda_sources)

    def test_cuda_wheel_rejects_changed_metadata_source_binary_and_record(self):
        prefix = "smsd_cuda-7.2.0.dist-info/"
        changes = [{"smsd_cuda/__init__.py": b"changed\n"}, {"smsd_cuda/libsmsd_cuda_screen.so": b"not an ELF library"},
                   {prefix + "METADATA": b"Name: smsd-cuda\nVersion: 7.2.1\nRequires-Python: >=3.10\n\nwrong description\n"},
                   {prefix + "WHEEL": b"Root-Is-Purelib: true\nTag: py3-none-any\n"},
                   {prefix + "licenses/LICENSE": b"changed licence"}, {prefix + "RECORD": b"invalid record\n"},
                   {"../unexpected.py": b"unsafe\n"}]
        for change in changes:
            with self.subTest(change=list(change)):
                self.cuda_wheel(change)
                with self.assertRaises(ValueError):
                    GATE.check_cuda_wheel(self.directory, "7.2.0", self.source, self.cuda_sources)

    def test_cuda_package_inputs_come_from_the_frozen_checkout_not_the_core_sdk(self):
        source_hash = self.source["sdist_sha256"]
        self.cuda_wheel()
        files = GATE.check_source(self.directory, "7.2.0", self.source)
        self.assertFalse(any(name.startswith("packaging/") for name in files))
        self.assertEqual(BUILD.sha256(self.sdist), source_hash)
        project = GATE.tomllib.loads((SCRIPTS.parent / "pyproject.toml").read_text())
        self.assertIn("packaging/**", project["tool"]["scikit-build"]["sdist"]["exclude"])
        listing = self.git(["git", "ls-tree"])
        with mock.patch.object(GATE, "git", return_value=listing) as git:
            cuda_files = GATE.check_cuda_source(self.checkout, self.source)
        git.assert_called_once_with(self.checkout, "ls-tree", "-r", "-z", "--full-tree", SHA)
        self.assertEqual(cuda_files, self.cuda_sources)
        GATE.check_cuda_wheel(self.directory, "7.2.0", self.source, cuda_files)

    def test_cuda_checkout_rejects_stale_missing_and_untracked_package_inputs(self):
        self.cuda_wheel()
        listing = self.git(["git", "ls-tree"])
        for name, original in self.cuda_sources.items():
            with self.subTest(stale=name), mock.patch.object(GATE, "git", return_value=listing):
                path = self.checkout / name
                path.write_bytes(original + b"stale\n")
                with self.assertRaisesRegex(ValueError, "differs from frozen source"):
                    GATE.check_cuda_source(self.checkout, self.source)
                path.write_bytes(original)
        name = "packaging/smsd-cuda/NOTICE"
        path = self.checkout / name
        path.unlink()
        with mock.patch.object(GATE, "git", return_value=listing):
            with self.assertRaisesRegex(ValueError, "Missing or unsafe"):
                GATE.check_cuda_source(self.checkout, self.source)
        path.write_bytes(self.cuda_sources[name])
        entries = [entry for entry in listing.split(b"\0") if entry and not entry.endswith(b"\t" + name.encode())]
        with mock.patch.object(GATE, "git", return_value=b"\0".join(entries) + b"\0"):
            with self.assertRaisesRegex(ValueError, "Frozen commit is missing"):
                GATE.check_cuda_source(self.checkout, self.source)

    def test_cuda_validation_stages_one_unchanged_wheel_and_rejects_stage_mutations(self):
        wheel, _ = self.cuda_wheel()
        (self.directory / "SHA256SUMS").write_text("".join(BUILD.sha256(path) + "  " + path.name + "\n" for path in self.directory.iterdir()))
        metadata = self.root / "release-metadata.json"
        metadata.write_text(json.dumps(self.published_metadata()))
        args = types.SimpleNamespace(tag="v7.2.0", package="smsd-cuda", release_dir=self.directory, metadata_file=metadata,
                                     source_checkout=self.checkout, dist_dir=self.root / "upload", receipt_file=self.root / "receipt.json")
        def git(checkout, *arguments):
            if arguments[0] == "rev-parse":
                return SHA.encode()
            if arguments[0] == "diff":
                return b""
            if arguments[0] == "ls-tree":
                return self.git(["git", "ls-tree"])
            raise AssertionError(arguments)
        with mock.patch.object(GATE, "git", side_effect=git), mock.patch.object(GATE.subprocess, "run", return_value=types.SimpleNamespace(returncode=0)) as run:
            GATE.validate(args)
        self.assertEqual(len(run.call_args_list), 1)
        self.assertIn("merge-base", run.call_args.args[0])
        self.assertEqual([path.name for path in args.dist_dir.iterdir()], [wheel.name])
        self.assertEqual((args.dist_dir / wheel.name).read_bytes(), wheel.read_bytes())
        GATE.verify_staged(args)
        receipt = json.loads(args.receipt_file.read_text())
        for missing in receipt:
            with self.subTest(missing=missing):
                args.receipt_file.write_text(json.dumps({key: value for key, value in receipt.items() if key != missing}))
                with self.assertRaises((ValueError, KeyError)):
                    GATE.verify_staged(args)
        args.receipt_file.write_text(json.dumps(receipt))
        args.package = "smsd"
        with self.assertRaisesRegex(ValueError, "package"):
            GATE.verify_staged(args)
        args.package = "smsd-cuda"
        (args.dist_dir / "unexpected.whl").write_bytes(b"unexpected")
        with self.assertRaisesRegex(ValueError, "contents"):
            GATE.verify_staged(args)
        (args.dist_dir / "unexpected.whl").unlink()
        (args.dist_dir / wheel.name).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "hash"):
            GATE.verify_staged(args)

    def test_cuda_source_description_and_paired_version_must_match(self):
        self.cuda_wheel()
        changed = dict(self.cuda_sources)
        changed["packaging/smsd-cuda/python/smsd_cuda/__init__.py"] = b'__version__ = "7.2.1"\n__smsd_version__ = "7.2.0"\n'
        with self.assertRaisesRegex(ValueError, "module version"):
            GATE.check_cuda_wheel(self.directory, "7.2.0", self.source, changed)
        with mock.patch.object(GATE, "git", side_effect=lambda checkout, *args: b"packaging/smsd-cuda/README.md\0" if args[0] == "diff" else SHA.encode()), mock.patch.object(GATE.subprocess, "run", return_value=types.SimpleNamespace(returncode=0)):
            with self.assertRaisesRegex(ValueError, "CUDA package description"):
                GATE.check_git(self.checkout, "v7.2.0", self.source, self.sources, "smsd-cuda")

    def test_cuda_download_requires_linux_before_any_service_call(self):
        for package, platform in (("smsd-cuda", "windows"), ("smsd-cuda", "all"), ("smsd", "linux")):
            with self.subTest(package=package, platform=platform), mock.patch.object(GATE.subprocess, "check_output") as service:
                with self.assertRaisesRegex(ValueError, "platform"):
                    GATE.download(types.SimpleNamespace(package=package, platform=platform))
                service.assert_not_called()

    def test_cuda_selection_cannot_start_a_core_build_matrix(self):
        workflow = (SCRIPTS.parent / ".github/workflows/python-publish.yml").read_text()
        snippet = workflow.split("          python - <<'PY'\n", 1)[1].split("          PY\n", 1)[0]
        with mock.patch.dict(os.environ, {"DISTRIBUTION": "smsd-cuda"}), mock.patch.object(subprocess, "check_output") as git:
            with self.assertRaisesRegex(SystemExit, "only publishes an existing"):
                exec(compile(textwrap.dedent(snippet), "workflow-cuda-selection", "exec"), {})
            git.assert_not_called()

    def test_staged_upload_requires_the_twenty_targets_and_original_sdist(self):
        self.all_wheels()
        self.collect()
        directory = self.root / "upload"
        directory.mkdir()
        names = [self.sdist.name, *self.targets.values()]
        for name in names:
            shutil.copyfile(self.directory / name, directory / name)
        receipt = self.root / "upload-receipt.json"
        receipt.write_text(json.dumps({"version": "7.2.0", "source_commit": SHA,
                                      "files": {name: BUILD.sha256(directory / name) for name in names}}))
        args = types.SimpleNamespace(receipt_file=receipt, dist_dir=directory)
        GATE.verify_staged(args)
        (directory / self.targets["cp310-manylinux_aarch64"]).unlink()
        with self.assertRaises(ValueError):
            GATE.verify_staged(args)

    def check_tagged_documentation(self, changed):
        tagged = "f" * 40
        def git(checkout, *arguments):
            if arguments[0] == "rev-parse":
                return (SHA if arguments[-1] == SHA + "^{commit}" else tagged).encode()
            if arguments[0] == "diff":
                return changed.encode() + b"\0"
            if arguments[0] == "ls-tree":
                return self.git(["git", "ls-tree"])
            raise AssertionError(arguments)
        with mock.patch.object(GATE, "git", side_effect=git), mock.patch.object(GATE.subprocess, "run", return_value=types.SimpleNamespace(returncode=0)):
            return GATE.check_git(self.checkout, "v7.2.0", self.source, self.sources)

    def test_publish_gate_rejects_changed_package_description_at_tag(self):
        for content in (b"Changed package description\n", self.sources["python/README.md"]):
            with self.subTest(checkout_content=content):
                (self.checkout / "python/README.md").write_bytes(content)
                with self.assertRaisesRegex(ValueError, "package description"):
                    self.check_tagged_documentation("python/README.md")

    def test_publish_gate_keeps_other_documentation_updates_allowed(self):
        collector = self.check_tagged_documentation("docs/EXAMPLES.md")
        self.assertEqual(collector, self.checkout / "scripts/collect-release-wheels.py")

    def test_collector_does_not_replace_a_different_manifest(self):
        self.wheel("macos")
        path = self.directory / "WHEEL_MANIFEST.json"
        original = b'{"existing": "different"}\n'
        path.write_bytes(original)
        with self.assertRaises(ValueError):
            self.collect("--allow-incomplete")
        self.assertEqual(path.read_bytes(), original)

    def test_explicit_provenance_path_supported(self):
        _, proof = self.wheel()
        moved = self.root / "receipt.json"
        proof.rename(moved)
        self.collect("--allow-incomplete", "--provenance", str(moved))
        self.assertTrue((self.directory / "WHEEL_MANIFEST.json").is_file())

    def build_args(self):
        return types.SimpleNamespace(sdist=self.sdist, source_manifest=self.manifest, platform="macos", arch="arm64",
                                     output_dir=self.root / "built", dry_run=False, manylinux_image=None, python_version="3.14")

    def test_build_records_success_only_after_checked_cibuildwheel_exit(self):
        args = self.build_args()
        def successful(command, **kwargs):
            self.assertTrue(kwargs["check"])
            WHEELS.make_wheel(args.output_dir, family="macos", repaired=False, sources=self.sources)
        with mock.patch.object(BUILD.importlib.metadata, "version", return_value="4.2.1"), mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.subprocess, "run", side_effect=successful):
            BUILD.build(args)
        result = json.loads((args.output_dir / "wheel-build-provenance.json").read_text())
        self.assertEqual(result["build_system"], "cibuildwheel")
        self.assertEqual(result["source_commit"], SHA)
        self.assertEqual(result["native_debug_tests"], "passed")
        self.assertEqual(result["installed_wheel_tests"], "passed")

    def test_failed_cibuildwheel_has_no_success_provenance(self):
        args = self.build_args()
        with mock.patch.object(BUILD.importlib.metadata, "version", return_value="4.2.1"), mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.subprocess, "run", side_effect=subprocess.CalledProcessError(1, ["cibuildwheel"])):
            with self.assertRaises(subprocess.CalledProcessError):
                BUILD.build(args)
        self.assertFalse((args.output_dir / "wheel-build-provenance.json").exists())

    def test_linux_aarch64_python310_build_keeps_the_requested_target(self):
        args = self.build_args()
        args.platform, args.arch, args.python_version = "linux", "aarch64", "3.10"
        def successful(command, **kwargs):
            self.assertEqual(command[command.index("--only") + 1], "cp310-manylinux_aarch64")
            self.assertNotIn("CIBW_TEST_SKIP", kwargs["env"])
            WHEELS.make_wheel(args.output_dir, family="linux", architecture="aarch64", python_version="3.10",
                              repaired=False, sources=self.sources)
        with mock.patch.object(BUILD.importlib.metadata, "version", return_value="4.2.1"), mock.patch.object(BUILD.subprocess, "run", side_effect=successful):
            BUILD.build(args)
        proof = json.loads((args.output_dir / "wheel-build-provenance.json").read_text())
        self.assertEqual((proof["build_identifier"], proof["architecture"], proof["python_version"]),
                         ("cp310-manylinux_aarch64", "aarch64", "3.10"))
        self.assertEqual(proof["sdist_sha256"], self.source["sdist_sha256"])

    def test_build_rejects_another_python_output_without_a_success_receipt(self):
        args = self.build_args()
        args.python_version = "3.10"
        def wrong_target(command, **kwargs):
            WHEELS.make_wheel(args.output_dir, family="macos", python_version="3.14", repaired=False, sources=self.sources)
        with mock.patch.object(BUILD.importlib.metadata, "version", return_value="4.2.1"), mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.subprocess, "run", side_effect=wrong_target):
            with self.assertRaisesRegex(ValueError, "requested target"):
                BUILD.build(args)
        self.assertFalse((args.output_dir / "wheel-build-provenance.json").exists())

    def test_workflow_builds_all_twenty_targets_from_the_one_source_artifact(self):
        workflow = (SCRIPTS.parent / ".github/workflows/python-publish.yml").read_text()
        snippet = workflow.split("          python - <<'PY'\n", 1)[1].split("          PY\n", 1)[0]
        output = self.root / "matrix-output"
        with mock.patch.dict(os.environ, {"RELEASE_REF": SHA, "BUILD_TARGET": "all", "PUBLISH": "false", "GITHUB_OUTPUT": str(output)}), mock.patch.object(subprocess, "check_output", return_value=SHA + "\n"):
            exec(compile(textwrap.dedent(snippet), "workflow-matrix", "exec"), {})
        matrix = json.loads(output.read_text().removeprefix("value="))["include"]
        self.assertEqual(len(matrix), 20)
        self.assertEqual({(row["platform"], row["arch"], row["python"]) for row in matrix},
                         {(platform, architecture, version) for platform, architecture in
                          (("windows", "AMD64"), ("linux", "x86_64"), ("linux", "aarch64"), ("macos", "arm64"))
                          for version in ("3.10", "3.11", "3.12", "3.13", "3.14")})
        self.assertIn("name: python-release-source", workflow)
        self.assertIn('--python-version "$WHEEL_PYTHON_VERSION"', workflow)

    def test_docker_workflow_reuses_the_canonical_sdk_and_rejects_another_release_ref(self):
        workflow = (SCRIPTS.parent / ".github/workflows/release.yml").read_text()
        self.assertIn("name: python-release-source", workflow)
        self.assertIn("name: python-wheel-linux-x86_64-3.14", workflow)
        self.assertEqual(workflow.count("run-id: ${{ inputs.python_build_run_id }}"), 2)
        self.assertNotIn("build_python_wheels.py prepare", workflow)
        snippet = workflow.split("          python - <<'PY'\n")[2].split("          PY\n", 1)[0]
        directory = self.root / "python-source"
        directory.mkdir()
        shutil.copyfile(self.sdist, directory / self.sdist.name)
        shutil.copyfile(self.manifest, directory / self.manifest.name)
        environment_file = self.root / "github-env"
        loader = types.SimpleNamespace(loader=types.SimpleNamespace(exec_module=lambda module: None))
        original_directory = Path.cwd()
        try:
            os.chdir(self.root)
            with mock.patch.dict(os.environ, {"GITHUB_ENV": str(environment_file)}), mock.patch.object(importlib.util, "spec_from_file_location", return_value=loader), mock.patch.object(importlib.util, "module_from_spec", return_value=BUILD), mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git):
                exec(compile(textwrap.dedent(snippet), "docker-canonical-source", "exec"), {})
                self.assertEqual(BUILD.sha256(self.root / "release-input" / self.sdist.name), self.source["sdist_sha256"])
                self.assertEqual(environment_file.read_text(), "SMSD_REVISION=" + SHA + "\n")
                shutil.rmtree(self.root / "release-input")
                changed = {**self.source, "source_commit": "0" * 40, "release_ref": "0" * 40}
                (directory / "SOURCE_MANIFEST.json").write_text(json.dumps(changed))
                with self.assertRaisesRegex(SystemExit, "selected release ref differ"):
                    exec(compile(textwrap.dedent(snippet), "docker-canonical-source", "exec"), {})
                self.assertFalse((self.root / "release-input").exists())
        finally:
            os.chdir(original_directory)

    def test_build_checks_every_source_identity_field(self):
        args = self.build_args()
        for field, value in (("source_commit", "wrong"), ("release_ref", "v7.2.0"),
                             ("source_file_count", 999), ("working_tree_clean", False), ("sdist", "other.tar.gz")):
            with self.subTest(field=field):
                self.manifest.write_text(json.dumps({**self.source, field: value}))
                with self.assertRaises(ValueError):
                    BUILD.build(args)

    def local_args(self):
        wheel, proof = self.wheel(python_version=".".join(map(str, sys.version_info[:2])))
        proof.unlink()
        native = self.root / "native-debug"
        native.mkdir()
        return types.SimpleNamespace(sdist=self.sdist, wheel=wheel, source_manifest=self.manifest,
                                     native_test_dir=native, source_dir=self.checkout)

    def test_local_receipt_requires_all_three_checked_test_commands(self):
        args = self.local_args()
        listing = json.dumps({"tests": [{"name": str(index)} for index in range(BUILD.NATIVE_TEST_FLOOR)]})
        calls = []
        def checked(command, **kwargs):
            self.assertTrue(kwargs["check"])
            self.assertNotIn("PYTHONPATH", kwargs["env"])
            self.assertNotIn("SMSD_DISABLE_GPU", kwargs["env"])
            self.assertFalse((args.wheel.parent / "wheel-build-provenance.json").exists())
            calls.append(command)
        with mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.platform, "machine", return_value="arm64"), mock.patch.object(BUILD.subprocess, "check_output", return_value=listing), mock.patch.object(BUILD.subprocess, "run", side_effect=checked), mock.patch.object(BUILD, "verify_installed_wheel") as installed:
            BUILD.verify_local(args)
        self.assertEqual(len(calls), 5)
        self.assertEqual(calls[0][:2], ["cmake", "-S"])
        self.assertIn(str(self.checkout.resolve() / "cpp"), calls[0])
        self.assertIn("-DCMAKE_BUILD_TYPE=Debug", calls[0])
        self.assertIn("-DSMSD_BUILD_METAL=OFF", calls[0])
        self.assertIn("-DSMSD_BUILD_CUDA=OFF", calls[0])
        self.assertEqual(calls[1][:2], ["cmake", "--build"])
        self.assertEqual(calls[2][0], "ctest")
        self.assertIn("--require-rdkit", calls[3])
        self.assertIn("--require-openmp", calls[3])
        self.assertIn(str(self.checkout.resolve() / "python/tests"), calls[4])
        self.assertIn("--import-mode=importlib", calls[4])
        installed.assert_called_once_with(args.wheel)
        proof = json.loads((args.wheel.parent / "wheel-build-provenance.json").read_text())
        self.assertEqual(proof["build_system"], "local")
        self.assertNotIn("cibuildwheel", proof)
        self.assertNotIn(str(self.root), json.dumps(proof))
        self.collect("--allow-incomplete")

    def test_failed_local_gate_never_records_passed(self):
        args = self.local_args()
        listing = json.dumps({"tests": [{}] * BUILD.NATIVE_TEST_FLOOR})
        for failing_gate in range(5):
            with self.subTest(failing_gate=failing_gate):
                calls = []
                def failed(command, **kwargs):
                    calls.append(command)
                    if len(calls) == failing_gate + 1:
                        raise subprocess.CalledProcessError(1, command)
                with mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.platform, "machine", return_value="arm64"), mock.patch.object(BUILD.subprocess, "check_output", return_value=listing), mock.patch.object(BUILD.subprocess, "run", side_effect=failed), mock.patch.object(BUILD, "verify_installed_wheel") as installed:
                    with self.assertRaises(subprocess.CalledProcessError):
                        BUILD.verify_local(args)
                self.assertEqual(len(calls), failing_gate + 1)
                installed.assert_not_called()
                self.assertFalse((args.wheel.parent / "wheel-build-provenance.json").exists())

    def test_local_checks_reject_different_sdk_or_stale_build(self):
        args = self.local_args()
        for change in ("sdk", "stale"):
            with self.subTest(change=change):
                (args.source_dir / "python/smsd/search.py").write_bytes(
                    b"different\n" if change == "sdk" else self.sources["python/smsd/search.py"])
                if change == "stale":
                    (args.native_test_dir / "CMakeCache.txt").write_text("unrelated stale build\n")
                with mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.platform, "machine", return_value="arm64"), mock.patch.object(BUILD.subprocess, "check_output", return_value='{"tests": []}'), mock.patch.object(BUILD.subprocess, "run") as run:
                    with self.assertRaises(ValueError):
                        BUILD.verify_local(args)
                    run.assert_not_called()

    def test_local_checks_reject_empty_or_disabled_native_suites(self):
        args = self.local_args()
        for tests in ([], [{"properties": [{"name": "DISABLED", "value": True}]}] * BUILD.NATIVE_TEST_FLOOR):
            with self.subTest(tests=tests):
                with mock.patch.object(BUILD.platform, "system", return_value="Darwin"), mock.patch.object(BUILD.platform, "machine", return_value="arm64"), mock.patch.object(BUILD.subprocess, "check_output", return_value=json.dumps({"tests": tests})), mock.patch.object(BUILD.subprocess, "run") as run:
                    with self.assertRaises(ValueError):
                        BUILD.verify_local(args)
                    self.assertEqual(run.call_count, 2)
                self.assertFalse((args.wheel.parent / "wheel-build-provenance.json").exists())

    def test_partial_manifest_is_not_silently_replaced_when_platforms_added(self):
        self.wheel("macos")
        self.collect("--allow-incomplete")
        path = self.directory / "WHEEL_MANIFEST.json"
        original = path.read_bytes()
        self.wheel("linux")
        with self.assertRaisesRegex(ValueError, "new release directory"):
            self.collect("--allow-incomplete")
        self.assertEqual(path.read_bytes(), original)
        self.assertFalse((self.directory / self.targets["cp314-manylinux_x86_64"]).exists())

    def test_installed_package_must_equal_the_selected_wheel(self):
        wheel, _ = self.wheel()
        installed = self.root / "installed"
        with zipfile.ZipFile(wheel) as archive:
            archive.extractall(installed)
        distribution = types.SimpleNamespace(locate_file=lambda name: installed / name)
        with mock.patch.object(BUILD.importlib.metadata, "distribution", return_value=distribution):
            BUILD.verify_installed_wheel(wheel)
            (installed / "smsd/search.py").write_text("different package\n")
            with self.assertRaises(ValueError):
                BUILD.verify_installed_wheel(wheel)

    def test_collector_rejects_conflicting_receipts(self):
        _, proof = self.wheel()
        conflicting = self.root / "other-receipt.json"
        other = json.loads(proof.read_text())
        other["build_system"] = "local"
        other.pop("cibuildwheel")
        conflicting.write_text(json.dumps(other))
        with self.assertRaises(ValueError):
            self.collect("--allow-incomplete", "--provenance", str(conflicting))

    def test_manifests_reject_duplicate_json_keys(self):
        self.manifest.write_text('{"version":"7.2.0","version":"7.2.0"}')
        with self.assertRaises(ValueError):
            BUILD.checked_source(self.sdist)


if __name__ == "__main__":
    unittest.main()

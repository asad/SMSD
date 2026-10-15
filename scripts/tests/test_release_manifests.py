"""Source and wheel provenance producer/consumer regression tests."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
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
            "pyproject.toml": b'[project]\nname = "smsd"\nversion = "7.2.0"\n',
            "pom.xml": b'<project xmlns="http://maven.apache.org/POM/4.0.0"><version>7.2.0</version></project>',
            "java/pom.xml": b'<project xmlns="http://maven.apache.org/POM/4.0.0"><version>7.2.0</version></project>',
            "cpp/CMakeLists.txt": b'project(smsd VERSION 7.2.0 LANGUAGES CXX)\n',
            "python/README.md": b'Fixture package description\n',
            "PKG-INFO": b'Name: smsd\nVersion: 7.2.0\n\nFixture package description\n',
        }
        for name in ("cpp/include/smsd/mcs.hpp", "cpp/bindings/pybind11/smsd_bindings.cpp",
                     "python/tests/test_smsd.py", "python/tests/test_aromaticity_parity.py",
                     "java/src/test/resources/com/bioinception/smsd/aromaticity_parity.json",
                     "scripts/check_python_wheel.py", "scripts/check_native_wheel_build.py",
                     "scripts/build_python_wheels.py", "scripts/collect-release-wheels.py", "scripts/cibuildwheel.toml",
                     "scripts/repair_windows_wheel.py", "licenses/msvc/README.md",
                     "licenses/msvc/LICENSE-2022.docx", "licenses/msvc/LICENSE-2026.docx"):
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
        self.families = {}

    def git(self, command, **kwargs):
        arguments = command[1:]
        if "status" in arguments:
            return ""
        if "rev-parse" in arguments:
            return SHA + "\n"
        if "ls-tree" in arguments:
            entries = []
            for name, data in self.sources.items():
                if name == "PKG-INFO":
                    continue
                blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
                entries.append(f"100644 blob {blob}\t{name}".encode())
            return b"\0".join(entries) + b"\0"
        raise AssertionError(command)

    def wheel(self, family="macos", **proof_changes):
        directory = self.root / family
        directory.mkdir(exist_ok=True)
        wheel = WHEELS.make_wheel(directory, family=family, repaired=False, sources=self.sources,
                                 changes={WHEELS.PREFIX + "METADATA": self.sources["PKG-INFO"],
                                          "include/smsd/mcs.hpp": self.sources["cpp/include/smsd/mcs.hpp"]})
        identifier = {"macos": "cp314-macosx_arm64", "linux": "cp314-manylinux_x86_64", "windows": "cp314-win_amd64"}[family]
        proof = {**self.source, "wheel": wheel.name, "wheel_sha256": BUILD.sha256(wheel),
                 "build_identifier": identifier, "build_system": "cibuildwheel", "cibuildwheel": "4.2.1",
                 "native_debug_tests": "passed", "installed_wheel_tests": "passed"}
        proof.update(proof_changes)
        path = directory / "wheel-build-provenance.json"
        path.write_text(json.dumps(proof))
        self.proofs.append(path)
        self.families[family] = wheel.name
        return wheel, path

    def collect(self, *extra):
        arguments = ["collect-release-wheels.py", "--release-dir", str(self.directory)]
        for family in self.families:
            arguments += ["--wheel-dir", str(self.root / family)]
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

    def test_prepare_rejects_sdist_file_not_in_frozen_tree(self):
        output = self.root / "different"
        args = types.SimpleNamespace(output_dir=output, release_ref="v7.2.0", no_isolation=True)
        def changed_source(command, **kwargs):
            archive(output, {**self.sources, "python/smsd/search.py": b'changed after freeze\n'})
        with mock.patch.object(BUILD.subprocess, "check_output", side_effect=self.git), mock.patch.object(BUILD.subprocess, "run", side_effect=changed_source):
            with self.assertRaises(ValueError):
                BUILD.prepare(args)
        self.assertFalse((output / "SOURCE_MANIFEST.json").exists())

    def test_collector_rejects_missing_provenance(self):
        _, proof = self.wheel()
        proof.unlink()
        with self.assertRaises(ValueError):
            self.collect("--allow-incomplete")
        self.assertFalse((self.directory / "WHEEL_MANIFEST.json").exists())

    def test_collector_emits_all_platforms_from_valid_receipts(self):
        for family in ("linux", "macos", "windows"):
            self.wheel(family)
        self.collect()
        result = json.loads((self.directory / "WHEEL_MANIFEST.json").read_text())
        self.assertTrue(result["complete"])
        self.assertEqual({entry["platform"] for entry in result["wheels"]}, set(self.families))
        GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.families, self.sources["PKG-INFO"])
        self.assertNotIn(str(self.root), json.dumps(result))

    def test_collector_rejects_mismatched_or_unpassed_receipts(self):
        wheel, proof = self.wheel()
        original = json.loads(proof.read_text())
        for field, value in (("wheel_sha256", "0" * 64), ("source_commit", "0" * 40),
                             ("sdist_sha256", "0" * 64), ("release_ref", "v7.2.0"),
                             ("native_debug_tests", "not run"), ("installed_wheel_tests", None),
                             ("build_identifier", "cp314-win_amd64")):
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
        self.assertEqual(result["missing_platforms"], ["linux", "windows"])
        self.assertEqual([entry["platform"] for entry in result["wheels"]], ["macos"])

    def test_check_only_accepts_published_manifest_without_sidecars(self):
        for family in ("linux", "macos", "windows"):
            self.wheel(family)
        self.collect()
        for proof in self.proofs:
            proof.unlink()
        before = {path.name: path.read_bytes() for path in self.directory.iterdir()}
        self.families.clear()
        self.collect("--check-only")
        self.assertEqual(before, {path.name: path.read_bytes() for path in self.directory.iterdir()})

    def test_check_only_rejects_hand_assembled_passed_strings(self):
        for family in ("linux", "macos", "windows"):
            self.wheel(family)
        self.collect()
        path = self.directory / "WHEEL_MANIFEST.json"
        manifest = json.loads(path.read_text())
        for entry in manifest["wheels"]:
            entry.pop("build_provenance", None)
        path.write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            self.collect("--check-only")
        with self.assertRaises((ValueError, SystemExit)):
            GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.families, self.sources["PKG-INFO"])

    def test_publish_gate_rejects_status_only_test_evidence(self):
        entries = []
        for family in ("linux", "macos", "windows"):
            wheel, _ = self.wheel(family)
            (self.directory / wheel.name).write_bytes(wheel.read_bytes())
            entries.append({"platform": family, "wheel": wheel.name,
                            "wheel_sha256": BUILD.sha256(wheel), "native_debug_tests": "passed",
                            "installed_wheel_tests": "passed"})
        manifest = {"schema_version": 1, "complete": True, "missing_platforms": [], "wheels": entries,
                    **{field: self.source[field] for field in ("version", "source_commit", "sdist", "sdist_sha256")}}
        (self.directory / "WHEEL_MANIFEST.json").write_text(json.dumps(manifest))
        with self.assertRaises(ValueError):
            GATE.check_wheel_manifest(self.directory, "7.2.0", self.source, self.families, self.sources["PKG-INFO"])

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
                                     output_dir=self.root / "built", dry_run=False, manylinux_image=None)

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

    def test_build_checks_every_source_identity_field(self):
        args = self.build_args()
        for field, value in (("source_commit", "wrong"), ("release_ref", "v7.2.0"),
                             ("source_file_count", 999), ("working_tree_clean", False), ("sdist", "other.tar.gz")):
            with self.subTest(field=field):
                self.manifest.write_text(json.dumps({**self.source, field: value}))
                with self.assertRaises(ValueError):
                    BUILD.build(args)

    def local_args(self):
        wheel, proof = self.wheel()
        proof.unlink()
        native = self.root / "native-debug"
        native.mkdir()
        return types.SimpleNamespace(sdist=self.sdist, wheel=wheel, source_manifest=self.manifest,
                                     native_test_dir=native, source_dir=self.checkout)

    def test_local_receipt_requires_all_three_checked_test_commands(self):
        args = self.local_args()
        listing = json.dumps({"tests": [{"name": str(index)} for index in range(12)]})
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
        listing = json.dumps({"tests": [{}] * 12})
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
        for tests in ([], [{"properties": [{"name": "DISABLED", "value": True}]}] * 12):
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
        self.assertFalse((self.directory / self.families["linux"]).exists())

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

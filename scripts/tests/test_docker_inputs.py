# SPDX-License-Identifier: Apache-2.0
"""Offline Docker input and callback-observation checks; no image or GPU is needed."""
import contextlib
import ctypes
import importlib.util
import io
import json
from pathlib import Path
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


PREPARE = load("docker_prepare", SCRIPTS / "prepare_docker_context.py")
INSTALL = load("docker_install", SCRIPTS / "docker_install.py")
SMOKE = load("docker_smoke", SCRIPTS / "docker_smoke.py")
FIXTURES = load("docker_release_fixtures", Path(__file__).with_name("test_release_manifests.py"))


class DockerInputsTest(unittest.TestCase):
    def setUp(self):
        fixture = FIXTURES.ReleaseManifestsTest()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.root = fixture.root
        self.release = fixture.directory
        sources = {name: data.replace(b"7.2.0", b"8.0.2") for name, data in fixture.sources.items()}
        old = FIXTURES.archive(self.release, sources)
        self.sdist = old.with_name("smsd-8.0.2.tar.gz")
        old.rename(self.sdist)
        source = {**PREPARE.BUILD.inspect_source(self.sdist), "source_commit": FIXTURES.SHA,
                  "release_ref": FIXTURES.SHA, "working_tree_clean": True}
        (self.release / "SOURCE_MANIFEST.json").write_text(json.dumps(source))
        entries = []
        for architecture in ("x86_64", "aarch64"):
            tag = "cp314-manylinux_" + architecture
            wheel = self.release / ("smsd-8.0.2-cp314-cp314-manylinux_2_28_" + architecture + ".whl")
            with zipfile.ZipFile(wheel, "w") as archive:
                archive.writestr("smsd-8.0.2.dist-info/METADATA", "Name: smsd\nVersion: 8.0.2\n")
            proof = {**source, **PREPARE.BUILD.BUILD_IDENTIFIERS[tag], "wheel": wheel.name,
                     "wheel_sha256": PREPARE.BUILD.sha256(wheel), "build_identifier": tag,
                     "build_system": "cibuildwheel", "cibuildwheel": "4.2.1",
                     "native_debug_tests": "passed", "installed_wheel_tests": "passed"}
            entries.append(PREPARE.BUILD.wheel_entry(wheel, tag, source, proof))
        manifest = {key: source[key] for key in ("version", "source_commit", "sdist", "sdist_sha256")}
        manifest.update(schema_version=2, wheels=entries, complete=False,
                        missing_targets=sorted(set(PREPARE.BUILD.BUILD_IDENTIFIERS) - {entry["build_identifier"] for entry in entries}))
        (self.release / "WHEEL_MANIFEST.json").write_text(json.dumps(manifest))
        self.args = types.SimpleNamespace(release_dir=self.release, output=self.root / "context",
            architecture=None, cuda_wheel=None, cuda_provenance=None, runtime_wheel=None)

    def prepare(self):
        with contextlib.redirect_stdout(io.StringIO()):
            PREPARE.prepare(self.args)
        return json.loads((self.args.output / "ARTIFACTS.json").read_text())

    def test_context_uses_frozen_files_and_each_matching_architecture(self):
        manifest = self.prepare()
        self.assertEqual((self.args.output / "Dockerfile").read_bytes(), b"fixture\n")
        self.assertEqual((self.args.output / "source.tar.gz").read_bytes(), self.sdist.read_bytes())
        for architecture, platform in (("amd64", "x86_64"), ("arm64", "aarch64")):
            entries = list(INSTALL.selected_files(self.args.output, manifest, architecture, False))
            self.assertEqual(len(entries), 1)
            self.assertTrue(entries[0][1].name.endswith(platform + ".whl"))
        with mock.patch("sys.argv", ["docker_install.py", "--root", str(self.args.output),
                "--architecture", "amd64", "--revision", FIXTURES.SHA, "--verify-only"]):
            INSTALL.main()

    def test_changed_wheel_is_rejected_before_context_creation(self):
        wheel = next(self.release.glob("*.whl"))
        wheel.write_bytes(wheel.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "SHA256"):
            self.prepare()
        self.assertFalse(self.args.output.exists())

    def test_core_manylinux_alias_order_keeps_required_architecture_and_floor(self):
        manifest = self.prepare()
        entry = manifest["cpu"]["amd64"]
        original = self.args.output / entry["path"]
        for tags, accepted in (("manylinux_2_17_x86_64.manylinux_2_28_x86_64", True),
                               ("manylinux_2_28_x86_64.manylinux_2_17_x86_64", True),
                               ("manylinux_2_17_x86_64", False),
                               ("manylinux_2_28_x86_64.manylinux_2_28_aarch64", False)):
            with self.subTest(tags=tags):
                entry["path"] = str(Path(entry["path"]).with_name("smsd-8.0.2-cp314-cp314-" + tags + ".whl"))
                target = self.args.output / entry["path"]
                target.write_bytes(original.read_bytes())
                if accepted:
                    self.assertEqual(len(list(INSTALL.selected_files(self.args.output, manifest, "amd64", False))), 1)
                else:
                    with self.assertRaisesRegex(ValueError, "CPython3.14"):
                        list(INSTALL.selected_files(self.args.output, manifest, "amd64", False))

    def test_source_filename_cannot_escape_release_directory(self):
        path = self.release / "SOURCE_MANIFEST.json"
        manifest = json.loads(path.read_text())
        manifest["sdist"] = "../smsd-8.0.2.tar.gz"
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "basename"):
            self.prepare()
        self.assertFalse(self.args.output.exists())

    def test_image_rejects_wrong_interpreter_and_missing_gpu_pair(self):
        manifest = self.prepare()
        original = manifest["cpu"]["amd64"]["path"]
        changed = original.replace("cp314-cp314", "cp310-cp310")
        (self.args.output / original).rename(self.args.output / changed)
        manifest["cpu"]["amd64"]["path"] = changed
        with self.assertRaisesRegex(ValueError, "CPython3.14"):
            list(INSTALL.selected_files(self.args.output, manifest, "amd64", False))
        with self.assertRaisesRegex(ValueError, "matching core"):
            list(INSTALL.selected_files(self.args.output, manifest, "amd64", True))
        with self.assertRaisesRegex(ValueError, "architecture"):
            list(INSTALL.selected_files(self.args.output, manifest, "arm64", True))

    def test_existing_context_is_not_overwritten(self):
        self.prepare()
        with self.assertRaises(FileExistsError):
            self.prepare()


class CallbackObserverTest(unittest.TestCase):
    def test_observation_preserves_minimum_and_delegates_real_arguments_status(self):
        calls = []
        @SMOKE.MinimumTargets
        def minimum(operation, value, error, size):
            calls.append(operation)
            value[0] = (1 << 64) - 1
            return 0
        @SMOKE.FingerprintScreen
        def screen(resident, count, words, cutoffs, stride, hits, capacity, offsets, error, size):
            calls.append((resident, count, stride, capacity))
            return 9
        table = SMOKE.AcceleratorTable(abi_version=1, struct_size=ctypes.sizeof(SMOKE.AcceleratorTable))
        for name, callback in (("minimum_targets", minimum), ("fingerprint_screen", screen)):
            setattr(table, name, ctypes.cast(callback, ctypes.c_void_p).value)
        table.create_fingerprints = table.destroy_resident = 1
        library = types.SimpleNamespace(smsd_cuda_screen_accelerator=lambda: ctypes.addressof(table))
        observer = SMOKE.CallbackObserver(library)
        self.assertEqual(observer.minimum.value, (1 << 64) - 1)
        self.assertEqual(observer.table.minimum_targets, table.minimum_targets)
        status = observer.callback(123, 2, None, None, 2049, None, 15, None, None, 0)
        self.assertEqual((status, observer.requests, observer.successes), (9, 1, 0))
        self.assertEqual(calls, [1, (123, 2, 2049, 15)])


if __name__ == "__main__":
    unittest.main()

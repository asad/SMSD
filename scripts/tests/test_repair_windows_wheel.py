"""Local tests for Windows runtime discovery, version checks and repair inputs."""

import importlib.util
import io
from pathlib import Path
import struct
import sys
import tempfile
import types
import unittest
from unittest import mock
import zipfile


SCRIPT = Path(__file__).resolve().parents[1] / "repair_windows_wheel.py"
SPEC = importlib.util.spec_from_file_location("repair_windows_wheel", SCRIPT)
REPAIR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = REPAIR
SPEC.loader.exec_module(REPAIR)


def pe_image(linker=(14, 44), machine=0x8664):
    image = bytearray(100)
    image[:2] = b"MZ"
    struct.pack_into("<I", image, 60, 64)
    image[64:68] = b"PE\0\0"
    struct.pack_into("<H", image, 68, machine)
    struct.pack_into("<H", image, 88, 0x20B)
    image[90:92] = bytes(linker)
    return bytes(image)


class WindowsRuntimeRepairTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="smsd-runtime-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        self.info = {}

    def dll(self, directory, name, version, company="Microsoft Corporation", flags=0, machine=0x8664):
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / name
        path.write_bytes(pe_image(machine=machine))
        self.info[path] = REPAIR.DllInfo(version, company, flags)
        return path

    def select(self, directories):
        return REPAIR.select_runtimes(directories, (14, 44), self.info.__getitem__)

    def wheel(self, extra=None):
        path = self.root / "smsd-7.2.0-cp314-cp314-win_amd64.whl"
        with zipfile.ZipFile(path, "w") as archive:
            archive.writestr("smsd/_smsd.cp314-win_amd64.pyd", pe_image())
            for name, data in (extra or {}).items():
                archive.writestr(name, data)
        return path

    def test_baseline_comes_from_built_extension(self):
        self.assertEqual(REPAIR.wheel_runtime_baseline(self.wheel()), (14, 44))

    def test_rejects_invalid_architecture_and_already_repaired_input(self):
        with self.assertRaisesRegex(ValueError, "AMD64"):
            REPAIR.pe_linker_version(pe_image(machine=0x14C))
        with self.assertRaisesRegex(ValueError, "PE image"):
            REPAIR.pe_linker_version(b"not a PE")
        with self.assertRaisesRegex(ValueError, "previously bundled"):
            REPAIR.wheel_runtime_baseline(self.wheel({"smsd.libs/msvcp140-" + "a" * 32 + ".dll": pe_image()}))

    def test_accepts_supported_family_without_comparing_compiler_patch(self):
        # MSVC 19.44.35229, VCTools 14.44.35207 and redist 14.44.35211
        # have different patch fields. The baseline is the 14.44 family.
        REPAIR.validate_runtime(REPAIR.DllInfo((14, 44, 35211, 0), "Microsoft Corporation"), (14, 44))
        REPAIR.validate_runtime(REPAIR.DllInfo((14, 51, 36231, 0), "Microsoft Corporation"), (14, 44))
        for version in ((14, 40, 33810, 0), (13, 99, 99999, 0), (15, 44, 0, 0)):
            with self.subTest(version=version), self.assertRaisesRegex(ValueError, "required family"):
                REPAIR.validate_runtime(REPAIR.DllInfo(version, "Microsoft Corporation"), (14, 44))

    def test_prefers_supported_visual_studio_redist_over_system32(self):
        redist, system32 = self.root / "redist", self.root / "System32"
        for name in REPAIR.RUNTIME_NAMES:
            self.dll(redist, name, (14, 44, 35211, 0))
            self.dll(system32, name, (14, 51, 36231, 0))
        selected = self.select([(system32, "Windows System32"), (redist, "Visual Studio redist")])
        self.assertTrue(all(path.parent == redist for path, _, _ in selected.values()))

    def test_selects_newest_redist_and_can_fall_back_per_file(self):
        old, new, system32 = (self.root / name for name in ("old", "new", "System32"))
        self.dll(old, "msvcp140.dll", (14, 40, 33810, 0))
        self.dll(new, "msvcp140.dll", (14, 44, 35211, 0))
        self.dll(system32, "vcomp140.dll", (14, 51, 36231, 0))
        selected = self.select([(old, "Visual Studio redist"), (new, "Visual Studio redist"),
                                (system32, "Windows System32")])
        self.assertEqual(selected["msvcp140.dll"][0].parent, new)
        self.assertEqual(selected["vcomp140.dll"][0].parent, system32)

    def test_rejects_non_microsoft_prerelease_debug_and_wrong_architecture(self):
        for kwargs in ({"company": "Other vendor"}, {"flags": 1}, {"flags": 2}, {"machine": 0xAA64}):
            with self.subTest(kwargs=kwargs):
                directory = self.root / "invalid"
                self.dll(directory, "msvcp140.dll", (14, 44, 35211, 0), **kwargs)
                with self.assertRaisesRegex(ValueError, "No supported AMD64 Microsoft"):
                    self.select([(directory, "Visual Studio redist")])

    def test_discovery_excludes_debug_x86_and_arbitrary_path(self):
        installation = self.root / "VS"
        release = installation / "VC/Redist/MSVC/14.44.35211"
        crt, omp = release / "x64/Microsoft.VC143.CRT", release / "x64/Microsoft.VC143.OpenMP"
        for directory in (crt, omp, release / "x86/Microsoft.VC143.CRT",
                          release / "debug_nonredist/x64/Microsoft.VC143.DebugCRT"):
            directory.mkdir(parents=True)
        system32 = self.root / "System32"
        directories = REPAIR.runtime_directories({"PATH": str(self.root / "Java"),
                                                 "VCToolsRedistDir": str(release)},
                                                [installation], system32)
        self.assertEqual(directories, [(crt, "Visual Studio redist"), (omp, "Visual Studio redist"),
                                       (system32, "Windows System32")])

    def test_developer_environment_and_vswhere_discovery(self):
        installation = self.root / "VS"
        program_files = self.root / "Program Files (x86)"
        vswhere = program_files / "Microsoft Visual Studio/Installer/vswhere.exe"
        vswhere.parent.mkdir(parents=True)
        vswhere.write_bytes(b"fixture")
        environment = {"VSINSTALLDIR": str(installation), "VCINSTALLDIR": str(installation / "VC"),
                       "VCToolsInstallDir": str(installation / "VC/Tools/MSVC/14.44.35207"),
                       "ProgramFiles(x86)": str(program_files)}
        with mock.patch.object(REPAIR.subprocess, "check_output", return_value=str(installation) + "\n") as call:
            self.assertEqual(REPAIR.visual_studio_installations(environment), [installation])
        self.assertIn("-requires", call.call_args.args[0])
        self.assertNotIn("-prerelease", call.call_args.args[0])

    def test_verifies_exact_selected_bytes_and_future_runtime_versions(self):
        original = {name: pe_image() for name in REPAIR.RUNTIME_NAMES}
        members = {"smsd.libs/" + name[:-4] + "-" + "a" * 32 + ".dll": data
                   for name, data in original.items()}
        REPAIR.verify_bundled_runtimes(self.wheel(members), original, (14, 44))
        members["smsd.libs/msvcp140_1-" + "b" * 32 + ".dll"] = pe_image()
        with self.assertRaisesRegex(ValueError, "required family"):
            REPAIR.verify_bundled_runtimes(self.wheel(members), original, (14, 44),
                                          lambda path: REPAIR.DllInfo((14, 40, 33810, 0), "Microsoft Corporation"))
        members["smsd.libs/msvcp140-" + "a" * 32 + ".dll"] += b"altered"
        with self.assertRaisesRegex(ValueError, "selected, unmodified"):
            REPAIR.verify_bundled_runtimes(self.wheel(members), original, (14, 44))

    def test_repair_stages_selected_files_before_ambient_path(self):
        selected = {}
        for name in REPAIR.RUNTIME_NAMES:
            path = self.dll(self.root / name[:-4], name, (14, 44, 35211, 0))
            selected[name] = (path, self.info[path], "Visual Studio redist")
        wheel, destination = self.wheel(), self.root / "output"

        def fake_delvewheel(command, check):
            self.assertTrue(check)
            self.assertEqual(command[:4], [sys.executable, "-m", "delvewheel", "repair"])
            staging = Path(command[command.index("--add-path") + 1])
            self.assertEqual(set(path.name for path in staging.iterdir()), set(REPAIR.RUNTIME_NAMES))
            with zipfile.ZipFile(destination / wheel.name, "w") as archive:
                for name, (original, _, _) in selected.items():
                    self.assertEqual((staging / name).read_bytes(), original.read_bytes())
                    archive.writestr("smsd.libs/" + name[:-4] + "-" + "a" * 32 + ".dll", original.read_bytes())

        with mock.patch.object(REPAIR, "os", types.SimpleNamespace(name="nt", environ={})), \
                mock.patch.object(REPAIR.importlib.metadata, "version", return_value="1.13.1"), \
                mock.patch.object(REPAIR, "visual_studio_installations", return_value=[]), \
                mock.patch.object(REPAIR, "system_directory", return_value=self.root / "System32"), \
                mock.patch.object(REPAIR, "select_runtimes", return_value=selected), \
                mock.patch.object(REPAIR.subprocess, "run", side_effect=fake_delvewheel), \
                mock.patch("sys.stdout", new_callable=io.StringIO) as output:
            REPAIR.repair(wheel, destination)
        self.assertIn('"version": "14.44.35211.0"', output.getvalue())
        self.assertIn("Verified both bundled Microsoft runtimes", output.getvalue())


if __name__ == "__main__":
    unittest.main()

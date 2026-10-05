"""Release-wheel integrity tests; run with python -m unittest discover -s scripts/tests."""

import base64
import csv
import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
import zipfile


SCRIPT = Path(__file__).resolve().parents[1] / "collect-release-wheels.py"
SPEC = importlib.util.spec_from_file_location("collect_release_wheels", SCRIPT)
COLLECTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COLLECTOR)

PREFIX = "smsd-7.2.0.dist-info/"
ORIGINAL = b'"""SMSD fixture."""\n\nimport os\n\n__version__ = "7.2.0"\n'
BOOTSTRAP = b"""# start delvewheel patch
def _delvewheel_patch_1_13_1():
    import os
    if os.path.isdir(libs_dir := os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'smsd.libs'))):
        os.add_dll_directory(libs_dir)


_delvewheel_patch_1_13_1()
del _delvewheel_patch_1_13_1
# end delvewheel patch
"""
REPAIRED = b'"""SMSD fixture."""\n\n\n' + BOOTSTRAP + b'\nimport os\n\n__version__ = "7.2.0"\n'
SOURCES = {
    "python/smsd/__init__.py": ORIGINAL,
    "python/smsd/mcs_engine.py": b"def find_mcs():\n    return {}\n",
    "cpp/include/smsd/smsd.hpp": b"#pragma once\n",
    "LICENSE": b"package license\n",
    "NOTICE": b"package notice\n",
    "licenses/libomp/LICENSE.TXT": b"libomp license\n",
    "licenses/libgomp/COPYING3": b"libgomp license\n",
    "licenses/libgomp/COPYING.RUNTIME": b"libgomp exception\n",
}


def make_wheel(directory, family="windows", changes=None, repaired=True, sources=None):
    sources = SOURCES if sources is None else sources
    if family == "windows":
        tag = "win_amd64"
        extension = "smsd/_smsd.cp314-win_amd64.pyd"
        binary = bytearray(70)
        binary[:2] = b"MZ"
        binary[60:64] = (64).to_bytes(4, "little")
        binary[64:70] = b"PE\0\0\x64\x86"
    elif family == "linux":
        tag = "manylinux_2_28_x86_64"
        extension = "smsd/_smsd.cpython-314-x86_64-linux-gnu.so"
        binary = bytearray(20)
        binary[:6] = b"\x7fELF\x02\x01"
        binary[18:20] = b"\x3e\0"
    else:
        tag = "macosx_26_0_arm64"
        extension = "smsd/_smsd.cpython-314-darwin.so"
        binary = b"\xcf\xfa\xed\xfe\x0c\0\0\x01"
    contents = {
        "smsd/__init__.py": REPAIRED if repaired else ORIGINAL,
        "smsd/mcs_engine.py": sources["python/smsd/mcs_engine.py"],
        "include/smsd/smsd.hpp": sources["cpp/include/smsd/smsd.hpp"],
        extension: bytes(binary),
        PREFIX + "METADATA": b"Name: smsd\nVersion: 7.2.0\n",
        PREFIX + "WHEEL": f"Tag: cp314-cp314-{tag}\n".encode(),
    }
    for name, data in sources.items():
        if not name.startswith(("python/", "cpp/")):
            contents[PREFIX + "licenses/" + name] = data
    if repaired:
        contents[PREFIX + "DELVEWHEEL"] = b"Version: 1.13.1\nArguments: ['delvewheel', 'repair']\n"
        contents["smsd.libs/runtime.dll"] = b"DLL fixture"
    for name, data in (changes or {}).items():
        if data is None:
            contents.pop(name, None)
        else:
            contents[name] = data
    # Regenerate valid hashes after mutations so source/template checks are tested.
    record = io.StringIO(newline="")
    writer = csv.writer(record, lineterminator="\n")
    for name, data in contents.items():
        if name.endswith("/"):
            continue
        checksum = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
        writer.writerow((name, "sha256=" + checksum, len(data)))
    writer.writerow((PREFIX + "RECORD", "", ""))
    contents[PREFIX + "RECORD"] = record.getvalue().encode()
    path = Path(directory) / f"smsd-7.2.0-cp314-cp314-{tag}.whl"
    with zipfile.ZipFile(path, "w") as archive:
        for name, data in contents.items():
            archive.writestr(name, data)
    return path


class CollectReleaseWheelsTest(unittest.TestCase):
    def check(self, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            wheel = make_wheel(directory, **kwargs)
            return COLLECTOR.check_wheel(wheel, "7.2.0", kwargs.get("sources", SOURCES))

    def test_supported_windows_repair(self):
        self.assertEqual(self.check(), "windows")

    def test_supported_insertion_preserves_utf8_and_future_imports(self):
        cases = (
            ('"""SMSD é fixture."""\n\nimport os\n'.encode(),
             '"""SMSD é fixture."""\n\n\n'.encode() + BOOTSTRAP + b"\nimport os\n"),
            (b'"""SMSD fixture."""\nfrom __future__ import annotations\n\nimport os\n',
             b'"""SMSD fixture."""\nfrom __future__ import annotations\n\n\n' + BOOTSTRAP + b"\nimport os\n"),
        )
        with tempfile.TemporaryDirectory() as directory:
            with zipfile.ZipFile(make_wheel(directory)) as wheel:
                for original, repaired in cases:
                    with self.subTest(original=original):
                        self.assertTrue(COLLECTOR.matches_windows_repaired_init(wheel, PREFIX, original, repaired))

    def test_original_source_stays_strict_on_all_platforms(self):
        for family in ("windows", "linux", "macos"):
            with self.subTest(family=family):
                self.assertEqual(self.check(family=family, repaired=False), family)

    def test_application_changes_rejected_with_valid_record(self):
        for changed in (REPAIRED.replace(b'__version__ = "7.2.0"', b'__version__ = "9.0.0"'),
                        REPAIRED.replace(b"\nimport os\n", b"\nimport  os\n"),
                        REPAIRED + b"print('extra application code')\n"):
            with self.subTest(changed=changed):
                with self.assertRaisesRegex(ValueError, "Wheel and source differ"):
                    self.check(changes={"smsd/__init__.py": changed})

    def test_altered_bootstrap_rejected_with_valid_record(self):
        for changed in (
            REPAIRED.replace(b"os.add_dll_directory(libs_dir)", b"exec('malicious code')"),
            REPAIRED.replace(b"'smsd.libs'", b"'other.libs'"),
            REPAIRED.replace(b"# end delvewheel patch", b"    exec('malicious code')\n# end delvewheel patch"),
            REPAIRED.replace(BOOTSTRAP, BOOTSTRAP + BOOTSTRAP),
            REPAIRED.replace(BOOTSTRAP, b"# start delvewheel patch\nexec('malicious code')\n# end delvewheel patch\n"),
            ORIGINAL + b"\n" + BOOTSTRAP,
        ):
            with self.subTest(changed=changed):
                with self.assertRaisesRegex(ValueError, "Wheel and source differ"):
                    self.check(changes={"smsd/__init__.py": changed})

    def test_repair_requires_supported_metadata_and_bundled_dll(self):
        for changes in (
            {PREFIX + "DELVEWHEEL": None},
            {PREFIX + "DELVEWHEEL": b"Version: 1.13.2\n"},
            {PREFIX + "DELVEWHEEL": b"Arguments: []\nVersion: 1.13.1\n"},
            {PREFIX + "DELVEWHEEL": b"Version: 1.13.1 \n"},
            {"smsd.libs/runtime.dll": None, "smsd.libs/": b""},
            {"smsd.libs/runtime.dll": None, "smsd.libs/nested/runtime.dll": b"DLL fixture"},
        ):
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(ValueError, "Wheel and source differ"):
                    self.check(changes=changes)

    def test_repair_exception_is_windows_only(self):
        for family in ("linux", "macos"):
            with self.subTest(family=family):
                with self.assertRaisesRegex(ValueError, "Wheel and source differ"):
                    self.check(family=family)

    def test_other_python_and_headers_stay_byte_exact(self):
        for name in ("smsd/mcs_engine.py", "include/smsd/smsd.hpp"):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, "Wheel and source differ"):
                    self.check(changes={name: b"modified contents\n"})

    def test_license_metadata_tags_and_architecture_checks_remain(self):
        cases = (
            ({PREFIX + "licenses/LICENSE": b"modified license\n"}, "Wrong license copy"),
            ({PREFIX + "METADATA": b"Name: smsd\nVersion: 9.0.0\n"}, "Wrong wheel metadata"),
            ({PREFIX + "WHEEL": b"Tag: cp314-cp314-win32\n"}, "Filename and wheel tags disagree"),
            ({"smsd/_smsd.cp314-win_amd64.pyd": b"wrong architecture"}, "Native extension does not match"),
        )
        for changes, message in cases:
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(ValueError, message):
                    self.check(changes=changes)

    def test_required_msvc_license_copies_cannot_be_missing_or_altered(self):
        msvc_licenses = {
            "licenses/msvc/README.md": b"MSVC runtime license notice\n",
            "licenses/msvc/LICENSE-2022.docx": b"MSVC 2022 license fixture",
            "licenses/msvc/LICENSE-2026.docx": b"MSVC 2026 license fixture",
        }
        sources = {**SOURCES, **msvc_licenses}
        self.assertEqual(self.check(sources=sources), "windows")
        for name in msvc_licenses:
            entry = PREFIX + "licenses/" + name
            with self.subTest(name=name, mutation="missing"):
                with self.assertRaises(KeyError):
                    self.check(sources=sources, changes={entry: None})
            with self.subTest(name=name, mutation="altered"):
                with self.assertRaisesRegex(ValueError, "Wrong license copy"):
                    self.check(sources=sources, changes={entry: b"modified license contents"})

    def test_record_hash_stays_strict(self):
        with tempfile.TemporaryDirectory() as directory:
            wheel = make_wheel(directory)
            with zipfile.ZipFile(wheel) as archive:
                contents = {name: archive.read(name) for name in archive.namelist()}
            contents["smsd/__init__.py"] += b"# changed without updating RECORD\n"
            with zipfile.ZipFile(wheel, "w") as archive:
                for name, data in contents.items():
                    archive.writestr(name, data)
            with self.assertRaisesRegex(ValueError, "Wrong RECORD hash or size"):
                COLLECTOR.check_wheel(wheel, "7.2.0", SOURCES)


if __name__ == "__main__":
    unittest.main()

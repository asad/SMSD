#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Check and collect matching CPython 3.14 wheels without publishing."""

import argparse
import ast
import base64
import csv
import hashlib
import importlib.util
import io
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tomllib
import zipfile
from email.parser import BytesParser


_BUILD_SPEC = importlib.util.spec_from_file_location(
    "smsd_release_build", Path(__file__).resolve().with_name("build_python_wheels.py"))
BUILD = importlib.util.module_from_spec(_BUILD_SPEC)
_BUILD_SPEC.loader.exec_module(BUILD)


def matches_windows_repaired_init(wheel, prefix, original, repaired):
    """Accept only delvewheel 1.13.1's known CPython 3.14 DLL bootstrap."""
    repair_metadata = prefix + "DELVEWHEEL"
    names = wheel.namelist()
    if (repair_metadata not in names
            or not wheel.read(repair_metadata).startswith(b"Version: 1.13.1\n")
            or not any(name.startswith("smsd.libs/") and name.count("/") == 1
                       and name.endswith(".dll") for name in names)):
        return False
    # This release uses LF source with a module docstring. Reconstruct the
    # repair's insertion instead of removing arbitrary marker-delimited code.
    if (b"\r" in original or not original.endswith(b"\n")
            or b"# start delvewheel patch" in original
            or b"# end delvewheel patch" in original):
        return False
    try:
        tree = ast.parse(original)
    except SyntaxError:
        return False
    if not ast.get_docstring(tree) or len(tree.body) < 2:
        return False
    patch = b"""# start delvewheel patch
def _delvewheel_patch_1_13_1():
    import os
    if os.path.isdir(libs_dir := os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'smsd.libs'))):
        os.add_dll_directory(libs_dir)


_delvewheel_patch_1_13_1()
del _delvewheel_patch_1_13_1
# end delvewheel patch
"""
    lines = original.splitlines(keepends=True)
    future_imports = [node for node in tree.body
                      if isinstance(node, ast.ImportFrom) and node.module == "__future__"]
    if future_imports:
        boundary = sum(map(len, lines[:future_imports[-1].end_lineno]))
    else:
        docstring = tree.body[0].value
        boundary = (sum(map(len, lines[:docstring.end_lineno - 1]))
                    + docstring.end_col_offset)
        # Upstream requires the closing docstring line to contain only whitespace.
        trailing = original[boundary:].split(b"\n", 1)[0]
        if trailing and not trailing.isspace():
            return False
    expected = (original[:boundary].rstrip() + b"\n\n\n" + patch
                + b"\n" + original[boundary:].lstrip())
    return repaired == expected


def check_wheel(path, version, sources):
    fields = path.name[:-4].split("-")
    if len(fields) != 5 or fields[:4] != ["smsd", version, "cp314", "cp314"]:
        raise ValueError(f"Unexpected package, version or interpreter: {path.name}")
    platform = fields[4]
    if platform == "win_amd64":
        family = "windows"
    elif platform == "macosx_26_0_arm64":
        family = "macos"
    elif ("manylinux_2_28_x86_64" in platform.split(".")
          and all(re.fullmatch(r"manylinux_2_\d+_x86_64", tag)
                  for tag in platform.split("."))):
        family = "linux"
    else:
        raise ValueError(f"Unexpected release platform: {platform}")
    with zipfile.ZipFile(path) as wheel:
        names = wheel.namelist()
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate wheel paths: {path.name}")
        for name in names:
            if name.startswith("/") or ".." in PurePosixPath(name).parts:
                raise ValueError(f"Unsafe wheel member: {name}")
        extensions = [name for name in names if name.startswith("smsd/_smsd.")
                      and name.endswith((".so", ".pyd")) and name.count("/") == 1]
        if len(extensions) != 1:
            raise ValueError(f"Expected one native extension: {path.name}")
        extension = extensions[0]
        binary = wheel.read(extension)
        if family == "windows":
            expected_extension = "smsd/_smsd.cp314-win_amd64.pyd"
            pe_offset = int.from_bytes(binary[60:64], "little")
            valid_binary = (binary[:2] == b"MZ" and pe_offset >= 64
                            and binary[pe_offset:pe_offset + 4] == b"PE\0\0"
                            and binary[pe_offset + 4:pe_offset + 6] == b"\x64\x86")
        elif family == "linux":
            expected_extension = "smsd/_smsd.cpython-314-x86_64-linux-gnu.so"
            valid_binary = (binary[:6] == b"\x7fELF\x02\x01"
                            and binary[18:20] == b"\x3e\0")
        else:
            expected_extension = "smsd/_smsd.cpython-314-darwin.so"
            valid_binary = (binary[:4] == b"\xcf\xfa\xed\xfe"
                            and binary[4:8] == b"\x0c\0\0\x01")
        if extension != expected_extension or not valid_binary:
            raise ValueError(f"Native extension does not match the target: {path.name}")
        prefix = f"smsd-{version}.dist-info/"
        metadata = BytesParser().parsebytes(wheel.read(prefix + "METADATA"))
        if metadata["Name"] != "smsd" or metadata["Version"] != version:
            raise ValueError(f"Wrong wheel metadata: {path.name}")
        wheel_metadata = BytesParser().parsebytes(wheel.read(prefix + "WHEEL"))
        expected_tags = {f"cp314-cp314-{tag}" for tag in platform.split(".")}
        if set(wheel_metadata.get_all("Tag", [])) != expected_tags:
            raise ValueError(f"Filename and wheel tags disagree: {path.name}")
        record_name = prefix + "RECORD"
        records = list(csv.reader(io.StringIO(wheel.read(record_name).decode("utf-8"))))
        if len({row[0] for row in records}) != len(records):
            raise ValueError(f"Duplicate RECORD entries: {path.name}")
        if {row[0] for row in records} != {name for name in names if not name.endswith("/")}:
            raise ValueError(f"Incomplete wheel RECORD: {path.name}")
        for name, checksum, size in records:
            if name == record_name:
                continue
            data = wheel.read(name)
            expected = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
            if checksum != "sha256=" + expected or size != str(len(data)):
                raise ValueError(f"Wrong RECORD hash or size: {path.name}: {name}")
        for name, data in sources.items():
            if name.startswith("python/smsd/") and name.endswith(".py"):
                suffix = name[len("python/"):]
            elif name.startswith("cpp/include/smsd/") and name.endswith(".hpp"):
                suffix = name[len("cpp/"):]
            else:
                continue
            matches = [entry for entry in names if entry == suffix or entry.endswith("/" + suffix)]
            if len(matches) != 1:
                raise ValueError(f"Wheel and source differ: {path.name}: {name}")
            installed = wheel.read(matches[0])
            if installed == data:
                continue
            if (family == "windows" and suffix == "smsd/__init__.py" and matches[0] == suffix
                    and matches_windows_repaired_init(wheel, prefix, data, installed)):
                continue
            raise ValueError(f"Wheel and source differ: {path.name}: {name}")
        licenses = ["LICENSE", "NOTICE", "licenses/libomp/LICENSE.TXT",
                    "licenses/libgomp/COPYING3", "licenses/libgomp/COPYING.RUNTIME"]
        licenses.extend(name for name in ("licenses/msvc/README.md",
                                          "licenses/msvc/LICENSE-2022.docx",
                                          "licenses/msvc/LICENSE-2026.docx")
                        if name in sources)
        for name in licenses:
            if wheel.read(prefix + "licenses/" + name) != sources[name]:
                raise ValueError(f"Wrong licence copy: {path.name}: {name}")
    return family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--wheel-dir", type=Path, action="append", default=[])
    parser.add_argument("--source-manifest", type=Path)
    parser.add_argument("--provenance", type=Path, action="append", default=[],
                        help="Checked build receipt; adjacent wheel-build-provenance.json files are found automatically")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="Collect local build results before all three platforms are available")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    directory = args.release_dir.resolve()
    archives = list(directory.glob("smsd-*.tar.gz"))
    archives = [p for p in archives if p.name.count("-") == 1]
    if len(archives) != 1:
        raise ValueError("Expected one Python source distribution in the release directory")
    version = archives[0].name[len("smsd-"):-len(".tar.gz")]
    source = BUILD.checked_source(archives[0], args.source_manifest)
    sources = BUILD.read_source_archive(archives[0])
    project = tomllib.loads(sources["pyproject.toml"].decode("utf-8"))["project"]
    if project["name"] != "smsd" or project["version"] != version:
        raise ValueError("Source filename and package metadata disagree")
    checkout = Path(__file__).resolve().parents[1]
    required = {"python/smsd/__init__.py", "python/smsd/search.py",
                "cpp/CMakeLists.txt", "cpp/bindings/pybind11/smsd_bindings.cpp"}
    required.update(p.relative_to(checkout).as_posix()
                    for p in (checkout / "cpp/include/smsd").rglob("*.hpp"))
    if "cpp/include/smsd/smsd.hpp" not in required:
        raise ValueError("Run collection from a complete SMSD source checkout")
    for name in required:
        if name not in sources or sources[name] != (checkout / name).read_bytes():
            raise ValueError(f"Source distribution and release checkout differ: {name}")
    candidates = list(directory.glob("*.whl"))
    for wheel_dir in args.wheel_dir:
        candidates.extend(wheel_dir.rglob("*.whl"))
    selected = {}
    for wheel in sorted({p.resolve() for p in candidates}):
        family = check_wheel(wheel, version, sources)
        previous = selected.get(family)
        if previous and BUILD.sha256(previous) != BUILD.sha256(wheel):
            raise ValueError(f"Multiple different {family} wheels; select one validated build")
        selected[family] = wheel
    missing = {"linux", "macos", "windows"} - set(selected)
    if missing and not args.allow_incomplete:
        raise ValueError("Missing release wheels: " + ", ".join(sorted(missing)))
    manifest_path = directory / "WHEEL_MANIFEST.json"
    if args.check_only:
        BUILD.checked_wheel_manifest(directory, source, selected, BUILD.read_manifest(manifest_path))
    else:
        receipts = {path.resolve() for path in args.provenance}
        receipts.update(path.parent / "wheel-build-provenance.json" for path in candidates
                        if (path.parent / "wheel-build-provenance.json").is_file())
        proofs = [BUILD.read_manifest(path) for path in sorted(receipts)]
        if manifest_path.exists():
            existing = BUILD.read_manifest(manifest_path)
            recorded_entries = existing.get("wheels")
            if (isinstance(recorded_entries, list)
                    and {entry.get("platform") for entry in recorded_entries if isinstance(entry, dict)} != set(selected)):
                raise ValueError("Existing wheel manifest has different platforms; collect the full set in a new release directory")
            BUILD.checked_wheel_manifest(directory, source, selected, existing)
            proofs.extend(entry["build_provenance"] for entry in existing["wheels"])
        entries = []
        for family, wheel in sorted(selected.items()):
            matching = [proof for proof in proofs if proof.get("wheel") == wheel.name]
            if not matching:
                raise ValueError("Missing checked-build provenance: " + wheel.name)
            entry = BUILD.wheel_entry(wheel, family, source, matching[0])
            if any(BUILD.wheel_entry(wheel, family, source, proof) != entry for proof in matching[1:]):
                raise ValueError("Conflicting checked-build provenance: " + wheel.name)
            entries.append(entry)
        for path in args.provenance:
            proof = BUILD.read_manifest(path)
            if proof.get("wheel") not in {wheel.name for wheel in selected.values()}:
                raise ValueError("Explicit provenance does not identify a selected wheel")
        manifest = {"schema_version": 1,
                    **{field: source[field] for field in ("version", "source_commit", "sdist", "sdist_sha256")},
                    "complete": not missing, "missing_platforms": sorted(missing), "wheels": entries}
        # Check both destinations before copying any release assets.
        for path, document in ((directory / "SOURCE_MANIFEST.json", source), (manifest_path, manifest)):
            if path.is_symlink() or (path.exists() and BUILD.read_manifest(path) != document):
                raise ValueError("Existing manifest differs: " + path.name)
        for wheel in selected.values():
            destination = directory / wheel.name
            if wheel != destination:
                shutil.copy2(wheel, destination)
        BUILD.write_manifest(directory / "SOURCE_MANIFEST.json", source)
        BUILD.write_manifest(manifest_path, manifest)
        files = sorted(p for p in directory.iterdir() if p.is_file() and p.name != "SHA256SUMS")
        (directory / "SHA256SUMS").write_text(
            "".join(f"{BUILD.sha256(p)}  {p.name}\n" for p in files), encoding="utf-8")
    print("Checked release wheels:", ", ".join(sorted(selected)) or "none")
    if missing:
        print("Still required before publication:", ", ".join(sorted(missing)))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, zipfile.BadZipFile, tarfile.TarError) as error:
        raise SystemExit(str(error)) from error

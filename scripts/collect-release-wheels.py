#!/usr/bin/env python3
"""Check and collect matching CPython 3.14 wheels without publishing."""

import argparse
import base64
import csv
import hashlib
import io
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
import tomllib
import zipfile
from email.parser import BytesParser


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_files(path):
    with tarfile.open(path, "r:gz") as archive:
        files = {}
        roots = set()
        for member in archive.getmembers():
            parts = PurePosixPath(member.name).parts
            if (not parts or member.issym() or member.islnk()
                    or ".." in parts or member.name.startswith("/")):
                raise ValueError(f"Unsafe source member: {member.name}")
            roots.add(parts[0])
            if member.isfile():
                name = "/".join(parts[1:])
                if name in files:
                    raise ValueError(f"Duplicate source member: {member.name}")
                files[name] = archive.extractfile(member).read()
        if len(roots) != 1:
            raise ValueError("Expected one source archive root")
        return files


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
            if len(matches) != 1 or wheel.read(matches[0]) != data:
                raise ValueError(f"Wheel and source differ: {path.name}: {name}")
        licenses = ["LICENSE", "NOTICE", "licenses/libomp/LICENSE.TXT",
                    "licenses/libgomp/COPYING3", "licenses/libgomp/COPYING.RUNTIME"]
        for name in licenses:
            if wheel.read(prefix + "licenses/" + name) != sources[name]:
                raise ValueError(f"Wrong license copy: {path.name}: {name}")
    return family


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--wheel-dir", type=Path, action="append", default=[])
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
    sources = source_files(archives[0])
    project = tomllib.loads(sources["pyproject.toml"].decode("utf-8"))["project"]
    if project["name"] != "smsd" or project["version"] != version:
        raise ValueError("Source filename and package metadata disagree")
    checkout = Path(__file__).resolve().parents[1]
    required = {"python/smsd/__init__.py", "python/smsd/mcs_engine.py",
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
        if previous and digest(previous) != digest(wheel):
            raise ValueError(f"Multiple different {family} wheels; select one validated build")
        selected[family] = wheel
    missing = {"linux", "macos", "windows"} - set(selected)
    if missing and not args.allow_incomplete:
        raise ValueError("Missing release wheels: " + ", ".join(sorted(missing)))
    if not args.check_only:
        for wheel in selected.values():
            destination = directory / wheel.name
            if wheel != destination:
                shutil.copy2(wheel, destination)
        files = sorted(p for p in directory.iterdir() if p.is_file() and p.name != "SHA256SUMS")
        (directory / "SHA256SUMS").write_text(
            "".join(f"{digest(p)}  {p.name}\n" for p in files), encoding="utf-8")
    print("Checked release wheels:", ", ".join(sorted(selected)) or "none")
    if missing:
        print("Still required before publication:", ", ".join(sorted(missing)))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, zipfile.BadZipFile, tarfile.TarError) as error:
        raise SystemExit(str(error)) from error

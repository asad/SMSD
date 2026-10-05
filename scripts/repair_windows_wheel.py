#!/usr/bin/env python3
"""Repair an AMD64 wheel using validated Microsoft release runtimes.

Run on native 64-bit Windows with delvewheel==1.13.1 installed. Runtime files
come from installed Visual Studio redistributables, or Windows System32; the
ambient PATH is never used to choose msvcp140.dll or vcomp140.dll.
"""

import argparse
import ctypes
from dataclasses import dataclass
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile
import zipfile


RUNTIME_NAMES = ("msvcp140.dll", "vcomp140.dll")
DELVEWHEEL_VERSION = "1.13.1"
VC_RUNTIME_PATTERN = re.compile(r"(?:msvcp|msvcr|vcruntime|vcomp|concrt|vccorlib|mfc)\d+\w*(?:-[0-9a-f]{32})?\.dll", re.I)


@dataclass(frozen=True)
class DllInfo:
    version: tuple
    company: str
    flags: int = 0


def pe_linker_version(data):
    """Read the architecture and MSVC linker family from a PE32+ image."""
    if len(data) < 64 or data[:2] != b"MZ":
        raise ValueError("Not a PE image")
    offset = struct.unpack_from("<I", data, 60)[0]
    if offset + 28 > len(data) or data[offset:offset + 4] != b"PE\0\0":
        raise ValueError("Invalid PE header")
    machine = struct.unpack_from("<H", data, offset + 4)[0]
    magic = struct.unpack_from("<H", data, offset + 24)[0]
    if machine != 0x8664 or magic != 0x20B:
        raise ValueError("Expected an AMD64 PE32+ image")
    return tuple(data[offset + 26:offset + 28])


def wheel_runtime_baseline(wheel):
    if not wheel.name.endswith("-win_amd64.whl"):
        raise ValueError("Expected a win_amd64 wheel")
    with zipfile.ZipFile(wheel) as archive:
        extensions = [name for name in archive.namelist() if name.endswith(".pyd")]
        if not extensions:
            raise ValueError("Wheel contains no native extension")
        versions = [pe_linker_version(archive.read(name)) for name in extensions]
        if any(version[0] != 14 or version[1] == 0 for version in versions):
            raise ValueError("Expected an MSVC v14 linker version in every extension")
        if any(re.search(r"/(msvcp140|vcomp140)(?:-[^/]+)?\.dll$", name, re.I)
               for name in archive.namelist()):
            raise ValueError("Repair requires a wheel without previously bundled Microsoft runtimes")
    # Microsoft's deployment rule matches the runtime major version and requires
    # the same or a newer minor version. Compiler/toolset and redistributable
    # patch/build numbers are distinct, so do not compare their full tuples.
    # https://learn.microsoft.com/en-us/cpp/windows/determining-which-dlls-to-redistribute
    return max(versions)


def read_dll_info(path):
    """Use Windows version resources rather than directory names or filenames."""
    from ctypes import wintypes

    version_api = ctypes.WinDLL("version", use_last_error=True)
    version_api.GetFileVersionInfoSizeW.argtypes = [wintypes.LPCWSTR, ctypes.POINTER(wintypes.DWORD)]
    version_api.GetFileVersionInfoSizeW.restype = wintypes.DWORD
    version_api.GetFileVersionInfoW.argtypes = [wintypes.LPCWSTR, wintypes.DWORD,
                                              wintypes.DWORD, wintypes.LPVOID]
    version_api.GetFileVersionInfoW.restype = wintypes.BOOL
    version_api.VerQueryValueW.argtypes = [wintypes.LPCVOID, wintypes.LPCWSTR,
                                         ctypes.POINTER(ctypes.c_void_p),
                                         ctypes.POINTER(wintypes.UINT)]
    version_api.VerQueryValueW.restype = wintypes.BOOL
    ignored = wintypes.DWORD()
    size = version_api.GetFileVersionInfoSizeW(str(path), ctypes.byref(ignored))
    if not size:
        raise ValueError("DLL has no readable version resource: " + str(path))
    buffer = ctypes.create_string_buffer(size)
    if not version_api.GetFileVersionInfoW(str(path), 0, size, buffer):
        raise ctypes.WinError(ctypes.get_last_error())

    def query(key):
        pointer, length = ctypes.c_void_p(), wintypes.UINT()
        if not version_api.VerQueryValueW(buffer, key, ctypes.byref(pointer), ctypes.byref(length)):
            raise ValueError("Missing DLL version resource " + key + ": " + str(path))
        return pointer, length.value

    pointer, length = query("\\")
    if length < 52:
        raise ValueError("Invalid fixed DLL version resource: " + str(path))
    fixed = struct.unpack("<13I", ctypes.string_at(pointer, 52))
    if fixed[0] != 0xFEEF04BD:
        raise ValueError("Invalid DLL version resource signature: " + str(path))
    version = (fixed[2] >> 16, fixed[2] & 0xFFFF, fixed[3] >> 16, fixed[3] & 0xFFFF)
    pointer, length = query("\\VarFileInfo\\Translation")
    translations = ctypes.string_at(pointer, length)
    companies = []
    for offset in range(0, length - 3, 4):
        language, codepage = struct.unpack_from("<HH", translations, offset)
        try:
            pointer, characters = query(f"\\StringFileInfo\\{language:04x}{codepage:04x}\\CompanyName")
            companies.append(ctypes.wstring_at(pointer, characters).rstrip("\0").strip())
        except ValueError:
            continue
    company = "Microsoft Corporation" if "Microsoft Corporation" in companies else "; ".join(companies)
    return DllInfo(version, company, fixed[7] & fixed[6])


def visual_studio_installations(environment):
    installations = []
    if environment.get("VSINSTALLDIR"):
        installations.append(Path(environment["VSINSTALLDIR"]))
    if environment.get("VCINSTALLDIR"):
        installations.append(Path(environment["VCINSTALLDIR"]).parent)
    if environment.get("VCToolsInstallDir"):
        installations.append(Path(environment["VCToolsInstallDir"]).parents[3])
    program_files = environment.get("ProgramFiles(x86)") or environment.get("ProgramFiles")
    if program_files:
        vswhere = Path(program_files) / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
        if vswhere.is_file():
            output = subprocess.check_output([str(vswhere), "-all", "-products", "*", "-requires",
                                              "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
                                              "-property", "installationPath", "-utf8"],
                                             encoding="utf-8-sig")
            installations.extend(Path(line.strip()) for line in output.splitlines() if line.strip())
    return list(dict.fromkeys(path.resolve() for path in installations))


def runtime_directories(environment, installations, system32):
    """Return only VS x64 CRT/OpenMP redist folders and the OS runtime folder."""
    roots = []
    if environment.get("VCToolsRedistDir"):
        roots.append(Path(environment["VCToolsRedistDir"]))
    for installation in installations:
        base = installation / "VC" / "Redist" / "MSVC"
        if base.is_dir():
            roots.extend(sorted(base.iterdir(), key=lambda path: path.name, reverse=True))
    directories = []
    for root in roots:
        for suffix in ("CRT", "OpenMP"):
            directories.extend((path.resolve(), "Visual Studio redist")
                               for path in sorted(root.glob(f"x64/Microsoft.VC*.{suffix}"))
                               if path.is_dir())
    directories.append((system32.resolve(), "Windows System32"))
    return list(dict.fromkeys(directories))


def system_directory():
    buffer = ctypes.create_unicode_buffer(32768)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GetSystemDirectoryW.argtypes = [ctypes.c_wchar_p, ctypes.c_uint]
    kernel32.GetSystemDirectoryW.restype = ctypes.c_uint
    size = kernel32.GetSystemDirectoryW(buffer, len(buffer))
    if not size or size >= len(buffer):
        raise ctypes.WinError(ctypes.get_last_error())
    return Path(buffer.value)


def select_runtimes(directories, baseline, metadata_reader=read_dll_info):
    selected = {}
    for name in RUNTIME_NAMES:
        candidates, rejected = [], []
        for directory, source in directories:
            path = directory / name
            if not path.is_file():
                continue
            try:
                pe_linker_version(path.read_bytes())
                info = metadata_reader(path)
                validate_runtime(info, baseline)
                candidates.append((source == "Visual Studio redist", info.version,
                                   str(path).casefold(), path, info, source))
            except (OSError, ValueError) as error:
                rejected.append(str(path) + ": " + str(error))
        if not candidates:
            detail = "\n".join(rejected) or "no DLL found in Visual Studio redist or Windows System32"
            raise ValueError("No supported AMD64 Microsoft " + name + ":\n" + detail)
        _, _, _, path, info, source = max(candidates, key=lambda candidate: candidate[:3])
        selected[name] = (path, info, source)
    return selected


def validate_runtime(info, baseline):
    if info.company != "Microsoft Corporation":
        raise ValueError("not a Microsoft runtime")
    if info.flags & 0x3:  # VS_FF_DEBUG or VS_FF_PRERELEASE
        raise ValueError("debug or prerelease runtime")
    if info.version[0] != baseline[0] or info.version[:2] < baseline:
        raise ValueError("runtime " + ".".join(map(str, info.version))
                         + " does not support required family " + ".".join(map(str, baseline)))


def verify_bundled_runtimes(wheel, original_bytes, baseline, metadata_reader=read_dll_info):
    with zipfile.ZipFile(wheel) as archive, tempfile.TemporaryDirectory(prefix="smsd-runtime-check-") as temporary:
        matched = set()
        for name, expected in original_bytes.items():
            stem = name[:-4]
            matches = [member for member in archive.namelist()
                       if re.fullmatch(r"smsd\.libs/" + stem + r"-[0-9a-f]{32}\.dll", member, re.I)]
            if len(matches) != 1 or archive.read(matches[0]) != expected:
                raise ValueError("Repaired wheel must contain the selected, unmodified " + name)
            matched.add(matches[0])
        # A future transitive VC dependency must obey the same release/version
        # contract. Other third-party DLLs remain handled by delvewheel normally.
        for member in archive.namelist():
            filename = member.rsplit("/", 1)[-1]
            if member not in matched and VC_RUNTIME_PATTERN.fullmatch(filename):
                data = archive.read(member)
                pe_linker_version(data)
                path = Path(temporary) / filename
                path.write_bytes(data)
                validate_runtime(metadata_reader(path), baseline)


def repair(wheel, destination):
    if os.name != "nt" or sys.maxsize <= 2**32:
        raise ValueError("Windows wheel repair requires native 64-bit Windows Python")
    if importlib.metadata.version("delvewheel") != DELVEWHEEL_VERSION:
        raise ValueError("Install delvewheel==" + DELVEWHEEL_VERSION + " before Windows repair")
    baseline = wheel_runtime_baseline(wheel)
    directories = runtime_directories(os.environ, visual_studio_installations(os.environ), system_directory())
    selected = select_runtimes(directories, baseline)
    destination.mkdir(parents=True, exist_ok=True)
    output = destination / wheel.name
    if output.exists():
        raise ValueError("Repaired output wheel already exists: " + str(output))
    original_bytes = {name: path.read_bytes() for name, (path, _, _) in selected.items()}
    print(json.dumps({"minimum_runtime_family": ".".join(map(str, baseline)),
                      "delvewheel": DELVEWHEEL_VERSION,
                      "runtimes": [{"name": name, "source": source, "path": str(path),
                                    "version": ".".join(map(str, info.version)),
                                    "sha256": hashlib.sha256(original_bytes[name]).hexdigest()}
                                   for name, (path, info, source) in selected.items()]}, indent=2), flush=True)
    with tempfile.TemporaryDirectory(prefix="smsd-msvc-redist-") as temporary:
        # A single directory prevents a rejected copy of one DLL in the selected
        # other DLL's source directory from winning delvewheel's search order.
        for name, data in original_bytes.items():
            (Path(temporary) / name).write_bytes(data)
        # delvewheel's documented --add-path directories precede ambient PATH.
        # https://github.com/adang1345/delvewheel#additional-options
        subprocess.run([sys.executable, "-m", "delvewheel", "repair", "--add-path", temporary,
                        "-w", str(destination), "-v", str(wheel)], check=True)
    verify_bundled_runtimes(output, original_bytes, baseline)
    print("Verified both bundled Microsoft runtimes match the selected release DLLs.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--dest-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        repair(args.wheel.resolve(), args.dest_dir.resolve())
    except (OSError, ValueError, importlib.metadata.PackageNotFoundError,
            subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

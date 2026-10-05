#!/usr/bin/env python3
"""Prepare an sdist or build repaired CPython 3.14 CPU wheels from one sdist.

Install build/twine for preparation and cibuildwheel==4.2.1 for wheel builds.
Linux builds use local Docker or Podman; Windows/macOS require that host OS.
"""

import argparse
from email.parser import BytesParser
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import sys
import tarfile
import tomllib
import zipfile


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect_source(path):
    with tarfile.open(path, "r:gz") as archive:
        files = {}
        roots = set()
        for member in archive.getmembers():
            name = PurePosixPath(member.name)
            if (not name.parts or name.is_absolute() or ".." in name.parts or "\\" in member.name
                    or ":" in member.name or not (member.isdir() or member.isfile())):
                raise ValueError("Unsafe source archive member: " + member.name)
            roots.add(name.parts[0])
            if member.isfile():
                relative = str(PurePosixPath(*name.parts[1:]))
                if relative in files:
                    raise ValueError("Duplicate source archive member: " + relative)
                files[relative] = archive.extractfile(member).read()
        if len(roots) != 1:
            raise ValueError("Expected one source archive root")
    required = {"pyproject.toml", "python/smsd/__init__.py", "cpp/CMakeLists.txt",
                "cpp/include/smsd/mcs.hpp", "cpp/bindings/pybind11/smsd_bindings.cpp",
                "python/tests/test_smsd.py", "scripts/check_python_wheel.py",
                "scripts/check_native_wheel_build.py", "scripts/build_python_wheels.py",
                "scripts/cibuildwheel.toml", "LICENSE", "NOTICE"}
    missing = sorted(required - files.keys())
    if missing:
        raise ValueError("Incomplete source distribution: " + ", ".join(missing))
    version = tomllib.loads(files["pyproject.toml"].decode())["project"]["version"]
    python_version = re.search(rb'^__version__ = "([^"]+)"', files["python/smsd/__init__.py"], re.M)
    cpp_version = re.search(rb'project\(smsd VERSION ([^ ]+)', files["cpp/CMakeLists.txt"])
    if not python_version or not cpp_version or python_version[1].decode() != version or cpp_version[1].decode() != version:
        raise ValueError("Python, CMake and package versions differ")
    return {"version": version, "sdist": path.name, "sdist_sha256": sha256(path),
            "source_file_count": len(files)}


def prepare(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if list(args.output_dir.glob("*.tar.gz")):
        raise ValueError("Prepare output directory already contains a source distribution")
    command = [sys.executable, "-m", "build", "--sdist", "--outdir", str(args.output_dir)]
    if args.no_isolation:
        command.append("--no-isolation")
    subprocess.run(command, check=True)
    archives = list(args.output_dir.glob("smsd-*.tar.gz"))
    if len(archives) != 1:
        raise ValueError("Expected exactly one SMSD source distribution")
    manifest = inspect_source(archives[0])
    if args.release_ref:
        if args.release_ref.startswith("v") and args.release_ref != "v" + manifest["version"]:
            raise ValueError("Release tag and source version differ")
        resolved = subprocess.check_output(["git", "rev-parse", args.release_ref + "^{commit}"], text=True).strip()
    else:
        resolved = None
    manifest["release_ref"] = args.release_ref
    manifest["source_commit"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    manifest["working_tree_clean"] = not bool(subprocess.check_output(["git", "status", "--porcelain"], text=True).strip())
    if resolved and resolved != manifest["source_commit"]:
        raise ValueError("Source checkout differs from the requested release ref")
    if args.release_ref and not manifest["working_tree_clean"]:
        raise ValueError("Release-ref preparation requires a clean source checkout")
    (args.output_dir / "source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def build(args):
    source = inspect_source(args.sdist)
    if args.source_manifest:
        recorded = json.loads(args.source_manifest.read_text())
        if recorded["sdist_sha256"] != source["sdist_sha256"] or recorded["version"] != source["version"]:
            raise ValueError("Source distribution does not match its manifest")
        source.update(recorded)
    architectures = {"linux": ("x86_64", "aarch64"), "macos": ("arm64", "x86_64"), "windows": ("AMD64",)}
    arch = args.arch or architectures[args.platform][0]
    if arch not in architectures[args.platform]:
        raise ValueError("Unsupported architecture for " + args.platform)
    if not args.dry_run and args.platform in ("windows", "macos"):
        expected_host = "Windows" if args.platform == "windows" else "Darwin"
        if platform.system() != expected_host:
            raise ValueError(args.platform + " builds require a native " + expected_host + " host")
    tag = {"linux": "manylinux_" + arch, "macos": "macosx_" + arch, "windows": "win_amd64"}[args.platform]
    identifier = "cp314-" + tag
    env = os.environ.copy()
    # Keep caller cache/container settings, but enforce this release's build and test contract.
    for key in list(env):
        if key.startswith("CIBW_") and key not in ("CIBW_CACHE_PATH", "CIBW_CONTAINER_ENGINE"):
            del env[key]
    for key in ("CMAKE_ARGS", "SKBUILD_CMAKE_ARGS", "SKBUILD_CMAKE_DEFINE", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    if args.manylinux_image:
        if args.platform != "linux":
            raise ValueError("--manylinux-image is only valid for Linux")
        env["CIBW_MANYLINUX_" + arch.upper() + "_IMAGE"] = args.manylinux_image
    command = [sys.executable, "-m", "cibuildwheel", str(args.sdist.resolve()),
               "--only", identifier, "--config-file", "{package}/scripts/cibuildwheel.toml",
               "--output-dir", str(args.output_dir.resolve())]
    print(json.dumps({"command": command, "source": source,
                      "manylinux_image": args.manylinux_image or ("manylinux_2_28" if args.platform == "linux" else None)}, indent=2))
    if args.dry_run:
        return
    if importlib.metadata.version("cibuildwheel") != "4.2.1":
        raise ValueError("Use cibuildwheel==4.2.1 for this release build")
    if args.output_dir.exists() and list(args.output_dir.glob("*.whl")):
        raise ValueError("Wheel output directory already contains wheels")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, env=env, check=True)
    wheels = list(args.output_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError("Expected exactly one repaired, tested wheel")
    wheel = wheels[0]
    with zipfile.ZipFile(wheel) as archive:
        metadata = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata) != 1 or BytesParser().parsebytes(archive.read(metadata[0]))["Version"] != source["version"]:
            raise ValueError("Built wheel and source versions differ")
        if args.platform == "linux":
            for license_name in ("COPYING3", "COPYING.RUNTIME"):
                if not any(name.endswith("/" + license_name) for name in archive.namelist()):
                    raise ValueError("Linux wheel is missing libgomp license notice " + license_name)
    provenance = {**source, "build_identifier": identifier, "wheel": wheel.name,
                  "wheel_sha256": sha256(wheel), "cibuildwheel": "4.2.1",
                  "host_system": platform.system(), "host_architecture": platform.machine(),
                  "gpu_backends": "disabled", "openmp": "required by installed-wheel check",
                  "manylinux_image": args.manylinux_image or ("manylinux_2_28" if args.platform == "linux" else None),
                  "native_debug_tests": "passed", "installed_wheel_tests": "passed"}
    (args.output_dir / "wheel-build-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    source = commands.add_parser("prepare", help="Build and validate one source distribution")
    source.add_argument("--output-dir", type=Path, required=True)
    source.add_argument("--release-ref")
    source.add_argument("--no-isolation", action="store_true", help="Use already installed local build dependencies")
    source.set_defaults(function=prepare)
    wheel = commands.add_parser("build", help="Repair and test a wheel from an existing sdist")
    wheel.add_argument("--sdist", type=Path, required=True)
    wheel.add_argument("--source-manifest", type=Path)
    wheel.add_argument("--platform", choices=("linux", "macos", "windows"), required=True)
    wheel.add_argument("--arch")
    wheel.add_argument("--manylinux-image")
    wheel.add_argument("--output-dir", type=Path, required=True)
    wheel.add_argument("--dry-run", action="store_true")
    wheel.set_defaults(function=build)
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

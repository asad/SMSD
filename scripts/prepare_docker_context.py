#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Prepare Docker inputs from one source distribution and its checked wheels."""
import argparse
from email.parser import BytesParser
import importlib.util
import json
from pathlib import Path
import re
import shutil
import zipfile

_SPEC = importlib.util.spec_from_file_location("smsd_wheel_build", Path(__file__).with_name("build_python_wheels.py"))
BUILD = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BUILD)
RUNTIME_VERSION = "12.8.90"
RUNTIME_WHEEL = "nvidia_cuda_runtime_cu12-12.8.90-py3-none-manylinux2014_x86_64.manylinux_2_17_x86_64.whl"
RUNTIME_SHA256 = "adade8dcbd0edf427b7204d480d6066d33902cab2a4707dcfc48a2d0fd44ab90"


def wheel_metadata(path, package, version):
    with zipfile.ZipFile(path) as archive:
        names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(names) != 1:
            raise ValueError("Expected one wheel metadata file")
        metadata = BytesParser().parsebytes(archive.read(names[0]))
        if metadata["Name"].replace("_", "-") != package or metadata["Version"] != version:
            raise ValueError("Wheel metadata differs: " + path.name)


def checked_inputs(args):
    source_manifest = BUILD.read_manifest(args.release_dir / "SOURCE_MANIFEST.json")
    if Path(source_manifest["sdist"]).name != source_manifest["sdist"]:
        raise ValueError("Source manifest filename is not a basename")
    sdist = args.release_dir / source_manifest["sdist"]
    source = BUILD.checked_source(sdist)
    if source["version"] != "8.1.0":
        raise ValueError("Expected SMSD 8.1.0")
    manifest = BUILD.read_manifest(args.release_dir / "WHEEL_MANIFEST.json")
    selected = {}
    for entry in manifest["wheels"]:
        name = entry["wheel"]
        if Path(name).name != name:
            raise ValueError("Wheel manifest filename is not a basename")
        path = args.release_dir / name
        identifier = BUILD.wheel_target(name, source["version"])
        if identifier in selected:
            raise ValueError("Duplicate wheel target")
        selected[identifier] = path
    BUILD.checked_wheel_manifest(args.release_dir, source, selected, manifest)
    files = BUILD.read_source_archive(sdist)
    required = ("Dockerfile", "Dockerfile.cuda", ".dockerignore", "scripts/docker_install.py", "scripts/docker_smoke.py",
                "java/cli/pom.xml", "java/core/pom.xml", "LICENSE", "NOTICE")
    for name in required:
        if name not in files:
            raise ValueError("Docker source missing from frozen distribution: " + name)
    cores = {}
    for architecture in args.architecture or ("amd64", "arm64"):
        platform = "x86_64" if architecture == "amd64" else "aarch64"
        identifier = "cp314-manylinux_" + platform
        if identifier not in selected:
            raise ValueError("Checked CPython3.14 wheel missing: " + identifier)
        cores[architecture] = selected[identifier]
    cuda = []
    license_text = cub_text = cccl_text = None
    if args.cuda_wheel:
        if "amd64" not in cores or not args.cuda_provenance or not args.runtime_wheel:
            raise ValueError("CUDA needs the amd64 core, CUDA provenance and runtime wheel")
        proof = BUILD.read_manifest(args.cuda_provenance)
        expected = {"source_commit": source["source_commit"], "wheel": args.cuda_wheel.name,
                    "wheel_sha256": BUILD.sha256(args.cuda_wheel), "smsd_version": source["version"],
                    "installed_wheel_tests": "passed"}
        if any(proof.get(name) != value for name, value in expected.items()):
            raise ValueError("CUDA wheel provenance differs or lacks installed checks")
        if not re.fullmatch(r"smsd_cuda-" + re.escape(source["version"])
                            + r"-py3-none-(?:manylinux[^-]*x86_64)\.whl", args.cuda_wheel.name):
            raise ValueError("Expected a repaired Linuxx86_64 CUDA wheel")
        wheel_metadata(args.cuda_wheel, "smsd-cuda", source["version"])
        with zipfile.ZipFile(args.cuda_wheel) as archive:
            cub = sorted(name for name in archive.namelist() if name.endswith("/licenses/cub/LICENSE"))
            cub_text = archive.read(cub[0]) if cub else b""
            cccl = sorted(name for name in archive.namelist() if name.endswith("/licenses/cccl/LICENSE"))
            cccl_text = archive.read(cccl[0]) if cccl else b""
        if not cub_text:
            raise ValueError("CUB licence missing from the CUDA wheel")
        if not cccl_text:
            raise ValueError("CCCL licence missing from the CUDA wheel")
        if args.runtime_wheel.name != RUNTIME_WHEEL or BUILD.sha256(args.runtime_wheel) != RUNTIME_SHA256:
            raise ValueError("Expected the pinned CUDA12.8.90 runtime wheel")
        wheel_metadata(args.runtime_wheel, "nvidia-cuda-runtime-cu12", RUNTIME_VERSION)
        cuda = [(args.cuda_wheel, "smsd-cuda", source["version"]),
                (args.runtime_wheel, "nvidia-cuda-runtime-cu12", RUNTIME_VERSION)]
        with zipfile.ZipFile(args.runtime_wheel) as archive:
            notices = [(name, archive.read(name)) for name in archive.namelist()
                       if Path(name).name.lower() in ("license", "license.txt", "licence", "licence.txt")]
        if not notices:
            raise ValueError("NVIDIA runtime wheel licence missing")
        license_text = b"\n\n".join(name.encode() + b"\n\n" + text for name, text in notices) + b"\n"
    elif args.cuda_provenance or args.runtime_wheel:
        raise ValueError("CUDA inputs must be supplied together")
    return source, sdist, files, cores, cuda, license_text, cub_text, cccl_text


def prepare(args):
    source, sdist, files, cores, cuda, license_text, cub_text, cccl_text = checked_inputs(args)
    args.output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(sdist, args.output / "source.tar.gz")
    shutil.copyfile(args.release_dir / "SOURCE_MANIFEST.json", args.output / "SOURCE_MANIFEST.json")
    for name in ("Dockerfile", "Dockerfile.cuda", ".dockerignore", "scripts/docker_install.py", "scripts/docker_smoke.py", "LICENSE", "NOTICE"):
        destination = args.output / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(files[name])
    manifest = {"version": source["version"], "source_commit": source["source_commit"],
                "sdist_sha256": source["sdist_sha256"], "cpu": {}, "cuda": [], "licenses": {}}
    def copy_wheel(path, directory, package, version):
        destination = args.output / "wheels" / directory / path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
        return {"path": destination.relative_to(args.output).as_posix(), "sha256": BUILD.sha256(path), "package": package, "version": version}
    for architecture, path in cores.items():
        manifest["cpu"][architecture] = copy_wheel(path, architecture, "smsd", source["version"])
    for path, package, version in cuda:
        manifest["cuda"].append(copy_wheel(path, "cuda", package, version))
    if license_text:
        destination = args.output / "licenses/NVIDIA-CUDA.txt"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(license_text)
    for name, text in (("licenses/cub/LICENSE", cub_text), ("licenses/cccl/LICENSE", cccl_text)):
        if text:
            destination = args.output / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(text)
    for name in ("LICENSE", "NOTICE", "licenses/NVIDIA-CUDA.txt", "licenses/cub/LICENSE", "licenses/cccl/LICENSE"):
        if (args.output / name).exists():
            manifest["licenses"][name] = BUILD.sha256(args.output / name)
    (args.output / "ARTIFACTS.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Prepared Docker context for", source["source_commit"], "at", args.output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--architecture", choices=("amd64", "arm64"), action="append")
    parser.add_argument("--cuda-wheel", type=Path)
    parser.add_argument("--cuda-provenance", type=Path)
    parser.add_argument("--runtime-wheel", type=Path)
    prepare(parser.parse_args())


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError, zipfile.BadZipFile) as error:
        raise SystemExit(str(error)) from error

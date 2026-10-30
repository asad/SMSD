#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Collect a CUDA wheel and its actual checked-build proof without publishing."""
import argparse
import importlib.util
from pathlib import Path
import shutil
import subprocess


SPEC = importlib.util.spec_from_file_location("pypi_checks", Path(__file__).with_name("check_pypi_release.py"))
CHECKS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKS)
BUILD = CHECKS.BUILD


def collect(release_dir, wheel, provenance, checkout):
    source = BUILD.read_manifest(release_dir / "SOURCE_MANIFEST.json")
    BUILD.checked_source(release_dir / source["sdist"], release_dir / "SOURCE_MANIFEST.json")
    proof = BUILD.read_manifest(provenance)
    expected = {
        "schema_version": 1, "source_commit": source["source_commit"],
        "sdist_sha256": source["sdist_sha256"], "smsd_version": source["version"],
        "wheel": CHECKS.cuda_wheel_name(source["version"]), "wheel_sha256": BUILD.sha256(wheel),
        "installed_wheel_tests": "passed", "native_gpu_tests": "passed",
    }
    if set(proof) != set(expected) or any(type(proof[key]) is not type(value) or proof[key] != value
                                         for key, value in expected.items()):
        raise ValueError("CUDA proof differs from the source, wheel or successful checks")
    if wheel.name != expected["wheel"]:
        raise ValueError("CUDA wheel filename differs from its proof")
    files = {}
    listing = subprocess.check_output(["git", "-C", str(checkout), "ls-tree", "-r", "--name-only",
                                      source["source_commit"], "--", "packaging/smsd-cuda/"])
    for name in listing.decode().splitlines():
        files[name] = subprocess.check_output(["git", "-C", str(checkout), "show",
                                               source["source_commit"] + ":" + name])
    destination = release_dir / wheel.name
    manifest = release_dir / "CUDA_WHEEL_MANIFEST.json"
    for path in (destination, manifest, release_dir / "SHA256SUMS"):
        if path.is_symlink():
            raise ValueError("Release output must not be a symbolic link")
    if destination.exists() and BUILD.sha256(destination) != expected["wheel_sha256"]:
        raise ValueError("Existing CUDA wheel differs")
    if manifest.exists():
        existing = BUILD.read_manifest(manifest)
        if set(existing) != set(proof) or any(type(existing[key]) is not type(value)
                                             or existing[key] != value for key, value in proof.items()):
            raise ValueError("Existing CUDA manifest differs")
    # Validate in temporary staging so rejected bytes never enter the release directory.
    import tempfile
    with tempfile.TemporaryDirectory() as temporary:
        staging = Path(temporary)
        shutil.copyfile(wheel, staging / wheel.name)
        BUILD.write_manifest(staging / manifest.name, proof)
        CHECKS.check_cuda_wheel(staging, source["version"], source, files)
    if wheel.resolve() != destination.resolve():
        shutil.copyfile(wheel, destination)
    BUILD.write_manifest(manifest, proof)
    assets = sorted(path for path in release_dir.iterdir() if path.is_file() and path.name != "SHA256SUMS")
    (release_dir / "SHA256SUMS").write_text("".join(f"{BUILD.sha256(path)}  {path.name}\n" for path in assets))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--checkout", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    collect(args.release_dir, args.wheel, args.provenance, args.checkout)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error

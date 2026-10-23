#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Check and install the wheel inputs of an SMSD Docker image."""
import argparse
from email.parser import BytesParser
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import sysconfig
import zipfile


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected_files(root, manifest, architecture, gpu):
    if architecture not in ("amd64", "arm64") or (gpu and architecture != "amd64"):
        raise ValueError("Unsupported image architecture")
    entries = [manifest["cpu"][architecture]]
    if gpu:
        entries.extend(manifest["cuda"])
    expected = {"smsd": "8.0.2"}
    if gpu:
        expected.update({"smsd-cuda": "8.0.2", "nvidia-cuda-runtime-cu12": "12.8.90"})
    if (len(entries) != len(expected)
            or {entry["package"]: entry["version"] for entry in entries} != expected):
        raise ValueError("Expected the matching core, accelerator and runtime packages")
    for entry in entries:
        relative = Path(entry["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe wheel path")
        path = root / relative
        if entry["package"] == "smsd":
            platform = "x86_64" if architecture == "amd64" else "aarch64"
            prefix = "smsd-8.0.2-cp314-cp314-"
            tags = path.name[len(prefix):-4].split(".")
            if (not path.name.startswith(prefix) or not path.name.endswith(".whl")
                    or "manylinux_2_28_" + platform not in tags
                    or not all(re.fullmatch(r"manylinux_2_\d+_" + platform, tag) for tag in tags)):
                raise ValueError("Expected the matching CPython3.14 manylinux core wheel")
        if path.is_symlink() or sha256(path) != entry["sha256"]:
            raise ValueError("Wheel SHA256 differs: " + str(relative))
        with zipfile.ZipFile(path) as wheel:
            names = [name for name in wheel.namelist() if name.endswith(".dist-info/METADATA")]
            if len(names) != 1:
                raise ValueError("Expected one wheel metadata file")
            metadata = BytesParser().parsebytes(wheel.read(names[0]))
            if metadata["Name"].replace("_", "-") != entry["package"] or metadata["Version"] != entry["version"]:
                raise ValueError("Wheel metadata differs: " + str(relative))
        yield entry, path


def verify_installed(entries):
    for entry, path in entries:
        distribution = importlib.metadata.distribution(entry["package"])
        if distribution.version != entry["version"]:
            raise ValueError("Installed version differs")
        with zipfile.ZipFile(path) as wheel:
            for name in wheel.namelist():
                if name.endswith("/") or ".data/" in name or name.endswith(".dist-info/RECORD"):
                    continue
                installed = distribution.locate_file(name)
                if not installed.is_file() or installed.read_bytes() != wheel.read(name):
                    raise ValueError("Installed wheel bytes differ: " + name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--architecture", required=True)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    manifest = json.loads((args.root / "ARTIFACTS.json").read_text())
    if (manifest["version"] != "8.0.2" or not re.fullmatch(r"[0-9a-f]{40}", args.revision)
            or manifest["source_commit"] != args.revision):
        raise ValueError("Image version or source revision differs")
    entries = list(selected_files(args.root, manifest, args.architecture, args.gpu))
    if args.verify_only:
        source = json.loads((args.root / "SOURCE_MANIFEST.json").read_text())
        if (sha256(args.root / "source.tar.gz") != manifest["sdist_sha256"]
                or source["sdist_sha256"] != manifest["sdist_sha256"]
                or source["source_commit"] != args.revision or source["working_tree_clean"] is not True):
            raise ValueError("Frozen source identity differs")
        for name, digest in manifest.get("licenses", {}).items():
            if name == "licenses/NVIDIA-CUDA.txt" and not args.gpu:
                continue
            if sha256(args.root / name) != digest:
                raise ValueError("Licence bytes differ: " + name)
        return
    if sys.version_info[:2] != (3, 14):
        raise ValueError("The image needs CPython3.14")
    requirements = args.root / "requirements.txt"
    requirements.write_text("".join(f"{entry['package']}=={entry['version']} --hash=sha256:{entry['sha256']}\n"
                                    for entry, path in entries))
    links = sorted({str(path.parent) for entry, path in entries})
    command = [sys.executable, "-m", "pip", "install", "--no-cache-dir", "--no-index", "--no-deps",
               "--only-binary=:all:", "--require-hashes", "-r", str(requirements)]
    for directory in links:
        command.extend(("--find-links", directory))
    subprocess.run(command, check=True)
    subprocess.run([sys.executable, "-m", "pip", "check"], check=True)
    verify_installed(entries)
    for name in ("ssl", "sqlite3", "ctypes", "bz2", "lzma", "readline", "dbm.gnu", "dbm.ndbm", "compression.zstd", "uuid"):
        importlib.import_module(name)
    os.environ["SMSD_DISABLE_GPU"] = "1"
    import smsd
    if smsd.__version__ != manifest["version"]:
        raise ValueError("Installed SMSD version differs")
    destination = Path("/usr/share/doc/smsd")
    destination.mkdir(parents=True, exist_ok=True)
    python_license = Path(sysconfig.get_path("stdlib")) / "LICENSE.txt"
    if not python_license.is_file() or not python_license.stat().st_size:
        raise ValueError("The Python distribution licence and notices are missing")
    licence_destination = destination / "licenses/python/LICENSE.txt"
    licence_destination.parent.mkdir(parents=True, exist_ok=True)
    licence_destination.write_bytes(python_license.read_bytes())
    (destination / "ARTIFACTS.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

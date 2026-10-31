#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Validate published SMSD distributions before uploading the unchanged files."""

import argparse
import base64
import csv
from email.parser import BytesParser
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile


_BUILD_SPEC = importlib.util.spec_from_file_location(
    "smsd_release_build", Path(__file__).resolve().with_name("build_python_wheels.py"))
BUILD = importlib.util.module_from_spec(_BUILD_SPEC)
_BUILD_SPEC.loader.exec_module(BUILD)
tomllib = BUILD.tomllib


def require(condition, message):
    if not condition:
        raise ValueError(message)


def release_version(tag):
    match = re.fullmatch(r"v(\d+\.\d+\.\d+)", tag)
    require(match is not None, "Publishing requires a stable release tag")
    return match[1]


def safe_name(name):
    return (isinstance(name, str) and name not in ("", ".", "..")
            and not any(char in name for char in "/\\\r\n\0"))


def cuda_wheel_name(version):
    return f"smsd_cuda-{version}-py3-none-manylinux_2_34_x86_64.whl"


def dependency(value):
    match = re.fullmatch(r"([A-Za-z0-9_.-]+)\s*(.*)", value)
    require(match is not None, "Malformed CUDA dependency")
    specifier = match[2].strip()
    if specifier.startswith("(") or specifier.endswith(")"):
        require(specifier.startswith("(") and specifier.endswith(")"), "Malformed CUDA dependency parentheses")
        specifier = specifier[1:-1]
    constraints = specifier.replace(" ", "").split(",")
    require(all(re.fullmatch(r"(?:===|==|!=|~=|<=|>=|<|>)[^,;()\s]+", item) for item in constraints),
            "Unsupported CUDA dependency constraint")
    return re.sub(r"[-_.]+", "-", match[1]).lower(), tuple(sorted(constraints))


def selected_assets(metadata, tag, package="smsd"):
    require(package in ("smsd", "smsd-cuda"), "Unknown release package")
    version = release_version(tag)
    require(metadata.get("tag_name") == tag, "GitHub release tag differs")
    require(metadata.get("draft") is False, "GitHub release is still a draft")
    require(metadata.get("prerelease") is False, "A stable GitHub release is required")
    require(isinstance(metadata.get("published_at"), str) and metadata["published_at"],
            "GitHub release has not been published")
    assets = {}
    for asset in metadata.get("assets", []):
        name = asset.get("name")
        require(safe_name(name) and name not in assets, "Unsafe or duplicate release asset")
        assets[name] = asset
    selected = {f"smsd-{version}.tar.gz", "SOURCE_MANIFEST.json", "SHA256SUMS",
                "WHEEL_MANIFEST.json" if package == "smsd" else "CUDA_WHEEL_MANIFEST.json"}
    targets = {}
    for name in assets:
        if name.endswith(".whl"):
            if name.startswith("smsd_cuda-"):
                fields = name[:-4].split("-")
                require(len(fields) == 5 and fields[:4] == ["smsd_cuda", version, "py3", "none"]
                        and all(tag == "linux_x86_64" or re.fullmatch(r"manylinux_2_\d+_x86_64", tag)
                                for tag in fields[4].split(".")), "Unexpected CUDA release wheel: " + name)
                if package == "smsd-cuda":
                    require(name == cuda_wheel_name(version) and "cuda" not in targets,
                            "Expected one repaired Linux glibc2.34 x86_64 CUDA wheel")
                    targets["cuda"] = name
                continue
            identifier = BUILD.wheel_target(name, version)
            if package == "smsd":
                require(identifier not in targets, "Multiple " + identifier + " release wheels")
                targets[identifier] = name
        elif name.startswith(("smsd_cuda-", "smsd-cuda-")) and package == "smsd-cuda":
            raise ValueError("Unexpected CUDA distribution: " + name)
    missing = (BUILD.BUILD_IDENTIFIERS.keys() if package == "smsd" else {"cuda"}) - targets.keys()
    require(not missing, "Missing release wheels: " + ", ".join(sorted(missing)))
    selected.update(targets.values())
    require(selected <= assets.keys(), "Published release is missing required metadata or source")
    for name in selected:
        asset = assets[name]
        require(asset.get("state") == "uploaded", "Release asset is not uploaded: " + name)
        require(type(asset.get("size")) is int and asset["size"] > 0,
                "Invalid release asset size: " + name)
        require(isinstance(asset.get("digest"), str)
                and re.fullmatch(r"sha256:[0-9a-f]{64}", asset["digest"]),
                "Expected GitHub SHA256 digest: " + name)
    return version, {name: assets[name] for name in selected}, targets


def download(args):
    package = getattr(args, "package", "smsd")
    require(package in ("smsd", "smsd-cuda"), "Unknown release package")
    require(args.platform == ("all" if package == "smsd" else "linux"),
            "Publishing requires platform=all for smsd or platform=linux for smsd-cuda")
    release_version(args.tag)
    require(re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository),
            "Expected owner/repository")
    require(not args.release_dir.exists() or not any(args.release_dir.iterdir()),
            "Release download directory must be empty")
    raw = subprocess.check_output(["gh", "api",
                                   f"repos/{args.repository}/releases/tags/{args.tag}"])
    args.metadata_file.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_file.write_bytes(raw)
    _, assets, _ = selected_assets(BUILD.read_manifest(args.metadata_file), args.tag, package)
    args.release_dir.mkdir(parents=True, exist_ok=True)
    command = ["gh", "release", "download", args.tag, "--repo", args.repository,
               "--dir", str(args.release_dir)]
    for name in sorted(assets):
        command.extend(["--pattern", name])
    subprocess.run(command, check=True)
    print("Downloaded one published SDK and 16 release wheels with their manifests" if package == "smsd"
          else "Downloaded the published CUDA wheel with its source and build manifests")


def check_source(directory, version, manifest):
    name = f"smsd-{version}.tar.gz"
    files = BUILD.read_source_archive(directory / name)
    require(manifest.get("version") == version and manifest.get("sdist") == name,
            "Source manifest version or filename differs")
    require(manifest.get("sdist_sha256") == BUILD.sha256(directory / name),
            "Source manifest SHA256 differs")
    require(type(manifest.get("source_file_count")) is int
            and manifest["source_file_count"] == len(files), "Source manifest file count differs")
    require(manifest.get("working_tree_clean") is True, "Frozen checkout was not clean")
    for field in ("source_commit", "release_ref"):
        require(isinstance(manifest.get(field), str)
                and re.fullmatch(r"[0-9a-f]{40}", manifest[field]),
                "Expected full frozen Git identity: " + field)
    required = {"pyproject.toml", "pom.xml", "java/pom.xml", "PKG-INFO",
                "python/README.md", "python/smsd/__init__.py", "cpp/CMakeLists.txt",
                "scripts/collect-release-wheels.py", "LICENSE", "NOTICE"}
    require(required <= files.keys(), "SDK is missing required package files")
    project = tomllib.loads(files["pyproject.toml"].decode("utf-8"))["project"]
    package = BytesParser().parsebytes(files["PKG-INFO"])
    require(project["name"] == package["Name"] == "smsd"
            and project["version"] == package["Version"] == version,
            "SDK package metadata differs")
    require(project.get("requires-python") == package["Requires-Python"] == ">=3.12",
            "SDK Python requirement differs from >=3.12")
    metadata_parts = files["PKG-INFO"].split(b"\n\n", 1)
    require(len(metadata_parts) == 2 and metadata_parts[1] == files["python/README.md"],
            "SDK description differs from its Python README")
    for pom in ("pom.xml", "java/pom.xml"):
        require(ET.fromstring(files[pom]).findtext("{http://maven.apache.org/POM/4.0.0}version")
                == version, "Java package version differs")
    python_version = re.search(rb'^__version__ = "([^"]+)"', files["python/smsd/__init__.py"], re.M)
    cpp_version = re.search(rb"project\(smsd VERSION ([^ ]+)", files["cpp/CMakeLists.txt"])
    require(python_version and cpp_version and python_version[1].decode() == version
            and cpp_version[1].decode() == version, "Python or C++ version differs")
    return files


def git(checkout, *arguments):
    return subprocess.check_output(["git", "-C", str(checkout), *arguments]).strip()


def check_git(checkout, tag, manifest, files, package="smsd"):
    frozen = git(checkout, "rev-parse", "--verify", manifest["release_ref"] + "^{commit}").decode()
    require(frozen == manifest["source_commit"], "Frozen source ref resolves to another commit")
    tagged = git(checkout, "rev-parse", "--verify", f"refs/tags/{tag}" + "^{commit}").decode()
    require(git(checkout, "rev-parse", "HEAD").decode() == tagged,
            "Collector checkout differs from the selected tag")
    related = any(subprocess.run(
        ["git", "-C", str(checkout), "merge-base", "--is-ancestor", earlier, later],
        check=False).returncode == 0 for earlier, later in ((frozen, tagged), (tagged, frozen)))
    require(related, "Frozen source and release tag have unrelated histories")
    changed = git(checkout, "diff", "--name-only", "-z", frozen, tagged).split(b"\0")
    require(all(not name or Path(name.decode()).suffix in (".md", ".rst")
                or name == b"CITATION.cff" for name in changed),
            "Tagged production files differ from the frozen SDK")
    require(b"python/README.md" not in changed, "Tagged package description differs from SDK")
    if package == "smsd-cuda":
        require(b"packaging/smsd-cuda/README.md" not in changed, "Tagged CUDA package description differs from SDK")
    blobs = BUILD.tree_blobs(git(checkout, "ls-tree", "-r", "-z", "--full-tree", frozen))
    for name, content in files.items():
        if name == "PKG-INFO":
            continue
        require(blobs.get(name) == BUILD.git_blob_id(content), "SDK file differs from frozen source: " + name)
    collector = "scripts/collect-release-wheels.py"
    require((checkout / collector).read_bytes() == files[collector], "Tagged collector differs from SDK")
    return checkout / collector


def check_wheel_manifest(directory, version, source, targets, package_metadata):
    manifest = BUILD.read_manifest(directory / "WHEEL_MANIFEST.json")
    BUILD.checked_wheel_manifest(directory, source,
                                 {identifier: directory / name for identifier, name in targets.items()}, manifest)
    require(manifest.get("schema_version") == 2, "Unexpected wheel manifest schema")
    for field in ("version", "source_commit", "sdist", "sdist_sha256"):
        require(manifest.get(field) == source[field], "Wheel manifest identity differs: " + field)
    entries = manifest.get("wheels")
    require(isinstance(entries, list) and len(entries) == len(BUILD.BUILD_IDENTIFIERS), "Expected 16 wheel manifest entries")
    recorded = {}
    for entry in entries:
        identifier = entry.get("build_identifier")
        require(identifier in targets and identifier not in recorded, "Duplicate or unknown manifest target")
        require(entry.get("wheel") == targets[identifier], "Wheel manifest filename differs")
        require(entry.get("wheel_sha256") == BUILD.sha256(directory / targets[identifier]),
                "Wheel manifest SHA256 differs")
        with zipfile.ZipFile(directory / targets[identifier]) as wheel:
            require(wheel.read(f"smsd-{version}.dist-info/METADATA") == package_metadata,
                    "Wheel package metadata or long description differs from the final SDK")
        require(entry.get("native_debug_tests") == "passed"
                and entry.get("installed_wheel_tests") == "passed",
                "Wheel manifest lacks successful native and installed tests")
        recorded[identifier] = entry
    require(set(recorded) == set(BUILD.BUILD_IDENTIFIERS), "Wheel manifest is incomplete")


def check_hashes(directory, assets):
    require({path.name for path in directory.iterdir()} == set(assets),
            "Release input directory does not contain the exact selected files")
    sums = {}
    for line in (directory / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        row = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        require(row is not None and safe_name(row[2]) and row[2] not in sums,
                "Unsafe, malformed or duplicate release checksum")
        sums[row[2]] = row[1]
    for name, asset in assets.items():
        path = directory / name
        require(path.is_file() and not path.is_symlink(), "Expected a regular release file")
        actual = BUILD.sha256(path)
        require(path.stat().st_size == asset["size"] and asset["digest"] == "sha256:" + actual,
                "GitHub asset digest or size differs: " + name)
        if name != "SHA256SUMS":
            require(sums.get(name) == actual, "Release SHA256SUMS differs: " + name)


def check_cuda_source(checkout, source):
    root = "packaging/smsd-cuda/"
    blobs = BUILD.tree_blobs(git(checkout, "ls-tree", "-r", "-z", "--full-tree", source["source_commit"]))
    names = {root + name for name in ("pyproject.toml", "README.md", "LICENSE", "NOTICE",
                                     "python/smsd_cuda/__init__.py")}
    names.update(name for name in blobs if name.startswith(root + "python/") and name.endswith(".py"))
    require(names <= blobs.keys(), "Frozen commit is missing CUDA package inputs")
    files = {}
    for name in sorted(names):
        path = checkout / name
        require(path.is_file() and not path.is_symlink(), "Missing or unsafe CUDA package input: " + name)
        content = path.read_bytes()
        require(BUILD.git_blob_id(content) == blobs[name], "CUDA package input differs from frozen source: " + name)
        files[name] = content
    return files


def check_cuda_wheel(directory, version, source, files):
    name = cuda_wheel_name(version)
    expected = {"schema_version": 1, "source_commit": source["source_commit"],
                "sdist_sha256": source["sdist_sha256"], "smsd_version": version,
                "wheel": name, "wheel_sha256": BUILD.sha256(directory / name),
                "installed_wheel_tests": "passed", "native_gpu_tests": "passed"}
    proof = BUILD.read_manifest(directory / "CUDA_WHEEL_MANIFEST.json")
    require(set(proof) == set(expected) and all(type(proof[key]) is type(value) and proof[key] == value
            for key, value in expected.items()), "CUDA wheel proof differs or lacks successful actual checks")
    root = "packaging/smsd-cuda/"
    require({root + path for path in ("pyproject.toml", "README.md", "LICENSE", "NOTICE", "python/smsd_cuda/__init__.py")} <= files.keys(),
            "Frozen source is missing CUDA package inputs")
    project = tomllib.loads(files[root + "pyproject.toml"].decode())["project"]
    dependencies = {dependency(value) for value in project["dependencies"]}
    require(project["name"] == "smsd-cuda" and project["version"] == version
            and project["requires-python"] == ">=3.12" and ("smsd", (f"=={version}",)) in dependencies,
            "Frozen CUDA package metadata or paired core version differs")
    module = files[root + "python/smsd_cuda/__init__.py"]
    for field in ("__version__", "__smsd_version__"):
        assigned = re.search(rb'^' + field.encode() + rb' = "([^"]+)"', module, re.M)
        require(assigned and assigned[1].decode() == version, "Frozen CUDA module version differs")
    prefix = f"smsd_cuda-{version}.dist-info/"
    with zipfile.ZipFile(directory / name) as wheel:
        names = wheel.namelist()
        require(len(names) == len(set(names)) and all(safe_name(path.replace("/", "")) and ":" not in path
                and not path.startswith("/") and ".." not in Path(path).parts for path in names),
                "Unsafe or duplicate CUDA wheel path")
        metadata_bytes = wheel.read(prefix + "METADATA")
        metadata = BytesParser().parsebytes(metadata_bytes)
        require(metadata["Name"] == "smsd-cuda" and metadata["Version"] == version
                and metadata["Requires-Python"] == ">=3.12"
                and {dependency(value) for value in metadata.get_all("Requires-Dist", [])} == dependencies,
                "CUDA wheel package metadata differs from frozen source")
        description = metadata_bytes.split(b"\n\n", 1)
        require(len(description) == 2 and description[1] == files[root + "README.md"], "CUDA wheel description differs")
        tags = BytesParser().parsebytes(wheel.read(prefix + "WHEEL"))
        require(tags.get_all("Tag") == ["py3-none-manylinux_2_34_x86_64"]
                and tags["Root-Is-Purelib"] == "false", "CUDA wheel tags differ")
        for path in names:
            if path.endswith(".py"):
                require(root + "python/" + path in files and wheel.read(path) == files[root + "python/" + path],
                        "CUDA wheel Python source differs from frozen source")
        require(wheel.read("smsd_cuda/__init__.py") == files[root + "python/smsd_cuda/__init__.py"], "CUDA module missing or changed")
        for license_name in ("LICENSE", "NOTICE"):
            require(wheel.read(prefix + "licenses/" + license_name) == files[root + license_name], "CUDA licence copy differs")
        binary = wheel.read("smsd_cuda/libsmsd_cuda_screen.so")
        require(binary[:6] == b"\x7fELF\x02\x01" and binary[18:20] == b"\x3e\x00", "CUDA library must be Linux x86_64 ELF")
        record_name = prefix + "RECORD"
        records = list(csv.reader(io.StringIO(wheel.read(record_name).decode())))
        require(all(len(row) == 3 for row in records) and len({row[0] for row in records}) == len(records)
                and {row[0] for row in records} == {path for path in names if not path.endswith("/")}, "CUDA RECORD inventory differs")
        for path, checksum, size in records:
            if path == record_name:
                require(checksum == size == "", "CUDA RECORD self entry differs")
                continue
            content = wheel.read(path)
            digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest()).rstrip(b"=").decode()
            require(checksum == "sha256=" + digest and size == str(len(content)), "CUDA RECORD hash or size differs")


def validate(args):
    package = getattr(args, "package", "smsd")
    version, assets, targets = selected_assets(BUILD.read_manifest(args.metadata_file), args.tag, package)
    check_hashes(args.release_dir, assets)
    source = BUILD.read_manifest(args.release_dir / "SOURCE_MANIFEST.json")
    files = check_source(args.release_dir, version, source)
    collector = check_git(args.source_checkout, args.tag, source, files, package)
    if package == "smsd":
        check_wheel_manifest(args.release_dir, version, source, targets, files["PKG-INFO"])
        subprocess.run([sys.executable, str(collector), "--release-dir", str(args.release_dir.resolve()),
                        "--check-only"], check=True)
        distributions = sorted({f"smsd-{version}.tar.gz", *targets.values()})
    else:
        cuda_files = check_cuda_source(args.source_checkout, source)
        check_cuda_wheel(args.release_dir, version, source, cuda_files)
        distributions = [targets["cuda"]]
    require(not args.dist_dir.exists() or not any(args.dist_dir.iterdir()), "Upload directory must be empty")
    args.dist_dir.mkdir(parents=True, exist_ok=True)
    for name in distributions:
        shutil.copyfile(args.release_dir / name, args.dist_dir / name)
    receipt = {"package": package, "version": version, "source_commit": source["source_commit"],
               "files": {name: BUILD.sha256(args.release_dir / name) for name in distributions}}
    args.receipt_file.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("Validated unchanged SDK and all 16 release wheels" if package == "smsd"
          else "Validated the unchanged published CUDA wheel")


def verify_staged(args):
    receipt = BUILD.read_manifest(args.receipt_file)
    version = release_version("v" + receipt["version"])
    files = receipt.get("files", {})
    package = getattr(args, "package", "smsd")
    require(package in ("smsd", "smsd-cuda"), "Unknown upload package")
    require(receipt.get("package", "smsd") == package, "Upload receipt package differs")
    if package == "smsd-cuda":
        require(set(receipt) == {"package", "version", "source_commit", "files"}
                and isinstance(receipt["source_commit"], str) and re.fullmatch(r"[0-9a-f]{40}", receipt["source_commit"]),
                "Invalid CUDA source receipt")
        require(isinstance(files, dict) and set(files) == {cuda_wheel_name(version)}, "Invalid CUDA distribution receipt")
    else:
        require(isinstance(files, dict) and len(files) == len(BUILD.BUILD_IDENTIFIERS) + 1
                and f"smsd-{version}.tar.gz" in files, "Invalid distribution receipt")
    require({path.name for path in args.dist_dir.iterdir()} == set(files), "Upload directory contents differ")
    targets = set()
    for name, expected in files.items():
        require(safe_name(name) and re.fullmatch(r"[0-9a-f]{64}", expected)
                and (args.dist_dir / name).is_file() and not (args.dist_dir / name).is_symlink()
                and BUILD.sha256(args.dist_dir / name) == expected, "Staged distribution hash differs")
        if name.endswith(".whl") and package == "smsd":
            targets.add(BUILD.wheel_target(name, version))
    if package == "smsd":
        require(targets == set(BUILD.BUILD_IDENTIFIERS), "Staged interpreter/platform targets differ")
    print("Final 21-file upload identity check passed" if package == "smsd"
          else "Final one-wheel CUDA upload identity check passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    fetch = commands.add_parser("download")
    fetch.add_argument("--tag", required=True)
    fetch.add_argument("--platform", required=True)
    fetch.add_argument("--repository", required=True)
    fetch.add_argument("--release-dir", type=Path, required=True)
    fetch.add_argument("--metadata-file", type=Path, required=True)
    fetch.set_defaults(function=download)
    check = commands.add_parser("validate", help="Offline validation of already downloaded release files")
    check.add_argument("--tag", required=True)
    check.add_argument("--release-dir", type=Path, required=True)
    check.add_argument("--metadata-file", type=Path, required=True)
    check.add_argument("--source-checkout", type=Path, required=True)
    check.add_argument("--dist-dir", type=Path, required=True)
    check.add_argument("--receipt-file", type=Path, required=True)
    check.set_defaults(function=validate)
    staged = commands.add_parser("verify-staged")
    staged.add_argument("--dist-dir", type=Path, required=True)
    staged.add_argument("--receipt-file", type=Path, required=True)
    staged.set_defaults(function=verify_staged)
    for command in (fetch, check, staged):
        command.add_argument("--package", choices=("smsd", "smsd-cuda"), default="smsd")
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, KeyError, OSError, tarfile.TarError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

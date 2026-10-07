#!/usr/bin/env python3
"""Validate published SMSD distributions before uploading the unchanged files."""

import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tarfile
import tomllib
import xml.etree.ElementTree as ET
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    def unique_pairs(pairs):
        value = {}
        for key, item in pairs:
            require(key not in value, "Duplicate JSON key: " + key)
            value[key] = item
        return value
    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_pairs)


def release_version(tag):
    match = re.fullmatch(r"v(\d+\.\d+\.\d+)", tag)
    require(match is not None, "Publishing requires a stable release tag")
    return match[1]


def safe_name(name):
    return (isinstance(name, str) and name not in ("", ".", "..")
            and not any(char in name for char in "/\\\r\n\0"))


def wheel_platform(name, version):
    prefix = f"smsd-{version}-cp314-cp314-"
    require(name.startswith(prefix) and name.endswith(".whl"),
            "Unexpected release wheel: " + name)
    tag = name[len(prefix):-4]
    if tag == "win_amd64":
        return "windows"
    if tag == "macosx_26_0_arm64":
        return "macos"
    if ("manylinux_2_28_x86_64" in tag.split(".")
            and all(re.fullmatch(r"manylinux_2_\d+_x86_64", item)
                    for item in tag.split("."))):
        return "linux"
    raise ValueError("Unexpected release wheel platform: " + name)


def selected_assets(metadata, tag):
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
    selected = {f"smsd-{version}.tar.gz", "SOURCE_MANIFEST.json",
                "WHEEL_MANIFEST.json", "SHA256SUMS"}
    families = {}
    for name in assets:
        if name.endswith(".whl"):
            family = wheel_platform(name, version)
            require(family not in families, "Multiple " + family + " release wheels")
            families[family] = name
    missing = {"linux", "macos", "windows"} - families.keys()
    require(not missing, "Missing release wheels: " + ", ".join(sorted(missing)))
    selected.update(families.values())
    require(selected <= assets.keys(), "Published release is missing required metadata or source")
    for name in selected:
        asset = assets[name]
        require(asset.get("state") == "uploaded", "Release asset is not uploaded: " + name)
        require(type(asset.get("size")) is int and asset["size"] > 0,
                "Invalid release asset size: " + name)
        require(isinstance(asset.get("digest"), str)
                and re.fullmatch(r"sha256:[0-9a-f]{64}", asset["digest"]),
                "Expected GitHub SHA256 digest: " + name)
    return version, {name: assets[name] for name in selected}, families


def download(args):
    require(args.platform == "all", "Publishing requires platform=all")
    release_version(args.tag)
    require(re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository),
            "Expected owner/repository")
    require(not args.release_dir.exists() or not any(args.release_dir.iterdir()),
            "Release download directory must be empty")
    raw = subprocess.check_output(["gh", "api",
                                   f"repos/{args.repository}/releases/tags/{args.tag}"])
    args.metadata_file.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_file.write_bytes(raw)
    _, assets, _ = selected_assets(read_json(args.metadata_file), args.tag)
    args.release_dir.mkdir(parents=True, exist_ok=True)
    command = ["gh", "release", "download", args.tag, "--repo", args.repository,
               "--dir", str(args.release_dir)]
    for name in sorted(assets):
        command.extend(["--pattern", name])
    subprocess.run(command, check=True)
    print("Downloaded one published SDK and three release wheels with their manifests")


def source_files(path):
    files, roots = {}, set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            parts = PurePosixPath(member.name).parts
            require(parts and not member.name.startswith("/") and ".." not in parts
                    and "\\" not in member.name and ":" not in member.name
                    and (member.isfile() or member.isdir()),
                    "Unsafe SDK member: " + member.name)
            roots.add(parts[0])
            if member.isfile():
                relative = PurePosixPath(*parts[1:]).as_posix()
                require(relative and relative not in files, "Duplicate SDK member: " + relative)
                files[relative] = archive.extractfile(member).read()
    require(len(roots) == 1, "Expected one SDK archive root")
    return files


def check_source(directory, version, manifest):
    name = f"smsd-{version}.tar.gz"
    files = source_files(directory / name)
    require(manifest.get("version") == version and manifest.get("sdist") == name,
            "Source manifest version or filename differs")
    require(manifest.get("sdist_sha256") == sha256(directory / name),
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


def check_git(checkout, tag, manifest, files):
    frozen = git(checkout, "rev-parse", "--verify", manifest["release_ref"] + "^{commit}").decode()
    require(frozen == manifest["source_commit"], "Frozen source ref resolves to another commit")
    tagged = git(checkout, "rev-parse", "--verify", f"refs/tags/{tag}" + "^{commit}").decode()
    require(git(checkout, "rev-parse", "HEAD").decode() == tagged,
            "Collector checkout differs from the selected tag")
    subprocess.run(["git", "-C", str(checkout), "merge-base", "--is-ancestor", frozen, tagged], check=True)
    changed = git(checkout, "diff", "--name-only", "-z", frozen, tagged).split(b"\0")
    require(all(not name or Path(name.decode()).suffix in (".md", ".rst")
                or name == b"CITATION.cff" for name in changed),
            "Tagged production files differ from the frozen SDK")
    tree = git(checkout, "ls-tree", "-r", "-z", "--full-tree", frozen)
    blobs = {}
    for entry in tree.split(b"\0"):
        if entry:
            metadata, name = entry.split(b"\t", 1)
            mode, kind, identity = metadata.split()
            if kind == b"blob" and mode in (b"100644", b"100755"):
                blobs[name.decode("utf-8")] = identity.decode("ascii")
    for name, content in files.items():
        if name == "PKG-INFO":
            continue
        identity = hashlib.sha1(b"blob " + str(len(content)).encode("ascii") + b"\0" + content).hexdigest()
        require(blobs.get(name) == identity, "SDK file differs from frozen source: " + name)
    collector = "scripts/collect-release-wheels.py"
    require((checkout / collector).read_bytes() == files[collector], "Tagged collector differs from SDK")
    return checkout / collector


def check_wheel_manifest(directory, version, source, families, package_metadata):
    manifest = read_json(directory / "WHEEL_MANIFEST.json")
    require(manifest.get("schema_version") == 1, "Unexpected wheel manifest schema")
    for field in ("version", "source_commit", "sdist", "sdist_sha256"):
        require(manifest.get(field) == source[field], "Wheel manifest identity differs: " + field)
    entries = manifest.get("wheels")
    require(isinstance(entries, list) and len(entries) == 3, "Expected three wheel manifest entries")
    recorded = {}
    for entry in entries:
        family = entry.get("platform")
        require(family in families and family not in recorded, "Duplicate or unknown manifest platform")
        require(entry.get("wheel") == families[family], "Wheel manifest filename differs")
        require(entry.get("wheel_sha256") == sha256(directory / families[family]),
                "Wheel manifest SHA256 differs")
        with zipfile.ZipFile(directory / families[family]) as wheel:
            require(wheel.read(f"smsd-{version}.dist-info/METADATA") == package_metadata,
                    "Wheel package metadata or long description differs from the final SDK")
        require(entry.get("native_debug_tests") == "passed"
                and entry.get("installed_wheel_tests") == "passed",
                "Wheel manifest lacks successful native and installed tests")
        recorded[family] = entry
    require(set(recorded) == {"linux", "macos", "windows"}, "Wheel manifest is incomplete")


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
        actual = sha256(path)
        require(path.stat().st_size == asset["size"] and asset["digest"] == "sha256:" + actual,
                "GitHub asset digest or size differs: " + name)
        if name != "SHA256SUMS":
            require(sums.get(name) == actual, "Release SHA256SUMS differs: " + name)


def validate(args):
    version, assets, families = selected_assets(read_json(args.metadata_file), args.tag)
    check_hashes(args.release_dir, assets)
    source = read_json(args.release_dir / "SOURCE_MANIFEST.json")
    files = check_source(args.release_dir, version, source)
    collector = check_git(args.source_checkout, args.tag, source, files)
    check_wheel_manifest(args.release_dir, version, source, families, files["PKG-INFO"])
    subprocess.run([sys.executable, str(collector), "--release-dir", str(args.release_dir.resolve()),
                    "--check-only"], check=True)
    distributions = sorted({f"smsd-{version}.tar.gz", *families.values()})
    require(not args.dist_dir.exists() or not any(args.dist_dir.iterdir()), "Upload directory must be empty")
    args.dist_dir.mkdir(parents=True, exist_ok=True)
    for name in distributions:
        shutil.copyfile(args.release_dir / name, args.dist_dir / name)
    receipt = {"version": version, "source_commit": source["source_commit"],
               "files": {name: sha256(args.release_dir / name) for name in distributions}}
    args.receipt_file.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("Validated unchanged SDK and all three release wheels")


def verify_staged(args):
    receipt = read_json(args.receipt_file)
    version = release_version("v" + receipt["version"])
    files = receipt.get("files", {})
    require(isinstance(files, dict) and len(files) == 4
            and f"smsd-{version}.tar.gz" in files, "Invalid distribution receipt")
    require({path.name for path in args.dist_dir.iterdir()} == set(files), "Upload directory contents differ")
    families = set()
    for name, expected in files.items():
        require(safe_name(name) and re.fullmatch(r"[0-9a-f]{64}", expected)
                and (args.dist_dir / name).is_file() and not (args.dist_dir / name).is_symlink()
                and sha256(args.dist_dir / name) == expected, "Staged distribution hash differs")
        if name.endswith(".whl"):
            families.add(wheel_platform(name, version))
    require(families == {"linux", "macos", "windows"}, "Staged platforms differ")
    print("Final four-file upload identity check passed")


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
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, KeyError, OSError, tarfile.TarError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

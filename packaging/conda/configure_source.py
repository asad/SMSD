# SPDX-License-Identifier: Apache-2.0
"""Prepare conda variant inputs from the exact frozen full-source archive bytes."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
from urllib.parse import urlsplit

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--url", required=True, help="The immutable HTTPS URL serving these exact bytes")
    parser.add_argument("--commit", required=True, help="The reviewed full source commit")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.commit):
        parser.error("--commit must be a full lowercase SHA-1")
    url = urlsplit(args.url)
    if url.scheme != "https" or not url.netloc or url.username or url.password or url.fragment:
        parser.error("--url must be a credential-free immutable HTTPS source URL")
    if args.commit not in url.path:
        parser.error("The source URL path must identify the full source commit")
    with tarfile.open(args.archive, "r:*") as archive:
        members = archive.getmembers()
        roots = {PurePosixPath(member.name).parts[0] for member in members if member.name}
        if len(roots) != 1:
            raise ValueError("Expected one top-level source directory")
        prefix = roots.pop() + "/"

        def read(relative):
            member = archive.getmember(prefix + relative)
            if not member.isfile() or member.size > 1024 * 1024:
                raise ValueError("Invalid source metadata: " + relative)
            return archive.extractfile(member).read().decode("utf-8")

        core = tomllib.loads(read("pyproject.toml"))["project"]
        cuda = tomllib.loads(read("packaging/smsd-cuda/pyproject.toml"))["project"]
        version = core["version"]
        native = re.search(r"^project\(smsd VERSION ([0-9.]+)\b", read("cpp/CMakeLists.txt"), re.M)
        if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
            raise ValueError("Expected a stable source version")
        if cuda["version"] != version or native is None or native[1] != version or f"smsd=={version}" not in cuda["dependencies"]:
            raise ValueError("Core, CUDA and native source versions must match exactly")
    digest = hashlib.sha256()
    with args.archive.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    variants = {"python": ["3.12", "3.13", "3.14", "3.15"],
                "smsd_version": [version], "smsd_source_commit": [args.commit],
                "smsd_source_url": [args.url], "smsd_source_sha256": [digest.hexdigest()]}
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(variants, stream, indent=2)
        stream.write("\n")
    print("Created YAML-compatible variant inputs; remote URL byte equality and Git identity require release review")


if __name__ == "__main__":
    main()

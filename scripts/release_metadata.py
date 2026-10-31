#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Update release versions and extract the matching GitHub release text locally."""
import argparse
from datetime import date
from pathlib import Path
import re


VERSION_FILES = (
    "pyproject.toml", "pom.xml", "java/pom.xml", "java/core/pom.xml", "java/cli/pom.xml",
    "cpp/CMakeLists.txt", "python/smsd/__init__.py", "Dockerfile", "Dockerfile.cuda",
    "packaging/smsd-cuda/pyproject.toml", "packaging/smsd-cuda/python/smsd_cuda/__init__.py",
    "scripts/docker_install.py", "scripts/docker_smoke.py", "scripts/prepare_docker_context.py",
    "packaging/mcp/Dockerfile", "examples/huggingface/Dockerfile",
)


def stable_version(version):
    if not re.fullmatch(r"(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)\.(?:0|[1-9][0-9]*)", version):
        raise ValueError("Expected a stable major.minor.patch version")
    return version


def set_version(root, version, released):
    stable_version(version)
    date.fromisoformat(released)
    project = root / "pyproject.toml"
    found = re.search(r'^version = "([^"]+)"$', project.read_text(), re.M)
    if not found:
        raise ValueError("Package version is missing")
    old = stable_version(found[1])
    updates = {}
    for name in VERSION_FILES:
        path = root / name
        text = path.read_text()
        if old not in text:
            raise ValueError("Current version is missing from " + name)
        updates[path] = text.replace(old, version)
    citation = root / "CITATION.cff"
    text, count = re.subn(r"^version: .+$", "version: " + version, citation.read_text(), flags=re.M)
    if count != 1:
        raise ValueError("Expected one citation version")
    if re.search(r"^date-released:", text, re.M):
        text = re.sub(r"^date-released:.*$", "date-released: " + released, text, flags=re.M)
    else:
        text = text.replace("version: " + version + "\n", "version: " + version + "\ndate-released: " + released + "\n", 1)
    updates[citation] = text
    for path, text in updates.items():
        path.write_text(text)


def release_notes(root, version):
    stable_version(version)
    text = (root / "docs/RELEASE_NOTES.md").read_text()
    heading = "# SMSD Pro " + version
    found = re.search(r"^" + re.escape(heading) + r"\s*$", text, re.M)
    if not found:
        raise ValueError("Release notes are missing for " + version)
    section = text[found.end():].lstrip("\n")
    section = re.split(r"^# ", section, maxsplit=1, flags=re.M)[0].rstrip()
    if not section:
        raise ValueError("Release notes are empty")
    return section + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    commands = parser.add_subparsers(dest="command", required=True)
    update = commands.add_parser("set-version")
    update.add_argument("version")
    update.add_argument("--date", default=date.today().isoformat())
    notes = commands.add_parser("notes")
    notes.add_argument("version")
    notes.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "set-version":
        set_version(args.root, args.version, args.date)
    else:
        text = release_notes(args.root, args.version)
        with args.output.open("x") as output:
            output.write(text)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error

#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Check that the OpenMP runtime can be bundled for the requested macOS target."""

import argparse
import os
from pathlib import Path
import re
import subprocess


def version_tuple(value):
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+){1,2}", value):
        raise ValueError("Expected a macOS major.minor[.patch] version")
    fields = tuple(int(field) for field in value.split("."))
    return fields + (0,) * (3 - len(fields))


def minimum_versions(load_commands):
    versions = []
    command = None
    for line in load_commands.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0] == "cmd":
            command = fields[1]
        elif len(fields) == 2 and ((command == "LC_BUILD_VERSION" and fields[0] == "minos")
                                  or (command == "LC_VERSION_MIN_MACOSX" and fields[0] == "version")):
            versions.append(version_tuple(fields[1]))
    if not versions:
        raise ValueError("OpenMP library has no macOS minimum-version load command")
    return versions


def check(prefix, deployment_target):
    target = version_tuple(deployment_target)
    library = prefix / "lib/libomp.dylib"
    if not library.is_file() or not (prefix / "include/omp.h").is_file():
        raise ValueError("OpenMP prefix must contain lib/libomp.dylib and include/omp.h")
    architectures = subprocess.check_output(["lipo", "-archs", str(library)], text=True).split()
    if "arm64" not in architectures:
        raise ValueError("OpenMP library does not support arm64")
    commands = subprocess.check_output(["otool", "-l", "-arch", "arm64", str(library)], text=True)
    if any(minimum > target for minimum in minimum_versions(commands)):
        raise ValueError("OpenMP requires a newer macOS version; set SMSD_MACOS_OPENMP_PREFIX "
                         "to an arm64 OpenMP installation built for macOS " + deployment_target + " or earlier")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--deployment-target", default=os.environ.get("MACOSX_DEPLOYMENT_TARGET", "14.0"))
    args = parser.parse_args()
    prefix = args.prefix or os.environ.get("SMSD_MACOS_OPENMP_PREFIX")
    if prefix is None:
        prefix = subprocess.check_output(["brew", "--prefix", "libomp"], text=True).strip()
    check(Path(prefix), args.deployment_target)
    print("Validated arm64 OpenMP runtime for macOS " + args.deployment_target)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error

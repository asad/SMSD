#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Run all native CTest suites with Debug assertions before building a wheel."""

import argparse
import importlib.util
from pathlib import Path
import tempfile


_BUILD_SPEC = importlib.util.spec_from_file_location(
    "smsd_release_build", Path(__file__).resolve().with_name("build_python_wheels.py"))
BUILD = importlib.util.module_from_spec(_BUILD_SPEC)
_BUILD_SPEC.loader.exec_module(BUILD)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    args = parser.parse_args()
    try:
        with tempfile.TemporaryDirectory(prefix="smsd-native-tests-") as temporary:
            count = BUILD.run_native_debug_tests(args.source.resolve(), Path(temporary) / "build", 2)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    print("Native Debug suites passed:", count)


if __name__ == "__main__":
    main()

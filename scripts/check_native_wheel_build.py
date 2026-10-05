#!/usr/bin/env python3
"""Run all native CTest suites with Debug assertions before building a wheel."""

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="smsd-native-tests-") as temporary:
        build = Path(temporary) / "build"
        subprocess.run(["cmake", "-S", str(args.source.resolve() / "cpp"), "-B", str(build),
                        "-DCMAKE_BUILD_TYPE=Debug", "-DSMSD_BUILD_TESTS=ON",
                        "-DSMSD_BUILD_PYTHON=OFF", "-DSMSD_BUILD_OPENMP=ON",
                        "-DSMSD_BUILD_METAL=OFF", "-DSMSD_BUILD_CUDA=OFF",
                        "-DSMSD_WITH_RDKIT=OFF"], check=True)
        subprocess.run(["cmake", "--build", str(build), "--config", "Debug", "--parallel", "2"], check=True)
        listing = subprocess.check_output(["ctest", "--test-dir", str(build),
                                           "--build-config", "Debug", "--show-only=json-v1"], text=True)
        count = len(json.loads(listing)["tests"])
        if count < 12:
            raise SystemExit("Expected all 12 native suites; found " + str(count))
        subprocess.run(["ctest", "--test-dir", str(build), "--build-config", "Debug",
                        "--output-on-failure"], check=True)
        print("Native Debug suites passed:", count)


if __name__ == "__main__":
    main()

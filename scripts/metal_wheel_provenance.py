#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Test a Metal wheel against the canonical source on Apple silicon, then collect it with its proof.

The Metal native tests are built from the source distribution, and the core and smsd-metal Python
suites run against the installed core and Metal wheels in a new environment. The proof is written
and the wheel collected only when every native and Metal case passed on the device and the
installed files are those of the two wheels.
"""

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROOF = _load("cuda_proof", "cuda_wheel_provenance.py")
COLLECTOR = _load("addon_collector", "collect_cuda_wheels.py")
BUILD = COLLECTOR.BUILD

NATIVE_TESTS = {
    *(f"smsd_metal_failure_{failure}_tests"
      for failure in ("pipeline", "buffer", "command", "encoder", "status", "disabled")),
    *(f"smsd_metal_fingerprint_{mode}_tests"
      for mode in ("contract", "concurrent", "device", "queue", "pipeline", "function", "buffer", "contents",
                   "exception", "command", "encoder", "status", "disabled")),
    "smsd_metal_resident_device_tests", "smsd_metal_threshold_regression_tests", "smsd_metal_screening_regression_tests",
}
# The checkout paths the Metal wheel is built from: its package and the native sources it compiles.
METAL_SOURCES = ("packaging/smsd-metal", "cpp/metal", "cpp/include", "cpp/CMakeLists.txt")
# Cases that show smsd found and attached the installed package on the device.
DEVICE_CASES = {"test_smsd_finds_the_package", "test_smsd_attaches_the_device"}

# Run by the test environment's interpreter: what it loads, from where, and the accelerator it attaches.
_IDENTITY = """
import json, sys, smsd, smsd_metal
from smsd import _accelerator
library = smsd_metal._load()
print(json.dumps({"prefix": sys.prefix, "core_module": smsd.__file__, "native_module": smsd._smsd.__file__,
                  "module": smsd_metal.__file__, "library": library and library._name,
                  "versions": [smsd.__version__, smsd._smsd.__version__, smsd_metal.__version__,
                               smsd_metal.__smsd_version__],
                  "accelerator": _accelerator.active()}))
"""


def native_tests(sdk, work):
    """Build the source distribution's tests with Metal and run the Metal ones; their JUnit results."""
    build, results = work / "native", work / "native.xml"
    subprocess.run(["cmake", "-S", str(sdk / "cpp"), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release",
                    "-DSMSD_BUILD_TESTS=ON", "-DSMSD_BUILD_METAL=ON", "-DSMSD_BUILD_CUDA=OFF",
                    "-DSMSD_BUILD_OPENMP=OFF"], check=True)
    subprocess.run(["cmake", "--build", str(build), "--parallel"], check=True)
    subprocess.run(["ctest", "--test-dir", str(build), "--output-on-failure", "--no-tests=error",
                    "-R", "^smsd_metal_", "--output-junit", str(results)], check=True)
    return results


def python_tests(python, sdk, checkout, cpu_wheel, wheel, work):
    """Install both wheels in a new environment and run the core and smsd-metal suites there."""
    subprocess.run([python, "-m", "venv", str(work / "venv")], check=True)
    venv, results = str(work / "venv/bin/python"), work / "python.xml"
    subprocess.run([venv, "-m", "pip", "install", "-q", "pytest==9.1.1", "pytest-timeout==2.4.0"], check=True)
    subprocess.run([venv, "-m", "pip", "install", "-q", "--no-deps", str(cpu_wheel), str(wheel)], check=True)
    subprocess.run([venv, "-m", "pytest", str(sdk / "python/tests"), str(checkout / "packaging/smsd-metal/tests"),
                    "-q", "-rs", "-p", "no:cacheprovider", "--junitxml=" + str(results),
                    "-c", str(sdk / "pyproject.toml"), "--rootdir", str(sdk)], cwd=work, check=True)
    return results, json.loads(subprocess.check_output([venv, "-c", _IDENTITY], cwd=work))


def proof_from_results(source, wheel, cpu_wheel, identity, native_xml, python_xml):
    """The Metal wheel proof, once the native and installed results are complete and successful."""
    version = source["version"]
    if identity["versions"] != [version] * 4 or not identity["accelerator"] or identity["accelerator"][0] != "Metal":
        raise ValueError("The tests did not load this smsd and smsd-metal with a Metal device")
    PROOF.installed_files(cpu_wheel, "smsd", identity["prefix"], identity["core_module"], identity["native_module"])
    PROOF.installed_files(wheel, "smsd_metal", identity["prefix"], identity["module"], identity["library"])
    native = PROOF.junit_results(native_xml)
    if (native["skipped"] or len(native["passed"]) != len(NATIVE_TESTS)
            or {name for _, name in native["passed"]} != NATIVE_TESTS):
        raise ValueError("Native XML does not contain the complete passed Metal suite")
    python = PROOF.junit_results(python_xml)
    metal = [case for case in python["cases"] if "smsd-metal" in case["classname"].split(".")]
    if any(case["outcome"] != "passed" for case in metal) or not DEVICE_CASES <= {case["name"] for case in metal}:
        raise ValueError("Python XML has a skipped or missing Metal case")
    if not any("test_smsd" in classname.split(".") for classname, _ in python["passed"]):
        raise ValueError("Python XML has no passed core suite")
    return {"schema_version": 1, "source_commit": source["source_commit"],
            "sdist_sha256": source["sdist_sha256"], "smsd_version": version,
            "wheel": wheel.name, "wheel_sha256": PROOF.sha256(wheel),
            "installed_wheel_tests": "passed", "native_gpu_tests": "passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True, help="the smsd-metal wheel")
    parser.add_argument("--cpu-wheel", type=Path, required=True, help="the smsd wheel of the same release for this Mac")
    parser.add_argument("--work", type=Path, required=True, help="a new directory for the builds and results")
    parser.add_argument("--python", default=sys.executable, help="the interpreter of the test environment")
    parser.add_argument("--checkout", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    work = args.work.resolve()
    if work.is_relative_to(args.release_dir.resolve()):
        raise ValueError("Private validation evidence must be outside the release directory")
    work.mkdir(parents=True)
    manifest = args.release_dir / "SOURCE_MANIFEST.json"
    source = BUILD.read_manifest(manifest)
    source = BUILD.checked_source(args.release_dir / source["sdist"], manifest)
    subprocess.run([sys.executable, str(args.checkout / "scripts/build_python_wheels.py"), "check-source",
                    "--sdist", str((args.release_dir / source["sdist"]).resolve()),
                    "--source-manifest", str(manifest.resolve()), "--release-commit", source["source_commit"]],
                   cwd=args.checkout, check=True, stdout=subprocess.DEVNULL)
    # The Metal suite runs from the checkout, which must hold the frozen package and the native sources
    # its library is compiled from unchanged.
    if (subprocess.run(["git", "-C", str(args.checkout), "diff", "--quiet", source["source_commit"], "--",
                        *METAL_SOURCES]).returncode
            or subprocess.check_output(["git", "-C", str(args.checkout), "status", "--porcelain",
                                        "--untracked-files=all", "--", *METAL_SOURCES])):
        raise ValueError("The checkout's smsd-metal package or Metal sources differ from the frozen source")
    sdk = work / "sdk"
    for name, data in BUILD.read_source_archive(args.release_dir / source["sdist"]).items():
        (sdk / name).parent.mkdir(parents=True, exist_ok=True)
        (sdk / name).write_bytes(data)
    native_xml = native_tests(sdk, work)
    python_xml, identity = python_tests(args.python, sdk, args.checkout.resolve(), args.cpu_wheel.resolve(),
                                        args.wheel.resolve(), work)
    proof_path = work / "metal-wheel-build-proof.json"
    BUILD.write_manifest(proof_path, proof_from_results(source, args.wheel, args.cpu_wheel, identity,
                                                        native_xml, python_xml))
    COLLECTOR.collect(args.release_dir, args.wheel, proof_path, args.checkout, "smsd-metal")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, ET.ParseError, zipfile.BadZipFile, subprocess.CalledProcessError) as error:
        raise SystemExit(str(error)) from error

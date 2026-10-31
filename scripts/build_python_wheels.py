#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Prepare one sdist or build repaired CPython 3.12–3.15 CPU wheels from it.

Install build/twine for preparation and cibuildwheel==4.2.1 for wheel builds.
Linux builds use local Docker or Podman; Windows/macOS require that host OS.
"""

import argparse
from email.parser import BytesParser
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import subprocess
import sys
import tarfile
import tempfile
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
import xml.etree.ElementTree as ET
import zipfile


SOURCE_FIELDS = ("version", "sdist", "sdist_sha256", "source_file_count",
                 "source_commit", "release_ref", "working_tree_clean")
PYTHON_VERSIONS = ("3.12", "3.13", "3.14", "3.15")
ARCHITECTURES = {"linux": ("x86_64", "aarch64"), "macos": ("arm64",), "windows": ("AMD64",)}
BUILD_IDENTIFIERS = {
    "cp" + version.replace(".", "") + "-" + tag:
        {"platform": family, "architecture": architecture, "python_version": version}
    for version in PYTHON_VERSIONS
    for family, architecture, tag in (("linux", "x86_64", "manylinux_x86_64"),
                                      ("linux", "aarch64", "manylinux_aarch64"),
                                      ("macos", "arm64", "macosx_arm64"),
                                      ("windows", "x86_64", "win_amd64"))
}
# CPU-only configurations register 35 native suites on Linux and Windows;
# macOS adds six Metal failure cases, and RDKit or CUDA builds add their own.
NATIVE_TEST_FLOOR = 35
SHARED_VALIDATION_CHECKS = ("cpp_release", "cpp_asan", "java8", "java25", "python", "examples", "check_source")
NATIVE_CMAKE_OPTIONS = ("-DCMAKE_BUILD_TYPE=Debug", "-DSMSD_BUILD_TESTS=ON", "-DSMSD_BUILD_PYTHON=OFF",
                        "-DSMSD_BUILD_OPENMP=ON", "-DSMSD_BUILD_METAL=OFF", "-DSMSD_BUILD_CUDA=OFF",
                        "-DSMSD_WITH_RDKIT=OFF")


def wheel_target(name, version):
    """Return the interpreter/platform target of a supported release wheel."""
    fields = name[:-4].split("-")
    if (not name.endswith(".whl") or len(fields) != 5 or fields[:2] != ["smsd", version]
            or fields[2] != fields[3] or fields[2] not in {"cp" + value.replace(".", "") for value in PYTHON_VERSIONS}):
        raise ValueError("Unexpected package, version or interpreter: " + name)
    platform_tag = fields[4]
    if platform_tag == "win_amd64":
        target_tag = "win_amd64"
    elif re.fullmatch(r"macosx_(?:1[1-9]|2[0-6])_0_arm64", platform_tag):
        target_tag = "macosx_arm64"
    else:
        target_tag = None
        for architecture in ARCHITECTURES["linux"]:
            tags = platform_tag.split(".")
            if ("manylinux_2_28_" + architecture in tags
                    and all(re.fullmatch(r"manylinux_2_\d+_" + architecture, tag) for tag in tags)):
                target_tag = "manylinux_" + architecture
                break
        if target_tag is None:
            raise ValueError("Unexpected release platform: " + platform_tag)
    return fields[2] + "-" + target_tag


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_blob_id(data):
    """Return the Git blob identity of file contents."""
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\0" + data).hexdigest()


def tree_blobs(listing):
    """Map each regular file in `git ls-tree -r -z --full-tree` output to its blob identity."""
    blobs = {}
    for entry in listing.split(b"\0"):
        if entry:
            metadata, name = entry.split(b"\t", 1)
            mode, kind, identity = metadata.split()
            if kind == b"blob" and mode in (b"100644", b"100755"):
                blobs[name.decode("utf-8")] = identity.decode("ascii")
    return blobs


def read_manifest(path):
    def unique_pairs(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("Duplicate JSON key: " + key)
            value[key] = item
        return value
    manifest = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_pairs)
    if not isinstance(manifest, dict):
        raise ValueError("Expected a manifest object: " + path.name)
    return manifest


def read_source_archive(path):
    """Return the files of a one-root source archive, keyed by their relative POSIX path."""
    files, roots = {}, set()
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            name = PurePosixPath(member.name)
            relative = "/".join(name.parts[1:])
            if (not name.parts or name.is_absolute() or ".." in name.parts or "\\" in member.name
                    or ":" in member.name or not (member.isdir() or member.isfile())
                    or (member.isfile() and not relative)):
                raise ValueError("Unsafe source archive member: " + member.name)
            roots.add(name.parts[0])
            if member.isfile():
                if relative in files:
                    raise ValueError("Duplicate source archive member: " + relative)
                files[relative] = archive.extractfile(member).read()
    if len(roots) != 1:
        raise ValueError("Expected one source archive root")
    return files


def inspect_source(path):
    files = read_source_archive(path)
    required = {"pyproject.toml", "pom.xml", "java/pom.xml", "PKG-INFO", "python/README.md",
                "python/smsd/__init__.py", "cpp/CMakeLists.txt",
                "cpp/include/smsd/mcs.hpp", "cpp/bindings/pybind11/smsd_bindings.cpp",
                "python/tests/test_smsd.py", "python/tests/test_aromaticity_parity.py",
                "java/src/test/resources/com/bioinception/smsd/aromaticity_parity.json",
                "scripts/check_python_wheel.py",
                "scripts/check_native_wheel_build.py", "scripts/build_python_wheels.py",
                "scripts/collect-release-wheels.py",
                "scripts/cibuildwheel.toml", "scripts/repair_windows_wheel.py",
                "licenses/msvc/README.md", "licenses/msvc/LICENSE-2022.docx",
                "licenses/msvc/LICENSE-2026.docx", "LICENSE", "NOTICE",
                "Dockerfile", "Dockerfile.cuda", ".dockerignore", "scripts/docker_install.py",
                "scripts/docker_smoke.py", "scripts/prepare_docker_context.py",
                "java/core/pom.xml", "java/cli/pom.xml",
                "java/src/main/java/com/bioinception/smsd/core/MCS.java",
                "java/src/main/java/com/bioinception/smsd/core/Substructure.java",
                "java/cli/src/main/java/com/bioinception/smsd/cli/SMSDcli.java",
                "licenses/cdk/LICENSE", "licenses/beam/LICENSE", "licenses/vecmath/LICENSE"}
    missing = sorted(required - files.keys())
    if missing:
        raise ValueError("Incomplete source distribution: " + ", ".join(missing))
    project = tomllib.loads(files["pyproject.toml"].decode())["project"]
    version = project["version"]
    metadata = BytesParser().parsebytes(files["PKG-INFO"])
    description = files["PKG-INFO"].split(b"\n\n", 1)
    if (project["name"] != "smsd" or metadata["Name"] != "smsd" or metadata["Version"] != version
            or len(description) != 2 or description[1] != files["python/README.md"]):
        raise ValueError("Source package metadata or description differs")
    if project.get("requires-python") != ">=3.12" or metadata["Requires-Python"] != project["requires-python"]:
        raise ValueError("Source Python requirement must agree with >=3.12")
    python_version = re.search(rb'^__version__ = "([^"]+)"', files["python/smsd/__init__.py"], re.M)
    cpp_version = re.search(rb'project\(smsd VERSION ([^ ]+)', files["cpp/CMakeLists.txt"])
    if not python_version or not cpp_version or python_version[1].decode() != version or cpp_version[1].decode() != version:
        raise ValueError("Python, CMake and package versions differ")
    for pom in ("pom.xml", "java/pom.xml"):
        if ET.fromstring(files[pom]).findtext("{http://maven.apache.org/POM/4.0.0}version") != version:
            raise ValueError("Java and package versions differ: " + pom)
    return {"version": version, "sdist": path.name, "sdist_sha256": sha256(path),
            "source_file_count": len(files)}


def write_manifest(path, manifest):
    """Create a manifest, retaining an existing identical document unchanged."""
    if path.is_symlink():
        raise ValueError("Manifest must be a regular file: " + path.name)
    if path.exists():
        if read_manifest(path) != manifest:
            raise ValueError("Existing manifest differs: " + path.name)
        return
    with path.open("x", encoding="utf-8") as output:
        output.write(json.dumps(manifest, indent=2) + "\n")


def checked_source(path, manifest_path=None):
    inspected = inspect_source(path)
    manifest_path = manifest_path or path.parent / "SOURCE_MANIFEST.json"
    recorded = read_manifest(manifest_path)
    for field, value in inspected.items():
        if recorded.get(field) != value:
            raise ValueError("Source distribution and manifest differ: " + field)
    for field in ("source_commit", "release_ref"):
        if not isinstance(recorded.get(field), str) or not re.fullmatch(r"[0-9a-f]{40}", recorded[field]):
            raise ValueError("Expected a full frozen Git identity: " + field)
    if recorded["release_ref"] != recorded["source_commit"] or recorded.get("working_tree_clean") is not True:
        raise ValueError("Source manifest must identify one clean frozen commit")
    if recorded.get("schema_version", 1) != 1:
        raise ValueError("Unexpected source manifest schema")
    return {"schema_version": 1, **{field: recorded[field] for field in SOURCE_FIELDS}}


def verify_source_tree(path, commit):
    """Require every packaged source file to belong to the frozen Git tree."""
    blobs = tree_blobs(subprocess.check_output(["git", "ls-tree", "-r", "-z", "--full-tree", commit]))
    files = read_source_archive(path)
    docker_java_inputs = {name for name in blobs if name.startswith(("java/src/main/", "java/cli/src/main/",
        "java/cli/src/scripts/", "licenses/cdk/", "licenses/beam/", "licenses/vecmath/"))}
    if not docker_java_inputs <= files.keys():
        raise ValueError("Source archive omits frozen Java/Docker inputs: " + ", ".join(sorted(docker_java_inputs - files.keys())))
    converted = []
    for relative, data in files.items():
        if relative == "PKG-INFO" or blobs.get(relative) == git_blob_id(data):
            continue
        normalised = data.replace(b"\r\n", b"\n")
        if (blobs.get(relative) != git_blob_id(normalised) or b"\r\n" in normalised
                or data != normalised.replace(b"\n", b"\r\n")):
            raise ValueError("Source archive differs from frozen Git tree: " + relative)
        converted.append(relative)
    if converted:
        info_attributes = Path(subprocess.check_output(
            ["git", "rev-parse", "--git-path", "info/attributes"], text=True).strip())
        if info_attributes.is_file() and info_attributes.read_bytes().strip():
            raise ValueError("Cannot verify CRLF source files with local Git attributes")
        environment = dict(os.environ, GIT_ATTR_NOSYSTEM="1")
        output = subprocess.check_output(["git", "-c", "core.attributesFile=" + os.devnull,
            "check-attr", "--source=" + commit, "-z", "--stdin", "text", "eol", "filter"],
            input=b"\0".join(name.encode("utf-8") for name in converted) + b"\0", env=environment)
        fields, attributes = output.rstrip(b"\0").split(b"\0"), {}
        for name, attribute, value in zip(fields[::3], fields[1::3], fields[2::3]):
            attributes.setdefault(name.decode("utf-8"), {})[attribute.decode("ascii")] = value.decode("utf-8")
        for relative in converted:
            declared = attributes.get(relative, {})
            if (declared.get("text") != "set" or declared.get("eol") != "crlf"
                    or declared.get("filter") not in ("unset", "unspecified")):
                raise ValueError("Source archive differs from frozen Git tree: " + relative)


def run_native_debug_tests(source_dir, build_dir, jobs, environment=None, cwd=None):
    """Build every native Debug test suite from an SDK checkout and run all but the slow ones."""
    def run(command):
        subprocess.run(command, env=environment, cwd=cwd, check=True)
    run(["cmake", "-S", str(source_dir / "cpp"), "-B", str(build_dir), *NATIVE_CMAKE_OPTIONS])
    run(["cmake", "--build", str(build_dir), "--config", "Debug", "--parallel", str(jobs)])
    listing = subprocess.check_output(["ctest", "--test-dir", str(build_dir), "--build-config", "Debug",
                                       "--show-only=json-v1"], env=environment, cwd=cwd, text=True)
    tests = json.loads(listing).get("tests", [])
    if len(tests) < NATIVE_TEST_FLOOR:
        raise ValueError("Expected at least " + str(NATIVE_TEST_FLOOR) + " native Debug test suites; found "
                         + str(len(tests)))
    if any(prop.get("name") == "DISABLED" and prop.get("value")
           for test in tests for prop in test.get("properties", [])):
        raise ValueError("Native Debug test suites must not be disabled")
    run(["ctest", "--test-dir", str(build_dir), "--build-config", "Debug", "--output-on-failure",
         "--label-exclude", "slow", "--parallel", str(jobs)])
    return len(tests)


def proof_digest(proof):
    data = json.dumps(proof, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def checked_shared_validation(source, validation):
    """Check the shared frozen-source results without retaining private evidence."""
    fields = {"source_commit", "sdist_sha256", *SHARED_VALIDATION_CHECKS}
    if not isinstance(validation, dict) or set(validation) != fields:
        raise ValueError("Expected the shared frozen-validation fields only")
    for field in ("source_commit", "sdist_sha256"):
        if validation[field] != source[field]:
            raise ValueError("Shared validation source identity differs: " + field)
    if any(validation[field] != "passed" for field in SHARED_VALIDATION_CHECKS):
        raise ValueError("Shared validation lacks successful frozen-source checks")
    return dict(validation)


def wheel_entry(wheel, identifier, source, proof, shared_validation=None):
    """Validate a checked build receipt and retain only public release fields."""
    for field in SOURCE_FIELDS:
        if proof.get(field) != source[field]:
            raise ValueError("Wheel provenance source identity differs: " + field)
    if proof.get("wheel") != wheel.name or proof.get("wheel_sha256") != sha256(wheel):
        raise ValueError("Wheel provenance filename or SHA256 differs")
    target = BUILD_IDENTIFIERS[identifier]
    if (proof.get("build_identifier") != identifier or wheel_target(wheel.name, source["version"]) != identifier
            or any(proof.get(field) != value for field, value in target.items())):
        raise ValueError("Wheel provenance build platform differs")
    if proof.get("build_system") not in ("local", "cibuildwheel"):
        raise ValueError("Expected provenance from a checked wheel build")
    if shared_validation is not None:
        validation = checked_shared_validation(source, shared_validation)
        if (proof.get("validation_mode") != "lean"
                or proof.get("shared_validation_sha256") != proof_digest(validation)):
            raise ValueError("Wheel provenance shared validation differs")
        if any(field in proof for field in ("native_debug_tests", "installed_wheel_tests")):
            raise ValueError("Lean provenance must not claim per-wheel full-suite results")
        if proof["build_system"] == "cibuildwheel":
            check, forbidden = "platform_recipe_checks", "installed_wheel_identity"
        else:
            check, forbidden = "installed_wheel_identity", "platform_recipe_checks"
        if forbidden in proof:
            raise ValueError("Wheel provenance claims checks outside its build system: " + forbidden)
        checks = (check, "installed_wheel_smoke")
        validation_fields = ("validation_mode", "shared_validation_sha256", *checks)
    else:
        if "validation_mode" in proof:
            raise ValueError("Lean provenance requires shared frozen validation")
        checks = validation_fields = ("native_debug_tests", "installed_wheel_tests")
    for field in checks:
        if proof.get(field) != "passed":
            raise ValueError("Wheel provenance lacks successful checks: " + field)
    fields = (*SOURCE_FIELDS, "wheel", "wheel_sha256", "build_identifier", "platform", "architecture", "python_version", "build_system",
              *validation_fields)
    retained = {field: proof[field] for field in fields}
    if proof["build_system"] == "cibuildwheel":
        if proof.get("cibuildwheel") != "4.2.1":
            raise ValueError("Expected cibuildwheel 4.2.1 provenance")
        retained["cibuildwheel"] = proof["cibuildwheel"]
    return {**target, "build_identifier": identifier, "wheel": wheel.name, "wheel_sha256": sha256(wheel),
            **{field: proof[field] for field in validation_fields},
            "build_provenance": retained, "provenance_sha256": proof_digest(retained)}


def checked_wheel_manifest(directory, source, selected, manifest):
    if manifest.get("schema_version") not in (2, 3):
        raise ValueError("Unexpected wheel manifest schema")
    shared_validation = None
    if manifest["schema_version"] == 3:
        shared_validation = checked_shared_validation(source, manifest.get("shared_validation"))
        if manifest.get("shared_validation_sha256") != proof_digest(shared_validation):
            raise ValueError("Wheel manifest shared validation digest differs")
    for field in ("version", "source_commit", "sdist", "sdist_sha256"):
        if manifest.get(field) != source[field]:
            raise ValueError("Wheel manifest source identity differs: " + field)
    entries = manifest.get("wheels")
    if not isinstance(entries, list) or len(entries) != len(selected):
        raise ValueError("Wheel manifest entries differ from selected wheels")
    recorded = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Expected a wheel manifest entry")
        identifier = entry.get("build_identifier")
        if identifier not in selected or identifier in recorded:
            raise ValueError("Duplicate or unknown wheel manifest target")
        proof = entry.get("build_provenance")
        if not isinstance(proof, dict) or wheel_entry(selected[identifier], identifier, source, proof, shared_validation) != entry:
            raise ValueError("Wheel manifest lacks matching checked-build provenance")
        recorded.add(identifier)
    missing = sorted(set(BUILD_IDENTIFIERS) - recorded)
    if manifest.get("complete") is not (not missing) or manifest.get("missing_targets") != missing:
        raise ValueError("Wheel manifest completeness differs from recorded targets")


def verify_installed_wheel(wheel):
    """Bind installed package checks to the wheel whose receipt is recorded."""
    distribution = importlib.metadata.distribution("smsd")
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if (name.endswith("/") or not (name.startswith("smsd/")
                    or name.endswith(".dist-info/METADATA"))):
                continue
            installed = distribution.locate_file(name)
            if not installed.is_file() or installed.read_bytes() != archive.read(name):
                raise ValueError("Installed package differs from the checked wheel: " + name)


def verify_local(args):
    source = checked_source(args.sdist, args.source_manifest)
    shared_validation = (checked_shared_validation(source, read_manifest(args.lean_validation))
                         if getattr(args, "lean_validation", None) else None)
    identifier = wheel_target(args.wheel.name, source["version"])
    target = BUILD_IDENTIFIERS[identifier]
    family = target["platform"]
    host = {"macos": "Darwin", "windows": "Windows", "linux": "Linux"}[family]
    architectures = ("arm64", "aarch64") if target["architecture"] in ("arm64", "aarch64") else ("amd64", "x86_64")
    if platform.system() != host or platform.machine().lower() not in architectures:
        raise ValueError("Local wheel checks require the wheel's native host architecture")
    if platform.python_implementation() != "CPython" or ".".join(map(str, sys.version_info[:2])) != target["python_version"]:
        raise ValueError("Local wheel checks require the wheel's CPython version")
    source_dir = args.source_dir.resolve()
    for relative, data in read_source_archive(args.sdist).items():
        path = source_dir / relative
        if not path.is_file() or path.read_bytes() != data:
            raise ValueError("Extracted SDK differs from the source distribution: " + relative)
    if shared_validation is None:
        if args.native_test_dir is None:
            raise ValueError("Full local checks require --native-test-dir")
        native_dir = args.native_test_dir.resolve()
        if native_dir.exists() and any(native_dir.iterdir()):
            raise ValueError("Native release checks require a new or empty build directory")
    jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "2")
    if not jobs.isdigit() or int(jobs) < 1:
        raise ValueError("CMAKE_BUILD_PARALLEL_LEVEL must be a positive integer")
    environment = os.environ.copy()
    for name in ("PYTHONPATH", "PYTHONHOME", "SMSD_DISABLE_GPU"):
        environment.pop(name, None)
    with tempfile.TemporaryDirectory(prefix="smsd-wheel-check-") as temporary:
        if shared_validation is None:
            run_native_debug_tests(source_dir, native_dir, jobs, environment, temporary)
        checks = [[sys.executable, str(source_dir / "scripts/check_python_wheel.py"),
                   "--source", str(source_dir), "--require-rdkit", "--require-openmp"]]
        if shared_validation is None:
            checks.append([sys.executable, "-m", "pytest", str(source_dir / "python/tests"), "-q",
                           "--import-mode=importlib"])
        for command in checks:
            subprocess.run(command, env=environment, cwd=temporary, check=True)
        verify_installed_wheel(args.wheel)
    record_build(source, args.wheel, args.wheel.parent, "local", identifier, shared_validation)


def record_build(source, wheel, output_dir, build_system, identifier, shared_validation=None):
    proof = {**source, **BUILD_IDENTIFIERS[identifier], "build_system": build_system, "build_identifier": identifier,
             "wheel": wheel.name, "wheel_sha256": sha256(wheel)}
    if shared_validation is None:
        proof.update(native_debug_tests="passed", installed_wheel_tests="passed")
    else:
        validation = checked_shared_validation(source, shared_validation)
        proof.update(validation_mode="lean", shared_validation_sha256=proof_digest(validation),
                     installed_wheel_smoke="passed")
        if build_system == "cibuildwheel":
            proof["platform_recipe_checks"] = "passed"
        elif build_system == "local":
            proof["installed_wheel_identity"] = "passed"
        else:
            raise ValueError("Expected a local or cibuildwheel validation path")
    if build_system == "cibuildwheel":
        proof["cibuildwheel"] = "4.2.1"
    write_manifest(output_dir / "wheel-build-provenance.json", proof)


def prepare(args):
    if not args.release_ref:
        raise ValueError("Source preparation requires --release-ref")
    resolved = subprocess.check_output(["git", "rev-parse", "--verify", args.release_ref + "^{commit}"], text=True).strip()
    commit = subprocess.check_output(["git", "rev-parse", "--verify", "HEAD"], text=True).strip()
    if not re.fullmatch(r"[0-9a-f]{40}", resolved) or resolved != commit:
        raise ValueError("Source checkout differs from the requested frozen release ref")
    if subprocess.check_output(["git", "status", "--porcelain"], text=True).strip():
        raise ValueError("Release-ref preparation requires a clean source checkout")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if list(args.output_dir.glob("*.tar.gz")):
        raise ValueError("Prepare output directory already contains a source distribution")
    command = [sys.executable, "-m", "build", "--sdist", "--outdir", str(args.output_dir)]
    if args.no_isolation:
        command.append("--no-isolation")
    subprocess.run(command, check=True)
    archives = list(args.output_dir.glob("smsd-*.tar.gz"))
    if len(archives) != 1:
        raise ValueError("Expected exactly one SMSD source distribution")
    manifest = inspect_source(archives[0])
    if args.release_ref.startswith("v") and args.release_ref != "v" + manifest["version"]:
        raise ValueError("Release tag and source version differ")
    verify_source_tree(archives[0], commit)
    if (subprocess.check_output(["git", "rev-parse", "--verify", "HEAD"], text=True).strip() != commit
            or subprocess.check_output(["git", "status", "--porcelain"], text=True).strip()):
        raise ValueError("Source checkout changed during source preparation")
    manifest.update(schema_version=1, release_ref=resolved, source_commit=commit, working_tree_clean=True)
    write_manifest(args.output_dir / "SOURCE_MANIFEST.json", manifest)
    print(json.dumps(manifest, indent=2))


def check_frozen_checkout(commit):
    head = subprocess.check_output(["git", "rev-parse", "--verify", "HEAD"], text=True).strip()
    if head != commit or subprocess.check_output(["git", "status", "--porcelain"], text=True).strip():
        raise ValueError("Release builds require a clean checkout of the frozen commit")


def check_source(args):
    """Check the canonical archive against a clean checkout of the frozen commit."""
    if not re.fullmatch(r"[0-9a-f]{40}", args.release_commit):
        raise ValueError("Release builds require a full frozen commit SHA")
    source = checked_source(args.sdist, args.source_manifest)
    if source["source_commit"] != args.release_commit:
        raise ValueError("Source distribution differs from the frozen release commit")
    check_frozen_checkout(args.release_commit)
    verify_source_tree(args.sdist, args.release_commit)
    check_frozen_checkout(args.release_commit)
    print(json.dumps(source, indent=2))


def build(args):
    source = checked_source(args.sdist, args.source_manifest)
    shared_validation = (checked_shared_validation(source, read_manifest(args.lean_validation))
                         if getattr(args, "lean_validation", None) else None)
    arch = args.arch or ARCHITECTURES[args.platform][0]
    if arch not in ARCHITECTURES[args.platform]:
        raise ValueError("Unsupported architecture for " + args.platform)
    if not args.dry_run and args.platform in ("windows", "macos"):
        expected_host = "Windows" if args.platform == "windows" else "Darwin"
        if platform.system() != expected_host:
            raise ValueError(args.platform + " builds require a native " + expected_host + " host")
    tag = {"linux": "manylinux_" + arch, "macos": "macosx_" + arch, "windows": "win_amd64"}[args.platform]
    python_version = args.python_version
    if python_version not in PYTHON_VERSIONS:
        raise ValueError("Unsupported release Python version")
    identifier = "cp" + python_version.replace(".", "") + "-" + tag
    env = os.environ.copy()
    # Keep caller cache/container settings, but enforce this release's build and test contract.
    for key in list(env):
        if key.startswith("CIBW_") and key not in ("CIBW_CACHE_PATH", "CIBW_CONTAINER_ENGINE"):
            del env[key]
    for key in ("CMAKE_ARGS", "SKBUILD_CMAKE_ARGS", "SKBUILD_CMAKE_DEFINE", "PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    if shared_validation is not None:
        env["CIBW_BEFORE_BUILD"] = ('python -m pip install --upgrade pip==26.2.2 && '
            'python -m pip install cmake==4.4.3 ninja==1.13.2')
        env["CIBW_TEST_COMMAND"] = ('python "{package}/scripts/check_python_wheel.py" '
            '--source "{package}" --require-rdkit --require-openmp')
    if args.manylinux_image:
        if args.platform != "linux":
            raise ValueError("--manylinux-image is only valid for Linux")
        env["CIBW_MANYLINUX_" + arch.upper() + "_IMAGE"] = args.manylinux_image
    command = [sys.executable, "-m", "cibuildwheel", str(args.sdist.resolve()),
               "--only", identifier, "--config-file", "{package}/scripts/cibuildwheel.toml",
               "--output-dir", str(args.output_dir.resolve())]
    print(json.dumps({"command": command, "source": source,
                      "manylinux_image": args.manylinux_image or ("manylinux_2_28" if args.platform == "linux" else None)}, indent=2))
    if args.dry_run:
        return
    if importlib.metadata.version("cibuildwheel") != "4.2.1":
        raise ValueError("Use cibuildwheel==4.2.1 for this release build")
    if args.output_dir.exists() and list(args.output_dir.glob("*.whl")):
        raise ValueError("Wheel output directory already contains wheels")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, env=env, check=True)
    wheels = list(args.output_dir.glob("*.whl"))
    if len(wheels) != 1:
        raise ValueError("Expected exactly one repaired, tested wheel")
    wheel = wheels[0]
    if wheel_target(wheel.name, source["version"]) != identifier:
        raise ValueError("Built wheel and requested target differ")
    with zipfile.ZipFile(wheel) as archive:
        metadata = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
        if len(metadata) != 1 or BytesParser().parsebytes(archive.read(metadata[0]))["Version"] != source["version"]:
            raise ValueError("Built wheel and source versions differ")
        if args.platform == "linux":
            for license_name in ("COPYING3", "COPYING.RUNTIME"):
                if not any(name.endswith("/" + license_name) for name in archive.namelist()):
                    raise ValueError("Linux wheel is missing libgomp licence notice " + license_name)
    record_build(source, wheel, args.output_dir, "cibuildwheel", identifier, shared_validation)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    source = commands.add_parser("prepare", help="Build and validate one source distribution")
    source.add_argument("--output-dir", type=Path, required=True)
    source.add_argument("--release-ref", required=True)
    source.add_argument("--no-isolation", action="store_true", help="Use already installed local build dependencies")
    source.set_defaults(function=prepare)
    frozen = commands.add_parser("check-source", help="Check the canonical sdist before release builds")
    frozen.add_argument("--sdist", type=Path, required=True)
    frozen.add_argument("--source-manifest", type=Path)
    frozen.add_argument("--release-commit", required=True)
    frozen.set_defaults(function=check_source)
    wheel = commands.add_parser("build", help="Repair and test a wheel from an existing sdist")
    wheel.add_argument("--sdist", type=Path, required=True)
    wheel.add_argument("--source-manifest", type=Path)
    wheel.add_argument("--platform", choices=("linux", "macos", "windows"), required=True)
    wheel.add_argument("--arch")
    wheel.add_argument("--python-version", choices=PYTHON_VERSIONS, default="3.14")
    wheel.add_argument("--manylinux-image")
    wheel.add_argument("--output-dir", type=Path, required=True)
    wheel.add_argument("--dry-run", action="store_true")
    wheel.add_argument("--lean-validation", type=Path,
                       help="Shared frozen-source checks; run platform recipes and installed smoke only")
    wheel.set_defaults(function=build)
    local = commands.add_parser("verify-local", help="Build Debug tests from the SDK and check the installed local wheel")
    local.add_argument("--sdist", type=Path, required=True)
    local.add_argument("--wheel", type=Path, required=True)
    local.add_argument("--source-manifest", type=Path)
    local.add_argument("--native-test-dir", type=Path)
    local.add_argument("--source-dir", type=Path, required=True)
    local.add_argument("--lean-validation", type=Path,
                       help="Shared frozen-source checks; run installed smoke without repeating full suites")
    local.set_defaults(function=verify_local)
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, OSError, KeyError, zipfile.BadZipFile, tarfile.TarError, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()

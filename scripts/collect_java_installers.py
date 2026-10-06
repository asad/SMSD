#!/usr/bin/env python3
"""Collect one tested MSI, DMG and DEB for the same Java CLI release."""

import argparse
import json
from pathlib import Path
from pathlib import PurePosixPath
import re
import shutil
import sys

from build_java_installer import JAVA_BUILD, JAVA_VERSION, digest, release_inputs, require


TARGETS = {"windows": ("amd64", "msi"), "macos": ("arm64", "dmg"), "linux": ("amd64", "deb")}
SOURCE_NAME = "OpenJDK25U-jdk-sources_25.0.4.1_1.tar.gz"
SOURCE_SHA256 = "cb9e50d2eb3de72ffd28466d7c511b4039b4bc4bd8f0d276d1e4cf81e30fdf09"
CLI_CASES = {"version", "help", "SMILES substructure", "SMARTS substructure", "aromatic substructure",
             "aromatic MCS", "SDF batch with UTF-8 path", "JSON export with UTF-8 path", "CML input", "PDB input"}
INSTALL_METHODS = {"windows": "native MSI install/run/remove", "macos": "DMG mount/copy/run/remove/detach",
                   "linux": "dpkg install/run/remove"}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def check_qa(qa, version):
    require(isinstance(qa, dict), "CLI check report is missing")
    checks = qa.get("checks")
    require(isinstance(checks, list) and len(checks) == len(CLI_CASES)
            and all(isinstance(case, dict) for case in checks), "Expected ten CLI checks")
    require({case.get("case") for case in checks} == CLI_CASES, "CLI cases are missing or duplicated")
    require(all(case.get("status") == "passed" for case in checks), "CLI checks contain a failure or limitation")
    require(type(qa.get("passed_cases")) is int and qa["passed_cases"] == len(CLI_CASES)
            and type(qa.get("known_limitations")) is int and qa["known_limitations"] == 0
            and qa.get("external_java_required") is False, "CLI check totals or standalone runtime status differ")
    cases = {case["case"]: case for case in checks}
    require(cases["version"].get("version") == version, "CLI version check differs")
    for name in ("SMILES substructure", "CML input", "PDB input"):
        require(cases[name].get("positive") is True and cases[name].get("negative") is False,
                "Positive/negative chemical checks differ: " + name)
    require(type(cases["SMILES substructure"].get("negative_exit_code")) is int
            and cases["SMILES substructure"]["negative_exit_code"] == 1, "Negative search exit code differs")
    for name in ("SMARTS substructure", "JSON export with UTF-8 path"):
        require(cases[name].get("exists") is True, "Chemical result or file export differs: " + name)
    for name in ("aromatic substructure", "aromatic MCS"):
        case = cases[name]
        require(type(case.get("mapped_atoms")) is int and case["mapped_atoms"] == 6
                and type(case.get("mapping_count")) is int and case["mapping_count"] > 0
                and case.get("mapping_injective") is True and case.get("ring_bonds_preserved") is True,
                "Aromatic mapping checks are incomplete: " + name)
    mcs = cases["aromatic MCS"]
    require(type(mcs.get("mcs_size")) is int and mcs["mcs_size"] == 6
            and mcs["mapping_count"] == 1 and mcs.get("fragment_exported") is True,
            "Aromatic MCS size or fragment export differs")
    batch = cases["SDF batch with UTF-8 path"]
    require(type(batch.get("target_count")) is int and batch["target_count"] == 2
            and batch.get("target_indices") == [0, 1]
            and all(type(index) is int for index in batch["target_indices"])
            and isinstance(batch.get("exists"), list) and len(batch["exists"]) == 2
            and batch["exists"][0] is True and batch["exists"][1] is False, "SDF batch results or indices differ")
    modules = qa.get("runtime_modules")
    require(isinstance(modules, list) and all(isinstance(name, str) for name in modules)
            and len(modules) == len(set(modules))
            and {"java.base", "java.xml", "java.desktop", "java.sql", "jdk.charsets", "jdk.unsupported"}.issubset(modules),
            "Standalone runtime modules are missing or duplicated")


def check_runtime(runtime, target):
    require(isinstance(runtime, dict) and runtime.get("IMPLEMENTOR") == "Eclipse Adoptium"
            and runtime.get("JAVA_VERSION") == JAVA_VERSION
            and runtime.get("JAVA_RUNTIME_VERSION") in (JAVA_BUILD, JAVA_BUILD + "-LTS")
            and runtime.get("IMPLEMENTOR_VERSION") == "Temurin-" + JAVA_BUILD,
            "Runtime differs from the pinned Temurin release")
    arch = ("aarch64",) if target == "macos" else ("amd64", "x86_64")
    require(runtime.get("OS_ARCH") in arch, "Runtime architecture differs")
    os_name = runtime.get("OS_NAME", "")
    require(isinstance(os_name, str) and ((target == "macos" and os_name == "Mac OS X")
            or (target == "linux" and os_name == "Linux")
            or (target == "windows" and os_name.startswith("Windows"))), "Runtime operating system differs")


def check_hashes(hashes, *, required):
    require(isinstance(hashes, dict) and (not required or bool(hashes)), "Runtime licence hashes are missing")
    for name, sha in hashes.items():
        require(isinstance(name, str) and not PurePosixPath(name).is_absolute()
                and ".." not in PurePosixPath(name).parts and "\\" not in name and ":" not in name
                and isinstance(sha, str) and re.fullmatch(r"[0-9a-f]{64}", sha) is not None,
                "Invalid runtime licence identity")
    if required:
        require("legal/java.base/LICENSE" in hashes, "Runtime GPL licence identity is missing")


def check_summary(record, inputs, path):
    require(isinstance(record, dict), "Collected installer record is missing")
    target = record.get("platform")
    require(target in TARGETS, "Unknown collected installer platform")
    architecture, extension = TARGETS[target]
    name = f"smsd-{inputs['version']}-{target}-{architecture}.{extension}"
    expected = {"version": inputs["version"], "platform": target, "architecture": architecture,
                "installer": name, "cli_jar_sha256": inputs["jar_sha256"], "runtime": JAVA_BUILD,
                "installation_method": INSTALL_METHODS[target], "cli_checks_passed": len(CLI_CASES),
                "installation": "passed", "removal": "passed"}
    require(all(record.get(key) == value for key, value in expected.items()), "Collected installer summary differs")
    require(path.name == name and digest(path) == record.get("installer_sha256")
            and type(record.get("size_bytes")) is int and path.stat().st_size == record["size_bytes"]
            and record["size_bytes"] > 0, "Collected installer bytes differ: " + name)
    check_runtime(record.get("runtime_details"), target)
    for field in ("build_qa", "installed_qa"):
        check_qa(record.get(field), inputs["version"])
    require(record["build_qa"] == record["installed_qa"], "Built and installed CLI evidence differs")
    launcher = record.get("launcher", {})
    require(isinstance(launcher, dict) and launcher.get("architecture") == architecture
            and launcher.get("format") == {"windows": "PE32+", "macos": "Mach-O64", "linux": "ELF64"}[target],
            "Collected native launcher identity differs")
    check_hashes(record.get("runtime_legal_hashes"), required=True)
    check_hashes(record.get("runtime_vendor_notice_hashes"), required=False)
    signing = record.get("signing", {})
    require(isinstance(signing, dict) and signing.get("installer") == "unsigned", "Installer signing status differs")
    environment = record.get("execution_environment")
    require(isinstance(environment, str) and environment.strip(), "Installer execution environment is missing")
    if target == "windows":
        require(launcher.get("subsystem") == 3 and record.get("authenticode_status") == "NotSigned",
                "Windows console or installer signing status differs")
        require(environment == "native Windows AMD64", "MSI needs native Windows execution")
        for field in ("install_exit_code", "remove_exit_code"):
            require(type(record.get(field)) is int and record[field] in (0, 3010), "MSI execution exit code differs")
    elif target == "macos":
        require(signing.get("application") == "ad hoc" and signing.get("developer_id") is False
                and signing.get("minimum_macos") == "11.0", "macOS signing or deployment target differs")


def check_report(directory, inputs):
    provenance = read_json(directory / "installer-provenance.json")
    installation = read_json(directory / "installation-check.json")
    require(isinstance(provenance, dict) and isinstance(installation, dict), "Installer reports must be JSON objects")
    target = provenance.get("platform")
    require(target in TARGETS, "Unknown installer platform")
    architecture, extension = TARGETS[target]
    name = f"smsd-{inputs['version']}-{target}-{architecture}.{extension}"
    path = directory / name
    require(path.is_file() and path.stat().st_size > 0, "Installer is missing: " + name)
    require(provenance.get("installer_size_bytes") == path.stat().st_size, "Installer size differs")
    expected = {"version": inputs["version"], "platform": target, "architecture": architecture,
                "installer": name, "installer_sha256": digest(path), "cli_jar_sha256": inputs["jar_sha256"]}
    for key, value in expected.items():
        require(provenance.get(key) == value and installation.get(key) == value,
                "Installer build/installation identity differs: " + key)
    require(installation.get("status") == "passed" and installation.get("install_verified") is True
            and installation.get("cleanup_verified") is True, "Installation/removal checks are incomplete")
    image = installation.get("installed_image", {})
    require(isinstance(image, dict), "Installed image report is missing")
    for key in ("version", "platform", "architecture", "cli_jar_sha256"):
        require(image.get(key) == expected[key], "Installed image identity differs: " + key)
    runtime = image.get("runtime", {})
    check_runtime(runtime, target)
    require(runtime == provenance.get("runtime") and image.get("launcher") == provenance.get("launcher"),
            "Built and installed runtime or launcher identity differs")
    require(image.get("runtime_legal_hashes") == provenance.get("runtime_legal_hashes")
            and image.get("runtime_legal_hashes"), "Installed runtime licences differ")
    require(image.get("runtime_vendor_notice_hashes") == provenance.get("runtime_vendor_notice_hashes"),
            "Installed runtime vendor notices differ")
    for report in (provenance, image):
        check_qa(report.get("qa"), inputs["version"])
    summary = {**expected, "size_bytes": path.stat().st_size, "runtime": JAVA_BUILD,
               "installation_method": installation["installation_method"],
               "execution_environment": installation["execution_environment"],
               "cli_checks_passed": image["qa"]["passed_cases"], "runtime_details": runtime,
               "launcher": image["launcher"], "build_qa": provenance["qa"], "installed_qa": image["qa"],
               "runtime_legal_hashes": image["runtime_legal_hashes"],
               "runtime_vendor_notice_hashes": image["runtime_vendor_notice_hashes"],
               "signing": provenance["signing"], "installation": "passed", "removal": "passed"}
    if target == "windows":
        for field in ("authenticode_status", "install_exit_code", "remove_exit_code"):
            summary[field] = installation.get(field)
    check_summary(summary, inputs, path)
    return path, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--installer-dir", type=Path, action="append", default=[])
    parser.add_argument("--runtime-source", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    try:
        inputs = release_inputs(args.release_dir)
        manifest_path = args.release_dir / "INSTALLER_MANIFEST.json"
        if args.check_only:
            require(not args.installer_dir and args.runtime_source is None, "Check-only does not collect files")
            manifest = read_json(manifest_path)
            require(isinstance(manifest, dict) and manifest.get("schema_version") == 2
                    and manifest.get("version") == inputs["version"] and manifest.get("cli_jar_sha256") == inputs["jar_sha256"],
                    "Installer manifest differs from the release")
            records = manifest.get("installers", [])
            require(isinstance(records, list) and len(records) == 3 and all(isinstance(item, dict) for item in records)
                    and {item.get("platform") for item in records} == set(TARGETS),
                    "Expected Windows, macOS and Linux installers")
            for record in records:
                path = args.release_dir / record["installer"]
                check_summary(record, inputs, path)
        else:
            require(len(args.installer_dir) == 3 and args.runtime_source is not None,
                    "Provide three tested installer directories and matching OpenJDK source")
            checked = [check_report(directory, inputs) for directory in args.installer_dir]
            require({record["platform"] for _, record in checked} == set(TARGETS), "Duplicate or missing installer platform")
            require(args.runtime_source.name == SOURCE_NAME and digest(args.runtime_source) == SOURCE_SHA256,
                    "OpenJDK source archive differs from the pinned runtime")
            # Validate the whole set before copying anything into the release.
            for path, _ in checked:
                destination = args.release_dir / path.name
                require(not destination.exists() or digest(destination) == digest(path), "Refusing to replace an installer")
            for path, _ in checked:
                shutil.copyfile(path, args.release_dir / path.name)
            shutil.copyfile(args.runtime_source, args.release_dir / SOURCE_NAME)
            manifest = {"schema_version": 2, "version": inputs["version"], "cli_jar_sha256": inputs["jar_sha256"],
                        "runtime_source": SOURCE_NAME, "runtime_source_sha256": SOURCE_SHA256,
                        "installers": sorted((record for _, record in checked), key=lambda item: item["platform"])}
            manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        require(manifest.get("runtime_source") == SOURCE_NAME and manifest.get("runtime_source_sha256") == SOURCE_SHA256
                and digest(args.release_dir / SOURCE_NAME) == SOURCE_SHA256, "Collected runtime source differs")
        print("Checked Java installers: Windows MSI, macOS DMG, Linux DEB")
    except (ValueError, OSError, KeyError, TypeError) as error:
        parser.exit(1, "Installer collection failed: " + str(error) + "\n")


if __name__ == "__main__":
    main()

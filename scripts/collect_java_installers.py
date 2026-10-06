#!/usr/bin/env python3
"""Collect one tested MSI, DMG and DEB for the same Java CLI release."""

import argparse
import json
from pathlib import Path
import shutil
import sys

from build_java_installer import JAVA_BUILD, JAVA_VERSION, digest, release_inputs, require


TARGETS = {"windows": ("amd64", "msi"), "macos": ("arm64", "dmg"), "linux": ("amd64", "deb")}
SOURCE_NAME = "OpenJDK25U-jdk-sources_25.0.4.1_1.tar.gz"
SOURCE_SHA256 = "cb9e50d2eb3de72ffd28466d7c511b4039b4bc4bd8f0d276d1e4cf81e30fdf09"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def check_report(directory, inputs):
    provenance = read_json(directory / "installer-provenance.json")
    installation = read_json(directory / "installation-check.json")
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
    for key in ("version", "platform", "architecture", "cli_jar_sha256"):
        require(image.get(key) == expected[key], "Installed image identity differs: " + key)
    runtime = image.get("runtime", {})
    require(runtime.get("IMPLEMENTOR") == "Eclipse Adoptium" and runtime.get("JAVA_VERSION") == JAVA_VERSION
            and runtime.get("JAVA_RUNTIME_VERSION") in (JAVA_BUILD, JAVA_BUILD + "-LTS"),
            "Installed runtime differs from the pinned Temurin release")
    require(image.get("runtime_legal_hashes") == provenance.get("runtime_legal_hashes")
            and image.get("runtime_legal_hashes"), "Installed runtime licences differ")
    require(image.get("runtime_vendor_notice_hashes") == provenance.get("runtime_vendor_notice_hashes"),
            "Installed runtime vendor notices differ")
    for report in (provenance, image):
        qa = report.get("qa", {})
        checks = qa.get("checks", [])
        require(checks and all(case.get("status") == "passed" for case in checks),
                "CLI checks contain a failure or known limitation")
        require(qa.get("known_limitations") == 0 and qa.get("external_java_required") is False,
                "Installer needs external Java or retains an input limitation")
        require({"CML input", "PDB input"}.issubset({case.get("case") for case in checks}),
                "CML/PDB installation checks are missing")
    if target == "windows":
        require(image.get("launcher", {}).get("subsystem") == 3, "Windows launcher must use the console")
        require(installation.get("authenticode_status") == "NotSigned", "Record MSI signing status explicitly")
    summary = {**expected, "size_bytes": path.stat().st_size, "runtime": JAVA_BUILD,
               "installation_method": installation["installation_method"],
               "execution_environment": installation["execution_environment"],
               "cli_checks_passed": image["qa"]["passed_cases"],
               "signing": provenance["signing"], "installation": "passed", "removal": "passed"}
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
            require(manifest.get("version") == inputs["version"] and manifest.get("cli_jar_sha256") == inputs["jar_sha256"],
                    "Installer manifest differs from the release")
            records = manifest.get("installers", [])
            require(len(records) == 3 and {item.get("platform") for item in records} == set(TARGETS),
                    "Expected Windows, macOS and Linux installers")
            for record in records:
                architecture, extension = TARGETS[record["platform"]]
                name = f"smsd-{inputs['version']}-{record['platform']}-{architecture}.{extension}"
                require(record.get("installer") == name and record.get("cli_jar_sha256") == inputs["jar_sha256"],
                        "Collected installer identity differs")
                path = args.release_dir / name
                require(digest(path) == record["installer_sha256"] and path.stat().st_size == record["size_bytes"],
                        "Collected installer bytes differ: " + name)
                require(record.get("installation") == record.get("removal") == "passed", "Installer execution is incomplete")
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
            manifest = {"version": inputs["version"], "cli_jar_sha256": inputs["jar_sha256"],
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

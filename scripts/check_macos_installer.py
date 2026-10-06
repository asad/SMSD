#!/usr/bin/env python3
"""Mount a release DMG, check its installed CLI, then remove and detach it."""

import argparse
import json
import os
from pathlib import Path
import platform
import plistlib
import re
import shutil
import subprocess
import tempfile

import build_java_installer as builder


MACH_O_MAGICS = {b"\xcf\xfa\xed\xfe", b"\xfe\xed\xfa\xcf", b"\xce\xfa\xed\xfe",
                 b"\xfe\xed\xfa\xce", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca",
                 b"\xca\xfe\xba\xbf", b"\xbf\xba\xfe\xca"}


def invoke(command, *, input_text=None):
    return subprocess.run([str(value) for value in command], input=input_text,
                          capture_output=True, text=True, encoding="utf-8", timeout=60)


def check_installer(release_dir, output_dir):
    _, target, architecture, _ = builder.host_target()
    builder.require(target == "macos" and architecture == "arm64", "Run this check on native macOS arm64")
    inputs = builder.release_inputs(release_dir)
    output_dir = output_dir.resolve()
    installers = sorted(output_dir.glob("*.dmg"))
    builder.require(len(installers) == 1, "Expected exactly one release DMG")
    installer = installers[0]
    provenance = json.loads((output_dir / "installer-provenance.json").read_text(encoding="utf-8"))
    builder.require(provenance.get("version") == inputs["version"], "Installer and release versions differ")
    builder.require(provenance.get("platform") == "macos" and provenance.get("architecture") == "arm64",
                    "Installer provenance has the wrong platform")
    builder.require(provenance.get("installer") == installer.name
                    and provenance.get("installer_sha256") == builder.digest(installer),
                    "Installer does not match its provenance")
    builder.require(provenance.get("cli_jar_sha256") == inputs["jar_sha256"],
                    "Installer provenance refers to another CLI JAR")

    work = Path(tempfile.mkdtemp(prefix="SMSD installer ")).resolve()
    mount = work / "Mounted image"
    mount.mkdir()
    copied_image = work / "Installation check é" / "SMSD.app"
    root = Path(__file__).resolve().parents[1]
    mounted = False
    report = {"version": inputs["version"], "platform": "macos", "architecture": "arm64",
              "installer": installer.name, "installer_sha256": builder.digest(installer),
              "cli_jar_sha256": inputs["jar_sha256"], "status": "failed",
              "installation_method": "DMG mount/copy/run/remove/detach",
              "install_verified": False, "cleanup_verified": False,
              "execution_environment": f"native macOS {platform.mac_ver()[0]} arm64",
              "quarantined_internet_download_tested": False, "minimum_os_execution_tested": False}

    def sanitise(value):
        for path, label in ((str(copied_image), "<installed-image>"), (str(work), "<installation-check>"),
                            (str(installer), installer.name), (str(output_dir), "<installer-output>"),
                            (str(release_dir.resolve()), "<release-directory>"),
                            (str(root), "<checkout>"), (str(Path.home()), "<user-home>")):
            value = value.replace(path, label)
        return re.sub(r"(?:/private)?/var/folders/[^\n]*", "<temporary-path>", value)

    error = None
    try:
        builder.require(not work.is_relative_to(root), "Installation checks must run outside the checkout")
        verified = invoke(["hdiutil", "verify", installer])
        builder.require(verified.returncode == 0, "DMG verification failed: " + sanitise(verified.stderr))
        report["dmg_verify"] = "passed"
        # A jpackage DMG displays the included Apache licence before its plist.
        attached = invoke(["hdiutil", "attach", "-readonly", "-nobrowse", "-plist",
                           "-mountpoint", mount, installer], input_text="Y\n")
        mounted = os.path.ismount(mount)
        builder.require(attached.returncode == 0, "DMG mount failed: " + sanitise(attached.stderr))
        start = attached.stdout.find("<?xml")
        builder.require(start >= 0, "DMG mount did not return a plist")
        entities = plistlib.loads(attached.stdout[start:].encode("utf-8"))["system-entities"]
        builder.require(any(item.get("mount-point") == str(mount) for item in entities),
                        "DMG was not mounted at the requested directory")
        builder.require(mounted, "Requested DMG mount is unavailable")
        report["mount_read_only"] = True
        source_image = mount / "SMSD.app"
        builder.require(source_image.is_dir() and not source_image.is_symlink(), "DMG application is missing")
        copied_image.parent.mkdir()
        copied = invoke(["ditto", source_image, copied_image])
        builder.require(copied.returncode == 0, "DMG application copy failed: " + sanitise(copied.stderr))
        report["installed_image"] = builder.verify_image(copied_image, release_dir)

        native_files = []
        for path in sorted(copied_image.rglob("*")):
            if not path.is_file() or path.is_symlink():
                continue
            with path.open("rb") as stream:
                magic = stream.read(4)
            if magic not in MACH_O_MAGICS:
                continue
            result = invoke(["lipo", "-archs", path])
            builder.require(result.returncode == 0 and result.stdout.strip() == "arm64",
                            "Installed native file does not target arm64: " + path.relative_to(copied_image).as_posix())
            native_files.append(path.relative_to(copied_image).as_posix())
        builder.require(native_files, "Installed application has no Mach-O binaries")
        report["native_files"] = native_files
        signature = invoke(["codesign", "--verify", "--deep", "--strict", copied_image])
        builder.require(signature.returncode == 0, "Installed signature check failed: " + sanitise(signature.stderr))
        description = invoke(["codesign", "--display", "--verbose=4", copied_image])
        builder.require(description.returncode == 0, "Installed signature details are unavailable")
        developer_id = "Authority=Developer ID Application" in description.stderr
        builder.require("Signature=adhoc" in description.stderr or developer_id,
                        "Installed application has an unexpected signing identity")
        report["signature_verify"] = "passed"
        report["signature_type"] = "Developer ID" if developer_id else "ad hoc"
        report["developer_id_signed"] = developer_id
        gatekeeper = invoke(["spctl", "--assess", "--verbose=2", "--type", "execute", copied_image])
        report["gatekeeper"] = {"accepted": gatekeeper.returncode == 0, "exit_code": gatekeeper.returncode,
                                "details": sanitise(gatekeeper.stdout + gatekeeper.stderr).strip()}
        report["notarised"] = gatekeeper.returncode == 0 and "source=Notarized Developer ID" in gatekeeper.stderr
        report["install_verified"] = True
    except Exception as exception:
        error = sanitise(str(exception))
    finally:
        mounted = mounted or os.path.ismount(mount)
        if copied_image.parent.exists():
            try:
                shutil.rmtree(copied_image.parent)
                report["installed_copy_removed"] = not copied_image.exists()
            except Exception as exception:
                error = (error + "; " if error else "") + sanitise(str(exception))
        if mounted:
            try:
                detached = invoke(["hdiutil", "detach", mount])
                if detached.returncode != 0:
                    detached = invoke(["hdiutil", "detach", "-force", mount])
                builder.require(detached.returncode == 0 and not os.path.ismount(mount),
                                "Test DMG did not detach: " + sanitise(detached.stderr))
                report["dmg_detach"] = "passed"
                mounted = False
            except Exception as exception:
                error = (error + "; " if error else "") + sanitise(str(exception))
        if not mounted and not os.path.ismount(mount):
            try:
                shutil.rmtree(work)
                report["cleanup_verified"] = not work.exists()
            except Exception as exception:
                error = (error + "; " if error else "") + sanitise(str(exception))

    if error:
        report["error"] = error
    elif report["install_verified"] and report["cleanup_verified"]:
        report["status"] = "passed"
    destination = output_dir / "installation-check.json"
    temporary_report = destination.with_suffix(".json.tmp")
    temporary_report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary_report.replace(destination)
    builder.require(report["status"] == "passed", "Mac installer check failed: " + (error or "Incomplete cleanup"))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    arguments = parser.parse_args()
    try:
        report = check_installer(arguments.release_dir, arguments.output_dir)
    except (OSError, ValueError, subprocess.SubprocessError) as exception:
        parser.exit(1, str(exception) + "\n")
    print(json.dumps({"installer": report["installer"], "status": report["status"],
                      "install_verified": report["install_verified"], "cleanup_verified": report["cleanup_verified"],
                      "developer_id_signed": report["developer_id_signed"], "notarised": report["notarised"]}, indent=2))


if __name__ == "__main__":
    main()

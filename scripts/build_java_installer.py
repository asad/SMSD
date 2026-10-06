#!/usr/bin/env python3
"""Package the validated Java CLI with a target-platform Java 25 runtime.

Build on the target operating system. Installers contain the published CLI
JAR; this helper does not rebuild Java sources or publish release assets.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import plistlib
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from urllib.parse import urlparse
import zipfile


MAIN_CLASS = "com.bioinception.smsd.cli.SMSDcli"
UPGRADE_UUID = "a8925984-4386-5a94-b379-a2b2140c5120"
JLINK_OPTIONS = "--bind-services --strip-debug --no-man-pages --no-header-files"
JAVA_VERSION = "25.0.4.1"
JAVA_BUILD = "25.0.4.1+1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def run(command, **kwargs):
    return subprocess.run([str(part) for part in command], check=True,
                          text=True, encoding="utf-8", **kwargs)


def release_inputs(directory):
    directory = directory.resolve()
    jars = sorted(directory.glob("smsd-*-jar-with-dependencies.jar"))
    require(len(jars) == 1, "Expected exactly one published CLI JAR")
    jar = jars[0]
    match = re.fullmatch(r"smsd-(\d+\.\d+\.\d+)-jar-with-dependencies\.jar", jar.name)
    require(match is not None, "Expected a three-part release version")
    checksums = {}
    for line in (directory / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        record = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        require(record is not None, "Invalid SHA256SUMS record")
        require(record[2] not in checksums, "Duplicate SHA256SUMS filename")
        checksums[record[2]] = record[1]
    require(checksums.get(jar.name) == digest(jar), "CLI JAR checksum does not match the release")
    with zipfile.ZipFile(jar) as archive:
        manifest = archive.read("META-INF/MANIFEST.MF").decode("utf-8")
        require("Main-Class: " + MAIN_CLASS in manifest, "Unexpected CLI main class")
        require("Java-Version: 25" in manifest, "CLI JAR must target Java 25")
        legal = {name: archive.read("META-INF/smsd/" + name) for name in ("LICENSE", "NOTICE")}
    root = Path(__file__).resolve().parents[1]
    for name, data in legal.items():
        require((root / name).read_bytes() == data, "Published JAR and root " + name + " differ")
    return {"version": match[1], "jar": jar, "jar_sha256": checksums[jar.name], "legal": legal}


def host_target():
    machine = platform.machine().lower()
    if sys.platform == "darwin":
        require(machine in ("arm64", "aarch64"), "The macOS release targets arm64")
        return "dmg", "macos", "arm64", "SMSD"
    if sys.platform == "win32":
        require(machine in ("amd64", "x86_64"), "The Windows release targets AMD64")
        return "msi", "windows", "amd64", "SMSD"
    require(sys.platform.startswith("linux"), "Unsupported installer operating system")
    require(machine in ("amd64", "x86_64"), "The Linux release targets AMD64")
    return "deb", "linux", "amd64", "smsd"


def java_release(path, java_executable=None):
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r'([A-Z_]+)="(.*)"', line)
        if match:
            result[match[1]] = match[2]
    if java_executable is not None:
        # jlink retains only JAVA_VERSION and MODULES in the release file.
        # Query the bundled JVM for its actual vendor, build and architecture.
        environment = os.environ.copy()
        for name in ("JAVA_HOME", "JRE_HOME", "JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS", "CLASSPATH"):
            environment.pop(name, None)
        output = run([java_executable, "-XshowSettings:properties", "-version"],
                     capture_output=True, env=environment, timeout=30).stderr
        properties = dict(re.findall(r"^\s+([a-z.]+) = (.+)$", output, re.M))
        require(result.get("JAVA_VERSION") == properties.get("java.version"), "Runtime version metadata differs from its JVM")
        for field, property_name in {"IMPLEMENTOR": "java.vendor", "IMPLEMENTOR_VERSION": "java.vendor.version",
                                     "JAVA_RUNTIME_VERSION": "java.runtime.version", "OS_ARCH": "os.arch",
                                     "OS_NAME": "os.name"}.items():
            result[field] = properties.get(property_name, "")
    require(result.get("JAVA_VERSION") == JAVA_VERSION, "Use the pinned Java " + JAVA_VERSION + " release")
    require(result.get("JAVA_RUNTIME_VERSION") in (JAVA_BUILD, JAVA_BUILD + "-LTS"),
            "Use the pinned Java build " + JAVA_BUILD)
    require(result.get("IMPLEMENTOR") == "Eclipse Adoptium", "Use the pinned Eclipse Temurin runtime")
    return result


def binary_architecture(path):
    data = path.read_bytes()
    if data[:4] == b"\x7fELF":
        require(data[4:6] == b"\x02\x01", "Expected ELF64 little endian")
        require(struct.unpack_from("<H", data, 18)[0] == 62, "Expected ELF AMD64")
        return {"format": "ELF64", "architecture": "amd64"}
    if data[:4] == b"\xcf\xfa\xed\xfe":
        require(struct.unpack_from("<I", data, 4)[0] == 0x0100000C, "Expected Mach-O arm64")
        return {"format": "Mach-O64", "architecture": "arm64"}
    require(data[:2] == b"MZ", "Unknown native executable format")
    offset = struct.unpack_from("<I", data, 0x3C)[0]
    require(data[offset:offset + 4] == b"PE\0\0", "Invalid PE signature")
    require(struct.unpack_from("<H", data, offset + 4)[0] == 0x8664, "Expected PE AMD64")
    require(struct.unpack_from("<H", data, offset + 24)[0] == 0x20B, "Expected PE32+")
    return {"format": "PE32+", "architecture": "amd64",
            "subsystem": struct.unpack_from("<H", data, offset + 24 + 68)[0]}


def image_paths(image, target):
    if target == "macos":
        return image / "Contents/app", image / "Contents/MacOS/SMSD", image / "Contents/runtime/Contents/Home"
    if target == "windows":
        return image / "app", image / "SMSD.exe", image / "runtime"
    return image / "lib/app", image / "bin/smsd", image / "lib/runtime"


def terminal_readme(version, target):
    launcher = {"macos": '"/Applications/SMSD.app/Contents/MacOS/SMSD"',
                "linux": "/opt/smsd/bin/smsd",
                "windows": '& "$env:LOCALAPPDATA\\SMSD\\SMSD.exe"'}[target]
    cml_note = ("The published 7.2.1 CLI has a known CML/PDB reader limitation. Use SMILES, MOL or SDF input.\n"
                if version == "7.2.1" else "")
    return (f"SMSD {version}: Java command-line application\n\n"
            "This installer includes Java 25. Run SMSD from a terminal; it has no graphical interface.\n"
            "Python packages are installed separately with pip.\n\n"
            f"{launcher} --help\n{launcher} --version\n"
            f'{launcher} --Q SMI --q CC --T SMI --t CCC --mode sub --json -\n\n'
            "The Windows per-user installation directory can be selected during installation.\n"
            + cml_note +
            "Runtime licences are retained in the bundled runtime legal directory.\n"
            f"Matching OpenJDK source: https://github.com/asad/SMSD/releases/download/v{version}/"
            "OpenJDK25U-jdk-sources_25.0.4.1_1.tar.gz\n"
            "https://github.com/asad/SMSD\n")


def molecule_record(name, elements, bonds):
    lines = [name, "  SMSD", "", f"{len(elements):3d}{len(bonds):3d}  0  0  0  0            999 V2000"]
    lines.extend(f"{float(index):10.4f}{0.0:10.4f}{0.0:10.4f} {element:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
                 for index, element in enumerate(elements))
    lines.extend(f"{first:3d}{second:3d}{order:3d}  0  0  0  0" for first, second, order in bonds)
    return "\n".join(lines + ["M  END", "$$$$", ""])


def smoke_cli(launcher, runtime, version):
    results = []
    with tempfile.TemporaryDirectory(prefix="smsd installer ") as temporary:
        work = Path(temporary) / "Molecule checks é"
        work.mkdir()
        environment = os.environ.copy()
        for name in ("JAVA_HOME", "JRE_HOME", "JAVA_TOOL_OPTIONS", "_JAVA_OPTIONS", "JDK_JAVA_OPTIONS", "CLASSPATH"):
            environment.pop(name, None)
        environment["PATH"] = str(work)

        def invoke(arguments, expected=0):
            completed = subprocess.run([str(launcher)] + arguments, cwd=work, env=environment,
                                       capture_output=True, text=True, encoding="utf-8", timeout=40)
            require(completed.returncode == expected,
                    "CLI check failed: " + " ".join(arguments[:2]) + "\n" + completed.stderr[:3000])
            return completed

        def search(query, target, *, query_type="SMI", target_type="SMI", mode="sub", expected=0):
            return json.loads(invoke(["--Q", query_type, "--q", query, "--T", target_type,
                                      "--t", target, "--mode", mode, "--json", "-"], expected).stdout)

        output = invoke(["--version"]).stdout
        require("SMSD Pro " + version in output, "Launcher version differs from the release")
        results.append({"case": "version", "status": "passed", "version": version})
        require("--Q" in invoke(["--help"]).stdout, "CLI help is incomplete")
        results.append({"case": "help", "status": "passed"})
        require(search("CC", "CCC").get("exists") is True, "Positive substructure result is wrong")
        require(search("N", "CCC", expected=1).get("exists") is False, "Negative substructure result is wrong")
        results.append({"case": "SMILES substructure", "status": "passed", "positive": True,
                        "negative": False, "negative_exit_code": 1})
        require(search("[#6]-[#6]", "CCC", query_type="SIG").get("exists") is True, "SMARTS result is wrong")
        results.append({"case": "SMARTS substructure", "status": "passed", "exists": True})
        target_edges = {frozenset((index, 1 if index == 6 else index + 1)) for index in range(1, 7)}

        def validate_ring(pairs):
            require(len(pairs) == 6, "Aromatic mapping must cover six atoms")
            require(len({pair[0] for pair in pairs}) == 6 and len({pair[1] for pair in pairs}) == 6,
                    "Aromatic mapping is not injective")
            require({pair[0] for pair in pairs} == set(range(6))
                    and {pair[1] for pair in pairs} == set(range(1, 7)), "Aromatic mapping contains incorrect elements")
            atom_map = dict(pairs)
            require(all(frozenset((atom_map[index], atom_map[(index + 1) % 6])) in target_edges
                        for index in range(6)), "Aromatic mapping does not preserve ring bonds")

        aromatic_sub = json.loads(invoke(["--Q", "SMI", "--q", "c1ccccc1", "--T", "SMI", "--t", "Oc1ccccc1",
                                         "--mode", "sub", "--mappings", "--json", "-"]).stdout)
        aromatic_maps = aromatic_sub.get("mappings", [])
        require(aromatic_maps, "Aromatic substructure mapping is missing")
        for item in aromatic_maps:
            validate_ring(item.get("pairs", []))
        results.append({"case": "aromatic substructure", "status": "passed", "mapped_atoms": 6,
                        "mapping_count": len(aromatic_maps), "mapping_injective": True,
                        "ring_bonds_preserved": True})
        mapping = search("c1ccccc1", "Oc1ccccc1", mode="mcs")
        pairs = mapping.get("pairs", [])
        require(mapping.get("mcs_size") == 6 and len(pairs) == 6, "Aromatic MCS size is wrong")
        validate_ring(pairs)
        require(mapping.get("mcs_smiles"), "MCS fragment export is missing")
        results.append({"case": "aromatic MCS", "status": "passed", "mcs_size": 6,
                        "mapped_atoms": 6, "mapping_count": 1, "mapping_injective": True,
                        "ring_bonds_preserved": True, "fragment_exported": True})
        sdf = work / "Targets é.sdf"
        sdf.write_text(molecule_record("Ethane", ["C", "C"], [(1, 2, 1)])
                       + molecule_record("Ammonia", ["N"], []), encoding="utf-8")
        batch = search("CC", str(sdf), target_type="SDF")
        rows = batch.get("results", [])
        require(batch.get("target_count") == 2 and len(rows) == 2, "SDF batch count is wrong")
        require([row.get("target_index") for row in rows] == [0, 1], "SDF batch indices changed")
        require([row.get("exists") for row in rows] == [True, False], "SDF batch results are wrong")
        require(not any("error" in row for row in rows), "SDF batch contains errors")
        results.append({"case": "SDF batch with UTF-8 path", "status": "passed", "target_count": 2,
                        "target_indices": [0, 1], "exists": [True, False]})
        export = work / "Search result é.json"
        invoke(["--Q", "SMI", "--q", "CC", "--T", "SMI", "--t", "CCC", "--json", str(export)])
        require(json.loads(export.read_text(encoding="utf-8")).get("exists") is True, "JSON file export failed")
        results.append({"case": "JSON export with UTF-8 path", "status": "passed", "exists": True})
        cml = work / "Ethane é.cml"
        cml.write_text('<?xml version="1.0"?><cml xmlns="http://www.xml-cml.org/schema"><molecule>'
                       '<atomArray><atom id="a1" elementType="C" hydrogenCount="3"/>'
                       '<atom id="a2" elementType="C" hydrogenCount="3"/></atomArray>'
                       '<bondArray><bond atomRefs2="a1 a2" order="1"/></bondArray></molecule></cml>',
                       encoding="utf-8")
        response = subprocess.run([str(launcher), "--Q", "SMI", "--q", "CC", "--T", "CML", "--t", str(cml),
                                   "--json", "-"], cwd=work, env=environment, capture_output=True,
                                  text=True, encoding="utf-8", timeout=40)
        if response.returncode == 0:
            require(json.loads(response.stdout).get("exists") is True, "CML result is wrong")
            require(search("N", str(cml), target_type="CML", expected=1).get("exists") is False,
                    "Negative CML substructure result is wrong")
            results.append({"case": "CML input", "status": "passed", "positive": True, "negative": False})
        else:
            require(version == "7.2.1" and response.returncode == 1
                    and "Only supported is reading of ChemFile objects." in response.stderr,
                    "Unexpected CML failure: " + response.stderr[:3000])
            results.append({"case": "CML input", "status": "known limitation",
                            "reason": "Published CLI passes an atom container to the ChemFile-only CML reader."})
        pdb = work / "Ethane é.pdb"
        pdb.write_text("HETATM    1  C1  ETH A   1       0.000   0.000   0.000  1.00  0.00           C  \n"
                       "HETATM    2  C2  ETH A   1       1.540   0.000   0.000  1.00  0.00           C  \n"
                       "CONECT    1    2\nCONECT    2    1\nEND\n", encoding="utf-8")
        response = subprocess.run([str(launcher), "--Q", "SMI", "--q", "CC", "--T", "PDB", "--t", str(pdb),
                                   "--json", "-"], cwd=work, env=environment, capture_output=True,
                                  text=True, encoding="utf-8", timeout=40)
        if response.returncode == 0:
            require(json.loads(response.stdout).get("exists") is True, "PDB result is wrong")
            require(search("N", str(pdb), target_type="PDB", expected=1).get("exists") is False,
                    "Negative PDB substructure result is wrong")
            results.append({"case": "PDB input", "status": "passed", "positive": True, "negative": False})
        else:
            require(version == "7.2.1" and response.returncode == 1
                    and "Only supported is reading of ChemFile objects." in response.stderr,
                    "Unexpected PDB failure: " + response.stderr[:3000])
            results.append({"case": "PDB input", "status": "known limitation",
                            "reason": "Published CLI passes an atom container to the ChemFile-only PDB reader."})
        java = runtime / "bin" / ("java.exe" if sys.platform == "win32" else "java")
        modules = run([java, "--list-modules"], cwd=work, env=environment, capture_output=True).stdout
        module_names = sorted(line.split("@")[0] for line in modules.splitlines())
        require({"java.base", "java.xml", "java.desktop", "jdk.charsets"}.issubset(module_names),
                "Bundled runtime is missing CLI or service provider modules")
    return {"checks": results, "passed_cases": sum(item["status"] == "passed" for item in results),
            "known_limitations": sum(item["status"] == "known limitation" for item in results),
            "external_java_required": False, "runtime_modules": module_names}


def verify_image(image, release_dir):
    inputs = release_inputs(release_dir)
    _, target, architecture, name = host_target()
    image = image.resolve()
    app, launcher, runtime = image_paths(image, target)
    require(image.is_dir() and launcher.is_file(), "Application image or launcher is missing")
    jars = sorted(app.rglob("*.jar"))
    require(jars == [app / inputs["jar"].name], "Image must contain exactly the published CLI JAR")
    require(not [path for path in image.rglob("*.jar")
                 if path not in jars and path != runtime / "lib/jrt-fs.jar"], "Unexpected JAR outside the application payload")
    require(digest(jars[0]) == inputs["jar_sha256"], "Image JAR differs from the release")
    for legal_name, content in inputs["legal"].items():
        require((app / legal_name).read_bytes() == content, "Image " + legal_name + " differs from the release")
    readme = app / "README.txt"
    require(readme.read_text(encoding="utf-8") == terminal_readme(inputs["version"], target),
            "Image terminal instructions differ from the expected content")
    expected = {inputs["jar"].name, "LICENSE", "NOTICE", "README.txt", name + ".cfg", ".jpackage.xml", ".package"}
    require(set(path.name for path in app.iterdir()).issubset(expected), "Unexpected application payload")
    require(all(path.is_file() for path in app.iterdir()), "Unexpected application payload directory")
    if (app / ".package").exists():
        require((app / ".package").read_bytes() == name.encode("utf-8"), "Unexpected installed-package marker")
    cfg = (app / (name + ".cfg")).read_text(encoding="utf-8")
    require(MAIN_CLASS in cfg and inputs["jar"].name in cfg, "Unexpected launcher configuration")
    java = runtime / "bin" / ("java.exe" if target == "windows" else "java")
    metadata = java_release(runtime / "release", java)
    expected_java_arch = "aarch64" if architecture == "arm64" else "x86_64"
    require(metadata.get("OS_ARCH") in (expected_java_arch, "amd64" if architecture == "amd64" else "aarch64"),
            "Runtime architecture differs from its target")
    native = binary_architecture(launcher)
    require(native["architecture"] == architecture, "Launcher architecture differs from its target")
    if target == "windows":
        require(native["subsystem"] == 3, "Windows needs a console launcher")
    require(binary_architecture(java)["architecture"] == architecture, "Bundled Java architecture differs")
    legal_dir = runtime / "legal"
    require(legal_dir.is_dir() and (legal_dir / "java.base/LICENSE").is_file(), "Runtime GPL licence is missing")
    legal_hashes = {str(path.relative_to(runtime)).replace(os.sep, "/"): digest(path)
                    for path in sorted(legal_dir.rglob("*")) if path.is_file()}
    require(legal_hashes, "Runtime licences are missing")
    require(all(path.is_relative_to(runtime) for path in (item.resolve() for item in legal_dir.rglob("*"))),
            "Runtime licence symlink escapes its image")
    signing = {"application": "not checked"}
    if target == "macos":
        plist = plistlib.loads((image / "Contents/Info.plist").read_bytes())
        require(plist.get("CFBundleShortVersionString") == inputs["version"], "macOS bundle version differs")
        require(plist.get("LSMinimumSystemVersion") == "11.0", "macOS runtime requires deployment target 11.0")
        run(["/usr/bin/codesign", "--verify", "--deep", "--strict", image], capture_output=True)
        signature = run(["/usr/bin/codesign", "--display", "--verbose=4", image], capture_output=True).stderr
        signing = {"application": "ad hoc" if "Signature=adhoc" in signature else "signed",
                   "developer_id": "Authority=Developer ID Application" in signature,
                   "notarisation": "not checked", "minimum_macos": "11.0"}
    qa = smoke_cli(launcher, runtime, inputs["version"])
    return {"version": inputs["version"], "platform": target, "architecture": architecture,
            "cli_jar": inputs["jar"].name, "cli_jar_sha256": inputs["jar_sha256"], "launcher": native,
            "runtime": {key: metadata[key] for key in ("IMPLEMENTOR", "IMPLEMENTOR_VERSION", "JAVA_VERSION",
                        "JAVA_RUNTIME_VERSION", "OS_ARCH", "OS_NAME", "SOURCE", "BUILD_SOURCE", "SOURCE_REPO",
                        "BUILD_SOURCE_REPO") if key in metadata},
            "runtime_legal_hashes": legal_hashes,
            "runtime_vendor_notice_hashes": {name: digest(runtime / name) for name in ("LICENSE", "NOTICE")
                                             if (runtime / name).is_file()},
            "signing": signing, "qa": qa}


def build(args):
    inputs = release_inputs(args.release_dir)
    package_format, target, architecture, name = host_target()
    require(args.format == package_format, "Build this installer on its target operating system")
    require(re.fullmatch(r"[0-9a-f]{64}", args.jdk_archive_sha256) is not None, "Invalid JDK archive SHA-256")
    source = urlparse(args.jdk_source_url)
    require(source.scheme == "https" and source.netloc == "github.com"
            and source.path.startswith("/adoptium/temurin25-binaries/releases/download/"),
            "Provide the pinned official Temurin 25 archive URL")
    java_home = args.java_home.resolve()
    metadata = java_release(java_home / "release")
    # Temurin 25 supports runtime linking without packaged JMODs. Exclude the
    # linker and jpackage, which depends on it, from the service provider set.
    limit_modules = sorted(set(metadata.get("MODULES", "").split()) - {"jdk.jlink", "jdk.jpackage"})
    require({"java.se", "jdk.unsupported"}.issubset(limit_modules), "Pinned JDK module metadata is incomplete")
    jlink_options = JLINK_OPTIONS + " --limit-modules " + ",".join(limit_modules)
    jpackage = java_home / "bin" / ("jpackage.exe" if target == "windows" else "jpackage")
    require(jpackage.is_file(), "JDK jpackage tool is missing")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    stage = out / "input"
    images = out / "app-image"
    packages = out / "package-output"
    for directory in (stage, images, packages):
        require(not directory.exists() or not any(directory.iterdir()), "Build directory is not empty: " + directory.name)
        directory.mkdir(exist_ok=True)
    shutil.copyfile(inputs["jar"], stage / inputs["jar"].name)
    for legal_name, data in inputs["legal"].items():
        (stage / legal_name).write_bytes(data)
    (stage / "README.txt").write_text(terminal_readme(inputs["version"], target), encoding="utf-8")
    command = [jpackage, "--type", "app-image", "--input", stage, "--dest", images, "--name", name,
               "--app-version", inputs["version"], "--main-jar", inputs["jar"].name,
               "--main-class", MAIN_CLASS, "--vendor", "BioInception PVT LTD",
               "--description", "Terminal application for chemical substructure and MCS search",
               "--copyright", "Copyright 2018-2026 Syed Asad Rahman - BioInception PVT LTD",
               "--add-modules", "java.se,jdk.unsupported", "--jlink-options", jlink_options,
               "--java-options", "-Xmx512m",
               "--java-options", "--add-opens=java.base/java.lang=ALL-UNNAMED"]
    if args.icon:
        require(args.icon.is_file(), "Installer icon is missing")
        command.extend(["--icon", args.icon.resolve()])
    if target == "windows":
        command.append("--win-console")
    elif target == "macos":
        command.extend(["--mac-package-identifier", "com.bioinceptionlabs.smsd"])
    run(command)
    image = images / (name + ".app" if target == "macos" else name)
    _, _, runtime = image_paths(image, target)
    for vendor_legal in ("NOTICE", "LICENSE"):
        source_legal = java_home / vendor_legal
        if source_legal.is_file():
            shutil.copyfile(source_legal, runtime / vendor_legal)
    if target == "macos":
        plist_path = image / "Contents/Info.plist"
        plist = plistlib.loads(plist_path.read_bytes())
        plist["LSMinimumSystemVersion"] = "11.0"
        plist_path.write_bytes(plistlib.dumps(plist, sort_keys=False))
        run(["/usr/bin/codesign", "--force", "--deep", "--sign", "-", image])
    verification = verify_image(image, args.release_dir)
    command = [jpackage, "--type", args.format, "--app-image", image, "--dest", packages,
               "--name", name, "--app-version", inputs["version"], "--vendor", "BioInception PVT LTD",
               "--description", "Terminal application for chemical substructure and MCS search",
               "--license-file", stage / "LICENSE", "--about-url", "https://github.com/asad/SMSD"]
    if target == "windows":
        command.extend(["--win-per-user-install", "--win-dir-chooser", "--win-upgrade-uuid", UPGRADE_UUID])
    elif target == "linux":
        command.extend(["--linux-package-name", "smsd", "--linux-deb-maintainer", "asad.rahman@bioinceptionlabs.com",
                        "--linux-app-category", "science", "--install-dir", "/opt"])
    run(command)
    installers = sorted(packages.glob("*." + args.format))
    require(len(installers) == 1, "Expected exactly one installer")
    final = out / f"smsd-{inputs['version']}-{target}-{architecture}.{args.format}"
    require(not final.exists(), "Installer output already exists")
    shutil.copyfile(installers[0], final)
    verification.update({"installer": final.name, "installer_sha256": digest(final),
                         "installer_size_bytes": final.stat().st_size, "jdk_archive_sha256": args.jdk_archive_sha256,
                         "jdk_archive_url": args.jdk_source_url, "jlink_options": jlink_options.split(),
                         "runtime_root_modules": ["java.se", "jdk.unsupported"],
                         "excluded_packaging_modules": ["jdk.jlink", "jdk.jpackage"],
                         "build_host": {"system": platform.system(), "release": platform.release(),
                                        "architecture": platform.machine()},
                         "package_installation": "pending target-platform installation check",
                         "jdk_version": metadata["JAVA_RUNTIME_VERSION"],
                         "jdk_source_build": {key: metadata[key] for key in ("SOURCE", "SOURCE_REPO", "BUILD_SOURCE",
                                              "BUILD_SOURCE_REPO") if key in metadata}})
    verification["signing"]["installer"] = "unsigned"
    (out / "installer-provenance.json").write_text(json.dumps(verification, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"installer": final.name, "sha256": digest(final),
                      "passed_cli_cases": verification["qa"]["passed_cases"],
                      "known_limitations": verification["qa"]["known_limitations"]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("build", help="Build and check a target-platform installer")
    prepare.add_argument("--release-dir", type=Path, required=True)
    prepare.add_argument("--format", choices=("dmg", "deb", "msi"), required=True)
    prepare.add_argument("--output-dir", type=Path, required=True)
    prepare.add_argument("--java-home", type=Path, required=True)
    prepare.add_argument("--jdk-archive-sha256", required=True)
    prepare.add_argument("--jdk-source-url", required=True, help="Official URL of the verified JDK archive")
    prepare.add_argument("--icon", type=Path)
    check = commands.add_parser("verify-image", help="Check an extracted or installed application image")
    check.add_argument("--image", type=Path, required=True)
    check.add_argument("--release-dir", type=Path, required=True)
    check.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "build":
            build(args)
        else:
            result = verify_image(args.image, args.release_dir)
            output = json.dumps(result, indent=2) + "\n"
            if args.output_json:
                args.output_json.parent.mkdir(parents=True, exist_ok=True)
                args.output_json.write_text(output, encoding="utf-8")
            print(output, end="")
    except (ValueError, OSError, subprocess.SubprocessError, zipfile.BadZipFile) as error:
        parser.exit(1, "Installer validation failed: " + str(error) + "\n")


if __name__ == "__main__":
    main()

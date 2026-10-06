# Java CLI installers

Each release provides one installer per platform, with Java 25 included:

The 7.2.2 JAR targets Java 8 or later. Installers include Java 25 LTS so users
do not need a separate Java installation.

| Platform | Package | Architecture |
|---|---|---|
| Windows | MSI | x86_64 |
| macOS | DMG | arm64 |
| Debian/Ubuntu Linux | DEB | x86_64 |

These packages run SMSD from a terminal. Python wheels are installed separately.
Verify downloads against `SHA256SUMS`. macOS and Windows installers are unsigned.

## Install and run

Windows: open the MSI and select the per-user installation directory. From PowerShell,
using the default directory:

```powershell
& "$env:LOCALAPPDATA\SMSD\SMSD.exe" --version
& "$env:LOCALAPPDATA\SMSD\SMSD.exe" --Q SMI --q CC --T SMI --t CCC --json -
```

macOS: open the DMG and copy `SMSD.app` to Applications. From Terminal:

```bash
"/Applications/SMSD.app/Contents/MacOS/SMSD" --version
"/Applications/SMSD.app/Contents/MacOS/SMSD" --Q SMI --q CC --T SMI --t CCC --json -
```

The application is not notarised. The bundled Java runtime requires macOS 11
or later; execution checks use the OS version listed in the release validation.
Python wheels have a separate macOS 26 minimum.

Debian/Ubuntu, from the directory containing the downloaded package:

```bash
sudo apt install ./smsd-7.2.2-linux-amd64.deb
/opt/smsd/bin/smsd --version
/opt/smsd/bin/smsd --Q SMI --q CC --T SMI --t CCC --json -
```

Remove the Windows package through Installed apps, delete `SMSD.app` on macOS,
or run `sudo apt remove smsd` on Linux. Other Linux distributions can use the
portable Java package with Java 8 or later installed; Java 25 LTS is preferred.

## Prepare release packages

Build on the target OS. All three installers must contain the exact same
validated CLI JAR. The shared builder checks its release checksum, licences,
native architecture and chemical searches before packaging it.

Linux can be prepared locally through Docker, including on an arm64 host:

```bash
bash scripts/build-linux-installer.sh dist/release-7.2.2 \
  build/platform-release/7.2.2/installers/linux
```

On macOS, use the verified Temurin 25.0.4.1+1 arm64 JDK:

```bash
python scripts/build_java_installer.py build \
  --release-dir dist/release-7.2.2 --format dmg \
  --output-dir build/platform-release/7.2.2/installers/macos/package \
  --java-home "$SMSD_JAVA_HOME" \
  --jdk-archive-sha256 61979887f7506a24a57439ff99adb8b3a7fc89977d9cfe3b8984f58a981b7b9d \
  --jdk-source-url 'https://github.com/adoptium/temurin25-binaries/releases/download/jdk-25.0.4.1%2B1/OpenJDK25U-jdk_aarch64_mac_hotspot_25.0.4.1_1.tar.gz'
python scripts/check_macos_installer.py \
  --release-dir dist/release-7.2.2 \
  --output-dir build/platform-release/7.2.2/installers/macos/package
```

Set `SMSD_JAVA_HOME` to the extracted JDK's `Contents/Home` directory. Verify
the archive checksum before extracting it. The builder accepts an optional
platform-native icon with `--icon`.

The manual `installers.yml` workflow prepares the Windows MSI and optionally
the Python wheel together. It downloads the validated JAR and source archive
from an existing draft or published release, verifies checksums, then installs,
runs and removes the MSI. It collects build artifacts without publishing them.
macOS and Linux builds remain local.

Before publication, collect all three tested packages:

```bash
python scripts/collect_java_installers.py \
  --release-dir dist/release-7.2.2 \
  --installer-dir build/platform-release/7.2.2/installers/macos/package \
  --installer-dir build/platform-release/7.2.2/installers/linux/package \
  --installer-dir build/platform-release/7.2.2/installers/windows \
  --runtime-source OpenJDK25U-jdk-sources_25.0.4.1_1.tar.gz
python scripts/collect_java_installers.py --release-dir dist/release-7.2.2 --check-only
```

Collection requires matching JAR/runtime versions, successful installation and
removal reports, and passing CML/PDB checks. It includes the matching OpenJDK
source archive. Runtime licences stay in each installed runtime's `legal/`
directory; SMSD's LICENSE and NOTICE stay beside the CLI JAR.

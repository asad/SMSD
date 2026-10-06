#!/usr/bin/env bash
# Build and check the x86_64 DEB locally, including on an arm64 Docker host.
set -euo pipefail

if [[ "${1:-}" == --help ]]; then
  printf 'Usage: %s [RELEASE_DIR [OUTPUT_DIR]]\n' "$0"
  printf 'Build and verify an x86_64 DEB locally with Docker and bundled Java 25.\n'
  exit 0
fi
[[ $# -le 2 ]] || { printf 'Expected at most two directory arguments.\n' >&2; exit 2; }

SMSD_INSTALLER_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SMSD_INSTALLER_VERSION="$(python3 - "$SMSD_INSTALLER_ROOT/java/pom.xml" <<'PYVERSION'
import sys
import xml.etree.ElementTree as ET
print(ET.parse(sys.argv[1]).findtext('{http://maven.apache.org/POM/4.0.0}version'))
PYVERSION
)"
SMSD_INSTALLER_RELEASE="${1:-$SMSD_INSTALLER_ROOT/dist/release-$SMSD_INSTALLER_VERSION}"
SMSD_INSTALLER_OUTPUT="${2:-$SMSD_INSTALLER_ROOT/build/platform-release/$SMSD_INSTALLER_VERSION/installers/linux}"
SMSD_INSTALLER_HOST_ARCH="$(uname -m)"
SMSD_INSTALLER_HOST_OS="$(uname -s)"
SMSD_INSTALLER_RUNTIME_CACHE="${SMSD_INSTALLER_RUNTIME_CACHE:-$SMSD_INSTALLER_ROOT/build/installer-runtime-cache/linux}"
mkdir -p "$SMSD_INSTALLER_OUTPUT" "$SMSD_INSTALLER_RUNTIME_CACHE"
SMSD_INSTALLER_RELEASE="$(cd "$SMSD_INSTALLER_RELEASE" && pwd)"
SMSD_INSTALLER_OUTPUT="$(cd "$SMSD_INSTALLER_OUTPUT" && pwd)"
SMSD_INSTALLER_RUNTIME_CACHE="$(cd "$SMSD_INSTALLER_RUNTIME_CACHE" && pwd)"

# Pin the amd64 manifest: some legacy Docker builders ignore --platform when
# the base is a multi-platform index.
SMSD_INSTALLER_BASE="maven@sha256:f109669454c37f3ce5e5c92897873abf6115f9f9c7dabc912fe37a9d8194c674"
SMSD_INSTALLER_IMAGE="smsd-java-installer-linux:25"
SMSD_INSTALLER_JDK_NAME="OpenJDK25U-jdk_x64_linux_hotspot_25.0.4.1_1.tar.gz"
SMSD_INSTALLER_JDK_URL="https://github.com/adoptium/temurin25-binaries/releases/download/jdk-25.0.4.1%2B1/$SMSD_INSTALLER_JDK_NAME"
SMSD_INSTALLER_JDK_SHA="dbb698396d478e7fa2b1e50f4103324b2a99b90569ee27c33f2261f9215cf41e"
SMSD_INSTALLER_JDK_ARCHIVE="$SMSD_INSTALLER_RUNTIME_CACHE/$SMSD_INSTALLER_JDK_NAME"

if [[ ! -f "$SMSD_INSTALLER_JDK_ARCHIVE" ]]; then
  SMSD_INSTALLER_DOWNLOAD="$(mktemp "$SMSD_INSTALLER_RUNTIME_CACHE/jdk-download.XXXXXX")"
  trap 'rm -f "$SMSD_INSTALLER_DOWNLOAD"' EXIT
  curl --fail --location --retry 3 "$SMSD_INSTALLER_JDK_URL" \
    --output "$SMSD_INSTALLER_DOWNLOAD"
  printf '%s  %s\n' "$SMSD_INSTALLER_JDK_SHA" "$SMSD_INSTALLER_DOWNLOAD" | shasum -a 256 -c -
  mv "$SMSD_INSTALLER_DOWNLOAD" "$SMSD_INSTALLER_JDK_ARCHIVE"
  trap - EXIT
fi
printf '%s  %s\n' "$SMSD_INSTALLER_JDK_SHA" "$SMSD_INSTALLER_JDK_ARCHIVE" | shasum -a 256 -c -

docker build --platform linux/amd64 --tag "$SMSD_INSTALLER_IMAGE" - <<DOCKERFILE
FROM $SMSD_INSTALLER_BASE
RUN apt-get update && apt-get install -y --no-install-recommends fakeroot binutils python3 && rm -rf /var/lib/apt/lists/*
DOCKERFILE
[[ "$(docker image inspect "$SMSD_INSTALLER_IMAGE" --format '{{.Architecture}}')" == amd64 ]]
SMSD_INSTALLER_IMAGE_ID="$(docker image inspect "$SMSD_INSTALLER_IMAGE" --format '{{.Id}}')"

docker run --rm --interactive --platform linux/amd64 \
  --volume "$SMSD_INSTALLER_ROOT:/workspace:ro" \
  --volume "$SMSD_INSTALLER_RELEASE:/release:ro" \
  --volume "$SMSD_INSTALLER_OUTPUT:/output" \
  --volume "$SMSD_INSTALLER_RUNTIME_CACHE:/runtime-cache:ro" \
  --workdir /workspace \
  --env "SMSD_INSTALLER_BASE=$SMSD_INSTALLER_BASE" \
  --env "SMSD_INSTALLER_IMAGE_ID=$SMSD_INSTALLER_IMAGE_ID" \
  --env "SMSD_INSTALLER_JDK_NAME=$SMSD_INSTALLER_JDK_NAME" \
  --env "SMSD_INSTALLER_JDK_URL=$SMSD_INSTALLER_JDK_URL" \
  --env "SMSD_INSTALLER_JDK_SHA=$SMSD_INSTALLER_JDK_SHA" \
  --env "SMSD_INSTALLER_HOST_ARCH=$SMSD_INSTALLER_HOST_ARCH" \
  --env "SMSD_INSTALLER_HOST_OS=$SMSD_INSTALLER_HOST_OS" \
  "$SMSD_INSTALLER_IMAGE" bash -s <<'CONTAINER'
set -euo pipefail
[[ "$(uname -m)" == x86_64 ]]
[[ "$(dpkg --print-architecture)" == amd64 ]]
printf '%s  %s\n' "$SMSD_INSTALLER_JDK_SHA" "/runtime-cache/$SMSD_INSTALLER_JDK_NAME" | sha256sum -c -
# Keep JDK and package staging on Linux storage. Some macOS Docker mounts
# cannot preserve the relative symlinks in the JDK's legal directory.
SMSD_INSTALLER_RUNTIME_DIR="$(mktemp -d /tmp/smsd-installer-runtime.XXXXXX)"
SMSD_INSTALLER_PACKAGE_DIR="$(mktemp -d /tmp/smsd-installer-package.XXXXXX)"
tar -xzf "/runtime-cache/$SMSD_INSTALLER_JDK_NAME" -C "$SMSD_INSTALLER_RUNTIME_DIR"
export JAVA_HOME="$SMSD_INSTALLER_RUNTIME_DIR/jdk-25.0.4.1+1"
export PATH="$JAVA_HOME/bin:$PATH"
# XZ remains compatible with older dpkg versions; level 1 keeps local
# emulated packaging practical without changing the installed payload.
export DPKG_DEB_COMPRESSOR_TYPE=xz
export DPKG_DEB_COMPRESSOR_LEVEL=1
export DPKG_DEB_THREADS_MAX=2
java -version
python3 scripts/build_java_installer.py build --release-dir /release \
  --format deb --output-dir "$SMSD_INSTALLER_PACKAGE_DIR" --java-home "$JAVA_HOME" \
  --jdk-archive-sha256 "$SMSD_INSTALLER_JDK_SHA" --jdk-source-url "$SMSD_INSTALLER_JDK_URL"

mapfile -t SMSD_INSTALLER_PACKAGES < <(find "$SMSD_INSTALLER_PACKAGE_DIR" -maxdepth 1 -name '*.deb' -type f)
[[ "${#SMSD_INSTALLER_PACKAGES[@]}" == 1 ]]
mkdir -p /output/package
cp "${SMSD_INSTALLER_PACKAGES[0]}" /output/package/
cp "$SMSD_INSTALLER_PACKAGE_DIR/installer-provenance.json" /output/package/
SMSD_INSTALLER_DEB="/output/package/$(basename "${SMSD_INSTALLER_PACKAGES[0]}")"
dpkg-deb --field "$SMSD_INSTALLER_DEB" > /output/deb-control.txt
dpkg-deb --contents "$SMSD_INSTALLER_DEB" > /output/deb-contents.txt
SMSD_INSTALLER_EXTRACTED="$(mktemp -d /tmp/smsd-installer-extracted.XXXXXX)"
dpkg-deb --extract "$SMSD_INSTALLER_DEB" "$SMSD_INSTALLER_EXTRACTED"
python3 scripts/build_java_installer.py verify-image --image "$SMSD_INSTALLER_EXTRACTED/opt/smsd" \
  --release-dir /release --output-json /output/extracted-image-check.json
dpkg --install "$SMSD_INSTALLER_DEB"
[[ -x /opt/smsd/bin/smsd ]]

python3 scripts/build_java_installer.py verify-image --image /opt/smsd \
  --release-dir /release --output-json /output/installed-image-check.json

python3 - "$SMSD_INSTALLER_DEB" <<'PY'
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

output = Path('/output')
deb = Path(sys.argv[1])
proof = {
    'platform': 'linux',
    'architecture': platform.machine(),
    'host_architecture': os.environ['SMSD_INSTALLER_HOST_ARCH'],
    'execution': 'QEMU amd64 on macOS arm64' if os.environ['SMSD_INSTALLER_HOST_OS'] == 'Darwin' and os.environ['SMSD_INSTALLER_HOST_ARCH'] in ('arm64', 'aarch64') else 'local Docker under x86_64 emulation' if os.environ['SMSD_INSTALLER_HOST_ARCH'] in ('arm64', 'aarch64') else 'local Docker on x86_64',
    'os_release': Path('/etc/os-release').read_text(),
    'build_base': os.environ['SMSD_INSTALLER_BASE'],
    'build_image_id': os.environ['SMSD_INSTALLER_IMAGE_ID'],
    'jdk_url': os.environ['SMSD_INSTALLER_JDK_URL'],
    'jdk_archive_sha256': os.environ['SMSD_INSTALLER_JDK_SHA'],
    'java_version': subprocess.check_output([os.environ['JAVA_HOME'] + '/bin/java', '-version'], stderr=subprocess.STDOUT, text=True),
    'deb_compression': {key: os.environ[key] for key in ('DPKG_DEB_COMPRESSOR_TYPE', 'DPKG_DEB_COMPRESSOR_LEVEL', 'DPKG_DEB_THREADS_MAX')},
    'package': deb.name,
    'package_sha256': hashlib.sha256(deb.read_bytes()).hexdigest(),
    'installed_image_check': json.loads((output / 'installed-image-check.json').read_text()),
    'tool_packages': subprocess.check_output(['dpkg-query', '-W', 'fakeroot', 'binutils', 'dpkg', 'python3'], text=True),
}
(output / 'linux-local-qa.json').write_text(json.dumps(proof, indent=2) + '\n')
print('Installed DEB passes CLI and chemical graph checks.')
PY
dpkg --remove smsd
[[ ! -e /opt/smsd/bin/smsd ]]
[[ ! -e /opt/smsd/lib/app ]]
[[ ! -e /opt/smsd/lib/runtime ]]
python3 - <<'PY'
import json
from pathlib import Path
p = Path('/output/linux-local-qa.json')
proof = json.loads(p.read_text())
proof['uninstall_check'] = 'passed; installed launcher, application and runtime removed'
p.write_text(json.dumps(proof, indent=2) + '\n')
installed = proof['installed_image_check']
report = {
    'version': installed['version'],
    'platform': installed['platform'],
    'architecture': installed['architecture'],
    'installer': proof['package'],
    'installer_sha256': proof['package_sha256'],
    'cli_jar_sha256': installed['cli_jar_sha256'],
    'status': 'passed',
    'installation_method': 'dpkg install/run/remove',
    'install_verified': True,
    'cleanup_verified': True,
    'cleanup_removed_paths': ['/opt/smsd/bin/smsd', '/opt/smsd/lib/app', '/opt/smsd/lib/runtime'],
    'execution_environment': proof['execution'],
    'installed_image': installed,
}
Path('/output/installation-check.json').write_text(json.dumps(report, indent=2) + '\n')
Path('/output/package/installation-check.json').write_text(json.dumps(report, indent=2) + '\n')
provenance_path = Path('/output/package/installer-provenance.json')
provenance = json.loads(provenance_path.read_text())
provenance['package_installation'] = 'passed; DEB installed, executed and removed'
provenance['installation_method'] = report['installation_method']
provenance['execution_environment'] = report['execution_environment']
provenance_path.write_text(json.dumps(provenance, indent=2) + '\n')
PY
CONTAINER

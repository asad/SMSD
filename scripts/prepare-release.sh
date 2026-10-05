#!/usr/bin/env bash
# Build, test, and assemble release assets locally. Does not publish or tag.
set -euo pipefail

SMSD_RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SMSD_RELEASE_ROOT"
SMSD_RELEASE_PYTHON="${SMSD_RELEASE_PYTHON:-python3}"
"$SMSD_RELEASE_PYTHON" -c 'import sys; assert sys.prefix != sys.base_prefix, "Run with a dedicated Python virtual environment"'
SMSD_RELEASE_VERSION="$("$SMSD_RELEASE_PYTHON" -c 'import xml.etree.ElementTree as E; print(E.parse("pom.xml").findtext("{http://maven.apache.org/POM/4.0.0}version"))')"
SMSD_RELEASE_FINAL_DIR="$SMSD_RELEASE_ROOT/dist/release-$SMSD_RELEASE_VERSION"
SMSD_RELEASE_JOBS="${SMSD_RELEASE_JOBS:-4}"
if [[ "$(uname -s)" == "Darwin" ]]; then
  # The release wheel targets macOS 26; bundled libraries must support this target.
  export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-26.0}"
fi
mkdir -p build dist
SMSD_RELEASE_DIR="$(mktemp -d "$SMSD_RELEASE_ROOT/build/release-assets.XXXXXX")"
SMSD_RELEASE_SOURCE="$(mktemp -d "$SMSD_RELEASE_ROOT/build/release-source.XXXXXX")"
trap 'rm -rf "$SMSD_RELEASE_SOURCE" "$SMSD_RELEASE_DIR"' EXIT

# Include the slow chemistry/benchmark suites and generate library documentation.
mvn -B -Dslow.tests.exclude=nothing clean verify \
  org.apache.maven.plugins:maven-source-plugin:3.3.1:jar-no-fork \
  org.apache.maven.plugins:maven-javadoc-plugin:3.6.3:jar
src/scripts/smsd --version
target/appassembler/bin/smsd --version

# Debug keeps C++ assert-based checks enabled. GPU builds are validated separately.
cmake -S cpp -B build/release-preflight -DCMAKE_BUILD_TYPE=Debug \
  -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_CUDA=OFF -DSMSD_BUILD_METAL=OFF
cmake --build build/release-preflight --parallel "$SMSD_RELEASE_JOBS"
ctest --test-dir build/release-preflight --output-on-failure

# Build from the source distribution to verify that it contains all build inputs.
"$SMSD_RELEASE_PYTHON" -m build --sdist --no-isolation --outdir "$SMSD_RELEASE_DIR"
tar -xzf "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" -C "$SMSD_RELEASE_SOURCE"
"$SMSD_RELEASE_PYTHON" -m build --wheel --no-isolation \
  -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF \
  --outdir "$SMSD_RELEASE_DIR" "$SMSD_RELEASE_SOURCE/smsd-$SMSD_RELEASE_VERSION"
if [[ "$(uname -s)" == "Darwin" ]]; then
  # Bundle external libraries such as OpenMP and validate macOS deployment tags.
  "$SMSD_RELEASE_PYTHON" -m delocate.cmd.delocate_wheel "$SMSD_RELEASE_DIR"/smsd-*.whl
fi
shopt -s nullglob
SMSD_RELEASE_WHEELS=("$SMSD_RELEASE_DIR"/smsd-*.whl)
[[ "${#SMSD_RELEASE_WHEELS[@]}" == 1 ]] || { echo "Expected exactly one release wheel" >&2; exit 1; }
"$SMSD_RELEASE_PYTHON" -m twine check --strict \
  "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" "${SMSD_RELEASE_WHEELS[0]}"
"$SMSD_RELEASE_PYTHON" -m pip install --no-deps --force-reinstall "${SMSD_RELEASE_WHEELS[0]}"
"$SMSD_RELEASE_PYTHON" - "$SMSD_RELEASE_VERSION" <<'PY'
import importlib.metadata
from pathlib import Path
import sys
import smsd
assert importlib.metadata.version("smsd") == smsd.__version__ == sys.argv[1]
assert Path(smsd.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
assert len(smsd.parse_smiles("c1ccccc1")) == 6
print("Verified installed wheel:", smsd.__file__)
PY
"$SMSD_RELEASE_PYTHON" -m pytest python/tests -q --import-mode=importlib

cp "target/smsd-$SMSD_RELEASE_VERSION.jar" \
  "target/smsd-$SMSD_RELEASE_VERSION-jar-with-dependencies.jar" \
  "target/smsd-$SMSD_RELEASE_VERSION-sources.jar" \
  "target/smsd-$SMSD_RELEASE_VERSION-javadoc.jar" "$SMSD_RELEASE_DIR/"
COPYFILE_DISABLE=1 tar -czf "$SMSD_RELEASE_DIR/smsd-cpp-$SMSD_RELEASE_VERSION-headers.tar.gz" \
  LICENSE NOTICE -C cpp/include smsd
COPYFILE_DISABLE=1 tar -czf "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION-cli.tar.gz" \
  LICENSE NOTICE -C target/appassembler bin repo
cp docs/RELEASE_NOTES.md "$SMSD_RELEASE_DIR/RELEASE_NOTES.md"
if [[ -f "docs/VALIDATION_$SMSD_RELEASE_VERSION.md" ]]; then
  cp "docs/VALIDATION_$SMSD_RELEASE_VERSION.md" "$SMSD_RELEASE_DIR/VALIDATION.md"
fi
if [[ -f "benchmarks/RESULTS_$SMSD_RELEASE_VERSION.md" ]]; then
  cp "benchmarks/RESULTS_$SMSD_RELEASE_VERSION.md" "$SMSD_RELEASE_DIR/BENCHMARK_RESULTS.md"
fi
if [[ -n "${SMSD_BENCHMARK_ARCHIVE:-}" ]]; then
  [[ -f "$SMSD_BENCHMARK_ARCHIVE" ]] || { echo "Benchmark archive not found" >&2; exit 1; }
  cp "$SMSD_BENCHMARK_ARCHIVE" "$SMSD_RELEASE_DIR/"
fi
"$SMSD_RELEASE_PYTHON" - "$SMSD_RELEASE_DIR" <<'PY'
from hashlib import sha256
from pathlib import Path
import sys
directory = Path(sys.argv[1])
assets = sorted(p for p in directory.iterdir() if p.is_file() and p.name != "SHA256SUMS")
(directory / "SHA256SUMS").write_text("".join(f"{sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in assets))
PY
rm -rf "$SMSD_RELEASE_FINAL_DIR"
mv "$SMSD_RELEASE_DIR" "$SMSD_RELEASE_FINAL_DIR"
echo "Verified local release assets: $SMSD_RELEASE_FINAL_DIR"

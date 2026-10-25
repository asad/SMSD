#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
# Build, test, and assemble release assets locally. Does not publish or tag.
set -euo pipefail

SMSD_RELEASE_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$SMSD_RELEASE_ROOT"
SMSD_RELEASE_PYTHON="${SMSD_RELEASE_PYTHON:-python3}"
"$SMSD_RELEASE_PYTHON" -c 'import sys; assert sys.prefix != sys.base_prefix, "Run with a dedicated Python virtual environment"'
mvn -version | "$SMSD_RELEASE_PYTHON" -c 'import re, sys; text = sys.stdin.read(); assert re.search(r"Java version: 25(?:[.,]|\s)", text), "Use JDK 25 for release builds; set JAVA_HOME accordingly"'
SMSD_RELEASE_VERSION="$("$SMSD_RELEASE_PYTHON" -c 'import xml.etree.ElementTree as E; print(E.parse("java/pom.xml").findtext("{http://maven.apache.org/POM/4.0.0}version"))')"
SMSD_RELEASE_FINAL_DIR="${SMSD_RELEASE_OUTPUT_DIR:-$SMSD_RELEASE_ROOT/dist/release-$SMSD_RELEASE_VERSION}"
if [[ -e "$SMSD_RELEASE_FINAL_DIR" ]]; then
  echo 'Release output already exists; set SMSD_RELEASE_OUTPUT_DIR to a new directory.' >&2
  exit 1
fi
SMSD_RELEASE_COMMIT="${SMSD_RELEASE_COMMIT:?Set SMSD_RELEASE_COMMIT to the full frozen commit SHA}"
SMSD_RELEASE_SDIST="${SMSD_RELEASE_SDIST:?Set SMSD_RELEASE_SDIST to the canonical source archive}"
SMSD_RELEASE_SOURCE_MANIFEST="${SMSD_RELEASE_SOURCE_MANIFEST:-$(dirname "$SMSD_RELEASE_SDIST")/SOURCE_MANIFEST.json}"
"$SMSD_RELEASE_PYTHON" scripts/build_python_wheels.py check-source \
  --sdist "$SMSD_RELEASE_SDIST" --source-manifest "$SMSD_RELEASE_SOURCE_MANIFEST" \
  --release-commit "$SMSD_RELEASE_COMMIT"
"$SMSD_RELEASE_PYTHON" -c 'import sys; assert sys.version_info[:2] == (3, 14), "Use CPython 3.14 for release preparation"'
"$SMSD_RELEASE_PYTHON" -c 'import importlib.metadata as M; required = {"pip": "26.2.2", "build": "1.6.1", "twine": "7.0.0"}; assert all(M.version(name) == version for name, version in required.items()), "Use the pinned release build tools described in docs/RELEASING.md"'
SMSD_RELEASE_JOBS="${SMSD_RELEASE_JOBS:-4}"
SMSD_RELEASE_JAVA8_HOME="${SMSD_RELEASE_JAVA8_HOME:-}"
[[ -x "$SMSD_RELEASE_JAVA8_HOME/bin/java" ]] || {
  echo 'Set SMSD_RELEASE_JAVA8_HOME to a Java 8 JDK for compatibility checks.' >&2
  exit 1
}
"$SMSD_RELEASE_JAVA8_HOME/bin/java" -version 2>&1 | \
  "$SMSD_RELEASE_PYTHON" -c 'import sys; assert "1.8.0_" in sys.stdin.read().splitlines()[0], "Compatibility runtime must be Java 8"'
if [[ "$(uname -s)" == "Darwin" ]]; then
  # The release wheel targets macOS 26; bundled libraries must support this target.
  export MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-26.0}"
fi
mkdir -p build dist
SMSD_RELEASE_DIR="$(mktemp -d "$SMSD_RELEASE_ROOT/build/release-assets.XXXXXX")"
SMSD_RELEASE_SOURCE="$(mktemp -d "$SMSD_RELEASE_ROOT/build/release-source.XXXXXX")"
trap 'printf "Retained release build directory: %s\n" "$SMSD_RELEASE_SOURCE" >&2; if [[ -d "$SMSD_RELEASE_DIR" ]]; then printf "Retained release assets: %s\n" "$SMSD_RELEASE_DIR" >&2; fi' EXIT
cp "$SMSD_RELEASE_SDIST" "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz"
cp "$SMSD_RELEASE_SOURCE_MANIFEST" "$SMSD_RELEASE_DIR/SOURCE_MANIFEST.json"
"$SMSD_RELEASE_PYTHON" scripts/build_python_wheels.py check-source \
  --sdist "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" \
  --source-manifest "$SMSD_RELEASE_DIR/SOURCE_MANIFEST.json" --release-commit "$SMSD_RELEASE_COMMIT"

# Include the slow chemistry suites and build the sources and Javadoc JARs of the
# library (java/core) and CLI (java/cli) modules. The source and Javadoc plugin
# versions come from the parent POM's plugin management.
"$SMSD_RELEASE_PYTHON" -m pytest scripts/tests -q
mvn -f java/pom.xml -B -Dslow.tests.exclude=nothing clean verify \
  source:jar-no-fork javadoc:jar
# The test phase finds the library classes for the CLI tests in the reactor;
# nothing is recompiled, so Java 8 runs the classes built above.
mvn -f java/pom.xml -B -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_RELEASE_JAVA8_HOME/bin/java" test
# The CLI reads its version from the JAR manifest, so each launcher must report the release.
for SMSD_VERSION_LINE in \
  "$("$SMSD_RELEASE_JAVA8_HOME/bin/java" -jar "java/cli/target/smsd-cli-$SMSD_RELEASE_VERSION-jar-with-dependencies.jar" --version)" \
  "$(java/cli/src/scripts/smsd --version)" \
  "$(java/cli/target/appassembler/bin/smsd --version)"; do
  [[ "$SMSD_VERSION_LINE" == "SMSD Pro $SMSD_RELEASE_VERSION "* ]] || {
    echo "CLI reports '$SMSD_VERSION_LINE' rather than $SMSD_RELEASE_VERSION" >&2
    exit 1
  }
done

# Build from the source distribution to verify that it contains all build inputs.
tar -xzf "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" -C "$SMSD_RELEASE_SOURCE"
"$SMSD_RELEASE_PYTHON" -m build --wheel --no-isolation \
  -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF \
  --outdir "$SMSD_RELEASE_DIR" "$SMSD_RELEASE_SOURCE/smsd-$SMSD_RELEASE_VERSION"
mkdir "$SMSD_RELEASE_SOURCE/unrepaired-wheels"
if [[ "$(uname -s)" == "Darwin" ]]; then
  # Bundle external libraries such as OpenMP and validate macOS deployment tags.
  cp "$SMSD_RELEASE_DIR"/smsd-*.whl "$SMSD_RELEASE_SOURCE/unrepaired-wheels/"
  "$SMSD_RELEASE_PYTHON" -m delocate.cmd.delocate_wheel "$SMSD_RELEASE_DIR"/smsd-*.whl
elif [[ "$(uname -s)" == "Linux" ]]; then
  SMSD_RELEASE_REPAIRED="$(mktemp -d "$SMSD_RELEASE_SOURCE/repaired.XXXXXX")"
  "$SMSD_RELEASE_PYTHON" -m auditwheel repair --plat manylinux_2_28_x86_64 \
    --wheel-dir "$SMSD_RELEASE_REPAIRED" "$SMSD_RELEASE_DIR"/smsd-*.whl
  mv "$SMSD_RELEASE_DIR"/smsd-*.whl "$SMSD_RELEASE_SOURCE/unrepaired-wheels/"
  mv "$SMSD_RELEASE_REPAIRED"/smsd-*.whl "$SMSD_RELEASE_DIR/"
else
  cp "$SMSD_RELEASE_DIR"/smsd-*.whl "$SMSD_RELEASE_SOURCE/unrepaired-wheels/"
fi
shopt -s nullglob
SMSD_RELEASE_WHEELS=("$SMSD_RELEASE_DIR"/smsd-*.whl)
[[ "${#SMSD_RELEASE_WHEELS[@]}" == 1 ]] || { echo "Expected exactly one release wheel" >&2; exit 1; }
"$SMSD_RELEASE_PYTHON" -m twine check --strict \
  "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" "${SMSD_RELEASE_WHEELS[0]}"
"$SMSD_RELEASE_PYTHON" -m pip install --no-deps --force-reinstall "${SMSD_RELEASE_WHEELS[0]}"
CMAKE_BUILD_PARALLEL_LEVEL="$SMSD_RELEASE_JOBS" \
"$SMSD_RELEASE_PYTHON" scripts/build_python_wheels.py verify-local \
  --sdist "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" \
  --wheel "${SMSD_RELEASE_WHEELS[0]}" \
  --native-test-dir "$SMSD_RELEASE_SOURCE/native-tests" \
  --source-dir "$SMSD_RELEASE_SOURCE/smsd-$SMSD_RELEASE_VERSION"
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
cp "java/core/target/smsd-$SMSD_RELEASE_VERSION.jar" \
  "java/core/target/smsd-$SMSD_RELEASE_VERSION-sources.jar" \
  "java/core/target/smsd-$SMSD_RELEASE_VERSION-javadoc.jar" \
  "java/cli/target/smsd-cli-$SMSD_RELEASE_VERSION-jar-with-dependencies.jar" "$SMSD_RELEASE_DIR/"
COPYFILE_DISABLE=1 tar -czf "$SMSD_RELEASE_DIR/smsd-cpp-$SMSD_RELEASE_VERSION-headers.tar.gz" \
  LICENSE NOTICE -C cpp/include smsd
COPYFILE_DISABLE=1 tar -czf "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION-cli.tar.gz" \
  LICENSE NOTICE licenses/cdk licenses/beam licenses/vecmath -C java/cli/target/appassembler bin repo
cp docs/RELEASE_NOTES.md "$SMSD_RELEASE_DIR/RELEASE_NOTES.md"
if [[ -f "docs/VALIDATION_$SMSD_RELEASE_VERSION.md" ]]; then
  cp "docs/VALIDATION_$SMSD_RELEASE_VERSION.md" "$SMSD_RELEASE_DIR/VALIDATION.md"
fi
SMSD_BENCHMARK_REPORT="${SMSD_BENCHMARK_REPORT:-}"
if [[ -n "$SMSD_BENCHMARK_REPORT" ]]; then
  [[ -f "$SMSD_BENCHMARK_REPORT" ]] || { echo "Benchmark report not found" >&2; exit 1; }
  cp "$SMSD_BENCHMARK_REPORT" "$SMSD_RELEASE_DIR/BENCHMARK_RESULTS.md"
fi
if [[ -n "${SMSD_BENCHMARK_ARCHIVE:-}" ]]; then
  [[ -f "$SMSD_BENCHMARK_ARCHIVE" ]] || { echo "Benchmark archive not found" >&2; exit 1; }
  cp "$SMSD_BENCHMARK_ARCHIVE" "$SMSD_RELEASE_DIR/"
fi
"$SMSD_RELEASE_PYTHON" scripts/build_python_wheels.py check-source \
  --sdist "$SMSD_RELEASE_DIR/smsd-$SMSD_RELEASE_VERSION.tar.gz" \
  --source-manifest "$SMSD_RELEASE_DIR/SOURCE_MANIFEST.json" --release-commit "$SMSD_RELEASE_COMMIT"
SMSD_RELEASE_COLLECT_ARGS=(--release-dir "$SMSD_RELEASE_DIR")
if [[ -n "${SMSD_PLATFORM_WHEELS_DIR:-}" ]]; then
  SMSD_RELEASE_COLLECT_ARGS+=(--wheel-dir "$SMSD_PLATFORM_WHEELS_DIR")
fi
"$SMSD_RELEASE_PYTHON" scripts/collect-release-wheels.py "${SMSD_RELEASE_COLLECT_ARGS[@]}"
if [[ -e "$SMSD_RELEASE_FINAL_DIR" ]]; then
  echo 'Release output was created during the build; it has been retained.' >&2
  exit 1
fi
mv "$SMSD_RELEASE_DIR" "$SMSD_RELEASE_FINAL_DIR"
echo "Verified complete release wheel matrix and local assets: $SMSD_RELEASE_FINAL_DIR"

# How to Build and Run

## Requirements
- Java 25+ (JDK 25 recommended)
- Maven 3.9+

## Build

```bash
mvn -U clean package
```

This produces `target/smsd-7.1.2-jar-with-dependencies.jar` (fat JAR with all dependencies, including CDK 2.13).

## Run Tests

```bash
mvn clean test
```

## Run the CLI

```bash
java -jar target/smsd-*-jar-with-dependencies.jar \
  --Q SMI --q "CCN" \
  --T SMI --t "CCCNC" \
  -m --json - --json-pretty
```

## Docker

```bash
docker build -t smsd .
docker run --rm smsd --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

## Notes
- The test suite exercises substructure and MCS, including recursive SMARTS,
  adversarial edge cases, and large molecules. Enable normally excluded suites
  with `mvn -Dslow.tests.exclude=nothing test`.
- Source launchers are at `src/scripts/smsd`, `smsd.bat`, and `smsd.ps1`.
  Maven also generates an isolated distribution under `target/appassembler/`.

## Configure Python bindings directly with CMake

Use CMake 3.18+, a C++17 compiler, and Python 3.9+ with pybind11 installed.
Select the same interpreter used to install pybind11:

```bash
cmake -S cpp -B build/python \
  -DSMSD_BUILD_PYTHON=ON -DSMSD_BUILD_TESTS=OFF \
  -DPython_EXECUTABLE="$PWD/.venv/bin/python" \
  -Dpybind11_DIR="$(.venv/bin/python -m pybind11 --cmakedir)"
cmake --build build/python --parallel 4
```

Python discovery uses `FindPython` and its `Development.Module` component.
Use `Python_EXECUTABLE` with this capitalization when selecting an interpreter.

## Prepare release assets locally

Use a dedicated Python environment with CMake, a C++17 compiler, and Java 25:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install build scikit-build-core pybind11 pytest pytest-timeout
# macOS dependency repair:
.venv/bin/python -m pip install delocate
SMSD_RELEASE_PYTHON=.venv/bin/python scripts/prepare-release.sh
```

Artifacts are assembled under `dist/release-7.1.2/` after validation succeeds.
On macOS, set `MACOSX_DEPLOYMENT_TARGET` to the supported deployment version;
delocate verifies dependencies against it. RDKit is optional for interoperability
tests. GPU test builds are separate from the CPU preflight. GitHub workflows
run only when manually dispatched; publishing and tagging are separate steps.

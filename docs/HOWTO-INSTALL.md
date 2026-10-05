# How to Build and Run

This source checkout targets version 7.2.0 with changes under Unreleased. Build
the checkout to use them. Current measured results and benchmark coverage are in
[the benchmark report](../benchmarks/RESULTS_7.2.0.md).
The currently published Maven Central and PyPI packages are 7.1.1; GitHub has
the 7.1.2 release. The commands below build the proposed 7.2.0 source and do
not assume its packages have already been published.

## Requirements
- Java 25+ (JDK 25 recommended)
- Maven 3.9+

## Build

```bash
mvn -U clean package
```

This produces `target/smsd-7.2.0-jar-with-dependencies.jar` (fat JAR with all dependencies, including CDK 2.13).

## Run Tests

```bash
mvn clean test
```

The full local correctness run includes the normally excluded algorithm and
stress suites:

```bash
mvn -Dslow.tests.exclude=nothing clean verify
```

Opt-in corpus benchmarks are separate. Their bounded defaults are one-second
pair budgets, no warmup and one measured trial; checkpoints are flushed after
each result. They retain their documented individual chemistry policies:

```bash
mvn test -Dslow.tests.exclude=nothing \
  '-Dtest=BenchmarkSuiteTest*,ExternalBenchmarkTest*,JavaCdkVsSmsdBenchmarkTest' \
  -Dbenchmark=true -Dsmsd.benchmark=true \
  -Dsmsd.benchmark.timeoutMs=1000 -Dsmsd.benchmark.rounds=1 \
  -Dsmsd.benchmark.warmup=0 \
  -Dsmsd.benchmark.outputDir=build/local-benchmarks/java
```

An elapsed-budget crossing does not identify cancellation: MCS may return a
valid incumbent when its budget ends. CDK diagnostic timings do not have the
same cancellation control as SMSD and do not establish a speed ranking.

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
  -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF \
  -DPython_EXECUTABLE="$PWD/.venv/bin/python" \
  -Dpybind11_DIR="$(.venv/bin/python -m pybind11 --cmakedir)"
cmake --build build/python --parallel 4
```

Python discovery uses `FindPython` and its `Development.Module` component.
Use `Python_EXECUTABLE` with this capitalization when selecting an interpreter.
The core extension does not require RDKit. Install RDKit in the same Python
environment for `from_rdkit` and the high-level RDKit molecule wrappers. Their
returned indices refer to the original RDKit inputs; raw native bindings use
`MolGraph` indices. See [the Python guide](PYTHON.md).

Build and install a CPU wheel through the packaging backend with the same
interpreter:

```bash
.venv/bin/python -m pip install build scikit-build-core pybind11
.venv/bin/python -m build --wheel \
  -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
.venv/bin/python -m pip install dist/smsd-7.2.0-*.whl
```

Use a wheel matching your Python version, operating system and architecture.
`SMSD_BUILD_METAL` and `SMSD_BUILD_CUDA` also accept `AUTO` or `ON` for optional
backends. Their availability does not change the CPU/OpenMP execution of core
batch matching.

## Prepare release assets locally

Use a dedicated Python environment with CMake, a C++17 compiler, and Java 25:

```bash
python3.14 -m venv .venv-release
.venv-release/bin/python -m pip install build scikit-build-core pybind11 pytest pytest-timeout twine
# macOS dependency repair:
.venv-release/bin/python -m pip install delocate
SMSD_RELEASE_PYTHON=.venv-release/bin/python scripts/prepare-release.sh
```

Artifacts are assembled under `dist/release-7.2.0/` after validation succeeds.
The macOS release target defaults to 26.0; delocate verifies bundled libraries
against it. RDKit is optional for interoperability
tests. GPU test builds are separate from the CPU preflight. GitHub workflows
run only when manually dispatched; publishing and tagging are separate steps.

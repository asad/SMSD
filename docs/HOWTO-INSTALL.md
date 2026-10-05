# How to Build and Run

This source checkout targets version 7.2.0 with changes under Unreleased. Build
the checkout to use them. Current measured results and benchmark coverage are in
[the benchmark report](../benchmarks/RESULTS_7.2.0.md).
The currently published Maven Central and PyPI packages are 7.1.1; GitHub has
the [7.1.2 release](https://github.com/asad/SMSD/releases/tag/v7.1.2). The commands
below use the proposed 7.2.0 source or locally prepared assets and do not assume
its packages have already been published.

## Choose a distribution

The Java CLI JAR and portable CLI archive contain platform-independent Java
code. Use the same files on Linux, macOS or Windows with Java 25 installed for
your machine's architecture; they do not include a Java runtime. The archive
includes Unix and Windows launchers. Native DMG, MSI and DEB installers are a
separate manual workflow and are not required to run the portable CLI.

Python wheels contain native code and must match the operating system,
architecture and Python interpreter. The current 7.2.0 local validation covers
CPython 3.14 on macOS arm64, executed on macOS 27.0.1 with a macOS 26+ deployment
target. Linux x86_64 passes 691 Python tests with 8 skips under local emulation;
native Windows x86_64 validation is pending. Intel
macOS is outside this compact wheel set; use a source build. A wheel tagged `cp314` is for ordinary CPython
3.14, not the free-threaded `cp314t` interpreter. Use a source build when a
matching wheel is unavailable. C++ headers are also available for source
builds on each operating system.

## Install the portable Java CLI

Download or copy `smsd-7.2.0-jar-with-dependencies.jar` from the release assets
after they are published, or build it with Maven below. Keep its exact filename
in these commands; shells differ in how they expand JAR filename wildcards.

Linux or macOS, from the directory containing the JAR:

```bash
java -version
java -jar smsd-7.2.0-jar-with-dependencies.jar --version
java -jar smsd-7.2.0-jar-with-dependencies.jar --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

Windows PowerShell, with Java 25 on `PATH`:

```powershell
java -version
java -jar .\smsd-7.2.0-jar-with-dependencies.jar --version
java -jar .\smsd-7.2.0-jar-with-dependencies.jar --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

If using the portable archive instead, extract it into its own directory so
`bin` and `repo` remain together:

```bash
mkdir smsd-cli
tar -xzf smsd-7.2.0-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

Windows PowerShell with the system `tar` command:

```powershell
New-Item -ItemType Directory -Path smsd-cli
tar -xzf .\smsd-7.2.0-cli.tar.gz -C .\smsd-cli
.\smsd-cli\bin\smsd.bat --version
```

The Windows archive launcher uses `java` from `PATH`; confirm that it is Java
25 even when `JAVA_HOME` is set. Retain the launcher's CRLF line endings. On
Unix, retain its executable permission or invoke it with `sh`. Paths containing
spaces should be quoted. The direct JAR command works independently of the
archive launchers.

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
java -jar target/smsd-7.2.0-jar-with-dependencies.jar \
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

## Build C++ locally on Linux, macOS or Windows

Install CMake 3.18+ and a C++17 compiler: GCC or Clang on Linux, Apple Clang
from Xcode Command Line Tools on macOS, or MSVC from Visual Studio Build Tools
with the Desktop development with C++ workload on Windows. On Windows, use a
Developer PowerShell prompt. The configuration flags below cover both
single-configuration and Visual Studio generators. All 12 native suites pass
on macOS and emulated Linux; native Windows execution remains pending:

```text
cmake -S cpp -B build/cpu -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpu --config Debug --parallel 4
ctest --test-dir build/cpu --build-config Debug --output-on-failure
```

OpenMP is detected when available; otherwise batch processing uses the
sequential fallback. The core headers do not require RDKit. The optional C++
RDKit adapter additionally requires an RDKit development installation and
C++20. See [the C++ guide](CPP.md) for integration details.

## Install or build Python locally

Install ordinary CPython 3.14 for the compact release wheel set. Native source
builds require the C++ compiler described above and Python development headers.
Start in the source checkout. On Linux or macOS:

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install build scikit-build-core pybind11 cmake ninja
.venv/bin/python -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

For the macOS release target, set `export MACOSX_DEPLOYMENT_TARGET=26.0` before
building. Install `libomp` separately if OpenMP support is desired. A local
source build is distinct from a repaired redistributable wheel that bundles
its external libraries.

On Windows, use Developer PowerShell with CPython 3.14 installed:

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install build scikit-build-core pybind11 cmake ninja
.\.venv\Scripts\python.exe -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Install the exact wheel matching your environment, from the release downloads
or the `dist` directory created by the build. Replace the filename below with
that wheel's full filename. Linux or macOS:

```bash
.venv/bin/python -m pip install /path/to/matching-smsd-wheel.whl
.venv/bin/python -c 'import smsd; print(smsd.__version__); assert len(smsd.parse_smiles("c1ccccc1")) == 6'
```

Windows PowerShell:

```powershell
.\.venv\Scripts\python.exe -m pip install "C:\path\to\matching-smsd-wheel.whl"
.\.venv\Scripts\python.exe -c "import smsd; print(smsd.__version__); assert len(smsd.parse_smiles('c1ccccc1')) == 6"
```

RDKit is optional. Install it in the same environment for the high-level RDKit
molecule wrappers. For local tests, install `pytest` and `pytest-timeout`, then
run the selected virtual environment's Python with
`-m pytest python/tests -q --import-mode=importlib`.

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

The direct CMake example uses Linux/macOS virtual environment paths. On Windows,
use `.venv/Scripts/python.exe` for `Python_EXECUTABLE`, and obtain `pybind11_DIR`
by running that interpreter with `-m pybind11 --cmakedir`.
`SMSD_BUILD_METAL` and `SMSD_BUILD_CUDA` also accept `AUTO` or `ON` for optional
backends. Their availability does not change the CPU/OpenMP execution of core
batch matching.

## Prepare release assets locally

The Bash preflight below was validated locally on macOS. Use a dedicated Python
environment with CMake, a C++17 compiler, and Java 25:

```bash
python3.14 -m venv .venv-release
.venv-release/bin/python -m pip install build scikit-build-core pybind11 pytest pytest-timeout twine
# macOS dependency repair:
.venv-release/bin/python -m pip install delocate
SMSD_RELEASE_PYTHON=.venv-release/bin/python scripts/prepare-release.sh
```

Artifacts are assembled under `dist/release-7.2.0/` after validation succeeds.
This prepares assets for the current native platform; it does not cross-build
Linux, macOS and Windows wheels in one invocation.
The macOS release target defaults to 26.0; delocate verifies bundled libraries
against it. RDKit is optional for interoperability
tests. GPU test builds are separate from the CPU preflight. GitHub workflows
run only when manually dispatched; publishing and tagging are separate steps.

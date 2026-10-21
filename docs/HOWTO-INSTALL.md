# Installation and source builds

This guide covers the SMSD Pro 8.0.0 packages, their platform requirements and
source builds of the Java, C++ and Python modules. Release files are published
on the [GitHub releases page](https://github.com/asad/SMSD/releases); the
[README](../README.md) lists where each package is distributed.

## Platform requirements

| Package | Platform | Requirements |
|---|---|---|
| Java library, CLI JAR and portable CLI | Windows, macOS or Linux | Java 8 or later; Java 25 LTS preferred |
| Windows installer (MSI) | Windows, x86_64 (AMD64) | Includes Java 25.0.4.1+1 |
| macOS installer (DMG) | macOS 11.0 or later, arm64 | Includes Java 25.0.4.1+1 |
| Linux installer (DEB) | Debian or Ubuntu, x86_64 (AMD64) | Includes Java 25.0.4.1+1 |
| Python wheel | Windows 10 or later, x86_64 (AMD64) | CPython 3.14 |
| Python wheel | Linux with glibc 2.28 or later, x86_64 (AMD64) | CPython 3.14 |
| Python wheel | macOS 26 or later, arm64 | CPython 3.14 |
| Python source build | Any platform with a C++17 compiler | Python 3.9 or later; CMake 3.18 or later |
| C++ headers | Any platform | C++17 compiler; core matching needs neither RDKit nor CDK |

Minimum operating-system versions are build targets. The versions on which a
release was tested are recorded in its validation record; see the
[8.0.0 validation record](VALIDATION_8.0.0.md).

## Verify downloads

Each GitHub release includes a `SHA256SUMS` file. Download it into the same
directory as the release files, then check each file you intend to use.

macOS or Linux:

```bash
grep ' smsd-cli-8.0.0-jar-with-dependencies.jar$' SHA256SUMS | shasum -a 256 -c
```

Windows PowerShell (the result is `True` when the checksum matches):

```powershell
$entry = Select-String -Path SHA256SUMS -SimpleMatch ' smsd-8.0.0-windows-amd64.msi'
(Get-FileHash .\smsd-8.0.0-windows-amd64.msi -Algorithm SHA256).Hash -eq $entry.Line.Split(' ')[0]
```

Replace the file name to check another download.

## Portable Java CLI

The shaded JAR `smsd-cli-8.0.0-jar-with-dependencies.jar` includes the SMSD
library, CDK 2.13 and the other dependencies of the CLI, and runs on Java 8 or
later. From the directory that contains it:

```bash
java -version
java -jar smsd-cli-8.0.0-jar-with-dependencies.jar --version
java -jar smsd-cli-8.0.0-jar-with-dependencies.jar \
  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

The last command runs a substructure search and writes every mapping as JSON
to standard output. In Windows PowerShell, use `.\smsd-cli-8.0.0-jar-with-dependencies.jar`
as the JAR path. The same JAR is published on Maven Central as
`com.bioinceptionlabs:smsd-cli:8.0.0` with the classifier `jar-with-dependencies`.

The archive `smsd-8.0.0-cli.tar.gz` contains Unix and Windows launchers. Extract
it into one directory so that `bin` and `repo` remain together:

```bash
mkdir smsd-cli
tar -xzf smsd-8.0.0-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

On Windows, run `smsd-cli\bin\smsd.bat`. The portable packages do not include a
Java runtime.

## Java CLI installers

The MSI, DMG and DEB packages install the command-line application with a
bundled Eclipse Temurin Java 25.0.4.1+1 runtime. They have no graphical
interface, and Python packages are installed separately. The packaging scripts
produce unsigned installer packages; the macOS application carries an ad hoc
signature and is not notarised.

Windows (`smsd-8.0.0-windows-amd64.msi`): open the MSI and choose the per-user
installation directory. With the default directory, from PowerShell:

```powershell
& "$env:LOCALAPPDATA\SMSD\SMSD.exe" --version
& "$env:LOCALAPPDATA\SMSD\SMSD.exe" --Q SMI --q CC --T SMI --t CCC --json -
```

macOS (`smsd-8.0.0-macos-arm64.dmg`): open the DMG and copy `SMSD.app` to
Applications. From Terminal:

```bash
"/Applications/SMSD.app/Contents/MacOS/SMSD" --version
"/Applications/SMSD.app/Contents/MacOS/SMSD" --Q SMI --q CC --T SMI --t CCC --json -
```

If Gatekeeper blocks a download whose checksum you have verified, follow
Apple's guidance on [opening apps safely](https://support.apple.com/en-gb/102445).

Debian or Ubuntu (`smsd-8.0.0-linux-amd64.deb`), from the download directory:

```bash
sudo apt install ./smsd-8.0.0-linux-amd64.deb
/opt/smsd/bin/smsd --version
/opt/smsd/bin/smsd --Q SMI --q CC --T SMI --t CCC --json -
```

To uninstall, remove SMSD under Settings > Apps on Windows, delete `SMSD.app`
on macOS, or run `sudo apt remove smsd` on Linux. On other Linux distributions,
use the portable Java CLI.

## Java source build

Use a JDK (Java 25 LTS preferred) and Maven 3.9 or later. From the repository
root:

```bash
mvn -f java/pom.xml clean package
java -jar java/cli/target/smsd-cli-8.0.0-jar-with-dependencies.jar --version
```

`java/pom.xml` is the parent of the library module, `java/core`, and the CLI
module, `java/cli`. The build runs the default test suites of both, writes the
library JAR to `java/core/target`, and writes the CLI JAR, the shaded CLI JAR
and the portable launchers to `java/cli/target`. To use the library as a Maven
dependency, run `mvn -f java/pom.xml install`, which installs
`com.bioinceptionlabs:smsd:8.0.0` and `com.bioinceptionlabs:smsd-cli:8.0.0` in
the local Maven repository. The [Java module README](../java/README.md#installation)
gives the Maven and Gradle dependencies of both. The root `pom.xml` also builds
the Java modules.

The algorithm, stress and benchmark suites are excluded by default. To include
the algorithm and stress suites:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing clean verify
```

The code is compiled for Java 8 (`--release 8`). To run the tests on a Java 8
runtime, set `SMSD_JAVA8_HOME` to a Java 8 JDK:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" test
```

The corpus benchmarks are opt-in and separate from the correctness tests:

```bash
mvn -f java/pom.xml -pl core test -Dslow.tests.exclude=nothing \
  '-Dtest=BenchmarkSuiteTest*,ExternalBenchmarkTest*,JavaCdkVsSmsdBenchmarkTest' \
  -Dbenchmark=true -Dsmsd.benchmark=true
```

The system properties `smsd.benchmark.timeoutMillis` (default `1000`),
`smsd.benchmark.rounds` (default `1`) and `smsd.benchmark.warmup` (default `0`)
set the protocol. Results are written to `build/local-benchmarks/java` unless
`smsd.benchmark.outputDir` names another directory. See the
[Java guide](JAVA.md) for the API and CLI reference.

### Maven Central deployment

The `release` profile of `java/pom.xml` signs the parent, library and CLI
artefacts and uploads them to the Central Portal as one deployment:

```bash
mvn -f java/pom.xml -P release clean deploy
```

Automatic publication is switched off. The command returns once Central has
validated the deployment, which then stays unpublished until a maintainer
publishes it by hand under Deployments at <https://central.sonatype.com>, or
drops it there. A release therefore reaches Maven Central only after that
manual confirmation.

## C++ source build

The C++ library is header-only. Release downloads provide the headers as
`smsd-cpp-8.0.0-headers.tar.gz`; in the source tree they are in `cpp/include`.

The library requires CMake 3.18 or later; the test commands below use
`ctest --test-dir`, which requires CMake 3.20 or later. On Windows, use
Developer PowerShell for Visual Studio with the Desktop development with C++
workload. From the repository root:

```text
cmake -S cpp -B build/cpu -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpu --config Debug --parallel 4
ctest --test-dir build/cpu --build-config Debug --output-on-failure
```

Batch operations use OpenMP when CMake finds it and otherwise run
sequentially. On macOS, install `libomp` (for example with Homebrew) and add
its prefix to `CMAKE_PREFIX_PATH` to enable OpenMP. `-DSMSD_WITH_RDKIT=ON`
builds the optional RDKit adapter, which requires the RDKit development files
and C++20. The [C++ module README](../cpp/README.md#install-and-use-with-cmake)
shows how to install the package and use it from CMake, and the
[C++ guide](CPP.md) is the API reference.

## Python

Install the release wheel for your platform from PyPI or the GitHub release:

```bash
python -m pip install "smsd==8.0.0"
python -c "import smsd; print(smsd.__version__); assert smsd.parse_smiles('CC').n == 2"
```

Release wheels are built for the standard CPython 3.14 interpreter (`cp314`),
not the free-threaded build (`cp314t`), and run batch operations on the CPU
with OpenMP. Other interpreters and architectures require a source build,
which supports Python 3.9 or later and needs a C++17 compiler, CMake and the
Python development headers. The root `pyproject.toml` builds the native extension and the
`python/smsd` package.

Linux or macOS, from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip build
.venv/bin/python -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Windows, from Developer PowerShell:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip build
.\.venv\Scripts\python.exe -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

The wheel is written to `dist/`; install it with the same environment's
Python. For development, install the checkout in editable mode with the test
dependencies, then run the tests:

```bash
python -m pip install -e ".[dev]"
python -c "import smsd; assert smsd.Substructure.find('CC').match('CCC').found"
python -m pytest python/tests -q
```

Metal (macOS) and CUDA are optional batch-screening backends. Source builds
enable each one when its toolkit is found; set `SMSD_BUILD_METAL` and
`SMSD_BUILD_CUDA` to `OFF`, as above, for a CPU-only build. RDKit is optional
and is needed only for RDKit molecule conversion and RDKit image depiction
(`python -m pip install "smsd[rdkit]==8.0.0"`). See the [Python guide](PYTHON.md) for
examples and build options.

## Docker

The repository `Dockerfile` builds an image that runs the Java CLI. From the
repository root:

```bash
docker build -t smsd .
docker run --rm smsd --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

The image holds `LICENSE`, `NOTICE` and the CDK licence files
(`licenses/cdk/`) in `/usr/share/doc/smsd`.

A manually dispatched release build can also attach the Java CLI image as a
compressed archive (`smsd-cli-image.tar.gz`), which `docker load` imports.

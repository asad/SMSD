# Installation and source builds

This guide covers installing SMSD Pro 8.0.1 and building it from source.
Release files are on the [GitHub releases page](https://github.com/asad/SMSD/releases);
the [README](../README.md) lists the Python and Maven packages.

## Platform requirements

| Package | Platform | Requirements |
|---|---|---|
| Java library, CLI JAR and portable CLI | Windows, macOS or Linux | Java 8 or later; Java 25 LTS preferred |
| Python wheel | Windows 10 or later, x86_64 (AMD64) | CPython 3.14 |
| Python wheel | Linux with glibc 2.28 or later, x86_64 (AMD64) | CPython 3.14 |
| Python wheel | macOS 26 or later, arm64 | CPython 3.14 |
| Python source build | Any platform with a C++17 compiler | Python 3.9 or later; CMake 3.18 or later |
| C++ headers | Any platform | C++17 compiler; neither RDKit nor CDK is needed |

Minimum operating-system versions are build targets.

## Verify downloads

Each GitHub release includes a `SHA256SUMS` file. Download it into the same
directory as the release files, then check each file you use.

macOS or Linux:

```bash
grep ' smsd-cli-8.0.1-jar-with-dependencies.jar$' SHA256SUMS | shasum -a 256 -c
```

Windows PowerShell (prints `True` when the checksum matches):

```powershell
$entry = Select-String -Path SHA256SUMS -SimpleMatch ' smsd-cli-8.0.1-jar-with-dependencies.jar'
(Get-FileHash .\smsd-cli-8.0.1-jar-with-dependencies.jar -Algorithm SHA256).Hash -eq $entry.Line.Split(' ')[0]
```

Replace the file name to check another download.

<a id="java-cli-installers"></a>

## Portable Java CLI

`smsd-cli-8.0.1-jar-with-dependencies.jar` contains the SMSD library, CDK 2.13
and the other CLI dependencies, and runs on Java 8 or later. From its directory:

```bash
java -jar smsd-cli-8.0.1-jar-with-dependencies.jar --version
java -jar smsd-cli-8.0.1-jar-with-dependencies.jar \
  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

The second command runs a substructure search and writes every mapping as JSON
to standard output. The same JAR is on Maven Central as
`com.bioinceptionlabs:smsd-cli:8.0.1` with the classifier `jar-with-dependencies`.

`smsd-8.0.1-cli.tar.gz` holds Unix and Windows launchers. Extract it into one
directory so that `bin` and `repo` stay together:

```bash
mkdir smsd-cli
tar -xzf smsd-8.0.1-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

On Windows, run `smsd-cli\bin\smsd.bat`. Neither package includes a Java runtime.

## Java source build

Use a JDK (Java 25 LTS preferred) and Maven 3.9 or later. From the repository
root:

```bash
mvn -f java/pom.xml clean package
java -jar java/cli/target/smsd-cli-8.0.1-jar-with-dependencies.jar --version
```

The build runs the default tests, writes the library JAR to `java/core/target`,
and writes the CLI JARs and the portable launchers to `java/cli/target`.
`mvn -f java/pom.xml install` also installs `com.bioinceptionlabs:smsd:8.0.1`
and `com.bioinceptionlabs:smsd-cli:8.0.1` in the local Maven repository; the
[Java module README](../java/README.md#installation) gives their Maven and
Gradle dependencies.

To run the slower test suites as well:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing clean verify
```

The code is compiled for Java 8 (`--release 8`). To run the tests on a Java 8
runtime, set `SMSD_JAVA8_HOME` to a Java 8 JDK:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" test
```

See the [Java guide](JAVA.md) for the API and CLI reference.

## C++ source build

The C++ library is header-only: use `cpp/include`, or the release archive
`smsd-cpp-8.0.1-headers.tar.gz`. The CMake build needs CMake 3.18 or later, and
`ctest --test-dir` needs CMake 3.20 or later. On Windows, use Developer
PowerShell for Visual Studio with the Desktop development with C++ workload.
From the repository root:

```text
cmake -S cpp -B build/cpu -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpu --config Debug --parallel 4
ctest --test-dir build/cpu --build-config Debug --output-on-failure
```

Batch operations use OpenMP when CMake finds it; on macOS, install `libomp`
(for example with Homebrew) and add its prefix to `CMAKE_PREFIX_PATH`.
`-DSMSD_WITH_RDKIT=ON` builds the optional RDKit adapter, which needs the RDKit
development files and C++20. The
[C++ module README](../cpp/README.md#install-and-use-with-cmake) shows how to
install the package and use it from CMake; the [C++ guide](CPP.md) is the API
reference.

## Python

Install the release wheel for your platform from PyPI or the GitHub release:

```bash
python -m pip install "smsd==8.0.1"
python -c "import smsd; print(smsd.__version__); assert smsd.parse_smiles('CC').n == 2"
```

Wheels are built for the standard CPython 3.14 interpreter (`cp314`), not the
free-threaded build (`cp314t`). Other interpreters and platforms need a source
build, with Python 3.9 or later, a C++17 compiler, CMake and the Python
development headers.

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

The wheel is written to `dist/`; install it with the same environment's Python.
For development, install the checkout in editable mode with the test
dependencies, then run the tests:

```bash
python -m pip install -e ".[dev]"
python -c "import smsd; assert smsd.Substructure.find('CC').match('CCC').found"
python -m pytest python/tests -q
```

Metal (macOS) and CUDA batch screening are optional: a source build enables
each one when its toolkit is found, and setting `SMSD_BUILD_METAL` and
`SMSD_BUILD_CUDA` to `OFF`, as above, gives a CPU-only build. RDKit is needed
only for RDKit molecule conversion and RDKit depiction
(`python -m pip install "smsd[rdkit]==8.0.1"`). See the
[Python guide](PYTHON.md) for examples and build options.

## Docker

The repository `Dockerfile` builds an image that runs the Java CLI:

```bash
docker build -t smsd .
docker run --rm smsd --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

The image holds `LICENSE`, `NOTICE` and the CDK licence files in
`/usr/share/doc/smsd`.

## Check the documented examples

After a source build, this runs every Python, Java and C++ code block in the
guides; examples that need RDKit are skipped when it is absent:

```bash
python scripts/check_examples.py \
  --jar java/cli/target/smsd-cli-8.0.1-jar-with-dependencies.jar \
  --report build/example-check.json
```

Add `--languages python` (or `java`, `cpp`) to check one language.

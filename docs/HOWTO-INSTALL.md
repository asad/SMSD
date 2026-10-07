# Build and run SMSD 7.2.2

Download packages from the [GitHub release](https://github.com/asad/SMSD/releases/tag/v7.2.2).
PyPI and Maven Central publication are pending. Verify downloads against
`SHA256SUMS` before installing.

| Distribution | Requirements |
|---|---|
| Java library and portable CLI | Java 8 or later on Windows, macOS or Linux; Java 25 LTS preferred |
| Java CLI installers | Windows x86_64 MSI, macOS arm64 DMG, Linux x86_64 DEB; Java 25.0.4.1+1 included |
| Python wheels | CPython 3.14: Windows 10+/x86_64, Linux x86_64/glibc 2.28+, macOS arm64/macOS 26+ |
| C++ headers | C++17 compiler; no RDKit or CDK dependency for core matching |

See [installer instructions](INSTALLERS.md) and the
[validation report](VALIDATION_7.2.2.md) for supported platforms and checks.

## Portable Java CLI

Run these commands in the directory containing the downloaded JAR:

```bash
java -version
java -jar smsd-7.2.2-jar-with-dependencies.jar --version
java -jar smsd-7.2.2-jar-with-dependencies.jar \
  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

Windows PowerShell:

```powershell
java -jar .\smsd-7.2.2-jar-with-dependencies.jar --version
java -jar .\smsd-7.2.2-jar-with-dependencies.jar --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

The portable archive also contains Unix and Windows launchers. Extract it into
one directory so `bin` and `repo` remain together:

```bash
mkdir smsd-cli
tar -xzf smsd-7.2.2-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

On Windows, use `smsd-cli\bin\smsd.bat`. The launchers use Java from `PATH`;
the portable packages do not include a runtime.

## Java source build

Use Maven 3.9+ and a JDK; Java 25 LTS is preferred. Run from the repository root:

```bash
mvn -f java/pom.xml -U clean package
java -jar java/target/smsd-7.2.2-jar-with-dependencies.jar --version
```

The shaded JAR includes CDK 2.13. The root `pom.xml` also builds the Java module.
To include the normally excluded algorithm and stress tests:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing clean verify
```

To test Java 8 compatibility, set `SMSD_JAVA8_HOME` to an installed Java 8 JDK:

```bash
mvn -f java/pom.xml -Dslow.tests.exclude=nothing \
  "-Djvm=$SMSD_JAVA8_HOME/bin/java" test
```

Opt-in corpus benchmarks are separate from correctness tests:

```bash
mvn -f java/pom.xml test -Dslow.tests.exclude=nothing \
  '-Dtest=BenchmarkSuiteTest*,ExternalBenchmarkTest*,JavaCdkVsSmsdBenchmarkTest' \
  -Dbenchmark=true -Dsmsd.benchmark=true \
  -Dsmsd.benchmark.timeoutMs=1000 -Dsmsd.benchmark.rounds=1 \
  -Dsmsd.benchmark.warmup=0 \
  -Dsmsd.benchmark.outputDir="$PWD/build/local-benchmarks/java"
```

## C++ source build

Use a C++17 compiler and CMake 3.20+ for these commands. The library itself
supports CMake 3.18+. On Windows, use Developer PowerShell with Visual Studio's
Desktop development with C++ tools installed.

```text
cmake -S cpp -B build/cpu -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpu --config Debug --parallel 4
ctest --test-dir build/cpu --build-config Debug --output-on-failure
```

OpenMP is used when available; otherwise batch operations run sequentially.
The optional C++ RDKit adapter requires RDKit development files and C++20.
See the [C++ guide](CPP.md) for integration examples.

## Python installation and source build

Install the wheel matching your interpreter and platform from the GitHub
release. When 7.2.2 is listed on PyPI:

```bash
python -m pip install smsd==7.2.2
python -c 'import smsd; print(smsd.__version__); assert smsd.is_substructure("CC", "CCC")'
```

Wheels target ordinary CPython 3.14, rather than the free-threaded `cp314t`
interpreter. Other interpreters and architectures require a source build.
Source metadata supports Python 3.9+; use a C++17 compiler, CMake and Python
development headers. The root `pyproject.toml` builds the C++ extension and
`python/smsd` package.

Linux or macOS, from the repository root:

```bash
python3.14 -m venv .venv
.venv/bin/python -m pip install --upgrade pip build
.venv/bin/python -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Windows Developer PowerShell:

```powershell
py -3.14 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip build
.\.venv\Scripts\python.exe -m build --wheel -Ccmake.define.SMSD_BUILD_METAL=OFF -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Install the resulting wheel with that environment's Python. For local tests,
install `pytest` and `pytest-timeout`, then run:

```bash
python -m pytest python/tests -q --import-mode=importlib
```

RDKit is optional for molecule conversion and drawing. Metal and CUDA are
optional source-build backends; core batch matching uses CPU/OpenMP.
See the [Python guide](PYTHON.md) for examples and build options.

## Docker

The GitHub release includes Java CLI image archives for Linux x86_64 and arm64,
with load/run commands in its `DOCKER.md` asset. To build locally:

```bash
docker build -t smsd .
docker run --rm smsd --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
```

For release preparation, platform collection and package publication, see
[publishing](PUBLISHING.md).

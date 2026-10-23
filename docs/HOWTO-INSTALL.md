# Installation and source builds

This guide covers installing SMSD Pro 8.0.2 and building it from source.
Release files are on the [GitHub releases page](https://github.com/asad/SMSD/releases);
the [README](../README.md) lists the Python and Maven packages.

<a id="platform-requirements"></a>

## Supported platforms

| Package | Platform | Interpreter or runtime | Compute |
|---|---|---|---|
| `smsd` wheel | Windows 10 or later, x86_64 (AMD64) | CPython 3.10–3.14 | CPU/OpenMP |
| `smsd` wheel | Linux with glibc 2.28 or later, x86_64 or aarch64 | CPython 3.10–3.14 | CPU/OpenMP |
| `smsd` wheel | macOS 26 or later, arm64 (Apple silicon) | CPython 3.10–3.14 | CPU/OpenMP |
| `smsd-cuda` wheel | Linux with glibc 2.34 or later, x86_64 | CPython 3.10–3.14; matching `smsd`; NVIDIA driver compatible with the installed CUDA 12 runtime | Automatic ECFP4 screening on Tesla T4 (compute capability 7.5) |
| CPU Docker image | Linux amd64 or arm64 | Python 3.14.9 and Java 25 included | CPU/OpenMP |
| GPU Docker image | Linux amd64; Windows through WSL2 with NVIDIA GPU support | Python 3.14.9 and CUDA 12.8.1 runtime included; compatible host driver and NVIDIA Container Toolkit | Optional NVIDIA CUDA screening |
| Java library, CLI JAR and portable CLI | Windows, macOS or Linux | Java 8 or later; Java 25 LTS preferred | CPU |
| Java CLI installers | Windows x86_64 (MSI), macOS arm64 (DMG), Linux amd64 (DEB) | Bundled Java 25 runtime | CPU |
| Python source build | A platform with a C++17 compiler | Python 3.10 or later; CMake 3.18 or later | CPU; optional CUDA package built separately |
| C++ headers | A platform with a C++17 compiler | No Python or Java runtime needed | CPU; native CUDA source build optional |

Minimum operating-system versions are build targets. The standard wheel matrix
has five CPython versions on four platforms. Free-threaded CPython and Intel
macOS need source builds. macOS GPU support and native Windows CUDA wheels are
not included. Search matching runs on the CPU in every package.

## Verify downloads

Each GitHub release includes a `SHA256SUMS` file. Download it into the same
directory as the release files, then check each file you use.

macOS or Linux:

```bash
grep ' smsd-cli-8.0.2-jar-with-dependencies.jar$' SHA256SUMS | shasum -a 256 -c
```

Windows PowerShell (prints `True` when the checksum matches):

```powershell
$entry = Select-String -Path SHA256SUMS -SimpleMatch ' smsd-cli-8.0.2-jar-with-dependencies.jar'
(Get-FileHash .\smsd-cli-8.0.2-jar-with-dependencies.jar -Algorithm SHA256).Hash -eq $entry.Line.Split(' ')[0]
```

Replace the file name to check another download.

<a id="java-cli-installers"></a>

## Portable Java CLI

`smsd-cli-8.0.2-jar-with-dependencies.jar` contains the SMSD library, CDK 2.13
and the other CLI dependencies, and runs on Java 8 or later. From its directory:

```bash
java -jar smsd-cli-8.0.2-jar-with-dependencies.jar --version
java -jar smsd-cli-8.0.2-jar-with-dependencies.jar \
  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

The second command runs a substructure search and writes every mapping as JSON
to standard output. The same JAR is on Maven Central as
`com.bioinceptionlabs:smsd-cli:8.0.2` with the classifier `jar-with-dependencies`.

`smsd-8.0.2-cli.tar.gz` holds Unix and Windows launchers. Extract it into one
directory so that `bin` and `repo` stay together:

```bash
mkdir smsd-cli
tar -xzf smsd-8.0.2-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

On Windows, run `smsd-cli\bin\smsd.bat`. Neither package includes a Java runtime.

The MSI, DMG and DEB installers include Java 25. Their architectures are listed
in the [supported-platforms table](#supported-platforms).

## Java source build

Use a JDK (Java 25 LTS preferred) and Maven 3.9 or later. From the repository
root:

```bash
mvn -f java/pom.xml clean package
java -jar java/cli/target/smsd-cli-8.0.2-jar-with-dependencies.jar --version
```

The build runs the default tests, writes the library JAR to `java/core/target`,
and writes the CLI JARs and the portable launchers to `java/cli/target`.
`mvn -f java/pom.xml install` also installs `com.bioinceptionlabs:smsd:8.0.2`
and `com.bioinceptionlabs:smsd-cli:8.0.2` in the local Maven repository; the
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
`smsd-cpp-8.0.2-headers.tar.gz`. The CMake build needs CMake 3.18 or later, and
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
python -m pip install "smsd==8.0.2"
python -c "import smsd; print(smsd.__version__); assert smsd.parse_smiles('CC').n == 2"
```

Wheels are built for standard CPython 3.10–3.14 (`cp310` through `cp314`),
not the free-threaded builds. Other interpreters and platforms need a source
build, with Python 3.10 or later, a C++17 compiler, CMake and the Python
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

The standard wheels use CPU/OpenMP. For optional NVIDIA CUDA screening on a
supported Linux system, install the matching package:

```bash
python -m pip install "smsd==8.0.2" "smsd-cuda==8.0.2"
```

The CUDA runtime is installed as a package dependency; the host still needs a
compatible NVIDIA driver. See [CUDA screening](PYTHON.md#optional-cuda-screening)
for dispatch and CPU fallback behaviour. RDKit is needed only for RDKit
molecule conversion and RDKit depiction
(`python -m pip install "smsd[rdkit]==8.0.2"`). See the
[Python guide](PYTHON.md) for examples and build options.

## Docker

The CPU image includes the Java CLI and the Python package. The default command
is the Java CLI; use `--entrypoint python` for Python. After loading the image
archive from the release:

```bash
docker run --rm smsd:8.0.2-cpu-amd64 \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
docker run --rm --entrypoint python smsd:8.0.2-cpu-amd64 \
  -c "import smsd; assert smsd.Substructure.find('CC').match('CCC').found"
```

On an arm64 host, use `smsd:8.0.2-cpu-arm64`. The GPU image includes the Python
package, matching `smsd-cuda` and CUDA runtime, with Python as its entry point:

```bash
docker run --rm --gpus all smsd:8.0.2-cuda12.8 \
  -c "import smsd; print(smsd.__version__); print(smsd.gpu_device_info())"
```

It needs the NVIDIA Container Toolkit and a compatible host driver, including
when used through WSL2 on Windows. The final image contains no CUDA toolkit.
Backend availability does not prove that a screening call used the GPU.

Source builds use a prepared context containing the matching release wheels
and their checksums. Given a completed release directory with its source and
wheel manifests, prepare a fresh context and use the `source_commit` recorded
in `SOURCE_MANIFEST.json` as the revision:

```bash
python scripts/prepare_docker_context.py \
  --release-dir release --output docker-context --architecture amd64
docker build --platform linux/amd64 \
  --build-arg SMSD_REVISION=SOURCE_COMMIT \
  -f docker-context/Dockerfile -t smsd:8.0.2-cpu-amd64 docker-context
```

Replace `SOURCE_COMMIT` with that recorded value. For arm64, use a separate
context with `--architecture arm64`, platform `linux/arm64` and tag
`smsd:8.0.2-cpu-arm64`. To prepare the GPU context, add `--cuda-wheel`,
`--cuda-provenance` and `--runtime-wheel` with the matching validated add-on,
its provenance receipt and the pinned runtime wheel. Build that context with
`Dockerfile.cuda`, platform `linux/amd64` and tag `smsd:8.0.2-cuda12.8`.
`python scripts/prepare_docker_context.py --help` lists the required inputs.
The CPU image holds `LICENSE`, `NOTICE` and the CLI dependency
licence files in `/usr/share/doc/smsd`; both images retain the notices of their
bundled runtimes.

## Check the documented examples

After a source build, this runs every Python, Java and C++ code block in the
guides; examples that need RDKit are skipped when it is absent:

```bash
python scripts/check_examples.py \
  --jar java/cli/target/smsd-cli-8.0.2-jar-with-dependencies.jar \
  --report build/example-check.json
```

Add `--languages python` (or `java`, `cpp`) to check one language.

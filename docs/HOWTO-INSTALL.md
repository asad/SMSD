# Installation and source builds

This guide covers installing SMSD Pro 8.1.0 and building it from source.
Release files are on the [GitHub releases page](https://github.com/asad/SMSD/releases);
the [README](../README.md) lists the Python and Maven packages and the conda
recipes.

## Install SMSD Pro 8.1.0

Install the Python package from PyPI:

```bash
python -m pip install "smsd==8.1.0"
```

For NVIDIA GPU screening on Linux x86_64, add the optional CUDA package of the
same version: `python -m pip install "smsd==8.1.0" "smsd-cuda==8.1.0"`. The
Java library and command-line tool stay at 8.0.3 on Maven Central, as
`com.bioinceptionlabs:smsd:8.0.3` and `com.bioinceptionlabs:smsd-cli:8.0.3`;
the [Java module README](../java/README.md#installation) gives their Maven and
Gradle dependencies. The [CPU Docker images](#docker) carry the 8.1.0
command-line tool, and a [Java source build](#java-source-build) produces the
8.1.0 Java library and tool. The C++ headers need a C++17 compiler and no
Python or Java runtime.

See [Python adapters and MCP](PYTHON_ADAPTERS.md),
[packed corpora and clustering](PACKED_CORPORA.md), and
[SMARTS export](SMARTS_EXPORT.md). PyTorch on macOS may require an OpenMP-free
SMSD source build to avoid incompatible bundled runtimes. Portable SMARTS
rejects unsupported matching constraints instead of silently dropping them.

The [conda guide](../packaging/conda/README.md) describes the recipes that
build `smsd` and, on Linux x86_64, `smsd-cuda` for NVIDIA GPU screening. It
also covers the conda, mamba and Pixi install routes.

<a id="platform-requirements"></a>

## Supported platforms

| Package | Platform | Interpreter or runtime | Compute |
|---|---|---|---|
| `smsd` wheel | Windows 10 or later, x86_64 (AMD64) | CPython 3.12–3.15 | CPU/OpenMP |
| `smsd` wheel | Linux with glibc 2.28 or later, x86_64 or aarch64 | CPython 3.12–3.15 | CPU/OpenMP |
| `smsd` wheel | macOS 14 or later, arm64 (Apple silicon) | CPython 3.12–3.15 | CPU/OpenMP |
| `smsd-cuda` wheel | Linux with glibc 2.34 or later, x86_64 | CPython 3.12–3.15; matching `smsd`; NVIDIA driver compatible with the installed CUDA 12 runtime | Automatic ECFP4 screening on Tesla T4 (compute capability 7.5) |
| Java library, CLI JAR and portable CLI | Windows, macOS or Linux | Java 8 or later; Java 25 LTS preferred | CPU |
| Java CLI installers (8.0.3 release) | Windows x86_64 (MSI), macOS arm64 (DMG), Linux amd64 (DEB) | Bundled Java 25 runtime | CPU |
| Python source build | A platform with a C++17 compiler | Python 3.12 or later; CMake 3.18 or later | CPU; optional CUDA package built separately |
| C++ headers | A platform with a C++17 compiler | No Python or Java runtime needed | CPU; native CUDA or macOS Metal source build optional |

Minimum operating-system versions are build targets. The standard wheel matrix
has four CPython versions on four platforms. Free-threaded CPython and Intel
macOS need source builds. Python wheels do not include macOS GPU support or
native Windows CUDA. Native C++ source builds on macOS include [optional Metal
support](CPP.md#optional-metal-fingerprint-screening), which a program uses only
when it links `smsd_metal`; Metal fingerprint screening runs only when requested
explicitly. Search matching runs on the CPU in every package. By default,
automatic ECFP4 screening of a prepared `TargetCorpus` uses the GPU only on a
Tesla T4 and only for a corpus of at least 100,000 targets; the
[CUDA screening](PYTHON.md#optional-cuda-screening) section gives the explicit
routes for other devices.

## Verify downloads

Each GitHub release includes a `SHA256SUMS` file. Download it into the same
directory as the release files, then check each file you use.

macOS or Linux:

```bash
grep ' smsd-cpp-8.1.0-headers.tar.gz$' SHA256SUMS | shasum -a 256 -c
```

Windows PowerShell (prints `True` when the checksum matches):

```powershell
$entry = Select-String -Path SHA256SUMS -SimpleMatch ' smsd-cpp-8.1.0-headers.tar.gz'
(Get-FileHash .\smsd-cpp-8.1.0-headers.tar.gz -Algorithm SHA256).Hash -eq $entry.Line.Split(' ')[0]
```

Replace the file name to check another download.

<a id="java-cli-installers"></a>

## Portable Java CLI

The standalone Java command-line tool stays at 8.0.3 in this release, so this
section uses the files of the [8.0.3 release](https://github.com/asad/SMSD/releases/tag/v8.0.3);
the [CPU Docker images](#docker) carry the 8.1.0 tool, and a
[Java source build](#java-source-build) produces it.

`smsd-cli-8.0.3-jar-with-dependencies.jar` contains the SMSD library, CDK 2.13
and the other CLI dependencies, and runs on Java 8 or later. From its directory:

```bash
java -jar smsd-cli-8.0.3-jar-with-dependencies.jar --version
java -jar smsd-cli-8.0.3-jar-with-dependencies.jar \
  --Q SMI --q "CCN" --T SMI --t "CCCNC" -m --json -
```

The second command runs a substructure search and writes every mapping as JSON
to standard output. The same JAR is on Maven Central as
`com.bioinceptionlabs:smsd-cli:8.0.3` with the classifier `jar-with-dependencies`.

`smsd-8.0.3-cli.tar.gz` holds Unix and Windows launchers. Extract it into one
directory so that `bin` and `repo` stay together:

```bash
mkdir smsd-cli
tar -xzf smsd-8.0.3-cli.tar.gz -C smsd-cli
./smsd-cli/bin/smsd --version
```

On Windows, run `smsd-cli\bin\smsd.bat`. Neither package includes a Java runtime.

The 8.0.3 MSI, DMG and DEB installers include Java 25. Their architectures are listed
in the [supported-platforms table](#supported-platforms).

## Java source build

Use a JDK (Java 25 LTS preferred) and Maven 3.9 or later. From the repository
root:

```bash
mvn -f java/pom.xml clean package
java -jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar --version
```

The build runs the default tests, writes the library JAR to `java/core/target`,
and writes the CLI JARs and the portable launchers to `java/cli/target`.
`mvn -f java/pom.xml install` also installs `com.bioinceptionlabs:smsd:8.1.0`
and `com.bioinceptionlabs:smsd-cli:8.1.0` in the local Maven repository. To
use them, set version 8.1.0 in the Maven or Gradle dependencies of the
[Java module README](../java/README.md#installation); Gradle also needs
`mavenLocal()` in its repositories.

To run the slow tests as well:

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
`smsd-cpp-8.1.0-headers.tar.gz`. The CMake build needs CMake 3.18 or later, and
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
python -m pip install "smsd==8.1.0"
python -c "import smsd; print(smsd.__version__); assert smsd.parse_smiles('CC').n == 2"
```

Wheels are built for standard CPython 3.12–3.15 (`cp312` through `cp315`),
not the free-threaded builds. Other interpreters and platforms need a source
build, with Python 3.12 or later, a C++17 compiler, CMake and the Python
development headers.

Linux or macOS, from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip build
.venv/bin/python -m build --wheel
```

Windows, from Developer PowerShell:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip build
.\.venv\Scripts\python.exe -m build --wheel
```

Source builds leave Metal and CUDA out by default. Python screening and search
do not run on them, so `-Ccmake.define.SMSD_BUILD_METAL=ON` or
`-Ccmake.define.SMSD_BUILD_CUDA=ON` (or `AUTO`, which includes each one when it
is available) changes only what `gpu_is_available()` and `gpu_device_info()`
report. NVIDIA screening comes from the separate `smsd-cuda` package.
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
python -m pip install "smsd==8.1.0" "smsd-cuda==8.1.0"
```

The CUDA runtime is installed as a package dependency; the host still needs a
compatible NVIDIA driver. See [CUDA screening](PYTHON.md#optional-cuda-screening)
for dispatch and CPU fallback behaviour. RDKit is needed only for RDKit
molecule conversion and RDKit depiction
(`python -m pip install "smsd[rdkit]==8.1.0"`). See the
[Python guide](PYTHON.md) for examples and build options.

## Docker

The CPU images include the 8.1.0 Java CLI and the Python package. The default
command is the Java CLI; use `--entrypoint python` for Python. After loading
the image archive from the release:

```bash
docker load --input smsd-8.1.0-cpu-linux-amd64.tar.xz
docker run --rm smsd:8.1.0-cpu-amd64 \
  --Q SMI --q "c1ccccc1" --T SMI --t "c1ccc(O)cc1" --json -
docker run --rm --entrypoint python smsd:8.1.0-cpu-amd64 \
  -c "import smsd; assert smsd.Substructure.find('CC').match('CCC').found"
```

On an arm64 host, load `smsd-8.1.0-cpu-linux-arm64.tar.xz` and use
`smsd:8.1.0-cpu-arm64`. The GPU image includes the Python
package, matching `smsd-cuda` and CUDA runtime, with Python as its entry point:

```bash
docker load --input smsd-8.1.0-cuda-linux-amd64.tar.xz
docker run --rm --gpus all smsd:8.1.0-cuda12.8 \
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
docker buildx build --load --platform linux/amd64 \
  --build-arg SMSD_REVISION=SOURCE_COMMIT \
  -f docker-context/Dockerfile -t smsd:8.1.0-cpu-amd64 docker-context
```

Replace `SOURCE_COMMIT` with that recorded value. For arm64, use a separate
context with `--architecture arm64`, platform `linux/arm64` and tag
`smsd:8.1.0-cpu-arm64`. To prepare the GPU context, add `--cuda-wheel`,
`--cuda-provenance` and `--runtime-wheel` with the matching `smsd-cuda` wheel,
its provenance record and the pinned runtime wheel. Build that context with
`Dockerfile.cuda`, platform `linux/amd64` and tag `smsd:8.1.0-cuda12.8`.
The GPU Dockerfile needs BuildKit with its bundled Dockerfile frontend 1.2
or later. Its wheel inputs are mounted for installation rather than copied
into final image layers; the same installed-file and runtime checks still run.
`python scripts/prepare_docker_context.py --help` lists the required inputs.
The CPU images hold `LICENSE`, `NOTICE` and the CLI dependency
licence files in `/usr/share/doc/smsd`; every image retains the notices of its
bundled runtimes.

To save a built image as a release archive:

```bash
python scripts/archive_docker_image.py --image smsd:8.1.0-cuda12.8 \
  --revision SOURCE_COMMIT --architecture amd64 \
  --output smsd-8.1.0-cuda-linux-amd64.tar.xz
docker load --input smsd-8.1.0-cuda-linux-amd64.tar.xz
```

If the Docker engine has no xz decompressor, decompress on the host before
loading:

```bash
xz --decompress --stdout smsd-8.1.0-cuda-linux-amd64.tar.xz | docker load
```

For a CPU image, use `--image smsd:8.1.0-cpu-amd64` and output
`smsd-8.1.0-cpu-linux-amd64.tar.xz`; replace both architecture values with
`arm64` for that platform. The helper saves the documented tag and checks that
its image identity, revision label and platform stay unchanged during export.
It requires the xz archive to be smaller than 2,000,000,000 bytes. An existing
output is never overwritten. Failed exports retain the partial archive,
stderr log and a JSON record with `accepted=false`.

## Check the documented examples

After a source build, this runs every Python, Java and C++ code block in the
main guides; examples that need RDKit are skipped when it is absent:

```bash
python scripts/check_examples.py \
  --jar java/cli/target/smsd-cli-8.1.0-jar-with-dependencies.jar \
  --report build/example-check.json
```

Add `--languages python` (or `java`, `cpp`) to check one language.

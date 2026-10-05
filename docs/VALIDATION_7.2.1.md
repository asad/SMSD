# SMSD 7.2.1 release validation

Version 7.2.1 is in preparation; publication is pending. This record tracks
fresh checks after separating the Java, C++ and Python source modules. It
does not inherit a passing result merely because a 7.2.0 artifact passed.

## Source and layout

Java sources and launchers are in `java/src/`, its Maven module is
`java/pom.xml`, and generated Java artifacts are in `java/target/`. The root
Maven aggregator supports `mvn verify`. C++ sources remain in `cpp/`, Python
sources/tests remain in `python/`, and shared scripts, documentation and
licenses stay at the root. The root `pyproject.toml` is the single Python
manifest because it builds the extension from the C++ tree.

The release commit, production-file fingerprint, source archive hash and
final asset checksums will be recorded after the source is frozen and the
checks below complete. Every wheel must correspond to that same source.

## Release gates

| Gate | 7.2.1 status | Required evidence |
|---|---|---|
| Source freeze | Pending | Exact commit, clean checkout and coherent Java/CMake/Python versions |
| Java 25 module | Pending | Clean full verification; library, CLI, sources and Javadoc; version/license checks and portable launchers |
| macOS arm64 CPU/OpenMP | Pending | Local native Debug suites and installed CPython 3.14 wheel tests; repaired dependencies and deployment tag |
| Linux x86_64 CPU/OpenMP | Pending | Local manylinux build with glibc 2.28 target; native Debug suites and installed CPython 3.14 tests |
| Windows x86_64 CPU/OpenMP | Pending | Native GitHub Windows build from the exact 7.2.1 source; MSVC Debug suites, repaired DLLs and installed CPython 3.14 tests |
| Three-wheel collection | Pending | Source/version agreement, binary architecture, runtime libraries, wrappers/headers, licenses and all RECORD hashes |
| Source and release assets | Pending | Complete source inputs, no private files, strict metadata checks and verified checksums |
| PyPI, Maven Central and GitHub | Pending | Publication followed by clean download/install checks |

The CPU configuration registers 12 native suites with assertions enabled.
Record the actual Python pass/skip counts, interpreter and RDKit versions for
each installed wheel after execution. OpenMP must be active in the release
wheel; a sequential fallback does not prove its runtime packaging.

The compact wheel set is ordinary CPython 3.14: Linux x86_64, macOS arm64
(macOS 26+ deployment target), and Windows x86_64. It excludes free-threaded
CPython, Intel macOS and Linux arm64 wheels. Other architectures can use
source builds, without a claim that this release tested them. CUDA and the
optional C++ RDKit adapter are outside the required CPU wheel gates; any new
GPU or adapter checks must be listed separately.

## Reproduce and assemble

From the repository root, verify the Java module and CPU native tests:

```bash
mvn -f java/pom.xml -B -Dslow.tests.exclude=nothing clean verify \
  org.apache.maven.plugins:maven-source-plugin:3.3.1:jar-no-fork \
  org.apache.maven.plugins:maven-javadoc-plugin:3.6.3:jar
java/src/scripts/smsd --version
java/target/appassembler/bin/smsd --version

cmake -S cpp -B build/release-preflight -DCMAKE_BUILD_TYPE=Debug \
  -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_CUDA=OFF -DSMSD_BUILD_METAL=OFF
cmake --build build/release-preflight --config Debug --parallel 4
ctest --test-dir build/release-preflight --build-config Debug --output-on-failure
```

Follow [publishing preparation](PUBLISHING.md) to build macOS and Linux
locally, validate Windows with the manual GitHub workflow, and retain the
source manifests and installed-wheel logs. The corrected 7.2.0 Windows run
([37293365203](https://github.com/asad/SMSD/actions/runs/37293365203)) passed and
is historical evidence, not a 7.2.1 test result.

Collect the verified wheels and check the complete set before publishing:

```bash
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 \
  --wheel-dir build/platform-release/linux \
  --wheel-dir build/platform-release/windows
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 --check-only
python -m twine check --strict \
  dist/release-7.2.1/smsd-7.2.1.tar.gz \
  dist/release-7.2.1/smsd-7.2.1-*.whl
```

Artifact inspection complements execution on the target operating system;
it cannot prove Windows runtime behavior from a macOS inspection alone.

## Historical search and benchmark evidence

The [7.2.0 validation](VALIDATION_7.2.0.md) and
[benchmark report](../benchmarks/RESULTS_7.2.0.md) retain their original source
versions, test counts, measured numbers, fingerprints and input hashes.
The retained raw archive is `smsd-7.2.0-benchmark-data.tar.gz`, SHA-256
`217e42f3b7f9cf7a2996120d57037826419ec075b23e82059afc2d44fb5a8ea5`.

The 7.2.1 version, layout and packaging changes have no new benchmark
measurements. These historical bounded-search and small-oracle results do
not prove global optimality for arbitrary molecular graphs, nor do macOS
timings establish performance on Linux or Windows.

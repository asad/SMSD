# SMSD 7.2.1 release validation

Version 7.2.1 is in preparation; publication is pending. This record tracks
fresh checks after separating the Java, C++ and Python source modules. It
does not inherit a passing result merely because a 7.2.0 artifact passed.

As checked on 2026-10-06, GitHub has published
[7.2.0](https://github.com/asad/SMSD/releases/tag/v7.2.0); PyPI and Maven Central
remain at 7.1.1. The native Windows 7.2.1 workflow passed, and all three
wheels have passed strict collection against the same source archive.

## Source and layout

Java sources and launchers are in `java/src/`, its Maven module is
`java/pom.xml`, and generated Java artifacts are in `java/target/`. The root
Maven aggregator supports `mvn verify`. C++ sources remain in `cpp/`, Python
sources/tests remain in `python/`, and shared scripts, documentation and
licenses stay at the root. The root `pyproject.toml` is the single Python
manifest because it builds the extension from the C++ tree.

The frozen build source is commit
`815efdfc9d7ecaf2f7957bff8f6e32882a503e44`, prepared from a clean checkout.
Its 215-file source distribution has SHA-256
`1ef98789ffeaf8f2fe22586a12eef46bdb2ec8af6226f98b12f4cc6f849e7a2c`.
Later documentation updates record results without changing the build inputs.
Every wheel must correspond to this frozen source; final asset checksums are
recorded in the release directory's `SHA256SUMS`.

## Release gates

| Gate | 7.2.1 status | Required evidence |
|---|---|---|
| Source freeze | Passed | Exact commit above, 215 source files checked, coherent Java/CMake/Python versions |
| Java 25 module | Passed across full and focused runs | 1,242 distinct passing cases, 15 opt-in skips; four JARs and both Unix launchers checked |
| macOS arm64 CPU/OpenMP | Passed from frozen archive | 12 native Debug suites, 691 Python passes and 8 skips on Python 3.14.8/RDKit 2026.03.6 |
| Docker CLI | Passed locally | Filtered source context; Linux arm64/Temurin 25.0.4.1; version/help and validated substructure/MCS results |
| Linux x86_64 CPU/OpenMP | Passed locally under emulation | Same frozen archive; 12 native Debug suites, 691 Python passes and 8 skips on Python 3.14.5/RDKit 2026.03.6 |
| Windows cross-compilation | Passed locally; no execution | All 12 native targets compile/link with MinGW GCC 16.2.0 and have AMD64 PE headers; sequential fallback |
| Windows x86_64 CPU/OpenMP | Passed on Windows Server 2022 | Same frozen archive; 12 MSVC Debug suites, 691 Python passes and 8 skips on Python 3.14.7/RDKit 2026.03.6; repaired DLLs checked |
| Three-wheel collection | Passed | Source/version agreement, binary architecture, runtime libraries, wrappers/headers, licenses and all RECORD hashes |
| Source and release assets | Passed locally | Complete source inputs, no private files, strict metadata checks, CLI/header packages and verified checksums |
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

## Recorded preparation checks

The first Java full run had 1,241 passes, 15 opt-in skips and one outer-guard
timeout. The affected fixture allowed a 30-second search and a 30-second test
guard, leaving no time for parsing or validation. After allowing 35 seconds
for 30-second drug searches and 12 seconds for 10-second pharmacophore
searches, all 30 affected fixtures passed. Search budgets and assertions
remain unchanged. Combined evidence covers 1,242 distinct passing cases;
it is not described as one clean full run. Root reactor release packaging
also passed with signing disabled. All four JARs contain the exact root
LICENSE/NOTICE; own class files target Java 25, source copies match, and both
Unix launchers report 7.2.1.

The final frozen-source macOS wheel passed 12 native Debug suites and 691 Python
tests with 8 optional skips. It used Python 3.14.8, RDKit 2026.03.6 and bundled
libomp 23.1.0. It executed on macOS 27.0.1, with a macOS 26 deployment tag;
execution on the minimum OS was not tested. Native CTest took 93.47 seconds;
installed-wheel pytest took 7.69 seconds. Delocate repair bundled libomp;
strict Twine metadata, architecture, all RECORD hashes, exact wrappers/headers
and legal copies passed. The wheel SHA-256 is
`535a08cb8baaebc2c9750fa5a915c0f40143f86a8b9f86e5716655cce5e2e903`.

The first complete Linux native run passed 11 suites and failed the existing
MCS deadline regression: a 5 ms budget exceeded its unchanged 100 ms guard.
Repeated probes on the same emulated x86_64 host measured about
146 ms for the 32-atom fixture. Individual seed extensions took at most
1.33 ms; throttled checks allowed many extensions after expiry. Immediate
checks at candidate boundaries reduced the focused fixture to about 5 ms,
with valid mappings and exact deadline restoration. The source also skips
already-expired seed and orientation setup. Fresh full native and installed
wheel runs now pass on macOS, Linux and Windows; the timeout assertion remains unchanged.
These timings describe that regression on an emulated host, not a general
performance comparison.

The final Linux wheel used the same frozen source archive in the local
manylinux 2.28 x86_64 container, running AlmaLinux 8.10/glibc 2.28 under QEMU
on an arm64 host. GCC 14.2.1 compiled the CPU/OpenMP build. All 12 Debug
suites passed in 402.06 seconds, including the MCS regression suite;
installed-wheel pytest passed 691 tests with 8 skips in 32.89 seconds.
The runtime was Python 3.14.5, RDKit 2026.03.6 and pytest 9.1.1, with active
OpenMP. Auditwheel repair and strict Twine metadata checks passed; collection
verified AMD64 ELF architecture, all RECORD hashes, exact source headers and
wrappers, and license copies. The bundled libgomp 8.5.0 runtime comes from the
AlmaLinux image; its executable section matches the original library, and the
extension's relative RPATH resolves the repaired copy. Auditwheel added both manylinux 2.27 and 2.28
tags; execution was tested on glibc 2.28 only. The repaired wheel SHA-256 is
`e2197aef8f2f310c77bf4637b32846bad2e58aef14a73efdc56c0558c8347a3d`.
Emulated execution establishes these checks, not native hardware timings.

Local MinGW cross-compilation from the frozen archive built all 12 native
Debug executables with assertions enabled and verified AMD64 PE headers.
That toolchain has no OpenMP runtime and used the sequential fallback. This
proves compilation/linking only; the native MSVC build below provides separate
Windows execution evidence.

The native Windows build passed in
[run 37394450131](https://github.com/asad/SMSD/actions/runs/37394450131), using
the byte-identical frozen source archive above. Windows Server 2022/AMD64 ran
all 12 MSVC Debug suites in 260.51 seconds and the installed-wheel Python suite
with 691 passes and 8 optional skips in 8.40 seconds. It used CPython 3.14.7,
RDKit 2026.03.6, pytest 9.1.1 and MSVC 19.44.35229.0; active OpenMP was checked.
The repaired wheel SHA-256 is
`36080c0929f8962a27bea33088d98dccd10cc17859f07d98ee21ad5df27e70fb`.

Both bundled Microsoft runtimes (`msvcp140.dll` and `vcomp140.dll`) are AMD64
release DLLs at version 14.44.35211.0, matching the extension's 14.44 linker
family. Their hashes match the explicitly selected Visual Studio
redistributables. Independent PE inspection resolved all 78 imports from those
DLLs and found no delay imports. Exact Microsoft license copies and the
supported delvewheel 1.13.1 loader passed inspection. Native execution also
checked UTF-8 file paths and the installed OpenMP backend. The workflow used
the updated checkout/upload/download actions and skipped publication.

Strict collection of the macOS, Linux and Windows wheels passed against the
single frozen archive. All wheel RECORD entries, installed Python wrappers,
27 C++ headers and legal copies were checked. Strict Twine metadata checks and
the complete release asset checksum list passed. These checks establish the
prepared artifacts; clean public download/install checks remain pending until
publication.

The Docker allowlist excludes generated Java API pages, test reports and
build artifacts. Inspection of the actual builder COPY layer found only the
35 required manifests, source/resource/launcher and legal files.

The production fingerprint covers 40 files: ten Java source files, 27 C++
headers, the native binding and two Python modules. SHA-256 over sorted
`SHA256  relative-path` records with newline delimiters is
`5818910469b616dd5aec65197bc98dab6211b539015b06b14ae91c9a8cb66980`.
This identifies the 7.2.1 code; it does not relabel historical benchmarks.

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
source manifests and installed-wheel logs. The fresh 7.2.1 Windows run is
[37394450131](https://github.com/asad/SMSD/actions/runs/37394450131). The corrected
7.2.0 Windows run remains separate historical evidence.

Collect the verified wheels and check the complete set before publishing:

```bash
python scripts/collect-release-wheels.py \
  --release-dir dist/release-7.2.1 \
  --wheel-dir build/platform-release/7.2.1/linux \
  --wheel-dir build/platform-release/7.2.1/windows-run-37394450131/wheel
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

The 7.2.1 deadline diagnostic is separate from the historical cross-solver
benchmarks. These historical bounded-search and small-oracle results do
not prove global optimality for arbitrary molecular graphs, nor do macOS
timings establish performance on Linux or Windows.

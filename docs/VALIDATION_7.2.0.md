# SMSD 7.2.0 local validation

This is an unreleased-source validation record, dated 2026-10-05. The baseline
is the 7.1.2 source snapshot `6807f31`. Earlier release artifacts are distinct
from the candidate tested here.

## Environment and scope

- macOS arm64, Apple M5; AppleClang 21.0.0.21000334; C++17.
- JDK 25, Maven 3.9.14, CDK 2.13.
- Search comparison: Python 3.13.14, pybind11 3.1.0, pytest 9.1.1.
- RDKit 2026.09.1, built locally from official `Release_2026_09_1` sources.
- NumPy 2.5.3 from conda-forge for compatible local Accelerate linkage.
- Release wheel: Python 3.14.8, RDKit 2026.03.6 and NumPy 2.5.3 from PyPI.
  This optional interoperability check is separate from the latest-source
  RDKit benchmark comparison.
- CPU-only Python wheel; OpenMP enabled. Metal checks run separately with
  access to the local hardware. CUDA was not run. Linux and Windows release
  status is recorded in the platform section below.

## Results

| Check | Result |
|---|---|
| Clean Java verification | 1,242 passed; 15 opt-in cases skipped; no failures/errors |
| Java artifacts | All four JARs and both CLI launchers report 7.2.0; JARs include exact SMSD LICENSE/NOTICE copies |
| Native Debug CPU | All 12 CTest suites passed; assertions enabled; includes platform portability checks |
| Installed Python comparison wheel | 691 passed; 8 optional/opt-in cases skipped; Python 3.13.14/RDKit 2026.09.1 |
| Repaired Python release wheel | 691 passed; 8 optional/opt-in cases skipped; Python 3.14.8/RDKit 2026.03.6; built from the source distribution |
| Optional C++ RDKit adapter | Enabled C++20 build, install and external consumer passed against RDKit 2026.09.1 |
| Selected Metal regressions | All 3 selected suites passed with hardware access |
| Independent native oracles | 84,096 small-model cases passed |
| Native recursion state | 1,024 additional connected/disconnected McGregor state-validity cases passed |
| ASan/UBSan | Focused objective, enumeration, symmetry, coverage, stereo, invalid-index, deadline and optional-array checks passed |
| Python guide snippets | All executable Python blocks passed |

The oracle count comprises 49,152 objective cases, 3,072 fragment cases,
8,192 direct McSplit/clique cases, 23,104 public MCS cases and 576 enumeration
cases. A separate 576-case coverage-validity oracle and the reproduced
missing-query-edge fixture passed sanitizer checks. Another 1,024 native
McGregor cases check recursive assignment state, separately from those
mathematical oracles. Java adds 1,024 independent
small-graph objective checks. These counts describe the tested models, rather
than a proof of arbitrary-molecule optimality.

The optional adapter check used a private SDK wrapper supplying Eigen and the
actual RDKit source headers, because the local RDKit build's generated export
referenced an uninstalled prefix. RDKit sources and configuration were
unchanged. This validates the adapter build and exported consumer target
within that SDK scope; the conversion limits are listed in [the C++ guide](CPP.md).

Full datasets, standalone benchmark programs and opt-in diagnostic suites are
reported in [the benchmark report](../benchmarks/RESULTS_7.2.0.md), with input
hashes, policies, exclusions and execution status. Their measurements are
separate from the default test-suite counts above.

## Reproduction

Run from the repository root with a dedicated Python environment:

```bash
mvn -B -Dslow.tests.exclude=nothing clean verify \
  org.apache.maven.plugins:maven-source-plugin:3.3.1:jar-no-fork \
  org.apache.maven.plugins:maven-javadoc-plugin:3.6.3:jar
src/scripts/smsd --version
target/appassembler/bin/smsd --version

cmake -S cpp -B build/validation-cpu -DCMAKE_BUILD_TYPE=Debug \
  -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/validation-cpu --parallel 2
ctest --test-dir build/validation-cpu --output-on-failure

python -m build --wheel --no-isolation \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
python -m pip install --no-deps --force-reinstall dist/smsd-7.2.0-*.whl
python -c 'import smsd; print(smsd.__file__, smsd.__version__, smsd.gpu_device_info())'
python -m pytest python/tests -q --import-mode=importlib
```

The explicit import mode prevents tests from loading source wrappers with an
unrelated installed extension. Check the installed path and version before
interpreting results. For the same RDKit comparison, install or build
2026.09.1; installing an older wheel does not reproduce that comparison.

On macOS, the selected hardware tests can be reproduced with:

```bash
cmake -S cpp -B build/validation-metal -DCMAKE_BUILD_TYPE=Debug \
  -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_METAL=ON -DSMSD_BUILD_CUDA=OFF
cmake --build build/validation-metal --parallel 2
ctest --test-dir build/validation-metal -R 'batch_gpu|gpu_domain' --output-on-failure
```

A GPU test skipped for unavailable hardware is not a hardware pass. Resource
budgets are cooperative; preparation and work units can overshoot. Public MCS
mappings carry no optimality certificate or cancellation flag. Native weighted
scores truncate to integer millipoints; Java uses double scores. Exact mapping
canonicalization explicitly rejects incomplete generators or resource limits.

## Release preparation

`scripts/prepare-release.sh` assembles locally validated artifacts and checksum
files. It does not tag, publish or run hosted jobs. The four README files link
to the measured report and remove unsupported historical speed and quality
claims. Local instruction/configuration files and private review logs are
excluded from public sources and artifacts.

The compact Python release targets CPython 3.14 wheels for Linux x86_64,
macOS arm64 and Windows x86_64, plus a source distribution. A platform target
does not imply a completed execution check. The macOS wheel's OpenMP library
is bundled; linkage and
all RECORD hashes were checked. Execution was validated on macOS 27.0.1;
the deployment tag does not certify a separate macOS 26 execution. The Maven
release profile passed local verification with signing skipped. Publication
still requires a PyPI token and interactive GPG signing; see
[publishing commands](PUBLISHING.md).

## Platform release checks

The portable Java 25 JAR and CLI archive contain Java code and Unix/Windows
launchers. Their macOS execution does not certify each Windows launcher or an
OS installer. Native installers are outside the compact asset set.

| Python target | Execution status |
|---|---|
| macOS 26+ arm64 | Installed wheel tested on macOS 27.0.1; counts above |
| Linux x86_64, glibc 2.28+ | All 12 native Debug suites; 691 Python tests passed, 8 skips |
| Windows x86_64 | Native build and execution checks pending |

Following the portability edits, all 12 native CPU suites pass on macOS and
the rebuilt CPython 3.14.8 wheel again passes 691 tests with 8 skips. Its SHA-256
is `5d974fe90f5b973d72d56fb108c9f500a320c2fa3ac386b59e2834ca1ab3054a`.
The additional native suite checks the standalone depiction header and
Unicode MOL/SDF round-trips. Release-wheel hashes differ from the comparison
wheels recorded in the benchmark report.

Linux checks use the local `manylinux_2_28_x86_64:2026.06.04-1` image under
x86_64 emulation on the arm64 Mac: glibc 2.28, GCC 14.2.1, CPython 3.14.5,
pybind11 3.1.0 and RDKit 2026.03.6. Auditwheel repair bundles `libgomp` and
emits both manylinux 2.27 and 2.28 compatibility tags; execution was checked
on the glibc 2.28 container. Twine checks, installed-module origin, CPU/OpenMP
functionality, full Python tests and artifact hashes pass. Timings from
emulation do not extend the macOS benchmark comparisons.
The final Linux wheel SHA-256 is
`9ec857dd14969eb05316428f4836b2604c87a3ad512a612a770f2a818ae22672`.

The initial Linux run exposed a batch test's one-second hardware speed
threshold. The test now checks all 1,000 expected substructure results and
fingerprint lengths/content/repeated-input consistency, while reporting
elapsed time. The corrected suite passed on Linux and macOS; the other 11
Linux suites passed in the initial run. Production search code was unchanged.

The platform review corrected a standalone depiction header's dependence on
`M_PI`. A regression test includes that header with the macro undefined and
checks ring geometry and SVG output. Native MOL/SDF APIs now interpret filenames
as UTF-8, and MSVC builds select UTF-8 source and executable character sets.
A MinGW-w64 GCC 16.2.0 Windows cross-build compiles and links all 12 native
targets and an installed CMake consumer; this is compiler evidence rather than Windows runtime
validation. Windows publication requires its native wheel and tests.

All three release wheels must match the source distribution's Python wrappers,
C++ headers and license copies. `scripts/collect-release-wheels.py` also checks
wheel RECORD hashes and requires all three platforms by default. The Linux
repair includes the GCC OpenMP runtime license and runtime-library exception;
macOS includes the LLVM OpenMP license.

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
  access to the local hardware. CUDA and other operating systems were not run.

## Results

| Check | Result |
|---|---|
| Clean Java verification | 1,242 passed; 15 opt-in cases skipped; no failures/errors |
| Java artifacts | All four JARs and both CLI launchers report 7.2.0; JARs include exact SMSD LICENSE/NOTICE copies |
| Native Debug CPU | All 11 CTest suites passed; assertions enabled |
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

The compact Python asset set contains one CPython 3.14/macOS 26+ arm64 wheel
and a source distribution. The wheel's OpenMP library is bundled; linkage and
all RECORD hashes were checked. Execution was validated on macOS 27.0.1;
the deployment tag does not certify a separate macOS 26 execution. The Maven
release profile passed local verification with signing skipped. Publication
still requires a PyPI token and interactive GPG signing; see
[publishing commands](PUBLISHING.md).

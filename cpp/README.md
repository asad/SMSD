# SMSD C++ module

The C++17 core is header-only and provides molecular graphs, substructure
and MCS search, fingerprints, parsing, stereo/CIP and depiction. Public headers
are in `include/smsd/`, native tests in `tests/`, and the Python extension
bindings in `bindings/pybind11/`. Java is not required for the native core.

See [the C++ guide](../docs/CPP.md) for APIs, matching contracts, installation
and RDKit integration. The current source targets 7.2.1; fresh release checks
passed all 12 native Debug suites on macOS arm64, emulated Linux x86_64 and
native Windows x86_64. These CPU/OpenMP builds use the same frozen source.
See [validation](../docs/VALIDATION_7.2.1.md) for toolchains and execution scope;
the headers are available in the GitHub release.

## Build and test

Run from the repository root with CMake 3.18+ and a C++17 compiler. The CPU
configuration works with single-configuration and Visual Studio generators:

```text
cmake -S cpp -B build/cpp -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_OPENMP=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF -DSMSD_WITH_RDKIT=OFF
cmake --build build/cpp --config Debug --parallel 4
ctest --test-dir build/cpp --build-config Debug --output-on-failure
```

Assertions stay enabled in the test targets. OpenMP is detected when
available; otherwise batch operations use the sequential fallback. Install
the headers and exported `smsd::smsd` CMake target with:

```bash
cmake --install build/cpp --config Debug --prefix "$PWD/build/cpp-install"
```

Include the public API through `smsd/smsd.hpp`. An installed consumer can
use `find_package(smsd 7.2 CONFIG REQUIRED)` and link `smsd::smsd`.

## Optional builds

- `SMSD_BUILD_PYTHON=ON` requires Python development headers and pybind11.
  Python package builds run from the repository root using its single
  `pyproject.toml`; see [the Python module](../python/README.md).
- `SMSD_WITH_RDKIT=ON` adds the `smsd::smsd_rdkit` adapter target and requires
  a compatible RDKit development installation and C++20. Its metadata
  conversion limits are documented in [the C++ guide](../docs/CPP.md).
- `SMSD_BUILD_METAL=ON` requires macOS and Metal; `SMSD_BUILD_CUDA=ON`
  requires the CUDA toolkit. Both also accept `AUTO`. These optional paths
  require separate hardware validation; the release wheels use CPU/OpenMP.

Bounded MCS searches can return a valid mapping without proving a global
optimum. The retained [7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md)
records its original source versions and measurement scope. The 7.2.1 seed
deadline regression is recorded separately; the full corpus comparison has
not been rerun for its patched source.

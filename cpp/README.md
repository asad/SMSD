# SMSD Pro C++ library

SMSD Pro provides a header-only C++17 library for molecular graphs,
substructure search, maximum common substructure (MCS) search, SMARTS,
fingerprints, MOL and SDF input and output, stereochemistry, CIP descriptors
and SVG depiction. The core needs only the C++ standard library and a threads
implementation; it does not need Java or RDKit.

| Path | Contents |
| --- | --- |
| `include/smsd/` | Public headers |
| `tests/` | Native test suites |
| `bindings/pybind11/` | Python extension bindings |
| `src/rdkit_adapter.cpp` | Optional RDKit adapter |

This page covers building and integration. The [C++ guide](../docs/CPP.md) is
the API reference, and the [changelog](../CHANGELOG.md) describes what the
release provides.

## Quick start

Save this program as `example.cpp` in the repository root:

```cpp
#include <iostream>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");

    const auto embedding = smsd::Substructure::find(query).withTimeout(1000).match(target);
    const auto common = smsd::MCS::find(query).withTimeout(1000).match(target);
    if (embedding.matches.empty() || common.matches.empty() || !common.optimal) return 1;

    std::cout << "Substructure atoms: " << embedding.matches.front().atomMapping.size() << '\n'
              << "MCS atoms: " << common.matches.front().atomMapping.size() << '\n';
}
```

Compile and run it on macOS or Linux:

```bash
c++ -std=c++17 -O2 -pthread -I cpp/include example.cpp -o example
./example
```

The program reports six substructure atoms and six MCS atoms. Mappings use
zero-based query-to-target atom indices.

`Substructure` finds full-query embeddings of a molecule or SMARTS query. `MCS`
defaults to the connected maximum common induced subgraph (MCIS), which
maximises mapped atoms and then bonds while preserving both edges and non-edges
between mapped atoms; the maximum common edge subgraph (MCES) objective
maximises common bonds instead. Both classes provide `match`, `matchAll` and
`matchMany`, and every result reports its stopping status and whether it is
optimal and complete. Check these fields before treating a bounded result as a
proved optimum or a complete enumeration.

## Build and test

Use a C++17 compiler and CMake 3.20 or later for these commands, which run from
the repository root (`ctest --test-dir` needs CMake 3.20); the library itself
supports CMake 3.18 or later. They work with single-configuration and Visual
Studio generators:

```text
cmake -S cpp -B build/cpp -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF
cmake --build build/cpp --config Debug --parallel 4
ctest --test-dir build/cpp --build-config Debug --output-on-failure
```

Assertions stay enabled in the test targets, including Release builds.

| CMake option | Default | Effect |
| --- | --- | --- |
| `SMSD_BUILD_TESTS` | `ON` | Builds the native test suites |
| `SMSD_HEADER_CHECK_WERROR` | `OFF` | Makes warnings errors in the header self-containment test (`-Werror`, or `/WX` with MSVC) |
| `SMSD_BUILD_OPENMP` | `ON` | Uses OpenMP, when found, for batch fingerprint and screening operations; without it they run sequentially |
| `SMSD_BUILD_PYTHON` | `OFF` | Builds the Python extension; requires Python development headers and pybind11 |
| `SMSD_WITH_RDKIT` | `OFF` | Builds the `smsd_rdkit` adapter; requires an RDKit development installation and C++20 |
| `SMSD_BUILD_METAL` | `AUTO` | Metal screening on macOS; `ON` fails configuration elsewhere |
| `SMSD_BUILD_CUDA` | `AUTO` | CUDA screening when the CUDA toolkit and compiler are found; `ON` requires them |

Python packages are built from the repository root with its `pyproject.toml`;
see the [Python module](../python/README.md). The RDKit adapter's conversion
rules are described under
[Optional native adapter](../docs/CPP.md#optional-native-adapter). The Metal
and CUDA libraries are not part of the installed package, so GPU screening is
available only to programs built in the same CMake project. The release Python
wheels are built without Metal or CUDA.

## Install and use with CMake

Install the headers and the exported `smsd::smsd` target:

```text
cmake --install build/cpp --config Debug --prefix build/cpp-install
```

No compiled library is needed for the core. The imported target supplies C++17,
threads and, when configured, OpenMP; with `SMSD_WITH_RDKIT=ON` the package also
exports `smsd::smsd_rdkit`. To use the package, save the quick-start program as
`example/main.cpp` and this file as `example/CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.18)
project(smsd_example LANGUAGES CXX)
find_package(smsd 8.0.0 CONFIG REQUIRED)
add_executable(smsd_example main.cpp)
target_link_libraries(smsd_example PRIVATE smsd::smsd)
```

Build from the repository root, replacing `<install-prefix>` with the absolute
path of `build/cpp-install`:

```text
cmake -S example -B build/example -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<install-prefix>"
cmake --build build/example --config Release
```

Run `build/example/smsd_example` on macOS or Linux, or
`build/example/Release/smsd_example.exe` with a Visual Studio generator. The
package uses same-major-version compatibility, so `find_package(smsd 8.0.0)`
accepts any later 8.x installation.

## Documentation

- [C++ guide](../docs/CPP.md): API reference, search and chemistry options,
  fingerprints, file input and depiction.
- [Installation guide](../docs/HOWTO-INSTALL.md): platform requirements and
  source builds.
- [Validation record](../docs/VALIDATION_8.0.0.md): release checks for 8.0.0.
- [Benchmark guide](../benchmarks/README.md): workloads shared by C++, Java and
  Python, and mapping validation.
- [Changelog](../CHANGELOG.md): what the release provides.

## Licence

Distributed under the [Apache License 2.0](../LICENSE). See [NOTICE](../NOTICE)
for attribution and trademark terms.

# SMSD C++

SMSD 7.2.2 provides a C++17 header-only core for molecular graphs,
substructure and MCS search, fingerprints, SMARTS, MOL/SDF I/O, stereo/CIP
and depiction. It uses the standard library and does not require Java or
RDKit. Public headers are in `include/smsd/`, native tests in `tests/`, and
Python extension bindings in `bindings/pybind11/`.

## Quick start

Save this program as `example.cpp` in the repository root:

```cpp
#include <iostream>
#include "smsd/smsd.hpp"

int main() {
    const auto query = smsd::parseSMILES("c1ccccc1");
    const auto target = smsd::parseSMILES("c1ccc(O)cc1");
    const smsd::ChemOptions chemistry;
    smsd::MCSOptions options;
    options.timeoutMs = 1000;

    const auto embedding = smsd::findSubstructure(query, target, chemistry, 1000);
    const auto mapping = smsd::findMCS(query, target, chemistry, options);
    if (embedding.size() != 6 || mapping.size() != 6 ||
        !smsd::validateMapping(query, target, mapping, chemistry).empty()) {
        return 1;
    }
    std::cout << "Substructure atoms: " << embedding.size() << '\n'
              << "MCS atoms: " << mapping.size() << '\n';
}
```

Compile and run on macOS or Linux:

```bash
c++ -std=c++17 -O2 -I cpp/include example.cpp -o example
./example
```

The program reports six substructure atoms and six MCS atoms. Mappings use
zero-based query-to-target atom indices. See [the C++ guide](../docs/CPP.md)
for chemistry options, mapping contracts and more examples.

## Build, test and install

Run from the repository root with CMake 3.20+ and a C++17 compiler. The CPU
configuration works with single-configuration and Visual Studio generators:

```text
cmake -S cpp -B build/cpp -DCMAKE_BUILD_TYPE=Debug -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_OPENMP=ON -DSMSD_BUILD_METAL=OFF -DSMSD_BUILD_CUDA=OFF -DSMSD_WITH_RDKIT=OFF
cmake --build build/cpp --config Debug --parallel 4
ctest --test-dir build/cpp --build-config Debug --output-on-failure
```

Assertions stay enabled in the test targets. OpenMP is used when detected;
otherwise batch operations run sequentially. Install the headers and
exported `smsd::smsd` CMake target with:

```text
cmake --install build/cpp --config Debug --prefix build/cpp-install
```

No compiled C++ core library is needed. To use the installed package, save
the quick-start program as `example/main.cpp` and this file as
`example/CMakeLists.txt`:

```cmake
cmake_minimum_required(VERSION 3.20)
project(smsd_example LANGUAGES CXX)
find_package(smsd 7.2.2 CONFIG REQUIRED)
add_executable(smsd_example main.cpp)
target_link_libraries(smsd_example PRIVATE smsd::smsd)
```

Build from the repository root, replacing `<install-prefix>` with the
absolute path to `build/cpp-install`:

```text
cmake -S example -B build/example -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<install-prefix>"
cmake --build build/example --config Release
```

Run `build/example/smsd_example` on macOS/Linux, or
`build/example/Release/smsd_example.exe` with a Visual Studio generator.
The imported target supplies C++17 and any configured OpenMP dependency.

## Optional builds

- `SMSD_BUILD_PYTHON=ON` requires Python development headers and pybind11.
  Python package builds run from the repository root using its single
  `pyproject.toml`; see [the Python module](../python/README.md).
- `SMSD_WITH_RDKIT=ON` adds the `smsd::smsd_rdkit` adapter target and requires
  a compatible RDKit development installation and C++20. Its metadata
  conversion limits are documented in [the C++ guide](../docs/CPP.md#rdkit-integration).
- `SMSD_BUILD_METAL=ON` requires macOS and Metal; `SMSD_BUILD_CUDA=ON`
  requires the CUDA toolkit. Both also accept `AUTO`. These optional paths
  require separate hardware validation; the release wheels use CPU/OpenMP.

GPU support requires a source build. The MSI, DMG and DEB packages provide
the Java CLI with a bundled runtime; they do not install the C++ headers
or Python package.

## Validation

The 7.2.2 release checks cover 12 native Debug suites and 691 installed Python
tests with 8 optional skips on macOS arm64, Linux x86_64 and Windows x86_64.
Platform versions, Linux emulation and reused evidence are detailed in
[release validation](../docs/VALIDATION_7.2.2.md).

Bounded MCS searches can return a valid mapping without proving a global
optimum. The [7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) retains
its measured source versions and scope. The full corpus comparison has not
been rerun for 7.2.2, so it does not establish current performance against
other toolkits.

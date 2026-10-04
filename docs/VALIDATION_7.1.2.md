# SMSD Pro 7.1.2 local validation

Validated on 4 October 2026 using macOS arm64, AppleClang 21, Java 25.0.2,
Maven 3.9.14, CPython 3.13, CDK 2.13 and RDKit 2026.3.6.

```bash
MACOSX_DEPLOYMENT_TARGET=26.0 SMSD_RELEASE_PYTHON=.venv/bin/python scripts/prepare-release.sh
```

The final preflight passed 1,192 Java tests (15 opt-in benchmark cases skipped),
all six C++ CPU suites with assertions enabled, and 603 installed-wheel Python
tests (6 optional cases skipped). The source distribution was extracted to a
fresh directory and used to build the tested wheel. Delocate bundles OpenMP;
the extension uses `@loader_path/.dylibs/libomp.dylib` rather than a local
Homebrew path. The wheel requires CPython 3.13, macOS 26+, and arm64.

The separate Metal/OpenMP batch test passed on the Apple M5 GPU. CUDA and
native Windows execution were not tested. Additional platform wheels, native
installers, PyPI publishing, and container publishing are outside this draft.

The 23 new Java regression cases all fail against the original sources and pass
with these fixes. Another 6,000 deterministic comparisons agree with exhaustive
injective-mapping oracles across matching profiles and search engines.

The general matcher agrees with an exhaustive oracle on all 33,868 graphs with
0–6 vertices and 3,072 seeded graphs with 7–12 vertices. These cases run in
`smsd_general_matching_tests`. Targeted AddressSanitizer/UndefinedBehaviorSanitizer
checks for matching, deadlines, kekulization, MCS, and canonical hashing pass.
Canonical SMILES and hashes are unchanged for a 13-molecule before/after corpus.
Relocated CMake installations compile and run independent consumers with OpenMP
enabled and disabled, while propagating C++17.

## C++ primitive measurements

These are three-run median synthetic measurements with `clang++ -std=c++17 -O3`
against original commit `52733cb` and the updated headers on the same machine.
They measure the individual primitives, not whole-application throughput.

| Primitive | Original | Updated | Ratio |
|---|---:|---:|---:|
| NLF3, 100-atom chain, 50,000 constructions | 11.10 ms | 6.84 ms | 1.6× |
| NLF3, 1,000-atom chain, 50,000 constructions | 45.95 ms | 8.22 ms | 5.6× |
| Clique, fixed 80-vertex graph, 20 searches | 64.95 ms | 24.16 ms | 2.7× |

```bash
clang++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_search_primitives.cpp -o build/search-primitives
build/search-primitives
```

Compile the same harness against an isolated checkout of `52733cb` for the
baseline. Clique results retain the same maximum size (9). The NLF change also
corrects its semantics to cumulative neighborhoods, so its label counts differ
where required for safe substructure pruning.

# SMSD Pro 7.1.2 historical local validation

This records the earlier 7.1.2 preparation against baseline `52733cb`. It is
not the validation result for current 7.2.0 Unreleased changes. Current
measurements and coverage are linked from the
[7.2.0 benchmark report](../benchmarks/RESULTS_7.2.0.md) and
[algorithm review](ALGORITHM_REVIEW.md).

Validated on 4 October 2026 using macOS arm64, AppleClang 21, Java 25.0.2,
Maven 3.9.14, CPython 3.13, CDK 2.13, Jackson 2.21.7 and RDKit 2026.3.6.

GitHub's 11 open advisories for Jackson 2.20.0 require patched versions up to
2.21.7; the release resolves Jackson databind/core to 2.21.7. Java tests and
CLI JSON output were revalidated after this dependency change. Alerts on the
default branch remain until the changes are merged.

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

The earlier review recorded targeted Java baseline failures and passing fixed
regressions. Its 6,000 deterministic comparisons agreed with exhaustive
injective-mapping oracles across matching profiles and search engines.

The general matcher agrees with an exhaustive oracle on all 33,868 graphs with
0–6 vertices and 3,072 seeded graphs with 7–12 vertices. These cases run in
`smsd_general_matching_tests`. Targeted AddressSanitizer/UndefinedBehaviorSanitizer
checks for matching, deadlines, kekulization, MCS, and canonical hashing pass.
Canonical SMILES and hashes are unchanged for a 13-molecule before/after corpus.
Relocated CMake installations compile and run independent consumers with OpenMP
enabled and disabled, while propagating C++17.

## Historical C++ primitive measurements

These are three-run median synthetic measurements with `clang++ -std=c++17 -O3`
against original commit `52733cb` and the updated headers on the same machine.
They measure the individual primitives, not whole-application throughput, and
are not the current 7.2.0 baseline/candidate comparison.

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

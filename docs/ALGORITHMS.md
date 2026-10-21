# SMSD search algorithms

The C++, Java and Python searches share one design: a prepared query owns its
chemistry and resource options, and every `match`, `matchAll` or `matchMany`
call runs an independent search against one target. The C++ and Java
common-bond (MCES) searches are checked against each other node for node; the
substructure and MCIS searches are compared across the languages for equal
proved optima and, where an enumeration completes, equal match sets; and the
Python package calls the native search.
The language guides ([Python](PYTHON.md), [C++](CPP.md) and
[Java](JAVA.md)) document the public API; this page summarises the search
stages, the meaning of a result and the differences between the languages.

## Search stages

Preparing a query snapshots its graph and resolves its chemistry options
once; SMARTS queries are compiled once and reused for every target.

Substructure search finds embeddings of the query in the target. Non-induced
matching requires every query bond in the target; induced matching also
forbids extra target bonds between mapped atoms. Enumeration returns every
embedding up to the requested limit, including symmetric self-mappings.

MCIS search maximises the mapped atoms, then bonds, or the supplied query atom
weights; MCES search maximises the common bonds, then the atoms they touch.
Every candidate is checked against the chemistry options before it is
reported, connected components are chosen by the requested objective, and
fragment limits constrain the reported components.

Canonical mapping selection chooses a representative mapping under the
symmetries of both molecules, and raises an explicit error when it cannot do
so within its limits.

## Results and budgets

Every search reports `status`, `optimal` and `complete` separately from its
matches. `Completed` means the search ran to its natural end; `TimedOut`,
`Cancelled` and `LimitReached` describe why it stopped early. A stopped search
may still return a valid incumbent without proving the optimum, and a capped
enumeration may prove the objective while leaving the enumeration incomplete.
An empty incomplete result does not prove absence. A negative timeout selects
an adaptive budget from the graph sizes; zero stops immediately; a node limit
of zero leaves the node count unrestricted.

The regression suites check small-graph results against independent
brute-force oracles rather than recorded outputs: every simple graph in the
smallest size classes for substructure, MCS and enumeration, seeded random
graphs for the Java solver and bounds, exhaustive matrices for assignment and
matching, and focused chemistry fixtures for connectivity, fragments, weights,
aromaticity, isotopes, chirality, rings and tautomers. Exhaustive oracles prove
the tested cases; a time or node limit on a large molecule can leave a valid
incumbent without a proof of global optimality.

## Cross-language differences

These differences are deliberate; each follows the conventions of its
language.

| Topic | C++ | Java | Python |
| --- | --- | --- | --- |
| Enum constants | `SearchStatus::Completed`, `MCS::Connectivity::Connected`; `ChemOptions` enums such as `BondOrderMode::STRICT` | `SearchResult.Status.Completed`, `MCS.Connectivity.Connected`; `ChemOptions` enums such as `BondOrderMode.STRICT` | `SearchStatus.COMPLETED`, `MCS.Connectivity.CONNECTED`, `BondOrderMode.STRICT` |
| Time budget and elapsed time | `timeoutMs` option, `elapsedMillis` result | `timeoutMillis` option, `elapsedMillis()` result | `timeout_ms` option, `elapsed_ms` result |
| Per-match score | `SearchResult::score` only | `Match.score()` as well as `SearchResult.score()` | `SearchResult.score` only |
| Witness accessors | Inspect `matches` directly | `found()` | `found` and `best` properties |
| Fluent option helpers | `withOptions`, `withTimeout` | `withOptions`, `withTimeout` | `with_options`, `with_timeout`, `with_objective`, `with_fragments`, `with_atom_weights`, `with_induced` |
| SMILES convenience | `MCS::find(std::string)` and `Substructure::find(std::string)` parse SMILES | Graphs and CDK containers only | `find(str)` parses SMILES; RDKit molecules are accepted |
| Batch searches | `matchMany(targets, numThreads = 0)`; zero uses every processor | `matchMany(List)` runs sequentially; `matchMany(List, numThreads)` uses up to `numThreads` threads, zero for one per processor | `match_many(targets, num_threads=0)` |
| Tuning switches | None | Deprecated `ChemOptions` tuning switches that never change a result | None |
| Cancellation | Shared `CancellationToken` | `BooleanSupplier`; thread interruption also cancels | Shared `CancellationToken` |

Batch results are identical whichever thread count is used, provided no
search reaches its time limit; completed and node-limited searches do not
depend on the thread count.

## Reproduce locally

Run from the repository root:

```sh
mvn -f java/pom.xml -B -Dslow.tests.exclude=nothing clean verify

cmake -S cpp -B build/algorithm-tests \
  -DCMAKE_BUILD_TYPE=Debug \
  -DSMSD_BUILD_TESTS=ON -DSMSD_BUILD_PYTHON=OFF \
  -DSMSD_BUILD_CUDA=OFF -DSMSD_BUILD_METAL=AUTO
cmake --build build/algorithm-tests --parallel 4
ctest --test-dir build/algorithm-tests --output-on-failure
```

The Metal target is conditional on macOS. On supported hardware its native
domain test runs; unavailable hardware produces a CTest skip. Release checks
are recorded in the [validation record](VALIDATION_8.0.0.md) as they complete;
a timing ratio does not establish a speed-up when mappings differ in validity,
objective, chemical policy or completion status.

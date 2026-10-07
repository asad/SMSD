# Search algorithm review

This review records the **7.2.0 source candidate** against the immutable 7.1.2
source snapshot `6807f31`. Those search fixes carry forward into 7.2.2; current release checks are in
[7.2.2 validation](VALIDATION_7.2.2.md). Latest cross-engine measurements, policies, input
hashes and quality/cancellation outcomes are in the
[benchmark report](../benchmarks/RESULTS_7.2.0.md). Earlier review results below
are historical and are identified separately.

## Reviewed objective, chemistry and cache fixes

Tautomer-aware Java and native matching preserves element identity when atom
matching is enabled. Tetrahedral matching compares normalized R/S descriptors
and fully mapped ligand permutations rather than raw traversal tags; the Java
graph retains CDK's exact ligand ordering. Opposite resolved configurations
exclude the center even in partial mappings. Unspecified stereo retains its
existing wildcard behavior.

Signed weights can make an identity mapping worse than a proper subset.
Small objective searches use positive-weight or bond-count upper bounds;
component/fragment selection and constrained target choice use the selected
objective. Java weight comparisons retain double precision. Native scores
retain their documented integer-millipoint range and input validation.
Constrained Java batches preserve original target chemistry and indices while
masking claimed atoms, support Builder graphs, and apply the per-pair timeout
without modifying caller options.

Canonical mapping selection explores generator-orbit closure rather than a
greedy local descent. Incomplete generators or exceeded resource bounds raise
an explicit error; internal enumeration uses raw keys if symmetry work cannot
complete. Weighted enumeration does not treat an unweighted graph automorphism
as evidence that user-weighted mappings are equivalent.

Java's molecule cache has weak values as well as weak keys: a cached MolGraph
otherwise retains its own CDK key. New focused Java coverage passed 77 cases,
including 25 new regressions and 1,024 independent three-vertex comparisons
across signed/positive weights, bond objectives, connectedness and induced
matching. Seven of the initial eight defect fixtures failed against the
baseline; valid positive-weight behavior was retained as an additional check.
The focused suite also covers CDK bond-storage reordering, completed ligand
permutations, Builder winding, scoped exclusions and constrained batch budgets.
Three of four additional chemical-symmetry fixtures failed before their fixes;
the valid methyl symmetry away from an annotated center remained available.
Large/skipped or node-capped generator capture is marked incomplete, discovered
permutations preserve chemical properties, and orbit pruning uses only proven
subgroup equivalences. Another 20 existing orbit, canonical-label and SMILES
tests passed after this change; canonical labels/hashes were left unchanged.

## Earlier correctness changes (review baseline `49d7303`)

Java and C++ MCS searches preserve the caller's query direction for
non-induced matching and keep atom weights attached to query indices.
Partial induced subgraphs use conservative chemical-label bounds: vertices
may omit neighbors, so their full-graph degrees cannot bound an MCS.
Candidates are validated and postfiltered before becoming incumbents or
triggering an upper-bound exit. Identity results obey connectivity and
fragment options. Connected components are selected by the requested
objective, including atom weights or conserved bonds.
After reversing a search, both languages reapply the original query's filters.
A fused ten-carbon query ring system against a six-carbon target ring cannot
retain a partial ring when complete rings are requested; its two-atom N–N
pendant remains valid. This is covered in Java, C++ and the Python binding.

The Java product graph includes permitted query non-edges and partial
vertices. McSplit allows sufficient skip depth, preserves non-induced
candidate partitions, and restricts symmetry pruning to valid contexts.
Recursive extension branches preserve their candidate buffers. MCS
enumeration applies default options and the same validity filters as the
single-result search. Similarity bounds honor chemical matching options.
Expired Java and C++ budgets stay expired during recursive unwinding. Java
clique search checks deadlines around costly coloring and pivot work;
orientation probes, recovery and ring retries share the caller's total budget.
Valid connected incumbents can extend even below the near-optimal threshold.
For graphs with 51–100 atoms, Java tries compatible single-atom anchors within
50 ms and one eighth of the remaining budget. This retains a valid macrolide
seed before expensive search; the 51/52-atom erythromycin/azithromycin literals
improve from six to 28 mapped atoms under the one-second request. A separate
49-atom connected witness establishes that the stress test's 25-atom lower
bound is valid. The bounded search still does not prove the large-pair optimum.

Python `find_mcs(strategy="auto")` retains the lightweight route for its
supported default policy, timeout, ring matching and strict/any bond mode.
Requests for induced, isotope, chirality, complete-ring, disconnected or
other unsupported lightweight constraints use the native solver. Explicit
`strategy="lightweight"` rejects unsupported options with `ValueError`.
Unknown native option names surface their binding error.

C++ substructure enumeration includes symmetric self mappings (two for a
three-carbon path, twelve for benzene). A disconnected cycle plus path no
longer enters the simple-path shortcut. Unrestricted bond order still
enforces strict aromaticity. Accelerator domains are conservative supersets
and are refined with the CPU isotope, chirality, ring and tautomer policies.
The portable GPU test double checks dispatch and CPU agreement; the native
test additionally requires successful kernel execution and checks its full
domain bitmatrix against an independent basic-compatibility oracle.

The native clique helpers enforce mapped query bonds, select the largest
connected component, normalize results before applying caps, and check an
already expired deadline. Rectangular assignment validates matrix shape and
finite values, and rejects reduced-cost overflow. C++ atom-weight validation
rejects nonfinite values and unrepresentable integer millipoint scores.

## Earlier independent validation

The regression suites use independent brute-force oracles rather than
expected results copied from the implementation:

| Area | Coverage |
|---|---|
| C++ substructure | 23,104 cases covering every simple graph with zero to four vertices, induced/non-induced matching, VF2/VF2++, exact mapping sets, existence and injectivity |
| C++ MCS | 23,104 cases over the same small graph family, induced/non-induced and connected/disconnected objectives, optimum size and complete mapping validity |
| C++ MCS enumeration | 576 exhaustive small-graph cases checking every returned mapping against the requested objective and validity constraints |
| Java MCS | 400 seeded random graph comparisons with three to six vertices, induced/non-induced matching, optimum size and mapping validity |
| Java bounds | All 4,096 pairs of four-vertex simple graphs checked against an exhaustive induced-MCS oracle |
| Native clique | All 33,868 simple graphs through six vertices, plus 256 random graphs with seven to ten vertices and result-cap/incumbent checks |
| Native helper search | 825 substructure graph pairs and 128 pipeline pairs in both orientations with three result caps |
| Rectangular assignment | 3,200 seeded matrices in both rectangular orientations, each tested with three unmatched penalties against brute force |
| General matching | All 33,868 simple graphs through six vertices and 3,072 seeded larger graphs checked against an independent matching oracle |
| Python API regressions | 30 cases covering solver selection, unsupported options, induced/isotope results, connected/full-fragment salt mappings and caller-direction complete rings |

Focused chemistry and option fixtures cover identity connectivity, fragment
precedence, unequal query direction, query weights, objective-aware component
selection, aromatic bonds, isotopes, chirality, ring policies and tautomers.
The padded nine-atom connected fixture exercises heuristic search beyond the
new eight-atom exact-search cutoff.
An independent eleven-atom biphenyl/diphenylmethane mapping verifies connected
seed extension. Chemistry fixtures distinguish strict keto/enol matching
(two atoms) from tautomer-aware matching (four), and connected alanine
enantiomer matching (three atoms) from disconnected matching (five).
Round-trip salt checks explicitly permit disconnected fragments; direction
stability checks use symmetric induced topology. Parallel Java reuse retains
its per-task guards and allows the combined five 15-second search budgets.
Outer wall-clock guards allow graph construction and assertion overhead in
addition to the unchanged one-second macrolide and ten-second PDE5 search budgets.
Baseline failures were reproduced before applying the fixes. The C++ MCS
baseline had 5,280 size mismatches in the small-graph oracle; the updated
implementation has zero. Focused AddressSanitizer/UndefinedBehaviorSanitizer
runs cover the native algorithm regressions.

Historical validation recorded on 2026-10-05, before the reviewed 7.2.0 changes:

| Check | Result |
|---|---|
| Maven `clean verify`, including slow suites | 1,216 passed, 15 skipped (1,231 total); zero failures or errors; JAR and CLI packaging completed |
| C++ CTest suites | All 11 portable suites passed; the native Metal domain test also passed with hardware access |
| Repaired macOS arm64 CPython 3.13 wheel, built from a fresh source distribution | 633 passed, 6 skipped; both Python and extension imports verified inside the extracted wheel |
| Focused native sanitizer checks | AddressSanitizer and UndefinedBehaviorSanitizer passed |

The native Metal test was rerun with hardware access after the sandboxed
CTest run reported the GPU unavailable. It verifies successful kernel
execution and the resulting compatibility matrix, so CPU fallback does not
count as GPU validation.

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
domain test runs; unavailable hardware produces a CTest skip. The portable
domain test always runs. CUDA, Windows and Linux native builds were not
validated in this local review. Hosted workflows require manual dispatch.

The existing [release preparation script](../scripts/prepare-release.sh)
builds source distributions, repairs macOS wheels and validates the installed
package. The earlier review validated candidate wheels under
`build/algorithm-ringfix-wheel-assets`. Current source and release preparation
target 7.2.2; neither those older packages nor a 7.1.2 version label identifies
the current changes. The reproduction commands use the new `java/` module
layout; the recorded counts and benchmark measurements retain their historical
source scope.

## Performance evidence

Use the [current benchmark report](../benchmarks/RESULTS_7.2.0.md) for measured
results. Primitive timings, cross-engine molecular searches, conversion costs,
subprocess startup and tautomer feature runs are distinct workloads. A timing
ratio does not establish a speedup when mappings differ in validity, objective,
chemical policy or completion status.

Reproducible primitive sources include
[assignment/matching](../benchmarks/benchmark_assignment_matching.cpp),
[matcher setup](../benchmarks/benchmark_substructure_setup.cpp),
[component filtering](../benchmarks/benchmark_mcs_components.cpp) and
[search primitives](../benchmarks/benchmark_search_primitives.cpp).
Rectangular assignment uses `O(min(m,n)² max(m,n))` time with `O(m+n)`
auxiliary space. Component filtering traverses graph neighbors instead of
scanning every pair of mapped vertices. These implementation properties do
not by themselves establish representative molecular-search throughput.

Large MCS searches retain bounded heuristic stages. Exhaustive oracles prove
the tested cases; a time or node limit can leave a valid incumbent without a
proof of global optimality.

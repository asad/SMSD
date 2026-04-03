# SMSD Optimization Roadmap

This note is the repo-specific performance plan for SMSD `6.9.x`.
It is intentionally tied to the current code layout rather than generic
"make Python faster" advice.

## Current Read

The maintained strict leaderboard already shows the split clearly:

- Python overhead matters for microsecond-scale screening paths.
- The worst outliers are native exact-MCS search-depth problems.
- Handwritten assembly is not the next lever.

Current evidence from the codebase:

- Heavy Python bindings already release the GIL in
  `cpp/bindings/pybind11/smsd_bindings.cpp`.
- Batch APIs already exist in `cpp/include/smsd/batch.hpp` and
  `python/smsd/__init__.py`.
- The exact MCS core already uses bitsets, `popcount`, orbit data, LFUB,
  McSplit, Bron-Kerbosch, and McGregor extension in
  `cpp/include/smsd/mcs.hpp`.
- VF2++ already uses flat bit domains and AC-3 pruning in
  `cpp/include/smsd/vf2pp.hpp`.

That means the next wins should come from better API shape, cache reuse, and
stronger pruning before any architecture-specific low-level work.

## Priority Order

1. Measure the real split between binding cost, setup cost, and kernel cost.
2. Expose better reusable Python objects for repeated workloads.
3. Reuse precomputed graph invariants across batches and threads.
4. Tighten exact-MCS pruning on the known hard pairs.
5. Reduce container churn in the recursive core.
6. Only then tune bitset kernels further with intrinsics if profiling justifies it.

## Phase 0: Instrumentation First

Before changing algorithms, split end-to-end time into:

- parse / object conversion
- lazy graph prewarm
- compatibility graph or domain build
- core search
- Python result materialization

Primary targets:

- `benchmarks/benchmark_python.py`
- `cpp/include/smsd/mcs.hpp`
- `cpp/include/smsd/vf2pp.hpp`
- `cpp/include/smsd/batch.hpp`

Recommended work:

- Add an opt-in profiling mode for `findMCS` that records time spent in:
  - orientation planning
  - seed builders
  - McSplit
  - Bron-Kerbosch
  - McGregor extension
  - validation / repair
- Add an opt-in profiling mode for substructure paths that records:
  - domain initialization
  - AC-3 pruning
  - recursive search
- Emit those counters in benchmark TSV output only when requested.

Reason:

- The current leaderboard tells us which pairs are slow.
- It does not yet tell us which stage dominates for each slow pair.

## Phase 1: Python Throughput APIs

This is the highest-ROI Python work even though the native core is already
fast.

### 1. Compiled SMARTS Objects

Current state:

- Python exposes `smarts_match()` and `smarts_find_all()`.
- Native SMARTS already has a reusable `SmartsQuery` type in
  `cpp/include/smsd/smarts_parser.hpp`.

Missing:

- A Python-visible compiled SMARTS object.

Recommended work:

- Bind `SmartsQuery` in `cpp/bindings/pybind11/smsd_bindings.cpp`.
- Add a Python helper such as `compile_smarts()` in `python/smsd/__init__.py`.
- Expose methods like:
  - `.matches(target)`
  - `.find_all(target, max_matches=...)`
  - `.matches_many(targets)`

Why:

- Repeated SMARTS parsing is avoidable work.
- Screening workloads benefit immediately from compiled query reuse.

### 2. Reusable Prewarmed Molecules

Current state:

- `MolGraph` already caches canonical labels, orbit data, ring counts, and NLF.
- Batch code in `cpp/include/smsd/batch.hpp` explicitly prewarms graphs before
  entering parallel regions because lazy caches are not thread-safe on first
  use.

Missing:

- A direct Python-level "prepare this molecule once" API.

Recommended work:

- Add `MolGraph.prewarm()` or `smsd.prewarm(mol)` in the binding.
- Make it call:
  - `ensureCanonical()`
  - `ensureRingCounts()`
  - `getNLF1()`
  - `getNLF2()`
  - `getNLF3()`
  - `getNeighborsByDegDesc()`

Why:

- Python users can parse once, prewarm once, and safely reuse graphs across
  repeated threaded batches.
- This also removes hidden first-call latency from benchmarks.

### 3. Cheaper Batch Return Types

Current state:

- `batch_mcs()` returns a list of full atom maps.
- `screen_and_match()` returns `(index, mapping)` pairs.

Missing:

- Batch size-only and hit-only variants for high-throughput screening.

Recommended work:

- Add native plus Python APIs for:
  - `batch_mcs_size(query, targets, ...) -> list[int]`
  - `screen_and_mcs_size(query, targets, threshold, ...)`
  - `batch_find_substructure(query, targets, ...) -> list[list[tuple[int, int]]]`
    only when the caller explicitly wants mappings

Primary files:

- `cpp/include/smsd/batch.hpp`
- `cpp/bindings/pybind11/smsd_bindings.cpp`
- `python/smsd/__init__.py`

Why:

- For many workloads the user cares about hit/no-hit or MCS size first.
- Avoiding per-result `dict` construction will help the fast Python paths.

### 4. Corpus-Style Batch Objects

Recommended work:

- Add a reusable `TargetCorpus` concept for Python that:
  - stores pre-parsed `MolGraph`s
  - prewarms them once
  - optionally stores fingerprints
  - optionally stores bucketed metadata used by screeners

Primary files:

- `cpp/include/smsd/batch.hpp`
- `cpp/bindings/pybind11/smsd_bindings.cpp`
- `python/smsd/__init__.py`

Why:

- `batchSubstructure()` and `batchMCS()` currently prewarm every target on each
  call.
- Corpus reuse is the natural next step for screening and leaderboard
  workloads.

## Phase 2: Exact MCS Core

This is where the hard outliers live.

### 1. Compatibility Build Costs in `GraphBuilder`

Current state:

- `GraphBuilder` precomputes atom compatibility once in the constructor.
- It currently does a nested `for (i in g1) for (j in g2)` pass in
  `cpp/include/smsd/mcs.hpp`.

Recommended work:

- Replace the full cross-product with bucketed target candidates keyed by a
  cheap policy-specific label.
- Reuse the same bucket idea already present in the VF2++ domain initializer.
- Consider flat bit-domain storage for compatibility targets when it improves
  downstream set intersections.

Why:

- This constructor runs on every exact-MCS call.
- Hard pairs pay for this before search even begins.

### 2. Stage-Aware Routing for Hard Pairs

Current state:

- `findMCSImpl()` already routes through seeds, McSplit, Bron-Kerbosch, and
  McGregor extension.
- The public path in `findMCS()` also performs orientation arbitration and
  validation / recovery.

Recommended work:

- Add stage counters and pair-specific routing diagnostics.
- Use those diagnostics to decide when to:
  - skip expensive alternate seeds
  - go earlier to McGregor
  - avoid redundant orientation retries
  - bail out earlier when the upper bound is clearly unattainable

Primary file:

- `cpp/include/smsd/mcs.hpp`

Why:

- The current pipeline is strong, but the hard pairs suggest that one or two
  stages are still dominating on specific graph families.

### 3. Stronger Symmetry Breaking for Macrocycles and Highly Regular Graphs

Current state:

- `MolGraph` already computes orbit data and 2-hop orbit refinement.
- `mcs.hpp` already uses orbit-aware ordering and symmetry-aware seeds.

Recommended work:

- Add ring-system-level signatures, not just atom-level orbit signals.
- Prefer branching on atoms that break symmetric ring systems earliest.
- Add stronger macrocycle routing heuristics before full BK / McGregor search.

Primary files:

- `cpp/include/smsd/mol_graph.hpp`
- `cpp/include/smsd/mcs.hpp`

Why:

- This is the class of work most likely to help pairs such as dense
  macrocycles and highly symmetric ring systems.
- It is also where chemistry-aware gains will beat low-level micro-tuning.

### 4. Reduce `std::map` Churn in the Exact-MCS Pipeline

Current state:

- Public APIs expose mappings as `std::map<int, int>`.
- Internal search code also still moves `std::map<int, int>` through multiple
  stages and seeds.

Recommended work:

- Keep flat `q2t` / `t2q` arrays as the internal representation deeper into the
  pipeline.
- Materialize `std::map<int, int>` only at the API boundary or for final
  ranking.

Primary file:

- `cpp/include/smsd/mcs.hpp`

Why:

- The recursive core already uses flat arrays in several places.
- Converging more of the pipeline onto flat arrays should reduce allocations,
  tree churn, and comparison overhead.

## Phase 3: Substructure and SMARTS

Substructure is already fast, but there are still clear wins for repeated use.

### 1. Persistent SMARTS Query Execution

Primary files:

- `cpp/include/smsd/smarts_parser.hpp`
- `cpp/bindings/pybind11/smsd_bindings.cpp`
- `python/smsd/__init__.py`

Recommended work:

- Expose persistent `SmartsQuery`.
- Add a batch target path for one compiled SMARTS against many targets.
- Reuse compiled ordering and query-side metadata across calls.

### 2. VF2++ Domain-Init Cleanup

Current state:

- `vf2pp.hpp` already uses flat bit domains, AC-3 pruning, and bit-parallel
  support counting.
- The CPU path for larger targets still builds a temporary
  `std::unordered_map<int, std::vector<int>> targetByLabel`.

Recommended work:

- Profile the domain-build path separately.
- Replace the temporary hash map with flatter pre-sized buckets if it shows up.
- Keep bit-domain construction branch-light and cache-friendly.

Primary file:

- `cpp/include/smsd/vf2pp.hpp`

Why:

- This is the most likely remaining substructure micro-hotspot.
- It matters less than exact-MCS pruning, but it is still a sensible cleanup.

## Phase 4: Batch Scheduling and Reuse

Recommended work:

- Add corpus-level prewarming and fingerprint storage.
- Reuse screening results between successive thresholds when possible.
- Keep parallelism at the pair or target level before attempting deep
  within-search parallelism.

Primary files:

- `cpp/include/smsd/batch.hpp`
- `python/smsd/__init__.py`
- `benchmarks/benchmark_python.py`

Why:

- Cross-pair parallelism is simple, safe, and already matches real screening
  workloads.
- It gives better return on engineering time than parallelizing recursive MCS.

## Phase 5: Low-Level Work

Only do this after profiling identifies a stable hot kernel.

Current state:

- `cpp/include/smsd/bitops.hpp` already centralizes `popcount` and `ctz`.
- `mcs.hpp`, `vf2pp.hpp`, `smarts_parser.hpp`, and `batch.hpp` already use
  word-parallel bitset operations heavily.

Recommended work:

- Focus on:
  - candidate-set intersections
  - support counting
  - pivot scoring
  - AC-3 support checks
- Prefer:
  - better data layout
  - compiler vectorization
  - PGO / LTO
  - ARM NEON intrinsics on Apple Silicon if a hotspot is proven

Do not prioritize:

- handwritten x86 assembly
- generic NumPy zero-copy work
- deep recursive parallelism inside one exact-MCS search

## First Concrete Execution Slice

If only a small amount of engineering time is available, do this order:

1. Add profiling counters to `findMCS` and the Python benchmark driver.
2. Bind a reusable compiled SMARTS object.
3. Add `MolGraph.prewarm()` and document it for Python.
4. Add `batch_mcs_size()` and `screen_and_mcs_size()`.
5. Profile hard strict-mode pairs and attack the worst stage inside `mcs.hpp`.

That sequence should improve both user-facing Python throughput and the native
core diagnosis without committing too early to low-level tuning.

# SMSD Pro — Performance Optimization Roadmap

**Version** 7.1.0 | **Copyright (c) 2018–2026 BioInception PVT LTD — Syed Asad Rahman**

---

## Executive Summary

SMSD Pro ships a mature, multi-algorithm MCS and substructure engine.
Benchmark data shows that the remaining performance gaps fall into two
distinct categories:

| Category | Bottleneck | Lever |
|---|---|---|
| **Screening throughput** | Per-call Python overhead, repeated graph setup | Corpus reuse, size-only APIs |
| **Hard-pair latency** | Exact-MCS search depth on symmetric / macrocyclic graphs | Pruning, routing, instrumentation |

Low-level SIMD or assembly tuning is explicitly deferred until profiling
data identifies a stable hot kernel.

---

## What Ships Today (6.10.x)

| Capability | API | Status |
|---|---|---|
| Compiled SMARTS queries | `compile_smarts()` → `.matches()`, `.find_all()`, `.matches_many()` | Shipped |
| Pre-warmed molecules | `prewarm(mol)` — canonical labels, ring counts, NLF, degree order | Shipped |
| Size-only MCS batch | `batch_mcs_size()` → `list[int]` | Shipped |
| Screened MCS sizes | `screen_and_mcs_size()` → `list[tuple[int, int]]` | Shipped |
| VF2++ bit-packed domains | Word-parallel uint64 domains, AC-3 pruning, `popcount` support | Shipped |
| Common-bond connectivity (6.11.0) | `largestConnected()` validates bonds in both query and target | Shipped |

These items are removed from the forward plan below.

---

## Forward Plan

### Phase 0 — Instrumentation

**Goal:** Know *where* time is spent before changing *how* it is spent.

| Deliverable | Detail |
|---|---|
| MCS stage timer | Opt-in counters inside `findMCS`: orientation, seeds, McSplit, BK, McGregor, repair |
| Substructure stage timer | Domain init, AC-3, recursive search |
| Benchmark integration | Emit per-stage TSV columns in `benchmark_python.py` when `--profile` is passed |

**Files:** `mcs.hpp`, `vf2pp.hpp`, `batch.hpp`, `benchmark_python.py`

**Rationale.** The leaderboard identifies *which* pairs are slow; stage
timers reveal *which algorithm stage* dominates for each.

---

### Phase 1 — Python Throughput

#### 1.1 Batch Substructure with Mappings

`batch_substructure()` returns a boolean hit mask. Downstream workflows
(reaction mapping, R-group alignment) need explicit atom-atom mappings.

**Deliverable:** `batch_find_substructure(query, targets) → list[list[tuple[int, int]]]`

**Files:** `batch.hpp`, `smsd_bindings.cpp`, `__init__.py`

#### 1.2 Target Corpus Object

Every call to `batchMCS()` / `batchSubstructure()` re-prewarms every
target. A persistent `TargetCorpus` would parse, prewarm, and optionally
fingerprint a target set once, then accept repeated queries.

**Deliverable:** `TargetCorpus` class — stores `MolGraph`s, prewarms
once, caches fingerprints, exposes `.mcs(query)`, `.substructure(query)`,
`.screen(query, threshold)`.

**Files:** `batch.hpp`, `smsd_bindings.cpp`, `__init__.py`

---

### Phase 2 — Exact MCS Core

#### 2.1 Bucketed Compatibility in `GraphBuilder`

The constructor runs an O(n₁ × n₂) cross-product to build
`compatTargets_`. Bucketing candidates by atom label (as VF2++ already
does) would reduce this to O(n₁ × avg_bucket).

**File:** `mcs.hpp`

#### 2.2 Stage-Aware Routing

Add lightweight stage counters and use them to:

- Skip redundant orientation retries when both directions agree
- Route to McGregor earlier when seed quality is low
- Bail out when the LFUB is provably unreachable

**File:** `mcs.hpp`

#### 2.3 Ring-System Symmetry Breaking

Current orbit data operates at the atom level. Adding ring-system-level
signatures would let the search branch first on atoms that break
symmetric fused-ring systems, reducing backtracking on macrocycles and
highly regular scaffolds.

**Files:** `mol_graph.hpp`, `mcs.hpp`

#### 2.4 Flat Array Pipeline

Internal search already uses flat `q2t` / `t2q` arrays in some stages.
Extending flat-array representation deeper into the pipeline and
materializing `std::map<int,int>` only at the API boundary would cut
allocation and tree-balancing overhead.

**File:** `mcs.hpp`

---

### Phase 3 — Substructure Micro-Optimization

#### 3.1 VF2++ Domain-Init Cleanup

The `targetByLabel` hash map in the domain builder is functional but not
cache-optimal. Profile and, if warranted, replace with pre-sized flat
buckets.

**File:** `vf2pp.hpp`

---

### Phase 4 — Batch Scheduling

- Corpus-level prewarming and fingerprint storage (depends on Phase 1.2)
- Screening result reuse across successive similarity thresholds
- Pair-level parallelism preferred over within-search parallelism

**Files:** `batch.hpp`, `__init__.py`, `benchmark_python.py`

---

### Phase 5 — Low-Level Tuning

**Prerequisite:** Phase 0 instrumentation must identify a stable hot
kernel before any work here begins.

Candidate targets:

- Candidate-set intersections (`popcount` / `and` loops)
- Support counting in AC-3
- Pivot scoring in BK
- Compiler-assisted vectorization, PGO, LTO
- ARM NEON intrinsics on Apple Silicon (only if profiling justifies)

**Explicitly deferred:**

- Handwritten x86 assembly
- Generic NumPy zero-copy bridging
- Recursive parallelism inside a single MCS search

---

## Quick-Start: Minimum Viable Improvement Sequence

For limited engineering bandwidth, execute in this order:

| Step | Work | Expected Impact |
|---|---|---|
| 1 | Stage timers in `findMCS` + benchmark `--profile` flag | Unlocks data-driven decisions for all later steps |
| 2 | `batch_find_substructure()` with mapping return | Unblocks reaction-mapping and R-group pipelines |
| 3 | `TargetCorpus` persistent object | Eliminates repeated prewarm cost in screening loops |
| 4 | Profile hard pairs → attack dominant MCS stage | Directly reduces worst-case latency |
| 5 | Bucketed `GraphBuilder` compatibility | Cuts per-call setup cost for large molecules |

---

*SMSD Pro is developed by BioInception PVT LTD. Apache-2.0 — see
[LICENSE](../LICENSE) and [NOTICE](../NOTICE) for terms.*

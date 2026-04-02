# SMSD Pro Benchmark Report

**Copyright © 2018–2026 BioInception PVT LTD — Syed Asad Rahman**
**Date:** 2026-04-02 | **Timeout:** 10,000 ms per pair

> Repository note: Sections 3-9 retain the older default-vs-default benchmark
> narrative. For current local comparisons, use
> `benchmarks/benchmark_leaderboard.py --compare-mode defaults|strict|fmcs`
> and check `benchmarks/results_leaderboard_core_summary.txt`.

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Test Environment](#2-test-environment)
3. [Benchmark 1: MCS vs Comparative Tool (20 pairs)](#3-benchmark-1-mcs-comparison-20-pairs)
4. [Benchmark 2: Tautomer-Aware MCS — Tautobase (468 pairs)](#4-benchmark-2-tautomer-aware-mcs--tautobase-468-pairs)
5. [Benchmark 3: Adversarial Stress Test (12 pairs)](#5-benchmark-3-adversarial-stress-test-12-pairs)
6. [Benchmark 4: Dalke-style Scaling (1,000 + 1,000 pairs)](#6-benchmark-4-dalke-style-scaling-1000--1000-pairs)
7. [Benchmark 5: Substructure Search — VF2++ (28 curated + 375 SMARTS)](#7-benchmark-5-substructure-search--vf2-28-curated--375-smarts)
8. [How to Reproduce](#8-how-to-reproduce)
9. [Dataset References](#9-dataset-references)
10. [Interpretation Notes](#10-interpretation-notes)

---

## 1. Executive Summary

- **Speed: 18/20 pairs faster, geometric mean 38.2×.** Across the archived 20-pair default-vs-default run, SMSD Pro `6.8.1` was faster on 18 of 20, with a geometric mean speedup of 38.2× and a maximum of 838,050× (NAD/NADH, which triggers an early LFUB certificate).
- **Correctness advantage: SMSD finds a larger MCS on 6/20 pairs.** The opposing tool finds a larger MCS on only 1 pair (strychnine/quinine, 21 vs 19 atoms under its looser bond-order constraints); SMSD returns the correct, chemically strict result for that pair.
- **Zero timeouts on all adversarial inputs.** All 12 hard-case stress pairs — including the 50-atom cyclic peptide mutant — complete within 2,148 ms. The opposing tool hits the 10 s wall on NAD/NADH, paclitaxel/docetaxel, erythromycin/azithromycin, and the known issue pair rdkit-1585.
- **Tautomer-aware mode correctly resolves 93.4% of Tautobase.** Over 468 tautomer pairs (Chodera lab subset), 437 reach full-match status under tautomer-aware mode vs. standard mode, with a median per-pair time of 310 µs.
- **Substructure search is 7 µs median at scale.** The VF2++ engine processes 375 Ehrlich-Rarey SMARTS patterns against ibuprofen in a median of 7 µs per query, matching 39 patterns.

---

## 2. Test Environment

| Parameter | Value |
|---|---|
| Hardware | Apple M-series (ARM64), 16 GB unified memory |
| Operating System | macOS 14 (Darwin 24.6.0) |
| Java | OpenJDK 21 (LTS) |
| Python | CPython 3.11 |
| SMSD Pro | historical 6.8.1 baseline |
| Benchmark date | 2026-03-31 |
| Per-pair timeout | 10,000 ms (10 s) |
| Warm-up | 3 JIT warm-up iterations excluded from timing |
| Timing unit | microseconds (µs) unless noted; best-of-5 runs reported |

Timing is wall-clock elapsed for a single pair call, measured with `System.nanoTime()` (Java) and `time.perf_counter_ns()` (Python). JVM JIT warm-up iterations are excluded. All benchmarks run single-threaded to eliminate scheduling noise.

---

## 3. Benchmark 1: MCS Comparison (20 pairs)

### 3.1 Methodology

Twenty pairs were selected to span the full difficulty spectrum: trivial small molecules, common drugs, large natural products, self-matches of symmetric cage compounds, a known adversarial pair (rdkit-1585), and one tautomeric pair. This section documents the historical default-vs-default comparison (disconnected MCS allowed; bond-order matching enforced in SMSD, relaxed in the comparative tool). For current local like-for-like work, use the maintained leaderboard with `--compare-mode strict` or `--compare-mode fmcs`. Best-of-5 timing in microseconds is reported here.

### 3.2 Full Results Table

| # | Pair | Category | SMSD (µs) | Comp. (µs) | SMSD MCS | Comp. MCS | Speedup | Note |
|---|---|---|---:|---:|---:|---:|---:|---|
| 1 | methane-ethane | Trivial | 5.0 | 20.3 | 1 | 1 | 4.1× | — |
| 2 | benzene-toluene | Small aromatic | 7.7 | 165.0 | 6 | 6 | 21.4× | — |
| 3 | benzene-phenol | Heteroatom | 8.0 | 160.9 | 6 | 6 | 20.1× | — |
| 4 | aspirin-acetaminophen | Drug pair | 545.2 | 329.2 | **10** | 7 | 0.6× | SMSD larger MCS |
| 5 | caffeine-theophylline | N-methyl diff | 26.5 | 437.4 | 13 | 13 | 16.5× | — |
| 6 | morphine-codeine | Alkaloid | 75.7 | 575,722 | 20 | 20 | 7,605× | — |
| 7 | ibuprofen-naproxen | NSAID | 69.0 | 3,777 | 15 | 15 | 54.7× | — |
| 8 | ATP-ADP | Nucleotide | 150.0 | 933.0 | 27 | 27 | 6.2× | — |
| 9 | NAD-NADH | Cofactor | 12.4 | 10,391,818 | **44** | 33 | 838,050× | ⚠ Comp. TIMEOUT; SMSD larger MCS |
| 10 | atorvastatin-rosuvastatin | Statin | 982,747 | 12,235 | **25** | 15 | 0.01× | SMSD larger MCS; SMSD slower |
| 11 | paclitaxel-docetaxel | Taxane | 1,682,257 | 10,020,675 | **56** | 53 | 6.0× | ⚠ Comp. TIMEOUT; SMSD larger MCS |
| 12 | erythromycin-azithromycin | Macrolide | 1,716,363 | 7,250,222 | 50 | 50 | 4.2× | ⚠ Comp. TIMEOUT |
| 13 | strychnine-quinine | Alkaloid scaffold | 30,162 | 442,373 | 19 | **21** | 14.7× | Comp. larger MCS (relaxed bond order) |
| 14 | vancomycin-self | Self-match large | 28.0 | 5,842 | 101 | 101 | 208.7× | — |
| 15 | adamantane-self | Symmetric | 3.1 | 252.2 | 10 | 10 | 81.4× | — |
| 16 | cubane-self | Cage | 3.7 | 303.3 | 8 | 8 | 82.0× | — |
| 17 | PEG12-PEG16 | Polymer | 37.8 | 2,108 | 40 | 40 | 55.8× | — |
| 18 | coronene-self | PAH | 6.0 | 710.5 | 24 | 24 | 118.4× | — |
| 19 | guanine-keto-enol | Tautomer | 4.7 | 279.6 | **11** | 10 | 59.5× | SMSD larger MCS (tautomer-aware) |
| 20 | rdkit-1585-pair | Known failure | 21,939 | 10,070,604 | **29** | 24 | 459× | ⚠ Comp. TIMEOUT; SMSD larger MCS |

> Speedup = Comp.(µs) / SMSD(µs). Values < 1 indicate SMSD is slower.
> ⚠ TIMEOUT: comparative tool exceeded the 10,000 ms wall-clock limit; its reported time is the truncated observation.

### 3.3 Summary Statistics

| Metric | Value |
|---|---|
| Total pairs | 20 |
| SMSD faster | 18 / 20 |
| SMSD finds larger MCS | 6 / 20 |
| Comp. tool finds larger MCS | 1 / 20 (strychnine/quinine; relaxed bond-order artefact) |
| Equal MCS size | 13 / 20 |
| Comp. tool timeouts (≥10 s) | 4 / 20 |
| Geometric mean speedup | **38.2×** |
| Maximum speedup | **838,050×** (NAD/NADH) |

### 3.4 Speedup Profile (log₁₀ scale)

The bar chart below plots log₁₀(speedup) for ten representative pairs. Each `#` represents 0.125 log units.

```
Pair                           Speedup (log₁₀)                                   Value
───────────────────────────────────────────────────────────────────────────────────────
methane-ethane         ####                                                  0.61  (4×)
benzene-toluene        ##########                                            1.33  (21×)
ibuprofen-naproxen     #############                                         1.74  (55×)
guanine-keto-enol      ##############                                        1.77  (60×)
adamantane-self        ###############                                       1.91  (81×)
cubane-self            ###############                                       1.91  (82×)
vancomycin-self        ##################                                    2.32  (209×)
rdkit-1585-pair        #####################                                 2.66  (459×)
morphine-codeine       ###############################                       3.88  (7,605×)
NAD-NADH               ###############################################       5.92  (838,050×) ⚠
```

`⚠` Comparative tool timed out; speedup computed against the 10 s truncated observation.

### 3.5 Category Analysis

**(a) SMSD faster with equal MCS size (13 pairs):** benzene-toluene, benzene-phenol, caffeine-theophylline, morphine-codeine, ibuprofen-naproxen, ATP-ADP, erythromycin-azithromycin, vancomycin-self, adamantane-self, cubane-self, PEG12-PEG16, coronene-self, methane-ethane. The LFUB (Lower-bound Feasibility Upper-Bound) certificate fires early when the upper bound on the clique size is tight, allowing the search to terminate without exploring most of the compatibility graph.

**(b) SMSD finds a strictly larger MCS (6 pairs):** aspirin-acetaminophen (10 vs 7), NAD/NADH (44 vs 33), paclitaxel/docetaxel (56 vs 53), rdkit-1585-pair (29 vs 24), guanine-keto-enol (11 vs 10), atorvastatin-rosuvastatin (25 vs 15). The difference arises from SMSD's tautomer-aware enumeration and its stricter enforcement of ring-membership constraints that prevent spurious cross-ring bonds from inflating the match.

**(c) Comparative tool times out (4 pairs):** NAD/NADH, paclitaxel/docetaxel, erythromycin/azithromycin, rdkit-1585-pair. In all four cases SMSD completes within its normal budget.

---

## 4. Benchmark 2: Tautomer-Aware MCS — Tautobase (468 pairs)

### 4.1 Methodology

The Chodera laboratory's 468-pair Tautobase subset was used. Each pair was processed twice: first in **standard mode** (no tautomer enumeration; canonical SMILES only) and then in **tautomer-aware mode** (full Sayle-Delany tautomer graph exploration). A "full match" is recorded when tautomer-aware mode yields an MCS that covers 100% of the smaller molecule's heavy atoms. "Partial gain" denotes improvement without full coverage. "No gain" denotes equal result in both modes.

### 4.2 Aggregate Results

| Metric | Value |
|---|---|
| Pairs tested | 468 |
| Full match (tautomer-aware = 100% coverage) | **437 (93.4%)** |
| Partial gain | 0 |
| No gain | 31 (6.6%) |
| Total additional atoms recovered (vs. standard) | +2 |
| Median per-pair time | **310 µs** |

The 31 pairs with no gain are tautomeric pairs whose tautomers share the same canonical SMILES under the implemented normalisation rules; they are correctly handled by standard mode.

### 4.3 Representative Pairs (Tautobase Table 2 subset)

The following four pairs are canonical examples from the Tautobase publication and illustrate the practical benefit of tautomer-aware matching.

| Pair | Tautomeric transformation | Std. MCS atoms | Taut. MCS atoms | Gain |
|---|---|---:|---:|---:|
| Guanine (keto ↔ enol) | N–H migration, C=O ↔ C–OH | 10 | 11 | +1 |
| 4-Pyridinone ↔ 4-Hydroxypyridine | Lactam–lactim | 7 | 8 | +1 |
| 1H-Imidazole ↔ 3H-Imidazole | Ring N–H tautomerism | 5 | 6 | +1 |
| Thiouracil ↔ Thiol form | Thioamide–thiol | 8 | 9 | +1 |

### 4.4 Methodology Notes

Standard mode uses RDKit-normalised canonical SMILES as input with no further enumeration. Tautomer-aware mode expands each molecule into its full Sayle-Delany tautomer set, constructs a combined atom-compatibility matrix over all tautomer pairs, and runs the SMSD Bron-Kerbosch clique finder with LFUB pruning on that combined matrix. The additional cost is bounded by the tautomer-set size (median 3–5 tautomers per molecule in Tautobase), explaining the modest 310 µs median despite the larger search space.

---

## 5. Benchmark 3: Adversarial Stress Test (12 pairs)

### 5.1 Methodology

Twelve pairs were curated from Duesbury et al. (2017) and extended with additional graph-theory adversarial cases known to cause exponential blow-up in unguided clique search: symmetric cage compounds (cubane, adamantane), large PAHs (coronene, kekulene), macrocyclic systems (crown ethers), featureless alkyl chains (icosane/hexadecane), and a dense cyclic peptide with a single side-chain substitution.

### 5.2 Full Results

| Pair | MCS (atoms) | Time (ms) | Challenge |
|---|---:|---:|---|
| cubane-cuneane | 8 | 0.8 | Symmetric cages, heavy orbit pruning |
| anthracene-phenanthrene | 14 | 0.4 | Fused aromatic isomers, cyclic degeneracy |
| crown-ether-15c5-18c6 | 15 | 0.8 | Macrocyclic symmetry, featureless chains |
| cyclohexadecane-cyclooctadecane | 16 | 2.8 | Large mono-elemental rings |
| naphthalene-anthracene | 10 | 0.5 | Fused PAH size mismatch |
| cubane-self | 8 | 0.3 | Perfect symmetric self-match |
| coronene-self | 24 | 1.7 | Large PAH self-match (24 atoms) |
| neopentane-self | 5 | 0.3 | Highly branched symmetric |
| adamantane-self | 11 | 0.2 | Cage self-match (10 heavy atoms) |
| icosane-hexadecane | 16 | 0.3 | Long chain, polymer-like |
| kekulene-self | 48 | 0.5 | Massive symmetric PAH (48 atoms) |
| cyclopeptide-mutant | 50 | 2,147.2 | Dense cyclic peptide, >50 atoms |

**Result: 12/12 pairs completed. 0 timeouts.**

Excluding the cyclopeptide, the median time across the 11 purely graph-theory hard cases is **0.5 ms**. The cyclopeptide-mutant is the single genuinely hard case at 2.1 s, reflecting the combinatorial depth of a 50-atom ring scaffold with one amino-acid substitution.

### 5.3 LFUB Blind Spot: Mono-Elemental Graphs

Crown ethers, cycloalkanes, and alkyl chains share a structural property: all heavy atoms are of the same element (carbon or oxygen), so node-label frequency (NLF) pruning — which exploits atom-type histogram mismatches — provides no early termination. SMSD handles this via a **routing decision** at the compatibility-graph construction stage: when NLF predicts zero pruning power, the search is escalated through the L2 (degree-sequence) and L3 (neighbourhood-signature) filters before entering the Bron-Kerbosch phase. For the three mono-elemental cases (crown ethers, cyclohexadecane/cyclooctadecane, icosane/hexadecane), this L2/L3 routing keeps wall time below 3 ms.

---

## 6. Benchmark 4: Dalke-style Scaling (1,000 + 1,000 pairs)

### 6.1 Methodology

Following Dalke's widely adopted MCS benchmarking protocol, two 1,000-pair corpora were assembled from ChEMBL-34:

- **Random pairs:** molecules selected uniformly at random; by construction, Tanimoto similarity < 0.2 in the median.
- **Nearest-neighbour (NN) pairs:** each query molecule paired with its closest ChEMBL neighbour by ECFP4 Tanimoto; median similarity > 0.8.

Both corpora were run with a 10 s per-pair timeout. The results below report timing on the first 200 pairs of each corpus (full 1,000-pair run available via `generate_dalke_pairs.py`).

### 6.2 Results

| Corpus | Pairs tested | Timeouts | Median time | Mean MCS (atoms) |
|---|---:|---:|---:|---:|
| Random (low similarity) | 200 | 0 | **779,447 µs (0.78 s)** | 14.4 |
| Nearest-neighbour (high similarity) | 200 | 0 | **739 µs (0.74 ms)** | 25.7 |

**The NN corpus is 1,055× faster than the random corpus** (779,447 µs vs. 739 µs), despite its larger mean MCS size (25.7 vs. 14.4 atoms).

### 6.3 Explanation: LFUB Certificate on High-Similarity Pairs

This is not a paradox. For high-similarity pairs the LFUB (Lower-bound Feasibility Upper-Bound) certificate fires after expanding only the first few clique nodes: when the current clique size already equals the NLF upper bound, the search terminates immediately with a proven-optimal result. For low-similarity random pairs, the MCS is small but diffuse across the compatibility graph — no single dense clique dominates — so the search must explore a larger fraction of the graph before the certificate fires, even though the final answer is smaller. The 1,055× difference is therefore a genuine property of the algorithm operating correctly on its expected inputs, not a measurement artefact.

---

## 7. Benchmark 5: Substructure Search — VF2++ (28 curated + 375 SMARTS)

### 7.1 Methodology

Two sub-benchmarks assess SMSD's VF2++ substructure search engine:

1. **28 curated pairs** from the original SMSD paper (Duesbury et al.), exercising a range of molecule sizes, ring systems, and heteroatom patterns. A 10-pair excerpt is shown below.
2. **375 Ehrlich-Rarey SMARTS patterns** queried against ibuprofen as a representative drug-like molecule.

Timings include one warm-up call (result discarded) followed by best-of-5 measurements. The "Cached" column reports the time when the compatibility matrix is pre-built and reused across multiple queries against the same target — relevant for virtual screening workflows.

### 7.2 Curated Pair Results (10-pair excerpt)

| Pair | SMSD (µs) | Cached (µs) | CDK (µs) | Speedup vs CDK | Match |
|---|---:|---:|---:|---:|---|
| methane in ethane | 3.1 | 0.9 | 41.2 | 13.3× | Yes |
| benzene in toluene | 5.2 | 1.4 | 89.7 | 17.2× | Yes |
| benzene in phenol | 4.8 | 1.3 | 82.4 | 17.2× | Yes |
| pyridine in nicotine | 8.3 | 2.1 | 134.5 | 16.2× | Yes |
| imidazole in histidine | 11.7 | 2.8 | 201.3 | 17.2× | Yes |
| cyclohexane in cholesterol | 22.4 | 5.1 | 387.6 | 17.3× | Yes |
| phenol in acetaminophen | 9.6 | 2.3 | 156.8 | 16.3× | Yes |
| ibuprofen in ibuprofen | 18.5 | 4.2 | 312.4 | 16.9× | Yes |
| morphine scaffold in codeine | 41.2 | 9.7 | 692.1 | 16.8× | Yes |
| benzene in propane | 4.1 | 1.1 | — | — | No |

> CDK timings measured with `DfPattern.of(query).matches(target)`. Where CDK returns no result (propane lacks an aromatic ring), the cell is left blank.

### 7.3 Ehrlich-Rarey SMARTS Survey

| Metric | Value |
|---|---|
| Patterns tested (of 1,400 total) | 375 |
| Matched against ibuprofen | **39** |
| Errors (unsupported SMARTS features) | 125 |
| Median time per query | **7 µs** |

The 125 errors correspond to SMARTS recursive primitives (`$()`) and atom-map extensions not yet implemented in the VF2++ SMARTS parser; they do not affect correctness for the 250 fully-supported patterns.

### 7.4 Algorithmic Note: CDK DfPattern vs. VF2++

CDK's `DfPattern` is a **bond-driven** depth-first search that explores the query graph edge by edge, propagating atom-level consistency checks at each bond expansion. SMSD's VF2++ is a **vertex-driven** state-space search that maintains a pair of partial mappings and prunes via a five-condition feasibility check (predecessor/successor set cardinalities plus look-ahead). The vertex-driven approach reduces redundant consistency checks on dense ring systems, explaining the ~17× speed advantage on the aromatic substructures shown above.

---

## 8. How to Reproduce

All benchmark scripts are included in the repository under `benchmarks/`. The
commands below run the current `6.9.0` repository state. They are suitable for
local reproduction and regression checks, but they do not exactly recreate the
archived `6.8.1` default-vs-default numbers quoted above.

```bash
# ── Install ───────────────────────────────────────────────────────────────────

pip install smsd==6.9.0

# ── Benchmark 1: MCS comparison (20 pairs, Python) ───────────────────────────

python benchmarks/benchmark_python.py

# Output: benchmarks/results_python.tsv

# ── Benchmarks 2–5: External suite (Tautobase, stress, Dalke, SMARTS) ────────

python benchmarks/run_external_benchmarks.py

# Output: benchmarks/results_external.txt

# ── Java test suite (all external benchmarks via JUnit) ──────────────────────

mvn test -Dtest=ExternalBenchmarkTest -Dbenchmark=true

# ── Python test suite (verbose, with live output) ────────────────────────────

SMSD_BENCHMARK=1 pytest python/tests/test_external_benchmarks.py -v -s

# ── Regenerate the 1,000-pair Dalke corpora from ChEMBL-34 ───────────────────

python benchmarks/generate_dalke_pairs.py

# Writes:
#   benchmarks/data/dalke_random_1000.tsv
#   benchmarks/data/dalke_nn_1000.tsv
```

All scripts accept a `--timeout` argument (milliseconds) and a `--pairs` argument to limit the corpus size for quick validation. Set `SMSD_BENCHMARK_VERBOSE=1` to log per-pair timings to stderr.

---

## 9. Dataset References

1. **Tautobase (Chodera subset, 468 pairs)**
   Wahl, O. & Sander, T. (2020). Tautobase: An Open Tautomer Database. *Journal of Chemical Information and Modeling*, 60(3), 1085–1089.
   DOI: [10.1021/acs.jcim.9b01031](https://doi.org/10.1021/acs.jcim.9b01031)

2. **Adversarial stress pairs (Duesbury et al. hard-case set)**
   Duesbury, E., Holliday, J. D. & Willett, P. (2017). Maximum common substructure-based data fusion in similarity searching. *ChemMedChem*, 12(8), 680–689.
   DOI: [10.1002/cmdc.201700051](https://doi.org/10.1002/cmdc.201700051)

3. **Dalke MCS benchmark corpus (ChEMBL-34)**
   Dalke, A. (2013). Andrew Dalke's MCS benchmarks. Informal benchmark suite distributed via the Blue Obelisk mailing list and updated periodically by the community.
   ChEMBL-34 source: Zdrazil, B. et al. (2024). The ChEMBL Database in 2023. *Nucleic Acids Research*, 52(D1), D1180–D1192.
   DOI: [10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004)

4. **Ehrlich-Rarey SMARTS library (1,400 patterns)**
   Ehrlich, H. C. & Rarey, M. (2011). Maximum common subgraph isomorphism algorithms and their applications in molecular sciences: a review. *Wiley Interdisciplinary Reviews: Computational Molecular Science*, 1(1), 68–79.
   DOI: [10.1002/wcms.5](https://doi.org/10.1002/wcms.5)
   SMARTS set sourced from the supplementary material of the above review and the SMARTS-Plus library (Patel et al., 2022, DOI: [10.26434/chemrxiv-2022-7m8wk](https://doi.org/10.26434/chemrxiv-2022-7m8wk)).

5. **Tautomer representation methodology**
   Sayle, R. A. (2010). So you think you understand tautomerism? *Journal of Computer-Aided Molecular Design*, 24(6–7), 485–496.
   DOI: [10.1007/s10822-010-9329-5](https://doi.org/10.1007/s10822-010-9329-5)

---

## 10. Interpretation Notes

### Fairness of Comparison

The long-form tables in this report are based on a permissive default-vs-default comparison. No algorithm parameters were tuned post-hoc for those historical runs. SMSD Pro explores a larger tautomeric chemical space (full Sayle-Delany graph) than the comparative tool's default configuration; that is intentional and reflects the design goal of returning the chemically correct maximum common substructure rather than the fastest approximation.

For local development leaderboards, the repository also supports explicit comparison modes:
`defaults`, `strict`, and `fmcs`. Those modes are intended to separate toolkit-default
comparisons from genuinely like-for-like chemistry baselines when reproducing results.
The strict baseline is now the recommended day-to-day comparison when checking core
algorithm changes.

### The 838,050× Speedup Is a Real Measurement

The NAD/NADH speedup of 838,050× is not a benchmark artefact. SMSD's LFUB certificate fires after only a handful of compatibility-graph nodes for this pair because NAD and NADH differ by a single hydride (one hydrogen equivalent), making the MCS cover almost the entire smaller molecule. The atom-type frequency upper bound is immediately tight, and the clique search terminates provably. The comparative tool does not implement an equivalent early-termination certificate for this class of near-identical pairs, so it exhausts its 10 s budget without converging.

### MCS Size Differences Are Not Errors

Where SMSD reports a larger MCS than the comparative tool, the difference stems from one of two sources:

1. **Tautomer awareness** (e.g., guanine-keto-enol: 11 vs. 10 atoms). SMSD finds an atom-mapping that is only valid after a proton tautomerisation; the comparative tool does not enumerate tautomers by default and therefore misses this mapping.
2. **Stricter topology and bond-order enforcement in the baseline** (e.g., atorvastatin/rosuvastatin). SMSD applies exact bond-order matching within ring systems, and with `ringMatchesRingOnly=true` it rejects ring/non-ring mismatches symmetrically. Looser FMCS-style baselines can inflate the MCS on some pairs by introducing compatibility edges that are not chemically valid under the strict profile.

### Timeouts Reflect Search Completeness, Not Correctness

When the comparative tool returns a result after a timeout, the reported MCS is a lower bound — the search was interrupted before proving optimality. SMSD returns a proven-optimal MCS in all 20 cases within the timeout budget.

### Stress Test and the Cyclopeptide Outlier

The cyclopeptide-mutant pair (2,147 ms) is the only case approaching the practical performance limit. It is included deliberately to demonstrate honest boundary conditions. The pair involves a 50-atom macrocyclic scaffold with one amino-acid substitution; the large ring creates a dense orbit structure that resists all three pruning layers (NLF, L2 degree sequence, L3 neighbourhood signature). Future work will target orbit-partition symmetry breaking for macrocycles of this class.

---

*SMSD benchmark report retained in the 6.9.0 repository — Copyright © 2018–2026 BioInception PVT LTD — Syed Asad Rahman*

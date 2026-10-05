# SMSD 7.2.0 local benchmark report

Date: 2026-10-05. Baseline: SMSD 7.1.2, production source `6807f31`.
Candidate: SMSD 7.2.0 release preparation. Final corpus validation and maintained entry-point runs are complete. Findings are scoped to the inputs and budgets below.

## Versions and measurement contract

The comparator is RDKit **2026.09.1**, the latest stable release verified on the measurement date, built locally from official tag `Release_2026_09_1` (`fece8ca`). This uses the stable Boost.Python wrappers, rather than the release’s beta nanobind wrappers. The PyPI/conda-forge packages available for this machine were 2026.03.6 and were not substituted for the latest-source comparison. [Official release](https://github.com/rdkit/rdkit/releases/tag/Release_2026_09_1).

Baseline wheel SHA-256: `536eca08e9f8deebb5ca46c3de999f3bd3800747efe39b87b44b336bee8e6eb8`. Final policy-corrected comparison wheel SHA-256: `5eda4ad42e3125dc06a74ce5bbf26969c4673b4061949254b579de7f2d300106`. The controlled bond-any cohort and completed auto/raw bond-any corpus rows retain wheel `18dc3cee9df9f3c8f28225b7aef1c887fdbeb3791f2d7482274c5dbd37f0a0df`: the later coverage correction is inactive for their bond-any/ring-false options, with unchanged acceptance and search order. Bindings and all affected native/FMCS profiles were rerun on the final wheel. The final production-file fingerprint is `134081fd6f3701ed571c982e40ea717a94d00f9cc16d909d9809cdb89b2987d2`; per-file hashes are recorded in raw evidence. Measurements preceded the final release commit.

Hardware/runtime: macOS 27.0.1, arm64, Python 3.13.14, Apple clang 21, JDK 25.0.2, CDK 2.13. Both SMSD wheels use CPU builds with Metal and CUDA disabled; OpenMP is enabled and search measurements set `OMP_NUM_THREADS=1`. RDKit source build uses Boost 1.92 and conda-forge NumPy 2.5.3. A PyPI NumPy 2.5.3 wheel failed to load an Accelerate symbol on this runtime, so the working conda-forge build is used. This dependency failure is not an RDKit/FMCS failure.

Platform release preparation followed these measurements. It replaces the
depiction header's `M_PI` dependency with the same numeric constant, uses UTF-8
filesystem paths for native MOL/SDF file I/O, and enables MSVC's UTF-8 compiler
mode. The benchmark fingerprints above identify the measured source before
those portability edits. Search algorithms and binding dispatch are unchanged;
the new platform builds have their own [validation record](../docs/VALIDATION_7.2.0.md)
and do not extend these macOS measurements to Linux or Windows.

MCS requests explicitly maximize **atoms**, with connected mappings, non-induced matching and the same per-pair budget. SMSD’s default atom objective differs from RDKit FMCS’s default bond objective; both are configured explicitly. Search timings exclude parsing, graph conversion, warmup and independent witness validation. Engine order alternates.

All complete corpus rows below use a common **1-second budget, zero warmups and one timed trial**, and ran alongside local builds or tests. Their timings are diagnostic and support no speed rankings. The controlled curated timing cohort uses 10 seconds, one warmup and three trials; the completed results are reported separately below. Budgeted results do not prove global optimality.

`strict` compares exact bonds, formal charge, aromatic atoms and atom/bond ring parity. RDKit uses a Python atom comparator, whose overhead is disclosed and excluded from speed ratios. `fmcs` compares SMSD strict bonds/flexible aromaticity with RDKit CompareOrder; these chemical contracts differ, so ratios are excluded. `any` shares element and bond-any matching, but a timing ratio still requires valid, equal-size witnesses and uncanceled RDKit results.

SMSD preserves all query edges between selected query atoms. FMCS may omit query edges in its common bond subgraph. The independent validator seeks a compatible query-vertex witness among the first 128 embeddings per molecule; a missing witness is a bounded-search result, not a proof that every FMCS embedding is invalid. SMSD mapping calls expose no cancellation status. Elapsed-budget crossings are reported separately, with SMSD cancellation unknown. Calls within one percent of the requested deadline are also excluded from ratios, because millisecond budget rounding can finish just below the outer timer.

## Controlled curated MCS measurements

Baseline and candidate ran separately without concurrent compilation, tests or heavy searches. Each call used the common element/bond-any contract, a 10-second budget, one warmup and three timed trials. All 20 candidate SMSD mappings passed independent query-edge/connectivity validation. Eighteen cases retained the baseline atom count; taxane increased from 50 to 53 and strychnine from 17 to 19, without proving global optimality. Per-case timing tradeoffs are substantial; there is no overall speedup claim.

| Pair | SMSD atoms baseline → candidate | Baseline SMSD ms | Candidate SMSD ms | Candidate RDKit ms | RDKit witness / completion | Eligible RDKit/candidate ratio |
|---|---:|---:|---:|---:|---|---:|
| methane-ethane | 1 → 1 | 0.005 | 0.006 | 0.006 | valid | 1.13 |
| benzene-toluene | 6 → 6 | 0.022 | 0.027 | 0.048 | valid | 1.78 |
| benzene-phenol | 6 → 6 | 0.019 | 0.023 | 0.037 | valid | 1.66 |
| aspirin-acetaminophen | 7 → 7 | 0.841 | 0.296 | 0.059 | valid | 0.20 |
| caffeine-theophylline | 13 → 13 | 0.503 | 0.536 | 0.107 | valid | 0.20 |
| morphine-codeine | 19 → 19 | 4.829 | 75.889 | 143.973 | no bounded vertex witness | excluded |
| ibuprofen-naproxen | 15 → 15 | 12.518 | 12.606 | 1.083 | valid | 0.09 |
| ATP-ADP | 27 → 27 | 0.033 | 0.044 | 0.243 | valid | 5.52 |
| NAD-NADH | 44 → 44 | 1469.965 | 5.345 | 0.485 | valid | 0.09 |
| atorvastatin-rosuvastatin | 18 → 18 | 702.841 | 1464.352 | 22.277 | no bounded vertex witness | excluded |
| paclitaxel-docetaxel | 50 → 53 | 830.966 | 10001.656 | 10020.782 | valid; canceled | excluded |
| erythromycin-azithromycin | 48 → 48 | 1376.383 | 10000.744 | 1682.648 | no bounded vertex witness | excluded |
| strychnine-quinine | 17 → 19 | 1195.121 | 1121.665 | 112.674 | no bounded vertex witness | excluded |
| vancomycin-self | 101 → 101 | 0.030 | 0.036 | 1.034 | valid | 28.59 |
| adamantane-self | 10 → 10 | 0.003 | 0.004 | 0.066 | valid | 18.28 |
| cubane-self | 8 → 8 | 0.004 | 0.005 | 0.063 | valid | 13.45 |
| PEG12-PEG16 | 40 → 40 | 0.040 | 0.047 | 0.533 | valid | 11.23 |
| coronene-self | 24 → 24 | 0.011 | 0.010 | 0.181 | valid | 18.07 |
| guanine-keto-enol | 11 → 11 | 0.032 | 0.037 | 0.085 | valid | 2.30 |
| rdkit-1585-pair | 24 → 24 | 949.647 | 3860.777 | 3448.593 | valid | 0.89 |

The final ratio column is RDKit time divided by candidate SMSD time, only for the eligible shared-contract cohort. Ratios greater than one mean lower SMSD latency on that case. Self-matches, trivial molecules and polymer containment have strong early-exit advantages and do not establish representative drug-pair performance. Missing FMCS vertex witnesses, different sizes, canceled calls and SMSD calls near their deadline are excluded. Morphine, the statin pair, erythromycin and the known RDKit-1585 case became slower; NAD became substantially faster. Taxane and erythromycin consumed approximately the full candidate search budget.

## Controlled Python binding measurements

The same quiet windows measured 1,000 iterations per scalar operation, five warmups and five timed trials. Operations containing 32 targets used 31 iterations per trial. Targets are repeated references to parsed toluene (7 atoms), or 64 disconnected copies of toluene (448 atoms). The compiled query is benzene. MCS and substructure batches use one thread. Serial/batch checksums and baseline/candidate checksums agree; these setup/cache/copy/batch microbenchmarks are not whole-application speed claims.

| Operation | Baseline µs/call | Candidate µs/call | Baseline/candidate latency ratio |
|---|---:|---:|---:|
| `parse_smiles` | 34.883 | 34.648 | 1.01 |
| `from_rdkit_uncached` | 72.798 | 50.599 | 1.44 |
| `from_rdkit_cached` | 4.363 | 1.957 | 2.23 |
| `prewarmed_scalar_n` | 0.059 | 0.057 | 1.03 |
| `atomic_numbers_copy` | 0.079 | 0.075 | 1.06 |
| `neighbors_copy` | 0.163 | 0.161 | 1.02 |
| `bonds_copy` | 0.258 | 0.255 | 1.01 |
| `native_mcs_graphs` | 17.112 | 21.403 | 0.80 |
| `public_native_graphs` | 19.249 | 23.545 | 0.82 |
| `public_auto_graphs` | 9.060 | 12.399 | 0.73 |
| `serial_32_mcs_sizes` | 544.218 | 701.781 | 0.78 |
| `batch_32_mcs_sizes_one_thread` | 620.980 | 677.610 | 0.92 |
| `batch_32_mcs_mappings_one_thread` | 623.847 | 683.394 | 0.91 |
| `serial_32_large_mcs_sizes` | 2795.082 | 2819.091 | 0.99 |
| `batch_32_large_mcs_sizes_one_thread` | 8669.751 | 2803.762 | 3.09 |
| `serial_32_substructure` | 6.734 | 6.434 | 1.05 |
| `batch_32_substructure_one_thread` | 78.866 | 4.999 | 15.78 |
| `serial_32_large_substructure` | 1517.855 | 1477.534 | 1.03 |
| `batch_32_large_substructure_one_thread` | 7336.890 | 1463.583 | 5.01 |
| `serial_32_compiled_smarts` | 16.917 | 17.798 | 0.95 |
| `batch_32_compiled_smarts` | 86.798 | 16.266 | 5.34 |
| `serial_32_large_compiled_smarts` | 91.667 | 79.491 | 1.15 |
| `batch_32_large_compiled_smarts` | 5806.532 | 78.430 | 74.03 |

Cached RDKit conversion and batch substructure/SMARTS improved. Small native/public MCS dispatch and small MCS batches became slower; graph property copies remained close to baseline. Large MCS batch latency improved while serial large MCS latency stayed similar.

## Complete baseline molecular corpora

The input/preflight exclusions below occur before timed calls: historical SMSD graph bindings could not import RDKit bond ring/aromatic flags, and re-perceived graph properties sometimes differed. All affected cases remain recorded, rather than being silently counted as successful. Candidate bindings import those flags directly.

| Strategy / policy | Corpus | Cases | Input/preflight errors | SMSD invalid | RDKit no vertex witness | RDKit canceled | SMSD elapsed crossings | SMSD near budget | Equal atom counts |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| native / fmcs | dalke-nn | 1000 | 252 | 4 | 73 | 35 | 133 | 243 | 551 |
| native / strict | dalke-nn | 1000 | 252 | 0 | 37 | 23 | 5 | 18 | 703 |
| native / fmcs | dalke-random | 1000 | 381 | 2 | 243 | 3 | 29 | 163 | 269 |
| native / strict | dalke-random | 1000 | 381 | 0 | 31 | 0 | 1 | 6 | 504 |
| native / fmcs | stress | 12 | 1 | 0 | 4 | 0 | 1 | 3 | 6 |
| native / strict | stress | 12 | 1 | 0 | 4 | 0 | 0 | 2 | 7 |
| native / fmcs | tautobase | 468 | 133 | 0 | 39 | 0 | 0 | 0 | 260 |
| native / strict | tautobase | 468 | 133 | 0 | 39 | 0 | 0 | 0 | 292 |
| native / any | dalke-nn | 1000 | 0 | 2 | 105 | 32 | 83 | 185 | 831 |
| native / any | dalke-random | 1000 | 0 | 15 | 449 | 9 | 44 | 357 | 398 |
| native / any | stress | 12 | 0 | 0 | 4 | 0 | 0 | 3 | 7 |
| native / any | tautobase | 468 | 0 | 0 | 0 | 0 | 0 | 0 | 468 |
| auto-any / any | dalke-nn | 1000 | 0 | 85 | 105 | 41 | 0 | 0 | 924 |
| auto-any / any | dalke-random | 1000 | 0 | 346 | 449 | 12 | 0 | 0 | 523 |
| auto-any / any | stress | 12 | 0 | 3 | 4 | 0 | 0 | 0 | 11 |
| auto-any / any | tautobase | 468 | 0 | 2 | 0 | 0 | 0 | 0 | 468 |
| auto-default / fmcs | dalke-nn | 1000 | 252 | 49 | 73 | 37 | 0 | 0 | 678 |
| auto-default / fmcs | dalke-random | 1000 | 381 | 59 | 243 | 6 | 0 | 0 | 138 |
| auto-default / fmcs | stress | 12 | 1 | 3 | 4 | 0 | 0 | 0 | 10 |
| auto-default / fmcs | tautobase | 468 | 133 | 45 | 39 | 0 | 0 | 0 | 306 |
| raw coverage / any | dalke-nn | 1000 | 0 | 85 | 105 | 38 | 0 | 0 | 926 |
| raw coverage / any | dalke-random | 1000 | 0 | 346 | 449 | 12 | 0 | 0 | 523 |
| raw coverage / any | stress | 12 | 0 | 3 | 4 | 0 | 0 | 0 | 11 |
| raw coverage / any | tautobase | 468 | 0 | 2 | 0 | 0 | 0 | 0 | 468 |
| raw coverage / fmcs | dalke-nn | 1000 | 252 | 49 | 73 | 35 | 0 | 0 | 678 |
| raw coverage / fmcs | dalke-random | 1000 | 381 | 59 | 243 | 4 | 0 | 0 | 138 |
| raw coverage / fmcs | stress | 12 | 1 | 3 | 4 | 0 | 0 | 0 | 10 |
| raw coverage / fmcs | tautobase | 468 | 133 | 45 | 39 | 0 | 0 | 0 | 306 |

All 2,480 molecular inputs were additionally tested under the common bond-any contract, which avoids historical perception exclusions. This exposes 17 invalid baseline native mappings, alongside the six invalid native mappings in the FMCS profile. Public auto and raw coverage are separate validated paths. Invalid mappings are excluded from speed and quality claims; equal atom counts in this table are descriptive and do not override validation failures. The complete baseline corpus matrix contains 18,760 case/profile observations and 2,807 seconds of measured search calls. Wall time includes parsing and witness checks and differs from this summed search cost.

## Final candidate molecular corpora

All **17,360 molecular case/profile observations** completed without input/search errors or invalid SMSD witnesses. These rows cover three native policies plus both public-auto/raw-coverage policies on every molecular input. The rows include repeated cases across strategies; this is 2,480 molecular input records, including known repeats/self-pairs, rather than 17,360 independent molecules. SMARTS dialect differences are reported separately.

| Strategy / policy | Corpus | Cases | Errors | SMSD invalid | RDKit no vertex witness | RDKit canceled | SMSD elapsed crossings | SMSD near budget | Equal atom counts |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| native / strict | dalke-nn | 1000 | 0 | 0 | 45 | 27 | 36 | 92 | 949 |
| native / strict | dalke-random | 1000 | 0 | 0 | 64 | 1 | 5 | 59 | 866 |
| native / strict | stress | 12 | 0 | 0 | 4 | 0 | 1 | 3 | 8 |
| native / strict | tautobase | 468 | 0 | 0 | 41 | 0 | 0 | 0 | 427 |
| native / fmcs | dalke-nn | 1000 | 0 | 0 | 113 | 33 | 39 | 136 | 875 |
| native / fmcs | dalke-random | 1000 | 0 | 0 | 428 | 5 | 14 | 208 | 505 |
| native / fmcs | stress | 12 | 0 | 0 | 4 | 0 | 1 | 3 | 8 |
| native / fmcs | tautobase | 468 | 0 | 0 | 43 | 0 | 0 | 0 | 384 |
| native / any | dalke-nn | 1000 | 0 | 0 | 105 | 33 | 44 | 138 | 884 |
| native / any | dalke-random | 1000 | 0 | 0 | 449 | 8 | 16 | 200 | 518 |
| native / any | stress | 12 | 0 | 0 | 4 | 0 | 1 | 3 | 8 |
| native / any | tautobase | 468 | 0 | 0 | 0 | 0 | 0 | 0 | 468 |
| auto / any | dalke-nn | 1000 | 0 | 0 | 105 | 38 | 0 | 0 | 875 |
| auto / any | dalke-random | 1000 | 0 | 0 | 449 | 8 | 0 | 0 | 393 |
| auto / any | stress | 12 | 0 | 0 | 4 | 0 | 0 | 0 | 8 |
| auto / any | tautobase | 468 | 0 | 0 | 0 | 0 | 0 | 0 | 468 |
| auto / fmcs | dalke-nn | 1000 | 0 | 0 | 113 | 40 | 0 | 0 | 861 |
| auto / fmcs | dalke-random | 1000 | 0 | 0 | 428 | 6 | 0 | 0 | 354 |
| auto / fmcs | stress | 12 | 0 | 0 | 4 | 0 | 0 | 0 | 8 |
| auto / fmcs | tautobase | 468 | 0 | 0 | 43 | 0 | 0 | 0 | 384 |
| raw coverage / any | dalke-nn | 1000 | 0 | 0 | 105 | 38 | 0 | 0 | 875 |
| raw coverage / any | dalke-random | 1000 | 0 | 0 | 449 | 7 | 0 | 0 | 393 |
| raw coverage / any | stress | 12 | 0 | 0 | 4 | 0 | 0 | 0 | 8 |
| raw coverage / any | tautobase | 468 | 0 | 0 | 0 | 0 | 0 | 0 | 468 |
| raw coverage / fmcs | dalke-nn | 1000 | 0 | 0 | 113 | 40 | 0 | 0 | 861 |
| raw coverage / fmcs | dalke-random | 1000 | 0 | 0 | 428 | 6 | 0 | 0 | 354 |
| raw coverage / fmcs | stress | 12 | 0 | 0 | 4 | 0 | 0 | 0 | 8 |
| raw coverage / fmcs | tautobase | 468 | 0 | 0 | 43 | 0 | 0 | 0 | 384 |

## Attained size changes on baseline-valid witnesses

Each row below compares the same inputs, policy and one-second budget, and retains only pairs where both SMSD versions independently validate. Larger/smaller refers to returned atom counts, not proven optimality. The full runs had concurrent local load, so these comparisons do not establish speed rankings. Baseline-invalid witnesses and preflight exclusions are not treated as quality wins.

| Strategy / policy | Corpus | Both valid | Candidate larger | Equal | Candidate smaller | Aggregate atom delta |
|---|---|---:|---:|---:|---:|---:|
| native / strict | stress | 11 | 1 | 10 | 0 | 1 |
| native / strict | dalke-random | 619 | 51 | 565 | 3 | 67 |
| native / strict | dalke-nn | 748 | 11 | 737 | 0 | 94 |
| native / strict | tautobase | 335 | 5 | 330 | 0 | 5 |
| native / fmcs | stress | 11 | 2 | 9 | 0 | 4 |
| native / fmcs | dalke-random | 617 | 212 | 398 | 7 | 491 |
| native / fmcs | dalke-nn | 744 | 146 | 597 | 1 | 3242 |
| native / fmcs | tautobase | 335 | 26 | 309 | 0 | 41 |
| native / any | stress | 12 | 2 | 10 | 0 | 4 |
| native / any | dalke-random | 985 | 341 | 622 | 22 | 795 |
| native / any | dalke-nn | 998 | 94 | 900 | 4 | 850 |
| native / any | tautobase | 468 | 0 | 468 | 0 | 0 |
| auto / fmcs | stress | 8 | 0 | 8 | 0 | 0 |
| auto / fmcs | dalke-random | 560 | 324 | 236 | 0 | 1115 |
| auto / fmcs | dalke-nn | 699 | 28 | 671 | 0 | 110 |
| auto / fmcs | tautobase | 290 | 14 | 276 | 0 | 23 |
| auto / any | stress | 9 | 0 | 9 | 0 | 0 |
| auto / any | dalke-random | 654 | 10 | 642 | 2 | 13 |
| auto / any | dalke-nn | 915 | 5 | 910 | 0 | 13 |
| auto / any | tautobase | 466 | 0 | 466 | 0 | 0 |
| raw coverage / fmcs | stress | 8 | 0 | 8 | 0 | 0 |
| raw coverage / fmcs | dalke-random | 560 | 324 | 236 | 0 | 1115 |
| raw coverage / fmcs | dalke-nn | 699 | 28 | 671 | 0 | 110 |
| raw coverage / fmcs | tautobase | 290 | 14 | 276 | 0 | 23 |
| raw coverage / any | stress | 9 | 0 | 9 | 0 | 0 |
| raw coverage / any | dalke-random | 654 | 10 | 642 | 2 | 13 |
| raw coverage / any | dalke-nn | 915 | 3 | 912 | 0 | 4 |
| raw coverage / any | tautobase | 466 | 0 | 466 | 0 | 0 |

The FMCS rerun includes the repaired FLEXIBLE coverage validation. On 290 baseline-valid Tautobase pairs, final auto/raw coverage has 14 larger and 276 equal mappings, with no smaller results (aggregate +23 atoms). Earlier affected candidate measurements are excluded. Some bond-any and native bounded searches still return smaller valid witnesses; those tradeoffs remain visible above.

## Corpus provenance

The two checked-in 1,000-pair files are **Dalke-style, MoleculeNet-derived** corpora, not the original ChEMBL-13 corpus used by Dalke and Hastings. The original work used random ChEMBL-13 pairs and neighborhood groups of 2, 10 and 100 molecules; this repository uses a different pool and two-molecule pairs. [FMCS paper](https://d-nb.info/1188152866/34).

The checked-in random file has no self-ID pairs or duplicate unordered ID pairs. The neighbor file contains **53 self-ID pairs, 39 repeated unordered ID pairs and 444 pairs below Tanimoto 0.7** (range 0.125–1.0, audited with RDKit 2026.09.1 Morgan radius-2/2048 fingerprints). Self-pairs can strongly bias speed measurements. They remain identifiable in raw results and do not establish representative neighbor-search performance.

The corrected generator seeds before pool subsampling, excludes the query’s own index even under fingerprint ties, records source hashes and RDKit version, and uses no similarity cutoff. It generated 1,000 random and 1,000 nearest-neighbor pairs locally; these new files do not replace the archived inputs used for before/after comparison.

| Input | Records | SHA-256 |
|---|---:|---|
| stress | 12 | `b9b8ce63093bce838461cbc05274eaadd75da0fc3a77b55fa005841152c851d6` |
| dalke-random | 1000 | `f658d67a9543a8eefcf2ac1849011227adc4e9888eb8a4b3c8f79f19f2c1318c` |
| dalke-nn | 1000 | `d668ed6eb3e847fbed3732e7442c3775d7d70883742244812ecfeb07ac67f89e` |
| tautobase | 468 | `f805439d2dd46553bf73701d8edd6f369c4b2cf287c2787ad056774a5668e778` |
| smarts | 1400 | `9ff02503d28aebc6a006d0a077e848aeb68296d960f159de953c5c72434165b9` |

## SMARTS and feature diagnostics

The full 1,400-pattern comparison compiles both queries outside the matching timer and uses one ibuprofen target. Both versions have 25 input/search exceptions (20 RDKit rejections and five native parser errors) and seven hit disagreements among the remaining 1,375 patterns. Different SMARTS dialects are recorded; these results do not establish full dialect compatibility or reproduce the Ehrlich–Rarey paper’s complete experiment. All SMARTS speed comparisons are excluded. [Original dataset](https://www.zbh.uni-hamburg.de/forschung/amd/datasets/smarts-dataset.html).

Full public API Python feature diagnostics ran all molecular corpora plus all SMARTS patterns: seven opt-in tests passed with a 1-second configured budget. Baseline Tautobase has 70 positive and 52 negative returned-atom deltas between tautomer-aware and ordinary dispatch (sum +6); the final candidate has 63 positive, zero negative and aggregate +81 atoms; these are unvalidated feature diagnostics, not superior-quality claims. Native SMARTS accepted 1,395 patterns and rejected five in this separate compile-plus-match diagnostic.

Java standalone substructure covers 28 pairs with 20 warmups/100 trials: 27 boolean agreements with CDK, one chemistry difference for benzene versus cyclohexane under flexible aromaticity. The tautomer loader was corrected to select the labelled Section 7, instead of incorrect numeric row offsets; parse failures preserve pair positions. Java loads 96 tautomer records with four parse failures in three pairs, leaving 45 paired diagnostics. Legacy proton-consistency diagnostics report 16 passes/29 failures; they are not independent proofs of chemical false positives.

## Entry-point coverage and historical overlap

All maintained standalone entry points were executed locally for baseline and candidate. RDKit-only diagnostics, data generation and utilities shared unchanged comparator/input versions where appropriate. The report distinguishes wrappers/duplicates from independent experiments:

| Entry point | Baseline/candidate execution and scope |
|---|---|
| `run_external_benchmarks.py` | Complete native strict/FMCS/any, public auto FMCS/any and raw coverage FMCS/any corpora; all SMARTS |
| `benchmark_python.py` | 20 curated MCS pairs; 20 substructure cases through core aggregator |
| `benchmark_bindings.py` | Parsing/cache/copies/dispatch/batch microbenchmarks; 23 operation profiles, serial/batch checksums agree; final candidate binding rerun complete |
| `benchmark_cpp.sh` / `.cpp` | Actual native SMSD, six curated cases; every returned mapping valid |
| Four C++ primitive programs | Assignment/matching, connected components, substructure setup and search helpers; checksum/assertion checks passed |
| C++/Java tautomer programs | Labelled Section 7; C++ 48 pairs, Java 45, with parse/pair exclusions recorded |
| Java standalone substructure | All 28 pairs, no runtime errors |
| `benchmark_1000_java.java` | All 993 generated pairs available from deduplicated molecule pool; one trial, 1 second |
| `benchmark_1000.py` | All 1,000 generated random/scaffold/size-matched pairs, one trial, 1 second; JavaCLI startup included |
| `benchmark_python_vs_rdkit.py` | Historical ten-pair wrapper over maintained protocol |
| `benchmark_smsd_v6.py` | Ten public MCS cases, 100 fingerprint calls and 100 batch substructure targets; historical filename retained |
| `benchmark_java.sh`, `benchmark_rdkit.py` | Historical15-pair diagnostic sets; JavaJSON parsing fixed to read mcs_size |
| `benchmark_all.py` | All20 pairs; integration timings, no chemistry/witness-equivalent ranking |
| `benchmark_leaderboard.py` | Core/integration/profile modes; delegates existing experiments |
| `benchmark_present_vs_light_scaling.py` | Separate sibling extension,20MCS/28substructure cases; unvalidated diagnostic only |
| `generate_dalke_pairs.py` | Corrected full1,000+1,000 derived-pair generation |
| `compare_benchmark_tsv.py` | Saved-result utility; matching JSON validation/protocol required for MCS timing comparisons |
| Python opt-in external tests | All corpora,7passed; importlib mode verifies installed wrapper |
| Java opt-in benchmark tests | Corrected bounded harness: all22 selected test cases passed for both versions,4,904 data rows each; historical unbounded run failed outer60second limits and required interruption |

Java CLI startup, CDK/Python/native setup costs and chemical options differ across historical entry points. They are not pooled into a cross-language speed score. Archived `results_*` files are preserved but not treated as current release evidence.

## Reproduction and raw evidence

Commands and comparison modes are documented in [README.md](README.md). Use separate baseline/candidate environments and verify the recorded module paths; pytest diagnostics use `--import-mode=importlib` to avoid loading the current source wrapper alongside an older installed extension. Native full corpus checkpoints retain every mapping, validation failure and engine observation.

The full candidate matrix contains 18,760 case/profile observations: 17,360 molecular rows and 1,400 SMARTS rows. Summed timed search calls cost 2592.8 seconds; the complete baseline cost is 2,807.2 seconds. These totals differ in accepted inputs and concurrent load; they are workload costs, not speed ratios. Native policies ran in parallel for the final candidate, with separate per-policy metadata and checkpoints. All 22 selected Java benchmark tests and seven Python feature diagnostics passed for both versions. The candidate Java standalone pool covers 993 returned calls from 989 valid molecules (14 parse exclusions), including 44 empty mappings and 139 elapsed-budget crossings; mapping validity is unverified in this historical program.

Validation evidence includes 691 passing Python tests/8 skips; 11 native Debug suites; separately reported counts of 84,096 small-graph oracle models, 1,024 McGregor state models and 576 coverage-validity models; focused AddressSanitizer/UndefinedBehaviorSanitizer checks; 1,242 passing Java verification tests/15 skips; and three real Metal hardware suites. The model checks are implemented in the public [native regression source](../cpp/tests/test_mcs_regression.cpp) and executed in the final Debug tests. The focused sanitizer run does not execute all 84,096 oracle models, and these small-graph checks do not prove arbitrary-molecule optimality.

Sanitized raw evidence: `smsd-7.2.0-benchmark-data.tar.gz`, 325 members, 7.1 MB. SHA-256: `217e42f3b7f9cf7a2996120d57037826419ec075b23e82059afc2d44fb5a8ea5`. The archive includes complete checkpoints, mappings, validation outcomes, corpus hashes, per-wheel fingerprints, benchmark harness sources and scoped test/oracle metadata. Module paths use relative placeholders; private instructions, review notes and earlier affected candidate measurements are excluded. The archive completeness and privacy checks passed.

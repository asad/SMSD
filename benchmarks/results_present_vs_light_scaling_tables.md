# Present vs Light SMSD: Comparison Tables

Python-to-Python comparison between:

- Present SMSD Pro extension: `/Users/bioinception/tools/github/SMSD/cpp/build-py313/_smsd.cpython-313-darwin.so`
- Lighter SMSD extension: `/Users/bioinception/tools/github/bioinception/python/smsd/_smsd.cpython-313-darwin.so`

Benchmark settings:

- Python: `3.13`
- MCS timeout: `1000 ms`
- Warmup: `1`
- Timed iterations: `3`

## Table 1. MCS Scaling Summary by Size Bucket

| Combined atoms | Cases | Present median (ms) | Light median (ms) | Median speedup |
|---|---:|---:|---:|---:|
| <=16 | 4 | 0.045 | 0.046 | 2.81x |
| 17-32 | 4 | 0.905 | 0.199 | 3.62x |
| 33-64 | 5 | 2.430 | 1.324 | 4.43x |
| 65-128 | 6 | 224.854 | 28.487 | 6.85x |
| 129+ | 1 | 0.072 | 0.013 | 5.73x |

Overall MCS agreement: `19/20`

## Table 2. Per-Pair MCS Comparison

| Pair | Combined atoms | Present (ms) | Light (ms) | Present MCS | Light MCS | Speedup | Agree |
|---|---:|---:|---:|---:|---:|---:|---|
| methane-ethane | 3 | 0.026 | 0.005 | 1 | 1 | 4.99x | Yes |
| benzene-toluene | 13 | 0.068 | 0.087 | 6 | 6 | 0.78x | Yes |
| benzene-phenol | 13 | 0.064 | 0.089 | 6 | 6 | 0.72x | Yes |
| aspirin-acetaminophen | 24 | 2.159 | 0.887 | 7 | 7 | 2.43x | Yes |
| caffeine-theophylline | 27 | 0.149 | 0.063 | 13 | 13 | 2.35x | Yes |
| morphine-codeine | 41 | 2.430 | 4.372 | 19 | 19 | 0.56x | Yes |
| ibuprofen-naproxen | 34 | 6.358 | 1.324 | 12 | 12 | 4.80x | Yes |
| ATP-ADP | 58 | 0.351 | 0.049 | 27 | 27 | 7.09x | Yes |
| NAD-NADH | 88 | 3363.006 | 1134.661 | 24 | 24 | 2.96x | Yes |
| atorvastatin-rosuvastatin | 74 | 156.000 | 23.967 | 10 | 10 | 6.51x | Yes |
| paclitaxel-docetaxel | 120 | 3011.960 | 1104.131 | 53 | 53 | 2.73x | Yes |
| erythromycin-azithromycin | 103 | 23.836 | 1.951 | 49 | 49 | 12.22x | Yes |
| strychnine-quinine | 49 | 3432.750 | 1294.847 | 0 | 9 | 2.65x | No |
| vancomycin-self | 202 | 0.072 | 0.013 | 101 | 101 | 5.73x | Yes |
| adamantane-self | 20 | 0.007 | 0.001 | 10 | 10 | 4.80x | Yes |
| cubane-self | 16 | 0.006 | 0.001 | 8 | 8 | 4.83x | Yes |
| PEG12-PEG16 | 92 | 0.328 | 0.046 | 40 | 40 | 7.19x | Yes |
| coronene-self | 48 | 0.013 | 0.003 | 24 | 24 | 4.43x | Yes |
| guanine-keto-enol | 22 | 1.660 | 0.334 | 10 | 10 | 4.97x | Yes |
| rdkit-1585-pair | 68 | 293.707 | 33.007 | 17 | 17 | 8.90x | Yes |

## Table 3. Substructure Scaling Summary by Size Bucket

| Combined atoms | Cases | Present median (ms) | Light median (ms) | Median speedup | Agreement |
|---|---:|---:|---:|---:|---:|
| <=16 | 19 | 0.0043 | 0.0008 | 5.04x | 19/19 |
| 17-32 | 8 | 0.0048 | 0.0009 | 5.20x | 8/8 |
| 33-64 | 1 | 0.0071 | 0.0014 | 5.10x | 1/1 |

Overall substructure agreement: `28/28`

## Table 4. Representative Substructure Comparison Cases

| Pair | Combined atoms | Present (ms) | Light (ms) | Present hit | Light hit | Speedup | Agree |
|---|---:|---:|---:|---|---|---:|---|
| Piperidine in Fentanyl | 31 | 0.0131 | 0.0027 | Yes | Yes | 4.93x | Yes |
| Glucose in Sucrose | 35 | 0.0071 | 0.0014 | Yes | Yes | 5.10x | Yes |
| Cyclobutane in Cyclopentane | 9 | 0.0068 | 0.0014 | No | No | 4.79x | Yes |
| Sulfonamide in Sulfisoxazole | 26 | 0.0061 | 0.0012 | Yes | Yes | 5.19x | Yes |
| Acetanilide in Acetaminophen | 21 | 0.0056 | 0.0016 | Yes | Yes | 3.49x | Yes |
| Pyridine in Niacin | 15 | 0.0053 | 0.0010 | Yes | Yes | 5.20x | Yes |
| Phenol in Acetaminophen | 18 | 0.0050 | 0.0010 | Yes | Yes | 5.20x | Yes |
| Benzene in Phenol | 13 | 0.0047 | 0.0007 | Yes | Yes | 6.52x | Yes |
| Piperazine in Phenylpiperazine | 18 | 0.0047 | 0.0009 | Yes | Yes | 5.38x | Yes |
| Methane in Ethane | 3 | 0.0047 | 0.0014 | Yes | Yes | 3.30x | Yes |

## Notes

- These tables are directly derived from:
  - `/Users/bioinception/tools/github/SMSD/benchmarks/results_present_vs_light_scaling_mcs.tsv`
  - `/Users/bioinception/tools/github/SMSD/benchmarks/results_present_vs_light_scaling_sub.tsv`
- The only MCS disagreement in this run was `strychnine-quinine`, where the present build returned `0` atoms and the lighter build returned `9` under the tighter `1000 ms` timeout.
- If you want manuscript-ready tables, the next step is to rerun the hard MCS tail with a longer timeout and export the same tables as LaTeX.

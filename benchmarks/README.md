# SMSD Cross-Platform Benchmarks

Publication-quality benchmark suite comparing MCS (Maximum Common Subgraph)
performance across three implementations:

| Engine | Language | Invocation |
|--------|----------|------------|
| **SMSD** | Java (CDK) | CLI jar with `--mode mcs --json -` |
| **RDKit FindMCS** | Python (C++ backend) | `rdFMCS.FindMCS()` |
| **SMSD C++** | C++17 header-only | Standalone binary or via `benchmark_cpp.sh` |

All benchmarks use the same **20 molecule pairs** covering trivial to
computationally challenging cases.

## Quick Start

```bash
# 1. Build SMSD Java jar (one-time)
cd /path/to/SMSD
mvn package -DskipTests

# 2. Install RDKit
pip install rdkit

# 3. Run the maintained local core leaderboard
python3 benchmarks/benchmark_leaderboard.py --mode core
```

Output files:
- `benchmarks/results_leaderboard_core_mcs.tsv` -- machine-readable MCS comparison
- `benchmarks/results_leaderboard_core_python.tsv` -- Python-core detailed timings
- `benchmarks/results_leaderboard_core_java_sub.tsv` -- Java cached-substructure timings
- `benchmarks/results_leaderboard_core_summary.txt` -- formatted summary

## Running Individual Benchmarks

### RDKit only

```bash
pip install rdkit
python3 benchmarks/benchmark_rdkit.py
```

### SMSD Java only

```bash
mvn package -DskipTests
bash benchmarks/benchmark_java.sh
```

### C++ only

```bash
bash benchmarks/benchmark_cpp.sh
```

This script handles CMake configuration, compilation in Release mode, and
execution. Requires CMake 3.16+ and a C++17 compiler.

### Unified (all engines)

```bash
python3 benchmarks/benchmark_all.py
```

Runs RDKit and SMSD Java head-to-head on all 20 pairs. If the C++ binary
is compiled (via `benchmark_cpp.sh`), it is automatically included.

### Core leaderboard

```bash
python3 benchmarks/benchmark_leaderboard.py --mode core --compare-mode strict
```

This is the maintained day-to-day benchmark for local core changes. It runs:

- in-process Python core vs RDKit for the 20 curated MCS pairs
- Java cached-path substructure on the 28 curated query/target pairs

For quicker turnaround while iterating on the native core, reduce the Python
MCS measured runs:

```bash
python3 benchmarks/benchmark_leaderboard.py --mode core --compare-mode strict --py-iters 5
```

Comparison modes:

- `--compare-mode defaults`: toolkit defaults vs toolkit defaults
- `--compare-mode strict`: exact/ring-parity baseline for the Python SMSD vs RDKit MCS comparison
- `--compare-mode fmcs`: loose FMCS-style baseline for the Python SMSD vs RDKit MCS comparison

The Java/CDK substructure benchmark is still a best-effort structural comparison rather
than a perfect flag-for-flag chemistry match, because CDK's matcher API does not expose
all of the same ring/topology controls as SMSD.

## Local Before/After Workflow

To benchmark your local Python/C++ changes rather than an already-installed
wheel, use the Python benchmark with an explicit output path and the local
import guard:

```bash
# 1. Build/install the current working tree in editable mode
python3 -m pip install -e .

# 2. Capture a baseline
python3 benchmarks/benchmark_python.py \
  --mcs-only \
  --require-local-smsd \
  --print-smsd-path \
  --output benchmarks/runs/mcs_before.tsv

# 3. Make changes, rebuild editable install if needed
python3 -m pip install -e .

# 4. Capture the candidate run
python3 benchmarks/benchmark_python.py \
  --mcs-only \
  --require-local-smsd \
  --print-smsd-path \
  --output benchmarks/runs/mcs_after.tsv

# 5. Compare the two runs
python3 benchmarks/compare_benchmark_tsv.py \
  benchmarks/runs/mcs_before.tsv \
  benchmarks/runs/mcs_after.tsv
```

Notes:

- `--require-local-smsd` fails fast if `import smsd` or `smsd._smsd` resolves
  outside this repository.
- `--output` prevents accidental overwriting of the default TSV.
- `--mcs-only` is the right mode for changes in `mcs.hpp`; substructure timings
  mainly reflect `vf2pp.hpp`.
- The default chemistry profile uses strict ring parity. If a benchmark is
  intended to mimic loose FMCS behavior, disable that explicitly with
  `ChemOptions.fmcsProfile()` / `ringMatchesRingOnly=false`.

## Molecule Pairs (20)

| # | Pair | Category | Notes |
|---|------|----------|-------|
| 1 | Methane / Ethane | Trivial | 1-2 atoms |
| 2 | Benzene / Toluene | Small aromatic | Methyl substituent |
| 3 | Benzene / Phenol | Heteroatom | OH substituent |
| 4 | Aspirin / Acetaminophen | Drug pair | Different scaffolds |
| 5 | Caffeine / Theophylline | N-methyl diff | Single N-CH3 removal |
| 6 | Morphine / Codeine | Alkaloid | Complex ring system |
| 7 | Ibuprofen / Naproxen | NSAID | Scaffold change |
| 8 | ATP / ADP | Nucleotide | Phosphate chain |
| 9 | NAD+ / NADH | Cofactor | Redox pair, large |
| 10 | Atorvastatin / Rosuvastatin | Statin | Different heterocycles |
| 11 | Paclitaxel / Docetaxel | Taxane | ~50 atoms, complex |
| 12 | Erythromycin / Azithromycin | Macrolide | Macrocyclic ring |
| 13 | Strychnine / Quinine | Alkaloid scaffold | Different frameworks |
| 14 | Vancomycin / self | Self-match large | ~100 atoms |
| 15 | Adamantane / self | Symmetric | Cage, high symmetry |
| 16 | Cubane / self | Cage | 8-atom cube |
| 17 | PEG-12 / PEG-16 | Polymer | Repetitive chain |
| 18 | Coronene / self | PAH | Polycyclic aromatic |
| 19 | Guanine keto / enol | Tautomer | Tautomeric forms |
| 20 | RDKit #1585 pair | Known failure | Stress test |

## Output Format

The unified benchmark (`benchmark_all.py`) produces a comparison table:

```
 #  Pair                          Category            RDKit ms  R.MCS   SMSD ms  S.MCS  Winner  Speed    Quality
------------------------------------------------------------------------------------------------------------------
 1  methane-ethane                 Trivial                0.012      1     0.008      1  SMSD      1.5x  equal
 2  benzene-toluene                Small aromatic         0.045      6     0.023      6  SMSD      2.0x  equal
...
```

Columns:
- **RDKit ms / SMSD ms**: Best time across 5 runs (lower is better)
- **R.MCS / S.MCS**: MCS size in atoms (higher is better)
- **Winner**: Which engine was faster
- **Speed**: Speedup factor
- **Quality**: Whether MCS sizes agree

## Prerequisites

| Benchmark | Requirements |
|-----------|-------------|
| Unified (`benchmark_all.py`) | Python 3.8+, `rdkit`, Java 11+, built SMSD jar |
| Java only | Java 11+, Maven (to build SMSD jar) |
| Python only | Python 3.8+, `pip install rdkit` |
| C++ only | CMake 3.16+, C++17 compiler |

## Methodology

- **Warmup**: 1 warmup run per pair (discarded) to eliminate JIT and cache effects
- **Timed runs**: 5 runs per pair; best and median times recorded
- **Timeout**: 10 seconds per pair (both engines)
- **Timing**: `time.perf_counter()` (Python), wall-clock subprocess timing (Java CLI)
- **MCS mode**: Atom-count-based MCS (not bond-count)

## Notes

- Java CLI timing includes JVM startup overhead per pair invocation. For
  in-process Java benchmarks, use the `--benchmark` CLI flag.
- The C++ benchmark covers a subset of pairs (those buildable without a SMILES
  parser). It measures raw graph algorithm speed.
- RDKit `FindMCS` returns a SMARTS pattern; SMSD returns atom-pair mappings.
  MCS sizes should agree for correct implementations.
- An asterisk (`*`) after RDKit time indicates the timeout was hit.

---

## Substructure Search Benchmark

A separate benchmark suite for **substructure search** (subgraph isomorphism),
comparing SMSD against CDK's `DfPattern` in a pure Java-to-Java setting.

### Files

| File | Description |
|------|-------------|
| `substructure_pairs.tsv` | 28 query/target SMILES pairs across 7 test sections |
| `benchmark_substructure_java.java` | Java benchmark: SMSD vs CDK DfPattern |
| `results_substructure.tsv` | Output with per-pair timings and speedup factors |

### Test Data (`substructure_pairs.tsv`)

28 query/target pairs organised into 7 sections:

| Section | Pairs | Description |
|---------|-------|-------------|
| Trivial / chain | 1-4 | Methane in ethane, propane chain, ethanol substructure |
| Ring containment | 5-8 | Benzene in toluene, naphthalene, biphenyl |
| Aromatic | 9-12 | Pyridine in nicotine, indole in tryptophan |
| Heterocyclic drugs | 13-16 | Imidazole in histamine, purine in caffeine |
| Drug fragments | 17-20 | Phenol in acetaminophen, piperidine in fentanyl |
| Edge cases | 21-24 | Self-match, single atom, symmetric cages |
| Negative controls | 25-28 | Pairs with no substructure match expected |

Format: `query_smiles<TAB>target_smiles<TAB>label`

All SMILES use kekulized form for CDK 2.11 compatibility.

### Running the Benchmark (`benchmark_substructure_java.java`)

Compares three approaches:

1. **SMSD IAtomContainer API** -- standard atom-container interface
2. **SMSD cached MolGraph API** -- pre-built graph, avoids repeated setup
3. **CDK DfPattern** -- CDK's built-in substructure matcher

```bash
# Set up Java (adjust path as needed)
export JAVA_HOME=/usr/local/Cellar/openjdk/25.0.2/libexec/openjdk.jdk/Contents/Home
export PATH="$JAVA_HOME/bin:$PATH"

# Build SMSD jar
mvn package -DskipTests

# Compile and run
javac -cp target/smsd-6.8.1-jar-with-dependencies.jar benchmarks/benchmark_substructure_java.java -d benchmarks/
java -cp target/smsd-6.8.1-jar-with-dependencies.jar:benchmarks/ benchmark_substructure_java
```

### Results (`results_substructure.tsv`)

Output columns:

| Column | Description |
|--------|-------------|
| `Pair` | Pair label from test data |
| `SMSD_us` | SMSD IAtomContainer time (microseconds) |
| `SMSD_cached_us` | SMSD cached MolGraph time (microseconds) |
| `CDK_us` | CDK DfPattern time (microseconds) |
| `Speedup` | CDK time / SMSD time |
| `Cached_Speedup` | CDK time / SMSD cached time |
| `SMSD_match` | Whether SMSD found a match (true/false) |
| `CDK_match` | Whether CDK found a match (true/false) |

**Key result**: SMSD cached MolGraph API beats CDK DfPattern on all 28/28
pairs, with speedups ranging from 2.5x to 13x.

---

## External Benchmark Datasets

Community-standard datasets for large-scale reproducible evaluation.
All data in `benchmarks/data/`. See [`data/README.md`](data/README.md) for download sources and licences.

### Datasets

| Dataset | File | Records | Source | Purpose |
|---------|------|--------:|--------|---------|
| Tautobase (Chodera subset) | `chodera_tautobase_subset.txt` | 468 pairs | Wahl & Sander 2020 | Tautomer-aware MCS |
| Tautobase (full SMIRKS) | `tautobase_smirks.txt` | 1,680 pairs | Wahl & Sander 2020 | Tautomer coverage |
| Ehrlich-Rarey SMARTS v2.0 | `ehrlich_rarey_smarts.txt` | 1,400 patterns | Ehrlich & Rarey 2012, ZBH Hamburg | Substructure search |
| Dalke random pairs | `dalke_random_pairs.tsv` | 1,000 pairs | MoleculeNet drug collections | Low-similarity MCS |
| Dalke nearest-neighbor pairs | `dalke_nn_pairs.tsv` | 1,000 pairs | MoleculeNet drug collections | High-similarity MCS |
| Stress pairs | `stress_pairs.tsv` | 12 pairs | Duesbury et al. 2017 | Adversarial robustness |
| Molecule pool | `chembl_mcs_benchmark.smi` | 5,590 molecules | MoleculeNet (BBBP, SIDER, ClinTox, BACE) | Pair generation |

### External Benchmark Results

#### Stress Test (12 adversarial hard cases)

| Pair | MCS | Time | Challenge |
|------|----:|-----:|-----------|
| Cubane / Cuneane | 8 | 2.9 ms | Symmetric cages, orbit pruning |
| Anthracene / Phenanthrene | 14 | 0.4 ms | Fused aromatic isomers |
| Crown ether 15-C-5 / 18-C-6 | 15 | 1.3 ms | Macrocyclic symmetry |
| Cyclohexadecane / Cyclooctadecane | 16 | 3.2 ms | Large mono-elemental rings |
| Naphthalene / Anthracene | 10 | 0.3 ms | Fused PAH size mismatch |
| Cubane (self) | 8 | 0.4 ms | Symmetric self-match |
| Coronene (self) | 24 | 1.9 ms | Large PAH (24 atoms) |
| Neopentane (self) | 5 | 0.5 ms | Highly branched |
| Adamantane (self) | 11 | 0.2 ms | Cage self-match |
| Icosane / Hexadecane | 16 | 0.7 ms | Long chain polymer-like |
| Kekulene (self) | 48 | 0.6 ms | Massive symmetric PAH (48 atoms) |
| Cyclopeptide / mutant | 50 | 2,026 ms | Dense cyclic peptide (>50 atoms) |

**12/12 pairs complete, 0 timeouts.** All self-match pairs return exact atom count.

#### Tautobase Tautomer Pairs (468 pairs)

| Metric | Value |
|--------|-------|
| Pairs tested | 468 |
| Full match (MCS = all atoms) | 437 (93.4%) |
| Median time per pair | 309 us |
| Timeouts | 0 |

**93.4% of tautomer pairs achieve a full heavy-atom match** without tautomer enumeration.

#### Dalke-style MCS Benchmark

| Track | Pairs | Timeouts | Median time | Mean MCS |
|-------|------:|:--------:|------------:|---------:|
| Random (low similarity) | 200 | 0 | 792 ms | 14.4 atoms |
| Nearest-neighbor (high similarity) | 200 | 0 | 774 us | 25.7 atoms |

High-similarity pairs are **1,000x faster** than random pairs (774 us vs 792 ms),
confirming that the LFUB certificate and polynomial fast-paths short-circuit the search.

#### Ehrlich-Rarey SMARTS Substructure (500 of 1,400 patterns)

| Metric | Value |
|--------|-------|
| Patterns tested | 375 |
| Matched against ibuprofen | 39 |
| Unsupported SMARTS (recursive) | 125 |
| Median query time | 6 us |

### Running External Benchmarks

```bash
# Full external benchmark suite (Python, ~10 min)
python benchmarks/run_external_benchmarks.py

# External benchmark test suite with assertions
SMSD_BENCHMARK=1 pytest python/tests/test_external_benchmarks.py -v -s

# Java external benchmarks
mvn test -Dtest=ExternalBenchmarkTest -Dbenchmark=true

# Regenerate Dalke-style pairs (requires RDKit)
python benchmarks/generate_dalke_pairs.py
```

### References

1. Wahl O, Sander T. Tautobase: An Open Tautomer Database. *J Chem Inf Model*. 2020;60(3):1085-1089.
2. Ehrlich HC, Rarey M. Systematic benchmark of substructure search in molecular graphs. *J Cheminform*. 2012;4:13.
3. Dalke A, Hastings J. FMCS: a novel algorithm for the multiple MCS problem. *J Cheminform*. 2013;5:6.
4. Duesbury E, Holliday J, Sheridan RP. Comparison of MCS isomorphism algorithms. *ChemMedChem*. 2017;12(21):1687-1694.
5. Wu Z, Ramsundar B, et al. MoleculeNet: a benchmark for molecular machine learning. *Chem Sci*. 2018;9(2):513-530.

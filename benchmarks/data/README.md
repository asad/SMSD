# External Benchmark Datasets

Community-standard datasets for reproducible evaluation of SMSD Pro.
Each dataset is associated with a specific published benchmark protocol.

---

## Dataset → Result File → Author / Paper

| # | Dataset | Input file | Result file | Pairs / Patterns | Author(s) | Reference |
|---|---------|-----------|-------------|-----------------|-----------|-----------|
| 1 | **Stress pairs** | `stress_pairs.tsv` | `results_stress_12pairs.txt` | 12 adversarial cases | Duesbury, Mayfield, Sheridan | *ChemMedChem* 2017, 12(21):1687–1694 |
| 2 | **Tautobase (Chodera subset)** | `chodera_tautobase_subset.txt` | `results_tautobase_468pairs.txt` | 468 tautomer pairs | Wahl & Sander | *J. Chem. Inf. Model.* 2020, 60(3):1085–1089 |
| 3 | **Dalke random pairs** | `dalke_random_pairs.tsv` | `results_dalke_random_1000pairs.txt` | 1,000 low-similarity pairs | Dalke & Hastings | *J. Cheminform.* 2013, 5(Suppl 1):O6 (FMCS) |
| 4 | **Dalke nearest-neighbor pairs** | `dalke_nn_pairs.tsv` | `results_dalke_nn_1000pairs.txt` | 1,000 high-similarity pairs | Dalke & Hastings | *J. Cheminform.* 2013, 5(Suppl 1):O6 (FMCS) |
| 5 | **Ehrlich-Rarey SMARTS v2.0** | `ehrlich_rarey_smarts.txt` | `results_ehrlich_rarey_1400smarts.txt` | 1,400 SMARTS patterns | Ehrlich & Rarey | *J. Cheminform.* 2012, 4:13 |
| — | Combined summary | — | `results_external_all_combined.txt` | All 5 above | — | — |

### Python MCS + Substructure (SMSD vs RDKit)

| Input | Result file | N pairs | Authors / Protocol |
|-------|------------|---------|-------------------|
| 20 curated SMILES pairs | `results_smsd_vs_rdkit_20pairs_mcs_sub.tsv` | 20 | Internal, spans trivial → vancomycin |

### Java CDK vs SMSD Substructure

Run: `mvn test -Dtest=JavaCdkVsSmsdBenchmarkTest -Dbenchmark=true`
Output: JUnit stdout (CDK 2.12 DfPattern vs SMSD Pro VF2++, 20 pairs)

---

## Input Data Files

| File | Lines | Contents | Source |
|------|-------|----------|--------|
| `stress_pairs.tsv` | 12 | Adversarial SMILES pairs (cubane, PAH, macrocycle) | BioInception curated |
| `chodera_tautobase_subset.txt` | 468 | Tautomer pair SMILES (name, smi1, smi2, ddG) | Chodera Lab / Wahl & Sander 2020 |
| `tautobase_smirks.txt` | 1,680 | Full Tautobase SMIRKS transforms | Wahl & Sander 2020 |
| `dalke_random_pairs.tsv` | 1,000 | Random low-similarity drug pairs | Generated from ChEMBL subset below |
| `dalke_nn_pairs.tsv` | 1,000 | Nearest-neighbor high-similarity pairs | Generated from ChEMBL subset below |
| `ehrlich_rarey_smarts.txt` | 1,400 | SMARTS patterns for substructure search | ZBH Hamburg (Ehrlich & Rarey) |
| `chembl_mcs_benchmark.smi` | 5,590 | Drug-like molecules (source pool) | MoleculeNet: BBBP + SIDER + ClinTox + BACE |
| `BBBP.csv` | — | Blood-brain barrier permeability | MoleculeNet |
| `bace.csv` | — | BACE-1 inhibitors | MoleculeNet |
| `clintox.csv` | — | Clinical toxicity | MoleculeNet |
| `sider.csv` | — | Side-effect resource | MoleculeNet |
| `bzr.sdf` | — | Benzodiazepine receptor ligands | DUD-E |
| `cdk2.sdf` | — | CDK2 kinase inhibitors | DUD-E |

---

## How to regenerate Dalke-style pairs

```bash
# Requires RDKit
python benchmarks/generate_dalke_pairs.py
# Reads: data/chembl_mcs_benchmark.smi
# Writes: data/dalke_random_pairs.tsv  (1,000 random pairs, seed=42)
#         data/dalke_nn_pairs.tsv       (1,000 NN pairs, Tanimoto >= 0.7)
```

---

## How to run all external benchmarks

### Python (SMSD vs RDKit, MCS + Substructure)
```bash
python benchmarks/benchmark_python.py
# Output: benchmarks/results_smsd_vs_rdkit_20pairs_mcs_sub.tsv
```

### Python (external community datasets)
```bash
python benchmarks/run_external_benchmarks.py
# Output: benchmarks/results_stress_12pairs.txt
#         benchmarks/results_tautobase_468pairs.txt
#         benchmarks/results_dalke_random_1000pairs.txt
#         benchmarks/results_dalke_nn_1000pairs.txt
#         benchmarks/results_ehrlich_rarey_1400smarts.txt
#         benchmarks/results_external_all_combined.txt
```

### Java (CDK 2.12 vs SMSD Pro VF2++, Substructure)
```bash
mvn test -Dtest=JavaCdkVsSmsdBenchmarkTest -Dbenchmark=true
```

### Java (External: Tautobase, Dalke, Ehrlich-Rarey, Stress)
```bash
mvn test -Dtest=ExternalBenchmarkTest -Dbenchmark=true
```

---

## References

1. **Wahl & Sander 2020** — Tautobase: An Open Tautomer Database.
   *J. Chem. Inf. Model.* 60(3):1085–1089. DOI: 10.1021/acs.jcim.0c00035

2. **Ehrlich & Rarey 2012** — Systematic benchmark of substructure search in molecular graphs.
   *J. Cheminform.* 4:13. DOI: 10.1186/1758-2946-4-13

3. **Dalke & Hastings 2013** — FMCS: a novel algorithm for the multiple MCS problem.
   *J. Cheminform.* 5(Suppl 1):O6. DOI: 10.1186/1758-2946-5-S1-O6

4. **Duesbury, Mayfield & Sheridan 2017** — Comparison of maximum common subgraph isomorphism algorithms.
   *ChemMedChem* 12(21):1687–1694. DOI: 10.1002/cmdc.201700316

5. **McCreesh, Prosser & Trimble 2017** — A partitioning algorithm for maximum common subgraph problems.
   *Proc. IJCAI-2017*, pp. 712–719.

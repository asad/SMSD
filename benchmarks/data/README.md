# Benchmark data

The [7.2.0 benchmark report](../RESULTS_7.2.0.md) records input hashes, parsing,
matching rules and measured results. These files are repository test corpora;
they do not all reproduce their named publications' original datasets.

## Checked-in records

Counts exclude comments and CSV/TSV headers and precede parsing.

| File | Records | Source |
|---|---:|---|
| `stress_pairs.tsv` | 12 | Curated cage, ring, chain and peptide examples |
| `chodera_tautobase_subset.txt` | 468 | Tautobase-labelled tautomer subset; extraction metadata is unavailable |
| `tautobase_smirks.txt` | 1,680 | Tautobase-labelled transforms and source references |
| `dalke_random_pairs.tsv` | 1,000 | Derived from the MoleculeNet pool below |
| `dalke_nn_pairs.tsv` | 1,000 | Legacy neighbouring pairs from the same pool |
| `ehrlich_rarey_smarts.txt` | 1,400 | Deduplicated SMARTS; attribution, source URL and licence are in its header |
| `chembl_mcs_benchmark.smi` | 5,590 | MoleculeNet-derived pool |
| `BBBP.csv`, `bace.csv`, `clintox.csv`, `sider.csv` | — | Source collections for that pool |
| `bzr.sdf` | 163 | Legacy SDF collection; extraction metadata is unavailable |
| `cdk2.sdf` | 47 | Legacy SDF collection; extraction metadata is unavailable |

The derived pairs are separate from the original Dalke/Hastings ChEMBL-13
corpus. The legacy neighbouring file contains 53 self-pairs, 39 duplicate
unordered source-ID pairs and 444 recorded Tanimoto scores below 0.7.
Compare newly generated pairs separately. SMARTS results on the supplied
ibuprofen target do not reproduce the full Ehrlich–Rarey experiment.

## Generate pairs

RDKit is required. This small example creates ten pairs from up to 50 molecules:

```bash
python benchmarks/generate_dalke_pairs.py \
  --input benchmarks/data/chembl_mcs_benchmark.smi \
  --output-dir build/local-benchmarks/generated-pairs \
  --seed 42 --pairs 10 --max-molecules 50
```

Use `--pairs 1000` and a larger `--max-molecules` for a larger set. The generator
records the RDKit version and source hash, excludes the query index from
neighbour selection, and uses radius-2/1,024-bit Morgan Tanimoto without a
similarity cutoff.

See the [benchmark guide](../README.md) for running comparisons. Java's CDK
comparison measures substructure search; CDK 2.13 is pinned in `java/pom.xml`.

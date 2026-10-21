# Benchmark data

These files are repository test corpora; they do not all reproduce their named
publications' original datasets. The [benchmark guide](../README.md) sets out
the protocol under which they are measured, and each published benchmark report
records the input hashes, parsing, matching rules and measured results of its
release.

## Checked-in records

Counts are taken before parsing and exclude comments and CSV/TSV headers.

| File | Records | Contents | Source and licence |
|---|---:|---|---|
| `mcs_pairs.tsv` | 20 | Curated MCS pairs from trivial to glycopeptide, shared by the Python and integration drivers | Curated for this repository; Apache-2.0 |
| `stress_pairs.tsv` | 12 | Curated cage, ring, chain and peptide examples | Curated for this repository after Duesbury et al. (2017); Apache-2.0 |
| `chodera_tautobase_subset.txt` | 468 | Tautobase-labelled tautomer subset | Subset of Tautobase (Wahl and Sander, 2020); extraction source and licence not recorded |
| `tautobase_smirks.txt` | 1,680 | Tautobase-labelled transforms and source references | Tautobase (Wahl and Sander, 2020), <https://github.com/WahlOya/Tautobase>; the source states no licence |
| `dalke_random_pairs.tsv` | 1,000 | Random pairs from the MoleculeNet pool below | Generated from `chembl_mcs_benchmark.smi`; terms as for that pool |
| `dalke_nn_pairs.tsv` | 1,000 | Nearest-neighbour pairs from the same pool | Generated from `chembl_mcs_benchmark.smi`; terms as for that pool |
| `ehrlich_rarey_smarts.txt` | 1,400 | Deduplicated SMARTS; attribution and source URL are in its header | Ehrlich and Rarey, J. Cheminform. 4:13 (2012), an article published under CC BY 2.0; the dataset page states no separate licence |
| `chembl_mcs_benchmark.smi` | 4,687 | MoleculeNet-derived pool | Structures from the MoleculeNet (Wu et al., 2018) BBBP, BACE and ClinTox collections, each occurring, as written, in at least one of them; licences not recorded |
| `BBBP.csv`, `bace.csv`, `clintox.csv` | — | Source collections for that pool | Processed copies of the MoleculeNet BBBP, BACE and ClinTox collections (Wu et al., 2018); licence not recorded. The SIDER collection is not distributed here |
| `bzr.sdf` | 163 | Benzodiazepine receptor ligands with activity values (SDF) | Identical to `Docs/Book/data/bzr.sdf` of RDKit, <https://github.com/rdkit/rdkit>, at its last revision before 10 October 2026; BSD 3-Clause, see `RDKIT_LICENSE.txt` |
| `cdk2.sdf` | 47 | CDK2 ligands with 3D coordinates and cluster labels (SDF) | Identical to `Docs/Book/data/cdk2.sdf` of RDKit, <https://github.com/rdkit/rdkit>; BSD 3-Clause, see `RDKIT_LICENSE.txt` |

These files are not part of the Python source distribution.

The derived pairs are separate from the original Dalke/Hastings ChEMBL-13
corpus. Both files draw on the 4,664 pool structures that have at least five
heavy atoms. The nearest-neighbour file has no self-pairs; it contains 53
pairs that recur in reverse order, where two queries are each other's nearest
neighbour, 35 pairs of different SMILES with a Tanimoto score of 1 and 461
scores below 0.7. SMARTS results on the supplied ibuprofen target do not
reproduce the full Ehrlich–Rarey experiment.

## Generate pairs

RDKit is required. This small example creates ten pairs from up to 50 molecules:

```bash
python benchmarks/generate_dalke_pairs.py \
  --input benchmarks/data/chembl_mcs_benchmark.smi \
  --output-dir build/local-benchmarks/generated-pairs \
  --seed 42 --pairs 10 --max-molecules 50
```

The checked-in pair files were generated from the default input with the
default arguments, `--seed 42 --pairs 1000 --max-molecules 5000`, under RDKit
2026.03.6; the same arguments and RDKit release reproduce them byte for byte.
The generator records the RDKit version and source hash, excludes the query
index from neighbour selection, and uses radius-2/1,024-bit Morgan Tanimoto
without a similarity cut-off.

See the [benchmark guide](../README.md) for running comparisons. Java's CDK
comparison measures substructure search; CDK 2.13 is pinned in `java/pom.xml`.

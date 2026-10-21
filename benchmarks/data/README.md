# Benchmark data

These files are the inputs of the [benchmarks](../README.md). Some are subsets
or derivatives of published datasets and do not reproduce those datasets in
full. They are not part of the Python source distribution.

| File | Contents | Source and licence |
|---|---|---|
| `mcs_pairs.tsv` | Curated MCS pairs, from trivial to glycopeptide | Curated for this repository; Apache-2.0 |
| `stress_pairs.tsv` | Curated cage, ring, chain and peptide examples | Curated for this repository after Duesbury et al. (2017); Apache-2.0 |
| `chodera_tautobase_subset.txt` | Tautobase-labelled tautomer pairs | Subset of Tautobase (Wahl and Sander, 2020); extraction source and licence not recorded |
| `tautobase_smirks.txt` | Tautobase-labelled transforms and source references | Tautobase (Wahl and Sander, 2020), <https://github.com/WahlOya/Tautobase>; the source states no licence |
| `dalke_random_pairs.tsv` | Random molecule pairs in the style of the Dalke benchmark | Generated from `chembl_mcs_benchmark.smi` by `generate_dalke_pairs.py`; not the original Dalke/Hastings ChEMBL corpus; terms as for that file |
| `dalke_nn_pairs.tsv` | Nearest-neighbour molecule pairs in the same style | As for `dalke_random_pairs.tsv` |
| `ehrlich_rarey_smarts.txt` | SMARTS patterns; attribution and source URL are in its header | Ehrlich and Rarey, J. Cheminform. 4:13 (2012), an article published under CC BY 2.0; the dataset page states no separate licence |
| `chembl_mcs_benchmark.smi` | Molecule pool for the generated pairs | Structures from the MoleculeNet (Wu et al., 2018) BBBP, BACE and ClinTox collections; licences not recorded |
| `BBBP.csv`, `bace.csv`, `clintox.csv` | Source collections for that pool | Processed copies of the MoleculeNet BBBP, BACE and ClinTox collections (Wu et al., 2018); licence not recorded |
| `bzr.sdf` | Benzodiazepine receptor ligands with activity values | `Docs/Book/data/bzr.sdf` of RDKit, <https://github.com/rdkit/rdkit>; BSD 3-Clause, see `RDKIT_LICENSE.txt` |
| `cdk2.sdf` | CDK2 ligands with 3D coordinates and cluster labels | `Docs/Book/data/cdk2.sdf` of RDKit, <https://github.com/rdkit/rdkit>; BSD 3-Clause, see `RDKIT_LICENSE.txt` |

`generate_dalke_pairs.py` reproduces `dalke_random_pairs.tsv` and
`dalke_nn_pairs.tsv` with its default arguments under RDKit 2026.03.6. This
small example (RDKit required) writes ten pairs from up to 50 molecules:

```bash
python benchmarks/generate_dalke_pairs.py \
  --input benchmarks/data/chembl_mcs_benchmark.smi \
  --output-dir build/local-benchmarks/generated-pairs \
  --seed 42 --pairs 10 --max-molecules 50
```

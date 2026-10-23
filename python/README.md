<p align="center">
  <a href="https://github.com/asad/SMSD" aria-label="SMSD Pro">
    <img src="https://raw.githubusercontent.com/asad/SMSD/master/icons/icon.svg" alt="SMSD Pro" width="140"/>
  </a>
</p>

# SMSD Pro™ for Python

[![PyPI](https://img.shields.io/pypi/v/smsd)](https://pypi.org/project/smsd/)
[![Licence](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/asad/SMSD/blob/master/LICENSE)

SMSD Pro provides substructure search, maximum common substructure (MCS)
search, fingerprints and similarity screening for chemical graphs. Core matching
needs no RDKit, CDK or Java; RDKit is optional, for RDKit molecules and the
RDKit image helpers.

## Installation

```bash
python -m pip install "smsd==8.0.2"
```

To include RDKit, install the `rdkit` extra:

```bash
python -m pip install "smsd[rdkit]==8.0.2"
```

The CPU/OpenMP wheels cover CPython 3.10–3.14 on macOS arm64, Linux x86_64
and aarch64, and Windows x86_64. See the
[supported-platforms table](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md#supported-platforms)
for operating-system requirements and Docker support. Free-threaded CPython
and Intel macOS need a source build.

### Optional CUDA package

`smsd-cuda` provides optional NVIDIA CUDA screening on Linux x86_64 and needs
the same version of `smsd`:

```bash
python -m pip install "smsd==8.0.2" "smsd-cuda==8.0.2"
```

The GPU computes integer overlaps; scores use the CPU expressions and hits
keep target order. Matching runs on the CPU. See the
[CUDA section](https://github.com/asad/SMSD/blob/master/docs/PYTHON.md#optional-cuda-screening)
for dispatch and CPU fallback behaviour.

## Quick start

```python
import smsd

query = "c1ccccc1"       # benzene
target = "c1ccc(O)cc1"   # phenol

assert smsd.Substructure.find(query).match(target).found
result = smsd.MCS.find(query).with_timeout(1000).match(target)
print(f"MCS: {result.best.atom_count} atoms")  # MCS: 6 atoms
```

`Substructure` and `MCS` prepare a query once and provide `match`,
`match_all` and `match_many`. Each match maps query to target atom indices and
records the selected bonds and fragments. A search that stops early can return
a valid but unproved answer, so check `result.optimal`, `result.complete` and
`result.status`; an empty incomplete result does not prove absence.

`MCS` defaults to a connected maximum common induced subgraph (MCIS), ranked by
atoms and then bonds. `MCS.Objective.MCES` selects the maximum common edge
subgraph (MCES), which maximises common bonds, then their atoms, and adds no
isolated atoms. `with_fragments()` allows disconnected common graphs.

## Batch search

Parse molecules once when you reuse them. Results follow the target order, and
each target has its own search budget:

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
targets = [smsd.parse_smiles(s) for s in ["c1ccc(O)cc1", "CCO"]]
results = smsd.Substructure.find(query).match_many(targets, num_threads=2)
assert [result.found for result in results] == [True, False]
assert [result.best.atom_count if result.found else 0 for result in results] == [6, 0]
results = smsd.MCS.find(query).with_timeout(1000).match_many(targets)
assert [result.best.atom_count for result in results] == [6, 2]
```

`match_all` returns up to 1,000 matches by default, symmetry-related mappings
included; check `complete` before treating an enumeration as exhaustive.
`TargetCorpus` holds targets for repeated queries. Batch calls use OpenMP, or
run sequentially without it.

For repeated screening, `TargetCorpus.screen_many` returns ECFP4 Tanimoto hit
indices for each query, and `screen_targets_many` returns similarity-bound hit
indices. See the
[corpus API](https://github.com/asad/SMSD/blob/master/docs/PYTHON.md#smarts-and-target-collections)
for preparation and return order.

## Fingerprints

Radius 2 gives ECFP4; `mode="fcfp"` gives functional-class fingerprints
([Rogers and Hahn, 2010](https://doi.org/10.1021/ci100050t)). `similarity()` is
a screening bound, not a fingerprint metric; use the coefficients below for
fingerprint similarity:

```python
import smsd

query_fp = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)
target_fp = smsd.fingerprint_from_smiles("c1ccc(O)cc1", radius=2, fp_size=2048)
score = smsd.tanimoto_coefficient(query_fp, target_fp)
assert 0.0 <= score <= 1.0
```

`FingerprintOptions` configures circular, path, torsion and structural-key
fingerprints, including the 166 MACCS keys, with the same bits as in C++ and
Java:

```python
import smsd

options = smsd.FingerprintOptions(family="ecfp", ring_size=True, cip_atoms=True)
fp = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=options)
keys = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=smsd.FingerprintOptions.maccs())
assert fp and keys
```

### Screening

Rank a library by ECFP4 similarity, then search it by substructure as a
prepared corpus:

```python
import smsd

library = ["CC(=O)Oc1ccccc1", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", "c1ccccc1"]
ranked = smsd.batch_similarity("CC(=O)Oc1ccccc1C(=O)O", library)  # (index, Tanimoto), most similar first
corpus = smsd.TargetCorpus.from_smiles(library)
corpus.prewarm()
acids = smsd.Substructure.smarts("C(=O)[OX2H1]").match_many(corpus, num_threads=0)
assert library[ranked[0][0]] == "CC(=O)Oc1ccccc1" and sum(r.found for r in acids) == 1
```

With OpenMP, as in the published wheels, `num_threads=0` uses every core
unless `OMP_NUM_THREADS` sets fewer. The
[Python examples](https://github.com/asad/SMSD/blob/master/docs/EXAMPLES.md#screening)
screen a library of 2,744 molecules and set aside targets that cannot reach an
MCS similarity threshold before computing the MCS; they also cover count
fingerprints and the other coefficients.

## RDKit molecules

With RDKit installed, pass its molecules directly. Mappings use the original
RDKit atom indices:

```python
from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("c1ccccc1")
target = Chem.MolFromSmiles("c1ccc(O)cc1")
result = smsd.MCS.find(query).with_timeout(1000).match(target)
assert result.best.atom_count == 6
```

Conversion passes through SMILES, so coordinates, data fields and query
predicates are not transferred. `from_openbabel` and `to_openbabel` do the same
for Open Babel, without its atom numbering. See
[input atom numbering](https://github.com/asad/SMSD/blob/master/docs/PYTHON.md#input-atom-numbering).

## Documentation

- [Python guide](https://github.com/asad/SMSD/blob/master/docs/PYTHON.md):
  searches, chemistry options, timeouts and the full API reference.
- [Python examples](https://github.com/asad/SMSD/blob/master/docs/EXAMPLES.md):
  fingerprints, drawing, layout, stereo, MOL and SDF files, R-groups and
  scaffolds.
- [Installation guide](https://github.com/asad/SMSD/blob/master/docs/HOWTO-INSTALL.md):
  platform requirements and source builds.
- [Benchmark guide](https://github.com/asad/SMSD/blob/master/benchmarks/README.md):
  workloads and mapping validation.
- [Changelog](https://github.com/asad/SMSD/blob/master/CHANGELOG.md): what the
  release provides.

## Build from source

From the repository root, with Python 3.10 or later, a C++17 compiler and CMake
3.18 or later:

```bash
python -m pip install build
python -m pip install -e ".[dev]"
python -m build
```

The standard wheels use CPU/OpenMP. The optional `smsd-cuda` package needs a
compatible NVIDIA driver; macOS uses the CPU. `gpu_device_info()` reports an
available backend, not the one a particular call used.

## Other languages

SMSD Pro is also available for
[Java](https://github.com/asad/SMSD/blob/master/docs/JAVA.md), with CDK
integration, and as a header-only
[C++17 library](https://github.com/asad/SMSD/blob/master/docs/CPP.md).

## Citation

If you use SMSD Pro in your research, please cite:

> Rahman SA.
> *SMSD Pro: Coverage-Driven, Tautomer-Aware Maximum Common Substructure Search.*
> ChemRxiv, 2026.
> DOI: [10.26434/chemrxiv.15001534/v1](https://doi.org/10.26434/chemrxiv.15001534/v1)

For the original SMSD toolkit, please also cite:

> Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
> *Small Molecule Subgraph Detector (SMSD) toolkit.*
> Journal of Cheminformatics, 1:12, 2009.
> DOI: [10.1186/1758-2946-1-12](https://doi.org/10.1186/1758-2946-1-12)

A machine-readable [CITATION.cff](https://github.com/asad/SMSD/blob/master/CITATION.cff)
is provided for citation software.

## Licence

Apache License 2.0. Copyright (c) 2018-2026 BioInception PVT LTD. Algorithm
Copyright (c) 2009-2026 Syed Asad Rahman. See
[LICENSE](https://github.com/asad/SMSD/blob/master/LICENSE) and
[NOTICE](https://github.com/asad/SMSD/blob/master/NOTICE).

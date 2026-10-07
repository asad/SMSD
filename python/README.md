<p align="center">
  <a href="https://github.com/asad/SMSD" aria-label="SMSD">
    <img src="https://raw.githubusercontent.com/asad/SMSD/master/icons/icon.svg" alt="SMSD" width="140"/>
  </a>
</p>

# SMSD 7.2.2 for Python

[![PyPI](https://img.shields.io/pypi/v/smsd)](https://pypi.org/project/smsd/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/asad/SMSD/blob/master/LICENSE)

Substructure search, maximum common substructure (MCS), fingerprints and
similarity screening. Core matching does not require RDKit, CDK or Java.
RDKit is optional for molecule conversion and drawing.

## Install

The 7.2.2 wheels use CPython 3.14 and include CPU/OpenMP support:

| Platform | Architecture | Requirement |
|---|---|---|
| Windows | x86_64 / AMD64 | Windows 10 or later |
| Linux | x86_64 | glibc 2.28+ |
| macOS | arm64 / Apple Silicon | macOS 26+ |

Install 7.2.2 from [PyPI](https://pypi.org/project/smsd/) when listed:

```bash
python -m pip install smsd==7.2.2
```

You can also download and install the matching wheel from the
[GitHub release](https://github.com/asad/SMSD/releases/tag/v7.2.2).
Intel macOS, Linux arm64 and other Python versions require a source build.

## Quick start

```python
import smsd

query = "c1ccccc1"       # benzene
target = "c1ccc(O)cc1"   # phenol

assert smsd.is_substructure(query, target)
mapping = smsd.find_mcs(query, target, timeout_ms=1000)
print(f"MCS: {len(mapping)} atoms")  # MCS: 6 atoms
```

Mappings link query atom indices to target atom indices. An empty mapping means
no match was found. A search timeout can leave a smaller MCS.

## Fingerprints

Radius 2 gives ECFP4; use `mode="fcfp"` for functional-class fingerprints
(Rogers and Hahn, 2010).
`similarity()` is a screening upper bound. Use fingerprint metrics for fingerprint
similarity:

```python
import smsd

query_fp = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)
target_fp = smsd.fingerprint_from_smiles("c1ccc(O)cc1", radius=2, fp_size=2048)
score = smsd.tanimoto_coefficient(query_fp, target_fp)
assert 0.0 <= score <= 1.0
```

For count fingerprints and other metrics, see the
[fingerprint examples](https://github.com/asad/SMSD/blob/master/docs/EXAMPLES.md#6-fingerprints-and-similarity).

## Batch search

Parse molecules once when reusing them. Batch results follow target order:

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
targets = [smsd.parse_smiles(s) for s in ["c1ccc(O)cc1", "CCO"]]
assert smsd.batch_substructure(query, targets, num_threads=2) == [True, False]
assert [len(m) for m in smsd.batch_find_substructure(query, targets)] == [6, 0]
sizes = smsd.batch_mcs_size(query, targets, timeout_ms=1000)
assert sizes == [6, 2]
```

Use `batch_mcs()` for mappings and `batch_mcs_size()` for atom counts.
`TargetCorpus` supports repeated queries against one collection.

## Using RDKit

Install RDKit separately to pass its molecules directly to SMSD. Returned
mappings use the original RDKit atom indices:

```python
from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("c1ccccc1")
target = Chem.MolFromSmiles("c1ccc(O)cc1")
mapping = smsd.find_mcs(query, target, timeout_ms=1000)
assert len(mapping) == 6
```

## More examples

The [Python guide](https://github.com/asad/SMSD/blob/master/docs/PYTHON.md) covers
chemistry options, SMARTS, stereo, tautomer matching, fingerprints, MOL/SDF I/O,
drawing and batch operations. See the
[examples](https://github.com/asad/SMSD/blob/master/docs/EXAMPLES.md) for complete workflows.

## Build from source

Run these commands at the repository root with a C++17 compiler and CMake
3.18 or later. Source metadata allows Python 3.9 or later:

```bash
python -m pip install build
python -m pip install -e ".[dev]"
python -m build
```

Release wheels use CPU/OpenMP. Metal and CUDA are optional source-build features
that need compatible tools and hardware. `gpu_device_info()` reports the active
backend; batch matching uses the CPU.

## Tests and benchmarks

Version 7.2.2 passes 691 Python tests with 8 optional skips per platform; see the
[test report](https://github.com/asad/SMSD/blob/master/docs/VALIDATION_7.2.2.md).
The [7.2.0 benchmark report](https://github.com/asad/SMSD/blob/master/benchmarks/RESULTS_7.2.0.md)
contains measured comparisons for its recorded versions and molecules.

## Other languages

[Java](https://github.com/asad/SMSD/tree/master/java) and
[C++](https://github.com/asad/SMSD/tree/master/cpp) are also available.
Java 7.2.2 is on GitHub; Maven Central remains at `com.bioinceptionlabs:smsd:7.1.1`
until 7.2.2 is published there.

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

A machine-readable [CITATION.cff](https://github.com/asad/SMSD/blob/master/CITATION.cff) is available for automated citation tools.

## Licence

Apache 2.0 — Copyright (c) 2018-2026 Syed Asad Rahman, BioInception PVT LTD.
See [LICENSE](https://github.com/asad/SMSD/blob/master/LICENSE) and [NOTICE](https://github.com/asad/SMSD/blob/master/NOTICE) for licensing and attribution.

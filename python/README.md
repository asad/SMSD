<p align="center">
  <a href="https://github.com/asad/SMSD" aria-label="SMSD">
    <img src="https://raw.githubusercontent.com/asad/SMSD/master/icons/icon.svg" alt="SMSD" width="140"/>
  </a>
</p>

# SMSD for Python

The examples use the current source API. Release 7.2.2 predates the
`MCS`/`Substructure` API change; build from source for these examples.

[![PyPI](https://img.shields.io/pypi/v/smsd)](https://pypi.org/project/smsd/)
[![Licence](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/asad/SMSD/blob/master/LICENSE)

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

Install 7.2.2 from [PyPI](https://pypi.org/project/smsd/):

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

assert smsd.Substructure.find(query).match(target).found
result = smsd.MCS.find(query).with_timeout(1000).match(target)
print(f"MCS: {result.best.atom_count} atoms")  # MCS: 6 atoms
```

Matches link query atom indices to target atom indices and retain selected
bonds and fragments. Check `result.optimal`, `result.complete` and
`result.status` before treating an answer as a completed search.

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
results = smsd.Substructure.find(query).match_many(targets, num_threads=2)
assert [result.found for result in results] == [True, False]
assert [result.best.atom_count if result.found else 0 for result in results] == [6, 0]
results = smsd.MCS.find(query).with_timeout(1000).match_many(targets)
assert [result.best.atom_count for result in results] == [6, 2]
```

Both queries expose `match`, `match_all` and `match_many`. `MCS` supports
MCIS and MCES objectives and disconnected fragment mapping. `TargetCorpus`
supports repeated queries against one collection.

## Using RDKit

Install RDKit separately to pass its molecules directly to SMSD. Returned
mappings use the original RDKit atom indices:

```python
from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("c1ccccc1")
target = Chem.MolFromSmiles("c1ccc(O)cc1")
result = smsd.MCS.find(query).with_timeout(1000).match(target)
assert result.best.atom_count == 6
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
For the current prepared-query API, use the shared workloads and mapping
validation in the [benchmark guide](https://github.com/asad/SMSD/blob/master/benchmarks/README.md).

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

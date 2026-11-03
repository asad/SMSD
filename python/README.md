# SMSD Pro™ for Python

[![PyPI](https://img.shields.io/pypi/v/smsd?label=PyPI)](https://pypi.org/project/smsd/)
[![Conda](https://img.shields.io/conda/vn/bioinception/smsd?label=conda)](https://anaconda.org/bioinception/smsd)
[![Python](https://img.shields.io/badge/Python-3.12%E2%80%933.15-3776ab)](https://asad.github.io/SMSD/HOWTO-INSTALL/)
[![Documentation](https://img.shields.io/badge/docs-user%20guide-0d3330)](https://asad.github.io/SMSD/guide/)
[![Licence](https://img.shields.io/badge/licence-Apache%202.0-blue)](https://github.com/asad/SMSD/blob/master/LICENSE)
[![Software DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23639199-007ec6)](https://doi.org/10.5281/zenodo.23639199)

Find chemical motifs, compare molecular structures and screen collections
with substructure search, maximum common substructure (MCS) and molecular
fingerprints. The Python package uses SMSD's native C++ engine. The examples
below need only `smsd`.

**[User guide](https://asad.github.io/SMSD/guide/)** ·
[Getting started](https://asad.github.io/SMSD/guide/GETTING_STARTED/) ·
[Notebooks](https://asad.github.io/SMSD/notebooks/) ·
[Source](https://github.com/asad/SMSD) ·
[Report a bug](https://github.com/asad/SMSD/issues)

## Install

### pip

```bash
python -m pip install "smsd==8.2.1"
```

### conda or mamba

```bash
conda install -c bioinception -c conda-forge smsd=8.2.1
```

Use the same channels and package name with `mamba`.

Python wheels are available for CPython 3.12–3.15:

| Operating system | Wheel architecture |
|---|---|
| Windows | x86_64 |
| Linux | x86_64, aarch64 |
| macOS 14 or later | arm64 (Apple silicon) |

See the [installation guide](https://asad.github.io/SMSD/HOWTO-INSTALL/)
for Linux and Windows requirements, conda packages and source builds. Intel
macOS and free-threaded CPython require a source build. Batch calls use OpenMP when available and otherwise run
sequentially.

### Optional GPU screening

Install the add-on for your hardware, keeping its version equal to `smsd`:

| Hardware | pip command |
|---|---|
| NVIDIA GPU, Linux x86_64 | `python -m pip install "smsd==8.2.1" "smsd-cuda==8.2.1"` |
| Apple silicon GPU, macOS 14+ | `python -m pip install "smsd==8.2.1" "smsd-metal==8.2.1"` |

These add-ons run supported screening operations on the GPU. Exact graph
matching runs on the CPU. The [CPU and GPU guide](https://asad.github.io/SMSD/guide/CPU_GPU/)
explains requirements and fallback behaviour; the
[Windows GPU guide](https://asad.github.io/SMSD/guide/GPU_ON_WINDOWS/)
covers NVIDIA screening through WSL2.

## Find a chemical motif

Search for a carboxylic acid with a SMARTS pattern. Batch results follow the
input order:

```python
import smsd

library = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ethanol": "CCO",
}
acid = smsd.Substructure.smarts("C(=O)[OX2H1]")
results = acid.match_many(list(library.values()), num_threads=1)

for name, hit in zip(library, results):
    print(name, hit.found)
# aspirin True
# paracetamol False
# ethanol False
```

For molecule queries, use `Substructure.find(query)`. See
[substructure search](https://asad.github.io/SMSD/guide/SUBSTRUCTURE/) and
[SMARTS patterns](https://asad.github.io/SMSD/guide/SMARTS/) for chemistry options.

## Compare two molecules

Find the common carbon skeleton of ethanol and ethylamine:

```python
import smsd

query, target = "CCO", "CCN"
result = smsd.MCS.find(query).with_timeout(1000).match(target)

if result.found:
    print(result.best.atom_count, result.best.bond_count)  # 2 1
    print(result.best.atom_mapping)  # query atom index -> target atom index
print(result.status.name, result.optimal, result.complete)
```

The timeout is in milliseconds per target. A stopped search may return a
valid match without proving it is the largest. Check `optimal`, `complete`
and `status`; an empty unfinished result does not prove that no match exists.
The default is a connected maximum common induced subgraph (MCIS), ranked by
atoms, then bonds. The [MCS guide](https://asad.github.io/SMSD/guide/MCS/)
covers edge-based MCES, disconnected fragments and matching rules.

### Draw the match

Draw the match as a native SVG. `mcs_svg` reuses the result's mapping, so it
runs no second search and needs no optional drawing package:

```python
from pathlib import Path

import smsd
from smsd.highlight import mcs_svg

query, target = "CCO", "CCN"
result = smsd.MCS.find(query).with_timeout(1000).match(target)
svg = mcs_svg(query, target, result)
Path("common-substructure.svg").write_text(svg, encoding="utf-8")
```

The SVG retains search status and highlights the mapped region. See
[drawing a result](https://asad.github.io/SMSD/guide/MCS/#draw-the-result)
for notebook display and selected-bond rendering limits.

## Configure a search

Set the time and node limits and the chemistry rules before preparing a query:

```python
import smsd

options = smsd.MCS.Options()
options.timeout_ms = 1000
options.node_limit = 100000
options.chemistry.match_bond_order = smsd.BondOrderMode.STRICT
options.chemistry.ring_matches_ring_only = True
options.chemistry.use_chirality = True
options.chemistry.use_bond_stereo = True

prepared = smsd.MCS.find("c1ccc(O)cc1").with_options(options)
configured = prepared.match("c1ccc(N)cc1")
atoms = configured.best.atom_count if configured.found else 0
print(atoms, configured.status.name, configured.optimal, configured.complete)  # 6 COMPLETED True True
```

| Setting | Controls |
|---|---|
| `timeout_ms`, `node_limit` | Time per target and search-node allowance |
| `chemistry.match_bond_order` | Equal bond orders (`STRICT`) or relaxed matching |
| `chemistry.aromaticity_mode` | Flexible or strict aromaticity matching |
| `chemistry.ring_matches_ring_only` | Ring atoms and bonds map only to rings |
| `chemistry.complete_rings_only` | Matched ring atoms and bonds lie in fully matched rings of both molecules |
| `chemistry.use_chirality`, `chemistry.use_bond_stereo` | Tagged atom and double-bond stereochemistry |
| `chemistry.match_formal_charge`, `chemistry.match_isotope` | Charge agreement; isotope agreement where both atoms specify a mass |
| `chemistry.tautomer_aware` | Recognised tautomeric relaxations; usually also needs a relaxed `match_bond_order` |
| `.with_objective(smsd.MCS.Objective.MCES)` | Maximise common bonds instead of common atoms |
| `.with_fragments()` | Allow a disconnected common graph |

`with_options` copies the settings into a new prepared query. Molecule-based
substructure searches use `Substructure.Options()` in the same way. SMARTS
patterns state their own atom and bond rules; among these chemistry settings,
only `complete_rings_only` applies to SMARTS. See the
[full configuration table](https://asad.github.io/SMSD/guide/SUBSTRUCTURE/#chemistry-options).

## Compare molecular fingerprints

Compute ECFP4 fingerprints and their Tanimoto similarity:

```python
import smsd

benzene = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)
phenol = smsd.fingerprint_from_smiles("c1ccc(O)cc1", radius=2, fp_size=2048)
score = smsd.tanimoto_coefficient(benzene, phenol)
print(f"Tanimoto: {score:.3f}")
```

The [fingerprint guide](https://asad.github.io/SMSD/guide/FINGERPRINTS/)
covers families and options; the [screening guide](https://asad.github.io/SMSD/guide/SCREENING/)
covers searching a collection. `smsd.similarity()` is a separate screening
bound; use `tanimoto_coefficient()` for the fingerprint comparison above.

## Screen a library, then compute the MCS

Rank a library by fingerprint similarity, or prepare a corpus for repeated
MCS-bound screening before graph search:

```python
import smsd

query = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
library = [query, "OC(=O)c1ccccc1O", "CCO"]
corpus = smsd.TargetCorpus.from_smiles(library)

for index, score in smsd.batch_similarity(query, library):
    print(index, f"ECFP4 Tanimoto: {score:.3f}")

candidates = smsd.screen_targets(query, corpus, 0.7)
matches = smsd.MCS.find(query).with_timeout(1000).match_many(
    [library[index] for index in candidates], num_threads=1
)
for index, match in zip(candidates, matches):
    atoms = match.best.atom_count if match.found else 0
    print(index, atoms, match.status.name, match.optimal, match.complete)
```

The `0.7` threshold applies to an **upper bound on atom-count MCS Tanimoto**,
under the same chemistry rules as the default search above. It is separate
from ECFP4 Tanimoto and is not a confirmed MCS score. A surviving candidate
still needs graph search and inspection of its result status.
The [screening examples](https://github.com/asad/SMSD/blob/master/docs/EXAMPLES.md#screening)
extend this workflow to larger libraries and substructure confirmation.

## Guides and examples

| Task | Documentation |
|---|---|
| Work through complete examples | [Jupyter and Colab notebooks](https://asad.github.io/SMSD/notebooks/) |
| Reuse a target library | [Corpus screening](https://asad.github.io/SMSD/guide/SCREENING/) |
| Store fingerprints and select top-k hits | [Packed fingerprint index](https://asad.github.io/SMSD/guide/PACKED_INDEX/) |
| Cluster molecules or select a diverse subset | [Clustering and diversity](https://asad.github.io/SMSD/guide/CLUSTERING/) |
| Use data frames, arrays and learning libraries | [Data adapters](https://asad.github.io/SMSD/guide/MACHINE_LEARNING/) |
| Look up Python methods and options | [Python reference](https://asad.github.io/SMSD/PYTHON/) |
| Use C++ or Java | [Other languages](https://asad.github.io/SMSD/guide/OTHER_LANGUAGES/) |

The [installation guide](https://asad.github.io/SMSD/HOWTO-INSTALL/) also
covers optional integrations and source-build dependencies. See the
[release notes](https://asad.github.io/SMSD/RELEASE_NOTES/) when upgrading.
Existing deprecated APIs remain supported throughout the 8.x line.

SMSD is maintained at the 8.2.0 feature baseline across C++, Python and Java.
Updates cover bug fixes, security, compatibility and documentation corrections.

## Cite SMSD Pro

If you use SMSD Pro in research, cite the paper and software release:

> Rahman SA. *SMSD Pro: Coverage-Driven, Tautomer-Aware Maximum Common
> Substructure Search.* ChemRxiv, 2026.
> [doi:10.26434/chemrxiv.15001534/v1](https://doi.org/10.26434/chemrxiv.15001534/v1)

Software: [doi:10.5281/zenodo.23639199](https://doi.org/10.5281/zenodo.23639199).
This DOI covers all releases; cite the version you used.

For the original toolkit:

> Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
> *Small Molecule Subgraph Detector (SMSD) toolkit.*
> Journal of Cheminformatics, 1:12, 2009.
> [doi:10.1186/1758-2946-1-12](https://doi.org/10.1186/1758-2946-1-12)

[CITATION.cff](https://github.com/asad/SMSD/blob/master/CITATION.cff)
provides machine-readable citation details.

## Licence

Apache License 2.0. Copyright (c) 2018-2026 BioInception PVT LTD. Algorithm
Copyright (c) 2009-2026 Syed Asad Rahman. See
[LICENSE](https://github.com/asad/SMSD/blob/master/LICENSE) and
[NOTICE](https://github.com/asad/SMSD/blob/master/NOTICE).

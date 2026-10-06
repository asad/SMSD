<p align="center">
  <a href="https://github.com/asad/SMSD" aria-label="SMSD Pro">
    <img src="https://raw.githubusercontent.com/asad/SMSD/master/icons/icon.svg" alt="SMSD Pro" width="140"/>
  </a>
</p>

# SMSD Pro for Python

[![PyPI](https://img.shields.io/pypi/v/smsd)](https://pypi.org/project/smsd/)
[![Downloads](https://img.shields.io/pypi/dm/smsd)](https://pypi.org/project/smsd/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/asad/SMSD/blob/master/LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/smsd)](https://pypi.org/project/smsd/)

Python bindings for SMSD native graph matching, including substructure search,
maximum common substructure (MCS), fingerprints, and molecular similarity.
RDKit and CDK are not required for the core SMSD path.

Performance depends on corpus, chemistry constraints, search budget and
result validity. See the [current local report](https://github.com/asad/SMSD/blob/master/benchmarks/RESULTS_7.2.0.md)
for the 7.1.2 baseline, 7.2.0 source candidate and RDKit 2026.09.1 comparison.
The 7.2.1 deadline regression is recorded separately. Full corpus comparisons
have not been rerun for its patched source, layout and packaging changes.
The MoleculeNet-derived Dalke-style pairs are not the original Dalke benchmark.

## Install

```bash
pip install smsd
```

Build from source from the repository root, which contains the C++ sources
and the canonical package metadata. `python/smsd/` contains the Python layer,
`cpp/` contains its native extension, and the root `pyproject.toml` builds both.
There is no separate Python manifest under `python/`:

```bash
python -m pip install -e ".[dev]"
python -m build
```

The package declares CPython `3.9` or later; wheel availability depends on
platform and architecture. The proposed 7.2.1 release targets Python 3.14 wheels
for Linux x86_64 (glibc 2.28+), macOS arm64 (26+) and Windows x86_64, plus a
source distribution. Each wheel requires an installed-package test on its
target operating system before publication. Intel macOS and Linux arm64
remain source-build targets. The controlled search review uses Python
`3.13.14` on macOS arm64.
Source builds default to Metal/CUDA auto-detection; release and comparison wheels disable
both explicitly. Core batch matching uses CPU/OpenMP. RDKit
remains optional for interop and depiction rather than a core dependency.

## Quick Start

```python
import smsd

# Substructure search
assert smsd.is_substructure(
    smsd.parse_smiles("c1ccccc1"),     # benzene
    smsd.parse_smiles("c1ccc(O)cc1"))  # phenol

# Maximum Common Substructure
mcs = smsd.find_mcs("c1ccccc1", "c1ccc2ccccc2c1")
print(f"MCS: {len(mcs)} atoms")  # 6

# Tautomer-aware MCS
mcs = smsd.find_mcs("CC(=O)C", "CC(O)=C", tautomer_aware=True)

# Circular fingerprint (ECFP4)
ecfp4 = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)

# Similarity
sim = smsd.similarity("c1ccccc1", "c1ccc(O)cc1")
```

## Features

### Search & Matching
| Feature | Description |
|---------|-------------|
| **Substructure search** | Native engine with 3-level NLF pruning, GPU-accelerated domain init |
| **MCS** | Multi-strategy native pipeline (chain/tree fast paths + general clique solver) |
| **SMARTS matching** | Connectivity, degree, valence, rings, stereo, recursive patterns and logical operators; dialect differences are recorded in the benchmark report |
| **Tautomer matching** | 30 transforms with pKa-informed weights, 6 solvents, pH-sensitive |
| **CIP R/S/E/Z** | `assign_rs()`, `assign_ez()` — native stereo descriptor assignment; see validation scope |
| **MCS SMILES** | `find_mcs_smiles()` — extract MCS as canonical SMILES string |
| **Multi-result MCS** | `find_mcs(mol1, mol2, max_results=N)` — top-N MCS enumeration |
| **SMARTS MCS** | `find_mcs_smarts()` — largest substructure matching a SMARTS pattern |
| **R-group decomposition** | `decompose_r_groups()` — scaffold + R-group extraction |

### Fingerprints
| Type | Description |
|------|-------------|
| **Circular ECFP** | Tautomer-aware structural invariants, configurable radius (2=ECFP4, 3=ECFP6, -1=whole molecule) |
| **Circular FCFP** | Pharmacophoric invariants (H-bond donor/acceptor, ionisable, aromatic, hydrophobic) |
| **Count-based ECFP/FCFP** | `circular_fingerprint_counts()` / `fcfp_counts()` — retain feature multiplicities |
| **Topological Torsion** | `topological_torsion()` — 4-atom path fingerprint |
| **Path fingerprint** | Graph-aware DFS path enumeration, tautomer-invariant |
| **MCS fingerprint** | MCS-aware, uses chemical matching rules for path compatibility |
| **Similarity metrics** | `overlap_coefficient()`, `tanimoto_coefficient()`, `dice()`, `cosine()`, `soergel()` — binary + count-vector |
| **Format conversions** | `to_hex()`, `to_binary_string()`, `from_hex()` — for database storage and REST APIs |
| **Subset check** | `fingerprint_subset()` — fast substructure pre-screening |

### Infrastructure
| Feature | Description |
|---------|-------------|
| **MatchResult** | `mcs_from_smiles()` — structured result: size, mapping, overlap_coefficient |
| **RDKit interop** | `mcs_rdkit_native()`, `batch_mcs_rdkit()`, `from_rdkit()` with correct indices |
| **Similarity screening** | RASCAL O(V+E) upper bound for fast pre-filtering |
| **Lenient parser** | Best-effort recovery from malformed SMILES |
| **Batch operations** | OpenMP-parallel `batch_substructure()`, `batch_mcs()`, `batch_mcs_rdkit()` |
| **Adaptive timeout** | `min(30s, 500+n1*n2*2)` based on molecule size |
| **GPU acceleration** | CUDA + Apple Metal for domain init and RASCAL screening |
| **Force-directed layout** | `force_directed_layout()` for bond-crossing minimisation |
| **SMACOF stress majorisation** | `stress_majorisation()` to reduce embedding stress |
| **Scaffold templates** | `match_template()` for 10 pre-computed common scaffolds |
| **Ring perception** | `compute_sssr()`, `layout_sssr()` — clean SSSR APIs |

## Performance and validation

The [benchmark report](https://github.com/asad/SMSD/blob/master/benchmarks/RESULTS_7.2.0.md)
records timings alongside atom/bond counts, mapping validity, budgets and
RDKit cancellation flags. Differences in MCS semantics can make raw size or
speed comparisons misleading. No universal speed or quality claim is made.

For repeated work, reuse parsed `MolGraph` objects and use `batch_mcs_size`
when only counts are needed. Core batch bindings retain graph references
instead of copying their caches into a temporary vector. Measurements of
conversion, cache hits and batch overhead are reported separately.

On the controlled local workload, cached RDKit conversion changed from
4.36 to 1.96 microseconds, a 32-target substructure batch from 78.9 to 5.0
microseconds, and compiled SMARTS matching on 32 repeated 448-atom targets
from 5,807 to 78.4 microseconds. Small native MCS dispatch increased from
17.1 to 21.4 microseconds. These checksum-matched measurements exclude setup;
they describe this workload rather than an application-wide speedup.

The historical 7.2.0 installed CPU wheels each passed **691 Python tests** with **8 skips**
on macOS arm64: Python 3.13.14 with RDKit 2026.09.1, and Python 3.14.8 with
the published RDKit 2026.03.6 wheel. The 3.14 release wheel bundles OpenMP
and targets macOS 26 or later. See
[7.2.0 validation](https://github.com/asad/SMSD/blob/master/docs/VALIDATION_7.2.0.md)
for scope and reproduction commands. Its Linux x86_64 wheel also passed 691
tests with 8 skips on CPython 3.14.5/RDKit 2026.03.6 under local emulation.
The corrected 7.2.0 Windows x86_64 wheel passed the same Python test counts
on Windows Server 2022 with CPython 3.14.7/RDKit 2026.03.6. All 12 native
Debug suites also passed
on Windows with MSVC. These checks do not extend the macOS benchmark timings
to other platforms.
These 7.2.0 results are separate from 7.2.1 validation.
The frozen-source 7.2.1 macOS arm64 wheel passed all 12 native Debug suites and
691 Python tests with 8 skips on Python 3.14.8/RDKit 2026.03.6, with bundled
OpenMP. Execution was on macOS 27.0.1; the wheel targets macOS 26+, without
a claim of testing the minimum OS. The same-source Linux x86_64 wheel also
passed all 12 native suites and 691 Python tests with 8 skips on Python
3.14.5/RDKit 2026.03.6, using glibc 2.28 under local emulation with bundled
OpenMP. The Windows Server 2022/AMD64 wheel passed the same test counts on
Python 3.14.7/RDKit 2026.03.6 with active OpenMP and checked Microsoft runtimes.
Strict collection of all three wheels passed against one source archive; see
[7.2.1 validation](https://github.com/asad/SMSD/blob/master/docs/VALIDATION_7.2.1.md).
Publication remains pending.

## Circular Fingerprints

Tautomer-aware Morgan/ECFP — includes tautomer class in the atom invariant,
with an optional tautomer-class invariant. Similarity depends on the selected
features and molecules.

**ECFP vs FCFP:** SMSD supports **both** fingerprint types (Rogers & Hahn 2010):

- **ECFP** (Extended Connectivity): atom invariant = atomic number, degree, charge, ring, aromaticity, tautomer class. Encodes atom environments.
- **FCFP** (Functional Class): atom invariant = pharmacophoric features (H-bond donor/acceptor, positive/negative ionisable, aromatic, hydrophobic). Encodes feature environments.

| Name | Radius | Type | SMSD call |
|------|--------|------|-----------|
| ECFP2 | 1 | Structural | `circular_fingerprint(mol, radius=1)` |
| ECFP4 | 2 | Structural | `circular_fingerprint(mol, radius=2)` |
| ECFP6 | 3 | Structural | `circular_fingerprint(mol, radius=3)` |
| FCFP2 | 1 | Pharmacophoric | `circular_fingerprint(mol, radius=1, mode="fcfp")` |
| FCFP4 | 2 | Pharmacophoric | `circular_fingerprint(mol, radius=2, mode="fcfp")` |
| FCFP6 | 3 | Pharmacophoric | `circular_fingerprint(mol, radius=3, mode="fcfp")` |
| Whole | -1 | Either | `circular_fingerprint(mol, radius=-1)` |

```python
import smsd

# ECFP4 (structural, recommended default)
ecfp4 = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)

# FCFP4 (pharmacophoric — H-bond donors/acceptors, ionisable, aromatic, hydrophobic)
fcfp4 = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048, mode="fcfp")

# ECFP6 (radius 3, captures larger environments)
ecfp6 = smsd.fingerprint_from_smiles("c1ccccc1", radius=3, fp_size=2048)

# ECFP2 (radius 1, smaller neighborhoods)
ecfp2 = smsd.fingerprint_from_smiles("c1ccccc1", radius=1, fp_size=2048)

# Whole molecule (radius -1 = expand until convergence)
whole = smsd.fingerprint_from_smiles("c1ccccc1", radius=-1, fp_size=2048)

# Tanimoto similarity (works with any fingerprint type)
sim = smsd.overlap_coefficient(
    smsd.fingerprint_from_smiles("c1ccccc1", radius=2),
    smsd.fingerprint_from_smiles("c1ccc(O)cc1", radius=2))
```

## Using with RDKit

SMSD works standalone or alongside RDKit. Use RDKit for parsing and drawing,
SMSD for native graph matching:

```python
from rdkit import Chem
import smsd

mol1 = Chem.MolFromSmiles("c1ccccc1")
mol2 = Chem.MolFromSmiles("c1ccc(O)cc1")

# MCS with RDKit molecules
result = smsd.mcs_rdkit(mol1, mol2)

# Depict MCS with highlighted atoms (works in Jupyter)
img = smsd.depict_mcs("c1ccccc1", "c1ccc(O)cc1")

# Export to SDF with the native writer
smsd.export_sdf([smsd.parse_smiles(s) for s in ["CCO", "c1ccccc1"]], "output.sdf")

# Convert between SMSD and RDKit
g = smsd.from_rdkit(mol1)
rdmol = smsd.to_rdkit(g)
```

> RDKit is optional — `pip install smsd` works without it.

## Solvent-Aware Tautomer Matching

```python
import smsd

mcs = smsd.find_mcs("CC(=O)C", "CC(O)=C", tautomer_aware=True)
```

Supported solvents: `AQUEOUS`, `DMSO`, `METHANOL`, `CHLOROFORM`, `ACETONITRILE`, `DIETHYL_ETHER`

## Native MOL/SDF I/O

```python
import smsd

g = smsd.parse_smiles("CCO")  # or smsd.read_molfile("input.mol")
mol_block = smsd.write_mol_block(g)
mol_block_v3000 = smsd.write_mol_block_v3000(g)
sdf_record = smsd.write_sdf_record(g)

smsd.write_molfile(g, "out_v2000.mol")
smsd.write_molfile(g, "out_v3000.mol", v3000=True)
smsd.write_molfile(g, "out.sdf", sdf=True)
```

The native writer preserves practical chemistry metadata in `7.1.1`:
- names, comments, and SDF properties
- charges, isotopes, atom classes, and atom maps
- `R#`/`R<n>` plus `M  RGP`
- practical V2000/V3000 stereo round-trip

## Platform & GPU Support

SMSD automatically dispatches to the best available compute backend:

| Platform | CPU | GPU |
|----------|-----|-----|
| macOS (Apple Silicon) | OpenMP | Metal (shared buffers) |
| macOS (Intel) | OpenMP | CPU fallback |
| Linux | OpenMP | CUDA (if available) |
| Windows | OpenMP | CUDA (if available) |

```python
import smsd

# Check GPU availability
if smsd.gpu_is_available():
    print(smsd.gpu_device_info())
    # e.g. "Metal GPU: Apple M2 Pro [OpenMP 5.0, 10 threads]"
    # e.g. "GPU: Tesla T4 [OpenMP 4.5, 8 threads]"

# Some screening operations can use an available compiled GPU backend.
# Core batch matching uses OpenMP on the CPU.
results = smsd.batch_substructure(query, targets, num_threads=2)
```

## API Reference

### Core Functions

```python
# Parsing
mol = smsd.parse_smiles("c1ccccc1")
smi = smsd.to_smiles(mol)

# Substructure
smsd.is_substructure(query, target)
mapping = smsd.find_substructure(query, target)

# MCS
mapping = smsd.find_mcs(mol1, mol2)
mapping = smsd.find_mcs("SMILES1", "SMILES2", tautomer_aware=True)

# Similarity
sim = smsd.similarity("SMILES1", "SMILES2")
ub = smsd.similarity_upper_bound(mol1, mol2)
hits = smsd.screen_targets(query, library, threshold=0.5)

# Fingerprints
fp = smsd.path_fingerprint(mol, path_length=7, fp_size=2048)
fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)           # ECFP4
fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048, mode="fcfp")  # FCFP4
sim = smsd.overlap_coefficient(fp1, fp2)
ok = smsd.fingerprint_subset(query_fp, target_fp)

# Format conversions (database storage, REST APIs)
hex_str = smsd.to_hex(fp, fp_size=2048)
fp_back = smsd.from_hex(hex_str)
bits    = smsd.to_binary_string(fp, fp_size=2048)

# Batch (OpenMP parallel)
results = smsd.batch_substructure(query, targets)
results = smsd.batch_mcs(query, targets)

# RASCAL pre-screen + exact MCS in one call
matches = smsd.screen_and_match(query, targets, threshold=0.5)

# Batch find substructure with atom-atom mappings (v7.1.1)
mappings = smsd.batch_find_substructure(query, targets)

# TargetCorpus — prewarm once, query many times (v7.1.1)
corpus = smsd.TargetCorpus.from_smiles(["c1ccccc1", "c1ccc(O)cc1", "CCO"])
corpus.prewarm()
hits = corpus.substructure(smsd.parse_smiles("c1ccccc1"))
sizes = corpus.mcs_size(smsd.parse_smiles("c1ccccc1"))
passing = corpus.screen(smsd.parse_smiles("c1ccccc1"), threshold=0.5)

# GPU
smsd.gpu_is_available()
smsd.gpu_device_info()
```

### Configuration

```python
opts = smsd.ChemOptions()
opts.match_atom_type = True
opts.tautomer_aware = True
opts.ring_fusion_mode = smsd.RingFusionMode.STRICT
opts.match_bond_order = smsd.BondOrderMode.LOOSE

# Profiles
opts = smsd.ChemOptions.tautomer_profile()
opts = smsd.ChemOptions.profile("strict")
```

## Scope and interoperability

SMSD provides molecular graph search, fingerprints, basic I/O and depiction.
RDKit integration is optional. Choose tools using your own molecules,
chemistry settings and validation requirements; feature lists alone do not
establish comparative accuracy or performance.

## Also Available

- **Java**: `com.bioinceptionlabs:smsd:7.1.1` on [Maven Central](https://central.sonatype.com/artifact/com.bioinceptionlabs/smsd)
- **C++**: Header-only, zero dependencies — [GitHub](https://github.com/asad/SMSD)

## Citation

If you use SMSD Pro in your research, please cite:

> Rahman SA.
> *SMSD Pro: Tautomer-Aware Maximum Common Substructure Search.*
> ChemRxiv, 2025.
> DOI: [10.26434/chemrxiv.15001534](https://doi.org/10.26434/chemrxiv.15001534/v1)

For the original SMSD toolkit, please also cite:

> Rahman SA, Bashton M, Holliday GL, Schrader R, Thornton JM.
> *Small Molecule Subgraph Detector (SMSD) toolkit.*
> Journal of Cheminformatics, 1:12, 2009.
> DOI: [10.1186/1758-2946-1-12](https://doi.org/10.1186/1758-2946-1-12)

A machine-readable [CITATION.cff](https://github.com/asad/SMSD/blob/master/CITATION.cff) is available for automated citation tools.

## License

Apache 2.0 — Copyright (c) 2018-2026 Syed Asad Rahman, BioInception PVT LTD.
See [NOTICE](https://github.com/asad/SMSD/blob/master/NOTICE) for attribution, trademark, and novel algorithm terms.

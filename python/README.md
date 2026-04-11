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

## Install

```bash
pip install smsd
```

Supported CPython versions: `3.9` through the latest stable release series.
Current default test target: `Python 3.12`.
The native SMSD path is CPU-first with optional GPU acceleration. RDKit
remains optional for interop and depiction rather than a core dependency.

## Quick Start

```python
import smsd

# Substructure search
assert smsd.is_substructure(
    smsd.parse_smiles("c1ccccc1"),     # benzene
    smsd.parse_smiles("c1ccc(O)cc1"))  # phenol

# Maximum Common Substructure
mcs = smsd.mcs("c1ccccc1", "c1ccc2ccccc2c1")
print(f"MCS: {len(mcs)} atoms")  # 6

# Tautomer-aware MCS
mcs = smsd.mcs("CC(=O)C", "CC(O)=C", tautomer_aware=True)

# Circular fingerprint (ECFP4, tautomer-aware)
ecfp4 = smsd.circular_fingerprint("c1ccccc1", radius=2, fp_size=2048)

# Similarity
sim = smsd.similarity("c1ccccc1", "c1ccc(O)cc1")
```

## Features

### Search & Matching
| Feature | Description |
|---------|-------------|
| **Substructure search** | VF2++ with 3-level NLF pruning, GPU-accelerated domain init |
| **MCS** | 11-level funnel: chain/tree DP → greedy → McSplit → BK → McGregor |
| **SMARTS matching** | Full SMARTS support including `X` (total connectivity), `D` (degree), `v` (valence), `R` (ring count), `r` (ring size), `x` (ring connectivity), `/` `\` (E/Z stereo), `$()` (recursive), logical AND/OR/NOT |
| **Tautomer matching** | 30 transforms with pKa-informed weights, 6 solvents, pH-sensitive |
| **Tautomer validation** | `validate_tautomer_consistency()` — proton conservation check |
| **CIP R/S/E/Z** | `assign_rs()`, `assign_ez()` — full digraph-based stereo descriptors (IUPAC 2013) |
| **MCS SMILES** | `find_mcs_smiles()` — extract MCS as canonical SMILES string |
| **findAllMCS** | `find_all_mcs()` — top-N MCS enumeration with canonical dedup |
| **SMARTS MCS** | `find_mcs_smarts()` — largest substructure matching a SMARTS pattern |
| **R-group decomposition** | `decompose_r_groups()` — scaffold + R-group extraction |

### Fingerprints
| Type | Description |
|------|-------------|
| **Circular ECFP** | Tautomer-aware structural invariants, configurable radius (2=ECFP4, 3=ECFP6, -1=whole molecule) |
| **Circular FCFP** | Pharmacophoric invariants (H-bond donor/acceptor, ionisable, aromatic, hydrophobic) |
| **Count-based ECFP/FCFP** | `ecfp_counts()` / `fcfp_counts()` — superior to binary for ML |
| **Topological Torsion** | `topological_torsion()` — 4-atom path fingerprint (SOTA on peptides) |
| **Path fingerprint** | Graph-aware DFS path enumeration, tautomer-invariant |
| **MCS fingerprint** | MCS-aware, uses chemical matching rules for path compatibility |
| **Similarity metrics** | `overlapCoefficient()`, `dice()`, `cosine()`, `soergel()` — binary + count-vector |
| **Format conversions** | `to_hex()`, `to_binary_string()`, `from_hex()` — for database storage and REST APIs |
| **Subset check** | `fingerprint_subset()` — fast substructure pre-screening |

### Infrastructure
| Feature | Description |
|---------|-------------|
| **MatchResult** | `mcs_result()` — structured result: size, mapping, overlapCoefficient |
| **RDKit interop** | `mcs_rdkit_native()`, `batch_mcs_rdkit()`, `from_rdkit()` with correct indices |
| **Similarity screening** | RASCAL O(V+E) upper bound for fast pre-filtering |
| **Lenient parser** | Best-effort recovery from malformed SMILES |
| **Batch operations** | OpenMP-parallel `batch_substructure()`, `batch_mcs()`, `batch_mcs_rdkit()` |
| **Adaptive timeout** | `min(30s, 500+n1*n2*2)` based on molecule size |
| **GPU acceleration** | CUDA + Apple Metal for domain init and RASCAL screening |
| **Force-directed layout** | `force_directed_layout()` for bond-crossing minimisation |
| **SMACOF stress majorisation** | `stress_majorisation()` for optimal 2D embedding |
| **Scaffold templates** | `match_template()` for 10 pre-computed common scaffolds |
| **Reaction-aware MCS** | `find_mcs(..., reaction_aware=True)` or `map_reaction_aware()` post-filter |
| **Ring perception** | `compute_sssr()`, `layout_sssr()` — clean SSSR APIs |

## Performance

Benchmarked alongside RDKit 2025.09.2 on the same machine, same Python process.
Both toolkits excel at different tasks — use whichever fits your workflow, or both together.

| Pair | SMSD | RDKit | Notes |
|------|------|-------|-------|
| Morphine / Codeine | 79 us | 579 ms | Complex ring system |
| Coronene self-match | 6 us | 712 us | Symmetric PAH |
| Caffeine / Theophylline | 17 us | 373 us | N-methyl difference |
| PEG-12 / PEG-16 | 39 us | 2.1 ms | Linear polymer |

Full data: [benchmarks/results_python.tsv](https://github.com/asad/SMSD/blob/master/benchmarks/results_python.tsv)

## Circular Fingerprint (Novel)

Tautomer-aware Morgan/ECFP — includes tautomer class in the atom invariant,
so tautomeric forms of the same molecule produce more similar fingerprints.

**ECFP vs FCFP:** SMSD supports **both** fingerprint types (Rogers & Hahn 2010):

- **ECFP** (Extended Connectivity): atom invariant = atomic number, degree, charge, ring, aromaticity, tautomer class. Best for structural similarity.
- **FCFP** (Functional Class): atom invariant = pharmacophoric features (H-bond donor/acceptor, positive/negative ionisable, aromatic, hydrophobic). Best for activity-based similarity and SAR.

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
ecfp4 = smsd.circular_fingerprint("c1ccccc1", radius=2, fp_size=2048)

# FCFP4 (pharmacophoric — H-bond donors/acceptors, ionisable, aromatic, hydrophobic)
fcfp4 = smsd.circular_fingerprint("c1ccccc1", radius=2, fp_size=2048, mode="fcfp")

# ECFP6 (radius 3, captures larger environments)
ecfp6 = smsd.circular_fingerprint("c1ccccc1", radius=3, fp_size=2048)

# ECFP2 (radius 1, fastest, less discriminating)
ecfp2 = smsd.circular_fingerprint("c1ccccc1", radius=1, fp_size=2048)

# Whole molecule (radius -1 = expand until convergence)
whole = smsd.circular_fingerprint("c1ccccc1", radius=-1, fp_size=2048)

# Tanimoto similarity (works with any fingerprint type)
sim = smsd.overlapCoefficient(
    smsd.circular_fingerprint("c1ccccc1", radius=2),
    smsd.circular_fingerprint("c1ccc(O)cc1", radius=2))
```

## Using with RDKit

SMSD works standalone or alongside RDKit. Use RDKit for parsing and drawing,
SMSD for fast matching:

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

opts = smsd.ChemOptions.tautomer_profile()
opts.solvent = smsd.Solvent.DMSO       # adjust for DMSO
opts.with_ph(5.0)                       # adjust for pH 5.0

mcs = smsd.mcs("CC(=O)C", "CC(O)=C", chem=opts)
```

Supported solvents: `AQUEOUS`, `DMSO`, `METHANOL`, `CHLOROFORM`, `ACETONITRILE`, `DIETHYL_ETHER`

## Native MOL/SDF I/O

```python
import smsd

g = smsd.read_molfile("input.mol")
mol_block = smsd.write_mol_block(g)
mol_block_v3000 = smsd.write_mol_block_v3000(g)
sdf_record = smsd.write_sdf_record(g)

smsd.write_molfile(g, "out_v2000.mol")
smsd.write_molfile(g, "out_v3000.mol", v3000=True)
smsd.write_molfile(g, "out.sdf", sdf=True)
```

The native writer preserves practical chemistry metadata in `6.12.0`:
- names, comments, and SDF properties
- charges, isotopes, atom classes, and atom maps
- `R#`/`R<n>` plus `M  RGP`
- practical V2000/V3000 stereo round-trip

## Platform & GPU Support

SMSD automatically dispatches to the best available compute backend:

| Platform | CPU | GPU |
|----------|-----|-----|
| macOS (Apple Silicon) | OpenMP | Metal (zero-copy unified memory) |
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

# GPU is used automatically for:
# - Domain initialization in VF2++ substructure search
# - RASCAL batch screening
# - NLF histogram computation

# CPU fallback is seamless — no code changes needed
results = smsd.batch_substructure(query, targets)  # uses GPU if available
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
mapping = smsd.mcs("SMILES1", "SMILES2", tautomer_aware=True)

# Similarity
sim = smsd.similarity("SMILES1", "SMILES2")
ub = smsd.similarity_upper_bound(mol1, mol2)
hits = smsd.screen_targets(query, library, threshold=0.5)

# Fingerprints
fp = smsd.path_fingerprint(mol, path_length=7, fp_size=2048)
fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)           # ECFP4
fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048, mode="fcfp")  # FCFP4
sim = smsd.overlapCoefficient(fp1, fp2)
ok = smsd.fingerprint_subset(query_fp, target_fp)

# Format conversions (database storage, REST APIs)
hex_str = smsd.to_hex(fp, fp_size=2048)
fp_back = smsd.from_hex(hex_str)
bits    = smsd.to_binary_string(fp, fp_size=2048)

# Batch (OpenMP parallel)
results = smsd.batch_substructure(query, targets)
results = smsd.batch_mcs(query, targets)
fps     = smsd.batch_fingerprint(mols, path_length=7, fp_size=2048)
hits    = smsd.batch_fingerprint_screen(query_fp, target_fps)

# RASCAL pre-screen + exact MCS in one call
matches = smsd.screen_and_match(query, targets, threshold=0.5)

# Batch find substructure with atom-atom mappings (v6.12.0)
mappings = smsd.batch_find_substructure(query, targets)

# TargetCorpus — prewarm once, query many times (v6.12.0)
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

## Complementary Strengths

SMSD is designed to work alongside existing toolkits, not replace them.
Each toolkit brings unique strengths to the cheminformatics ecosystem:

| Capability | SMSD Pro | RDKit | CDK |
|---|---|---|---|
| Tautomer-aware circular FP | Yes | — | — |
| pH/solvent-sensitive matching | Yes | — | — |
| Multi-level MCS pipeline | 11 levels | FMCS | MCSPlus |
| GPU acceleration | CUDA + Metal | — | — |
| Descriptor calculation | — | Extensive | Extensive |
| Reaction handling | Basic | Comprehensive | Comprehensive |
| 3D conformers | — | Yes | Yes |
| Header-only C++ | Yes | — | — |

**Recommended workflow:** Use RDKit or CDK for parsing, descriptors, and 3D — use SMSD for MCS and substructure matching.

## Also Available

- **Java**: `com.bioinceptionlabs:smsd:6.12.0` on [Maven Central](https://central.sonatype.com/artifact/com.bioinceptionlabs/smsd)
- **C++**: Header-only, zero dependencies — [GitHub](https://github.com/asad/SMSD)

## Citation

If you use SMSD Pro in your research, please cite:

> Rahman SA.
> *SMSD Pro: Coverage-Driven, Tautomer-Aware Maximum Common Substructure Search.*
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

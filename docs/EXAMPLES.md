# SMSD Pro — Examples, How-To, and Cautions

**Version 6.11.2** | Copyright (c) 2018-2026 BioInception PVT LTD

This document provides worked examples for every major SMSD Pro feature. Each section
includes runnable code, expected output, practical cautions, and performance notes.

---

## Table of Contents

- [1. Quick Start](#1-quick-start)
- [2. Substructure Search](#2-substructure-search)
- [3. Maximum Common Substructure (MCS)](#3-maximum-common-substructure-mcs)
- [4. MCS Variants](#4-mcs-variants)
- [5. Tautomer-Aware MCS](#5-tautomer-aware-mcs)
- [6. Fingerprints and Similarity](#6-fingerprints-and-similarity)
- [7. Depiction (SVG)](#7-depiction-svg)
- [8. 2D / 3D Layout](#8-2d--3d-layout)
- [9. Stereo and CIP Assignment](#9-stereo-and-cip-assignment)
- [10. File I/O (MOL, SDF)](#10-file-io-mol-sdf)
- [11. R-Group Decomposition](#11-r-group-decomposition)
- [12. Reaction Atom Mapping](#12-reaction-atom-mapping)
- [13. Batch Operations](#13-batch-operations)
- [14. SMARTS Matching](#14-smarts-matching)
- [15. Scaffold Analysis](#15-scaffold-analysis)
- [16. Graph Utilities](#16-graph-utilities)
- [17. Performance Tuning](#17-performance-tuning)
- [18. Common Cautions](#18-common-cautions)

---

## 1. Quick Start

### Python

```python
import smsd

# Substructure: is benzene a substructure of phenol?
assert smsd.is_substructure(
    smsd.parse_smiles("c1ccccc1"),    # query
    smsd.parse_smiles("c1ccc(O)cc1")  # target
)

# MCS: what do aspirin and salicylic acid share?
mapping = smsd.find_mcs(
    smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O"),  # aspirin
    smsd.parse_smiles("Oc1ccccc1C(=O)O")          # salicylic acid
)
print(f"MCS size: {len(mapping)} atoms")  # 10

# Similarity
fp1 = smsd.circular_fingerprint("CC(=O)Oc1ccccc1C(=O)O", radius=2, fp_size=2048)
fp2 = smsd.circular_fingerprint("Oc1ccccc1C(=O)O", radius=2, fp_size=2048)
print(f"Tanimoto: {smsd.tanimoto(fp1, fp2):.3f}")
```

### Java

```java
import com.bioinception.smsd.core.*;

SMSD smsd = new SMSD(mol1, mol2, new ChemOptions());
boolean isSub = smsd.isSubstructure();
var mapping = smsd.findMCS();
```

### C++

```cpp
#include "smsd/smsd.hpp"

auto mol1 = smsd::parseSMILES("c1ccccc1");
auto mol2 = smsd::parseSMILES("c1ccc(O)cc1");

bool isSub = smsd::isSubstructure(mol1, mol2, smsd::ChemOptions{});
auto mcs = smsd::findMCS(mol1, mol2, smsd::ChemOptions{}, smsd::MCSOptions{});
```

---

## 2. Substructure Search

### Basic substructure check

```python
import smsd

query  = smsd.parse_smiles("c1ccccc1")      # benzene ring
target = smsd.parse_smiles("c1ccc(N)cc1")   # aniline

hit = smsd.is_substructure(query, target)
print(hit)  # True

# Get the atom mapping
mapping = smsd.find_substructure(query, target)
print(mapping)  # {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 6} (example)
```

### All substructure matches

```python
# Find up to 10 distinct mappings
all_matches = smsd.find_all_substructures(query, target, max_matches=10)
print(f"Found {len(all_matches)} distinct mappings")
```

### Caution: ring-matches-ring

By default, `ringMatchesRingOnly=True` — a ring atom in the query will only
match ring atoms in the target. This prevents false positives where a benzene
ring matches a cyclohexane chain.

```python
# This will NOT match with default options (ring vs. non-ring)
ring   = smsd.parse_smiles("c1ccccc1")  # aromatic ring
chain  = smsd.parse_smiles("C1CCCCC1")  # non-aromatic ring
print(smsd.is_substructure(ring, chain))  # False with default options
```

---

## 3. Maximum Common Substructure (MCS)

### Basic MCS

```python
import smsd

mol1 = smsd.parse_smiles("c1ccc2c(c1)cc1ccccc1c2")  # phenanthrene
mol2 = smsd.parse_smiles("c1cc2ccc3cccc4ccc(c1)c2c34")  # pyrene

mapping = smsd.find_mcs(mol1, mol2)
print(f"MCS size: {len(mapping)} atoms")

# Extract MCS as SMILES
mcs_smiles = smsd.mcs_to_smiles(mol1, mapping)
print(f"MCS SMILES: {mcs_smiles}")
```

### Structured result

```python
result = smsd.mcs_result("c1ccccc1", "c1ccc(O)cc1")
print(f"Size: {result.size}")
print(f"Tanimoto: {result.tanimoto:.3f}")
print(f"MCS SMILES: {result.mcs_smiles}")
print(f"Mapping: {result.mapping}")
```

### With timeout

```python
# For very large molecules, set a timeout to prevent hanging
mapping = smsd.find_mcs(mol1, mol2, timeout_ms=5000)  # 5 seconds max
```

### Caution: MCS size vs. similarity

MCS size alone does not indicate similarity — always normalise. A 6-atom MCS
between two 50-atom molecules is a poor match, but a 6-atom MCS between a
6-atom and a 7-atom molecule is excellent.

```python
# Always compute Tanimoto for meaningful comparison
tanimoto = len(mapping) / (len(mol1) + len(mol2) - len(mapping))
```

---

## 4. MCS Variants

```python
import smsd

mol1 = smsd.parse_smiles("c1ccccc1")
mol2 = smsd.parse_smiles("c1ccc(O)cc1")

# Connected MCS (default) — all matched atoms form a single connected component
mcs = smsd.mcs(mol1, mol2)

# Disconnected MCS — matched atoms may be in separate fragments
mcs = smsd.mcs(mol1, mol2, connected_only=False)

# Induced MCS — preserved bond orders between matched atoms
mcs = smsd.mcs(mol1, mol2, induced=True)

# Edge MCS (MCES) — maximise matched bonds rather than atoms
mcs = smsd.mcs(mol1, mol2, maximize_bonds=True)

# Find top-5 distinct MCS solutions
all_mcs = smsd.find_all_mcs(mol1, mol2, max_results=5)
for i, m in enumerate(all_mcs):
    print(f"Solution {i+1}: {len(m)} atoms")
```

### Multi-Molecule MCS (N-MCS)

```python
# Find the common substructure shared by 3+ molecules
molecules = [
    smsd.parse_smiles("c1ccc(O)cc1"),   # phenol
    smsd.parse_smiles("c1ccc(N)cc1"),   # aniline
    smsd.parse_smiles("c1ccc(Cl)cc1"),  # chlorobenzene
]
nmcs = smsd.find_nmcs(molecules, threshold=0.5, timeout_ms=10000)
print(f"N-MCS size: {len(nmcs)} atoms")
```

### Caution: disconnected MCS

Disconnected MCS can return fragments that have no chemical meaning when
considered independently. Always check fragment connectivity for downstream use.

---

## 5. Tautomer-Aware MCS

```python
import smsd

# Keto-enol tautomers: acetone and prop-1-en-2-ol
keto = smsd.parse_smiles("CC(=O)C")
enol = smsd.parse_smiles("CC(O)=C")

# Without tautomer awareness — misses the equivalence
mcs_strict = smsd.mcs(keto, enol)
print(f"Strict MCS: {len(mcs_strict)} atoms")

# With tautomer awareness — recognises keto-enol equivalence
mcs_tauto = smsd.mcs(keto, enol, tautomer_aware=True)
print(f"Tautomer MCS: {len(mcs_tauto)} atoms")  # larger

# Solvent-aware tautomer equilibrium
mcs = smsd.mcs(keto, enol, tautomer_aware=True, solvent="DMSO", pH=7.4)
```

### Caution: tautomer-aware is slower

Tautomer-aware MCS explores up to 15 tautomeric transforms per atom. This is
typically 2-5x slower than strict MCS. For large-scale screening, use strict
mode first, then tautomer-aware on the top hits only.

---

## 6. Fingerprints and Similarity

### ECFP4 (Morgan) fingerprint

```python
import smsd

# Binary ECFP4 (2048 bits, radius 2)
fp = smsd.circular_fingerprint("c1ccc(O)cc1", radius=2, fp_size=2048)

# Count-based ECFP4 (better for ML)
counts = smsd.ecfp_counts("c1ccc(O)cc1", radius=2, fp_size=2048)

# FCFP4 (pharmacophore-aware)
fcfp = smsd.circular_fingerprint("c1ccc(O)cc1", radius=2, fp_size=2048, mode="fcfp")

# Topological torsion
torsion = smsd.topological_torsion("c1ccc(O)cc1", fp_size=2048)
```

### Similarity metrics

```python
fp1 = smsd.circular_fingerprint("CCO", radius=2)
fp2 = smsd.circular_fingerprint("CCCO", radius=2)

tan  = smsd.tanimoto(fp1, fp2)          # Jaccard index
dice = smsd.dice(fp1, fp2)              # Dice coefficient
cos  = smsd.cosine(fp1, fp2)            # Cosine similarity

# Count-based metrics (use count vectors, not binary)
c1 = smsd.ecfp_counts("CCO", radius=2)
c2 = smsd.ecfp_counts("CCCO", radius=2)
ct = smsd.count_tanimoto(c1, c2)
cd = smsd.count_dice(c1, c2)
```

### RASCAL upper bound (fast pre-filter)

```python
# Fast similarity upper bound — avoids full MCS on dissimilar pairs
ub = smsd.similarity_upper_bound(mol1, mol2)
if ub > 0.3:
    # Only compute full MCS if upper bound suggests a reasonable match
    mcs = smsd.find_mcs(mol1, mol2)
```

### Caution: fingerprint choice matters

- **ECFP4** (radius=2): best for small-molecule similarity
- **FCFP4**: better for pharmacophore-level comparison
- **Topological torsion**: best for peptides and linear chains
- **Count-based**: always prefer over binary for machine learning models

---

## 7. Depiction (SVG)

SMSD includes a zero-dependency SVG renderer conforming to ACS 1996 standard
(Nature, Science, JACS, Springer). No external tools required.

### Single molecule

```python
import smsd

# From SMILES (auto-layout)
svg = smsd.depict_svg("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
smsd.save_svg(svg, "aspirin.svg")

# From MolGraph
mol = smsd.parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C")  # caffeine
svg = smsd.depict_svg(mol)
smsd.save_svg(svg, "caffeine.svg")
```

### MCS comparison (side-by-side)

```python
mol1 = smsd.parse_smiles("c1ccccc1")
mol2 = smsd.parse_smiles("c1ccc(O)cc1")
mapping = smsd.find_mcs(mol1, mol2)

svg = smsd.depict_pair(mol1, mol2, mapping)
smsd.save_svg(svg, "benzene_vs_phenol.svg")
```

Both molecules are rendered with matched atoms highlighted in green, matched
bonds drawn bold, and blue superscript numbers showing the atom-atom correspondence.

### Substructure highlighting

```python
query  = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(NC(=O)C)cc1")  # acetanilide
mapping = smsd.find_substructure(query, target)

svg = smsd.depict_mapping(target, mapping)
smsd.save_svg(svg, "acetanilide_highlight.svg")
```

### Custom styling

```python
# All ACS 1996 proportions auto-scale from bond_length
svg = smsd.depict_svg("c1ccc2c(c1)cc1ccccc1c2",  # phenanthrene
    bond_length=50,           # larger for poster / slide
    width=800,                # fixed canvas size
    height=400,
    padding=40,
    show_atom_indices=True,   # debug: show atom numbers
    font_family="Times New Roman, serif"
)

# Using DepictOptions for batch consistency
opts = smsd.DepictOptions()
opts.bond_length = 35
opts.show_map_numbers = False

for smi in ["CCO", "c1ccccc1", "CC(=O)O"]:
    svg = smsd.depict_svg(smi, opts=opts)
    smsd.save_svg(svg, f"{smi.replace('(', '').replace(')', '')}.svg")
```

### Converting SVG to PNG (external tools)

```bash
# Inkscape (highest quality for publication)
inkscape molecule.svg --export-type=png --export-dpi=600

# ImageMagick
convert -density 600 molecule.svg molecule.png

# cairosvg (Python)
pip install cairosvg
python -c "import cairosvg; cairosvg.svg2png(url='molecule.svg', write_to='molecule.png', dpi=600)"
```

### Caution: font rendering

SVG text rendering depends on the viewer's installed fonts. Arial and Helvetica
are near-universal. If using a custom font, embed it as a base64 web font in
the SVG or ensure the font is installed on the target system.

### Caution: large molecules

For molecules with >100 atoms, the built-in auto-layout may produce suboptimal
results. Use `generate_coords_2d()` first for better placement, then pass
pre-computed coordinates to the renderer.

---

## 8. 2D / 3D Layout

### Generate 2D coordinates

```python
import smsd

mol = smsd.parse_smiles("c1ccc2c(c1)cc1ccccc1c2")  # phenanthrene
coords = smsd.generate_coords_2d(mol, target_bond_length=1.5)
# coords is a list of (x, y) tuples, one per atom

# Check layout quality
quality = smsd.layout_quality(mol, coords)
print(f"Layout quality score: {quality:.3f}")  # lower is better
```

### Generate 3D coordinates

```python
mol = smsd.parse_smiles("c1ccccc1")
coords_3d = smsd.generate_coords_3d(mol, target_bond_length=1.5)
# coords_3d is a list of (x, y, z) tuples
```

### Layout algorithms

```python
# Force-directed (Fruchterman-Reingold + crossing penalty)
coords = smsd.force_directed_layout(mol, coords, max_iter=500, target_bond_length=1.5)

# Stress majorisation (SMACOF — globally optimal distances)
coords = smsd.stress_majorisation(mol, max_iter=300, target_bond_length=1.5)

# Simulated annealing crossing reduction
coords = smsd.reduce_crossings(mol, coords, max_iter=2000)
```

### Coordinate transforms

```python
# All transforms operate on (x, y) coordinate lists
coords = smsd.translate_2d(mol, coords, dx=10.0, dy=5.0)
coords = smsd.rotate_2d(mol, coords, angle_degrees=45.0)
coords = smsd.scale_2d(mol, coords, factor=2.0)
coords = smsd.mirror_x(mol, coords)
coords = smsd.mirror_y(mol, coords)
coords = smsd.center_2d(mol, coords)
coords = smsd.normalise_bond_length(mol, coords, target=1.5)
coords = smsd.canonical_orientation(mol, coords)

# Align one molecule onto another (Kabsch rotation)
coords = smsd.align_2d(mol, coords, ref_coords, ref_indices, mol_indices)

# Bounding box
bbox = smsd.bounding_box_2d(mol, coords)  # (min_x, min_y, max_x, max_y)
```

### Caution: layout is non-deterministic

Force-directed and SMACOF layouts use random initialisation. Results may vary
between runs. For reproducible output, use `generate_coords_2d()` which uses
template matching and deterministic placement.

---

## 9. Stereo and CIP Assignment

```python
import smsd

# R/S tetrahedral chirality
mol = smsd.parse_smiles("N[C@@H](C)C(=O)O")  # L-alanine
stereo = smsd.assign_rs(mol)
print(stereo)  # {1: 'S'}

# E/Z double bond geometry
mol = smsd.parse_smiles("C/C=C/C")  # E-2-butene
ez = smsd.assign_ez(mol)
print(ez)  # {(1,2): 'E'}

# Full CIP assignment from SMILES (convenience)
result = smsd.assign_cip_from_smiles("N[C@@H](C)C(=O)O")
```

### Caution: CIP requires 3D or stereo flags

CIP assignment operates on stored stereo flags (`@`, `@@`, `/`, `\` in SMILES).
If your molecule was parsed from a SMILES without stereo notation, no
assignments will be returned. Ensure input SMILES contain explicit stereo
markers when stereo analysis is required.

---

## 10. File I/O (MOL, SDF)

### Read/write MOL files

```python
import smsd

# Read from file
mol = smsd.read_mol_file("molecule.mol")

# Read from string
mol_block = """
  molecule
     BioInception

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0
    3.0000    0.0000    0.0000 O   0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
"""
mol = smsd.read_mol_block(mol_block)

# Write V2000
block = smsd.write_mol_block(mol)
print(block)

# Write V3000
block_v3 = smsd.write_mol_block_v3000(mol)
```

### Read/write SDF files

```python
# Read an SDF (returns list of MolGraph)
mols = smsd.read_sdf("compounds.sdf")

# Write molecules to SDF
smsd.write_sdf(mols, "output.sdf")
```

### Caution: V2000 atom limit

V2000 format has a hard limit of 999 atoms. For larger molecules, use V3000
format or SMILES.

---

## 11. R-Group Decomposition

```python
import smsd

core = smsd.parse_smiles("c1ccccc1")  # benzene core
targets = [
    smsd.parse_smiles("c1ccc(O)cc1"),    # phenol
    smsd.parse_smiles("c1ccc(N)cc1"),    # aniline
    smsd.parse_smiles("c1ccc(Cl)cc1"),   # chlorobenzene
]

results = smsd.decompose_rgroups(core, targets)
for r in results:
    print(f"Core: {r['core']}, R-groups: {r['rgroups']}")
```

### Caution: core must be a valid substructure

R-group decomposition requires the core to be a substructure of every target.
If a target does not contain the core, it will be skipped silently.

---

## 12. Reaction Atom Mapping

```python
import smsd

# Map atoms between reactant and product
reactant = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
product  = smsd.parse_smiles("Oc1ccccc1C(=O)O")         # salicylic acid

mapping = smsd.map_reaction(reactant, product)
print(mapping)  # {reactant_atom: product_atom, ...}
```

### Caution: reaction mapping is NP-hard

For complex multi-step reactions with many bond changes, the mapping may
time out. Use `timeout_ms` to bound the computation and accept a
best-effort result.

---

## 13. Batch Operations

```python
import smsd

# Parse a library
library = [smsd.parse_smiles(s) for s in [
    "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1",
    "c1ccc(Cl)cc1", "c1ccc(F)cc1"
]]

# Pre-warm for repeated use (10-30% speedup)
for mol in library:
    smsd.prewarm_graph(mol)

# Batch substructure screening
query = smsd.parse_smiles("c1ccccc1")
hits = smsd.batch_substructure(query, library)
print(f"Hits: {sum(hits)}/{len(library)}")

# Batch MCS
mcs_results = smsd.batch_mcs(library[0], library[1:])

# Batch fingerprint screening (RASCAL upper bound)
fps = [smsd.circular_fingerprint(mol, radius=2) for mol in library]
```

### Caution: prewarm for batch workflows

Always call `prewarm_graph()` on molecules used in batch operations.
This pre-computes canonical hashes and NLF invariants, avoiding redundant
recomputation on every query.

---

## 14. SMARTS Matching

```python
import smsd

# Match a SMARTS pattern
mol = smsd.parse_smiles("c1ccc(O)cc1")

# Hydroxyl group
matches = smsd.smarts_match("[OH]", mol)
print(f"Found hydroxyl: {len(matches) > 0}")

# Find all occurrences of a pattern
all_matches = smsd.smarts_find_all("[#6]~[#6]", mol, max_matches=50)
print(f"C-C bonds: {len(all_matches)}")

# SMARTS-based MCS
mcs = smsd.find_mcs_smarts("[#6]~[#7]", mol)
```

### Caution: SMARTS performance

Complex SMARTS patterns with many wildcards can be expensive. Keep patterns
specific — use `[#6]` (any carbon) rather than `[*]` (any atom) when possible.

---

## 15. Scaffold Analysis

```python
import smsd

# Murcko scaffold (ring systems + linkers, side chains removed)
mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
scaffold = smsd.murcko_scaffold(mol)
print(f"Scaffold: {smsd.to_smiles(scaffold)}")

# Scaffold MCS — find shared scaffold between two molecules
s = smsd.find_scaffold_mcs(
    smsd.parse_smiles("c1ccc2c(c1)cccc2"),    # naphthalene
    smsd.parse_smiles("c1ccc2c(c1)cc1ccccc1c2") # phenanthrene
)
```

---

## 16. Graph Utilities

```python
import smsd

mol = smsd.parse_smiles("c1ccccc1.CCO")  # benzene + ethanol (disconnected)

# Count connected components
n = smsd.count_components(mol)
print(f"Components: {n}")  # 2

# Split into separate molecules
components = smsd.split_components(mol)
for c in components:
    print(f"  Fragment: {smsd.to_smiles(c)}")

# Check canonical equivalence
mol1 = smsd.parse_smiles("c1ccccc1")
mol2 = smsd.parse_smiles("C1=CC=CC=C1")
print(smsd.same_canonical_graph(mol1, mol2))  # True
```

---

## 17. Performance Tuning

### Pre-warm molecules for batch use

```python
smsd.prewarm_graph(mol)  # 10-30% latency reduction on repeated ops
```

### Use timeouts for safety

```python
# Always set timeouts on large molecules to prevent hanging
mapping = smsd.find_mcs(mol1, mol2, timeout_ms=10000)  # 10 second limit
```

### RASCAL pre-screening

```python
# For large-scale screening, use RASCAL upper bound first
ub = smsd.similarity_upper_bound(query, target)
if ub < 0.2:
    pass  # skip — guaranteed dissimilar
else:
    mcs = smsd.find_mcs(query, target)  # only compute MCS on promising pairs
```

### GPU acceleration

```python
if smsd.gpu_is_available():
    print(smsd.gpu_device_info())
    # GPU automatically used for RASCAL batch screening and domain init
```

---

## 18. Common Cautions

### Never use regex to strip atom maps

```python
# WRONG — regex can corrupt ring closures like "C1CC1" → "CCC"
import re
smi = re.sub(r':\d+', '', mapped_smiles)  # DO NOT DO THIS

# CORRECT — use the SMSD parser which preserves ring closures
mol = smsd.parse_smiles(mapped_smiles)
clean_smi = smsd.to_smiles(mol)  # atom maps stripped safely
```

### Timeout on large molecules

SMSD will never hang on any molecule. However, MCS computation is NP-hard.
For molecules with >80 heavy atoms, always set a timeout:

```python
mapping = smsd.find_mcs(mol1, mol2, timeout_ms=5000)
```

### Ring-matches-ring (default)

The default `ringMatchesRingOnly=True` prevents a ring atom from matching a
chain atom. This is correct for most chemical applications. Use `fmcsProfile()`
only when you explicitly want loose FMCS-style topology.

### Tautomer-aware is 2-5x slower

Only enable tautomer-aware mode when chemical equivalence of tautomers is
required. For pure topological comparison, strict mode is faster and sufficient.

### Fingerprint choice for ML

Always prefer **count-based** fingerprints (`ecfp_counts`, `count_tanimoto`)
over binary fingerprints for machine learning models. Count vectors preserve
frequency information that binary vectors discard.

### V2000 atom limit

V2000 MOL format is limited to 999 atoms. For larger molecules, use V3000 or
SMILES format.

### Font rendering in SVG

SVG depiction uses Arial/Helvetica by default. If the viewing system lacks
these fonts, atom labels may render in a fallback font. For publication,
convert SVG to PDF/PNG using Inkscape at 600 DPI.

---

## Language-Specific Guides

For in-depth API reference and language-specific examples:

| Language | Guide |
|----------|-------|
| Python | [docs/PYTHON.md](PYTHON.md) |
| Java | [docs/JAVA.md](JAVA.md) |
| C++ | [docs/CPP.md](CPP.md) |

---

*SMSD Pro by BioInception PVT LTD. Algorithm Copyright 2009-2026 Syed Asad Rahman. Apache-2.0.*

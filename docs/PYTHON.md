# SMSD Pro Python Guide

SMSD Pro for Python exposes the native C++ engines for substructure, MCS,
fingerprinting, SMARTS matching, stereo/CIP assignment, layout, and chemistry
I/O. The default path is native SMSD with no RDKit dependency. RDKit can still
be used alongside SMSD for interop, cross-validation, or depiction workflows.

## Support

- CPython `3.9` through the latest stable release series
- Recommended default runtime: `Python 3.12`
- CPU-only execution by default
- Optional GPU acceleration:
  - macOS: Metal
  - Linux/Windows: CUDA
- CPU fallback remains available when no GPU backend is compiled

## Install

```bash
pip install smsd
```

Optional extras:

```bash
pip install "smsd[rdkit]"
pip install "smsd[gpu]"
```

## API Quick Reference

| Category | Functions |
|----------|-----------|
| **Parsing/I/O** | `parse_smiles`, `to_smiles`, `read_mol_block`, `write_mol_block`, `write_mol_block_v3000`, `read_mol_file`, `read_sdf`, `write_sdf`, `write_sdf_record` |
| **Substructure** | `is_substructure`, `find_substructure`, `find_all_substructures` |
| **MCS** | `find_mcs`, `find_all_mcs`, `find_mcs_progressive`, `find_mcs_smiles`, `mcs_size`, `mcs_to_smiles`, `find_nmcs`, `find_scaffold_mcs` |
| **Mapping** | `validate_mapping`, `is_mapping_maximal`, `canonicalize_mapping`, `translate_to_atom_ids` |
| **Reaction** | `map_reaction`, `reaction_aware_mcs`, `map_reaction_aware`, `bond_change_score` |
| **Scaffold/R-group** | `murcko_scaffold`, `decompose_rgroups`, `extract_subgraph` |
| **Graph Utils** | `count_components`, `split_components`, `same_canonical_graph`, `prewarm_graph` |
| **Fingerprints** | `circular_fingerprint`, `circular_fingerprint_counts`, `topological_torsion`, `topological_torsion_counts`, `path_fingerprint`, `mcs_fingerprint`, `counts_to_array` |
| **Similarity** | `tanimoto`, `dice`, `cosine`, `soergel`, `count_tanimoto`, `count_dice`, `count_cosine`, `similarity_upper_bound` |
| **Batch** | `batch_substructure`, `batch_mcs`, `batch_mcs_size`, `batch_mcs_constrained`, `batch_fingerprint`, `batch_fingerprint_screen`, `screen_targets`, `screen_and_match`, `screen_and_mcs_size` |
| **Stereo/CIP** | `assign_rs`, `assign_ez`, `assign_cip`, `assign_cip_from_smiles` |
| **SMARTS** | `smarts_match`, `smarts_find_all`, `find_mcs_smarts`, `compile_smarts`, `to_smarts` |
| **Layout** | `generate_coords_2d`, `generate_coords_3d`, `force_directed_layout`, `stress_majorisation`, `reduce_crossings`, `match_template`, `compute_sssr`, `layout_sssr`, `count_crossings`, `is_degenerate_layout`, `resolve_overlaps`, `layout_quality` |
| **Transforms** | `translate_2d`, `rotate_2d`, `scale_2d`, `mirror_x`, `mirror_y`, `center_2d`, `align_2d`, `bounding_box_2d`, `normalise_bond_length`, `canonical_orientation` |
| **Chemistry** | `classify_pharmacophore`, `implicit_h`, `validate_tautomer_consistency` |
| **Depiction** | `depict_svg`, `depict_pair`, `depict_mapping`, `save_svg`, `DepictOptions` |
| **GPU** | `gpu_is_available`, `gpu_device_info` |
| **Clique Solver** *(6.12.1)* | `find_max_cliques`, `find_mcs_clique`, `match_substructure`, `match_substructure_from_elements`, `score_mapping` |

## Options Reference

### ChemOptions — what counts as a "match"

Controls atom/bond compatibility during substructure and MCS search.

```python
import smsd

chem = smsd.ChemOptions()

# Atom-level matching
chem.match_atom_type    = True     # carbon only matches carbon
chem.match_formal_charge = False   # ignore charge differences (good for reactions)
chem.match_isotope      = False    # ignore mass number
chem.use_chirality      = False    # ignore R/S
chem.use_bond_stereo    = False    # ignore E/Z

# Bond-level matching
chem.match_bond_order = smsd.BondOrderMode.STRICT    # single≠double≠triple
                      # smsd.BondOrderMode.LOOSE     # aromatic ≈ single/double
                      # smsd.BondOrderMode.ANY       # any bond matches any bond

# Ring constraints
chem.ring_matches_ring_only = True   # ring atom only matches ring atom
chem.complete_rings_only    = False  # MCS must include full rings (not partial)
chem.ring_fusion_mode = smsd.RingFusionMode.IGNORE      # don't care about fusion
                      # smsd.RingFusionMode.PERMISSIVE  # fused ≈ non-fused ok
                      # smsd.RingFusionMode.STRICT      # fused must match fused

# Aromaticity — two separate concepts:
#   AromaticityModel = which algorithm detects aromaticity (perception)
#   AromaticityMode  = how strict aromatic matching is during search
chem.aromaticity_model = smsd.AromaticityModel.DAYLIGHT_LIKE  # only model for now
chem.aromaticity_mode  = smsd.AromaticityMode.FLEXIBLE  # aromatic ≈ Kekulé ok
                       # smsd.AromaticityMode.STRICT    # aromatic must match aromatic

# Tautomer handling
chem.tautomer_aware = False     # set True for keto/enol, amide, etc.
chem.pH             = 7.4       # affects tautomer equilibrium weights
chem.solvent        = smsd.Solvent.AQUEOUS  # DMSO, METHANOL, CHLOROFORM, etc.

# Engine
chem.matcher_engine = smsd.MatcherEngine.VF2PP  # fastest (default)
                    # smsd.MatcherEngine.VF2    # classic
                    # smsd.MatcherEngine.VF3    # experimental

# Pruning (advanced — usually leave defaults)
chem.use_two_hop_nlf   = True    # 2-hop neighbourhood label frequency
chem.use_three_hop_nlf = False   # 3-hop (slower, stronger pruning)
chem.use_bit_parallel_feasibility = True   # bitset-accelerated AC-3
chem.induced           = False   # induced subgraph isomorphism

# Profiles (shortcuts)
strict = smsd.ChemOptions.profile("strict")      # exact matching
taut   = smsd.ChemOptions.tautomer_profile()      # tautomer-aware defaults
```

### MCSOptions — how to search

Controls the MCS search strategy and stopping criteria.

```python
mcs_opts = smsd.MCSOptions()

# Basic
mcs_opts.timeout_ms     = 1000     # wall-clock limit (-1 = no limit)
mcs_opts.induced        = False    # induced MCS (no extra edges)
mcs_opts.connected_only = True     # connected MCS only
mcs_opts.disconnected_mcs = False  # allow disconnected MCS fragments
mcs_opts.maximize_bonds = False    # edge MCS (MCES) instead of atom MCS

# Fragment control (for disconnected MCS)
mcs_opts.min_fragment_size = 1     # minimum atoms per fragment
mcs_opts.max_fragments     = 100   # cap number of fragments

# Pipeline depth — higher = more thorough, slower
mcs_opts.max_stage = 5   # 0=identity only, 1=substructure (fast for reactions),
                          # 2=McSplit, 3=Bron-Kerbosch, 4=McGregor, 5=full pipeline

# Near-MCS exploration (reaction-aware)
mcs_opts.near_mcs_delta      = 2    # try K-1, K-2 variants
mcs_opts.near_mcs_candidates = 20   # how many near-MCS to evaluate

# Chemistry-aware features
mcs_opts.reaction_aware   = False   # prefer heteroatom mappings
mcs_opts.bond_change_aware = False  # penalise implausible bond changes
mcs_opts.extra_seeds       = True   # try additional seed strategies

# Seed search tuning (advanced)
mcs_opts.seed_neighborhood_radius = 2    # local neighbourhood for seed extension
mcs_opts.seed_max_anchors         = 12   # max anchor points per seed
mcs_opts.template_fuzzy_atoms     = 0    # fuzzy atom tolerance

# NLF pruning during extension (advanced)
mcs_opts.use_two_hop_nlf_in_extension   = True
mcs_opts.use_three_hop_nlf_in_extension = False
```

### Aromaticity: Perception vs Matching

SMSD separates *detecting* aromaticity from *comparing* it during search:

```python
mol = smsd.parse_smiles("c1ccccc1")

# Perception — which algorithm detects aromatic atoms/bonds
# Currently only DAYLIGHT_LIKE (Hückel 4n+2 on planar rings)
mol.perceive_aromaticity(smsd.AromaticityModel.DAYLIGHT_LIKE)

# Matching — how strict aromatic comparison is in MCS/substructure
chem = smsd.ChemOptions()
chem.aromaticity_mode = smsd.AromaticityMode.STRICT    # aromatic must match aromatic
chem.aromaticity_mode = smsd.AromaticityMode.FLEXIBLE  # aromatic ≈ Kekulé single/double

# Both can be set on ChemOptions:
chem.aromaticity_model = smsd.AromaticityModel.DAYLIGHT_LIKE  # perception model
chem.aromaticity_mode  = smsd.AromaticityMode.FLEXIBLE        # matching strictness
```

## Lightweight MCS Engine (6.12.1)

High-level coverage-driven MCS with automatic LFUB termination.
Accepts SMILES, MolGraph, or RDKit Mol. Uses the C++ clique solver
for heavy combinatorial work; Python handles chemistry.

```python
from smsd.mcs_engine import find_mcs_lightweight

# From SMILES — zero boilerplate
result = find_mcs_lightweight("c1ccc(O)cc1", "c1ccc(N)cc1")
print(result.size)          # 6
print(result.mapping)       # [(1,1), (2,2), ...]
print(result.candidates)    # all candidate mappings
print(result.lfub)          # label-frequency upper bound
print(result.elapsed_ms)    # wall-clock time

# From MolGraph
import smsd
g1 = smsd.parse_smiles("c1ccc(O)cc1")
g2 = smsd.parse_smiles("c1ccc(N)cc1")
result = find_mcs_lightweight(g1, g2, timeout=2.0, ring_matches_ring=True)

# From RDKit Mol
from rdkit import Chem
result = find_mcs_lightweight(
    Chem.MolFromSmiles("c1ccccc1"),
    Chem.MolFromSmiles("c1ccc(O)cc1"),
    bond_any=True,
)
```

### Pipeline stages (automatic escalation)

| Stage | Algorithm | Typical coverage |
|-------|-----------|-----------------|
| L0.75 | Greedy atom-by-atom (Morgan rank + NLF) | 45% of pairs solved here |
| L1 | Substructure containment (C++ VF2 or RDKit) | +20% |
| L1.5 | Seed-and-extend from heteroatom anchors | +15% |
| L3 | C++ clique solver (BK + Tomita pivoting) | +15% |
| L4 | McGregor bond-grow backtracking | +5% |

Each stage checks against the label-frequency upper bound (LFUB).
If the current best reaches LFUB, search stops immediately.

## Clique Solver — Low-Level API (6.12.1)

Direct access to the C++ clique-based MCS and substructure engine.
Chemistry stays in Python; only the graph search runs in C++.

```python
import smsd._smsd as _smsd

# Substructure from element lists (C++ builds compat internally)
matches = _smsd.match_substructure_from_elements(
    ["C","C","C","C","C","C"],
    ["C","C","C","C","C","C","O"],
    {(0,1):1, (1,2):1, (2,3):1, (3,4):1, (4,5):1, (0,5):1},
    {(0,1):1, (1,2):1, (2,3):1, (3,4):1, (4,5):1, (0,5):1, (3,6):1},
    True, 500, 8,
)

# Substructure from pre-built compat pairs
matches = _smsd.match_substructure(
    n_query, n_target, compat,
    bonds_query, bonds_target,
    False, 500, 8,
)

# MCS via clique solver
result = _smsd.find_mcs_clique(
    compat, bonds_a, bonds_b, n_a, n_b,
    True, 1000, 8,
)
# result.best_size, result.candidates, result.elapsed_us

# Score a mapping
score = _smsd.score_mapping(
    mapping, elements_a, elements_b, bonds_a, bonds_b)
```

## Core Search

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(O)cc1")

assert smsd.is_substructure(query, target)
mapping = smsd.find_substructure(query, target)
mcs = smsd.find_mcs(query, target)
```

Convenience wrappers accept `MolGraph`, SMILES, or RDKit molecules:

```python
mcs = smsd.find_mcs("CC(=O)O", "CCC(=O)O")
hit = smsd.find_substructure("c1ccccc1", "c1ccc(O)cc1")
```

## SMARTS

```python
import smsd

assert smsd.smarts_match("[#6]~[#7]", "CCN")
matches = smsd.smarts_find_all("[OH]", "c1ccc(O)cc1")
largest = smsd.find_mcs_smarts("[CX3](=O)[NX3]", "CC(=O)Nc1ccc(O)cc1")
```

SMSD supports ring primitives (`R`, `r`, `x`), recursive SMARTS, logical
operators, aromatic atoms, charges, isotopes, and stereo-aware query handling
through the native matcher.

## Fingerprints

```python
import smsd

path_fp = smsd.path_fingerprint("c1ccccc1", path_length=7, fp_size=2048)
mcs_fp = smsd.mcs_fingerprint("c1ccc(O)cc1", path_length=7, fp_size=2048)
ecfp4 = smsd.fingerprint_from_smiles("c1ccccc1", radius=2, fp_size=2048)
torsion = smsd.topological_torsion("CC(N)C(=O)O")

assert smsd.fingerprint_subset(path_fp, path_fp)
score = smsd.tanimoto(ecfp4, ecfp4)
```

## Stereo and CIP

```python
import smsd

g = smsd.parse_smiles("N[C@@H](C)C(=O)O")
rs = smsd.assign_rs(g)
ez = smsd.assign_ez(smsd.parse_smiles("C/C=C/C"))
```

## Native MOL/SDF I/O

Use SMSD’s native readers/writers by default:

```python
import smsd

g = smsd.read_molfile("molecule.mol")
mol_block = smsd.write_mol_block(g)
mol_block_v3000 = smsd.write_mol_block_v3000(g)
sdf_record = smsd.write_sdf_record(g)

smsd.write_molfile(g, "out_v2000.mol")
smsd.write_molfile(g, "out_v3000.mol", v3000=True)
smsd.write_molfile(g, "out.sdf", sdf=True)
smsd.export_sdf([g], "library.sdf")
```

Covered natively in `6.11.0`:
- MOL V2000 and V3000 graph round-trip
- names, program line, comments, and SDF properties
- charges, isotopes, atom classes, atom maps
- `R#`/`R<n>` and `M  RGP`
- practical stereo flags

Not yet a complete native replacement for every exotic MDL query feature:
- full atom lists and not-lists
- link nodes
- variable attachment semantics
- the complete Markush query feature set

## R-group Fragments

```python
import smsd

frag = smsd.parse_smiles("[R1]c1ccccc1")
assert smsd.to_smiles(frag).startswith("[R1]")

rgroups = smsd.decompose_r_groups(
    "c1ccccc1",
    ["c1ccc(O)cc1", "c1ccc(N)cc1"],
)
```

## Multi-Molecule MCS (v6.11.0)

Find the common substructure shared across N molecules:

```python
import smsd

# Scaffold shared by a congeneric series
mols = [smsd.parse_smiles(s) for s in [
    "c1ccc(O)cc1",    # phenol
    "c1ccc(N)cc1",    # aniline
    "c1ccc(F)cc1",    # fluorobenzene
]]
nmcs = smsd.find_nmcs(mols, threshold=1.0, timeout_ms=30000)
print(f"Common scaffold: {len(nmcs)} atoms")
# Result: 6 (benzene ring shared by all three)
```

**Parameters:**
- `threshold` (0.0-1.0): fraction of molecules that must contain the MCS.
  Use 0.8 to tolerate one outlier in a set of 5.
- `timeout_ms`: per-pairwise-MCS timeout. Total time scales with N.

**Caution:** Runtime is O(N) pairwise MCS calls. For large sets (>50),
pre-filter with `screen_targets()` first.

**Impact:** Enables scaffold-hopping analysis, SAR series alignment, and
library clustering in a single call.

## Scaffold MCS (v6.11.0)

MCS on Murcko scaffolds only (strip side chains first):

```python
atorvastatin = smsd.parse_smiles("CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O")
rosuvastatin = smsd.parse_smiles("CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O")

scaffold_mcs = smsd.find_scaffold_mcs(atorvastatin, rosuvastatin)
print(f"Scaffold overlap: {len(scaffold_mcs)} atoms")
```

**Impact:** Focuses on pharmacophore-relevant ring frameworks, ignoring
substituent decoration. Better for series-to-series comparisons.

## Mapping Validation (v6.11.0)

Verify that an atom mapping is chemically correct:

```python
mol1 = smsd.parse_smiles("c1ccccc1")
mol2 = smsd.parse_smiles("c1ccc(O)cc1")
mapping = smsd.find_mcs(mol1, mol2)

# Validate: checks injectivity, atom/bond compatibility
errors = smsd.validate_mapping(mol1, mol2, mapping)
assert len(errors) == 0, f"Bad mapping: {errors}"

# Check maximality: can any more atom pairs be added?
is_max = smsd.is_mapping_maximal(mol1, mol2, mapping)
print(f"Maximal: {is_max}")
```

**Caution:** `validate_mapping` returns a list of human-readable error strings.
An empty list means the mapping is valid. Always validate mappings from external
sources (e.g., reaction databases) before downstream analysis.

**Impact:** Essential QA step for atom-atom mapping pipelines, especially when
integrating with external tools or reaction databases.

## R-Group Decomposition (v6.11.0)

Decompose molecules into a shared core plus variable R-groups:

```python
core = smsd.parse_smiles("c1ccccc1")   # benzene as core
mols = [
    smsd.parse_smiles("Cc1ccccc1"),     # toluene
    smsd.parse_smiles("c1ccc(O)cc1"),   # phenol
    smsd.parse_smiles("c1ccc(N)cc1"),   # aniline
]

results = smsd.decompose_rgroups(core, mols, timeout_ms=10000)
for i, r in enumerate(results):
    rg = r["rgroups"]  # dict: "R1" -> MolGraph, "R2" -> MolGraph, ...
    r_smiles = {k: smsd.to_smiles(v) for k, v in rg.items()}
    print(f"Mol {i}: R-groups = {r_smiles}")
# Output:
#   Mol 0: R-groups = {'R1': 'C'}
#   Mol 1: R-groups = {'R1': 'O'}
#   Mol 2: R-groups = {'R1': 'N'}
```

**Caution:** The core must be a valid substructure of each molecule. If a
molecule doesn't contain the core, the result has an empty core and no R-groups.
Check `len(r["core"]) > 0` before processing.

**Impact:** Enables SAR table generation, matched-molecular-pair analysis,
and automated medicinal chemistry reporting.

## Reaction Atom Mapping (v6.11.0)

Map atoms between reactants and products using disconnected MCS:

```python
reactant = smsd.parse_smiles("CCO")      # ethanol
product  = smsd.parse_smiles("CC=O")     # acetaldehyde

mapping = smsd.map_reaction(reactant, product, timeout_ms=10000)
# mapping: {reactant_atom_idx: product_atom_idx}
print(f"Mapped {len(mapping)} atoms")
```

**Caution:** Uses disconnected MCS internally, which may over-map for complex
multi-center reactions. For reaction-class-aware mapping with bond-change
scoring, use `smsd.reaction_aware_mcs()` and `smsd.bond_change_score()` instead.

## All Substructure Matches (v6.11.0)

Find every possible way a pattern embeds into a target:

```python
benzene = smsd.parse_smiles("c1ccccc1")
phenol = smsd.parse_smiles("c1ccc(O)cc1")

all_maps = smsd.find_all_substructures(benzene, phenol)
print(f"Found {len(all_maps)} embeddings")
# benzene in phenol: 12 (6 rotations x 2 reflections)
```

**Caution:** Returns up to 10,000 mappings. Highly symmetric molecules can
produce combinatorial explosion (e.g., C60 fullerene). For most workflows,
`find_substructure()` (single match) is sufficient and faster.

**Impact:** Required for exhaustive enumeration in reaction-site detection,
symmetry analysis, and crystallographic mapping.

## Graph Utilities (v6.11.0)

```python
import smsd

# Murcko scaffold: strip side chains, keep ring framework
ibuprofen = smsd.parse_smiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")
scaffold = smsd.murcko_scaffold(ibuprofen)
print(f"Scaffold: {smsd.to_smiles(scaffold)}")  # benzene ring

# Extract subgraph by atom indices
sub = smsd.extract_subgraph(ibuprofen, [4, 5, 6, 7, 8, 9])
print(f"Subgraph: {len(sub)} atoms")

# Connected components (separate salts/mixtures)
salt = smsd.parse_smiles("CC(=O)[O-].[Na+]")
n = smsd.count_components(salt)
parts = smsd.split_components(salt)
print(f"{n} components: {[smsd.to_smiles(p) for p in parts]}")

# Fast graph identity check (canonical SMILES comparison)
a = smsd.parse_smiles("c1ccccc1")
b = smsd.parse_smiles("C1=CC=CC=C1")
print(smsd.same_canonical_graph(a, b))  # True
```

**Caution:** `extract_subgraph` re-indexes atoms from 0. If you need to track
original indices, maintain a separate index map.

## File I/O (v6.11.0)

```python
import smsd

# Read single MOL file
mol = smsd.read_mol_file("molecule.mol")

# Read entire SDF library
library = smsd.read_sdf("compounds.sdf")
print(f"Loaded {len(library)} molecules")

# Write SDF
smsd.write_sdf(library, "output.sdf")
```

**Caution:** `read_sdf()` loads the entire file into memory. For very large
files (>100K molecules), process in chunks using `read_mol_block()` on
individual records instead.

## Chemistry Utilities (v6.11.0)

```python
import smsd

mol = smsd.parse_smiles("c1ccc(O)cc1")

# Pharmacophore classification (FCFP invariant)
# Returns bitmask: bit 0=donor, 1=acceptor, 2=pos, 3=neg, 4=aromatic, 5=hydrophobic
for i in range(len(mol)):
    features = smsd.classify_pharmacophore(mol, i)
    labels = []
    if features & 1:  labels.append("donor")
    if features & 2:  labels.append("acceptor")
    if features & 4:  labels.append("positive")
    if features & 8:  labels.append("negative")
    if features & 16: labels.append("aromatic")
    if features & 32: labels.append("hydrophobic")
    print(f"Atom {i}: {labels}")

# Implicit hydrogen count (MDL valence model)
print(smsd.implicit_h(6, 3, 0))  # Carbon, 3 bonds, no charge -> 1 H
print(smsd.implicit_h(7, 2, 0))  # Nitrogen, 2 bonds, no charge -> 1 H
print(smsd.implicit_h(8, 2, 0))  # Oxygen, 2 bonds, no charge -> 0 H
```

## 2D/3D Coordinate Generation (v6.11.0)

Generate publication-quality 2D or initial 3D coordinates:

```python
import smsd

mol = smsd.parse_smiles("c1ccc2c(c1)cc1ccccc1n2")  # acridine

# Full pipeline: template match -> ring layout -> chain zig-zag
# -> force refinement -> overlap resolution -> crossing reduction
# -> canonical orientation -> bond length normalisation
coords_2d = smsd.generate_coords_2d(mol, target_bond_length=1.5)
print(f"Quality: {smsd.layout_quality(mol, coords_2d):.4f}")  # 0 = perfect
print(f"Crossings: {smsd.count_crossings(mol, coords_2d)}")

# 3D coordinates via distance geometry + force-field refinement
coords_3d = smsd.generate_coords_3d(mol, target_bond_length=1.5)
```

**Caution:** `generate_coords_3d()` uses simplified distance geometry. For
production conformer generation (drug design, docking), use external tools
(ETKDG, MMFF, UFF). SMSD 3D provides a reasonable starting geometry.

## Coordinate Transforms (v6.11.0)

```python
import smsd, math

mol = smsd.parse_smiles("c1ccccc1")
coords = smsd.generate_coords_2d(mol)

# Translate, rotate, scale, mirror
moved   = smsd.translate_2d(coords, dx=5.0, dy=3.0)
rotated = smsd.rotate_2d(coords, angle=math.pi/4)
scaled  = smsd.scale_2d(coords, factor=2.0)
flipped = smsd.mirror_x(coords)    # flip vertically
flipped = smsd.mirror_y(coords)    # flip horizontally
centred = smsd.center_2d(coords)

# Align to a reference layout (minimise RMSD via optimal rotation)
rmsd, aligned = smsd.align_2d(coords, reference_coords)
print(f"RMSD after alignment: {rmsd:.4f}")

# Bounding box
min_x, min_y, max_x, max_y = smsd.bounding_box_2d(coords)

# Normalise bond lengths to target
normalised = smsd.normalise_bond_length(mol, coords, target=1.5)

# Canonical orientation (30-degree grid alignment, IUPAC convention)
canonical = smsd.canonical_orientation(mol, coords)
```

**Impact:** Full transform pipeline enables integration with any rendering
or molecular graphics framework. `align_2d()` is essential for superimposing
MCS mappings or comparing conformations.

## Layout Quality & Overlap Resolution (v6.11.0)

```python
import smsd

mol = smsd.parse_smiles("CC(C)(C)c1ccc(O)cc1")
coords = smsd.generate_coords_2d(mol)

# Quality score: 0.0 = perfect (bond uniformity + overlap + crossing count)
score = smsd.layout_quality(mol, coords)

# Count bond crossings
crossings = smsd.count_crossings(mol, coords)

# Check for degenerate layout
degenerate = smsd.is_degenerate_layout(coords)

# Resolve atom-atom overlaps
remaining, fixed = smsd.resolve_overlaps(coords, threshold=0.3, max_iter=100)
```

**Impact:** Use `layout_quality()` to benchmark different layout strategies.
`resolve_overlaps()` is a post-processing step for crowded molecules.

## Fingerprint Utilities (v6.11.0)

```python
import smsd

# Count fingerprint -> dense array (for ML pipelines)
counts = smsd.circular_fingerprint_counts("c1ccc(O)cc1", radius=2, fp_size=2048)
dense = smsd.counts_to_array(counts, 2048)
print(f"Feature vector length: {len(dense)}")  # 2048

# Pre-warm graph caches for batch operations
mols = [smsd.parse_smiles(s) for s in ["CCO", "CC=O", "CC(=O)O"]]
for m in mols:
    smsd.prewarm_graph(m)  # 10-30% latency reduction on repeated ops

# Then batch operations are faster:
results = smsd.batch_mcs(mols[0], mols[1:])
```

**Impact:** `prewarm_graph()` pre-computes canonical hashes and NLF invariants.
Call it once per molecule before repeated MCS/substructure operations to avoid
redundant recomputation. Measurable benefit on batch workflows (>10 molecules).

## Publication-Quality Depiction (v6.11.0)

SMSD includes a zero-dependency SVG renderer conforming to ACS 1996 standard
(the same specification used by Nature, Science, JACS, and Springer journals).
No RDKit or external tools required.

### Quick Start

```python
import smsd

# Render any molecule from SMILES
svg = smsd.depict_svg("c1ccc(O)cc1")
smsd.save_svg(svg, "phenol.svg")

# Render from MolGraph
mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
svg = smsd.depict_svg(mol)
```

### MCS Comparison (Side-by-Side)

```python
mol1 = smsd.parse_smiles("c1ccccc1")     # benzene
mol2 = smsd.parse_smiles("c1ccc(O)cc1")  # phenol
mapping = smsd.find_mcs(mol1, mol2)

svg = smsd.depict_pair(mol1, mol2, mapping)
smsd.save_svg(svg, "benzene_vs_phenol.svg")
```

Matched atoms and bonds are highlighted in green with numbered labels showing
the atom-atom correspondence.

### Substructure Highlighting

```python
query = smsd.parse_smiles("c1ccccc1")   # benzene ring
target = smsd.parse_smiles("c1ccc(N)cc1")  # aniline
mapping = smsd.find_mcs(query, target)

svg = smsd.depict_mapping(target, mapping)
smsd.save_svg(svg, "aniline_highlight.svg")
```

### Customizing Options

All ACS 1996 proportions auto-scale from `bond_length`. You can override any
parameter via the `DepictOptions` object or keyword arguments:

```python
# Method 1: Keyword arguments
svg = smsd.depict_svg("c1ccccc1",
    bond_length=40,        # larger molecules
    width=600,             # fixed SVG width
    height=400,            # fixed SVG height
    show_atom_indices=True # debugging aid
)

# Method 2: DepictOptions object (for reuse across many renders)
opts = smsd.DepictOptions()
opts.bond_length = 50           # publication: 30-50 typical
opts.font_size = 12             # 0 = auto from ACS ratio
opts.font_family = "Times New Roman, serif"
opts.show_carbon_labels = False # standard skeletal formula
opts.padding = 40

svg = smsd.depict_svg("Cn1cnc2c1c(=O)n(c(=O)n2C)C", opts=opts)  # caffeine
```

### Available Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `bond_length` | 30 | Reference unit (pixels). All ACS proportions scale from this. |
| `line_width` | 0 (auto) | Bond line width. Auto = ACS ratio (0.6/14.4 * bond_length). |
| `bold_width` | 0 (auto) | Wedge bond max width. Auto = ACS ratio (2.0/14.4 * bond_length). |
| `bond_spacing` | 0 (auto) | Double bond offset. Auto = 18% of bond_length. |
| `font_size` | 0 (auto) | Atom label font size. Auto = ACS ratio. |
| `subscript_scale` | 0.70 | H-count and charge subscript/superscript relative size. |
| `show_carbon_labels` | False | Show "C" at carbon positions. |
| `show_atom_indices` | False | Show atom index numbers (for debugging). |
| `show_map_numbers` | True | Show atom-atom mapping numbers. |
| `highlight_radius` | 0 (auto) | Atom highlight circle radius. Auto = 40% of bond_length. |
| `highlight_opacity` | 0.30 | Highlight fill opacity. |
| `match_width` | 0 (auto) | Highlighted bond width. Auto = 2.2x line_width. |
| `padding` | 30 | SVG padding in pixels. |
| `width` | 0 (auto) | Fixed SVG width. 0 = auto-fit to molecule. |
| `height` | 0 (auto) | Fixed SVG height. 0 = auto-fit to molecule. |
| `font_family` | Arial, Helvetica, sans-serif | CSS font family string. |

### ACS 1996 Standard

The renderer implements the ACS 1996 standard used by most chemistry journals:

- **Bond length:** 14.4 pt (all proportions are ratios of this)
- **Line width:** 0.6 pt
- **Wedge width:** 2.0 pt
- **Double bond spacing:** 18% of bond length
- **Chain angle:** 120 degrees
- **Font:** Arial/Helvetica, 10 pt (atom labels)
- **Subscripts:** 70% of label font size

### Rendering Features

- **Skeletal formula:** Carbon labels suppressed (standard for organic chemistry)
- **Heteroatom labels:** Element symbol + H-count subscript + charge superscript
- **H placement:** Automatic left/right based on bond directions (e.g., HO-, H₂N-)
- **Single bonds:** Round-capped lines at ACS line width
- **Double bonds:** Asymmetric offset toward ring interior (ring bonds), symmetric (chain bonds)
- **Triple bonds:** Center line + two thinner offset lines
- **Wedge bonds:** Solid filled triangle (up) or dashed stripes (down) for stereocenters
- **Aromatic circles:** Solid inner circle in aromatic rings
- **Bond-to-label clipping:** Bonds stop at label boundaries (no overlap)
- **White label backgrounds:** Clean masking of bonds behind atom labels
- **Jmol/CPK element colors:** Standard publication colors (N=blue, O=red, S=amber, etc.)
- **Atom map numbers:** Blue superscript numbering for MCS correspondences
- **SVG precision:** `geometricPrecision` rendering for sharp vector output

### Cautions

- SVG output is resolution-independent (vector). For raster (PNG), convert externally:
  ```bash
  # Using Inkscape (recommended for publication)
  inkscape phenol.svg --export-type=png --export-dpi=600

  # Using ImageMagick
  convert -density 600 phenol.svg phenol.png

  # Using cairosvg (Python)
  pip install cairosvg
  python -c "import cairosvg; cairosvg.svg2png(url='phenol.svg', write_to='phenol.png', dpi=600)"
  ```
- Auto-layout uses SMSD's internal layout engine. For large molecules (>100 atoms),
  consider using `generate_coords_2d()` first and passing pre-computed coordinates.
- Font rendering depends on the viewer's installed fonts. Arial/Helvetica are near-universal.

## GPU and CPU Backends

```python
import smsd

print(smsd.gpu_is_available())
print(smsd.gpu_device_info())
```

GPU support is opportunistic. If a GPU backend is unavailable, the same API
continues on CPU with no code changes.

## Third-party Interop

SMSD does not require RDKit or CDK. RDKit remains optional for:
- drawing and notebook depiction
- cross-validation
- mixed-toolkit workflows

Open Babel is also supported as an optional interop layer for users who prefer
that stack:

```python
import smsd
from openbabel import pybel

obmol = pybel.readstring("smi", "c1ccc(O)cc1")
g = smsd.from_openbabel(obmol)
back = smsd.to_openbabel(g)
```

When needed:

```python
import smsd
from rdkit import Chem

rdmol = Chem.MolFromSmiles("c1ccccc1")
g = smsd.from_rdkit(rdmol)
back = smsd.to_rdkit(g)
```

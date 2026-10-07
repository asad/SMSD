# SMSD 7.2.2 Python examples

Worked examples for substructure search, MCS, fingerprints, drawing and molecular
I/O. Each Python block runs independently after installing SMSD. Examples that
write files use the current directory. See the [Python guide](PYTHON.md) for
CPython 3.14 wheel platforms and matching options.

Java and C++ examples are in the [Java guide](JAVA.md) and [C++ guide](CPP.md).

## Contents

- [1. Quick start](#1-quick-start)
- [2. Substructure search](#2-substructure-search)
- [3. Maximum common substructure](#3-maximum-common-substructure)
- [4. MCS variants](#4-mcs-variants)
- [5. Tautomer and solvent settings](#5-tautomer-and-solvent-settings)
- [6. Fingerprints and similarity](#6-fingerprints-and-similarity)
- [7. SVG drawing](#7-svg-drawing)
- [8. Coordinates and layout](#8-coordinates-and-layout)
- [9. Stereo assignment](#9-stereo-assignment)
- [10. MOL and SDF files](#10-mol-and-sdf-files)
- [11. R-group decomposition](#11-r-group-decomposition)
- [12. Related molecules and constrained batches](#12-related-molecules-and-constrained-batches)
- [13. Batch operations](#13-batch-operations)
- [14. SMARTS](#14-smarts)
- [15. Scaffolds](#15-scaffolds)
- [16. Graph utilities](#16-graph-utilities)
- [17. Repeated searches](#17-repeated-searches)
- [18. Atom maps](#18-atom-maps)

## 1. Quick start

```python
import smsd

query = "c1ccccc1"       # benzene
target = "c1ccc(O)cc1"   # phenol
assert smsd.is_substructure(query, target)
mapping = smsd.find_mcs(query, target, timeout_ms=1000)
assert len(mapping) == 6
print(f"MCS: {len(mapping)} atoms")  # MCS: 6 atoms
```

Mappings use zero-based input atom indices. SMILES use parser order; RDKit
inputs use their original atom order. Exact pairs can vary for symmetric
molecules even when the chemical result is equivalent.

## 2. Substructure search

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(N)cc1")   # aniline
assert smsd.is_substructure(query, target)
mapping = smsd.find_substructure(query, target)
assert len(mapping) == 6
embeddings = smsd.find_substructure(query, target, max_results=10)
assert embeddings and all(len(m) == 6 for m in embeddings)

# Require strict aromaticity for benzene versus saturated cyclohexane.
strict = smsd.ChemOptions.profile("strict")
assert not smsd._smsd.is_substructure(query, smsd.parse_smiles("C1CCCCC1"), strict, 1000)
```

`max_results > 1` returns a list of mappings. Aromaticity and bond order are
separate from ring membership: benzene and cyclohexane are both rings. The default
aromaticity policy is flexible; the native call above uses an explicit strict profile.
Ring-only matching is off by default; enable it explicitly for MCS when needed:

```python
import smsd

chain = "CCC"
ring = "C1CC1"
assert len(smsd.find_mcs(chain, ring, timeout_ms=1000)) == 3
assert not smsd.find_mcs(chain, ring, ring_matches_ring_only=True, timeout_ms=1000)
```

## 3. Maximum common substructure

```python
import smsd

query = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
target = smsd.parse_smiles("Oc1ccccc1C(=O)O")         # salicylic acid
mapping = smsd.find_mcs(query, target, timeout_ms=1000)
assert len(mapping) == 10
assert smsd.validate_mapping(query, target, mapping) == []
fragment = smsd.mcs_to_smiles(query, mapping)
assert len(smsd.parse_smiles(fragment)) == len(mapping)
print(fragment)
```

A timeout bounds search effort, not every preparation step. The returned mapping
has no cancellation flag or optimality certificate; a timed search can return a
smaller valid result. `validate_mapping()` checks atom and bond compatibility,
not optimality or every search-objective restriction.

### Structured result

```python
import smsd

result = smsd.mcs_result("c1ccccc1", "c1ccc(O)cc1", timeout_ms=1000)
assert result.size == 6
assert result.overlap == 1.0
assert abs(result.tanimoto - 6 / 7) < 1e-12
print(result.mcs_smiles)
```

Overlap is `size / min(query_atoms, target_atoms)`; Tanimoto is
`size / (query_atoms + target_atoms - size)`. Choose the measure that answers
your question instead of comparing atom counts alone.

## 4. MCS variants

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(O)cc1")
connected = smsd.find_mcs(query, target, timeout_ms=1000)
disconnected = smsd.find_mcs(query, target, connected_only=False, timeout_ms=1000)
induced = smsd.find_mcs(query, target, induced=True, timeout_ms=1000)
edge_objective = smsd.find_mcs(query, target, maximize_bonds=True, timeout_ms=1000)
assert all(len(m) == 6 for m in (connected, disconnected, induced, edge_objective))

# Up to five MCS mappings, rather than a ranking of five different sizes.
mappings = smsd.find_mcs(query, target, max_results=5, timeout_ms=1000)
assert 1 <= len(mappings) <= 5 and all(len(m) == 6 for m in mappings)
```

Connected MCS retains one query fragment. Disconnected MCS can retain several.
Induced matching also preserves nonedges between mapped atoms.
`maximize_bonds=True` ranks matched bonds before atom count.

### Multiple molecules

```python
import smsd

molecules = [smsd.parse_smiles(s) for s in
             ["c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(Cl)cc1"]]
mapping = smsd.find_nmcs(molecules, threshold=1.0, timeout_ms=1000)
assert len(mapping) == 6
```

`find_nmcs()` uses sequential pairwise reduction. Its mapping relates the
smallest molecule's atom indices to common-fragment positions. This workflow
is order-dependent and does not prove a globally optimal multi-molecule result.

## 5. Tautomer and solvent settings

Tautomer matching and bond-order matching are separate choices. For this
keto/enol pair, loose bond order permits all four heavy atoms to match:

```python
import smsd

keto = "CC(=O)C"
enol = "CC(O)=C"
strict = smsd.find_mcs(keto, enol, timeout_ms=1000)
loose = smsd.find_mcs(keto, enol, tautomer_aware=True,
                      match_bond_order="loose", timeout_ms=1000)
assert len(strict) == 2 and len(loose) == 4
```

For native functions, solvent and pH settings belong to `ChemOptions`:

```python
import smsd

chem = smsd.ChemOptions.tautomer_profile().with_solvent(smsd.Solvent.DMSO)
chem.pH = 7.0
chem.match_bond_order = smsd.BondOrderMode.LOOSE
search = smsd.MCSOptions()
search.timeout_ms = 1000
fragment = smsd.find_mcs_smiles(smsd.parse_smiles("CC(=O)C"),
                               smsd.parse_smiles("CC(O)=C"), chem=chem, opts=search)
assert len(smsd.parse_smiles(fragment)) == 4
```

These settings control matching; the result is not a calculated equilibrium
population or a guarantee that every tautomer is equivalent.

## 6. Fingerprints and similarity

### Circular and torsion fingerprints

```python
import smsd

mol = smsd.parse_smiles("c1ccc(O)cc1")
fp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048)  # ECFP4
counts = smsd.circular_fingerprint_counts(mol, radius=2, fp_size=2048)
fcfp = smsd.circular_fingerprint(mol, radius=2, fp_size=2048, mode="fcfp")
torsion = smsd.topological_torsion(mol, fp_size=2048)
assert fp and counts and fcfp and torsion
```

ECFP describes structural environments; FCFP describes functional classes
(Rogers and Hahn, 2010). Binary fingerprints contain set bit positions; count
fingerprints contain `(position, count)` pairs. Validate the representation and
thresholds on your own task; no fingerprint is universally best for similarity
or machine learning.

### Binary and count metrics

```python
import smsd

mol1 = smsd.parse_smiles("CCO")
mol2 = smsd.parse_smiles("CCCO")
fp1 = smsd.circular_fingerprint(mol1, radius=2)
fp2 = smsd.circular_fingerprint(mol2, radius=2)
scores = [smsd.tanimoto_coefficient(fp1, fp2), smsd.dice(fp1, fp2),
          smsd.cosine(fp1, fp2)]
assert all(0.0 <= score <= 1.0 for score in scores)

c1 = smsd.circular_fingerprint_counts(mol1, radius=2)
c2 = smsd.circular_fingerprint_counts(mol2, radius=2)
count_tanimoto = smsd.count_tanimoto_coefficient(c1, c2)
dense1 = smsd.counts_to_array(c1, 2048)
dense2 = smsd.counts_to_array(c2, 2048)
count_dice = smsd.count_dice(dense1, dense2)
assert 0.0 <= count_tanimoto <= 1.0 and 0.0 <= count_dice <= 1.0

# Hexadecimal storage for a binary fingerprint.
hex_string = smsd.to_hex(fp1, fp_size=2048)
assert set(smsd.from_hex(hex_string)) == set(fp1)
```

`count_dice()` and `count_cosine()` require dense vectors. The sparse count
helpers are `count_overlap_coefficient()` and `count_tanimoto_coefficient()`.
Binary `tanimoto_coefficient()` compares set positions, even for sparse input;
it does not retain feature multiplicities.

## 7. SVG drawing

The native SVG renderer works without RDKit. Check the drawing before using it
in a report; layout and font rendering depend on the molecule and viewer.

```python
import xml.etree.ElementTree as ET
import smsd

svg = smsd.depict_svg("CC(=O)Oc1ccccc1C(=O)O")
assert ET.fromstring(svg).tag.endswith("svg")
smsd.save_svg(svg, "aspirin.svg")

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(O)cc1")
mapping = smsd.find_mcs(query, target, timeout_ms=1000)
smsd.save_svg(smsd.depict_pair(query, target, mapping), "mcs_pair.svg")

# Highlight target atoms: keys must be indices in the rendered molecule.
target_mapping = {b: a for a, b in mapping.items()}
smsd.save_svg(smsd.depict_mapping(target, target_mapping), "target.svg")
```

### Styling

```python
import smsd

options = smsd.DepictOptions()
options.bond_length = 35
options.show_map_numbers = False
options.font_family = "Arial, sans-serif"
for index, smiles in enumerate(["CCO", "c1ccccc1", "CC(=O)O"]):
    smsd.save_svg(smsd.depict_svg(smiles, opts=options), f"molecule_{index}.svg")

svg = smsd.depict_svg("c1ccc2c(c1)cc1ccccc1c2", bond_length=50,
                      width=800, height=400, padding=40, show_atom_indices=True)
smsd.save_svg(svg, "large.svg")
```

### Optional SVG conversion

After creating `aspirin.svg` above, choose an external converter. These commands
require Inkscape or CairoSVG to be installed separately:

```bash
inkscape aspirin.svg --export-type=png --export-filename=aspirin.png --export-dpi=600
python -c "import cairosvg; cairosvg.svg2png(url='aspirin.svg', write_to='aspirin.png', dpi=600)"
```

## 8. Coordinates and layout

### Two and three dimensions

```python
import smsd

mol = smsd.parse_smiles("c1ccccc1")
coords = smsd.generate_coords_2d(mol, target_bond_length=1.5)
assert len(coords) == len(mol) and all(len(point) == 2 for point in coords)
quality = smsd.layout_quality(mol, coords)
assert quality >= 0.0
coords_3d = smsd.generate_coords_3d(mol, target_bond_length=1.5)
assert len(coords_3d) == len(mol) and all(len(point) == 3 for point in coords_3d)
```

Three-dimensional coordinates provide a starting geometry, not a validated
conformer ensemble or energy minimum.

### Refinement

```python
import smsd

mol = smsd.parse_smiles("c1ccccc1")
coords = smsd.generate_coords_2d(mol)
_, coords = smsd.force_directed_layout(mol, coords, max_iter=100)
_, coords = smsd.stress_majorisation(mol, coords, max_iter=100)
crossings, coords = smsd.reduce_crossings(mol, coords, max_iter=100)
assert len(coords) == len(mol)
```

Compare layouts with `layout_quality()` and inspect the drawing. Refinement
does not guarantee a better result for every molecule.

### Transform coordinates

```python
import math
import smsd

mol = smsd.parse_smiles("c1ccccc1")
coords = smsd.generate_coords_2d(mol)
coords = smsd.translate_2d(coords, dx=10.0, dy=5.0)
coords = smsd.rotate_2d(coords, angle=math.radians(45))
coords = smsd.scale_2d(coords, factor=2.0)
coords = smsd.mirror_x(coords)
coords = smsd.mirror_y(coords)
coords = smsd.center_2d(coords)
coords = smsd.normalise_bond_length(mol, coords, target=1.5)
coords = smsd.canonical_orientation(mol, coords)
reference = smsd.generate_coords_2d(mol)
rmsd, aligned = smsd.align_2d(coords, reference)
box = smsd.bounding_box_2d(aligned)
assert len(aligned) == len(mol) and len(box) == 4 and rmsd >= 0.0
```

Rotation angles are in radians. Coordinate generation and refinement can change
with algorithms or settings; save coordinates when exact reproducibility matters.

## 9. Stereo assignment

```python
import smsd

alanine = smsd.parse_smiles("N[C@@H](C)C(=O)O")
assert smsd.assign_rs(alanine) == {1: "S"}
butene = smsd.parse_smiles("C/C=C/C")
assert smsd.assign_ez(butene) == {(1, 2): "E"}
result = smsd.assign_cip(alanine)
```

Include explicit SMILES stereo markers when the assignment should reflect
input stereochemistry. Unspecified stereo is not a resolved configuration.

## 10. MOL and SDF files

These examples create their own input files:

```python
import smsd

mol = smsd.parse_smiles("CCO")
smsd.write_molfile(mol, "molecule.mol")
loaded = smsd.read_mol_file("molecule.mol")
assert len(loaded) == 3
block = smsd.write_mol_block(loaded)
block_v3000 = smsd.write_mol_block_v3000(loaded)
assert len(smsd.read_mol_block(block)) == 3
assert len(smsd.read_mol_block(block_v3000)) == 3
```

```python
import smsd

molecules = [smsd.parse_smiles(s) for s in ["CCO", "c1ccccc1"]]
smsd.write_sdf(molecules, "compounds.sdf")
loaded = smsd.read_sdf("compounds.sdf")
assert [len(mol) for mol in loaded] == [3, 6]
smsd.export_sdf(loaded, "output.sdf")
```

`read_sdf()` loads the whole file; malformed records may produce empty graphs.
Check each record before use. For large files, iterate records and pass each to
`read_mol_block()`. V2000 supports at most 999 atoms; use V3000 for larger graphs.

## 11. R-group decomposition

```python
import smsd

core = "c1ccccc1"
targets = ["c1ccc(O)cc1", "c1ccc(N)cc1", "CCO"]
results = smsd.decompose_r_groups(core, targets, timeout_ms=1000)
assert len(results) == len(targets)
assert len(results[0]["core"]) == 6 and len(results[1]["core"]) == 6
assert results[2] == {}
```

Results retain target order. A missing core produces an empty dictionary.
R-group atom lists describe attached substituents; they are not a complete
chemical standardisation or reaction-mapping workflow.

## 12. Related molecules and constrained batches

```python
import smsd

parent = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
child = smsd.parse_smiles("Oc1ccccc1C(=O)O")
mapping = smsd.find_mcs(parent, child, timeout_ms=1000)
assert len(mapping) == 10

assigned = smsd.batch_mcs_constrained(["NO"], ["C", "CNO"],
                                      return_target_indices=True, timeout_ms=1000)
target_index, atom_mapping = assigned[0]
assert target_index == 1 and len(atom_mapping) == 2
```

The constrained batch assigns queries greedily without overlapping target atoms.
It does not prove a globally optimal reaction mapping or account for reaction
mechanisms. Preserve the returned target index when annotating several products.

## 13. Batch operations

```python
import smsd

library = [smsd.parse_smiles(s) for s in
           ["c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "CCO"]]
query = smsd.parse_smiles("c1ccccc1")
assert smsd.batch_substructure(query, library, num_threads=2) == [True, True, True, False]
mappings = smsd.batch_mcs(query, library, timeout_ms=1000, num_threads=2)
sizes = smsd.batch_mcs_size(query, library, timeout_ms=1000, num_threads=2)
assert sizes == [len(mapping) for mapping in mappings]
```

Each query/target pair receives its own time budget. Batch calls prepare graphs
before starting workers. Before separate concurrent calls share graphs, prewarm
them and avoid mutation while matching is running.

## 14. SMARTS

```python
import smsd

pattern = smsd.compile_smarts("[#6]-[#8]")
targets = [smsd.parse_smiles(s) for s in ["CCO", "CCN"]]
assert pattern.matches_many(targets) == [True, False]
phenol = smsd.parse_smiles("c1ccc(O)cc1")
assert smsd.smarts_match("[OH]", phenol)
embeddings = smsd.smarts_find_all("[#6]~[#6]", phenol, max_matches=50)
assert embeddings
assert not smsd.find_mcs_smarts("[#6]~[#7]", phenol)
```

SMARTS results are embeddings; their count is not necessarily a count of unique
bonds or functional groups. `find_mcs_smarts()` returns a complete pattern
embedding or an empty mapping. Check dialect-specific patterns on your inputs.

## 15. Scaffolds

```python
import smsd

mol = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
scaffold = smsd.murcko_scaffold(mol)
assert len(scaffold) == 6
mapping = smsd.find_scaffold_mcs(smsd.parse_smiles("c1ccc2c(c1)cccc2"),
                                smsd.parse_smiles("c1ccc2c(c1)cc1ccccc1c2"))
assert mapping
```

Scaffold extraction removes side chains. Scaffold-MCS mappings use the extracted
scaffold indices, not the original molecule indices.

## 16. Graph utilities

```python
import smsd

mol = smsd.parse_smiles("c1ccccc1.CCO")
assert smsd.count_components(mol) == 2
components = smsd.split_components(mol)
assert sorted(len(part) for part in components) == [3, 6]
assert smsd.same_canonical_graph(smsd.parse_smiles("c1ccccc1"),
                                 smsd.parse_smiles("C1=CC=CC=C1"))
```

Canonical graph comparison concerns structure. Use the requested chemistry and
stereo settings when matching stereoisomers.

## 17. Repeated searches

```python
import smsd

query = smsd.parse_smiles("c1ccccc1")
targets = [smsd.parse_smiles(s) for s in ["c1ccc(O)cc1", "CCO"]]
for graph in [query, *targets]:
    smsd.prewarm_graph(graph)
selected = smsd.screen_targets(query, targets, threshold=0.5)
assert selected == [0]
mappings = [smsd.find_mcs(query, targets[index], timeout_ms=1000) for index in selected]
assert len(mappings[0]) == 6
print(smsd.gpu_device_info())
```

Prewarming helps repeated work; it is not required before every batch call.
`similarity_upper_bound()` and `screen_targets()` provide RASCAL screening bounds.
Release wheels use CPU/OpenMP; optional GPU screening requires a source build
and compatible hardware. Benchmark your own molecules and settings.

## 18. Atom maps

```python
import smsd

mapped = "[CH3:1][C:2](=[O:3])[OH:4]"
clean = smsd.strip_atom_maps(mapped)
assert clean == "[CH3][C](=[O])[OH]"
assert len(smsd.parse_smiles(clean)) == 4
```

Use `strip_atom_maps()` rather than a regular expression. Colons can also belong
to aromatic ring closures, so a text-only replacement can alter the structure.

## Further reading

- [Python guide](PYTHON.md): options, atom indices and interoperability.
- [Validation](VALIDATION_7.2.2.md): platform tests and their scope.
- [Benchmark report](../benchmarks/RESULTS_7.2.0.md): measured results for its recorded versions.
- [Citation information](../CITATION.cff), [LICENSE](../LICENSE) and [NOTICE](../NOTICE).

Copyright (c) 2009-2026 Syed Asad Rahman, BioInception PVT LTD. Apache-2.0.

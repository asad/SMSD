# SMSD Pro Python examples

These runnable examples cover the `smsd` Python package beyond molecule search:
fingerprints, drawing, layout, stereo, MOL and SDF files, R-groups, scaffolds,
graph utilities and atom maps. Each block is self-contained; blocks that write
files create them in the current directory. Substructure search, maximum common
substructure (MCS) search, chemistry options, timeouts and the full API
reference are in the [Python guide](PYTHON.md).

## Fingerprints and similarity

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

ECFP describes structural environments and FCFP describes functional classes
([Rogers and Hahn, 2010](https://doi.org/10.1021/ci100050t)). Binary
fingerprints contain set bit positions; count fingerprints contain
`(position, count)` pairs. No fingerprint suits every similarity or
machine-learning task, so validate the representation and thresholds on your own data.

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

`count_dice()` and `count_cosine()` require dense vectors; the sparse count
helpers are `count_tanimoto_coefficient()` and `count_overlap_coefficient()`.
Binary `tanimoto_coefficient()` compares set positions, even for count input,
and does not retain feature multiplicities. Two empty fingerprints score 1.0
with `tanimoto_coefficient()` and `overlap_coefficient()`, as identical empty
sets, in C++, Java and Python.

## SVG drawing

The native SVG renderer does not need RDKit. Check each drawing before using it
in a report, as layout and font rendering depend on the molecule and the viewer.

```python
import xml.etree.ElementTree as ET
import smsd

svg = smsd.depict_svg("CC(=O)Oc1ccccc1C(=O)O")
assert ET.fromstring(svg).tag.endswith("svg")
smsd.save_svg(svg, "aspirin.svg")

query = smsd.parse_smiles("c1ccccc1")
target = smsd.parse_smiles("c1ccc(O)cc1")
result = smsd.MCS.find(query).with_timeout(1000).match(target)
mapping = dict(result.best.atom_mapping)
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

Keyword arguments such as `bond_length=50` override `DepictOptions`
attributes; when `opts` is also passed, they are set on that object.

### Optional conversion to PNG

After creating `aspirin.svg` as above, use an external converter. Inkscape and
CairoSVG are installed separately:

```bash
inkscape aspirin.svg --export-type=png --export-filename=aspirin.png --export-dpi=600
python -c "import cairosvg; cairosvg.svg2png(url='aspirin.svg', write_to='aspirin.png', dpi=600)"
```

## Coordinates and layout

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

A `layout_quality()` score of 0.0 is ideal; lower is better. Three-dimensional
coordinates provide a starting geometry, not a validated conformer ensemble or
energy minimum.

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

Refinement does not improve every layout. Compare the results with
`layout_quality()` and inspect the drawing.

### Transforms

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

Rotation angles are in radians. Generated and refined coordinates can change
with the algorithm or its settings, so save them when exact reproducibility
matters.

## Stereo assignment

```python
import smsd

alanine = smsd.parse_smiles("N[C@@H](C)C(=O)O")
assert smsd.assign_rs(alanine) == {1: "S"}
butene = smsd.parse_smiles("C/C=C/C")
assert smsd.assign_ez(butene) == {(1, 2): "E"}
descriptors = smsd.assign_cip(alanine)
assert descriptors.rs_labels[1] == smsd.RSLabel.S
```

Write explicit stereo markers in the SMILES when the assignment should reflect
the input stereochemistry. Unspecified stereo is not a resolved configuration.

## MOL and SDF files

Each example creates its own input files.

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
smsd.write_sdf(loaded, "output.sdf")
```

`read_sdf()` loads the whole file, and malformed records become empty graphs,
so check each record before use. For large files, split the records and pass
each to `read_mol_block()`. The V2000 format holds at most 999 atoms and 999
bonds, so `write_mol_block()` and the SDF writers switch to V3000 for larger
graphs.

## R-group decomposition

```python
import smsd

options = smsd.Substructure.Options()
options.timeout_ms = 1000
results = smsd.Substructure.decompose_r_groups("c1ccccc1", ["Oc1ccccc1", "Nc1ccccc1", "CCO"], options)
assert [result["search"].found for result in results] == [True, True, False]
assert results[0]["core"].n == 6 and results[0]["rgroups"]["R1"].n == 1
```

Each entry holds the core graph, the substituent graphs and the search result
for the core, which shows whether a missing core was proved absent.

## Scaffolds

```python
import smsd

parent = smsd.parse_smiles("Cc1ccccc1")
scaffold = smsd.murcko_scaffold(parent)
assert scaffold.n == 6
assert smsd.Substructure.find(scaffold).match(parent).found
result = smsd.MCS.find(scaffold).match("Oc1ccccc1")
assert result.best.atom_count == 6
```

The scaffold retains ring systems, the paths that connect them and every atom
joined to either by a double or triple bond, so that carbonyl groups stay as in
the Bemis and Murcko scaffold. An acyclic input is returned unchanged rather
than as an empty scaffold.

## Graph utilities

```python
import smsd

mol = smsd.parse_smiles("c1ccccc1.CCO")
assert smsd.count_components(mol) == 2
components = smsd.split_components(mol)
assert sorted(len(part) for part in components) == [3, 6]
assert smsd.same_canonical_graph(smsd.parse_smiles("c1ccccc1"),
                                 smsd.parse_smiles("C1=CC=CC=C1"))
```

Canonical graph comparison is a structural check. To compare stereoisomers in a
search, enable the stereo options described in the
[Python guide](PYTHON.md#chemistry-options).

## Atom maps

```python
import smsd

mapped = "[CH3:1][C:2](=[O:3])[OH:4]"
clean = smsd.strip_atom_maps(mapped)
assert clean == "[CH3][C](=[O])[OH]"
assert len(smsd.parse_smiles(clean)) == 4
```

Use `strip_atom_maps()` rather than a regular expression. Outside brackets, a
colon is an aromatic bond and can precede a ring-closure digit, so a text-only
replacement can alter the structure. `parse_mapped_smiles()` strips the maps and
parses in one call.

## Further reading

See the [Python guide](PYTHON.md), the [installation guide](HOWTO-INSTALL.md),
the [C++ guide](CPP.md) and the [Java guide](JAVA.md).

# SMSD Pro Python examples

Runnable examples for the `smsd` Python package beyond molecule search. Each
block is self-contained, and blocks that write files create them in the current
directory. Searches, options and the full API reference are in the
[Python guide](PYTHON.md). Interpreter and platform requirements are in the
[supported-platforms table](HOWTO-INSTALL.md#supported-platforms).

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

ECFP describes structural environments and FCFP functional classes
([Rogers and Hahn, 2010](https://doi.org/10.1021/ci100050t)). Binary
fingerprints hold set bit positions; count fingerprints hold
`(position, count)` pairs. Choose the fingerprint and threshold on your own
data.

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

`count_dice()` and `count_cosine()` accept dense vectors, `(position, count)`
pairs or dictionaries; `count_tanimoto_coefficient()` and
`count_overlap_coefficient()` take pairs. Binary `tanimoto_coefficient()`
compares set positions only, even for count input. Two empty fingerprints
score 1.0 with `tanimoto_coefficient()` and `overlap_coefficient()`.

### Configurable fingerprints

```python
import smsd

aspirin = "CC(=O)Oc1ccccc1C(=O)O"
# ECFP4 with ring-size and stereo terms; the default ECFP options give the
# same bits as circular_fingerprint().
options = smsd.FingerprintOptions(family="ecfp", ring_size=True, cip_atoms=True, cip_bonds=True)
fp = smsd.fingerprint(aspirin, options=options)
keys = smsd.fingerprint(aspirin, options=smsd.FingerprintOptions.maccs())  # 166 MACCS keys
assert fp and keys
```

The same `FingerprintOptions` give the same bits in C++ and Java.

## Screening

### Rank a small set by similarity

```python
import smsd

query = "CC(=O)Oc1ccccc1C(=O)O"                    # aspirin
library = ["OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1", "c1ccccc1", "CC(=O)Oc1ccccc1"]
ranked = smsd.batch_similarity(query, library)    # (index, Tanimoto), most similar first
assert library[ranked[0][0]] == "CC(=O)Oc1ccccc1"
for index, score in ranked:
    print(f"{score:.2f}  {library[index]}")
```

### Screen a large library, then confirm by substructure

```python
import itertools
import smsd

groups = ["C", "CC", "OC", "N", "O", "F", "Cl", "Br", "C#N", "C(F)(F)F",
          "C(=O)O", "C(=O)N", "S(=O)(=O)N", "OCC"]
library = [f"c1c({a})cc({b})cc1{c}" for a, b, c in itertools.product(groups, repeat=3)]
assert len(library) == 2744

corpus = smsd.TargetCorpus.from_smiles(library)   # parse once, search many times
corpus.prewarm()

query = "c1c(C)cc(O)cc1C(=O)O"
nearest = smsd.batch_similarity(query, library)[:10]

acid = smsd.Substructure.smarts("C(=O)[OX2H1]")
results = acid.match_many(corpus, num_threads=0)  # every processor core
hits = [library[i] for i, result in enumerate(results) if result.found]
assert library[nearest[0][0]] == query and query in hits
```

With OpenMP, as in the published wheels, `num_threads=0` uses every core unless
`OMP_NUM_THREADS` sets fewer; results keep the order of the targets.

### Screen by MCS similarity, then compute the MCS

```python
import itertools
import smsd

groups = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = [f"c1c({a})cc({b})cc1{c}" for a, b, c in itertools.product(groups, repeat=3)]

corpus = smsd.TargetCorpus.from_smiles(library)   # prepared once for repeated screens

query = "c1c(C)cc(O)cc1C(=O)O"
candidates = smsd.screen_targets(query, corpus, 0.8)   # targets that can reach 0.8
results = smsd.MCS.find(query).with_timeout(1000).match_many([library[i] for i in candidates])
assert candidates and all(result.found for result in results)
```

`screen_targets` bounds the atom-count Tanimoto similarity of the MCS from
above, so under the same chemistry options it never rejects a target able to
reach the threshold, and the MCS runs only where it can matter. For repeated
screens of the same library, a prepared `TargetCorpus` such as this one saves
converting the targets on every call; a list of molecules or SMILES gives the
same indices.

This threshold applies to atom-count Tanimoto, not bond Tanimoto or weighted MCS scores.
`MCS.match_many` itself searches every supplied target without a similarity threshold.

With `smsd-cuda` and an NVIDIA GPU,
[`gpu_screening.py`](../python/examples/gpu_screening.py) runs packed and
corpus screens on the device and checks every result against the CPU; without
them it runs on the CPU.

## SVG drawing

The native SVG renderer needs no RDKit. Check each drawing before publishing
it, as fonts depend on the viewer.

```python
import xml.etree.ElementTree as ET
import smsd

svg = smsd.depict_svg("CC(=O)Oc1ccccc1C(=O)O")
assert ET.fromstring(svg).tag.endswith("svg")
smsd.save_svg(svg, "aspirin.svg")

query = smsd.parse_smiles("CC(=O)NC")   # N-methylacetamide
target = smsd.parse_smiles("CCC(=O)N")  # propanamide
result = smsd.MCS.find(query).with_timeout(1000).match(target)
mapping = dict(result.best.atom_mapping)
assert mapping == {0: 1, 1: 2, 2: 3, 3: 4}
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

Keyword arguments such as `bond_length=50` set `DepictOptions` attributes, on
`opts` when it is passed.

### Optional conversion to PNG

Convert `aspirin.svg` from above with Inkscape or CairoSVG, installed
separately:

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

A `layout_quality()` score of 0.0 is ideal. Three-dimensional coordinates are
a starting geometry, not a validated conformer or energy minimum.

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

Refinement does not improve every layout; compare `layout_quality()` and
inspect the drawing.

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

Angles are in radians. Save coordinates when exact reproducibility matters.

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

Write explicit stereo markers in the SMILES; unspecified stereo is not
assigned.

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

`read_sdf()` loads the whole file and turns malformed records into empty
graphs, so check each record. For large files, pass records one at a time to
`read_mol_block()`. `write_mol_block()` and the SDF writers switch to V3000
above 999 atoms or bonds.

## R-group decomposition

```python
import smsd

options = smsd.Substructure.Options()
options.timeout_ms = 1000
results = smsd.Substructure.decompose_r_groups("c1ccccc1", ["Oc1ccccc1", "Nc1ccccc1", "CCO"], options)
assert [result["search"].found for result in results] == [True, True, False]
assert results[0]["core"].n == 6 and results[0]["rgroups"]["R1"].n == 1
```

Each entry holds the core graph, the substituent graphs and the core's search
result, which shows whether a missing core was proved absent.

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

The Bemis–Murcko scaffold keeps ring systems, their linkers and atoms double-
or triple-bonded to them, such as carbonyl oxygens. An acyclic input is
returned unchanged.

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

To compare stereoisomers in a search, enable the stereo
[chemistry options](PYTHON.md#chemistry-options).

## Atom maps

```python
import smsd

mapped = "[CH3:1][C:2](=[O:3])[OH:4]"
clean = smsd.strip_atom_maps(mapped)
assert clean == "[CH3][C](=[O])[OH]"
assert len(smsd.parse_smiles(clean)) == 4
```

Use `strip_atom_maps()` rather than a regular expression: outside brackets, a
colon is an aromatic bond, so text replacement can alter the structure.
`parse_mapped_smiles()` strips and parses in one call.

## AI assistants

An AI assistant can call the molecular tools through the local MCP server
(`smsd-mcp`) or through tool definitions in your own Python loop. The tools run
on your machine. Each search result reports its `status` and whether it is
`optimal` and `complete`. See [Connecting AI assistants](AI_ASSISTANTS.md) for
set-up steps and worked examples.
[`agent_workflow.py`](../python/examples/agent_workflow.py) chains screening,
MCS and depiction through `call_tool` without a model, and
[`mcp_stdio_client.py`](../python/examples/mcp_stdio_client.py) checks the MCP
server from Python.

## Further reading

See the [Python guide](PYTHON.md), the [installation guide](HOWTO-INSTALL.md),
the [C++ guide](CPP.md) and the [Java guide](JAVA.md). A C++ program that
screens packed fingerprints on a Metal GPU, with the CPU screen as the fallback,
is the [Metal screening example](../cpp/examples/metal_fingerprint_screen.cpp).

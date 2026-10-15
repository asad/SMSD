# SMSD Python examples

These standalone examples use the current source API. Release 7.2.2 predates
this API; see [building from source](PYTHON.md#build-from-source).

## 1. Quick start

```python
import smsd

query, target = "c1ccccc1", "c1ccc(O)cc1"
assert smsd.Substructure.find(query).match(target).found
result = smsd.MCS.find(query).with_timeout(1000).match(target)
assert result.optimal and result.best.atom_count == 6
print(result.best.atom_mapping)
```

## 2. Substructure search

```python
import smsd

query = smsd.Substructure.find("CO")
assert query.match("CCO").found
embeddings = query.match_all("OCO", limit=10)
assert len(embeddings.matches) == 2 and embeddings.complete
assert not smsd.Substructure.find("CCC").with_induced().match("C1CC1").found
```

`match` and `match_all` return a `SearchResult`. Empty results prove absence
only when the search completed. Molecule queries use SMILES or graph inputs;
SMARTS predicate queries use `Substructure.smarts`.

## 3. Maximum common substructure

```python
import smsd

query, target = "c1ccc(O)cc1", "c1ccc(N)cc1"
result = smsd.MCS.find(query).match(target)
assert result.best.atom_count == 6
assert result.best.bond_count == 6
assert smsd.MCS.validate_mapping(query, target, result.best) == []
fragment = result.best.to_smiles(query)
assert smsd.parse_smiles(fragment).n == 6
```

Results retain selected bonds, fragment membership, completion status and an
optimality flag. A stopped search may return a valid incumbent without a proof.
Mappings use each input's original atom order.

## 4. MCIS, MCES and fragments

```python
import smsd

query, target = "C1CC1", "CCC"
induced = smsd.MCS.find(query).match(target)
edges = smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES).match(target)
assert (induced.best.atom_count, induced.best.bond_count) == (2, 1)
assert (edges.best.atom_count, edges.best.bond_count) == (3, 2)
assert len(smsd.parse_smiles(edges.best.to_smiles(query)).bonds()) == 2
fragments = smsd.MCS.find("CO.CN").with_fragments().match("CO.CN")
assert fragments.best.atom_count == 4 and len(fragments.best.fragments) == 2
all_matches = smsd.MCS.find(query).match_all(target, limit=100)
assert all_matches.optimal and all_matches.complete
```

MCIS maximises atoms and preserves edges and nonedges. MCES maximises selected
common edges and then incident atoms, allowing input edges to be omitted.
Connectivity and fragment limits apply to selected common edges.

```python
import smsd

result = smsd.MCS.consensus(["CCO", "CCN", "CCCl"], threshold=1.0)
assert result.found and result.support_count == 3
assert result.graph.n == 2 and not result.optimal
```

Consensus uses progressive reduction and verifies support across the inputs.
It does not prove the global multi-molecule optimum.

## 5. Chemistry and budgets

```python
import smsd

options = smsd.MCS.Options()
options.chemistry.tautomer_aware = True
options.chemistry.match_bond_order = smsd.BondOrderMode.ANY
options.timeout_ms = 1000
result = smsd.MCS.find("CC(=O)C").with_options(options).match("CC(O)=C")
assert result.best.atom_count == 4
options.chemistry.match_bond_order = smsd.BondOrderMode.STRICT
strict = smsd.MCS.find("CC(=O)C").with_options(options).match("CC(O)=C")
assert strict.best.atom_count == 2
```

Chemistry includes atom type, charge, isotope, aromaticity, rings, stereo and
tautomer settings. Complete-ring matching checks both inputs. A timeout of
`-1` selects an adaptive budget, zero requests an immediate stop and positive
values apply per pair. Node limits and a shared `CancellationToken` also stop
cooperatively. Check `status`, `optimal` and `complete`.

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

options = smsd.Substructure.Options()
options.timeout_ms = 1000
results = smsd.Substructure.decompose_r_groups("c1ccccc1", ["Oc1ccccc1", "Nc1ccccc1", "CCO"], options)
assert [result["search"].found for result in results] == [True, True, False]
assert results[0]["core"].n == 6 and results[0]["rgroups"]["R1"].n == 1
```

The decomposition returns core/substituent graphs and its search result,
including whether a missing core was proved absent.

## 12. Constrained mapping

```python
import smsd

results = smsd.MCS.match_constrained(["CO", "CN"], "OCCN")
assert all(result.found for result in results)
first, second = [set(result.best.atom_mapping.values()) for result in results]
assert not first.intersection(second)
```

Queries claim atoms in the original target without re-perceiving a residual
molecule. Results retain query order and use independent per-pair budgets.
The sequential assignment does not prove a globally optimal joint mapping.

## 13. Batch operations

```python
import smsd

query = smsd.Substructure.find("c1ccccc1")
targets = ["c1ccc(O)cc1", "Cc1ccccc1", "CCO"]
results = query.match_many(targets, num_threads=2)
assert [result.found for result in results] == [True, True, False]
```

Both query classes expose `match_many`, preserving target order and per-pair
status. Prepared queries can be reused concurrently.

## 14. SMARTS

```python
import smsd

pattern = smsd.Substructure.smarts("[#6]-[#8]")
assert pattern.match("CCO").found
assert not pattern.match("CCN").found
assert smsd.Substructure.smarts("[OH]").match("Oc1ccccc1").found
assert smsd.Substructure.smarts("[#6]~[#6]").match_all("c1ccccc1", limit=50).matches
```

SMARTS searches complete predicate-query embeddings. Compile a prepared
query once for repeated targets.

## 15. Scaffolds

```python
import smsd

parent = smsd.parse_smiles("Cc1ccccc1")
scaffold = smsd.murcko_scaffold(parent)
assert scaffold.n == 6
assert smsd.Substructure.find(scaffold).match(parent).found
result = smsd.MCS.find(scaffold).match("Oc1ccccc1")
assert result.best.atom_count == 6
```

Scaffold extraction retains rings and connecting paths. Acyclic inputs are
returned unchanged.

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

corpus = smsd.TargetCorpus.from_smiles(["Oc1ccccc1", "CCO"])
corpus.prewarm(num_threads=1)
results = smsd.Substructure.find("c1ccccc1").match_many(corpus, num_threads=2)
assert [result.found for result in results] == [True, False]
```

The corpus captures input graphs and atom correspondence when added. It can
precompute ECFP fingerprints for `screen(query, threshold)`. Fingerprint
screening and graph-search completion are separate results.

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

See the [Python guide](PYTHON.md), [C++ guide](CPP.md) and [Java guide](JAVA.md).

# Fingerprints and similarity

SMSD computes hashed fingerprints, substructure-key fingerprints and the
coefficients that compare them. A binary fingerprint is a sorted list of the
positions of its set bits. A count fingerprint is a list of `(position, count)`
pairs for the bins that are not zero. Both are plain Python lists, so they can
be stored, compared and turned into arrays without any further package.

This chapter starts with one molecule, then covers each family, count
fingerprints, configurable specifications, substructure keys, batch
generation, the similarity coefficients, encodings and NumPy arrays. The
Options section lists every signature with its defaults. Each Python block is a
complete program: run it as it stands with SMSD 8.1.0 and it prints the output
shown below it. The block under NumPy arrays also needs NumPy.

Compare fingerprints only when they were made with the same family, radius or
path length, size and switches. Nothing in a fingerprint records how it was
made, so a mismatch gives a meaningless score without an error.

## Fingerprint one molecule

`smsd.fingerprint` takes a SMILES string, an SMSD `MolGraph` or an RDKit
molecule and returns the sorted set positions. The `kind` argument chooses the
family. The default kind differs between entry points: `fingerprint` uses the
labelled path fingerprint `"mcs"`, while `fingerprint_from_smiles` uses ECFP.
Name the kind, or pass a `FingerprintOptions`, whenever the fingerprints will
be compared with others.

```python
import smsd

aspirin = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin

default = smsd.fingerprint(aspirin)                  # kind="mcs"
from_smiles = smsd.fingerprint_from_smiles(aspirin)  # mode="ecfp"
ecfp4 = smsd.fingerprint(aspirin, kind="ecfp")       # radius 2, 2048 bits

print("fingerprint() set bits:", len(default))
print("fingerprint_from_smiles() set bits:", len(from_smiles))
print("fingerprint_from_smiles() equals kind='ecfp':", from_smiles == ecfp4)
print("first five ECFP4 positions:", ecfp4[:5])
assert len(default) == 127 and default == smsd.fingerprint(aspirin, kind="mcs")
assert len(from_smiles) == 36 and from_smiles == ecfp4
assert ecfp4[:5] == [23, 44, 61, 101, 172]
```

The program prints:

```text
fingerprint() set bits: 127
fingerprint_from_smiles() set bits: 36
fingerprint_from_smiles() equals kind='ecfp': True
first five ECFP4 positions: [23, 44, 61, 101, 172]
```

The default fingerprint of each entry point:

| Entry point | Default fingerprint |
| --- | --- |
| `fingerprint`, `fingerprints` | `kind="mcs"` |
| `smsd.adapters.fingerprints_numpy`, `fingerprints_torch` | `kind="mcs"` |
| `fingerprint_from_smiles` | `mode="ecfp"` |
| `FingerprintCorpus.from_smiles` | `kind="ecfp"` |
| `smsd.adapters.fingerprints_scipy`, `molecules_pandas`, `molecules_polars`, `FingerprintTransformer`, `SimilarityTransformer` | `kind="ecfp"` |
| `batch_similarity` | ECFP with `radius=2` and `fp_size=2048`, both adjustable |
| `TargetCorpus.screen` and `screen_many`, `warm_start`, `stream_screen_smiles`, `stream_screen_many_smiles` | Fixed: binary ECFP, radius 2, 2,048 bits |

## Families

Each kind has a fixed function that computes it from a parsed molecule:

| Kind | What it hashes | Fixed function | Fixed function takes SMILES |
| --- | --- | --- | --- |
| `"ecfp"` (also `"circular"`) | Circular atom environments: element, degree, valence, formal charge, isotope, hydrogen count, tautomer class, aromaticity, bond order and ring membership | `circular_fingerprint(graph, radius, fp_size)` | No |
| `"fcfp"` | Circular environments of pharmacophore classes, with degree, bond order and ring membership | `circular_fingerprint(graph, radius, fp_size, mode="fcfp")` | No |
| `"path"` | Heavy-atom paths of up to `path_length` bonds, by atomic number only; bonds are not hashed | `path_fingerprint(graph, path_length, fp_size)` | No |
| `"mcs"` | Heavy-atom paths of up to `path_length` bonds, labelled by element, ring membership, aromaticity, tautomer class and degree, and by bond order, ring and aromatic flags | `mcs_fingerprint(graph, path_length, fp_size)` | No |
| `"torsion"` | Linear four-atom paths, typed by element, heavy-atom degree, pi electrons and ring membership | `topological_torsion(mol, fp_size)` | Yes |
| keys | Substructure keys, including the 166 MACCS keys; set through `FingerprintOptions` only | `fingerprint(mol, options=FingerprintOptions.maccs())` | Yes |

The functions that take a graph accept an SMSD `MolGraph` only, so pass
`smsd.parse_smiles(smiles)`; a SMILES string raises `TypeError`. `fingerprint`,
`fingerprints`, `fingerprint_from_smiles` and `topological_torsion` also take
SMILES strings and RDKit molecules. Each kind gives exactly the bits of its
fixed function, and the classic `FingerprintOptions` of the ECFP, FCFP, path
and torsion families give them again:

```python
import smsd

aspirin = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
graph = smsd.parse_smiles(aspirin)

fixed = {
    "ecfp": smsd.circular_fingerprint(graph, radius=2, fp_size=2048),
    "fcfp": smsd.circular_fingerprint(graph, radius=2, fp_size=2048, mode="fcfp"),
    "path": smsd.path_fingerprint(graph, path_length=7, fp_size=2048),
    "mcs": smsd.mcs_fingerprint(graph, path_length=7, fp_size=2048),
    "torsion": smsd.topological_torsion(graph, fp_size=2048),
}
classic = {
    "ecfp": smsd.FingerprintOptions.ecfp(),
    "fcfp": smsd.FingerprintOptions.fcfp(),
    "path": smsd.FingerprintOptions.path(),
    "torsion": smsd.FingerprintOptions.torsion(),
}
for kind, bits in fixed.items():
    assert smsd.fingerprint(aspirin, kind=kind) == bits
    if kind in classic:
        assert smsd.fingerprint(aspirin, options=classic[kind]) == bits
    print(f"{kind:8} {len(bits):3d} set bits")
maccs = smsd.fingerprint(aspirin, options=smsd.FingerprintOptions.maccs())
print(f"{'maccs':8} {len(maccs):3d} set bits")
assert [len(bits) for bits in fixed.values()] == [36, 29, 29, 127, 16]
assert len(maccs) == 21

try:
    smsd.circular_fingerprint(aspirin)
except TypeError:
    print("circular_fingerprint() needs a MolGraph")
else:
    raise AssertionError("a SMILES string was accepted")

try:
    smsd.fingerprint(aspirin, kind="maccs")
except ValueError as error:
    print(error)
    assert str(error).startswith("Unknown fingerprint kind: 'maccs'")
else:
    raise AssertionError("an unknown kind was accepted")
```

The program prints:

```text
ecfp      36 set bits
fcfp      29 set bits
path      29 set bits
mcs      127 set bits
torsion   16 set bits
maccs     21 set bits
circular_fingerprint() needs a MolGraph
Unknown fingerprint kind: 'maccs'. Use 'path', 'mcs', 'ecfp', 'circular', 'fcfp', or 'torsion'.
```

MACCS keys are not a `kind`; they come from `FingerprintOptions.maccs()`, as
shown under [Substructure keys and MACCS](#substructure-keys-and-maccs).

The FCFP pharmacophore class of each atom is available on its own.
`classify_pharmacophore(graph, atom_index)` returns a bit mask: 1 for a
hydrogen-bond donor, 2 for an acceptor, 4 for a positive ionisable atom, 8 for
a negative ionisable atom, 16 for an aromatic atom and 32 for a hydrophobic
atom.

```python
import smsd

paracetamol = smsd.parse_smiles("CC(=O)Nc1ccc(O)cc1")  # paracetamol
classes = [smsd.classify_pharmacophore(paracetamol, atom)
           for atom in range(paracetamol.n)]
print(classes)
assert classes == [32, 0, 2, 3, 16, 16, 16, 16, 3, 16, 16]
```

The program prints:

```text
[32, 0, 2, 3, 16, 16, 16, 16, 3, 16, 16]
```

## Count fingerprints

A count fingerprint records how often each folded bin is hit. Its positions
are the set positions of the binary fingerprint of the same specification.
`circular_fingerprint_counts` and `topological_torsion_counts` compute counts
from a graph; for any family, `FingerprintOptions(..., counts=True)` does the
same through `fingerprint` and `fingerprints`.

```python
import smsd

hexane = smsd.parse_smiles("CCCCCC")  # hexane
counts = smsd.circular_fingerprint_counts(hexane, radius=2, fp_size=2048)
print("ECFP4 counts:", counts)
assert len(counts) == 8 and counts[0] == (592, 2) and counts[-1] == (1861, 4)
assert [position for position, _ in counts] == smsd.circular_fingerprint(hexane)
counted = smsd.FingerprintOptions("ecfp", counts=True)
assert smsd.fingerprint("CCCCCC", options=counted) == counts

path = smsd.FingerprintOptions("path", counts=True)
path_counts = smsd.fingerprint("CCCCCC", options=path)
print("path counts:", path_counts)
assert all(count == 1 for _, count in path_counts)

aspirin = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
maccs = smsd.FingerprintOptions.maccs()
bits = smsd.fingerprint(aspirin, options=maccs)
maccs.counts = True
key_counts = dict(smsd.fingerprint(aspirin, options=maccs))
print("MACCS key 125 count:", key_counts[124], "| bit set:", 124 in bits)
assert key_counts[124] == 1 and 124 not in bits
```

The program prints:

```text
ECFP4 counts: [(592, 2), (906, 2), (1015, 2), (1147, 2), (1512, 2), (1702, 2), (1828, 2), (1861, 4)]
path counts: [(145, 1), (273, 1), (277, 1), (401, 1), (405, 1), (533, 1)]
MACCS key 125 count: 1 | bit set: False
```

Path counts are always 1. Key counts are the raw number of matches and do not
apply a key's minimum: MACCS key 125 needs two aromatic rings, so aspirin
reports a count of 1 for it while the binary bit stays clear.

## Configurable fingerprints

`FingerprintOptions(family="ecfp", **switches)` describes a fingerprint
completely. It starts from the classic configuration of the family, which gives
the bits of the fixed function, and each keyword then switches a hashed term on
or off. The class methods `ecfp`, `fcfp`, `path`, `torsion` and `maccs` build
the classic configurations. Pass the value as `options` to
`fingerprint`, `fingerprints`, `fingerprint_from_smiles`,
`FingerprintCorpus.from_smiles` or the NumPy, PyTorch, SciPy, pandas and
Polars adapters; the same options give the same bits in C++ and Java.

Stereo and ring terms are off in the classic configurations:

```python
import smsd

fcfp4 = smsd.FingerprintOptions.fcfp()
print(fcfp4)
assert repr(fcfp4).endswith("terms=[degree, bond_order, ring_membership])")

first, second = "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"  # alanine enantiomers
classic = smsd.FingerprintOptions()
stereo = smsd.FingerprintOptions(cip_atoms=True)
same = smsd.fingerprint(first, options=classic) == smsd.fingerprint(
    second, options=classic)
apart = smsd.fingerprint(first, options=stereo) != smsd.fingerprint(
    second, options=stereo)
print("classic ECFP4 separates the enantiomers:", not same)
print("cip_atoms separates them:", apart)
assert same and apart

e_butene, z_butene = "C/C=C/C", "C/C=C\\C"  # E- and Z-2-butene
geometry = smsd.FingerprintOptions(cip_bonds=True)
apart = smsd.fingerprint(e_butene, options=geometry) != smsd.fingerprint(
    z_butene, options=geometry)
print("cip_bonds separates E and Z:", apart)
assert apart

cyclohexane = "C1CCCCC1"  # cyclohexane
ring = smsd.FingerprintOptions(ring_size=True)
changed = smsd.fingerprint(cyclohexane, options=ring) != smsd.fingerprint(
    cyclohexane, options=classic)
print("ring_size changes the bits:", changed)
assert changed
```

The program prints:

```text
FingerprintOptions(family='fcfp', size=2048, radius=2, terms=[degree, bond_order, ring_membership])
classic ECFP4 separates the enantiomers: False
cip_atoms separates them: True
cip_bonds separates E and Z: True
ring_size changes the bits: True
```

The `scaffold`, `fragments` and `radius` keywords change what part of the
molecule is fingerprinted. `scaffold="murcko"` fingerprints the ring framework
with its linkers. `fragments="largest"` keeps the component with the most heavy
atoms, which removes counter-ions, and `"per_fragment"` tags each atom term
with its component. A radius of -1 grows each environment to the whole
component.

```python
import smsd

ecfp4 = smsd.FingerprintOptions()

murcko = smsd.FingerprintOptions(scaffold="murcko")
toluene = smsd.fingerprint("Cc1ccccc1", options=murcko)  # toluene
benzene = smsd.fingerprint("c1ccccc1", options=ecfp4)    # benzene
print("Murcko framework of toluene equals benzene:", toluene == benzene)
assert toluene == benzene

salt = "CC(=O)Oc1ccccc1C(=O)O.[Na+].[Cl-]"  # aspirin with sodium chloride
parent = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=ecfp4)  # aspirin
whole = smsd.fingerprint(salt, options=ecfp4)
largest = smsd.fingerprint(salt, options=smsd.FingerprintOptions(fragments="largest"))
print("whole salt equals aspirin:", whole == parent)
print("largest fragment equals aspirin:", largest == parent)
assert whole != parent and largest == parent

decane = "CCCCCCCCCC"  # decane
full = smsd.fingerprint(decane, options=smsd.FingerprintOptions(radius=-1))
print("decane set bits, radius -1 and radius 2:",
      len(full), len(smsd.fingerprint(decane, options=ecfp4)))
assert len(full) == 44 and len(smsd.fingerprint(decane, options=ecfp4)) == 9
```

The program prints:

```text
Murcko framework of toluene equals benzene: True
whole salt equals aspirin: False
largest fragment equals aspirin: True
decane set bits, radius -1 and radius 2: 44 9
```

Assigning `family` after construction keeps the current switches, so ECFP
options changed to `"fcfp"` are not classic FCFP; build
`FingerprintOptions("fcfp")` or `FingerprintOptions.fcfp()` instead. Switch
values must be `True` or `False`, and an unknown keyword raises `TypeError`.
Values that are out of range raise `ValueError` when the fingerprint is
computed. Some messages use the C++ field names, such as `ringSizeCap` for
`ring_size_cap` and `pathLength` for `path_length`.

```python
import smsd

aspirin = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
changed = smsd.FingerprintOptions.ecfp()
changed.family = "fcfp"  # keeps the ECFP switches
classic = smsd.FingerprintOptions.fcfp()
same = smsd.fingerprint(aspirin, options=changed) == smsd.fingerprint(
    aspirin, options=classic)
print("assigned family gives classic FCFP:", same)
assert not same

expected = {
    "radius=-2": "radius must be >= -1 (use -1 for unlimited)",
    "ring_size_cap=2": "ringSizeCap must be at least 3",
}
cases = {
    "radius=-2": smsd.FingerprintOptions(radius=-2),
    "ring_size_cap=2": smsd.FingerprintOptions(ring_size=True, ring_size_cap=2),
}
for label, options in cases.items():
    try:
        smsd.fingerprint(aspirin, options=options)
    except ValueError as error:
        print(f"{label}: {error}")
        assert str(error) == expected[label]
    else:
        raise AssertionError(f"{label} was accepted")

try:
    smsd.FingerprintOptions(element=1)
except TypeError as error:
    print(f"element=1: {error}")
    assert str(error) == "element must be True or False"
else:
    raise AssertionError("a non-boolean switch was accepted")
```

The program prints:

```text
assigned family gives classic FCFP: False
radius=-2: radius must be >= -1 (use -1 for unlimited)
ring_size_cap=2: ringSizeCap must be at least 3
element=1: element must be True or False
```

## Substructure keys and MACCS

The `"keys"` family sets one bit per substructure key, in the order given. A
key is a SMARTS string, a `(smarts, minimum)` pair or a
`FingerprintKey(smarts, minimum=1)`; its bit is set when at least `minimum`
distinct atom sets match. The size of a key fingerprint is the number of keys.

`FingerprintOptions.maccs()` gives the 166 MACCS keys. Key k is at position
k − 1, so add 1 to a position to read the key number. `maccs_keys()` lists the
definitions as `FingerprintKey` values with `smarts`, `minimum` and `kind`;
`kind` is `"smarts"` for most keys and `"isotope"`, `"aromatic_rings"` or
`"fragments"` for the few keys that count something else.

```python
import smsd

keys = smsd.FingerprintOptions("keys", keys=[
    "[#8]",                                   # key 0: an oxygen atom
    ("[#8]", 3),                              # key 1: three or more oxygens
    smsd.FingerprintKey("c1ccccc1"),          # key 2: a benzene ring
    smsd.FingerprintKey("[NX3]", minimum=2),  # key 3: two trivalent nitrogens or more
])
aspirin = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=keys)  # aspirin
paracetamol = smsd.fingerprint("CC(=O)Nc1ccc(O)cc1", options=keys)  # paracetamol
print("size:", keys.size)
print("aspirin:", aspirin)
print("paracetamol:", paracetamol)
assert keys.size == 4 and aspirin == [0, 1, 2] and paracetamol == [0, 2]

maccs = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O",
                         options=smsd.FingerprintOptions.maccs())
key_numbers = [position + 1 for position in maccs]
print("MACCS keys set for aspirin:", key_numbers)
assert len(smsd.FingerprintOptions.maccs_keys()) == 166
assert len(key_numbers) == 21
assert key_numbers[:3] == [89, 113, 123] and key_numbers[-1] == 165
```

The program prints:

```text
size: 4
aspirin: [0, 1, 2]
paracetamol: [0, 2]
MACCS keys set for aspirin: [89, 113, 123, 126, 127, 136, 139, 140, 143, 144, 146, 150, 152, 154, 157, 159, 160, 162, 163, 164, 165]
```

Each key search has fixed time and embedding limits. A key that cannot be
decided within them raises `RuntimeError` instead of returning a bit that could
depend on the machine.

## Many molecules at once

`fingerprints(mols, ...)` computes the fingerprint of every molecule in one
native call. It takes the same `kind`, `path_length`, `fp_size`, `radius` and
`options` arguments as `fingerprint`, and returns the results in input order.
The results equal one-at-a-time results whatever the thread count.
`num_threads=0` uses the OpenMP default, which `OMP_NUM_THREADS` controls. The
input may mix SMILES strings, `MolGraph` objects and RDKit molecules.

```python
import smsd

names = ["aspirin", "salicylic acid", "paracetamol", "ibuprofen", "naproxen",
         "caffeine", "phenol"]
library = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
           "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
           "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Oc1ccccc1"]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

batch = smsd.fingerprints(library, options=ecfp4, num_threads=2)
assert batch == [smsd.fingerprint(smiles, options=ecfp4) for smiles in library]
assert batch == smsd.fingerprints(library, options=ecfp4, num_threads=1)
mixed = [library[0], smsd.parse_smiles(library[1])]
assert smsd.fingerprints(mixed, options=ecfp4) == batch[:2]

for name, bits in zip(names, batch):
    print(f"{name:15} {len(bits):3d} set bits")
assert [len(bits) for bits in batch] == [36, 28, 24, 32, 46, 33, 14]

try:
    smsd.fingerprints(["CCO", "C1CC"], options=ecfp4)
except ValueError as error:
    print("ValueError:", error)
    assert str(error) == "Unclosed ring(s): 1"
else:
    raise AssertionError("an invalid SMILES string was accepted")
```

The program prints:

```text
aspirin          36 set bits
salicylic acid   28 set bits
paracetamol      24 set bits
ibuprofen        32 set bits
naproxen         46 set bits
caffeine         33 set bits
phenol           14 set bits
ValueError: Unclosed ring(s): 1
```

The first molecule that fails raises its own error, as `fingerprint` would.
A parse error does not name the input position; a ring-perception or key-search
refusal adds the zero-based input index to its message.

## Similarity coefficients

The binary coefficients compare set positions A and B:

- `tanimoto_coefficient`: |A ∩ B| / |A ∪ B|.
- `dice`: 2|A ∩ B| / (|A| + |B|).
- `cosine`: |A ∩ B| / √(|A| |B|).
- `overlap_coefficient`: |A ∩ B| / min(|A|, |B|), which is 1.0 when one
  fingerprint is contained in the other.
- `soergel`: 1 − Tanimoto, a distance rather than a similarity.

`dice`, `cosine` and `soergel` take the fingerprint width as `fp_size` (default
2048) and ignore positions at or above it; always pass the real width.
`tanimoto_coefficient` and `overlap_coefficient` do not read `fp_size`.

```python
import smsd

ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
a = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=ecfp4)  # aspirin
b = smsd.fingerprint("OC(=O)c1ccccc1O", options=ecfp4)        # salicylic acid
print("set bits:", len(a), len(b), "| in common:", len(set(a) & set(b)))
assert (len(a), len(b), len(set(a) & set(b))) == (36, 28, 10)

scores = {
    "tanimoto_coefficient": smsd.tanimoto_coefficient(a, b),
    "dice": smsd.dice(a, b, fp_size=2048),
    "cosine": smsd.cosine(a, b, fp_size=2048),
    "overlap_coefficient": smsd.overlap_coefficient(a, b),
    "soergel": smsd.soergel(a, b, fp_size=2048),
}
for name, score in scores.items():
    print(f"{name:21} {score:.3f}")
rounded = [round(score, 3) for score in scores.values()]
assert rounded == [0.185, 0.312, 0.315, 0.357, 0.815]
```

The program prints:

```text
set bits: 36 28 | in common: 10
tanimoto_coefficient  0.185
dice                  0.312
cosine                0.315
overlap_coefficient   0.357
soergel               0.815
```

The count coefficients compare count fingerprints a and b position by
position:

- `count_tanimoto_coefficient`: Σ min(a, b) / Σ max(a, b).
- `count_dice`: 2 Σ min(a, b) / (Σ a + Σ b).
- `count_cosine`: Σ ab / (√Σ a² √Σ b²).
- `count_overlap_coefficient`: Σ min(a, b) / min(Σ a, Σ b).

Pass `(position, count)` pairs, as count fingerprints are returned.
`count_dice` and `count_cosine` also take dense lists of equal length or
`{position: count}` dictionaries. `count_tanimoto_coefficient` and
`count_overlap_coefficient` read a dictionary by its keys only, as positions
with a count of 1. The binary `tanimoto_coefficient` accepts count pairs but
compares their positions only.

```python
import smsd

counted = smsd.FingerprintOptions("ecfp", counts=True)
a = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=counted)  # aspirin
b = smsd.fingerprint("OC(=O)c1ccccc1O", options=counted)        # salicylic acid
scores = {
    "count_tanimoto_coefficient": smsd.count_tanimoto_coefficient(a, b),
    "count_dice": smsd.count_dice(a, b),
    "count_cosine": smsd.count_cosine(a, b),
    "count_overlap_coefficient": smsd.count_overlap_coefficient(a, b),
    "tanimoto_coefficient": smsd.tanimoto_coefficient(a, b),  # positions only
}
for name, score in scores.items():
    print(f"{name:27} {score:.3f}")
rounded = [round(score, 3) for score in scores.values()]
assert rounded == [0.211, 0.348, 0.462, 0.4, 0.185]

pairs = smsd.count_tanimoto_coefficient([(1, 2)], [(1, 1)])
mapping = smsd.count_tanimoto_coefficient({1: 2}, {1: 1})  # keys only
print("pairs:", pairs, "| dictionaries:", mapping)
assert pairs == 0.5 and mapping == 1.0
```

The program prints:

```text
count_tanimoto_coefficient  0.211
count_dice                  0.348
count_cosine                0.462
count_overlap_coefficient   0.400
tanimoto_coefficient        0.185
pairs: 0.5 | dictionaries: 1.0
```

The width matters for `dice`, `cosine` and `soergel`. With 4,096-bit
fingerprints, the default `fp_size=2048` drops every position from 2,048
upwards and gives a different score without an error:

```python
import smsd

wide = smsd.FingerprintOptions.ecfp(radius=2, size=4096)
a = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=wide)  # aspirin
b = smsd.fingerprint("OC(=O)c1ccccc1O", options=wide)        # salicylic acid
right = smsd.dice(a, b, fp_size=4096)
wrong = smsd.dice(a, b)  # reads positions below 2048 only
print(f"dice with fp_size=4096: {right:.3f}")
print(f"dice with the default:  {wrong:.3f}")
assert round(right, 3) == 0.312 and round(wrong, 3) == 0.333
assert smsd.tanimoto_coefficient(a, b) == smsd.tanimoto_coefficient(
    a, b, fp_size=4096)
```

The program prints:

```text
dice with fp_size=4096: 0.312
dice with the default:  0.333
```

SMSD has no Tversky function. The asymmetric Tversky index is a few lines over
the set positions; with both weights 1 it equals Tanimoto, and with both
weights 0.5 it equals Dice.

```python
import smsd


def tversky(query, target, alpha=0.9, beta=0.1):
    q, t = set(query), set(target)
    common = len(q & t)
    denominator = common + alpha * len(q - t) + beta * len(t - q)
    return 1.0 if denominator == 0 else common / denominator


ecfp4 = smsd.FingerprintOptions.ecfp()
phenol = smsd.fingerprint("Oc1ccccc1", options=ecfp4)                # phenol
paracetamol = smsd.fingerprint("CC(=O)Nc1ccc(O)cc1", options=ecfp4)  # paracetamol
forward = tversky(phenol, paracetamol)
reverse = tversky(paracetamol, phenol)
print(f"phenol in paracetamol: {forward:.3f}")
print(f"paracetamol in phenol: {reverse:.3f}")
assert round(forward, 3) == 0.533 and round(reverse, 3) == 0.348
tanimoto = smsd.tanimoto_coefficient(phenol, paracetamol)
assert abs(tversky(phenol, paracetamol, 1, 1) - tanimoto) < 1e-12
dice = smsd.dice(phenol, paracetamol, fp_size=2048)
assert abs(tversky(phenol, paracetamol, 0.5, 0.5) - dice) < 1e-12
```

The program prints:

```text
phenol in paracetamol: 0.533
paracetamol in phenol: 0.348
```

`smsd.similarity(mol1, mol2)` is not a fingerprint score. It is the same as
`similarity_upper_bound`: an upper bound on the atom-count similarity of the
maximum common substructure, described under
[MCS similarity](MCS.md#mcs-similarity).

```python
import smsd

phenol, aspirin = "Oc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"  # phenol, aspirin
bound = smsd.similarity(phenol, aspirin)
ecfp4 = smsd.FingerprintOptions.ecfp()
tanimoto = smsd.tanimoto_coefficient(smsd.fingerprint(phenol, options=ecfp4),
                                     smsd.fingerprint(aspirin, options=ecfp4))
print(f"smsd.similarity: {bound:.3f}")
print(f"ECFP4 Tanimoto:  {tanimoto:.3f}")
assert bound == smsd.similarity_upper_bound(phenol, aspirin)
assert round(bound, 3) == 0.538 and round(tanimoto, 3) == 0.020
```

The program prints:

```text
smsd.similarity: 0.538
ECFP4 Tanimoto:  0.020
```

## Empty fingerprints

A molecule can have an empty fingerprint: methanol and ethanol, for example,
have no four-atom torsion. The coefficients treat empty input as follows:

| Coefficient | Both empty | One empty |
| --- | --- | --- |
| `tanimoto_coefficient` | 1.0 | 0.0 |
| `overlap_coefficient` | 1.0 | 0.0 |
| `dice` | 0.0 | 0.0 |
| `cosine` | 0.0 | 0.0 |
| `soergel` | 0.0 | 1.0 |
| `count_tanimoto_coefficient`, `count_overlap_coefficient`, `count_dice`, `count_cosine` | 0.0 | 0.0 |

`tanimoto_coefficient` also returns 0.0 when either argument is `None`.

```python
import smsd

rules = [
    ("tanimoto_coefficient", smsd.tanimoto_coefficient, [1], 1.0, 0.0),
    ("overlap_coefficient", smsd.overlap_coefficient, [1], 1.0, 0.0),
    ("dice", smsd.dice, [1], 0.0, 0.0),
    ("cosine", smsd.cosine, [1], 0.0, 0.0),
    ("soergel", smsd.soergel, [1], 0.0, 1.0),
    ("count_tanimoto_coefficient", smsd.count_tanimoto_coefficient,
     [(1, 1)], 0.0, 0.0),
    ("count_overlap_coefficient", smsd.count_overlap_coefficient,
     [(1, 1)], 0.0, 0.0),
    ("count_dice", smsd.count_dice, [(1, 1)], 0.0, 0.0),
    ("count_cosine", smsd.count_cosine, [(1, 1)], 0.0, 0.0),
]
for name, metric, one, both_empty, one_empty in rules:
    print(f"{name:27} {metric([], []):.1f} {metric([], one):.1f}")
    assert metric([], []) == both_empty and metric([], one) == one_empty

methanol = smsd.topological_torsion("CO")  # methanol
ethanol = smsd.topological_torsion("CCO")  # ethanol
score = smsd.tanimoto_coefficient(methanol, ethanol)
print("torsion bits:", methanol, ethanol, "| Tanimoto:", score)
assert methanol == [] and ethanol == [] and score == 1.0
```

The program prints:

```text
tanimoto_coefficient        1.0 0.0
overlap_coefficient         1.0 0.0
dice                        0.0 0.0
cosine                      0.0 0.0
soergel                     0.0 1.0
count_tanimoto_coefficient  0.0 0.0
count_overlap_coefficient   0.0 0.0
count_dice                  0.0 0.0
count_cosine                0.0 0.0
torsion bits: [] [] | Tanimoto: 1.0
```

Two molecules too small for a family therefore score 1.0 by Tanimoto. Filter
empty fingerprints out before ranking when that matters.

## Ranking a small set

`batch_similarity(query, targets, fp_size=2048, radius=2)` fingerprints the
query and each target with ECFP of the given radius and width, and returns
`(index, Tanimoto)` pairs from most to least similar. Targets with equal scores
keep their input order. For other families, or for libraries that are screened
repeatedly, see
[Similarity screen in memory](SCREENING.md#similarity-screen-in-memory) and the
[packed fingerprint index](PACKED_INDEX.md#build-an-index).

```python
import smsd

query = "CC(=O)Oc1ccccc1C(=O)O"  # aspirin
library = [
    "OC(=O)c1ccccc1O",               # salicylic acid
    "CC(=O)Nc1ccc(O)cc1",            # paracetamol
    "c1ccccc1",                      # benzene
    "CC(=O)Oc1ccccc1",               # phenyl acetate
    "CC(=O)Oc1ccccc1C(=O)OC",        # methyl acetylsalicylate
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # caffeine
    "OC(=O)c1ccccc1",                # benzoic acid
]
ranking = smsd.batch_similarity(query, library, fp_size=2048, radius=2)
for index, score in ranking[:3]:
    print(f"{score:.3f} {library[index]}")
top = [(index, round(score, 3)) for index, score in ranking[:3]]
assert top == [(4, 0.63), (3, 0.439), (6, 0.256)]

ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
q = smsd.fingerprint(query, options=ecfp4)
scores = [smsd.tanimoto_coefficient(q, t)
          for t in smsd.fingerprints(library, options=ecfp4)]
by_hand = sorted(enumerate(scores), key=lambda pair: -pair[1])
assert by_hand == ranking
```

The program prints:

```text
0.630 CC(=O)Oc1ccccc1C(=O)OC
0.439 CC(=O)Oc1ccccc1
0.256 OC(=O)c1ccccc1
```

## Encodings

`to_hex(fp, fp_size)` writes a binary fingerprint as hexadecimal text and
`from_hex(text)` reads it back. The text holds the fingerprint in 64-bit words
in ascending order, each written as 16 lower-case hexadecimal digits with the
most significant digit first, so bit 0 is the last digit of the first word. A
width that is not a multiple of 64 is padded to whole words. `from_hex` also
accepts upper case. This layout is not the byte layout of a packed fingerprint
file.

`to_binary_string(fp, fp_size)` writes one character per bit, character i
being bit i. `counts_to_array(counts, fp_size)` turns `(position, count)` pairs
or a `{position: count}` dictionary into a dense list, ignoring positions
outside the width.

```python
import smsd

phenol = smsd.fingerprint("Oc1ccccc1", kind="ecfp")  # phenol
text = smsd.to_hex(phenol, fp_size=2048)
print("hex characters:", len(text))
print("round trip:", smsd.from_hex(text) == phenol)
assert len(text) == 512 and smsd.from_hex(text) == phenol

for positions, width in (([0], 64), ([63], 64), ([0], 8)):
    print(f"to_hex({positions}, {width}):", smsd.to_hex(positions, fp_size=width))
assert smsd.to_hex([0], fp_size=64) == "0000000000000001"
assert smsd.to_hex([63], fp_size=64) == "8000000000000000"
assert smsd.to_hex([0], fp_size=8) == "0000000000000001"
assert smsd.from_hex("8000000000000000".upper()) == [63]

bits = smsd.to_binary_string([0, 3], fp_size=8)
print("to_binary_string([0, 3], 8):", bits)
assert bits == "10010000"

dense = smsd.counts_to_array([(1, 2), (9, 5)], 4)
print("counts_to_array([(1, 2), (9, 5)], 4):", dense)
assert dense == [0, 2, 0, 0]

try:
    smsd.from_hex("0x01")
except ValueError as error:
    print("ValueError:", error)
    assert str(error) == "hex_str must contain only hexadecimal digits"
else:
    raise AssertionError("non-hexadecimal text was accepted")
```

The program prints:

```text
hex characters: 512
round trip: True
to_hex([0], 64): 0000000000000001
to_hex([63], 64): 8000000000000000
to_hex([0], 8): 0000000000000001
to_binary_string([0, 3], 8): 10010000
counts_to_array([(1, 2), (9, 5)], 4): [0, 2, 0, 0]
ValueError: hex_str must contain only hexadecimal digits
```

## NumPy arrays

A fingerprint becomes a NumPy row by setting the listed positions to 1.
`smsd.adapters.fingerprints_numpy(molecules, ...)` builds the whole matrix
natively: one row per molecule, `uint8` for bits and `int64` for counts. Like
`fingerprints`, it defaults to `kind="mcs"`, so pass the kind or the options.
A Tanimoto matrix then follows from one matrix product. These functions need
NumPy (`python -m pip install "smsd[numpy]==8.1.0"`).

```python
import numpy as np
import smsd
from smsd.adapters import fingerprints_numpy

names = ["aspirin", "salicylic acid", "paracetamol", "phenol"]
smiles = ["CC(=O)Oc1ccccc1C(=O)O", "OC(=O)c1ccccc1O", "CC(=O)Nc1ccc(O)cc1",
          "Oc1ccccc1"]
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)

fps = smsd.fingerprints(smiles, options=ecfp4)
by_hand = np.zeros((len(fps), 2048), dtype=np.uint8)
for row, positions in enumerate(fps):
    by_hand[row, positions] = 1

matrix = fingerprints_numpy(smiles, options=ecfp4)
print("shape:", matrix.shape, "| dtype:", matrix.dtype)
assert matrix.shape == (4, 2048) and matrix.dtype == np.uint8
assert (matrix == by_hand).all()

a = matrix.astype(np.int32)
common = a @ a.T
sizes = a.sum(axis=1)
tanimoto = common / (sizes[:, None] + sizes[None, :] - common)
expected = [[smsd.tanimoto_coefficient(p, q) for q in fps] for p in fps]
assert np.allclose(tanimoto, expected)
for name, row in zip(names, tanimoto.round(3).tolist()):
    print(f"{name:15}", row)
assert tanimoto.round(3).tolist()[0] == [1.0, 0.185, 0.034, 0.02]
```

The program prints:

```text
shape: (4, 2048) | dtype: uint8
aspirin         [1.0, 0.185, 0.034, 0.02]
salicylic acid  [0.185, 1.0, 0.13, 0.273]
paracetamol     [0.034, 0.13, 1.0, 0.267]
phenol          [0.02, 0.273, 0.267, 1.0]
```

SciPy sparse matrices, pandas and Polars columns, Arrow buffers, PyTorch
tensors and scikit-learn transformers are covered in the
[Python feature adapters](../PYTHON_ADAPTERS.md) page.

## Fingerprint quality

`analyze_fp_quality(fp, fp_size=2048)` reports how full a binary fingerprint
is: `set_bits`, `fp_size`, `density`, `expected_collision_rate` (the estimated
share of hashed features that fell on a bit already set) and `is_saturated`
(density above 0.75). A saturated fingerprint has lost most of its ability to
separate molecules; use a larger width or a shorter path length.

```python
import smsd

diazepam = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"  # diazepam
for size in (2048, 256, 64):
    bits = smsd.fingerprint(diazepam, kind="mcs", fp_size=size)
    quality = smsd.analyze_fp_quality(bits, fp_size=size)
    row = (quality["set_bits"], round(quality["density"], 3),
           round(quality["expected_collision_rate"], 3), quality["is_saturated"])
    print(f"fp_size {size}: set {row[0]}, density {row[1]}, "
          f"collision rate {row[2]}, saturated {row[3]}")
    assert row == {2048: (297, 0.145, 0.074, False),
                   256: (182, 0.711, 0.426, False),
                   64: (64, 1.0, 1.0, True)}[size]
```

The program prints:

```text
fp_size 2048: set 297, density 0.145, collision rate 0.074, saturated False
fp_size 256: set 182, density 0.711, collision rate 0.426, saturated False
fp_size 64: set 64, density 1.0, collision rate 1.0, saturated True
```

`fingerprint_subset(fp1, fp2)` is true when every position of `fp1` is set in
`fp2`. Both lists must be sorted; unsorted input gives a wrong answer without an
error. It compares two sets of positions and is not a substructure test: with
the similarity fingerprints, a target that contains the query can fail it, and a
target that passes it need not contain the query.

```python
import smsd

ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
query = smsd.fingerprint("Oc1ccccc1", options=ecfp4)          # phenol
search = smsd.Substructure.find("Oc1ccccc1")
rows = []
for name, smiles in (("paracetamol", "CC(=O)Nc1ccc(O)cc1"),   # paracetamol
                     ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")):  # caffeine
    subset = smsd.fingerprint_subset(query, smsd.fingerprint(smiles, options=ecfp4))
    found = search.match(smiles).found
    print(f"{name:12} ECFP4 subset {subset}, contains phenol {found}")
    rows.append((subset, found))
assert rows == [(False, True), (False, False)]
assert smsd.fingerprint_subset([1, 5], [1, 2, 5]) and not smsd.fingerprint_subset([1, 6], [1, 2, 5])
```

The program prints:

```text
paracetamol  ECFP4 subset False, contains phenol True
caffeine     ECFP4 subset False, contains phenol False
```

For substructure screening of a library, use `Substructure.match_corpus`
(see [Search a prepared corpus](SUBSTRUCTURE.md#search-a-prepared-corpus)),
which rules out targets that cannot contain the query before matching, with a
proof.

## Options

Fingerprint functions, with their defaults. "Graph" means an SMSD `MolGraph`
only; "molecule" means a SMILES string, a `MolGraph` or an RDKit molecule.

| Function | Takes | Returns |
| --- | --- | --- |
| `fingerprint(mol, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None)` | molecule | Sorted positions, or `(position, count)` pairs when `options.counts` is true |
| `fingerprints(mols, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None, num_threads=0)` | iterable of molecules | One fingerprint per molecule, in input order |
| `fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp", *, path_length=7, options=None)` | non-empty SMILES string | As `fingerprint` |
| `circular_fingerprint(mol, radius=2, fp_size=2048, mode="ecfp")` | graph | Sorted positions; `mode` is `"ecfp"` or `"fcfp"` |
| `circular_fingerprint_counts(mol, radius=2, fp_size=2048, mode="ecfp")` | graph | `(position, count)` pairs |
| `path_fingerprint(mol, path_length=7, fp_size=2048)` | graph | Sorted positions |
| `mcs_fingerprint(mol, path_length=7, fp_size=2048, tautomer_aware=False)` | graph | Sorted positions; `tautomer_aware` hashes atoms of a tautomeric group without their element and tautomer class |
| `topological_torsion(mol, fp_size=2048)` | molecule | Sorted positions |
| `topological_torsion_counts(mol, fp_size=2048)` | graph | `(position, count)` pairs |
| `classify_pharmacophore(mol, atom_index)` | graph | Pharmacophore class bit mask of one atom |
| `tanimoto_coefficient(fp1, fp2, fp_size=2048)`, `overlap_coefficient(fp1, fp2, fp_size=2048)` | positions or pairs | `float`; `fp_size` is not read |
| `dice(fp1, fp2, fp_size=2048)`, `cosine(fp1, fp2, fp_size=2048)`, `soergel(fp1, fp2, fp_size=2048)` | positions | `float`; positions at or above `fp_size` are ignored |
| `count_tanimoto_coefficient(fp1, fp2)`, `count_overlap_coefficient(fp1, fp2)` | pairs or positions | `float` |
| `count_dice(fp1, fp2)`, `count_cosine(fp1, fp2)` | pairs, dictionaries or dense lists | `float`; dense lists of unequal length give 0.0 |
| `batch_similarity(query, targets, fp_size=2048, radius=2)` | molecules | `(index, Tanimoto)` pairs, most similar first |
| `to_hex(fp, fp_size=2048)`, `from_hex(hex_str)` | positions; text | Text; sorted positions |
| `to_binary_string(fp, fp_size=2048)` | positions | Text of `fp_size` characters |
| `counts_to_array(counts_dict, fp_size)` | pairs or dictionary | Dense list of length `fp_size` |
| `fingerprint_subset(fp1, fp2)` | sorted positions | `bool` |
| `analyze_fp_quality(fp, fp_size=2048)` | positions | `dict` |
| `smsd.adapters.fingerprints_numpy(molecules, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None, num_threads=0, dtype=None)` | molecules | NumPy matrix |

When `options` is given, it defines the fingerprint completely and `kind`,
`path_length`, `fp_size` and `radius` are not read. `fp_size` and
`path_length` must be at least 1, and `num_threads` must not be negative. With
`kind`, a radius of -1 means the whole component and a more negative radius is
treated as -1, while `FingerprintOptions` raises `ValueError` for it.

`FingerprintOptions(family="ecfp", **keywords)` takes the family `"ecfp"` (or
`"circular"`), `"fcfp"`, `"path"`, `"torsion"` or `"keys"`, as a string in any
case or as a `FingerprintFamily` value. The keywords, with their defaults for
the classic ECFP, FCFP, path and torsion configurations (a dash marks a value
the family does not read):

| Keyword | Type | ECFP | FCFP | Path | Torsion | Meaning |
| --- | --- | --- | --- | --- | --- | --- |
| `radius` | `int` | 2 | 2 | – | – | Circular radius; 2 gives ECFP4 or FCFP4, -1 the whole component |
| `path_length` | `int` | – | – | 7 | – | Longest path, in bonds |
| `size` | `int` | 2048 | 2048 | 2048 | 2048 | Bits or bins; for keys, the number of keys |
| `counts` | `bool` | False | False | False | False | Return `(position, count)` pairs |
| `element` | `bool` | True | False | True | True | Atomic number |
| `degree` | `bool` | True | True | False | True | Number of explicit neighbours |
| `valence` | `bool` | True | False | False | True | Bond order sum; pi electrons for torsion |
| `formal_charge`, `isotope`, `hydrogen_count`, `tautomer_class`, `aromaticity` | `bool` | True | False | False | False | The named atom property |
| `bond_order` | `bool` | True | True | False | False | Bond order of neighbour or path bonds |
| `aromatic_bond`, `ring_bond` | `bool` | False | False | False | False | Aromatic and ring flags of bonds |
| `ring_membership` | `bool` | True | True | False | True | Whether the atom is in a ring |
| `ring_size` | `bool` | False | False | False | False | Smallest ring through the atom |
| `ring_size_cap` | `int` | 8 | 8 | 8 | 8 | Rings of this size or more share a class; at least 3 |
| `cip_atoms` | `bool` | False | False | False | False | R/S (and r/s) labels of stereocentres |
| `cip_bonds` | `bool` | False | False | False | False | E/Z labels of double bonds |
| `fragments` | `str` or `FragmentMode` | `"whole"` | `"whole"` | `"whole"` | `"whole"` | `"whole"`, `"largest"` or `"per_fragment"` |
| `scaffold` | `str` or `ScaffoldMode` | `"none"` | `"none"` | `"none"` | `"none"` | `"none"` or `"murcko"` |
| `keys` | sequence | `[]` | `[]` | `[]` | `[]` | Keys of the `"keys"` family: SMARTS strings, `(smarts, minimum)` pairs or `FingerprintKey` values |

Every keyword is also a readable and writable attribute, and `copy()` returns
an independent copy. The class methods are `ecfp(radius=2, size=2048)`,
`fcfp(radius=2, size=2048)`, `path(path_length=7, size=2048)`,
`torsion(size=2048)`, `maccs()` and `maccs_keys()`.
`FingerprintKey(smarts, minimum=1)` has the read-only attributes `smarts`,
`minimum` and `kind`.

## Notes

- `FingerprintOptions` calls the width `size`; the functions call it
  `fp_size`.
- The FCFP pharmacophore class is always hashed and is not listed among the
  terms in the `repr`.
- Stereo terms use whole-molecule CIP labels. Radius 0 does not separate
  stereoisomers, and double bonds in rings of fewer than eight atoms carry no
  E/Z term.
- A whole-molecule fingerprint includes counter-ions and solvent; use
  `fragments="largest"` to fingerprint the parent.
- `mcs_fingerprint` is not configurable through `FingerprintOptions`.
- Ring-size terms, double-bond stereo terms, aromatic-ring keys and some
  SMARTS keys need a complete ring basis. When ring perception cannot complete
  it, they raise `smsd.RingBasisUnavailable` naming the cause.
- Fingerprints are always computed on CPU threads. A GPU applies only to some
  comparison steps; see [What runs where](CPU_GPU.md#what-runs-where).
- The same `FingerprintOptions` gives the same bits in C++ and Java. The C++
  `batchFingerprint` is a different all-atom path fingerprint and does not
  match Python `path_fingerprint`.
- `overlapCoefficient` and `count_overlapCoefficient` are deprecated aliases
  of `overlap_coefficient` and `count_overlap_coefficient`.

## See also

- [Screening a library](SCREENING.md#choose-a-route), for similarity and
  substructure screens over many targets.
- [Packed fingerprint index](PACKED_INDEX.md#build-an-index), for stored
  fingerprints of any specification, neighbours and clustering.
- The Python guide sections on
  [configurable fingerprints](../PYTHON.md#configurable-fingerprints) and
  [coefficients](../PYTHON.md#similarity-coefficients-and-encodings).
- [Python feature adapters](../PYTHON_ADAPTERS.md), for NumPy, SciPy, pandas,
  Polars, Arrow, PyTorch and scikit-learn.
- [Fingerprints and similarity](../EXAMPLES.md#fingerprints-and-similarity)
  in the Python examples.
- The [C++ guide](../CPP.md#configurable-fingerprints) and the
  [Java guide](../JAVA.md#configurable-fingerprints).
- The fingerprints and similarity
  [notebook](../../python/examples/fingerprints_similarity.ipynb), which runs
  the same steps on a small set of drugs.

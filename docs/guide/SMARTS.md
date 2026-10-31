# SMARTS queries

A SMARTS pattern describes the atoms and bonds a match must have, rather than
a whole molecule. SMSD patterns follow standard SMARTS semantics. This chapter
lists the primitives and operators SMSD accepts, shows each with a target it
matches and one it does not, and covers recursion, components,
stereochemistry, the chemistry options that apply, limits, error messages and
common traps. Searching with a pattern works as described in
[Substructure search](SUBSTRUCTURE.md).

## Compile and run a pattern

`Substructure.smarts(pattern)` compiles a pattern once. The compiled query
offers the same calls as a molecule query: `match`, `match_all` and
`match_many`. An empty or non-string pattern, or one that does not parse,
raises `ValueError` with the parser's message.

```python
import smsd

hydroxyl = smsd.Substructure.smarts("[OX2H]")   # hydroxyl oxygen
targets = ["CCO",       # ethanol
           "COC",       # dimethyl ether
           "CC(=O)O"]   # acetic acid

print("ethanol:", hydroxyl.match("CCO").found)
assert hydroxyl.match("CCO").found
found = [result.found for result in hydroxyl.match_many(targets, num_threads=1)]
print("ethanol, dimethyl ether, acetic acid:", found)
assert found == [True, False, True]

every = hydroxyl.with_timeout(1000).match_all("OCCO", limit=10)   # ethylene glycol
print("ethylene glycol:", len(every.matches), every.status.name)
assert len(every.matches) == 2 and every.complete

errors = {"": "SMARTS must be a non-empty string",
          "CC-": "Trailing bond in SMARTS",
          "[C": "Expected ']' in bracket atom at pos 2",
          "C(": "Unmatched '(' in SMARTS"}
for pattern, message in errors.items():
    try:
        smsd.Substructure.smarts(pattern)
    except ValueError as error:
        print(f"{pattern!r:6} -> {error}")
        assert str(error) == message
    else:
        raise AssertionError(pattern)
```

The program prints:

```text
ethanol: True
ethanol, dimethyl ether, acetic acid: [True, False, True]
ethylene glycol: 2 COMPLETED
''     -> SMARTS must be a non-empty string
'CC-'  -> Trailing bond in SMARTS
'[C'   -> Expected ']' in bracket atom at pos 2
'C('   -> Unmatched '(' in SMARTS
```

## Atom primitives

Inside square brackets, an atom is described by one or more primitives. An
element symbol outside brackets, such as `C`, `c` or `Cl`, is the same as the
symbol in brackets.

| Primitive | Meaning | Example | Matches | Does not match |
| --- | --- | --- | --- | --- |
| `*` | Any atom | `C*O` | `CCO` | `CO` |
| `#n` | Atomic number, aromatic or not | `[#7]` | `c1ccncc1` | `c1ccccc1` |
| Upper-case symbol | Aliphatic atom of that element | `[C]` | `CC` | `c1ccccc1` |
| Lower-case symbol | Aromatic atom of that element | `[c]` | `c1ccccc1` | `C1CCCCC1` |
| `a` | Aromatic atom | `[a]` | `c1ccncc1` | `C1CCCCC1` |
| `A` | Aliphatic atom | `[N;A]` | `CN` | `c1ccncc1` |
| `D<n>` | Explicit connections, counting explicit hydrogen atoms | `[D3]` | `CC(C)C` | `CCCC` |
| `X<n>` | Total connections, implicit hydrogens included | `[CX4]` | `CC` | `C=C` |
| `d<n>` | Heavy-atom neighbours | `[C;d3]` | `CC(C)C` | `CCC` |
| `v<n>` | Valence, counting bond orders and hydrogens | `[Cv4]` | `CC` | `[CH3]` |
| `H<n>` | Attached hydrogens | `[CH3]` | `CC` | `C=C` |
| `R` | In a ring | `[R]` | `C1CC1` | `CC` |
| `R<n>` | Number of rings that contain the atom; `R0` is not in a ring | `[R2]` | `c1ccc2ccccc2c1` | `c1ccccc1` |
| `r<n>` | The smallest ring that contains the atom has `n` atoms | `[c;r5]` | `c1ccoc1` | `c1ccccc1` |
| `x<n>` | Ring bonds | `[Cx2]` | `C1CCCCC1` | `CCC` |
| `z<n>` | Neighbours that are neither carbon nor hydrogen | `[Cz2]` | `OC(N)C` | `OCC` |
| `Z<n>` | Aliphatic neighbours that are neither carbon nor hydrogen | `[Z1]` | `CCO` | `Cc1ccccn1` |
| `^<n>` | Hybridisation: 0 to 5 for s, sp, sp2, sp3, sp3d and sp3d2 | `[C^1]` | `CC#N` | `CC=C` |
| `+`, `-`, `+<n>`, `-<n>` | Formal charge | `[#7;+]` | `C[N+](C)(C)C` | `CN` |
| Mass number | Isotope | `[13C]` | `[13CH4]` | `C` |
| `[H]` | A hydrogen atom, not a hydrogen count | `[H]` | `[H]C` | `C` |
| `:n` | Atom class; a label that does not restrict the match | `[C:1]O` | `CO` | `CN` |

Recursive SMARTS, `$(...)`, and the chirality primitives `@`, `@@`, `@?` and
`@@?` are described below. The program below checks every row of the table.

```python
import smsd

rows = [
    ("C*O", "CCO", "CO"),
    ("[#7]", "c1ccncc1", "c1ccccc1"),
    ("[C]", "CC", "c1ccccc1"),
    ("[c]", "c1ccccc1", "C1CCCCC1"),
    ("[a]", "c1ccncc1", "C1CCCCC1"),
    ("[N;A]", "CN", "c1ccncc1"),
    ("[D3]", "CC(C)C", "CCCC"),
    ("[CX4]", "CC", "C=C"),
    ("[C;d3]", "CC(C)C", "CCC"),
    ("[Cv4]", "CC", "[CH3]"),
    ("[CH3]", "CC", "C=C"),
    ("[R]", "C1CC1", "CC"),
    ("[R2]", "c1ccc2ccccc2c1", "c1ccccc1"),
    ("[c;r5]", "c1ccoc1", "c1ccccc1"),
    ("[Cx2]", "C1CCCCC1", "CCC"),
    ("[Cz2]", "OC(N)C", "OCC"),
    ("[Z1]", "CCO", "Cc1ccccn1"),
    ("[C^1]", "CC#N", "CC=C"),
    ("[#7;+]", "C[N+](C)(C)C", "CN"),
    ("[13C]", "[13CH4]", "C"),
    ("[H]", "[H]C", "C"),
    ("[C:1]O", "CO", "CN"),
]
for pattern, hit, miss in rows:
    query = smsd.Substructure.smarts(pattern)
    outcome = (query.match(hit).found, query.match(miss).found)
    print(f"{pattern:8} {hit:15} {outcome[0]!s:5} {miss:10} {outcome[1]}")
    assert outcome == (True, False)
```

The program prints:

```text
C*O      CCO             True  CO         False
[#7]     c1ccncc1        True  c1ccccc1   False
[C]      CC              True  c1ccccc1   False
[c]      c1ccccc1        True  C1CCCCC1   False
[a]      c1ccncc1        True  C1CCCCC1   False
[N;A]    CN              True  c1ccncc1   False
[D3]     CC(C)C          True  CCCC       False
[CX4]    CC              True  C=C        False
[C;d3]   CC(C)C          True  CCC        False
[Cv4]    CC              True  [CH3]      False
[CH3]    CC              True  C=C        False
[R]      C1CC1           True  CC         False
[R2]     c1ccc2ccccc2c1  True  c1ccccc1   False
[c;r5]   c1ccoc1         True  c1ccccc1   False
[Cx2]    C1CCCCC1        True  CCC        False
[Cz2]    OC(N)C          True  OCC        False
[Z1]     CCO             True  Cc1ccccn1  False
[C^1]    CC#N            True  CC=C       False
[#7;+]   C[N+](C)(C)C    True  CN         False
[13C]    [13CH4]         True  C          False
[H]      [H]C            True  C          False
[C:1]O   CO              True  CN         False
```

## Bond primitives and operators

| Bond | Meaning | Example | Matches | Does not match |
| --- | --- | --- | --- | --- |
| `-` | Single | `C-C` | `CC` | `C=C` |
| `=` | Double | `C=O` | `CC=O` | `CCO` |
| `#` | Triple | `C#N` | `CC#N` | `CC=N` |
| `:` | Aromatic | `[#6]:[#7]` | `c1ccncc1` | `Nc1ccccc1` |
| `~` | Any bond | `[#6]~[#7]` | `CC#N` | `CCO` |
| `@` | Ring bond | `C@C` | `C1CCCCC1` | `CCCCCC` |
| Unwritten | Single or aromatic | `CC` | `CC` | `C=C` |
| `/`, `\` | Double-bond configuration | See [Stereochemistry](#stereochemistry) | | |
| `/?`, `\?` | That configuration or none | See [Stereochemistry](#stereochemistry) | | |

Operators combine primitives in atom and bond expressions. In order of
precedence, from the highest:

| Operator | Meaning | Example | Matches | Does not match |
| --- | --- | --- | --- | --- |
| `!` | Not | `[!#6;!#1]` | `CCO` | `CC` |
| `&` or nothing | And | `[N,O&H1]` | `CN(C)C` | `CC` |
| `,` | Or | `[C,N]=O` | `CN=O` | `CS=O` |
| `;` | And, with low precedence | `[N,O;H1]` | `CNC` | `CN(C)C` |

`[N,O&H1]` reads as nitrogen, or oxygen with one hydrogen, so it matches the
nitrogen of trimethylamine. `[N,O;H1]` reads as nitrogen or oxygen, with one
hydrogen, so it does not.

```python
import smsd

rows = [
    ("C-C", "CC", "C=C"),
    ("C=O", "CC=O", "CCO"),
    ("C#N", "CC#N", "CC=N"),
    ("[#6]:[#7]", "c1ccncc1", "Nc1ccccc1"),
    ("[#6]~[#7]", "CC#N", "CCO"),
    ("C@C", "C1CCCCC1", "CCCCCC"),
    ("CC", "CC", "C=C"),
    ("[!#6;!#1]", "CCO", "CC"),
    ("[N,O&H1]", "CN(C)C", "CC"),
    ("[C,N]=O", "CN=O", "CS=O"),
    ("[N,O;H1]", "CNC", "CN(C)C"),
]
for pattern, hit, miss in rows:
    query = smsd.Substructure.smarts(pattern)
    outcome = (query.match(hit).found, query.match(miss).found)
    print(f"{pattern:10} {hit:9} {outcome[0]!s:5} {miss:10} {outcome[1]}")
    assert outcome == (True, False)
```

The program prints:

```text
C-C        CC        True  C=C        False
C=O        CC=O      True  CCO        False
C#N        CC#N      True  CC=N       False
[#6]:[#7]  c1ccncc1  True  Nc1ccccc1  False
[#6]~[#7]  CC#N      True  CCO        False
C@C        C1CCCCC1  True  CCCCCC     False
CC         CC        True  C=C        False
[!#6;!#1]  CCO       True  CC         False
[N,O&H1]   CN(C)C    True  CC         False
[C,N]=O    CN=O      True  CS=O       False
[N,O;H1]   CNC       True  CN(C)C     False
```

## Recursive SMARTS

`$(...)` matches an atom that is the first atom of the enclosed pattern, so it
tests the atom's surroundings without adding them to the match. Alternatives
are joined with `,`. Recursive patterns may be nested up to 20 levels deep;
Python cannot change this limit, and a deeper pattern raises `ValueError`.

```python
import smsd

# The nitrogen of an amide: bonded to a carbon that bears a double-bonded oxygen
amide_n = smsd.Substructure.smarts("[$(C=O)]N")
print("acetamide:", amide_n.match("CC(=O)N").found,
      "ethylamine:", amide_n.match("CCN").found)
assert amide_n.match("CC(=O)N").found and not amide_n.match("CCN").found

# A carboxyl carbon or a nitrile carbon
either = smsd.Substructure.smarts("[$(C(=O)O),$(C#N)]")
names = {"CC(=O)O": "acetic acid", "CC#N": "acetonitrile", "CC=O": "acetaldehyde"}
found = {name: either.match(smiles).found for smiles, name in names.items()}
print(found)
assert found == {"acetic acid": True, "acetonitrile": True, "acetaldehyde": False}


def nested(depth):
    pattern = "C"
    for _ in range(depth):
        pattern = f"[$({pattern})]"
    return pattern


deepest = smsd.Substructure.smarts(nested(20)).match("C").found
print("depth 20:", deepest)
assert deepest
try:
    smsd.Substructure.smarts(nested(21))
except ValueError as error:
    print("depth 21:", error)
    assert "maximum depth of 20" in str(error)
```

The program prints:

```text
acetamide: True ethylamine: False
{'acetic acid': True, 'acetonitrile': True, 'acetaldehyde': False}
depth 20: True
depth 21: Recursive SMARTS nesting exceeds maximum depth of 20
```

## Components

A dot separates pattern components. Each component may match anywhere, with
no shared atoms. Parentheses around the whole pattern or around components
group them by target component:

- `(C.C)` requires both atoms in the same target component.
- `(C).(C)` requires them in different target components.

```python
import smsd

rows = [("C.C", ["CC", "C.C"], [True, True]),
        ("(C.C)", ["CC.CC", "C.C"], [True, False]),
        ("(C).(C)", ["CC", "C.C"], [False, True])]
for pattern, targets, expected in rows:
    query = smsd.Substructure.smarts(pattern)
    found = [query.match(target).found for target in targets]
    print(f"{pattern:8}", dict(zip(targets, found)))
    assert found == expected
```

The program prints:

```text
C.C      {'CC': True, 'C.C': True}
(C.C)    {'CC.CC': True, 'C.C': False}
(C).(C)  {'CC': False, 'C.C': True}
```

## Stereochemistry

Stereochemistry written in a pattern is always enforced; `use_chirality` and
`use_bond_stereo` are not needed. A pattern with a configuration does not match
a target that leaves it unspecified, unless the pattern uses `@?`, `@@?`, `/?`
or `\?`, which also accept no configuration. A double bond with a direction
written on one side only places no constraint.

A tetrahedral tag is read relative to the written order of the neighbours, as
in SMILES. To compare a pattern with a target, write both with the neighbours
in the same order.

```python
import smsd

trans = smsd.Substructure.smarts("F/C=C/F")
targets = {"trans": "F/C=C/F", "cis": "F/C=C\\F", "unspecified": "FC=CF"}
found = {name: trans.match(smiles).found for name, smiles in targets.items()}
print("F/C=C/F:", found)
assert found == {"trans": True, "cis": False, "unspecified": False}

either = smsd.Substructure.smarts("F/C=C/?F")
found = {name: either.match(smiles).found for name, smiles in targets.items()}
print("F/C=C/?F:", found)
assert found == {"trans": True, "cis": False, "unspecified": True}

one_side = smsd.Substructure.smarts("F/C=CF")
found = {name: one_side.match(smiles).found for name, smiles in targets.items()}
print("F/C=CF:", found)
assert found == {"trans": True, "cis": True, "unspecified": True}

centres = {"@": "[C@H](F)(Cl)Br", "@@": "[C@@H](F)(Cl)Br", "none": "C(F)(Cl)Br"}
for pattern in ("[C@H](F)(Cl)Br", "[C@?H](F)(Cl)Br"):
    query = smsd.Substructure.smarts(pattern)
    found = {name: query.match(smiles).found for name, smiles in centres.items()}
    print(f"{pattern}:", found)
    assert found["@"] and not found["@@"] and found["none"] == ("?" in pattern)
```

The program prints:

```text
F/C=C/F: {'trans': True, 'cis': False, 'unspecified': False}
F/C=C/?F: {'trans': True, 'cis': False, 'unspecified': True}
F/C=CF: {'trans': True, 'cis': True, 'unspecified': True}
[C@H](F)(Cl)Br: {'@': True, '@@': False, 'none': False}
[C@?H](F)(Cl)Br: {'@': True, '@@': False, 'none': True}
```

## Which options apply

A pattern states its own element, charge, hydrogen count, aromaticity, bond
order and stereochemistry. Of the chemistry options, a SMARTS query therefore
reads only `complete_rings_only`, together with `Substructure.Options.induced`.
The other options have no effect on a pattern: put such constraints into the
pattern itself. The time limit, node limit and cancellation token apply as
usual.

```python
import smsd


def found(pattern, target, induced=False, **fields):
    options = smsd.Substructure.Options()
    options.induced = induced
    for name, value in fields.items():
        setattr(options.chemistry, name, value)
    return smsd.Substructure.smarts(pattern).with_options(options).match(target).found


rows = [
    ("C-C in ethene", "C-C", "C=C",
     dict(match_bond_order=smsd.BondOrderMode.ANY), (False, False)),
    ("C[N] in methylammonium", "C[N]", "C[NH3+]",
     dict(match_formal_charge=True), (True, True)),
    ("CCC in cyclohexane", "CCC", "C1CCCCC1",
     dict(ring_matches_ring_only=True), (True, True)),
    ("CCC in cyclohexane", "CCC", "C1CCCCC1",
     dict(complete_rings_only=True), (True, False)),
    ("CCC in cyclopropane", "CCC", "C1CC1", dict(induced=True), (True, False)),
]
for label, pattern, target, fields, expected in rows:
    outcome = (found(pattern, target), found(pattern, target, **fields))
    setting = ", ".join(fields)
    print(f"{label:23} {setting:23} {outcome[0]!s:5} -> {outcome[1]}")
    assert outcome == expected
```

The program prints:

```text
C-C in ethene           match_bond_order        False -> False
C[N] in methylammonium  match_formal_charge     True  -> True
CCC in cyclohexane      ring_matches_ring_only  True  -> True
CCC in cyclohexane      complete_rings_only     True  -> False
CCC in cyclopropane     induced                 True  -> False
```

## Writing SMARTS from a molecule

`smsd.to_smarts(molecule, ...)` writes a molecule as SMARTS. Its switches
choose which properties the pattern fixes:

| Argument | Default | Adds |
| --- | --- | --- |
| `include_aromaticity` | `True` | `a` or `A` on each atom |
| `include_charge` | `True` | The formal charge of each charged atom |
| `include_ring_member` | `False` | `R` or `!R` on each atom |
| `include_h_count` | `False` | `H<n>` on each atom |
| `include_isotope` | `False` | The mass number of each labelled atom |

`Match.to_smarts(query, chemistry=None)` writes only the matched atoms and
bonds of a result. The pattern encodes the matching rules of `chemistry`, or of
the default options when it is omitted, and a rule that SMARTS cannot express,
such as tautomer matching, raises `ValueError`. An MCS can be turned into a
pattern in one step with `MCS.match_with_smarts`; see
[Maximum common substructure](MCS.md#export-the-common-graph) and
[SMARTS export](../SMARTS_EXPORT.md).

```python
import smsd

phenol = "c1ccccc1O"
written = [
    (smsd.to_smarts(phenol),
     "[#6;a]:1(:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]:1)-[#8;A]"),
    (smsd.to_smarts(phenol, include_ring_member=True),
     "[#6;a;R]:1(:[#6;a;R]:[#6;a;R]:[#6;a;R]:[#6;a;R]:[#6;a;R]:1)-[#8;A;!R]"),
    (smsd.to_smarts("C[NH3+]"), "[#6;A]-[#7;A;+]"),
    (smsd.to_smarts("C[NH3+]", include_charge=False), "[#6;A]-[#7;A]"),
    (smsd.to_smarts("[13CH4]"), "[#6;A]"),
    (smsd.to_smarts("[13CH4]", include_isotope=True), "[13#6;A]"),
]
for text, expected in written:
    print(text)
    assert text == expected

# Acetaldehyde found in propanoic acid, written as SMARTS
result = smsd.Substructure.find("CC=O").match("CCC(=O)O")
pattern = result.best.to_smarts("CC=O")
print(pattern)
assert pattern == "[#6]-,:[#6]=,:[#8]"
assert smsd.Substructure.smarts(pattern).match("CCC(=O)O").found
```

The program prints:

```text
[#6;a]:1(:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]:1)-[#8;A]
[#6;a;R]:1(:[#6;a;R]:[#6;a;R]:[#6;a;R]:[#6;a;R]:[#6;a;R]:1)-[#8;A;!R]
[#6;A]-[#7;A;+]
[#6;A]-[#7;A]
[#6;A]
[13#6;A]
[#6]-,:[#6]=,:[#8]
```

Under the default flexible aromaticity, aromatic bonds also match single or
double bonds, so the exported bonds read `-,:` and `=,:`. To
fix hydrogen counts, see
[Hydrogen counts](SUBSTRUCTURE.md#hydrogen-counts).

## Limits

| Limit | Value |
| --- | --- |
| Pattern length, in characters | 100,000 |
| Atoms, counting those inside recursive SMARTS | 10,000 |
| Expression nodes in the whole pattern | 100,000 |
| Expression nodes in one atom or bond expression | 10,000 |
| Nesting depth of recursive SMARTS | 20 |

A pattern beyond a limit raises `ValueError`, for example
`SMARTS exceeds the maximum length of 100000 characters`.

## Unsupported syntax

- The implicit-hydrogen primitive `h` is not supported. `[Ch2]` raises
  `ValueError` ("Unknown aromatic element in SMARTS: h"); write `H<n>` or
  `X<n>` instead.
- Named predicates such as `$isHalogen` are not accepted in Python. Write the
  definition as SMARTS, using `$(...)` where it describes an atom's
  surroundings. The C++ and Java libraries can expand named predicates; see
  [C++ SMARTS](../CPP.md#smarts) and
  [Java prepared searches](../JAVA.md#prepared-searches).

```python
import smsd

errors = {"[Ch2]": "Unknown aromatic element in SMARTS: h",
          "[$isHalogen]": "Expected '(' after '$' in recursive SMARTS"}
for pattern, message in errors.items():
    try:
        smsd.Substructure.smarts(pattern)
    except ValueError as error:
        print(f"{pattern:13} -> {error}")
        assert str(error) == message
    else:
        raise AssertionError(pattern)

halogen = smsd.Substructure.smarts("[F,Cl,Br,I]")
found = [halogen.match(smiles).found for smiles in ("CCl", "CBr", "CO")]
print("halogen in chloromethane, bromomethane, methanol:", found)
assert found == [True, True, False]
```

The program prints:

```text
[Ch2]         -> Unknown aromatic element in SMARTS: h
[$isHalogen]  -> Expected '(' after '$' in recursive SMARTS
halogen in chloromethane, bromomethane, methanol: [True, True, False]
```

## Notes

- SMARTS `C` is an aliphatic carbon. SMSD perceives aromaticity when it reads
  a molecule, so Kekulé input can become aromatic: caffeine written as
  `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` has aromatic carbonyl carbons, and
  `C(=O)[NX3]` misses them. Write `[#6]` when aromaticity should not matter.
- Inside brackets, a two-letter element symbol takes precedence: `[Cd3]` is
  read as cadmium, not as carbon with `d3`. Write `[C;d3]` for a carbon with
  three heavy-atom neighbours.
- A molecule written as SMARTS is stricter than the same molecule used as a
  query. The pattern's lower-case atoms must be aromatic, while a molecule
  query with the default flexible aromaticity also matches non-aromatic rings.
- `D<n>` counts explicit hydrogen atoms, such as those written `[H]`, while
  `d<n>` counts heavy atoms only.

The program below shows the first and third notes on shared example data.

```python
import smsd

caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
print("read as:", smsd.to_smiles(caffeine))
assert smsd.to_smiles(caffeine) == "Cn1c2c(c(n(C)c1=O)=O)n(C)cn2"
aliphatic = smsd.Substructure.smarts("C(=O)[NX3]").match(caffeine).found
any_carbon = smsd.Substructure.smarts("[#6](=O)[#7X3]").match(caffeine).found
print("C(=O)[NX3]:", aliphatic, "[#6](=O)[#7X3]:", any_carbon)
assert not aliphatic and any_carbon

benzodiazepines = {
    "diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    "nordazepam": "O=C1CN=C(c2ccccc2)c2cc(Cl)ccc2N1",
    "temazepam": "CN1C(=O)C(O)N=C(c2ccccc2)c2cc(Cl)ccc21",
    "oxazepam": "OC1N=C(c2ccccc2)c2cc(Cl)ccc2NC1=O",
    "lorazepam": "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O",
    "clonazepam": "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1",
    "nitrazepam": "O=C1CN=C(c2ccccc2)c2cc([N+](=O)[O-])ccc2N1",
    "flunitrazepam": "CN1C(=O)CN=C(c2ccccc2F)c2cc([N+](=O)[O-])ccc21",
    "bromazepam": "O=C1CN=C(c2ccccn2)c2cc(Br)ccc2N1",
    "tetrazepam": "CN1C(=O)CN=C(C2=CCCCC2)c2cc(Cl)ccc21",
    "alprazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1)=NC2",
    "triazolam": "Cc1nnc2n1-c1ccc(Cl)cc1C(c1ccccc1Cl)=NC2",
    "midazolam": "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2",
    "clobazam": "CN1C(=O)CC(=O)N(c2ccccc2)c2cc(Cl)ccc21",
}
core = "O=C1CN=C(c2ccccc2)c2ccccc2N1"   # 5-phenyl-1,4-benzodiazepin-2-one
names, targets = list(benzodiazepines), list(benzodiazepines.values())
hits = {}
for label, query in (("molecule query", smsd.Substructure.find(core)),
                     ("SMARTS", smsd.Substructure.smarts(core))):
    results = query.match_many(targets, num_threads=1)
    assert all(result.status == smsd.SearchStatus.COMPLETED for result in results)
    hits[label] = [name for name, result in zip(names, results) if result.found]
    print(f"{label}: {len(hits[label])} hits")
missing = sorted(set(hits["molecule query"]) - set(hits["SMARTS"]))
print("found only by the molecule query:", missing)
assert len(hits["molecule query"]) == 9 and missing == ["tetrazepam"]
```

The program prints:

```text
read as: Cn1c2c(c(n(C)c1=O)=O)n(C)cn2
C(=O)[NX3]: False [#6](=O)[#7X3]: True
molecule query: 9 hits
SMARTS: 8 hits
found only by the molecule query: ['tetrazepam']
```

## See also

- [Substructure search](SUBSTRUCTURE.md)
- [Maximum common substructure](MCS.md)
- [SMARTS export](../SMARTS_EXPORT.md)
- [C++ SMARTS](../CPP.md#smarts)
- [Java prepared searches](../JAVA.md#prepared-searches)

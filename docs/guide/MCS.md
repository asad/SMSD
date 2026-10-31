# Maximum common substructure

A maximum common substructure (MCS) search finds the largest graph that two
molecules share under a set of chemistry rules. SMSD Pro offers two
objectives: the maximum common induced subgraph (MCIS), which maximises common
atoms, and the maximum common edge subgraph (MCES), which maximises common
bonds. Either can be connected or split into fragments. Every result states
whether it is proved optimal, so a stopped search never passes for a complete
one. This chapter covers the objectives, limits and partial results, chemistry
options, atom weights, reserved atoms, common graphs of a series, validation,
MCS similarity, export, drawing and tautomers.

## Prepare and run a search

`smsd.MCS.find(query)` prepares a query once from a SMILES string, a
`MolGraph` or an RDKit molecule. The prepared query is immutable and can be
shared between threads; `with_timeout`, `with_options`, `with_objective`,
`with_fragments` and `with_atom_weights` each return a new query and leave the
original unchanged. `match(target)` returns a `SearchResult` with at most one
`Match`.

```python
import smsd

query, target = "c1ccc(O)cc1", "c1ccc(N)cc1"  # phenol, aniline
prepared = smsd.MCS.find(query).with_timeout(1000)
result = prepared.match(target)
best = result.best

assert result.found and result.status == smsd.SearchStatus.COMPLETED
assert result.optimal and result.complete and result.score == 6.0
assert (best.atom_count, best.bond_count) == (6, 6)
mapping = sorted(best.atom_mapping.items())
assert mapping == [(0, 0), (1, 1), (2, 2), (3, 3), (5, 5), (6, 6)]
assert best.bond_mapping[0] == (0, 1, 0, 1) and len(best.bond_mapping) == 6
assert best.fragments == ((0, 1, 2, 3, 5, 6),)
assert best.to_smiles(query) == "c1ccccc1"

print("status:", result.status.name, "| optimal:", result.optimal,
      "| complete:", result.complete)
print("score:", result.score, "| atoms:", best.atom_count,
      "| bonds:", best.bond_count)
print("atom mapping:", mapping)
print("first bond:", best.bond_mapping[0], "of", len(best.bond_mapping))
print("fragments:", best.fragments)
print("common graph:", best.to_smiles(query))
```

The program prints:

```text
status: COMPLETED | optimal: True | complete: True
score: 6.0 | atoms: 6 | bonds: 6
atom mapping: [(0, 0), (1, 1), (2, 2), (3, 3), (5, 5), (6, 6)]
first bond: (0, 1, 0, 1) of 6
fragments: ((0, 1, 2, 3, 5, 6),)
common graph: c1ccccc1
```

The result fields are:

| Field | Meaning |
| --- | --- |
| `found`, `best`, `matches` | Whether a match exists; the first `Match` or `None`; all matches as a tuple |
| `status` | `COMPLETED`, `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED` |
| `optimal` | The optimum is proved |
| `complete` | The search, or the enumeration of `match_all`, finished |
| `score` | Atoms for MCIS, bonds for MCES, weight units with atom weights |
| `nodes_visited`, `elapsed_ms` | Search effort and elapsed milliseconds, for diagnostics only |

Each `Match` holds `atom_mapping` (a read-only map from query atom to target
atom), `bond_mapping` (one `(query_first, query_second, target_first,
target_second)` tuple per selected bond), `fragments` (the components of the
common graph, as query atom indices), `atom_count` and `bond_count`. Indices
follow the caller's numbering, including the atom order of an RDKit molecule.
The mapping of an RDKit input is not sorted by key, so sort it before printing.
[Reading a result](GETTING_STARTED.md#reading-a-result) explains the three
status fields in more depth.

## MCIS and MCES

MCIS, the default, ranks common graphs by atoms and then by bonds. The mapped
atoms keep every bond and every non-bond between them, so the common graph is
an induced subgraph of both molecules. MCES ranks by bonds and then by atoms.
It may leave out a bond between two mapped atoms, and it never keeps an atom
without a selected bond.

```python
import smsd

query, target = "C1CC1", "CC(C)C"  # cyclopropane, isobutane
induced = smsd.MCS.find(query).match(target)
edges = smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES).match(target)

rows = []
for label, result in (("MCIS", induced), ("MCES", edges)):
    assert result.optimal
    best = result.best
    rows.append((label, best.atom_count, best.bond_count, result.score,
                 best.to_smiles(query)))
assert rows == [("MCIS", 2, 1, 2.0, "CC"), ("MCES", 3, 2, 2.0, "C(C)C")]
for label, atoms, bonds, score, smiles in rows:
    print(f"{label}: atoms {atoms}, bonds {bonds}, score {score}, {smiles}")
```

The program prints:

```text
MCIS: atoms 2, bonds 1, score 2.0, CC
MCES: atoms 3, bonds 2, score 2.0, C(C)C
```

Both scores are 2.0, but they count different things: two atoms for MCIS and
two bonds for MCES. Compare scores only between searches with the same
objective.

## Connected or fragments

By default the common graph is connected. `with_fragments(min_atoms=1,
max_fragments=2147483647)` allows several components, each of at least
`min_atoms` atoms, up to `max_fragments` components. With `min_atoms=1`, MCIS
can keep isolated atoms; MCES never does.

```python
import smsd

query, target = "CCOCC", "CCNCC"  # diethyl ether, diethylamine
mces = smsd.MCS.Objective.MCES
connected = smsd.MCS.find(query).match(target)
split = smsd.MCS.find(query).with_fragments().match(target)
large = smsd.MCS.find(query).with_fragments(min_atoms=3).match(target)
edges = smsd.MCS.find(query).with_objective(mces).with_fragments().match(target)

assert connected.best.atom_count == 2
assert split.best.atom_count == 4 and split.best.fragments == ((0, 1), (3, 4))
assert split.best.to_smiles(query) == "CC.CC"
assert not large.found and large.status == smsd.SearchStatus.COMPLETED
assert large.optimal
assert edges.best.bond_count == 2 and edges.score == 2.0
print("connected MCIS:", connected.best.atom_count, "atoms")
print("fragments MCIS:", split.best.atom_count, "atoms",
      split.best.fragments, split.best.to_smiles(query))
print("fragments of 3 atoms or more:", large.found, large.status.name,
      "| optimal:", large.optimal)
print("fragments MCES:", edges.best.bond_count, "bonds, score", edges.score)

single = smsd.MCS.find("C.O").with_fragments().match("C.O")  # methane, water
none = smsd.MCS.find("C.O").with_objective(mces).with_fragments().match("C.O")
assert single.best.fragments == ((0,), (1,)) and not none.found
print("C.O, MCIS fragments:", single.best.fragments,
      "| MCES found:", none.found)

try:
    smsd.MCS.find(query).with_fragments(min_atoms=0)
except ValueError as error:
    message = str(error)
else:
    raise AssertionError("min_atoms=0 was accepted")
assert message == "Fragment limits must be positive"
print("min_atoms=0:", message)
```

The program prints:

```text
connected MCIS: 2 atoms
fragments MCIS: 4 atoms ((0, 1), (3, 4)) CC.CC
fragments of 3 atoms or more: False COMPLETED | optimal: True
fragments MCES: 2 bonds, score 2.0
C.O, MCIS fragments: ((0,), (1,)) | MCES found: False
min_atoms=0: Fragment limits must be positive
```

A result with no match, `COMPLETED` status and `optimal=True` proves that no
qualifying common graph exists: here, no pair of fragments of three or more
atoms each.

## Proof, limits and partial results

Three limits stop a search early:

- `timeout_ms`: by default the time limit adapts to the size of the pair; `0`
  stops at once, and a positive value applies to each query and target pair;
- `node_limit`: `0`, the default, sets no limit; reaching one reports
  `LIMIT_REACHED`;
- `cancellation`: a `CancellationToken` whose `cancel()` can be called from any
  thread; the search reports `CANCELLED`. A token cannot be reset, so use a new
  one for each search.

A stopped search can still return a valid common graph, but it is not proved
to be the largest. Describe a result as the maximum only when `optimal` is
true.

```python
import smsd

query = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # ibuprofen
target = "COc1ccc2cc(ccc2c1)C(C)C(=O)O"  # naproxen

full = smsd.MCS.find(query).with_timeout(2000).match(target)
assert full.status == smsd.SearchStatus.COMPLETED and full.optimal
assert full.best.atom_count == 14
print("no node limit:", full.status.name, "| optimal:", full.optimal,
      "| atoms:", full.best.atom_count)

rows = []
for limit in (2, 5, 20, 100):
    options = smsd.MCS.Options()
    options.timeout_ms = 2000
    options.node_limit = limit
    result = smsd.MCS.find(query).with_options(options).match(target)
    valid = smsd.MCS.validate_mapping(query, target, result.best, options) == []
    rows.append((limit, result.status.name, result.optimal,
                 result.best.atom_count, valid))
assert rows == [
    (2, "LIMIT_REACHED", False, 2, True),
    (5, "LIMIT_REACHED", False, 5, True),
    (20, "LIMIT_REACHED", False, 13, True),
    (100, "COMPLETED", True, 14, True),
]
for limit, status, optimal, atoms, valid in rows:
    print(f"node_limit {limit}: {status} | optimal: {optimal} "
          f"| atoms: {atoms} | valid: {valid}")

stopped = smsd.MCS.find(query).with_timeout(0).match(target)
token = smsd.CancellationToken()
token.cancel()
options = smsd.MCS.Options()
options.cancellation = token
cancelled = smsd.MCS.find(query).with_options(options).match(target)
for label, result in (("timeout 0", stopped), ("cancelled", cancelled)):
    assert not (result.found or result.optimal or result.complete)
    print(f"{label}: {result.status.name} | found: {result.found}")
assert stopped.status == smsd.SearchStatus.TIMED_OUT
assert cancelled.status == smsd.SearchStatus.CANCELLED
```

The program prints:

```text
no node limit: COMPLETED | optimal: True | atoms: 14
node_limit 2: LIMIT_REACHED | optimal: False | atoms: 2 | valid: True
node_limit 5: LIMIT_REACHED | optimal: False | atoms: 5 | valid: True
node_limit 20: LIMIT_REACHED | optimal: False | atoms: 13 | valid: True
node_limit 100: COMPLETED | optimal: True | atoms: 14 | valid: True
timeout 0: TIMED_OUT | found: False
cancelled: CANCELLED | found: False
```

An empty result proves absence only when it is `COMPLETED` and optimal. An
empty result from a stopped search proves nothing.

## All optimal mappings

`match_all(target, limit=1000)` returns every optimal mapping, including those
related by symmetry, up to `limit`. When the limit cuts the enumeration short,
the status is `LIMIT_REACHED` with `complete=False`, while `optimal` can stay
true because the optimum itself is proved. `MCS.canonicalize_mapping` picks one
representative under the symmetries of both molecules, which removes the
duplicates.

```python
import smsd

query, target = "C1CC1", "CC(C)C"  # cyclopropane, isobutane
prepared = smsd.MCS.find(query).with_timeout(1000)
everything = prepared.match_all(target, limit=100)
capped = prepared.match_all(target, limit=2)

assert len(everything.matches) == 18
assert everything.status == smsd.SearchStatus.COMPLETED
assert everything.optimal and everything.complete
assert len(capped.matches) == 2
assert capped.status == smsd.SearchStatus.LIMIT_REACHED
assert capped.optimal and not capped.complete
print("limit 100:", len(everything.matches), "mappings,",
      everything.status.name, "| complete:", everything.complete)
print("limit 2:", len(capped.matches), "mappings,", capped.status.name,
      "| optimal:", capped.optimal, "| complete:", capped.complete)

canonical = {
    tuple(sorted(smsd.MCS.canonicalize_mapping(
        query, target, dict(match.atom_mapping)).items()))
    for match in everything.matches
}
assert canonical == {((0, 0), (1, 1))}
print("distinct up to symmetry:", len(canonical), sorted(canonical))
```

The program prints:

```text
limit 100: 18 mappings, COMPLETED | complete: True
limit 2: 2 mappings, LIMIT_REACHED | optimal: True | complete: False
distinct up to symmetry: 1 [((0, 0), (1, 1))]
```

## Chemistry options in an MCS

The rules that decide which atoms and bonds may correspond live in
`options.chemistry`, a `ChemOptions` object edited in place. Option objects
take no keyword arguments; set each field as an attribute. The full table of
fields and named profiles is in
[Chemistry options](SUBSTRUCTURE.md#chemistry-options); this section shows how
the main fields change an MCS.

```python
import smsd


def common_atoms(query, target, **fields):
    options = smsd.MCS.Options()
    options.timeout_ms = 1000
    for name, value in fields.items():
        setattr(options.chemistry, name, value)
    result = smsd.MCS.find(query).with_options(options).match(target)
    assert result.status == smsd.SearchStatus.COMPLETED and result.optimal
    return result.best.atom_count if result.found else 0


bond = smsd.BondOrderMode
rows = [
    ("propene, propane: STRICT", common_atoms("C=CC", "CCC")),
    ("propene, propane: LOOSE",
     common_atoms("C=CC", "CCC", match_bond_order=bond.LOOSE)),
    ("propene, propane: ANY",
     common_atoms("C=CC", "CCC", match_bond_order=bond.ANY)),
    ("benzene, cyclohexa-1,3-diene: FLEXIBLE",
     common_atoms("c1ccccc1", "C1=CCCC=C1")),
    ("benzene, cyclohexa-1,3-diene: STRICT",
     common_atoms("c1ccccc1", "C1=CCCC=C1",
                  aromaticity_mode=smsd.AromaticityMode.STRICT)),
    ("acetate, acetic acid: default",
     common_atoms("CC(=O)[O-]", "CC(=O)O")),
    ("acetate, acetic acid: match_formal_charge",
     common_atoms("CC(=O)[O-]", "CC(=O)O", match_formal_charge=True)),
    ("13C and 12C ethanol: default", common_atoms("[13CH3]CO", "[12CH3]CO")),
    ("13C and 12C ethanol: match_isotope",
     common_atoms("[13CH3]CO", "[12CH3]CO", match_isotope=True)),
]
assert [atoms for _, atoms in rows] == [2, 3, 3, 6, 0, 4, 3, 3, 2]
for label, atoms in rows:
    print(f"{label}: {atoms} atoms")
```

The program prints:

```text
propene, propane: STRICT: 2 atoms
propene, propane: LOOSE: 3 atoms
propene, propane: ANY: 3 atoms
benzene, cyclohexa-1,3-diene: FLEXIBLE: 6 atoms
benzene, cyclohexa-1,3-diene: STRICT: 0 atoms
acetate, acetic acid: default: 4 atoms
acetate, acetic acid: match_formal_charge: 3 atoms
13C and 12C ethanol: default: 3 atoms
13C and 12C ethanol: match_isotope: 2 atoms
```

`BondOrderMode.STRICT`, the default, needs equal bond orders; `LOOSE` and
`ANY` relax that rule. `AromaticityMode.FLEXIBLE`, the default, lets an
aromatic bond match a single or double bond, which is why benzene shares six
atoms with cyclohexa-1,3-diene until aromaticity is `STRICT`. Formal charge
and isotope are ignored unless requested; an isotope is compared only where
both atoms carry a mass number.

Three ring options restrict how rings may map, and each does so differently:

```python
import smsd


def common_graph(query, target, **fields):
    options = smsd.MCS.Options()
    options.timeout_ms = 1000
    for name, value in fields.items():
        setattr(options.chemistry, name, value)
    result = smsd.MCS.find(query).with_options(options).match(target)
    assert result.status == smsd.SearchStatus.COMPLETED and result.optimal
    if not result.found:
        return 0, ""
    return result.best.atom_count, result.best.to_smiles(query)


naphthalene, styrene = "c1ccc2ccccc2c1", "C=Cc1ccccc1"
rows = [
    ("naphthalene, styrene: default", common_graph(naphthalene, styrene)),
    ("naphthalene, styrene: complete_rings_only",
     common_graph(naphthalene, styrene, complete_rings_only=True)),
    ("naphthalene, styrene: ring_matches_ring_only",
     common_graph(naphthalene, styrene, ring_matches_ring_only=True)),
    ("naphthalene, styrene: ring_fusion_mode STRICT",
     common_graph(naphthalene, styrene,
                  ring_fusion_mode=smsd.RingFusionMode.STRICT)),
    ("cyclohexane, hexane: default", common_graph("C1CCCCC1", "CCCCCC")),
    ("cyclohexane, hexane: ring_matches_ring_only",
     common_graph("C1CCCCC1", "CCCCCC", ring_matches_ring_only=True)),
]
assert [value for _, value in rows] == [
    (8, "C(C)c1ccccc1"), (6, "c1ccccc1"), (6, "c1ccccc1"),
    (6, "C(CCC)CC"), (5, "C(CC)CC"), (0, ""),
]
for label, (atoms, smiles) in rows:
    print(f"{label}: {atoms} atoms {smiles or '(none)'}")
```

The program prints:

```text
naphthalene, styrene: default: 8 atoms C(C)c1ccccc1
naphthalene, styrene: complete_rings_only: 6 atoms c1ccccc1
naphthalene, styrene: ring_matches_ring_only: 6 atoms c1ccccc1
naphthalene, styrene: ring_fusion_mode STRICT: 6 atoms C(CCC)CC
cyclohexane, hexane: default: 5 atoms C(CC)CC
cyclohexane, hexane: ring_matches_ring_only: 0 atoms (none)
```

- `complete_rings_only` keeps a ring only when the whole ring maps.
- `ring_matches_ring_only` lets ring atoms and bonds match only ring atoms and
  bonds.
- `ring_fusion_mode` `STRICT` needs matched ring atoms to lie in the same
  number of rings. It does not require complete rings, so the result can be a
  chain cut from the rings, as above.

Stereochemistry is a hard constraint. With `use_chirality` or
`use_bond_stereo`, an MCIS search leaves out atoms rather than mapping centres
of opposite configuration, and an MCES search drops the fewest bonds. A centre
without a stated configuration matches either form.

```python
import smsd


def search(query, target, objective=smsd.MCS.Objective.MCIS, **fields):
    options = smsd.MCS.Options()
    options.timeout_ms = 1000
    options.objective = objective
    for name, value in fields.items():
        setattr(options.chemistry, name, value)
    result = smsd.MCS.find(query).with_options(options).match(target)
    assert result.status == smsd.SearchStatus.COMPLETED and result.optimal
    return result.best.atom_count, result.best.bond_count


l_ala, d_ala = "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"  # L- and D-alanine
trans, cis = "C/C=C/C", "C/C=C\\C"  # trans- and cis-2-butene
r_thf, s_thf = "C[C@@H]1CCCO1", "C[C@H]1CCCO1"  # 2-methyltetrahydrofuran
mces = smsd.MCS.Objective.MCES
rows = [
    ("L-, D-alanine: default", search(l_ala, d_ala)),
    ("L-, D-alanine: use_chirality",
     search(l_ala, d_ala, use_chirality=True)),
    ("trans-, cis-2-butene: default", search(trans, cis)),
    ("trans-, cis-2-butene: use_bond_stereo",
     search(trans, cis, use_bond_stereo=True)),
    ("2-methyltetrahydrofuran enantiomers: MCIS",
     search(r_thf, s_thf, use_chirality=True)),
    ("2-methyltetrahydrofuran enantiomers: MCES",
     search(r_thf, s_thf, mces, use_chirality=True)),
    ("L-alanine, unspecified alanine",
     search(l_ala, "NC(C)C(=O)O", use_chirality=True)),
]
assert [value for _, value in rows] == [
    (6, 5), (5, 4), (4, 3), (3, 2), (5, 5), (6, 5), (6, 5),
]
for label, (atoms, bonds) in rows:
    print(f"{label}: {atoms} atoms, {bonds} bonds")
```

The program prints:

```text
L-, D-alanine: default: 6 atoms, 5 bonds
L-, D-alanine: use_chirality: 5 atoms, 4 bonds
trans-, cis-2-butene: default: 4 atoms, 3 bonds
trans-, cis-2-butene: use_bond_stereo: 3 atoms, 2 bonds
2-methyltetrahydrofuran enantiomers: MCIS: 5 atoms, 5 bonds
2-methyltetrahydrofuran enantiomers: MCES: 6 atoms, 5 bonds
L-alanine, unspecified alanine: 6 atoms, 5 bonds
```

Tautomer matching is covered under [Tautomers](#tautomers).

## Atom weights

`with_atom_weights(weights)` ranks MCIS results by the sum of the weights of
the matched query atoms, and then by atoms and bonds. Give exactly one finite
weight per query atom, in the query's atom order. `score` is then in the units
of the weights, and a negative weight can make a smaller common graph win.
Weights apply to MCIS only.

```python
import math

import smsd

result = smsd.MCS.find("CCO").with_atom_weights([0.5, 1.0, 3.0]).match("CCN")
assert result.score == 1.5 and result.best.atom_count == 2
print("ethanol, ethylamine weighted score:", result.score)

query, target = "OCCCN", "OCC.CCCN"  # 3-aminopropanol; ethanol and propylamine
plain = smsd.MCS.find(query).match(target)
heavy = smsd.MCS.find(query).with_atom_weights([5.0, 1.0, 1.0, 1.0, 1.0])
weighted = heavy.match(target)
negative = smsd.MCS.find("CCO").with_atom_weights([1.0, 1.0, -5.0]).match("CCO")
assert (plain.score, plain.best.to_smiles(query)) == (4.0, "C(CC)N")
assert (weighted.score, weighted.best.to_smiles(query)) == (7.0, "C(C)O")
assert (negative.score, negative.best.to_smiles("CCO")) == (2.0, "CC")
print("unweighted:", plain.score, plain.best.to_smiles(query))
print("oxygen weighted 5:", weighted.score, weighted.best.to_smiles(query))
print("oxygen weighted -5:", negative.score, negative.best.to_smiles("CCO"))

mces = smsd.MCS.Objective.MCES
messages = []
for weights, objective in (([1, 1, 1], mces),
                           ([1, 1], smsd.MCS.Objective.MCIS),
                           ([1, math.nan, 1], smsd.MCS.Objective.MCIS)):
    try:
        smsd.MCS.find("CCO").with_objective(objective).with_atom_weights(weights)
    except ValueError as error:
        messages.append(str(error))
assert messages == [
    "Atom weights are supported only for MCIS",
    "Atom weights must cover every query atom",
    "Atom weights must be finite",
]
for message in messages:
    print("ValueError:", message)
```

The program prints:

```text
ethanol, ethylamine weighted score: 1.5
unweighted: 4.0 C(CC)N
oxygen weighted 5: 7.0 C(C)O
oxygen weighted -5: 2.0 CC
ValueError: Atom weights are supported only for MCIS
ValueError: Atom weights must cover every query atom
ValueError: Atom weights must be finite
```

Weights whose positive or negative total is beyond about 2.1 million also
raise `ValueError`. Atom weights cannot be combined with `MCS.consensus`.

## Reserved atoms and constrained placement

`MCS.Options.excluded_target_atoms` makes target atoms unavailable to a search
without editing the target. The indices follow the caller's numbering.
`MCS.match_constrained(queries, target, options=None)` places several queries
on one target in turn: each search excludes the target atoms that earlier
queries used, and the indices in every result refer to the original target.

```python
import smsd

query, target = "CO", "COCCO"  # methanol; 2-methoxyethanol
free = smsd.MCS.find(query).match(target)
options = smsd.MCS.Options()
options.excluded_target_atoms = [3, 4]
reserved = smsd.MCS.find(query).with_options(options).match(target)
assert sorted(free.best.atom_mapping.items()) == [(0, 3), (1, 4)]
assert sorted(reserved.best.atom_mapping.items()) == [(0, 0), (1, 1)]
assert reserved.optimal
print("no exclusions:", sorted(free.best.atom_mapping.items()))
print("atoms 3 and 4 reserved:", sorted(reserved.best.atom_mapping.items()))

placed = smsd.MCS.match_constrained(["CO", "CN"], "CO.CN")
mappings = [sorted(result.best.atom_mapping.items()) for result in placed]
assert mappings == [[(0, 0), (1, 1)], [(0, 2), (1, 3)]]
print("CO, then CN on CO.CN:", mappings)

placed = smsd.MCS.match_constrained(["CC", "CC", "CC"], "CCCC")
mappings = [sorted(result.best.atom_mapping.items()) if result.found else None
            for result in placed]
assert mappings == [[(0, 0), (1, 1)], [(0, 3), (1, 2)], None]
assert all(result.status == smsd.SearchStatus.COMPLETED for result in placed)
print("three CC on butane:", mappings)
```

The program prints:

```text
no exclusions: [(0, 3), (1, 4)]
atoms 3 and 4 reserved: [(0, 0), (1, 1)]
CO, then CN on CO.CN: [[(0, 0), (1, 1)], [(0, 2), (1, 3)]]
three CC on butane: [[(0, 0), (1, 1)], [(0, 3), (1, 2)], None]
```

Each result is proved optimal for its own query, given the atoms that earlier
queries took. In the last case butane has no free pair of bonded carbons left
for the third query. The joint placement is not proved optimal and can depend
on the order of the queries.

## Common graph of a series

`MCS.consensus(molecules, options=None, threshold=1.0)` looks for a common
graph present in at least `ceil(threshold × len(molecules))` of the inputs,
under one shared time limit. The result is an approximate common graph, not
proved optimal: `optimal` is always false, and `complete` says only that the
pass finished. `threshold` must lie in (0, 1]; atom weights and excluded
atoms raise `ValueError`.

```python
import smsd

series = [
    "c1ccccc1CCN",  # phenethylamine
    "c1ccccc1CCO",  # 2-phenylethanol
    "c1ccccc1CC(=O)O",  # phenylacetic acid
    "Cc1ccccc1CCN",  # 2-methylphenethylamine
]
common = smsd.MCS.consensus(series)
assert common.found and common.to_smiles() == "C(C)c1ccccc1"
assert common.support_count == 4 and common.graph.n == 8
assert not common.optimal and common.complete
assert sorted(common.mappings) == [0, 1, 2, 3]
print("common graph:", common.to_smiles(), "| support:", common.support_count,
      "| optimal:", common.optimal, "| complete:", common.complete)
print("inputs with a mapping:", sorted(common.mappings))

partial = smsd.MCS.consensus(["CCO", "CCN", "O"], threshold=2 / 3)
assert partial.to_smiles() == "CC" and partial.support_count == 2
assert sorted(partial.mappings) == [0, 1]
print("two of three:", partial.to_smiles(), "| supported by", sorted(partial.mappings))

empty = smsd.MCS.consensus([])
assert not empty.found and empty.source_index == -1
print("no inputs:", empty.found, empty.source_index)
```

The program prints:

```text
common graph: C(C)c1ccccc1 | support: 4 | optimal: False | complete: True
inputs with a mapping: [0, 1, 2, 3]
two of three: CC | supported by [0, 1]
no inputs: False -1
```

`mappings` maps each supporting input to a `Match` from the atoms of `graph` to
that input's atoms. `fragment` selects the same graph in the input numbered
`source_index`.

## Validate a mapping

`MCS.validate_mapping(query, target, match, options=None)` checks a `Match` or
a plain atom mapping against the objective, connectivity, fragment limits,
excluded atoms and chemistry, and returns a list of errors; an empty list
means the mapping is valid. Without `options` it applies the default rules,
connected MCIS, so validate a result with the options used for the search.
Pass the `Match` itself for an MCES result, because a plain atom mapping
selects every compatible bond between the mapped atoms.

```python
import smsd

query, target = "c1ccc(O)cc1", "c1ccc(N)cc1"  # phenol, aniline
best = smsd.MCS.find(query).match(target).best
checks = [
    ("search result", smsd.MCS.validate_mapping(query, target, best)),
    ("oxygen on nitrogen", smsd.MCS.validate_mapping(query, target, {4: 4})),
    ("target atom twice",
     smsd.MCS.validate_mapping(query, target, {0: 0, 1: 0})),
    ("index out of range", smsd.MCS.validate_mapping(query, target, {0: 99})),
    ("propane on cyclopropane",
     smsd.MCS.validate_mapping("CCC", "C1CC1", {0: 0, 1: 1, 2: 2})),
    ("two ethyl groups",
     smsd.MCS.validate_mapping("CCOCC", "CCNCC", {0: 0, 1: 1, 3: 3, 4: 4})),
]
fragments = smsd.MCS.Options()
fragments.connectivity = smsd.MCS.Connectivity.FRAGMENTS
checks.append(("two ethyl groups, fragments", smsd.MCS.validate_mapping(
    "CCOCC", "CCNCC", {0: 0, 1: 1, 3: 3, 4: 4}, fragments)))

mces = smsd.MCS.Options()
mces.objective = smsd.MCS.Objective.MCES
edges = smsd.MCS.find("C1CC1").with_options(mces).match("CC(C)C").best
checks.append(("MCES result, MCES rules",
               smsd.MCS.validate_mapping("C1CC1", "CC(C)C", edges, mces)))
checks.append(("MCES result, default rules",
               smsd.MCS.validate_mapping("C1CC1", "CC(C)C", edges)))
assert [errors for _, errors in checks] == [
    [], ["Mapped atoms are incompatible"],
    ["Target atom occurs more than once"], ["Mapped atom index out of range"],
    ["Mapped topology is not induced"], ["Selected graph is disconnected"],
    [], [],
    ["Mapped topology is not induced", "MCIS must select every induced bond"],
]
for label, errors in checks:
    print(f"{label}: {errors}")

maximal = smsd.MCS.is_mapping_maximal("CCCO", "CCC.CO", {2: 3, 3: 4})
optimum = smsd.MCS.find("CCCO").match("CCC.CO").best.atom_count
assert maximal and optimum == 3
print("2-atom mapping is maximal:", maximal, "| optimum:", optimum, "atoms")
```

The program prints:

```text
search result: []
oxygen on nitrogen: ['Mapped atoms are incompatible']
target atom twice: ['Target atom occurs more than once']
index out of range: ['Mapped atom index out of range']
propane on cyclopropane: ['Mapped topology is not induced']
two ethyl groups: ['Selected graph is disconnected']
two ethyl groups, fragments: []
MCES result, MCES rules: []
MCES result, default rules: ['Mapped topology is not induced', 'MCIS must select every induced bond']
2-atom mapping is maximal: True | optimum: 3 atoms
```

The possible errors are `Mapped atom index out of range`, `Target atom occurs
more than once`, `Mapped target atom is excluded`, `Mapped atoms are
incompatible`, `Selected bond is outside the atom mapping`, `Selected bond
occurs more than once`, `Selected bonds are incompatible`, `Mapped topology is
not induced`, `MCIS must select every induced bond`, `MCES cannot retain
isolated atoms`, `Selected graph is disconnected`, `Too many fragments`,
`Fragment is too small`, `Selected stereochemistry is incompatible` and
`Selected ring is incomplete`. An excluded atom index outside the target
raises `IndexError`.

`MCS.is_mapping_maximal(query, target, match, options=None)` returns true when
a valid mapping cannot be extended by any further compatible atom. That is
weaker than optimal: the 2-atom mapping above cannot grow, yet the optimum has
three atoms. It returns false for an invalid mapping and raises `RuntimeError`
when its proof cannot finish.

## MCS similarity

The MCS similarity of two molecules is the atom-count Tanimoto of their common
graph: common atoms divided by query atoms plus target atoms minus common
atoms. No function returns this exact value; compute it from
`best.atom_count`. `smsd.similarity_upper_bound(mol1, mol2)` returns a fast
upper bound on it, and `smsd.screen_targets(query, targets, threshold)` keeps
every target whose bound reaches the threshold, so it never rejects a target
that qualifies. Confirm the targets it keeps with an MCS search.

```python
import smsd

query = "c1ccccc1C(=O)O"  # benzoic acid
library = [
    "c1ccccc1C(=O)OC",  # methyl benzoate
    "c1ccc(O)cc1",  # phenol
    "CC(=O)O",  # acetic acid
    "c1ccncc1C(=O)O",  # nicotinic acid
    "CCCC",  # butane
]
prepared = smsd.MCS.find(query).with_timeout(1000)
results = prepared.match_many(library, num_threads=1)
query_atoms = smsd.parse_smiles(query).n
exact, bound = [], []
for smiles, result in zip(library, results):
    assert result.optimal
    common = result.best.atom_count if result.found else 0
    target_atoms = smsd.parse_smiles(smiles).n
    exact.append(round(common / (query_atoms + target_atoms - common), 3))
    bound.append(round(smsd.similarity_upper_bound(query, smiles), 3))
assert exact == [0.9, 0.6, 0.444, 0.8, 0.444]
assert bound == [0.9, 0.778, 0.444, 0.8, 0.444]
for smiles, value, limit in zip(library, exact, bound):
    print(f"{smiles:16} MCS similarity {value:.3f} | bound {limit:.3f}")

kept = smsd.screen_targets(query, library, 0.7, num_threads=1)
confirmed = [index for index in kept if exact[index] >= 0.7]
assert kept == [0, 1, 3] and confirmed == [0, 3]
print("bound keeps:", kept, "| MCS confirms:", confirmed)
```

The program prints:

```text
c1ccccc1C(=O)OC  MCS similarity 0.900 | bound 0.900
c1ccc(O)cc1      MCS similarity 0.600 | bound 0.778
CC(=O)O          MCS similarity 0.444 | bound 0.444
c1ccncc1C(=O)O   MCS similarity 0.800 | bound 0.800
CCCC             MCS similarity 0.444 | bound 0.444
bound keeps: [0, 1, 3] | MCS confirms: [0, 3]
```

`smsd.similarity(mol1, mol2)` is the same upper bound, not a fingerprint
score. For phenol and aspirin the bound happens to equal the exact MCS
similarity, while a fingerprint Tanimoto of the same pair is far lower:

```python
import smsd

phenol, aspirin = "c1ccccc1O", "CC(=O)Oc1ccccc1C(=O)O"
bound = smsd.similarity(phenol, aspirin)
assert bound == smsd.similarity_upper_bound(phenol, aspirin)
result = smsd.MCS.find(phenol).with_timeout(1000).match(aspirin)
assert result.optimal
common = result.best.atom_count
total = smsd.parse_smiles(phenol).n + smsd.parse_smiles(aspirin).n
ecfp4 = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
fingerprint = smsd.tanimoto_coefficient(
    smsd.fingerprint(phenol, options=ecfp4),
    smsd.fingerprint(aspirin, options=ecfp4))
values = (round(bound, 3), round(common / (total - common), 3),
          round(fingerprint, 3))
assert values == (0.538, 0.538, 0.02)
print("smsd.similarity: %.3f | exact MCS similarity: %.3f "
      "| ECFP4 Tanimoto: %.3f" % values)
```

The program prints:

```text
smsd.similarity: 0.538 | exact MCS similarity: 0.538 | ECFP4 Tanimoto: 0.020
```

For a screen of a whole library, `MCS.match_corpus` decides for each target
whether its MCS reaches a minimum atom count, and targets that cannot reach it
are ruled out before searching, with a proof. See
[MCS-similarity screen](SCREENING.md#mcs-similarity-screen).

## Export the common graph

`Match.to_smiles(query)` writes the common graph as SMILES with the cut
valences capped, and `Match.to_graph(query)` returns it as a `MolGraph` that
keeps the original hydrogen counts. The graph can be written as a MOL block
with `smsd.write_mol_block` or drawn with `smsd.depict_svg`.

```python
import smsd

query, target = "C1CC1", "CC(C)C"  # cyclopropane, isobutane
edges = smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES)
best = edges.match(target).best
graph = best.to_graph(query)
counts = smsd.write_mol_block(graph).splitlines()[3]
assert best.to_smiles(query) == "C(C)C"
assert smsd.to_smiles(graph) == "[CH2]C[CH2]"
assert counts.split()[:2] == ["3", "2"]
print("to_smiles:", best.to_smiles(query))
print("to_graph:", smsd.to_smiles(graph))
print("MOL block atoms and bonds:", counts.split()[:2])
```

The program prints:

```text
to_smiles: C(C)C
to_graph: [CH2]C[CH2]
MOL block atoms and bonds: ['3', '2']
```

`match_with_smarts(target)` runs the search and returns a `SmartsResult` with
the result and the common graph as a SMARTS pattern under the query's
chemistry; the pattern is empty when nothing is found. The pattern can then
search a library. Under the default `FLEXIBLE` aromaticity its bonds also
match single and double bonds, so a benzene pattern finds cyclohexane; choose
`AromaticityMode.STRICT` for an aromatic-only pattern. Stereo, tautomer,
complete-ring and strict ring-fusion rules cannot be written as SMARTS, so
they raise `ValueError` rather than being dropped.

```python
import smsd

library = ["c1ccccc1O", "c1ccccc1N", "Cc1ccccc1", "CCO", "C1CCCCC1"]
exported = smsd.MCS.find("c1ccccc1O").with_timeout(1000).match_with_smarts(
    "c1ccccc1N")  # phenol, aniline
pattern = exported.smarts
assert exported.result.optimal
assert pattern == "[#6]-,=,:1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]1"
search = smsd.Substructure.smarts(pattern).with_timeout(1000)
hits = [smiles for smiles, result
        in zip(library, search.match_many(library, num_threads=1))
        if result.found]
assert hits == ["c1ccccc1O", "c1ccccc1N", "Cc1ccccc1", "C1CCCCC1"]
print("pattern:", pattern)
print("library hits:", hits)

strict = smsd.MCS.Options()
strict.chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
aromatic = smsd.MCS.find("c1ccccc1O").with_options(strict).match_with_smarts(
    "c1ccccc1N").smarts
assert aromatic == "[#6;a]:1:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]1"
print("STRICT aromaticity:", aromatic)

stereo = smsd.MCS.Options()
stereo.chemistry.use_chirality = True
try:
    smsd.MCS.find("N[C@@H](C)C(=O)O").with_options(stereo).match_with_smarts(
        "N[C@@H](C)C(=O)O")
except ValueError as error:
    message = str(error)
else:
    raise AssertionError("a chiral SMARTS export was accepted")
assert message == ("SMARTS export cannot represent these stereo, tautomer or "
                   "ring constraints")
print("use_chirality:", message)

nothing = smsd.MCS.find("O").match_with_smarts("C")
assert nothing.smarts == "" and not nothing.result.found
print("no common atom:", repr(nothing.smarts))
```

The program prints:

```text
pattern: [#6]-,=,:1-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]-,=,:[#6]1
library hits: ['c1ccccc1O', 'c1ccccc1N', 'Cc1ccccc1', 'C1CCCCC1']
STRICT aromaticity: [#6;a]:1:[#6;a]:[#6;a]:[#6;a]:[#6;a]:[#6;a]1
use_chirality: SMARTS export cannot represent these stereo, tautomer or ring constraints
no common atom: ''
```

The pattern carries neither the induced rule of MCIS nor the proof of the
search. Writing patterns by hand is covered in [SMARTS queries](SMARTS.md).

## Draw the result

`smsd.highlight.mcs_svg(query, target, result, *, opts=None, **kwargs)` draws
both molecules with the common graph highlighted, without RDKit and without a
new search. It takes the `SearchResult`, not a `Match`, and records the
status, `optimal` and `complete` in the SVG, both as attributes of the root
element and in its title, so a figure cannot hide a stopped search. Import it
from `smsd.highlight`; `import smsd` alone does not provide it. Keyword
arguments are the fields of `smsd.DepictOptions`, such as `width` and
`height`.

```python
import re
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd
from smsd.highlight import mcs_svg

query, target = "c1ccccc1O", "c1ccccc1N"  # phenol, aniline
result = smsd.MCS.find(query).with_timeout(1000).match(target)
svg = mcs_svg(query, target, result, width=600, height=300)
title = re.search(r"<title>(.*?)</title>", svg).group(1)
flags = re.findall(r'data-(search-status|optimal|complete)="(\w+)"', svg)
assert title == ("MCS: 6 mapped atoms; status COMPLETED; optimal true; "
                 "complete true")
assert flags == [("search-status", "COMPLETED"), ("optimal", "true"),
                 ("complete", "true")]
print(title)
print(flags)

with TemporaryDirectory() as directory:
    path = Path(directory) / "phenol_aniline.svg"
    smsd.save_svg(svg, path)
    assert path.read_text(encoding="utf-8") == svg
    print("saved", path.name)

try:
    mcs_svg(query, target, result.best)
except TypeError as error:
    message = str(error)
else:
    raise AssertionError("a Match was accepted")
assert message == "result must be a SearchResult"
print("TypeError:", message)
```

The program prints:

```text
MCS: 6 mapped atoms; status COMPLETED; optimal true; complete true
[('search-status', 'COMPLETED'), ('optimal', 'true'), ('complete', 'true')]
saved phenol_aniline.svg
TypeError: result must be a SearchResult
```

An MCES result can leave out a bond between two mapped atoms, which the pair
drawing cannot show, so `mcs_svg` raises `ValueError` for it. Draw the
selected graph from `to_graph` instead:

```python
import smsd
from smsd.highlight import mcs_svg

query, target = "C1CC1", "CC(C)C"  # cyclopropane, isobutane
edges = smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES)
result = edges.match(target)
try:
    mcs_svg(query, target, result)
except ValueError:
    rejected = True
else:
    rejected = False
svg = smsd.depict_svg(result.best.to_graph(query), width=300, height=200)
assert rejected and "<svg" in svg
print("mcs_svg rejects the MCES result:", rejected)
print("selected graph drawn:", "<svg" in svg)
```

The program prints:

```text
mcs_svg rejects the MCES result: True
selected graph drawn: True
```

`smsd.depict_pair(mol1, mol2, mapping)` draws two molecules with any mapping
from atoms of `mol1` to atoms of `mol2` highlighted, and
`smsd.depict_mapping(mol, mapping)` highlights the keys of a mapping on one
molecule, so invert an MCS mapping to highlight target atoms. Both return
native SVG. `smsd.depict_mcs(mol1, mol2, *, tautomer_aware=False,
timeout_ms=10000, size=(600, 300))` returns an image drawn with RDKit. It runs
its own MCIS search and shows no status, so prefer `mcs_svg` for any figure
that should state whether the result is proved. Its `highlight_color` argument
has no effect and is deprecated. The drawing options are listed under
[Depiction](../PYTHON.md#depiction).

## Tautomers

With `tautomer_aware`, atoms and bonds that differ only by the position of a
mobile hydrogen can correspond. Such matches usually also need a relaxed bond
order, so `ChemOptions.tautomer_profile()` sets `tautomer_aware` with
`BondOrderMode.LOOSE`.

```python
import smsd

query = smsd.parse_smiles("ONC(C)=O")  # acetohydroxamic acid
target = smsd.parse_smiles("ON=C(C)O")  # its imidic acid tautomer
options = smsd.MCS.Options()
options.timeout_ms = 1000
options.chemistry = smsd.ChemOptions.tautomer_profile()
result = smsd.MCS.find(query).with_options(options).match(target)
match = result.best
assert result.optimal and match.atom_count == 5
assert smsd.MCS.validate_tautomer_consistency(query, target, match)
print("common atoms:", match.atom_count, "| optimal:", result.optimal)
print("protons conserved:",
      smsd.MCS.validate_tautomer_consistency(query, target, match))

scores = []
for pH in (2.0, 7.4, 9.0, 12.0):
    chemistry = smsd.ChemOptions.tautomer_profile().with_ph(pH)
    score = smsd.MCS.compute_tautomer_confidence(query, target, match, chemistry)
    scores.append((pH, round(score, 3)))
assert scores == [(2.0, 0.85), (7.4, 0.829), (9.0, 0.425), (12.0, 0.1)]
for pH, score in scores:
    print(f"pH {pH}: confidence {score}")

keto = smsd.MCS.Options()
keto.chemistry.tautomer_aware = True
keto.chemistry.match_bond_order = smsd.BondOrderMode.ANY
default = smsd.MCS.find("CC(=O)C").match("CC(O)=C")  # acetone, its enol
relaxed = smsd.MCS.find("CC(=O)C").with_options(keto).match("CC(O)=C")
assert (default.best.atom_count, relaxed.best.atom_count) == (2, 4)
print("acetone, prop-1-en-2-ol:", default.best.atom_count, "atoms by default,",
      relaxed.best.atom_count, "with tautomer_aware and ANY")
```

The program prints:

```text
common atoms: 5 | optimal: True
protons conserved: True
pH 2.0: confidence 0.85
pH 7.4: confidence 0.829
pH 9.0: confidence 0.425
pH 12.0: confidence 0.1
acetone, prop-1-en-2-ol: 2 atoms by default, 4 with tautomer_aware and ANY
```

`MCS.validate_tautomer_consistency(query, target, match)` reports whether a
`Match` or atom mapping conserves the mobile N–H, O–H and S–H protons of each
mapped tautomeric group. `MCS.compute_tautomer_confidence(query, target,
match, chemistry)` scores the tautomeric relaxations a `Match` uses at the
`pH` and `solvent` of `chemistry`. The score is a heuristic in [0, 1], not a
probability, and it changes neither the molecules nor the search. It needs a
`Match`; a plain atom mapping raises `TypeError`.

## Options

`MCS` methods:

| Method | Returns |
| --- | --- |
| `MCS.find(query)` | A prepared query with `MCS.Options()` |
| `with_options(options)` | A new query holding a copy of `options` |
| `with_timeout(milliseconds)` | A new query with that time limit; negative selects the adaptive limit |
| `with_objective(objective)` | A new query with `MCS.Objective.MCIS` or `MCES` |
| `with_fragments(*, min_atoms=1, max_fragments=2147483647)` | A new query allowing fragments |
| `with_atom_weights(weights)` | A new MCIS query ranked by atom weights |
| `match(target)` | A `SearchResult` with at most one `Match` |
| `match_all(target, limit=1000)` | Every optimal mapping up to `limit` |
| `match_many(targets, *, num_threads=0)` | One `SearchResult` per target, in order, each with its own time limit |
| `match_with_smarts(target)` | A `SmartsResult(result, smarts)` |
| `match_corpus(corpus, minimum_atoms, *, gpu_prefilter=False, num_threads=0)` | One `CorpusSearchResult` per target of a `TargetCorpus` |
| `MCS.validate_mapping(query, target, match, options=None)` | A list of errors, empty when valid |
| `MCS.is_mapping_maximal(query, target, match, options=None)` | Whether no compatible extension exists |
| `MCS.canonicalize_mapping(query, target, atom_mapping)` | The representative mapping under both molecules' symmetries |
| `MCS.validate_tautomer_consistency(query, target, match)` | Whether mobile protons are conserved |
| `MCS.compute_tautomer_confidence(query, target, match, chemistry)` | A heuristic score in [0, 1] |
| `MCS.match_constrained(queries, target, options=None)` | One `SearchResult` per query, in order |
| `MCS.consensus(molecules, options=None, threshold=1.0)` | A `ConsensusResult` |

`num_threads=0` uses the OpenMP default, which `OMP_NUM_THREADS` sets; a
negative value raises `ValueError`. `limit` must be positive.

`MCS.Options` fields:

| Field | Type | Default | Meaning |
| --- | --- | --- | --- |
| `chemistry` | `ChemOptions` | `ChemOptions()` | Atom and bond matching rules |
| `timeout_ms` | `int` | `-1` | Time limit per pair in milliseconds; negative selects min(30,000, 500 + 2 × query atoms × target atoms); `0` stops at once |
| `node_limit` | `int` | `0` | Search nodes allowed; `0` sets no limit |
| `cancellation` | `CancellationToken` or `None` | `None` | Stops the search with `CANCELLED` once cancelled |
| `objective` | `MCS.Objective` | `MCIS` | `MCIS` maximises atoms, `MCES` maximises bonds |
| `connectivity` | `MCS.Connectivity` | `CONNECTED` | `CONNECTED` or `FRAGMENTS` |
| `min_fragment_atoms` | `int` | `1` | Smallest fragment allowed with `FRAGMENTS` |
| `max_fragments` | `int` | `2147483647` | Most fragments allowed with `FRAGMENTS` |
| `atom_weights` | `list[float]` | `[]` | One weight per query atom; MCIS only |
| `excluded_target_atoms` | `list[int]` | `[]` | Target atoms unavailable to the search |

The default time limit is 1.3 s for two 20-atom molecules, 5.5 s for two
50-atom molecules and 30 s from 122 atoms each. Parsing and preparation fall
outside it.

## Notes

- Report `status`, `optimal` and `complete` together with any result. A
  stopped search can return a valid but smaller common graph, or none.
- Scores mean atoms for MCIS, bonds for MCES and weight units with atom
  weights; never compare them across objectives.
- Tied optimal mappings are not ranked chemically. A stopped MCIS search can
  return a different mapping in another language or version, so do not rely
  on the exact mapping of a stopped search.
- MCES declines very large pairs before searching, with `LIMIT_REACHED` and no
  match; see [Size limits](../PYTHON.md#size-limits).
- `with_options` copies the options when it is called, so later edits to the
  options object do not reach the prepared query. A `CancellationToken` is
  shared, however.
- Other entry points use their own bounded defaults: the `find_mcs` molecular
  tool uses 1,000 ms and 100,000 nodes, the command-line `--timeout` is
  10,000 ms per pair, and `depict_mcs` uses 10,000 ms.
- The MCS search always runs on the CPU. `match_corpus(...,
  gpu_prefilter=True)` affects only the step that rules out targets before
  searching; see
  [Request the GPU for a corpus pre-screen](CPU_GPU.md#request-the-gpu-for-a-corpus-pre-screen).
- R-group decomposition needs the whole core in every target, so it is a
  substructure operation; see
  [R-group decomposition](SUBSTRUCTURE.md#r-group-decomposition).

## See also

- [Substructure search](SUBSTRUCTURE.md) and [SMARTS queries](SMARTS.md)
- [Screening a library](SCREENING.md) and
  [Fingerprints and similarity](FINGERPRINTS.md)
- [MCS on the command line](COMMAND_LINE.md#mcs) and
  [C++ and Java](OTHER_LANGUAGES.md)
- The [Python reference](../PYTHON.md), the [C++ guide](../CPP.md), the
  [Java guide](../JAVA.md) and the [examples](../EXAMPLES.md)
- The notebook [mcs_search.ipynb](../../python/examples/mcs_search.ipynb)

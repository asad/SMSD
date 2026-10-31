# Substructure search

A substructure search asks whether a query occurs whole inside a target, and
where. SMSD prepares the query once and then searches one target, every
embedding in one target, a list of targets or a captured corpus. Every call
returns a `SearchResult` that holds the matches together with the reason the
search ended, so an empty result proves absence only when its status is
`COMPLETED` and `optimal` is true.

This chapter covers molecule queries in Python. The SMARTS syntax is described
in [SMARTS queries](SMARTS.md), and the exhaustive API list is in the
[Python guide](../PYTHON.md#searches).

## Prepare a query

- `Substructure.find(query)` prepares a molecule query from a SMILES string, a
  `MolGraph` or an RDKit molecule. The query graph is copied, so later changes
  to the input do not reach it.
- `Substructure.smarts(pattern)` compiles a SMARTS pattern once.

A prepared query is immutable. `with_options`, `with_timeout` and
`with_induced` return a new query with a copy of the options and leave the
original unchanged. Searches release the global interpreter lock, so several
threads can share one prepared query.

```python
from concurrent.futures import ThreadPoolExecutor

import smsd

# The same phenol query, prepared from SMILES, from a MolGraph and as SMARTS
queries = {
    "SMILES": smsd.Substructure.find("c1ccccc1O"),
    "MolGraph": smsd.Substructure.find(smsd.parse_smiles("c1ccccc1O")),
    "SMARTS": smsd.Substructure.smarts("c1ccccc1[OX2H]"),
}
bounded = queries["SMILES"].with_timeout(500)
assert bounded is not queries["SMILES"]

targets = ["CC(=O)Nc1ccc(O)cc1",          # paracetamol
           "CC(C)Cc1ccc(cc1)C(C)C(=O)O",  # ibuprofen
           "OC(=O)c1ccccc1O"]             # salicylic acid

# One prepared query serves several threads at once
with ThreadPoolExecutor(max_workers=2) as pool:
    for name, query in queries.items():
        found = [result.found for result in pool.map(query.match, targets)]
        print(f"{name:8}", found)
        assert found == [True, False, True]
```

The program prints:

```text
SMILES   [True, False, True]
MolGraph [True, False, True]
SMARTS   [True, False, True]
```

## Find one embedding

`match(target)` returns at most one embedding. `atom_mapping` maps each query
atom index to a target atom index, in each input's own numbering (see
[Getting started](GETTING_STARTED.md#atom-numbering)). `score` is the number of
query atoms when a match is found and `0.0` otherwise.

```python
import smsd

query = smsd.Substructure.find("c1ccccc1O")   # phenol
result = query.match("CC(C)(C)c1ccc(O)cc1")   # 4-tert-butylphenol

print("found:", result.found)
print("status:", result.status.name, "optimal:", result.optimal,
      "complete:", result.complete)
print("score:", result.score)
print("mapping:", sorted(result.best.atom_mapping.items()))
print("matched query atoms as SMILES:", result.best.to_smiles("c1ccccc1O"))
assert result.found and result.status == smsd.SearchStatus.COMPLETED
assert result.optimal and result.complete and result.score == 7.0
assert sorted(result.best.atom_mapping.items()) == [
    (0, 6), (1, 5), (2, 4), (3, 10), (4, 9), (5, 7), (6, 8)]
assert result.best.to_smiles("c1ccccc1O") == "c1(ccccc1)O"

# A completed search that finds nothing proves absence
absent = query.match("Nc1ccccc1")             # aniline
print("aniline:", absent.found, absent.status.name, absent.optimal,
      absent.complete)
assert not absent.found and absent.status == smsd.SearchStatus.COMPLETED
assert absent.optimal and absent.complete
```

The program prints:

```text
found: True
status: COMPLETED optimal: True complete: True
score: 7.0
mapping: [(0, 6), (1, 5), (2, 4), (3, 10), (4, 9), (5, 7), (6, 8)]
matched query atoms as SMILES: c1(ccccc1)O
aniline: False COMPLETED True True
```

## Enumerate embeddings

`match_all(target, limit=1000)` lists every embedding up to `limit`.
Symmetry-related embeddings are listed separately: benzene sits in toluene in
twelve ways, all on the same six ring atoms. To count occurrences, collect the
target atom sets with `frozenset`.

When the limit is reached, the status is `LIMIT_REACHED` and `complete` is
false; `optimal` stays true, because a found embedding already proves that the
query occurs. A limit that is not a positive integer raises `ValueError`. Give
every enumeration on a large target both a limit and a time limit.

```python
import smsd

benzene = smsd.Substructure.find("c1ccccc1").with_timeout(1000)
toluene = "Cc1ccccc1"

every = benzene.match_all(toluene, limit=100)
occurrences = {frozenset(match.atom_mapping.values()) for match in every.matches}
print("embeddings:", len(every.matches), every.status.name,
      "complete:", every.complete)
print("distinct ring occurrences:", len(occurrences))
assert len(every.matches) == 12 and every.complete and len(occurrences) == 1

capped = benzene.match_all(toluene, limit=5)
print("limit 5:", len(capped.matches), capped.status.name,
      "optimal:", capped.optimal, "complete:", capped.complete)
assert capped.status == smsd.SearchStatus.LIMIT_REACHED
assert capped.optimal and not capped.complete

# Two carbonyl groups in aspirin
carbonyl = smsd.Substructure.smarts("C=O").with_timeout(1000)
groups = carbonyl.match_all("CC(=O)Oc1ccccc1C(=O)O", limit=100)
atoms = sorted(sorted(match.atom_mapping.values()) for match in groups.matches)
print("carbonyl atoms in aspirin:", atoms)
assert atoms == [[1, 2], [10, 11]]

try:
    benzene.match_all(toluene, limit=0)
except ValueError as error:
    print("limit 0:", error)
    assert str(error) == "Match limit must be positive"
```

The program prints:

```text
embeddings: 12 COMPLETED complete: True
distinct ring occurrences: 1
limit 5: 5 LIMIT_REACHED optimal: True complete: False
carbonyl atoms in aspirin: [[1, 2], [10, 11]]
limit 0: Match limit must be positive
```

## Search many targets

`match_many(targets, *, num_threads=0)` searches each target independently and
returns one `SearchResult` per target, in input order. SMILES strings,
`MolGraph` objects and RDKit molecules can be mixed in one list. Each target
gets its own time limit and keeps its own atom numbering. `num_threads=0`
uses the OpenMP default, which `OMP_NUM_THREADS` sets.

`match_many` keeps at most one embedding per target. To count occurrences in
each target, call `match_all` on that target.

```python
import smsd

acid = smsd.Substructure.smarts("C(=O)[OX2H1]")      # carboxylic acid
targets = ["CC(=O)Oc1ccccc1C(=O)O",                  # aspirin
           smsd.parse_smiles("CC(=O)Nc1ccc(O)cc1"),  # paracetamol
           "CC(C)Cc1ccc(cc1)C(C)C(=O)O"]             # ibuprofen
results = acid.match_many(targets, num_threads=2)
found = [result.found for result in results]
completed = all(result.status == smsd.SearchStatus.COMPLETED for result in results)
print("carboxylic acid:", found, "all completed:", completed)
assert found == [True, False, True] and completed

benzene = smsd.Substructure.find("c1ccccc1").with_timeout(1000)
pair = ["Cc1ccccc1", "c1ccc2ccccc2c1"]               # toluene, naphthalene
kept = [len(result.matches) for result in benzene.match_many(pair, num_threads=2)]
listed = [len(benzene.match_all(target, limit=100).matches) for target in pair]
print("match_many keeps:", kept, "match_all lists:", listed)
assert kept == [1, 1] and listed == [12, 24]
```

The program prints:

```text
carboxylic acid: [True, False, True] all completed: True
match_many keeps: [1, 1] match_all lists: [12, 24]
```

## Search a prepared corpus

`TargetCorpus.from_smiles(targets)` captures a library once, with each
target's atom numbering, so that many queries can search it. A substructure
search needs no `prewarm` and no NumPy. `match_many` accepts a corpus in place
of a list.

`match_corpus(corpus, *, gpu_prefilter=False, num_threads=0)` returns one
`CorpusSearchResult` for each corpus slot, in order. Targets that cannot
contain the query are ruled out before matching, with a proof; every other
target is searched as usual.

| Field | Meaning |
| --- | --- |
| `target_index` | The slot in the corpus |
| `atom_upper_bound` | Always `None` for a substructure search |
| `proven_below_threshold` | True when the query is proved absent |
| `threshold_reached` | True when an embedding was found |
| `result` | The `SearchResult`, or `None` for a slot ruled out before matching |

When both flags are false, the outcome is unknown, for example after a search
was stopped. Every slot is searched, with no slot ruled out, for a SMARTS
query, with `match_atom_type=False`, with a `node_limit` other than 0, with
`timeout_ms` set to 0 or with a cancellation token. With
`gpu_prefilter=True`, an attached device runs the step that rules targets out,
and the CPU runs it when no device is usable; the decisions are the same either
way. See
[CPU and GPU screening](CPU_GPU.md#request-the-gpu-for-a-corpus-pre-screen).

```python
import smsd

names = ["aspirin", "paracetamol", "ibuprofen", "caffeine", "salicylic acid"]
corpus = smsd.TargetCorpus.from_smiles([
    "CC(=O)Oc1ccccc1C(=O)O",
    "CC(=O)Nc1ccc(O)cc1",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "OC(=O)c1ccccc1O",
])

searches = [("benzene ring", smsd.Substructure.find("c1ccccc1"),
             ["aspirin", "paracetamol", "ibuprofen", "salicylic acid"]),
            ("aromatic hydroxyl", smsd.Substructure.smarts("c[OX2H]"),
             ["paracetamol", "salicylic acid"])]
for label, query, expected in searches:
    results = query.match_many(corpus, num_threads=2)
    hits = [name for name, result in zip(names, results) if result.found]
    print(f"{label}:", hits)
    assert all(result.status == smsd.SearchStatus.COMPLETED for result in results)
    assert hits == expected

decisions = smsd.Substructure.find("c1ccccc1O").match_corpus(corpus, num_threads=2)
outcomes = []
for name, decision in zip(names, decisions):
    if decision.result is None:
        outcome = "ruled out before matching"
    elif decision.threshold_reached:
        outcome = "found"
    elif decision.proven_below_threshold:
        outcome = "searched, proved absent"
    else:
        outcome = "unknown"
    print(f"{name:15}", outcome)
    assert decision.atom_upper_bound is None
    outcomes.append(outcome)
assert outcomes == ["found", "found", "searched, proved absent",
                    "ruled out before matching", "found"]
```

The program prints:

```text
benzene ring: ['aspirin', 'paracetamol', 'ibuprofen', 'salicylic acid']
aromatic hydroxyl: ['paracetamol', 'salicylic acid']
aspirin         found
paracetamol     found
ibuprofen       searched, proved absent
caffeine        ruled out before matching
salicylic acid  found
```

## Induced matching

By default a search is non-induced: query bonds must map to target bonds, but
two query atoms that are not bonded may map to target atoms that are.
`with_induced()` also requires every pair of unbonded query atoms to map to
unbonded target atoms. Propane is then no longer found in cyclopropane, while
cyclohexane is still found in decalin.

```python
import smsd

pairs = [("CCC", "C1CC1", "propane in cyclopropane"),
         ("C1CCCCC1", "C1CCC2CCCCC2C1", "cyclohexane in decalin")]
for query, target, label in pairs:
    plain = smsd.Substructure.find(query).match(target).found
    induced = smsd.Substructure.find(query).with_induced().match(target).found
    print(f"{label:24} default: {plain!s:5} induced: {induced}")
    assert plain and induced == (target != "C1CC1")
```

The program prints:

```text
propane in cyclopropane  default: True  induced: False
cyclohexane in decalin   default: True  induced: True
```

## Chemistry options

`options.chemistry` is a `ChemOptions` object that decides which atoms and
bonds may match. For a SMARTS query, the pattern itself states charge,
hydrogen count, aromaticity, bond order and stereochemistry, so only
`complete_rings_only` applies (see
[SMARTS queries](SMARTS.md#which-options-apply)).

| Attribute | Default | Effect on a molecule query | SMARTS |
| --- | --- | --- | --- |
| `match_atom_type` | `True` | Matched atoms have the same element | No |
| `match_formal_charge` | `False` | Matched atoms have the same formal charge | No |
| `match_isotope` | `False` | Matched atoms have the same mass number where both atoms state one | No |
| `use_chirality` | `False` | A tagged query centre needs the same configuration in the target; an untagged query centre matches any; an untagged target centre does not satisfy a tagged query | No |
| `use_bond_stereo` | `False` | The same rule for double-bond configuration | No |
| `ring_matches_ring_only` | `False` | Ring atoms and bonds match only ring atoms and bonds | No |
| `complete_rings_only` | `False` | Matched ring atoms and bonds lie in fully matched rings of the target | Yes |
| `ring_fusion_mode` | `RingFusionMode.IGNORE` | `STRICT` requires the same number of rings on each matched atom; `PERMISSIVE` behaves as `IGNORE` | No |
| `match_bond_order` | `BondOrderMode.STRICT` | Bond orders agree; `LOOSE` and `ANY` both accept any order and give identical results | No |
| `aromaticity_mode` | `AromaticityMode.FLEXIBLE` | Aromatic atoms and bonds may also match non-aromatic ones; `STRICT` requires agreement | No |
| `tautomer_aware` | `False` | Recognised tautomeric bonds may match; usually combined with a relaxed bond order | No |
| `pH`, `solvent` | `7.4`, `Solvent.AQUEOUS` | Used only by `MCS.compute_tautomer_confidence` | No |

`LOOSE` and `ANY` still respect `aromaticity_mode=STRICT` and
`ring_matches_ring_only`. Options objects take no keyword arguments, so set
each field as an attribute. The program below shows the effect of each option
on a small pair of molecules.

```python
import smsd


def found(query, target, **fields):
    options = smsd.Substructure.Options()
    for name, value in fields.items():
        setattr(options.chemistry, name, value)
    return smsd.Substructure.find(query).with_options(options).match(target).found


BOND = smsd.BondOrderMode
AROMATIC = smsd.AromaticityMode
FUSION = smsd.RingFusionMode
rows = [
    ("C=C in ethane", "C=C", "CC", dict(match_bond_order=BOND.ANY)),
    ("C=C in ethane", "C=C", "CC", dict(match_bond_order=BOND.LOOSE)),
    ("phenol in cyclohexanol", "c1ccccc1O", "OC1CCCCC1",
     dict(aromaticity_mode=AROMATIC.STRICT)),
    ("C[NH3+] in methylamine", "C[NH3+]", "CN", dict(match_formal_charge=True)),
    ("[13CH4] in methane", "[13CH4]", "C", dict(match_isotope=True)),
    ("methanol in methylamine", "CO", "CN", dict(match_atom_type=False)),
    ("L-alanine in D-alanine", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O",
     dict(use_chirality=True)),
    ("L-alanine in racemic alanine", "N[C@@H](C)C(=O)O", "NC(C)C(=O)O",
     dict(use_chirality=True)),
    ("racemic alanine in L-alanine", "NC(C)C(=O)O", "N[C@@H](C)C(=O)O",
     dict(use_chirality=True)),
    ("F/C=C/F in F/C=C\\F", "F/C=C/F", "F/C=C\\F", dict(use_bond_stereo=True)),
    ("propane in cyclohexane", "CCC", "C1CCCCC1",
     dict(ring_matches_ring_only=True)),
    ("propane in cyclohexane", "CCC", "C1CCCCC1", dict(complete_rings_only=True)),
    ("propane in propylcyclohexane", "CCC", "CCCC1CCCCC1",
     dict(complete_rings_only=True)),
    ("cyclohexane in decalin", "C1CCCCC1", "C1CCC2CCCCC2C1",
     dict(ring_fusion_mode=FUSION.STRICT)),
    ("acetone in its enol", "CC(=O)C", "CC(O)=C", dict(tautomer_aware=True)),
    ("acetone in its enol", "CC(=O)C", "CC(O)=C",
     dict(tautomer_aware=True, match_bond_order=BOND.ANY)),
]
expected = [(False, True), (False, True), (True, False), (True, False),
            (True, True), (False, True), (True, False), (True, False),
            (True, True), (True, False), (True, False), (True, False),
            (True, True), (True, False), (False, False), (False, True)]
for (label, query, target, fields), outcome in zip(rows, expected):
    setting = ", ".join(f"{name}={getattr(value, 'name', value)}"
                        for name, value in fields.items())
    default, changed = found(query, target), found(query, target, **fields)
    print(f"{label:29} {setting:40} {default!s:5} -> {changed}")
    assert (default, changed) == outcome
```

The program prints:

```text
C=C in ethane                 match_bond_order=ANY                     False -> True
C=C in ethane                 match_bond_order=LOOSE                   False -> True
phenol in cyclohexanol        aromaticity_mode=STRICT                  True  -> False
C[NH3+] in methylamine        match_formal_charge=True                 True  -> False
[13CH4] in methane            match_isotope=True                       True  -> True
methanol in methylamine       match_atom_type=False                    False -> True
L-alanine in D-alanine        use_chirality=True                       True  -> False
L-alanine in racemic alanine  use_chirality=True                       True  -> False
racemic alanine in L-alanine  use_chirality=True                       True  -> True
F/C=C/F in F/C=C\F            use_bond_stereo=True                     True  -> False
propane in cyclohexane        ring_matches_ring_only=True              True  -> False
propane in cyclohexane        complete_rings_only=True                 True  -> False
propane in propylcyclohexane  complete_rings_only=True                 True  -> True
cyclohexane in decalin        ring_fusion_mode=STRICT                  True  -> False
acetone in its enol           tautomer_aware=True                      False -> False
acetone in its enol           tautomer_aware=True, match_bond_order=ANY False -> True
```

Three rows deserve attention. Under the default `FLEXIBLE` mode, phenol is
found in cyclohexanol; `STRICT` aromaticity prevents this. `match_isotope`
compares mass numbers only where both atoms state one, so the labelled query
still matches unlabelled methane. Tautomer matching usually also needs a
relaxed bond order.

`ChemOptions.profile(name)` returns a named set of options, and the name is
not case-sensitive. Each profile starts from the defaults and changes only the
fields shown below. `strict` changes only the aromaticity mode, so charge and
isotope matching still need `match_formal_charge` or `match_isotope`. An
unknown name raises `ValueError`.

```python
import smsd

FIELDS = ("match_atom_type", "match_formal_charge", "match_isotope",
          "use_chirality", "use_bond_stereo", "ring_matches_ring_only",
          "complete_rings_only", "ring_fusion_mode", "match_bond_order",
          "aromaticity_mode", "tautomer_aware")
default = smsd.ChemOptions()


def changes(chemistry):
    changed = []
    for name in FIELDS:
        value = getattr(chemistry, name)
        if value != getattr(default, name):
            changed.append(f"{name}={getattr(value, 'name', value)}")
    return changed


expected = {
    "default": [],
    "strict": ["aromaticity_mode=STRICT"],
    "pharma": ["match_formal_charge=True", "ring_matches_ring_only=True",
               "complete_rings_only=True"],
    "reaction": ["match_bond_order=LOOSE"],
    "compat-fmcs": ["match_bond_order=LOOSE"],
    "compat-substruct": ["match_bond_order=LOOSE"],
}
for name, fields in expected.items():
    print(f"{name:17}", changes(smsd.ChemOptions.profile(name)))
    assert changes(smsd.ChemOptions.profile(name.upper())) == fields
tautomer = changes(smsd.ChemOptions.tautomer_profile())
print(f"{'tautomer_profile':17}", tautomer)
assert tautomer == ["match_bond_order=LOOSE", "tautomer_aware=True"]

try:
    smsd.ChemOptions.profile("fast")
except ValueError as error:
    print(error)
    assert str(error) == "Unknown ChemOptions profile: fast"
```

The program prints:

```text
default           []
strict            ['aromaticity_mode=STRICT']
pharma            ['match_formal_charge=True', 'ring_matches_ring_only=True', 'complete_rings_only=True']
reaction          ['match_bond_order=LOOSE']
compat-fmcs       ['match_bond_order=LOOSE']
compat-substruct  ['match_bond_order=LOOSE']
tautomer_profile  ['match_bond_order=LOOSE', 'tautomer_aware=True']
Unknown ChemOptions profile: fast
```

## Hydrogen counts

A molecule query does not compare hydrogen counts. Phenol as a molecule query
therefore matches aspirin, whose ester oxygen bears no hydrogen. When hydrogen
counts matter, write SMARTS with an `H` primitive, such as `c[OX2H]`, or turn
the drawn query into SMARTS with `to_smarts(..., include_h_count=True)`. The
second route fixes the hydrogen count of every atom, including the ring
positions, so the pattern below matches unsubstituted phenol only.

```python
import smsd

aspirin = "CC(=O)Oc1ccccc1C(=O)O"
phenol = "Oc1ccccc1"
butylphenol = "CC(C)(C)c1ccc(O)cc1"   # 4-tert-butylphenol

molecule = smsd.Substructure.find("c1ccccc1O").match(aspirin).found
hydroxyl = smsd.Substructure.smarts("c[OX2H]").match(aspirin).found
print("phenol as a molecule query, in aspirin:", molecule)
print("c[OX2H] in aspirin:", hydroxyl)
assert molecule and not hydroxyl

pattern = smsd.to_smarts("c1ccccc1O", include_h_count=True)
print(pattern)
assert pattern.endswith("-[#8;A;H1]") and pattern.count("H1") == 6
fixed = smsd.Substructure.smarts(pattern)
outcome = [fixed.match(target).found for target in (aspirin, phenol, butylphenol)]
print("aspirin, phenol, 4-tert-butylphenol:", outcome)
assert outcome == [False, True, False]
```

The program prints:

```text
phenol as a molecule query, in aspirin: True
c[OX2H] in aspirin: False
[#6;a;H0]:1(:[#6;a;H1]:[#6;a;H1]:[#6;a;H1]:[#6;a;H1]:[#6;a;H1]:1)-[#8;A;H1]
aspirin, phenol, 4-tert-butylphenol: [False, True, False]
```

## Time limits, node limits and cancellation

Each search stops at the first of these limits:

- `timeout_ms`: the time limit per target in milliseconds. The default `-1`
  selects 10,000 ms for a substructure search; `0` stops at once.
- `node_limit`: the number of search nodes. The default `0` sets no limit;
  reaching one gives `LIMIT_REACHED`.
- `cancellation`: a `CancellationToken`. Calling `cancel()` from any thread
  stops the search with `CANCELLED`.

A stopped search with no matches proves nothing: `optimal` and `complete` are
both false. A token stays cancelled and is shared by every options object and
prepared query that holds it, so use a new token for each request.

```python
import smsd

query = smsd.Substructure.find("c1ccccc1O")   # phenol
target = "CC(C)(C)c1ccc(O)cc1"                # 4-tert-butylphenol


def report(label, result):
    print(f"{label:14}", result.status.name, "found:", result.found,
          "optimal:", result.optimal, "complete:", result.complete)


stopped = query.with_timeout(0).match(target)
report("timeout 0", stopped)
assert stopped.status == smsd.SearchStatus.TIMED_OUT and not stopped.optimal

options = smsd.Substructure.Options()
options.node_limit = 3
capped = query.with_options(options).match(target)
report("node limit 3", capped)
assert capped.status == smsd.SearchStatus.LIMIT_REACHED and not capped.complete

token = smsd.CancellationToken()
options = smsd.Substructure.Options()
options.cancellation = token
options.timeout_ms = 1000
cancellable = query.with_options(options)
before = cancellable.match(target)
report("before cancel", before)
assert before.found and before.status == smsd.SearchStatus.COMPLETED
token.cancel()
cancelled = cancellable.match(target)
report("after cancel", cancelled)
assert cancelled.status == smsd.SearchStatus.CANCELLED and token.is_cancelled
```

The program prints:

```text
timeout 0      TIMED_OUT found: False optimal: False complete: False
node limit 3   LIMIT_REACHED found: False optimal: False complete: False
before cancel  COMPLETED found: True optimal: True complete: True
after cancel   CANCELLED found: False optimal: False complete: False
```

To stop a long search from another thread, raise the token there. Keep a
limit and a time limit as a safety net. The number of embeddings found before
the stop depends on timing, so the program below does not print it.

```python
import threading

import smsd

# A carbon ladder of 181 atoms holds a very large number of 20-carbon chains
ladder = ("".join(f"C%{ring}" for ring in range(10, 100)) + "C"
          + "".join(f"C%{ring}" for ring in range(99, 9, -1)))
chain = smsd.Substructure.find("C" * 20)

token = smsd.CancellationToken()
options = smsd.Substructure.Options()
options.cancellation = token
options.timeout_ms = 5000                     # safety net
timer = threading.Timer(0.1, token.cancel)
timer.start()
result = chain.with_options(options).match_all(ladder, limit=1_000_000)
timer.join()

print("status:", result.status.name, "complete:", result.complete)
assert result.status == smsd.SearchStatus.CANCELLED and not result.complete
```

The program prints:

```text
status: CANCELLED complete: False
```

## Draw a hit

`depict_mapping(molecule, mapping)` draws a molecule natively, with no RDKit,
and highlights the atoms in `mapping`. Its keys are atoms of the drawn
molecule, so invert the search mapping to draw the target. This route works
for molecule and SMARTS queries alike.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

paracetamol = "CC(=O)Nc1ccc(O)cc1"
hit = smsd.Substructure.smarts("c[OX2H]").match(paracetamol)
assert hit.found

highlight = {target: query for query, target in hit.best.atom_mapping.items()}
print("highlighted target atoms:", sorted(highlight.items()))
assert highlight == {7: 0, 8: 1}

svg = smsd.depict_mapping(paracetamol, highlight)
with TemporaryDirectory() as directory:
    path = Path(directory) / "paracetamol_hit.svg"
    smsd.save_svg(svg, path)
    written = path.read_text(encoding="utf-8")
print("SVG written:", "<svg" in written)
assert "<svg" in written
```

The program prints:

```text
highlighted target atoms: [(7, 0), (8, 1)]
SVG written: True
```

With RDKit installed, `depict_substructure(query, target, *, timeout_ms=10000,
size=(600, 300))` runs its own molecule-query search and returns an RDKit
image of the query and the target with the matched atoms highlighted. It does
not accept SMARTS and reports no search status: when the search stops without
a match, it draws the pair without a highlight.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("c1ccccc1O")              # phenol
target = Chem.MolFromSmiles("CC(C)(C)c1ccc(O)cc1")   # 4-tert-butylphenol
image = smsd.depict_substructure(query, target)
print("image size:", image.size)
assert image.size == (600, 300)

with TemporaryDirectory() as directory:
    path = Path(directory) / "phenol_hit.png"
    image.save(path)
    assert path.stat().st_size > 0
```

The program prints:

```text
image size: (600, 300)
```

## R-group decomposition

`Substructure.decompose_r_groups(core, targets, options=None)` searches each
target for the core and returns a dictionary with the matched `core` graph,
the substituents as `rgroups` (`R1`, `R2` and so on) and the `search` result.
The core keeps its original hydrogen counts, so its attachment atoms are
written as bracket atoms.

```python
import smsd

targets = ["Cc1ccccc1O",   # 2-methylphenol
           "CCN"]          # ethylamine
parts = smsd.Substructure.decompose_r_groups("c1ccccc1", targets)   # benzene core
summary = []
for target, part in zip(targets, parts):
    found = part["search"].found
    core = smsd.to_smiles(part["core"]) if found else None
    groups = {name: smsd.to_smiles(graph) for name, graph in part["rgroups"].items()}
    print(f"{target:11}", "found:", found, "core:", core, "R-groups:", groups)
    summary.append((found, core, groups))
assert summary == [(True, "[c]1[c]cccc1", {"R1": "[CH3]", "R2": "[OH]"}),
                   (False, None, {})]
```

The program prints:

```text
Cc1ccccc1O  found: True core: [c]1[c]cccc1 R-groups: {'R1': '[CH3]', 'R2': '[OH]'}
CCN         found: False core: None R-groups: {}
```

## Options

| Method | Default | Meaning |
| --- | --- | --- |
| `Substructure.find(query)` | | Prepare a molecule query from SMILES, a `MolGraph` or an RDKit molecule |
| `Substructure.smarts(pattern)` | | Compile a SMARTS pattern; an empty or non-string pattern raises `ValueError` |
| `with_options(options)` | | A new query with a copy of `options` |
| `with_timeout(milliseconds)` | | A new query with `timeout_ms` set |
| `with_induced(enabled=True)` | `enabled=True` | A new query with induced matching on or off |
| `match(target)` | | One embedding, with its status |
| `match_all(target, limit=1000)` | `limit=1000` | Every embedding up to `limit`; `limit` must be a positive integer |
| `match_many(targets, *, num_threads=0)` | `num_threads=0` | One result per target, in order; `num_threads` must be a non-negative integer |
| `match_corpus(corpus, *, gpu_prefilter=False, num_threads=0)` | `gpu_prefilter=False`, `num_threads=0` | One decision per corpus slot; `corpus` must be a `TargetCorpus` |
| `Substructure.decompose_r_groups(core, targets, options=None)` | `options=None` (default options) | Core, R-groups and search result for each target |

`Substructure.Options()` holds these fields. Copy one with
`Substructure.Options(other)`.

| Field | Type | Default | Meaning |
| --- | --- | --- | --- |
| `chemistry` | `ChemOptions` | `ChemOptions()` | Matching rules (see [Chemistry options](#chemistry-options)) |
| `timeout_ms` | `int` | `-1` | Time limit per target; a negative value selects 10,000 ms; `0` stops at once |
| `node_limit` | `int` | `0` | Search node limit; `0` sets no limit |
| `cancellation` | `CancellationToken` or `None` | `None` | A token that stops the search when cancelled |
| `induced` | `bool` | `False` | Require an induced embedding |

## Notes

- `Substructure.Options` and `ChemOptions` take no keyword arguments. Set each
  field as an attribute; a misspelt attribute raises `AttributeError`.
  `options.chemistry` can be edited in place.
- `with_options` copies the options object, so later edits do not reach the
  query. A `CancellationToken` held in the options is shared, not copied.
- A query with several components, such as `C.O`, maps each component to
  distinct target atoms: `C.O` is found in ethanol, `O.O` is not.
- An empty SMILES string raises `ValueError`. An empty `MolGraph` query, from
  `parse_smiles("")`, is found in every target with one empty match.
- `Match.atom_mapping` is read-only, and `Match` and `SearchResult` are frozen.
- `elapsed_ms` is in whole milliseconds and `nodes_visited` measures effort;
  both vary between runs and builds.
- `match` and `match_all` count the time spent building Python results
  towards the time limit. `match_many` stops building results only when
  cancelled.
- A molecule query with an explicit `[H]` atom needs an explicit hydrogen atom
  in the target; a bracket hydrogen count is not an atom (see
  [Input atom numbering](../PYTHON.md#input-atom-numbering)).
- `num_threads=0` means the OpenMP default in Python and one thread per
  processor in C++. In Java, `matchMany(targets)` without a thread count runs
  on one thread.
- The substructure tools of the local molecular tools server (`smsd-mcp`)
  use their own defaults: a time limit of 1,000 ms and a node limit of
  100,000; see [Molecular tools server](COMMAND_LINE.md#molecular-tools-server).

## See also

- [SMARTS queries](SMARTS.md)
- [Maximum common substructure](MCS.md)
- [Library screening](SCREENING.md#substructure-screen-with-exact-confirmation)
- [Chemistry options in the Python guide](../PYTHON.md#chemistry-options)
- [Search guarantees](../ALGORITHMS.md#what-each-search-guarantees)
- [C++ prepared queries](../CPP.md#prepared-queries)
- [Java prepared searches](../JAVA.md#prepared-searches)

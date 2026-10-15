# SMSD Python guide

The source API uses `MCS` and `Substructure` queries. It replaces the search
functions in release 7.2.2; build the current source to run these examples.
Molecular I/O, fingerprints, layout and optional molecule conversion remain
available from `smsd`.

## Build from source

Use Python 3.9 or later, a C++17 compiler and CMake 3.18 or later:

```bash
python -m pip install -e ".[dev]" \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

CPU matching requires no external chemistry toolkit. OpenMP is optional.
Optional conversion and drawing require their own Python dependencies.

## Search and results

```python
import smsd

query = smsd.Substructure.find("CO")
assert query.match("CCO").found
assert not query.match("CCN").found
embeddings = query.match_all("OCO", limit=10)
assert len(embeddings.matches) == 2 and embeddings.complete

result = smsd.MCS.find("c1ccc(O)cc1").with_timeout(1000).match("c1ccc(N)cc1")
assert result.best.atom_count == 6 and result.optimal
print(result.best.atom_mapping)
```

`match` requests one best match. `match_all` requests all optimal atom and
selected-bond embeddings, up to a positive limit. Both return `SearchResult`:

| Field | Meaning |
|---|---|
| `matches` | Immutable tuple of `Match` objects |
| `best` / `found` | First match or `None`; whether a witness exists |
| `status` | `COMPLETED`, `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED` |
| `optimal` | The objective is proved optimal |
| `complete` | The requested search or enumeration has completed |
| `score` | Primary objective value, in input units for weights |
| `nodes_visited` / `elapsed_ms` | Search effort and elapsed milliseconds |

A timeout or node limit can return a valid incumbent without proving an
optimum. An empty incomplete result does not prove absence. An enumeration
limit can leave `optimal=True` and `complete=False`.

Each `Match` contains an injective `atom_mapping`, selected `bond_mapping`
entries `(query_first, query_second, target_first, target_second)`, and
`fragments` as query atom indices. `atom_count` and `bond_count` describe
that selected graph. `match.to_smiles(query)` writes those atoms and bonds,
capping cut valences where needed. `to_graph(query)` instead retains the
original chemical context for further matching.

## MCIS, MCES and fragments

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
```

MCIS is the default: maximise atoms, then common bonds, preserving edges and
nonedges between mapped vertices. MCES maximises selected compatible common
edges, then incident atoms. It can omit input edges and excludes isolated
atom padding. Connectivity refers to selected common edges. The default
requires one component; `with_fragments(min_atoms=1, max_fragments=...)`
permits multiple components.

## Chemistry, weights and budgets

```python
import smsd

options = smsd.MCS.Options()
options.chemistry.match_formal_charge = True
options.chemistry.use_chirality = True
options.chemistry.use_bond_stereo = True
options.chemistry.complete_rings_only = True
options.timeout_ms = 1000
search = smsd.MCS.find("c1ccccc1").with_options(options)
assert search.match("Cc1ccccc1").best.atom_count == 6

sub_options = smsd.Substructure.Options()
sub_options.induced = True
assert not smsd.Substructure.find("CCC").with_options(sub_options).match("C1CC1").found
```

`ChemOptions` controls element, charge, isotope, aromaticity, ring, mapped
stereo and tautomer predicates. Defaults match elements and strict bond order;
charge, isotope and stereo matching are opt-in. `tautomer_aware` does not
change bond-order mode. Set `match_bond_order` to `BondOrderMode.LOOSE` or
`ANY` explicitly when needed. `ChemOptions.profile(name)` supplies named
profiles; unknown names raise `ValueError`.

`complete_rings_only` applies to both inputs. A retained ring atom or bond
must belong to a fully selected relevant ring in its own input. A fusion atom
may survive through one complete ring.

`timeout_ms=-1` selects an adaptive budget; zero requests an immediate stop.
Positive values apply to one query/target call. `node_limit=0` is unlimited.
Set `options.cancellation` to a `CancellationToken` and call `cancel()` from
another thread to stop cooperatively. Preparation and individual steps can
extend elapsed time beyond the deadline.

For MCIS, `with_atom_weights(weights)` ranks finite signed query atom weights
before atom and bond tie-breaks. Negative values can favour a smaller result.
Weights use bounded integer millipoints. MCES rejects atom weights.

`MCS.compute_tautomer_confidence(query, target, match, chemistry)` reports a
heuristic confidence for selected tautomeric chemistry relaxations. Its pH
and solvent context is applied when called, without changing molecules or
the search objective. Identical selected chemistry returns 1. The value is
not a probability.

Prepared queries own their query and options. Fluent configuration returns
a new query; separate calls have independent search state and may run
concurrently. Do not mutate a target while a call uses it.

## SMARTS and target collections

```python
import smsd

query = smsd.Substructure.smarts("[#6]-[#8]")
assert [result.found for result in query.match_many(["CCO", "CCN"])] == [True, False]
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN"])
corpus.prewarm(num_threads=1)
assert [result.found for result in query.match_many(corpus, num_threads=2)] == [True, False]
```

SMARTS compiles once and searches complete pattern embeddings using the same
result, budget and enumeration API. `match_many` returns one `SearchResult`
per target in input order, with an independent per-pair budget.
`num_threads=0` uses the configured OpenMP worker count; without OpenMP it
runs sequentially. `TargetCorpus` captures molecules and their input atom
numbering when added, and provides reusable fingerprints and graph invariants.

## Input atom numbering

SMILES strings use parser order. `MolGraph` inputs use graph order. Optional
external molecule inputs use their original order, independently on each side.
Weights use the query input's order. Explicit hydrogen vertices remain vertices;
bracket hydrogen counts do not add vertices.

```python
from rdkit import Chem
import smsd

query = Chem.RenumberAtoms(Chem.MolFromSmiles("NCCO"), [3, 1, 0, 2])
target = Chem.MolFromSmiles("CC(O)CN")
result = smsd.MCS.find(query).match(target)
assert result.best.atom_count == 4
for first, second in result.best.atom_mapping.items():
    assert query.GetAtomWithIdx(first).GetAtomicNum() == target.GetAtomWithIdx(second).GetAtomicNum()
```

Conversions are reused while the external molecule remains unchanged.
`from_rdkit(molecule, use_cache=False)` bypasses the cache. `get_index_map(graph)`
returns a copy of the correspondence; `clear_cache()` retains metadata for
live graphs. Search results retain their own correspondence.

## Utilities

`MCS.validate_mapping(query, target, match, options)` checks the selected graph
against an objective and chemistry policy. Pass a `Match` for MCES so omitted
edges remain explicit. `MCS.is_mapping_maximal` checks extension; maximality
is weaker than an optimality proof. `MCS.canonicalize_mapping` chooses an atom
mapping representative under symmetry and raises when the calculation cannot
complete.

Other APIs include fingerprints and metrics (`fingerprint`, `circular_fingerprint`,
`tanimoto_coefficient`), MOL/SDF I/O, CIP assignment, graph extraction and SVG
drawing. Canonical hashes are screening keys, not proofs of identity.
See [examples](EXAMPLES.md), [C++](CPP.md) and [Java](JAVA.md).

# SMSD Pro Python guide

This guide is the reference for the `smsd` Python package, version 8.0.0. It
describes the search model, every option and its default, and each public name
that `smsd` exports. [Python examples](EXAMPLES.md) gives runnable code for
fingerprints, drawing, layout, stereo, MOL and SDF files, R-groups, scaffolds
and graph utilities.

## Installation

```bash
python -m pip install "smsd==8.0.0"
```

Wheels are built for CPython 3.14; the [installation guide](HOWTO-INSTALL.md)
lists the supported platforms. `smsd.__version__` reports the installed version.

Core matching needs no other chemistry toolkit. RDKit is optional and is needed
only to pass RDKit molecules or to use `to_rdkit`, `to_svg`, `depict_mcs` and
`depict_substructure`; install it separately or with the `rdkit` extra,
`python -m pip install "smsd[rdkit]==8.0.0"`. Open Babel is optional and is needed only
for `from_openbabel` and `to_openbabel`.

### Build from source

Use Python 3.9 or later, a C++17 compiler and CMake 3.18 or later. From the
repository root:

```bash
python -m pip install -e ".[dev]" \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

These settings match the published CPU/OpenMP wheels. Without them, CMake builds
a Metal or CUDA backend when it finds one (`AUTO`). A GPU backend can compute
the initial atom-compatibility filter; search and proof always run on the CPU.
Setting the environment variable `SMSD_DISABLE_GPU` to a value other than an
empty string or `0` disables the GPU at run time. OpenMP is optional; without
it, batch calls run sequentially.

## Searches

`Substructure` finds embeddings of a whole query in a target. `MCS` finds a
maximum common substructure (MCS) of a query and a target. Each prepares its
query once for repeated use:

| Prepare | Query |
|---|---|
| `Substructure.find(query)` | A molecule; SMILES strings are parsed as molecules |
| `Substructure.smarts(query)` | A SMARTS pattern, compiled once; an empty or non-string pattern raises `ValueError` |
| `MCS.find(query)` | A molecule |

Queries and targets may be SMILES strings, `MolGraph` objects or RDKit
molecules. A prepared query provides:

| Method | Returns |
|---|---|
| `match(target)` | A `SearchResult` with at most one match: one optimum for `MCS`, one embedding for `Substructure` |
| `match_all(target, limit=1000)` | A `SearchResult` with up to `limit` matches: all optimal `MCS` mappings, or all `Substructure` embeddings |
| `match_many(targets, *, num_threads=0)` | A list with one `SearchResult` per target, in input order |
| `with_options(options)`, `with_timeout(milliseconds)` | A new prepared query with a copy of the options |
| `MCS.with_objective(objective)`, `MCS.with_fragments(*, min_atoms=1, max_fragments=2147483647)`, `MCS.with_atom_weights(weights)` | A new `MCS` query with that setting |
| `Substructure.with_induced(enabled=True)` | A new `Substructure` query with induced matching |

`limit` must be a positive integer and `num_threads` a non-negative integer;
other values raise `ValueError`. Symmetry-related mappings are reported as
separate matches.

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

### Results

Every search returns a `SearchResult`:

| Field | Meaning |
|---|---|
| `matches` | Tuple of `Match` objects |
| `best`, `found` | The first match or `None`; whether any match exists |
| `status` | `SearchStatus.COMPLETED`, `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED` |
| `optimal` | The MCS optimum, or the substructure existence or absence answer, is proved |
| `complete` | The requested search or enumeration finished |
| `score` | Primary objective value; with atom weights, in the supplied weight units |
| `nodes_visited`, `elapsed_ms` | Search effort and elapsed milliseconds |

A timeout or node limit can return a valid incumbent without proving an
optimum, and `match_all` may then hold only the best matches found so far. An
empty incomplete result does not prove absence. Check `status`, `optimal` and
`complete` before treating a result as final. An enumeration limit can leave
`optimal=True` and `complete=False`; reaching exactly the limit does not imply
truncation when no further matches exist. The order of tied matches is not a
chemical ranking.

```python
import smsd

result = smsd.Substructure.find("CC").match_all("CCC", limit=1)
assert len(result.matches) == 1 and result.optimal and not result.complete
assert result.status == smsd.SearchStatus.LIMIT_REACHED
```

Each `Match` holds an injective `atom_mapping` (a read-only mapping from query
atom index to target atom index), `bond_mapping` entries
`(query_first, query_second, target_first, target_second)` for the selected
bonds, and `fragments`, the connected components of the selected bonds as
query atom indices. `atom_count` and `bond_count` describe that selected graph.
`match.to_smiles(query)` writes the selected atoms and bonds as SMILES, capping
cut valences where needed; `match.to_graph(query)` returns them as a `MolGraph`
that keeps the original chemical perception and hydrogen counts for further
matching.

### Concurrency

Prepared queries own their query and options, and every configuration method
returns a new query. Calls keep independent search state and release the
Python global interpreter lock while they search, so one prepared query can
serve several threads. Do not modify a target while a call is using it.

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
all_matches = smsd.MCS.find(query).match_all(target, limit=100)
assert all_matches.optimal and all_matches.complete
```

`MCS` defaults to a connected maximum common induced subgraph (MCIS): it
maximises mapped atoms, then common bonds, and preserves both edges and
non-edges between mapped atoms. The maximum common edge subgraph (MCES)
objective maximises selected compatible common bonds, then their incident
atoms; it may omit input bonds and never adds isolated atoms.

Connectivity refers to the selected common bonds. The default requires one
connected component. `with_fragments(min_atoms=1, max_fragments=2147483647)`
permits several components, each with at least `min_atoms` atoms; with a
minimum of one, MCIS may include isolated mapped atoms. Fragment limits below
one raise `ValueError`.

`Substructure.find` is non-induced by default. `with_induced()`, or
`induced = True` on `Substructure.Options`, also preserves non-edges.

## Chemistry options

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

`ChemOptions`, available as `options.chemistry`, decides which atoms and bonds
are compatible:

| Attribute | Default | Effect |
|---|---|---|
| `match_atom_type` | `True` | Matched atoms must be the same element |
| `match_bond_order` | `BondOrderMode.STRICT` | `STRICT` requires equal bond orders, subject to `aromaticity_mode`; `LOOSE` and `ANY` relax the comparison |
| `aromaticity_mode` | `AromaticityMode.FLEXIBLE` | `FLEXIBLE` lets aromatic bonds match single or double bonds; `STRICT` requires atoms and bonds to agree on aromaticity |
| `aromaticity_model` | `AromaticityModel.DAYLIGHT_LIKE` | Model used where aromaticity is perceived again |
| `match_formal_charge` | `False` | Matched atoms must carry the same formal charge |
| `match_isotope` | `False` | Mass numbers must agree where both atoms specify one |
| `use_chirality` | `False` | Check tetrahedral stereo |
| `use_bond_stereo` | `False` | Check double-bond stereo |
| `ring_matches_ring_only` | `False` | Ring atoms and bonds match only ring atoms and bonds |
| `complete_rings_only` | `False` | Retain only complete rings |
| `ring_fusion_mode` | `RingFusionMode.IGNORE` | `STRICT` requires matched ring atoms to belong to the same number of relevant cycles found by the ring search; `PERMISSIVE` currently behaves as `IGNORE` |
| `tautomer_aware` | `False` | Allow recognised tautomeric relaxations |
| `pH` | `7.4` | Context for `MCS.compute_tautomer_confidence` only |
| `solvent` | `Solvent.AQUEOUS` | Context for `MCS.compute_tautomer_confidence` only |

Stereo is checked where the selected mapping determines the configuration;
omitting ligands can leave a stereocentre of a partial MCS undetermined. A
`Substructure` search embeds the query, so with `use_chirality` a tagged query
centre needs a tagged target centre of the same winding over the mapped
ligands, and with `use_bond_stereo` a configured query double bond needs a
configured target bond of the same relation, as a SMARTS chirality primitive or
pair of directional bonds does; an untagged or unconfigured query atom or bond
matches any target. An `MCS` search keeps an unspecified target atom or bond as
a wildcard. For MCS, `complete_rings_only` applies to both inputs: a retained
ring atom or bond must belong to a fully selected relevant ring in its own
input, and a fusion atom may be retained through one complete ring.

`tautomer_aware` enables recognised tautomeric relaxations without changing
the configured bond-order mode; set `match_bond_order` to `LOOSE` or `ANY` when
that broader policy is needed:

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

`ChemOptions.profile(name="default")` returns a named preset; names are not
case-sensitive and an unknown name raises `ValueError`:

| Preset | Differs from the defaults |
|---|---|
| `default` | No differences |
| `strict` | `aromaticity_mode=STRICT` |
| `pharma` | `match_formal_charge`, `ring_matches_ring_only` and `complete_rings_only` are `True` |
| `reaction`, `compat-fmcs`, `compat-substruct` | `match_bond_order=LOOSE` |
| `ChemOptions.tautomer_profile()` | `tautomer_aware=True` and `match_bond_order=LOOSE` |

No preset enables the isotope or stereo checks, and only `pharma` enables the
charge check. `ChemOptions(options)`, `copy.copy` and `copy.deepcopy` give an
independent copy. The setters `with_bond_order_mode`, `with_aromaticity_mode`,
`with_complete_rings_only`, `with_match_isotope`, `with_ring_fusion_mode`,
`with_tautomer_aware`, `with_solvent` and `with_ph` each set a field and return
the same options, as in C++ and Java.
`ChemOptions.aromaticity_model` is deprecated: no search reads it,
and `MolGraphBuilder.build()` and `perceive_aromaticity()` take the model
directly.

## Timeouts, node limits and cancellation

The search budget of a call is its timeout together with its node limit.

- `timeout_ms` defaults to `-1`. Any negative value selects the default
  timeout: 10 seconds for `Substructure`, and for `MCS`
  500 + 2 × (query atoms × target atoms) milliseconds, capped at 30 seconds.
  Zero requests an immediate stop. A positive value applies to each query and
  target pair.
- `node_limit` defaults to `0`, which sets no node limit; the timeout still
  applies. A node stop reports `SearchStatus.LIMIT_REACHED`.
- `cancellation` defaults to `None`. Set it to a `CancellationToken` and call
  `cancel()` from another thread to stop the search cooperatively; the result
  reports `SearchStatus.CANCELLED`. A cancelled token stays cancelled, so use a
  new token for a later search.

Parsing and other preparation can extend the total elapsed time beyond the
search deadline. For `match` and `match_all`, converting matches into Python
objects counts towards the timeout, as native result building does: once the
timeout and a further eighth of it have elapsed, or the token is cancelled, the
matches converted so far are returned with `TIMED_OUT` or `CANCELLED` and
`complete=False`. In `match_many`, only cancellation stops that conversion.

## Atom weights

For MCIS, `with_atom_weights(weights)` ranks results by the sum of finite,
signed query atom weights before the atom and bond tie-breaks; negative weights
can favour a smaller result. Supply one weight for every query atom, in the
query's input order. The search scores weights in whole thousandths, truncated
towards zero, so differences below 0.001 can tie; `result.score` reports the
weighted score in the supplied units. Non-finite weights, weight totals beyond
the supported range of about ±2.1 million, a list that does not cover every
query atom, and weights combined with MCES all raise `ValueError`.

## Tautomer confidence

`MCS.compute_tautomer_confidence(query, target, match, chemistry)` returns a
heuristic confidence for the tautomeric relaxations that a `Match` uses, read
in the `pH` and `solvent` context of `chemistry`. It changes neither the
molecules nor the search objective. Identical selected chemistry returns 1.
The value is not a probability. `match` must be a `Match`; a plain atom mapping
raises `TypeError`. `ChemOptions.with_ph(pH)` sets a finite pH, clamped to
[0, 14], and `ChemOptions.with_solvent(solvent)` sets the solvent.

The confidence is the geometric mean, over the relaxed atoms with a tautomer
class in either molecule, of the smaller of the two atom weights. Each atom of
a class takes the weight of the first motif in the table below whose pattern
covers the atom, and the ionisable site and aqueous pKa that the matches of
that motif covering the atom carry. Each motif takes every atom of all its
matches before the next is tried, so the motif chosen for an atom does not
depend on the order of the atoms or of their neighbours. The atoms of one match
share a class, and matches that share an atom join their classes, so atoms that
interconvert through any chain of motifs form one class, in the manner of the
mobile-hydrogen groups of InChI but without requiring a mobile hydrogen, and
the classes hold exactly the atoms that some motif covers. The partition
depends on the graph alone: `ONC(C)=O` and `CC(=O)NO` both hold one class of
all five atoms, succinimide one class of all seven, and 2-hydroxypyridine,
written `Oc1ccccn1` or `n1ccccc1O`, one class of the oxygen, the ring nitrogen,
the carbon between them and its para carbon. `MolGraph.tautomer_class` numbers
the classes from 0 in increasing order of a key, the sorted colours of their
atoms, refined in the manner of Weisfeiler and Leman from element, charge, mass
number, hydrogen count, aromaticity, ring membership, degree and motif cover;
of classes that still share a key, which the refinement cannot tell apart (in
practice those a symmetry of the graph exchanges), the one with the lowest atom
index takes the next number and the refinement resumes. The numbers are
therefore canonical up to Weisfeiler–Leman refinement, the same for every atom
order except among classes that the refinement cannot tell apart: the two
carbonyl classes of `CC(=O)CCC(C)=O` are numbered 0 and 1, the class of the
lower atoms first. An atom outside every class has class −1, weight 1 and no
site. The weight is the relevance of the dominant neutral tautomer of the motif
and does not depend on pH.

| Order | Motif | Weight | Site | pKa |
| --- | --- | --- | --- | --- |
| 1 | Hydroxamic acid | 0.85 | acid | 9 |
| 2 | Imide | 0.92 | acid | 9.5 |
| 3 | Barbituric: a ring N between two ring carbonyls, whose N and carbons the imide motif takes first, so it takes the carbonyl oxygens | 0.80 | — | — |
| 4 | Urea | 0.96 | acid | 14 |
| 5 | 1,3-Diketone | 0.72 | acid | 9.5 |
| 6 | β-Enaminone | 0.85 | — | — |
| 7 | 1,5-Keto-enol | 0.75 | — | — |
| 8 | Pyridone/hydroxypyridine | 0.95 | acid | 11.7, or 8.6 for a pyrimidin-4(3H)-one |
| 9 | Pyridinone | 0.95 | acid | 11.7, or 8.6 for a pyrimidin-4(3H)-one |
| 10 | Hydroxypyrimidine, the pattern of row 8, so it takes no atom | 0.90 | — | — |
| 11 | Amide/imidic acid | 0.97 | acid | 16 |
| 12 | Keto/enol | 0.98 | acid | 20 |
| 13 | Thioamide | 0.94 | acid | 13 |
| 14 | Thione/thiol | 0.94 | acid | 10.5 |
| 15 | Thioamide/iminothiol | 0.95 | acid | 13 |
| 16 | Nitroso/oxime | 0.95 | — | — |
| 17 | Phenol/quinone methide | 0.92 | — | — |
| 18 | Phenol/quinone | 0.88 | acid | 10 |
| 19 | Tetrazole, 1,2,3-triazole, benzimidazole, pyrazole and imidazole N–H shifts | 0.50 | — | — |
| 20 | Guanidine, amidine | 0.50 | — | — |
| 21 | Phosphonate | 0.50 | — | — |
| 22 | Nitro/aci-nitro | 0.95 | — | — |
| 23 | Furanose/pyranose ring–chain | 0.60 | — | — |
| 24 | Lactol ring–chain, the pattern of row 23, so it takes no atom | 0.65 | — | — |
| 25 | Allyl shift | 0.70 | — | — |
| 26 | Selenoketone/selenol | 0.94 | — | — |

The confidence is computed by the native core, as in C++ and Java. The solvent
shift is applied first; the pH then scales the weight of every atom with a site
by the neutral fraction of its class. With pKa_a the lowest acid pKa and pKa_b
the highest base pKa (that of BH+) among the sites of the class:

```text
f0 = 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH))
weight = min(0.99, max(0.1, solvent-corrected weight × f0))
```

A class with acid sites only drops the base term, which leaves
1 / (1 + 10^(pH − pKa_a)), the fraction of the neutral acid HA; a class with
base sites only drops the acid term, which leaves 1 / (1 + 10^(pKa_b − pH)),
the fraction of the neutral base B. The library assigns no pKa to the motifs
marked — in the table, and the pH does not change their weights.

A match of a motif with a site carries it only when the atom of the match that
would release the proton bears a hydrogen: the N–O oxygen of a hydroxamic acid,
so an N-alkoxy amide such as a Weinreb amide has none; the nitrogen of an
imide, urea, amide, thioamide or thioamide/iminothiol; the alpha carbon of a
ketone or thione and the bridging carbon of a 1,3-diketone; the ring nitrogen
or the oxygen of a pyridone/hydroxypyridine or pyridinone; and a phenol oxygen
with no second heavy-atom neighbour, so anisoles, aryl esters and lactam or
lactone carbonyls have none. When that atom belongs to an earlier motif, the
match carries the site of that motif instead, and an atom that several matches
of its motif cover takes the lowest pKa among them. The
pyridone/hydroxypyridine and pyridinone motifs take the N–H pKa of 2-pyridone,
11.7, or 8.6 when the other ring neighbour of the nitrogen is a carbon bonded
to a second aromatic ring nitrogen, as in pyrimidin-4(3H)-one. No motif in the
table has a base site. The pKa values are aqueous and are applied in every
solvent, and a class with several sites uses its lowest acid pKa and its
highest base pKa.

For example, against its iminol `ON=C(C)O` the hydroxamic acid `ONC(C)=O`
scores 0.85 at pH 2, 0.829 at pH 7.4, 0.773 at pH 8, 0.425 at pH 9 and the 0.1
floor at pH 10 and 12, and the same with its atoms written `CC(=O)NO`, while
the N-alkoxy amide `CONC(C)=O`, which has no O–H, scores 0.85 at every pH
against `CON=C(C)O`. The 2-pyridone `[nH]1ccccc1=O`, whose class also holds the
carbonyl and para carbons that the phenol motif reaches without an O–H, takes
the N–H pKa of 11.7 and scores 0.914 at pH 2 and 7.4, 0.905 at pH 10 and 0.317
at pH 12 against `n1ccccc1O`; `O=c1cccc[nH]1` scores the same against
`Oc1ccccn1`, whose para carbon takes that pKa through the O–H of the
pyridone/hydroxypyridine motif. 4(3H)-Pyrimidinone `O=c1ccnc[nH]1`, of pKa 8.6,
scores 0.95 at pH 2, 0.894 at pH 7.4, 0.271 at pH 9 and the 0.1 floor at pH 10
and 12 against `Oc1ccncn1`. The amidine `CC(=N)N` against `CC(N)=N` scores 0.5
at every pH. The first example in Python:

```python
import smsd

query, target = smsd.parse_smiles("ONC(C)=O"), smsd.parse_smiles("ON=C(C)O")
options = smsd.MCS.Options()
options.chemistry = smsd.ChemOptions.tautomer_profile()
match = smsd.MCS.find(query).with_options(options).match(target).best
for pH in (2.0, 7.4, 9.0, 12.0):
    chemistry = smsd.ChemOptions.tautomer_profile().with_ph(pH)
    print(pH, round(smsd.MCS.compute_tautomer_confidence(query, target, match, chemistry), 3))
```

prints the confidences 0.85, 0.829, 0.425 and 0.1 after the four pH values.

`MCS.validate_tautomer_consistency(query, target, match)` reports whether a
mapping conserves the mobile N–H, O–H and S–H protons within each mapped
tautomeric group; `match` may be a `Match` or a plain atom mapping.

## SMARTS and target collections

```python
import smsd

query = smsd.Substructure.smarts("[#6]-[#8]")
assert [result.found for result in query.match_many(["CCO", "CCN"])] == [True, False]
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN"])
corpus.prewarm(num_threads=1)
assert [result.found for result in query.match_many(corpus, num_threads=2)] == [True, False]
```

A SMARTS query is compiled once and finds complete pattern embeddings with the
same result, timeout and enumeration rules as a molecule query. `match_many`
gives each target its own search budget. `num_threads=0` uses the OpenMP
default thread count, which `OMP_NUM_THREADS` can set; without OpenMP, the call
runs sequentially.

`TargetCorpus` captures molecules and their input atom numbering when they are
added, so later changes to the inputs do not affect it. `prewarm()` computes
graph invariants and ECFP4 fingerprints (2,048 bits) once. `screen(query,
threshold)` returns the indices of targets whose ECFP4 Tanimoto similarity to
the query is at least `threshold`; it returns an empty list until `prewarm()`
has run after the last addition. A screening threshold does not prove that an
MCS or substructure match is absent.

## Constrained and consensus searches

```python
import smsd

results = smsd.MCS.match_constrained(["CO", "CN"], "CO.CN")
assert [result.best.atom_count for result in results] == [2, 2]
assert set(results[0].best.atom_mapping.values()).isdisjoint(
    results[1].best.atom_mapping.values()
)

consensus = smsd.MCS.consensus(["CCO", "CCN", "O"], threshold=2 / 3)
assert consensus.found and consensus.support_count == 2
assert consensus.graph.n == 2 and not consensus.optimal
assert consensus.to_smiles() == "CC"
```

`MCS.match_constrained(queries, target, options=None)` searches the queries in
order and excludes target atoms claimed by an earlier result. Each pair has its
own budget and `SearchResult`. The sequential assignment does not prove a
globally optimal joint mapping. `MCS.Options.excluded_target_atoms` reserves
further original target indices without deleting atoms or changing their
chemical context.

`MCS.consensus(molecules, options=None, threshold=1.0)` progressively reduces
candidate graphs and verifies support in at least
`ceil(threshold * len(molecules))` inputs; `threshold` must lie in (0, 1]. One
budget covers the whole call. The result reports `optimal=False` even when
`complete=True`, because completing this procedure does not prove a globally
optimal common graph. Atom weights and target exclusions raise `ValueError`.
`ConsensusResult.mappings` maps each supporting input index to a `Match` from
the atoms of `consensus.graph`, numbered from 0; `fragment` and `source_index`
identify the selected graph in its source input.

## Validation and symmetry

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

`MCS.validate_mapping(query, target, match, options=None)` checks a selected
graph against an objective and chemistry policy and returns a list of errors;
an empty list means the mapping is valid. Pass a `Match` for MCES so that
omitted bonds remain explicit; a plain atom mapping selects every compatible
mapped bond. `MCS.is_mapping_maximal(query, target, match, options=None)`
checks whether a valid common graph admits a compatible extension; maximality
is weaker than an optimality proof, and the call raises an exception if its
extension search stops before a proof. `MCS.canonicalize_mapping(query, target,
atom_mapping)` returns the lexicographically smallest atom mapping equivalent
under the symmetries of both molecules, and raises an exception when that
calculation cannot complete.

## R-group decomposition

`Substructure.decompose_r_groups(core, targets, options=None)` returns one
dictionary per target with `core` (a `MolGraph`), `rgroups` (labels such as
`"R1"` mapped to graphs) and `search` (the `SearchResult` for the core). Check
the search result before using a decomposition; see the
[R-group example](EXAMPLES.md#r-group-decomposition).

## Input atom numbering

SMILES strings use parser order and `MolGraph` inputs use graph order. RDKit
molecule inputs keep their original order, independently on each side, and
search results report RDKit atom indices. Atom weights use the query input's
order. Explicit hydrogen atoms are vertices; bracket hydrogen counts do not add
vertices.

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

Conversion passes through SMILES, so coordinates, data fields and query
predicates are not transferred. Conversions are reused while the RDKit molecule
is unchanged; `from_rdkit(molecule, use_cache=False)` bypasses the cache.
`get_index_map(graph)` returns a copy of the atom correspondence, and
`clear_cache()` keeps that correspondence for live graphs. Search results hold
their own correspondence.

`from_openbabel` accepts Pybel molecules and `OBMol` objects and
`to_openbabel` returns a Pybel molecule, both through SMILES. This bridge does
not retain the original Open Babel atom numbering. For conversion in C++, see
the [optional native adapter](CPP.md#optional-native-adapter).

## API reference

The tables list every public name in `smsd.__all__`. In the descriptions, a
*molecule* is a SMILES string, `MolGraph` or RDKit molecule; a *graph* is a
`MolGraph` only. Parse SMILES with `parse_smiles` before calling a function
that takes a graph.

### Search types and enumerations

| Name | Description |
|---|---|
| `MCS`, `Substructure` | Prepared queries; see [Searches](#searches) |
| `MCS.Options`, `Substructure.Options` | Search options; `MCS.Options(other)` and `Substructure.Options(other)` copy existing options |
| `MCS.Objective` | `MCIS` (default) or `MCES` |
| `MCS.Connectivity` | `CONNECTED` (default) or `FRAGMENTS` |
| `SearchResult`, `Match` | Search results; see [Results](#results) |
| `ConsensusResult` | `fragment`, `source_index`, `graph`, `mappings`, `support_count`, `found`, `status`, `optimal`, `complete`, `nodes_visited`, `elapsed_ms` and `to_smiles()` |
| `SearchStatus` | `COMPLETED`, `TIMED_OUT`, `CANCELLED`, `LIMIT_REACHED` |
| `CancellationToken` | `cancel()` and the `is_cancelled` property |
| `ChemOptions` | Chemistry options; see [Chemistry options](#chemistry-options) |
| `BondOrderMode` | `STRICT`, `LOOSE`, `ANY` |
| `AromaticityMode` | `STRICT`, `FLEXIBLE` |
| `AromaticityModel` | `DAYLIGHT_LIKE` |
| `RingFusionMode` | `IGNORE`, `PERMISSIVE`, `STRICT` |
| `Solvent` | `AQUEOUS`, `DMSO`, `METHANOL`, `CHLOROFORM`, `ACETONITRILE`, `DIETHYL_ETHER` |

Search option fields and defaults:

| Field | Default | Applies to |
|---|---|---|
| `chemistry` | `ChemOptions()` | Both |
| `timeout_ms` | `-1`, the default timeout | Both |
| `node_limit` | `0`, no node limit | Both |
| `cancellation` | `None` | Both |
| `induced` | `False` | `Substructure` |
| `objective` | `MCS.Objective.MCIS` | `MCS` |
| `connectivity` | `MCS.Connectivity.CONNECTED` | `MCS` |
| `min_fragment_atoms` | `1` | `MCS` |
| `max_fragments` | `2147483647` | `MCS` |
| `atom_weights` | `[]` | `MCS` with MCIS |
| `excluded_target_atoms` | `[]` | `MCS` |

### Molecules and SMILES

| Name | Description |
|---|---|
| `MolGraph` | Molecular graph used throughout the package; members are listed below |
| `MolGraphBuilder` | Builds a `MolGraph` from per-atom lists and adjacency lists |
| `parse_smiles(smiles, *, lenient=False)` | Parse SMILES into a graph; an empty string gives an empty graph, and `lenient=True` recovers from malformed text that would otherwise raise `ValueError`. A bracket atom outside the OpenSMILES limits raises `ValueError`: a hydrogen count has one digit, a charge two digits or repeated signs within -15 to +15, and an isotope three digits |
| `parse_mapped_smiles(smiles)` | Remove atom-map numbers, then parse |
| `strip_atom_maps(smiles)` | Remove atom-map numbers from SMILES text, preserving ring closures and aromatic bonds |
| `to_smiles(mol)` | Canonical SMILES of a graph, SMILES string or RDKit molecule; canonical order may differ between languages. `canonical_smiles(mol)` is a deprecated alias |
| `to_smarts(mol, include_aromaticity=True, include_charge=True, include_ring_member=False, include_h_count=False, include_isotope=False)` | SMARTS for a graph with the selected predicates |
| `canonical_hash(mol)` | Integer hash of a graph covering connectivity, bond order, aromaticity, ring membership, charge, mass number and stereo; a screening key, not a proof of identity |
| `prewarm(mol)` | Compute and cache common invariants of a molecule; returns the graph |
| `prewarm_graph(mol)` | Deprecated alias of `prewarm` that returns `None` |
| `extract_subgraph(mol, atom_indices)` | Induced subgraph of a graph, numbered from 0 in the supplied order |
| `murcko_scaffold(mol)` | Bemis-Murcko scaffold of a graph: its ring systems, the shortest paths between them and every atom joined to either by a double or triple bond, so that exocyclic C=O, C=N and C=S groups stay; each stripped substituent is replaced by hydrogen, a stereocentre keeps its configuration when all its neighbours stay and a double bond when each end keeps a ligand, and an acyclic graph is returned unchanged |
| `count_components(mol)`, `split_components(mol)` | Number of connected components of a graph, or one graph per component |
| `same_canonical_graph(mol1, mol2)` | Whether two graphs have the same atom count, canonical hash and canonical SMILES |
| `translate_to_atom_ids(mapping, g1, g2)` | Translate a mapping into external atom IDs, such as those set with `MolGraphBuilder.atom_ids`; returned unchanged when neither graph has IDs |

`MolGraph` members:

| Member | Description |
|---|---|
| `n`, `len(graph)` | Atom count |
| `atomic_num`, `formal_charge`, `mass_number`, `hydrogen_count`, `atom_class` | Per-atom values |
| `aromatic`, `ring`, `degree`, `neighbors` | Per-atom flags, degrees and adjacency lists |
| `bonds()` | All bonds as `(i, j, order)` tuples |
| `bond_order(i, j)`, `has_bond(i, j)`, `bond_in_ring(i, j)`, `bond_aromatic(i, j)` | Bond queries; an index outside the graph raises `IndexError` |
| `name`, `program_line`, `comment`, `properties` | Writable MOL header lines and SDF data fields |
| `num_rings()`, `ensure_ring_counts()`, `ring_system_of(atom)`, `ring_system_count()`, `ring_system_sig(rs_id)` | Ring counts and ring systems |
| `get_canonical_hash()`, `ensure_canonical()`, `morgan_rank`, `canonical_label`, `orbit` | Canonical ranks, labels and symmetry classes; `canonical_label` and `orbit` compute the labelling on first access, as the C++ getters do, and `morgan_rank` is filled by the same step |
| `get_automorphism_generators()`, `automorphism_generators_truncated()` | Automorphism generators, and whether they were capped |
| `label`, `tautomer_class` | Per-atom matching labels and tautomer classes |
| `prewarm()` | Compute and cache common invariants |
| `perceive_aromaticity(model=AromaticityModel.DAYLIGHT_LIKE)`, `kekulize()`, `dearomatize()` | Change the graph in place |

`MolGraphBuilder` setters return the builder, so calls can be chained:
`atom_count(n)`, `atomic_numbers(nums)`, `formal_charges(charges)`,
`mass_numbers(masses)`, `hydrogen_counts(counts)`, `ring_flags(flags)`,
`aromatic_flags(flags)`, `neighbors(adj)`, `bond_orders(orders)` (parallel to
the adjacency lists), `bond_ring_flags(flags)`, `bond_aromatic_flags(flags)`,
`atom_ids(ids)`, `name(value)`, `program_line(value)`, `comment(value)` and
`properties(values)`. `build(perceive_aromaticity=True,
aromaticity_model=AromaticityModel.DAYLIGHT_LIKE)` returns the graph.

### Toolkit conversion

| Name | Description |
|---|---|
| `from_rdkit(mol, use_cache=True)` | RDKit molecule to graph through SMILES, recording the atom correspondence; empty and hydrogen-only molecules raise `ValueError` |
| `to_rdkit(mol_or_smiles)` | Graph or SMILES string to an RDKit molecule through canonical SMILES; atom order follows that SMILES |
| `get_index_map(g)` | Copy of the correspondence `result[smsd_index] = rdkit_index` for a graph from `from_rdkit`, otherwise `None` |
| `translate_mapping(mapping, g_query=None, g_target=None)` | Translate a mapping between graphs from `from_rdkit` into RDKit indices; an index outside its molecule raises `IndexError` |
| `clear_cache()` | Clear reusable conversions; live graphs keep their correspondence |
| `from_openbabel(mol)` | Pybel molecule or `OBMol` to graph through SMILES |
| `to_openbabel(mol)` | Molecule to Pybel molecule through SMILES |

Search calls translate RDKit numbering themselves; `translate_mapping` is needed
only when you search converted graphs directly.

### Stereo, rings and atom chemistry

| Name | Description |
|---|---|
| `assign_rs(mol)` | `{atom: "R" or "S"}` for tetrahedral centres with resolvable priorities |
| `assign_ez(mol)` | `{(atom, atom): "E" or "Z"}` for stereogenic double bonds with resolvable priorities |
| `assign_cip(mol)` | `CIPDescriptors` with `rs_labels`, one `RSLabel` per atom, and `ez_bonds`, a list of `(atom, atom, EZLabel)` tuples |
| `RSLabel`, `EZLabel`, `CIPDescriptors` | `RSLabel.NONE`, `R`, `S` and the pseudoasymmetric `r`, `s`; `EZLabel.NONE`, `E`, `Z`; the descriptor set |
| `compute_sssr(mol)` | Smallest set of smallest rings, as atom-index lists in ascending size |
| `perceive_aromaticity(mol, *, model=AromaticityModel.DAYLIGHT_LIKE)` | Recompute ring membership and aromaticity |
| `kekulize(mol)` | Replace aromatic bonds with an explicit Kekulé assignment; raises `ValueError` when none exists |
| `dearomatize(mol)` | Deprecated alias of `kekulize`; issues a `DeprecationWarning` |
| `classify_pharmacophore(mol, atom_index)` | Pharmacophore bit mask of a graph atom: bit 0 hydrogen-bond donor, 1 acceptor, 2 positive ionisable, 3 negative ionisable, 4 aromatic, 5 hydrophobic |
| `implicit_h(atomic_num, bond_order_sum, formal_charge, aromatic=False)` | Implicit hydrogen count under the SMILES valence model; 0 when no default valence fits |

CIP priorities treat duplicate atoms as IUPAC 2013 Rule 1a states: at a phosphorus centre the oxo oxygen, which carries a duplicate phosphorus, ranks above a hydroxy or alkoxy oxygen, so phosphinates and phosphonates may receive the opposite descriptor from toolkits that rank it lower.

Every function in this table except `classify_pharmacophore` and `implicit_h`
accepts a molecule. `perceive_aromaticity`, `kekulize` and `dearomatize` change
a `MolGraph` argument in place and return it; SMILES and RDKit inputs are
converted to a new graph.

### Fingerprints

Binary fingerprints are sorted lists of set bit positions. Count fingerprints
are lists of `(position, count)` pairs for non-zero positions.

| Name | Description |
|---|---|
| `fingerprint(mol, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None)` | Fingerprint of a molecule; `kind` is `"mcs"`, `"path"`, `"ecfp"` or `"circular"` (ECFP4), `"fcfp"` (FCFP4) or `"torsion"`; a `FingerprintOptions` value as `options` defines the fingerprint instead |
| `fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp", *, path_length=7, options=None)` | Parse SMILES and fingerprint it; `mode` is `"ecfp"`, `"fcfp"`, `"path"` or `"mcs"`, or `options` defines the fingerprint |
| `FingerprintOptions(family="ecfp", **switches)` | A fingerprint specification shared with C++ and Java; see below |
| `circular_fingerprint(mol, radius=2, fp_size=2048, mode="ecfp")` | ECFP or FCFP of a graph; radius 2 gives ECFP4 |
| `circular_fingerprint_counts(mol, radius=2, fp_size=2048, mode="ecfp")` | Count form of `circular_fingerprint` |
| `path_fingerprint(mol, path_length=7, fp_size=2048)` | Heavy-atom paths of up to `path_length` bonds in a graph, hashed by atomic number |
| `mcs_fingerprint(mol, path_length=7, fp_size=2048, tautomer_aware=False)` | MCS-aware path fingerprint of a graph |
| `topological_torsion(mol, fp_size=2048)` | Four-atom torsion fingerprint of a molecule (Nilakantan et al., 1987) |
| `topological_torsion_counts(mol, fp_size=2048)` | Count form of the torsion fingerprint, for a graph |
| `fingerprint_subset(fp1, fp2)` | Whether every position of `fp1` is set in `fp2`; both sorted |
| `analyze_fp_quality(fp, fp_size=2048)` | Dictionary with `set_bits`, `fp_size`, `density`, `expected_collision_rate` and `is_saturated` (density above 0.75); a position outside `fp_size` raises `ValueError` |

### Configurable fingerprints

`FingerprintOptions` is a specification shared with C++ and Java: the same
options give the same bits in all three. `FingerprintOptions(family=...)`
starts from the classic configuration of the family, which reproduces the bits
of `circular_fingerprint`, `path_fingerprint` and `topological_torsion` exactly,
and applies the keyword arguments; `FingerprintOptions()` is the classic ECFP4,
and the class methods `ecfp(radius=2, size=2048)`, `fcfp(radius=2, size=2048)`,
`path(path_length=7, size=2048)`, `torsion(size=2048)` and `maccs()` build the
same configurations. Pass the value as the `options` of `fingerprint` or
`fingerprint_from_smiles`; the result is the set positions, or the
`(position, count)` pairs of the non-zero bins when `counts=True`.

```python
import smsd

chiral = smsd.FingerprintOptions(family="ecfp", radius=2, size=2048, cip_atoms=True, ring_size=True)
left = smsd.fingerprint("C[C@H](N)C(=O)O", options=chiral)
right = smsd.fingerprint("C[C@@H](N)C(=O)O", options=chiral)
assert left != right
assert smsd.fingerprint("C[C@H](N)C(=O)O", options=smsd.FingerprintOptions()) == smsd.fingerprint(
    "C[C@@H](N)C(=O)O", options=smsd.FingerprintOptions())

framework = smsd.FingerprintOptions(scaffold="murcko")
assert smsd.fingerprint("Cc1ccccc1", options=framework) == smsd.fingerprint("c1ccccc1", options=smsd.FingerprintOptions())

keys = smsd.fingerprint_from_smiles("CC(=O)Oc1ccccc1C(=O)O", options=smsd.FingerprintOptions.maccs())
assert 162 in keys  # MACCS key 163, a six-membered ring, is position 162
counts = smsd.fingerprint("CC(=O)Oc1ccccc1C(=O)O", options=smsd.FingerprintOptions("keys", keys=["[#8]", ("c1ccccc1", 1)], counts=True))
assert counts == [(0, 4), (1, 1)]
```

| Keyword | Meaning | ecfp | fcfp | path | torsion |
| --- | --- | --- | --- | --- | --- |
| `family` | `"ecfp"`, `"fcfp"`, `"path"`, `"torsion"` or `"keys"`, or a `FingerprintFamily` member | | | | |
| `radius`, `path_length`, `size` | Radius (−1 for the whole component), path length in bonds, and bits or bins | 2, –, 2048 | 2, –, 2048 | –, 7, 2048 | –, –, 2048 |
| `counts` | Return `(position, count)` pairs instead of positions | False | False | False | False |
| `element` | Atomic number | on | off | on | on |
| `degree` | Number of explicit neighbours | on | on | off | on |
| `valence` | Bond order sum; pi-electron count for `"torsion"` | on | off | off | on |
| `formal_charge`, `isotope`, `hydrogen_count`, `tautomer_class`, `aromaticity` | Formal charge, mass number, attached hydrogens, tautomer class when assigned, aromatic atom | on | off | off | off |
| `bond_order` | Bond order of the neighbour term, or of each path bond | on | on | off | off |
| `aromatic_bond`, `ring_bond` | Aromatic and ring bond terms | off | off | off | off |
| `ring_membership` | Ring atom | on | on | off | on |
| `ring_size`, `ring_size_cap` | Smallest SSSR ring through the atom, as a class; rings of the cap size or more share one class | off, 8 | off, 8 | off, 8 | off, 8 |
| `cip_atoms`, `cip_bonds` | CIP R/S (and r/s) of tetrahedral centres and E/Z of double bonds, as `assign_rs` and `assign_ez` give them | off | off | off | off |
| `fragments` | `"whole"`, `"largest"` (the component with the most heavy atoms, ties resolved without regard to atom order) or `"per_fragment"` (every atom term carries its component's class), or a `FragmentMode` member | `"whole"` | | | |
| `scaffold` | `"none"`, or `"murcko"`: the Bemis-Murcko framework, that is ring atoms, the linkers between rings and the atoms double- or triple-bonded to them, capped with hydrogen where side chains were; a molecule with no ring is used whole | `"none"` | | | |
| `keys` | The keys of the `"keys"` family, one bit each, as SMARTS strings, `(smarts, minimum)` pairs or `FingerprintKey(smarts, minimum=1)` values; a key is set when the molecule holds at least `minimum` distinct matched atom sets | | | | |

FCFP always hashes the pharmacophore feature class first. The stereo terms are
the CIP descriptors of the whole molecule, so a remote substituent that changes
a centre's priorities changes its label. They join the circular hashes after
the first iteration, so radius 0 does not separate stereoisomers and the bits
of a stereoisomer pair differ only within radius − 1 bonds of the centre or
bond; the path family folds R/S into its atom terms, so single-atom paths
already differ, and E/Z into its bond terms, and the torsion family folds both
into the atom types. Double bonds in rings of fewer than eight atoms carry no
E/Z term. The framework keeps exocyclic carbonyls and imines, as the usual
convention does, whereas `murcko_scaffold` removes every side chain; a
tetrahedral centre or double bond that stays stereogenic after capping keeps
its configuration. The scaffold is taken before the fragment choice. `FingerprintOptions.maccs()` is the key
family of the 166 MACCS structural keys, key k at position k − 1, where keys 1,
125 and 166 count isotopes, aromatic rings and fragments;
`FingerprintOptions.maccs_keys()` lists the definitions and each key's `kind`
names what it counts. Every attribute can be read and assigned; `copy()` returns
an independent copy. An unknown keyword raises `TypeError`, an unknown family,
fragment or scaffold name `ValueError`, and inconsistent options or a SMARTS key
that does not parse `ValueError` when the fingerprint is computed. A key search
that does not complete within its time budget, or that reaches 100,000
embeddings before its count is settled, raises `RuntimeError` rather than
returning bits that depend on the machine. The MCS fingerprint is a separate
function that these options do not configure.

### Similarity coefficients and encodings

| Name | Description |
|---|---|
| `tanimoto_coefficient(fp1, fp2, fp_size=2048)` | \|A ∩ B\| / \|A ∪ B\| over set positions; two empty fingerprints score 1.0 |
| `overlap_coefficient(fp1, fp2, fp_size=2048)` | \|A ∩ B\| / min(\|A\|, \|B\|); two empty fingerprints score 1.0, one empty fingerprint 0.0 |
| `overlapCoefficient` | Deprecated alias of `overlap_coefficient` |
| `dice(fp1, fp2, fp_size=2048)` | 2\|A ∩ B\| / (\|A\| + \|B\|) |
| `cosine(fp1, fp2, fp_size=2048)` | \|A ∩ B\| / √(\|A\| \|B\|) |
| `soergel(fp1, fp2, fp_size=2048)` | Distance, 1 − Tanimoto |
| `count_tanimoto_coefficient(fp1, fp2)` | Σ min(a, b) / Σ max(a, b) |
| `count_overlap_coefficient(fp1, fp2)` | Σ min(a, b) / min(Σ a, Σ b) |
| `count_overlapCoefficient` | Deprecated alias of `count_overlap_coefficient` |
| `count_dice(fp1, fp2)`, `count_cosine(fp1, fp2)` | Dice and cosine of count vectors given densely, as `(position, count)` pairs or as dictionaries |
| `counts_to_array(counts_dict, fp_size)` | Dense count list of length `fp_size` from a dictionary or `(position, count)` pairs; positions outside the range are ignored |
| `to_hex(fp, fp_size=2048)`, `from_hex(hex_str)` | Hexadecimal text and back; a non-hexadecimal character raises `ValueError` |
| `to_binary_string(fp, fp_size=2048)` | Text of `fp_size` `0` and `1` characters |

`tanimoto_coefficient` and `overlap_coefficient` also accept count
fingerprints and compare their positions only; their `fp_size` argument is
accepted for compatibility and does not change the result. The two sparse count
helpers accept `(position, count)` pairs or position lists, and return 0.0 when
both fingerprints are empty. `count_dice` and `count_cosine` need dense vectors,
for example from `counts_to_array`.

### Screening and target corpora

| Name | Description |
|---|---|
| `similarity(mol1, mol2)` | Upper bound, from atom-label counts, on the MCS similarity of two molecules, in [0, 1]; a screening value, not a fingerprint metric |
| `similarity_upper_bound(g1, g2, chemistry=None)` | The same bound for two molecules under the given chemistry, the default options when `None`; the `chem` keyword is a deprecated spelling |
| `screen_targets(query, targets, threshold, chemistry=None)` | Indices of the targets whose bound is at least `threshold`; the `chem` keyword is a deprecated spelling |
| `batch_similarity(query, targets, fp_size=2048, radius=2)` | `(index, Tanimoto)` pairs for a query molecule and target molecules on ECFP fingerprints, most similar first |
| `gpu_is_available()` | Whether a GPU backend is available at run time |
| `gpu_device_info()` | Description of the compute backend; it does not show which backend a particular call used |

`TargetCorpus` members:

| Member | Description |
|---|---|
| `TargetCorpus()`, `TargetCorpus.from_smiles(targets)` | An empty corpus, or one built from molecules |
| `add_target(target)`, `add_targets(targets)` | Capture molecules and their current atom numbering |
| `prewarm(*, num_threads=0)`, `is_prewarmed` | Compute invariants and ECFP4 fingerprints; whether they cover every target added |
| `screen(query, threshold, *, num_threads=0)` | Indices of targets with ECFP4 Tanimoto similarity of at least `threshold` |
| `len(corpus)` | Number of targets |

### MOL and SDF

| Name | Description |
|---|---|
| `read_mol_block(mol_block)` | Parse MOL V2000 or V3000 text |
| `read_mol_file(filename)` | Read a MOL file (V2000 or V3000) from a path; a file that cannot be opened raises `ValueError`. `read_molfile(filename)` is a deprecated alias that reads the file as UTF-8 text and calls `read_mol_block` |
| `read_sdf(filename)` | Read every record at once; malformed records become empty graphs |
| `write_mol_block(mol)`, `write_mol_block_v3000(mol)` | MOL text for a graph: V2000, or V3000 above 999 atoms or bonds; always V3000 |
| `write_sdf_record(mol)` | One SDF record for a graph, in the format of `write_mol_block` |
| `write_molfile(mol, filename, *, v3000=False, sdf=False)` | Write a molecule as a MOL block, a V3000 MOL block or one SDF record |
| `write_sdf(molecules, filename)` | Write molecules as SDF records, in the format of `write_mol_block`, to a path |
| `export_sdf(molecules, filename)` | Deprecated alias of `write_sdf` |

### Coordinates and layout

Functions return coordinates as lists of `(x, y)` or `(x, y, z)` tuples, one
per atom, and accept lists of `[x, y]` or `(x, y)` pairs. Refinement and
transform functions return new lists and leave their input unchanged.

| Name | Description |
|---|---|
| `generate_coords_2d(mol, target_bond_length=1.5)` | Two-dimensional coordinates for a graph |
| `generate_coords_3d(mol, target_bond_length=1.5)` | A three-dimensional starting geometry, not a validated conformer |
| `layout_quality(mol, coords, target_bond_length=1.5)` | Score combining bond-length uniformity, overlaps and crossings; 0.0 is ideal and lower is better |
| `count_crossings(mol, coords)` | Number of crossing bond pairs |
| `is_degenerate_layout(coords)` | Whether the layout has zero width or height |
| `layout_sssr(mol)` | Smallest set of smallest rings, ordered for layout |
| `match_template(mol, target_bond_length=1.5)` | Template coordinates for benzene, naphthalene, indole, purine, steroid, morphinan, quinoline, biphenyl, cyclohexane and piperidine; otherwise an empty list |
| `force_directed_layout(mol, coords, max_iter=500, target_bond_length=1.5)` | `(stress, coords)` after force-directed refinement |
| `stress_majorisation(mol, coords, max_iter=300, target_bond_length=1.5)` | `(stress, coords)` after stress majorisation |
| `reduce_crossings(mol, coords, max_iter=1000)` | `(crossings, coords)` after flipping ring systems |
| `resolve_overlaps(coords, threshold=0.3, max_iter=100)` | `(remaining_overlaps, coords)` |
| `translate_2d(coords, dx, dy)` | Translate |
| `rotate_2d(coords, angle)`, `scale_2d(coords, factor)` | Rotate by `angle` radians or scale, about the centroid |
| `mirror_x(coords)`, `mirror_y(coords)` | Reflect in the x axis or the y axis |
| `center_2d(coords)` | Centre on the origin |
| `align_2d(coords, reference)` | `(rmsd, aligned)` after optimal rotation onto `reference` |
| `bounding_box_2d(coords)` | `(min_x, min_y, max_x, max_y)` |
| `normalise_bond_length(mol, coords, target=1.5)` | Scale to the target bond length |
| `canonical_orientation(mol, coords)` | Rotate to a canonical orientation on a 30-degree grid |
| `Point2D(x, y)`, `Point3D(x, y, z)` | Coordinate value types with writable coordinate attributes |

`layout_sssr`, `match_template`, `force_directed_layout`, `stress_majorisation`
and `reduce_crossings` accept molecules; the other functions that take `mol`
need a graph.

### Depiction

| Name | Description |
|---|---|
| `depict_svg(mol_or_smiles, opts=None, **kwargs)` | SVG text for a molecule |
| `depict_mapping(mol_or_smiles, mapping, opts=None, **kwargs)` | SVG with mapped atoms highlighted and numbered; mapping keys are atom indices in the drawn molecule |
| `depict_pair(mol1, mol2, mapping, opts=None, **kwargs)` | Side-by-side SVG with the mapping from `mol1` to `mol2` highlighted |
| `save_svg(svg_str, path)` | Write SVG text to a file |
| `DepictOptions` | Native drawing options; see below |
| `depict_mcs(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000, size=(600, 300), highlight_color=(0.5, 0.8, 1.0))` | RDKit image of the MCS of two molecules with matched atoms highlighted; `highlight_color` is deprecated and has no effect, and a value other than the default warns |
| `depict_substructure(query, target, *, timeout_ms=10000, size=(600, 300))` | RDKit image of a substructure match between two molecules |
| `to_svg(mol_or_smiles, *, size=(300, 200))` | Deprecated; SVG drawn by RDKit. Use `depict_svg`, which needs no RDKit |

The native functions do not need RDKit. Their keyword arguments set
`DepictOptions` attributes; when `opts` is supplied, the attributes are set on
that object.

| `DepictOptions` attribute | Default |
|---|---|
| `bond_length`, `padding` | 30 pixels |
| `line_width`, `bold_width`, `font_size`, `match_width` | 0, derived automatically |
| `bond_spacing` | 0, giving 18% of `bond_length` |
| `highlight_radius` | 0, giving 40% of `bond_length` |
| `highlight_opacity` | 0.30 |
| `subscript_scale` | 0.70 |
| `show_carbon_labels`, `show_atom_indices` | `False` |
| `show_map_numbers` | `True` |
| `width`, `height` | 0, fitted to the molecule |
| `font_family` | `"Arial, Helvetica, sans-serif"` |

## Further reading

See the [Python package README](../python/README.md), the
[Python examples](EXAMPLES.md), the [installation guide](HOWTO-INSTALL.md), the
[C++ guide](CPP.md), the [Java guide](JAVA.md) and the
[changelog](../CHANGELOG.md).

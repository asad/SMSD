# SMSD Pro Python guide

Reference for the `smsd` Python package, version 8.0.3. The
[Python examples](EXAMPLES.md) cover fingerprints, drawing, layout, stereo,
files, R-groups and scaffolds.

## Installation

```bash
python -m pip install "smsd==8.0.3"
```

The CPU/OpenMP wheels cover CPython 3.10–3.14; the
[supported-platforms table](HOWTO-INSTALL.md#supported-platforms) lists the
operating systems and architectures. `smsd.__version__` gives the installed version. Core
matching needs no other toolkit. RDKit is needed only for RDKit molecules,
`to_rdkit`, `to_svg`, `depict_mcs` and `depict_substructure`
(`python -m pip install "smsd[rdkit]==8.0.3"`); Open Babel only for
`from_openbabel` and `to_openbabel`.

### Build from source

With Python 3.10 or later, a C++17 compiler and CMake 3.18 or later, from the
repository root:

```bash
python -m pip install -e ".[dev]" \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

This builds the CPU/OpenMP package. Without OpenMP, batch calls run sequentially.

### Optional CUDA screening

`smsd-cuda` provides optional NVIDIA CUDA screening on Linux x86_64 and needs
the same version of `smsd`:

```bash
python -m pip install "smsd==8.0.3" "smsd-cuda==8.0.3"
```

The GPU computes integer overlaps; scores use the CPU expressions and hits
keep target order. Matching runs on the CPU. The backend keeps prepared
corpus fingerprints on the device for reuse. ECFP4 Tanimoto
screening and the atom-count MCS similarity bound are distinct operations.

Automatic ECFP4 screening uses CUDA on a Tesla T4 with compute capability 7.5
when the prepared corpus has at least 100,000 targets. Smaller corpora and
other devices use the CPU. Label-bound and substructure candidate screening
use the CPU. `gpu_is_available()` and `gpu_device_info()` report backend
availability, not whether a call used the device. Missing devices, allocation
failure and device loss use the CPU.
Setting `SMSD_DISABLE_GPU` to anything but an empty string or `0` disables the
GPU at run time. macOS uses the CPU; the
[supported-platforms table](HOWTO-INSTALL.md#supported-platforms) also covers
the GPU Docker image and Windows through WSL2.

## Searches

A *molecule* is a SMILES string, `MolGraph` or RDKit molecule; a *graph* is a
`MolGraph` (parse SMILES with `parse_smiles`). Prepare a query once:

- `Substructure.find(query)`: a molecule to embed whole in each target.
- `Substructure.smarts(query)`: a SMARTS pattern; an empty or non-string
  pattern raises `ValueError`.
- `MCS.find(query)`: a molecule for a maximum common substructure (MCS) with
  each target.

A prepared query provides:

| Method | Returns |
|---|---|
| `match(target)` | A `SearchResult` with at most one match |
| `match_all(target, limit=1000)` | Up to `limit` matches: all optimal `MCS` mappings or all `Substructure` embeddings, symmetry-related ones included |
| `match_many(targets, *, num_threads=0)` | One `SearchResult` per target, in order, each with its own budget; `0` uses the OpenMP default (`OMP_NUM_THREADS`) |
| `with_options(options)`, `with_timeout(milliseconds)` | A new query with a copy of the options |
| `MCS.with_objective(objective)`, `MCS.with_fragments(*, min_atoms=1, max_fragments=2147483647)`, `MCS.with_atom_weights(weights)` | A new `MCS` query with that setting |
| `Substructure.with_induced(enabled=True)` | A new `Substructure` query with induced matching |

`limit` must be a positive integer and `num_threads` a non-negative integer,
or `ValueError` is raised. Searches release the global interpreter lock, so a
prepared query can serve several threads; do not modify a target during a call.

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

| `SearchResult` field | Meaning |
|---|---|
| `matches`, `best`, `found` | Tuple of `Match` objects; the first or `None`; whether any exists |
| `status` | `SearchStatus.COMPLETED`, `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED` |
| `optimal` | The MCS optimum, or the substructure answer, is proved |
| `complete` | The search or enumeration finished |
| `score` | Objective value, in weight units with atom weights |
| `nodes_visited`, `elapsed_ms` | Search effort and elapsed milliseconds |

`COMPLETED` means the search ran to its end; `TIMED_OUT` and `CANCELLED` mean
the time limit or the cancellation token stopped it. `LIMIT_REACHED` means the
node limit stopped it, an enumeration reached its match limit, or an MCES pair
beyond the [size limits](#size-limits) was declined before any search; a
declined pair returns no match and proves nothing.

With `complete_rings_only` or `RingFusionMode.STRICT`, truncated ring perception
in either input leaves the mappings selected using the perceived rings, but
cannot prove the result. The status is `LIMIT_REACHED`, with `optimal` and
`complete` both `False`; `TIMED_OUT` and `CANCELLED` take precedence.

A stopped search can return a valid but unproved match, and an empty
incomplete result does not prove absence, so check `status`, `optimal` and
`complete`. An enumeration limit can leave `optimal=True` with
`complete=False`. Tied matches are not ranked chemically.

```python
import smsd

result = smsd.Substructure.find("CC").match_all("CCC", limit=1)
assert len(result.matches) == 1 and result.optimal and not result.complete
assert result.status == smsd.SearchStatus.LIMIT_REACHED
```

Each `Match` has:

- `atom_mapping`: a read-only, one-to-one map from query to target atom index;
- `bond_mapping`: `(query_first, query_second, target_first, target_second)`
  for each selected bond;
- `fragments`: the components of the selected bonds, as query atom indices;
- `atom_count` and `bond_count` of the selected graph;
- `to_smiles(query)` and `to_graph(query)`: the selection as SMILES with cut
  valences capped, or as a `MolGraph` with the original chemistry and hydrogen
  counts.

### Size limits

An MCES search declines a pair before searching it, with `LIMIT_REACHED`, no
match, and `optimal` and `complete` both `False`, when:

- query atoms × target atoms exceed 2³⁰ (1,073,741,824), as for two molecules
  of more than 32,768 atoms each;
- the two molecules are too large in bonds, as for two molecules of more than
  16,384 bonds each;
- rarely, a molecule has an atom with 512 or more bonds and the two molecules
  share many kinds of bond.

MCIS searches are never declined; the time and node limits bound them. C++
and Java use the same limits.

## MCIS, MCES and fragments

```python
import smsd

query, target = "C1CC1", "CC(C)C"  # cyclopropane, isobutane
induced = smsd.MCS.find(query).match(target)
edges = smsd.MCS.find(query).with_objective(smsd.MCS.Objective.MCES).match(target)
assert (induced.best.atom_count, induced.best.bond_count) == (2, 1)
assert (edges.best.atom_count, edges.best.bond_count) == (3, 2)
assert len(smsd.parse_smiles(edges.best.to_smiles(query)).bonds()) == 2
fragments = smsd.MCS.find("CCOCC").with_fragments().match("CCNCC")
assert fragments.best.atom_count == 4 and len(fragments.best.fragments) == 2
all_matches = smsd.MCS.find(query).match_all(target, limit=100)
assert all_matches.optimal and all_matches.complete
```

- `MCS` defaults to a connected maximum common induced subgraph (MCIS): most
  atoms, then most bonds, keeping bonds and non-bonds between mapped atoms.
- `MCS.Objective.MCES`, the maximum common edge subgraph, maximises common
  bonds, then their atoms; it may omit input bonds and adds no isolated atoms.
- `with_fragments()` allows several components of at least `min_atoms` atoms
  each, so MCIS may include isolated atoms when `min_atoms=1`. Limits below
  one raise `ValueError`.
- `Substructure` is non-induced unless `with_induced()` or
  `Substructure.Options.induced` is set.

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

`ChemOptions`, as `options.chemistry`, decides which atoms and bonds match:

| Attribute | Default | Effect |
|---|---|---|
| `match_atom_type` | `True` | Same element |
| `match_bond_order` | `BondOrderMode.STRICT` | Equal orders, subject to `aromaticity_mode`; `LOOSE` and `ANY` relax this |
| `aromaticity_mode` | `AromaticityMode.FLEXIBLE` | Aromatic bonds also match single or double bonds; `STRICT` needs agreement |
| `aromaticity_model` | `AromaticityModel.DAYLIGHT_LIKE` | Deprecated; no search reads it |
| `match_formal_charge` | `False` | Same formal charge |
| `match_isotope` | `False` | Same mass number where both atoms give one |
| `use_chirality`, `use_bond_stereo` | `False` | Mapped tetrahedral centres or double bonds agree in configuration, by the parity of their mapped neighbours, not CIP labels |
| `ring_matches_ring_only` | `False` | Ring atoms and bonds match only ring atoms and bonds |
| `complete_rings_only` | `False` | Complete rings only, in each input for `MCS` |
| `ring_fusion_mode` | `RingFusionMode.IGNORE` | `STRICT`: matched ring atoms lie in the same number of rings; `PERMISSIVE` currently equals `IGNORE` |
| `tautomer_aware` | `False` | Allow recognised tautomeric relaxations; `match_bond_order` is unchanged |
| `pH`, `solvent` | `7.4`, `Solvent.AQUEOUS` | Context for `MCS.compute_tautomer_confidence` only |

Stereo is checked only where the mapping determines it. In `Substructure`, a
tagged query centre or double bond needs the same configuration in the target,
and an untagged one matches any; in `MCS`, an unspecified target atom or bond
matches either configuration. The comparison uses the parity of the mapped
neighbours, not CIP labels, whose priorities depend on atoms outside the
mapping. In `MCS` it is a hard constraint: an MCIS is the largest common
subgraph whose mapped stereo agrees, so it can be smaller than the MCIS found
without these options, and an MCES leaves out the fewest common bonds that
resolve each disagreement.

Tautomer matching usually also needs a relaxed bond order:

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

`ChemOptions.profile(name="default")` returns a preset; names are not
case-sensitive, and an unknown name raises `ValueError`:

| Preset | Differs from the defaults |
|---|---|
| `default` | No differences |
| `strict` | `aromaticity_mode=STRICT` |
| `pharma` | `match_formal_charge`, `ring_matches_ring_only`, `complete_rings_only` |
| `reaction`, `compat-fmcs`, `compat-substruct` | `match_bond_order=LOOSE` |
| `ChemOptions.tautomer_profile()` | `tautomer_aware=True`, `match_bond_order=LOOSE` |

`ChemOptions(options)`, `copy.copy` and `copy.deepcopy` copy options. The
setters `with_bond_order_mode`, `with_aromaticity_mode`,
`with_complete_rings_only`, `with_match_isotope`, `with_ring_fusion_mode`,
`with_tautomer_aware`, `with_solvent` and `with_ph` set one field and return the
same object; `with_ph` clamps a finite pH to [0, 14].

## Timeouts, node limits and cancellation

- `timeout_ms`, default `-1`: a negative value selects 10 s for `Substructure`
  and 500 + 2 × (query atoms × target atoms) ms, at most 30 s, for `MCS`: 1.3 s
  for two 20-atom molecules, 5.5 s for two 50-atom molecules and 30 s from 122
  atoms each; `0` stops at once; a positive value applies to each query–target
  pair.
- `node_limit`, default `0` (none): reaching it reports `LIMIT_REACHED`.
- `cancellation`, default `None`: a `CancellationToken`; `cancel()` from any
  thread stops the search with `CANCELLED`. Tokens stay cancelled, so use a new
  one per search.

Parsing and preparation can take the total time past the deadline. Converting
matches to Python objects counts towards it in `match` and `match_all`, which
then return the matches converted so far with `TIMED_OUT` or `CANCELLED` and
`complete=False`; `match_many` stops conversion only on cancellation.

## Atom weights

`with_atom_weights(weights)`, for MCIS only, ranks results by the sum of query
atom weights before atom and bond counts; negative weights can favour a smaller
result. Give one finite weight per query atom, in input order. Weights are
compared to 0.001; `result.score` is in the supplied units. Non-finite weights,
totals beyond about ±2.1 million, a list missing a query atom, and use with
MCES raise `ValueError`.

## Tautomer confidence

`MCS.compute_tautomer_confidence(query, target, match, chemistry)` scores the
tautomeric relaxations a `Match` uses, in the `pH` and `solvent` of
`chemistry`; identical selected chemistry scores 1. The score is a heuristic,
not a probability, and changes neither molecules nor search. A plain atom
mapping raises `TypeError`.

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

`MCS.validate_tautomer_consistency(query, target, match)` reports whether a
`Match` or atom mapping conserves the mobile N–H, O–H and S–H protons of each
mapped tautomeric group.

## SMARTS and target collections

```python
import smsd

query = smsd.Substructure.smarts("[#6]-[#8]")
assert [result.found for result in query.match_many(["CCO", "CCN"])] == [True, False]
corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN"])
corpus.prewarm(num_threads=1)
assert [result.found for result in query.match_many(corpus, num_threads=2)] == [True, False]
```

SMARTS queries follow the result, timeout and enumeration rules of molecule
queries. `TargetCorpus` holds targets for repeated queries:

- `TargetCorpus()`, `TargetCorpus.from_smiles(targets)`, `add_target(target)`,
  `add_targets(targets)`: capture molecules and their atom numbering,
  unaffected by later changes to the inputs; `len(corpus)` counts them.
- `prewarm(*, num_threads=0)` computes invariants and ECFP4 fingerprints;
  `is_prewarmed` says whether they cover every target.
- `screen(query, threshold, *, num_threads=0)`: indices of targets with ECFP4
  Tanimoto of at least `threshold`, empty until `prewarm()` follows the last
  addition. ECFP similarity is not a substructure test.
- `screen_many(queries, threshold, *, num_threads=0)`: the same ECFP4 screen
  for each query, returning a list of hit-index lists in query order. Each
  hit list keeps target order and needs the same `prewarm()`.
- `screen_targets_many(queries, threshold, chemistry=None, *, num_threads=0)`:
  the atom-count MCS similarity-bound screen for each query, with one chemistry
  setting for the batch. No `prewarm()` is needed. Each hit list keeps target
  order, as `screen_targets(query, corpus, threshold, chemistry)` does.

```python
import smsd

corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "CCC"])
corpus.prewarm(num_threads=1)
assert corpus.screen_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
assert corpus.screen_targets_many(["CCO", "CCN"], 1.0, num_threads=1) == [[0], [1]]
```

Corpus substructure searches can reject impossible molecule candidates with a
conservative element-path screen before CPU matching. SMARTS queries, untyped
queries and finite node limits bypass this screen. Target results keep their
original slots, atom numbering and search status.

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

- `MCS.match_constrained(queries, target, options=None)`: one `SearchResult`
  per query, in order, each excluding target atoms claimed earlier; the joint
  mapping is not proved optimal.
- `MCS.Options.excluded_target_atoms` reserves further target atoms without
  changing the target.
- `MCS.consensus(molecules, options=None, threshold=1.0)`: a common graph in at
  least `ceil(threshold * len(molecules))` inputs, `threshold` in (0, 1], under
  one budget; always `optimal=False`. Atom weights and exclusions raise
  `ValueError`. With complete rings or strict ring fusion, truncated ring
  perception in any input makes it `LIMIT_REACHED` and incomplete, unless it is
  already `TIMED_OUT` or `CANCELLED`; its fragment and supporting mappings are
  retained.
- `ConsensusResult` has `graph`, `to_smiles()`, `support_count`, `mappings`
  (supporting input index to a `Match` from `graph`), `fragment` and
  `source_index` (the graph in its source), `found`, `status`, `optimal`,
  `complete`, `nodes_visited` and `elapsed_ms`.

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

- `MCS.validate_mapping(query, target, match, options=None)`: a list of errors,
  empty when valid. A plain atom mapping selects every compatible mapped bond,
  so pass a `Match` for MCES.
- `MCS.is_mapping_maximal(query, target, match, options=None)`: whether no
  compatible extension exists, which is weaker than optimality. An invalid
  mapping gives `False`; a valid mapping raises `RuntimeError` if the proof
  cannot finish or required ring perception was truncated.
- `MCS.canonicalize_mapping(query, target, atom_mapping)`: the
  lexicographically smallest equivalent mapping under both molecules'
  symmetries; raises if it cannot finish.

## R-group decomposition

`Substructure.decompose_r_groups(core, targets, options=None)` returns, per
target, a dictionary with `core` (a `MolGraph`), `rgroups` (labels such as
`"R1"` mapped to graphs) and `search` (the core's `SearchResult`, to check
first); see the [R-group example](EXAMPLES.md#r-group-decomposition).

## Input atom numbering

SMILES use parser order, `MolGraph` objects graph order and RDKit molecules
their own order, separately on each side; search results and atom weights use
these indices. Other functions number an RDKit molecule as `from_rdkit` does;
`get_index_map(from_rdkit(mol))` gives the RDKit index of each atom. Explicit
hydrogen atoms are vertices; bracket hydrogen counts are not.

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
predicates are lost. For C++, see the
[optional native adapter](CPP.md#optional-native-adapter).

- `from_rdkit(mol, use_cache=True)`: a graph with the atom correspondence,
  reused while the molecule is unchanged; empty or hydrogen-only molecules
  raise `ValueError`.
- `to_rdkit(mol_or_smiles)`: an RDKit molecule in canonical SMILES order.
- `get_index_map(g)`: a copy of `result[smsd_index] = rdkit_index` for a
  `from_rdkit` graph, else `None`.
- `translate_mapping(mapping, g_query=None, g_target=None)`: a mapping between
  `from_rdkit` graphs in RDKit indices; an out-of-range index raises
  `IndexError`.
- `clear_cache()`: clears reusable conversions; live graphs keep their
  correspondence.
- `from_openbabel(mol)`, `to_openbabel(mol)`: Pybel molecule or `OBMol` to
  graph, and back to Pybel, without the Open Babel numbering.

## API reference

Every name in `smsd.__all__` appears in this guide.

### Enumerations and search options

| Name | Members |
|---|---|
| `MCS.Objective`, `MCS.Connectivity` | `MCIS` (default), `MCES`; `CONNECTED` (default), `FRAGMENTS` |
| `SearchStatus` | `COMPLETED`, `TIMED_OUT`, `CANCELLED`, `LIMIT_REACHED` |
| `BondOrderMode`, `AromaticityMode`, `AromaticityModel` | `STRICT`, `LOOSE`, `ANY`; `STRICT`, `FLEXIBLE`; `DAYLIGHT_LIKE` |
| `RingFusionMode` | `IGNORE`, `PERMISSIVE`, `STRICT` |
| `Solvent` | `AQUEOUS`, `DMSO`, `METHANOL`, `CHLOROFORM`, `ACETONITRILE`, `DIETHYL_ETHER` |
| `CancellationToken` | `cancel()`, `is_cancelled` |

`MCS.Options` and `Substructure.Options` (copied by `MCS.Options(other)` and
`Substructure.Options(other)`):

| Field | Default | Applies to |
|---|---|---|
| `chemistry`, `timeout_ms`, `node_limit`, `cancellation` | `ChemOptions()`, `-1`, `0`, `None` | Both |
| `induced` | `False` | `Substructure` |
| `objective`, `connectivity` | `MCS.Objective.MCIS`, `MCS.Connectivity.CONNECTED` | `MCS` |
| `min_fragment_atoms`, `max_fragments` | `1`, `2147483647` | `MCS` |
| `atom_weights`, `excluded_target_atoms` | `[]`, `[]` | `MCS`; weights with MCIS only |

### Molecules and SMILES

- `parse_smiles(smiles, *, lenient=False)`: a graph; `""` gives an empty one.
  Malformed text raises `ValueError` unless `lenient=True` recovers it; bracket
  atoms beyond the OpenSMILES limits (one-digit hydrogen count, charge within
  ±15, three-digit isotope) always raise it.
- `parse_mapped_smiles(smiles)`: parse without atom maps;
  `strip_atom_maps(smiles)` removes them from the text, keeping ring closures
  and aromatic bonds.
- `to_smiles(mol)`: canonical SMILES, whose order may differ between
  languages; `canonical_smiles` is a deprecated alias. `/` and `\` go on single
  bonds wherever a configured double bond has one at an end; an end with only
  aromatic ring bonds carries them on an aromatic bond, as RDKit writes it,
  which strict OpenSMILES readers take as a single bond.
- `to_smarts(mol, include_aromaticity=True, include_charge=True, include_ring_member=False, include_h_count=False, include_isotope=False)`:
  SMARTS for a molecule with the selected predicates.
- `canonical_hash(mol)`: a hash of connectivity, bond orders, aromaticity,
  rings, charges, mass numbers and stereo; a screening key, not proof of
  identity. `same_canonical_graph(mol1, mol2)` compares atom count, canonical
  hash and canonical SMILES.
- `prewarm(mol)`: caches invariants and returns the graph; `prewarm_graph` is
  a deprecated alias returning `None`.
- `extract_subgraph(mol, atom_indices)`: the induced subgraph, numbered in the
  supplied order.
- `murcko_scaffold(mol)`: the Bemis–Murcko scaffold (ring systems, linkers and
  atoms double- or triple-bonded to them), other substituents replaced by
  hydrogen and stereo kept where defined; acyclic graphs are returned
  unchanged.
- `count_components(mol)`, `split_components(mol)`: the number of connected
  components, or one graph for each.
- `translate_to_atom_ids(mapping, g1, g2)`: a mapping in the IDs set by
  `MolGraphBuilder.atom_ids`; unchanged without IDs.

`MolGraph` members:

- `n`, `len(graph)`: atom count.
- `atomic_num`, `formal_charge`, `mass_number`, `hydrogen_count`,
  `atom_class`, `aromatic`, `ring`, `degree`, `neighbors`, `label`: per-atom
  values, adjacency lists and matching labels.
- `bonds()` as `(i, j, order)` tuples; `bond_order(i, j)`, `has_bond(i, j)`,
  `bond_in_ring(i, j)`, `bond_aromatic(i, j)`, raising `IndexError` out of
  range.
- `name`, `program_line`, `comment`, `properties`: writable MOL header lines
  and SDF data fields.
- `num_rings()`, `ensure_ring_counts()`, `ring_system_of(atom)`,
  `ring_system_count()`, `ring_system_sig(rs_id)`: rings and ring systems.
- `get_canonical_hash()`, `ensure_canonical()`, `morgan_rank`,
  `canonical_label`, `orbit`: canonical ranks, labels and symmetry classes,
  computed on first access.
- `get_automorphism_generators()`, `automorphism_generators_truncated()`:
  automorphism generators, and whether they were capped.
- `tautomer_class`: tautomer class per atom from 0, or −1; independent of atom
  order except between symmetric classes.
- `prewarm()` caches invariants; `perceive_aromaticity(model=AromaticityModel.DAYLIGHT_LIKE)`,
  `kekulize()` and `dearomatize()` change the graph in place.

`MolGraphBuilder` builds a `MolGraph` through chainable setters:
`atom_count(n)`, `atomic_numbers(nums)`, `formal_charges(charges)`,
`mass_numbers(masses)`, `hydrogen_counts(counts)`, `ring_flags(flags)`,
`aromatic_flags(flags)`, `neighbors(adj)`, `bond_orders(orders)` (parallel to
`neighbors`), `bond_ring_flags(flags)`, `bond_aromatic_flags(flags)`,
`atom_ids(ids)`, `name(value)`, `program_line(value)`, `comment(value)`,
`properties(values)`, then `build(perceive_aromaticity=True,
aromaticity_model=AromaticityModel.DAYLIGHT_LIKE)`.

### Stereo, rings and atom chemistry

- `assign_rs(mol)`, `assign_ez(mol)`: `{atom: "R" or "S"}` and
  `{(atom, atom): "E" or "Z"}` where priorities resolve.
- `assign_cip(mol)` returns `CIPDescriptors` with `rs_labels`, one `RSLabel`
  per atom (`NONE`, `R`, `S`, pseudoasymmetric `r`, `s`), and `ez_bonds`,
  `(atom, atom, EZLabel)` tuples (`NONE`, `E`, `Z`).
- `compute_sssr(mol)`: the smallest set of smallest rings, as atom lists by
  ascending size.
- `rings_truncated(mol)`: whether the ring search behind `compute_sssr` stopped
  at its work bound, as for a very large ring system, so that the rings may be
  fewer than the graph holds.
- `perceive_aromaticity(mol, *, model=AromaticityModel.DAYLIGHT_LIKE)`:
  recompute ring membership and aromaticity.
- `kekulize(mol)`: replace aromatic bonds with Kekulé bonds, or raise
  `ValueError`; `dearomatize(mol)` is a deprecated alias issuing a
  `DeprecationWarning`.
- `classify_pharmacophore(mol, atom_index)`: bit mask of a graph atom: 0 donor,
  1 acceptor, 2 positive and 3 negative ionisable, 4 aromatic, 5 hydrophobic.
- `implicit_h(atomic_num, bond_order_sum, formal_charge, aromatic=False)`:
  implicit hydrogens by SMILES valence; 0 when none fits.

The other functions here accept a molecule; `perceive_aromaticity`,
`kekulize` and `dearomatize` change a `MolGraph` in place and return it. CIP
labels follow IUPAC 2013: at phosphorus, oxo outranks hydroxy and alkoxy
oxygen, so phosphinates and phosphonates can differ from toolkits that rank it
lower.

### Fingerprints

Binary fingerprints are sorted set positions; count fingerprints are
`(position, count)` pairs.

- `fingerprint(mol, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None)`:
  `kind` is `"mcs"`, `"path"`, `"ecfp"` or `"circular"` (ECFP4), `"fcfp"`
  (FCFP4) or `"torsion"`; `options` overrides it.
- `fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode="ecfp", *, path_length=7, options=None)`:
  `mode` takes the `kind` values above; `options` overrides it.
- `fingerprints(mols, *, kind="mcs", path_length=7, fp_size=2048, radius=2, options=None, num_threads=0)`:
  the fingerprint of every molecule, each as `fingerprint` gives it, in input
  order. SMILES strings are parsed and the fingerprints computed natively on
  `num_threads` threads (0 uses the OpenMP default); the first molecule that
  fails raises the error `fingerprint` gives for it.
- `circular_fingerprint(mol, radius=2, fp_size=2048, mode="ecfp")`,
  `circular_fingerprint_counts(mol, radius=2, fp_size=2048, mode="ecfp")`: ECFP
  or FCFP of a graph; radius 2 is ECFP4.
- `path_fingerprint(mol, path_length=7, fp_size=2048)`: heavy-atom paths of up
  to `path_length` bonds in a graph.
- `mcs_fingerprint(mol, path_length=7, fp_size=2048, tautomer_aware=False)`:
  an MCS-aware path fingerprint of a graph.
- `topological_torsion(mol, fp_size=2048)`,
  `topological_torsion_counts(mol, fp_size=2048)`: four-atom torsions
  (Nilakantan et al., 1987) of a molecule, or counts for a graph.
- `fingerprint_subset(fp1, fp2)`: whether every position of `fp1` is in `fp2`;
  both sorted.
- `analyze_fp_quality(fp, fp_size=2048)`: `set_bits` (distinct positions),
  `fp_size`, `density` (set bits over `fp_size`), `expected_collision_rate`
  and `is_saturated` (density over 0.75); a position beyond `fp_size` raises
  `ValueError`. With s of N bits set, the fingerprint holds about
  n = ln(1 − s/N) / ln(1 − 1/N) uniformly hashed features, and the collision
  rate is the share 1 − s/n of them that fell on a bit already set: 0 for an
  empty fingerprint and 1 for a full one. Path fingerprints of a power-of-two
  size collide more often than uniform hashing predicts, so the rate
  understates their collisions. Java `analyzeFpQuality` uses the same
  definitions.

### Configurable fingerprints

`FingerprintOptions(family="ecfp", **switches)` gives the same bits in Python,
C++ and Java. Each family starts from the classic configuration that
reproduces `circular_fingerprint`, `path_fingerprint` and `topological_torsion`;
`FingerprintOptions()` is ECFP4. The class methods `ecfp(radius=2, size=2048)`,
`fcfp(radius=2, size=2048)`, `path(path_length=7, size=2048)`,
`torsion(size=2048)` and `maccs()` build the same configurations. Pass the
value as `options` to `fingerprint`, `fingerprints` or
`fingerprint_from_smiles`.

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

Keywords, with their defaults for ecfp, fcfp, path and torsion:

- `family`: `"ecfp"` (default), `"fcfp"`, `"path"`, `"torsion"`, `"keys"`, or
  a `FingerprintFamily` (`ECFP`, `FCFP`, `PATH`, `TORSION`, `KEYS`).
- `radius` (2, 2, –, –; −1 for the whole component), `path_length` in bonds
  (–, –, 7, –), `size` in bits or bins (2048), and `counts` to return
  `(position, count)` pairs (False).
- Atom terms: `element` (on, off, on, on); `degree`, explicit neighbours (on,
  on, off, on); `valence`, the bond order sum or pi electrons for torsion (on,
  off, off, on); `formal_charge`, `isotope`, `hydrogen_count`,
  `tautomer_class` and `aromaticity` (on, off, off, off); `ring_membership`
  (on, on, off, on); `ring_size`, the smallest ring through the atom, with
  rings at or above `ring_size_cap` (8) sharing a class (off).
- Bond terms: `bond_order` of neighbour or path bonds (on, on, off, off);
  `aromatic_bond` and `ring_bond` (off).
- `cip_atoms`, `cip_bonds`: R/S (and r/s) and E/Z as from `assign_rs` and
  `assign_ez` (off).
- `fragments` (ecfp, `"whole"`): `"whole"`, `"largest"` (most heavy atoms,
  independent of atom order), `"per_fragment"` (terms carry the component), or
  a `FragmentMode` (`WHOLE`, `LARGEST`, `PER_FRAGMENT`).
- `scaffold` (ecfp, `"none"`): `"none"`, `"murcko"` (the Bemis–Murcko
  framework, taken before `fragments`; ring-free molecules stay whole), or a
  `ScaffoldMode` (`NONE`, `MURCKO`).
- `keys`, for `"keys"`: SMARTS strings, `(smarts, minimum)` pairs or
  `FingerprintKey(smarts, minimum=1)`; a bit is set when at least `minimum`
  distinct atom sets match.

FCFP always includes the pharmacophore class. Stereo terms use whole-molecule
CIP labels; radius 0 does not separate stereoisomers, and double bonds in rings
under eight atoms carry no E/Z term. `maccs()` gives the 166 MACCS keys, key k
at position k − 1; `maccs_keys()` lists them as `FingerprintKey` values with
`smarts`, `minimum` and `kind` (what the key counts). Attributes are readable
and writable, and `copy()` returns an independent copy. An unknown keyword
raises `TypeError`; an unknown family, fragment or scaffold name raises
`ValueError`, as do inconsistent options and unparsable keys when computing. A
key count not settled in time raises `RuntimeError`. `mcs_fingerprint` ignores
these options.

### Similarity coefficients and encodings

- `tanimoto_coefficient(fp1, fp2, fp_size=2048)`: |A ∩ B| / |A ∪ B|; 1.0 for
  two empty fingerprints.
- `overlap_coefficient(fp1, fp2, fp_size=2048)`: |A ∩ B| / min(|A|, |B|); 1.0
  for two empty, 0.0 for one. `overlapCoefficient` is a deprecated alias.
- `dice(fp1, fp2, fp_size=2048)`, `cosine(fp1, fp2, fp_size=2048)`:
  2|A ∩ B| / (|A| + |B|) and |A ∩ B| / √(|A| |B|).
- `soergel(fp1, fp2, fp_size=2048)`: 1 − Tanimoto.
- `count_tanimoto_coefficient(fp1, fp2)`: Σ min(a, b) / Σ max(a, b).
- `count_overlap_coefficient(fp1, fp2)`: Σ min(a, b) / min(Σ a, Σ b);
  `count_overlapCoefficient` is a deprecated alias.
- `count_dice(fp1, fp2)`, `count_cosine(fp1, fp2)`: Dice and cosine of dense
  count vectors of equal length, `(position, count)` pairs or dictionaries.
- `counts_to_array(counts_dict, fp_size)`: a dense list of length `fp_size`
  from a dictionary or pairs, ignoring positions outside it.
- `to_hex(fp, fp_size=2048)`, `from_hex(hex_str)`: hexadecimal text and back;
  a non-hexadecimal character raises `ValueError`.
- `to_binary_string(fp, fp_size=2048)`: `fp_size` characters of `0` and `1`.

`tanimoto_coefficient` and `overlap_coefficient` compare positions only, even
for counts, and ignore `fp_size`. The two sparse count coefficients take pairs
or position lists and return 0.0 when both are empty.

### Screening

- `similarity(mol1, mol2)`: an upper bound in [0, 1] on MCS similarity, for
  screening; not a fingerprint metric.
- `similarity_upper_bound(g1, g2, chemistry=None)`: the same bound under
  `chemistry`; `screen_targets(query, targets, threshold, chemistry=None, *,
  num_threads=0)`: indices of targets whose bound reaches `threshold`, in
  target order; `0` uses the OpenMP default (`OMP_NUM_THREADS`), as in
  `match_many`. For both, `chem` is a deprecated keyword. To screen the same
  targets repeatedly, pass a prepared `TargetCorpus` as `targets`, so that
  they are not converted again on each call.
- `batch_similarity(query, targets, fp_size=2048, radius=2)`:
  `(index, Tanimoto)` pairs on ECFP, most similar first.
- `gpu_is_available()`, `gpu_device_info()`: whether a GPU backend is
  available, and a description of the run-time compute backend (the CPU when
  there is no GPU), not of the backend a call used.

### MOL and SDF

- `read_mol_block(mol_block)`, `read_mol_file(filename)`: V2000 or V3000 text
  or file; an unopenable file or a block that cannot be read, such as a V3000
  block with more than 10,000 atoms, raises `ValueError`; V2000 numeric fields
  are read leniently. `read_molfile` is a deprecated alias.
- `read_sdf(filename)`: every record at once; malformed records become empty
  graphs.
- `write_mol_block(mol)` (V2000, or V3000 above 999 atoms or bonds),
  `write_mol_block_v3000(mol)` and `write_sdf_record(mol)`: MOL text or one SDF
  record.
- `write_molfile(mol, filename, *, v3000=False, sdf=False)`: write a MOL block,
  V3000 block or SDF record.
- `write_sdf(molecules, filename)`: write SDF records; `export_sdf` is a
  deprecated alias.

### Coordinates and layout

Coordinates are lists of `(x, y)` or `(x, y, z)` tuples, one per atom; `[x, y]`
lists are accepted, and new lists are returned. Every function taking `mol`
accepts a molecule.

- `generate_coords_2d(mol, target_bond_length=1.5)`,
  `generate_coords_3d(mol, target_bond_length=1.5)`: 2D coordinates, or a 3D
  starting geometry that is not a validated conformer.
- `layout_quality(mol, coords, target_bond_length=1.5)`: a bond-length,
  overlap and crossing score; 0.0 is ideal.
- `count_crossings(mol, coords)`, `is_degenerate_layout(coords)`: crossing
  bond pairs; whether width or height is zero.
- `layout_sssr(mol)`: the smallest set of smallest rings, ordered for layout.
- `match_template(mol, target_bond_length=1.5)`: coordinates for benzene,
  naphthalene, indole, purine, steroid, morphinan, quinoline, biphenyl,
  cyclohexane or piperidine, else `[]`.
- `force_directed_layout(mol, coords, max_iter=500, target_bond_length=1.5)`,
  `stress_majorisation(mol, coords, max_iter=300, target_bond_length=1.5)`:
  `(stress, coords)`.
- `reduce_crossings(mol, coords, max_iter=1000)`: `(crossings, coords)` after
  flipping ring systems; `resolve_overlaps(coords, threshold=0.3, max_iter=100)`:
  `(remaining_overlaps, coords)`.
- `translate_2d(coords, dx, dy)` translates; `rotate_2d(coords, angle)`, in
  radians, and `scale_2d(coords, factor)` act about the centroid.
- `mirror_x(coords)`, `mirror_y(coords)`: reflect in the x or y axis;
  `center_2d(coords)` centres on the origin.
- `align_2d(coords, reference)`: `(rmsd, aligned)`;
  `bounding_box_2d(coords)`: `(min_x, min_y, max_x, max_y)`.
- `normalise_bond_length(mol, coords, target=1.5)` scales to the bond length;
  `canonical_orientation(mol, coords)` rotates onto a 30-degree grid.
- `Point2D(x, y)`, `Point3D(x, y, z)`: value types with writable coordinates.

### Depiction

- `depict_svg(mol_or_smiles, opts=None, **kwargs)`: SVG text.
- `depict_mapping(mol_or_smiles, mapping, opts=None, **kwargs)`: SVG
  highlighting and numbering mapped atoms, keyed by drawn-atom index.
- `depict_pair(mol1, mol2, mapping, opts=None, **kwargs)`: a side-by-side SVG
  of a mapping.
- `save_svg(svg_str, path)`: write SVG text to a file.
- `depict_mcs(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000, size=(600, 300), highlight_color=(0.5, 0.8, 1.0))`:
  an RDKit image of the MCS; `highlight_color` is deprecated, ignored, and
  warns when changed.
- `depict_substructure(query, target, *, timeout_ms=10000, size=(600, 300))`:
  an RDKit image of a substructure match.
- `to_svg(mol_or_smiles, *, size=(300, 200))`: deprecated RDKit SVG; use
  `depict_svg`.

The native SVG functions need no RDKit; keyword arguments set `DepictOptions`
attributes, on `opts` when given. Defaults: `bond_length` and `padding` 30 px;
`line_width`, `bold_width`, `font_size` and `match_width` 0 (automatic);
`width` and `height` 0 (fitted to the molecule); `bond_spacing` and
`highlight_radius` 0 (18% and 40% of `bond_length`); `highlight_opacity` 0.30;
`subscript_scale` 0.70; `show_carbon_labels` and `show_atom_indices` `False`;
`show_map_numbers` `True`; `font_family` `"Arial, Helvetica, sans-serif"`.

## Further reading

The [Python package README](../python/README.md),
[installation guide](HOWTO-INSTALL.md), [C++ guide](CPP.md),
[Java guide](JAVA.md) and [changelog](../CHANGELOG.md).

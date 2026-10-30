# SMSD Pro 8.1.0 API reference

This reference describes the 8.1.0 source API. Use a matching wheel or build;
an earlier installed native extension cannot provide new Python APIs.

The Python package provides molecular graphs, exact searches, fingerprints,
screening, depiction and optional data adapters. Java and C++ expose related
graph operations; their names and available options are documented in the
[Java guide](../java/README.md) and [C++ guide](CPP.md). Identical fingerprint
parameters do not establish bit identity across toolkits.

## MCS and substructure searches

```python
import smsd

options = smsd.MCS.Options()
options.timeout_ms = 1000
options.node_limit = 100000
result = smsd.MCS.find("CCO").with_options(options).match("CCN")
assert result.found and result.best.atom_count == 2
assert smsd.MCS.validate_mapping("CCO", "CCN", result.best, options) == []
```

`MCS.Objective.MCIS` maximises common atoms while preserving edges and
non-edges. `MCES` maximises selected common bonds; its witness includes the
chosen bonds as well as atoms. The objectives are different and should not be
compared as interchangeable scores. `MCS.Connectivity.CONNECTED` requires a
connected common graph; `FRAGMENTS` permits disconnected selected components.
`with_fragments(min_atoms=..., max_fragments=...)` selects the latter mode.
`with_atom_weights` assigns weights to query atoms in MCIS; weighted scores
are not atom counts. Use the matching options when validating a witness.

Prepared searches copy options and support `match`, `match_all(limit=...)`
and `match_many(num_threads=...)`. A prepared query can be reused; avoid
mutating its source molecule or shared mutable inputs during calls. Options
are copied, but a shared cancellation token intentionally remains a signal
to all searches using it. `num_threads=0` selects the native OpenMP default;
set an explicit positive thread count for a controlled batch.

`SearchResult` contains `matches`, `status`, `optimal`, `complete`, `score`,
`nodes_visited` and `elapsed_ms`; `found` and `best` are convenience properties.
`best` is `None` when there is no witness. `Match` contains `atom_mapping`,
selected `bond_mapping` and `fragments`, with `atom_count` and `bond_count`.
For SMILES inputs, atom indices are zero-based in parsing order. Foreign
molecule inputs retain their caller-visible numbering through conversion.

Inspect `found`, `status`, `optimal` and `complete` separately. A time-, node-
or cancellation-stopped call can retain a useful witness without proving an
optimum. An unfinished result with no witness does not prove absence.
`complete` also records enumeration completeness; a requested result limit
can stop enumeration after an optimum is found. Ring-perception work limits
can prevent a proof even when a mapping is feasible. Preserve stopped cases
when analysing outcomes; stopped latency is not completed-work throughput.

Search budgets apply after input parsing and preparation. General API options
permit native defaults and unlimited node settings; set explicit finite limits
when accepting untrusted or interactive requests. The MCP tools impose their
own stricter limits. `timeout_ms=0` requests an immediate timeout. A
`CancellationToken` is one-way: call `cancel()` and create a fresh token for a
later independent request.

`ChemOptions` controls atom type, charge, isotope, bond order, aromatic
compatibility, stereo, ring and optional tautomer handling. The default
profile does not require matching formal charges or isotopes. Use explicit
options or `ChemOptions.profile("strict")` when those restrictions matter;
do not infer them from a tool or method name. `complete_rings_only`,
`ring_matches_ring_only` and `ring_fusion_mode` impose different restrictions.
The graph aromaticity model is a preparation choice, selected through
`perceive_aromaticity` or the builder. The deprecated chemistry
`aromaticity_model` field is not read by search.

```python
query = smsd.Substructure.smarts("[OX2H]")
assert query.match("CCO").found
assert not query.match("CCN").found

embeddings = smsd.Substructure.smarts("[#8]").match_all("OCCO", limit=8)
assert len(embeddings.matches) == 2
```

`Substructure.find` takes a molecule query; `Substructure.smarts` compiles a
predicate query. `with_induced(True)` additionally requires query non-bonds
to remain non-bonds. SMARTS predicates can specify hydrogen, isotope, charge,
aromaticity and topology directly. A parsed SMARTS string and a proved search
result are separate facts.

`MCS.match_constrained` processes queries in order while excluding already
used target atoms. `MCS.consensus` builds a progressive common graph with a
support threshold and supporting mappings. It is an approximation and always
reports `optimal=False`, including empty and single-input calls. `complete`
reports whether the progressive pass finished. Neither it nor `COMPLETED`
proves a global multi-molecule optimum or absence of a qualifying common graph.
`is_mapping_maximal` tests extendibility of a given valid mapping, which is different from proving a maximum. Symmetry
canonicalisation and tautomer-confidence helpers also have separate contracts;
confidence is a heuristic, not a search proof or stability prediction.

## Result graphs, SMARTS and images

`Match.to_graph(query)` and `Match.to_smiles(query)` extract selected atoms
and bonds. A partial selected graph can have capped valences and should not
be mistaken for an unchanged fragment of the original complete molecule.
`Match.to_smarts(query, chemistry)` and `MCS.match_with_smarts(target)` export
portable predicates for represented selected graphs. Unsupported stereo,
tautomer and ring-export constraints raise rather than being silently dropped.
A SMARTS predicate does not carry an optimum proof, enumeration status or all
induced non-edge requirements; retain the associated `SearchResult`.

```python
exported = smsd.MCS.find("CCO").match_with_smarts("CCN")
assert exported.result.found
assert smsd.Substructure.smarts(exported.smarts).match("CCN").found

from smsd.highlight import mcs_svg
svg = mcs_svg("CCO", "CCN", exported.result)
assert '<svg' in svg and 'data-search-status=' in svg
```

Native `depict_svg`, `depict_mapping` and `depict_pair` return SVG strings and
use `DepictOptions`. The pair renderer highlights mapped atoms and input
bonds between them. `smsd.highlight.mcs_svg` carries search flags in the SVG
and rejects selected-edge witnesses that this renderer would over-highlight;
render `best.to_graph(query)` separately in that case. The MCP depiction tool
instead labels this highlighting scope and returns the selected bond mapping
alongside the image.

The older `depict_mcs`, `depict_substructure` and deprecated `to_svg` helpers
use the optional RDKit drawing dependency. Their convenience images do not
replace the proof-bearing result. The native SVG functions need no external
drawing toolkit. `save_svg` writes only when explicitly called.

## Molecule I/O and chemistry

`parse_smiles` accepts OpenSMILES text; an empty string constructs an empty
graph. `lenient=True` explicitly requests recovery from malformed input.
`to_smiles` writes the native canonical representation, including supported
stereo and atom maps. Canonical ordering can differ between language ports.
`canonical_hash` and graph canonical hashes are compact screening keys, not
proofs of molecular identity or difference. `same_canonical_graph` compares
native canonical outputs; use a validated mapping when exact identity is
required by an application.

`read_mol_block`, `read_mol_file`, `write_mol_block`,
`write_mol_block_v3000` and `write_molfile` cover MOL V2000/V3000 text or
files. `read_sdf`, `write_sdf` and `write_sdf_record` cover SDF records.
`read_sdf` loads the whole file and represents malformed records as empty
graphs; inspect them rather than counting every record as accepted chemistry.

```python
from pathlib import Path
from tempfile import TemporaryDirectory

with TemporaryDirectory() as directory:
    destination = Path(directory) / "molecules.sdf"
    smsd.write_sdf(["CCO", "CCN"], destination)
    assert [graph.n for graph in smsd.read_sdf(destination)] == [3, 3]
```

`MolGraphBuilder` takes explicit atom properties, adjacency and parallel bond
properties. Its `build` can run aromaticity perception; select that behaviour
deliberately when carrying an externally prepared representation. `MolGraph`
exposes atom and bond properties and lazy invariant getters. `prewarm` prepares
common invariants; `perceive_aromaticity` and `kekulize` explicitly change a
graph. Do not mutate it concurrently with reads/searches. Molecule names,
comments, program lines and SDF properties are ordinary metadata.

`compute_sssr`, `rings_truncated`, `relevant_cycles_truncated` and graph
ring-system methods expose ring perception. `rings_truncated` keeps its 8.0.x
meaning: the ring search behind `compute_sssr` stopped before finding a
complete set of rings. That search works within one fixed, deterministic work
limit; see [Results and budgets](ALGORITHMS.md#results-and-budgets).
`relevant_cycles_truncated` separately reports incomplete relevant-cycle
enumeration. Strict ring fusion compares per-atom ring counts; complete-ring
matching uses relevant cycles. Each checks completeness of
the perception it uses. `assign_cip`, `assign_rs` and
`assign_ez` expose represented stereo descriptors, including unresolved cases;
missing annotations are not inferred experimental configurations.
`count_components`, `split_components`, `extract_subgraph` and
`murcko_scaffold` return graph topology or transformed graph objects, rather
than implicitly standardising every input.

`compute_relevant_cycle_families` groups the relevant cycles into families,
each with one member and its exact number of members as a Python integer, so
they can be counted without listing them; `truncated` reports an incomplete
enumeration. `relevant_cycle_family_of` finds the family of a ring and
`relevant_cycle_family_members` lists members up to a limit. When the
enumeration is complete, the relevant cycles and their number do not depend on
the input atom order; the grouping into families can.

`from_rdkit`/`to_rdkit` and `from_openbabel`/`to_openbabel` use optional toolkit
bindings. They pass through molecular text and are not unrestricted object
serialisation. `get_index_map` and `translate_mapping` preserve RDKit atom
numbering. Cache clearing releases reusable conversions without invalidating
index maps retained by live graphs. `translate_to_atom_ids` uses explicitly
supplied external graph atom IDs. Keep toolkit versions and conversion choices
with cross-tool results.

## Fingerprints and coefficients

`fingerprint` and `fingerprints` accept molecules or SMILES and return sorted
set-bit positions, or non-zero `(position, count)` pairs when counts are
requested. `fingerprint_from_smiles` is a text convenience helper. Families
include MCS-aware paths, element paths, ECFP/circular, FCFP and topological
torsions; `FingerprintOptions` adds configurable invariants, fragment/scaffold
choices, stereo and SMARTS keys. `FingerprintOptions.maccs()` supplies 166
MACCS keys, with key number `k` at position `k - 1`. Fingerprint family names
and options are not guarantees of cross-tool bit identity.

```python
specification = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
first, second = smsd.fingerprints(["CCO", "CCN"], options=specification, num_threads=1)
assert 0 <= smsd.tanimoto_coefficient(first, second) <= 1

counted = smsd.FingerprintOptions(family="ecfp", size=128, counts=True)
counts = smsd.fingerprint("CCCCCC", options=counted)
assert any(count > 1 for position, count in counts)
```

Bit Tanimoto is intersection over union; bit Dice and cosine use bit counts.
Count Tanimoto uses `sum(min)/sum(max)`, count Dice uses
`2*sum(min)/(sum(a)+sum(b))`, and count cosine uses the dot product and L2
norms. These are different definitions. Use the count-specific helpers when
retaining multiplicity. Empty-input conventions are documented on individual
functions. Pass the actual fingerprint width to bit Dice/cosine/distance
helpers whose default is 2,048 bits. `to_hex`, `from_hex`,
`to_binary_string` and `counts_to_array` provide explicit serialisation forms.
Fingerprint quality estimates concern hashed-vector density/collisions, not
biological accuracy.

`batch_similarity` ranks by ECFP Tanimoto with stable input-order ties.
`similarity_upper_bound`, its `similarity` alias and `screen_targets` are
graph-label bounds used for screening; they are neither the same coefficient
nor an exact MCS search. `TargetCorpus` captures molecules once for repeated
graph-bound or ECFP screens and exact prepared-query batches.

## Corpus search decisions

`Substructure.match_corpus(corpus, *, gpu_prefilter=False, num_threads=0)`
decides for each target whether it contains the query, and
`MCS.match_corpus(corpus, minimum_atoms, *, gpu_prefilter=False, num_threads=0)`
whether its unweighted MCIS with the query has at least `minimum_atoms` atoms.
A conservative pre-screen rules out targets that cannot qualify; the CPU
searches the rest, including targets whose bound equals the threshold.
Existing `match_many` is unchanged.

Each `CorpusSearchResult` has `target_index`, `atom_upper_bound` (integer or
None), `proven_below_threshold`, `threshold_reached` and `result` (the actual
`SearchResult`, or None for a target ruled out by the pre-screen). Both
decision flags false means unknown; a search stopped below the threshold
proves nothing. Substructure has no atom bound. MCES and weighted searches
return CPU results, but their scores never prove that a target falls below the
atom threshold.

`gpu_prefilter=False` runs the pre-screen on the CPU. `True` asks the optional
GPU provider to run it and falls back to the CPU when the device cannot. The
GPU is never selected automatically, and graph matching always runs on the
CPU. With unsupported query options, node caps, immediate timeouts or
cancellation tokens, every target is searched. Original target order and atom
numbering are retained. See the
[Python examples and budget scope](PYTHON.md#explicit-corpus-prefilters).

## Packed corpora, clustering and optional data adapters

`FingerprintCorpus` is a read-only memory-mapped packed-bit file with explicit
fingerprint width and row count. `create`/`from_smiles` construct it;
`fingerprint`, `iter_packed_chunks`, `screen_many`, `similarities`,
`threshold_pairs` and `top_k_pairs` read it. Use its context manager or `close`.
The underlying file must remain unchanged while open. Chunk bounds limit
working arrays, not returned-hit storage, mapped pages or total comparison
work. Counts of representable rows are not throughput or memory measurements.

`stream_screen_smiles` and `stream_screen_many_smiles` process target
iterables in chunks. `butina` clusters by similarity to a chosen centre;
members need not be mutually close. `maxmin` makes a greedy diverse selection
with deterministic input-index ties. Full threshold pair/clustering work can
be quadratic even when it avoids a dense distance matrix.

The `smsd.adapters` module imports optional dependencies only when requested:

| API | Output and dependency |
| --- | --- |
| `to_numpy`, `fingerprints_numpy` | Dense binary/count arrays; NumPy. |
| `to_scipy`, `fingerprints_scipy` | Sparse matrices; SciPy and NumPy. |
| `molecules_pandas`, `molecules_polars` | Molecule/fingerprint columns; pandas or Polars. |
| `to_arrow` | Fixed-size-list fingerprint table; PyArrow and NumPy. |
| `to_torch`, `fingerprints_torch` | Tensors; PyTorch and NumPy. |
| `to_pyg` | Molecular graph `Data`; PyTorch Geometric and PyTorch. |
| `FingerprintTransformer`, `SimilarityTransformer` | Fit/transform estimators for scikit-learn pipelines. |

Arrow can share a native-endian, C-contiguous 2-D integer NumPy buffer and
retains its owner. Shared input must not be mutated. Other inputs/layouts can
require copies; `zero_copy_only=True` rejects those cases and table metadata
records the ownership/copy result. A tensor/graph adapter's device argument
does not itself establish a GPU chemistry kernel.

## Accelerators and molecular tools

Exact MCS/substructure proof remains on the CPU. GPU support is operation-,
build-, device- and package-dependent; availability alone does not prove a
particular call dispatched to a device. `gpu_is_available` and
`gpu_device_info` describe discovery. `warm_start` explicitly prepares a
reusable corpus/query path. Packed corpus APIs expose explicit backends and
retain deterministic CPU-defined score/hit semantics. Native Metal packed
fingerprint screening is available only in C++ and runs only when requested
explicitly; Python screening and search never run on Metal. See
[CUDA operation scope](PYTHON.md#optional-cuda-screening),
[native Metal scope](CPP.md#optional-metal-fingerprint-screening) and the optional
[`smsd-cuda` package](../packaging/smsd-cuda/README.md) for build and device
requirements; the [supported-platforms table](HOWTO-INSTALL.md#supported-platforms)
lists the GPU that automatic screening uses by default.

The eight [MCP tools](MCP_SERVER.md) expose finite text, atom, target, result
and search limits plus status/proof/timing metadata. Native default
standardisation changes representation only; opt-in chemical transformations
record their optional backend. The tool registry also exports local function
definitions in the `input-schema`, `function` and `chat-function` formats, and
LangChain and LlamaIndex tools, without creating a model client.

## Compatibility names

The 8.x API retains deprecated names such as `canonical_smiles`,
`prewarm_graph`, `read_molfile`, `export_sdf`, `dearomatize`,
`overlapCoefficient`, `count_overlapCoefficient`, `to_svg` and `chem=`.
Use their documented replacements. Retaining these names is not permission
to remove or reinterpret them in a minor release; warnings and old return
values remain part of the compatibility behaviour.

## Native graph and option objects

| Object | Public operations |
| --- | --- |
| `MolGraph` | `bonds()`, `has_bond(i, j)`, `bond_order(i, j)`, `bond_in_ring(i, j)`, `bond_aromatic(i, j)` read represented topology. `prewarm()`, `ensure_canonical()` and `ensure_ring_counts()` prepare lazy invariants. `get_canonical_hash()`, `get_automorphism_generators()` and `automorphism_generators_truncated()` expose bounded canonical/symmetry results. `num_rings()`, `ring_system_of(atom)`, `ring_system_sig(rs_id)` and `ring_system_count()` describe perceived rings. `perceive_aromaticity(model)`, `kekulize()` and `dearomatize()` explicitly mutate representation. |
| `MolGraphBuilder` | `atom_count`, `atomic_numbers`, `formal_charges`, `mass_numbers`, `hydrogen_counts`, `ring_flags`, `aromatic_flags`, `neighbors`, `bond_orders`, `bond_ring_flags`, `bond_aromatic_flags`, `atom_ids`, `name`, `program_line`, `comment` and `properties` set builder data and return the builder. `build(perceive_aromaticity=True, aromaticity_model=DAYLIGHT_LIKE)` constructs the graph. Bond-property rows run parallel to neighbour rows. |
| `ChemOptions` | `profile(name)` selects a named profile and `tautomer_profile()` selects the tautomer-aware profile. Fluent setters `with_match_isotope`, `with_complete_rings_only`, `with_tautomer_aware`, `with_solvent`, `with_bond_order_mode`, `with_aromaticity_mode`, `with_ph` and `with_ring_fusion_mode` change and return that options object. Constructing `ChemOptions(existing)` copies it. |
| `MCS.Options` | Search chemistry, `timeout_ms`, `node_limit`, `cancellation`, `objective`, `connectivity`, `min_fragment_atoms`, `max_fragments`, `atom_weights` and `excluded_target_atoms`. |
| `Substructure.Options` | Search chemistry, `timeout_ms`, `node_limit`, `cancellation` and `induced`. |
| `CancellationToken` | `cancel()` sets the one-way signal; `is_cancelled` reads it. |
| `FingerprintOptions` | Constructor keyword options and mutable fields describe fingerprint terms. `copy()` makes an independent copy; factories `ecfp(radius=2, size=2048)`, `fcfp(radius=2, size=2048)`, `path(path_length=7, size=2048)`, `torsion(size=2048)` and `maccs()` select classic configurations. `maccs_keys()` returns the key definitions. |
| `FingerprintKey` | `FingerprintKey(smarts, minimum=1)` defines a SMARTS occurrence key; read-only `smarts`, `minimum` and `kind` describe it. |
| `CIPDescriptors` | `rs_labels` is a label per atom; `ez_bonds` contains `(atom1, atom2, label)` entries. Produced by `assign_cip`. |
| `RelevantCycleFamilies` | Read-only `families`, `relevant_cycle_count` (the exact number of members of the families listed) and `truncated`. Produced by `compute_relevant_cycle_families`. |
| `RelevantCycleFamily` | Read-only `size`, `root`, `first`, `second`, `opposite`, `prototype` (one member as atom indices) and the exact `member_count`. |
| `Point2D`, `Point3D` | Empty or coordinate constructors, with writable `x`, `y` and, for 3-D, `z`. |
| `DepictOptions` | Native SVG dimensions, bond lengths/widths/spacing, atom labels, map/index visibility, highlight styling, padding and font family. |

Native enums are `SearchStatus`, `BondOrderMode`, `AromaticityMode`,
`AromaticityModel`, `RingFusionMode`, `Solvent`, `FingerprintFamily`,
`FragmentMode`, `ScaffoldMode`, `RSLabel`, `EZLabel`, and the nested
`MCS.Objective`/`MCS.Connectivity`. Use their named values rather than numeric
ordinals. Graph atom arrays are in graph index order; `canonical_label` and
`orbit` are derived data and are not permission to replace caller numbering.
`properties` supports ordinary string-key/value mapping operations.

`smsd.RingBasisUnavailable` is the native `RuntimeError` subclass raised when
a result needs the minimum cycle basis and ring perception reached its work
limit before completing it; batch calls add the zero-based input index.

## Callable index

The source-linked index below includes the top-level purpose modules, their
public prepared-search/corpus methods, the separate adapters/highlight modules
and molecular tool functions. Native graph/options fields and builder methods
are grouped above and in their runtime docstrings. Optional dependencies are
required only for the corresponding call.

<!-- CALLABLE_INDEX -->

### Molecules and I/O

Source: [_molecule.py](../python/smsd/_molecule.py).

| Callable | Contract |
| --- | --- |
| `smsd.canonical_hash(mol)` | The canonical hash of a molecule as a Python int. |
| `smsd.canonical_smiles(mol)` | Deprecated alias of :func:to_smiles. |
| `smsd.count_components(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.export_sdf(molecules, filename)` | Deprecated alias of :func:write_sdf. |
| `smsd.clear_cache()` | Clear reusable toolkit conversions; live graphs keep their atom correspondence. |
| `smsd.extract_subgraph(mol, atom_indices)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.murcko_scaffold(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.parse_mapped_smiles(smiles)` | Parse a SMILES string that contains atom map numbers. |
| `smsd.parse_smiles(smiles, *, lenient=False)` | Parse a SMILES string into a MolGraph. |
| `smsd.prewarm(mol)` | Compute and cache the common lazy invariants of a molecule once; returns its graph. |
| `smsd.prewarm_graph(mol)` | Deprecated alias of :func:prewarm that returns None. |
| `smsd.read_mol_block(mol_block)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.read_mol_file(filename)` | Read a MOL file (V2000 or V3000) into a MolGraph; one that cannot be opened or read raises ValueError. |
| `smsd.read_molfile(filename)` | Deprecated alias of :func:read_mol_file that reads the file as UTF-8 text. |
| `smsd.read_sdf(filename)` | Read every record of an SDF file at once; malformed records become empty graphs. |
| `smsd.same_canonical_graph(mol1, mol2)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.split_components(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.strip_atom_maps(smiles)` | Remove atom map numbers from a SMILES string. |
| `smsd.to_smarts(mol, include_aromaticity=True, include_charge=True, include_ring_member=False, include_h_count=False, include_isotope=False)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.to_smiles(mol)` | Canonical SMILES of a MolGraph, SMILES string or RDKit Mol. |
| `smsd.translate_to_atom_ids(mapping, g1, g2)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.write_mol_block(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.write_mol_block_v3000(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.write_molfile(mol, filename, *, v3000=False, sdf=False)` | Write a molecule to disk as MOL V2000, MOL V3000 or a single SDF record. |
| `smsd.write_sdf(molecules, filename)` | Write molecules as SDF records, each a V2000 MOL block followed by $$$$. |
| `smsd.write_sdf_record(mol)` | Native callable; its runtime docstring includes an ordinary example. |

### Chemical graph helpers

Source: [_chemistry.py](../python/smsd/_chemistry.py).

| Callable | Contract |
| --- | --- |
| `smsd.CIPDescriptors (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.EZLabel (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.RSLabel (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.RelevantCycleFamilies (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.RelevantCycleFamily (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.assign_cip(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.assign_ez(mol)` | Assign CIP E/Z descriptors to all stereogenic double bonds. |
| `smsd.assign_rs(mol)` | Assign CIP R/S descriptors to all tetrahedral stereocentres. |
| `smsd.classify_pharmacophore(mol, atom_index)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.compute_relevant_cycle_families(mol)` | The relevant cycles grouped into families, each with one member and its exact number of members. |
| `smsd.compute_sssr(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.dearomatize(mol)` | Deprecated alias of :func:kekulize. |
| `smsd.implicit_h(atomic_num, bond_order_sum, formal_charge, aromatic=False)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.kekulize(mol)` | Convert an aromatic MolGraph into an explicit Kekule form. |
| `smsd.perceive_aromaticity(mol, *, model=_smsd.AromaticityModel.DAYLIGHT_LIKE)` | Recompute ring membership and aromaticity on a MolGraph. |
| `smsd.rings_truncated(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.relevant_cycles_truncated(mol)` | Whether full relevant-cycle enumeration is incomplete, independently of the minimum-basis flag. |
| `smsd.relevant_cycle_family_members(mol, family, limit)` | Up to `limit` members of a family as atom-index lists, starting with its prototype. |
| `smsd.relevant_cycle_family_of(mol, families, ring)` | The index in `families.families` of the family that has `ring` as a member, or -1. |

### Prepared searches and results

Source: [search.py](../python/smsd/search.py).

| Callable | Contract |
| --- | --- |
| `smsd.ConsensusResult (class)` | A progressive common graph and supporting embeddings; `optimal` is always false and `complete` records progressive-pass completion. |
| `smsd.ConsensusResult.to_smiles()` | The common graph as SMILES. |
| `smsd.MCS (class)` | Find common induced or edge subgraphs, including disconnected fragments. |
| `smsd.MCS.match_with_smarts(target)` | Return the match and its SMARTS under this query's chemistry. |
| `smsd.MCS.find(query)` | Prepare a molecule query for repeated and concurrent searches. |
| `smsd.MCS.with_objective(objective)` | This query maximising objective: common atoms (MCIS) or common bonds (MCES). |
| `smsd.MCS.with_fragments(*, min_atoms=1, max_fragments=2147483647)` | This query allowing disconnected common fragments of at least min_atoms atoms each. |
| `smsd.MCS.with_atom_weights(weights)` | This query scoring matched query atoms by weights, one per input atom. |
| `smsd.MCS.validate_mapping(query, target, match, options=None)` | Return selected-graph validation errors; an empty list is valid. |
| `smsd.MCS.is_mapping_maximal(query, target, match, options=None)` | Check whether a valid common graph admits a compatible extension. |
| `smsd.MCS.validate_tautomer_consistency(query, target, match)` | Whether a mapping conserves mobile protons across tautomeric groups. |
| `smsd.MCS.canonicalize_mapping(query, target, atom_mapping)` | Choose a representative under the two molecules' symmetries. |
| `smsd.MCS.compute_tautomer_confidence(query, target, match, chemistry)` | Score selected chemistry relaxations using a read-time pH/solvent heuristic. |
| `smsd.MCS.match_constrained(queries, target, options=None)` | Match queries in order without reusing original target atoms. |
| `smsd.MCS.consensus(molecules, options=None, threshold=1.0)` | Find an approximate progressive common graph supported by a fraction of inputs; never a proof of optimality or absence. |
| `smsd.Match (class)` | Atom and selected-bond correspondences, with common-edge fragments. |
| `smsd.Match.atom_count (property)` | The number of matched atoms. |
| `smsd.Match.bond_count (property)` | The number of matched bonds. |
| `smsd.Match.to_smiles(query)` | Write the selected query atoms and bonds as a capped molecule. |
| `smsd.Match.to_graph(query)` | Extract selected edges while retaining original matching context. |
| `smsd.Match.to_smarts(query, chemistry=None)` | Export selected bonds as SMARTS; reject rules SMARTS cannot represent. |
| `smsd.SearchResult (class)` | Matches and the proof, enumeration and stopping status of a search. |
| `smsd.SearchResult.found (property)` | Whether the search found at least one match. |
| `smsd.SearchResult.best (property)` | The first match, or None when nothing was found. |
| `smsd.SmartsResult (class)` | A search result with a portable SMARTS query for its best selected graph. |
| `smsd.Substructure (class)` | Find molecule or SMARTS embeddings using one public query interface. |
| `smsd.Substructure.find(query)` | Prepare a molecule query; SMILES strings are parsed as molecules. |
| `smsd.Substructure.smarts(query)` | Compile a SMARTS predicate query once for repeated searches. |
| `smsd.Substructure.with_induced(enabled=True)` | This query requiring an induced embedding, so query non-bonds are kept apart too. |
| `smsd.Substructure.decompose_r_groups(core, targets, options=None)` | Return core and substituent graphs alongside each search result. |
| `smsd.MCS.with_options(options)` | This query with a copy of options. |
| `smsd.Substructure.with_options(options)` | This query with a copy of options. |
| `smsd.MCS.with_timeout(milliseconds)` | This query with a time limit per target, in milliseconds; a negative value selects the default limit and zero stops at once. |
| `smsd.Substructure.with_timeout(milliseconds)` | This query with a time limit per target, in milliseconds; a negative value selects the default limit and zero stops at once. |
| `smsd.MCS.match(target)` | Return one best match, with its proof and stopping status. |
| `smsd.Substructure.match(target)` | Return one best match, with its proof and stopping status. |
| `smsd.MCS.match_all(target, limit=1000)` | Return optimal matches up to the limit, with enumeration status. |
| `smsd.Substructure.match_all(target, limit=1000)` | Return optimal matches up to the limit, with enumeration status. |
| `smsd.CorpusSearchResult (class)` | Original corpus slot, bound decision and optional actual CPU result. |
| `smsd.MCS.match_corpus(corpus, minimum_atoms, *, gpu_prefilter=False, num_threads=0)` | Atom-count decisions for unweighted MCIS, with unsupported contracts kept on the CPU. |
| `smsd.Substructure.match_corpus(corpus, *, gpu_prefilter=False, num_threads=0)` | Containment decisions per target after an optional pre-screen, with ordinary CPU embeddings. |
| `smsd.MCS.match_many(targets, *, num_threads=0)` | Search targets independently and preserve each target's atom indices. |
| `smsd.Substructure.match_many(targets, *, num_threads=0)` | Search targets independently and preserve each target's atom indices. |

### Fingerprints and metrics

Source: [_fingerprints.py](../python/smsd/_fingerprints.py).

| Callable | Contract |
| --- | --- |
| `smsd.analyze_fp_quality(fp, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.circular_fingerprint(mol, radius=2, fp_size=2048, mode="ecfp")` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.circular_fingerprint_counts(mol, radius=2, fp_size=2048, mode="ecfp")` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.cosine(fp1, fp2, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.count_cosine(fp1, fp2)` | Count-vector cosine: dot(a,b) / (\|a\| * \|b\|). |
| `smsd.count_dice(fp1, fp2)` | Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)). |
| `smsd.count_overlapCoefficient(fp1, fp2)` | Deprecated alias of :func:count_overlap_coefficient. |
| `smsd.count_overlap_coefficient(fp1, fp2)` | Count-vector Szymkiewicz-Simpson overlap: sum(min(a,b)) / min(sum(a), sum(b)). |
| `smsd.count_tanimoto_coefficient(fp1, fp2)` | Count-vector Tanimoto (generalised Jaccard): sum(min(a,b)) / sum(max(a,b)). |
| `smsd.counts_to_array(counts_dict, fp_size)` | Convert a sparse count dict {bit_position: count} into a dense list of length fp_size. |
| `smsd.dice(fp1, fp2, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.fingerprint(mol, *, kind='mcs', path_length=7, fp_size=2048, radius=2, options=None)` | Compute a molecular fingerprint. |
| `smsd.fingerprint_from_smiles(smiles, radius=2, fp_size=2048, mode='ecfp', *, path_length=7, options=None)` | Parse a SMILES string and compute its fingerprint in one call. |
| `smsd.fingerprint_subset(fp1, fp2)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.fingerprints(mols, *, kind='mcs', path_length=7, fp_size=2048, radius=2, options=None, num_threads=0)` | Compute the fingerprints of many molecules in one call. |
| `smsd.from_hex(hex_str)` | Convert a hex string back to a fingerprint (list of set-bit positions). |
| `smsd.mcs_fingerprint(mol, path_length=7, fp_size=2048, tautomer_aware=False)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.overlapCoefficient(fp1, fp2, fp_size=2048)` | Deprecated alias of :func:overlap_coefficient. |
| `smsd.overlap_coefficient(fp1, fp2, fp_size=2048)` | Overlap coefficient (Szymkiewicz-Simpson): \|A∩B\| / min(\|A\|, \|B\|). |
| `smsd.path_fingerprint(mol, path_length=7, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.soergel(fp1, fp2, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.tanimoto_coefficient(fp1, fp2, fp_size=2048)` | Tanimoto coefficient (Jaccard index): \|A∩B\| / \|A∪B\|. |
| `smsd.to_binary_string(fp, fp_size=2048)` | Convert a fingerprint to a binary string of '0' and '1' characters. |
| `smsd.to_hex(fp, fp_size=2048)` | Convert a fingerprint (list of set-bit positions) to a hex string. |
| `smsd.topological_torsion(mol, fp_size=2048)` | Compute a topological torsion fingerprint (Nilakantan et al. 1987). |
| `smsd.topological_torsion_counts(mol, fp_size=2048)` | Native callable; its runtime docstring includes an ordinary example. |

### Graph and fingerprint screening

Source: [_screening.py](../python/smsd/_screening.py).

| Callable | Contract |
| --- | --- |
| `smsd.TargetCorpus (class)` | Targets captured once, with their input atom numbering, for repeated searches and screening. |
| `smsd.TargetCorpus.from_smiles(targets)` | Create a corpus from molecule objects or SMILES strings. |
| `smsd.TargetCorpus.add_target(target)` | Capture a molecule and its current input atom numbering. |
| `smsd.TargetCorpus.add_targets(targets)` | Capture targets together; later input changes do not alter them. |
| `smsd.TargetCorpus.prewarm(*, num_threads=0)` | Compute graph invariants and ECFP fingerprints once; later searches read the prewarmed graphs. |
| `smsd.TargetCorpus.is_prewarmed (property)` | Whether :meth:prewarm has run since the last targets were added. |
| `smsd.TargetCorpus.screen(query, threshold, *, num_threads=0)` | Return target indices meeting the ECFP Tanimoto threshold. |
| `smsd.TargetCorpus.screen_many(queries, threshold, *, num_threads=0)` | Return ECFP Tanimoto hit indices for each query in input order. |
| `smsd.TargetCorpus.screen_targets_many(queries, threshold, chemistry=None, *, num_threads=0)` | Return similarity-bound hit indices for each query in input order. |
| `smsd.batch_similarity(query, targets, fp_size=2048, radius=2)` | Rank targets by ECFP Tanimoto similarity to a query. |
| `smsd.gpu_device_info()` | Human-readable compute backend string, including a CPU-only backend. |
| `smsd.gpu_is_available()` | True when a GPU is available at run time: through an installed smsd-cuda package that finds a device, or through a native module built with CUDA or Metal. Python screening and search never run on Metal; a Metal-enabled source build changes only what this function and `gpu_device_info()` report. |
| `smsd.screen_targets(query, targets, threshold, chemistry=None, *, chem=None, num_threads=0)` | Indices of the targets whose similarity bound to the query reaches threshold. |
| `smsd.similarity(mol1, mol2)` | Upper bound on the similarity of two molecules, in [0.0, 1.0]. |
| `smsd.similarity_upper_bound(g1, g2, chemistry=None, *, chem=None)` | Linear-time upper bound on the similarity of two molecules, in [0.0, 1.0]. |

### Coordinates and layout

Source: [_layout.py](../python/smsd/_layout.py).

| Callable | Contract |
| --- | --- |
| `smsd.Point2D (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.Point3D (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.align_2d(coords, reference)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.bounding_box_2d(coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.canonical_orientation(mol, coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.center_2d(coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.count_crossings(mol, coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.force_directed_layout(mol, coords, max_iter=500, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.generate_coords_2d(mol, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.generate_coords_3d(mol, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.is_degenerate_layout(coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.layout_quality(mol, coords, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.layout_sssr(mol)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.match_template(mol, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.mirror_x(coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.mirror_y(coords)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.normalise_bond_length(mol, coords, target=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.reduce_crossings(mol, coords, max_iter=1000)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.resolve_overlaps(coords, threshold=0.3, max_iter=100)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.rotate_2d(coords, angle)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.scale_2d(coords, factor)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.stress_majorisation(mol, coords, max_iter=300, target_bond_length=1.5)` | Native callable; its runtime docstring includes an ordinary example. |
| `smsd.translate_2d(coords, dx, dy)` | Native callable; its runtime docstring includes an ordinary example. |

### Depiction

Source: [_depiction.py](../python/smsd/_depiction.py).

| Callable | Contract |
| --- | --- |
| `smsd.DepictOptions (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.depict_mapping(mol_or_smiles, mapping, opts=None, **kwargs)` | Render a molecule with MCS or substructure atoms highlighted and numbered. |
| `smsd.depict_mcs(mol1, mol2, *, tautomer_aware=False, timeout_ms=10000, size=(600, 300), highlight_color=(0.5, 0.8, 1.0))` | Find the MCS and return an RDKit image with matched atoms highlighted; `highlight_color` is deprecated and has no effect. |
| `smsd.depict_pair(mol1, mol2, mapping, opts=None, **kwargs)` | Render two molecules side by side with an MCS mapping highlighted. |
| `smsd.depict_substructure(query, target, *, timeout_ms=10000, size=(600, 300))` | Find a substructure and return an RDKit image with the matched atoms highlighted. |
| `smsd.depict_svg(mol_or_smiles, opts=None, **kwargs)` | Render a molecule as an SVG string. |
| `smsd.save_svg(svg_str, path)` | Write SVG text to a file as UTF-8. |
| `smsd.to_svg(mol_or_smiles, *, size=(300, 200))` | Deprecated RDKit rendering of a molecule as SVG. |

### Proof-bearing highlights

Source: [highlight.py](../python/smsd/highlight.py).

| Callable | Contract |
| --- | --- |
| `smsd.highlight.mcs_svg(query, target, result, *, opts=None, **kwargs)` | Render the two molecules with the best result's atom correspondence. |

### Packed corpus and streaming

Source: [corpus.py](../python/smsd/corpus.py).

| Callable | Contract |
| --- | --- |
| `smsd.FingerprintCorpus (class)` | A read-only packed fingerprint file, indexed in input order. |
| `smsd.FingerprintCorpus.create(path, fingerprints, *, fp_size=2048, chunk_size=4096, overwrite=False)` | Stream binary fingerprints to a new file and open it read-only. |
| `smsd.FingerprintCorpus.from_smiles(path, molecules, *, kind='ecfp', fp_size=2048, radius=2, path_length=7, options=None, chunk_size=4096, num_threads=0, overwrite=False)` | Fingerprint molecule inputs in bounded batches, then write packed rows. |
| `smsd.FingerprintCorpus.fingerprint(index)` | Return the sorted set-bit positions of one row, as an owned list. |
| `smsd.FingerprintCorpus.iter_packed_chunks(*, chunk_size=4096)` | Yield (start_index, uint8_array) copies in input order. |
| `smsd.FingerprintCorpus.cuda_calls (property)` | Number of successfully completed packed CUDA calls made by this corpus. |
| `smsd.FingerprintCorpus.screen_many(queries, *, threshold=0.0, chunk_size=4096, backend='cpu')` | Yield (query_index, target_index, Tanimoto) from at most 32 queries. |
| `smsd.FingerprintCorpus.similarities(query, *, threshold=0.0, chunk_size=4096, backend='cpu')` | Yield (target_index, Tanimoto) hits in input order for a bit set. |
| `smsd.FingerprintCorpus.threshold_pairs(threshold, *, include_self=False, chunk_size=4096, backend='cpu')` | Yield undirected (first, second, Tanimoto) pairs meeting the threshold. |
| `smsd.FingerprintCorpus.top_k_pairs(k, *, threshold=0.0, include_self=False, chunk_size=4096, backend='cpu')` | Yield at most k directed neighbours per row, ranked by similarity. |
| `smsd.FingerprintCorpus.close()` | Wait for an active mapped copy, then release; repeated closes are harmless. |
| `smsd.stream_screen_smiles(query, targets, threshold, *, chunk_size=4096, num_threads=0)` | Yield input indices passing an ECFP Tanimoto screen, one target chunk at a time. |
| `smsd.stream_screen_many_smiles(queries, targets, threshold, *, target_chunk_size=4096, query_chunk_size=64, max_queries=4096, num_threads=0)` | Yield (query_index, target_index) ECFP screen hits, one bounded chunk at a time. |

### Clustering and selection

Source: [clustering.py](../python/smsd/clustering.py).

| Callable | Contract |
| --- | --- |
| `smsd.butina(corpus, threshold=0.7, *, reordering=False, chunk_size=4096, backend='cpu')` | Cluster by a binary Tanimoto similarity threshold (Butina, 1999). |
| `smsd.maxmin(corpus, count, *, seed=0, chunk_size=4096, backend='cpu')` | Select diverse row indices by greedy MaxMin Tanimoto distance. |

### Optional data adapters

Source: [adapters.py](../python/smsd/adapters.py).

| Callable | Contract |
| --- | --- |
| `smsd.adapters.to_numpy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None)` | Expand fingerprints to a NumPy array of shape (n_molecules, fp_size). |
| `smsd.adapters.to_torch(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None, device='cpu')` | Expand fingerprints to a PyTorch tensor without requiring NumPy. |
| `smsd.adapters.fingerprints_numpy(molecules, *, kind='mcs', path_length=7, fp_size=2048, radius=2, options=None, num_threads=0, dtype=None)` | Compute native fingerprints and expand them to a NumPy feature matrix. |
| `smsd.adapters.fingerprints_torch(molecules, *, kind='mcs', path_length=7, fp_size=2048, radius=2, options=None, num_threads=0, dtype=None, device='cpu')` | Compute native fingerprints and expand them to a PyTorch feature tensor. |
| `smsd.adapters.FingerprintTransformer (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.adapters.SimilarityTransformer (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.adapters.to_scipy(sparse_fingerprints, *, fp_size=2048, counts=False, dtype=None)` | Return a canonical SciPy CSR matrix without expanding dense rows. |
| `smsd.adapters.fingerprints_scipy(molecules, *, kind='ecfp', path_length=7, fp_size=2048, radius=2, options=None, num_threads=0, dtype=None)` | Compute native fingerprints as an integer CSR feature matrix. |
| `smsd.adapters.molecules_pandas(molecules, *, kind='ecfp', fp_size=2048, radius=2, path_length=7, options=None, num_threads=0)` | Return a pandas DataFrame with molecule and fingerprint object columns. |
| `smsd.adapters.molecules_polars(molecules, *, kind='ecfp', fp_size=2048, radius=2, path_length=7, options=None, num_threads=0)` | Return a Polars DataFrame with explicit Object molecules and typed lists. |
| `smsd.adapters.to_arrow(fingerprints, *, fp_size=2048, counts=False, dtype=None, copy=False, zero_copy_only=False)` | Export a fixed-size-list Arrow fingerprint column with ownership metadata. |
| `smsd.adapters.to_pyg(molecule, *, device='cpu')` | Export native graph features as a PyTorch Geometric Data object. |

### Warm start

Source: [_accelerator.py](../python/smsd/_accelerator.py).

| Callable | Contract |
| --- | --- |
| `smsd.warm_start(corpus, query, threshold=0.0, *, num_threads=0)` | Prepare and screen a real target corpus before its next query. |

### Molecular tools

Source: [mcp_server.py](../python/smsd/mcp_server.py).

| Callable | Contract |
| --- | --- |
| `smsd.mcp_server.create_server()` | Create the official SDK stdio server without starting it. |
| `smsd.mcp_server.find_mcs(query_smiles: str, target_smiles: str, objective: str='mcis', connected: bool=True, timeout_ms: int=1000, node_limit: int=100000, match_stereo: bool=False, complete_rings_only: bool=False, ring_matches_ring_only: bool=False)` | Find one common graph and retain native proof and stopping flags. |
| `smsd.mcp_server.find_substructure(query: str, target_smiles: str, query_format: str='smarts', induced: bool=False, timeout_ms: int=1000, node_limit: int=100000, match_stereo: bool=False, complete_rings_only: bool=False, ring_matches_ring_only: bool=False)` | Find one SMARTS or molecule embedding, retaining native stopping status. |
| `smsd.mcp_server.molecular_fingerprint(smiles: str, kind: str='ecfp', fp_size: int=2048, radius: int=2, path_length: int=7, counts: bool=False)` | Compute native set-bit positions or sparse position/count pairs, not an MCS proof. |
| `smsd.mcp_server.search_smarts(query_smarts: str, target_smiles: str, max_matches: int=16, timeout_ms: int=1000, node_limit: int=100000, match_stereo: bool=False, complete_rings_only: bool=False, ring_matches_ring_only: bool=False)` | Enumerate bounded SMARTS embeddings without turning incomplete results into absence. |
| `smsd.mcp_server.molecular_similarity(query_smiles: str, target_smiles: str, metric: str='tanimoto', kind: str='ecfp', fp_size: int=2048, radius: int=2, path_length: int=7, counts: bool=False)` | Compare native fingerprints; the coefficient is not an MCS or biological-accuracy proof. |
| `smsd.mcp_server.screen_molecules(query_smiles: str, targets_smiles: list[str], threshold: float=0.7, max_hits: int=32, offset: int=0, kind: str='ecfp', fp_size: int=2048, radius: int=2, path_length: int=7, counts: bool=False)` | Screen a bounded in-memory list by Tanimoto, with stable ranked ties and explicit pagination. |
| `smsd.mcp_server.standardise_molecule(smiles: str, remove_atom_maps: bool=False, cleanup: bool=False, largest_fragment: bool=False, neutralise: bool=False, canonical_tautomer: bool=False)` | Canonicalise representation; chemical transformations require explicit options and RDKit. |
| `smsd.mcp_server.depict_mcs(query_smiles: str, target_smiles: str, objective: str='mcis', connected: bool=True, timeout_ms: int=1000, node_limit: int=100000, match_stereo: bool=False, complete_rings_only: bool=False, ring_matches_ring_only: bool=False, width: int=800, height: int=360)` | Return a native mapping-highlight SVG and the unchanged search proof/status. |
| `smsd.mcp_server.main()` | Serve over standard input/output, with diagnostics on standard error. |

### Tool definitions

Source: [tool_definitions.py](../python/smsd/tool_definitions.py).

| Callable | Contract |
| --- | --- |
| `smsd.tool_definitions.ToolDefinition (class)` | A registered molecular callable with a fresh JSON schema on each access. |
| `smsd.tool_definitions.ToolDefinition.name (property)` | The public callable name. |
| `smsd.tool_definitions.ToolDefinition.description (property)` | The callable's opening documentation paragraph. |
| `smsd.tool_definitions.ToolDefinition.input_schema (property)` | Object schema with finite single-argument bounds and Python defaults. |
| `smsd.tool_definitions.TOOL_REGISTRY (class/enum)` | See the native graph/options table or the adapter contract above. |
| `smsd.tool_definitions.get_tool_definitions(provider: str='mcp')` | Export definitions in the `mcp`, `input-schema`, `function` or `chat-function` format. |
| `smsd.tool_definitions.call_tool(name: str, arguments: dict)` | Invoke a registered tool with checked JSON arguments and ordinary defaults. |
| `smsd.tool_definitions.langchain_tools()` | Make LangChain StructuredTools with lazy optional imports and the same bounds. |
| `smsd.tool_definitions.llamaindex_tools()` | Make LlamaIndex FunctionTools with lazy optional imports and the same bounds. |

### CUDA package entry points

Source: [__init__.py](../packaging/smsd-cuda/python/smsd_cuda/__init__.py).

| Callable | Contract |
| --- | --- |
| `smsd_cuda.capsule()` | A capsule holding the accelerator's table, which smsd attaches. |
| `smsd_cuda.device_info()` | A description of the device, such as 'Tesla T4 (compute 7.5, 14930 MB)'. |
| `smsd_cuda.is_available()` | Whether a CUDA device can be used. |
| `smsd_cuda.preload()` | Initialise the current CUDA device and load screening kernels without screening data. |
| `smsd_cuda.unavailable_reason()` | Why no CUDA device can be used, or None when one can. |

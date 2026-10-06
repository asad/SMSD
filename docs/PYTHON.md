# SMSD Python Guide

SMSD exposes native C++ molecular graph search through pybind11. Core parsing,
MCS, substructure, fingerprints and SVG depiction work without RDKit. RDKit
is optional for molecule conversion, independent checks and drawing.

The proposed 7.2.1 changes are unreleased. The current PyPI release is
7.1.1; the GitHub release is 7.2.0. See the [local benchmark report](../benchmarks/RESULTS_7.2.0.md) for
versions, settings, measurements and limitations.

## Install and build

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install smsd
# Optional interoperability:
python -m pip install rdkit
```

The package declares Python 3.9 or later; wheel availability depends on Python,
platform and architecture. Release preparation targets CPython 3.14 wheels
for Linux x86_64, macOS arm64 and Windows x86_64, plus a source distribution.
Local macOS/Linux builds and native GitHub Windows checks must validate the
same 7.2.1 source before publication. See [7.2.1 validation](VALIDATION_7.2.1.md)
for pending results. The historical 7.2.0 search comparison runs Python 3.13.14
on macOS arm64 so both versions use the same interpreter and RDKit.

The root `pyproject.toml` is the canonical package manifest. It combines the
native extension in `cpp/` with the Python package in `python/smsd/`; build
from the repository root. Java sources and Maven artifacts are separate under
`java/`.
Source builds enable OpenMP when available. Metal and CUDA detection default
to `AUTO`; the local comparison uses CPU-only builds explicitly:

```bash
python -m pip install build scikit-build-core pybind11
python -m build --wheel --no-isolation \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Check the installed backend with `smsd.gpu_device_info()`. GPU screening is
optional; MCS search does not become a GPU solver merely because a backend is
available.

## Search and indices

```python
import smsd

mapping = smsd.find_mcs("c1ccc(O)cc1", "c1ccc(N)cc1", timeout_ms=1000)
assert len(mapping) == 6
assert smsd.is_substructure("c1ccccc1", "c1ccc(O)cc1")
embeddings = smsd.find_substructure("CO", "CCO", max_results=10)
```

`find_mcs` and `find_substructure` return a dictionary for one result and a
list of dictionaries for `max_results > 1`. Keys are query atom indices and
values are target atom indices. Indices are zero-based: SMILES inputs use
parser order, `MolGraph` inputs use graph order, and RDKit inputs use the
original RDKit molecule order. Mixed input types follow the same rule for
each side. Explicit hydrogen vertices retained by an RDKit molecule are retained by
the conversion. Bracket hydrogen counts do not create separate vertices.

Mappings are injective. SMSD preserves every query edge whose endpoints are
mapped, with the requested chemistry constraints. `induced=True` also rejects
extra target edges. RDKit FMCS can omit query edges from its selected common
subgraph, so equal atom counts do not necessarily mean equivalent results.

The default `strategy="auto"` uses the native coverage pipeline for supported
settings and the full native solver for advanced options. `strategy="native"`
selects the full solver explicitly. `strategy="lightweight"` supports a single
connected atom objective, timeout, ring matching and strict/any bond order;
unsupported options raise `ValueError`. Its separate historical wrapper,
`find_mcs_lightweight`, returns **one-based** pairs in `LightMCSResult.mapping`.

A time budget produces the best valid mapping found by the pipeline. The
public mapping API has no cancellation flag or optimality certificate.
Larger searches use heuristics; reaching the budget does not prove maximum
size. A false/empty substructure result can mean no witness was found within
the budget, rather than a completed proof of absence. Deadlines are cooperative: graph preparation and a work unit may
overshoot the requested wall-clock limit.

## Chemistry and search options

```python
mapping = smsd.find_mcs(
    "CCC", "C1CC1", strategy="native", timeout_ms=1000,
    induced=True, connected_only=True, match_bond_order="strict",
)
assert len(mapping) == 2

weighted = smsd.find_mcs(
    "CCC", "CCC", strategy="native", atom_weights=[10.0, -30.0, 1.0],
)
assert set(weighted) == {0}
```

| Keyword | Default | Meaning |
|---|---|---|
| `match_atom_type` | `True` | Match elements; `False` explicitly relaxes this |
| `match_formal_charge` | `False` | Require equal formal charges |
| `match_isotope` | `False` | Compare specified isotope labels; unspecified labels act as wildcards |
| `match_bond_order` | `"strict"` | `"loose"` and `"any"` both relax bond order; other chemistry constraints still apply |
| `ring_matches_ring_only` | `False` | Match ring atoms/bonds only to ring atoms/bonds |
| `complete_rings_only` | `False` | Preserve complete query rings in a result |
| `use_chirality` | `False` | Compare resolved R/S labels and mapped ligand parity; unspecified target tags remain permissive |
| `use_bond_stereo` | `False` | Compare annotated double-bond stereo |
| `tautomer_aware` | `False` | Relax eligible tautomer bonds while preserving elements unless explicitly disabled |
| `connected_only` | `True` | Return one connected query fragment |
| `induced` | `False` | Preserve query nonedges as well as edges |
| `maximize_bonds` | `False` | Rank mapped bonds before atom count |
| `timeout_ms` | `10000` | Per-search budget; native `MCSOptions.timeout_ms=-1` requests adaptive sizing |
| `atom_weights` | `[]` | Finite per-query scores; negative scores may prefer a smaller or empty result |
| `min_fragment_size` | `1` | Minimum retained query fragment size |
| `max_fragments` | native limit | Maximum retained fragments |
| `max_stage` | `5` | Effort setting; reducing it can reduce result size |

Native weights scale sums to integer millipoints by truncation toward zero and must fit
the supported integer score range. Weight vectors follow the query input's atom ordering, including original
RDKit indices; high-level wrappers translate them into native graph order.
Raw native calls use graph order. Weighted constrained batches with RDKit
inputs currently support one query; multiple queries must be converted to
graphs with weights supplied in graph order.

Raw functions in `smsd._smsd` accept `ChemOptions` and `MCSOptions` objects:

```python
q = smsd.parse_smiles("CC")
t = smsd.parse_smiles("CCC")
chem = smsd.ChemOptions()
opts = smsd.MCSOptions()
opts.timeout_ms = 1000
mapping = smsd._smsd.find_mcs(q, t, chem, opts)
```

Raw bindings return graph indices and do not translate RDKit atom ordering.
`validate_mapping` checks a witness with the supplied options;
`is_mapping_maximal` checks whether a compatible pair can be added. Maximality
is weaker than a proof of the maximum objective.

## Batch search and thread safety

```python
query = smsd.parse_smiles("C1CC1")
targets = [smsd.parse_smiles(s) for s in ["CC1CC1", "C1CC1", "CCC"]]
mappings = smsd.batch_mcs(query, targets, timeout_ms=1000, num_threads=2)
sizes = smsd.batch_mcs_size(query, targets, timeout_ms=1000, num_threads=2)
assert sizes == [len(m) for m in mappings]
hits = smsd.batch_substructure(query, targets, timeout_ms=1000, num_threads=2)
assert hits == [True, True, False]
```

The four core batch functions accept SMILES, graph or RDKit inputs.
`batch_mcs` and `batch_mcs_size` preserve chemistry/search keywords.
`num_threads=0` uses OpenMP defaults; `1` runs sequentially. `timeout_ms`
applies to each pair, independently of the worker count.

Bindings retain Python graph owners and pass native references through these
core batches and `SmartsQuery.matches_many`, avoiding graph copies. Native
search releases the GIL; progress callbacks reacquire it. Batch setup warms
canonical, ring, fingerprint, neighborhood and pharmacophore caches before
workers start. Owning collections such as `TargetCorpus` still copy graphs.

For concurrent calls that share graphs, call `smsd.prewarm_graph(graph)` once
before starting threads. Do not mutate graph properties or shared option
objects while a native call is running. Prewarming establishes initialized
caches; it does not promise a fixed percentage improvement.

Constrained reaction batches also expose which target was selected:

```python
results = smsd.batch_mcs_constrained(
    ["NO"], ["C", "CNO"], return_target_indices=True, timeout_ms=1000,
)
target_index, mapping = results[0]
assert target_index == 1 and len(mapping) == 2
```

Results retain query order. Queries are processed by decreasing size and
claim target atoms without overlap. This greedy assignment across queries
is not a proof of the globally best reaction mapping. Without
`return_target_indices=True`, the existing API returns mappings only; empty
results use target index `-1`.

## RDKit conversion and caches

### A molecule-first workflow

Keep RDKit molecules as the inputs when their atom order is needed for
annotations, reaction maps or drawing. The high-level functions translate
both sides back to those original indices:

```python
from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("CO")
targets = [Chem.MolFromSmiles(s) for s in ("CCO", "CCN")]
hits = [smsd.substructure_rdkit(query, target) for target in targets]
assert len(hits[0]) == 2 and hits[1] == {}
mapping = smsd.mcs_rdkit_native(query, targets[0], timeout_ms=1000)
query_atoms = sorted(mapping)
target_atoms = [mapping[index] for index in query_atoms]
assert all(query.GetAtomWithIdx(a).GetAtomicNum() ==
           targets[0].GetAtomWithIdx(b).GetAtomicNum()
           for a, b in zip(query_atoms, target_atoms))
```

| RDKit convention | SMSD counterpart | Result contract |
|---|---|---|
| `Chem.MolFromSmiles` | `parse_smiles`, or pass a RDKit molecule directly | Native parser order for strings; original order for RDKit inputs |
| `GetSubstructMatch` | `substructure_rdkit` / `find_substructure` | Dictionary from query index to target index |
| `GetSubstructMatches` | `find_substructure(..., max_results=N)` | List of dictionaries when `N > 1` |
| `rdFMCS.FindMCS` | `mcs_rdkit_native` / `find_mcs` | A selected atom mapping; no cancellation or optimality certificate |
| `MolFromSmarts` | `compile_smarts` | Native SMARTS query; dialect differences require checking |

SMSD and RDKit fingerprints are separate implementations. Equal names or
radius settings do not establish identical bit vectors or interchangeable
similarity thresholds. Likewise, compare MCS engines only after aligning
chemistry, connectivity, objective and witness validity.

### Conversion lifetime and atom ordering

```python
from rdkit import Chem
import smsd

query = Chem.MolFromSmiles("NCCO")
target = Chem.MolFromSmiles("CC(O)CN")
query = Chem.RenumberAtoms(query, [3, 1, 0, 2])
mapping = smsd.mcs_rdkit_native(query, target, timeout_ms=1000)
for a, b in mapping.items():
    assert query.GetAtomWithIdx(a).GetAtomicNum() == target.GetAtomWithIdx(b).GetAtomicNum()

graph = smsd.from_rdkit(query)
index_map = smsd.get_index_map(graph)  # graph index -> original RDKit index
smsd.clear_cache()
assert smsd.get_index_map(graph) == index_map
```

`mcs_rdkit`, `mcs_rdkit_native`, `batch_mcs_rdkit` and `substructure_rdkit`
use the original RDKit indices. Translation preserves the selected native
mapping; it does not rematch a derived SMARTS or discard explicitly relaxed
atom matches.

RDKit conversions are cached per live molecule using weak keys. A binary
structure signature invalidates a cached conversion after molecule edits.
`use_cache=False` bypasses reuse. Clearing conversion caches preserves index
metadata for graphs the caller still holds. Returned index maps are copies.
A failed index-order reconstruction raises an error instead of inventing an
identity mapping. Empty or hydrogen-only RDKit molecules are rejected.

### Reuse compiled SMARTS

```python
pattern = smsd.compile_smarts("[#6]-[#8]")
graphs = [smsd.parse_smiles(s) for s in ("CCO", "CCN")]
assert pattern.matches_many(graphs) == [True, False]
```

Compile once for repeated queries and reuse parsed target graphs. This keeps
compilation and conversion out of the matching loop. The binding benchmark
reports those setup costs separately from search latency.

## Progress callbacks

```python
seen = []
mapping = smsd.find_mcs_progressive(
    "CCO", "CCN", timeout_ms=1000,
    on_progress=lambda result, size, elapsed_ms: seen.append((size, elapsed_ms)),
)
assert seen[-1][0] == len(mapping)
```

The current native implementation performs one search and invokes the callback
with the final mapping. The API name is retained for compatibility; it does
not currently emit intermediate stage results. Caller options and atom-index
translation are preserved. Callback exceptions propagate to Python.

## Other APIs

- SMARTS: `compile_smarts`, `smarts_match`, `smarts_find_all`, `find_mcs_smarts`.
- Fingerprints: `fingerprint`, `circular_fingerprint`, `topological_torsion`
  and count variants; `tanimoto_coefficient` and count similarity functions.
- I/O: `parse_smiles`, `to_smiles`, `read_mol_block`, `write_mol_block`,
  `read_sdf`, `write_sdf`, `write_mol_block_v3000`.
- Chemistry: `assign_rs`, `assign_ez`, `assign_cip`, `murcko_scaffold`.
- Graph operations: `extract_subgraph`, `split_components`, `count_components`.
- Mapping enumeration: `find_all_mcs`, `canonicalize_mapping` and validation.
- Depiction: `depict_svg`, `depict_pair`, `generate_coords_2d` and layout helpers.

Exact symmetry canonicalization can fail when molecular automorphism
generators are incomplete or the orbit exceeds resource limits. Native
`length_error` is translated to Python `ValueError`; expiry is `RuntimeError`.
Internal enumeration conservatively retains raw mappings when symmetry proof
is unavailable, so chemically equivalent embeddings may remain in the output.

See [examples](EXAMPLES.md), the [C++ guide](CPP.md) and the
[benchmark protocol](../benchmarks/README.md) for detailed usage and scope.

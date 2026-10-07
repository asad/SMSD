# SMSD 7.2.2 Python guide

Substructure search, maximum common substructure (MCS), fingerprints, molecular
I/O and SVG drawing. Core functions do not require RDKit or Java. RDKit is
optional for molecule conversion and drawing.

## Install

Release wheels use CPython 3.14 with CPU/OpenMP support:

| Platform | Architecture | Requirement |
|---|---|---|
| Windows | x86_64 / AMD64 | Windows 10 or later |
| Linux | x86_64 | glibc 2.28+ |
| macOS | arm64 / Apple Silicon | macOS 26+ |

Install from PyPI when 7.2.2 is listed, or download the matching wheel from the
[GitHub release](https://github.com/asad/SMSD/releases/tag/v7.2.2):

```bash
python -m pip install smsd==7.2.2
# Optional RDKit interoperability and drawing:
python -m pip install rdkit
```

## Build from source

Source metadata allows Python 3.9 or later. Other interpreters and architectures
require a source build. Run at the repository root with a C++17 compiler and
CMake 3.18 or later. The root `pyproject.toml` builds `python/smsd/` together
with the native extension in `cpp/`:

```bash
python -m pip install build
python -m build --wheel \
  -Ccmake.define.SMSD_BUILD_METAL=OFF \
  -Ccmake.define.SMSD_BUILD_CUDA=OFF
```

Metal and CUDA are optional source-build features requiring compatible tools
and hardware. Release wheels use CPU/OpenMP. Check the active backend with
`smsd.gpu_device_info()`; core MCS and batch matching run on the CPU.

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
extra target edges between mapped atoms. Equal atom counts from different
engines do not establish that their atom and bond mappings are equivalent.

The default `strategy="auto"` selects a solver that supports the requested
options. `strategy="native"` selects the full native solver.
`strategy="lightweight"` supports one connected atom objective, a timeout,
ring matching and strict/any bond order; unsupported options raise `ValueError`.

A time budget produces the best valid mapping found by the pipeline. The
public mapping API has no cancellation flag or optimality certificate.
Reaching the budget does not prove maximum size. A false/empty substructure result can mean no witness was found within
the budget, rather than a completed proof of absence. Timeouts are cooperative;
preparation and a search step can take the elapsed time beyond the budget.

## Chemistry and search options

```python
import smsd

mapping = smsd.find_mcs(
    "CCC", "C1CC1", strategy="native", timeout_ms=1000,
    induced=True, connected_only=True, match_bond_order="strict",
)
assert len(mapping) == 2

weighted = smsd.find_mcs(
    "CCC", "CCC", strategy="native", timeout_ms=1000, atom_weights=[10.0, -30.0, 1.0],
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
| `tautomer_aware` | `False` | Relax eligible tautomer bonds while preserving elements; bond order remains an explicit setting |
| `connected_only` | `True` | Return one connected query fragment |
| `induced` | `False` | Preserve query nonedges as well as edges |
| `maximize_bonds` | `False` | Rank mapped bonds before atom count |
| `timeout_ms` | `10000` | Per-search budget; native `MCSOptions.timeout_ms=-1` requests adaptive sizing |
| `atom_weights` | `[]` | Finite per-query scores; negative scores may prefer a smaller or empty result |
| `min_fragment_size` | `1` | Minimum retained query fragment size |
| `max_fragments` | native limit | Maximum retained fragments |
| `max_stage` | `5` | Effort setting; reducing it can reduce result size |

Weights follow the query input's atom order, including original RDKit indices.
Scores are converted to integer millipoints and must fit the supported range.
Weighted constrained batches with RDKit inputs support one query; for multiple
queries, convert them to graphs and provide weights in graph order.

Raw functions in `smsd._smsd` accept `ChemOptions` and `MCSOptions` objects:

```python
import smsd

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
import smsd

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

Before concurrent calls share a graph, call `smsd.prewarm_graph(graph)` once.
Do not mutate shared graphs or option objects while matching is running. Ordinary
batch calls prepare their input graphs before starting workers.

Constrained reaction batches also expose which target was selected:

```python
import smsd

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

## Working with RDKit

### Preserve input atom indices

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

Conversions are reused while the RDKit molecule remains unchanged.
`from_rdkit(mol, use_cache=False)` bypasses that reuse. Clearing caches preserves
index metadata for graphs still held by the caller. Returned index maps are
copies; empty or hydrogen-only RDKit molecules are rejected.

### Reuse compiled SMARTS

```python
import smsd

pattern = smsd.compile_smarts("[#6]-[#8]")
graphs = [smsd.parse_smiles(s) for s in ("CCO", "CCN")]
assert pattern.matches_many(graphs) == [True, False]
```

Compile once for repeated queries and reuse parsed target graphs.

## Progress callbacks

```python
import smsd

seen = []
mapping = smsd.find_mcs_progressive(
    "CCO", "CCN", timeout_ms=1000,
    on_progress=lambda result, size, elapsed_ms: seen.append((size, elapsed_ms)),
)
assert seen[-1][0] == len(mapping)
```

The callback receives the final mapping once; it does not report intermediate
stages. Callback exceptions propagate to the caller.

## Other APIs

- SMARTS: `compile_smarts`, `smarts_match`, `smarts_find_all`, `find_mcs_smarts`.
- Fingerprints: `fingerprint`, `circular_fingerprint`, `topological_torsion`
  and count variants; `tanimoto_coefficient` and count similarity functions.
- Fingerprint storage: `to_hex`, `from_hex`, `to_binary_string`, `counts_to_array`.
- I/O: `parse_smiles`, `to_smiles`, `read_mol_block`, `write_mol_block`,
  `read_sdf`, `write_sdf`, `write_mol_block_v3000`.
- Chemistry: `assign_rs`, `assign_ez`, `assign_cip`, `murcko_scaffold`.
- Graph operations: `extract_subgraph`, `split_components`, `count_components`.
- Mapping enumeration: `find_mcs(..., max_results=10)`, `canonicalize_mapping`
  and `validate_mapping`.
- Depiction: `depict_svg`, `depict_pair`, `generate_coords_2d` and layout helpers.

Exact symmetry canonicalisation can exceed the time or resource budget.
`ValueError` or `RuntimeError` reports that the canonical representative was
not established. Mapping enumeration can retain chemically equivalent results
when symmetry reduction is incomplete.

See [examples](EXAMPLES.md), the [C++ guide](CPP.md) and the
[benchmark protocol](../benchmarks/README.md) for detailed usage and scope.

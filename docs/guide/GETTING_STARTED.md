# Getting started

This chapter installs SMSD Pro 8.1.0, runs a first substructure search and a
first maximum common substructure (MCS) search, and sets out the conventions
that every other chapter relies on: what a result proves, how molecules are
given and numbered, the default chemistry, and time limits and threads.

## Install

Install the Python package from PyPI:

```bash
python -m pip install "smsd==8.1.0"
python -c "import smsd; print(smsd.__version__); assert smsd.parse_smiles('CC').n == 2"
```

The package has no required dependencies. Searches, fingerprints, prepared
target corpora, streamed screens and corpus searches all work without the
optional packages. Extras add the packages that particular features need, for
example `python -m pip install "smsd[rdkit,numpy]==8.1.0"`:

| Extra | Installs | Needed for |
| --- | --- | --- |
| `rdkit` | `rdkit` | RDKit molecules as input, and RDKit depiction |
| `numpy` | `numpy` | Packed fingerprint indexes, `butina` and `maxmin` clustering, NumPy feature arrays |
| `pandas` | `pandas` | Molecule columns in pandas data frames |
| `polars` | `polars` | Molecule columns in Polars data frames |
| `scipy` | `numpy`, `scipy` | Sparse fingerprint matrices |
| `sklearn` | `numpy`, `scikit-learn` | Fingerprint features in scikit-learn pipelines |
| `torch` | `numpy`, `torch` | Fingerprint tensors |
| `pyg` | `torch`, `torch-geometric` | Molecular graphs for graph learning |
| `arrow` | `numpy`, `pyarrow` | Arrow fingerprint columns |
| `mcp` | `mcp` | The local molecular tools server, `smsd-mcp` |

For NVIDIA GPU screening on Linux x86_64, add the optional CUDA package of the
same version, `python -m pip install "smsd==8.1.0" "smsd-cuda==8.1.0"`; the
[CPU and GPU screening](CPU_GPU.md) chapter covers it. Interpreters and
platforms are listed in the
[supported-platforms table](../HOWTO-INSTALL.md#supported-platforms), which also
leads to the source builds.

## A first search

A prepared query is built once and then matched against any number of targets.
This program finds phenol in 4-tert-butylphenol and the MCS of
N-methylacetamide and propanamide:

```python
import smsd

# Substructure: phenol in 4-tert-butylphenol
result = smsd.Substructure.find("Oc1ccccc1").match("CC(C)(C)c1ccc(O)cc1")
assert result.found and result.optimal and result.complete
assert result.status == smsd.SearchStatus.COMPLETED
mapping = dict(sorted(result.best.atom_mapping.items()))
assert mapping == {0: 8, 1: 7, 2: 6, 3: 5, 4: 4, 5: 10, 6: 9}
print("substructure:", result.found, result.status.name, result.optimal, result.complete)
print("atom mapping:", mapping)

# MCS: N-methylacetamide and propanamide
common = smsd.MCS.find("CC(=O)NC").with_timeout(1000).match("CCC(=O)N")
assert common.found and common.optimal and common.complete
assert common.status == smsd.SearchStatus.COMPLETED
common_mapping = dict(sorted(common.best.atom_mapping.items()))
assert common.best.atom_count == 4 and common_mapping == {0: 1, 1: 2, 2: 3, 3: 4}
assert common.best.to_smiles("CC(=O)NC") == "C(C)(N)=O"
print("MCS:", common.found, common.status.name, common.optimal, common.complete)
print("MCS atoms:", common.best.atom_count, common_mapping)
print("common graph:", common.best.to_smiles("CC(=O)NC"))
```

The program prints:

```text
substructure: True COMPLETED True True
atom mapping: {0: 8, 1: 7, 2: 6, 3: 5, 4: 4, 5: 10, 6: 9}
MCS: True COMPLETED True True
MCS atoms: 4 {0: 1, 1: 2, 2: 3, 3: 4}
common graph: C(C)(N)=O
```

The atom mapping runs from query atom to target atom: the phenol oxygen, atom
0 of the query, is atom 8 of the target. `to_smiles(query)` writes the matched
part of the query as SMILES, with the cut valences capped.

## Reading a result

Each `match` and `match_all` call returns a `SearchResult`, and `match_many`
and `MCS.match_constrained` return one per query or target. The other entry
points wrap one or report their own fields:

- `match_corpus` returns one `CorpusSearchResult` per target, holding a
  `SearchResult`, or `None` when the target is proved not to qualify (see
  [Search a prepared corpus](SUBSTRUCTURE.md#search-a-prepared-corpus)).
- `MCS.match_with_smarts` returns a `SmartsResult` with the `SearchResult`
  and a SMARTS string (see
  [Export the common graph](MCS.md#export-the-common-graph)).
- `Substructure.decompose_r_groups` returns one dictionary per target, with
  the `SearchResult` under `search`.
- `MCS.consensus` returns a `ConsensusResult` with its own `status`,
  `optimal`, `complete` and `found` fields, and `optimal` is always False
  (see [Common graph of a series](MCS.md#common-graph-of-a-series)).

The fields of a `SearchResult`:

| Field | Type | Meaning |
| --- | --- | --- |
| `matches` | tuple of `Match` | The matches found, at most one from `match` |
| `status` | `SearchStatus` | `COMPLETED`, `TIMED_OUT`, `CANCELLED` or `LIMIT_REACHED` |
| `optimal` | bool | The answer is proved: for a substructure search, a match exists or a completed search proved absence; for an MCS, the maximum is proved |
| `complete` | bool | The search, or the enumeration of `match_all`, ran to its end |
| `score` | float | Substructure: the number of query atoms when found, otherwise 0.0; MCIS: atoms; MCES: bonds; weighted MCIS: the weight sum |
| `nodes_visited` | int | Search effort, for diagnostics |
| `elapsed_ms` | int | Elapsed milliseconds, which vary from run to run |
| `found` | bool | Whether `matches` is not empty |
| `best` | `Match` or `None` | The first match |

A `Match` has a read-only `atom_mapping` from query atom to target atom,
`bond_mapping` as `(query_first, query_second, target_first, target_second)`
tuples, `fragments` (the components of the matched bonds, as query atoms),
`atom_count`, `bond_count`, and the methods `to_smiles(query)`,
`to_graph(query)` and `to_smarts(query)`.

The rule that the rest of the guide relies on:

- Any match returned is a valid embedding or common subgraph, whatever the
  status.
- An empty result proves absence only when `status` is `COMPLETED` and
  `optimal` is `True`.
- An MCS is proved maximum only when `optimal` is `True`.
- A search stopped by its time limit, node limit or cancellation token proves
  nothing about what it did not find.

A time limit of 0 stops a search at once, which shows the stopped state:

```python
import smsd

anthracene = "c1ccc2cc3ccccc3cc2c1"                     # anthracene
benzene = smsd.Substructure.find("c1ccccc1")            # benzene

stopped = benzene.with_timeout(0).match(anthracene)
assert stopped.status == smsd.SearchStatus.TIMED_OUT
assert not (stopped.found or stopped.optimal or stopped.complete)
print(stopped.status.name, stopped.found, stopped.optimal, stopped.complete)

finished = benzene.match(anthracene)
assert finished.status == smsd.SearchStatus.COMPLETED
assert finished.found and finished.optimal and finished.complete and finished.score == 6.0
print(finished.status.name, finished.found, finished.optimal, finished.complete, finished.score)
```

The program prints:

```text
TIMED_OUT False False False
COMPLETED True True True 6.0
```

## Molecule input

Searches and most other functions accept a molecule in any of these forms; a
few low-level functions, such as `circular_fingerprint`, take a `MolGraph`
only.

- A SMILES string. It is parsed on every call, so when a molecule is used
  repeatedly, parse it once with `smsd.parse_smiles(smiles, *, lenient=False)`.
  `lenient=True` recovers from malformed text, such as an unclosed ring, that
  otherwise raises `ValueError`.
- A `MolGraph`, from `parse_smiles`, `read_mol_file`, `read_mol_block` or
  `read_sdf`.
- An RDKit molecule, when RDKit is installed.

`read_sdf(filename)` reads every record of a file at once and
`write_sdf(molecules, filename)` writes V2000 records. Each `MolGraph` carries
an editable `name` (the record title) and `properties` (the SD fields), and
`write_sdf` writes both. A malformed record is read as an empty graph with no
warning, so keep the graphs with `n > 0`; a missing file raises `ValueError`.

```python
import smsd
from pathlib import Path
from tempfile import TemporaryDirectory

aspirin = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")      # aspirin
aspirin.name = "aspirin"
aspirin.properties = {"SOURCE": "inline"}
paracetamol = smsd.parse_smiles("CC(=O)Nc1ccc(O)cc1")     # paracetamol
paracetamol.name = "paracetamol"
paracetamol.properties = {"SOURCE": "inline"}

with TemporaryDirectory() as directory:
    path = Path(directory) / "drugs.sdf"
    smsd.write_sdf([aspirin, paracetamol], path)
    with path.open("a", encoding="utf-8") as handle:    # a damaged third record
        handle.write("damaged\nnot a MOL block\n\nM  END\n$$$$\n")
    records = smsd.read_sdf(path)

sizes = [graph.n for graph in records]
assert sizes == [13, 11, 0]
molecules = [graph for graph in records if graph.n > 0]
names = [graph.name for graph in molecules]
fields = dict(molecules[0].properties)
assert names == ["aspirin", "paracetamol"] and fields == {"SOURCE": "inline"}
print("atoms per record:", sizes)
print("kept:", names)
print("fields of the first:", fields)
```

The program prints:

```text
atoms per record: [13, 11, 0]
kept: ['aspirin', 'paracetamol']
fields of the first: {'SOURCE': 'inline'}
```

A line of a SMILES file that carries a title, or even trailing whitespace, is
not valid SMILES, so split off the first column before parsing.
`parse_smiles("")` returns an empty graph, while searches, fingerprints and
target corpora raise `ValueError` for an empty string:

```python
import smsd

lines = ["CC(=O)Oc1ccccc1C(=O)O aspirin", "CC(=O)Nc1ccc(O)cc1 paracetamol"]
try:
    smsd.parse_smiles(lines[0])
except ValueError as error:
    message = str(error)
assert message == "Unexpected character at position 21: ' '"
print("whole line:", message)

graphs = []
for line in lines:
    smiles, title = line.split(maxsplit=1)
    graph = smsd.parse_smiles(smiles)
    graph.name = title
    graphs.append(graph)
assert [(graph.name, graph.n) for graph in graphs] == [("aspirin", 13), ("paracetamol", 11)]
print("parsed:", [(graph.name, graph.n) for graph in graphs])

recovered = smsd.to_smiles(smsd.parse_smiles("C1CC", lenient=True))
assert recovered == "C(C)C"
print("lenient C1CC:", recovered)
```

The program prints:

```text
whole line: Unexpected character at position 21: ' '
parsed: [('aspirin', 13), ('paracetamol', 11)]
lenient C1CC: C(C)C
```

Explicit hydrogens, such as `[H]` in SMILES or hydrogen atoms in a MOL block,
are atoms of the graph, whereas bracket hydrogen counts are not. A molecule
query with an explicit hydrogen therefore needs an explicit hydrogen in the
target, and an explicit-hydrogen record gives different fingerprints from its
hydrogen-suppressed form. Normalise hydrogens the same way for every molecule
before comparing them.

## Atom numbering

Search results and atom weights use each input's own numbering, separately for
the query and the target: parser order for SMILES, graph order for a
`MolGraph`, and RDKit order for an RDKit molecule. Sort the items of
`atom_mapping` before printing them, because their order is not sorted for
RDKit input:

```python
from rdkit import Chem
import smsd

target = Chem.MolFromSmiles("OCC1=CC=CC=C1")              # benzyl alcohol, RDKit order
result = smsd.Substructure.find("c1ccccc1CO").match(target)
mapping = dict(sorted(result.best.atom_mapping.items()))
assert mapping == {0: 3, 1: 4, 2: 5, 3: 6, 4: 7, 5: 2, 6: 1, 7: 0}
oxygen = target.GetAtomWithIdx(mapping[7]).GetSymbol()
assert oxygen == "O"
print("RDKit-numbered mapping:", mapping)
print("query atom 7 is target atom", mapping[7], oxygen)
```

The program prints:

```text
RDKit-numbered mapping: {0: 3, 1: 4, 2: 5, 3: 6, 4: 7, 5: 2, 6: 1, 7: 0}
query atom 7 is target atom 0 O
```

The [Python guide](../PYTHON.md#input-atom-numbering) describes the RDKit
index helpers.

## Default chemistry

Every search starts from the same chemistry options:

| Option | Default | Effect |
| --- | --- | --- |
| `match_atom_type` | `True` | Elements must agree |
| `match_bond_order` | `BondOrderMode.STRICT` | Bond orders must agree |
| `aromaticity_mode` | `AromaticityMode.FLEXIBLE` | An aromatic bond also matches a single or double bond |
| `match_formal_charge`, `match_isotope` | `False` | Charges and isotopes are not compared |
| `use_chirality`, `use_bond_stereo` | `False` | Stereochemistry is not compared |
| `ring_matches_ring_only` | `False` | Ring atoms may match chain atoms |
| `complete_rings_only` | `False` | An MCS may include part of a ring |
| `ring_fusion_mode` | `RingFusionMode.IGNORE` | Ring fusion is not compared |
| `tautomer_aware` | `False` | Tautomeric forms are not treated as equivalent |

Two consequences surprise new users:

- Because aromatic bonds also match single and double bonds, benzene is found
  in cyclohexanol.
- A molecule query does not compare hydrogen counts, so the phenol graph also
  embeds in the ester oxygen of aspirin. A SMARTS pattern such as `c[OX2H]`
  states the hydrogen and does not.

```python
import smsd

benzene = smsd.Substructure.find("c1ccccc1")              # benzene
cyclohexanol = "OC1CCCCC1"                                 # cyclohexanol
aspirin = "CC(=O)Oc1ccccc1C(=O)O"                          # aspirin
flexible = benzene.match(cyclohexanol).found
assert flexible
print("benzene in cyclohexanol:", flexible)

options = smsd.Substructure.Options()
options.chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
strict = benzene.with_options(options).match(cyclohexanol).found
assert not strict
print("with STRICT aromaticity:", strict)

molecule_query = smsd.Substructure.find("Oc1ccccc1").match(aspirin).found      # phenol
smarts_query = smsd.Substructure.smarts("c[OX2H]").match(aspirin).found
assert molecule_query and not smarts_query
print("phenol in aspirin:", molecule_query, "| c[OX2H] in aspirin:", smarts_query)
```

The program prints:

```text
benzene in cyclohexanol: True
with STRICT aromaticity: False
phenol in aspirin: True | c[OX2H] in aspirin: False
```

The [substructure chapter](SUBSTRUCTURE.md#chemistry-options) describes every
option and the named profiles.

## Time limits and threads

Each search has its own limits, set on an options object:

- `timeout_ms`, default `-1`, selects 10,000 ms per target for a substructure
  search, and 500 + 2 × query atoms × target atoms ms, at most 30,000 ms, for
  an MCS. A positive value applies to each query and target pair; `0` stops at
  once.
- `node_limit`, default `0` (no cap): reaching it reports `LIMIT_REACHED`.
- `cancellation`, default `None`: a `CancellationToken` whose `cancel()`, from
  any thread, stops the search with `CANCELLED`. A token stays cancelled, so
  use a new one for each request.

Options objects take no keyword arguments: create one and set its attributes.
A misspelt attribute raises `AttributeError`. `with_options` and
`with_timeout` return a new query with a copy of the options, apart from the
cancellation token, which stays shared.

Batch calls such as `match_many` take `num_threads`; `0`, the default, uses the
OpenMP default, which follows the `OMP_NUM_THREADS` environment variable.
Results do not depend on the thread count. Pass a non-negative integer:
`TargetCorpus.screen` and `TargetCorpus.prewarm` do not check it.

```python
import smsd

anthracene = "c1ccc2cc3ccccc3cc2c1"                       # anthracene
benzene = smsd.Substructure.find("c1ccccc1")              # benzene

options = smsd.Substructure.Options()
assert (options.timeout_ms, options.node_limit, options.cancellation) == (-1, 0, None)
try:
    smsd.Substructure.Options(timeout_ms=2000)
except TypeError:
    print("keyword arguments: TypeError")
else:
    raise AssertionError("keyword arguments were accepted")
options.timeout_ms = 2000                                 # per target
options.node_limit = 3
capped = benzene.with_options(options).match(anthracene)
assert capped.status == smsd.SearchStatus.LIMIT_REACHED
assert not (capped.found or capped.complete)
print("node limit:", capped.status.name, capped.found, capped.complete)

token = smsd.CancellationToken()
options = smsd.Substructure.Options()
options.cancellation = token
cancellable = benzene.with_options(options)
token.cancel()                                            # the query shares the token
cancelled = cancellable.match(anthracene)
assert cancelled.status == smsd.SearchStatus.CANCELLED and not cancelled.found
print("cancelled:", cancelled.status.name, cancelled.found)

targets = ["Oc1ccccc1", "OC1CCCCC1", "CCO", "Cc1ccccc1"]   # phenol, cyclohexanol, ethanol, toluene
one_thread = [result.found for result in benzene.match_many(targets, num_threads=1)]
default = [result.found for result in benzene.match_many(targets)]
assert one_thread == default == [True, True, False, True]
print("match_many:", default)
```

The program prints:

```text
keyword arguments: TypeError
node limit: LIMIT_REACHED False False
cancelled: CANCELLED False
match_many: [True, True, False, True]
```

## Next steps

- [Substructure search](SUBSTRUCTURE.md): embeddings, many targets, prepared
  corpora and chemistry options.
- [SMARTS queries](SMARTS.md): patterns, primitives and limits.
- [Maximum common substructure](MCS.md): MCIS and MCES, fragments, proofs and
  export.
- [Fingerprints and similarity](FINGERPRINTS.md): fingerprint families and
  similarity coefficients.
- [Screening a library](SCREENING.md): similarity and substructure screens with
  exact confirmation.
- [Packed fingerprint index](PACKED_INDEX.md): building, reopening and
  searching an index file.
- [CPU and GPU screening](CPU_GPU.md): where work runs and how to choose.
- [Command line and containers](COMMAND_LINE.md): the command-line tool and the
  Docker images.
- [C++ and Java](OTHER_LANGUAGES.md): the same tasks in the other languages.

The [Python guide](../PYTHON.md) lists every function and option, and the
[Python examples](../EXAMPLES.md) cover depiction, coordinates, stereo
assignment and other tasks beyond searching.

# Screening a library

Screening narrows a library down to the molecules worth a closer look. SMSD Pro
offers these routes: a ranked list for a handful of molecules, a fingerprint
threshold screen held in memory or streamed from a file, a saved packed index,
an exact substructure screen and an MCS-similarity screen. This chapter shows
which route suits which task, how to prepare the library, how to run each
screen and how to save the hits. Every example uses the same inline library of
533 molecules: 14 benzodiazepines, 7 common drugs and 512 generated
trisubstituted benzenes, so each program runs as it stands, with no data files.

## Choose a route

| Task | Use | Fingerprint | Prepare first | Needs NumPy | GPU |
| --- | --- | --- | --- | --- | --- |
| Rank a small list | `batch_similarity` | ECFP; radius 2 and 2048 bits by default | Nothing | No | No |
| Threshold screen in memory | `TargetCorpus.screen` or `screen_many` | ECFP radius 2, 2048 bits, fixed | `prewarm()` after every addition | No | With `smsd-cuda`, automatic on a Tesla T4 for 100,000 targets or more, or on any supported device from `SMSD_CUDA_MIN_TARGETS` targets; see the [device rule](CPU_GPU.md#automatic-use-in-prepared-screens) |
| Screen a large file once | `stream_screen_smiles` or `stream_screen_many_smiles` | ECFP radius 2, 2048 bits, fixed | Nothing | No | As above, chunk by chunk |
| Search a saved index many times | `FingerprintCorpus`; see the [packed index](PACKED_INDEX.md) | Any binary specification | Build the file once | Yes | On request, with `backend="cuda"` or `"auto"` |
| Exact substructure | `Substructure.match_corpus` | None: exact search | A `TargetCorpus` | No | Optional pre-screen, with `gpu_prefilter=True` |
| MCS similarity | `screen_targets`, then `MCS.match_many`; or `MCS.match_corpus` | None: an atom-count bound, then the exact MCS | A list or a `TargetCorpus` | No | `screen_targets` runs on the CPU; `MCS.match_corpus` takes `gpu_prefilter=True` |

The in-memory, streamed and packed-index screens give the same hits for the
same query and threshold when the index is built with ECFP radius 2 and 2048
bits, because they then compare the same fingerprints. Graph matching always
runs on the CPU, and every GPU route returns the same results as the CPU.
[CPU and GPU screening](CPU_GPU.md) explains when a device is used.

## Prepare the library

A screen is only as consistent as its input. Four points matter before any
corpus is built:

- Malformed SD records. `smsd.read_sdf` reads the whole file and returns a
  malformed record as an empty graph, with no warning. An empty graph still
  takes part in a screen and passes one with a threshold of 0, so keep only
  graphs with `graph.n > 0`.
- SMILES files. A line that holds a name after the SMILES does not parse as
  SMILES, so split off the first column. There is no SMILES-file reader; a
  generator over the lines serves, and keeps memory flat for the streamed
  screens.
- Identifiers. A `TargetCorpus` or `FingerprintCorpus` stores no names and
  cannot be indexed or iterated. Keep the names in a list in the same order;
  every hit index refers to that list.
- Hydrogens. Explicit hydrogens are graph atoms, so a record written with them
  has quite different fingerprints from the same molecule without them. Write
  or export every record in the same form, with implicit hydrogens.

The program below defines the library, writes it as an SD file with a title and
a `SOURCE` field per record and as a SMILES file with `SMILES name` lines, adds
a malformed record, then reads both files back:

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

with TemporaryDirectory() as directory:
    sdf_path = Path(directory) / "library.sdf"
    smi_path = Path(directory) / "library.smi"

    # Write an SD file with a title and a SOURCE field per record ...
    graphs = []
    for name, text in library:
        graph = smsd.parse_smiles(text)
        graph.name = name
        graph.properties = {"SOURCE": "named" if name in NAMED else "generated"}
        graphs.append(graph)
    smsd.write_sdf(graphs, sdf_path)
    with sdf_path.open("a", encoding="utf-8") as stream:
        stream.write("garbage\nnot a molfile\n\nM  END\n$$$$\n")  # malformed
    # ... and a SMILES file with "SMILES name" lines.
    smi_path.write_text("".join(f"{text} {name}\n" for name, text in library),
                        encoding="utf-8")

    # A malformed SD record comes back as an empty graph: drop it.
    records = smsd.read_sdf(sdf_path)
    kept = [graph for graph in records if graph.n > 0]
    print("SD records:", len(records), "| with atoms:", len(kept))
    print("first record:", kept[0].name, dict(kept[0].properties))

    # A SMILES line holds the SMILES first; the rest of the line is the name.
    lines = smi_path.read_text(encoding="utf-8").splitlines()
    rows = [line.split(maxsplit=1) for line in lines if line.strip()]
    file_smiles = [row[0] for row in rows]
    file_names = [row[1] if len(row) > 1 else "" for row in rows]
    print("SMILES lines:", len(rows), "| a name with a space:", file_names[15])

assert len(records) == 534 and len(kept) == 533
assert [graph.name for graph in kept] == names == file_names
assert dict(kept[0].properties) == {"SOURCE": "named"}
assert file_smiles == smiles
```

The program prints:

```text
SD records: 534 | with atoms: 533
first record: diazepam {'SOURCE': 'named'}
SMILES lines: 533 | a name with a space: salicylic acid
```

Duplicates are common in assembled libraries. `smsd.canonical_hash` groups
candidate duplicates quickly, but a hash is a screening key rather than proof
of identity, so confirm each candidate with `smsd.same_canonical_graph`. The
generated benzenes carry their three substituents in equivalent positions, so
the 512 SMILES describe only 120 distinct structures, and the whole library
holds 141:

```python
import itertools

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

# Group records by canonical hash, then confirm each candidate with a graph
# comparison, because a hash alone is not proof of identity.
graphs = [smsd.parse_smiles(text) for text in smiles]
seen = {}
unique, repeats = [], []
for index, graph in enumerate(graphs):
    bucket = seen.setdefault(smsd.canonical_hash(graph), [])
    same = [kept for kept in bucket
            if smsd.same_canonical_graph(graphs[kept], graph)]
    if same:
        repeats.append((index, same[0]))
    else:
        bucket.append(index)
        unique.append(index)
print("records:", len(graphs), "| distinct structures:", len(unique))
print("first repeat:", names[repeats[0][0]], "is", names[repeats[0][1]])
assert len(unique) == 141 and len(repeats) == 392
assert [names[i] for i in repeats[0]] == ["GEN-008", "GEN-001"]

# Explicit hydrogens are graph atoms, so they change every fingerprint.
plain = smsd.parse_smiles(NAMED["paracetamol"])
explicit = smsd.parse_smiles(  # paracetamol with explicit hydrogens
    "[H]Oc1c([H])c([H])c(N([H])C(=O)C([H])([H])[H])c([H])c1[H]")
score = smsd.tanimoto_coefficient(smsd.fingerprint(plain, kind="ecfp"),
                                  smsd.fingerprint(explicit, kind="ecfp"))
print("paracetamol atoms:", plain.n, "and", explicit.n,
      "| ECFP4 Tanimoto:", f"{score:.3f}")
assert (plain.n, explicit.n) == (11, 20)
assert f"{score:.3f}" == "0.018"
```

The program prints:

```text
records: 533 | distinct structures: 141
first repeat: GEN-008 is GEN-001
paracetamol atoms: 11 and 20 | ECFP4 Tanimoto: 0.018
```

The same paracetamol with and without explicit hydrogens shares almost no ECFP
bits, which is why a library should be consistent before it is screened.

## Similarity screen in memory

`TargetCorpus.from_smiles` parses the library once and keeps it in memory.
`prewarm()` computes an ECFP radius 2, 2048-bit fingerprint for every target;
`screen` and `screen_many` then compare one or more queries with them and
return the indices of the targets whose Tanimoto score is at least the
threshold, in library order. They return indices only. For scores, use
`batch_similarity` on the hits, as in [Save the hits](#save-the-hits), or a
[packed index](PACKED_INDEX.md).

A screen run before `prewarm()`, or after any target is added, returns an empty
list rather than raising, so call `prewarm()` again after every addition.
`smsd.warm_start(corpus, query, threshold)` prepares the corpus when it needs
it and then screens.

```python
import itertools

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

diazepam = NAMED["diazepam"]
alprazolam = NAMED["alprazolam"]

corpus = smsd.TargetCorpus.from_smiles(smiles[:21])  # the named drugs
corpus.prewarm(num_threads=2)
corpus.add_targets(smiles[21:])                      # then the generated set
# An addition clears the preparation, and an unprepared screen finds nothing.
assert not corpus.is_prewarmed
assert corpus.screen(diazepam, 0.5, num_threads=2) == []

corpus.prewarm(num_threads=2)
hits = corpus.screen(diazepam, 0.5, num_threads=2)
print("diazepam, at least 0.5:", [names[i] for i in hits])

per_query = corpus.screen_many([diazepam, alprazolam], 0.5, num_threads=2)
print("alprazolam, at least 0.5:", [names[i] for i in per_query[1]])

# The threshold is inclusive: 1.0 keeps identical fingerprints.
exact = corpus.screen(diazepam, 1.0, num_threads=2)
print("diazepam, at least 1.0:", [names[i] for i in exact])

# warm_start prepares a corpus when needed, then screens it.
fresh = smsd.TargetCorpus.from_smiles(smiles)
assert smsd.warm_start(fresh, diazepam, 0.5, num_threads=2) == hits
assert fresh.is_prewarmed

assert [names[i] for i in hits] == [
    "diazepam", "nordazepam", "temazepam", "tetrazepam", "clobazam"]
assert per_query[0] == hits
assert [names[i] for i in per_query[1]] == [
    "alprazolam", "triazolam", "midazolam"]
assert exact == [0]
```

The program prints:

```text
diazepam, at least 0.5: ['diazepam', 'nordazepam', 'temazepam', 'tetrazepam', 'clobazam']
alprazolam, at least 0.5: ['alprazolam', 'triazolam', 'midazolam']
diazepam, at least 1.0: ['diazepam']
```

`screen` and `screen_many` do not check the threshold: a value above 1 or a NaN
returns no hits, and a negative value returns every target, so pass a number
from 0 to 1. `warm_start` does check it, but defaults to 0.0, which returns
every target; always pass a threshold. Thread counts must be non-negative
integers; `screen` and `prewarm` do not check this either.

## Streamed similarity screen

A library that is screened once does not need to sit in memory.
`smsd.stream_screen_smiles(query, targets, threshold)` takes any one-shot
iterable of molecules, such as a generator over the lines of a file, builds one
chunk of targets at a time, screens it with the same fixed ECFP radius 2,
2048-bit fingerprints and yields the input index of each hit, in input order.
Only one chunk is held at a time, and NumPy is not needed.

`smsd.stream_screen_many_smiles` screens several queries in bounded tiles of
queries by targets and yields `(query_index, target_index)` pairs. The pairs
come tile by tile, then by query, then by target: they are neither grouped by
query nor ranked, so sort them when order matters.

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533


def smiles_column(path):
    """Yield the first column of a SMILES file, one line at a time."""
    with open(path, encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                yield line.split()[0]


diazepam = NAMED["diazepam"]
alprazolam = NAMED["alprazolam"]
with TemporaryDirectory() as directory:
    path = Path(directory) / "library.smi"
    path.write_text("".join(f"{text} {name}\n" for name, text in library),
                    encoding="utf-8")

    hits = list(smsd.stream_screen_smiles(
        diazepam, smiles_column(path), 0.5, chunk_size=100, num_threads=2))
    print("diazepam, at least 0.5, input indices:", hits)

    pairs = list(smsd.stream_screen_many_smiles(
        [diazepam, alprazolam], smiles_column(path), 0.5,
        target_chunk_size=4, query_chunk_size=2, num_threads=2))
    print("(query, input index) as streamed:", pairs)
    print("(query, input index) sorted:", sorted(pairs))

assert [names[i] for i in hits] == [
    "diazepam", "nordazepam", "temazepam", "tetrazepam", "clobazam"]
assert pairs == [(0, 0), (0, 1), (0, 2), (0, 9), (1, 10), (1, 11), (0, 13),
                 (1, 12)]
assert [index for query, index in sorted(pairs) if query == 0] == hits
```

The program prints:

```text
diazepam, at least 0.5, input indices: [0, 1, 2, 9, 13]
(query, input index) as streamed: [(0, 0), (0, 1), (0, 2), (0, 9), (1, 10), (1, 11), (0, 13), (1, 12)]
(query, input index) sorted: [(0, 0), (0, 1), (0, 2), (0, 9), (0, 13), (1, 10), (1, 11), (1, 12)]
```

With tiles of four targets, the diazepam hit at index 13 arrives after the
alprazolam hits at indices 10 and 11, because it lies in a later tile; the
sorted pairs match the in-memory results. An invalid molecule later in the input
raises `ValueError` after the hits of earlier chunks have been yielded, and the
message does not name the failing input, so check a doubtful file before a
long run.

## Substructure screen with exact confirmation

`Substructure.match_corpus(corpus)` searches a `TargetCorpus` and returns one
`CorpusSearchResult` per target, in library order. Targets that cannot contain
the query are ruled out before matching, with a proof; every other target is
searched exactly on the CPU. The decisions are the same as those of
`match_many` over the whole library, and no `prewarm()` is needed.

| Field | Meaning |
| --- | --- |
| `target_index` | The position of the target in the corpus |
| `threshold_reached` | The query was found; `result` holds the match |
| `proven_below_threshold` | The query is proved absent, before matching or by a completed search |
| `result` | The `SearchResult` of the search, or `None` for a target ruled out before matching |
| `atom_upper_bound` | Always `None` for a substructure search |

When both flags are false the outcome is unknown, for example after a search
that stopped at its time limit; that is not proof of absence.

```python
import itertools

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

corpus = smsd.TargetCorpus.from_smiles(smiles)
# The 5-phenyl-1,4-benzodiazepin-2-one core
core = smsd.Substructure.find("O=C1CN=C(c2ccccc2)c2ccccc2N1")

decisions = core.match_corpus(corpus, num_threads=2)
found = [names[d.target_index] for d in decisions if d.threshold_reached]
ruled_out = sum(d.result is None for d in decisions)
absent = [names[d.target_index] for d in decisions
          if d.result is not None and d.proven_below_threshold]
unknown = sum(not d.threshold_reached and not d.proven_below_threshold
              for d in decisions)
print("contain the core:", found)
print("ruled out before matching:", ruled_out, "| searched and absent:",
      absent, "| unknown:", unknown)

# STRICT aromaticity: an aromatic query bond needs an aromatic target bond.
options = smsd.Substructure.Options()
options.chemistry = smsd.ChemOptions().with_aromaticity_mode(
    smsd.AromaticityMode.STRICT)
strict = core.with_options(options).match_corpus(corpus, num_threads=2)
strict_found = [names[d.target_index] for d in strict if d.threshold_reached]
print("STRICT drops:", sorted(set(found) - set(strict_found)))

# A SMARTS query searches every target.
pattern = smsd.Substructure.smarts("O=C1CN=C(c2ccccc2)c2ccccc2N1")
by_smarts = pattern.match_corpus(corpus, num_threads=2)
print("SMARTS: ruled out", sum(d.result is None for d in by_smarts),
      "| found", sum(d.threshold_reached for d in by_smarts))

assert len(found) == 9 and len(strict_found) == 8
assert set(found) - set(strict_found) == {"tetrazepam"}
assert "bromazepam" not in found and absent == ["bromazepam"]
assert ruled_out + 10 == len(decisions) and unknown == 0
assert all(d.result is not None for d in by_smarts)
assert [d.threshold_reached for d in by_smarts] == [
    d.threshold_reached for d in strict]
assert [r.found for r in core.match_many(corpus, num_threads=2)] == [
    d.threshold_reached for d in decisions]
```

The program prints:

```text
contain the core: ['diazepam', 'nordazepam', 'temazepam', 'oxazepam', 'lorazepam', 'clonazepam', 'nitrazepam', 'flunitrazepam', 'tetrazepam']
ruled out before matching: 523 | searched and absent: ['bromazepam'] | unknown: 0
STRICT drops: ['tetrazepam']
SMARTS: ruled out 0 | found 8
```

Bromazepam is searched and found not to contain the core, because its pyridine
nitrogen sits where the core has a phenyl carbon. Tetrazepam contains the core
under the default chemistry, because the default `FLEXIBLE` aromaticity lets an
aromatic query bond match a single or double bond, so its cyclohexene ring
stands in for the phenyl ring. `STRICT` aromaticity requires an aromatic bond
and drops it. The SMARTS form of the same pattern behaves like `STRICT` here,
because a lower-case `c` in SMARTS is an aromatic carbon.

SMARTS queries, and searches with a node limit, a cancellation token or an
immediate time limit, search every target without ruling any out first. ECFP
similarity bits cannot rule out a substructure, so do not filter a substructure
screen by fingerprint similarity. With `smsd-cuda` installed,
`gpu_prefilter=True` asks the GPU to rule targets out; the decisions are the
same, and the step falls back to the CPU on any failure (see
[Request the GPU for a corpus pre-screen](CPU_GPU.md#request-the-gpu-for-a-corpus-pre-screen)).

## MCS-similarity screen

MCS similarity here is the atom-count Tanimoto of the maximum common
substructure: `common / (query atoms + target atoms - common)`. Computing the
MCS for every target is the expensive step, so screen first.
`smsd.screen_targets(query, targets, threshold)` computes a fast upper bound on
that score for every target and returns the indices that reach the threshold.
Under the same chemistry options it never rejects a target that can reach the
threshold. It needs no `prewarm()`. The bound can keep targets that fall short,
so run the exact MCS on the survivors and keep the results that are optimal and
reach the threshold.

`MCS.match_corpus(corpus, minimum_atoms)` answers a related question for each
target: whether the MCS with the query has at least `minimum_atoms` atoms. It
returns the same `CorpusSearchResult` fields as the substructure screen, with
`atom_upper_bound` holding an upper bound on the MCS atom count where one was
computed.

```python
import itertools

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

corpus = smsd.TargetCorpus.from_smiles(smiles)
query = smsd.parse_smiles(NAMED["diazepam"])

# 1. A fast upper bound that never drops a target able to reach 0.8.
candidates = smsd.screen_targets(query, corpus, 0.8, num_threads=2)

# 2. The exact MCS on the survivors, scored as an atom-count Tanimoto.
search = smsd.MCS.find(query).with_timeout(1000)
results = search.match_many([smiles[i] for i in candidates], num_threads=2)
confirmed = []
for index, result in zip(candidates, results):
    target_atoms = smsd.parse_smiles(smiles[index]).n
    common = result.best.atom_count if result.found else 0
    score = common / (query.n + target_atoms - common)
    if result.optimal and score >= 0.8:
        confirmed.append((names[index], common, f"{score:.3f}"))
print("bound kept:", [names[i] for i in candidates])
print("exact MCS, at least 0.8:")
for name, common, score in confirmed:
    print(f"  {name:<12} {common} common atoms  {score}")

# 3. Or one atom-count decision per target: an MCS of at least 18 atoms.
decisions = smsd.MCS.find(query).with_timeout(1000).match_corpus(
    corpus, 18, num_threads=2)
print("MCS of 18 atoms or more:",
      [names[d.target_index] for d in decisions if d.threshold_reached])

# Default chemistry lets an aromatic query bond pair with a single or double
# bond, so the cyclohexene ring of tetrazepam stands in for the phenyl ring.
strict = smsd.MCS.Options()
strict.chemistry = smsd.ChemOptions().with_aromaticity_mode(
    smsd.AromaticityMode.STRICT)
strict.timeout_ms = 1000
result = smsd.MCS.find(query).with_options(strict).match(NAMED["tetrazepam"])
target_atoms = smsd.parse_smiles(NAMED["tetrazepam"]).n
common = result.best.atom_count
score = common / (query.n + target_atoms - common)
print("tetrazepam with STRICT aromaticity:", common, "common atoms",
      f"{score:.3f}")

assert result.optimal and (common, f"{score:.3f}") == (14, "0.538")
assert len(candidates) == 9 and len(confirmed) == 7
assert ("tetrazepam", 20, "1.000") in confirmed
assert sum(d.threshold_reached for d in decisions) == 12
assert not any(not d.threshold_reached and not d.proven_below_threshold
               for d in decisions)
```

The program prints:

```text
bound kept: ['diazepam', 'nordazepam', 'temazepam', 'oxazepam', 'lorazepam', 'clonazepam', 'tetrazepam', 'alprazolam', 'clobazam']
exact MCS, at least 0.8:
  diazepam     20 common atoms  1.000
  nordazepam   19 common atoms  0.950
  temazepam    20 common atoms  0.952
  oxazepam     19 common atoms  0.905
  lorazepam    19 common atoms  0.864
  tetrazepam   20 common atoms  1.000
  alprazolam   19 common atoms  0.826
MCS of 18 atoms or more: ['diazepam', 'nordazepam', 'temazepam', 'oxazepam', 'lorazepam', 'clonazepam', 'nitrazepam', 'flunitrazepam', 'tetrazepam', 'alprazolam', 'triazolam', 'midazolam']
tetrazepam with STRICT aromaticity: 14 common atoms 0.538
```

The bound kept nine targets, and the exact MCS confirmed seven: clonazepam and
clobazam fall short of 0.8. Tetrazepam reaches 1.000 for the same reason as in
the substructure screen; with `STRICT` aromaticity it shares 14 atoms with
diazepam. The guarantee of `screen_targets` holds for the atom-count Tanimoto
only, not for bond counts, MCES or weighted scores, and only when the exact
search uses the same `ChemOptions`. A result that is not optimal, because its
search stopped at the time limit, gives a lower bound on the common atoms.

## Save the hits

`MolGraph.name` and `MolGraph.properties` can be set on any graph, and
`smsd.write_sdf` writes the name as the record title and every property as an
SD field. `batch_similarity` uses the same ECFP radius 2, 2048-bit fingerprints
as `TargetCorpus.screen`, so it supplies the scores for the hits. Format scores
as strings with three decimal places.

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

diazepam = NAMED["diazepam"]
corpus = smsd.TargetCorpus.from_smiles(smiles)
corpus.prewarm(num_threads=2)
hits = corpus.screen(diazepam, 0.5, num_threads=2)  # indices, no scores

# batch_similarity uses the same ECFP4 bits and returns the scores.
scored = smsd.batch_similarity(diazepam, [smiles[i] for i in hits])
graphs = []
for position, score in scored:  # most similar first
    graph = smsd.parse_smiles(smiles[hits[position]])
    graph.name = names[hits[position]]
    graph.properties = {"TANIMOTO": f"{score:.3f}"}
    graphs.append(graph)

with TemporaryDirectory() as directory:
    path = Path(directory) / "hits.sdf"
    smsd.write_sdf(graphs, path)
    text = path.read_text(encoding="utf-8")
    saved = smsd.read_sdf(path)
for graph in saved:
    print(f"{graph.name:<12} TANIMOTO {graph.properties['TANIMOTO']}")

assert "> <TANIMOTO>" in text
assert hits == [index for index, score in
                sorted(smsd.batch_similarity(diazepam, smiles)) if score >= 0.5]
assert scored[-1][1] == 0.5  # clobazam: kept, as the threshold is inclusive
assert [graph.name for graph in saved] == [
    "diazepam", "temazepam", "tetrazepam", "nordazepam", "clobazam"]
assert all(float(graph.properties["TANIMOTO"]) >= 0.5 for graph in saved)
```

The program prints:

```text
diazepam     TANIMOTO 1.000
temazepam    TANIMOTO 0.582
tetrazepam   TANIMOTO 0.542
nordazepam   TANIMOTO 0.528
clobazam     TANIMOTO 0.500
```

Clobazam scores exactly 0.500 and is kept, because the threshold is inclusive.

## End-to-end example

The program below runs the similarity, substructure and MCS-similarity screens
over the library and writes the molecules that pass all three to an SD file,
with both scores as fields:

```python
import itertools
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

NAMED = {
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
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "salicylic acid": "OC(=O)c1ccccc1O",
    "paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "phenol": "Oc1ccccc1",
}
GROUPS = ["C", "CC", "OC", "N", "O", "F", "Cl", "C(=O)O"]
library = list(NAMED.items()) + [
    (f"GEN-{number:03d}", f"c1c({a})cc({b})cc1{c}")
    for number, (a, b, c) in enumerate(itertools.product(GROUPS, repeat=3))
]
names = [name for name, _ in library]
smiles = [text for _, text in library]
assert len(library) == 533

QUERY = NAMED["diazepam"]
CORE = "O=C1CN=C(c2ccccc2)c2ccccc2N1"  # benzodiazepinone core

# 1. Prepare the library once; keep the names beside it.
corpus = smsd.TargetCorpus.from_smiles(smiles)
corpus.prewarm(num_threads=2)

# 2. Similarity: ECFP4 Tanimoto of at least 0.5.
similar = corpus.screen(QUERY, 0.5, num_threads=2)

# 3. Substructure: exact decisions for the core.
decisions = smsd.Substructure.find(CORE).match_corpus(corpus, num_threads=2)
with_core = [d.target_index for d in decisions if d.threshold_reached]

# 4. MCS similarity: bound first, then the exact MCS on the survivors.
query = smsd.parse_smiles(QUERY)
candidates = smsd.screen_targets(query, corpus, 0.8, num_threads=2)
results = smsd.MCS.find(query).with_timeout(1000).match_many(
    [smiles[i] for i in candidates], num_threads=2)
mcs_scores = {}
for index, result in zip(candidates, results):
    common = result.best.atom_count if result.found else 0
    score = common / (query.n + smsd.parse_smiles(smiles[index]).n - common)
    if result.optimal and score >= 0.8:
        mcs_scores[index] = score

# 5. Write the molecules that pass all three screens, with their scores.
fingerprint_scores = dict(smsd.batch_similarity(QUERY, smiles))
selected = [i for i in similar if i in with_core and i in mcs_scores]
graphs = []
for index in selected:
    graph = smsd.parse_smiles(smiles[index])
    graph.name = names[index]
    graph.properties = {"TANIMOTO": f"{fingerprint_scores[index]:.3f}",
                        "MCS_TANIMOTO": f"{mcs_scores[index]:.3f}"}
    graphs.append(graph)
with TemporaryDirectory() as directory:
    path = Path(directory) / "hits.sdf"
    smsd.write_sdf(graphs, path)
    saved = smsd.read_sdf(path)

print("similar:", len(similar), "| core:", len(with_core),
      "| MCS:", len(mcs_scores), "| all three:", len(saved))
for graph in saved:
    print(f"  {graph.name:<12} TANIMOTO {graph.properties['TANIMOTO']}"
          f"  MCS_TANIMOTO {graph.properties['MCS_TANIMOTO']}")

counts = (len(similar), len(with_core), len(mcs_scores), len(saved))
assert counts == (5, 9, 7, 4)
assert [graph.name for graph in saved] == [names[i] for i in selected]
assert names[selected[0]] == "diazepam"
```

The program prints:

```text
similar: 5 | core: 9 | MCS: 7 | all three: 4
  diazepam     TANIMOTO 1.000  MCS_TANIMOTO 1.000
  nordazepam   TANIMOTO 0.528  MCS_TANIMOTO 0.950
  temazepam    TANIMOTO 0.582  MCS_TANIMOTO 0.952
  tetrazepam   TANIMOTO 0.542  MCS_TANIMOTO 1.000
```

Clobazam passes the fingerprint screen but has neither the core nor an MCS
similarity of 0.8; oxazepam and lorazepam have the core and a high MCS
similarity, but their fingerprints fall below 0.5. The
[library screening notebook](../../python/examples/library_screening.ipynb)
runs the same steps cell by cell, with a packed index and clustering.

## Options

### TargetCorpus

| Member | Returns | Meaning |
| --- | --- | --- |
| `TargetCorpus.from_smiles(targets)` | `TargetCorpus` | Captures SMILES strings, `MolGraph` objects or RDKit molecules from any iterable |
| `add_target(target)`, `add_targets(targets)` | `None` | Adds targets; the corpus then needs `prewarm()` again |
| `prewarm(*, num_threads=0)` | `None` | Computes the screening fingerprints and graph invariants |
| `is_prewarmed` | `bool` | Whether `prewarm()` has run since the last addition |
| `len(corpus)` | `int` | The number of targets |
| `screen(query, threshold, *, num_threads=0)` | `list[int]` | ECFP radius 2, 2048-bit Tanimoto hits, at least `threshold`, in target order |
| `screen_many(queries, threshold, *, num_threads=0)` | `list[list[int]]` | One hit list per query, in query order |
| `screen_targets_many(queries, threshold, chemistry=None, *, num_threads=0)` | `list[list[int]]` | The MCS-similarity bound for each query; no `prewarm()` needed |
| `smsd.warm_start(corpus, query, threshold=0.0, *, num_threads=0)` | `list[int]` | Prepares the corpus if needed, then screens; the threshold is checked |

A `TargetCorpus` cannot be saved or pickled. `num_threads=0` uses the OpenMP
default, which `OMP_NUM_THREADS` sets.

### Streamed screens

`stream_screen_smiles(query, targets, threshold, *, chunk_size=4096,
num_threads=0)` and `stream_screen_many_smiles(queries, targets, threshold, *,
target_chunk_size=4096, query_chunk_size=64, max_queries=4096, num_threads=0)`
take these arguments:

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `query`, `queries` | molecule, iterable of molecules | required | SMILES strings, `MolGraph` objects or RDKit molecules |
| `targets` | iterable | required | Read once; a generator is suitable |
| `threshold` | float | required | Inclusive; finite and from 0 to 1, otherwise `ValueError` |
| `chunk_size` | int | `4096` | Targets per chunk; values above 1,048,576 are reduced to 1,048,576 |
| `target_chunk_size` | int | `4096` | Targets per tile, from 1 to 1,048,576 |
| `query_chunk_size` | int | `64` | Queries per tile, from 1 to 1,048,576; a tile may hold at most 1,048,576 query and target pairs |
| `max_queries` | int | `4096` | Most queries accepted, from 1 to 1,048,576; more raise `ValueError` |
| `num_threads` | int | `0` | Threads per chunk; 0 uses the OpenMP default |

With no queries, `stream_screen_many_smiles` yields nothing and does not read
the targets.

### screen_targets

`smsd.screen_targets(query, targets, threshold, chemistry=None, *,
num_threads=0)` returns the indices of the targets whose bound reaches the
threshold, in target order.

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `query` | molecule | required | SMILES string, `MolGraph` or RDKit molecule |
| `targets` | list or `TargetCorpus` | required | A `TargetCorpus` is converted once and reused |
| `threshold` | float | required | Inclusive atom-count Tanimoto threshold; not checked |
| `chemistry` | `ChemOptions` or `None` | `None` | `None` means `ChemOptions()`; use the options of the exact search |
| `num_threads` | int | `0` | Non-negative; 0 uses the OpenMP default |

### batch_similarity

`smsd.batch_similarity(query, targets, fp_size=2048, radius=2)` returns
`(index, Tanimoto)` pairs, most similar first; equal scores keep their input
order.

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `query` | molecule | required | SMILES string or `MolGraph` |
| `targets` | list | required | Every target is fingerprinted on each call, so keep the list small |
| `fp_size` | int | `2048` | Fingerprint width in bits |
| `radius` | int | `2` | ECFP radius |

## Notes

- Hits do not depend on the thread count or the chunk sizes. Only the order of
  the pairs from `stream_screen_many_smiles` and
  `FingerprintCorpus.screen_many` depends on the tile sizes.
- `smsd.similarity(a, b)` is the MCS-similarity bound of `screen_targets`, the
  same value as `smsd.similarity_upper_bound`, not a fingerprint score.
- `smsd.prewarm(mol)` caches the invariants of a molecule that is searched many
  times, such as a query run against several libraries.
- A parse error in `TargetCorpus.from_smiles`, the streamed screens or
  `FingerprintCorpus.from_smiles` does not name the failing input.
- `smsd.to_smiles` raises `ValueError` for the few stereo records it cannot
  write; keep those records in SD form. `smsd.canonical_smiles` is a
  deprecated alias of `smsd.to_smiles`.
- A substructure or MCS screen needs the whole `TargetCorpus` in memory; only
  the fingerprint screens stream.
- Speed depends on the data, the hardware and the thread count. Time a
  representative sample on your own machine, and compare the results as well
  as the times.

## See also

- [Packed fingerprint index](PACKED_INDEX.md): build, reopen and search a saved
  index, with clustering.
- [CPU and GPU screening](CPU_GPU.md): which screens can use a GPU, and how to
  check.
- [Fingerprints and similarity](FINGERPRINTS.md): fingerprint families and
  similarity coefficients.
- [Substructure search](SUBSTRUCTURE.md#search-a-prepared-corpus) and
  [maximum common substructure](MCS.md#mcs-similarity).
- [Streamed native screening](../PACKED_CORPORA.md#streamed-native-screening)
  and [explicit corpus prefilters](../PYTHON.md#explicit-corpus-prefilters) in
  the reference pages.
- [Screening examples](../EXAMPLES.md#screening) and the
  [library screening notebook](../../python/examples/library_screening.ipynb).

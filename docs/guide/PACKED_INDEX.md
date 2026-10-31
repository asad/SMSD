# Packed fingerprint index

A packed fingerprint index is a file of binary fingerprints that is built once,
reopened in any later session and scanned exactly, without holding the library
in memory. `smsd.FingerprintCorpus` (also `smsd.corpus.FingerprintCorpus`)
builds and reads it; `smsd.butina` and `smsd.maxmin` cluster and select from
it. Every operation on the index, including opening it, needs NumPy, which the
`numpy` extra installs:

```bash
python -m pip install "smsd[numpy]==8.1.0"
```

The examples use the same inline library as [Screening a library](SCREENING.md):
14 benzodiazepines, 7 common drugs and 512 generated trisubstituted benzenes.
Each program writes its index in a temporary directory.

## Build an index

`FingerprintCorpus.from_smiles(path, molecules, options=...)` fingerprints
SMILES strings, `MolGraph` objects (for example from `smsd.read_sdf`) or RDKit
molecules in bounded batches and writes the file. Any iterable works, including
a generator over a large file. `FingerprintCorpus.create(path, fingerprints,
fp_size=...)` writes fingerprints that already exist as lists of set bits;
`from_fingerprints` is the same method under another name. Both return the new
index, already open.

The file is written in full under a temporary name in the destination
directory and then published, so a failed build leaves no file behind. An
existing file is kept, and the build raises `FileExistsError`, unless
`overwrite=True` replaces it. Without `overwrite`, the filesystem must support
hard links. There is no append or merge: to change an index, build it again.

The file records the width and the row count only, not the fingerprint
specification or any identifiers. Store those beside it, for example in a small
JSON sidecar file:

```python
import itertools
import json
import struct
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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)  # ECFP4
with TemporaryDirectory() as directory:
    index_path = Path(directory) / "library.fp"
    with smsd.FingerprintCorpus.from_smiles(
            index_path, smiles, options=options, num_threads=2) as corpus:
        print("rows:", len(corpus), "| bits:", corpus.fp_size)
        print("set bits in row 0:", len(corpus.fingerprint(0)))
        assert len(corpus.fingerprint(0)) == 42

    # The file holds bits only: keep the specification and the names with it.
    spec = {"family": "ecfp", "radius": options.radius, "size": options.size,
            "names": names}
    index_path.with_suffix(".json").write_text(json.dumps(spec),
                                               encoding="utf-8")
    size = index_path.stat().st_size
    header = struct.unpack("<8sQQ", index_path.read_bytes()[:24])
    print(f"file size: {size:,} bytes | header: {header}")

    # An existing file is kept unless overwrite=True.
    try:
        smsd.FingerprintCorpus.from_smiles(index_path, smiles[:21],
                                           options=options)
    except FileExistsError:
        print("library.fp already exists and is unchanged")
    with smsd.FingerprintCorpus.from_smiles(
            index_path, smiles[:21], options=options, overwrite=True) as corpus:
        print("replaced with", len(corpus), "rows")
        assert len(corpus) == 21

    # Bit lists from any source: MACCS keys, 166 bits wide.
    maccs = smsd.FingerprintOptions.maccs()
    keys = smsd.fingerprints(smiles[:21], options=maccs)
    with smsd.FingerprintCorpus.create(Path(directory) / "maccs.fp", keys,
                                       fp_size=maccs.size) as corpus:
        print("MACCS index:", len(corpus), "rows of", corpus.fp_size, "bits")
        assert (len(corpus), corpus.fp_size) == (21, 166)
        assert corpus.fingerprint(16) == keys[16]

assert size == 24 + 533 * 256
assert header == (b"SMSDFP1\x00", 2048, 533)
```

The program prints:

```text
rows: 533 | bits: 2048
set bits in row 0: 42
file size: 136,472 bytes | header: (b'SMSDFP1\x00', 2048, 533)
library.fp already exists and is unchanged
replaced with 21 rows
MACCS index: 21 rows of 166 bits
```

A 2048-bit row takes 256 bytes, so the 533-row file is 24 header bytes plus
533 rows of 256 bytes.

## Reopen an index

`FingerprintCorpus(path)` opens an existing file read-only. It checks the
header and the file length and maps the rows into memory without reading them
all. Use it as a context manager, or call `close()`; `len(corpus)` and
`corpus.fp_size` still work after closing. `corpus.fingerprint(i)` returns the
set bits of row `i`.

Queries must be made with the same specification as the index. A query made
another way, for example with the default kind of `smsd.fingerprint` or with
another width, still runs and returns scores, but they answer a different
question and nothing reports the mismatch. Rebuild the options from the
sidecar file:

```python
import itertools
import json
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
    index_path = Path(directory) / "library.fp"
    sidecar = Path(directory) / "library.json"

    # Session 1: build the index and record how it was built.
    build = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
    smsd.FingerprintCorpus.from_smiles(index_path, smiles,
                                       options=build).close()
    sidecar.write_text(json.dumps({"family": "ecfp", "radius": 2, "size": 2048,
                                   "names": names}), encoding="utf-8")

    # Session 2: rebuild the options from the sidecar, then open the index.
    spec = json.loads(sidecar.read_text(encoding="utf-8"))
    options = smsd.FingerprintOptions.ecfp(radius=spec["radius"],
                                           size=spec["size"])
    row_names = spec["names"]
    diazepam = NAMED["diazepam"]
    with smsd.FingerprintCorpus(index_path) as corpus:
        assert corpus.fp_size == spec["size"] and len(corpus) == len(row_names)
        assert corpus.fingerprint(9) == smsd.fingerprint(smiles[9],
                                                         options=options)
        query = smsd.fingerprint(diazepam, options=options)
        right = [row_names[row]
                 for row, _ in corpus.similarities(query, threshold=0.5)]
        # Queries made another way still run, but answer a different question.
        default_kind = smsd.fingerprint(diazepam)  # kind "mcs"
        narrow = smsd.fingerprint(diazepam, kind="ecfp", fp_size=1024)
        wrong = [len(list(corpus.similarities(bits, threshold=0.5)))
                 for bits in (default_kind, narrow)]
    print("rows:", len(corpus), "| bits:", corpus.fp_size)
    print("same specification:", right)
    print("default kind, 1024-bit ECFP:", wrong, "rows")

assert right == ["diazepam", "nordazepam", "temazepam", "tetrazepam",
                 "clobazam"]
assert wrong == [0, 0]
```

The program prints:

```text
rows: 533 | bits: 2048
same specification: ['diazepam', 'nordazepam', 'temazepam', 'tetrazepam', 'clobazam']
default kind, 1024-bit ECFP: [0, 0] rows
```

## Similarity search with a threshold

`corpus.similarities(query, threshold=...)` scans every row and yields
`(row, Tanimoto)` for the rows whose score is at least the threshold, in row
order. The scores are exactly those of `smsd.tanimoto_coefficient` on the same
bits. An empty query scores 1.0 against an empty row and 0.0 against any other.

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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
query = smsd.fingerprint(NAMED["diazepam"], options=options)
with TemporaryDirectory() as directory:
    with smsd.FingerprintCorpus.from_smiles(
            Path(directory) / "library.fp", smiles,
            options=options) as corpus:
        hits = list(corpus.similarities(query, threshold=0.5))
for row, score in hits:
    print(f"{row:3d}  {names[row]:<12} {score:.3f}")

# The scores are tanimoto_coefficient on the same bits ...
for row, score in hits:
    bits = smsd.fingerprint(smiles[row], options=options)
    assert score == smsd.tanimoto_coefficient(query, bits)
# ... and the rows are those of the in-memory ECFP4 screen.
targets = smsd.TargetCorpus.from_smiles(smiles)
targets.prewarm(num_threads=2)
assert [row for row, _ in hits] == targets.screen(NAMED["diazepam"], 0.5,
                                                  num_threads=2)
```

The program prints:

```text
  0  diazepam     1.000
  1  nordazepam   0.528
  2  temazepam    0.582
  9  tetrazepam   0.542
 13  clobazam     0.500
```

With ECFP radius 2 and 2048 bits, the rows are the same as those of
`TargetCorpus.screen` at the same threshold.

## The N most similar

There is no top-N method for a query from outside the index. With the default
threshold of 0, `similarities` yields every row, and `heapq.nlargest` keeps the
best N while the scan streams past. The key `(score, -row)` makes equal scores
prefer the lower row, so the result is repeatable.

```python
import heapq
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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
query = smsd.fingerprint(NAMED["oxazepam"], options=options)
with TemporaryDirectory() as directory:
    with smsd.FingerprintCorpus.from_smiles(
            Path(directory) / "library.fp", smiles, options=options) as corpus:
        # Every row with its score; equal scores prefer the lower row.
        best = heapq.nlargest(5, corpus.similarities(query),
                              key=lambda hit: (hit[1], -hit[0]))
for row, score in best:
    print(f"{names[row]:<12} {score:.3f}")

assert [names[row] for row, _ in best] == [
    "oxazepam", "lorazepam", "nordazepam", "temazepam", "diazepam"]
assert [score for _, score in best] == sorted(
    (score for _, score in best), reverse=True)
```

The program prints:

```text
oxazepam     1.000
lorazepam    0.755
nordazepam   0.558
temazepam    0.554
diazepam     0.377
```

## Several queries in one pass

`corpus.screen_many(queries, threshold=...)` scans the file once for up to 32
queries and yields `(query, row, Tanimoto)`. The output comes chunk by chunk,
then by query, then by row: it is neither grouped by query nor ranked. Sort it
when order matters. For more than 32 queries, call it again for each group of
32.

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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
queries = [smsd.fingerprint(NAMED[name], options=options)
           for name in ("diazepam", "alprazolam", "aspirin")]
with TemporaryDirectory() as directory:
    with smsd.FingerprintCorpus.from_smiles(
            Path(directory) / "library.fp", smiles, options=options) as corpus:
        # Small chunks show the order: by chunk, then query, then row.
        hits = list(corpus.screen_many(queries, threshold=0.5, chunk_size=4))
        whole = list(corpus.screen_many(queries, threshold=0.5))
        try:
            list(corpus.screen_many(queries * 11, threshold=0.5))
        except ValueError as error:
            message = str(error)
            print("33 queries:", message)
print("as yielded:", [(query, row) for query, row, _ in hits])
for query, row, score in sorted(hits):
    print(f"  query {query}  {names[row]:<12} {score:.3f}")

assert message == "screen_many supports at most 32 queries"
assert [(query, row) for query, row, _ in hits] == [
    (0, 0), (0, 1), (0, 2), (0, 9), (1, 10), (1, 11), (0, 13), (1, 12), (2, 14)]
assert sorted(hits) == sorted(whole)
```

The program prints:

```text
33 queries: screen_many supports at most 32 queries
as yielded: [(0, 0), (0, 1), (0, 2), (0, 9), (1, 10), (1, 11), (0, 13), (1, 12), (2, 14)]
  query 0  diazepam     1.000
  query 0  nordazepam   0.528
  query 0  temazepam    0.582
  query 0  tetrazepam   0.542
  query 0  clobazam     0.500
  query 1  alprazolam   1.000
  query 1  triazolam    0.760
  query 1  midazolam    0.508
  query 2  aspirin      1.000
```

With chunks of four rows, the diazepam hit in row 13 arrives after the
alprazolam hits in rows 10 and 11, because row 13 lies in a later chunk. With
the default chunk size, all 533 rows here form a single chunk.

## Neighbours and pairs within the index

`corpus.top_k_pairs(k)` yields the k nearest neighbours of every row as
`(row, neighbour, Tanimoto)`, grouped by row and best first; equal scores
prefer the lower neighbour. `corpus.threshold_pairs(threshold)` yields every
pair of rows whose score is at least the threshold, once each, as
`(first, second, Tanimoto)` with `first < second`. Both scan every pair
exactly, so the work grows with the square of the row count, and neither
builds a similarity matrix. Consume them as they stream: collecting every pair
can need memory that grows with the square of the row count, and collecting all
neighbours memory for n times k results.

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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    named = Path(directory) / "named.fp"
    with smsd.FingerprintCorpus.from_smiles(  # the 21 named drugs
            named, smiles[:21], options=options) as corpus:
        # The nearest neighbour of every row, best first, ties to the lower row.
        nearest = list(corpus.top_k_pairs(1))
        # Every pair at 0.5 or above, each once, first < second.
        pairs = list(corpus.threshold_pairs(0.5))
for query, target, score in nearest[:4]:
    print(f"{names[query]:<12} -> {names[target]:<12} {score:.3f}")
print("pairs at 0.5 or above:")
for first, second, score in pairs:
    print(f"  {names[first]:<12} {names[second]:<12} {score:.3f}")

assert len(nearest) == 21 and len(pairs) == 13
assert all(first < second for first, second, _ in pairs)
```

The program prints:

```text
diazepam     -> temazepam    0.582
nordazepam   -> oxazepam     0.558
temazepam    -> diazepam     0.582
oxazepam     -> lorazepam    0.755
pairs at 0.5 or above:
  diazepam     nordazepam   0.528
  diazepam     temazepam    0.582
  diazepam     tetrazepam   0.542
  diazepam     clobazam     0.500
  nordazepam   oxazepam     0.558
  nordazepam   nitrazepam   0.509
  nordazepam   bromazepam   0.519
  temazepam    oxazepam     0.554
  oxazepam     lorazepam    0.755
  clonazepam   nitrazepam   0.750
  alprazolam   triazolam    0.760
  alprazolam   midazolam    0.508
  triazolam    midazolam    0.541
```

## Clustering and diverse selection

`smsd.butina(corpus, threshold=0.7)` groups rows around cluster centres and
returns a list of tuples, each with its centre first and the other members in
row order. Every member reaches the threshold with its centre, but not
necessarily with the other members. `reordering=True` updates the neighbour
counts after each cluster is taken. `smsd.maxmin(corpus, count)` picks `count`
diverse rows: it starts at row `seed` (0 by default), then repeatedly takes the
row that is least similar to its closest chosen row. Both are deterministic,
and ties go to the lower row.

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

options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    named = Path(directory) / "named.fp"
    with smsd.FingerprintCorpus.from_smiles(  # the 21 named drugs
            named, smiles[:21], options=options) as corpus:
        clusters = smsd.butina(corpus, threshold=0.4)
        diverse = smsd.maxmin(corpus, 5)
print("Butina clusters at 0.4:", len(clusters))
for cluster in clusters:
    if len(cluster) > 1:
        print("  centre", names[cluster[0]], "with",
              [names[member] for member in cluster[1:]])
print("MaxMin, 5 diverse:", [names[row] for row in diverse])

assert sorted(row for cluster in clusters for row in cluster) == list(range(21))
assert len(clusters) == 13 and len(set(diverse)) == 5
assert [names[row] for row in diverse] == [
    "diazepam", "phenol", "caffeine", "naproxen", "aspirin"]
```

The program prints:

```text
Butina clusters at 0.4: 13
  centre nordazepam with ['diazepam', 'oxazepam', 'lorazepam', 'clonazepam', 'nitrazepam', 'bromazepam']
  centre alprazolam with ['triazolam', 'midazolam']
MaxMin, 5 diverse: ['diazepam', 'phenol', 'caffeine', 'naproxen', 'aspirin']
```

## Substructure searches

The index holds similarity fingerprints only, so it cannot answer a
substructure question. Load the library into a `TargetCorpus` and use
`Substructure.match_corpus`, which rules out targets that cannot contain the
query before matching and confirms every other target exactly; see
[Substructure screen with exact confirmation](SCREENING.md#substructure-screen-with-exact-confirmation).

## File format

Version 1 of the format is a 24-byte header followed by fixed-width rows:

| Offset | Size | Content |
| --- | --- | --- |
| 0 | 8 bytes | The magic bytes `SMSDFP1` followed by a zero byte |
| 8 | 8 bytes | `fp_size`, the width in bits, an unsigned little-endian integer from 1 to 16,777,216 |
| 16 | 8 bytes | The row count, an unsigned little-endian integer |
| 24 | rows times `ceil(fp_size / 8)` bytes | One row per fingerprint, in input order |

Bit `p` of a row is stored in byte `p // 8` of that row, at bit `p % 8`, so bit
0 is the least significant bit of the first byte. Unused bits in the last byte
of a row must be zero. The file length is exactly 24 bytes plus the rows. There
is no other metadata and no version field beyond the magic bytes.

Opening a file checks the magic bytes, the width and the length; searches also
reject non-zero padding bits. Any program can write the format. Within SMSD
Pro, only the Python package reads it. The program below writes a 12-bit,
three-row file by hand and opens it:

```python
import struct
from pathlib import Path
from tempfile import TemporaryDirectory

import smsd

fp_size = 12                      # any width from 1 to 16,777,216 bits
rows = [[0, 3, 11], [3], []]
row_bytes = (fp_size + 7) // 8    # 2 bytes per row
payload = bytearray()
for bits in rows:
    packed = bytearray(row_bytes)
    for bit in bits:
        packed[bit // 8] |= 1 << (bit % 8)  # bit 0: lowest bit of byte 0
    payload += packed

with TemporaryDirectory() as directory:
    path = Path(directory) / "by_hand.fp"
    header = struct.pack("<8sQQ", b"SMSDFP1\0", fp_size, len(rows))
    path.write_bytes(header + payload)
    data = path.read_bytes()
    print("bytes:", data.hex(" "))
    with smsd.FingerprintCorpus(path) as corpus:
        stored = [corpus.fingerprint(row) for row in range(len(corpus))]
        print("rows read back:", stored)
        hits = list(corpus.similarities([3], threshold=0.3))
        print("similar to [3]:",
              [(row, round(score, 3)) for row, score in hits])

    # Bits past fp_size in the last byte of a row must be zero.
    broken = Path(directory) / "padding.fp"
    broken.write_bytes(data[:25] + bytes([data[25] | 0x80]) + data[26:])
    with smsd.FingerprintCorpus(broken) as corpus:
        try:
            list(corpus.similarities([3]))
        except ValueError as error:
            message = str(error)
            print("ValueError:", message)

assert len(data) == 24 + len(rows) * row_bytes
assert stored == rows
assert [(row, round(score, 3)) for row, score in hits] == [(0, 0.333), (1, 1.0)]
assert message == "fingerprint corpus contains non-zero padding bits"
```

The program prints:

```text
bytes: 53 4d 53 44 46 50 31 00 0c 00 00 00 00 00 00 00 03 00 00 00 00 00 00 00 09 08 08 00 00 00
rows read back: [[0, 3, 11], [3], []]
similar to [3]: [(0, 0.333), (1, 1.0)]
ValueError: fingerprint corpus contains non-zero padding bits
```

## Build options

`FingerprintCorpus.from_smiles(path, molecules, *, kind="ecfp", fp_size=2048,
radius=2, path_length=7, options=None, chunk_size=4096, num_threads=0,
overwrite=False)` and `FingerprintCorpus.create(path, fingerprints, *,
fp_size=2048, chunk_size=4096, overwrite=False)` take these arguments:

| Name | Type | Default | Meaning |
| --- | --- | --- | --- |
| `path` | str or path | required | The file to create |
| `molecules` | iterable | required | `from_smiles`: SMILES strings, `MolGraph` objects or RDKit molecules |
| `fingerprints` | iterable of bit lists | required | `create`: set-bit positions, each from 0 to `fp_size - 1`; duplicates are folded, count pairs are rejected |
| `options` | `FingerprintOptions` or `None` | `None` | The full specification; when given, its `size` sets the width; count options are rejected |
| `kind` | str | `"ecfp"` | Used when `options` is `None`: `"path"`, `"mcs"`, `"ecfp"`, `"circular"`, `"fcfp"` or `"torsion"` |
| `fp_size` | int | `2048` | Width in bits, from 1 to 16,777,216 |
| `radius` | int | `2` | Circular radius, used when `options` is `None` |
| `path_length` | int | `7` | Path length, used when `options` is `None` |
| `chunk_size` | int | `4096` | Molecules or rows per batch; reduced when needed to keep each batch's arrays within about 64 MiB |
| `num_threads` | int | `0` | Fingerprinting threads; 0 uses the OpenMP default |
| `overwrite` | bool | `False` | `True` replaces an existing file; `False` raises `FileExistsError` |

The default kind here is `"ecfp"`, whereas `smsd.fingerprint` defaults to
`"mcs"`, so pass the same `FingerprintOptions` object to both.

## Search options

| Method | Arguments and defaults | Yields |
| --- | --- | --- |
| `similarities` | `query, *, threshold=0.0, chunk_size=4096, backend="cpu"` | `(row, Tanimoto)` in row order |
| `screen_many` | `queries, *, threshold=0.0, chunk_size=4096, backend="cpu"` | `(query, row, Tanimoto)` by chunk, then query, then row; at most 32 queries |
| `threshold_pairs` | `threshold, *, include_self=False, chunk_size=4096, backend="cpu"` | `(first, second, Tanimoto)` with `first < second` |
| `top_k_pairs` | `k, *, threshold=0.0, include_self=False, chunk_size=4096, backend="cpu"` | `(row, neighbour, Tanimoto)`, at most `k` per row, best first |
| `iter_packed_chunks` | `*, chunk_size=4096` | `(first_row, array)`: copies of the packed bytes as `uint8` arrays |
| `smsd.butina` | `corpus, threshold=0.7, *, reordering=False, chunk_size=4096, backend="cpu"` | Returns a list of tuples, centre first |
| `smsd.maxmin` | `corpus, count, *, seed=0, chunk_size=4096, backend="cpu"` | Returns a list of rows in selection order |

Every threshold is inclusive and must be a finite number from 0 to 1;
`threshold=True` raises `TypeError` and a value outside the range raises
`ValueError`. A query bit outside 0 to `fp_size - 1` raises `ValueError`.
`chunk_size` changes how many rows are compared at a time, and so the working
memory and the order of `screen_many` output, never the hits. `backend` is
`"cpu"` (the default), `"cuda"` or `"auto"`; any other value raises
`ValueError`. [Choose a backend for packed indexes](CPU_GPU.md#choose-a-backend-for-packed-indexes)
describes the GPU choices and `corpus.cuda_calls`.

## Notes

- The working arrays of a search or build are kept to about 64 MiB per chunk.
  Mapped pages, queries and the results you keep are additional; see
  [Ordering, work and memory](../PACKED_CORPORA.md#ordering-work-and-memory).
- Keep the file unchanged while it is open. Opening checks only the header and
  the length.
- After `close()`, reading a row raises `ValueError`, but `len(corpus)` still
  works. Arrays from `iter_packed_chunks` are copies and stay valid.
- `TargetCorpus`, `stream_screen_smiles`, `stream_screen_many_smiles` and
  `Substructure.match_corpus` need no NumPy.

## See also

- [Screening a library](SCREENING.md): every screening route, from a ranked
  list to an MCS-similarity screen.
- [CPU and GPU screening](CPU_GPU.md): GPU backends for the index and for
  clustering.
- [Fingerprints and similarity](FINGERPRINTS.md): fingerprint families and
  `FingerprintOptions`.
- [Packed fingerprint corpora](../PACKED_CORPORA.md): the reference contract,
  including [CPU and CUDA backends](../PACKED_CORPORA.md#cpu-and-cuda-backends).
- The [library screening notebook](../../python/examples/library_screening.ipynb)
  builds, reopens and clusters an index.

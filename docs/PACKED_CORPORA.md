# Packed fingerprint corpora

`smsd.corpus.FingerprintCorpus` stores binary fingerprints as fixed-width packed
rows in a read-only, memory-mapped file. NumPy is required for this API and for
`smsd.clustering`; importing the modules does not import NumPy.

For a task-based walk-through of building, reopening and searching an index,
see the [packed fingerprint index](guide/PACKED_INDEX.md) chapter of the user
guide.

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from smsd.corpus import FingerprintCorpus
from smsd.clustering import butina, maxmin

with TemporaryDirectory() as directory:
    path = Path(directory) / "library.fp"
    with FingerprintCorpus.create(path, [[1], [1], [2]], fp_size=8) as corpus:
        assert len(corpus) == 3
        assert corpus.fingerprint(0) == [1]
        assert list(corpus.similarities([1], threshold=1.0)) == [(0, 1.0), (1, 1.0)]
        assert list(corpus.screen_many([[1], [2]], threshold=1.0)) == [
            (0, 0, 1.0), (0, 1, 1.0), (1, 2, 1.0)
        ]
        assert list(corpus.threshold_pairs(1.0)) == [(0, 1, 1.0)]
        assert list(corpus.top_k_pairs(1)) == [
            (0, 1, 1.0), (1, 0, 1.0), (2, 0, 0.0)
        ]
        assert butina(corpus, threshold=1.0) == [(0, 1), (2,)]
        assert maxmin(corpus, count=3, seed=0) == [0, 2, 1]
        chunks = list(corpus.iter_packed_chunks(chunk_size=2))
        assert chunks[0][0] == 0
        assert chunks[0][1].tolist() == [[2], [2]]
        assert corpus.cuda_calls == 0  # The default backend is CPU.
    # Chunk arrays own their data and remain usable after the mapping is closed.
    assert chunks[1][1].tolist() == [[4]]
    with FingerprintCorpus(path) as reopened:
        assert reopened.fingerprint(2) == [2]
```

`from_fingerprints` is an alias for `create`. The fingerprint input can be a
generator. Construction publishes a complete file atomically, preserves an
existing destination by default and removes its temporary file after an input
failure. Use `overwrite=True` for replacement. The containing filesystem must
support hard links for the default create mode.

To create packed rows directly from molecule inputs, use `from_smiles`. Queries
must use the same fingerprint specification:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
import smsd
from smsd.corpus import FingerprintCorpus

molecules = ["CCO", "CCN", "c1ccccc1", "CC(=O)O"]
options = smsd.FingerprintOptions.ecfp(radius=2, size=2048)
with TemporaryDirectory() as directory:
    with FingerprintCorpus.from_smiles(
        Path(directory) / "molecules.fp", iter(molecules), options=options,
        chunk_size=1024
    ) as corpus:
        query = smsd.fingerprint("CCO", options=options)
        for index, score in corpus.similarities(query, threshold=0.3):
            print(index, score)
```

The file stores its bit width and row count, rather than the chemistry
specification. Keep that specification with the application. Count fingerprints
are rejected; duplicate binary positions are folded. Two empty fingerprints
score 1, and an empty fingerprint compared with a non-empty one scores 0, as
with `smsd.tanimoto_coefficient`. Thresholds are inclusive and must be finite
values from 0 to 1.

## CPU and CUDA backends

`screen_many`, `similarities`, `threshold_pairs`, `top_k_pairs`, `butina` and
`maxmin` accept the same explicit `backend` argument:

| Backend | Behaviour |
| --- | --- |
| `"cpu"` | Default. NumPy computes packed intersections; no CUDA provider is requested. |
| `"cuda"` | Requires a compatible `smsd-cuda` package and device. Unavailable providers or unsupported widths raise. |
| `"auto"` | Uses the CUDA provider when available and supported; otherwise uses CPU comparisons. |

The core and CUDA package versions must match. Packed CUDA supports fingerprint
widths up to 65,536 bits and a little-endian host. Wider packed fingerprints
remain supported on the CPU; `auto` uses that CPU path. Errors during a device
call, result validation or resident cleanup propagate. They do not silently
retry the operation on the CPU. These packed operations do not use the native
`TargetCorpus` operation crossovers.

A CUDA-capable installation can request the device explicitly:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from smsd.corpus import FingerprintCorpus
from smsd.clustering import butina, maxmin

with TemporaryDirectory() as directory:
    with FingerprintCorpus.create(
        Path(directory) / "bits.fp", [[1], [1, 2], [3]], fp_size=8
    ) as corpus:
        before = corpus.cuda_calls
        hits = list(corpus.screen_many([[1], [3]], threshold=0.5, backend="cuda"))
        assert hits == [(0, 0, 1.0), (0, 1, 0.5), (1, 2, 1.0)]
        assert corpus.cuda_calls > before
        pairs = list(corpus.threshold_pairs(0.5, backend="cuda"))
        neighbours = list(corpus.top_k_pairs(1, backend="cuda"))
        clusters = butina(corpus, threshold=0.5, backend="cuda")
        selected = maxmin(corpus, count=2, backend="cuda")
```

CUDA computes exact integer intersections. `screen_many` applies integer cutoff
tables on device, chosen to reproduce the CPU's binary64 threshold comparison;
the host forms the returned Tanimoto scores from integer counts. `top_k_pairs`
performs threshold selection and per-tile top-k on device, then merges those
lists on the CPU with the same descending score and ascending target-index tie
break. Pair searches and clustering use the same exact scores. GPU acceleration
does not make those algorithms approximate.

Butina centre selection, count updates and assignment bookkeeping remain on the
CPU. MaxMin's minimum-distance updates and next-row choice also remain on the
CPU. Both have deterministic lower-index tie breaks. The `cuda_calls` property
is a cumulative count of completed packed CUDA calls, including clustering
passes. A positive increase is evidence that this adapter completed CUDA work;
it is not a count of queries, hits or GPU kernels, and does not establish a
speedup or the fraction of an application executed on the GPU.

## Ordering, work and memory

`screen_many` accepts at most 32 packed queries and consumes at most one extra
query input to detect an excess before reading targets. It yields target-chunk
order, then query index and target index. This is not a global query-major or
similarity ranking. `similarities` yields target input order for one query.

CUDA tiles contain at most 1,048,576 query/target pairs. A conservative 64 MiB
host-work estimate also limits each tile, so the effective target chunk can be
smaller than the requested `chunk_size`. CPU packed comparison temporaries are
similarly bounded to approximately 64 MiB; `chunk_size` can lower the bound.
These are array/work estimates, not hard caps on process memory. Mapped pages,
query storage, Python objects and application results are additional costs.
Molecule construction and parsing memory are also input-dependent.

The native resident workspace separately retains at most 256 MiB of device
scratch. That cap is not a total GPU-memory limit: resident target fingerprints,
other temporary allocations, CUDA runtime resources and host pinned staging are
separate. The packed adapter creates and releases a target resident handle for
each tile; an open memory-mapped file does not imply that the entire corpus is
resident on the GPU or that targets persist on device between calls. The native
`TargetCorpus` API has its own retained-library lifetime.

`top_k_pairs` yields directed neighbours, ordered by descending similarity and
then target index, with at most `k` heap entries in the global merge.
`threshold_pairs` yields each undirected pair once in index order. Both exclude
self pairs by default. Neither creates a dense all-against-all matrix. Retaining
every yielded threshold pair can nevertheless require quadratic application
memory, and retaining all top-k results can require O(n * k) memory.

The exact all-against-all scans take quadratic comparison work. Butina keeps
linear auxiliary arrays and returns tuples with the centre first. Its default
uses initial neighbour counts; `reordering=True` updates those counts after
each cluster. A member must meet the similarity threshold with its centre,
but need not meet it with every other member. MaxMin keeps a linear array of
minimum Tanimoto distances and uses a deterministic starting row; comparison
work is O(n * selected_count).

The format has a 24-byte little-endian header and `ceil(fp_size / 8)` bytes per
row. The supported width is 1 to 16,777,216 bits. A 2048-bit corpus requires
256 bytes per row plus the header. The format and bounded tiles do not establish
validation at a particular corpus size or any throughput result.

Keep an open corpus file unchanged, and close it with a context manager.
Concurrent reads copy their mapped rows under a shared lifetime lock. `close()`
waits for an active copy before unmapping; repeated closes are harmless.
Computation on copied rows and already yielded results remains valid. A running
iterator raises `ValueError` when it next requests a mapped row after close.
Closing does not cancel a whole search. Explicit lifecycle calls are available
when a context manager is unsuitable:

```python
from pathlib import Path
from tempfile import TemporaryDirectory
from smsd.corpus import FingerprintCorpus

with TemporaryDirectory() as directory:
    corpus = FingerprintCorpus.create(Path(directory) / "bits.fp", [[1]], fp_size=8)
    try:
        assert corpus.fingerprint(0) == [1]
    finally:
        corpus.close()
    corpus.close()
```

## Streamed native screening

`stream_screen_smiles` screens a target iterable with the public native corpus
API and releases each target chunk before loading the next. It uses the fixed
ECFP radius 2, 2048-bit specification of `TargetCorpus.screen`. It requires no
NumPy installation. This native graph path is separate from the packed backend
choices above.

```python
from smsd.corpus import stream_screen_smiles

targets = iter(["CCO", "CCN", "c1ccccc1", "CCO"])
for index in stream_screen_smiles("CCO", targets, 0.5, chunk_size=4096):
    print(index)
```

Each chunk holds at most `chunk_size` target inputs, their prepared native
graphs and one chunk's hit indices; the maximum accepted chunk is 1,048,576
targets. Graph memory depends on the input molecules. Hit indices retain their
global input offsets. An invalid later molecule raises after any earlier
chunks' hits have been yielded.

The native screen uses the optional `smsd-cuda` backend when eligible and
retains its ordinary CPU fallback. The operation's crossover applies to each
chunk:
a small chunk may run on the CPU even when a GPU is available. This wrapper
does not force GPU execution. For repeated queries over one retained library,
use `TargetCorpus` directly so its resident device allocation can be reused.

For several query inputs, `stream_screen_many_smiles` uses bounded query and
target tiles through the same native operation:

```python
from smsd.corpus import stream_screen_many_smiles

molecules = ["CCO", "CCN", "c1ccccc1", "CC(=O)O"]
for query_index, target_index in stream_screen_many_smiles(
    iter(["CCO", "CCN"]), iter(molecules), 0.5,
    query_chunk_size=64, target_chunk_size=4096, max_queries=4096
):
    print(query_index, target_index)
```

The query inventory is limited by `max_queries`; one extra query input may be
consumed to reject an excess before reading targets. Queries are prepared once,
and the target iterable is consumed once. Native result storage covers at most
one query tile by one target tile. All three controls must be positive and at
most 1,048,576, and the query/target tile product cannot exceed 1,048,576.
These are bounds on molecule counts and result pairs; graph memory depends on
the input molecules. Hits follow target-tile order, then query index, then
target index. They are not globally query-major or ranked by similarity.
The packed corpus's 32-query limit does not apply to this native tiled wrapper.
Empty queries yield no hits without consuming targets. Invalid later targets
can raise after earlier tiles have yielded hits.

`smsd.warm_start(corpus, query, threshold=0.5)` discovers available accelerators,
prewarms that real corpus if needed and returns its ordinary screen hit indices.
The corpus retains any prepared fingerprints and resident device allocation.
The same dispatch rules apply; a successful call alone does not establish that
a GPU kernel ran.

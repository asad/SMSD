# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Deterministic Butina clustering and MaxMin selection of binary fingerprints.

Both functions default to ``backend="cpu"``. ``"cuda"`` uses exact device
intersections, while the CPU still chooses centres or selections and updates
bookkeeping arrays. ``"auto"`` may use CPU comparisons if the CUDA provider or
fingerprint width is unsupported; device-call failures are reported. Calls
completed by the CUDA adapter increase the input corpus's ``cuda_calls``.

>>> import tempfile
>>> from pathlib import Path
>>> with tempfile.TemporaryDirectory() as directory:
...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1], [2]], fp_size=8) as corpus:
...         butina(corpus, threshold=1.0, backend="cpu")
[(0, 1), (2,)]
"""

from .corpus import FingerprintCorpus, _chunk_size, _integer, _numpy, _threshold

__all__ = ["butina", "maxmin"]


def _checked_corpus(corpus, chunk_size, backend):
    """Validate the shared corpus and bounded scan controls."""
    from ._packed_accelerator import _backend
    _backend(backend)
    if not isinstance(corpus, FingerprintCorpus):
        raise TypeError("corpus must be a FingerprintCorpus")
    corpus._check_open()
    _chunk_size(chunk_size, corpus._row_bytes)


def butina(corpus, threshold=0.7, *, reordering=False, chunk_size=4096, backend="cpu"):
    """Cluster by a binary Tanimoto similarity threshold (Butina, 1999).

    Return tuples with the centre first and other members in input order. The
    centre has the most neighbours; equal counts prefer the lower input index.
    ``reordering=False`` uses initial neighbour counts. With ``reordering=True``,
    counts are updated after each cluster to include only unassigned neighbours.

    Exact comparisons use bounded chunks and O(n) auxiliary arrays, without a
    distance matrix or adjacency lists. The full comparison work is quadratic.
    A cluster is defined by proximity to its centre; its other members need not
    be within the threshold of each other.

    ``backend="cpu"`` is the default. CUDA computes exact integer intersections
    for widths up to 65,536 bits; centre selection, threshold comparisons and
    count updates stay on the CPU. Auto may fall back only for unavailable or
    unsupported providers; device-call failures propagate. ``corpus.cuda_calls``
    counts completed CUDA calls, including all clustering passes. Packed CUDA
    tiles obey the 1,048,576-pair and estimated 64 MiB host-work bounds; the O(n)
    bookkeeping and returned clusters are additional memory.

    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1], [2]], fp_size=8) as corpus:
    ...         butina(corpus, threshold=1.0, reordering=True, backend="cpu")
    [(0, 1), (2,)]
    """
    _checked_corpus(corpus, chunk_size, backend)
    threshold = _threshold(threshold)
    if not isinstance(reordering, bool):
        raise TypeError("reordering must be a bool")
    np = _numpy()
    counts = np.zeros(len(corpus), dtype=np.int64)
    active = np.ones(len(corpus), dtype=bool)
    for first in range(len(corpus)):
        query = corpus._packed_row(first)
        for start, scores in corpus._score_chunks(query, start=first + 1, chunk_size=chunk_size, backend=backend):
            neighbours = np.flatnonzero(scores >= threshold)
            counts[first] += len(neighbours)
            counts[start + neighbours] += 1
    clusters = []
    remaining = len(corpus)
    while remaining:
        centre = int(np.argmax(counts))
        members = [centre]
        query = corpus._packed_row(centre)
        for start, scores in corpus._score_chunks(query, chunk_size=chunk_size, backend=backend):
            for local in np.flatnonzero((scores >= threshold) & active[start:start + len(scores)]):
                index = start + int(local)
                if index != centre:
                    members.append(index)
        active[members] = False
        counts[members] = -1
        remaining -= len(members)
        clusters.append(tuple(members))
        if reordering and remaining:
            for removed in members:
                query = corpus._packed_row(removed)
                for start, scores in corpus._score_chunks(query, chunk_size=chunk_size, backend=backend):
                    neighbours = np.flatnonzero((scores >= threshold) & active[start:start + len(scores)])
                    counts[start + neighbours] -= 1
    return clusters


def maxmin(corpus, count, *, seed=0, chunk_size=4096, backend="cpu"):
    """Select diverse row indices by greedy MaxMin Tanimoto distance.

    Start at ``seed``, then maximise the minimum distance ``1 - Tanimoto`` to
    selected rows. Ties prefer the smaller input index. Return selection order;
    no randomness is used. Memory is O(n) plus one bounded comparison chunk,
    and comparison work is O(n * count). ``count=0`` returns an empty list.

    ``backend="cpu"`` is the default. CUDA computes exact integer intersections
    for widths up to 65,536 bits, with the same binary64 scores. Minimum-distance
    updates and the deterministic next-row choice remain CPU bookkeeping. Auto
    may use CPU comparisons for an unavailable provider or unsupported width;
    errors after device execution starts propagate. Successful CUDA calls
    increment ``corpus.cuda_calls``. Packed CUDA tiles obey the 1,048,576-pair
    and estimated 64 MiB host-work bounds, separate from the O(n) arrays.

    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1], [2]], fp_size=8) as corpus:
    ...         maxmin(corpus, count=3, seed=0, backend="cpu")
    [0, 2, 1]
    """
    _checked_corpus(corpus, chunk_size, backend)
    count = _integer(count, "count", 0, len(corpus))
    seed = _integer(seed, "seed", 0)
    if not count:
        if seed and seed >= len(corpus):
            raise IndexError("seed index out of range")
        return []
    seed = corpus._index(seed)
    np = _numpy()
    distances = np.ones(len(corpus), dtype=np.float64)
    selected = []
    next_index = seed
    for _ in range(count):
        selected.append(next_index)
        query = corpus._packed_row(next_index)
        for start, scores in corpus._score_chunks(query, chunk_size=chunk_size, backend=backend):
            np.minimum(distances[start:start + len(scores)], 1.0 - scores,
                       out=distances[start:start + len(scores)])
        distances[selected] = -1.0
        next_index = int(np.argmax(distances))
    return selected

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Read-only, memory-mapped binary fingerprints and streamed Tanimoto searches.

NumPy is imported when a packed corpus is opened or created. Fingerprints are
binary sets of bit positions, with the same empty-set convention as
``smsd.tanimoto_coefficient``. Packed searches default to ``backend="cpu"``;
``"cuda"`` requires the compatible CUDA package, and ``"auto"`` uses the CPU
when that provider is unavailable or the fingerprint width is unsupported.
Device-call failures are reported rather than retried on the CPU.

>>> import tempfile
>>> from pathlib import Path
>>> with tempfile.TemporaryDirectory() as directory:
...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [2]], fp_size=8) as corpus:
...         list(corpus.similarities([1], threshold=1.0, backend="cpu"))
[(0, 1.0)]
"""

import heapq
import math
import numbers
import operator
import os
import struct
import threading
from itertools import islice

__all__ = ["FingerprintCorpus", "stream_screen_smiles", "stream_screen_many_smiles"]

_HEADER = struct.Struct("<8sQQ")
_MAGIC = b"SMSDFP1\0"
_MAX_FP_SIZE = 1 << 24
_WORKING_BYTES = 64 * 1024 * 1024
_MAX_NATIVE_INPUTS = 1 << 20
_MAX_NATIVE_PAIRS = 1 << 20


def _numpy():
    """Load the optional array dependency only when it is needed."""
    try:
        import numpy
    except ImportError as error:
        raise ImportError("FingerprintCorpus requires NumPy; install numpy>=1.23") from error
    return numpy


def _integer(value, name, minimum=0, maximum=None):
    """Validate integral controls without accepting booleans or rounding floats."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise TypeError(f"{name} must be an integer")
    try:
        result = operator.index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if result < minimum or (maximum is not None and result > maximum):
        limit = f" between {minimum} and {maximum}" if maximum is not None else f" at least {minimum}"
        raise ValueError(f"{name} must be{limit}")
    return result


def _threshold(value):
    """A finite binary Tanimoto threshold."""
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError("threshold must be a real number")
    result = float(value)
    if not math.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ValueError("threshold must be finite and in [0, 1]")
    return result


def _chunk_size(value, row_bytes):
    """Bound packed comparison temporaries as well as the requested row count."""
    requested = _integer(value, "chunk_size", 1)
    return min(requested, max(1, _WORKING_BYTES // (3 * row_bytes + 64)))


def _pack(fingerprint, fp_size):
    """Pack a binary set, rejecting count pairs and positions outside its extent."""
    packed = bytearray((fp_size + 7) // 8)
    for bit in fingerprint:
        position = _integer(bit, "fingerprint position", 0, fp_size - 1)
        packed[position // 8] |= 1 << (position % 8)
    return packed


def stream_screen_smiles(query, targets, threshold, *, chunk_size=4096, num_threads=0):
    """Yield input indices passing an ECFP Tanimoto screen, one target chunk at a time.

    Inputs have the public ``TargetCorpus.from_smiles`` semantics. Each chunk
    is prewarmed, screened through ``TargetCorpus.screen_many`` and released.
    The native operation uses the optional smsd-cuda (CUDA) backend when
    eligible, including its operation-specific crossover and CPU fallback.
    Small chunks may therefore run on the CPU even when a GPU is available.

    Uses ECFP radius 2 and 2048 bits, as ``TargetCorpus.screen`` does. At most
    ``chunk_size`` input molecules, their graphs and one chunk's hits are held;
    molecular graph memory depends on the inputs. No NumPy dependency is needed.
    Earlier hits can have been yielded when an invalid later input raises.

    >>> list(stream_screen_smiles("CCO", iter(["CCO", "N", "CCO"]), 1.0, chunk_size=2))
    [0, 2]
    """
    from ._molecule import _ensure_mol
    from ._screening import TargetCorpus

    threshold = _threshold(threshold)
    chunk_size = _integer(chunk_size, "chunk_size", 1)
    chunk_size = min(chunk_size, 1 << 20)
    num_threads = _integer(num_threads, "num_threads", 0)
    query = _ensure_mol(query)
    targets = iter(targets)
    offset = 0
    while True:
        batch = list(islice(targets, chunk_size))
        if not batch:
            return
        corpus = TargetCorpus.from_smiles(batch)
        corpus.prewarm(num_threads=num_threads)
        hits = corpus.screen_many([query], threshold, num_threads=num_threads)[0]
        yield from (offset + index for index in hits)
        offset += len(batch)
        del corpus, batch, hits


def stream_screen_many_smiles(queries, targets, threshold, *, target_chunk_size=4096,
                              query_chunk_size=64, max_queries=4096, num_threads=0):
    """Yield ``(query_index, target_index)`` ECFP screen hits from bounded native tiles.

    Queries and targets may be one-shot iterables. At most ``max_queries``
    query graphs are retained; one extra query input is consumed to detect an
    excess before any targets are read. Each native call screens at most
    ``query_chunk_size`` by ``target_chunk_size`` pairs. Positive controls may
    not exceed 1,048,576 and their tile product may not exceed 1,048,576.

    Results follow target-tile order, then query index, then target index.
    This is not a global query-major or similarity ranking. At most one target
    tile's graphs and one result tile are retained, alongside the bounded query
    inventory. Molecular graph size remains input-dependent. Native CUDA
    eligibility, operation crossovers and CPU fallback apply to every tile.
    This native graph API is separate from the packed corpus's 32-query limit.

    >>> list(stream_screen_many_smiles(
    ...     iter(["CCO", "N"]), iter(["CCO", "N"]), 1.0,
    ...     query_chunk_size=1, target_chunk_size=1))
    [(0, 0), (1, 1)]
    """
    from ._molecule import _ensure_mol
    from ._screening import TargetCorpus

    threshold = _threshold(threshold)
    target_chunk_size = _integer(target_chunk_size, "target_chunk_size", 1, _MAX_NATIVE_INPUTS)
    query_chunk_size = _integer(query_chunk_size, "query_chunk_size", 1, _MAX_NATIVE_INPUTS)
    max_queries = _integer(max_queries, "max_queries", 1, _MAX_NATIVE_INPUTS)
    num_threads = _integer(num_threads, "num_threads", 0)
    if target_chunk_size * query_chunk_size > _MAX_NATIVE_PAIRS:
        raise ValueError("query/target tile product must not exceed 1048576")
    inputs = list(islice(iter(queries), max_queries + 1))
    if len(inputs) > max_queries:
        raise ValueError("query count exceeds max_queries")
    graphs = [_ensure_mol(query) for query in inputs]
    del inputs
    if not graphs:
        return
    targets = iter(targets)
    offset = 0
    while True:
        batch = list(islice(targets, target_chunk_size))
        if not batch:
            return
        corpus = TargetCorpus.from_smiles(batch)
        corpus.prewarm(num_threads=num_threads)
        for query_start in range(0, len(graphs), query_chunk_size):
            hits = corpus.screen_many(graphs[query_start:query_start + query_chunk_size], threshold,
                                      num_threads=num_threads)
            for local_query, target_indices in enumerate(hits):
                for local_target in target_indices:
                    yield query_start + local_query, offset + local_target
            del hits
        offset += len(batch)
        del corpus, batch


class FingerprintCorpus:
    """A read-only packed fingerprint file, indexed in input order.

    The version 1 format has a 24-byte little-endian header followed by fixed
    width rows, with bit zero in the least significant bit of each byte.
    ``fp_size`` is in [1, 16777216]. Opening validates the header and exact file
    length without reading all rows. The file must remain unchanged while open.

    Packed searches accept ``backend="cpu"`` (the default), ``"cuda"`` or
    ``"auto"``. ``"cuda"`` needs the optional smsd-cuda package of the same
    version, a usable NVIDIA device, a little-endian host and a width of at
    most 65,536 bits, and raises an error otherwise; a non-empty
    ``SMSD_DISABLE_GPU`` other than ``0`` makes the device unusable.
    ``"auto"`` uses CUDA under the same conditions and the CPU otherwise.
    Metal is never used. Once a device call starts, failures propagate. GPU work
    uses exact integer intersections; threshold decisions and returned binary64
    Tanimoto scores follow the CPU convention. Pair ordering, global top-k merges
    and clustering decisions remain deterministic CPU bookkeeping.

    Construction and packed comparisons use an approximately 64 MiB temporary
    array/host-work budget; smaller ``chunk_size`` values lower tile sizes. CUDA
    tiles contain at most 32 queries and 1,048,576 pairs. The native retained
    workspace's separate 256 MiB device-scratch cap is not a total-memory limit:
    resident targets, other device allocations, mapped pages, query arrays and
    returned data are additional. This adapter releases resident targets after
    each CUDA tile. The format alone establishes no molecule-count validation.

    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[], [1]], fp_size=8) as corpus:
    ...         len(corpus), corpus.fingerprint(1)
    (2, [1])
    """

    def __init__(self, path):
        """Open an existing packed file without reading its complete row payload.

        Use a context manager to release its mapping. Opening the same path does
        not recreate the fingerprints or retain an earlier CUDA resident handle.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     path = Path(directory) / "bits.fp"
        ...     FingerprintCorpus.create(path, [[1]], fp_size=8).close()
        ...     with FingerprintCorpus(path) as corpus:
        ...         corpus.fingerprint(0)
        [1]
        """
        from pathlib import Path

        np = _numpy()
        self.path = Path(path)
        self._lock = threading.RLock()
        self._packed = None
        self._closed = True
        self._cuda_calls = 0
        with self.path.open("rb") as stream:
            header = stream.read(_HEADER.size)
            if len(header) != _HEADER.size:
                raise ValueError("truncated fingerprint corpus header")
            magic, fp_size, count = _HEADER.unpack(header)
            if magic != _MAGIC:
                raise ValueError("unsupported fingerprint corpus format")
            self.fp_size = _integer(fp_size, "fp_size", 1, _MAX_FP_SIZE)
            self._row_bytes = (self.fp_size + 7) // 8
            expected = _HEADER.size + count * self._row_bytes
            if os.fstat(stream.fileno()).st_size != expected:
                raise ValueError("fingerprint corpus size does not match its header")
            self._count = count
            if count:
                self._packed = np.memmap(stream, mode="r", dtype=np.uint8,
                                         offset=_HEADER.size, shape=(count, self._row_bytes))
            else:
                self._packed = np.empty((0, self._row_bytes), dtype=np.uint8)
                self._packed.flags.writeable = False
        self._closed = False

    @classmethod
    def create(cls, path, fingerprints, *, fp_size=2048, chunk_size=4096, overwrite=False):
        """Stream binary fingerprints to a new file and open it read-only.

        Duplicate positions are folded. Invalid positions and count pairs are
        rejected. An incomplete construction leaves no destination file. The
        default preserves an existing destination; ``overwrite=True`` replaces
        it atomically after successful construction. The containing filesystem
        must support hard links for the default create mode. ``from_fingerprints``
        is an alias with the same arguments and behaviour.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     path = Path(directory) / "bits.fp"
        ...     with FingerprintCorpus.from_fingerprints(path, iter([[2, 1, 2], []]), fp_size=8) as corpus:
        ...         corpus.fingerprint(0), corpus.fingerprint(1)
        ([1, 2], [])
        """
        from pathlib import Path
        import tempfile

        _numpy()
        fp_size = _integer(fp_size, "fp_size", 1, _MAX_FP_SIZE)
        row_count = _chunk_size(chunk_size, (fp_size + 7) // 8)
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a bool")
        destination = Path(path)
        if not overwrite and destination.exists():
            raise FileExistsError(destination)
        descriptor, temporary = tempfile.mkstemp(prefix=".smsd-fp-", dir=destination.parent)
        try:
            with os.fdopen(descriptor, "wb") as stream:
                stream.write(_HEADER.pack(_MAGIC, fp_size, 0))
                count = 0
                buffer = bytearray()
                for fingerprint in fingerprints:
                    buffer.extend(_pack(fingerprint, fp_size))
                    count += 1
                    if count % row_count == 0:
                        stream.write(buffer)
                        buffer.clear()
                stream.write(buffer)
                stream.seek(0)
                stream.write(_HEADER.pack(_MAGIC, fp_size, count))
                stream.flush()
                os.fsync(stream.fileno())
            if overwrite:
                os.replace(temporary, destination)
            else:
                # A hard link publishes a complete file without replacing a concurrent writer.
                os.link(temporary, destination)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
        # Open only after the temporary name is gone: Windows cannot remove a name of a mapped file.
        return cls(destination)

    from_fingerprints = create

    @classmethod
    def from_smiles(cls, path, molecules, *, kind="ecfp", fp_size=2048, radius=2,
                    path_length=7, options=None, chunk_size=4096, num_threads=0, overwrite=False):
        """Fingerprint molecule inputs in bounded batches, then write packed rows.

        The fingerprint arguments have the semantics of ``smsd.fingerprints``.
        Count fingerprints are rejected. All subsequent queries must use the
        same fingerprint specification; the file stores its bit width, not a
        chemistry specification. Molecule and native parsing memory is separate
        from the packed-array bound.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.from_smiles(Path(directory) / "molecules.fp", iter(["CCO", "CCN"])) as corpus:
        ...         len(corpus), corpus.fp_size
        (2, 2048)
        """
        from . import _smsd
        from ._fingerprints import fingerprints

        _integer(num_threads, "num_threads", 0)
        if options is not None:
            if not isinstance(options, _smsd.FingerprintOptions):
                raise TypeError("options must be a FingerprintOptions")
            if options.counts:
                raise ValueError("FingerprintCorpus stores binary fingerprints, not counts")
            fp_size = options.size
        fp_size = _integer(fp_size, "fp_size", 1, _MAX_FP_SIZE)
        batch_size = _chunk_size(chunk_size, (fp_size + 7) // 8)
        fingerprints([], kind=kind, fp_size=fp_size, radius=radius,
                     path_length=path_length, options=options, num_threads=num_threads)

        def rows():
            inputs = iter(molecules)
            while True:
                batch = list(islice(inputs, batch_size))
                if not batch:
                    return
                yield from fingerprints(batch, kind=kind, fp_size=fp_size, radius=radius,
                                        path_length=path_length, options=options, num_threads=num_threads)

        return cls.create(path, rows(), fp_size=fp_size, chunk_size=batch_size, overwrite=overwrite)

    def __len__(self):
        """Return the number of stored rows, including after the mapping is closed.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[], [1]], fp_size=8) as corpus:
        ...         len(corpus)
        2
        """
        return self._count

    def _check_open(self):
        """Reject access after the mapped file has been released."""
        with self._lock:
            if self._closed:
                raise ValueError("fingerprint corpus is closed")

    def _index(self, index):
        self._check_open()
        index = _integer(index, "index", 0)
        if index >= len(self):
            raise IndexError("fingerprint index out of range")
        return index

    def _packed_row(self, index):
        """Copy one row so iterators do not retain borrowed mapped views."""
        with self._lock:
            return self._packed[self._index(index)].copy()

    def _packed_chunk(self, start, stop):
        """Keep the mapping alive until the bounded read has its own storage."""
        with self._lock:
            self._check_open()
            return self._packed[start:stop].copy()

    def fingerprint(self, index):
        """Return the sorted set-bit positions of one row, as an owned list.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[3, 1, 3]], fp_size=8) as corpus:
        ...         corpus.fingerprint(0)
        [1, 3]
        """
        np = _numpy()
        row = self._packed_row(index)
        return np.flatnonzero(np.unpackbits(row, bitorder="little")[:self.fp_size]).tolist()

    def iter_packed_chunks(self, *, chunk_size=4096):
        """Yield ``(start_index, uint8_array)`` copies in input order.

        Each array owns its memory and remains usable after the corpus is closed.
        Its columns hold packed bytes, rather than one column per bit.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[0, 7], [1]], fp_size=8) as corpus:
        ...         chunks = list(corpus.iter_packed_chunks(chunk_size=1))
        ...     [(start, packed.tolist()) for start, packed in chunks]
        [(0, [[129]]), (1, [[2]])]
        """
        self._check_open()
        rows = _chunk_size(chunk_size, self._row_bytes)
        for start in range(0, len(self), rows):
            yield start, self._packed_chunk(start, start + rows)

    @property
    def cuda_calls(self):
        """Number of successfully completed packed CUDA calls made by this corpus.

        The cumulative count includes pair searches and clustering, not CPU work.
        One call may process a whole tile, so this is not a query, hit or kernel
        count. A positive increase establishes that this adapter completed CUDA
        work; it does not measure speed or the fraction of an application on GPU.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1]], fp_size=8) as corpus:
        ...         hits = list(corpus.similarities([1], backend="cpu"))
        ...         corpus.cuda_calls
        0
        """
        with self._lock:
            return self._cuda_calls

    def _record_cuda_call(self):
        with self._lock:
            self._cuda_calls += 1

    def _score_chunks(self, query, *, chunk_size=4096, start=0, backend="cpu"):
        """Yield score arrays for one packed query without a corpus-sized array."""
        np = _numpy()
        self._check_open()
        from ._packed_accelerator import _bit_counts, _cuda_backend, _cuda_rows
        provider = _cuda_backend(backend)
        if provider is not None and self.fp_size > 65536:
            if backend == "cuda":
                raise ValueError("CUDA packed screening supports at most 65536 bits")
            provider = None
        rows = _chunk_size(chunk_size, self._row_bytes)
        rows = _cuda_rows(self._row_bytes, self.fp_size, 1, rows) if provider is not None else rows
        if self.fp_size % 8 and int(query[-1]) >> (self.fp_size % 8):
            raise ValueError("fingerprint corpus contains non-zero padding bits")
        query_count = int(_bit_counts(query, np.uint64, np))
        for offset in range(start, len(self), rows):
            targets = self._packed_chunk(offset, offset + rows)
            if self.fp_size % 8 and np.any(targets[:, -1] >> (self.fp_size % 8)):
                raise ValueError("fingerprint corpus contains non-zero padding bits")
            if provider is None:
                intersection = _bit_counts(np.bitwise_and(targets, query), np.uint64, np)
                union = _bit_counts(targets, np.uint64, np) + query_count - intersection
                scores = np.ones(len(targets), dtype=np.float64)
                np.divide(intersection, union, out=scores, where=union != 0)
            else:
                scores = provider.scores(targets, query, self.fp_size)
                self._record_cuda_call()
            yield offset, scores

    def screen_many(self, queries, *, threshold=0.0, chunk_size=4096, backend="cpu"):
        """Yield ``(query_index, target_index, Tanimoto)`` from at most 32 queries.

        Consume at most 33 query inputs to reject excess before reading targets.
        Output follows target-chunk order, then query index and target index.
        Each CUDA tile has at most 1,048,576 pairs and an estimated 64 MiB host-work
        budget, which may reduce ``chunk_size``. CUDA supports widths up to 65,536.

        ``backend="cpu"`` is the default. ``"cuda"`` requires a compatible provider;
        ``"auto"`` permits CPU fallback for an unavailable provider or unsupported
        width. Device execution or validation errors propagate. CUDA computes
        integer intersections and threshold selection on device; host code forms
        the binary64 scores. Successful calls increment ``cuda_calls``. Returned
        hits are streamed; collecting all of them needs output-proportional memory.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[], [1]], fp_size=8) as corpus:
        ...         list(corpus.screen_many([[], [1]], threshold=1.0, backend="cpu"))
        [(0, 0, 1.0), (1, 1, 1.0)]
        """
        np = _numpy()
        from ._packed_accelerator import _backend, _bit_counts, _cuda_backend, _cuda_rows
        _backend(backend)
        threshold = _threshold(threshold)
        rows = _chunk_size(chunk_size, self._row_bytes)
        self._check_open()
        inputs = list(islice(iter(queries), 33))
        if len(inputs) > 32:
            raise ValueError("screen_many supports at most 32 queries")
        if not inputs:
            return
        packed = np.array([np.frombuffer(_pack(query, self.fp_size), dtype=np.uint8)
                           for query in inputs], dtype=np.uint8)
        del inputs
        if self.fp_size > 65536 and backend == "cuda":
            raise ValueError("CUDA packed screening supports at most 65536 bits")
        provider = None if self.fp_size > 65536 else _cuda_backend(backend)
        if provider is not None:
            rows = _cuda_rows(self._row_bytes, self.fp_size, len(packed), rows)
        query_counts = _bit_counts(packed, np.uint64, np)
        for start in range(0, len(self), rows):
            targets = self._packed_chunk(start, start + rows)
            if self.fp_size % 8 and np.any(targets[:, -1] >> (self.fp_size % 8)):
                raise ValueError("fingerprint corpus contains non-zero padding bits")
            if provider is not None:
                selected = provider.select(targets, packed, self.fp_size, threshold)
                self._record_cuda_call()
                for query_index, (indices, _, scores) in enumerate(selected):
                    for index, score in zip(indices, scores):
                        yield query_index, start + int(index), float(score)
            else:
                target_counts = _bit_counts(targets, np.uint64, np)
                for query_index, query in enumerate(packed):
                    common = _bit_counts(np.bitwise_and(targets, query), np.uint64, np)
                    union = target_counts + query_counts[query_index] - common
                    scores = np.ones(len(targets), dtype=np.float64)
                    np.divide(common, union, out=scores, where=union != 0)
                    for local in np.flatnonzero(scores >= threshold):
                        yield query_index, start + int(local), float(scores[local])

    def similarities(self, query, *, threshold=0.0, chunk_size=4096, backend="cpu"):
        """Yield ``(target_index, Tanimoto)`` hits in input order for a bit set.

        The threshold is inclusive. Backend selection and CUDA bounds are the
        same as ``screen_many``; CPU is the default. Empty versus empty scores 1,
        and empty versus non-empty scores 0.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[], [1]], fp_size=8) as corpus:
        ...         list(corpus.similarities([1], threshold=1.0, backend="cpu"))
        [(1, 1.0)]
        """
        for _, index, score in self.screen_many([query], threshold=threshold,
                                                chunk_size=chunk_size, backend=backend):
            yield index, score

    def threshold_pairs(self, threshold, *, include_self=False, chunk_size=4096, backend="cpu"):
        """Yield undirected ``(first, second, Tanimoto)`` pairs meeting the threshold.

        Pairs are in index order, with ``first < second`` unless self pairs are
        requested. Computation is quadratic; storage is a query and one bounded
        chunk. Consuming the generator into a list can require quadratic memory.

        CPU is the default backend. With CUDA, exact integer intersections are
        computed on device and inclusive binary64 threshold comparisons use the
        same scores on the host. Auto fallback and width limits follow the class
        contract; successful device calls increment ``cuda_calls``.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1, 2], []], fp_size=8) as corpus:
        ...         list(corpus.threshold_pairs(0.5, backend="cpu"))
        [(0, 1, 0.5)]
        """
        np = _numpy()
        from ._packed_accelerator import _backend
        _backend(backend)
        threshold = _threshold(threshold)
        if not isinstance(include_self, bool):
            raise TypeError("include_self must be a bool")
        _chunk_size(chunk_size, self._row_bytes)
        self._check_open()
        for first in range(len(self)):
            query = self._packed_row(first)
            for start, scores in self._score_chunks(query, chunk_size=chunk_size,
                                                   start=first if include_self else first + 1,
                                                   backend=backend):
                for local in np.flatnonzero(scores >= threshold):
                    yield first, start + int(local), float(scores[local])

    def top_k_pairs(self, k, *, threshold=0.0, include_self=False, chunk_size=4096, backend="cpu"):
        """Yield at most ``k`` directed neighbours per row, ranked by similarity.

        Output is ``(query_index, target_index, Tanimoto)``, grouped by query.
        Equal scores prefer the smaller target index. At most ``k`` heap entries
        and one comparison chunk are kept. Every candidate is compared; this is
        an exact quadratic scan rather than an approximate neighbour index.

        CPU is the default backend. CUDA performs integer intersections, threshold
        selection and per-tile top-k on device; the CPU merges those exact lists
        with the same global score/index tie break. It requires the CUDA package's
        top-k entry point and does not retry failed device calls on the CPU.
        Successful device calls increment ``cuda_calls``. Collecting all results
        adds up to O(n * k) application memory.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     with FingerprintCorpus.create(Path(directory) / "bits.fp", [[1], [1], [2]], fp_size=8) as corpus:
        ...         list(corpus.top_k_pairs(1, backend="cpu"))
        [(0, 1, 1.0), (1, 0, 1.0), (2, 0, 0.0)]
        """
        from ._packed_accelerator import _cuda_backend, _cuda_rows
        from ._packed_accelerator import _backend
        _backend(backend)
        k = _integer(k, "k", 0)
        threshold = _threshold(threshold)
        if not isinstance(include_self, bool):
            raise TypeError("include_self must be a bool")
        _chunk_size(chunk_size, self._row_bytes)
        self._check_open()
        k = min(k, len(self) if include_self else max(0, len(self) - 1))
        if not k:
            return
        if self.fp_size > 65536 and backend == "cuda":
            raise ValueError("CUDA packed screening supports at most 65536 bits")
        provider = None if self.fp_size > 65536 else _cuda_backend(backend)
        for query_index in range(len(self)):
            best = []
            query = self._packed_row(query_index)
            rows = _chunk_size(chunk_size, self._row_bytes)
            rows = _cuda_rows(self._row_bytes, self.fp_size, 1, rows) if provider is not None else rows
            if provider is None:
                chunks = ((start, enumerate(scores)) for start, scores in
                          self._score_chunks(query, chunk_size=rows))
            else:
                def ranked_chunks():
                    for start in range(0, len(self), rows):
                        targets = self._packed_chunk(start, start + rows)
                        extra = int(not include_self and start <= query_index < start + len(targets))
                        indices, _, scores = provider.select(targets, query, self.fp_size,
                                                             threshold, top_k=min(k + extra, len(targets)))[0]
                        self._record_cuda_call()
                        yield start, zip(indices, scores)
                chunks = ranked_chunks()
            for start, candidates in chunks:
                for local, score in candidates:
                    target_index = start + int(local)
                    if (not include_self and target_index == query_index) or score < threshold:
                        continue
                    candidate = (float(score), -target_index)
                    if len(best) < k:
                        heapq.heappush(best, candidate)
                    elif candidate > best[0]:
                        heapq.heapreplace(best, candidate)
            for score, negative_index in sorted(best, key=lambda item: (-item[0], -item[1])):
                yield query_index, -negative_index, score

    def close(self):
        """Wait for an active mapped copy, then release; repeated closes are harmless.

        Array computations and already copied results may finish after close.
        Iterators raise on their next mapped read once the corpus is closed.
        Closing does not cancel an entire search already using copied rows.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     corpus = FingerprintCorpus.create(Path(directory) / "bits.fp", [[1]], fp_size=8)
        ...     row = corpus.fingerprint(0)
        ...     corpus.close()
        ...     corpus.close()
        ...     row
        [1]
        """
        with self._lock:
            if not self._closed:
                self._closed = True
                if hasattr(self._packed, "_mmap"):
                    self._packed._mmap.close()
                self._packed = None

    def __enter__(self):
        """Return this open corpus for use in a ``with`` statement.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     corpus = FingerprintCorpus.create(Path(directory) / "bits.fp", [[1]], fp_size=8)
        ...     corpus.__enter__() is corpus
        ...     corpus.close()
        True
        """
        self._check_open()
        return self

    def __exit__(self, *exception):
        """Close the mapping on leaving a context, without suppressing exceptions.

        >>> import tempfile
        >>> from pathlib import Path
        >>> with tempfile.TemporaryDirectory() as directory:
        ...     corpus = FingerprintCorpus.create(Path(directory) / "bits.fp", [[1]], fp_size=8)
        ...     corpus.__exit__(None, None, None) is None
        True
        """
        self.close()

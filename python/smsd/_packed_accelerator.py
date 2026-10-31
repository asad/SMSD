# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Bounded packed-fingerprint calls to the optional CUDA package.

The caller owns mapped-file lifetime and result ordering. This module neither
imports NumPy nor initialises CUDA until a CUDA backend is requested.
"""

import ctypes
import os
import sys

_MAX_PAIRS = 1 << 20
_HOST_BYTES = 64 << 20
# NumPy 2 counts set bits with bitwise_count; older NumPy, or False here, uses a byte table.
_BITWISE_COUNT = True
_byte_bits = None


class _Hit(ctypes.Structure):
    _fields_ = [("target_index", ctypes.c_int32), ("common_bits", ctypes.c_int32)]


_Create = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int64, ctypes.c_int32,
    ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32),
    ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_size_t)
_Screen = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_int32,
    ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
    ctypes.POINTER(_Hit), ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64),
    ctypes.c_void_p, ctypes.c_size_t)
_Destroy = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t)


class _Table(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
        ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
        ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
        ("minimum_targets", ctypes.c_void_p), ("create_fingerprints", _Create),
        ("fingerprint_screen", _Screen), ("create_labels", ctypes.c_void_p),
        ("resident_label_overlaps", ctypes.c_void_p), ("destroy_resident", _Destroy)]


def _backend(value):
    if value not in ("cpu", "cuda", "auto"):
        raise ValueError("backend must be 'cpu', 'cuda' or 'auto'")
    return value


def _cuda_backend(value):
    """An optional provider, or None only for an explicitly allowed CPU path."""
    _backend(value)
    if value == "cpu":
        return None
    try:
        return _PackedCuda()
    except (ImportError, OSError, RuntimeError) as error:
        if value == "auto":
            return None
        raise RuntimeError(f"The requested CUDA packed-screening backend is unavailable: {error}") from error


def _cuda_rows(row_bytes, fp_size, query_count, requested):
    """Conservative per-call host-work and pair limits, independent of hit density."""
    words = (row_bytes + 7) // 8
    fixed = 4 * query_count * (fp_size + 1) + 24 * words * query_count
    per_row = 4 * words * 8 + 32 * query_count
    if fixed + per_row > _HOST_BYTES:
        raise ValueError("Packed CUDA query tile exceeds its host-work bound")
    return min(requested, _MAX_PAIRS // query_count, (_HOST_BYTES - fixed) // per_row)


def _bit_counts(rows, dtype, np):
    """Set bits of uint8 rows, summed along the last axis in ``dtype``.

    Both counting paths give the same integers in the same dtype.
    """
    global _byte_bits
    count = getattr(np, "bitwise_count", None) if _BITWISE_COUNT else None
    if count is None:
        if _byte_bits is None:
            _byte_bits = np.array([value.bit_count() for value in range(256)], dtype=np.uint8)
        return _byte_bits[rows].sum(axis=-1, dtype=dtype)
    if rows.dtype == np.uint8 and rows.size and rows.shape[-1] % 8 == 0 and rows.flags.c_contiguous:
        rows = rows.view(np.uint64)
    return count(rows).sum(axis=-1, dtype=dtype)


def _cutoffs(query_counts, maximum_target_count, threshold, np):
    """The first intersection whose binary64 score reaches the threshold."""
    rows = np.zeros((len(query_counts), maximum_target_count + 1), dtype=np.int32)
    if threshold == 0.0:
        return rows
    target_counts = np.arange(maximum_target_count + 1, dtype=np.int64)
    for index, query_count in enumerate(query_counts):
        lower = np.zeros(len(target_counts), dtype=np.int64)
        upper = np.minimum(int(query_count), target_counts) + 1
        while np.any(lower < upper):
            active = lower < upper
            middle = lower + (upper - lower) // 2
            union = int(query_count) + target_counts - middle
            scores = np.ones(len(union), dtype=np.float64)
            np.divide(middle, union, out=scores, where=union != 0)
            accepted = scores >= threshold
            upper = np.where(active & accepted, middle, upper)
            lower = np.where(active & ~accepted, middle + 1, lower)
        rows[index] = lower
    return rows


class _PackedCuda:
    """One loaded provider; resident target ownership is bounded to each call."""

    def __init__(self):
        if sys.byteorder != "little":
            raise RuntimeError("CUDA packed screening requires a little-endian host")
        if os.environ.get("SMSD_DISABLE_GPU", "") not in ("", "0"):
            raise RuntimeError("CUDA screening is disabled by SMSD_DISABLE_GPU")
        import smsd_cuda
        from . import __version__

        if getattr(smsd_cuda, "__smsd_version__", None) != __version__:
            raise RuntimeError("smsd and smsd-cuda must have the same version")
        if not smsd_cuda.is_available():
            raise RuntimeError(smsd_cuda.unavailable_reason())
        self.library = smsd_cuda._load()
        address = self.library.smsd_cuda_screen_accelerator()
        if not address:
            raise RuntimeError("The CUDA screening interface is null")
        self.table = _Table.from_address(address)
        if (self.table.abi_version != 1 or self.table.struct_size < ctypes.sizeof(_Table)
                or self.table.smsd_version.decode() != __version__):
            raise RuntimeError("The CUDA resident-screening interface is incompatible")
        self.top_k = getattr(self.library, "smsd_cuda_fingerprint_top_k", None)
        if self.top_k is not None:
            self.top_k.restype = ctypes.c_int
            self.top_k.argtypes = [ctypes.c_void_p, ctypes.c_int32,
                ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32),
                ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(_Hit), ctypes.c_uint64,
                ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p, ctypes.c_size_t]
        self.calls = 0

    def select(self, targets, queries, fp_size, threshold=0.0, top_k=None):
        """Return owned index/intersection/score arrays, one tuple per query."""
        import numpy as np

        targets = np.asarray(targets, dtype=np.uint8)
        queries = np.atleast_2d(np.asarray(queries, dtype=np.uint8))
        row_bytes = (fp_size + 7) // 8
        if (targets.ndim != 2 or queries.ndim != 2
                or targets.shape[1] != row_bytes or queries.shape[1] != row_bytes):
            raise ValueError("Packed fingerprint dimensions differ")
        target_count, query_count = len(targets), len(queries)
        if fp_size > 65536 or not 0 <= query_count <= 32:
            raise ValueError("CUDA packed calls support at most 65536 bits and 32 queries")
        if target_count * query_count > _MAX_PAIRS:
            raise ValueError("CUDA packed query/target tile exceeds 1048576 pairs")
        if fp_size % 8 and (np.any(targets[:, -1] >> (fp_size % 8))
                            or np.any(queries[:, -1] >> (fp_size % 8))):
            raise ValueError("fingerprint corpus contains non-zero padding bits")
        if top_k is not None and self.top_k is None:
            raise RuntimeError("The CUDA package has no on-device top-k operation")
        if not target_count or not query_count:
            return [(np.empty(0, np.int32), np.empty(0, np.int32), np.empty(0, np.float64))
                    for _ in range(query_count)]
        target_counts = _bit_counts(targets, np.int32, np)
        query_counts = _bit_counts(queries, np.int32, np)
        maximum_count = int(target_counts.max())
        words = (row_bytes + 7) // 8
        estimate = (4 * target_count * words * 8 + query_count * target_count * 32
                    + 4 * query_count * (maximum_count + 1) + 24 * words * query_count)
        if estimate > _HOST_BYTES:
            raise ValueError("Reduce the packed CUDA tile to fit its 64 MiB host-work bound")
        packed_targets = np.zeros((target_count, words * 8), dtype=np.uint8)
        packed_targets[:, :row_bytes] = targets
        word_major = np.ascontiguousarray(packed_targets.view("<u8").T)
        del packed_targets
        packed_queries = np.zeros((query_count, words * 8), dtype=np.uint8)
        packed_queries[:, :row_bytes] = queries
        query_words = packed_queries.view("<u8")
        cutoffs = _cutoffs(query_counts, maximum_count, threshold, np)
        capacity = query_count * (target_count if top_k is None else min(top_k, target_count))
        hits = (_Hit * capacity)()
        offsets = (ctypes.c_uint64 * (query_count + 1))()
        handle, error = ctypes.c_void_p(), ctypes.create_string_buffer(512)
        status = self.table.create_fingerprints(target_count, words,
            word_major.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
            target_counts.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            ctypes.byref(handle), error, len(error))
        if status != 0 or not handle.value:
            raise RuntimeError("CUDA resident creation failed: " + error.value.decode(errors="replace"))
        try:
            arguments = [handle, query_count,
                query_words.ctypes.data_as(ctypes.POINTER(ctypes.c_uint64)),
                cutoffs.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)), maximum_count + 1]
            if top_k is None:
                function = self.table.fingerprint_screen
            else:
                function = self.top_k
                arguments.append(min(top_k, target_count))
            status = function(*arguments, hits, capacity, offsets, error, len(error))
            if status != 0:
                raise RuntimeError("CUDA packed screening failed: " + error.value.decode(errors="replace"))
            self.calls += 1
            if offsets[0] != 0 or not all(0 <= offsets[i] <= offsets[i + 1] <= capacity
                                         for i in range(query_count)):
                raise RuntimeError("Invalid CUDA hit offsets")
            records = np.frombuffer(hits, dtype=np.int32).reshape(-1, 2)
            result = []
            for query_index in range(query_count):
                row = records[offsets[query_index]:offsets[query_index + 1]]
                indices, common = row[:, 0].copy(), row[:, 1].copy()
                if np.any(indices < 0) or np.any(indices >= target_count):
                    raise RuntimeError("Invalid CUDA target index")
                if np.any(common < 0) or np.any(common > np.minimum(query_counts[query_index], target_counts[indices])):
                    raise RuntimeError("Invalid CUDA integer intersection")
                union = target_counts[indices].astype(np.int64) + int(query_counts[query_index]) - common
                scores = np.ones(len(indices), dtype=np.float64)
                np.divide(common, union, out=scores, where=union != 0)
                if (np.any(scores < threshold) or len(np.unique(indices)) != len(indices)
                        or (top_k is None and np.any(indices[1:] <= indices[:-1]))):
                    raise RuntimeError("Invalid CUDA hit order or threshold")
                if top_k is not None and any((-scores[i], indices[i]) > (-scores[i + 1], indices[i + 1])
                                            for i in range(len(indices) - 1)):
                    raise RuntimeError("Invalid CUDA top-k order")
                result.append((indices, common, scores))
            return result
        finally:
            status = self.table.destroy_resident(handle, error, len(error))
            if status != 0 and sys.exc_info()[0] is None:
                raise RuntimeError("CUDA resident cleanup failed: " + error.value.decode(errors="replace"))

    def scores(self, targets, query, fp_size):
        """All target scores, with no ranking or CPU intersection fallback."""
        import numpy as np

        indices, _, scores = self.select(targets, query, fp_size)[0]
        if not np.array_equal(indices, np.arange(len(targets), dtype=np.int32)):
            raise RuntimeError("CUDA omitted a zero-threshold target")
        return scores

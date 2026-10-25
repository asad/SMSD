# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Integer resident screening through the actual CUDA callback, without CPU fallback."""
import ctypes
import math

import pytest

smsd_cuda = pytest.importorskip("smsd_cuda")
pytestmark = pytest.mark.skipif(not smsd_cuda.is_available(), reason=str(smsd_cuda.unavailable_reason()))


class Hit(ctypes.Structure):
    _fields_ = [("target_index", ctypes.c_int32), ("common_bits", ctypes.c_int32)]


Create = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int64, ctypes.c_int32,
                         ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32),
                         ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_size_t)
Screen = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_int32,
                         ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32,
                         ctypes.POINTER(Hit), ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64),
                         ctypes.c_void_p, ctypes.c_size_t)
Destroy = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t)


class Table(ctypes.Structure):
    _fields_ = [("abi_version", ctypes.c_uint32), ("struct_size", ctypes.c_uint32),
                ("smsd_version", ctypes.c_char_p), ("name", ctypes.c_char_p),
                ("device_info", ctypes.c_void_p), ("label_overlaps", ctypes.c_void_p),
                ("minimum_targets", ctypes.c_void_p), ("create_fingerprints", Create),
                ("fingerprint_screen", Screen), ("create_labels", ctypes.c_void_p),
                ("resident_label_overlaps", ctypes.c_void_p), ("destroy_resident", Destroy)]


def _score(common, query_bits, target_bits):
    union = query_bits + target_bits - common
    return common / union if union else 1.0


def _cutoff(query_bits, target_bits, threshold):
    lower, upper = 0, min(query_bits, target_bits) + 1
    while lower < upper:
        common = lower + (upper - lower) // 2
        if _score(common, query_bits, target_bits) >= threshold:
            upper = common
        else:
            lower = common + 1
    return lower


def _fingerprint(index, words):
    if index % 7 == 0:
        return [0] * words
    if words > 32:
        return [1 << ((index + word) % 64) for word in range(words)]
    patterns = [(1 << 64) - 1, 0xAAAAAAAAAAAAAAAA, 0x5555555555555555]
    return [patterns[index % 3] if index % 7 < 4 else 1 << ((index + word) % 64)
            for word in range(words)]


@pytest.mark.parametrize("words,query_count", [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
                                             (32, 5), (32, 9), (32, 65), (1024, 5), (1025, 5)])
def test_resident_integer_thresholds_and_order(words, query_count):
    table = Table.from_address(smsd_cuda._library.smsd_cuda_screen_accelerator())
    assert table.abi_version == 1 and table.struct_size >= ctypes.sizeof(Table)
    targets = [_fingerprint(index, words) for index in range(257)]
    queries = [_fingerprint(index, words) for index in range(query_count)]
    target_counts = [sum(word.bit_count() for word in target) for target in targets]
    query_counts = [sum(word.bit_count() for word in query) for query in queries]
    packed_targets = (ctypes.c_uint64 * (words * len(targets)))(
        *(targets[target][word] for word in range(words) for target in range(len(targets))))
    packed_queries = (ctypes.c_uint64 * (words * query_count))(*(word for query in queries for word in query))
    counts = (ctypes.c_int32 * len(targets))(*target_counts)
    handle = ctypes.c_void_p()
    error = ctypes.create_string_buffer(256)
    assert table.create_fingerprints(len(targets), words, packed_targets, counts,
                                     ctypes.byref(handle), error, len(error)) == 0, error.value
    assert handle.value
    try:
        stride = max(target_counts) + 1
        boundaries = [-1.0, 0.0, math.nextafter(1 / 3, 0.0), 1 / 3,
                      math.nextafter(1 / 3, 1.0), 1.0, math.nextafter(1.0, 2.0), math.nan]
        for threshold in boundaries:
            cutoffs = (ctypes.c_int32 * (stride * query_count))(
                *(_cutoff(query_bits, target_bits, threshold)
                  for query_bits in query_counts for target_bits in range(stride)))
            hits = (Hit * (query_count * len(targets)))()
            offsets = (ctypes.c_uint64 * (query_count + 1))()
            assert table.fingerprint_screen(handle, query_count, packed_queries, cutoffs, stride,
                                             hits, len(hits), offsets, error, len(error)) == 0, error.value
            assert offsets[0] == 0 and offsets[-1] <= len(hits)
            for query_index, query in enumerate(queries):
                expected = []
                for target_index, target in enumerate(targets):
                    common = sum((left & right).bit_count() for left, right in zip(query, target))
                    if _score(common, query_counts[query_index], target_counts[target_index]) >= threshold:
                        expected.append((target_index, common))
                actual = [(hits[index].target_index, hits[index].common_bits)
                          for index in range(offsets[query_index], offsets[query_index + 1])]
                assert actual == expected
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value

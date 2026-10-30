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


def _create_resident_fingerprints(targets, words):
    table = Table.from_address(smsd_cuda._library.smsd_cuda_screen_accelerator())
    assert table.abi_version == 1 and table.struct_size >= ctypes.sizeof(Table)
    packed_targets = (ctypes.c_uint64 * (words * len(targets)))(
        *(target[word] for word in range(words) for target in targets))
    target_counts = (ctypes.c_int32 * len(targets))(
        *(sum(word.bit_count() for word in target) for target in targets))
    handle = ctypes.c_void_p()
    error = ctypes.create_string_buffer(256)
    assert table.create_fingerprints(len(targets), words, packed_targets, target_counts,
                                     ctypes.byref(handle), error, len(error)) == 0, error.value
    assert handle.value
    return table, handle, error


def _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets):
    assert offsets[0] == 0
    assert all(0 <= offsets[index] <= offsets[index + 1] <= len(hits)
               for index in range(len(queries)))
    expected_total = 0
    for query_index, query in enumerate(queries):
        expected = []
        for target_index, target in enumerate(targets):
            common = sum((left & right).bit_count() for left, right in zip(query, target))
            target_bits = sum(word.bit_count() for word in target)
            if common >= cutoffs[query_index * stride + target_bits]:
                expected.append((target_index, common))
        actual = [(hits[index].target_index, hits[index].common_bits)
                  for index in range(offsets[query_index], offsets[query_index + 1])]
        assert actual == expected
        expected_total += len(expected)
    assert offsets[-1] == expected_total


@pytest.mark.parametrize("target_count,query_count,words", [
    (0, 4, 2), (1, 1, 1), (255, 3, 2), (256, 4, 32), (257, 5, 32),
    (257, 8, 32), (255, 16, 1), (256, 32, 33), (257, 0, 2),
    (17, 7, 32), (17, 15, 32),
    (3, 17, 255), (3, 17, 256), (3, 17, 257), (3, 17, 511),
    (3, 17, 512), (3, 17, 513), (3, 17, 1024), (3, 17, 1025),
])
def test_resident_integer_cutoffs_and_batch_tails(target_count, query_count, words):
    targets = [_fingerprint(index, words) for index in range(target_count)]
    queries = [_fingerprint(index, words) for index in range(query_count)]
    packed_queries = (ctypes.c_uint64 * (words * query_count))(
        *(word for query in queries for word in query))
    maximum_bits = max((sum(word.bit_count() for word in target) for target in targets), default=0)
    stride = maximum_bits + 1
    table, handle, error = _create_resident_fingerprints(targets, words)
    try:
        # Different queries accept or reject the same target popcount in each batch.
        for pass_index in range(5):
            cutoffs = (ctypes.c_int32 * (stride * query_count))(
                *([-1, 0, target_bits, target_bits + 1, (1 << 31) - 1]
                  [(query_index + pass_index) % 5]
                  for query_index in range(query_count) for target_bits in range(stride)))
            hits = (Hit * (query_count * target_count))()
            offsets = (ctypes.c_uint64 * (query_count + 1))(
                *([(1 << 64) - 1] * (query_count + 1)))
            assert table.fingerprint_screen(handle, query_count, packed_queries, cutoffs, stride,
                                             hits, len(hits), offsets, error, len(error)) == 0, error.value
            _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


@pytest.mark.parametrize("invalid_argument", ["short-capacity", "short-stride"])
def test_resident_rejects_short_buffers_and_keeps_handle(invalid_argument):
    targets = [[0], [1], [3]]
    queries = [[0], [1], [3], [7], [15]]
    packed_queries = (ctypes.c_uint64 * len(queries))(*(query[0] for query in queries))
    stride = 3
    cutoffs = (ctypes.c_int32 * (stride * len(queries)))()
    hits = (Hit * (len(queries) * len(targets)))()
    offsets = (ctypes.c_uint64 * (len(queries) + 1))()
    table, handle, error = _create_resident_fingerprints(targets, 1)
    try:
        capacity = len(hits) - 1 if invalid_argument == "short-capacity" else len(hits)
        supplied_stride = stride - 1 if invalid_argument == "short-stride" else stride
        assert table.fingerprint_screen(handle, len(queries), packed_queries, cutoffs, supplied_stride,
                                         hits, capacity, offsets, error, len(error)) != 0
        assert error.value
        assert table.fingerprint_screen(handle, len(queries), packed_queries, cutoffs, stride,
                                         hits, len(hits), offsets, error, len(error)) == 0, error.value
        _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def _top_k_function():
    function = smsd_cuda._library.smsd_cuda_fingerprint_top_k
    function.restype = ctypes.c_int
    function.argtypes = [ctypes.c_void_p, ctypes.c_int32, ctypes.POINTER(ctypes.c_uint64),
                         ctypes.POINTER(ctypes.c_int32), ctypes.c_int32, ctypes.c_int32,
                         ctypes.POINTER(Hit), ctypes.c_uint64, ctypes.POINTER(ctypes.c_uint64),
                         ctypes.c_void_p, ctypes.c_size_t]
    return function


def _top_k_result(handle, targets, queries, words, threshold, k):
    query_counts = [sum(word.bit_count() for word in query) for query in queries]
    target_counts = [sum(word.bit_count() for word in target) for target in targets]
    stride = max(target_counts, default=0) + 1
    packed_queries = (ctypes.c_uint64 * (words * len(queries)))(
        *(word for query in queries for word in query))
    cutoffs = (ctypes.c_int32 * (stride * len(queries)))(
        *(_cutoff(query_bits, target_bits, threshold)
          for query_bits in query_counts for target_bits in range(stride)))
    hits = (Hit * (len(queries) * min(k, len(targets))))()
    offsets = (ctypes.c_uint64 * (len(queries) + 1))()
    error = ctypes.create_string_buffer(256)
    assert _top_k_function()(handle, len(queries), packed_queries, cutoffs, stride, k,
                             hits, len(hits), offsets, error, len(error)) == 0, error.value
    result = []
    for query_index, query in enumerate(queries):
        expected = []
        for target_index, target in enumerate(targets):
            common = sum((left & right).bit_count() for left, right in zip(query, target))
            score = _score(common, query_counts[query_index], target_counts[target_index])
            if score >= threshold:
                expected.append((target_index, common, score))
        expected.sort(key=lambda hit: (-hit[2], hit[0]))
        expected = [(index, common) for index, common, score in expected[:k]]
        actual = [(hits[index].target_index, hits[index].common_bits)
                  for index in range(offsets[query_index], offsets[query_index + 1])]
        assert actual == expected
        result.append(actual)
    assert offsets[0] == 0 and offsets[-1] == sum(map(len, result))
    return result


@pytest.mark.parametrize("words,query_count,target_count", [(1, 1, 0), (1, 0, 9),
    (1, 1, 17), (3, 17, 257), (4, 16, 256), (5, 17, 9), (33, 7, 17)])
def test_resident_top_k_double_thresholds_ties_and_tails(words, query_count, target_count):
    targets = [_fingerprint(index, words) for index in range(target_count)]
    queries = [_fingerprint(index, words) for index in range(query_count)]
    table, handle, error = _create_resident_fingerprints(targets, words)
    try:
        for threshold in [0.0, math.nextafter(1 / 3, 0.0), 1 / 3,
                          math.nextafter(1 / 3, 1.0), 1.0, math.nan]:
            for k in [0, 1, 3, target_count + 1]:
                _top_k_result(handle, targets, queries, words, threshold, k)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def test_resident_top_k_rejects_invalid_arguments_and_keeps_handle():
    targets, queries = [[0], [1], [3]], [[0], [1]]
    table, handle, error = _create_resident_fingerprints(targets, 1)
    try:
        function = _top_k_function()
        packed = (ctypes.c_uint64 * 2)(0, 1)
        cutoffs = (ctypes.c_int32 * 6)()
        hits, offsets = (Hit * 6)(), (ctypes.c_uint64 * 3)()
        for stride, k, capacity in [(2, 2, 6), (3, -1, 6), (3, 2, 3)]:
            for hit in hits:
                hit.target_index, hit.common_bits = -1234567, -7654321
            before = bytes(hits)
            assert function(handle, 2, packed, cutoffs, stride, k, hits, capacity,
                            offsets, error, len(error)) != 0
            assert error.value
            assert bytes(hits) == before
            _top_k_result(handle, targets, queries, 1, 0.0, 2)
        assert function(handle, 0, None, None, 0, 3, None, 0, offsets, error, len(error)) == 0
        assert offsets[0] == 0
        assert function(handle, 2, None, None, 0, 0, None, 0, offsets, error, len(error)) == 0
        assert list(offsets) == [0, 0, 0]
        _top_k_result(handle, targets, queries, 1, 0.0, 2)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def test_resident_top_k_concurrent_reuse():
    from concurrent.futures import ThreadPoolExecutor

    targets = [_fingerprint(index, 3) for index in range(257)]
    queries = [_fingerprint(index, 3) for index in range(17)]
    table, handle, error = _create_resident_fingerprints(targets, 3)
    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(_top_k_result, handle, targets, queries, 3, threshold, 7)
                       for threshold in [0.0, 1 / 3]]
            for future in futures:
                future.result()
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


@pytest.mark.parametrize("words", [1, 3, 4, 5, 33])
def test_opt_in_binary_mma_integer_intersections(monkeypatch, words):
    available = smsd_cuda._library.smsd_cuda_binary_mma_available
    available.restype = ctypes.c_int
    available.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    error = ctypes.create_string_buffer(256)
    status = available(error, len(error))
    assert status >= 0, error.value
    if status == 0:
        pytest.skip("binary MMA prototype needs an enabled build and an SM80-SM90 device")
    monkeypatch.setenv("SMSD_CUDA_BINARY_MMA", "1")
    smsd_cuda.preload()
    targets = [_fingerprint(index, words) for index in range(9)]
    queries = [_fingerprint(index, words) for index in range(17)]
    table, handle, error = _create_resident_fingerprints(targets, words)
    try:
        stride = max(sum(word.bit_count() for word in target) for target in targets) + 1
        packed = (ctypes.c_uint64 * (words * len(queries)))(*(word for query in queries for word in query))
        for threshold in [0.0, math.nextafter(1 / 3, 0.0), 1 / 3, math.nextafter(1 / 3, 1.0), 1.0]:
            cutoffs = (ctypes.c_int32 * (stride * len(queries)))(
                *(_cutoff(sum(word.bit_count() for word in query), target_bits, threshold)
                  for query in queries for target_bits in range(stride)))
            hits, offsets = (Hit * (len(targets) * len(queries)))(), (ctypes.c_uint64 * (len(queries) + 1))()
            assert table.fingerprint_screen(handle, len(queries), packed, cutoffs, stride,
                                             hits, len(hits), offsets, error, len(error)) == 0, error.value
            _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def test_binary_mma_request_keeps_unsupported_device_on_portable_path(monkeypatch):
    available = smsd_cuda._library.smsd_cuda_binary_mma_available
    available.restype = ctypes.c_int
    available.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    error = ctypes.create_string_buffer(256)
    status = available(error, len(error))
    assert status >= 0, error.value
    if status != 0:
        pytest.skip("this control requires the portable fallback device or build")
    monkeypatch.setenv("SMSD_CUDA_BINARY_MMA", "1")
    targets = [[0, 0, 0], [1, 2, 4], [3, 6, 12], [1, 2, 4], [7, 7, 7],
               [8, 16, 32], [0, 0, 0], [15, 15, 15], [1, 0, 0]]
    queries = [targets[index % len(targets)] for index in range(17)]
    table, handle, error = _create_resident_fingerprints(targets, 3)
    try:
        stride = max(sum(word.bit_count() for word in target) for target in targets) + 1
        packed = (ctypes.c_uint64 * (3 * len(queries)))(*(word for query in queries for word in query))
        cutoffs = (ctypes.c_int32 * (stride * len(queries)))()
        hits, offsets = (Hit * (len(targets) * len(queries)))(), (ctypes.c_uint64 * (len(queries) + 1))()
        assert table.fingerprint_screen(handle, len(queries), packed, cutoffs, stride,
                                         hits, len(hits), offsets, error, len(error)) == 0, error.value
        _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


@pytest.mark.parametrize("failure", [1, 2])
def test_resident_top_k_injected_failure_keeps_payload_and_handle(failure):
    inject = getattr(smsd_cuda._library, "smsd_cuda_test_fail_next", None)
    if inject is None:
        pytest.skip("failure injection requires the separate SMSD_CUDA_TESTING library")
    inject.argtypes, inject.restype = [ctypes.c_int], None
    targets, queries = [[0], [1], [3]], [[0], [1]]
    table, handle, error = _create_resident_fingerprints(targets, 1)
    try:
        packed, cutoffs = (ctypes.c_uint64 * 2)(0, 1), (ctypes.c_int32 * 6)()
        hits, offsets = (Hit * 4)(), (ctypes.c_uint64 * 3)()
        for hit in hits:
            hit.target_index, hit.common_bits = -1234567, -7654321
        before = bytes(hits)
        inject(failure)
        assert _top_k_function()(handle, 2, packed, cutoffs, 3, 2, hits, len(hits),
                                 offsets, error, len(error)) != 0
        assert error.value and bytes(hits) == before
        _top_k_result(handle, targets, queries, 1, 0.0, 2)
    finally:
        inject(0)
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def _screen_prepared(handle, table, targets, queries, words, threshold, require_success=True):
    stride = max((sum(word.bit_count() for word in target) for target in targets), default=0) + 1
    packed = (ctypes.c_uint64 * (words * len(queries)))(*(word for query in queries for word in query))
    cutoffs = (ctypes.c_int32 * (stride * len(queries)))(
        *(_cutoff(sum(word.bit_count() for word in query), target_bits, threshold)
          for query in queries for target_bits in range(stride)))
    hits, offsets = (Hit * (len(targets) * len(queries)))(), (ctypes.c_uint64 * (len(queries) + 1))()
    for hit in hits:
        hit.target_index, hit.common_bits = -1234567, -7654321
    before = bytes(hits)
    error = ctypes.create_string_buffer(256)
    status = table.fingerprint_screen(handle, len(queries), packed, cutoffs, stride,
                                      hits, len(hits), offsets, error, len(error))
    if require_success:
        assert status == 0, error.value
        _assert_integer_hits(targets, queries, cutoffs, stride, hits, offsets)
    else:
        assert status != 0 and error.value
        assert bytes(hits) == before


def test_resident_graph_replay_changed_inputs_and_shapes(monkeypatch):
    monkeypatch.setenv("SMSD_CUDA_GRAPHS", "1")
    monkeypatch.setenv("SMSD_CUDA_POOL", "1")
    monkeypatch.delenv("SMSD_CUDA_BINARY_MMA", raising=False)
    targets = [_fingerprint(index, 3) for index in range(257)]
    table, handle, error = _create_resident_fingerprints(targets, 3)
    try:
        for query_count in [4, 16, 8, 17, 1, 17]:
            for shift, threshold in enumerate([0.0, 1 / 3, 1.0]):
                queries = [_fingerprint(index + shift, 3) for index in range(query_count)]
                _screen_prepared(handle, table, targets, queries, 3, threshold)
        # Disabling the lease must preserve the same ordinary callback contract.
        monkeypatch.setenv("SMSD_CUDA_POOL", "0")
        _screen_prepared(handle, table, targets, [[1, 2, 4]], 3, 1 / 3)
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


def test_resident_graph_concurrent_reuse(monkeypatch):
    from concurrent.futures import ThreadPoolExecutor

    monkeypatch.setenv("SMSD_CUDA_GRAPHS", "1")
    monkeypatch.setenv("SMSD_CUDA_POOL", "1")
    monkeypatch.delenv("SMSD_CUDA_BINARY_MMA", raising=False)
    targets = [_fingerprint(index, 3) for index in range(257)]
    table, handle, error = _create_resident_fingerprints(targets, 3)
    try:
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(_screen_prepared, handle, table, targets,
                                       [_fingerprint(index + shift, 3) for index in range(query_count)],
                                       3, threshold)
                       for query_count, shift, threshold in [(8, 0, 0.0), (17, 3, 1 / 3)]]
            for future in futures:
                future.result()
    finally:
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value


@pytest.mark.parametrize("failure", [1, 2])
def test_resident_graph_failure_recovery_preserves_payload(monkeypatch, failure):
    inject = getattr(smsd_cuda._library, "smsd_cuda_test_fail_next", None)
    if inject is None:
        pytest.skip("failure injection requires the separate SMSD_CUDA_TESTING library")
    inject.argtypes, inject.restype = [ctypes.c_int], None
    monkeypatch.setenv("SMSD_CUDA_GRAPHS", "1")
    monkeypatch.setenv("SMSD_CUDA_POOL", "1")
    monkeypatch.delenv("SMSD_CUDA_BINARY_MMA", raising=False)
    targets = [_fingerprint(index, 3) for index in range(17)]
    table, handle, error = _create_resident_fingerprints(targets, 3)
    try:
        _screen_prepared(handle, table, targets, targets[:4], 3, 0.0)
        # A larger batch requires new scratch, rather than relying on a cached allocation hook.
        inject(failure)
        _screen_prepared(handle, table, targets, targets, 3, 0.0, require_success=False)
        _screen_prepared(handle, table, targets, targets, 3, 0.0)
        _screen_prepared(handle, table, targets, targets[:4], 3, 1 / 3)
    finally:
        inject(0)
        assert table.destroy_resident(handle, error, len(error)) == 0, error.value

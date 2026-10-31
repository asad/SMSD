# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Packed corpus searches agree with ordinary binary fingerprint similarity."""

import struct
from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pytest
import smsd

np = pytest.importorskip("numpy")

from smsd.corpus import FingerprintCorpus


@pytest.fixture
def fingerprints():
    return [[], [], [0, 2, 12], [0, 2], [1, 4], [0, 2, 12]]


def test_stream_round_trip_and_owned_chunks(tmp_path, fingerprints):
    path = tmp_path / "library.fp"
    with FingerprintCorpus.create(path, (iter(row) for row in fingerprints), fp_size=13,
                                  chunk_size=2) as corpus:
        assert len(corpus) == len(fingerprints)
        assert [corpus.fingerprint(index) for index in range(len(corpus))] == fingerprints
        chunks = list(corpus.iter_packed_chunks(chunk_size=2))
        assert [(start, packed.shape) for start, packed in chunks] == [(0, (2, 2)), (2, (2, 2)), (4, (2, 2))]
        assert all(packed.dtype == np.uint8 for _, packed in chunks)
    assert chunks[1][1].tolist() == [[5, 16], [5, 0]]
    with FingerprintCorpus(path) as reopened:
        assert reopened.fingerprint(2) == [0, 2, 12]
    with pytest.raises(ValueError, match="closed"):
        reopened.fingerprint(0)


@pytest.mark.parametrize("chunk_size", [1, 2, 100])
@pytest.mark.parametrize("threshold", [0.0, 2 / 3, 1.0])
def test_threshold_pairs_match_scalar(tmp_path, fingerprints, chunk_size, threshold):
    with FingerprintCorpus.create(tmp_path / "library.fp", fingerprints, fp_size=13) as corpus:
        expected = [(first, second, smsd.tanimoto_coefficient(query, target))
                    for first, query in enumerate(fingerprints)
                    for second, target in enumerate(fingerprints) if first < second
                    and smsd.tanimoto_coefficient(query, target) >= threshold]
        assert list(corpus.threshold_pairs(threshold, chunk_size=chunk_size)) == expected
        expected_hits = [(index, smsd.tanimoto_coefficient(fingerprints[2], target))
                         for index, target in enumerate(fingerprints)
                         if smsd.tanimoto_coefficient(fingerprints[2], target) >= threshold]
        assert list(corpus.similarities(fingerprints[2], threshold=threshold,
                                      chunk_size=chunk_size)) == expected_hits


@pytest.mark.parametrize("include_self", [False, True])
@pytest.mark.parametrize("k", [0, 1, 3, 100])
def test_top_k_exact_ranking_and_input_order_ties(tmp_path, fingerprints, include_self, k):
    with FingerprintCorpus.create(tmp_path / "library.fp", fingerprints, fp_size=13) as corpus:
        expected = []
        for query_index, query in enumerate(fingerprints):
            scores = [(index, smsd.tanimoto_coefficient(query, target))
                      for index, target in enumerate(fingerprints)
                      if include_self or index != query_index]
            scores.sort(key=lambda item: (-item[1], item[0]))
            expected.extend((query_index, index, score) for index, score in scores[:k])
        assert list(corpus.top_k_pairs(k, include_self=include_self, chunk_size=2)) == expected
        thresholded = [pair for pair in expected if pair[2] >= 0.5]
        assert list(corpus.top_k_pairs(k, include_self=include_self,
                                     threshold=0.5, chunk_size=1)) == thresholded


def test_from_smiles_uses_bounded_native_batches(tmp_path, monkeypatch):
    from smsd import _fingerprints

    molecules = ["CCO", "CCN", "c1ccccc1", "CCO", "CC(=O)O"]
    options = smsd.FingerprintOptions.ecfp(radius=2, size=257)
    expected = [smsd.fingerprint(molecule, options=options) for molecule in molecules]
    original = _fingerprints.fingerprints
    batch_sizes = []

    def checked_batch(batch, **kwargs):
        batch_sizes.append(len(batch))
        return original(batch, **kwargs)

    monkeypatch.setattr(_fingerprints, "fingerprints", checked_batch)
    with FingerprintCorpus.from_smiles(tmp_path / "library.fp", iter(molecules), options=options,
                                       chunk_size=2, num_threads=1) as corpus:
        assert corpus.fp_size == 257
        assert [corpus.fingerprint(index) for index in range(len(corpus))] == expected
    assert batch_sizes == [0, 2, 2, 1]


def test_empty_duplicate_and_self_conventions(tmp_path):
    with FingerprintCorpus.from_fingerprints(tmp_path / "empty.fp", [], fp_size=9) as corpus:
        assert len(corpus) == 0
        assert list(corpus.threshold_pairs(0)) == []
        assert list(corpus.top_k_pairs(2)) == []
        assert list(corpus.similarities([])) == []
    with FingerprintCorpus.create(tmp_path / "single.fp", [[8, 8, 0]], fp_size=9) as corpus:
        assert corpus.fingerprint(0) == [0, 8]
        assert list(corpus.threshold_pairs(1, include_self=True)) == [(0, 0, 1.0)]


def test_failed_build_preserves_destination_and_cleans_temporary(tmp_path):
    path = tmp_path / "library.fp"
    with FingerprintCorpus.create(path, [[1]], fp_size=8):
        pass
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        FingerprintCorpus.create(path, [[2]], fp_size=8)
    with pytest.raises(ValueError, match="position"):
        FingerprintCorpus.create(path, [[1], [8]], fp_size=8, overwrite=True)
    assert path.read_bytes() == original
    assert sorted(item.name for item in tmp_path.iterdir()) == ["library.fp"]
    with FingerprintCorpus.create(path, [[2]], fp_size=8, overwrite=True) as corpus:
        assert corpus.fingerprint(0) == [2]


@pytest.mark.parametrize("arguments, error", [
    ({"fp_size": 0}, ValueError), ({"fp_size": True}, TypeError),
    ({"fp_size": 2**24 + 1}, ValueError), ({"chunk_size": 0}, ValueError),
    ({"chunk_size": 1.5}, TypeError), ({"overwrite": 1}, TypeError),
])
def test_invalid_build_controls_do_not_create_files(tmp_path, arguments, error):
    with pytest.raises(error):
        FingerprintCorpus.create(tmp_path / "invalid.fp", [], **arguments)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("row, error", [([-1], ValueError), ([8], ValueError),
                                         ([1.5], TypeError), ([True], TypeError),
                                         ([(1, 2)], TypeError)])
def test_invalid_binary_fingerprints(tmp_path, row, error):
    with pytest.raises(error):
        FingerprintCorpus.create(tmp_path / "invalid.fp", [row], fp_size=8)
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("threshold, error", [(float("nan"), ValueError),
                                               (float("inf"), ValueError),
                                               (-0.1, ValueError), (1.1, ValueError),
                                               (True, TypeError), ("0.5", TypeError)])
def test_invalid_search_controls(tmp_path, threshold, error):
    with FingerprintCorpus.create(tmp_path / "library.fp", [[]], fp_size=8) as corpus:
        with pytest.raises(error):
            list(corpus.threshold_pairs(threshold))
        with pytest.raises(error):
            list(corpus.top_k_pairs(1, threshold=threshold))
        with pytest.raises(error):
            list(corpus.similarities([], threshold=threshold))


def test_reject_truncation_trailing_data_and_padding(tmp_path):
    path = tmp_path / "invalid.fp"
    for contents in [b"", b"wrong format" * 2,
                     struct.pack("<8sQQ", b"SMSDFP1\0", 8, 2) + b"\0",
                     struct.pack("<8sQQ", b"SMSDFP1\0", 8, 0) + b"\0"]:
        path.write_bytes(contents)
        with pytest.raises(ValueError):
            FingerprintCorpus(path)
    path.write_bytes(struct.pack("<8sQQ", b"SMSDFP1\0", 9, 1) + b"\0\x80")
    with FingerprintCorpus(path) as corpus:
        with pytest.raises(ValueError, match="padding"):
            list(corpus.similarities([]))


def test_threshold_pairs_rejects_padding_in_first_query_row(tmp_path):
    path = tmp_path / "invalid.fp"
    path.write_bytes(struct.pack("<8sQQ", b"SMSDFP1\0", 9, 2) + b"\x01\x80\x01\0")
    with FingerprintCorpus(path) as corpus:
        with pytest.raises(ValueError, match="padding"):
            list(corpus.threshold_pairs(0.0))


@pytest.mark.parametrize("read", ["fingerprint", "chunks", "scores"])
def test_close_waits_for_active_mapped_read(tmp_path, read):
    corpus = FingerprintCorpus.create(tmp_path / "library.fp", [[0, 8], [1]], fp_size=9)
    mapped = corpus._packed
    entered, release, closing, closed = Event(), Event(), Event(), Event()

    class PausedMapping:
        _mmap = mapped._mmap

        def __getitem__(self, index):
            entered.set()
            assert release.wait(2)
            return mapped[index]

    corpus._packed = PausedMapping()

    def mapped_read():
        if read == "fingerprint":
            return corpus.fingerprint(0)
        if read == "chunks":
            return next(corpus.iter_packed_chunks(chunk_size=1))[1].tolist()
        return next(corpus.similarities([0, 8], chunk_size=1))

    def close():
        closing.set()
        corpus.close()
        closed.set()

    try:
        with ThreadPoolExecutor(max_workers=2) as workers:
            result = workers.submit(mapped_read)
            assert entered.wait(2)
            close_result = workers.submit(close)
            assert closing.wait(2)
            try:
                assert not closed.wait(0.05)
                assert not mapped._mmap.closed
            finally:
                release.set()
            assert result.result(timeout=2) == {"fingerprint": [0, 8], "chunks": [[1, 1]], "scores": (0, 1.0)}[read]
            close_result.result(timeout=2)
        assert mapped._mmap.closed
        with pytest.raises(ValueError, match="closed"):
            corpus.fingerprint(0)
    finally:
        release.set()
        corpus.close()


@pytest.mark.parametrize("row_bytes", [1, 3, 7, 8, 9, 16, 33])
def test_bit_counts_match_byte_table(monkeypatch, row_bytes):
    from smsd import _packed_accelerator

    rng = np.random.default_rng(row_bytes)
    rows = rng.integers(0, 256, size=(37, row_bytes), dtype=np.uint8)
    rows[0] = 0
    rows[1] = 255
    rows[2] = 0x55
    rows[3] = 0xAA
    expected = np.unpackbits(rows, axis=1).sum(axis=1)
    samples = [rows, rows[4], rows[:0], rows[:, ::-1], np.bitwise_and(rows, rows[5])]
    for dtype in (np.uint64, np.int32):
        for sample in samples:
            monkeypatch.setattr(_packed_accelerator, "_BITWISE_COUNT", False)
            table = _packed_accelerator._bit_counts(sample, dtype, np)
            monkeypatch.setattr(_packed_accelerator, "_BITWISE_COUNT", True)
            counted = _packed_accelerator._bit_counts(sample, dtype, np)
            assert np.asarray(counted).dtype == np.asarray(table).dtype == dtype
            assert np.array_equal(counted, table)
        assert np.array_equal(_packed_accelerator._bit_counts(rows, dtype, np), expected)


@pytest.mark.parametrize("fp_size", [13, 64, 197, 1024])
def test_corpus_searches_do_not_depend_on_bit_counting(tmp_path, monkeypatch, fp_size):
    from smsd import _packed_accelerator
    from smsd.clustering import butina, maxmin

    rng = np.random.default_rng(fp_size)
    rows = [sorted(set(rng.integers(0, fp_size, size=int(rng.integers(0, fp_size // 3 + 2))).tolist()))
            for _ in range(40)] + [[], list(range(fp_size))]
    with FingerprintCorpus.create(tmp_path / "library.fp", rows, fp_size=fp_size) as corpus:
        def searches():
            return (list(corpus.threshold_pairs(0.2, chunk_size=7)),
                    list(corpus.similarities(rows[3], chunk_size=5)),
                    list(corpus.screen_many(rows[:4] + [[]], threshold=0.1, chunk_size=9)),
                    butina(corpus, threshold=0.3, chunk_size=6), maxmin(corpus, 5, chunk_size=8))

        monkeypatch.setattr(_packed_accelerator, "_BITWISE_COUNT", False)
        table = searches()
        monkeypatch.setattr(_packed_accelerator, "_BITWISE_COUNT", True)
        assert searches() == table
    assert table[1] == [(index, smsd.tanimoto_coefficient(rows[3], target)) for index, target in enumerate(rows)]

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Packed CPU and explicit CUDA results retain exact set semantics."""

import builtins
import math
import os
import sys
import types

import pytest

np = pytest.importorskip("numpy")
from smsd.corpus import FingerprintCorpus
from smsd.clustering import butina, maxmin
from smsd import _packed_accelerator as accelerator


def _score(first, second):
    first, second = set(first), set(second)
    return len(first & second) / len(first | second) if first | second else 1.0


def _rows(size):
    return [[], [0], list(range(min(size, 3))), [size - 1], [0], [],
            list(range(0, size, 3)), list(range(1, size, 3))]


@pytest.mark.parametrize("value", ["false", "no", "FALSE", "1", "yes"])
def test_disable_flag_matches_existing_contract_before_import(monkeypatch, value):
    monkeypatch.setenv("SMSD_DISABLE_GPU", value)
    original = builtins.__import__

    def guarded(name, *args, **kwargs):
        assert name != "smsd_cuda", "disabled path must not import the CUDA provider"
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", guarded)
    with pytest.raises(RuntimeError, match="disabled"):
        accelerator._PackedCuda()


@pytest.mark.parametrize("value", ["", "0"])
def test_enabled_flag_reaches_provider_import(monkeypatch, value):
    monkeypatch.setenv("SMSD_DISABLE_GPU", value)
    original = builtins.__import__

    def missing(name, *args, **kwargs):
        if name == "smsd_cuda":
            raise ImportError("provider reached")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", missing)
    with pytest.raises(ImportError, match="provider reached"):
        accelerator._PackedCuda()


def test_old_device_fails_explicit_cuda_and_leaves_auto_on_the_cpu(tmp_path, monkeypatch):
    from smsd import __version__

    reason = "the CUDA device Tesla P4 has compute capability 6.1; CUDA screening needs 7.5 or later"
    provider = types.ModuleType("smsd_cuda")
    provider.__smsd_version__ = __version__
    provider.is_available = lambda: False
    provider.unavailable_reason = lambda: reason
    monkeypatch.delenv("SMSD_DISABLE_GPU", raising=False)
    monkeypatch.setitem(sys.modules, "smsd_cuda", provider)
    assert accelerator._cuda_backend("auto") is None
    with pytest.raises(RuntimeError, match="compute capability 6.1") as failure:
        accelerator._cuda_backend("cuda")
    assert str(failure.value.__cause__) == reason
    with FingerprintCorpus.create(tmp_path / "bits.fp", [[0], [1]], fp_size=8) as corpus:
        assert list(corpus.similarities([0], backend="auto")) == [(0, 1.0), (1, 0.0)]
        with pytest.raises(RuntimeError, match="compute capability 6.1"):
            list(corpus.similarities([0], backend="cuda"))
        assert corpus.cuda_calls == 0


@pytest.mark.parametrize("threshold", [0.0, math.nextafter(1 / 3, 0), 1 / 3,
                                       math.nextafter(1 / 3, 1), 1.0])
def test_integer_cutoffs_match_binary64_acceptance(threshold):
    counts = np.array([0, 1, 3, 8], dtype=np.int32)
    cutoffs = accelerator._cutoffs(counts, 9, threshold, np)
    for row, query_count in enumerate(counts):
        for target_count in range(10):
            for common in range(min(query_count, target_count) + 1):
                union = int(query_count) + target_count - common
                score = common / union if union else 1.0
                assert (common >= cutoffs[row, target_count]) == (score >= threshold)


def test_bounded_batch_order_empty_queries_and_excess(tmp_path):
    rows = _rows(13)
    with FingerprintCorpus.create(tmp_path / "bits.fp", rows, fp_size=13) as corpus:
        expected = [(query, target, _score(rows[query], rows[target]))
                    for start in range(0, len(rows), 3)
                    for query in range(len(rows))
                    for target in range(start, min(start + 3, len(rows)))
                    if _score(rows[query], rows[target]) >= 1 / 3]
        assert list(corpus.screen_many(iter(rows), threshold=1 / 3, chunk_size=3)) == expected
        assert list(corpus.screen_many([])) == []
        consumed = []

        def excessive():
            for index in range(100):
                consumed.append(index)
                yield []

        with pytest.raises(ValueError, match="32 queries"):
            list(corpus.screen_many(excessive()))
        assert len(consumed) == 33
        assert corpus.cuda_calls == 0


def test_strict_backend_and_mid_call_error_never_fall_back(tmp_path, monkeypatch):
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    with FingerprintCorpus.create(tmp_path / "bits.fp", [[0], [1]], fp_size=8) as corpus:
        assert list(corpus.similarities([0], backend="auto")) == [(0, 1.0), (1, 0.0)]
        with pytest.raises(RuntimeError, match="unavailable"):
            list(corpus.similarities([0], backend="cuda"))

        class Broken:
            def select(self, *args):
                raise RuntimeError("device operation failed")

        monkeypatch.setattr(accelerator, "_cuda_backend", lambda backend: Broken())
        with pytest.raises(RuntimeError, match="device operation failed"):
            list(corpus.similarities([0], backend="auto"))
        assert corpus.cuda_calls == 0


def test_invalid_controls_are_rejected_before_provider_lookup(tmp_path, monkeypatch):
    def forbidden(*args):
        raise AssertionError("invalid control touched the optional runtime")

    monkeypatch.setattr(accelerator, "_cuda_backend", forbidden)
    with FingerprintCorpus.create(tmp_path / "empty.fp", [], fp_size=8) as corpus:
        for operation in [lambda: list(corpus.screen_many([], backend="bogus")),
                          lambda: list(corpus.threshold_pairs(0, backend="bogus")),
                          lambda: list(corpus.top_k_pairs(0, backend="bogus")),
                          lambda: butina(corpus, backend="bogus"),
                          lambda: maxmin(corpus, 0, backend="bogus")]:
            with pytest.raises(ValueError, match="backend"):
                operation()
        with pytest.raises(TypeError, match="k must be an integer"):
            list(corpus.top_k_pairs(True, backend="cuda"))
        assert list(corpus.top_k_pairs(0, backend="cuda")) == []


@pytest.mark.skipif(os.environ.get("SMSD_TEST_PACKED_CUDA") != "1",
                    reason="requires an explicitly owned actual CUDA test session")
@pytest.mark.parametrize("size,query_count", [(1, 1), (13, 8), (65, 17), (257, 32), (2048, 8)])
def test_actual_cuda_batch_pair_ranking_and_clustering(tmp_path, size, query_count):
    rows = _rows(size)
    queries = [rows[index % len(rows)] for index in range(query_count)]
    with FingerprintCorpus.create(tmp_path / "bits.fp", rows, fp_size=size) as corpus:
        for threshold in [0.0, math.nextafter(1 / 3, 0), 1 / 3,
                          math.nextafter(1 / 3, 1), 1.0]:
            expected = list(corpus.screen_many(queries, threshold=threshold, chunk_size=3))
            actual = list(corpus.screen_many(queries, threshold=threshold, chunk_size=3, backend="cuda"))
            assert actual == expected
        for include_self in [False, True]:
            for k in [1, 3, 100]:
                assert list(corpus.top_k_pairs(k, threshold=1 / 3, include_self=include_self,
                                              chunk_size=3, backend="cuda")) == list(
                    corpus.top_k_pairs(k, threshold=1 / 3, include_self=include_self, chunk_size=3))
            assert list(corpus.threshold_pairs(1 / 3, include_self=include_self,
                                               chunk_size=3, backend="cuda")) == list(
                corpus.threshold_pairs(1 / 3, include_self=include_self, chunk_size=3))
        for reordering in [False, True]:
            assert butina(corpus, 1 / 3, reordering=reordering, chunk_size=3, backend="cuda") == butina(
                corpus, 1 / 3, reordering=reordering, chunk_size=3)
        assert maxmin(corpus, 5, seed=2, chunk_size=3, backend="cuda") == maxmin(
            corpus, 5, seed=2, chunk_size=3)
        assert corpus.cuda_calls > 0

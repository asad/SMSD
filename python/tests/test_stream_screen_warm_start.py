# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Streaming screens preserve public hit indices, preparation and input failures."""

import pytest
import smsd

from smsd import _accelerator
from smsd.corpus import stream_screen_many_smiles, stream_screen_smiles


@pytest.mark.parametrize("chunk_size", [1, 2, 100])
@pytest.mark.parametrize("threshold", [0.0, 0.3, 1.0])
def test_streamed_screen_matches_scalar_fingerprints(chunk_size, threshold):
    query = "CCO"
    targets = ["CCO", "CCN", "c1ccccc1", "CCO", "CC(=O)O"]
    query_fp = smsd.fingerprint(query, kind="ecfp")
    expected = [index for index, target in enumerate(targets)
                if smsd.tanimoto_coefficient(query_fp, smsd.fingerprint(target, kind="ecfp")) >= threshold]
    assert list(stream_screen_smiles(query, iter(targets), threshold, chunk_size=chunk_size,
                                    num_threads=1)) == expected


def test_streamed_hits_arrive_before_later_input_is_consumed():
    consumed = []

    def targets():
        for index, smiles in enumerate(["CCO", "CCN", "CCO", "invalid smiles!"]):
            consumed.append(index)
            yield smiles

    hits = stream_screen_smiles("CCO", targets(), 1.0, chunk_size=2, num_threads=1)
    assert next(hits) == 0
    assert consumed == [0, 1]
    with pytest.raises(ValueError):
        list(hits)
    assert consumed == [0, 1, 2, 3]


def test_streamed_empty_and_invalid_controls():
    assert list(stream_screen_smiles("CCO", [], 0.5)) == []
    with pytest.raises(ValueError):
        list(stream_screen_smiles("CCO", [], 0.5, chunk_size=0))
    with pytest.raises(ValueError):
        list(stream_screen_smiles("CCO", [], float("nan")))


def test_warm_start_uses_and_keeps_prepared_real_corpus(monkeypatch):
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    corpus = smsd.TargetCorpus.from_smiles(["CCO", "CCN", "CCO"])
    assert not corpus.is_prewarmed
    assert _accelerator.warm_start(corpus, "CCO", 1.0, num_threads=1) == [0, 2]
    assert corpus.is_prewarmed
    assert _accelerator.warm_start(corpus, "CCN", 1.0, num_threads=1) == [1]
    assert corpus.screen_many(["CCO"], 1.0, num_threads=1) == [[0, 2]]


def test_warm_start_validates_before_preparing():
    corpus = smsd.TargetCorpus.from_smiles(["CCO"])
    with pytest.raises(TypeError):
        _accelerator.warm_start([], "CCO")
    with pytest.raises(ValueError):
        _accelerator.warm_start(corpus, "CCO", float("inf"))
    with pytest.raises(ValueError):
        _accelerator.warm_start(corpus, "CCO", num_threads=-1)
    assert not corpus.is_prewarmed


def test_multiquery_tiles_match_scalar_and_bound_each_native_call(monkeypatch):
    from smsd import _screening

    queries = ["CCO", "CCN", "c1ccccc1"]
    targets = ["CCO", "CCN", "CCO", "c1ccccc1", "CCN"]
    original = _screening.TargetCorpus.screen_many
    tiles = []

    def screen_many(corpus, batch, threshold, **kwargs):
        tiles.append((len(batch), len(corpus)))
        return original(corpus, batch, threshold, **kwargs)

    monkeypatch.setattr(_screening.TargetCorpus, "screen_many", screen_many)
    expected = []
    for start in range(0, len(targets), 2):
        for query_index, query in enumerate(queries):
            for target_index in range(start, min(start + 2, len(targets))):
                score = smsd.tanimoto_coefficient(smsd.fingerprint(query, kind="ecfp"),
                                                  smsd.fingerprint(targets[target_index], kind="ecfp"))
                if score >= 0.3:
                    expected.append((query_index, target_index))
    assert list(stream_screen_many_smiles(iter(queries), iter(targets), 0.3, query_chunk_size=2,
                                         target_chunk_size=2, max_queries=3, num_threads=1)) == expected
    assert tiles == [(2, 2), (1, 2), (2, 2), (1, 2), (2, 1), (1, 1)]


def test_multiquery_hits_precede_later_target_consumption():
    consumed = []

    def targets():
        for index, target in enumerate(["CCO", "CCN", "CCO", "invalid smiles!"]):
            consumed.append(index)
            yield target

    hits = stream_screen_many_smiles(["CCO", "CCN"], targets(), 1, query_chunk_size=1,
                                     target_chunk_size=2, num_threads=1)
    assert next(hits) == (0, 0)
    assert consumed == [0, 1]
    assert next(hits) == (1, 1)
    with pytest.raises(ValueError):
        list(hits)


@pytest.mark.parametrize("controls", [{"max_queries": 0}, {"query_chunk_size": 0},
                                       {"target_chunk_size": 2**63},
                                       {"query_chunk_size": 2048, "target_chunk_size": 2048}])
def test_multiquery_rejects_invalid_bounds_before_consuming_inputs(controls):
    def untouched():
        raise AssertionError("input was consumed")
        yield

    with pytest.raises(ValueError):
        list(stream_screen_many_smiles(untouched(), untouched(), 0.5, **controls))


def test_multiquery_count_limit_and_empty_query_do_not_consume_targets():
    def untouched():
        raise AssertionError("targets were consumed")
        yield

    with pytest.raises(ValueError, match="max_queries"):
        list(stream_screen_many_smiles(iter(["CCO", "CCN"]), untouched(), 0.5, max_queries=1))
    assert list(stream_screen_many_smiles([], untouched(), 0.5)) == []

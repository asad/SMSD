# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Batched corpus screening preserves the ordered results of individual calls."""
import math
import sys

import pytest
import smsd
from smsd import _accelerator


@pytest.fixture(autouse=True)
def cpu_only(monkeypatch):
    monkeypatch.setenv("SMSD_DISABLE_GPU", "1")
    monkeypatch.setitem(sys.modules, "smsd_cuda", None)
    _accelerator._reset()
    yield
    _accelerator._reset()


@pytest.fixture
def graphs():
    return [smsd.parse_smiles(text) for text in ("C", "CCO", "CC(=O)O", "c1ccccc1", "c1ccncc1")]


@pytest.mark.parametrize("threshold", [0.0, 0.5, 1.0, -math.inf, math.inf, math.nan])
def test_batches_equal_ordered_single_screens(graphs, threshold):
    corpus = smsd.TargetCorpus()
    corpus.add_targets(graphs * 3)
    assert corpus.screen_many(graphs, threshold) == [[] for _ in graphs]
    corpus.prewarm()
    assert corpus.screen_many(graphs, threshold) == [corpus.screen(query, threshold) for query in graphs]
    assert corpus.screen_targets_many(graphs, threshold) == [
        smsd.screen_targets(query, corpus, threshold) for query in graphs
    ]


def test_batches_keep_exact_score_boundaries(graphs):
    corpus = smsd.TargetCorpus()
    corpus.add_targets(graphs)
    corpus.prewarm()
    fingerprints = [smsd.circular_fingerprint(graph, 2, 2048) for graph in graphs]
    score = smsd.tanimoto_coefficient(fingerprints[0], fingerprints[1])
    for threshold in (math.nextafter(score, -math.inf), score, math.nextafter(score, math.inf)):
        assert corpus.screen_many(graphs, threshold) == [corpus.screen(query, threshold) for query in graphs]


def test_empty_batches_and_corpora(graphs):
    corpus = smsd.TargetCorpus()
    assert corpus.screen_many([], 0.0) == []
    assert corpus.screen_targets_many([], 0.0) == []
    assert corpus.screen_many(graphs, 0.0) == [[] for _ in graphs]
    assert corpus.screen_targets_many(graphs, 0.0) == [[] for _ in graphs]


def test_batch_label_chemistry_matches_singles(graphs):
    corpus = smsd.TargetCorpus()
    corpus.add_targets(graphs)
    chemistry = smsd.ChemOptions()
    chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
    chemistry.ring_matches_ring_only = True
    assert corpus.screen_targets_many(graphs, 0.5, chemistry) == [
        smsd.screen_targets(query, corpus, 0.5, chemistry) for query in graphs
    ]

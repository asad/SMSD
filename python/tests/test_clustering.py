# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Cluster membership and diversity selection follow scalar Tanimoto distances."""

import pytest
import smsd

pytest.importorskip("numpy")

from smsd.clustering import butina, maxmin
from smsd.corpus import FingerprintCorpus


def scalar_butina(fingerprints, threshold, reordering):
    neighbours = [{second for second, target in enumerate(fingerprints)
                   if first != second and smsd.tanimoto_coefficient(query, target) >= threshold}
                  for first, query in enumerate(fingerprints)]
    active = set(range(len(fingerprints)))
    clusters = []
    while active:
        centre = min(active, key=lambda index: (-len(neighbours[index] & active if reordering
                                                   else neighbours[index]), index))
        members = (centre, *sorted(neighbours[centre] & active))
        clusters.append(members)
        active.difference_update(members)
    return clusters


@pytest.mark.parametrize("reordering", [False, True])
@pytest.mark.parametrize("threshold", [0.0, 0.4, 0.6, 1.0])
def test_butina_matches_scalar_and_centre_contract(tmp_path, reordering, threshold):
    fingerprints = [[], [], [0, 1], [0, 1, 2], [1, 2], [2, 3], [3], [3, 4], [6]]
    with FingerprintCorpus.create(tmp_path / "library.fp", fingerprints, fp_size=8) as corpus:
        actual = butina(corpus, threshold, reordering=reordering, chunk_size=2)
        assert actual == scalar_butina(fingerprints, threshold, reordering)
        assert sorted(index for cluster in actual for index in cluster) == list(range(len(fingerprints)))
        for centre, *members in actual:
            assert all(smsd.tanimoto_coefficient(fingerprints[centre], fingerprints[index]) >= threshold
                       for index in members)


@pytest.mark.parametrize("seed", [0, 2, 4])
@pytest.mark.parametrize("count", [0, 1, 3, 6])
def test_maxmin_matches_scalar_selection(tmp_path, seed, count):
    fingerprints = [[], [0, 1], [0], [2, 3], [2], []]
    expected = []
    while len(expected) < count:
        if not expected:
            expected.append(seed)
        else:
            candidates = set(range(len(fingerprints))) - set(expected)
            expected.append(min(candidates, key=lambda index: (
                -min(1 - smsd.tanimoto_coefficient(fingerprints[index], fingerprints[selected])
                     for selected in expected), index)))
    with FingerprintCorpus.create(tmp_path / "library.fp", fingerprints, fp_size=8) as corpus:
        assert maxmin(corpus, count, seed=seed, chunk_size=1) == expected


def test_empty_and_identical_inputs_have_deterministic_ties(tmp_path):
    with FingerprintCorpus.create(tmp_path / "empty.fp", [], fp_size=8) as corpus:
        assert butina(corpus) == []
        assert maxmin(corpus, 0) == []
        with pytest.raises(ValueError):
            maxmin(corpus, 1)
    with FingerprintCorpus.create(tmp_path / "identical.fp", [[0]] * 4, fp_size=8) as corpus:
        assert butina(corpus, 1) == [(0, 1, 2, 3)]
        assert maxmin(corpus, 4, seed=2) == [2, 0, 1, 3]
        with pytest.raises(TypeError):
            maxmin(corpus, True)
        with pytest.raises(IndexError):
            maxmin(corpus, 1, seed=4)
        with pytest.raises(TypeError):
            butina(corpus, reordering=1)

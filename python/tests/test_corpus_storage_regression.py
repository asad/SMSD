# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Captured target slots keep their values, numbering and order."""
from concurrent.futures import ThreadPoolExecutor
import threading

import smsd
from smsd import _smsd


def fields(result):
    return (result.matches, result.status, result.optimal, result.complete,
            result.score, result.nodes_visited)


def test_repeated_targets_are_captured_before_source_changes():
    carbon = smsd.parse_smiles("CC")
    carbon.properties["tag"] = "before"
    oxygen = smsd.parse_smiles("O")
    corpus = smsd.TargetCorpus.from_smiles([carbon, oxygen, carbon])
    carbon.properties["tag"] = "after"
    snapshot = corpus._snapshot()
    assert len(snapshot) == 3
    assert [snapshot[index][0].properties["tag"] for index in (0, 2)] == ["before", "before"]
    query = smsd.Substructure.find("CC")
    expected = [fields(query.match(text)) for text in ("CC", "O", "CC")]
    assert [fields(result) for result in query.match_many(corpus, num_threads=2)] == expected
    corpus.add_target(carbon)
    assert len(snapshot) == 3 and len(corpus) == 4
    assert corpus._snapshot()[3][0].properties["tag"] == "after"
    corpus.prewarm(num_threads=1)
    assert [fields(result) for result in query.match_many(corpus, num_threads=2)] == expected + expected[:1]


def test_native_value_snapshots_remain_independent():
    graph = smsd.parse_smiles("CCO")
    graph.properties["tag"] = "original"
    corpus = _smsd.TargetCorpus()
    corpus.add_targets([graph, graph])
    first, second = corpus.targets()
    first.properties["tag"] = "changed"
    assert second.properties["tag"] == "original"
    assert [target.properties["tag"] for target in corpus.targets()] == ["original", "original"]


def test_native_iterator_captures_each_yielded_value():
    graph = smsd.parse_smiles("CCO")

    def targets():
        graph.properties["tag"] = "first"
        yield graph
        graph.properties["tag"] = "second"
        yield graph

    corpus = _smsd.TargetCorpus()
    corpus.add_targets(targets())
    assert [target.properties["tag"] for target in corpus.targets()] == ["first", "second"]


def test_native_search_snapshots_are_detached_and_keep_repeated_slots():
    carbon = smsd.parse_smiles("CC")
    carbon.properties["tag"] = "carbon"
    oxygen = smsd.parse_smiles("O")
    oxygen.properties["tag"] = "oxygen"
    corpus = _smsd.TargetCorpus()
    corpus.add_targets((carbon, oxygen, carbon))

    first = corpus._search_targets()
    second = corpus._search_targets()
    assert [graph.properties["tag"] for graph in first] == ["carbon", "oxygen", "carbon"]
    first[0].properties["tag"] = "changed snapshot"
    carbon.properties["tag"] = "changed source"
    assert [graph.properties["tag"] for graph in second] == ["carbon", "oxygen", "carbon"]
    assert [graph.properties["tag"] for graph in corpus.targets()] == ["carbon", "oxygen", "carbon"]
    query = smsd.Substructure.find("CC")
    assert [result.best.atom_count if result.found else 0 for result in query.match_many(second)] == [2, 0, 2]


def test_native_search_snapshots_survive_concurrent_additions():
    carbon = smsd.parse_smiles("CC")
    carbon.properties["tag"] = "carbon"
    oxygen = smsd.parse_smiles("O")
    oxygen.properties["tag"] = "oxygen"
    nitrogen = smsd.parse_smiles("N")
    nitrogen.properties["tag"] = "nitrogen"
    corpus = _smsd.TargetCorpus()
    corpus.add_targets([carbon, oxygen, carbon])
    start = threading.Barrier(2)

    def snapshots():
        start.wait()
        return [corpus._search_targets() for _ in range(4)]

    def append():
        start.wait()
        corpus.add_target(nitrogen)

    with ThreadPoolExecutor(max_workers=2) as workers:
        captured = workers.submit(snapshots)
        added = workers.submit(append)
        results = captured.result()
        added.result()
    for snapshot in results:
        assert len(snapshot) in (3, 4)
        assert [graph.properties["tag"] for graph in snapshot[:3]] == ["carbon", "oxygen", "carbon"]
        if len(snapshot) == 4:
            assert snapshot[3].properties["tag"] == "nitrogen"
    assert len(corpus) == 4

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Screening and bounded MCS retain valid common subgraphs."""

import math

import pytest

import smsd
from .search_helpers import configured_mcs, atom_map


@pytest.mark.parametrize("query,target", [
    ("c1ccccc1CC(=O)O", "C1CCCCC1CC(=O)O"),
    ("C1CCCCC1", "CCCCCC"),
])
def test_default_screening_relaxes_aromaticity_and_ring_membership(query, target):
    query, target = smsd.parse_smiles(query), smsd.parse_smiles(target)
    assert smsd.similarity_upper_bound(query, target) == 1.0
    assert smsd.screen_targets(query, [target], 0.5) == [0]


@pytest.mark.parametrize("operation", ["mapping", "size"])
def test_screened_mcs_keeps_the_full_relaxed_aromatic_match(operation):
    query = smsd.parse_smiles("c1ccccc1CC(=O)O")
    targets = [smsd.parse_smiles("C1CCCCC1CC(=O)O")]
    indices = smsd.screen_targets(query, targets, 0.5)
    assert indices == [0]
    result = smsd.MCS.find(query).match_many([targets[index] for index in indices])[0]
    assert (len(result.best.atom_mapping) if operation == "mapping" else result.best.atom_count) == 10


@pytest.mark.parametrize("as_graph", [False, True])
def test_auto_refines_the_partial_common_alkaloid_graph(as_graph):
    query = "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(=O)CC4"
    target = "COc1ccc2[nH]cc(CCN(C)C)c2c1"
    if as_graph:
        query, target = smsd.parse_smiles(query), smsd.parse_smiles(target)
    mapping = atom_map(configured_mcs(query, timeout_ms=1000).match(target))
    assert len(mapping) >= 14
    query = smsd.parse_smiles(query) if isinstance(query, str) else query
    target = smsd.parse_smiles(target) if isinstance(target, str) else target
    assert smsd.MCS.validate_mapping(query, target, mapping) == []


def test_positive_budget_can_prove_the_optimum():
    result = smsd.MCS.find("CCO").with_timeout(1000).match("CCCO")
    assert result.best.atom_count == 3
    assert result.optimal and result.complete
    assert result.status == smsd.SearchStatus.COMPLETED


def test_node_limit_does_not_claim_an_unproved_optimum():
    options = smsd.MCS.Options()
    options.node_limit = 1
    result = smsd.MCS.find("CCO").with_options(options).match("CCCO")
    assert result.status == smsd.SearchStatus.LIMIT_REACHED
    assert not result.optimal and not result.complete
    for match in result.matches:
        assert smsd.MCS.validate_mapping(smsd.parse_smiles("CCO"), smsd.parse_smiles("CCCO"), match) == []


def test_expired_deadline_reports_incompletion():
    result = smsd.MCS.find("CCO").with_timeout(0).match("CCCO")
    assert result.status == smsd.SearchStatus.TIMED_OUT
    assert not result.optimal and not result.complete


def test_search_never_returns_an_incompatible_candidate():
    mapping = atom_map(configured_mcs('CCO').match('CNN'))
    assert len(mapping) == 1
    assert smsd.MCS.validate_mapping(smsd.parse_smiles("CCO"), smsd.parse_smiles("CNN"), mapping) == []


@pytest.mark.parametrize("mode", ["aromatic", "ring", "untyped"])
def test_screening_uses_the_requested_atom_policy(mode):
    chemistry = smsd.ChemOptions()
    query, target = "C1CCCCC1", "CCCCCC"
    if mode == "aromatic":
        query, target = "c1ccccc1", "C1CCCCC1"
        chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
    elif mode == "ring":
        chemistry.ring_matches_ring_only = True
    else:
        query, target = "C", "N"
        chemistry.match_atom_type = False
    query, target = smsd.parse_smiles(query), smsd.parse_smiles(target)
    expected = 1.0 if mode == "untyped" else 0.0
    assert smsd.similarity_upper_bound(query, target, chemistry) == expected
    assert smsd.screen_targets(query, [target], 0.5, chemistry) == ([0] if expected else [])


SCREEN_LIBRARY = ["c1ccccc1", "CCO", "CC(=O)Oc1ccccc1C(=O)O", "C1CCCCC1", "c1ccncc1", "CCN(CC)CC",
                  "O=C(O)c1ccccc1O", "C[N+](C)(C)C", "ClC(Cl)Cl", "c1ccc2ccccc2c1", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"]


@pytest.mark.parametrize("ring_only", [False, True])
def test_screen_targets_gives_the_serial_indices_on_every_input_and_thread_count(ring_only):
    chemistry = smsd.ChemOptions()
    chemistry.ring_matches_ring_only = ring_only
    smiles = [SCREEN_LIBRARY[(index * 7) % len(SCREEN_LIBRARY)] for index in range(3000)]
    graphs = [smsd.parse_smiles(text) for text in smiles]
    corpus = smsd.TargetCorpus.from_smiles(smiles)
    query = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")
    scores = [smsd.similarity_upper_bound(query, graph, chemistry) for graph in graphs]
    # Thresholds on every score and on the next double above it show that each screened score
    # equals the pairwise one exactly.
    distinct = sorted(set(scores))
    for threshold in distinct + [math.nextafter(score, math.inf) for score in distinct] + [0.0, 1.0]:
        expected = [index for index, score in enumerate(scores) if score >= threshold]
        for threads in (1, 2, 0):
            for targets in (smiles, graphs, corpus, iter(smiles)):
                assert smsd.screen_targets(query, targets, threshold, chemistry, num_threads=threads) == expected


def _label_rules():
    """One chemistry for each label rule of the bound: untyped, typed, strict aromaticity, ring only, both."""
    rules = []
    for typed, strict, ring_only in ((False, False, False), (True, False, False), (True, True, False),
                                     (True, False, True), (True, True, True)):
        chemistry = smsd.ChemOptions()
        chemistry.match_atom_type = typed
        chemistry.ring_matches_ring_only = ring_only
        if strict:
            chemistry.aromaticity_mode = smsd.AromaticityMode.STRICT
        rules.append(chemistry)
    return rules


def test_corpus_screen_keeps_the_exact_scores_across_label_rules_and_additions():
    graphs = [smsd.parse_smiles(SCREEN_LIBRARY[(index * 7) % len(SCREEN_LIBRARY)]) for index in range(3000)]
    corpus = smsd.TargetCorpus()
    corpus.add_targets(graphs)
    query = smsd.parse_smiles("CC(=O)Oc1ccccc1C(=O)O")

    def check(chemistry):
        scores = [smsd.similarity_upper_bound(query, graph, chemistry) for graph in graphs]
        distinct = sorted(set(scores))
        for threshold in distinct + [math.nextafter(score, math.inf) for score in distinct]:
            expected = [index for index, score in enumerate(scores) if score >= threshold]
            for threads in (1, 2, 0):
                assert smsd.screen_targets(query, corpus, threshold, chemistry, num_threads=threads) == expected

    # Each screen follows one under another label rule, in both orders.
    rules = _label_rules()
    for chemistry in rules + rules[::-1]:
        check(chemistry)
    corpus.add_target("c1ccccc1CCN")
    graphs.append(smsd.parse_smiles("c1ccccc1CCN"))
    for chemistry in rules:
        check(chemistry)
    added = [smsd.parse_smiles(text) for text in ("CC(=O)Oc1ccccc1C(=O)O", "C1CCCCC1CCCl")]
    corpus.add_targets(added)
    graphs.extend(added)
    for chemistry in rules:
        check(chemistry)


def test_screen_targets_reads_mixed_and_repeated_inputs():
    graphs = [smsd.parse_smiles(text) for text in SCREEN_LIBRARY]
    mixed = [graph if index % 2 else text for index, (graph, text) in enumerate(zip(graphs, SCREEN_LIBRARY))]
    expected = smsd.screen_targets("c1ccccc1O", graphs, 0.2)
    assert expected and smsd.screen_targets("c1ccccc1O", mixed, 0.2) == expected
    assert smsd.screen_targets("c1ccccc1O", tuple(SCREEN_LIBRARY * 100), 0.2) == [
        index + len(SCREEN_LIBRARY) * repeat for repeat in range(100) for index in expected]
    assert smsd.screen_targets("c1ccccc1O", graphs[:1] * 1000, 0.0, num_threads=0) == list(range(1000))


@pytest.mark.parametrize("targets,error", [
    (["CCO", "C1CC", 5], ValueError),
    (["CCO", 5, "C1CC"], TypeError),
    (["CCO", "", "C1CC"], ValueError),
    (["CCO"] * 600 + ["C1CC"] + [None] * 600, ValueError),
    (["CCO"] * 600 + [None] + ["C1CC"] * 600, TypeError),
])
def test_screen_targets_raises_the_error_of_the_first_failing_target(targets, error):
    for threads in (1, 2, 0):
        with pytest.raises(error):
            smsd.screen_targets("CCO", targets, 0.5, num_threads=threads)


def test_screen_targets_reports_target_errors_before_argument_errors():
    with pytest.raises(ValueError):
        smsd.screen_targets("CCO", ["CCO", "C1CC"], "high")
    with pytest.raises(TypeError, match="threshold"):
        smsd.screen_targets("CCO", ["CCO"], "high")
    with pytest.raises(TypeError, match="chem"):
        smsd.screen_targets("CCO", ["CCO"], 0.5, 5)
    with pytest.raises(TypeError):
        smsd.screen_targets("CCO", smsd.TargetCorpus.from_smiles(["CCO"]), "high")
    for threads in (-1, 1.5, True):
        with pytest.raises(ValueError):
            smsd.screen_targets("CCO", ["CCO"], 0.5, num_threads=threads)


def test_screen_targets_reports_a_converter_that_returns_no_graph():
    query = smsd.parse_smiles("CCO")
    with pytest.raises(TypeError, match="convert must return a MolGraph"):
        smsd._smsd.screen_targets(query, ["CCO", 3], 0.1, smsd.ChemOptions(), 0, lambda item: 5)

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Screening and bounded MCS retain valid common subgraphs."""

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

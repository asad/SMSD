# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Selected-edge search preserves chemistry, raw embeddings and fixed mappings."""
import smsd


def _key(match):
    return tuple(sorted(match.atom_mapping.items())), tuple(sorted(match.bond_mapping))


def test_wildcard_fragments_keep_all_full_embeddings():
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    options.connectivity = smsd.MCS.Connectivity.FRAGMENTS
    options.chemistry.match_atom_type = False
    options.min_fragment_atoms = 2
    options.max_fragments = 2
    result = smsd.MCS.find("CCO.CN").with_options(options).match_all("CCN.CO", limit=100)
    assert result.optimal and result.complete
    assert result.score == 3 and len(result.matches) == 4
    assert len({_key(match) for match in result.matches}) == 4
    for match in result.matches:
        assert (match.atom_count, match.bond_count) == (5, 3)
        assert sorted(map(len, match.fragments)) == [2, 3]
        assert smsd.MCS.validate_mapping("CCO.CN", "CCN.CO", match, options) == []


def test_typed_fragments_do_not_retain_an_unmatched_element():
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    options.connectivity = smsd.MCS.Connectivity.FRAGMENTS
    result = smsd.MCS.find("CO.CN").with_options(options).match_all("NCCN", limit=100)
    assert result.optimal and result.complete
    assert result.score == 1 and len(result.matches) == 2
    assert {tuple(sorted(match.atom_mapping.items())) for match in result.matches} == {
        ((2, 1), (3, 0)), ((2, 2), (3, 3))}


def test_one_fragment_enumerates_the_same_selected_graphs_as_connected_search():
    query = smsd.MCS.find("CC.CC").with_objective(smsd.MCS.Objective.MCES)
    connected = query.match_all("CCCC", limit=100)
    one_fragment = query.with_fragments(max_fragments=1).match_all("CCCC", limit=100)
    assert connected.optimal and connected.complete and len(connected.matches) == 12
    assert one_fragment.optimal and one_fragment.complete
    assert {_key(match) for match in connected.matches} == {_key(match) for match in one_fragment.matches}


def test_maximality_preserves_the_fixed_edge_orientation():
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    extendible = smsd.Match({0: 0, 1: 1}, ((0, 1, 0, 1),), ((0, 1),))
    blocked = smsd.Match({0: 1, 1: 0}, ((0, 1, 1, 0),), ((0, 1),))
    assert not smsd.MCS.is_mapping_maximal("CCO", "CCO", extendible, options)
    assert smsd.MCS.is_mapping_maximal("CCO", "CCO", blocked, options)
    unequal = smsd.Match({1: 0, 2: 1}, ((1, 2, 0, 1),), ((1, 2),))
    assert smsd.MCS.is_mapping_maximal("CCCO", "CCC", unequal, options)

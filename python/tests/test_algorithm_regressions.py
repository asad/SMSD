# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Regression checks for observable MCS options and prepared-query reuse."""

import pytest
import smsd


def test_mcis_and_mces_preserve_requested_topology():
    edge = smsd.MCS.find("CCC").with_objective(smsd.MCS.Objective.MCES).match("C1CC1")
    induced = smsd.MCS.find("CCC").match("C1CC1")
    assert edge.best.atom_count == 3 and edge.best.bond_count == 2
    assert induced.best.atom_count == 2
    assert edge.optimal and induced.optimal


def test_isotope_matching_rejects_incompatible_labelled_atoms():
    assert smsd.MCS.find("[13C]").match("[12C]").best.atom_count == 1
    options = smsd.MCS.Options()
    options.chemistry.match_isotope = True
    assert not smsd.MCS.find("[13C]").with_options(options).match("[12C]").found


@pytest.mark.parametrize("as_graph", [False, True])
def test_salt_self_mcs_respects_connected_and_full_fragment_modes(as_graph):
    molecule = "[Na+].[O-]C(=O)C"
    if as_graph:
        molecule = smsd.parse_smiles(molecule)
    query = smsd.MCS.find(molecule)
    assert query.match(molecule).best.atom_count == 4
    fragments = query.with_fragments().match(molecule)
    assert fragments.best.atom_count == 5
    assert sorted(map(len, fragments.best.fragments)) == [1, 4]


@pytest.mark.parametrize("query,target,field,value,expected", [
    ("C", "N", "match_atom_type", False, 1),
    ("[NH4+]", "N", "match_formal_charge", True, 0),
    ("C=C", "CC", "match_bond_order", smsd.BondOrderMode.LOOSE, 2),
    ("CC", "CCC", "ring_matches_ring_only", True, 2),
    ("CC", "CCC", "complete_rings_only", True, 2),
])
def test_chemistry_options_preserved_across_repeated_searches(query, target, field, value, expected):
    options = smsd.MCS.Options()
    setattr(options.chemistry, field, value)
    prepared = smsd.MCS.find(query).with_options(options)
    for result in [prepared.match(target), *prepared.match_many([target, target], num_threads=2)]:
        assert (result.best.atom_count if result.found else 0) == expected
        assert result.optimal


@pytest.mark.parametrize("field,value", [
    ("match_isotope", True), ("complete_rings_only", True),
    ("use_chirality", True), ("use_bond_stereo", True),
    ("tautomer_aware", True), ("match_atom_type", False),
    ("match_formal_charge", True), ("ring_matches_ring_only", True),
    ("match_bond_order", smsd.BondOrderMode.STRICT),
    ("match_bond_order", smsd.BondOrderMode.LOOSE),
    ("match_bond_order", smsd.BondOrderMode.ANY),
    ("aromaticity_mode", smsd.AromaticityMode.STRICT),
    ("aromaticity_mode", smsd.AromaticityMode.FLEXIBLE),
    ("ring_fusion_mode", smsd.RingFusionMode.STRICT),
])
def test_prepared_query_copies_chemistry_options(field, value):
    options = smsd.MCS.Options()
    setattr(options.chemistry, field, value)
    prepared = smsd.MCS.find("CC").with_options(options)
    options.chemistry.match_atom_type = False
    assert prepared.match("CCC").best.atom_count == 2
    if field != "match_atom_type":
        assert not prepared.match("NNN").found


def test_unknown_options_are_rejected():
    with pytest.raises(AttributeError, match="unknown_option"):
        smsd.MCS.Options().unknown_option = True


def test_default_bond_and_ring_options_keep_ordinary_matches():
    assert smsd.MCS.find("CC").match("CCC").best.atom_count == 2
    options = smsd.MCS.Options()
    options.chemistry.match_bond_order = smsd.BondOrderMode.ANY
    assert smsd.MCS.find("C=C").with_options(options).with_timeout(1000).match("CCC").best.atom_count == 2
    options.chemistry.ring_matches_ring_only = True
    assert smsd.MCS.find("C1CC1").with_options(options).match("CC1CC1").best.atom_count == 3


def test_prepared_queries_keep_independent_objectives():
    prepared = smsd.MCS.find("CCC")
    edge = prepared.with_objective(smsd.MCS.Objective.MCES)
    assert prepared.match("C1CC1").best.atom_count == 2
    assert edge.match("C1CC1").best.atom_count == 3
    assert prepared.match("C1CC1").best.atom_count == 2


@pytest.mark.parametrize("field", ["objective", "match_bond_order"])
def test_invalid_objective_or_bond_mode_has_clear_error(field):
    with pytest.raises(TypeError):
        if field == "objective":
            smsd.MCS.Options().objective = "typo"
        else:
            smsd.ChemOptions().match_bond_order = "typo"


@pytest.mark.parametrize("objective", [smsd.MCS.Objective.MCIS, smsd.MCS.Objective.MCES])
def test_complete_selected_rings_survive_reverse_search(objective):
    # A fusion atom may be retained when one complete relevant ring covers it.
    # Here a full six-membered ring and the N-N pendant are common.
    options = smsd.MCS.Options()
    options.chemistry.complete_rings_only = True
    options.objective = objective
    molecules = [smsd.parse_smiles("NNc1ccc2ccccc2c1"), smsd.parse_smiles("NNc1ccccc1")]
    for query, target in [molecules, list(reversed(molecules))]:
        result = smsd.MCS.find(query).with_options(options).with_timeout(1000).match(target)
        assert result.best.atom_count == 8 and result.best.bond_count == 8
        assert result.optimal and smsd.MCS.validate_mapping(query, target, result.best, options) == []
        for first, second, target_first, target_second in result.best.bond_mapping:
            assert query.has_bond(first, second) and target.has_bond(target_first, target_second)
        for offset in [0, 2]:
            adjacency = {}
            for bond in result.best.bond_mapping:
                first, second = bond[offset:offset + 2]
                adjacency.setdefault(first, set()).add(second)
                adjacency.setdefault(second, set()).add(first)
            def closes_cycle(path):
                if len(path) == 6:
                    return path[0] in adjacency[path[-1]]
                return any(closes_cycle(path + [neighbor]) for neighbor in adjacency[path[-1]] if neighbor not in path)
            assert any(closes_cycle([atom]) for atom in adjacency)

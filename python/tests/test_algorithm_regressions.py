# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
"""Regression checks for MCS option preservation across Python strategies."""

import pytest

import smsd


def test_auto_induced_mcs_preserves_requested_topology():
    # A three-atom path embeds in a triangle only with non-induced matching.
    assert len(smsd.find_mcs("CCC", "C1CC1")) == 3
    induced = smsd.find_mcs("CCC", "C1CC1", induced=True)
    assert len(induced) == 2
    assert induced == smsd.find_mcs("CCC", "C1CC1", induced=True, strategy="native")


def test_auto_isotope_matching_rejects_incompatible_labelled_atoms():
    assert len(smsd.find_mcs("[13C]", "[12C]")) == 1
    assert smsd.find_mcs("[13C]", "[12C]", match_isotope=True) == {}
    assert smsd.find_mcs("[13C]", "[12C]", match_isotope=True, strategy="native") == {}


@pytest.mark.parametrize("as_graph", [False, True])
def test_salt_self_mcs_respects_connected_and_full_fragment_modes(as_graph):
    molecule = "[Na+].[O-]C(=O)C"
    if as_graph:
        molecule = smsd.parse_smiles(molecule)
    assert len(smsd.find_mcs(molecule, molecule)) == 4
    assert len(smsd.find_mcs(molecule, molecule, connected_only=False)) == 5


@pytest.mark.parametrize(
    "query,target,options",
    [
        ("C", "N", {"match_atom_type": False}),
        ("[NH4+]", "N", {"match_formal_charge": True}),
        ("C=C", "CC", {"match_bond_order": "loose"}),
        ("CC", "CCC", {"max_stage": 1}),
        ("CC", "CCC", {"extra_seeds": False}),
    ],
)
def test_auto_advanced_options_agree_with_native(query, target, options):
    assert smsd.find_mcs(query, target, **options) == smsd.find_mcs(
        query, target, strategy="native", **options
    )


@pytest.mark.parametrize(
    "options,name",
    [
        ({"induced": True}, "induced"),
        ({"match_isotope": True}, "match_isotope"),
        ({"complete_rings_only": True}, "complete_rings_only"),
        ({"use_chirality": True}, "use_chirality"),
        ({"use_bond_stereo": True}, "use_bond_stereo"),
        ({"tautomer_aware": True}, "tautomer_aware"),
        ({"connected_only": False}, "connected_only"),
        ({"maximize_bonds": True}, "maximize_bonds"),
        ({"max_stage": 1}, "max_stage"),
        ({"match_atom_type": False}, "match_atom_type"),
        ({"match_formal_charge": True}, "match_formal_charge"),
        ({"match_bond_order": "loose"}, "match_bond_order"),
        ({"max_results": 2}, "max_results"),
        ({"extra_seeds": False}, "extra_seeds"),
    ],
)
def test_explicit_lightweight_rejects_unsupported_options(options, name):
    with pytest.raises(ValueError, match=name):
        smsd.find_mcs("CC", "CCC", strategy="lightweight", **options)


def test_unknown_native_keyword_cannot_be_hidden_by_auto_fast_return():
    with pytest.raises(AttributeError, match="unknown_option"):
        smsd.find_mcs("CC", "CCC", unknown_option=True)


def test_default_and_supported_options_keep_lightweight_fast_path(monkeypatch):
    def unexpected_native(*args, **kwargs):
        raise AssertionError("supported lightweight options unexpectedly used native")

    monkeypatch.setattr(smsd, "_native_find_mcs", unexpected_native)
    assert len(smsd.find_mcs("CC", "CCC")) == 2
    assert len(smsd.find_mcs("C=C", "CCC", match_bond_order="Any", timeout_ms=1000)) == 2
    assert len(smsd.find_mcs("C1CC1", "CC1CC1", ring_matches_ring_only=True)) == 3


def test_auto_advanced_search_skips_lightweight_candidate(monkeypatch):
    calls = []

    def unexpected_lightweight(*args, **kwargs):
        calls.append((args, kwargs))
        raise AssertionError("advanced options must not run the lightweight engine")

    monkeypatch.setattr(smsd, "_find_mcs_light", unexpected_lightweight)
    assert len(smsd.find_mcs("CCC", "C1CC1", induced=True)) == 2
    assert calls == []


@pytest.mark.parametrize("options,name", [({"strategy": "typo"}, "strategy"),
                                          ({"match_bond_order": "typo"}, "match_bond_order")])
def test_invalid_strategy_or_bond_mode_has_clear_error(options, name):
    with pytest.raises(ValueError, match=name):
        smsd.find_mcs("CC", "CCC", **options)


@pytest.mark.parametrize("induced", [False, True])
def test_complete_query_rings_survive_reverse_search(induced):
    # Completing either overlapping naphthalene ring requires ten query
    # carbons. The target has only six, so only the N-N pendant can remain.
    mapping = smsd.find_mcs(
        "NNc1ccc2ccccc2c1", "NNc1ccccc1",
        complete_rings_only=True, induced=induced, timeout_ms=1000,
    )
    assert set(mapping) == {0, 1}
    assert set(mapping.values()) == {0, 1}

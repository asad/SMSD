# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Search entry points apply the same chemistry and fingerprint defaults."""
import pytest
import smsd
from .search_helpers import configured_mcs, atom_map


@pytest.mark.parametrize("entry", ["single", "all", "batch", "corpus", "selected_output", "copied_options"])
def test_tautomer_option_keeps_strict_bond_default(entry):
    query, target = smsd.parse_smiles("CC(=O)C"), smsd.parse_smiles("CC(O)=C")
    expected = len(atom_map(configured_mcs(query, tautomer_aware=True).match(target)))
    assert expected == 2
    search = configured_mcs(query, tautomer_aware=True)
    if entry == "all":
        result = search.match_all(target)
    elif entry == "batch":
        result = search.match_many([target])[0]
    elif entry == "corpus":
        result = search.match_many(smsd.TargetCorpus.from_smiles([target]))[0]
    elif entry == "copied_options":
        result = search.with_timeout(1000).match(target)
    else:
        result = search.match(target)
    size = smsd.parse_smiles(result.best.to_smiles(query)).n if entry == "selected_output" else result.best.atom_count
    assert size == expected


@pytest.mark.parametrize("entry", ["result", "selected_output"])
def test_result_helpers_accept_chemistry_options(entry):
    result = configured_mcs("CC(=O)C", tautomer_aware=True, match_bond_order="any").match("CC(O)=C")
    assert (smsd.parse_smiles(result.best.to_smiles("CC(=O)C")).n if entry == "selected_output" else result.best.atom_count) == 4


@pytest.mark.parametrize("name", ["strict", "STRICT", "compat-substruct", "Compat-FMCS"])
def test_named_profiles_match_documented_contract(name):
    options = smsd.ChemOptions.profile(name)
    assert not options.match_formal_charge
    assert not options.ring_matches_ring_only
    expected = smsd.BondOrderMode.STRICT if name.lower() == "strict" else smsd.BondOrderMode.LOOSE
    assert options.match_bond_order == expected
    aromaticity = smsd.AromaticityMode.STRICT if name.lower() == "strict" else smsd.AromaticityMode.FLEXIBLE
    assert options.aromaticity_mode == aromaticity


def test_unknown_profile_raises():
    with pytest.raises(ValueError):
        smsd.ChemOptions.profile("strcit")


def test_mcs_fingerprint_default_label_golden():
    assert smsd.mcs_fingerprint(smsd.parse_smiles("C")) == [1973]


def test_batch_similarity_uses_circular_fingerprints():
    query, targets = "c1ccccc1", ["c1ccc(O)cc1", "CCCC", "c1ccncc1"]
    query_fp = smsd.circular_fingerprint(smsd.parse_smiles(query))
    expected = sorted([(index, smsd.tanimoto_coefficient(query_fp, smsd.circular_fingerprint(smsd.parse_smiles(text))))
                       for index, text in enumerate(targets)], key=lambda item: item[1], reverse=True)
    assert smsd.batch_similarity(query, targets) == expected


@pytest.mark.parametrize("mode", ["path", "mcs"])
def test_path_fingerprint_convenience_defaults(mode):
    molecule = smsd.parse_smiles("CCCCCCCC")
    expected = smsd.fingerprint(molecule, kind=mode)
    assert smsd.fingerprint_from_smiles("CCCCCCCC", mode=mode) == expected

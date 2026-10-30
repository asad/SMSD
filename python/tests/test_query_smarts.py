# SPDX-License-Identifier: Apache-2.0
"""A portable MCS query retains its original search result and selected bonds."""
import pytest
import smsd


def test_mcs_result_carries_query_and_proof():
    exported = smsd.MCS.find("CCO").match_with_smarts("CCN")
    assert exported.result.found and exported.result.optimal
    assert smsd.Substructure.smarts(exported.smarts).match("CCCl").found


def test_immediate_stop_has_no_pattern_or_optimum():
    exported = smsd.MCS.find("CCO").with_timeout(0).match_with_smarts("CCN")
    assert exported.smarts == ""
    assert not exported.result.found and not exported.result.optimal


def test_export_keeps_prepared_charge_options():
    options = smsd.MCS.Options()
    options.chemistry.match_formal_charge = True
    prepared = smsd.MCS.find("C").with_options(options)
    exported = prepared.match_with_smarts("C")
    assert smsd.Substructure.smarts(exported.smarts).match("C").found
    assert not smsd.Substructure.smarts(exported.smarts).match("[CH3+]").found


def test_unsupported_export_rules_are_explicit():
    options = smsd.MCS.Options()
    options.chemistry.use_chirality = True
    with pytest.raises(ValueError, match="SMARTS export"):
        smsd.MCS.find("CCO").with_options(options).match_with_smarts("CCO")

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""MCS and substructure SVG highlights keep correspondence and proof information."""

from xml.etree import ElementTree

import pytest
import smsd
from smsd.highlight import mcs_svg, substructure_svg


def _halo_fills(root):
    return {node.get("fill") for layer in root.iter() if layer.get("class") == "highlight"
            for node in layer if node.get("fill") != "none"}


def test_mcs_svg_highlights_both_molecules_and_keeps_status():
    query, target = "c1ccccc1O", "c1ccccc1N"
    result = smsd.MCS.find(query).match(target)
    root = ElementTree.fromstring(mcs_svg(query, target, result, width=600, height=300))
    assert root.attrib["data-search-status"] == result.status.name
    assert root.attrib["data-optimal"] == "true"
    assert root.attrib["data-complete"] == "true"
    assert "6 mapped atoms" in root.find("{*}title").text
    labels = [node for node in root.findall(".//{*}text") if node.get("class") == "map"]
    assert len(labels) == 12


def test_timed_out_highlight_does_not_claim_a_proved_mcs():
    result = smsd.MCS.find("CCO").with_timeout(0).match("CCCO")
    root = ElementTree.fromstring(mcs_svg("CCO", "CCCO", result))
    assert root.attrib["data-search-status"] == "TIMED_OUT"
    assert root.attrib["data-optimal"] == root.attrib["data-complete"] == "false"


def test_selected_edge_highlight_draws_only_carried_bonds():
    # Propane maps onto two of cyclopropane's bonds; the third has no partner bond and stays plain.
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    result = smsd.MCS.find("CCC").with_options(options).match("C1CC1")
    assert result.best.atom_count == 3 and result.best.bond_count == 2
    svg = mcs_svg("CCC", "C1CC1", result)
    assert svg.count('class="bond-highlight"') == 4


def test_selected_edge_highlight_does_not_colour_an_extra_bond():
    # Cyclopropene's double bond is left out of the match, yet its atoms map onto a bond of
    # cyclopropane, so the renderer would highlight it.
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    result = smsd.MCS.find("C1=CC1").with_options(options).match("C1CC1")
    assert result.best.atom_count == 3 and result.best.bond_count == 2
    with pytest.raises(ValueError, match="selected-edge MCS"):
        mcs_svg("C1=CC1", "C1CC1", result)


def test_mcs_svg_draws_a_mint_halo_on_request():
    result = smsd.MCS.find("c1ccccc1O").match("c1ccccc1N")
    assert _halo_fills(ElementTree.fromstring(mcs_svg("c1ccccc1O", "c1ccccc1N", result))) == {"#F9E6B3", "#F0C246"}
    root = ElementTree.fromstring(mcs_svg("c1ccccc1O", "c1ccccc1N", result, highlight="mint"))
    assert _halo_fills(root) == {"#B8F6E9", "#4CE1C1"}


@pytest.mark.parametrize("query, search", [
    ("S(=O)(=O)N", smsd.Substructure.find), ("[SX4](=O)(=O)[NX3]", smsd.Substructure.smarts)])
def test_substructure_svg_highlights_the_hit_in_mint_and_keeps_status(query, search):
    target = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    result = search(query).match(target)
    root = ElementTree.fromstring(substructure_svg(query, target, result, width=400, height=300))
    assert root.attrib["data-search-status"] == result.status.name
    assert root.attrib["data-optimal"] == root.attrib["data-complete"] == "true"
    assert "4 hit atoms" in root.find("{*}title").text
    assert _halo_fills(root) == {"#B8F6E9", "#4CE1C1"}
    labels = [node for node in root.findall(".//{*}text") if node.get("class") == "map"]
    assert len(labels) == 4
    gold = ElementTree.fromstring(substructure_svg(query, target, result, highlight="gold"))
    assert _halo_fills(gold) == {"#F9E6B3", "#F0C246"}


def test_substructure_svg_without_a_hit_draws_the_plain_target():
    result = smsd.Substructure.find("CCN").match("CCC")
    root = ElementTree.fromstring(substructure_svg("CCN", "CCC", result))
    assert "0 hit atoms" in root.find("{*}title").text
    assert not _halo_fills(root)


def test_substructure_svg_does_not_colour_a_bond_the_query_lacks():
    result = smsd.Substructure.find("CCCCCC").match("C1CCCCC1")
    with pytest.raises(ValueError, match="share a bond that the query lacks"):
        substructure_svg("CCCCCC", "C1CCCCC1", result)


def test_substructure_svg_rejects_mismatched_inputs():
    result = smsd.Substructure.find("CO").match("CCO")
    with pytest.raises(TypeError, match="SearchResult"):
        substructure_svg("CO", "CCO", result.best)
    with pytest.raises(ValueError, match="does not belong"):
        substructure_svg("CO", "C", result)
    with pytest.raises(ValueError):
        substructure_svg("C(", "CCO", result)
    with pytest.raises(ValueError, match="highlight must be 'gold' or 'mint'"):
        substructure_svg("CO", "CCO", result, highlight="red")

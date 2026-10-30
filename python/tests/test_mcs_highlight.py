# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""MCS SVG highlights keep correspondence and proof information."""

from xml.etree import ElementTree

import pytest
import smsd
from smsd.highlight import mcs_svg


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


def test_selected_edge_highlight_does_not_colour_an_extra_bond():
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    result = smsd.MCS.find("CCC").with_options(options).match("C1CC1")
    assert result.best.atom_count == 3 and result.best.bond_count == 2
    with pytest.raises(ValueError, match="selected-edge MCS"):
        mcs_svg("CCC", "C1CC1", result)

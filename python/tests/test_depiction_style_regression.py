# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""SVG structure and mapped-atom labelling regression controls."""
import xml.etree.ElementTree as ET
import pytest
import smsd


def _render(pair, **options):
    if pair:
        return smsd.depict_pair("CCO", "CCO", {0: 0, 1: 1, 2: 2}, **options)
    return smsd.depict_svg("CCO", **options)


@pytest.mark.parametrize("pair", [False, True])
@pytest.mark.parametrize("family", ["a}</style><script>alert(1)</script><style>", "Arial; fill:red", "\x01<>{};"])
def test_font_family_preserves_svg_and_css_structure(pair, family):
    root = ET.fromstring(_render(pair, font_family=family))
    assert not root.findall(".//{*}script")
    style = root.find(".//{*}style").text
    declaration = style.split("font-family:", 1)[1].split("}", 1)[0]
    assert declaration.count(";") == 1
    assert declaration.strip() != ";"


@pytest.mark.parametrize("pair", [False, True])
def test_ordinary_font_family_list_is_preserved(pair):
    family = "'Noto Sans', Helvetica, sans-serif"
    assert family in _render(pair, font_family=family)


@pytest.mark.parametrize("pair", [False, True])
def test_map_numbers_cover_carbon_and_heteroatoms(pair):
    if pair:
        svg = smsd.depict_pair("c1ccccc1", "c1ccccc1O", {atom: atom for atom in range(6)})
        expected = 12
    else:
        svg = smsd.depict_mapping("c1ccccc1O", {atom: atom for atom in range(7)})
        expected = 7
    root = ET.fromstring(svg)
    labels = [node for node in root.findall(".//{*}text") if node.get("class") == "map"]
    assert len(labels) == expected
    assert sorted(int(label.text) for label in labels) == sorted(list(range(1, 7)) * 2 if pair else range(1, 8))

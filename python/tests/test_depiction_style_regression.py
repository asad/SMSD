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


@pytest.mark.parametrize("highlight, tint, edge", [
    (None, "#F9E6B3", "#F0C246"), ("gold", "#F9E6B3", "#F0C246"), ("mint", "#B8F6E9", "#4CE1C1")])
@pytest.mark.parametrize("pair", [False, True])
def test_matched_region_is_one_flat_halo(pair, highlight, tint, edge):
    options = {} if highlight is None else {"highlight": highlight}
    if pair:
        svg = smsd.depict_pair("c1ccccc1O", "c1ccccc1N", {atom: atom for atom in range(6)}, **options)
    else:
        svg = smsd.depict_mapping("c1ccccc1O", {atom: atom for atom in range(7)}, **options)
    root = ET.fromstring(svg)
    layers = [node for node in root.iter() if node.get("class") == "highlight"]
    assert len(layers) == (2 if pair else 1)
    colours = {node.get("stroke") if node.get("fill") == "none" else node.get("fill")
               for layer in layers for node in layer}
    assert colours == {tint, edge}
    assert "rgba(" not in svg and "opacity" not in svg


def test_highlight_keeps_keyword_overrides_and_given_options():
    options = smsd.DepictOptions()
    svg = smsd.depict_mapping("c1ccccc1O", {0: 0}, options, highlight="mint", width=320, height=200)
    assert 'width="320"' in svg and "#B8F6E9" in svg and "#F9E6B3" not in svg
    assert "#F9E6B3" in smsd.depict_mapping("c1ccccc1O", {0: 0}, options)


@pytest.mark.parametrize("draw, default, other", [
    ("depict_mcs", (242, 204, 102), "mint"), ("depict_substructure", (78, 232, 199), "gold")])
def test_rdkit_helpers_highlight_in_gold_or_mint(draw, default, other):
    pytest.importorskip("rdkit.Chem.Draw")
    colours = {"gold": (242, 204, 102), "mint": (78, 232, 199)}

    def drawn(**options):
        image = getattr(smsd, draw)("c1ccccc1", "c1ccc(O)cc1", **options).convert("RGB")
        return {colour for _, colour in image.getcolors(image.width * image.height)}

    assert default in drawn() and colours[other] not in drawn()
    assert colours[other] in drawn(highlight=other) and default not in drawn(highlight=other)


@pytest.mark.parametrize("highlight", ["Gold", "teal", "", None, 1, ("gold",)])
@pytest.mark.parametrize("draw", [
    lambda highlight: smsd.depict_mapping("CCO", {0: 0}, highlight=highlight),
    lambda highlight: smsd.depict_pair("CCO", "CCN", {0: 0}, highlight=highlight),
    lambda highlight: smsd.depict_mcs("CCO", "CCN", highlight=highlight),
    lambda highlight: smsd.depict_substructure("CO", "CCO", highlight=highlight),
])
def test_unknown_highlight_names_the_allowed_values(draw, highlight):
    with pytest.raises(ValueError, match="highlight must be 'gold' or 'mint'"):
        draw(highlight)


def test_rdkit_helpers_highlight_only_carried_bonds():
    Chem = pytest.importorskip("rdkit.Chem")
    from smsd._depiction import _carried_bonds
    hexane, cyclohexane = Chem.MolFromSmiles("CCCCCC"), Chem.MolFromSmiles("C1CCCCC1")
    first, second = _carried_bonds(hexane, cyclohexane, {atom: atom for atom in range(6)})
    assert len(first) == len(second) == 5
    assert cyclohexane.GetBondBetweenAtoms(0, 5).GetIdx() not in second

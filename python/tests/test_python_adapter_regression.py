# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
"""Conversions preserve source molecules and depiction atom correspondence."""
import pytest
import smsd
from .search_helpers import configured_mcs, atom_map


@pytest.mark.parametrize("operation", ["kekulize", "dearomatize", "perceive_aromaticity"])
def test_aromaticity_operations_preserve_cached_input(operation):
    chemistry = pytest.importorskip("rdkit.Chem")
    molecule = chemistry.MolFromSmiles("Oc1ccccc1")
    cached = smsd.from_rdkit(molecule)
    before = smsd.to_smiles(cached)
    converted = getattr(smsd, operation)(molecule)
    assert converted is not cached
    assert smsd.to_smiles(smsd.from_rdkit(molecule)) == before
    assert smsd.Substructure.smarts('c1ccccc1').match(molecule).found


@pytest.mark.parametrize("operation", ["kekulize", "dearomatize", "perceive_aromaticity"])
def test_aromaticity_operations_modify_native_graph(operation):
    molecule = smsd.parse_smiles("Oc1ccccc1")
    assert getattr(smsd, operation)(molecule) is molecule


@pytest.mark.parametrize("drawing", ["depict_mcs", "depict_substructure"])
@pytest.mark.parametrize("input_kind", ["graph", "external", "smiles"])
def test_depiction_highlights_the_mapped_elements(monkeypatch, drawing, input_kind):
    chemistry = pytest.importorskip("rdkit.Chem")
    from rdkit.Chem import Draw
    first, second = "OCC", "OCCN"
    if input_kind == "graph":
        first, second = smsd.parse_smiles(first), smsd.parse_smiles(second)
    elif input_kind == "external":
        first, second = chemistry.MolFromSmiles(first), chemistry.MolFromSmiles(second)
    captured = {}
    def capture(molecules, **options):
        captured["elements"] = [sorted(molecule.GetAtomWithIdx(index).GetAtomicNum() for index in atoms)
                                for molecule, atoms in zip(molecules, options["highlightAtomLists"])]
        return captured
    monkeypatch.setattr(Draw, "MolsToGridImage", capture)
    result = getattr(smsd, drawing)(first, second)
    assert result["elements"] == [[6, 6, 8], [6, 6, 8]]


def test_adaptive_search_uses_native_budget():
    result = smsd.MCS.find("CCCC").with_timeout(-1).match("CCCCC")
    assert result.best.atom_count == 4 and result.optimal
    assert result.status == smsd.SearchStatus.COMPLETED


def test_zero_timeout_reports_an_incomplete_search():
    result = smsd.MCS.find("CCCC").with_timeout(0).match("CCCCC")
    assert result.status == smsd.SearchStatus.TIMED_OUT
    assert not result.optimal and not result.complete

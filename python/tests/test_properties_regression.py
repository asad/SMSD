# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Molecule metadata mutations remain visible in graph and SDF output."""
import gc

import pytest
import smsd


def test_property_item_assignment_and_deletion_update_the_graph():
    graph = smsd.parse_smiles("CCO")
    graph.properties = {"Identifier": "initial"}
    view = graph.properties
    view["Identifier"] = "updated"
    view["Source"] = "sample"
    assert dict(graph.properties) == {"Identifier": "updated", "Source": "sample"}
    assert "updated" in smsd.write_sdf_record(graph)
    del view["Source"]
    assert "Source" not in graph.properties
    with pytest.raises(KeyError):
        del view["missing"]


def test_property_view_retains_its_molecule():
    graph = smsd.parse_smiles("C")
    graph.properties = {"Identifier": "sample"}
    view = graph.properties
    del graph
    gc.collect()
    assert view["Identifier"] == "sample"
    view.update({"Source": "archive"})
    assert view.get("Source") == "archive"
    assert view.get("missing") is None
    assert set(view) == {"Identifier", "Source"}
    assert view.copy() == {"Identifier": "sample", "Source": "archive"}


def test_property_replacement_updates_existing_views():
    graph = smsd.parse_smiles("C")
    view = graph.properties
    graph.properties = {"Identifier": "replacement"}
    assert view == {"Identifier": "replacement"}
    assert view == graph.properties
    view.clear()
    assert dict(graph.properties) == {}

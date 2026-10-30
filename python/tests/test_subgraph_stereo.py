# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Stereo kept by extract_subgraph and R-group decomposition, from the fixture also read by the C++ and Java tests."""

import json
from pathlib import Path

import pytest

import smsd


FIXTURE_PATH = (
    Path(__file__).resolve().parents[2]
    / "java/src/test/resources/com/bioinception/smsd/subgraph_stereo.json"
)
with FIXTURE_PATH.open(encoding="utf-8") as handle:
    FIXTURE = json.load(handle)


def has_stereo(smiles):
    return any(mark in smiles for mark in "@/\\")


@pytest.mark.parametrize("case", FIXTURE["subgraphs"], ids=lambda case: case["name"])
def test_extract_subgraph_keeps_defined_stereo(case):
    parent = smsd.parse_smiles(case["smiles"])
    assert has_stereo(smsd.to_smiles(parent))
    fragment = smsd.to_smiles(smsd.extract_subgraph(parent, case["atoms"]))
    assert fragment == smsd.to_smiles(smsd.parse_smiles(case["expected"]))
    assert has_stereo(fragment) == case["stereo"]


@pytest.mark.parametrize("case", FIXTURE["rgroups"], ids=lambda case: case["name"])
def test_r_groups_keep_defined_stereo(case):
    results = smsd.Substructure.decompose_r_groups(case["core"], [case["smiles"]])
    assert len(results) == 1 and results[0]["search"].found
    groups = sorted(smsd.to_smiles(group) for group in results[0]["rgroups"].values())
    assert groups == sorted(smsd.to_smiles(smsd.parse_smiles(smiles)) for smiles in case["expected"])
    assert any(has_stereo(group) for group in groups) == case["stereo"]


def test_atom_order_does_not_change_a_centre():
    parent = smsd.parse_smiles("N[C@@H](C)C(=O)O")
    in_order = smsd.extract_subgraph(parent, [0, 1, 2, 3, 4])
    permuted = smsd.extract_subgraph(parent, [3, 4, 1, 0, 2])
    assert smsd.to_smiles(permuted) == smsd.to_smiles(in_order)
    assert "@" in smsd.to_smiles(permuted)

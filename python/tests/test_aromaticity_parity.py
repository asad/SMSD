# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Shared aromaticity parity corpus for native SMSD."""

import json
from pathlib import Path

import pytest

from smsd import (
    AromaticityModel,
    dearomatize,
    kekulize,
    parse_smiles,
    perceive_aromaticity,
)


CORPUS_PATH = (
    Path(__file__).resolve().parents[2]
    / "src/test/resources/com/bioinception/smsd/aromaticity_parity.json"
)
with CORPUS_PATH.open(encoding="utf-8") as handle:
    AROMATICITY_CASES = json.load(handle)


def aromatic_atom_count(mol):
    return sum(bool(flag) for flag in mol.aromatic)


def ring_atom_count(mol):
    return sum(bool(flag) for flag in mol.ring)


def aromatic_bond_count(mol):
    return sum(
        1
        for i in range(mol.n)
        for j in range(i + 1, mol.n)
        if mol.bond_aromatic(i, j)
    )


def bond_signature(mol):
    return [
        (i, j, mol.bond_order(i, j), mol.bond_aromatic(i, j))
        for i in range(mol.n)
        for j in range(i + 1, mol.n)
        if mol.bond_order(i, j)
    ]


def aromatic_cases():
    return [
        entry
        for entry in AROMATICITY_CASES
        if entry["aromatic_atom_count"] > 0
    ]


@pytest.mark.parametrize("entry", AROMATICITY_CASES, ids=lambda entry: entry["name"])
def test_daylight_like_parity_corpus(entry):
    mol = parse_smiles(entry["kekule"])

    assert aromatic_atom_count(mol) == entry["aromatic_atom_count"]
    assert aromatic_bond_count(mol) == entry["aromatic_bond_count"]
    assert ring_atom_count(mol) == entry["ring_atom_count"]


@pytest.mark.parametrize("entry", aromatic_cases(), ids=lambda entry: entry["name"])
def test_kekulize_round_trip_restores_daylight_like_perception(entry):
    mol = parse_smiles(entry["kekule"])

    kekulized = kekulize(mol)
    assert kekulized is mol
    assert aromatic_atom_count(mol) == 0
    assert aromatic_bond_count(mol) == 0
    assert all(order != 4 for _, _, order, _ in bond_signature(mol))

    perceived = perceive_aromaticity(mol, model=AromaticityModel.DAYLIGHT_LIKE)
    assert perceived is mol
    assert aromatic_atom_count(mol) == entry["aromatic_atom_count"]
    assert aromatic_bond_count(mol) == entry["aromatic_bond_count"]
    assert ring_atom_count(mol) == entry["ring_atom_count"]


@pytest.mark.parametrize("entry", aromatic_cases(), ids=lambda entry: entry["name"])
def test_dearomatize_alias_matches_kekulize(entry):
    kekule_graph = parse_smiles(entry["kekule"])
    dearomatized_graph = parse_smiles(entry["kekule"])

    kekulize(kekule_graph)
    dearomatize(dearomatized_graph)

    assert list(kekule_graph.aromatic) == list(dearomatized_graph.aromatic)
    assert bond_signature(kekule_graph) == bond_signature(dearomatized_graph)

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""MOL blocks give E/Z only to double bonds that can carry it.

A drawn molecule then has the canonical SMILES of the same molecule read as SMILES.
The cases are shared with the C++ and Java tests.
"""

import json
import math
from pathlib import Path

import pytest

import smsd


FIXTURE_PATH = (
    Path(__file__).resolve().parents[2]
    / "java/cli/src/test/resources/com/bioinception/smsd/cli/mol_ez_eligibility.json"
)
with FIXTURE_PATH.open(encoding="utf-8") as handle:
    FIXTURE = json.load(handle)


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"]) for case in FIXTURE["cases"]])
def test_mol_block_gives_the_smiles_configuration(case):
    drawn = smsd.read_mol_block("\n".join(case["molblock"]) + "\n")
    written = smsd.parse_smiles(case["smiles"])
    assert len(smsd.assign_ez(drawn)) == case["ez_bonds"]
    assert smsd.to_smiles(drawn) == smsd.to_smiles(written)


def _v2000(atoms, bonds, tail=""):
    lines = ["drawn", "  SMSD              2D", "",
             f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for symbol, x, y in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{0.:10.4f} {symbol:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for first, second, order in bonds:
        lines.append(f"{first:3d}{second:3d}{order:3d}  0")
    return "\n".join(lines) + "\n" + tail + "M  END\n"


@pytest.mark.parametrize("size,expected", [(5, 0), (6, 0), (7, 0), (8, 1)])
def test_ring_double_bonds_need_eight_atoms(size, expected):
    radius = 0.75 / math.sin(math.pi / size)
    atoms = [("C", radius * math.cos(2 * math.pi * k / size), radius * math.sin(2 * math.pi * k / size))
             for k in range(size)]
    bonds = [(k + 1, (k + 1) % size + 1, 2 if k == 0 else 1) for k in range(size)]
    assert len(smsd.assign_ez(smsd.read_mol_block(_v2000(atoms, bonds)))) == expected


def test_small_ring_configuration_does_not_survive_a_mol_round_trip():
    tagged = smsd.parse_smiles(r"C1/C=C\CCC1")
    for writer in (smsd.write_mol_block, smsd.write_mol_block_v3000):
        restored = smsd.read_mol_block(writer(tagged))
        assert smsd.assign_ez(restored) == {}
        assert smsd.to_smiles(restored) == smsd.to_smiles(smsd.parse_smiles("C1CC=CCC1"))


PROPENE = [("C", 0.0, 0.0), ("C", 1.5, 0.0), ("C", 2.25, 1.299), ("H", -0.75, 1.299), ("H", -0.75, -1.299)]
PROPENE_BONDS = [(1, 2, 2), (2, 3, 1), (1, 4, 1), (1, 5, 1)]


@pytest.mark.parametrize("tail,expected", [
    ("", 0),                              # =CH2 drawn with both hydrogens
    ("M  ISO  1   5   2\n", 1),           # =CHD: the hydrogens differ
    ("M  ISO  2   4   2   5   2\n", 0),   # =CD2: the hydrogens are alike
])
def test_like_hydrogens_at_one_end_give_no_configuration(tail, expected):
    block = _v2000(PROPENE, PROPENE_BONDS, tail)
    graph = smsd.read_mol_block(block)
    assert len(smsd.assign_ez(graph)) == expected
    # assign_ez leaves out bonds whose ends tie, so compare with the bond drawn crossed (either).
    crossed = smsd.read_mol_block(block.replace("  1  2  2  0", "  1  2  2  3", 1))
    unchanged = smsd.to_smiles(graph) == smsd.to_smiles(crossed)
    assert unchanged == (expected == 0)
    if expected == 0:
        assert smsd.canonical_hash(graph) == smsd.canonical_hash(crossed)

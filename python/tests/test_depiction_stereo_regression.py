# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Depictions distinguish the handedness of ordinary stereocentres."""
import pytest
import smsd


@pytest.mark.parametrize("first,second", [
    ("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
    ("F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I"),
    ("[C@H](F)(Cl)Br", "[C@@H](F)(Cl)Br"),
])
def test_depiction_distinguishes_enantiomers(first, second):
    assert smsd.depict_svg(first) != smsd.depict_svg(second)


@pytest.mark.parametrize("text", [
    "c1ccc2ccccc2c1", "c1ccc2[nH]ccc2c1", "Cc1ccc2ccccc2c1C",
    "c1ccccc1CCc1ccccc1", "C1CCCCCCCCCCCCC1", "CC(=O)[O-].[Na+]",
    "[Na+].[Cl-]", "CCCCCC.CCCCCC", "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
])
def test_ring_and_fragment_layout_preserves_geometry(text):
    import math
    molecule = smsd.parse_smiles(text)
    coordinates = smsd.generate_coords_2d(molecule)
    for first in range(molecule.n):
        for second in range(first + 1, molecule.n):
            distance = math.dist(coordinates[first], coordinates[second]) / 1.5
            assert distance > 0.4
            if molecule.bond_order(first, second):
                assert 0.8 < distance < 1.2
    assert smsd.count_crossings(molecule, coordinates) == 0


def test_canonical_orientation_aligns_principal_axis():
    import math
    coordinates = [(0.0, 0.0), (math.cos(0.35), math.sin(0.35))]
    aligned = smsd.canonical_orientation(smsd.parse_smiles("CC"), coordinates)
    assert aligned[1][1] == pytest.approx(aligned[0][1], abs=1e-9)

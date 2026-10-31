# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Relevant-cycle families: the API surface, exact counts, membership and member order."""
import pytest
import smsd


def cyclodextrin(units):
    return "C%101C(O)C(O)C(OC1CO)O" + "C1C(O)C(O)C(OC1CO)O" * (units - 1) + "%10"


def test_surface_is_exported():
    for name in ("RelevantCycleFamily", "RelevantCycleFamilies", "compute_relevant_cycle_families",
                 "relevant_cycle_family_of", "relevant_cycle_family_members"):
        assert name in smsd.__all__ and hasattr(smsd, name)


def test_cubane_has_six_faces():
    result = smsd.compute_relevant_cycle_families("C12C3C4C1C5C2C3C45")
    assert len(result.families) == 6 and result.relevant_cycle_count == 6 and not result.truncated
    assert all(family.size == 4 and family.member_count == 1 for family in result.families)


@pytest.mark.parametrize("units", [6, 7, 34, 97])
def test_cyclodextrin_counts_are_exact_python_integers(units):
    result = smsd.compute_relevant_cycle_families(cyclodextrin(units))
    assert type(result.relevant_cycle_count) is int
    assert result.relevant_cycle_count == 2 ** units + units and not result.truncated
    assert sum(family.member_count for family in result.families) == result.relevant_cycle_count
    assert sum(1 for family in result.families if family.size == 6) == units


def test_odd_and_even_simple_rings():
    five = smsd.compute_relevant_cycle_families("C1CCCC1").families
    six = smsd.compute_relevant_cycle_families("c1ccccc1").families
    assert [(f.size, f.root, f.first, f.second, f.opposite) for f in five] == [(5, 4, 1, 2, -1)]
    assert [(f.size, f.root, f.first, f.second, f.opposite) for f in six] == [(6, 5, 1, 3, 2)]


def test_members_come_in_a_fixed_order_and_belong_to_their_family():
    mol = smsd.parse_smiles(cyclodextrin(6))
    families = smsd.compute_relevant_cycle_families(mol)
    family = max(families.families, key=lambda f: f.member_count)
    members = smsd.relevant_cycle_family_members(mol, family, 1000)
    assert len(members) == family.member_count and members[0] == family.prototype
    assert members == smsd.relevant_cycle_family_members(mol, family, 1000)
    assert smsd.relevant_cycle_family_members(mol, family, 3) == members[:3]
    for member in members:
        index = smsd.relevant_cycle_family_of(mol, families, member)
        assert index >= 0
        found = families.families[index]
        assert (found.size, found.root, found.first, found.second, found.opposite) == \
            (family.size, family.root, family.first, family.second, family.opposite)


def test_a_ring_whose_arc_skips_a_level_is_no_member():
    # r = 7; a1, a2, a3, b1, b2, b3, y = 0 .. 6. The ring r-a1-b1-a2-b2-y-b3-a3 passes a2,
    # which is next to r, at position three of its arc, so it belongs to no family.
    bonds = [(7, 0), (7, 1), (7, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 6), (4, 6), (5, 6)]
    neighbours = [[] for _ in range(8)]
    for a, b in bonds:
        neighbours[a].append(b)
        neighbours[b].append(a)
    mol = smsd.MolGraphBuilder().atom_count(8).atomic_numbers([6] * 8).neighbors(neighbours).build(False)
    families = smsd.compute_relevant_cycle_families(mol)
    assert smsd.relevant_cycle_family_of(mol, families, [7, 0, 3, 1, 4, 6, 5, 2]) == -1
    for index, family in enumerate(families.families):
        for member in smsd.relevant_cycle_family_members(mol, family, 100):
            assert smsd.relevant_cycle_family_of(mol, families, member) == index


def test_invalid_rings_are_no_members():
    mol = smsd.parse_smiles("c1ccccc1")
    families = smsd.compute_relevant_cycle_families(mol)
    for ring in ([], [0, 1], [0, 1, 3], [0, 1, 2, 3, 4, 5, 0], [0, 1, 2, 3, 4, 99]):
        assert smsd.relevant_cycle_family_of(mol, families, ring) == -1

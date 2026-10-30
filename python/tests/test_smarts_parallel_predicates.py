# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Every predicate on a SMARTS bond constrains its embedding."""
import pytest
import smsd


@pytest.mark.parametrize("pattern,target,expected,atoms", [
    ("C1=C-1", "C1=CC1", False, 2),
    ("[$(C1=C-1)]", "C1=CC1", False, 1),
    ("[!$(C1=C-1)]", "C1=CC1", True, 1),
    ("[#6,$(C1=C-1)]", "C1=CC1", True, 1),
    ("C11", "C1=CC1", False, 1),
    ("[!$(C11)]", "C1=CC1", True, 1),
    ("C1-C-1", "CC", True, 2),
    (r"F/?&\?C=C/F", "FC=CF", True, 4),
    (r"F/?&\?C=C/F", "F/C=C/F", False, 4),
    ("F/?C=C/F", "FC=CF", True, 4),
    ("F/C=C/F", "FC=CF", False, 4),
])
def test_parallel_and_optional_bond_predicates(pattern, target, expected, atoms):
    prepared = smsd.Substructure.smarts(pattern).with_timeout(10000)
    molecule = smsd.parse_smiles(target)
    for result in (prepared.match(molecule), prepared.match_all(molecule)):
        assert result.status == smsd.SearchStatus.COMPLETED
        assert result.optimal and result.complete
        assert result.found is expected
        assert result.score == (atoms if expected else 0)


def test_parallel_predicates_select_one_edge_per_embedding():
    result = smsd.Substructure.smarts("C1-C-1").with_timeout(10000).match_all("CC")
    assert result.complete and result.optimal
    assert len(result.matches) == 2
    assert all(match.atom_count == 2 and match.bond_count == 1 for match in result.matches)

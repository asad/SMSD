# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Stereo matching follows mapped ligand order rather than whole-molecule labels."""

import pytest
import smsd
from .search_helpers import configured_mcs, atom_map


@pytest.mark.parametrize("query_text,target_text,chirality,expected", [
    ("C[C@H](O)CC", "ClC[C@H](O)CC", True, True),
    ("C[C@@H](O)CC", "ClC[C@H](O)CC", True, False),
    ("N[C@@H](C)C(=O)O", "N[C@@H](CS)C(=O)O", True, True),
    ("N[C@@H](C)C(=O)O", "N[C@@](F)(C)C(=O)O", True, True),
    ("N[C@@H](C)C(=O)O", "N[C@](F)(C)C(=O)O", True, False),
    ("C/C=C/C", "Cl/C(C)=C/C", False, False),
    ("C/C=C\\C", "Cl/C(C)=C/C", False, True),
    ("C/C(Cl)=C/C", "Cl/C(C)=C\\C", False, True),
    ("F/C=C/Cl", "F/C=C\\Cl", False, False),
    # An embedding needs specified target stereo for specified query stereo, as SMARTS does.
    ("C[C@H](O)CC", "CC(O)CC", True, False),
    ("CC(O)CC", "C[C@H](O)CC", True, True),
    ("N[C@@H](C)C(=O)O", "NC(C)C(=O)O", True, False),
    ("F/C=C/F", "FC=CF", False, False),
    ("FC=CF", "F/C=C/F", False, True),
])
def test_molecular_stereo_entry_points(query_text, target_text, chirality, expected):
    query, target = smsd.parse_smiles(query_text), smsd.parse_smiles(target_text)
    options = smsd.Substructure.Options()
    options.chemistry.use_chirality = chirality
    options.chemistry.use_bond_stereo = not chirality
    search = smsd.Substructure.find(query).with_options(options)
    assert search.match(target).found is expected
    assert bool(search.match_all(target).matches) is expected


def test_mcs_uses_mapped_stereochemistry():
    assert len(atom_map(configured_mcs('C[C@H](O)CC', use_chirality=True).match('ClC[C@H](O)CC'))) == 5
    # A common-subgraph search keeps an unspecified target as a wildcard.
    assert len(atom_map(configured_mcs('C[C@H](O)CC', use_chirality=True).match('CC(O)CC'))) == 5
    assert len(atom_map(configured_mcs('C/C=C/C', use_bond_stereo=True).match('Cl/C(C)=C/C'))) < 4
    assert len(atom_map(configured_mcs('C/C=C\\C', use_bond_stereo=True).match('Cl/C(C)=C/C'))) == 4


@pytest.mark.parametrize("pattern,target_text,expected", [
    ("[C@@H](F)(Cl)Br", "[C@H](Cl)(F)Br", True),
    ("[C@@H](F)(Cl)Br", "[C@@H](Cl)(F)Br", False),
    ("C[C@H]1CCCCO1", "C[C@@H]1OCCCC1", True),
    ("[C@@H,N](F)(Cl)Br", "[C@H](Cl)(F)Br", True),
    ("[C@@H,N](F)(Cl)Br", "N(F)(Cl)Br", True),
    ("[C;!@@;H1](F)(Cl)Br", "[C@H](F)(Cl)Br", True),
    ("[C;!@@;H1](F)(Cl)Br", "[C@@H](F)(Cl)Br", False),
    ("[C;!@@;H1](F)(Cl)Br", "[C@@H](Cl)(F)Br", True),
    ("[C;!@@;H1](F)(Cl)Br", "[C@H](Cl)(F)Br", False),
    ("[C@?H](F)(Cl)Br", "C(F)(Cl)Br", True),
    ("[C@@H](F)(Cl)Br", "C(F)(Cl)Br", False),
    ("[C;!@@;H1](F)(Cl)Br", "C(F)(Cl)Br", True),
    ("F/C=C/F", "F/C=C/F", True),
    ("F/C=C/F", "F/C=C\\F", False),
    ("F/C=C/F", "FC=CF", False),
    ("C/C", "CC", True),
    ("C/C", "C=C", False),
    ("C/C=C/C", "Cl/C(C)=C/C", False),
    ("C/C=C\\C", "Cl/C(C)=C/C", True),
    ("F/,-C=C/F", "F/C=C/F", True),
    ("F/,-C=C/F", "F/C=C\\F", True),
    ("F/,-C=C/F", "FC=CF", True),
    ("F!/C=C/F", "FC=CF", False),
    ("F!/C=C/F", "F/C=C\\F", True),
    ("F!/C=C/F", "F/C=C/F", False),
    ("[$([C@@H](F)(Cl)Br)]", "[C@H](Cl)(F)Br", True),
    ("[$([C@@H](F)(Cl)Br)]", "[C@@H](Cl)(F)Br", False),
    ("C/C=C(/F)\\C(Cl)=C/C", "C\\C=C(\\F)/C(Cl)=C\\C", True),
    ("C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C", False),
])
def test_compiled_smarts_stereo_entry_points(pattern, target_text, expected):
    query = smsd.Substructure.smarts(pattern)
    target = smsd.parse_smiles(target_text)
    assert query.match(target).found is expected
    assert bool(query.match_all(target).matches) is expected

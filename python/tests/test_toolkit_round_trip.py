# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Chemical identity and correspondence across the optional molecule bridges."""

import pytest
import smsd


MOLECULES = [
    "CC(=O)N", "Oc1ccccc1", "c1cc[nH]c1", "[NH3+]CC(=O)[O-]",
    "[13CH3]CO", "N[C@@H](C)C(=O)O", "C/C=C/C", "C/C=C\\C",
    "[Na+].CC(=O)[O-]",
]


@pytest.mark.parametrize("smiles", MOLECULES)
def test_rdkit_round_trip_preserves_chemical_identity(smiles):
    chemistry = pytest.importorskip("rdkit.Chem")
    source = chemistry.MolFromSmiles(smiles)
    expected = chemistry.MolToSmiles(source, isomericSmiles=True)
    source.SetProp("sample", "source molecule")
    graph = smsd.from_rdkit(source, use_cache=False)
    restored = smsd.to_rdkit(graph)

    assert chemistry.MolToSmiles(restored, isomericSmiles=True) == expected
    assert restored.GetNumAtoms() == source.GetNumAtoms()
    assert restored.GetNumBonds() == source.GetNumBonds()
    assert chemistry.GetFormalCharge(restored) == chemistry.GetFormalCharge(source)
    assert len(chemistry.GetMolFrags(restored)) == len(chemistry.GetMolFrags(source))
    assert sorted((atom.GetAtomicNum(), atom.GetIsotope(), atom.GetIsAromatic())
                  for atom in restored.GetAtoms()) == sorted(
        (atom.GetAtomicNum(), atom.GetIsotope(), atom.GetIsAromatic())
        for atom in source.GetAtoms())
    assert chemistry.MolToSmiles(source, isomericSmiles=True) == expected
    assert source.GetProp("sample") == "source molecule"
    assert not restored.HasProp("sample")

    native_source = smsd.parse_smiles(smiles)
    native_restored = smsd.from_rdkit(smsd.to_rdkit(native_source), use_cache=False)
    assert smsd.same_canonical_graph(native_source, native_restored)


def test_rdkit_round_trip_preserves_input_correspondence():
    chemistry = pytest.importorskip("rdkit.Chem")
    source = chemistry.RenumberAtoms(chemistry.MolFromSmiles("NCCO"), [3, 1, 0, 2])
    graph = smsd.from_rdkit(source, use_cache=False)
    restored = smsd.to_rdkit(graph)
    result = smsd.Substructure.find(source).match(restored)
    assert result.complete and result.optimal and result.found
    assert len(result.best.atom_mapping) == source.GetNumAtoms()
    for query_atom, target_atom in result.best.atom_mapping.items():
        assert source.GetAtomWithIdx(query_atom).GetAtomicNum() == restored.GetAtomWithIdx(target_atom).GetAtomicNum()
    for bond in source.GetBonds():
        target_bond = restored.GetBondBetweenAtoms(
            result.best.atom_mapping[bond.GetBeginAtomIdx()],
            result.best.atom_mapping[bond.GetEndAtomIdx()])
        assert target_bond is not None and target_bond.GetBondType() == bond.GetBondType()

    index_map = smsd.get_index_map(graph)
    assert sorted(index_map) == list(range(source.GetNumAtoms()))
    for graph_index, source_index in enumerate(index_map):
        assert graph.atomic_num[graph_index] == source.GetAtomWithIdx(source_index).GetAtomicNum()


def test_rdkit_explicit_hydrogen_import_and_export_normalisation():
    chemistry = pytest.importorskip("rdkit.Chem")
    source = chemistry.AddHs(chemistry.MolFromSmiles("C[C@H](O)F"))
    graph = smsd.from_rdkit(source, use_cache=False)
    assert graph.n == source.GetNumAtoms()
    index_map = smsd.get_index_map(graph)
    assert sorted(index_map) == list(range(source.GetNumAtoms()))
    for graph_index, source_index in enumerate(index_map):
        assert graph.atomic_num[graph_index] == source.GetAtomWithIdx(source_index).GetAtomicNum()
    restored = smsd.to_rdkit(graph)
    expected = chemistry.RemoveHs(source)
    assert restored.GetNumAtoms() == expected.GetNumAtoms()
    assert chemistry.MolToSmiles(restored) == chemistry.MolToSmiles(expected)
    assert source.GetNumAtoms() > restored.GetNumAtoms()


@pytest.mark.parametrize("smiles", MOLECULES)
@pytest.mark.parametrize("input_kind", ["pybel", "obmol"])
def test_openbabel_round_trip_preserves_chemical_identity(smiles, input_kind):
    pybel = pytest.importorskip("openbabel.pybel")
    source = pybel.readstring("smi", smiles)
    expected = source.write("can").split()[0]
    source.title = "source molecule"
    graph = smsd.from_openbabel(source if input_kind == "pybel" else source.OBMol)
    restored = smsd.to_openbabel(graph)

    assert restored.write("can").split()[0] == expected
    assert restored.OBMol.NumAtoms() == source.OBMol.NumAtoms()
    assert restored.OBMol.NumBonds() == source.OBMol.NumBonds()
    assert restored.OBMol.GetTotalCharge() == source.OBMol.GetTotalCharge()
    assert restored.formula == source.formula
    assert sorted((atom.atomicnum, atom.OBAtom.GetIsotope(), atom.OBAtom.IsAromatic())
                  for atom in restored.atoms) == sorted(
        (atom.atomicnum, atom.OBAtom.GetIsotope(), atom.OBAtom.IsAromatic())
        for atom in source.atoms)
    assert source.write("can").split()[0] == expected
    assert source.title == "source molecule"
    assert restored.title != source.title
    assert smsd.get_index_map(graph) is None

    native_source = smsd.parse_smiles(smiles)
    native_restored = smsd.from_openbabel(smsd.to_openbabel(native_source))
    assert smsd.same_canonical_graph(native_source, native_restored)


@pytest.mark.parametrize("input_kind", ["pybel", "obmol"])
def test_openbabel_explicit_hydrogens_preserve_chemistry(input_kind):
    pybel = pytest.importorskip("openbabel.pybel")
    source = pybel.readstring("smi", "C[C@H](O)F")
    source.addh()
    expected = source.write("can").split()[0]
    atom_count = source.OBMol.NumAtoms()
    graph = smsd.from_openbabel(source if input_kind == "pybel" else source.OBMol)
    restored = smsd.to_openbabel(graph)
    assert restored.write("can").split()[0] == expected
    assert restored.formula == source.formula
    assert restored.OBMol.NumHvyAtoms() == source.OBMol.NumHvyAtoms()
    assert source.OBMol.NumAtoms() == atom_count


@pytest.mark.parametrize("first,second", [
    ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
    ("C/C=C/C", "C/C=C\\C"),
    ("[13CH3]CO", "CCO"),
    ("CC(=O)[O-]", "CC(=O)O"),
    ("c1ccccc1", "C1CCCCC1"),
    ("CCO", "CCN"),
])
@pytest.mark.parametrize("toolkit", ["rdkit", "openbabel"])
def test_round_trip_does_not_collapse_distinct_chemistry(first, second, toolkit):
    if toolkit == "rdkit":
        chemistry = pytest.importorskip("rdkit.Chem")

        def round_trip(smiles):
            source = chemistry.MolFromSmiles(smiles)
            restored = smsd.to_rdkit(smsd.from_rdkit(source, use_cache=False))
            return chemistry.MolToSmiles(restored, isomericSmiles=True)
    else:
        pybel = pytest.importorskip("openbabel.pybel")

        def round_trip(smiles):
            source = pybel.readstring("smi", smiles)
            restored = smsd.to_openbabel(smsd.from_openbabel(source))
            return restored.write("can").split()[0]

    assert round_trip(first) != round_trip(second)

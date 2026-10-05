"""Benchmark validity and corpus-generation regressions, independent of timings."""
from pathlib import Path
import random
import sys

import pytest

pytest.importorskip("rdkit")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "benchmarks"))
from rdkit import Chem
from generate_dalke_pairs import main as generate_pairs, nearest_neighbor_index, random_pair_indices
from mcs_protocol import graph_from_rdkit, rdkit_parameters, summarize, validate_mapping


def test_seed_applies_before_pool_subsampling(tmp_path):
    source = tmp_path / "pool.smi"
    source.write_text("\n".join("C"*n + f" molecule_{n}" for n in range(5, 15)))
    args = ["--input", str(source), "--max-molecules", "3", "--pairs", "2"]
    generate_pairs(args + ["--output-dir", str(tmp_path / "first")])
    # Global RNG activity must not alter either the subsample or its derived pairs.
    random.random(); random.sample(range(100), 10)
    generate_pairs(args + ["--output-dir", str(tmp_path / "second")])
    for name in ("dalke_random_pairs.tsv", "dalke_nn_pairs.tsv"):
        assert (tmp_path / "first" / name).read_bytes() == (tmp_path / "second" / name).read_bytes()


def test_nearest_neighbor_explicitly_excludes_query_under_fingerprint_ties():
    assert nearest_neighbor_index(1, [1.0, 1.0, 0.5]) == 0
    assert nearest_neighbor_index(0, [1.0, 1.0, 1.0]) == 1


def test_pair_generation_rejects_impossible_size_and_deduplicates():
    with pytest.raises(ValueError): random_pair_indices(3, 4, random.Random(42))
    pairs = random_pair_indices(3, 3, random.Random(42))
    assert len({tuple(sorted(pair)) for pair in pairs}) == 3


def test_graph_conversion_preserves_indices_and_rdkit_aromaticity():
    mol = Chem.MolFromSmiles("[O-]c1ccccc1")
    graph = graph_from_rdkit(mol)
    assert graph.n == mol.GetNumAtoms()
    assert graph.atomic_num == [a.GetAtomicNum() for a in mol.GetAtoms()]
    assert graph.formal_charge == [a.GetFormalCharge() for a in mol.GetAtoms()]
    assert [bool(x) for x in graph.aromatic] == [a.GetIsAromatic() for a in mol.GetAtoms()]
    for bond in mol.GetBonds():
        assert graph.has_bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        assert graph.bond_in_ring(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) == bond.IsInRing()
        assert graph.bond_aromatic(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) == bond.GetIsAromatic()


def test_nonaromatic_ring_bonds_are_imported():
    mol = Chem.MolFromSmiles("O=C1CCCCC1")
    graph = graph_from_rdkit(mol)
    for bond in mol.GetBonds():
        assert graph.bond_in_ring(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) == bond.IsInRing()


def test_nonstandard_bonds_are_excluded_from_exact_bond_contracts():
    mol = Chem.MolFromSmiles("N->[Cu]")
    assert mol.GetBondWithIdx(0).GetBondType() == Chem.BondType.DATIVE
    for policy in ("strict", "fmcs"):
        with pytest.raises(ValueError, match="unsupported bond type"):
            graph_from_rdkit(mol, policy)
    graph = graph_from_rdkit(mol, "any")
    assert graph.n == 2 and graph.has_bond(0, 1)


def test_validation_rejects_missing_query_edge_disconnected_and_reused_atoms():
    ring = Chem.MolFromSmiles("C1CC1")
    path = Chem.MolFromSmiles("CCC")
    assert not validate_mapping(ring, path, {0:0, 1:1, 2:2}, "any")[0]
    assert not validate_mapping(path, path, {0:0, 2:2}, "any")[0]
    assert not validate_mapping(path, path, {0:0, 1:0}, "any")[0]
    assert validate_mapping(path, ring, {0:0, 1:1, 2:2}, "any")[0]


def test_flexible_validation_includes_aromatic_endpoint_resonance():
    left = Chem.MolFromSmiles("c1ccccc1-c1ccccc1")
    right = Chem.RWMol(left)
    right.GetBondBetweenAtoms(5, 6).SetBondType(Chem.BondType.DOUBLE)
    # Match graph semantics without re-perceiving this deliberately constructed
    # aromatic-endpoint bond; SMSD FLEXIBLE permits single/double resonance.
    mapping = {i: i for i in range(left.GetNumAtoms())}
    assert validate_mapping(left, right, mapping, "fmcs")[0]
    assert not validate_mapping(left, right, mapping, "strict")[0]


def test_objective_is_explicitly_atoms_and_timeouts_cannot_be_speed_wins():
    assert not rdkit_parameters("any", 10).MaximizeBonds
    rows = [{"engine": e, "elapsed_us": 10.0, "atoms": 3, "valid": True,
             "timeout": e == "rdkit", "budget_reached": False} for e in ("smsd", "rdkit")]
    result = summarize(rows, "any")
    assert result["smsd_timeouts"] is None
    assert not result["speed_comparable"] and result["rdkit_over_smsd_ratio"] is None
    rows[1]["timeout"] = False
    rows[0]["valid"] = False
    assert not summarize(rows, "any")["speed_comparable"]
    rows[0]["valid"] = True
    rows[0]["near_budget"] = True
    assert not summarize(rows, "any")["speed_comparable"]

# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Consistency tests for the public Python API.

These tests focus on wrapper parity and stable behavior across input types
instead of single-path implementation details.
"""

import pytest

import smsd
from smsd import (
    AromaticityModel,
    MolGraphBuilder,
    assign_rs,
    batch_mcs_size,
    compile_smarts,
    export_sdf,
    from_hex,
    fingerprint,
    find_mcs,
    parse_smiles,
    perceive_aromaticity,
    prewarm,
    read_molfile,
    smarts_match,
    smarts_find_all,
    find_substructure,
    to_hex,
    to_smiles,
    write_molfile,
)


try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - optional dependency
    Chem = None

try:
    from openbabel import pybel  # type: ignore[import]
except ImportError:  # pragma: no cover - optional dependency
    try:
        import pybel  # type: ignore[import]
    except ImportError:  # pragma: no cover - optional dependency
        pybel = None


MCS_CASES = [
    ("c1ccccc1", "c1ccc(O)cc1", 6),
    ("CC(=O)O", "CCC(=O)O", 4),
    ("C1CCCCC1", "C1CCCCC1", 6),
]


SUBSTRUCTURE_CASES = [
    ("c1ccccc1", "c1ccc(O)cc1", True),
    ("CC(=O)O", "CCC(=O)O", True),
    ("c1ccc(O)cc1", "c1ccccc1", False),
]


@pytest.mark.parametrize("query,target,expected", SUBSTRUCTURE_CASES)
def test_substructure_wrappers_stay_in_sync(query, target, expected):
    """SMILES and MolGraph wrappers should agree on substructure behavior."""
    gq = parse_smiles(query)
    gt = parse_smiles(target)

    assert smsd.is_substructure(gq, gt) is expected
    assert find_substructure(query, target) == find_substructure(gq, gt)

    mapping = find_substructure(gq, gt)
    assert bool(mapping) is expected
    if expected:
        assert len(mapping) == len(gq)


@pytest.mark.parametrize("query,target,min_size", MCS_CASES)
def test_mcs_wrappers_stay_in_sync(query, target, min_size):
    """SMILES and MolGraph wrappers should return the same MCS size."""
    gq = parse_smiles(query)
    gt = parse_smiles(target)

    mapping_smiles = smsd.find_mcs(query, target)
    mapping_graph = smsd.find_mcs(gq, gt)

    assert len(mapping_smiles) >= min_size
    assert len(mapping_graph) == len(mapping_smiles)

    if Chem is not None:
        rdkit_q = Chem.MolFromSmiles(query)
        rdkit_t = Chem.MolFromSmiles(target)
        native = smsd.mcs_rdkit_native(rdkit_q, rdkit_t)
        assert len(native) >= min_size
        for q_idx, t_idx in native.items():
            assert rdkit_q.GetAtomWithIdx(q_idx).GetAtomicNum() == \
                rdkit_t.GetAtomWithIdx(t_idx).GetAtomicNum()


@pytest.mark.parametrize("smiles", [
    "c1ccccc1",
    "CC(=O)O",
    "C1CCCCC1",
])
def test_lazy_canonical_fields_are_populated_by_smiles_rendering(smiles):
    """`to_smiles()` should populate the lazy ranking fields consistently."""
    mol = parse_smiles(smiles)
    assert mol.morgan_rank == []
    assert mol.canonical_label == []

    rendered = to_smiles(mol)
    assert isinstance(rendered, str)
    assert len(rendered) > 0
    assert len(mol.morgan_rank) == mol.n
    assert len(mol.canonical_label) == mol.n


def test_fingerprint_serialisation_round_trip():
    """Fingerprint serialisation helpers should round-trip cleanly."""
    fp = fingerprint("c1ccccc1")
    hex_fp = to_hex(fp)
    assert from_hex(hex_fp) == fp


def test_builder_api_constructs_expected_graph():
    """The builder should assemble a small graph using the public fluent API."""
    builder = MolGraphBuilder()
    builder.atom_count(3)
    builder.atomic_numbers([6, 6, 8])
    builder.formal_charges([0, 0, 0])
    builder.ring_flags([0, 0, 0])
    builder.aromatic_flags([0, 0, 0])
    builder.neighbors([[1], [0, 2], [1]])
    builder.bond_orders([[1], [1, 1], [1]])

    mol = builder.build()
    assert mol.n == 3
    assert mol.atomic_num == [6, 6, 8]
    assert mol.has_bond(0, 1) is True
    assert mol.has_bond(1, 2) is True


def test_gpu_helpers_return_stable_types():
    """GPU helpers should be safe to call on CPU-only machines."""
    assert isinstance(smsd.gpu_is_available(), bool)
    assert isinstance(smsd.gpu_device_info(), str)


@pytest.mark.skipif(pybel is None, reason="Open Babel bindings not installed")
def test_openbabel_interop_round_trip():
    """Optional Open Babel interop should round-trip through canonical SMILES."""
    obmol = pybel.readstring("smi", "c1ccc(O)cc1")
    g = smsd.from_openbabel(obmol)
    back = smsd.to_openbabel(g)

    assert to_smiles(g) == "Oc1ccccc1" or to_smiles(g) == "c1ccc(O)cc1"
    assert back.write("smi").strip()


@pytest.mark.parametrize("v3000", [False, True])
def test_native_molfile_helpers_round_trip(tmp_path, v3000):
    """Native file helpers should preserve chemistry and metadata."""
    mol = parse_smiles("N[C@@H](C)C(=O)O")
    mol.name = "alanine"
    mol.comment = "native molfile helper"
    mol.properties = {"ID": "AA-001"}

    path = tmp_path / ("alanine_v3000.mol" if v3000 else "alanine.mol")
    write_molfile(mol, path, v3000=v3000)
    reread = read_molfile(path)

    assert reread.n == mol.n
    assert reread.name == "alanine"
    assert reread.comment == "native molfile helper"
    assert reread.properties["ID"] == "AA-001"
    assert assign_rs(reread) == assign_rs(mol)


def test_native_sdf_export_has_no_rdkit_dependency(tmp_path):
    """export_sdf() should use the native writer for ordinary MolGraph inputs."""
    out = tmp_path / "library.sdf"
    mols = [parse_smiles("CCO"), parse_smiles("c1ccccc1"), parse_smiles("[R1]c1ccccc1")]
    for idx, mol in enumerate(mols, start=1):
        mol.name = f"mol-{idx}"
        mol.properties = {"ID": str(idx)}

    export_sdf(mols, out)
    text = out.read_text(encoding="utf-8")

    assert text.count("$$$$") == 3
    assert "> <ID>" in text
    assert "mol-1" in text


def test_rgroup_fragment_smarts_mcs_and_fingerprint_pipeline():
    """R-group fragments should stay usable across native search and fingerprint APIs."""
    core = parse_smiles("[R1]c1ccccc1")
    target = parse_smiles("[R2]c1ccc(O)cc1")

    assert to_smiles(core).startswith("[R1]")
    assert smarts_match("[*]c1ccccc1", target) is True
    assert len(fingerprint(core, kind="path")) > 0
    assert len(find_mcs(core, target)) >= 6


def test_compiled_smarts_queries_can_be_reused_across_targets():
    """Compiled SMARTS objects should avoid reparsing while matching repeatedly."""
    compiled = compile_smarts("[#6]~[#7]")
    targets = [parse_smiles("CCN"), parse_smiles("CCC"), parse_smiles("NCCN")]

    assert compiled.matches(targets[0]) is True
    assert compiled.matches(targets[1]) is False
    assert compiled.matches_many(targets) == [True, False, True]
    assert smarts_match(compiled, targets[0]) is True
    assert smarts_find_all(compiled, "CCN") == [{0: 1, 1: 2}]


def test_prewarm_helper_populates_lazy_graph_fields():
    """Explicit prewarm should fill canonical/orbit caches ahead of batch use."""
    mol = parse_smiles("c1ccccc1")
    assert mol.morgan_rank == []
    assert mol.canonical_label == []
    assert mol.orbit == []

    prepared = prewarm(mol)
    assert prepared is mol
    assert len(mol.morgan_rank) == mol.n
    assert len(mol.canonical_label) == mol.n
    assert len(mol.orbit) == mol.n


def test_explicit_aromaticity_perception_normalizes_raw_kekule_builder_graph():
    """Python should expose aromaticity perception as an explicit MolGraph operation."""
    builder = MolGraphBuilder()
    builder.atom_count(6)
    builder.atomic_numbers([6, 6, 6, 6, 6, 6])
    builder.formal_charges([0, 0, 0, 0, 0, 0])
    builder.ring_flags([0, 0, 0, 0, 0, 0])
    builder.aromatic_flags([0, 0, 0, 0, 0, 0])
    builder.neighbors([[1, 5], [0, 2], [1, 3], [2, 4], [3, 5], [4, 0]])
    builder.bond_orders([[2, 1], [2, 1], [1, 2], [2, 1], [1, 2], [2, 1]])

    mol = builder.build(
        perceive_aromaticity=False,
        aromaticity_model=AromaticityModel.DAYLIGHT_LIKE,
    )
    assert sum(bool(x) for x in mol.aromatic) == 0
    assert sum(bool(x) for x in mol.ring) == 0

    perceived = perceive_aromaticity(mol, model=AromaticityModel.DAYLIGHT_LIKE)
    assert perceived is mol
    assert sum(bool(x) for x in mol.aromatic) == 6
    assert sum(bool(x) for x in mol.ring) == 6
    assert to_smiles(mol) == "c1ccccc1"


def test_size_first_batch_helpers_keep_hit_order_and_sizes():
    """Size-oriented batch helpers should expose screening-friendly results."""
    query = parse_smiles("c1ccccc1")
    targets = [parse_smiles("c1ccccc1"), parse_smiles("c1ccc(O)cc1")]

    assert batch_mcs_size(query, targets, timeout_ms=1000) == [6, 6]

    mapping_hits = smsd.screen_and_match(query, targets, 0.95, timeout_ms=1000)
    size_hits = smsd.screen_and_mcs_size(query, targets, 0.95, timeout_ms=1000)

    assert len(mapping_hits) == 1
    assert mapping_hits[0][0] == 0
    assert len(mapping_hits[0][1]) == 6
    assert size_hits == [(0, 6)]


def test_stereo_smarts_substructure_and_mcs_pipeline():
    """Stereo-bearing molecules should remain stable across native matching APIs."""
    query = parse_smiles("N[C@@H](C)C(=O)O")
    target = parse_smiles("CC(N)C(=O)O")

    assert smarts_match("[NX3][CH]([#6])C(=O)O", query) is True
    assert find_substructure("NC(C)C(=O)O", target)
    assert len(find_mcs(query, target)) >= 5
    assert assign_rs(query)

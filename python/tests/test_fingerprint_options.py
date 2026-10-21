# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Configurable fingerprints: classic configurations keep the published bits, switches extend them."""

from pathlib import Path

import pytest

import smsd
from smsd import FingerprintFamily, FingerprintKey, FingerprintOptions, FragmentMode, ScaffoldMode

DATA_DIR = Path(__file__).resolve().parents[2] / "benchmarks/data"
MOLECULES = [
    "C", "CCO", "CC(=O)O", "c1ccccc1", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "C[C@H](N)C(=O)O", "C/C=C/C(=O)O", "[13CH4]", "[Na+].CC(=O)[O-]",
    "c1ccccc1.CCO", "C12C3C4C1C5C2C3C45", "O=C1CCCCC1", "NS(=O)(=O)c1ccc(Cl)cc1", "[H]OC", "C1CC1.CC.[He]",
]


def corpus(limit=3000):
    path = DATA_DIR / "chembl_mcs_benchmark.smi"
    if not path.is_file():
        return MOLECULES
    with path.open(encoding="utf-8") as handle:
        return MOLECULES + [line.split()[0] for line in handle if line.strip()][:limit]


def test_classic_configurations_keep_the_published_bits():
    for smiles in corpus():
        graph = smsd.parse_smiles(smiles)
        for family in ("ecfp", "fcfp"):
            for radius in (0, 2, -1):
                options = FingerprintOptions(family=family, radius=radius, size=1024)
                assert smsd.fingerprint(graph, options=options) == smsd.circular_fingerprint(graph, radius, 1024, family)
                options.counts = True
                assert smsd.fingerprint(graph, options=options) == smsd.circular_fingerprint_counts(graph, radius, 1024, family)
        assert smsd.fingerprint(graph, options=FingerprintOptions("path")) == smsd.path_fingerprint(graph)
        assert smsd.fingerprint(graph, options=FingerprintOptions.path(3, 512)) == smsd.path_fingerprint(graph, 3, 512)
        assert smsd.fingerprint(graph, options=FingerprintOptions("torsion")) == smsd.topological_torsion(graph)
        assert smsd.fingerprint(graph, options=FingerprintOptions(family="torsion", counts=True)) == smsd.topological_torsion_counts(graph)
    assert smsd.fingerprint("CCO", options=FingerprintOptions()) == smsd.circular_fingerprint(smsd.parse_smiles("CCO"))


def test_options_accept_names_and_members():
    assert FingerprintOptions(family=FingerprintFamily.PATH).family == FingerprintFamily.PATH
    assert FingerprintOptions("Torsion").family == FingerprintFamily.TORSION
    options = FingerprintOptions(fragments="largest", scaffold=ScaffoldMode.MURCKO)
    assert options.fragments == FragmentMode.LARGEST and options.scaffold == ScaffoldMode.MURCKO
    options.fragments = "per_fragment"
    assert options.fragments == FragmentMode.PER_FRAGMENT
    assert FingerprintOptions.ecfp().radius == 2 and FingerprintOptions.fcfp(3, 512).size == 512
    assert FingerprintOptions.path().path_length == 7 and FingerprintOptions.torsion(256).size == 256
    assert "path_length=5" in repr(FingerprintOptions("path", path_length=5))
    copy = options.copy()
    copy.radius = 3
    assert options.radius == 2


def test_classic_switches_follow_the_family():
    path = FingerprintOptions("path")
    assert path.element and not path.degree and not path.bond_order and not path.ring_membership
    fcfp = FingerprintOptions("fcfp")
    assert fcfp.degree and fcfp.ring_membership and not fcfp.element and not fcfp.hydrogen_count
    torsion = FingerprintOptions("torsion")
    assert torsion.element and torsion.degree and torsion.valence and torsion.ring_membership and not torsion.aromaticity
    ecfp = FingerprintOptions()
    assert all([ecfp.element, ecfp.degree, ecfp.valence, ecfp.formal_charge, ecfp.isotope, ecfp.hydrogen_count,
                ecfp.tautomer_class, ecfp.aromaticity, ecfp.bond_order, ecfp.ring_membership])
    assert not any([ecfp.aromatic_bond, ecfp.ring_bond, ecfp.ring_size, ecfp.cip_atoms, ecfp.cip_bonds, ecfp.counts])


def test_options_reject_bad_arguments():
    with pytest.raises(TypeError):
        FingerprintOptions(radii=3)
    with pytest.raises(TypeError):
        FingerprintOptions(radius=True)
    with pytest.raises(TypeError):
        FingerprintOptions(cip_atoms="yes")
    with pytest.raises(ValueError):
        FingerprintOptions(family="morgan")
    with pytest.raises(ValueError):
        smsd.fingerprint("CCO", options=FingerprintOptions(size=0))
    with pytest.raises(ValueError):
        smsd.fingerprint("CCO", options=FingerprintOptions(radius=-2))
    with pytest.raises(ValueError):
        smsd.fingerprint("CCO", options=FingerprintOptions("keys"))
    with pytest.raises(ValueError):
        smsd.fingerprint("CCO", options=FingerprintOptions("keys", keys=[("C", 0)]))
    with pytest.raises(TypeError):
        FingerprintOptions("keys", keys="c1ccccc1")


def test_stereo_terms_separate_isomers_only_where_a_centre_lies():
    plain, chiral = FingerprintOptions(size=4096), FingerprintOptions(size=4096, cip_atoms=True)
    left, right = "C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"
    assert smsd.fingerprint(left, options=plain) == smsd.fingerprint(right, options=plain)
    assert smsd.fingerprint(left, options=chiral) != smsd.fingerprint(right, options=chiral)
    assert smsd.fingerprint(left, options=FingerprintOptions(radius=0, cip_atoms=True)) == smsd.fingerprint(
        right, options=FingerprintOptions(radius=0, cip_atoms=True))
    differing = set(smsd.fingerprint(left, options=chiral)) ^ set(smsd.fingerprint(right, options=chiral))
    assert len(differing) <= 10
    geometric = FingerprintOptions(cip_bonds=True)
    assert smsd.fingerprint("C/C=C/C(=O)O", options=geometric) != smsd.fingerprint("C/C=C\\C(=O)O", options=geometric)
    for family in ("path", "torsion"):
        options = FingerprintOptions(family, cip_atoms=True)
        assert smsd.fingerprint(left, options=options) != smsd.fingerprint(right, options=options)


def test_ring_size_fragments_and_scaffold():
    ring_size = FingerprintOptions(ring_size=True)
    assert smsd.fingerprint("C1CCCCCC1", options=ring_size) != smsd.fingerprint("C1CCCCCCC1", options=ring_size)
    assert smsd.fingerprint("C1CCCCCCC1", options=ring_size) == smsd.fingerprint("C1CCCCCCCC1", options=ring_size)
    assert smsd.fingerprint("C1CCCCCCC1", options=FingerprintOptions(ring_size=True, ring_size_cap=12)) != smsd.fingerprint(
        "C1CCCCCCCC1", options=FingerprintOptions(ring_size=True, ring_size_cap=12))
    assert smsd.fingerprint("CCO.CC(=O)O", options=FingerprintOptions(fragments="largest")) == smsd.fingerprint(
        "CC(=O)O", options=FingerprintOptions())
    assert smsd.fingerprint("c1ccccc1.CCO", options=FingerprintOptions(fragments="per_fragment")) != smsd.fingerprint(
        "c1ccccc1.CCO", options=FingerprintOptions())
    benzene = smsd.fingerprint("c1ccccc1", options=FingerprintOptions())
    framework = FingerprintOptions(scaffold="murcko")
    assert smsd.fingerprint("Cc1ccccc1", options=framework) == benzene
    assert smsd.fingerprint("CC(=O)Nc1ccc(O)cc1", options=framework) == benzene
    assert smsd.fingerprint("CC(=O)c1ccccc1", options=framework) == benzene
    assert smsd.fingerprint("O=C1CCCCC1", options=framework) == smsd.fingerprint("O=C1CCCCC1", options=FingerprintOptions())
    assert smsd.fingerprint("Cn1cnc2c1c(=O)n(C)c(=O)n2C", options=framework) == smsd.fingerprint(
        "[nH]1cnc2c1c(=O)[nH]c(=O)[nH]2", options=FingerprintOptions())
    assert smsd.fingerprint("CCO", options=framework) == smsd.fingerprint("CCO", options=FingerprintOptions())
    chiral, framed = FingerprintOptions(cip_atoms=True), FingerprintOptions(cip_atoms=True, scaffold="murcko")
    assert smsd.fingerprint("C[C@]12CCCC[C@@H]1CCCO2", options=framed) == smsd.fingerprint("[C@H]12CCCC[C@@H]1CCCO2", options=chiral)
    assert smsd.fingerprint("C[C@@]12CCCC[C@@H]1CCCO2", options=framed) == smsd.fingerprint("[C@@H]12CCCC[C@@H]1CCCO2", options=chiral)
    geometric, framed_geometric = FingerprintOptions(cip_bonds=True), FingerprintOptions(cip_bonds=True, scaffold="murcko")
    assert smsd.fingerprint("C/C(c1ccccc1)=C/c1ccccc1", options=framed_geometric) == smsd.fingerprint(
        "c1ccccc1/C=C\\c1ccccc1", options=geometric)
    assert smsd.fingerprint("C1C/C=C\\CC1", options=geometric) == smsd.fingerprint("C1CC=CCC1", options=geometric)
    assert smsd.fingerprint("C1CCC/C=C/CC1", options=geometric) != smsd.fingerprint("C1CCC/C=C\\CC1", options=geometric)


def test_substructure_keys():
    maccs = FingerprintOptions.maccs()
    assert maccs.family == FingerprintFamily.KEYS and len(maccs.keys) == 166 and len(FingerprintOptions.maccs_keys()) == 166
    kinds = [key.kind for key in maccs.keys]
    assert kinds[0] == "isotope" and kinds[124] == "aromatic_rings" and kinds[165] == "fragments"
    assert all(kind == "smarts" for index, kind in enumerate(kinds) if index not in (0, 124, 165))
    assert maccs.keys[162].smarts == "*1~*~*~*~*~*~1" and maccs.keys[139].minimum == 4
    benzene = smsd.fingerprint("c1ccccc1", options=maccs)
    assert {162, 164, 161} <= set(benzene) and 163 not in benzene
    assert 165 in smsd.fingerprint("CC(=O)O.O", options=maccs)
    assert 0 in smsd.fingerprint("[13CH4]", options=maccs)
    assert 124 in smsd.fingerprint("c1ccc2ccccc2c1", options=maccs)
    assert smsd.fingerprint(smsd.MolGraph(), options=maccs) == []
    custom = FingerprintOptions("keys", keys=["c1ccccc1", ("[#8]", 2), FingerprintKey("[#8]", 5)])
    assert custom.size == 3 and [key.minimum for key in custom.keys] == [1, 2, 5]
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    assert smsd.fingerprint(aspirin, options=custom) == [0, 1]
    custom.counts = True
    assert smsd.fingerprint(aspirin, options=custom) == [(0, 1), (1, 4), (2, 4)]
    assert smsd.fingerprint_from_smiles(aspirin, options=maccs) == smsd.fingerprint(aspirin, options=maccs)


@pytest.mark.parametrize("first, second, switches", [
    ("C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", {"cip_atoms": True}),
    ("C/C=C/C(=O)O", "OC(=O)/C=C/C", {"cip_bonds": True}),
    ("c1ccccc1.CCO", "OCC.c1ccccc1", {"fragments": "per_fragment"}),
    ("CCO.CCN", "CCN.CCO", {"fragments": "largest"}),
    ("Cc1ccccc1CC(C)c1ccccc1O", "Oc1ccccc1C(C)Cc1ccccc1C", {"scaffold": "murcko"}),
    ("C1CC1c1ccccc1", "c1ccc(cc1)C1CC1", {"ring_size": True}),
])
def test_bits_do_not_depend_on_the_atom_order(first, second, switches):
    for family in ("ecfp", "path", "torsion"):
        options = FingerprintOptions(family, **switches)
        assert smsd.fingerprint(first, options=options) == smsd.fingerprint(second, options=options)


def test_options_take_precedence_over_the_classic_arguments():
    options = FingerprintOptions.torsion(512)
    assert smsd.fingerprint("CCCC", kind="path", fp_size=64, options=options) == smsd.topological_torsion("CCCC", 512)
    assert smsd.fingerprint_from_smiles("CCCC", 3, 64, "fcfp", options=options) == smsd.topological_torsion("CCCC", 512)
    with pytest.raises(ValueError):
        smsd.fingerprint_from_smiles("", options=options)

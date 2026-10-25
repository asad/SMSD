# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""fingerprints() gives the fingerprints of the single-molecule calls, in input order, on any thread count."""

from pathlib import Path

import pytest

import smsd
from smsd import FingerprintOptions

DATA_DIR = Path(__file__).resolve().parents[2] / "benchmarks/data"
MOLECULES = [
    "C", "CCO", "CC(=O)O", "c1ccccc1", "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1",
    "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "C[C@H](N)C(=O)O", "C/C=C/C(=O)O", "[13CH4]", "[Na+].CC(=O)[O-]",
    "c1ccccc1.CCO", "C12C3C4C1C5C2C3C45", "O=C1CCCCC1", "NS(=O)(=O)c1ccc(Cl)cc1", "[H]OC", "C1CC1.CC.[He]",
]


def library(limit=300):
    path = DATA_DIR / "chembl_mcs_benchmark.smi"
    if not path.is_file():
        return MOLECULES * 20
    with path.open(encoding="utf-8") as handle:
        return MOLECULES + [line.split()[0] for line in handle if line.strip()][:limit]


OPTIONS = {
    "ecfp4": FingerprintOptions.ecfp(2, 2048),
    "ecfp6": FingerprintOptions.ecfp(3, 2048),
    "fcfp4": FingerprintOptions.fcfp(2, 2048),
    "path": FingerprintOptions.path(7, 2048),
    "torsion": FingerprintOptions.torsion(2048),
    "maccs": FingerprintOptions.maccs(),
    "ecfp-cip-counts": FingerprintOptions("ecfp", radius=2, size=1024, cip_atoms=True, ring_size=True, counts=True),
    "torsion-counts": FingerprintOptions.torsion(512),
}
OPTIONS["torsion-counts"].counts = True


@pytest.mark.parametrize("name", list(OPTIONS))
def test_fingerprints_match_the_single_molecule_function(name):
    options = OPTIONS[name]
    smiles = library()
    graphs = [smsd.parse_smiles(text) for text in smiles]
    mixed = [graph if index % 3 else text for index, (graph, text) in enumerate(zip(graphs, smiles))]
    for molecules in (smiles, graphs, mixed):
        expected = [smsd.fingerprint(molecule, options=options) for molecule in molecules]
        for threads in (1, 2, 0):
            assert smsd.fingerprints(molecules, options=options, num_threads=threads) == expected


@pytest.mark.parametrize("arguments", [
    {}, {"kind": "ecfp"}, {"kind": "circular", "radius": 3}, {"kind": "fcfp", "fp_size": 1024},
    {"kind": "path", "path_length": 5}, {"kind": "mcs", "path_length": 3, "fp_size": 512}, {"kind": "torsion"},
])
def test_fingerprints_match_each_kind(arguments):
    smiles = library(100)
    expected = [smsd.fingerprint(text, **arguments) for text in smiles]
    for threads in (1, 2, 0):
        assert smsd.fingerprints(smiles, num_threads=threads, **arguments) == expected


def test_fingerprints_read_any_iterable_and_repeated_items():
    graph = smsd.parse_smiles("c1ccccc1O")
    expected = smsd.fingerprint(graph, kind="ecfp")
    assert smsd.fingerprints([graph] * 1000, kind="ecfp") == [expected] * 1000
    assert smsd.fingerprints((text for text in ["c1ccccc1O", "CCO"]), kind="ecfp") == [
        expected, smsd.fingerprint("CCO", kind="ecfp")]
    assert smsd.fingerprints([]) == []


def test_fingerprints_convert_rdkit_molecules():
    chem = pytest.importorskip("rdkit.Chem")
    molecules = [chem.MolFromSmiles(text) for text in ("CC(=O)Oc1ccccc1C(=O)O", "c1ccncc1")] * 300
    maccs = FingerprintOptions.maccs()
    assert smsd.fingerprints(molecules, options=maccs) == [smsd.fingerprint(m, options=maccs) for m in molecules]


def first_error(call):
    try:
        call()
    except Exception as error:
        return type(error), str(error)
    return None


@pytest.mark.parametrize("molecules", [
    ["CCO", "C1CC", 5],
    ["CCO", 5, "C1CC"],
    ["CCO", "", "C1CC"],
    ["CCO"] * 600 + ["C1CC"] + [None] * 600,
    ["CCO"] * 600 + [None] + ["C1CC"] * 600,
    ["CCO"] * 600 + ["C1CC", "C(C"] * 300,
])
def test_fingerprints_raise_the_error_of_the_first_failing_molecule(molecules):
    expected = first_error(lambda: [smsd.fingerprint(molecule, kind="ecfp") for molecule in molecules])
    assert expected is not None
    for threads in (1, 2, 0):
        assert first_error(lambda: smsd.fingerprints(molecules, kind="ecfp", num_threads=threads)) == expected


def test_fingerprints_check_their_arguments():
    for threads in (-1, 1.5, True):
        with pytest.raises(ValueError):
            smsd.fingerprints(["CCO"], num_threads=threads)
    with pytest.raises(ValueError, match="Unknown fingerprint kind"):
        smsd.fingerprints(["CCO"], kind="morgan")
    with pytest.raises(ValueError):
        smsd.fingerprints(["CCO"], fp_size=0)
    with pytest.raises(TypeError, match="FingerprintOptions"):
        smsd.fingerprints(["CCO"], options="ecfp")
    with pytest.raises(ValueError):
        smsd.fingerprints(["CCO"], options=FingerprintOptions("ecfp", radius=-2))

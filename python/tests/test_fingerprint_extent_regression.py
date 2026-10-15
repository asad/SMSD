"""Whole-graph fingerprint extent, aromatic torsions and size boundaries."""
import pytest
import smsd


@pytest.mark.parametrize("smiles,extent", [("CCO", 2), ("c1ccccc1", 3),
    ("CC(C)(C)C", 2), ("C1CC1.CC.[He]", 1), ("[He].[Ne]", 0)])
@pytest.mark.parametrize("mode", ["ecfp", "fcfp"])
@pytest.mark.parametrize("size", [1, 8, 17, 1000, 2048])
def test_whole_graph_radius_stops_at_the_largest_component_extent(smiles, extent, mode, size):
    graph = smsd.parse_smiles(smiles)
    assert smsd.circular_fingerprint(graph, -1, size, mode) == smsd.circular_fingerprint(graph, extent, size, mode)
    assert smsd.circular_fingerprint_counts(graph, -1, size, mode) == smsd.circular_fingerprint_counts(graph, extent, size, mode)


def test_benzene_whole_graph_counts_have_three_expansion_rounds():
    graph = smsd.parse_smiles("c1ccccc1")
    assert sum(count for _, count in smsd.circular_fingerprint_counts(graph, -1)) == 24
    assert smsd.circular_fingerprint(graph, -1) != smsd.circular_fingerprint(graph, 4)


def test_aromatic_and_saturated_torsions_are_distinct():
    aromatic = smsd.parse_smiles("c1ccccc1")
    saturated = smsd.parse_smiles("C1CCCCC1")
    assert smsd.topological_torsion(aromatic) != smsd.topological_torsion(saturated)
    assert smsd.topological_torsion_counts(aromatic) != smsd.topological_torsion_counts(saturated)
    assert smsd.topological_torsion(aromatic) == smsd.topological_torsion("C1=CC=CC=C1")


@pytest.mark.parametrize("size", [1, 8, 17, 1000])
def test_quality_handles_small_and_nonuniform_sizes(size):
    quality = smsd.analyze_fp_quality([size - 1, 0, size - 1], size)
    assert quality["set_bits"] == len({0, size - 1})
    assert quality["density"] == len({0, size - 1}) / size
    assert 0 <= quality["expected_collision_rate"] <= 1


@pytest.mark.parametrize("bits,size", [([], 0), ([], -1), ([-1], 8), ([8], 8)])
def test_quality_rejects_invalid_sizes_and_bit_positions(bits, size):
    with pytest.raises(ValueError):
        smsd.analyze_fp_quality(bits, size)


@pytest.mark.parametrize("length", [0, -1, 2**31 - 1])
@pytest.mark.parametrize("function", [smsd.path_fingerprint, smsd.mcs_fingerprint])
def test_path_fingerprints_reject_unrepresentable_lengths(function, length):
    with pytest.raises(ValueError):
        function(smsd.parse_smiles("CCO"), path_length=length)


def test_fixed_radius_hashes_are_preserved():
    graph = smsd.parse_smiles("CCO")
    assert smsd.circular_fingerprint(graph, 0) == [452, 1828, 1861]
    assert smsd.circular_fingerprint(graph, 2) == [273, 452, 501, 1237, 1670, 1702, 1746, 1828, 1861]
    assert smsd.circular_fingerprint(graph, 2, mode="fcfp") == [357, 395, 814, 1172, 1716, 1838, 1933, 1942, 2027]
    assert smsd.topological_torsion("C1CCCCC1") == [1177]

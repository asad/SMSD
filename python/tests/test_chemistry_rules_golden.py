# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
"""Shared chemistry-rule fixture, also read by the C++ and Java tests."""

import json
import re
from pathlib import Path

import pytest

import smsd
from smsd import _smsd


RULES_PATH = (
    Path(__file__).resolve().parents[2]
    / "java/src/test/resources/com/bioinception/smsd/chemistry_rules.json"
)
with RULES_PATH.open(encoding="utf-8") as handle:
    RULES = json.load(handle)


def cases(section):
    return [pytest.param(case, id=case.get("name", case.get("note", case.get("query"))))
            for case in RULES[section]["cases"]]


@pytest.mark.parametrize("case", cases("similarity"))
def test_similarity_metrics_share_their_empty_fingerprint_rules(case):
    first, second, size = case["first"], case["second"], case["fp_size"]
    counts = ([(bit, 1) for bit in first], [(bit, 1) for bit in second])
    for metric, function in (("tanimoto", smsd.tanimoto_coefficient), ("overlap", smsd.overlap_coefficient)):
        assert function(first, second, size) == pytest.approx(case[metric], abs=1e-12)
        assert function(*counts) == pytest.approx(case[metric], abs=1e-12)
    for metric in ("dice", "cosine", "soergel"):
        assert getattr(smsd, metric)(first, second, size) == pytest.approx(case[metric], abs=1e-12)


@pytest.mark.parametrize("case", cases("fingerprint_quality"))
def test_fingerprint_quality_shares_its_definitions(case):
    tolerance = RULES["fingerprint_quality"]["tolerance"]
    size, positions = case["fp_size"], list(range(case["set_bits"]))
    # Repeated and unordered positions give the same quality.
    for bits in (positions, positions[::-1] + positions):
        quality = smsd.analyze_fp_quality(bits, size)
        assert list(quality) == ["set_bits", "fp_size", "density", "expected_collision_rate", "is_saturated"]
        assert quality["set_bits"] == case["set_bits"] and quality["fp_size"] == size
        assert quality["density"] == pytest.approx(case["density"], abs=tolerance)
        assert quality["expected_collision_rate"] == pytest.approx(case["expected_collision_rate"], abs=tolerance)
        assert quality["is_saturated"] is case["is_saturated"]


@pytest.mark.parametrize("case", cases("circular_similarity"))
def test_batch_similarity_of_empty_molecules(case):
    query = smsd.parse_smiles(case["query"])
    targets = [smsd.parse_smiles(target) for target in case["targets"]]
    scores = dict(smsd.batch_similarity(query, targets))
    assert [scores[index] for index in range(len(targets))] == pytest.approx(case["tanimoto"], abs=1e-12)


@pytest.mark.parametrize("case", cases("implicit_hydrogens"))
def test_implicit_hydrogens_follow_the_shared_valence_table(case):
    assert smsd.implicit_h(case["atomic_number"], case["bond_order_sum"], case["charge"],
                           aromatic=case["aromatic"]) == case["hydrogens"]


def matched_atoms(smarts, graph):
    result = smsd.Substructure.smarts(smarts).match_all(graph)
    return sorted(match.atom_mapping[0] for match in result.matches)


def without_hydrogen_counts(graph):
    """The same graph built without hydrogen counts, so every count follows the valence model."""
    orders = [[graph.bond_order(atom, neighbor) for neighbor in graph.neighbors[atom]] for atom in range(graph.n)]
    return (smsd.MolGraphBuilder().atom_count(graph.n).atomic_numbers(list(graph.atomic_num))
            .formal_charges(list(graph.formal_charge)).ring_flags(list(graph.ring))
            .aromatic_flags(list(graph.aromatic)).neighbors(list(graph.neighbors)).bond_orders(orders).build())


@pytest.mark.parametrize("case", cases("smarts"))
def test_smarts_hydrogen_connectivity_and_valence(case):
    assert matched_atoms(case["smarts"], smsd.parse_smiles(case["smiles"])) == case["atoms"]


@pytest.mark.parametrize("case", cases("smarts_without_hydrogen_counts"))
def test_smarts_on_graphs_without_hydrogen_counts(case):
    graph = without_hydrogen_counts(smsd.parse_smiles(case["smiles"]))
    assert matched_atoms(case["smarts"], graph) == case["atoms"]


def tautomer_confidence_inputs(case):
    query, target = smsd.parse_smiles(case["query"]), smsd.parse_smiles(case["target"])
    atoms = {first: second for first, second in case["atoms"]}
    chemistry = smsd.ChemOptions.tautomer_profile()
    chemistry.solvent = getattr(smsd.Solvent, case["solvent"])
    chemistry.pH = case["pH"]
    bonds = tuple((atom, neighbor, image, atoms[neighbor])
                  for atom, image in sorted(atoms.items()) if 0 <= atom < query.n
                  for neighbor in query.neighbors[atom]
                  if neighbor > atom and neighbor in atoms and target.has_bond(image, atoms[neighbor]))
    return query, target, smsd.Match(atom_mapping=atoms, bond_mapping=bonds, fragments=()), chemistry


@pytest.mark.parametrize("case", cases("tautomer_confidence"))
def test_tautomer_confidence(case):
    query, target, match, chemistry = tautomer_confidence_inputs(case)
    if case.get("error"):
        with pytest.raises(IndexError):
            smsd.MCS.compute_tautomer_confidence(query, target, match, chemistry)
        return
    tolerance = RULES["tautomer_confidence"]["tolerance"]
    confidence = smsd.MCS.compute_tautomer_confidence(query, target, match, chemistry)
    assert confidence == pytest.approx(case["confidence"], abs=tolerance)
    if "same_as" in case:
        # The named case differs from this one only in the order of the atoms.
        (reference,) = [other for other in RULES["tautomer_confidence"]["cases"]
                        if other["name"] == case["same_as"]]
        assert confidence == pytest.approx(
            smsd.MCS.compute_tautomer_confidence(*tautomer_confidence_inputs(reference)), abs=tolerance)


@pytest.mark.parametrize("case", cases("path_fingerprints"))
def test_path_fingerprint_bits_are_unchanged(case):
    graph = smsd.parse_smiles(case["smiles"])
    length, size = case["path_length"], case["path_fp_size"]
    assert smsd.path_fingerprint(graph, length, size) == case["path_fingerprint"]
    words = _smsd.batch_fingerprint([graph], length, case["batch_fp_size"], 1)[0]
    assert [64 * index + bit for index, word in enumerate(words) for bit in range(64)
            if word >> bit & 1] == case["batch_fingerprint"]
    assert smsd.mcs_fingerprint(graph, length, size, case["tautomer_aware"]) == case["mcs_fingerprint"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"])
                                  for case in RULES["tetrahedral_centres"]["matches"]])
def test_chiral_substructure_follows_the_shared_ligand_rule(case):
    options = smsd.Substructure.Options()
    options.chemistry.use_chirality = True
    query, target = smsd.parse_smiles(case["query"]), smsd.parse_smiles(case["target"])
    assert smsd.Substructure.find(query).with_options(options).match(target).found is case["found"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["smiles"])
                                  for case in RULES["tetrahedral_centres"]["labels"]])
def test_cip_labels_follow_the_shared_ligand_rule(case):
    assert {str(atom): label for atom, label in smsd.assign_rs(case["smiles"]).items()} == case["labels"]


@pytest.mark.parametrize("case", cases("cip_mancude"))
def test_mancude_cip_labels_do_not_depend_on_the_kekule_structure(case):
    assert {str(atom): label for atom, label in smsd.assign_rs(case["smiles"]).items()} == case["labels"]


@pytest.mark.parametrize("case", cases("cip_auxiliary_descriptors"))
def test_auxiliary_descriptors_follow_sequence_rules_4_and_5(case):
    labels = smsd.assign_cip(case["smiles"]).rs_labels
    assert {str(atom): label.name for atom, label in enumerate(labels)
            if label != smsd.RSLabel.NONE} == case["labels"]


@pytest.mark.parametrize("case", cases("aromatic_input"))
def test_aromaticity_is_perceived_from_the_structure(case):
    graph = smsd.parse_smiles(case["smiles"])
    assert sum(bool(flag) for flag in graph.aromatic) == case["aromatic_atoms"]


@pytest.mark.parametrize("case", cases("cip_hierarchical_ranking"))
def test_first_pass_ranks_by_rule_1a_and_leaves_the_rest_to_the_hierarchical_rules(case):
    descriptors = smsd.assign_cip(case["smiles"])
    assert {str(atom): label.name for atom, label in enumerate(descriptors.rs_labels)
            if label != smsd.RSLabel.NONE} == case["labels"]
    assert {f"{min(first, second)}-{max(first, second)}": label.name
            for first, second, label in descriptors.ez_bonds} == case["bonds"]


def fingerprint_options_of(spec):
    """The options of a fixture case: the classic configuration of its family, then its keyword arguments."""
    spec = dict(spec)
    if spec.pop("maccs", False):
        options = smsd.FingerprintOptions.maccs()
        for key, value in spec.items():
            if key != "family":
                setattr(options, key, value)
        return options
    return smsd.FingerprintOptions(**spec)


@pytest.mark.parametrize("case", cases("fingerprint_options"))
def test_fingerprint_options_bits_are_shared_across_languages(case):
    graph = smsd.parse_smiles(case["smiles"]) if case["smiles"] else smsd.MolGraph()
    options = fingerprint_options_of(case["options"])
    assert smsd.fingerprint(graph, options=options) == case["bits"]
    options.counts = True
    assert [list(pair) for pair in smsd.fingerprint(graph, options=options)] == case["counts"]
    if case.get("same_bits_as"):
        other = next(item for item in RULES["fingerprint_options"]["cases"] if item["name"] == case["same_bits_as"])
        assert other["bits"] == case["bits"]
    if not case.get("classic"):
        return
    # The classic configuration of a family gives the bits of its fixed function.
    family, size = case["options"]["family"], case["options"]["size"]
    if family in ("ecfp", "fcfp"):
        radius = case["options"].get("radius", 2)
        assert smsd.circular_fingerprint(graph, radius, size, family) == case["bits"]
        assert [list(pair) for pair in smsd.circular_fingerprint_counts(graph, radius, size, family)] == case["counts"]
    elif family == "path":
        assert smsd.path_fingerprint(graph, case["options"].get("path_length", 7), size) == case["bits"]
    elif family == "torsion":
        assert smsd.topological_torsion(graph, size) == case["bits"]
        assert [list(pair) for pair in smsd.topological_torsion_counts(graph, size)] == case["counts"]


@pytest.mark.parametrize("case", cases("smarts_daylight"))
def test_smarts_follows_daylight_semantics(case):
    target = smsd.parse_smiles(case["smiles"])
    assert len(smsd.Substructure.smarts(case["smarts"]).match_all(target).matches) == case["embeddings"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"]) for case in RULES["smarts_daylight"]["rejected"]])
def test_smarts_rejects_conflicting_ring_closure_bonds(case):
    with pytest.raises(ValueError):
        smsd.Substructure.smarts(case["smarts"])


CHIRAL_SMARTS = re.compile(r"\[[^\]]*@")


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"]) for case in RULES["smarts_daylight"]["cases"]
                                  if CHIRAL_SMARTS.search(case["smarts"])])
def test_smarts_chirality_agrees_with_rdkit(case):
    """@ and @@ as RDKit matches them with chirality; @? and @@? as @ and @@ on a target with
    configured atoms, and as the pattern without chirality on a target without any."""
    chem = pytest.importorskip("rdkit.Chem")
    target = chem.MolFromSmiles(case["smiles"])
    query = chem.MolFromSmarts(case["smarts"].replace("?", ""))
    configured = any(atom.GetChiralTag() != chem.ChiralType.CHI_UNSPECIFIED for atom in target.GetAtoms())
    chirality = "?" not in case["smarts"] or configured
    expected = len(target.GetSubstructMatches(query, useChirality=chirality, uniquify=False))
    assert expected == case["embeddings"]


def profile_fields(options):
    """Every field of a profile as the fixture names it."""
    return {"match_atom_type": options.match_atom_type, "match_formal_charge": options.match_formal_charge,
            "use_chirality": options.use_chirality, "use_bond_stereo": options.use_bond_stereo,
            "ring_matches_ring_only": options.ring_matches_ring_only,
            "complete_rings_only": options.complete_rings_only, "match_isotope": options.match_isotope,
            "ring_fusion_mode": options.ring_fusion_mode.name, "tautomer_aware": options.tautomer_aware,
            "pH": options.pH, "solvent": options.solvent.name, "match_bond_order": options.match_bond_order.name,
            "aromaticity_mode": options.aromaticity_mode.name}


@pytest.mark.parametrize("case", cases("chemistry_profiles"))
def test_chemistry_profiles_change_only_their_listed_fields(case):
    expected = dict(RULES["chemistry_profiles"]["defaults"], **case["changes"])
    options = (smsd.ChemOptions.tautomer_profile() if case["profile"] == "tautomer_profile"
               else smsd.ChemOptions.profile(case["profile"]))
    assert profile_fields(options) == expected


def chemistry_of(fields):
    """The chemistry of a fixture: the defaults with the listed fields changed."""
    chemistry = smsd.ChemOptions()
    for name, value in fields.items():
        setattr(chemistry, name, getattr(smsd.AromaticityMode, value) if name == "aromaticity_mode" else value)
    return chemistry


@pytest.mark.parametrize("case", cases("similarity_upper_bound"))
def test_similarity_upper_bound_uses_the_shared_atom_labels(case):
    tolerance = RULES["similarity_upper_bound"]["tolerance"]
    chemistry = chemistry_of(case["chemistry"])
    query, target = smsd.parse_smiles(case["query"]), smsd.parse_smiles(case["target"])
    assert smsd.similarity_upper_bound(query, target, chemistry) == pytest.approx(case["bound"], abs=tolerance)
    assert smsd.similarity_upper_bound(target, query, chemistry) == pytest.approx(case["bound"], abs=tolerance)
    assert smsd.screen_targets(query, [target], case["bound"], chemistry) == [0]


@pytest.mark.parametrize("case", cases("consensus"))
def test_consensus_reduces_each_seed_in_input_order(case):
    options = smsd.MCS.Options()
    if "objective" in case:
        options.objective = getattr(smsd.MCS.Objective, case["objective"])
    result = smsd.MCS.consensus(case["molecules"], options, case["threshold"])
    assert result.found == case["found"] and result.complete
    if not result.found:
        return
    assert result.fragment.atom_count == case["atoms"] and result.fragment.bond_count == case["bonds"]
    assert result.support_count == case["support"] and result.source_index == case["source_index"]
    assert sorted(result.mappings) == case["mappings"] and result.graph.n == case["atoms"]
@pytest.mark.parametrize("case", [pytest.param(case, id=case["smiles"])
                                  for case in RULES["aromatic_bond_hydrogens"]["cases"]])
def test_explicit_aromatic_bonds_take_the_shared_hydrogen_counts(case):
    assert list(smsd.parse_smiles(case["smiles"]).hydrogen_count) == case["hydrogens"]


def stereo_chemistry(case, options):
    options.chemistry.use_chirality = case["chirality"]
    options.chemistry.use_bond_stereo = case["bond_stereo"]
    return options


@pytest.mark.parametrize("case", cases("stereo_embeddings"))
def test_substructure_embeddings_need_specified_stereo(case):
    options = stereo_chemistry(case, smsd.Substructure.Options())
    prepared = smsd.Substructure.find(case["query"]).with_options(options)
    assert prepared.match(case["target"]).found is case["found"]
    assert bool(prepared.match_all(case["target"]).matches) is case["found"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"])
                                  for case in RULES["stereo_embeddings"]["common_subgraphs"]])
def test_common_subgraphs_keep_unspecified_stereo_as_a_wildcard(case):
    options = stereo_chemistry(case, smsd.MCS.Options())
    result = smsd.MCS.find(case["query"]).with_options(options).match(case["target"])
    assert (result.best.atom_count if result.found else 0) == case["atoms"]


@pytest.mark.parametrize("case", cases("tautomer_consistency"))
def test_tautomer_consistency_counts_stored_hydrogens(case):
    mapping = {query_atom: target_atom for query_atom, target_atom in case["atoms"]}
    assert smsd.MCS.validate_tautomer_consistency(case["query"], case["target"], mapping) is case["consistent"]


@pytest.mark.parametrize("case", cases("tautomer_classes"))
def test_tautomer_classes_do_not_depend_on_the_atom_order(case):
    assert list(smsd.parse_smiles(case["smiles"]).tautomer_class) == case["classes"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["smiles"])
                                  for case in RULES["murcko_scaffolds"]["cases"]])
def test_murcko_scaffolds_cap_stripped_atoms_with_hydrogen(case):
    scaffold = smsd.murcko_scaffold(smsd.parse_smiles(case["smiles"]))
    assert list(scaffold.hydrogen_count) == case["hydrogens"]
    assert smsd.to_smiles(scaffold) == case["scaffold"]


@pytest.mark.parametrize("case", cases("consensus_native_order"))
def test_consensus_follows_the_native_reduction_order(case):
    options = smsd.MCS.Options()
    options.objective = getattr(smsd.MCS.Objective, case["objective"])
    options.connectivity = getattr(smsd.MCS.Connectivity, case["connectivity"].upper())
    result = smsd.MCS.consensus(case["molecules"], options, case["threshold"])
    assert result.found is case["found"]
    assert result.source_index == case["source_index"]
    atoms = sorted(result.fragment.atom_mapping) if result.found else []
    assert atoms == case["atoms"]
    assert (len(result.fragment.bond_mapping) if result.found else 0) == case["bonds"]
    assert sorted(result.mappings) == case["support"]
    assert result.support_count == len(case["support"])
    assert result.to_smiles() == case["smiles"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"])
                                  for case in RULES["strict_ring_fusion"]["ring_counts"]])
def test_ring_search_reports_its_bound(case):
    assert smsd.rings_truncated(case["smiles"]) is case["truncated"]


@pytest.mark.parametrize("case", [pytest.param(case, id=case["name"])
                                  for case in RULES["strict_ring_fusion"]["common_bonds"]])
def test_strict_ring_fusion_does_not_depend_on_atom_order(case):
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    options.chemistry.ring_fusion_mode = smsd.RingFusionMode.STRICT
    for query, target in ((case["query"], case["target"]), (case["target"], case["query"])):
        result = smsd.MCS.find(query).with_options(options).match(target)
        assert result.optimal
        assert result.best.bond_count == case["bonds"] and result.best.atom_count == case["atoms"]


@pytest.mark.parametrize("case", cases("double_bond_ligands"))
def test_double_bond_ligands_are_joined_to_their_end(case):
    options = smsd.MCS.Options()
    options.objective = smsd.MCS.Objective.MCES
    options.connectivity = getattr(smsd.MCS.Connectivity, case["connectivity"].upper())
    options.chemistry.use_bond_stereo = True
    result = smsd.MCS.find(case["query"]).with_options(options).match(case["target"])
    assert result.optimal
    assert result.best.bond_count == case["bonds"] and result.best.atom_count == case["atoms"]
    assert smsd.MCS.validate_mapping(case["query"], case["target"], result.best, options) == []

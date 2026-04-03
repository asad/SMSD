/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Consolidated parser test suite -- SMARTS matching, hydrogen handling,
 * and comprehensive SMILES parser tests (441 tests, gated by ifdef).
 *
 * Merges: test_smarts.cpp, test_hydrogen.cpp, test_smiles_comprehensive.cpp
 *
 * Compile:
 *   clang++ -std=c++17 -O2 -I include tests/test_parsers.cpp -o test_parsers
 *   clang++ -std=c++17 -O2 -I include -DSMSD_TEST_SMILES_COMPREHENSIVE \
 *           tests/test_parsers.cpp -o test_parsers_full
 */

#include "smsd/smsd.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/smiles_parser.hpp"

// Only include SMARTS parser when the guard is defined (always for this suite)
#include "smsd/smarts_parser.hpp"

#include <cassert>
#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <set>

// ============================================================================
// Common test harness
// ============================================================================

static int g_passed = 0;
static int g_failed = 0;

#define RUN_TEST(name) do { \
    std::cout << "  [" #name "] "; \
    try { test_##name(); g_passed++; std::cout << "PASS\n"; } \
    catch (const std::exception& e) { g_failed++; std::cout << "FAIL: " << e.what() << "\n"; } \
    catch (...) { g_failed++; std::cout << "FAIL (unknown)\n"; } \
} while(0)

// ============================================================================
// SECTION 1: SMARTS parser and matcher tests
// ============================================================================

// Convenience helpers
static bool smartsMatches(const std::string& smarts, const std::string& smiles) {
    auto q = smsd::parseSMARTS(smarts);
    auto mol = smsd::parseSMILES(smiles);
    return q.matches(mol);
}

static int smartsCountAll(const std::string& smarts, const std::string& smiles, int max = 1000) {
    auto q = smsd::parseSMARTS(smarts);
    auto mol = smsd::parseSMILES(smiles);
    return static_cast<int>(q.findAll(mol, max).size());
}

// --- Atom primitives ---

void test_smarts_atom_atomic_num_carbon() {
    assert(smartsMatches("[#6]", "C"));
    assert(smartsMatches("[#6]", "CC"));
    assert(!smartsMatches("[#6]", "[He]"));
}

void test_smarts_atom_atomic_num_nitrogen() {
    assert(smartsMatches("[#7]", "N"));
    assert(!smartsMatches("[#7]", "C"));
}

void test_smarts_atom_not_hydrogen() {
    assert(smartsMatches("[!#1]", "C"));
    assert(smartsMatches("[!#1]", "N"));
    assert(!smartsMatches("[!#1]", "[H]"));
}

void test_smarts_atom_degree_2() {
    assert(smartsMatches("[D2]", "CCC"));
    assert(!smartsMatches("[D2]", "C"));
}

void test_smarts_atom_degree_3() {
    assert(smartsMatches("[D3]", "CC(C)C"));
}

void test_smarts_atom_ring() {
    assert(smartsMatches("[R]", "C1CCC1"));
    assert(!smartsMatches("[R]", "CCCC"));
}

void test_smarts_atom_ring_size_6() {
    assert(smartsMatches("[r6]", "C1CCCCC1"));
    assert(!smartsMatches("[r6]", "C1CCC1"));
}

void test_smarts_atom_valence_4() {
    assert(smartsMatches("[v4]", "C(F)(F)(F)F"));
}

void test_smarts_atom_ring_connectivity_2() {
    assert(smartsMatches("[x2]", "C1CCC1"));
    assert(!smartsMatches("[x2]", "CCC"));
}

void test_smarts_atom_hcount_1() {
    assert(smartsMatches("[H1]", "C(Cl)(Cl)Cl"));
}

void test_smarts_atom_positive_charge() {
    assert(smartsMatches("[+1]", "[NH4+]"));
    assert(!smartsMatches("[+1]", "C"));
}

void test_smarts_atom_negative_charge() {
    assert(smartsMatches("[-1]", "[O-]"));
    assert(!smartsMatches("[-1]", "O"));
}

void test_smarts_atom_isotope() {
    assert(smartsMatches("[13C]", "[13CH4]"));
    assert(!smartsMatches("[13C]", "C"));
}

void test_smarts_atom_aromatic() {
    assert(smartsMatches("[a]", "c1ccccc1"));
    assert(!smartsMatches("[a]", "C1CCCCC1"));
}

void test_smarts_atom_aliphatic() {
    assert(smartsMatches("[A]", "C"));
    assert(smartsMatches("[A]", "CCC"));
}

void test_smarts_atom_wildcard_bracket() {
    assert(smartsMatches("[*]", "C"));
    assert(smartsMatches("[*]", "N"));
    assert(smartsMatches("[*]", "[He]"));
}

// --- Bond primitives ---

void test_smarts_bond_single() {
    assert(smartsMatches("C-C", "CC"));
    assert(!smartsMatches("C-C", "C=C"));
}

void test_smarts_bond_double() {
    assert(smartsMatches("C=C", "C=C"));
    assert(!smartsMatches("C=C", "CC"));
}

void test_smarts_bond_triple() {
    assert(smartsMatches("C#N", "C#N"));
    assert(!smartsMatches("C#N", "C=N"));
}

void test_smarts_bond_aromatic() {
    assert(smartsMatches("c:c", "c1ccccc1"));
}

void test_smarts_bond_any() {
    assert(smartsMatches("C~C", "CC"));
    assert(smartsMatches("C~C", "C=C"));
    assert(smartsMatches("C~C", "C#C"));
}

void test_smarts_bond_ring() {
    assert(smartsMatches("[#6]@[#6]", "C1CCC1"));
    assert(!smartsMatches("[#6]@[#6]", "CCC"));
}

void test_smarts_bond_wedge() {
    assert(smartsMatches("C/C", "CC"));
    assert(smartsMatches("C/C", "C=C"));
}

void test_smarts_bond_dash_stereo() {
    assert(smartsMatches("C\\C", "CC"));
    assert(smartsMatches("C\\C", "C=C"));
}

// --- Logical operators ---

void test_smarts_logic_or() {
    assert(smartsMatches("[#6,#7]", "C"));
    assert(smartsMatches("[#6,#7]", "N"));
    assert(!smartsMatches("[#6,#7]", "O"));
}

void test_smarts_logic_and_low() {
    assert(smartsMatches("[#6;R]", "C1CCC1"));
    assert(!smartsMatches("[#6;R]", "CCC"));
}

void test_smarts_logic_and_high() {
    assert(smartsMatches("[#6&R]", "C1CCC1"));
    assert(!smartsMatches("[#6&R]", "CCC"));
}

void test_smarts_logic_not() {
    assert(!smartsMatches("[!#6]", "C"));
    assert(smartsMatches("[!#6]", "N"));
}

void test_smarts_logic_complex_or_and() {
    assert(smartsMatches("[#6,#7;R]", "C1CCC1"));
    assert(!smartsMatches("[#6,#7;R]", "CCC"));
}

// --- Recursive SMARTS ---

void test_smarts_recursive_sp3_nitrogen() {
    assert(smartsMatches("[$([NX3])]", "CCN(C)C"));
}

void test_smarts_recursive_hydroxyl() {
    assert(smartsMatches("[$([OH])]", "CO"));
    assert(!smartsMatches("[$([OH])]", "C=O"));
}

void test_smarts_writer_roundtrip_ring_charge_and_hcount() {
    auto mol = smsd::parseSMILES("[NH3+]C1=CC=CC=C1");
    smsd::SmartsWriteOptions opts;
    opts.includeAromaticity = true;
    opts.includeCharge = true;
    opts.includeRingMember = true;
    opts.includeHCount = true;
    std::string smarts = smsd::toSMARTS(mol, opts);
    auto q = smsd::parseSMARTS(smarts);
    assert(q.matches(mol));
}

void test_smarts_writer_preserves_atom_class_queries() {
    auto mol = smsd::parseSMILES("[CH3:7][*:1]");
    std::string smarts = smsd::toSMARTS(mol);
    assert(smarts.find(":7") != std::string::npos);
    assert(smarts.find(":1") != std::string::npos);
    auto q = smsd::parseSMARTS(smarts);
    assert(q.matches(mol));
}

void test_mol_block_preserves_headers_and_properties() {
    std::string molBlock =
        "sample-name\n"
        "  SMSD  2D\n"
        "sample-comment\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n"
        "> <ID>\n"
        "cmpd-001\n"
        "\n"
        "> <NOTE>\n"
        "first line\n"
        "second line\n"
        "\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.name == "sample-name");
    assert(g.programLine == "  SMSD  2D");
    assert(g.comment == "sample-comment");
    assert(g.properties.at("ID") == "cmpd-001");
    assert(g.properties.at("NOTE") == "first line\nsecond line");

    auto written = smsd::writeMolBlock(g);
    assert(written.find("sample-name\n") == 0);
    assert(written.find("> <ID>\ncmpd-001\n\n") != std::string::npos);
    assert(written.find("> <NOTE>\nfirst line\nsecond line\n\n") != std::string::npos);
}

void test_v3000_roundtrip_metadata_and_atom_class() {
    auto mol = smsd::parseSMILES("[CH3:7][R1]");
    mol.name = "v3000-demo";
    mol.comment = "metadata";
    mol.properties["ID"] = "v30-1";
    std::string block = smsd::writeMolBlockV3000(mol);
    auto reparsed = smsd::readMolBlock(block);
    assert(reparsed.name == "v3000-demo");
    assert(reparsed.comment == "metadata");
    assert(reparsed.properties.at("ID") == "v30-1");
    assert(reparsed.n == mol.n);
    assert(reparsed.atomClass[0] == 7);
    assert(reparsed.atomClass[1] == 1);
}

void test_v3000_stereo_roundtrip() {
    auto mol = smsd::parseSMILES("N[C@@H](C)C(=O)O");
    std::string block = smsd::writeMolBlockV3000(mol);
    auto reparsed = smsd::readMolBlock(block);
    auto desc = smsd::cip::assignAll(reparsed);
    assert(!desc.rsLabels.empty());
    bool found = false;
    for (auto lab : desc.rsLabels) {
        if (lab != smsd::cip::RSLabel::NONE) { found = true; break; }
    }
    assert(found);
}

void test_v2000_rgroup_roundtrip() {
    auto mol = smsd::parseSMILES("[R1]c1ccccc1");
    std::string block = smsd::writeMolBlock(mol);
    assert(block.find("R#") != std::string::npos);
    assert(block.find("M  RGP") != std::string::npos);
    auto reparsed = smsd::readMolBlock(block);
    assert(reparsed.atomicNum[0] == 0);
    assert(reparsed.atomClass[0] == 1);
}

void test_v2000_alias_rgroup_reads() {
    std::string molBlock =
        "alias-rgroup\n"
        "  SMSD\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "A    1\n"
        "R2\n"
        "M  END\n";
    auto g = smsd::readMolBlock(molBlock);
    assert(g.atomicNum[0] == 0);
    assert(g.atomClass[0] == 2);
}

void test_v2000_multi_attachment_roundtrip() {
    auto mol = smsd::parseSMILES("[R1]c1ccc([R2])cc1");
    std::string block = smsd::writeMolBlock(mol);
    auto reparsed = smsd::readMolBlock(block);
    int labels = 0;
    for (int cls : reparsed.atomClass) {
        if (cls == 1 || cls == 2) labels++;
    }
    assert(labels == 2);
}

void test_v3000_double_bond_ez_roundtrip() {
    auto mol = smsd::parseSMILES("C/C=C/C");
    std::string block = smsd::writeMolBlockV3000(mol);
    auto reparsed = smsd::readMolBlock(block);
    auto ez = smsd::cip::assignEZ(reparsed, 1, 2);
    assert(ez == smsd::cip::EZLabel::E);
}

void test_v3000_multiple_stereocentres_roundtrip() {
    auto mol = smsd::parseSMILES("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1");
    std::string block = smsd::writeMolBlockV3000(mol);
    auto reparsed = smsd::readMolBlock(block);
    auto desc = smsd::cip::assignAll(reparsed);
    int count = 0;
    for (auto lab : desc.rsLabels) {
        if (lab != smsd::cip::RSLabel::NONE) count++;
    }
    assert(count >= 4);
}

void test_smarts_recursive_carbon_bonded_to_oxygen() {
    assert(smartsMatches("[$([#6][#8])]", "CCO"));
    assert(!smartsMatches("[$([#6][#8])]", "CCC"));
}

// --- Named predicate registry ---

void test_smarts_predicate_registry_basic() {
    auto reg = smsd::SmartsPredicateRegistry::builtinRegistry();
    std::string expanded = reg.expandPredicates("[$isHalogen]");
    assert(expanded == "[$([F,Cl,Br,I])]");
}

void test_smarts_predicate_registry_custom() {
    smsd::SmartsPredicateRegistry reg;
    reg.registerPredicate("myPred", "[#6]");
    std::string expanded = reg.expandPredicates("[$myPred]");
    assert(expanded == "[$([#6])]");
}

void test_smarts_predicate_no_clobber_recursive() {
    smsd::SmartsPredicateRegistry reg;
    reg.registerPredicate("foo", "[#7]");
    std::string input = "[$([#6])]";
    std::string expanded = reg.expandPredicates(input);
    assert(expanded == input);
}

void test_smarts_predicate_isKetone_expand() {
    auto reg = smsd::SmartsPredicateRegistry::builtinRegistry();
    std::string expanded = reg.expandPredicates("[$isKetone]");
    assert(expanded == "[$([#6][CX3](=O)[#6])]");
}

void test_smarts_predicate_isHalogen_match() {
    auto reg = smsd::SmartsPredicateRegistry::builtinRegistry();
    std::string expanded = reg.expandPredicates("[$isHalogen]");
    auto q = smsd::parseSMARTS(expanded);
    auto mol = smsd::parseSMILES("CF");
    assert(q.matches(mol));
}

void test_smarts_predicate_isPositive_match() {
    auto reg = smsd::SmartsPredicateRegistry::builtinRegistry();
    std::string expanded = reg.expandPredicates("[$isPositive]");
    auto q = smsd::parseSMARTS(expanded);
    assert(q.matches(smsd::parseSMILES("[NH4+]")));
    assert(!q.matches(smsd::parseSMILES("C")));
}

void test_smarts_predicate_isNegative_match() {
    auto reg = smsd::SmartsPredicateRegistry::builtinRegistry();
    std::string expanded = reg.expandPredicates("[$isNegative]");
    auto q = smsd::parseSMARTS(expanded);
    assert(q.matches(smsd::parseSMILES("[O-]C")));
}

// --- Complex patterns ---

void test_smarts_complex_carboxylic_acid() {
    assert(smartsMatches("[CX3](=O)[OX2H1]", "CC(=O)O"));
    assert(!smartsMatches("[CX3](=O)[OX2H1]", "CC"));
}

void test_smarts_complex_aromatic_ring() {
    assert(smartsMatches("c1ccccc1", "c1ccccc1"));
}

void test_smarts_complex_primary_amine_not_amide() {
    assert(smartsMatches("[NX3;H2,H1;!$(NC=O)]", "CCN"));
}

void test_smarts_complex_carbonyl() {
    assert(smartsMatches("[CX3]=[OX1]", "CC=O"));
    assert(smartsMatches("[CX3]=[OX1]", "CC(=O)C"));
}

// --- Edge cases ---

void test_smarts_edge_empty_pattern() {
    auto q = smsd::parseSMARTS("");
    auto mol = smsd::parseSMILES("C");
    assert(q.matches(mol));
}

void test_smarts_edge_single_bracket_atom() {
    assert(smartsMatches("[C]", "C"));
    assert(smartsMatches("[C]", "CC"));
}

void test_smarts_edge_wildcard_star() {
    assert(smartsMatches("*", "C"));
    assert(smartsMatches("*", "N"));
    assert(smartsMatches("*", "[He]"));
    assert(smartsMatches("*", "c1ccccc1"));
}

void test_smarts_edge_invalid_smarts_throws() {
    bool threw = false;
    try {
        smsd::parseSMARTS("[#6");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
}

void test_smarts_edge_multi_component() {
    auto q = smsd::parseSMARTS("C.N");
    auto mol = smsd::parseSMILES("CN");
    assert(q.matches(mol));
}

void test_smarts_edge_findAll_count() {
    int count = smartsCountAll("[#6]", "CCCC");
    assert(count == 4);
}

// --- Drug molecule matching ---

void test_smarts_drug_benzene_in_aspirin() {
    assert(smartsMatches("c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"));
}

void test_smarts_drug_amide_in_acetaminophen() {
    assert(smartsMatches("C(=O)N", "CC(=O)Nc1ccc(O)cc1"));
}

void test_smarts_drug_carboxyl_in_aspirin() {
    assert(smartsMatches("C(=O)O", "CC(=O)Oc1ccccc1C(=O)O"));
}

void test_smarts_drug_hydroxyl_in_acetaminophen() {
    assert(smartsMatches("[OH]", "CC(=O)Nc1ccc(O)cc1"));
}

void test_smarts_drug_ester_in_aspirin() {
    assert(smartsMatches("[CX3](=O)[OX2][#6]", "CC(=O)Oc1ccccc1C(=O)O"));
}

// --- RDKit-parity SMARTS extensions ---

void test_smarts_ext_hetero_neighbors_z2_pyrimidine() {
    assert(smartsMatches("[z2]", "c1ncncc1"));
}

void test_smarts_ext_hetero_neighbors_z0() {
    assert(smartsMatches("[z0]", "CC"));
    assert(!smartsMatches("[z1]", "CC"));
}

void test_smarts_ext_hetero_neighbors_z1_cn() {
    assert(smartsMatches("[z1]", "CN"));
}

void test_smarts_ext_hetero_neighbors_default() {
    assert(smartsMatches("[z]", "CN"));
    assert(!smartsMatches("[z]", "CC"));
}

void test_smarts_ext_aliphatic_hetero_Z1() {
    assert(smartsMatches("[Z1]", "CN"));
}

void test_smarts_ext_aliphatic_hetero_aromatic() {
    assert(!smartsMatches("[Z1]", "c1ccncc1"));
}

void test_smarts_ext_heavy_degree_d2() {
    assert(smartsMatches("[d2]", "CCC"));
    assert(!smartsMatches("[d2]", "C=C"));
}

void test_smarts_ext_heavy_degree_d1() {
    assert(smartsMatches("[d1]", "CC"));
}

void test_smarts_ext_heavy_degree_d3() {
    assert(smartsMatches("[d3]", "CC(C)C"));
}

void test_smarts_ext_hybridization_sp2_aromatic() {
    assert(smartsMatches("[^2]", "c1ccccc1"));
}

void test_smarts_ext_hybridization_sp3() {
    assert(smartsMatches("[^3]", "C"));
    assert(smartsMatches("[^3]", "CC"));
}

void test_smarts_ext_hybridization_sp2_double_bond() {
    assert(smartsMatches("[^2]", "C=C"));
}

void test_smarts_ext_hybridization_sp() {
    assert(smartsMatches("[^1]", "C#C"));
}

void test_smarts_ext_hybridization_sp3_not_sp2() {
    assert(!smartsMatches("[^2]", "C"));
}

void test_smarts_ext_range_degree_2_to_4() {
    assert(smartsMatches("[D{2-4}]", "CC(C)(C)C"));
    assert(!smartsMatches("[D{2-4}]", "C"));
}

void test_smarts_ext_range_min_only() {
    assert(smartsMatches("[D{2-}]", "CCC"));
    assert(smartsMatches("[D{2-}]", "CC(C)(C)C"));
    assert(!smartsMatches("[D{2-}]", "C"));
}

void test_smarts_ext_range_max_only() {
    assert(smartsMatches("[D{-2}]", "CCC"));
    assert(smartsMatches("[D{-2}]", "C"));
}

void test_smarts_ext_range_hcount() {
    assert(smartsMatches("[H{1-3}]", "CCC"));
    assert(!smartsMatches("[H{1-3}]", "C(Cl)(Cl)(Cl)Cl"));
}

void test_smarts_ext_range_exact_brace() {
    assert(smartsMatches("[D{3}]", "CC(C)C"));
    assert(!smartsMatches("[D{3}]", "CCC"));
}

void test_smarts_ext_combined_hybridization_and_element() {
    assert(smartsMatches("[#6^2]", "C=C"));
    assert(!smartsMatches("[#6^2]", "CC"));
}

void test_smarts_ext_range_valence() {
    assert(smartsMatches("[v{2-4}]", "C=C"));
    assert(!smartsMatches("[v{2-4}]", "C"));
}

void test_smarts_ext_range_hetero_neighbors() {
    assert(smartsMatches("[z{1-2}]", "CN"));
    assert(!smartsMatches("[z{1-2}]", "CC"));
}

// ============================================================================
// SECTION 2: Hydrogen handling tests
// ============================================================================

static smsd::McsOptions defaultMcsOpts() {
    smsd::McsOptions o;
    o.timeoutMs = 10000;
    return o;
}

static int mcsSize(const std::string& smi1, const std::string& smi2,
                   smsd::ChemOptions chem = smsd::ChemOptions{}) {
    auto g1 = smsd::parseSMILES(smi1);
    auto g2 = smsd::parseSMILES(smi2);
    auto m  = smsd::findMCS(g1, g2, chem, defaultMcsOpts());
    return static_cast<int>(m.size());
}

static bool isSub(const std::string& query, const std::string& target,
                  smsd::ChemOptions chem = smsd::ChemOptions{}) {
    auto q = smsd::parseSMILES(query);
    auto t = smsd::parseSMILES(target);
    return smsd::isSubstructure(q, t, chem, 10000);
}

// -- Implicit-H normalisation --

void test_h_methane_implicit_vs_bracket() {
    assert(isSub("[CH4]", "C"));
    assert(isSub("C", "[CH4]"));
    int sz = mcsSize("C", "[CH4]");
    std::cout << "mcs=" << sz << " ";
    assert(sz == 1);
}

void test_h_benzene_bracketed_aromatic_h() {
    int sz = mcsSize("c1ccccc1", "[cH]1[cH][cH][cH][cH][cH]1");
    std::cout << "mcs=" << sz << " ";
    assert(sz == 6);
}

void test_h_ethanol_bracket_vs_implicit() {
    int sz = mcsSize("CCO", "[CH3][CH2][OH]");
    std::cout << "mcs=" << sz << " ";
    assert(sz == 3);
}

void test_h_aspirin_bracket_notation() {
    int sz = mcsSize(
        "CC(=O)Oc1ccccc1C(=O)O",
        "[CH3]C(=O)Oc1ccccc1C(=O)O"
    );
    std::cout << "mcs=" << sz << " ";
    assert(sz == 13);
}

void test_h_water_OH2_vs_implicit() {
    assert(isSub("O",     "[OH2]"));
    assert(isSub("[OH2]", "O"));
    std::cout << "both_ways=ok ";
}

// -- Explicit [H] as graph nodes --

void test_h_explicit_H_methane_atom_count() {
    auto mol = smsd::parseSMILES("[H]C([H])([H])[H]");
    std::cout << "n=" << mol.n << " ";
    assert(mol.n == 5);
}

void test_h_explicit_H_query_no_match_implicit_target() {
    bool sub = isSub("[H]C([H])([H])[H]", "C");
    std::cout << "sub=" << (sub ? "true" : "false") << " ";
    assert(!sub);
}

void test_h_explicit_H_water_mcs_heavy_only() {
    int sz = mcsSize("[H]O[H]", "O");
    std::cout << "mcs=" << sz << " ";
    assert(sz == 1);
}

void test_h_H2_self_match() {
    assert(isSub("[H][H]", "[H][H]"));
    std::cout << "ok ";
}

// -- H-count notation and ring/aromaticity --

void test_h_cyclohexane_CH2_notation() {
    int sz = mcsSize("C1CCCCC1", "[CH2]1[CH2][CH2][CH2][CH2][CH2]1");
    std::cout << "mcs=" << sz << " ";
    assert(sz == 6);
}

void test_h_caffeine_theophylline_nH() {
    int sz = mcsSize(
        "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "Cn1cnc2c1c(=O)[nH]c(=O)n2C"
    );
    std::cout << "mcs=" << sz << " ";
    assert(sz >= 11);
}

// -- Drug-like molecules --

void test_h_morphine_stereo_h_notation() {
    int sz = mcsSize(
        "CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
        "C[N]1CC[C@@]23[C@@H]4[C@H]1C[C@@H]5=C([C@@H]2[C@H](C=C4)O3)C=C(C=C5)O"
    );
    std::cout << "mcs=" << sz << " ";
    assert(sz >= 17);
}

void test_h_ibuprofen_stereo_bracket() {
    int sz = mcsSize(
        "CC(C)Cc1ccc(CC(C)C(=O)O)cc1",
        "[CH3][C@@H]([CH3])Cc1ccc(CC(C)C(=O)O)cc1"
    );
    std::cout << "mcs=" << sz << " ";
    assert(sz >= 12);
}

// -- Formal charge orthogonality --

void test_h_ammonium_vs_ammonia_no_charge() {
    smsd::ChemOptions chem;
    chem.matchFormalCharge = false;
    int sz = mcsSize("N", "[NH4+]", chem);
    std::cout << "mcs=" << sz << " ";
    assert(sz == 1);
}

void test_h_carboxylate_vs_acid_no_charge() {
    smsd::ChemOptions chem;
    chem.matchFormalCharge = false;
    int sz = mcsSize("CC(=O)[O-]", "CC(=O)O", chem);
    std::cout << "mcs=" << sz << " ";
    assert(sz == 4);
}

// ============================================================================
// SECTION 3: Comprehensive SMILES parser tests (441 tests, gated by ifdef)
//
// When -DSMSD_TEST_SMILES_COMPREHENSIVE is set, the test_smiles_comprehensive.cpp
// translation unit is compiled alongside this file. It exports
// smiles_comprehensive_run(int& pass, int& fail).
// ============================================================================

#ifdef SMSD_TEST_SMILES_COMPREHENSIVE
extern void smiles_comprehensive_run(int& outPassed, int& outFailed);
#endif

// ============================================================================
// main
// ============================================================================

int main() {
    std::cout << "=== SMSD Parser Test Suite ===\n\n";

    // --- SMARTS tests ---
    std::cout << "-- SMARTS atom primitives --\n";
    RUN_TEST(smarts_atom_atomic_num_carbon);
    RUN_TEST(smarts_atom_atomic_num_nitrogen);
    RUN_TEST(smarts_atom_not_hydrogen);
    RUN_TEST(smarts_atom_degree_2);
    RUN_TEST(smarts_atom_degree_3);
    RUN_TEST(smarts_atom_ring);
    RUN_TEST(smarts_atom_ring_size_6);
    RUN_TEST(smarts_atom_valence_4);
    RUN_TEST(smarts_atom_ring_connectivity_2);
    RUN_TEST(smarts_atom_hcount_1);
    RUN_TEST(smarts_atom_positive_charge);
    RUN_TEST(smarts_atom_negative_charge);
    RUN_TEST(smarts_atom_isotope);
    RUN_TEST(smarts_atom_aromatic);
    RUN_TEST(smarts_atom_aliphatic);
    RUN_TEST(smarts_atom_wildcard_bracket);

    std::cout << "\n-- SMARTS bond primitives --\n";
    RUN_TEST(smarts_bond_single);
    RUN_TEST(smarts_bond_double);
    RUN_TEST(smarts_bond_triple);
    RUN_TEST(smarts_bond_aromatic);
    RUN_TEST(smarts_bond_any);
    RUN_TEST(smarts_bond_ring);
    RUN_TEST(smarts_bond_wedge);
    RUN_TEST(smarts_bond_dash_stereo);

    std::cout << "\n-- SMARTS logical operators --\n";
    RUN_TEST(smarts_logic_or);
    RUN_TEST(smarts_logic_and_low);
    RUN_TEST(smarts_logic_and_high);
    RUN_TEST(smarts_logic_not);
    RUN_TEST(smarts_logic_complex_or_and);

    std::cout << "\n-- SMARTS recursive --\n";
    RUN_TEST(smarts_recursive_sp3_nitrogen);
    RUN_TEST(smarts_recursive_hydroxyl);
    RUN_TEST(smarts_recursive_carbon_bonded_to_oxygen);
    RUN_TEST(smarts_writer_roundtrip_ring_charge_and_hcount);
    RUN_TEST(smarts_writer_preserves_atom_class_queries);

    std::cout << "\n-- SMARTS named predicates --\n";
    RUN_TEST(smarts_predicate_registry_basic);
    RUN_TEST(smarts_predicate_registry_custom);
    RUN_TEST(smarts_predicate_no_clobber_recursive);
    RUN_TEST(smarts_predicate_isKetone_expand);
    RUN_TEST(smarts_predicate_isHalogen_match);
    RUN_TEST(smarts_predicate_isPositive_match);
    RUN_TEST(smarts_predicate_isNegative_match);

    std::cout << "\n-- SMARTS complex patterns --\n";
    RUN_TEST(smarts_complex_carboxylic_acid);
    RUN_TEST(smarts_complex_aromatic_ring);
    RUN_TEST(smarts_complex_primary_amine_not_amide);
    RUN_TEST(smarts_complex_carbonyl);

    std::cout << "\n-- SMARTS edge cases --\n";
    RUN_TEST(smarts_edge_empty_pattern);
    RUN_TEST(smarts_edge_single_bracket_atom);
    RUN_TEST(smarts_edge_wildcard_star);
    RUN_TEST(smarts_edge_invalid_smarts_throws);
    RUN_TEST(smarts_edge_multi_component);
    RUN_TEST(smarts_edge_findAll_count);

    std::cout << "\n-- SMARTS drug matching --\n";
    RUN_TEST(smarts_drug_benzene_in_aspirin);
    RUN_TEST(smarts_drug_amide_in_acetaminophen);
    RUN_TEST(smarts_drug_carboxyl_in_aspirin);
    RUN_TEST(smarts_drug_hydroxyl_in_acetaminophen);
    RUN_TEST(smarts_drug_ester_in_aspirin);

    std::cout << "\n-- SMARTS RDKit-parity extensions --\n";
    RUN_TEST(smarts_ext_hetero_neighbors_z2_pyrimidine);
    RUN_TEST(smarts_ext_hetero_neighbors_z0);
    RUN_TEST(smarts_ext_hetero_neighbors_z1_cn);
    RUN_TEST(smarts_ext_hetero_neighbors_default);
    RUN_TEST(smarts_ext_aliphatic_hetero_Z1);
    RUN_TEST(smarts_ext_aliphatic_hetero_aromatic);
    RUN_TEST(smarts_ext_heavy_degree_d2);
    RUN_TEST(smarts_ext_heavy_degree_d1);
    RUN_TEST(smarts_ext_heavy_degree_d3);
    RUN_TEST(smarts_ext_hybridization_sp2_aromatic);
    RUN_TEST(smarts_ext_hybridization_sp3);
    RUN_TEST(smarts_ext_hybridization_sp2_double_bond);
    RUN_TEST(smarts_ext_hybridization_sp);
    RUN_TEST(smarts_ext_hybridization_sp3_not_sp2);
    RUN_TEST(smarts_ext_range_degree_2_to_4);
    RUN_TEST(smarts_ext_range_min_only);
    RUN_TEST(smarts_ext_range_max_only);
    RUN_TEST(smarts_ext_range_hcount);
    RUN_TEST(smarts_ext_range_exact_brace);
    RUN_TEST(smarts_ext_combined_hybridization_and_element);
    RUN_TEST(smarts_ext_range_valence);
    RUN_TEST(smarts_ext_range_hetero_neighbors);

    // --- Hydrogen handling tests ---
    std::cout << "\n-- Hydrogen: implicit-H normalisation --\n";
    RUN_TEST(h_methane_implicit_vs_bracket);
    RUN_TEST(h_benzene_bracketed_aromatic_h);
    RUN_TEST(h_ethanol_bracket_vs_implicit);
    RUN_TEST(h_aspirin_bracket_notation);
    RUN_TEST(h_water_OH2_vs_implicit);

    std::cout << "\n-- Hydrogen: explicit [H] as graph nodes --\n";
    RUN_TEST(h_explicit_H_methane_atom_count);
    RUN_TEST(h_explicit_H_query_no_match_implicit_target);
    RUN_TEST(h_explicit_H_water_mcs_heavy_only);
    RUN_TEST(h_H2_self_match);

    std::cout << "\n-- Hydrogen: ring/aromaticity unaffected by H notation --\n";
    RUN_TEST(h_cyclohexane_CH2_notation);
    RUN_TEST(h_caffeine_theophylline_nH);

    std::cout << "\n-- Hydrogen: drug-like molecules --\n";
    RUN_TEST(h_morphine_stereo_h_notation);
    RUN_TEST(h_ibuprofen_stereo_bracket);

    std::cout << "\n-- Hydrogen: formal charge orthogonality --\n";
    RUN_TEST(h_ammonium_vs_ammonia_no_charge);
    RUN_TEST(h_carboxylate_vs_acid_no_charge);

    std::cout << "\n-- MOL/SDF metadata --\n";
    RUN_TEST(mol_block_preserves_headers_and_properties);
    RUN_TEST(v2000_rgroup_roundtrip);
    RUN_TEST(v2000_alias_rgroup_reads);
    RUN_TEST(v2000_multi_attachment_roundtrip);
    RUN_TEST(v3000_roundtrip_metadata_and_atom_class);
    RUN_TEST(v3000_stereo_roundtrip);
    RUN_TEST(v3000_double_bond_ez_roundtrip);
    RUN_TEST(v3000_multiple_stereocentres_roundtrip);

#ifdef SMSD_TEST_SMILES_COMPREHENSIVE
    {
        std::cout << "\n-- Comprehensive SMILES parser tests (441) --\n";
        int cp = 0, cf = 0;
        smiles_comprehensive_run(cp, cf);
        g_passed += cp;
        g_failed += cf;
    }
#endif

    std::cout << "\n================================================\n";
    std::cout << g_passed << " passed, " << g_failed << " failed\n";
    return g_failed > 0 ? 1 : 0;
}

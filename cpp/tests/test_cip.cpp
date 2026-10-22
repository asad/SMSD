/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * CIP stereo descriptor test suite.
 * Tests R/S assignment for tetrahedral centres and E/Z for double bonds.
 */
#include "test_support.hpp"

#include "smsd/smsd.hpp"
#include "smsd/cip.hpp"
#include "smsd/smiles_parser.hpp"
#include <cassert>
#include <cstdio>
#include <iostream>
#include <string>

// ============================================================================
// Test harness
// ============================================================================

static int g_pass = 0, g_fail = 0;

#define RUN_TEST(name) do { \
    std::cout << "  [" #name "] "; \
    try { test_##name(); std::cout << "PASS\n"; } \
    catch (const std::exception& e) { g_fail++; std::cout << "FAIL: " << e.what() << "\n"; } \
    catch (...) { g_fail++; std::cout << "FAIL (unknown)\n"; } \
} while(0)

#define CIP_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "  FAIL: %s  [%s:%d]\n", msg, __FILE__, __LINE__); \
        throw std::runtime_error(msg); \
    } else { \
        g_pass++; \
    } \
} while(0)

using namespace smsd;
using namespace smsd::cip;

// ============================================================================
// R/S tests
// ============================================================================

void test_L_alanine() {
    // L-alanine: N[C@@H](C)C(=O)O -> S
    auto g = parseSMILES("N[C@@H](C)C(=O)O");
    RSLabel rs = assignRS(g, 1);
    CIP_ASSERT(rs == RSLabel::S, "L-alanine should be S");
}

void test_component_root_ring_stereo() {
    using smsd::cip::RSLabel;
    auto s = smsd::parseSMILES("CC1.[C@H]1(N)O");
    auto r = smsd::parseSMILES("CC1.[C@@H]1(N)O");
    CIP_ASSERT(smsd::cip::assignRS(s, 2) == RSLabel::S, "component-root ring closure should preserve S");
    CIP_ASSERT(smsd::cip::assignRS(r, 2) == RSLabel::R, "component-root ring closure should preserve R");
    auto reference = smsd::parseSMILES("CC[C@@H](N)O");
    smsd::ChemOptions chemistry;
    chemistry.useChirality = true;
    CIP_ASSERT(smsd_test::substructureExists(s, reference, chemistry), "equivalent stereo traversals should embed");
    CIP_ASSERT(!smsd_test::substructureExists(r, reference, chemistry), "opposite stereoisomers should not embed");
    CIP_ASSERT(smsd_test::commonAtoms(s, reference, chemistry, {}).size() == 5, "equivalent stereo MCS should include all atoms");
}

void test_D_alanine() {
    // D-alanine: N[C@H](C)C(=O)O -> R
    auto g = parseSMILES("N[C@H](C)C(=O)O");
    RSLabel rs = assignRS(g, 1);
    CIP_ASSERT(rs == RSLabel::R, "D-alanine should be R");
}

void test_bromochlorofluoromethane() {
    // Bromochlorofluoromethane: [C@@H](Br)(Cl)F -> S
    // Priorities: Br(35) > Cl(17) > F(9) > H(1)
    auto g = parseSMILES("[C@@H](Br)(Cl)F");
    RSLabel rs = assignRS(g, 0);
    CIP_ASSERT(rs == RSLabel::S, "Bromochlorofluoromethane @@ should be S");
}

void test_bromochlorofluoromethane_R() {
    // [C@H](Br)(Cl)F -> R
    auto g = parseSMILES("[C@H](Br)(Cl)F");
    RSLabel rs = assignRS(g, 0);
    CIP_ASSERT(rs == RSLabel::R, "Bromochlorofluoromethane @ should be R");
}

void test_L_cysteine() {
    // L-cysteine: N[C@@H](CS)C(=O)O -> R
    // Note: the thiol group (CS) has higher CIP priority than COOH
    // because S (Z=16) > O (Z=8) at the first point of difference
    auto g = parseSMILES("N[C@@H](CS)C(=O)O");
    RSLabel rs = assignRS(g, 1);
    CIP_ASSERT(rs == RSLabel::R, "L-cysteine should be R");
}

void test_R_glyceraldehyde() {
    // (R)-glyceraldehyde: OC[C@@H](O)C=O -> R
    auto g = parseSMILES("OC[C@@H](O)C=O");
    RSLabel rs = assignRS(g, 2);
    CIP_ASSERT(rs == RSLabel::R, "(R)-glyceraldehyde should be R");
}

void test_no_stereocentre_methane() {
    // CH4 has no stereocentre
    auto g = parseSMILES("C");
    RSLabel rs = assignRS(g, 0);
    CIP_ASSERT(rs == RSLabel::NONE, "Methane has no stereocentre");
}

void test_no_stereocentre_symmetric() {
    // CC(C)C - isobutane, central carbon has two identical methyls
    auto g = parseSMILES("CC(C)C");
    // No chirality annotation, so should return NONE
    RSLabel rs = assignRS(g, 1);
    CIP_ASSERT(rs == RSLabel::NONE, "Isobutane central carbon has no chirality annotation");
}

void test_citric_acid_no_stereocentre() {
    // Citric acid: OC(=O)CC(O)(CC(=O)O)C(=O)O
    // The central carbon has OH + COOH + CH2COOH + CH2COOH
    // Two identical CH2COOH groups -> not a stereocentre
    auto g = parseSMILES("OC(=O)CC(O)(CC(=O)O)C(=O)O");
    // No chirality annotation in the SMILES, should be NONE
    for (int i = 0; i < g.n; ++i) {
        RSLabel rs = assignRS(g, i);
        CIP_ASSERT(rs == RSLabel::NONE,
                   "Citric acid should have no stereocentres");
    }
}

// ============================================================================
// E/Z tests
// ============================================================================

void test_E_2_butene() {
    // (E)-2-butene: C/C=C/C -> E (trans)
    auto g = parseSMILES("C/C=C/C");
    // Double bond is between atoms 1 and 2
    EZLabel ez = assignEZ(g, 1, 2);
    CIP_ASSERT(ez == EZLabel::E, "(E)-2-butene should be E");
}

void test_Z_2_butene() {
    // (Z)-2-butene: C/C=C\C -> Z (cis)
    auto g = parseSMILES("C/C=C\\C");
    EZLabel ez = assignEZ(g, 1, 2);
    CIP_ASSERT(ez == EZLabel::Z, "(Z)-2-butene should be Z");
}

void test_no_EZ_ethene() {
    // Ethene: C=C has no E/Z (identical H substituents on each end)
    auto g = parseSMILES("C=C");
    EZLabel ez = assignEZ(g, 0, 1);
    // No / \ annotation in SMILES, so dbStereo is 0
    CIP_ASSERT(ez == EZLabel::NONE, "Ethene has no E/Z (no annotation)");
}

void test_no_EZ_single_bond() {
    // Single bond has no E/Z
    auto g = parseSMILES("CC");
    EZLabel ez = assignEZ(g, 0, 1);
    CIP_ASSERT(ez == EZLabel::NONE, "Single bond has no E/Z");
}

// ============================================================================
// assignAll / assignFromSMILES tests
// ============================================================================

void test_assignAll_L_alanine() {
    auto g = parseSMILES("N[C@@H](C)C(=O)O");
    auto desc = assignAll(g);
    CIP_ASSERT(desc.rsLabels.size() == static_cast<size_t>(g.n),
               "rsLabels size matches atom count");
    CIP_ASSERT(desc.rsLabels[1] == RSLabel::S,
               "assignAll: L-alanine C1 should be S");
    // No stereogenic double bonds in alanine
    CIP_ASSERT(desc.ezBonds.empty(),
               "assignAll: L-alanine has no E/Z bonds");
}

void test_assignFromSMILES() {
    auto desc = assignFromSMILES("N[C@@H](C)C(=O)O");
    CIP_ASSERT(desc.rsLabels[1] == RSLabel::S,
               "assignFromSMILES: L-alanine should be S");
}

void test_assignAll_E_2_butene() {
    auto g = parseSMILES("C/C=C/C");
    auto desc = assignAll(g);
    CIP_ASSERT(desc.ezBonds.size() == 1,
               "E-2-butene should have one E/Z bond");
    auto [a1, a2, ez] = desc.ezBonds[0];
    CIP_ASSERT(ez == EZLabel::E,
               "assignAll: E-2-butene should be E");
}

// ============================================================================
// Priority computation tests
// ============================================================================

void test_priority_basic() {
    // In bromochlorofluoromethane, priorities should be:
    // H(1) < F(9) < Cl(17) < Br(35)
    auto g = parseSMILES("[C@@H](Br)(Cl)F");
    auto prios = cip::detail::computePriorities(g, 0);
    // prios: {atomIdx, rank}
    // Find ranks for each ligand
    int rankBr = -1, rankCl = -1, rankF = -1, rankH = -1;
    for (auto& [idx, rank] : prios) {
        if (idx == 1) rankBr = rank; // Br
        if (idx == 2) rankCl = rank; // Cl
        if (idx == 3) rankF = rank;  // F
        if (idx == -1) rankH = rank; // implicit H
    }
    CIP_ASSERT(rankH < rankF, "H priority < F priority");
    CIP_ASSERT(rankF < rankCl, "F priority < Cl priority");
    CIP_ASSERT(rankCl < rankBr, "Cl priority < Br priority");
}

void test_digraph_basic() {
    // Simple test: build digraph for methane
    auto g = parseSMILES("C");
    auto dg = cip::detail::buildDigraph(g, 0);
    CIP_ASSERT(!dg.empty(), "Digraph should not be empty");
    CIP_ASSERT(dg[0].atomIdx == 0, "Root should be atom 0");
}

void test_to_string() {
    CIP_ASSERT(to_string(RSLabel::R) == "R", "R to_string");
    CIP_ASSERT(to_string(RSLabel::S) == "S", "S to_string");
    CIP_ASSERT(to_string(RSLabel::NONE) == "", "NONE to_string");
    CIP_ASSERT(to_string(EZLabel::E) == "E", "E to_string");
    CIP_ASSERT(to_string(EZLabel::Z) == "Z", "Z to_string");
    CIP_ASSERT(to_string(EZLabel::NONE) == "", "NONE EZ to_string");
}

// Sequence-rule cases shared with the Java and Python suites: duplicated
// substituents, isotopes, mancude rings, sulfoxides and double-bond ligands.
void test_cip_priority_sequence() {
    struct Example { const char* smiles; int centre; RSLabel expected; };
    const Example examples[] = {
        {"F[C@H](C(OC)C)C(O)CO", 1, RSLabel::S},
        {"[2H]C([2H])([C@H](O)CCF)C", 3, RSLabel::S},
        {"O[C@H](C[13CH2]C)CCCF", 1, RSLabel::R},
        {"F[C@H](C(CO)C(CBr)(CBr)CBr)C(C(O)C)C(C)C", 1, RSLabel::R},
        {"O[C@@H](c1ccccc1)C(C)(C)C", 1, RSLabel::R},
        {"O[C@@H](C1=CC=CC=C1)C(C)(C)C", 1, RSLabel::R},
        {"O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1", 10, RSLabel::S},
        {"C[S@](=O)c1ccccc1", 1, RSLabel::S},
        {"COc1ccc2[nH]c([S@@](=O)Cc3ncc(C)c(OC)c3C)nc2c1", 8, RSLabel::S},
        {"C/C=C/[C@@H](F)/C=C\\C", 3, RSLabel::R},
        {"C/C([C@H](F)/C(C)=C/C)=C/C", 2, RSLabel::R},
        {"C/C=C(C)/[C@@H](F)/C(C)=C\\C", 4, RSLabel::R},
        {"C/C(Cl)=C(C)/[C@@H](F)/C(C)=C(Cl)\\C", 5, RSLabel::S},
        {"N[C@@H](C)C(=O)O", 1, RSLabel::S},
        {"N[C@H](C)C(=O)O", 1, RSLabel::R}
    };
    for (const auto& example : examples) {
        const auto graph = parseSMILES(example.smiles);
        const std::string message = std::string("sequence rules for ") + example.smiles
            + " expected " + to_string(example.expected);
        CIP_ASSERT(assignRS(graph, example.centre) == example.expected, message.c_str());
        CIP_ASSERT(assignAll(graph).rsLabels[example.centre] == example.expected, message.c_str());
    }
}

// ============================================================================
// Rules 4 and 5: auxiliary descriptors and pseudoasymmetric centres
// ============================================================================

void test_pseudoasymmetric_1_4_disubstituted_cyclohexane() {
    // Both ring carbons of a 1,4-disubstituted cyclohexane are told apart only by
    // Rule 5, so they are pseudoasymmetric; assignRS() reads every rule, as assignAll() does.
    auto g = parseSMILES("C[C@H]1CC[C@H](O)CC1");
    CIP_ASSERT(assignRS(g, 1) == RSLabel::r, "assignRS gives C1 of 4-methylcyclohexanol its Rule 5 descriptor");
    auto all = assignAll(g);
    CIP_ASSERT(all.rsLabels[1] == RSLabel::r, "(1r,4r)-4-methylcyclohexanol C1 should be r");
    CIP_ASSERT(all.rsLabels[4] == RSLabel::r, "(1r,4r)-4-methylcyclohexanol C4 should be r");
    auto other = assignAll(parseSMILES("C[C@H]1CC[C@@H](O)CC1"));
    CIP_ASSERT(other.rsLabels[1] == RSLabel::s && other.rsLabels[4] == RSLabel::s,
               "the other diastereomer of 4-methylcyclohexanol should be (1s,4s)");
}

void test_pseudoasymmetric_acyclic_centre() {
    // 2,4-dichloropentan-3-ol: C3 carries an R and an S ligand of the same constitution.
    auto s = assignAll(parseSMILES("C[C@H](Cl)[C@H](O)[C@@H](C)Cl"));
    CIP_ASSERT(s.rsLabels[1] == RSLabel::S && s.rsLabels[5] == RSLabel::R, "the flanking centres keep R and S");
    CIP_ASSERT(s.rsLabels[3] == RSLabel::s, "2,4-dichloropentan-3-ol C3 should be s");
    auto r = assignAll(parseSMILES("C[C@H](Cl)[C@@H](O)[C@@H](C)Cl"));
    CIP_ASSERT(r.rsLabels[3] == RSLabel::r, "the epimer at C3 should be r");
}

void test_rule_4b_like_pairs_give_a_true_centre() {
    // Nateglinide: the ring carbon bearing the isopropyl group is resolved by Rule 4b
    // (like against unlike descriptor pairs), so it is a true centre, while the carbon
    // bearing the amide is pseudoasymmetric.
    auto all = assignAll(parseSMILES("CC(C)[C@H]1CC[C@H](C(=O)N[C@H](Cc2ccccc2)C(=O)[O-])CC1"));
    CIP_ASSERT(all.rsLabels[3] == RSLabel::R, "nateglinide C4 (isopropyl) should be R by Rule 4b");
    CIP_ASSERT(all.rsLabels[6] == RSLabel::r, "nateglinide C1 (amide) should be r");
    CIP_ASSERT(all.rsLabels[10] == RSLabel::R, "nateglinide phenylalanine centre should be R");
}

void test_tropane_c3_reads_the_bridgehead_descriptors() {
    auto ester = assignAll(parseSMILES("C[NH+]1[C@@H]2CC[C@H]1C[C@@H](OC(=O)C(CO)c1ccccc1)C2"));
    CIP_ASSERT(ester.rsLabels[2] == RSLabel::R && ester.rsLabels[5] == RSLabel::S, "tropane bridgeheads are R and S");
    CIP_ASSERT(ester.rsLabels[7] == RSLabel::r, "tropane C3 should be r");
    auto epimer = assignAll(parseSMILES("CN1[C@@H]2CC[C@H]1C[C@H](OC(=O)C(CO)c1ccccc1)C2"));
    CIP_ASSERT(epimer.rsLabels[7] == RSLabel::s, "the C3 epimer should be s");
}

void test_hexachlorocyclohexane_mixes_true_and_pseudoasymmetric_centres() {
    auto all = assignAll(parseSMILES("Cl[C@H]1[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H]1Cl"));
    const RSLabel expected[] = {RSLabel::S, RSLabel::r, RSLabel::R, RSLabel::S, RSLabel::r, RSLabel::R};
    const int atoms[] = {1, 2, 4, 6, 8, 10};
    for (int i = 0; i < 6; ++i)
        CIP_ASSERT(all.rsLabels[atoms[i]] == expected[i], "hexachlorocyclohexane ring carbon label");
}

void test_symmetric_centre_stays_unlabelled() {
    // Two identical ligands: no rule tells them apart, and the tag is left unlabelled.
    auto all = assignAll(parseSMILES("C[C@H](C)O"));
    CIP_ASSERT(all.rsLabels[1] == RSLabel::NONE, "propan-2-ol carries no descriptor");
    CIP_ASSERT(assignRS(parseSMILES("C[C@H](C)O"), 1, {{0, 1}}) == RSLabel::NONE,
               "the full rules leave propan-2-ol unlabelled as well");
}

void test_rule_2_unlabelled_atom_counts_its_standard_atomic_weight() {
    // Sequence Rule 2 reads an unlabelled atom at its standard atomic weight, so an isotope
    // lighter than that weight ranks below the unlabelled element: 18F (18) below F (18.998).
    CIP_ASSERT(assignRS(parseSMILES("C[C@H](F)[C@H](Cl)[C@H](C)[18F]"), 3) == RSLabel::R,
               "the 18F branch ranks below the F branch at C3");
    // Both passes read the rule the same way, so the pseudoasymmetric centre between an
    // 18F-bearing and an F-bearing ligand agrees with the descriptors of those ligands.
    auto all = assignAll(parseSMILES("[18F][C@H](F)[C@H](O)[C@H](F)[18F]"));
    CIP_ASSERT(all.rsLabels[1] == RSLabel::R && all.rsLabels[5] == RSLabel::S, "the flanking centres are R and S");
    CIP_ASSERT(all.rsLabels[3] == RSLabel::r, "the centre between them is r");
}

void test_assign_rs_applies_every_rule() {
    // The first pass ranks by Rule 1a alone; a ring-closure duplicate (Rule 1b) outranks the 13C
    // label (Rule 2) further out, so the bicyclic centre is R, as the hierarchical rules read it.
    CIP_ASSERT(assignRS(parseSMILES("C1=[13CH]C[C@@H]2C[C@H]12"), 3) == RSLabel::R,
               "Rule 1b decides before the isotope");
    // assignRS() gives the descriptor assignAll() does, pseudoasymmetric centres included.
    const auto graph = parseSMILES("Cl[C@H]1[C@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@H](Cl)[C@H]1Cl");
    const auto all = assignAll(graph);
    for (int atom = 0; atom < graph.n; ++atom)
        CIP_ASSERT(assignRS(graph, atom) == all.rsLabels[atom], "assignRS agrees with assignAll");
}

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << "=== CIP Stereo Descriptor Tests ===\n";

    std::cout << "\n--- R/S assignment ---\n";
    RUN_TEST(L_alanine);
    RUN_TEST(component_root_ring_stereo);
    RUN_TEST(cip_priority_sequence);
    RUN_TEST(D_alanine);
    RUN_TEST(bromochlorofluoromethane);
    RUN_TEST(bromochlorofluoromethane_R);
    RUN_TEST(L_cysteine);
    RUN_TEST(R_glyceraldehyde);
    RUN_TEST(no_stereocentre_methane);
    RUN_TEST(no_stereocentre_symmetric);
    RUN_TEST(citric_acid_no_stereocentre);

    std::cout << "\n--- Rules 4 and 5 ---\n";
    RUN_TEST(pseudoasymmetric_1_4_disubstituted_cyclohexane);
    RUN_TEST(pseudoasymmetric_acyclic_centre);
    RUN_TEST(rule_4b_like_pairs_give_a_true_centre);
    RUN_TEST(tropane_c3_reads_the_bridgehead_descriptors);
    RUN_TEST(hexachlorocyclohexane_mixes_true_and_pseudoasymmetric_centres);
    RUN_TEST(symmetric_centre_stays_unlabelled);
    RUN_TEST(rule_2_unlabelled_atom_counts_its_standard_atomic_weight);
    RUN_TEST(assign_rs_applies_every_rule);

    std::cout << "\n--- E/Z assignment ---\n";
    RUN_TEST(E_2_butene);
    RUN_TEST(Z_2_butene);
    RUN_TEST(no_EZ_ethene);
    RUN_TEST(no_EZ_single_bond);

    std::cout << "\n--- Convenience API ---\n";
    RUN_TEST(assignAll_L_alanine);
    RUN_TEST(assignFromSMILES);
    RUN_TEST(assignAll_E_2_butene);

    std::cout << "\n--- Internal helpers ---\n";
    RUN_TEST(priority_basic);
    RUN_TEST(digraph_basic);
    RUN_TEST(to_string);

    std::cout << "\n=== Results: " << g_pass << " passed, "
              << g_fail << " failed ===\n";
    return g_fail > 0 ? 1 : 0;
}

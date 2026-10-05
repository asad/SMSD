/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * CIP stereo descriptor test suite.
 * Tests R/S assignment for tetrahedral centres and E/Z for double bonds.
 */

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

// ============================================================================
// Main
// ============================================================================

int main() {
    std::cout << "=== CIP Stereo Descriptor Tests ===\n";

    std::cout << "\n--- R/S assignment ---\n";
    RUN_TEST(L_alanine);
    RUN_TEST(D_alanine);
    RUN_TEST(bromochlorofluoromethane);
    RUN_TEST(bromochlorofluoromethane_R);
    RUN_TEST(L_cysteine);
    RUN_TEST(R_glyceraldehyde);
    RUN_TEST(no_stereocentre_methane);
    RUN_TEST(no_stereocentre_symmetric);
    RUN_TEST(citric_acid_no_stereocentre);

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

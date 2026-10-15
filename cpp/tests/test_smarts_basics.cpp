/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/substructure.hpp"
#include <cassert>
#include <iostream>

namespace smsd_smarts_test {

inline void testWildcard() {
    auto q = smsd::parseSMARTS("*");
    auto mol = smsd::parseSMILES("C");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testWildcard\n";
}

inline void testElementMatch() {
    auto q = smsd::parseSMARTS("[#6]");
    auto mol = smsd::parseSMILES("C");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto q2 = smsd::parseSMARTS("[#7]");
    assert(smsd::Substructure::smarts(q2).match(mol).matches.empty());
    std::cout << "  [PASS] testElementMatch\n";
}

inline void testOrganicShorthand() {
    auto q = smsd::parseSMARTS("C");
    auto mol = smsd::parseSMILES("CC");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto q2 = smsd::parseSMARTS("N");
    assert(smsd::Substructure::smarts(q2).match(mol).matches.empty());
    std::cout << "  [PASS] testOrganicShorthand\n";
}

inline void testBondSingle() {
    auto q = smsd::parseSMARTS("C-C");
    auto mol = smsd::parseSMILES("CC");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("C=C");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testBondSingle\n";
}

inline void testBondDouble() {
    auto q = smsd::parseSMARTS("C=O");
    auto mol = smsd::parseSMILES("CC=O");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testBondDouble\n";
}

inline void testBondTriple() {
    auto q = smsd::parseSMARTS("C#N");
    auto mol = smsd::parseSMILES("CC#N");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("CC=N");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testBondTriple\n";
}

inline void testBondAny() {
    auto q = smsd::parseSMARTS("C~C");
    auto mol1 = smsd::parseSMILES("CC");
    auto mol2 = smsd::parseSMILES("C=C");
    assert(!smsd::Substructure::smarts(q).match(mol1).matches.empty());
    assert(!smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testBondAny\n";
}

inline void testDegree() {
    // D2 = degree 2
    auto q = smsd::parseSMARTS("[D2]");
    auto mol = smsd::parseSMILES("CCC");  // middle C has degree 2
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("C");   // single C has degree 0
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testDegree\n";
}

inline void testRing() {
    auto q = smsd::parseSMARTS("[R]");
    auto mol = smsd::parseSMILES("C1CCC1");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("CCC");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testRing\n";
}

inline void testNotRing() {
    auto q = smsd::parseSMARTS("[!R]");
    auto mol = smsd::parseSMILES("CCC");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testNotRing\n";
}

inline void testOrExpression() {
    // [#6,#7] = carbon or nitrogen
    auto q = smsd::parseSMARTS("[#6,#7]");
    auto molC = smsd::parseSMILES("C");
    auto molN = smsd::parseSMILES("N");
    auto molO = smsd::parseSMILES("O");
    assert(!smsd::Substructure::smarts(q).match(molC).matches.empty());
    assert(!smsd::Substructure::smarts(q).match(molN).matches.empty());
    assert(smsd::Substructure::smarts(q).match(molO).matches.empty());
    std::cout << "  [PASS] testOrExpression\n";
}

inline void testAndExpression() {
    // [#6&D2] = carbon with degree 2
    auto q = smsd::parseSMARTS("[#6&D2]");
    auto mol = smsd::parseSMILES("CCC");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());  // middle C is #6 and D2
    std::cout << "  [PASS] testAndExpression\n";
}

inline void testCharge() {
    auto q = smsd::parseSMARTS("[+1]");
    auto mol = smsd::parseSMILES("[NH4+]");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("C");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testCharge\n";
}

inline void testTwoLetterElement() {
    auto q = smsd::parseSMARTS("[Cl]");
    auto mol = smsd::parseSMILES("CCl");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testTwoLetterElement\n";
}

inline void testPattern_CC() {
    auto q = smsd::parseSMARTS("CC");
    auto mol = smsd::parseSMILES("CCC");  // has C-C substructure
    auto matches = smsd::Substructure::smarts(q).matchAll(mol, 100).matches;
    assert(!matches.empty());
    std::cout << "  [PASS] testPattern_CC (" << matches.size() << " matches)\n";
}

inline void testAromatic() {
    auto q = smsd::parseSMARTS("[a]");
    auto mol = smsd::parseSMILES("c1ccccc1");  // benzene
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("C1CCCCC1");  // cyclohexane
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testAromatic\n";
}

inline void testRingBond() {
    auto q = smsd::parseSMARTS("[#6]@[#6]");
    auto mol = smsd::parseSMILES("C1CCC1");  // cyclobutane: ring bond
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("CCC");  // no ring
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testRingBond\n";
}

inline void testBranches() {
    // C(=O)O = carboxylic acid
    auto q = smsd::parseSMARTS("C(=O)O");
    auto mol = smsd::parseSMILES("CC(=O)O");  // acetic acid
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testBranches (carboxylic acid)\n";
}

inline void testNotOperator() {
    // [!#6] = not carbon
    auto q = smsd::parseSMARTS("[!#6]");
    auto molC = smsd::parseSMILES("C");
    auto molN = smsd::parseSMILES("N");
    assert(smsd::Substructure::smarts(q).match(molC).matches.empty());
    assert(!smsd::Substructure::smarts(q).match(molN).matches.empty());
    std::cout << "  [PASS] testNotOperator\n";
}

inline void testRingClosure() {
    // c1ccccc1 = aromatic 6-ring
    auto q = smsd::parseSMARTS("c1ccccc1");
    auto mol = smsd::parseSMILES("c1ccccc1");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testRingClosure\n";
}

inline void testEmptySmarts() {
    auto q = smsd::parseSMARTS("");
    auto mol = smsd::parseSMILES("C");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());  // empty pattern matches anything
    std::cout << "  [PASS] testEmptySmarts\n";
}

inline void testInvalidSmarts() {
    bool threw = false;
    try {
        smsd::parseSMARTS("[#6");  // missing ]
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
    std::cout << "  [PASS] testInvalidSmarts\n";
}

inline void testSemicolonAnd() {
    // [#6;R] = carbon AND in ring (semicolon = low-priority AND)
    auto q = smsd::parseSMARTS("[#6;R]");
    auto mol = smsd::parseSMILES("C1CCC1");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("CCC");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testSemicolonAnd\n";
}

inline void testMultipleMatches() {
    auto q = smsd::parseSMARTS("[#6]");
    auto mol = smsd::parseSMILES("CCCC");
    auto matches = smsd::Substructure::smarts(q).matchAll(mol, 100).matches;
    assert(static_cast<int>(matches.size()) == 4);  // 4 carbons
    std::cout << "  [PASS] testMultipleMatches (" << matches.size() << " matches)\n";
}

inline void testRecursiveSmarts() {
    // [$([#6][#8])] = carbon bonded to oxygen
    auto q = smsd::parseSMARTS("[$([#6][#8])]");
    auto mol = smsd::parseSMILES("CCO");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());  // C bonded to O

    auto mol2 = smsd::parseSMILES("CCC");
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());
    std::cout << "  [PASS] testRecursiveSmarts\n";
}

inline void testHeteroNeighbors() {
    // z1 in methylamine: C has 1 N neighbour
    auto q1 = smsd::parseSMARTS("[z1]");
    auto mol1 = smsd::parseSMILES("CN");
    assert(!smsd::Substructure::smarts(q1).match(mol1).matches.empty());

    // z0 in ethane: only carbon neighbours
    auto q2 = smsd::parseSMARTS("[z0]");
    auto mol2 = smsd::parseSMILES("CC");
    assert(!smsd::Substructure::smarts(q2).match(mol2).matches.empty());

    std::cout << "  [PASS] testHeteroNeighbors\n";
}

inline void testHeavyDegree() {
    // d2 in propane: middle C has 2 heavy neighbours
    auto q = smsd::parseSMARTS("[d2]");
    auto mol = smsd::parseSMILES("CCC");
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());
    std::cout << "  [PASS] testHeavyDegree\n";
}

inline void testHybridization() {
    // ^2 matches SP2 (aromatic carbon)
    auto q1 = smsd::parseSMARTS("[^2]");
    auto mol1 = smsd::parseSMILES("c1ccccc1");
    assert(!smsd::Substructure::smarts(q1).match(mol1).matches.empty());

    // ^3 matches SP3 (methane)
    auto q2 = smsd::parseSMARTS("[^3]");
    auto mol2 = smsd::parseSMILES("C");
    assert(!smsd::Substructure::smarts(q2).match(mol2).matches.empty());

    // ^1 matches SP (acetylene)
    auto q3 = smsd::parseSMARTS("[^1]");
    auto mol3 = smsd::parseSMILES("C#C");
    assert(!smsd::Substructure::smarts(q3).match(mol3).matches.empty());

    std::cout << "  [PASS] testHybridization\n";
}

inline void testRangeQuery() {
    // D{2-4} matches degree 2 through 4
    auto q = smsd::parseSMARTS("[D{2-4}]");
    auto mol = smsd::parseSMILES("CC(C)(C)C"); // central C has D4
    assert(!smsd::Substructure::smarts(q).match(mol).matches.empty());

    auto mol2 = smsd::parseSMILES("C"); // D0
    assert(smsd::Substructure::smarts(q).match(mol2).matches.empty());

    // D{2-} means >= 2
    auto q2 = smsd::parseSMARTS("[D{2-}]");
    auto mol3 = smsd::parseSMILES("CCC"); // middle C has D2
    assert(!smsd::Substructure::smarts(q2).match(mol3).matches.empty());

    std::cout << "  [PASS] testRangeQuery\n";
}

inline void runAllTests() {
    std::cout << "smarts_parser tests:\n";
    testWildcard();
    testElementMatch();
    testOrganicShorthand();
    testBondSingle();
    testBondDouble();
    testBondTriple();
    testBondAny();
    testDegree();
    testRing();
    testNotRing();
    testOrExpression();
    testAndExpression();
    testCharge();
    testTwoLetterElement();
    testPattern_CC();
    testAromatic();
    testRingBond();
    testBranches();
    testNotOperator();
    testRingClosure();
    testEmptySmarts();
    testInvalidSmarts();
    testSemicolonAnd();
    testMultipleMatches();
    testRecursiveSmarts();
    testHeteroNeighbors();
    testHeavyDegree();
    testHybridization();
    testRangeQuery();
    std::cout << "All smarts_parser tests passed.\n";
}

} // namespace smsd_smarts_test


int main() { smsd_smarts_test::runAllTests(); }

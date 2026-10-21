/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/substructure.hpp"
#include <cassert>
#include <functional>
#include <iostream>
#include <map>
#include <set>

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

// An uppercase element in brackets is aliphatic and a lowercase one aromatic, as outside
// brackets; an atomic number matches either form.
inline void testBracketElementAromaticity() {
    auto mappings = [](const std::string& smarts, const std::string& smiles) {
        std::set<std::map<int, int>> found;
        for (const auto& match : smsd::Substructure::smarts(smarts).matchAll(smiles).matches) found.insert(match.atomMapping);
        return found;
    };
    const std::string toluene = "Cc1ccccc1";
    assert(mappings("[C]", toluene) == (std::set<std::map<int, int>>{{{0, 0}}}));
    assert(mappings("[c]", toluene).size() == 6 && !mappings("[c]", toluene).count({{0, 0}}));
    assert(mappings("[#6]", toluene).size() == 7);
    assert(mappings("[C,c]", toluene).size() == 7);
    assert(mappings("[C;a]", toluene).empty());
    for (const char* smiles : {"Cc1ccccc1", "C1CCCCC1", "c1ccc2[nH]ccc2c1", "CC(=O)Nc1ccc(O)cc1", "c1cc[se]c1CC[Se]C"}) {
        assert(mappings("[C]", smiles) == mappings("C", smiles));
        assert(mappings("[c]", smiles) == mappings("c", smiles));
        assert(mappings("[$([C])]", smiles) == mappings("C", smiles));
    }
    assert(smsd::Substructure::smarts("[NX3]").match("c1cc[nH]c1").matches.empty());
    assert(!smsd::Substructure::smarts("[nX3]").match("c1cc[nH]c1").matches.empty());
    assert(!smsd::Substructure::smarts("[#7X3]").match("c1cc[nH]c1").matches.empty());
    assert(mappings("[Se]", "c1cc[se]c1CC[Se]C") == (std::set<std::map<int, int>>{{{0, 7}}}));
    assert(mappings("[se]", "c1cc[se]c1CC[Se]C") == (std::set<std::map<int, int>>{{{0, 3}}}));
    std::cout << "  [PASS] testBracketElementAromaticity\n";
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

inline std::set<std::map<int, int>> mappingSet(const smsd::SearchResult& result) {
    std::set<std::map<int, int>> mappings;
    for (const auto& match : result.matches) mappings.insert(match.atomMapping);
    return mappings;
}

/** Compiled-pattern screening, shared domains and recursive caches keep search results. */
inline void testPreparedPatternParity() {
    // Element and aromaticity counts reject only targets without an embedding.
    assert(smsd::Substructure::smarts("N~*~N").match("CCN").matches.empty());
    assert(!smsd::Substructure::smarts("N~*~N").match("NCN").matches.empty());
    assert(smsd::Substructure::smarts("[#7].[#7]").match("CN").matches.empty());
    assert(!smsd::Substructure::smarts("[#7].[#7]").match("NCCN").matches.empty());
    assert(smsd::Substructure::smarts("c1ccccc1").match("C1CCCCC1").matches.empty());
    assert(smsd::Substructure::smarts("c1ccccc1").matchAll("Cc1ccccc1").matches.size() == 12);

    // Atom classes never restrict matching, so numbered and plain atoms share domains.
    const auto toluene = smsd::parseSMILES("Cc1ccccc1");
    assert(mappingSet(smsd::Substructure::smarts("[cH:1][cH:2]").matchAll(toluene))
           == mappingSet(smsd::Substructure::smarts("[cH][cH]").matchAll(toluene)));

    // A recursive predicate parsed for each evaluation agrees with the compiled one.
    auto compiled = smsd::parseSMARTS("[C;$(C(=O)O)]~[$([#8;$(O[H]),$(O=*)])]");
    auto transient = compiled;
    std::function<void(smsd::smarts_detail::AtomExpr&)> clear = [&](smsd::smarts_detail::AtomExpr& expression) {
        expression.cachedRecursive.reset();
        for (auto& child : expression.children) clear(child);
    };
    for (auto& atom : transient.atoms) clear(atom.expr);
    for (const char* smiles : {"CC(=O)O", "CC(=O)N", "OC(=O)CC(=O)O", "C(=O)OC=O"}) {
        const auto target = smsd::parseSMILES(smiles);
        assert(mappingSet(smsd::Substructure::smarts(compiled).matchAll(target))
               == mappingSet(smsd::Substructure::smarts(transient).matchAll(target)));
    }

    // Separately compiled copies of one subpattern are each evaluated correctly.
    assert(!smsd::Substructure::smarts("[$(C=O)]~[$(C=O)]").match("O=CC=O").matches.empty());
    assert(smsd::Substructure::smarts("[$(C=O)]~[$(C=O)]").match("O=CCC=O").matches.empty());

    // A truncated enumeration still reports sorted, valid embeddings.
    const auto naphthalene = smsd::parseSMILES("c1ccc2ccccc2c1");
    const auto ring = smsd::Substructure::smarts("*1*****1");
    const auto limited = ring.matchAll(naphthalene, 5);
    const auto all = mappingSet(ring.matchAll(naphthalene));
    assert(limited.status == smsd::SearchStatus::LimitReached && limited.matches.size() == 5 && !limited.complete);
    for (size_t index = 0; index < limited.matches.size(); ++index) {
        assert(all.count(limited.matches[index].atomMapping));
        assert(index == 0 || limited.matches[index - 1].atomMapping < limited.matches[index].atomMapping);
        assert(limited.matches[index].fragments.size() == 1 && limited.matches[index].bondMapping.size() == 6);
    }

    // Induced embeddings exclude extra target bonds.
    smsd::SubstructureOptions induced;
    induced.induced = true;
    assert(!smsd::Substructure::smarts("C~C~C").match("C1CC1").matches.empty());
    assert(smsd::Substructure::smarts("C~C~C").withOptions(induced).match("C1CC1").matches.empty());
    std::cout << "  [PASS] testPreparedPatternParity\n";
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
    testBracketElementAromaticity();
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
    testPreparedPatternParity();
    std::cout << "All smarts_parser tests passed.\n";
}

} // namespace smsd_smarts_test


int main() { smsd_smarts_test::runAllTests(); }

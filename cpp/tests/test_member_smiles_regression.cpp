/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mol_graph.hpp"
#include "smsd/smiles_parser.hpp"

#include <iostream>
#include <stdexcept>

using smsd_test::require;

static smsd::MolGraph atom(int number, int isotope, int charge, int hydrogens, int atomClass = 0) {
    return smsd::MolGraph::Builder().atomCount(1).atomicNumbers({number})
        .massNumbers({isotope}).formalCharges({charge}).hydrogenCounts({hydrogens})
        .atomClasses({atomClass}).setNeighbors({{}}).build();
}

static void testAtomProperties() {
    require(atom(6, 13, 0, 4).toCanonicalSmiles() == "[13CH4]", "member writer lost isotope or H");
    require(atom(7, 0, 1, 4).toCanonicalSmiles() == "[NH4+]", "member writer lost charged atom H");
    require(atom(6, 0, 0, 0).toCanonicalSmiles() == "[C]", "member writer changed supplied zero H");
    require(atom(6, 0, 0, 2).toCanonicalSmiles() == "[CH2]", "member writer changed radical H");
    require(atom(6, 0, 0, 3, 7).toCanonicalSmiles() == "[CH3:7]", "member writer lost atom class");
}

static void testTetrahedralStereo() {
    const auto graph = smsd::MolGraph::Builder().atomCount(4).atomicNumbers({6, 6, 9, 17})
        .hydrogenCounts({3, 1, 0, 0}).setNeighbors({{1}, {0, 2, 3}, {1}, {1}})
        .setBondOrders({{1}, {1, 1, 1}, {1}, {1}}).tetrahedralChirality({0, 1, 0, 0}).build();
    const auto text = graph.toCanonicalSmiles();
    require(text.find('@') != std::string::npos && text.find('H') != std::string::npos,
            "member writer lost tetrahedral hydrogen");
}

static void testDoubleBondStereo() {
    std::vector<std::vector<int>> configurations(4, std::vector<int>(4, 0));
    configurations[1][2] = configurations[2][1] = 2;
    const auto graph = smsd::MolGraph::Builder().atomCount(4).atomicNumbers({6, 6, 6, 6})
        .hydrogenCounts({3, 1, 1, 3}).setNeighbors({{1}, {0, 2}, {1, 3}, {2}})
        .setBondOrders({{1}, {1, 2}, {2, 1}, {1}}).doubleBondStereo(configurations).build();
    const auto text = graph.toCanonicalSmiles();
    require(text.find('/') != std::string::npos || text.find('\\') != std::string::npos,
            "member writer lost alkene configuration");
}

static void testParserWriterParity() {
    for (const char* text : {"[13CH4]", "[NH4+]", "[C]", "[CH2]", "[CH3:7]",
            "C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", "C/C=C/C", "C/C=C\\C",
            "C/C=C(/F)\\C(Cl)=C/C", "C1[C@H](F)CCC[C@@H]1F",
            "C[N@](CC)CCC", "c1ccc(-c2ccccc2)cc1", "[nH]1cccc1", "CCO.[Na+]"}) {
        const auto graph = smsd::parseSMILES(text);
        const auto member = graph.toCanonicalSmiles();
        require(member == smsd::toSMILES(graph), std::string("member and free writer differ: ") + text);
        require(smsd::parseSMILES(member).getCanonicalHash() == graph.getCanonicalHash(),
                std::string("member writer changed molecule: ") + text);
    }
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"atom-properties", testAtomProperties}, {"tetrahedral", testTetrahedralStereo},
        {"double-bond", testDoubleBondStereo}, {"parity", testParserWriterParity}});
}

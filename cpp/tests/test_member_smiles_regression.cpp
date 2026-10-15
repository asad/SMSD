/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/mol_graph.hpp"

#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

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

int main(int argc, char** argv) {
    int failures = 0;
    const std::pair<const char*, void (*)()> tests[] = {
        {"atom-properties", testAtomProperties}, {"tetrahedral", testTetrahedralStereo},
        {"double-bond", testDoubleBondStereo}
    };
    for (const auto& test : tests) {
        if (argc > 1 && argv[1] != std::string(test.first)) continue;
        try { test.second(); std::cout << "PASS " << test.first << '\n'; }
        catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}

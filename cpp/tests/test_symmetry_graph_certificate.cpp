/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/detail/search_graph_validation.hpp"
#include "smsd/smiles_parser.hpp"

#include <cassert>
#include <iostream>

static bool certified(const smsd::MolGraph& graph) {
    return smsd::search_detail::symmetryGraphCoherent(graph, [] { return false; });
}

int main() {
    const auto ordinary = smsd::parseSMILES("CC(C)CO");
    assert(certified(ordinary));
    auto repeated = ordinary;
    repeated.neighbors[0].push_back(1);
    assert(!certified(repeated));
    repeated.degree[0] = static_cast<int>(repeated.neighbors[0].size());
    assert(!certified(repeated));
    auto selfListing = ordinary;
    selfListing.neighbors[0].push_back(0);
    ++selfListing.degree[0];
    assert(!certified(selfListing));
    auto wrongDegree = ordinary;
    wrongDegree.degree[0] = -1;
    assert(!certified(wrongDegree));
    wrongDegree.degree[0] = 0;
    assert(!certified(wrongDegree));
    auto asymmetric = ordinary;
    asymmetric.bondPacked[0][1] ^= smsd::MolGraph::BOND_AROM_BIT;
    assert(!certified(asymmetric));
    asymmetric = ordinary;
    asymmetric.bondPacked[0][1] = 2;
    assert(!certified(asymmetric));
    asymmetric = ordinary;
    asymmetric.bondPacked[0][1] ^= smsd::MolGraph::BOND_RING_BIT;
    assert(!certified(asymmetric));
    auto hidden = ordinary;
    hidden.bondPacked[0][2] = hidden.bondPacked[2][0] = 1;
    assert(!certified(hidden));
    auto missing = ordinary;
    missing.bondPacked[0][1] = missing.bondPacked[1][0] = 0;
    assert(!certified(missing));
    auto shortRow = ordinary;
    shortRow.bondPacked[1].clear();
    assert(!certified(shortRow));
    auto outOfRange = ordinary;
    outOfRange.neighbors[0][0] = ordinary.n;
    assert(!certified(outOfRange));
    auto missingAtoms = ordinary;
    missingAtoms.atomicNum.clear();
    assert(!certified(missingAtoms));
    auto shortStereo = ordinary;
    shortStereo.hasDbStereo = true;
    shortStereo.dbStereoConf.clear();
    assert(!certified(shortStereo));
    auto sparse = ordinary;
    sparse.useDense = false;
    sparse.sparseBondProps.clear();
    for (int atom = 0; atom < ordinary.n; ++atom)
        for (int neighbour : ordinary.neighbors[atom])
            if (atom < neighbour)
                sparse.sparseBondProps[smsd::MolGraph::bondKey(atom, neighbour)] = {
                    ordinary.bondOrder(atom, neighbour), int(ordinary.bondInRing(atom, neighbour)),
                    int(ordinary.bondAromatic(atom, neighbour))};
    assert(certified(sparse));
    auto negativeOrder = sparse;
    negativeOrder.sparseBondProps[smsd::MolGraph::bondKey(0, 1)][0] = -1;
    assert(!certified(negativeOrder));
    auto invalidKey = sparse;
    invalidKey.sparseBondProps[-1] = {1, 0, 0};
    assert(!certified(invalidKey));
    sparse.sparseBondProps[smsd::MolGraph::bondKey(0, 2)] = {1, 0, 0};
    assert(!certified(sparse));
    auto changed = smsd::parseSMILES("CC.CC");
    changed.setBondOrder(1, 2, 1);
    assert(certified(changed));
    changed.setBondOrder(1, 2, 0);
    assert(certified(changed));
    int completePolls = 0;
    assert(smsd::search_detail::symmetryGraphCoherent(ordinary, [&] { ++completePolls; return false; }));
    for (int stopAt = 1; stopAt <= completePolls; ++stopAt) {
        int observed = 0;
        assert(!smsd::search_detail::symmetryGraphCoherent(ordinary, [&] { return ++observed >= stopAt; }));
    }
    int polls = 0;
    assert(!smsd::search_detail::symmetryGraphCoherent(ordinary, [&] { return ++polls == 3; }));
    assert(certified(ordinary));
    std::cout << "symmetry certificate controls passed\n";
}

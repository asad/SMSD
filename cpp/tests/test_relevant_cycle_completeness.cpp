/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/substructure.hpp"

#include <stdexcept>

#include <string>
#include <vector>

static smsd::MolGraph cyclobutaneMacrocycle(int units) {
    const int chain = 2 * units, atoms = 3 * units + 1 + chain;
    std::vector<std::vector<int>> neighbours(static_cast<size_t>(atoms)), orders(static_cast<size_t>(atoms));
    auto bond = [&](int first, int second) {
        neighbours[first].push_back(second);
        neighbours[second].push_back(first);
        orders[first].push_back(1);
        orders[second].push_back(1);
    };
    for (int unit = 0; unit < units; ++unit) {
        const int corner = 3 * unit, opposite = 3 * (unit + 1);
        bond(corner, corner + 1);
        bond(corner + 1, opposite);
        bond(corner, corner + 2);
        bond(corner + 2, opposite);
    }
    int previous = 3 * units;
    for (int atom = 3 * units + 1; atom < atoms; ++atom) {
        bond(previous, atom);
        previous = atom;
    }
    bond(previous, 0);
    return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(static_cast<size_t>(atoms), 6))
        .setNeighbors(neighbours).setBondOrders(orders).build();
}

int main() {
    struct Example { const char* text; size_t rings; };
    for (const auto& example : {Example{"C1CCCCC1", 1}, {"C12C3C4C1C5C4C3C25", 6},
                               {"C1C2CC3CC1CC(C2)C3", 4}}) {
        const auto graph = smsd::parseSMILES(example.text);
        const auto perceived = smsd::detail_rf::perceiveRingSets(graph);
        smsd_test::require(!perceived.truncated && perceived.rcb.size() == example.rings,
                          std::string("ordinary relevant cycles changed: ") + example.text);
        bool truncated = true;
        const auto member = graph.computeRelevantCyclesReportingTruncation(truncated);
        smsd_test::require(!truncated && member.size() == example.rings,
                          std::string("ordinary member relevant cycles changed: ") + example.text);
    }
    const auto exact = cyclobutaneMacrocycle(5);
    smsd_test::require(!smsd::relevantCyclesTruncated(exact)
                      && smsd::computeRelevantCycles(exact).size() == 5 + (size_t(1) << 5),
                      "small diamond chains retain every relevant cycle");
    const auto graph = cyclobutaneMacrocycle(13);
    const size_t expected = 13 + (size_t(1) << 13);
    const auto perceived = smsd::detail_rf::perceiveRingSets(graph);
    smsd_test::require(perceived.relevantTruncated && !perceived.basisTruncated && perceived.rcb.size() <= expected,
                      "the bounded large chain retains a complete basis and reports a partial relevant set");
    bool truncated = false;
    const auto member = graph.computeRelevantCyclesReportingTruncation(truncated);
    smsd_test::require(truncated || member.size() == expected,
                      "a partial member relevant-cycle list must report incomplete perception");
    const auto numeric = smsd::Substructure::smarts("[R1]").match(graph);
    smsd_test::require(numeric.status == smsd::SearchStatus::Completed && numeric.optimal && numeric.complete,
                      "numeric ring membership remains exact when only relevant enumeration is bounded");
    const auto ordinaryPattern = smsd::Substructure::smarts("[R1]").matchAll("C1CCCCC1");
    smsd_test::require(ordinaryPattern.complete && ordinaryPattern.optimal
                      && ordinaryPattern.matches.size() == 6,
                      "Numeric SMARTS ring membership must retain ordinary complete results");
    return 0;
}

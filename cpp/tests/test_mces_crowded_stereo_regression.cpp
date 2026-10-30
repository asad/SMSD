/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static smsd::MolGraph configuredGraph(const std::vector<int>& elements,
                                      const std::vector<std::vector<int>>& neighbours,
                                      int configuration = 1) {
    const int atomCount = static_cast<int>(elements.size());
    std::vector<std::vector<int>> orders(atomCount, std::vector<int>(atomCount, 0));
    std::vector<std::vector<int>> configurations(atomCount, std::vector<int>(atomCount, 0));
    for (int atom = 0; atom < atomCount; ++atom)
        for (int neighbour : neighbours[atom])
            orders[atom][neighbour] = (atom == 0 && neighbour == 1) || (atom == 1 && neighbour == 0) ? 2 : 1;
    configurations[0][1] = configurations[1][0] = configuration;
    return smsd::MolGraph::Builder().atomCount(atomCount).atomicNumbers(elements)
        .hydrogenCounts(std::vector<int>(atomCount, 0)).setNeighbors(neighbours)
        .setBondOrders(orders).doubleBondStereo(configurations).build(false);
}

/** A crowded endpoint must retain every common bond when its full selection agrees. */
static void testCrowdedDoubleBondFullSelection() {
    const auto query = configuredGraph({15, 7, 6, 16, 9}, {{1, 2, 3}, {0, 4}, {0, 3}, {0, 2}, {1}});
    const auto target = configuredGraph({15, 7, 17, 6, 16, 9},
        {{1, 2, 3, 4}, {0, 5}, {0}, {0, 4}, {0, 3}, {1}}, 2);
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.connectivity = smsd::MCSConnectivity::Connected;
    options.nodeLimit = 0;
    options.timeoutMs = 10000;
    options.chemistry.matchAtomType = true;
    options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    options.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    options.chemistry.useBondStereo = true;
    options.chemistry.useChirality = false;
    options.chemistry.tautomerAware = false;
    options.chemistry.ringMatchesRingOnly = false;
    options.chemistry.completeRingsOnly = false;
    options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::IGNORE;

    smsd::Match expected;
    expected.atomMapping = {{0, 0}, {1, 1}, {2, 3}, {3, 4}, {4, 5}};
    expected.bondMapping = {{0, 1, 0, 1}, {0, 2, 0, 3}, {0, 3, 0, 4}, {1, 4, 1, 5}, {2, 3, 3, 4}};
    require(smsd::MCS::validateMapping(query, target, expected, options).empty(),
            "The full five-bond, five-atom witness must satisfy the stereo contract");

    const auto result = smsd::MCS::find(query).withOptions(options).match(target);
    require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
            "The full-selection control must complete and prove its objective");
    require(!result.matches.empty(), "The full-selection control must have a common subgraph");
    for (const auto& match : result.matches) {
        require(smsd::MCS::validateMapping(query, target, match, options).empty(),
                "The full-selection witness must be valid");
        require(result.score == 5 && match.bondMapping.size() == 5 && match.atomMapping.size() == 5,
                "An agreeing crowded selection must retain all five bonds and atoms");
    }
    std::cout << "MCES crowded stereo full-set control passed" << std::endl;
}

/** Dropping the first ligand can retain a connected MCES with a different stereo reference. */
static void testCrowdedDoubleBondSelection() {
    const auto query = configuredGraph({15, 7, 6, 16, 9}, {{1, 2, 3}, {0, 4}, {0, 3}, {0, 2}, {1}});
    const auto target = configuredGraph({15, 7, 17, 6, 16, 9},
        {{1, 2, 3, 4}, {0, 5}, {0}, {0, 4}, {0, 3}, {1}});
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.connectivity = smsd::MCSConnectivity::Connected;
    options.nodeLimit = 0;
    options.timeoutMs = 10000;
    options.chemistry.matchAtomType = true;
    options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    options.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    options.chemistry.useBondStereo = true;
    options.chemistry.useChirality = false;
    options.chemistry.tautomerAware = false;
    options.chemistry.ringMatchesRingOnly = false;
    options.chemistry.completeRingsOnly = false;
    options.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::IGNORE;

    smsd::Match expected;
    expected.atomMapping = {{0, 0}, {1, 1}, {2, 3}, {3, 4}, {4, 5}};
    expected.bondMapping = {{0, 1, 0, 1}, {0, 3, 0, 4}, {1, 4, 1, 5}, {2, 3, 3, 4}};
    require(smsd::MCS::validateMapping(query, target, expected, options).empty(),
            "The four-bond, five-atom connected witness must satisfy the stereo contract");

    const auto result = smsd::MCS::find(query).withOptions(options).match(target);
    require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
            "The small configured-graph search must complete and prove its objective");
    require(!result.matches.empty(), "The configured graphs must have a common subgraph");
    for (const auto& match : result.matches) {
        require(smsd::MCS::validateMapping(query, target, match, options).empty(),
                "The returned selected-bond witness must be valid");
        require(result.score == 4 && match.bondMapping.size() == 4 && match.atomMapping.size() == 5,
                "MCES must retain four bonds and all five atoms after dropping the first ligand");
    }
}

int main() {
    testCrowdedDoubleBondFullSelection();
    testCrowdedDoubleBondSelection();
    std::cout << "MCES crowded stereo regression passed\n";
}

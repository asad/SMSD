/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */

#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

static void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

static void removeBond() {
    auto target = smsd::parseSMILES("CC");
    target.setBondOrder(0, 1, 0);
    smsd::MCSOptions options;
    options.connectivity = smsd::MCSConnectivity::Fragments;
    options.nodeLimit = 10000;
    options.timeoutMs = 10000;
    const auto query = smsd::parseSMILES("C.C.N");
    const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 10);
    require(result.score == 2 && result.optimal && result.complete && result.matches.size() == 2,
            "Removing a bond must permit both optimal induced isolated-atom mappings");
    for (const auto& match : result.matches)
        require(smsd::MCS::validateMapping(query, target, match, options).empty(), "Invalid removed-bond witness");
}

static void queryTopology() {
    auto query = smsd::parseSMILES("C1CC1");
    query.setBondOrder(0, 1, 0);
    const auto result = smsd::Substructure::find(query).matchAll(smsd::parseSMILES("CCC"), 10);
    require(result.complete && result.matches.size() == 2,
            "Removing a query bond must update degree domains and neighbour anchors");
}

static void restoreBond() {
    auto target = smsd::parseSMILES("CC");
    target.setBondOrder(0, 1, 0);
    target.setBondOrder(0, 1, 1);
    const auto result = smsd::Substructure::find("CC").matchAll(target, 10);
    require(result.complete && result.matches.size() == 2, "Restoring a bond must restore its topology");
    auto added = smsd::parseSMILES("C.C");
    added.setBondOrder(0, 1, 1);
    require(smsd::Substructure::find("CC").matchAll(added, 10).matches.size() == 2,
            "Adding a dense bond must create its neighbour entries");
}

static void ringsAndPositiveOrder() {
    auto triangle = smsd::parseSMILES("C1CC1");
    triangle.ensureRingCounts();
    require(!triangle.computeRings().empty(), "Triangle cycle missing");
    triangle.setBondOrder(0, 1, 0);
    triangle.ensureRingCounts();
    require(triangle.computeRings().empty()
            && std::all_of(triangle.ringCount.begin(), triangle.ringCount.end(), [](int n) { return n == 0; }),
            "Removed ring bond must invalidate stored cycle and count values");
    auto ordinary = smsd::parseSMILES("CC");
    const auto neighbors = ordinary.neighbors;
    const auto adjacency = ordinary.adjLong;
    ordinary.setBondOrder(0, 1, 2);
    require(ordinary.neighbors == neighbors && ordinary.adjLong == adjacency
            && !smsd::Substructure::find("C=C").match(ordinary).matches.empty(),
            "Positive order changes must preserve topology and update bond matching");
}

static void sparseTopology() {
    const int count = smsd::MolGraph::SPARSE_THRESHOLD + 1;
    std::vector<int> labels(count, 7);
    labels[0] = labels[1] = labels[2] = 6;
    std::vector<std::vector<int>> triangle(count), chain(count);
    triangle[0] = {1, 2}; triangle[1] = {0, 2}; triangle[2] = {0, 1};
    chain[0] = {1}; chain[1] = {0, 2}; chain[2] = {1};
    auto query = smsd::MolGraph::Builder().atomCount(count).atomicNumbers(labels)
        .setNeighbors(triangle).build(false);
    const auto target = smsd::MolGraph::Builder().atomCount(count).atomicNumbers(labels)
        .setNeighbors(chain).build(false);
    require(!query.useDense, "Sparse boundary fixture is not sparse");
    query.setBondOrder(0, 1, 0);
    require(!smsd::Substructure::find(query).withTimeout(10000).match(target).matches.empty(),
            "Sparse removal must update full-query degree domains");
    query.setBondOrder(0, 1, 1);
    const auto restored = smsd::MolGraph::Builder().atomCount(count).atomicNumbers(labels)
        .setNeighbors(triangle).build(false);
    require(!smsd::Substructure::find(query).withTimeout(10000).match(restored).matches.empty(),
            "Restoring an existing sparse entry must restore its topology");
    query.setBondOrder(3, 4, 1);
    require(!query.hasBond(3, 4), "A missing sparse entry must retain its no-op setter contract");
}

int main() {
    int failures = 0;
    for (const auto& test : {&removeBond, &queryTopology, &restoreBond, &ringsAndPositiveOrder, &sparseTopology}) {
        try { test(); }
        catch (const std::exception& error) { ++failures; std::cerr << error.what() << '\n'; }
    }
    std::cout << "groups=5 failures=" << failures << '\n';
    return failures == 0 ? 0 : 1;
}

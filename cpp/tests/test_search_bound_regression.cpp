/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>

using smsd_test::require;

static smsd::MolGraph graph(int mask) {
    std::vector<std::vector<int>> adjacent(3);
    int bit = 0;
    for (int first = 0; first < 3; ++first)
        for (int second = first + 1; second < 3; ++second, ++bit)
            if ((mask & (1 << bit)) != 0) {
                adjacent[first].push_back(second); adjacent[second].push_back(first);
            }
    return smsd::MolGraph::Builder().atomCount(3).atomicNumbers({6, 6, 6}).setNeighbors(adjacent).build(false);
}

using Score = std::tuple<int, int, int>;
using Mapping = std::map<int, int>;

static std::pair<Score, std::set<Mapping>> oracle(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                               const smsd::MCSOptions& options) {
    Score best{0, 0, 0};
    std::set<Mapping> matches;
    Mapping mapping;
    std::set<int> used;
    std::function<void(int)> enumerate = [&](int atom) {
        if (atom != query.n) {
            enumerate(atom + 1);
            for (int candidate = 0; candidate < target.n; ++candidate) {
                if (used.count(candidate) || std::find(options.excludedTargetAtoms.begin(), options.excludedTargetAtoms.end(), candidate)
                    != options.excludedTargetAtoms.end()) continue;
                bool compatible = true;
                for (const auto& previous : mapping)
                    if (query.hasBond(atom, previous.first) != target.hasBond(candidate, previous.second)) compatible = false;
                if (!compatible) continue;
                mapping.emplace(atom, candidate); used.insert(candidate);
                enumerate(atom + 1);
                used.erase(candidate); mapping.erase(atom);
            }
            return;
        }
        if (mapping.empty()) return;
        int bonds = 0, weight = 0;
        for (const auto& first : mapping) {
            if (!options.atomWeights.empty()) weight += static_cast<int>(options.atomWeights[first.first] * 4);
            for (const auto& second : mapping)
                if (first.first < second.first && query.hasBond(first.first, second.first)) ++bonds;
        }
        if (options.connectivity == smsd::MCSConnectivity::Connected) {
            std::set<int> reached{mapping.begin()->first};
            for (;;) {
                const size_t before = reached.size();
                for (const auto& retained : mapping)
                    for (int neighbor : query.neighbors[retained.first])
                        if (mapping.count(neighbor) && reached.count(neighbor)) reached.insert(retained.first);
                if (before == reached.size()) break;
            }
            if (reached.size() != mapping.size()) return;
        }
        const int atoms = static_cast<int>(mapping.size());
        const Score score{options.atomWeights.empty() ? atoms : weight * 250, atoms, bonds};
        if (score < best) return;
        if (best < score) { best = score; matches.clear(); }
        matches.insert(mapping);
    };
    enumerate(0);
    return {best, matches};
}

static void exhaustiveSmallGraphs() {
    for (int source = 0; source < 8; ++source) for (int destination = 0; destination < 8; ++destination) {
        const auto query = graph(source), target = graph(destination);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (const auto& weights : {std::vector<double>{}, {-1.0, 0.25, 1.5}, {0.0, 0.0, 0.0}})
                for (const auto& excluded : {std::vector<int>{}, {2}}) {
                    smsd::MCSOptions options;
                    options.connectivity = connectivity; options.atomWeights = weights;
                    options.excludedTargetAtoms = excluded; options.timeoutMs = 10000;
                    const auto expected = oracle(query, target, options);
                    const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
                    require(result.optimal && result.complete, "Small graph search must prove and enumerate the optimum");
                    const double score = weights.empty() ? std::get<0>(expected.first) : std::get<0>(expected.first) / 1000.0;
                    require(result.score == score, "Objective differs from the independent graph oracle");
                    std::set<Mapping> actual;
                    for (const auto& match : result.matches) {
                        require(smsd::MCS::validateMapping(query, target, match, options).empty(), "Invalid oracle mapping");
                        require(match.atomMapping.size() == static_cast<size_t>(std::get<1>(expected.first))
                            && match.bondMapping.size() == static_cast<size_t>(std::get<2>(expected.first)),
                            "Secondary objective differs from the independent oracle");
                        require(actual.insert(match.atomMapping).second, "Repeated raw embedding");
                    }
                    require(actual == expected.second, "Raw tied embeddings differ from the independent oracle");
                    auto singleOptions = options;
                    singleOptions.nodeLimit = 1000;
                    const auto single = smsd::MCS::find(query).withOptions(singleOptions).match(target);
                    require(single.optimal && single.complete && single.score == score,
                            "Single optimum differs from the independent graph oracle");
                    require(single.matches.size() == (expected.second.empty() ? 0u : 1u),
                            "Single search must retain one optimal witness when one exists");
                    if (!single.matches.empty()) {
                        const auto& match = single.matches.front();
                        require(expected.second.count(match.atomMapping) == 1
                            && match.bondMapping.size() == static_cast<size_t>(std::get<2>(expected.first))
                            && smsd::MCS::validateMapping(query, target, match, options).empty(),
                            "Single witness differs from the independent optimum set");
                    }
                }
    }
}

static void connectedNaturalProductProof() {
    const auto query = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C");
    const auto target = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O");
    smsd::MCSOptions options; options.nodeLimit = 10000; options.timeoutMs = 10000;
    const auto result = smsd::MCS::find(query).withOptions(options).match(target);
    require(result.optimal, "The connected common graph must be proved within the node budget");
    require(result.score == 53 && result.matches.front().bondMapping.size() == 58,
            "The connected natural-product common graph differs from the reference optimum");
    require(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty(), "Invalid natural-product mapping");
}

static void optionalAdjacencyWords() {
    const auto query = graph(7);
    auto target = graph(7);
    target.adjLong.clear();
    smsd::MCSOptions options; options.timeoutMs = 10000;
    const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
    require(result.optimal && result.complete && result.matches.size() == 6,
            "Derived adjacency words must not change the six triangle embeddings");
}

static void candidatesAcrossWordBoundaries() {
    std::vector<int> elements(130, 6);
    elements[63] = elements[64] = elements[65] = 7;
    std::vector<std::vector<int>> adjacent(130);
    for (int atom = 1; atom < 130; ++atom) {
        adjacent[atom - 1].push_back(atom); adjacent[atom].push_back(atom - 1);
    }
    const auto target = smsd::MolGraph::Builder().atomCount(130).atomicNumbers(elements).setNeighbors(adjacent).build(false);
    const auto query = smsd::parseSMILES("NNN");
    smsd::MCSOptions options; options.timeoutMs = 10000;
    auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
    require(result.optimal && result.complete && result.score == 3 && result.matches.size() == 2,
            "Both full nitrogen-chain orientations must survive target-word boundaries");
    options.excludedTargetAtoms = {64};
    result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
    require(result.optimal && result.complete && result.score == 1 && result.matches.size() == 6,
            "An excluded bridge must separate the available connected target atoms");
    options.connectivity = smsd::MCSConnectivity::Fragments;
    result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
    require(result.optimal && result.complete && result.score == 2 && result.matches.size() == 2,
            "Fragment search must retain both separated available atoms");

    smsd::SubstructureOptions substructureOptions; substructureOptions.timeoutMs = 10000;
    result = smsd::Substructure::find(query).withOptions(substructureOptions).matchAll(target, 1000);
    std::set<Mapping> embeddings;
    for (const auto& match : result.matches) embeddings.insert(match.atomMapping);
    const std::set<Mapping> expected{{{0, 63}, {1, 64}, {2, 65}}, {{0, 65}, {1, 64}, {2, 63}}};
    require(result.optimal && result.complete && result.matches.size() == 2 && embeddings == expected,
            "Substructure enumeration must retain both full orientations across target-word boundaries");
    elements[64] = 6;
    const auto interrupted = smsd::MolGraph::Builder().atomCount(130).atomicNumbers(elements).setNeighbors(adjacent).build(false);
    result = smsd::Substructure::find(query).withOptions(substructureOptions).matchAll(interrupted, 1000);
    require(result.optimal && result.complete && result.matches.empty(),
            "A carbon between the nitrogen atoms must prevent a full substructure embedding");
}

using SelectedKey = std::pair<Mapping, std::vector<std::pair<int, int>>>;

static std::pair<Score, std::set<SelectedKey>> edgeOracle(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                                       const smsd::MCSOptions& options) {
    Score best{0, 0, 0};
    std::set<SelectedKey> matches;
    Mapping mapping;
    std::set<int> used;
    std::function<void(int)> enumerate = [&](int atom) {
        if (atom != query.n) {
            enumerate(atom + 1);
            for (int candidate = 0; candidate < target.n; ++candidate) {
                if (used.count(candidate) || std::find(options.excludedTargetAtoms.begin(), options.excludedTargetAtoms.end(), candidate)
                    != options.excludedTargetAtoms.end()) continue;
                mapping.emplace(atom, candidate); used.insert(candidate);
                enumerate(atom + 1);
                used.erase(candidate); mapping.erase(atom);
            }
            return;
        }
        std::vector<std::pair<int, int>> edges;
        for (const auto& first : mapping) for (const auto& second : mapping)
            if (first.first < second.first && query.hasBond(first.first, second.first)
                && target.hasBond(first.second, second.second)) edges.emplace_back(first.first, second.first);
        for (int subset = 1; subset < (1 << edges.size()); ++subset) {
            Mapping incident;
            std::vector<std::pair<int, int>> selected;
            for (size_t edge = 0; edge < edges.size(); ++edge) if ((subset & (1 << edge)) != 0) {
                selected.push_back(edges[edge]);
                incident.emplace(edges[edge].first, mapping.at(edges[edge].first));
                incident.emplace(edges[edge].second, mapping.at(edges[edge].second));
            }
            if (options.connectivity == smsd::MCSConnectivity::Connected || options.maxFragments == 1) {
                std::set<int> reached{incident.begin()->first};
                for (;;) {
                    const size_t before = reached.size();
                    for (const auto& edge : selected)
                        if (reached.count(edge.first) || reached.count(edge.second)) {
                            reached.insert(edge.first); reached.insert(edge.second);
                        }
                    if (before == reached.size()) break;
                }
                if (reached.size() != incident.size()) continue;
            }
            const Score score{static_cast<int>(selected.size()), static_cast<int>(incident.size()), 0};
            if (score < best) continue;
            if (best < score) { best = score; matches.clear(); }
            matches.emplace(std::move(incident), std::move(selected));
        }
    };
    enumerate(0);
    return {best, matches};
}

static void exhaustiveCommonEdges() {
    for (int source = 0; source < 8; ++source) for (int destination = 0; destination < 8; ++destination) {
        const auto query = graph(source), target = graph(destination);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (const auto& excluded : {std::vector<int>{}, {2}}) {
                smsd::MCSOptions options; options.objective = smsd::MCSObjective::MCES;
                options.connectivity = connectivity; options.excludedTargetAtoms = excluded; options.timeoutMs = 10000;
                const auto expected = edgeOracle(query, target, options);
                const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
                require(result.optimal && result.complete && result.score == std::get<0>(expected.first),
                        "MCES differs from the independent selected-edge oracle");
                std::set<SelectedKey> actual;
                for (const auto& match : result.matches) {
                    require(smsd::MCS::validateMapping(query, target, match, options).empty(), "Invalid MCES oracle mapping");
                    require(match.atomMapping.size() == static_cast<size_t>(std::get<1>(expected.first)),
                            "MCES incident-atom tie differs from the oracle");
                    std::vector<std::pair<int, int>> selected;
                    for (const auto& bond : match.bondMapping) selected.emplace_back(bond[0], bond[1]);
                    std::sort(selected.begin(), selected.end());
                    require(actual.emplace(match.atomMapping, selected).second, "Repeated MCES embedding");
                }
                require(actual == expected.second, "MCES raw tied embeddings differ from the independent oracle");
            }
    }
}

static void oneFragmentIsConnected() {
    const auto query = smsd::parseSMILES("CC.CC"), target = smsd::parseSMILES("CCCC");
    for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
        smsd::MCSOptions options; options.objective = objective; options.timeoutMs = 10000;
        const auto connected = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
        options.connectivity = smsd::MCSConnectivity::Fragments; options.maxFragments = 1;
        const auto one = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
        require(one.optimal && one.complete && one.score == connected.score
                && one.matches.size() == connected.matches.size(), "One allowed fragment must have connected-search semantics");
        for (const auto& match : one.matches) require(match.fragments.size() == 1, "One-fragment result is disconnected");
    }
    smsd::MCSOptions options; options.objective = smsd::MCSObjective::MCES;
    options.connectivity = smsd::MCSConnectivity::Fragments; options.maxFragments = 2; options.timeoutMs = 10000;
    const auto two = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
    require(two.optimal && two.complete && two.score == 2 && two.matches.size() == 8,
            "Two disconnected common edges must retain all eight incident-atom embeddings");
}

static void exhaustedFullQueryProof() {
    const auto query = smsd::parseSMILES("C1CCCCCCCCCCC1");
    const auto target = smsd::parseSMILES("C1CCCCCCCCCCCCCC1");
    for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
        smsd::MCSOptions options;
        options.connectivity = connectivity; options.nodeLimit = 600; options.timeoutMs = 10000;
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        require(result.optimal && result.complete && result.status == smsd::SearchStatus::Completed,
                "An exhausted full-query search must retain its near-complete optimum proof");
        require(result.score == 11 && result.matches.front().bondMapping.size() == 10,
                "Unequal cycles must retain the eleven-atom common path");
        require(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty(),
                "Invalid near-complete cycle mapping");
    }
    for (const auto& pair : {std::pair<const char*, const char*>{"CCCC", "C1CCC1"}, {"CC.C", "CCC"}}) {
        const auto source = smsd::parseSMILES(pair.first), destination = smsd::parseSMILES(pair.second);
        const auto all = smsd::MCS::find(source).withTimeout(10000).matchAll(destination, 1000);
        const size_t expected = source.n == 4 ? 16 : 4;
        require(all.optimal && all.complete && all.matches.size() == expected,
                "A near-complete proof must preserve every tied embedding");
        for (const auto& match : all.matches)
            require(smsd::MCS::validateMapping(source, destination, match).empty(),
                    "Invalid near-complete tied mapping");
    }
}

static void heterocycleCommonEdgeProof() {
    const auto query = smsd::parseSMILES("C1CCCNC1");
    const auto target = smsd::parseSMILES("C1CCCCCCCC1");
    for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
        smsd::MCSOptions options;
        options.objective = smsd::MCSObjective::MCES;
        options.connectivity = connectivity;
        options.nodeLimit = 128;
        options.timeoutMs = 10000;
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        require(result.optimal && result.complete && result.status == smsd::SearchStatus::Completed,
                "The heterocycle common-edge optimum must be proved within the node budget");
        require(result.score == 4 && result.matches.front().atomMapping.size() == 5,
                "The common-edge result must retain the five-carbon path");
        require(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty(),
                "Invalid common-edge heterocycle mapping");
    }
}

static void twinCandidateOrder() {
    auto nitrogenStar = [](const std::vector<int>& centreNeighbours) {
        const int atoms = static_cast<int>(centreNeighbours.size()) + 1;
        std::vector<int> elements(atoms, 6);
        elements[0] = 7;
        std::vector<std::vector<int>> neighbours(atoms), orders(atoms);
        neighbours[0] = centreNeighbours;
        orders[0].assign(centreNeighbours.size(), 1);
        for (int atom = 1; atom < atoms; ++atom) {
            neighbours[atom] = {0};
            orders[atom] = {1};
        }
        return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(elements)
            .setNeighbors(neighbours).setBondOrders(orders).build(false);
    };
    const auto query = nitrogenStar({1, 2, 3});
    std::set<Mapping> expected;
    for (int firstLeaf = 1; firstLeaf <= 3; ++firstLeaf)
        for (int secondLeaf = firstLeaf + 1; secondLeaf <= 3; ++secondLeaf) {
            expected.insert(Mapping{{0, 0}, {firstLeaf, 1}, {secondLeaf, 2}});
            expected.insert(Mapping{{0, 0}, {firstLeaf, 2}, {secondLeaf, 1}});
        }
    for (const auto& centreNeighbours : {std::vector<int>{1, 2}, std::vector<int>{2, 1}}) {
        const auto target = nitrogenStar(centreNeighbours);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
            smsd::MCSOptions options;
            options.connectivity = connectivity;
            options.nodeLimit = 1000;
            options.timeoutMs = 10000;
            const auto search = smsd::MCS::find(query).withOptions(options);
            const auto all = search.matchAll(target, 1000);
            require(all.status == smsd::SearchStatus::Completed && all.optimal && all.complete && all.score == 3,
                    "Neighbour order must preserve the three-atom optimum");
            std::set<Mapping> actual;
            for (const auto& match : all.matches) {
                require(match.bondMapping.size() == 2 && match.fragments.size() == 1
                    && smsd::MCS::validateMapping(query, target, match, options).empty(),
                    "Invalid star optimum mapping");
                require(actual.insert(match.atomMapping).second, "Repeated star optimum mapping");
            }
            require(actual == expected, "Every one of the six labelled star optima must be retained");
            const auto single = search.match(target);
            require(single.status == smsd::SearchStatus::Completed && single.optimal && single.complete
                && single.score == 3 && single.matches.size() == 1
                && expected.count(single.matches.front().atomMapping) == 1
                && single.matches.front().bondMapping.size() == 2
                && smsd::MCS::validateMapping(query, target, single.matches.front(), options).empty(),
                "A single star optimum must belong to the complete labelled optimum set");
            options.excludedTargetAtoms = {1};
            std::set<Mapping> excludedExpected;
            for (int leaf = 1; leaf <= 3; ++leaf)
                excludedExpected.insert(Mapping{{0, 0}, {leaf, 2}});
            const auto excludedSearch = smsd::MCS::find(query).withOptions(options);
            const auto excludedAll = excludedSearch.matchAll(target, 1000);
            require(excludedAll.optimal && excludedAll.complete && excludedAll.score == 2,
                    "Excluded target twins must retain the two-atom optimum");
            std::set<Mapping> excludedActual;
            for (const auto& match : excludedAll.matches) {
                require(smsd::MCS::validateMapping(query, target, match, options).empty(),
                        "Excluded-twin search returned an invalid mapping");
                require(excludedActual.insert(match.atomMapping).second,
                        "Excluded-twin search duplicated a mapping");
            }
            require(excludedActual == excludedExpected,
                    "Excluded-twin search omitted a labelled optimum");
            const auto excludedSingle = excludedSearch.match(target);
            require(excludedSingle.optimal && excludedSingle.complete && excludedSingle.score == 2
                && excludedSingle.matches.size() == 1
                && excludedExpected.count(excludedSingle.matches.front().atomMapping) == 1
                && smsd::MCS::validateMapping(query, target, excludedSingle.matches.front(), options).empty(),
                "An excluded class minimum must not hide a feasible twin");
        }
    }
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"oracle", exhaustiveSmallGraphs}, {"proof", connectedNaturalProductProof},
        {"adjacency", optionalAdjacencyWords}, {"word-boundaries", candidatesAcrossWordBoundaries},
        {"edge-oracle", exhaustiveCommonEdges}, {"one-fragment", oneFragmentIsConnected},
        {"full-query-proof", exhaustedFullQueryProof}, {"edge-proof", heterocycleCommonEdgeProof},
        {"twin-order", twinCandidateOrder}});
}

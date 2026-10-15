/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/mcs.hpp"
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <tuple>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

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

int main(int argc, char** argv) {
    const std::string selected = argc > 1 ? argv[1] : "all";
    int failures = 0;
    for (const auto& group : {std::pair<const char*, void(*)()>{"oracle", exhaustiveSmallGraphs},
                              {"proof", connectedNaturalProductProof}, {"adjacency", optionalAdjacencyWords},
                              {"word-boundaries", candidatesAcrossWordBoundaries}, {"edge-oracle", exhaustiveCommonEdges},
                              {"one-fragment", oneFragmentIsConnected}}) {
        if (selected != "all" && selected != group.first) continue;
        try { group.second(); std::cout << "PASS " << group.first << '\n'; }
        catch (const std::exception& error) { std::cerr << "FAIL " << group.first << ": " << error.what() << '\n'; ++failures; }
    }
    return failures == 0 ? 0 : 1;
}

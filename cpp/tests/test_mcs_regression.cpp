/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/cip.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace {

struct Graph {
    int atoms;
    std::vector<std::vector<int>> edges;
    smsd::MolGraph molecule;
};

Graph makeGraph(int atoms, const std::vector<std::pair<int, int>>& bonds,
                std::vector<int> atomicNumbers = {}) {
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    std::vector<std::vector<int>> edges(atoms, std::vector<int>(atoms));
    for (const auto& [a, b] : bonds) {
        edges[a][b] = edges[b][a] = 1;
        neighbors[a].push_back(b);
        neighbors[b].push_back(a);
        orders[a].push_back(1);
        orders[b].push_back(1);
    }
    if (atomicNumbers.empty()) atomicNumbers.assign(atoms, 6);
    auto molecule = smsd::MolGraph::Builder()
        .atomCount(atoms).atomicNumbers(atomicNumbers)
        .setNeighbors(neighbors).setBondOrders(orders).build();
    return {atoms, std::move(edges), std::move(molecule)};
}

Graph graphFromMask(int atoms, int mask) {
    std::vector<std::pair<int, int>> bonds;
    int bit = 0;
    for (int a = 0; a < atoms; ++a)
        for (int b = a + 1; b < atoms; ++b, ++bit)
            if (mask & (1 << bit)) bonds.emplace_back(a, b);
    return makeGraph(atoms, bonds);
}

bool connectedMapping(const Graph& query, const std::vector<int>& mapping) {
    int mappedCount = 0, root = -1;
    for (int atom = 0; atom < query.atoms; ++atom) {
        if (mapping[atom] < 0) continue;
        ++mappedCount;
        root = atom;
    }
    if (root < 0) return true;
    std::vector<int> seen(query.atoms), frontier{root};
    seen[root] = 1;
    for (size_t head = 0; head < frontier.size(); ++head) {
        for (int atom = 0; atom < query.atoms; ++atom) {
            if (mapping[atom] < 0 || seen[atom]
                || !query.edges[frontier[head]][atom]) continue;
            seen[atom] = 1;
            frontier.push_back(atom);
        }
    }
    return static_cast<int>(frontier.size()) == mappedCount;
}


using Objective = smsd::MCS::Objective;
using Connectivity = smsd::MCS::Connectivity;

smsd::MCS::Options options(bool induced, bool connected) {
    smsd::MCS::Options result;
    result.objective = induced ? Objective::MCIS : Objective::MCES;
    result.connectivity = connected ? Connectivity::Connected : Connectivity::Fragments;
    result.timeoutMs = 1000;
    return result;
}

/** Enumerate partial injections using only local adjacency matrices.
 * MCES selects common edges; MCIS also requires equal induced adjacency. */
double maximumScore(const Graph& query, const Graph& target, const smsd::MCS::Options& opts) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    double best = 0;
    std::function<void(int)> visit = [&](int queryAtom) {
        if (queryAtom == query.atoms) {
            std::vector<std::vector<int>> adjacency(query.atoms);
            int atoms = 0, bonds = 0;
            double weight = 0;
            for (int first = 0; first < query.atoms; ++first) if (mapping[first] >= 0) {
                ++atoms;
                if (!opts.atomWeights.empty()) weight += opts.atomWeights[first];
                for (int second = first + 1; second < query.atoms; ++second)
                    if (mapping[second] >= 0 && query.edges[first][second] && target.edges[mapping[first]][mapping[second]]) {
                        adjacency[first].push_back(second); adjacency[second].push_back(first); ++bonds;
                    }
            }
            if (opts.objective == Objective::MCES)
                for (int atom = 0; atom < query.atoms; ++atom)
                    if (mapping[atom] >= 0 && adjacency[atom].empty()) return;
            int fragments = 0;
            std::vector<bool> seen(query.atoms);
            for (int root = 0; root < query.atoms; ++root) if (mapping[root] >= 0 && !seen[root]) {
                std::vector<int> component{root}; seen[root] = true;
                for (size_t head = 0; head < component.size(); ++head)
                    for (int neighbor : adjacency[component[head]]) if (!seen[neighbor]) {
                        seen[neighbor] = true; component.push_back(neighbor);
                    }
                if (static_cast<int>(component.size()) < opts.minFragmentAtoms) return;
                ++fragments;
            }
            if (fragments > opts.maxFragments || (opts.connectivity == Connectivity::Connected && fragments > 1)) return;
            best = std::max(best, opts.objective == Objective::MCES ? double(bonds)
                : opts.atomWeights.empty() ? double(atoms) : weight);
            return;
        }
        visit(queryAtom + 1);
        for (int targetAtom = 0; targetAtom < target.atoms; ++targetAtom) {
            if (used[targetAtom] || query.molecule.atomicNum[queryAtom] != target.molecule.atomicNum[targetAtom]) continue;
            bool valid = true;
            if (opts.objective == Objective::MCIS)
                for (int previous = 0; previous < queryAtom; ++previous)
                    if (mapping[previous] >= 0 && query.edges[queryAtom][previous] != target.edges[targetAtom][mapping[previous]]) valid = false;
            if (!valid) continue;
            mapping[queryAtom] = targetAtom; used[targetAtom] = 1;
            visit(queryAtom + 1);
            mapping[queryAtom] = -1; used[targetAtom] = 0;
        }
    };
    visit(0);
    return best;
}

void validate(const Graph& query, const Graph& target, const smsd::Match& match,
              const smsd::MCS::Options& opts) {
    std::set<int> used;
    std::map<int, std::set<int>> adjacency;
    for (const auto& atom : match.atomMapping) {
        assert(atom.first >= 0 && atom.first < query.atoms && atom.second >= 0 && atom.second < target.atoms);
        assert(used.insert(atom.second).second);
        assert(query.molecule.atomicNum[atom.first] == target.molecule.atomicNum[atom.second]);
        adjacency[atom.first];
    }
    std::set<std::pair<int, int>> selected;
    for (const auto& bond : match.bondMapping) {
        assert(match.atomMapping.at(bond[0]) == bond[2] && match.atomMapping.at(bond[1]) == bond[3]);
        assert(query.edges[bond[0]][bond[1]] && target.edges[bond[2]][bond[3]]);
        assert(selected.insert(std::minmax(bond[0], bond[1])).second);
        adjacency[bond[0]].insert(bond[1]); adjacency[bond[1]].insert(bond[0]);
    }
    if (opts.objective == Objective::MCES)
        for (const auto& atom : match.atomMapping) assert(!adjacency[atom.first].empty());
    if (opts.objective == Objective::MCIS)
        for (const auto& first : match.atomMapping) for (const auto& second : match.atomMapping) if (first.first < second.first) {
            assert(query.edges[first.first][second.first] == target.edges[first.second][second.second]);
            assert(bool(query.edges[first.first][second.first]) == bool(selected.count({first.first, second.first})));
        }
    std::set<int> seen;
    std::set<std::set<int>> expected;
    for (const auto& atom : match.atomMapping) if (seen.insert(atom.first).second) {
        std::vector<int> component{atom.first};
        for (size_t head = 0; head < component.size(); ++head)
            for (int neighbor : adjacency[component[head]]) if (seen.insert(neighbor).second) component.push_back(neighbor);
        assert(static_cast<int>(component.size()) >= opts.minFragmentAtoms);
        expected.emplace(component.begin(), component.end());
    }
    std::set<std::set<int>> actual;
    for (const auto& fragment : match.fragments) actual.emplace(fragment.begin(), fragment.end());
    assert(actual == expected);
    assert(static_cast<int>(actual.size()) <= opts.maxFragments);
    assert(opts.connectivity != Connectivity::Connected || actual.size() <= 1);
}

smsd::SearchResult search(const Graph& query, const Graph& target, const smsd::MCS::Options& opts) {
    auto result = smsd::MCS::find(query.molecule).withOptions(opts).match(target.molecule);
    for (const auto& match : result.matches) validate(query, target, match, opts);
    return result;
}

void testOptionRegressions() {
    auto fragments = makeGraph(5, {{0, 1}, {2, 3}, {3, 4}});
    auto result = search(fragments, fragments, options(true, true));
    assert(result.matches.front().atomMapping.size() == 3);
    auto disconnected = options(true, false);
    assert(search(fragments, fragments, disconnected).matches.front().atomMapping.size() == 5);
    disconnected.minFragmentAtoms = 3; disconnected.maxFragments = 1;
    assert(search(fragments, fragments, disconnected).matches.front().atomMapping.size() == 3);
    auto path = graphFromMask(3, 5), triangle = graphFromMask(3, 7);
    assert(search(path, triangle, options(true, true)).matches.front().atomMapping.size() == 2);
    auto edgeResult = search(path, triangle, options(false, true));
    assert(edgeResult.matches.front().atomMapping.size() == 3 && edgeResult.score == 2);
    auto isolates = graphFromMask(3, 0), edge = graphFromMask(3, 1);
    for (const auto& pair : {std::pair<const Graph*, const Graph*>{&isolates, &edge}, {&edge, &isolates}})
        assert(search(*pair.first, *pair.second, options(true, false)).matches.front().atomMapping.size() == 2);
    auto star = graphFromMask(4, 7), edgeAndIsolates = graphFromMask(4, 1), path4 = graphFromMask(4, 13);
    assert(search(edgeAndIsolates, star, options(true, true)).matches.front().atomMapping.size() == 2);
    assert(search(edgeAndIsolates, path4, options(true, false)).matches.front().atomMapping.size() == 3);
    auto twoIsolates = graphFromMask(2, 0);
    assert(search(edge, twoIsolates, options(true, false)).matches.front().atomMapping.size() == 2);
    auto paddedEdge = makeGraph(9, {{0, 1}}), paddedStar = makeGraph(9, {{0, 1}, {0, 2}, {0, 3}});
    assert(search(paddedEdge, paddedStar, options(true, true)).matches.front().atomMapping.size() == 2);
}

void testWeightRegressions() {
    auto clique = graphFromMask(4, 63), path = makeGraph(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    auto weighted = options(true, true); weighted.atomWeights = {1, 1, 1, 10};
    auto result = search(clique, path, weighted);
    assert(result.matches.front().atomMapping.size() == 2 && result.matches.front().atomMapping.count(3));
    auto fragments = makeGraph(6, {{0, 1}, {1, 2}, {2, 3}, {4, 5}});
    weighted.atomWeights = {1, 1, 1, 1, 5, 5};
    auto heavy = search(fragments, fragments, weighted);
    assert(heavy.score == 10 && heavy.matches.front().atomMapping.count(4) && heavy.matches.front().atomMapping.count(5));
    auto bonded = makeGraph(9, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {5, 6}, {5, 7}, {5, 8}, {6, 7}, {6, 8}, {7, 8}});
    auto bonds = options(false, true);
    auto dense = search(bonded, bonded, bonds);
    assert(dense.score == 6 && dense.matches.front().atomMapping.size() == 4);
    auto reject = [&](std::vector<double> weights) {
        auto invalid = options(true, true); invalid.atomWeights = std::move(weights);
        for (bool all : {false, true}) {
            bool rejected = false;
            try {
                const auto prepared = smsd::MCS::find(clique.molecule).withOptions(invalid);
                if (all) (void)prepared.matchAll(clique.molecule); else (void)prepared.match(clique.molecule);
            } catch (const std::invalid_argument&) { rejected = true; }
            assert(rejected);
        }
    };
    reject({1});
    reject({std::numeric_limits<double>::quiet_NaN(), 0, 0, 0});
    reject({std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({-std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({std::numeric_limits<double>::max(), 0, 0, 0});
    reject({1.2e6, 1.2e6, -1.2e6, -1.2e6});
    reject({0, 0, -1.2e6, -1.2e6});
    auto atom = graphFromMask(1, 0);
    auto boundary = options(true, true); boundary.atomWeights = {(double(INT_MAX) + 0.5) / 1000.0};
    assert(search(atom, atom, boundary).score == double(INT_MAX) / 1000.0);
    boundary.atomWeights = {(double(INT_MIN) - 0.5) / 1000.0};
    assert(search(atom, atom, boundary).score == 0);
    auto signedPath = makeGraph(3, {{0, 1}, {1, 2}});
    auto signedOptions = options(true, true); signedOptions.atomWeights = {1, -5, 1};
    assert(search(signedPath, signedPath, signedOptions).score == 1);
    auto all = smsd::MCS::find(signedPath.molecule).withOptions(signedOptions).matchAll(signedPath.molecule, 10);
    assert(all.matches.size() == 6 && all.complete);
    for (const auto& match : all.matches) { assert(match.atomMapping.size() == 1); validate(signedPath, signedPath, match, signedOptions); }
    signedOptions.atomWeights = {-1, -2, -3};
    assert(search(signedPath, signedPath, signedOptions).matches.empty());
    weighted.connectivity = Connectivity::Fragments; weighted.maxFragments = 1;
    assert(search(fragments, fragments, weighted).score == 10);
    bonds.connectivity = Connectivity::Fragments; bonds.maxFragments = 1;
    assert(search(bonded, bonded, bonds).score == 6);
    auto large = makeGraph(41, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    auto triangle = graphFromMask(3, 7);
    auto singleton = options(true, true); singleton.atomWeights.assign(41, 1); singleton.atomWeights[40] = 10;
    auto selected = search(large, triangle, singleton);
    assert(selected.matches.front().atomMapping.size() == 1 && selected.matches.front().atomMapping.count(40));
}

void testObjectiveOracle() {
    size_t cases = 0;
    for (int queryMask = 0; queryMask < 64; ++queryMask) for (int targetMask = 0; targetMask < 64; ++targetMask) {
        auto query = graphFromMask(4, queryMask), target = graphFromMask(4, targetMask);
        for (bool connected : {false, true}) for (int objective = 0; objective < 4; ++objective) {
            auto opts = options(objective != 3, connected);
            if (objective == 0) opts.atomWeights = {1, 2, 5, 10};
            if (objective == 1) opts.atomWeights = {2, -5, 3, -1};
            auto result = search(query, target, opts);
            assert(result.optimal && result.score == maximumScore(query, target, opts)); ++cases;
        }
    }
    assert(cases == 32768);
    std::mt19937 random(9123);
    size_t fragmentCases = 0;
    for (int sample = 0; sample < 128; ++sample) {
        auto query = graphFromMask(4, random() % 64), target = graphFromMask(4, random() % 64);
        for (int objective = 0; objective < 4; ++objective) for (int minimum : {1, 2}) for (int maximum : {1, 2}) {
            auto opts = options(objective != 3, false); opts.minFragmentAtoms = minimum; opts.maxFragments = maximum;
            if (objective == 0) opts.atomWeights = {1, 2, 5, 10};
            if (objective == 1) opts.atomWeights = {2, -5, 3, -1};
            auto result = search(query, target, opts);
            assert(result.optimal && result.score == maximumScore(query, target, opts)); ++fragmentCases;
        }
    }
    assert(fragmentCases == 2048);
    std::cout << "Validated " << cases << " objective and " << fragmentCases << " constrained-fragment cases\n";
}

void testConstrainedBatchObjectives() {
    auto query = makeGraph(4, {{0, 1}, {1, 2}}, {6, 6, 6, 7});
    auto target = makeGraph(4, {{0, 1}, {1, 2}}, {6, 6, 6, 7});
    auto weighted = options(true, true); weighted.atomWeights = {1, 1, 1, 10};
    auto results = smsd::MCS::matchConstrained({query.molecule, query.molecule}, target.molecule, weighted);
    assert(results.size() == 2 && results[0].score == 10 && results[1].score == 3);
    assert((results[0].matches.front().atomMapping == std::map<int, int>{{3, 3}}));
    std::set<int> used;
    for (const auto& result : results) for (const auto& atom : result.matches.front().atomMapping) assert(used.insert(atom.second).second);
    auto bondQuery = makeGraph(9, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {5, 6}, {5, 7}, {5, 8}, {6, 7}, {6, 8}, {7, 8}});
    auto clique = graphFromMask(4, 63);
    auto bonds = options(false, true);
    assert(smsd::MCS::matchConstrained({bondQuery.molecule}, clique.molecule, bonds).front().score == 6);
    auto nitrogen = makeGraph(1, {}, {7}), carbon = makeGraph(3, {{0, 1}, {1, 2}});
    assert(smsd::MCS::matchConstrained({nitrogen.molecule}, carbon.molecule).front().matches.empty());
}

void testCoverageMappingValidity() {
    auto query = smsd::parseSMILES("O=C([O-])Cc1ccccc1"), target = smsd::parseSMILES("O=C([O-])CCCc1ccccc1");
    for (bool any : {false, true}) {
        auto opts = options(true, true);
        opts.chemistry.matchBondOrder = any ? smsd::ChemOptions::BondOrderMode::ANY : smsd::ChemOptions::BondOrderMode::STRICT;
        auto result = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(!result.matches.empty() && smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
    }
    auto flexibleQuery = smsd::parseSMILES("O=C(Cc1ccccc1)c1ccncc1");
    auto flexibleTarget = smsd::parseSMILES("O/C(/c1ccncc1)=C\\c1ccccc1");
    const std::map<int, int> witness{{1,11}, {2,10}, {3,9}, {4,8}, {5,1},
                                    {6,2}, {7,3}, {9,12}, {10,13}, {11,14}};
    for (int aromaticOrder : {1, 4}) {
        for (auto* graph : {&flexibleQuery, &flexibleTarget}) for (int atom = 0; atom < graph->n; ++atom)
            for (int other : graph->neighbors[atom]) if (other > atom && graph->bondAromatic(atom, other)) graph->setBondOrder(atom, other, aromaticOrder);
        auto opts = options(false, true);
        opts.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
        assert(smsd::MCS::validateMapping(flexibleQuery, flexibleTarget, witness, opts).empty());
        auto result = smsd::MCS::find(flexibleQuery).withOptions(opts).match(flexibleTarget);
        assert(result.matches.front().atomMapping.size() >= witness.size());
        assert(smsd::MCS::validateMapping(flexibleQuery, flexibleTarget, result.matches.front(), opts).empty());
    }
    auto single = smsd::parseSMILES("CO"), doubled = smsd::parseSMILES("C=O");
    auto opts = options(true, true);
    assert(smsd::MCS::find(single).withOptions(opts).match(doubled).matches.front().atomMapping.size() == 1);
    opts.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    assert(smsd::MCS::find(single).withOptions(opts).match(doubled).matches.front().atomMapping.size() == 2);
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 3; ++atoms) for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask) graphs.push_back(graphFromMask(atoms, mask));
    for (const auto& first : graphs) for (const auto& second : graphs) for (bool any : {false, true}) for (bool ring : {false, true}) {
        auto opts = options(true, true); opts.chemistry.ringMatchesRingOnly = ring;
        opts.chemistry.matchBondOrder = any ? smsd::ChemOptions::BondOrderMode::ANY : smsd::ChemOptions::BondOrderMode::STRICT;
        auto result = search(first, second, opts);
        for (const auto& match : result.matches) assert(smsd::MCS::validateMapping(first.molecule, second.molecule, match, opts).empty());
    }
}

void testNativeStages() {
    size_t cases = 0;
    for (int queryMask = 0; queryMask < 64; ++queryMask) for (int targetMask = 0; targetMask < 64; ++targetMask) {
        auto query = graphFromMask(4, queryMask), target = graphFromMask(4, targetMask);
        for (bool induced : {false, true}) {
            auto opts = options(induced, false);
            const auto prepared = smsd::MCS::find(query.molecule).withOptions(opts);
            auto first = prepared.match(target.molecule), repeated = prepared.match(target.molecule);
            assert(first.optimal && repeated.optimal && first.score == maximumScore(query, target, opts) && repeated.score == first.score);
            for (const auto& match : first.matches) validate(query, target, match, opts);
            ++cases;
        }
    }
    assert(cases == 8192);
}

void testSymmetryCanonicalization() {
    auto cycle = makeGraph(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    std::map<int, int> mapping{{0, 1}, {1, 0}, {2, 3}}, expected{{0, 0}, {1, 1}, {2, 2}};
    assert(smsd::MCS::canonicalizeMapping(cycle.molecule, cycle.molecule, mapping) == expected);
    for (int querySign : {-1, 1}) for (int queryShift = 0; queryShift < 4; ++queryShift)
        for (int targetSign : {-1, 1}) for (int targetShift = 0; targetShift < 4; ++targetShift) {
            std::map<int, int> transformed;
            for (const auto& atom : mapping) transformed[(querySign * atom.first + queryShift + 4) % 4] = (targetSign * atom.second + targetShift + 4) % 4;
            assert(smsd::MCS::canonicalizeMapping(cycle.molecule, cycle.molecule, transformed) == expected);
        }
}

void testMolecularAutomorphismProperties() {
    for (const char* smiles : {"[C+].C", "[12C].[13C]", "C=CC", "C[C@H](F)[C@H](F)C"}) {
        auto graph = smsd::parseSMILES(smiles); graph.ensureCanonical();
        assert(!graph.automorphismGeneratorsTruncated());
        std::map<int, int> mapping{{0, 0}};
        auto canonical = smsd::MCS::canonicalizeMapping(graph, graph, mapping);
        assert(smsd::MCS::canonicalizeMapping(graph, graph, canonical) == canonical);
    }
    auto capped = makeGraph(8, {});
    assert(capped.molecule.automorphismGeneratorsTruncated());
}

void testSharedDeadline() {
    std::mt19937 random(781);
    for (int atoms : {24, 32}) {
        std::vector<std::pair<int, int>> queryEdges, targetEdges;
        for (int first = 0; first < atoms; ++first) for (int second = first + 1; second < atoms; ++second) {
            if (random() % 100 < 22) queryEdges.emplace_back(first, second);
            if (random() % 100 < 22) targetEdges.emplace_back(first, second);
        }
        auto query = makeGraph(atoms, queryEdges), target = makeGraph(atoms, targetEdges);
        auto opts = options(true, true); opts.timeoutMs = 5;
        auto started = std::chrono::steady_clock::now();
        auto result = search(query, target, opts);
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - started).count();
        assert(elapsed < 100 && result.status == smsd::SearchStatus::TimedOut && !result.complete);
        auto atom = graphFromMask(1, 0);
        assert(search(atom, atom, options(true, true)).complete);
    }
    auto atom = graphFromMask(1, 0); auto invalid = options(true, true);
    invalid.atomWeights = {std::numeric_limits<double>::quiet_NaN()};
    bool rejected = false;
    try { (void)search(atom, atom, invalid); } catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected && search(atom, atom, options(true, true)).complete);
}

void testMappingValidation() {
    auto edge = graphFromMask(2, 1), isolates = graphFromMask(2, 0);
    for (int invalid : {-1, 2, INT_MAX}) {
        std::map<int, int> mapping{{0, 0}, {1, invalid}};
        assert(!smsd::MCS::validateMapping(edge.molecule, edge.molecule, mapping).empty());
    }
    std::map<int, int> full{{0, 0}, {1, 1}};
    auto edgeOptions = options(false, false), induced = options(true, false);
    assert(!smsd::MCS::validateMapping(isolates.molecule, edge.molecule, full, edgeOptions).empty());
    assert(smsd::MCS::validateMapping(edge.molecule, edge.molecule, full, edgeOptions).empty());
    assert(!smsd::MCS::validateMapping(isolates.molecule, edge.molecule, full, induced).empty());
    assert(!smsd::MCS::isMappingMaximal(isolates.molecule, edge.molecule, std::map<int, int>{{0, 0}}, edgeOptions));
    assert(smsd::MCS::isMappingMaximal(isolates.molecule, edge.molecule, std::map<int, int>{{0, 0}}, induced));
    assert(!smsd::MCS::isMappingMaximal(edge.molecule, edge.molecule, std::map<int, int>{{0, INT_MAX}}, induced));
    auto legacy = makeGraph(3, {{0, 1}}).molecule;
    legacy.formalCharge.clear(); legacy.massNumber.clear(); legacy.tetraChirality.clear();
    smsd::Substructure::Options optional; optional.chemistry.matchFormalCharge = true;
    optional.chemistry.matchIsotope = true; optional.chemistry.useChirality = true;
    auto result = smsd::Substructure::find(legacy).withOptions(optional).withTimeout(1000).match(legacy);
    assert(result.complete && result.matches.front().atomMapping.size() == 3);
    for (int atom = 0; atom < legacy.n; ++atom) assert(smsd::cip::assignRS(legacy, atom) == smsd::cip::RSLabel::NONE);
    const auto split = smsd::Substructure::smarts("C").match(legacy);
    assert(!split.matches.empty());
    const auto fragments = smsd::MCS::find(legacy).withOptions(options(true, false)).match(legacy);
    assert(fragments.matches.front().fragments.size() == 2);
    const auto graph = fragments.matches.front().toGraph(legacy);
    for (int atom = 0; atom < graph.n; ++atom) {
        assert(graph.formalCharge[atom] == 0 && graph.massNumber[atom] == 0 && graph.tetraChirality[atom] == 0);
    }
}

void testTautomerBound() {
    auto query = graphFromMask(4, 7), target = graphFromMask(4, 7);
    query.molecule.atomicNum.assign(4, 7); target.molecule.atomicNum.assign(4, 8);
    query.molecule.tautomerClass.assign(4, 0); target.molecule.tautomerClass.assign(4, 0);
    query.molecule.refreshAtomLabels(); target.molecule.refreshAtomLabels();
    auto opts = options(true, true); opts.chemistry.tautomerAware = true;
    assert(search(query, target, opts).matches.empty());
    smsd::Substructure::Options sub; sub.chemistry = opts.chemistry;
    assert(smsd::Substructure::find(query.molecule).withOptions(sub).match(target.molecule).matches.empty());
}

void testSparseConnectedSelection() {
    auto query = makeGraph(2048, {{1500, 1501}, {1700, 1701}, {1701, 1702}});
    auto target = makeGraph(2048, {{1500, 1501}, {1700, 1701}});
    // The uniquely labelled sparse component is checked without a dense product search.
    for (auto* graph : {&query.molecule, &target.molecule}) {
        graph->atomicNum[1500] = 7; graph->atomicNum[1501] = 8;
        graph->refreshAtomLabels();
    }
    auto small = makeGraph(2, {{0, 1}}, {7, 8});
    auto result = search(small, target, options(true, true));
    assert(result.matches.front().atomMapping.at(0) == 1500 && result.matches.front().atomMapping.at(1) == 1501);
    assert(result.matches.front().fragments.size() == 1);
    assert(smsd::Substructure::find(small.molecule).match(query.molecule).matches.front().atomMapping.size() == 2);
}

void testFlexibleAromaticExtension() {
    auto query = smsd::parseSMILES("c1ccc(-c2ccccc2)cc1"), target = smsd::parseSMILES("c1ccc(Cc2ccccc2)cc1");
    const std::map<int, int> witness{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {10, 11}, {11, 12}};
    auto opts = options(true, true);
    assert(smsd::MCS::validateMapping(query, target, witness, opts).empty());
    auto result = smsd::MCS::find(query).withOptions(opts).match(target);
    assert(result.matches.front().atomMapping.size() >= witness.size());
    assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
}

void testConnectedStatinSeedRecovery() {
    auto query = smsd::parseSMILES("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4");
    auto target = smsd::parseSMILES("CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C");
    const std::map<int, int> witness{{4, 3}, {5, 8}, {6, 7}, {7, 6}, {8, 5}, {19, 20}, {20, 21}, {21, 22}, {22, 23}, {23, 24}, {24, 25}, {25, 26}, {26, 9}, {27, 10}, {28, 11}, {29, 12}, {30, 13}, {32, 1}};
    auto opts = options(false, true); opts.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    assert(smsd::MCS::validateMapping(query, target, witness, opts).empty());
    const auto prepared = smsd::MCS::find(query).withOptions(opts);
    for (int repeat = 0; repeat < 2; ++repeat) {
        auto result = prepared.match(target);
        assert(result.matches.front().atomMapping.size() >= witness.size());
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        assert(result.matches.front().fragments.size() == 1);
    }
}

void testMcGregorStateIsolation() {
    size_t cases = 0;
    for (int queryMask = 0; queryMask < 64; ++queryMask) for (int targetMask : {0, 7, 18, 23, 32, 41, 56, 63}) {
        auto query = graphFromMask(4, queryMask), core = graphFromMask(4, targetMask);
        std::vector<std::pair<int, int>> bonds;
        for (int first = 0; first < 4; ++first) for (int second = first + 1; second < 4; ++second)
            if (core.edges[first][second]) bonds.emplace_back(first, second);
        auto target = makeGraph(20, bonds);
        for (bool connected : {false, true}) {
            auto opts = options(true, connected); opts.timeoutMs = 100;
            const auto prepared = smsd::MCS::find(query.molecule).withOptions(opts);
            auto result = prepared.match(target.molecule);
            for (const auto& match : result.matches) validate(query, target, match, opts);
            if (connected) assert(result.score <= maximumScore(query, core, opts));
            assert(!result.matches.empty());
            assert(prepared.match(query.molecule).score >= result.score);
            ++cases;
        }
    }
    assert(cases == 1024);
}

void testLargeConnectedCoverageSeed() {
    auto query = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C");
    auto target = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O");
    assert(std::min(query.n, target.n) > 50);
    auto opts = options(false, true); opts.timeoutMs = 2000;
    opts.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    auto result = smsd::MCS::find(query).withOptions(opts).match(target);
    assert(result.matches.front().atomMapping.size() >= 50);
    assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
    assert(result.matches.front().fragments.size() == 1);
}

void testReverseCompleteRingFilter() {
    auto first = smsd::parseSMILES("NNc1ccc2ccccc2c1"), second = smsd::parseSMILES("NNc1ccccc1");
    for (bool induced : {false, true}) for (bool reverse : {false, true}) {
        const auto& query = reverse ? second : first;
        const auto& target = reverse ? first : second;
        auto opts = options(induced, true); opts.chemistry.completeRingsOnly = true;
        auto result = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(result.matches.front().atomMapping.size() == 8 && result.matches.front().bondMapping.size() == 8);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        const auto graph = result.matches.front().toGraph(query);
        assert(graph.computeRelevantCycles().size() == 1 && graph.computeRelevantCycles().front().size() == 6);
    }
}

void testExhaustiveSmallGraphs() {
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 4; ++atoms) for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask) graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) for (const auto& target : graphs)
        for (bool induced : {false, true}) for (bool connected : {false, true}) {
            auto opts = options(induced, connected); auto result = search(query, target, opts);
            assert(result.optimal && result.score == maximumScore(query, target, opts)); ++cases;
        }
    assert(cases == 23104);
    std::cout << "Validated " << cases << " exhaustive small-graph MCS cases\n";
}

void testEnumeratedMappings() {
    auto isolates = graphFromMask(2, 0), path = graphFromMask(3, 5);
    auto result = smsd::MCS::find(isolates.molecule).withOptions(options(true, false)).matchAll(path.molecule, 10);
    assert(!result.matches.empty());
    for (const auto& match : result.matches) { assert(match.atomMapping.size() == 2); validate(isolates, path, match, options(true, false)); }
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 3; ++atoms) for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask) graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) for (const auto& target : graphs)
        for (bool induced : {false, true}) for (bool connected : {false, true}) {
            auto opts = options(induced, connected);
            auto all = smsd::MCS::find(query.molecule).withOptions(opts).matchAll(target.molecule, 3);
            double expected = maximumScore(query, target, opts);
            assert(all.optimal && all.score == expected && (expected == 0 || !all.matches.empty()));
            for (const auto& match : all.matches) validate(query, target, match, opts);
            ++cases;
        }
    assert(cases == 576);
}

smsd::MCSOptions commonBondOptions(smsd::MCSConnectivity connectivity, smsd::ChemOptions::BondOrderMode bondOrder) {
    smsd::MCSOptions result;
    result.objective = smsd::MCSObjective::MCES;
    result.connectivity = connectivity;
    result.chemistry.matchBondOrder = bondOrder;
    result.timeoutMs = 30000;
    return result;
}

std::string mappingText(const smsd::Match& match) {
    std::string text = "{";
    for (const auto& atom : match.atomMapping)
        text += (text.size() > 1 ? ", " : "") + std::to_string(atom.first) + "=" + std::to_string(atom.second);
    return text + "}";
}

void testCommonBondHardFragmentProofs() {
    // ChEMBL random pairs 3471, 7381 and 3168 with their optimal common bond counts.
    const std::vector<std::tuple<const char*, const char*, int>> pairs{
        {"N#Cc1ccc(C(Nc2ccc3cc4ccc(NC(c5ccc(C#N)cc5)=O)cc4nc3c2)=O)cc1",
         "OCCCNP(O)(=O)OC1CC(n2cnc3c(N)ncnc23)OC1COP(O)(=O)OC1CC(n2ccc(N)nc2=O)OC1COP(O)(O)=O", 24},
        {"COc1c(Cl)cc(Nc2nc3ccc(C(=O)c4ccccc4)cc3n3c(C)nnc32)c(OC)c1",
         "CCCCc1n(Cc2ccc(-c3ccccc3-c3nn[nH]n3)cc2)c(SCc2ccccc2-c2nn[nH]n2)nn1", 27},
        {"CCCCCCCCCCCCCCCNC(OCCOCCOC(CCCCC[n+]1ccsc1)=O)=O",
         "OP(O)(=O)C(F)(F)c1ccc(CC(n2nnc3ccccc23)(c2scc(-c3ccccc3)n2)Cc2ccc(C(F)(F)P(O)(O)=O)cc2)cc1", 27}};
    for (const auto& [querySmiles, targetSmiles, bonds] : pairs) {
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        auto opts = commonBondOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::ANY);
        const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(result.score == bonds && result.matches.front().bondMapping.size() == static_cast<size_t>(bonds));
        assert(result.nodesVisited < 20000);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        // The node limit caps the pairs the search maps: with as many as the proof took, the
        // search proves the same optimum; with fewer, it stops with a valid mapping and no proof.
        opts.nodeLimit = result.nodesVisited;
        const auto enough = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(enough.status == smsd::SearchStatus::Completed && enough.optimal);
        assert(enough.matches.front().bondMapping.size() == static_cast<size_t>(bonds));
        opts.nodeLimit = std::max<uint64_t>(1, result.nodesVisited / 2);
        const auto capped = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(capped.status == smsd::SearchStatus::LimitReached && !capped.optimal);
        assert(capped.nodesVisited <= opts.nodeLimit);
        assert(capped.matches.front().bondMapping.size() <= static_cast<size_t>(bonds));
        assert(smsd::MCS::validateMapping(query, target, capped.matches.front(), opts).empty());
    }
}

void testCommonBondStoppedSearches() {
    // ChEMBL random pair 8031, which takes far longer than these limits to prove.
    const auto query = smsd::parseSMILES("Cn1nnnc1Sc1ncnc2scc(-c3c(Cl)cccc3)c21");
    const auto target = smsd::parseSMILES("NCCCCC(N)C(=O)NC(CCCCN)C(NC(CCCCN)C(NC(C(NC(C(NC(C(NC(C(NC(C=O)CCCCN)=O)"
                                          "CCCCN)=O)CCCCN)=O)CCCCN)=O)CCCCN)=O)=O");
    auto opts = commonBondOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
    opts.timeoutMs = 200;
    const auto timedOut = smsd::MCS::find(query).withOptions(opts).match(target);
    opts.timeoutMs = 30000;
    opts.cancellation = std::make_shared<smsd::CancellationToken>();
    std::thread canceller([token = opts.cancellation] {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        token->cancel();
    });
    const auto cancelled = smsd::MCS::find(query).withOptions(opts).match(target);
    canceller.join();
    assert(timedOut.status == smsd::SearchStatus::TimedOut && cancelled.status == smsd::SearchStatus::Cancelled);
    assert(cancelled.elapsedMillis < 10000);
    for (const auto* result : {&timedOut, &cancelled}) {
        assert(!result->optimal && !result->complete && !result->matches.empty());
        assert(result->score == static_cast<double>(result->matches.front().bondMapping.size()));
        assert(smsd::MCS::validateMapping(query, target, result->matches.front(), opts).empty());
    }
}

void testCommonBondScoresIgnoreUnreachableNodeLimits() {
    const std::vector<std::pair<const char*, const char*>> pairs{
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"},
        {"CC(C)Cc1ccc(C(C)C(=O)O)cc1", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"},
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"},
        {"c1ccc2ccccc2c1", "c1ccccc1CC"},
        {"C1CCC2CCCCC2C1", "C12C3C4C1C5C2C3C45"},
        {"NCC(=O)O.CCO", "OCC(N)=O"}};
    for (const auto& [querySmiles, targetSmiles] : pairs) {
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (auto bondOrder : {smsd::ChemOptions::BondOrderMode::STRICT, smsd::ChemOptions::BondOrderMode::ANY})
                for (bool exclude : {false, true}) {
                    auto opts = commonBondOptions(connectivity, bondOrder);
                    if (exclude) opts.excludedTargetAtoms = {0, target.n - 1};
                    const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
                    // An unreachable node limit leaves the search and its objective unchanged.
                    auto generalOptions = opts;
                    generalOptions.nodeLimit = std::numeric_limits<uint64_t>::max();
                    const auto general = smsd::MCS::find(query).withOptions(generalOptions).match(target);
                    assert(result.optimal && general.optimal && result.score == general.score);
                    assert(result.matches.empty() == general.matches.empty());
                    if (result.matches.empty()) continue;
                    const auto& match = result.matches.front();
                    assert(match.bondMapping.size() == general.matches.front().bondMapping.size());
                    assert(match.atomMapping.size() == general.matches.front().atomMapping.size());
                    assert(smsd::MCS::validateMapping(query, target, match, opts).empty());
                    for (const auto& atom : match.atomMapping)
                        assert(!exclude || (atom.second != 0 && atom.second != target.n - 1));
                }
    }
}

// Every MCES option goes to the common-bond search, which gives the same tied mapping as Java.
void testCommonBondOptionsKeepResults() {
    struct Case {
        const char* mode;
        smsd::MCSConnectivity connectivity;
        double score;
        size_t matches;
        const char* firstMapping;
    };
    const auto connected = smsd::MCSConnectivity::Connected, fragments = smsd::MCSConnectivity::Fragments;
    const std::vector<Case> cases{
        {"all", connected, 7, 2, "{3=8, 4=7, 5=6, 6=5, 7=4, 8=10, 9=9}"},
        {"all", fragments, 9, 2, "{0=0, 1=1, 2=2, 3=8, 4=7, 5=6, 6=5, 7=4, 8=10, 9=9}"},
        {"chirality", connected, 2, 1, "{3=1, 4=2, 5=3}"},
        {"chirality", fragments, 2, 1, "{3=1, 4=2, 5=3}"},
        {"bondStereo", connected, 3, 1, "{0=5, 1=4, 2=3, 3=2}"},
        {"bondStereo", fragments, 4, 1, "{0=0, 1=1, 2=2, 3=3, 4=4, 5=5}"},
        {"rings", connected, 6, 1, "{3=0, 4=1, 5=2, 6=3, 7=4, 8=5}"},
        {"rings", fragments, 6, 1, "{3=0, 4=1, 5=2, 6=3, 7=4, 8=5}"},
        {"weights", connected, 9, 1, "{3=1, 4=2}"},
        {"weights", fragments, 14, 1, "{1=3, 2=4, 3=1, 4=2}"},
        {"minAtoms", connected, 3, 1, "{0=0, 1=1, 2=2, 3=3}"},
        {"minAtoms", fragments, 5, 1, "{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}"}};
    for (const auto& test : cases) {
        const std::string mode = test.mode;
        const char* querySmiles = mode == "all" ? "CC(=O)Oc1ccccc1C(=O)O" : mode == "chirality" ? "CCC[C@H](F)Cl"
            : mode == "bondStereo" ? "CC/C=C/CC" : mode == "rings" ? "c1ccc2ccccc2c1" : mode == "weights" ? "CCO.CN" : "CCCC.OCCO";
        const char* targetSmiles = mode == "all" ? "CC(=O)Nc1ccc(O)cc1" : mode == "chirality" ? "C[C@@H](F)Cl"
            : mode == "bondStereo" ? "CC/C=C\\CC" : mode == "rings" ? "c1ccccc1CC" : mode == "weights" ? "CCN.CO" : "CCCCC.OCCCO";
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        auto opts = commonBondOptions(test.connectivity, smsd::ChemOptions::BondOrderMode::STRICT);
        opts.chemistry.useChirality = mode == "chirality";
        opts.chemistry.useBondStereo = mode == "bondStereo";
        opts.chemistry.completeRingsOnly = mode == "rings";
        if (mode == "weights") {
            opts.objective = smsd::MCSObjective::MCIS;
            opts.atomWeights = {1, 2, 3, 4, 5};
        }
        if (mode == "minAtoms") opts.minFragmentAtoms = 3;
        const auto prepared = smsd::MCS::find(query).withOptions(opts);
        const auto result = mode == "all" ? prepared.matchAll(target, 1000) : prepared.match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal);
        assert(result.score == test.score && result.matches.size() == test.matches);
        assert(mappingText(result.matches.front()) == test.firstMapping);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
    }
}

/** Draws from a fixed 64-bit linear congruential sequence, the same as in the Java tests. */
struct Draws {
    uint64_t state;
    int next(int bound) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<int>((state >> 33) % static_cast<uint64_t>(bound));
    }
};

void testCommonBondSeededOptionsIgnoreUnreachableNodeLimits() {
    // The same pairs and draws are in the Java tests.
    const std::vector<std::pair<const char*, const char*>> pairs{
        {"Oc1ccccc1", "Oc1ccncc1"}, {"CC(=O)[O-]", "CC(=O)OC"}, {"[13CH3]C(=O)O", "CC(=O)O"},
        {"c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"}, {"O=c1cc[nH]cc1", "Oc1ccncc1"},
        {"C1CC2CCC1C2", "CC1CCCCC1"}, {"C[N+](C)(C)CC(=O)O", "CN(C)CC(N)=O"},
        {"NCC(=O)O.CCO", "OCC(N)=O"}, {"c1ccc(cc1)-c1ccccn1", "C=Cc1ccccc1"},
        {"[15NH2]c1ccccc1", "Nc1ccc(O)cc1"}};
    using Chem = smsd::ChemOptions;
    const Chem::BondOrderMode bondOrders[] = {Chem::BondOrderMode::LOOSE, Chem::BondOrderMode::ANY};
    const Chem::AromaticityMode aromaticities[] = {Chem::AromaticityMode::STRICT, Chem::AromaticityMode::FLEXIBLE};
    const Chem::RingFusionMode fusions[] = {Chem::RingFusionMode::IGNORE, Chem::RingFusionMode::PERMISSIVE,
                                            Chem::RingFusionMode::STRICT};
    Draws draws{20261019U};
    int cases = 0;
    for (const auto& [querySmiles, targetSmiles] : pairs)
        for (int round = 0; round < 8; ++round, ++cases) {
            // Fresh graphs, so that validation prepares what the options compare by itself.
            const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
            const auto connectivity = draws.next(2) == 0 ? smsd::MCSConnectivity::Connected
                                                         : smsd::MCSConnectivity::Fragments;
            auto opts = commonBondOptions(connectivity, bondOrders[draws.next(2)]);
            opts.chemistry.aromaticityMode = aromaticities[draws.next(2)];
            opts.chemistry.matchIsotope = draws.next(2) == 1;
            opts.chemistry.matchFormalCharge = draws.next(2) == 1;
            opts.chemistry.ringFusionMode = fusions[draws.next(3)];
            opts.chemistry.matchAtomType = draws.next(4) != 0;
            for (int excluded = draws.next(3); excluded > 0; --excluded) {
                const int atom = draws.next(target.n);
                if (std::find(opts.excludedTargetAtoms.begin(), opts.excludedTargetAtoms.end(), atom)
                    == opts.excludedTargetAtoms.end())
                    opts.excludedTargetAtoms.push_back(atom);
            }
            const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
            // An unreachable node limit leaves the search and its objective unchanged.
            auto generalOptions = opts;
            generalOptions.nodeLimit = std::numeric_limits<uint64_t>::max();
            const auto general = smsd::MCS::find(query).withOptions(generalOptions).match(target);
            assert(result.optimal && general.optimal && result.score == general.score);
            assert(result.matches.empty() == general.matches.empty());
            if (result.matches.empty()) continue;
            for (const auto* found : {&result, &general}) {
                const auto& match = found->matches.front();
                assert(smsd::MCS::validateMapping(query, target, match, opts).empty());
                for (const auto& atom : match.atomMapping)
                    assert(std::find(opts.excludedTargetAtoms.begin(), opts.excludedTargetAtoms.end(), atom.second)
                           == opts.excludedTargetAtoms.end());
            }
            assert(result.matches.front().bondMapping.size() == general.matches.front().bondMapping.size());
            assert(result.matches.front().atomMapping.size() == general.matches.front().atomMapping.size());
        }
    assert(cases == 80);
}

/** Three aromatic carbons joined by bonds of the given orders that are not aromatic. */
smsd::MolGraph aromaticChain(int firstOrder, int secondOrder) {
    return smsd::MolGraph::Builder().atomCount(3).atomicNumbers({6, 6, 6}).aromaticFlags({1, 1, 1})
        .setNeighbors({{1}, {0, 2}, {1}}).setBondOrders({{firstOrder}, {firstOrder, secondOrder}, {secondOrder}})
        .build(false);
}

void testFlexibleAromaticityResonanceBonds() {
    // The same graphs and expectations are in the Java tests.
    const auto query = aromaticChain(1, 1), target = aromaticChain(2, 1);
    for (auto mode : {smsd::ChemOptions::AromaticityMode::STRICT, smsd::ChemOptions::AromaticityMode::FLEXIBLE}) {
        const bool flexible = mode == smsd::ChemOptions::AromaticityMode::FLEXIBLE;
        smsd::Substructure::Options substructure;
        substructure.chemistry.aromaticityMode = mode;
        assert(smsd::Substructure::find(query).withOptions(substructure).match(target).matches.empty() != flexible);
        for (uint64_t nodeLimit : {uint64_t(0), std::numeric_limits<uint64_t>::max()}) {
            auto opts = commonBondOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
            opts.chemistry.aromaticityMode = mode;
            opts.nodeLimit = nodeLimit;
            const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
            assert(result.optimal && !result.matches.empty());
            assert(result.matches.front().bondMapping.size() == (flexible ? 2U : 1U));
            assert(result.matches.front().atomMapping.size() == (flexible ? 3U : 2U));
            assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        }
    }
}

smsd::MCSOptions commonAtomOptions(smsd::MCSConnectivity connectivity, smsd::ChemOptions::BondOrderMode bondOrder) {
    smsd::MCSOptions result;
    result.objective = smsd::MCSObjective::MCIS;
    result.connectivity = connectivity;
    result.chemistry.matchBondOrder = bondOrder;
    result.timeoutMs = 30000;
    return result;
}

void testCommonAtomHardFragmentProofs() {
    // ChEMBL random pairs 5776, 7041 and 8361 with their optimal atoms and bonds. The same
    // mappings are in the Java tests.
    const std::vector<std::tuple<const char*, const char*, int, int, const char*>> pairs{
        {"CC(OC1C(=O)C2(C)C(O)CC3OCC3(OC(C)=O)C2C(OC(c2cc(N=[N+]=[N-])ccc2)=O)C2(O)CC(OC(C(O)C(NC(OC(C)(C)C)=O)C(F)F)=O)C(C)=C1C2(C)C)=O",
         "CC1N(C(=O)c2ccc(C(F)(F)F)cc2)CC(=O)OC1(Cn1cncn1)c1c(F)cc(F)cc1", 26, 23,
         "{0=15, 1=16, 2=18, 3=19, 4=1, 6=0, 22=3, 23=5, 24=14, 25=13, 29=8, 30=7, 31=6, 32=4, 33=32, "
         "35=30, 36=29, 41=22, 42=23, 43=24, 51=11, 52=10, 54=27, 56=26, 57=33, 60=17}"},
        {"COc1ccc(CCNC(=O)Cn2cc(S(Cc3ccc(F)cc3)(=O)=O)c3ccccc23)cc1OC",
         "CC(OC1OC2OC3(C)CCC4C(C)CCC(C1C)C24OO3)c1ccc(C(=O)O)cc1", 25, 20,
         "{0=5, 1=6, 2=7, 3=9, 4=10, 5=11, 6=12, 7=13, 10=27, 13=0, 14=1, 16=3, 17=17, 18=18, 22=15, "
         "23=16, 24=28, 25=20, 26=22, 27=23, 28=24, 29=25, 30=29, 31=30, 33=8}"},
        {"C=CCNC(=O)C(C#N)=Cc1c(=O)n2cccc(C)c2nc1N1CCN(c2ccc(OC)cc2)CC1",
         "CS(=O)c1c(-c2nc(C(=O)Nc3ccccc3)c(N)nc2)cccc1", 21, 18,
         "{5=2, 10=7, 11=8, 12=9, 13=10, 14=11, 15=12, 16=13, 17=14, 18=15, 20=18, 21=17, 22=19, "
         "24=0, 26=4, 27=21, 28=22, 29=23, 32=24, 33=3, 35=20}"}};
    for (const auto& [querySmiles, targetSmiles, atoms, bonds, mapping] : pairs) {
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::ANY);
        const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(result.matches.front().atomMapping.size() == static_cast<size_t>(atoms));
        assert(result.matches.front().bondMapping.size() == static_cast<size_t>(bonds));
        assert(result.score == atoms && mappingText(result.matches.front()) == mapping);
        assert(result.nodesVisited < 20000);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        // The finite vertex search may finish before using its node allowance.
        opts.nodeLimit = 20000;
        const auto general = smsd::MCS::find(query).withOptions(opts).match(target);
        std::cerr << "Capped common-atom result: expected_atoms=" << atoms
                  << ", expected_bonds=" << bonds << ", status=" << static_cast<int>(general.status)
                  << ", score=" << general.score << ", optimal=" << general.optimal
                  << ", complete=" << general.complete << ", nodes=" << general.nodesVisited
                  << ", atoms=" << (general.matches.empty() ? 0 : general.matches.front().atomMapping.size())
                  << ", bonds=" << (general.matches.empty() ? 0 : general.matches.front().bondMapping.size()) << '\n';
        assert(!general.matches.empty());
        assert(general.nodesVisited <= opts.nodeLimit);
        assert(general.status == smsd::SearchStatus::Completed
               || general.status == smsd::SearchStatus::LimitReached);
        assert(general.score == static_cast<double>(general.matches.front().atomMapping.size()));
        assert(smsd::MCS::validateMapping(query, target, general.matches.front(), opts).empty());
        assert(general.matches.front().atomMapping.size() <= static_cast<size_t>(atoms));
        if (general.matches.front().atomMapping.size() == static_cast<size_t>(atoms))
            assert(general.matches.front().bondMapping.size() <= static_cast<size_t>(bonds));
        if (general.status == smsd::SearchStatus::Completed)
            assert(general.optimal && general.complete);
        else
            assert(!general.complete);
        if (general.optimal) {
            assert(general.score == atoms);
            assert(general.matches.front().atomMapping.size() == static_cast<size_t>(atoms));
            assert(general.matches.front().bondMapping.size() == static_cast<size_t>(bonds));
        } else {
            assert(general.status == smsd::SearchStatus::LimitReached);
        }
    }
}

void testCommonAtomStoppedSearches() {
    // ChEMBL random pair 7175, which takes far longer than these limits to prove.
    const auto query = smsd::parseSMILES("COc1cc2c(c(OC)c1OC)-c1c(cc(=O)c(SC)cc1)C(NC(=O)c1ccc([N+]([O-])=O)"
                                         "c(CO[N+]([O-])=O)c1)CC2");
    const auto target = smsd::parseSMILES("CC(C)=CCc1cc(C2CC(=O)c3c(c(CC=C(C)C)c(O)c(CC=C(C)C)c3O)O2)cc(O)c1O");
    auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
    opts.timeoutMs = 200;
    const auto timedOut = smsd::MCS::find(query).withOptions(opts).match(target);
    opts.timeoutMs = 30000;
    opts.cancellation = std::make_shared<smsd::CancellationToken>();
    std::thread canceller([token = opts.cancellation] {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        token->cancel();
    });
    const auto cancelled = smsd::MCS::find(query).withOptions(opts).match(target);
    canceller.join();
    assert(timedOut.status == smsd::SearchStatus::TimedOut && cancelled.status == smsd::SearchStatus::Cancelled);
    assert(cancelled.elapsedMillis < 10000);
    for (const auto* result : {&timedOut, &cancelled}) {
        assert(!result->optimal && !result->complete && !result->matches.empty());
        assert(result->score == static_cast<double>(result->matches.front().atomMapping.size()));
        assert(smsd::MCS::validateMapping(query, target, result->matches.front(), opts).empty());
    }
}

void testCommonAtomScoresAgreeWithGeneralSearch() {
    const std::vector<std::pair<const char*, const char*>> pairs{
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"},
        {"CC(C)Cc1ccc(C)cc1", "COc1ccc(C(C)C=O)cc1"},
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"},
        {"c1ccc2ccccc2c1", "c1ccccc1CC"},
        {"C1CCC2CCCCC2C1", "C12C3C4C1C5C2C3C45"},
        {"NCC(=O)O.CCO", "OCC(N)=O"},
        {"CC(C)(C)c1ccccc1.O.N", "O.Nc1ccccc1"},
        {"C1CC1C1CC1.C", "CC.C1CC1"}};
    for (const auto& [querySmiles, targetSmiles] : pairs) {
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        for (int minAtoms = 1; minAtoms <= 2; ++minAtoms)
            for (auto bondOrder : {smsd::ChemOptions::BondOrderMode::STRICT, smsd::ChemOptions::BondOrderMode::ANY})
                for (bool exclude : {false, true}) {
                    auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, bondOrder);
                    opts.minFragmentAtoms = minAtoms;
                    if (exclude) opts.excludedTargetAtoms = {0, target.n - 1};
                    const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
                    // An unreachable node limit keeps the general search without changing the objective.
                    auto generalOptions = opts;
                    generalOptions.nodeLimit = std::numeric_limits<uint64_t>::max();
                    const auto general = smsd::MCS::find(query).withOptions(generalOptions).match(target);
                    assert(result.optimal && general.optimal && result.score == general.score);
                    assert(result.matches.empty() == general.matches.empty());
                    if (result.matches.empty()) continue;
                    const auto& match = result.matches.front();
                    assert(match.atomMapping.size() == general.matches.front().atomMapping.size());
                    assert(match.bondMapping.size() == general.matches.front().bondMapping.size());
                    assert(smsd::MCS::validateMapping(query, target, match, opts).empty());
                    for (const auto& atom : match.atomMapping)
                        assert(!exclude || (atom.second != 0 && atom.second != target.n - 1));
                }
    }
}

void testCommonAtomIneligibleOptionsKeepResults() {
    struct Case {
        const char* mode;
        smsd::MCSConnectivity connectivity;
        double score;
        size_t matches;
        const char* firstMapping;
    };
    const auto connected = smsd::MCSConnectivity::Connected, fragments = smsd::MCSConnectivity::Fragments;
    const std::vector<Case> cases{
        {"all", fragments, 9, 12, "{0=0, 1=1, 2=2, 4=4, 5=5, 6=6, 7=7, 8=9, 9=10}"},
        {"chirality", fragments, 3, 1, "{3=1, 4=2, 5=3}"},
        {"bondStereo", fragments, 5, 1, "{0=0, 1=1, 2=2, 3=3, 5=5}"},
        {"rings", fragments, 6, 1, "{0=2, 1=1, 2=0, 3=5, 8=4, 9=3}"},
        {"weights", fragments, 14, 1, "{1=3, 2=4, 3=1, 4=2}"},
        {"minAtoms", fragments, 7, 1, "{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}"},
        {"maxFragments", fragments, 7, 1, "{0=0, 1=1, 2=2, 3=3, 4=5, 5=6, 6=7}"},
        {"plain", connected, 7, 1, "{3=8, 4=7, 5=9, 6=10, 7=4, 8=5, 9=6}"}};
    for (const auto& test : cases) {
        const std::string mode = test.mode;
        const bool drugs = mode == "all" || mode == "plain";
        const char* querySmiles = drugs ? "CC(=O)Oc1ccccc1C(=O)O" : mode == "chirality" ? "CCC[C@H](F)Cl"
            : mode == "bondStereo" ? "CC/C=C/CC" : mode == "rings" ? "c1ccc2ccccc2c1" : mode == "weights" ? "CCO.CN" : "CCCC.OCCO";
        const char* targetSmiles = drugs ? "CC(=O)Nc1ccc(O)cc1" : mode == "chirality" ? "C[C@@H](F)Cl"
            : mode == "bondStereo" ? "CC/C=C\\CC" : mode == "rings" ? "c1ccccc1CC" : mode == "weights" ? "CCN.CO" : "CCCCC.OCCCO";
        const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
        auto opts = commonAtomOptions(test.connectivity, smsd::ChemOptions::BondOrderMode::STRICT);
        opts.chemistry.useChirality = mode == "chirality";
        opts.chemistry.useBondStereo = mode == "bondStereo";
        opts.chemistry.completeRingsOnly = mode == "rings";
        if (mode == "weights") opts.atomWeights = {1, 2, 3, 4, 5};
        if (mode == "minAtoms") opts.minFragmentAtoms = 3;
        if (mode == "maxFragments") opts.maxFragments = 2;
        const auto prepared = smsd::MCS::find(query).withOptions(opts);
        const auto result = mode == "all" ? prepared.matchAll(target, 1000) : prepared.match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal);
        assert(result.score == test.score && result.matches.size() == test.matches);
        assert(mappingText(result.matches.front()) == test.firstMapping);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
    }
}

void testCommonAtomSeededOptionsAgreeWithGeneralSearch() {
    // The same pairs and draws are in the Java tests.
    const std::vector<std::pair<const char*, const char*>> pairs{
        {"Oc1ccccc1", "Oc1ccncc1"}, {"CC(=O)[O-]", "CC(=O)OC"}, {"[13CH3]C(=O)O", "CC(=O)O"},
        {"c1ccc2ccccc2c1", "C1CCC2CCCCC2C1"}, {"O=c1cc[nH]cc1", "Oc1ccncc1"},
        {"C1CC2CCC1C2", "CC1CCCCC1"}, {"C[N+](C)(C)CC(=O)O", "CN(C)CC(N)=O"},
        {"NCC(=O)O.CCO", "OCC(N)=O"}, {"c1ccc(cc1)-c1ccccn1", "C=Cc1ccccc1"},
        {"[15NH2]c1ccccc1", "Nc1ccc(O)cc1"}};
    using Chem = smsd::ChemOptions;
    const Chem::BondOrderMode bondOrders[] = {Chem::BondOrderMode::STRICT, Chem::BondOrderMode::LOOSE,
                                              Chem::BondOrderMode::ANY};
    const Chem::AromaticityMode aromaticities[] = {Chem::AromaticityMode::STRICT, Chem::AromaticityMode::FLEXIBLE};
    const Chem::RingFusionMode fusions[] = {Chem::RingFusionMode::IGNORE, Chem::RingFusionMode::PERMISSIVE,
                                            Chem::RingFusionMode::STRICT};
    Draws draws{20261020U};
    int cases = 0;
    for (const auto& [querySmiles, targetSmiles] : pairs)
        for (int round = 0; round < 8; ++round, ++cases) {
            // Fresh graphs, so that validation prepares what the options compare by itself.
            const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
            auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, bondOrders[draws.next(3)]);
            opts.minFragmentAtoms = 1 + draws.next(2);
            opts.chemistry.aromaticityMode = aromaticities[draws.next(2)];
            opts.chemistry.matchIsotope = draws.next(2) == 1;
            opts.chemistry.matchFormalCharge = draws.next(2) == 1;
            opts.chemistry.ringFusionMode = fusions[draws.next(3)];
            opts.chemistry.matchAtomType = draws.next(4) != 0;
            for (int excluded = draws.next(3); excluded > 0; --excluded) {
                const int atom = draws.next(target.n);
                if (std::find(opts.excludedTargetAtoms.begin(), opts.excludedTargetAtoms.end(), atom)
                    == opts.excludedTargetAtoms.end())
                    opts.excludedTargetAtoms.push_back(atom);
            }
            const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
            // An unreachable node limit keeps the general search without changing the objective.
            auto generalOptions = opts;
            generalOptions.nodeLimit = std::numeric_limits<uint64_t>::max();
            const auto general = smsd::MCS::find(query).withOptions(generalOptions).match(target);
            assert(result.optimal && general.optimal && result.score == general.score);
            assert(result.matches.empty() == general.matches.empty());
            if (result.matches.empty()) continue;
            for (const auto* found : {&result, &general}) {
                const auto& match = found->matches.front();
                assert(smsd::MCS::validateMapping(query, target, match, opts).empty());
                for (const auto& atom : match.atomMapping)
                    assert(std::find(opts.excludedTargetAtoms.begin(), opts.excludedTargetAtoms.end(), atom.second)
                           == opts.excludedTargetAtoms.end());
            }
            assert(result.matches.front().atomMapping.size() == general.matches.front().atomMapping.size());
            assert(result.matches.front().bondMapping.size() == general.matches.front().bondMapping.size());
        }
    assert(cases == 80);
}

void testCommonAtomResonanceBonds() {
    // The same graphs and expectations are in the Java tests.
    const auto query = aromaticChain(1, 1), target = aromaticChain(2, 1);
    for (auto mode : {smsd::ChemOptions::AromaticityMode::STRICT, smsd::ChemOptions::AromaticityMode::FLEXIBLE}) {
        const bool flexible = mode == smsd::ChemOptions::AromaticityMode::FLEXIBLE;
        for (uint64_t nodeLimit : {uint64_t(0), std::numeric_limits<uint64_t>::max()}) {
            auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
            opts.chemistry.aromaticityMode = mode;
            opts.nodeLimit = nodeLimit;
            const auto result = smsd::MCS::find(query).withOptions(opts).match(target);
            assert(result.optimal && !result.matches.empty());
            assert(result.matches.front().atomMapping.size() == (flexible ? 3U : 2U));
            assert(result.matches.front().bondMapping.size() == (flexible ? 2U : 1U));
            assert(smsd::MCS::validateMapping(query, target, result.matches.front(), opts).empty());
        }
    }
}


smsd::MolGraph isolatedAtoms(int count, int element) {
    return smsd::MolGraph::Builder().atomCount(count).atomicNumbers(std::vector<int>(count, element))
        .setNeighbors(std::vector<std::vector<int>>(count)).setBondOrders(std::vector<std::vector<int>>(count)).build(false);
}

/** A graph of carbons in a cycle, each bonded to the next `step` atoms around it. */
smsd::MolGraph circulant(int atoms, int step) {
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    for (int atom = 0; atom < atoms; ++atom)
        for (int offset = 1; offset <= step; ++offset) {
            const int other = (atom + offset) % atoms;
            neighbors[atom].push_back(other);
            neighbors[other].push_back(atom);
            orders[atom].push_back(1);
            orders[other].push_back(1);
        }
    return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 6))
        .setNeighbors(neighbors).setBondOrders(orders).build(false);
}

/** A carbon bonded to `leaves` atoms: an oxygen first when `oxygen`, and carbons otherwise. */
smsd::MolGraph star(int leaves, bool oxygen) {
    std::string smiles = oxygen ? "C(O)" : "C(C)";
    for (int leaf = 2; leaf < leaves; ++leaf) smiles += "(C)";
    return smsd::parseSMILES(smiles + "C");
}

void testVertexSearchHonoursTimeoutsOnSparseGraphs() {
    // The first query atom of a connected MCIS of isolated atoms has every target atom as a
    // candidate. Each candidate is compared with the mapped atoms only, and the clock is read
    // while the candidates are listed, so the time limit holds.
    const auto graph = isolatedAtoms(32000, 8);
    smsd::MCSOptions opts;
    opts.timeoutMs = 300;
    const auto result = smsd::MCS::find(graph).withOptions(opts).match(graph);
    assert(result.status == smsd::SearchStatus::TimedOut && !result.optimal && !result.complete);
    // Comparing each candidate with every query atom would take about 800 ms here.
    assert(result.elapsedMillis < 550);
}

void testCommonBondPairLimit() {
    namespace detail = smsd::common_bond_detail;
    // 8192 target bonds take 128 words a query bond, so 32768 query bonds fill the limit exactly.
    const auto target = circulant(8192, 1), fits = circulant(32768, 1), over = circulant(32769, 1);
    assert(int64_t(32768) * 128 == detail::MAX_BOND_PAIR_WORDS);
    assert(detail::CommonBondGraphs::supports(fits, target) && !detail::CommonBondGraphs::supports(over, target));
    assert(detail::CommonBondGraphs::supports(target, fits));
    // A dense pair beyond the limit is declined at once: no match, with the limit reached.
    const auto query = circulant(600, 40), dense = circulant(640, 40);
    assert(!detail::CommonBondGraphs::supports(query, dense));
    for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
        auto opts = commonBondOptions(connectivity, smsd::ChemOptions::BondOrderMode::STRICT);
        opts.timeoutMs = 200;
        const auto result = smsd::MCS::find(query).withOptions(opts).match(dense);
        assert(result.status == smsd::SearchStatus::LimitReached && !result.optimal && result.matches.empty());
        assert(result.elapsedMillis < 5000);
    }
}

/** A common-atom or common-bond search of two molecules that the clock stops after `polls` polls; -1 runs it to its end. */
struct PolledSearch {
    bool finished = false;
    int64_t polls = 0;
    smsd::Match match;
};

PolledSearch polledSearch(const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& opts,
                          int64_t polls) {
    const auto& chemistry = opts.chemistry;
    auto atoms = [&](int queryAtom, int targetAtom) {
        return smsd::search_detail::atomCompatible(query, queryAtom, target, targetAtom, chemistry);
    };
    auto bonds = [&](int queryFirst, int querySecond, int targetFirst, int targetSecond) {
        return smsd::search_detail::bondCompatible(query, queryFirst, querySecond, target, targetFirst, targetSecond, chemistry);
    };
    auto key = [](const smsd::MolGraph& graph, int first, int second) {
        return smsd::search_detail::bondCompatibilityKey(graph, first, second);
    };
    PolledSearch result;
    auto stopped = [&] { return ++result.polls > polls && polls >= 0; };
    std::vector<int> mapping;
    std::vector<std::array<int, 4>> found;
    if (opts.objective == smsd::MCSObjective::MCIS) {
        auto outcome = smsd::common_atom_detail::searchCommonAtoms(query, target, atoms, bonds, key, {}, opts.minFragmentAtoms, stopped);
        result.finished = outcome.finished;
        mapping = std::move(outcome.mapping);
        found = std::move(outcome.bonds);
    } else {
        auto outcome = smsd::common_bond_detail::searchCommonBonds(query, target, atoms, bonds, key, {},
                                                                    opts.connectivity == smsd::MCSConnectivity::Connected, stopped);
        result.finished = outcome.finished;
        mapping = std::move(outcome.mapping);
        found = std::move(outcome.bonds);
    }
    if (!mapping.empty()) result.match = smsd::search_detail::makeMatch(mapping, std::move(found));
    return result;
}

void testStarSeedsGrowQuickly() {
    // A seed passes the taken targets next to a hub once, and looks for pairs only at the
    // frontier atoms that might beat the best pair so far, so seeding a star is not cubic in its
    // leaves. Two thousand leaves would take about two million polls if each frontier atom of
    // each step scanned the taken targets again.
    for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
        const bool bonds = objective == smsd::MCSObjective::MCES;
        const auto query = star(2000, bonds), target = star(2000, false);
        auto opts = bonds ? commonBondOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT)
                          : commonAtomOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
        const auto result = polledSearch(query, target, opts, -1);
        assert(result.finished && result.polls < 100000);
        assert(result.match.atomMapping.size() == (bonds ? 2000U : 2001U));
        assert(result.match.bondMapping.size() == (bonds ? 1999U : 2000U));
        assert(smsd::MCS::validateMapping(query, target, result.match, opts).empty());
    }
}

void testIsolatedAtomSeedsGrowIncrementally() {
    // Each seed round of isolated atoms keeps one more atom. A round takes its roots from runs of
    // target cells and passes each target that cannot root a seed only once, so seeding four
    // thousand atoms takes about forty thousand polls, where scanning every pair of atoms in each
    // round would take about fifty million.
    const auto graph = isolatedAtoms(4000, 8);
    const auto opts = commonAtomOptions(smsd::MCSConnectivity::Fragments, smsd::ChemOptions::BondOrderMode::STRICT);
    const auto whole = polledSearch(graph, graph, opts, -1);
    assert(whole.finished && whole.polls < 200000 && whole.match.atomMapping.size() == 4000U);
    assert(smsd::MCS::validateMapping(graph, graph, whole.match, opts).empty());
    // Stopped part way through seeding, the search returns the atoms kept so far.
    const auto stopped = polledSearch(graph, graph, opts, 36000);
    assert(!stopped.finished && stopped.match.atomMapping.size() >= 2000U);
    assert(smsd::MCS::validateMapping(graph, graph, stopped.match, opts).empty());
}

void testSearchesStoppedWhileSeedingKeepTheirSeed() {
    // Stopped a few polls into seeding, which takes most of these searches, each search returns
    // the seed it was growing, which is a valid mapping.
    const auto leaves = star(400, false), oxygen = star(400, true);
    for (auto connectivity : {smsd::MCSConnectivity::Fragments, smsd::MCSConnectivity::Connected})
        for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
            const bool bonds = objective == smsd::MCSObjective::MCES;
            if (!bonds && connectivity == smsd::MCSConnectivity::Connected) continue;
            const auto& query = bonds ? oxygen : leaves;
            auto opts = bonds ? commonBondOptions(connectivity, smsd::ChemOptions::BondOrderMode::STRICT)
                              : commonAtomOptions(connectivity, smsd::ChemOptions::BondOrderMode::STRICT);
            const auto whole = polledSearch(query, leaves, opts, -1);
            assert(whole.finished && whole.polls > 800);
            const auto stopped = polledSearch(query, leaves, opts, 800);
            assert(!stopped.finished && !stopped.match.atomMapping.empty() && !stopped.match.bondMapping.empty());
            assert(smsd::MCS::validateMapping(query, leaves, stopped.match, opts).empty());
        }
}
} // namespace

int main() {
    testOptionRegressions();
    testWeightRegressions();
    testObjectiveOracle();
    testConstrainedBatchObjectives();
    testCoverageMappingValidity();
    testNativeStages();
    testSymmetryCanonicalization();
    testMolecularAutomorphismProperties();
    testSharedDeadline();
    testMappingValidation();
    testTautomerBound();
    testSparseConnectedSelection();
    testFlexibleAromaticExtension();
    testConnectedStatinSeedRecovery();
    testMcGregorStateIsolation();
    testLargeConnectedCoverageSeed();
    testReverseCompleteRingFilter();
    testExhaustiveSmallGraphs();
    testEnumeratedMappings();
    testCommonBondHardFragmentProofs();
    testCommonBondStoppedSearches();
    testCommonBondScoresIgnoreUnreachableNodeLimits();
    testCommonBondOptionsKeepResults();
    testCommonBondSeededOptionsIgnoreUnreachableNodeLimits();
    testFlexibleAromaticityResonanceBonds();
    testCommonAtomHardFragmentProofs();
    testCommonAtomStoppedSearches();
    testCommonAtomScoresAgreeWithGeneralSearch();
    testCommonAtomIneligibleOptionsKeepResults();
    testCommonAtomSeededOptionsAgreeWithGeneralSearch();
    testCommonAtomResonanceBonds();
    testVertexSearchHonoursTimeoutsOnSparseGraphs();
    testCommonBondPairLimit();
    testStarSeedsGrowQuickly();
    testSearchesStoppedWhileSeedingKeepTheirSeed();
    testIsolatedAtomSeedsGrowIncrementally();
    std::cout << "MCS regressions passed\n";
}

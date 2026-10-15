/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/cip.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
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
    std::cout << "MCS regressions passed\n";
}

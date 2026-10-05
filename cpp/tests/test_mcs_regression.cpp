/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/mcs.hpp"

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

// Independent oracle: enumerate every partial injective vertex mapping using
// only these adjacency matrices. No SMSD compatibility, bound, or search helper
// participates in the expected result or in returned-mapping validation.
int maximumAtoms(const Graph& query, const Graph& target,
                 bool induced, bool connected) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    int best = 0;
    std::function<void(int, int)> visit = [&](int atom, int count) {
        if (atom == query.atoms) {
            if (!connected || connectedMapping(query, mapping))
                best = std::max(best, count);
            return;
        }
        for (int other = 0; other < target.atoms; ++other) {
            if (used[other]) continue;
            bool valid = true;
            for (int previous = 0; previous < atom; ++previous) {
                if (mapping[previous] < 0) continue;
                bool qEdge = query.edges[atom][previous] != 0;
                bool tEdge = target.edges[other][mapping[previous]] != 0;
                if ((qEdge && !tEdge) || (induced && qEdge != tEdge)) valid = false;
            }
            if (!valid) continue;
            mapping[atom] = other;
            used[other] = 1;
            visit(atom + 1, count + 1);
            used[other] = 0;
            mapping[atom] = -1;
        }
        visit(atom + 1, count);
    };
    visit(0, 0);
    return best;
}

int maximumScore(const Graph& query, const Graph& target, const smsd::MCSOptions& opts) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    int best = 0;
    std::function<void(int)> visit = [&](int atom) {
        if (atom == query.atoms) {
            if (opts.connectedOnly && !opts.disconnectedMCS && !connectedMapping(query, mapping)) return;
            std::vector<int> seen(query.atoms);
            int fragments = 0;
            for (int root = 0; root < query.atoms; ++root) {
                if (mapping[root] < 0 || seen[root]) continue;
                std::vector<int> component{root};
                seen[root] = 1;
                for (size_t head = 0; head < component.size(); ++head)
                    for (int neighbor = 0; neighbor < query.atoms; ++neighbor)
                        if (mapping[neighbor] >= 0 && !seen[neighbor] && query.edges[component[head]][neighbor]) {
                            seen[neighbor] = 1;
                            component.push_back(neighbor);
                        }
                if (opts.disconnectedMCS && static_cast<int>(component.size()) < opts.minFragmentSize) return;
                ++fragments;
            }
            if (opts.disconnectedMCS && fragments > opts.maxFragments) return;
            int count = 0, bonds = 0;
            double weight = 0.0;
            for (int q = 0; q < query.atoms; ++q) {
                if (mapping[q] < 0) continue;
                ++count;
                if (!opts.atomWeights.empty()) weight += opts.atomWeights[q];
                for (int k = q + 1; k < query.atoms; ++k)
                    if (mapping[k] >= 0 && query.edges[q][k]) ++bonds;
            }
            best = std::max(best, !opts.atomWeights.empty() ? static_cast<int>(weight * 1000.0)
                : opts.maximizeBonds ? bonds : count);
            return;
        }
        visit(atom + 1);
        for (int other = 0; other < target.atoms; ++other) {
            if (used[other]) continue;
            bool valid = true;
            for (int previous = 0; previous < atom; ++previous) {
                if (mapping[previous] < 0) continue;
                bool qEdge = query.edges[atom][previous], tEdge = target.edges[other][mapping[previous]];
                if ((qEdge && !tEdge) || (opts.induced && qEdge != tEdge)) { valid = false; break; }
            }
            if (!valid) continue;
            mapping[atom] = other; used[other] = 1;
            visit(atom + 1);
            mapping[atom] = -1; used[other] = 0;
        }
    };
    visit(0);
    return best;
}

void validate(const Graph& query, const Graph& target,
              const std::map<int, int>& result, bool induced, bool connected) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    for (const auto& [atom, other] : result) {
        assert(atom >= 0 && atom < query.atoms);
        assert(other >= 0 && other < target.atoms);
        assert(!used[other]);
        used[other] = 1;
        mapping[atom] = other;
    }
    for (const auto& [atom, other] : result) {
        for (const auto& [previous, mappedPrevious] : result) {
            bool qEdge = query.edges[atom][previous] != 0;
            bool tEdge = target.edges[other][mappedPrevious] != 0;
            assert(!qEdge || tEdge);
            assert(!induced || qEdge == tEdge);
        }
    }
    assert(!connected || connectedMapping(query, mapping));
}

smsd::MCSOptions options(bool induced, bool connected) {
    smsd::MCSOptions result;
    result.induced = induced;
    result.connectedOnly = connected;
    result.disconnectedMCS = !connected;
    result.timeoutMs = 1000;
    return result;
}

void testOptionRegressions() {
    auto fragments = makeGraph(5, {{0, 1}, {2, 3}, {3, 4}});
    auto connected = smsd::findMCS(
        fragments.molecule, fragments.molecule, {}, options(false, true));
    assert(connected.size() == 3); // Identity still obeys connectedOnly.
    validate(fragments, fragments, connected, false, true);

    auto dOpts = options(false, false);
    dOpts.connectedOnly = true; // Explicit dMCS takes precedence.
    assert(smsd::findMCS(fragments.molecule, fragments.molecule, {}, dOpts).size() == 5);
    dOpts.minFragmentSize = 3;
    dOpts.maxFragments = 1;
    auto fragment = smsd::findMCS(fragments.molecule, fragments.molecule, {}, dOpts);
    assert(fragment.size() == 3);
    validate(fragments, fragments, fragment, false, true);

    auto path = graphFromMask(3, 5), triangle = graphFromMask(3, 7);
    auto induced = smsd::findMCS(path.molecule, triangle.molecule, {}, options(true, true));
    assert(induced.size() == 2); // A cycle is not a linear-path shortcut.
    validate(path, triangle, induced, true, true);

    auto isolates = graphFromMask(3, 0), edge = graphFromMask(3, 1);
    auto forward = smsd::findMCS(isolates.molecule, edge.molecule, {}, options(false, false));
    auto reverse = smsd::findMCS(edge.molecule, isolates.molecule, {}, options(false, false));
    assert(forward.size() == 3 && reverse.size() == 2);
    validate(isolates, edge, forward, false, false);
    validate(edge, isolates, reverse, false, false);

    auto partialInduced = smsd::findMCS(
        edge.molecule, isolates.molecule, {}, options(true, false));
    assert(partialInduced.size() == 2); // Omitted neighbors do not impose degrees.
    validate(edge, isolates, partialInduced, true, false);

    auto star = graphFromMask(4, 7), edgeAndIsolates = graphFromMask(4, 1);
    auto commonEdge = smsd::findMCS(
        edgeAndIsolates.molecule, star.molecule, {}, options(false, true));
    assert(commonEdge.size() == 2); // The largest raw seed can have a worse component.
    validate(edgeAndIsolates, star, commonEdge, false, true);

    auto path4 = graphFromMask(4, 13);
    auto inducedPartial = smsd::findMCS(
        edgeAndIsolates.molecule, path4.molecule, {}, options(true, false));
    assert(inducedPartial.size() == 3); // Invalid augmented seeds must not stop exact search.
    validate(edgeAndIsolates, path4, inducedPartial, true, false);

    auto twoIsolates = graphFromMask(2, 0);
    auto unequal = smsd::findMCS(
        edge.molecule, twoIsolates.molecule, {}, options(false, false));
    assert(unequal.size() == 2); // Reverse containment is not a caller-direction proof.
    validate(edge, twoIsolates, unequal, false, false);

    // Padding bypasses the tiny exact solver and checks heuristic promotions.
    auto paddedEdge = makeGraph(9, {{0, 1}});
    auto paddedStar = makeGraph(9, {{0, 1}, {0, 2}, {0, 3}});
    auto padded = smsd::findMCS(
        paddedEdge.molecule, paddedStar.molecule, {}, options(false, true));
    assert(padded.size() == 2);
    validate(paddedEdge, paddedStar, padded, false, true);
}

void testWeightRegressions() {
    auto clique = graphFromMask(4, 63);
    auto path = makeGraph(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    auto weighted = options(false, true);
    weighted.atomWeights = {1, 1, 1, 10};
    auto result = smsd::findMCS(clique.molecule, path.molecule, {}, weighted);
    assert(result.size() == 2 && result.count(3) == 1);
    validate(clique, path, result, false, true);

    auto fragments = makeGraph(6, {{0, 1}, {1, 2}, {2, 3}, {4, 5}});
    weighted.atomWeights = {1, 1, 1, 1, 5, 5};
    auto heavy = smsd::findMCS(fragments.molecule, fragments.molecule, {}, weighted);
    assert(heavy.size() == 2 && heavy.count(4) && heavy.count(5));
    validate(fragments, fragments, heavy, false, true);

    auto bonded = makeGraph(9, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4},
        {5, 6}, {5, 7}, {5, 8}, {6, 7}, {6, 8}, {7, 8}});
    auto bondOptions = options(false, true);
    bondOptions.maximizeBonds = true;
    auto dense = smsd::findMCS(bonded.molecule, bonded.molecule, {}, bondOptions);
    assert(dense.size() == 4 && dense.count(5) && dense.count(8));
    validate(bonded, bonded, dense, false, true);

    weighted.atomWeights = {1};
    bool threw = false;
    try {
        (void)smsd::findMCS(clique.molecule, clique.molecule, {}, weighted);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw); // Identity must validate per-query options too.

    auto reject = [&](std::vector<double> weights) {
        auto invalid = options(false, true);
        invalid.atomWeights = std::move(weights);
        bool rejected = false;
        try {
            (void)smsd::findMCS(clique.molecule, clique.molecule, {}, invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    };
    reject({std::numeric_limits<double>::quiet_NaN(), 0, 0, 0});
    reject({std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({-std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({std::numeric_limits<double>::max(), 0, 0, 0});
    reject({1.2e6, 1.2e6, -1.2e6, -1.2e6}); // Full score cancels, partial score overflows.
    reject({0, 0, -1.2e6, -1.2e6});

    auto atom = graphFromMask(1, 0);
    auto boundary = options(false, true);
    boundary.atomWeights = {(static_cast<double>(INT_MAX) + 0.5) / 1000.0};
    auto full = smsd::findMCS(atom.molecule, atom.molecule, {}, boundary);
    assert(smsd::detail::mcsScore(atom.molecule, full, boundary) == INT_MAX);
    boundary.atomWeights = {(static_cast<double>(INT_MIN) - 0.5) / 1000.0};
    assert(smsd::detail::mcsScore(atom.molecule, full, boundary) == INT_MIN);
    int q2t[] = {0};
    assert(smsd::detail::mcsScoreFlat(atom.molecule, q2t, 1, 1, boundary) == INT_MIN);

    // Score helpers also protect direct/internal callers that bypass API validation.
    boundary.atomWeights = {std::numeric_limits<double>::max()};
    for (bool flat : {false, true}) {
        bool rejected = false;
        try {
            if (flat) (void)smsd::detail::mcsScoreFlat(atom.molecule, q2t, 1, 1, boundary);
            else (void)smsd::detail::mcsScore(atom.molecule, full, boundary);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    }

    auto signedPath = makeGraph(3, {{0, 1}, {1, 2}});
    auto signedOptions = options(false, true);
    signedOptions.atomWeights = {1, -5, 1};
    auto signedResult = smsd::findMCS(signedPath.molecule, signedPath.molecule, {}, signedOptions);
    assert(smsd::detail::mcsScore(signedPath.molecule, signedResult, signedOptions) == 1000);
    auto signedMappings = smsd::findAllMCS(signedPath.molecule, signedPath.molecule, {}, signedOptions, 10);
    assert(signedMappings.size() == 6);
    for (const auto& mapping : signedMappings) {
        assert(mapping.size() == 1);
        assert(smsd::detail::mcsScore(signedPath.molecule, mapping, signedOptions) == 1000);
    }
    signedOptions.atomWeights = {-1, -2, -3};
    assert(smsd::findMCS(signedPath.molecule, signedPath.molecule, {}, signedOptions).empty());

    auto limited = options(false, false);
    limited.maxFragments = 1;
    limited.atomWeights = {1, 1, 1, 1, 5, 5};
    auto limitedResult = smsd::findMCS(fragments.molecule, fragments.molecule, {}, limited);
    assert(smsd::detail::mcsScore(fragments.molecule, limitedResult, limited) == 10000);
    bondOptions.disconnectedMCS = true;
    bondOptions.maxFragments = 1;
    auto bondFragment = smsd::findMCS(bonded.molecule, bonded.molecule, {}, bondOptions);
    assert(smsd::detail::mcsScore(bonded.molecule, bondFragment, bondOptions) == 6);

    auto large = makeGraph(41, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}});
    auto triangle = graphFromMask(3, 7);
    auto heavySingleton = options(false, true);
    heavySingleton.atomWeights.assign(41, 1);
    heavySingleton.atomWeights[40] = 10;
    auto singleton = smsd::findMCS(large.molecule, triangle.molecule, {}, heavySingleton);
    assert(singleton.size() == 1 && singleton.count(40));
    validate(large, triangle, singleton, false, true);
}

void testObjectiveOracle() {
    size_t cases = 0;
    for (int qMask = 0; qMask < 64; ++qMask) {
        auto query = graphFromMask(4, qMask);
        for (int tMask = 0; tMask < 64; ++tMask) {
            auto target = graphFromMask(4, tMask);
            for (bool induced : {false, true}) for (bool connected : {false, true}) {
                for (int objective = 0; objective < 3; ++objective) {
                    auto opts = options(induced, connected);
                    if (objective == 0) opts.atomWeights = {1, 2, 5, 10};
                    if (objective == 1) opts.atomWeights = {2, -5, 3, -1};
                    if (objective == 2) opts.maximizeBonds = true;
                    auto result = smsd::findMCS(query.molecule, target.molecule, {}, opts);
                    validate(query, target, result, induced, connected);
                    assert(smsd::detail::mcsScore(query.molecule, result, opts) == maximumScore(query, target, opts));
                    ++cases;
                }
            }
        }
    }
    std::cout << "Validated " << cases << " objective oracle cases\n";

    std::mt19937 random(9123);
    size_t fragmentCases = 0;
    for (int sample = 0; sample < 128; ++sample) {
        auto query = graphFromMask(4, random() % 64), target = graphFromMask(4, random() % 64);
        for (bool induced : {false, true}) for (int objective = 0; objective < 3; ++objective)
            for (int minimum : {1, 2}) for (int maximum : {1, 2}) {
                auto opts = options(induced, false);
                opts.minFragmentSize = minimum; opts.maxFragments = maximum;
                if (objective == 0) opts.atomWeights = {1, 2, 5, 10};
                if (objective == 1) opts.atomWeights = {2, -5, 3, -1};
                if (objective == 2) opts.maximizeBonds = true;
                auto result = smsd::findMCS(query.molecule, target.molecule, {}, opts);
                validate(query, target, result, induced, false);
                assert(smsd::detail::mcsScore(query.molecule, result, opts) == maximumScore(query, target, opts));
                ++fragmentCases;
            }
    }
    std::cout << "Validated " << fragmentCases << " constrained-fragment oracle cases\n";
}

void testConstrainedBatchObjectives() {
    auto query = makeGraph(4, {{0, 1}, {1, 2}}, {6, 6, 6, 7});
    auto carbon = makeGraph(3, {{0, 1}, {1, 2}});
    auto nitrogen = makeGraph(1, {}, {7});
    auto weighted = options(false, true);
    weighted.atomWeights = {1, 1, 1, 10};
    std::vector<int> selected;
    auto mappings = smsd::batchMCSConstrained({query.molecule, query.molecule},
        {carbon.molecule, nitrogen.molecule}, {}, weighted, &selected);
    assert(selected == std::vector<int>({1, 0}));
    assert((mappings[0] == std::map<int, int>({{3, 0}})));
    assert(mappings[1].size() == 3);

    auto bondQuery = makeGraph(9, {{0, 1}, {1, 2}, {2, 3}, {3, 4},
        {5, 6}, {5, 7}, {5, 8}, {6, 7}, {6, 8}, {7, 8}});
    auto path = makeGraph(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    auto clique = graphFromMask(4, 63);
    auto bonds = options(false, true);
    bonds.maximizeBonds = true;
    mappings = smsd::batchMCSConstrained({bondQuery.molecule}, {path.molecule, clique.molecule}, {}, bonds, &selected);
    assert(selected == std::vector<int>({1}));
    assert(smsd::detail::mcsScore(bondQuery.molecule, mappings[0], bonds) == 6);
    mappings = smsd::batchMCSConstrained({nitrogen.molecule}, {carbon.molecule}, {}, {}, &selected);
    assert(mappings[0].empty() && selected == std::vector<int>({-1}));
}

void testCoverageMappingValidity() {
    auto q = smsd::parseSMILES("O=C([O-])Cc1ccccc1");
    auto t = smsd::parseSMILES("O=C([O-])CCCc1ccccc1");
    for (bool any : {false, true}) {
        auto result = smsd::findMCSCoverage(q, t, false, any, 1000);
        smsd::ChemOptions chemistry;
        chemistry.matchBondOrder = any ? smsd::ChemOptions::BondOrderMode::ANY : smsd::ChemOptions::BondOrderMode::STRICT;
        assert(smsd::validateMapping(q, t, result, chemistry).empty());
        assert(!result.empty());
    }
    auto flexQuery = smsd::parseSMILES("O=C(Cc1ccccc1)c1ccncc1");
    auto flexTarget = smsd::parseSMILES("O/C(/c1ccncc1)=C\\c1ccccc1");
    const std::map<int,int> witness{{1,11}, {2,10}, {3,9}, {4,8}, {5,1},
                                  {6,2}, {7,3}, {9,12}, {10,13}, {11,14}};
    smsd::ChemOptions flexible;
    flexible.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    flexible.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    for (int aromaticOrder : {1, 4}) {
        for (auto* graph : {&flexQuery, &flexTarget})
            for (int atom = 0; atom < graph->n; ++atom)
                for (int other : graph->neighbors[atom])
                    if (other > atom && graph->bondAromatic(atom, other))
                        graph->setBondOrder(atom, other, aromaticOrder);
        assert(smsd::isValidMCSMapping(flexQuery, flexTarget, witness, flexible));
        auto result = smsd::findMCSCoverage(flexQuery, flexTarget, false, false, 1000);
        assert(result.size() >= witness.size());
        assert(smsd::isValidMCSMapping(flexQuery, flexTarget, result, flexible));
        assert(smsd::detail::largestConnected(flexQuery, result, &flexTarget).size() == result.size());
    }
    auto single = smsd::parseSMILES("CO");
    auto doubled = smsd::parseSMILES("C=O");
    assert(smsd::findMCSCoverage(single, doubled, false, false, 1000).size() == 1);
    assert(smsd::findMCSCoverage(single, doubled, false, true, 1000).size() == 2);
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 3; ++atoms)
        for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask) graphs.push_back(graphFromMask(atoms, mask));
    for (const auto& query : graphs) for (const auto& target : graphs)
        for (bool any : {false, true}) for (bool ring : {false, true}) {
            auto result = smsd::findMCSCoverage(query.molecule, target.molecule, ring, any, 100);
            validate(query, target, result, false, true);
        }
}

void testNativeStages() {
    size_t cases = 0;
    for (int qMask = 0; qMask < 64; ++qMask) {
        auto query = graphFromMask(4, qMask);
        for (int tMask = 0; tMask < 64; ++tMask) {
            auto target = graphFromMask(4, tMask);
            for (bool induced : {false, true}) {
                smsd::detail::GraphBuilder builder(query.molecule, target.molecule, {}, induced);
                smsd::detail::TimeBudget splitBudget(1000), cliqueBudget(1000);
                int64_t nodes = 0;
                auto split = builder.mcSplitSeed(splitBudget, nodes);
                auto clique = builder.maximumCliqueSeed(cliqueBudget);
                int expected = maximumAtoms(query, target, induced, false);
                validate(query, target, split, induced, false);
                validate(query, target, clique, induced, false);
                assert(static_cast<int>(split.size()) == expected);
                assert(static_cast<int>(clique.size()) == expected);
                ++cases;
            }
        }
    }
    std::cout << "Validated " << cases << " native seed oracle cases\n";
}

void testSymmetryCanonicalization() {
    auto cycle = makeGraph(4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
    std::map<int, int> mapping{{0, 1}, {1, 0}, {2, 3}}, expected{{0, 0}, {1, 1}, {2, 2}};
    assert(smsd::canonicalizeMapping(cycle.molecule, cycle.molecule, mapping) == expected);
    for (int qSign : {-1, 1}) for (int qShift = 0; qShift < 4; ++qShift)
        for (int tSign : {-1, 1}) for (int tShift = 0; tShift < 4; ++tShift) {
            std::map<int, int> transformed;
            for (const auto& [q, t] : mapping)
                transformed[(qSign * q + qShift + 4) % 4] = (tSign * t + tShift + 4) % 4;
            assert(smsd::canonicalizeMapping(cycle.molecule, cycle.molecule, transformed) == expected);
            assert(smsd::areMappingsEquivalent(cycle.molecule, cycle.molecule, transformed, expected));
        }
}

void testMolecularAutomorphismProperties() {
    for (const char* smiles : {"[C+].C", "[12C].[13C]", "C=CC", "C[C@H](F)[C@H](F)C"}) {
        auto graph = smsd::parseSMILES(smiles);
        graph.ensureCanonical();
        assert(graph.automorphismGeneratorsTruncated());
        std::map<int, int> mapping{{0, 0}};
        bool rejected = false;
        try { (void)smsd::canonicalizeMapping(graph, graph, mapping); }
        catch (const std::length_error&) { rejected = true; }
        assert(rejected);
        assert(smsd::detail::enumerationMappingKey(graph, graph, mapping) == mapping);
    }
    auto capped = makeGraph(8, {});
    assert(capped.molecule.automorphismGeneratorsTruncated());
}

void testSharedDeadline() {
    std::mt19937 random(781);
    for (int atoms : {24, 32}) {
        std::vector<std::pair<int, int>> queryEdges, targetEdges;
        for (int a = 0; a < atoms; ++a) for (int b = a + 1; b < atoms; ++b) {
            if (random() % 100 < 22) queryEdges.emplace_back(a, b);
            if (random() % 100 < 22) targetEdges.emplace_back(a, b);
        }
        auto query = makeGraph(atoms, queryEdges), target = makeGraph(atoms, targetEdges);
        query.molecule.ensureCanonical(); target.molecule.ensureCanonical();
        auto opts = options(false, true);
        opts.timeoutMs = 5;
        auto started = std::chrono::steady_clock::now();
        auto result = smsd::findMCS(query.molecule, target.molecule, {}, opts);
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - started).count();
        assert(elapsed < 100);
        validate(query, target, result, false, true);
        assert(!smsd::detail::global_deadline::active);
    }
    smsd::detail::global_deadline::set(1000);
    auto previous = smsd::detail::global_deadline::deadline;
    auto atom = graphFromMask(1, 0);
    auto invalid = options(false, true);
    invalid.atomWeights = {std::numeric_limits<double>::quiet_NaN()};
    bool rejected = false;
    try { (void)smsd::findMCS(atom.molecule, atom.molecule, {}, invalid); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected && smsd::detail::global_deadline::active);
    assert(smsd::detail::global_deadline::deadline == previous);
    smsd::detail::global_deadline::clear();
}

void testMappingValidation() {
    auto edge = graphFromMask(2, 1), isolates = graphFromMask(2, 0);
    for (int invalid : {-1, 2, INT_MAX}) {
        std::map<int, int> mapping{{0, 0}, {1, invalid}};
        assert(!smsd::isValidMCSMapping(edge.molecule, edge.molecule, mapping, {}));
        assert(!smsd::validateMapping(edge.molecule, edge.molecule, mapping, {}).empty());
    }
    std::map<int, int> full{{0, 0}, {1, 1}};
    smsd::ChemOptions chemistry;
    assert(smsd::isValidMCSMapping(isolates.molecule, edge.molecule, full, chemistry));
    chemistry.induced = true;
    assert(!smsd::isValidMCSMapping(isolates.molecule, edge.molecule, full, chemistry));
    assert(!smsd::validateMapping(isolates.molecule, edge.molecule, full, chemistry).empty());

    int q2t[] = {0, -1}, t2q[] = {0, -1};
    std::vector<uint64_t> mapped{1};
    assert(!smsd::detail::mappedBondCompatQuery(edge.molecule, isolates.molecule, {}, false, 1, 1, q2t, mapped));
    assert(!smsd::detail::mappedBondCompatTarget(edge.molecule, isolates.molecule, {}, false, 1, 1, t2q, mapped));
    assert(smsd::detail::mappedBondCompatQuery(isolates.molecule, edge.molecule, {}, false, 1, 1, q2t, mapped));
    assert(!smsd::detail::mappedBondCompatQuery(isolates.molecule, edge.molecule, {}, true, 1, 1, q2t, mapped));
    assert(!smsd::isMappingMaximal(isolates.molecule, edge.molecule, {{0, 0}}, {}));
    assert(smsd::isMappingMaximal(isolates.molecule, edge.molecule, {{0, 0}}, chemistry));
    assert(!smsd::isMappingMaximal(edge.molecule, edge.molecule, {{0, INT_MAX}}, {}));

    auto legacy = makeGraph(3, {{0, 1}}).molecule;
    legacy.formalCharge.clear();
    legacy.massNumber.clear();
    legacy.tetraChirality.clear();
    smsd::ChemOptions optionalProperties;
    optionalProperties.matchFormalCharge = true;
    optionalProperties.matchIsotope = true;
    optionalProperties.useChirality = true;
    assert(smsd::isSubstructure(legacy, legacy, optionalProperties, 1000));
    assert(smsd::findSubstructure(legacy, legacy, optionalProperties, 1000).size() == 3);
    for (int atom = 0; atom < legacy.n; ++atom)
        assert(smsd::cip::assignRS(legacy, atom) == smsd::cip::RSLabel::NONE);
    auto components = smsd::detail::splitComponents(legacy);
    assert(components.size() == 2);
    for (const auto& component : components)
        for (int atom = 0; atom < component.n; ++atom) {
            assert(component.formalCharge[atom] == 0);
            assert(component.massNumber[atom] == 0);
            assert(component.tetraChirality[atom] == 0);
        }
}

void testTautomerBound() {
    auto query = graphFromMask(4, 7), target = graphFromMask(4, 7);
    query.molecule.atomicNum.assign(4, 7);
    target.molecule.atomicNum.assign(4, 8);
    query.molecule.tautomerClass.assign(4, 0);
    target.molecule.tautomerClass.assign(4, 0);
    query.molecule.refreshAtomLabels();
    target.molecule.refreshAtomLabels();
    smsd::ChemOptions chemistry;
    chemistry.tautomerAware = true;
    // Tautomer relaxation preserves element identity; the frequency bound may stay loose.
    for (int atom = 0; atom < 4; ++atom)
        assert(!smsd::detail::atomsCompatFast(query.molecule, atom, target.molecule, atom, chemistry));
    assert(smsd::detail::labelFrequencyUpperBound(query.molecule, target.molecule, chemistry) == 4);
    assert(smsd::detail::labelFrequencyUpperBoundDirected(query.molecule, target.molecule, chemistry) == 4);
    assert(smsd::findMCS(query.molecule, target.molecule, chemistry, options(false, true)).empty());
}

void testSparseConnectedSelection() {
    auto query = makeGraph(2048, {{1500, 1501}, {1700, 1701}, {1701, 1702}});
    std::map<int, int> mapping{{1500, 1500}, {1501, 1501},
                              {1700, 1700}, {1701, 1701}, {1702, 1702}};
    auto result = smsd::detail::largestConnected(query.molecule, mapping, &query.molecule);
    assert(result.size() == 3 && result.count(1702));
    validate(query, query, result, false, true);

    auto target = makeGraph(2048, {{1500, 1501}, {1700, 1701}});
    auto common = smsd::detail::largestConnected(query.molecule, mapping, &target.molecule);
    assert(common.size() == 2 && common.count(1500));
    validate(query, target, common, false, true);
}

void testFlexibleAromaticExtension() {
    auto queryMolecule = smsd::parseSMILES("c1ccc(-c2ccccc2)cc1");
    auto targetMolecule = smsd::parseSMILES("c1ccc(Cc2ccccc2)cc1");
    auto toGraph = [](smsd::MolGraph molecule) {
        int atoms = molecule.n;
        std::vector<std::vector<int>> edges(atoms, std::vector<int>(atoms));
        for (int a = 0; a < atoms; ++a)
            for (int b = 0; b < atoms; ++b)
                edges[a][b] = molecule.bondOrder(a, b) != 0;
        return Graph{atoms, std::move(edges), std::move(molecule)};
    };
    auto query = toGraph(std::move(queryMolecule));
    auto target = toGraph(std::move(targetMolecule));
    auto checkChemistry = [&](const std::map<int, int>& mapping) {
        validate(query, target, mapping, false, true);
        for (const auto& [a, b] : mapping) {
            assert(query.molecule.atomicNum[a] == target.molecule.atomicNum[b]);
            assert(query.molecule.formalCharge[a] == target.molecule.formalCharge[b]);
            for (const auto& [k, l] : mapping) {
                int qOrder = query.molecule.bondOrder(a, k);
                if (!qOrder) continue;
                int tOrder = target.molecule.bondOrder(b, l);
                bool flexible = (query.molecule.bondAromatic(a, k) && (tOrder == 1 || tOrder == 2))
                    || (target.molecule.bondAromatic(b, l) && (qOrder == 1 || qOrder == 2));
                assert(qOrder == tOrder || flexible);
            }
        }
    };
    // A complete phenyl ring plus its attached query ring atom maps onto a
    // phenyl ring plus target methylene. FLEXIBLE atom aromaticity permits it.
    std::map<int, int> seven{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {10, 11}, {11, 12}};
    checkChemistry(seven);
    auto result = smsd::findMCS(query.molecule, target.molecule, {}, options(false, true));
    assert(result.size() >= seven.size());
    checkChemistry(result);
}

void testConnectedStatinSeedRecovery() {
    auto query = smsd::parseSMILES("CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4");
    auto target = smsd::parseSMILES("CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C");
    smsd::ChemOptions chemistry;
    chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    const std::map<int, int> witness{{4, 3}, {5, 8}, {6, 7}, {7, 6}, {8, 5},
        {19, 20}, {20, 21}, {21, 22}, {22, 23}, {23, 24}, {24, 25}, {25, 26},
        {26, 9}, {27, 10}, {28, 11}, {29, 12}, {30, 13}, {32, 1}};
    assert(smsd::isValidMCSMapping(query, target, witness, chemistry));
    auto opts = options(false, true);
    opts.timeoutMs = 1000;
    for (int stage : {1, 5}) {
        opts.maxStage = stage;
        auto result = smsd::findMCS(query, target, chemistry, opts);
        assert(result.size() >= witness.size());
        assert(smsd::isValidMCSMapping(query, target, result, chemistry));
        assert(smsd::detail::largestConnected(query, result, &target).size() == result.size());
    }
}

void testMcGregorStateIsolation() {
    size_t cases = 0;
    for (int queryMask = 0; queryMask < 64; ++queryMask) {
        auto query = graphFromMask(4, queryMask);
        for (int targetMask : {0, 7, 18, 23, 32, 41, 56, 63}) {
            auto core = graphFromMask(4, targetMask);
            std::vector<std::pair<int, int>> bonds;
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    if (core.edges[a][b]) bonds.emplace_back(a, b);
            auto target = makeGraph(20, bonds);
            smsd::ChemOptions chemistry;
            chemistry.induced = true;
            for (bool connected : {false, true}) {
                smsd::detail::TimeBudget budget(100);
                auto result = smsd::detail::mcGregorExtend(query.molecule, target.molecule,
                    {{0, 0}}, chemistry, budget, 100, false, false, connected);
                validate(query, target, result, true, connected);
                assert(!result.empty());
                if (connected)
                    assert(static_cast<int>(result.size()) <= maximumAtoms(query, core, true, true));
                ++cases;
            }
        }
    }
    assert(cases == 1024);
}

void testLargeConnectedCoverageSeed() {
    auto query = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C");
    auto target = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O");
    assert(std::min(query.n, target.n) > 50);
    smsd::ChemOptions chemistry;
    chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    auto opts = options(false, true);
    opts.timeoutMs = 2000;
    opts.maxStage = 1;
    auto result = smsd::findMCS(query, target, chemistry, opts);
    assert(result.size() >= 50);
    assert(smsd::isValidMCSMapping(query, target, result, chemistry));
    assert(smsd::detail::largestConnected(query, result, &target).size() == result.size());
}

void testReverseCompleteRingFilter() {
    std::vector<int> queryElements(12, 6), targetElements(8, 6);
    queryElements[10] = queryElements[11] = 7;
    targetElements[6] = targetElements[7] = 7;
    auto query = makeGraph(12, {
        {0, 1}, {1, 2}, {2, 3}, {3, 8}, {8, 9}, {9, 0},
        {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {0, 10}, {10, 11}}, queryElements);
    auto target = makeGraph(8, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 6}, {6, 7}}, targetElements);
    smsd::ChemOptions chemistry;
    chemistry.completeRingsOnly = true;
    auto result = smsd::findMCS(query.molecule, target.molecule, chemistry, options(false, true));
    // The shared atoms of the fused query rings require all ten ring carbons.
    // The target has only six carbons, so the two pendant nitrogens are optimal.
    assert(result.size() == 2 && result.count(10) && result.count(11));
    assert(result.at(10) >= 6 && result.at(11) >= 6);
    validate(query, target, result, false, true);
}

void testExhaustiveSmallGraphs() {
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 4; ++atoms)
        for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask)
            graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) {
        for (const auto& target : graphs) {
            for (bool induced : {false, true}) {
                for (bool connected : {false, true}) {
                    auto result = smsd::findMCS(
                        query.molecule, target.molecule, {}, options(induced, connected));
                    validate(query, target, result, induced, connected);
                    assert(static_cast<int>(result.size()) == maximumAtoms(query, target, induced, connected));
                    ++cases;
                }
            }
        }
    }
    assert(cases == 23104);
    std::cout << "Validated " << cases << " exhaustive small-graph MCS cases\n";
}

void testEnumeratedMappings() {
    auto isolates = graphFromMask(2, 0);
    auto path = graphFromMask(3, 5);
    auto result = smsd::findAllMCS(
        isolates.molecule, path.molecule, {}, options(true, false), 10);
    assert(!result.empty());
    for (const auto& mapping : result) {
        assert(mapping.size() == 2);
        validate(isolates, path, mapping, true, false);
    }

    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 3; ++atoms)
        for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask)
            graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) {
        for (const auto& target : graphs) {
            for (bool induced : {false, true}) {
                for (bool connected : {false, true}) {
                    auto mappings = smsd::findAllMCS(
                        query.molecule, target.molecule, {}, options(induced, connected), 3);
                    int expected = maximumAtoms(query, target, induced, connected);
                    assert(expected == 0 || !mappings.empty());
                    for (const auto& mapping : mappings) {
                        assert(static_cast<int>(mapping.size()) == expected);
                        validate(query, target, mapping, induced, connected);
                    }
                    ++cases;
                }
            }
        }
    }
    assert(cases == 576);
    std::cout << "Validated " << cases << " exhaustive enumeration cases\n";
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
    return 0;
}

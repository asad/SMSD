/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Independent graph oracles for complete common subgraphs, result limits
 * and molecular substructure/MCS mappings.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <set>
#include <map>
#include <functional>

#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

using Mapping = std::vector<std::pair<int, int>>;
using Bonds = std::map<std::pair<int, int>, int>;

using smsd_test::require;

struct Graph {
    int n;
    std::vector<std::vector<int>> adj;
};

static Graph graphFromMask(int n, std::uint64_t mask) {
    Graph graph{n, std::vector<std::vector<int>>(n)};
    int bit = 0;
    for (int first = 0; first < n; ++first) for (int second = first + 1; second < n; ++second, ++bit)
        if ((mask & (std::uint64_t(1) << bit)) != 0) {
            graph.adj[first].push_back(second); graph.adj[second].push_back(first);
        }
    return graph;
}

static smsd::MolGraph molecule(int n, const Bonds& bonds, bool uniqueLabels = false) {
    std::vector<std::vector<int>> adjacent(n), orders(n);
    std::vector<int> labels(n, 6);
    if (uniqueLabels) std::iota(labels.begin(), labels.end(), 6);
    for (const auto& bond : bonds) {
        const auto [first, second] = bond.first;
        adjacent[first].push_back(second); adjacent[second].push_back(first);
        orders[first].push_back(bond.second); orders[second].push_back(bond.second);
    }
    return smsd::MolGraph::Builder().atomCount(n).atomicNumbers(labels).setNeighbors(adjacent).setBondOrders(orders).build(false);
}

// Small searches here must prove their answers; the limit only guards against a hang, so it
// allows for slow Debug builds sharing a runner with other tests.
static constexpr int kProofTimeoutMs = 30000;

static int completeSubgraphOracle(const Graph& graph) {
    int maximum = 0;
    for (unsigned mask = 1; mask < (1u << graph.n); ++mask) {
        std::vector<int> vertices;
        bool complete = true;
        for (int atom = 0; atom < graph.n && complete; ++atom) {
            if ((mask & (1u << atom)) == 0) continue;
            for (int previous : vertices)
                if (!std::binary_search(graph.adj[atom].begin(), graph.adj[atom].end(), previous)) { complete = false; break; }
            vertices.push_back(atom);
        }
        if (complete) maximum = std::max(maximum, static_cast<int>(vertices.size()));
    }
    return maximum;
}

static void verifyCompleteSubgraph(const Graph& graph) {
    Bonds source, complete;
    for (int atom = 0; atom < graph.n; ++atom) {
        for (int neighbor : graph.adj[atom]) if (atom < neighbor) source[{atom, neighbor}] = 1;
        for (int neighbor = atom + 1; neighbor < graph.n; ++neighbor) complete[{atom, neighbor}] = 1;
    }
    const auto result = smsd::MCS::find(molecule(graph.n, source)).withTimeout(kProofTimeoutMs).match(molecule(graph.n, complete));
    std::string context;
    if (!result.optimal || !result.complete) {
        std::uint64_t mask = 0;
        int bit = 0;
        for (int first = 0; first < graph.n; ++first)
            for (int second = first + 1; second < graph.n; ++second, ++bit)
                if (source.count({first, second})) mask |= std::uint64_t(1) << bit;
        std::ostringstream details;
        details << ": n=" << graph.n << " edges={";
        for (const auto& bond : source)
            details << '(' << bond.first.first << ',' << bond.first.second << ')';
        details << "} mask=" << mask << " status=" << static_cast<int>(result.status)
                << " optimal=" << result.optimal << " complete=" << result.complete
                << " score=" << result.score << " nodes=" << result.nodesVisited;
        context = details.str();
    }
    require(result.optimal && result.complete, "small complete-subgraph search did not prove its answer" + context);
    require(result.score == completeSubgraphOracle(graph), "MCIS complete-subgraph size differs from exhaustive oracle");
}

static void testCompleteSubgraphs() {
    std::size_t count = 0;
    for (int atoms = 0; atoms <= 6; ++atoms) {
        const auto masks = std::uint64_t(1) << (atoms * (atoms - 1) / 2);
        for (std::uint64_t mask = 0; mask < masks; ++mask) { verifyCompleteSubgraph(graphFromMask(atoms, mask)); ++count; }
    }
    std::mt19937 random(0xC1192026);
    for (int trial = 0; trial < 256; ++trial) {
        const int atoms = 7 + trial % 4;
        std::uint64_t mask = 0;
        for (int bit = 0; bit < atoms * (atoms - 1) / 2; ++bit)
            if (random() % 100 < static_cast<unsigned>(trial % 101)) mask |= std::uint64_t(1) << bit;
        verifyCompleteSubgraph(graphFromMask(atoms, mask));
    }
    std::cout << "Complete-subgraph oracle passed: " << count << " exhaustive graphs and 256 seeded graphs.\n";
}

static void testCompleteInputConnectivity() {
    const auto path = molecule(3, {{{0, 1}, 1}, {{1, 2}, 1}});
    const auto triangle = molecule(3, {{{0, 1}, 1}, {{0, 2}, 1}, {{1, 2}, 1}});
    const auto isolated = molecule(3, {});
    for (bool reverse : {false, true}) {
        const auto& query = reverse ? triangle : path;
        const auto& target = reverse ? path : triangle;
        std::set<Mapping> expected;
        for (int queryFirst = 0; queryFirst < 3; ++queryFirst)
            for (int querySecond = queryFirst + 1; querySecond < 3; ++querySecond) {
                if (!query.hasBond(queryFirst, querySecond)) continue;
                for (int targetFirst = 0; targetFirst < 3; ++targetFirst)
                    for (int targetSecond = 0; targetSecond < 3; ++targetSecond)
                        if (targetFirst != targetSecond && target.hasBond(targetFirst, targetSecond))
                            expected.insert({{queryFirst, targetFirst}, {querySecond, targetSecond}});
            }
        require(expected.size() == 12, "two-atom clique oracle must have twelve mappings");
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
            smsd::MCSOptions options;
            options.connectivity = connectivity;
            options.timeoutMs = kProofTimeoutMs;
            const auto single = smsd::MCS::find(query).withOptions(options).match(target);
            require(single.status == smsd::SearchStatus::Completed && single.optimal && single.complete
                    && single.score == 2.0 && single.matches.size() == 1,
                    "a complete input must give a proved two-atom clique");
            require(expected.count(Mapping(single.matches[0].atomMapping.begin(), single.matches[0].atomMapping.end())) == 1,
                    "single clique witness differs from the injection oracle");
            require(smsd::MCS::validateMapping(query, target, single.matches[0], options).empty(),
                    "complete-input search returned an invalid single witness");
            const auto all = smsd::MCS::find(query).withOptions(options).matchAll(target, 100);
            require(all.status == smsd::SearchStatus::Completed && all.optimal && all.complete && all.score == 2.0,
                    "all clique mappings must be proved");
            std::set<Mapping> found;
            for (const auto& match : all.matches) {
                require(smsd::MCS::validateMapping(query, target, match, options).empty(),
                        "complete-input search returned an invalid witness");
                require(match.fragments.size() == 1, "a nonempty induced match of a complete graph is connected");
                require(found.emplace(match.atomMapping.begin(), match.atomMapping.end()).second,
                        "complete-input search duplicated a mapping");
            }
            require(found == expected, "complete-input search omitted tied optimal mappings");
            options.timeoutMs = 0;
            const auto timedOut = smsd::MCS::find(query).withOptions(options).match(target);
            require(timedOut.status == smsd::SearchStatus::TimedOut && !timedOut.optimal && !timedOut.complete
                    && timedOut.matches.empty(), "expired complete-input searches must do no work");
            options.timeoutMs = kProofTimeoutMs;
            options.cancellation = std::make_shared<smsd::CancellationToken>();
            options.cancellation->cancel();
            const auto cancelled = smsd::MCS::find(query).withOptions(options).match(target);
            require(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete
                    && cancelled.matches.empty(), "complete-input searches must honour cancellation");
            options.cancellation.reset();
            options.nodeLimit = 1;
            const auto capped = smsd::MCS::find(query).withOptions(options).match(target);
            require(capped.status == smsd::SearchStatus::LimitReached && !capped.optimal && !capped.complete,
                    "complete-input searches must honour the node cap");
            options.nodeLimit = 0;
            options.atomWeights = {1.0, 2.0, 4.0};
            const auto weighted = smsd::MCS::find(query).withOptions(options).match(target);
            require(weighted.status == smsd::SearchStatus::Completed && weighted.optimal && weighted.complete
                    && weighted.score == 6.0, "complete-input search changed the weighted objective");
            options.atomWeights.clear();
            options.excludedTargetAtoms = {0, 1};
            options.minFragmentAtoms = 2;
            const auto tooSmall = smsd::MCS::find(query).withOptions(options).match(target);
            require(tooSmall.status == smsd::SearchStatus::Completed && tooSmall.optimal && tooSmall.complete
                    && tooSmall.score == 0.0 && tooSmall.matches.empty(),
                    "complete-input search must retain exclusions and minimum fragment size");
        }
        const auto& noncompleteQuery = reverse ? isolated : path;
        const auto& noncompleteTarget = reverse ? path : isolated;
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
            smsd::MCSOptions options;
            options.connectivity = connectivity;
            options.timeoutMs = kProofTimeoutMs;
            const auto result = smsd::MCS::find(noncompleteQuery).withOptions(options).match(noncompleteTarget);
            require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete
                    && result.score == (connectivity == smsd::MCSConnectivity::Connected ? 1.0 : 2.0),
                    "noncomplete inputs must preserve the connectivity distinction");
            for (const auto& match : result.matches)
                require(smsd::MCS::validateMapping(noncompleteQuery, noncompleteTarget, match, options).empty(),
                        "noncomplete-input search returned an invalid witness");
        }
    }
}

static bool validMapping(const Mapping& mapping, const Bonds& query,
                         const Bonds& target, bool anyOrder) {
    std::map<int, int> mapped;
    std::set<int> used;
    for (const auto& [q, t] : mapping) {
        if (!mapped.emplace(q, t).second || !used.insert(t).second) return false;
    }
    for (const auto& [edge, order] : query) {
        if (!mapped.count(edge.first) || !mapped.count(edge.second)) continue;
        const int a = mapped.at(edge.first), b = mapped.at(edge.second);
        const auto other = target.find({std::min(a, b), std::max(a, b)});
        if (other == target.end() || (!anyOrder && other->second != order)) return false;
    }
    return true;
}

static Bonds bondsFromMask(int n, std::uint64_t mask) {
    Bonds bonds;
    int bit = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j, ++bit)
            if (mask & (std::uint64_t(1) << bit)) bonds[{i, j}] = 1;
    return bonds;
}

static std::set<Mapping> embeddingOracle(int nQuery, int nTarget,
                                        const Bonds& query, const Bonds& target) {
    std::set<Mapping> embeddings;
    Mapping current;
    std::vector<bool> used(nTarget, false);
    std::function<void(int)> enumerate = [&](int q) {
        if (q == nQuery) {
            if (validMapping(current, query, target, false)) embeddings.insert(current);
            return;
        }
        for (int t = 0; t < nTarget; ++t) {
            if (used[t]) continue;
            current.emplace_back(q, t);
            used[t] = true;
            enumerate(q + 1);
            used[t] = false;
            current.pop_back();
        }
    };
    enumerate(0);
    return embeddings;
}

static void testSubstructureOracle() {
    std::size_t count = 0;
    for (int nQuery = 1; nQuery <= 3; ++nQuery) {
        for (int nTarget = 1; nTarget <= 4; ++nTarget) {
            std::vector<std::pair<int, int>> compat;
            for (int q = 0; q < nQuery; ++q)
                for (int t = 0; t < nTarget; ++t) compat.emplace_back(q, t);
            for (int qm = 0; qm < (1 << (nQuery * (nQuery - 1) / 2)); ++qm) {
                const auto query = bondsFromMask(nQuery, qm);
                for (int tm = 0; tm < (1 << (nTarget * (nTarget - 1) / 2)); ++tm) {
                    const auto target = bondsFromMask(nTarget, tm);
                    const auto expected = embeddingOracle(nQuery, nTarget, query, target);
                    for (int cap : {0, 1, 100}) {
                        if (cap == 0) {
                            bool rejected = false;
                            try { smsd::Substructure::find(molecule(nQuery, query)).matchAll(molecule(nTarget, target), 0); }
                            catch (const std::invalid_argument&) { rejected = true; }
                            require(rejected, "zero result cap must be rejected");
                            continue;
                        }
                        const auto result = smsd::Substructure::find(molecule(nQuery, query)).withTimeout(kProofTimeoutMs).matchAll(molecule(nTarget, target), cap);
                        std::vector<Mapping> found;
                        for (const auto& match : result.matches) found.emplace_back(match.atomMapping.begin(), match.atomMapping.end());
                        require(found.size() == std::min(expected.size(), static_cast<std::size_t>(cap)),
                                "standalone substructure differs from injection oracle");
                        std::set<Mapping> seen;
                        for (const auto& mapping : found) {
                            require(expected.count(mapping) == 1, "standalone oracle found an invalid embedding");
                            require(seen.insert(mapping).second, "standalone oracle found duplicate embeddings");
                        }
                    }
                    ++count;
                }
            }
        }
    }
    std::cout << "Substructure injection oracle passed: " << count << " graph pairs.\n";
}

static void testMolecularValidity() {
    const Bonds query{{{0, 1}, 1}, {{1, 2}, 1}}, wrongOrders{{{0, 1}, 2}, {{1, 2}, 2}};
    const auto source = molecule(3, query, true);
    for (const auto& target : {Bonds{}, wrongOrders}) {
        const auto destination = molecule(3, target, true);
        const auto result = smsd::MCS::find(source).withTimeout(kProofTimeoutMs).match(destination);
        require(result.optimal && result.score <= 2.0, "strict bonds falsely reach the full graph");
        for (const auto& match : result.matches)
            require(smsd::MCS::validateMapping(source, destination, match).empty(), "search returned incompatible bonds");
    }
    smsd::MCSOptions options; options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    require(smsd::MCS::find(source).withOptions(options).match(molecule(3, wrongOrders, true)).score == 3.0,
            "any-order matching must accept changed bonds");
    const Bonds disconnected{{{1, 2}, 1}, {{2, 3}, 1}};
    require(smsd::MCS::find(molecule(4, disconnected)).match(molecule(4, disconnected)).score == 3.0,
            "connected search must retain the largest component");
    const auto stopped = smsd::MCS::find(source).withTimeout(0).match(source);
    require(stopped.status == smsd::SearchStatus::TimedOut && !stopped.optimal && stopped.matches.empty(),
            "expired searches must not perform work");
}

static void testRandomMolecularValidity() {
    std::mt19937 random(0x9C52026);
    for (int trial = 0; trial < 128; ++trial) {
        constexpr int nQuery = 4, nTarget = 5;
        Bonds query, target;
        for (int i = 0; i < nQuery; ++i) for (int j = i + 1; j < nQuery; ++j)
            if (random() % 3 != 0) query[{i, j}] = 1 + random() % 2;
        for (int i = 0; i < nTarget; ++i) for (int j = i + 1; j < nTarget; ++j)
            if (random() % 3 != 0) target[{i, j}] = 1 + random() % 2;
        std::vector<std::pair<int, int>> compat;
        for (int q = 0; q < nQuery; ++q) for (int t = 0; t < nTarget; ++t)
            if (random() % 4 != 0) compat.emplace_back(q, t);
        for (bool anyOrder : {false, true}) {
            for (int cap : {0, 1, 8}) {
                if (cap == 0) continue;
                const auto source = molecule(nQuery, query), destination = molecule(nTarget, target);
                smsd::MCSOptions options;
                options.chemistry.matchBondOrder = anyOrder ? smsd::ChemOptions::BondOrderMode::ANY : smsd::ChemOptions::BondOrderMode::STRICT;
                const auto result = smsd::MCS::find(source).withOptions(options).matchAll(destination, cap);
                std::vector<Mapping> candidates;
                for (const auto& match : result.matches) {
                    require(smsd::MCS::validateMapping(source, destination, match, options).empty(), "selected graph is invalid");
                    candidates.emplace_back(match.atomMapping.begin(), match.atomMapping.end());
                }
                require(candidates.size() <= static_cast<std::size_t>(cap),
                        "random search exceeds its candidate cap");
                std::set<Mapping> seen;
                for (const auto& mapping : candidates) {
                    require(validMapping(mapping, query, target, anyOrder),
                            "random search returned incompatible bonds");
                    require(seen.insert(mapping).second, "random search duplicates a candidate");

                }
                if (!candidates.empty())
                    require(candidates[0].size() == static_cast<std::size_t>(result.score),
                            "candidate cap removed the best search mapping");
            }
        }
    }
    std::cout << "Molecular validity passed: 128 seeded graph pairs across orders/caps.\n";
}

static void testMolecular() {
    testMolecularValidity();
    testRandomMolecularValidity();
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"complete", testCompleteSubgraphs}, {"complete-input", testCompleteInputConnectivity},
        {"substructure", testSubstructureOracle}, {"molecular", testMolecular}});
}

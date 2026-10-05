/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Independent small-graph oracles for capped clique enumeration and the
 * compatibility-table substructure/MCS helpers.
 */
#include "smsd/clique_solver.hpp"

#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using Mapping = std::vector<std::pair<int, int>>;
using Bonds = std::map<std::pair<int, int>, int>;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static smsd::clique::ProductGraph graphFromMask(int n, std::uint64_t mask) {
    std::vector<std::pair<int, int>> edges;
    int bit = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j, ++bit)
            if (mask & (std::uint64_t(1) << bit)) edges.emplace_back(i, j);
    smsd::clique::ProductGraph graph;
    graph.build(std::vector<smsd::clique::ProductVertex>(n), edges);
    return graph;
}

static std::set<std::vector<int>> cliqueOracle(const smsd::clique::ProductGraph& graph) {
    std::set<std::vector<int>> result;
    std::size_t maximum = 0;
    for (unsigned mask = 1; mask < (1u << graph.n); ++mask) {
        std::vector<int> vertices;
        bool clique = true;
        for (int i = 0; i < graph.n && clique; ++i) {
            if (!(mask & (1u << i))) continue;
            for (int j : vertices)
                if (!std::binary_search(graph.adj[i].begin(), graph.adj[i].end(), j)) {
                    clique = false;
                    break;
                }
            vertices.push_back(i);
        }
        if (!clique || vertices.size() < maximum) continue;
        if (vertices.size() > maximum) {
            maximum = vertices.size();
            result.clear();
        }
        result.insert(vertices);
    }
    return result;
}

static void verifyClique(const smsd::clique::ProductGraph& graph) {
    const auto expected = cliqueOracle(graph);
    const int maximum = expected.empty() ? 0 : static_cast<int>(expected.begin()->size());
    for (int incumbent : {0, std::max(0, maximum - 1), maximum}) {
        for (int cap : {-1, 0, 1, 2, 8, 100}) {
            const auto found = smsd::clique::findMaxCliques(graph, cap, 1000, incumbent);
            require(!found.timed_out, "small clique oracle unexpectedly timed out");
            require(found.max_size == maximum, "clique maximum differs from exhaustive oracle");
            require(found.cliques.size() == std::min(expected.size(),
                static_cast<std::size_t>(std::max(0, cap))), "clique cap/tie count is incorrect");
            std::set<std::vector<int>> seen;
            for (auto clique : found.cliques) {
                std::sort(clique.begin(), clique.end());
                require(expected.count(clique) == 1, "returned clique is not an oracle maximum");
                require(seen.insert(clique).second, "returned maximum clique is duplicated");
            }
        }
    }
}

static void testCliqueEnumeration() {
    std::size_t count = 0;
    for (int n = 0; n <= 6; ++n) {
        const auto masks = std::uint64_t(1) << (n * (n - 1) / 2);
        for (std::uint64_t mask = 0; mask < masks; ++mask) {
            verifyClique(graphFromMask(n, mask));
            ++count;
        }
    }
    std::mt19937 rng(0xC1192026);
    for (int trial = 0; trial < 256; ++trial) {
        const int n = 7 + trial % 4;
        std::uint64_t mask = 0;
        for (int bit = 0; bit < n * (n - 1) / 2; ++bit)
            if (rng() % 100 < static_cast<unsigned>(trial % 101)) mask |= std::uint64_t(1) << bit;
        verifyClique(graphFromMask(n, mask));
    }
    std::cout << "Clique oracle passed: " << count
              << " exhaustive graphs and 256 seeded graphs across caps/incumbents.\n";
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

static void testStandaloneSubstructure() {
    const auto one = smsd::clique::substructureMatch(1, 1, {{0, 0}}, {}, {});
    require(one.size() == 1 && one[0] == Mapping{{0, 0}},
            "standalone substructure helper loses the empty-start one-atom embedding");

    const std::vector<std::pair<int, int>> compat{{0, 0}, {0, 1}, {0, 2},
                                                {1, 0}, {1, 1}, {1, 2}};
    const Bonds query{{{0, 1}, 1}};
    const Bonds target{{{0, 1}, 1}, {{1, 2}, 2}};
    for (bool anyOrder : {false, true}) {
        std::set<Mapping> expected;
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            Mapping mapping{{0, a}, {1, b}};
            if (validMapping(mapping, query, target, anyOrder)) expected.insert(mapping);
        }
        for (int cap : {0, 1, 8}) {
            const auto found = smsd::clique::substructureMatch(
                2, 3, compat, query, target, anyOrder, 1000, cap);
            require(found.size() == std::min(expected.size(), static_cast<std::size_t>(cap)),
                    "standalone substructure cap/bond matching is incorrect");
            std::set<Mapping> seen;
            for (const auto& mapping : found) {
                require(expected.count(mapping) == 1, "standalone helper returned an invalid embedding");
                require(seen.insert(mapping).second, "standalone helper returned a duplicate embedding");
            }
        }
    }
    require(smsd::clique::substructureMatch(2, 3, {{0, 0}}, query, target).empty(),
            "query with no compatible target for an atom should not match");
    require(smsd::clique::substructureMatch(1, 1, {{0, -1}, {0, 1}}, {}, {}).empty(),
            "standalone substructure helper accepted out-of-range target atoms");
    std::cout << "Standalone substructure regressions passed.\n";
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
                        const auto found = smsd::clique::substructureMatch(
                            nQuery, nTarget, compat, query, target, false, 1000, cap);
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

static void testPipelineValidity() {
    const Bonds query{{{0, 1}, 1}, {{1, 2}, 1}};
    // Query and target have the same connectivity but atom compatibilities force
    // every target bond order to disagree with its query counterpart.
    const Bonds wrongOrders{{{0, 1}, 2}, {{1, 2}, 2}};
    const std::vector<std::pair<int, int>> compat{{0, 0}, {1, 1}, {2, 2}};
    for (const auto& target : {Bonds{}, wrongOrders}) {
        const auto result = smsd::clique::findMCSPipeline(
            compat, query, target, 3, 3, false, 1000, 8);
        for (const auto& mapping : result.candidates)
            require(validMapping(mapping, query, target, false),
                    "pipeline greedy result contains absent or incompatible bonds");
        require(result.best_size <= 2, "pipeline falsely reached full-size upper bound");
    }
    const auto anyOrder = smsd::clique::findMCSPipeline(
        compat, query, wrongOrders, 3, 3, true, 1000, 8);
    require(anyOrder.best_size == 3, "any-order pipeline should accept changed bond orders");
    for (const auto& mapping : anyOrder.candidates)
        require(validMapping(mapping, query, wrongOrders, true), "any-order pipeline mapping is invalid");
    const auto capped = smsd::clique::findMCSPipeline(
        compat, query, wrongOrders, 3, 3, true, 1000, 0);
    require(capped.candidates.empty(), "pipeline zero result cap still returns candidates");

    const Bonds disconnected{{{1, 2}, 1}, {{2, 3}, 1}};
    const std::vector<std::pair<int, int>> disconnectedCompat{{0, 0}, {1, 1}, {2, 2}, {3, 3}};
    const auto largest = smsd::clique::findMCSPipeline(
        disconnectedCompat, disconnected, disconnected, 4, 4, false, 1000, 1);
    require(largest.best_size == 3 && largest.candidates.size() == 1,
            "pipeline must keep the largest connected component rather than the first one");

    std::vector<std::vector<int>> adjacency{{1}, {0, 2}, {1}};
    std::vector<std::vector<int>> byQuery{{0}, {1}, {2}};
    const auto expired = smsd::clique::mcgregorDFSExtend(
        {{0, 0}}, compat, byQuery, query, query, adjacency,
        3, 3, false, std::chrono::steady_clock::now() - std::chrono::seconds(1));
    require(expired == Mapping{{0, 0}}, "expired DFS extension must return its seed without search");
    std::cout << "Pipeline validity regressions passed.\n";
}

static void testRandomPipelineValidity() {
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
                const auto result = smsd::clique::findMCSPipeline(
                    compat, query, target, nQuery, nTarget, anyOrder, 1000, cap);
                require(result.candidates.size() <= static_cast<std::size_t>(cap),
                        "random pipeline exceeds its candidate cap");
                std::set<Mapping> seen;
                for (const auto& mapping : result.candidates) {
                    require(validMapping(mapping, query, target, anyOrder),
                            "random pipeline returned incompatible bonds");
                    require(seen.insert(mapping).second, "random pipeline duplicates a candidate");
                    for (const auto& pair : mapping)
                        require(std::find(compat.begin(), compat.end(), pair) != compat.end(),
                                "random pipeline uses an incompatible atom pair");
                }
                if (!result.candidates.empty())
                    require(result.candidates[0].size() == static_cast<std::size_t>(result.best_size),
                            "candidate cap removed the best pipeline mapping");
            }
        }
    }
    std::cout << "Pipeline validity passed: 128 seeded graph pairs across orders/caps.\n";
}

int main(int argc, char** argv) {
    try {
        const std::string category = argc > 1 ? argv[1] : "all";
        if (category == "all" || category == "clique") testCliqueEnumeration();
        if (category == "all" || category == "substructure") {
            testStandaloneSubstructure();
            testSubstructureOracle();
        }
        if (category == "all" || category == "pipeline") {
            testPipelineValidity();
            testRandomPipelineValidity();
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Clique/helper regression failed: " << error.what() << '\n';
        return 1;
    }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * General matching regression tests with an independent exhaustive oracle.
 */
#include "smsd/general_matching.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

using Graph = std::vector<std::vector<int>>;

static void addEdge(Graph& graph, int a, int b) {
    graph[a].push_back(b);
    graph[b].push_back(a);
}

// Consider every possible partner for the first available vertex, or leave it
// unmatched. Unlike the implementation under test, this enumerates matchings
// directly and does not use augmenting paths or blossom contraction.
static int oracle(const Graph& graph, std::uint32_t available,
                  std::vector<int>& memo) {
    if (available == 0) return 0;
    int& cached = memo[available];
    if (cached != -1) return cached;
    int first = 0;
    while ((available & (std::uint32_t(1) << first)) == 0) ++first;
    const auto remaining = available & ~(std::uint32_t(1) << first);
    cached = oracle(graph, remaining, memo);
    for (int neighbor : graph[first]) {
        if ((remaining & (std::uint32_t(1) << neighbor)) != 0)
            cached = std::max(cached, 1 + oracle(
                graph, remaining & ~(std::uint32_t(1) << neighbor), memo));
    }
    return cached;
}

static int validate(const Graph& graph) {
    const auto mate = smsd::detail::maximumMatching(graph);
    if (mate.size() != graph.size())
        throw std::runtime_error("matching has the wrong number of vertices");
    int edges = 0;
    for (std::size_t v = 0; v < mate.size(); ++v) {
        const int u = mate[v];
        if (u == -1) continue;
        if (u < 0 || static_cast<std::size_t>(u) >= mate.size()
            || static_cast<std::size_t>(u) == v || mate[u] != static_cast<int>(v))
            throw std::runtime_error("matching is not a set of disjoint pairs");
        if (std::find(graph[v].begin(), graph[v].end(), u) == graph[v].end())
            throw std::runtime_error("matching uses a nonexistent edge");
        if (static_cast<int>(v) < u) ++edges;
    }
    return edges;
}

static void compareWithOracle(const Graph& graph) {
    std::vector<int> memo(std::size_t(1) << graph.size(), -1);
    const int expected = oracle(graph,
        (std::uint32_t(1) << graph.size()) - 1, memo);
    if (validate(graph) != expected)
        throw std::runtime_error("matching cardinality differs from exhaustive oracle");
}

static void checkExamples() {
    Graph triangleWithLeaves(6);
    for (const auto edge : {std::pair<int, int>{0, 1}, {1, 2}, {2, 0},
                           {0, 3}, {1, 4}, {2, 5}})
        addEdge(triangleWithLeaves, edge.first, edge.second);
    if (validate(triangleWithLeaves) != 3)
        throw std::runtime_error("blossom with three leaves needs a perfect matching");

    // Azulene: fused five- and seven-membered carbon rings sharing edge 3--4.
    Graph azulene(10);
    for (const auto edge : {std::pair<int, int>{0, 1}, {1, 2}, {2, 3}, {3, 4},
                           {4, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 3}})
        addEdge(azulene, edge.first, edge.second);
    if (validate(azulene) != 5)
        throw std::runtime_error("azulene needs a perfect matching");
    compareWithOracle(azulene);

    Graph petersen(10);
    for (int i = 0; i < 5; ++i) {
        addEdge(petersen, i, (i + 1) % 5);
        addEdge(petersen, i, i + 5);
        addEdge(petersen, i + 5, (i + 2) % 5 + 5);
    }
    if (validate(petersen) != 5)
        throw std::runtime_error("Petersen graph needs a perfect matching");
    compareWithOracle(petersen);

    for (int n : {127, 128}) {
        Graph complete(n);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) addEdge(complete, i, j);
        if (validate(complete) != n / 2)
            throw std::runtime_error("complete graph matching cardinality is incorrect");
    }
}

static std::size_t checkAllSmallGraphs() {
    std::size_t count = 0;
    for (int n = 0; n <= 6; ++n) {
        const int possibleEdges = n * (n - 1) / 2;
        const auto graphCount = std::uint32_t(1) << possibleEdges;
        for (std::uint32_t mask = 0; mask < graphCount; ++mask) {
            Graph graph(n);
            int edge = 0;
            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j, ++edge)
                    if ((mask & (std::uint32_t(1) << edge)) != 0) addEdge(graph, i, j);
            compareWithOracle(graph);
            ++count;
        }
    }
    return count;
}

static std::size_t checkRandomGraphs() {
    std::mt19937 rng(0x5A5D2026);
    std::size_t count = 0;
    for (int n = 7; n <= 12; ++n) {
        for (int threshold : {10, 30, 50, 80}) {
            for (int trial = 0; trial < 128; ++trial) {
                Graph graph(n);
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        if (rng() % 100 >= static_cast<unsigned>(threshold)) continue;
                        addEdge(graph, i, j);
                        if (trial % 2 == 0) addEdge(graph, i, j);
                    }
                    if (trial % 2 == 0) graph[i].push_back(i);
                    std::shuffle(graph[i].begin(), graph[i].end(), rng);
                }
                compareWithOracle(graph);
                ++count;
            }
        }
    }
    return count;
}

static void checkInvalidIndices() {
    for (int invalid : {-1, 1}) {
        bool rejected = false;
        try {
            smsd::detail::maximumMatching(Graph{{invalid}});
        } catch (const std::out_of_range&) {
            rejected = true;
        }
        if (!rejected) throw std::runtime_error("invalid vertex index was accepted");
    }
}

int main() {
    try {
        checkExamples();
        const auto exhaustive = checkAllSmallGraphs();
        const auto randomized = checkRandomGraphs();
        checkInvalidIndices();
        std::cout << "General matching passed: " << exhaustive
                  << " exhaustive graphs, " << randomized
                  << " seeded random graphs, blossom/azulene/Petersen examples.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "General matching failed: " << error.what() << '\n';
        return 1;
    }
}

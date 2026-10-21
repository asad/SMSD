/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Synthetic native primitive benchmarks. These do not measure MCS or
 * molecular substructure throughput, whole-application speed or GPUs: the
 * unbalanced matrices, complete graphs, long paths and cold targets expose
 * particular set-up and augmentation costs. Input construction is outside
 * the timer; validation, allocation, solving and checksums are inside it.
 *
 * Build from the repository root and select one workload:
 *   c++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_primitives.cpp \
 *       -o build/benchmark_primitives
 *   build/benchmark_primitives assignment [repeat-multiplier 1..1000]
 *   build/benchmark_primitives primitives
 *   build/benchmark_primitives setup
 *
 * assignment  Rectangular assignment and general matching. The multiplier
 *             repeats each measurement; the checksums must equal the recorded
 *             values scaled by that multiplier.
 * primitives  NLF3 construction on long paths and a complete-subgraph MCS.
 * setup       Substructure cold-target domain initialisation.
 *
 * Compare identical selectors, multipliers and compiler flags on the same
 * machine and require identical checksums.
 */
#include "smsd/smsd.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

long long microsecondsSince(Clock::time_point start) {
    return std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - start).count();
}

void benchmarkAssignments(int multiplier) {
    std::mt19937 random(91724);
    const double expectedChecksums[] = {820, 1380, 176220};
    int example = 0;
    for (const auto dims : {std::pair<int, int>{8, 512}, {512, 8}, {128, 128}}) {
        std::vector<std::vector<double>> cost(dims.first, std::vector<double>(dims.second));
        for (auto& row : cost)
            for (double& value : row) value = random() % 10000;
        const int repeats = 10 * multiplier;
        double checksum = 0;
        const auto start = Clock::now();
        for (int trial = 0; trial < repeats; ++trial)
            checksum += smsd::optimalAssign(cost).totalCost;
        const auto elapsed = microsecondsSince(start);
        if (checksum != expectedChecksums[example++] * multiplier)
            throw std::runtime_error("assignment objective checksum differs from baseline");
        std::cout << "primitive=assignment rows=" << dims.first << " cols=" << dims.second
                  << " repeats=" << repeats << " elapsed_us=" << elapsed
                  << " total_cost_checksum=" << checksum << '\n';
    }
}

void benchmarkMatchings(int multiplier) {
    for (int n : {128, 512}) {
        std::vector<std::vector<int>> graph(n);
        for (int i = 0; i < n; ++i) for (int j = i + 1; j < n; ++j) {
            graph[i].push_back(j);
            graph[j].push_back(i);
        }
        const int repeats = 5 * multiplier;
        int matched = 0;
        const auto start = Clock::now();
        for (int trial = 0; trial < repeats; ++trial) {
            const auto mate = smsd::detail::maximumMatching(graph);
            for (int v : mate) matched += v != -1;
        }
        const auto elapsed = microsecondsSince(start);
        if (matched != n * repeats)
            throw std::runtime_error("complete graph matching cardinality is incorrect");
        std::cout << "primitive=general_matching graph=complete vertices=" << n
                  << " repeats=" << repeats << " elapsed_us=" << elapsed
                  << " matched_vertex_checksum=" << matched << '\n';
    }
}

int runAssignment(int argc, char** argv) {
    int multiplier = 1;
    if (argc > 3) throw std::invalid_argument("usage: benchmark_primitives assignment [repeat-multiplier 1..1000]");
    if (argc == 3) {
        std::size_t consumed = 0;
        multiplier = std::stoi(argv[2], &consumed);
        if (consumed != std::string(argv[2]).size() || multiplier < 1 || multiplier > 1000)
            throw std::invalid_argument("repeat multiplier must be an integer in [1,1000]");
    }
    benchmarkAssignments(multiplier);
    benchmarkMatchings(multiplier);
    return 0;
}

int runPrimitives() {
    for (int n : {100, 1000}) {
        std::vector<std::vector<int>> neighbors(n), orders(n);
        for (int i = 1; i < n; ++i) {
            neighbors[i - 1].push_back(i);
            neighbors[i].push_back(i - 1);
            orders[i - 1].push_back(1);
            orders[i].push_back(1);
        }
        auto graph = smsd::MolGraph::Builder().atomCount(n)
            .atomicNumbers(std::vector<int>(n, 6))
            .setNeighbors(neighbors).setBondOrders(orders).build();
        long checksum = 0;
        const int repeats = n == 100 ? 500 : 50;
        const auto start = Clock::now();
        for (int repeat = 0; repeat < repeats; ++repeat)
            for (int i = 0; i < n; ++i)
                checksum += smsd::MolGraph::buildNLF3(graph, i).size();
        std::cout << "NLF3 n=" << n << " constructions=" << repeats * n
                  << " us=" << microsecondsSince(start) << " checksum=" << checksum << '\n';
    }

    std::mt19937 random(216491);
    std::vector<std::pair<int, int>> edges;
    constexpr int n = 80;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (random() % 100 < 55) edges.emplace_back(i, j);
    std::vector<std::vector<int>> neighbors(n), completeNeighbors(n);
    for (const auto& edge : edges) {
        neighbors[edge.first].push_back(edge.second); neighbors[edge.second].push_back(edge.first);
    }
    for (int first = 0; first < n; ++first) for (int second = 0; second < n; ++second)
        if (first != second) completeNeighbors[first].push_back(second);
    const auto graph = smsd::MolGraph::Builder().atomCount(n).atomicNumbers(std::vector<int>(n, 6)).setNeighbors(neighbors).build(false);
    const auto complete = smsd::MolGraph::Builder().atomCount(n).atomicNumbers(std::vector<int>(n, 6)).setNeighbors(completeNeighbors).build(false);
    const auto query = smsd::MCS::find(complete).withTimeout(10000);
    const auto start = Clock::now();
    int checksum = 0;
    for (int repeat = 0; repeat < 20; ++repeat) {
        const auto result = query.match(graph);
        if (!result.optimal) return 2;
        checksum += static_cast<int>(result.score);
    }
    std::cout << "complete-subgraph n=80 repeats=20 us=" << microsecondsSince(start)
              << " checksum=" << checksum << '\n';
    return 0;
}

smsd::MolGraph cycle(int n) {
    std::vector<std::vector<int>> neighbors(n), orders(n);
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        neighbors[i].push_back(j); neighbors[j].push_back(i);
        orders[i].push_back(1); orders[j].push_back(1);
    }
    return smsd::MolGraph::Builder().atomCount(n)
        .atomicNumbers(std::vector<int>(n, 6))
        .setNeighbors(neighbors).setBondOrders(orders).build(false);
}

int runSetup() {
    const auto query = cycle(24), target = cycle(128);
    smsd::ChemOptions c;
    int checksum = 0;
    const auto start = Clock::now();
    for (int i = 0; i < 20; ++i) {
        auto coldTarget = target;
        smsd::SubstructureOptions options; options.chemistry = c;
        smsd::search_detail::Budget budget(options, query.n, coldTarget.n, false);
        smsd::search_detail::Collector collector(nullptr, 1, false);
        const smsd::search_detail::QueryView view(query);
        smsd::search_detail::VertexSearch matcher(view, coldTarget, c, nullptr, false, false, budget, collector);
        checksum += query.n * coldTarget.n;
    }
    std::cout << "Substructure cold-target setup repeats=20 us=" << microsecondsSince(start)
              << " pair_checksum=" << checksum << '\n';
    return checksum == 61440 ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::string selector = argc > 1 ? argv[1] : "";
        if (selector == "assignment") return runAssignment(argc, argv);
        if (selector == "primitives" && argc == 2) return runPrimitives();
        if (selector == "setup" && argc == 2) return runSetup();
        throw std::invalid_argument(
            "usage: benchmark_primitives assignment [repeat-multiplier 1..1000] | primitives | setup");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

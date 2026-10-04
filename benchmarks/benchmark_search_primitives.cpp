/* SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Synthetic setup/search primitives; does not measure whole-application speed.
 */
#include "smsd/smsd.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

int main() {
    using Clock = std::chrono::steady_clock;
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
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - start).count();
        std::cout << "NLF3 n=" << n << " constructions=" << repeats * n
                  << " us=" << elapsed << " checksum=" << checksum << '\n';
    }

    std::mt19937 random(216491);
    std::vector<std::pair<int, int>> edges;
    constexpr int n = 80;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (random() % 100 < 55) edges.emplace_back(i, j);
    smsd::clique::ProductGraph graph;
    graph.build(std::vector<smsd::clique::ProductVertex>(n), edges);
    const auto start = Clock::now();
    int checksum = 0;
    for (int repeat = 0; repeat < 20; ++repeat) {
        const auto result = smsd::clique::findMaxCliques(graph, 200000, 10000);
        if (result.timed_out) return 2;
        checksum += result.max_size;
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(Clock::now() - start).count();
    std::cout << "clique n=80 repeats=20 us=" << elapsed << " checksum=" << checksum << '\n';
}

/* SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Measures shared Substructure domain initialisation on cold targets.
 * c++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_substructure_setup.cpp -o benchmark_setup
 */
#include "smsd/substructure.hpp"
#include <chrono>
#include <iostream>

static smsd::MolGraph cycle(int n) {
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

int main() {
    const auto query = cycle(24), target = cycle(128);
    smsd::ChemOptions c;
    int checksum = 0;
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 20; ++i) {
        auto coldTarget = target;
        smsd::SubstructureOptions options; options.chemistry = c;
        smsd::search_detail::Budget budget(options, query.n, coldTarget.n, false);
        smsd::search_detail::Collector collector(nullptr, 1, false);
        const smsd::search_detail::QueryView view(query);
        smsd::search_detail::VertexSearch matcher(view, coldTarget, c, nullptr, false, false, budget, collector);
        checksum += query.n * coldTarget.n;
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << "Substructure cold-target setup repeats=20 us=" << elapsed
              << " pair_checksum=" << checksum << '\n';
    return checksum == 61440 ? 0 : 1;
}

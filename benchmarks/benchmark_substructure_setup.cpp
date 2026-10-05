/* SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Measures cold-target VF2++ setup, not complete search throughput.
 *
 * From the repository root, build and run current headers:
 *   c++ -std=c++17 -O3 -Icpp/include \
 *       benchmarks/benchmark_substructure_setup.cpp -o /tmp/smsd_setup_after
 *   /tmp/smsd_setup_after
 *
 * Extract baseline headers to an empty directory and run the same harness:
 *   mkdir -p /tmp/smsd_algorithm_baseline
 *   git archive 49d7303 cpp/include | tar -x -C /tmp/smsd_algorithm_baseline
 *   c++ -std=c++17 -O3 -I/tmp/smsd_algorithm_baseline/cpp/include \
 *       benchmarks/benchmark_substructure_setup.cpp -o /tmp/smsd_setup_before
 *   /tmp/smsd_setup_before
 *
 * Use identical compiler flags and require domain_checksum=2560. Highly
 * symmetric cold targets expose previously unused canonicalization work;
 * prewarmed graphs and typical molecular search need separate measurements.
 */
#include "smsd/vf2pp.hpp"
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
        smsd::detail::VF2PPMatcher matcher(query, coldTarget, c, 10000);
        checksum += matcher.domainSupport(0);
    }
    const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << "VF2++ cold-target setup repeats=20 us=" << elapsed
              << " domain_checksum=" << checksum << '\n';
    return checksum == 2560 ? 0 : 1;
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 *
 * Synthetic assignment/matching primitive benchmarks. These do not measure
 * MCS, molecular substructure throughput, whole-application speed, or GPUs.
 * Input construction is outside the timer; validation, solver allocations,
 * solving, and checksum accumulation are inside it. Defaults reproduce the
 * original 10-assignment/5-matching measurements. An optional integer argument
 * multiplies repetitions for less noisy timing, e.g. ./benchmark 3.
 *
 * Build current headers from the repository root:
 *   c++ -std=c++17 -O2 -Icpp/include \
 *       benchmarks/benchmark_assignment_matching.cpp -o /tmp/smsd_assign_match
 *   /tmp/smsd_assign_match
 *
 * Build the identical harness against baseline headers (49d7303):
 *   git show 49d7303:cpp/include/smsd/hungarian.hpp > /tmp/smsd_assign_before.hpp
 *   git show 49d7303:cpp/include/smsd/general_matching.hpp > /tmp/smsd_match_before.hpp
 *   c++ -std=c++17 -O2 -Icpp/include \
 *       -DSMSD_ASSIGN_HEADER='"/tmp/smsd_assign_before.hpp"' \
 *       -DSMSD_MATCH_HEADER='"/tmp/smsd_match_before.hpp"' \
 *       benchmarks/benchmark_assignment_matching.cpp -o /tmp/smsd_assign_match_before
 *   /tmp/smsd_assign_match_before
 *
 * Compare identical repetition counts/compiler flags on the same machine and
 * require identical checksums. The intentionally unbalanced matrices and
 * complete graphs expose these particular setup/augmentation costs; their
 * speed ratios are not representative chemistry-workload claims.
 */
#ifndef SMSD_ASSIGN_HEADER
#define SMSD_ASSIGN_HEADER "smsd/hungarian.hpp"
#endif
#ifndef SMSD_MATCH_HEADER
#define SMSD_MATCH_HEADER "smsd/general_matching.hpp"
#endif
#include SMSD_ASSIGN_HEADER
#include SMSD_MATCH_HEADER

#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;

static void benchmarkAssignments(int multiplier) {
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
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
            Clock::now() - start).count();
        if (checksum != expectedChecksums[example++] * multiplier)
            throw std::runtime_error("assignment objective checksum differs from baseline");
        std::cout << "primitive=assignment rows=" << dims.first << " cols=" << dims.second
                  << " repeats=" << repeats << " elapsed_us=" << elapsed
                  << " total_cost_checksum=" << checksum << '\n';
    }
}

static void benchmarkMatchings(int multiplier) {
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
        const auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
            Clock::now() - start).count();
        if (matched != n * repeats)
            throw std::runtime_error("complete graph matching cardinality is incorrect");
        std::cout << "primitive=general_matching graph=complete vertices=" << n
                  << " repeats=" << repeats << " elapsed_us=" << elapsed
                  << " matched_vertex_checksum=" << matched << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        int multiplier = 1;
        if (argc > 2) throw std::invalid_argument("usage: benchmark [repeat-multiplier 1..1000]");
        if (argc == 2) {
            std::size_t consumed = 0;
            multiplier = std::stoi(argv[1], &consumed);
            if (consumed != std::string(argv[1]).size() || multiplier < 1 || multiplier > 1000)
                throw std::invalid_argument("repeat multiplier must be an integer in [1,1000]");
        }
        benchmarkAssignments(multiplier);
        benchmarkMatchings(multiplier);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

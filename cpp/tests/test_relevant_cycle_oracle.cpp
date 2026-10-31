/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/ring_finder.hpp"
#include <array>
#include <functional>
#include <map>
#include <random>
#include <set>

namespace {
uint64_t reduce(uint64_t bits, const std::array<uint64_t, 64>& rows) {
    for (int pivot = 63; pivot >= 0; --pivot)
        if ((bits & (uint64_t(1) << pivot)) != 0) {
            if (rows[pivot] == 0) return bits;
            bits ^= rows[pivot];
        }
    return 0;
}
void insert(uint64_t bits, std::array<uint64_t, 64>& rows) {
    bits = reduce(bits, rows);
    for (int pivot = 63; pivot >= 0; --pivot)
        if ((bits & (uint64_t(1) << pivot)) != 0) { rows[pivot] = bits; return; }
}
uint64_t cycleMask(const std::vector<int>& cycle, const std::vector<std::vector<int>>& edgeIndex) {
    uint64_t bits = 0;
    for (size_t atom = 0; atom < cycle.size(); ++atom) {
        const int edge = edgeIndex[cycle[atom]][cycle[(atom + 1) % cycle.size()]];
        smsd_test::require(edge >= 0, "Every perceived cycle bond belongs to the input graph");
        bits |= uint64_t(1) << edge;
    }
    smsd_test::require(smsd::popcount64(bits) == static_cast<int>(cycle.size()),
                       "Perceived cycles have distinct bonds");
    return bits;
}
}

int main() {
    std::mt19937 random(0x5eedc0deu);
    for (int sample = 0; sample < 128; ++sample) {
        const int atoms = 3 + sample % 6;
        std::vector<std::vector<int>> neighbors(atoms), orders(atoms), edgeIndex(atoms, std::vector<int>(atoms, -1));
        int edges = 0;
        for (int first = 0; first < atoms; ++first)
            for (int second = first + 1; second < atoms; ++second)
                if (random() % 8 < 3) {
                    neighbors[first].push_back(second); neighbors[second].push_back(first);
                    orders[first].push_back(1); orders[second].push_back(1);
                    edgeIndex[first][second] = edgeIndex[second][first] = edges++;
                }
        std::set<uint64_t> all;
        std::vector<uint8_t> used(atoms, 0);
        std::vector<int> path;
        std::function<void(int, int)> enumerate = [&](int root, int atom) {
            for (int next : neighbors[atom]) {
                if (next == root && path.size() >= 3) all.insert(cycleMask(path, edgeIndex));
                else if (next > root && !used[next]) {
                    used[next] = 1; path.push_back(next);
                    enumerate(root, next);
                    path.pop_back(); used[next] = 0;
                }
            }
        };
        int components = 0;
        std::vector<uint8_t> reached(atoms, 0);
        for (int root = 0; root < atoms; ++root) {
            used[root] = 1; path.assign(1, root); enumerate(root, root); used[root] = 0;
            if (reached[root]) continue;
            ++components;
            std::vector<int> queue{root}; reached[root] = 1;
            for (size_t head = 0; head < queue.size(); ++head)
                for (int next : neighbors[queue[head]])
                    if (!reached[next]) { reached[next] = 1; queue.push_back(next); }
        }
        std::map<int, std::vector<uint64_t>> byLength;
        for (uint64_t cycle : all) byLength[smsd::popcount64(cycle)].push_back(cycle);
        std::array<uint64_t, 64> rows{};
        std::set<uint64_t> expected;
        size_t basisSize = 0, minimumWeight = 0;
        for (const auto& group : byLength) {
            // The first pass sees only strictly shorter cycles. The second extends their span.
            for (uint64_t cycle : group.second) if (reduce(cycle, rows) != 0) expected.insert(cycle);
            for (uint64_t cycle : group.second)
                if (reduce(cycle, rows) != 0) {
                    insert(cycle, rows); ++basisSize; minimumWeight += group.first;
                }
        }
        smsd_test::require(basisSize == static_cast<size_t>(edges - atoms + components),
                           "The independent oracle spans the cycle space");
        const auto graph = smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 6))
            .setNeighbors(neighbors).setBondOrders(orders).build(false);
        smsd_test::require(!smsd::ringsTruncated(graph) && !smsd::relevantCyclesTruncated(graph),
                           "Tiny exhaustive controls fit both ring allowances");
        std::set<uint64_t> actual;
        for (const auto& cycle : smsd::computeRelevantCycles(graph))
            smsd_test::require(actual.insert(cycleMask(cycle, edgeIndex)).second,
                               "Relevant cycles are not repeated");
        smsd_test::require(actual == expected, "Relevant cycles equal the independent complete simple-cycle oracle");
        const auto basis = smsd::computeRings(graph);
        size_t weight = 0;
        std::array<uint64_t, 64> actualRows{};
        for (const auto& cycle : basis) {
            const uint64_t bits = cycleMask(cycle, edgeIndex);
            smsd_test::require(reduce(bits, actualRows) != 0, "The returned basis is independent");
            insert(bits, actualRows); weight += cycle.size();
        }
        smsd_test::require(basis.size() == basisSize && weight == minimumWeight,
                           "The returned basis has full rank and minimum total length");
    }
    std::cout << "PASS bounded exhaustive relevant-cycle oracle\n";
    return 0;
}

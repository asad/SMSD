/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Rectangular assignment regressions with an exhaustive injection oracle.
 */
#include "smsd/hungarian.hpp"

#include <functional>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>

using Costs = std::vector<std::vector<double>>;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static double assignmentOracle(const Costs& cost) {
    const int rows = static_cast<int>(cost.size());
    const int cols = static_cast<int>(cost[0].size());
    const bool transpose = rows > cols;
    const int small = std::min(rows, cols), large = std::max(rows, cols);
    double best = std::numeric_limits<double>::infinity();
    std::vector<bool> used(large, false);
    std::function<void(int, double)> enumerate = [&](int a, double total) {
        if (a == small) {
            best = std::min(best, total);
            return;
        }
        for (int b = 0; b < large; ++b) {
            if (used[b]) continue;
            used[b] = true;
            enumerate(a + 1, total + (transpose ? cost[b][a] : cost[a][b]));
            used[b] = false;
        }
    };
    enumerate(0, 0);
    return best;
}

static void checkCost(const Costs& cost, double penalty) {
    const auto found = smsd::optimalAssign(cost, penalty);
    const std::size_t rows = cost.size(), cols = cost[0].size();
    require(found.assignment.size() == std::min(rows, cols), "assignment cardinality is incorrect");
    require(std::is_sorted(found.assignment.begin(), found.assignment.end()), "assignment is not deterministic by row");
    std::set<int> usedRows, usedCols;
    double sum = 0;
    for (const auto& [row, col] : found.assignment) {
        require(row >= 0 && static_cast<std::size_t>(row) < rows
                && col >= 0 && static_cast<std::size_t>(col) < cols, "assignment index is out of range");
        require(usedRows.insert(row).second && usedCols.insert(col).second, "assignment reuses a row or column");
        sum += cost[row][col];
    }
    require(sum == found.totalCost, "reported total differs from assigned real costs");
    require(sum == assignmentOracle(cost), "assignment differs from exhaustive injection oracle");
}

static void testFiniteCosts() {
    require(smsd::optimalAssign({}).assignment.empty(), "empty matrix should have an empty assignment");
    require(smsd::optimalAssign(Costs(3)).assignment.empty(), "zero-column matrix should have an empty assignment");
    std::mt19937 random(0xA5512026);
    std::size_t count = 0;
    for (int rows = 1; rows <= 5; ++rows) {
        for (int cols = 1; cols <= 5; ++cols) {
            for (int trial = 0; trial < 128; ++trial) {
                Costs cost(rows, std::vector<double>(cols));
                for (auto& row : cost)
                    for (double& value : row) value = static_cast<int>(random() % 21) - 10;
                for (double penalty : {-13.0, 0.0, 7.0}) checkCost(cost, penalty);
                ++count;
            }
        }
    }
    checkCost({{7, -3, 1, 8, -5}}, 3);
    checkCost({{7}, {-3}, {1}, {8}, {-5}}, 3);
    checkCost({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, 0);
    std::cout << "Assignment oracle passed: " << count << " matrices across three penalties.\n";
}

static void testInvalidInputs() {
    for (const auto& cost : {Costs{{1, 2}, {3}}, Costs{{}, {1}},
                            Costs{{std::numeric_limits<double>::infinity()}},
                            Costs{{std::numeric_limits<double>::quiet_NaN()}}}) {
        bool rejected = false;
        try {
            smsd::optimalAssign(cost);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        require(rejected, "ragged/nonfinite cost input was not rejected");
    }
    bool rejectedPenalty = false;
    try {
        smsd::optimalAssign({{1, 2}}, std::numeric_limits<double>::infinity());
    } catch (const std::invalid_argument&) {
        rejectedPenalty = true;
    }
    require(rejectedPenalty, "nonfinite unmatched penalty was not rejected");
    bool rejectedOverflow = false;
    try {
        const double maximum = std::numeric_limits<double>::max();
        smsd::optimalAssign({{-maximum, maximum}, {-maximum, maximum}});
    } catch (const std::overflow_error&) {
        rejectedOverflow = true;
    }
    require(rejectedOverflow, "reduced-cost overflow should fail instead of looping");
    std::cout << "Assignment invalid-input regressions passed.\n";
}

int main(int argc, char** argv) {
    try {
        const std::string category = argc > 1 ? argv[1] : "all";
        if (category == "all" || category == "finite") testFiniteCosts();
        if (category == "all" || category == "invalid") testInvalidInputs();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Assignment regression failed: " << error.what() << '\n';
        return 1;
    }
}

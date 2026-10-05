/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// Optimal assignment solver
//
// Computes a minimum-cost assignment for a rectangular m x n cost matrix.
// Assigns every vertex of the smaller side. Uniform unmatched penalties
// contribute the same constant to every assignment and need no dummy matrix.
//
// Time complexity : O(min(m,n)^2 * max(m,n))
// Auxiliary space: O(m+n)
// ============================================================================

struct AssignmentResult {
    std::vector<std::pair<int,int>> assignment;  // (row, col) pairs
    double totalCost;
};

/**
 * Solve the rectangular assignment problem.
 *
 * @param cost     m x n cost matrix (may be rectangular)
 * @param penalty  cost for unmatched rows/columns in rectangular problems (default 1.0)
 * @return         Optimal (row, col) pairs and the sum of their real costs;
 *                 totalCost excludes unmatched penalties, as in earlier releases
 * @throws         std::invalid_argument for ragged matrices or nonfinite costs/penalty
 * @throws         std::overflow_error if reduced-cost arithmetic overflows
 */
inline AssignmentResult optimalAssign(
        const std::vector<std::vector<double>>& cost,
        double penalty = 1.0) {

    if (cost.empty()) return {{}, 0.0};
    const std::size_t columns = cost[0].size();
    if (!std::isfinite(penalty))
        throw std::invalid_argument("optimalAssign: penalty must be finite");
    for (const auto& row : cost) {
        if (row.size() != columns)
            throw std::invalid_argument("optimalAssign: matrix must be rectangular");
        for (double value : row)
            if (!std::isfinite(value))
                throw std::invalid_argument("optimalAssign: costs must be finite");
    }
    if (columns == 0) return {{}, 0.0};
    const auto indexLimit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (cost.size() >= indexLimit || columns >= indexLimit)
        throw std::length_error("optimalAssign: matrix dimensions exceed integer indexing");

    const int m = static_cast<int>(cost.size()), n = static_cast<int>(columns);
    const bool transpose = m > n;
    const int rows = std::min(m, n), cols = std::max(m, n);
    auto realCost = [&](int i, int j) {
        return transpose ? cost[j][i] : cost[i][j];
    };

    // u[i] = potential for row i, v[j] = potential for column j
    std::vector<double> u(rows + 1, 0.0), v(cols + 1, 0.0);
    // p[j] = row assigned to column j (1-indexed, 0 = unassigned)
    std::vector<int> p(cols + 1, 0), way(cols + 1, 0);
    std::vector<double> minv(cols + 1);
    std::vector<bool> used(cols + 1);

    for (int i = 1; i <= rows; i++) {
        p[0] = i;
        int j0 = 0;
        std::fill(minv.begin(), minv.end(), std::numeric_limits<double>::infinity());
        std::fill(used.begin(), used.end(), false);

        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            double delta = std::numeric_limits<double>::infinity();

            for (int j = 1; j <= cols; j++) {
                if (used[j]) continue;
                double cur = realCost(i0 - 1, j - 1) - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }

            // Finite input can still overflow reduced-cost arithmetic. Fail
            // explicitly instead of repeatedly choosing the sentinel column.
            if (!std::isfinite(delta) || j1 == 0)
                throw std::overflow_error("optimalAssign: reduced cost overflow");
            for (int j = 0; j <= cols; j++) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
        } while (p[j0] != 0);

        // Augmenting path
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }

    // Extract assignment (only real rows and columns)
    AssignmentResult result;
    result.totalCost = 0.0;
    result.assignment.reserve(rows);
    for (int j = 1; j <= cols; j++) {
        if (p[j] == 0) continue;
        const int row = transpose ? j - 1 : p[j] - 1;
        const int col = transpose ? p[j] - 1 : j - 1;
        result.assignment.emplace_back(row, col);
        result.totalCost += cost[row][col];
    }

    // Sort by row for deterministic output
    std::sort(result.assignment.begin(), result.assignment.end());
    return result;
}

// Backward-compatible aliases
using HungarianResult = AssignmentResult;
inline AssignmentResult hungarianSolve(
        const std::vector<std::vector<double>>& cost,
        double penalty = 1.0) {
    return optimalAssign(cost, penalty);
}

} // namespace smsd

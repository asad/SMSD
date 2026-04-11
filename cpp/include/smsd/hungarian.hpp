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
#include <vector>

namespace smsd {

// ============================================================================
// Optimal assignment solver
//
// Computes a minimum-cost assignment for a rectangular m x n cost matrix.
// For unbalanced problems (m != n), virtual rows or columns with penalty
// cost are added to make the problem square.
//
// Time complexity : O(max(m,n)^3)
// Space complexity: O(max(m,n)^2)
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
 * @return         AssignmentResult with optimal (row, col) assignment pairs and total cost
 */
inline AssignmentResult optimalAssign(
        const std::vector<std::vector<double>>& cost,
        double penalty = 1.0) {

    int m = static_cast<int>(cost.size());
    if (m == 0) return {{}, 0.0};
    int n = static_cast<int>(cost[0].size());
    if (n == 0) return {{}, 0.0};

    // Pad to square matrix
    int sz = std::max(m, n);
    std::vector<std::vector<double>> c(sz, std::vector<double>(sz, penalty));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            c[i][j] = cost[i][j];

    // u[i] = potential for row i, v[j] = potential for column j
    std::vector<double> u(sz + 1, 0.0), v(sz + 1, 0.0);
    // p[j] = row assigned to column j (1-indexed, 0 = unassigned)
    std::vector<int> p(sz + 1, 0), way(sz + 1, 0);

    for (int i = 1; i <= sz; i++) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(sz + 1, std::numeric_limits<double>::infinity());
        std::vector<bool> used(sz + 1, false);

        do {
            used[j0] = true;
            int i0 = p[j0], j1 = 0;
            double delta = std::numeric_limits<double>::infinity();

            for (int j = 1; j <= sz; j++) {
                if (used[j]) continue;
                double cur = c[i0 - 1][j - 1] - u[i0] - v[j];
                if (cur < minv[j]) {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if (minv[j] < delta) {
                    delta = minv[j];
                    j1 = j;
                }
            }

            for (int j = 0; j <= sz; j++) {
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
    for (int j = 1; j <= sz; j++) {
        int row = p[j] - 1;
        int col = j - 1;
        if (row < m && col < n) {
            result.assignment.emplace_back(row, col);
            result.totalCost += cost[row][col];
        }
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

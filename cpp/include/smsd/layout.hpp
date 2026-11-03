/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only C++17 2D layout post-processing utilities.
 * Provides:
 *   - Bond crossing reduction via simulated annealing on ring orientations.
 *   - Force-directed layout minimisation (attractive/repulsive + crossing penalty).
 *   - Stress majorisation (SMACOF) for graph-distance-proportional layouts.
 *   - Scaffold template library for common ring systems.
 *
 * Zero external dependencies -- only requires smsd/mol_graph.hpp and
 * smsd/ring_finder.hpp.
 */
#pragma once
#ifndef SMSD_LAYOUT_HPP
#define SMSD_LAYOUT_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/ring_finder.hpp"

// MSVC does not define M_PI by default
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace smsd {

// ── 2D point for layout coordinates ─────────────────────────────────────────

struct Point2D {
    double x, y;
};

inline std::vector<Point2D> matchTemplate(const MolGraph& g,
                                          double targetBondLength);

namespace detail_layout {

inline void normalizeCoords(const MolGraph& g, std::vector<Point2D>& coords,
                            double spacing = 1.5) {
    if (coords.size() < static_cast<size_t>(g.n)) {
        size_t oldSize = coords.size();
        coords.resize(static_cast<size_t>(g.n));
        for (size_t i = oldSize; i < coords.size(); ++i) {
            coords[i].x = static_cast<double>(i) * spacing;
            coords[i].y = 0.0;
        }
    } else if (coords.size() > static_cast<size_t>(g.n)) {
        coords.resize(static_cast<size_t>(g.n));
    }
}

inline bool isDegenerateLayout(const std::vector<Point2D>& coords) {
    if (coords.size() < 3) return true;
    double minX = coords[0].x, maxX = coords[0].x;
    double minY = coords[0].y, maxY = coords[0].y;
    for (const auto& p : coords) {
        minX = std::min(minX, p.x);
        maxX = std::max(maxX, p.x);
        minY = std::min(minY, p.y);
        maxY = std::max(maxY, p.y);
    }
    return (maxX - minX < 1e-6) || (maxY - minY < 1e-6);
}

// --------------------------------------------------------------------------
// Line-segment intersection test (for bond crossing detection)
// --------------------------------------------------------------------------

/// Test whether segments (p1,p2) and (p3,p4) properly cross.
/// Returns true if the interior of the two segments intersect (not just touch).
inline bool segmentsCross(const Point2D& p1, const Point2D& p2,
                          const Point2D& p3, const Point2D& p4) {
    auto cross2d = [](double ax, double ay, double bx, double by) {
        return ax * by - ay * bx;
    };
    double d1x = p2.x - p1.x, d1y = p2.y - p1.y;
    double d2x = p4.x - p3.x, d2y = p4.y - p3.y;

    double denom = cross2d(d1x, d1y, d2x, d2y);
    if (std::abs(denom) < 1e-12) return false; // parallel or collinear

    double dx = p3.x - p1.x, dy = p3.y - p1.y;
    double t = cross2d(dx, dy, d2x, d2y) / denom;
    double u = cross2d(dx, dy, d1x, d1y) / denom;

    // Strict interior intersection (exclude endpoints to avoid counting
    // shared-vertex bonds as crossings)
    constexpr double eps = 1e-9;
    return (t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps);
}

// --------------------------------------------------------------------------
// Count total bond crossings in the current coordinate set
// --------------------------------------------------------------------------

inline int countCrossings(const MolGraph& g, const std::vector<Point2D>& coords) {
    // Collect all bonds
    struct Bond { int a, b; };
    std::vector<Bond> bonds;
    bonds.reserve(g.n * 2);
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) bonds.push_back({i, j});
        }
    }

    int crossings = 0;
    int nb = static_cast<int>(bonds.size());
    for (int i = 0; i < nb; i++) {
        for (int j = i + 1; j < nb; j++) {
            int a1 = bonds[i].a, a2 = bonds[i].b;
            int b1 = bonds[j].a, b2 = bonds[j].b;
            // Skip bonds sharing a vertex
            if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
            if (segmentsCross(coords[a1], coords[a2], coords[b1], coords[b2]))
                crossings++;
        }
    }
    return crossings;
}

// --------------------------------------------------------------------------
// Identify ring systems (connected components of ring adjacency)
// --------------------------------------------------------------------------

struct RingSystem {
    std::vector<std::vector<int>> rings; // each ring = list of atom indices
    std::set<int> atoms;                 // all atoms in this system
};

// An individual ring within a multi-ring system, split into exclusive
// (flippable) atoms and shared (pivot) atoms.
struct FlippableRing {
    std::vector<int> exclusive; // atoms only in this ring — get mirrored
    std::vector<int> shared;    // atoms shared with other rings — stay fixed
};

inline std::vector<RingSystem> findRingSystems(const MolGraph& g) {
    auto sssr = computeSSSR(g);
    if (sssr.empty()) return {};

    const auto groups = detail_rf::groupRingsBySharedBond(sssr);
    std::vector<RingSystem> systems(groups.systemCount);
    for (size_t i = 0; i < sssr.size(); i++) {
        systems[groups.system[i]].rings.push_back(sssr[i]);
        for (int a : sssr[i]) systems[groups.system[i]].atoms.insert(a);
    }
    return systems;
}

} // namespace detail_layout


// ==========================================================================
// reduceCrossings — bond crossing reduction via simulated annealing
// ==========================================================================

/**
 * Reduce bond crossings in a 2D molecular layout via two-phase
 * simulated annealing.
 *
 * Phase 1 — System-level flipping:
 *   Flip entire ring systems (connected components of fused rings)
 *   about their centroid to find the best global orientation.
 *
 * Phase 2 — Individual ring flipping:
 *   Within multi-ring systems, flip individual SSSR rings while
 *   keeping shared (fusion) atoms fixed as pivots.  This resolves
 *   crossings that system-level moves cannot reach.
 *
 * Both phases use Metropolis criterion with exponential cooling.
 *
 * @param g        The molecular graph.
 * @param coords   Input/output 2D coordinates (coords[atomIdx] = {x, y}).
 *                 Modified in-place with optimised positions.
 * @param maxIter  Maximum total SA iterations across both phases (default 1000).
 * @return         Number of bond crossings remaining after optimisation.
 */
inline int reduceCrossings(const MolGraph& g, std::vector<Point2D>& coords,
                           int maxIter = 1000) {
    detail_layout::normalizeCoords(g, coords);
    if (g.n < 4) return 0;

    auto systems = detail_layout::findRingSystems(g);
    if (systems.empty()) return detail_layout::countCrossings(g, coords);

    int nSys = static_cast<int>(systems.size());
    int currentCrossings = detail_layout::countCrossings(g, coords);
    if (currentCrossings == 0) return 0;

    // Simulated annealing parameters (shared across both phases)
    double T0 = 2.0;          // initial temperature
    double Tmin = 0.01;       // minimum temperature
    double alpha = 0.995;      // cooling rate

    std::mt19937 rng(42);      // deterministic seed for reproducibility
    std::uniform_real_distribution<double> probDist(0.0, 1.0);
    std::uniform_int_distribution<int> axisDist(0, 1);

    double T = T0;

    // ---- Phase 1: system-level flipping ----
    int phase1Iters = maxIter / 2;
    {
        std::uniform_int_distribution<int> sysDist(0, nSys - 1);
        for (int iter = 0; iter < phase1Iters && currentCrossings > 0; iter++) {
            int si = sysDist(rng);
            auto& sys = systems[si];
            if (sys.atoms.size() < 3) continue;

            double cx = 0, cy = 0;
            for (int a : sys.atoms) { cx += coords[a].x; cy += coords[a].y; }
            cx /= static_cast<double>(sys.atoms.size());
            cy /= static_cast<double>(sys.atoms.size());

            std::vector<std::pair<int, Point2D>> saved;
            saved.reserve(sys.atoms.size());
            for (int a : sys.atoms) saved.push_back({a, coords[a]});

            int axis = axisDist(rng);
            for (int a : sys.atoms) {
                if (axis == 0) coords[a].x = 2.0 * cx - coords[a].x;
                else           coords[a].y = 2.0 * cy - coords[a].y;
            }

            int newCrossings = detail_layout::countCrossings(g, coords);
            int delta = newCrossings - currentCrossings;

            bool accept = (delta <= 0);
            if (!accept && T > Tmin) {
                accept = (probDist(rng) < std::exp(-static_cast<double>(delta) / T));
            }

            if (accept) { currentCrossings = newCrossings; }
            else { for (auto& [a, p] : saved) coords[a] = p; }

            T *= alpha;
            if (T < Tmin) T = Tmin;
        }
    }

    if (currentCrossings == 0) return 0;

    // ---- Phase 2: individual ring flipping within multi-ring systems ----
    // Build flippable rings: for each ring in a multi-ring system, identify
    // exclusive atoms (only in this ring) and shared atoms (in ≥2 rings).
    std::vector<detail_layout::FlippableRing> flippable;
    for (auto& sys : systems) {
        if (sys.rings.size() < 2) continue; // single ring — handled in phase 1

        // Count how many rings each atom belongs to within this system
        std::unordered_map<int, int> atomRingCount;
        for (auto& ring : sys.rings) {
            for (int a : ring) atomRingCount[a]++;
        }

        for (auto& ring : sys.rings) {
            detail_layout::FlippableRing fr;
            for (int a : ring) {
                if (atomRingCount[a] >= 2) fr.shared.push_back(a);
                else                       fr.exclusive.push_back(a);
            }
            // Need both exclusive and shared atoms for a meaningful flip
            if (!fr.exclusive.empty() && !fr.shared.empty()) {
                flippable.push_back(std::move(fr));
            }
        }
    }

    if (!flippable.empty()) {
        int nFlip = static_cast<int>(flippable.size());
        std::uniform_int_distribution<int> flipDist(0, nFlip - 1);

        // Reset temperature for phase 2
        T = T0;
        int phase2Iters = maxIter - phase1Iters;

        for (int iter = 0; iter < phase2Iters && currentCrossings > 0; iter++) {
            int fi = flipDist(rng);
            auto& fr = flippable[fi];

            // Pivot = centroid of shared (fusion) atoms
            double cx = 0, cy = 0;
            for (int a : fr.shared) { cx += coords[a].x; cy += coords[a].y; }
            cx /= static_cast<double>(fr.shared.size());
            cy /= static_cast<double>(fr.shared.size());

            // Save exclusive atom coordinates
            std::vector<std::pair<int, Point2D>> saved;
            saved.reserve(fr.exclusive.size());
            for (int a : fr.exclusive) saved.push_back({a, coords[a]});

            // Mirror exclusive atoms about the shared-atom centroid
            int axis = axisDist(rng);
            for (int a : fr.exclusive) {
                if (axis == 0) coords[a].x = 2.0 * cx - coords[a].x;
                else           coords[a].y = 2.0 * cy - coords[a].y;
            }

            int newCrossings = detail_layout::countCrossings(g, coords);
            int delta = newCrossings - currentCrossings;

            bool accept = (delta <= 0);
            if (!accept && T > Tmin) {
                accept = (probDist(rng) < std::exp(-static_cast<double>(delta) / T));
            }

            if (accept) { currentCrossings = newCrossings; }
            else { for (auto& [a, p] : saved) coords[a] = p; }

            T *= alpha;
            if (T < Tmin) T = Tmin;
        }
    }

    return currentCrossings;
}


// ==========================================================================
// detail_layout — BFS graph-distance matrix
// ==========================================================================

namespace detail_layout {

/// Compute all-pairs shortest path distances via BFS from each atom.
inline std::vector<std::vector<int>> bfsDistanceMatrix(const MolGraph& g) {
    int n = g.n;
    std::vector<std::vector<int>> D(n, std::vector<int>(n, n + 1));
    for (int src = 0; src < n; src++) {
        D[src][src] = 0;
        std::deque<int> q;
        q.push_back(src);
        while (!q.empty()) {
            int u = q.front(); q.pop_front();
            for (int v : g.neighbors[u]) {
                if (D[src][v] > D[src][u] + 1) {
                    D[src][v] = D[src][u] + 1;
                    q.push_back(v);
                }
            }
        }
    }
    return D;
}

/// Compute segment–segment intersection point.
/// Returns true if segments (p1,p2) and (p3,p4) cross, and writes the
/// intersection point into 'out'.
inline bool segmentIntersectionPoint(const Point2D& p1, const Point2D& p2,
                                     const Point2D& p3, const Point2D& p4,
                                     Point2D& out) {
    double d1x = p2.x - p1.x, d1y = p2.y - p1.y;
    double d2x = p4.x - p3.x, d2y = p4.y - p3.y;
    double denom = d1x * d2y - d1y * d2x;
    if (std::abs(denom) < 1e-12) return false;
    double dx = p3.x - p1.x, dy = p3.y - p1.y;
    double t = (dx * d2y - dy * d2x) / denom;
    double u = (dx * d1y - dy * d1x) / denom;
    constexpr double eps = 1e-9;
    if (t <= eps || t >= 1.0 - eps || u <= eps || u >= 1.0 - eps) return false;
    out.x = p1.x + t * d1x;
    out.y = p1.y + t * d1y;
    return true;
}

} // namespace detail_layout


// ==========================================================================
// forceDirectedLayout — force-directed 2D layout minimisation
// ==========================================================================

/**
 * Force-directed 2D layout minimisation.
 * Iteratively moves atoms to minimise:
 *   1. Bond-length stress (attractive: keep bonded atoms at targetBondLength)
 *   2. Non-bonded repulsion (push overlapping atoms apart)
 *   3. Crossing penalty (push crossing bond midpoints apart)
 *
 * @param g                The molecular graph.
 * @param coords           Input/output 2D coordinates.
 * @param maxIter          Maximum number of iterations (default 500).
 * @param targetBondLength Desired bond length in coordinate units (default 1.5).
 * @return                 Final stress energy.
 */
inline double forceDirectedLayout(const MolGraph& g, std::vector<Point2D>& coords,
                                   int maxIter = 500, double targetBondLength = 1.5) {
    int n = g.n;
    if (n < 2) return 0.0;
    detail_layout::normalizeCoords(g, coords, targetBondLength);
    if (detail_layout::isDegenerateLayout(coords)) {
        auto templ = matchTemplate(g, targetBondLength);
        if (templ.size() == static_cast<size_t>(n)) coords = std::move(templ);
    }

    // Force constants
    constexpr double k_attract  = 1.0;
    constexpr double k_repel    = 0.5;
    constexpr double k_crossing = 2.0;

    // Step size with cooling
    double step = 0.1;
    constexpr double cooling     = 0.98;
    constexpr double convergence = 0.01;

    // Build adjacency lookup for O(1) bonded checks
    std::vector<std::vector<bool>> bonded(n, std::vector<bool>(n, false));
    for (int i = 0; i < n; i++)
        for (int j : g.neighbors[i])
            bonded[i][j] = true;

    // Collect bonds for crossing detection
    struct Bond { int a, b; };
    std::vector<Bond> bonds;
    bonds.reserve(n * 2);
    for (int i = 0; i < n; i++)
        for (int j : g.neighbors[i])
            if (j > i) bonds.push_back({i, j});
    int nb = static_cast<int>(bonds.size());

    double totalStress = 0.0;

    for (int iter = 0; iter < maxIter; iter++) {
        // Net forces on each atom
        std::vector<double> fx(n, 0.0), fy(n, 0.0);

        // 1. Attractive forces (bonded pairs)
        for (auto& bond : bonds) {
            int i = bond.a, j = bond.b;
            double dx = coords[j].x - coords[i].x;
            double dy = coords[j].y - coords[i].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist < 1e-10) dist = 1e-10;
            double force = k_attract * (dist - targetBondLength);
            double ux = dx / dist, uy = dy / dist;
            fx[i] += force * ux;
            fy[i] += force * uy;
            fx[j] -= force * ux;
            fy[j] -= force * uy;
        }

        // 2. Repulsive forces (all non-bonded pairs)
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (bonded[i][j]) continue;
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                double distSq = dx * dx + dy * dy;
                if (distSq < 1e-6) distSq = 1e-6;
                double force = -k_repel / distSq;
                double dist = std::sqrt(distSq);
                double ux = dx / dist, uy = dy / dist;
                fx[i] += force * ux;
                fy[i] += force * uy;
                fx[j] -= force * ux;
                fy[j] -= force * uy;
            }
        }

        // 3. Crossing penalty forces
        for (int bi = 0; bi < nb; bi++) {
            for (int bj = bi + 1; bj < nb; bj++) {
                int a1 = bonds[bi].a, a2 = bonds[bi].b;
                int b1 = bonds[bj].a, b2 = bonds[bj].b;
                if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;

                Point2D crossPt;
                if (detail_layout::segmentIntersectionPoint(
                        coords[a1], coords[a2], coords[b1], coords[b2], crossPt)) {
                    // Compute midpoints of each bond
                    double m1x = (coords[a1].x + coords[a2].x) * 0.5;
                    double m1y = (coords[a1].y + coords[a2].y) * 0.5;
                    double m2x = (coords[b1].x + coords[b2].x) * 0.5;
                    double m2y = (coords[b1].y + coords[b2].y) * 0.5;
                    // Push midpoints apart
                    double dx = m2x - m1x;
                    double dy = m2y - m1y;
                    double dist = std::sqrt(dx * dx + dy * dy);
                    if (dist < 1e-8) { dx = 0.1; dy = 0.1; dist = std::sqrt(0.02); }
                    double force = -k_crossing / (dist + 0.1);
                    double ux = dx / dist, uy = dy / dist;
                    // Distribute force to bond endpoints
                    fx[a1] += 0.5 * force * ux;
                    fy[a1] += 0.5 * force * uy;
                    fx[a2] += 0.5 * force * ux;
                    fy[a2] += 0.5 * force * uy;
                    fx[b1] -= 0.5 * force * ux;
                    fy[b1] -= 0.5 * force * uy;
                    fx[b2] -= 0.5 * force * ux;
                    fy[b2] -= 0.5 * force * uy;
                }
            }
        }

        // Update positions with PrEd-style clamping (Bertault 2000):
        // Limit each atom's displacement to avoid introducing NEW crossings.
        // If moving atom i by (dx,dy) would cause any of its bonds to cross
        // a non-crossing bond, halve the displacement iteratively.
        int crossingsBefore = detail_layout::countCrossings(g, coords);
        double maxDisp = 0.0;
        for (int i = 0; i < n; i++) {
            double dx = step * fx[i];
            double dy = step * fy[i];
            // Track original unclamped displacement for convergence
            double origDisp = std::sqrt(dx * dx + dy * dy);
            Point2D orig = coords[i];
            coords[i].x += dx;
            coords[i].y += dy;
            int crossingsAfter = detail_layout::countCrossings(g, coords);
            if (crossingsAfter > crossingsBefore) {
                // PrEd: halve displacement until no new crossings
                for (int s = 0; s < 4; s++) {
                    dx *= 0.5; dy *= 0.5;
                    coords[i] = {orig.x + dx, orig.y + dy};
                    if (detail_layout::countCrossings(g, coords) <= crossingsBefore) break;
                }
                if (detail_layout::countCrossings(g, coords) > crossingsBefore)
                    coords[i] = orig;
            }
            if (origDisp > maxDisp) maxDisp = origDisp;
        }

        // Cooling
        step *= cooling;

        // Convergence check
        if (maxDisp < convergence) break;
    }

    // Compute final stress: sum of bond-length deviations
    totalStress = 0.0;
    for (auto& bond : bonds) {
        double dx = coords[bond.b].x - coords[bond.a].x;
        double dy = coords[bond.b].y - coords[bond.a].y;
        double dist = std::sqrt(dx * dx + dy * dy);
        double delta = dist - targetBondLength;
        totalStress += delta * delta;
    }
    return totalStress;
}


// ==========================================================================
// stressMajorisation — SMACOF algorithm for 2D layout
// ==========================================================================

/**
 * Stress majorisation (SMACOF algorithm) for 2D layout.
 * Minimises: sum_{i<j} w_ij * (d_ij - D_ij)^2
 * where D_ij = graph_distance(i,j) * targetBondLength
 *
 * Uses iterative weighted averaging for guaranteed monotone convergence.
 *
 * @param g                The molecular graph.
 * @param coords           Input/output 2D coordinates.
 * @param maxIter          Maximum iterations (default 300).
 * @param targetBondLength Desired bond length (default 1.5).
 * @return                 Final normalised stress value.
 */
inline double stressMajorisation(const MolGraph& g, std::vector<Point2D>& coords,
                                  int maxIter = 300, double targetBondLength = 1.5,
                                  int nInit = 3) {
    int n = g.n;
    if (n < 2) return 0.0;
    detail_layout::normalizeCoords(g, coords, targetBondLength);
    if (detail_layout::isDegenerateLayout(coords)) {
        auto templ = matchTemplate(g, targetBondLength);
        if (templ.size() == static_cast<size_t>(n)) coords = std::move(templ);
    }

    // Compute target distances from BFS graph distances
    auto graphDist = detail_layout::bfsDistanceMatrix(g);
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            D[i][j] = static_cast<double>(graphDist[i][j]) * targetBondLength;

    // Weights: w_ij = 1 / D_ij^2 (standard stress majorisation weighting)
    std::vector<std::vector<double>> W(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j && D[i][j] > 1e-10)
                W[i][j] = 1.0 / (D[i][j] * D[i][j]);

    // Weighted degree: sum of all weights for atom i
    std::vector<double> Wsum(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            Wsum[i] += W[i][j];

    // Multiple random initialisations (SMACOF best practice, de Leeuw 1977).
    // Keep the result with lowest final stress.
    std::vector<Point2D> bestCoords = coords;
    double bestStress = std::numeric_limits<double>::max();
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> perturb(-0.5, 0.5);

    for (int init = 0; init < std::max(1, nInit); ++init) {
        // Perturb initial coordinates (except first run which uses input)
        auto runCoords = coords;
        if (init > 0) {
            for (int i = 0; i < n; i++) {
                runCoords[i].x += perturb(rng) * targetBondLength;
                runCoords[i].y += perturb(rng) * targetBondLength;
            }
        }

    // Local aliases for this run
    auto& lCoords = runCoords;
    double prevStress = [&]() -> double {
        double s = 0.0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                double dx = lCoords[i].x - lCoords[j].x;
                double dy = lCoords[i].y - lCoords[j].y;
                double dij = std::sqrt(dx * dx + dy * dy);
                s += W[i][j] * (dij - D[i][j]) * (dij - D[i][j]);
            }
        return s;
    }();
    constexpr double eps = 1e-6;

    for (int iter = 0; iter < maxIter; iter++) {
        // Build the Guttman transform: X^(k+1) = V^{-1} B(X^k) X^k
        // For diagonal V, this simplifies to coordinate-wise weighted average.
        std::vector<Point2D> newCoords(n);
        for (int i = 0; i < n; i++) {
            double sx = 0.0, sy = 0.0;
            for (int j = 0; j < n; j++) {
                if (i == j || W[i][j] < 1e-15) continue;
                double dx = lCoords[i].x - lCoords[j].x;
                double dy = lCoords[i].y - lCoords[j].y;
                double dij = std::sqrt(dx * dx + dy * dy);
                if (dij < 1e-10) dij = 1e-10;
                double ratio = D[i][j] / dij;
                sx += W[i][j] * (lCoords[j].x + ratio * (lCoords[i].x - lCoords[j].x));
                sy += W[i][j] * (lCoords[j].y + ratio * (lCoords[i].y - lCoords[j].y));
            }
            if (Wsum[i] > 1e-15) {
                newCoords[i].x = sx / Wsum[i];
                newCoords[i].y = sy / Wsum[i];
            } else {
                newCoords[i] = lCoords[i];
            }
        }
        lCoords = newCoords;

        double stress = [&]() -> double {
            double s = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++) {
                    double dx = lCoords[i].x - lCoords[j].x;
                    double dy = lCoords[i].y - lCoords[j].y;
                    double dij = std::sqrt(dx * dx + dy * dy);
                    s += W[i][j] * (dij - D[i][j]) * (dij - D[i][j]);
                }
            return s;
        }();
        if (prevStress > 1e-10 && (prevStress - stress) / prevStress < eps)
            break;
        prevStress = stress;
    }

    if (prevStress < bestStress) {
        bestStress = prevStress;
        bestCoords = lCoords;
    }
    } // end nInit loop

    coords = bestCoords;
    return bestStress;
}


// ==========================================================================
// Template library — pre-computed 2D coordinates for common scaffolds
// ==========================================================================

namespace detail_layout {

/// Scaffold template descriptor: ring sizes, atom count, coordinates.
struct ScaffoldTemplate {
    const char* name;
    int atomCount;
    /// Ring sizes list (terminated by 0).
    int ringSizes[8];
    /// Pre-computed normalised coordinates (unit bond length).
    /// Maximum 30 atoms per template.
    Point2D coords[30];
};

/// Canonical ring-size signature from SSSR for matching.
inline std::vector<int> ringSizeSignature(const MolGraph& g) {
    auto sssr = computeSSSR(g);
    std::vector<int> sizes;
    sizes.reserve(sssr.size());
    for (auto& r : sssr) sizes.push_back(static_cast<int>(r.size()));
    std::sort(sizes.begin(), sizes.end());
    return sizes;
}

// --------------------------------------------------------------------------
// The 10 most common scaffolds with pre-computed regular-polygon coordinates.
// Coordinates use a canonical atom ordering; template matching checks topology.
// --------------------------------------------------------------------------

/// Regular polygon coordinates centred at origin with unit edge length.
inline void regularPolygon(int n, Point2D* out) {
    double r = 0.5 / std::sin(M_PI / n); // circumradius for unit edge
    for (int i = 0; i < n; i++) {
        double angle = 2.0 * M_PI * i / n - M_PI / 2.0;
        out[i].x = r * std::cos(angle);
        out[i].y = r * std::sin(angle);
    }
}

// Pre-compute the 10 scaffold templates.
// Each entry stores normalised coordinates for the scaffold atoms.

inline std::vector<ScaffoldTemplate> buildTemplates() {
    std::vector<ScaffoldTemplate> templates;

    // Helper: fused-ring system builder.
    // Place ring A as a regular polygon, then ring B sharing an edge.
    auto fusedRings = [](int /*sizeA*/, int sizeB, int sharedEdgeA0, int sharedEdgeA1,
                         Point2D* coords, int& atomIdx) {
        // Ring A is already placed. Ring B shares edge (sharedEdgeA0, sharedEdgeA1).
        // Place ring B vertices on the opposite side of that edge.
        Point2D pA = coords[sharedEdgeA0];
        Point2D pB = coords[sharedEdgeA1];
        double mx = (pA.x + pB.x) * 0.5;
        double my = (pA.y + pB.y) * 0.5;
        double edx = pB.x - pA.x;
        double edy = pB.y - pA.y;
        // Normal to the shared edge (pointing away from ring A centre)
        double nx = -edy, ny = edx;
        // Find centre of ring A to determine outward direction
        double cax = 0, cay = 0;
        for (int i = 0; i < atomIdx; i++) { cax += coords[i].x; cay += coords[i].y; }
        cax /= atomIdx; cay /= atomIdx;
        // If normal points towards ring A centre, flip it
        double dot = nx * (cax - mx) + ny * (cay - my);
        if (dot > 0) { nx = -nx; ny = -ny; }
        // Place ring B as regular polygon, then transform
        double rB = 0.5 / std::sin(M_PI / sizeB);
        double edgeLen = std::sqrt(edx * edx + edy * edy);
        // Centre of ring B
        double centreDistFromMidpoint = rB * std::cos(M_PI / sizeB);
        double cx = mx + centreDistFromMidpoint * nx / std::sqrt(nx * nx + ny * ny);
        double cy = my + centreDistFromMidpoint * ny / std::sqrt(nx * nx + ny * ny);
        (void)edgeLen;
        const double startAngle = std::atan2(pB.y - cy, pB.x - cx);
        const double direction = ((pA.x - cx) * (pB.y - cy)
                                - (pA.y - cy) * (pB.x - cx)) > 0 ? 1.0 : -1.0;
        for (int k = 0; k < sizeB - 2; k++) {
            double angle = startAngle + direction * 2.0 * M_PI * (k + 1) / sizeB;
            coords[atomIdx].x = cx + rB * std::cos(angle);
            coords[atomIdx].y = cy + rB * std::sin(angle);
            atomIdx++;
        }
    };

    // 1. Benzene (6-ring)
    {
        ScaffoldTemplate t;
        t.name = "Benzene";
        t.atomCount = 6;
        t.ringSizes[0] = 6; t.ringSizes[1] = 0;
        regularPolygon(6, t.coords);
        templates.push_back(t);
    }

    // 2. Naphthalene (fused 6-6, 10 atoms)
    {
        ScaffoldTemplate t;
        t.name = "Naphthalene";
        t.atomCount = 10;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 3. Indole (fused 5-6, 9 atoms)
    {
        ScaffoldTemplate t;
        t.name = "Indole";
        t.atomCount = 9;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 5, 2, 3, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 4. Purine (fused 5-6 with N, 9 atoms)
    {
        ScaffoldTemplate t;
        t.name = "Purine";
        t.atomCount = 9;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 5, 3, 4, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 5. Steroid core (4 fused rings: 6-6-6-5, ABCD pattern, 17 atoms)
    {
        ScaffoldTemplate t;
        t.name = "Steroid";
        t.atomCount = 17;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 6;
        t.ringSizes[3] = 6; t.ringSizes[4] = 0;
        // Ring A (6-membered)
        regularPolygon(6, t.coords);
        int idx = 6;
        // Ring B fused to A
        fusedRings(6, 6, 1, 2, t.coords, idx);
        // Ring C fused to B (share edge between atoms 6 and 7)
        fusedRings(6, 6, 6, 7, t.coords, idx);
        // Ring D (5-membered) fused to C
        fusedRings(6, 5, 10, 11, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 6. Morphinan (bridged 3-ring system)
    {
        ScaffoldTemplate t;
        t.name = "Morphinan";
        t.atomCount = 14;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6;
        t.ringSizes[3] = 0;
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 3, 4, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 7. Quinoline (fused 6-6 with N, same geometry as naphthalene)
    {
        ScaffoldTemplate t;
        t.name = "Quinoline";
        t.atomCount = 10;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 2, 3, t.coords, idx);
        t.atomCount = idx;
        templates.push_back(t);
    }

    // 8. Biphenyl (two 6-rings connected by single bond, 12 atoms)
    {
        ScaffoldTemplate t;
        t.name = "Biphenyl";
        t.atomCount = 12;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        // First ring
        regularPolygon(6, t.coords);
        regularPolygon(6, t.coords + 6);
        for (int i = 6; i < 12; ++i) t.coords[i].y += 3.0;
        templates.push_back(t);
    }

    // 9. Cyclohexane (6-ring, flat projection)
    {
        ScaffoldTemplate t;
        t.name = "Cyclohexane";
        t.atomCount = 6;
        t.ringSizes[0] = 6; t.ringSizes[1] = 0;
        regularPolygon(6, t.coords);
        templates.push_back(t);
    }

    // 10. Piperidine (6-ring with N, same geometry)
    {
        ScaffoldTemplate t;
        t.name = "Piperidine";
        t.atomCount = 6;
        t.ringSizes[0] = 6; t.ringSizes[1] = 0;
        regularPolygon(6, t.coords);
        templates.push_back(t);
    }

    return templates;
}

inline const std::vector<ScaffoldTemplate>& getTemplates() {
    static const std::vector<ScaffoldTemplate> templates = buildTemplates();
    return templates;
}

} // namespace detail_layout


// ==========================================================================
// matchTemplate — scaffold template matching for common ring systems
// ==========================================================================

/**
 * Check if a molecule matches a known scaffold template.
 * Returns pre-computed coordinates scaled to targetBondLength if a match
 * is found, empty vector otherwise.
 *
 * Ring-size signatures select candidates; a topology mapping assigns each
 * template atom to the corresponding graph atom before coordinates are copied.
 *
 * @param g                The molecular graph.
 * @param targetBondLength Desired bond length for scaling (default 1.5).
 * @return                 Pre-computed coordinates, or empty if no match.
 */
inline std::vector<Point2D> matchTemplate(const MolGraph& g,
                                           double targetBondLength = 1.5) {
    const auto signature = detail_layout::ringSizeSignature(g);
    for (const auto& scaffold : detail_layout::getTemplates()) {
        if (scaffold.atomCount != g.n) continue;
        std::vector<int> scaffoldSignature;
        for (int ringIndex = 0; ringIndex < 8 && scaffold.ringSizes[ringIndex]; ++ringIndex)
            scaffoldSignature.push_back(scaffold.ringSizes[ringIndex]);
        std::sort(scaffoldSignature.begin(), scaffoldSignature.end());
        if (scaffoldSignature != signature) continue;
        std::vector<std::vector<bool>> edges(g.n, std::vector<bool>(g.n, false));
        std::vector<int> degrees(g.n, 0);
        for (int first = 0; first < g.n; ++first)
            for (int second = first + 1; second < g.n; ++second) {
                const double distance = std::hypot(scaffold.coords[first].x - scaffold.coords[second].x,
                                                   scaffold.coords[first].y - scaffold.coords[second].y);
                if (std::abs(distance - 1.0) < 1e-6) {
                    edges[first][second] = edges[second][first] = true;
                    ++degrees[first]; ++degrees[second];
                }
            }
        std::vector<int> atomOrder(g.n), mapping(g.n, -1);
        std::iota(atomOrder.begin(), atomOrder.end(), 0);
        std::stable_sort(atomOrder.begin(), atomOrder.end(), [&](int first, int second) {
            return degrees[first] > degrees[second];
        });
        std::vector<bool> used(g.n, false);
        auto match = [&](auto&& self, int depth) -> bool {
            if (depth == g.n) return true;
            const int templateAtom = atomOrder[depth];
            for (int graphAtom = 0; graphAtom < g.n; ++graphAtom) {
                if (used[graphAtom] || degrees[templateAtom] != static_cast<int>(g.neighbors[graphAtom].size()))
                    continue;
                bool compatible = true;
                for (int previous = 0; previous < depth; ++previous) {
                    const int other = atomOrder[previous];
                    if (edges[templateAtom][other] != (g.bondOrder(graphAtom, mapping[other]) != 0)) {
                        compatible = false; break;
                    }
                }
                if (!compatible) continue;
                mapping[templateAtom] = graphAtom;
                used[graphAtom] = true;
                if (self(self, depth + 1)) return true;
                used[graphAtom] = false;
            }
            return false;
        };
        if (!match(match, 0)) continue;
        std::vector<Point2D> result(g.n);
        for (int templateAtom = 0; templateAtom < g.n; ++templateAtom)
            result[mapping[templateAtom]] = {scaffold.coords[templateAtom].x * targetBondLength,
                                             scaffold.coords[templateAtom].y * targetBondLength};
        return result;
    }
    return {};
}


// ==========================================================================
// Point3D — 3D coordinate type for optional 3D layout
// ==========================================================================

struct Point3D {
    double x = 0.0, y = 0.0, z = 0.0;
};

// ==========================================================================
// Coordinate transforms — 2D and 3D (translate, rotate, scale, mirror, align)
// ==========================================================================

namespace transform {

/// Translate all 2D coordinates by (dx, dy).
inline void translate2D(std::vector<Point2D>& coords, double dx, double dy) {
    for (auto& p : coords) { p.x += dx; p.y += dy; }
}

/// Translate all 3D coordinates by (dx, dy, dz).
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void translate3D(std::vector<Point3D>& coords, double dx, double dy, double dz) {
    for (auto& p : coords) { p.x += dx; p.y += dy; p.z += dz; }
}

/// Rotate 2D coordinates by angle (radians) about the centroid.
inline void rotate2D(std::vector<Point2D>& coords, double angle) {
    if (coords.empty()) return;
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());
    double cosA = std::cos(angle), sinA = std::sin(angle);
    for (auto& p : coords) {
        double rx = p.x - cx, ry = p.y - cy;
        p.x = cx + rx * cosA - ry * sinA;
        p.y = cy + rx * sinA + ry * cosA;
    }
}

/// Rotate 2D coordinates about a specific centre point.
inline void rotate2DAbout(std::vector<Point2D>& coords, double angle,
                          double cx, double cy) {
    double cosA = std::cos(angle), sinA = std::sin(angle);
    for (auto& p : coords) {
        double rx = p.x - cx, ry = p.y - cy;
        p.x = cx + rx * cosA - ry * sinA;
        p.y = cy + rx * sinA + ry * cosA;
    }
}

/// Rotate 3D coordinates about Z-axis by angle (radians) around centroid.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void rotate3DZ(std::vector<Point3D>& coords, double angle) {
    if (coords.empty()) return;
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());
    double cosA = std::cos(angle), sinA = std::sin(angle);
    for (auto& p : coords) {
        double rx = p.x - cx, ry = p.y - cy;
        p.x = cx + rx * cosA - ry * sinA;
        p.y = cy + rx * sinA + ry * cosA;
    }
}

/// Rotate 3D coordinates about X-axis by angle (radians) around centroid.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void rotate3DX(std::vector<Point3D>& coords, double angle) {
    if (coords.empty()) return;
    double cy = 0, cz = 0;
    for (auto& p : coords) { cy += p.y; cz += p.z; }
    cy /= static_cast<double>(coords.size());
    cz /= static_cast<double>(coords.size());
    double cosA = std::cos(angle), sinA = std::sin(angle);
    for (auto& p : coords) {
        double ry = p.y - cy, rz = p.z - cz;
        p.y = cy + ry * cosA - rz * sinA;
        p.z = cz + ry * sinA + rz * cosA;
    }
}

/// Rotate 3D coordinates about Y-axis by angle (radians) around centroid.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void rotate3DY(std::vector<Point3D>& coords, double angle) {
    if (coords.empty()) return;
    double cx = 0, cz = 0;
    for (auto& p : coords) { cx += p.x; cz += p.z; }
    cx /= static_cast<double>(coords.size());
    cz /= static_cast<double>(coords.size());
    double cosA = std::cos(angle), sinA = std::sin(angle);
    for (auto& p : coords) {
        double rx = p.x - cx, rz = p.z - cz;
        p.x = cx + rx * cosA + rz * sinA;
        p.z = cz - rx * sinA + rz * cosA;
    }
}

/// Uniform scale about centroid.
inline void scale2D(std::vector<Point2D>& coords, double factor) {
    if (coords.empty()) return;
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());
    for (auto& p : coords) {
        p.x = cx + (p.x - cx) * factor;
        p.y = cy + (p.y - cy) * factor;
    }
}

/// Uniform 3D scale about centroid.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void scale3D(std::vector<Point3D>& coords, double factor) {
    if (coords.empty()) return;
    double cx = 0, cy = 0, cz = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; cz += p.z; }
    double n = static_cast<double>(coords.size());
    cx /= n; cy /= n; cz /= n;
    for (auto& p : coords) {
        p.x = cx + (p.x - cx) * factor;
        p.y = cy + (p.y - cy) * factor;
        p.z = cz + (p.z - cz) * factor;
    }
}

/// Mirror 2D coordinates about X-axis (flip vertically).
inline void mirrorX(std::vector<Point2D>& coords) {
    double cy = 0;
    for (auto& p : coords) cy += p.y;
    cy /= static_cast<double>(coords.size());
    for (auto& p : coords) p.y = 2.0 * cy - p.y;
}

/// Mirror 2D coordinates about Y-axis (flip horizontally).
inline void mirrorY(std::vector<Point2D>& coords) {
    double cx = 0;
    for (auto& p : coords) cx += p.x;
    cx /= static_cast<double>(coords.size());
    for (auto& p : coords) p.x = 2.0 * cx - p.x;
}

/// Centre coordinates at origin.
inline void center2D(std::vector<Point2D>& coords) {
    if (coords.empty()) return;
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());
    for (auto& p : coords) { p.x -= cx; p.y -= cy; }
}

/// Centre 3D coordinates at origin.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline void center3D(std::vector<Point3D>& coords) {
    if (coords.empty()) return;
    double cx = 0, cy = 0, cz = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; cz += p.z; }
    double n = static_cast<double>(coords.size());
    cx /= n; cy /= n; cz /= n;
    for (auto& p : coords) { p.x -= cx; p.y -= cy; p.z -= cz; }
}

/// Normalise bond lengths to target.
inline void normaliseBondLength(const MolGraph& g, std::vector<Point2D>& coords,
                                double target = 1.5) {
    if (coords.size() < static_cast<size_t>(g.n))
        throw std::invalid_argument("Coordinates must cover every atom");
    if (g.n < 2) return;
    double sumLen = 0.0;
    int count = 0;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                sumLen += std::sqrt(dx * dx + dy * dy);
                count++;
            }
        }
    }
    if (count > 0 && sumLen > 1e-10) {
        double avg = sumLen / count;
        double factor = target / avg;
        scale2D(coords, factor);
    }
}

/// Canonical orientation: rotate to align the longest axis with the horizontal,
/// then by at most 15 degrees more, so that the bonds lie as near as they will
/// go to multiples of 30 degrees.
inline void canonicalOrientation(const MolGraph& g, std::vector<Point2D>& coords) {
    if (coords.size() < 2) return;

    // Find principal axis via PCA (covariance matrix eigenvectors)
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());

    double sxx = 0, sxy = 0, syy = 0;
    for (auto& p : coords) {
        double dx = p.x - cx, dy = p.y - cy;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }

    // Principal axis angle from the 2x2 covariance matrix.  A symmetric shape
    // (benzene, a large ring) has no principal axis, and terms at rounding
    // level are taken as zero so that builds which round differently agree.
    const double noise = 1e-9 * (sxx + syy);
    const double twiceCovariance = std::abs(2.0 * sxy) > noise ? 2.0 * sxy : 0.0;
    const double spread = std::abs(sxx - syy) > noise ? sxx - syy : 0.0;
    double angle = twiceCovariance == 0.0 && spread == 0.0 ? 0.0 : 0.5 * std::atan2(twiceCovariance, spread);

    // Turn by at most 15 degrees more, so the bonds lie as near as they can to
    // the 30-degree grid: the mean of their directions taken twelvefold.
    const size_t atomCount = std::min(coords.size(), static_cast<size_t>(g.n));
    double gridX = 0, gridY = 0;
    int bonds = 0;
    for (size_t atom = 0; atom < atomCount; ++atom)
        for (int neighbour : g.neighbors[atom]) {
            if (neighbour <= static_cast<int>(atom) || static_cast<size_t>(neighbour) >= atomCount) continue;
            const double direction = std::atan2(coords[neighbour].y - coords[atom].y,
                                                coords[neighbour].x - coords[atom].x) - angle;
            gridX += std::cos(12.0 * direction);
            gridY += std::sin(12.0 * direction);
            ++bonds;
        }
    const double resultant = std::hypot(gridX, gridY);
    if (resultant > 1e-6 * bonds) {
        if (std::abs(gridY) <= 1e-9 * resultant) gridY = 0.0;
        angle += std::atan2(gridY, gridX) / 12.0;
    }

    rotate2DAbout(coords, -angle, cx, cy);
}

/// Align 2D coordinates to a reference via Kabsch-like rotation minimising RMSD.
/// Both vectors must have the same length.
inline double align2D(std::vector<Point2D>& coords,
                      const std::vector<Point2D>& reference) {
    size_t n = std::min(coords.size(), reference.size());
    if (n < 2) return 0.0;

    // Centroid
    double cx1 = 0, cy1 = 0, cx2 = 0, cy2 = 0;
    for (size_t i = 0; i < n; i++) {
        cx1 += coords[i].x; cy1 += coords[i].y;
        cx2 += reference[i].x; cy2 += reference[i].y;
    }
    cx1 /= n; cy1 /= n; cx2 /= n; cy2 /= n;

    // Cross-covariance for 2D rotation
    double sumCross = 0, sumDot = 0;
    for (size_t i = 0; i < n; i++) {
        double ax = coords[i].x - cx1, ay = coords[i].y - cy1;
        double bx = reference[i].x - cx2, by = reference[i].y - cy2;
        sumDot   += ax * bx + ay * by;
        sumCross += ax * by - ay * bx;
    }
    double theta = std::atan2(sumCross, sumDot);

    // Apply rotation and translation
    double cosT = std::cos(theta), sinT = std::sin(theta);
    double rmsd = 0.0;
    for (size_t i = 0; i < n; i++) {
        double rx = coords[i].x - cx1, ry = coords[i].y - cy1;
        coords[i].x = cx2 + rx * cosT - ry * sinT;
        coords[i].y = cy2 + rx * sinT + ry * cosT;
        double dx = coords[i].x - reference[i].x;
        double dy = coords[i].y - reference[i].y;
        rmsd += dx * dx + dy * dy;
    }
    return std::sqrt(rmsd / n);
}

/// Project 3D coordinates to 2D (XY plane).
inline std::vector<Point2D> projectTo2D(const std::vector<Point3D>& coords3d) {
    std::vector<Point2D> result;
    result.reserve(coords3d.size());
    for (auto& p : coords3d) result.push_back({p.x, p.y});
    return result;
}

/// Lift 2D coordinates to 3D (Z = 0).
inline std::vector<Point3D> liftTo3D(const std::vector<Point2D>& coords2d) {
    std::vector<Point3D> result;
    result.reserve(coords2d.size());
    for (auto& p : coords2d) result.push_back({p.x, p.y, 0.0});
    return result;
}

/// Compute bounding box of 2D coordinates.
/// Returns (minX, minY, maxX, maxY).
inline std::array<double, 4> boundingBox2D(const std::vector<Point2D>& coords) {
    if (coords.empty()) return {0, 0, 0, 0};
    double minX = coords[0].x, minY = coords[0].y;
    double maxX = coords[0].x, maxY = coords[0].y;
    for (auto& p : coords) {
        minX = std::min(minX, p.x); minY = std::min(minY, p.y);
        maxX = std::max(maxX, p.x); maxY = std::max(maxY, p.y);
    }
    return {minX, minY, maxX, maxY};
}

/// Compute RMSD between two coordinate sets.
[[deprecated("Unused; the 3-D transforms and rmsd2D are withdrawn.")]]
inline double rmsd2D(const std::vector<Point2D>& a, const std::vector<Point2D>& b) {
    size_t n = std::min(a.size(), b.size());
    if (n == 0) return 0.0;
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
        double dx = a[i].x - b[i].x, dy = a[i].y - b[i].y;
        sum += dx * dx + dy * dy;
    }
    return std::sqrt(sum / n);
}

} // namespace transform


// ==========================================================================
// Overlap resolution — resolve atom-atom overlaps in 2D layout
// ==========================================================================

namespace detail_layout {

/// Count atom-atom overlaps (distance < threshold).
inline int countOverlaps(const std::vector<Point2D>& coords, double threshold = 0.3) {
    int n = static_cast<int>(coords.size());
    int count = 0;
    double thr2 = threshold * threshold;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = coords[i].x - coords[j].x;
            double dy = coords[i].y - coords[j].y;
            if (dx * dx + dy * dy < thr2) count++;
        }
    }
    return count;
}

} // namespace detail_layout

/// Resolve atom-atom overlaps by pushing overlapping atoms apart.
/// Uses iterative repulsion with displacement limiting.
/// @param coords   In/out 2D coordinates.
/// @param threshold Minimum allowed inter-atom distance (default 0.3).
/// @param maxIter  Maximum iterations (default 100).
/// @return Number of remaining overlaps (0 = fully resolved).
inline int resolveOverlaps(std::vector<Point2D>& coords,
                           double threshold = 0.3, int maxIter = 100) {
    int n = static_cast<int>(coords.size());
    if (n < 2) return 0;
    double thr2 = threshold * threshold;

    for (int iter = 0; iter < maxIter; iter++) {
        bool anyOverlap = false;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                double d2 = dx * dx + dy * dy;
                if (d2 >= thr2) continue;
                anyOverlap = true;
                double d = std::sqrt(d2);
                if (d < 1e-10) { dx = 0.1; dy = 0.1; d = std::sqrt(0.02); }
                double push = (threshold - d) * 0.55; // slightly more than half
                double ux = dx / d, uy = dy / d;
                coords[i].x -= push * ux;
                coords[i].y -= push * uy;
                coords[j].x += push * ux;
                coords[j].y += push * uy;
            }
        }
        if (!anyOverlap) return 0;
    }
    return detail_layout::countOverlaps(coords, threshold);
}


// ==========================================================================
// Expanded scaffold template library — 40+ pharmaceutical scaffolds
// ==========================================================================

namespace detail_layout {

/// Extended template library with comprehensive pharmaceutical scaffolds.
/// Coordinates computed analytically using regular polygon fusion.
inline std::vector<ScaffoldTemplate> buildExtendedTemplates() {
    std::vector<ScaffoldTemplate> templates;

    // Reuse fusedRings helper
    auto fusedRings = [](int /*sizeA*/, int sizeB, int sharedEdgeA0, int sharedEdgeA1,
                         Point2D* coords, int& atomIdx) {
        Point2D pA = coords[sharedEdgeA0];
        Point2D pB = coords[sharedEdgeA1];
        double mx = (pA.x + pB.x) * 0.5;
        double my = (pA.y + pB.y) * 0.5;
        double edx = pB.x - pA.x, edy = pB.y - pA.y;
        double nx = -edy, ny = edx;
        double cax = 0, cay = 0;
        for (int i = 0; i < atomIdx; i++) { cax += coords[i].x; cay += coords[i].y; }
        cax /= atomIdx; cay /= atomIdx;
        double dot = nx * (cax - mx) + ny * (cay - my);
        if (dot > 0) { nx = -nx; ny = -ny; }
        double rB = 0.5 / std::sin(M_PI / sizeB);
        double nmag = std::sqrt(nx * nx + ny * ny);
        double centreDistFromMidpoint = rB * std::cos(M_PI / sizeB);
        double cx = mx + centreDistFromMidpoint * nx / nmag;
        double cy = my + centreDistFromMidpoint * ny / nmag;
        const double startAngle = std::atan2(pB.y - cy, pB.x - cx);
        const double direction = ((pA.x - cx) * (pB.y - cy)
                                - (pA.y - cy) * (pB.x - cx)) > 0 ? 1.0 : -1.0;
        for (int k = 0; k < sizeB - 2; k++) {
            double angle = startAngle + direction * 2.0 * M_PI * (k + 1) / sizeB;
            coords[atomIdx].x = cx + rB * std::cos(angle);
            coords[atomIdx].y = cy + rB * std::sin(angle);
            atomIdx++;
        }
    };

    // Helper: create simple ring template
    auto addRing = [&](const char* name, int size) {
        ScaffoldTemplate t;
        t.name = name;
        t.atomCount = size;
        t.ringSizes[0] = size; t.ringSizes[1] = 0;
        regularPolygon(size, t.coords);
        templates.push_back(t);
    };

    // Helper: create fused 2-ring template
    auto addFused2 = [&](const char* name, int sA, int sB, int e0, int e1) {
        ScaffoldTemplate t;
        t.name = name;
        regularPolygon(sA, t.coords);
        int idx = sA;
        fusedRings(sA, sB, e0, e1, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = sA; t.ringSizes[1] = sB; t.ringSizes[2] = 0;
        templates.push_back(t);
    };

    // --- Single rings (3-8 membered) ---
    addRing("cyclopropane", 3);      // cyclopropane / aziridine
    addRing("cyclobutane", 4);       // cyclobutane / azetidine / oxetane
    addRing("cyclopentane", 5);      // cyclopentane / furan / thiophene / pyrrole
    addRing("cyclohexane", 6);       // cyclohexane / benzene / pyridine / piperidine
    addRing("cycloheptane", 7);      // tropane / azepane
    addRing("cyclooctane", 8);       // cyclooctane

    // --- Fused 2-ring systems ---
    addFused2("naphthalene", 6, 6, 1, 2);         // 6-6 fused
    addFused2("indole", 6, 5, 2, 3);              // 6-5 fused (indole, benzofuran, benzothiophene)
    addFused2("indane", 6, 5, 1, 2);              // 6-5 fused (indane, dihydrobenzofuran)
    addFused2("quinoline", 6, 6, 2, 3);           // 6-6 fused (quinoline, isoquinoline)
    addFused2("benzimidazole", 6, 5, 3, 4);       // 6-5 fused (benzimidazole, purine)
    addFused2("isobenzofuran", 6, 5, 0, 1);       // 6-5 fused (isobenzofuran, phthalimide)
    addFused2("pentalene", 5, 5, 1, 2);           // 5-5 fused
    addFused2("azulene", 5, 7, 1, 2);             // 5-7 fused (non-alternant aromatic)
    addFused2("tetrahydroisoquinoline", 6, 6, 3, 4); // 6-6 THiQ scaffold

    // --- Fused 3-ring systems ---
    // Anthracene: linear 6-6-6
    {
        ScaffoldTemplate t;
        t.name = "anthracene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 6, 7, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // Phenanthrene: angular 6-6-6
    {
        ScaffoldTemplate t;
        t.name = "phenanthrene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 7, 8, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // Acridine: angular 6-6-6 (different fusion points)
    {
        ScaffoldTemplate t;
        t.name = "acridine";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 2, 3, t.coords, idx);
        fusedRings(6, 6, 4, 5, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // Carbazole: 6-5-6
    {
        ScaffoldTemplate t;
        t.name = "carbazole";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 5, 2, 3, t.coords, idx);
        fusedRings(5, 6, 6, 7, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // Fluorene: 6-5-6 (cyclopenta[d,e,f]phenanthrene variant)
    {
        ScaffoldTemplate t;
        t.name = "fluorene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 5, 1, 2, t.coords, idx);
        fusedRings(5, 6, 6, 7, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // Xanthene: 6-6-6 (O-bridged)
    {
        ScaffoldTemplate t;
        t.name = "xanthene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 2, 3, t.coords, idx);
        fusedRings(6, 6, 8, 9, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        templates.push_back(t);
    }

    // --- Fused 4-ring systems ---

    // Steroid core: 6-6-6-5 (ABCD)
    {
        ScaffoldTemplate t;
        t.name = "steroid";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 6, 7, t.coords, idx);
        fusedRings(6, 5, 10, 11, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 6;
        t.ringSizes[3] = 6; t.ringSizes[4] = 0;
        templates.push_back(t);
    }

    // Pyrene: 4 fused 6-rings (peri-fused PAH)
    {
        ScaffoldTemplate t;
        t.name = "pyrene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 3, 4, t.coords, idx);
        fusedRings(6, 6, 6, 7, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6;
        t.ringSizes[3] = 6; t.ringSizes[4] = 0;
        templates.push_back(t);
    }

    // Naphthacene (tetracene): linear 6-6-6-6
    {
        ScaffoldTemplate t;
        t.name = "naphthacene";
        regularPolygon(6, t.coords);
        int idx = 6;
        fusedRings(6, 6, 1, 2, t.coords, idx);
        fusedRings(6, 6, 6, 7, t.coords, idx);
        fusedRings(6, 6, 10, 11, t.coords, idx);
        t.atomCount = idx;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6;
        t.ringSizes[3] = 6; t.ringSizes[4] = 0;
        templates.push_back(t);
    }

    // --- Spiro compounds ---
    // Spiro[4.5]decane: 5-ring spiro-fused to 6-ring
    {
        ScaffoldTemplate t;
        t.name = "spiro45";
        regularPolygon(5, t.coords);
        int idx = 5;
        // Place 6-ring sharing only the spiro atom (atom 0), on opposite side
        double r6 = 0.5 / std::sin(M_PI / 6.0);
        // Spiro atom is coords[0], place 6-ring centred opposite to 5-ring
        double cx5 = 0, cy5 = 0;
        for (int i = 0; i < 5; i++) { cx5 += t.coords[i].x; cy5 += t.coords[i].y; }
        cx5 /= 5; cy5 /= 5;
        double dx = t.coords[0].x - cx5, dy = t.coords[0].y - cy5;
        double mag = std::sqrt(dx * dx + dy * dy);
        if (mag < 1e-10) { dx = 1; dy = 0; mag = 1; }
        double cx6 = t.coords[0].x + r6 * dx / mag;
        double cy6 = t.coords[0].y + r6 * dy / mag;
        double baseAngle = std::atan2(dy, dx);
        for (int k = 1; k < 6; k++) {
            double angle = baseAngle + M_PI + (2.0 * M_PI * k) / 6.0;
            t.coords[idx].x = cx6 + r6 * std::cos(angle);
            t.coords[idx].y = cy6 + r6 * std::sin(angle);
            idx++;
        }
        t.atomCount = idx;
        t.ringSizes[0] = 5; t.ringSizes[1] = 6; t.ringSizes[2] = 0;
        templates.push_back(t);
    }

    // --- Bridged systems ---
    // Norbornane (bicyclo[2.2.1]heptane): 7 atoms
    {
        ScaffoldTemplate t;
        t.name = "norbornane";
        t.atomCount = 7;
        t.ringSizes[0] = 5; t.ringSizes[1] = 5; t.ringSizes[2] = 0;
        // Manual coordinates for norbornane boat-like 2D projection
        double bl = 1.0; // unit bond length
        t.coords[0] = { 0.0,       0.0};
        t.coords[1] = { bl,        0.0};
        t.coords[2] = { 1.5 * bl, -0.87 * bl};
        t.coords[3] = { bl,       -1.73 * bl};
        t.coords[4] = { 0.0,      -1.73 * bl};
        t.coords[5] = {-0.5 * bl, -0.87 * bl};
        t.coords[6] = { 0.5 * bl, -0.60 * bl}; // bridgehead
        templates.push_back(t);
    }

    // Adamantane (tricyclo[3.3.1.1]decane): 10 atoms
    {
        ScaffoldTemplate t;
        t.name = "adamantane";
        t.atomCount = 10;
        t.ringSizes[0] = 6; t.ringSizes[1] = 6; t.ringSizes[2] = 6; t.ringSizes[3] = 0;
        // Diamond lattice projection
        double s = 0.8;
        t.coords[0] = { 0.0,       1.5 * s};
        t.coords[1] = {-1.0 * s,   0.5 * s};
        t.coords[2] = { 1.0 * s,   0.5 * s};
        t.coords[3] = { 0.0,      -0.5 * s};
        t.coords[4] = {-1.5 * s,  -0.5 * s};
        t.coords[5] = { 1.5 * s,  -0.5 * s};
        t.coords[6] = {-0.5 * s,  -1.5 * s};
        t.coords[7] = { 0.5 * s,  -1.5 * s};
        t.coords[8] = { 0.0,       0.0};
        t.coords[9] = { 0.0,      -2.0 * s};
        templates.push_back(t);
    }

    // --- Common drug scaffolds ---
    // Benzodiazepine: 6-7 fused
    addFused2("benzodiazepine", 6, 7, 2, 3);

    // Dihydropyridine: just a 6-ring
    addRing("dihydropyridine", 6);

    // Morpholine: 6-ring
    addRing("morpholine", 6);

    // Piperazine: 6-ring
    addRing("piperazine", 6);

    // Imidazole: 5-ring
    addRing("imidazole", 5);

    // Oxazole: 5-ring
    addRing("oxazole", 5);

    // Thiazole: 5-ring
    addRing("thiazole", 5);

    // Pyrimidine: 6-ring
    addRing("pyrimidine", 6);

    // Triazine: 6-ring
    addRing("triazine", 6);

    // Tetrazole: 5-ring
    addRing("tetrazole", 5);

    return templates;
}

inline const std::vector<ScaffoldTemplate>& getExtendedTemplates() {
    static const std::vector<ScaffoldTemplate> templates = buildExtendedTemplates();
    return templates;
}

} // namespace detail_layout


// ==========================================================================
// Clash removal and refinement that keep ring templates exact
// ==========================================================================

namespace detail_layout {

/// Bonds within this distance of the unit length count as exact, which allows
/// for rounding to the grid of onGrid.
constexpr double kBondLengthTolerance = 1e-5;
/// Non-bonded atoms closer than this, in bond lengths, clash.
constexpr double kMinNonBondedDistance = 0.8;
/// Non-bonded atoms closer than this, in bond lengths, hide one another.
constexpr double kHiddenDistance = 0.4;
/// Refinement parts non-bonded atoms to at least this distance, clear of kHiddenDistance.
constexpr double kPartedDistance = 0.5;

/// True when `score` beats `best` by more than rounding.  Candidates that tie
/// keep the earlier one, so builds that round differently choose alike.
inline bool clearlyHigher(double score, double best) {
    return score - best > 1e-9 * std::max(1.0, std::abs(best));
}

/// Direction, from `atom`, of the bisector of a gap between its bonds to
/// `bonded` whose unit-length point lies furthest from the atoms in `others`.
/// A gap inside a ring (less than a half-turn between two bonds that
/// `ringAngles` pairs at this atom) is passed over while any other gap is
/// open, so at a ring fusion the new bond takes the exterior gap.  With
/// `terminal` set, as for a bond to an atom with no other neighbour, an
/// exterior gap counts as open only while its unit-length point keeps clear
/// of the atoms in `others`, so that a 9-alpha fluorine can take a ring's
/// interior gap when the angular methyl already fills the bay.
inline double clearestGapBisector(const std::vector<Point2D>& coords, int atom,
                                  const std::vector<int>& bonded, const std::vector<int>& others,
                                  const std::vector<std::pair<int, int>>& ringAngles = {},
                                  bool terminal = false) {
    if (bonded.empty()) return 0;
    auto direction = [&](int neighbour) {
        return std::atan2(coords[neighbour].y - coords[atom].y, coords[neighbour].x - coords[atom].x);
    };
    // Gaps are taken anticlockwise from the first bond, so their order never
    // turns on which side of the -x axis a bond falls.
    const double reference = direction(bonded.front());
    std::vector<std::pair<double, int>> bonds{{0.0, bonded.front()}};
    for (size_t k = 1; k < bonded.size(); ++k) {
        const double angle = std::remainder(direction(bonded[k]) - reference, 2 * M_PI);
        bonds.push_back({angle < 0 ? angle + 2 * M_PI : angle, bonded[k]});
    }
    std::sort(bonds.begin(), bonds.end());
    auto gapEnd = [&](size_t k) { return k + 1 < bonds.size() ? bonds[k + 1].first : 2 * M_PI; };
    auto insideRing = [&](size_t k) {
        const int first = bonds[k].second, second = bonds[(k + 1) % bonds.size()].second;
        if (gapEnd(k) - bonds[k].first >= M_PI) return false;
        return std::any_of(ringAngles.begin(), ringAngles.end(), [&](const std::pair<int, int>& pair) {
            return (pair.first == first && pair.second == second) || (pair.first == second && pair.second == first);
        });
    };
    auto pointIn = [&](size_t k) {
        const double bisector = reference + 0.5 * (bonds[k].first + gapEnd(k));
        return Point2D{coords[atom].x + std::cos(bisector), coords[atom].y + std::sin(bisector)};
    };
    auto keepsClear = [&](size_t k) {
        const Point2D point = pointIn(k);
        return std::none_of(others.begin(), others.end(), [&](int other) {
            return other != atom && std::hypot(point.x - coords[other].x, point.y - coords[other].y)
                                        < kMinNonBondedDistance;
        });
    };
    bool exteriorOpen = false;
    for (size_t k = 0; k < bonds.size(); ++k)
        exteriorOpen = exteriorOpen || (!insideRing(k) && (!terminal || keepsClear(k)));
    double best = 0, bestScore = 0;
    bool found = false;
    for (size_t k = 0; k < bonds.size(); ++k) {
        if (exteriorOpen && insideRing(k)) continue;
        const double bisector = reference + 0.5 * (bonds[k].first + gapEnd(k));
        const Point2D point = pointIn(k);
        double score = 0;
        for (int other : others) {
            if (other == atom) continue;
            const double dx = point.x - coords[other].x, dy = point.y - coords[other].y;
            score -= 1.0 / std::max(0.01, dx * dx + dy * dy);
        }
        if (!found || clearlyHigher(score, bestScore)) { bestScore = score; best = bisector; found = true; }
    }
    return best;
}

/// True when a bonded pair is one bond length apart and a non-bonded pair does not clash.
inline bool isExactPair(const MolGraph& g, const std::vector<Point2D>& coords, int first, int second) {
    const double distance = std::hypot(coords[first].x - coords[second].x,
                                       coords[first].y - coords[second].y);
    return g.bondOrder(first, second) ? std::abs(distance - 1.0) <= kBondLengthTolerance
                                      : distance >= kMinNonBondedDistance;
}

/// True when every bond has unit length, no bonds cross and no atoms clash.
inline bool isCleanLayout(const MolGraph& g, const std::vector<Point2D>& coords) {
    for (int first = 0; first < g.n; ++first)
        for (int second = first + 1; second < g.n; ++second)
            if (!isExactPair(g, coords, first, second)) return false;
    return countCrossings(g, coords) == 0;
}

/// Members of each rigid group, indexed by group id.
inline std::vector<std::vector<int>> groupMembers(const std::vector<int>& group) {
    std::vector<std::vector<int>> members;
    for (int atom = 0; atom < static_cast<int>(group.size()); ++atom) {
        if (group[atom] < 0) continue;
        if (group[atom] >= static_cast<int>(members.size())) members.resize(group[atom] + 1);
        members[group[atom]].push_back(atom);
    }
    return members;
}

/// Rigid group of every atom (-1 for none): ring systems whose own geometry is
/// exact.  A system spiro-fused to a group joins it while the two keep clear of
/// each other; otherwise it moves on its own, without the shared atom, so that
/// refinement can turn it clear while its bonds to that atom settle.  Where
/// three or more systems meet at one atom (a chelating metal, say) and clash,
/// regular rings cannot all fit round it, and those systems stay free.
inline std::vector<int> rigidRingGroups(const MolGraph& g, const std::vector<Point2D>& coords,
                                        const std::vector<RingSystem>& systems) {
    auto allExact = [&](const std::vector<int>& atoms) {
        for (size_t first = 0; first < atoms.size(); ++first)
            for (size_t second = first + 1; second < atoms.size(); ++second)
                if (!isExactPair(g, coords, atoms[first], atoms[second])) return false;
        return true;
    };
    std::vector<bool> loose(systems.size(), false);
    for (int centre = 0; centre < g.n; ++centre) {
        std::vector<size_t> meeting;
        std::set<int> atoms;
        for (size_t index = 0; index < systems.size(); ++index)
            if (systems[index].atoms.count(centre)) {
                meeting.push_back(index);
                atoms.insert(systems[index].atoms.begin(), systems[index].atoms.end());
            }
        if (meeting.size() >= 3 && !allExact(std::vector<int>(atoms.begin(), atoms.end())))
            for (size_t index : meeting) loose[index] = true;
    }
    std::vector<int> group(g.n, -1);
    int groups = 0;
    for (size_t index = 0; index < systems.size(); ++index) {
        const auto& system = systems[index];
        const std::vector<int> atoms(system.atoms.begin(), system.atoms.end());
        if (loose[index] || !allExact(atoms)) continue;
        std::set<int> joined;
        for (int atom : atoms)
            if (group[atom] >= 0) joined.insert(group[atom]);
        bool clear = true;
        for (int other = 0; other < g.n && clear; ++other)
            if (group[other] >= 0 && joined.count(group[other]) && !system.atoms.count(other))
                for (size_t first = 0; first < atoms.size() && clear; ++first)
                    clear = isExactPair(g, coords, atoms[first], other);
        const int id = groups++;
        for (int atom : atoms) {
            const int previous = group[atom];
            if (previous >= 0 && !clear) continue;
            if (previous >= 0) std::replace(group.begin(), group.end(), previous, id);
            group[atom] = id;
        }
    }
    // A terminal substituent of a rigid system (an angular methyl, a ring halogen)
    // moves with it, so that refinement cannot fold its bond onto a ring bond,
    // unless it already clashes with the system or another such substituent.
    const auto members = groupMembers(group);
    std::vector<bool> joinedGroup(g.n, false);
    for (int atom = 0; atom < g.n; ++atom) {
        if (group[atom] >= 0 || g.neighbors[atom].size() != 1) continue;
        const int root = g.neighbors[atom].front();
        if (group[root] < 0) continue;
        bool clear = true;
        for (int member : members[group[root]])
            clear = clear && isExactPair(g, coords, atom, member);
        for (int other = 0; other < atom && clear; ++other)
            if (joinedGroup[other] && group[other] == group[root]) clear = isExactPair(g, coords, atom, other);
        if (!clear) continue;
        group[atom] = group[root];
        joinedGroup[atom] = true;
    }
    return group;
}

/// Crossing bonds and clashing atom pairs, compared by crossings, then by the
/// pairs close enough to hide an atom, then by clashes.
struct LayoutDefects {
    int hidden = 0;
    int crossings = 0;
    int clashes = 0;
    std::vector<bool> involved; // atoms on a crossing bond or in a clashing pair

    bool none() const { return crossings == 0 && clashes == 0; }
    bool fewerThan(const LayoutDefects& other) const {
        if (crossings != other.crossings) return crossings < other.crossings;
        return hidden != other.hidden ? hidden < other.hidden : clashes < other.clashes;
    }
    bool noneAdded(const LayoutDefects& before) const {
        return crossings <= before.crossings && hidden <= before.hidden && clashes <= before.clashes;
    }
};

/**
 * The defects of a layout.  Atoms and bonds are binned in cells one bond
 * length wide, so that only nearby pairs are tested: a clashing pair lies in
 * the same or neighbouring cells, and two bonds can cross only where their
 * bounding boxes meet, which happens first in one cell that both reach.
 */
inline LayoutDefects findDefects(const MolGraph& g, const std::vector<Point2D>& coords) {
    LayoutDefects defects;
    defects.involved.assign(g.n, false);
    auto clash = [&](int first, int second) {
        if (g.bondOrder(first, second) || isExactPair(g, coords, first, second)) return;
        ++defects.clashes;
        defects.hidden += std::hypot(coords[first].x - coords[second].x,
                                     coords[first].y - coords[second].y) < kHiddenDistance;
        defects.involved[first] = defects.involved[second] = true;
    };
    auto cross = [&](const std::pair<int, int>& one, const std::pair<int, int>& other) {
        const auto [a1, a2] = one;
        const auto [b1, b2] = other;
        if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) return;
        if (!segmentsCross(coords[a1], coords[a2], coords[b1], coords[b2])) return;
        ++defects.crossings;
        defects.involved[a1] = defects.involved[a2] = defects.involved[b1] = defects.involved[b2] = true;
    };
    std::vector<std::pair<int, int>> bonds;
    for (int first = 0; first < g.n; ++first)
        for (int second : g.neighbors[first])
            if (second > first) bonds.push_back({first, second});
    std::sort(bonds.begin(), bonds.end());
    const bool binnable = std::all_of(coords.begin(), coords.end(), [](const Point2D& point) {
        return std::abs(point.x) < 1e6 && std::abs(point.y) < 1e6;
    });
    if (!binnable) {
        for (int first = 0; first < g.n; ++first)
            for (int second = first + 1; second < g.n; ++second) clash(first, second);
        for (size_t i = 0; i < bonds.size(); ++i)
            for (size_t j = i + 1; j < bonds.size(); ++j) cross(bonds[i], bonds[j]);
        return defects;
    }
    auto cellIndex = [](double value) { return static_cast<int64_t>(std::floor(value)); };
    auto cell = [](int64_t column, int64_t row) {
        return (static_cast<uint64_t>(column) << 32) ^ (static_cast<uint64_t>(row) & 0xffffffffu);
    };
    std::unordered_map<uint64_t, std::vector<int>> atomCells;
    for (int atom = 0; atom < g.n; ++atom) atomCells[cell(cellIndex(coords[atom].x), cellIndex(coords[atom].y))].push_back(atom);
    for (int first = 0; first < g.n; ++first)
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy) {
                const auto found = atomCells.find(cell(cellIndex(coords[first].x) + dx, cellIndex(coords[first].y) + dy));
                if (found == atomCells.end()) continue;
                for (int second : found->second)
                    if (second > first) clash(first, second);
            }
    // Bounding boxes, a little enlarged, and the cells each reaches.
    std::vector<std::array<double, 4>> boxes;
    std::unordered_map<uint64_t, std::vector<size_t>> bondCells;
    for (size_t k = 0; k < bonds.size(); ++k) {
        const Point2D& a = coords[bonds[k].first];
        const Point2D& b = coords[bonds[k].second];
        boxes.push_back({std::min(a.x, b.x) - 1e-6, std::min(a.y, b.y) - 1e-6,
                         std::max(a.x, b.x) + 1e-6, std::max(a.y, b.y) + 1e-6});
        const auto& box = boxes.back();
        for (int64_t column = cellIndex(box[0]); column <= cellIndex(box[2]); ++column)
            for (int64_t row = cellIndex(box[1]); row <= cellIndex(box[3]); ++row) bondCells[cell(column, row)].push_back(k);
    }
    for (const auto& [key, members] : bondCells)
        for (size_t i = 0; i < members.size(); ++i)
            for (size_t j = i + 1; j < members.size(); ++j) {
                const auto& one = boxes[members[i]];
                const auto& other = boxes[members[j]];
                const double left = std::max(one[0], other[0]), top = std::max(one[1], other[1]);
                if (left > std::min(one[2], other[2]) || top > std::min(one[3], other[3])) continue;
                if (cell(cellIndex(left), cellIndex(top)) != key) continue;  // tested where the boxes first meet
                cross(bonds[std::min(members[i], members[j])], bonds[std::max(members[i], members[j])]);
            }
    return defects;
}

/// Atoms reached from `root` without the bond to `pivot`; empty when that bond is in a ring.
inline std::vector<int> branchFrom(const MolGraph& g, int root, int pivot) {
    std::vector<int> atoms{root};
    std::vector<bool> seen(g.n, false);
    seen[root] = true;
    for (size_t k = 0; k < atoms.size(); ++k)
        for (int neighbour : g.neighbors[atoms[k]]) {
            if (neighbour == pivot) {
                if (atoms[k] == root) continue;
                return {};
            }
            if (!seen[neighbour]) { seen[neighbour] = true; atoms.push_back(neighbour); }
        }
    return atoms;
}

/**
 * A coordinate rounded to a grid of 2^-20 bond lengths.  The choices made
 * after Phase 3 rest on rounded coordinates, so that builds which round
 * differently (with or without fused multiply-add, say) seldom reach a
 * different choice, and refinement's many small steps do not drift apart.
 */
inline double onGrid(double value) {
    constexpr double kGrid = 1048576.0; // a power of two, so the scaling is exact
    return std::round(value * kGrid) / kGrid;
}

inline void snapToGrid(std::vector<Point2D>& coords) {
    for (auto& point : coords) point = {onGrid(point.x), onGrid(point.y)};
}

/// Mirror `atoms` in the line through `pivot` and `root`, onto the grid.
inline void reflectInBondAxis(std::vector<Point2D>& coords, const std::vector<int>& atoms, int pivot, int root,
                              double gridScale = 1.0) {
    const Point2D origin = coords[pivot];
    const double axisX = coords[root].x - origin.x, axisY = coords[root].y - origin.y;
    const double lengthSquared = axisX * axisX + axisY * axisY;
    if (lengthSquared < 1e-20) return;
    for (int atom : atoms) {
        const double relativeX = coords[atom].x - origin.x, relativeY = coords[atom].y - origin.y;
        const double projection = (relativeX * axisX + relativeY * axisY) / lengthSquared;
        coords[atom] = {gridScale * onGrid((origin.x + 2 * projection * axisX - relativeX) / gridScale),
                        gridScale * onGrid((origin.y + 2 * projection * axisY - relativeY) / gridScale)};
    }
}

/**
 * Draw every acyclic double bond with stored E/Z as that isomer.  Where the
 * first ligands at its two ends lie on the wrong sides, the smaller side of
 * the bond is mirrored in the bond's axis.  The mirror is an isometry of that
 * side, so it keeps the drawn relation of every other double bond.
 */
inline void honourDoubleBondStereo(const MolGraph& g, std::vector<Point2D>& coords) {
    auto firstLigand = [&](int end, int partner) {
        for (int neighbour : g.neighbors[end])
            if (neighbour != partner) return neighbour;
        return -1;
    };
    for (int first = 0; first < g.n; ++first)
        for (int second : g.neighbors[first]) {
            if (second < first || g.bondOrder(first, second) != 2) continue;
            const int stored = g.dbStereo(first, second);
            const int firstReference = firstLigand(first, second), secondReference = firstLigand(second, first);
            if (stored == 0 || firstReference < 0 || secondReference < 0) continue;
            const double axisX = coords[second].x - coords[first].x, axisY = coords[second].y - coords[first].y;
            const double firstSide = axisX * (coords[firstReference].y - coords[first].y)
                                   - axisY * (coords[firstReference].x - coords[first].x);
            const double secondSide = axisX * (coords[secondReference].y - coords[second].y)
                                    - axisY * (coords[secondReference].x - coords[second].x);
            // A ligand on the axis has no side, and no mirror gives it one.
            if (std::abs(firstSide * secondSide) < 1e-10) continue;
            if ((firstSide * secondSide > 0 ? 1 : 2) == stored) continue;
            auto secondAtoms = branchFrom(g, second, first);
            if (secondAtoms.empty()) continue;
            auto firstAtoms = branchFrom(g, first, second);
            if (firstAtoms.size() < secondAtoms.size()) reflectInBondAxis(coords, firstAtoms, second, first);
            else reflectInBondAxis(coords, secondAtoms, first, second);
        }
}

/**
 * Reflect branches across acyclic single bonds while each reflection lowers
 * the defects: crossings, then atoms hidden by others, then clashes.  The
 * smaller side of the bond is mirrored in the bond's own axis, which is an
 * isometry of that side: ring templates, bond lengths and bond angles all
 * stay exact, and so does the drawn E/Z geometry of every double bond.  The
 * loop ends because each accepted reflection strictly lowers the defects.
 */
inline void reflectClashingBranches(const MolGraph& g, std::vector<Point2D>& coords) {
    auto defects = findDefects(g, coords);
    if (defects.none()) return;

    struct Branch { int pivot, root; std::vector<int> atoms; };
    std::vector<Branch> branches;
    for (int first = 0; first < g.n; ++first)
        for (int second : g.neighbors[first]) {
            if (second < first || g.bondOrder(first, second) != 1) continue;
            auto secondSide = branchFrom(g, second, first);
            if (secondSide.empty()) continue;
            auto firstSide = branchFrom(g, first, second);
            Branch branch = firstSide.size() < secondSide.size() ? Branch{second, first, std::move(firstSide)}
                                                                 : Branch{first, second, std::move(secondSide)};
            // A lone terminal atom lies on the axis, so reflecting it changes nothing.
            if (branch.atoms.size() > 1) branches.push_back(std::move(branch));
        }

    auto reflect = [&](const Branch& branch) { reflectInBondAxis(coords, branch.atoms, branch.pivot, branch.root); };
    for (bool improved = true; improved && !defects.none();) {
        improved = false;
        for (const auto& branch : branches) {
            // A branch free of defects can only gain them by moving.
            if (std::none_of(branch.atoms.begin(), branch.atoms.end(),
                             [&](int atom) { return defects.involved[atom]; })) continue;
            std::vector<Point2D> saved;
            for (int atom : branch.atoms) saved.push_back(coords[atom]);
            reflect(branch);
            auto next = findDefects(g, coords);
            if (next.fewerThan(defects)) {
                defects = std::move(next);
                improved = true;
                if (defects.none()) break;
            } else {
                for (size_t k = 0; k < saved.size(); ++k) coords[branch.atoms[k]] = saved[k];
            }
        }
    }
}

/// Two bonds closer than this to a straight line, 170 degrees, leave stereo unreadable.
constexpr double kStraightAngle = M_PI * 17 / 18;

/// The partner of `atom` in a double bond outside any ring with stored E/Z, or -1.
inline int stereoPartner(const MolGraph& g, int atom) {
    for (int neighbour : g.neighbors[atom])
        if (g.bondOrder(atom, neighbour) == 2 && g.dbStereo(atom, neighbour) != 0
                && !branchFrom(g, neighbour, atom).empty())
            return neighbour;
    return -1;
}

/**
 * True when the drawn bonds at `atom` leave its stored stereo unreadable, each
 * to within kStraightAngle: two bonds in a line at a tetrahedral centre drawn
 * with three, whose hydrogen then has no side (IUPAC ST-1.1); all four bonds
 * of a centre within a half-turn; or a substituent of a stereo double bond in
 * line with that bond (ST-2).
 */
inline bool stereoUnreadable(const MolGraph& g, const std::vector<Point2D>& coords, int atom) {
    auto direction = [&](int neighbour) {
        return std::atan2(coords[neighbour].y - coords[atom].y, coords[neighbour].x - coords[atom].x);
    };
    auto inLine = [&](int first, int second) {
        return std::abs(std::remainder(direction(first) - direction(second), 2 * M_PI)) >= kStraightAngle;
    };
    const auto& neighbours = g.neighbors[atom];
    const bool centre = atom < static_cast<int>(g.tetraChirality.size()) && g.tetraChirality[atom] != 0;
    if (centre && neighbours.size() == 3) {
        if (inLine(neighbours[0], neighbours[1]) || inLine(neighbours[0], neighbours[2]) || inLine(neighbours[1], neighbours[2]))
            return true;
    } else if (centre && neighbours.size() > 3) {
        std::vector<double> angles;
        for (int neighbour : neighbours) angles.push_back(direction(neighbour));
        std::sort(angles.begin(), angles.end());
        double widest = angles.front() + 2 * M_PI - angles.back();
        for (size_t k = 1; k < angles.size(); ++k) widest = std::max(widest, angles[k] - angles[k - 1]);
        if (widest >= kStraightAngle) return true;
    }
    const int partner = stereoPartner(g, atom);
    if (partner < 0) return false;
    return std::any_of(neighbours.begin(), neighbours.end(),
                       [&](int neighbour) { return neighbour != partner && inLine(neighbour, partner); });
}

/**
 * Turn a branch about each atom whose stereo stereoUnreadable finds unreadable,
 * so that its bond bisects the widest gap the other bonds leave, or stands at
 * 120 degrees to a double bond that is the only other bond, or a third or two
 * thirds of the way there.  The branch turns as a rigid body, so no other bond
 * angle changes.  A turn is kept only when it makes the atom readable, keeps
 * the stored E/Z of its double bond and adds no defects; of several, the one
 * with the fewest defects is kept, and of those the furthest turn.
 */
inline void spreadStereoBonds(const MolGraph& g, std::vector<Point2D>& coords) {
    auto direction = [&](int from, int to) {
        return std::atan2(coords[to].y - coords[from].y, coords[to].x - coords[from].x);
    };
    auto firstLigand = [&](int end, int partner) {
        for (int neighbour : g.neighbors[end])
            if (neighbour != partner) return neighbour;
        return -1;
    };
    // As honourDoubleBondStereo reads it; a ligand on the axis is left to its own end.
    auto keepsStoredIsomer = [&](int end, int partner) {
        if (partner < 0) return true;
        const int endLigand = firstLigand(end, partner), partnerLigand = firstLigand(partner, end);
        if (endLigand < 0 || partnerLigand < 0) return true;
        const double axisX = coords[partner].x - coords[end].x, axisY = coords[partner].y - coords[end].y;
        const double endSide = axisX * (coords[endLigand].y - coords[end].y) - axisY * (coords[endLigand].x - coords[end].x);
        const double partnerSide = axisX * (coords[partnerLigand].y - coords[partner].y)
                                 - axisY * (coords[partnerLigand].x - coords[partner].x);
        if (std::abs(endSide * partnerSide) < 1e-10) return true;
        return (endSide * partnerSide > 0 ? 1 : 2) == g.dbStereo(end, partner);
    };
    for (int atom = 0; atom < g.n; ++atom) {
        if (!stereoUnreadable(g, coords, atom)) continue;
        const int partner = stereoPartner(g, atom);
        const auto defects = findDefects(g, coords);
        std::vector<Point2D> best;
        LayoutDefects bestDefects;
        for (int moved : g.neighbors[atom]) {
            if (moved == partner) continue;
            const auto branch = branchFrom(g, moved, atom);
            if (branch.empty()) continue;
            std::vector<double> others;
            for (int other : g.neighbors[atom])
                if (other != moved) others.push_back(direction(atom, other));
            std::sort(others.begin(), others.end());
            std::vector<double> targets;
            if (others.size() == 1) {
                targets = {others.front() + 2 * M_PI / 3, others.front() - 2 * M_PI / 3};
            } else {
                double widest = others.front() + 2 * M_PI - others.back(), target = others.back() + 0.5 * widest;
                for (size_t k = 1; k < others.size(); ++k)
                    if (others[k] - others[k - 1] > widest) {
                        widest = others[k] - others[k - 1];
                        target = others[k - 1] + 0.5 * widest;
                    }
                targets = {target};
            }
            std::vector<double> turns;
            for (double target : targets) {
                const double full = std::remainder(target - direction(atom, moved), 2 * M_PI);
                for (double fraction : {1.0, 2.0 / 3, 1.0 / 3}) turns.push_back(full * fraction);
            }
            for (double turn : turns) {
                const std::vector<Point2D> saved = coords;
                const double cosine = std::cos(turn), sine = std::sin(turn);
                const Point2D pivot = coords[atom];
                for (int member : branch) {
                    const double dx = coords[member].x - pivot.x, dy = coords[member].y - pivot.y;
                    coords[member] = {onGrid(pivot.x + cosine * dx - sine * dy), onGrid(pivot.y + sine * dx + cosine * dy)};
                }
                if (!stereoUnreadable(g, coords, atom) && keepsStoredIsomer(atom, partner)) {
                    auto next = findDefects(g, coords);
                    if (next.noneAdded(defects) && (best.empty() || next.fewerThan(bestDefects))) {
                        best = coords;
                        bestDefects = std::move(next);
                    }
                }
                coords = saved;
            }
        }
        if (!best.empty()) coords = std::move(best);
    }
}

/** Whether the first ligands show the stored relation without collinear bonds. */
inline bool doubleBondStereoMatches(const MolGraph& g, const std::vector<Point2D>& coords,
                                   int first, int second) {
    const double axisX = coords[second].x - coords[first].x;
    const double axisY = coords[second].y - coords[first].y;
    const double axisSquared = axisX * axisX + axisY * axisY;
    if (!std::isfinite(axisSquared) || axisSquared < 1e-20) return false;
    auto referenceSide = [&](int end, int partner) {
        double reference = 0;
        for (int ligand : g.neighbors[end]) {
            if (ligand == partner) continue;
            const double armX = coords[ligand].x - coords[end].x;
            const double armY = coords[ligand].y - coords[end].y;
            const double length = std::sqrt(axisSquared * (armX * armX + armY * armY));
            if (!std::isfinite(length) || length < 1e-20) return 0.0;
            const double side = (axisX * armY - axisY * armX) / length;
            if (std::abs(side) < 1e-6) return 0.0;
            if (reference == 0) reference = side;
            else if (reference * side >= 0) return 0.0;
        }
        return reference;
    };
    const double relation = referenceSide(first, second) * referenceSide(second, first);
    return relation != 0 && (relation > 0 ? 1 : 2) == g.dbStereo(first, second);
}

/**
 * Repair ring E/Z by reflecting a component between two fixed ring pivots.
 * Every moved bond is internal or ends on a pivot, so its length is retained.
 * A move must correct more configured bonds, preserve the ones already correct
 * and add no layout defects. An unsuccessful repair reports failure rather
 * than returning coordinates for the other isomer.
 */
inline void honourRingDoubleBondStereo(const MolGraph& g, std::vector<Point2D>& coords,
                                      const std::vector<RingSystem>& systems,
                                      const std::vector<bool>& fixed = {}) {
    std::vector<std::pair<int, int>> bonds;
    std::vector<size_t> ringBonds;
    double totalLength = 0;
    int bondCount = 0;
    for (int first = 0; first < g.n; ++first)
        for (int second : g.neighbors[first]) {
            if (second <= first) continue;
            totalLength += std::hypot(coords[first].x - coords[second].x,
                                      coords[first].y - coords[second].y);
            ++bondCount;
            if (g.bondOrder(first, second) != 2 || g.dbStereo(first, second) == 0) continue;
            if (branchFrom(g, second, first).empty()) ringBonds.push_back(bonds.size());
            bonds.emplace_back(first, second);
        }
    if (ringBonds.empty()) return;
    const double bondLength = totalLength / bondCount;
    if (!std::isfinite(bondLength) || bondLength < 1e-10)
        throw std::invalid_argument("Cannot lay out stored ring double-bond stereo");
    auto matches = [&]() {
        std::vector<bool> result;
        for (const auto& bond : bonds)
            result.push_back(doubleBondStereoMatches(g, coords, bond.first, bond.second));
        return result;
    };
    auto defects = [&]() {
        auto unit = coords;
        for (auto& point : unit) { point.x /= bondLength; point.y /= bondLength; }
        return findDefects(g, unit);
    };
    auto correct = matches();
    for (bool improved = true; improved;) {
        improved = false;
        for (size_t bondIndex : ringBonds) {
            if (correct[bondIndex]) continue;
            const auto before = coords;
            const auto beforeDefects = defects();
            std::vector<Point2D> best;
            std::vector<bool> bestMatches;
            LayoutDefects bestDefects;
            int bestCorrect = static_cast<int>(std::count(correct.begin(), correct.end(), true));
            double bestDisplacement = std::numeric_limits<double>::infinity();
            const auto bond = bonds[bondIndex];
            for (const auto& ends : {bond, std::make_pair(bond.second, bond.first)})
                for (const auto& system : systems) {
                    if (!system.atoms.count(ends.first) || !system.atoms.count(ends.second)) continue;
                    for (int pivot : system.atoms) {
                        if (pivot == ends.first || pivot == ends.second) continue;
                        std::vector<int> moved{ends.second};
                        std::vector<bool> moving(g.n, false);
                        moving[ends.second] = true;
                        for (size_t at = 0; at < moved.size(); ++at)
                            for (int neighbour : g.neighbors[moved[at]])
                                if (neighbour != ends.first && neighbour != pivot && !moving[neighbour]) {
                                    moving[neighbour] = true;
                                    moved.push_back(neighbour);
                                }
                        if (!fixed.empty() && std::any_of(moved.begin(), moved.end(),
                                [&](int atom) { return fixed[atom]; })) continue;
                        bool splitsRing = false;
                        for (const auto& protectedSystem : systems)
                            for (const auto& ring : protectedSystem.rings) {
                                if (ring.size() >= 8 && !std::all_of(ring.begin(), ring.end(),
                                        [&](int atom) { return g.aromatic[atom]; })) continue;
                                int inside = 0, outside = 0;
                                for (int atom : ring) {
                                    if (atom == ends.first || atom == pivot) continue;
                                    if (moving[atom]) ++inside; else ++outside;
                                }
                                if (inside && outside) splitsRing = true;
                            }
                        if (splitsRing) continue;
                        reflectInBondAxis(coords, moved, ends.first, pivot, bondLength);
                        const auto nextMatches = matches();
                        bool keepsStereo = nextMatches[bondIndex];
                        for (int atom = 0; atom < g.n && keepsStereo; ++atom)
                            for (int neighbour : g.neighbors[atom]) {
                                if (neighbour <= atom) continue;
                                const double previousLength = std::hypot(before[atom].x - before[neighbour].x,
                                                                         before[atom].y - before[neighbour].y);
                                const double nextLength = std::hypot(coords[atom].x - coords[neighbour].x,
                                                                     coords[atom].y - coords[neighbour].y);
                                if (!std::isfinite(nextLength)
                                        || std::abs(nextLength - previousLength) > 1e-4 * bondLength)
                                    keepsStereo = false;
                            }
                        for (size_t index = 0; index < bonds.size(); ++index)
                            if (correct[index] && !nextMatches[index]) keepsStereo = false;
                        for (int atom = 0; atom < g.n && keepsStereo; ++atom)
                            if (atom < static_cast<int>(g.tetraChirality.size()) && g.tetraChirality[atom] != 0
                                    && !stereoUnreadable(g, before, atom)
                                    && stereoUnreadable(g, coords, atom)) keepsStereo = false;
                        if (keepsStereo) {
                            const auto nextDefects = defects();
                            const int nextCorrect = static_cast<int>(std::count(nextMatches.begin(), nextMatches.end(), true));
                            double displacement = 0;
                            for (int atom : moved) {
                                const double dx = coords[atom].x - before[atom].x;
                                const double dy = coords[atom].y - before[atom].y;
                                displacement += dx * dx + dy * dy;
                            }
                            if (nextDefects.noneAdded(beforeDefects) && (nextCorrect > bestCorrect
                                    || (nextCorrect == bestCorrect && !best.empty()
                                        && (nextDefects.fewerThan(bestDefects)
                                            || (!bestDefects.fewerThan(nextDefects) && displacement < bestDisplacement))))) {
                                best = coords;
                                bestMatches = nextMatches;
                                bestDefects = nextDefects;
                                bestCorrect = nextCorrect;
                                bestDisplacement = displacement;
                            }
                        }
                        coords = before;
                    }
                }
            if (!best.empty()) {
                coords = std::move(best);
                correct = std::move(bestMatches);
                improved = true;
            }
        }
    }
    if (std::any_of(ringBonds.begin(), ringBonds.end(), [&](size_t index) { return !correct[index]; }))
        throw std::invalid_argument("Cannot lay out stored ring double-bond stereo");
}

/**
 * Relax a layout in which each rigid group moves as one body.  Bonds pull
 * towards unit length and non-bonded atoms push apart; a free atom follows
 * its net force, while a group shifts with the net force on its atoms and
 * turns with their net torque, so its own geometry changes only by the
 * rounding to the grid.
 */
inline void relaxAroundRigidGroups(const MolGraph& g, std::vector<Point2D>& coords,
                                   const std::vector<int>& group, int iterations) {
    constexpr double kRepulsion = 0.5, kMaxShift = 0.2, kCooling = 0.97;
    const auto members = groupMembers(group);
    double step = 0.1;
    for (int iteration = 0; iteration < iterations; ++iteration, step *= kCooling) {
        snapToGrid(coords);
        std::vector<Point2D> force(g.n, {0, 0});
        for (int first = 0; first < g.n; ++first)
            for (int second = first + 1; second < g.n; ++second) {
                if (group[first] >= 0 && group[first] == group[second]) continue;
                const double dx = coords[second].x - coords[first].x;
                const double dy = coords[second].y - coords[first].y;
                const double distance = std::max(1e-6, std::hypot(dx, dy));
                const double pull = g.bondOrder(first, second) ? distance - 1.0
                                                               : -kRepulsion / (distance * distance);
                force[first].x += pull * dx / distance; force[first].y += pull * dy / distance;
                force[second].x -= pull * dx / distance; force[second].y -= pull * dy / distance;
            }
        auto clamp = [&](double dx, double dy) {
            const double length = std::hypot(dx, dy);
            return length > kMaxShift ? Point2D{dx * kMaxShift / length, dy * kMaxShift / length} : Point2D{dx, dy};
        };
        for (int atom = 0; atom < g.n; ++atom)
            if (group[atom] < 0) {
                const auto shift = clamp(step * force[atom].x, step * force[atom].y);
                coords[atom].x += shift.x; coords[atom].y += shift.y;
            }
        for (const auto& atoms : members) {
            if (atoms.empty()) continue;
            Point2D centre{0, 0}, net{0, 0};
            for (int atom : atoms) {
                centre.x += coords[atom].x; centre.y += coords[atom].y;
                net.x += force[atom].x; net.y += force[atom].y;
            }
            centre.x /= atoms.size(); centre.y /= atoms.size();
            double torque = 0, inertia = 0;
            for (int atom : atoms) {
                const double rx = coords[atom].x - centre.x, ry = coords[atom].y - centre.y;
                torque += rx * force[atom].y - ry * force[atom].x;
                inertia += rx * rx + ry * ry;
            }
            const auto shift = clamp(step * net.x / atoms.size(), step * net.y / atoms.size());
            const double turn = inertia > 0 ? std::clamp(step * torque / inertia, -0.1, 0.1) : 0.0;
            const double cosine = std::cos(turn), sine = std::sin(turn);
            for (int atom : atoms) {
                const double rx = coords[atom].x - centre.x, ry = coords[atom].y - centre.y;
                coords[atom] = {centre.x + shift.x + cosine * rx - sine * ry,
                                centre.y + shift.y + sine * rx + cosine * ry};
            }
        }
    }
}

/**
 * Pull bonds back to unit length after refinement.  Each sweep first parts
 * non-bonded atoms closer than kPartedDistance, so that a bond angle cannot
 * fold shut while the bonds settle.  A rigid group moves as a whole, and each
 * end of a pair gives way in proportion to the other's size.
 */
inline void restoreUnitBonds(const MolGraph& g, std::vector<Point2D>& coords,
                             const std::vector<int>& group) {
    const auto members = groupMembers(group);
    auto size = [&](int atom) { return group[atom] < 0 ? 1.0 : static_cast<double>(members[group[atom]].size()); };
    auto shift = [&](int atom, double dx, double dy) {
        if (group[atom] < 0) { coords[atom].x += dx; coords[atom].y += dy; return; }
        for (int member : members[group[atom]]) { coords[member].x += dx; coords[member].y += dy; }
    };
    auto settle = [&](int first, int second, double target) {
        const double dx = coords[second].x - coords[first].x;
        const double dy = coords[second].y - coords[first].y;
        const double distance = std::hypot(dx, dy);
        if (distance < 1e-10) return;
        const double correction = (distance - target) / distance / (size(first) + size(second));
        shift(first, dx * correction * size(second), dy * correction * size(second));
        shift(second, -dx * correction * size(first), -dy * correction * size(first));
    };
    auto sameGroup = [&](int first, int second) { return group[first] >= 0 && group[first] == group[second]; };
    for (int iteration = 0; iteration < 60; ++iteration) {
        snapToGrid(coords);
        for (int first = 0; first < g.n; ++first)
            for (int second = first + 1; second < g.n; ++second) {
                if (g.bondOrder(first, second) || sameGroup(first, second)) continue;
                if (std::hypot(coords[second].x - coords[first].x, coords[second].y - coords[first].y)
                        < kPartedDistance)
                    settle(first, second, kPartedDistance);
            }
        for (int atom = 0; atom < g.n; ++atom)
            for (int neighbour : g.neighbors[atom])
                if (neighbour > atom && !sameGroup(atom, neighbour)) settle(atom, neighbour, 1.0);
    }
}

} // namespace detail_layout


// ==========================================================================
// generateCoords2D — comprehensive structure diagram generation pipeline
// ==========================================================================

namespace detail_layout {

inline std::vector<Point2D> generateCoords2DImpl(const MolGraph& g, double targetBondLength,
                                               bool requireRingStereo) {
    int n = g.n;
    if (n == 0) return {};
    if (n == 1) return {{0.0, 0.0}};

    std::vector<Point2D> coords(n);
    std::vector<bool> placed(n, false);

    // --- Phase 1: Template matching ---
    auto templ = matchTemplate(g, 1.0); // unit bond length first
    if (templ.size() == static_cast<size_t>(n)) {
        coords = std::move(templ);
        for (int i = 0; i < n; i++) placed[i] = true;
    }

    auto systems = detail_layout::findRingSystems(g);
    std::vector<int> atomSystem(n, -1);
    // The two ring bonds on either side of each atom of every ring.
    std::vector<std::vector<std::pair<int, int>>> ringAngles(n);
    for (int systemIndex = 0; systemIndex < static_cast<int>(systems.size()); ++systemIndex) {
        for (int atom : systems[systemIndex].atoms) atomSystem[atom] = systemIndex;
        for (const auto& ring : systems[systemIndex].rings)
            for (size_t k = 0; k < ring.size(); ++k)
                ringAngles[ring[k]].push_back({ring[(k + ring.size() - 1) % ring.size()], ring[(k + 1) % ring.size()]});
    }
    // --- Phase 2: Ring-first layout ---
    // Each ring system is drawn in a frame of its own, so that an atom two
    // systems share at a spiro centre keeps its true place in both.
    std::vector<std::vector<std::pair<int, Point2D>>> frames(systems.size());
    auto recordFrame = [&](int systemIndex) {
        for (int atom : systems[systemIndex].atoms) frames[systemIndex].push_back({atom, coords[atom]});
    };
    if (placed[0]) {
        for (int systemIndex = 0; systemIndex < static_cast<int>(systems.size()); ++systemIndex)
            recordFrame(systemIndex);
    } else { // not fully placed by template

        // Place ring systems
        for (int systemIndex = 0; systemIndex < static_cast<int>(systems.size()); ++systemIndex) {
            auto& sys = systems[systemIndex];
            if (sys.rings.empty()) continue;
            for (int atom : sys.atoms) placed[atom] = false;

            // Pick the largest ring as seed (prefer 6-membered, then largest)
            int seedIdx = 0;
            int bestScore = 0;
            for (int i = 0; i < static_cast<int>(sys.rings.size()); i++) {
                int sz = static_cast<int>(sys.rings[i].size());
                int score = sz * 10;
                if (sz == 6) score += 100; // prefer hexagonal
                if (score > bestScore) { bestScore = score; seedIdx = i; }
            }

            // Place seed ring as regular polygon
            auto& seedRing = sys.rings[seedIdx];
            int rSize = static_cast<int>(seedRing.size());
            std::vector<Point2D> polyCoords(rSize);
            detail_layout::regularPolygon(rSize, polyCoords.data());
            for (int i = 0; i < rSize; i++) {
                coords[seedRing[i]] = polyCoords[i];
                placed[seedRing[i]] = true;
            }

            // Fuse remaining rings to already-placed atoms
            std::vector<bool> ringPlaced(sys.rings.size(), false);
            ringPlaced[seedIdx] = true;
            bool progress = true;
            while (progress) {
                progress = false;
                for (int ri = 0; ri < static_cast<int>(sys.rings.size()); ri++) {
                    if (ringPlaced[ri]) continue;
                    auto& ring = sys.rings[ri];

                    // Find shared edge with an already-placed ring
                    int shA = -1, shB = -1, sharedIndex = -1;
                    for (size_t k = 0; k < ring.size(); k++) {
                        int a = ring[k], b = ring[(k + 1) % ring.size()];
                        if (placed[a] && placed[b]) {
                            shA = a; shB = b; sharedIndex = static_cast<int>(k);
                            break;
                        }
                    }
                    if (shA < 0) continue; // no shared edge yet

                    // Place this ring fused to shared edge
                    int ringSize = static_cast<int>(ring.size());
                    double rB = 0.5 / std::sin(M_PI / ringSize);
                    Point2D pA = coords[shA], pB = coords[shB];
                    double mx = (pA.x + pB.x) * 0.5, my = (pA.y + pB.y) * 0.5;
                    double edx = pB.x - pA.x, edy = pB.y - pA.y;
                    double nx = -edy, ny = edx;

                    // Place the new polygon opposite the existing neighbours of the shared edge.
                    double cax = 0, cay = 0;
                    int placedCount = 0;
                    for (int sharedAtom : {shA, shB})
                        for (int neighbour : g.neighbors[sharedAtom])
                            if (placed[neighbour] && std::find(ring.begin(), ring.end(), neighbour) == ring.end()) {
                                cax += coords[neighbour].x; cay += coords[neighbour].y; ++placedCount;
                            }
                    if (placedCount && nx * (cax / placedCount - mx) + ny * (cay / placedCount - my) > 0) {
                        nx = -nx; ny = -ny;
                    }
                    const double nmag = std::hypot(nx, ny);
                    const double apothem = rB * std::cos(M_PI / ringSize);
                    const double cx = mx + apothem * nx / nmag;
                    const double cy = my + apothem * ny / nmag;
                    const double startAngle = std::atan2(pB.y - cy, pB.x - cx);
                    const double direction = ((pA.x - cx) * (pB.y - cy)
                                            - (pA.y - cy) * (pB.x - cx)) > 0 ? 1.0 : -1.0;
                    for (int offset = 1; offset <= ringSize - 2; ++offset) {
                        const int atom = ring[(sharedIndex + 1 + offset) % ringSize];
                        if (placed[atom]) continue;
                        const double angle = startAngle + direction * 2.0 * M_PI * offset / ringSize;
                        coords[atom] = {cx + rB * std::cos(angle), cy + rB * std::sin(angle)};
                        placed[atom] = true;
                    }

                    ringPlaced[ri] = true;
                    progress = true;
                }
            }
            recordFrame(systemIndex);
        }
    }

    // --- Phase 3: Attach ring systems and chains along the molecular graph. ---
    std::fill(placed.begin(), placed.end(), false);
    std::deque<int> bfsQ;
    // Places `atom` at `position` with the rest of its ring system, if it has
    // one, turned so that the system's outer side at `atom` faces back along `angle`.
    auto attach = [&](int atom, Point2D position, double angle, int systemIndex) {
        if (systemIndex < 0) {
            coords[atom] = position;
            placed[atom] = true;
            bfsQ.push_back(atom);
            return;
        }
        const auto& atoms = systems[systemIndex].atoms;
        std::vector<Point2D> frame(n);
        for (const auto& [member, point] : frames[systemIndex]) frame[member] = point;
        std::vector<int> ringNeighbours;
        for (int neighbour : g.neighbors[atom])
            if (atoms.count(neighbour)) ringNeighbours.push_back(neighbour);
        const double outward = detail_layout::clearestGapBisector(
            frame, atom, ringNeighbours, std::vector<int>(atoms.begin(), atoms.end()), ringAngles[atom]);
        const double rotation = angle + M_PI - outward;
        for (int member : atoms) {
            const double dx = frame[member].x - frame[atom].x;
            const double dy = frame[member].y - frame[atom].y;
            coords[member] = {position.x + dx * std::cos(rotation) - dy * std::sin(rotation),
                              position.y + dx * std::sin(rotation) + dy * std::cos(rotation)};
            placed[member] = true;
            bfsQ.push_back(member);
        }
    };
    auto placedAtoms = [&]() {
        std::vector<int> atoms;
        for (int existing = 0; existing < n; ++existing)
            if (placed[existing]) atoms.push_back(existing);
        return atoms;
    };
    // Turns the ring system just attached at the spiro atom `pivot` by up to 30 degrees either
    // way, in 5-degree steps, to where its atoms lie furthest from `others`, scored as in
    // clearestGapBisector.
    auto turnClear = [&](int pivot, int systemIndex, const std::vector<int>& others) {
        const auto& members = systems[systemIndex].atoms;
        std::vector<int> fixed;
        for (int other : others)
            if (!members.count(other)) fixed.push_back(other);
        const std::vector<Point2D> start = coords;
        auto crowding = [&](double turn) {
            const double cosine = std::cos(turn), sine = std::sin(turn);
            double total = 0;
            for (int member : members) {
                const double dx = start[member].x - start[pivot].x, dy = start[member].y - start[pivot].y;
                coords[member] = {start[pivot].x + cosine * dx - sine * dy, start[pivot].y + sine * dx + cosine * dy};
                for (int other : fixed) {
                    const double ox = coords[member].x - coords[other].x, oy = coords[member].y - coords[other].y;
                    total += 1.0 / std::max(0.01, ox * ox + oy * oy);
                }
            }
            return total;
        };
        double bestTurn = 0, least = crowding(0);
        for (int step = -6; step <= 6; ++step) {
            const double crowd = crowding(step * M_PI / 36);
            if (detail_layout::clearlyHigher(least, crowd)) { least = crowd; bestTurn = step * M_PI / 36; }
        }
        crowding(bestTurn);
    };
    double componentRight = -3.0;
    for (int seed = 0; seed < n; ++seed) {
        if (placed[seed]) continue;
        attach(seed, {componentRight + 3.0, 0}, 0, atomSystem[seed]);
        while (!bfsQ.empty()) {
            const int atom = bfsQ.front(); bfsQ.pop_front();
            for (int neighbour : g.neighbors[atom]) {
                if (placed[neighbour]) continue;
                std::vector<int> placedNeighbours;
                for (int existing : g.neighbors[atom])
                    if (placed[existing]) placedNeighbours.push_back(existing);
                const int system = atomSystem[neighbour];
                if (system >= 0 && systems[system].atoms.count(atom)) {
                    // A spiro centre: the second ring system opens into the clearest gap, then
                    // turns within it to where it lies clearest.
                    const auto others = placedAtoms();
                    const double angle = detail_layout::clearestGapBisector(coords, atom, placedNeighbours,
                                                                            others, ringAngles[atom]);
                    attach(atom, coords[atom], angle, system);
                    turnClear(atom, system, others);
                    continue;
                }
                if (placedNeighbours.size() >= 3) {
                    // With three bonds already drawn, 120-degree angles are out of reach:
                    // the new bond bisects the clearest gap (an angular methyl, say).
                    const double angle = detail_layout::clearestGapBisector(coords, atom, placedNeighbours,
                                                                            placedAtoms(), ringAngles[atom],
                                                                            g.neighbors[neighbour].size() == 1);
                    attach(neighbour, {coords[atom].x + std::cos(angle), coords[atom].y + std::sin(angle)}, angle,
                           system);
                    continue;
                }
                double bestAngle = 0, bestScore = 0;
                double reference = 0;
                if (!placedNeighbours.empty())
                    reference = std::atan2(coords[placedNeighbours.front()].y - coords[atom].y,
                                           coords[placedNeighbours.front()].x - coords[atom].x);
                for (int directionIndex = 0; directionIndex < 12; ++directionIndex) {
                    const double angle = reference + directionIndex * M_PI / 6.0;
                    const Point2D candidate{coords[atom].x + std::cos(angle), coords[atom].y + std::sin(angle)};
                    double score = 0;
                    for (int existing = 0; existing < n; ++existing) {
                        if (!placed[existing] || existing == atom) continue;
                        const double distance = std::hypot(candidate.x - coords[existing].x,
                                                            candidate.y - coords[existing].y);
                        score -= 1.0 / std::max(0.01, distance * distance);
                    }
                    for (int existing : g.neighbors[atom])
                        if (placed[existing]) {
                            const double existingAngle = std::atan2(coords[existing].y - coords[atom].y,
                                                                    coords[existing].x - coords[atom].x);
                            const double separation = std::abs(std::remainder(angle - existingAngle, 2 * M_PI));
                            score -= 4.0 * std::abs(separation - 2 * M_PI / 3);
                        }
                    if (directionIndex == 0 || detail_layout::clearlyHigher(score, bestScore)) {
                        bestScore = score;
                        bestAngle = angle;
                    }
                }
                attach(neighbour, {coords[atom].x + std::cos(bestAngle), coords[atom].y + std::sin(bestAngle)},
                       bestAngle, system);
            }
        }
        for (int atom = 0; atom < n; ++atom)
            if (placed[atom]) componentRight = std::max(componentRight, coords[atom].x);
    }

    // --- Phase 4: Stored E/Z, then rigid clash removal ---
    detail_layout::snapToGrid(coords);
    detail_layout::honourDoubleBondStereo(g, coords);
    detail_layout::reflectClashingBranches(g, coords);

    if (!detail_layout::isCleanLayout(g, coords)) {
        // --- Phase 5: Refinement around rigid ring systems ---
        const auto rigidGroup = detail_layout::rigidRingGroups(g, coords, systems);
        detail_layout::relaxAroundRigidGroups(g, coords, rigidGroup, 100);
        // A reflection that removes a crossing can bring two atoms together, so the bonds
        // settle again, up to twice more, while an atom lies hidden.
        for (int pass = 0; pass < 3; ++pass) {
            detail_layout::restoreUnitBonds(g, coords, rigidGroup);
            detail_layout::snapToGrid(coords);
            detail_layout::honourDoubleBondStereo(g, coords);
            detail_layout::reflectClashingBranches(g, coords);
            if (detail_layout::findDefects(g, coords).hidden == 0) break;
        }
    }
    detail_layout::spreadStereoBonds(g, coords);
    if (requireRingStereo) detail_layout::honourRingDoubleBondStereo(g, coords, systems);

    // --- Phase 6: Canonical orientation ---
    transform::canonicalOrientation(g, coords);

    // --- Phase 7: Bond length normalisation ---
    transform::normaliseBondLength(g, coords, targetBondLength);
    if (std::all_of(g.neighbors.begin(), g.neighbors.end(),
                    [](const auto& neighbours) { return neighbours.empty(); }))
        transform::scale2D(coords, targetBondLength);

    return coords;
}

} // namespace detail_layout

/**
 * Generate publication-quality 2D coordinates for a molecule.
 *
 * Multi-phase pipeline:
 *
 * Phase 1: Template matching — try extended template library for known scaffolds.
 * Phase 2: Ring-first layout — place ring systems as regular polygons, fuse
 *          adjacent rings by shared edges.
 * Phase 3: Chain layout — place acyclic branches in zig-zag convention.
 * Phase 4: Stored E/Z, then rigid clash removal — reflect branches across
 *          acyclic single bonds.
 * Phase 5: Refinement, only when crossings, clashes or uneven bonds remain —
 *          exact ring systems move as rigid bodies, so they stay regular.
 *          Afterwards, in every layout, branches at stereocentres and stereo
 *          double bonds turn where two bonds would lie in a line, so that the
 *          drawn stereo stays readable.
 * Phase 6: Canonical orientation — principal axis horizontal, then bonds
 *          turned onto the 30-degree grid as nearly as they will go.
 * Phase 7: Bond length normalisation — scale to target bond length.
 *
 * @param g                The molecular graph.
 * @param targetBondLength Desired bond length in coordinate units (default 1.5).
 * @return                 2D coordinates for each atom, crossing-free where possible.
 * @throws std::invalid_argument if no layout showing the stored ring double-bond stereo is found.
 */
inline std::vector<Point2D> generateCoords2D(const MolGraph& g,
                                          double targetBondLength = 1.5) {
    return detail_layout::generateCoords2DImpl(g, targetBondLength, true);
}


// ==========================================================================
// generateCoords3D — basic 3D coordinate generation via distance geometry
// ==========================================================================

/**
 * Generate initial 3D coordinates using distance geometry embedding.
 *
 * Uses a simplified distance-bounds + eigendecomposition approach:
 * 1. Compute graph distances as lower bound estimates.
 * 2. Set upper bounds using van der Waals radii.
 * 3. Sample random distances within bounds.
 * 4. Embed in 3D via classical MDS (eigenvalue decomposition of Gram matrix).
 * 5. Refine with a short force-field pass.
 *
 * For production 3D conformer generation, external tools (ETKDG, MMFF)
 * are recommended. This provides a reasonable starting geometry.
 *
 * @param g                The molecular graph.
 * @param targetBondLength Desired bond length (default 1.5 Angstroms).
 * @return                 3D coordinates for each atom.
 */
inline std::vector<Point3D> generateCoords3D(const MolGraph& g,
                                              double targetBondLength = 1.5) {
    int n = g.n;
    if (n == 0) return {};
    if (n == 1) return {{0.0, 0.0, 0.0}};

    // Distance geometry: use graph distances * bond length as target
    auto graphDist = detail_layout::bfsDistanceMatrix(g);

    // Target distance matrix
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            D[i][j] = static_cast<double>(graphDist[i][j]) * targetBondLength;

    // Classical MDS: embed in 3D from distance matrix
    // 1. Double centring: B = -0.5 * H * D^2 * H where H = I - (1/n)*11'
    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
    {
        std::vector<double> rowMean(n, 0.0);
        double grandMean = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double d2 = D[i][j] * D[i][j];
                rowMean[i] += d2;
            }
            rowMean[i] /= n;
            grandMean += rowMean[i];
        }
        grandMean /= n;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                B[i][j] = -0.5 * (D[i][j] * D[i][j] - rowMean[i] - rowMean[j] + grandMean);
    }

    // 2. Power iteration to find top 3 eigenvectors
    std::vector<Point3D> coords(n);
    // Store unscaled unit eigenvectors for correct deflation
    std::vector<std::vector<double>> eigenvecs;
    std::mt19937 rng(42);
    std::normal_distribution<double> normal(0.0, 1.0);

    auto powerIteration = [&](int dim) -> std::pair<double, std::vector<double>> {
        std::vector<double> v(n);
        for (int i = 0; i < n; i++) v[i] = normal(rng);

        double eigenval = 0.0;
        for (int iter = 0; iter < 100; iter++) {
            // Matrix-vector multiply
            std::vector<double> Bv(n, 0.0);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    Bv[i] += B[i][j] * v[j];

            // Deflate against previous unscaled eigenvectors
            for (int d = 0; d < dim; d++) {
                double dotPrev = 0.0;
                for (int i = 0; i < n; i++)
                    dotPrev += Bv[i] * eigenvecs[d][i];
                for (int i = 0; i < n; i++)
                    Bv[i] -= dotPrev * eigenvecs[d][i];
            }

            // Normalise
            double norm = 0.0;
            for (int i = 0; i < n; i++) norm += Bv[i] * Bv[i];
            norm = std::sqrt(norm);
            if (norm < 1e-15) break;
            eigenval = norm;
            for (int i = 0; i < n; i++) v[i] = Bv[i] / norm;
        }
        return {eigenval, v};
    };

    // Extract top 3 dimensions
    for (int dim = 0; dim < 3; dim++) {
        auto [eigenval, eigenvec] = powerIteration(dim);
        eigenvecs.push_back(eigenvec);  // store unscaled for deflation
        double scale = (eigenval > 0) ? std::sqrt(eigenval) : 0.0;
        for (int i = 0; i < n; i++) {
            double val = eigenvec[i] * scale;
            if (dim == 0)      coords[i].x = val;
            else if (dim == 1) coords[i].y = val;
            else               coords[i].z = val;
        }
    }

    // 3. Short force-field refinement in 3D
    for (int iter = 0; iter < 100; iter++) {
        std::vector<double> fx(n, 0.0), fy(n, 0.0), fz(n, 0.0);

        // Bond-length forces
        for (int i = 0; i < n; i++) {
            for (int j : g.neighbors[i]) {
                if (j <= i) continue;
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                double dz = coords[j].z - coords[i].z;
                double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                if (dist < 1e-10) dist = 1e-10;
                double force = (dist - targetBondLength) * 0.5;
                double ux = dx / dist, uy = dy / dist, uz = dz / dist;
                fx[i] += force * ux; fy[i] += force * uy; fz[i] += force * uz;
                fx[j] -= force * ux; fy[j] -= force * uy; fz[j] -= force * uz;
            }
        }

        // Non-bonded repulsion
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                double dz = coords[j].z - coords[i].z;
                double d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < 1e-6) d2 = 1e-6;
                if (d2 < targetBondLength * targetBondLength * 4) {
                    double force = -0.3 / d2;
                    double dist = std::sqrt(d2);
                    double ux = dx / dist, uy = dy / dist, uz = dz / dist;
                    fx[i] += force * ux; fy[i] += force * uy; fz[i] += force * uz;
                    fx[j] -= force * ux; fy[j] -= force * uy; fz[j] -= force * uz;
                }
            }
        }

        double step = 0.05 * std::pow(0.98, iter);
        for (int i = 0; i < n; i++) {
            coords[i].x += step * fx[i];
            coords[i].y += step * fy[i];
            coords[i].z += step * fz[i];
        }
    }

    return coords;
}


// ==========================================================================
// Layout quality metrics
// ==========================================================================

/// Compute layout quality score (0.0 = perfect, higher = worse).
/// Combines: bond length uniformity, overlap count, crossing count.
inline double layoutQuality(const MolGraph& g, const std::vector<Point2D>& coords,
                            double targetBondLength = 1.5) {
    if (g.n < 2) return 0.0;
    int n = g.n;

    // Bond length uniformity: stddev / target
    double sumLen = 0.0, sumLen2 = 0.0;
    int nBonds = 0;
    for (int i = 0; i < n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                double dx = coords[j].x - coords[i].x;
                double dy = coords[j].y - coords[i].y;
                double len = std::sqrt(dx * dx + dy * dy);
                sumLen += len;
                sumLen2 += len * len;
                nBonds++;
            }
        }
    }
    double bondScore = 0.0;
    if (nBonds > 0) {
        double mean = sumLen / nBonds;
        double variance = sumLen2 / nBonds - mean * mean;
        if (variance < 0) variance = 0;
        bondScore = std::sqrt(variance) / targetBondLength;
    }

    // Overlap penalty
    double overlapScore = static_cast<double>(detail_layout::countOverlaps(coords, 0.3 * targetBondLength));

    // Crossing penalty
    double crossingScore = static_cast<double>(detail_layout::countCrossings(g, coords));

    return bondScore + overlapScore * 5.0 + crossingScore * 10.0;
}


} // namespace smsd

#endif // SMSD_LAYOUT_HPP

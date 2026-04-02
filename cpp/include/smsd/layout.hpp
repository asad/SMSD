/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
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

    int nr = static_cast<int>(sssr.size());

    // Build edge sets per ring
    std::vector<std::set<int64_t>> edgeSets(nr);
    for (int i = 0; i < nr; i++) {
        auto& r = sssr[i];
        for (size_t j = 0; j < r.size(); j++) {
            int a = r[j], b = r[(j + 1) % r.size()];
            int lo = std::min(a, b), hi = std::max(a, b);
            edgeSets[i].insert((static_cast<int64_t>(lo) << 32) |
                               static_cast<int64_t>(static_cast<uint32_t>(hi)));
        }
    }

    // Ring adjacency
    std::vector<std::vector<int>> adj(nr);
    for (int i = 0; i < nr; i++) {
        for (int j = i + 1; j < nr; j++) {
            for (auto& e : edgeSets[i]) {
                if (edgeSets[j].count(e)) {
                    adj[i].push_back(j);
                    adj[j].push_back(i);
                    break;
                }
            }
        }
    }

    // Connected components
    std::vector<int> comp(nr, -1);
    int nComp = 0;
    for (int i = 0; i < nr; i++) {
        if (comp[i] >= 0) continue;
        int c = nComp++;
        comp[i] = c;
        std::vector<int> stk = {i};
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            for (int v : adj[u]) {
                if (comp[v] < 0) { comp[v] = c; stk.push_back(v); }
            }
        }
    }

    std::vector<RingSystem> systems(nComp);
    for (int i = 0; i < nr; i++) {
        systems[comp[i]].rings.push_back(sssr[i]);
        for (int a : sssr[i]) systems[comp[i]].atoms.insert(a);
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
 *                 Modified in-place with optimized positions.
 * @param maxIter  Maximum total SA iterations across both phases (default 1000).
 * @return         Number of bond crossings remaining after optimization.
 * @since 6.5.1
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
 * @since 6.3.3
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
        // Update positions with PrEd-style clamping (Bertault 2000):
        // Limit displacement to avoid introducing NEW crossings.
        int crossingsBefore = detail_layout::countCrossings(g, coords);
        double maxDisp = 0.0;
        for (int i = 0; i < n; i++) {
            double dx = step * fx[i];
            double dy = step * fy[i];
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
            double disp = std::sqrt(dx * dx + dy * dy);
            if (disp > maxDisp) maxDisp = disp;
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
 * @since 6.4.0
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

    // Compute stress
    auto computeStress = [&]() -> double {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double dx = coords[i].x - coords[j].x;
                double dy = coords[i].y - coords[j].y;
                double dij = std::sqrt(dx * dx + dy * dy);
                double delta = dij - D[i][j];
                s += W[i][j] * delta * delta;
            }
        }
        return s;
    };

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
// Coordinates are for a canonical atom ordering; real matching requires
// subgraph isomorphism which is handled via the existing VF2PP engine.
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

inline const std::vector<ScaffoldTemplate>& getTemplates() {
    static std::vector<ScaffoldTemplate> templates;
    static bool initialised = false;
    if (initialised) return templates;
    initialised = true;

    // Helper: fused-ring system builder.
    // Place ring A as a regular polygon, then ring B sharing an edge.
    auto fusedRings = [](int sizeA, int sizeB, int sharedEdgeA0, int sharedEdgeA1,
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
        // If normal points toward ring A centre, flip it
        double dot = nx * (cax - mx) + ny * (cay - my);
        if (dot > 0) { nx = -nx; ny = -ny; }
        // Place ring B as regular polygon, then transform
        double rB = 0.5 / std::sin(M_PI / sizeB);
        double edgeLen = std::sqrt(edx * edx + edy * edy);
        // Rotation angle of shared edge
        double theta = std::atan2(edy, edx);
        // Centre of ring B
        double centreDistFromMidpoint = rB * std::cos(M_PI / sizeB);
        double cx = mx + centreDistFromMidpoint * nx / std::sqrt(nx * nx + ny * ny);
        double cy = my + centreDistFromMidpoint * ny / std::sqrt(nx * nx + ny * ny);
        (void)edgeLen;
        // Place remaining sizeB - 2 atoms of ring B
        for (int k = 0; k < sizeB - 2; k++) {
            // Angle offset: start from just past the shared edge
            double angle = theta + M_PI + (2.0 * M_PI * (k + 1)) / sizeB;
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
        // Second ring offset by bond length
        double offsetX = t.coords[0].x - t.coords[3].x + 1.0;
        for (int i = 0; i < 6; i++) {
            double angle = 2.0 * M_PI * i / 6 - M_PI / 2.0;
            double r = 0.5 / std::sin(M_PI / 6.0);
            t.coords[6 + i].x = r * std::cos(angle) + offsetX + 1.5;
            t.coords[6 + i].y = r * std::sin(angle);
        }
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

} // namespace detail_layout


// ==========================================================================
// matchTemplate — scaffold template matching for common ring systems
// ==========================================================================

/**
 * Check if a molecule matches a known scaffold template.
 * Returns pre-computed coordinates scaled to targetBondLength if a match
 * is found, empty vector otherwise.
 *
 * Matching is based on ring-size signature (sorted list of SSSR ring sizes)
 * and atom count. For molecules with matching topology, the template
 * coordinates provide optimal crossing-free layout.
 *
 * @param g                The molecular graph.
 * @param targetBondLength Desired bond length for scaling (default 1.5).
 * @return                 Pre-computed coordinates, or empty if no match.
 */
inline std::vector<Point2D> matchTemplate(const MolGraph& g,
                                           double targetBondLength = 1.5) {
    auto sig = detail_layout::ringSizeSignature(g);
    auto& templates = detail_layout::getTemplates();

    for (auto& t : templates) {
        // Compare ring-size signature
        std::vector<int> tSig;
        for (int k = 0; k < 8 && t.ringSizes[k] != 0; k++)
            tSig.push_back(t.ringSizes[k]);
        std::sort(tSig.begin(), tSig.end());
        if (tSig != sig) continue;
        if (t.atomCount != g.n) continue;

        // Match found — return scaled coordinates
        std::vector<Point2D> result(t.atomCount);
        for (int i = 0; i < t.atomCount; i++) {
            result[i].x = t.coords[i].x * targetBondLength;
            result[i].y = t.coords[i].y * targetBondLength;
        }
        return result;
    }
    return {};
}


} // namespace smsd

#endif // SMSD_LAYOUT_HPP

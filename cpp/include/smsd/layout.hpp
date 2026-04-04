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

/// Center coordinates at origin.
inline void center2D(std::vector<Point2D>& coords) {
    if (coords.empty()) return;
    double cx = 0, cy = 0;
    for (auto& p : coords) { cx += p.x; cy += p.y; }
    cx /= static_cast<double>(coords.size());
    cy /= static_cast<double>(coords.size());
    for (auto& p : coords) { p.x -= cx; p.y -= cy; }
}

/// Center 3D coordinates at origin.
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

/// Canonical orientation: rotate to align longest axis with horizontal,
/// then snap to nearest 30-degree increment (IUPAC convention).
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

    // Principal axis angle from the 2x2 covariance matrix
    double angle = 0.5 * std::atan2(2.0 * sxy, sxx - syy);

    // Rotate so principal axis is horizontal
    rotate2DAbout(coords, -angle, cx, cy);

    // Snap to nearest 30-degree increment (IUPAC: M_PI/6)
    double snapAngle = std::round(angle / (M_PI / 6.0)) * (M_PI / 6.0);
    double correction = snapAngle - angle;
    if (std::abs(correction) > 1e-6) {
        rotate2DAbout(coords, correction, cx, cy);
    }
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
// Chain layout — zig-zag convention for acyclic chains
// ==========================================================================

namespace detail_layout {

/// Place acyclic chain atoms in zig-zag pattern.
/// sp3 bonds use 109.5 angle, sp2 bonds use 120 angle, sp use 180.
/// Direction alternates up/down for each bond.
inline void zigzagChain(const MolGraph& g, const std::vector<int>& chain,
                        std::vector<Point2D>& coords, double bondLen,
                        double startAngle = 0.0) {
    if (chain.empty()) return;
    if (chain.size() == 1) return;

    double angle = startAngle;
    int direction = 1; // +1 or -1 for zig-zag

    for (size_t k = 1; k < chain.size(); k++) {
        int prev = chain[k - 1];
        int cur  = chain[k];

        // Determine bond angle based on hybridisation
        int bo = g.bondOrder(prev, cur);
        double bendAngle;
        if (bo == 3) {
            bendAngle = 0.0;  // sp: linear 180
        } else if (bo == 2 || g.aromatic[prev]) {
            bendAngle = M_PI / 6.0;  // sp2: 120 = pi - pi/3 -> bend pi/6 from straight
        } else {
            bendAngle = (M_PI - 109.5 * M_PI / 180.0) / 2.0; // sp3: 109.5 zig-zag
        }

        angle += direction * bendAngle;
        coords[cur].x = coords[prev].x + bondLen * std::cos(angle);
        coords[cur].y = coords[prev].y + bondLen * std::sin(angle);
        direction = -direction;
    }
}

} // namespace detail_layout


// ==========================================================================
// Expanded scaffold template library — 40+ pharmaceutical scaffolds
// ==========================================================================

namespace detail_layout {

/// Extended template library with comprehensive pharmaceutical scaffolds.
/// Coordinates computed analytically using regular polygon fusion.
inline const std::vector<ScaffoldTemplate>& getExtendedTemplates() {
    static std::vector<ScaffoldTemplate> templates;
    static bool initialised = false;
    if (initialised) return templates;
    initialised = true;

    // Reuse fusedRings helper
    auto fusedRings = [](int sizeA, int sizeB, int sharedEdgeA0, int sharedEdgeA1,
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
        double theta = std::atan2(edy, edx);
        double nmag = std::sqrt(nx * nx + ny * ny);
        double centreDistFromMidpoint = rB * std::cos(M_PI / sizeB);
        double cx = mx + centreDistFromMidpoint * nx / nmag;
        double cy = my + centreDistFromMidpoint * ny / nmag;
        for (int k = 0; k < sizeB - 2; k++) {
            double angle = theta + M_PI + (2.0 * M_PI * (k + 1)) / sizeB;
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

} // namespace detail_layout


// ==========================================================================
// generateCoords2D — comprehensive structure diagram generation pipeline
// ==========================================================================

/**
 * Generate publication-quality 2D coordinates for a molecule.
 *
 * Multi-phase pipeline inspired by CoordGen and CDK SDG:
 *
 * Phase 1: Template matching — try extended template library for known scaffolds.
 * Phase 2: Ring-first layout — place ring systems as regular polygons, fuse
 *          adjacent rings by shared edges.
 * Phase 3: Chain layout — place acyclic branches in zig-zag convention.
 * Phase 4: Force refinement — short force-directed pass for local optimisation.
 * Phase 5: Overlap resolution — push apart any remaining clashes.
 * Phase 6: Crossing reduction — simulated annealing ring flips.
 * Phase 7: Canonical orientation — align to 30-degree grid.
 * Phase 8: Bond length normalisation — scale to target bond length.
 *
 * @param g                The molecular graph.
 * @param targetBondLength Desired bond length in coordinate units (default 1.5).
 * @return                 2D coordinates for each atom, crossing-free where possible.
 * @since 6.11.0
 */
inline std::vector<Point2D> generateCoords2D(const MolGraph& g,
                                              double targetBondLength = 1.5) {
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

    // --- Phase 2: Ring-first layout ---
    if (!placed[0]) { // not fully placed by template
        auto systems = detail_layout::findRingSystems(g);

        // Place ring systems
        for (auto& sys : systems) {
            if (sys.rings.empty()) continue;

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
                    int shA = -1, shB = -1;
                    for (size_t k = 0; k < ring.size(); k++) {
                        int a = ring[k], b = ring[(k + 1) % ring.size()];
                        if (placed[a] && placed[b]) {
                            shA = a; shB = b;
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

                    // Determine outward direction
                    double cax = 0, cay = 0;
                    int placedCount = 0;
                    for (int a : ring) {
                        if (placed[a]) { cax += coords[a].x; cay += coords[a].y; placedCount++; }
                    }
                    if (placedCount > 0) {
                        cax /= placedCount; cay /= placedCount;
                        double dot = nx * (cax - mx) + ny * (cay - my);
                        if (dot > 0) { nx = -nx; ny = -ny; }
                    }

                    double nmag = std::sqrt(nx * nx + ny * ny);
                    if (nmag < 1e-10) nmag = 1.0;
                    double centreDistFromMidpoint = rB * std::cos(M_PI / ringSize);
                    double cx = mx + centreDistFromMidpoint * nx / nmag;
                    double cy = my + centreDistFromMidpoint * ny / nmag;
                    double theta = std::atan2(edy, edx);

                    // Place unplaced atoms
                    int atomSlot = 0;
                    for (size_t k = 0; k < ring.size(); k++) {
                        if (!placed[ring[k]]) {
                            atomSlot++;
                            double angle = theta + M_PI + (2.0 * M_PI * atomSlot) / ringSize;
                            coords[ring[k]].x = cx + rB * std::cos(angle);
                            coords[ring[k]].y = cy + rB * std::sin(angle);
                            placed[ring[k]] = true;
                        }
                    }

                    ringPlaced[ri] = true;
                    progress = true;
                }
            }
        }
    }

    // --- Phase 3: Chain layout (BFS from placed atoms) ---
    // BFS to place remaining atoms attached to already-placed atoms
    std::deque<int> bfsQ;
    for (int i = 0; i < n; i++) {
        if (placed[i]) bfsQ.push_back(i);
    }
    // Handle disconnected molecules: ensure at least one atom is placed
    if (bfsQ.empty() && n > 0) {
        coords[0] = {0.0, 0.0};
        placed[0] = true;
        bfsQ.push_back(0);
    }

    while (!bfsQ.empty()) {
        int u = bfsQ.front(); bfsQ.pop_front();
        for (int v : g.neighbors[u]) {
            if (placed[v]) continue;

            // Determine bond angle: look at already-placed neighbors of u
            // and pick an angle that avoids collisions
            double baseAngle = 0.0;
            int placedNbCount = 0;
            double sumAngle = 0.0;
            for (int w : g.neighbors[u]) {
                if (w != v && placed[w]) {
                    double dx = coords[w].x - coords[u].x;
                    double dy = coords[w].y - coords[u].y;
                    sumAngle += std::atan2(dy, dx);
                    placedNbCount++;
                }
            }

            if (placedNbCount == 0) {
                baseAngle = 0.0;
            } else if (placedNbCount == 1) {
                // Place at 120 degrees from existing neighbor
                double existAngle = sumAngle;
                // Determine bend direction based on hybridisation
                int bo = g.bondOrder(u, v);
                double bendAngle = (bo >= 2) ? (2.0 * M_PI / 3.0) : (109.5 * M_PI / 180.0);
                // Try both directions, pick the one farther from existing atoms
                double a1 = existAngle + bendAngle;
                double a2 = existAngle - bendAngle;
                Point2D p1 = {coords[u].x + std::cos(a1), coords[u].y + std::sin(a1)};
                Point2D p2 = {coords[u].x + std::cos(a2), coords[u].y + std::sin(a2)};
                double minDist1 = 1e10, minDist2 = 1e10;
                for (int w = 0; w < n; w++) {
                    if (!placed[w] || w == u) continue;
                    double d1 = (p1.x - coords[w].x) * (p1.x - coords[w].x)
                              + (p1.y - coords[w].y) * (p1.y - coords[w].y);
                    double d2 = (p2.x - coords[w].x) * (p2.x - coords[w].x)
                              + (p2.y - coords[w].y) * (p2.y - coords[w].y);
                    if (d1 < minDist1) minDist1 = d1;
                    if (d2 < minDist2) minDist2 = d2;
                }
                baseAngle = (minDist1 >= minDist2) ? a1 : a2;
            } else {
                // Multiple placed neighbors: place opposite to centroid of neighbors
                double avgAngle = sumAngle / placedNbCount;
                baseAngle = avgAngle + M_PI;
            }

            coords[v].x = coords[u].x + std::cos(baseAngle);
            coords[v].y = coords[u].y + std::sin(baseAngle);
            placed[v] = true;
            bfsQ.push_back(v);
        }
    }

    // --- Phase 4: Short force-directed refinement (50 iterations) ---
    forceDirectedLayout(g, coords, 50, 1.0);

    // --- Phase 5: Overlap resolution ---
    resolveOverlaps(coords, 0.3, 50);

    // --- Phase 6: Crossing reduction ---
    reduceCrossings(g, coords, 500);

    // --- Phase 7: Canonical orientation ---
    transform::canonicalOrientation(g, coords);

    // --- Phase 8: Bond length normalisation ---
    transform::normaliseBondLength(g, coords, targetBondLength);

    return coords;
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
 * @since 6.11.0
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
    // 1. Double centering: B = -0.5 * H * D^2 * H where H = I - (1/n)*11'
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

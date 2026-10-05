/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only C++17 URF ring finder.
 * Implements:
 *   - Horton candidate generation
 *   - 2-phase GF(2) elimination (Vismara's rule)
 *   - SSSR (Minimum Cycle Basis)
 *   - Relevant Cycle Basis (union of all MCBs)
 *   - Unique Ring Families (Kolodzik et al. 2012)
 *
 * Zero external dependencies -- only requires smsd/mol_graph.hpp.
 */
#pragma once
#ifndef SMSD_RING_FINDER_HPP
#define SMSD_RING_FINDER_HPP

#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// --------------------------------------------------------------------------
// Portable bit intrinsics
// --------------------------------------------------------------------------
namespace smsd { namespace detail_rf {

// MSVC intrinsics header (included at top level for portability)
#if defined(_MSC_VER)
#include <intrin.h>
#endif

inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#elif defined(_MSC_VER)
    return static_cast<int>(__popcnt64(x));
#else
    // Fallback: Hamming weight
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    return static_cast<int>((((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL) * 0x0101010101010101ULL) >> 56);
#endif
}

inline int ctz64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#elif defined(_MSC_VER)
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return static_cast<int>(idx);
#else
    // Fallback: De Bruijn
    static const int debruijn[64] = {
        0,1,2,7,3,13,8,19,4,25,14,28,9,34,20,40,
        5,17,26,38,15,46,29,48,10,31,35,54,21,50,41,57,
        63,6,12,18,24,27,33,39,16,37,45,47,30,53,49,56,
        62,11,23,32,36,44,52,55,61,22,43,51,60,42,59,58
    };
    return debruijn[((x & (~x + 1)) * 0x0218A392CD3D5DBFULL) >> 58];
#endif
}

}} // namespace smsd::detail_rf


namespace smsd {

// ==========================================================================
// BFS shortest path avoiding a specific edge
// ==========================================================================

/**
 * BFS shortest path from @p start to @p target in the graph @p g,
 * excluding the undirected edge (avoidU, avoidV).
 * Returns the atom path [start, ..., target], or empty if unreachable.
 */
inline std::vector<int> bfsPath(const MolGraph& g, int start, int target,
                                int avoidU, int avoidV) {
    if (start == target) return { start };
    int n = g.n;
    std::vector<int> parent(n, -1);
    std::vector<bool> vis(n, false);
    vis[start] = true;
    std::deque<int> q;
    q.push_back(start);
    while (!q.empty()) {
        int u = q.front(); q.pop_front();
        for (int v : g.neighbors[u]) {
            if ((u == avoidU && v == avoidV) || (u == avoidV && v == avoidU))
                continue;
            if (vis[v]) continue;
            vis[v] = true;
            parent[v] = u;
            if (v == target) {
                // Reconstruct path
                std::vector<int> path;
                for (int c = target; c != -1; c = parent[c]) path.push_back(c);
                // parent[start] is -1, so loop ends; but start was never assigned parent
                // since vis[start]=true from the beginning. So we need to add start.
                if (path.back() != start) path.push_back(start);
                std::reverse(path.begin(), path.end());
                return path;
            }
            q.push_back(v);
        }
    }
    return {}; // unreachable
}

// ==========================================================================
// Internal helpers
// ==========================================================================
namespace detail_rf {

// Pack two atom indices into a single int64 for edge lookup
inline int64_t edgeKey(int i, int j) {
    int lo = std::min(i, j), hi = std::max(i, j);
    return (static_cast<int64_t>(lo) << 32) | static_cast<int64_t>(static_cast<uint32_t>(hi));
}

// Find the lowest set bit position in a GF(2) vector of uint64_t words
inline int findPivot(const std::vector<uint64_t>& vec) {
    for (size_t w = 0; w < vec.size(); w++) {
        if (vec[w] != 0) return static_cast<int>((w << 6) + ctz64(vec[w]));
    }
    return -1;
}

// Build edge-incidence vector (GF(2)) for a cycle
inline std::vector<uint64_t> cycleToEdgeVec(
        const std::vector<int>& cycle,
        const std::unordered_map<int64_t, int>& edgeIndex,
        int longWords) {
    std::vector<uint64_t> vec(longWords, 0);
    for (size_t i = 0; i < cycle.size(); i++) {
        int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
        auto it = edgeIndex.find(edgeKey(a, b));
        if (it != edgeIndex.end()) {
            int idx = it->second;
            vec[idx >> 6] ^= 1ULL << (idx & 63);
        }
    }
    return vec;
}

// Canonical edge-set key for cycle deduplication
inline std::vector<int64_t> canonicalCycleKey(const std::vector<int>& cycle) {
    std::vector<int64_t> keys;
    keys.reserve(cycle.size());
    for (size_t i = 0; i < cycle.size(); i++) {
        int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
        int lo = std::min(a, b), hi = std::max(a, b);
        keys.push_back(static_cast<int64_t>(lo) * 100000 + hi);
    }
    std::sort(keys.begin(), keys.end());
    return keys;
}

// Build edge index for the molecule
inline std::unordered_map<int64_t, int> buildEdgeIndex(const MolGraph& g) {
    std::unordered_map<int64_t, int> edgeIndex;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                int64_t key = edgeKey(i, j);
                if (edgeIndex.find(key) == edgeIndex.end()) {
                    int idx = static_cast<int>(edgeIndex.size());
                    edgeIndex[key] = idx;
                }
            }
        }
    }
    return edgeIndex;
}

// Build edge list from graph
inline std::vector<std::pair<int,int>> buildEdgeList(const MolGraph& g) {
    std::vector<std::pair<int,int>> edges;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) edges.push_back({i, j});
        }
    }
    return edges;
}

// Count connected components
inline int countComponents(const MolGraph& g) {
    int n = g.n;
    std::vector<bool> visited(n, false);
    int components = 0;
    for (int i = 0; i < n; i++) {
        if (visited[i]) continue;
        components++;
        std::deque<int> bfsQ;
        bfsQ.push_back(i);
        visited[i] = true;
        while (!bfsQ.empty()) {
            int u = bfsQ.front(); bfsQ.pop_front();
            for (int v : g.neighbors[u]) {
                if (!visited[v]) { visited[v] = true; bfsQ.push_back(v); }
            }
        }
    }
    return components;
}

// Count edges
inline int countEdges(const MolGraph& g) {
    int edgeCount = 0;
    for (int i = 0; i < g.n; i++)
        edgeCount += static_cast<int>(g.neighbors[i].size());
    return edgeCount / 2;
}

// Compute cycle rank = |E| - |V| + components
inline int cycleRank(const MolGraph& g) {
    return countEdges(g) - g.n + countComponents(g);
}

// Union-Find
inline int ufFind(std::vector<int>& parent, int i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
}
inline void ufUnion(std::vector<int>& parent, int a, int b) {
    a = ufFind(parent, a); b = ufFind(parent, b);
    if (a != b) parent[a] = b;
}

} // namespace detail_rf


// ==========================================================================
// Shared Horton + GF(2) pipeline (returns both SSSR and RCB in one pass)
// ==========================================================================

struct RingResult {
    std::vector<std::vector<int>> sssr;  // Minimum Cycle Basis
    std::vector<std::vector<int>> rcb;   // Relevant Cycle Basis
};

inline RingResult runHortonPipeline(const MolGraph& g) {
    RingResult result;
    if (g.n == 0) return result;

    int cRank = detail_rf::cycleRank(g);
    if (cRank <= 0) return result;

    auto edgeIndex = detail_rf::buildEdgeIndex(g);
    auto edgeList  = detail_rf::buildEdgeList(g);
    int totalEdges = static_cast<int>(edgeIndex.size());
    int longWords  = (totalEdges + 63) >> 6;

    // --- Horton candidate generation ---
    std::set<std::vector<int64_t>> seenMasks;
    std::vector<std::vector<uint64_t>> candidateVecs;
    std::vector<std::vector<int>>      candidateAtoms;
    std::vector<int>                   candidateSizes;

    for (int v = 0; v < g.n; v++) {
        for (auto& [eU, eV] : edgeList) {
            auto path1 = bfsPath(g, v, eU, eU, eV);
            if (path1.empty()) continue;
            auto path2 = bfsPath(g, v, eV, eU, eV);
            if (path2.empty()) continue;

            std::unordered_set<int> p1Interior;
            for (size_t i = 1; i < path1.size(); i++) p1Interior.insert(path1[i]);
            bool disjoint = true;
            for (size_t i = 1; i < path2.size() - 1; i++) {
                if (p1Interior.count(path2[i])) { disjoint = false; break; }
            }
            if (!disjoint) continue;

            int ringSize = static_cast<int>(path1.size() + path2.size()) - 1;
            if (ringSize < 3 || ringSize > 20) continue;

            std::vector<int> atoms;
            atoms.reserve(ringSize);
            for (int x : path1) atoms.push_back(x);
            for (int i = static_cast<int>(path2.size()) - 1; i >= 1; i--)
                atoms.push_back(path2[i]);

            std::vector<uint64_t> vec(longWords, 0);
            {
                auto it = edgeIndex.find(detail_rf::edgeKey(eU, eV));
                if (it != edgeIndex.end()) vec[it->second >> 6] ^= 1ULL << (it->second & 63);
            }
            bool valid = true;
            for (size_t i = 0; i + 1 < path1.size() && valid; i++) {
                auto it = edgeIndex.find(detail_rf::edgeKey(path1[i], path1[i+1]));
                if (it == edgeIndex.end()) valid = false;
                else vec[it->second >> 6] ^= 1ULL << (it->second & 63);
            }
            for (size_t i = 0; i + 1 < path2.size() && valid; i++) {
                auto it = edgeIndex.find(detail_rf::edgeKey(path2[i], path2[i+1]));
                if (it == edgeIndex.end()) valid = false;
                else vec[it->second >> 6] ^= 1ULL << (it->second & 63);
            }
            if (!valid) continue;

            auto maskKey = detail_rf::canonicalCycleKey(atoms);
            if (!seenMasks.insert(maskKey).second) continue;

            candidateVecs.push_back(std::move(vec));
            candidateAtoms.push_back(std::move(atoms));
            candidateSizes.push_back(ringSize);
        }
    }

    // Sort candidates by size (smallest first)
    std::vector<int> order(candidateVecs.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [&](int a, int b) { return candidateSizes[a] < candidateSizes[b]; });

    // --- 2-phase GF(2) interleaved pipeline ---
    std::set<std::vector<int64_t>> rcbSeen;
    std::vector<std::vector<uint64_t>> basisMatrix;
    std::vector<int> basisSizes;
    std::vector<int> basisByLeadingBit(totalEdges, -1);

    for (int idx : order) {
        auto& mask = candidateVecs[idx];
        int cSize = candidateSizes[idx];

        // Phase 1: reduce against STRICTLY SHORTER basis vectors only
        std::vector<uint64_t> phase1Mask(mask);
        for (int b = 0; b < totalEdges; b++) {
            if ((phase1Mask[b >> 6] >> (b & 63)) & 1ULL) {
                int basisIdx = basisByLeadingBit[b];
                if (basisIdx != -1 && basisSizes[basisIdx] < cSize) {
                    for (int w = 0; w < longWords; w++)
                        phase1Mask[w] ^= basisMatrix[basisIdx][w];
                }
            }
        }

        bool relevant = false;
        for (auto w : phase1Mask) if (w != 0) { relevant = true; break; }
        if (!relevant) continue;

        // Add to RCB (deduplicated)
        auto rcbKey = detail_rf::canonicalCycleKey(candidateAtoms[idx]);
        if (rcbSeen.insert(rcbKey).second) {
            result.rcb.push_back(candidateAtoms[idx]);
        }

        // Phase 2: reduce against ALL basis vectors for SSSR
        std::vector<uint64_t> phase2Mask(phase1Mask);
        for (int b = 0; b < totalEdges; b++) {
            if ((phase2Mask[b >> 6] >> (b & 63)) & 1ULL) {
                int basisIdx = basisByLeadingBit[b];
                if (basisIdx != -1) {
                    for (int w = 0; w < longWords; w++)
                        phase2Mask[w] ^= basisMatrix[basisIdx][w];
                }
            }
        }

        bool independent = false;
        int newLeadingBit = -1;
        for (int b = 0; b < totalEdges; b++) {
            if ((phase2Mask[b >> 6] >> (b & 63)) & 1ULL) {
                independent = true; newLeadingBit = b; break;
            }
        }

        if (independent && static_cast<int>(result.sssr.size()) < cRank) {
            result.sssr.push_back(candidateAtoms[idx]);
            basisMatrix.push_back(phase2Mask);
            basisSizes.push_back(cSize);
            basisByLeadingBit[newLeadingBit] = static_cast<int>(basisMatrix.size()) - 1;
        }
    }

    return result;
}

// ==========================================================================
// computeRings — SSSR (thin wrapper over shared pipeline)
// ==========================================================================

/**
 * Computes the Smallest Set of Smallest Rings (SSSR / MCB).
 *
 * Algorithm:
 *   1. Horton candidate generation: for every vertex v and every edge (eU,eV),
 *      BFS from v to eU and v to eV avoiding edge (eU,eV), form cycle from
 *      the two vertex-disjoint paths + the edge.
 *   2. Sort candidates by size.
 *   3. 2-phase GF(2) interleaved pipeline:
 *      - Phase 1: XOR against STRICTLY SHORTER basis vectors only.
 *        If non-zero after phase 1, the cycle is "relevant" (goes to RCB).
 *      - Phase 2: XOR against ALL basis vectors.
 *        If non-zero after phase 2, it's linearly independent (goes to SSSR).
 *
 * @param g   The molecular graph.
 * @return    SSSR rings, each ring as an ordered vector of atom indices.
 */
inline std::vector<std::vector<int>> computeRings(const MolGraph& g) {
    return runHortonPipeline(g).sssr;
}

// (Old duplicated computeRings body removed — now uses shared runHortonPipeline)


// ==========================================================================
// computeRelevantCycles — Relevant Cycle Basis (Vismara's 2-phase rule)
// ==========================================================================

/**
 * Computes the set of relevant cycles (union of all minimum cycle bases).
 * A cycle is relevant if it cannot be expressed as the XOR (symmetric
 * difference) of strictly shorter cycles.
 *
 * Uses the same Horton + 2-phase GF(2) pipeline as computeRings, but
 * collects all cycles that pass phase-1 (relevance test) instead of
 * only the linearly independent ones.
 *
 * @param g   The molecular graph.
 * @return    Relevant cycles, each as an ordered vector of atom indices.
 */
inline std::vector<std::vector<int>> computeRelevantCycles(const MolGraph& g) {
    return runHortonPipeline(g).rcb;
}


// ==========================================================================
// computeURFs — Unique Ring Families (Kolodzik et al. 2012)
// ==========================================================================

/**
 * Computes Unique Ring Families. Two relevant cycles belong to the same URF
 * if they share the same orbit signature (lexmin rotation of orbit labels
 * in both directions) and all atoms are vertex-transitive (same orbit).
 *
 * @param g   The molecular graph.
 * @return    URFs, each URF is a vector of cycles (each cycle is an atom vector).
 */
inline std::vector<std::vector<std::vector<int>>> computeURFs(const MolGraph& g) {
    g.ensureCanonical();   // orbit[] is lazy-computed
    auto rc = computeRelevantCycles(g);
    if (rc.empty()) return {};

    // Build edge index + reverse lookup
    auto edgeIndex = detail_rf::buildEdgeIndex(g);
    int numEdges = static_cast<int>(edgeIndex.size());
    int longWords = (numEdges + 63) >> 6;

    std::vector<std::pair<int,int>> edgeByIndex(numEdges);
    for (auto& [key, idx] : edgeIndex) {
        edgeByIndex[idx] = { static_cast<int>(key >> 32),
                             static_cast<int>(key & 0xFFFFFFFF) };
    }

    // Build edge-incidence vectors for each relevant cycle
    int nrc = static_cast<int>(rc.size());
    std::vector<std::vector<uint64_t>> vectors(nrc, std::vector<uint64_t>(longWords, 0));
    for (int ci = 0; ci < nrc; ci++) {
        auto& c = rc[ci];
        for (size_t i = 0; i < c.size(); i++) {
            int a = c[i], b = c[(i + 1) % c.size()];
            auto it = edgeIndex.find(detail_rf::edgeKey(a, b));
            if (it != edgeIndex.end())
                vectors[ci][it->second >> 6] ^= 1ULL << (it->second & 63);
        }
    }

    // Union-Find
    std::vector<int> family(nrc);
    std::iota(family.begin(), family.end(), 0);

    // Generate orbit-based signatures for each cycle
    // Orbit signature = lexmin sequence of orbit[atom] values considering both directions
    auto getCycleOrbitSig = [&](const std::vector<int>& c) -> std::vector<int> {
        int len = static_cast<int>(c.size());
        std::vector<int> best;
        for (int dir = 0; dir <= 1; dir++) {
            for (int start = 0; start < len; start++) {
                std::vector<int> seq(len);
                for (int i = 0; i < len; i++) {
                    int idx = dir == 0 ? (start + i) % len : (start - i + len) % len;
                    seq[i] = g.orbit[c[idx]];
                }
                if (best.empty() || seq < best) best = seq;
            }
        }
        return best;
    };

    std::vector<std::vector<int>> sigs(nrc);
    for (int i = 0; i < nrc; i++) sigs[i] = getCycleOrbitSig(rc[i]);

    // Group cycles into URFs using two criteria (in order):
    //   1. Orbit-signature match: only valid when every atom in both cycles
    //      shares the same orbit value (vertex-transitive ring systems, e.g.
    //      cubane where all 8 atoms have orbit=0). A matching sequence produced
    //      by cycles whose atoms span multiple orbit classes just means the
    //      cycles look the same locally — not that they are interchangeable in
    //      the molecule, so the guard prevents incorrect merges in fused systems
    //      like naphthalene.
    //   2. Kolodzik single-edge exchange: symmetric difference has exactly 2
    //      edges that share a common vertex.
    for (int i = 0; i < nrc; i++) {
        for (int j = i + 1; j < nrc; j++) {
            if (rc[i].size() != rc[j].size()) continue;

            // 1. Orbit-signature match — guarded: apply only when all atoms
            //    in both cycles belong to the same orbit (vertex-transitive).
            if (sigs[i] == sigs[j]) {
                int orb0 = g.orbit[rc[i][0]];
                bool allSame = true;
                for (int a : rc[i]) { if (g.orbit[a] != orb0) { allSame = false; break; } }
                if (allSame) {
                    for (int a : rc[j]) { if (g.orbit[a] != orb0) { allSame = false; break; } }
                }
                if (allSame) {
                    detail_rf::ufUnion(family, i, j);
                    continue;
                }
            }

            // Note: Kolodzik's single-edge-exchange rule (XOR popcount == 2)
            // is mathematically impossible on simple chemical graphs — the
            // symmetric difference of two cycles is always Eulerian, and an
            // Eulerian subgraph with exactly 2 edges requires parallel edges.
            // The orbit-signature check above already handles URF deduplication.
        }
    }

    // Group by family
    std::unordered_map<int, std::vector<int>> groups;
    for (int i = 0; i < nrc; i++) groups[detail_rf::ufFind(family, i)].push_back(i);

    std::vector<std::vector<std::vector<int>>> result;
    for (auto& [rep, members] : groups) {
        std::vector<std::vector<int>> fam;
        for (int m : members) fam.push_back(rc[m]);
        result.push_back(std::move(fam));
    }
    return result;
}

// ==========================================================================
// computeSSSR — explicit SSSR API (alias for computeRings, pre-sorted)
// ==========================================================================

/**
 * Return the Smallest Set of Smallest Rings (minimum cycle basis).
 * Rings are sorted ascending by size. No redundant candidates.
 *
 * This is equivalent to computeRings() but with an explicit SSSR name
 * and a guaranteed ascending size sort.
 *
 * @param g   The molecular graph.
 * @return    SSSR rings, each ring as an ordered vector of atom indices,
 *            sorted ascending by ring size.
 */
inline std::vector<std::vector<int>> computeSSSR(const MolGraph& g) {
    auto rings = computeRings(g);
    // computeRings already returns rings in roughly ascending order from the
    // Horton pipeline, but sort explicitly to guarantee the contract.
    std::sort(rings.begin(), rings.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() < b.size();
        });
    return rings;
}


// ==========================================================================
// layoutSSSR — layout-optimized ring ordering
// ==========================================================================

/**
 * Return SSSR rings optimized for 2D coordinate generation:
 *   - Largest ring system first
 *   - Fused rings ordered by shared-edge adjacency within each system
 *   - Within each system, rings sorted by size (smallest first for placement)
 *
 * A "ring system" is a connected component of the ring adjacency graph
 * (two rings are adjacent if they share at least one edge).
 *
 * @param g   The molecular graph.
 * @return    Layout-ordered SSSR rings.
 */
inline std::vector<std::vector<int>> layoutSSSR(const MolGraph& g) {
    auto rings = computeSSSR(g);
    if (rings.size() <= 1) return rings;

    int nr = static_cast<int>(rings.size());

    // Build edge sets for each ring
    std::vector<std::set<int64_t>> ringEdgeSets(nr);
    for (int i = 0; i < nr; i++) {
        auto& r = rings[i];
        for (size_t j = 0; j < r.size(); j++) {
            int a = r[j], b = r[(j + 1) % r.size()];
            ringEdgeSets[i].insert(detail_rf::edgeKey(a, b));
        }
    }

    // Build ring adjacency graph (shared-edge adjacency)
    std::vector<std::vector<int>> ringAdj(nr);
    for (int i = 0; i < nr; i++) {
        for (int j = i + 1; j < nr; j++) {
            bool shared = false;
            for (auto& e : ringEdgeSets[i]) {
                if (ringEdgeSets[j].count(e)) { shared = true; break; }
            }
            if (shared) {
                ringAdj[i].push_back(j);
                ringAdj[j].push_back(i);
            }
        }
    }

    // Find connected components (ring systems)
    std::vector<int> comp(nr, -1);
    int nComp = 0;
    for (int i = 0; i < nr; i++) {
        if (comp[i] >= 0) continue;
        int c = nComp++;
        comp[i] = c;
        std::deque<int> bfs;
        bfs.push_back(i);
        while (!bfs.empty()) {
            int u = bfs.front(); bfs.pop_front();
            for (int v : ringAdj[u]) {
                if (comp[v] < 0) { comp[v] = c; bfs.push_back(v); }
            }
        }
    }

    // Group rings by component and compute total atom count per system
    std::vector<std::vector<int>> systems(nComp);
    for (int i = 0; i < nr; i++) systems[comp[i]].push_back(i);

    // Count total unique atoms in each system
    std::vector<int> systemSize(nComp, 0);
    for (int c = 0; c < nComp; c++) {
        std::set<int> atoms;
        for (int ri : systems[c])
            for (int a : rings[ri]) atoms.insert(a);
        systemSize[c] = static_cast<int>(atoms.size());
    }

    // Sort systems: largest first
    std::vector<int> sysOrder(nComp);
    std::iota(sysOrder.begin(), sysOrder.end(), 0);
    std::sort(sysOrder.begin(), sysOrder.end(),
        [&](int a, int b) { return systemSize[a] > systemSize[b]; });

    // Within each system, BFS from the smallest ring to order by adjacency
    std::vector<std::vector<int>> result;
    result.reserve(nr);

    for (int ci : sysOrder) {
        auto& sysRings = systems[ci];
        if (sysRings.size() == 1) {
            result.push_back(rings[sysRings[0]]);
            continue;
        }

        // Sort by size to find the smallest ring as seed
        std::sort(sysRings.begin(), sysRings.end(),
            [&](int a, int b) { return rings[a].size() < rings[b].size(); });

        // BFS from smallest ring through ring adjacency
        std::vector<bool> visited(nr, false);
        std::deque<int> bfs;
        bfs.push_back(sysRings[0]);
        visited[sysRings[0]] = true;
        std::vector<int> ordered;
        while (!bfs.empty()) {
            int u = bfs.front(); bfs.pop_front();
            ordered.push_back(u);
            // Sort neighbors by size before enqueuing (smallest first)
            std::vector<int> nbrs;
            for (int v : ringAdj[u]) {
                if (!visited[v] && comp[v] == ci) nbrs.push_back(v);
            }
            std::sort(nbrs.begin(), nbrs.end(),
                [&](int a, int b) { return rings[a].size() < rings[b].size(); });
            for (int v : nbrs) {
                if (!visited[v]) { visited[v] = true; bfs.push_back(v); }
            }
        }

        for (int ri : ordered) result.push_back(rings[ri]);
    }

    return result;
}

} // namespace smsd

#endif // SMSD_RING_FINDER_HPP

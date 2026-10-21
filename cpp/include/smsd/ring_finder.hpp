/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only C++17 URF ring finder.
 * Implements:
 *   - Horton candidate generation
 *   - 2-phase GF(2) elimination (Vismara's rule)
 *   - SSSR (Minimum Cycle Basis)
 *   - Relevant Cycle Basis (union of all MCBs)
 *   - Unique Ring Families (Kolodzik et al. 2012)
 *
 * Zero external dependencies -- only requires smsd/mol_graph.hpp and smsd/bitops.hpp.
 */
#pragma once
#ifndef SMSD_RING_FINDER_HPP
#define SMSD_RING_FINDER_HPP

#include "smsd/bitops.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <deque>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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

// Canonical edge-set key for cycle deduplication
inline std::vector<int64_t> canonicalCycleKey(const std::vector<int>& cycle) {
    std::vector<int64_t> keys;
    keys.reserve(cycle.size());
    for (size_t i = 0; i < cycle.size(); i++)
        keys.push_back(edgeKey(cycle[i], cycle[(i + 1) % cycle.size()]));
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

/**
 * Rings grouped into ring systems: two rings are adjacent when they share a bond, and
 * each system is a connected component of that adjacency, numbered in ring order.
 */
struct RingSystemGroups {
    std::vector<std::vector<int>> adjacency;  // rings sharing a bond with each ring
    std::vector<int> system;                  // system of each ring
    int systemCount = 0;
};

inline RingSystemGroups groupRingsBySharedBond(const std::vector<std::vector<int>>& rings) {
    const int count = static_cast<int>(rings.size());
    std::vector<std::set<int64_t>> edgeSets(count);
    for (int i = 0; i < count; i++) {
        const auto& ring = rings[i];
        for (size_t j = 0; j < ring.size(); j++)
            edgeSets[i].insert(edgeKey(ring[j], ring[(j + 1) % ring.size()]));
    }
    RingSystemGroups groups;
    groups.adjacency.resize(count);
    for (int i = 0; i < count; i++) {
        for (int j = i + 1; j < count; j++) {
            for (const auto& edge : edgeSets[i]) {
                if (edgeSets[j].count(edge)) {
                    groups.adjacency[i].push_back(j);
                    groups.adjacency[j].push_back(i);
                    break;
                }
            }
        }
    }
    groups.system.assign(count, -1);
    for (int i = 0; i < count; i++) {
        if (groups.system[i] >= 0) continue;
        const int system = groups.systemCount++;
        groups.system[i] = system;
        std::vector<int> pending = {i};
        while (!pending.empty()) {
            const int ring = pending.back(); pending.pop_back();
            for (int other : groups.adjacency[ring])
                if (groups.system[other] < 0) { groups.system[other] = system; pending.push_back(other); }
        }
    }
    return groups;
}

} // namespace detail_rf


// ==========================================================================
// Shared Horton + GF(2) pipeline (returns both SSSR and RCB in one pass)
// ==========================================================================

struct RingResult {
    std::vector<std::vector<int>> sssr;  // Minimum Cycle Basis
    std::vector<std::vector<int>> rcb;   // Relevant Cycle Basis
};

/**
 * Work bound of the exhaustive ring candidate search below, counted as the Java implementation
 * counts it, in atoms and neighbours read, path atoms and the bonds of each candidate, so that
 * both implementations give the search up at the same point: it pairs atoms with the bonds of their ring
 * systems through two breadth-first searches, which grows with the cube of the ring system
 * size. Ordinary molecules, including those of several hundred atoms, finish well within it.
 * A search that exhausts it, which takes a fraction of a second, is abandoned; the graph,
 * such as a very large ring system, then takes its candidates from
 * MolGraph::computeRelevantCycles() and MolGraph::computeRings(), whose work is bounded, and
 * the same selection rules apply.
 */
constexpr uint64_t RING_CANDIDATE_WORK_LIMIT = uint64_t(1) << 26;

namespace detail_rf {

/**
 * The two-phase selection of the ring candidate search, over given candidates in ascending
 * size: a candidate independent of strictly shorter selected cycles is relevant, and a
 * relevant one independent of every selected cycle joins the basis until it has cycleRank
 * cycles. Each basis row keeps its lowest bond as pivot, so a reduction reads only the rows
 * its own lowest bonds select.
 */
inline RingResult selectRings(const MolGraph& g, std::vector<std::vector<int>> candidates, int cycleRank) {
    RingResult result;
    const auto edgeIndex = buildEdgeIndex(g);
    const int totalEdges = static_cast<int>(edgeIndex.size());
    const int longWords = (totalEdges + 63) >> 6;
    std::stable_sort(candidates.begin(), candidates.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); });
    std::vector<std::vector<uint64_t>> rows;
    std::vector<size_t> rowSizes;
    std::vector<int> rowByPivot(totalEdges, -1);
    std::vector<uint64_t> bits(longWords);
    // Reduces bits by the rows shorter than size (all rows when size is 0); returns the
    // lowest bond left without a row, or -1 when the bits are in their span.
    auto reduce = [&](size_t size) {
        for (int word = 0; word < longWords; ++word) {
            while (bits[word] != 0) {
                const int pivot = (word << 6) + ctz64(bits[word]);
                const int row = rowByPivot[pivot];
                if (row < 0 || (size != 0 && rowSizes[row] >= size)) return pivot;
                for (int other = word; other < longWords; ++other) bits[other] ^= rows[row][other];
            }
        }
        return -1;
    };
    for (auto& cycle : candidates) {
        std::fill(bits.begin(), bits.end(), 0);
        for (size_t i = 0; i < cycle.size(); i++) {
            const auto it = edgeIndex.find(edgeKey(cycle[i], cycle[(i + 1) % cycle.size()]));
            if (it != edgeIndex.end()) bits[it->second >> 6] ^= 1ULL << (it->second & 63);
        }
        if (reduce(cycle.size()) < 0) continue;
        const int pivot = reduce(0);
        const bool independent = pivot >= 0 && static_cast<int>(result.sssr.size()) < cycleRank;
        if (independent) {
            rowByPivot[pivot] = static_cast<int>(rows.size());
            rows.push_back(bits);
            rowSizes.push_back(cycle.size());
            result.sssr.push_back(cycle);
        }
        result.rcb.push_back(std::move(cycle));
    }
    return result;
}

/**
 * The exhaustive ring perception: Horton candidate cycles from every atom and bond of a
 * ring system, selected into the relevant cycles and a minimum cycle basis by the GF(2)
 * rules above. Beyond RING_CANDIDATE_WORK_LIMIT the candidates come from the graph's own
 * bounded ring perception instead.
 */
inline RingResult perceiveRingSets(const MolGraph& g) {
    RingResult result;
    if (g.n == 0) return result;

    int cRank = detail::cycleRank(g);
    if (cRank <= 0) return result;

    // Relevant cycles first, then any basis ring they omit, each cycle once.
    auto boundedPipeline = [&] {
        std::vector<std::vector<int>> candidates = g.computeRelevantCycles();
        std::set<std::vector<int64_t>> seen;
        for (const auto& cycle : candidates) seen.insert(detail_rf::canonicalCycleKey(cycle));
        for (const auto& cycle : g.computeRings())
            if (seen.insert(detail_rf::canonicalCycleKey(cycle)).second) candidates.push_back(cycle);
        return detail_rf::selectRings(g, std::move(candidates), cRank);
    };
    uint64_t work = 0;

    auto edgeIndex = detail_rf::buildEdgeIndex(g);
    auto edgeList  = detail_rf::buildEdgeList(g);
    int totalEdges = static_cast<int>(edgeIndex.size());
    int longWords  = (totalEdges + 63) >> 6;

    // --- Horton candidate generation ---
    std::set<std::vector<int64_t>> seenMasks;
    std::vector<std::vector<uint64_t>> candidateVecs;
    std::vector<std::vector<int>>      candidateAtoms;
    std::vector<int>                   candidateSizes;

    // Only pairs of an atom and a bond that can add a candidate are searched, in the same
    // order, as in the Java implementation. A bridge leaves one of its atoms unreachable. From outside
    // the ring system of a bond, both paths enter the system through the same atom, so they
    // overlap and close no ring. In a system that is one simple ring, every pair inside it
    // repeats the ring that its lowest atom and first bond give.
    const auto systems = g.cycleSystems();
    const bool bySystem = !systems.ofAtom.empty();
    std::vector<std::vector<int>> systemEdges;
    std::vector<int> lowestAtom;
    if (bySystem) {
        systemEdges.resize(systems.atoms.size());
        lowestAtom.assign(systems.atoms.size(), -1);
        for (size_t edge = 0; edge < edgeList.size(); ++edge) {
            const auto [first, second] = edgeList[edge];
            const int system = systems.ofAtom[first];
            if (system < 0 || systems.ofAtom[second] != system) continue;
            systemEdges[system].push_back(static_cast<int>(edge));
        }
        for (int atom = g.n - 1; atom >= 0; --atom)
            if (systems.ofAtom[atom] >= 0) lowestAtom[systems.ofAtom[atom]] = atom;
    }
    std::vector<int> pairs;
    // bfsPath() with reused buffers. With ring systems known, both targets lie in the system of
    // the bond and the search enters no other system; atoms of that system are reached in the
    // same order and with the same parents, so the path is the one bfsPath() finds.
    std::vector<int> parent(g.n, -1), queue;
    std::vector<uint32_t> visited(g.n, 0);
    uint32_t stamp = 0;
    auto shortestPath = [&](int start, int target, int avoidU, int avoidV) {
        std::vector<int> path;
        if (start == target) { path.push_back(start); return path; }
        if (++stamp == 0) { std::fill(visited.begin(), visited.end(), 0u); stamp = 1; }
        const int system = bySystem ? systems.ofAtom[target] : -1;
        visited[start] = stamp;
        parent[start] = -1;
        queue.assign(1, start);
        for (size_t head = 0; head < queue.size(); ++head) {
            const int u = queue[head];
            work += 1 + g.neighbors[u].size();
            for (int next : g.neighbors[u]) {
                if ((u == avoidU && next == avoidV) || (u == avoidV && next == avoidU)) continue;
                if (visited[next] == stamp || (bySystem && systems.ofAtom[next] != system)) continue;
                visited[next] = stamp;
                parent[next] = u;
                if (next == target) {
                    for (int c = target; c != -1; c = parent[c]) path.push_back(c);
                    std::reverse(path.begin(), path.end());
                    return path;
                }
                queue.push_back(next);
            }
        }
        return path;
    };

    for (int v = 0; v < g.n; v++) {
        pairs.clear();
        if (!bySystem) {
            for (size_t edge = 0; edge < edgeList.size(); ++edge) pairs.push_back(static_cast<int>(edge));
        } else {
            const int system = systems.ofAtom[v];
            if (system >= 0) {
                const auto& edges = systemEdges[system];
                if (!systems.simpleRing[system]) pairs.insert(pairs.end(), edges.begin(), edges.end());
                else if (v == lowestAtom[system] && !edges.empty()) pairs.push_back(edges.front());
            }
        }
        for (int edge : pairs) {
            const int eU = edgeList[edge].first, eV = edgeList[edge].second;
            if (work > RING_CANDIDATE_WORK_LIMIT) return boundedPipeline();
            auto path1 = shortestPath(v, eU, eU, eV);
            if (path1.empty()) continue;
            auto path2 = shortestPath(v, eV, eU, eV);
            if (path2.empty()) continue;
            work += path1.size() + path2.size();

            // The paths share only v. A first path through the far end of the bond walks
            // that end twice and closes no simple ring.
            std::unordered_set<int> p1Interior;
            for (size_t i = 1; i < path1.size(); i++) p1Interior.insert(path1[i]);
            bool disjoint = true;
            for (size_t i = 1; i < path2.size(); i++) {
                if (p1Interior.count(path2[i])) { disjoint = false; break; }
            }
            if (!disjoint) continue;

            int ringSize = static_cast<int>(path1.size() + path2.size()) - 1;
            if (ringSize < 3) continue;

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
            // The bonds of the candidate, read as the Java implementation reads its bond mask.
            work += static_cast<uint64_t>(totalEdges);

            auto maskKey = detail_rf::canonicalCycleKey(atoms);
            if (!seenMasks.insert(maskKey).second) continue;

            candidateVecs.push_back(std::move(vec));
            candidateAtoms.push_back(std::move(atoms));
            candidateSizes.push_back(ringSize);
        }
    }

    // Sort candidates by size (smallest first). Equal sizes keep the order in which they were
    // found, as in the Java implementation, so every standard library selects the same rings.
    std::vector<int> order(candidateVecs.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
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

} // namespace detail_rf

/// @deprecated Use detail_rf::perceiveRingSets(), or the public computeRings() and
/// computeRelevantCycles() functions below.
[[deprecated("Use detail_rf::perceiveRingSets(), computeRings() or computeRelevantCycles().")]]
inline RingResult runHortonPipeline(const MolGraph& g) {
    return detail_rf::perceiveRingSets(g);
}

// ==========================================================================
// computeRings — SSSR (thin wrapper over shared pipeline)
// ==========================================================================

/**
 * Computes the Smallest Set of Smallest Rings (SSSR / MCB).
 *
 * This free form is the canonical ring perception. MolGraph::computeRings() is the
 * graph's bounded shortest-path basis: it finds the same relevant cycles and the same
 * ring sizes, but may choose differently among equal-size rings of a symmetric system.
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
    return detail_rf::perceiveRingSets(g).sssr;
}

namespace detail_rf {

/**
 * A smallest set of smallest rings with as many rings as the cycle rank: the graph's bounded
 * basis, MolGraph::computeRings(), when it is complete, else the rings of the candidate search,
 * which the Java implementation reads. The bounded basis keeps one shortest ring through each bond, so it
 * can omit a ring whose every bond lies in an equally small ring, as the central ring of
 * coronene.
 */
inline std::vector<std::vector<int>> completeBasis(const MolGraph& g) {
    const auto& bounded = g.computeRings();
    if (static_cast<int>(bounded.size()) >= detail::cycleRank(g)) return bounded;
    auto rings = computeRings(g);
    return rings.size() > bounded.size() ? rings : bounded;
}

} // namespace detail_rf

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
    return detail_rf::perceiveRingSets(g).rcb;
}

/**
 * Per-atom ring counts by default: the relevant cycles through each atom, as the Java
 * implementation counts them, found over each atom's neighbours in index order. The relevant cycles do not
 * depend on which minimum cycle basis is chosen, so cages such as adamantane and cubane count
 * every equivalent ring whatever the atom order. Of several equally short paths the candidate
 * search keeps the first it reaches; reading the neighbours in index order makes that choice
 * follow the atom order alone, in both implementations, and not the order in which a parser listed the
 * bonds. Where a ring system has more relevant cycles of one size than the search reaches,
 * as the macrocycle of a cyclodextrin, each of whose glucose units offers two equally short
 * paths, the counts still depend on the atom order.
 */
inline void MolGraph::computeRingCounts() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    bool indexOrder = true;
    for (int atom = 0; atom < n && indexOrder; ++atom)
        indexOrder = std::is_sorted(neighbors[atom].begin(), neighbors[atom].end());
    if (indexOrder) {
        computeRingCounts(detail_rf::perceiveRingSets(*this).rcb);
        return;
    }
    std::vector<std::vector<int>> ordered(neighbors);
    for (int atom = 0; atom < n; ++atom) {
        auto& row = ordered[atom];
        std::sort(row.begin(), row.end());
        // A graph that lists an atom as its own neighbour, or a neighbour twice, has no ring
        // systems and no valid copy; it is read as given.
        if (std::adjacent_find(row.begin(), row.end()) != row.end()
            || std::binary_search(row.begin(), row.end(), atom)) {
            computeRingCounts(detail_rf::perceiveRingSets(*this).rcb);
            return;
        }
    }
    const MolGraph topology = MolGraph::Builder().atomCount(n).atomicNumbers(atomicNum)
        .setNeighbors(std::move(ordered)).build(false);
    computeRingCounts(detail_rf::perceiveRingSets(topology).rcb);
}


// ==========================================================================
// computeURFs — Unique Ring Families (Kolodzik et al. 2012)
// ==========================================================================

/** Compute Unique Ring Families using shared edges and strictly shorter cycle spans. */
inline std::vector<std::vector<std::vector<int>>> computeURFs(const MolGraph& g) {
    return MolGraph::groupUniqueRingFamilies(computeRelevantCycles(g), detail_rf::buildEdgeIndex(g));
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
    // Horton pipeline, but sort explicitly to guarantee the contract; rings of one
    // size keep their order.
    std::stable_sort(rings.begin(), rings.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() < b.size();
        });
    return rings;
}


// ==========================================================================
// layoutSSSR — layout-optimised ring ordering
// ==========================================================================

/**
 * Return SSSR rings optimised for 2D coordinate generation:
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

    // Ring systems: components of the shared-bond adjacency between rings
    const auto groups = detail_rf::groupRingsBySharedBond(rings);
    const auto& ringAdj = groups.adjacency;
    const auto& comp = groups.system;
    const int nComp = groups.systemCount;

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

    // Sort systems: largest first, systems of one size in ring order
    std::vector<int> sysOrder(nComp);
    std::iota(sysOrder.begin(), sysOrder.end(), 0);
    std::stable_sort(sysOrder.begin(), sysOrder.end(),
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
        std::stable_sort(sysRings.begin(), sysRings.end(),
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
            // Sort neighbours by size before enqueuing (smallest first)
            std::vector<int> nbrs;
            for (int v : ringAdj[u]) {
                if (!visited[v] && comp[v] == ci) nbrs.push_back(v);
            }
            std::stable_sort(nbrs.begin(), nbrs.end(),
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

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Test support: the pairwise ring candidate search of 8.1.0, frozen, for differential tests of
 * detail_rf::perceiveRingSets(), as written there but for one section comment. Only the main loop
 * is frozen; it shares the unchanged helpers (selectRings, certifiedMinimumBasis, the bounded
 * relevant enumeration). The replay below runs the same pair searches and returns each pair's
 * work and paths.
 */
#pragma once

#include "smsd/ring_finder.hpp"

#include <algorithm>
#include <cstdint>
#include <set>
#include <unordered_set>
#include <vector>

namespace smsd {
namespace ring_reference {

using detail_rf::certifiedMinimumBasis;

/** detail_rf::perceiveRingSets() as 8.1.0 wrote it. */
inline RingResult perceiveRingSets(const MolGraph& g, bool enumerateRelevant = true,
                                   uint64_t workLimit = RING_CANDIDATE_WORK_LIMIT) {
    RingResult result;
    if (g.n == 0) return result;

    int cRank = detail::cycleRank(g);
    if (cRank <= 0) return result;

    // Relevant cycles first, then any basis ring they omit, each cycle once.
    auto boundedPipeline = [&] {
        bool legacyTruncated = false;
        std::vector<std::vector<int>> candidates = g.computeRelevantCyclesFromBasis(
            g.computeRings(), legacyTruncated, false);
        std::set<std::vector<int64_t>> seen;
        for (const auto& cycle : candidates) seen.insert(detail_rf::canonicalCycleKey(cycle));
        for (const auto& cycle : g.computeRings())
            if (seen.insert(detail_rf::canonicalCycleKey(cycle)).second) candidates.push_back(cycle);
        RingResult bounded = detail_rf::selectRings(g, std::move(candidates), cRank);
        bounded.legacyRcb = bounded.rcb;
        // The work bound stays recorded: fusion counts still come from the bounded candidates.
        bounded.truncated = true;
        bounded.basisTruncated = !certifiedMinimumBasis(g, bounded.sssr, cRank, workLimit);
        if (enumerateRelevant) {
            bounded.rcb = g.computeRelevantCyclesFromBasis(
                bounded.basisTruncated ? g.computeRings() : bounded.sssr, bounded.relevantTruncated);
            bounded.relevantTruncated = bounded.relevantTruncated || bounded.basisTruncated;
        } else {
            bounded.relevantTruncated = bounded.basisTruncated;
        }
        return bounded;
    };
    uint64_t work = 0;

    auto edgeIndex = detail_rf::buildEdgeIndex(g);
    auto edgeList  = detail_rf::buildEdgeList(g);
    int totalEdges = static_cast<int>(edgeIndex.size());
    int longWords  = (totalEdges + 63) >> 6;

    // --- Candidate generation ---
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
    // Shortest paths with reused buffers. With ring systems known, both targets lie in the system of
    // the bond and the search enters no other system; atoms of that system are reached in the
    // same order and with the same parents, so the path is the one bfsPath() finds.
    std::vector<int> parent(g.n, -1), queue;
    std::vector<uint32_t> visited(g.n, 0);
    uint32_t stamp = 0;
    auto shortestPaths = [&](int start, int firstTarget, int secondTarget, int avoidU, int avoidV,
                             std::vector<int>& firstPath, std::vector<int>& secondPath) {
        bool firstReached = start == firstTarget, secondReached = start == secondTarget;
        if (firstReached && secondReached) {
            firstPath.push_back(start);
            secondPath.push_back(start);
            return;
        }
        if (++stamp == 0) { std::fill(visited.begin(), visited.end(), 0u); stamp = 1; }
        const int system = bySystem ? systems.ofAtom[firstTarget] : -1;
        visited[start] = stamp;
        parent[start] = -1;
        queue.assign(1, start);
        for (size_t head = 0; head < queue.size() && !(firstReached && secondReached); ++head) {
            const int u = queue[head];
            work += 1 + g.neighbors[u].size();
            for (int next : g.neighbors[u]) {
                if ((u == avoidU && next == avoidV) || (u == avoidV && next == avoidU)) continue;
                if (visited[next] == stamp || (bySystem && systems.ofAtom[next] != system)) continue;
                visited[next] = stamp;
                parent[next] = u;
                if (next == firstTarget) firstReached = true;
                if (next == secondTarget) secondReached = true;
                if (firstReached && secondReached) break;
                queue.push_back(next);
            }
        }
        if (!firstReached) return;
        for (int c = firstTarget; c != -1; c = parent[c]) firstPath.push_back(c);
        std::reverse(firstPath.begin(), firstPath.end());
        if (!secondReached) return;
        for (int c = secondTarget; c != -1; c = parent[c]) secondPath.push_back(c);
        std::reverse(secondPath.begin(), secondPath.end());
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
            if (work > workLimit) return boundedPipeline();
            std::vector<int> path1, path2;
            shortestPaths(v, eU, eV, eU, eV, path1, path2);
            if (path1.empty()) continue;
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
        result.legacyRcb.push_back(candidateAtoms[idx]);

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

    result.basisTruncated = static_cast<int>(result.sssr.size()) != cRank;
    result.truncated = result.basisTruncated;
    if (enumerateRelevant) {
        result.rcb = g.computeRelevantCyclesFromBasis(result.sssr, result.relevantTruncated);
        result.relevantTruncated = result.relevantTruncated || result.basisTruncated;
    } else {
        result.rcb = result.legacyRcb;
        result.relevantTruncated = result.basisTruncated;
    }
    return result;
}

/** One pair of the pairwise search: its root and bond, its paths and the work it counts. */
struct PairRecord {
    int root = -1, edge = -1;
    uint64_t search = 0;          // atoms and neighbours its breadth-first search reads
    std::vector<int> first, second;
    bool valid = false;           // a candidate: both paths found, sharing only the root
    uint64_t charge = 0;          // everything the pair adds to the work count
};

/** Every pair of the pairwise search in its order, for graphs with ring systems. */
inline std::vector<PairRecord> replayPairs(const MolGraph& g) {
    std::vector<PairRecord> records;
    const auto systems = g.cycleSystems();
    if (systems.ofAtom.empty()) return records;
    const auto edgeList = detail_rf::buildEdgeList(g);
    const uint64_t totalEdges = detail_rf::buildEdgeIndex(g).size();
    std::vector<std::vector<int>> systemEdges(systems.atoms.size());
    std::vector<int> lowestAtom(systems.atoms.size(), -1);
    for (size_t edge = 0; edge < edgeList.size(); ++edge) {
        const int system = systems.ofAtom[edgeList[edge].first];
        if (system < 0 || systems.ofAtom[edgeList[edge].second] != system) continue;
        systemEdges[system].push_back(static_cast<int>(edge));
    }
    for (int atom = g.n - 1; atom >= 0; --atom)
        if (systems.ofAtom[atom] >= 0) lowestAtom[systems.ofAtom[atom]] = atom;
    std::vector<int> parent(g.n, -1), queue;
    std::vector<char> visited(g.n, 0);
    for (int v = 0; v < g.n; ++v) {
        const int system = systems.ofAtom[v];
        if (system < 0) continue;
        std::vector<int> pairs;
        if (!systems.simpleRing[system]) pairs = systemEdges[system];
        else if (v == lowestAtom[system] && !systemEdges[system].empty()) pairs.push_back(systemEdges[system].front());
        for (int edge : pairs) {
            PairRecord record;
            record.root = v;
            record.edge = edge;
            const int eU = edgeList[edge].first, eV = edgeList[edge].second;
            bool firstReached = v == eU, secondReached = v == eV;
            std::fill(visited.begin(), visited.end(), 0);
            visited[v] = 1;
            parent[v] = -1;
            queue.assign(1, v);
            for (size_t head = 0; head < queue.size() && !(firstReached && secondReached); ++head) {
                const int u = queue[head];
                record.search += 1 + g.neighbors[u].size();
                for (int next : g.neighbors[u]) {
                    if ((u == eU && next == eV) || (u == eV && next == eU)) continue;
                    if (visited[next] || systems.ofAtom[next] != system) continue;
                    visited[next] = 1;
                    parent[next] = u;
                    if (next == eU) firstReached = true;
                    if (next == eV) secondReached = true;
                    if (firstReached && secondReached) break;
                    queue.push_back(next);
                }
            }
            record.charge = record.search;
            if (firstReached && secondReached) {
                for (int c = eU; c != -1; c = parent[c]) record.first.push_back(c);
                for (int c = eV; c != -1; c = parent[c]) record.second.push_back(c);
                std::reverse(record.first.begin(), record.first.end());
                std::reverse(record.second.begin(), record.second.end());
                record.charge += record.first.size() + record.second.size();
                std::unordered_set<int> interior(record.first.begin() + 1, record.first.end());
                record.valid = record.first.size() + record.second.size() >= 4;
                for (size_t i = 1; i < record.second.size() && record.valid; ++i)
                    if (interior.count(record.second[i])) record.valid = false;
                if (record.valid) record.charge += totalEdges;
            }
            records.push_back(std::move(record));
        }
    }
    return records;
}

/** The work of all pairs before the last, which decides the fallback, and the last pair's. */
struct WorkTotals { uint64_t pairs = 0, beforeLast = 0, last = 0; };

inline WorkTotals workTotals(const std::vector<PairRecord>& records) {
    WorkTotals totals;
    totals.pairs = records.size();
    for (size_t k = 0; k < records.size(); ++k) {
        if (k + 1 < records.size()) totals.beforeLast += records[k].charge;
        else totals.last = records[k].charge;
    }
    return totals;
}

} // namespace ring_reference
} // namespace smsd

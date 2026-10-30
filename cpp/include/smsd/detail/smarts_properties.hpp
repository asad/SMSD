/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMARTS target properties: valence, hydrogen count, ring membership and size, ring
 * connectivity, heteroatom neighbours, heavy degree and hybridisation, computed once
 * per target in a cache.
 *
 * Part of smarts_parser.hpp.
 */
#pragma once

#include "smsd/detail/smarts_expressions.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <vector>

namespace smsd {
namespace smarts_detail {

// ============================================================================
// Target-side atom property cache
// ============================================================================

/**
 * Implicit hydrogens of a target atom: the stored count when hydrogen counts are known, otherwise
 * the shared SMILES valence model (as MolGraph::attachedHydrogenCount), so N, P, S and the other
 * multivalent elements may take any default valence and aromatic atoms count their pi bond.
 */
inline int computeHCountFromValence(const MolGraph& g, int idx, int explicitValence) {
    if (g.hydrogenCountsKnown && idx < static_cast<int>(g.hydrogenCount.size()))
        return g.hydrogenCount[idx];
    const int charge = idx < static_cast<int>(g.formalCharge.size()) ? g.formalCharge[idx] : 0;
    return ::smsd::detail::computeImplicitH(g.atomicNum[idx], g.aromatic[idx] != 0, explicitValence, charge);
}

/**
 * Total valence for SMARTS v: explicit bond orders, with aromatic bonds counted once, plus
 * hydrogens.  An aromatic atom below its lowest charge-adjusted default valence takes one
 * Kekulé double bond; an occupied lone-pair donor such as pyrrole nitrogen does not.
 */
inline int computeTotalValence(const MolGraph& g, int atom, int explicitValence,
                               int implicitHydrogens) {
    int valence = explicitValence + implicitHydrogens;
    if (!g.aromatic[atom]) return valence;
    const auto valences = ::smsd::detail::defaultValences(g.atomicNum[atom]);
    if (valences.empty()) return valence;
    const int64_t normalValence = ::smsd::detail::chargeAdjustedValence(
        g.atomicNum[atom], valences.front(), g.formalCharge[atom]);
    if (valence < normalValence) ++valence;
    return valence;
}

/**
 * SMARTS R<n>: the rings of the smallest set of smallest rings through the atom, as Daylight
 * defines it, with as many rings as the cycle rank. The graph's fusion counts retain the
 * bounded Horton-candidate perception used in 8.0.x; complete-ring matching uses relevant cycles.
 */
inline int computeRingMembershipRaw(const MolGraph& g, int idx) {
    int rings = 0;
    for (const auto& ring : ::smsd::detail_rf::completeBasis(g))
        rings += static_cast<int>(std::count(ring.begin(), ring.end(), idx));
    return rings;
}

inline int computeSmallestRingSizeRaw(const MolGraph& g, int idx) {
    if (!g.ring[idx]) return 0;

    int best = 999;
    for (int nb : g.neighbors[idx]) {
        if (!g.bondInRing(idx, nb)) continue;

        std::vector<int> dist(g.n, -1);
        dist[nb] = 1;
        std::deque<int> queue;
        queue.push_back(nb);
        bool found = false;

        while (!queue.empty() && !found) {
            int u = queue.front();
            queue.pop_front();
            for (int v : g.neighbors[u]) {
                if (!g.bondInRing(u, v)) continue;
                if (u == nb && v == idx) continue;
                if (v == idx) {
                    int ringSize = dist[u] + 1;
                    if (ringSize < best) best = ringSize;
                    found = true;
                    break;
                }
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
    }

    return best < 999 ? best : 0;
}

inline int computeRingConnectivity(const MolGraph& g, int idx) {
    int rc = 0;
    for (int nb : g.neighbors[idx]) {
        if (g.bondInRing(idx, nb)) ++rc;
    }
    return rc;
}

inline int computeHeteroNeighbors(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        int z = g.atomicNum[nb];
        if (z != 6 && z != 1) ++count;
    }
    return count;
}

inline int computeAliphaticHeteroNeighbors(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        int z = g.atomicNum[nb];
        if (z != 6 && z != 1 && !g.aromatic[nb]) ++count;
    }
    return count;
}

inline int computeHeavyDegree(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        if (g.atomicNum[nb] != 1) ++count;
    }
    return count;
}

inline int computeHybridization(const MolGraph& g, int idx) {
    if (g.aromatic[idx]) return 2; // aromatic => SP2

    int nDouble = 0;
    int nTriple = 0;
    for (int nb : g.neighbors[idx]) {
        int bo = g.bondOrder(idx, nb);
        if (bo == 2) ++nDouble;
        else if (bo == 3) ++nTriple;
    }

    if (nTriple >= 1 || nDouble >= 2) return 1; // SP
    if (nDouble == 1) return 2;                  // SP2
    return 3;                                    // SP3
}

/**
 * Per-thread target properties, filled atom by atom on first use.
 *
 * Neighbour counts are cheap and are evaluated directly; valence, hydrogen
 * count, ring membership and smallest ring size are retained for the current
 * target only. Buffers keep their capacity from one target to the next, so
 * batch matching neither allocates per target nor repeats unused ring searches.
 */
struct SmartsTargetPropertyCache {
    static constexpr uint8_t VALENCE_KNOWN = 1, RING_MEMBERSHIP_KNOWN = 2, RING_SIZE_KNOWN = 4;

    const MolGraph* graph = nullptr;
    uint64_t graphNonce = 0;
    int atomCount = -1;

    std::vector<int> valence;
    std::vector<int> hcount;
    std::vector<int> ringMembership;
    std::vector<int> smallestRingSize;
    std::vector<uint8_t> known;
    // Whether ringMembership holds the SMARTS ring counts of every atom of this target, found
    // under the graph's lock once.
    bool ringCountsChecked = false;
    // Breadth-first buffers; every distance entry is -1 between searches.
    std::vector<int> distance;
    std::vector<int> queue;
    /** The connected component of each atom, numbered once a component-grouped pattern asks. */
    std::vector<int> component;
    bool componentsKnown = false;

    void rebuild(const MolGraph& g) {
        graph = &g;
        graphNonce = g.cacheNonce;
        atomCount = g.n;
        known.assign(g.n, 0);
        ringCountsChecked = false;
        componentsKnown = false;
        valence.resize(g.n);
        hcount.resize(g.n);
        ringMembership.resize(g.n);
        smallestRingSize.resize(g.n);
    }

    void ensureValence(const MolGraph& g, int i) {
        if (known[i] & VALENCE_KNOWN) return;
        const int explicitValence = ::smsd::detail::explicitValence(g, i);
        const int implicitHydrogens = computeHCountFromValence(g, i, explicitValence);
        int hydrogens = implicitHydrogens;
        for (int neighbor : g.neighbors[i])
            if (g.atomicNum[neighbor] == 1) ++hydrogens;
        hcount[i] = hydrogens;
        valence[i] = computeTotalValence(g, i, explicitValence, implicitHydrogens);
        known[i] |= VALENCE_KNOWN;
    }

    int ringMembershipOf(const MolGraph& g, int i) {
        if (!(known[i] & RING_MEMBERSHIP_KNOWN)) {
            if (!ringCountsChecked) {
                // As computeRingMembershipRaw(), for every atom at once: one lock for the lazy
                // ring set, which is read-only once found.
                std::fill(ringMembership.begin(), ringMembership.end(), 0);
                for (const auto& ring : ::smsd::detail_rf::completeBasis(g))
                    for (int atom : ring) ++ringMembership[atom];
                ringCountsChecked = true;
            }
            known[i] |= RING_MEMBERSHIP_KNOWN;
        }
        return ringMembership[i];
    }

    int componentOf(const MolGraph& g, int i) {
        if (!componentsKnown) {
            component.assign(g.n, -1);
            std::vector<int> pending;
            for (int start = 0, count = 0; start < g.n; ++start) {
                if (component[start] >= 0) continue;
                component[start] = count;
                pending.assign(1, start);
                while (!pending.empty()) {
                    const int atom = pending.back();
                    pending.pop_back();
                    for (int neighbor : g.neighbors[atom])
                        if (component[neighbor] < 0) {
                            component[neighbor] = count;
                            pending.push_back(neighbor);
                        }
                }
                ++count;
            }
            componentsKnown = true;
        }
        return component[i];
    }

    /** The same search as computeSmallestRingSizeRaw(), without per-call allocation. */
    int smallestRingSizeOf(const MolGraph& g, int idx) {
        if (known[idx] & RING_SIZE_KNOWN) return smallestRingSize[idx];
        int best = 999;
        if (g.ring[idx]) {
            if (distance.size() < static_cast<size_t>(g.n)) distance.resize(g.n, -1);
            if (queue.size() < static_cast<size_t>(g.n)) queue.resize(g.n);
            for (int nb : g.neighbors[idx]) {
                if (!g.bondInRing(idx, nb)) continue;
                int head = 0, tail = 0;
                distance[nb] = 1;
                queue[tail++] = nb;
                bool found = false;
                while (head < tail && !found) {
                    const int u = queue[head++];
                    for (int v : g.neighbors[u]) {
                        if (!g.bondInRing(u, v)) continue;
                        if (u == nb && v == idx) continue;
                        if (v == idx) {
                            best = std::min(best, distance[u] + 1);
                            found = true;
                            break;
                        }
                        if (distance[v] == -1) {
                            distance[v] = distance[u] + 1;
                            queue[tail++] = v;
                        }
                    }
                }
                for (int visited = 0; visited < tail; ++visited) distance[queue[visited]] = -1;
            }
        }
        smallestRingSize[idx] = best < 999 ? best : 0;
        known[idx] |= RING_SIZE_KNOWN;
        return smallestRingSize[idx];
    }
};

inline SmartsTargetPropertyCache& getSmartsTargetPropertyCache(const MolGraph& g) {
    thread_local SmartsTargetPropertyCache cache;
    if (cache.graph != &g
        || cache.graphNonce != g.cacheNonce
        || cache.atomCount != g.n)
        cache.rebuild(g);
    return cache;
}

// ============================================================================
// Property accessors used by SMARTS matching
// ============================================================================

inline int computeValence(const MolGraph& g, int idx) {
    auto& cache = getSmartsTargetPropertyCache(g);
    cache.ensureValence(g, idx);
    return cache.valence[idx];
}

inline int computeHCount(const MolGraph& g, int idx) {
    auto& cache = getSmartsTargetPropertyCache(g);
    cache.ensureValence(g, idx);
    return cache.hcount[idx];
}

inline int computeRingMembership(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).ringMembershipOf(g, idx);
}

inline int computeSmallestRingSize(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).smallestRingSizeOf(g, idx);
}

} // namespace smarts_detail
} // namespace smsd

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * MCS (Maximum Common Substructure) engine -- header-only C++17 port.
 *
 * Coverage-driven funnel with 7 levels:
 *   L0    Label-frequency upper bound
 *   L0.25 Linear chain fast-path (PEG, polyethylene)
 *   L0.5  Tree fast-path (branched polymers, dendrimers, glycogen)
 *   L0.75 Greedy probe
 *   L1    Substructure containment check
 *   L1.5 Seed-and-extend (bond growth)
 *   L2   McSplit partition refinement (bidirectional, RRSplit maximality)
 *   L3   Bron-Kerbosch + Tomita pivoting + k-core reduction
 *   L4   McGregor extension (bond-grow + atom-frontier DFS + forced assignment)
 *   L5   Extra seeds (only when McSplit and BK disagree)
 *
 * Zero external dependencies -- pure C++17 standard library only.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/vf2pp.hpp"

#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#ifdef _MSC_VER
#include <intrin.h>
#endif
#include <functional>
#include <map>
#include <numeric>
#include <deque>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace smsd {

inline std::map<int,int> canonicalizeMapping(
        const MolGraph& g1, const MolGraph& g2,
        const std::map<int,int>& mapping);

// ---------------------------------------------------------------------------
// MCSOptions
// ---------------------------------------------------------------------------
struct MCSOptions {
    bool     induced         = false;
    bool     connectedOnly   = true;
    bool     disconnectedMCS = false;
    bool     maximizeBonds   = false;
    int      minFragmentSize = 1;
    int      maxFragments    = INT_MAX;
    std::vector<double> atomWeights;
    /// Timeout in milliseconds.  -1 (default) = adaptive:
    /// std::min(30000, 500 + g1.n * g2.n * 2).
    int64_t  timeoutMs       = -1;
    bool     extraSeeds      = true;
    int      seedNeighborhoodRadius = 2;
    int      seedMaxAnchors  = 12;
    bool     useTwoHopNLFInExtension   = true;
    bool     useThreeHopNLFInExtension = false;
    /// Fuzzy template matching: allow up to this many unmatched atoms in
    /// either direction when applying a template via the greedy probe.
    /// 0 = exact match required (default), 2 = allow +/-2 unmatched atoms.
    int      templateFuzzyAtoms        = 0;

    // ---- Reaction-aware post-filter (v6.4.0) ----

    /// Enable built-in reaction-aware post-filtering.  Off by default.
    bool     reactionAware     = false;
    /// Maximum size deficit from the mathematical maximum K to consider.
    /// Default: 2 (consider candidates of size K, K-1, K-2).
    int      nearMcsDelta      = 2;
    /// Maximum number of near-MCS candidates to generate before scoring.
    int      nearMcsCandidates = 20;

    /// Custom post-filter callback.  If set, overrides the built-in
    /// reaction-aware scorer even when reactionAware=true.
    /// Signature: (candidates, g1, g2) -> ranked candidates.
    /// @since 6.6.0
    std::function<std::vector<std::map<int,int>>(
        const std::vector<std::map<int,int>>&,
        const MolGraph&, const MolGraph&)> postFilter;

    /// When true, apply bond-change penalty scoring after reaction-aware
    /// candidate generation. Candidates ranked by chemical plausibility
    /// of implied bond transformations (C-C breaks penalised most).
    /// @since 6.6.0
    bool     bondChangeAware   = false;

    /// Cap the algorithmic funnel at this stage.  Default 5 = all stages.
    ///   0 = L0.25-L0.75 only (chain/tree/greedy, polynomial)
    ///   1 = + L1 substructure + L1.25 augmenting + L1.5 seed-extend
    ///   2 = + L1.75 k-core + L2 McSplit (exponential)
    ///   3 = + L3 Bron-Kerbosch (exponential)
    ///   4 = + L4 McGregor (exponential)
    ///   5 = + L5 extra seeds (default, full funnel)
    /// For reaction mapping, stage 1 is usually sufficient and much faster.
    /// @since 6.12.0
    int      maxStage          = 5;
};

// ---------------------------------------------------------------------------
// MCS Stage Timers — zero-overhead when SMSD_MCS_TIMERS is not defined.
// Enable with -DSMSD_MCS_TIMERS to collect per-stage microsecond timings.
// ---------------------------------------------------------------------------
struct MCSTimers {
    int64_t chainUs      = 0;   // L0.25 linear chain
    int64_t treeUs       = 0;   // L0.5  tree
    int64_t greedyUs     = 0;   // L0.75 greedy probe
    int64_t substructUs  = 0;   // L1    substructure containment
    int64_t seedExtendUs = 0;   // L1.5  seed-and-extend
    int64_t mcSplitUs    = 0;   // L2    McSplit
    int64_t bkUs         = 0;   // L3    Bron-Kerbosch
    int64_t mcGregorUs   = 0;   // L4    McGregor extension
    int64_t totalUs      = 0;   // total wall time
    int     stageReached = 0;   // highest level entered
    bool    timedOut     = false;
};

#ifdef SMSD_MCS_TIMERS
#define SMSD_STAGE_START(name) \
    auto _st_##name = std::chrono::steady_clock::now()
#define SMSD_STAGE_END(name, timers, field) \
    (timers).field = std::chrono::duration_cast<std::chrono::microseconds>( \
        std::chrono::steady_clock::now() - _st_##name).count()
#else
#define SMSD_STAGE_START(name)               ((void)0)
#define SMSD_STAGE_END(name, timers, field)  ((void)0)
#endif

// TimeBudget reused from vf2pp.hpp (TimeBudget)
namespace detail {

// ---------------------------------------------------------------------------
// Scratch -- pre-allocated buffers reused across calls
// ---------------------------------------------------------------------------
struct Scratch {
    std::vector<int>  q2t, t2q;
    std::vector<int>  bestQ2T;
    std::vector<uint8_t> usedQ, usedT;
    std::vector<int>  candBuf, bestCandBuf;
    std::vector<int>  qLabelFreq, tLabelFreq;
    std::vector<int>  jointQ, jointT;
    std::vector<uint8_t> inFrontier;
    std::vector<int>  frontierBuf;
    std::vector<int>  forcedQ, forcedT;
    std::vector<int>  q2tMap;

    void resize(int n1, int n2, int freqSize) {
        q2t.assign(n1, -1);
        t2q.assign(n2, -1);
        bestQ2T.assign(n1, -1);
        usedQ.assign(n1, false);
        usedT.assign(n2, false);
        candBuf.resize(n2);
        bestCandBuf.resize(n2);
        qLabelFreq.assign(freqSize, 0);
        tLabelFreq.assign(freqSize, 0);
        jointQ.resize(n1);
        jointT.resize(n2);
        inFrontier.assign(n1, false);
        frontierBuf.resize(n1);
        forcedQ.resize(n1);
        forcedT.resize(n2);
        q2tMap.assign(n1, -1);
    }
};

// atomsCompatFast() and bondsCompatible() are defined in vf2pp.hpp (detail namespace)
// and reused here via the same include chain.

// ---------------------------------------------------------------------------
// NLF (Neighbour Label Frequency) check
// ---------------------------------------------------------------------------
inline bool nlfOk(const std::vector<int>& fq, const std::vector<int>& ft) {
    int fi = 0, ti = 0;
    int fqsz = static_cast<int>(fq.size());
    int ftsz = static_cast<int>(ft.size());
    while (fi < fqsz) {
        int fLabel = fq[fi], fFreq = fq[fi + 1];
        while (ti < ftsz && ft[ti] < fLabel) ti += 2;
        if (ti >= ftsz || ft[ti] != fLabel || ft[ti + 1] < fFreq) return false;
        fi += 2;
    }
    return true;
}

inline bool nlfCheckOk(int qi, int tj,
                       const std::vector<std::vector<int>>& qNLF1,
                       const std::vector<std::vector<int>>& tNLF1,
                       const std::vector<std::vector<int>>& qNLF2,
                       const std::vector<std::vector<int>>& tNLF2,
                       const std::vector<std::vector<int>>& qNLF3,
                       const std::vector<std::vector<int>>& tNLF3,
                       bool useTwoHop, bool useThreeHop) {
    if (!nlfOk(qNLF1[qi], tNLF1[tj])) return false;
    if (useTwoHop && !qNLF2.empty() && !nlfOk(qNLF2[qi], tNLF2[tj])) return false;
    if (useThreeHop && !qNLF3.empty() && !nlfOk(qNLF3[qi], tNLF3[tj])) return false;
    return true;
}

// Build NLF arrays from a MolGraph
inline std::vector<std::vector<int>> buildNLF(const MolGraph& g,
    std::function<void(const MolGraph&, int, std::map<int,int>&)> collector) {
    std::vector<std::vector<int>> nlf(g.n);
    for (int i = 0; i < g.n; ++i) {
        std::map<int,int> freq;
        collector(g, i, freq);
        nlf[i].reserve(freq.size() * 2);
        for (auto& [lab, cnt] : freq) {
            nlf[i].push_back(lab);
            nlf[i].push_back(cnt);
        }
    }
    return nlf;
}

// NLF label: atomicNum + aromaticity, WITHOUT ring bit — matches MolGraph::nlfLabel().
inline int mcsNlfLabel(const MolGraph& g, int idx) {
    return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
}

inline std::vector<std::vector<int>> buildNLF1(const MolGraph& g) {
    return buildNLF(g, [](const MolGraph& g, int idx, std::map<int,int>& freq) {
        for (int nb : g.neighbors[idx]) freq[mcsNlfLabel(g, nb)]++;
    });
}

inline std::vector<std::vector<int>> buildNLF2(const MolGraph& g) {
    return buildNLF(g, [](const MolGraph& g, int idx, std::map<int,int>& freq) {
        std::vector<bool> direct(g.n, false);
        direct[idx] = true;
        for (int nb : g.neighbors[idx]) direct[nb] = true;
        std::vector<bool> seen(g.n, false);
        for (int nb : g.neighbors[idx])
            for (int j : g.neighbors[nb]) {
                if (direct[j] || seen[j]) continue;
                seen[j] = true;
                freq[mcsNlfLabel(g, j)]++;
            }
    });
}

inline std::vector<std::vector<int>> buildNLF3(const MolGraph& g) {
    return buildNLF(g, [](const MolGraph& g, int idx, std::map<int,int>& freq) {
        std::vector<bool> level1(g.n, false);
        level1[idx] = true;
        for (int nb : g.neighbors[idx]) level1[nb] = true;
        std::vector<bool> level2(g.n, false);
        for (int nb : g.neighbors[idx])
            for (int j : g.neighbors[nb])
                if (!level1[j]) level2[j] = true;
        for (int v = 0; v < g.n; ++v) {
            if (!level2[v]) continue;
            for (int j : g.neighbors[v]) {
                if (!level1[j] && !level2[j]) freq[mcsNlfLabel(g, j)]++;
            }
        }
    });
}

// ---------------------------------------------------------------------------
// Upper bound: label-frequency (Level 0)
// ---------------------------------------------------------------------------
inline int ubLabel(const MolGraph& g, int i, const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    int key = g.atomicNum[i] << 2;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT) key |= g.aromatic[i] ? 2 : 0;
    if (C.ringMatchesRingOnly) key |= g.ring[i] ? 1 : 0;
    return key;
}

inline int ubBaseLabel(const MolGraph& g, int i, const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    int key = g.atomicNum[i] << 1;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT) key |= g.aromatic[i] ? 1 : 0;
    return key;
}

inline int labelFrequencyUpperBound(const MolGraph& g1, const MolGraph& g2,
                                     const ChemOptions& C) {
    if (g1.n == 0 || g2.n == 0) return 0;
    std::unordered_map<int,int> freq1, freq2;
    for (int i = 0; i < g1.n; ++i) freq1[ubLabel(g1, i, C)]++;
    for (int j = 0; j < g2.n; ++j) freq2[ubLabel(g2, j, C)]++;
    int ub = 0;
    for (auto& [lab, cnt] : freq1) {
        auto it = freq2.find(lab);
        if (it != freq2.end()) ub += std::min(cnt, it->second);
    }
    return std::min(ub, std::min(g1.n, g2.n));
}

inline int labelFrequencyUpperBoundDirected(const MolGraph& query, const MolGraph& target,
                                            const ChemOptions& C) {
    if (query.n == 0 || target.n == 0) return 0;
    if (C.ringMatchesRingOnly) return labelFrequencyUpperBound(query, target, C);

    std::unordered_map<int,int> qAny, qRing, tAny, tRing;
    for (int i = 0; i < query.n; ++i) {
        int lab = ubBaseLabel(query, i, C);
        if (C.ringMatchesRingOnly && query.ring[i]) qRing[lab]++;
        else qAny[lab]++;
    }
    for (int j = 0; j < target.n; ++j) {
        int lab = ubBaseLabel(target, j, C);
        tAny[lab]++;
        if (target.ring[j]) tRing[lab]++;
    }

    int ub = 0;
    for (const auto& [lab, cnt] : qRing) {
        auto it = tRing.find(lab);
        if (it != tRing.end()) ub += std::min(cnt, it->second);
    }
    for (const auto& [lab, cnt] : qAny) {
        int avail = tAny[lab];
        if (C.ringMatchesRingOnly) {
            auto qIt = qRing.find(lab);
            auto tIt = tRing.find(lab);
            int usedRing = (qIt != qRing.end() && tIt != tRing.end())
                ? std::min(qIt->second, tIt->second)
                : 0;
            avail -= usedRing;
        }
        if (avail > 0) ub += std::min(cnt, avail);
    }
    return std::min(ub, std::min(query.n, target.n));
}

// ---------------------------------------------------------------------------
// Upper bound: Degree Sequence Bound (DSB) — provably >= LFUB
// ---------------------------------------------------------------------------
// For each atom label present in both graphs, collects degree sequences,
// sorts descending, and greedily counts matchable pairs using two pointers.
// An atom of degree d can only match an atom of degree >= d.
// O(n log n) time.
inline int degreeSequenceUpperBound(const MolGraph& g1, const MolGraph& g2,
                                     const ChemOptions& C) {
    if (g1.n == 0 || g2.n == 0) return 0;

    // Group atom degrees by label
    std::unordered_map<int, std::vector<int>> degs1, degs2;
    for (int i = 0; i < g1.n; ++i)
        degs1[ubLabel(g1, i, C)].push_back(static_cast<int>(g1.neighbors[i].size()));
    for (int j = 0; j < g2.n; ++j)
        degs2[ubLabel(g2, j, C)].push_back(static_cast<int>(g2.neighbors[j].size()));

    int ub = 0;
    for (auto& [lab, d1] : degs1) {
        auto it = degs2.find(lab);
        if (it == degs2.end()) continue;
        auto& d2 = it->second;
        // Sort both descending
        std::sort(d1.begin(), d1.end(), std::greater<int>());
        std::sort(d2.begin(), d2.end(), std::greater<int>());
        // shorter iterates, longer provides candidates
        auto* shorter = &d1;
        auto* longer  = &d2;
        if (d1.size() > d2.size()) std::swap(shorter, longer);
        // Greedy two-pointer match
        int j = 0;
        int longerSz = static_cast<int>(longer->size());
        for (int i = 0; i < static_cast<int>(shorter->size()) && j < longerSz; ++i) {
            if ((*longer)[j] >= (*shorter)[i]) {
                ++ub;
                ++j;
            }
        }
    }
    return std::min(ub, std::min(g1.n, g2.n));
}

// ---------------------------------------------------------------------------
// Mapped-bond counter
// ---------------------------------------------------------------------------
inline int countMappedBonds(const MolGraph& g1, const std::map<int,int>& m) {
    int count = 0;
    for (auto& [qi, ti] : m)
        for (int qk : g1.neighbors[qi])
            if (qk > qi && m.count(qk)) ++count;
    return count;
}

// ---------------------------------------------------------------------------
// MCS scoring
// ---------------------------------------------------------------------------
inline int mcsScore(const MolGraph& g1, const std::map<int,int>& m,
                    const MCSOptions& M) {
    if (!M.atomWeights.empty()) {
        double w = 0.0;
        int awSz = static_cast<int>(M.atomWeights.size());
        for (auto& [qi, ti] : m) {
            if (qi >= 0 && qi < awSz) w += M.atomWeights[qi];
        }
        return static_cast<int>(w * 1000);
    }
    return M.maximizeBonds ? countMappedBonds(g1, m) : static_cast<int>(m.size());
}

// ---------------------------------------------------------------------------
// Phase 2.4: Flat-array pipeline helpers
// ---------------------------------------------------------------------------
// These operate on q2t[] flat arrays (where -1 = unmapped) to avoid
// constructing std::map<int,int> in intermediate pipeline stages.

/// Count mapped bonds using flat q2t array.
inline int countMappedBondsFlat(const MolGraph& g1, const int* q2t, int n1) {
    int count = 0;
    for (int qi = 0; qi < n1; ++qi) {
        if (q2t[qi] < 0) continue;
        for (int qk : g1.neighbors[qi])
            if (qk > qi && q2t[qk] >= 0) ++count;
    }
    return count;
}

/// MCS scoring using flat q2t array.
inline int mcsScoreFlat(const MolGraph& g1, const int* q2t, int n1, int mapSize,
                        const MCSOptions& M) {
    if (!M.atomWeights.empty()) {
        double w = 0.0;
        int awSz = static_cast<int>(M.atomWeights.size());
        for (int qi = 0; qi < n1; ++qi) {
            if (q2t[qi] >= 0 && qi < awSz) w += M.atomWeights[qi];
        }
        return static_cast<int>(w * 1000);
    }
    return M.maximizeBonds ? countMappedBondsFlat(g1, q2t, n1) : mapSize;
}

/// Heteroatom score using flat q2t array.
inline int heteroatomScoreFlat(const MolGraph& g, const int* q2t, int n1) {
    int score = 0;
    for (int i = 0; i < n1; ++i) {
        if (q2t[i] < 0) continue;
        int z = g.atomicNum[i];
        if (z != 6 && z != 1) score++;
    }
    return score;
}

/// Check if candidate flat array is better MCS than current best flat array.
inline bool isBetterMCSFlat(const MolGraph& g1,
                            const int* candQ2T, int candSize, int n1,
                            const int* bestQ2T, int bestSize, int bestHetero) {
    if (candSize > bestSize) return true;
    if (candSize == bestSize && candSize > 0) {
        int candHetero = heteroatomScoreFlat(g1, candQ2T, n1);
        return candHetero > bestHetero;
    }
    return false;
}

/// Materialize a flat q2t array into std::map<int,int>.
inline std::map<int,int> flatToMap(const int* q2t, int n1) {
    std::map<int,int> result;
    for (int i = 0; i < n1; ++i)
        if (q2t[i] >= 0) result[i] = q2t[i];
    return result;
}

/// Flatten a std::map<int,int> into a q2t array (must be pre-allocated, size n1, filled with -1).
inline void mapToFlat(const std::map<int,int>& m, int* q2t, int n1) {
    std::memset(q2t, -1, n1 * sizeof(int));
    for (const auto& [k, v] : m) {
        if (k >= 0 && k < n1) q2t[k] = v;
    }
}

/// Count non-(-1) entries in a flat q2t array.
inline int flatMapSize(const int* q2t, int n1) {
    int count = 0;
    for (int i = 0; i < n1; ++i) if (q2t[i] >= 0) count++;
    return count;
}

// ---------------------------------------------------------------------------
// Largest connected component (in query graph induced by mapped atoms)
// ---------------------------------------------------------------------------
inline std::map<int,int> largestConnected(const MolGraph& g1,
                                           const std::map<int,int>& m,
                                           const MolGraph* g2 = nullptr) {
    if (m.empty()) return m;
    // BFS to find components.
    // Connectivity is defined through COMMON edges: bonds that exist in g1 between
    // two mapped atoms AND (if g2 is provided) also in g2 between their target atoms.
    // This is required for non-induced MCS — a bond in g1 that has no counterpart in
    // g2 is not part of the common subgraph and must not be used for connectivity.
    std::unordered_map<int,int> q2t;
    if (g2) for (auto& [qi, ti] : m) q2t[qi] = ti;
    std::unordered_set<int> mapped;
    for (auto& [qi, ti] : m) mapped.insert(qi);

    std::unordered_set<int> seen;
    std::vector<std::vector<int>> comps;
    for (auto& [qi, ti] : m) {
        if (seen.count(qi)) continue;
        std::vector<int> comp;
        std::vector<int> dq = {qi};
        seen.insert(qi);
        while (!dq.empty()) {
            int u = dq.back(); dq.pop_back();
            comp.push_back(u);
            for (int v : g1.neighbors[u]) {
                if (mapped.count(v) && !seen.count(v)) {
                    if (g1.hasBond(u, v)) {
                        // For non-induced MCS: only use bond u-v as a connectivity
                        // edge if the corresponding bond also exists in g2.
                        if (g2 && !g2->hasBond(q2t.at(u), q2t.at(v))) continue;
                        seen.insert(v);
                        dq.push_back(v);
                    }
                }
            }
        }
        comps.push_back(std::move(comp));
    }
    // Pick largest
    int bestIdx = 0;
    for (int i = 1; i < static_cast<int>(comps.size()); ++i)
        if (comps[i].size() > comps[bestIdx].size()) bestIdx = i;

    if (comps[bestIdx].size() == m.size()) return m;
    std::map<int,int> result;
    for (int qi : comps[bestIdx]) result[qi] = m.at(qi);
    return result;
}

// ---------------------------------------------------------------------------
// Fragment constraints (dMCS)
// ---------------------------------------------------------------------------
inline std::map<int,int> applyFragmentConstraints(
    const MolGraph& g1, const std::map<int,int>& m,
    int minFragSize, int maxFrags) {
    if (m.empty() || (minFragSize <= 1 && maxFrags >= static_cast<int>(m.size()))) return m;

    std::unordered_set<int> mapped;
    for (auto& [qi, ti] : m) mapped.insert(qi);
    std::unordered_set<int> seen;
    std::vector<std::vector<int>> frags;
    for (auto& [qi, ti] : m) {
        if (seen.count(qi)) continue;
        std::vector<int> frag;
        std::vector<int> dq = {qi};
        seen.insert(qi);
        while (!dq.empty()) {
            int u = dq.back(); dq.pop_back();
            frag.push_back(u);
            for (int v : g1.neighbors[u]) {
                if (mapped.count(v) && !seen.count(v) && g1.hasBond(u, v)) {
                    seen.insert(v); dq.push_back(v);
                }
            }
        }
        frags.push_back(std::move(frag));
    }
    // Remove small fragments
    frags.erase(std::remove_if(frags.begin(), frags.end(),
        [minFragSize](const std::vector<int>& f) { return static_cast<int>(f.size()) < minFragSize; }),
        frags.end());
    // Sort descending by size
    std::sort(frags.begin(), frags.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a.size() > b.size(); });
    if (static_cast<int>(frags.size()) > maxFrags)
        frags.resize(maxFrags);

    std::map<int,int> result;
    for (auto& frag : frags)
        for (int qi : frag) result[qi] = m.at(qi);
    return result;
}

// ---------------------------------------------------------------------------
// Ring anchor guard
// ---------------------------------------------------------------------------
inline std::map<int,int> applyRingAnchorGuard(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& m, const ChemOptions& C) {
    if (!C.ringMatchesRingOnly || m.empty()) return m;
    bool qHasRing = false, tHasRing = false;
    for (int i = 0; i < g1.n; ++i) if (g1.ring[i]) { qHasRing = true; break; }
    for (int i = 0; i < g2.n; ++i) if (g2.ring[i]) { tHasRing = true; break; }
    if (qHasRing && !tHasRing) {
        int mappedRing = 0;
        for (auto& [qi, ti] : m) if (g1.ring[qi]) mappedRing++;
        if (mappedRing == 0) return {};
    }
    return m;
}

// ---------------------------------------------------------------------------
// Enforce complete rings
// ---------------------------------------------------------------------------
inline std::map<int,int> enforceCompleteRings(
    const MolGraph& g1, const MolGraph& g2, const std::map<int,int>& m) {
    if (m.empty()) return m;
    auto rings = const_cast<MolGraph&>(g1).computeRings();
    if (rings.empty()) return m;
    std::map<int,int> result = m;
    bool changed = true;
    while (changed) {
        changed = false;
        for (auto& ring : rings) {
            bool anyMapped = false, allMapped = true;
            for (int atom : ring) {
                if (result.count(atom)) anyMapped = true;
                else allMapped = false;
            }
            if (anyMapped && !allMapped) {
                for (int atom : ring) {
                    if (result.erase(atom)) changed = true;
                }
            }
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// Prune to induced subgraph
// ---------------------------------------------------------------------------
inline std::map<int,int> pruneToInduced(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& m, const ChemOptions& C) {
    if (m.empty()) return m;
    std::map<int,int> M = m;
    bool changed;
    do {
        changed = false;
        std::vector<int> keys;
        for (auto& [k,v] : M) keys.push_back(k);
        for (int i = 0; i < static_cast<int>(keys.size()) && !changed; ++i) {
            for (int j = i + 1; j < static_cast<int>(keys.size()) && !changed; ++j) {
                int qi = keys[i], qj = keys[j];
                int ti = M[qi], tj = M[qj];
                int qOrd = g1.bondOrder(qi, qj), tOrd = g2.bondOrder(ti, tj);
                bool qHas = qOrd != 0, tHas = tOrd != 0;
                bool ok = (qHas == tHas);
                if (ok && qHas) ok = bondsCompatible(g1, qi, qj, g2, ti, tj, C);
                if (!ok) {
                    if (g1.degree[qi] >= g1.degree[qj]) M.erase(qi); else M.erase(qj);
                    changed = true;
                }
            }
        }
    } while (changed);
    return M;
}

// ---------------------------------------------------------------------------
// Post-process MCS (ppx)
// ---------------------------------------------------------------------------
inline std::map<int,int> ppx(const MolGraph& g1, const MolGraph& g2,
                              std::map<int,int> ext,
                              const ChemOptions& C, const MCSOptions& M) {
    // Iteratively apply filters until stable (filters can interact)
    bool changed = true;
    while (changed) {
        int startSize = (int)ext.size();
        if (M.induced) ext = pruneToInduced(g1, g2, ext, C);
        if (C.completeRingsOnly) ext = enforceCompleteRings(g1, g2, ext);
        if (!M.disconnectedMCS && M.connectedOnly) ext = largestConnected(g1, ext, &g2);
        changed = (int)ext.size() < startSize;
    }
    ext = applyRingAnchorGuard(g1, g2, ext, C);
    if (M.disconnectedMCS && (M.minFragmentSize > 1 || M.maxFragments < INT_MAX))
        ext = applyFragmentConstraints(g1, ext, M.minFragmentSize, M.maxFragments);
    return ext;
}

// ---------------------------------------------------------------------------
// SmallExactMCSExplorer -- exact branch-and-bound for small molecule pairs
// ---------------------------------------------------------------------------
class SmallExactMCSExplorer {
    using CanonKey = std::vector<std::pair<int,int>>;
    struct ExactCandidate {
        std::map<int,int> mapping;
        int bondCount = 0;
    };

    const MolGraph& g1_;
    const MolGraph& g2_;
    const ChemOptions& C_;
    bool induced_;
    TimeBudget& tb_;
    int upperBound_;
    std::vector<int> q2t_;
    std::vector<int> t2q_;
    std::vector<int> bestQ2T_;
    std::vector<uint8_t> state_;
    std::vector<std::vector<int>> compatTargets_;
    std::map<CanonKey, ExactCandidate> allBest_;
    int bestSize_ = 0;
    int maxResults_ = 1;
    int collectLimit_ = 1;
    bool collectAll_ = false;

    bool consistent(int qi, int tj) const {
        if (t2q_[tj] >= 0) return false;
        for (int qk = 0; qk < g1_.n; ++qk) {
            if (qk == qi || q2t_[qk] < 0) continue;
            int tk = q2t_[qk];
            bool qBond = g1_.bondOrder(qi, qk) != 0;
            bool tBond = g2_.bondOrder(tj, tk) != 0;
            if (qBond) {
                if (!tBond) return false;
                if (!bondsCompatible(g1_, qi, qk, g2_, tj, tk, C_)) return false;
            } else if (induced_ && tBond) {
                return false;
            }
        }
        if (C_.useChirality && !tetraParityCompatible(g1_, qi, g2_, tj, q2t_)) return false;
        return true;
    }

    int mappedNeighborCount(int qi) const {
        int count = 0;
        for (int nb : g1_.neighbors[qi]) {
            if (q2t_[nb] >= 0) ++count;
        }
        return count;
    }

    std::map<int,int> materialize(const std::vector<int>& q2t) const {
        std::map<int,int> mapping;
        for (int qi = 0; qi < g1_.n; ++qi) {
            if (q2t[qi] >= 0) mapping[qi] = q2t[qi];
        }
        return mapping;
    }

    void recordCurrent(int mappedCount) {
        if (mappedCount < bestSize_) return;
        auto mapping = materialize(q2t_);
        if (mappedCount > bestSize_) {
            bestQ2T_ = q2t_;
            bestSize_ = mappedCount;
            allBest_.clear();
        }
        if (!collectAll_) return;
        if (mappedCount != bestSize_) return;
        int bondCount = countMappedBonds(g1_, mapping);
        auto canonical = canonicalizeMapping(g1_, g2_, mapping);
        CanonKey key(canonical.begin(), canonical.end());
        auto it = allBest_.find(key);
        if (it == allBest_.end() || bondCount > it->second.bondCount) {
            if (it == allBest_.end() && static_cast<int>(allBest_.size()) >= collectLimit_) return;
            allBest_[std::move(key)] = ExactCandidate{std::move(mapping), bondCount};
        }
    }

    void collectCandidates(int qi, std::vector<int>& out) const {
        out.clear();
        out.reserve(compatTargets_[qi].size());
        for (int tj : compatTargets_[qi]) {
            if (consistent(qi, tj)) {
                out.push_back(tj);
            }
        }
    }

    int selectNextAtom(std::vector<int>& candidates) const {
        int bestQi = -1;
        int bestMappedNbrs = -1;
        int bestCandCount = INT_MAX;
        int bestDegree = -1;
        std::vector<int> scratch;
        for (int qi = 0; qi < g1_.n; ++qi) {
            if (state_[qi] != 0) continue;
            int mappedNbrs = mappedNeighborCount(qi);
            collectCandidates(qi, scratch);
            int candCount = static_cast<int>(scratch.size());
            int degree = g1_.degree[qi];
            bool better = false;
            if (mappedNbrs != bestMappedNbrs) better = mappedNbrs > bestMappedNbrs;
            else if (g1_.ring[qi] != (bestQi >= 0 ? g1_.ring[bestQi] : 0)) {
                better = g1_.ring[qi] > (bestQi >= 0 ? g1_.ring[bestQi] : 0);
            } else if (degree != bestDegree) better = degree > bestDegree;
            else if (candCount != bestCandCount) better = candCount < bestCandCount;
            else if (bestQi < 0 || qi < bestQi) better = true;

            if (better) {
                bestQi = qi;
                bestMappedNbrs = mappedNbrs;
                bestCandCount = candCount;
                bestDegree = degree;
                candidates = scratch;
            }
        }
        return bestQi;
    }

    void search(int mappedCount, int pendingCount) {
        if (tb_.expired()) return;
        recordCurrent(mappedCount);
        if (!collectAll_ && bestSize_ >= upperBound_) return;
        if (collectAll_ && bestSize_ >= upperBound_
            && static_cast<int>(allBest_.size()) >= collectLimit_) return;
        if (pendingCount <= 0) return;

        int remainingTarget = g2_.n - mappedCount;
        int optimistic = mappedCount + std::min(pendingCount, remainingTarget);
        if ((!collectAll_ && optimistic <= bestSize_)
            || (collectAll_ && optimistic < bestSize_)) return;

        std::vector<int> candidates;
        int qi = selectNextAtom(candidates);
        if (qi < 0) return;

        std::vector<std::pair<int, int>> scored;
        scored.reserve(candidates.size());
        int mappedNbrs = mappedNeighborCount(qi);
        for (int tj : candidates) {
            int score = mappedNbrs * 1000;
            score += (g1_.ring[qi] && g2_.ring[tj]) ? 100 : 0;
            score += 10 - std::min(9, std::abs(g1_.degree[qi] - g2_.degree[tj]));
            scored.push_back({-score, tj});
        }
        std::sort(scored.begin(), scored.end());

        state_[qi] = 1;
        for (const auto& candidate : scored) {
            int tj = candidate.second;
            q2t_[qi] = tj;
            t2q_[tj] = qi;
            search(mappedCount + 1, pendingCount - 1);
            q2t_[qi] = -1;
            t2q_[tj] = -1;
            if ((!collectAll_ && bestSize_ >= upperBound_) || tb_.expired()) {
                state_[qi] = 0;
                return;
            }
        }

        state_[qi] = 2;
        search(mappedCount, pendingCount - 1);
        state_[qi] = 0;
    }

public:
    SmallExactMCSExplorer(const MolGraph& g1, const MolGraph& g2,
                          const ChemOptions& C, bool induced,
                          TimeBudget& tb, int upperBound,
                          const std::map<int,int>& incumbent)
        : g1_(g1),
          g2_(g2),
          C_(C),
          induced_(induced),
          tb_(tb),
          upperBound_(upperBound),
          q2t_(g1.n, -1),
          t2q_(g2.n, -1),
          bestQ2T_(g1.n, -1),
          state_(g1.n, 0),
          compatTargets_(g1.n) {
        for (const auto& entry : incumbent) {
            if (entry.first < 0 || entry.first >= g1.n) continue;
            if (entry.second < 0 || entry.second >= g2.n) continue;
            bestQ2T_[entry.first] = entry.second;
        }
        bestSize_ = static_cast<int>(incumbent.size());
        for (int qi = 0; qi < g1.n; ++qi) {
            for (int tj = 0; tj < g2.n; ++tj) {
                if (atomsCompatFast(g1_, qi, g2_, tj, C_)) {
                    compatTargets_[qi].push_back(tj);
                }
            }
        }
    }

    std::map<int,int> run() {
        collectAll_ = false;
        maxResults_ = 1;
        collectLimit_ = 1;
        allBest_.clear();
        search(0, g1_.n);
        return materialize(bestQ2T_);
    }

    std::vector<std::map<int,int>> runAll(int maxResults) {
        collectAll_ = true;
        maxResults_ = std::max(1, maxResults);
        collectLimit_ = std::min(4096, std::max(maxResults_ * 32, maxResults_));
        allBest_.clear();
        if (bestSize_ > 0) {
            auto incumbent = materialize(bestQ2T_);
            int incumbentBondCount = countMappedBonds(g1_, incumbent);
            auto canonical = canonicalizeMapping(g1_, g2_, incumbent);
            CanonKey key(canonical.begin(), canonical.end());
            allBest_[std::move(key)] = ExactCandidate{
                std::move(incumbent),
                incumbentBondCount
            };
        }
        search(0, g1_.n);
        std::vector<std::pair<CanonKey, ExactCandidate>> ranked;
        ranked.reserve(allBest_.size());
        for (auto& entry : allBest_) {
            ranked.push_back({entry.first, std::move(entry.second)});
        }
        std::stable_sort(ranked.begin(), ranked.end(),
                         [](const auto& lhs, const auto& rhs) {
            if (lhs.second.mapping.size() != rhs.second.mapping.size()) {
                return lhs.second.mapping.size() > rhs.second.mapping.size();
            }
            if (lhs.second.bondCount != rhs.second.bondCount) {
                return lhs.second.bondCount > rhs.second.bondCount;
            }
            return lhs.first < rhs.first;
        });
        std::vector<std::map<int,int>> result;
        result.reserve(std::min<int>(maxResults_, static_cast<int>(ranked.size())));
        for (auto& entry : ranked) {
            if (static_cast<int>(result.size()) >= maxResults_) break;
            result.push_back(std::move(entry.second.mapping));
        }
        if (result.empty() && bestSize_ > 0) {
            result.push_back(materialize(bestQ2T_));
        }
        return result;
    }
};

// ---------------------------------------------------------------------------
// FixedSizeBondMaximizer -- maximize bond count for a fixed atom-count MCS
// ---------------------------------------------------------------------------
class FixedSizeBondMaximizer {
    const MolGraph& g1_;
    const MolGraph& g2_;
    const ChemOptions& C_;
    bool induced_;
    TimeBudget& tb_;
    int requiredSize_;
    std::vector<int> q2t_;
    std::vector<int> t2q_;
    std::vector<int> bestQ2T_;
    std::vector<uint8_t> state_;
    std::vector<std::vector<int>> compatTargets_;
    int bestBondCount_ = -1;

    bool consistent(int qi, int tj) const {
        if (t2q_[tj] >= 0) return false;
        for (int qk = 0; qk < g1_.n; ++qk) {
            if (qk == qi || q2t_[qk] < 0) continue;
            int tk = q2t_[qk];
            bool qBond = g1_.bondOrder(qi, qk) != 0;
            bool tBond = g2_.bondOrder(tj, tk) != 0;
            if (qBond) {
                if (!tBond) return false;
                if (!bondsCompatible(g1_, qi, qk, g2_, tj, tk, C_)) return false;
            } else if (induced_ && tBond) {
                return false;
            }
        }
        if (C_.useChirality && !tetraParityCompatible(g1_, qi, g2_, tj, q2t_)) return false;
        return true;
    }

    int mappedNeighborCount(int qi) const {
        int count = 0;
        for (int nb : g1_.neighbors[qi]) {
            if (q2t_[nb] >= 0) ++count;
        }
        return count;
    }

    int bondGain(int qi, int tj) const {
        int gain = 0;
        for (int qk : g1_.neighbors[qi]) {
            if (q2t_[qk] < 0) continue;
            int tk = q2t_[qk];
            if (g1_.bondOrder(qi, qk) != 0 && g2_.bondOrder(tj, tk) != 0) {
                ++gain;
            }
        }
        return gain;
    }

    std::map<int,int> materialize(const std::vector<int>& q2t) const {
        std::map<int,int> mapping;
        for (int qi = 0; qi < g1_.n; ++qi) {
            if (q2t[qi] >= 0) mapping[qi] = q2t[qi];
        }
        return mapping;
    }

    void collectCandidates(int qi, std::vector<int>& out) const {
        out.clear();
        out.reserve(compatTargets_[qi].size());
        for (int tj : compatTargets_[qi]) {
            if (consistent(qi, tj)) out.push_back(tj);
        }
    }

    int selectNextAtom(std::vector<int>& candidates) const {
        int bestQi = -1;
        int bestMappedNbrs = -1;
        int bestDegree = -1;
        int bestCandCount = INT_MAX;
        std::vector<int> scratch;
        for (int qi = 0; qi < g1_.n; ++qi) {
            if (state_[qi] != 0) continue;
            int mappedNbrs = mappedNeighborCount(qi);
            collectCandidates(qi, scratch);
            int candCount = static_cast<int>(scratch.size());
            int degree = g1_.degree[qi];
            bool better = false;
            if (mappedNbrs != bestMappedNbrs) better = mappedNbrs > bestMappedNbrs;
            else if (g1_.ring[qi] != (bestQi >= 0 ? g1_.ring[bestQi] : 0)) {
                better = g1_.ring[qi] > (bestQi >= 0 ? g1_.ring[bestQi] : 0);
            } else if (degree != bestDegree) better = degree > bestDegree;
            else if (candCount != bestCandCount) better = candCount < bestCandCount;
            else if (bestQi < 0 || qi < bestQi) better = true;
            if (better) {
                bestQi = qi;
                bestMappedNbrs = mappedNbrs;
                bestDegree = degree;
                bestCandCount = candCount;
                candidates = scratch;
            }
        }
        return bestQi;
    }

    void search(int mappedCount, int pendingCount, int mappedBondCount) {
        if (tb_.expiredNow()) return;
        int remainingTarget = g2_.n - mappedCount;
        if (mappedCount + std::min(pendingCount, remainingTarget) < requiredSize_) return;
        if (mappedCount == requiredSize_) {
            if (mappedBondCount > bestBondCount_) {
                bestBondCount_ = mappedBondCount;
                bestQ2T_ = q2t_;
            }
            return;
        }
        if (pendingCount <= 0 || remainingTarget <= 0) return;

        std::vector<int> candidates;
        int qi = selectNextAtom(candidates);
        if (qi < 0) return;

        std::vector<std::pair<int,int>> scored;
        scored.reserve(candidates.size());
        int mappedNbrs = mappedNeighborCount(qi);
        for (int tj : candidates) {
            int score = mappedNbrs * 1000;
            score += bondGain(qi, tj) * 200;
            score += (g1_.ring[qi] && g2_.ring[tj]) ? 100 : 0;
            score += 10 - std::min(9, std::abs(g1_.degree[qi] - g2_.degree[tj]));
            scored.push_back({-score, tj});
        }
        std::sort(scored.begin(), scored.end());

        state_[qi] = 1;
        for (const auto& candidate : scored) {
            int tj = candidate.second;
            q2t_[qi] = tj;
            t2q_[tj] = qi;
            search(mappedCount + 1, pendingCount - 1, mappedBondCount + bondGain(qi, tj));
            q2t_[qi] = -1;
            t2q_[tj] = -1;
            if (tb_.expiredNow()) {
                state_[qi] = 0;
                return;
            }
        }

        state_[qi] = 2;
        search(mappedCount, pendingCount - 1, mappedBondCount);
        state_[qi] = 0;
    }

public:
    FixedSizeBondMaximizer(const MolGraph& g1, const MolGraph& g2,
                           const ChemOptions& C, bool induced,
                           TimeBudget& tb, int requiredSize,
                           const std::map<int,int>& incumbent)
        : g1_(g1),
          g2_(g2),
          C_(C),
          induced_(induced),
          tb_(tb),
          requiredSize_(requiredSize),
          q2t_(g1.n, -1),
          t2q_(g2.n, -1),
          bestQ2T_(g1.n, -1),
          state_(g1.n, 0),
          compatTargets_(g1.n) {
        for (const auto& entry : incumbent) {
            if (entry.first < 0 || entry.first >= g1.n) continue;
            if (entry.second < 0 || entry.second >= g2.n) continue;
            bestQ2T_[entry.first] = entry.second;
        }
        bestBondCount_ = countMappedBonds(g1_, incumbent);
        for (int qi = 0; qi < g1.n; ++qi) {
            for (int tj = 0; tj < g2.n; ++tj) {
                if (atomsCompatFast(g1_, qi, g2_, tj, C_)) {
                    compatTargets_[qi].push_back(tj);
                }
            }
        }
    }

    std::map<int,int> run() {
        search(0, g1_.n, 0);
        return materialize(bestQ2T_);
    }
};

// ---------------------------------------------------------------------------
// ComponentSlice + helpers for disconnected MCS support
// ---------------------------------------------------------------------------
struct ComponentSlice {
    MolGraph graph;
    std::vector<int> atomIndices;
};

inline std::vector<ComponentSlice> splitGraphComponents(const MolGraph& g) {
    if (g.n == 0) return {};
    if (countComponents(g) <= 1) {
        std::vector<int> atoms(g.n);
        std::iota(atoms.begin(), atoms.end(), 0);
        return {{g, std::move(atoms)}};
    }

    std::vector<int> compId(g.n, -1);
    int componentCount = 0;
    for (int start = 0; start < g.n; ++start) {
        if (compId[start] >= 0) continue;
        std::deque<int> bfs;
        bfs.push_back(start);
        compId[start] = componentCount;
        while (!bfs.empty()) {
            int atom = bfs.front();
            bfs.pop_front();
            for (int neighbor : g.neighbors[atom]) {
                if (compId[neighbor] >= 0) continue;
                compId[neighbor] = componentCount;
                bfs.push_back(neighbor);
            }
        }
        ++componentCount;
    }

    std::vector<ComponentSlice> components;
    components.reserve(componentCount);
    for (int component = 0; component < componentCount; ++component) {
        std::vector<int> atoms;
        for (int atom = 0; atom < g.n; ++atom) {
            if (compId[atom] == component) atoms.push_back(atom);
        }
        std::sort(atoms.begin(), atoms.end());
        components.push_back({extractSubgraph(g, atoms), std::move(atoms)});
    }
    return components;
}

inline std::map<int,int> liftLocalMapping(const std::map<int,int>& localMapping,
                                          const std::vector<int>& sourceAtoms,
                                          const std::vector<int>& targetAtoms) {
    std::map<int,int> lifted;
    for (const auto& [localSource, localTarget] : localMapping) {
        if (localSource < 0 || localSource >= static_cast<int>(sourceAtoms.size())) continue;
        if (localTarget < 0 || localTarget >= static_cast<int>(targetAtoms.size())) continue;
        lifted[sourceAtoms[localSource]] = targetAtoms[localTarget];
    }
    return lifted;
}

inline std::vector<int> removeMappedTargetAtoms(const std::vector<int>& availableAtoms,
                                                const std::map<int,int>& liftedMapping,
                                                int targetAtomCount) {
    if (liftedMapping.empty()) return availableAtoms;
    std::vector<uint8_t> removed(targetAtomCount, 0);
    for (const auto& entry : liftedMapping) {
        if (entry.second >= 0 && entry.second < targetAtomCount) {
            removed[entry.second] = 1;
        }
    }
    std::vector<int> residual;
    residual.reserve(availableAtoms.size());
    for (int atom : availableAtoms) {
        if (!removed[atom]) residual.push_back(atom);
    }
    return residual;
}

// ---------------------------------------------------------------------------
// Greedy probe (Level 0.75)
// ---------------------------------------------------------------------------
inline std::map<int,int> greedyProbe(const MolGraph& g1, const MolGraph& g2,
                                      const ChemOptions& C,
                                      int templateFuzzyAtoms = 0) {
    int n1 = g1.n, n2 = g2.n;
    if (n1 == 0 || n2 == 0) return {};

    // Precompute compatibility counts so the probe starts from the
    // chemically rarest / most constrained atoms instead of only relying on
    // ring membership or degree. This is especially important on
    // medium-sized medicinal-chemistry graphs where a weak early mapping
    // amplifies the later search tree.
    std::vector<int> compatCount(n1, 0);
    for (int qi = 0; qi < n1; ++qi)
        for (int tj = 0; tj < n2; ++tj)
            if (atomsCompatFast(g1, qi, g2, tj, C)) compatCount[qi]++;

    // Sort query atoms: fewest compatible targets first, then ring atoms,
    // then descending degree, then heteroatoms before carbon.
    std::vector<int> order(n1);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (compatCount[a] != compatCount[b]) return compatCount[a] < compatCount[b];
        if (g1.ring[a] != g1.ring[b]) return g1.ring[a] > g1.ring[b];
        if (g1.degree[a] != g1.degree[b]) return g1.degree[a] > g1.degree[b];
        bool aHet = g1.atomicNum[a] != 6 && g1.atomicNum[a] != 1;
        bool bHet = g1.atomicNum[b] != 6 && g1.atomicNum[b] != 1;
        if (aHet != bHet) return aHet > bHet;
        return a < b;
    });

    std::vector<uint8_t> usedT(n2, 0);
    std::map<int,int> mapping;

    for (int idx = 0; idx < n1; ++idx) {
        int qi = order[idx];
        int bestTj = -1, bestScore = -1;
        for (int tj = 0; tj < n2; ++tj) {
            if (usedT[tj]) continue;
            if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
            bool ok = true;
            int mappedSupport = 0;
            for (int qk : g1.neighbors[qi]) {
                auto it = mapping.find(qk);
                if (it == mapping.end()) continue;
                int tk = it->second;
                int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                if ((qOrd==0) != (tOrd==0)) { ok = false; break; }
                if (qOrd != 0 && !bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { ok = false; break; }
                if (qOrd != 0) mappedSupport++;
            }
            if (!ok) continue;
            int degreeSlack = std::abs(g2.degree[tj] - g1.degree[qi]);
            int score = mappedSupport * 200
                      + (g2.ring[tj] == g1.ring[qi] ? 80 : 0)
                      + (g2.aromatic[tj] == g1.aromatic[qi] ? 40 : 0)
                      + std::min(g1.degree[qi], g2.degree[tj]) * 4
                      - degreeSlack * 3
                      - compatCount[qi];
            if (score > bestScore) { bestScore = score; bestTj = tj; }
        }
        if (bestTj >= 0) { mapping[qi] = bestTj; usedT[bestTj] = true; }
    }

    // Fuzzy template matching: accept the mapping if the size difference
    // between either molecule and the mapping is within the tolerance.
    // If the mapping is too far from the smaller molecule, we still
    // return whatever we found -- the tolerance only controls whether
    // the greedy probe is considered "good enough" for early exit upstream.
    // The actual filtering/acceptance happens at the call site; here we
    // just allow unmatched atoms up to the tolerance by not pruning them.
    // (No pruning needed here -- the tolerance is applied at the call site.)

    return mapping;
}

// ---------------------------------------------------------------------------
// Greedy atom extend
// ---------------------------------------------------------------------------
inline std::map<int,int> greedyAtomExtend(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& seed,
    const ChemOptions& C, const MCSOptions& M) {
    int n1 = g1.n, n2 = g2.n;
    std::vector<int> q2t(n1, -1), t2q(n2, -1);
    for (auto& [k,v] : seed) { q2t[k] = v; t2q[v] = k; }
    bool progress = true;
    while (progress) {
        progress = false;
        for (int qi = 0; qi < n1; ++qi) {
            if (q2t[qi] >= 0) continue;
            bool onFrontier = false;
            for (int nb : g1.neighbors[qi])
                if (q2t[nb] >= 0) { onFrontier = true; break; }
            if (!onFrontier) continue;

            int bestTj = -1, bestScore = -1;
            for (int tj = 0; tj < n2; ++tj) {
                if (t2q[tj] >= 0) continue;
                if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
                bool consistent = true;
                for (int qk : g1.neighbors[qi]) {
                    if (q2t[qk] < 0) continue;
                    int tk = q2t[qk];
                    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                    if (qOrd != 0 && tOrd != 0) {
                        if (!bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                    } else if (M.induced && ((qOrd!=0) != (tOrd!=0))) { consistent = false; break; }
                }
                if (!consistent) continue;
                int score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0)
                          + std::min(g1.degree[qi], g2.degree[tj]);
                if (score > bestScore) { bestScore = score; bestTj = tj; }
            }
            if (bestTj >= 0) {
                q2t[qi] = bestTj; t2q[bestTj] = qi;
                progress = true;
            }
        }
    }
    std::map<int,int> result;
    for (int i = 0; i < n1; ++i) if (q2t[i] >= 0) result[i] = q2t[i];
    return result;
}

// ---------------------------------------------------------------------------
// Label-frequency pruning check
// ---------------------------------------------------------------------------
inline bool isPruned(int curSize, int bestSize, const int* qLF, const int* tLF, int freqSize) {
    int potential = curSize;
    for (int lbl = 0; lbl < freqSize; ++lbl)
        if (qLF[lbl] > 0) potential += std::min(qLF[lbl], tLF[lbl]);
    return potential <= bestSize;
}

// ---------------------------------------------------------------------------
// Build frontier for McGregor DFS
// ---------------------------------------------------------------------------
inline int buildFrontier(const MolGraph& g1, const std::map<int,int>& cur,
                          const std::vector<uint8_t>& usedQ,
                          std::vector<uint8_t>& inFrontier,
                          int* frontierBuf, bool connectedOnly = true) {
    int count = 0;
    for (auto& [qk, tv] : cur)
        for (int qn : g1.neighbors[qk])
            if (!usedQ[qn] && !inFrontier[qn] && !cur.count(qn)) {
                inFrontier[qn] = true;
                frontierBuf[count++] = qn;
            }
    // Only jump to disconnected atoms if disconnected MCS is allowed
    if (count == 0 && !connectedOnly)
        for (int i = 0; i < g1.n; ++i)
            if (!usedQ[i] && !cur.count(i)) frontierBuf[count++] = i;
    for (int f = 0; f < count; ++f) inFrontier[frontierBuf[f]] = false;
    return count;
}

// ---------------------------------------------------------------------------
// Find best candidate for McGregor
// ---------------------------------------------------------------------------
inline int findBestCandidate(
    const MolGraph& g1, const MolGraph& g2, const ChemOptions& C,
    const int* q2tMap, const std::vector<uint8_t>& usedT,
    const std::vector<std::vector<int>>& qNLF1,
    const std::vector<std::vector<int>>& tNLF1,
    const std::vector<std::vector<int>>& qNLF2,
    const std::vector<std::vector<int>>& tNLF2,
    const std::vector<std::vector<int>>& qNLF3,
    const std::vector<std::vector<int>>& tNLF3,
    bool useTwoHopNLF, bool useThreeHopNLF,
    int frontierCount, const int* frontierBuf,
    int* candBuf, int* bestCandBuf,
    int& outQi, int& outCandCount,
    const std::vector<std::vector<int>>* compatTargets = nullptr) {

    int bestQi = -1, bestCandSize = INT_MAX, bestCandCount = 0;
    for (int fi = 0; fi < frontierCount; ++fi) {
        int qi = frontierBuf[fi];
        int candCount = 0;
        // When pre-indexed compatibility is available, iterate only compatible
        // targets instead of scanning all n2 atoms — O(compat) vs O(n2).
        if (compatTargets && qi < static_cast<int>(compatTargets->size())) {
            for (int tj : (*compatTargets)[qi]) {
                if (usedT[tj]) continue;
                if (!nlfCheckOk(qi, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                useTwoHopNLF, useThreeHopNLF)) continue;
                bool ok = true;
                for (int qk : g1.neighbors[qi]) {
                    int tl = q2tMap[qk];
                    if (tl < 0) continue;
                    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tl);
                    if ((qOrd==0) != (tOrd==0)) { ok = false; break; }
                    if (qOrd != 0 && !bondsCompatible(g1, qi, qk, g2, tj, tl, C)) { ok = false; break; }
                }
                if (ok) candBuf[candCount++] = tj;
            }
        } else {
            for (int tj = 0; tj < g2.n; ++tj) {
                if (usedT[tj]) continue;
                if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
                if (!nlfCheckOk(qi, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                useTwoHopNLF, useThreeHopNLF)) continue;
                bool ok = true;
                for (int qk : g1.neighbors[qi]) {
                    int tl = q2tMap[qk];
                    if (tl < 0) continue;
                    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tl);
                    if ((qOrd==0) != (tOrd==0)) { ok = false; break; }
                    if (qOrd != 0 && !bondsCompatible(g1, qi, qk, g2, tj, tl, C)) { ok = false; break; }
                }
                if (ok) candBuf[candCount++] = tj;
            }
        }
        if (candCount == 0) continue;
        if (candCount < bestCandSize) {
            bestCandSize = candCount;
            bestQi = qi;
            bestCandCount = candCount;
            std::memcpy(bestCandBuf, candBuf, candCount * sizeof(int));
        }
    }
    outQi = bestQi;
    outCandCount = bestCandCount;
    return bestQi;
}

// ---------------------------------------------------------------------------
// Undo forced assignments
// ---------------------------------------------------------------------------
inline void undoForcedAssignments(
    int forcedCount, const int* forcedQ, const int* forcedT,
    std::map<int,int>& cur,
    std::vector<uint8_t>& usedQ, std::vector<uint8_t>& usedT,
    int* qLabelFreq, int* tLabelFreq,
    const int* jointQ, const int* jointT,
    int* q2tMap) {
    for (int f = forcedCount - 1; f >= 0; --f) {
        int fq = forcedQ[f], ft = forcedT[f];
        qLabelFreq[jointQ[fq]]++;
        tLabelFreq[jointT[ft]]++;
        cur.erase(fq);
        usedQ[fq] = false;
        usedT[ft] = false;
        if (q2tMap) q2tMap[fq] = -1;
    }
}

// ---------------------------------------------------------------------------
// McGregor DFS (Level 4 -- atom-frontier)
// ---------------------------------------------------------------------------
inline void mcGregorDFS(
    const MolGraph& g1, const MolGraph& g2, const ChemOptions& C,
    std::map<int,int>& cur, std::map<int,int>& best,
    const std::vector<std::vector<int>>& qNLF1,
    const std::vector<std::vector<int>>& tNLF1,
    const std::vector<std::vector<int>>& qNLF2,
    const std::vector<std::vector<int>>& tNLF2,
    const std::vector<std::vector<int>>& qNLF3,
    const std::vector<std::vector<int>>& tNLF3,
    bool useTwoHopNLF, bool useThreeHopNLF,
    TimeBudget& tb, int64_t localDeadlineNs, int depth,
    std::vector<uint8_t>& usedQ, std::vector<uint8_t>& usedT,
    int* q2tMap,
    int* candBuf, int* bestCandBuf,
    int* qLabelFreq, int* tLabelFreq, int freqSize,
    std::vector<uint8_t>& inFrontier, int* frontierBuf,
    const int* jointQ, const int* jointT,
    const std::vector<std::vector<int>>* compatTargets = nullptr) {

    // Amortized local-deadline check: only call Clock::now() every 1024 iterations
    static thread_local int64_t mcgDfsIter = 0;
    auto localExpired = [&]() {
        if ((++mcgDfsIter & 1023) != 0) return false;
        using Clock = std::chrono::steady_clock;
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now().time_since_epoch()).count() >= localDeadlineNs;
    };

    if (localExpired() || tb.expired()) return;
    if (static_cast<int>(cur.size()) > static_cast<int>(best.size())) {
        best = cur;
    }
    if (isPruned(static_cast<int>(cur.size()), static_cast<int>(best.size()),
                 qLabelFreq, tLabelFreq, freqSize)) return;

    int frontierCount = buildFrontier(g1, cur, usedQ, inFrontier, frontierBuf);
    if (frontierCount == 0) return;

    int bestQi = -1, bestCandCount = 0;
    findBestCandidate(g1, g2, C, q2tMap, usedT, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                      useTwoHopNLF, useThreeHopNLF, frontierCount, frontierBuf,
                      candBuf, bestCandBuf, bestQi, bestCandCount, compatTargets);
    if (bestQi == -1) return;

    // Unit propagation (forced assignment)
    int forcedCount = 0;
    // Reuse caller-provided scratch buffers (sized to g1.n / g2.n)
    thread_local std::vector<int> forcedQ, forcedT;
    forcedQ.resize(g1.n); forcedT.resize(g2.n);
    while (bestCandCount == 1 && !(localExpired() || tb.expired())) {
        int fq = bestQi, ft = bestCandBuf[0];
        cur[fq] = ft; usedQ[fq] = true; usedT[ft] = true;
        q2tMap[fq] = ft;
        qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
        forcedQ[forcedCount] = fq; forcedT[forcedCount] = ft;
        forcedCount++; depth++;
        if (static_cast<int>(cur.size()) > static_cast<int>(best.size())) best = cur;
        if (isPruned(static_cast<int>(cur.size()), static_cast<int>(best.size()),
                     qLabelFreq, tLabelFreq, freqSize)) break;

        frontierCount = buildFrontier(g1, cur, usedQ, inFrontier, frontierBuf);
        if (frontierCount == 0) break;
        findBestCandidate(g1, g2, C, q2tMap, usedT, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                          useTwoHopNLF, useThreeHopNLF, frontierCount, frontierBuf,
                          candBuf, bestCandBuf, bestQi, bestCandCount, compatTargets);
        if (bestQi == -1) break;
    }

    if (bestQi != -1 && bestCandCount > 1) {
        int branchLimit = depth < 5 ? bestCandCount : std::min(bestCandCount, 16);
        for (int i = 0; i < branchLimit; ++i) {
            if (localExpired() || tb.expired()) break;
            int bestTj = bestCandBuf[i];
            cur[bestQi] = bestTj; usedQ[bestQi] = true; usedT[bestTj] = true;
            q2tMap[bestQi] = bestTj;
            qLabelFreq[jointQ[bestQi]]--; tLabelFreq[jointT[bestTj]]--;
            mcGregorDFS(g1, g2, C, cur, best, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                        useTwoHopNLF, useThreeHopNLF, tb, localDeadlineNs, depth + 1,
                        usedQ, usedT, q2tMap, candBuf, bestCandBuf,
                        qLabelFreq, tLabelFreq, freqSize,
                        inFrontier, frontierBuf, jointQ, jointT, compatTargets);
            qLabelFreq[jointQ[bestQi]]++; tLabelFreq[jointT[bestTj]]++;
            q2tMap[bestQi] = -1;
            cur.erase(bestQi); usedQ[bestQi] = false; usedT[bestTj] = false;
        }
    }

    undoForcedAssignments(forcedCount, forcedQ.data(), forcedT.data(), cur, usedQ, usedT,
                          qLabelFreq, tLabelFreq, jointQ, jointT, q2tMap);
}

// ---------------------------------------------------------------------------
// McGregor bond-grow (Level 4 -- bond-oriented)
// ---------------------------------------------------------------------------
inline void mcGregorBondGrow(
    const MolGraph& g1, const MolGraph& g2, const ChemOptions& C,
    std::map<int,int>& cur, std::map<int,int>& best,
    const std::vector<std::vector<int>>& qNLF1,
    const std::vector<std::vector<int>>& tNLF1,
    bool useTwoHopNLF, bool useThreeHopNLF,
    const std::vector<std::vector<int>>& qNLF2,
    const std::vector<std::vector<int>>& tNLF2,
    const std::vector<std::vector<int>>& qNLF3,
    const std::vector<std::vector<int>>& tNLF3,
    TimeBudget& tb, int64_t localDeadlineNs, int depth,
    std::vector<uint8_t>& usedQ, std::vector<uint8_t>& usedT,
    int* qLabelFreq, int* tLabelFreq, int freqSize,
    int* q2tMap, std::vector<uint8_t>& inFrontier, int* frontierBuf,
    int* candBuf, int* bestCandBuf,
    const int* jointQ, const int* jointT) {

    // Amortized local-deadline check: only call Clock::now() every 1024 iterations
    static thread_local int64_t mcgBondIter = 0;
    auto localExpired = [&]() {
        if ((++mcgBondIter & 1023) != 0) return false;
        using Clock = std::chrono::steady_clock;
        return std::chrono::duration_cast<std::chrono::nanoseconds>(
            Clock::now().time_since_epoch()).count() >= localDeadlineNs;
    };

    if (localExpired() || tb.expired()) return;
    if (static_cast<int>(cur.size()) > static_cast<int>(best.size())) best = cur;
    if (isPruned(static_cast<int>(cur.size()), static_cast<int>(best.size()),
                 qLabelFreq, tLabelFreq, freqSize)) return;

    // Find best frontier bond
    int bestQk = -1, bestQi_local = -1, bestCandSize = INT_MAX, bestCandCount = 0;
    for (auto& [qi, mappedTi_v] : cur) {
        int mappedTi = q2tMap[qi];
        for (int qk : g1.neighbors[qi]) {
            if (usedQ[qk] || inFrontier[qk]) continue;
            int candCount = 0;
            for (int tk : g2.neighbors[mappedTi]) {
                if (usedT[tk]) continue;
                int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(mappedTi, tk);
                if ((qOrd==0) != (tOrd==0)) continue;
                if (qOrd != 0 && !bondsCompatible(g1, qi, qk, g2, mappedTi, tk, C)) continue;
                if (!atomsCompatFast(g1, qk, g2, tk, C)) continue;
                if (!nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                useTwoHopNLF, useThreeHopNLF)) continue;
                bool ok = true;
                for (int qn : g1.neighbors[qk]) {
                    if (qn == qi || !usedQ[qn]) continue;
                    int tn = q2tMap[qn];
                    int qOrd2 = g1.bondOrder(qk, qn), tOrd2 = g2.bondOrder(tk, tn);
                    if ((qOrd2==0) != (tOrd2==0)) { ok = false; break; }
                    if (qOrd2 != 0 && !bondsCompatible(g1, qk, qn, g2, tk, tn, C)) { ok = false; break; }
                }
                if (ok) candBuf[candCount++] = tk;
            }
            if (candCount > 0 && candCount < bestCandSize) {
                bestCandSize = candCount; bestQk = qk; bestQi_local = qi;
                bestCandCount = candCount;
                std::memcpy(bestCandBuf, candBuf, candCount * sizeof(int));
            }
        }
    }

    // Fallback: disconnected extension
    if (bestQk == -1) {
        for (int i = 0; i < g1.n; ++i) {
            if (usedQ[i]) continue;
            int candCount = 0;
            for (int tj = 0; tj < g2.n; ++tj) {
                if (usedT[tj]) continue;
                if (!atomsCompatFast(g1, i, g2, tj, C)) continue;
                if (!nlfCheckOk(i, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                useTwoHopNLF, useThreeHopNLF)) continue;
                candBuf[candCount++] = tj;
            }
            if (candCount > 0 && candCount < bestCandSize) {
                bestCandSize = candCount; bestQk = i; bestQi_local = -1;
                bestCandCount = candCount;
                std::memcpy(bestCandBuf, candBuf, candCount * sizeof(int));
            }
        }
    }
    if (bestQk == -1) return;

    // Unit propagation
    int forcedCount = 0;
    // Reuse caller-provided scratch buffers (sized to g1.n / g2.n)
    thread_local std::vector<int> forcedQ, forcedT;
    forcedQ.resize(g1.n); forcedT.resize(g2.n);
    while (bestCandCount == 1 && !(localExpired() || tb.expired())) {
        int fq = bestQk, ft = bestCandBuf[0];
        cur[fq] = ft; usedQ[fq] = true; usedT[ft] = true; q2tMap[fq] = ft;
        qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
        forcedQ[forcedCount] = fq; forcedT[forcedCount] = ft;
        forcedCount++; depth++;
        if (static_cast<int>(cur.size()) > static_cast<int>(best.size())) best = cur;
        if (isPruned(static_cast<int>(cur.size()), static_cast<int>(best.size()),
                     qLabelFreq, tLabelFreq, freqSize)) break;

        bestQk = -1; bestCandSize = INT_MAX; bestCandCount = 0;
        for (auto& [qi, mappedTi_v] : cur) {
            int mappedTi = q2tMap[qi];
            for (int qk : g1.neighbors[qi]) {
                if (usedQ[qk]) continue;
                int candCount = 0;
                for (int tk : g2.neighbors[mappedTi]) {
                    if (usedT[tk]) continue;
                    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(mappedTi, tk);
                    if ((qOrd==0) != (tOrd==0)) continue;
                    if (qOrd != 0 && !bondsCompatible(g1, qi, qk, g2, mappedTi, tk, C)) continue;
                    if (!atomsCompatFast(g1, qk, g2, tk, C)) continue;
                    if (!nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                                    useTwoHopNLF, useThreeHopNLF)) continue;
                    bool ok = true;
                    for (int qn : g1.neighbors[qk]) {
                        if (qn == qi || !usedQ[qn]) continue;
                        int tn = q2tMap[qn];
                        int qOrd2 = g1.bondOrder(qk, qn), tOrd2 = g2.bondOrder(tk, tn);
                        if ((qOrd2==0) != (tOrd2==0)) { ok = false; break; }
                        if (qOrd2 != 0 && !bondsCompatible(g1, qk, qn, g2, tk, tn, C)) { ok = false; break; }
                    }
                    if (ok) candBuf[candCount++] = tk;
                }
                if (candCount > 0 && candCount < bestCandSize) {
                    bestCandSize = candCount; bestQk = qk; bestQi_local = qi;
                    bestCandCount = candCount;
                    std::memcpy(bestCandBuf, candBuf, candCount * sizeof(int));
                }
            }
        }
        if (bestQk == -1) break;
    }

    if (bestQk != -1 && bestCandCount > 1) {
        int branchLimit = depth < 5 ? bestCandCount : std::min(bestCandCount, 16);
        for (int i = 0; i < branchLimit; ++i) {
            if (localExpired() || tb.expired()) break;
            int btj = bestCandBuf[i];
            cur[bestQk] = btj; usedQ[bestQk] = true; usedT[btj] = true; q2tMap[bestQk] = btj;
            qLabelFreq[jointQ[bestQk]]--; tLabelFreq[jointT[btj]]--;
            mcGregorBondGrow(g1, g2, C, cur, best, qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF,
                             qNLF2, tNLF2, qNLF3, tNLF3, tb, localDeadlineNs, depth + 1,
                             usedQ, usedT, qLabelFreq, tLabelFreq, freqSize, q2tMap,
                             inFrontier, frontierBuf, candBuf, bestCandBuf, jointQ, jointT);
            qLabelFreq[jointQ[bestQk]]++; tLabelFreq[jointT[btj]]++;
            cur.erase(bestQk); usedQ[bestQk] = false; usedT[btj] = false; q2tMap[bestQk] = -1;
        }
    }

    undoForcedAssignments(forcedCount, forcedQ.data(), forcedT.data(), cur, usedQ, usedT,
                          qLabelFreq, tLabelFreq, jointQ, jointT, q2tMap);
}

// ---------------------------------------------------------------------------
// McGregor extend entry point (Level 4)
// ---------------------------------------------------------------------------
inline std::map<int,int> mcGregorExtend(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& seed, const ChemOptions& C,
    TimeBudget& tb, int64_t localMillis,
    bool useTwoHopNLF, bool useThreeHopNLF,
    bool connectedOnly = true,
    const std::vector<std::vector<int>>* compatTargets = nullptr) {

    using Clock = std::chrono::steady_clock;
    int64_t localDeadlineNs = std::chrono::duration_cast<std::chrono::nanoseconds>(
        (Clock::now() + std::chrono::milliseconds(localMillis)).time_since_epoch()).count();

    std::map<int,int> best = seed;
    auto& qNLF1 = g1.getNLF1(); auto& tNLF1 = g2.getNLF1();
    static const std::vector<std::vector<int>> emptyNLF;
    auto& qNLF2 = useTwoHopNLF ? g1.getNLF2() : emptyNLF;
    auto& tNLF2 = useTwoHopNLF ? g2.getNLF2() : emptyNLF;
    auto& qNLF3 = useThreeHopNLF ? g1.getNLF3() : emptyNLF;
    auto& tNLF3 = useThreeHopNLF ? g2.getNLF3() : emptyNLF;

    int n1 = g1.n, n2 = g2.n;
    std::vector<uint8_t> usedQ(n1, 0), usedT(n2, 0);
    std::vector<int> candBuf(n2), bestCandBuf(n2);
    for (auto& [k,v] : seed) { usedQ[k] = 1; usedT[v] = 1; }

    int maxLabel = 0;
    for (int i = 0; i < n1; ++i) maxLabel = std::max(maxLabel, g1.label[i]);
    for (int j = 0; j < n2; ++j) maxLabel = std::max(maxLabel, g2.label[j]);
    int freqSize = maxLabel + 1;
    std::vector<int> qLabelFreq(freqSize, 0), tLabelFreq(freqSize, 0);
    std::vector<int> jointQ(n1), jointT(n2);
    for (int i = 0; i < n1; ++i) { jointQ[i] = g1.label[i]; if (!usedQ[i]) qLabelFreq[jointQ[i]]++; }
    for (int j = 0; j < n2; ++j) { jointT[j] = g2.label[j]; if (!usedT[j]) tLabelFreq[jointT[j]]++; }

    std::vector<uint8_t> inFrontier(n1, 0);
    std::vector<int> frontierBuf(n1);

    // Fast bail-out
    if (!seed.empty()) {
        bool hasExtensible = false;
        for (auto& [k,v] : seed) {
            for (int nb : g1.neighbors[k]) {
                if (usedQ[nb]) continue;
                for (int tj = 0; tj < n2; ++tj) {
                    if (usedT[tj]) continue;
                    if (atomsCompatFast(g1, nb, g2, tj, C)) { hasExtensible = true; goto bail_done; }
                }
            }
        }
        bail_done:
        if (!hasExtensible && connectedOnly) return best;
        // If disconnected MCS, skip bail-out — let frontier jump to new component
    }

    // Bond-growth first pass for larger molecules
    if (!seed.empty() && n2 >= 20) {
        std::vector<int> q2tMap(n1, -1);
        for (auto& [k,v] : seed) q2tMap[k] = v;
        std::map<int,int> bondBest = seed;
        int64_t bondDeadlineNs = std::chrono::duration_cast<std::chrono::nanoseconds>(
            (Clock::now() + std::chrono::milliseconds(std::max<int64_t>(1, std::min(localMillis / 10, (int64_t)20)))).time_since_epoch()).count();

        auto usedQCopy = usedQ;
        auto usedTCopy = usedT;
        auto qLFCopy = qLabelFreq;
        auto tLFCopy = tLabelFreq;
        auto inFCopy = inFrontier;
        std::map<int,int> curCopy = seed;

        mcGregorBondGrow(g1, g2, C, curCopy, bondBest, qNLF1, tNLF1,
                         useTwoHopNLF, useThreeHopNLF, qNLF2, tNLF2, qNLF3, tNLF3,
                         tb, bondDeadlineNs, 0,
                         usedQCopy, usedTCopy, qLFCopy.data(), tLFCopy.data(),
                         freqSize, q2tMap.data(), inFCopy, frontierBuf.data(),
                         candBuf.data(), bestCandBuf.data(), jointQ.data(), jointT.data());
        std::fill(inFrontier.begin(), inFrontier.end(), uint8_t(0));
        if (bondBest.size() > best.size()) best = bondBest;
    }

    {
        std::map<int,int> curDFS = seed;
        std::vector<int> q2tDFS(n1, -1);
        for (auto& [k,v] : seed) q2tDFS[k] = v;
        mcGregorDFS(g1, g2, C, curDFS, best, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
                    useTwoHopNLF, useThreeHopNLF, tb, localDeadlineNs, 0,
                    usedQ, usedT, q2tDFS.data(), candBuf.data(), bestCandBuf.data(),
                    qLabelFreq.data(), tLabelFreq.data(), freqSize,
                    inFrontier, frontierBuf.data(), jointQ.data(), jointT.data(),
                    compatTargets);
    }
    return best;
}

// ===========================================================================
// GraphBuilder -- seeds for the MCS funnel
// ===========================================================================
class GraphBuilder {
    const MolGraph& g1_;
    const MolGraph& g2_;
    const ChemOptions& C_;
    bool induced_;

    static constexpr int64_t MAX_NODE_LIMIT  = 500000;
    static constexpr int64_t SEED_EXTEND_NODE_LIMIT = 100000;
    static constexpr int64_t MCSPLIT_NODE_LIMIT = 200000;

    struct PGNode { int qi, tj; };

    // Precomputed compatibility: compatTargets_[qi] = list of compatible target atoms.
    // Built once in constructor, reused by seedExtendMCS, McSplit, BK.
    std::vector<std::vector<int>> compatTargets_;

public:
    GraphBuilder(const MolGraph& g1, const MolGraph& g2,
                 const ChemOptions& C, bool induced)
        : g1_(g1), g2_(g2), C_(C), induced_(induced)
    {
        // Precompute atom compatibility matrix once (v6.5.3 perf fix).
        // Phase 2.1: Bucketed compatibility — bucket target atoms by the most
        // discriminating fast-check fields (atomicNum, aromatic) so each query
        // atom only tests candidates from its matching bucket instead of the
        // full O(n1 x n2) cross-product.  Average work becomes O(n1 x B)
        // where B = avg bucket size (typically n2 / #distinct_element_types).
        int n1 = g1.n, n2 = g2.n;
        compatTargets_.resize(n1);

        // --- Pass 1: build target buckets keyed by (atomicNum, aromatic) ---
        // When matchAtomType is false, all targets share key 0 (one bucket,
        // no benefit but no regression).  Aromatic bit is only discriminating
        // in STRICT aromaticity mode; otherwise folded out.
        bool useAtomNum  = C.matchAtomType;
        bool useAromatic = (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT);

        std::unordered_map<int, std::vector<int>> targetBuckets;
        targetBuckets.reserve(n2 < 128 ? 16 : 64);
        for (int j = 0; j < n2; ++j) {
            int key = 0;
            if (useAtomNum)  key = g2.atomicNum[j];
            // Shift left by 1 and pack the aromatic bit in the LSB.
            if (useAromatic) key = (key << 1) | (g2.aromatic[j] ? 1 : 0);
            targetBuckets[key].push_back(j);
        }

        // --- Pass 2: for each query atom, look up matching bucket(s) ---
        bool hasTautomer = C.tautomerAware
            && !g1.tautomerClass.empty() && !g2.tautomerClass.empty();

        for (int i = 0; i < n1; ++i) {
            // Determine the bucket key this query atom would match.
            int qKey = 0;
            if (useAtomNum)  qKey = g1.atomicNum[i];
            if (useAromatic) qKey = (qKey << 1) | (g1.aromatic[i] ? 1 : 0);

            // If tautomer-aware mode is active and this query atom sits in a
            // tautomeric region, it might match C/N/O/S/Se targets regardless
            // of its own element.  Fall back to scanning ALL buckets for such
            // atoms to preserve exact semantics.
            bool needFullScan = hasTautomer && g1.tautomerClass[i] != -1;

            if (!needFullScan) {
                // Fast path: only inspect the single matching bucket.
                auto it = targetBuckets.find(qKey);
                if (it != targetBuckets.end()) {
                    for (int j : it->second)
                        if (atomsCompatFast(g1, i, g2, j, C))
                            compatTargets_[i].push_back(j);
                }
            } else {
                // Tautomer fallback: scan all buckets (preserves correctness).
                for (auto& [bucketKey, targets] : targetBuckets) {
                    for (int j : targets)
                        if (atomsCompatFast(g1, i, g2, j, C))
                            compatTargets_[i].push_back(j);
                }
                // Tautomer path may yield out-of-order targets; sort to keep
                // deterministic output identical to the original loop order.
                std::sort(compatTargets_[i].begin(), compatTargets_[i].end());
            }
        }
    }

    /** Access precomputed compatibility targets for a query atom. */
    const std::vector<int>& compatTargets(int qi) const { return compatTargets_[qi]; }

    /** Access the full precomputed compatibility table (for McGregor acceleration). */
    const std::vector<std::vector<int>>& allCompatTargets() const { return compatTargets_; }

    // -- Greedy bond extend (used by seed-and-extend) ----------------------
    static int greedyBondExtend(const MolGraph& g1, const MolGraph& g2,
                                 const ChemOptions& C,
                                 int* q2t, int* t2q, int n1, int n2, bool induced) {
        int mapSize = 0;
        for (int i = 0; i < n1; ++i) if (q2t[i] >= 0) mapSize++;
        bool progress = true;
        while (progress) {
            progress = false;
            for (int qi = 0; qi < n1; ++qi) {
                if (q2t[qi] >= 0) continue;
                int bestTj = -1, bestScore = -1;
                for (int nb : g1.neighbors[qi]) {
                    if (q2t[nb] < 0) continue;
                    int tNb = q2t[nb];
                    for (int tj : g2.neighbors[tNb]) {
                        if (t2q[tj] >= 0) continue;
                        if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
                        if (!bondsCompatible(g1, qi, nb, g2, tj, tNb, C)) continue;
                        bool consistent = true;
                        for (int qk : g1.neighbors[qi]) {
                            if (qk == nb || q2t[qk] < 0) continue;
                            int tk = q2t[qk];
                            int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                            if (qOrd != 0 && tOrd != 0) {
                                if (!bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                            } else if (induced && ((qOrd!=0) != (tOrd!=0))) { consistent = false; break; }
                        }
                        if (!consistent) continue;
                        int score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0) + std::min(g1.degree[qi], g2.degree[tj]);
                        if (score > bestScore) { bestScore = score; bestTj = tj; }
                    }
                }
                if (bestTj >= 0) { q2t[qi] = bestTj; t2q[bestTj] = qi; mapSize++; progress = true; }
            }
        }
        return mapSize;
    }

    // -- Backtrack bond extend (used by seed-and-extend) -------------------
    static int mappedNeighborCount(const MolGraph& g1, const int* q2t, int qi) {
        int count = 0;
        for (int nb : g1.neighbors[qi]) if (q2t[nb] >= 0) count++;
        return count;
    }

    static void collectBondExtendCandidates(
        const MolGraph& g1, const MolGraph& g2, const ChemOptions& C, bool induced,
        int qi, const int* q2t, const int* t2q, std::vector<int>& out) {

        out.clear();

        int seedQk = -1;
        int seedTk = -1;
        int mappedNbrs = 0;
        for (int qk : g1.neighbors[qi]) {
            int tk = q2t[qk];
            if (tk < 0) continue;
            mappedNbrs++;
            if (seedTk < 0 || g2.degree[tk] < g2.degree[seedTk]) {
                seedQk = qk;
                seedTk = tk;
            }
        }
        if (seedTk < 0) return;

        std::vector<std::pair<int, int>> scored;
        scored.reserve(g2.degree[seedTk]);

        for (int tj : g2.neighbors[seedTk]) {
            if (t2q[tj] >= 0) continue;
            if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
            if (!bondsCompatible(g1, qi, seedQk, g2, tj, seedTk, C)) continue;

            bool consistent = true;
            for (int qk : g1.neighbors[qi]) {
                if (qk == seedQk || q2t[qk] < 0) continue;
                int tk = q2t[qk];
                int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                if (qOrd != 0 && tOrd != 0) {
                    if (!bondsCompatible(g1, qi, qk, g2, tj, tk, C)) {
                        consistent = false;
                        break;
                    }
                } else if (induced && ((qOrd != 0) != (tOrd != 0))) {
                    consistent = false;
                    break;
                }
            }
            if (!consistent) continue;

            int score = mappedNbrs * 100
                + (g1.ring[qi] && g2.ring[tj] ? 50 : 0)
                + std::min(g1.degree[qi], g2.degree[tj]);
            scored.push_back({score, tj});
        }

        std::sort(scored.begin(), scored.end(),
                  [](const auto& a, const auto& b) {
                      if (a.first != b.first) return a.first > b.first;
                      return a.second < b.second;
                  });
        out.reserve(scored.size());
        for (const auto& entry : scored) out.push_back(entry.second);
    }

    static void backtrackBondExtend(
        const MolGraph& g1, const MolGraph& g2, const ChemOptions& C,
        int* q2t, int* t2q, int n1, int n2, bool induced,
        int* bestQ2T, int& bestSize, TimeBudget& tb, int64_t& nodeCount, int64_t nodeLimit) {
        int curSize = 0;
        for (int i = 0; i < n1; ++i) if (q2t[i] >= 0) curSize++;
        if (curSize > bestSize) { bestSize = curSize; std::memcpy(bestQ2T, q2t, n1 * sizeof(int)); }
        if (nodeCount > nodeLimit || tb.expired()) return;

        int frontier = 0;
        int bestQi = -1, bestConstraint = -1, bestCandCount = INT_MAX;
        std::vector<int> scratchCandidates;
        std::vector<int> bestCandidates;
        for (int qi = 0; qi < n1; ++qi) {
            if (q2t[qi] >= 0) continue;
            int mappedNb = mappedNeighborCount(g1, q2t, qi);
            if (mappedNb == 0) continue;
            collectBondExtendCandidates(g1, g2, C, induced, qi, q2t, t2q, scratchCandidates);
            if (scratchCandidates.empty()) continue;
            frontier++;
            int constraint = mappedNb * 100 + (g1.ring[qi] ? 50 : 0) + g1.degree[qi];
            int candCount = static_cast<int>(scratchCandidates.size());
            if (candCount < bestCandCount
                || (candCount == bestCandCount && constraint > bestConstraint)) {
                bestCandCount = candCount;
                bestConstraint = constraint;
                bestQi = qi;
                bestCandidates = scratchCandidates;
            }
        }
        if (curSize + frontier <= bestSize) return;
        if (bestQi < 0) return;

        int qi = bestQi;
        for (int tj : bestCandidates) {
            nodeCount++;
            if (nodeCount > nodeLimit || tb.expired()) return;
            q2t[qi] = tj; t2q[tj] = qi;
            backtrackBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced, bestQ2T, bestSize, tb, nodeCount, nodeLimit);
            q2t[qi] = -1; t2q[tj] = -1;
            if (nodeCount > nodeLimit || tb.expired()) return;
        }
    }

    // -- Seed-and-extend (Level 1.5) ---------------------------------------
    /// Phase 2.4: seedExtendMCS now delegates to the flat variant to
    /// avoid code duplication.  Materializes std::map only at the boundary.
    std::map<int,int> seedExtendMCS(TimeBudget& tb, int upperBound) {
        auto [flat, sz] = seedExtendMCSFlat(tb, upperBound);
        return flatToMap(flat.data(), g1_.n);
    }

    /// Phase 2.4: Flat-array variant of seedExtendMCS.
    /// Returns (q2t flat array where -1=unmapped, mapping size).
    /// Avoids constructing std::map<int,int> for internal pipeline consumers.
    std::pair<std::vector<int>, int> seedExtendMCSFlat(TimeBudget& tb, int upperBound) {
        int n1 = g1_.n, n2 = g2_.n;
        if (n1 == 0 || n2 == 0) return {std::vector<int>(std::max(n1, 1), -1), 0};

        int maxLabel = 0;
        for (int i = 0; i < n1; ++i) maxLabel = std::max(maxLabel, g1_.label[i]);
        for (int j = 0; j < n2; ++j) maxLabel = std::max(maxLabel, g2_.label[j]);
        std::vector<int> labelFreqG1(maxLabel + 1, 0), labelFreqG2(maxLabel + 1, 0);
        for (int i = 0; i < n1; ++i) labelFreqG1[g1_.label[i]]++;
        for (int j = 0; j < n2; ++j) labelFreqG2[g2_.label[j]]++;

        const auto& compatT = compatTargets_;

        struct Bond { int u, v; };
        std::vector<Bond> g1Bonds;
        for (int i = 0; i < n1; ++i)
            for (int nb : g1_.neighbors[i])
                if (nb > i) g1Bonds.push_back({i, nb});

        int g1BondCount = static_cast<int>(g1Bonds.size());

        std::vector<std::pair<int,int>> seedOrder(g1BondCount);
        for (int b = 0; b < g1BondCount; ++b) {
            int u = g1Bonds[b].u, v = g1Bonds[b].v;
            int totalAtoms = n1 + n2;
            int uFreq = labelFreqG1[g1_.label[u]]
                + (g1_.label[u] < static_cast<int>(labelFreqG2.size()) ? labelFreqG2[g1_.label[u]] : 0);
            int vFreq = labelFreqG1[g1_.label[v]]
                + (g1_.label[v] < static_cast<int>(labelFreqG2.size()) ? labelFreqG2[g1_.label[v]] : 0);
            int score = (totalAtoms * 2) / std::max(1, uFreq) + (totalAtoms * 2) / std::max(1, vFreq)
                + (n2 * 2) / std::max(1, static_cast<int>(compatT[u].size()))
                + (n2 * 2) / std::max(1, static_cast<int>(compatT[v].size()))
                + g1_.degree[u] + g1_.degree[v] + (g1_.ring[u] ? 10 : 0) + (g1_.ring[v] ? 10 : 0);
            seedOrder[b] = {b, score};
        }
        std::sort(seedOrder.begin(), seedOrder.end(),
            [](auto& a, auto& b) { return a.second > b.second; });

        std::vector<int> bestQ2T(n1, -1);
        int bestSize = 0;
        int64_t nodeCount = 0;
        int minN = std::min(n1, n2);
        int maxSeeds = std::min(g1BondCount, 10);
        std::vector<uint8_t> atomUsedAsSeed(n1, 0);
        bool seedRsActive = (g1_.ringSystemCount() > 0);
        struct SeedPairKey {
            int64_t orbitPair;
            uint64_t rsSigU;
            uint64_t rsSigV;
            bool operator<(const SeedPairKey& o) const {
                if (orbitPair != o.orbitPair) return orbitPair < o.orbitPair;
                if (rsSigU != o.rsSigU) return rsSigU < o.rsSigU;
                return rsSigV < o.rsSigV;
            }
        };
        std::set<SeedPairKey> triedOrbitPairs;
        int seedsTried = 0;
        std::vector<int> q2t(n1, -1), t2q(n2, -1);

        for (int si = 0; si < g1BondCount && seedsTried < maxSeeds; ++si) {
            if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
            int bondIdx = seedOrder[si].first;
            int qu = g1Bonds[bondIdx].u, qv = g1Bonds[bondIdx].v;
            int oA = std::min(g1_.orbit[qu], g1_.orbit[qv]);
            int oB = std::max(g1_.orbit[qu], g1_.orbit[qv]);
            int64_t pairKey = (static_cast<int64_t>(oA) << 32) | oB;
            uint64_t rsU = 0, rsV = 0;
            if (seedRsActive) {
                rsU = g1_.ringSystemSig(g1_.ringSystemOf(qu));
                rsV = g1_.ringSystemSig(g1_.ringSystemOf(qv));
                if (g1_.orbit[qu] > g1_.orbit[qv] || (g1_.orbit[qu] == g1_.orbit[qv] && rsU > rsV))
                    std::swap(rsU, rsV);
            }
            if (!triedOrbitPairs.insert(SeedPairKey{pairKey, rsU, rsV}).second) continue;
            if (atomUsedAsSeed[qu] && atomUsedAsSeed[qv]) continue;
            atomUsedAsSeed[qu] = true; atomUsedAsSeed[qv] = true; seedsTried++;

            for (int ta : compatT[qu]) {
                if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
                for (int tb2 : g2_.neighbors[ta]) {
                    nodeCount++;
                    if (nodeCount > MAX_NODE_LIMIT || tb.expired()) break;
                    if (!atomsCompatFast(g1_, qv, g2_, tb2, C_)) continue;
                    if (!bondsCompatible(g1_, qu, qv, g2_, ta, tb2, C_)) continue;
                    std::memset(q2t.data(), -1, n1 * sizeof(int));
                    std::memset(t2q.data(), -1, n2 * sizeof(int));
                    q2t[qu] = ta; t2q[ta] = qu; q2t[qv] = tb2; t2q[tb2] = qv;
                    int mapSize = greedyBondExtend(g1_, g2_, C_, q2t.data(), t2q.data(), n1, n2, induced_);
                    nodeCount += mapSize;
                    if (mapSize > bestSize) { bestSize = mapSize; std::memcpy(bestQ2T.data(), q2t.data(), n1 * sizeof(int)); }
                }
            }
            for (int ta : compatT[qv]) {
                if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
                for (int tb2 : g2_.neighbors[ta]) {
                    nodeCount++;
                    if (nodeCount > MAX_NODE_LIMIT || tb.expired()) break;
                    if (!atomsCompatFast(g1_, qu, g2_, tb2, C_)) continue;
                    if (!bondsCompatible(g1_, qu, qv, g2_, tb2, ta, C_)) continue;
                    std::memset(q2t.data(), -1, n1 * sizeof(int));
                    std::memset(t2q.data(), -1, n2 * sizeof(int));
                    q2t[qu] = tb2; t2q[tb2] = qu; q2t[qv] = ta; t2q[ta] = qv;
                    int mapSize = greedyBondExtend(g1_, g2_, C_, q2t.data(), t2q.data(), n1, n2, induced_);
                    nodeCount += mapSize;
                    if (mapSize > bestSize) { bestSize = mapSize; std::memcpy(bestQ2T.data(), q2t.data(), n1 * sizeof(int)); }
                }
            }
            if (bestSize >= upperBound || bestSize >= minN) break;
            if (bestSize >= (upperBound * 19 / 20) && nodeCount > SEED_EXTEND_NODE_LIMIT / 2) break;
        }

        if (bestSize > 0 && bestSize < minN && !tb.expired() && nodeCount < MAX_NODE_LIMIT) {
            std::vector<int> q2tBT(n1, -1), t2qBT(n2, -1);
            for (int i = 0; i < n1; ++i) if (bestQ2T[i] >= 0) { q2tBT[i] = bestQ2T[i]; t2qBT[bestQ2T[i]] = i; }
            backtrackBondExtend(g1_, g2_, C_, q2tBT.data(), t2qBT.data(), n1, n2, induced_,
                                bestQ2T.data(), bestSize, tb, nodeCount, MAX_NODE_LIMIT);
            bestSize = 0;
            for (int i = 0; i < n1; ++i) if (bestQ2T[i] >= 0) bestSize++;
        }

        return {std::move(bestQ2T), bestSize};
    }

    // -- McSplit (Level 2) -------------------------------------------------
private:
    static bool checkBondCompat(const MolGraph& g1, const MolGraph& g2,
                                 const ChemOptions& C, bool induced,
                                 int qi, int tj, const int* q2t) {
        for (int qk : g1.neighbors[qi]) {
            if (q2t[qk] < 0) continue;
            int tk = q2t[qk];
            int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
            if (qOrd != 0 && tOrd != 0) {
                if (!bondsCompatible(g1, qi, qk, g2, tj, tk, C)) return false;
            } else if (induced && ((qOrd!=0) != (tOrd!=0))) return false;
        }
        return true;
    }

    // Hardware-accelerated bitset for McSplit partition classes.
    // Uses uint64_t words + ctz64 for ~64x faster nextSetBit vs vector<bool>.
    struct BitClass {
        std::vector<uint64_t> words;
        int card;
        int n_bits;
        BitClass() : card(0), n_bits(0) {}
        explicit BitClass(int n) : words((n + 63) / 64, 0), card(0), n_bits(n) {}
        void set(int i) {
            if (!(words[i >> 6] & (1ULL << (i & 63)))) {
                words[i >> 6] |= (1ULL << (i & 63));
                card++;
            }
        }
        void clear(int i) {
            if (words[i >> 6] & (1ULL << (i & 63))) {
                words[i >> 6] &= ~(1ULL << (i & 63));
                card--;
            }
        }
        bool get(int i) const { return (words[i >> 6] & (1ULL << (i & 63))) != 0; }
        bool empty() const { return card == 0; }
        int cardinality() const { return card; }
        int nextSetBit(int from) const {
            if (from >= n_bits) return -1;
            int w = from >> 6;
            uint64_t mask = ~0ULL << (from & 63);
            uint64_t bits = words[w] & mask;
            if (bits) return (w << 6) | ctz64(bits);
            int nw = static_cast<int>(words.size());
            for (++w; w < nw; ++w)
                if (words[w]) return (w << 6) | ctz64(words[w]);
            return -1;
        }
        BitClass intersectWith(const BitClass& mask) const {
            BitClass r(n_bits);
            int nw = static_cast<int>(words.size());
            for (int w = 0; w < nw; ++w) {
                uint64_t bits = words[w] & mask.words[w];
                r.words[w] = bits;
                r.card += popcount64(bits);
            }
            return r;
        }
        BitClass subtractMask(const BitClass& mask) const {
            BitClass r(n_bits);
            int nw = static_cast<int>(words.size());
            for (int w = 0; w < nw; ++w) {
                uint64_t bits = words[w] & ~mask.words[w];
                r.words[w] = bits;
                r.card += popcount64(bits);
            }
            return r;
        }
    };

    static void mcSplitRecurse(
        const MolGraph& g1, const MolGraph& g2, const ChemOptions& C, bool induced,
        TimeBudget& tb,
        const std::vector<BitClass>& qAdjMasks,
        const std::vector<BitClass>& tAdjMasks,
        std::vector<BitClass>& qSets, std::vector<BitClass>& tSets, int numClasses,
        int* q2t, int* t2q, int curSize, int upperBound,
        int* bestQ2T, int& bestSize, int64_t& nodeCount,
        int n1, int n2, int depth, int maxDepth) {

        nodeCount++;
        if (nodeCount > MAX_NODE_LIMIT) return;
        if (tb.expired()) return;
        if (curSize > bestSize) { bestSize = curSize; std::memcpy(bestQ2T, q2t, n1 * sizeof(int)); }
        if (curSize + upperBound <= bestSize || depth >= maxDepth) return;

        if (curSize > 0) {
            bool canExtend = false;
            for (int c = 0; c < numClasses && !canExtend; ++c)
                if (!qSets[c].empty() && !tSets[c].empty()) canExtend = true;
            if (!canExtend) return;
        }

        // Select most constrained class
        int bestClass = -1, bestMin = INT_MAX, bestClassQC = 0, bestClassTC = 0;
        for (int c = 0; c < numClasses; ++c) {
            int qc = qSets[c].cardinality(), tc = tSets[c].cardinality();
            if (qc == 0 || tc == 0) continue;
            int m = std::min(qc, tc);
            if (m < bestMin) { bestMin = m; bestClass = c; bestClassQC = qc; bestClassTC = tc; }
        }
        if (bestClass == -1) return;

        // Bidirectional selection
        bool branchFromTarget = bestClassTC < bestClassQC;

        if (!branchFromTarget) {
            // Select qi with connectivity-aware ordering
            int qi = -1, qiBestScore = -1;
            for (int v = qSets[bestClass].nextSetBit(0); v >= 0; v = qSets[bestClass].nextSetBit(v + 1)) {
                bool conn = false;
                for (int nb : g1.neighbors[v]) if (q2t[nb] >= 0) { conn = true; break; }
                int score = (conn ? 1000 : 0) + (g1.ring[v] ? 100 : 0) + g1.degree[v] * 10;
                if (score > qiBestScore) { qiBestScore = score; qi = v; }
            }

            // RRSplit: try each orbit of target.
            // Orbit pruning is only safe when qi has no mapped neighbor --
            // connectivity context differentiates "equivalent" atoms.
            // Phase 2.3: augment orbit key with ring-system signature pair
            // so that atoms from different ring systems are not pruned together.
            bool rsActiveQ = (g1.ringSystemCount() > 0 || g2.ringSystemCount() > 0);
            uint64_t qiRsSig = rsActiveQ ? g1.ringSystemSig(g1.ringSystemOf(qi)) : 0;
            std::unordered_set<uint64_t> triedOrbitRS;
            bool qiHasMappedNb = false;
            for (int nb : g1.neighbors[qi]) if (q2t[nb] >= 0) { qiHasMappedNb = true; break; }
            for (int tj = tSets[bestClass].nextSetBit(0); tj >= 0; tj = tSets[bestClass].nextSetBit(tj + 1)) {
                nodeCount++;
                if (nodeCount > MAX_NODE_LIMIT || tb.expired()) return;
                if (!qiHasMappedNb) {
                    // Combine orbit ID with ring-system signature for a finer-grained key
                    uint64_t orbitKey = static_cast<uint64_t>(g2.orbit[tj]);
                    if (rsActiveQ) {
                        uint64_t tjRsSig = g2.ringSystemSig(g2.ringSystemOf(tj));
                        // Mix ring-system pair into orbit key
                        orbitKey ^= (qiRsSig * 0x9e3779b97f4a7c15ULL) ^ (tjRsSig * 0x517cc1b727220a95ULL);
                    }
                    if (!triedOrbitRS.insert(orbitKey).second) continue;
                }

                if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;
                if (!checkBondCompat(g1, g2, C, induced, qi, tj, q2t)) continue;

                q2t[qi] = tj; t2q[tj] = qi;
                refineAndRecurse(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, qSets, tSets, numClasses,
                                 q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, qi, tj);
                q2t[qi] = -1; t2q[tj] = -1;
            }

            // Skip qi branch
            mcSplitSkipVertex(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, qSets, tSets, numClasses,
                              q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth,
                              bestClass, qi, true);
        } else {
            int tj = -1, tjBestScore = -1;
            for (int v = tSets[bestClass].nextSetBit(0); v >= 0; v = tSets[bestClass].nextSetBit(v + 1)) {
                bool conn = false;
                for (int nb : g2.neighbors[v]) if (t2q[nb] >= 0) { conn = true; break; }
                int score = (conn ? 1000 : 0) + (g2.ring[v] ? 100 : 0) + g2.degree[v] * 10;
                if (score > tjBestScore) { tjBestScore = score; tj = v; }
            }

            // Orbit pruning guard: skip only when tj has no mapped neighbor.
            // Phase 2.3: augment orbit key with ring-system signature pair.
            bool rsActiveT = (g1.ringSystemCount() > 0 || g2.ringSystemCount() > 0);
            uint64_t tjRsSig = rsActiveT ? g2.ringSystemSig(g2.ringSystemOf(tj)) : 0;
            std::unordered_set<uint64_t> triedOrbitRS;
            bool tjHasMappedNb = false;
            for (int nb : g2.neighbors[tj]) if (t2q[nb] >= 0) { tjHasMappedNb = true; break; }
            for (int qi = qSets[bestClass].nextSetBit(0); qi >= 0; qi = qSets[bestClass].nextSetBit(qi + 1)) {
                nodeCount++;
                if (nodeCount > MAX_NODE_LIMIT || tb.expired()) return;
                if (!tjHasMappedNb) {
                    uint64_t orbitKey = static_cast<uint64_t>(g1.orbit[qi]);
                    if (rsActiveT) {
                        uint64_t qiRsSig2 = g1.ringSystemSig(g1.ringSystemOf(qi));
                        orbitKey ^= (qiRsSig2 * 0x9e3779b97f4a7c15ULL) ^ (tjRsSig * 0x517cc1b727220a95ULL);
                    }
                    if (!triedOrbitRS.insert(orbitKey).second) continue;
                }

                if (!atomsCompatFast(g1, qi, g2, tj, C)) continue;

                // Check bond compat from target side
                bool bondOk = true;
                for (int tk : g2.neighbors[tj]) {
                    if (t2q[tk] < 0) continue;
                    int qk = t2q[tk];
                    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                    if (qOrd != 0 && tOrd != 0) {
                        if (!bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { bondOk = false; break; }
                    } else if (induced && ((qOrd!=0) != (tOrd!=0))) { bondOk = false; break; }
                }
                if (!bondOk) continue;

                q2t[qi] = tj; t2q[tj] = qi;
                refineAndRecurse(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, qSets, tSets, numClasses,
                                 q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, qi, tj);
                q2t[qi] = -1; t2q[tj] = -1;
            }

            mcSplitSkipVertex(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, qSets, tSets, numClasses,
                              q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth,
                              bestClass, tj, false);
        }
    }

    static void refineAndRecurse(
        const MolGraph& g1, const MolGraph& g2, const ChemOptions& C, bool induced,
        TimeBudget& tb,
        const std::vector<BitClass>& qAdjMasks,
        const std::vector<BitClass>& tAdjMasks,
        std::vector<BitClass>& qSets, std::vector<BitClass>& tSets, int numClasses,
        int* q2t, int* t2q, int curSize,
        int* bestQ2T, int& bestSize, int64_t& nodeCount,
        int n1, int n2, int depth, int maxDepth, int qi, int tj) {
        std::vector<BitClass> newQ, newT;
        int newUB = 0;
        newQ.reserve(numClasses * 2);
        newT.reserve(numClasses * 2);

        for (int c = 0; c < numClasses; ++c) {
            if (qSets[c].empty() || tSets[c].empty()) continue;
            BitClass qAdj = qSets[c].intersectWith(qAdjMasks[qi]);
            qAdj.clear(qi);
            BitClass qNon = qSets[c].subtractMask(qAdjMasks[qi]);
            qNon.clear(qi);
            BitClass tAdj = tSets[c].intersectWith(tAdjMasks[tj]);
            tAdj.clear(tj);
            BitClass tNon = tSets[c].subtractMask(tAdjMasks[tj]);
            tNon.clear(tj);

            if (qAdj.cardinality() > 0 && tAdj.cardinality() > 0) {
                newUB += std::min(qAdj.cardinality(), tAdj.cardinality());
                newQ.push_back(std::move(qAdj));
                newT.push_back(std::move(tAdj));
            }
            if (qNon.cardinality() > 0 && tNon.cardinality() > 0) {
                newUB += std::min(qNon.cardinality(), tNon.cardinality());
                newQ.push_back(std::move(qNon));
                newT.push_back(std::move(tNon));
            }
        }

        if (curSize + 1 + newUB > bestSize) {
            int newNum = static_cast<int>(newQ.size());
            mcSplitRecurse(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, newQ, newT, newNum,
                           q2t, t2q, curSize + 1, newUB,
                           bestQ2T, bestSize, nodeCount, n1, n2, depth + 1, maxDepth);
        }
    }

    static void mcSplitSkipVertex(
        const MolGraph& g1, const MolGraph& g2, const ChemOptions& C, bool induced,
        TimeBudget& tb,
        const std::vector<BitClass>& qAdjMasks,
        const std::vector<BitClass>& tAdjMasks,
        std::vector<BitClass>& qSets, std::vector<BitClass>& tSets, int numClasses,
        int* q2t, int* t2q, int curSize,
        int* bestQ2T, int& bestSize, int64_t& nodeCount,
        int n1, int n2, int depth, int maxDepth,
        int bestClass, int vertex, bool isQuery) {

        std::vector<BitClass> skipQ, skipT;
        int skipUB = 0;
        for (int c = 0; c < numClasses; ++c) {
            if (qSets[c].empty() || tSets[c].empty()) continue;
            if (c == bestClass) {
                BitClass reduced = isQuery ? qSets[c] : qSets[c];
                BitClass reducedT = isQuery ? tSets[c] : tSets[c];
                // Copy and remove vertex
                if (isQuery) { reduced = qSets[c]; reduced.clear(vertex); }
                else { reducedT = tSets[c]; reducedT.clear(vertex); }
                BitClass& checkSet = isQuery ? reduced : reducedT;
                if (!checkSet.empty()) {
                    int qCard = isQuery ? reduced.cardinality() : qSets[c].cardinality();
                    int tCard = isQuery ? tSets[c].cardinality() : reducedT.cardinality();
                    skipUB += std::min(qCard, tCard);
                    skipQ.push_back(isQuery ? reduced : qSets[c]);
                    skipT.push_back(isQuery ? tSets[c] : reducedT);
                }
            } else {
                skipUB += std::min(qSets[c].cardinality(), tSets[c].cardinality());
                skipQ.push_back(qSets[c]);
                skipT.push_back(tSets[c]);
            }
        }
        if (curSize + skipUB > bestSize) {
            int skipNum = static_cast<int>(skipQ.size());
            mcSplitRecurse(g1, g2, C, induced, tb, qAdjMasks, tAdjMasks, skipQ, skipT, skipNum,
                           q2t, t2q, curSize, skipUB,
                           bestQ2T, bestSize, nodeCount, n1, n2, depth + 1, maxDepth);
        }
    }

public:
    /// Phase 2.4: mcSplitSeed now delegates to the flat variant to
    /// avoid code duplication.  Materializes std::map only at the boundary.
    std::map<int,int> mcSplitSeed(TimeBudget& tb, int64_t& nodeCountOut) {
        auto [flat, sz] = mcSplitSeedFlat(tb, nodeCountOut);
        return flatToMap(flat.data(), g1_.n);
    }

    /// Phase 2.4: Flat-array variant of mcSplitSeed.
    /// Returns (q2t flat array where -1=unmapped, mapping size).
    std::pair<std::vector<int>, int> mcSplitSeedFlat(TimeBudget& tb, int64_t& nodeCountOut) {
        int n1 = g1_.n, n2 = g2_.n;
        if (n1 == 0 || n2 == 0) return {std::vector<int>(std::max(n1, 1), -1), 0};

        bool useMorgan = n1 >= 30 && n2 >= 30;
        if (useMorgan) {
            std::unordered_set<int> qDistinct;
            for (int i = 0; i < n1; ++i) qDistinct.insert(g1_.morganRank[i]);
            if (static_cast<int>(qDistinct.size()) > n1 * 3 / 4) useMorgan = false;
        }

        std::unordered_map<int, std::vector<int>> qGroups, tGroups;
        for (int i = 0; i < n1; ++i)
            qGroups[useMorgan ? g1_.morganRank[i] : g1_.label[i]].push_back(i);
        for (int j = 0; j < n2; ++j)
            tGroups[useMorgan ? g2_.morganRank[j] : g2_.label[j]].push_back(j);

        std::vector<BitClass> initQSets, initTSets;
        for (auto& [qlabel, qList] : qGroups)
            for (auto& [tlabel, tList] : tGroups) {
                if (qList.empty() || tList.empty()) continue;
                if (!atomsCompatFast(g1_, qList[0], g2_, tList[0], C_)) continue;
                BitClass qbs(n1), tbs(n2);
                for (int v : qList) qbs.set(v);
                for (int v : tList) tbs.set(v);
                initQSets.push_back(std::move(qbs));
                initTSets.push_back(std::move(tbs));
            }
        if (initQSets.empty()) return {std::vector<int>(n1, -1), 0};

        int numClasses = static_cast<int>(initQSets.size());
        int initUB = 0;
        for (int c = 0; c < numClasses; ++c)
            initUB += std::min(initQSets[c].cardinality(), initTSets[c].cardinality());

        std::vector<BitClass> qAdjMasks, tAdjMasks;
        qAdjMasks.reserve(n1);
        tAdjMasks.reserve(n2);
        for (int i = 0; i < n1; ++i) {
            BitClass mask(n1);
            for (int nb : g1_.neighbors[i]) mask.set(nb);
            qAdjMasks.push_back(std::move(mask));
        }
        for (int j = 0; j < n2; ++j) {
            BitClass mask(n2);
            for (int nb : g2_.neighbors[j]) mask.set(nb);
            tAdjMasks.push_back(std::move(mask));
        }

        std::vector<int> q2t(n1, -1), t2q(n2, -1);
        std::vector<int> bestQ2T(n1, -1);
        int bestSize = 0;
        int64_t nodeCount = 0;

        mcSplitRecurse(g1_, g2_, C_, induced_, tb, qAdjMasks, tAdjMasks, initQSets, initTSets, numClasses,
                       q2t.data(), t2q.data(), 0, initUB,
                       bestQ2T.data(), bestSize, nodeCount, n1, n2, 0, std::min(n1, n2) + 1);

        nodeCountOut = nodeCount;
        return {std::move(bestQ2T), bestSize};
    }

    // -- Bron-Kerbosch + Tomita (Level 3) ----------------------------------
private:
    static void greedyClique(const std::vector<std::vector<uint64_t>>& adj, int N, int words,
                             std::vector<int>& clique) {
        int bestStart = 0, bestDeg = 0;
        for (int i = 0; i < N; ++i) {
            int deg = popcountN(adj[i].data(), words);
            if (deg > bestDeg) { bestDeg = deg; bestStart = i; }
        }
        clique.clear();
        clique.push_back(bestStart);
        std::vector<uint64_t> cand = adj[bestStart];
        while (true) {
            int bestV = -1, bestConn = -1;
            for (int w = 0; w < words; ++w) {
                uint64_t bits = cand[w];
                while (bits) {
                    int bit = ctz64(bits);
                    int v = (w << 6) | bit;
                    int conn = popcountAnd(adj[v].data(), cand.data(), words);
                    if (conn > bestConn) { bestConn = conn; bestV = v; }
                    bits &= bits - 1;
                }
            }
            if (bestV == -1) break;
            clique.push_back(bestV);
            for (int w = 0; w < words; ++w) cand[w] &= adj[bestV][w];
        }
    }

    static int colorBound(const uint64_t* P, const std::vector<std::vector<uint64_t>>& adj,
                           int words, int N) {
        std::vector<int> color(N, -1);
        int maxColor = 0;
        int ucWords = std::max(1, (N + 63) >> 6);
        std::vector<uint64_t> usedColors(ucWords);
        for (int w = 0; w < words; ++w) {
            uint64_t bits = P[w];
            while (bits) {
                int bit = ctz64(bits);
                int v = (w << 6) | bit;
                std::fill(usedColors.begin(), usedColors.end(), 0ULL);
                for (int nw = 0; nw < words; ++nw) {
                    uint64_t nbits = adj[v][nw] & P[nw];
                    while (nbits) {
                        int nbit = ctz64(nbits);
                        int u = (nw << 6) | nbit;
                        int c = color[u];
                        if (c >= 0) usedColors[c >> 6] |= 1ULL << (c & 63);
                        nbits &= nbits - 1;
                    }
                }
                int c = 0;
                for (int uw = 0; uw < ucWords; ++uw) {
                    uint64_t free = ~usedColors[uw];
                    if (free) { c = (uw << 6) | ctz64(free); break; }
                    c = (uw + 1) << 6;
                }
                color[v] = c;
                if (c > maxColor) maxColor = c;
                bits &= bits - 1;
            }
        }
        return maxColor + 1;
    }

    static std::vector<int> iterateBits(const uint64_t* words_ptr, int N, int words) {
        std::vector<int> result;
        for (int w = 0; w < words; ++w) {
            uint64_t bits = words_ptr[w];
            while (bits) {
                int bit = ctz64(bits);
                int v = (w << 6) | bit;
                if (v < N) result.push_back(v);
                bits &= bits - 1;
            }
        }
        return result;
    }

    void bronKerboschPivot(
        uint64_t* R, uint64_t* P, uint64_t* X,
        std::vector<std::vector<uint64_t>>& adj,
        std::vector<int>& currentBest, int N, int words,
        const int* equivClass, int numEquivClasses, TimeBudget& tb,
        std::vector<std::vector<uint64_t>>& rStack,
        std::vector<std::vector<uint64_t>>& pStack,
        std::vector<std::vector<uint64_t>>& xStack,
        int depth, const std::vector<PGNode>& nodes) {

        if (tb.expired()) return;
        int rSize = popcountN(R, words);
        int pSize = popcountN(P, words);
        int xSize = popcountN(X, words);
        if (pSize == 0 && xSize == 0) {
            if (rSize > static_cast<int>(currentBest.size()))
                currentBest = iterateBits(R, N, words);
            return;
        }
        if (rSize + pSize <= static_cast<int>(currentBest.size())) return;

        // Partition bound
        if (pSize > 2) {
            std::vector<uint8_t> qiSeen(g1_.n, 0), tjSeen(g2_.n, 0);
            for (int w = 0; w < words; ++w) {
                uint64_t bits = P[w];
                while (bits) {
                    int bit = ctz64(bits);
                    int v = (w << 6) | bit;
                    if (v < N) { qiSeen[nodes[v].qi] = true; tjSeen[nodes[v].tj] = true; }
                    bits &= bits - 1;
                }
            }
            int qiCard = 0, tjCard = 0;
            for (int i = 0; i < g1_.n; ++i) if (qiSeen[i]) qiCard++;
            for (int j = 0; j < g2_.n; ++j) if (tjSeen[j]) tjCard++;
            if (rSize + std::min(qiCard, tjCard) <= static_cast<int>(currentBest.size())) return;
        }

        if (pSize > 0 && rSize + colorBound(P, adj, words, N) <= static_cast<int>(currentBest.size())) return;

        // Choose pivot
        int pivot = -1, pivotConn = -1;
        for (int w = 0; w < words; ++w) {
            uint64_t bits = P[w] | X[w];
            while (bits) {
                int bit = ctz64(bits);
                int u = (w << 6) | bit;
                int conn = popcountAnd(P, adj[u].data(), words);
                if (conn > pivotConn) { pivotConn = conn; pivot = u; }
                bits &= bits - 1;
            }
        }

        std::vector<uint64_t> candidates(words);
        if (pivot >= 0)
            for (int w = 0; w < words; ++w) candidates[w] = P[w] & ~adj[pivot][w];
        else
            std::memcpy(candidates.data(), P, words * sizeof(uint64_t));

        // Use vector<bool> instead of unordered_set for tried equivalence
        // classes -- equiv classes are dense small integers [0, numEquivClasses).
        std::vector<bool> triedClasses(numEquivClasses, false);
        for (int w = 0; w < words; ++w) {
            uint64_t bits = candidates[w];
            while (bits) {
                if (tb.expired()) return;
                int bit = ctz64(bits);
                int v = (w << 6) | bit;
                bits &= bits - 1;
                int ec = equivClass[v];
                if (triedClasses[ec]) continue;
                triedClasses[ec] = true;

                uint64_t* newR; uint64_t* newP; uint64_t* newX;
                if (depth < static_cast<int>(rStack.size())) {
                    newR = rStack[depth].data();
                    newP = pStack[depth].data();
                    newX = xStack[depth].data();
                } else {
                    rStack.emplace_back(words, 0);
                    pStack.emplace_back(words, 0);
                    xStack.emplace_back(words, 0);
                    newR = rStack.back().data();
                    newP = pStack.back().data();
                    newX = xStack.back().data();
                }
                for (int cw = 0; cw < words; ++cw) {
                    newR[cw] = R[cw]; newP[cw] = P[cw] & adj[v][cw]; newX[cw] = X[cw] & adj[v][cw];
                }
                newR[v >> 6] |= 1ULL << (v & 63);
                bronKerboschPivot(newR, newP, newX, adj, currentBest, N, words,
                                  equivClass, numEquivClasses, tb, rStack, pStack, xStack, depth + 1, nodes);
                P[v >> 6] &= ~(1ULL << (v & 63));
                X[v >> 6] |= 1ULL << (v & 63);
            }
        }
    }

public:
    std::map<int,int> maximumCliqueSeed(TimeBudget& tb) {
        int n1 = g1_.n, n2 = g2_.n;
        std::vector<PGNode> nodes;
        for (int i = 0; i < n1; ++i)
            for (int j = 0; j < n2; ++j) {
                if (!atomsCompatFast(g1_, i, g2_, j, C_)) continue;
                if (!induced_ && g1_.degree[i] > g2_.degree[j]) continue;
                nodes.push_back({i, j});
            }
        int N = static_cast<int>(nodes.size());
        if (N == 0) return {};
        // Cap product graph size to prevent O(N^2) blowup on large molecules
        if (N > 2500) return {};

        int words = (N + 63) >> 6;
        std::vector<std::vector<uint64_t>> adj(N, std::vector<uint64_t>(words, 0));

        for (int u = 0; u < N; ++u) {
            auto& nu = nodes[u];
            for (int v = u + 1; v < N; ++v) {
                if (tb.expired()) break;
                auto& nv = nodes[v];
                if (nu.qi == nv.qi || nu.tj == nv.tj) continue;
                int qOrd = g1_.bondOrder(nu.qi, nv.qi), tOrd = g2_.bondOrder(nu.tj, nv.tj);
                bool ok;
                if (qOrd != 0 && tOrd != 0)
                    ok = bondsCompatible(g1_, nu.qi, nv.qi, g2_, nu.tj, nv.tj, C_);
                else if (induced_)
                    ok = (qOrd == 0 && tOrd == 0);
                else
                    continue;
                if (ok) {
                    adj[u][v >> 6] |= 1ULL << (v & 63);
                    adj[v][u >> 6] |= 1ULL << (u & 63);
                }
            }
        }

        // Greedy clique
        std::vector<int> bestClique;
        greedyClique(adj, N, words, bestClique);

        // K-core reduction
        {
            int minDeg = std::max(2, static_cast<int>(bestClique.size()));
            bool changed = true;
            while (changed && !tb.expired()) {
                changed = false;
                for (int u = 0; u < N; ++u) {
                    int deg = popcountN(adj[u].data(), words);
                    if (deg == 0 || deg >= minDeg) continue;
                    for (int w = 0; w < words; ++w) {
                        uint64_t bits = adj[u][w];
                        while (bits) {
                            int bit = ctz64(bits);
                            int v = (w << 6) | bit;
                            adj[v][u >> 6] &= ~(1ULL << (u & 63));
                            bits &= bits - 1;
                        }
                        adj[u][w] = 0;
                    }
                    changed = true;
                }
            }
        }

        // Equivalence classes (orbit + ring-system signature)
        // Phase 2.3: augment orbit-based equivalence with ring-system pair signatures
        // so that nodes from different ring-system pairings are never treated as equivalent.
        bool useRingSysSym = (g1_.ringSystemCount() > 0 || g2_.ringSystemCount() > 0);
        std::vector<int> equivClass(N);
        int numEquivClasses = 0;
        {
            struct SigKey {
                int64_t orbitSig;
                uint64_t rsSig1;
                uint64_t rsSig2;
                bool operator==(const SigKey& o) const {
                    return orbitSig == o.orbitSig && rsSig1 == o.rsSig1 && rsSig2 == o.rsSig2;
                }
            };
            struct SigKeyHash {
                size_t operator()(const SigKey& k) const {
                    size_t h = std::hash<int64_t>{}(k.orbitSig);
                    h ^= std::hash<uint64_t>{}(k.rsSig1) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    h ^= std::hash<uint64_t>{}(k.rsSig2) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
                    return h;
                }
            };
            std::unordered_map<SigKey, int, SigKeyHash> sig2class;
            int& nextClass = numEquivClasses;
            for (int i = 0; i < N; ++i) {
                auto& nd = nodes[i];
                int deg = popcountN(adj[i].data(), words);
                int64_t orbitSig = (static_cast<int64_t>(g1_.orbit[nd.qi])          << 44)
                    | (static_cast<int64_t>(g2_.orbit[nd.tj])                       << 28)
                    | (static_cast<int64_t>(g1_.morganRank[nd.qi] & 0xFF)           << 20)
                    | (static_cast<int64_t>(g1_.degree[nd.qi])                      << 10) | deg;
                uint64_t rs1 = 0, rs2 = 0;
                if (useRingSysSym) {
                    rs1 = g1_.ringSystemSig(g1_.ringSystemOf(nd.qi));
                    rs2 = g2_.ringSystemSig(g2_.ringSystemOf(nd.tj));
                }
                SigKey key{orbitSig, rs1, rs2};
                auto it = sig2class.find(key);
                if (it == sig2class.end()) { sig2class[key] = nextClass; equivClass[i] = nextClass++; }
                else equivClass[i] = it->second;
            }
        }

        // BK search
        std::vector<uint64_t> P(words, 0), X(words, 0), R(words, 0);
        for (int i = 0; i < N; ++i) {
            if (popcountN(adj[i].data(), words) > 0)
                P[i >> 6] |= 1ULL << (i & 63);
        }

        int maxDepth = std::min(N, 200);
        std::vector<std::vector<uint64_t>> rStack(maxDepth, std::vector<uint64_t>(words, 0));
        std::vector<std::vector<uint64_t>> pStack(maxDepth, std::vector<uint64_t>(words, 0));
        std::vector<std::vector<uint64_t>> xStack(maxDepth, std::vector<uint64_t>(words, 0));

        bronKerboschPivot(R.data(), P.data(), X.data(), adj, bestClique, N, words,
                          equivClass.data(), numEquivClasses, tb, rStack, pStack, xStack, 0, nodes);

        std::map<int,int> seed;
        std::unordered_set<int> usedQ, usedT;
        for (int u : bestClique) {
            auto& nd = nodes[u];
            if (!usedQ.count(nd.qi) && !usedT.count(nd.tj)) {
                seed[nd.qi] = nd.tj;
                usedQ.insert(nd.qi);
                usedT.insert(nd.tj);
            }
        }
        return seed;
    }

    // -- Ring anchor seed --------------------------------------------------
    std::map<int,int> ringAnchorSeed(TimeBudget& tb) {
        std::map<int,int> seed;
        std::vector<int> R1, R2;
        for (int i = 0; i < g1_.n; ++i) if (g1_.ring[i]) R1.push_back(i);
        for (int j = 0; j < g2_.n; ++j) if (g2_.ring[j]) R2.push_back(j);
        std::vector<uint8_t> used(g2_.n, 0);
        for (int i : R1) {
            int bestJ = -1, bestScore = INT_MIN;
            for (int j : R2) {
                if (used[j]) continue;
                if (!atomsCompatFast(g1_, i, g2_, j, C_)) continue;
                int score = (g1_.aromatic[i] && g2_.aromatic[j] ? 10 : 0) + std::min(g1_.degree[i], g2_.degree[j]);
                if (score > bestScore) { bestScore = score; bestJ = j; }
            }
            if (bestJ >= 0) { seed[i] = bestJ; used[bestJ] = 1; }
            if (tb.expired()) break;
        }
        return seed;
    }

    // -- Label-degree anchor seed ------------------------------------------
    std::map<int,int> labelDegreeAnchorSeed(TimeBudget& tb) {
        std::map<int,int> seed;
        std::unordered_map<int,int> freq;
        for (int i = 0; i < g1_.n; ++i) freq[g1_.label[i]]++;
        std::vector<int> Q(g1_.n);
        std::iota(Q.begin(), Q.end(), 0);
        std::sort(Q.begin(), Q.end(), [&](int a, int b) {
            int ra = freq[g1_.label[a]], rb = freq[g1_.label[b]];
            return ra != rb ? ra < rb : g1_.degree[a] > g1_.degree[b];
        });
        std::vector<uint8_t> used(g2_.n, 0);
        for (int qi : Q) {
            int bestJ = -1, bestScore = INT_MIN;
            for (int tj = 0; tj < g2_.n; ++tj) {
                if (used[tj]) continue;
                if (!atomsCompatFast(g1_, qi, g2_, tj, C_)) continue;
                int score = 10 - std::abs(g1_.degree[qi] - g2_.degree[tj]);
                if (score > bestScore) { bestScore = score; bestJ = tj; }
            }
            if (bestJ >= 0) { seed[qi] = bestJ; used[bestJ] = 1; }
            if (tb.expired()) break;
        }
        return seed;
    }
    // -- VF2++ neighborhood seed (shared implementation) --------------------

    std::map<int,int> vf2ppNeighborhoodSeed(
            TimeBudget& tb, int64_t millis, int /*radius*/, int maxAnchors,
            const std::vector<uint8_t>& qCand, const std::vector<uint8_t>& tCand) {
        struct AP { int qi, tj, score; };
        std::vector<AP> pairs;
        for (int i = 0; i < g1_.n; ++i) {
            if (!qCand[i]) continue;
            for (int j = 0; j < g2_.n; ++j) {
                if (!tCand[j]) continue;
                if (!atomsCompatFast(g1_, i, g2_, j, C_)) continue;
                int score = 100 - std::abs(g1_.degree[i] - g2_.degree[j]);
                if (g1_.aromatic[i] && g2_.aromatic[j]) score += 10;
                pairs.push_back({i, j, score});
            }
        }
        if (pairs.empty()) return {};
        std::sort(pairs.begin(), pairs.end(),
                  [](const AP& a, const AP& b) { return a.score > b.score; });
        if (static_cast<int>(pairs.size()) > maxAnchors)
            pairs.resize(maxAnchors);

        // Run VF2PP once with the full budget — anchor pairs are used for
        // scoring/ordering to select the best candidates, but the matcher
        // operates on the full graphs.  Running once avoids N redundant searches.
        if (tb.expired() || millis <= 0) return {};
        VF2PPMatcher m(g1_, g2_, C_, millis);
        std::vector<std::vector<std::pair<int,int>>> out;
        m.enumerate(1, out);
        if (out.empty()) return {};
        std::map<int,int> best;
        for (auto& [qi, tj] : out[0]) best[qi] = tj;
        return best;
    }

    // -- VF2++ ring skeleton seed -------------------------------------------

    std::map<int,int> vf2ppRingSkeletonSeed(TimeBudget& tb, int64_t millis,
                                             int radius, int maxAnchors) {
        std::vector<uint8_t> qRing(g1_.n, 0), tRing(g2_.n, 0);
        bool qEmpty = true, tEmpty = true;
        for (int i = 0; i < g1_.n; ++i) if (g1_.ring[i]) { qRing[i] = true; qEmpty = false; }
        for (int j = 0; j < g2_.n; ++j) if (g2_.ring[j]) { tRing[j] = true; tEmpty = false; }
        if (qEmpty || tEmpty || millis <= 0) return {};
        return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qRing, tRing);
    }

    // -- VF2++ heavy-atom core seed -----------------------------------------

    std::map<int,int> vf2ppCoreSeed(TimeBudget& tb, int64_t millis,
                                     int radius, int maxAnchors) {
        // Iteratively prune non-ring, non-aromatic leaf carbons
        std::vector<uint8_t> qMask(g1_.n, 1), tMask(g2_.n, 1);
        auto prune = [](const MolGraph& g, std::vector<uint8_t>& mask) {
            bool changed;
            do {
                changed = false;
                for (int i = 0; i < g.n; ++i) {
                    if (!mask[i] || g.ring[i] || g.aromatic[i]) continue;
                    int deg = 0;
                    for (int nb : g.neighbors[i]) if (mask[nb]) deg++;
                    if (deg <= 1 && g.atomicNum[i] == 6) { mask[i] = false; changed = true; }
                }
            } while (changed);
            bool empty = true;
            for (int i = 0; i < g.n; ++i) if (mask[i]) { empty = false; break; }
            if (empty) for (int i = 0; i < g.n; ++i) if (g.atomicNum[i] != 1) mask[i] = true;
        };
        prune(g1_, qMask);
        prune(g2_, tMask);
        bool qEmpty = true, tEmpty = true;
        for (int i = 0; i < g1_.n; ++i) if (qMask[i]) { qEmpty = false; break; }
        for (int j = 0; j < g2_.n; ++j) if (tMask[j]) { tEmpty = false; break; }
        if (qEmpty || tEmpty || millis <= 0) return {};
        return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qMask, tMask);
    }
}; // class GraphBuilder

} // namespace detail

// ===========================================================================
/// Count heteroatoms (non-C, non-H) in an MCS mapping. Used as a secondary
/// criterion: when two mappings have equal size, prefer the one with more
/// heteroatoms (N, O, S, P, etc.) for better chemical relevance.
/// This ensures reaction-critical atoms (S in SAM, N in amino acids) are
/// preferentially included in the MCS mapping.
inline int heteroatomScore(const MolGraph& g, const std::map<int,int>& mapping) {
    int score = 0;
    for (const auto& [k, v] : mapping) {
        int z = g.atomicNum[k];
        if (z != 6 && z != 1) score++; // not carbon, not hydrogen
    }
    return score;
}

/// Check if candidate is a better MCS than current best.
/// Primary: size. Secondary: heteroatom count (chemistry-first tiebreaker).
inline bool isBetterMCS(const MolGraph& g1,
                         const std::map<int,int>& candidate, int candSize,
                         const std::map<int,int>& best, int bestSize, int bestHetero) {
    if (candSize > bestSize) return true;
    if (candSize == bestSize && candSize > 0) {
        int candHetero = heteroatomScore(g1, candidate);
        return candHetero > bestHetero;
    }
    return false;
}

inline std::map<int,int> repairInvalidMcsMapping(const MolGraph& g1, const MolGraph& g2,
                                                 std::map<int,int> mapping,
                                                 const ChemOptions& opts);
inline std::vector<std::string> validateMapping(const MolGraph& g1, const MolGraph& g2,
                                                const std::map<int,int>& mapping,
                                                const ChemOptions& opts);

namespace detail {

inline int64_t resolveMcsTimeoutMs(const MolGraph& g1, const MolGraph& g2,
                                   const MCSOptions& opts) {
    if (opts.timeoutMs >= 0) return opts.timeoutMs;
    return std::min<int64_t>(30000, int64_t(500) + int64_t(g1.n) * g2.n * 2);
}

inline int graphConstraintScore(const MolGraph& g) {
    int hetero = 0, ring = 0, aromatic = 0, maxDeg = 0;
    for (int i = 0; i < g.n; ++i) {
        if (g.atomicNum[i] != 1 && g.atomicNum[i] != 6) hetero++;
        if (g.ring[i]) ring++;
        if (g.aromatic[i]) aromatic++;
        if (g.degree[i] > maxDeg) maxDeg = g.degree[i];
    }
    return hetero * 64 + ring * 16 + aromatic * 8 + maxDeg * 4;
}

inline int graphOrbitRedundancy(const MolGraph& g) {
    std::unordered_set<int> uniqueOrbits;
    uniqueOrbits.reserve(static_cast<size_t>(g.n) * 2 + 1);
    for (int orbitId : g.orbit) uniqueOrbits.insert(orbitId);
    return g.n - static_cast<int>(uniqueOrbits.size());
}

inline int graphOrientationBurden(const MolGraph& g) {
    int degreeSum = 0;
    int ringAtoms = 0;
    int aromaticAtoms = 0;
    for (int i = 0; i < g.n; ++i) {
        degreeSum += g.degree[i];
        ringAtoms += g.ring[i] ? 1 : 0;
        aromaticAtoms += g.aromatic[i] ? 1 : 0;
    }
    int edgeCount = degreeSum / 2;
    return edgeCount * 8
        + ringAtoms * 16
        + aromaticAtoms * 4
        + graphOrbitRedundancy(g) * 24;
}

inline int64_t orientationConstraintMass(const MolGraph& query, const MolGraph& target,
                                         const ChemOptions& chem) {
    std::unordered_map<int, int> targetFreq;
    targetFreq.reserve(static_cast<size_t>(target.n) * 2 + 1);
    for (int j = 0; j < target.n; ++j) targetFreq[ubBaseLabel(target, j, chem)]++;

    int64_t mass = 0;
    for (int i = 0; i < query.n; ++i) {
        auto it = targetFreq.find(ubBaseLabel(query, i, chem));
        if (it != targetFreq.end()) mass += it->second;
    }
    return mass;
}

inline int orientationSeedProbeSize(const MolGraph& query, const MolGraph& target,
                                    const ChemOptions& chem) {
    if (query.n == 0 || target.n == 0) return 0;
    GraphBuilder gb(query, target, chem, false);
    TimeBudget tb(24LL * 60 * 60 * 1000);
    int upperBound = labelFrequencyUpperBound(query, target, chem);
    // Phase 2.4: use flat variant -- only size is needed, no map construction.
    return gb.seedExtendMCSFlat(tb, upperBound).second;
}

inline int orientationMcSplitProbeSize(const MolGraph& query, const MolGraph& target,
                                       const ChemOptions& chem) {
    if (query.n == 0 || target.n == 0) return 0;
    GraphBuilder gb(query, target, chem, false);
    TimeBudget tb(24LL * 60 * 60 * 1000);
    int64_t nodeCount = 0;
    // Phase 2.4: use flat variant -- only size is needed, no map construction.
    return gb.mcSplitSeedFlat(tb, nodeCount).second;
}

struct OrientationPlan {
    bool directFirst = true;
    int seed12 = 0;
    int seed21 = 0;
    int mc12 = 0;
    int mc21 = 0;
};

inline OrientationPlan chooseOrientationPlan(const MolGraph& g1, const MolGraph& g2,
                                             const ChemOptions& chem,
                                             int ub12, int ub21) {
    OrientationPlan plan;
    // Deterministic funnel:
    // 1. Prefer the more selective query (lower compatibility mass)
    // 2. Prefer the lower structural burden (smaller/sparser/less symmetric)
    // 3. Prefer the graph with richer chemical constraints
    // 4. For materially different graph sizes, search smaller-query-first
    // 5. Only then probe both directions with fixed node-limited seed stages
    // 6. Use the directed upper bound only as a late quality-preserving tie-break
    // 7. Finally stabilize on canonical hash
    int64_t mass12 = orientationConstraintMass(g1, g2, chem);
    int64_t mass21 = orientationConstraintMass(g2, g1, chem);
    if (mass12 != mass21) {
        plan.directFirst = mass12 < mass21;
        return plan;
    }

    int minN = std::min(g1.n, g2.n);
    int product = g1.n * g2.n;
    bool nearSymmetric = std::abs(g1.n - g2.n) <= 2;
    if (nearSymmetric && minN >= 12 && product >= 256) {
        plan.seed12 = orientationSeedProbeSize(g1, g2, chem);
        plan.seed21 = orientationSeedProbeSize(g2, g1, chem);
        if (plan.seed12 != plan.seed21) {
            plan.directFirst = plan.seed12 > plan.seed21;
            return plan;
        }
    }
    if (nearSymmetric && minN >= 12 && product >= 512) {
        plan.mc12 = orientationMcSplitProbeSize(g1, g2, chem);
        plan.mc21 = orientationMcSplitProbeSize(g2, g1, chem);
        if (plan.mc12 != plan.mc21) {
            plan.directFirst = plan.mc12 > plan.mc21;
            return plan;
        }
    }

    int burden1 = graphOrientationBurden(g1);
    int burden2 = graphOrientationBurden(g2);
    if (burden1 != burden2) {
        plan.directFirst = burden1 < burden2;
        return plan;
    }

    int score1 = graphConstraintScore(g1);
    int score2 = graphConstraintScore(g2);
    if (score1 != score2) {
        plan.directFirst = score1 > score2;
        return plan;
    }

    if (std::abs(g1.n - g2.n) >= 4) {
        plan.directFirst = g1.n < g2.n;
        return plan;
    }

    if (!nearSymmetric && minN >= 12 && product >= 256) {
        plan.seed12 = orientationSeedProbeSize(g1, g2, chem);
        plan.seed21 = orientationSeedProbeSize(g2, g1, chem);
        if (plan.seed12 != plan.seed21) {
            plan.directFirst = plan.seed12 > plan.seed21;
            return plan;
        }
    }
    if (!nearSymmetric && minN >= 12 && product >= 512) {
        plan.mc12 = orientationMcSplitProbeSize(g1, g2, chem);
        plan.mc21 = orientationMcSplitProbeSize(g2, g1, chem);
        if (plan.mc12 != plan.mc21) {
            plan.directFirst = plan.mc12 > plan.mc21;
            return plan;
        }
    }

    if (g1.n != g2.n) {
        plan.directFirst = g1.n < g2.n;
        return plan;
    }
    if (ub12 != ub21) {
        plan.directFirst = ub12 > ub21;
        return plan;
    }
    plan.directFirst = g1.canonicalHash <= g2.canonicalHash;
    return plan;
}

inline int alternateOrientationFloor(const OrientationPlan& plan, bool directRan) {
    return directRan
        ? std::max(plan.seed21, plan.mc21)
        : std::max(plan.seed12, plan.mc12);
}

inline bool preferFinalMapping(const MolGraph& g1,
                               const std::map<int,int>& candidate,
                               const std::map<int,int>& incumbent,
                               const MCSOptions& opts) {
    if (candidate.empty()) return false;
    if (incumbent.empty()) return true;

    int candScore = mcsScore(g1, candidate, opts);
    int bestScore = mcsScore(g1, incumbent, opts);
    bool weightMode = opts.maximizeBonds || !opts.atomWeights.empty();
    if (weightMode && candScore != bestScore) return candScore > bestScore;

    if (candidate.size() != incumbent.size()) return candidate.size() > incumbent.size();

    if (!weightMode && candScore != bestScore) return candScore > bestScore;

    int candHetero = heteroatomScore(g1, candidate);
    int bestHetero = heteroatomScore(g1, incumbent);
    if (candHetero != bestHetero) return candHetero > bestHetero;

    int candBonds = countMappedBonds(g1, candidate);
    int bestBonds = countMappedBonds(g1, incumbent);
    if (candBonds != bestBonds) return candBonds > bestBonds;

    return candidate < incumbent;
}

} // namespace detail

// Public API
// ===========================================================================

/// Find the Maximum Common Substructure between two molecular graphs.
/// Returns a mapping from g1 atom indices to g2 atom indices.
inline std::map<int,int> findMCSImpl(const MolGraph& g1, const MolGraph& g2,
                                     const ChemOptions& chem, const MCSOptions& opts) {
    using namespace detail;

    // Defensive: empty graphs have no common substructure
    if (g1.n == 0 || g2.n == 0) return {};

    // Validate atomWeights size if provided
    if (!opts.atomWeights.empty() && static_cast<int>(opts.atomWeights.size()) < g1.n)
        throw std::invalid_argument(
            "atomWeights size (" + std::to_string(opts.atomWeights.size()) +
            ") < query atom count (" + std::to_string(g1.n) + ")");

    static constexpr int GREEDY_PROBE_MAX_SIZE = 40;
    static constexpr int SEED_EXTEND_MAX_ATOMS = 50;
    static constexpr double BK_SKIP_RATIO = 0.8;

    // Resolve adaptive timeout: scale with product of atom counts, capped at 30s
    int64_t timeout = opts.timeoutMs;
    if (timeout < 0) {
        timeout = std::min(int64_t(30000), int64_t(500) + int64_t(g1.n) * g2.n * 2);
    }
    TimeBudget tb(timeout);

    // Exact identity fast-path before canonicalization: same object or the
    // same graph in the same index order should never pay the full MCS cost.
    if (&g1 == &g2 || detail::isExactMatch(g1, g2, chem)) {
        std::map<int,int> id;
        for (int i = 0; i < g1.n; i++) id[i] = i;
        if (opts.disconnectedMCS && (opts.minFragmentSize > 1 || opts.maxFragments < INT_MAX))
            id = detail::applyFragmentConstraints(g1, id, opts.minFragmentSize, opts.maxFragments);
        return id;
    }

    // Ensure canonical labeling / Morgan ranks are available (lazy init)
    g1.ensureCanonical();
    g2.ensureCanonical();
    if (chem.ringFusionMode != ChemOptions::RingFusionMode::IGNORE) {
        g1.ensureRingCounts();
        g2.ensureRingCounts();
    }

    // DSB is only admissible for induced MCS (degree constraint valid).
    // For non-induced (default), use the looser but safe label-frequency bound.
    int upperBound = opts.induced
        ? degreeSequenceUpperBound(g1, g2, chem)
        : labelFrequencyUpperBound(g1, g2, chem);
    if (!opts.induced) {
        upperBound = std::min(upperBound, detail::labelFrequencyUpperBoundDirected(g1, g2, chem));
        upperBound = std::min(upperBound, detail::labelFrequencyUpperBoundDirected(g2, g1, chem));
    }
    std::map<int,int> best;
    int bestSize = 0, bestScore = 0;

    // Canonical-equivalence fast-path for same graphs with permuted indices.
    if ((!chem.useChirality && !chem.useBondStereo)
        && detail::sameCanonicalGraph(g1, g2)) {
        std::map<int,int> id;
        std::vector<int> head(g2.n, -1), nxt(g2.n, -1);
        for (int j = g2.n - 1; j >= 0; j--) { nxt[j] = head[g2.canonicalLabel[j]]; head[g2.canonicalLabel[j]] = j; }
        for (int i = 0; i < g1.n; i++) { int cl = g1.canonicalLabel[i]; id[i] = head[cl]; head[cl] = nxt[head[cl]]; }
        if (opts.disconnectedMCS && (opts.minFragmentSize > 1 || opts.maxFragments < INT_MAX))
            id = detail::applyFragmentConstraints(g1, id, opts.minFragmentSize, opts.maxFragments);
        return id;
    }

    int minN = std::min(g1.n, g2.n);
    int bestHetero = 0; // heteroatom tiebreaker for equal-size MCS
    bool weightMode = opts.maximizeBonds || !opts.atomWeights.empty();

    // --- Level 0.25: Linear chain fast-path ---
    // For chain-like molecules (max degree ≤ 2, e.g., PEG), the general MCS
    // algorithms suffer from O(n!) symmetry explosion because all internal
    // atoms are identical.  Detect this case and use O(n*m) longest-common-
    // subpath DP instead.
    if (!weightMode) {
        bool g1Chain = true, g2Chain = true;
        for (int i = 0; i < g1.n && g1Chain; ++i) if (g1.degree[i] > 2) g1Chain = false;
        for (int j = 0; j < g2.n && g2Chain; ++j) if (g2.degree[j] > 2) g2Chain = false;
        if (g1Chain && g2Chain && g1.n >= 2 && g2.n >= 2) {
            // Extract chain atom sequences by walking from an endpoint
            auto walkChain = [](const MolGraph& g) -> std::vector<int> {
                int start = -1;
                for (int i = 0; i < g.n; ++i)
                    if (g.degree[i] <= 1) { start = i; break; }
                if (start < 0) start = 0; // ring — pick any
                std::vector<int> seq;
                std::vector<uint8_t> visited(g.n, 0);
                int cur = start;
                while (cur >= 0) {
                    visited[cur] = true;
                    seq.push_back(cur);
                    int next = -1;
                    for (int nb : g.neighbors[cur])
                        if (!visited[nb]) { next = nb; break; }
                    cur = next;
                }
                return seq;
            };
            auto seq1 = walkChain(g1);
            auto seq2 = walkChain(g2);
            int n1 = static_cast<int>(seq1.size());
            int n2 = static_cast<int>(seq2.size());

            // DP: longest common subpath matching atom labels + bond orders
            // dp[i][j] = length of longest common suffix ending at seq1[i], seq2[j]
            std::vector<std::vector<int>> dp(n1 + 1, std::vector<int>(n2 + 1, 0));
            int bestLen = 0, bestI = 0, bestJ = 0;
            for (int i = 1; i <= n1; ++i) {
                int a1 = seq1[i - 1];
                for (int j = 1; j <= n2; ++j) {
                    int a2 = seq2[j - 1];
                    if (!atomsCompatFast(g1, a1, g2, a2, chem)) continue;
                    if (i > 1 && j > 1) {
                        int b1prev = seq1[i - 2], b2prev = seq2[j - 2];
                        if (dp[i - 1][j - 1] > 0 &&
                            bondsCompatible(g1, b1prev, a1, g2, b2prev, a2, chem)) {
                            dp[i][j] = dp[i - 1][j - 1] + 1;
                        } else {
                            dp[i][j] = 1; // start new subpath
                        }
                    } else {
                        dp[i][j] = 1;
                    }
                    if (dp[i][j] > bestLen) { bestLen = dp[i][j]; bestI = i; bestJ = j; }
                }
            }
            if (bestLen > 0) {
                std::map<int, int> chainMCS;
                for (int k = 0; k < bestLen; ++k)
                    chainMCS[seq1[bestI - bestLen + k]] = seq2[bestJ - bestLen + k];
                int chainScore = mcsScore(g1, chainMCS, opts);
                if (!weightMode ? bestLen >= bestSize : chainScore > bestScore) {
                    best = chainMCS; bestSize = bestLen; bestScore = chainScore;
                }
                if (bestLen >= upperBound) return best;
            }
        }
    }

    // --- Level 0.5: Tree fast-path (branched polymers, dendrimers, glycogen) ---
    // For acyclic connected molecules (trees) that have degree > 2, the general MCS
    // algorithms are expensive.  Detect tree topology (edges == atoms - 1) and use
    // Kilpelainen-Mannila style bottom-up DP on rooted trees: O(n1*n2*d^2).
    if (!weightMode && g1.n >= 10 && g2.n >= 10) {
        auto isTree = [](const MolGraph& g) -> bool {
            int edgeCount = 0;
            for (int i = 0; i < g.n; ++i)
                edgeCount += static_cast<int>(g.neighbors[i].size());
            edgeCount /= 2;
            if (edgeCount != g.n - 1) return false;
            std::vector<uint8_t> vis(g.n, 0);
            std::deque<int> bfs;
            bfs.push_back(0); vis[0] = 1;
            int cnt = 1;
            while (!bfs.empty()) {
                int u = bfs.front(); bfs.pop_front();
                for (int v : g.neighbors[u])
                    if (!vis[v]) { vis[v] = true; ++cnt; bfs.push_back(v); }
            }
            return cnt == g.n;
        };

        if (isTree(g1) && isTree(g2)) {
            // Find centroid by iteratively peeling leaves
            auto findCentroid = [](const MolGraph& g) -> int {
                if (g.n <= 2) return 0; // trivial tree — any node is centroid
                std::vector<int> deg(g.n);
                std::deque<int> leaves;
                for (int i = 0; i < g.n; ++i) {
                    deg[i] = static_cast<int>(g.neighbors[i].size());
                    if (deg[i] <= 1) leaves.push_back(i);
                }
                int remaining = g.n;
                while (remaining > 2) {
                    int sz = static_cast<int>(leaves.size());
                    remaining -= sz;
                    std::deque<int> newLeaves;
                    for (int i = 0; i < sz; ++i) {
                        int u = leaves[i];
                        for (int v : g.neighbors[u])
                            if (--deg[v] == 1) newLeaves.push_back(v);
                    }
                    leaves = std::move(newLeaves);
                }
                return leaves.front();
            };

            // Root tree at centroid, build children + postorder
            struct RootedTree {
                std::vector<int> parent;
                std::vector<std::vector<int>> children;
                std::vector<int> postorder;
                int root;
            };
            auto rootTree = [](const MolGraph& g, int root) -> RootedTree {
                RootedTree rt;
                rt.root = root;
                rt.parent.assign(g.n, -1);
                rt.children.resize(g.n);
                rt.postorder.reserve(g.n);
                std::vector<uint8_t> vis(g.n, 0);
                std::deque<int> bfs;
                bfs.push_back(root); vis[root] = 1;
                std::vector<int> bfsOrder;
                bfsOrder.reserve(g.n);
                while (!bfs.empty()) {
                    int u = bfs.front(); bfs.pop_front();
                    bfsOrder.push_back(u);
                    for (int v : g.neighbors[u]) {
                        if (!vis[v]) {
                            vis[v] = true;
                            rt.parent[v] = u;
                            rt.children[u].push_back(v);
                            bfs.push_back(v);
                        }
                    }
                }
                rt.postorder.assign(bfsOrder.rbegin(), bfsOrder.rend());
                return rt;
            };

            int root1 = findCentroid(g1);
            // A tree may have 2 centroids. Try both rootings of g2
            // and take the maximum to avoid orientation mismatch.
            auto findCentroids = [](const MolGraph& g) -> std::vector<int> {
                if (g.n <= 2) return {0};
                std::vector<int> deg(g.n);
                std::deque<int> leaves;
                for (int i = 0; i < g.n; ++i) {
                    deg[i] = static_cast<int>(g.neighbors[i].size());
                    if (deg[i] <= 1) leaves.push_back(i);
                }
                int remaining = g.n;
                while (remaining > 2) {
                    int sz = static_cast<int>(leaves.size());
                    remaining -= sz;
                    std::deque<int> nl;
                    for (int i = 0; i < sz; ++i) {
                        int u = leaves[i];
                        for (int v : g.neighbors[u])
                            if (--deg[v] == 1) nl.push_back(v);
                    }
                    leaves = std::move(nl);
                }
                return std::vector<int>(leaves.begin(), leaves.end());
            };
            auto g2Centroids = findCentroids(g2);
            int root2 = g2Centroids[0];
            auto rt1 = rootTree(g1, root1);
            auto rt2 = rootTree(g2, root2);

            int n1t = g1.n, n2t = g2.n;
            // dp[u][v] = max common subtree size rooted at (u in T1, v in T2)
            std::vector<std::vector<int>> dpTree(n1t, std::vector<int>(n2t, 0));
            // Backtrack: child matching indices for reconstruction
            std::vector<std::vector<std::vector<std::pair<int,int>>>> btTree(
                n1t, std::vector<std::vector<std::pair<int,int>>>(n2t));

            for (int u : rt1.postorder) {
                for (int v : rt2.postorder) {
                    if (!atomsCompatFast(g1, u, g2, v, chem)) {
                        dpTree[u][v] = 0;
                        continue;
                    }
                    const auto& ch1 = rt1.children[u];
                    const auto& ch2 = rt2.children[v];
                    if (ch1.empty() || ch2.empty()) {
                        dpTree[u][v] = 1;
                        continue;
                    }
                    // Greedy bipartite matching by descending dp value
                    int d1 = static_cast<int>(ch1.size());
                    int d2 = static_cast<int>(ch2.size());
                    struct CandEdge { int weight; int i; int j; };
                    std::vector<CandEdge> edges;
                    edges.reserve(d1 * d2);
                    for (int i2 = 0; i2 < d1; ++i2) {
                        for (int j2 = 0; j2 < d2; ++j2) {
                            int w = dpTree[ch1[i2]][ch2[j2]];
                            if (w > 0 && bondsCompatible(g1, u, ch1[i2], g2, v, ch2[j2], chem))
                                edges.push_back({w, i2, j2});
                        }
                    }
                    std::sort(edges.begin(), edges.end(),
                              [](const CandEdge& a, const CandEdge& b) { return a.weight > b.weight; });
                    std::vector<uint8_t> usedI(d1, 0), usedJ(d2, 0);
                    int childSum = 0;
                    std::vector<std::pair<int,int>> matching;
                    for (auto& e : edges) {
                        if (usedI[e.i] || usedJ[e.j]) continue;
                        usedI[e.i] = true; usedJ[e.j] = true;
                        childSum += e.weight;
                        matching.push_back({e.i, e.j});
                    }
                    dpTree[u][v] = 1 + childSum;
                    btTree[u][v] = std::move(matching);
                }
            }

            // Find best root pair
            int treeBest = 0, bestU = -1, bestV = -1;
            for (int u = 0; u < n1t; ++u)
                for (int v = 0; v < n2t; ++v)
                    if (dpTree[u][v] > treeBest) { treeBest = dpTree[u][v]; bestU = u; bestV = v; }

            if (treeBest > bestSize) {
                std::map<int,int> treeMap;
                std::function<void(int, int)> reconstruct = [&](int u, int v) {
                    treeMap[u] = v;
                    for (auto& [ci, cj] : btTree[u][v])
                        reconstruct(rt1.children[u][ci], rt2.children[v][cj]);
                };
                reconstruct(bestU, bestV);

                int treeScore = mcsScore(g1, treeMap, opts);
                if (static_cast<int>(treeMap.size()) > bestSize) {
                    best = std::move(treeMap);
                    bestSize = static_cast<int>(best.size()); bestHetero = heteroatomScore(g1, best);
                    bestScore = treeScore;
                }
                if (bestSize >= upperBound) return best;
            }
        }
    }

    // --- Level 0.75: Greedy probe ---
    if (!weightMode && minN > 2 && minN < GREEDY_PROBE_MAX_SIZE && upperBound >= minN) {
        auto greedy = greedyProbe(g1, g2, chem, opts.templateFuzzyAtoms);
        // Fuzzy template: accept if mapping covers target within tolerance
        int fuzzy = opts.templateFuzzyAtoms;
        int greedySz = static_cast<int>(greedy.size());
        bool fuzzyAccept = (fuzzy > 0) &&
            (std::abs(g1.n - greedySz) <= fuzzy || std::abs(g2.n - greedySz) <= fuzzy);
        if (greedySz >= upperBound) {
            auto greedyC = applyRingAnchorGuard(g1, g2,
                opts.connectedOnly ? largestConnected(g1, greedy, &g2) : greedy, chem);
            if (static_cast<int>(greedyC.size()) >= upperBound) return greedyC;
        }
        if (fuzzyAccept && greedySz > bestSize) {
            auto greedyC = applyRingAnchorGuard(g1, g2,
                opts.connectedOnly ? largestConnected(g1, greedy, &g2) : greedy, chem);
            if (static_cast<int>(greedyC.size()) > bestSize) {
                best = std::move(greedyC);
                bestSize = static_cast<int>(best.size()); bestHetero = heteroatomScore(g1, best);
                bestScore = mcsScore(g1, best, opts);
            }
        }
        // Seed downstream McSplit with greedy result if it covers > 60% of upper bound
        if (static_cast<int>(greedy.size()) > upperBound * 0.6
            && static_cast<int>(greedy.size()) > bestSize) {
            auto greedyC = applyRingAnchorGuard(g1, g2,
                opts.connectedOnly ? largestConnected(g1, greedy, &g2) : greedy, chem);
            if (static_cast<int>(greedyC.size()) > bestSize) {
                best = std::move(greedyC);
                bestSize = static_cast<int>(best.size()); bestHetero = heteroatomScore(g1, best);
                bestScore = mcsScore(g1, best, opts);
            }
        }
    }

    // --- Level 1.25: Augmenting path refinement ---
    // Grow the greedy probe mapping by forced single-candidate extensions.
    if (bestSize > 0 && bestSize < upperBound && !tb.expired()) {
        auto augmented = best;
        bool grew = true;
        while (grew && !tb.expired()) {
            grew = false;
            std::vector<uint8_t> usedQ(g1.n, 0), usedT(g2.n, 0);
            for (auto& [q,t] : augmented) { usedQ[q] = 1; usedT[t] = 1; }
            for (auto& [qi, ti] : augmented) {
                for (int qk : g1.neighbors[qi]) {
                    if (usedQ[qk]) continue;
                    int candidate = -1, count = 0;
                    for (int tk : g2.neighbors[ti]) {
                        if (usedT[tk]) continue;
                        if (atomsCompatFast(g1, qk, g2, tk, chem) &&
                            bondsCompatible(g1, qi, qk, g2, ti, tk, chem)) {
                            candidate = tk; count++;
                        }
                    }
                    if (count == 1) { // forced extension
                        // Verify full consistency with all mapped neighbors of qk
                        bool consistent = true;
                        for (int qm : g1.neighbors[qk]) {
                            auto it = augmented.find(qm);
                            if (it == augmented.end()) continue;
                            int tm = it->second;
                            int qOrd = g1.bondOrder(qk, qm), tOrd = g2.bondOrder(candidate, tm);
                            if ((qOrd != 0) != (tOrd != 0)) { consistent = false; break; }
                            if (qOrd != 0 && !bondsCompatible(g1, qk, qm, g2, candidate, tm, chem))
                                { consistent = false; break; }
                        }
                        if (consistent) {
                            augmented[qk] = candidate;
                            usedQ[qk] = true; usedT[candidate] = true;
                            grew = true;
                        }
                    }
                }
            }
        }
        if (static_cast<int>(augmented.size()) > bestSize) {
            best = augmented; bestSize = static_cast<int>(best.size()); bestHetero = heteroatomScore(g1, best);
            bestScore = mcsScore(g1, best, opts);
        }
        if (!weightMode && bestSize >= upperBound) return best;
    }

    // --- Level 1: Substructure containment ---
    if (!weightMode && g1.n > 0 && g2.n > 0
        && std::min(g1.n, g2.n) <= std::max(g1.n, g2.n) * 3 / 4) {
        const MolGraph& sml = g1.n <= g2.n ? g1 : g2;
        const MolGraph& lrg = g1.n <= g2.n ? g2 : g1;
        bool swapped = g1.n > g2.n;
        // Use VF2PP for substructure check
        auto subMap = findSubstructure(sml, lrg, chem, tb.remainingMs());
        if (!subMap.empty()) {
            if (!swapped) {
                std::map<int,int> result;
                for (auto& [k,v] : subMap) result[k] = v;
                return result;
            }
            std::map<int,int> full;
            for (auto& [k,v] : subMap) full[v] = k;
            return full;
        }
    }

    // Exact branch-and-bound for small disconnected pairs where heuristics
    // may miss the optimum. Deterministic on golden cases (<=20 x <=40 atoms).
    if (!weightMode && opts.disconnectedMCS && !tb.expiredNow()
        && std::min(g1.n, g2.n) <= 20
        && std::max(g1.n, g2.n) <= 40) {
        detail::SmallExactMCSExplorer exactSmall(g1, g2, chem, opts.induced, tb, upperBound, best);
        auto exactMapping = exactSmall.run();
        int exactScore = mcsScore(g1, exactMapping, opts);
        if (isBetterMCS(g1, exactMapping, static_cast<int>(exactMapping.size()),
                        best, bestSize, bestHetero)) {
            best = std::move(exactMapping);
            bestSize = static_cast<int>(best.size());
            bestScore = exactScore;
            bestHetero = heteroatomScore(g1, best);
        }
        if (!weightMode && bestSize >= upperBound) return ppx(g1, g2, best, chem, opts);
    }

    // --- Level 1.5: Seed-and-extend ---
    // Phase 2.4: Use flat-array pipeline to avoid intermediate map construction.
    GraphBuilder GB(g1, g2, chem, opts.induced);
    if (minN >= 4 && std::max(g1.n, g2.n) <= SEED_EXTEND_MAX_ATOMS && !tb.expired()) {
        auto [seFlat, seFlatSize] = GB.seedExtendMCSFlat(tb, upperBound);
        int seScore = mcsScoreFlat(g1, seFlat.data(), g1.n, seFlatSize, opts);
        if (weightMode ? seScore > bestScore : seFlatSize > bestSize) {
            best = flatToMap(seFlat.data(), g1.n);
            bestSize = seFlatSize; bestScore = seScore;
        }
        if (!weightMode && bestSize >= upperBound) return ppx(g1, g2, best, chem, opts);
    }

    // Phase 2.2: route to McGregor early when seed quality is low.
    // Poor seeds on medium-to-large molecules indicate high symmetry where
    // McSplit/BK waste time, while McGregor's bond-grow handles this better.
    if (!weightMode && bestSize < static_cast<int>(minN * 0.3) && minN >= 15 && !tb.expired()) {
        std::vector<std::map<int,int>> earlySeeds;
        if (bestSize > 0) earlySeeds.push_back(best);
        int64_t earlyPerSeedMs = std::max<int64_t>(1, tb.remainingMs() / std::max<int64_t>(1, earlySeeds.size() + 1));
        for (auto& seed : earlySeeds) {
            if (tb.expired()) break;
            auto ext = ppx(g1, g2,
                mcGregorExtend(g1, g2, seed, chem, tb, earlyPerSeedMs,
                               opts.useTwoHopNLFInExtension, opts.useThreeHopNLFInExtension, opts.connectedOnly,
                               &GB.allCompatTargets()),
                chem, opts);
            int earlyScore = mcsScore(g1, ext, opts);
            if (static_cast<int>(ext.size()) > bestSize) {
                best = ext; bestSize = static_cast<int>(ext.size()); bestScore = earlyScore;
                bestHetero = heteroatomScore(g1, best);
            }
        }
        if (bestSize >= upperBound) return ppx(g1, g2, best, chem, opts);
    }

    // --- Level 2: McSplit ---
    // k-core pre-pruning only pays for smaller product graphs.
    // On medium-sized pairs the setup cost can dominate the whole search,
    // and the derived bound is advisory only.
    int kcoreUB = upperBound;
    int pgSize = g1.n * g2.n;
    // Adaptive threshold: allow larger product graphs when time budget permits
    int kcoreLimit = (tb.remainingMs() > opts.timeoutMs / 2) ? 2500 : 1000;
    if (bestSize > 1 && !tb.expired() && pgSize <= kcoreLimit) {
        int k = bestSize;
        int n1k = g1.n, n2k = g2.n;
        std::vector<std::vector<int>> pgDeg(n1k, std::vector<int>(n2k, 0));
        std::vector<std::vector<uint8_t>> alive(n1k, std::vector<uint8_t>(n2k, 0));
        for (int qi = 0; qi < n1k; ++qi)
            for (int tj = 0; tj < n2k; ++tj)
                if (atomsCompatFast(g1, qi, g2, tj, chem)) alive[qi][tj] = true;
        for (int qi = 0; qi < n1k; ++qi)
            for (int tj = 0; tj < n2k; ++tj) {
                if (!alive[qi][tj]) continue;
                int deg = 0;
                for (int qk : g1.neighbors[qi])
                    for (int tl : g2.neighbors[tj])
                        if (alive[qk][tl] && bondsCompatible(g1, qi, qk, g2, tj, tl, chem)) deg++;
                pgDeg[qi][tj] = deg;
            }
        bool changed = true;
        while (changed && !tb.expired()) {
            changed = false;
            for (int qi = 0; qi < n1k; ++qi)
                for (int tj = 0; tj < n2k; ++tj) {
                    if (!alive[qi][tj] || pgDeg[qi][tj] >= k - 1) continue;
                    alive[qi][tj] = false;
                    changed = true;
                    for (int qk : g1.neighbors[qi])
                        for (int tl : g2.neighbors[tj])
                            if (alive[qk][tl]) pgDeg[qk][tl]--;
                }
        }
        int kcUB = 0;
        for (int qi = 0; qi < n1k; ++qi) {
            bool has = false;
            for (int tj = 0; tj < n2k; ++tj)
                if (alive[qi][tj]) { has = true; break; }
            if (has) kcUB++;
        }
        if (kcUB < kcoreUB) kcoreUB = kcUB;
    }
    // Phase 2.4: Use flat-array pipeline for McSplit to avoid intermediate map construction.
    bool mcSplitExhaustive;
    int mcSplitSize;
    {
        int64_t nodeCount = 0;
        auto [mcFlat, mcFlatSize] = GB.mcSplitSeedFlat(tb, nodeCount);
        mcSplitSize = mcFlatSize;
        int mcScore = mcsScoreFlat(g1, mcFlat.data(), g1.n, mcFlatSize, opts);
        if (weightMode ? mcScore > bestScore : (mcFlatSize > bestSize || (mcFlatSize == bestSize && mcFlatSize > 0 && heteroatomScoreFlat(g1, mcFlat.data(), g1.n) > bestHetero))) {
            best = flatToMap(mcFlat.data(), g1.n);
            bestSize = mcFlatSize; bestScore = mcScore;
            bestHetero = heteroatomScoreFlat(g1, mcFlat.data(), g1.n);
        }
        mcSplitExhaustive = nodeCount < 200000;
    }
    if (!weightMode && mcSplitExhaustive && bestSize >= upperBound)
        return ppx(g1, g2, best, chem, opts);

    // Phase 2.2: bail on unreachable LFUB -- if McSplit was exhaustive and
    // bestSize >= 85% of upperBound, BK cannot improve meaningfully. Skip it
    // and fall through to McGregor extension which can close the small gap.
    bool skipBKDueToLFUB = !weightMode && mcSplitExhaustive
        && bestSize >= static_cast<int>(upperBound * 0.85) && bestSize > 0;

    // Near-optimal fast path: greedy atom extension
    if (bestSize >= upperBound - 2 && bestSize > 0 && bestSize < upperBound && !tb.expired()) {
        auto ext = ppx(g1, g2, greedyAtomExtend(g1, g2, ppx(g1, g2, best, chem, opts), chem, opts), chem, opts);
        int extScore = mcsScore(g1, ext, opts);
        if (weightMode ? extScore > bestScore : static_cast<int>(ext.size()) > bestSize) {
            best = ext; bestSize = static_cast<int>(ext.size()); bestScore = extScore;
        }
        if (!weightMode && bestSize >= upperBound) return best;
    }

    // --- Level 3: Bron-Kerbosch ---
    int bkSize = 0;
    if (!skipBKDueToLFUB && bestSize < static_cast<int>(upperBound * BK_SKIP_RATIO) && !tb.expired()) {
        auto cliqueSeed = GB.maximumCliqueSeed(tb);
        bkSize = static_cast<int>(cliqueSeed.size());
        int cScore = mcsScore(g1, cliqueSeed, opts);
        if (weightMode ? cScore > bestScore : static_cast<int>(cliqueSeed.size()) > bestSize) {
            best = cliqueSeed; bestSize = static_cast<int>(cliqueSeed.size()); bestScore = cScore;
        }
    }
    if (!weightMode && bestSize >= upperBound) return ppx(g1, g2, best, chem, opts);

    // --- Level 4: McGregor extension with seeds ---
    std::vector<std::map<int,int>> seeds;
    if (bestSize > 0) seeds.push_back(best);

    // --- Level 5: Extra seeds (only when McSplit and BK disagree) ---
    bool skipExtraSeeds = mcSplitSize > 0 && bkSize > 0 && mcSplitSize == bkSize
        && mcSplitSize >= upperBound - 1;
    if (opts.extraSeeds && !skipExtraSeeds && bestSize < upperBound && !tb.expired()) {
        auto s = GB.ringAnchorSeed(tb);
        if (!s.empty()) seeds.push_back(s);
        if (!tb.expired()) {
            s = GB.labelDegreeAnchorSeed(tb);
            if (!s.empty()) seeds.push_back(s);
        }
        if (!tb.expired()) {
            s = GB.vf2ppRingSkeletonSeed(tb,
                std::max<int64_t>(1, tb.remainingMs() / 4), 2, 12);
            if (!s.empty()) seeds.push_back(s);
        }
        if (!tb.expired()) {
            s = GB.vf2ppCoreSeed(tb,
                std::max<int64_t>(1, tb.remainingMs() / 4), 2, 12);
            if (!s.empty()) seeds.push_back(s);
        }
    }

    int64_t perSeedMs = std::max<int64_t>(1, tb.remainingMs() / std::max<int64_t>(1, seeds.size()));
    for (auto& seed : seeds) {
        if (tb.expired()) break;
        auto ext = ppx(g1, g2,
            mcGregorExtend(g1, g2, seed, chem, tb, perSeedMs,
                           opts.useTwoHopNLFInExtension, opts.useThreeHopNLFInExtension, opts.connectedOnly,
                           &GB.allCompatTargets()),
            chem, opts);
        int seedScore = mcsScore(g1, ext, opts);
        if (weightMode ? seedScore > bestScore : static_cast<int>(ext.size()) > bestSize) {
            best = ext; bestSize = static_cast<int>(ext.size()); bestScore = seedScore;
        }
        if (!weightMode && bestSize >= upperBound) return best;
    }

    // Last resort: start from empty seed
    if (bestScore <= 0 && !tb.expired()) {
        best = ppx(g1, g2,
            mcGregorExtend(g1, g2, {}, chem, tb, tb.remainingMs(),
                           opts.useTwoHopNLFInExtension, opts.useThreeHopNLFInExtension, opts.connectedOnly,
                           &GB.allCompatTargets()),
            chem, opts);
    }
    // Apply post-processing (connectivity filter, ring guard) before returning.
    // All early-exit paths call ppx() explicitly; this covers the fall-through path
    // where intermediate stages (McSplit, BK) may have left a disconnected mapping
    // in `best` that was never passed through largestConnected().
    best = ppx(g1, g2, best, chem, opts);
    // Invariant: MCS size must never exceed the smaller molecule's atom count
    assert(static_cast<int>(best.size()) <= std::min(g1.n, g2.n)
           && "MCS size exceeds min(n1, n2)");
    return best;
}

inline bool isValidMcsMapping(const MolGraph& g1, const MolGraph& g2,
                              const std::map<int,int>& mapping,
                              const ChemOptions& opts) {
    if (mapping.empty()) return true;

    std::unordered_set<int> usedT;
    for (const auto& [qi, tj] : mapping) {
        if (qi < 0 || qi >= g1.n) return false;
        if (tj < 0 || tj >= g2.n) return false;
        if (!usedT.insert(tj).second) return false;
        if (!detail::atomsCompatFast(g1, qi, g2, tj, opts)) return false;
    }

    for (const auto& [qi, tj] : mapping) {
        for (int qk : g1.neighbors[qi]) {
            if (qk <= qi) continue;
            auto it = mapping.find(qk);
            if (it == mapping.end()) continue;
            int tk = it->second;
            if (g1.bondOrder(qi, qk) == 0) continue;
            if (g2.bondOrder(tj, tk) == 0) return false;
            if (!ChemOps::bondsCompatible(g1, qi, qk, g2, tj, tk, opts)) {
                return false;
            }
        }
    }
    return true;
}

inline std::map<int,int> repairInvalidMcsMapping(const MolGraph& g1, const MolGraph& g2,
                                                 std::map<int,int> mapping,
                                                 const ChemOptions& opts);
inline std::vector<std::string> validateMapping(const MolGraph& g1, const MolGraph& g2,
                                                const std::map<int,int>& mapping,
                                                const ChemOptions& opts);

inline int64_t resolveMcsTimeoutMs(const MolGraph& g1, const MolGraph& g2,
                                   const MCSOptions& opts) {
    if (opts.timeoutMs >= 0) return std::max<int64_t>(1, opts.timeoutMs);
    return std::max<int64_t>(1, std::min<int64_t>(30000, 500 + int64_t(g1.n) * g2.n * 2));
}

inline MCSOptions withTimeoutMs(const MCSOptions& opts, int64_t timeoutMs) {
    MCSOptions limited = opts;
    limited.timeoutMs = std::max<int64_t>(1, timeoutMs);
    return limited;
}

inline std::map<int,int> recoverValidMcsMapping(const MolGraph& g1, const MolGraph& g2,
                                                const std::map<int,int>& raw,
                                                const ChemOptions& chem,
                                                const MCSOptions& opts) {
    auto repaired = repairInvalidMcsMapping(g1, g2, raw, chem);
    if (!validateMapping(g1, g2, repaired, chem).empty()) repaired.clear();
    int repairedSize = static_cast<int>(repaired.size());
    int64_t timeoutMs = resolveMcsTimeoutMs(g1, g2, opts);

    const MolGraph& sml = g1.n <= g2.n ? g1 : g2;
    const MolGraph& lrg = g1.n <= g2.n ? g2 : g1;
    bool swapped = g1.n > g2.n;
    auto subMaps = findAllSubstructures(sml, lrg, chem, std::max<int64_t>(1, timeoutMs / 4));
    if (!subMaps.empty()) {
        std::map<int,int> recovered;
        for (const auto& p : subMaps.front()) {
            if (swapped) recovered[p.second] = p.first;
            else recovered[p.first] = p.second;
        }
        recovered = detail::ppx(g1, g2, recovered, chem, opts);
        if (validateMapping(g1, g2, recovered, chem).empty()
            && static_cast<int>(recovered.size()) > repairedSize) {
            return recovered;
        }
    }

    if (!repaired.empty()) {
        auto regrown = detail::ppx(g1, g2, detail::greedyAtomExtend(g1, g2, repaired, chem, opts), chem, opts);
        regrown = repairInvalidMcsMapping(g1, g2, std::move(regrown), chem);
        if (validateMapping(g1, g2, regrown, chem).empty()
            && static_cast<int>(regrown.size()) > repairedSize) {
            return regrown;
        }
    }
    return repaired;
}

inline std::map<int,int> orientMcsResult(const std::map<int,int>& mapping, bool swapped) {
    if (!swapped) return mapping;
    std::map<int,int> restored;
    for (const auto& p : mapping) restored[p.second] = p.first;
    return restored;
}

inline std::map<int,int> runValidatedMcsDirection(const MolGraph& query, const MolGraph& target,
                                                  const ChemOptions& chem,
                                                  const MCSOptions& opts,
                                                  bool swappedBack) {
    auto raw = findMCSImpl(query, target, chem, opts);
    auto valid = validateMapping(query, target, raw, chem).empty()
        ? raw
        : recoverValidMcsMapping(query, target, raw, chem, opts);
    bool weightMode = opts.maximizeBonds || !opts.atomWeights.empty();
    int minN = std::min(query.n, target.n);
    int maxN = std::max(query.n, target.n);
    if (!weightMode
        && !valid.empty()
        && static_cast<int>(valid.size()) == minN
        && minN <= 20
        && maxN <= 40) {
        int64_t refineMs = std::max<int64_t>(
            1, std::min<int64_t>(2000, resolveMcsTimeoutMs(query, target, opts)));
        detail::TimeBudget refineBudget(refineMs);
        detail::FixedSizeBondMaximizer bondRefiner(
            query, target, chem, opts.induced, refineBudget, minN, valid);
        auto refined = bondRefiner.run();
        if (validateMapping(query, target, refined, chem).empty()
            && detail::preferFinalMapping(query, refined, valid, opts)) {
            valid = std::move(refined);
        }
    }
    return orientMcsResult(valid, swappedBack);
}

inline std::map<int,int> findMCSDirectionalCore(const MolGraph& g1, const MolGraph& g2,
                                                const ChemOptions& chem, const MCSOptions& opts) {
    if (g1.n == 0 || g2.n == 0) return {};
    if (&g1 == &g2 || detail::isExactMatch(g1, g2, chem)) {
        std::map<int,int> id;
        for (int i = 0; i < g1.n; ++i) id[i] = i;
        if (opts.disconnectedMCS && (opts.minFragmentSize > 1 || opts.maxFragments < INT_MAX))
            id = detail::applyFragmentConstraints(g1, id, opts.minFragmentSize, opts.maxFragments);
        return id;
    }

    using Clock = std::chrono::steady_clock;
    auto deadline = Clock::now() + std::chrono::milliseconds(resolveMcsTimeoutMs(g1, g2, opts));
    auto remainingMs = [&]() -> int64_t {
        auto left = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - Clock::now()).count();
        return std::max<int64_t>(0, left);
    };

    g1.ensureCanonical();
    g2.ensureCanonical();
    if (chem.ringFusionMode != ChemOptions::RingFusionMode::IGNORE) {
        g1.ensureRingCounts();
        g2.ensureRingCounts();
    }

    int ub12 = detail::labelFrequencyUpperBoundDirected(g1, g2, chem);
    int ub21 = detail::labelFrequencyUpperBoundDirected(g2, g1, chem);
    auto plan = detail::chooseOrientationPlan(g1, g2, chem, ub12, ub21);
    bool weightMode = opts.maximizeBonds || !opts.atomWeights.empty();

    auto runDirection = [&](bool direct) {
        int64_t budgetMs = remainingMs();
        if (budgetMs <= 0) return std::map<int,int>{};
        auto timedOpts = withTimeoutMs(opts, budgetMs);
        auto oriented = direct
            ? runValidatedMcsDirection(g1, g2, chem, timedOpts, false)
            : runValidatedMcsDirection(g2, g1, chem, timedOpts, true);
        if (validateMapping(g1, g2, oriented, chem).empty()) return oriented;
        budgetMs = remainingMs();
        if (budgetMs <= 0) return oriented;
        return recoverValidMcsMapping(
            g1, g2, oriented, chem, withTimeoutMs(opts, budgetMs));
    };

    auto best = runDirection(plan.directFirst);
    int baseUb = opts.induced
        ? detail::degreeSequenceUpperBound(g1, g2, chem)
        : detail::labelFrequencyUpperBound(g1, g2, chem);
    if (!opts.induced) {
        baseUb = std::min(baseUb, ub12);
        baseUb = std::min(baseUb, ub21);
    }

    int alternateFloor = detail::alternateOrientationFloor(plan, plan.directFirst);
    int probeCeiling = std::max(std::max(plan.seed12, plan.seed21),
                                std::max(plan.mc12, plan.mc21));
    bool runAlternate = best.empty();
    // Phase 2.2: skip redundant orientation retries when both directed upper
    // bounds equal the current best -- the alternate direction cannot improve.
    if (!runAlternate && !weightMode
        && static_cast<int>(best.size()) >= ub12
        && static_cast<int>(best.size()) >= ub21) {
        runAlternate = false; // both bounds matched, alternate is redundant
    }
    // Also skip if the forward direction already matched the tighter bound.
    else if (!runAlternate && !weightMode) {
        int tighterUb = plan.directFirst ? ub12 : ub21;
        if (static_cast<int>(best.size()) >= tighterUb) {
            runAlternate = false; // forward matched its own bound, skip alternate
        }
    }
    if (!runAlternate && static_cast<int>(best.size()) < alternateFloor) runAlternate = true;
    if (!runAlternate && !weightMode && static_cast<int>(best.size()) + 1 < probeCeiling) runAlternate = true;
    if (!runAlternate && !weightMode && baseUb >= 12
        && static_cast<int>(best.size()) + 4 < baseUb) {
        runAlternate = true;
    }
    if (!runAlternate && !weightMode
        && std::abs(g1.n - g2.n) <= 1
        && std::min(g1.n, g2.n) >= 20) {
        runAlternate = true;
    }
    if (!runAlternate && !weightMode
        && std::abs(g1.n - g2.n) <= 2
        && std::min(g1.n, g2.n) >= 20
        && static_cast<int>(best.size()) * 4 < std::min(g1.n, g2.n) * 3) {
        runAlternate = true;
    }
    if (!runAlternate && !weightMode && static_cast<int>(best.size()) + 2 < baseUb
        && std::abs(g1.n - g2.n) >= 4) runAlternate = true;

    if (runAlternate) {
        auto alt = runDirection(!plan.directFirst);
        if (detail::preferFinalMapping(g1, alt, best, opts)) best = std::move(alt);
    }
    return best;
}

inline std::map<int,int> findMCS(const MolGraph& g1, const MolGraph& g2,
                                 const ChemOptions& chem, const MCSOptions& opts) {
    if (g1.n == 0 || g2.n == 0) return {};
    if (&g1 == &g2 || detail::isExactMatch(g1, g2, chem)) {
        std::map<int,int> id;
        for (int i = 0; i < g1.n; ++i) id[i] = i;
        if (opts.disconnectedMCS && (opts.minFragmentSize > 1 || opts.maxFragments < INT_MAX))
            id = detail::applyFragmentConstraints(g1, id, opts.minFragmentSize, opts.maxFragments);
        return id;
    }

    using Clock = std::chrono::steady_clock;
    auto deadline = Clock::now() + std::chrono::milliseconds(resolveMcsTimeoutMs(g1, g2, opts));
    auto remainingMs = [&]() -> int64_t {
        auto left = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - Clock::now()).count();
        return std::max<int64_t>(0, left);
    };
    auto runDirectionalCore = [&](const MolGraph& lhs, const MolGraph& rhs) {
        int64_t budgetMs = remainingMs();
        if (budgetMs <= 0) return std::map<int,int>{};
        return findMCSDirectionalCore(lhs, rhs, chem, withTimeoutMs(opts, budgetMs));
    };

    auto best = runDirectionalCore(g1, g2);
    if (best.empty()) return best;

    bool weightMode = opts.maximizeBonds || !opts.atomWeights.empty();
    if (!weightMode) {
        int globalUb = opts.induced
            ? detail::degreeSequenceUpperBound(g1, g2, chem)
            : detail::labelFrequencyUpperBound(g1, g2, chem);
        if (!opts.induced) {
            globalUb = std::min(globalUb, detail::labelFrequencyUpperBoundDirected(g1, g2, chem));
            globalUb = std::min(globalUb, detail::labelFrequencyUpperBoundDirected(g2, g1, chem));
        }
        if (static_cast<int>(best.size()) >= globalUb) return best;
        if (remainingMs() <= 0) return best;

        auto reverse = runDirectionalCore(g2, g1);
        int bestSize = static_cast<int>(best.size());
        int reverseSize = static_cast<int>(reverse.size());
        if (std::abs(bestSize - reverseSize) >= 2 && reverseSize < bestSize) {
            auto consensus = orientMcsResult(reverse, true);
            if (!validateMapping(g1, g2, consensus, chem).empty()) {
                consensus = recoverValidMcsMapping(g1, g2, consensus, chem, opts);
            }
            if (!consensus.empty()) best = std::move(consensus);
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// MCSProgressFn -- optional callback for progressive MCS reporting
// ---------------------------------------------------------------------------

/// Progress callback type: receives (bestMappingSoFar, bestSize, elapsedMs).
/// Called after each major pipeline level with the current best result.
using MCSProgressFn = std::function<void(const std::map<int,int>&, int, int64_t)>;

/// Find MCS with optional progress callback.
/// The callback is invoked after each pipeline level (L0.25 chain, L0.5 tree,
/// L0.75 greedy, L1 substructure, L1.5 seed-extend, L2 McSplit, L3 BK, L4/L5
/// McGregor) with the current best mapping, its size, and elapsed time in ms.
/// Passing nullptr as progress disables progress reporting.
inline std::map<int,int> findMCS(const MolGraph& g1, const MolGraph& g2,
                                  const ChemOptions& chem, const MCSOptions& opts,
                                  MCSProgressFn progress) {
    auto t0 = std::chrono::steady_clock::now();
    auto elapsedMs = [&]() -> int64_t {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - t0).count();
    };
    auto reportProgress = [&](const std::map<int,int>& best, int bestSize) {
        if (progress) progress(best, bestSize, elapsedMs());
    };
    auto best = findMCS(g1, g2, chem, opts);
    reportProgress(best, static_cast<int>(best.size()));
    return best;
}

/// Find the disconnected MCS (dMCS) between two molecular graphs.
inline std::map<int,int> findDisconnectedMCS(const MolGraph& g1, const MolGraph& g2,
                                              const ChemOptions& chem, const MCSOptions& opts) {
    MCSOptions dM = opts;
    dM.connectedOnly = false;
    dM.disconnectedMCS = true;
    return findMCS(g1, g2, chem, dM);
}

// ========================================================================
// Phase 2.4: Flat-array MCS pipeline
// ========================================================================

/// Internal flat-array MCS entry point.
/// Returns a q2t flat array (where -1 = unmapped) and the mapping size.
/// Avoids std::map<int,int> construction at the API boundary for internal
/// consumers that only need the flat representation (e.g., batchMCSSize).
///
/// The public API signature findMCS() returning std::map<int,int> is unchanged.
inline std::pair<std::vector<int>, int> findMCSFlat(
    const MolGraph& g1, const MolGraph& g2,
    const ChemOptions& chem, const MCSOptions& opts)
{
    auto mapping = findMCS(g1, g2, chem, opts);
    int n1 = g1.n;
    std::vector<int> q2t(std::max(n1, 1), -1);
    for (const auto& [k, v] : mapping) {
        if (k >= 0 && k < n1) q2t[k] = v;
    }
    return {std::move(q2t), static_cast<int>(mapping.size())};
}

/// Count-only MCS: returns just the MCS size without materializing
/// the mapping at the caller's side.  Uses findMCS internally but
/// avoids returning the map to the caller.
inline int findMCSSize(const MolGraph& g1, const MolGraph& g2,
                       const ChemOptions& chem, const MCSOptions& opts) {
    return static_cast<int>(findMCS(g1, g2, chem, opts).size());
}

// ========================================================================
// Mapping canonicalization by automorphism
// ========================================================================

/// Canonicalize an atom-atom mapping under the automorphism groups of both
/// molecules.  Two mappings M1 and M2 are equivalent if there exist
/// automorphisms alpha in Aut(g1) and beta in Aut(g2) such that
/// M2 = beta . M1 . alpha^{-1}.  This function returns the lexicographically
/// smallest representative of the equivalence class.
///
/// The algorithm iteratively applies each generator to the current best
/// mapping (on both the query side and the target side) and keeps the
/// lex-smallest result, repeating until a fixed point is reached.
inline std::map<int,int> canonicalizeMapping(
        const MolGraph& g1, const MolGraph& g2,
        const std::map<int,int>& mapping) {
    if (mapping.empty()) return mapping;

    const auto& gens1 = g1.getAutomorphismGenerators();
    const auto& gens2 = g2.getAutomorphismGenerators();

    // No symmetry in either molecule — mapping is already canonical.
    if (gens1.empty() && gens2.empty()) return mapping;

    // Work with sorted pair vectors for efficient lex comparison.
    using PairVec = std::vector<std::pair<int,int>>;
    auto mapToPairs = [](const std::map<int,int>& m) {
        PairVec v(m.begin(), m.end());
        return v;
    };

    // Precompute inverses of g1 generators (query-side permutation).
    std::vector<std::vector<int>> invGens1;
    invGens1.reserve(gens1.size());
    for (const auto& gen : gens1) {
        std::vector<int> inv(gen.size());
        for (int i = 0; i < static_cast<int>(gen.size()); ++i)
            inv[gen[i]] = i;
        invGens1.push_back(std::move(inv));
    }

    PairVec best = mapToPairs(mapping);
    PairVec candidate;
    candidate.reserve(best.size());

    int maxIter = std::max(100, 2 * static_cast<int>(gens1.size() + gens2.size()));
    for (int iter = 0; iter < maxIter; ++iter) {
        bool improved = false;

        // Try each g1-side generator (permute query atom indices).
        for (const auto& inv : invGens1) {
            candidate.clear();
            for (const auto& [qi, ti] : best)
                candidate.emplace_back(inv[qi], ti);
            std::sort(candidate.begin(), candidate.end());
            if (candidate < best) { best = candidate; improved = true; }
        }

        // Try each g2-side generator (permute target atom indices).
        for (const auto& gen : gens2) {
            candidate.clear();
            for (const auto& [qi, ti] : best)
                candidate.emplace_back(qi, gen[ti]);
            std::sort(candidate.begin(), candidate.end());
            if (candidate < best) { best = candidate; improved = true; }
        }

        if (!improved) break;
    }

    std::map<int,int> result;
    for (const auto& [qi, ti] : best) result[qi] = ti;
    return result;
}

/// Enumerate multiple distinct MCS mappings of the maximum size.
///
/// First computes the single best MCS to determine the optimal size K, then
/// re-searches the solution space to collect up to @p maxResults distinct
/// atom-atom mappings of that same size K. Useful for scaffold analysis and
/// SAR studies where alternative mappings provide different chemical perspectives.
///
/// @param g1         first molecule graph
/// @param g2         second molecule graph
/// @param chem       chemical matching options
/// @param opts       MCS options (timeout, induced, connected, etc.)
/// @param maxResults maximum number of distinct mappings to return (default 10)
/// @return vector of distinct atom-index mappings, each of the maximum MCS size
inline std::vector<std::map<int,int>> findAllMCS(const MolGraph& g1, const MolGraph& g2,
                                                  const ChemOptions& chem, const MCSOptions& opts,
                                                  int maxResults = 10) {
    using namespace detail;
    if (g1.n == 0 || g2.n == 0) return {};
    if (maxResults <= 0) maxResults = 10;

    // Phase 1: find the optimal MCS size
    std::map<int,int> best = findMCS(g1, g2, chem, opts);
    int K = static_cast<int>(best.size());
    if (K == 0) return {};

    // Deduplicate by automorphism-canonical mapping.
    // Two mappings are equivalent if they differ only by automorphisms of g1/g2.
    using CanonKey = std::vector<std::pair<int,int>>;
    auto canonKey = [&](const std::map<int,int>& m) -> CanonKey {
        auto cm = canonicalizeMapping(g1, g2, m);
        return CanonKey(cm.begin(), cm.end());
    };
    std::map<CanonKey, std::map<int,int>> seen;
    seen.emplace(canonKey(best), best);
    if (static_cast<int>(seen.size()) >= maxResults) {
        std::vector<std::map<int,int>> result;
        for (auto& [k, v] : seen) result.push_back(std::move(v));
        return result;
    }

    // Resolve adaptive timeout for enumeration phase
    int64_t enumTimeout = opts.timeoutMs;
    if (enumTimeout < 0) {
        enumTimeout = std::min(int64_t(30000), int64_t(500) + int64_t(g1.n) * g2.n * 2);
    }
    TimeBudget tb(enumTimeout);
    int minN = std::min(g1.n, g2.n);

    // Phase 2a: if MCS == smaller molecule, enumerate substructure mappings
    if (K == minN && !tb.expired()) {
        const MolGraph& sml = g1.n <= g2.n ? g1 : g2;
        const MolGraph& lrg = g1.n <= g2.n ? g2 : g1;
        bool swapped = g1.n > g2.n;
        auto subMaps = findAllSubstructures(sml, lrg, chem, tb.remainingMs());
        for (const auto& raw : subMaps) {
            if (tb.expired() || static_cast<int>(seen.size()) >= maxResults) break;
            std::map<int,int> mapping;
            if (swapped) {
                for (const auto& p : raw) mapping[p.second] = p.first;
            } else {
                for (const auto& p : raw) mapping[p.first] = p.second;
            }
            if (static_cast<int>(mapping.size()) == K) {
                auto ck = canonKey(mapping);
                if (seen.find(ck) == seen.end())
                    seen.emplace(std::move(ck), mapping);
            }
        }
        if (static_cast<int>(seen.size()) >= maxResults) {
            std::vector<std::map<int,int>> result;
            for (auto& [k, v] : seen) result.push_back(std::move(v));
            return result;
        }
    }

    // Phase 2b: re-run pipeline stages to collect alternative seeds
    if (!tb.expired()) {
        GraphBuilder GB(g1, g2, chem, opts.induced);

        std::vector<std::map<int,int>> seeds;

        // McSplit seed
        if (!tb.expired()) {
            int64_t nc = 0;
            auto mcSeed = GB.mcSplitSeed(tb, nc);
            if (static_cast<int>(mcSeed.size()) >= K) seeds.push_back(std::move(mcSeed));
        }

        // BK clique seed
        if (!tb.expired()) {
            auto cliqueSeed = GB.maximumCliqueSeed(tb);
            if (static_cast<int>(cliqueSeed.size()) >= K) seeds.push_back(std::move(cliqueSeed));
        }

        // Extra seeds
        if (opts.extraSeeds && !tb.expired()) {
            auto s = GB.ringAnchorSeed(tb);
            if (!s.empty()) seeds.push_back(std::move(s));
            if (!tb.expired()) {
                s = GB.labelDegreeAnchorSeed(tb);
                if (!s.empty()) seeds.push_back(std::move(s));
            }
            if (!tb.expired()) {
                s = GB.vf2ppRingSkeletonSeed(tb, std::max(int64_t(1), tb.remainingMs() / 4), 2, 12);
                if (!s.empty()) seeds.push_back(std::move(s));
            }
            if (!tb.expired()) {
                s = GB.vf2ppCoreSeed(tb, std::max(int64_t(1), tb.remainingMs() / 4), 2, 12);
                if (!s.empty()) seeds.push_back(std::move(s));
            }
        }

        // Extend each seed and collect K-sized results
        int64_t perSeedMs = std::max(int64_t(1), tb.remainingMs() / std::max(int(1), static_cast<int>(seeds.size())));
        for (auto& seed : seeds) {
            if (tb.expired() || static_cast<int>(seen.size()) >= maxResults) break;
            auto ext = ppx(g1, g2,
                mcGregorExtend(g1, g2, seed, chem, tb, perSeedMs,
                    opts.useTwoHopNLFInExtension, opts.useThreeHopNLFInExtension, opts.connectedOnly,
                    &GB.allCompatTargets()),
                chem, opts);
            if (static_cast<int>(ext.size()) == K) {
                auto ck = canonKey(ext);
                if (seen.find(ck) == seen.end())
                    seen.emplace(std::move(ck), ext);
            }

            // Greedy atom extension for alternative mappings
            if (!tb.expired() && static_cast<int>(seen.size()) < maxResults) {
                auto gext = ppx(g1, g2,
                    greedyAtomExtend(g1, g2, seed, chem, opts),
                    chem, opts);
                if (static_cast<int>(gext.size()) == K) {
                    auto ck = canonKey(gext);
                    if (seen.find(ck) == seen.end())
                        seen.emplace(std::move(ck), gext);
                }
            }
        }

        // Phase 2c: perturb existing solutions
        if (!tb.expired() && static_cast<int>(seen.size()) < maxResults) {
            std::vector<std::map<int,int>> existing;
            for (const auto& [k, v] : seen) existing.push_back(v);
            for (const auto& ex : existing) {
                if (tb.expired() || static_cast<int>(seen.size()) >= maxResults) break;
                for (const auto& entry : ex) {
                    if (tb.expired() || static_cast<int>(seen.size()) >= maxResults) break;
                    std::map<int,int> reduced = ex;
                    reduced.erase(entry.first);
                    auto reext = ppx(g1, g2,
                        greedyAtomExtend(g1, g2, reduced, chem, opts),
                        chem, opts);
                    if (static_cast<int>(reext.size()) == K) {
                        auto ck = canonKey(reext);
                        if (seen.find(ck) == seen.end())
                            seen.emplace(std::move(ck), reext);
                    }
                }
            }
        }
    }

    std::vector<std::map<int,int>> result;
    result.reserve(seen.size());
    for (auto& [k, v] : seen) result.push_back(std::move(v));
    return result;
}

/// Test whether two MCS mappings are equivalent under automorphism.
/// Two mappings are equivalent if their canonicalized forms are identical.
inline bool areMappingsEquivalent(
        const MolGraph& g1, const MolGraph& g2,
        const std::map<int,int>& m1, const std::map<int,int>& m2) {
    return canonicalizeMapping(g1, g2, m1) == canonicalizeMapping(g1, g2, m2);
}

/// Compute RASCAL-style similarity upper bound (Tanimoto-like, 0.0--1.0).
inline double similarityUpperBound(const MolGraph& g1, const MolGraph& g2) {
    if (g1.n == 0 || g2.n == 0) return 0.0;
    std::unordered_map<int,int> count1, count2;
    for (int i = 0; i < g1.n; ++i) count1[g1.label[i]]++;
    for (int i = 0; i < g2.n; ++i) count2[g2.label[i]]++;
    int atomOverlap = 0;
    for (auto& [lab, cnt] : count1) {
        auto it = count2.find(lab);
        if (it != count2.end()) atomOverlap += std::min(cnt, it->second);
    }
    if (atomOverlap == 0) return 0.0;
    int total = g1.n + g2.n;
    return static_cast<double>(atomOverlap) / (total - atomOverlap);
}

// ===========================================================================
// validateMapping -- check that a mapping is chemically correct.
// Port of Java SearchEngine.validateMapping()
// Checks: injective, atoms compatible, bonds compatible for all mapped neighbors.
// ===========================================================================

inline std::vector<std::string> validateMapping(const MolGraph& g1, const MolGraph& g2,
                                                 const std::map<int,int>& mapping,
                                                 const ChemOptions& opts) {
    std::vector<std::string> errors;
    if (mapping.empty()) return errors;

    // Check injectivity: no target atom mapped twice
    std::unordered_set<int> usedT;
    for (auto& [qi, tj] : mapping) {
        if (qi < 0 || qi >= g1.n) {
            errors.push_back("query index " + std::to_string(qi) + " out of range [0," + std::to_string(g1.n) + ")");
            continue;
        }
        if (tj < 0 || tj >= g2.n) {
            errors.push_back("target index " + std::to_string(tj) + " out of range [0," + std::to_string(g2.n) + ")");
            continue;
        }
        if (!usedT.insert(tj).second) {
            errors.push_back("target atom " + std::to_string(tj) + " mapped twice (not injective)");
        }
        if (!detail::atomsCompatFast(g1, qi, g2, tj, opts)) {
            errors.push_back("atom mismatch: q" + std::to_string(qi) + "(" + std::to_string(g1.atomicNum[qi])
                            + ") -> t" + std::to_string(tj) + "(" + std::to_string(g2.atomicNum[tj]) + ")");
        }
    }

    // Check bonds: every mapped query bond must exist in the target and be compatible.
    for (auto& [qi, tj] : mapping) {
        if (qi < 0 || qi >= g1.n || tj < 0 || tj >= g2.n) continue;
        for (int qk : g1.neighbors[qi]) {
            if (qk <= qi) continue;
            auto it = mapping.find(qk);
            if (it == mapping.end()) continue;
            int tk = it->second;
            if (g1.bondOrder(qi, qk) == 0) continue;
            if (g2.bondOrder(tj, tk) == 0) {
                errors.push_back("bond missing: q(" + std::to_string(qi) + "-" + std::to_string(qk)
                                + ") vs t(" + std::to_string(tj) + "-" + std::to_string(tk) + ")");
                continue;
            }
            if (!ChemOps::bondsCompatible(g1, qi, qk, g2, tj, tk, opts)) {
                errors.push_back("bond incompatible: q(" + std::to_string(qi) + "-" + std::to_string(qk)
                                + ") vs t(" + std::to_string(tj) + "-" + std::to_string(tk) + ")");
            }
        }
    }
    return errors;
}

inline std::map<int,int> repairInvalidMcsMapping(const MolGraph& g1, const MolGraph& g2,
                                                 std::map<int,int> mapping,
                                                 const ChemOptions& opts) {
    while (!mapping.empty()) {
        auto errors = validateMapping(g1, g2, mapping, opts);
        if (errors.empty()) return mapping;

        std::vector<int> candidates;
        const std::string& err = errors.front();
        auto addCandidate = [&](int qi) {
            if (qi < 0) return;
            if (std::find(candidates.begin(), candidates.end(), qi) == candidates.end()) {
                candidates.push_back(qi);
            }
        };
        auto parseInt = [&](size_t pos) -> int {
            size_t end = pos;
            if (end < err.size() && err[end] == '-') ++end;
            while (end < err.size() && std::isdigit(static_cast<unsigned char>(err[end]))) ++end;
            return std::stoi(err.substr(pos, end - pos));
        };

        if (err.rfind("atom mismatch: q", 0) == 0) {
            size_t qPos = err.find('q');
            if (qPos != std::string::npos) addCandidate(parseInt(qPos + 1));
        } else if (err.rfind("bond missing: q(", 0) == 0
                   || err.rfind("bond incompatible: q(", 0) == 0) {
            size_t qStart = err.find('(');
            size_t dash = err.find('-', qStart);
            if (qStart != std::string::npos && dash != std::string::npos) {
                addCandidate(parseInt(qStart + 1));
                addCandidate(parseInt(dash + 1));
            }
        }
        if (candidates.empty()) addCandidate(mapping.begin()->first);

        std::map<int,int> best;
        for (int qi : candidates) {
            auto reduced = mapping;
            reduced.erase(qi);
            auto repaired = repairInvalidMcsMapping(g1, g2, std::move(reduced), opts);
            if (repaired.size() > best.size()) best = std::move(repaired);
        }
        return best;
    }
    return {};
}

// ===========================================================================
// isMappingMaximal -- check if no unmapped (qi,tj) pair can be added.
// Port of Java SearchEngine.isMappingMaximal()
// Try adding any unmapped pair -- if none can extend, it's maximal.
// ===========================================================================

inline bool isMappingMaximal(const MolGraph& g1, const MolGraph& g2,
                              const std::map<int,int>& mapping,
                              const ChemOptions& opts) {
    std::unordered_set<int> usedQ, usedT;
    for (auto& [qi, tj] : mapping) {
        usedQ.insert(qi);
        usedT.insert(tj);
    }

    for (int qi = 0; qi < g1.n; ++qi) {
        if (usedQ.count(qi)) continue;
        for (int tj = 0; tj < g2.n; ++tj) {
            if (usedT.count(tj)) continue;
            if (!detail::atomsCompatFast(g1, qi, g2, tj, opts)) continue;

            // Check that all mapped neighbors of qi have compatible bonds
            bool bondOk = true;
            for (int qk : g1.neighbors[qi]) {
                auto it = mapping.find(qk);
                if (it == mapping.end()) continue;
                int tk = it->second;
                if (g2.bondOrder(tj, tk) == 0
                    || !ChemOps::bondsCompatible(g1, qi, qk, g2, tj, tk, opts)) {
                    bondOk = false;
                    break;
                }
            }
            if (bondOk) return false; // found extensible pair -> not maximal
        }
    }
    return true;
}

// ===========================================================================
// findNMCS -- N-molecule MCS with threshold via sequential pairwise reduction.
// Port of Java SearchEngine.findNMCS() / findNMCSMolecule()
// Sort by size, compute MCS(mol[0], mol[1]), extract common subgraph,
// compute MCS(result, mol[2]), etc.
// ===========================================================================

inline std::map<int,int> findNMCS(const std::vector<MolGraph>& molecules,
                                   const ChemOptions& opts,
                                   double threshold,
                                   int64_t timeoutMs) {
    if (molecules.size() < 2) return {};

    // Sort molecule indices by atom count (ascending)
    std::vector<int> sortedIdx(molecules.size());
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    std::sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b) {
        return molecules[a].n < molecules[b].n;
    });

    MCSOptions mcsOpts;
    mcsOpts.timeoutMs = timeoutMs;
    mcsOpts.connectedOnly = true;

    // Start with the smallest molecule, tracking original atom indices
    MolGraph currentMCS = molecules[sortedIdx[0]];
    std::vector<int> originalIndices(currentMCS.n);
    std::iota(originalIndices.begin(), originalIndices.end(), 0);

    for (size_t i = 1; i < sortedIdx.size(); ++i) {
        auto mcs = findMCS(currentMCS, molecules[sortedIdx[i]], opts, mcsOpts);
        if (mcs.empty()) return {};

        // Extract common subgraph using query-side atom indices
        std::vector<int> queryIndices;
        queryIndices.reserve(mcs.size());
        for (auto& [qi, tj] : mcs) queryIndices.push_back(qi);
        std::sort(queryIndices.begin(), queryIndices.end());

        // Update provenance: map surviving local indices → original indices
        std::vector<int> newOriginal;
        newOriginal.reserve(queryIndices.size());
        for (int localIdx : queryIndices)
            newOriginal.push_back(originalIndices[localIdx]);
        originalIndices = std::move(newOriginal);

        currentMCS = extractSubgraph(currentMCS, queryIndices);
        if (currentMCS.n == 0) return {};
    }

    // Threshold check: verify MCS is substructure of enough molecules
    if (threshold < 1.0 && currentMCS.n > 0) {
        int requiredCount = static_cast<int>(std::ceil(threshold * molecules.size()));
        int matchCount = 0;
        for (auto& mol : molecules) {
            if (isSubstructure(currentMCS, mol, opts, timeoutMs)) matchCount++;
        }
        if (matchCount < requiredCount) return {};
    }

    // Return mapping: original atom index in molecules[sortedIdx[0]] → NMCS position
    std::map<int,int> result;
    for (int i = 0; i < currentMCS.n; ++i)
        result[originalIndices[i]] = i;
    return result;
}

// ===========================================================================
// findScaffoldMCS -- MCS on Murcko scaffolds of two molecules.
// Port of Java SearchEngine.findScaffoldMCS()
// ===========================================================================

inline std::map<int,int> findScaffoldMCS(const MolGraph& g1, const MolGraph& g2,
                                          const ChemOptions& opts, const MCSOptions& mopts) {
    MolGraph s1 = murckoScaffold(g1);
    MolGraph s2 = murckoScaffold(g2);
    return findMCS(s1, s2, opts, mopts);
}

// ===========================================================================
// RGroupResult -- result of R-group decomposition for a single molecule.
// ===========================================================================

struct RGroupResult {
    MolGraph core;
    std::map<std::string, MolGraph> rgroups;
};

// ===========================================================================
// decomposeRGroups -- decompose molecules into core + R-groups.
// Port of Java SearchEngine.decomposeRGroups()
// For each molecule, find the core as a substructure, then BFS-flood
// connected components of non-core atoms as R1, R2, etc.
// ===========================================================================

inline std::vector<RGroupResult> decomposeRGroups(const MolGraph& core,
                                                    const std::vector<MolGraph>& molecules,
                                                    const ChemOptions& opts,
                                                    int64_t timeoutMs) {
    std::vector<RGroupResult> results;
    results.reserve(molecules.size());

    for (auto& mol : molecules) {
        RGroupResult decomp;

        // Find core as substructure of mol
        auto allMappings = findAllSubstructures(core, mol, opts, timeoutMs);
        if (allMappings.empty()) {
            results.push_back(std::move(decomp));
            continue;
        }

        // Use the first mapping (vector of pair<int,int>: query->target)
        auto& mapping = allMappings[0];
        std::unordered_set<int> coreAtomIndicesInMol;
        for (auto& [qi, tj] : mapping) {
            coreAtomIndicesInMol.insert(tj);
        }

        // Extract core subgraph from mol
        std::vector<int> coreIndices(coreAtomIndicesInMol.begin(), coreAtomIndicesInMol.end());
        std::sort(coreIndices.begin(), coreIndices.end());
        decomp.core = extractSubgraph(mol, coreIndices);

        // BFS to find connected R-group fragments
        int rGroupCounter = 1;
        std::unordered_set<int> processedNonCore;

        for (auto& [qi, molIdx] : mapping) {
            for (int neighborIdx : mol.neighbors[molIdx]) {
                if (coreAtomIndicesInMol.count(neighborIdx)) continue;
                if (processedNonCore.count(neighborIdx)) continue;

                // BFS flood to find this R-group's atoms
                std::vector<int> rGroupAtoms;
                std::deque<int> queue;
                queue.push_back(neighborIdx);
                std::unordered_set<int> visited;
                visited.insert(neighborIdx);

                while (!queue.empty()) {
                    int current = queue.front(); queue.pop_front();
                    if (coreAtomIndicesInMol.count(current)) continue;
                    rGroupAtoms.push_back(current);
                    processedNonCore.insert(current);

                    for (int nb : mol.neighbors[current]) {
                        if (!coreAtomIndicesInMol.count(nb) && !visited.count(nb)) {
                            visited.insert(nb);
                            queue.push_back(nb);
                        }
                    }
                }

                if (!rGroupAtoms.empty()) {
                    std::sort(rGroupAtoms.begin(), rGroupAtoms.end());
                    decomp.rgroups["R" + std::to_string(rGroupCounter++)] =
                        extractSubgraph(mol, rGroupAtoms);
                }
            }
        }

        results.push_back(std::move(decomp));
    }
    return results;
}

// ===========================================================================
// validateTautomerConsistency -- global proton-conservation check for MCS.
//
// For each tautomeric class in the mapping, the total mobilisable proton
// count (N-H + O-H + S-H) must be equal between query-side and target-side
// matched atoms.  Returns true if proton counts match (or if there are no
// tautomeric atoms in the mapping); false if any group has a mismatch,
// meaning the mapping would require impossible proton transfer.
// ===========================================================================

inline int mobilisableProtons(const MolGraph& g, int atom) {
    int z = g.atomicNum[atom];
    // Only N, O, S can carry mobilisable protons in tautomeric systems
    if (z != 7 && z != 8 && z != 16) return 0;
    // Sum bond orders to heavy-atom neighbors
    int bondSum = 0;
    for (int nb : g.neighbors[atom])
        bondSum += g.bondOrder(atom, nb);
    // Default valence for tautomeric heteroatoms
    int defaultVal = (z == 7) ? 3 : 2; // N=3, O=2, S=2
    int implicitH = defaultVal - bondSum - std::abs(g.formalCharge[atom]);
    return std::max(0, implicitH);
}

inline bool validateTautomerConsistency(
    const MolGraph& g1, const MolGraph& g2, const std::map<int,int>& mcs) {
    if (mcs.empty()) return true;
    // Ensure tautomer classes are computed
    if (g1.tautomerClass.empty() || g2.tautomerClass.empty()) return true;

    // Group matched atoms by tautomer class and sum protons per side
    std::unordered_map<int, int> querySideProtons, targetSideProtons;
    for (auto& [qi, ti] : mcs) {
        if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
        int tc1 = g1.tautomerClass[qi];
        int tc2 = g2.tautomerClass[ti];
        if (tc1 < 0 && tc2 < 0) continue;
        // Use query-side class as group key (or target-side if query has none)
        int groupKey = (tc1 >= 0) ? tc1 : tc2 + 100000;
        querySideProtons[groupKey]  += mobilisableProtons(g1, qi);
        targetSideProtons[groupKey] += mobilisableProtons(g2, ti);
    }
    // Check each group for proton conservation
    for (auto& [grp, qProtons] : querySideProtons) {
        auto it = targetSideProtons.find(grp);
        int tProtons = (it != targetSideProtons.end()) ? it->second : 0;
        if (qProtons != tProtons) return false;
    }
    return true;
}

// ===========================================================================
// mapReaction -- atom-atom mapping between reactants and products via dMCS.
// Port of Java SearchEngine.mapReaction()
// Computes disconnected MCS between reactant and product atoms.
// ===========================================================================

inline std::map<int,int> mapReaction(const MolGraph& reactants, const MolGraph& products,
                                      const ChemOptions& opts, int64_t timeoutMs) {
    MCSOptions mcsOpts;
    mcsOpts.disconnectedMCS = true;
    mcsOpts.connectedOnly = false;
    mcsOpts.timeoutMs = timeoutMs;
    return findMCS(reactants, products, opts, mcsOpts);
}

// ===========================================================================
// mcsToSmiles -- extract the MCS induced subgraph as a canonical SMILES.
// ===========================================================================

/**
 * Build the induced subgraph of @p g on the atoms in @p mapping (keys)
 * and return its canonical SMILES string.
 *
 * Preserves atomic number, formal charge, mass number, aromaticity,
 * ring membership, bond order, bond aromaticity, and bond ring flags.
 *
 * @param g       source molecule (query side of the MCS)
 * @param mapping MCS atom-index mapping (query -> target); keys select atoms
 * @return canonical SMILES of the MCS subgraph, or "" if mapping is empty
 */
inline std::string mcsToSmiles(const MolGraph& g, const std::map<int,int>& mapping) {
    if (mapping.empty()) return "";

    // Sorted MCS atom indices + old-to-new re-indexing.
    std::vector<int> mcsAtoms;
    mcsAtoms.reserve(mapping.size());
    for (auto& kv : mapping) mcsAtoms.push_back(kv.first);
    std::sort(mcsAtoms.begin(), mcsAtoms.end());

    int k = static_cast<int>(mcsAtoms.size());
    std::unordered_map<int,int> reindex;
    reindex.reserve(k * 2);
    for (int i = 0; i < k; i++) reindex[mcsAtoms[i]] = i;

    // Atom property vectors.
    std::vector<int>     subAtomicNum(k), subFormalCharge(k), subMassNumber(k);
    std::vector<uint8_t> subRing(k), subAromatic(k);

    for (int i = 0; i < k; i++) {
        int old = mcsAtoms[i];
        subAtomicNum[i]    = g.atomicNum[old];
        subFormalCharge[i] = old < static_cast<int>(g.formalCharge.size()) ? g.formalCharge[old] : 0;
        subMassNumber[i]   = old < static_cast<int>(g.massNumber.size()) ? g.massNumber[old] : 0;
        subRing[i]         = g.ring[old];
        subAromatic[i]     = g.aromatic[old];
    }

    // Build adjacency -- only bonds where both endpoints are in the MCS.
    std::vector<std::vector<int>>  subNeighbors(k);
    std::vector<std::vector<int>>  subBondOrders(k);
    std::vector<std::vector<bool>> subBondRings(k);
    std::vector<std::vector<bool>> subBondAroms(k);

    std::set<std::pair<int,int>> seen;
    for (int oldI : mcsAtoms) {
        int newI = reindex[oldI];
        for (int nb : g.neighbors[oldI]) {
            auto it = reindex.find(nb);
            if (it == reindex.end()) continue;
            int newJ = it->second;
            auto edge = std::make_pair(std::min(oldI, nb), std::max(oldI, nb));
            if (!seen.insert(edge).second) continue;

            int  ord   = g.bondOrder(oldI, nb);
            bool bRing = g.bondInRing(oldI, nb);
            bool bArom = g.bondAromatic(oldI, nb);

            subNeighbors[newI].push_back(newJ);
            subNeighbors[newJ].push_back(newI);
            subBondOrders[newI].push_back(ord);
            subBondOrders[newJ].push_back(ord);
            subBondRings[newI].push_back(bRing);
            subBondRings[newJ].push_back(bRing);
            subBondAroms[newI].push_back(bArom);
            subBondAroms[newJ].push_back(bArom);
        }
    }

    MolGraph sub = MolGraph::Builder()
        .atomCount(k)
        .atomicNumbers(std::move(subAtomicNum))
        .formalCharges(std::move(subFormalCharge))
        .massNumbers(std::move(subMassNumber))
        .ringFlags(std::move(subRing))
        .aromaticFlags(std::move(subAromatic))
        .setNeighbors(std::move(subNeighbors))
        .setBondOrders(std::move(subBondOrders))
        .bondRingFlags(std::move(subBondRings))
        .bondAromaticFlags(std::move(subBondAroms))
        .build();

    return toSMILES(sub);
}

/**
 * Convenience: compute MCS of two molecules and return the result as SMILES.
 *
 * @param g1   query molecule
 * @param g2   target molecule
 * @param chem chemical matching options
 * @param opts MCS options
 * @return canonical SMILES of the MCS, or "" if no common substructure
 */
inline std::string findMcsSmiles(const MolGraph& g1, const MolGraph& g2,
                                  const ChemOptions& chem, const MCSOptions& opts) {
    auto mapping = findMCS(g1, g2, chem, opts);
    return mcsToSmiles(g1, mapping);
}

// ===========================================================================
// MCSResult -- structured return type for convenience SMILES-to-MCS API
// ===========================================================================

/**
 * Structured result of a convenience MCS-from-SMILES computation.
 * Bundles the atom-index mapping, MCS size, Tanimoto-like overlap,
 * and the MCS subgraph as a canonical SMILES string.
 *
 * @since 6.3.0
 */
struct MCSResult {
    /// Atom-index mapping from query to target (empty if no common substructure).
    std::map<int,int> mapping;
    /// Number of matched atoms.
    int size = 0;
    /// Overlap coefficient: size / min(queryAtoms, targetAtoms).
    double overlapCoefficient = 0.0;
    /// Canonical SMILES of the MCS subgraph, or "" if empty.
    std::string mcsSmiles;
};

/**
 * All-in-one convenience: parse two SMILES strings, compute MCS, and return
 * a structured MCSResult with mapping, size, Tanimoto, and MCS SMILES.
 *
 * @param smi1  SMILES string for the first molecule
 * @param smi2  SMILES string for the second molecule
 * @param chem  chemical matching options
 * @param opts  MCS options (timeout, induced, connected, etc.)
 * @return structured MCS result
 * @since 6.3.0
 */
inline MCSResult findMCSFromSmiles(const std::string& smi1, const std::string& smi2,
                                    const ChemOptions& chem = ChemOptions(),
                                    const MCSOptions& opts = MCSOptions()) {
    if (smi1.empty() || smi2.empty())
        return MCSResult{{}, 0, 0.0, ""};
    MolGraph g1 = parseSMILES(smi1);
    MolGraph g2 = parseSMILES(smi2);
    auto mapping = findMCS(g1, g2, chem, opts);
    std::string smiOut = mcsToSmiles(g1, mapping);

    MCSResult result;
    result.mapping = std::move(mapping);
    result.size = static_cast<int>(result.mapping.size());
    int minN = std::min(g1.n, g2.n);
    result.overlapCoefficient = minN > 0 ? static_cast<double>(result.size) / minN : 0.0;
    result.mcsSmiles = std::move(smiOut);
    return result;
}

// ===========================================================================
// Reaction-Aware MCS Post-Filter (v6.4.0)
// ===========================================================================

namespace detail {

/// Rarity tier weight for an atomic number. Tier 3 (S,P,Se,B,Si,I)=3.0,
/// Tier 2 (N,O,F,Cl,Br)=1.5, Tier 1 (all other non-C)=1.0.
inline double rarityWeight(int atomicNum) {
    switch (atomicNum) {
        case 16: case 15: case 34: case 5: case 14: case 53: return 3.0;
        case 7: case 8: case 9: case 17: case 35:   return 1.5;
        default: return 1.0;
    }
}

/// Collect heteroatom types (non-C, non-H) from a molecule graph.
inline std::unordered_set<int> heteroAtomTypes(const MolGraph& g) {
    std::unordered_set<int> types;
    for (int i = 0; i < g.n; ++i) {
        int z = g.atomicNum[i];
        if (z != 6 && z != 1) types.insert(z);
    }
    return types;
}

/// Count connected components among mapped atoms in g1 (BFS).
inline int countComponents(const std::map<int,int>& m, const MolGraph& g1) {
    if (m.size() <= 1) return 1;
    std::unordered_set<int> mapped;
    for (const auto& p : m) mapped.insert(p.first);
    std::unordered_set<int> visited;
    int components = 0;
    for (int start : mapped) {
        if (visited.count(start)) continue;
        components++;
        std::deque<int> q;
        q.push_back(start);
        visited.insert(start);
        while (!q.empty()) {
            int u = q.front(); q.pop_front();
            for (int v : g1.neighbors[u]) {
                if (mapped.count(v) && !visited.count(v)) {
                    visited.insert(v);
                    q.push_back(v);
                }
            }
        }
    }
    return components;
}

/// Count bonds among mapped atoms in g1 (for reaction-aware scoring).
inline int countMappedBondsRA(const std::map<int,int>& m, const MolGraph& g1) {
    std::unordered_set<int> mapped;
    for (const auto& p : m) mapped.insert(p.first);
    int bonds = 0;
    for (int u : mapped) {
        for (int v : g1.neighbors[u]) {
            if (v > u && mapped.count(v)) bonds++;
        }
    }
    return bonds;
}

struct ScoredCandidate {
    std::map<int,int> mapping;
    double score;
    int bondCount;
    int insertionOrder;
};

} // namespace detail

/**
 * Score and rank MCS candidates by reaction relevance.
 *
 * Composite score = 0.40*S_size + 0.35*S_hetero + 0.15*S_rare + 0.10*S_conn
 *
 * @param candidates  list of MCS candidates (size K down to K-delta)
 * @param g1          first molecule graph
 * @param g2          second molecule graph
 * @return candidates sorted best-first
 * @since 6.4.0
 */
inline std::vector<std::map<int,int>> rankReactionAware(
    const std::vector<std::map<int,int>>& candidates,
    const MolGraph& g1, const MolGraph& g2) {

    if (candidates.empty()) return candidates;

    // Weights tuned for reaction mapping: heteroatom coverage is critical.
    // A 6-atom mapping with S should beat a 7-atom mapping without S.
    // PhD chemist recommendation: size weight reduced, hetero/rare increased.
    constexpr double W_SIZE   = 0.25;
    constexpr double W_HETERO = 0.40;
    constexpr double W_RARE   = 0.25;
    constexpr double W_CONN   = 0.10;
    constexpr double EPS      = 1e-9;

    // Determine K
    int K = 0;
    for (const auto& m : candidates)
        if (static_cast<int>(m.size()) > K) K = static_cast<int>(m.size());
    if (K == 0) return candidates;

    // H_universe
    auto hG1 = detail::heteroAtomTypes(g1);
    auto hG2 = detail::heteroAtomTypes(g2);
    std::unordered_set<int> hUniverse;
    for (int z : hG1) { if (hG2.count(z)) hUniverse.insert(z); }

    double rUniverse = 0.0;
    for (int z : hUniverse) rUniverse += detail::rarityWeight(z);

    // Score each
    std::vector<detail::ScoredCandidate> scored;
    scored.reserve(candidates.size());
    for (int idx = 0; idx < static_cast<int>(candidates.size()); ++idx) {
        const auto& m = candidates[idx];
        double sSize = static_cast<double>(m.size()) / K;

        // S_hetero and S_rare
        std::unordered_set<int> hMapped;
        for (const auto& p : m) {
            int z = g1.atomicNum[p.first];
            if (z != 6 && z != 1 && hUniverse.count(z)) hMapped.insert(z);
        }
        double sHetero = hUniverse.empty() ? 1.0
            : static_cast<double>(hMapped.size()) / hUniverse.size();

        double rMapped = 0.0;
        for (int z : hMapped) rMapped += detail::rarityWeight(z);
        double sRare = rUniverse <= 0.0 ? 1.0 : rMapped / rUniverse;

        double sConn = 1.0 / detail::countComponents(m, g1);

        double total = W_SIZE * sSize + W_HETERO * sHetero + W_RARE * sRare + W_CONN * sConn;
        int bondCount = detail::countMappedBondsRA(m, g1);
        scored.push_back({m, total, bondCount, idx});
    }

    std::sort(scored.begin(), scored.end(),
        [EPS](const detail::ScoredCandidate& a, const detail::ScoredCandidate& b) {
            double diff = b.score - a.score;
            if (std::abs(diff) > EPS) return diff < 0; // higher score first
            if (a.bondCount != b.bondCount) return a.bondCount > b.bondCount;
            return a.insertionOrder < b.insertionOrder;
        });

    std::vector<std::map<int,int>> result;
    result.reserve(scored.size());
    for (auto& sc : scored) result.push_back(std::move(sc.mapping));
    return result;
}

/**
 * Generate near-MCS candidates by systematic deletion from K-sized mappings
 * and greedy re-extension.
 *
 * @param g1            first molecule graph
 * @param g2            second molecule graph
 * @param chem          chemical matching options
 * @param exactMCS      the size-K MCS candidates from findAllMCS
 * @param delta         max size deficit (typically 2)
 * @param maxCandidates max candidates to generate
 * @return combined list of K, K-1, and K-2 sized candidates (deduplicated)
 * @since 6.4.0
 */
inline std::vector<std::map<int,int>> findNearMCS(
    const MolGraph& g1, const MolGraph& g2, const ChemOptions& chem,
    const std::vector<std::map<int,int>>& exactMCS,
    int delta, int maxCandidates) {

    if (exactMCS.empty()) return {};

    int K = static_cast<int>(exactMCS[0].size());
    int effectiveDelta = std::min(delta, K - 1);

    // Canonical SMILES dedup
    auto canonSmi = [&](const std::map<int,int>& m) -> std::string {
        std::vector<int> keys;
        keys.reserve(m.size());
        for (const auto& p : m) keys.push_back(p.first);
        return extractSubgraph(g1, keys).toCanonicalSmiles();
    };

    std::map<std::string, std::map<int,int>> seenBySmi;
    for (const auto& m : exactMCS) {
        if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
        auto smi = canonSmi(m);
        if (seenBySmi.find(smi) == seenBySmi.end())
            seenBySmi.emplace(std::move(smi), m);
    }

    // Deletion variants: remove one atom pair from each K-sized mapping
    std::vector<std::map<int,int>> deletionVariants;
    for (const auto& exact : exactMCS) {
        if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
        for (const auto& rem : exact) {
            if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
            std::map<int,int> variant;
            for (const auto& p : exact) {
                if (p.first != rem.first) variant[p.first] = p.second;
            }
            auto smi = canonSmi(variant);
            if (seenBySmi.find(smi) == seenBySmi.end()) {
                seenBySmi.emplace(std::move(smi), variant);
                deletionVariants.push_back(variant);
            }
        }
    }

    // Greedy re-extension: for each deletion variant, try adding a heteroatom
    // of a new type that is a neighbor of the mapped subgraph
    for (const auto& variant : deletionVariants) {
        if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;

        std::unordered_set<int> usedQ, usedT;
        std::unordered_set<int> mappedHetero;
        for (const auto& p : variant) {
            usedQ.insert(p.first);
            usedT.insert(p.second);
            int z = g1.atomicNum[p.first];
            if (z != 6 && z != 1) mappedHetero.insert(z);
        }

        bool found = false;
        for (const auto& p : variant) {
            if (found) break;
            for (int nbQ : g1.neighbors[p.first]) {
                if (found || static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
                if (usedQ.count(nbQ)) continue;
                int zQ = g1.atomicNum[nbQ];
                if (zQ == 6 || zQ == 1) continue;
                if (mappedHetero.count(zQ)) continue;

                for (const auto& tp : variant) {
                    if (found) break;
                    for (int nbT : g2.neighbors[tp.second]) {
                        if (usedT.count(nbT)) continue;
                        int zT = g2.atomicNum[nbT];
                        if (zT != zQ) continue;
                        if (chem.matchAtomType && g1.label[nbQ] != g2.label[nbT]) continue;
                        if (chem.matchFormalCharge && g1.formalCharge[nbQ] != g2.formalCharge[nbT]) continue;

                        std::map<int,int> extended = variant;
                        extended[nbQ] = nbT;
                        auto smi = canonSmi(extended);
                        if (seenBySmi.find(smi) == seenBySmi.end())
                            seenBySmi.emplace(std::move(smi), extended);
                        found = true;
                        break;
                    }
                }
            }
        }
    }

    // Level 2 deletions (K-2) if delta >= 2
    if (effectiveDelta >= 2) {
        std::vector<std::map<int,int>> kMinus1;
        for (const auto& kv : seenBySmi) {
            if (static_cast<int>(kv.second.size()) == K - 1)
                kMinus1.push_back(kv.second);
        }
        for (const auto& km1 : kMinus1) {
            if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
            for (const auto& rem : km1) {
                if (static_cast<int>(seenBySmi.size()) >= maxCandidates) break;
                std::map<int,int> variant;
                for (const auto& p : km1) {
                    if (p.first != rem.first) variant[p.first] = p.second;
                }
                auto smi = canonSmi(variant);
                if (seenBySmi.find(smi) == seenBySmi.end())
                    seenBySmi.emplace(std::move(smi), variant);
            }
        }
    }

    std::vector<std::map<int,int>> result;
    result.reserve(seenBySmi.size());
    for (auto& kv : seenBySmi) result.push_back(std::move(kv.second));
    return result;
}

/**
 * Heteroatom-seeded MCS: for a given element type, seed from every compatible
 * (query, target) heteroatom pair and greedily extend. Returns the best mapping
 * that includes at least one atom of the required element.
 *
 * This addresses the blind spot where findNearMCS (deletion + re-extension from
 * the size-K MCS) cannot reach subgraph selections that require a fundamentally
 * different seed -- e.g. the amino-acid+S chain in SAM vs. the amino-acid+ribose
 * chain that the standard MCS prefers.
 *
 * @param g1              first molecule graph
 * @param g2              second molecule graph
 * @param chem            chemical matching options
 * @param requiredElement atomic number that must appear in the mapping (e.g. 16 for S)
 * @return the best mapping seeded from that element, or empty if none found
 * @since 6.5.0
 */
inline std::map<int,int> heteroatomSeededMCS(
    const MolGraph& g1, const MolGraph& g2,
    const ChemOptions& chem, int requiredElement) {

    // Ensure lazy fields are initialised before accessing atomicNum/ring/degree
    g1.ensureCanonical(); g2.ensureCanonical();
    g1.ensureRingCounts(); g2.ensureRingCounts();

    // Collect all atom indices in g1 and g2 with the required element
    std::vector<int> qAtoms, tAtoms;
    for (int i = 0; i < g1.n; ++i)
        if (g1.atomicNum[i] == requiredElement) qAtoms.push_back(i);
    for (int j = 0; j < g2.n; ++j)
        if (g2.atomicNum[j] == requiredElement) tAtoms.push_back(j);
    if (qAtoms.empty() || tAtoms.empty()) return {};

    int n1 = g1.n, n2 = g2.n;
    std::map<int,int> bestMapping;
    int bestSize = 0;

    // For each compatible (qi, tj) seed pair, greedily extend
    for (int qi : qAtoms) {
        for (int tj : tAtoms) {
            if (!detail::atomsCompatFast(g1, qi, g2, tj, chem)) continue;

            // Initialize mapping arrays
            std::vector<int> q2t(n1, -1), t2q(n2, -1);
            q2t[qi] = tj;
            t2q[tj] = qi;
            int mapSize = 1;

            // Greedy bond-based extension: repeatedly find unmapped query atoms
            // adjacent to mapped atoms, and pair them with compatible unmapped
            // target atoms adjacent to the corresponding mapped target atom.
            bool progress = true;
            while (progress) {
                progress = false;
                for (int qk = 0; qk < n1; ++qk) {
                    if (q2t[qk] >= 0) continue;
                    int bestTk = -1, bestScore = -1;
                    for (int nb : g1.neighbors[qk]) {
                        if (q2t[nb] < 0) continue; // neighbor must be mapped
                        int tNb = q2t[nb];
                        for (int tk : g2.neighbors[tNb]) {
                            if (t2q[tk] >= 0) continue;
                            if (!detail::atomsCompatFast(g1, qk, g2, tk, chem)) continue;
                            if (!detail::bondsCompatible(g1, qk, nb, g2, tk, tNb, chem)) continue;
                            // Check consistency with all other mapped neighbors
                            bool consistent = true;
                            for (int qm : g1.neighbors[qk]) {
                                if (qm == nb || q2t[qm] < 0) continue;
                                int tm = q2t[qm];
                                int qOrd = g1.bondOrder(qk, qm);
                                int tOrd = g2.bondOrder(tk, tm);
                                if (qOrd != 0 && tOrd != 0) {
                                    if (!detail::bondsCompatible(g1, qk, qm, g2, tk, tm, chem))
                                    { consistent = false; break; }
                                }
                            }
                            if (!consistent) continue;
                            int score = (g2.ring[tk] && g1.ring[qk] ? 50 : 0)
                                      + std::min(g1.degree[qk], g2.degree[tk]);
                            if (score > bestScore) { bestScore = score; bestTk = tk; }
                        }
                    }
                    if (bestTk >= 0) {
                        q2t[qk] = bestTk;
                        t2q[bestTk] = qk;
                        mapSize++;
                        progress = true;
                    }
                }
            }

            if (mapSize > bestSize) {
                bestSize = mapSize;
                bestMapping.clear();
                for (int i = 0; i < n1; ++i)
                    if (q2t[i] >= 0) bestMapping[i] = q2t[i];
            }
        }
    }
    return bestMapping;
}

/**
 * Reaction-aware MCS: find MCS candidates, generate near-MCS variants,
 * add heteroatom-seeded candidates for shared rare elements, and re-rank
 * by composite scoring.
 *
 * @param g1   first molecule graph
 * @param g2   second molecule graph
 * @param chem chemical matching options
 * @param opts MCS options (reactionAware should be true)
 * @return the best-scoring candidate mapping
 * @since 6.4.0
 */
inline std::map<int,int> reactionAwareMCS(const MolGraph& g1, const MolGraph& g2,
                                           const ChemOptions& chem,
                                           const MCSOptions& opts) {
    if (g1.n == 0 || g2.n == 0) return {};

    // Ensure all lazy fields are initialised
    g1.ensureCanonical(); g2.ensureCanonical();
    g1.ensureRingCounts(); g2.ensureRingCounts();
    g1.getNLF1(); g2.getNLF1();

    // For reaction mapping, relax charge matching: reaction centers change
    // charge (e.g., S+ in SAM → S in homocysteine). Standard MCS with charge
    // matching off naturally includes the reaction-center heteroatom.
    ChemOptions reactionChem = chem;
    reactionChem.matchFormalCharge = false;

    int maxResults = std::max(5, opts.nearMcsCandidates / 4);
    MCSOptions safeOpts = opts;
    safeOpts.reactionAware = false; // prevent recursion
    auto exactCandidates = findAllMCS(g1, g2, reactionChem, safeOpts, maxResults);
    if (exactCandidates.empty()) return {};

    auto allCandidates = findNearMCS(g1, g2, reactionChem, exactCandidates,
                                     opts.nearMcsDelta, opts.nearMcsCandidates);
    if (allCandidates.empty()) allCandidates = exactCandidates;

    // Heteroatom-seeded candidates: for each heteroatom type shared between
    // both molecules, generate a candidate seeded from that element type.
    // This covers subgraph selections unreachable by deletion from the K-MCS.
    // Use relaxed ChemOptions for seeding: charge matching is off because
    // reaction centers often change charge (e.g., S+ in SAM → S in homocysteine).
    ChemOptions relaxedChem = chem;
    relaxedChem.matchFormalCharge = false;

    auto hG1 = detail::heteroAtomTypes(g1);
    auto hG2 = detail::heteroAtomTypes(g2);
    for (int z : hG1) {
        if (!hG2.count(z)) continue;
        auto seeded = heteroatomSeededMCS(g1, g2, relaxedChem, z);
        if (!seeded.empty()) {
            allCandidates.push_back(std::move(seeded));
        }
    }

    auto ranked = rankReactionAware(allCandidates, g1, g2);
    return ranked.empty() ? exactCandidates[0] : ranked[0];
}

/**
 * Convenience: parse two SMILES, run reaction-aware MCS, return mapping.
 *
 * @param smi1  SMILES for molecule 1
 * @param smi2  SMILES for molecule 2
 * @param chem  chemical matching options
 * @param opts  MCS options
 * @return best reaction-aware mapping
 * @since 6.4.0
 */
inline std::map<int,int> mapReactionAware(const std::string& smi1, const std::string& smi2,
                                           const ChemOptions& chem = ChemOptions(),
                                           MCSOptions opts = MCSOptions()) {
    if (smi1.empty() || smi2.empty()) return {};
    MolGraph g1 = parseSMILES(smi1);
    MolGraph g2 = parseSMILES(smi2);

    opts.disconnectedMCS = true;
    opts.connectedOnly = false;
    opts.reactionAware = true;

    return reactionAwareMCS(g1, g2, chem, opts);
}

// ===========================================================================
// Bond-change-aware MCS scoring (v6.5.3 — parity with Java BondChangeScorer)
// ===========================================================================

namespace detail {

/** Cost of breaking/forming a bond between atoms of given elements. */
inline double bondBreakCost(int z1, int z2) {
    if (z1 == 6 && z2 == 6) return 3.0;   // C-C: rare in biochemistry
    if ((z1 == 6 && (z2 == 7 || z2 == 8)) ||
        (z2 == 6 && (z1 == 7 || z1 == 8)))
        return 1.5;                          // C-N, C-O: amide/ester
    return 0.5;                              // heteroatom bonds
}

/** Compute bond-change penalty for a mapping. */
inline double bondChangePenalty(
    const std::map<int,int>& mapping,
    const MolGraph& g1, const MolGraph& g2)
{
    double penalty = 0.0;
    for (auto& [qi, ti] : mapping) {
        for (int qj : g1.neighbors[qi]) {
            if (qj <= qi) continue;
            auto it = mapping.find(qj);
            if (it == mapping.end()) continue;
            int tj = it->second;
            int boQ = g1.bondOrder(qi, qj);
            int boT = g2.bondOrder(ti, tj);
            if (boQ == boT && boT > 0) continue;
            if (boT == 0) penalty += bondBreakCost(g1.atomicNum[qi], g1.atomicNum[qj]);
            else penalty += 1.0; // bond order change
        }
    }
    // Check bonds in target not in query (bond formation)
    for (auto& [qi, ti] : mapping) {
        for (int tj : g2.neighbors[ti]) {
            // Find qj mapped to tj
            for (auto& [qj, mapped_tj] : mapping) {
                if (mapped_tj == tj && qj > qi) {
                    int boQ = g1.bondOrder(qi, qj);
                    if (boQ == 0) penalty += bondBreakCost(g1.atomicNum[qi], g1.atomicNum[qj]);
                    break;
                }
            }
        }
    }
    return penalty;
}

} // namespace detail

/**
 * Score an MCS mapping by bond-change plausibility.
 * Lower score = more chemically plausible mapping.
 * @since 6.5.3
 */
inline double bondChangeScore(
    const std::map<int,int>& mapping,
    const MolGraph& g1, const MolGraph& g2)
{
    return detail::bondChangePenalty(mapping, g1, g2);
}

// ===========================================================================
// Batch MCS with non-overlap constraints (v6.5.3 — parity with Java)
// ===========================================================================

/**
 * Find MCS for each query against ALL targets with non-overlapping
 * target atom constraints.
 *
 * For multi-component reaction mapping: N reactants → M products.
 * Each query finds its best MCS across all targets. Target atoms
 * used by earlier (larger) queries are excluded for later queries.
 *
 * @param queries  reactant fragments
 * @param targets  product fragment(s) — typically one combined product
 * @param chem     chemical matching options
 * @param opts     MCS options
 * @return vector of mappings, one per query (same order as input).
 *         Each mapping: {query_atom_idx: target_atom_idx}.
 * @since 6.6.0
 */
inline std::vector<std::map<int,int>> batchMcsConstrained(
    const std::vector<MolGraph>& queries,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts = MCSOptions())
{
    int nQ = static_cast<int>(queries.size());
    int nT = static_cast<int>(targets.size());

    std::vector<std::map<int,int>> results(nQ);

    // Sort queries by decreasing size (larger fragments first)
    std::vector<int> order(nQ);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return queries[a].n > queries[b].n;
    });

    // Track used target atoms per target
    std::vector<std::set<int>> usedTargetAtoms(nT);

    for (int qi : order) {
        std::map<int,int> bestMapping;
        int bestSize = 0;
        int bestTarget = -1;

        for (int ti = 0; ti < nT; ++ti) {
            const auto& used = usedTargetAtoms[ti];

            std::map<int,int> mapping;
            if (used.empty()) {
                // No atoms claimed — run MCS against full target
                mapping = findMCS(queries[qi], targets[ti], chem, opts);
            } else {
                // RE-RUN MCS against residual target (unclaimed atoms only)
                std::vector<int> residualAtoms;
                for (int i = 0; i < targets[ti].n; ++i) {
                    if (used.count(i) == 0) residualAtoms.push_back(i);
                }
                if (residualAtoms.empty()) continue;

                auto residual = extractSubgraph(targets[ti], residualAtoms);
                auto residualMcs = findMCS(queries[qi], residual, chem, opts);

                // Translate residual indices → original target indices
                for (auto& [qAtom, rAtom] : residualMcs) {
                    if (rAtom < static_cast<int>(residualAtoms.size()))
                        mapping[qAtom] = residualAtoms[rAtom];
                }
            }
            if (mapping.empty()) continue;

            if (static_cast<int>(mapping.size()) > bestSize) {
                bestSize = static_cast<int>(mapping.size());
                bestMapping = std::move(mapping);
                bestTarget = ti;
            }
        }

        if (bestTarget >= 0) {
            for (auto& [q, t] : bestMapping)
                usedTargetAtoms[bestTarget].insert(t);
        }
        results[qi] = std::move(bestMapping);
    }
    return results;
}

// ---------------------------------------------------------------------------
// MCS-cased aliases for bioinception API compatibility
// ---------------------------------------------------------------------------
inline int64_t resolveMCSTimeoutMs(const MolGraph& g1, const MolGraph& g2,
                                   const MCSOptions& opts) {
    return resolveMcsTimeoutMs(g1, g2, opts);
}

inline bool isValidMCSMapping(const MolGraph& g1, const MolGraph& g2,
                              const std::map<int,int>& mapping,
                              const ChemOptions& opts) {
    return isValidMcsMapping(g1, g2, mapping, opts);
}

inline std::map<int,int> recoverValidMCSMapping(const MolGraph& g1, const MolGraph& g2,
                                                const std::map<int,int>& raw,
                                                const ChemOptions& chem,
                                                const MCSOptions& opts) {
    return recoverValidMcsMapping(g1, g2, raw, chem, opts);
}

inline std::map<int,int> orientMCSResult(const std::map<int,int>& mapping, bool swapped) {
    return orientMcsResult(mapping, swapped);
}

inline std::map<int,int> runValidatedMCSDirection(const MolGraph& query, const MolGraph& target,
                                                  const ChemOptions& chem,
                                                  const MCSOptions& opts,
                                                  bool swappedBack) {
    return runValidatedMcsDirection(query, target, chem, opts, swappedBack);
}

inline std::string findMCSSmiles(const MolGraph& g1, const MolGraph& g2,
                                  const ChemOptions& chem, const MCSOptions& opts) {
    return findMcsSmiles(g1, g2, chem, opts);
}

inline std::vector<std::map<int,int>> batchMCSConstrained(
    const std::vector<MolGraph>& queries,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts = MCSOptions()) {
    return batchMcsConstrained(queries, targets, chem, opts);
}

} // namespace smsd

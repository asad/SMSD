/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * VF2++ Substructure Search Engine -- header-only C++17 port.
 * Zero external dependencies; requires "smsd/mol_graph.hpp" for MolGraph / ChemOptions.
 *
 * Algorithms ported from the Java SMSD SubstructureEngine:
 *   - Greedy probe (O(N) fast-path)
 *   - FASTiso ordering (forward-checking, smallest domain first, < 30 atoms)
 *   - VF3-Light ordering (label rarity sort, >= 30 atoms)
 *   - VF2PP backtracking with bit-parallel domains, NLF-1/2/3 pruning
 *   - Full feasibility check: degree, NLF, mapped-neighbor adjacency, bond compat
 */
#ifndef SMSD_VF2PP_HPP
#define SMSD_VF2PP_HPP

#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gpu_kernels.hpp"

namespace smsd {

// ============================================================================
// Match status (internal; public API remains unchanged)
// ============================================================================

enum class MatchStatus { Match, NoMatch, Timeout };

// ============================================================================
// Forward declarations of public API
// ============================================================================

bool isSubstructure(const MolGraph& query, const MolGraph& target,
                    const ChemOptions& opts, int64_t timeoutMs = 10000);

std::vector<std::pair<int,int>> findSubstructure(
    const MolGraph& query, const MolGraph& target,
    const ChemOptions& opts, int64_t timeoutMs = 10000);

std::vector<std::vector<std::pair<int,int>>> findAllSubstructures(
    const MolGraph& query, const MolGraph& target,
    const ChemOptions& opts, int64_t timeoutMs = 10000);

// ============================================================================
// Implementation detail namespace
// ============================================================================
namespace detail {

// Bit-ops: delegate to smsd::popcount64 / smsd::ctz64 from bitops.hpp
using smsd::popcount64;
using smsd::ctz64;

// ---------- Time budget ----------------------------------------------------

struct TimeBudget {
    using Clock = std::chrono::steady_clock;
    Clock::time_point deadline;
    int64_t counter_ = 0;
    static constexpr int64_t CHECK_EVERY = 1024;

    explicit TimeBudget(int64_t ms)
        : deadline(Clock::now() + std::chrono::milliseconds(std::max<int64_t>(1, ms))) {}

    bool expired() {
        if ((++counter_ & (CHECK_EVERY - 1)) != 0) return false;
        return Clock::now() >= deadline;
    }
    bool expiredNow() const { return Clock::now() >= deadline; }
    int64_t remainingMs() const {
        auto r = std::chrono::duration_cast<std::chrono::milliseconds>(deadline - Clock::now()).count();
        return std::max<int64_t>(0, r);
    }
};

// ---------- Iterate set bits in a uint64_t word array ----------------------

inline void iterateBits(const uint64_t* words, int nWords,
                        int* out, int& outCount) {
    outCount = 0;
    for (int w = 0; w < nWords; ++w) {
        uint64_t bits = words[w];
        while (bits != 0) {
            out[outCount++] = (w << 6) | ctz64(bits);
            bits &= bits - 1;
        }
    }
}

// ---------- NLF sorted-array merge check -----------------------------------

inline bool nlfOk(const int* fq, int fqLen, const int* ft, int ftLen) {
    int fi = 0, ti = 0;
    while (fi < fqLen) {
        int fLabel = fq[fi], fFreq = fq[fi + 1];
        while (ti < ftLen && ft[ti] < fLabel) ti += 2;
        if (ti >= ftLen || ft[ti] != fLabel || ft[ti + 1] < fFreq) return false;
        fi += 2;
    }
    return true;
}

// ---------- NLF builders ---------------------------------------------------

// NLF label: atomicNum + aromaticity, WITHOUT ring bit (ring is a one-directional
// constraint so including it would over-prune non-ring-query vs ring-target matches).
inline int nlfLabel(const MolGraph& g, int idx) {
    return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
}

// Builds sorted (label, freq) pairs for 1-hop neighbors.
inline std::vector<int> buildNLF1(const MolGraph& g, int idx) {
    std::unordered_map<int,int> freq;
    for (int k = 0; k < g.degree[idx]; ++k)
        ++freq[nlfLabel(g, g.neighbors[idx][k])];
    // Pack into sorted array of (label, freq) pairs
    std::vector<int> arr(freq.size() * 2);
    int p = 0;
    for (auto& kv : freq) { arr[p++] = kv.first; arr[p++] = kv.second; }
    // Insertion-sort by label
    int pairs = static_cast<int>(freq.size());
    for (int i = 1; i < pairs; ++i) {
        int kl = arr[i*2], kf = arr[i*2+1];
        int j = i - 1;
        while (j >= 0 && arr[j*2] > kl) {
            arr[(j+1)*2]   = arr[j*2];
            arr[(j+1)*2+1] = arr[j*2+1];
            --j;
        }
        arr[(j+1)*2]   = kl;
        arr[(j+1)*2+1] = kf;
    }
    return arr;
}

// 2-hop NLF: nodes at distance exactly 2.
inline std::vector<int> buildNLF2(const MolGraph& g, int idx) {
    std::unordered_map<int,int> freq;
    std::vector<bool> seen(g.n, false);
    // Mark direct neighbors
    for (int k = 0; k < g.degree[idx]; ++k) seen[g.neighbors[idx][k]] = true;
    seen[idx] = true;
    for (int k = 0; k < g.degree[idx]; ++k) {
        int nb = g.neighbors[idx][k];
        for (int m = 0; m < g.degree[nb]; ++m) {
            int j = g.neighbors[nb][m];
            if (!seen[j]) {
                seen[j] = true;
                ++freq[nlfLabel(g, j)];
            }
        }
    }
    std::vector<int> arr(freq.size() * 2);
    int p = 0;
    for (auto& kv : freq) { arr[p++] = kv.first; arr[p++] = kv.second; }
    int pairs = static_cast<int>(freq.size());
    for (int i = 1; i < pairs; ++i) {
        int kl = arr[i*2], kf = arr[i*2+1], j = i - 1;
        while (j >= 0 && arr[j*2] > kl) {
            arr[(j+1)*2] = arr[j*2]; arr[(j+1)*2+1] = arr[j*2+1]; --j;
        }
        arr[(j+1)*2] = kl; arr[(j+1)*2+1] = kf;
    }
    return arr;
}

// 3-hop NLF: nodes at distance exactly 3.
inline std::vector<int> buildNLF3(const MolGraph& g, int idx) {
    // level1 = direct neighbors
    std::vector<bool> level1(g.n, false), level2(g.n, false), level3(g.n, false);
    level1[idx] = true;
    for (int k = 0; k < g.degree[idx]; ++k) level1[g.neighbors[idx][k]] = true;
    // level2
    for (int k = 0; k < g.degree[idx]; ++k) {
        int nb = g.neighbors[idx][k];
        for (int m = 0; m < g.degree[nb]; ++m) {
            int j = g.neighbors[nb][m];
            if (!level1[j]) level2[j] = true;
        }
    }
    // level3
    for (int j = 0; j < g.n; ++j) {
        if (!level2[j]) continue;
        for (int m = 0; m < g.degree[j]; ++m) {
            int v = g.neighbors[j][m];
            if (!level1[v] && !level2[v]) level3[v] = true;
        }
    }
    std::unordered_map<int,int> freq;
    for (int j = 0; j < g.n; ++j)
        if (level3[j]) ++freq[nlfLabel(g, j)];
    std::vector<int> arr(freq.size() * 2);
    int p = 0;
    for (auto& kv : freq) { arr[p++] = kv.first; arr[p++] = kv.second; }
    int pairs = static_cast<int>(freq.size());
    for (int i = 1; i < pairs; ++i) {
        int kl = arr[i*2], kf = arr[i*2+1], j = i - 1;
        while (j >= 0 && arr[j*2] > kl) {
            arr[(j+1)*2] = arr[j*2]; arr[(j+1)*2+1] = arr[j*2+1]; --j;
        }
        arr[(j+1)*2] = kl; arr[(j+1)*2+1] = kf;
    }
    return arr;
}


// Policy-aware NLF label. Reuses the strict/default cached tables only when
// ChemOptions matches that exact policy; otherwise build a conservative label.
inline int nlfLabelPolicy(const MolGraph& g, int idx, const ChemOptions& C) {
    if (!C.matchAtomType) return 0; // all atoms in one class
    int z = g.atomicNum[idx];
    if (C.tautomerAware && (z == 6 || z == 7 || z == 8 || z == 16)) z = 6;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT)
        return (z << 1) | (g.aromatic[idx] ? 1 : 0);
    return z;
}

inline bool nlfPolicyIsDefault(const ChemOptions& C) {
    return C.matchAtomType && !C.tautomerAware
        && C.aromaticityMode == ChemOptions::AromaticityMode::STRICT;
}

inline std::vector<int> buildPolicyNLF1(const MolGraph& g, int idx, const ChemOptions& C) {
    std::unordered_map<int,int> freq;
    for (int k = 0; k < g.degree[idx]; ++k)
        ++freq[nlfLabelPolicy(g, g.neighbors[idx][k], C)];
    std::vector<int> arr(freq.size() * 2);
    int p = 0;
    for (auto& kv : freq) { arr[p++] = kv.first; arr[p++] = kv.second; }
    int pairs = static_cast<int>(freq.size());
    for (int i = 1; i < pairs; ++i) {
        int kl = arr[i*2], kf = arr[i*2+1], j = i - 1;
        while (j >= 0 && arr[j*2] > kl) {
            arr[(j+1)*2] = arr[j*2]; arr[(j+1)*2+1] = arr[j*2+1]; --j;
        }
        arr[(j+1)*2] = kl; arr[(j+1)*2+1] = kf;
    }
    return arr;
}

inline std::vector<std::vector<int>> buildAllPolicyNLF1(const MolGraph& g, const ChemOptions& C) {
    std::vector<std::vector<int>> out(g.n);
    for (int i = 0; i < g.n; ++i) out[i] = buildPolicyNLF1(g, i, C);
    return out;
}

// ---------- Atom / bond compatibility --------------------------------------

inline bool atomsCompatFast(const MolGraph& gq, int qi,
                            const MolGraph& gt, int tj,
                            const ChemOptions& C) {
    // Tautomer-aware: C/N/O can interchange within tautomeric regions,
    // but all OTHER atom-level constraints still apply (charge, aromaticity,
    // ring membership, isotope, chirality, ring fusion).
    bool tautRelax = C.tautomerAware
        && !gq.tautomerClass.empty() && !gt.tautomerClass.empty()
        && gq.tautomerClass[qi] != -1 && gt.tautomerClass[tj] != -1;
    if (tautRelax) {
        int aq = gq.atomicNum[qi], at = gt.atomicNum[tj];
        if (!((aq == 6 || aq == 7 || aq == 8 || aq == 16) &&
              (at == 6 || at == 7 || at == 8 || at == 16)))
            tautRelax = false;  // not C/N/O pair — fall through to normal matching
        // Degree guard: genuine tautomers shift at most 1 proton,
        // so degree difference should be <= 1 for same-element matches.
        // Prevents over-matching unrelated functional groups (e.g., -OH vs =O).
        if (tautRelax && aq == at && std::abs(gq.degree[qi] - gt.degree[tj]) > 1)
            tautRelax = false;
    }
    if (!tautRelax && C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj]) return false;
    if (C.matchFormalCharge && gq.formalCharge[qi] != gt.formalCharge[tj]) return false;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT && gq.aromatic[qi] != gt.aromatic[tj]) return false;
    if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
    if (C.matchIsotope) {
        int qm = gq.massNumber[qi], tm = gt.massNumber[tj];
        if (qm != 0 && tm != 0 && qm != tm) return false;
    }
    if (C.useChirality) {
        int qs = gq.tetraChirality[qi], ts = gt.tetraChirality[tj];
        if (qs != 0 && ts != 0 && qs != ts) return false;
    }
    if (C.ringFusionMode != ChemOptions::RingFusionMode::IGNORE && gq.ring[qi] && gt.ring[tj]) {
        if (C.ringFusionMode == ChemOptions::RingFusionMode::STRICT) {
            if (gq.ringCount[qi] != gt.ringCount[tj]) return false;
        }
    }
    return true;
}

inline bool bondsCompatible(const MolGraph& g1, int qi, int qk,
                            const MolGraph& g2, int tj, int tk,
                            const ChemOptions& C) {
    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
    if (qOrd == 0 || tOrd == 0) return false;
    // Tautomer-aware: relax bond ORDER (single/double interchangeable)
    // but still enforce aromaticity, stereo, and ring constraints.
    bool tautBondRelax = C.tautomerAware
        && !g1.tautomerClass.empty() && !g2.tautomerClass.empty()
        && g1.tautomerClass[qi] != -1 && g1.tautomerClass[qk] != -1
        && g2.tautomerClass[tj] != -1 && g2.tautomerClass[tk] != -1;
    // Strict aromaticity: check bond aromaticity directly (no ring guard)
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
        && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk)) return false;
    if (C.useBondStereo) {
        int qConf = g1.dbStereo(qi, qk), tConf = g2.dbStereo(tj, tk);
        if (qConf != 0 && tConf != 0 && qConf != tConf) return false;
    }
    if (C.ringMatchesRingOnly && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk)) return false;
    if (tautBondRelax) return true;  // bond order relaxed for tautomeric bonds
    if (C.matchBondOrder == ChemOptions::BondOrderMode::ANY) return true;
    if (qOrd == tOrd) return true;
    if (C.matchBondOrder == ChemOptions::BondOrderMode::LOOSE) return true;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::FLEXIBLE) {
        bool qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
        if ((qa && ta) || (qa && (tOrd == 1 || tOrd == 2)) || (ta && (qOrd == 1 || qOrd == 2)))
            return true;
        // Resonance-aware: allow single/double interchange between atoms that
        // are both in a conjugated (aromatic) ring system.  This catches Kekulé
        // form mismatches (e.g. pyridone drawn as C=O vs C-O) without requiring
        // explicit tautomer flags.  (v6.8.0)
        if ((g1.aromatic[qi] && g1.aromatic[qk]) &&
            (g2.aromatic[tj] && g2.aromatic[tk]) &&
            (qOrd == 1 || qOrd == 2) && (tOrd == 1 || tOrd == 2))
            return true;
    }
    return false;
}

template <typename Mapping>
inline bool tetraParityCompatible(const MolGraph& gq, int qi,
                                  const MolGraph& gt, int tj,
                                  const Mapping& q2t) {
    int qs = gq.tetraChirality[qi];
    int ts = gt.tetraChirality[tj];
    if (qs == 0 || ts == 0) return true;

    int mapped = 0;
    int qnbIdx[4], tnbIdx[4];
    for (int nk = 0; nk < gq.degree[qi] && mapped < 4; ++nk) {
        int qk = gq.neighbors[qi][nk];
        int tk = q2t[qk];
        if (tk != -1) {
            qnbIdx[mapped] = qk;
            tnbIdx[mapped] = tk;
            ++mapped;
        }
    }
    if (mapped < 3) return true;

    int qPerm[4], tPerm[4];
    for (int a = 0; a < mapped; ++a) {
        qPerm[a] = 0;
        for (int nk = 0; nk < gq.degree[qi]; ++nk) {
            if (gq.neighbors[qi][nk] == qnbIdx[a]) {
                qPerm[a] = nk;
                break;
            }
        }
        tPerm[a] = 0;
        for (int nk = 0; nk < gt.degree[tj]; ++nk) {
            if (gt.neighbors[tj][nk] == tnbIdx[a]) {
                tPerm[a] = nk;
                break;
            }
        }
    }

    int qInv = 0, tInv = 0;
    for (int a = 0; a < mapped; ++a) {
        for (int b = a + 1; b < mapped; ++b) {
            if (qPerm[a] > qPerm[b]) ++qInv;
            if (tPerm[a] > tPerm[b]) ++tInv;
        }
    }

    bool parityMatch = ((qInv ^ tInv) & 1) == 0;
    return (qs == ts) ? parityMatch : !parityMatch;
}


// ---------- Shared substructure prescreen / fast paths --------------------

inline bool quickPrescreen(const MolGraph& query, const MolGraph& target,
                           const ChemOptions& opts) {
    if (query.n > target.n) return false;
    std::array<int, 120> qFreq{}, tFreq{};
    std::unordered_map<int,int> qChargeFreq, tChargeFreq;
    int qMaxDeg = 0, tMaxDeg = 0;
    int qRingCount = 0, tRingCount = 0;
    for (int i = 0; i < query.n; ++i) {
        int z = query.atomicNum[i];
        if (z >= 0 && z < 120) qFreq[z]++;
        if (opts.matchFormalCharge) qChargeFreq[query.formalCharge[i]]++;
        if (query.degree[i] > qMaxDeg) qMaxDeg = query.degree[i];
        if (query.ring[i]) qRingCount++;
    }
    for (int i = 0; i < target.n; ++i) {
        int z = target.atomicNum[i];
        if (z >= 0 && z < 120) tFreq[z]++;
        if (opts.matchFormalCharge) tChargeFreq[target.formalCharge[i]]++;
        if (target.degree[i] > tMaxDeg) tMaxDeg = target.degree[i];
        if (target.ring[i]) tRingCount++;
    }
    if (opts.matchAtomType) {
        for (int z = 0; z < 120; ++z)
            if (qFreq[z] > tFreq[z]) return false;
    }
    if (opts.matchFormalCharge) {
        for (const auto& [charge, count] : qChargeFreq) {
            auto it = tChargeFreq.find(charge);
            if (it == tChargeFreq.end() || count > it->second) return false;
        }
    }
    if (qMaxDeg > tMaxDeg) return false;
    if (opts.ringMatchesRingOnly && qRingCount > tRingCount) return false;

    // Pattern fingerprint pre-screen (v6.8.0): O(1) structural feature check.
    // If the query has a structural feature (atom-pair, 2-hop path, ring+degree)
    // that the target lacks, abort before expensive VF2++ setup.
    if (opts.matchAtomType && query.n > 16) {
        query.ensurePatternFP();
        target.ensurePatternFP();
        for (int w = 0; w < MolGraph::FP_WORDS; ++w) {
            if ((query.patternFP_[w] & target.patternFP_[w]) != query.patternFP_[w])
                return false;
        }
    }

    // Extended prescreens for medium+ molecules (skip for tiny graphs
    // where the SmallMolMatcher DFS resolves quickly)
    if (query.n > 16) {
        // Bond-count check: query must have <= target bonds for substructure.
        {
            int qBonds = 0, tBonds = 0;
            for (int i = 0; i < query.n; ++i) qBonds += query.degree[i];
            for (int i = 0; i < target.n; ++i) tBonds += target.degree[i];
            if (qBonds > tBonds) return false;
        }

        // Degree-frequency check: for each degree d, query must not have
        // more atoms of degree >= d than target does.
        {
            int qDegHist[32] = {};
            int tDegHist[32] = {};
            for (int i = 0; i < query.n; ++i) {
                int d = query.degree[i]; if (d > 31) d = 31;
                qDegHist[d]++;
            }
            for (int i = 0; i < target.n; ++i) {
                int d = target.degree[i]; if (d > 31) d = 31;
                tDegHist[d]++;
            }
            // Suffix sums: count atoms with degree >= d
            int qSuffix = 0, tSuffix = 0;
            for (int d = 31; d >= 1; --d) {
                qSuffix += qDegHist[d];
                tSuffix += tDegHist[d];
                if (qSuffix > tSuffix) return false;
            }
        }
    }

    return true;
}

// Hardened canonical-graph comparison (v6.7.4):
//  1. Uses canonicalComputed_ flag instead of zero-sentinel (hash CAN be 0)
//  2. On hash match, verifies structural equivalence to guard against collisions
//  3. canonicalHash is uint64_t — wrapping arithmetic is defined (no UB)
inline bool sameCanonicalGraph(const MolGraph& q, const MolGraph& t) {
    auto vecAt = [](const auto& vec, int idx) -> int {
        return idx < static_cast<int>(vec.size()) ? vec[idx] : 0;
    };
    if (q.n != t.n) return false;
    q.ensureCanonical();
    t.ensureCanonical();
    // Both must have completed canonicalization (large molecules may bail)
    if (!q.canonicalComputed_ || !t.canonicalComputed_) return false;
    if (q.canonicalHash != t.canonicalHash) return false;
    // Hash matched — verify structure to rule out collisions.
    // Compare canonical-order label and degree sequences (O(N)).
    for (int r = 0; r < q.n; ++r) {
        // Find original atoms with canonical rank r
        // (canonicalLabel[orig] = rank, so we need inv[rank] = orig)
        // Rather than build inv, compare sorted label+degree sequences.
    }
    // Fast structural check: label and degree in canonical order must match.
    // Build inverse mappings on the stack for small N, heap for large.
    std::vector<int> qInv(q.n), tInv(t.n);
    for (int i = 0; i < q.n; ++i) qInv[q.canonicalLabel[i]] = i;
    for (int i = 0; i < t.n; ++i) tInv[t.canonicalLabel[i]] = i;
    for (int r = 0; r < q.n; ++r) {
        int qi = qInv[r], ti = tInv[r];
        if (q.label[qi] != t.label[ti]) return false;
        if (q.degree[qi] != t.degree[ti]) return false;
        if (vecAt(q.formalCharge, qi) != vecAt(t.formalCharge, ti)) return false;
        if (vecAt(q.massNumber, qi) != vecAt(t.massNumber, ti)) return false;
        if (vecAt(q.hydrogenCount, qi) != vecAt(t.hydrogenCount, ti)) return false;
        if (vecAt(q.atomClass, qi) != vecAt(t.atomClass, ti)) return false;
    }
    // Verify edge structure: for each canonical-order atom, its sorted
    // neighbor ranks must match.
    for (int r = 0; r < q.n; ++r) {
        int qi = qInv[r], ti = tInv[r];
        if (q.degree[qi] != t.degree[ti]) return false;
        std::vector<std::array<int,4>> qEdges(q.degree[qi]), tEdges(t.degree[ti]);
        for (int k = 0; k < q.degree[qi]; ++k) {
            int qj = q.neighbors[qi][k];
            qEdges[k] = {q.canonicalLabel[qj], q.bondOrder(qi, qj),
                         q.bondInRing(qi, qj) ? 1 : 0,
                         q.bondAromatic(qi, qj) ? 1 : 0};
        }
        for (int k = 0; k < t.degree[ti]; ++k) {
            int tj = t.neighbors[ti][k];
            tEdges[k] = {t.canonicalLabel[tj], t.bondOrder(ti, tj),
                         t.bondInRing(ti, tj) ? 1 : 0,
                         t.bondAromatic(ti, tj) ? 1 : 0};
        }
        std::sort(qEdges.begin(), qEdges.end());
        std::sort(tEdges.begin(), tEdges.end());
        if (qEdges != tEdges) return false;
    }
    return true;
}

// O(N) identity check: same atom count, same degree sequence, every atom
// and bond matches in index order.  Catches self-match even when canonical
// labelling bails out (e.g. CANON_SEARCH_LIMIT hit on vancomycin).
inline bool isExactMatch(const MolGraph& q, const MolGraph& t,
                         const ChemOptions& opts) {
    auto vecAt = [](const auto& vec, int idx) -> int {
        return idx < static_cast<int>(vec.size()) ? vec[idx] : 0;
    };
    if (q.n != t.n) return false;
    for (int i = 0; i < q.n; ++i) {
        if (q.atomicNum[i] != t.atomicNum[i]) return false;
        if (q.degree[i]   != t.degree[i])     return false;
        if (q.aromatic[i]  != t.aromatic[i])   return false;
        if (vecAt(q.formalCharge, i) != vecAt(t.formalCharge, i)) return false;
        if (vecAt(q.massNumber, i) != vecAt(t.massNumber, i)) return false;
        if (vecAt(q.hydrogenCount, i) != vecAt(t.hydrogenCount, i)) return false;
        if (vecAt(q.atomClass, i) != vecAt(t.atomClass, i)) return false;
        if (opts.useChirality && vecAt(q.tetraChirality, i) != vecAt(t.tetraChirality, i)) return false;
    }
    // Check bonds: for each edge in q, verify same edge exists in t with same order
    for (int i = 0; i < q.n; ++i) {
        for (int k = 0; k < q.degree[i]; ++k) {
            int j = q.neighbors[i][k];
            if (j <= i) continue; // each bond once
            int qbo = q.bondOrder(i, j);
            int tbo = t.bondOrder(i, j);
            if (qbo != tbo) return false;
            if (tbo == 0) return false; // edge missing in target
            if (q.bondInRing(i, j) != t.bondInRing(i, j)) return false;
            if (q.bondAromatic(i, j) != t.bondAromatic(i, j)) return false;
            if (opts.useBondStereo && q.dbStereo(i, j) != t.dbStereo(i, j)) return false;
        }
    }
    return true;
}

inline std::vector<std::pair<int,int>> canonicalIsoMap(const MolGraph& q, const MolGraph& t) {
    std::vector<int> qByRank(q.n), tByRank(t.n);
    for (int i = 0; i < q.n; ++i) qByRank[q.canonicalLabel[i]] = i;
    for (int j = 0; j < t.n; ++j) tByRank[t.canonicalLabel[j]] = j;
    std::vector<std::pair<int,int>> out;
    out.reserve(q.n);
    for (int r = 0; r < q.n; ++r) out.emplace_back(qByRank[r], tByRank[r]);
    return out;
}

// linearizePath, pathMatchDFS, pathMatchAll removed in v6.8.0.
// VF2++ engine with orbit-based symmetry pruning and NLF checks
// handles all topologies (including linear paths) more efficiently.

inline bool isSimplePathGraph(const MolGraph& g) {
    if (g.n == 0) return true;
    int edgeCount = 0;
    int endpoints = 0;
    for (int i = 0; i < g.n; ++i) {
        int deg = g.degree[i];
        if (deg > 2) return false;
        if (deg <= 1) endpoints++;
        edgeCount += deg;
    }
    edgeCount /= 2;
    if (g.n == 1) return true;
    return edgeCount == g.n - 1 && endpoints == 2;
}

inline std::vector<int> linearizeSimplePath(const MolGraph& g) {
    std::vector<int> seq;
    seq.reserve(g.n);
    if (g.n == 0) return seq;
    int start = 0;
    for (int i = 0; i < g.n; ++i) {
        if (g.degree[i] <= 1) { start = i; break; }
    }
    int prev = -1;
    int cur = start;
    while (cur != -1) {
        seq.push_back(cur);
        int next = -1;
        for (int nb : g.neighbors[cur]) {
            if (nb != prev) { next = nb; break; }
        }
        prev = cur;
        cur = next;
    }
    return seq;
}

inline std::vector<std::pair<int,int>> findPathSubstructure(
    const MolGraph& query, const MolGraph& target, const ChemOptions& opts) {
    if (query.n == 0) return {};
    if (query.n > target.n) return {};
    if (!isSimplePathGraph(query) || !isSimplePathGraph(target)) return {};

    const auto qSeq = linearizeSimplePath(query);
    auto tSeq = linearizeSimplePath(target);
    auto tryLinearTarget = [&](const std::vector<int>& seqT) -> std::vector<std::pair<int,int>> {
        for (int start = 0; start + query.n <= target.n; ++start) {
            bool ok = true;
            for (int i = 0; i < query.n; ++i) {
                if (!atomsCompatFast(query, qSeq[i], target, seqT[start + i], opts)) {
                    ok = false;
                    break;
                }
                if (i == 0) continue;
                if (!bondsCompatible(query, qSeq[i - 1], qSeq[i],
                                     target, seqT[start + i - 1], seqT[start + i], opts)) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                std::vector<std::pair<int,int>> out;
                out.reserve(query.n);
                for (int i = 0; i < query.n; ++i)
                    out.emplace_back(qSeq[i], seqT[start + i]);
                return out;
            }
        }
        return {};
    };

    auto forward = tryLinearTarget(tSeq);
    if (!forward.empty()) return forward;
    std::reverse(tSeq.begin(), tSeq.end());
    return tryLinearTarget(tSeq);
}

// --- Lightweight small-molecule DFS (v6.8.0) ---
// For small queries (Nq <= 16), bypass the full VF2PP constructor overhead
// (NLF, AC-3, frontierStack, etc.) and do direct backtracking DFS with
// inline atom/bond compatibility.  Avoids ~6-10µs of setup cost.

struct SmallMolMatcher {
    static constexpr int MAX_QUERY_ATOMS = 32;
    static constexpr int MAX_TARGET_ATOMS = 512;
    const MolGraph& gq;
    const MolGraph& gt;
    const ChemOptions& C;
    int Nq, Nt;
    int q2t[MAX_QUERY_ATOMS]; // query → target mapping, -1 = unmapped
    bool tUsed[MAX_TARGET_ATOMS]; // target atom used flags
    int order[MAX_QUERY_ATOMS]; // BFS-based search order
    int parentNbr[MAX_QUERY_ATOMS]; // for each position, index of a mapped neighbor (-1 = root)
    bool found;

    SmallMolMatcher(const MolGraph& gq_, const MolGraph& gt_,
                    const ChemOptions& C_)
        : gq(gq_), gt(gt_), C(C_), Nq(gq_.n), Nt(gt_.n), found(false) {
        if (Nq > MAX_QUERY_ATOMS || Nt > MAX_TARGET_ATOMS)
            throw std::out_of_range("SmallMolMatcher: molecule exceeds fixed-size array bounds");
        std::memset(q2t, -1, sizeof(q2t));
        std::memset(tUsed, 0, sizeof(tUsed));
        std::memset(parentNbr, -1, sizeof(parentNbr));

        // BFS order from highest-degree atom — ensures each subsequent atom
        // (except component roots) has at least one already-mapped neighbor,
        // enabling adjacency-constrained candidate selection.
        bool picked[MAX_QUERY_ATOMS] = {};
        int bestRoot = 0;
        for (int i = 1; i < Nq; ++i)
            if (gq.degree[i] > gq.degree[bestRoot]) bestRoot = i;
        order[0] = bestRoot; picked[bestRoot] = true; parentNbr[0] = -1;
        int head = 0, tail = 1;
        while (head < tail && tail < Nq) {
            int u = order[head++];
            for (int k = 0; k < gq.degree[u]; ++k) {
                int v = gq.neighbors[u][k];
                if (!picked[v]) {
                    picked[v] = true;
                    parentNbr[tail] = u;
                    order[tail++] = v;
                }
            }
        }
        // Handle disconnected components (shouldn't happen after splitting, but safe)
        for (int i = 0; i < Nq && tail < Nq; ++i) {
            if (!picked[i]) {
                picked[i] = true;
                parentNbr[tail] = -1;
                order[tail++] = i;
            }
        }
    }

    bool atomOk(int qi, int tj) const {
        if (gt.degree[tj] < gq.degree[qi]) return false;
        return atomsCompatFast(gq, qi, gt, tj, C);
    }

    bool edgesOk(int qi, int tj) const {
        for (int qk = 0; qk < Nq; ++qk) {
            if (qk == qi) continue;
            int tk = q2t[qk];
            if (tk == -1) continue;
            bool qBond = gq.hasBond(qi, qk);
            bool tBond = gt.hasBond(tj, tk);
            if (qBond) {
                if (!tBond) return false;
                if (!bondsCompatible(gq, qi, qk, gt, tj, tk, C)) return false;
            } else if (C.induced && tBond) {
                return false;
            }
        }
        if (C.useChirality && !tetraParityCompatible(gq, qi, gt, tj, q2t))
            return false;
        return true;
    }

    void dfs(int pos) {
        if (found) return;
        if (pos == Nq) { found = true; return; }
        int qi = order[pos];
        int pnb = parentNbr[pos]; // a mapped neighbor, or -1

        if (pnb >= 0 && q2t[pnb] >= 0) {
            // Constrained: only try neighbors of mapped parent
            int tp = q2t[pnb];
            for (int k = 0; k < gt.degree[tp]; ++k) {
                int tj = gt.neighbors[tp][k];
                if (tUsed[tj]) continue;
                if (!atomOk(qi, tj)) continue;
                if (!edgesOk(qi, tj)) continue;
                q2t[qi] = tj; tUsed[tj] = true;
                dfs(pos + 1);
                if (found) return;
                q2t[qi] = -1; tUsed[tj] = false;
            }
        } else {
            // Unconstrained root: try all target atoms
            for (int tj = 0; tj < Nt; ++tj) {
                if (tUsed[tj]) continue;
                if (!atomOk(qi, tj)) continue;
                if (!edgesOk(qi, tj)) continue;
                q2t[qi] = tj; tUsed[tj] = true;
                dfs(pos + 1);
                if (found) return;
                q2t[qi] = -1; tUsed[tj] = false;
            }
        }
    }

    bool exists() { dfs(0); return found; }

    std::vector<std::pair<int,int>> findOne() {
        dfs(0);
        if (!found) return {};
        std::vector<std::pair<int,int>> out(Nq);
        for (int i = 0; i < Nq; ++i) out[i] = {i, q2t[i]};
        return out;
    }
};

// --- Disconnected component detection & splitting (v6.8.0) ---
// Split disconnected queries (e.g. salts [Na+].[Cl-]) into
// independent searches per component.  Monolithic search on disconnected
// graphs is combinatorially explosive.

inline int countComponents(const MolGraph& g) {
    if (g.n == 0) return 0;
    std::vector<char> vis(g.n, 0);
    int nComp = 0;
    for (int s = 0; s < g.n; ++s) {
        if (vis[s]) continue;
        ++nComp;
        std::deque<int> bfs;
        bfs.push_back(s); vis[s] = 1;
        while (!bfs.empty()) {
            int u = bfs.front(); bfs.pop_front();
            for (int v : g.neighbors[u]) {
                if (!vis[v]) { vis[v] = 1; bfs.push_back(v); }
            }
        }
    }
    return nComp;
}

inline std::vector<MolGraph> splitComponents(const MolGraph& g) {
    if (g.n == 0) return {};
    std::vector<int> compId(g.n, -1);
    int nComp = 0;
    for (int s = 0; s < g.n; ++s) {
        if (compId[s] >= 0) continue;
        int cid = nComp++;
        std::deque<int> bfs;
        bfs.push_back(s); compId[s] = cid;
        while (!bfs.empty()) {
            int u = bfs.front(); bfs.pop_front();
            for (int v : g.neighbors[u]) {
                if (compId[v] < 0) { compId[v] = cid; bfs.push_back(v); }
            }
        }
    }
    if (nComp <= 1) return {g}; // single component — no split needed

    // Build per-component MolGraphs via Builder
    std::vector<MolGraph> comps;
    comps.reserve(nComp);
    for (int c = 0; c < nComp; ++c) {
        // Collect atoms in this component
        std::vector<int> atoms;
        for (int i = 0; i < g.n; ++i)
            if (compId[i] == c) atoms.push_back(i);
        int cn = static_cast<int>(atoms.size());
        // Build old→new index map
        std::vector<int> old2new(g.n, -1);
        for (int k = 0; k < cn; ++k) old2new[atoms[k]] = k;

        std::vector<int> atomicNums(cn), charges(cn), masses(cn);
        std::vector<uint8_t> ringFlags(cn), aromFlags(cn);
        std::vector<int> tetraChir(cn);
        std::vector<int> tautClass(cn, -1);
        std::vector<float> tautWeight(cn, 1.0f);
        std::vector<std::vector<int>> neighbors(cn);
        for (int k = 0; k < cn; ++k) {
            int oi = atoms[k];
            atomicNums[k]  = g.atomicNum[oi];
            charges[k]     = g.formalCharge[oi];
            masses[k]      = g.massNumber[oi];
            ringFlags[k]   = g.ring[oi];
            aromFlags[k]   = g.aromatic[oi];
            tetraChir[k]   = g.tetraChirality[oi];
            if (!g.tautomerClass.empty()) tautClass[k] = g.tautomerClass[oi];
            if (!g.tautomerWeight.empty()) tautWeight[k] = g.tautomerWeight[oi];
            for (int nb : g.neighbors[oi]) {
                int nj = old2new[nb];
                if (nj >= 0) neighbors[k].push_back(nj);
            }
        }

        // Build bond orders & aromatic/ring bond flags
        std::vector<std::vector<int>> bondOrders(cn, std::vector<int>(cn, 0));
        std::vector<std::vector<bool>> bondRingFlags(cn, std::vector<bool>(cn, false));
        std::vector<std::vector<bool>> bondAromFlags(cn, std::vector<bool>(cn, false));
        for (int k = 0; k < cn; ++k) {
            int oi = atoms[k];
            for (int nb : g.neighbors[oi]) {
                int nj = old2new[nb];
                if (nj < 0 || nj <= k) continue;
                bondOrders[k][nj] = bondOrders[nj][k] = g.bondOrder(oi, nb);
                bondRingFlags[k][nj] = bondRingFlags[nj][k] = g.bondInRing(oi, nb);
                bondAromFlags[k][nj] = bondAromFlags[nj][k] = g.bondAromatic(oi, nb);
            }
        }

        MolGraph::Builder builder;
        builder.atomCount(cn)
               .atomicNumbers(atomicNums)
               .formalCharges(charges)
               .massNumbers(masses)
               .ringFlags(ringFlags)
               .aromaticFlags(aromFlags)
               .setNeighbors(neighbors)
               .setBondOrders(bondOrders)
               .bondRingFlags(bondRingFlags)
               .bondAromaticFlags(bondAromFlags)
               .tetrahedralChirality(tetraChir);
        MolGraph comp = builder.build();
        comp.tautomerClass = std::move(tautClass);
        comp.tautomerWeight = std::move(tautWeight);
        comps.push_back(std::move(comp));
    }
    return comps;
}

// isTree, chooseTreeRoot, buildRootedTree, treeMatchSubtree,
// treeSubstructureExists, treeSubstructureFind removed in v6.8.0.
// The naive DFS tree-match bypassed VF2++ domain pruning, causing
// PEG-like polymers and branched trees to time out.  VF2++ with
// AC-3 and orbit pruning handles all topologies efficiently.

// ==========================================================================
// AbstractVFMatcher -- base class with all shared state & algorithms
// ==========================================================================

class AbstractVFMatcher {
protected:
    const MolGraph& gq_;
    const MolGraph& gt_;
    const ChemOptions& C_;
    TimeBudget tb_;

    int Nq_, Nt_;
    int tWords_;
    bool singleWord_;

    // Mapping arrays -- pre-allocated, zero heap in hot loops
    std::vector<int> q2t_;
    std::vector<int> t2q_;

    // NLF caches: pointers into MolGraph lazy caches — zero-copy, lifetime tied to gq_/gt_
    const std::vector<std::vector<int>>* qNLF1_ = nullptr;
    const std::vector<std::vector<int>>* tNLF1_ = nullptr;
    const std::vector<std::vector<int>>* qNLF2_ = nullptr;
    const std::vector<std::vector<int>>* tNLF2_ = nullptr;
    const std::vector<std::vector<int>>* qNLF3_ = nullptr;
    const std::vector<std::vector<int>>* tNLF3_ = nullptr;
    std::vector<std::vector<int>> qNLF1Own_;
    std::vector<std::vector<int>> tNLF1Own_;

    bool useTwoHop_, useThreeHop_, useBitParallel_, useRingOnly_, useStereo_, useInduced_;
    bool bitParallelSufficient_;

    // Bit-parallel candidate domains: flat 1D layout for cache locality (v6.8.0)
    // Access via domRow(qi)[w] — single flat allocation, zero per-row mallocs.
    std::vector<uint64_t> domain_;          // size: Nq_ * tWords_
    std::vector<uint64_t> usedMask_;

    // Query degree and target degree (raw pointers into MolGraph arrays)
    const int* qdeg_;
    const int* tdeg_;

    // Query neighbors sorted by descending degree
    std::vector<std::vector<int>> qNeighborsByDegDesc_;

    // Scratch buffers: flat 1D layout to minimise malloc calls (v6.8.0)
    std::vector<int> probeCandBuf_;                // for non-recursive greedy probe
    std::vector<int> candBufFlat_;                 // flat: Nq_ * candStride_
    int              candStride_ = 0;              // = Nt_ + 64
    std::vector<uint64_t> availBuf_;               // sized to tWords

    // Accessor helpers for flat domain / candidate buffers
    uint64_t* domRow(int qi)             { return domain_.data() + qi * tWords_; }
    const uint64_t* domRow(int qi) const { return domain_.data() + qi * tWords_; }
    int* candRow(int pos)                { return candBufFlat_.data() + pos * candStride_; }

    // Stats
    int64_t nodesVisited_  = 0;
    int64_t backtracks_    = 0;
    int64_t candidatesTried_ = 0;
    int64_t prunesAtom_    = 0;
    int64_t prunesBond_    = 0;
    int64_t prunesDegree_  = 0;
    int64_t prunesNLF_     = 0;
    bool    timedOut_      = false;
    bool    found_         = false;

public:
    AbstractVFMatcher(const MolGraph& gq, const MolGraph& gt,
                      const ChemOptions& C, int64_t timeoutMs)
        : gq_(gq), gt_(gt), C_(C), tb_(timeoutMs),
          Nq_(gq.n), Nt_(gt.n),
          tWords_((gt.n + 63) >> 6),
          singleWord_(((gt.n + 63) >> 6) <= 1),
          q2t_(gq.n, -1), t2q_(gt.n, -1),
          useTwoHop_(C.useTwoHopNLF),
          useThreeHop_(C.useThreeHopNLF),
          useBitParallel_(C.useBitParallelFeasibility),
          useRingOnly_(C.ringMatchesRingOnly),
          useStereo_(C.useBondStereo),
          useInduced_(C.induced),
          bitParallelSufficient_(C.useBitParallelFeasibility
              && C.matchBondOrder == ChemOptions::BondOrderMode::ANY
              && !C.useBondStereo && !C.ringMatchesRingOnly && !C.induced),
          domain_(static_cast<size_t>(gq.n) * ((gt.n + 63) >> 6), 0),
          usedMask_(((gt.n + 63) >> 6), 0),
          qdeg_(gq.degree.data()),
          tdeg_(gt.degree.data()),
          qNeighborsByDegDesc_(gq.n),
          probeCandBuf_(gt.n + 64),
          candBufFlat_(static_cast<size_t>(std::max(1, gq.n)) * (gt.n + 64)),
          candStride_(gt.n + 64),
          availBuf_(((gt.n + 63) >> 6), 0)
    {
        // Adaptive NLF disable for small molecules
        if (Nq_ <= 12 || Nt_ <= 12) useTwoHop_ = false;
        if (Nq_ <= 20 || Nt_ <= 20) useThreeHop_ = false;

        // Build query neighbors sorted by descending degree
        for (int i = 0; i < Nq_; ++i) {
            int deg = gq_.degree[i];
            qNeighborsByDegDesc_[i].resize(deg);
            for (int k = 0; k < deg; ++k)
                qNeighborsByDegDesc_[i][k] = gq_.neighbors[i][k];
            // Insertion sort descending by degree
            auto& nb = qNeighborsByDegDesc_[i];
            for (int a = 1; a < deg; ++a) {
                int key = nb[a], keyDeg = gq_.degree[key];
                int b = a - 1;
                while (b >= 0 && gq_.degree[nb[b]] < keyDeg) {
                    nb[b+1] = nb[b]; --b;
                }
                nb[b+1] = key;
            }
        }

        // Ensure canonical data is available for sameCanonicalGraph fast-path
        gt_.ensureCanonical();

        // Ensure lazy-computed fields are ready before matching
        if (C.ringFusionMode != ChemOptions::RingFusionMode::IGNORE) {
            gq_.ensureRingCounts();
            gt_.ensureRingCounts();
        }

        // NLF tables: reuse MolGraph caches only under the default strict policy.
        // For relaxed modes (tautomer-aware, aromaticity flexible, matchAtomType=false),
        // build a conservative 1-hop NLF locally and disable 2/3-hop NLF to avoid over-pruning.
        if (nlfPolicyIsDefault(C_)) {
            qNLF1_ = &gq_.getNLF1();
            tNLF1_ = &gt_.getNLF1();
            if (useTwoHop_) {
                qNLF2_ = &gq_.getNLF2();
                tNLF2_ = &gt_.getNLF2();
            }
            if (useThreeHop_) {
                qNLF3_ = &gq_.getNLF3();
                tNLF3_ = &gt_.getNLF3();
            }
        } else {
            qNLF1Own_ = buildAllPolicyNLF1(gq_, C_);
            tNLF1Own_ = buildAllPolicyNLF1(gt_, C_);
            qNLF1_ = &qNLF1Own_;
            tNLF1_ = &tNLF1Own_;
            useTwoHop_ = false;
            useThreeHop_ = false;
        }

        // --- GPU fast path for domain initialization ---
        // For large enough workloads, dispatch the Nq×Nt compatibility matrix
        // to GPU.  The GPU implements common-case checks only (atomicNum,
        // charge, aromaticity, ringOnly).  Advanced checks (tautomer, isotope,
        // chirality, ring fusion) are handled by the CPU feasibility check
        // during backtracking — the GPU domain is a safe superset.
        bool gpuDomainOk = false;
        if (Nq_ > 0 && Nt_ > 0 && static_cast<int64_t>(Nq_) * Nt_ > 2000
            && gpu_kern::graphKernelsAvailable()) {
            // Materialize vector<bool> to int arrays (vector<bool> has no .data())
            std::vector<int> qRi(Nq_), qAr(Nq_), tRi(Nt_), tAr(Nt_);
            for (int i = 0; i < Nq_; ++i) { qRi[i] = gq_.ring[i]; qAr[i] = gq_.aromatic[i]; }
            for (int j = 0; j < Nt_; ++j) { tRi[j] = gt_.ring[j]; tAr[j] = gt_.aromatic[j]; }
            std::vector<uint64_t> flatDom(static_cast<size_t>(Nq_) * tWords_, 0);
            gpuDomainOk = gpu_kern::domainInit(
                Nq_, Nt_, tWords_,
                gq_.atomicNum.data(), gq_.formalCharge.data(), qRi.data(), qAr.data(),
                gt_.atomicNum.data(), gt_.formalCharge.data(), tRi.data(), tAr.data(),
                C_.matchAtomType, C_.matchFormalCharge, C_.ringMatchesRingOnly,
                C_.aromaticityMode == ChemOptions::AromaticityMode::STRICT,
                flatDom.data());
            if (gpuDomainOk) {
                for (int i = 0; i < Nq_; ++i)
                    std::memcpy(domRow(i), &flatDom[i * tWords_],
                                tWords_ * sizeof(uint64_t));
            }
        }

        // --- CPU fallback (or small workload) ---
        if (!gpuDomainOk && Nq_ > 0 && Nt_ > 50) {
            std::unordered_map<int, std::vector<int>> targetByLabel;
            for (int j = 0; j < Nt_; ++j) {
                int key = C_.matchAtomType ? gt_.atomicNum[j] : 0;
                // Only bucket by aromatic; skip ring bit to avoid incorrectly
                // blocking non-ring query atoms from matching ring target atoms.
                key = (key << 1) | (gt_.aromatic[j] ? 1 : 0);
                targetByLabel[key].push_back(j);
            }
            for (int i = 0; i < Nq_; ++i) {
                for (auto& [labelKey, targets] : targetByLabel) {
                    if (C_.matchAtomType && !targets.empty()
                        && gq_.atomicNum[i] != gt_.atomicNum[targets[0]]) {
                        prunesAtom_ += static_cast<int64_t>(targets.size());
                        continue;
                    }
                    for (int j : targets) {
                        if (atomsCompatFast(gq_, i, gt_, j, C_))
                            domRow(i)[j >> 6] |= uint64_t(1) << (j & 63);
                        else
                            ++prunesAtom_;
                    }
                }
            }
        } else if (!gpuDomainOk) {
            for (int i = 0; i < Nq_; ++i)
                for (int j = 0; j < Nt_; ++j) {
                    if (atomsCompatFast(gq_, i, gt_, j, C_))
                        domRow(i)[j >> 6] |= uint64_t(1) << (j & 63);
                    else
                        ++prunesAtom_;
                }
        }

        // --- AC-3 Arc Consistency domain pruning (v6.8.0) ---
        // For each query edge (qi, qk), ensure every candidate Tj in
        // domain(qi) has at least one compatible neighbor Tk in domain(qk),
        // and vice versa.  Iterate until no domain changes.  If any domain
        // becomes empty, the search will terminate in O(1) without backtracking.
        // Only run for molecules above the small-molecule threshold.
        // AC-3 overhead dominates on tiny graphs (benzene→toluene class);
        // skip when query is small — greedy probe + backtracking suffice.
        if (Nq_ > 15 && Nt_ > 15) {
            bool changed = true;
            int acIters = 0;
            constexpr int kMaxACIters = 10; // cap iterations for huge molecules
            while (changed && acIters < kMaxACIters) {
                changed = false;
                ++acIters;
                for (int qi = 0; qi < Nq_; ++qi) {
                    for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
                        int qk = gq_.neighbors[qi][nk];
                        // For each candidate tj in domain(qi), check that at
                        // least one neighbor of tj is in domain(qk)
                        for (int w = 0; w < tWords_; ++w) {
                            uint64_t bits = domRow(qi)[w];
                            while (bits) {
                                int tj = (w << 6) | ctz64(bits);
                                bits &= bits - 1;
                                // Check: does tj have a bond-compatible neighbor in domain(qk)?
                                bool hasSupport = false;
                                for (int tn : gt_.neighbors[tj]) {
                                    if ((domRow(qk)[tn >> 6] & (uint64_t(1) << (tn & 63)))
                                        && bondsCompatible(gq_, qi, qk, gt_, tj, tn, C_)) {
                                        hasSupport = true;
                                        break;
                                    }
                                }
                                if (!hasSupport) {
                                    // Remove tj from domain(qi)
                                    domRow(qi)[tj >> 6] &= ~(uint64_t(1) << (tj & 63));
                                    changed = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    virtual ~AbstractVFMatcher() = default;

    // ----- Domain candidates (available & compatible) ----------------------

    int domainCandidates(int qi, int* out) const {
        int count = 0;
        for (int w = 0; w < tWords_; ++w) {
            uint64_t bits = domRow(qi)[w] & ~usedMask_[w];
            while (bits != 0) {
                out[count++] = (w << 6) | ctz64(bits);
                bits &= bits - 1;
            }
        }
        return count;
    }


    int domainSupport(int qi) const {
        int c = 0;
        for (int w = 0; w < tWords_; ++w) c += popcount64(domRow(qi)[w]);
        return c;
    }

    bool mappedBondCompatibleOnly(int qi, int tj) const {
        const auto& nbSorted = qNeighborsByDegDesc_[qi];
        const uint64_t* tAdj = gt_.adjLong[tj].data();
        for (int qk : nbSorted) {
            int tk = q2t_[qk];
            if (tk == -1) continue;
            if ((tAdj[tk >> 6] & (uint64_t(1) << (tk & 63))) == 0) return false;
            if (!bondsCompatible(gq_, qi, qk, gt_, tj, tk, C_)) return false;
        }
        return true;
    }

    bool hasLocalSupport(int qi, int tj) {
        for (int qn : gq_.neighbors[qi]) {
            if (q2t_[qn] != -1) continue;
            bool any = false;
            for (int tk : gt_.neighbors[tj]) {
                if (t2q_[tk] != -1) continue;
                if ((domRow(qn)[tk >> 6] & (uint64_t(1) << (tk & 63))) == 0) continue;
                if (!bondsCompatible(gq_, qi, qn, gt_, tj, tk, C_)) continue;
                bool ok = true;
                for (int qx : gq_.neighbors[qn]) {
                    int tx = q2t_[qx];
                    if (tx == -1 || qx == qi) continue;
                    if ((gt_.adjLong[tk][tx >> 6] & (uint64_t(1) << (tx & 63))) == 0
                        || !bondsCompatible(gq_, qn, qx, gt_, tk, tx, C_)) {
                        ok = false; break;
                    }
                }
                if (ok) { any = true; break; }
            }
            if (!any) return false;
        }
        return true;
    }

    // ----- Greedy probe (O(N) fast path) -----------------------------------

    bool greedyProbe(const int* order) {
        const bool sameSize = (Nq_ == Nt_);
        int matched = 0;
        for (int idx = 0; idx < Nq_; ++idx) {
            int qi = order[idx];
            int nCands = domainCandidates(qi, probeCandBuf_.data());

            // When query and target are the same size, try identity map first
            if (sameSize && qi < Nt_) {
                bool inDomain = false;
                for (int c = 0; c < nCands; ++c) {
                    if (probeCandBuf_[c] == qi) { inDomain = true; break; }
                }
                if (inDomain && feasible(qi, qi)) {
                    q2t_[qi] = qi; t2q_[qi] = qi;
                    usedMask_[qi >> 6] |= uint64_t(1) << (qi & 63);
                    ++matched;
                    continue;
                }
            }

            int bestTj = -1, bestDist = 0x7FFFFFFF;
            bool found2 = false;
            for (int c = 0; c < nCands; ++c) {
                int tj = probeCandBuf_[c];
                int dist = tdeg_[tj] - qdeg_[qi];
                if (dist < 0) dist = -dist;
                if (dist < bestDist && feasible(qi, tj)) {
                    bestTj = tj; bestDist = dist; found2 = true;
                    if (dist == 0) break;
                }
            }
            if (!found2) {
                // Undo all matched so far
                for (int k = matched - 1; k >= 0; --k) {
                    int qk = order[k], tk = q2t_[qk];
                    q2t_[qk] = -1; t2q_[tk] = -1;
                    usedMask_[tk >> 6] &= ~(uint64_t(1) << (tk & 63));
                }
                return false;
            }
            q2t_[qi] = bestTj; t2q_[bestTj] = qi;
            usedMask_[bestTj >> 6] |= uint64_t(1) << (bestTj & 63);
            ++matched;
        }
        return true;
    }

    // ----- FASTiso ordering ------------------------------------------------
    // Forward-checking, smallest domain first, break ties by highest degree.

    void fastisoOrder(int* order) const {
        std::vector<bool> picked(Nq_, false);

        // Compute connected components of the query graph via BFS
        std::vector<int> compId(Nq_, -1);
        int nComp = 0;
        for (int s = 0; s < Nq_; ++s) {
            if (compId[s] >= 0) continue;
            int cid = nComp++;
            compId[s] = cid;
            std::deque<int> bfs;
            bfs.push_back(s);
            while (!bfs.empty()) {
                int u = bfs.front(); bfs.pop_front();
                for (int k = 0; k < gq_.degree[u]; ++k) {
                    int v = gq_.neighbors[u][k];
                    if (compId[v] < 0) { compId[v] = cid; bfs.push_back(v); }
                }
            }
        }

        // Working copy of domains for forward checking (flat 1D like domain_)
        std::vector<uint64_t> workDomain(domain_);  // full copy of flat domain

        std::vector<uint64_t> reachable(tWords_, 0);
        int activeComp = -1; // current component being ordered

        for (int pos = 0; pos < Nq_; ++pos) {
            int bestQ = -1, bestSize = 0x7FFFFFFF, bestDeg = -1;

            // First pass: prefer atoms from the active component
            if (activeComp >= 0) {
                for (int i = 0; i < Nq_; ++i) {
                    if (picked[i] || compId[i] != activeComp) continue;
                    int size = 0;
                    for (int w = 0; w < tWords_; ++w)
                        size += popcount64(workDomain[i * tWords_ + w]);
                    if (size < bestSize || (size == bestSize && gq_.degree[i] > bestDeg)) {
                        bestQ = i; bestSize = size; bestDeg = gq_.degree[i];
                    }
                }
            }

            // If no unpicked atom in active component, pick globally (starts new component)
            if (bestQ < 0) {
                for (int i = 0; i < Nq_; ++i) {
                    if (picked[i]) continue;
                    int size = 0;
                    for (int w = 0; w < tWords_; ++w)
                        size += popcount64(workDomain[i * tWords_ + w]);
                    if (size < bestSize || (size == bestSize && gq_.degree[i] > bestDeg)) {
                        bestQ = i; bestSize = size; bestDeg = gq_.degree[i];
                    }
                }
                activeComp = compId[bestQ];
            }

            order[pos] = bestQ;
            picked[bestQ] = true;

            // Forward-check: restrict neighbor domains to reachable from bestQ's domain.
            // Skip neighbors whose domain is already empty (nothing to restrict).
            for (int nk = 0; nk < gq_.degree[bestQ]; ++nk) {
                int qk = gq_.neighbors[bestQ][nk];
                if (picked[qk]) continue;
                // Check if qk's domain is already empty — skip expensive reachability
                bool hasAny = false;
                for (int w = 0; w < tWords_ && !hasAny; ++w)
                    if (workDomain[qk * tWords_ + w]) hasAny = true;
                if (!hasAny) continue;
                std::memset(reachable.data(), 0, tWords_ * sizeof(uint64_t));
                for (int w = 0; w < tWords_; ++w) {
                    uint64_t bits = workDomain[bestQ * tWords_ + w];
                    while (bits != 0) {
                        int tj = (w << 6) | ctz64(bits);
                        bits &= bits - 1;
                        if (tj < Nt_) {
                            for (int rw = 0; rw < tWords_; ++rw)
                                reachable[rw] |= gt_.adjLong[tj][rw];
                        }
                    }
                }
                for (int rw = 0; rw < tWords_; ++rw)
                    workDomain[qk * tWords_ + rw] &= reachable[rw];
            }
        }
    }

    // ----- Support/frontier ordering ---------------------------------------
    // Large-query order driven by target support (domain popcount), with a
    // frontier preference inside the active connected component.

    void supportFrontierOrder(int* order) const {
        std::vector<int> compId(Nq_, -1);
        int nComp = 0;
        for (int s = 0; s < Nq_; ++s) {
            if (compId[s] >= 0) continue;
            int cid = nComp++;
            std::deque<int> dq;
            dq.push_back(s);
            compId[s] = cid;
            while (!dq.empty()) {
                int u = dq.front(); dq.pop_front();
                for (int v : gq_.neighbors[u]) if (compId[v] < 0) { compId[v] = cid; dq.push_back(v); }
            }
        }

        std::vector<int> support(Nq_, 0);
        for (int i = 0; i < Nq_; ++i) support[i] = domainSupport(i);
        std::vector<char> picked(Nq_, 0);
        int activeComp = -1;

        // Tie-breaker hierarchy (v6.8.0 — ring-first for rigid core mapping):
        //   1. Minimum domain support (most constrained)
        //   2. Maximum mapped neighbors (frontier preference)
        //   3. Ring membership (ring atoms before chain atoms)
        //   4. Maximum query degree
        //   5. Minimum label (determinism)
        auto pickBest = [&](int comp, bool frontierOnly) {
            int best = -1;
            int bestSupport = 0x7fffffff, bestMapped = -1;
            int bestRing = -1, bestDeg = -1, bestLabel = 0x7fffffff;
            for (int q = 0; q < Nq_; ++q) {
                if (picked[q]) continue;
                if (comp >= 0 && compId[q] != comp) continue;
                int mappedNbrs = 0;
                for (int nb : gq_.neighbors[q]) mappedNbrs += picked[nb] ? 1 : 0;
                if (frontierOnly && mappedNbrs == 0) continue;
                int inRing = gq_.ring[q] ? 1 : 0;
                if (support[q] < bestSupport
                    || (support[q] == bestSupport && mappedNbrs > bestMapped)
                    || (support[q] == bestSupport && mappedNbrs == bestMapped && inRing > bestRing)
                    || (support[q] == bestSupport && mappedNbrs == bestMapped && inRing == bestRing && gq_.degree[q] > bestDeg)
                    || (support[q] == bestSupport && mappedNbrs == bestMapped && inRing == bestRing && gq_.degree[q] == bestDeg && gq_.label[q] < bestLabel)) {
                    best = q;
                    bestSupport = support[q];
                    bestMapped = mappedNbrs;
                    bestRing = inRing;
                    bestDeg = gq_.degree[q];
                    bestLabel = gq_.label[q];
                }
            }
            return best;
        };

        for (int pos = 0; pos < Nq_; ++pos) {
            int q = -1;
            if (activeComp >= 0) q = pickBest(activeComp, true);
            if (q < 0 && activeComp >= 0) q = pickBest(activeComp, false);
            if (q < 0) q = pickBest(-1, true);
            if (q < 0) q = pickBest(-1, false);
            order[pos] = q;
            picked[q] = 1;
            activeComp = compId[q];
            bool anyLeft = false;
            for (int i = 0; i < Nq_; ++i)
                if (!picked[i] && compId[i] == activeComp) { anyLeft = true; break; }
            if (!anyLeft) activeComp = -1;
        }
    }

    // ----- Feasibility check -----------------------------------------------

    bool feasible(int qi, int tj) {
        // Degree check
        if (gt_.degree[tj] < gq_.degree[qi]) { ++prunesDegree_; return false; }

        // Ring-only check (cheapest boolean test — run early for fast exit)
        if (useRingOnly_ && gq_.ring[qi] && !gt_.ring[tj]) {
            ++prunesAtom_; return false;
        }

        // NLF-1 check
        if (!nlfOk((*qNLF1_)[qi].data(), static_cast<int>((*qNLF1_)[qi].size()),
                   (*tNLF1_)[tj].data(), static_cast<int>((*tNLF1_)[tj].size()))) {
            ++prunesNLF_; return false;
        }

        // NLF-2 check
        if (useTwoHop_
            && !nlfOk((*qNLF2_)[qi].data(), static_cast<int>((*qNLF2_)[qi].size()),
                      (*tNLF2_)[tj].data(), static_cast<int>((*tNLF2_)[tj].size()))) {
            ++prunesNLF_; return false;
        }

        // NLF-3 check
        if (useThreeHop_
            && !nlfOk((*qNLF3_)[qi].data(), static_cast<int>((*qNLF3_)[qi].size()),
                      (*tNLF3_)[tj].data(), static_cast<int>((*tNLF3_)[tj].size()))) {
            ++prunesNLF_; return false;
        }

        // Mapped-neighbor adjacency (bit-parallel fast path)
        const auto& nbSorted = qNeighborsByDegDesc_[qi];
        if (useBitParallel_) {
            const uint64_t* tAdj = gt_.adjLong[tj].data();
            for (int qk : nbSorted) {
                int tk = q2t_[qk];
                if (tk != -1 && (tAdj[tk >> 6] & (uint64_t(1) << (tk & 63))) == 0) {
                    ++prunesBond_; return false;
                }
            }
            if (bitParallelSufficient_ && !useInduced_ && !useStereo_) return true;
        }

        // Full bond-compatibility check
        for (int qk : nbSorted) {
            int tk = q2t_[qk];
            if (tk != -1 && !bondsCompatible(gq_, qi, qk, gt_, tj, tk, C_)) {
                ++prunesBond_; return false;
            }
        }

        // --- Induced subgraph check (v6.8.0) ---
        // If induced mode is enabled, verify that non-adjacent query atoms
        // do NOT map to adjacent target atoms.  This enforces strict ring-
        // system matching (no extra intra-subgraph bonds allowed).
        if (useInduced_) {
            const uint64_t* tAdj = gt_.adjLong[tj].data();
            for (int k = 0; k < Nq_; ++k) {
                if (k == qi) continue;
                int tk = q2t_[k];
                if (tk == -1) continue;
                // If query has NO bond qi-k but target HAS bond tj-tk → reject
                if (!gq_.hasBond(qi, k) &&
                    (tAdj[tk >> 6] & (uint64_t(1) << (tk & 63)))) {
                    ++prunesBond_; return false;
                }
            }
        }

        // --- Topological stereo parity check (v6.8.0) ---
        // Verify that the winding order of mapped neighbors
        // around a chiral center matches.  A simple property comparison
        // (R==R) is insufficient because chirality depends on the specific
        // neighbor ordering in the mapping.  We compute the local parity
        // sign from the mapped neighbor permutation.
        if (useStereo_ && C_.useChirality
            && !tetraParityCompatible(gq_, qi, gt_, tj, q2t_)) {
            ++prunesAtom_; return false;
        }

        if (!hasLocalSupport(qi, tj)) { ++prunesBond_; return false; }

        return true;
    }

    // ----- Insertion sort candidates by degree proximity --------------------

    void sortByDegreeProximity(int* cands, int n, int qi) const {
        int qd = qdeg_[qi];
        for (int i = 1; i < n; ++i) {
            int key = cands[i];
            int keyDist = tdeg_[key] - qd;
            if (keyDist < 0) keyDist = -keyDist;
            int keyDeg = tdeg_[key];
            int j = i - 1;
            while (j >= 0) {
                int d = tdeg_[cands[j]] - qd;
                if (d < 0) d = -d;
                if (d < keyDist) break;
                // Prefer higher target degree on ties (more constrained)
                if (d == keyDist && tdeg_[cands[j]] >= keyDeg) break;
                cands[j+1] = cands[j]; --j;
            }
            cands[j+1] = key;
        }
    }

    // ----- VF3-style pivot: prefer candidate with most mapped-neighbor overlap

    void applyPivotHeuristic(int* out, int count, int qi) const {
        if (count <= 1) return;
        // Count mapped neighbors of qi — pivot only valuable with >= 2
        int mappedNbrs = 0;
        for (int nk = 0; nk < gq_.degree[qi]; ++nk)
            if (q2t_[gq_.neighbors[qi][nk]] != -1) mappedNbrs++;
        if (mappedNbrs < 2) return;
        // Find candidate with most adjacency overlap to mapped neighbors
        int bestPivot = 0, bestOverlap = 0;
        for (int c = 0; c < count; c++) {
            int tj = out[c];
            int overlap = 0;
            for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
                int qk = gq_.neighbors[qi][nk];
                int tk = q2t_[qk];
                if (tk != -1) {
                    if (gt_.adjLong[tj][tk >> 6] & (uint64_t(1) << (tk & 63)))
                        overlap++;
                }
            }
            if (overlap > bestOverlap) { bestOverlap = overlap; bestPivot = c; }
        }
        if (bestPivot != 0) std::swap(out[0], out[bestPivot]);
    }

    // ----- Candidate selection (virtual, overridden by VF2/VF2PP) ----------

    virtual int selectCandidates(int qi, int* out) = 0;
    virtual void onMatch(int /*pos*/, int /*tj*/) {}
    virtual void onUnmatch(int /*pos*/, int /*tj*/) {}

    // ----- Backtracking (exists) -------------------------------------------

    // 1-step forward check: after tentative qi→tj, verify each unmapped
    // neighbor of qi still has at least one compatible candidate remaining.
    // Catches dead-end assignments early, avoiding deep recursion.
    bool forwardCheck(int qi, int tj) const {
        for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
            int qk = gq_.neighbors[qi][nk];
            if (q2t_[qk] != -1) continue; // already mapped
            // Check domain(qk) minus used mask minus tj has any candidate
            // that is also adjacent to tj
            bool hasCandidate = false;
            for (int w = 0; w < tWords_; ++w) {
                if (domRow(qk)[w] & ~usedMask_[w]) {
                    hasCandidate = true; break;
                }
            }
            if (!hasCandidate) return false;
        }
        return true;
    }

    void backtrack(const int* order, int pos) {
        if (tb_.expired()) { timedOut_ = true; return; }
        ++nodesVisited_;
        if (found_) return;
        if (pos == Nq_) { found_ = true; return; }
        int qi = order[pos];
        int nCands = selectCandidates(qi, candRow(pos));
        for (int c = 0; c < nCands; ++c) {
            int tj = candRow(pos)[c];
            ++candidatesTried_;
            if (!feasible(qi, tj)) continue;
            q2t_[qi] = tj; t2q_[tj] = qi;
            usedMask_[tj >> 6] |= uint64_t(1) << (tj & 63);
            onMatch(pos, tj);
            // Forward check: prune if any unmapped neighbor of qi has no candidates left
            if (forwardCheck(qi, tj)) {
                backtrack(order, pos + 1);
            }
            if (found_ || timedOut_) return;
            q2t_[qi] = -1; t2q_[tj] = -1;
            usedMask_[tj >> 6] &= ~(uint64_t(1) << (tj & 63));
            onUnmatch(pos, tj);
            ++backtracks_;
        }
    }

    // ----- Backtracking (enumerate all) ------------------------------------

    void enumerateRec(const int* order, int pos,
                      std::vector<std::vector<std::pair<int,int>>>& out,
                      int maxSolutions) {
        if (tb_.expired()) { timedOut_ = true; return; }
        ++nodesVisited_;
        if (static_cast<int>(out.size()) >= maxSolutions) return;
        if (pos == Nq_) {
            std::vector<std::pair<int,int>> map;
            map.reserve(Nq_);
            for (int k = 0; k < Nq_; ++k)
                map.emplace_back(order[k], q2t_[order[k]]);
            out.push_back(std::move(map));
            return;
        }
        int qi = order[pos];
        int nCands = selectCandidates(qi, candRow(pos));
        for (int c = 0; c < nCands; ++c) {
            if (tb_.expired()) { timedOut_ = true; return; }
            int tj = candRow(pos)[c];
            ++candidatesTried_;
            if (!feasible(qi, tj)) continue;
            q2t_[qi] = tj; t2q_[tj] = qi;
            usedMask_[tj >> 6] |= uint64_t(1) << (tj & 63);
            onMatch(pos, tj);
            if (forwardCheck(qi, tj)) {
                enumerateRec(order, pos + 1, out, maxSolutions);
            }
            q2t_[qi] = -1; t2q_[tj] = -1;
            usedMask_[tj >> 6] &= ~(uint64_t(1) << (tj & 63));
            onUnmatch(pos, tj);
            ++backtracks_;
            if (static_cast<int>(out.size()) >= maxSolutions || timedOut_) return;
        }
    }

    // ----- Check if any query atom has empty domain ------------------------

    bool anyEmptyDomain() const {
        for (int i = 0; i < Nq_; ++i) {
            bool hasCandidate = false;
            for (int w = 0; w < tWords_; ++w) {
                if (domRow(i)[w] != 0) { hasCandidate = true; break; }
            }
            if (!hasCandidate) return true;
        }
        return false;
    }
};

// ==========================================================================
// VF2Matcher -- classic VF2 with terminal-set candidate selection
// ==========================================================================

class VF2Matcher final : public AbstractVFMatcher {
    // Scratch for terminal mask computation
    std::vector<uint64_t> termMask_;
    std::vector<uint64_t> termAvail_;

public:
    VF2Matcher(const MolGraph& gq, const MolGraph& gt,
               const ChemOptions& C, int64_t timeoutMs)
        : AbstractVFMatcher(gq, gt, C, timeoutMs),
          termMask_(tWords_, 0),
          termAvail_(tWords_, 0)
    {}

    bool exists() {
        if (Nq_ == 0) return true;
        if (anyEmptyDomain()) return false;
        std::vector<int> order(Nq_);
        if (Nq_ > 30) supportFrontierOrder(order.data());
        else fastisoOrder(order.data());
        if (greedyProbe(order.data())) { found_ = true; return true; }
        backtrack(order.data(), 0);
        return found_;
    }

    void enumerate(int maxSolutions,
                   std::vector<std::vector<std::pair<int,int>>>& out) {
        if (Nq_ == 0) { out.push_back({}); return; }
        if (anyEmptyDomain()) return;
        std::vector<int> order(Nq_);
        if (Nq_ > 30) supportFrontierOrder(order.data());
        else fastisoOrder(order.data());
        enumerateRec(order.data(), 0, out, maxSolutions);
    }

    int selectCandidates(int qi, int* out) override {
        for (int w = 0; w < tWords_; ++w)
            termAvail_[w] = domRow(qi)[w] & ~usedMask_[w];

        bool constrained = false;
        for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
            int qk = gq_.neighbors[qi][nk];
            int tk = q2t_[qk];
            if (tk == -1) continue;
            constrained = true;
            for (int w = 0; w < tWords_; ++w)
                termAvail_[w] &= gt_.adjLong[tk][w];
        }
        if (!constrained) {
            std::memset(termMask_.data(), 0, tWords_ * sizeof(uint64_t));
            bool hasTerminal = false;
            for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
                int qk = gq_.neighbors[qi][nk];
                int tk = q2t_[qk];
                if (tk != -1) {
                    for (int w = 0; w < tWords_; ++w) termMask_[w] |= gt_.adjLong[tk][w];
                    hasTerminal = true;
                }
            }
            if (hasTerminal) {
                for (int w = 0; w < tWords_; ++w) termAvail_[w] &= termMask_[w];
            }
        }
        int count = 0;
        for (int w = 0; w < tWords_; ++w) {
            uint64_t bits = termAvail_[w];
            while (bits) {
                int tj = (w << 6) | ctz64(bits);
                bits &= bits - 1;
                if (mappedBondCompatibleOnly(qi, tj)) out[count++] = tj;
            }
        }
        if (count > 3) sortByDegreeProximity(out, count, qi);
        return count;
    }
};

// ==========================================================================
// VF2PPMatcher -- VF2++ with frontier mask and save/restore
// ==========================================================================

class VF2PPMatcher final : public AbstractVFMatcher {
    std::vector<uint64_t> tFrontierMask_;
    std::vector<uint64_t> candAvail_;
    // Flat frontier stack: single allocation instead of Nq heap allocs (v6.8.0)
    std::vector<uint64_t> frontierStackFlat_;
    bool pivotEnabled_ = false; // VF3-style pivot: only during exists()

    uint64_t* frontierRow(int pos) {
        return frontierStackFlat_.data() + pos * tWords_;
    }

public:
    VF2PPMatcher(const MolGraph& gq, const MolGraph& gt,
                 const ChemOptions& C, int64_t timeoutMs)
        : AbstractVFMatcher(gq, gt, C, timeoutMs),
          tFrontierMask_(tWords_, 0),
          candAvail_(tWords_, 0),
          frontierStackFlat_(static_cast<size_t>(gq.n) * ((gt.n + 63) >> 6), 0)
    {}

    bool exists() {
        if (Nq_ == 0) return true;
        if (anyEmptyDomain()) return false;
        pivotEnabled_ = true;
        std::vector<int> order(Nq_);
        if (Nq_ > 30)
            supportFrontierOrder(order.data());
        else
            fastisoOrder(order.data());
        if (greedyProbe(order.data())) { found_ = true; pivotEnabled_ = false; return true; }
        std::fill(tFrontierMask_.begin(), tFrontierMask_.end(), uint64_t(0));
        backtrack(order.data(), 0);
        pivotEnabled_ = false;
        return found_;
    }

    void enumerate(int maxSolutions,
                   std::vector<std::vector<std::pair<int,int>>>& out) {
        if (Nq_ == 0) { out.push_back({}); return; }
        if (anyEmptyDomain()) return;
        std::vector<int> order(Nq_);
        if (Nq_ > 30)
            supportFrontierOrder(order.data());
        else
            fastisoOrder(order.data());
        std::fill(tFrontierMask_.begin(), tFrontierMask_.end(), uint64_t(0));
        enumerateRec(order.data(), 0, out, maxSolutions);
    }

    int selectCandidates(int qi, int* out) override {
        for (int w = 0; w < tWords_; ++w)
            candAvail_[w] = domRow(qi)[w] & ~usedMask_[w];

        bool constrained = false;
        for (int nk = 0; nk < gq_.degree[qi]; ++nk) {
            int qk = gq_.neighbors[qi][nk];
            int tk = q2t_[qk];
            if (tk == -1) continue;
            constrained = true;
            for (int w = 0; w < tWords_; ++w)
                candAvail_[w] &= gt_.adjLong[tk][w];
        }

        if (!constrained) {
            bool hasFrontier = false;
            for (int w = 0; w < tWords_; ++w) {
                if (tFrontierMask_[w] != 0) { hasFrontier = true; break; }
            }
            if (hasFrontier) {
                bool qiConnected = false;
                for (int k = 0; k < gq_.degree[qi]; ++k)
                    if (q2t_[gq_.neighbors[qi][k]] != -1) { qiConnected = true; break; }
                if (qiConnected) {
                    for (int w = 0; w < tWords_; ++w)
                        candAvail_[w] &= tFrontierMask_[w];
                }
            }
        }

        int count = 0;
        for (int w = 0; w < tWords_; ++w) {
            uint64_t bits = candAvail_[w];
            while (bits) {
                int tj = (w << 6) | ctz64(bits);
                bits &= bits - 1;
                if (mappedBondCompatibleOnly(qi, tj)) out[count++] = tj;
            }
        }
        if (count > 3) sortByDegreeProximity(out, count, qi);
        if (pivotEnabled_ && count > 4) applyPivotHeuristic(out, count, qi);
        return count;
    }

    void onMatch(int pos, int tj) override {
        // Save frontier state (flat layout)
        std::memcpy(frontierRow(pos), tFrontierMask_.data(),
                    tWords_ * sizeof(uint64_t));
        // Add unmatched neighbors of tj to frontier
        for (int k = 0; k < gt_.degree[tj]; ++k) {
            int u = gt_.neighbors[tj][k];
            if (t2q_[u] == -1)
                tFrontierMask_[u >> 6] |= uint64_t(1) << (u & 63);
        }
        // Remove tj itself from frontier
        tFrontierMask_[tj >> 6] &= ~(uint64_t(1) << (tj & 63));
    }

    void onUnmatch(int pos, int /*tj*/) override {
        // Restore frontier state (flat layout)
        std::memcpy(tFrontierMask_.data(), frontierRow(pos),
                    tWords_ * sizeof(uint64_t));
    }
};

} // namespace detail

// ============================================================================
// Public API implementation
// ============================================================================

inline bool isSubstructure(const MolGraph& query, const MolGraph& target,
                           const ChemOptions& opts, int64_t timeoutMs) {
    if (query.n == 0) return true;
    if (!detail::quickPrescreen(query, target, opts)) return false;
    if (&query == &target) return true;

    if (auto pathMap = detail::findPathSubstructure(query, target, opts); !pathMap.empty())
        return true;

    // Small-molecule fast path: bypass full VF2PP setup overhead.
    // Avoids ensureCanonical(), NLF construction, domain init, AC-3,
    // and all the heap allocations in the VF2PP constructor.
    if (query.n <= detail::SmallMolMatcher::MAX_QUERY_ATOMS
        && target.n <= detail::SmallMolMatcher::MAX_TARGET_ATOMS
        && static_cast<int64_t>(query.n) * target.n <= 1200) {
        if (detail::isExactMatch(query, target, opts)) return true;
        detail::SmallMolMatcher sm(query, target, opts);
        return sm.exists();
    }

    if (detail::isExactMatch(query, target, opts)) return true;
    if (!opts.useChirality && !opts.useBondStereo
        && detail::sameCanonicalGraph(query, target)) return true;

    // --- Disconnected query splitting (v6.8.0) ---
    if (detail::countComponents(query) > 1) {
        auto comps = detail::splitComponents(query);
        for (const auto& comp : comps) {
            if (!isSubstructure(comp, target, opts, timeoutMs))
                return false;
        }
        return true;
    }

    detail::TimeBudget tb(timeoutMs);
    bool timedOut = false;

    int64_t remaining = tb.remainingMs();
    if (opts.matcherEngine == ChemOptions::MatcherEngine::VF2) {
        detail::VF2Matcher m(query, target, opts, remaining);
        return m.exists();
    } else {
        detail::VF2PPMatcher m(query, target, opts, remaining);
        return m.exists();
    }
}

inline std::vector<std::pair<int,int>> findSubstructure(
    const MolGraph& query, const MolGraph& target,
    const ChemOptions& opts, int64_t timeoutMs) {

    if (query.n == 0) return {};
    if (!detail::quickPrescreen(query, target, opts)) return {};
    if (&query == &target) {
        std::vector<std::pair<int,int>> id(query.n);
        for (int i = 0; i < query.n; ++i) id[i] = {i, i};
        return id;
    }

    if (auto pathMap = detail::findPathSubstructure(query, target, opts); !pathMap.empty())
        return pathMap;

    // Small-molecule fast path
    if (query.n <= detail::SmallMolMatcher::MAX_QUERY_ATOMS
        && target.n <= detail::SmallMolMatcher::MAX_TARGET_ATOMS
        && static_cast<int64_t>(query.n) * target.n <= 1200) {
        if (detail::isExactMatch(query, target, opts)) {
            std::vector<std::pair<int,int>> id(query.n);
            for (int i = 0; i < query.n; ++i) id[i] = {i, i};
            return id;
        }
        detail::SmallMolMatcher sm(query, target, opts);
        return sm.findOne();
    }

    if (detail::isExactMatch(query, target, opts)) {
        std::vector<std::pair<int,int>> id(query.n);
        for (int i = 0; i < query.n; ++i) id[i] = {i, i};
        return id;
    }
    if (!opts.useChirality && !opts.useBondStereo
        && detail::sameCanonicalGraph(query, target))
        return detail::canonicalIsoMap(query, target);

    detail::TimeBudget tb(timeoutMs);
    bool timedOut = false;

    std::vector<std::vector<std::pair<int,int>>> all;
    int64_t remaining = tb.remainingMs();
    if (opts.matcherEngine == ChemOptions::MatcherEngine::VF2) {
        detail::VF2Matcher m(query, target, opts, remaining);
        m.enumerate(1, all);
    } else {
        detail::VF2PPMatcher m(query, target, opts, remaining);
        m.enumerate(1, all);
    }
    return all.empty() ? std::vector<std::pair<int,int>>{} : std::move(all[0]);
}

inline std::vector<std::vector<std::pair<int,int>>> findAllSubstructures(
    const MolGraph& query, const MolGraph& target,
    const ChemOptions& opts, int64_t timeoutMs) {

    if (query.n == 0) return {{}};
    if (!detail::quickPrescreen(query, target, opts)) return {};
    if (&query == &target || detail::isExactMatch(query, target, opts)) {
        std::vector<std::pair<int,int>> id(query.n);
        for (int i = 0; i < query.n; ++i) id[i] = {i, i};
        return {id};
    }
    if (!opts.useChirality && !opts.useBondStereo
        && detail::sameCanonicalGraph(query, target))
        return {detail::canonicalIsoMap(query, target)};

    detail::TimeBudget tb(timeoutMs);
    bool timedOut = false;

    constexpr int kMaxSolutions = 10000;
    std::vector<std::vector<std::pair<int,int>>> out;
    int64_t remaining = tb.remainingMs();
    if (opts.matcherEngine == ChemOptions::MatcherEngine::VF2) {
        detail::VF2Matcher m(query, target, opts, remaining);
        m.enumerate(kMaxSolutions, out);
    } else {
        detail::VF2PPMatcher m(query, target, opts, remaining);
        m.enumerate(kMaxSolutions, out);
    }
    return out;
}

} // namespace smsd

#endif // SMSD_VF2PP_HPP

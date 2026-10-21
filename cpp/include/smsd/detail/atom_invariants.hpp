/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph per-atom invariants used to prune searches: neighbour-label frequencies,
 * pharmacophore features, hydrogen counts, and prewarm() for shared graphs.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// Pattern fingerprint (deprecated; nothing reads it)
// ============================================================================

inline void MolGraph::ensurePatternFP() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (patternFPComputed_) return;
    for (int w = 0; w < FP_WORDS; ++w) patternFP_[w] = 0;
    // Encode atom features: element + aromaticity → 1 bit per type
    for (int i = 0; i < n; ++i) {
        uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 31 + (aromatic[i] ? 1 : 0);
        patternFP_[0] |= uint64_t(1) << (h & 63);
    }
    // Encode bond-pair features: (Z1, bondOrder, Z2) for each edge
    for (int i = 0; i < n; ++i) {
        for (int j : neighbors[i]) {
            if (j <= i) continue;  // each edge once
            int bo = bondOrder(i, j);
            uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 997
                       + static_cast<uint32_t>(atomicNum[j]) * 31
                       + static_cast<uint32_t>(bo) * 7;
            patternFP_[1] |= uint64_t(1) << (h & 63);
            // Also encode reverse pair for asymmetric element pairs
            uint32_t h2 = static_cast<uint32_t>(atomicNum[j]) * 997
                        + static_cast<uint32_t>(atomicNum[i]) * 31
                        + static_cast<uint32_t>(bo) * 7;
            patternFP_[1] |= uint64_t(1) << (h2 & 63);
        }
    }
    // Encode 2-hop path features: (Z1, Z2, Z3) for paths of length 2
    for (int i = 0; i < n; ++i) {
        for (int j : neighbors[i]) {
            for (int k : neighbors[j]) {
                if (k == i) continue;
                uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 997
                           + static_cast<uint32_t>(atomicNum[j]) * 31
                           + static_cast<uint32_t>(atomicNum[k]);
                patternFP_[2] |= uint64_t(1) << (h & 63);
            }
        }
    }
    // Encode ring membership + degree features
    for (int i = 0; i < n; ++i) {
        uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 31
                   + static_cast<uint32_t>(degree[i]) * 7
                   + (ring[i] ? 3u : 0u);
        patternFP_[3] |= uint64_t(1) << (h & 63);
    }
    patternFPComputed_ = true;
}

// ============================================================================
// NLF (Neighbour Label Frequency) helpers
// ============================================================================

inline std::vector<int> MolGraph::freqMapToSortedArray(
        const std::unordered_map<int,int>& freq) {
    if (freq.empty()) return {};
    std::vector<std::pair<int,int>> pairs(freq.begin(), freq.end());
    std::sort(pairs.begin(), pairs.end());
    std::vector<int> arr(pairs.size() * 2);
    for (size_t i = 0; i < pairs.size(); i++) {
        arr[i * 2]     = pairs[i].first;
        arr[i * 2 + 1] = pairs[i].second;
    }
    return arr;
}

inline int MolGraph::nlfLabel(const MolGraph& g, int idx) {
    return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
}

inline std::vector<int> MolGraph::buildNLF1(const MolGraph& g, int idx) {
    std::unordered_map<int,int> freq;
    for (int nb : g.neighbors[idx]) freq[nlfLabel(g, nb)]++;
    return freqMapToSortedArray(freq);
}

inline std::vector<int> MolGraph::buildRadiusNLF(const MolGraph& g, int idx, int radius) {
    std::unordered_map<int,int> freq;
    std::vector<uint8_t> seen(g.n, 0);
    std::vector<int> queue{idx};
    seen[idx] = 1;
    size_t begin = 0;
    for (int depth = 0; depth < radius && begin < queue.size(); ++depth) {
        const size_t end = queue.size();
        for (; begin < end; ++begin) {
            for (int nb : g.neighbors[queue[begin]]) {
                if (seen[nb]) continue;
                seen[nb] = 1;
                queue.push_back(nb);
                ++freq[nlfLabel(g, nb)];
            }
        }
    }
    return freqMapToSortedArray(freq);
}

inline std::vector<int> MolGraph::buildNLF2(const MolGraph& g, int idx) {
    return buildRadiusNLF(g, idx, 2);
}

inline std::vector<int> MolGraph::buildNLF3(const MolGraph& g, int idx) {
    return buildRadiusNLF(g, idx, 3);
}

inline std::vector<std::vector<int>> MolGraph::buildAllNLF(
        const MolGraph& g, NLFBuilder builder) {
    std::vector<std::vector<int>> nlf(g.n);
    for (int i = 0; i < g.n; i++) nlf[i] = builder(g, i);
    return nlf;
}

inline const std::vector<std::vector<int>>& MolGraph::getNLF1() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (cachedNLF1.empty() && n > 0)
        cachedNLF1 = buildAllNLF(*this, buildNLF1);
    return cachedNLF1;
}

inline const std::vector<std::vector<int>>& MolGraph::getNLF2() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (cachedNLF2.empty() && n > 0)
        cachedNLF2 = buildAllNLF(*this, buildNLF2);
    return cachedNLF2;
}

inline const std::vector<std::vector<int>>& MolGraph::getNLF3() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (cachedNLF3.empty() && n > 0)
        cachedNLF3 = buildAllNLF(*this, buildNLF3);
    return cachedNLF3;
}

inline bool MolGraph::nlfOk(const std::vector<int>& fq, const std::vector<int>& ft) {
    size_t fi = 0, ti = 0;
    while (fi < fq.size()) {
        int fLabel = fq[fi], fFreq = fq[fi + 1];
        while (ti < ft.size() && ft[ti] < fLabel) ti += 2;
        if (ti >= ft.size() || ft[ti] != fLabel || ft[ti + 1] < fFreq)
            return false;
        fi += 2;
    }
    return true;
}

inline bool MolGraph::nlfCheckOk(
        int qi, int tj,
        const std::vector<std::vector<int>>& qNLF1,
        const std::vector<std::vector<int>>& tNLF1,
        const std::vector<std::vector<int>>& qNLF2,
        const std::vector<std::vector<int>>& tNLF2,
        const std::vector<std::vector<int>>& qNLF3,
        const std::vector<std::vector<int>>& tNLF3,
        bool useTwoHop, bool useThreeHop) {
    if (!nlfOk(qNLF1[qi], tNLF1[tj])) return false;
    if (useTwoHop && !nlfOk(qNLF2[qi], tNLF2[tj])) return false;
    return !useThreeHop || nlfOk(qNLF3[qi], tNLF3[tj]);
}

// ============================================================================
// Cached atom properties and prewarm()
// ============================================================================

inline const std::vector<std::vector<int>>& MolGraph::getNeighborsByDegDesc() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (cachedNeighborsByDegDesc.empty() && n > 0) {
        cachedNeighborsByDegDesc.resize(n);
        for (int i = 0; i < n; ++i) {
            auto nb = neighbors[i]; // copy
            std::sort(nb.begin(), nb.end(), [this](int a, int b) {
                return degree[a] > degree[b]; // descending
            });
            cachedNeighborsByDegDesc[i] = std::move(nb);
        }
    }
    return cachedNeighborsByDegDesc;
}

namespace detail {

/** Explicit valence of an atom for the valence model: its bond orders, with an aromatic bond counted once. */
inline int explicitValence(const MolGraph& graph, int atom) {
    int64_t bondOrderSum = 0;
    for (int neighbor : graph.neighbors[atom]) {
        const int order = graph.bondOrder(atom, neighbor);
        bondOrderSum += order == 4 ? 1 : order;
    }
    return static_cast<int>(std::min(bondOrderSum, int64_t(std::numeric_limits<int>::max()) - 1));
}

} // namespace detail

inline int MolGraph::implicitHydrogenCount(int atom) const {
    if (atom < 0 || atom >= n) throw std::out_of_range("Atom index out of range");
    if (hydrogenCountsKnown && atom < static_cast<int>(hydrogenCount.size())) return hydrogenCount[atom];
    const int charge = atom < static_cast<int>(formalCharge.size()) ? formalCharge[atom] : 0;
    return detail::computeImplicitH(atomicNum[atom], aromatic[atom] != 0, detail::explicitValence(*this, atom), charge);
}

inline int64_t MolGraph::attachedHydrogenCount(int atom) const {
    if (atom < 0 || atom >= n) throw std::out_of_range("Atom index out of range");
    int64_t explicitHydrogens = 0;
    for (int neighbor : neighbors[atom]) explicitHydrogens += atomicNum[neighbor] == 1;
    return implicitHydrogenCount(atom) + explicitHydrogens;
}

inline const std::vector<int>& MolGraph::getPharmacophoreFeatures() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (!pharmComputed_) {
        cachedPharmacophoreFeatures.resize(n);
        for (int i = 0; i < n; ++i) {
            int z = atomicNum[i];
            int charge = formalCharge[i];
            bool arom = (aromatic[i] != 0);
            int deg = degree[i];
            int features = 0;
            const int64_t hydrogens = attachedHydrogenCount(i);

            // H-bond donor: N-H, O-H, S-H
            if (z == 7 || z == 8 || z == 16) {
                if (hydrogens > 0) features |= 1;
            }

            // H-bond acceptor: N (not pyrrole), O, F, S (not thiophene)
            // Pyrrole N: aromatic + has H (lone pair donated to pi) — NOT acceptor
            // Pyridine N: aromatic + no H (lone pair in plane) — IS acceptor
            bool isPyrroleN = false;
            if (z == 7 && arom) {
                isPyrroleN = (hydrogens > 0);
            }
            if ((z == 7 && charge <= 0 && !isPyrroleN) || z == 8 || z == 9 ||
                (z == 16 && !arom))
                features |= 2;

            // Positive ionisable: basic amines (not amide, not aniline)
            if (z == 7 && charge > 0) features |= 4;
            if (z == 7 && !arom && deg <= 3 && charge == 0) {
                bool isAmide = false, isAniline = false, hasMultipleBond = false;
                for (int nb : neighbors[i]) {
                    if (bondOrder(i, nb) != 1) hasMultipleBond = true;
                    if (atomicNum[nb] == 6) {
                        if (aromatic[nb]) { isAniline = true; break; }
                        for (int nb2 : neighbors[nb]) {
                            if (nb2 != i && atomicNum[nb2] == 8 && bondOrder(nb, nb2) == 2) {
                                isAmide = true; break;
                            }
                        }
                    }
                    if (isAmide) break;
                }
                if (!isAmide && !isAniline && !hasMultipleBond) features |= 4;
            }

            // Negative ionisable: carboxylic, phosphoric, sulfonic
            if (z == 8 && charge < 0) features |= 8;
            if (z == 8 && charge == 0 && deg == 1) {
                for (int nb : neighbors[i]) {
                    int nbZ = atomicNum[nb];
                    if (nbZ == 6 && bondOrder(i, nb) == 1) {
                        for (int nb2 : neighbors[nb]) {
                            if (nb2 != i && atomicNum[nb2] == 8 && bondOrder(nb, nb2) == 2)
                                features |= 8;
                        }
                    }
                    if (nbZ == 15 || nbZ == 16) features |= 8;
                }
            }

            // Aromatic
            if (arom) features |= 16;

            // Hydrophobic: non-aromatic C (no hetero neighbours), or halogens
            if (z == 6 && !arom) {
                bool hasHetero = false;
                for (int nb : neighbors[i]) {
                    int nbZ = atomicNum[nb];
                    if (nbZ != 6 && nbZ != 1) { hasHetero = true; break; }
                }
                if (!hasHetero) features |= 32;
            }
            if (z == 17 || z == 35 || z == 53) features |= 32;

            cachedPharmacophoreFeatures[i] = features;
        }
        pharmComputed_ = true;
    }
    return cachedPharmacophoreFeatures;
}

inline void MolGraph::prewarm() const {
    ensureCanonical();
    ensureRingCounts();
    ensureRingSystems();
    getPharmacophoreFeatures();
}

} // namespace smsd

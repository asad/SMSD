/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * RASCAL — O(V+E) similarity upper-bound screening.
 * Raymond, Gardiner & Willett, JCAMD 2002.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace smsd {

namespace screening_detail {
inline int ubLabel(const MolGraph& g, int i, const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    int z = g.atomicNum[i];
    if (C.tautomerAware && (z == 6 || z == 7 || z == 8 || z == 16 || z == 34)) z = 6;
    int key = z << 2;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT) key |= g.aromatic[i] ? 2 : 0;
    if (C.ringMatchesRingOnly) key |= g.ring[i] ? 1 : 0;
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

} // namespace screening_detail

inline double similarityUpperBound(const MolGraph& g1, const MolGraph& g2, const ChemOptions& chemistry) {
    if (g1.n == 0 || g2.n == 0) return 0.0;
    int atomOverlap = screening_detail::labelFrequencyUpperBound(g1, g2, chemistry);
    if (atomOverlap == 0) return 0.0;
    int total = g1.n + g2.n;
    return static_cast<double>(atomOverlap) / (total - atomOverlap);
}

inline double similarityUpperBound(const MolGraph& g1, const MolGraph& g2) {
    return similarityUpperBound(g1, g2, ChemOptions());
}



/// Batch RASCAL screening: find all targets with similarity >= threshold.
/// Returns indices of passing targets.
inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold, const ChemOptions& chemistry) {
    std::vector<int> hits;
    hits.reserve(targets.size() / 10); // expect ~10% pass
    for (int i = 0; i < static_cast<int>(targets.size()); ++i) {
        if (similarityUpperBound(query, targets[i], chemistry) >= threshold) {
            hits.push_back(i);
        }
    }
    return hits;
}

inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold) {
    return screenTargets(query, targets, threshold, ChemOptions());
}

} // namespace smsd

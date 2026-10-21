/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Label-histogram similarity upper bound for screening, O(V+E).
 * After Raymond, Gardiner & Willett, J. Comput.-Aided Mol. Des. 2002.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace smsd {

namespace screening_detail {
/**
 * The label two atoms must share to pair in the bound: the element, the aromatic flag under
 * strict aromaticity and the ring flag when ring atoms match ring atoms only. Tautomerism
 * changes proton positions and bond orders, not elements, so atomCompatible keeps elements apart
 * under tautomer-aware matching and the label does too. The Java MCS.similarityUpperBound uses
 * the same label.
 */
inline int ubLabel(const MolGraph& g, int i, const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    const int z = g.atomicNum[i];
    // Shifted as unsigned, so that an out-of-range negative atomic number keeps its label without undefined behaviour.
    int key = static_cast<int>(static_cast<unsigned>(z) << 2);
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

/**
 * Query label counts held in a dense table, so one query can be bounded against many targets
 * without building hash maps per target. Gives the same value as labelFrequencyUpperBound.
 * The query graph must outlive this object.
 */
class QueryLabelCounts {
public:
    /** Labels at or above this limit, or below zero, use the hashed path. */
    static constexpr int denseLimit = 1 << 10;

    QueryLabelCounts(const MolGraph& query, const ChemOptions& chemistry)
        : query_(query), chemistry_(chemistry), counts_(denseLimit, 0) {
        for (int atom = 0; atom < query.n && dense_; ++atom) {
            const int label = ubLabel(query, atom, chemistry);
            if (label < 0 || label >= denseLimit) dense_ = false;
            else ++counts_[label];
        }
    }

    /** scratch must hold denseLimit zeros; it is left that way for the next target. */
    int upperBound(const MolGraph& target, std::vector<int>& scratch) const {
        if (query_.n == 0 || target.n == 0) return 0;
        if (!dense_) return labelFrequencyUpperBound(query_, target, chemistry_);
        int overlap = 0, atom = 0;
        for (; atom < target.n; ++atom) {
            const int label = ubLabel(target, atom, chemistry_);
            if (label < 0 || label >= denseLimit) break;
            // The first min(query, target) atoms of each label pair up.
            if (++scratch[label] <= counts_[label]) ++overlap;
        }
        const bool complete = atom == target.n;
        for (int counted = 0; counted < atom; ++counted) scratch[ubLabel(target, counted, chemistry_)] = 0;
        if (!complete) return labelFrequencyUpperBound(query_, target, chemistry_);
        return std::min(overlap, std::min(query_.n, target.n));
    }

    /** Same value as similarityUpperBound(query, target, chemistry). */
    double similarity(const MolGraph& target, std::vector<int>& scratch) const {
        if (query_.n == 0 || target.n == 0) return 0.0;
        const int atomOverlap = upperBound(target, scratch);
        if (atomOverlap == 0) return 0.0;
        const int total = query_.n + target.n;
        return static_cast<double>(atomOverlap) / (total - atomOverlap);
    }

private:
    const MolGraph& query_;
    ChemOptions chemistry_;
    std::vector<int> counts_;
    bool dense_ = true;
};

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



/// Batch similarity-bound screening: find all targets with similarity >= threshold.
/// Returns indices of passing targets.
inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold, const ChemOptions& chemistry) {
    std::vector<int> hits;
    hits.reserve(targets.size() / 10); // expect ~10% pass
    const screening_detail::QueryLabelCounts queryCounts(query, chemistry);
    std::vector<int> scratch(screening_detail::QueryLabelCounts::denseLimit, 0);
    for (int i = 0; i < static_cast<int>(targets.size()); ++i) {
        if (queryCounts.similarity(targets[i], scratch) >= threshold) {
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

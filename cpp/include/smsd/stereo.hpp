/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#pragma once
#ifndef SMSD_STEREO_HPP
#define SMSD_STEREO_HPP

#include "mol_graph.hpp"
#include <algorithm>
#include <array>

namespace smsd {
namespace stereo_detail {

/** Ligand order used by graph tetrahedral tags, with -1 for the implicit ligand. */
inline bool tetrahedralLigands(const std::vector<int>& neighbors, int atom,
                               std::array<int, 4>& ligands) {
    if (neighbors.size() != 3 && neighbors.size() != 4) return false;
    const int hydrogenPosition = !neighbors.empty() && neighbors[0] < atom ? 1 : 0;
    int position = 0;
    for (int slot = 0; slot < 4; ++slot)
        ligands[slot] = neighbors.size() == 3 && slot == hydrogenPosition ? -1 : neighbors[position++];
    return true;
}

/** Compare winding once three explicit ligands determine the mapped permutation. */
template <typename Mapping>
inline bool tetrahedralParityMatches(int queryTag, const std::array<int, 4>& queryLigands,
                                     int targetTag, const std::array<int, 4>& targetLigands,
                                     const Mapping& mapping) {
    std::array<int, 4> permutation{{-1, -1, -1, -1}};
    std::array<bool, 4> used{};
    int mappedCount = 0;
    for (int slot = 0; slot < 4; ++slot) {
        const int queryLigand = queryLigands[slot];
        if (queryLigand < 0 || mapping[queryLigand] < 0) continue;
        const auto found = std::find(targetLigands.begin(), targetLigands.end(), mapping[queryLigand]);
        if (found == targetLigands.end()) return false;
        const int targetSlot = static_cast<int>(found - targetLigands.begin());
        if (used[targetSlot]) return false;
        used[targetSlot] = true;
        permutation[slot] = targetSlot;
        ++mappedCount;
    }
    if (mappedCount < 3) return true;
    for (int slot = 0; slot < 4; ++slot) {
        if (permutation[slot] >= 0) continue;
        for (int targetSlot = 0; targetSlot < 4; ++targetSlot) {
            if (!used[targetSlot]) {
                permutation[slot] = targetSlot;
                used[targetSlot] = true;
                break;
            }
        }
    }
    int inversions = 0;
    for (int first = 0; first < 4; ++first)
        for (int second = first + 1; second < 4; ++second)
            inversions += permutation[first] > permutation[second];
    return (queryTag == targetTag) == ((inversions & 1) == 0);
}

inline int doubleBondReference(const MolGraph& graph, int atom, int otherEnd) {
    for (int neighbor : graph.neighbors[atom]) if (neighbor != otherEnd) return neighbor;
    return -1;
}

/** Compare E/Z relative to mapped ligands, rather than either graph's own references. */
template <typename Mapping>
inline bool doubleBondParityMatches(const MolGraph& query, int first, int second,
                                     const MolGraph& target, const Mapping& mapping) {
    const int mappedFirst = mapping[first], mappedSecond = mapping[second];
    if (mappedFirst < 0 || mappedSecond < 0) return true;
    int queryConfiguration = query.dbStereo(first, second);
    const int targetConfiguration = target.dbStereo(mappedFirst, mappedSecond);
    if (queryConfiguration == 0 || targetConfiguration == 0) return true;
    for (int end = 0; end < 2; ++end) {
        const int queryEnd = end == 0 ? first : second;
        const int otherQueryEnd = end == 0 ? second : first;
        const int targetEnd = mapping[queryEnd], otherTargetEnd = mapping[otherQueryEnd];
        const int queryReference = doubleBondReference(query, queryEnd, otherQueryEnd);
        const int targetReference = doubleBondReference(target, targetEnd, otherTargetEnd);
        int mappedLigand = -1;
        for (int neighbor : query.neighbors[queryEnd]) {
            if (neighbor != otherQueryEnd && mapping[neighbor] >= 0) { mappedLigand = neighbor; break; }
        }
        if (mappedLigand < 0 || queryReference < 0 || targetReference < 0) return true;
        if (!target.hasBond(targetEnd, mapping[mappedLigand])) return false;
        if ((mappedLigand != queryReference) != (mapping[mappedLigand] != targetReference))
            queryConfiguration = 3 - queryConfiguration;
    }
    return queryConfiguration == targetConfiguration;
}

} // namespace stereo_detail
} // namespace smsd
#endif

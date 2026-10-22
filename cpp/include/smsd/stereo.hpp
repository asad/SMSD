/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once
#ifndef SMSD_STEREO_HPP
#define SMSD_STEREO_HPP

#include "smsd/mol_graph.hpp"
#include <algorithm>
#include <array>

namespace smsd {
namespace stereo_detail {

/**
 * Whether an atom of this element with three neighbours and no hydrogen is a tetrahedral centre
 * whose lone pair is the fourth ligand: N, P, S, As, Se, Sb and Te.
 */
inline bool lonePairCentre(int atomicNumber) {
    return atomicNumber == 7 || atomicNumber == 15 || atomicNumber == 16 || atomicNumber == 33
        || atomicNumber == 34 || atomicNumber == 51 || atomicNumber == 52;
}

/**
 * Whether an atom with `degree` neighbours and `hydrogens` implicit hydrogens is a tetrahedral
 * centre: four neighbours and no hydrogen, or three neighbours and either one hydrogen or, for
 * an element of lonePairCentre(), none.
 */
inline bool tetrahedralCentre(size_t degree, int hydrogens, int atomicNumber) {
    return degree == 4 ? hydrogens == 0
         : degree == 3 && (hydrogens == 1 || (hydrogens == 0 && lonePairCentre(atomicNumber)));
}

/**
 * Ligand order of a tetrahedral tag over `neighbors`, the neighbours of `atom` in order, with -1
 * for the implicit ligand of three neighbours (a hydrogen or a lone pair). The implicit ligand
 * comes second when the first neighbour precedes the atom, as after the atom it was written from
 * in SMILES, and first otherwise. False unless there are three or four neighbours; a SMARTS atom,
 * which has no hydrogen count of its own, takes this order as it is.
 */
inline bool tetrahedralLigands(const std::vector<int>& neighbors, int atom,
                               std::array<int, 4>& ligands) {
    if (neighbors.size() != 3 && neighbors.size() != 4) return false;
    const int hydrogenPosition = !neighbors.empty() && neighbors[0] < atom ? 1 : 0;
    int position = 0;
    for (int slot = 0; slot < 4; ++slot)
        ligands[slot] = neighbors.size() == 3 && slot == hydrogenPosition ? -1 : neighbors[position++];
    return true;
}

/**
 * Ligand order of the tetrahedral tag of `atom` in `graph`, as above, when the atom is a
 * tetrahedralCentre() by its implicit hydrogens. Java's MolGraph.tetrahedralLigands() applies the
 * same rule.
 */
inline bool tetrahedralLigands(const MolGraph& graph, int atom, std::array<int, 4>& ligands) {
    const size_t degree = graph.neighbors[atom].size();
    if (degree != 3 && degree != 4) return false;
    return tetrahedralCentre(degree, graph.implicitHydrogenCount(atom), graph.atomicNum[atom])
        && tetrahedralLigands(graph.neighbors[atom], atom, ligands);
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

/** The ligand test of doubleBondParityMatches() that takes every mapped neighbour of an end. */
struct AnyLigand {
    bool operator()(int, int) const { return true; }
};

/**
 * Compare E/Z relative to mapped ligands, rather than either graph's own references. A query
 * bond without a configuration matches any target bond. With `embedding`, as in a substructure
 * search, a configured query bond needs a configured target bond, as a SMARTS pair of
 * directional bonds does; otherwise a target bond without a configuration matches as well. A
 * mapped neighbour of an end is its ligand when `joined(end, neighbour)` holds, as when the bond
 * between them is selected; Java's Substructure applies the same rule.
 */
template <typename Mapping, typename Joined = AnyLigand>
inline bool doubleBondParityMatches(const MolGraph& query, int first, int second,
                                     const MolGraph& target, const Mapping& mapping, bool embedding = false,
                                     const Joined& joined = Joined{}) {
    const int mappedFirst = mapping[first], mappedSecond = mapping[second];
    if (mappedFirst < 0 || mappedSecond < 0) return true;
    int queryConfiguration = query.dbStereo(first, second);
    if (queryConfiguration == 0) return true;
    const int targetConfiguration = target.dbStereo(mappedFirst, mappedSecond);
    if (targetConfiguration == 0) return !embedding;
    for (int end = 0; end < 2; ++end) {
        const int queryEnd = end == 0 ? first : second;
        const int otherQueryEnd = end == 0 ? second : first;
        const int targetEnd = mapping[queryEnd], otherTargetEnd = mapping[otherQueryEnd];
        const int queryReference = doubleBondReference(query, queryEnd, otherQueryEnd);
        const int targetReference = doubleBondReference(target, targetEnd, otherTargetEnd);
        int mappedLigand = -1;
        for (int neighbor : query.neighbors[queryEnd]) {
            if (neighbor != otherQueryEnd && mapping[neighbor] >= 0 && joined(queryEnd, neighbor)) {
                mappedLigand = neighbor;
                break;
            }
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

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SmartsQuery, the parsed SMARTS pattern: adjacency, atom and bond tests, stereo
 * and component checks of a mapping, and recursive SMARTS evaluation.
 *
 * Part of smarts_parser.hpp.
 */
#pragma once

#include "smsd/stereo.hpp"
#include "smsd/detail/smarts_evaluation.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// SmartsQuery -- the parsed SMARTS pattern
// ============================================================================


struct SmartsQuery {
    std::vector<smarts_detail::SmartsAtom> atoms;
    std::vector<smarts_detail::SmartsBond> bonds;

    // Adjacency for the query graph
    std::vector<std::vector<int>> neighbors;   // neighbour atom indices
    std::vector<std::vector<int>> bondIndices; // parallel to neighbours: bond index

    /**
     * Expression nodes of each atom and bond as parseSMARTS() counted them, which weigh how
     * often a search reads the clock; a pattern assembled by hand may leave them empty, and
     * then each expression counts one node.
     */
    std::vector<uint32_t> atomNodes, bondNodes;

    int atomCount() const { return static_cast<int>(atoms.size()); }
    size_t atomNodesOf(int atom) const {
        return static_cast<size_t>(atom) < atomNodes.size() ? atomNodes[static_cast<size_t>(atom)] : 1;
    }
    size_t bondNodesOf(int bond) const {
        return static_cast<size_t>(bond) < bondNodes.size() ? bondNodes[static_cast<size_t>(bond)] : 1;
    }

    // Build adjacency from atoms/bonds
    void buildAdjacency();

    // Primitive compatibility helpers
    bool atomCompat(int qi, const MolGraph& target, int ti) const;
    bool bondCompat(int bi, const MolGraph& target, int ti, int tj) const;

    /**
     * A plain breadth-first order from atom 0, for debugging and simple callers.
     *
     * The optimised matcher below uses a target-aware, anchor-preserving order
     * instead.
     */
    std::vector<int> matchingOrder() const;


    bool hasStereochemistry() const;
    bool stereoCompatible(const MolGraph& target, const std::vector<int>& mapping) const;

    /**
     * The component group of each atom, -1 outside one, or empty when the pattern has none, see
     * SmartsParser::componentGroups().
     */
    std::vector<int> componentGroups;
    /** Whether a complete mapping is checked beyond its atoms and bonds: stereochemistry or component groups. */
    bool checksMapping() const { return !componentGroups.empty() || hasStereochemistry(); }
    /** Whether the mapped atoms of each component group lie in one target component, and those of different groups in different ones. */
    bool componentsCompatible(const MolGraph& target, const std::vector<int>& mapping) const;
    /** stereoCompatible() and componentsCompatible() together. */
    bool mappingCompatible(const MolGraph& target, const std::vector<int>& mapping) const {
        return componentsCompatible(target, mapping) && stereoCompatible(target, mapping);
    }

};

inline void SmartsQuery::buildAdjacency() {
    int n = atomCount();
    neighbors.assign(n, {});
    bondIndices.assign(n, {});
    for (int bi = 0; bi < static_cast<int>(bonds.size()); ++bi) {
        int a = bonds[bi].from, b = bonds[bi].to;
        neighbors[a].push_back(b);
        neighbors[b].push_back(a);
        bondIndices[a].push_back(bi);
        bondIndices[b].push_back(bi);
    }
}

inline bool SmartsQuery::hasStereochemistry() const {
    for (const auto& atom : atoms)
        if (smarts_detail::hasChiralityPrimitive(atom.expr)) return true;
    for (const auto& bond : bonds)
        if (smarts_detail::hasDirectionalPrimitive(bond.expr)) return true;
    return false;
}

inline bool SmartsQuery::stereoCompatible(const MolGraph& target, const std::vector<int>& mapping) const {
    using smarts_detail::StereoTruth;
    for (int atom = 0; atom < atomCount(); ++atom) {
        if (!smarts_detail::hasChiralityPrimitive(atoms[atom].expr)) continue;
        const int targetAtom = mapping[atom];
        std::array<int, 4> queryLigands{}, targetLigands{};
        const bool validLigands = stereo_detail::tetrahedralLigands(neighbors[atom], atom, queryLigands)
            && stereo_detail::tetrahedralLigands(target, targetAtom, targetLigands);
        const int targetTag = target.tetraChirality.empty() ? 0 : target.tetraChirality[targetAtom];
        const auto matched = smarts_detail::evalAtomExprWithStereo(atoms[atom].expr, target, targetAtom, [&](int queryTag) {
            if (queryTag == 0) return StereoTruth::TRUE_VALUE;
            if (targetTag == 0) return StereoTruth::FALSE_VALUE;
            const bool compatible = validLigands
                ? stereo_detail::tetrahedralParityMatches(queryTag, queryLigands, targetTag, targetLigands, mapping)
                : queryTag == targetTag;
            return compatible ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
        });
        if (matched != StereoTruth::TRUE_VALUE) return false;
    }

    std::vector<int> variableForBond(bonds.size(), -1), allowedDirections;
    for (int bondIndex = 0; bondIndex < static_cast<int>(bonds.size()); ++bondIndex) {
        const auto& bond = bonds[bondIndex];
        if (!smarts_detail::hasDirectionalPrimitive(bond.expr)) continue;
        const int first = mapping[bond.from], second = mapping[bond.to];
        if (target.bondOrder(first, second) != 1 || target.bondAromatic(first, second)) continue;
        // Bits of the directions the bond may take: 1 up, 2 down and 4 none, for a double bond
        // without a configuration.
        int allowed = 0;
        if (smarts_detail::evalBondExprWithDirection(bond.expr, target, first, second, 1)) allowed |= 1;
        if (smarts_detail::evalBondExprWithDirection(bond.expr, target, first, second, -1)) allowed |= 2;
        if (smarts_detail::evalBondExprWithoutDirection(bond.expr, target, first, second) == smarts_detail::StereoTruth::TRUE_VALUE)
            allowed |= 4;
        if (allowed == 0) return false;
        variableForBond[bondIndex] = static_cast<int>(allowedDirections.size());
        allowedDirections.push_back(allowed);
    }
    if (allowedDirections.empty()) return true;
    std::vector<std::vector<std::pair<int, int>>> constraints(allowedDirections.size());
    auto addConstraint = [&](int firstBond, int secondBond, int relation) {
        const int firstVariable = variableForBond[firstBond], secondVariable = variableForBond[secondBond];
        constraints[firstVariable].push_back({secondVariable, relation});
        constraints[secondVariable].push_back({firstVariable, relation});
    };
    auto outwardOrientation = [&](int bondIndex, int atom) { return bonds[bondIndex].from == atom ? 1 : -1; };
    // A pair of directional bonds constrains the bond between their atoms, as Daylight SMARTS and
    // Java read it: a double bond takes the configuration they give it, and an aromatic bond or a
    // double bond without a configuration matches only when either bond is written with ?.
    for (const auto& focus : bonds) {
        const int first = focus.from, second = focus.to;
        const int mappedFirst = mapping[first], mappedSecond = mapping[second];
        const int order = target.bondOrder(mappedFirst, mappedSecond);
        const bool aromatic = order == 4 || target.bondAromatic(mappedFirst, mappedSecond);
        if ((order != 2 && !aromatic) || smarts_detail::hasDirectionalPrimitive(focus.expr)) continue;
        std::array<std::vector<int>, 2> carriers;
        for (int end = 0; end < 2; ++end) {
            const int atom = end == 0 ? first : second, other = end == 0 ? second : first;
            for (size_t position = 0; position < neighbors[atom].size(); ++position) {
                const int bondIndex = bondIndices[atom][position];
                if (neighbors[atom][position] != other && variableForBond[bondIndex] >= 0)
                    carriers[end].push_back(bondIndex);
            }
            for (size_t previous = 0; previous < carriers[end].size(); ++previous)
                for (size_t next = previous + 1; next < carriers[end].size(); ++next)
                    addConstraint(carriers[end][previous], carriers[end][next],
                        -outwardOrientation(carriers[end][previous], atom) * outwardOrientation(carriers[end][next], atom));
        }
        if (carriers[0].empty() || carriers[1].empty()) continue;
        const int configuration = aromatic ? 0 : target.dbStereo(mappedFirst, mappedSecond);
        if (configuration == 0) {
            bool unspecified = false;
            for (const auto& end : carriers)
                for (int bondIndex : end) unspecified |= (allowedDirections[variableForBond[bondIndex]] & 4) != 0;
            if (!unspecified) return false;
            continue;
        }
        const int firstReference = stereo_detail::doubleBondReference(target, mappedFirst, mappedSecond);
        const int secondReference = stereo_detail::doubleBondReference(target, mappedSecond, mappedFirst);
        for (int firstBond : carriers[0]) {
            const auto& firstCarrier = bonds[firstBond];
            const int firstLigand = firstCarrier.from == first ? firstCarrier.to : firstCarrier.from;
            for (int secondBond : carriers[1]) {
                const auto& secondCarrier = bonds[secondBond];
                const int secondLigand = secondCarrier.from == second ? secondCarrier.to : secondCarrier.from;
                const int relation = (configuration == 1 ? 1 : -1)
                    * (mapping[firstLigand] == firstReference ? 1 : -1)
                    * (mapping[secondLigand] == secondReference ? 1 : -1)
                    * outwardOrientation(firstBond, first) * outwardOrientation(secondBond, second);
                addConstraint(firstBond, secondBond, relation);
            }
        }
    }
    std::vector<int> relativeDirections(allowedDirections.size(), 0);
    for (int start = 0; start < static_cast<int>(allowedDirections.size()); ++start) {
        if (relativeDirections[start] != 0) continue;
        relativeDirections[start] = 1;
        int allowedRoot = 3;
        std::deque<int> pending{start};
        while (!pending.empty()) {
            const int current = pending.front(); pending.pop_front();
            const int allowed = allowedDirections[current] & 3;
            allowedRoot &= relativeDirections[current] == 1 ? allowed : ((allowed & 1) << 1) | ((allowed & 2) >> 1);
            if (allowedRoot == 0) return false;
            for (const auto& constraint : constraints[current]) {
                const int required = relativeDirections[current] * constraint.second;
                if (relativeDirections[constraint.first] == 0) {
                    relativeDirections[constraint.first] = required;
                    pending.push_back(constraint.first);
                } else if (relativeDirections[constraint.first] != required) return false;
            }
        }
    }
    return true;
}

inline bool SmartsQuery::componentsCompatible(const MolGraph& target, const std::vector<int>& mapping) const {
    if (componentGroups.empty()) return true;
    auto& cache = smarts_detail::getSmartsTargetPropertyCache(target);
    std::vector<int> groupComponent;
    for (int atom = 0; atom < atomCount() && atom < static_cast<int>(componentGroups.size()); ++atom) {
        const int group = componentGroups[atom];
        if (group < 0 || mapping[atom] < 0) continue;
        if (group >= static_cast<int>(groupComponent.size())) groupComponent.resize(group + 1, -1);
        const int component = cache.componentOf(target, mapping[atom]);
        if (groupComponent[group] < 0) groupComponent[group] = component;
        else if (groupComponent[group] != component) return false;
    }
    std::vector<int> used;
    for (int component : groupComponent) if (component >= 0) used.push_back(component);
    std::sort(used.begin(), used.end());
    return std::adjacent_find(used.begin(), used.end()) == used.end();
}

inline bool SmartsQuery::atomCompat(int qi, const MolGraph& target, int ti) const {
    return smarts_detail::evalAtomExpr(atoms[qi].expr, target, ti);
}

inline bool SmartsQuery::bondCompat(int bi, const MolGraph& target, int ti, int tj) const {
    const auto& expression = bonds[bi].expr;
    if (smarts_detail::hasDirectionalPrimitive(expression))
        return smarts_detail::evalBondExprWithDirection(expression, target, ti, tj, 1)
            || smarts_detail::evalBondExprWithDirection(expression, target, ti, tj, -1);
    return smarts_detail::evalBondExpr(expression, target, ti, tj);
}

inline std::vector<int> SmartsQuery::matchingOrder() const {
    int n = atomCount();
    if (n == 0) return {};
    std::vector<int> order;
    order.reserve(n);
    std::vector<bool> visited(n, false);
    std::deque<int> queue;
    queue.push_back(0);
    visited[0] = true;
    while (!queue.empty()) {
        int u = queue.front(); queue.pop_front();
        order.push_back(u);
        for (int v : neighbors[u]) {
            if (!visited[v]) {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            visited[i] = true;
            queue.push_back(i);
            while (!queue.empty()) {
                int u = queue.front(); queue.pop_front();
                order.push_back(u);
                for (int v : neighbors[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
        }
    }
    return order;
}

// Deferred recursive SMARTS evaluation (needs complete SmartsQuery type)
inline bool smarts_detail::evalRecursiveSmarts(const smarts_detail::AtomExpr& expr, const MolGraph& g, int idx) {
    if (expr.recursiveSmarts.empty()) return false;

    if (smarts_detail::recursiveMatcher) {
        if (expr.cachedRecursive) return smarts_detail::recursiveMatcher(*expr.cachedRecursive, g, idx);
        const SmartsQuery transient = parseSMARTS(expr.recursiveSmarts);
        struct Depth {
            Depth() { ++smarts_detail::transientRecursionDepth; }
            ~Depth() { --smarts_detail::transientRecursionDepth; }
        } depth;
        return smarts_detail::recursiveMatcher(transient, g, idx);
    }
    throw std::logic_error("Recursive SMARTS evaluation requires a Substructure search context");
}

} // namespace smsd

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Graphs of selected matches: Match::toGraph(), Match::toSmiles() and the
 * connected components of a graph.
 *
 * Included by search_types.hpp.
 */
#pragma once

#include "smsd/search_types.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/query_smarts.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {

inline MolGraph Match::toGraph(const MolGraph& query) const {
    const std::lock_guard<std::recursive_mutex> lock(query.lazyCacheMutex());
    std::vector<int> originals;
    std::vector<int> reindex(query.n, -1);
    for (const auto& atom : atomMapping) {
        if (atom.first < 0 || atom.first >= query.n)
            throw std::out_of_range("Selected atom index out of range");
        reindex[atom.first] = static_cast<int>(originals.size());
        originals.push_back(atom.first);
    }
    std::set<int64_t> selected;
    for (const auto& bond : bondMapping) {
        if (bond[0] < 0 || bond[1] < 0 || bond[0] >= query.n || bond[1] >= query.n)
            throw std::out_of_range("Selected bond atom index out of range");
        if (reindex[bond[0]] < 0 || reindex[bond[1]] < 0 || !query.hasBond(bond[0], bond[1])
            || atomMapping.at(bond[0]) != bond[2] || atomMapping.at(bond[1]) != bond[3])
            throw std::invalid_argument("Selected bond is not in the atom mapping");
        if (!selected.insert(MolGraph::bondKey(bond[0], bond[1])).second)
            throw std::invalid_argument("Selected bond occurs more than once");
    }
    const int count = static_cast<int>(originals.size());
    std::vector<int> numbers(count), charges(count), masses(count), hydrogens(count), classes(count), chirality(count), ids;
    std::vector<uint8_t> rings(count), aromatic(count);
    std::vector<std::vector<int>> adjacent(count), orders(count);
    std::vector<std::vector<bool>> bondRings(count), bondAromatic(count);
    std::vector<bool> allLigands(count, true);
    if (!query.atomId.empty()) ids.resize(count);
    auto property = [](const auto& values, int atom, int fallback = 0) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : fallback;
    };
    for (int atom = 0; atom < count; ++atom) {
        const int original = originals[atom];
        numbers[atom] = query.atomicNum[original]; charges[atom] = property(query.formalCharge, original);
        masses[atom] = property(query.massNumber, original); hydrogens[atom] = property(query.hydrogenCount, original);
        classes[atom] = property(query.atomClass, original); rings[atom] = query.ring[original]; aromatic[atom] = query.aromatic[original];
        if (!ids.empty()) ids[atom] = query.atomId[original];
        for (int neighbor : query.neighbors[original]) {
            if (reindex[neighbor] < 0 || selected.count(MolGraph::bondKey(original, neighbor)) == 0) {
                allLigands[atom] = false;
                continue;
            }
            adjacent[atom].push_back(reindex[neighbor]); orders[atom].push_back(query.bondOrder(original, neighbor));
            bondRings[atom].push_back(query.bondInRing(original, neighbor));
            bondAromatic[atom].push_back(query.bondAromatic(original, neighbor));
        }
        if (allLigands[atom]) chirality[atom] = property(query.tetraChirality, original);
    }
    std::vector<std::tuple<int, int, int>> stereo;
    for (const auto& bond : bondMapping) {
        const int first = reindex[bond[0]], second = reindex[bond[1]];
        if (!allLigands[first] || !allLigands[second]) continue;
        const int configuration = query.dbStereo(bond[0], bond[1]);
        if (configuration != 0) stereo.emplace_back(first, second, configuration);
    }
    auto builder = MolGraph::Builder().atomCount(count).atomicNumbers(std::move(numbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses)).hydrogenCounts(std::move(hydrogens))
        .atomClasses(std::move(classes)).ringFlags(std::move(rings)).aromaticFlags(std::move(aromatic))
        .tetrahedralChirality(std::move(chirality)).setNeighbors(std::move(adjacent)).setBondOrders(std::move(orders))
        .bondRingFlags(std::move(bondRings)).bondAromaticFlags(std::move(bondAromatic))
        .doubleBondStereoSparse(std::move(stereo)).solvent(query.solvent_).pH(query.pH_);
    if (!ids.empty()) builder.atomIds(std::move(ids));
    MolGraph graph = builder.build(false);
    graph.hydrogenCountsKnown = query.hydrogenCountsKnown;
    query.ensureRingCounts();
    for (int atom = 0; atom < count; ++atom) {
        const int original = originals[atom];
        graph.ringCount[atom] = query.ringCount[original];
        graph.tautomerClass[atom] = property(query.tautomerClass, original, -1);
        graph.tautomerWeight[atom] = original < static_cast<int>(query.tautomerWeight.size()) ? query.tautomerWeight[original] : 1.0f;
        graph.tautomerSite[atom] = static_cast<int8_t>(property(query.tautomerSite, original, MolGraph::TAUTOMER_SITE_NONE));
        graph.tautomerPKa[atom] = original < static_cast<int>(query.tautomerPKa.size()) ? query.tautomerPKa[original] : 0.0f;
    }
    graph.ringCountsComputed_ = true;
    return graph;
}

inline std::string Match::toSmarts(const MolGraph& query, const ChemOptions& chemistry) const {
    return toQuerySMARTS(toGraph(query), chemistry);
}

inline std::string Match::toSmiles(const MolGraph& g) const {
    const auto& mapping = atomMapping;
    for (const auto& atom : mapping)
        if (atom.first < 0 || atom.first >= g.n) throw std::out_of_range("Selected atom index out of range");
    std::set<int64_t> selectedBonds;
    for (const auto& bond : bondMapping) {
        if (bond[0] < 0 || bond[1] < 0 || bond[0] >= g.n || bond[1] >= g.n)
            throw std::out_of_range("Selected bond atom index out of range");
        if (mapping.count(bond[0]) == 0 || mapping.count(bond[1]) == 0
            || !g.hasBond(bond[0], bond[1]) || mapping.at(bond[0]) != bond[2] || mapping.at(bond[1]) != bond[3])
            throw std::invalid_argument("Selected bond is not in the query mapping");
        if (!selectedBonds.insert(MolGraph::bondKey(bond[0], bond[1])).second)
            throw std::invalid_argument("Selected bond occurs more than once");
    }
    if (mapping.empty()) return "";
    for (const auto& pair : mapping)
        if (pair.first < 0 || pair.first >= g.n)
            throw std::out_of_range("MCS atom index out of range");

    std::vector<int> atoms;
    atoms.reserve(mapping.size());
    for (const auto& pair : mapping) atoms.push_back(pair.first);
    const int atomCount = static_cast<int>(atoms.size());
    std::unordered_map<int, int> reindex;
    for (int atom = 0; atom < atomCount; ++atom) reindex[atoms[atom]] = atom;

    std::vector<int> atomicNumbers(atomCount), charges(atomCount), masses(atomCount), hydrogens(atomCount);
    std::vector<int> chirality(atomCount, 0), classes(atomCount, 0);
    std::vector<uint8_t> aromatic(atomCount, 0), ring(atomCount, 0);
    std::vector<std::vector<int>> neighbors(atomCount), orders(atomCount);
    std::vector<std::vector<bool>> bondRings(atomCount), bondAromatic(atomCount);
    std::vector<bool> allLigandsRetained(atomCount, true);
    auto property = [](const auto& values, int atom) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : 0;
    };
    for (int atom = 0; atom < atomCount; ++atom) {
        const int original = atoms[atom];
        atomicNumbers[atom] = g.atomicNum[original];
        charges[atom] = property(g.formalCharge, original);
        classes[atom] = property(g.atomClass, original);
        masses[atom] = property(g.massNumber, original);
        hydrogens[atom] = detail::hydrogenCountForWrite(g, original);
        for (int neighbor : g.neighbors[original]) {
            const auto retained = reindex.find(neighbor);
            if (retained == reindex.end() || selectedBonds.count(MolGraph::bondKey(original, neighbor)) == 0) {
                allLigandsRetained[atom] = false;
                continue;
            }
            neighbors[atom].push_back(retained->second);
            orders[atom].push_back(g.bondOrder(original, neighbor));
        }
        if (allLigandsRetained[atom]) chirality[atom] = property(g.tetraChirality, original);
    }

    std::vector<bool> retainedRingAtoms;
    std::vector<std::set<int>> retainedRingBonds;
    detail::detectRings(atomCount, neighbors, retainedRingAtoms, retainedRingBonds);
    for (int atom = 0; atom < atomCount; ++atom) {
        ring[atom] = retainedRingAtoms[atom];
        for (int bondIndex = 0; bondIndex < static_cast<int>(neighbors[atom].size()); ++bondIndex) {
            const int neighbor = neighbors[atom][bondIndex];
            const bool inRing = retainedRingBonds[atom].count(neighbor) != 0;
            const bool isAromatic = inRing && g.bondAromatic(atoms[atom], atoms[neighbor]);
            bondRings[atom].push_back(inRing);
            bondAromatic[atom].push_back(isAromatic);
            if (isAromatic) aromatic[atom] = 1;
            else if (orders[atom][bondIndex] == 4) orders[atom][bondIndex] = 1;
        }
        if (!allLigandsRetained[atom] || aromatic[atom] != g.aromatic[atoms[atom]]) {
            int bondOrderSum = 0;
            for (int order : orders[atom]) bondOrderSum += order == 4 ? 1 : order;
            if (aromatic[atom] && (atomicNumbers[atom] == 7 || atomicNumbers[atom] == 15))
                hydrogens[atom] += static_cast<int>(g.neighbors[atoms[atom]].size() - neighbors[atom].size());
            else
                hydrogens[atom] = atomicNumbers[atom] <= 1 ? 0
                    : detail::computeImplicitH(atomicNumbers[atom], aromatic[atom], bondOrderSum, charges[atom]);
        }
    }

    std::vector<std::tuple<int, int, int>> doubleBondStereo;
    for (int atom = 0; atom < atomCount; ++atom) {
        if (!allLigandsRetained[atom]) continue;
        for (int neighbor : neighbors[atom]) {
            if (atom >= neighbor || !allLigandsRetained[neighbor]) continue;
            const int configuration = g.dbStereo(atoms[atom], atoms[neighbor]);
            if (configuration == 0) continue;
            doubleBondStereo.emplace_back(atom, neighbor, configuration);
        }
    }
    MolGraph sub = MolGraph::Builder()
        .atomCount(atomCount).atomicNumbers(std::move(atomicNumbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses)).atomClasses(std::move(classes))
        .hydrogenCounts(std::move(hydrogens)).tetrahedralChirality(std::move(chirality))
        .ringFlags(std::move(ring)).aromaticFlags(std::move(aromatic))
        .setNeighbors(std::move(neighbors)).setBondOrders(std::move(orders))
        .bondRingFlags(std::move(bondRings)).bondAromaticFlags(std::move(bondAromatic))
        .doubleBondStereoSparse(std::move(doubleBondStereo)).build();
    if (!sub.kekulize()) {
        for (int atom = 0; atom < sub.n; ++atom) {
            if (!sub.aromatic[atom]) continue;
            sub.aromatic[atom] = 0;
            int bondOrderSum = 0;
            for (int neighbor : sub.neighbors[atom]) {
                if (sub.bondOrder(atom, neighbor) == 4) sub.setBondOrder(atom, neighbor, 1);
                sub.setBondAromatic(atom, neighbor, false);
                bondOrderSum += sub.bondOrder(atom, neighbor);
            }
            sub.hydrogenCount[atom] = detail::computeImplicitH(
                sub.atomicNum[atom], false, bondOrderSum, sub.formalCharge[atom]);
        }
    }
    sub.perceiveAromaticity();

    return toSMILES(sub);
}

namespace search_detail {

inline std::vector<std::vector<int>> componentAtoms(const MolGraph& graph) {
    std::vector<bool> visited(graph.n, false);
    std::vector<std::vector<int>> components;
    for (int atom = 0; atom < graph.n; ++atom) {
        if (visited[atom]) continue;
        visited[atom] = true;
        std::vector<int> component{atom};
        for (size_t head = 0; head < component.size(); ++head)
            for (int neighbor : graph.neighbors[component[head]])
                if (!visited[neighbor]) { visited[neighbor] = true; component.push_back(neighbor); }
        std::sort(component.begin(), component.end());
        components.push_back(std::move(component));
    }
    return components;
}

inline int countComponents(const MolGraph& graph) {
    return static_cast<int>(componentAtoms(graph).size());
}

inline std::vector<MolGraph> splitComponents(const MolGraph& graph) {
    std::vector<MolGraph> result;
    for (const auto& component : componentAtoms(graph)) {
        Match selected;
        for (int atom : component) selected.atomMapping.emplace(atom, atom);
        for (int atom : component) for (int neighbor : graph.neighbors[atom])
            if (atom < neighbor) selected.bondMapping.push_back({atom, neighbor, atom, neighbor});
        result.push_back(selected.toGraph(graph));
    }
    return result;
}

inline bool sameCanonicalGraph(const MolGraph& first, const MolGraph& second) {
    return first.n == second.n && first.getCanonicalHash() == second.getCanonicalHash()
        && first.toCanonicalSmiles() == second.toCanonicalSmiles();
}

} // namespace search_detail

} // namespace smsd

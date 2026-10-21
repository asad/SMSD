/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph aromaticity: candidate cycles, pi-electron counts, ring and aromatic
 * flags, and Kekule assignment by matching the double-bond demands of atoms.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// Aromaticity perception
// ============================================================================

inline bool MolGraph::isAromaticDonorElement(int z) {
    return z == 7 || z == 8 || z == 15 || z == 16 || z == 33 || z == 34
        || z == 51 || z == 52;
}

inline std::pair<bool, bool> MolGraph::exocyclicMultipleBondInfo(
    int atom,
    int prev,
    int next) const {
    bool hasExocyclicMultipleBond = false;
    bool hasExocyclicMultipleBondToHetero = false;
    for (int nb : neighbors[atom]) {
        if (nb == prev || nb == next) continue;
        int ord = bondOrder(atom, nb);
        if (ord < 2 || ord == 4) continue;
        hasExocyclicMultipleBond = true;
        if (atomicNum[nb] != 6) hasExocyclicMultipleBondToHetero = true;
    }
    return {hasExocyclicMultipleBond, hasExocyclicMultipleBondToHetero};
}

inline bool MolGraph::isSimpleCycle(const std::vector<int>& cycle) {
    if (cycle.size() < 3) return false;
    std::unordered_set<int> seen;
    seen.reserve(cycle.size());
    for (int atom : cycle) {
        if (!seen.insert(atom).second) return false;
    }
    return true;
}

inline int MolGraph::aromaticPiContribution(const std::vector<int>& cycle, int pos) const {
    const int len = static_cast<int>(cycle.size());
    const int atom = cycle[pos];
    const int prev = cycle[(pos + len - 1) % len];
    const int next = cycle[(pos + 1) % len];

    int cyclePiBonds = 0;
    bool hasExplicitAromaticBond = false;
    for (int nb : {prev, next}) {
        int ord = bondOrder(atom, nb);
        if (ord == 3) return -1;
        if (ord == 4) {
            hasExplicitAromaticBond = true;
        } else if (ord == 2) {
            ++cyclePiBonds;
        }
    }

    if (hasExplicitAromaticBond) {
        if (atomicNum[atom] == 6) {
            for (int neighbor : neighbors[atom]) {
                if (bondOrder(atom, neighbor) == 2 && !bondInRing(atom, neighbor))
                    return atomicNum[neighbor] == 6 ? 1 : 0;
            }
        }
        return isAromaticDonorElement(atomicNum[atom]) && aromaticKekuleDemand(atom) == 0 ? 2 : 1;
    }
    if (cyclePiBonds == 1) return 1;
    if (cyclePiBonds > 1) return -1;

    // A fused-ring double bond remains cyclic across candidate-cycle boundaries.
    for (int neighbor : neighbors[atom]) {
        if (neighbor != prev && neighbor != next && bondOrder(atom, neighbor) == 2
            && bondInRing(atom, neighbor)) return 1;
    }

    auto [hasExocyclicMultipleBond, hasExocyclicMultipleBondToHetero]
        = exocyclicMultipleBondInfo(atom, prev, next);

    const int z = atomicNum[atom];
    const int q = formalCharge.empty() ? 0 : formalCharge[atom];

    if (z == 6) {
        if (q < 0 && !hasExocyclicMultipleBond) return 2;
        if (q == 0 && hasExocyclicMultipleBondToHetero) return 0;
        if (q == 0)
            for (int neighbor : neighbors[atom]) {
                if (neighbor != prev && neighbor != next && atomicNum[neighbor] == 6
                    && bondOrder(atom, neighbor) == 2) return 1;
            }
        if (q > 0) return 0;
        return -1;
    }

    if (z == 5 || z == 13) {
        if (q < 0 && !hasExocyclicMultipleBond) return 2;
        return -1;
    }

    if (isAromaticDonorElement(z)) {
        if (q > 0) {
            if (hasExocyclicMultipleBond) return 0;
            const int hydrogens = hydrogenCount.empty() ? 0 : hydrogenCount[atom];
            const int64_t sigmaBonds = int64_t(neighbors[atom].size()) + hydrogens;
            return sigmaBonds <= 2 ? 0 : -1;
        }
        if (hasExocyclicMultipleBond) return -1;
        return 2;
    }

    return -1;
}

inline bool MolGraph::cycleLooksAromatic(const std::vector<int>& cycle) const {
    const int len = static_cast<int>(cycle.size());
    if (len < 3) return false;

    int piElectrons = 0;
    bool sawCyclePi = false;
    for (int pos = 0; pos < len; ++pos) {
        int contrib = aromaticPiContribution(cycle, pos);
        if (contrib < 0) return false;
        if (contrib == 1 || contrib == 2) sawCyclePi = true;
        piElectrons += contrib;
    }

    if (!sawCyclePi || piElectrons < 2) return false;
    return ((piElectrons - 2) % 4) == 0;
}

inline void MolGraph::enumerateRingCyclesDfs(
    int root,
    int current,
    const std::vector<std::vector<int>>& ringAdj,
    const std::vector<uint8_t>& inComponent,
    std::vector<uint8_t>& inPath,
    std::vector<int>& path,
    std::set<std::vector<int64_t>>& seenCycles,
    std::vector<std::vector<int>>& cycles,
    size_t maxCycles,
    uint64_t& steps) const {
    if (cycles.size() >= maxCycles || steps >= AROMATIC_CYCLE_STEP_LIMIT) return;
    struct Frame { int atom; size_t nextNeighbor; };
    std::vector<Frame> pending{{current, 0}};
    const size_t initialPathSize = path.size();
    while (!pending.empty() && cycles.size() < maxCycles && ++steps < AROMATIC_CYCLE_STEP_LIMIT) {
        auto& frame = pending.back();
        if (frame.nextNeighbor == ringAdj[frame.atom].size()) {
            pending.pop_back();
            if (path.size() > initialPathSize) {
                inPath[path.back()] = 0;
                path.pop_back();
            }
            continue;
        }
        const int nb = ringAdj[frame.atom][frame.nextNeighbor++];
        if (!inComponent[nb] || nb < root) continue;
        if (nb == root) {
            if (path.size() >= 3) {
                steps += path.size();
                auto key = canonicalCycleKey(path);
                if (seenCycles.insert(key).second) cycles.push_back(path);
            }
            continue;
        }
        if (inPath[nb]) continue;

        inPath[nb] = 1;
        path.push_back(nb);
        pending.push_back({nb, 0});
    }
    while (path.size() > initialPathSize) {
        inPath[path.back()] = 0;
        path.pop_back();
    }
}

inline std::vector<std::vector<int>> MolGraph::collectAromaticCandidateCycles() const {
    std::vector<std::vector<int>> cycles;
    std::set<std::vector<int64_t>> seenCycles;
    uint64_t steps = 0;

    auto addCycle = [&](const std::vector<int>& cycle) {
        if (!isSimpleCycle(cycle)) return;
        if (cycle.size() < 3) return;
        auto key = canonicalCycleKey(cycle);
        if (seenCycles.insert(key).second) cycles.push_back(cycle);
    };

    for (const auto& cycle : computeRelevantCycles()) addCycle(cycle);

    std::vector<std::vector<int>> ringAdj(n);
    for (int i = 0; i < n; ++i) {
        if (!ring[i]) continue;
        for (int nb : neighbors[i]) {
            if (i < nb && bondInRing(i, nb)) {
                ringAdj[i].push_back(nb);
                ringAdj[nb].push_back(i);
            }
        }
    }

    std::vector<uint8_t> visited(n, uint8_t(0));
    std::vector<uint8_t> inComponent(n, uint8_t(0));
    std::vector<uint8_t> inPath(n, uint8_t(0));
    for (int start = 0; start < n; ++start) {
        if (!ring[start] || visited[start] || ringAdj[start].empty()) continue;

        std::vector<int> component;
        std::deque<int> queue;
        queue.push_back(start);
        visited[start] = 1;
        while (!queue.empty()) {
            int atom = queue.front();
            queue.pop_front();
            component.push_back(atom);
            for (int nb : ringAdj[atom]) {
                if (visited[nb]) continue;
                visited[nb] = 1;
                queue.push_back(nb);
            }
        }

        int edgeCount = 0;
        for (int atom : component) edgeCount += static_cast<int>(ringAdj[atom].size());
        edgeCount /= 2;
        const int cycleRank = edgeCount - static_cast<int>(component.size()) + 1;
        if (cycleRank <= 1) continue;

        std::fill(inComponent.begin(), inComponent.end(), uint8_t(0));
        for (int atom : component) inComponent[atom] = 1;

        std::sort(component.begin(), component.end());
        for (int root : component) {
            if (ringAdj[root].size() < 2) continue;
            std::fill(inPath.begin(), inPath.end(), uint8_t(0));
            steps += static_cast<uint64_t>(n);
            std::vector<int> path = { root };
            inPath[root] = 1;

            for (int nb : ringAdj[root]) {
                if (!inComponent[nb] || nb <= root) continue;
                inPath[nb] = 1;
                path.push_back(nb);
                enumerateRingCyclesDfs(
                    root, nb, ringAdj, inComponent, inPath, path,
                    seenCycles, cycles, 4096, steps);
                path.pop_back();
                inPath[nb] = 0;
                if (cycles.size() >= 4096 || steps >= AROMATIC_CYCLE_STEP_LIMIT) break;
            }
            inPath[root] = 0;
            if (cycles.size() >= 4096 || steps >= AROMATIC_CYCLE_STEP_LIMIT) break;
        }
        if (cycles.size() >= 4096 || steps >= AROMATIC_CYCLE_STEP_LIMIT) break;
    }

    return cycles;
}

inline void MolGraph::clearRingFlagsFromTopology() {
    std::fill(ring.begin(), ring.end(), uint8_t(0));
    if (useDense) {
        for (auto& row : bondPacked) {
            for (auto& b : row) b &= ~BOND_RING_BIT;
        }
    } else {
        for (auto& kv : sparseBondProps) {
            kv.second[1] = 0;
        }
    }
}

inline void MolGraph::markRingEdge(int i, int j) {
    ring[i] = 1;
    ring[j] = 1;
    if (useDense) {
        bondPacked[i][j] |= BOND_RING_BIT;
        bondPacked[j][i] |= BOND_RING_BIT;
    } else {
        auto it = sparseBondProps.find(bondKey(i, j));
        if (it != sparseBondProps.end()) it->second[1] = 1;
    }
}

inline void MolGraph::markAromaticEdge(int i, int j) {
    if (useDense) {
        bondPacked[i][j] = (bondPacked[i][j] & BOND_RING_BIT) | 4 | BOND_AROM_BIT;
        bondPacked[j][i] = bondPacked[i][j];
    } else {
        auto it = sparseBondProps.find(bondKey(i, j));
        if (it != sparseBondProps.end()) {
            it->second[0] = 4;
            it->second[2] = 1;
        }
    }
}

inline void MolGraph::markKekuleEdge(int i, int j, int order) {
    if (useDense) {
        bondPacked[i][j] = (bondPacked[i][j] & BOND_RING_BIT) | static_cast<uint8_t>(order & 0x07);
        bondPacked[j][i] = bondPacked[i][j];
    } else {
        auto it = sparseBondProps.find(bondKey(i, j));
        if (it != sparseBondProps.end()) {
            it->second[0] = order;
            it->second[2] = 0;
        }
    }
}

inline void MolGraph::refreshAtomLabels() {
    label.resize(n);
    for (int i = 0; i < n; ++i) {
        label[i] = (atomicNum[i] << 2)
                 | (aromatic[i] ? 2 : 0)
                 | (ring[i] ? 1 : 0);
    }
}

inline void MolGraph::normalizeRingAndAromaticity(
    AromaticityModel model) {
    if (n == 0) return;
    (void) model;

    clearRingFlagsFromTopology();
    const auto& sssr = computeRings();
    for (const auto& cycle : sssr) {
        for (size_t i = 0; i < cycle.size(); ++i) {
            int a = cycle[i];
            int b = cycle[(i + 1) % cycle.size()];
            markRingEdge(a, b);
        }
    }

    std::vector<uint8_t> perceivedAromatic = aromatic;
    std::vector<std::pair<int, int>> aromaticEdges;
    for (int i = 0; i < n; ++i) {
        if (!ring[i]) continue;
        for (int nb : neighbors[i]) {
            if (i < nb && bondInRing(i, nb) && bondAromatic(i, nb)) {
                aromaticEdges.emplace_back(i, nb);
            }
        }
    }

    auto candidateCycles = collectAromaticCandidateCycles();
    for (const auto& cycle : candidateCycles) {
        if (!cycleLooksAromatic(cycle)) continue;
        for (int atom : cycle) perceivedAromatic[atom] = 1;
        for (size_t i = 0; i < cycle.size(); ++i) {
            int a = cycle[i];
            int b = cycle[(i + 1) % cycle.size()];
            aromaticEdges.emplace_back(std::min(a, b), std::max(a, b));
        }
    }

    for (int i = 0; i < n; ++i) {
        if (perceivedAromatic[i]) aromatic[i] = 1;
    }

    std::sort(aromaticEdges.begin(), aromaticEdges.end());
    aromaticEdges.erase(std::unique(aromaticEdges.begin(), aromaticEdges.end()), aromaticEdges.end());
    for (const auto& edge : aromaticEdges) {
        markAromaticEdge(edge.first, edge.second);
    }
}

inline void MolGraph::perceiveAromaticity(
    AromaticityModel model) {
    normalizeRingAndAromaticity(model);
    refreshAtomLabels();
    ringCount.assign(n, 0);
    ringCountsComputed_ = false;
    invalidateBondDerivedState();
}

// ============================================================================
// Kekule assignment
// ============================================================================

inline int MolGraph::aromaticKekuleDemand(int atom) const {
    if (!aromatic[atom]) return 0;
    const int z = atomicNum[atom];
    const int q = formalCharge.empty() ? 0 : formalCharge[atom];
    const int h = hydrogenCount.empty() ? 0 : hydrogenCount[atom];
    const bool hasExocyclicMultipleBond = exocyclicMultipleBondInfo(atom).first;

    if (z == 6) {
        if (q == 0 && hasExocyclicMultipleBond) return 0;
        return q == 0 ? 1 : 0;
    }
    if (z == 5 || z == 13) return q >= 0 ? 1 : 0;
    if (z == 7 || z == 15) {
        if (h > 0 && q <= 0) return 0;
        if (q == 0 && neighbors[atom].size() >= 3) return 0;
        if (q < 0) return 0;
        return 1;
    }
    if (z == 8 || z == 16 || z == 34 || z == 52)
        return q == 1 && int64_t(neighbors[atom].size()) + h == 2 ? 1 : 0;
    return 1;
}

inline bool MolGraph::findKekuleAssignment(const std::vector<int>& component,
                         std::vector<std::pair<int, int>>& singleBonds,
                         std::vector<std::pair<int, int>>& doubleBonds) const {
    if (component.empty()) return true;

    std::unordered_set<int> inComponent(component.begin(), component.end());
    std::vector<std::pair<int, int>> aromaticEdges;
    aromaticEdges.reserve(component.size() * 2);
    for (int atom : component) {
        for (int nb : neighbors[atom]) {
            if (atom < nb && inComponent.count(nb)
                && (bondAromatic(atom, nb) || bondOrder(atom, nb) == 4)) {
                aromaticEdges.emplace_back(atom, nb);
            }
        }
    }
    if (aromaticEdges.empty()) {
        return true;
    }

    std::vector<std::pair<int, int>> reservedDoubleBonds;
    std::unordered_set<int> reservedAtoms;
    for (const auto& edge : aromaticEdges) {
        if (dbStereo(edge.first, edge.second) == 0) continue;
        if (aromaticKekuleDemand(edge.first) <= 0 || aromaticKekuleDemand(edge.second) <= 0
            || reservedAtoms.count(edge.first) || reservedAtoms.count(edge.second)) return false;
        reservedAtoms.insert(edge.first);
        reservedAtoms.insert(edge.second);
        reservedDoubleBonds.push_back(edge);
    }

    std::vector<int> demandAtoms;
    demandAtoms.reserve(component.size());
    std::unordered_map<int, int> localIndex;
    for (int atom : component) {
        if (aromaticKekuleDemand(atom) > 0 && !reservedAtoms.count(atom)) {
            localIndex[atom] = static_cast<int>(demandAtoms.size());
            demandAtoms.push_back(atom);
        }
    }

    if (!demandAtoms.empty()) {
        if (demandAtoms.size() % 2 != 0) return false;
        std::vector<std::vector<int>> demandAdj(demandAtoms.size());
        for (const auto& edge : aromaticEdges) {
            auto ita = localIndex.find(edge.first);
            auto itb = localIndex.find(edge.second);
            if (ita == localIndex.end() || itb == localIndex.end()) continue;
            demandAdj[ita->second].push_back(itb->second);
            demandAdj[itb->second].push_back(ita->second);
        }

        std::vector<int> color(demandAtoms.size(), -1);
        bool bipartite = true;
        std::deque<int> queue;
        for (int start = 0; start < static_cast<int>(demandAtoms.size()); ++start) {
            if (color[start] != -1) continue;
            color[start] = 0;
            queue.push_back(start);
            while (!queue.empty()) {
                int u = queue.front();
                queue.pop_front();
                for (int v : demandAdj[u]) {
                    if (color[v] == -1) {
                        color[v] = color[u] ^ 1;
                        queue.push_back(v);
                    } else if (color[v] == color[u]) {
                        bipartite = false;
                    }
                }
            }
        }

        // Non-alternant aromatic systems such as azulene contain odd
        // cycles. They can still have a perfect matching of double-bond
        // demands; use blossom matching rather than rejecting the system.
        if (!bipartite) {
            const auto mate = detail::maximumMatching(demandAdj);
            if (std::any_of(mate.begin(), mate.end(), [](int v) { return v < 0; }))
                return false;
            singleBonds.insert(singleBonds.end(), aromaticEdges.begin(), aromaticEdges.end());
            doubleBonds.insert(doubleBonds.end(), reservedDoubleBonds.begin(), reservedDoubleBonds.end());
            for (int i = 0; i < static_cast<int>(mate.size()); ++i) {
                if (i < mate[i]) doubleBonds.emplace_back(demandAtoms[i], demandAtoms[mate[i]]);
            }
            return true;
        }

        std::vector<int> leftNodes, rightNodes;
        leftNodes.reserve(demandAtoms.size());
        rightNodes.reserve(demandAtoms.size());
        for (int i = 0; i < static_cast<int>(demandAtoms.size()); ++i) {
            if (color[i] == 0) leftNodes.push_back(i);
            else rightNodes.push_back(i);
        }
        if (leftNodes.size() != rightNodes.size()) return false;

        std::unordered_map<int, int> rightPos;
        rightPos.reserve(rightNodes.size());
        for (int i = 0; i < static_cast<int>(rightNodes.size()); ++i) {
            rightPos[rightNodes[i]] = i;
        }

        std::vector<std::vector<int>> matchAdj(leftNodes.size());
        for (int li = 0; li < static_cast<int>(leftNodes.size()); ++li) {
            int u = leftNodes[li];
            for (int v : demandAdj[u]) {
                auto it = rightPos.find(v);
                if (it != rightPos.end()) matchAdj[li].push_back(it->second);
            }
            if (matchAdj[li].empty()) return false;
        }

        std::vector<int> matchRight(rightNodes.size(), -1);
        std::function<bool(int, std::vector<uint8_t>&)> augment =
            [&](int li, std::vector<uint8_t>& seen) {
                for (int ri : matchAdj[li]) {
                    if (seen[ri]) continue;
                    seen[ri] = 1;
                    if (matchRight[ri] == -1 || augment(matchRight[ri], seen)) {
                        matchRight[ri] = li;
                        return true;
                    }
                }
                return false;
            };

        for (int li = 0; li < static_cast<int>(leftNodes.size()); ++li) {
            std::vector<uint8_t> seen(rightNodes.size(), uint8_t(0));
            if (!augment(li, seen)) return false;
        }

        singleBonds.insert(singleBonds.end(), aromaticEdges.begin(), aromaticEdges.end());
        doubleBonds.insert(doubleBonds.end(), reservedDoubleBonds.begin(), reservedDoubleBonds.end());

        for (int ri = 0; ri < static_cast<int>(matchRight.size()); ++ri) {
            int li = matchRight[ri];
            if (li < 0) return false;
            int a = demandAtoms[leftNodes[li]];
            int b = demandAtoms[rightNodes[ri]];
            doubleBonds.emplace_back(a, b);
        }
        return true;
    }

    singleBonds.insert(singleBonds.end(), aromaticEdges.begin(), aromaticEdges.end());
    doubleBonds.insert(doubleBonds.end(), reservedDoubleBonds.begin(), reservedDoubleBonds.end());
    return true;
}

inline bool MolGraph::kekulizeAromaticComponent(const std::vector<int>& component) {
    std::vector<std::pair<int, int>> singleBonds, doubleBonds;
    if (!findKekuleAssignment(component, singleBonds, doubleBonds)) return false;
    for (const auto& edge : singleBonds) markKekuleEdge(edge.first, edge.second, 1);
    for (const auto& edge : doubleBonds) markKekuleEdge(edge.first, edge.second, 2);
    for (int atom : component) aromatic[atom] = 0;
    return true;
}

inline bool MolGraph::kekulize() {
    std::vector<uint8_t> visited(n, uint8_t(0));
    std::vector<std::pair<int, int>> singleBonds, doubleBonds;
    for (int start = 0; start < n; ++start) {
        if (!aromatic[start] || visited[start]) continue;

        std::vector<int> component;
        std::deque<int> queue;
        queue.push_back(start);
        visited[start] = 1;
        while (!queue.empty()) {
            int atom = queue.front();
            queue.pop_front();
            component.push_back(atom);
            for (int nb : neighbors[atom]) {
                if (visited[nb] || !aromatic[nb]) continue;
                if (!(bondAromatic(atom, nb) || bondOrder(atom, nb) == 4)) continue;
                visited[nb] = 1;
                queue.push_back(nb);
            }
        }

        if (!findKekuleAssignment(component, singleBonds, doubleBonds)) return false;
    }

    for (const auto& edge : singleBonds) markKekuleEdge(edge.first, edge.second, 1);
    for (const auto& edge : doubleBonds) markKekuleEdge(edge.first, edge.second, 2);
    for (int atom = 0; atom < n; ++atom)
        if (visited[atom]) aromatic[atom] = 0;
    refreshAtomLabels();
    invalidateBondDerivedState();
    return true;
}

inline bool MolGraph::dearomatize() {
    return kekulize();
}

} // namespace smsd

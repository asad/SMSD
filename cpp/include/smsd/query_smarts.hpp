/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace smsd {

/** Write a selected graph as SMARTS under the representable molecular matching rules.
 * Stereo, tautomer, complete-ring and strict ring-fusion constraints need the original
 * search and are rejected rather than weakened. Atom numbering determines traversal, and
 * ring numbers are reused once closed; a graph that needs more than 99 simultaneously open
 * ring closures is rejected. Source and target bonds must use normalised molecular aromatic
 * annotations.
 */
inline std::string toQuerySMARTS(const MolGraph& graph, const ChemOptions& chemistry = {}) {
    if (chemistry.useChirality || chemistry.useBondStereo || chemistry.tautomerAware
        || chemistry.completeRingsOnly || chemistry.ringFusionMode == ChemOptions::RingFusionMode::STRICT)
        throw std::invalid_argument("SMARTS export cannot represent these stereo, tautomer or ring constraints");
    const std::lock_guard<std::recursive_mutex> lock(graph.lazyCacheMutex());
    auto atomText = [&](int atom) {
        const std::string element = chemistry.matchAtomType ? "#" + std::to_string(graph.atomicNum[atom]) : "*";
        std::string text = element;
        const int isotope = atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0;
        if (chemistry.matchIsotope && isotope != 0)
            text = "0" + element + "," + std::to_string(isotope) + element;
        if (chemistry.aromaticityMode == ChemOptions::AromaticityMode::STRICT)
            text += graph.aromatic[atom] ? ";a" : ";A";
        if (chemistry.matchFormalCharge) {
            const int charge = atom < static_cast<int>(graph.formalCharge.size()) ? graph.formalCharge[atom] : 0;
            if (charge < -15 || charge > 15) throw std::invalid_argument("SMARTS charge is outside the supported range");
            text += ";" + std::string(charge < 0 ? "-" : "+") + std::to_string(charge < 0 ? -charge : charge);
        }
        if (chemistry.ringMatchesRingOnly) text += graph.ring[atom] ? ";R" : ";!R";
        return "[" + text + "]";
    };
    auto bondText = [&](int first, int second) {
        const int order = graph.bondOrder(first, second);
        const bool aromatic = graph.bondAromatic(first, second);
        if (aromatic != (order == 4))
            throw std::invalid_argument("SMARTS export needs normalised aromatic bond annotations");
        const bool strict = chemistry.aromaticityMode == ChemOptions::AromaticityMode::STRICT;
        const bool any = chemistry.matchBondOrder != ChemOptions::BondOrderMode::STRICT;
        std::string text;
        if (strict && aromatic) text = ":";
        else if (any) text = strict ? "!:" : "~";
        else {
            if (!strict && !aromatic && graph.aromatic[first] && graph.aromatic[second]
                && (order == 1 || order == 2))
                throw std::invalid_argument("SMARTS export cannot represent the aromatic resonance bond rule");
            if (aromatic) text = "-,=,:";
            else {
                if (order < 1 || order > 3) throw std::invalid_argument("Unsupported SMARTS bond order");
                text = order == 1 ? "-" : order == 2 ? "=" : "#";
                if (!strict && order <= 2) text += ",:";
                if (strict) text += ";!:";
            }
        }
        if (chemistry.ringMatchesRingOnly) text += graph.bondInRing(first, second) ? ";@" : ";!@";
        return text;
    };
    std::vector<std::vector<int>> adjacent = graph.neighbors;
    for (auto& neighbors : adjacent) std::sort(neighbors.begin(), neighbors.end());
    std::vector<int> parent(graph.n, -1), visit(graph.n, -1), roots;
    std::vector<std::vector<int>> children(graph.n);
    // Ring number of each closure bond while it is open; numbers are reused once closed.
    std::map<int64_t, int> closures;
    int visited = 0;
    struct TraversalFrame { int atom; size_t next; };
    std::vector<TraversalFrame> traversal;
    for (int atom = 0; atom < graph.n; ++atom) if (visit[atom] < 0) {
        roots.push_back(atom);
        visit[atom] = visited++;
        traversal.push_back({atom, 0});
        while (!traversal.empty()) {
            auto& frame = traversal.back();
            if (frame.next == adjacent[frame.atom].size()) {
                traversal.pop_back();
                continue;
            }
            const int current = frame.atom;
            const int neighbor = adjacent[current][frame.next++];
            if (neighbor == parent[current]) continue;
            if (visit[neighbor] < 0) {
                parent[neighbor] = current;
                children[current].push_back(neighbor);
                visit[neighbor] = visited++;
                traversal.push_back({neighbor, 0});
            } else if (visit[neighbor] < visit[current]) {
                closures.emplace(MolGraph::bondKey(current, neighbor), 0);
            }
        }
    }
    std::string result;
    std::vector<bool> openRings(100, false);
    auto appendAtom = [&](int atom) {
        result += atomText(atom);
        for (int neighbor : adjacent[atom]) {
            auto closure = closures.find(MolGraph::bondKey(atom, neighbor));
            if (closure == closures.end()) continue;
            if (visit[atom] < visit[neighbor]) {
                int ring = 1;
                while (ring < 100 && openRings[ring]) ++ring;
                if (ring == 100)
                    throw std::invalid_argument("SMARTS output needs more than 99 simultaneously open ring closures");
                openRings[ring] = true;
                closure->second = ring;
                result += bondText(atom, neighbor);
            } else {
                openRings[closure->second] = false;
            }
            const int ring = closure->second;
            result += (ring < 10 ? "" : "%") + std::to_string(ring);
        }
    };
    struct EmitFrame { int atom; size_t next; bool closeBranch; };
    std::vector<EmitFrame> emission;
    for (int root : roots) {
        if (!result.empty()) result += ".";
        appendAtom(root);
        emission.push_back({root, 0, false});
        while (!emission.empty()) {
            auto& frame = emission.back();
            if (frame.next == children[frame.atom].size()) {
                if (frame.closeBranch) result += ")";
                emission.pop_back();
                continue;
            }
            const int current = frame.atom;
            const int child = children[current][frame.next++];
            const bool branch = frame.next < children[current].size();
            if (branch) result += "(";
            result += bondText(current, child);
            appendAtom(child);
            emission.push_back({child, 0, branch});
        }
    }
    return result;
}
} // namespace smsd

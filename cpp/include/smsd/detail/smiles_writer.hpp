/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Canonical SMILES writer and MolGraph::toCanonicalSmiles(). Part of mol_graph.hpp,
 * which includes this header; including it first includes the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include "smsd/smiles_utils.hpp"

#include <algorithm>
#include <deque>
#include <functional>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace smsd {
namespace detail {

// --------------------------------------------------------------------------
// Canonical SMILES writer internals
// --------------------------------------------------------------------------

// Get the aromatic lowercase symbol for an element (returns "" if not applicable)
inline std::string aromaticSymbol(int z) {
    switch (z) {
        case 5:  return "b";
        case 6:  return "c";
        case 7:  return "n";
        case 8:  return "o";
        case 15: return "p";
        case 16: return "s";
        case 34: return "se";
        case 52: return "te";
        default: return "";
    }
}

// Determine if atom needs brackets
inline bool needsBrackets(const MolGraph& g, int idx) {
    int z = g.atomicNum[idx];
    if (!isOrganicSubset(z)) return true;
    if (g.formalCharge[idx] != 0) return true;
    if (g.massNumber[idx] != 0) return true;
    if (!g.atomClass.empty() && g.atomClass[idx] != 0) return true;
    if (g.tetraChirality[idx] != 0) return true;
    for (int neighbor : g.neighbors[idx])
        if (g.bondOrder(idx, neighbor) == 4 && g.dbStereo(idx, neighbor) != 0) return true;
    // If aromatic, check if it's one of the standard aromatic atoms
    if (g.aromatic[idx]) {
        std::string asym = aromaticSymbol(z);
        if (asym.empty()) return true;
    }
    return false;
}

/** Use supplied hydrogen counts, including zero, before the valence fallback. */
inline int hydrogenCountForWrite(const MolGraph& graph, int atom) {
    const int atomicNumber = graph.atomicNum[atom];
    const bool aromaticDonor = graph.aromatic[atom] && graph.formalCharge[atom] == 0
        && (atomicNumber == 7 || atomicNumber == 15);
    if (aromaticDonor && graph.neighbors[atom].size() >= 3) return 0;
    if (graph.hydrogenCountsKnown && atom < static_cast<int>(graph.hydrogenCount.size()))
        return graph.hydrogenCount[atom];
    if (atomicNumber <= 1) return 0;
    int hydrogens = computeImplicitH(atomicNumber, graph.aromatic[atom], explicitValence(graph, atom), graph.formalCharge[atom]);
    if (aromaticDonor && graph.neighbors[atom].size() == 2 && hydrogens == 0) {
        for (const auto& ring : graph.computeRings()) {
            if (ring.size() != 5 || std::find(ring.begin(), ring.end(), atom) == ring.end()) continue;
            int heteroatoms = 0;
            for (int member : ring)
                if (graph.aromatic[member] && graph.atomicNum[member] != 6) ++heteroatoms;
            if (heteroatoms == 1) return 1;
        }
    }
    return hydrogens;
}

// Write atom SMILES
inline std::string writeAtom(const MolGraph& g, int idx, int implicitH, int chirality = -1) {
    if (chirality < 0) chirality = g.tetraChirality[idx];
    bool brackets = needsBrackets(g, idx);
    if (!brackets) {
        const int defaultHydrogens = computeImplicitH(
            g.atomicNum[idx], g.aromatic[idx], explicitValence(g, idx), g.formalCharge[idx]);
        brackets = implicitH != defaultHydrogens
            || (implicitH > 0 && g.aromatic[idx] && (g.atomicNum[idx] == 7 || g.atomicNum[idx] == 15));
    }
    std::ostringstream oss;

    if (brackets) {
        oss << '[';
        if (g.massNumber[idx] != 0) oss << g.massNumber[idx];

        if (g.atomicNum[idx] == 0 && !g.atomClass.empty() && g.atomClass[idx] > 0
            && g.formalCharge[idx] == 0 && g.massNumber[idx] == 0
            && g.tetraChirality[idx] == 0 && implicitH == 0) {
            oss << 'R' << g.atomClass[idx];
        } else if (g.aromatic[idx]) {
            std::string asym = aromaticSymbol(g.atomicNum[idx]);
            if (!asym.empty())
                oss << asym;
            else
                oss << atomicNumToSymbol(g.atomicNum[idx]);
        } else {
            oss << atomicNumToSymbol(g.atomicNum[idx]);
        }

        if (chirality == 1) oss << '@';
        else if (chirality == 2) oss << "@@";

        if (implicitH > 0) {
            oss << 'H';
            if (implicitH > 1) oss << implicitH;
        }

        if (g.formalCharge[idx] > 0) {
            oss << '+';
            if (g.formalCharge[idx] > 1) oss << g.formalCharge[idx];
        } else if (g.formalCharge[idx] < 0) {
            oss << '-';
            if (g.formalCharge[idx] < -1) oss << -g.formalCharge[idx];
        }

        if (!g.atomClass.empty() && g.atomClass[idx] != 0 && g.atomicNum[idx] != 0) {
            oss << ':' << g.atomClass[idx];
        } else if (!g.atomClass.empty() && g.atomClass[idx] != 0 && g.atomicNum[idx] == 0
                   && !(g.formalCharge[idx] == 0 && g.massNumber[idx] == 0
                        && g.tetraChirality[idx] == 0 && implicitH == 0)) {
            oss << ':' << g.atomClass[idx];
        }

        oss << ']';
    } else {
        if (g.aromatic[idx]) {
            oss << aromaticSymbol(g.atomicNum[idx]);
        } else {
            oss << atomicNumToSymbol(g.atomicNum[idx]);
        }
    }

    return oss.str();
}

// Bond order symbol
inline std::string bondSymbol(int order, bool fromAromatic, bool toAromatic) {
    if (order == 4) return ""; // aromatic implicit
    if (order == 1) {
        if (fromAromatic && toAromatic) return "-"; // explicit single between aromatics
        return ""; // implicit single
    }
    if (order == 2) return "=";
    if (order == 3) return "#";
    return "";
}

// --------------------------------------------------------------------------
// Canonical depth-first traversal shared by the SMILES and SMARTS writers
// --------------------------------------------------------------------------

/** A ring-closure digit at an atom: its number, the bond order and the partner atom. */
struct RingDigit { int rn; int bo; int partner; };

/** Connected components of the graph, each in breadth-first discovery order. */
inline std::vector<std::vector<int>> connectedComponents(const MolGraph& g) {
    std::vector<bool> seen(g.n, false);
    std::vector<std::vector<int>> components;
    for (int i = 0; i < g.n; i++) {
        if (seen[i]) continue;
        std::vector<int> component;
        std::deque<int> queue;
        queue.push_back(i);
        seen[i] = true;
        while (!queue.empty()) {
            int u = queue.front(); queue.pop_front();
            component.push_back(u);
            for (int v : g.neighbors[u])
                if (!seen[v]) { seen[v] = true; queue.push_back(v); }
        }
        components.push_back(std::move(component));
    }
    return components;
}

/**
 * Writes one connected component depth first from its lowest canonical label, visiting
 * neighbours in canonical order. A first pass assigns ring-closure numbers to the back
 * edges; the second pass emits, for every atom, atomToken(atom, parent, children, ringDigits),
 * then each ring digit preceded by ringBond(atom, digit, opening), then the children,
 * each preceded by bond(atom, child) and all but the last wrapped in a branch.
 */
template <class AtomToken, class RingBond, class Bond>
inline std::string writeCanonicalComponent(const MolGraph& g, const std::vector<int>& canon,
                                           const std::vector<int>& component,
                                           AtomToken&& atomToken, RingBond&& ringBond, Bond&& bond) {
    const int n = g.n;
    auto sortedNeighbors = [&](int atom) {
        std::vector<int> nbs = g.neighbors[atom];
        std::sort(nbs.begin(), nbs.end(), [&](int a, int b) { return canon[a] < canon[b]; });
        return nbs;
    };

    int root = component[0];
    for (int a : component)
        if (canon[a] < canon[root]) root = a;

    // Pass 1: the depth-first tree, and one ring-closure number for each back edge,
    // recorded at both of its atoms.
    std::vector<bool> vis(n, false);
    std::set<int64_t> treeEdges;
    std::set<int64_t> assignedRings;
    std::vector<std::vector<RingDigit>> ringDigits(n);
    int ringCounter = 1;
    std::function<void(int, int)> dfs1 = [&](int atom, int parent) {
        vis[atom] = true;
        bool parentSkipped = false;
        for (int nb : sortedNeighbors(atom)) {
            if (nb == parent && !parentSkipped) { parentSkipped = true; continue; }
            if (vis[nb]) {
                if (assignedRings.insert(MolGraph::bondKey(atom, nb)).second) {
                    const int rn = ringCounter++;
                    const int bo = g.bondOrder(atom, nb);
                    ringDigits[atom].push_back({rn, bo, nb});
                    ringDigits[nb].push_back({rn, bo, atom});
                }
            } else {
                treeEdges.insert(MolGraph::bondKey(atom, nb));
                dfs1(nb, atom);
            }
        }
    };
    dfs1(root, -1);

    // Pass 2: emit the tokens.
    vis.assign(n, false);
    std::string text;
    std::unordered_set<int> emittedRings;
    std::function<void(int, int)> dfs2 = [&](int atom, int parent) {
        vis[atom] = true;
        std::vector<int> children;
        for (int nb : sortedNeighbors(atom))
            if (!vis[nb] && treeEdges.count(MolGraph::bondKey(atom, nb))) children.push_back(nb);
        text += atomToken(atom, parent, children, ringDigits[atom]);
        for (const auto& digit : ringDigits[atom]) {
            const bool opening = emittedRings.insert(digit.rn).second;
            text += ringBond(atom, digit, opening);
            if (digit.rn < 10) text += std::to_string(digit.rn);
            else text += "%" + std::to_string(digit.rn);
        }
        for (size_t i = 0; i < children.size(); i++) {
            const int nb = children[i];
            if (vis[nb]) continue;
            const std::string symbol = bond(atom, nb);
            if (i < children.size() - 1) {
                text += "(" + symbol;
                dfs2(nb, atom);
                text += ")";
            } else {
                text += symbol;
                dfs2(nb, atom);
            }
        }
    };
    dfs2(root, -1);
    return text;
}

// Canonical SMILES writer using DFS with canonical ordering.
//
// Strategy:
//   1. Run a preliminary DFS to identify all back-edges (ring closures)
//      and assign ring numbers.  For each ring bond, record the ring
//      number at both endpoint atoms.
//   2. Run a second DFS that emits atom symbols, ring digits, bonds,
//      and branches.
inline std::string writeCanonicalSMILES(const MolGraph& g) {
    if (g.n == 0) return "";

    g.ensureCanonical();   // canonicalLabel is lazy
    int n = g.n;
    const auto& canon = g.canonicalLabel;

    std::vector<int> implicitH(n, 0);
    for (int atom = 0; atom < n; ++atom)
        implicitH[atom] = hydrogenCountForWrite(g, atom);

    // Directional bonds form a small parity-constraint graph. A shared
    // single or aromatic carrier must keep one direction.
    struct DirectionConstraint { int first, second, relation; };
    std::vector<DirectionConstraint> constraints;
    std::unordered_map<int64_t, int> directionIndex;
    auto directionBond = [&](int first, int second) {
        const int64_t key = MolGraph::bondKey(first, second);
        const auto inserted = directionIndex.emplace(key, static_cast<int>(directionIndex.size()));
        return inserted.first->second;
    };
    auto isDirectionalCarrier = [&](int first, int second) {
        const int order = g.bondOrder(first, second);
        return order == 1 || (order == 4 && g.dbStereo(first, second) == 0);
    };
    for (int atom = 0; atom < n; ++atom) {
        for (int other : g.neighbors[atom]) {
            const int order = g.bondOrder(atom, other);
            if (other <= atom || (order != 2 && order != 4) || g.dbStereo(atom, other) == 0) continue;
            int firstLigand = -1, secondLigand = -1;
            for (int neighbor : g.neighbors[atom])
                if (neighbor != other && (firstLigand < 0 || canon[neighbor] < canon[firstLigand])) firstLigand = neighbor;
            for (int neighbor : g.neighbors[other])
                if (neighbor != atom && (secondLigand < 0 || canon[neighbor] < canon[secondLigand])) secondLigand = neighbor;
            if (firstLigand < 0 || secondLigand < 0 ||
                !isDirectionalCarrier(atom, firstLigand) || !isDirectionalCarrier(other, secondLigand))
                throw std::invalid_argument("Double-bond stereo requires a single or aromatic reference bond at each end");
            const int firstOrientation = atom < firstLigand ? 1 : -1;
            const int secondOrientation = other < secondLigand ? 1 : -1;
            const int configuration = MolGraph::canonicalDoubleBondTag(g, atom, other, canon);
            const int relation = (configuration == 1 ? 1 : -1)
                                 * firstOrientation * secondOrientation;
            constraints.push_back({directionBond(atom, firstLigand), directionBond(other, secondLigand), relation});
        }
    }
    for (int atom = 0; atom < n; ++atom) {
        std::vector<int> markedLigands;
        bool doubleBondEnd = false;
        for (int neighbor : g.neighbors[atom]) {
            const int order = g.bondOrder(atom, neighbor);
            if (order == 2 || (order == 4 && g.dbStereo(atom, neighbor) != 0)) doubleBondEnd = true;
            if (directionIndex.count(MolGraph::bondKey(atom, neighbor))) markedLigands.push_back(neighbor);
        }
        if (!doubleBondEnd || markedLigands.size() != 2) continue;
        const int first = markedLigands[0], second = markedLigands[1];
        const int relation = -(atom < first ? 1 : -1) * (atom < second ? 1 : -1);
        constraints.push_back({directionBond(atom, first), directionBond(atom, second), relation});
    }
    std::vector<std::vector<std::pair<int, int>>> directionEdges(directionIndex.size());
    for (const auto& constraint : constraints) {
        directionEdges[constraint.first].push_back({constraint.second, constraint.relation});
        directionEdges[constraint.second].push_back({constraint.first, constraint.relation});
    }
    std::vector<int> directions(directionIndex.size(), 0);
    std::vector<int64_t> directionBonds(directionIndex.size());
    for (const auto& bond : directionIndex) directionBonds[bond.second] = bond.first;
    auto canonicalEnds = [&](int64_t bond) {
        const int first = canon[static_cast<int>(bond >> 32)], second = canon[static_cast<int>(bond)];
        return std::make_pair(std::min(first, second), std::max(first, second));
    };
    std::vector<int> directionOrder(directionIndex.size());
    std::iota(directionOrder.begin(), directionOrder.end(), 0);
    std::sort(directionOrder.begin(), directionOrder.end(), [&](int first, int second) {
        return canonicalEnds(directionBonds[first]) < canonicalEnds(directionBonds[second]);
    });
    for (int start : directionOrder) {
        if (directions[start] != 0) continue;
        const int64_t bond = directionBonds[start];
        directions[start] = canon[static_cast<int>(bond >> 32)] < canon[static_cast<int>(bond)] ? 1 : -1;
        std::deque<int> pending{static_cast<int>(start)};
        while (!pending.empty()) {
            const int current = pending.front(); pending.pop_front();
            for (const auto& edge : directionEdges[current]) {
                const int required = directions[current] * edge.second;
                if (directions[edge.first] == 0) {
                    directions[edge.first] = required;
                    pending.push_back(edge.first);
                } else if (directions[edge.first] != required) {
                    throw std::invalid_argument("Inconsistent double-bond stereo annotations");
                }
            }
        }
    }
    auto writeBond = [&](int first, int second) {
        if (g.bondOrder(first, second) == 4 && g.dbStereo(first, second) != 0) return std::string("=");
        const auto found = directionIndex.find(MolGraph::bondKey(first, second));
        if (found != directionIndex.end()) {
            const int direction = directions[found->second] * (first < second ? 1 : -1);
            return std::string(1, direction == 1 ? '/' : '\\');
        }
        return bondSymbol(g.bondOrder(first, second), g.aromatic[first], g.aromatic[second]);
    };

    // Each component is written from its lowest canonical label; the SMILES of a
    // tetrahedral centre carries the parity of the emitted ligand order against the
    // stored one.
    std::vector<std::string> componentSmiles;
    for (const auto& component : connectedComponents(g)) {
        componentSmiles.push_back(writeCanonicalComponent(g, canon, component,
            [&](int atom, int parent, const std::vector<int>& children, const std::vector<RingDigit>& digits) {
                int chirality = g.tetraChirality[atom];
                if (chirality != 0) {
                    std::vector<int> storedOrder = g.neighbors[atom];
                    std::vector<int> emittedOrder;
                    if (parent >= 0) emittedOrder.push_back(parent);
                    if (storedOrder.size() == 3) {
                        const size_t hydrogenPosition = atom > 0 && storedOrder[0] < atom ? 1 : 0;
                        storedOrder.insert(storedOrder.begin() + hydrogenPosition, -1);
                        emittedOrder.push_back(-1);
                    }
                    for (const auto& ring : digits) emittedOrder.push_back(ring.partner);
                    emittedOrder.insert(emittedOrder.end(), children.begin(), children.end());
                    int inversions = 0;
                    if (storedOrder.size() == emittedOrder.size()) {
                        for (size_t first = 0; first < emittedOrder.size(); ++first) {
                            const auto firstPosition = std::find(storedOrder.begin(), storedOrder.end(), emittedOrder[first]);
                            for (size_t second = first + 1; second < emittedOrder.size(); ++second) {
                                const auto secondPosition = std::find(storedOrder.begin(), storedOrder.end(), emittedOrder[second]);
                                if (firstPosition > secondPosition) ++inversions;
                            }
                        }
                    }
                    if (inversions % 2 != 0) chirality = 3 - chirality;
                }
                return writeAtom(g, atom, implicitH[atom], chirality);
            },
            [&](int atom, const RingDigit& digit, bool opening) {
                return opening ? writeBond(atom, digit.partner) : std::string();
            },
            [&](int atom, int nb) { return writeBond(atom, nb); }));
    }

    std::sort(componentSmiles.begin(), componentSmiles.end());
    std::ostringstream result;
    for (size_t component = 0; component < componentSmiles.size(); ++component) {
        if (component > 0) result << '.';
        result << componentSmiles[component];
    }
    return result.str();
}

} // namespace detail

inline std::string MolGraph::toCanonicalSmiles() const {
    return detail::writeCanonicalSMILES(*this);
}

} // namespace smsd

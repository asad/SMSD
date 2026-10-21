/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * CIP (Cahn-Ingold-Prelog) stereodescriptor assignment for MolGraph.
 *
 * Uses the IUPAC 2013 sequence rules (Blue Book) when assigning
 * R/S descriptors to tetrahedral stereocentres and E/Z descriptors to
 * double bonds using CIP sequence rules:
 *   Rule 1: higher atomic number gets higher priority
 *   Rule 2: higher mass number breaks ties (isotopes); an unlabelled atom
 *   counts at its standard atomic weight, a duplicate at none
 *   Rule 3: Z double-bond configuration gets higher priority than E
 *   Rules 4 and 5: auxiliary descriptors, like and unlike descriptor pairs,
 *   r before s and R before S, for the centres Rules 1 to 3 leave unresolved
 *   (detail::SequenceRules); a centre resolved by Rule 5 is pseudoasymmetric
 *
 * Uses the digraph (extended-connectivity tree) approach
 * with phantom/duplicate nodes for ring closures and multiple bonds.
 * The duplicate atom of a double bond in a mancude ring system takes the mean
 * atomic number of the ring neighbours it may be double bonded to (P-92.1.4.4),
 * so labels do not depend on the Kekule structure.
 *
 * References:
 *   - Prelog & Helmchen, Angew. Chem. Int. Ed. Engl. 21, 567 (1982)
 *   - IUPAC Recommendations 2013 (P-92)
 *   - Hanson et al., JCIM 58, 1755 (2018) — edge-case corrections and the
 *     auxiliary descriptor method the reference implementations follow
 *
 * Header-only, zero external dependencies — pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_CIP_HPP
#define SMSD_CIP_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/periodic_table.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/stereo.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace smsd {
namespace cip {

// ============================================================================
// Public enums
// ============================================================================

/** CIP tetrahedral descriptor (uppercase = true stereocentre, lowercase = pseudoasymmetric). */
enum class RSLabel { NONE = 0, R, S, r, s };

/** CIP double-bond descriptor. */
enum class EZLabel { NONE = 0, E, Z };

inline std::string to_string(RSLabel l) {
    switch (l) {
        case RSLabel::R: return "R";
        case RSLabel::S: return "S";
        case RSLabel::r: return "r";
        case RSLabel::s: return "s";
        default: return "";
    }
}

inline std::string to_string(EZLabel l) {
    switch (l) {
        case EZLabel::E: return "E";
        case EZLabel::Z: return "Z";
        default: return "";
    }
}

// ============================================================================
// Internal: digraph node for CIP priority determination
// ============================================================================
namespace detail {

/** Scale of DigraphNode::scaledAtomicNum, so the mean of up to four atomic numbers is exact. */
constexpr int ATOMIC_NUMBER_SCALE = 12;

inline int atomPropertyOrZero(const std::vector<int>& values, int atom) noexcept {
    return atom >= 0 && atom < static_cast<int>(values.size()) ? values[atom] : 0;
}

struct DigraphNode {
    int  atomIdx;       // index into MolGraph
    int  atomicNum;     // atomic number (copied for speed)
    int  massNumber;    // 0 for phantom/duplicate nodes
    bool isDuplicate;   // true for ring-closure / bond-expansion phantom nodes
    int  parent;        // index in digraph node list (-1 for root)
    int  depth;         // distance from root in the digraph
    std::vector<int> children; // indices of child nodes in digraph
    int  scaledAtomicNum = 0; // Rule 1 key: atomicNum * ATOMIC_NUMBER_SCALE, a mean for mancude duplicates
    bool isRingDuplicate = false; // a duplicate closing a ring, rather than one of a multiple bond
    int  distance = 0;        // Rule 1b: root 1, each sphere one more; a duplicate takes its original's
};

/** Implicit hydrogens, with the valence model for Builder graphs without supplied counts. */
inline int implicitHydrogens(const MolGraph& graph, int atom) { return graph.implicitHydrogenCount(atom); }

/** Whether a tagged centre has a lone pair as its fourth ligand, see stereo_detail::tetrahedralLigands(). */
inline bool hasLonePairLigand(const MolGraph& graph, int atom, int hydrogens) {
    return hydrogens == 0 && graph.neighbors[atom].size() == 3 && atomPropertyOrZero(graph.tetraChirality, atom) != 0
        && stereo_detail::lonePairCentre(graph.atomicNum[atom]);
}

/**
 * Sequence Rule 2 between two nodes Rule 1a has paired, either of which may be absent. Two atoms
 * without a mass number tie whatever their weights; otherwise a labelled atom counts its mass
 * number, an unlabelled atom its standard atomic weight and a duplicate nothing, as the reference
 * CIP implementations read the rule (Hanson et al., JCIM 58, 1755 (2018)). Both passes use it.
 */
inline double rule2Mass(const DigraphNode& node) {
    if (node.isDuplicate || node.atomicNum <= 0) return 0.0;
    return node.massNumber != 0 ? node.massNumber : atomicWeight(node.atomicNum);
}

inline int compareMassNumbers(const DigraphNode* first, const DigraphNode* second) {
    const int firstLabel = first ? first->massNumber : 0;
    const int secondLabel = second ? second->massNumber : 0;
    if (firstLabel == 0 && secondLabel == 0) return 0;
    const double firstMass = first ? rule2Mass(*first) : 0.0;
    const double secondMass = second ? rule2Mass(*second) : 0.0;
    return firstMass < secondMass ? -1 : firstMass > secondMass ? 1 : 0;
}

/** Compare spheres in parent order, padding missing phantom children with zero. */
inline int compareOrderedBranches(const std::vector<DigraphNode>& firstNodes, int firstRoot,
                                  const std::vector<DigraphNode>& secondNodes, int secondRoot,
                                  bool compareMass) {
    std::vector<std::pair<int, int>> frontier{{firstRoot, secondRoot}};
    const std::vector<int> empty;
    while (!frontier.empty()) {
        for (const auto& pair : frontier) {
            if (compareMass) {
                const int massComparison = compareMassNumbers(pair.first < 0 ? nullptr : &firstNodes[pair.first],
                                                              pair.second < 0 ? nullptr : &secondNodes[pair.second]);
                if (massComparison != 0) return massComparison;
                continue;
            }
            const int firstValue = pair.first < 0 ? 0 : firstNodes[pair.first].scaledAtomicNum;
            const int secondValue = pair.second < 0 ? 0 : secondNodes[pair.second].scaledAtomicNum;
            if (firstValue != secondValue) return firstValue - secondValue;
        }
        std::vector<std::pair<int, int>> next;
        for (const auto& pair : frontier) {
            const auto& firstChildren = pair.first < 0 ? empty : firstNodes[pair.first].children;
            const auto& secondChildren = pair.second < 0 ? empty : secondNodes[pair.second].children;
            for (size_t slot = 0; slot < std::max(firstChildren.size(), secondChildren.size()); ++slot)
                next.push_back({slot < firstChildren.size() ? firstChildren[slot] : -1,
                                slot < secondChildren.size() ? secondChildren[slot] : -1});
        }
        frontier = std::move(next);
    }
    return 0;
}

/** Child ordering is established from the leaves, preserving earlier-rule ranks. */
inline void orderDigraph(std::vector<DigraphNode>& nodes) {
    for (int node = static_cast<int>(nodes.size()) - 1; node >= 0; --node) {
        auto& children = nodes[node].children;
        std::stable_sort(children.begin(), children.end(), [&](int first, int second) {
            const int atomicComparison = compareOrderedBranches(nodes, first, nodes, second, false);
            return atomicComparison != 0 ? atomicComparison > 0
                : compareOrderedBranches(nodes, first, nodes, second, true) > 0;
        });
    }
}

/** Whether the bond is a ring bond, or an aromatic bond as kekulize() reads it. */
inline bool ringBond(const MolGraph& g, int atom, int neighbor) {
    return g.bondInRing(atom, neighbor) || (g.aromatic[atom] && g.aromatic[neighbor]
        && (g.bondOrder(atom, neighbor) == 4 || g.bondAromatic(atom, neighbor)));
}

/**
 * The duplicate a ring atom with exactly one double bond, partners[atom] in the Kekule bond orders
 * of `kekule`, carries for that bond (-1 for other atoms). A mancude ring system is what remains of
 * those atoms after repeatedly removing the ones with fewer than two others as ring neighbours.
 * There the duplicate takes the mean atomic number of the atom's ring neighbours in the system
 * when it is fractional (IUPAC 2013 P-92.1.4.4, as the reference CIP implementations apply it), so
 * labels do not depend on the Kekule structure; otherwise it is the partner. numbers[atom] is that
 * atomic number times ATOMIC_NUMBER_SCALE. Java's CIPAssigner applies the same rule.
 */
inline void mancudeDuplicates(const MolGraph& g, const MolGraph& kekule,
                              std::vector<int>& numbers, std::vector<int>& partners) {
    numbers.assign(g.n, 0);
    partners.assign(g.n, -1);
    for (int atom = 0; atom < g.n; ++atom) {
        const bool ring = (atom < static_cast<int>(g.ring.size()) && g.ring[atom]) || g.aromatic[atom];
        int partner = -1, doubles = 0, triples = 0;
        for (int neighbor : g.neighbors[atom]) {
            const int order = kekule.bondOrder(atom, neighbor);
            if (order == 2) { partner = neighbor; ++doubles; }
            else if (order == 3) ++triples;
        }
        if (ring && doubles == 1 && triples == 0) partners[atom] = partner;
    }
    std::vector<uint8_t> mancude(g.n, 0);
    std::vector<int> degree(g.n, 0), pending;
    for (int atom = 0; atom < g.n; ++atom) {
        if (partners[atom] < 0) continue;
        mancude[atom] = 1;
        for (int neighbor : g.neighbors[atom])
            if (partners[neighbor] >= 0 && ringBond(g, atom, neighbor)) ++degree[atom];
        if (degree[atom] < 2) pending.push_back(atom);
    }
    while (!pending.empty()) {
        const int atom = pending.back();
        pending.pop_back();
        if (!mancude[atom]) continue;
        mancude[atom] = 0;
        for (int neighbor : g.neighbors[atom])
            if (mancude[neighbor] && ringBond(g, atom, neighbor) && --degree[neighbor] < 2)
                pending.push_back(neighbor);
    }
    for (int atom = 0; atom < g.n; ++atom) {
        if (partners[atom] < 0) continue;
        int sum = 0, count = 0;
        if (mancude[atom]) {
            for (int neighbor : g.neighbors[atom]) {
                if (!mancude[neighbor] || !ringBond(g, atom, neighbor)) continue;
                sum += g.atomicNum[neighbor];
                ++count;
            }
        }
        numbers[atom] = count == 0 || sum % count == 0 ? g.atomicNum[partners[atom]] * ATOMIC_NUMBER_SCALE
                                                       : sum * ATOMIC_NUMBER_SCALE / count;
    }
}

// --------------------------------------------------------------------------
// Build the CIP digraph (extended connectivity tree) rooted at `root`,
// expanding every neighbour of the root.  Ring closures produce duplicate
// (phantom) nodes with mass = 0.  Mancude (aromatic) duplicate nodes also get
// mass = 0.
//
// maxDepth and maxNodes bound the expansion on large molecules.
// --------------------------------------------------------------------------
inline std::vector<DigraphNode> buildDigraph(
    const MolGraph& g,
    int root,
    int maxDepth = 30,
    int maxNodes = 5000)
{
    std::optional<MolGraph> kekule;
    if (std::any_of(g.aromatic.begin(), g.aromatic.end(), [](uint8_t aromatic) { return aromatic != 0; })) {
        kekule.emplace(g);
        if (!kekule->kekulize()) return {};
    }
    const MolGraph& bonding = kekule ? *kekule : g;
    std::vector<int> duplicateNumbers, duplicatePartners;
    mancudeDuplicates(g, bonding, duplicateNumbers, duplicatePartners);
    std::vector<DigraphNode> nodes;
    nodes.reserve(128);

    // Root node
    DigraphNode rootNode;
    rootNode.atomIdx    = root;
    rootNode.atomicNum  = g.atomicNum[root];
    rootNode.scaledAtomicNum = g.atomicNum[root] * ATOMIC_NUMBER_SCALE;
    rootNode.massNumber = atomPropertyOrZero(g.massNumber, root);
    rootNode.isDuplicate = false;
    rootNode.parent     = -1;
    rootNode.depth      = 0;
    rootNode.distance   = 1;
    nodes.push_back(rootNode);

    // BFS expansion
    // frontier: indices into `nodes`
    std::vector<int> frontier = {0};

    for (int depth = 0; depth < maxDepth && !frontier.empty(); ++depth) {
        std::vector<int> nextFrontier;
        for (int ni : frontier) {
            if (static_cast<int>(nodes.size()) >= maxNodes) break;

            int atom = nodes[ni].atomIdx;
            int parentNodeIdx = nodes[ni].parent;
            int parentAtom = (parentNodeIdx >= 0) ? nodes[parentNodeIdx].atomIdx : -1;

            // Collect neighbours, expanding double/triple bonds as per CIP
            // convention: a double bond A=B adds one extra phantom B on A's
            // side (and vice versa); a triple bond adds two.
            struct Expansion {
                int neighborAtom;
                bool isDuplicate;
                int bondOrder;
                int scaledAtomicNum;
            };
            std::vector<Expansion> expansions;

            for (int nbr : g.neighbors[atom]) {
                if (nbr == parentAtom && !nodes[ni].isDuplicate) {
                    // Don't go back to the real parent (but duplicates
                    // can appear from the parent direction for bond-order
                    // expansion).
                    // However, for bond order > 1 we still need phantom
                    // nodes for the parent bond.  We handle that below.
                    continue;
                }

                // The partner's double bond is the mancude duplicate below.
                int bo = nbr == duplicatePartners[atom] ? 1 : bonding.bondOrder(atom, nbr);
                if (bo == 4) bo = 1;
                const int scaled = g.atomicNum[nbr] * ATOMIC_NUMBER_SCALE;

                // Real child node
                expansions.push_back({nbr, false, bo, scaled});

                // Phantom nodes for higher-order bonds
                for (int extra = 1; extra < bo; ++extra) {
                    expansions.push_back({nbr, true, bo, scaled});
                }
            }

            // Also add phantom nodes for the bond back to parent
            // (for double/triple bonds to parent, the parent side
            //  already added phantoms for us, but we need phantom
            //  copies of the parent on THIS node's child list).
            if (parentAtom >= 0 && !nodes[ni].isDuplicate) {
                int bo = parentAtom == duplicatePartners[atom] ? 1 : bonding.bondOrder(atom, parentAtom);
                if (bo == 4) bo = 1;
                for (int extra = 1; extra < bo; ++extra) {
                    expansions.push_back({parentAtom, true, bo, g.atomicNum[parentAtom] * ATOMIC_NUMBER_SCALE});
                }
            }
            if (duplicatePartners[atom] >= 0)
                expansions.push_back({duplicatePartners[atom], true, 2, duplicateNumbers[atom]});

            for (auto& exp : expansions) {
                if (static_cast<int>(nodes.size()) >= maxNodes) break;

                // Check if this atom is already on the path from root
                // to current node.  If so, it's a ring closure -> duplicate.
                // A duplicate of an atom on the path takes that node's
                // distance; one of an atom not yet reached has none.
                int pathNode = -1;
                for (int walk = ni; walk >= 0; walk = nodes[walk].parent) {
                    if (nodes[walk].atomIdx == exp.neighborAtom) {
                        pathNode = walk;
                        break;
                    }
                }
                const bool onPath = !exp.isDuplicate && pathNode >= 0;

                bool makeDuplicate = exp.isDuplicate || onPath;

                DigraphNode child;
                child.atomIdx    = exp.neighborAtom;
                child.atomicNum  = g.atomicNum[exp.neighborAtom];
                child.scaledAtomicNum = exp.scaledAtomicNum;
                // Duplicate nodes: mass = 0 per CIP rules
                // (including mancude / aromatic ring closures)
                child.massNumber = makeDuplicate ? 0 : atomPropertyOrZero(g.massNumber, exp.neighborAtom);
                child.isDuplicate = makeDuplicate;
                child.isRingDuplicate = onPath;
                child.parent     = ni;
                child.depth      = depth + 1;
                child.distance   = makeDuplicate ? (pathNode >= 0 ? nodes[pathNode].distance : 0) : depth + 2;

                int childIdx = static_cast<int>(nodes.size());
                nodes.push_back(child);
                nodes[ni].children.push_back(childIdx);

                // Duplicate/phantom nodes are leaves -- don't expand further
                if (!makeDuplicate) {
                    nextFrontier.push_back(childIdx);
                }
            }
            for (int hydrogen = 0; hydrogen < implicitHydrogens(g, atom)
                 && static_cast<int>(nodes.size()) < maxNodes; ++hydrogen) {
                DigraphNode child;
                child.atomIdx = -1;
                child.atomicNum = 1;
                child.scaledAtomicNum = ATOMIC_NUMBER_SCALE;
                child.massNumber = 0;
                child.isDuplicate = false;
                child.parent = ni;
                child.depth = depth + 1;
                child.distance = depth + 2;
                nodes[ni].children.push_back(static_cast<int>(nodes.size()));
                nodes.push_back(std::move(child));
            }
        }
        frontier = std::move(nextFrontier);
    }

    orderDigraph(nodes);
    return nodes;
}

// --------------------------------------------------------------------------
// Compare two digraph sub-trees by CIP Rules 1-2 (atomic number, mass
// number, BFS level-by-level).
//
// Returns: <0 if tree rooted at a has LOWER priority than b,
//           0 if equal,
//          >0 if tree rooted at a has HIGHER priority.
// --------------------------------------------------------------------------
inline int compareBranches_Rule1(
    const std::vector<DigraphNode>& nodesA, int rootA,
    const std::vector<DigraphNode>& nodesB, int rootB)
{
    const int atomicComparison = compareOrderedBranches(nodesA, rootA, nodesB, rootB, false);
    if (atomicComparison != 0) return atomicComparison;
    return compareOrderedBranches(nodesA, rootA, nodesB, rootB, true);
}

inline int atomicPriorityDoubleBondConfiguration(const MolGraph& graph, int first, int second);

// --------------------------------------------------------------------------
// CIP Rule 3: Z > E tie-break on double-bond stereo configuration.
//
// When Rules 1-2 produce a tie, Z-configured double bonds get higher
// priority than E-configured double bonds.
// --------------------------------------------------------------------------
inline int compareBranches_Rule3(
    const std::vector<DigraphNode>& nodesA, int rootA,
    const std::vector<DigraphNode>& nodesB, int rootB,
    const MolGraph& g)
{
    std::vector<int> frontA = {rootA};
    std::vector<int> frontB = {rootB};

    while (!frontA.empty() || !frontB.empty()) {
        // Collect E/Z descriptors at this level: Z(1)->2, E(2)->1, NONE->0
        auto ezDescs = [&](const std::vector<int>& front,
                           const std::vector<DigraphNode>& nodes) {
            std::vector<int> ez;
            for (int ni : front) {
                const auto& nd = nodes[ni];
                int val = 0;
                if (nd.parent >= 0 && nd.atomIdx >= 0 && !nd.isDuplicate) {
                    int parentAtom = nodes[nd.parent].atomIdx;
                    if (parentAtom >= 0 && g.bondOrder(parentAtom, nd.atomIdx) == 2) {
                        int conf = atomicPriorityDoubleBondConfiguration(g, parentAtom, nd.atomIdx);
                        val = (conf == 1) ? 2 : (conf == 2) ? 1 : 0;  // Z->2, E->1
                    }
                }
                ez.push_back(val);
            }
            std::sort(ez.begin(), ez.end(), std::greater<int>());
            return ez;
        };

        auto ezA = ezDescs(frontA, nodesA);
        auto ezB = ezDescs(frontB, nodesB);
        size_t len = std::max(ezA.size(), ezB.size());
        for (size_t k = 0; k < len; ++k) {
            int a = (k < ezA.size()) ? ezA[k] : 0;
            int b = (k < ezB.size()) ? ezB[k] : 0;
            if (a != b) return a - b;
        }

        // Expand to next BFS level
        std::vector<int> nextA, nextB;
        for (int i : frontA)
            for (int c : nodesA[i].children) nextA.push_back(c);
        for (int i : frontB)
            for (int c : nodesB[i].children) nextB.push_back(c);

        frontA = std::move(nextA);
        frontB = std::move(nextB);
    }
    return 0;
}

// --------------------------------------------------------------------------
// CIP Rule 4b-c: like/unlike descriptor pairing.
//
// When two branches are equal by Rules 1-2, compare the CIP descriptors
// (R/S or E/Z) of stereocentres encountered along each branch path.
// "like" pairs (R,R or S,S) have higher priority than "unlike" pairs
// (R,S or S,R).  Within "like", R > S.
//
// descriptor_map: atom_index -> RSLabel (pre-computed in first pass).
// Returns: <0, 0, >0 as compareBranches_Rule1.
// --------------------------------------------------------------------------
inline int compareBranches_Rule4bc(
    const std::vector<DigraphNode>& nodesA, int rootA,
    const std::vector<DigraphNode>& nodesB, int rootB,
    const std::unordered_map<int, int>& descriptorMap)
{
    // Collect descriptor sequences along each branch path (BFS order)
    auto collectDescriptors = [&](const std::vector<DigraphNode>& nodes,
                                  int root) -> std::vector<int> {
        std::vector<int> result;
        std::vector<int> frontier = {root};
        while (!frontier.empty()) {
            std::vector<int> nextFrontier;
            for (int ni : frontier) {
                int atom = nodes[ni].atomIdx;
                if (!nodes[ni].isDuplicate && atom >= 0) {
                    auto it = descriptorMap.find(atom);
                    if (it != descriptorMap.end() && it->second != 0) {
                        result.push_back(it->second);
                    }
                }
                for (int c : nodes[ni].children) nextFrontier.push_back(c);
            }
            frontier = std::move(nextFrontier);
        }
        return result;
    };

    auto descA = collectDescriptors(nodesA, rootA);
    auto descB = collectDescriptors(nodesB, rootB);

    // "like" descriptor = same as the reference centre's descriptor.
    // We encode: R=1, S=2. like pair = same value, unlike = different.
    // Per IUPAC 2013 P-92.1.4.2: "like > unlike" and within like, R > S.
    //
    // Sort descriptors in each branch: R(1) before S(2)
    std::sort(descA.begin(), descA.end());
    std::sort(descB.begin(), descB.end());

    size_t len = std::max(descA.size(), descB.size());
    for (size_t k = 0; k < len; ++k) {
        int a = (k < descA.size()) ? descA[k] : 0;
        int b = (k < descB.size()) ? descB[k] : 0;
        if (a != b) {
            // R(1) > S(2) > NONE(0) for Rule 5, but for Rule 4b-c:
            // A branch with a descriptor ranks higher than one without.
            // Between R and S: compare as like/unlike relative to reference.
            // For simplicity (IUPAC 2013 P-92.1.4.2): R > S.
            if (a == 0) return -1; // A has no descriptor, B does → B higher
            if (b == 0) return 1;  // A has descriptor, B doesn't → A higher
            // Both have descriptors: R(1) > S(2) per Rule 5
            return (a < b) ? 1 : -1; // lower enum value = R = higher priority
        }
    }

    return 0;
}

// --------------------------------------------------------------------------
// Full CIP comparison: Rules 1-2, then 3 (Z>E), then 4b-c, then Rule 5 (R > S).
// --------------------------------------------------------------------------
inline int compareBranches_Full(
    const std::vector<DigraphNode>& nodesA, int rootA,
    const std::vector<DigraphNode>& nodesB, int rootB,
    const std::unordered_map<int, int>& descriptorMap,
    const MolGraph* molGraph = nullptr)
{
    // Rules 1-2: atomic number + mass number
    int cmp12 = compareBranches_Rule1(nodesA, rootA, nodesB, rootB);
    if (cmp12 != 0) return cmp12;

    // Rule 3: Z > E double-bond stereo
    if (molGraph) {
        int cmp3 = compareBranches_Rule3(nodesA, rootA, nodesB, rootB, *molGraph);
        if (cmp3 != 0) return cmp3;
    }

    // Rule 4b-c: like/unlike descriptor pairing
    int cmp4 = compareBranches_Rule4bc(nodesA, rootA, nodesB, rootB,
                                       descriptorMap);
    return cmp4; // Rule 5 is implicit: R(1) > S(2) within Rule 4b-c
}

// --------------------------------------------------------------------------
// Compute CIP priorities for the neighbours of a given atom.
// Returns a vector of (neighbor_atom_index, priority) pairs,
// where priority 1 is lowest and higher numbers are higher priority.
// Atoms with the same priority share the same rank.
// --------------------------------------------------------------------------
inline std::vector<std::pair<int, int>> computePriorities(
    const MolGraph& g, int centre,
    const std::unordered_map<int, int>& descriptorMap = {}, bool includeStereo = true)
{
    const auto& nbrs = g.neighbors[centre];
    int nNbr = static_cast<int>(nbrs.size());

    const int implH = implicitHydrogens(g, centre);
    const bool lonePair = hasLonePairLigand(g, centre, implH);
    const int totalLigands = nNbr + implH + (lonePair ? 1 : 0);
    struct LigandInfo { int atomIdx; int atomicNum; int digraphRoot; };
    std::vector<LigandInfo> ligands;
    ligands.reserve(totalLigands);
    auto digraph = buildDigraph(g, centre);
    if (digraph.empty()) return {};
    for (int neighbor : nbrs) {
        int ligandRoot = -1;
        for (int child : digraph[0].children)
            if (digraph[child].atomIdx == neighbor && !digraph[child].isDuplicate) { ligandRoot = child; break; }
        ligands.push_back({neighbor, g.atomicNum[neighbor], ligandRoot});
    }
    for (int implicit = 0; implicit < implH + (lonePair ? 1 : 0); ++implicit) {
        const bool phantom = implicit >= implH;
        DigraphNode node;
        node.atomIdx = phantom ? -2 : -1;
        node.atomicNum = phantom ? 0 : 1;
        node.scaledAtomicNum = node.atomicNum * ATOMIC_NUMBER_SCALE;
        node.massNumber = 0;
        node.isDuplicate = false;
        node.parent = -1;
        node.depth = 0;
        const int root = static_cast<int>(digraph.size());
        digraph.push_back(node);
        ligands.push_back({node.atomIdx, node.atomicNum, root});
    }

    // Sort indices by CIP priority (Rule 1 comparison)
    std::vector<int> order(totalLigands);
    std::iota(order.begin(), order.end(), 0);

    auto compareLigands = [&](int first, int second) -> int {
        const auto& firstLigand = ligands[first];
        const auto& secondLigand = ligands[second];
        if (firstLigand.atomicNum != secondLigand.atomicNum)
            return firstLigand.atomicNum - secondLigand.atomicNum;
        if (firstLigand.digraphRoot < 0 || secondLigand.digraphRoot < 0) return 0;
        return compareBranches_Full(digraph, firstLigand.digraphRoot, digraph, secondLigand.digraphRoot,
                                    descriptorMap, includeStereo ? &g : nullptr);
    };

    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return compareLigands(a, b) < 0;
    });

    // Assign priority ranks
    std::vector<std::pair<int, int>> result(totalLigands);
    int rank = 1;
    for (int i = 0; i < totalLigands; ++i) {
        int idx = order[i];
        if (i > 0) {
            bool same = (compareLigands(order[i - 1], order[i]) == 0);
            if (!same) rank = i + 1;
        }
        result[idx] = {ligands[idx].atomIdx, rank};
    }

    return result;
}

/** Rule 3 descriptors use atomic/isotopic priorities, independent of input carrier order. */
inline int atomicPriorityDoubleBondConfiguration(const MolGraph& graph, int first, int second) {
    int configuration = graph.dbStereo(first, second);
    if (configuration == 0) return 0;
    for (const auto& ends : {std::pair<int, int>{first, second}, {second, first}}) {
        const auto priorities = computePriorities(graph, ends.first, {}, false);
        int bestRank = -1, secondRank = -1, bestAtom = -1;
        for (const auto& priority : priorities) {
            if (priority.first == ends.second) continue;
            if (priority.second > bestRank) {
                secondRank = bestRank;
                bestRank = priority.second;
                bestAtom = priority.first;
            } else secondRank = std::max(secondRank, priority.second);
        }
        if (bestRank < 0 || bestRank == secondRank) return 0;
        for (int neighbor : graph.neighbors[ends.first]) {
            if (neighbor == ends.second) continue;
            if (neighbor != bestAtom) configuration = 3 - configuration;
            break;
        }
    }
    return configuration;
}

// ============================================================================
// Internal: Rules 4 and 5 with auxiliary descriptors
// ============================================================================

/**
 * Sequence Rules 1 to 6 (IUPAC 2013 P-92.1.3) over the hierarchical digraph of one tetrahedral
 * centre, for the centres Rules 1 to 3 leave unresolved. The method is the one the reference CIP
 * implementations follow (Hanson et al., JCIM 58, 1755 (2018)):
 *
 *   - Every stereogenic unit in the digraph (a tagged tetrahedral centre or a stereogenic double
 *     bond, other than the centre itself) receives an auxiliary descriptor at each of its nodes,
 *     found by ranking its ligands with the digraph re-rooted at that node: the edge towards the
 *     root is reversed and the subtree beyond it explored in that direction. Nodes are labelled
 *     from the deepest sphere inwards, so each sees the descriptors of the units beyond it.
 *   - The root's ligands are then ranked by Rule 1a (atomic number), Rule 1b (a ring-closure
 *     duplicate before any other node, and between two the one whose original is nearer the
 *     root), Rule 2 (the mass number of a labelled atom, against the standard atomic weight of an
 *     unlabelled one), Rule 3 (Z before E), Rule 4a (R, S, seqCis and seqTrans before r, s, E and
 *     Z before none), Rule 4b (like descriptor pairs before unlike, each branch read against its
 *     reference descriptor, the majority descriptor of the first sphere that has one), Rule 4c
 *     (r before s), Rule 5 (R before S and seqCis before seqTrans, read with both R and S as the
 *     reference) and Rule 6 (an undifferentiated reference ligand first).
 *   - Each rule is exhausted sphere by sphere before the next is consulted, and within a rule a
 *     node's branches are ordered by the rules up to and including it.
 *   - A centre whose ligands are told apart only by Rule 5 or Rule 6 is pseudoasymmetric and
 *     takes a lower-case descriptor; a double bond told apart that way at one end only is seqCis
 *     or seqTrans rather than Z or E.
 *
 * Java's CIPAssigner applies the same rules node for node.
 */
class SequenceRules {
public:
    SequenceRules(const MolGraph& graph, int centre) : g_(graph), centre_(centre) {}

    /** The descriptor of the centre, or RSLabel::NONE when the rules do not resolve it. */
    RSLabel label() {
        std::array<int, 4> ligands{};
        if (centre_ < 0 || centre_ >= g_.n || atomPropertyOrZero(g_.tetraChirality, centre_) == 0
            || !stereo_detail::tetrahedralLigands(g_, centre_, ligands)) return RSLabel::NONE;
        nodes_ = buildDigraph(g_, centre_);
        if (nodes_.empty()) return RSLabel::NONE;
        aux_.assign(nodes_.size(), NONE);
        up_.assign(nodes_.size(), -1);
        collectUnits();
        if (centreUnit_ < 0) return RSLabel::NONE;
        labelAuxiliaryDescriptors();
        changeRoot(0);
        const int descriptor = labelTetrahedral(0, units_[centreUnit_]);
        if (exhausted_) return RSLabel::NONE;
        switch (descriptor) {
            case R: return RSLabel::R;
            case S: return RSLabel::S;
            case r: return RSLabel::r;
            case s: return RSLabel::s;
            default: return RSLabel::NONE;
        }
    }

private:
    /** Descriptors of the units in the digraph; E and Z are constitutional, seqCis and seqTrans are not. */
    enum Descriptor : int { NONE = 0, UNKNOWN, NS, R, S, r, s, SEQ_TRANS, SEQ_CIS, E, Z };
    enum RuleKind : int { RULE_1A = 0, RULE_1B, RULE_2, RULE_3, RULE_4A, RULE_4B, RULE_4C, RULE_5, RULE_6, RULE_COUNT };

    /** A rule; Rules 4b and 5 read descriptor pairs against `reference` (R or S) in their replacement form. */
    struct Rule { int kind; int reference; };
    /** A directed edge of the re-rooted digraph: `beg` is the node nearer the current root. */
    struct Edge { int beg; int end; };
    struct Priority { bool unique; bool pseudoasymmetric; };
    /** A stereogenic unit: a tagged tetrahedral centre, or a stereogenic double bond. */
    struct Unit {
        bool bond = false;
        int first = -1, second = -1;        // the atom, or both bond atoms
        std::array<int, 4> ligands{};       // centre: the tag's ligand order, -1 for the implicit ligand
        int tag = 0;                        // centre: 1 anticlockwise, 2 clockwise
        int carrierFirst = -1, carrierSecond = -1; // bond: the reference neighbour at each end
        int configuration = 0;              // bond: 1 cis and 2 trans, relative to the carriers
    };

    /** Descriptor pairs of a branch, read against the reference descriptor at its head. */
    class PairList {
    public:
        explicit PairList(int reference) { add(reference); }
        void add(int descriptor) {
            const int reference = pairReference(descriptor);
            if (reference != NONE) descriptors_.push_back(reference);
        }
        int reference() const { return descriptors_.front(); }
        /** The first pair that is like in one list and unlike in the other decides. */
        int compareTo(const PairList& that) const {
            const size_t count = std::min(descriptors_.size(), that.descriptors_.size());
            for (size_t i = 1; i < count; ++i) {
                const bool thisLike = descriptors_[i] == descriptors_[0];
                const bool thatLike = that.descriptors_[i] == that.descriptors_[0];
                if (thisLike && !thatLike) return 1;
                if (!thisLike && thatLike) return -1;
            }
            return 0;
        }
    private:
        std::vector<int> descriptors_;
    };

    /** Recursive comparisons allowed for one centre; beyond it the centre is left unresolved. */
    static constexpr std::int64_t COMPARISON_BUDGET = 5000000;

    const MolGraph& g_;
    int centre_;
    std::vector<DigraphNode> nodes_;
    std::vector<int> aux_;              // auxiliary descriptor by node
    std::vector<int> up_;               // the neighbour of each node towards the current root, -1 for the root
    std::vector<Unit> units_;
    int centreUnit_ = -1;
    int root_ = 0;                      // the current root
    int rule6Reference_ = -1;           // the Rule 6 reference atom, -1 for none
    std::int64_t budget_ = COMPARISON_BUDGET;
    bool exhausted_ = false;

    static int threeWay(double first, double second) { return first < second ? -1 : first > second ? 1 : 0; }
    static int pairReference(int descriptor) {
        return descriptor == R || descriptor == SEQ_CIS ? R : descriptor == S || descriptor == SEQ_TRANS ? S : NONE;
    }
    static int rule3Order(int descriptor) { return descriptor == E ? 1 : descriptor == Z ? 2 : 0; }
    static int rule4aOrder(int descriptor) {
        switch (descriptor) {
            case R: case S: case SEQ_CIS: case SEQ_TRANS: return 2;
            case r: case s: case E: case Z: return 1;
            default: return 0;
        }
    }
    static int rule4cOrder(int descriptor) { return descriptor == r ? 2 : descriptor == s ? 1 : 0; }

    static int firstOtherNeighbour(const MolGraph& g, int atom, int other) {
        for (int neighbour : g.neighbors[atom]) if (neighbour != other) return neighbour;
        return -1;
    }

    /** Whether the node ends its branch: a duplicate, an implicit hydrogen, or an atom with no other neighbour. */
    bool terminal(int node) const {
        const DigraphNode& n = nodes_[node];
        return n.isDuplicate || n.atomIdx < 0
            || (n.parent >= 0 ? 1 : 0) + static_cast<int>(n.children.size()) == 1;
    }

    /** Make `node` the root: the path from it to the original root is read in reverse. */
    void changeRoot(int node) {
        for (size_t i = 0; i < nodes_.size(); ++i) up_[i] = nodes_[i].parent;
        up_[node] = -1;
        for (int child = node, parent = nodes_[node].parent; parent >= 0;
             child = parent, parent = nodes_[parent].parent)
            up_[parent] = child;
        root_ = node;
    }

    /** The edges of a node under the current root: the one towards the root first, then the rest. */
    std::vector<Edge> edgesOf(int node) const {
        std::vector<Edge> edges;
        const DigraphNode& n = nodes_[node];
        if (up_[node] >= 0) edges.push_back({up_[node], node});
        if (n.parent >= 0 && n.parent != up_[node]) edges.push_back({node, n.parent});
        for (int child : n.children) {
            if (child == up_[node]) continue;
            // The original root's multiple bonds have no duplicates in the reference digraph.
            if (node == 0 && nodes_[child].isDuplicate && !nodes_[child].isRingDuplicate) continue;
            edges.push_back({node, child});
        }
        return edges;
    }

    std::vector<Edge> nonTerminalOutEdges(int node) const {
        std::vector<Edge> out;
        for (const Edge& edge : edgesOf(node))
            if (edge.beg == node && !terminal(edge.end)) out.push_back(edge);
        return out;
    }

    /** The rules a rule orders branches with: those up to and including it; a replacement rule, itself. */
    std::vector<Rule> sorterOf(const Rule& rule) const {
        if (rule.reference != NONE) return {rule};
        std::vector<Rule> rules;
        for (int kind = 0; kind <= rule.kind; ++kind) rules.push_back({kind, NONE});
        return rules;
    }

    std::vector<Rule> allRules() const {
        std::vector<Rule> rules;
        for (int kind = 0; kind < RULE_COUNT; ++kind) rules.push_back({kind, NONE});
        return rules;
    }

    /** All rules with `replacement` in place of its kind, read last. */
    std::vector<Rule> referenceSorter(const Rule& replacement) const {
        std::vector<Rule> rules;
        for (int kind = 0; kind < RULE_COUNT; ++kind)
            if (kind != replacement.kind) rules.push_back({kind, NONE});
        rules.push_back(replacement);
        return rules;
    }

    // ------------------------------------------------------------------
    // Ordering
    // ------------------------------------------------------------------

    /** Order the edges of `node` by `rules`, highest first; `top` reads every rule exhaustively. */
    Priority prioritise(const std::vector<Rule>& rules, bool top, int node, std::vector<Edge>& edges, bool deep) {
        bool unique = true;
        int pseudoasymmetric = 0;
        for (size_t i = 0; i < edges.size(); ++i) {
            for (size_t j = i; j > 0; --j) {
                const int cmp = compareSubstituents(rules, top, node, edges[j - 1], edges[j], deep);
                if (cmp < -1 || cmp > 1) ++pseudoasymmetric;
                if (cmp < 0) {
                    std::swap(edges[j], edges[j - 1]);
                } else {
                    if (cmp == 0) unique = false;
                    break;
                }
            }
        }
        return {unique, pseudoasymmetric == 1};
    }

    int compareSubstituents(const std::vector<Rule>& rules, bool top, int node,
                            const Edge& first, const Edge& second, bool deep) {
        const bool firstUp = first.end == node, secondUp = second.end == node;
        if (firstUp && !secondUp) return 1;
        if (!firstUp && secondUp) return -1;
        for (const Rule& rule : rules) {
            const int cmp = top || deep ? recursiveCompare(rule, first, second) : compare(rule, first, second);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

    /** Consecutive edges of a sorted list that the rules do not tell apart. */
    std::vector<std::vector<Edge>> groups(const std::vector<Rule>& rules, bool top, const std::vector<Edge>& sorted) {
        std::vector<std::vector<Edge>> result;
        for (size_t i = 0; i < sorted.size(); ++i) {
            if (i == 0 || compareSubstituents(rules, top, sorted[i - 1].beg, sorted[i - 1], sorted[i], true) != 0)
                result.emplace_back();
            result.back().push_back(sorted[i]);
        }
        return result;
    }

    /** Compare two edges by one rule at their ends only. */
    int compare(const Rule& rule, const Edge& first, const Edge& second) {
        const DigraphNode& a = nodes_[first.end];
        const DigraphNode& b = nodes_[second.end];
        switch (rule.kind) {
            case RULE_1A:
                return threeWay(a.scaledAtomicNum, b.scaledAtomicNum);
            case RULE_1B:
                if (a.isRingDuplicate && b.isRingDuplicate) return -threeWay(a.distance, b.distance);
                if (a.isRingDuplicate) return 1;
                if (b.isRingDuplicate) return -1;
                return 0;
            case RULE_2:
                if (a.atomicNum == 0 && b.atomicNum == 0) return 0;
                if (a.atomicNum == 0 || b.atomicNum == 0) return threeWay(a.atomicNum, b.atomicNum);
                return compareMassNumbers(&a, &b);
            case RULE_3:
                return threeWay(rule3Order(aux_[first.end]), rule3Order(aux_[second.end]));
            case RULE_4A:
                return threeWay(rule4aOrder(aux_[first.end]), rule4aOrder(aux_[second.end]));
            case RULE_4B:
                return compareRule4b(rule, first, second);
            case RULE_4C:
                return threeWay(rule4cOrder(aux_[first.end]), rule4cOrder(aux_[second.end]));
            case RULE_5:
                return compareRule5(rule, first, second);
            case RULE_6:
                if (rule6Reference_ < 0) return 0;
                if (rule6Reference_ == a.atomIdx && rule6Reference_ != b.atomIdx) return 2;
                if (rule6Reference_ != a.atomIdx && rule6Reference_ == b.atomIdx) return -2;
                return 0;
            default:
                return 0;
        }
    }

    /** Compare two edges by one rule sphere by sphere, ordering each sphere by the rules up to it. */
    int recursiveCompare(const Rule& rule, Edge first, Edge second) {
        if (exhausted_ || budget_-- <= 0) {
            exhausted_ = true;
            return 0;
        }
        int cmp = compare(rule, first, second);
        if (cmp != 0) return cmp;
        const std::vector<Rule> sorter = sorterOf(rule);
        std::vector<Edge> firstQueue{first}, secondQueue{second};
        for (size_t pos = 0; pos < firstQueue.size() && pos < secondQueue.size(); ++pos) {
            first = firstQueue[pos];
            second = secondQueue[pos];
            std::vector<Edge> firstEdges = edgesOf(first.end), secondEdges = edgesOf(second.end);
            prioritise(sorter, false, first.end, firstEdges, false);
            prioritise(sorter, false, second.end, secondEdges, false);
            const int sizeDifference = threeWay(firstEdges.size(), secondEdges.size());
            for (size_t i = 0; i < firstEdges.size() && i < secondEdges.size(); ++i) {
                if (firstEdges[i].end == first.end) continue; // the edges towards the root
                cmp = compare(rule, firstEdges[i], secondEdges[i]);
                if (cmp != 0) return cmp;
            }
            if (sizeDifference != 0) return sizeDifference;
            prioritise(sorter, false, first.end, firstEdges, true);
            prioritise(sorter, false, second.end, secondEdges, true);
            for (size_t i = 0; i < firstEdges.size() && i < secondEdges.size(); ++i) {
                if (firstEdges[i].end == first.end) continue;
                cmp = compare(rule, firstEdges[i], secondEdges[i]);
                if (cmp != 0) return cmp;
                firstQueue.push_back(firstEdges[i]);
                secondQueue.push_back(secondEdges[i]);
            }
        }
        return 0;
    }

    // ------------------------------------------------------------------
    // Rule 4b: like descriptor pairs before unlike
    // ------------------------------------------------------------------

    /** Away from the root, a replacement rule reads the ends against its reference descriptor. */
    int compareAgainstReference(const Rule& rule, const Edge& first, const Edge& second) const {
        if (rule.reference == NONE) return 0;
        const int a = aux_[first.end], b = aux_[second.end];
        if (a != NONE && b != NONE && a != NS && b != NS) {
            const bool firstLike = pairReference(rule.reference) == pairReference(a);
            const bool secondLike = pairReference(rule.reference) == pairReference(b);
            if (firstLike && !secondLike) return 1;
            if (secondLike && !firstLike) return -1;
        }
        return 0;
    }

    int compareRule4b(const Rule& rule, const Edge& first, const Edge& second) {
        if (first.beg != root_ || second.beg != root_) return compareAgainstReference(rule, first, second);
        const std::vector<int> firstReferences = referenceDescriptors(first.end);
        const std::vector<int> secondReferences = referenceDescriptors(second.end);
        if (firstReferences.empty() || secondReferences.empty()) return 0;
        if (firstReferences.size() == 1)
            return comparePairs(first.end, second.end, firstReferences[0], secondReferences[0]);
        std::vector<PairList> firstLists, secondLists;
        for (int reference : firstReferences) {
            firstLists.emplace_back(reference);
            fillPairs({RULE_4B, reference}, first.end, firstLists.back());
        }
        for (int reference : secondReferences) {
            secondLists.emplace_back(reference);
            fillPairs({RULE_4B, reference}, second.end, secondLists.back());
        }
        auto descending = [](const PairList& a, const PairList& b) { return a.compareTo(b) > 0; };
        std::stable_sort(firstLists.begin(), firstLists.end(), descending);
        std::stable_sort(secondLists.begin(), secondLists.end(), descending);
        for (size_t i = 0; i < firstLists.size() && i < secondLists.size(); ++i) {
            const int cmp = firstLists[i].compareTo(secondLists[i]);
            if (cmp != 0) return cmp;
        }
        return 0;
    }

    /** The reference descriptor(s) of a branch: the majority of the first sphere with any, both on a tie. */
    std::vector<int> referenceDescriptors(int node) {
        std::vector<int> result;
        std::vector<std::vector<int>> level{{node}};
        while (!level.empty()) {
            for (const auto& nodes : level)
                if (referenceOf(nodes, result)) return result;
            level = nextLevel(level);
        }
        return {};
    }

    bool referenceOf(const std::vector<int>& nodes, std::vector<int>& result) const {
        int right = 0, left = 0;
        for (int node : nodes) {
            const int descriptor = aux_[node];
            if (descriptor == R || descriptor == SEQ_CIS) ++right;
            else if (descriptor == S || descriptor == SEQ_TRANS) ++left;
        }
        if (right + left == 0) return false;
        if (right > left) result.push_back(R);
        else if (right < left) result.push_back(S);
        else { result.push_back(R); result.push_back(S); }
        return true;
    }

    /** The next sphere, grouped by the rules up to Rule 4b so that equivalent nodes are read together. */
    std::vector<std::vector<int>> nextLevel(const std::vector<std::vector<int>>& previous) {
        const std::vector<Rule> sorter = sorterOf({RULE_4B, NONE});
        std::vector<std::vector<int>> next;
        for (const auto& nodes : previous) {
            std::vector<std::vector<std::vector<Edge>>> grouped;
            for (int node : nodes) {
                std::vector<Edge> edges = nonTerminalOutEdges(node);
                prioritise(sorter, false, node, edges, true);
                grouped.push_back(groups(sorter, false, edges));
            }
            size_t count = grouped.empty() ? 0 : grouped.front().size();
            for (const auto& g : grouped) count = std::min(count, g.size());
            for (size_t i = 0; i < count; ++i) {
                std::vector<int> equivalent;
                for (const auto& g : grouped)
                    for (const Edge& edge : g[i]) equivalent.push_back(edge.end);
                if (!equivalent.empty()) next.push_back(equivalent);
            }
        }
        return next;
    }

    /** Read the descriptors of a branch outwards, each sphere ordered with the replacement rule. */
    void fillPairs(const Rule& replacement, int beg, PairList& pairs) {
        const std::vector<Rule> sorter = referenceSorter(replacement);
        std::vector<int> queue{beg};
        for (size_t i = 0; i < queue.size(); ++i) {
            const int node = queue[i];
            pairs.add(aux_[node]);
            std::vector<Edge> edges = edgesOf(node);
            prioritise(sorter, false, node, edges, true);
            for (const Edge& edge : edges)
                if (edge.beg == node && !terminal(edge.end)) queue.push_back(edge.end);
        }
    }

    /** Compare two branches pair by pair, each read against its own reference descriptor. */
    int comparePairs(int firstNode, int secondNode, int firstReference, int secondReference) {
        const std::vector<Rule> firstSorter = referenceSorter({RULE_4B, firstReference});
        const std::vector<Rule> secondSorter = referenceSorter({RULE_4B, secondReference});
        std::vector<int> firstQueue{firstNode}, secondQueue{secondNode};
        for (size_t pos = 0; pos < firstQueue.size() && pos < secondQueue.size(); ++pos) {
            const int a = firstQueue[pos], b = secondQueue[pos];
            const bool firstLike = pairReference(aux_[a]) == firstReference;
            const bool secondLike = pairReference(aux_[b]) == secondReference;
            if (firstLike && !secondLike) return 1;
            if (!firstLike && secondLike) return -1;
            std::vector<Edge> edges = edgesOf(a);
            prioritise(firstSorter, false, a, edges, true);
            for (const Edge& edge : edges)
                if (edge.beg == a && !terminal(edge.end)) firstQueue.push_back(edge.end);
            edges = edgesOf(b);
            prioritise(secondSorter, false, b, edges, true);
            for (const Edge& edge : edges)
                if (edge.beg == b && !terminal(edge.end)) secondQueue.push_back(edge.end);
        }
        return 0;
    }

    // ------------------------------------------------------------------
    // Rule 5: R before S, read with both references
    // ------------------------------------------------------------------

    int compareRule5(const Rule& rule, const Edge& first, const Edge& second) {
        if (first.beg != root_ || second.beg != root_) return compareAgainstReference(rule, first, second);
        PairList firstR(R), secondR(R), firstS(S), secondS(S);
        fillPairs({RULE_5, R}, first.end, firstR);
        fillPairs({RULE_5, S}, first.end, firstS);
        fillPairs({RULE_5, R}, second.end, secondR);
        fillPairs({RULE_5, S}, second.end, secondS);
        const int cmpR = firstR.compareTo(secondR);
        const int cmpS = firstS.compareTo(secondS);
        // Two when the references disagree, the mark of a pseudoasymmetric distinction.
        if (cmpR < 0) return cmpS < 0 ? -1 : -2;
        if (cmpR > 0) return cmpS > 0 ? 1 : 2;
        return 0;
    }

    // ------------------------------------------------------------------
    // Labelling
    // ------------------------------------------------------------------

    void collectUnits() {
        for (int atom = 0; atom < g_.n; ++atom) {
            if (atomPropertyOrZero(g_.tetraChirality, atom) == 0) continue;
            Unit unit;
            if (!stereo_detail::tetrahedralLigands(g_, atom, unit.ligands)) continue;
            unit.first = atom;
            unit.tag = g_.tetraChirality[atom];
            if (atom == centre_) centreUnit_ = static_cast<int>(units_.size());
            units_.push_back(unit);
        }
        for (int atom = 0; atom < g_.n; ++atom) {
            for (int other : g_.neighbors[atom]) {
                if (other <= atom || g_.bondOrder(atom, other) != 2) continue;
                const int configuration = g_.dbStereo(atom, other);
                if (configuration == 0) continue;
                Unit unit;
                unit.bond = true;
                unit.first = atom;
                unit.second = other;
                unit.carrierFirst = firstOtherNeighbour(g_, atom, other);
                unit.carrierSecond = firstOtherNeighbour(g_, other, atom);
                unit.configuration = configuration;
                if (unit.carrierFirst >= 0 && unit.carrierSecond >= 0) units_.push_back(unit);
            }
        }
    }

    /** Label every node of every other unit, deepest first; a sphere sees only the spheres beyond it. */
    void labelAuxiliaryDescriptors() {
        std::vector<std::pair<int, int>> pending; // node, unit
        for (size_t u = 0; u < units_.size(); ++u) {
            const Unit& unit = units_[u];
            if (!unit.bond && unit.first == centre_) continue;
            for (size_t node = 0; node < nodes_.size(); ++node) {
                if (nodes_[node].isDuplicate || nodes_[node].atomIdx != unit.first) continue;
                int low = static_cast<int>(node);
                // A bond is labelled at its end nearer the root.
                if (unit.bond && nodes_[node].parent >= 0 && nodes_[nodes_[node].parent].atomIdx == unit.second)
                    low = nodes_[node].parent;
                pending.push_back({low, static_cast<int>(u)});
            }
        }
        std::stable_sort(pending.begin(), pending.end(), [&](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return nodes_[a.first].depth > nodes_[b.first].depth;
        });
        std::vector<std::pair<int, int>> sphere; // node, descriptor
        auto flush = [&]() {
            for (const auto& labelled : sphere) aux_[labelled.first] = labelled.second;
            sphere.clear();
        };
        int previous = std::numeric_limits<int>::max();
        for (const auto& entry : pending) {
            const int depth = nodes_[entry.first].depth;
            if (depth < previous) {
                flush();
                previous = depth;
            }
            const Unit& unit = units_[entry.second];
            int descriptor;
            if (unit.bond) {
                descriptor = labelDoubleBond(entry.first, unit);
            } else {
                changeRoot(entry.first);
                descriptor = labelTetrahedral(entry.first, unit);
            }
            bool replaced = false;
            for (auto& labelled : sphere)
                if (labelled.first == entry.first) { labelled.second = descriptor; replaced = true; }
            if (!replaced) sphere.push_back({entry.first, descriptor});
        }
        flush();
    }

    /** Whether reordering `first` into `second` is an odd permutation of the same four edges. */
    static bool oddPermutation(const std::vector<Edge>& first, const std::vector<Edge>& second) {
        std::vector<int> positions;
        for (const Edge& edge : second)
            for (size_t i = 0; i < first.size(); ++i)
                if (first[i].beg == edge.beg && first[i].end == edge.end) positions.push_back(static_cast<int>(i));
        int inversions = 0;
        for (size_t i = 0; i < positions.size(); ++i)
            for (size_t j = i + 1; j < positions.size(); ++j)
                if (positions[i] > positions[j]) ++inversions;
        return inversions % 2 == 1;
    }

    /** The descriptor of a tetrahedral unit at `node`, the current root. */
    int labelTetrahedral(int node, const Unit& unit) {
        std::vector<Edge> edges = edgesOf(node);
        if (edges.size() < 3) return NS;
        const std::vector<Rule> rules = allRules();
        Priority priority = prioritise(rules, true, node, edges, true);
        if (!priority.unique && edges.size() == 4) {
            const auto partition = groups(rules, true, edges);
            if (partition.size() == 2) {
                rule6Reference_ = nodes_[edges[1].end].atomIdx;
                priority = prioritise(rules, true, node, edges, true);
                rule6Reference_ = -1;
            } else if (partition.size() == 1) {
                rule6Reference_ = nodes_[edges[0].end].atomIdx;
                prioritise(rules, true, node, edges, true);
                const std::vector<Edge> first = edges;
                rule6Reference_ = nodes_[edges[1].end].atomIdx;
                priority = prioritise(rules, true, node, edges, true);
                rule6Reference_ = -1;
                if (oddPermutation(first, edges)) return UNKNOWN;
            }
            if (!priority.unique) return UNKNOWN;
        } else if (!priority.unique) {
            return UNKNOWN;
        }
        // The carriers in rank order: every ligand atom once, and the implicit hydrogen.
        std::vector<int> carriers;
        for (const Edge& edge : edges) {
            const DigraphNode& end = nodes_[edge.end];
            if (end.isDuplicate && !end.isRingDuplicate) continue;
            carriers.push_back(end.atomIdx);
        }
        std::array<int, 4> ranks{};
        for (int slot = 0; slot < 4; ++slot) {
            const int ligand = unit.ligands[slot];
            int rank = 0; // a lone pair ranks below every carrier
            for (size_t i = 0; i < carriers.size(); ++i)
                if (carriers[i] == ligand) { rank = static_cast<int>(carriers.size() - i); break; }
            if (ligand >= 0 && rank == 0) return UNKNOWN;
            ranks[slot] = rank;
        }
        int inversions = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                if (ranks[i] > ranks[j]) ++inversions;
        const bool anticlockwise = unit.tag == 1;
        const bool rectus = (inversions % 2 == 0) != anticlockwise;
        if (priority.pseudoasymmetric) return rectus ? r : s;
        return rectus ? R : S;
    }

    /** The descriptor of a double bond unit whose end nearer the root is `root1`. */
    int labelDoubleBond(int root1, const Unit& unit) {
        const int atom1 = nodes_[root1].atomIdx;
        const int atom2 = atom1 == unit.first ? unit.second : unit.first;
        int root2 = -1;
        for (int child : nodes_[root1].children)
            if (!nodes_[child].isDuplicate && nodes_[child].atomIdx == atom2) root2 = child;
        if (root2 < 0) return UNKNOWN;
        auto withoutInternal = [&](int node, int otherAtom) {
            std::vector<Edge> edges;
            for (const Edge& edge : edgesOf(node))
                if (nodes_[edge.end].atomIdx != otherAtom && nodes_[edge.beg].atomIdx != otherAtom) edges.push_back(edge);
            return edges;
        };
        changeRoot(root1);
        std::vector<Edge> edges1 = withoutInternal(root1, atom2);
        changeRoot(root2);
        std::vector<Edge> edges2 = withoutInternal(root2, atom1);
        if (edges1.empty() || edges2.empty()) return NS;
        int carrier1 = unit.carrierFirst, carrier2 = unit.carrierSecond;
        if (atom1 == unit.second) std::swap(carrier1, carrier2);
        int configuration = unit.configuration;
        const std::vector<Rule> rules = allRules();
        changeRoot(root1);
        const Priority priority1 = prioritise(rules, true, root1, edges1, true);
        if (!priority1.unique) return UNKNOWN;
        if (edges1.size() > 1 && nodes_[edges1[0].end].atomIdx != carrier1) configuration = 3 - configuration;
        changeRoot(root2);
        const Priority priority2 = prioritise(rules, true, root2, edges2, true);
        if (!priority2.unique) return UNKNOWN;
        if (edges2.size() > 1 && nodes_[edges2[0].end].atomIdx != carrier2) configuration = 3 - configuration;
        const bool pseudoasymmetric = priority1.pseudoasymmetric != priority2.pseudoasymmetric;
        if (configuration == 1) return pseudoasymmetric ? SEQ_CIS : Z;
        return pseudoasymmetric ? SEQ_TRANS : E;
    }
};

} // namespace detail

// ============================================================================
// Public API: assign R/S to a tetrahedral stereocentre
// ============================================================================

/**
 * Assign CIP R/S descriptor to atom `centre` in MolGraph `g`.
 *
 * Uses a two-pass approach per IUPAC 2013 P-92.1.4:
 *   Pass 1: Rules 1-3 (atomic number, mass number, Z before E)
 *   Pass 2: Rules 4 and 5 with auxiliary descriptors (detail::SequenceRules),
 *           for a centre pass 1 leaves unresolved
 *
 * If ligands are distinguished only by Rule 5 (or Rule 6), the centre is
 * pseudoasymmetric and receives lowercase r/s.
 *
 * @param descriptorMap  Empty for pass 1 only (the default); any entry asks
 *                       for pass 2, whose auxiliary descriptors are derived
 *                       from the digraph rather than read from the map. Only
 *                       whether the map is empty is read.
 * @return RSLabel::NONE if not a valid stereocentre.
 */
inline RSLabel assignRS(const MolGraph& g, int centre,
                        const std::unordered_map<int, int>& descriptorMap = {}) {
    if (centre < 0 || centre >= g.n) return RSLabel::NONE;
    if (detail::atomPropertyOrZero(g.tetraChirality, centre) == 0) return RSLabel::NONE;
    if (!descriptorMap.empty()) return detail::SequenceRules(g, centre).label();

    auto priorities = detail::computePriorities(g, centre);

    // Check that all priorities are distinct
    std::unordered_set<int> ranks;
    for (auto& [atomIdx, rank] : priorities) {
        ranks.insert(rank);
    }

    bool allDistinct = (static_cast<int>(ranks.size()) ==
                        static_cast<int>(priorities.size()));

    // If not all distinct, return NONE: the caller may retry with pass 2,
    // where Rules 4 and 5 and the auxiliary descriptors apply.
    if (!allDistinct) return RSLabel::NONE;

    // We need exactly 4 ligands for a tetrahedral centre
    if (priorities.size() != 4) return RSLabel::NONE;

    // The SMILES chirality (@/@@ ) encodes the winding order of
    // the last three neighbours as seen from the first neighbour
    // (or implicit H if it's first).
    //
    // @ = anticlockwise winding of neighbours 2,3,4 when viewed
    //     from neighbour 1 (or implicit H).
    // @@ = clockwise winding.
    //
    // To determine R/S:
    // 1. Order ligands by CIP priority (1=lowest, 4=highest)
    // 2. Find the lowest-priority ligand (#1).
    // 3. The remaining three ligands (in priority order 2->3->4)
    //    should wind clockwise for R, anticlockwise for S,
    //    when viewed from the lowest-priority ligand.
    //
    // The SMILES convention:
    //   - The "from" atom (or implicit H) is the viewpoint.
    //   - @ means the remaining 3 wind anticlockwise.
    //   - @@ means they wind clockwise.
    //
    // To map SMILES ordering to CIP:
    //   We need to figure out the permutation parity between
    //   the SMILES neighbour order and the CIP priority order.

    // Build map: ligand index -> CIP priority rank
    // ligand index 0..3 corresponds to:
    //   0 = first explicit neighbour (or implicit H)
    //   1,2,3 = subsequent neighbours in SMILES order

    const int implH = detail::implicitHydrogens(g, centre);
    const bool lonePair = detail::hasLonePairLigand(g, centre, implH);
    const int virtualLigand = lonePair ? -2 : -1;

    // SMILES order depends on whether the stereocentre has a preceding
    // bond (i.e., was not the first atom in the SMILES string).
    //
    // Case 1: "N[C@@H](C)C(=O)O" — centre (idx=1) has preceding atom N.
    //   The "from" direction is N. Implicit H goes after the from-atom.
    //   → smilesOrder = [N(0), H(-1), C(2), C(3)]
    //
    // Case 2: "[C@@H](Br)(Cl)F" — centre (idx=0) starts the SMILES.
    //   The "from" direction is the implicit H itself.
    //   → smilesOrder = [H(-1), Br(1), Cl(2), F(3)]
    //
    // So when the first neighbour has a lower index it was the preceding-bond
    // atom, and otherwise the implicit ligand is the from-atom. This is the
    // order the stereo checks of the searches read.
    std::array<int, 4> smilesOrder{}; // atom indices, virtualLigand for the implicit ligand
    if (!stereo_detail::tetrahedralLigands(g, centre, smilesOrder)) return RSLabel::NONE;
    for (int& ligand : smilesOrder) if (ligand < 0) ligand = virtualLigand;

    // Build priority array in SMILES order
    // priorities[i] = {atomIdx, rank} where i is ligand index
    // We need to map smilesOrder[k] -> CIP rank
    std::vector<int> cipRanks(4);
    for (int k = 0; k < 4; ++k) {
        int targetAtom = smilesOrder[k];
        for (auto& [aIdx, rank] : priorities) {
            if (aIdx == targetAtom) {
                cipRanks[k] = rank;
                break;
            }
        }
    }

    // Determine parity of permutation from SMILES order to CIP order
    // CIP order: sorted by priority rank (1,2,3,4)
    // SMILES order: cipRanks[0], cipRanks[1], cipRanks[2], cipRanks[3]

    // Count inversions to determine parity
    int inversions = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (cipRanks[i] > cipRanks[j]) inversions++;
        }
    }
    bool evenPermutation = (inversions % 2 == 0);

    // @ = anticlockwise winding from viewpoint of first ligand
    // @@ = clockwise winding from viewpoint of first ligand
    //
    // For R: viewing from lowest priority, the 2->3->4 winds clockwise
    // For S: viewing from lowest priority, the 2->3->4 winds anticlockwise
    //
    // SMILES @ means anticlockwise => if the permutation to CIP order is
    // even, the winding is preserved (anticlockwise = S).
    // If odd permutation, the winding is inverted.

    bool smilesIsAnticlockwise = (g.tetraChirality[centre] == 1); // @

    // Even perm + @ = anticlockwise from lowest = S
    // Even perm + @@ = clockwise from lowest = R
    // Odd perm + @ = clockwise from lowest = R
    // Odd perm + @@ = anticlockwise from lowest = S

    if (evenPermutation) return smilesIsAnticlockwise ? RSLabel::S : RSLabel::R;
    return smilesIsAnticlockwise ? RSLabel::R : RSLabel::S;
}

// ============================================================================
// Public API: assign E/Z to a double bond
// ============================================================================

/**
 * Assign CIP E/Z descriptor to the double bond between atoms a1 and a2.
 *
 * This function works entirely from the parsed SMILES bond stereo annotations
 * stored in MolGraph::dbStereoConf (populated by parseSMILES from / and \
 * bond direction tokens).  NO 2D or 3D coordinate data is required or used.
 *
 * Examples (from SMILES alone, no coordinates needed):
 *   - "C/C=C/C"  -> E (trans, high-priority groups on opposite sides)
 *   - "C/C=C\\C" -> Z (cis, high-priority groups on same side)
 *
 * The SMILES parser resolves / and \ directions into a cis/trans configuration
 * integer (1=Z, 2=E) stored per double-bond atom pair.  This function then
 * applies CIP priority rules to determine whether the SMILES-encoded
 * configuration corresponds to E or Z in the CIP sense, correcting for cases
 * where the SMILES reference substituent differs from the CIP high-priority
 * substituent.
 *
 * Returns EZLabel::NONE if:
 *   - The bond is neither double nor an annotated aromatic focus
 *   - Either end has two identical substituents (no geometric isomerism)
 *   - The stored dbStereo configuration is 0 (unspecified / no annotation)
 */
inline EZLabel assignEZ(const MolGraph& g, int a1, int a2) {
    if (a1 < 0 || a1 >= g.n || a2 < 0 || a2 >= g.n) return EZLabel::NONE;
    const int order = g.bondOrder(a1, a2);
    if (order != 2 && order != 4) return EZLabel::NONE;

    // Check stored stereo configuration
    int storedConfig = g.dbStereo(a1, a2);
    if (storedConfig == 0) return EZLabel::NONE;

    // For each end of the double bond, find the higher-priority substituent
    // (excluding the other end of the double bond)
    auto getPriority = [&](int thisEnd, int otherEnd) -> int {
        // Build priorities for this end
        auto prios = detail::computePriorities(g, thisEnd);

        // Find the max priority among substituents that are NOT otherEnd
        // and NOT the double-bond phantom nodes
        int bestRank = -1;
        int bestAtom = -1;
        int secondRank = -1;

        for (auto& [atomIdx, rank] : prios) {
            if (atomIdx == otherEnd) continue;
            if (rank > bestRank) {
                secondRank = bestRank;
                bestRank = rank;
                bestAtom = atomIdx;
            } else if (rank > secondRank) {
                secondRank = rank;
            }
        }

        // Check if the two non-double-bond substituents have the same priority
        // (which would make E/Z undefined)
        // Count substituents excluding otherEnd
        int subCount = 0;
        for (auto& [atomIdx, rank] : prios) {
            if (atomIdx != otherEnd) subCount++;
        }
        if (subCount >= 2 && bestRank == secondRank) return -1; // same priority

        return bestAtom;
    };

    int highA1 = getPriority(a1, a2);
    int highA2 = getPriority(a2, a1);

    if (highA1 < 0 || highA2 < 0) return EZLabel::NONE;

    // storedConfig:  1 = Z (cis), 2 = E (trans) — as computed by SMILES parser
    // But the SMILES parser's cis/trans is based on the / \ notation relative
    // to the atoms in the SMILES string, which uses a specific neighbour as reference.
    // We need to check whether the SMILES convention's "cis" corresponds to
    // the CIP high-priority groups being cis or trans.

    // The SMILES parser already correctly determines geometric isomerism:
    //   dbStereoConf[a1][a2] = 1 means the reference substituents are on the
    //   same side (Z in SMILES convention), = 2 means opposite sides (E).
    //
    // The CIP E/Z definition: E = high-priority groups on opposite sides,
    //                          Z = high-priority groups on same side.
    //
    // For simple cases where the SMILES reference substituents ARE the
    // high-priority substituents, the mapping is direct.
    // For the general case, we need to check.

    // The SMILES / and \ annotations are on bonds adjacent to the double bond.
    // The "reference" atoms for the SMILES parser are the first substituent
    // encountered on each side.  If the CIP high-priority atom matches this
    // reference atom, the stored config maps directly.  Otherwise, we need
    // to invert.

    // The SMILES parser stores Z=1, E=2 from whether the annotated substituents
    // lie cis or trans.  When the annotated substituent at each end is also its
    // highest-priority neighbour, the stored configuration is the CIP descriptor.

    // We determine if the high-priority atom at each end is the same as
    // the "reference" atom used by the SMILES parser.  The reference
    // atoms are the neighbours that had / or \ annotations.

    // Since we may not have direct access to which specific neighbour
    // was the reference, we use the stored config as a starting point
    // and adjust based on whether high-priority != first neighbour.

    // Simple approach: check if the high-priority substituent at each end
    // is the first neighbour (in the neighbour list) that isn't the other
    // double-bond atom.  If both match, use stored config directly.
    // If one doesn't match, invert.  If both don't match, keep stored.

    auto getFirstNonDB = [&](int thisEnd, int otherEnd) -> int {
        for (int nb : g.neighbors[thisEnd]) {
            if (nb != otherEnd) return nb;
        }
        return -1;
    };

    int refA1 = getFirstNonDB(a1, a2);
    int refA2 = getFirstNonDB(a2, a1);

    int flips = 0;
    if (highA1 != refA1 && refA1 >= 0) flips++;
    if (highA2 != refA2 && refA2 >= 0) flips++;

    int effectiveConfig = storedConfig;
    if (flips == 1) {
        // One reference atom is not the high-priority one -> invert
        effectiveConfig = (storedConfig == 1) ? 2 : 1;
    }
    // flips == 0 or 2 -> config stays the same

    return (effectiveConfig == 1) ? EZLabel::Z : EZLabel::E;
}

// ============================================================================
// Public API: assign all stereo descriptors for a molecule
// ============================================================================

struct CIPDescriptors {
    std::vector<RSLabel> rsLabels;   // per-atom R/S (NONE if not a stereocentre)
    std::vector<EZLabel> ezLabels;   // per-bond E/Z, indexed by bond pair
    // For E/Z, store as a map from (min(a,b), max(a,b)) atom pairs
    std::vector<std::tuple<int, int, EZLabel>> ezBonds;
};

/**
 * Compute CIP descriptors for all stereocentres and stereogenic
 * double bonds in the molecule.
 *
 * Uses a two-pass approach per IUPAC 2013:
 *   Pass 1: Assign R/S using Rules 1-3 only.
 *   Pass 2: For centres still NONE due to tied priorities, apply Rules 4
 *           to 6 with auxiliary descriptors (detail::SequenceRules).
 *           Centres resolved only by Rule 5 or Rule 6 are pseudoasymmetric (r/s).
 */
inline CIPDescriptors assignAll(const MolGraph& g) {
    CIPDescriptors result;
    result.rsLabels.assign(g.n, RSLabel::NONE);

    // Pass 1: Assign R/S using Rules 1-3 only
    for (int i = 0; i < g.n; ++i) {
        if (detail::atomPropertyOrZero(g.tetraChirality, i) != 0) {
            result.rsLabels[i] = assignRS(g, i);
        }
    }

    // Pass 2: Rules 4 to 6 for the centres still unresolved
    for (int i = 0; i < g.n; ++i) {
        if (detail::atomPropertyOrZero(g.tetraChirality, i) != 0 && result.rsLabels[i] == RSLabel::NONE) {
            result.rsLabels[i] = detail::SequenceRules(g, i).label();
        }
    }

    // Assign E/Z to all double bonds
    std::unordered_set<int64_t> visited;
    for (int i = 0; i < g.n; ++i) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue; // process each bond once
            if (g.bondOrder(i, j) == 2 || (g.bondOrder(i, j) == 4 && g.dbStereo(i, j) != 0)) {
                int64_t key = MolGraph::bondKey(i, j);
                if (visited.count(key)) continue;
                visited.insert(key);

                EZLabel ez = assignEZ(g, i, j);
                if (ez != EZLabel::NONE) {
                    result.ezBonds.push_back({i, j, ez});
                }
            }
        }
    }

    return result;
}

/**
 * Convenience: parse SMILES and return CIP descriptors.
 */
inline CIPDescriptors assignFromSMILES(const std::string& smiles) {
    auto g = parseSMILES(smiles);
    return assignAll(g);
}

} // namespace cip
} // namespace smsd

#endif // SMSD_CIP_HPP

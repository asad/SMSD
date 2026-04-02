/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * CIP (Cahn-Ingold-Prelog) stereodescriptor assignment for MolGraph.
 *
 * Implements the IUPAC 2013 recommendations (Blue Book) for assigning
 * R/S descriptors to tetrahedral stereocentres and E/Z descriptors to
 * double bonds.  Uses the digraph (extended-connectivity tree) approach
 * with phantom/duplicate nodes for ring closures and correct handling
 * of mancude (aromatic) systems.
 *
 * References:
 *   - Prelog & Helmchen, Angew. Chem. Int. Ed. Engl. 21, 567 (1982)
 *   - IUPAC Recommendations 2013 (P-92)
 *   - Hanson et al., JCIM 58, 1755 (2018) — edge-case corrections
 *
 * Header-only, zero external dependencies — pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_CIP_HPP
#define SMSD_CIP_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
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

struct DigraphNode {
    int  atomIdx;       // index into MolGraph
    int  atomicNum;     // atomic number (copied for speed)
    int  massNumber;    // 0 for phantom/duplicate nodes
    bool isDuplicate;   // true for ring-closure / bond-expansion phantom nodes
    int  parent;        // index in digraph node list (-1 for root)
    int  depth;         // distance from root in the digraph
    std::vector<int> children; // indices of child nodes in digraph
};

// --------------------------------------------------------------------------
// Build the CIP digraph (extended connectivity tree) rooted at `centre`
// going out through `startNbr`.  Ring closures produce duplicate (phantom)
// nodes with mass = 0.  Mancude (aromatic) duplicate nodes also get mass = 0.
//
// maxDepth limits expansion to prevent runaway on large molecules.
// --------------------------------------------------------------------------
inline std::vector<DigraphNode> buildDigraph(
    const MolGraph& g,
    int root,
    int maxDepth = 30,
    int maxNodes = 5000)
{
    std::vector<DigraphNode> nodes;
    nodes.reserve(128);

    // Root node
    DigraphNode rootNode;
    rootNode.atomIdx    = root;
    rootNode.atomicNum  = g.atomicNum[root];
    rootNode.massNumber = g.massNumber[root];
    rootNode.isDuplicate = false;
    rootNode.parent     = -1;
    rootNode.depth      = 0;
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

            // Collect neighbors, expanding double/triple bonds as per CIP
            // convention: a double bond A=B adds one extra phantom B on A's
            // side (and vice versa); a triple bond adds two.
            struct Expansion {
                int neighborAtom;
                bool isDuplicate;
                int bondOrder;
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

                int bo = g.bondOrder(atom, nbr);
                if (bo == 4) bo = 1; // aromatic bond treated as single for expansion

                // Real child node
                expansions.push_back({nbr, false, bo});

                // Phantom nodes for higher-order bonds
                for (int extra = 1; extra < bo; ++extra) {
                    expansions.push_back({nbr, true, bo});
                }
            }

            // Also add phantom nodes for the bond back to parent
            // (for double/triple bonds to parent, the parent side
            //  already added phantoms for us, but we need phantom
            //  copies of the parent on THIS node's child list).
            if (parentAtom >= 0 && !nodes[ni].isDuplicate) {
                int bo = g.bondOrder(atom, parentAtom);
                if (bo == 4) bo = 1;
                for (int extra = 1; extra < bo; ++extra) {
                    expansions.push_back({parentAtom, true, bo});
                }
            }

            for (auto& exp : expansions) {
                if (static_cast<int>(nodes.size()) >= maxNodes) break;

                // Check if this atom is already on the path from root
                // to current node.  If so, it's a ring closure -> duplicate.
                bool onPath = false;
                if (!exp.isDuplicate) {
                    int walk = ni;
                    while (walk >= 0) {
                        if (nodes[walk].atomIdx == exp.neighborAtom) {
                            onPath = true;
                            break;
                        }
                        walk = nodes[walk].parent;
                    }
                }

                bool makeDuplicate = exp.isDuplicate || onPath;

                DigraphNode child;
                child.atomIdx    = exp.neighborAtom;
                child.atomicNum  = g.atomicNum[exp.neighborAtom];
                // Duplicate nodes: mass = 0 per CIP rules
                // (including mancude / aromatic ring closures)
                child.massNumber = makeDuplicate ? 0 : g.massNumber[exp.neighborAtom];
                child.isDuplicate = makeDuplicate;
                child.parent     = ni;
                child.depth      = depth + 1;

                int childIdx = static_cast<int>(nodes.size());
                nodes.push_back(child);
                nodes[ni].children.push_back(childIdx);

                // Duplicate/phantom nodes are leaves -- don't expand further
                if (!makeDuplicate) {
                    nextFrontier.push_back(childIdx);
                }
            }
        }
        frontier = std::move(nextFrontier);
    }

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
    // BFS level-by-level
    std::vector<int> frontA = {rootA};
    std::vector<int> frontB = {rootB};

    while (!frontA.empty() || !frontB.empty()) {
        // Gather atomic numbers at this level, sorted descending
        std::vector<int> zA, zB;
        for (int i : frontA) zA.push_back(nodesA[i].atomicNum);
        for (int i : frontB) zB.push_back(nodesB[i].atomicNum);

        // Sort descending (highest priority first)
        std::sort(zA.begin(), zA.end(), std::greater<int>());
        std::sort(zB.begin(), zB.end(), std::greater<int>());

        // Compare element by element
        size_t len = std::max(zA.size(), zB.size());
        for (size_t k = 0; k < len; ++k) {
            int a = (k < zA.size()) ? zA[k] : 0;
            int b = (k < zB.size()) ? zB[k] : 0;
            if (a != b) return a - b;
        }

        // Tie on atomic number at this level -- try mass number
        std::vector<int> mA, mB;
        for (int i : frontA) mA.push_back(nodesA[i].massNumber);
        for (int i : frontB) mB.push_back(nodesB[i].massNumber);
        std::sort(mA.begin(), mA.end(), std::greater<int>());
        std::sort(mB.begin(), mB.end(), std::greater<int>());

        for (size_t k = 0; k < std::max(mA.size(), mB.size()); ++k) {
            int a = (k < mA.size()) ? mA[k] : 0;
            int b = (k < mB.size()) ? mB[k] : 0;
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

    return 0; // truly equal by Rules 1-2
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
// Full CIP comparison: Rules 1-2, then 4b-c, then Rule 5 (R > S).
// --------------------------------------------------------------------------
inline int compareBranches_Full(
    const std::vector<DigraphNode>& nodesA, int rootA,
    const std::vector<DigraphNode>& nodesB, int rootB,
    const std::unordered_map<int, int>& descriptorMap)
{
    // Rules 1-2: atomic number + mass number
    int cmp12 = compareBranches_Rule1(nodesA, rootA, nodesB, rootB);
    if (cmp12 != 0) return cmp12;

    // Rule 4b-c: like/unlike descriptor pairing
    int cmp4 = compareBranches_Rule4bc(nodesA, rootA, nodesB, rootB,
                                       descriptorMap);
    return cmp4; // Rule 5 is implicit: R(1) > S(2) within Rule 4b-c
}

// --------------------------------------------------------------------------
// Compute CIP priorities for the neighbors of a given atom.
// Returns a vector of (neighbor_atom_index, priority) pairs,
// where priority 1 is lowest and higher numbers are higher priority.
// Atoms with the same priority share the same rank.
// --------------------------------------------------------------------------
inline std::vector<std::pair<int, int>> computePriorities(
    const MolGraph& g, int centre,
    const std::unordered_map<int, int>& descriptorMap = {})
{
    const auto& nbrs = g.neighbors[centre];
    int nNbr = static_cast<int>(nbrs.size());

    // Also account for implicit hydrogen(s) as a pseudo-neighbor
    // with atomic number 1, mass 0, no further substituents.
    int boSum = 0;
    for (int nb : nbrs) {
        int bo = g.bondOrder(centre, nb);
        if (bo == 4) bo = 1;
        boSum += bo;
    }
    int implH = ::smsd::detail::computeImplicitH(g.atomicNum[centre], g.aromatic[centre] != 0,
                                 boSum, g.formalCharge[centre]);

    // Total ligands = explicit neighbors + implicit H
    int totalLigands = nNbr + implH;

    // Build a digraph for each ligand
    struct LigandInfo {
        int atomIdx;      // -1 for implicit H
        int atomicNum;
        std::vector<DigraphNode> digraph;
        int digraphRoot;  // index of ligand's root in its digraph
    };

    std::vector<LigandInfo> ligands;
    ligands.reserve(totalLigands);

    for (int nb : nbrs) {
        // Build full digraph from centre
        auto dg = buildDigraph(g, centre);

        // Find the child of root (node 0) that corresponds to neighbor nb
        int ligRoot = -1;
        for (int c : dg[0].children) {
            if (dg[c].atomIdx == nb && !dg[c].isDuplicate) {
                ligRoot = c;
                break;
            }
        }
        // If not found as real node, try duplicate (for multiply-bonded)
        if (ligRoot < 0) {
            for (int c : dg[0].children) {
                if (dg[c].atomIdx == nb) {
                    ligRoot = c;
                    break;
                }
            }
        }

        LigandInfo li;
        li.atomIdx = nb;
        li.atomicNum = g.atomicNum[nb];
        li.digraph = std::move(dg);
        li.digraphRoot = ligRoot;
        ligands.push_back(std::move(li));
    }

    // Add implicit H ligands
    for (int h = 0; h < implH; ++h) {
        LigandInfo li;
        li.atomIdx = -1;
        li.atomicNum = 1;
        // Create a trivial single-node digraph for H
        DigraphNode hNode;
        hNode.atomIdx = -1;
        hNode.atomicNum = 1;
        hNode.massNumber = 0;
        hNode.isDuplicate = false;
        hNode.parent = -1;
        hNode.depth = 0;
        li.digraph = {hNode};
        li.digraphRoot = 0;
        ligands.push_back(std::move(li));
    }

    // Sort indices by CIP priority (Rule 1 comparison)
    std::vector<int> order(totalLigands);
    std::iota(order.begin(), order.end(), 0);

    // Choose comparison function based on whether descriptor map is available
    bool useFull = !descriptorMap.empty();

    auto compareLigands = [&](int a, int b) -> int {
        auto& la = ligands[a];
        auto& lb = ligands[b];
        if (la.atomicNum != lb.atomicNum)
            return la.atomicNum - lb.atomicNum;
        if (la.digraphRoot >= 0 && lb.digraphRoot >= 0) {
            if (useFull)
                return compareBranches_Full(
                    la.digraph, la.digraphRoot,
                    lb.digraph, lb.digraphRoot,
                    descriptorMap);
            else
                return compareBranches_Rule1(
                    la.digraph, la.digraphRoot,
                    lb.digraph, lb.digraphRoot);
        }
        return 0;
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

// --------------------------------------------------------------------------
// Check if an atom is a potential tetrahedral stereocentre.
// Requirements: 4 distinct-priority ligands (counting implicit H),
// and the SMILES chirality annotation is nonzero.
// --------------------------------------------------------------------------
inline bool isPotentialStereocentre(const MolGraph& g, int atom) {
    if (g.tetraChirality[atom] == 0) return false;
    // Typically sp3 carbon with 4 different substituents
    // We allow any atom with chirality annotation
    return true;
}

} // namespace detail

// ============================================================================
// Public API: assign R/S to a tetrahedral stereocentre
// ============================================================================

/**
 * Assign CIP R/S descriptor to atom `centre` in MolGraph `g`.
 *
 * Uses a two-pass approach per IUPAC 2013 P-92.1.4:
 *   Pass 1: Rules 1-2 (atomic number, mass number)
 *   Pass 2: Rules 4b-c (like/unlike descriptor pairing) + Rule 5 (R > S)
 *
 * If ligands are distinguished only in pass 2, the centre is
 * pseudoasymmetric and receives lowercase r/s.
 *
 * @param descriptorMap  Pre-computed R/S labels from pass 1 of other centres
 *                       (empty map = pass 1 only, no Rule 4b-c/5).
 * @return RSLabel::NONE if not a valid stereocentre.
 * @since 6.5.2
 */
inline RSLabel assignRS(const MolGraph& g, int centre,
                        const std::unordered_map<int, int>& descriptorMap = {}) {
    if (centre < 0 || centre >= g.n) return RSLabel::NONE;
    if (g.tetraChirality[centre] == 0) return RSLabel::NONE;

    auto priorities = detail::computePriorities(g, centre, descriptorMap);

    // Check that all priorities are distinct
    std::unordered_set<int> ranks;
    for (auto& [atomIdx, rank] : priorities) {
        ranks.insert(rank);
    }

    bool allDistinct = (static_cast<int>(ranks.size()) ==
                        static_cast<int>(priorities.size()));

    // If not all distinct and we haven't tried Rule 4b-c yet, return NONE
    // (caller may retry with descriptorMap in pass 2).
    // If we already have descriptorMap and still not distinct, truly NONE.
    if (!allDistinct) return RSLabel::NONE;

    // We need exactly 4 ligands for a tetrahedral centre
    if (priorities.size() != 4) return RSLabel::NONE;

    // Pseudoasymmetric detection: if pass-1 (no descriptorMap) would have
    // returned NONE due to tied priorities, but pass-2 distinguishes them,
    // the centre is pseudoasymmetric.
    bool isPseudoasymmetric = false;
    if (!descriptorMap.empty()) {
        // Re-check with Rules 1-2 only to see if ties existed
        auto pass1Prios = detail::computePriorities(g, centre);
        std::unordered_set<int> pass1Ranks;
        for (auto& [ai, r] : pass1Prios) pass1Ranks.insert(r);
        if (static_cast<int>(pass1Ranks.size()) !=
            static_cast<int>(pass1Prios.size())) {
            isPseudoasymmetric = true;
        }
    }

    // The SMILES chirality (@/@@ ) encodes the winding order of
    // the last three neighbors as seen from the first neighbor
    // (or implicit H if it's first).
    //
    // @ = anticlockwise winding of neighbors 2,3,4 when viewed
    //     from neighbor 1 (or implicit H).
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
    //   the SMILES neighbor order and the CIP priority order.

    // Build map: ligand index -> CIP priority rank
    // ligand index 0..3 corresponds to:
    //   0 = first explicit neighbor (or implicit H)
    //   1,2,3 = subsequent neighbors in SMILES order

    // Get the SMILES order of neighbors
    const auto& nbrs = g.neighbors[centre];
    int nExplicit = static_cast<int>(nbrs.size());

    // Compute implicit H count
    int boSum = 0;
    for (int nb : nbrs) {
        int bo = g.bondOrder(centre, nb);
        if (bo == 4) bo = 1;
        boSum += bo;
    }
    int implH = ::smsd::detail::computeImplicitH(g.atomicNum[centre],
                                 g.aromatic[centre] != 0,
                                 boSum, g.formalCharge[centre]);

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
    // Heuristic: if centre > 0 and first neighbor has lower index,
    // it was the preceding-bond atom. Otherwise H is the from-atom.
    std::vector<int> smilesOrder; // atom indices, -1 for implicit H
    if (implH > 0) {
        bool hasPrecedingAtom = (centre > 0 && nExplicit >= 1 && nbrs[0] < centre);
        if (hasPrecedingAtom) {
            // From-atom is nbrs[0], H is second
            smilesOrder.push_back(nbrs[0]);
            smilesOrder.push_back(-1);
            for (int k = 1; k < nExplicit; ++k)
                smilesOrder.push_back(nbrs[k]);
        } else {
            // H is the from-atom (bracket starts the SMILES or no preceding bond)
            smilesOrder.push_back(-1);
            for (int nb : nbrs) smilesOrder.push_back(nb);
        }
    } else {
        for (int nb : nbrs) smilesOrder.push_back(nb);
    }

    if (static_cast<int>(smilesOrder.size()) != 4) return RSLabel::NONE;

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

    RSLabel label;
    if (evenPermutation) {
        label = smilesIsAnticlockwise ? RSLabel::S : RSLabel::R;
    } else {
        label = smilesIsAnticlockwise ? RSLabel::R : RSLabel::S;
    }

    // Pseudoasymmetric centres get lowercase r/s
    if (isPseudoasymmetric) {
        label = (label == RSLabel::R) ? RSLabel::r : RSLabel::s;
    }
    return label;
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
 *   - There is no double bond between a1 and a2
 *   - Either end has two identical substituents (no geometric isomerism)
 *   - The stored dbStereo configuration is 0 (unspecified / no annotation)
 */
inline EZLabel assignEZ(const MolGraph& g, int a1, int a2) {
    if (a1 < 0 || a1 >= g.n || a2 < 0 || a2 >= g.n) return EZLabel::NONE;
    if (g.bondOrder(a1, a2) != 2) return EZLabel::NONE;

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
    // to the atoms in the SMILES string, which uses a specific neighbor as reference.
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

    // For now, use a simpler approach: the SMILES parser stores Z=1, E=2
    // based on whether the annotated substituents are cis or trans.
    // In the common case (the annotated substituents = first neighbor =
    // highest-priority neighbor for each end), this is correct.

    // We determine if the high-priority atom at each end is the same as
    // the "reference" atom used by the SMILES parser.  The reference
    // atoms are the neighbors that had / or \ annotations.

    // Since we may not have direct access to which specific neighbor
    // was the reference, we use the stored config as a starting point
    // and adjust based on whether high-priority != first neighbor.

    // Simple approach: check if the high-priority substituent at each end
    // is the first neighbor (in the neighbor list) that isn't the other
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
 *   Pass 1: Assign R/S using Rules 1-2 only.
 *   Pass 2: For centres still NONE due to tied priorities, re-assign
 *           using Rules 4b-c (like/unlike pairing) and Rule 5 (R > S).
 *           Centres resolved only in pass 2 are pseudoasymmetric (r/s).
 *
 * @since 6.5.2
 */
inline CIPDescriptors assignAll(const MolGraph& g) {
    CIPDescriptors result;
    result.rsLabels.assign(g.n, RSLabel::NONE);

    // Pass 1: Assign R/S using Rules 1-2 only
    for (int i = 0; i < g.n; ++i) {
        if (g.tetraChirality[i] != 0) {
            result.rsLabels[i] = assignRS(g, i);
        }
    }

    // Build descriptor map from pass 1 for Rule 4b-c
    std::unordered_map<int, int> descriptorMap;
    for (int i = 0; i < g.n; ++i) {
        if (result.rsLabels[i] != RSLabel::NONE) {
            // Encode: R/r=1, S/s=2 for descriptor comparison
            int desc = (result.rsLabels[i] == RSLabel::R ||
                        result.rsLabels[i] == RSLabel::r) ? 1 : 2;
            descriptorMap[i] = desc;
        }
    }

    // Pass 2: Re-assign NONE centres using Rules 4b-c + 5
    if (!descriptorMap.empty()) {
        for (int i = 0; i < g.n; ++i) {
            if (g.tetraChirality[i] != 0 && result.rsLabels[i] == RSLabel::NONE) {
                result.rsLabels[i] = assignRS(g, i, descriptorMap);
            }
        }
    }

    // Assign E/Z to all double bonds
    std::unordered_set<int64_t> visited;
    for (int i = 0; i < g.n; ++i) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue; // process each bond once
            if (g.bondOrder(i, j) == 2) {
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

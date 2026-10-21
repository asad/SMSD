/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Operations on whole MolGraphs: induced subgraphs, ring-system scaffolds and the
 * translation of mappings to external atom IDs, with the deprecated ChemOps rules.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// ChemOps -- atom and bond compatibility checks (static utility, deprecated)
//
// The searches apply their own atom and bond rules (search_detail::atomCompatible
// and bondCompatible); nothing in the library calls these.
// ============================================================================

struct [[deprecated("Unused by the searches, whose rules differ; MCS::validateMapping checks a mapping.")]] ChemOps {

    /// Atom compatibility check for substructure / MCS matching.
    static bool atomsCompatible(
            const MolGraph& gq, int qi,
            const MolGraph& gt, int tj,
            const ChemOptions& C) {
        if (C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj])
            return false;
        if (C.matchFormalCharge && gq.formalCharge[qi] != gt.formalCharge[tj])
            return false;
        if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
            && gq.aromatic[qi] != gt.aromatic[tj])
            return false;
        if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj])
            return false;
        if (C.matchIsotope) {
            int qm = gq.massNumber[qi], tm = gt.massNumber[tj];
            if (qm != 0 && tm != 0 && qm != tm) return false;
        }
        // Relative tetrahedral tags require ligand-order checks by the matcher.
        if (C.ringFusionMode == ChemOptions::RingFusionMode::STRICT
            && gq.ring[qi] && gt.ring[tj]) {
            if (gq.ringCount[qi] != gt.ringCount[tj]) return false;
        }
        return true;
    }

    /// Bond compatibility check for substructure / MCS matching.
    static bool bondsCompatible(
            const MolGraph& g1, int qi, int qk,
            const MolGraph& g2, int tj, int tk,
            const ChemOptions& C) {
        int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
        if (qOrd == 0 || tOrd == 0) return false;

        // Tautomer-aware relaxation affects bond order only.
        bool tautomerBond = false;
        if (C.tautomerAware
            && !g1.tautomerClass.empty() && !g2.tautomerClass.empty()) {
            bool qBothTaut = g1.tautomerClass[qi] != -1
                          && g1.tautomerClass[qk] != -1;
            bool tBothTaut = g2.tautomerClass[tj] != -1
                          && g2.tautomerClass[tk] != -1;
            tautomerBond = qBothTaut && tBothTaut;
        }

        // Strict aromaticity: check bond aromaticity directly (no ring guard)
        if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
            && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk))
            return false;



        if (C.ringMatchesRingOnly
            && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk))
            return false;

        if (tautomerBond || C.matchBondOrder == ChemOptions::BondOrderMode::ANY)
            return true;
        if (qOrd == tOrd) return true;
        if (C.matchBondOrder == ChemOptions::BondOrderMode::LOOSE)
            return true;

        if (C.aromaticityMode == ChemOptions::AromaticityMode::FLEXIBLE) {
            bool qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
            if ((qa && ta)
                || (qa && (tOrd == 1 || tOrd == 2))
                || (ta && (qOrd == 1 || qOrd == 2)))
                return true;
        }
        return false;
    }
};

// ===========================================================================
// extractSubgraph -- extract a subgraph containing only the specified atom indices.
// Extract an induced graph in the supplied atom order.
// ===========================================================================

inline MolGraph extractSubgraph(const MolGraph& mol, const std::vector<int>& atomIndices) {
    for (int atom : atomIndices)
        if (atom < 0 || atom >= mol.n)
            throw std::out_of_range("Subgraph atom index out of range");
    int subN = static_cast<int>(atomIndices.size());
    if (subN == 0) return MolGraph();

    // Build old-to-new index mapping
    std::unordered_map<int, int> oldToNew;
    oldToNew.reserve(subN);
    for (int k = 0; k < subN; ++k)
        oldToNew[atomIndices[k]] = k;

    std::vector<int> atomicNum(subN), formalCharge(subN, 0), massNumber(subN, 0), hydrogenCount(subN, 0), atomClass(subN, 0);
    std::vector<int> atomIds;
    if (!mol.atomId.empty()) atomIds.resize(subN, 0);
    std::vector<uint8_t> ring(subN, uint8_t(0)), aromatic(subN, uint8_t(0));
    std::vector<std::vector<int>> neighbors(subN);
    std::vector<std::vector<int>> bondOrders(subN);
    std::vector<std::vector<bool>> bondRings(subN);
    std::vector<std::vector<bool>> bondAroms(subN);

    for (int k = 0; k < subN; ++k) {
        int oldIdx = atomIndices[k];
        atomicNum[k]    = mol.atomicNum[oldIdx];
        formalCharge[k]  = mol.formalCharge[oldIdx];
        massNumber[k]    = mol.massNumber[oldIdx];
        if (oldIdx < static_cast<int>(mol.hydrogenCount.size()))
            hydrogenCount[k] = mol.hydrogenCount[oldIdx];
        if (oldIdx < static_cast<int>(mol.atomClass.size()))
            atomClass[k] = mol.atomClass[oldIdx];
        if (!atomIds.empty() && oldIdx < static_cast<int>(mol.atomId.size()))
            atomIds[k] = mol.atomId[oldIdx];
        // Do NOT copy ring/aromatic flags from parent — they may be invalid
        // for the subgraph (e.g., partial ring extraction). These will be
        // re-perceived below after bond connectivity is established.
        ring[k]      = 0;
        aromatic[k]  = 0;
    }

    // Copy bonds where both endpoints are in the subgraph
    for (int k = 0; k < subN; ++k) {
        int oldI = atomIndices[k];
        for (int oldJ : mol.neighbors[oldI]) {
            auto it = oldToNew.find(oldJ);
            if (it == oldToNew.end()) continue;
            int newJ = it->second;
            neighbors[k].push_back(newJ);
            bondOrders[k].push_back(mol.bondOrder(oldI, oldJ));
            bondRings[k].push_back(mol.bondInRing(oldI, oldJ));
            bondAroms[k].push_back(mol.bondAromatic(oldI, oldJ));
        }
    }

    // Build the subgraph WITHOUT ring/aromatic flags first
    auto builder = MolGraph::Builder()
        .atomCount(subN)
        .atomicNumbers(std::move(atomicNum))
        .formalCharges(std::move(formalCharge))
        .massNumbers(std::move(massNumber))
        .hydrogenCounts(std::move(hydrogenCount))
        .atomClasses(std::move(atomClass))
        .ringFlags(std::move(ring))        // all zeros — re-perceived below
        .aromaticFlags(std::move(aromatic)) // all zeros — re-perceived below
        .setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(bondOrders))
        .bondRingFlags(std::move(bondRings))
        .bondAromaticFlags(std::move(bondAroms))
        .name(mol.name)
        .programLine(mol.programLine)
        .comment(mol.comment)
        .properties(mol.properties);
    if (!atomIds.empty()) builder.atomIds(std::move(atomIds));
    auto sub = builder.build();
    sub.hydrogenCountsKnown = mol.hydrogenCountsKnown;

    // Re-perceive ring membership and aromaticity on the actual subgraph topology.
    // This prevents aromatic flag leakage from partial ring extraction (e.g., 4 atoms
    // of a 6-membered ring are NOT aromatic as a standalone fragment).
    sub.ensureRingCounts();
    for (int i = 0; i < sub.n; ++i) {
        sub.ring[i] = (!sub.ringCount.empty() && sub.ringCount[i] > 0) ? 1 : 0;
        // Simple aromaticity: atom is aromatic only if it is in a ring AND was
        // aromatic in the parent AND all its ring neighbours in the subgraph are present
        if (sub.ring[i]) {
            int oldIdx = atomIndices[i];
            bool allRingNbPresent = true;
            for (int oldNb : mol.neighbors[oldIdx]) {
                if (mol.ring[oldNb] && mol.bondInRing(oldIdx, oldNb)) {
                    if (oldToNew.find(oldNb) == oldToNew.end()) {
                        allRingNbPresent = false;
                        break;
                    }
                }
            }
            sub.aromatic[i] = (mol.aromatic[oldIdx] && allRingNbPresent) ? 1 : 0;
        }
    }

    // Kekule fallback: atoms aromatic in the parent but not in the subgraph (a broken ring,
    // such as part of naphthalene) get single bonds (order 1) and lose the aromatic flag.
    // This prevents invalid aromatic SMILES output for partial ring fragments.
    for (int i = 0; i < sub.n; ++i) {
        int oldIdx = atomIndices[i];
        if (mol.aromatic[oldIdx] && !sub.aromatic[i]) {
            for (int nb : sub.neighbors[i]) {
                int oldNb = atomIndices[nb];
                if (mol.aromatic[oldNb] && !sub.aromatic[nb]) {
                    sub.setBondOrder(i, nb, 1);
                    sub.setBondAromatic(i, nb, false);
                }
            }
        }
    }

    return sub;
}

// ===========================================================================
// murckoScaffold -- extract the Bemis-Murcko scaffold (ring systems + linkers).
// Keeps the ring atoms, the atoms on paths between rings, and every atom joined
// to one of these by a double or triple bond; all other atoms are stripped.
// Double-bonded atoms are kept as in Bemis and Murcko (1996) and RDKit.
// ===========================================================================

inline MolGraph murckoScaffold(const MolGraph& mol) {
    if (mol.n == 0) return mol;

    // Identify ring atoms
    bool hasRing = false;
    for (int i = 0; i < mol.n; ++i) {
        if (mol.ring[i]) { hasRing = true; break; }
    }
    if (!hasRing) return mol;

    // Mark atoms to keep: all ring atoms + atoms on shortest paths between ring atoms
    std::vector<bool> keep(mol.n, false);
    std::vector<int> ringAtoms;
    for (int i = 0; i < mol.n; ++i) {
        if (mol.ring[i]) {
            keep[i] = true;
            ringAtoms.push_back(i);
        }
    }

    // BFS from each ring atom to find shortest paths to other ring atoms
    // through non-ring atoms (linkers)
    for (int src : ringAtoms) {
        std::vector<int> prev(mol.n, -1);
        std::vector<bool> seen(mol.n, false);
        std::deque<int> queue;
        queue.push_back(src);
        seen[src] = true;

        while (!queue.empty()) {
            int u = queue.front(); queue.pop_front();
            for (int v : mol.neighbors[u]) {
                if (seen[v]) continue;
                seen[v] = true;
                prev[v] = u;
                if (mol.ring[v] && v != src) {
                    // Trace back the path and mark all atoms on it
                    int cur = v;
                    while (cur != src && cur != -1) {
                        keep[cur] = true;
                        cur = prev[cur];
                    }
                }
                // Only continue BFS through non-ring atoms (looking for linkers)
                if (!mol.ring[v]) queue.push_back(v);
            }
        }
    }

    // An atom joined to a ring or linker atom by a double or triple bond stays with that bond,
    // so that a ring carrying an exocyclic C=O, C=N or C=S group can still be kekulised.
    const std::vector<bool> framework = keep;
    for (int atom = 0; atom < mol.n; ++atom) {
        if (!framework[atom]) continue;
        for (int neighbor : mol.neighbors[atom]) {
            const int order = mol.bondOrder(atom, neighbor);
            if (!framework[neighbor] && (order == 2 || order == 3)) keep[neighbor] = true;
        }
    }

    // Collect kept atom indices
    std::vector<int> keptIndices;
    for (int i = 0; i < mol.n; ++i) {
        if (keep[i]) keptIndices.push_back(i);
    }

    if (static_cast<int>(keptIndices.size()) == mol.n) return mol;

    // Every ring stays whole and only acyclic atoms are stripped, so the ring and aromaticity
    // perception of the input still holds and is copied rather than perceived again. As in the
    // Java scaffold, a tetrahedral centre keeps its configuration when all its neighbours stay, and
    // a double bond when each end keeps a ligand besides the other end.
    const int count = static_cast<int>(keptIndices.size());
    std::vector<int> newIndex(mol.n, -1);
    for (int k = 0; k < count; ++k) newIndex[keptIndices[k]] = k;
    auto property = [](const std::vector<int>& values, int atom) {
        return atom < static_cast<int>(values.size()) ? values[atom] : 0;
    };
    std::vector<int> atomicNum(count), formalCharge(count), massNumber(count), hydrogenCount(count), atomClass(count);
    std::vector<int> chirality(count, 0), atomIds;
    if (!mol.atomId.empty()) atomIds.resize(count, 0);
    std::vector<uint8_t> ring(count), aromatic(count);
    std::vector<std::vector<int>> neighbors(count), bondOrders(count);
    std::vector<std::vector<bool>> bondRings(count), bondAroms(count);
    std::vector<bool> complete(count, true);
    for (int k = 0; k < count; ++k) {
        const int original = keptIndices[k];
        atomicNum[k] = mol.atomicNum[original];
        formalCharge[k] = property(mol.formalCharge, original);
        massNumber[k] = property(mol.massNumber, original);
        hydrogenCount[k] = property(mol.hydrogenCount, original);
        atomClass[k] = property(mol.atomClass, original);
        if (!atomIds.empty()) atomIds[k] = property(mol.atomId, original);
        ring[k] = mol.ring[original];
        aromatic[k] = mol.aromatic[original];
        for (int neighbor : mol.neighbors[original]) {
            if (newIndex[neighbor] < 0) {
                complete[k] = false;
                continue;
            }
            neighbors[k].push_back(newIndex[neighbor]);
            bondOrders[k].push_back(mol.bondOrder(original, neighbor));
            bondRings[k].push_back(mol.bondInRing(original, neighbor));
            bondAroms[k].push_back(mol.bondAromatic(original, neighbor));
        }
        if (complete[k]) chirality[k] = property(mol.tetraChirality, original);
    }
    // A double-bond configuration is stored against the first ligand of each end besides the
    // other end, so it flips when that ligand is stripped and the second ligand of the end stays.
    auto sameReference = [&](int end, int other) {
        int first = -1, kept = -1;
        for (int ligand : mol.neighbors[end]) {
            if (ligand == other) continue;
            if (first < 0) first = ligand;
            if (newIndex[ligand] >= 0) {
                kept = ligand;
                break;
            }
        }
        if (kept < 0) return 0;
        if (kept == first) return 1;
        return mol.neighbors[end].size() == 3 ? -1 : 0;
    };
    std::vector<std::tuple<int, int, int>> stereoBonds;
    for (int k = 0; k < count; ++k)
        for (int other : neighbors[k]) {
            if (other <= k) continue;
            const int first = keptIndices[k], second = keptIndices[other];
            const int configuration = mol.dbStereo(first, second);
            if (configuration == 0) continue;
            const int relation = sameReference(first, second) * sameReference(second, first);
            if (relation != 0) stereoBonds.emplace_back(k, other, relation > 0 ? configuration : 3 - configuration);
        }
    auto builder = MolGraph::Builder()
        .atomCount(count)
        .atomicNumbers(std::move(atomicNum))
        .formalCharges(std::move(formalCharge))
        .massNumbers(std::move(massNumber))
        .hydrogenCounts(std::move(hydrogenCount))
        .atomClasses(std::move(atomClass))
        .ringFlags(std::move(ring))
        .aromaticFlags(std::move(aromatic))
        .tetrahedralChirality(std::move(chirality))
        .setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(bondOrders))
        .bondRingFlags(std::move(bondRings))
        .bondAromaticFlags(std::move(bondAroms))
        .name(mol.name)
        .programLine(mol.programLine)
        .comment(mol.comment)
        .properties(mol.properties);
    if (!atomIds.empty()) builder.atomIds(std::move(atomIds));
    if (!stereoBonds.empty()) builder.doubleBondStereoSparse(std::move(stereoBonds));
    MolGraph scaffold = builder.build(false);
    scaffold.hydrogenCountsKnown = mol.hydrogenCountsKnown;
    if (!mol.hydrogenCountsKnown) return scaffold;
    // A stripped substituent is replaced by hydrogen, as in a common fragment written as SMILES:
    // an atom that lost neighbours takes the valence model's count for its remaining bonds, and an
    // aromatic N or P keeps its count plus one per stripped neighbour. The Java scaffold caps its
    // atoms the same way.
    for (int atom = 0; atom < scaffold.n; ++atom) {
        const int original = keptIndices[atom];
        const size_t stripped = mol.neighbors[original].size() - scaffold.neighbors[atom].size();
        if (stripped == 0) continue;
        const int element = scaffold.atomicNum[atom];
        if (scaffold.aromatic[atom] && (element == 7 || element == 15)) {
            scaffold.hydrogenCount[atom] = mol.implicitHydrogenCount(original) + static_cast<int>(stripped);
            continue;
        }
        int bondOrderSum = 0;
        for (int neighbor : scaffold.neighbors[atom]) {
            const int order = scaffold.bondOrder(atom, neighbor);
            bondOrderSum += order == 4 ? 1 : order;
        }
        scaffold.hydrogenCount[atom] = element <= 1 ? 0
            : detail::computeImplicitH(element, scaffold.aromatic[atom] != 0, bondOrderSum, scaffold.formalCharge[atom]);
    }
    return scaffold;
}

// ===========================================================================
// translateToAtomIds -- convert MCS mapping from internal indices to external IDs
//
// If either graph has a non-empty atomId vector, the corresponding keys or
// values in the returned map are translated. If neither graph has external
// IDs, the input mapping is returned unchanged.
// ===========================================================================

inline std::map<int,int> translateToAtomIds(
        const std::map<int,int>& mapping,
        const MolGraph& g1,
        const MolGraph& g2) {
    if (mapping.empty()) return mapping;
    if (g1.atomId.empty() && g2.atomId.empty()) return mapping;
    std::map<int,int> translated;
    for (const auto& [k, v] : mapping) {
        if (!g1.atomId.empty() && (k < 0 || k >= static_cast<int>(g1.atomId.size()))) continue;
        if (!g2.atomId.empty() && (v < 0 || v >= static_cast<int>(g2.atomId.size()))) continue;
        int key = g1.atomId.empty() ? k : g1.atomId[k];
        int val = g2.atomId.empty() ? v : g2.atomId[v];
        translated[key] = val;
    }
    return translated;
}

} // namespace smsd

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Atom and bond compatibility under ChemOptions, stereo parity of mappings, and the
 * QueryView that answers both for a molecule or a SMARTS query.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include "smsd/detail/search_graph_validation.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/stereo.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {
namespace search_detail {

inline bool atomCompatible(const MolGraph& gq, int qi,
                            const MolGraph& gt, int tj,
                            const ChemOptions& C) {
    // Tautomerism changes proton positions and bond orders, not elements.
    if (C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj]) return false;
    auto property = [](const auto& values, int atom) { return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : 0; };
    if (C.matchFormalCharge && property(gq.formalCharge, qi) != property(gt.formalCharge, tj)) return false;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT && gq.aromatic[qi] != gt.aromatic[tj]) return false;
    if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
    if (C.matchIsotope) {
        int qm = property(gq.massNumber, qi), tm = property(gt.massNumber, tj);
        if (qm != 0 && tm != 0 && qm != tm) return false;
    }
    if (C.useChirality) {
        // Tags of atoms that are not tetrahedral centres, see stereo_detail::tetrahedralLigands(), are
        // compared as they are; those of centres by their ligands once mapped.
        const int queryTag = property(gq.tetraChirality, qi), targetTag = property(gt.tetraChirality, tj);
        if (queryTag != 0 && targetTag != 0 && queryTag != targetTag) {
            std::array<int, 4> ligands{};
            if (!stereo_detail::tetrahedralLigands(gq, qi, ligands) || !stereo_detail::tetrahedralLigands(gt, tj, ligands))
                return false;
        }
    }
    if (C.ringFusionMode != ChemOptions::RingFusionMode::IGNORE && gq.ring[qi] && gt.ring[tj]) {
        if (C.ringFusionMode == ChemOptions::RingFusionMode::STRICT) {
            if (property(gq.ringCount, qi) != property(gt.ringCount, tj)) return false;
        }
    }
    return true;
}

/**
 * Whether atomCompatible compares atomic numbers alone under these options. Keep it, and
 * atomCompatibilityTerms, in step with atomCompatible; the substructure regression tests check all three
 * against each other for every option that atomCompatible reads.
 */
inline bool atomCompatibleByElementOnly(const ChemOptions& C) {
    return C.matchAtomType && !C.matchFormalCharge && C.aromaticityMode != ChemOptions::AromaticityMode::STRICT
        && !C.ringMatchesRingOnly && !C.matchIsotope && !C.useChirality
        && C.ringFusionMode != ChemOptions::RingFusionMode::STRICT;
}

/**
 * Every atom property atomCompatible reads, followed by the neighbour count. Atoms with equal
 * terms are compatible with the same atoms under any options.
 */
inline std::array<int, 10> atomCompatibilityTerms(const MolGraph& graph, int atom) {
    auto property = [atom](const auto& values) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : 0;
    };
    const int tag = property(graph.tetraChirality);
    std::array<int, 4> ligands{};
    const int centre = tag != 0 && stereo_detail::tetrahedralLigands(graph, atom, ligands) ? 1 : 0;
    return {graph.atomicNum[atom], property(graph.formalCharge), graph.aromatic[atom], graph.ring[atom],
            property(graph.massNumber), tag, centre, graph.degree[atom],
            property(graph.ringCount), static_cast<int>(graph.neighbors[atom].size())};
}

inline bool bondCompatible(const MolGraph& g1, int qi, int qk,
                            const MolGraph& g2, int tj, int tk,
                            const ChemOptions& C) {
    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
    if (qOrd == 0 || tOrd == 0) return false;
    // Tautomer-aware: relax bond ORDER (single/double interchangeable)
    // but still enforce aromaticity, stereo, and ring constraints.
    bool tautBondRelax = C.tautomerAware
        && !g1.tautomerClass.empty() && !g2.tautomerClass.empty()
        && g1.tautomerClass[qi] != -1 && g1.tautomerClass[qk] != -1
        && g2.tautomerClass[tj] != -1 && g2.tautomerClass[tk] != -1;
    // Strict aromaticity: check bond aromaticity directly (no ring guard)
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
        && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk)) return false;
    if (C.ringMatchesRingOnly && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk)) return false;
    if (tautBondRelax) return true;  // bond order relaxed for tautomeric bonds
    if (C.matchBondOrder == ChemOptions::BondOrderMode::ANY) return true;
    if (qOrd == tOrd) return true;
    if (C.matchBondOrder == ChemOptions::BondOrderMode::LOOSE) return true;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::FLEXIBLE) {
        bool qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
        if ((qa && ta) || (qa && (tOrd == 1 || tOrd == 2)) || (ta && (qOrd == 1 || qOrd == 2)))
            return true;
        // Resonance-aware: allow single/double interchange between atoms that
        // are both in a conjugated (aromatic) ring system.  This catches Kekulé
        // form mismatches (e.g. pyridone drawn as C=O vs C-O) without requiring
        // explicit tautomer flags.
        if ((g1.aromatic[qi] && g1.aromatic[qk]) &&
            (g2.aromatic[tj] && g2.aromatic[tk]) &&
            (qOrd == 1 || qOrd == 2) && (tOrd == 1 || tOrd == 2))
            return true;
    }
    return false;
}

/**
 * Everything bondCompatible reads about one bond when tautomers are not used: its order,
 * aromaticity and ring membership, and whether both its atoms are aromatic. Bonds with equal
 * keys are compatible with the same bonds. Keep it in step with bondCompatible.
 */
inline int64_t bondCompatibilityKey(const MolGraph& graph, int first, int second) {
    return int64_t(graph.bondOrder(first, second)) * 8 + (graph.bondAromatic(first, second) ? 4 : 0)
        + (graph.bondInRing(first, second) ? 2 : 0) + (graph.aromatic[first] && graph.aromatic[second] ? 1 : 0);
}

/**
 * The key of bondCompatibilityKey() and, with `tautomers`, whether both atoms of the bond are in
 * tautomer classes, which bondCompatible then also reads.
 */
inline int64_t bondCompatibilityKey(const MolGraph& graph, int first, int second, bool tautomers) {
    const int64_t key = bondCompatibilityKey(graph, first, second);
    const bool tautomeric = tautomers && !graph.tautomerClass.empty()
        && graph.tautomerClass[first] != -1 && graph.tautomerClass[second] != -1;
    return key / 8 * 16 + (tautomeric ? 8 : 0) + key % 8;
}

/**
 * Whether the tetrahedral tags of a mapped atom pair agree over the mapped ligands. An untagged
 * query atom matches any target atom. With `embedding`, as in a substructure search, a tagged
 * query atom needs a tagged target atom, as a SMARTS chirality primitive does; otherwise an
 * untagged target atom matches as well, as in a common-subgraph search.
 */
template <typename Mapping>
inline bool tetraParityCompatible(const MolGraph& gq, int qi,
                                  const MolGraph& gt, int tj,
                                  const Mapping& q2t, bool embedding = false) {
    int qs = qi < static_cast<int>(gq.tetraChirality.size()) ? gq.tetraChirality[qi] : 0;
    int ts = tj < static_cast<int>(gt.tetraChirality.size()) ? gt.tetraChirality[tj] : 0;
    if (qs == 0) return true;
    if (ts == 0) return !embedding;

    std::array<int, 4> queryLigands{}, targetLigands{};
    if (!stereo_detail::tetrahedralLigands(gq, qi, queryLigands)
        || !stereo_detail::tetrahedralLigands(gt, tj, targetLigands)) return true;
    return stereo_detail::tetrahedralParityMatches(qs, queryLigands, ts, targetLigands, q2t);
}

/** The stereo test of a complete mapping; `embedding` is the substructure rule of tetraParityCompatible(). */
template <typename Mapping>
inline bool mappingStereoCompatible(const MolGraph& gq, const MolGraph& gt,
                                    const ChemOptions& opts, const Mapping& q2t, bool embedding = false) {
    if (!opts.useChirality && !opts.useBondStereo) return true;
    for (int queryAtom = 0; queryAtom < gq.n; ++queryAtom) {
        const int targetAtom = q2t[queryAtom];
        if (targetAtom < 0) continue;
        if (opts.useChirality && !tetraParityCompatible(gq, queryAtom, gt, targetAtom, q2t, embedding)) return false;
        if (opts.useBondStereo) {
            for (int neighbor : gq.neighbors[queryAtom])
                if (neighbor > queryAtom
                    && !stereo_detail::doubleBondParityMatches(gq, queryAtom, neighbor, gt, q2t, embedding)) return false;
        }
    }
    return true;
}

/**
 * Spare vectors kept per thread, so the many short searches of a batch reuse their buffers
 * instead of allocating; on many cores the shared allocator, not the search, otherwise limits
 * cheap batches. Only a bounded number of moderately sized buffers is kept.
 */
template<class T>
class SpareVectors {
    static constexpr size_t maximumSpares = 64;
    static constexpr size_t maximumCapacity = size_t(1) << 14;
    static std::vector<std::vector<T>>& spares() {
        thread_local std::vector<std::vector<T>> pool;
        return pool;
    }
public:
    static std::vector<T> take() {
        auto& pool = spares();
        if (pool.empty()) return {};
        std::vector<T> spare = std::move(pool.back());
        pool.pop_back();
        spare.clear();
        return spare;
    }
    static std::vector<T> take(size_t size, const T& value) {
        std::vector<T> spare = take();
        spare.assign(size, value);
        return spare;
    }
    static void give(std::vector<T>& spare) noexcept {
        if (spare.capacity() == 0 || spare.capacity() > maximumCapacity) return;
        auto& pool = spares();
        if (pool.size() >= maximumSpares) return;
        if (pool.capacity() < maximumSpares) {
            try { pool.reserve(maximumSpares); } catch (...) { return; }
        }
        pool.push_back(std::move(spare));
    }
    /** A spare vector that is handed back when the lease ends. */
    class Lease {
    public:
        std::vector<T> items;
        Lease() : items(take()) {}
        Lease(size_t size, const T& value) : items(take(size, value)) {}
        Lease(const Lease&) = delete;
        Lease& operator=(const Lease&) = delete;
        ~Lease() { give(items); }
    };
};

struct QueryView {
    const MolGraph* molecule = nullptr;
    const SmartsQuery* pattern = nullptr;
    /** Filled only when the pattern's adjacency lists do not describe its bonds. */
    std::unordered_map<int64_t, int> patternBondIndex;
    /** Scanning the short adjacency lists finds the same lowest bond index without a per-call map. */
    bool scanPatternAdjacency = false;
    explicit QueryView(const MolGraph& graph) : molecule(&graph) {}
    explicit QueryView(const SmartsQuery& query) : pattern(&query), scanPatternAdjacency(adjacencyListsBonds(query)) {
        if (scanPatternAdjacency) return;
        for (size_t bond = 0; bond < query.bonds.size(); ++bond)
            patternBondIndex.emplace(MolGraph::bondKey(query.bonds[bond].from, query.bonds[bond].to), static_cast<int>(bond));
    }
    /** True when each bond is listed exactly once at each of its atoms, as buildAdjacency does. */
    static bool adjacencyListsBonds(const SmartsQuery& query) {
        const size_t atoms = static_cast<size_t>(query.atomCount());
        if (query.neighbors.size() != atoms || query.bondIndices.size() != atoms) return false;
        // Per bond: entries at its first atom in the low half, at its second atom in the high half.
        std::vector<uint32_t> listed = SpareVectors<uint32_t>::take(query.bonds.size(), 0u);
        bool consistent = true;
        for (size_t atom = 0; atom < atoms && consistent; ++atom) {
            const auto& adjacent = query.neighbors[atom];
            const auto& bonds = query.bondIndices[atom];
            consistent = adjacent.size() == bonds.size();
            for (size_t position = 0; position < adjacent.size() && consistent; ++position) {
                const int bond = bonds[position];
                if (bond < 0 || static_cast<size_t>(bond) >= query.bonds.size()) { consistent = false; break; }
                const auto& entry = query.bonds[bond];
                const int self = static_cast<int>(atom), other = adjacent[position];
                consistent = (entry.from == self && entry.to == other) || (entry.to == self && entry.from == other);
                listed[bond] += (entry.from == self ? 1u : 0u) + (entry.to == self ? 0x10000u : 0u);
            }
        }
        for (size_t bond = 0; bond < listed.size() && consistent; ++bond)
            consistent = listed[bond] == (query.bonds[bond].from == query.bonds[bond].to ? 0x20002u : 0x10001u);
        SpareVectors<uint32_t>::give(listed);
        return consistent;
    }
    /** Lowest index of a bond between the two pattern atoms, or -1. */
    int patternBond(int first, int second) const {
        if (scanPatternAdjacency) {
            const auto& adjacent = pattern->neighbors[first];
            int found = -1;
            for (size_t position = 0; position < adjacent.size(); ++position)
                if (adjacent[position] == second) {
                    const int bond = pattern->bondIndices[first][position];
                    if (found < 0 || bond < found) found = bond;
                }
            return found;
        }
        const auto bond = patternBondIndex.find(MolGraph::bondKey(first, second));
        return bond == patternBondIndex.end() ? -1 : bond->second;
    }
    int size() const { return molecule ? molecule->n : pattern->atomCount(); }
    const std::vector<int>& neighbors(int atom) const {
        return molecule ? molecule->neighbors[atom] : pattern->neighbors[atom];
    }
    bool adjacent(int first, int second) const {
        return molecule ? molecule->hasBond(first, second) : patternBond(first, second) >= 0;
    }
    bool atomMatches(int queryAtom, const MolGraph& target, int targetAtom, const ChemOptions& chemistry) const {
        if (molecule) return atomCompatible(*molecule, queryAtom, target, targetAtom, chemistry);
        const auto truth = smarts_detail::evalAtomExprWithStereo(pattern->atoms[queryAtom].expr, target, targetAtom,
                                                                 [](const smarts_detail::AtomExpr&) {
            return smarts_detail::StereoTruth::UNKNOWN;
        });
        return truth != smarts_detail::StereoTruth::FALSE_VALUE;
    }
    bool bondMatches(int first, int second, const MolGraph& target,
                     int mappedFirst, int mappedSecond, const ChemOptions& chemistry) const {
        if (first > second) {
            std::swap(first, second);
            std::swap(mappedFirst, mappedSecond);
        }
        if (molecule) return bondCompatible(*molecule, first, second, target, mappedFirst, mappedSecond, chemistry);
        const int bond = patternBond(first, second);
        return bond >= 0 && pattern->bondCompat(bond, target, mappedFirst, mappedSecond);
    }
};

} // namespace search_detail
} // namespace smsd

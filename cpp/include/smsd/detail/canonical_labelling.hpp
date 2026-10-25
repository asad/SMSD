/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph canonical labelling: neighbourhood ranks, partition refinement with
 * individualisation, automorphism orbits and the canonical hash.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// Neighbourhood rank computation
// ============================================================================

inline std::vector<int> MolGraph::computeMorganRanks(
        int n, const std::vector<int>& initLabel,
        const std::vector<std::vector<int>>& neighbors) {
    std::vector<int> rank(initLabel.begin(), initLabel.end());
    std::vector<int> rankNew(n);
    std::vector<int> nbSorted;

    for (int iter = 0; iter < 8; iter++) {
        std::unordered_set<int> prevSet(rank.begin(), rank.end());
        int prevDistinct = static_cast<int>(prevSet.size());

        for (int v = 0; v < n; v++) {
            const auto& nbs = neighbors[v];
            int deg = static_cast<int>(nbs.size());
            nbSorted.resize(deg);
            for (int k = 0; k < deg; k++) nbSorted[k] = rank[nbs[k]];
            std::sort(nbSorted.begin(), nbSorted.end());
            // Hash modulo 2^32 explicitly; preserve the signed rank ordering
            // without relying on undefined signed-overflow behaviour.
            uint32_t h = static_cast<uint32_t>(rank[v]) * static_cast<uint32_t>(HASH_PRIME);
            for (int k = 0; k < deg; k++)
                h = h * 31u + static_cast<uint32_t>(nbSorted[k]);
            int32_t signedHash;
            std::memcpy(&signedHash, &h, sizeof(h));
            rankNew[v] = signedHash;
        }
        rank = rankNew;

        std::unordered_set<int> newSet(rank.begin(), rank.end());
        if (static_cast<int>(newSet.size()) <= prevDistinct) break;
    }
    return rank;
}

// ============================================================================
// Orbits from the union-find of automorphic atoms
// ============================================================================

inline std::vector<int> MolGraph::buildOrbits(std::vector<int>& uf, int n) {
    std::vector<int> orb(n);
    for (int i = 0; i < n; i++) orb[i] = detail::unionFindRoot(uf, i);
    std::vector<int> orbitMin(n, n);
    for (int i = 0; i < n; i++) {
        int r = orb[i];
        if (i < orbitMin[r]) orbitMin[r] = i;
    }
    for (int i = 0; i < n; i++) orb[i] = orbitMin[orb[i]];
    return orb;
}

// ============================================================================
// 2-hop neighbourhood signature for orbit disambiguation
// ============================================================================

inline uint64_t MolGraph::twoHopSignature(
        int v, const std::vector<int>& label,
        const std::vector<int>& degree,
        const std::vector<std::vector<int>>& neighbors) {
    // 1-hop: sorted (label, degree) of direct neighbours
    std::vector<uint64_t> sig1;
    sig1.reserve(neighbors[v].size());
    for (int u : neighbors[v]) {
        sig1.push_back(static_cast<uint64_t>(label[u]) * 1000003LL + degree[u]);
    }
    std::sort(sig1.begin(), sig1.end());

    // 2-hop: for each direct neighbour, the sorted (label, degree) of THEIR neighbours
    std::vector<uint64_t> sig2;
    for (int u : neighbors[v]) {
        for (int w : neighbors[u]) {
            if (w == v) continue;  // skip back-link to v itself
            sig2.push_back(static_cast<uint64_t>(label[w]) * 1000003LL + degree[w]);
        }
    }
    std::sort(sig2.begin(), sig2.end());

    // Combine modulo 2^64; orbit refinement only compares signatures.
    uint64_t h = static_cast<uint64_t>(label[v]) * HASH_PRIME + degree[v];
    h = h * HASH_PRIME + static_cast<uint64_t>(sig1.size());
    for (uint64_t x : sig1) h = h * 31 + x;
    h = h * HASH_PRIME + static_cast<uint64_t>(sig2.size());
    for (uint64_t x : sig2) h = h * 31 + x;
    return h;
}

inline std::vector<int> MolGraph::refineOrbitsBy2Hop(
        const std::vector<int>& orb, int n,
        const std::vector<int>& label,
        const std::vector<int>& degree,
        const std::vector<std::vector<int>>& neighbors) {
    // Group atoms by current orbit
    std::unordered_map<int, std::vector<int>> groups;
    for (int i = 0; i < n; i++) groups[orb[i]].push_back(i);

    std::vector<int> refined(n);
    for (auto& [orbitId, members] : groups) {
        if (members.size() <= 1) {
            refined[members[0]] = members[0];
            continue;
        }
        // Compute 2-hop signature for each member
        std::unordered_map<uint64_t, int> sigToRep;
        for (int v : members) {
            uint64_t sig = twoHopSignature(v, label, degree, neighbors);
            auto it = sigToRep.find(sig);
            if (it == sigToRep.end()) {
                sigToRep[sig] = v;
                refined[v] = v;
            } else {
                refined[v] = it->second;
            }
        }
    }
    // Normalise to minimum representative in each refined orbit
    std::unordered_map<int, int> repMin;
    for (int i = 0; i < n; i++) {
        int r = refined[i];
        auto it = repMin.find(r);
        if (it == repMin.end() || i < it->second) repMin[r] = i;
    }
    for (int i = 0; i < n; i++) refined[i] = repMin[refined[i]];
    return refined;
}

// ============================================================================
// Canonical labelling (partition refinement and individualisation)
// ============================================================================

inline int MolGraph::canonicalBondLabel(const MolGraph* graph, int first, int second) {
    if (!graph) return 1;
    return (graph->bondOrder(first, second) << 3)
        | (graph->bondInRing(first, second) ? 4 : 0)
        | (graph->bondAromatic(first, second) ? 2 : 0)
        | (graph->dbStereo(first, second) != 0 ? 1 : 0);
}

inline int MolGraph::canonicalTetrahedralTag(const MolGraph& graph, int atom,
                                    const std::vector<int>& ranking) {
    const int tag = atom < static_cast<int>(graph.tetraChirality.size()) ? graph.tetraChirality[atom] : 0;
    if (tag == 0) return 0;
    const auto& neighbors = graph.neighbors[atom];
    if (neighbors.size() != 3 && neighbors.size() != 4) return tag;
    std::array<int, 4> ligandOrder{};
    const size_t hydrogenPosition = atom > 0 && neighbors[0] < atom ? 1 : 0;
    size_t neighborPosition = 0;
    for (size_t position = 0; position < ligandOrder.size(); ++position)
        ligandOrder[position] = neighbors.size() == 3 && position == hydrogenPosition
            ? -1 : neighbors[neighborPosition++];
    auto rank = [&](int ligand) { return ligand < 0 || ranking.empty() ? ligand : ranking[ligand]; };
    int inversions = 0;
    for (size_t first = 0; first < ligandOrder.size(); ++first)
        for (size_t second = first + 1; second < ligandOrder.size(); ++second)
            inversions += rank(ligandOrder[first]) > rank(ligandOrder[second]);
    return inversions % 2 == 0 ? tag : 3 - tag;
}

inline int MolGraph::canonicalDoubleBondTag(const MolGraph& graph, int first, int second,
                                   const std::vector<int>& ranking) {
    int configuration = graph.dbStereo(first, second);
    if (configuration == 0) return 0;
    for (const auto& end : {std::make_pair(first, second), std::make_pair(second, first)}) {
        int reference = -1, canonicalReference = -1;
        for (int neighbor : graph.neighbors[end.first]) {
            if (neighbor == end.second) continue;
            if (reference < 0) reference = neighbor;
            if (canonicalReference < 0 || (ranking.empty() ? neighbor < canonicalReference
                    : ranking[neighbor] < ranking[canonicalReference])) canonicalReference = neighbor;
        }
        if (reference != canonicalReference) configuration = 3 - configuration;
    }
    return configuration;
}

inline std::vector<int> MolGraph::canonicalAtomLabels() const {
    auto valueAt = [](const std::vector<int>& values, int atom) {
        return atom < static_cast<int>(values.size()) ? values[atom] : 0;
    };
    std::vector<std::array<int, 6>> invariants(n);
    for (int atom = 0; atom < n; ++atom)
        invariants[atom] = {{label[atom], valueAt(formalCharge, atom), valueAt(massNumber, atom),
            valueAt(hydrogenCount, atom), valueAt(atomClass, atom), valueAt(tetraChirality, atom) != 0 ? 1 : 0}};
    std::vector<int> order(n), labels(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int first, int second) { return invariants[first] < invariants[second]; });
    int rank = 0;
    for (int position = 0; position < n; ++position) {
        if (position > 0 && invariants[order[position]] != invariants[order[position - 1]]) ++rank;
        labels[order[position]] = rank;
    }
    return labels;
}

inline int MolGraph::findSmallestNonSingleton(int n, const std::vector<bool>& cellEnd) {
    int bestSize = n + 1, bestEnd = -1, cellStart = 0;
    for (int i = 0; i < n; i++) {
        if (cellEnd[i]) {
            int size = i - cellStart + 1;
            if (size > 1 && size < bestSize) { bestSize = size; bestEnd = i; }
            cellStart = i + 1;
        }
    }
    return bestEnd;
}

inline void MolGraph::refinePartitionImpl(
        int n, std::vector<int>& perm, std::vector<bool>& cellEnd,
        const std::vector<std::vector<int>>& neighbors,
        const std::vector<int>& label,
        const std::vector<int>& degree, int /*targetPos*/, const MolGraph* graph) {
    std::vector<int> inv(n);
    for (int i = 0; i < n; i++) inv[perm[i]] = i;

    bool changed = true;
    int maxIter = n * 2 + 10;
    while (changed && maxIter-- > 0) {
        changed = false;
        int sStart = 0;
        for (int sEnd = 0; sEnd < n; sEnd++) {
            if (!cellEnd[sEnd]) continue;
            int cStart = 0;
            for (int cEnd = 0; cEnd < n; cEnd++) {
                if (!cellEnd[cEnd]) continue;
                int cSize = cEnd - cStart + 1;
                if (cSize <= 1) { cStart = cEnd + 1; continue; }

                std::vector<std::vector<int>> count(cSize);
                for (int ci = 0; ci < cSize; ci++) {
                    int v = perm[cStart + ci];
                    for (int nb : neighbors[v]) {
                        int nbPos = inv[nb];
                        if (nbPos >= sStart && nbPos <= sEnd)
                            count[ci].push_back(canonicalBondLabel(graph, v, nb));
                    }
                    std::sort(count[ci].begin(), count[ci].end());
                }

                bool allSame = true;
                for (int ci = 1; ci < cSize; ci++) {
                    if (count[ci] != count[0]) { allSame = false; break; }
                }
                if (allSame) { cStart = cEnd + 1; continue; }

                // Sort indices within the cell. The comparison ends with the atom index, so it is
                // a total order and every standard library sorts the cell alike.
                std::vector<int> idx(cSize);
                std::iota(idx.begin(), idx.end(), 0);
                const int fcs = cStart;
                std::sort(idx.begin(), idx.end(), [&](int a, int b) {
                    if (count[a] != count[b]) return count[a] < count[b];
                    int va = perm[fcs + a], vb = perm[fcs + b];
                    if (label[va] != label[vb]) return label[va] < label[vb];
                    // Use robust 2-hop signature for topological tie-breaking
                    // instead of raw node index — ensures deterministic canonical
                    // hashes regardless of input atom ordering.
                    uint64_t sigA = twoHopSignature(va, label, degree, neighbors);
                    uint64_t sigB = twoHopSignature(vb, label, degree, neighbors);
                    // Keep the signed order of the signatures while computing
                    // their hashes with defined unsigned arithmetic.
                    if (sigA != sigB)
                        return (sigA ^ (uint64_t(1) << 63)) < (sigB ^ (uint64_t(1) << 63));
                    return va < vb;
                });

                std::vector<int> tmp(cSize);
                for (int ci = 0; ci < cSize; ci++) tmp[ci] = perm[cStart + idx[ci]];
                std::copy(tmp.begin(), tmp.end(), perm.begin() + cStart);
                for (int ci = 0; ci < cSize; ci++) inv[perm[cStart + ci]] = cStart + ci;

                for (int ci = 0; ci < cSize - 1; ci++) {
                    bool wasBoundary = (cStart + ci == cEnd);
                    bool newBoundary = (count[idx[ci]] != count[idx[ci + 1]]);
                    cellEnd[cStart + ci] = newBoundary || wasBoundary;
                    if (newBoundary && !wasBoundary) changed = true;
                }
                cellEnd[cEnd] = true;
                cStart = cEnd + 1;
            }
            sStart = sEnd + 1;
        }
    }
}

inline void MolGraph::refinePartition(
        int n, std::vector<int>& perm, std::vector<bool>& cellEnd,
        const std::vector<std::vector<int>>& neighbors,
        const std::vector<int>& label, const std::vector<int>& degree, int targetPos) {
    refinePartitionImpl(n, perm, cellEnd, neighbors, label, degree, targetPos, nullptr);
}

inline int MolGraph::comparePermImpl(
        const std::vector<int>& perm1, const std::vector<int>& perm2,
        int n, const std::vector<int>& label,
        const std::vector<std::vector<int>>& neighbors, const MolGraph* graph) {
    std::vector<int> inv1(n), inv2(n);
    for (int i = 0; i < n; i++) { inv1[perm1[i]] = i; inv2[perm2[i]] = i; }
    for (int pos = 0; pos < n; pos++) {
        int v1 = perm1[pos], v2 = perm2[pos];
        int c = label[v1] - label[v2];
        if (c != 0) return c;
        if (graph) {
            c = canonicalTetrahedralTag(*graph, v1, inv1) - canonicalTetrahedralTag(*graph, v2, inv2);
            if (c != 0) return c;
        }
        const auto& nb1 = neighbors[v1];
        const auto& nb2 = neighbors[v2];
        c = static_cast<int>(nb1.size()) - static_cast<int>(nb2.size());
        if (c != 0) return c;
        std::vector<std::array<int, 3>> cn1(nb1.size()), cn2(nb2.size());
        for (size_t k = 0; k < nb1.size(); k++)
            cn1[k] = {{inv1[nb1[k]], canonicalBondLabel(graph, v1, nb1[k]),
                graph ? canonicalDoubleBondTag(*graph, v1, nb1[k], inv1) : 0}};
        for (size_t k = 0; k < nb2.size(); k++)
            cn2[k] = {{inv2[nb2[k]], canonicalBondLabel(graph, v2, nb2[k]),
                graph ? canonicalDoubleBondTag(*graph, v2, nb2[k], inv2) : 0}};
        std::sort(cn1.begin(), cn1.end());
        std::sort(cn2.begin(), cn2.end());
        for (size_t k = 0; k < cn1.size(); k++) {
            if (cn1[k] != cn2[k]) return cn1[k] < cn2[k] ? -1 : 1;
        }
    }
    return 0;
}

inline int MolGraph::comparePerm(
        const std::vector<int>& perm1, const std::vector<int>& perm2,
        int n, const std::vector<int>& label, const std::vector<std::vector<int>>& neighbors) {
    return comparePermImpl(perm1, perm2, n, label, neighbors, nullptr);
}

inline int MolGraph::partialCompare(
        const std::vector<int>& perm, const std::vector<int>& bestPerm,
        int n, const std::vector<bool>& cellEnd,
        const std::vector<int>& label,
        const std::vector<std::vector<int>>& /*neighbours*/) {
    int cellStart = 0;
    for (int pos = 0; pos < n; pos++) {
        if (cellEnd[pos]) {
            if (pos == cellStart) {
                int c = label[perm[pos]] - label[bestPerm[pos]];
                if (c != 0) return c;
            }
            cellStart = pos + 1;
        }
    }
    return 0;
}

inline MolGraph::CanonResult MolGraph::computeCanonicalLabelingImpl(
        int n, const std::vector<int>& label,
        const std::vector<int>& degree,
        const std::vector<std::vector<int>>& neighbors, const MolGraph* graph) {
    if (n == 0) return {{}, {}, {}, false};

    // Morgan-like refinement for initial ordering
    std::vector<int> mRank(label.begin(), label.end());
    std::vector<int> mNew(n);
    std::vector<std::pair<int, int>> scratch;

    for (int iter = 0; iter < 8; iter++) {
        std::unordered_set<int> s(mRank.begin(), mRank.end());
        int prevDistinct = static_cast<int>(s.size());
        for (int v = 0; v < n; v++) {
            const auto& nbs = neighbors[v];
            int d = static_cast<int>(nbs.size());
            scratch.resize(d);
            for (int k = 0; k < d; k++)
                scratch[k] = {mRank[nbs[k]], canonicalBondLabel(graph, v, nbs[k])};
            std::sort(scratch.begin(), scratch.end());
            // Keep wrapping arithmetic defined, as in computeMorganRanks.
            uint32_t h = static_cast<uint32_t>(mRank[v]) * static_cast<uint32_t>(HASH_PRIME);
            for (int k = 0; k < d; k++) {
                h = h * 31u + static_cast<uint32_t>(scratch[k].first);
                if (graph) h = h * 31u + static_cast<uint32_t>(scratch[k].second);
            }
            int32_t signedHash;
            std::memcpy(&signedHash, &h, sizeof(h));
            mNew[v] = signedHash;
        }
        mRank = mNew;
        std::unordered_set<int> ns(mRank.begin(), mRank.end());
        if (static_cast<int>(ns.size()) <= prevDistinct) break;
    }

    // Union-find for orbit computation
    std::vector<int> uf(n);
    std::iota(uf.begin(), uf.end(), 0);

    // Initial partition: sort atoms by (label, degree, mRank). Atoms equal in all three keep
    // their index order, as in the Java implementation, so the partition, and the labels beyond
    // CANON_SEARCH_LIMIT, are the same under every standard library.
    std::vector<int> sorted(n);
    std::iota(sorted.begin(), sorted.end(), 0);
    std::stable_sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        if (label[a] != label[b]) return label[a] < label[b];
        if (degree[a] != degree[b]) return degree[a] < degree[b];
        return mRank[a] < mRank[b];
    });

    std::vector<int> initPerm(sorted);
    std::vector<bool> initCellEnd(n, false);
    for (int i = 0; i < n - 1; i++) {
        int vi = sorted[i], vj = sorted[i + 1];
        initCellEnd[i] = (label[vi] != label[vj]
                       || degree[vi] != degree[vj]
                       || mRank[vi] != mRank[vj]);
    }
    initCellEnd[n - 1] = true;

    // For large molecules, skip individualisation search
    if (n > CANON_SEARCH_LIMIT) {
        std::vector<int> cl(n);
        for (int pos = 0; pos < n; pos++) cl[initPerm[pos]] = pos;
        auto orb = buildOrbits(uf, n);
        // Refine orbits using 2-hop neighbourhood signatures
        orb = refineOrbitsBy2Hop(orb, n, label, degree, neighbors);
        return { cl, orb, {}, true };
    }

    refinePartitionImpl(n, initPerm, initCellEnd, neighbors, label, degree, -1, graph);

    std::vector<int> bestPerm;
    std::vector<std::vector<int>> generators;
    static constexpr int MAX_GENERATORS = 500;
    bool generatorsTruncated = false;
    int nodeCount = 0;
    bool budgetExceeded = false;

    std::deque<std::vector<int>>  permStack;
    std::deque<std::vector<bool>> cellEndStack;

    int targetCell = findSmallestNonSingleton(n, initCellEnd);
    if (targetCell < 0) {
        bestPerm = initPerm;
    } else {
        int cellStart = 0;
        for (int i = 0; i < n; i++) {
            if (i > 0 && initCellEnd[i - 1]) cellStart = i;
            if (i == targetCell) break;
        }
        for (int i = targetCell; i >= cellStart; i--) {
            auto p = initPerm;
            auto ce = initCellEnd;
            if (i != cellStart) {
                int tmp = p[i];
                std::copy_backward(p.begin() + cellStart,
                                   p.begin() + i,
                                   p.begin() + i + 1);
                p[cellStart] = tmp;
            }
            ce[cellStart] = true;
            refinePartitionImpl(n, p, ce, neighbors, label, degree, cellStart, graph);
            permStack.push_front(std::move(p));
            cellEndStack.push_front(std::move(ce));
            nodeCount++;
        }
    }

    while (!permStack.empty() && !budgetExceeded) {
        auto perm = std::move(permStack.front()); permStack.pop_front();
        auto cellEnd = std::move(cellEndStack.front()); cellEndStack.pop_front();

        int nc = findSmallestNonSingleton(n, cellEnd);
        if (nc < 0) {
            if (bestPerm.empty()
                || comparePermImpl(perm, bestPerm, n, label, neighbors, graph) < 0) {
                bestPerm = perm;
            } else if (comparePermImpl(perm, bestPerm, n, label, neighbors, graph) == 0) {
                // Capture the automorphism generator: sigma maps
                // bestPerm[i] -> perm[i] for each position i.
                if (static_cast<int>(generators.size()) < MAX_GENERATORS) {
                    std::vector<int> gen(n);
                    for (int i = 0; i < n; i++) {
                        gen[bestPerm[i]] = perm[i];
                        detail::unionFindJoin(uf, perm[i], bestPerm[i]);
                    }
                    generators.push_back(std::move(gen));
                } else {
                    generatorsTruncated = true;
                }
            }
            continue;
        }
        if (!bestPerm.empty()
            && partialCompare(perm, bestPerm, n, cellEnd, label, neighbors) > 0)
            continue;

        int cellStart = 0;
        for (int i = 0; i < n; i++) {
            if (i > 0 && cellEnd[i - 1]) cellStart = i;
            if (i == nc) break;
        }
        int cellSize = nc - cellStart + 1;
        if (nodeCount + cellSize > MAX_SEARCH_NODES) {
            budgetExceeded = true; break;
        }

        for (int i = nc; i >= cellStart; i--) {
            auto p = perm;
            auto ce = cellEnd;
            if (i != cellStart) {
                int tmp = p[i];
                std::copy_backward(p.begin() + cellStart,
                                   p.begin() + i,
                                   p.begin() + i + 1);
                p[cellStart] = tmp;
            }
            ce[cellStart] = true;
            refinePartitionImpl(n, p, ce, neighbors, label, degree, cellStart, graph);
            permStack.push_front(std::move(p));
            cellEndStack.push_front(std::move(ce));
            nodeCount++;
        }
    }

    if (bestPerm.empty()) bestPerm = initPerm;

    std::vector<int> cl(n);
    for (int pos = 0; pos < n; pos++) cl[bestPerm[pos]] = pos;

    auto orb = buildOrbits(uf, n);
    // Refine orbits using 2-hop neighbourhood signatures
    orb = refineOrbitsBy2Hop(orb, n, label, degree, neighbors);
    return { cl, orb, std::move(generators), generatorsTruncated || budgetExceeded };
}

inline MolGraph::CanonResult MolGraph::computeCanonicalLabeling(
        int n, const std::vector<int>& label, const std::vector<int>& degree,
        const std::vector<std::vector<int>>& neighbors) {
    return computeCanonicalLabelingImpl(n, label, degree, neighbors, nullptr);
}

// ============================================================================
// Canonical hash
// ============================================================================

inline uint64_t MolGraph::computeCanonicalHash(
        const MolGraph& g, const std::vector<int>& canonLabel) {
    int n = g.n;
    auto vecAt = [](const auto& vec, int idx) -> int {
        return idx < static_cast<int>(vec.size()) ? vec[idx] : 0;
    };
    std::vector<int> inv(n);
    for (int i = 0; i < n; i++) inv[canonLabel[i]] = i;
    // Use uint64_t throughout: wrapping arithmetic is defined behaviour
    // (unlike signed int64_t overflow, which is technically UB).
    uint64_t hash = 0;
    for (int cp = 0; cp < n; cp++) {
        int v = inv[cp];
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(g.label[v]);
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(static_cast<int64_t>(vecAt(g.formalCharge, v)) + 512);
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(vecAt(g.massNumber, v));
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(vecAt(g.hydrogenCount, v));
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(vecAt(g.atomClass, v));
        // include tetrahedral chirality — enantiomers must hash differently
        hash = hash * static_cast<uint64_t>(HASH_PRIME)
             + static_cast<uint64_t>(canonicalTetrahedralTag(g, v, canonLabel));
    }
    for (int cp = 0; cp < n; cp++) {
        int v = inv[cp];
        const auto& nbs = g.neighbors[v];
        std::vector<uint64_t> edgeSigs;
        edgeSigs.reserve(nbs.size());
        for (int nb : nbs) {
            int cn = canonLabel[nb];
            uint64_t edgeKey = static_cast<uint64_t>(std::min(cp, cn)) * static_cast<uint64_t>(n)
                             + static_cast<uint64_t>(std::max(cp, cn));
            uint64_t sig = edgeKey;
            sig = sig * static_cast<uint64_t>(HASH_PRIME)
                + static_cast<uint64_t>(g.bondOrder(v, nb));
            sig = sig * 2ULL + static_cast<uint64_t>(g.bondInRing(v, nb) ? 1 : 0);
            sig = sig * 2ULL + static_cast<uint64_t>(g.bondAromatic(v, nb) ? 1 : 0);
            // include E/Z double-bond stereo (use canonical-position guard to avoid double-counting)
            if (cp < cn) {
                sig = sig * static_cast<uint64_t>(HASH_PRIME)
                    + static_cast<uint64_t>(canonicalDoubleBondTag(g, v, nb, canonLabel));
            }
            edgeSigs.push_back(sig);
        }
        std::sort(edgeSigs.begin(), edgeSigs.end());
        for (uint64_t sig : edgeSigs) {
            hash = hash * static_cast<uint64_t>(HASH_PRIME)
                 + sig;
        }
    }
    return hash;
}

// ============================================================================
// Lazy canonical properties
// ============================================================================

inline void MolGraph::ensureCanonical() const {
    if (lazyCaches_.ready(LazyCaches::Canonical)) return;
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (canonicalComputed_) { lazyCaches_.publish(LazyCaches::Canonical); return; }
    morganRank = computeMorganRanks(n, label, neighbors);
    const auto invariantLabels = canonicalAtomLabels();
    auto clResult = computeCanonicalLabelingImpl(n, invariantLabels, degree, neighbors, this);
    canonicalLabel = std::move(clResult.canonLabel);
    orbit = std::move(clResult.orbit);
    autGenerators_ = std::move(clResult.autGenerators);
    autGeneratorsTruncated_ = clResult.generatorsTruncated;
    auto preservesProperties = [&](const std::vector<int>& permutation) {
        if (permutation.size() != static_cast<size_t>(n)) return false;
        std::vector<uint8_t> used(n, 0);
        auto equal = [](const auto& values, int a, int b) {
            return (a < static_cast<int>(values.size()) ? static_cast<int>(values[a]) : 0)
                == (b < static_cast<int>(values.size()) ? static_cast<int>(values[b]) : 0);
        };
        for (int atom = 0; atom < n; ++atom) {
            int image = permutation[atom];
            if (image < 0 || image >= n || used[image]++) return false;
            if (!equal(atomicNum, atom, image) || !equal(formalCharge, atom, image)
                || !equal(massNumber, atom, image) || !equal(hydrogenCount, atom, image)
                || !equal(atomClass, atom, image) || !equal(ring, atom, image)
                || !equal(aromatic, atom, image)
                || canonicalTetrahedralTag(*this, atom, permutation) != canonicalTetrahedralTag(*this, image, {})) return false;
            for (int neighbor : neighbors[atom]) {
                if (neighbor <= atom) continue;
                int other = permutation[neighbor];
                if (other < 0 || other >= n
                    || bondOrder(atom, neighbor) != bondOrder(image, other)
                    || bondInRing(atom, neighbor) != bondInRing(image, other)
                    || bondAromatic(atom, neighbor) != bondAromatic(image, other)
                    || canonicalDoubleBondTag(*this, atom, neighbor, permutation)
                        != canonicalDoubleBondTag(*this, image, other, {})) return false;
            }
        }
        return true;
    };
    autGenerators_.erase(std::remove_if(autGenerators_.begin(), autGenerators_.end(), [&](const auto& generator) {
        if (preservesProperties(generator)) return false;
        autGeneratorsTruncated_ = true;
        return true;
    }), autGenerators_.end());
    std::vector<int> verifiedOrbits(n);
    std::iota(verifiedOrbits.begin(), verifiedOrbits.end(), 0);
    for (const auto& generator : autGenerators_)
        for (int atom = 0; atom < n; ++atom) detail::unionFindJoin(verifiedOrbits, atom, generator[atom]);
    orbit = buildOrbits(verifiedOrbits, n);
    canonicalHash = computeCanonicalHash(*this, canonicalLabel);
    canonicalComputed_ = true;
    lazyCaches_.publish(LazyCaches::Canonical);
}

} // namespace smsd

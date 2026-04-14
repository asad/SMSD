/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Topological Torsion Fingerprint (Nilakantan et al. 1987).
 * Enumerates all 4-atom linear paths A-B-C-D, hashes atom types
 * with canonical ordering (min of forward/reverse hash).
 */
#pragma once
#ifndef FP_MOL_TORSION_HPP
#define FP_MOL_TORSION_HPP

#include "fp/common.hpp"
#include "smsd/mol_graph.hpp"

#include <cstdint>
#include <vector>

namespace smsd {
namespace fp {
namespace mol {

/// Compute topological torsion atom type: (atomicNum, numHeavyNeighbors,
/// numPiElectrons, isRing) as a single 32-bit integer.
inline int ttAtomType(const MolGraph& g, int idx) {
    int z = g.atomicNum[idx];
    int heavyDeg = g.degree[idx];
    int nPi = 0;
    for (int nb : g.neighbors[idx]) {
        int bo = g.bondOrder(idx, nb);
        if (bo == 2) nPi += 1;
        else if (bo == 3) nPi += 2;
    }
    int inRing = g.ring[idx] ? 1 : 0;
    return (z << 16) | (heavyDeg << 8) | (nPi << 4) | inRing;
}

/// Topological torsion fingerprint (binary).
/// @param fpSize  Fingerprint size in bits
inline std::vector<uint64_t> computeTopologicalTorsion(
    const MolGraph& mol, int fpSize)
{
    auto fp = makeBitVector(fpSize);
    int n = mol.n;
    if (n < 4) return fp;

    std::vector<int> atomType(n);
    for (int i = 0; i < n; ++i)
        atomType[i] = ttAtomType(mol, i);

    for (int a = 0; a < n; ++a) {
        for (int b : mol.neighbors[a]) {
            for (int c : mol.neighbors[b]) {
                if (c == a) continue;
                for (int d : mol.neighbors[c]) {
                    if (d == b || d == a) continue;
                    uint64_t fwd = FNV1A_SEED;
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[a]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[b]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[c]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[d]));
                    uint64_t rev = FNV1A_SEED;
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[d]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[c]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[b]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[a]));
                    setBit(fp, fpSize, (fwd <= rev) ? fwd : rev);
                }
            }
        }
    }
    return fp;
}

/// Topological torsion count fingerprint.
/// @param fpSize  Number of counter bins
inline std::vector<int> computeTopologicalTorsionCounts(
    const MolGraph& mol, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    int n = mol.n;
    if (n < 4) return counts;

    std::vector<int> atomType(n);
    for (int i = 0; i < n; ++i)
        atomType[i] = ttAtomType(mol, i);

    for (int a = 0; a < n; ++a) {
        for (int b : mol.neighbors[a]) {
            for (int c : mol.neighbors[b]) {
                if (c == a) continue;
                for (int d : mol.neighbors[c]) {
                    if (d == b || d == a) continue;
                    if (a > d) continue; // canonical: count each path once
                    uint64_t fwd = FNV1A_SEED;
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[a]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[b]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[c]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[d]));
                    uint64_t rev = FNV1A_SEED;
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[d]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[c]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[b]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[a]));
                    uint64_t h = (fwd <= rev) ? fwd : rev;
                    counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
                }
            }
        }
    }
    return counts;
}

} // namespace mol
} // namespace fp
} // namespace smsd

#endif // FP_MOL_TORSION_HPP

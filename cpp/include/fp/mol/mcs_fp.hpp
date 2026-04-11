/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * MCS-aware path fingerprint.
 * Encodes element, ring, aromatic, tautomer class, degree per atom;
 * bond order, ring, aromatic per bond. Paths are hashed canonically
 * (direction-independent). Compatible with the Java SMSD MCS fingerprint.
 */
#pragma once
#ifndef FP_MOL_MCS_FP_HPP
#define FP_MOL_MCS_FP_HPP

#include "fp/common.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace smsd {
namespace fp {
namespace mol {

namespace mcs_fp_detail {

inline uint32_t atomHash(const MolGraph& g, bool tautAware, int atom) {
    uint32_t h = 17;
    uint32_t label = static_cast<uint32_t>(g.atomicNum[atom]);
    bool hasTautClass = !g.tautomerClass.empty() && g.tautomerClass[atom] >= 0;
    if (tautAware && hasTautClass) label = 999;
    h = h * 37 + label;
    h = h * 37 + (g.ring[atom] ? 1u : 0u);
    h = h * 37 + (g.aromatic[atom] ? 1u : 0u);
    int tc = (!g.tautomerClass.empty()) ? g.tautomerClass[atom] : -1;
    h = h * 37 + static_cast<uint32_t>(tautAware ? 0 : tc);
    h = h * 37 + static_cast<uint32_t>(g.degree[atom]);
    return h;
}

inline uint32_t bondHash(const MolGraph& g, int from, int to) {
    uint32_t h = 17;
    h = h * 37 + static_cast<uint32_t>(g.bondOrder(from, to));
    h = h * 37 + (g.bondInRing(from, to) ? 1u : 0u);
    h = h * 37 + (g.bondAromatic(from, to) ? 1u : 0u);
    return h;
}

inline int hashPath(const MolGraph& g, bool tautAware,
                    const int* path, int len) {
    uint32_t fwd = 17, rev = 17;
    for (int i = 0; i < len; ++i) {
        fwd = fwd * 31 + atomHash(g, tautAware, path[i]);
        rev = rev * 31 + atomHash(g, tautAware, path[len - 1 - i]);
        if (i < len - 1) {
            fwd = fwd * 31 + bondHash(g, path[i], path[i + 1]);
            rev = rev * 31 + bondHash(g, path[len - 1 - i], path[len - 2 - i]);
        }
    }
    return static_cast<int>(std::min(fwd & 0x7FFFFFFFu, rev & 0x7FFFFFFFu));
}

inline void setBitMCS(std::vector<uint64_t>& fp, int hash, int fpSize) {
    int bit = hash % fpSize;
    fp[bit / 64] |= 1ULL << (bit % 64);
}

inline void enumeratePaths(const MolGraph& g, bool tautAware,
                           const std::vector<uint8_t>& isHeavy,
                           std::vector<uint64_t>& fp, int fpSize,
                           std::vector<uint8_t>& visited, int* path,
                           int depth, int maxDepth) {
    int cur = path[depth - 1];
    for (int nb : g.neighbors[cur]) {
        if (visited[nb] || !isHeavy[nb]) continue;
        path[depth] = nb;
        visited[nb] = true;
        setBitMCS(fp, hashPath(g, tautAware, path, depth + 1), fpSize);
        if (depth < maxDepth)
            enumeratePaths(g, tautAware, isHeavy, fp, fpSize, visited,
                           path, depth + 1, maxDepth);
        visited[nb] = false;
    }
}

} // namespace mcs_fp_detail

/// MCS-aware path fingerprint.
/// @param pathLength  Maximum path length in bonds (default 7)
/// @param fpSize      Fingerprint size in bits (default 2048)
inline std::vector<uint64_t> computeMCSFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(numWords, 0ULL);
    int n = mol.n;
    if (n == 0) return fp;

    std::vector<uint8_t> isHeavy(n);
    int heavyCount = 0;
    for (int i = 0; i < n; ++i) {
        isHeavy[i] = (mol.atomicNum[i] != 1) ? 1 : 0;
        if (isHeavy[i]) ++heavyCount;
    }
    if (heavyCount == 0) return fp;

    bool tautAware = !mol.tautomerClass.empty();
    std::vector<int> pathBuf(pathLength + 1);
    std::vector<uint8_t> visited(n, 0);
    for (int start = 0; start < n; ++start) {
        if (!isHeavy[start]) continue;
        pathBuf[0] = start;
        visited[start] = true;
        mcs_fp_detail::setBitMCS(fp, mcs_fp_detail::hashPath(mol, tautAware, pathBuf.data(), 1), fpSize);
        mcs_fp_detail::enumeratePaths(mol, tautAware, isHeavy, fp, fpSize, visited,
                                      pathBuf.data(), 1, pathLength);
        visited[start] = false;
    }
    return fp;
}

} // namespace mol
} // namespace fp
} // namespace smsd

#endif // FP_MOL_MCS_FP_HPP

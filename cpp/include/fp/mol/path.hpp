/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Path-based molecular fingerprint.
 * Enumerates all simple paths up to a given length and hashes them
 * using FNV-1a 64-bit into a bit-vector.
 */
#pragma once
#ifndef FP_MOL_PATH_HPP
#define FP_MOL_PATH_HPP

#include "fp/common.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace smsd {
namespace fp {
namespace mol {

/// Simple path-based fingerprint generation for a single molecule.
/// Enumerates all paths up to pathLength and hashes them into fpSize bits.
inline std::vector<uint64_t> computePathFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    auto fp = makeBitVector(fpSize);
    if (mol.n == 0) return fp;
    pathLength = std::min(pathLength, 7);

    struct Frame {
        int atom;
        int depth;
        int nbIdx;
        int path[8];
        int bonds[7];
    };

    std::vector<uint8_t> visited(mol.n, 0);
    Frame stack[16];

    for (int start = 0; start < mol.n; ++start) {
        {
            int p[1] = { mol.label[start] };
            setBit(fp, fpSize, fnv1a(p, 1));
        }

        int top = 0;
        stack[0].atom = start;
        stack[0].depth = 0;
        stack[0].nbIdx = 0;
        stack[0].path[0] = mol.label[start];
        visited[start] = true;

        while (top >= 0) {
            Frame& f = stack[top];
            bool expanded = false;

            if (f.depth < pathLength) {
                const auto& nbs = mol.neighbors[f.atom];
                int nbCount = static_cast<int>(nbs.size());
                for (; f.nbIdx < nbCount; ++f.nbIdx) {
                    int nb = nbs[f.nbIdx];
                    if (visited[nb]) continue;

                    int newDepth = f.depth + 1;
                    int bondOrd = mol.bondOrder(f.atom, nb);

                    int pathData[16];
                    int pathLen = 0;
                    for (int d = 0; d <= f.depth; ++d) {
                        pathData[pathLen++] = f.path[d];
                        if (d < f.depth)
                            pathData[pathLen++] = f.bonds[d] + 1000;
                    }
                    pathData[pathLen++] = bondOrd + 1000;
                    pathData[pathLen++] = mol.label[nb];

                    uint64_t hFwd = fnv1a(pathData, pathLen);
                    int revData[16];
                    for (int r = 0; r < pathLen; ++r)
                        revData[r] = pathData[pathLen - 1 - r];
                    uint64_t hRev = fnv1a(revData, pathLen);
                    setBit(fp, fpSize, std::min(hFwd, hRev));

                    if (newDepth < pathLength && top + 1 < 15) {
                        visited[nb] = true;
                        f.nbIdx++;
                        top++;
                        stack[top].atom = nb;
                        stack[top].depth = newDepth;
                        stack[top].nbIdx = 0;
                        for (int d = 0; d <= stack[top-1].depth; ++d)
                            stack[top].path[d] = stack[top-1].path[d];
                        for (int d = 0; d < stack[top-1].depth; ++d)
                            stack[top].bonds[d] = stack[top-1].bonds[d];
                        stack[top].bonds[stack[top-1].depth] = bondOrd;
                        stack[top].path[newDepth] = mol.label[nb];
                        expanded = true;
                        break;
                    }
                }
            }

            if (!expanded) {
                visited[stack[top].atom] = false;
                top--;
            }
        }
        visited[start] = false;
    }
    return fp;
}

} // namespace mol
} // namespace fp
} // namespace smsd

#endif // FP_MOL_PATH_HPP

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Circular (Morgan/ECFP/FCFP) molecular fingerprints.
 * Rogers & Hahn (2010) with SMSD extensions for tautomer-invariance.
 *
 * Binary and count-based variants for both ECFP and FCFP modes.
 */
#pragma once
#ifndef FP_MOL_CIRCULAR_HPP
#define FP_MOL_CIRCULAR_HPP

#include "fp/common.hpp"
#include "fp/mol/pharmacophore.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <vector>

namespace smsd {
namespace fp {
namespace mol {

// ── Valence helpers (OpenSMILES multi-valence model) ─────────────────────────

/// Returns the list of allowed valences for common organic-subset elements.
inline std::pair<const int*, int> defaultValences(int z) {
    static const int v_H[]  = {1};
    static const int v_B[]  = {3};
    static const int v_C[]  = {4};
    static const int v_N[]  = {3, 5};
    static const int v_O[]  = {2};
    static const int v_F[]  = {1};
    static const int v_Al[] = {3};
    static const int v_Si[] = {4};
    static const int v_P[]  = {3, 5};
    static const int v_S[]  = {2, 4, 6};
    static const int v_Cl[] = {1};
    static const int v_As[] = {3, 5};
    static const int v_Se[] = {2, 4, 6};
    static const int v_Br[] = {1};
    static const int v_Sb[] = {3, 5};
    static const int v_Te[] = {2, 4, 6};
    static const int v_I[]  = {1, 3, 5, 7};
    #define VSPAN(arr) {arr, sizeof(arr)/sizeof(arr[0])}
    switch (z) {
        case 1:  return VSPAN(v_H);
        case 5:  return VSPAN(v_B);
        case 6:  return VSPAN(v_C);
        case 7:  return VSPAN(v_N);
        case 8:  return VSPAN(v_O);
        case 9:  return VSPAN(v_F);
        case 13: return VSPAN(v_Al);
        case 14: return VSPAN(v_Si);
        case 15: return VSPAN(v_P);
        case 16: return VSPAN(v_S);
        case 17: return VSPAN(v_Cl);
        case 33: return VSPAN(v_As);
        case 34: return VSPAN(v_Se);
        case 35: return VSPAN(v_Br);
        case 51: return VSPAN(v_Sb);
        case 52: return VSPAN(v_Te);
        case 53: return VSPAN(v_I);
        default: return {nullptr, 0};
    }
    #undef VSPAN
}

/// Compute implicit hydrogen count using OpenSMILES multi-valence model.
inline int implicitH(int atomicNum, int bondOrderSum, int formalCharge) {
    auto [vals, count] = defaultValences(atomicNum);
    if (!vals) return 0;
    int absCharge = std::abs(formalCharge);
    for (int i = 0; i < count; i++) {
        int capacity = vals[i] - absCharge;
        if (capacity >= bondOrderSum) return capacity - bondOrderSum;
    }
    return 0;
}

// ── ECFP (Extended Connectivity Fingerprint) ─────────────────────────────────

/// ECFP binary fingerprint.
/// @param radius  Maximum radius (-1 = unlimited, 2 = ECFP4, 3 = ECFP6)
/// @param fpSize  Fingerprint size in bits (default 2048)
inline std::vector<uint64_t> computeCircularFingerprintECFP(
    const MolGraph& mol, int radius, int fpSize)
{
    auto fp = makeBitVector(fpSize);
    if (mol.n == 0) return fp;

    int n = mol.n;
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, mol.atomicNum[i]);
        h = fnvMix(h, mol.degree[i]);
        int bondOrdSum = 0;
        for (int nb : mol.neighbors[i]) bondOrdSum += mol.bondOrder(i, nb);
        h = fnvMix(h, static_cast<uint64_t>(bondOrdSum));
        int hCount = implicitH(mol.atomicNum[i], bondOrdSum, mol.formalCharge[i]);
        h = fnvMix(h, static_cast<uint64_t>(hCount));
        h = fnvMix(h, mol.ring[i] ? 1 : 0);
        h = fnvMix(h, mol.aromatic[i] ? 1 : 0);
        h = fnvMix(h, static_cast<uint64_t>(mol.formalCharge[i] + 4));
        if (!mol.massNumber.empty() && mol.massNumber[i] > 0)
            h = fnvMix(h, static_cast<uint64_t>(mol.massNumber[i]));
        if (!mol.tautomerClass.empty() && mol.tautomerClass[i] >= 0)
            h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(mol.tautomerClass[i]));
        atomHash[i] = h;
        setBit(fp, fpSize, h);
    }

    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
    std::unordered_set<uint64_t> seenHashes;
    if (radius < 0) {
        for (int i = 0; i < n; ++i) seenHashes.insert(atomHash[i]);
    }

    for (int r = 1; r <= maxRadius; ++r) {
        bool anyNew = false;
        for (int i = 0; i < n; ++i) {
            int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                int nb = mol.neighbors[i][k];
                nbHashes[k] = prevHash[nb] ^ (static_cast<uint64_t>(mol.bondOrder(i, nb)) * 31ULL);
            }
            std::sort(nbHashes.begin(), nbHashes.end());

            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[i]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            nextHash[i] = h;

            setBit(fp, fpSize, h);
            if (radius < 0 && seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return fp;
}

/// ECFP count fingerprint.
/// @param radius  Maximum radius (-1 = unlimited, 2 = ECFP4, 3 = ECFP6)
/// @param fpSize  Number of counter bins (default 2048)
inline std::vector<int> computeCircularFingerprintECFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    if (mol.n == 0) return counts;

    int n = mol.n;
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, mol.atomicNum[i]);
        h = fnvMix(h, mol.degree[i]);
        int bondOrdSum = 0;
        for (int nb : mol.neighbors[i]) bondOrdSum += mol.bondOrder(i, nb);
        h = fnvMix(h, static_cast<uint64_t>(bondOrdSum));
        int hCount = implicitH(mol.atomicNum[i], bondOrdSum, mol.formalCharge[i]);
        h = fnvMix(h, static_cast<uint64_t>(hCount));
        h = fnvMix(h, mol.ring[i] ? 1 : 0);
        h = fnvMix(h, mol.aromatic[i] ? 1 : 0);
        h = fnvMix(h, static_cast<uint64_t>(mol.formalCharge[i] + 4));
        if (!mol.massNumber.empty() && mol.massNumber[i] > 0)
            h = fnvMix(h, static_cast<uint64_t>(mol.massNumber[i]));
        if (!mol.tautomerClass.empty() && mol.tautomerClass[i] >= 0)
            h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(mol.tautomerClass[i]));
        atomHash[i] = h;
        counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
    }

    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
    std::unordered_set<uint64_t> seenHashes;

    for (int r = 1; r <= maxRadius; ++r) {
        bool anyNew = false;
        for (int i = 0; i < n; ++i) {
            int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                int nb = mol.neighbors[i][k];
                nbHashes[k] = prevHash[nb] ^ (static_cast<uint64_t>(mol.bondOrder(i, nb)) * 31ULL);
            }
            std::sort(nbHashes.begin(), nbHashes.end());

            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[i]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            nextHash[i] = h;

            counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
            if (seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return counts;
}

// ── FCFP (Feature-Class Fingerprint) ─────────────────────────────────────────

/// Fingerprint invariant mode.
enum class FingerprintMode { ECFP, FCFP };

/// FCFP binary fingerprint using pharmacophoric feature classes.
inline std::vector<uint64_t> computeCircularFingerprintFCFP(
    const MolGraph& mol, int radius, int fpSize)
{
    auto fp = makeBitVector(fpSize);
    if (mol.n == 0) return fp;

    int n = mol.n;
    const auto& pharmFeatures = mol.getPharmacophoreFeatures();
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, static_cast<uint64_t>(pharmFeatures[i]));
        h = fnvMix(h, static_cast<uint64_t>(mol.degree[i]));
        h = fnvMix(h, mol.ring[i] ? 1ULL : 0ULL);
        atomHash[i] = h;
        setBit(fp, fpSize, h);
    }
    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
    std::unordered_set<uint64_t> seenHashes;
    if (radius < 0) {
        for (int i = 0; i < n; ++i) seenHashes.insert(atomHash[i]);
    }
    for (int r = 1; r <= maxRadius; ++r) {
        bool anyNew = false;
        for (int i = 0; i < n; ++i) {
            int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                int nb = mol.neighbors[i][k];
                nbHashes[k] = prevHash[nb] ^ (static_cast<uint64_t>(mol.bondOrder(i, nb)) * 31ULL);
            }
            std::sort(nbHashes.begin(), nbHashes.end());
            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[i]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            nextHash[i] = h;
            setBit(fp, fpSize, h);
            if (radius < 0 && seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return fp;
}

/// FCFP count fingerprint.
inline std::vector<int> computeCircularFingerprintFCFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    if (mol.n == 0) return counts;

    int n = mol.n;
    const auto& pharmFeatures = mol.getPharmacophoreFeatures();
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, static_cast<uint64_t>(pharmFeatures[i]));
        h = fnvMix(h, static_cast<uint64_t>(mol.degree[i]));
        h = fnvMix(h, mol.ring[i] ? 1ULL : 0ULL);
        atomHash[i] = h;
        counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
    }
    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
    std::unordered_set<uint64_t> seenHashes;

    for (int r = 1; r <= maxRadius; ++r) {
        bool anyNew = false;
        for (int i = 0; i < n; ++i) {
            int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                int nb = mol.neighbors[i][k];
                nbHashes[k] = prevHash[nb] ^ (static_cast<uint64_t>(mol.bondOrder(i, nb)) * 31ULL);
            }
            std::sort(nbHashes.begin(), nbHashes.end());
            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[i]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            nextHash[i] = h;
            counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
            if (seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return counts;
}

/// All-in-one convenience: parse a SMILES string and compute its ECFP
/// fingerprint as a bit-packed uint64_t vector.
inline std::vector<uint64_t> fingerprintFromSmiles(const std::string& smi,
                                                    int radius = 2,
                                                    int fpSize = 2048) {
    MolGraph g = smsd::parseSMILES(smi);
    return computeCircularFingerprintECFP(g, radius, fpSize);
}

} // namespace mol
} // namespace fp
} // namespace smsd

#endif // FP_MOL_CIRCULAR_HPP

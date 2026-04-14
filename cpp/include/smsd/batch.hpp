/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only OpenMP batch processing for SMSD.
 *
 * Provides multi-CPU parallel batch operations:
 *   - batchSubstructure: 1 query vs N targets, returns hit mask
 *   - batchFindSubstructure: 1 query vs N targets, returns atom mappings
 *   - batchMCS: 1 query vs N targets, returns all MCS mappings
 *   - screenAndMatch: RASCAL pre-screen + exact MCS on hits
 *   - batchFingerprint: parallel fingerprint generation
 *   - TargetCorpus: pre-warmed target collection for repeated queries
 *
 * Falls back to sequential execution when OpenMP is unavailable.
 * Thread-safe: each thread uses its own scratch buffers.
 *
 * Compile with: -fopenmp (gcc/clang) or /openmp (MSVC)
 * On macOS with clang: brew install libomp, then -Xpreprocessor -fopenmp -lomp
 */
#pragma once
#ifndef SMSD_BATCH_HPP
#define SMSD_BATCH_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/vf2pp.hpp"
#include "smsd/mcs.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <stdexcept>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace smsd {
namespace batch {

/// FNV-1a 64-bit offset basis.
constexpr uint64_t FNV1A_SEED  = 0xCBF29CE484222325ULL;
/// FNV-1a 64-bit prime.
constexpr uint64_t FNV1A_PRIME = 0x100000001B3ULL;

// ============================================================================
// Internal: thread count resolution
// ============================================================================
namespace detail {

/// Resolve thread count: 0 = auto (use all available), else clamp to [1, max].
inline int resolveThreads(int numThreads) {
#ifdef _OPENMP
    if (numThreads <= 0) {
        return omp_get_max_threads();
    }
    return std::max(1, std::min(numThreads, omp_get_max_threads()));
#else
    (void)numThreads;
    return 1;
#endif
}

/// Precompute common lazy MolGraph invariants from a single thread before
/// entering parallel matching code.
inline void prewarmGraph(const MolGraph& g) {
    g.ensureCanonical();
    g.ensureRingCounts();
    g.ensurePatternFP();
    g.getPharmacophoreFeatures();
    g.getNLF1();
    g.getNLF2();
    g.getNLF3();
    g.getNeighborsByDegDesc();
}

/// Simple path-based fingerprint generation for a single molecule.
/// Enumerates all paths up to pathLength and hashes them into fpSize bits.
/// Returns a vector of uint64_t words representing the fingerprint.
inline std::vector<uint64_t> computePathFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(numWords, 0ULL);

    if (mol.n == 0) return fp;
    pathLength = std::min(pathLength, 7); // clamp to path[8] buffer size

    // Hash helper: FNV-1a 64-bit
    auto fnv1a = [](const int* data, int len) -> uint64_t {
        uint64_t h = FNV1A_SEED;
        for (int i = 0; i < len; ++i) {
            h ^= static_cast<uint64_t>(static_cast<uint32_t>(data[i]));
            h *= FNV1A_PRIME;
        }
        return h;
    };

    // DFS-based path enumeration from each atom
    // Stack-based to avoid recursion overhead
    struct Frame {
        int atom;
        int depth;
        int nbIdx;   // neighbor iteration index — resume after backtrack
        int path[8]; // atom labels, max path length + 1
        int bonds[7]; // bond orders between consecutive path atoms
    };

    std::vector<uint8_t> visited(mol.n, 0);
    Frame stack[16]; // max DFS depth = pathLength (clamped to 7) + 1

    for (int start = 0; start < mol.n; ++start) {
        // Single atom path
        {
            int p[1] = { mol.label[start] };
            uint64_t h = fnv1a(p, 1);
            int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
            fp[bit / 64] |= (1ULL << (bit % 64));
        }

        // Multi-atom paths via iterative DFS
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

                    // Build path including all bond orders interleaved
                    int newDepth = f.depth + 1;
                    int bondOrd = mol.bondOrder(f.atom, nb);

                    // path = [label0, bond01, label1, bond12, label2, ...]
                    int pathData[16]; // max 2*pathLength + 1
                    int pathLen = 0;
                    for (int d = 0; d <= f.depth; ++d) {
                        pathData[pathLen++] = f.path[d];
                        if (d < f.depth) {
                            pathData[pathLen++] = f.bonds[d] + 1000;
                        }
                    }
                    pathData[pathLen++] = bondOrd + 1000;
                    pathData[pathLen++] = mol.label[nb];

                    uint64_t hFwd = fnv1a(pathData, pathLen);

                    // Canonical hash: use min(fwd, rev) so A→B and B→A set the
                    // same single bit, preserving correct bit density.
                    int revData[16];
                    for (int r = 0; r < pathLen; ++r)
                        revData[r] = pathData[pathLen - 1 - r];
                    uint64_t hRev = fnv1a(revData, pathLen);

                    uint64_t h = std::min(hFwd, hRev);
                    int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
                    fp[bit / 64] |= (1ULL << (bit % 64));

                    // Push to stack if we can go deeper
                    if (newDepth < pathLength && top + 1 < 15) {
                        visited[nb] = true;
                        f.nbIdx++; // advance BEFORE push so we resume at next neighbor on return
                        top++;
                        stack[top].atom = nb;
                        stack[top].depth = newDepth;
                        stack[top].nbIdx = 0; // fresh start for new depth
                        for (int d = 0; d <= stack[top-1].depth; ++d)
                            stack[top].path[d] = stack[top-1].path[d];
                        for (int d = 0; d < stack[top-1].depth; ++d)
                            stack[top].bonds[d] = stack[top-1].bonds[d];
                        stack[top].bonds[stack[top-1].depth] = bondOrd;
                        stack[top].path[newDepth] = mol.label[nb];
                        expanded = true;
                        break; // DFS: go deeper first
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

/// Returns the list of allowed valences for common organic-subset elements,
/// following the OpenSMILES / Daylight convention.  Multi-valent elements
/// (N, P, S, As, Se, Sb, Te, I) list every standard valence in ascending
/// order so that implicitH() can pick the smallest fitting one.
/// Returns an empty span for elements with no known default valence.
inline std::pair<const int*, int> defaultValences(int z) {
    // Static tables -- one per element, ascending order
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

/// Computes implicit hydrogen count for an atom using the OpenSMILES
/// multi-valence model: pick the smallest allowed valence v such that
/// (v - |charge|) >= bondOrderSum, then implicitH = (v - |charge|) - bondOrderSum.
/// Returns 0 for elements with no known default valence or when no
/// valence can accommodate the observed bond order sum.
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

/// Circular (Morgan/ECFP-style) fingerprint for a single molecule.
/// Each atom's environment is hashed at increasing radii using canonical
/// invariants (atomicNum, degree, implicitH, ring, aromatic, charge, tautomerClass).
/// Neighbour hashes are sorted for automorphism invariance.
///
/// @param radius  Maximum radius (-1 = unlimited / whole molecule, 2 = ECFP4, 3 = ECFP6)
/// @param fpSize  Fingerprint size in bits (default 2048)
inline std::vector<uint64_t> computeCircularFingerprintECFP(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(numWords, 0ULL);
    if (mol.n == 0) return fp;

    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };

    int n = mol.n;
    // Initial atom invariants (Rogers & Hahn 2010, Table 1)
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, mol.atomicNum[i]);                                  // R&H #1: atomic number
        h = fnvMix(h, mol.degree[i]);                                     // R&H #2: heavy atom degree
        int bondOrdSum = 0;
        for (int nb : mol.neighbors[i]) bondOrdSum += mol.bondOrder(i, nb);
        h = fnvMix(h, static_cast<uint64_t>(bondOrdSum));                 // R&H #3: bond order sum (valence)
        int massNum = (static_cast<int>(mol.massNumber.size()) > i)
                      ? mol.massNumber[i] : 0;
        h = fnvMix(h, static_cast<uint64_t>(massNum));                    // R&H #4: atomic mass number
        h = fnvMix(h, static_cast<uint64_t>(mol.formalCharge[i] + 4));    // R&H #5: formal charge
        int hCount = implicitH(mol.atomicNum[i], bondOrdSum, mol.formalCharge[i]);
        h = fnvMix(h, static_cast<uint64_t>(hCount));                     // R&H #6: attached H count
        h = fnvMix(h, mol.ring[i] ? 1 : 0);                              // R&H #7: ring membership
        h = fnvMix(h, mol.aromatic[i] ? 1 : 0);                          // Daylight extension: aromaticity
        if (!mol.tautomerClass.empty() && mol.tautomerClass[i] >= 0)
            h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(mol.tautomerClass[i]));
        atomHash[i] = h;
        int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
        fp[bit / 64] |= (1ULL << (bit % 64));
    }

    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;

    for (int r = 1; r <= maxRadius; ++r) {
        bool anyNew = false;
        for (int i = 0; i < n; ++i) {
            int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; } // isolated atom — no expansion
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

            int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
            uint64_t mask = 1ULL << (bit % 64);
            if (!(fp[bit / 64] & mask)) anyNew = true;
            fp[bit / 64] |= mask;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return fp;
}

/// Fingerprint invariant mode.
enum class FingerprintMode { ECFP, FCFP };

/// Classify atom into pharmacophoric feature classes (FCFP invariant).
/// Returns bitmask: bit 0=donor, 1=acceptor, 2=positive, 3=negative, 4=aromatic, 5=hydrophobic.
inline int classifyPharmacophore(const MolGraph& g, int idx) {
    int z = g.atomicNum[idx];
    int charge = g.formalCharge[idx];
    bool arom = g.aromatic[idx];
    int deg = g.degree[idx];
    int features = 0;

    // H-bond donor: N-H, O-H, S-H (Rogers & Hahn 2010)
    if (z == 7 || z == 8 || z == 16) {
        int typVal = (z == 7) ? 3 : 2; // S typical valence 2
        if (deg < typVal && charge >= 0) features |= 1;
    }

    // H-bond acceptor: N (not pyrrole-type), O, F, S (not thiophene-type)
    // Pyrrole N: aromatic + has H (lone pair donated to pi) — NOT acceptor
    // Pyridine N: aromatic + no H (lone pair in plane, available) — IS acceptor
    bool isPyrroleTypeN = false;
    if (z == 7 && arom) {
        isPyrroleTypeN = (g.hydrogenCount[idx] > 0);
    }
    bool isAcceptorN = (z == 7 && charge <= 0 && !isPyrroleTypeN);
    bool isAcceptorS = (z == 16 && !arom);
    if (isAcceptorN || z == 8 || z == 9 || isAcceptorS) features |= 2;

    // Positive ionisable: charged N, or basic amine (exclude amide N, aniline N)
    if (z == 7 && charge > 0) features |= 4;
    if (z == 7 && !arom && deg <= 3 && charge == 0) {
        bool isAmide = false;
        bool isAniline = false;
        for (int nb : g.neighbors[idx]) {
            if (g.atomicNum[nb] == 6) {
                if (g.aromatic[nb]) { isAniline = true; break; }
                for (int nb2 : g.neighbors[nb]) {
                    if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2) {
                        isAmide = true; break;
                    }
                }
            }
            if (isAmide) break;
        }
        if (!isAmide && !isAniline) features |= 4;
    }

    // Negative ionisable: carboxylic acid (require C=O pattern), phosphoric, sulfonic
    if (z == 8 && charge < 0) features |= 8;
    if (z == 8 && charge == 0 && deg == 1) {
        for (int nb : g.neighbors[idx]) {
            int nbZ = g.atomicNum[nb];
            if (nbZ == 6 && g.bondOrder(idx, nb) == 1) { // this O is single-bonded (OH)
                for (int nb2 : g.neighbors[nb])
                    if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2)
                        features |= 8;
            }
            if (nbZ == 15 || nbZ == 16) features |= 8;
        }
    }

    // Aromatic
    if (arom) features |= 16;

    // Hydrophobic: non-aromatic C with no heteroatom neighbors, or halogens
    if (z == 6 && !arom) {
        bool hasHetero = false;
        for (int nb : g.neighbors[idx]) {
            int nbZ = g.atomicNum[nb];
            if (nbZ != 6 && nbZ != 1) { hasHetero = true; break; }
        }
        if (!hasHetero) features |= 32;
    }
    if (z == 17 || z == 35 || z == 53) features |= 32; // Cl, Br, I

    return features;
}

/// FCFP-style circular fingerprint using pharmacophoric feature classes.
inline std::vector<uint64_t> computeCircularFingerprintFCFP(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(numWords, 0ULL);
    if (mol.n == 0) return fp;
    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };
    int n = mol.n;
    // Use cached pharmacophore features (v6.5.3 perf fix — parity with Java)
    const auto& pharmFeatures = mol.getPharmacophoreFeatures();
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, static_cast<uint64_t>(pharmFeatures[i]));
        h = fnvMix(h, static_cast<uint64_t>(mol.degree[i]));
        h = fnvMix(h, mol.ring[i] ? 1ULL : 0ULL);
        atomHash[i] = h;
        int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
        fp[bit / 64] |= (1ULL << (bit % 64));
    }
    int maxRadius = (radius < 0) ? n : radius;
    std::vector<uint64_t> prevHash = atomHash;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
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
            int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
            uint64_t mask = 1ULL << (bit % 64);
            if (!(fp[bit / 64] & mask)) anyNew = true;
            fp[bit / 64] |= mask;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return fp;
}

// ============================================================================
// Count-based Circular Fingerprints (ECFP/FCFP counts)
// ============================================================================

/// ECFP count fingerprint: each element of the returned vector holds the number
/// of substructure hashes mapped to that bin position. Preserves multiplicity
/// information lost by binary (OR) folding.
///
/// @param radius  Maximum radius (-1 = unlimited, 2 = ECFP4, 3 = ECFP6)
/// @param fpSize  Number of counter bins (default 2048)
/// @since 6.0.1
inline std::vector<int> computeCircularFingerprintECFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    if (mol.n == 0) return counts;

    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };

    int n = mol.n;
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, mol.atomicNum[i]);                                  // R&H #1
        h = fnvMix(h, mol.degree[i]);                                     // R&H #2
        int bondOrdSum = 0;
        for (int nb : mol.neighbors[i]) bondOrdSum += mol.bondOrder(i, nb);
        h = fnvMix(h, static_cast<uint64_t>(bondOrdSum));                 // R&H #3
        int massNum = (static_cast<int>(mol.massNumber.size()) > i)
                      ? mol.massNumber[i] : 0;
        h = fnvMix(h, static_cast<uint64_t>(massNum));                    // R&H #4
        h = fnvMix(h, static_cast<uint64_t>(mol.formalCharge[i] + 4));    // R&H #5
        int hCount = implicitH(mol.atomicNum[i], bondOrdSum, mol.formalCharge[i]);
        h = fnvMix(h, static_cast<uint64_t>(hCount));                     // R&H #6
        h = fnvMix(h, mol.ring[i] ? 1 : 0);                              // R&H #7
        h = fnvMix(h, mol.aromatic[i] ? 1 : 0);                          // Daylight extension
        if (!mol.tautomerClass.empty() && mol.tautomerClass[i] >= 0)
            h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(mol.tautomerClass[i]));
        atomHash[i] = h;
        int idx = static_cast<int>(h % static_cast<uint64_t>(fpSize));
        counts[idx]++;
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

            int idx = static_cast<int>(h % static_cast<uint64_t>(fpSize));
            counts[idx]++;
            if (seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return counts;
}

/// FCFP count fingerprint using pharmacophoric feature classes. Each element of
/// the returned vector holds the number of substructure hashes mapped to that
/// bin position.
///
/// @param radius  Maximum radius (-1 = unlimited, 2 = FCFP4, 3 = FCFP6)
/// @param fpSize  Number of counter bins (default 2048)
/// @since 6.0.1
inline std::vector<int> computeCircularFingerprintFCFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    if (mol.n == 0) return counts;

    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };

    int n = mol.n;
    const auto& pharmFeatures = mol.getPharmacophoreFeatures();
    std::vector<uint64_t> atomHash(n);
    for (int i = 0; i < n; ++i) {
        uint64_t h = FNV1A_SEED;
        h = fnvMix(h, static_cast<uint64_t>(pharmFeatures[i]));
        h = fnvMix(h, static_cast<uint64_t>(mol.degree[i]));
        h = fnvMix(h, mol.ring[i] ? 1ULL : 0ULL);
        atomHash[i] = h;
        int idx = static_cast<int>(h % static_cast<uint64_t>(fpSize));
        counts[idx]++;
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

            int idx = static_cast<int>(h % static_cast<uint64_t>(fpSize));
            counts[idx]++;
            if (seenHashes.insert(h).second) anyNew = true;
        }
        std::swap(prevHash, nextHash);
        if (radius < 0 && !anyNew) break;
    }
    return counts;
}

// ============================================================================
// MCS-aware path fingerprint
// ============================================================================

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

inline void setBit(std::vector<uint64_t>& fp, int hash, int fpSize) {
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
        setBit(fp, hashPath(g, tautAware, path, depth + 1), fpSize);
        if (depth < maxDepth)
            enumeratePaths(g, tautAware, isHeavy, fp, fpSize, visited,
                           path, depth + 1, maxDepth);
        visited[nb] = false;
    }
}

} // namespace mcs_fp_detail

/// MCS-aware path fingerprint: encodes element, ring, aromatic, tautomer class,
/// degree per atom; bond order, ring, aromatic per bond. Paths are hashed
/// canonically (direction-independent). Compatible with the Java SMSD MCS fingerprint.
///
/// @param pathLength  Maximum path length in bonds (default 7)
/// @param fpSize      Fingerprint size in bits (default 2048)
inline std::vector<uint64_t> computeMcsFingerprint(
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
        mcs_fp_detail::setBit(fp, mcs_fp_detail::hashPath(mol, tautAware, pathBuf.data(), 1), fpSize);
        mcs_fp_detail::enumeratePaths(mol, tautAware, isHeavy, fp, fpSize, visited,
                                      pathBuf.data(), 1, pathLength);
        visited[start] = false;
    }
    return fp;
}

// ============================================================================
// Topological Torsion Fingerprint (Nilakantan et al. 1987)
// ============================================================================

/// Compute topological torsion atom type: encodes (atomicNum, numHeavyNeighbors,
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

/// Topological torsion fingerprint (binary). Enumerates all 4-atom linear paths
/// A-B-C-D, hashes the atom types using FNV-1a with canonical ordering
/// (min of forward and reverse hash), and folds into fpSize bits.
///
/// Reference: Nilakantan et al., J. Chem. Inf. Comput. Sci. 1987, 27, 82-85.
///
/// @param fpSize  Fingerprint size in bits
/// @since 6.0.1
inline std::vector<uint64_t> computeTopologicalTorsion(
    const MolGraph& mol, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    int numWords = (fpSize + 63) / 64;
    std::vector<uint64_t> fp(numWords, 0ULL);
    int n = mol.n;
    if (n < 4) return fp;

    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };

    // Precompute atom types
    std::vector<int> atomType(n);
    for (int i = 0; i < n; ++i)
        atomType[i] = ttAtomType(mol, i);

    // Enumerate all 4-atom linear paths: A -> B -> C -> D
    for (int a = 0; a < n; ++a) {
        for (int b : mol.neighbors[a]) {
            for (int c : mol.neighbors[b]) {
                if (c == a) continue;
                for (int d : mol.neighbors[c]) {
                    if (d == b || d == a) continue;
                    // Hash forward: A-B-C-D
                    uint64_t fwd = FNV1A_SEED;
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[a]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[b]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[c]));
                    fwd = fnvMix(fwd, static_cast<uint64_t>(atomType[d]));
                    // Hash reverse: D-C-B-A
                    uint64_t rev = FNV1A_SEED;
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[d]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[c]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[b]));
                    rev = fnvMix(rev, static_cast<uint64_t>(atomType[a]));
                    // Canonical: take minimum
                    uint64_t h = (fwd <= rev) ? fwd : rev;
                    int bit = static_cast<int>(h % static_cast<uint64_t>(fpSize));
                    fp[bit / 64] |= (1ULL << (bit % 64));
                }
            }
        }
    }
    return fp;
}

/// Topological torsion count fingerprint. Each element of the returned vector
/// holds the number of distinct 4-atom torsion paths that hash to that bin.
///
/// @param fpSize  Number of counter bins
/// @since 6.0.1
inline std::vector<int> computeTopologicalTorsionCounts(
    const MolGraph& mol, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    int n = mol.n;
    if (n < 4) return counts;

    auto fnvMix = [](uint64_t h, uint64_t val) -> uint64_t {
        h ^= val; h *= FNV1A_PRIME; return h;
    };

    std::vector<int> atomType(n);
    for (int i = 0; i < n; ++i)
        atomType[i] = ttAtomType(mol, i);

    // Enumerate canonical paths only (a < d) to avoid double-counting
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
                    int idx = static_cast<int>(h % static_cast<uint64_t>(fpSize));
                    counts[idx]++;
                }
            }
        }
    }
    return counts;
}

} // namespace detail

// ============================================================================
// Batch Substructure: 1 query vs N targets
// ============================================================================

/// Check if query is a substructure of each target molecule.
/// Returns a boolean vector: true if query is substructure of targets[i].
///
/// Thread-safe: each thread runs its own VF2PP instance with independent state.
/// @param numThreads  0 = auto (all cores), >0 = specific thread count.
inline std::vector<bool> batchSubstructure(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& opts,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());

    if (N == 0) return {};
    if (query.n == 0) return std::vector<bool>(N, true);

    // Eagerly initialize ALL lazy fields before entering the parallel region.
    // MolGraph lazy caches are not thread-safe on first construction.
    detail::prewarmGraph(query);
    for (const MolGraph& t : targets) {
        if (t.n >= query.n) detail::prewarmGraph(t);
    }

    // Use uint8_t instead of vector<bool> — vector<bool> packs bits into
    // shared bytes, causing data races under OpenMP concurrent writes.
    std::vector<uint8_t> buf(N, 0);

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads; // suppress unused warning when no OpenMP

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n >= query.n) {
            buf[i] = isSubstructure(query, targets[i], opts) ? 1 : 0;
        }
    }

    return std::vector<bool>(buf.begin(), buf.end());
}

// ============================================================================
// Batch Find Substructure: 1 query vs N targets (with atom mappings)
// ============================================================================

/// Find substructure atom-atom mappings for query against each target molecule.
/// Returns a vector of atom mappings (query atom, target atom) for each target.
/// Empty inner vector means no substructure match.
///
/// Thread-safe: each thread runs its own VF2PP instance with independent state.
/// @param numThreads  0 = auto (all cores), >0 = specific thread count.
inline std::vector<std::vector<std::pair<int,int>>> batchFindSubstructure(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& opts,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());
    std::vector<std::vector<std::pair<int,int>>> results(N);

    if (N == 0) return results;
    if (query.n == 0) return results;  // empty query: no meaningful mapping

    // Eagerly initialize ALL lazy fields before entering the parallel region.
    detail::prewarmGraph(query);
    for (const MolGraph& t : targets) {
        if (t.n >= query.n) detail::prewarmGraph(t);
    }

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads; // suppress unused warning when no OpenMP

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n >= query.n) {
            results[i] = findSubstructure(query, targets[i], opts);
        }
    }

    return results;
}

// ============================================================================
// Batch MCS: 1 query vs N targets
// ============================================================================

/// Find the MCS between query and each target molecule.
/// Returns a vector of mappings (query atom -> target atom).
/// Empty map means no significant common substructure found.
///
/// Thread-safe: each thread creates its own MCS scratch buffers internally.
/// @param numThreads  0 = auto, >0 = specific thread count.
inline std::vector<std::map<int,int>> batchMCS(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());
    std::vector<std::map<int,int>> results(N);

    if (N == 0 || query.n == 0) return results;

    // Eagerly initialize all lazy fields before parallel region.
    detail::prewarmGraph(query);
    for (const MolGraph& t : targets) detail::prewarmGraph(t);

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n > 0) {
            results[i] = findMCS(query, targets[i], chem, opts);
        }
    }

    return results;
}

/// Find the MCS size between query and each target molecule.
/// Returns atom counts only, avoiding mapping materialization at the API
/// boundary when callers only need screening sizes.
/// Phase 2.4: Uses findMCSSize() to avoid constructing std::map<int,int>
/// at the caller boundary.
inline std::vector<int> batchMCSSize(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());
    std::vector<int> results(N, 0);

    if (N == 0 || query.n == 0) return results;

    detail::prewarmGraph(query);
    for (const MolGraph& t : targets) detail::prewarmGraph(t);

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n > 0) {
            results[i] = findMCSSize(query, targets[i], chem, opts);
        }
    }

    return results;
}

// ============================================================================
// Screen and Match: RASCAL pre-screen + exact MCS on hits
// ============================================================================

/// Two-phase pipeline:
///   Phase 1: RASCAL similarity upper-bound screening (parallel)
///   Phase 2: Exact MCS computation on hits only (parallel)
///
/// Returns vector of (target_index, MCS_mapping) for targets that pass
/// the RASCAL threshold AND have a non-empty MCS.
///
/// This is the recommended workflow for large databases: RASCAL screening
/// eliminates ~90% of targets cheaply, then exact MCS runs only on candidates.
inline std::vector<std::pair<int, std::map<int,int>>> screenAndMatch(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts,
    double rascalThreshold = 0.3,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());
    if (N == 0 || query.n == 0) return {};

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

    // --- Phase 1: RASCAL screening (parallel) ---
    std::vector<double> scores(N, 0.0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n > 0) {
            scores[i] = similarityUpperBound(query, targets[i]);
        }
    }

    // Collect hits
    std::vector<int> hits;
    hits.reserve(N / 10);
    for (int i = 0; i < N; ++i) {
        if (scores[i] >= rascalThreshold) {
            hits.push_back(i);
        }
    }

    if (hits.empty()) return {};

    detail::prewarmGraph(query);
    for (int idx : hits) detail::prewarmGraph(targets[idx]);

    // --- Phase 2: Exact MCS on hits only (parallel) ---
    const int H = static_cast<int>(hits.size());
    std::vector<std::map<int,int>> mcsMaps(H);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int h = 0; h < H; ++h) {
        int idx = hits[h];
        mcsMaps[h] = findMCS(query, targets[idx], chem, opts);
    }

    // Collect non-empty results
    std::vector<std::pair<int, std::map<int,int>>> results;
    results.reserve(H);
    for (int h = 0; h < H; ++h) {
        if (!mcsMaps[h].empty()) {
            results.emplace_back(hits[h], std::move(mcsMaps[h]));
        }
    }

    return results;
}

/// Two-phase RASCAL screening followed by exact MCS size computation on hits.
/// Returns (target_index, mcs_size) for each screened-in target with a non-zero
/// exact MCS.
inline std::vector<std::pair<int, int>> screenAndMCSSize(
    const MolGraph& query,
    const std::vector<MolGraph>& targets,
    const ChemOptions& chem,
    const MCSOptions& opts,
    double rascalThreshold = 0.3,
    int numThreads = 0)
{
    const int N = static_cast<int>(targets.size());
    if (N == 0 || query.n == 0) return {};

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

    std::vector<double> scores(N, 0.0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        if (targets[i].n > 0) {
            scores[i] = similarityUpperBound(query, targets[i]);
        }
    }

    std::vector<int> hits;
    hits.reserve(N / 10);
    for (int i = 0; i < N; ++i) {
        if (scores[i] >= rascalThreshold) {
            hits.push_back(i);
        }
    }

    if (hits.empty()) return {};

    detail::prewarmGraph(query);
    for (int idx : hits) detail::prewarmGraph(targets[idx]);

    const int H = static_cast<int>(hits.size());
    std::vector<int> mcsSizes(H, 0);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
    for (int h = 0; h < H; ++h) {
        int idx = hits[h];
        mcsSizes[h] = static_cast<int>(findMCS(query, targets[idx], chem, opts).size());
    }

    std::vector<std::pair<int, int>> results;
    results.reserve(H);
    for (int h = 0; h < H; ++h) {
        if (mcsSizes[h] > 0) {
            results.emplace_back(hits[h], mcsSizes[h]);
        }
    }

    return results;
}

// ============================================================================
// Batch Fingerprint Generation
// ============================================================================

/// Generate path-based fingerprints for all molecules in parallel.
/// Each fingerprint is a vector of uint64_t words (fpSize / 64 words).
///
/// @param pathLength  Maximum path length to enumerate (default 7).
/// @param fpSize      Fingerprint size in bits (default 1024, must be multiple of 64).
/// @param numThreads  0 = auto, >0 = specific thread count.
inline std::vector<std::vector<uint64_t>> batchFingerprint(
    const std::vector<MolGraph>& mols,
    int pathLength = 7,
    int fpSize = 1024,
    int numThreads = 0)
{
    const int N = static_cast<int>(mols.size());
    std::vector<std::vector<uint64_t>> results(N);

    if (N == 0) return results;

    // Ensure fpSize is a multiple of 64
    fpSize = ((fpSize + 63) / 64) * 64;

    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 8) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        results[i] = detail::computePathFingerprint(mols[i], pathLength, fpSize);
    }

    return results;
}

// ============================================================================
// Utility: Tanimoto similarity between two fingerprints
// ============================================================================

/// Compute Tanimoto coefficient between two fingerprints.
inline double fingerprintTanimoto(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;

    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));

    int andBits = 0, orBits = 0;

    // Use the portable popcount64 from bitops.hpp (no need for inline lambda)
    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        orBits  += popcount64(a[i] | b[i]);
    }

    // Extra words in the longer fingerprint contribute only to OR
    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i) {
        orBits += popcount64(longer[i]);
    }

    return (orBits == 0) ? 1.0 : static_cast<double>(andBits) /
                                  static_cast<double>(orBits);
}

/// Compute Dice coefficient between two fingerprints: 2*|A&B| / (|A| + |B|).
/// Dice weights shared features more heavily than Tanimoto, suitable for
/// scaffold-hopping searches and count-based fingerprint comparisons.
inline double fingerprintDice(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;

    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));

    int andBits = 0, aBits = 0, bBits = 0;

    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        aBits   += popcount64(a[i]);
        bBits   += popcount64(b[i]);
    }

    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i) {
        if (&longer == &a) aBits += popcount64(longer[i]);
        else               bBits += popcount64(longer[i]);
    }

    int sum = aBits + bBits;
    return (sum == 0) ? 0.0 : (2.0 * andBits) / sum;
}

/// Compute Cosine similarity between two fingerprints: |A&B| / sqrt(|A| * |B|).
/// Normalises for fingerprint density differences between molecules of varying size.
inline double fingerprintCosine(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 0.0;

    int minWords = static_cast<int>(std::min(a.size(), b.size()));
    int maxWords = static_cast<int>(std::max(a.size(), b.size()));

    int andBits = 0, aBits = 0, bBits = 0;

    for (int i = 0; i < minWords; ++i) {
        andBits += popcount64(a[i] & b[i]);
        aBits   += popcount64(a[i]);
        bBits   += popcount64(b[i]);
    }

    const auto& longer = (a.size() > b.size()) ? a : b;
    for (int i = minWords; i < maxWords; ++i) {
        if (&longer == &a) aBits += popcount64(longer[i]);
        else               bBits += popcount64(longer[i]);
    }

    double denom = std::sqrt(static_cast<double>(aBits) * bBits);
    return (denom == 0.0) ? 0.0 : static_cast<double>(andBits) / denom;
}

/// Compute Soergel distance between two fingerprints: 1 - Tanimoto.
/// A proper metric distance useful for clustering and k-NN algorithms.
inline double fingerprintSoergel(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    if (a.empty() || b.empty()) return 1.0; // max distance for empty inputs
    return 1.0 - fingerprintTanimoto(a, b);
}

// ---- Count-vector similarity metrics ----

/// Count-vector Tanimoto: sum(min(a,b)) / sum(max(a,b)).
/// Generalisation of Tanimoto for integer count vectors (e.g. ECFP counts).
inline double countTanimoto(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long minSum = 0, maxSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        minSum += std::min(a[i], b[i]);
        maxSum += std::max(a[i], b[i]);
    }
    return (maxSum == 0) ? 0.0 : static_cast<double>(minSum) / maxSum;
}

/// Count-vector Dice: 2*sum(min(a,b)) / (sum(a) + sum(b)).
inline double countDice(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long minSum = 0, aSum = 0, bSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        minSum += std::min(a[i], b[i]);
        aSum   += a[i];
        bSum   += b[i];
    }
    long long denom = aSum + bSum;
    return (denom == 0) ? 0.0 : (2.0 * minSum) / denom;
}

/// Count-vector Cosine: dot(a,b) / (|a| * |b|).
inline double countCosine(
    const std::vector<int>& a,
    const std::vector<int>& b)
{
    if (a.empty() || b.empty() || a.size() != b.size()) return 0.0;
    long long dot = 0, aSqSum = 0, bSqSum = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        dot    += static_cast<long long>(a[i]) * b[i];
        aSqSum += static_cast<long long>(a[i]) * a[i];
        bSqSum += static_cast<long long>(b[i]) * b[i];
    }
    double denom = std::sqrt(static_cast<double>(aSqSum) * bSqSum);
    return (denom == 0.0) ? 0.0 : static_cast<double>(dot) / denom;
}

/// Check if fingerprint a is a subset of fingerprint b (all ON bits of a are ON in b).
/// Used for substructure pre-screening: if query is substructure of target,
/// then query FP must be a subset of target FP.
inline bool fingerprintSubset(
    const std::vector<uint64_t>& query,
    const std::vector<uint64_t>& target)
{
    if (query.empty()) return true;

    int qWords = static_cast<int>(query.size());
    int tWords = static_cast<int>(target.size());

    for (int i = 0; i < qWords; ++i) {
        uint64_t qw = query[i];
        uint64_t tw = (i < tWords) ? target[i] : 0ULL;
        if ((qw & tw) != qw) return false;
    }
    return true;
}

// ============================================================================
// Batch Fingerprint Subset Screen (CPU)
// ============================================================================

/// CPU-based fingerprint subset screening: 1 query vs N targets.
/// Returns indices of targets where query FP is a subset of target FP.
inline std::vector<int> batchFingerprintScreen(
    const std::vector<uint64_t>& queryFP,
    const std::vector<std::vector<uint64_t>>& targetFPs,
    int numThreads = 0)
{
    const int N = static_cast<int>(targetFPs.size());
    if (N == 0 || queryFP.empty()) return {};

    std::vector<int> mask(N, 0);
    int nThreads = detail::resolveThreads(numThreads);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 64) num_threads(nThreads)
#endif
    for (int i = 0; i < N; ++i) {
        mask[i] = fingerprintSubset(queryFP, targetFPs[i]) ? 1 : 0;
    }

    std::vector<int> hits;
    for (int i = 0; i < N; ++i) {
        if (mask[i]) hits.push_back(i);
    }
    return hits;
}

// ============================================================================
// Fingerprint Format Conversions
// ============================================================================

/// Convert a fingerprint to a hexadecimal string (lowercase, zero-padded).
inline std::string toHex(const std::vector<uint64_t>& fp) {
    std::string out;
    out.reserve(fp.size() * 16);
    for (uint64_t w : fp) {
        char buf[17];
        std::snprintf(buf, sizeof(buf), "%016" PRIx64, w);
        out.append(buf, 16);
    }
    return out;
}

/// Parse a hexadecimal string back into a fingerprint.
inline std::vector<uint64_t> fromHex(const std::string& hex) {
    int words = (static_cast<int>(hex.size()) + 15) / 16;
    std::vector<uint64_t> fp(words, 0ULL);
    for (int i = 0; i < words; ++i) {
        int start = i * 16;
        int len = std::min(16, static_cast<int>(hex.size()) - start);
        std::string chunk = hex.substr(start, len);
        fp[i] = std::strtoull(chunk.c_str(), nullptr, 16);
    }
    return fp;
}

/// Convert a fingerprint to a binary string of '0' and '1' characters.
inline std::string toBinaryString(const std::vector<uint64_t>& fp, int fpSize) {
    std::string out(fpSize, '0');
    for (int i = 0; i < fpSize; ++i) {
        if (fp[i / 64] & (1ULL << (i % 64))) out[i] = '1';
    }
    return out;
}

// ===========================================================================
// fingerprintFromSmiles -- convenience: SMILES -> ECFP fingerprint
// ===========================================================================

/**
 * All-in-one convenience: parse a SMILES string and compute its ECFP
 * fingerprint as a bit-packed uint64_t vector.
 *
 * @param smi     SMILES string
 * @param radius  ECFP radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
 * @param fpSize  fingerprint size in bits (default 2048)
 * @return bit-packed fingerprint vector
 * @since 6.3.0
 */
inline std::vector<uint64_t> fingerprintFromSmiles(const std::string& smi,
                                                    int radius = 2,
                                                    int fpSize = 2048) {
    MolGraph g = smsd::parseSMILES(smi);
    return detail::computeCircularFingerprintECFP(g, radius, fpSize);
}

// ===========================================================================
// countsToArray -- sparse count map -> dense array
// ===========================================================================

/**
 * Convert a sparse count map {bitPosition -> count} into a dense integer
 * array of length fpSize. Useful for serialising count-based fingerprints
 * to fixed-length vectors suitable for database storage.
 *
 * Positions outside [0, fpSize) are silently ignored.
 *
 * @param counts  sparse count map (bitPosition to count)
 * @param fpSize  target array length
 * @return dense count vector of length fpSize
 * @since 6.3.0
 */
inline std::vector<int> countsToArray(const std::map<int,int>& counts, int fpSize) {
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> arr(fpSize, 0);
    for (auto& kv : counts) {
        if (kv.first >= 0 && kv.first < fpSize) {
            arr[kv.first] = kv.second;
        }
    }
    return arr;
}

// ============================================================================
// TargetCorpus: pre-warmed, fingerprinted target collection
// ============================================================================

/// A reusable collection of target molecules with pre-computed invariants
/// and fingerprints. Amortises the cost of prewarm + fingerprint computation
/// across many queries.
///
/// Usage:
///   TargetCorpus corpus;
///   corpus.addTargets(targets);
///   corpus.prewarm();  // one-time cost
///   auto hits = corpus.substructure(query, opts);       // reuse many times
///   auto maps = corpus.findSubstructure(query, opts);
///   auto sizes = corpus.mcsSize(query, chem, mcsOpts);
///   auto idxs  = corpus.screen(query, 0.7);
class TargetCorpus {
    std::vector<MolGraph> targets_;
    std::vector<std::vector<uint64_t>> fingerprints_;
    bool prewarmed_ = false;

public:
    /// Add a single target molecule.
    void addTarget(MolGraph g) {
        targets_.push_back(std::move(g));
        prewarmed_ = false;
    }

    /// Add multiple target molecules.
    void addTargets(std::vector<MolGraph> gs) {
        targets_.reserve(targets_.size() + gs.size());
        for (auto& g : gs) {
            targets_.push_back(std::move(g));
        }
        prewarmed_ = false;
    }

    /// Pre-compute lazy graph invariants and ECFP fingerprints for all targets.
    /// Call once after loading all targets. Thread-safe to call from one thread;
    /// subsequent query methods are safe to call concurrently.
    void prewarm(int numThreads = 0) {
        const int N = static_cast<int>(targets_.size());
        if (N == 0) { prewarmed_ = true; return; }

        int nThreads = detail::resolveThreads(numThreads);
        (void)nThreads;

        // Prewarm graph invariants in parallel
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads)
#endif
        for (int i = 0; i < N; ++i) {
            detail::prewarmGraph(targets_[i]);
        }

        // Compute ECFP fingerprints in parallel
        fingerprints_.resize(N);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 8) num_threads(nThreads)
#endif
        for (int i = 0; i < N; ++i) {
            fingerprints_[i] = detail::computeCircularFingerprintECFP(
                targets_[i], /*radius=*/2, /*fpSize=*/2048);
        }

        prewarmed_ = true;
    }

    /// Number of targets in the corpus.
    size_t size() const { return targets_.size(); }

    /// Whether prewarm() has been called since the last addTarget/addTargets.
    bool isPrewarmed() const { return prewarmed_; }

    /// Check if query is a substructure of each target.
    /// Returns boolean hit mask. Automatically prewarms query.
    std::vector<bool> substructure(
        const MolGraph& query, const ChemOptions& opts, int nThreads = 0) const
    {
        const int N = static_cast<int>(targets_.size());
        if (N == 0) return {};
        if (query.n == 0) return std::vector<bool>(N, true);

        detail::prewarmGraph(query);
        // If not prewarmed, prewarm targets inline (const-safe: lazy caches are mutable)
        if (!prewarmed_) {
            for (const MolGraph& t : targets_) {
                if (t.n >= query.n) detail::prewarmGraph(t);
            }
        }

        std::vector<uint8_t> buf(N, 0);
        int nt = detail::resolveThreads(nThreads);
        (void)nt;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nt)
#endif
        for (int i = 0; i < N; ++i) {
            if (targets_[i].n >= query.n) {
                buf[i] = isSubstructure(query, targets_[i], opts) ? 1 : 0;
            }
        }
        return std::vector<bool>(buf.begin(), buf.end());
    }

    /// Find substructure atom-atom mappings for query against each target.
    /// Empty inner vector = no match.
    std::vector<std::vector<std::pair<int,int>>> findSubstructure(
        const MolGraph& query, const ChemOptions& opts, int nThreads = 0) const
    {
        const int N = static_cast<int>(targets_.size());
        std::vector<std::vector<std::pair<int,int>>> results(N);
        if (N == 0 || query.n == 0) return results;

        detail::prewarmGraph(query);
        if (!prewarmed_) {
            for (const MolGraph& t : targets_) {
                if (t.n >= query.n) detail::prewarmGraph(t);
            }
        }

        int nt = detail::resolveThreads(nThreads);
        (void)nt;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nt)
#endif
        for (int i = 0; i < N; ++i) {
            if (targets_[i].n >= query.n) {
                results[i] = smsd::findSubstructure(query, targets_[i], opts);
            }
        }
        return results;
    }

    /// Find MCS sizes between query and each target.
    std::vector<int> mcsSize(
        const MolGraph& query, const ChemOptions& chem,
        const MCSOptions& opts, int nThreads = 0) const
    {
        const int N = static_cast<int>(targets_.size());
        std::vector<int> results(N, 0);
        if (N == 0 || query.n == 0) return results;

        detail::prewarmGraph(query);
        if (!prewarmed_) {
            for (const MolGraph& t : targets_) detail::prewarmGraph(t);
        }

        int nt = detail::resolveThreads(nThreads);
        (void)nt;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nt)
#endif
        for (int i = 0; i < N; ++i) {
            if (targets_[i].n > 0) {
                results[i] = static_cast<int>(
                    findMCS(query, targets_[i], chem, opts).size());
            }
        }
        return results;
    }

    /// Fingerprint-based Tanimoto screening: returns indices of targets
    /// whose Tanimoto similarity to the query exceeds the threshold.
    /// Requires prewarm() to have been called (for fingerprints).
    std::vector<int> screen(
        const MolGraph& query, double threshold, int nThreads = 0) const
    {
        const int N = static_cast<int>(targets_.size());
        if (N == 0) return {};

        // Compute query fingerprint
        auto queryFP = detail::computeCircularFingerprintECFP(
            query, /*radius=*/2, /*fpSize=*/2048);

        // If fingerprints are not computed yet, fall back to empty
        if (fingerprints_.empty() || static_cast<int>(fingerprints_.size()) != N) {
            return {};
        }

        std::vector<uint8_t> mask(N, 0);
        int nt = detail::resolveThreads(nThreads);
        (void)nt;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 64) num_threads(nt)
#endif
        for (int i = 0; i < N; ++i) {
            if (fingerprintTanimoto(queryFP, fingerprints_[i]) >= threshold) {
                mask[i] = 1;
            }
        }

        std::vector<int> hits;
        for (int i = 0; i < N; ++i) {
            if (mask[i]) hits.push_back(i);
        }
        return hits;
    }
};

} // namespace batch
} // namespace smsd

#endif // SMSD_BATCH_HPP

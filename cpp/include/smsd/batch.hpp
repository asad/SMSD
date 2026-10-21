/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only OpenMP batch processing for SMSD.
 *
 * Provides parallel fingerprint operations and reusable target storage.
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

#include "smsd/bitops.hpp"
#include "smsd/mol_graph.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <stdexcept>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <limits>
#include <mutex>
#include <shared_mutex>
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

/// Non-owning targets for bindings; the caller retains every graph's lifetime.
class MolGraphView {
    std::vector<std::reference_wrapper<const MolGraph>> graphs_;
public:
    void reserve(size_t size) { graphs_.reserve(size); }
    void push_back(const MolGraph& graph) { graphs_.emplace_back(graph); }
    size_t size() const { return graphs_.size(); }
    const MolGraph& operator[](size_t index) const { return graphs_[index].get(); }
    auto begin() const { return graphs_.begin(); }
    auto end() const { return graphs_.end(); }
};

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

/// Precompute the lazy MolGraph invariants used by matching and fingerprints.
inline void prewarmGraph(const MolGraph& g) {
    g.prewarm();
}

// ============================================================================
// Path fingerprints: one adjacency builder and one path walk
// ============================================================================

namespace path_fp_detail {

inline uint64_t mixTerm(uint64_t hash, uint32_t term) {
    hash ^= static_cast<uint64_t>(term);
    return hash * FNV1A_PRIME;
}

/// Atoms kept for a path fingerprint, their neighbours in graph order and one term per bond,
/// stored contiguously.  Atoms that are not kept have no neighbours.
struct PathAdjacency {
    std::vector<int> offset, atom;
    std::vector<uint32_t> bondTerm;
    std::vector<uint8_t> kept;
};

/// Builds the adjacency of the atoms for which `keep(atom)` holds, with `bondTerm(from, to)` per bond.
template<class Keep, class BondTerm>
inline PathAdjacency pathAdjacency(const MolGraph& mol, Keep keep, BondTerm bondTerm) {
    PathAdjacency adjacency;
    adjacency.kept.assign(mol.n, 0);
    for (int i = 0; i < mol.n; ++i) adjacency.kept[i] = keep(i) ? 1 : 0;
    adjacency.offset.assign(mol.n + 1, 0);
    for (int i = 0; i < mol.n; ++i) {
        adjacency.offset[i] = static_cast<int>(adjacency.atom.size());
        if (!adjacency.kept[i]) continue;
        for (int nb : mol.neighbors[i]) {
            if (!adjacency.kept[nb]) continue;
            adjacency.atom.push_back(nb);
            adjacency.bondTerm.push_back(bondTerm(i, nb));
        }
    }
    adjacency.offset[mol.n] = static_cast<int>(adjacency.atom.size());
    return adjacency;
}

/// Rolling form of the direction-independent path hash.  The forward hash folds the
/// interleaved atom and bond terms from the first atom (h = 31h + term); the reverse
/// hash folds the same terms from the last atom, so appending a term adds term * 31^k.
/// Both are exact in 32-bit modular arithmetic, so a path costs O(1) instead of O(length).
struct RollingPathHash {
    uint32_t forward = 17, reverseSum = 0, power = 1;   // reverse = 17 * power + reverseSum

    RollingPathHash append(uint32_t term) const {
        return {forward * 31u + term, reverseSum + term * power, power * 31u};
    }
    int value() const {
        const uint32_t reverse = 17u * power + reverseSum;
        return static_cast<int>(std::min(forward & 0x7FFFFFFFu, reverse & 0x7FFFFFFFu));
    }
};

/// Rolling hash of the atom terms alone; bond terms are not hashed.
struct AtomPathHash {
    RollingPathHash hash;
    static AtomPathHash start(uint32_t atom, uint32_t*) { return {RollingPathHash{}.append(atom)}; }
    AtomPathHash extend(uint32_t, uint32_t atom) const { return {hash.append(atom)}; }
    int bit(int fpSize) const { return hash.value() % fpSize; }
};

/// Rolling hash of the interleaved atom and bond terms.
struct LabelledPathHash {
    RollingPathHash hash;
    static LabelledPathHash start(uint32_t atom, uint32_t*) { return {RollingPathHash{}.append(atom)}; }
    LabelledPathHash extend(uint32_t bond, uint32_t atom) const { return {hash.append(bond).append(atom)}; }
    int bit(int fpSize) const { return hash.value() % fpSize; }
};

/// FNV-1a hash of the interleaved atom and bond terms [atom0, bond01, atom1, ...] in the
/// canonical direction.  The forward hash is carried down the walk; the reverse hash is
/// folded from the end of the shared term buffer when a bit is set.
struct FoldedPathHash {
    uint32_t* terms;
    int length;
    uint64_t forward;

    static FoldedPathHash start(uint32_t atom, uint32_t* terms) {
        terms[0] = atom;
        return {terms, 1, mixTerm(FNV1A_SEED, atom)};
    }
    FoldedPathHash extend(uint32_t bond, uint32_t atom) const {
        terms[length] = bond;
        terms[length + 1] = atom;
        return {terms, length + 2, mixTerm(mixTerm(forward, bond), atom)};
    }
    int bit(int fpSize) const {
        uint64_t reverse = FNV1A_SEED;
        for (int k = length - 1; k >= 0; --k) reverse = mixTerm(reverse, terms[k]);
        return static_cast<int>(std::min(forward, reverse) % static_cast<uint64_t>(fpSize));
    }
};

inline void setBit(std::vector<uint64_t>& fp, int bit) {
    fp[bit >> 6] |= uint64_t(1) << (bit & 63);
}

/// Sets a bit for every simple path extending the current one.  Each path of two or more
/// atoms is reached from both ends with the same canonical hash, so only the copy ending
/// above its start atom sets a bit.
template<class Hash>
inline void extendPaths(const PathAdjacency& adjacency, const std::vector<uint32_t>& atomTerm,
                        std::vector<uint64_t>& fp, int fpSize, std::vector<uint8_t>& visited,
                        int start, int current, const Hash& hash, int bonds, int maxBonds) {
    for (int k = adjacency.offset[current]; k < adjacency.offset[current + 1]; ++k) {
        const int next = adjacency.atom[k];
        if (visited[next]) continue;
        const Hash extended = hash.extend(adjacency.bondTerm[k], atomTerm[next]);
        if (next > start) setBit(fp, extended.bit(fpSize));
        if (bonds + 1 < maxBonds) {
            visited[next] = 1;
            extendPaths(adjacency, atomTerm, fp, fpSize, visited, start, next, extended, bonds + 1, maxBonds);
            visited[next] = 0;
        }
    }
}

/// Folds every simple path of up to `maxBonds` bonds between kept atoms into `fpSize` bits.
template<class Hash>
inline std::vector<uint64_t> pathFingerprint(const PathAdjacency& adjacency,
                                             const std::vector<uint32_t>& atomTerm,
                                             int maxBonds, int fpSize) {
    std::vector<uint64_t> fp((fpSize - 1) / 64 + 1, 0ULL);
    const int n = static_cast<int>(adjacency.kept.size());
    std::vector<uint8_t> visited(n, 0);
    std::vector<uint32_t> terms(2 * static_cast<size_t>(std::min(maxBonds, n)) + 1);
    for (int start = 0; start < n; ++start) {
        if (!adjacency.kept[start]) continue;
        const Hash hash = Hash::start(atomTerm[start], terms.data());
        setBit(fp, hash.bit(fpSize));
        visited[start] = 1;
        extendPaths(adjacency, atomTerm, fp, fpSize, visited, start, start, hash, 0, maxBonds);
        visited[start] = 0;
    }
    return fp;
}

} // namespace path_fp_detail

/// Path fingerprint over all atoms, including explicit hydrogens.  Enumerates all paths of
/// up to pathLength bonds (at most 7), hashes each path's interleaved atom labels and bond
/// orders in canonical direction, and folds the hashes into fpSize bits.  Used by
/// batchFingerprint.  Returns a vector of uint64_t words.
inline std::vector<uint64_t> computePathFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    if (pathLength <= 0) throw std::invalid_argument("pathLength must be positive");
    pathLength = std::min(pathLength, 7); // paths have at most seven bonds
    const auto adjacency = path_fp_detail::pathAdjacency(mol, [](int) { return true; },
        [&](int from, int to) { return static_cast<uint32_t>(mol.bondOrder(from, to) + 1000); });
    std::vector<uint32_t> atomTerm(mol.n);
    for (int i = 0; i < mol.n; ++i) atomTerm[i] = static_cast<uint32_t>(mol.label[i]);
    return path_fp_detail::pathFingerprint<path_fp_detail::FoldedPathHash>(adjacency, atomTerm,
                                                                           pathLength, fpSize);
}

/// Element path fingerprint over heavy atoms (the Python path_fingerprint): each path of up to
/// pathLength bonds is hashed by its atomic-number sequence in canonical direction.  Bonds
/// are not hashed.  Matches the Java FingerprintEngine.pathFingerprint bits.
inline std::vector<uint64_t> computeElementPathFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    if (pathLength <= 0 || pathLength == std::numeric_limits<int>::max())
        throw std::invalid_argument("pathLength must be positive and representable");
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    const auto adjacency = path_fp_detail::pathAdjacency(mol,
        [&](int atom) { return mol.atomicNum[atom] != 1; }, [](int, int) { return 0u; });
    std::vector<uint32_t> atomTerm(mol.n);
    for (int i = 0; i < mol.n; ++i) atomTerm[i] = static_cast<uint32_t>(mol.atomicNum[i]);
    return path_fp_detail::pathFingerprint<path_fp_detail::AtomPathHash>(adjacency, atomTerm,
                                                                         pathLength, fpSize);
}

/** Largest shortest-path distance within any connected component. */
inline int maximumComponentExtent(const MolGraph& graph) {
    int extent = 0;
    std::vector<int> distance(graph.n), queue(graph.n);
    for (int start = 0; start < graph.n; ++start) {
        std::fill(distance.begin(), distance.end(), -1);
        size_t head = 0, tail = 0;
        queue[tail++] = start;
        distance[start] = 0;
        while (head < tail) {
            const int atom = queue[head++];
            for (int neighbor : graph.neighbors[atom]) {
                if (distance[neighbor] >= 0) continue;
                distance[neighbor] = distance[atom] + 1;
                extent = std::max(extent, distance[neighbor]);
                queue[tail++] = neighbor;
            }
        }
    }
    return extent;
}

// ============================================================================
// Circular (Morgan) fingerprints: one iteration, two atom invariants, two sinks
// ============================================================================

inline uint64_t fnvMix(uint64_t hash, uint64_t value) {
    hash ^= value;
    return hash * FNV1A_PRIME;
}

/// Initial ECFP atom invariant (Rogers and Hahn 2010, Table 1, with the Daylight aromaticity
/// extension and the tautomer class when one is assigned).
inline uint64_t rogersHahnInvariant(const MolGraph& mol, int i) {
    uint64_t h = FNV1A_SEED;
    h = fnvMix(h, mol.atomicNum[i]);                                  // R&H #1: atomic number
    h = fnvMix(h, mol.degree[i]);                                     // R&H #2: heavy atom degree
    int bondOrdSum = 0;
    for (int nb : mol.neighbors[i]) bondOrdSum += mol.bondOrder(i, nb);
    h = fnvMix(h, static_cast<uint64_t>(bondOrdSum));                 // R&H #3: bond order sum (valence)
    int massNum = (static_cast<int>(mol.massNumber.size()) > i)
                  ? mol.massNumber[i] : 0;
    h = fnvMix(h, static_cast<uint64_t>(massNum));                    // R&H #4: atomic mass number
    h = fnvMix(h, static_cast<uint64_t>(int64_t(mol.formalCharge[i]) + 4)); // R&H #5: formal charge
    const int64_t hCount = mol.attachedHydrogenCount(i);
    h = fnvMix(h, static_cast<uint64_t>(hCount));                     // R&H #6: attached H count
    h = fnvMix(h, mol.ring[i] ? 1 : 0);                              // R&H #7: ring membership
    h = fnvMix(h, mol.aromatic[i] ? 1 : 0);                          // Daylight extension: aromaticity
    if (!mol.tautomerClass.empty() && mol.tautomerClass[i] >= 0)
        h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(mol.tautomerClass[i]));
    return h;
}

/// Initial FCFP atom invariant: the pharmacophore feature class, degree and ring membership.
inline uint64_t pharmacophoreInvariant(const MolGraph& mol, const std::vector<int>& features, int i) {
    uint64_t h = FNV1A_SEED;
    h = fnvMix(h, static_cast<uint64_t>(features[i]));
    h = fnvMix(h, static_cast<uint64_t>(mol.degree[i]));
    h = fnvMix(h, mol.ring[i] ? 1ULL : 0ULL);
    return h;
}

/**
 * The Morgan iteration: every atom starts from invariant(atom), and at each radius takes the
 * hash of its radius, its previous hash and the sorted neighbour hashes (each mixed with its
 * bond order), so the result is invariant under atom renumbering. Every hash produced, at
 * radius 0 and beyond, is handed to sink.
 *
 * @param radius  Maximum radius; -1 iterates to the largest component extent.
 */
template <class Invariant, class Sink>
inline void morganHashes(const MolGraph& mol, int radius, Invariant&& invariant, Sink&& sink) {
    const int n = mol.n;
    std::vector<uint64_t> prevHash(n);
    for (int i = 0; i < n; ++i) {
        prevHash[i] = invariant(i);
        sink(prevHash[i]);
    }
    const int maxRadius = radius < 0 ? maximumComponentExtent(mol) : radius;
    std::vector<uint64_t> nextHash(n);
    std::vector<uint64_t> nbHashes;
    for (int r = 1; r <= maxRadius; ++r) {
        for (int i = 0; i < n; ++i) {
            const int deg = static_cast<int>(mol.neighbors[i].size());
            if (deg == 0) { nextHash[i] = prevHash[i]; continue; } // isolated atom: no expansion
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                const int nb = mol.neighbors[i][k];
                nbHashes[k] = prevHash[nb] ^ (static_cast<uint64_t>(mol.bondOrder(i, nb)) * 31ULL);
            }
            std::sort(nbHashes.begin(), nbHashes.end());

            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[i]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            nextHash[i] = h;
            sink(h);
        }
        std::swap(prevHash, nextHash);
    }
}

/// Folds every Morgan hash into a bit of an fpSize-bit fingerprint.
template <class Invariant>
inline std::vector<uint64_t> morganBits(const MolGraph& mol, int radius, int fpSize, Invariant&& invariant) {
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<uint64_t> fp((fpSize - 1) / 64 + 1, 0ULL);
    if (mol.n == 0) return fp;
    morganHashes(mol, radius, invariant, [&](uint64_t h) {
        path_fp_detail::setBit(fp, static_cast<int>(h % static_cast<uint64_t>(fpSize)));
    });
    return fp;
}

/// Counts the Morgan hashes falling into each of fpSize bins.
template <class Invariant>
inline std::vector<int> morganCounts(const MolGraph& mol, int radius, int fpSize, Invariant&& invariant) {
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    if (mol.n == 0) return counts;
    morganHashes(mol, radius, invariant, [&](uint64_t h) {
        counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
    });
    return counts;
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
    return morganBits(mol, radius, fpSize, [&](int atom) { return rogersHahnInvariant(mol, atom); });
}

/// Fingerprint invariant mode.
enum class FingerprintMode { ECFP, FCFP };

/// Classify atom into pharmacophoric feature classes (FCFP invariant).
/// Returns bitmask: bit 0=donor, 1=acceptor, 2=positive, 3=negative, 4=aromatic, 5=hydrophobic.
inline int classifyPharmacophore(const MolGraph& graph, int atom) {
    return graph.getPharmacophoreFeatures().at(atom);
}

/// FCFP-style circular fingerprint using pharmacophoric feature classes.
inline std::vector<uint64_t> computeCircularFingerprintFCFP(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    if (mol.n == 0) return std::vector<uint64_t>((fpSize - 1) / 64 + 1, 0ULL);
    const auto& features = mol.getPharmacophoreFeatures();
    return morganBits(mol, radius, fpSize, [&](int atom) { return pharmacophoreInvariant(mol, features, atom); });
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
inline std::vector<int> computeCircularFingerprintECFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    return morganCounts(mol, radius, fpSize, [&](int atom) { return rogersHahnInvariant(mol, atom); });
}

/// FCFP count fingerprint using pharmacophoric feature classes. Each element of
/// the returned vector holds the number of substructure hashes mapped to that
/// bin position.
///
/// @param radius  Maximum radius (-1 = unlimited, 2 = FCFP4, 3 = FCFP6)
/// @param fpSize  Number of counter bins (default 2048)
inline std::vector<int> computeCircularFingerprintFCFPCounts(
    const MolGraph& mol, int radius, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    if (mol.n == 0) return std::vector<int>(fpSize, 0);
    const auto& features = mol.getPharmacophoreFeatures();
    return morganCounts(mol, radius, fpSize, [&](int atom) { return pharmacophoreInvariant(mol, features, atom); });
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

} // namespace mcs_fp_detail

/// MCS-aware path fingerprint: encodes element, ring, aromatic, tautomer class,
/// degree per atom; bond order, ring, aromatic per bond. Paths are hashed
/// canonically (direction-independent). Compatible with the Java SMSD MCS fingerprint.
///
/// @param pathLength  Maximum path length in bonds (default 7)
/// @param fpSize      Fingerprint size in bits (default 2048)
inline std::vector<uint64_t> computeMCSFingerprintForChemistry(
    const MolGraph& mol, int pathLength, int fpSize, const ChemOptions& chemistry)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    if (pathLength <= 0 || pathLength == std::numeric_limits<int>::max())
        throw std::invalid_argument("pathLength must be positive and representable");
    const bool tautAware = chemistry.tautomerAware;
    // Atom and bond terms depend only on the graph, so compute each once.
    const auto adjacency = path_fp_detail::pathAdjacency(mol,
        [&](int atom) { return mol.atomicNum[atom] != 1; },
        [&](int from, int to) { return mcs_fp_detail::bondHash(mol, from, to); });
    std::vector<uint32_t> atomTerm(mol.n, 0);
    for (int i = 0; i < mol.n; ++i)
        if (adjacency.kept[i]) atomTerm[i] = mcs_fp_detail::atomHash(mol, tautAware, i);
    return path_fp_detail::pathFingerprint<path_fp_detail::LabelledPathHash>(adjacency, atomTerm,
                                                                             pathLength, fpSize);
}

inline std::vector<uint64_t> computeMCSFingerprint(
    const MolGraph& mol, int pathLength, int fpSize)
{
    return computeMCSFingerprintForChemistry(mol, pathLength, fpSize, ChemOptions{});
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
    if (g.aromatic[idx]) ++nPi;
    int inRing = g.ring[idx] ? 1 : 0;
    return (z << 16) | (heavyDeg << 8) | (nPi << 4) | inRing;
}

/**
 * Enumerates every 4-atom linear path A-B-C-D once and hands its canonical hash (the
 * lesser of the forward and reverse FNV-1a hashes of the atom types) to sink.
 */
template <class Sink>
inline void torsionHashes(const MolGraph& mol, Sink&& sink) {
    const int n = mol.n;
    if (n < 4) return;
    std::vector<int> atomType(n);
    for (int i = 0; i < n; ++i)
        atomType[i] = ttAtomType(mol, i);

    // The canonical hash is the same in both directions, so each path is hashed once,
    // from its lower end atom.
    for (int a = 0; a < n; ++a) {
        for (int b : mol.neighbors[a]) {
            for (int c : mol.neighbors[b]) {
                if (c == a) continue;
                for (int d : mol.neighbors[c]) {
                    if (d == b || d <= a) continue;
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
                    sink(fwd <= rev ? fwd : rev);
                }
            }
        }
    }
}

/// Topological torsion fingerprint (binary). Enumerates all 4-atom linear paths
/// A-B-C-D, hashes the atom types using FNV-1a with canonical ordering
/// (min of forward and reverse hash), and folds into fpSize bits.
///
/// Reference: Nilakantan et al., J. Chem. Inf. Comput. Sci. 1987, 27, 82-85.
///
/// @param fpSize  Fingerprint size in bits
inline std::vector<uint64_t> computeTopologicalTorsion(
    const MolGraph& mol, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<uint64_t> fp((fpSize - 1) / 64 + 1, 0ULL);
    torsionHashes(mol, [&](uint64_t h) {
        path_fp_detail::setBit(fp, static_cast<int>(h % static_cast<uint64_t>(fpSize)));
    });
    return fp;
}

/// Topological torsion count fingerprint. Each element of the returned vector
/// holds the number of distinct 4-atom torsion paths that hash to that bin.
///
/// @param fpSize  Number of counter bins
inline std::vector<int> computeTopologicalTorsionCounts(
    const MolGraph& mol, int fpSize)
{
    if (fpSize <= 0) throw std::invalid_argument("fpSize must be positive");
    std::vector<int> counts(fpSize, 0);
    torsionHashes(mol, [&](uint64_t h) {
        counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++;
    });
    return counts;
}

} // namespace detail

// ============================================================================
// Batch Fingerprint Generation
// ============================================================================

namespace detail {
/// Workers for a loop of count items when each worker should receive at least minimumPerWorker
/// of them; small loops then run on the calling thread instead of waking a team.
inline int workersFor(int count, int numThreads, int minimumPerWorker) {
    return std::max(1, std::min(resolveThreads(numThreads), count / std::max(1, minimumPerWorker)));
}
} // namespace detail

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
    if (fpSize <= 0 || fpSize > std::numeric_limits<int>::max() - 63)
        throw std::invalid_argument("fpSize must be positive and representable after rounding");
    if (pathLength <= 0) throw std::invalid_argument("pathLength must be positive");

    const int N = static_cast<int>(mols.size());
    std::vector<std::vector<uint64_t>> results(N);

    if (N == 0) return results;

    // Ensure fpSize is a multiple of 64
    fpSize = ((fpSize + 63) / 64) * 64;

    const int nThreads = detail::workersFor(N, numThreads, 1);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 8) num_threads(nThreads) if(nThreads > 1)
#endif
    for (int i = 0; i < N; ++i) {
        results[i] = detail::computePathFingerprint(mols[i], pathLength, fpSize);
    }

    return results;
}

// ============================================================================
// Fingerprint similarity: the one home of each metric and its empty-set rule
// ============================================================================

/// Set-bit counts of two fingerprints: bits in both, bits in the first and bits in the second.
struct BitCounts {
    int64_t common = 0, first = 0, second = 0;
};

/// Counts the set bits of two word fingerprints; words missing from the shorter one are empty.
inline BitCounts bitCounts(const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
    BitCounts counts;
    const size_t shared = std::min(a.size(), b.size());
    for (size_t i = 0; i < shared; ++i) {
        counts.common += popcount64(a[i] & b[i]);
        counts.first  += popcount64(a[i]);
        counts.second += popcount64(b[i]);
    }
    for (size_t i = shared; i < a.size(); ++i) counts.first += popcount64(a[i]);
    for (size_t i = shared; i < b.size(); ++i) counts.second += popcount64(b[i]);
    return counts;
}

/// Counts two sets of bit positions; a repeated position counts once.
inline BitCounts bitCounts(std::vector<int> a, std::vector<int> b) {
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    std::sort(b.begin(), b.end());
    b.erase(std::unique(b.begin(), b.end()), b.end());
    BitCounts counts;
    counts.first = static_cast<int64_t>(a.size());
    counts.second = static_cast<int64_t>(b.size());
    for (size_t i = 0, j = 0; i < a.size() && j < b.size();) {
        if (a[i] < b[j]) ++i;
        else if (b[j] < a[i]) ++j;
        else { ++counts.common; ++i; ++j; }
    }
    return counts;
}

/// Tanimoto |A & B| / |A | B|.  Two empty fingerprints are identical sets and score 1.
inline double tanimoto(const BitCounts& counts) {
    const int64_t either = counts.first + counts.second - counts.common;
    return either == 0 ? 1.0 : static_cast<double>(counts.common) / static_cast<double>(either);
}

/// Overlap |A & B| / min(|A|, |B|).  Two empty fingerprints score 1; one empty fingerprint scores 0.
inline double overlap(const BitCounts& counts) {
    if (counts.first == 0 && counts.second == 0) return 1.0;
    const int64_t smaller = std::min(counts.first, counts.second);
    return smaller == 0 ? 0.0 : static_cast<double>(counts.common) / static_cast<double>(smaller);
}

/// Dice 2|A & B| / (|A| + |B|).  Two empty fingerprints share no feature to weight and score 0.
inline double dice(const BitCounts& counts) {
    const int64_t sum = counts.first + counts.second;
    return sum == 0 ? 0.0 : (2.0 * static_cast<double>(counts.common)) / static_cast<double>(sum);
}

/// Cosine |A & B| / sqrt(|A| |B|).  An empty fingerprint has no direction and scores 0.
inline double cosine(const BitCounts& counts) {
    const double denominator = std::sqrt(static_cast<double>(counts.first) * static_cast<double>(counts.second));
    return denominator == 0.0 ? 0.0 : static_cast<double>(counts.common) / denominator;
}

/// Soergel distance 1 - Tanimoto, so two empty fingerprints are at distance 0.
inline double soergel(const BitCounts& counts) {
    return 1.0 - tanimoto(counts);
}

/// Tanimoto coefficient between two fingerprints; see tanimoto(const BitCounts&).
inline double fingerprintTanimoto(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    return tanimoto(bitCounts(a, b));
}

/// Overlap coefficient between two fingerprints; see overlap(const BitCounts&).
inline double fingerprintOverlap(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    return overlap(bitCounts(a, b));
}

/// Compute Dice coefficient between two fingerprints: 2*|A&B| / (|A| + |B|).
/// Dice weights shared features more heavily than Tanimoto, suitable for
/// scaffold-hopping searches and count-based fingerprint comparisons.
inline double fingerprintDice(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    return dice(bitCounts(a, b));
}

/// Compute Cosine similarity between two fingerprints: |A&B| / sqrt(|A| * |B|).
/// Normalises for fingerprint density differences between molecules of varying size.
inline double fingerprintCosine(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    return cosine(bitCounts(a, b));
}

/// Compute Soergel distance between two fingerprints: 1 - Tanimoto.
/// A proper metric distance useful for clustering and k-NN algorithms.
inline double fingerprintSoergel(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b)
{
    return soergel(bitCounts(a, b));
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
    // A subset test takes nanoseconds, so claims are large and short screens stay serial.
    const int nThreads = detail::workersFor(N, numThreads, 2048);
    (void)nThreads;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1024) num_threads(nThreads) if(nThreads > 1)
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
///   auto matches = Substructure::find(query).matchMany(corpus.targets());  // reuse many times
///   auto idxs = corpus.screen(query, 0.7);
class TargetCorpus {
    std::vector<MolGraph> targets_;
    std::vector<std::vector<uint64_t>> fingerprints_;
    bool prewarmed_ = false;

    std::shared_mutex& corpusMutex() const {
        static std::array<std::shared_mutex, 64> mutexes;
        return mutexes[(reinterpret_cast<uintptr_t>(this) >> 4) % mutexes.size()];
    }

public:
    /// Add a single target molecule.
    void addTarget(MolGraph g) {
        const std::unique_lock<std::shared_mutex> corpusLock(corpusMutex());
        targets_.push_back(std::move(g));
        prewarmed_ = false;
    }

    /// Add multiple target molecules.
    void addTargets(std::vector<MolGraph> gs) {
        const std::unique_lock<std::shared_mutex> corpusLock(corpusMutex());
        targets_.reserve(targets_.size() + gs.size());
        for (auto& g : gs) {
            targets_.push_back(std::move(g));
        }
        prewarmed_ = false;
    }

    /// Pre-compute lazy graph invariants and ECFP fingerprints for all targets.
    /// Query methods may run concurrently; additions and prewarming are exclusive.
    void prewarm(int numThreads = 0) {
        const std::unique_lock<std::shared_mutex> corpusLock(corpusMutex());
        const int N = static_cast<int>(targets_.size());
        if (N == 0) { prewarmed_ = true; return; }

        const int nThreads = detail::workersFor(N, numThreads, 1);
        (void)nThreads;

        // Prewarm graph invariants in parallel
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 1) num_threads(nThreads) if(nThreads > 1)
#endif
        for (int i = 0; i < N; ++i) {
            detail::prewarmGraph(targets_[i]);
        }

        // Compute ECFP fingerprints in parallel
        fingerprints_.resize(N);
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 8) num_threads(nThreads) if(nThreads > 1)
#endif
        for (int i = 0; i < N; ++i) {
            fingerprints_[i] = detail::computeCircularFingerprintECFP(
                targets_[i], /*radius=*/2, /*fpSize=*/2048);
        }

        prewarmed_ = true;
    }

    /// Number of targets in the corpus.
    size_t size() const {
        const std::shared_lock<std::shared_mutex> corpusLock(corpusMutex());
        return targets_.size();
    }

    /// Whether prewarm() has been called since the last addTarget/addTargets.
    bool isPrewarmed() const {
        const std::shared_lock<std::shared_mutex> corpusLock(corpusMutex());
        return prewarmed_;
    }

    /** Snapshot original target graphs for immutable prepared queries. */
    std::vector<MolGraph> targets() const {
        const std::shared_lock<std::shared_mutex> corpusLock(corpusMutex());
        std::vector<MolGraph> result;
        result.reserve(targets_.size());
        for (const auto& target : targets_) {
            const std::lock_guard<std::recursive_mutex> graphLock(target.lazyCacheMutex());
            result.push_back(target);
        }
        return result;
    }

    /// Fingerprint-based Tanimoto screening: returns indices of targets
    /// whose Tanimoto similarity to the query exceeds the threshold.
    /// Requires prewarm() to have been called (for fingerprints).
    std::vector<int> screen(
        const MolGraph& query, double threshold, int nThreads = 0) const
    {
        const std::shared_lock<std::shared_mutex> corpusLock(corpusMutex());
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
        // Each comparison is short, so claims are large and small corpora stay serial.
        const int nt = detail::workersFor(N, nThreads, 1024);
        (void)nt;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 512) num_threads(nt) if(nt > 1)
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

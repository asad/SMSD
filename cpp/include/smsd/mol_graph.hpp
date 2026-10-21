/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph, the header-only C++17 molecular graph that the searches, parsers and
 * fingerprints share. Zero external dependencies -- pure C++17 standard library.
 *
 * This header declares MolGraph and keeps its data and short accessors. Member
 * definitions live in purpose headers under detail/, included at the end:
 *   chem_options.hpp               ChemOptions and the aromaticity model
 *   detail/aromaticity.hpp         aromaticity perception and Kekule assignment
 *   detail/ring_perception.hpp     ring sets, relevant cycles, ring families and systems
 *   detail/tautomers.hpp           tautomer classes, weights and confidence
 *   detail/canonical_labelling.hpp canonical ranks, labels, orbits and hash
 *   detail/atom_invariants.hpp     neighbour labels, features, hydrogens and prewarm()
 *   detail/mol_graph_builder.hpp   MolGraph::Builder::build()
 *   detail/graph_operations.hpp    subgraphs, scaffolds and atom ID translation
 *   detail/smiles_writer.hpp       canonical SMILES
 * Each of them includes this header, so any one of them may be included first.
 */
#pragma once
#ifndef SMSD_MOL_GRAPH_HPP
#define SMSD_MOL_GRAPH_HPP

#include "smsd/bitops.hpp"
#include "smsd/chem_options.hpp"
#include "smsd/detail/union_find.hpp"
#include "smsd/general_matching.hpp"
#include "smsd/smiles_utils.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Branch prediction hint: dense bond storage is the common path.
#if defined(__GNUC__) || defined(__clang__)
#define SMSD_LIKELY(x) __builtin_expect(!!(x), 1)
#else
#define SMSD_LIKELY(x) (x)
#endif

namespace smsd {

// ============================================================================
// MolGraph -- immutable molecular graph representation
//
// Thread safety: lazy getters support concurrent readers. Copies and explicit
// mutations must not overlap a lazy fill or another mutation.
// ============================================================================

struct MolGraph {

    /** Synchronise lazy fills and snapshots of the cached public atom arrays. */
    std::recursive_mutex& lazyCacheMutex() const {
        static std::array<std::recursive_mutex, 64> mutexes;
        return mutexes[(reinterpret_cast<uintptr_t>(this) >> 4) % mutexes.size()];
    }

    static uint64_t nextInstanceNonce() {
        static std::atomic<uint64_t> counter{1};
        return counter.fetch_add(1, std::memory_order_relaxed);
    }

    static constexpr int SPARSE_THRESHOLD   = 200;
    static constexpr int HASH_PRIME         = 1000003;
    static constexpr int MAX_SEARCH_NODES   = 50000;
    static constexpr int CANON_SEARCH_LIMIT = 5000;

    // --- Core atom arrays ---
    int n     = 0;   // atom count
    int words = 0;   // number of uint64_t words for bit-parallel adjacency
    uint64_t cacheNonce = nextInstanceNonce();

    std::vector<int>  atomicNum;
    std::vector<int>  formalCharge;
    std::vector<int>  massNumber;
    std::vector<int>  hydrogenCount;
    bool hydrogenCountsKnown = false;
    std::vector<int>  atomClass;
    std::vector<int>  label;
    std::vector<uint8_t> ring;
    std::vector<uint8_t> aromatic;
    std::vector<int>  degree;
    mutable std::vector<int>  ringCount;
    mutable bool              ringCountsComputed_ = false;

    std::vector<std::vector<int>> neighbors;

    /// Optional external atom IDs for translating between SmsdGraph 0-based
    /// contiguous indices and non-contiguous IDs used by editors (e.g. BIME).
    /// When empty, internal indices are used directly (identity mapping).
    /// When populated, must have size == n. atomId[i] is the external ID of atom i.
    std::vector<int> atomId;

    // --- Molecule-level metadata ---
    std::string name;
    std::string programLine;
    std::string comment;
    std::map<std::string, std::string> properties;

    // --- Morgan / canonical / orbit (lazy — call ensureCanonical() before accessing) ---
    mutable std::vector<int>  morganRank;
    mutable std::vector<int>  canonicalLabel;
    mutable std::vector<int>  orbit;
    mutable std::vector<std::vector<int>> autGenerators_;   // automorphism generators (permutations)
    mutable bool              autGeneratorsTruncated_ = false;
    mutable uint64_t          canonicalHash = 0;  // uint64_t: wrapping is defined (no UB)
    mutable bool              canonicalComputed_ = false;

    // --- Deprecated pattern fingerprint; nothing reads it ---
    static constexpr int FP_WORDS = 4;  // 4 × 64 = 256 bits
    mutable uint64_t patternFP_[FP_WORDS] = {};
    mutable bool     patternFPComputed_ = false;

    [[deprecated("Nothing reads this fingerprint, and prewarm() does not compute it.")]]
    void ensurePatternFP() const;

    // --- Tautomer classes and per-atom pKa-informed relevance weights ---
    /**
     * The tautomer class of each atom, -1 outside every class. The atoms of one motif match share a
     * class and matches that share an atom join their classes, so the partition does not depend on
     * the atom order; the numbers, from 0, depend on it only for classes that Weisfeiler-Leman
     * refinement cannot tell apart (in practice those a symmetry of the graph exchanges).
     */
    std::vector<int>   tautomerClass;
    std::vector<float> tautomerWeight;  // 1.0 = non-tautomeric; <1.0 = relevance in the graph's solvent and pH
    /**
     * The ionisable site of each atom's tautomer motif, which adjustWeightsForPH() reads:
     * TAUTOMER_SITE_ACID when the neutral acid HA tautomerises, TAUTOMER_SITE_BASE when the
     * neutral base B does, and TAUTOMER_SITE_NONE when the library assigns the motif no pKa or
     * no match of the motif that covers the atom holds a hydrogen the site would release. The
     * motif of an atom in a class is the first, in a fixed priority order, whose pattern covers
     * the atom, so the weight, site and pKa do not depend on the order of the atoms. No motif
     * carries a base site.
     */
    std::vector<int8_t> tautomerSite;
    /** The aqueous pKa of that site (of the conjugate acid BH+ for a base); 0 when there is none. */
    std::vector<float>  tautomerPKa;
    ChemOptions::Solvent solvent_ = ChemOptions::Solvent::AQUEOUS;
    double pH_ = 7.4;

    /**
     * Scale the weight of every atom with an ionisable site by the neutral fraction of its class,
     * clamped to [0.1, 0.99]. With pKa_a the lowest acid pKa and pKa_b the highest base pKa among
     * the sites of the class, the fraction is 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH)), and a
     * class with sites of one kind drops the other term. Atoms without a site keep their weight
     * at every pH. Call after applySolventCorrection().
     */
    void adjustWeightsForPH(double pH);

    static constexpr int8_t TAUTOMER_SITE_NONE = 0;
    static constexpr int8_t TAUTOMER_SITE_ACID = 1;
    static constexpr int8_t TAUTOMER_SITE_BASE = -1;

    // Relevance weights of the dominant neutral tautomer (Sitzmann 2010); pH enters only through
    // the neutral fraction that adjustWeightsForPH() applies.
    static constexpr float TW_KETO_ENOL          = 0.98f; // keto form dominant (pKa ~20)
    static constexpr float TW_AMIDE_IMIDIC        = 0.97f; // amide dominant (pKa ~15-17)
    static constexpr float TW_LACTAM_LACTIM       = 0.97f; // lactam dominant
    static constexpr float TW_UREA                = 0.96f; // urea dominant (pKa ~14)
    static constexpr float TW_PYRIDINONE          = 0.95f; // lactam form dominant (N-H pKa ~11.7)
    static constexpr float TW_THIOAMIDE           = 0.94f; // thioamide dominant
    static constexpr float TW_THIONE_THIOL        = 0.94f; // thione dominant (pKa ~10-11)
    static constexpr float TW_ENAMINONE           = 0.85f; // vinylogous amide
    static constexpr float TW_HYDROXAMIC          = 0.85f; // hydroxamic acid (pKa ~9)
    static constexpr float TW_PHENOL_QUINONE      = 0.88f; // phenol dominant (pKa ~10)
    static constexpr float TW_HYDROXYPYRIMIDINE   = 0.90f; // hydroxy form dominant
    static constexpr float TW_PURINE_NH           = 0.78f; // purine N-H
    static constexpr float TW_DIKETONE_ENOL       = 0.72f; // 1,3-diketone enolisation: pH-sensitive, ~28% enol at pH 7.4 (pKa ~8-11)
    static constexpr float TW_NITROSO_OXIME       = 0.95f; // oxime strongly dominant over nitroso (Karabina 2014)
    static constexpr float TW_CYANAMIDE           = 0.80f; // cyanamide
    static constexpr float TW_IMIDE               = 0.92f; // imide NH (pKa ~9-10)
    static constexpr float TW_AMIDINE             = 0.50f; // symmetric tautomers
    static constexpr float TW_GUANIDINE           = 0.50f; // symmetric tautomers
    static constexpr float TW_IMIDAZOLE_NH        = 0.50f; // N1H/N3H symmetric
    static constexpr float TW_TRIAZOLE_NH         = 0.50f; // symmetric
    static constexpr float TW_TETRAZOLE_NH        = 0.50f; // 1H/2H symmetric
    static constexpr float TW_NITRO_ACI           = 0.95f; // nitro form strongly dominant at pH 7.4 (pKa ~3.5)

    // Tier 1b: additional transforms (T16-T30), Dhaked & Nicklaus 2024
    static constexpr float TW_THIOAMIDE_IMINOTHIOL = 0.95f; // thioamide form dominant (pKa ~13)
    static constexpr float TW_PHOSPHONATE          = 0.50f; // symmetric P(=O)(OH) ↔ P(OH)(=O)
    /// @deprecated Sulfoxide S=O is a covalent bond rather than a tautomeric form; no class uses this weight.
    [[deprecated("No tautomer class uses this weight.")]]
    static constexpr float TW_SULFOXIDE            = 0.0f;
    [[deprecated("Use TW_NITRO_ACI instead — unified nitro/aci-nitro weight.")]]
    static constexpr float TW_NITRO_ACI_NITRO      = TW_NITRO_ACI;
    /// @deprecated Nitrile and isonitrile are isomers rather than tautomers; no class uses this weight.
    [[deprecated("No tautomer class uses this weight.")]]
    static constexpr float TW_NITRILE_ISONITRILE   = 0.0f;
    static constexpr float TW_15_KETO_ENOL         = 0.75f; // 1,5-shift through conjugation
    static constexpr float TW_FURANOSE_PYRANOSE    = 0.60f; // ring-chain: furanose/pyranose
    static constexpr float TW_LACTOL               = 0.65f; // ring-chain: lactol
    static constexpr float TW_PHENOL_QUINONE_METHIDE = 0.92f; // phenol -> para-quinone methide
    static constexpr float TW_TRIAZOLE_NH_SHIFT    = 0.50f; // 1,2,3-triazole NH shift (symmetric)
    static constexpr float TW_BENZIMIDAZOLE_NH     = 0.50f; // NH-1 <-> NH-3 shift (symmetric)
    static constexpr float TW_PYRIDONE_HYDROXYPYRIDINE = TW_PYRIDINONE; // unified — both are pyridone/hydroxypyridine equilibria
    static constexpr float TW_BARBITURIC           = 0.80f; // tri-keto <-> enol forms
    static constexpr float TW_ALLYL_SHIFT          = 0.70f; // X-C=C-C <-> C=C-C-X (X=OH,SH,NH)
    static constexpr float TW_SELENOL_SELENOKETONE = 0.94f; // C=Se <-> C(-SeH)

    // Aqueous pKa of the ionisable site of each motif that has one, from the ranges above
    // (midpoints where a range is given), applied in every solvent. Every site is an acid, which a
    // match of its motif carries only when the atom that would release the proton bears a hydrogen.
    // The pKa ~3.5 of the nitro/aci-nitro motif belongs to the minor aci-nitro tautomer, not to the
    // dominant nitro form from which the motif arises, so that motif has no site.
    static constexpr float TPKA_KETO_ENOL            = 20.0f; // ketone alpha C-H
    static constexpr float TPKA_AMIDE_IMIDIC         = 16.0f; // amide N-H
    static constexpr float TPKA_UREA                 = 14.0f; // urea N-H
    static constexpr float TPKA_PYRIDINONE           = 11.7f; // 2-pyridone N-H (8.6 for pyrimidin-4(3H)-one)
    static constexpr float TPKA_PYRIDONE_HYDROXYPYRIDINE = TPKA_PYRIDINONE; // same equilibrium
    static constexpr float TPKA_THIOAMIDE            = 13.0f; // thioamide N-H
    static constexpr float TPKA_THIONE_THIOL         = 10.5f; // thione/thiol
    static constexpr float TPKA_HYDROXAMIC           = 9.0f;  // hydroxamic acid O-H
    static constexpr float TPKA_PHENOL_QUINONE       = 10.0f; // phenol O-H
    static constexpr float TPKA_DIKETONE_ENOL        = 9.5f;  // 1,3-diketone C-H
    static constexpr float TPKA_IMIDE                = 9.5f;  // imide N-H
    static constexpr float TPKA_THIOAMIDE_IMINOTHIOL = 13.0f; // thioamide N-H

    // --- Bit-parallel adjacency (n x words) ---
    std::vector<std::vector<uint64_t>> adjLong;

    // --- Bond storage: packed dense for small molecules, sparse for large ---
    // Dense mode: single uint8_t matrix packing order (bits 0-2), ring (bit 3), aromatic (bit 4).
    // Saves 66% memory vs 3 separate matrices.
    std::vector<std::vector<uint8_t>> bondPacked;    // dense n x n, 1 byte per pair
    bool useDense = false;
    static constexpr uint8_t BOND_ORDER_MASK = 0x07; // bits 0-2
    static constexpr uint8_t BOND_RING_BIT   = 0x08; // bit 3
    static constexpr uint8_t BOND_AROM_BIT   = 0x10; // bit 4

    // Sparse bond props: key = bondKey(i,j), value = {order, inRing, aromatic}
    std::unordered_map<int64_t, std::array<int,3>> sparseBondProps;
    std::unordered_map<int64_t, int> sparseDbStereo;

    // --- Stereochemistry ---
    std::vector<int>              tetraChirality;
    std::vector<std::vector<int>> dbStereoConf;   // dense, only for n <= SPARSE_THRESHOLD
    bool hasDbStereo = false;

    // --- SSSR rings (lazy, mutable) ---
    mutable std::vector<std::vector<int>> sssrRings;
    mutable bool sssrComputed = false;

    // --- Ring-system decomposition (Phase 2.3: symmetry breaking) ---
    mutable std::vector<int>      ringSystemId_;       // atom -> ring system ID (-1 if acyclic)
    mutable std::vector<uint64_t> ringSysSignature_;   // per ring system: hash signature
    mutable int                   ringSysCount_ = 0;   // number of distinct ring systems
    mutable bool                  ringSysComputed_ = false;

    // --- NLF caches (lazy, mutable for const-correct access) ---
    mutable std::vector<std::vector<int>> cachedNLF1, cachedNLF2, cachedNLF3;

    // --- Lazily computed per-atom caches ---
    mutable std::vector<std::vector<int>> cachedNeighborsByDegDesc;
    mutable std::vector<int> cachedPharmacophoreFeatures;
    mutable bool pharmComputed_ = false;

    // ========================================================================
    // Bond key: pack two indices into one int64
    // ========================================================================

    static int64_t bondKey(int i, int j) {
        int lo = std::min(i, j), hi = std::max(i, j);
        return (static_cast<int64_t>(lo) << 32)
             | static_cast<int64_t>(static_cast<uint32_t>(hi));
    }

    // ========================================================================
    // Bond accessors
    // ========================================================================

    int bondOrder(int i, int j) const {
        if (SMSD_LIKELY(useDense)) return bondPacked[i][j] & BOND_ORDER_MASK;
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() ? it->second[0] : 0;
    }

    bool bondInRing(int i, int j) const {
        if (SMSD_LIKELY(useDense)) return (bondPacked[i][j] & BOND_RING_BIT) != 0;
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() && it->second[1] != 0;
    }

    bool bondAromatic(int i, int j) const {
        if (SMSD_LIKELY(useDense)) return (bondPacked[i][j] & BOND_AROM_BIT) != 0;
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() && it->second[2] != 0;
    }

    bool hasBond(int i, int j) const { return bondOrder(i, j) != 0; }

    void invalidateBondDerivedState() {
        cacheNonce = nextInstanceNonce();
        canonicalComputed_ = false;
        canonicalHash = 0;
        morganRank.clear();
        canonicalLabel.clear();
        orbit.clear();
        autGenerators_.clear();
        autGeneratorsTruncated_ = false;
        ringSysComputed_ = false;
        ringSystemId_.clear();
        ringSysSignature_.clear();
        ringSysCount_ = 0;
        patternFPComputed_ = false;
        for (int w = 0; w < FP_WORDS; ++w) patternFP_[w] = 0;
        cachedPharmacophoreFeatures.clear();
        pharmComputed_ = false;
        cachedNLF1.clear();
        cachedNLF2.clear();
        cachedNLF3.clear();
        tautomerClass.clear();
        tautomerWeight.clear();
        tautomerSite.clear();
        tautomerPKa.clear();
        if (static_cast<int>(formalCharge.size()) == n
            && static_cast<int>(ring.size()) == n
            && static_cast<int>(aromatic.size()) == n) {
            computeTautomerClasses();
            applySolventCorrection(solvent_);
            adjustWeightsForPH(pH_);
        } else {
            tautomerClass.assign(n, -1);
            tautomerWeight.assign(n, 1.0f);
            tautomerSite.assign(n, TAUTOMER_SITE_NONE);
            tautomerPKa.assign(n, 0.0f);
        }
    }

    void setBondOrder(int i, int j, int ord) {
        if (useDense) {
            bondPacked[i][j] = (bondPacked[i][j] & ~BOND_ORDER_MASK) | static_cast<uint8_t>(ord & 0x07);
            bondPacked[j][i] = bondPacked[i][j];
        }
        else {
            auto key = bondKey(i, j);
            auto it = sparseBondProps.find(key);
            if (it != sparseBondProps.end()) it->second[0] = ord;
        }
        invalidateBondDerivedState();
    }

    void setBondAromatic(int i, int j, bool arom) {
        if (useDense) {
            if (arom) { bondPacked[i][j] |= BOND_AROM_BIT; bondPacked[j][i] |= BOND_AROM_BIT; }
            else      { bondPacked[i][j] &= ~BOND_AROM_BIT; bondPacked[j][i] &= ~BOND_AROM_BIT; }
        }
        else {
            auto key = bondKey(i, j);
            auto it = sparseBondProps.find(key);
            if (it != sparseBondProps.end()) it->second[2] = arom ? 1 : 0;
        }
        invalidateBondDerivedState();
    }

    int atomCount() const { return n; }

    int dbStereo(int i, int j) const {
        if (!hasDbStereo) return 0;
        if (useDense) return dbStereoConf[i][j];
        auto it = sparseDbStereo.find(bondKey(i, j));
        return it != sparseDbStereo.end() ? it->second : 0;
    }

    const std::vector<int>& getOrbit() const          { ensureCanonical(); return orbit; }
    const std::vector<int>& getCanonicalLabeling() const { ensureCanonical(); return canonicalLabel; }
    uint64_t getCanonicalHash() const                   { ensureCanonical(); return canonicalHash; }

    /// Return automorphism generators discovered during canonical labelling.
    /// Each generator is a permutation vector of length n where gen[i] is the
    /// image of atom i. Their closure is the full molecular automorphism group
    /// only when automorphismGeneratorsTruncated() is false.
    const std::vector<std::vector<int>>& getAutomorphismGenerators() const {
        ensureCanonical();
        return autGenerators_;
    }

    /// True when search limits or molecular-property checks prevent a complete
    /// molecular automorphism generator list. The orbit partition is advisory.
    bool automorphismGeneratorsTruncated() const {
        ensureCanonical();
        return autGeneratorsTruncated_;
    }

    // ========================================================================
    // Ring counts per atom (defined in detail/ring_perception.hpp)
    // ========================================================================

    /** Compute ring counts from a given ring set (SSSR or RCB). */
    void computeRingCounts(const std::vector<std::vector<int>>& rings) const;

    /** Default: compute ring counts from the relevant cycles of the ring candidate search over
     *  each atom's neighbours in index order, as the Java library counts them, so that the counts
     *  of cages such as adamantane and cubane do not depend on atom order. */
    void computeRingCounts() const;

    // ========================================================================
    // Aromaticity perception and Kekule assignment (detail/aromaticity.hpp)
    // ========================================================================

    static bool isAromaticDonorElement(int z);

    std::pair<bool, bool> exocyclicMultipleBondInfo(
        int atom,
        int prev = -1,
        int next = -1) const;

    static bool isSimpleCycle(const std::vector<int>& cycle);

    int aromaticPiContribution(const std::vector<int>& cycle, int pos) const;

    bool cycleLooksAromatic(const std::vector<int>& cycle) const;

    void enumerateRingCyclesDfs(
        int root,
        int current,
        const std::vector<std::vector<int>>& ringAdj,
        const std::vector<uint8_t>& inComponent,
        std::vector<uint8_t>& inPath,
        std::vector<int>& path,
        std::set<std::vector<int64_t>>& seenCycles,
        std::vector<std::vector<int>>& cycles,
        size_t maxCycles,
        uint64_t& steps) const;

    /**
     * Steps allowed to the simple-cycle walk that adds aromatic candidates beyond the relevant
     * cycles, counted in neighbour reads and candidate atoms. The walk also stops at 4096
     * cycles; in long fused ladders such as acenes it can otherwise branch exponentially before
     * reaching that count. Ordinary polycyclic molecules use a small fraction of it.
     */
    static constexpr uint64_t AROMATIC_CYCLE_STEP_LIMIT = uint64_t(1) << 24;

    std::vector<std::vector<int>> collectAromaticCandidateCycles() const;

    void clearRingFlagsFromTopology();

    void markRingEdge(int i, int j);

    void markAromaticEdge(int i, int j);

    void markKekuleEdge(int i, int j, int order);

    void refreshAtomLabels();

    void normalizeRingAndAromaticity(
        AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE);

    void perceiveAromaticity(
        AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE);

    int aromaticKekuleDemand(int atom) const;

    /** Determine bond assignments without changing the graph on failure. */
    bool findKekuleAssignment(const std::vector<int>& component,
                             std::vector<std::pair<int, int>>& singleBonds,
                             std::vector<std::pair<int, int>>& doubleBonds) const;

    bool kekulizeAromaticComponent(const std::vector<int>& component);

    bool kekulize();

    bool dearomatize();

    // ========================================================================
    // Ring perception (detail/ring_perception.hpp)
    // ========================================================================

    /**
     * Ring systems found from the bridges in linear time. Atoms joined by bonds that lie on a
     * cycle share a system; an atom whose bonds are all bridges has none (-1). Every cycle lies
     * inside one system, a bond is a bridge exactly when its atoms are not in the same system,
     * and a system with as many bonds as atoms is one simple ring. Both lists are empty when an
     * atom lists itself or a neighbour twice, and ring searches then read the whole graph.
     */
    struct CycleSystems {
        std::vector<int> ofAtom;
        std::vector<int> atoms;
        std::vector<uint8_t> simpleRing;
    };

    CycleSystems cycleSystems() const;

    // ========================================================================
    // Smallest set of smallest rings (shortest-path cycle basis)
    // ========================================================================

    const std::vector<std::vector<int>>& computeRings() const;

    // ========================================================================
    // Relevant cycles (union of all minimum cycle bases)
    // ========================================================================

private:

    friend std::vector<std::vector<std::vector<int>>> computeURFs(const MolGraph& graph);

    /** Group equal-length cycles whose shared-edge difference spans shorter cycles. */
    static std::vector<std::vector<std::vector<int>>> groupUniqueRingFamilies(
            const std::vector<std::vector<int>>& cycles,
            const std::unordered_map<int64_t, int>& edgeIndex);

    /** Find the lowest set bit position in a GF(2) vector. */
    static int findPivotBit(const std::vector<uint64_t>& vec);

    /** Reduce a GF(2) vector against a row echelon basis. Returns true if zero. */
    static bool reduceVectorGF2(std::vector<uint64_t>& vec,
                                const std::vector<std::vector<uint64_t>>& basis,
                                int basisSize);

    /** Add a GF(2) vector to row echelon basis if linearly independent. */
    static int addToBasisGF2(std::vector<uint64_t>& vec,
                             std::vector<std::vector<uint64_t>>& basis,
                             int basisSize);

    /** Canonical edge-set key for cycle deduplication, written into a reusable buffer. */
    static void canonicalCycleKeyInto(const std::vector<int>& cycle, std::vector<int64_t>& keys);

    /** Canonical edge-set key for cycle deduplication. */
    static std::vector<int64_t> canonicalCycleKey(const std::vector<int>& cycle);

    /** A parent list position on the explicit stack of shortest-path enumeration. */
    struct PathFrame { int node; size_t nextParent; };

    /**
     * All shortest paths from root to target in a BFS tree (at most 50), into reusable storage:
     * paths[0, count) receive them, and spare entries keep their capacity. Parents are followed
     * depth first in list order on an explicit stack, so a path as long as the molecule needs no
     * call stack; current holds the path from the target to the deepest node. Returns how many
     * paths were written.
     */
    static size_t allShortestPathsCpp(int root, int target, const std::vector<int>& dist,
                                      const std::vector<std::vector<int>>& parents, std::vector<int>& current,
                                      std::vector<PathFrame>& stack, std::vector<std::vector<int>>& paths);

public:

    /**
     * Work allowance of one relevant-cycle enumeration, counted in atoms read: breadth-first
     * visits, shortest-path atoms, overlap tests and candidate cycle atoms. Ordinary molecules,
     * macrocycles of a few thousand atoms included, use a small fraction of it. When an
     * enumeration exhausts it, which takes a fraction of a second, the remaining roots are not
     * searched and the result is the relevant cycles among the minimum cycle basis and the
     * candidates found so far, so ring perception stays bounded for very large ring systems.
     */
    static constexpr uint64_t RELEVANT_CYCLE_WORK_LIMIT = uint64_t(1) << 26;

    /**
     * Computes the set of relevant cycles (union of all minimum cycle bases).
     * A cycle is relevant if it cannot be expressed as the XOR (symmetric
     * difference) of strictly shorter cycles. The enumeration is bounded by
     * RELEVANT_CYCLE_WORK_LIMIT.
     */
    std::vector<std::vector<int>> computeRelevantCycles() const;

    /** Compute Unique Ring Families from relevant cycles and strictly shorter cycle spans. */
    std::vector<std::vector<std::vector<int>>> computeURFs();

    // ========================================================================
    // Tautomer class detection: pKa-informed transforms, Tier 1 (detail/tautomers.hpp)
    // ========================================================================

    void computeTautomerClasses();

    // ========================================================================
    // Tier 2: solvent-dependent tautomer weight corrections
    // ========================================================================

    /** Adjust tautomer weights for solvent effects. Call after computeTautomerClasses(). */
    void applySolventCorrection(ChemOptions::Solvent solvent);

    // ========================================================================
    // Neighbourhood rank computation (detail/canonical_labelling.hpp)
    // ========================================================================

    static std::vector<int> computeMorganRanks(
            int n, const std::vector<int>& initLabel,
            const std::vector<std::vector<int>>& neighbors);

    // ========================================================================
    // Orbits from the union-find of automorphic atoms
    // ========================================================================

    static std::vector<int> buildOrbits(std::vector<int>& uf, int n);

    // ========================================================================
    // 2-hop neighbourhood signature for orbit disambiguation
    // ========================================================================

    /**
     * Compute a deterministic 2-hop neighbourhood signature for atom v.
     *
     * The signature encodes the sorted multiset of (neighbor_label, neighbor_degree)
     * pairs for all atoms within 2 hops. Two atoms in the same automorphism orbit
     * must have identical 2-hop signatures, so splitting orbits by signature is safe.
     *
     * This refines orbits beyond what Morgan ranks alone achieve: for example,
     * in spiro compounds or fused ring systems, Morgan may converge prematurely
     * but 2-hop signatures can still distinguish structurally distinct positions.
     */
    static uint64_t twoHopSignature(
            int v, const std::vector<int>& label,
            const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors);

    /**
     * Refine an orbit partition using 2-hop neighbourhood signatures.
     *
     * Atoms that share the same orbit representative but have different 2-hop
     * signatures are split into separate orbits. This is a safe refinement:
     * true automorphism-equivalent atoms always have identical k-hop signatures
     * for any k.
     */
    static std::vector<int> refineOrbitsBy2Hop(
            const std::vector<int>& orb, int n,
            const std::vector<int>& label,
            const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors);

    // ========================================================================
    // Canonical labelling (partition refinement and individualisation)
    // ========================================================================

    static int canonicalBondLabel(const MolGraph* graph, int first, int second);

    static int canonicalTetrahedralTag(const MolGraph& graph, int atom,
                                        const std::vector<int>& ranking);

    static int canonicalDoubleBondTag(const MolGraph& graph, int first, int second,
                                       const std::vector<int>& ranking);

    std::vector<int> canonicalAtomLabels() const;

    static int findSmallestNonSingleton(int n, const std::vector<bool>& cellEnd);

    static void refinePartitionImpl(
            int n, std::vector<int>& perm, std::vector<bool>& cellEnd,
            const std::vector<std::vector<int>>& neighbors,
            const std::vector<int>& label,
            const std::vector<int>& degree, int /*targetPos*/, const MolGraph* graph);

    /// @deprecated Use getCanonicalLabeling() on the graph.
    [[deprecated("Use getCanonicalLabeling() on the graph.")]]
    static void refinePartition(
            int n, std::vector<int>& perm, std::vector<bool>& cellEnd,
            const std::vector<std::vector<int>>& neighbors,
            const std::vector<int>& label, const std::vector<int>& degree, int targetPos);

    static int comparePermImpl(
            const std::vector<int>& perm1, const std::vector<int>& perm2,
            int n, const std::vector<int>& label,
            const std::vector<std::vector<int>>& neighbors, const MolGraph* graph);

    /// @deprecated Use getCanonicalLabeling() on the graph.
    [[deprecated("Use getCanonicalLabeling() on the graph.")]]
    static int comparePerm(
            const std::vector<int>& perm1, const std::vector<int>& perm2,
            int n, const std::vector<int>& label, const std::vector<std::vector<int>>& neighbors);

    static int partialCompare(
            const std::vector<int>& perm, const std::vector<int>& bestPerm,
            int n, const std::vector<bool>& cellEnd,
            const std::vector<int>& label,
            const std::vector<std::vector<int>>& /*neighbours*/);

    struct CanonResult {
        std::vector<int> canonLabel;
        std::vector<int> orbit;
        std::vector<std::vector<int>> autGenerators;
        bool generatorsTruncated = false;
    };

    static CanonResult computeCanonicalLabelingImpl(
            int n, const std::vector<int>& label,
            const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors, const MolGraph* graph);

    /// @deprecated Use getCanonicalLabeling() on the graph.
    [[deprecated("Use getCanonicalLabeling() on the graph.")]]
    static CanonResult computeCanonicalLabeling(
            int n, const std::vector<int>& label, const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors);

    // ========================================================================
    // Canonical hash
    // ========================================================================

    static uint64_t computeCanonicalHash(
            const MolGraph& g, const std::vector<int>& canonLabel);

    // ========================================================================
    // Canonical SMILES generation (detail/smiles_writer.hpp)
    // ========================================================================

    /** Write canonical SMILES, retaining atom properties and stereochemistry. */
    std::string toCanonicalSmiles() const;

    // ========================================================================
    // NLF (Neighbour Label Frequency) helpers (detail/atom_invariants.hpp)
    // ========================================================================

    static std::vector<int> freqMapToSortedArray(
            const std::unordered_map<int,int>& freq);

    // NLF labels use atomicNum + aromaticity but NOT ring membership.
    // Ring compatibility is enforced separately in atom/bond checks so the
    // cached NLF tables remain reusable across option profiles.
    static int nlfLabel(const MolGraph& g, int idx);

    static std::vector<int> buildNLF1(const MolGraph& g, int idx);

    // Count distinct atoms within the radius, excluding the root. Exact
    // distance shells are unsafe for substructure pruning: additional target
    // bonds can shorten distances without invalidating a query embedding.
    static std::vector<int> buildRadiusNLF(const MolGraph& g, int idx, int radius);

    static std::vector<int> buildNLF2(const MolGraph& g, int idx);

    static std::vector<int> buildNLF3(const MolGraph& g, int idx);

    using NLFBuilder = std::vector<int>(*)(const MolGraph&, int);

    static std::vector<std::vector<int>> buildAllNLF(
            const MolGraph& g, NLFBuilder builder);

    /// @deprecated The searches do not read the neighbour-label tables.
    [[deprecated("The searches do not read the neighbour-label tables.")]]
    const std::vector<std::vector<int>>& getNLF1() const;

    /// @deprecated The searches do not read the neighbour-label tables.
    [[deprecated("The searches do not read the neighbour-label tables.")]]
    const std::vector<std::vector<int>>& getNLF2() const;

    /// @deprecated The searches do not read the neighbour-label tables.
    [[deprecated("The searches do not read the neighbour-label tables.")]]
    const std::vector<std::vector<int>>& getNLF3() const;

    static bool nlfOk(const std::vector<int>& fq, const std::vector<int>& ft);

    /// @deprecated The searches do not read the neighbour-label tables.
    [[deprecated("The searches do not read the neighbour-label tables.")]]
    static bool nlfCheckOk(
            int qi, int tj,
            const std::vector<std::vector<int>>& qNLF1,
            const std::vector<std::vector<int>>& tNLF1,
            const std::vector<std::vector<int>>& qNLF2,
            const std::vector<std::vector<int>>& tNLF2,
            const std::vector<std::vector<int>>& qNLF3,
            const std::vector<std::vector<int>>& tNLF3,
            bool useTwoHop, bool useThreeHop);

    // ========================================================================
    // Derived data and lazy caches: initDerivedFields() is defined in
    // detail/mol_graph_builder.hpp, each lazy getter with the property it fills
    // ========================================================================

    void initDerivedFields();

    /** Lazily compute Morgan ranks, canonical labelling, orbit, and canonical hash. */
    void ensureCanonical() const;

    /** Lazily compute per-atom ring counts. */
    void ensureRingCounts() const;

    /** Lazily compute ring-system decomposition via union-find on ring atoms.
     *  Two ring atoms belong to the same ring system if they share an edge where
     *  both endpoints are ring atoms (i.e. the edge is in a ring or connects two
     *  ring atoms in fused/spiro systems).
     *  Each ring system gets a canonical hash signature based on its member
     *  element types, bond orders, size, and aromatic atom count. */
    void ensureRingSystems() const;

    /** Ring system ID for an atom (-1 if acyclic). */
    int ringSystemOf(int atom) const;

    /** Hash signature for a ring system (by ID). Returns 0 for invalid IDs. */
    uint64_t ringSystemSig(int rsId) const;

    /** Number of distinct ring systems in this molecule. */
    int ringSystemCount() const { ensureRingSystems(); return ringSysCount_; }

    /**
     * The neighbours of each atom sorted by descending degree, computed on the first call and
     * cached with the graph.
     */
    [[deprecated("Nothing reads this cache, and prewarm() does not compute it.")]]
    const std::vector<std::vector<int>>& getNeighborsByDegDesc() const;

    /** Stored implicit H, or the SMILES valence model's count when hydrogen counts were not supplied. */
    int implicitHydrogenCount(int atom) const;

    /** Stored implicit H plus explicit hydrogen neighbours, with valence fallback when omitted. */
    int64_t attachedHydrogenCount(int atom) const;

    /**
     * Cached pharmacophore feature classification per atom for FCFP.
     * Computed once per MolGraph, so that FCFP does not repeat the O(n×degree²) classification.
     */
    const std::vector<int>& getPharmacophoreFeatures() const;

    /** Initialise the lazy properties used by matching before copying a shared graph. */
    void prewarm() const;

    // ========================================================================
    // Builder -- molecule construction via fluent API (build() is in
    // detail/mol_graph_builder.hpp)
    // ========================================================================

    class Builder {
    public:
        enum class RowEncoding : uint8_t { ABSENT = 0, PARALLEL = 1, DENSE = 2 };

        int n_ = 0;
        std::vector<int>  atomicNum_;
        std::vector<int>  formalCharge_;
        std::vector<int>  massNumber_;
        std::vector<int>  hydrogenCount_;
        std::vector<int>  atomClass_;
        std::vector<uint8_t> ring_;
        std::vector<uint8_t> aromatic_;
        std::vector<std::vector<int>>  neighbors_;
        std::vector<std::vector<int>>  bondOrders_;
        std::vector<std::vector<bool>> bondRings_;
        std::vector<std::vector<bool>> bondAroms_;
        std::vector<int>               tetraChirality_;
        std::vector<std::vector<int>>  dbStereoConf_;
        std::vector<std::tuple<int, int, int>> sparseDbStereoConf_;
        std::vector<int>               atomIds_;
        std::string                    name_;
        std::string                    programLine_;
        std::string                    comment_;
        std::map<std::string, std::string> properties_;
        ChemOptions::Solvent           solvent_ = ChemOptions::Solvent::AQUEOUS;
        double                         pH_ = 7.4;

        Builder& atomCount(int n)                                         { n_ = n; return *this; }
        Builder& atomicNumbers(std::vector<int> v)                        { atomicNum_ = std::move(v); return *this; }
        Builder& formalCharges(std::vector<int> v)                        { formalCharge_ = std::move(v); return *this; }
        Builder& massNumbers(std::vector<int> v)                          { massNumber_ = std::move(v); return *this; }
        Builder& hydrogenCounts(std::vector<int> v)                       { hydrogenCount_ = std::move(v); return *this; }
        Builder& atomClasses(std::vector<int> v)                          { atomClass_ = std::move(v); return *this; }
        Builder& ringFlags(std::vector<uint8_t> v)                          { ring_ = std::move(v); return *this; }
        Builder& aromaticFlags(std::vector<uint8_t> v)                    { aromatic_ = std::move(v); return *this; }
        Builder& setNeighbors(std::vector<std::vector<int>> v)            { neighbors_ = std::move(v); return *this; }
        Builder& setBondOrders(std::vector<std::vector<int>> v)           { bondOrders_ = std::move(v); return *this; }
        Builder& bondRingFlags(std::vector<std::vector<bool>> v)          { bondRings_ = std::move(v); return *this; }
        Builder& bondAromaticFlags(std::vector<std::vector<bool>> v)      { bondAroms_ = std::move(v); return *this; }
        /**
         * Tetrahedral tags: 0 for none, 1 for anticlockwise (@) and 2 for clockwise (@@). A tag
         * winds the neighbours in their order, with an implicit hydrogen or lone pair second when
         * the first neighbour has a lower index than the atom and first otherwise, see
         * stereo_detail::tetrahedralLigands().
         */
        Builder& tetrahedralChirality(std::vector<int> v)                 { tetraChirality_ = std::move(v); return *this; }
        Builder& doubleBondStereo(std::vector<std::vector<int>> v)        { dbStereoConf_ = std::move(v); sparseDbStereoConf_.clear(); return *this; }
        /** Supply (atom, atom, configuration) tuples; configuration is 1 (Z) or 2 (E). */
        Builder& doubleBondStereoSparse(std::vector<std::tuple<int, int, int>> v) {
            sparseDbStereoConf_ = std::move(v);
            dbStereoConf_.clear();
            return *this;
        }
        /** Optional external atom IDs (1-based, gaps allowed). Size must equal atomCount if provided. */
        Builder& atomIds(std::vector<int> v)                               { atomIds_ = std::move(v); return *this; }
        Builder& name(std::string v)                                        { name_ = std::move(v); return *this; }
        Builder& programLine(std::string v)                                 { programLine_ = std::move(v); return *this; }
        Builder& comment(std::string v)                                     { comment_ = std::move(v); return *this; }
        Builder& properties(std::map<std::string, std::string> v)           { properties_ = std::move(v); return *this; }
        Builder& solvent(ChemOptions::Solvent s)                            { solvent_ = s; return *this; }
        /** pH context for the built graph's baseline tautomer weights, clamped to [0, 14]. */
        Builder& pH(double v) {
            if (!std::isfinite(v)) throw std::invalid_argument("pH must be finite");
            pH_ = std::max(0.0, std::min(14.0, v));
            return *this;
        }

        MolGraph build(
            bool perceiveAromaticity = true,
            AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE) const;
    };
};

} // namespace smsd

#include "smsd/detail/aromaticity.hpp"
#include "smsd/detail/ring_perception.hpp"
#include "smsd/detail/tautomers.hpp"
#include "smsd/detail/canonical_labelling.hpp"
#include "smsd/detail/atom_invariants.hpp"
#include "smsd/detail/mol_graph_builder.hpp"
#include "smsd/detail/graph_operations.hpp"
#include "smsd/detail/smiles_writer.hpp"

#endif // SMSD_MOL_GRAPH_HPP

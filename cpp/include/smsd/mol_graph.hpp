/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only C++17 port of MolGraph / ChemOptions / ChemOps.
 * Zero external dependencies -- pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_MOL_GRAPH_HPP
#define SMSD_MOL_GRAPH_HPP

#include "smsd/bitops.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <deque>
#include <functional>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace smsd {

enum class AromaticityModel { DAYLIGHT_LIKE };

// ============================================================================
// ChemOptions -- configuration for chemical matching constraints
// ============================================================================

struct ChemOptions {
    enum class BondOrderMode   { STRICT, LOOSE, ANY };
    enum class AromaticityMode { STRICT, FLEXIBLE };
    enum class MatcherEngine   { VF2, VF2PP, VF3 };
    enum class RingFusionMode  { IGNORE, PERMISSIVE, STRICT };
    enum class Solvent         { AQUEOUS, DMSO, METHANOL, CHLOROFORM, ACETONITRILE, DIETHYL_ETHER };

    // Atom matching
    bool matchAtomType            = true;
    bool matchFormalCharge        = true;
    bool useChirality             = false;
    bool useBondStereo            = false;
    bool ringMatchesRingOnly      = true;

    // Ring completeness
    bool completeRingsOnly        = false;

    // Isotope matching
    bool matchIsotope             = false;

    // Ring fusion
    RingFusionMode ringFusionMode = RingFusionMode::IGNORE;

    // Tautomer matching
    bool   tautomerAware = false;
    /** Simulation pH for pKa-informed tautomer relevance scoring. Default 7.4 (physiological). */
    double pH            = 7.4;

    /** Set simulation pH, clamped to [0, 14]. */
    ChemOptions& withPH(double v) { pH = std::max(0.0, std::min(14.0, v)); return *this; }

    // Solvent for Tier 2 pKa-informed tautomer weight corrections
    Solvent solvent = Solvent::AQUEOUS;

    /** Set solvent for tautomer weight corrections (fluent). */
    ChemOptions& withSolvent(Solvent s) { solvent = s; return *this; }

    // Bond matching
    BondOrderMode  matchBondOrder  = BondOrderMode::STRICT;
    AromaticityMode aromaticityMode = AromaticityMode::FLEXIBLE;

    // Engine
    MatcherEngine matcherEngine   = MatcherEngine::VF2PP;

    // Induced subgraph isomorphism: when true, if query atoms qi and qk
    // have NO bond between them, the corresponding mapped target atoms
    // must also have NO bond. Default false (standard non-induced matching).
    bool induced = false;

    // Pruning
    bool useTwoHopNLF             = true;
    bool useThreeHopNLF           = false;
    bool useBitParallelFeasibility = true;

    // Fluent setters
    ChemOptions& withAromaticityMode(AromaticityMode m) { aromaticityMode = m; return *this; }
    ChemOptions& withBondOrderMode(BondOrderMode m)     { matchBondOrder = m; return *this; }
    ChemOptions& withTwoHopNLF(bool on)                 { useTwoHopNLF = on; return *this; }
    ChemOptions& withThreeHopNLF(bool on)               { useThreeHopNLF = on; return *this; }
    ChemOptions& withBitParallelFeasibility(bool on)     { useBitParallelFeasibility = on; return *this; }
    ChemOptions& withCompleteRingsOnly(bool on)          { completeRingsOnly = on; return *this; }
    ChemOptions& withMatchIsotope(bool on)               { matchIsotope = on; return *this; }
    ChemOptions& withRingFusionMode(RingFusionMode m)    { ringFusionMode = m; return *this; }
    ChemOptions& withTautomerAware(bool on)              { tautomerAware = on; return *this; }
    ChemOptions& withInduced(bool on)                    { induced = on; return *this; }

    // Named profiles
    static ChemOptions profile(const std::string& name) {
        ChemOptions opts;
        if (name == "strict") {
            opts.matchBondOrder  = BondOrderMode::STRICT;
            opts.aromaticityMode = AromaticityMode::STRICT;
        } else {
            // "compat-fmcs", "compat-substruct", default
            opts.matchBondOrder      = BondOrderMode::LOOSE;
            opts.aromaticityMode     = AromaticityMode::FLEXIBLE;
            opts.ringMatchesRingOnly = false;
        }
        return opts;
    }

    static ChemOptions tautomerProfile() {
        ChemOptions c;
        c.tautomerAware      = true;
        c.matchBondOrder      = BondOrderMode::LOOSE;
        c.aromaticityMode     = AromaticityMode::FLEXIBLE;
        c.matchFormalCharge   = false;
        c.ringMatchesRingOnly = false;
        return c;
    }

    /** Relaxed profile for loose FMCS-style matching — ring atoms may map to chain atoms,
     *  partial rings are accepted, bond orders and charges are not enforced.
     *  Produces numerically larger but topologically looser MCS than the strict default.
     *  Use for interoperability benchmarks or when loose topology matching is required. */
    static ChemOptions fmcsProfile() {
        ChemOptions c;
        c.matchBondOrder      = BondOrderMode::LOOSE;
        c.ringMatchesRingOnly = false;
        c.completeRingsOnly   = false;
        c.matchFormalCharge   = false;
        return c;
    }
};

// ============================================================================
// MolGraph -- immutable molecular graph representation
//
// Thread safety: MolGraph is safe to share across threads AFTER calling
// ensureCanonical() and ensureRingCounts() once from a single thread.
// Those methods lazily populate mutable cached fields and are NOT
// internally synchronised. Two threads calling them concurrently on the
// same instance will race.
// ============================================================================

struct MolGraph {

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

    // --- Pattern fingerprint for O(1) substructure pre-screening (v6.8.0) ---
    // 256-bit structural bitset encoding local atom-pair and path features.
    // Computed lazily; if (query_fp & target_fp) != query_fp, the target is
    // guaranteed to be missing a required feature → abort before VF2++.
    static constexpr int FP_WORDS = 4;  // 4 × 64 = 256 bits
    mutable uint64_t patternFP_[FP_WORDS] = {};
    mutable bool     patternFPComputed_ = false;

    void ensurePatternFP() const {
        if (patternFPComputed_) return;
        for (int w = 0; w < FP_WORDS; ++w) patternFP_[w] = 0;
        // Encode atom features: element + aromaticity → 1 bit per type
        for (int i = 0; i < n; ++i) {
            uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 31 + (aromatic[i] ? 1 : 0);
            patternFP_[0] |= uint64_t(1) << (h & 63);
        }
        // Encode bond-pair features: (Z1, bondOrder, Z2) for each edge
        for (int i = 0; i < n; ++i) {
            for (int j : neighbors[i]) {
                if (j <= i) continue;  // each edge once
                int bo = bondOrder(i, j);
                uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 997
                           + static_cast<uint32_t>(atomicNum[j]) * 31
                           + static_cast<uint32_t>(bo) * 7;
                patternFP_[1] |= uint64_t(1) << (h & 63);
                // Also encode reverse pair for asymmetric element pairs
                uint32_t h2 = static_cast<uint32_t>(atomicNum[j]) * 997
                            + static_cast<uint32_t>(atomicNum[i]) * 31
                            + static_cast<uint32_t>(bo) * 7;
                patternFP_[1] |= uint64_t(1) << (h2 & 63);
            }
        }
        // Encode 2-hop path features: (Z1, Z2, Z3) for paths of length 2
        for (int i = 0; i < n; ++i) {
            for (int j : neighbors[i]) {
                for (int k : neighbors[j]) {
                    if (k == i) continue;
                    uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 997
                               + static_cast<uint32_t>(atomicNum[j]) * 31
                               + static_cast<uint32_t>(atomicNum[k]);
                    patternFP_[2] |= uint64_t(1) << (h & 63);
                }
            }
        }
        // Encode ring membership + degree features
        for (int i = 0; i < n; ++i) {
            uint32_t h = static_cast<uint32_t>(atomicNum[i]) * 31
                       + static_cast<uint32_t>(degree[i]) * 7
                       + (ring[i] ? 3u : 0u);
            patternFP_[3] |= uint64_t(1) << (h & 63);
        }
        patternFPComputed_ = true;
    }

    // --- Tautomer classes and per-atom pKa-informed relevance weights ---
    std::vector<int>   tautomerClass;   // -1 = not tautomeric
    std::vector<float> tautomerWeight;  // 1.0 = non-tautomeric; <1.0 = relevance at pH 7.4
    ChemOptions::Solvent solvent_ = ChemOptions::Solvent::AQUEOUS;
    double pH_ = 7.4;

    /** Adjust tautomer weights for pH deviation from 7.4 (Henderson-Hasselbalch scaling). */
    void adjustWeightsForPH(double pH) {
        if (std::abs(pH - 7.4) < 0.01) return;
        // Dampen weights toward 0.5 as pH deviates from 7.4 (never cross 0.5)
        float dampen = static_cast<float>(0.1 * std::abs(pH - 7.4));
        for (int i = 0; i < n; ++i) {
            if (tautomerClass[i] == -1 || tautomerWeight[i] == 1.0f) continue;
            float w = tautomerWeight[i];
            float shift = (w > 0.5f) ? -dampen : dampen;
            float adjusted = w + shift;
            if (w > 0.5f) adjusted = std::max(0.5f, adjusted);
            else           adjusted = std::min(0.5f, adjusted);
            tautomerWeight[i] = std::max(0.1f, std::min(0.99f, adjusted));
        }
    }

    // pKa-informed relevance weights at pH 7.4 (Sitzmann 2010)
    static constexpr float TW_KETO_ENOL          = 0.98f; // keto form dominant (pKa ~20)
    static constexpr float TW_AMIDE_IMIDIC        = 0.97f; // amide dominant (pKa ~25)
    static constexpr float TW_LACTAM_LACTIM       = 0.97f; // lactam dominant
    static constexpr float TW_UREA                = 0.96f; // urea dominant (pKa ~14)
    static constexpr float TW_PYRIDINONE          = 0.95f; // lactam form dominant (pKa ~1)
    static constexpr float TW_THIOAMIDE           = 0.94f; // thioamide dominant
    static constexpr float TW_THIONE_THIOL        = 0.94f; // thione dominant (pKa ~10-11)
    static constexpr float TW_ENAMINONE           = 0.85f; // vinylogous amide
    static constexpr float TW_HYDROXAMIC          = 0.85f; // hydroxamic acid (pKa ~8)
    static constexpr float TW_PHENOL_QUINONE      = 0.88f; // phenol dominant (pKa ~10)
    static constexpr float TW_HYDROXYPYRIMIDINE   = 0.90f; // hydroxy form dominant
    static constexpr float TW_PURINE_NH           = 0.78f; // purine N-H
    static constexpr float TW_DIKETONE_ENOL       = 0.72f; // pH-sensitive at 7.4 (pKa ~8-11)
    static constexpr float TW_NITROSO_OXIME       = 0.70f; // oxime slightly more stable
    static constexpr float TW_CYANAMIDE           = 0.80f; // cyanamide
    static constexpr float TW_IMIDE               = 0.92f; // imide NH (pKa ~9-10)
    static constexpr float TW_AMIDINE             = 0.50f; // symmetric tautomers
    static constexpr float TW_GUANIDINE           = 0.50f; // symmetric tautomers
    static constexpr float TW_IMIDAZOLE_NH        = 0.50f; // N1H/N3H symmetric
    static constexpr float TW_TRIAZOLE_NH         = 0.50f; // symmetric
    static constexpr float TW_TETRAZOLE_NH        = 0.50f; // 1H/2H symmetric
    static constexpr float TW_NITRO_ACI           = 0.25f; // aci-nitro minor at pH 7.4

    // Tier 1b: additional transforms (T16-T30), Dhaked & Nicklaus 2024
    static constexpr float TW_THIOAMIDE_IMINOTHIOL = 0.95f; // thioamide form dominant (pKa ~13)
    static constexpr float TW_PHOSPHONATE          = 0.50f; // symmetric P(=O)(OH) ↔ P(OH)(=O)
    // T18 DISABLED: sulfoxide S=O is a true covalent bond, not tautomerism
    static constexpr float TW_SULFOXIDE            = 0.0f;  // DISABLED
    static constexpr float TW_NITRO_ACI_NITRO      = 0.95f; // nitro form strongly dominant
    // T20 DISABLED: nitrile/isonitrile is isomerism, not tautomerism (no proton shift)
    static constexpr float TW_NITRILE_ISONITRILE   = 0.0f;  // DISABLED
    static constexpr float TW_15_KETO_ENOL         = 0.75f; // 1,5-shift through conjugation
    static constexpr float TW_FURANOSE_PYRANOSE    = 0.60f; // ring-chain: furanose/pyranose
    static constexpr float TW_LACTOL               = 0.65f; // ring-chain: lactol
    static constexpr float TW_PHENOL_QUINONE_METHIDE = 0.92f; // phenol -> para-quinone methide
    static constexpr float TW_TRIAZOLE_NH_SHIFT    = 0.50f; // 1,2,3-triazole NH shift (symmetric)
    static constexpr float TW_BENZIMIDAZOLE_NH     = 0.50f; // NH-1 <-> NH-3 shift (symmetric)
    static constexpr float TW_PYRIDONE_HYDROXYPYRIDINE = 0.85f; // 2-pyridinol <-> 2-pyridone
    static constexpr float TW_BARBITURIC           = 0.80f; // tri-keto <-> enol forms
    static constexpr float TW_ALLYL_SHIFT          = 0.70f; // X-C=C-C <-> C=C-C-X (X=OH,SH,NH)
    static constexpr float TW_SELENOL_SELENOKETONE = 0.94f; // C=Se <-> C(-SeH)

    // --- Bit-parallel adjacency (n x words) ---
    std::vector<std::vector<uint64_t>> adjLong;

    // --- Bond storage: dense for small molecules, sparse for large ---
    std::vector<std::vector<int>>  bondOrdMatrix;   // dense n x n
    std::vector<std::vector<bool>> bondRingMatrix;
    std::vector<std::vector<bool>> bondAromMatrix;
    bool useDense = false;

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

    // --- NLF caches (lazy, mutable for const-correct access) ---
    mutable std::vector<std::vector<int>> cachedNLF1, cachedNLF2, cachedNLF3;

    // --- Performance caches (v6.5.3 parity with Java) ---
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
        if (useDense) return bondOrdMatrix[i][j];
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() ? it->second[0] : 0;
    }

    bool bondInRing(int i, int j) const {
        if (useDense) return bondRingMatrix[i][j];
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() && it->second[1] != 0;
    }

    bool bondAromatic(int i, int j) const {
        if (useDense) return bondAromMatrix[i][j];
        auto it = sparseBondProps.find(bondKey(i, j));
        return it != sparseBondProps.end() && it->second[2] != 0;
    }

    bool hasBond(int i, int j) const { return bondOrder(i, j) != 0; }

    void invalidateBondDerivedState() {
        canonicalComputed_ = false;
        canonicalHash = 0;
        morganRank.clear();
        canonicalLabel.clear();
        orbit.clear();
        autGenerators_.clear();
        autGeneratorsTruncated_ = false;
        patternFPComputed_ = false;
        for (int w = 0; w < FP_WORDS; ++w) patternFP_[w] = 0;
        cachedPharmacophoreFeatures.clear();
        pharmComputed_ = false;
        tautomerClass.clear();
        tautomerWeight.clear();
        if (static_cast<int>(formalCharge.size()) == n
            && static_cast<int>(ring.size()) == n
            && static_cast<int>(aromatic.size()) == n) {
            computeTautomerClasses();
            applySolventCorrection(solvent_);
            adjustWeightsForPH(pH_);
        } else {
            tautomerClass.assign(n, -1);
            tautomerWeight.assign(n, 1.0f);
        }
    }

    void setBondOrder(int i, int j, int ord) {
        if (useDense) { bondOrdMatrix[i][j] = ord; bondOrdMatrix[j][i] = ord; }
        else {
            auto key = bondKey(i, j);
            auto it = sparseBondProps.find(key);
            if (it != sparseBondProps.end()) it->second[0] = ord;
        }
        invalidateBondDerivedState();
    }

    void setBondAromatic(int i, int j, bool arom) {
        if (useDense) { bondAromMatrix[i][j] = arom; bondAromMatrix[j][i] = arom; }
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

    const std::vector<int>& getOrbit() const          { return orbit; }
    const std::vector<int>& getCanonicalLabeling() const { return canonicalLabel; }
    uint64_t getCanonicalHash() const                   { return canonicalHash; }

    /// Return automorphism generators discovered during canonical labeling.
    /// Each generator is a permutation vector of length n where gen[i] is the
    /// image of atom i.  The full automorphism group is the closure of these
    /// generators.
    const std::vector<std::vector<int>>& getAutomorphismGenerators() const {
        ensureCanonical();
        return autGenerators_;
    }

    /// True when the generator list was capped during canonical search.
    /// The orbit partition is still exact; only the explicit generator list
    /// may be incomplete.
    bool automorphismGeneratorsTruncated() const {
        ensureCanonical();
        return autGeneratorsTruncated_;
    }

    // ========================================================================
    // Ring-count per atom
    // ========================================================================

    /** Compute ring counts from a given ring set (SSSR or RCB). */
    void computeRingCounts(const std::vector<std::vector<int>>& rings) const {
        std::fill(ringCount.begin(), ringCount.end(), 0);
        for (const auto& ring : rings) {
            for (int atom : ring) ringCount[atom]++;
        }
    }

    /** Default: compute ring counts from internal SSSR.
     *  For deterministic results on symmetric molecules, call
     *  computeRingCounts(smsd::computeRelevantCycles(*this)) after including ring_finder.hpp. */
    void computeRingCounts() const {
        computeRingCounts(computeRings());
    }

    static bool isAromaticDonorElement(int z) {
        return z == 7 || z == 8 || z == 15 || z == 16 || z == 33 || z == 34
            || z == 51 || z == 52;
    }

    std::pair<bool, bool> exocyclicMultipleBondInfo(
        int atom,
        int prev = -1,
        int next = -1) const {
        bool hasExocyclicMultipleBond = false;
        bool hasExocyclicMultipleBondToHetero = false;
        for (int nb : neighbors[atom]) {
            if (nb == prev || nb == next) continue;
            int ord = bondOrder(atom, nb);
            if (ord < 2 || ord == 4) continue;
            hasExocyclicMultipleBond = true;
            if (atomicNum[nb] != 6) hasExocyclicMultipleBondToHetero = true;
        }
        return {hasExocyclicMultipleBond, hasExocyclicMultipleBondToHetero};
    }

    static bool isSimpleCycle(const std::vector<int>& cycle) {
        if (cycle.size() < 3) return false;
        std::unordered_set<int> seen;
        seen.reserve(cycle.size());
        for (int atom : cycle) {
            if (!seen.insert(atom).second) return false;
        }
        return true;
    }

    int aromaticPiContribution(const std::vector<int>& cycle, int pos) const {
        const int len = static_cast<int>(cycle.size());
        const int atom = cycle[pos];
        const int prev = cycle[(pos + len - 1) % len];
        const int next = cycle[(pos + 1) % len];

        int cyclePiBonds = 0;
        bool hasExplicitAromaticBond = false;
        for (int nb : {prev, next}) {
            int ord = bondOrder(atom, nb);
            if (ord == 3) return -1;
            if (ord == 4) {
                hasExplicitAromaticBond = true;
            } else if (ord == 2) {
                ++cyclePiBonds;
            }
        }

        if (hasExplicitAromaticBond) return 1;
        if (cyclePiBonds == 1) return 1;
        if (cyclePiBonds > 1) return -1;

        auto [hasExocyclicMultipleBond, hasExocyclicMultipleBondToHetero]
            = exocyclicMultipleBondInfo(atom, prev, next);

        const int z = atomicNum[atom];
        const int q = formalCharge.empty() ? 0 : formalCharge[atom];

        if (z == 6) {
            if (q < 0 && !hasExocyclicMultipleBond) return 2;
            if (q == 0 && hasExocyclicMultipleBondToHetero) return 0;
            if (q > 0) return 0;
            return -1;
        }

        if (z == 5 || z == 13) {
            if (q < 0 && !hasExocyclicMultipleBond) return 2;
            return -1;
        }

        if (isAromaticDonorElement(z)) {
            if (q > 0) return 0;
            if (hasExocyclicMultipleBond) return -1;
            return 2;
        }

        return -1;
    }

    bool cycleLooksAromatic(const std::vector<int>& cycle) const {
        const int len = static_cast<int>(cycle.size());
        if (len < 3 || len > 24) return false;

        int piElectrons = 0;
        bool sawCyclePi = false;
        for (int pos = 0; pos < len; ++pos) {
            int contrib = aromaticPiContribution(cycle, pos);
            if (contrib < 0) return false;
            if (contrib == 1 || contrib == 2) sawCyclePi = true;
            piElectrons += contrib;
        }

        if (!sawCyclePi || piElectrons < 2) return false;
        return ((piElectrons - 2) % 4) == 0;
    }

    void enumerateRingCyclesDfs(
        int root,
        int current,
        const std::vector<std::vector<int>>& ringAdj,
        const std::vector<uint8_t>& inComponent,
        std::vector<uint8_t>& inPath,
        std::vector<int>& path,
        std::set<std::vector<int64_t>>& seenCycles,
        std::vector<std::vector<int>>& cycles,
        size_t maxCycles) const {
        if (cycles.size() >= maxCycles || path.size() > 24) return;

        for (int nb : ringAdj[current]) {
            if (!inComponent[nb] || nb < root) continue;
            if (nb == root) {
                if (path.size() >= 3) {
                    auto key = canonicalCycleKey(path);
                    if (seenCycles.insert(key).second) cycles.push_back(path);
                }
                continue;
            }
            if (inPath[nb] || path.size() >= 24) continue;

            inPath[nb] = 1;
            path.push_back(nb);
            enumerateRingCyclesDfs(
                root, nb, ringAdj, inComponent, inPath, path,
                seenCycles, cycles, maxCycles);
            path.pop_back();
            inPath[nb] = 0;

            if (cycles.size() >= maxCycles) return;
        }
    }

    std::vector<std::vector<int>> collectAromaticCandidateCycles() const {
        std::vector<std::vector<int>> cycles;
        std::set<std::vector<int64_t>> seenCycles;

        auto addCycle = [&](const std::vector<int>& cycle) {
            if (!isSimpleCycle(cycle)) return;
            if (cycle.size() < 3 || cycle.size() > 24) return;
            auto key = canonicalCycleKey(cycle);
            if (seenCycles.insert(key).second) cycles.push_back(cycle);
        };

        for (const auto& cycle : computeRelevantCycles()) addCycle(cycle);

        std::vector<std::vector<int>> ringAdj(n);
        for (int i = 0; i < n; ++i) {
            if (!ring[i]) continue;
            for (int nb : neighbors[i]) {
                if (i < nb && bondInRing(i, nb)) {
                    ringAdj[i].push_back(nb);
                    ringAdj[nb].push_back(i);
                }
            }
        }

        std::vector<uint8_t> visited(n, uint8_t(0));
        std::vector<uint8_t> inComponent(n, uint8_t(0));
        std::vector<uint8_t> inPath(n, uint8_t(0));
        for (int start = 0; start < n; ++start) {
            if (!ring[start] || visited[start] || ringAdj[start].empty()) continue;

            std::vector<int> component;
            std::deque<int> queue;
            queue.push_back(start);
            visited[start] = 1;
            while (!queue.empty()) {
                int atom = queue.front();
                queue.pop_front();
                component.push_back(atom);
                for (int nb : ringAdj[atom]) {
                    if (visited[nb]) continue;
                    visited[nb] = 1;
                    queue.push_back(nb);
                }
            }

            int edgeCount = 0;
            for (int atom : component) edgeCount += static_cast<int>(ringAdj[atom].size());
            edgeCount /= 2;
            const int cycleRank = edgeCount - static_cast<int>(component.size()) + 1;
            if (cycleRank <= 1 || component.size() > 24 || cycleRank > 8) continue;

            std::fill(inComponent.begin(), inComponent.end(), uint8_t(0));
            for (int atom : component) inComponent[atom] = 1;

            std::sort(component.begin(), component.end());
            for (int root : component) {
                if (ringAdj[root].size() < 2) continue;
                std::fill(inPath.begin(), inPath.end(), uint8_t(0));
                std::vector<int> path = { root };
                inPath[root] = 1;

                for (int nb : ringAdj[root]) {
                    if (!inComponent[nb] || nb <= root) continue;
                    inPath[nb] = 1;
                    path.push_back(nb);
                    enumerateRingCyclesDfs(
                        root, nb, ringAdj, inComponent, inPath, path,
                        seenCycles, cycles, 4096);
                    path.pop_back();
                    inPath[nb] = 0;
                    if (cycles.size() >= 4096) break;
                }
                inPath[root] = 0;
                if (cycles.size() >= 4096) break;
            }
            if (cycles.size() >= 4096) break;
        }

        return cycles;
    }

    void clearRingFlagsFromTopology() {
        std::fill(ring.begin(), ring.end(), uint8_t(0));
        if (useDense) {
            for (auto& row : bondRingMatrix) {
                std::fill(row.begin(), row.end(), false);
            }
        } else {
            for (auto& kv : sparseBondProps) {
                kv.second[1] = 0;
            }
        }
    }

    void markRingEdge(int i, int j) {
        ring[i] = 1;
        ring[j] = 1;
        if (useDense) {
            bondRingMatrix[i][j] = true;
            bondRingMatrix[j][i] = true;
        } else {
            auto it = sparseBondProps.find(bondKey(i, j));
            if (it != sparseBondProps.end()) it->second[1] = 1;
        }
    }

    void markAromaticEdge(int i, int j) {
        if (useDense) {
            bondOrdMatrix[i][j] = 4;
            bondOrdMatrix[j][i] = 4;
            bondAromMatrix[i][j] = true;
            bondAromMatrix[j][i] = true;
        } else {
            auto it = sparseBondProps.find(bondKey(i, j));
            if (it != sparseBondProps.end()) {
                it->second[0] = 4;
                it->second[2] = 1;
            }
        }
    }

    void markKekuleEdge(int i, int j, int order) {
        if (useDense) {
            bondOrdMatrix[i][j] = order;
            bondOrdMatrix[j][i] = order;
            bondAromMatrix[i][j] = false;
            bondAromMatrix[j][i] = false;
        } else {
            auto it = sparseBondProps.find(bondKey(i, j));
            if (it != sparseBondProps.end()) {
                it->second[0] = order;
                it->second[2] = 0;
            }
        }
    }

    void refreshAtomLabels() {
        label.resize(n);
        for (int i = 0; i < n; ++i) {
            label[i] = (atomicNum[i] << 2)
                     | (aromatic[i] ? 2 : 0)
                     | (ring[i] ? 1 : 0);
        }
    }

    void normalizeRingAndAromaticity(
        AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE) {
        if (n == 0) return;
        (void) model;

        clearRingFlagsFromTopology();
        const auto& sssr = computeRings();
        for (const auto& cycle : sssr) {
            for (size_t i = 0; i < cycle.size(); ++i) {
                int a = cycle[i];
                int b = cycle[(i + 1) % cycle.size()];
                markRingEdge(a, b);
            }
        }

        std::vector<uint8_t> perceivedAromatic = aromatic;
        std::vector<std::pair<int, int>> aromaticEdges;
        for (int i = 0; i < n; ++i) {
            if (!ring[i]) continue;
            for (int nb : neighbors[i]) {
                if (i < nb && bondInRing(i, nb) && bondAromatic(i, nb)) {
                    aromaticEdges.emplace_back(i, nb);
                }
            }
        }

        auto candidateCycles = collectAromaticCandidateCycles();
        for (const auto& cycle : candidateCycles) {
            if (!cycleLooksAromatic(cycle)) continue;
            for (int atom : cycle) perceivedAromatic[atom] = 1;
            for (size_t i = 0; i < cycle.size(); ++i) {
                int a = cycle[i];
                int b = cycle[(i + 1) % cycle.size()];
                aromaticEdges.emplace_back(std::min(a, b), std::max(a, b));
            }
        }

        for (int i = 0; i < n; ++i) {
            if (perceivedAromatic[i]) aromatic[i] = 1;
        }

        std::sort(aromaticEdges.begin(), aromaticEdges.end());
        aromaticEdges.erase(std::unique(aromaticEdges.begin(), aromaticEdges.end()), aromaticEdges.end());
        for (const auto& edge : aromaticEdges) {
            markAromaticEdge(edge.first, edge.second);
        }
    }

    void perceiveAromaticity(
        AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE) {
        normalizeRingAndAromaticity(model);
        refreshAtomLabels();
        ringCount.assign(n, 0);
        ringCountsComputed_ = false;
        invalidateBondDerivedState();
    }

    int aromaticKekuleDemand(int atom) const {
        if (!aromatic[atom]) return 0;
        const int z = atomicNum[atom];
        const int q = formalCharge.empty() ? 0 : formalCharge[atom];
        const int h = hydrogenCount.empty() ? 0 : hydrogenCount[atom];
        auto [hasExocyclicMultipleBond, hasExocyclicMultipleBondToHetero]
            = exocyclicMultipleBondInfo(atom);
        (void) hasExocyclicMultipleBond;

        if (z == 6) {
            if (q == 0 && hasExocyclicMultipleBondToHetero) return 0;
            return q == 0 ? 1 : 0;
        }
        if (z == 5 || z == 13) return q >= 0 ? 1 : 0;
        if (z == 7 || z == 15) {
            if (h > 0 && q <= 0) return 0;
            if (q < 0) return 0;
            return 1;
        }
        if (z == 8 || z == 16 || z == 34 || z == 52) return 0;
        return 1;
    }

    bool kekulizeAromaticComponent(const std::vector<int>& component) {
        if (component.empty()) return true;

        std::unordered_set<int> inComponent(component.begin(), component.end());
        std::vector<std::pair<int, int>> aromaticEdges;
        aromaticEdges.reserve(component.size() * 2);
        for (int atom : component) {
            for (int nb : neighbors[atom]) {
                if (atom < nb && inComponent.count(nb)
                    && (bondAromatic(atom, nb) || bondOrder(atom, nb) == 4)) {
                    aromaticEdges.emplace_back(atom, nb);
                }
            }
        }
        if (aromaticEdges.empty()) {
            for (int atom : component) aromatic[atom] = 0;
            return true;
        }

        std::vector<int> demandAtoms;
        demandAtoms.reserve(component.size());
        std::unordered_map<int, int> localIndex;
        for (int atom : component) {
            if (aromaticKekuleDemand(atom) > 0) {
                localIndex[atom] = static_cast<int>(demandAtoms.size());
                demandAtoms.push_back(atom);
            }
        }

        if (!demandAtoms.empty()) {
            std::vector<std::vector<int>> demandAdj(demandAtoms.size());
            for (const auto& edge : aromaticEdges) {
                auto ita = localIndex.find(edge.first);
                auto itb = localIndex.find(edge.second);
                if (ita == localIndex.end() || itb == localIndex.end()) continue;
                demandAdj[ita->second].push_back(itb->second);
                demandAdj[itb->second].push_back(ita->second);
            }

            std::vector<int> color(demandAtoms.size(), -1);
            std::deque<int> queue;
            for (int start = 0; start < static_cast<int>(demandAtoms.size()); ++start) {
                if (color[start] != -1) continue;
                color[start] = 0;
                queue.push_back(start);
                while (!queue.empty()) {
                    int u = queue.front();
                    queue.pop_front();
                    for (int v : demandAdj[u]) {
                        if (color[v] == -1) {
                            color[v] = color[u] ^ 1;
                            queue.push_back(v);
                        } else if (color[v] == color[u]) {
                            return false;
                        }
                    }
                }
            }

            std::vector<int> leftNodes, rightNodes;
            leftNodes.reserve(demandAtoms.size());
            rightNodes.reserve(demandAtoms.size());
            for (int i = 0; i < static_cast<int>(demandAtoms.size()); ++i) {
                if (color[i] == 0) leftNodes.push_back(i);
                else rightNodes.push_back(i);
            }
            if (leftNodes.size() != rightNodes.size()) return false;

            std::unordered_map<int, int> rightPos;
            rightPos.reserve(rightNodes.size());
            for (int i = 0; i < static_cast<int>(rightNodes.size()); ++i) {
                rightPos[rightNodes[i]] = i;
            }

            std::vector<std::vector<int>> matchAdj(leftNodes.size());
            for (int li = 0; li < static_cast<int>(leftNodes.size()); ++li) {
                int u = leftNodes[li];
                for (int v : demandAdj[u]) {
                    auto it = rightPos.find(v);
                    if (it != rightPos.end()) matchAdj[li].push_back(it->second);
                }
                if (matchAdj[li].empty()) return false;
            }

            std::vector<int> matchRight(rightNodes.size(), -1);
            std::function<bool(int, std::vector<uint8_t>&)> augment =
                [&](int li, std::vector<uint8_t>& seen) {
                    for (int ri : matchAdj[li]) {
                        if (seen[ri]) continue;
                        seen[ri] = 1;
                        if (matchRight[ri] == -1 || augment(matchRight[ri], seen)) {
                            matchRight[ri] = li;
                            return true;
                        }
                    }
                    return false;
                };

            for (int li = 0; li < static_cast<int>(leftNodes.size()); ++li) {
                std::vector<uint8_t> seen(rightNodes.size(), uint8_t(0));
                if (!augment(li, seen)) return false;
            }

            for (const auto& edge : aromaticEdges) {
                markKekuleEdge(edge.first, edge.second, 1);
            }
            for (int atom : component) aromatic[atom] = 0;

            for (int ri = 0; ri < static_cast<int>(matchRight.size()); ++ri) {
                int li = matchRight[ri];
                if (li < 0) return false;
                int a = demandAtoms[leftNodes[li]];
                int b = demandAtoms[rightNodes[ri]];
                markKekuleEdge(a, b, 2);
            }
            return true;
        }

        for (const auto& edge : aromaticEdges) {
            markKekuleEdge(edge.first, edge.second, 1);
        }
        for (int atom : component) aromatic[atom] = 0;
        return true;
    }

    bool kekulize() {
        std::vector<uint8_t> visited(n, uint8_t(0));
        for (int start = 0; start < n; ++start) {
            if (!aromatic[start] || visited[start]) continue;

            std::vector<int> component;
            std::deque<int> queue;
            queue.push_back(start);
            visited[start] = 1;
            while (!queue.empty()) {
                int atom = queue.front();
                queue.pop_front();
                component.push_back(atom);
                for (int nb : neighbors[atom]) {
                    if (visited[nb] || !aromatic[nb]) continue;
                    if (!(bondAromatic(atom, nb) || bondOrder(atom, nb) == 4)) continue;
                    visited[nb] = 1;
                    queue.push_back(nb);
                }
            }

            if (!kekulizeAromaticComponent(component)) return false;
        }

        refreshAtomLabels();
        invalidateBondDerivedState();
        return true;
    }

    bool dearomatize() {
        return kekulize();
    }

    // ========================================================================
    // SSSR ring computation (Horton-style shortest-path basis)
    // ========================================================================

    const std::vector<std::vector<int>>& computeRings() const {
        if (sssrComputed) return sssrRings;
        sssrComputed = true;
        if (n == 0) return sssrRings;

        // Edge count
        int edgeCount = 0;
        for (int i = 0; i < n; i++)
            edgeCount += static_cast<int>(neighbors[i].size());
        edgeCount /= 2;

        // Connected components via BFS
        std::vector<bool> visited(n, false);
        int components = 0;
        for (int i = 0; i < n; i++) {
            if (visited[i]) continue;
            components++;
            std::deque<int> bfsQ;
            bfsQ.push_back(i);
            visited[i] = true;
            while (!bfsQ.empty()) {
                int u = bfsQ.front(); bfsQ.pop_front();
                for (int v : neighbors[u])
                    if (!visited[v]) { visited[v] = true; bfsQ.push_back(v); }
            }
        }

        int cycleRank = edgeCount - n + components;
        if (cycleRank <= 0) return sssrRings;

        // For each edge (start, nb) with nb > start, find shortest path
        // from start to nb that avoids that direct edge.
        std::vector<std::vector<int>> candidates;
        for (int start = 0; start < n; start++) {
            for (int nb : neighbors[start]) {
                if (nb <= start) continue;
                std::vector<int> parent(n, -2);
                parent[start] = -1;
                std::deque<int> queue;
                queue.push_back(start);
                bool found = false;
                while (!queue.empty() && !found) {
                    int u = queue.front(); queue.pop_front();
                    for (int v : neighbors[u]) {
                        if (u == start && v == nb) continue;
                        if (v == start && u == nb) continue;
                        if (parent[v] != -2) continue;
                        parent[v] = u;
                        if (v == nb) { found = true; break; }
                        queue.push_back(v);
                    }
                }
                if (found) {
                    std::vector<int> path;
                    int cur = nb;
                    while (cur != start && cur != -1) { path.push_back(cur); cur = parent[cur]; }
                    if (cur == start) {
                        path.push_back(start);
                        if (path.size() >= 3 && path.size() <= 20)
                            candidates.push_back(std::move(path));
                    }
                }
            }
        }

        std::sort(candidates.begin(), candidates.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
                return a.size() < b.size();
            });

        // Greedy basis: pick cycles that contribute at least one new edge
        std::vector<std::vector<int>> basis;
        std::vector<std::unordered_set<int64_t>> basisEdges;

        for (auto& cycle : candidates) {
            if (static_cast<int>(basis.size()) >= cycleRank) break;
            std::unordered_set<int64_t> edges;
            for (size_t i = 0; i < cycle.size(); i++) {
                int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
                edges.insert(bondKey(std::min(a, b), std::max(a, b)));
            }
            std::unordered_set<int64_t> allBasis;
            for (auto& be : basisEdges) allBasis.insert(be.begin(), be.end());
            bool hasNewEdge = false;
            for (int64_t e : edges) {
                if (allBasis.find(e) == allBasis.end()) { hasNewEdge = true; break; }
            }
            if (hasNewEdge || basis.empty()) {
                basis.push_back(cycle);
                basisEdges.push_back(std::move(edges));
            }
        }

        sssrRings = std::move(basis);
        return sssrRings;
    }

    // ========================================================================
    // Relevant cycles (union of all minimum cycle bases)
    // ========================================================================

private:

    /** Find the lowest set bit position in a GF(2) vector. */
    static int findPivotBit(const std::vector<uint64_t>& vec) {
        for (size_t w = 0; w < vec.size(); w++) {
            if (vec[w] != 0) return static_cast<int>((w << 6) + smsd::ctz64(vec[w]));
        }
        return -1;
    }

    /** Reduce a GF(2) vector against a row echelon basis. Returns true if zero. */
    static bool reduceVectorGF2(std::vector<uint64_t>& vec,
                                const std::vector<std::vector<uint64_t>>& basis,
                                int basisSize) {
        for (int i = 0; i < basisSize; i++) {
            int pivot = findPivotBit(basis[i]);
            if (pivot >= 0 && ((vec[pivot >> 6] >> (pivot & 63)) & 1ULL) != 0) {
                for (size_t w = 0; w < vec.size(); w++) vec[w] ^= basis[i][w];
            }
        }
        for (auto w : vec) if (w != 0) return false;
        return true;
    }

    /** Add a GF(2) vector to row echelon basis if linearly independent. */
    static int addToBasisGF2(std::vector<uint64_t>& vec,
                             std::vector<std::vector<uint64_t>>& basis,
                             int basisSize) {
        for (int i = 0; i < basisSize; i++) {
            int pivot = findPivotBit(basis[i]);
            if (pivot >= 0 && ((vec[pivot >> 6] >> (pivot & 63)) & 1ULL) != 0) {
                for (size_t w = 0; w < vec.size(); w++) vec[w] ^= basis[i][w];
            }
        }
        bool nonZero = false;
        for (auto w : vec) if (w != 0) { nonZero = true; break; }
        if (nonZero) {
            if (basisSize >= static_cast<int>(basis.size())) basis.push_back(vec);
            else basis[basisSize] = vec;
            return basisSize + 1;
        }
        return basisSize;
    }

    /** Canonical edge-set key for cycle deduplication. */
    static std::vector<int64_t> canonicalCycleKey(const std::vector<int>& cycle) {
        std::vector<int64_t> keys;
        keys.reserve(cycle.size());
        for (size_t i = 0; i < cycle.size(); i++) {
            int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
            int lo = std::min(a, b), hi = std::max(a, b);
            keys.push_back((static_cast<int64_t>(lo) << 32) | static_cast<int64_t>(hi));
        }
        std::sort(keys.begin(), keys.end());
        return keys;
    }

    /** Enumerate all shortest paths from root to target in a BFS tree. */
    static void buildPathsHelper(int root, int node, const std::vector<int>& dist,
                                 const std::vector<std::vector<int>>& parents,
                                 std::vector<int>& current,
                                 std::vector<std::vector<int>>& result, int limit) {
        if (static_cast<int>(result.size()) >= limit) return;
        if (node == root) {
            std::vector<int> path(current.rbegin(), current.rend());
            result.push_back(std::move(path));
            return;
        }
        for (int p : parents[node]) {
            if (dist[p] != dist[node] - 1) continue;
            current.push_back(p);
            buildPathsHelper(root, p, dist, parents, current, result, limit);
            current.pop_back();
        }
    }

    static std::vector<std::vector<int>> allShortestPathsCpp(int root, int target,
            const std::vector<int>& dist,
            const std::vector<std::vector<int>>& parents) {
        std::vector<std::vector<int>> result;
        std::vector<int> current = { target };
        buildPathsHelper(root, target, dist, parents, current, result, 50);
        return result;
    }

    // Union-Find helpers
    static int ufFind(std::vector<int>& parent, int i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    }
    static void ufUnion(std::vector<int>& parent, int a, int b) {
        a = ufFind(parent, a); b = ufFind(parent, b);
        if (a != b) parent[a] = b;
    }

public:

    /**
     * Computes the set of relevant cycles (union of all minimum cycle bases).
     * A cycle is relevant if it cannot be expressed as the XOR (symmetric
     * difference) of strictly shorter cycles.
     */
    std::vector<std::vector<int>> computeRelevantCycles() const {
        auto& mcb = computeRings();
        if (mcb.empty()) return {};

        // Build edge index
        std::unordered_map<int64_t, int> edgeIndex;
        for (int i = 0; i < n; i++) {
            for (int j : neighbors[i]) {
                if (j > i) {
                    int64_t key = bondKey(std::min(i, j), std::max(i, j));
                    if (edgeIndex.find(key) == edgeIndex.end()) {
                        int idx = static_cast<int>(edgeIndex.size());
                        edgeIndex[key] = idx;
                    }
                }
            }
        }
        int numEdges = static_cast<int>(edgeIndex.size());
        int longWords = (numEdges + 63) >> 6;

        // Find max cycle length in MCB
        int maxLen = 0;
        for (auto& c : mcb) if (static_cast<int>(c.size()) > maxLen) maxLen = static_cast<int>(c.size());

        // Collect all short cycles
        std::set<std::vector<int64_t>> uniqueCycleSet;
        std::vector<std::vector<int>> allCycles;

        // Add MCB cycles
        for (auto& c : mcb) {
            auto key = canonicalCycleKey(c);
            if (uniqueCycleSet.insert(key).second) allCycles.push_back(c);
        }

        // Enumerate short cycles via BFS from each vertex
        for (int root = 0; root < n; root++) {
            std::vector<int> dist(n, -1);
            dist[root] = 0;
            std::vector<std::vector<int>> parents(n);
            std::deque<int> queue;
            queue.push_back(root);
            while (!queue.empty()) {
                int u = queue.front(); queue.pop_front();
                if (dist[u] >= maxLen / 2) continue;
                for (int v : neighbors[u]) {
                    if (dist[v] == -1) {
                        dist[v] = dist[u] + 1;
                        parents[v].push_back(u);
                        queue.push_back(v);
                    } else if (dist[v] == dist[u] + 1) {
                        parents[v].push_back(u);
                    }
                }
            }

            // Odd-length cycles from edges
            for (int u = 0; u < n; u++) {
                if (dist[u] == -1) continue;
                for (int v : neighbors[u]) {
                    if (v <= u || dist[v] == -1) continue;
                    int cycleLen = dist[u] + dist[v] + 1;
                    if (cycleLen < 3 || cycleLen > maxLen) continue;

                    auto pathsU = allShortestPathsCpp(root, u, dist, parents);
                    auto pathsV = allShortestPathsCpp(root, v, dist, parents);

                    for (auto& pu : pathsU) {
                        for (auto& pv : pathsV) {
                            std::unordered_set<int> puSet(pu.begin(), pu.end());
                            bool overlap = false;
                            for (size_t k = 1; k < pv.size(); k++) {
                                if (puSet.count(pv[k])) { overlap = true; break; }
                            }
                            if (overlap) continue;

                            std::vector<int> cycle;
                            cycle.reserve(pu.size() + pv.size() - 1);
                            for (int x : pu) cycle.push_back(x);
                            for (int k = static_cast<int>(pv.size()) - 2; k >= 0; k--)
                                cycle.push_back(pv[k]);

                            if (static_cast<int>(cycle.size()) >= 3 &&
                                static_cast<int>(cycle.size()) <= maxLen) {
                                auto key = canonicalCycleKey(cycle);
                                if (uniqueCycleSet.insert(key).second)
                                    allCycles.push_back(std::move(cycle));
                            }
                        }
                    }
                }
            }

            // Even-length cycles from two disjoint paths to same vertex
            for (int u = 0; u < n; u++) {
                if (dist[u] == -1 || static_cast<int>(parents[u].size()) < 2) continue;
                int cycleLen = 2 * dist[u];
                if (cycleLen < 3 || cycleLen > maxLen) continue;

                auto paths = allShortestPathsCpp(root, u, dist, parents);
                for (size_t pi = 0; pi < paths.size(); pi++) {
                    for (size_t pj = pi + 1; pj < paths.size(); pj++) {
                        auto& p1 = paths[pi]; auto& p2 = paths[pj];
                        std::unordered_set<int> p1Int;
                        for (size_t k = 1; k + 1 < p1.size(); k++) p1Int.insert(p1[k]);
                        bool disjoint = true;
                        for (size_t k = 1; k + 1 < p2.size(); k++) {
                            if (p1Int.count(p2[k])) { disjoint = false; break; }
                        }
                        if (!disjoint) continue;

                        std::vector<int> cycle;
                        cycle.reserve(cycleLen);
                        for (size_t k = 0; k < p1.size(); k++) cycle.push_back(p1[k]);
                        for (int k = static_cast<int>(p2.size()) - 2; k >= 1; k--)
                            cycle.push_back(p2[k]);

                        if (static_cast<int>(cycle.size()) >= 3 &&
                            static_cast<int>(cycle.size()) <= maxLen) {
                            auto key = canonicalCycleKey(cycle);
                            if (uniqueCycleSet.insert(key).second)
                                allCycles.push_back(std::move(cycle));
                        }
                    }
                }
            }
        }

        // Sort by length
        std::sort(allCycles.begin(), allCycles.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
                return a.size() < b.size();
            });

        // Build GF(2) edge-incidence vectors
        std::vector<std::vector<uint64_t>> vectors(allCycles.size(),
            std::vector<uint64_t>(longWords, 0));
        for (size_t ci = 0; ci < allCycles.size(); ci++) {
            auto& c = allCycles[ci];
            for (size_t i = 0; i < c.size(); i++) {
                int a = c[i], b = c[(i + 1) % c.size()];
                int64_t key = bondKey(std::min(a, b), std::max(a, b));
                int idx = edgeIndex[key];
                vectors[ci][idx >> 6] ^= 1ULL << (idx & 63);
            }
        }

        // Filter relevant cycles via Gaussian elimination on GF(2)
        std::vector<std::vector<int>> relevant;
        std::vector<std::vector<uint64_t>> reducedBasis(numEdges,
            std::vector<uint64_t>(longWords, 0));
        int basisSize = 0;

        size_t ci = 0;
        while (ci < allCycles.size()) {
            int len = static_cast<int>(allCycles[ci].size());
            size_t cj = ci;
            while (cj < allCycles.size() && static_cast<int>(allCycles[cj].size()) == len) cj++;

            // Test each cycle against basis of strictly shorter cycles
            for (size_t k = ci; k < cj; k++) {
                auto vec = vectors[k];
                bool inSpan = reduceVectorGF2(vec, reducedBasis, basisSize);
                if (!inSpan) relevant.push_back(allCycles[k]);
            }

            // Add all cycles of this length to basis
            for (size_t k = ci; k < cj; k++) {
                auto vec = vectors[k];
                basisSize = addToBasisGF2(vec, reducedBasis, basisSize);
            }

            ci = cj;
        }

        return relevant;
    }

    /**
     * Computes Unique Ring Families (URFs). Two relevant cycles belong to the
     * same URF if one can be obtained from the other by a single edge exchange
     * (their symmetric difference has exactly two edges sharing a vertex).
     */
    std::vector<std::vector<std::vector<int>>> computeURFs() {
        auto rc = computeRelevantCycles();
        if (rc.empty()) return {};

        // Build edge index
        std::unordered_map<int64_t, int> edgeIndex;
        for (int i = 0; i < n; i++) {
            for (int j : neighbors[i]) {
                if (j > i) {
                    int64_t key = bondKey(std::min(i, j), std::max(i, j));
                    if (edgeIndex.find(key) == edgeIndex.end()) {
                        int idx = static_cast<int>(edgeIndex.size());
                        edgeIndex[key] = idx;
                    }
                }
            }
        }
        int numEdges = static_cast<int>(edgeIndex.size());
        int longWords = (numEdges + 63) >> 6;

        // Build reverse lookup
        std::vector<std::pair<int,int>> edgeByIndex(numEdges);
        for (auto& [key, idx] : edgeIndex) {
            edgeByIndex[idx] = { static_cast<int>(key >> 32),
                                 static_cast<int>(key & 0xFFFFFFFF) };
        }

        // Build edge-incidence vectors
        int nrc = static_cast<int>(rc.size());
        std::vector<std::vector<uint64_t>> vectors(nrc, std::vector<uint64_t>(longWords, 0));
        for (int ci = 0; ci < nrc; ci++) {
            auto& c = rc[ci];
            for (size_t i = 0; i < c.size(); i++) {
                int a = c[i], b = c[(i + 1) % c.size()];
                int64_t key = bondKey(std::min(a, b), std::max(a, b));
                int idx = edgeIndex[key];
                vectors[ci][idx >> 6] ^= 1ULL << (idx & 63);
            }
        }

        // Union-Find grouping
        std::vector<int> family(nrc);
        std::iota(family.begin(), family.end(), 0);

        for (int i = 0; i < nrc; i++) {
            for (int j = i + 1; j < nrc; j++) {
                if (rc[i].size() != rc[j].size()) continue;
                int popcount = 0;
                std::vector<uint64_t> xorVec(longWords);
                for (int w = 0; w < longWords; w++) {
                    xorVec[w] = vectors[i][w] ^ vectors[j][w];
                    popcount += smsd::popcount64(xorVec[w]);
                }
                if (popcount == 2) {
                    int diffEdges[2]; int di = 0;
                    for (int w = 0; w < longWords && di < 2; w++) {
                        uint64_t bits = xorVec[w];
                        while (bits != 0 && di < 2) {
                            int bit = smsd::ctz64(bits);
                            diffEdges[di++] = (w << 6) + bit;
                            bits &= bits - 1;
                        }
                    }
                    auto [u1, v1] = edgeByIndex[diffEdges[0]];
                    auto [u2, v2] = edgeByIndex[diffEdges[1]];
                    if (u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2) {
                        ufUnion(family, i, j);
                    }
                }
            }
        }

        // Group by family
        std::unordered_map<int, std::vector<int>> groups;
        for (int i = 0; i < nrc; i++) groups[ufFind(family, i)].push_back(i);

        std::vector<std::vector<std::vector<int>>> result;
        for (auto& [rep, members] : groups) {
            std::vector<std::vector<int>> fam;
            for (int m : members) fam.push_back(rc[m]);
            result.push_back(std::move(fam));
        }
        return result;
    }

    // ========================================================================
    // Tautomer class detection — 15 pKa-informed transforms (Tier 1)
    // ========================================================================

    void computeTautomerClasses() {
        tautomerClass.assign(n, -1);
        tautomerWeight.assign(n, 1.0f);
        int classId = 0;

        // Helper: assign cls + weight to a set of atoms.
        // If any atom already belongs to an existing class, merge by reusing
        // that class for all atoms in the group (prevents class fracture when
        // tautomeric motifs overlap on shared atoms).
        auto tc = [&](std::initializer_list<int> atoms, int cls, float w) {
            // Check if any atom already has a class — if so, merge into it
            for (int a : atoms)
                if (tautomerClass[a] != -1) { cls = tautomerClass[a]; break; }
            for (int a : atoms) {
                tautomerClass[a] = cls;
                if (tautomerWeight[a] == 1.0f) tautomerWeight[a] = w;
            }
        };

        for (int i = 0; i < n; i++) {
            if (tautomerClass[i] != -1) continue;

            // ------------------------------------------------------------------
            // T1: Keto/enol  C=O ↔ C-OH
            //     Amide carbonyl (N on C) → TW_AMIDE_IMIDIC; else TW_KETO_ENOL
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i]) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        bool hasN = false;
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 7) { hasN = true; break; }
                        }
                        float w = hasN ? TW_AMIDE_IMIDIC : TW_KETO_ENOL;
                        int cls = -1;
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 6 && !aromatic[k]) {
                                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k] == 1.0f) tautomerWeight[k] = w; }
                            }
                        }
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 7) {
                                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k] == 1.0f) tautomerWeight[k] = w; }
                            }
                        }
                        if (cls != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T2: Thione/thiol  C=S ↔ C-SH
            //     Thioamide N-C(=S) → TW_THIOAMIDE; else TW_THIONE_THIOL
            // ------------------------------------------------------------------
            if (atomicNum[i] == 16 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        bool hasN = false;
                        for (int k : neighbors[j]) { if (k != i && atomicNum[k] == 7) { hasN = true; break; } }
                        float w = hasN ? TW_THIOAMIDE : TW_THIONE_THIOL;
                        int cls = -1;
                        for (int k : neighbors[j]) {
                            if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) {
                                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k] == 1.0f) tautomerWeight[k] = w; }
                            }
                        }
                        if (cls != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T3: Nitroso/oxime  C=N-OH ↔ C(=O)-NH
            //     Requires N=C double bond to avoid matching simple amide N-C=O
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && !ring[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8) {
                                int cls = classId++;
                                tc({i, k}, cls, TW_NITROSO_OXIME);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T4: Phenol/quinone  arom-C-OH ↔ para-C=O
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
                        for (int a : neighbors[j]) {
                            if (a != i && atomicNum[a] == 6 && aromatic[a] && ring[a]) {
                                for (int b : neighbors[a]) {
                                    if (b != j && atomicNum[b] == 6 && aromatic[b] && ring[b]) {
                                        for (int c : neighbors[b]) {
                                            if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) {
                                                int cls = classId++;
                                                tc({i, j, c}, cls, TW_PHENOL_QUINONE);
                                                for (int d : neighbors[c]) {
                                                    if (atomicNum[d] == 8 && !aromatic[d] && tautomerClass[d] == -1) {
                                                        tautomerClass[d] = cls; tautomerWeight[d] = TW_PHENOL_QUINONE;
                                                    }
                                                }
                                                break;
                                            }
                                        }
                                        if (tautomerClass[i] != -1) break;
                                    }
                                }
                                if (tautomerClass[i] != -1) break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T5: 1,3-diketone  C(=O)-C-C(=O)  enolisation  (pKa ~8-11)
            //     Requires exactly one carbon bridge between the two carbonyls
            //     (true 1,3-pattern: O=C-CH-C=O, not fused/directly bonded carbonyls)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        for (int bridge : neighbors[j]) {
                            if (bridge != i && atomicNum[bridge] == 6 && !aromatic[bridge]) {
                                for (int m : neighbors[bridge]) {
                                    if (m != j && m != i && atomicNum[m] == 6) {
                                        // Verify true 1,3-separation: m must NOT be
                                        // a direct neighbor of j (which would mean
                                        // j and m are 1,2 not 1,3).
                                        bool directlyBonded = false;
                                        for (int nb : neighbors[j]) {
                                            if (nb == m) { directlyBonded = true; break; }
                                        }
                                        if (directlyBonded) continue;
                                        for (int o2 : neighbors[m]) {
                                            if (o2 != bridge && atomicNum[o2] == 8 && !aromatic[o2]
                                                && bondOrder(m, o2) == 2) {
                                                int cls = (tautomerClass[o2] != -1) ? tautomerClass[o2] : classId++;
                                                float w = TW_DIKETONE_ENOL;
                                                if (tautomerClass[i]      == -1) { tautomerClass[i]      = cls; tautomerWeight[i]      = w; }
                                                if (tautomerClass[bridge] == -1) { tautomerClass[bridge] = cls; tautomerWeight[bridge] = w; }
                                                if (tautomerClass[o2]     == -1) { tautomerClass[o2]     = cls; tautomerWeight[o2]     = w; }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T6: Pyridinone / lactam  aromatic N adjacent to exocyclic C=O
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j]) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8 && !aromatic[k] && tautomerClass[k] == -1) {
                                int cls = classId++;
                                tc({i, k}, cls, TW_PYRIDINONE);
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T7: Imidazole / pyrazole  aromatic N-C-N or direct N-N in ring
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                // Direct N-N bond (pyrazole / indazole)
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 7 && aromatic[j] && ring[j] && tautomerClass[j] == -1) {
                        int cls = classId++;
                        tc({i, j}, cls, TW_IMIDAZOLE_NH);
                        break;
                    }
                }
                // N-C-N through aromatic C (imidazole / benzimidazole)
                // Only match when both N atoms share a 5-membered ring (skip pyrimidine)
                if (tautomerClass[i] == -1) {
                    const auto& rings = computeRings();
                    for (int j : neighbors[i]) {
                        if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
                            for (int k : neighbors[j]) {
                                if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
                                    // Verify both N atoms share a 5-membered ring
                                    bool in5 = false;
                                    for (const auto& r : rings) {
                                        if (r.size() != 5) continue;
                                        bool hasI = false, hasK = false;
                                        for (int a : r) { if (a == i) hasI = true; if (a == k) hasK = true; }
                                        if (hasI && hasK) { in5 = true; break; }
                                    }
                                    if (!in5) continue;
                                    int cls = classId++;
                                    tc({i, k}, cls, TW_IMIDAZOLE_NH);
                                    break;
                                }
                            }
                            if (tautomerClass[i] != -1) break;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T8: Amidine  N=C-N ↔ N-C=N  (non-aromatic, symmetric, weight 0.50)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 7 && !aromatic[k]) {
                                int cls = classId++;
                                tc({i, j, k}, cls, TW_AMIDINE);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T9: Guanidine  C connected to ≥3 N  (symmetric, weight 0.50)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 6 && tautomerClass[i] == -1) {
                std::vector<int> nNbrs;
                for (int j : neighbors[i]) { if (atomicNum[j] == 7 && !aromatic[j]) nNbrs.push_back(j); }
                if (nNbrs.size() >= 3) {
                    int cls = classId++;
                    tautomerClass[i] = cls; tautomerWeight[i] = TW_GUANIDINE;
                    for (int ni : nNbrs) {
                        if (tautomerClass[ni] == -1) { tautomerClass[ni] = cls; tautomerWeight[ni] = TW_GUANIDINE; }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T10: Urea  O=C(-N)(-N)  — two N on same carbonyl C
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        std::vector<int> nn;
                        for (int k : neighbors[j]) { if (k != i && atomicNum[k] == 7) nn.push_back(k); }
                        if (nn.size() >= 2) {
                            int cls = classId++;
                            tautomerClass[i] = cls; tautomerWeight[i] = TW_UREA;
                            tautomerClass[j] = cls; tautomerWeight[j] = TW_UREA;
                            for (int ni : nn) {
                                if (tautomerClass[ni] == -1) { tautomerClass[ni] = cls; tautomerWeight[ni] = TW_UREA; }
                            }
                        }
                        break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T11: β-Enaminone  N-C=C-C=O  (vinylogous amide, weight 0.85)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
                bool found11 = false;
                for (int j : neighbors[i]) {
                    if (found11) break;
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 1) {
                        for (int k : neighbors[j]) {
                            if (found11) break;
                            if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 2) {
                                for (int m : neighbors[k]) {
                                    if (found11) break;
                                    if (m != j && atomicNum[m] == 6) {
                                        for (int o : neighbors[m]) {
                                            if (o != k && atomicNum[o] == 8 && bondOrder(m, o) == 2) {
                                                int cls = classId++;
                                                tc({i, j, k, m, o}, cls, TW_ENAMINONE);
                                                found11 = true; break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T12: Imide  C(=O)-N-C(=O)  — N flanked by two carbonyls (pKa ~9-10)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && tautomerClass[i] == -1) {
                std::vector<int> carbonyls;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) {
                                carbonyls.push_back(j); break;
                            }
                        }
                    }
                }
                if (carbonyls.size() >= 2) {
                    int cls = classId++;
                    tautomerClass[i] = cls; tautomerWeight[i] = TW_IMIDE;
                    for (int c : carbonyls) {
                        if (tautomerClass[c] == -1) { tautomerClass[c] = cls; tautomerWeight[c] = TW_IMIDE; }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T13: Hydroxamic acid  N(-OH)-C=O ↔ NH-C(=O)-OH  (pKa ~8)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && tautomerClass[i] == -1) {
                int oh = -1;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 8 && !aromatic[j] && bondOrder(i, j) == 1) { oh = j; break; }
                }
                if (oh != -1) {
                    for (int j : neighbors[i]) {
                        if (atomicNum[j] == 6) {
                            for (int k : neighbors[j]) {
                                if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) {
                                    int cls = classId++;
                                    tc({i, oh, j, k}, cls, TW_HYDROXAMIC);
                                    break;
                                }
                            }
                            if (tautomerClass[i] != -1) break;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T14: Hydroxypyrimidine / purine  exocyclic OH on arom-C with ring-N
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j] && ring[j] && bondOrder(i, j) == 1) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k]) {
                                int cls = classId++;
                                tc({i, j, k}, cls, TW_HYDROXYPYRIMIDINE);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T15: Tetrazole / triazole N-H  (aromatic ring with ≥2 adjacent N)
            //      1H/2H symmetric, weight 0.50
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                std::vector<int> adjN;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 7 && aromatic[j] && ring[j]) adjN.push_back(j);
                }
                if (adjN.size() >= 2) {
                    for (int j : adjN) {
                        if (tautomerClass[j] == -1) {
                            int cls = classId++;
                            tc({i, j}, cls, TW_TETRAZOLE_NH);
                            break;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T16: Thioamide/iminothiol  RC(=S)NR <-> RC(SH)=NR  (weight 0.95)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 16 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 7) {
                                int cls = classId++;
                                tc({i, j, k}, cls, TW_THIOAMIDE_IMINOTHIOL);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T17: Phosphonate ester  P(=O)(OH) <-> P(OH)(=O)  (symmetric, 0.50)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 15 && tautomerClass[i] == -1) {
                std::vector<int> oxygens;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 8) oxygens.push_back(j);
                }
                if (oxygens.size() >= 2) {
                    int cls = classId++;
                    tautomerClass[i] = cls; tautomerWeight[i] = TW_PHOSPHONATE;
                    for (int o : oxygens) {
                        if (tautomerClass[o] == -1) { tautomerClass[o] = cls; tautomerWeight[o] = TW_PHOSPHONATE; }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T18: DISABLED — sulfoxide S=O is a true covalent bond, not tautomerism
            // ------------------------------------------------------------------
            if (false && atomicNum[i] == 16 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 8 && bondOrder(i, j) == 2) {
                        int cls = classId++;
                        tc({i, j}, cls, TW_SULFOXIDE);
                        break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T19: Nitro/aci-nitro  C-N(=O)=O <-> C=N(=O)-OH  (weight 0.95)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
                std::vector<int> oNbrs;
                int cNbr = -1;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 8) oNbrs.push_back(j);
                    else if (atomicNum[j] == 6) cNbr = j;
                }
                if (oNbrs.size() >= 2 && cNbr != -1) {
                    int cls = classId++;
                    tautomerClass[i] = cls; tautomerWeight[i] = TW_NITRO_ACI_NITRO;
                    tautomerClass[cNbr] = cls; tautomerWeight[cNbr] = TW_NITRO_ACI_NITRO;
                    for (int o : oNbrs) {
                        if (tautomerClass[o] == -1) { tautomerClass[o] = cls; tautomerWeight[o] = TW_NITRO_ACI_NITRO; }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T20: DISABLED — nitrile/isonitrile is isomerism, not tautomerism
            // ------------------------------------------------------------------
            if (false && atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 3) {
                        int cls = classId++;
                        tc({i, j}, cls, TW_NITRILE_ISONITRILE);
                        break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T21: 1,5-keto-enol  C(=O)-C=C-C=C  (weight 0.75)
            //      5-atom keto-enol shift through conjugation
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                bool found21 = false;
                for (int j : neighbors[i]) {
                    if (found21) break;
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        for (int k : neighbors[j]) {
                            if (found21) break;
                            if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 1) {
                                for (int m : neighbors[k]) {
                                    if (found21) break;
                                    if (m != j && atomicNum[m] == 6 && bondOrder(k, m) == 2) {
                                        for (int p : neighbors[m]) {
                                            if (p != k && atomicNum[p] == 6 && !aromatic[p]) {
                                                int cls = classId++;
                                                tc({i, j, k, m, p}, cls, TW_15_KETO_ENOL);
                                                found21 = true; break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T22: Ring-chain: furanose/pyranose  (weight 0.60)
            //      O in ring bonded to C with OH neighbor (hemiacetal)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && ring[i] && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && ring[j]) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8 && !ring[k] && bondOrder(j, k) == 1) {
                                int cls = classId++;
                                tc({i, j, k}, cls, TW_FURANOSE_PYRANOSE);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T23: Ring-chain: lactol  (weight 0.65)
            //      C in ring bonded to both ring-O and exocyclic OH
            // ------------------------------------------------------------------
            if (atomicNum[i] == 6 && ring[i] && tautomerClass[i] == -1) {
                int ringO = -1, exoOH = -1;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 8 && ring[j] && !aromatic[j]) ringO = j;
                    else if (atomicNum[j] == 8 && !ring[j] && !aromatic[j] && bondOrder(i, j) == 1) exoOH = j;
                }
                if (ringO != -1 && exoOH != -1) {
                    int cls = classId++;
                    tc({ringO, i, exoOH}, cls, TW_LACTOL);
                }
            }

            // ------------------------------------------------------------------
            // T24: Phenol/quinone methide  (weight 0.92)
            //      Aromatic C-OH where para-C has exocyclic =CH2
            // ------------------------------------------------------------------
            if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j] && ring[j] && bondOrder(i, j) == 1) {
                        bool found24 = false;
                        for (int a : neighbors[j]) {
                            if (found24) break;
                            if (a != i && atomicNum[a] == 6 && aromatic[a] && ring[a]) {
                                for (int b : neighbors[a]) {
                                    if (found24) break;
                                    if (b != j && atomicNum[b] == 6 && aromatic[b] && ring[b]) {
                                        for (int c : neighbors[b]) {
                                            if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) {
                                                for (int d : neighbors[c]) {
                                                    if (atomicNum[d] == 6 && !ring[d] && bondOrder(c, d) == 2) {
                                                        int cls = classId++;
                                                        tc({i, j, c, d}, cls, TW_PHENOL_QUINONE_METHIDE);
                                                        found24 = true; break;
                                                    }
                                                }
                                                if (found24) break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T25: Triazole NH shift  1,2,3-triazole NH  (symmetric, 0.50)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                std::vector<int> adjNring;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 7 && aromatic[j] && ring[j]) adjNring.push_back(j);
                }
                if (adjNring.size() == 1) {
                    int j = adjNring[0];
                    for (int k : neighbors[j]) {
                        if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
                            int cls = classId++;
                            tc({i, j, k}, cls, TW_TRIAZOLE_NH_SHIFT);
                            break;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T26: Benzimidazole NH shift  NH-1 <-> NH-3  (symmetric, 0.50)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
                        int ringNbCount = 0;
                        for (int k : neighbors[j]) { if (ring[k] && aromatic[k]) ringNbCount++; }
                        if (ringNbCount >= 3) {
                            for (int k : neighbors[j]) {
                                if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
                                    int cls = classId++;
                                    tc({i, j, k}, cls, TW_BENZIMIDAZOLE_NH);
                                    break;
                                }
                            }
                            if (tautomerClass[i] != -1) break;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T27: Pyridone/hydroxypyridine  2-pyridinol <-> 2-pyridone  (0.85)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8 && !aromatic[k] && bondOrder(j, k) == 1) {
                                int cls = classId++;
                                tc({i, j, k}, cls, TW_PYRIDONE_HYDROXYPYRIDINE);
                                break;
                            }
                        }
                        if (tautomerClass[i] != -1) break;
                    }
                }
            }

            // ------------------------------------------------------------------
            // T28: Barbituric acid  tri-keto <-> enol  (weight 0.80)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 7 && ring[i] && tautomerClass[i] == -1) {
                std::vector<int> carbonylC;
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && ring[j]) {
                        for (int k : neighbors[j]) {
                            if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) {
                                carbonylC.push_back(j);
                                break;
                            }
                        }
                    }
                }
                if (carbonylC.size() >= 2) {
                    int cls = classId++;
                    tautomerClass[i] = cls; tautomerWeight[i] = TW_BARBITURIC;
                    for (int c : carbonylC) {
                        if (tautomerClass[c] == -1) { tautomerClass[c] = cls; tautomerWeight[c] = TW_BARBITURIC; }
                        for (int k : neighbors[c]) {
                            if (atomicNum[k] == 8 && bondOrder(c, k) == 2 && tautomerClass[k] == -1) {
                                tautomerClass[k] = cls; tautomerWeight[k] = TW_BARBITURIC;
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T29: Allyl shift  X-C=C-C <-> C=C-C-X (X=OH,SH,NH)  (0.70)
            //      Validates contiguous bond path: X(-H)--C==C--C(sp3)
            //      X must carry at least one H (checked via valence/degree)
            //      and k-m must be a single bond (not another double bond system)
            // ------------------------------------------------------------------
            if ((atomicNum[i] == 8 || atomicNum[i] == 16 || atomicNum[i] == 7)
                && !aromatic[i] && !ring[i] && tautomerClass[i] == -1) {
                // X must bear at least one H for the proton shift to occur.
                // Infer from degree: O(valence 2) with 1 heavy-atom bond has 1H,
                // S(valence 2) with 1 heavy-atom bond has 1H,
                // N(valence 3) with 1-2 heavy-atom bonds has H.
                int sumBO = 0;
                for (int nb : neighbors[i]) sumBO += bondOrder(i, nb);
                int maxValence = (atomicNum[i] == 7) ? 3 : 2; // N=3, O=2, S=2
                bool xHasH = (sumBO < maxValence + std::abs(formalCharge[i]));
                if (!xHasH) continue;  // no H on heteroatom, shift not possible
                bool found29 = false;
                for (int j : neighbors[i]) {
                    if (found29) break;
                    if (atomicNum[j] == 6 && !aromatic[j] && bondOrder(i, j) == 1) {
                        for (int k : neighbors[j]) {
                            if (found29) break;
                            if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 2) {
                                for (int m : neighbors[k]) {
                                    // Require single bond k-m (sp3 terminal carbon)
                                    if (m != j && atomicNum[m] == 6 && !aromatic[m]
                                        && bondOrder(k, m) == 1) {
                                        int cls = classId++;
                                        tc({i, j, k, m}, cls, TW_ALLYL_SHIFT);
                                        found29 = true; break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // ------------------------------------------------------------------
            // T30: Selenol/selenoketone  C=Se <-> C(-SeH)  (weight 0.94)
            // ------------------------------------------------------------------
            if (atomicNum[i] == 34 && !aromatic[i] && tautomerClass[i] == -1) {
                for (int j : neighbors[i]) {
                    if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
                        int cls = classId++;
                        tc({i, j}, cls, TW_SELENOL_SELENOKETONE);
                        for (int k : neighbors[j]) {
                            if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7) && tautomerClass[k] == -1) {
                                tautomerClass[k] = cls;
                                if (tautomerWeight[k] == 1.0f) tautomerWeight[k] = TW_SELENOL_SELENOKETONE;
                            }
                        }
                        break;
                    }
                }
            }
        }
    }

    // ========================================================================
    // Tier 2: Solvent-dependent tautomer weight corrections
    //
    // Aprotic/non-aqueous solvents shift keto-enol and related equilibria.
    // Literature pKa shifts (Bordwell 1988; Olmstead 1977):
    //   DMSO:       keto-enol +5 pKa units  -> enol more prevalent (shift -0.15)
    //               amide/lactam +3          -> imidic slightly more prevalent
    //   Methanol:   keto-enol +2             -> small shift (-0.05)
    //   Chloroform: keto-enol +8 (aprotic)   -> much more enol (-0.25)
    // ========================================================================

    /** Adjust tautomer weights for solvent effects. Call after computeTautomerClasses(). */
    void applySolventCorrection(ChemOptions::Solvent solvent) {
        if (solvent == ChemOptions::Solvent::AQUEOUS) return;

        float ketoShift  = 0.0f;
        float amideShift = 0.0f;
        float generalShift = 0.0f;

        if (solvent == ChemOptions::Solvent::DMSO) {
            ketoShift    = -0.15f;
            amideShift   = -0.08f;
            generalShift = -0.05f;
        } else if (solvent == ChemOptions::Solvent::CHLOROFORM) {
            ketoShift    = -0.25f;
            amideShift   = -0.12f;
            generalShift = -0.08f;
        } else if (solvent == ChemOptions::Solvent::METHANOL) {
            ketoShift    = -0.05f;
            amideShift   = -0.02f;
            generalShift = -0.01f;
        } else if (solvent == ChemOptions::Solvent::ACETONITRILE) {
            ketoShift    = -0.10f;
            amideShift   = -0.05f;
            generalShift = -0.03f;
        } else if (solvent == ChemOptions::Solvent::DIETHYL_ETHER) {
            ketoShift    = -0.20f;
            amideShift   = -0.10f;
            generalShift = -0.06f;
        }

        for (int i = 0; i < n; i++) {
            if (tautomerClass[i] == -1) continue;
            float w = tautomerWeight[i];
            if (w >= 1.0f) continue;

            float shift;
            if (w >= 0.96f) {
                shift = (w >= 0.97f) ? amideShift : ketoShift;
                if (std::abs(w - TW_KETO_ENOL) < 0.005f) shift = ketoShift;
            } else if (w >= 0.93f) {
                shift = amideShift;
            } else if (w >= 0.84f) {
                shift = generalShift;
            } else {
                shift = generalShift;
            }

            tautomerWeight[i] = std::max(0.05f, std::min(0.99f, w + shift));
        }
    }

    // ========================================================================
    // Morgan rank computation
    // ========================================================================

    static std::vector<int> computeMorganRanks(
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
                int h = rank[v] * HASH_PRIME;
                for (int k = 0; k < deg; k++) h = h * 31 + nbSorted[k];
                rankNew[v] = h;
            }
            rank = rankNew;

            std::unordered_set<int> newSet(rank.begin(), rank.end());
            if (static_cast<int>(newSet.size()) <= prevDistinct) break;
        }
        return rank;
    }

    // ========================================================================
    // Union-Find helpers (reuse ufFind/ufUnion defined above)
    // ========================================================================

    static std::vector<int> buildOrbits(std::vector<int>& uf, int n) {
        std::vector<int> orb(n);
        for (int i = 0; i < n; i++) orb[i] = ufFind(uf, i);
        std::vector<int> orbitMin(n, n);
        for (int i = 0; i < n; i++) {
            int r = orb[i];
            if (i < orbitMin[r]) orbitMin[r] = i;
        }
        for (int i = 0; i < n; i++) orb[i] = orbitMin[orb[i]];
        return orb;
    }

    // ========================================================================
    // 2-hop neighborhood signature for orbit disambiguation
    // ========================================================================

    /**
     * Compute a deterministic 2-hop neighborhood signature for atom v.
     *
     * The signature encodes the sorted multiset of (neighbor_label, neighbor_degree)
     * pairs for all atoms within 2 hops. Two atoms in the same automorphism orbit
     * must have identical 2-hop signatures, so splitting orbits by signature is safe.
     *
     * This refines orbits beyond what Morgan ranks alone achieve: for example,
     * in spiro compounds or fused ring systems, Morgan may converge prematurely
     * but 2-hop signatures can still distinguish structurally distinct positions.
     */
    static int64_t twoHopSignature(
            int v, const std::vector<int>& label,
            const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors) {
        // 1-hop: sorted (label, degree) of direct neighbors
        std::vector<int64_t> sig1;
        sig1.reserve(neighbors[v].size());
        for (int u : neighbors[v]) {
            sig1.push_back(static_cast<int64_t>(label[u]) * 1000003LL + degree[u]);
        }
        std::sort(sig1.begin(), sig1.end());

        // 2-hop: for each direct neighbor, the sorted (label, degree) of THEIR neighbors
        std::vector<int64_t> sig2;
        for (int u : neighbors[v]) {
            for (int w : neighbors[u]) {
                if (w == v) continue;  // skip back-link to v itself
                sig2.push_back(static_cast<int64_t>(label[w]) * 1000003LL + degree[w]);
            }
        }
        std::sort(sig2.begin(), sig2.end());

        // Combine into a single hash
        int64_t h = static_cast<int64_t>(label[v]) * HASH_PRIME + degree[v];
        h = h * HASH_PRIME + static_cast<int64_t>(sig1.size());
        for (int64_t x : sig1) h = h * 31 + x;
        h = h * HASH_PRIME + static_cast<int64_t>(sig2.size());
        for (int64_t x : sig2) h = h * 31 + x;
        return h;
    }

    /**
     * Refine an orbit partition using 2-hop neighborhood signatures.
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
            std::unordered_map<int64_t, int> sigToRep;
            for (int v : members) {
                int64_t sig = twoHopSignature(v, label, degree, neighbors);
                auto it = sigToRep.find(sig);
                if (it == sigToRep.end()) {
                    sigToRep[sig] = v;
                    refined[v] = v;
                } else {
                    refined[v] = it->second;
                }
            }
        }
        // Normalize to minimum representative in each refined orbit
        std::unordered_map<int, int> repMin;
        for (int i = 0; i < n; i++) {
            int r = refined[i];
            auto it = repMin.find(r);
            if (it == repMin.end() || i < it->second) repMin[r] = i;
        }
        for (int i = 0; i < n; i++) refined[i] = repMin[refined[i]];
        return refined;
    }

    // ========================================================================
    // Canonical labeling (Bliss-style partition refinement + individualization)
    // ========================================================================

    static int findSmallestNonSingleton(int n, const std::vector<bool>& cellEnd) {
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

    static void refinePartition(
            int n, std::vector<int>& perm, std::vector<bool>& cellEnd,
            const std::vector<std::vector<int>>& neighbors,
            const std::vector<int>& label,
            const std::vector<int>& degree, int /*targetPos*/) {
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

                    std::vector<int> count(cSize);
                    for (int ci = 0; ci < cSize; ci++) {
                        int v = perm[cStart + ci], cnt = 0;
                        for (int nb : neighbors[v]) {
                            int nbPos = inv[nb];
                            if (nbPos >= sStart && nbPos <= sEnd) cnt++;
                        }
                        count[ci] = cnt;
                    }

                    bool allSame = true;
                    for (int ci = 1; ci < cSize; ci++) {
                        if (count[ci] != count[0]) { allSame = false; break; }
                    }
                    if (allSame) { cStart = cEnd + 1; continue; }

                    // Sort indices within the cell
                    std::vector<int> idx(cSize);
                    std::iota(idx.begin(), idx.end(), 0);
                    const int fcs = cStart;
                    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
                        if (count[a] != count[b]) return count[a] < count[b];
                        int va = perm[fcs + a], vb = perm[fcs + b];
                        if (label[va] != label[vb]) return label[va] < label[vb];
                        // Use robust 2-hop signature for topological tie-breaking
                        // instead of raw node index — ensures deterministic canonical
                        // hashes regardless of input atom ordering.  (v6.8.0)
                        int64_t sigA = twoHopSignature(va, label, degree, neighbors);
                        int64_t sigB = twoHopSignature(vb, label, degree, neighbors);
                        if (sigA != sigB) return sigA < sigB;
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

    static int comparePerm(
            const std::vector<int>& perm1, const std::vector<int>& perm2,
            int n, const std::vector<int>& label,
            const std::vector<std::vector<int>>& neighbors) {
        std::vector<int> inv1(n), inv2(n);
        for (int i = 0; i < n; i++) { inv1[perm1[i]] = i; inv2[perm2[i]] = i; }
        for (int pos = 0; pos < n; pos++) {
            int v1 = perm1[pos], v2 = perm2[pos];
            int c = label[v1] - label[v2];
            if (c != 0) return c;
            const auto& nb1 = neighbors[v1];
            const auto& nb2 = neighbors[v2];
            c = static_cast<int>(nb1.size()) - static_cast<int>(nb2.size());
            if (c != 0) return c;
            std::vector<int> cn1(nb1.size()), cn2(nb2.size());
            for (size_t k = 0; k < nb1.size(); k++) cn1[k] = inv1[nb1[k]];
            for (size_t k = 0; k < nb2.size(); k++) cn2[k] = inv2[nb2[k]];
            std::sort(cn1.begin(), cn1.end());
            std::sort(cn2.begin(), cn2.end());
            for (size_t k = 0; k < cn1.size(); k++) {
                c = cn1[k] - cn2[k];
                if (c != 0) return c;
            }
        }
        return 0;
    }

    static int partialCompare(
            const std::vector<int>& perm, const std::vector<int>& bestPerm,
            int n, const std::vector<bool>& cellEnd,
            const std::vector<int>& label,
            const std::vector<std::vector<int>>& /*neighbors*/) {
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

    struct CanonResult {
        std::vector<int> canonLabel;
        std::vector<int> orbit;
        std::vector<std::vector<int>> autGenerators;
        bool generatorsTruncated = false;
    };

    static CanonResult computeCanonicalLabeling(
            int n, const std::vector<int>& label,
            const std::vector<int>& degree,
            const std::vector<std::vector<int>>& neighbors) {
        if (n == 0) return {{}, {}};

        // Morgan-like refinement for initial ordering
        std::vector<int> mRank(label.begin(), label.end());
        std::vector<int> mNew(n);
        std::vector<int> scratch;

        for (int iter = 0; iter < 8; iter++) {
            std::unordered_set<int> s(mRank.begin(), mRank.end());
            int prevDistinct = static_cast<int>(s.size());
            for (int v = 0; v < n; v++) {
                const auto& nbs = neighbors[v];
                int d = static_cast<int>(nbs.size());
                scratch.resize(d);
                for (int k = 0; k < d; k++) scratch[k] = mRank[nbs[k]];
                std::sort(scratch.begin(), scratch.end());
                int h = mRank[v] * HASH_PRIME;
                for (int k = 0; k < d; k++) h = h * 31 + scratch[k];
                mNew[v] = h;
            }
            mRank = mNew;
            std::unordered_set<int> ns(mRank.begin(), mRank.end());
            if (static_cast<int>(ns.size()) <= prevDistinct) break;
        }

        // Union-find for orbit computation
        std::vector<int> uf(n);
        std::iota(uf.begin(), uf.end(), 0);

        // Initial partition: sort atoms by (label, degree, mRank)
        std::vector<int> sorted(n);
        std::iota(sorted.begin(), sorted.end(), 0);
        std::sort(sorted.begin(), sorted.end(), [&](int a, int b) {
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

        // For large molecules, skip individualization search
        if (n > CANON_SEARCH_LIMIT) {
            int cs = 0;
            for (int i = 0; i < n; i++) {
                if (initCellEnd[i]) {
                    for (int j = cs + 1; j <= i; j++)
                        ufUnion(uf, initPerm[cs], initPerm[j]);
                    cs = i + 1;
                }
            }
            std::vector<int> cl(n);
            for (int pos = 0; pos < n; pos++) cl[initPerm[pos]] = pos;
            auto orb = buildOrbits(uf, n);
            // Refine orbits using 2-hop neighborhood signatures
            orb = refineOrbitsBy2Hop(orb, n, label, degree, neighbors);
            return { cl, orb, {}, false };
        }

        refinePartition(n, initPerm, initCellEnd, neighbors, label, degree, -1);

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
                refinePartition(n, p, ce, neighbors, label, degree, cellStart);
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
                    || comparePerm(perm, bestPerm, n, label, neighbors) < 0) {
                    bestPerm = perm;
                } else if (comparePerm(perm, bestPerm, n, label, neighbors) == 0) {
                    for (int i = 0; i < n; i++)
                        ufUnion(uf, perm[i], bestPerm[i]);
                    // Capture the automorphism generator: sigma maps
                    // bestPerm[i] -> perm[i] for each position i.
                    if (static_cast<int>(generators.size()) < MAX_GENERATORS) {
                        std::vector<int> gen(n);
                        for (int i = 0; i < n; i++) gen[bestPerm[i]] = perm[i];
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
                refinePartition(n, p, ce, neighbors, label, degree, cellStart);
                permStack.push_front(std::move(p));
                cellEndStack.push_front(std::move(ce));
                nodeCount++;
            }
        }

        if (bestPerm.empty()) bestPerm = initPerm;

        std::vector<int> cl(n);
        for (int pos = 0; pos < n; pos++) cl[bestPerm[pos]] = pos;

        // Merge initial partition orbits
        int cs = 0;
        for (int i = 0; i < n; i++) {
            if (initCellEnd[i]) {
                for (int j = cs + 1; j <= i; j++)
                    ufUnion(uf, initPerm[cs], initPerm[j]);
                cs = i + 1;
            }
        }
        auto orb = buildOrbits(uf, n);
        // Refine orbits using 2-hop neighborhood signatures
        orb = refineOrbitsBy2Hop(orb, n, label, degree, neighbors);
        return { cl, orb, std::move(generators), generatorsTruncated };
    }

    // ========================================================================
    // Canonical hash
    // ========================================================================

    static uint64_t computeCanonicalHash(
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
                 + static_cast<uint64_t>(vecAt(g.tetraChirality, v));
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
                        + static_cast<uint64_t>(g.dbStereo(v, nb));
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

    // ========================================================================
    // Canonical SMILES generation
    // ========================================================================

    static const std::unordered_map<int, std::string>& organicSymbols() {
        static const std::unordered_map<int, std::string> m = {
            {5,"B"},{6,"C"},{7,"N"},{8,"O"},{15,"P"},{16,"S"},
            {9,"F"},{17,"Cl"},{35,"Br"},{53,"I"}
        };
        return m;
    }

    static const char* elementSymbol(int z) {
        static const char* const syms[] = {
            "?","H","He","Li","Be","B","C","N","O","F","Ne",
            "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
            "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
            "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
            "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
            "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
            "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
            "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
            "Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
            "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
            "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
            "Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
        };
        constexpr int count = static_cast<int>(sizeof(syms) / sizeof(syms[0]));
        return (z >= 0 && z < count) ? syms[z] : "?";
    }

    std::vector<int> sortedNeighborsByCanonical(int atom) const {
        auto sorted = neighbors[atom];
        std::sort(sorted.begin(), sorted.end(), [this](int a, int b) {
            return canonicalLabel[a] < canonicalLabel[b];
        });
        return sorted;
    }

    // Stereo-aware atom writer. Emits @/@@  inside bracket when tetraChirality is set.
    // parent and sortedChildren are reserved for future parity-aware @ ordering.
    void writeAtom(int atom, int /*parent*/,
                   const std::vector<int>& /*sortedChildren*/,
                   std::string& sb) const {
        int z = atomicNum[atom];
        int charge = !formalCharge.empty() ? formalCharge[atom] : 0;
        bool arom = !aromatic.empty() && aromatic[atom];
        int chiral = !tetraChirality.empty() ? tetraChirality[atom] : 0;
        bool hasStereo = (chiral != 0);
        auto& org = organicSymbols();
        auto it = org.find(z);
        if (it != org.end() && charge == 0 && !hasStereo) {
            const std::string& sym = it->second;
            if (arom)
                sb += static_cast<char>(
                    std::tolower(static_cast<unsigned char>(sym[0])));
            else
                sb += sym;
        } else {
            sb += '[';
            const char* elemSym = elementSymbol(z);
            if (arom) {
                sb += static_cast<char>(
                    std::tolower(static_cast<unsigned char>(elemSym[0])));
                if (elemSym[1] != '\0') sb += &elemSym[1];
            } else {
                sb += elemSym;
            }
            if (hasStereo) {
                // chiral==1 → CLOCKWISE → @@, chiral==2 → ANTICLOCKWISE → @
                sb += (chiral == 1 ? "@@" : "@");
            }
            if (charge > 0) {
                sb += '+';
                if (charge > 1) sb += std::to_string(charge);
            } else if (charge < 0) {
                sb += '-';
                if (charge < -1) sb += std::to_string(-charge);
            }
            sb += ']';
        }
    }

    // Stereo-aware bond writer. Emits /\ for single bonds adjacent to stereo double bonds.
    void writeBondSymbol(int from, int to, std::string& sb) const {
        int order = bondOrder(from, to);
        bool fromArom = !aromatic.empty() && aromatic[from];
        bool toArom   = !aromatic.empty() && aromatic[to];
        if (order == 4 || (fromArom && toArom)) return;
        if (order == 2) { sb += '='; return; }
        if (order == 3) { sb += '#'; return; }
        if (order == 1 && (fromArom != toArom)) { sb += '-'; return; }
        // single bond: check E/Z stereo annotation
        if (order == 1 && !dbStereoConf.empty()) {
            int conf = dbStereo(from, to);
            if (conf != 0) {
                sb += (conf == 1 ? '/' : '\\');
            }
        }
    }

    static void writeRingDigit(int rid, std::string& sb) {
        if (rid < 10) {
            sb += static_cast<char>('0' + rid);
        } else {
            sb += '%';
            sb += std::to_string(rid);
        }
    }

    // Ring closure info: {ringId, otherAtom}
    using RingClosureList = std::vector<std::vector<std::array<int,2>>>;

    void preScanDfs(int atom, int parent, std::vector<int>& color,
                     RingClosureList& rc,
                     std::unordered_set<int64_t>& ringEdges,
                     int& nextRing) const {
        color[atom] = 1;
        auto sorted = sortedNeighborsByCanonical(atom);
        int skippedParent = 0;
        for (int nb : sorted) {
            if (nb == parent && skippedParent == 0) {
                skippedParent = 1; continue;
            }
            if (color[nb] == 0) {
                preScanDfs(nb, atom, color, rc, ringEdges, nextRing);
            } else if (color[nb] == 1) {
                int rid = nextRing++;
                rc[nb].push_back({rid, atom});
                rc[atom].push_back({rid, nb});
                ringEdges.insert(bondKey(atom, nb));
            }
        }
        color[atom] = 2;
    }

    void emitSmiles(int atom, int parent, std::vector<bool>& visited,
                     std::string& sb, const RingClosureList& rc,
                     const std::unordered_set<int64_t>& ringEdges) const {
        visited[atom] = true;
        auto sorted = sortedNeighborsByCanonical(atom);
        writeAtom(atom, parent, sorted, sb);
        for (const auto& entry : rc[atom]) {
            if (!visited[entry[1]]) writeBondSymbol(atom, entry[1], sb);
            writeRingDigit(entry[0], sb);
        }
        int childCount = 0;
        for (int nb : sorted) {
            if (!visited[nb]
                && ringEdges.find(bondKey(atom, nb)) == ringEdges.end())
                childCount++;
        }
        int childIdx = 0;
        for (int nb : sorted) {
            if (visited[nb]
                || ringEdges.find(bondKey(atom, nb)) != ringEdges.end())
                continue;
            childIdx++;
            bool needParen = childCount > 1 && childIdx < childCount;
            if (needParen) sb += '(';
            writeBondSymbol(atom, nb, sb);
            emitSmiles(nb, atom, visited, sb, rc, ringEdges);
            if (needParen) sb += ')';
        }
    }

    std::string toCanonicalSmiles() const {
        if (n == 0) return "";
        ensureCanonical(); // ensure canonicalLabel is computed
        std::vector<int> inv(n);
        for (int i = 0; i < n; i++) inv[canonicalLabel[i]] = i;

        std::vector<int> color(n, 0);
        RingClosureList rc(n);
        std::unordered_set<int64_t> ringEdges;
        int nextRing = 1;

        for (int ci = 0; ci < n; ci++) {
            int start = inv[ci];
            if (color[start] == 0)
                preScanDfs(start, -1, color, rc, ringEdges, nextRing);
        }

        std::vector<bool> visited(n, false);
        std::string sb;
        for (int ci = 0; ci < n; ci++) {
            int start = inv[ci];
            if (visited[start]) continue;
            if (!sb.empty()) sb += '.';
            emitSmiles(start, -1, visited, sb, rc, ringEdges);
        }
        return sb;
    }

    // ========================================================================
    // NLF (Neighbor Label Frequency) helpers
    // ========================================================================

    static std::vector<int> freqMapToSortedArray(
            const std::unordered_map<int,int>& freq) {
        if (freq.empty()) return {};
        std::vector<std::pair<int,int>> pairs(freq.begin(), freq.end());
        std::sort(pairs.begin(), pairs.end());
        std::vector<int> arr(pairs.size() * 2);
        for (size_t i = 0; i < pairs.size(); i++) {
            arr[i * 2]     = pairs[i].first;
            arr[i * 2 + 1] = pairs[i].second;
        }
        return arr;
    }

    // NLF labels use atomicNum + aromaticity but NOT ring membership.
    // Ring compatibility is enforced separately in atom/bond checks so the
    // cached NLF tables remain reusable across option profiles.
    static int nlfLabel(const MolGraph& g, int idx) {
        return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
    }

    static std::vector<int> buildNLF1(const MolGraph& g, int idx) {
        std::unordered_map<int,int> freq;
        for (int nb : g.neighbors[idx]) freq[nlfLabel(g, nb)]++;
        return freqMapToSortedArray(freq);
    }

    static std::vector<int> buildNLF2(const MolGraph& g, int idx) {
        std::unordered_map<int,int> freq;
        std::vector<bool> direct(g.n, false);
        for (int nb : g.neighbors[idx]) direct[nb] = true;
        std::vector<bool> seen(g.n, false);
        for (int nb : g.neighbors[idx]) {
            for (int j : g.neighbors[nb]) {
                if (j == idx || direct[j] || seen[j]) continue;
                seen[j] = true;
                freq[nlfLabel(g, j)]++;
            }
        }
        return freqMapToSortedArray(freq);
    }

    static std::vector<int> buildNLF3(const MolGraph& g, int idx) {
        std::unordered_map<int,int> freq;
        std::vector<bool> level1(g.n, false);
        for (int nb : g.neighbors[idx]) level1[nb] = true;

        std::vector<bool> level2(g.n, false);
        for (int i = 0; i < g.n; i++) {
            if (!level1[i]) continue;
            for (int j : g.neighbors[i])
                if (j != idx && !level1[j]) level2[j] = true;
        }
        std::vector<bool> level3(g.n, false);
        for (int i = 0; i < g.n; i++) {
            if (!level2[i]) continue;
            for (int j : g.neighbors[i])
                if (j != idx && !level1[j] && !level2[j]) level3[j] = true;
        }
        for (int j = 0; j < g.n; j++)
            if (level3[j]) freq[nlfLabel(g, j)]++;
        return freqMapToSortedArray(freq);
    }

    using NLFBuilder = std::vector<int>(*)(const MolGraph&, int);

    static std::vector<std::vector<int>> buildAllNLF(
            const MolGraph& g, NLFBuilder builder) {
        std::vector<std::vector<int>> nlf(g.n);
        for (int i = 0; i < g.n; i++) nlf[i] = builder(g, i);
        return nlf;
    }

    const std::vector<std::vector<int>>& getNLF1() const {
        if (cachedNLF1.empty() && n > 0)
            cachedNLF1 = buildAllNLF(*this, buildNLF1);
        return cachedNLF1;
    }

    const std::vector<std::vector<int>>& getNLF2() const {
        if (cachedNLF2.empty() && n > 0)
            cachedNLF2 = buildAllNLF(*this, buildNLF2);
        return cachedNLF2;
    }

    const std::vector<std::vector<int>>& getNLF3() const {
        if (cachedNLF3.empty() && n > 0)
            cachedNLF3 = buildAllNLF(*this, buildNLF3);
        return cachedNLF3;
    }

    static bool nlfOk(const std::vector<int>& fq, const std::vector<int>& ft) {
        size_t fi = 0, ti = 0;
        while (fi < fq.size()) {
            int fLabel = fq[fi], fFreq = fq[fi + 1];
            while (ti < ft.size() && ft[ti] < fLabel) ti += 2;
            if (ti >= ft.size() || ft[ti] != fLabel || ft[ti + 1] < fFreq)
                return false;
            fi += 2;
        }
        return true;
    }

    static bool nlfCheckOk(
            int qi, int tj,
            const std::vector<std::vector<int>>& qNLF1,
            const std::vector<std::vector<int>>& tNLF1,
            const std::vector<std::vector<int>>& qNLF2,
            const std::vector<std::vector<int>>& tNLF2,
            const std::vector<std::vector<int>>& qNLF3,
            const std::vector<std::vector<int>>& tNLF3,
            bool useTwoHop, bool useThreeHop) {
        if (!nlfOk(qNLF1[qi], tNLF1[tj])) return false;
        if (useTwoHop && !nlfOk(qNLF2[qi], tNLF2[tj])) return false;
        return !useThreeHop || nlfOk(qNLF3[qi], tNLF3[tj]);
    }

    // ========================================================================
    // initDerivedFields -- called after construction to compute all
    // derived data (adjLong, morgan, canonical label, orbits, tautomers, etc.)
    // ========================================================================

    void initDerivedFields() {
        words = (n + 63) >> 6;
        adjLong.assign(n, std::vector<uint64_t>(words, 0));
        for (int i = 0; i < n; i++) {
            for (int k : neighbors[i])
                adjLong[i][k >> 6] |= uint64_t(1) << (k & 63);
        }
        computeTautomerClasses();
        applySolventCorrection(solvent_);
        adjustWeightsForPH(pH_);
        // Canonical labeling and ring counts are deferred; call
        // ensureCanonical() / ensureRingCounts() before accessing those fields.
    }

    /** Lazily compute Morgan ranks, canonical labeling, orbit, and canonical hash.
     *  NOT thread-safe -- callers sharing a MolGraph across threads must
     *  call ensureCanonical() once before spawning threads. */
    void ensureCanonical() const {
        if (canonicalComputed_) return;
        morganRank = computeMorganRanks(n, label, neighbors);
        auto clResult = computeCanonicalLabeling(n, label, degree, neighbors);
        canonicalLabel = std::move(clResult.canonLabel);
        orbit = std::move(clResult.orbit);
        autGenerators_ = std::move(clResult.autGenerators);
        autGeneratorsTruncated_ = clResult.generatorsTruncated;
        canonicalHash = computeCanonicalHash(*this, canonicalLabel);
        canonicalComputed_ = true;
    }

    /** Lazily compute per-atom ring counts.
     *  NOT thread-safe -- callers sharing a MolGraph across threads must
     *  call ensureRingCounts() once before spawning threads. */
    void ensureRingCounts() const {
        if (ringCountsComputed_) return;
        computeRingCounts();
        ringCountsComputed_ = true;
    }

    /**
     * Cached neighbors sorted by descending degree for VF2++ query ordering.
     * Computed once, reused across all matcher constructions.
     * @since 6.5.3
     */
    const std::vector<std::vector<int>>& getNeighborsByDegDesc() const {
        if (cachedNeighborsByDegDesc.empty() && n > 0) {
            cachedNeighborsByDegDesc.resize(n);
            for (int i = 0; i < n; ++i) {
                auto nb = neighbors[i]; // copy
                std::sort(nb.begin(), nb.end(), [this](int a, int b) {
                    return degree[a] > degree[b]; // descending
                });
                cachedNeighborsByDegDesc[i] = std::move(nb);
            }
        }
        return cachedNeighborsByDegDesc;
    }

    /**
     * Cached pharmacophore feature classification per atom for FCFP.
     * Computed once per MolGraph, eliminates O(n×degree²) recomputation.
     * @since 6.5.3
     */
    const std::vector<int>& getPharmacophoreFeatures() const {
        if (!pharmComputed_) {
            cachedPharmacophoreFeatures.resize(n);
            for (int i = 0; i < n; ++i) {
                int z = atomicNum[i];
                int charge = formalCharge[i];
                bool arom = (aromatic[i] != 0);
                int deg = degree[i];
                int features = 0;

                // H-bond donor: N-H, O-H, S-H
                if (z == 7 || z == 8 || z == 16) {
                    int typVal = (z == 7) ? 3 : 2;
                    if (deg < typVal && charge >= 0) features |= 1;
                }

                // H-bond acceptor: N (not pyrrole), O, F, S (not thiophene)
                bool isPyrroleN = false;
                if (z == 7 && arom) {
                    int boSum = 0;
                    for (int nb : neighbors[i]) boSum += bondOrder(i, nb);
                    int hc = std::max(0, 3 - boSum - std::abs(charge));
                    isPyrroleN = (hc > 0);
                }
                if ((z == 7 && charge <= 0 && !isPyrroleN) || z == 8 || z == 9 ||
                    (z == 16 && !arom))
                    features |= 2;

                // Positive ionisable: basic amines (not amide, not aniline)
                if (z == 7 && charge > 0) features |= 4;
                if (z == 7 && !arom && deg <= 3 && charge == 0) {
                    bool isAmide = false, isAniline = false;
                    for (int nb : neighbors[i]) {
                        if (atomicNum[nb] == 6) {
                            if (aromatic[nb]) { isAniline = true; break; }
                            for (int nb2 : neighbors[nb]) {
                                if (nb2 != i && atomicNum[nb2] == 8 && bondOrder(nb, nb2) == 2) {
                                    isAmide = true; break;
                                }
                            }
                        }
                        if (isAmide) break;
                    }
                    if (!isAmide && !isAniline) features |= 4;
                }

                // Negative ionisable: carboxylic, phosphoric, sulfonic
                if (z == 8 && charge < 0) features |= 8;
                if (z == 8 && charge == 0 && deg == 1) {
                    for (int nb : neighbors[i]) {
                        int nbZ = atomicNum[nb];
                        if (nbZ == 6 && bondOrder(i, nb) == 1) {
                            for (int nb2 : neighbors[nb]) {
                                if (nb2 != i && atomicNum[nb2] == 8 && bondOrder(nb, nb2) == 2)
                                    features |= 8;
                            }
                        }
                        if (nbZ == 15 || nbZ == 16) features |= 8;
                    }
                }

                // Aromatic
                if (arom) features |= 16;

                // Hydrophobic: non-aromatic C (no hetero neighbors), or halogens
                if (z == 6 && !arom) {
                    bool hasHetero = false;
                    for (int nb : neighbors[i]) {
                        int nbZ = atomicNum[nb];
                        if (nbZ != 6 && nbZ != 1) { hasHetero = true; break; }
                    }
                    if (!hasHetero) features |= 32;
                }
                if (z == 17 || z == 35 || z == 53) features |= 32;

                cachedPharmacophoreFeatures[i] = features;
            }
            pharmComputed_ = true;
        }
        return cachedPharmacophoreFeatures;
    }

    // ========================================================================
    // Builder -- CDK-free molecule construction via fluent API
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
        Builder& tetrahedralChirality(std::vector<int> v)                 { tetraChirality_ = std::move(v); return *this; }
        Builder& doubleBondStereo(std::vector<std::vector<int>> v)        { dbStereoConf_ = std::move(v); return *this; }
        /** Optional external atom IDs (1-based, gaps allowed). Size must equal atomCount if provided. */
        Builder& atomIds(std::vector<int> v)                               { atomIds_ = std::move(v); return *this; }
        Builder& name(std::string v)                                        { name_ = std::move(v); return *this; }
        Builder& programLine(std::string v)                                 { programLine_ = std::move(v); return *this; }
        Builder& comment(std::string v)                                     { comment_ = std::move(v); return *this; }
        Builder& properties(std::map<std::string, std::string> v)           { properties_ = std::move(v); return *this; }
        Builder& solvent(ChemOptions::Solvent s)                            { solvent_ = s; return *this; }
        Builder& pH(double v)                                                { pH_ = v; return *this; }

        MolGraph build(
            bool perceiveAromaticity = true,
            AromaticityModel model = AromaticityModel::DAYLIGHT_LIKE) const {
            if (n_ < 0)
                throw std::invalid_argument("atomCount must be >= 0");
            if (static_cast<int>(atomicNum_.size()) != n_)
                throw std::invalid_argument("atomicNumbers must have length == atomCount");
            if (static_cast<int>(neighbors_.size()) != n_)
                throw std::invalid_argument("neighbors must have length == atomCount");
            if (!formalCharge_.empty() && static_cast<int>(formalCharge_.size()) != n_)
                throw std::invalid_argument("formalCharges must have length == atomCount");
            if (!massNumber_.empty() && static_cast<int>(massNumber_.size()) != n_)
                throw std::invalid_argument("massNumbers must have length == atomCount");
            if (!hydrogenCount_.empty() && static_cast<int>(hydrogenCount_.size()) != n_)
                throw std::invalid_argument("hydrogenCounts must have length == atomCount");
            if (!atomClass_.empty() && static_cast<int>(atomClass_.size()) != n_)
                throw std::invalid_argument("atomClasses must have length == atomCount");
            if (!ring_.empty() && static_cast<int>(ring_.size()) != n_)
                throw std::invalid_argument("ringFlags must have length == atomCount");
            if (!aromatic_.empty() && static_cast<int>(aromatic_.size()) != n_)
                throw std::invalid_argument("aromaticFlags must have length == atomCount");
            if (!tetraChirality_.empty() && static_cast<int>(tetraChirality_.size()) != n_)
                throw std::invalid_argument("tetrahedralChirality must have length == atomCount");
            if (!bondOrders_.empty() && static_cast<int>(bondOrders_.size()) != n_)
                throw std::invalid_argument("bondOrders must have length == atomCount");
            if (!bondRings_.empty() && static_cast<int>(bondRings_.size()) != n_)
                throw std::invalid_argument("bondRingFlags must have length == atomCount");
            if (!bondAroms_.empty() && static_cast<int>(bondAroms_.size()) != n_)
                throw std::invalid_argument("bondAromaticFlags must have length == atomCount");
            if (!dbStereoConf_.empty() && static_cast<int>(dbStereoConf_.size()) != n_)
                throw std::invalid_argument("doubleBondStereo must have length == atomCount");

            for (int i = 0; i < n_; ++i) {
                if (atomicNum_[i] < 0 || atomicNum_[i] > 118)
                    throw std::invalid_argument("atomicNumbers must be in [0, 118]");
                if (!massNumber_.empty() && massNumber_[i] < 0)
                    throw std::invalid_argument("massNumbers must be >= 0");
                if (!hydrogenCount_.empty() && hydrogenCount_[i] < 0)
                    throw std::invalid_argument("hydrogenCounts must be >= 0");
                if (!atomClass_.empty() && atomClass_[i] < 0)
                    throw std::invalid_argument("atomClasses must be >= 0");
                if (!ring_.empty() && ring_[i] > 1)
                    throw std::invalid_argument("ringFlags entries must be 0 or 1");
                if (!aromatic_.empty() && aromatic_[i] > 1)
                    throw std::invalid_argument("aromaticFlags entries must be 0 or 1");
                if (!tetraChirality_.empty()
                    && tetraChirality_[i] != 0
                    && tetraChirality_[i] != 1
                    && tetraChirality_[i] != 2)
                    throw std::invalid_argument("tetrahedralChirality entries must be 0, 1, or 2");
            }

            std::vector<std::unordered_map<int, int>> neighborIndex(n_);
            for (int i = 0; i < n_; ++i) {
                neighborIndex[i].reserve(neighbors_[i].size());
                for (int k = 0; k < static_cast<int>(neighbors_[i].size()); ++k) {
                    int nb = neighbors_[i][k];
                    if (nb < 0 || nb >= n_)
                        throw std::invalid_argument("neighbor index out of range");
                    if (nb == i)
                        throw std::invalid_argument("self-loop detected in neighbors");
                    if (!neighborIndex[i].emplace(nb, k).second)
                        throw std::invalid_argument("duplicate neighbor detected in adjacency list");
                }
            }

            auto classifyIntRows = [&](const std::vector<std::vector<int>>& rows,
                                       const char* name) {
                std::vector<RowEncoding> enc(n_, RowEncoding::ABSENT);
                if (rows.empty()) return enc;
                for (int i = 0; i < n_; ++i) {
                    const auto& row = rows[i];
                    const auto deg = neighbors_[i].size();
                    if (row.size() == deg) enc[i] = RowEncoding::PARALLEL;
                    else if (row.size() == static_cast<size_t>(n_)) enc[i] = RowEncoding::DENSE;
                    else
                        throw std::invalid_argument(std::string(name)
                            + " rows must have length == neighbor count or atomCount");
                }
                return enc;
            };

            auto classifyBoolRows = [&](const std::vector<std::vector<bool>>& rows,
                                        const char* name) {
                std::vector<RowEncoding> enc(n_, RowEncoding::ABSENT);
                if (rows.empty()) return enc;
                for (int i = 0; i < n_; ++i) {
                    const auto& row = rows[i];
                    const auto deg = neighbors_[i].size();
                    if (row.size() == deg) enc[i] = RowEncoding::PARALLEL;
                    else if (row.size() == static_cast<size_t>(n_)) enc[i] = RowEncoding::DENSE;
                    else
                        throw std::invalid_argument(std::string(name)
                            + " rows must have length == neighbor count or atomCount");
                }
                return enc;
            };

            auto bondOrderEnc = classifyIntRows(bondOrders_, "bondOrders");
            auto bondRingEnc = classifyBoolRows(bondRings_, "bondRingFlags");
            auto bondAromEnc = classifyBoolRows(bondAroms_, "bondAromaticFlags");

            auto intPropAt = [&](const std::vector<std::vector<int>>& rows,
                                 const std::vector<RowEncoding>& enc,
                                 int i, int k, int j, int fallback) {
                if (rows.empty()) return fallback;
                return enc[i] == RowEncoding::DENSE ? rows[i][j] : rows[i][k];
            };
            auto boolPropAt = [&](const std::vector<std::vector<bool>>& rows,
                                  const std::vector<RowEncoding>& enc,
                                  int i, int k, int j, bool fallback) -> bool {
                if (rows.empty()) return fallback;
                return enc[i] == RowEncoding::DENSE ? rows[i][j] : rows[i][k];
            };

            for (int i = 0; i < n_; ++i) {
                for (const auto& [nb, k] : neighborIndex[i]) {
                    auto rev = neighborIndex[nb].find(i);
                    if (rev == neighborIndex[nb].end())
                        throw std::invalid_argument("adjacency must be symmetric");
                    if (i > nb) continue;
                    int rk = rev->second;
                    if (!bondOrders_.empty()) {
                        int forward = intPropAt(bondOrders_, bondOrderEnc, i, k, nb, 1);
                        int backward = intPropAt(bondOrders_, bondOrderEnc, nb, rk, i, 1);
                        if (forward <= 0 || backward <= 0)
                            throw std::invalid_argument("bondOrders for present edges must be > 0");
                        if (forward != backward)
                            throw std::invalid_argument("bondOrders must be symmetric");
                    }
                    if (!bondRings_.empty()) {
                        bool forward = boolPropAt(bondRings_, bondRingEnc, i, k, nb, false);
                        bool backward = boolPropAt(bondRings_, bondRingEnc, nb, rk, i, false);
                        if (forward != backward)
                            throw std::invalid_argument("bondRingFlags must be symmetric");
                    }
                    if (!bondAroms_.empty()) {
                        bool forward = boolPropAt(bondAroms_, bondAromEnc, i, k, nb, false);
                        bool backward = boolPropAt(bondAroms_, bondAromEnc, nb, rk, i, false);
                        if (forward != backward)
                            throw std::invalid_argument("bondAromaticFlags must be symmetric");
                    }
                }
            }

            if (!dbStereoConf_.empty()) {
                for (int i = 0; i < n_; ++i) {
                    if (static_cast<int>(dbStereoConf_[i].size()) != n_)
                        throw std::invalid_argument("doubleBondStereo rows must have length == atomCount");
                }
                for (int i = 0; i < n_; ++i) {
                    for (int j = 0; j < n_; ++j) {
                        int conf = dbStereoConf_[i][j];
                        if (conf < 0 || conf > 2)
                            throw std::invalid_argument("doubleBondStereo entries must be 0, 1, or 2");
                        if (i == j && conf != 0)
                            throw std::invalid_argument("doubleBondStereo diagonal must be zero");
                        if (j <= i) continue;
                        if (conf != dbStereoConf_[j][i])
                            throw std::invalid_argument("doubleBondStereo must be symmetric");
                        if (conf != 0) {
                            auto it = neighborIndex[i].find(j);
                            if (it == neighborIndex[i].end())
                                throw std::invalid_argument("doubleBondStereo requires a bond");
                            int ord = intPropAt(bondOrders_, bondOrderEnc, i, it->second, j, 1);
                            if (ord != 2)
                                throw std::invalid_argument("doubleBondStereo requires bond order 2");
                        }
                    }
                }
            }

            if (!atomIds_.empty()) {
                if (static_cast<int>(atomIds_.size()) != n_)
                    throw std::invalid_argument("atomIds must have length == atomCount");
                std::unordered_set<int> seenIds;
                seenIds.reserve(atomIds_.size());
                for (int id : atomIds_) {
                    if (id <= 0)
                        throw std::invalid_argument("atomIds must be positive");
                    if (!seenIds.insert(id).second)
                        throw std::invalid_argument("atomIds must be unique");
                }
            }

            MolGraph g;
            g.n = n_;
            g.words = (n_ + 63) >> 6;

            // Copy atom arrays
            g.atomicNum    = atomicNum_;
            g.formalCharge = formalCharge_.empty()
                             ? std::vector<int>(n_, 0) : formalCharge_;
            g.massNumber   = massNumber_.empty()
                             ? std::vector<int>(n_, 0) : massNumber_;
            g.hydrogenCount = hydrogenCount_.empty()
                             ? std::vector<int>(n_, 0) : hydrogenCount_;
            g.atomClass    = atomClass_.empty()
                             ? std::vector<int>(n_, 0) : atomClass_;

            if (ring_.empty())     g.ring.assign(n_, uint8_t(0));     else g.ring = ring_;
            if (aromatic_.empty()) g.aromatic.assign(n_, uint8_t(0)); else g.aromatic = aromatic_;

            // Neighbors + degree
            g.neighbors = neighbors_;
            g.degree.resize(n_);
            for (int i = 0; i < n_; i++)
                g.degree[i] = static_cast<int>(g.neighbors[i].size());

            // Bond storage -- dense or sparse
            if (n_ <= SPARSE_THRESHOLD) {
                g.useDense = true;
                g.bondOrdMatrix.assign(n_, std::vector<int>(n_, 0));
                g.bondRingMatrix.assign(n_, std::vector<bool>(n_, false));
                g.bondAromMatrix.assign(n_, std::vector<bool>(n_, false));
                for (int i = 0; i < n_; i++) {
                    const auto& nbs = g.neighbors[i];
                    for (int k = 0; k < static_cast<int>(nbs.size()); k++) {
                        int j = nbs[k];
                        int ord = intPropAt(bondOrders_, bondOrderEnc, i, k, j, 1);
                        bool inRing = boolPropAt(bondRings_, bondRingEnc, i, k, j, false);
                        bool arom = boolPropAt(bondAroms_, bondAromEnc, i, k, j, false);
                        g.bondOrdMatrix[i][j]  = ord;
                        g.bondRingMatrix[i][j] = inRing;
                        g.bondAromMatrix[i][j] = arom;
                    }
                }
            } else {
                g.useDense = false;
                for (int i = 0; i < n_; i++) {
                    const auto& nbs = g.neighbors[i];
                    for (int k = 0; k < static_cast<int>(nbs.size()); k++) {
                        int j = nbs[k];
                        if (i >= j) continue;   // store each edge once
                        int ord = intPropAt(bondOrders_, bondOrderEnc, i, k, j, 1);
                        int inRing = boolPropAt(bondRings_, bondRingEnc, i, k, j, false) ? 1 : 0;
                        int arom = boolPropAt(bondAroms_, bondAromEnc, i, k, j, false) ? 1 : 0;
                        g.sparseBondProps[bondKey(i, j)] = {ord, inRing, arom};
                    }
                }
            }

            // Stereochemistry
            g.tetraChirality = tetraChirality_.empty()
                               ? std::vector<int>(n_, 0) : tetraChirality_;
            if (n_ <= SPARSE_THRESHOLD) {
                g.hasDbStereo = true;
                g.dbStereoConf.assign(n_, std::vector<int>(n_, 0));
                if (!dbStereoConf_.empty()) {
                    for (int i = 0; i < n_
                         && i < static_cast<int>(dbStereoConf_.size()); i++) {
                        for (int j = 0; j < n_
                             && j < static_cast<int>(dbStereoConf_[i].size()); j++) {
                            g.dbStereoConf[i][j] = dbStereoConf_[i][j];
                        }
                    }
                }
            } else if (!dbStereoConf_.empty()) {
                g.hasDbStereo = true;
                for (int i = 0; i < n_; ++i) {
                    for (int j = i + 1; j < n_; ++j) {
                        int conf = dbStereoConf_[i][j];
                        if (conf != 0)
                            g.sparseDbStereo[bondKey(i, j)] = conf;
                    }
                }
            }

            // External atom IDs (optional)
            if (!atomIds_.empty()) g.atomId = atomIds_;

            g.name = name_;
            g.programLine = programLine_;
            g.comment = comment_;
            g.properties = properties_;

            if (perceiveAromaticity) g.normalizeRingAndAromaticity(model);
            g.refreshAtomLabels();

            g.ringCount.assign(n_, 0);
            g.solvent_ = solvent_;
            g.initDerivedFields();
            return g;
        }
    };
};

// ============================================================================
// ChemOps -- atom and bond compatibility checks (static utility)
// ============================================================================

struct ChemOps {

    /// Atom compatibility check for substructure / MCS matching.
    static bool atomsCompatible(
            const MolGraph& gq, int qi,
            const MolGraph& gt, int tj,
            const ChemOptions& C) {
        // Tautomer-aware: relaxed matching for C/N/O/S in tautomeric positions
        if (C.tautomerAware
            && !gq.tautomerClass.empty() && !gt.tautomerClass.empty()
            && gq.tautomerClass[qi] != -1 && gt.tautomerClass[tj] != -1) {
            int aq = gq.atomicNum[qi], at = gt.atomicNum[tj];
            auto isTautElem = [](int z) { return z==6||z==7||z==8||z==16; };
            if (isTautElem(aq) && isTautElem(at)) {
                if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj])
                    return false;
                return true;
            }
        }
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
        if (C.useChirality) {
            int qs = gq.tetraChirality[qi], ts = gt.tetraChirality[tj];
            if (qs != 0 && ts != 0 && qs != ts) return false;
        }
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

        // Tautomer-aware: both endpoints tautomeric => any order matches
        if (C.tautomerAware
            && !g1.tautomerClass.empty() && !g2.tautomerClass.empty()) {
            bool qBothTaut = g1.tautomerClass[qi] != -1
                          && g1.tautomerClass[qk] != -1;
            bool tBothTaut = g2.tautomerClass[tj] != -1
                          && g2.tautomerClass[tk] != -1;
            if (qBothTaut && tBothTaut) return true;
        }

        // Strict aromaticity: check bond aromaticity directly (no ring guard)
        if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
            && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk))
            return false;

        if (C.useBondStereo) {
            int qConf = g1.dbStereo(qi, qk), tConf = g2.dbStereo(tj, tk);
            if (qConf != 0 && tConf != 0 && qConf != tConf) return false;
        }

        if (C.ringMatchesRingOnly
            && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk))
            return false;

        if (C.matchBondOrder == ChemOptions::BondOrderMode::ANY)
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
// computeTautomerConfidence -- post-MCS tautomer match quality score ∈ (0,1].
//
// For each atom pair (qi→ti) in the MCS that is a cross-element tautomeric
// match (atomicNum differs but both atoms are in a tautomeric class), the
// confidence is multiplied by min(w_query, w_target).  Non-tautomeric pairs
// contribute 1.0 (no penalty).  Returns 1.0 for an empty mapping or when
// neither graph has tautomeric atoms.
// ===========================================================================

inline double computeTautomerConfidence(
        const MolGraph& gq,
        const MolGraph& gt,
        const std::unordered_map<int,int>& mcs) {
    if (mcs.empty()) return 1.0;
    double conf = 1.0;
    for (const auto& kv : mcs) {
        int qi = kv.first, ti = kv.second;
        if (!gq.tautomerClass.empty() && !gt.tautomerClass.empty()
            && gq.tautomerClass[qi] >= 0 && gt.tautomerClass[ti] >= 0
            && gq.atomicNum[qi] != gt.atomicNum[ti]) {
            float wq = gq.tautomerWeight.empty() ? 1.0f : gq.tautomerWeight[qi];
            float wt = gt.tautomerWeight.empty() ? 1.0f : gt.tautomerWeight[ti];
            conf *= static_cast<double>(std::min(wq, wt));
        }
    }
    return conf;
}

// ===========================================================================
// extractSubgraph -- extract a subgraph containing only the specified atom indices.
// Port of Java SearchEngine.extractSubgraph()
// ===========================================================================

inline MolGraph extractSubgraph(const MolGraph& mol, const std::vector<int>& atomIndices) {
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

    // Re-perceive ring membership and aromaticity on the actual subgraph topology.
    // This prevents aromatic flag leakage from partial ring extraction (e.g., 4 atoms
    // of a 6-membered ring are NOT aromatic as a standalone fragment).
    sub.ensureRingCounts();
    for (int i = 0; i < sub.n; ++i) {
        sub.ring[i] = (!sub.ringCount.empty() && sub.ringCount[i] > 0) ? 1 : 0;
        // Simple aromaticity: atom is aromatic only if it is in a ring AND was
        // aromatic in the parent AND all its ring neighbors in the subgraph are present
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

    // Kekulisation fallback: for atoms that were aromatic in the parent but are
    // NOT aromatic in the subgraph (broken ring, e.g. partial naphthalene
    // Kekule fallback: atoms that lost aromaticity in the subgraph (broken ring)
    // get their bonds downgraded to single (order 1) and aromatic flag cleared.
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
// murckoScaffold -- extract Murcko scaffold (ring systems + linkers).
// Port of Java SearchEngine.murckoScaffold()
// Strips all atoms not in rings and not on shortest paths between rings.
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

    // Collect kept atom indices
    std::vector<int> keptIndices;
    for (int i = 0; i < mol.n; ++i) {
        if (keep[i]) keptIndices.push_back(i);
    }

    if (static_cast<int>(keptIndices.size()) == mol.n) return mol;

    return extractSubgraph(mol, keptIndices);
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

#endif // SMSD_MOL_GRAPH_HPP

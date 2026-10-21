/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Configurable fingerprints: one specification shared by the circular (ECFP and FCFP),
 * path, topological torsion and substructure-key families, with the same bits in
 * C++, Java and Python.
 *
 * The classic configuration of each family, FingerprintOptions::ecfp(), fcfp(), path()
 * and torsion(), reproduces the bits of the fixed functions in batch.hpp exactly; every
 * further switch extends the hashed terms without moving the classic ones.
 */
#pragma once
#ifndef SMSD_FINGERPRINT_OPTIONS_HPP
#define SMSD_FINGERPRINT_OPTIONS_HPP

#include "smsd/batch.hpp"
#include "smsd/cip.hpp"
#include "smsd/stereo.hpp"
#include "smsd/substructure.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {
namespace batch {

/// The hashed families and the substructure-key family.
enum class FingerprintFamily { ECFP, FCFP, PATH, TORSION, KEYS };

/// How the connected components of a molecule are fingerprinted.
enum class FragmentMode {
    WHOLE,          ///< The whole molecule, as the fixed functions do.
    LARGEST,        ///< The component with the most heavy atoms only.
    PER_FRAGMENT    ///< Every atom term carries the class of its component, so bits are component-specific.
};

/// Whether the molecule is reduced before it is fingerprinted.
enum class ScaffoldMode {
    NONE,     ///< The molecule as given.
    MURCKO    ///< Its Bemis-Murcko scaffold, as murckoScaffold keeps it: ring atoms, the linkers between rings and the atoms double or triple bonded to either, capped with hydrogen.
};

/// What a substructure key counts.
enum class KeyKind {
    SMARTS,          ///< Distinct atom sets matched by a SMARTS pattern.
    ISOTOPE,         ///< Atoms with a mass number.
    AROMATIC_RINGS,  ///< Rings of the smallest set of smallest rings whose bonds are all aromatic.
    FRAGMENTS        ///< Connected components.
};

/// One substructure key: its bit is set when the molecule holds at least `minimum` of what it
/// counts, which for a SMARTS pattern is its distinct matched atom sets.
struct FingerprintKey {
    KeyKind kind = KeyKind::SMARTS;
    std::string smarts;
    int minimum = 1;

    FingerprintKey() = default;
    FingerprintKey(std::string pattern, int minimumCount = 1)
        : smarts(std::move(pattern)), minimum(minimumCount) {}
    static FingerprintKey special(KeyKind what, int minimumCount = 1) {
        FingerprintKey key;
        key.kind = what;
        key.minimum = minimumCount;
        return key;
    }
};

/**
 * A fingerprint specification. The default-constructed value is the classic ECFP
 * configuration; use the factories for the classic configuration of another family, then
 * change the switches you need. The same options give the same bits in Java and Python.
 */
struct FingerprintOptions {
    FingerprintFamily family = FingerprintFamily::ECFP;
    /// Circular families: the largest radius; -1 iterates to the largest component extent.
    int radius = 2;
    /// Path family: the largest number of bonds in a path.
    int pathLength = 7;
    /// Bits of a bit fingerprint or bins of a count fingerprint; the key family uses one per key.
    int size = 2048;
    /// Count fingerprint rather than bits, where one entry point serves both (Python).
    bool counts = false;

    // Atom terms.
    bool element = true;         ///< Atomic number.
    bool degree = true;          ///< Number of explicit neighbours (heavy atoms of a prepared graph).
    bool valence = true;         ///< Bond order sum; the torsion family reads its pi-electron count here.
    bool formalCharge = true;
    bool isotope = true;         ///< Mass number.
    bool hydrogenCount = true;   ///< Attached hydrogens, implicit and explicit.
    bool tautomerClass = true;   ///< The tautomer class, when one is assigned.
    bool aromaticity = true;

    // Bond terms.
    bool bondOrder = true;       ///< Bond order; aromatic bonds have order 4.
    bool aromaticBond = false;
    bool ringBond = false;

    // Ring terms.
    bool ringMembership = true;
    bool ringSize = false;       ///< Smallest ring through the atom, from the SSSR, as a class up to ringSizeCap.
    int ringSizeCap = 8;         ///< Rings of this size or more share one class.

    // Stereo terms, from the library's own CIP assignment.
    bool cipAtoms = false;       ///< R/S (and r/s) of tetrahedral centres.
    bool cipBonds = false;       ///< E/Z of stereogenic double bonds.

    FragmentMode fragments = FragmentMode::WHOLE;
    ScaffoldMode scaffold = ScaffoldMode::NONE;

    /// The keys of the key family, in bit order.
    std::vector<FingerprintKey> keys;

    /// The classic configuration of a family.
    static FingerprintOptions of(FingerprintFamily which) {
        FingerprintOptions options;
        options.family = which;
        switch (which) {
        case FingerprintFamily::ECFP:
            break;
        case FingerprintFamily::FCFP:
            options.element = options.valence = options.formalCharge = options.isotope = false;
            options.hydrogenCount = options.tautomerClass = options.aromaticity = false;
            break;
        case FingerprintFamily::PATH:
            options.degree = options.valence = options.formalCharge = options.isotope = false;
            options.hydrogenCount = options.tautomerClass = options.aromaticity = false;
            options.ringMembership = false;
            options.bondOrder = false;
            break;
        case FingerprintFamily::TORSION:
            options.formalCharge = options.isotope = options.hydrogenCount = false;
            options.tautomerClass = options.aromaticity = false;
            options.bondOrder = false;
            break;
        case FingerprintFamily::KEYS:
            options.size = 0;
            break;
        }
        return options;
    }
    static FingerprintOptions ecfp(int radius = 2, int size = 2048) {
        auto options = of(FingerprintFamily::ECFP);
        options.radius = radius;
        options.size = size;
        return options;
    }
    static FingerprintOptions fcfp(int radius = 2, int size = 2048) {
        auto options = of(FingerprintFamily::FCFP);
        options.radius = radius;
        options.size = size;
        return options;
    }
    static FingerprintOptions path(int pathLength = 7, int size = 2048) {
        auto options = of(FingerprintFamily::PATH);
        options.pathLength = pathLength;
        options.size = size;
        return options;
    }
    static FingerprintOptions torsion(int size = 2048) {
        auto options = of(FingerprintFamily::TORSION);
        options.size = size;
        return options;
    }
    /// A key fingerprint over the given keys, one bit each.
    static FingerprintOptions withKeys(std::vector<FingerprintKey> definitions) {
        auto options = of(FingerprintFamily::KEYS);
        options.keys = std::move(definitions);
        options.size = static_cast<int>(options.keys.size());
        return options;
    }
    /// The 166 MACCS structural keys; key k is bit k - 1.
    static FingerprintOptions maccs();
    /// The MACCS key definitions, in key order.
    static const std::vector<FingerprintKey>& maccsKeys();
};

// ============================================================================
// MACCS keys
// ============================================================================

namespace fingerprint_detail {

struct MaccsDefinition { const char* smarts; int minimum; };

/// The public MACCS key definitions as SMARTS with a minimum count of distinct matches;
/// keys 1, 125 and 166 are the isotope, aromatic ring and fragment counts.
inline const std::array<MaccsDefinition, 166>& maccsDefinitions() {
    static const std::array<MaccsDefinition, 166> definitions = {{
        {"", 1},                                                    // 1 isotope
        {"[#104]", 1},                                              // 2
        {"[#32,#33,#34,#50,#51,#52,#82,#83,#84]", 1},               // 3
        {"[Ac,Th,Pa,U,Np,Pu,Am,Cm,Bk,Cf,Es,Fm,Md,No,Lr]", 1},       // 4
        {"[Sc,Ti,Y,Zr,Hf]", 1},                                     // 5
        {"[La,Ce,Pr,Nd,Pm,Sm,Eu,Gd,Tb,Dy,Ho,Er,Tm,Yb,Lu]", 1},      // 6
        {"[V,Cr,Mn,Nb,Mo,Tc,Ta,W,Re]", 1},                          // 7
        {"[!#6;!#1]1~*~*~*~1", 1},                                  // 8
        {"[Fe,Co,Ni,Ru,Rh,Pd,Os,Ir,Pt]", 1},                        // 9
        {"[Be,Mg,Ca,Sr,Ba,Ra]", 1},                                 // 10
        {"*1~*~*~*~1", 1},                                          // 11
        {"[Cu,Zn,Ag,Cd,Au,Hg]", 1},                                 // 12
        {"[#8]~[#7](~[#6])~[#6]", 1},                               // 13
        {"[#16]-[#16]", 1},                                         // 14
        {"[#8]~[#6](~[#8])~[#8]", 1},                               // 15
        {"[!#6;!#1]1~*~*~1", 1},                                    // 16
        {"[#6]#[#6]", 1},                                           // 17
        {"[#5,#13,#31,#49,#81]", 1},                                // 18
        {"*1~*~*~*~*~*~*~1", 1},                                    // 19
        {"[#14]", 1},                                               // 20
        {"[#6]=[#6](~[!#6;!#1])~[!#6;!#1]", 1},                     // 21
        {"*1~*~*~1", 1},                                            // 22
        {"[#7]~[#6](~[#8])~[#8]", 1},                               // 23
        {"[#7]-[#8]", 1},                                           // 24
        {"[#7]~[#6](~[#7])~[#7]", 1},                               // 25
        {"[#6]=;@[#6](@*)@*", 1},                                   // 26
        {"[I]", 1},                                                 // 27
        {"[!#6;!#1]~[CH2]~[!#6;!#1]", 1},                           // 28
        {"[#15]", 1},                                               // 29
        {"[#6]~[!#6;!#1](~[#6])(~[#6])~*", 1},                      // 30
        {"[!#6;!#1]~[F,Cl,Br,I]", 1},                               // 31
        {"[#6]~[#16]~[#7]", 1},                                     // 32
        {"[#7]~[#16]", 1},                                          // 33
        {"[CH2]=*", 1},                                             // 34
        {"[Li,Na,K,Rb,Cs,Fr]", 1},                                  // 35
        {"[#16R]", 1},                                              // 36
        {"[#7]~[#6](~[#8])~[#7]", 1},                               // 37
        {"[#7]~[#6](~[#6])~[#7]", 1},                               // 38
        {"[#8]~[#16](~[#8])~[#8]", 1},                              // 39
        {"[#16]-[#8]", 1},                                          // 40
        {"[#6]#[#7]", 1},                                           // 41
        {"F", 1},                                                   // 42
        {"[!#6;!#1;!H0]~*~[!#6;!#1;!H0]", 1},                       // 43
        {"[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]", 1}, // 44
        {"[#6]=[#6]~[#7]", 1},                                      // 45
        {"Br", 1},                                                  // 46
        {"[#16]~*~[#7]", 1},                                        // 47
        {"[#8]~[!#6;!#1](~[#8])(~[#8])", 1},                        // 48
        {"[!+0]", 1},                                               // 49
        {"[#6]=[#6](~[#6])~[#6]", 1},                               // 50
        {"[#6]~[#16]~[#8]", 1},                                     // 51
        {"[#7]~[#7]", 1},                                           // 52
        {"[!#6;!#1;!H0]~*~*~*~[!#6;!#1;!H0]", 1},                   // 53
        {"[!#6;!#1;!H0]~*~*~[!#6;!#1;!H0]", 1},                     // 54
        {"[#8]~[#16]~[#8]", 1},                                     // 55
        {"[#8]~[#7](~[#8])~[#6]", 1},                               // 56
        {"[#8R]", 1},                                               // 57
        {"[!#6;!#1]~[#16]~[!#6;!#1]", 1},                           // 58
        {"[#16]!:*:*", 1},                                          // 59
        {"[#16]=[#8]", 1},                                          // 60
        {"*~[#16](~*)~*", 1},                                       // 61
        {"*@*!@*@*", 1},                                            // 62
        {"[#7]=[#8]", 1},                                           // 63
        {"*@*!@[#16]", 1},                                          // 64
        {"c:n", 1},                                                 // 65
        {"[#6]~[#6](~[#6])(~[#6])~*", 1},                           // 66
        {"[!#6;!#1]~[#16]", 1},                                     // 67
        {"[!#6;!#1;!H0]~[!#6;!#1;!H0]", 1},                         // 68
        {"[!#6;!#1]~[!#6;!#1;!H0]", 1},                             // 69
        {"[!#6;!#1]~[#7]~[!#6;!#1]", 1},                            // 70
        {"[#7]~[#8]", 1},                                           // 71
        {"[#8]~*~*~[#8]", 1},                                       // 72
        {"[#16]=*", 1},                                             // 73
        {"[CH3]~*~[CH3]", 1},                                       // 74
        {"*!@[#7]@*", 1},                                           // 75
        {"[#6]=[#6](~*)~*", 1},                                     // 76
        {"[#7]~*~[#7]", 1},                                         // 77
        {"[#6]=[#7]", 1},                                           // 78
        {"[#7]~*~*~[#7]", 1},                                       // 79
        {"[#7]~*~*~*~[#7]", 1},                                     // 80
        {"[#16]~*(~*)~*", 1},                                       // 81
        {"*~[CH2]~[!#6;!#1;!H0]", 1},                               // 82
        {"[!#6;!#1]1~*~*~*~*~1", 1},                                // 83
        {"[NH2]", 1},                                               // 84
        {"[#6]~[#7](~[#6])~[#6]", 1},                               // 85
        {"[C;H2,H3][!#6;!#1][C;H2,H3]", 1},                         // 86
        {"[F,Cl,Br,I]!@*@*", 1},                                    // 87
        {"[#16]", 1},                                               // 88
        {"[#8]~*~*~*~[#8]", 1},                                     // 89
        {"[$([!#6;!#1;!H0]~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[CH2;R]1),"
         "$([!#6;!#1;!H0]~[R]1@[R]@[CH2;R]1)]", 1},                 // 90
        {"[$([!#6;!#1;!H0]~*~*~*~[CH2]~*),$([!#6;!#1;!H0;R]1@[R]@[R]@[R]@[CH2;R]1),"
         "$([!#6;!#1;!H0]~[R]1@[R]@[R]@[CH2;R]1),$([!#6;!#1;!H0]~*~[R]1@[R]@[CH2;R]1)]", 1}, // 91
        {"[#8]~[#6](~[#7])~[#6]", 1},                               // 92
        {"[!#6;!#1]~[CH3]", 1},                                     // 93
        {"[!#6;!#1]~[#7]", 1},                                      // 94
        {"[#7]~*~*~[#8]", 1},                                       // 95
        {"*1~*~*~*~*~1", 1},                                        // 96
        {"[#7]~*~*~*~[#8]", 1},                                     // 97
        {"[!#6;!#1]1~*~*~*~*~*~1", 1},                              // 98
        {"[#6]=[#6]", 1},                                           // 99
        {"*~[CH2]~[#7]", 1},                                        // 100
        {"[$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
         "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
         "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
         "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1),"
         "$([R]@1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]1)]", 1}, // 101
        {"[!#6;!#1]~[#8]", 1},                                      // 102
        {"Cl", 1},                                                  // 103
        {"[!#6;!#1;!H0]~*~[CH2]~*", 1},                             // 104
        {"*@*(@*)@*", 1},                                           // 105
        {"[!#6;!#1]~*(~[!#6;!#1])~[!#6;!#1]", 1},                   // 106
        {"[F,Cl,Br,I]~*(~*)~*", 1},                                 // 107
        {"[CH3]~*~*~*~[CH2]~*", 1},                                 // 108
        {"*~[CH2]~[#8]", 1},                                        // 109
        {"[#7]~[#6]~[#8]", 1},                                      // 110
        {"[#7]~*~[CH2]~*", 1},                                      // 111
        {"*~*(~*)(~*)~*", 1},                                       // 112
        {"[#8]!:*:*", 1},                                           // 113
        {"[CH3]~[CH2]~*", 1},                                       // 114
        {"[CH3]~*~[CH2]~*", 1},                                     // 115
        {"[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]", 1},         // 116
        {"[#7]~*~[#8]", 1},                                         // 117
        {"[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", 2},             // 118
        {"[#7]=*", 1},                                              // 119
        {"[!#6;R]", 2},                                             // 120
        {"[#7;R]", 1},                                              // 121
        {"*~[#7](~*)~*", 1},                                        // 122
        {"[#8]~[#6]~[#8]", 1},                                      // 123
        {"[!#6;!#1]~[!#6;!#1]", 1},                                 // 124
        {"", 2},                                                    // 125 aromatic rings
        {"*!@[#8]!@*", 1},                                          // 126
        {"*@*!@[#8]", 2},                                           // 127
        {"[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2;R]@[R]@[R]@[R]@[CH2;R]1),"
         "$(*~[CH2]~[R]1@[R]@[R]@[CH2;R]1),$(*~[CH2]~*~[R]1@[R]@[CH2;R]1)]", 1}, // 128
        {"[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2;R]1),$(*~[CH2]~[R]1@[R]@[CH2;R]1)]", 1}, // 129
        {"[!#6;!#1]~[!#6;!#1]", 2},                                 // 130
        {"[!#6;!#1;!H0]", 2},                                       // 131
        {"[#8]~*~[CH2]~*", 1},                                      // 132
        {"*@*!@[#7]", 1},                                           // 133
        {"[F,Cl,Br,I]", 1},                                         // 134
        {"[#7]!:*:*", 1},                                           // 135
        {"[#8]=*", 2},                                              // 136
        {"[!C;!c;R]", 1},                                           // 137
        {"[!#6;!#1]~[CH2]~*", 2},                                   // 138
        {"[O;!H0]", 1},                                             // 139
        {"[#8]", 4},                                                // 140
        {"[CH3]", 3},                                               // 141
        {"[#7]", 2},                                                // 142
        {"*@*!@[#8]", 1},                                           // 143
        {"*!:*:*!:*", 1},                                           // 144
        {"*1~*~*~*~*~*~1", 2},                                      // 145
        {"[#8]", 3},                                                // 146
        {"[$(*~[CH2]~[CH2]~*),$([R]1@[CH2;R]@[CH2;R]1)]", 1},       // 147
        {"*~[!#6;!#1](~*)~*", 1},                                   // 148
        {"[C;H3,H4]", 2},                                           // 149
        {"*!@*@*!@*", 1},                                           // 150
        {"[#7;!H0]", 1},                                            // 151
        {"[#8]~[#6](~[#6])~[#6]", 1},                               // 152
        {"[!#6;!#1]~[CH2]~*", 1},                                   // 153
        {"[#6]=[#8]", 1},                                           // 154
        {"*!@[CH2]!@*", 1},                                         // 155
        {"[#7]~*(~*)~*", 1},                                        // 156
        {"[#6]-[#8]", 1},                                           // 157
        {"[#6]-[#7]", 1},                                           // 158
        {"[#8]", 2},                                                // 159
        {"[C;H3,H4]", 1},                                           // 160
        {"[#7]", 1},                                                // 161
        {"a", 1},                                                   // 162
        {"*1~*~*~*~*~*~1", 1},                                      // 163
        {"[#8]", 1},                                                // 164
        {"[R]", 1},                                                 // 165
        {"", 2},                                                    // 166 fragments
    }};
    return definitions;
}

} // namespace fingerprint_detail

inline const std::vector<FingerprintKey>& FingerprintOptions::maccsKeys() {
    static const std::vector<FingerprintKey> keys = [] {
        std::vector<FingerprintKey> built;
        const auto& definitions = fingerprint_detail::maccsDefinitions();
        for (size_t index = 0; index < definitions.size(); ++index) {
            const auto& definition = definitions[index];
            if (index == 0) built.push_back(FingerprintKey::special(KeyKind::ISOTOPE, definition.minimum));
            else if (index == 124) built.push_back(FingerprintKey::special(KeyKind::AROMATIC_RINGS, definition.minimum));
            else if (index == 165) built.push_back(FingerprintKey::special(KeyKind::FRAGMENTS, definition.minimum));
            else built.emplace_back(definition.smarts, definition.minimum);
        }
        return built;
    }();
    return keys;
}

inline FingerprintOptions FingerprintOptions::maccs() {
    return withKeys(maccsKeys());
}

// ============================================================================
// The engine
// ============================================================================

namespace fingerprint_detail {

using detail::fnvMix;
using detail::maximumComponentExtent;
namespace path_fp_detail = detail::path_fp_detail;

/// Mixing tags of the stereo descriptors, so a descriptor never collides with a small term.
constexpr uint64_t CIP_ATOM_TAG = 0xC1A000ULL;
constexpr uint64_t CIP_BOND_TAG = 0xE2E000ULL;
/// Neighbour-term marks of the aromatic and ring bond terms of the circular families.
constexpr uint64_t AROMATIC_BOND_MARK = 0x9E3779B97F4A7C15ULL;
constexpr uint64_t RING_BOND_MARK = 0xC2B2AE3D27D4EB4FULL;
/// Embeddings examined per key; a key whose count is not settled within them is an error, never a guess.
constexpr size_t KEY_MATCH_LIMIT = 100000;
/// Time allowed to one key search, in milliseconds, so that it normally completes; a search that
/// does time out is an error rather than a fingerprint that depends on the machine.
constexpr int64_t KEY_SEARCH_MILLIS = 3600000;

inline std::shared_ptr<const Substructure> preparedKey(const std::string& pattern);

inline void validate(const FingerprintOptions& options) {
    if (options.family == FingerprintFamily::KEYS) {
        if (options.keys.empty()) throw std::invalid_argument("keys must name at least one key");
        for (const auto& key : options.keys) {
            if (key.minimum < 1) throw std::invalid_argument("a key minimum must be positive");
            if (key.kind == KeyKind::SMARTS && key.smarts.empty())
                throw std::invalid_argument("a SMARTS key must have a pattern");
            if (key.kind == KeyKind::SMARTS) preparedKey(key.smarts);   // a pattern that does not parse fails here
        }
    } else {
        if (options.size <= 0) throw std::invalid_argument("size must be positive");
    }
    if (options.radius < -1) throw std::invalid_argument("radius must be >= -1 (use -1 for unlimited)");
    if (options.pathLength <= 0 || options.pathLength == std::numeric_limits<int>::max())
        throw std::invalid_argument("pathLength must be positive and representable");
    if (options.ringSizeCap < 3) throw std::invalid_argument("ringSizeCap must be at least 3");
}

/// Connected components in graph order; returns the component of each atom.
inline std::vector<int> components(const MolGraph& graph, int& count) {
    std::vector<int> component(graph.n, -1), queue;
    count = 0;
    for (int start = 0; start < graph.n; ++start) {
        if (component[start] >= 0) continue;
        queue.assign(1, start);
        component[start] = count;
        for (size_t head = 0; head < queue.size(); ++head)
            for (int neighbor : graph.neighbors[queue[head]])
                if (component[neighbor] < 0) { component[neighbor] = count; queue.push_back(neighbor); }
        ++count;
    }
    return component;
}

/// Whether the two orders of the same four ligands differ by an odd permutation.
inline bool oddPermutation(const std::array<int, 4>& actual, const std::array<int, 4>& derived) {
    std::array<int, 4> position{};
    for (int i = 0; i < 4; ++i) {
        position[i] = -1;
        for (int j = 0; j < 4; ++j) if (derived[j] == actual[i]) { position[i] = j; break; }
        if (position[i] < 0) return false;
    }
    int inversions = 0;
    for (int i = 0; i < 4; ++i) for (int j = i + 1; j < 4; ++j) if (position[i] > position[j]) ++inversions;
    return (inversions & 1) == 1;
}

/**
 * The tetrahedral tag of atom `atom` as the reduced graph reads it, in which it is atom `k`
 * with the kept neighbours `neighbors`. A removed neighbour becomes the implicit hydrogen that
 * caps the atom, so the ligand order the tag winds is re-read against the order the reduced
 * graph derives, and the tag flips when the two differ by an odd permutation. An atom that is
 * no longer a tetrahedral centre has no tag.
 */
inline int reducedTag(const MolGraph& graph, int atom, int k, const std::vector<int>& newIndex,
                      const std::vector<int>& neighbors) {
    const int tag = graph.tetraChirality[atom];
    std::array<int, 4> ligands{};
    if (tag == 0 || !stereo_detail::tetrahedralLigands(graph, atom, ligands)) return 0;
    std::array<int, 4> actual{};
    int implicitCount = 0;
    for (int slot = 0; slot < 4; ++slot) {
        actual[slot] = ligands[slot] < 0 || newIndex[ligands[slot]] < 0 ? -1 : newIndex[ligands[slot]];
        if (actual[slot] < 0) ++implicitCount;
    }
    const int degree = static_cast<int>(neighbors.size());
    if (!((degree == 4 && implicitCount == 0) || (degree == 3 && implicitCount == 1))) return 0;
    std::array<int, 4> derived{};
    const int implicit = degree == 4 ? -1 : (neighbors[0] < k ? 1 : 0);
    for (int slot = 0, position = 0; slot < 4; ++slot) derived[slot] = slot == implicit ? -1 : neighbors[position++];
    return oddPermutation(actual, derived) ? 3 - tag : tag;
}

/**
 * The configuration of the double bond `atom`-`neighbor` as the reduced graph reads it, or 0
 * when it is dropped. The stored relation refers to the first other neighbour at each end; when
 * that neighbour is removed the next kept one takes its place and the relation inverts, and
 * when none is left the bond is no longer stereogenic.
 */
inline int reducedConfiguration(const MolGraph& graph, int atom, int neighbor, const std::vector<int>& newIndex) {
    int configuration = graph.dbStereo(atom, neighbor);
    if (configuration == 0) return 0;
    for (const auto& end : {std::pair<int, int>{atom, neighbor}, std::pair<int, int>{neighbor, atom}}) {
        const int reference = stereo_detail::doubleBondReference(graph, end.first, end.second);
        if (reference < 0) return 0;
        if (newIndex[reference] >= 0) continue;
        int next = -1;
        for (int other : graph.neighbors[end.first])
            if (other != end.second && newIndex[other] >= 0) { next = other; break; }
        if (next < 0) return 0;
        configuration = 3 - configuration;
    }
    return configuration;
}

/**
 * The induced graph of `atoms`, in index order. Each kept atom receives `extraHydrogens[atom]`
 * hydrogens on top of its implicit count. Tetrahedral tags and double-bond configurations are
 * carried over, re-read against the reduced graph's neighbour order, see reducedTag() and
 * reducedConfiguration().
 */
inline MolGraph inducedGraph(const MolGraph& graph, const std::vector<int>& atoms,
                             const std::vector<int>& extraHydrogens) {
    const int count = static_cast<int>(atoms.size());
    std::vector<int> newIndex(graph.n, -1);
    for (int k = 0; k < count; ++k) newIndex[atoms[k]] = k;
    std::vector<int> atomicNum(count), formalCharge(count), massNumber(count), hydrogenCount(count);
    std::vector<int> atomClass(count, 0), tags;
    std::vector<uint8_t> ring(count), aromatic(count);
    std::vector<std::vector<int>> neighbors(count), orders(count);
    std::vector<std::vector<bool>> ringBonds(count), aromaticBonds(count);
    std::vector<std::tuple<int, int, int>> stereoBonds;
    const bool hasTags = !graph.tetraChirality.empty();
    if (hasTags) tags.assign(count, 0);
    for (int k = 0; k < count; ++k) {
        const int atom = atoms[k];
        atomicNum[k] = graph.atomicNum[atom];
        formalCharge[k] = graph.formalCharge[atom];
        massNumber[k] = atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0;
        hydrogenCount[k] = graph.implicitHydrogenCount(atom) + extraHydrogens[atom];
        if (atom < static_cast<int>(graph.atomClass.size())) atomClass[k] = graph.atomClass[atom];
        ring[k] = graph.ring[atom];
        aromatic[k] = graph.aromatic[atom];
        for (int neighbor : graph.neighbors[atom]) {
            if (newIndex[neighbor] < 0) continue;
            neighbors[k].push_back(newIndex[neighbor]);
            orders[k].push_back(graph.bondOrder(atom, neighbor));
            ringBonds[k].push_back(graph.bondInRing(atom, neighbor));
            aromaticBonds[k].push_back(graph.bondAromatic(atom, neighbor));
        }
        if (hasTags) tags[k] = reducedTag(graph, atom, k, newIndex, neighbors[k]);
        if (graph.hasDbStereo)
            for (int neighbor : graph.neighbors[atom]) {
                if (neighbor <= atom || newIndex[neighbor] < 0 || graph.bondOrder(atom, neighbor) != 2) continue;
                const int configuration = reducedConfiguration(graph, atom, neighbor, newIndex);
                if (configuration != 0) stereoBonds.emplace_back(k, newIndex[neighbor], configuration);
            }
    }
    MolGraph::Builder builder;
    builder.atomCount(count).atomicNumbers(std::move(atomicNum)).formalCharges(std::move(formalCharge))
        .massNumbers(std::move(massNumber)).hydrogenCounts(std::move(hydrogenCount))
        .atomClasses(std::move(atomClass)).ringFlags(std::move(ring)).aromaticFlags(std::move(aromatic))
        .setNeighbors(std::move(neighbors)).setBondOrders(std::move(orders))
        .bondRingFlags(std::move(ringBonds)).bondAromaticFlags(std::move(aromaticBonds));
    if (hasTags) builder.tetrahedralChirality(std::move(tags));
    if (!stereoBonds.empty()) builder.doubleBondStereoSparse(std::move(stereoBonds));
    return builder.build(false);
}

/**
 * The atoms of the Bemis-Murcko framework: ring atoms, the atoms on paths between rings, and
 * the atoms joined to those by a double or triple bond, which the framework keeps as the usual
 * convention does (an exocyclic carbonyl stays; a side chain is cut at its single bond). Chain
 * atoms are pruned from the leaves inwards in linear time; a chain atom survives exactly when
 * two directions from it lead to ring atoms, that is, when it lies on a path between rings.
 */
inline std::vector<int> frameworkAtoms(const MolGraph& graph) {
    std::vector<int> kept;
    bool hasRing = false;
    for (int atom = 0; atom < graph.n && !hasRing; ++atom) hasRing = graph.ring[atom] != 0;
    if (!hasRing) return kept;
    std::vector<int> remaining(graph.n), queue;
    std::vector<uint8_t> removed(graph.n, 0);
    for (int atom = 0; atom < graph.n; ++atom) {
        remaining[atom] = static_cast<int>(graph.neighbors[atom].size());
        if (!graph.ring[atom] && remaining[atom] <= 1) queue.push_back(atom);
    }
    for (size_t head = 0; head < queue.size(); ++head) {
        const int atom = queue[head];
        if (removed[atom]) continue;
        removed[atom] = 1;
        for (int neighbor : graph.neighbors[atom]) {
            if (removed[neighbor]) continue;
            if (--remaining[neighbor] <= 1 && !graph.ring[neighbor]) queue.push_back(neighbor);
        }
    }
    std::vector<uint8_t> keep(graph.n, 0);
    for (int atom = 0; atom < graph.n; ++atom) keep[atom] = removed[atom] ? 0 : 1;
    for (int atom = 0; atom < graph.n; ++atom) {
        if (removed[atom]) continue;
        for (int neighbor : graph.neighbors[atom]) {
            const int order = graph.bondOrder(atom, neighbor);
            if (removed[neighbor] && (order == 2 || order == 3)) keep[neighbor] = 1;
        }
    }
    for (int atom = 0; atom < graph.n; ++atom) if (keep[atom]) kept.push_back(atom);
    return kept;
}

/// The framework of a molecule with a ring, capped with hydrogen where side chains were; a
/// molecule without a ring, or one that is its own framework, is returned unchanged.
inline MolGraph frameworkGraph(const MolGraph& graph) {
    const auto kept = frameworkAtoms(graph);
    if (kept.empty() || static_cast<int>(kept.size()) == graph.n) return graph;
    std::vector<uint8_t> keep(graph.n, 0);
    for (int atom : kept) keep[atom] = 1;
    std::vector<int> extraHydrogens(graph.n, 0);
    for (int atom : kept)
        for (int neighbor : graph.neighbors[atom])
            if (!keep[neighbor]) {
                const int order = graph.bondOrder(atom, neighbor);
                extraHydrogens[atom] += order >= 1 && order <= 3 ? order : 1;
            }
    return inducedGraph(graph, kept, extraHydrogens);
}

/// A renumbering-invariant signature of an atom: the classic ECFP atom invariant without its tautomer class.
inline uint64_t atomSignature(const MolGraph& graph, int atom) {
    uint64_t h = FNV1A_SEED;
    h = fnvMix(h, static_cast<uint64_t>(graph.atomicNum[atom]));
    h = fnvMix(h, static_cast<uint64_t>(graph.degree[atom]));
    int bondOrderSum = 0;
    for (int neighbor : graph.neighbors[atom]) bondOrderSum += graph.bondOrder(atom, neighbor);
    h = fnvMix(h, static_cast<uint64_t>(bondOrderSum));
    h = fnvMix(h, static_cast<uint64_t>(atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0));
    h = fnvMix(h, static_cast<uint64_t>(int64_t(graph.formalCharge[atom]) + 4));
    h = fnvMix(h, static_cast<uint64_t>(graph.attachedHydrogenCount(atom)));
    h = fnvMix(h, graph.ring[atom] ? 1 : 0);
    h = fnvMix(h, graph.aromatic[atom] ? 1 : 0);
    return h;
}

/// A renumbering-invariant signature of a set of atoms: its atom and bond counts and sorted atom signatures.
inline uint64_t componentSignature(const MolGraph& graph, const std::vector<int>& atoms) {
    std::vector<uint64_t> terms;
    terms.reserve(atoms.size());
    int ends = 0;
    for (int atom : atoms) {
        terms.push_back(atomSignature(graph, atom));
        ends += static_cast<int>(graph.neighbors[atom].size());
    }
    std::sort(terms.begin(), terms.end());
    uint64_t s = FNV1A_SEED;
    s = fnvMix(s, static_cast<uint64_t>(atoms.size()));
    s = fnvMix(s, static_cast<uint64_t>(ends / 2));
    for (uint64_t term : terms) s = fnvMix(s, term);
    return s;
}

/**
 * The component with the most heavy atoms, then the most atoms; among components alike in both,
 * the one with the greater signature, so that the choice does not depend on the atom order, and
 * among identical components the first in atom order.
 */
inline MolGraph largestComponent(const MolGraph& graph) {
    int count = 0;
    const auto component = components(graph, count);
    if (count <= 1) return graph;
    std::vector<int> heavy(count, 0), atoms(count, 0);
    std::vector<std::vector<int>> members(count);
    for (int atom = 0; atom < graph.n; ++atom) {
        ++atoms[component[atom]];
        members[component[atom]].push_back(atom);
        if (graph.atomicNum[atom] != 1) ++heavy[component[atom]];
    }
    std::vector<uint64_t> signature(count, 0);
    int best = 0;
    for (int c = 1; c < count; ++c) {
        if (heavy[c] != heavy[best]) { if (heavy[c] > heavy[best]) best = c; continue; }
        if (atoms[c] != atoms[best]) { if (atoms[c] > atoms[best]) best = c; continue; }
        if (signature[best] == 0) signature[best] = componentSignature(graph, members[best]);
        signature[c] = componentSignature(graph, members[c]);
        if (signature[c] > signature[best]) best = c;
    }
    std::vector<int> kept;
    for (int atom = 0; atom < graph.n; ++atom) if (component[atom] == best) kept.push_back(atom);
    return inducedGraph(graph, kept, std::vector<int>(graph.n, 0));
}

/// The graph a fingerprint reads: the molecule after its scaffold and fragment handling.
struct Subject {
    const MolGraph* graph;
    std::unique_ptr<MolGraph> owned;
    static Subject of(const MolGraph& molecule, const FingerprintOptions& options) {
        Subject subject{&molecule, nullptr};
        if (options.scaffold == ScaffoldMode::MURCKO) {
            subject.owned = std::make_unique<MolGraph>(frameworkGraph(*subject.graph));
            subject.graph = subject.owned.get();
        }
        if (options.fragments == FragmentMode::LARGEST) {
            auto largest = std::make_unique<MolGraph>(largestComponent(*subject.graph));
            subject.owned = std::move(largest);
            subject.graph = subject.owned.get();
        }
        return subject;
    }
};

/// The smallest SSSR ring through each atom, as a class: 0 for an acyclic atom, else the size up to the cap.
inline std::vector<int> ringSizeClasses(const MolGraph& graph, int cap) {
    std::vector<int> classes(graph.n, 0);
    for (const auto& ring : ::smsd::detail_rf::completeBasis(graph)) {
        const int size = static_cast<int>(ring.size());
        for (int atom : ring)
            if (classes[atom] == 0 || size < classes[atom]) classes[atom] = size;
    }
    for (int& value : classes) value = std::min(value, cap);
    return classes;
}

/// The CIP descriptors as codes: R 1, S 2, r 3, s 4 by atom; E 1, Z 2 by double bond, sorted at each end.
struct StereoCodes {
    std::vector<int> atom;
    std::vector<std::vector<int>> bondsAt;
    std::map<std::pair<int, int>, int> bond;
};

inline StereoCodes stereoCodes(const MolGraph& graph, bool atoms, bool bonds) {
    StereoCodes codes;
    codes.atom.assign(graph.n, 0);
    codes.bondsAt.assign(graph.n, {});
    if (!atoms && !bonds) return codes;
    const auto descriptors = cip::assignAll(graph);
    if (atoms)
        for (int atom = 0; atom < graph.n; ++atom) {
            switch (descriptors.rsLabels[atom]) {
            case cip::RSLabel::R: codes.atom[atom] = 1; break;
            case cip::RSLabel::S: codes.atom[atom] = 2; break;
            case cip::RSLabel::r: codes.atom[atom] = 3; break;
            case cip::RSLabel::s: codes.atom[atom] = 4; break;
            default: break;
            }
        }
    if (bonds) {
        // A double bond in a ring of fewer than eight atoms has its geometry fixed by the ring
        // and carries no configuration, whatever marks the input wrote on it.
        std::set<std::pair<int, int>> smallRingBonds;
        for (const auto& ring : ::smsd::detail_rf::completeBasis(graph)) {
            if (ring.size() >= 8) continue;
            for (size_t k = 0; k < ring.size(); ++k) {
                const int a = ring[k], b = ring[(k + 1) % ring.size()];
                smallRingBonds.insert({std::min(a, b), std::max(a, b)});
            }
        }
        for (const auto& entry : descriptors.ezBonds) {
            const int first = std::min(std::get<0>(entry), std::get<1>(entry));
            const int second = std::max(std::get<0>(entry), std::get<1>(entry));
            // Only true double bonds carry a configuration in every language.
            if (graph.bondOrder(first, second) != 2) continue;
            if (smallRingBonds.count({first, second})) continue;
            const int code = std::get<2>(entry) == cip::EZLabel::E ? 1 : 2;
            codes.bond[{first, second}] = code;
            codes.bondsAt[first].push_back(code);
            codes.bondsAt[second].push_back(code);
        }
        for (auto& list : codes.bondsAt) std::sort(list.begin(), list.end());
    }
    return codes;
}

/// Everything the atom and bond terms of a hashed family read.
struct Terms {
    const MolGraph& graph;
    const FingerprintOptions& options;
    const std::vector<int>* features = nullptr;
    std::vector<int> ringSize;
    StereoCodes stereo;
    std::vector<int> component;
    int componentCount = 0;

    Terms(const MolGraph& g, const FingerprintOptions& o) : graph(g), options(o) {
        if (o.family == FingerprintFamily::FCFP) features = &g.getPharmacophoreFeatures();
        if (o.ringSize) ringSize = ringSizeClasses(g, o.ringSizeCap);
        stereo = stereoCodes(g, o.cipAtoms, o.cipBonds);
        if (o.fragments == FragmentMode::PER_FRAGMENT) component = components(g, componentCount);
    }

    int bondOrderSum(int atom) const {
        int sum = 0;
        for (int neighbor : graph.neighbors[atom]) sum += graph.bondOrder(atom, neighbor);
        return sum;
    }
    int piElectrons(int atom) const {
        int count = 0;
        for (int neighbor : graph.neighbors[atom]) {
            const int order = graph.bondOrder(atom, neighbor);
            if (order == 2) count += 1;
            else if (order == 3) count += 2;
        }
        if (graph.aromatic[atom]) ++count;
        return count;
    }
    int tautomer(int atom) const {
        return !graph.tautomerClass.empty() ? graph.tautomerClass[atom] : -1;
    }
    int bondCode(int from, int to) const {
        if (!options.cipBonds) return 0;
        const auto found = stereo.bond.find({std::min(from, to), std::max(from, to)});
        return found == stereo.bond.end() ? 0 : found->second;
    }

    /// The circular invariant of an atom before the fragment term.
    uint64_t circularBase(int atom) const {
        uint64_t h = FNV1A_SEED;
        if (features) h = fnvMix(h, static_cast<uint64_t>((*features)[atom]));
        if (options.element) h = fnvMix(h, static_cast<uint64_t>(graph.atomicNum[atom]));
        if (options.degree) h = fnvMix(h, static_cast<uint64_t>(graph.degree[atom]));
        if (options.valence) h = fnvMix(h, static_cast<uint64_t>(bondOrderSum(atom)));
        if (options.isotope) {
            const int mass = atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0;
            h = fnvMix(h, static_cast<uint64_t>(mass));
        }
        if (options.formalCharge) h = fnvMix(h, static_cast<uint64_t>(int64_t(graph.formalCharge[atom]) + 4));
        if (options.hydrogenCount) h = fnvMix(h, static_cast<uint64_t>(graph.attachedHydrogenCount(atom)));
        if (options.ringMembership) h = fnvMix(h, graph.ring[atom] ? 1 : 0);
        if (options.aromaticity) h = fnvMix(h, graph.aromatic[atom] ? 1 : 0);
        if (options.tautomerClass && tautomer(atom) >= 0)
            h = fnvMix(h, 0xCAFE00ULL + static_cast<uint64_t>(tautomer(atom)));
        if (options.ringSize) h = fnvMix(h, static_cast<uint64_t>(ringSize[atom]));
        return h;
    }

    /// The path atom term before the fragment term: the atomic number alone in the classic form.
    uint32_t pathBase(int atom) const {
        uint32_t t = options.element ? static_cast<uint32_t>(graph.atomicNum[atom]) : 0u;
        if (options.degree) t = t * 37u + static_cast<uint32_t>(graph.degree[atom]);
        if (options.valence) t = t * 37u + static_cast<uint32_t>(bondOrderSum(atom));
        if (options.isotope) {
            const int mass = atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0;
            t = t * 37u + static_cast<uint32_t>(mass);
        }
        if (options.formalCharge) t = t * 37u + static_cast<uint32_t>(graph.formalCharge[atom] + 4);
        if (options.hydrogenCount) t = t * 37u + static_cast<uint32_t>(graph.attachedHydrogenCount(atom));
        if (options.ringMembership) t = t * 37u + (graph.ring[atom] ? 1u : 0u);
        if (options.aromaticity) t = t * 37u + (graph.aromatic[atom] ? 1u : 0u);
        if (options.tautomerClass) t = t * 37u + static_cast<uint32_t>(tautomer(atom) + 1);
        if (options.ringSize) t = t * 37u + static_cast<uint32_t>(ringSize[atom]);
        if (options.cipAtoms) t = t * 37u + static_cast<uint32_t>(stereo.atom[atom]);
        return t;
    }
    bool pathHashesBonds() const {
        return options.bondOrder || options.aromaticBond || options.ringBond || options.cipBonds;
    }
    uint32_t pathBond(int from, int to) const {
        uint32_t b = 17u;
        if (options.bondOrder) b = b * 37u + static_cast<uint32_t>(graph.bondOrder(from, to));
        if (options.ringBond) b = b * 37u + (graph.bondInRing(from, to) ? 1u : 0u);
        if (options.aromaticBond) b = b * 37u + (graph.bondAromatic(from, to) ? 1u : 0u);
        if (options.cipBonds) b = b * 37u + static_cast<uint32_t>(bondCode(from, to));
        return b;
    }

    /// The torsion atom type before the fragment term: the classic packed type, extended by mixing.
    uint64_t torsionBase(int atom) const {
        const int z = options.element ? graph.atomicNum[atom] : 0;
        const int deg = options.degree ? graph.degree[atom] : 0;
        const int pi = options.valence ? piElectrons(atom) : 0;
        const int inRing = options.ringMembership && graph.ring[atom] ? 1 : 0;
        uint64_t t = static_cast<uint64_t>(static_cast<int>((z << 16) | (deg << 8) | (pi << 4) | inRing));
        if (options.isotope) {
            const int mass = atom < static_cast<int>(graph.massNumber.size()) ? graph.massNumber[atom] : 0;
            t = fnvMix(t, static_cast<uint64_t>(mass));
        }
        if (options.formalCharge) t = fnvMix(t, static_cast<uint64_t>(int64_t(graph.formalCharge[atom]) + 4));
        if (options.hydrogenCount) t = fnvMix(t, static_cast<uint64_t>(graph.attachedHydrogenCount(atom)));
        if (options.aromaticity) t = fnvMix(t, graph.aromatic[atom] ? 1 : 0);
        if (options.tautomerClass) t = fnvMix(t, static_cast<uint64_t>(tautomer(atom) + 1));
        if (options.ringSize) t = fnvMix(t, static_cast<uint64_t>(ringSize[atom]));
        if (options.cipAtoms && stereo.atom[atom] != 0) t = fnvMix(t, CIP_ATOM_TAG + static_cast<uint64_t>(stereo.atom[atom]));
        if (options.cipBonds)
            for (int code : stereo.bondsAt[atom]) t = fnvMix(t, CIP_BOND_TAG + static_cast<uint64_t>(code));
        return t;
    }

    /// The class of each component: a hash of its atom and bond counts and its sorted base terms.
    std::vector<uint64_t> fragmentTerms(const std::vector<uint64_t>& base) const {
        std::vector<std::vector<uint64_t>> perComponent(componentCount);
        std::vector<int> bonds(componentCount, 0);
        for (int atom = 0; atom < graph.n; ++atom) {
            perComponent[component[atom]].push_back(base[atom]);
            bonds[component[atom]] += static_cast<int>(graph.neighbors[atom].size());
        }
        std::vector<uint64_t> signature(componentCount);
        for (int c = 0; c < componentCount; ++c) {
            auto& terms = perComponent[c];
            std::sort(terms.begin(), terms.end());
            uint64_t s = FNV1A_SEED;
            s = fnvMix(s, static_cast<uint64_t>(terms.size()));
            s = fnvMix(s, static_cast<uint64_t>(bonds[c] / 2));
            for (uint64_t term : terms) s = fnvMix(s, term);
            signature[c] = s;
        }
        std::vector<uint64_t> result(graph.n);
        for (int atom = 0; atom < graph.n; ++atom) result[atom] = signature[component[atom]];
        return result;
    }
};

/// Hands every hash of a circular fingerprint to `sink`, see morganHashes().
template <class Sink>
inline void circularHashes(const Terms& terms, Sink&& sink) {
    const MolGraph& graph = terms.graph;
    const FingerprintOptions& options = terms.options;
    const int n = graph.n;
    std::vector<uint64_t> prevHash(n);
    for (int atom = 0; atom < n; ++atom) prevHash[atom] = terms.circularBase(atom);
    if (options.fragments == FragmentMode::PER_FRAGMENT) {
        const auto fragment = terms.fragmentTerms(prevHash);
        for (int atom = 0; atom < n; ++atom) prevHash[atom] = fnvMix(prevHash[atom], fragment[atom]);
    }
    for (int atom = 0; atom < n; ++atom) sink(prevHash[atom]);
    const int maxRadius = options.radius < 0 ? maximumComponentExtent(graph) : options.radius;
    std::vector<uint64_t> nextHash(n), nbHashes;
    for (int r = 1; r <= maxRadius; ++r) {
        for (int atom = 0; atom < n; ++atom) {
            const int deg = static_cast<int>(graph.neighbors[atom].size());
            if (deg == 0) { nextHash[atom] = prevHash[atom]; continue; }
            nbHashes.resize(deg);
            for (int k = 0; k < deg; ++k) {
                const int neighbor = graph.neighbors[atom][k];
                uint64_t term = prevHash[neighbor];
                if (options.bondOrder) term ^= static_cast<uint64_t>(graph.bondOrder(atom, neighbor)) * 31ULL;
                if (options.aromaticBond && graph.bondAromatic(atom, neighbor)) term ^= AROMATIC_BOND_MARK;
                if (options.ringBond && graph.bondInRing(atom, neighbor)) term ^= RING_BOND_MARK;
                nbHashes[k] = term;
            }
            std::sort(nbHashes.begin(), nbHashes.end());
            uint64_t h = FNV1A_SEED;
            h = fnvMix(h, static_cast<uint64_t>(r));
            h = fnvMix(h, prevHash[atom]);
            for (uint64_t nh : nbHashes) h = fnvMix(h, nh);
            if (r == 1) {
                // The descriptors join after the first iteration, so a centre changes the bits of
                // its own environments from radius 1 and of its neighbours' from radius 2.
                if (options.cipAtoms && terms.stereo.atom[atom] != 0)
                    h = fnvMix(h, CIP_ATOM_TAG + static_cast<uint64_t>(terms.stereo.atom[atom]));
                if (options.cipBonds)
                    for (int code : terms.stereo.bondsAt[atom]) h = fnvMix(h, CIP_BOND_TAG + static_cast<uint64_t>(code));
            }
            nextHash[atom] = h;
            sink(h);
        }
        std::swap(prevHash, nextHash);
    }
}

/// The path fingerprint of the terms, over heavy atoms, as the classic element path fingerprint.
inline std::vector<uint64_t> pathBits(const Terms& terms, int fpSize) {
    const MolGraph& graph = terms.graph;
    std::vector<uint32_t> atomTerm(graph.n);
    for (int atom = 0; atom < graph.n; ++atom) atomTerm[atom] = terms.pathBase(atom);
    if (terms.options.fragments == FragmentMode::PER_FRAGMENT) {
        std::vector<uint64_t> base(atomTerm.begin(), atomTerm.end());
        const auto fragment = terms.fragmentTerms(base);
        for (int atom = 0; atom < graph.n; ++atom)
            atomTerm[atom] = atomTerm[atom] * 37u + static_cast<uint32_t>(fragment[atom] ^ (fragment[atom] >> 32));
    }
    const auto keep = [&](int atom) { return graph.atomicNum[atom] != 1; };
    if (!terms.pathHashesBonds()) {
        const auto adjacency = path_fp_detail::pathAdjacency(graph, keep, [](int, int) { return 0u; });
        return path_fp_detail::pathFingerprint<path_fp_detail::AtomPathHash>(adjacency, atomTerm,
                                                                             terms.options.pathLength, fpSize);
    }
    const auto adjacency = path_fp_detail::pathAdjacency(graph, keep,
        [&](int from, int to) { return terms.pathBond(from, to); });
    return path_fp_detail::pathFingerprint<path_fp_detail::LabelledPathHash>(adjacency, atomTerm,
                                                                             terms.options.pathLength, fpSize);
}

/// Hands the canonical hash of every four-atom path to `sink`, see torsionHashes().
template <class Sink>
inline void torsionHashesOf(const Terms& terms, Sink&& sink) {
    const MolGraph& graph = terms.graph;
    const int n = graph.n;
    std::vector<uint64_t> atomType(n);
    for (int atom = 0; atom < n; ++atom) atomType[atom] = terms.torsionBase(atom);
    if (terms.options.fragments == FragmentMode::PER_FRAGMENT) {
        const auto fragment = terms.fragmentTerms(atomType);
        for (int atom = 0; atom < n; ++atom) atomType[atom] = fnvMix(atomType[atom], fragment[atom]);
    }
    if (n < 4) return;
    for (int a = 0; a < n; ++a)
        for (int b : graph.neighbors[a])
            for (int c : graph.neighbors[b]) {
                if (c == a) continue;
                for (int d : graph.neighbors[c]) {
                    if (d == b || d <= a) continue;
                    uint64_t fwd = FNV1A_SEED;
                    fwd = fnvMix(fwd, atomType[a]);
                    fwd = fnvMix(fwd, atomType[b]);
                    fwd = fnvMix(fwd, atomType[c]);
                    fwd = fnvMix(fwd, atomType[d]);
                    uint64_t rev = FNV1A_SEED;
                    rev = fnvMix(rev, atomType[d]);
                    rev = fnvMix(rev, atomType[c]);
                    rev = fnvMix(rev, atomType[b]);
                    rev = fnvMix(rev, atomType[a]);
                    sink(fwd <= rev ? fwd : rev);
                }
            }
}

// ---- Substructure keys ----

/// Prepared SMARTS searches by pattern, shared by every key fingerprint.
inline std::shared_ptr<const Substructure> preparedKey(const std::string& pattern) {
    static std::mutex guard;
    static std::unordered_map<std::string, std::shared_ptr<const Substructure>> prepared;
    {
        const std::lock_guard<std::mutex> lock(guard);
        const auto found = prepared.find(pattern);
        if (found != prepared.end()) return found->second;
    }
    SubstructureOptions options;
    options.timeoutMs = KEY_SEARCH_MILLIS;
    auto query = std::make_shared<const Substructure>(Substructure::smarts(pattern).withOptions(options));
    const std::lock_guard<std::mutex> lock(guard);
    if (prepared.size() >= 4096) prepared.clear();
    return prepared.emplace(pattern, query).first->second;
}

/// A key search must run to its end: a search cut short by time would make the bits depend on the machine.
inline void requireSettled(const SearchResult& result) {
    if (result.status == SearchStatus::TimedOut || result.status == SearchStatus::Cancelled)
        throw std::runtime_error("a substructure key search did not complete within its time budget");
}

/**
 * How many distinct atom sets a SMARTS pattern matches, or whether it matches at all when
 * `atMost` is one. Beyond the embedding limit the count is settled only when it has reached
 * `atMost`; otherwise the key is reported as an error rather than guessed.
 */
inline int keyMatches(const Substructure& query, const MolGraph& graph, int atMost) {
    if (atMost <= 1) {
        const SearchResult first = query.match(graph);
        requireSettled(first);
        return first.matches.empty() ? 0 : 1;
    }
    const SearchResult all = query.matchAll(graph, KEY_MATCH_LIMIT);
    requireSettled(all);
    std::set<std::vector<int>> distinct;
    for (const auto& match : all.matches) {
        std::vector<int> atoms;
        atoms.reserve(match.atomMapping.size());
        for (const auto& pair : match.atomMapping) atoms.push_back(pair.second);
        std::sort(atoms.begin(), atoms.end());
        distinct.insert(std::move(atoms));
    }
    const int count = static_cast<int>(distinct.size());
    if (all.status == SearchStatus::LimitReached && count < atMost)
        throw std::runtime_error("a substructure key matched more than 100000 embeddings before its count was settled");
    return count;
}

inline int aromaticRingCount(const MolGraph& graph) {
    int count = 0;
    for (const auto& ring : ::smsd::detail_rf::completeBasis(graph)) {
        bool aromatic = !ring.empty();
        for (size_t k = 0; k < ring.size() && aromatic; ++k)
            aromatic = graph.bondAromatic(ring[k], ring[(k + 1) % ring.size()]);
        if (aromatic) ++count;
    }
    return count;
}

/// The count of each key; `exact` asks for the full count rather than whether the minimum is met.
inline std::vector<int> keyCounts(const MolGraph& graph, const FingerprintOptions& options, bool exact) {
    std::vector<int> counts(options.keys.size(), 0);
    int componentCount = -1;
    for (size_t k = 0; k < options.keys.size(); ++k) {
        const auto& key = options.keys[k];
        switch (key.kind) {
        case KeyKind::SMARTS:
            counts[k] = keyMatches(*preparedKey(key.smarts), graph, exact ? std::numeric_limits<int>::max() : key.minimum);
            break;
        case KeyKind::ISOTOPE:
            for (int atom = 0; atom < graph.n; ++atom)
                if (atom < static_cast<int>(graph.massNumber.size()) && graph.massNumber[atom] != 0) ++counts[k];
            break;
        case KeyKind::AROMATIC_RINGS:
            counts[k] = aromaticRingCount(graph);
            break;
        case KeyKind::FRAGMENTS:
            if (componentCount < 0) components(graph, componentCount);
            counts[k] = componentCount;
            break;
        }
    }
    return counts;
}

} // namespace fingerprint_detail

// ============================================================================
// Public entry points
// ============================================================================

/// The bit fingerprint of a molecule under `options`, as words of 64 bits.
inline std::vector<uint64_t> computeFingerprint(const MolGraph& molecule, const FingerprintOptions& options) {
    using namespace fingerprint_detail;
    validate(options);
    const Subject subject = Subject::of(molecule, options);
    const MolGraph& graph = *subject.graph;
    if (options.family == FingerprintFamily::KEYS) {
        const int size = static_cast<int>(options.keys.size());
        std::vector<uint64_t> fp((size - 1) / 64 + 1, 0ULL);
        const auto counts = keyCounts(graph, options, false);
        for (int k = 0; k < size; ++k)
            if (counts[k] >= options.keys[k].minimum) detail::path_fp_detail::setBit(fp, k);
        return fp;
    }
    const int fpSize = options.size;
    std::vector<uint64_t> fp((fpSize - 1) / 64 + 1, 0ULL);
    if (graph.n == 0) return fp;
    const Terms terms(graph, options);
    const auto setBit = [&](uint64_t h) {
        detail::path_fp_detail::setBit(fp, static_cast<int>(h % static_cast<uint64_t>(fpSize)));
    };
    switch (options.family) {
    case FingerprintFamily::ECFP:
    case FingerprintFamily::FCFP:
        circularHashes(terms, setBit);
        break;
    case FingerprintFamily::PATH:
        fp = pathBits(terms, fpSize);
        break;
    case FingerprintFamily::TORSION:
        torsionHashesOf(terms, setBit);
        break;
    case FingerprintFamily::KEYS:
        break;
    }
    return fp;
}

/// The count fingerprint of a molecule under `options`: one bin per bit position, or per key.
inline std::vector<int> computeFingerprintCounts(const MolGraph& molecule, const FingerprintOptions& options) {
    using namespace fingerprint_detail;
    validate(options);
    const Subject subject = Subject::of(molecule, options);
    const MolGraph& graph = *subject.graph;
    if (options.family == FingerprintFamily::KEYS) return keyCounts(graph, options, true);
    const int fpSize = options.size;
    std::vector<int> counts(fpSize, 0);
    if (graph.n == 0) return counts;
    const Terms terms(graph, options);
    const auto count = [&](uint64_t h) { counts[static_cast<int>(h % static_cast<uint64_t>(fpSize))]++; };
    switch (options.family) {
    case FingerprintFamily::ECFP:
    case FingerprintFamily::FCFP:
        circularHashes(terms, count);
        break;
    case FingerprintFamily::PATH: {
        // Path bits are set once per path, so the count form counts the set positions.
        const auto bits = pathBits(terms, fpSize);
        for (int bit = 0; bit < fpSize; ++bit)
            if (bits[bit >> 6] & (uint64_t(1) << (bit & 63))) counts[bit] = 1;
        break;
    }
    case FingerprintFamily::TORSION:
        torsionHashesOf(terms, count);
        break;
    case FingerprintFamily::KEYS:
        break;
    }
    return counts;
}

/// Parses a SMILES string and computes its fingerprint under `options`.
inline std::vector<uint64_t> fingerprintFromSmiles(const std::string& smiles, const FingerprintOptions& options) {
    return computeFingerprint(smsd::parseSMILES(smiles), options);
}

/// The fingerprints of every molecule under `options`, in input order, on up to numThreads threads.
inline std::vector<std::vector<uint64_t>> batchFingerprint(const std::vector<MolGraph>& mols,
                                                           const FingerprintOptions& options,
                                                           int numThreads = 0) {
    fingerprint_detail::validate(options);
    const int N = static_cast<int>(mols.size());
    std::vector<std::vector<uint64_t>> results(N);
    if (N == 0) return results;
    const int nThreads = detail::workersFor(N, numThreads, 1);
    (void)nThreads;
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 8) num_threads(nThreads) if(nThreads > 1)
#endif
    for (int i = 0; i < N; ++i) results[i] = computeFingerprint(mols[i], options);
    return results;
}

} // namespace batch
} // namespace smsd

#endif // SMSD_FINGERPRINT_OPTIONS_HPP

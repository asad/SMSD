/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace smsd {
namespace detail {

// --------------------------------------------------------------------------
// Element table: atomic number -> symbol, and the reverse lookup built from it
// --------------------------------------------------------------------------
inline const char* atomicNumToSymbol(int z) {
    static const char* const syms[] = {
        "*","H","He",
        "Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar",
        "K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
        "Ag","Cd","In","Sn","Sb","Te","I","Xe",
        "Cs","Ba",
        "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra",
        "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
        "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    };
    if (z < 0 || z > 118) return "?";
    return syms[z];
}

/// Atomic number of an element symbol, including the aromatic lower-case SMILES
/// forms; -1 for an unknown symbol.
inline int elementToAtomicNum(const std::string& sym) {
    static const std::unordered_map<std::string, int> table = [] {
        std::unordered_map<std::string, int> symbols;
        symbols.reserve(160);
        for (int z = 0; z <= 118; ++z) symbols.emplace(atomicNumToSymbol(z), z);
        for (const auto& aromatic : {std::pair<const char*, int>{"b", 5}, {"c", 6}, {"n", 7}, {"o", 8},
                                     {"p", 15}, {"s", 16}, {"se", 34}, {"te", 52}})
            symbols.emplace(aromatic.first, aromatic.second);
        return symbols;
    }();
    auto it = table.find(sym);
    return it != table.end() ? it->second : -1;
}

// --------------------------------------------------------------------------
// Organic subset: atoms that can appear without brackets
// --------------------------------------------------------------------------
inline bool isOrganicSubset(int z) {
    return z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16 ||
           z == 9 || z == 17 || z == 35 || z == 53;
}

// --------------------------------------------------------------------------
// Default valences for implicit H computation: the one valence table shared by
// the SMILES, SMARTS, MOL, CIP and fingerprint code (Java mirrors it in
// MolGraph.DEFAULT_VALENCES).
// --------------------------------------------------------------------------

/// The default valences of an element in ascending order: a view over a static
/// table, so looking them up on a hot path allocates nothing.
struct DefaultValences {
    const int* values;
    int count;
    bool empty() const { return count == 0; }
    int size() const { return count; }
    int operator[](int index) const { return values[index]; }
    int front() const { return values[0]; }
    int back() const { return values[count - 1]; }
    const int* begin() const { return values; }
    const int* end() const { return values + count; }
};

inline DefaultValences defaultValences(int z) {
    static constexpr int one[] = {1}, two[] = {2}, three[] = {3}, four[] = {4};
    static constexpr int threeFive[] = {3, 5}, twoFourSix[] = {2, 4, 6}, oddToSeven[] = {1, 3, 5, 7};
    switch (z) {
        // Period 1
        case 1:  return {one, 1};             // H
        // Period 2
        case 5:  return {three, 1};           // B
        case 6:  return {four, 1};            // C
        case 7:  return {threeFive, 2};       // N
        case 8:  return {two, 1};             // O
        case 9:  return {one, 1};             // F
        // Period 3
        case 13: return {three, 1};           // Al
        case 14: return {four, 1};            // Si
        case 15: return {threeFive, 2};       // P
        case 16: return {twoFourSix, 3};      // S
        case 17: return {one, 1};             // Cl
        // Period 4
        case 33: return {threeFive, 2};       // As
        case 34: return {twoFourSix, 3};      // Se
        case 35: return {one, 1};             // Br
        // Period 5
        case 51: return {threeFive, 2};       // Sb
        case 52: return {twoFourSix, 3};      // Te
        case 53: return {oddToSeven, 4};      // I
        default: return {nullptr, 0};
    }
}

// --------------------------------------------------------------------------
// Valence target of a charged atom
//   Cations of group-15/16 elements expand into the next valence shell
//   (N+->4, P+->4, O+->3, S+->3, As+->4, Se+->3); group-13 anions (B-, Al-)
//   become isoelectronic with group 14 ([BH4-] has valence 4).  Every other
//   charged atom loses |charge| from its available valence.
// --------------------------------------------------------------------------
inline int64_t chargeAdjustedValence(int atomicNum, int valence, int charge) {
    const int64_t absoluteCharge = charge < 0 ? -int64_t(charge) : int64_t(charge);
    const bool expands = (charge > 0 && (atomicNum == 7  || atomicNum == 15 ||
                                         atomicNum == 8  || atomicNum == 16 ||
                                         atomicNum == 33 || atomicNum == 34))
                      || (charge < 0 && (atomicNum == 5 || atomicNum == 13));
    return expands ? int64_t(valence) + absoluteCharge : int64_t(valence) - absoluteCharge;
}

// --------------------------------------------------------------------------
// Compute implicit H count for an organic-subset atom
// Following OpenSMILES valence model:
//   effectiveBondOrder = sum(bond orders) [aromatic bonds count as 1]
//                        + 1 if atom is aromatic (pi contribution)
//   implicitH = smallestFittingValence - effectiveBondOrder
// where smallestFittingValence is the smallest default valence v such that
//   (v - charge) >= effectiveBondOrder   [charge adjustment per group]
// --------------------------------------------------------------------------
inline int computeImplicitH(int atomicNum, bool isAromatic, int totalBondOrder, int charge) {
    const DefaultValences valences = defaultValences(atomicNum);
    if (valences.empty()) return 0;

    // Aromatic lone-pair donors do not gain H from a higher valence state.
    if (isAromatic && charge == 0 && totalBondOrder >= valences.front())
        return 0;

    // Effective bond order: sum of explicit bond orders + aromatic pi
    int ebo = totalBondOrder;
    if (isAromatic) ebo += 1;

    // Per OpenSMILES / Daylight convention:
    //   For organic-subset atoms (always charge 0), implicitH = v - ebo
    //   for the smallest fitting default valence v.
    //
    //   For charged bracket atoms, MolGraph doesn't store the original
    //   explicit H count.  We use the heuristic:  target = v + absCharge
    //   for cations of group-15/16 elements (N+->4, S+->3, etc.) and
    //   target = v - absCharge otherwise.  This covers the vast majority
    //   of drug-like molecules.
    // Early out: if ebo exceeds the largest conceivable target, no valence can fit
    int maxVal = valences.back();
    const int64_t absoluteCharge = charge < 0 ? -int64_t(charge) : int64_t(charge);
    const int64_t maxTarget = maxVal + absoluteCharge;
    if (ebo > maxTarget) return 0;

    for (int v : valences) {
        const int64_t hc = chargeAdjustedValence(atomicNum, v, charge) - ebo;
        if (hc >= 0) return static_cast<int>(std::min(hc, int64_t(std::numeric_limits<int>::max())));
    }

    return 0;
}

} // namespace detail
} // namespace smsd

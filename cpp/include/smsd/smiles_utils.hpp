/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
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
// Element table: symbol -> atomic number
// --------------------------------------------------------------------------
inline int elementToAtomicNum(const std::string& sym) {
    // This covers elements up to Oganesson (118).  Only the ones that
    // are likely to appear in SMILES strings are listed; the rest can
    // be added trivially.
    static const std::unordered_map<std::string, int> table = {
        {"*",0},
        {"H",1},{"He",2},
        {"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
        {"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},
        {"K",19},{"Ca",20},
        {"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},
        {"Cu",29},{"Zn",30},{"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},
        {"Rb",37},{"Sr",38},
        {"Y",39},{"Zr",40},{"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},
        {"Ag",47},{"Cd",48},{"In",49},{"Sn",50},{"Sb",51},{"Te",52},{"I",53},{"Xe",54},
        {"Cs",55},{"Ba",56},
        {"La",57},{"Ce",58},{"Pr",59},{"Nd",60},{"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},
        {"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},{"Lu",71},
        {"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},
        {"Au",79},{"Hg",80},{"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},
        {"Fr",87},{"Ra",88},
        {"Ac",89},{"Th",90},{"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Am",95},{"Cm",96},
        {"Bk",97},{"Cf",98},{"Es",99},{"Fm",100},{"Md",101},{"No",102},{"Lr",103},
        {"Rf",104},{"Db",105},{"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},
        {"Rg",111},{"Cn",112},{"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118},
        // Aromatic lowercase forms
        {"b",5},{"c",6},{"n",7},{"o",8},{"p",15},{"s",16},{"se",34},{"te",52},
    };
    auto it = table.find(sym);
    return it != table.end() ? it->second : -1;
}

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

// --------------------------------------------------------------------------
// Organic subset: atoms that can appear without brackets
// --------------------------------------------------------------------------
inline bool isOrganicSubset(int z) {
    return z == 5 || z == 6 || z == 7 || z == 8 || z == 15 || z == 16 ||
           z == 9 || z == 17 || z == 35 || z == 53;
}

// --------------------------------------------------------------------------
// Default valences for implicit H computation (organic subset)
// --------------------------------------------------------------------------
inline std::vector<int> defaultValences(int z) {
    switch (z) {
        // Period 1
        case 1:  return {1};             // H
        // Period 2
        case 5:  return {3};             // B
        case 6:  return {4};             // C
        case 7:  return {3, 5};          // N
        case 8:  return {2};             // O
        case 9:  return {1};             // F
        // Period 3
        case 13: return {3};             // Al
        case 14: return {4};             // Si
        case 15: return {3, 5};          // P
        case 16: return {2, 4, 6};       // S
        case 17: return {1};             // Cl
        // Period 4
        case 33: return {3, 5};          // As
        case 34: return {2, 4, 6};       // Se
        case 35: return {1};             // Br
        // Period 5
        case 51: return {3, 5};          // Sb
        case 52: return {2, 4, 6};       // Te
        case 53: return {1, 3, 5, 7};    // I
        default: return {};
    }
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
    auto valences = defaultValences(atomicNum);
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
        int64_t target;
        // Cations of group-15/16 elements expand into the next valence shell
        // (N+->4, P+->4, O+->3, S+->3, As+->4, Se+->3).
        // All other cases (anions, metals, etc.) simply lose available valence.
        if (charge > 0 && (atomicNum == 7  || atomicNum == 15 ||
                           atomicNum == 8  || atomicNum == 16 ||
                           atomicNum == 33 || atomicNum == 34)) {
            target = int64_t(v) + charge;
        } else if (charge < 0 && (atomicNum == 5 || atomicNum == 13)) {
            // Group 13 anions (B⁻, Al⁻) gain valence capacity,
            // becoming isoelectronic with Group 14 (e.g., [BH4-] has valence 4)
            target = v + absoluteCharge;
        } else {
            target = v - absoluteCharge;
        }
        const int64_t hc = target - ebo;
        if (hc >= 0) return static_cast<int>(std::min(hc, int64_t(std::numeric_limits<int>::max())));
    }

    return 0;
}

} // namespace detail
} // namespace smsd

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only OpenSMILES parser and writer for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_SMILES_PARSER_HPP
#define SMSD_SMILES_PARSER_HPP

#include "mol_graph.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <deque>
#include <functional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// Parse options: strict (default) vs lenient mode
// ============================================================================
struct ParseOptions {
    bool lenient = false;  // false = throw on errors; true = best-effort recovery
};

// ============================================================================
// Internal implementation details
// ============================================================================
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
// Parsed atom struct (internal)
// --------------------------------------------------------------------------
struct ParsedAtom {
    int atomicNum     = 6;
    int charge        = 0;
    int isotope       = 0;
    int hcount        = -1;  // -1 = compute from valence; >=0 = explicit bracket hcount
    int chirality     = 0;   // 0=none, 1=@, 2=@@
    int atomClass     = 0;
    bool isAromatic   = false;
    bool isBracket    = false;
};

// --------------------------------------------------------------------------
// Parsed bond struct (internal)
// --------------------------------------------------------------------------
struct ParsedBond {
    int from = -1;
    int to   = -1;
    int order = 1;       // 1=single, 2=double, 3=triple, 4=aromatic(1.5)
    int stereo = 0;      // 0=none, 1=up(/), 2=down(\)
};

// --------------------------------------------------------------------------
// SMILES tokeniser / parser state
// --------------------------------------------------------------------------
class SmilesParser {
public:
    explicit SmilesParser(const std::string& smiles, bool lenient = false)
        : smi_(smiles), pos_(0), lenient_(lenient) {}

    void parse() {
        if (smi_.empty()) throw std::invalid_argument("Empty SMILES string");
        parseChain(-1, 0, 0);
        if (pos_ < smi_.size()) {
            if (!lenient_)
                throw std::invalid_argument(
                    std::string("Unexpected character at position ") +
                    std::to_string(pos_) + ": '" + smi_[pos_] + "'");
            // lenient: skip extra ')' and continue parsing remaining input
            while (pos_ < smi_.size()) {
                if (smi_[pos_] == ')') {
                    pos_++; // silently ignore extra closing paren
                    // Continue parsing any remaining atoms after the stray ')'
                    if (pos_ < smi_.size() && smi_[pos_] != ')')
                        parseChain(atoms_.empty() ? -1 : static_cast<int>(atoms_.size()) - 1, 0, 0);
                } else {
                    pos_++; // skip truly unrecognized characters
                }
            }
        }
        closeRings();
    }

    const std::vector<ParsedAtom>& atoms() const { return atoms_; }
    const std::vector<ParsedBond>& bonds() const { return bonds_; }

private:
    std::string smi_;
    size_t pos_;
    bool lenient_;
    std::vector<ParsedAtom> atoms_;
    std::vector<ParsedBond> bonds_;

    // Ring closure map: ring number -> (atom index, bond order, stereo)
    std::unordered_map<int, std::tuple<int, int, int>> ringOpenings_;

    void parseChain(int prevAtom, int pendingBondOrder, int pendingStereo) {
        while (pos_ < smi_.size()) {
            char c = smi_[pos_];

            // Dot disconnection
            if (c == '.') {
                pos_++;
                prevAtom = -1;
                pendingBondOrder = 0;
                pendingStereo = 0;
                continue;
            }

            // Branch open
            if (c == '(') {
                pos_++;
                int bo = 0, st = 0;
                parseBondSymbol(bo, st);
                parseChain(prevAtom, bo, st);
                if (pos_ >= smi_.size() || smi_[pos_] != ')') {
                    if (!lenient_)
                        throw std::invalid_argument("Unmatched opening parenthesis");
                    // lenient: silently close the open branch
                } else {
                    pos_++;
                }
                continue;
            }

            // Branch close
            if (c == ')') {
                return;
            }

            // Bond symbol (not followed by atom yet, just record)
            int bondOrd = 0, bondStereo = 0;
            if (parseBondSymbol(bondOrd, bondStereo)) {
                pendingBondOrder = bondOrd;
                pendingStereo = bondStereo;
            }

            // Try to parse an atom
            int atomIdx = -1;
            if (pos_ < smi_.size() && smi_[pos_] == '[') {
                atomIdx = parseBracketAtom();
            } else if (pos_ < smi_.size()) {
                atomIdx = parseOrganicAtom();
            }

            if (atomIdx < 0) {
                // No atom parsed, and we have no more bond symbols
                if (pendingBondOrder != 0 || pendingStereo != 0) {
                    if (!lenient_)
                        throw std::invalid_argument("Bond symbol not followed by atom");
                    // lenient: default to single bond, discard pending bond
                    pendingBondOrder = 0;
                    pendingStereo = 0;
                }
                return;
            }

            // Connect to previous atom
            if (prevAtom >= 0) {
                int bo = pendingBondOrder;
                int st = pendingStereo;
                if (bo == 0) {
                    // Implicit bond: aromatic if both atoms aromatic, else single
                    if (atoms_[prevAtom].isAromatic && atoms_[atomIdx].isAromatic)
                        bo = 4; // aromatic
                    else
                        bo = 1;
                }
                bonds_.push_back({prevAtom, atomIdx, bo, st});
            }
            pendingBondOrder = 0;
            pendingStereo = 0;

            // Ring closures
            parseRingClosures(atomIdx);

            prevAtom = atomIdx;
        }
    }

    bool parseBondSymbol(int& order, int& stereo) {
        if (pos_ >= smi_.size()) return false;
        char c = smi_[pos_];
        switch (c) {
            case '-': order = 1; stereo = 0; pos_++; return true;
            case '=': order = 2; stereo = 0; pos_++; return true;
            case '#': order = 3; stereo = 0; pos_++; return true;
            case ':': order = 4; stereo = 0; pos_++; return true;
            case '/': order = 0; stereo = 1; pos_++; return true;  // order determined later
            case '\\': order = 0; stereo = 2; pos_++; return true;
            default: return false;
        }
    }

    int parseBracketAtom() {
        if (pos_ >= smi_.size() || smi_[pos_] != '[') return -1;
        pos_++; // skip '['

        ParsedAtom a;
        a.isBracket = true;

        // Isotope (optional digits)
        int iso = 0;
        bool hasIso = false;
        while (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
            iso = iso * 10 + (smi_[pos_] - '0');
            pos_++;
            hasIso = true;
        }
        if (hasIso) a.isotope = iso;

        // Element symbol
        if (pos_ >= smi_.size()) {
            if (!lenient_)
                throw std::invalid_argument("Unexpected end in bracket atom");
            return -1; // lenient: skip incomplete bracket atom
        }

        std::string sym;
        char c = smi_[pos_];

        // Check for aromatic lowercase atoms inside brackets
        if (c >= 'a' && c <= 'z') {
            // Try two-char aromatic: se, te
            if (pos_ + 1 < smi_.size()) {
                std::string two{c, smi_[pos_ + 1]};
                if (two == "se" || two == "te") {
                    sym = two;
                    pos_ += 2;
                    a.isAromatic = true;
                } else {
                    sym = std::string(1, c);
                    pos_++;
                    a.isAromatic = true;
                }
            } else {
                sym = std::string(1, c);
                pos_++;
                a.isAromatic = true;
            }
        } else if (c >= 'A' && c <= 'Z') {
            if (c == 'R' && pos_ + 1 < smi_.size()
                && std::isdigit(static_cast<unsigned char>(smi_[pos_ + 1]))) {
                // Patent-style R-group placeholder, e.g. [R1] -> wildcard + class 1.
                sym = "*";
                pos_++;
                int cls = 0;
                while (pos_ < smi_.size() && std::isdigit(static_cast<unsigned char>(smi_[pos_]))) {
                    cls = cls * 10 + (smi_[pos_] - '0');
                    pos_++;
                }
                a.atomClass = cls;
            } else {
                sym = std::string(1, c);
                pos_++;
                // Second letter (lowercase)
                if (pos_ < smi_.size() && smi_[pos_] >= 'a' && smi_[pos_] <= 'z') {
                    sym += smi_[pos_];
                    pos_++;
                }
            }
        } else if (c == '*') {
            sym = "*";
            pos_++;
        } else {
            if (!lenient_)
                throw std::invalid_argument(
                    std::string("Invalid element in bracket atom: '") + c + "'");
            // lenient: skip to closing ']' and return -1
            while (pos_ < smi_.size() && smi_[pos_] != ']') pos_++;
            if (pos_ < smi_.size()) pos_++; // skip ']'
            return -1;
        }

        int z = elementToAtomicNum(sym);
        if (z < 0) {
            // Try uppercase first letter only
            if (sym.size() == 2) {
                std::string s1(1, sym[0]);
                int z1 = elementToAtomicNum(s1);
                if (z1 >= 0) {
                    z = z1;
                    pos_--; // put back the second char
                } else {
                    if (!lenient_)
                        throw std::invalid_argument("Unknown element: " + sym);
                    // lenient: skip to closing ']' and return -1
                    while (pos_ < smi_.size() && smi_[pos_] != ']') pos_++;
                    if (pos_ < smi_.size()) pos_++; // skip ']'
                    return -1;
                }
            } else {
                if (!lenient_)
                    throw std::invalid_argument("Unknown element: " + sym);
                // lenient: skip to closing ']' and return -1
                while (pos_ < smi_.size() && smi_[pos_] != ']') pos_++;
                if (pos_ < smi_.size()) pos_++; // skip ']'
                return -1;
            }
        }
        a.atomicNum = z;

        // Chirality (optional)
        if (pos_ < smi_.size() && smi_[pos_] == '@') {
            pos_++;
            if (pos_ < smi_.size() && smi_[pos_] == '@') {
                a.chirality = 2;
                pos_++;
            } else {
                a.chirality = 1;
            }
            // Skip known CIP chirality class prefixes: TH, AL, SP, TB, OH.
            // Only consume uppercase letters + trailing digits when we see one
            // of these exact two-letter prefixes. This avoids consuming 'H'
            // (hydrogen count) or other valid bracket-atom tokens.
            if (pos_ < smi_.size() && (pos_ + 1) < smi_.size()) {
                char c0 = smi_[pos_], c1 = smi_[pos_ + 1];
                if ((c0 == 'T' && c1 == 'H') || (c0 == 'A' && c1 == 'L') ||
                    (c0 == 'S' && c1 == 'P') || (c0 == 'T' && c1 == 'B') ||
                    (c0 == 'O' && c1 == 'H')) {
                    while (pos_ < smi_.size() && std::isalpha(smi_[pos_])) pos_++;
                    while (pos_ < smi_.size() && std::isdigit(smi_[pos_])) pos_++;
                }
            }
        }

        // Hydrogen count (optional)
        if (pos_ < smi_.size() && smi_[pos_] == 'H') {
            pos_++;
            int hc = 1;
            if (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                hc = 0;
                while (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                    hc = hc * 10 + (smi_[pos_] - '0');
                    pos_++;
                }
            }
            a.hcount = hc;
        } else {
            a.hcount = 0;  // bracket atom with no H specified means 0
        }

        // Charge (optional)
        if (pos_ < smi_.size() && (smi_[pos_] == '+' || smi_[pos_] == '-')) {
            char sign = smi_[pos_];
            pos_++;
            int chg = 1;
            if (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                chg = 0;
                while (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                    chg = chg * 10 + (smi_[pos_] - '0');
                    pos_++;
                }
            } else {
                // Count repeated signs: ++ = +2, +++ = +3
                while (pos_ < smi_.size() && smi_[pos_] == sign) {
                    chg++;
                    pos_++;
                }
            }
            a.charge = (sign == '+') ? chg : -chg;
        }

        // Atom class (optional :nn)
        if (pos_ < smi_.size() && smi_[pos_] == ':') {
            pos_++;
            int cls = 0;
            while (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                cls = cls * 10 + (smi_[pos_] - '0');
                pos_++;
            }
            a.atomClass = cls;
        }

        if (pos_ >= smi_.size() || smi_[pos_] != ']') {
            if (!lenient_)
                throw std::invalid_argument("Missing closing ']' in bracket atom");
            // lenient: assume bracket is closed
        } else {
            pos_++;
        }

        atoms_.push_back(a);
        return static_cast<int>(atoms_.size()) - 1;
    }

    int parseOrganicAtom() {
        if (pos_ >= smi_.size()) return -1;
        char c = smi_[pos_];

        // Aromatic atoms: b, c, n, o, p, s
        if (c == 'b' || c == 'c' || c == 'n' || c == 'o' || c == 'p' || c == 's') {
            ParsedAtom a;
            a.isAromatic = true;
            a.atomicNum = elementToAtomicNum(std::string(1, c));
            a.hcount = -1; // compute from valence
            pos_++;
            atoms_.push_back(a);
            return static_cast<int>(atoms_.size()) - 1;
        }

        // Two-letter organic: Cl, Br
        if (c == 'C' && pos_ + 1 < smi_.size() && smi_[pos_ + 1] == 'l') {
            ParsedAtom a;
            a.atomicNum = 17;
            a.hcount = -1;
            pos_ += 2;
            atoms_.push_back(a);
            return static_cast<int>(atoms_.size()) - 1;
        }
        if (c == 'B' && pos_ + 1 < smi_.size() && smi_[pos_ + 1] == 'r') {
            ParsedAtom a;
            a.atomicNum = 35;
            a.hcount = -1;
            pos_ += 2;
            atoms_.push_back(a);
            return static_cast<int>(atoms_.size()) - 1;
        }

        // Single-letter organic: B, C, N, O, P, S, F, I
        if (c == 'B' || c == 'C' || c == 'N' || c == 'O' ||
            c == 'P' || c == 'S' || c == 'F' || c == 'I') {
            ParsedAtom a;
            a.atomicNum = elementToAtomicNum(std::string(1, c));
            a.hcount = -1;
            pos_++;
            atoms_.push_back(a);
            return static_cast<int>(atoms_.size()) - 1;
        }

        return -1;
    }

    void parseRingClosures(int atomIdx) {
        while (pos_ < smi_.size()) {
            int ringBondOrder = 0, ringBondStereo = 0;
            // Optional bond symbol before ring digit
            size_t savedPos = pos_;
            parseBondSymbol(ringBondOrder, ringBondStereo);

            int ringNum = -1;
            if (pos_ < smi_.size() && smi_[pos_] == '%') {
                pos_++;
                if (pos_ + 1 < smi_.size() && std::isdigit(smi_[pos_]) &&
                    std::isdigit(smi_[pos_ + 1])) {
                    ringNum = (smi_[pos_] - '0') * 10 + (smi_[pos_ + 1] - '0');
                    pos_ += 2;
                } else {
                    if (!lenient_)
                        throw std::invalid_argument("Invalid %nn ring closure");
                    // lenient: rewind past '%' and break out of ring parsing
                    pos_ = savedPos;
                    break;
                }
            } else if (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                ringNum = smi_[pos_] - '0';
                pos_++;
            } else {
                // No ring closure; rewind bond symbol if any
                pos_ = savedPos;
                break;
            }

            // Process ring closure
            auto it = ringOpenings_.find(ringNum);
            if (it == ringOpenings_.end()) {
                // Open ring
                ringOpenings_[ringNum] = std::make_tuple(atomIdx, ringBondOrder, ringBondStereo);
            } else {
                // Close ring
                auto [openAtom, openBo, openStereo] = it->second;
                ringOpenings_.erase(it);

                int bo = 0, st = 0;
                // Reconcile bond orders from open and close
                if (ringBondOrder != 0 && openBo != 0) {
                    if (ringBondOrder != openBo) {
                        if (!lenient_)
                            throw std::invalid_argument("Conflicting ring closure bond orders");
                        // lenient: default to single bond
                        bo = 1;
                    } else {
                        bo = ringBondOrder;
                    }
                } else if (ringBondOrder != 0) {
                    bo = ringBondOrder;
                } else if (openBo != 0) {
                    bo = openBo;
                }

                // Reconcile stereo
                if (ringBondStereo != 0) st = ringBondStereo;
                else if (openStereo != 0) st = openStereo;

                if (bo == 0) {
                    // Default: aromatic if both aromatic, else single
                    if (atoms_[openAtom].isAromatic && atoms_[atomIdx].isAromatic)
                        bo = 4;
                    else
                        bo = 1;
                }

                bonds_.push_back({openAtom, atomIdx, bo, st});
            }
        }
    }

    void closeRings() {
        if (!ringOpenings_.empty()) {
            if (!lenient_)
                throw std::invalid_argument(
                    "Unclosed ring(s): " + std::to_string(ringOpenings_.begin()->first));
            // lenient: silently discard unclosed ring openings
            ringOpenings_.clear();
        }
    }
};

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
    int maxTarget = maxVal + std::abs(charge); // upper bound across all heuristics
    if (ebo > maxTarget) return 0;

    for (int v : valences) {
        int target;
        // Cations of group-15/16 elements expand into the next valence shell
        // (N+->4, P+->4, O+->3, S+->3, As+->4, Se+->3).
        // All other cases (anions, metals, etc.) simply lose available valence.
        if (charge > 0 && (atomicNum == 7  || atomicNum == 15 ||
                           atomicNum == 8  || atomicNum == 16 ||
                           atomicNum == 33 || atomicNum == 34)) {
            target = v + charge;
        } else if (charge < 0 && (atomicNum == 5 || atomicNum == 13)) {
            // Group 13 anions (B⁻, Al⁻) gain valence capacity,
            // becoming isoelectronic with Group 14 (e.g., [BH4-] has valence 4)
            target = v + std::abs(charge);
        } else {
            target = v - std::abs(charge);
        }
        int hc = target - ebo;
        if (hc >= 0) return hc;
    }

    return 0;
}

// --------------------------------------------------------------------------
// Detect ring membership via DFS
// --------------------------------------------------------------------------
inline void detectRings(int n,
                        const std::vector<std::vector<int>>& adj,
                        std::vector<bool>& atomInRing,
                        std::vector<std::set<int>>& bondInRing) {
    atomInRing.assign(n, false);
    bondInRing.resize(n);

    if (n == 0) return;

    // DFS to find back edges, then trace ring
    std::vector<int> parent(n, -1);
    std::vector<int> depth(n, -1);
    std::vector<bool> visited(n, false);

    // Use iterative DFS with explicit stack
    // For each back-edge, trace the ring path
    for (int start = 0; start < n; start++) {
        if (visited[start]) continue;

        struct Frame {
            int node;
            int parentNode;
            int neighborIdx;
        };
        std::vector<Frame> stack;
        stack.push_back({start, -1, 0});
        visited[start] = true;
        depth[start] = 0;
        parent[start] = -1;

        while (!stack.empty()) {
            auto& f = stack.back();
            if (f.neighborIdx >= static_cast<int>(adj[f.node].size())) {
                stack.pop_back();
                continue;
            }

            int nb = adj[f.node][f.neighborIdx];
            f.neighborIdx++;

            if (nb == f.parentNode) continue;

            if (visited[nb]) {
                // Back edge: ring found
                if (depth[nb] < depth[f.node]) {
                    // Trace ring from f.node back to nb
                    int cur = f.node;
                    while (cur != nb) {
                        atomInRing[cur] = true;
                        int p = parent[cur];
                        if (p >= 0) {
                            bondInRing[cur].insert(p);
                            bondInRing[p].insert(cur);
                        }
                        cur = p;
                    }
                    atomInRing[nb] = true;
                    // Also mark the back-edge itself
                    bondInRing[f.node].insert(nb);
                    bondInRing[nb].insert(f.node);
                }
            } else {
                visited[nb] = true;
                depth[nb] = depth[f.node] + 1;
                parent[nb] = f.node;
                stack.push_back({nb, f.node, 0});
            }
        }
    }
}

// --------------------------------------------------------------------------
// Build MolGraph from parsed SMILES
// --------------------------------------------------------------------------
inline MolGraph buildMolGraph(const std::vector<ParsedAtom>& atoms,
                               const std::vector<ParsedBond>& bonds,
                               bool lenient = false) {
    int n = static_cast<int>(atoms.size());
    if (n == 0) {
        if (!lenient) throw std::invalid_argument("SMILES produced no atoms");
        // lenient: return empty MolGraph
        MolGraph g;
        g.n = 0;
        return g;
    }

    // Build adjacency lists and bond properties
    std::vector<std::vector<int>> neighbors(n);
    std::vector<std::vector<int>> bondOrders(n);

    // Also build a simple adjacency for ring detection
    std::vector<std::vector<int>> simpleAdj(n);

    // Track bond order sums for implicit H calculation
    std::vector<int> bondOrderSum(n, 0);

    // Double bond stereo tracking
    // We'll store stereo info on bonds: dbStereo[i][j] = 1 (/) or 2 (\)
    // Actually dbStereoConf in MolGraph stores config per atom pair for
    // the C=C bond endpoints. We need to track / and \ on single bonds.

    for (auto& b : bonds) {
        int i = b.from, j = b.to;
        int ord = b.order;
        int effectiveOrd = (ord == 4) ? 1 : ord;  // aromatic = 1 for bond order sum

        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        bondOrders[i].push_back(ord);
        bondOrders[j].push_back(ord);

        simpleAdj[i].push_back(j);
        simpleAdj[j].push_back(i);

        bondOrderSum[i] += effectiveOrd;
        bondOrderSum[j] += effectiveOrd;
    }

    // Ring detection
    std::vector<bool> atomInRing;
    std::vector<std::set<int>> bondInRingSet;
    detectRings(n, simpleAdj, atomInRing, bondInRingSet);

    // Build ring/aromatic flags per bond
    std::vector<std::vector<bool>> bondRingFlags(n);
    std::vector<std::vector<bool>> bondAromFlags(n);

    for (int i = 0; i < n; i++) {
        bondRingFlags[i].resize(neighbors[i].size(), false);
        bondAromFlags[i].resize(neighbors[i].size(), false);
        for (int k = 0; k < static_cast<int>(neighbors[i].size()); k++) {
            int j = neighbors[i][k];
            bool inRing = bondInRingSet[i].count(j) > 0;
            bondRingFlags[i][k] = inRing;
            bool arom = (bondOrders[i][k] == 4);
            bondAromFlags[i][k] = arom;
        }
    }

    // Compute implicit H counts and build atom arrays
    std::vector<int> atomicNums(n), charges(n, 0), isotopes(n, 0), hydrogenCounts(n, 0), atomClasses(n, 0);
    std::vector<uint8_t> aromaticFlags(n, uint8_t(0)), ringFlags(n, uint8_t(0));
    std::vector<int> chirality(n, 0);

    for (int i = 0; i < n; i++) {
        const auto& a = atoms[i];
        atomicNums[i] = a.atomicNum;
        charges[i] = a.charge;
        isotopes[i] = a.isotope;
        atomClasses[i] = a.atomClass;
        aromaticFlags[i] = a.isAromatic;
        ringFlags[i] = atomInRing[i];
        chirality[i] = a.chirality;

        // Implicit H: bracket atoms have explicit hcount, organic atoms compute from valence
        if (a.hcount == -1) {
            // Organic subset atom: compute implicit H
            // Bond order sum: for aromatic bonds count 1, for others count the order
            // Actually we already have bondOrderSum computed correctly
            // No need to adjust -- computeImplicitH handles aromatic flag
            hydrogenCounts[i] = computeImplicitH(
                a.atomicNum, a.isAromatic, bondOrderSum[i], a.charge);
        } else {
            hydrogenCounts[i] = a.hcount;
        }
    }

    // Build double bond stereo configuration matrix
    // For / and \ bonds on C=C systems, determine E/Z
    // dbStereoConf[i][j] = 1 for Z (cis), 2 for E (trans) on the C=C bond
    std::vector<std::vector<int>> dbStereoMat;
    bool hasAnyStereo = false;

    // Collect / and \ annotations
    struct StereoAnnotation {
        int from, to, stereo;
    };
    std::vector<StereoAnnotation> stereoAnnotations;
    for (auto& b : bonds) {
        if (b.stereo != 0) {
            stereoAnnotations.push_back({b.from, b.to, b.stereo});
            hasAnyStereo = true;
        }
    }

    if (hasAnyStereo) {
        dbStereoMat.assign(n, std::vector<int>(n, 0));

        // Find C=C bonds and determine configuration from / \ on adjacent single bonds
        for (auto& b : bonds) {
            if (b.order == 2) {
                int a1 = b.from, a2 = b.to;
                // Find stereo bonds adjacent to a1 and a2
                int s1 = 0, s2 = 0; // stereo values for neighbors of a1, a2
                int n1 = -1, n2 = -1; // the neighbor atoms

                for (auto& sa : stereoAnnotations) {
                    if (sa.from == a1 && sa.to != a2) { s1 = sa.stereo; n1 = sa.to; }
                    if (sa.to == a1 && sa.from != a2) { s1 = sa.stereo; n1 = sa.from; }
                    if (sa.from == a2 && sa.to != a1) { s2 = sa.stereo; n2 = sa.to; }
                    if (sa.to == a2 && sa.from != a1) { s2 = sa.stereo; n2 = sa.from; }
                }

                if (s1 != 0 && s2 != 0) {
                    // Same direction = Z (cis), different = E (trans)
                    // / means the bond goes up from left to right
                    // \ means the bond goes down from left to right
                    // Two / or two \ = trans (E), / + \ = cis (Z)
                    if (s1 == s2)
                        dbStereoMat[a1][a2] = dbStereoMat[a2][a1] = 2; // E (trans)
                    else
                        dbStereoMat[a1][a2] = dbStereoMat[a2][a1] = 1; // Z (cis)
                }
            }
        }
    }

    // Use the Builder to construct MolGraph
    MolGraph::Builder builder;
    builder.atomCount(n)
           .atomicNumbers(atomicNums)
           .formalCharges(charges)
           .massNumbers(isotopes)
           .hydrogenCounts(hydrogenCounts)
           .atomClasses(atomClasses)
           .ringFlags(ringFlags)
           .aromaticFlags(aromaticFlags)
           .setNeighbors(neighbors)
           .setBondOrders(bondOrders)
           .bondRingFlags(bondRingFlags)
           .bondAromaticFlags(bondAromFlags)
           .tetrahedralChirality(chirality);

    if (hasAnyStereo && !dbStereoMat.empty()) {
        builder.doubleBondStereo(dbStereoMat);
    }

    return builder.build();
}

// --------------------------------------------------------------------------
// Canonical SMILES writer internals
// --------------------------------------------------------------------------

// Get the aromatic lowercase symbol for an element (returns "" if not applicable)
inline std::string aromaticSymbol(int z) {
    switch (z) {
        case 5:  return "b";
        case 6:  return "c";
        case 7:  return "n";
        case 8:  return "o";
        case 15: return "p";
        case 16: return "s";
        case 34: return "se";
        case 52: return "te";
        default: return "";
    }
}

// Determine if atom needs brackets
inline bool needsBrackets(const MolGraph& g, int idx) {
    int z = g.atomicNum[idx];
    if (!isOrganicSubset(z)) return true;
    if (g.formalCharge[idx] != 0) return true;
    if (g.massNumber[idx] != 0) return true;
    if (!g.atomClass.empty() && g.atomClass[idx] != 0) return true;
    if (g.tetraChirality[idx] != 0) return true;
    // If aromatic, check if it's one of the standard aromatic atoms
    if (g.aromatic[idx]) {
        std::string asym = aromaticSymbol(z);
        if (asym.empty()) return true;
    }
    return false;
}

// Write atom SMILES
inline std::string writeAtom(const MolGraph& g, int idx, int implicitH) {
    bool brackets = needsBrackets(g, idx);
    // v7.1.1: force bracket form for aromatic N / P carrying implicit H,
    // so pyrrole-type centres emit `[nH]` instead of `n`.  The organic-
    // subset symbol `n` defaults to "pyridine-type, 0 H", and RDKit
    // refuses to re-kekulize 5-ring aromatic-N systems when the H count
    // is missing.  Benchmarked on CMNPD 1 000 sample: 53 round-trip
    // kekulize failures before the fix, 0 after.
    if (!brackets && implicitH > 0 && g.aromatic[idx]
        && (g.atomicNum[idx] == 7 || g.atomicNum[idx] == 15)) {
        brackets = true;
    }
    std::ostringstream oss;

    if (brackets) {
        oss << '[';
        if (g.massNumber[idx] != 0) oss << g.massNumber[idx];

        if (g.atomicNum[idx] == 0 && !g.atomClass.empty() && g.atomClass[idx] > 0
            && g.formalCharge[idx] == 0 && g.massNumber[idx] == 0
            && g.tetraChirality[idx] == 0 && implicitH == 0) {
            oss << 'R' << g.atomClass[idx];
        } else if (g.aromatic[idx]) {
            std::string asym = aromaticSymbol(g.atomicNum[idx]);
            if (!asym.empty())
                oss << asym;
            else
                oss << atomicNumToSymbol(g.atomicNum[idx]);
        } else {
            oss << atomicNumToSymbol(g.atomicNum[idx]);
        }

        if (g.tetraChirality[idx] == 1) oss << '@';
        else if (g.tetraChirality[idx] == 2) oss << "@@";

        if (implicitH > 0) {
            oss << 'H';
            if (implicitH > 1) oss << implicitH;
        }

        if (g.formalCharge[idx] > 0) {
            oss << '+';
            if (g.formalCharge[idx] > 1) oss << g.formalCharge[idx];
        } else if (g.formalCharge[idx] < 0) {
            oss << '-';
            if (g.formalCharge[idx] < -1) oss << -g.formalCharge[idx];
        }

        if (!g.atomClass.empty() && g.atomClass[idx] != 0 && g.atomicNum[idx] != 0) {
            oss << ':' << g.atomClass[idx];
        } else if (!g.atomClass.empty() && g.atomClass[idx] != 0 && g.atomicNum[idx] == 0
                   && !(g.formalCharge[idx] == 0 && g.massNumber[idx] == 0
                        && g.tetraChirality[idx] == 0 && implicitH == 0)) {
            oss << ':' << g.atomClass[idx];
        }

        oss << ']';
    } else {
        if (g.aromatic[idx]) {
            oss << aromaticSymbol(g.atomicNum[idx]);
        } else {
            oss << atomicNumToSymbol(g.atomicNum[idx]);
        }
    }

    return oss.str();
}

// Bond order symbol
inline std::string bondSymbol(int order, bool fromAromatic, bool toAromatic) {
    if (order == 4) return ""; // aromatic implicit
    if (order == 1) {
        if (fromAromatic && toAromatic) return "-"; // explicit single between aromatics
        return ""; // implicit single
    }
    if (order == 2) return "=";
    if (order == 3) return "#";
    return "";
}

// Canonical SMILES writer using DFS with canonical ordering.
//
// Strategy:
//   1. Run a preliminary DFS to identify all back-edges (ring closures)
//      and assign ring numbers.  For each ring bond, record the ring
//      number at both endpoint atoms.
//   2. Run a second DFS that emits atom symbols, ring digits, bonds,
//      and branches.
inline std::string writeCanonicalSMILES(const MolGraph& g) {
    if (g.n == 0) return "";

    g.ensureCanonical();   // canonicalLabel is lazy
    int n = g.n;
    const auto& canon = g.canonicalLabel;

    // Compute implicit H for each atom
    std::vector<int> implicitH(n, 0);
    for (int i = 0; i < n; i++) {
        if (!isOrganicSubset(g.atomicNum[i]) && g.formalCharge[i] == 0 &&
            g.massNumber[i] == 0 && g.tetraChirality[i] == 0 && !g.aromatic[i])
            continue;
        int boSum = 0;
        for (int j : g.neighbors[i]) {
            int bo = g.bondOrder(i, j);
            if (bo == 4) bo = 1;
            boSum += bo;
        }
        implicitH[i] = computeImplicitH(
            g.atomicNum[i], g.aromatic[i], boSum, g.formalCharge[i]);
    }

    // v7.1.1 canonical-SMILES pyrrole-type N inference.
    //
    // computeImplicitH above uses a valence model that cannot distinguish
    // pyrrole-type aromatic N / P (contributes 2 pi electrons from its
    // lone pair, carries 1 implicit H) from pyridine-type aromatic N / P
    // (contributes 1 pi electron, carries 0 implicit H).  For aromatic
    // N with 2 ring σ-bonds the model lands on the pyridine interpretation
    // (ebo = 2 + 1 = 3 = default valence, hc = 0), which writes the
    // canonical SMILES as `n` even for true pyrroles.  RDKit then refuses
    // to re-kekulize the 5-ring aromatic system (4 C + 1 pyridine-N gives
    // 5 pi electrons, not Hückel 4n+2 = 6).
    //
    // Two-step fix:
    //
    //   (1) If the MolGraph has an explicit per-atom hydrogenCount
    //       stored from the input source (SMILES [nH], SDF M  HCOUNT,
    //       or CDK IAtom.getImplicitHydrogenCount), honour it for
    //       aromatic N / P.  This covers imidazole, pyrazole,
    //       N-methyl pyrroles, and every case where the input
    //       supplied the H count directly.
    //
    //   (2) For aromatic N / P that still have implicitH == 0 AND
    //       whose ring degree is exactly 2 (unsubstituted, so there
    //       is room for an implicit H), check the SSSR.  If the atom
    //       appears in any 5-membered ring whose only aromatic
    //       heteroatom is itself, it must be pyrrole-type; set
    //       implicitH[i] = 1 so the writer emits `[nH]`.  This
    //       handles pyrrole, indole, carbazole, 7-azaindole and
    //       fused benzo-pyrroles — the dominant failure mode on
    //       CMNPD 2D SDF entries that lack an explicit H count.
    //
    // For neutral aromatic N / P the generic valence fallback inside
    // computeImplicitH is unreliable: when σ-bond sum exceeds the
    // default valence of 3, it falls through to v = 5 and returns a
    // spurious "implicit H = 1" for N-substituted pyrroles (e.g.
    // N-methyl pyrrole `Cn1cccc1` → `C[nH]1cccc1`).  Override the
    // computed value unconditionally:
    //
    //   * degree 3 (or more) → 0 H (σ slots saturated; the lone pair
    //     donates the 2 pi electrons, no H needed).
    //   * degree 2 and in a 5-ring whose only aromatic heteroatom is
    //     this one → 1 H (pyrrole-type: pyrrole, indole, carbazole,
    //     7-azaindole, benzo-pyrroles).
    //   * anything else → leave whatever computeImplicitH produced.
    //
    // Explicit input H count (`g.hydrogenCount`, populated by SMILES
    // `[nH]`, SDF `M  HCOUNT`, or CDK `IAtom.getImplicitHydrogenCount`)
    // wins over both the valence model and the ring-topology fallback.
    for (int i = 0; i < n; i++) {
        if (!g.aromatic[i]) continue;
        int z = g.atomicNum[i];
        if (z != 7 && z != 15) continue;
        if (!g.formalCharge.empty() && g.formalCharge[i] != 0) continue;
        int deg = static_cast<int>(g.neighbors[i].size());
        // (1) Degree ≥ 3: no room for implicit H.  Override regardless
        //     of what the stored hydrogenCount says (some parsers over-
        //     eagerly infer 1 H for organic-subset `n` in aromatic
        //     5-rings even when the N is substituted).
        if (deg >= 3) {
            implicitH[i] = 0;
            continue;
        }
        // (2) Degree == 2: honour an explicit input H count if the
        //     caller supplied one (SMILES `[nH]`, SDF M  HCOUNT,
        //     CDK IAtom.getImplicitHydrogenCount).
        if (!g.hydrogenCount.empty() && g.hydrogenCount[i] > 0) {
            implicitH[i] = g.hydrogenCount[i];
            continue;
        }
        // (3) Degree == 2 with no stored H: fall back to the 5-ring /
        //     single-heteroatom SSSR inference for unsubstituted
        //     aromatic N / P (pyrrole, indole, carbazole,
        //     benzo-pyrroles).
        implicitH[i] = 0;
        const auto& rings = g.computeRings();
        for (const auto& r : rings) {
            if (r.size() != 5) continue;
            bool containsI = false;
            int aromaticHetero = 0;
            for (int ra : r) {
                if (ra == i) containsI = true;
                if (ra < n && g.aromatic[ra] && g.atomicNum[ra] != 6) {
                    ++aromaticHetero;
                }
            }
            if (containsI && aromaticHetero == 1) {
                implicitH[i] = 1;
                break;
            }
        }
    }

    // Find connected components
    std::vector<bool> compVis(n, false);
    std::vector<std::vector<int>> components;
    for (int i = 0; i < n; i++) {
        if (compVis[i]) continue;
        std::vector<int> comp;
        std::deque<int> q;
        q.push_back(i);
        compVis[i] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop_front();
            comp.push_back(u);
            for (int v : g.neighbors[u]) {
                if (!compVis[v]) { compVis[v] = true; q.push_back(v); }
            }
        }
        components.push_back(std::move(comp));
    }

    // Per-atom ring-digit annotations:
    // {ringNumber, bondOrder, partner}
    struct RingDigit { int rn; int bo; int partner; };
    std::vector<std::vector<RingDigit>> ringDigits(n);

    // DFS parent map for tree identification
    std::vector<int> dfsParent(n, -1);

    std::ostringstream result;

    auto sortedNeighbors = [&](int atom) {
        std::vector<int> nbs = g.neighbors[atom];
        std::sort(nbs.begin(), nbs.end(), [&](int a, int b) {
            return canon[a] < canon[b];
        });
        return nbs;
    };

    for (size_t ci = 0; ci < components.size(); ci++) {
        if (ci > 0) result << '.';
        auto& comp = components[ci];

        int root = comp[0];
        for (int a : comp) {
            if (canon[a] < canon[root]) root = a;
        }

        // --- Phase 1: DFS to find tree and assign ring numbers ---
        std::vector<bool> vis(n, false);
        std::set<int64_t> treeEdges; // set of bondKey for tree edges
        std::set<int64_t> assignedRings; // prevent double-counting same bond
        int ringCounter = 1;

        std::function<void(int, int)> dfs1 = [&](int atom, int parent) {
            vis[atom] = true;
            dfsParent[atom] = parent;
            auto nbs = sortedNeighbors(atom);
            bool parentSkipped = false;
            for (int nb : nbs) {
                if (nb == parent && !parentSkipped) { parentSkipped = true; continue; }
                if (vis[nb]) {
                    // Back-edge: each bond gets exactly one ring-closure number.
                    // Guard with a bond-key set to avoid assigning from both ends.
                    int64_t bk = MolGraph::bondKey(atom, nb);
                    if (assignedRings.insert(bk).second) {
                        int rn = ringCounter++;
                        int bo = g.bondOrder(atom, nb);
                        ringDigits[atom].push_back({rn, bo, nb});
                        ringDigits[nb].push_back({rn, bo, atom});
                    }
                } else {
                    treeEdges.insert(MolGraph::bondKey(atom, nb));
                    dfs1(nb, atom);
                }
            }
        };
        dfs1(root, -1);

        // --- Phase 2: DFS to emit SMILES ---
        vis.assign(n, false);
        std::string smiStr;
        // Track which ring digits have been emitted (first occurrence opens,
        // second occurrence closes).
        std::unordered_set<int> emittedRings;

        std::function<void(int, int)> dfs2 = [&](int atom, int parent) {
            vis[atom] = true;
            smiStr += writeAtom(g, atom, implicitH[atom]);

            // Emit ring digits for this atom
            for (auto& rd : ringDigits[atom]) {
                bool isOpen = (emittedRings.find(rd.rn) == emittedRings.end());
                if (isOpen) {
                    emittedRings.insert(rd.rn);
                }
                // Write bond symbol if needed
                std::string bs = bondSymbol(rd.bo, g.aromatic[atom], g.aromatic[rd.partner]);
                smiStr += bs;
                if (rd.rn < 10) smiStr += std::to_string(rd.rn);
                else smiStr += "%" + std::to_string(rd.rn);
            }

            // Collect tree children (unvisited neighbors connected by tree edges)
            auto nbs = sortedNeighbors(atom);
            std::vector<int> children;
            for (int nb : nbs) {
                if (vis[nb]) continue;
                int64_t key = MolGraph::bondKey(atom, nb);
                if (treeEdges.count(key)) children.push_back(nb);
            }

            for (size_t i = 0; i < children.size(); i++) {
                int nb = children[i];
                if (vis[nb]) continue;

                int bo = g.bondOrder(atom, nb);
                std::string bs = bondSymbol(bo, g.aromatic[atom], g.aromatic[nb]);

                if (i < children.size() - 1) {
                    smiStr += "(" + bs;
                    dfs2(nb, atom);
                    smiStr += ")";
                } else {
                    smiStr += bs;
                    dfs2(nb, atom);
                }
            }
        };

        dfs2(root, -1);
        result << smiStr;
    }

    return result.str();
}

// --------------------------------------------------------------------------
// Canonical SMARTS writer internals
// --------------------------------------------------------------------------

// Write a single atom as a SMARTS primitive.
// Always uses [#Z] (atomic number) as the base — unambiguous across toolkits.
// Aromaticity, charge, ring membership, H count, and isotope are added as
// semicolon-separated AND primitives when the corresponding option is set.
inline std::string writeAtomSMARTS(const MolGraph& g, int idx, int implH,
                                   bool inclArom, bool inclCharge,
                                   bool inclRing, bool inclHCount,
                                   bool inclIsotope) {
    int z = g.atomicNum[idx];
    std::ostringstream oss;
    oss << '[';

    // Isotope prefix (e.g. [2#1] for deuterium)
    if (inclIsotope && g.massNumber[idx] != 0)
        oss << g.massNumber[idx];

    oss << '#' << z;

    // Aromaticity: ;a = aromatic, ;A = aliphatic
    if (inclArom) {
        if (g.aromatic[idx]) oss << ";a";
        else                 oss << ";A";
    }

    // H count: ;H0, ;H1, ;H2, ...
    if (inclHCount && implH >= 0) {
        oss << ";H" << implH;
    }

    // Formal charge: ;+, ;+2, ;-, ;-2, ...
    if (inclCharge && g.formalCharge[idx] != 0) {
        int q = g.formalCharge[idx];
        if (q > 0) { oss << ";+"; if (q > 1) oss << q; }
        else       { oss << ";-"; if (-q > 1) oss << -q; }
    }

    // Ring membership: ;R (in ring) or ;!R (chain)
    if (inclRing) {
        oss << (g.ring[idx] ? ";R" : ";!R");
    }

    if (!g.atomClass.empty() && g.atomClass[idx] != 0) {
        oss << ':' << g.atomClass[idx];
    }

    oss << ']';
    return oss.str();
}

// Bond symbol for SMARTS — always written explicitly so output is unambiguous.
// Bond order 4 is the internal code for aromatic bonds.
inline std::string bondSymbolSMARTS(int bo) {
    switch (bo) {
        case 1:  return "-";
        case 2:  return "=";
        case 3:  return "#";
        case 4:  return ":";   // aromatic
        default: return "~";   // any — fallback for unusual bond types
    }
}

// Canonical SMARTS writer. Uses the same Morgan-based canonical labeling and
// DFS traversal order as writeCanonicalSMILES, so two isomorphic MolGraphs
// always produce the same SMARTS string (canonical + invariant).
inline std::string writeCanonicalSMARTS(const MolGraph& g,
                                        bool inclArom, bool inclCharge,
                                        bool inclRing, bool inclHCount,
                                        bool inclIsotope) {
    if (g.n == 0) return "";

    g.ensureCanonical();   // canonicalLabel is lazy
    int n = g.n;
    const auto& canon = g.canonicalLabel;

    // Implicit H per atom (needed if inclHCount is true; computed regardless
    // to keep the logic simple — it is cheap)
    std::vector<int> implicitH(n, 0);
    for (int i = 0; i < n; i++) {
        int boSum = 0;
        for (int j : g.neighbors[i]) {
            int bo = g.bondOrder(i, j);
            if (bo == 4) bo = 1;
            boSum += bo;
        }
        implicitH[i] = computeImplicitH(
            g.atomicNum[i], g.aromatic[i], boSum, g.formalCharge[i]);
    }

    // Connected components
    std::vector<bool> compVis(n, false);
    std::vector<std::vector<int>> components;
    for (int i = 0; i < n; i++) {
        if (compVis[i]) continue;
        std::vector<int> comp;
        std::deque<int> bfsQ;
        bfsQ.push_back(i); compVis[i] = true;
        while (!bfsQ.empty()) {
            int u = bfsQ.front(); bfsQ.pop_front();
            comp.push_back(u);
            for (int v : g.neighbors[u])
                if (!compVis[v]) { compVis[v] = true; bfsQ.push_back(v); }
        }
        components.push_back(std::move(comp));
    }

    // Neighbour iterator in canonical order
    auto sortedNbs = [&](int atom) {
        std::vector<int> nbs = g.neighbors[atom];
        std::sort(nbs.begin(), nbs.end(),
                  [&](int a, int b){ return canon[a] < canon[b]; });
        return nbs;
    };

    struct RingDigit { int rn; int bo; };
    std::vector<std::vector<RingDigit>> ringDigits(n);
    std::ostringstream result;

    for (size_t ci = 0; ci < components.size(); ci++) {
        if (ci > 0) result << '.';
        auto& comp = components[ci];

        // Root = atom with lowest canonical label in this component
        int root = comp[0];
        for (int a : comp)
            if (canon[a] < canon[root]) root = a;

        // Phase 1: DFS to discover tree/back edges and assign ring-closure numbers.
        // Each bond gets exactly one ring-closure number; guard with assignedRings
        // to avoid recording the same back-edge from both its endpoints.
        std::vector<bool> vis(n, false);
        std::set<int64_t> treeEdges;
        std::set<int64_t> assignedRings;
        int ringCounter = 1;

        std::function<void(int, int)> dfs1 = [&](int atom, int parent) {
            vis[atom] = true;
            auto nbs = sortedNbs(atom);
            bool parentSkipped = false;
            for (int nb : nbs) {
                if (nb == parent && !parentSkipped) { parentSkipped = true; continue; }
                if (vis[nb]) {
                    int64_t bk = MolGraph::bondKey(atom, nb);
                    if (assignedRings.insert(bk).second) {
                        int rn = ringCounter++;
                        int bo = g.bondOrder(atom, nb);
                        ringDigits[atom].push_back({rn, bo});
                        ringDigits[nb].push_back({rn, bo});
                    }
                } else {
                    treeEdges.insert(MolGraph::bondKey(atom, nb));
                    dfs1(nb, atom);
                }
            }
        };
        dfs1(root, -1);

        // Phase 2: DFS to emit SMARTS tokens
        vis.assign(n, false);
        std::string smStr;
        std::unordered_set<int> emittedRings;

        std::function<void(int, int)> dfs2 = [&](int atom, int parent) {
            vis[atom] = true;
            smStr += writeAtomSMARTS(g, atom, implicitH[atom],
                                     inclArom, inclCharge, inclRing,
                                     inclHCount, inclIsotope);

            for (auto& rd : ringDigits[atom]) {
                if (emittedRings.find(rd.rn) == emittedRings.end())
                    emittedRings.insert(rd.rn);
                smStr += bondSymbolSMARTS(rd.bo);
                if (rd.rn < 10) smStr += std::to_string(rd.rn);
                else            smStr += "%" + std::to_string(rd.rn);
            }

            auto nbs = sortedNbs(atom);
            std::vector<int> children;
            for (int nb : nbs) {
                if (!vis[nb] && treeEdges.count(MolGraph::bondKey(atom, nb)))
                    children.push_back(nb);
            }

            for (size_t i = 0; i < children.size(); i++) {
                int nb = children[i];
                if (vis[nb]) continue;
                std::string bs = bondSymbolSMARTS(g.bondOrder(atom, nb));
                if (i < children.size() - 1) {
                    smStr += "(" + bs;
                    dfs2(nb, atom);
                    smStr += ")";
                } else {
                    smStr += bs;
                    dfs2(nb, atom);
                }
            }
        };
        dfs2(root, -1);
        result << smStr;
    }

    return result.str();
}

} // namespace detail

// ============================================================================
// Public API
// ============================================================================

/**
 * Options for toSMARTS() — controls which atom primitives are written.
 *
 * Defaults are tuned for scaffold/MCS queries: element + aromaticity +
 * formal charge. H count and ring membership are off by default because
 * they over-constrain most scaffold searches.
 *
 * Example — tighten for fingerprint-style queries:
 * @code
 * SmartsWriteOptions opts;
 * opts.includeHCount     = true;
 * opts.includeRingMember = true;
 * std::string smarts = smsd::toSMARTS(mol, opts);
 * @endcode
 */
struct SmartsWriteOptions {
    /** Emit ;a (aromatic) or ;A (aliphatic) primitive. Default: true. */
    bool includeAromaticity  = true;
    /** Emit formal-charge primitive (+/-). Default: true. */
    bool includeCharge       = true;
    /** Emit ;R (ring atom) or ;!R (chain atom) primitive. Default: false. */
    bool includeRingMember   = false;
    /** Emit ;H<n> implicit-H-count primitive. Default: false. */
    bool includeHCount       = false;
    /** Emit isotope mass prefix ([2#1] for deuterium). Default: false. */
    bool includeIsotope      = false;
};

/**
 * Parse a SMILES string into a MolGraph.
 *
 * Supports: organic subset atoms, bracket atoms (isotope, chirality,
 * hcount, charge, class), all bond types, ring closures (1-9, %10-%99),
 * branches, aromatic atoms, stereochemistry (@ @@ / \), dot disconnection.
 *
 * @param smiles  The SMILES string to parse.
 * @return        A fully initialized MolGraph with Morgan ranks,
 *                canonical labels, bit-parallel adjacency, etc.
 * @throws std::invalid_argument on malformed SMILES.
 */
/// Maximum allowed atom count to prevent memory exhaustion from crafted input.
constexpr int MAX_ATOMS = 10000;

inline MolGraph parseSMILES(const std::string& smiles) {
    detail::SmilesParser parser(smiles);
    parser.parse();
    if (static_cast<int>(parser.atoms().size()) > MAX_ATOMS)
        throw std::invalid_argument("Molecule exceeds MAX_ATOMS limit ("
            + std::to_string(MAX_ATOMS) + ")");
    return detail::buildMolGraph(parser.atoms(), parser.bonds());
}

/**
 * Parse a SMILES string into a MolGraph with configurable options.
 *
 * When opts.lenient is true, the parser performs best-effort recovery:
 *   - Unclosed rings are silently ignored (no bond created)
 *   - Unbalanced parentheses are silently handled
 *   - Unknown elements are skipped
 *   - Conflicting bond orders default to single bond
 *
 * @param smiles  The SMILES string to parse.
 * @param opts    Parsing options (default: strict mode).
 * @return        A MolGraph (may be empty in lenient mode if no atoms parsed).
 * @throws std::invalid_argument on malformed SMILES in strict mode.
 */
inline MolGraph parseSMILES(const std::string& smiles, const ParseOptions& opts) {
    detail::SmilesParser parser(smiles, opts.lenient);
    parser.parse();
    if (static_cast<int>(parser.atoms().size()) > MAX_ATOMS)
        throw std::invalid_argument("Molecule exceeds MAX_ATOMS limit ("
            + std::to_string(MAX_ATOMS) + ")");
    return detail::buildMolGraph(parser.atoms(), parser.bonds(), opts.lenient);
}

/**
 * Generate a canonical SMILES string from a MolGraph.
 *
 * Uses the canonical labeling computed by MolGraph::initDerivedFields()
 * to produce a deterministic SMILES output.
 *
 * @param g  The molecular graph.
 * @return   A canonical SMILES string.
 */
inline std::string toSMILES(const MolGraph& g) {
    return detail::writeCanonicalSMILES(g);
}

/**
 * Generate a canonical SMARTS string from a MolGraph.
 *
 * The output is invariant under atom renumbering: the Morgan-based canonical
 * labeling drives the DFS traversal in the same way as toSMILES(), so two
 * isomorphic graphs always produce the same SMARTS string.
 *
 * Atom primitives use the [#Z] form (atomic number) as the base, with
 * optional aromaticity, charge, H count, ring membership, and isotope
 * primitives controlled by SmartsWriteOptions. Bond symbols are always
 * written explicitly (-, =, #, :, ~) so there is no ambiguity.
 *
 * Typical use — write the scaffold of an MCS result:
 * @code
 * auto mcs = smsd.findMCS();
 * std::vector<int> atoms;
 * for (auto& [qi, ti] : mcs) atoms.push_back(qi);
 * auto sub = smsd::extractSubgraph(mol1, atoms);
 * std::string scaffold = smsd::toSMARTS(sub);   // e.g. [#6;a]:[#6;a]...
 * @endcode
 *
 * @param g     The molecular graph.
 * @param opts  Atom primitive options (default: element + aromaticity + charge).
 * @return      A canonical SMARTS string.
 */
inline std::string toSMARTS(const MolGraph& g,
                             SmartsWriteOptions opts = SmartsWriteOptions{}) {
    return detail::writeCanonicalSMARTS(g,
        opts.includeAromaticity,
        opts.includeCharge,
        opts.includeRingMember,
        opts.includeHCount,
        opts.includeIsotope);
}

} // namespace smsd

// ============================================================================
// Tests (compile with -DSMSD_TEST_SMILES)
// ============================================================================
#ifdef SMSD_TEST_SMILES

#include <iostream>

namespace smsd_test {

static int testsPassed = 0;
static int testsFailed = 0;

#define SMSD_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "FAIL: " << msg << " [" << __FILE__ << ":" << __LINE__ << "]" << std::endl; \
            testsFailed++; \
        } else { \
            testsPassed++; \
        } \
    } while(0)

#define SMSD_ASSERT_EQ(actual, expected, msg) \
    do { \
        if ((actual) != (expected)) { \
            std::cerr << "FAIL: " << msg << " (expected " << (expected) \
                      << ", got " << (actual) << ") [" << __FILE__ << ":" \
                      << __LINE__ << "]" << std::endl; \
            testsFailed++; \
        } else { \
            testsPassed++; \
        } \
    } while(0)

#define SMSD_ASSERT_THROWS(expr, msg) \
    do { \
        bool threw = false; \
        try { expr; } catch (const std::invalid_argument&) { threw = true; } \
        if (!threw) { \
            std::cerr << "FAIL: " << msg << " (expected exception) [" \
                      << __FILE__ << ":" << __LINE__ << "]" << std::endl; \
            testsFailed++; \
        } else { \
            testsPassed++; \
        } \
    } while(0)

inline int countBonds(const smsd::MolGraph& g) {
    int count = 0;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) count++;
        }
    }
    return count;
}

inline void testMethane() {
    auto g = smsd::parseSMILES("C");
    SMSD_ASSERT_EQ(g.n, 1, "C: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 6, "C: carbon");
    SMSD_ASSERT_EQ(g.degree[0], 0, "C: degree 0");
    SMSD_ASSERT_EQ(countBonds(g), 0, "C: 0 bonds");
}

inline void testEthane() {
    auto g = smsd::parseSMILES("CC");
    SMSD_ASSERT_EQ(g.n, 2, "CC: 2 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[0], 6, "CC: C0");
    SMSD_ASSERT_EQ(g.atomicNum[1], 6, "CC: C1");
    SMSD_ASSERT_EQ(countBonds(g), 1, "CC: 1 bond");
    SMSD_ASSERT_EQ(g.bondOrder(0, 1), 1, "CC: single bond");
}

inline void testDoubleBond() {
    auto g = smsd::parseSMILES("C=C");
    SMSD_ASSERT_EQ(g.n, 2, "C=C: 2 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 1, "C=C: 1 bond");
    SMSD_ASSERT_EQ(g.bondOrder(0, 1), 2, "C=C: double bond");
}

inline void testTripleBond() {
    auto g = smsd::parseSMILES("C#N");
    SMSD_ASSERT_EQ(g.n, 2, "C#N: 2 atoms");
    SMSD_ASSERT_EQ(g.bondOrder(0, 1), 3, "C#N: triple bond");
    SMSD_ASSERT_EQ(g.atomicNum[1], 7, "C#N: nitrogen");
}

inline void testBenzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    SMSD_ASSERT_EQ(g.n, 6, "benzene: 6 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 6, "benzene: 6 bonds");
    for (int i = 0; i < 6; i++) {
        SMSD_ASSERT(g.aromatic[i], "benzene: atom " + std::to_string(i) + " aromatic");
        SMSD_ASSERT_EQ(g.atomicNum[i], 6, "benzene: atom " + std::to_string(i) + " is carbon");
        SMSD_ASSERT(g.ring[i], "benzene: atom " + std::to_string(i) + " in ring");
    }
}

inline void testPhenol() {
    auto g = smsd::parseSMILES("c1ccc(O)cc1");
    SMSD_ASSERT_EQ(g.n, 7, "phenol: 7 atoms");
    // 6 ring bonds + 1 C-O bond = 7 bonds
    SMSD_ASSERT_EQ(countBonds(g), 7, "phenol: 7 bonds");
    // Find the oxygen
    int oIdx = -1;
    for (int i = 0; i < g.n; i++) {
        if (g.atomicNum[i] == 8) oIdx = i;
    }
    SMSD_ASSERT(oIdx >= 0, "phenol: has oxygen");
    SMSD_ASSERT(!g.aromatic[oIdx], "phenol: O not aromatic");
}

inline void testAmmonium() {
    auto g = smsd::parseSMILES("[NH4+]");
    SMSD_ASSERT_EQ(g.n, 1, "[NH4+]: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 7, "[NH4+]: nitrogen");
    SMSD_ASSERT_EQ(g.formalCharge[0], 1, "[NH4+]: charge +1");
}

inline void testIsotope() {
    auto g = smsd::parseSMILES("[13CH4]");
    SMSD_ASSERT_EQ(g.n, 1, "[13CH4]: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 6, "[13CH4]: carbon");
    SMSD_ASSERT_EQ(g.massNumber[0], 13, "[13CH4]: isotope 13");
}

inline void testAceticAcid() {
    auto g = smsd::parseSMILES("CC(=O)O");
    SMSD_ASSERT_EQ(g.n, 4, "acetic acid: 4 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 3, "acetic acid: 3 bonds");
    // Atom 0: C, Atom 1: C, Atom 2: O (=O), Atom 3: O
    SMSD_ASSERT_EQ(g.atomicNum[0], 6, "acetic acid: C0");
    SMSD_ASSERT_EQ(g.atomicNum[1], 6, "acetic acid: C1");
    SMSD_ASSERT_EQ(g.atomicNum[2], 8, "acetic acid: O2");
    SMSD_ASSERT_EQ(g.atomicNum[3], 8, "acetic acid: O3");
    SMSD_ASSERT_EQ(g.bondOrder(1, 2), 2, "acetic acid: C=O double bond");
    SMSD_ASSERT_EQ(g.bondOrder(1, 3), 1, "acetic acid: C-O single bond");
}

inline void testCyclopropane() {
    auto g = smsd::parseSMILES("C1CC1");
    SMSD_ASSERT_EQ(g.n, 3, "cyclopropane: 3 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 3, "cyclopropane: 3 bonds");
    for (int i = 0; i < 3; i++) {
        SMSD_ASSERT(g.ring[i], "cyclopropane: atom " + std::to_string(i) + " in ring");
    }
}

inline void testDotDisconnection() {
    auto g = smsd::parseSMILES("C.C");
    SMSD_ASSERT_EQ(g.n, 2, "C.C: 2 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 0, "C.C: 0 bonds (disconnected)");
}

inline void testEButene() {
    auto g = smsd::parseSMILES("C/C=C/C");
    SMSD_ASSERT_EQ(g.n, 4, "E-butene: 4 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 3, "E-butene: 3 bonds");
    SMSD_ASSERT_EQ(g.bondOrder(1, 2), 2, "E-butene: C=C double bond");
}

inline void testNaphthalene() {
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    SMSD_ASSERT_EQ(g.n, 10, "naphthalene: 10 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 11, "naphthalene: 11 bonds");
    for (int i = 0; i < 10; i++) {
        SMSD_ASSERT(g.aromatic[i], "naphthalene: atom " + std::to_string(i) + " aromatic");
        SMSD_ASSERT(g.ring[i], "naphthalene: atom " + std::to_string(i) + " in ring");
    }
}

inline void testAspirin() {
    auto g = smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O");
    SMSD_ASSERT_EQ(g.n, 13, "aspirin: 13 atoms");
    // Count carbons, oxygens
    int nC = 0, nO = 0;
    for (int i = 0; i < g.n; i++) {
        if (g.atomicNum[i] == 6) nC++;
        if (g.atomicNum[i] == 8) nO++;
    }
    SMSD_ASSERT_EQ(nC, 9, "aspirin: 9 carbons");
    SMSD_ASSERT_EQ(nO, 4, "aspirin: 4 oxygens");
}

inline void testIronBracket() {
    auto g = smsd::parseSMILES("[Fe+2]");
    SMSD_ASSERT_EQ(g.n, 1, "[Fe+2]: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 26, "[Fe+2]: iron");
    SMSD_ASSERT_EQ(g.formalCharge[0], 2, "[Fe+2]: charge +2");
}

inline void testAromaticNH() {
    auto g = smsd::parseSMILES("[nH]1cccc1");
    SMSD_ASSERT_EQ(g.n, 5, "pyrrole: 5 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[0], 7, "pyrrole: nitrogen");
    SMSD_ASSERT(g.aromatic[0], "pyrrole: N aromatic");
}

inline void testChirality() {
    auto g = smsd::parseSMILES("[C@@H](F)(Cl)Br");
    SMSD_ASSERT_EQ(g.n, 4, "chiral: 4 atoms");
    SMSD_ASSERT_EQ(g.tetraChirality[0], 2, "chiral: @@ = 2");
}

inline void testPercentRingClosure() {
    // Ring closure %10
    auto g = smsd::parseSMILES("C%10CCCCCCCCC%10");
    SMSD_ASSERT_EQ(g.n, 10, "%10 ring: 10 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 10, "%10 ring: 10 bonds (ring)");
}

inline void testBranching() {
    // Isobutane: C(C)(C)C
    auto g = smsd::parseSMILES("C(C)(C)C");
    SMSD_ASSERT_EQ(g.n, 4, "isobutane: 4 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 3, "isobutane: 3 bonds");
    SMSD_ASSERT_EQ(g.degree[0], 3, "isobutane: central C degree 3");
}

inline void testHalogens() {
    auto g = smsd::parseSMILES("ClBrFI");
    SMSD_ASSERT_EQ(g.n, 4, "ClBrFI: 4 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[0], 17, "Cl");
    SMSD_ASSERT_EQ(g.atomicNum[1], 35, "Br");
    SMSD_ASSERT_EQ(g.atomicNum[2], 9, "F");
    SMSD_ASSERT_EQ(g.atomicNum[3], 53, "I");
}

inline void testChloromethane() {
    auto g = smsd::parseSMILES("CCl");
    SMSD_ASSERT_EQ(g.n, 2, "CCl: 2 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[0], 6, "CCl: carbon");
    SMSD_ASSERT_EQ(g.atomicNum[1], 17, "CCl: chlorine");
}

inline void testNegativeCharge() {
    auto g = smsd::parseSMILES("[O-]");
    SMSD_ASSERT_EQ(g.n, 1, "[O-]: 1 atom");
    SMSD_ASSERT_EQ(g.formalCharge[0], -1, "[O-]: charge -1");
}

inline void testDoubleNegCharge() {
    auto g = smsd::parseSMILES("[O-2]");
    SMSD_ASSERT_EQ(g.n, 1, "[O-2]: 1 atom");
    SMSD_ASSERT_EQ(g.formalCharge[0], -2, "[O-2]: charge -2");
}

inline void testRepeatedCharge() {
    auto g = smsd::parseSMILES("[Fe++]");
    SMSD_ASSERT_EQ(g.formalCharge[0], 2, "[Fe++]: charge +2 (repeated)");
}

inline void testBoron() {
    auto g = smsd::parseSMILES("B");
    SMSD_ASSERT_EQ(g.n, 1, "B: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 5, "B: boron");
}

inline void testPhosphorus() {
    auto g = smsd::parseSMILES("P");
    SMSD_ASSERT_EQ(g.n, 1, "P: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 15, "P: phosphorus");
}

inline void testSulfur() {
    auto g = smsd::parseSMILES("S");
    SMSD_ASSERT_EQ(g.n, 1, "S: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 16, "S: sulfur");
}

inline void testWildcard() {
    auto g = smsd::parseSMILES("[*]");
    SMSD_ASSERT_EQ(g.n, 1, "[*]: 1 atom");
    SMSD_ASSERT_EQ(g.atomicNum[0], 0, "[*]: wildcard Z=0");
}

inline void testMultipleRings() {
    // Cubane: C12C3C4C1C5C4C3C25
    auto g = smsd::parseSMILES("C12C3C4C1C5C4C3C25");
    SMSD_ASSERT_EQ(g.n, 8, "cubane: 8 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 12, "cubane: 12 bonds");
}

inline void testLongChain() {
    auto g = smsd::parseSMILES("CCCCCCCCCC");
    SMSD_ASSERT_EQ(g.n, 10, "decane: 10 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 9, "decane: 9 bonds");
}

inline void testThiophene() {
    auto g = smsd::parseSMILES("c1ccsc1");
    SMSD_ASSERT_EQ(g.n, 5, "thiophene: 5 atoms");
    int sIdx = -1;
    for (int i = 0; i < g.n; i++) {
        if (g.atomicNum[i] == 16) sIdx = i;
    }
    SMSD_ASSERT(sIdx >= 0, "thiophene: has sulfur");
    SMSD_ASSERT(g.aromatic[sIdx], "thiophene: S aromatic");
}

inline void testPyridine() {
    auto g = smsd::parseSMILES("c1ccncc1");
    SMSD_ASSERT_EQ(g.n, 6, "pyridine: 6 atoms");
    int nIdx = -1;
    for (int i = 0; i < g.n; i++) {
        if (g.atomicNum[i] == 7) nIdx = i;
    }
    SMSD_ASSERT(nIdx >= 0, "pyridine: has nitrogen");
    SMSD_ASSERT(g.aromatic[nIdx], "pyridine: N aromatic");
}

inline void testInvalidSMILES() {
    SMSD_ASSERT_THROWS(smsd::parseSMILES(""), "empty SMILES");
    SMSD_ASSERT_THROWS(smsd::parseSMILES("C1CC"), "unclosed ring");
    SMSD_ASSERT_THROWS(smsd::parseSMILES("C(C"), "unclosed branch");
    SMSD_ASSERT_THROWS(smsd::parseSMILES("["), "unclosed bracket");
    SMSD_ASSERT_THROWS(smsd::parseSMILES("[Zq]"), "unknown element");
}

inline void testMorganRanks() {
    auto g = smsd::parseSMILES("c1ccccc1");
    g.ensureCanonical();   // morganRank is lazy
    // All atoms in benzene should have same morgan rank (by symmetry)
    SMSD_ASSERT(!g.morganRank.empty(), "benzene: morgan ranks computed");
    // They should all be equal for benzene
    for (int i = 1; i < g.n; i++) {
        SMSD_ASSERT_EQ(g.morganRank[i], g.morganRank[0],
                        "benzene: all morgan ranks equal");
    }
}

inline void testCanonicalLabels() {
    auto g = smsd::parseSMILES("c1ccccc1");
    g.ensureCanonical();   // canonicalLabel is lazy
    SMSD_ASSERT(!g.canonicalLabel.empty(), "benzene: canonical labels computed");
    SMSD_ASSERT_EQ(static_cast<int>(g.canonicalLabel.size()), 6,
                    "benzene: 6 canonical labels");
}

inline void testBitParallelAdjacency() {
    auto g = smsd::parseSMILES("CC");
    SMSD_ASSERT(!g.adjLong.empty(), "CC: adjLong computed");
    SMSD_ASSERT(g.adjLong[0][0] & (uint64_t(1) << 1), "CC: adj[0] has bit 1");
    SMSD_ASSERT(g.adjLong[1][0] & (uint64_t(1) << 0), "CC: adj[1] has bit 0");
}

inline void testSMILESWriter() {
    // Round-trip: parse then write, then re-parse and verify same structure
    auto g1 = smsd::parseSMILES("CC");
    std::string s1 = smsd::toSMILES(g1);
    auto g2 = smsd::parseSMILES(s1);
    SMSD_ASSERT_EQ(g2.n, 2, "writer round-trip CC: 2 atoms");
    SMSD_ASSERT_EQ(countBonds(g2), 1, "writer round-trip CC: 1 bond");

    // Benzene round-trip
    auto gb = smsd::parseSMILES("c1ccccc1");
    std::string sb = smsd::toSMILES(gb);
    auto gb2 = smsd::parseSMILES(sb);
    SMSD_ASSERT_EQ(gb2.n, 6, "writer round-trip benzene: 6 atoms");
    SMSD_ASSERT_EQ(countBonds(gb2), 6, "writer round-trip benzene: 6 bonds");
}

inline void testAtomClass() {
    auto g = smsd::parseSMILES("[CH3:1]C");
    SMSD_ASSERT_EQ(g.n, 2, "atom class: 2 atoms");
    SMSD_ASSERT_EQ(g.atomClass[0], 1, "atom class preserved on first atom");
    SMSD_ASSERT(smsd::toSMILES(g).find(":1") != std::string::npos,
                "atom class writes back with attachment label");
}

inline void testPatentRGroupPlaceholder() {
    auto g = smsd::parseSMILES("[R1]c1ccccc1");
    SMSD_ASSERT_EQ(g.n, 7, "R-group placeholder: 7 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[0], 0, "R-group placeholder stored as wildcard");
    SMSD_ASSERT_EQ(g.atomClass[0], 1, "R-group placeholder stored as class 1");
    SMSD_ASSERT_EQ(smsd::toSMILES(g), "[R1]c1ccccc1", "R-group placeholder round-trip");
}

inline void testCyclohexane() {
    auto g = smsd::parseSMILES("C1CCCCC1");
    SMSD_ASSERT_EQ(g.n, 6, "cyclohexane: 6 atoms");
    SMSD_ASSERT_EQ(countBonds(g), 6, "cyclohexane: 6 bonds");
    for (int i = 0; i < 6; i++) {
        SMSD_ASSERT(g.ring[i], "cyclohexane: atom in ring");
    }
}

inline void testDMSO() {
    // CS(=O)C  dimethyl sulfoxide
    auto g = smsd::parseSMILES("CS(=O)C");
    SMSD_ASSERT_EQ(g.n, 4, "DMSO: 4 atoms");
    SMSD_ASSERT_EQ(g.atomicNum[1], 16, "DMSO: sulfur");
}

inline void testCaffeine() {
    auto g = smsd::parseSMILES("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
    SMSD_ASSERT(g.n > 0, "caffeine: parsed OK");
    SMSD_ASSERT_EQ(g.n, 14, "caffeine: 14 heavy atoms");
}

// Benchmark suite molecules
inline void testBenchmarkMolecules() {
    // Methane
    auto g1 = smsd::parseSMILES("C");
    SMSD_ASSERT_EQ(g1.n, 1, "benchmark methane");

    // Water
    auto g2 = smsd::parseSMILES("O");
    SMSD_ASSERT_EQ(g2.n, 1, "benchmark water");

    // Ethanol
    auto g3 = smsd::parseSMILES("CCO");
    SMSD_ASSERT_EQ(g3.n, 3, "benchmark ethanol");

    // Benzene
    auto g4 = smsd::parseSMILES("c1ccccc1");
    SMSD_ASSERT_EQ(g4.n, 6, "benchmark benzene");

    // Naphthalene
    auto g5 = smsd::parseSMILES("c1ccc2ccccc2c1");
    SMSD_ASSERT_EQ(g5.n, 10, "benchmark naphthalene");

    // Aspirin
    auto g6 = smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O");
    SMSD_ASSERT_EQ(g6.n, 13, "benchmark aspirin");

    // Caffeine
    auto g7 = smsd::parseSMILES("Cn1cnc2c1c(=O)n(c(=O)n2C)C");
    SMSD_ASSERT_EQ(g7.n, 14, "benchmark caffeine");

    // Ibuprofen
    auto g8 = smsd::parseSMILES("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
    SMSD_ASSERT_EQ(g8.n, 15, "benchmark ibuprofen");

    // Glucose (open chain)
    auto g9 = smsd::parseSMILES("OCC(O)C(O)C(O)C(O)C=O");
    SMSD_ASSERT_EQ(g9.n, 12, "benchmark glucose");

    // Cholesterol
    auto g10 = smsd::parseSMILES("CC(CCCC(C)C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C");
    SMSD_ASSERT(g10.n > 20, "benchmark cholesterol parsed");

    // Strychnine
    auto g11 = smsd::parseSMILES("O=C1CC2OCC=C3CN4CCC5=CC=CC1C5C4CC23");
    SMSD_ASSERT(g11.n > 15, "benchmark strychnine parsed");

    // ATP (simplified)
    auto g12 = smsd::parseSMILES("c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N");
    SMSD_ASSERT(g12.n > 20, "benchmark ATP parsed");

    // Paclitaxel (taxol)
    auto g13 = smsd::parseSMILES("CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(c5ccccc5)NC(=O)c6ccccc6)O)O)OC(=O)C7CCCCC7)(CO4)OC(=O)C)O)C)OC(=O)C");
    SMSD_ASSERT(g13.n > 40, "benchmark paclitaxel parsed");

    // Vancomycin (large drug)
    auto g14 = smsd::parseSMILES("CC1C(C(CC(O1)OC2C(C(C(OC2OC3=CC=C4C(=C3Cl)OC5=C(C=C(C(=C5C(=O)NC(CC(=O)N)C(=O)NC4C(=O)NC6C7=CC(=C(C(=C7)OC8=CC=C(C(=C8)Cl)C9C(C(=O)NC(C(=O)N1CCCC1C(=O)NC(C(=O)NC6CO)CC(=O)N)CC(C)C)NC(=O)C(CC(C)C)NC9=O)O)O)O)C(=O)O)O)O)CO)O)O)(C)N)O");
    SMSD_ASSERT(g14.n > 80, "benchmark vancomycin parsed");
}

inline void runAllTests() {
    std::cout << "Running SMILES parser tests..." << std::endl;

    testMethane();
    testEthane();
    testDoubleBond();
    testTripleBond();
    testBenzene();
    testPhenol();
    testAmmonium();
    testIsotope();
    testAceticAcid();
    testCyclopropane();
    testDotDisconnection();
    testEButene();
    testNaphthalene();
    testAspirin();
    testIronBracket();
    testAromaticNH();
    testChirality();
    testPercentRingClosure();
    testBranching();
    testHalogens();
    testChloromethane();
    testNegativeCharge();
    testDoubleNegCharge();
    testRepeatedCharge();
    testBoron();
    testPhosphorus();
    testSulfur();
    testWildcard();
    testMultipleRings();
    testLongChain();
    testThiophene();
    testPyridine();
    testInvalidSMILES();
    testMorganRanks();
    testCanonicalLabels();
    testBitParallelAdjacency();
    testSMILESWriter();
    testAtomClass();
    testPatentRGroupPlaceholder();
    testCyclohexane();
    testDMSO();
    testCaffeine();
    testBenchmarkMolecules();

    std::cout << std::endl;
    std::cout << "Passed: " << testsPassed << std::endl;
    std::cout << "Failed: " << testsFailed << std::endl;
    std::cout << "Total:  " << testsPassed + testsFailed << std::endl;

    if (testsFailed > 0) {
        std::cout << "\n*** SOME TESTS FAILED ***" << std::endl;
    } else {
        std::cout << "\nAll tests passed." << std::endl;
    }
}

} // namespace smsd_test

int main() {
    smsd_test::runAllTests();
    return smsd_test::testsFailed > 0 ? 1 : 0;
}

#endif // SMSD_TEST_SMILES

#endif // SMSD_SMILES_PARSER_HPP

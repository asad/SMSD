/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only OpenSMILES parser and writer for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_SMILES_PARSER_HPP
#define SMSD_SMILES_PARSER_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/parser_utils.hpp"
#include "smsd/smiles_utils.hpp"
#include "smsd/stereo.hpp"

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
    bool hasIncomingBond = false;
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
                    pos_++; // skip truly unrecognised characters
                }
            }
        }
        closeRings();
    }

    const std::vector<ParsedAtom>& atoms() const { return atoms_; }
    const std::vector<ParsedBond>& bonds() const { return bonds_; }
    const std::vector<std::vector<int>>& neighborOrder() const { return neighborOrder_; }

private:
    std::string smi_;
    size_t pos_;
    bool lenient_;
    std::vector<ParsedAtom> atoms_;
    std::vector<ParsedBond> bonds_;
    std::vector<std::vector<int>> neighborOrder_;

    // Ring number -> atom, bond order, stereo, and neighbour-order placeholder.
    std::unordered_map<int, std::tuple<int, int, int, int>> ringOpenings_;

    void parseChain(int prevAtom, int pendingBondOrder, int pendingStereo) {
        std::vector<std::tuple<int, int, int>> branches;
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
                branches.emplace_back(prevAtom, pendingBondOrder, pendingStereo);
                pos_++;
                pendingBondOrder = 0;
                pendingStereo = 0;
                parseBondSymbol(pendingBondOrder, pendingStereo);
                continue;
            }

            // Branch close
            if (c == ')') {
                if (branches.empty()) return;
                std::tie(prevAtom, pendingBondOrder, pendingStereo) = branches.back();
                branches.pop_back();
                ++pos_;
                continue;
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
                if (!branches.empty() && !lenient_)
                    throw std::invalid_argument("Unmatched opening parenthesis");
                return;
            }

            // Connect to previous atom
            atoms_[atomIdx].hasIncomingBond = prevAtom >= 0;
            neighborOrder_.resize(atoms_.size());
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
                neighborOrder_[prevAtom].push_back(atomIdx);
                neighborOrder_[atomIdx].push_back(prevAtom);
            }
            pendingBondOrder = 0;
            pendingStereo = 0;

            // Ring closures
            parseRingClosures(atomIdx);

            prevAtom = atomIdx;
        }
        if (!branches.empty() && !lenient_)
            throw std::invalid_argument("Unmatched opening parenthesis");
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

    /**
     * Reads a numeric field of a bracket atom of at most `digits` digits, or returns -1 when no
     * digit follows. OpenSMILES allows three digits for an isotope, one for a hydrogen count and
     * two for a charge; a longer field is rejected, in lenient mode as well.
     */
    int readBracketField(size_t digits, const char* field, const char* limit) {
        const size_t start = pos_;
        const int value = parser_detail::readDecimal(smi_, pos_);
        if (pos_ - start > digits)
            throw std::invalid_argument(std::string(field) + " in bracket atom exceeds " + limit);
        return value;
    }

    int parseBracketAtom() {
        if (pos_ >= smi_.size() || smi_[pos_] != '[') return -1;
        pos_++; // skip '['

        ParsedAtom a;
        a.isBracket = true;

        // Isotope (optional digits, at most three)
        const int isotope = readBracketField(3, "Isotope", "three digits");
        if (isotope >= 0) a.isotope = isotope;

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
                a.atomClass = parser_detail::readDecimal(smi_, pos_);
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
            if (pos_ + 1 < smi_.size()) {
                const std::string stereoClass = smi_.substr(pos_, 2);
                if (stereoClass == "TH") {
                    pos_ += 2;
                    if (pos_ >= smi_.size() || (smi_[pos_] != '1' && smi_[pos_] != '2'))
                        throw std::invalid_argument("Tetrahedral chirality must be @TH1 or @TH2");
                    a.chirality = smi_[pos_++] - '0';
                } else if (stereoClass == "AL" || stereoClass == "SP" ||
                           stereoClass == "TB" || stereoClass == "OH") {
                    throw std::invalid_argument("Unsupported SMILES chirality class: " + stereoClass);
                }
            }
        }

        // Hydrogen count (optional): H alone is one hydrogen, and a single digit gives the count.
        if (pos_ < smi_.size() && smi_[pos_] == 'H') {
            pos_++;
            int hc = 1;
            if (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                hc = readBracketField(1, "Hydrogen count", "one digit");
            }
            a.hcount = hc;
        } else {
            a.hcount = 0;  // bracket atom with no H specified means 0
        }

        // Charge (optional): a sign with up to two digits, or repeated signs, within -15 to +15.
        if (pos_ < smi_.size() && (smi_[pos_] == '+' || smi_[pos_] == '-')) {
            char sign = smi_[pos_];
            pos_++;
            int chg = 1;
            if (pos_ < smi_.size() && std::isdigit(smi_[pos_])) {
                chg = readBracketField(2, "Charge", "two digits");
            } else {
                // Count repeated signs: ++ = +2, +++ = +3
                while (pos_ < smi_.size() && smi_[pos_] == sign) {
                    chg++;
                    pos_++;
                }
            }
            if (chg > 15)
                throw std::invalid_argument("Charge in bracket atom exceeds the OpenSMILES range of -15 to +15");
            a.charge = (sign == '+') ? chg : -chg;
        }

        // Atom class (optional :nn)
        if (pos_ < smi_.size() && smi_[pos_] == ':') {
            pos_++;
            const int atomClass = parser_detail::readDecimal(smi_, pos_);
            a.atomClass = atomClass < 0 ? 0 : atomClass;
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
                const int slot = static_cast<int>(neighborOrder_[atomIdx].size());
                neighborOrder_[atomIdx].push_back(-1);
                ringOpenings_[ringNum] = std::make_tuple(atomIdx, ringBondOrder, ringBondStereo, slot);
            } else {
                // Close ring
                auto [openAtom, openBo, openStereo, openSlot] = it->second;
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

                // Closing directions point back towards the opening atom.
                if (ringBondStereo != 0) st = 3 - ringBondStereo;
                else if (openStereo != 0) st = openStereo;

                if (bo == 0) {
                    // Default: aromatic if both aromatic, else single
                    if (atoms_[openAtom].isAromatic && atoms_[atomIdx].isAromatic)
                        bo = 4;
                    else
                        bo = 1;
                }

                bonds_.push_back({openAtom, atomIdx, bo, st});
                neighborOrder_[openAtom][openSlot] = atomIdx;
                neighborOrder_[atomIdx].push_back(openAtom);
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
                               bool lenient = false,
                               const std::vector<std::vector<int>>& neighborOrder = {}) {
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

    // Double bond stereo tracking: / and \ are read from the single bonds, and
    // dbStereoConf in MolGraph holds the configuration by the atom pair of each C=C bond.

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

    // Ring digits count at their first appearance in SMILES stereo ordering.
    for (int atom = 0; atom < n && atom < static_cast<int>(neighborOrder.size()); ++atom) {
        if (atoms[atom].chirality == 0) continue;
        std::vector<int> orderedNeighbors, orderedBonds;
        for (int neighbor : neighborOrder[atom]) {
            const auto found = std::find(neighbors[atom].begin(), neighbors[atom].end(), neighbor);
            if (found == neighbors[atom].end()) continue;
            orderedNeighbors.push_back(neighbor);
            orderedBonds.push_back(bondOrders[atom][found - neighbors[atom].begin()]);
        }
        if (orderedNeighbors.size() == neighbors[atom].size()) {
            neighbors[atom] = std::move(orderedNeighbors);
            bondOrders[atom] = std::move(orderedBonds);
        }
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
            if (bondOrders[i][k] == 4 && !inRing) bondOrders[i][k] = 1;
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
        // A ring closure can connect a new component to an earlier atom.
        // Normalise its H position to the graph's ligand-order convention.
        if (chirality[i] != 0 && neighbors[i].size() == 3
            && a.hasIncomingBond != (neighbors[i][0] < i)) chirality[i] = 3 - chirality[i];

        // Implicit H: bracket atoms have explicit hcount, organic atoms compute from valence
        if (a.hcount == -1) {
            // Organic subset atom: implicit H from bondOrderSum, which counts an aromatic
            // bond once; computeImplicitH adds the pi contribution of an aromatic atom.
            hydrogenCounts[i] = computeImplicitH(
                a.atomicNum, a.isAromatic, bondOrderSum[i], a.charge);
        } else {
            hydrogenCounts[i] = a.hcount;
        }
        // A tag on an atom that is not a tetrahedral centre is dropped, as Java parsing does.
        if (chirality[i] != 0 && !stereo_detail::tetrahedralCentre(neighbors[i].size(), hydrogenCounts[i], a.atomicNum))
            chirality[i] = 0;
        // Target properties use int counts, including explicit bonds and H atoms.
        const int64_t totalValence = int64_t(hydrogenCounts[i]) + bondOrderSum[i];
        const int64_t connectivity = int64_t(hydrogenCounts[i]) + neighbors[i].size();
        if (totalValence + (a.isAromatic ? 1 : 0) > std::numeric_limits<int>::max()
            || connectivity > std::numeric_limits<int>::max())
            throw std::invalid_argument("Hydrogen count exceeds the supported atom property range");
    }

    // Collect double bond configurations without a dense intermediate matrix.
    // For / and \ bonds on C=C systems, determine E/Z
    // dbStereoConf[i][j] = 1 for Z (cis), 2 for E (trans) on the C=C bond
    std::vector<std::tuple<int, int, int>> doubleBondConfigurations;
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
        // Find C=C bonds and determine configuration from / \ on adjacent single bonds
        for (auto& b : bonds) {
            if (b.order == 2) {
                int a1 = b.from, a2 = b.to;
                int firstSide = 0, secondSide = 0;
                int firstLigand = -1, secondLigand = -1;
                for (const auto& annotation : stereoAnnotations) {
                    const int direction = annotation.stereo == 1 ? 1 : -1;
                    if (annotation.from == a1 && annotation.to != a2) {
                        firstSide = direction; firstLigand = annotation.to;
                    } else if (annotation.to == a1 && annotation.from != a2) {
                        firstSide = -direction; firstLigand = annotation.from;
                    }
                    if (annotation.from == a2 && annotation.to != a1) {
                        secondSide = direction; secondLigand = annotation.to;
                    } else if (annotation.to == a2 && annotation.from != a1) {
                        secondSide = -direction; secondLigand = annotation.from;
                    }
                }
                if (firstSide != 0 && secondSide != 0) {
                    bool together = firstSide == secondSide;
                    // Store the relation against the first non-double-bond ligand.
                    // This reference is also used by CIP assignment and the writer.
                    for (int neighbor : neighbors[a1]) {
                        if (neighbor == a2) continue;
                        if (neighbor != firstLigand) together = !together;
                        break;
                    }
                    for (int neighbor : neighbors[a2]) {
                        if (neighbor == a1) continue;
                        if (neighbor != secondLigand) together = !together;
                        break;
                    }
                    doubleBondConfigurations.emplace_back(a1, a2, together ? 1 : 2);
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

    if (!doubleBondConfigurations.empty()) {
        builder.doubleBondStereoSparse(std::move(doubleBondConfigurations));
    }

    return builder.build();
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

// Canonical SMARTS writer. Uses the same Morgan-based canonical labelling and
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

    std::vector<int> implicitH(n, 0);
    for (int atom = 0; atom < n; ++atom)
        implicitH[atom] = hydrogenCountForWrite(g, atom);

    std::ostringstream result;
    const auto components = connectedComponents(g);
    for (size_t ci = 0; ci < components.size(); ci++) {
        if (ci > 0) result << '.';
        result << writeCanonicalComponent(g, canon, components[ci],
            [&](int atom, int, const std::vector<int>&, const std::vector<RingDigit>&) {
                return writeAtomSMARTS(g, atom, implicitH[atom], inclArom, inclCharge, inclRing, inclHCount, inclIsotope);
            },
            [&](int, const RingDigit& digit, bool) { return bondSymbolSMARTS(digit.bo); },
            [&](int atom, int nb) { return bondSymbolSMARTS(g.bondOrder(atom, nb)); });
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
 * @return        A fully initialised MolGraph with Morgan ranks,
 *                canonical labels, bit-parallel adjacency, etc.
 * @throws std::invalid_argument on malformed SMILES, or when the molecule has more
 *         than MAX_ATOMS atoms (parser_utils.hpp).
 */
inline MolGraph parseSMILES(const std::string& smiles) {
    detail::SmilesParser parser(smiles);
    parser.parse();
    if (static_cast<int>(parser.atoms().size()) > MAX_ATOMS)
        throw std::invalid_argument("Molecule exceeds MAX_ATOMS limit ("
            + std::to_string(MAX_ATOMS) + ")");
    return detail::buildMolGraph(parser.atoms(), parser.bonds(), false, parser.neighborOrder());
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
    return detail::buildMolGraph(parser.atoms(), parser.bonds(), opts.lenient, parser.neighborOrder());
}

/**
 * Generate a canonical SMILES string from a MolGraph.
 *
 * Uses the canonical labelling computed by MolGraph::initDerivedFields()
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
 * labelling drives the DFS traversal in the same way as toSMILES(), so two
 * isomorphic graphs always produce the same SMARTS string.
 *
 * Atom primitives use the [#Z] form (atomic number) as the base, with
 * optional aromaticity, charge, H count, ring membership, and isotope
 * primitives controlled by SmartsWriteOptions. Bond symbols are always
 * written explicitly (-, =, #, :, ~) so there is no ambiguity.
 *
 * Typical use, writing the scaffold of an MCS result:
 * @code
 * auto common = smsd::MCS::find(query).match(target);
 * std::vector<int> atoms;
 * for (const auto& [qi, ti] : common.matches.front().atomMapping) atoms.push_back(qi);
 * auto sub = smsd::extractSubgraph(query, atoms);
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

#endif // SMSD_SMILES_PARSER_HPP

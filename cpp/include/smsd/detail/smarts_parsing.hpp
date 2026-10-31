/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMARTS parsing: the recursive-descent parser of atoms, bonds, branches, ring
 * closures, component groups and recursive SMARTS.
 *
 * Part of smarts_parser.hpp.
 */
#pragma once

#include "smsd/parser_utils.hpp"
#include "smsd/detail/smarts_evaluation.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {
namespace smarts_detail {

// ============================================================================
// SMARTS parser
// ============================================================================

class SmartsParser {
    const std::string& smarts_;
    int pos_ = 0;
    std::vector<std::vector<int>> neighborOrder_;
    std::vector<bool> hasIncomingBond_;
    /** The component group of each atom, -1 outside one, and the group being read. */
    std::vector<int> componentGroup_;
    int group_ = -1;

    void recordAtomOrder(int atom, int parent) {
        componentGroup_.push_back(group_);
        neighborOrder_.emplace_back();
        hasIncomingBond_.push_back(parent >= 0);
        if (parent >= 0) {
            neighborOrder_[parent].push_back(atom);
            neighborOrder_[atom].push_back(parent);
        }
    }

    char peek() const {
        if (pos_ >= static_cast<int>(smarts_.size())) return '\0';
        return smarts_[pos_];
    }

    char advance() {
        if (pos_ >= static_cast<int>(smarts_.size())) return '\0';
        return smarts_[pos_++];
    }

    bool match(char c) {
        if (peek() == c) { pos_++; return true; }
        return false;
    }

    int parseDigits() {
        return parser_detail::readDecimal(smarts_, pos_);
    }

    // Try to parse a range query {min-max}, {min-}, or {-max}.
    // Returns true if a range was parsed; fills min/max (-1 = unbounded).
    bool tryParseRange(int& rmin, int& rmax) {
        if (peek() != '{') return false;
        int saved = pos_;
        advance(); // consume '{'
        rmin = -1;
        rmax = -1;
        int d = parseDigits();
        if (d >= 0) rmin = d;
        if (peek() == '-') {
            advance(); // consume '-'
            d = parseDigits();
            if (d >= 0) rmax = d;
        } else {
            // {n} without dash means exact match: min == max
            if (rmin >= 0) rmax = rmin;
        }
        if (peek() == '}') {
            advance(); // consume '}'
            return true;
        }
        // Not a valid range, rewind
        pos_ = saved;
        return false;
    }

    // After parsing a numeric primitive, check for optional range {n-m}.
    // If range found, set expr fields accordingly. Otherwise keep exact value.
    void applyOptionalRange(AtomExpr& expr) {
        int rmin, rmax;
        if (tryParseRange(rmin, rmax)) {
            expr.rangeQuery = true;
            expr.rangeMin = rmin;
            expr.rangeMax = rmax;
            expr.hasValue = true;
        }
    }

    // Check if current position starts a two-letter element (e.g., Al, Br, Cl)
    bool isTwoLetterElement(char upper) const {
        if (pos_ + 1 >= static_cast<int>(smarts_.size())) return false;
        char next = smarts_[pos_ + 1];
        if (!parser_detail::isAsciiLower(next)) return false;
        std::string twoChar{upper, next};
        return elementToAtomicNum(twoChar) >= 0;
    }

    // Parse an element symbol (uppercase or lowercase, 1 or 2 chars). As outside brackets, an
    // uppercase symbol is aliphatic and a lowercase one aromatic; only #n matches either form.
    AtomExpr parseElementSymbol() {
        char c = peek();
        if (parser_detail::isAsciiUpper(c)) {
            advance();
            std::string sym(1, c);
            if (pos_ < static_cast<int>(smarts_.size()) && parser_detail::isAsciiLower(peek())) {
                std::string twoChar = sym + peek();
                if (elementToAtomicNum(twoChar) >= 0) {
                    sym = twoChar;
                    advance();
                }
            }
            int z = elementToAtomicNum(sym);
            if (z < 0) throw std::invalid_argument("Unknown element in SMARTS: " + sym);
            return AtomExpr::makeElement(z, false);
        }
        if (parser_detail::isAsciiLower(c)) {
            advance();
            std::string sym(1, c);
            if (pos_ < static_cast<int>(smarts_.size()) && parser_detail::isAsciiLower(peek())) {
                std::string twoChar = sym + peek();
                if (elementToAtomicNum(twoChar) >= 0) {
                    sym = twoChar;
                    advance();
                }
            }
            int z = elementToAtomicNum(sym);
            if (z < 0) throw std::invalid_argument("Unknown aromatic element in SMARTS: " + sym);
            return AtomExpr::makeElement(z, true);
        }
        throw std::invalid_argument(
            std::string("Expected element symbol in SMARTS at '") + c + "'");
    }

    // Parse a single atom primitive inside brackets
    AtomExpr parsePrimitive() {
        char c = peek();

        // # followed by digits = atomic number
        if (c == '#') {
            advance();
            int num = parseDigits();
            if (num < 0) throw std::invalid_argument("Expected number after # in SMARTS");
            return AtomExpr::makePrim(AtomPrimType::ATOMIC_NUM, num, true);
        }

        // * = any atom
        if (c == '*') {
            advance();
            return AtomExpr::makePrim(AtomPrimType::TRUE_);
        }

        // A = aliphatic (but Al, Ag, As, Au, Ac, Am, Ar, At = elements)
        if (c == 'A') {
            if (isTwoLetterElement('A')) return parseElementSymbol();
            advance();
            return AtomExpr::makePrim(AtomPrimType::ALIPHATIC);
        }

        // a = aromatic atom primitive
        if (c == 'a') {
            advance();
            return AtomExpr::makePrim(AtomPrimType::AROMATIC);
        }

        // X or Q = total connectivity (degree + implicit H; but Xe = xenon)
        // Q is Ehrlich-Rarey/Abolmaali synonym for X (e.g. CQ1H2 = CX1H2)
        if (c == 'X' || c == 'Q') {
            if (c == 'X' && isTwoLetterElement('X')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;
            auto expr = AtomExpr::makePrim(AtomPrimType::TOTAL_CONNECTIVITY, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // D = degree (but Dy, Db, Ds = elements)
        if (c == 'D') {
            if (isTwoLetterElement('D')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // D without number = D1
            auto expr = AtomExpr::makePrim(AtomPrimType::DEGREE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // d = heavy degree (non-hydrogen degree) — inside brackets only
        // Note: 'd' could be deuterium in some contexts but we treat it as
        // a primitive inside bracket atoms (per OpenSMARTS convention).
        if (c == 'd') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // d without number = d1
            auto expr = AtomExpr::makePrim(AtomPrimType::HEAVY_DEGREE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // v = valence
        if (c == 'v') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;
            auto expr = AtomExpr::makePrim(AtomPrimType::VALENCE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // R = ring membership (but Rb, Re, Rf, Rg, Rh, Rn, Ru = elements)
        if (c == 'R') {
            if (isTwoLetterElement('R')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::IN_RING, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::IN_RING, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // r = ring size
        if (c == 'r') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::RING_SIZE, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::RING_SIZE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // x = ring connectivity
        if (c == 'x') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::RING_CONNECTIVITY, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::RING_CONNECTIVITY, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // z = heteroatom neighbour count
        if (c == 'z') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                // z without number = at least 1 heteroatom neighbour
                auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // Z = aliphatic heteroatom neighbour count (but Zn, Zr = elements)
        // Only treat as primitive inside brackets when followed by digit, '{', or
        // at end/delimiter position (not followed by lowercase letter forming element).
        if (c == 'Z') {
            if (isTwoLetterElement('Z')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::ALIPHATIC_HETERO_NEIGHBORS, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::ALIPHATIC_HETERO_NEIGHBORS, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // ^ = hybridisation (^0=S, ^1=SP, ^2=SP2, ^3=SP3, ^4=SP3D, ^5=SP3D2)
        if (c == '^') {
            advance();
            int num = parseDigits();
            if (num < 0 || num > 5) {
                throw std::invalid_argument("Expected digit 0-5 after ^ in SMARTS hybridisation");
            }
            return AtomExpr::makePrim(AtomPrimType::HYBRIDIZATION, num, true);
        }

        // + = positive charge
        if (c == '+') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                num = 1;
                while (peek() == '+') { advance(); num++; }
            }
            return AtomExpr::makePrim(AtomPrimType::CHARGE, num, true);
        }

        // - = negative charge (inside brackets only, not a bond)
        if (c == '-') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                num = 1;
                while (peek() == '-') { advance(); num++; }
            }
            return AtomExpr::makePrim(AtomPrimType::CHARGE, -num, true);
        }

        // H = hydrogen count (but Hf, Hg, Ho, Hs = elements)
        if (c == 'H') {
            if (pos_ + 1 < static_cast<int>(smarts_.size())) {
                char next = smarts_[pos_ + 1];
                if (next == 'f' || next == 'g' || next == 'o' || next == 's') {
                    return parseElementSymbol();
                }
            }
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // H without number = H1
            auto expr = AtomExpr::makePrim(AtomPrimType::HCOUNT, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // @ = anticlockwise, @@ = clockwise (tags 1 and 2, as SMILES stores them); a trailing ?
        // (@? or @@?) means that configuration or none, as Daylight SMARTS defines it
        if (c == '@') {
            advance();
            int tag = 1; // @, anticlockwise
            if (pos_ < static_cast<int>(smarts_.size()) && peek() == '@') {
                advance(); // @@, clockwise
                tag = 2;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::CHIRALITY, tag, true);
            if (pos_ < static_cast<int>(smarts_.size()) && peek() == '?') {
                advance(); // @? or @@?
                expr.orUnspecified = true;
            }
            return expr;
        }

        // :n = atom class / map number
        if (c == ':') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 0;
            return AtomExpr::makePrim(AtomPrimType::ATOM_CLASS, num, true);
        }

        // $ = recursive SMARTS
        if (c == '$') {
            advance();
            if (peek() != '(') throw std::invalid_argument("Expected '(' after '$' in recursive SMARTS");
            advance(); // consume '('
            int depth = 1;
            int start = pos_;
            while (pos_ < static_cast<int>(smarts_.size()) && depth > 0) {
                if (smarts_[pos_] == '(') depth++;
                else if (smarts_[pos_] == ')') depth--;
                if (depth > 0) pos_++;
            }
            if (depth != 0) throw std::invalid_argument("Unmatched '(' in recursive SMARTS");
            std::string subSmarts = smarts_.substr(start, pos_ - start);
            advance(); // consume ')'

            AtomExpr expr;
            expr.op = ExprOp::PRIM;
            expr.primType = AtomPrimType::RECURSIVE;
            expr.recursiveSmarts = subSmarts;
            return expr;
        }

        // Element symbol (uppercase or lowercase)
        if (parser_detail::isAsciiLetter(c)) {
            return parseElementSymbol();
        }

        // Isotope: digits at the start of a bracket atom
        if (parser_detail::isAsciiDigit(c)) {
            int num = parseDigits();
            return AtomExpr::makePrim(AtomPrimType::MASS, num, true);
        }

        throw std::invalid_argument(
            std::string("Unexpected character in SMARTS atom: '") + c + "'");
    }

    // Parse NOT-level: !prim
    AtomExpr parseNotExpr() {
        bool negated = false;
        while (peek() == '!') {
            advance();
            negated = !negated;
        }
        AtomExpr expression = parsePrimitive();
        return negated ? AtomExpr::makeNot(std::move(expression)) : std::move(expression);
    }

    // Parse AND-high: prim & prim & ... (implicit AND between consecutive primitives)
    // Operands are collected only when an operator follows the first one.
    AtomExpr parseAndHigh() {
        AtomExpr first = parseNotExpr();
        std::vector<AtomExpr> parts;
        while (pos_ < static_cast<int>(smarts_.size())) {
            char c = peek();
            if (c == ',' || c == ';' || c == ']' || c == '\0') break;
            // Explicit '&' or implicit AND between adjacent primitives
            if (c == '&') advance();
            if (parts.empty()) {
                parts.reserve(4);
                parts.push_back(std::move(first));
            }
            parts.push_back(parseNotExpr());
        }
        if (parts.empty()) return first;
        return AtomExpr::makeAnd(std::move(parts));
    }

    // Parse OR: ... , ...
    AtomExpr parseOr() {
        AtomExpr first = parseAndHigh();
        if (peek() != ',') return first;
        std::vector<AtomExpr> parts;
        parts.reserve(4);
        parts.push_back(std::move(first));
        while (peek() == ',') {
            advance();
            parts.push_back(parseAndHigh());
        }
        return AtomExpr::makeOr(std::move(parts));
    }

    // Parse AND-low (semicolon): ... ; ...
    AtomExpr parseAndLow() {
        AtomExpr first = parseOr();
        if (peek() != ';') return first;
        std::vector<AtomExpr> parts;
        parts.reserve(4);
        parts.push_back(std::move(first));
        while (peek() == ';') {
            advance();
            parts.push_back(parseOr());
        }
        return AtomExpr::makeAnd(std::move(parts));
    }

    /**
     * A bracket atom that is a hydrogen atom, as Daylight SMARTS reads one: an optional isotope, H,
     * an optional charge and an optional atom class, and nothing else, as in [H], [2H], [H+] and
     * [H:1]. Elsewhere H is a hydrogen count, as in [CH] or [H1]. Returns false, consuming
     * nothing, for any other bracket atom.
     */
    bool parseHydrogenAtom(AtomExpr& expression) {
        const int size = static_cast<int>(smarts_.size());
        int at = pos_, isotope = -1, charge = 0;
        if (at < size && parser_detail::isAsciiDigit(smarts_[at])) {
            isotope = 0;
            while (at < size && parser_detail::isAsciiDigit(smarts_[at])) {
                if (isotope > 100000) return false;
                isotope = isotope * 10 + (smarts_[at++] - '0');
            }
        }
        if (at >= size || smarts_[at] != 'H') return false;
        ++at;
        if (at < size && (smarts_[at] == '+' || smarts_[at] == '-')) {
            const char sign = smarts_[at];
            int magnitude = 0;
            while (at < size && smarts_[at] == sign) { ++magnitude; ++at; }
            if (magnitude == 1 && at < size && parser_detail::isAsciiDigit(smarts_[at])) {
                magnitude = 0;
                while (at < size && parser_detail::isAsciiDigit(smarts_[at])) {
                    if (magnitude > 100000) return false;
                    magnitude = magnitude * 10 + (smarts_[at++] - '0');
                }
            }
            charge = sign == '+' ? magnitude : -magnitude;
        }
        int atomClass = -1;
        if (at < size && smarts_[at] == ':') {
            ++at;
            if (at >= size || !parser_detail::isAsciiDigit(smarts_[at])) return false;
            atomClass = 0;
            while (at < size && parser_detail::isAsciiDigit(smarts_[at])) {
                if (atomClass > 100000) return false;
                atomClass = atomClass * 10 + (smarts_[at++] - '0');
            }
        }
        if (at >= size || smarts_[at] != ']') return false;
        std::vector<AtomExpr> parts;
        parts.push_back(AtomExpr::makeElement(1, false));
        if (isotope >= 0) parts.push_back(AtomExpr::makePrim(AtomPrimType::MASS, isotope, true));
        if (charge != 0) parts.push_back(AtomExpr::makePrim(AtomPrimType::CHARGE, charge, true));
        if (atomClass >= 0) parts.push_back(AtomExpr::makePrim(AtomPrimType::ATOM_CLASS, atomClass, true));
        expression = AtomExpr::makeAnd(std::move(parts));
        pos_ = at;
        return true;
    }

    // Parse bracket atom expression: everything between [ and ]
    AtomExpr parseBracketAtomExpr() {
        AtomExpr hydrogen;
        if (parseHydrogenAtom(hydrogen)) return hydrogen;
        return parseAndLow();
    }

    // Bond precedence: negation, high AND, OR, low AND.
    BondExpr parseBondExpr() {
        BondExpr first = parseBondOr();
        if (peek() != ';') return first;
        std::vector<BondExpr> parts;
        parts.push_back(std::move(first));
        while (peek() == ';') {
            advance();
            parts.push_back(parseBondOr());
        }
        return BondExpr::makeAnd(std::move(parts));
    }

    BondExpr parseBondOr() {
        BondExpr first = parseBondAndHigh();
        if (peek() != ',') return first;
        std::vector<BondExpr> parts;
        parts.push_back(std::move(first));
        while (peek() == ',') {
            advance();
            parts.push_back(parseBondAndHigh());
        }
        return BondExpr::makeOr(std::move(parts));
    }

    BondExpr parseBondAndHigh() {
        BondExpr first = parseSingleBond();
        if (peek() != '&' && !isBondChar(peek())) return first;
        std::vector<BondExpr> parts;
        parts.push_back(std::move(first));
        while (peek() == '&' || isBondChar(peek())) {
            if (peek() == '&') advance();
            parts.push_back(parseSingleBond());
        }
        return BondExpr::makeAnd(std::move(parts));
    }

    BondExpr parseSingleBond() {
        bool negated = false;
        while (peek() == '!') {
            advance();
            negated = !negated;
        }
        BondPrimType type;
        switch (peek()) {
            case '-': type = BondPrimType::SINGLE; break;
            case '=': type = BondPrimType::DOUBLE; break;
            case '#': type = BondPrimType::TRIPLE; break;
            case ':': type = BondPrimType::AROMATIC; break;
            case '~': type = BondPrimType::ANY; break;
            case '@': type = BondPrimType::RING; break;
            case '/': type = BondPrimType::WEDGE; break;
            case '\\': type = BondPrimType::DASH; break;
            default: throw std::invalid_argument("Expected bond primitive in SMARTS");
        }
        advance();
        // /? and \? also accept a double bond without a configuration.
        if ((type == BondPrimType::WEDGE || type == BondPrimType::DASH)
            && pos_ < static_cast<int>(smarts_.size()) && peek() == '?') {
            advance();
            type = type == BondPrimType::WEDGE ? BondPrimType::WEDGE_OR_UNSPECIFIED : BondPrimType::DASH_OR_UNSPECIFIED;
        }
        BondExpr expression = BondExpr::makePrim(type);
        return negated ? BondExpr::makeNot(std::move(expression)) : std::move(expression);
    }

    bool isBondChar(char c) const {
        return c == '-' || c == '=' || c == '#' || c == ':' || c == '~'
            || c == '@' || c == '!' || c == '/' || c == '\\';
    }

    bool isAtomStart(char c) const {
        return c == '[' || c == '*' || parser_detail::isAsciiLetter(c);
    }

public:
    SmartsParser(const std::string& smarts) : smarts_(smarts) {}

    /**
     * The component group of each atom: atoms within one top-level ( ) group must match atoms of
     * one target component, and atoms of different groups atoms of different components, as in
     * (C.C) and (C).(C). Empty when the pattern has no groups.
     */
    std::vector<int> componentGroups() const {
        for (int group : componentGroup_) if (group >= 0) return componentGroup_;
        return {};
    }

    /** Normalise parsed tetrahedral primitives to the query adjacency order. */
    void normaliseChirality(std::vector<SmartsAtom>& atoms,
                            const std::vector<std::vector<int>>& neighbors) const {
        std::function<void(AtomExpr&)> flip = [&](AtomExpr& expression) {
            if (expression.op == ExprOp::PRIM && expression.primType == AtomPrimType::CHIRALITY
                && expression.intVal != 0) expression.intVal = 3 - expression.intVal;
            for (auto& child : expression.children) flip(child);
        };
        for (int atom = 0; atom < static_cast<int>(atoms.size()); ++atom) {
            if (!hasChiralityPrimitive(atoms[atom].expr)) continue;
            std::array<int, 4> storedOrder{};
            if (!stereo_detail::tetrahedralLigands(neighbors[atom], atom, storedOrder)) continue;
            auto parsedOrder = neighborOrder_[atom];
            if (parsedOrder.size() == 3)
                parsedOrder.insert(parsedOrder.begin() + (hasIncomingBond_[atom] ? 1 : 0), -1);
            int inversions = 0;
            for (size_t first = 0; first < parsedOrder.size(); ++first) {
                const auto firstPosition = std::find(storedOrder.begin(), storedOrder.end(), parsedOrder[first]);
                for (size_t second = first + 1; second < parsedOrder.size(); ++second) {
                    const auto secondPosition = std::find(storedOrder.begin(), storedOrder.end(), parsedOrder[second]);
                    inversions += firstPosition > secondPosition;
                }
            }
            if (inversions % 2 != 0) flip(atoms[atom].expr);
        }
    }

    // Parse the full SMARTS string into atoms and bonds
    void parse(std::vector<SmartsAtom>& atoms, std::vector<SmartsBond>& bonds) {
        if (smarts_.empty()) return;

        std::vector<int> branchStack;   // stack of parent atom indices
        int currentAtom = -1, groups = 0;
        struct RingClosure { int atom; BondExpr expression; int neighborSlot; };
        std::unordered_map<int, RingClosure> ringClosures;
        BondExpr pendingBond;
        bool hasPendingBond = false;

        while (pos_ < static_cast<int>(smarts_.size())) {
            char c = peek();

            if (c == '(') {
                advance();
                // A ( that opens the pattern or follows a dot at the top level groups components.
                if (currentAtom < 0 && branchStack.empty() && !hasPendingBond) group_ = groups++;
                branchStack.push_back(currentAtom);
                hasPendingBond = false;
                continue;
            }

            if (c == ')') {
                advance();
                if (branchStack.empty()) {
                    throw std::invalid_argument("Unmatched ')' in SMARTS");
                }
                currentAtom = branchStack.back();
                branchStack.pop_back();
                if (branchStack.empty() && currentAtom < 0) group_ = -1;
                hasPendingBond = false;
                continue;
            }

            // Bond symbol?
            if (isBondChar(c)) {
                pendingBond = parseBondExpr();
                hasPendingBond = true;
                continue;
            }

            // Ring closure digit or %nn
            if (parser_detail::isAsciiDigit(c) || c == '%') {
                int ringNum;
                if (c == '%') {
                    advance();
                    ringNum = parseDigits();
                    if (ringNum < 0) throw std::invalid_argument("Expected digits after % in ring closure");
                } else {
                    ringNum = advance() - '0';
                }

                if (currentAtom < 0) {
                    throw std::invalid_argument("Ring closure before any atom");
                }

                auto it = ringClosures.find(ringNum);
                if (it != ringClosures.end()) {
                    // Close the ring
                    SmartsBond bond;
                    const int openAtom = it->second.atom;
                    bond.from = hasPendingBond ? currentAtom : openAtom;
                    bond.to = hasPendingBond ? openAtom : currentAtom;
                    neighborOrder_[openAtom][it->second.neighborSlot] = currentAtom;
                    neighborOrder_[currentAtom].push_back(openAtom);
                    // Use the bond from open if it had one, or the current pending bond; a bond
                    // written at both ends must be the same expression at each, as Java reads it.
                    // A direction at the closing digit is read from the closing atom, so one
                    // written at both ends is rejected, as Java rejects it.
                    const bool openWritten = it->second.expression.primType != BondPrimType::DEFAULT_
                        || it->second.expression.op != ExprOp::PRIM;
                    if (hasPendingBond && openWritten
                        && (bondExprHasDirection(pendingBond) || bondExprHasDirection(it->second.expression)))
                        throw std::invalid_argument("Ring closure " + std::to_string(ringNum)
                                                    + " has a bond direction written at both ends");
                    if (hasPendingBond && openWritten && !sameBondExpr(pendingBond, it->second.expression))
                        throw std::invalid_argument("Ring closure " + std::to_string(ringNum)
                                                    + " has different bond expressions at its two ends");
                    if (hasPendingBond) {
                        bond.expr = std::move(pendingBond);
                        hasPendingBond = false;
                    } else if (openWritten) {
                        bond.expr = std::move(it->second.expression);
                    } else {
                        bond.expr = BondExpr::makePrim(BondPrimType::DEFAULT_);
                    }
                    bonds.push_back(std::move(bond));
                    ringClosures.erase(it);
                } else {
                    // Open the ring
                    BondExpr ringBondExpr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    const int slot = static_cast<int>(neighborOrder_[currentAtom].size());
                    neighborOrder_[currentAtom].push_back(-1);
                    ringClosures[ringNum] = {currentAtom, std::move(ringBondExpr), slot};
                }
                continue;
            }

            // Atom: [ bracket atom ] or organic shorthand
            if (c == '[') {
                advance(); // consume [

                SmartsAtom atom;
                atom.expr = parseBracketAtomExpr();

                if (peek() != ']') {
                    throw std::invalid_argument(
                        "Expected ']' in bracket atom at pos " + std::to_string(pos_));
                }
                advance(); // consume ]

                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            // Organic shorthand atoms: B C N O P S F Cl Br I * and aromatic b c n o p s
            if (c == '*') {
                advance();
                SmartsAtom atom;
                atom.expr = AtomExpr::makePrim(AtomPrimType::TRUE_);
                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            if (parser_detail::isAsciiLetter(c)) {
                // Outside brackets, a and A are always atom wildcards.
                if (c == 'a' || c == 'A') {
                    advance();
                    SmartsAtom atom;
                    atom.isAromatic = (c == 'a');
                    atom.expr = AtomExpr::makePrim(c == 'a' ? AtomPrimType::AROMATIC
                                                           : AtomPrimType::ALIPHATIC, 0, true);
                    int newIdx = static_cast<int>(atoms.size());
                    atoms.push_back(std::move(atom));
                    recordAtomOrder(newIdx, currentAtom);
                    if (currentAtom >= 0) {
                        SmartsBond bond;
                        bond.from = currentAtom; bond.to = newIdx;
                        bond.expr = hasPendingBond ? std::move(pendingBond)
                            : BondExpr::makePrim(BondPrimType::DEFAULT_);
                        hasPendingBond = false;
                        bonds.push_back(std::move(bond));
                    }
                    currentAtom = newIdx;
                    continue;
                }
                bool aromatic = parser_detail::isAsciiLower(c);
                std::string sym(1, c);
                advance();

                if ((sym == "C" && peek() == 'l') || (sym == "B" && peek() == 'r')) {
                    sym += advance();
                }

                int z = elementToAtomicNum(sym);
                if (z < 0) {
                    throw std::invalid_argument("Unknown element in SMARTS: " + sym);
                }

                SmartsAtom atom;
                atom.isAromatic = aromatic;
                atom.expr = AtomExpr::makeElement(z, aromatic);

                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            // Skip whitespace
            if (parser_detail::isAsciiSpace(c)) {
                advance();
                continue;
            }

            // Dot (disconnected)
            if (c == '.') {
                advance();
                currentAtom = -1;
                hasPendingBond = false;
                continue;
            }

            throw std::invalid_argument(
                std::string("Unexpected character in SMARTS: '") + c
                + "' at position " + std::to_string(pos_));
        }

        if (!branchStack.empty()) {
            throw std::invalid_argument("Unmatched '(' in SMARTS");
        }
        if (!ringClosures.empty()) {
            throw std::invalid_argument("Unclosed ring in SMARTS");
        }
        if (hasPendingBond) {
            throw std::invalid_argument("Trailing bond in SMARTS");
        }
    }
};

} // namespace smarts_detail
} // namespace smsd

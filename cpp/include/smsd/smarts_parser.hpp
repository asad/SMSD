/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only SMARTS pattern parser for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 *
 * Supports:
 *   - Atom primitives: #n, A, a, D<n>, v<n>, R, r<n>, +/-, H<n>, *, elements
 *   - Extended primitives: z<n> (heteroatom neighbours), Z<n> (aliphatic
 *     heteroatom neighbours), d<n> (heavy degree), ^<n> (hybridisation)
 *   - Range queries: {n-m}, {n-}, {-m} on numeric primitives
 *   - Logical operators: , (OR), & (AND high), ; (AND low), ! (NOT)
 *   - Bond primitives: - (single), = (double), # (triple), : (aromatic),
 *                       ~ (any), @ (ring bond), / (wedge), \ (dash)
 *   - Ring closures, branches ()
 *   - Recursive SMARTS: $(...)
 *   - Named predicate registry (SmartsPredicateRegistry) with 10 built-ins
 *
 * The parts are in detail/: smarts_expressions.hpp, smarts_properties.hpp,
 * smarts_evaluation.hpp, smarts_parsing.hpp and smarts_query.hpp. This header
 * holds the named predicate registry and parseSMARTS().
 */
#pragma once
#ifndef SMSD_SMARTS_PARSER_HPP
#define SMSD_SMARTS_PARSER_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/parser_utils.hpp"
#include "smsd/stereo.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "smsd/detail/smarts_expressions.hpp"
#include "smsd/detail/smarts_properties.hpp"
#include "smsd/detail/smarts_evaluation.hpp"
#include "smsd/detail/smarts_parsing.hpp"
#include "smsd/detail/smarts_query.hpp"

namespace smsd {

// ============================================================================
// Named predicate registry (mirrors Java Standardiser.expandNamedPredicates)
// ============================================================================

class SmartsPredicateRegistry {
    std::unordered_map<std::string, std::string> predicates_;

public:
    SmartsPredicateRegistry() = default;

    /// Register a named predicate: name -> SMARTS definition.
    void registerPredicate(const std::string& name, const std::string& smarts) {
        predicates_[name] = smarts;
    }

    /// Expand all `$name` tokens (inside brackets) with `$(definition)`.
    /// Scans for `$` followed by an identifier (alphanumeric + underscore);
    /// if the identifier is a registered predicate, replaces `$name` with
    /// `$(definition)`.  Tokens already of the form `$(...)` are left alone.
    std::string expandPredicates(const std::string& input) const {
        std::string result;
        result.reserve(input.size());
        size_t i = 0;
        while (i < input.size()) {
            if (input[i] == '$' && i + 1 < input.size() && input[i + 1] != '(') {
                // Collect the identifier after '$'
                size_t start = i + 1;
                size_t end = start;
                while (end < input.size() &&
                       (parser_detail::isAsciiAlnum(input[end]) || input[end] == '_')) {
                    ++end;
                }
                std::string name = input.substr(start, end - start);
                auto it = predicates_.find(name);
                if (it != predicates_.end()) {
                    result += "$(";
                    result += it->second;
                    result += ')';
                    i = end;
                } else {
                    result += input[i];
                    ++i;
                }
            } else {
                result += input[i];
                ++i;
            }
        }
        return result;
    }

    /// Create a registry pre-loaded with the 10 built-in predicates
    /// (same set as Java Standardiser).
    static SmartsPredicateRegistry builtinRegistry() {
        SmartsPredicateRegistry reg;
        reg.registerPredicate("isAmideN",
            "[$([NX3;H2,H1;!$(NC=O)]),$([NX3;$(NC=O)])]");
        reg.registerPredicate("isCarboxylC",
            "[CX3](=O)[OX2H1,OX1-]");
        reg.registerPredicate("isSulfonamideN",
            "[NX3;$(NS(=O)=O)]");
        reg.registerPredicate("isEster",
            "[#6][CX3](=O)[OX2][#6]");
        reg.registerPredicate("isKetone",
            "[#6][CX3](=O)[#6]");
        reg.registerPredicate("isHalogen",
            "[F,Cl,Br,I]");
        reg.registerPredicate("isPositive",
            "[*+]");
        reg.registerPredicate("isNegative",
            "[*-]");
        reg.registerPredicate("isHBDonor",
            "[N!H0,O!H0,S!H0]");
        reg.registerPredicate("isHBAcceptor",
            "[N,O,S;H0;!+0]");
        return reg;
    }
};

// ============================================================================
// parseSMARTS -- top-level API
// ============================================================================

namespace smarts_detail {

/** Nodes of an atom expression as MAX_SMARTS_EXPRESSION_NODES counts them; a recursive SMARTS counts one. */
inline size_t expressionNodes(const AtomExpr& expression) {
    size_t nodes = expression.op == ExprOp::PRIM || expression.op == ExprOp::NOT
        ? 1 : expression.children.size() - 1;
    for (const auto& child : expression.children) nodes += expressionNodes(child);
    return nodes;
}

inline size_t expressionNodes(const BondExpr& expression) {
    size_t nodes = expression.op == ExprOp::PRIM || expression.op == ExprOp::NOT
        ? 1 : expression.children.size() - 1;
    for (const auto& child : expression.children) nodes += expressionNodes(child);
    return nodes;
}

/** parseSMARTS() of one pattern level, adding its atoms and nodes, and those of its recursive SMARTS, to the totals. */
inline SmartsQuery parsePattern(const std::string& smarts, int maxRecursionDepth, size_t& atoms, size_t& nodes);

} // namespace smarts_detail

/**
 * Parse a SMARTS pattern. A pattern longer than MAX_SMARTS_LENGTH characters, with more than
 * MAX_ATOMS atoms or more than MAX_SMARTS_EXPRESSION_NODES expression nodes, counting its
 * recursive SMARTS, with an atom or bond expression of more than
 * MAX_SMARTS_NODES_PER_EXPRESSION nodes, or with recursive SMARTS nested deeper than
 * maxRecursionDepth, is rejected with std::invalid_argument, as is malformed SMARTS.
 */
inline SmartsQuery parseSMARTS(const std::string& smarts, int maxRecursionDepth) {
    if (smarts.size() > MAX_SMARTS_LENGTH)
        throw std::invalid_argument("SMARTS exceeds the maximum length of "
            + std::to_string(MAX_SMARTS_LENGTH) + " characters");
    size_t atoms = 0, nodes = 0;
    return smarts_detail::parsePattern(smarts, maxRecursionDepth, atoms, nodes);
}

inline SmartsQuery smarts_detail::parsePattern(const std::string& smarts, int maxRecursionDepth,
                                               size_t& atoms, size_t& nodes) {
    if (smarts.empty()) return SmartsQuery{};

    // Guard against unbounded recursive SMARTS nesting by tracking actual
    // nested $(...) frames rather than counting raw occurrences.
    int parenDepth = 0;
    std::vector<int> recursiveFrames;
    recursiveFrames.reserve(8);
    for (size_t i = 0; i < smarts.size(); ++i) {
        char c = smarts[i];
        if (c == '(') {
            ++parenDepth;
            if (i > 0 && smarts[i - 1] == '$') {
                recursiveFrames.push_back(parenDepth);
                if (static_cast<int>(recursiveFrames.size()) > maxRecursionDepth) {
                    throw std::invalid_argument(
                        "Recursive SMARTS nesting exceeds maximum depth of "
                        + std::to_string(maxRecursionDepth));
                }
            }
        } else if (c == ')') {
            if (!recursiveFrames.empty() && recursiveFrames.back() == parenDepth)
                recursiveFrames.pop_back();
            if (parenDepth > 0) --parenDepth;
        }
    }

    SmartsQuery query;
    // Reserve for the atoms written at this level and their ring closures.
    size_t atomEstimate = 0, closureDigits = 0;
    int bracketDepth = 0;
    for (char c : smarts) {
        if (c == '[') { if (bracketDepth++ == 0) ++atomEstimate; }
        else if (c == ']') { if (bracketDepth > 0) --bracketDepth; }
        else if (bracketDepth == 0) {
            if (parser_detail::isAsciiLetter(c) || c == '*') ++atomEstimate;
            else if (parser_detail::isAsciiDigit(c)) ++closureDigits;
        }
    }
    query.atoms.reserve(atomEstimate);
    query.bonds.reserve(atomEstimate + closureDigits / 2);
    smarts_detail::SmartsParser parser(smarts);
    parser.parse(query.atoms, query.bonds);
    query.componentGroups = parser.componentGroups();
    atoms += query.atoms.size();
    if (atoms > static_cast<size_t>(MAX_ATOMS))
        throw std::invalid_argument("SMARTS exceeds the maximum of " + std::to_string(MAX_ATOMS) + " atoms");
    auto expression = [&](size_t count) {
        if (count > MAX_SMARTS_NODES_PER_EXPRESSION)
            throw std::invalid_argument("SMARTS atom or bond expression exceeds the maximum of "
                + std::to_string(MAX_SMARTS_NODES_PER_EXPRESSION) + " nodes");
        nodes += count;
        return static_cast<uint32_t>(count);
    };
    query.atomNodes.reserve(query.atoms.size());
    query.bondNodes.reserve(query.bonds.size());
    for (const auto& atom : query.atoms) query.atomNodes.push_back(expression(expressionNodes(atom.expr)));
    for (const auto& bond : query.bonds) query.bondNodes.push_back(expression(expressionNodes(bond.expr)));
    if (nodes > MAX_SMARTS_EXPRESSION_NODES)
        throw std::invalid_argument("SMARTS exceeds the maximum of "
            + std::to_string(MAX_SMARTS_EXPRESSION_NODES) + " expression nodes");
    query.buildAdjacency();
    parser.normaliseChirality(query.atoms, query.neighbors);
    std::function<void(smarts_detail::AtomExpr&)> compileRecursive = [&](smarts_detail::AtomExpr& expression) {
        if (expression.primType == smarts_detail::AtomPrimType::RECURSIVE)
            expression.cachedRecursive = std::make_shared<SmartsQuery>(
                parsePattern(expression.recursiveSmarts, maxRecursionDepth, atoms, nodes));
        for (auto& child : expression.children) compileRecursive(child);
    };
    for (auto& atom : query.atoms) compileRecursive(atom.expr);
    return query;
}

} // namespace smsd

#endif // SMSD_SMARTS_PARSER_HPP

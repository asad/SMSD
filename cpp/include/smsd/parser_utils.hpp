/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace smsd {

// ============================================================================
// Input limits shared by the SMILES and SMARTS parsers. They lie far above any
// real chemistry and keep the memory and work of crafted input bounded; input
// beyond them is rejected with std::invalid_argument (ValueError in Python).
// ============================================================================

/// Maximum atom count of a parsed SMILES molecule or SMARTS pattern; a SMARTS
/// pattern counts the atoms of its recursive SMARTS as well.
constexpr int MAX_ATOMS = 10000;

/// Maximum length of a SMARTS pattern, in characters, after named predicates are expanded.
constexpr size_t MAX_SMARTS_LENGTH = 100000;

/// Maximum number of atom and bond expression nodes in a SMARTS pattern, counting those of
/// its recursive SMARTS: each primitive test and each logical operator joining two operands
/// counts once, and repeated negations count as at most one.
constexpr size_t MAX_SMARTS_EXPRESSION_NODES = 100000;

/// Maximum number of nodes in one atom or bond expression of a SMARTS pattern; a recursive
/// SMARTS counts one node in the expression that holds it.
constexpr size_t MAX_SMARTS_NODES_PER_EXPRESSION = 10000;

namespace parser_detail {

/** Read a non-negative int, returning -1 when no digits are present. */
template<class Position>
inline int readDecimal(const std::string& text, Position& position) {
    int value = 0;
    const Position start = position;
    while (static_cast<size_t>(position) < text.size()
           && text[position] >= '0' && text[position] <= '9') {
        const int digit = text[position] - '0';
        if (value > (std::numeric_limits<int>::max() - digit) / 10)
            throw std::invalid_argument("Numeric field exceeds the supported integer range");
        value = value * 10 + digit;
        ++position;
    }
    return position == start ? -1 : value;
}

} // namespace parser_detail
} // namespace smsd

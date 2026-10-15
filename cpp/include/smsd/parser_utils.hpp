/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#pragma once

#include <limits>
#include <stdexcept>
#include <string>

namespace smsd {
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

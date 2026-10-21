/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Disjoint sets over item indices, stored as a parent array.
 */
#pragma once

#include <vector>

namespace smsd {
namespace detail {

/** The root of the set of `item`, halving the path to it on the way. */
inline int unionFindRoot(std::vector<int>& parent, int item) {
    while (parent[item] != item) item = parent[item] = parent[parent[item]];
    return item;
}

/** Joins the sets of `first` and `second`; the root of `first` is placed below the root of `second`. */
inline void unionFindJoin(std::vector<int>& parent, int first, int second) {
    first = unionFindRoot(parent, first);
    second = unionFindRoot(parent, second);
    if (first != second) parent[first] = second;
}

} // namespace detail
} // namespace smsd

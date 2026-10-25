/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once
#include <algorithm>
#include <cstdint>

namespace smsd {
namespace screening_detail {
/**
 * The bound for a label overlap, the sum over labels of the smaller of the two counts, between
 * graphs of n1 and n2 atoms: the overlap is capped by the smaller graph and gives the similarity
 * overlap / (n1 + n2 - overlap). Every CPU screen scores through this, so their scores agree bit for bit.
 */
inline double overlapSimilarity(int overlap, int n1, int n2) {
    if (n1 == 0 || n2 == 0) return 0.0;
    const int atomOverlap = std::min(overlap, std::min(n1, n2));
    if (atomOverlap == 0) return 0.0;
    const int64_t atomUnion = int64_t(n1) + n2 - atomOverlap;
    return static_cast<double>(atomOverlap) / static_cast<double>(atomUnion);
}
}
}

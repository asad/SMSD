/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * RASCAL — O(V+E) similarity upper-bound screening.
 * Raymond, Gardiner & Willett, JCAMD 2002.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace smsd {

// similarityUpperBound() defined in mcs.hpp

/// Batch RASCAL screening: find all targets with similarity >= threshold.
/// Returns indices of passing targets.
inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold) {
    std::vector<int> hits;
    hits.reserve(targets.size() / 10); // expect ~10% pass
    for (int i = 0; i < static_cast<int>(targets.size()); ++i) {
        if (similarityUpperBound(query, targets[i]) >= threshold) {
            hits.push_back(i);
        }
    }
    return hits;
}

} // namespace smsd

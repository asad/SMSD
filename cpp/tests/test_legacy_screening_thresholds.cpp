/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#ifndef SMSD_TEST_NATIVE_GPU
#define SMSD_ENABLE_METAL
#endif
#include "test_support.hpp"
#include "smsd/gpu.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>

#ifndef SMSD_TEST_NATIVE_GPU
namespace smsd::metal {
bool metalAvailable() noexcept { return true; }
namespace detail {
std::vector<int> batchSimilarityOverlaps(
    const MetalMolecule& query, const std::vector<MetalMolecule>& targets) {
    std::vector<int> overlaps(targets.size());
    for (size_t target = 0; target < targets.size(); ++target)
        for (int label = 0; label < METAL_LABEL_BINS; ++label)
            overlaps[target] += std::min(query.labelCounts[label], targets[target].labelCounts[label]);
    return overlaps;
}
}
std::vector<MetalScreenResult> batchSimilarityScreen(
    const MetalMolecule& query, const std::vector<MetalMolecule>& targets, float threshold) {
    const auto overlaps = detail::batchSimilarityOverlaps(query, targets);
    std::vector<MetalScreenResult> hits;
    for (size_t target = 0; target < targets.size(); ++target) {
        const int64_t atomUnion = int64_t(query.n) + targets[target].n - overlaps[target];
        const float score = query.n == 0 || targets[target].n == 0 || atomUnion <= 0
            ? 0.0f : static_cast<float>(overlaps[target]) / static_cast<float>(atomUnion);
        if (score >= threshold) hits.push_back({static_cast<int>(target), score});
    }
    return hits;
}
}
#endif

int main() {
    try {
#ifdef SMSD_TEST_NATIVE_GPU
        if (!smsd::gpu::isAvailable()) return 77;
#endif
        const auto graph = [](const std::vector<int>& elements) {
            return smsd::MolGraph::Builder().atomCount(static_cast<int>(elements.size()))
                .atomicNumbers(elements).setNeighbors(std::vector<std::vector<int>>(elements.size())).build(false);
        };
        const auto query = graph({6, 7, 7, 7, 7});
        const auto target = graph({6, 8, 8, 8, 8, 8});
        const std::vector<smsd::MolGraph> targets{target, target, query, graph({})};
        const auto scores = smsd::screening_detail::similarityScores(query, targets, smsd::ChemOptions(), 1);
        const double boundary = scores[0];
        smsd_test::require(boundary == 0.1, "Decimal fixture has one common atom out of ten union atoms");
        for (double threshold : {std::nextafter(boundary, 0.0), boundary,
                std::nextafter(boundary, 1.0), 0.0, -0.0, 1.0,
                -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::quiet_NaN()}) {
            const auto expected = smsd::screening_detail::indicesReaching(scores, threshold);
            const auto actual = smsd::gpu::batchSimilarityScreenAuto(query, targets, threshold);
            smsd_test::require(actual.size() == expected.size(), "Threshold membership must follow the CPU double score");
            for (size_t hit = 0; hit < actual.size(); ++hit) {
                smsd_test::require(actual[hit].index == expected[hit], "Tied targets retain input order");
                smsd_test::require(std::memcmp(&actual[hit].score, &scores[expected[hit]], sizeof(double)) == 0,
                                   "Host scores must agree bit for bit");
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

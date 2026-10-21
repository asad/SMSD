/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/batch.hpp"
#include "smsd/gpu.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/smiles_parser.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>

int main() {
    int failures = 0;
    const auto require = [&](bool condition, const char* message) {
        if (!condition) { std::cerr << message << '\n'; ++failures; }
    };
    const auto query = smsd::parseSMILES("c1ccccc1CC(=O)O");
    const auto target = smsd::parseSMILES("C1CCCCC1CC(=O)O");
    const std::vector<smsd::MolGraph> targets{target};
    smsd::ChemOptions chemistry;
    smsd::MCSOptions options;
    options.timeoutMs = 1000;
    require(std::abs(smsd::similarityUpperBound(query, target) - 1.0) < 1e-8,
            "Default screening must allow aromatic atoms to match aliphatic atoms");
    require(smsd::screenTargets(query, targets, 0.5) == std::vector<int>{0},
            "Default screening must retain the complete common graph");
    const auto mappings = smsd_test::screenedCommon(query, targets, chemistry, options, 0.5, 1);
    require(mappings.size() == 1 && mappings[0].second.size() == 10,
            "Screened MCS must retain all ten common atoms");
    const auto sizes = smsd_test::screenedCommonSizes(query, targets, chemistry, options, 0.5, 1);
    require(sizes == std::vector<std::pair<int, int>>{{0, 10}},
            "Screened MCS sizes must retain all ten common atoms");
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    require(smsd_test::screenedCommonSizes(query, targets, chemistry, options, 0.5, 1).empty(),
            "Strict aromaticity must still constrain the screening bound");
    chemistry = smsd::ChemOptions();
    chemistry.matchAtomType = false;
    require(smsd_test::screenedCommonSizes(smsd::parseSMILES("C"), {smsd::parseSMILES("N")},
                chemistry, options, 0.5, 1) == std::vector<std::pair<int, int>>{{0, 1}},
            "Untyped matching must not be rejected by element screening");
    const auto ring = smsd::parseSMILES("C1CCCCC1");
    const auto chain = smsd::parseSMILES("CCCCCC");
    require(std::abs(smsd::similarityUpperBound(ring, chain) - 1.0) < 1e-8,
            "Default screening must allow ring atoms to match chain atoms");
    require(smsd::screenTargets(smsd::MolGraph(), targets, 0.5).empty(),
            "An empty query must not pass a positive screening threshold");

    const auto gpuHits = smsd::gpu::batchSimilarityScreenAuto(query, targets, 0.5);
    require(gpuHits.size() == 1 && std::abs(gpuHits[0].score - 1.0) < 1e-6,
            "Automatic screening must use the same default bound as CPU screening");
    const auto overlapHits = smsd::gpu::batchSimilarityScreenAuto(
        smsd::parseSMILES("CCCO"), {smsd::parseSMILES("CCCN")}, 0.0);
    require(overlapHits.size() == 1 && std::abs(overlapHits[0].score - 0.6) < 1e-6,
            "Automatic screening scores must use the Tanimoto denominator");
    const auto cisplatin = smsd::parseSMILES("N[Pt](N)(Cl)Cl");
    const auto heavyHits = smsd::gpu::batchSimilarityScreenAuto(cisplatin, {cisplatin}, 0.9);
    require(heavyHits.size() == 1 && std::abs(heavyHits[0].score - 1.0) < 1e-6,
            "Heavy elements must remain in screening histograms");
    chemistry = smsd::ChemOptions();
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    require(std::abs(smsd::similarityUpperBound(query, target, chemistry) - 0.25) < 1e-8,
            "The explicit chemistry bound must preserve strict aromaticity");
    require(smsd::gpu::batchSimilarityScreenAuto(query, targets, chemistry, 0.5).empty(),
            "Automatic screening must preserve strict aromaticity");
    chemistry = smsd::ChemOptions();
    chemistry.ringMatchesRingOnly = true;
    require(smsd::screenTargets(ring, {chain}, 0.5, chemistry).empty()
            && smsd::gpu::batchSimilarityScreenAuto(ring, {chain}, chemistry, 0.5).empty(),
            "Explicit ring-only screening must reject chain atoms");
    chemistry = smsd::ChemOptions();
    chemistry.matchAtomType = false;
    const auto untypedHits = smsd::gpu::batchSimilarityScreenAuto(
        smsd::parseSMILES("C"), {smsd::parseSMILES("N")}, chemistry, 0.5);
    require(untypedHits.size() == 1 && std::abs(untypedHits[0].score - 1.0) < 1e-6,
            "Automatic untyped screening must allow different elements");
    require(smsd::gpu::batchSimilarityScreenAuto(smsd::MolGraph(), targets, 0.5).empty(),
            "Automatic screening must reject empty queries at positive thresholds");

    // Labels outside the dense table, at or above its limit or below zero, take the hashed path
    // for the query or for one target, and leave the scratch table clear for the next target.
    {
        auto labelled = [](std::vector<int> elements) {
            auto graph = smsd::MolGraph::Builder().atomCount(static_cast<int>(elements.size()))
                .atomicNumbers(std::vector<int>(elements.size(), 6))
                .setNeighbors(std::vector<std::vector<int>>(elements.size())).build(false);
            graph.atomicNum = std::move(elements);
            return graph;
        };
        const int beyond = smsd::screening_detail::QueryLabelCounts::denseLimit / 4;
        const std::vector<smsd::MolGraph> graphs{labelled({6, 6, 8}), labelled({beyond, 6}), labelled({-1, 6, 6}),
                                                 labelled({6, 8, beyond, 6, -1}), labelled({8, 8, 6}), labelled({beyond + 1, beyond})};
        const smsd::ChemOptions typed;
        std::vector<int> scratch(smsd::screening_detail::QueryLabelCounts::denseLimit, 0);
        bool consistent = true;
        for (const auto& first : graphs) {
            const smsd::screening_detail::QueryLabelCounts counts(first, typed);
            for (const auto& second : graphs) {
                consistent = consistent
                    && counts.upperBound(second, scratch) == smsd::screening_detail::labelFrequencyUpperBound(first, second, typed)
                    && std::all_of(scratch.begin(), scratch.end(), [](int count) { return count == 0; })
                    && std::abs(counts.similarity(second, scratch) - smsd::similarityUpperBound(first, second, typed)) < 1e-12;
            }
        }
        require(consistent, "Dense query label counts must agree with the hashed bound for every label");
    }
#ifdef SMSD_ENABLE_METAL
    const auto packed = smsd::metal::toMetalMolecule(cisplatin.n, cisplatin.label.data());
    require(std::accumulate(std::begin(packed.labelCounts), std::end(packed.labelCounts), 0) == cisplatin.n,
            "Metal packing must preserve every atom, including platinum");
    std::cout << (smsd::metal::metalAvailable() ? "Metal device available\n" : "Metal device unavailable; CPU fallback exercised\n");
#endif
    return failures == 0 ? 0 : 1;
}

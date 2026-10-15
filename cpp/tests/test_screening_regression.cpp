#include "search_test_helpers.hpp"
/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/batch.hpp"
#include "smsd/gpu.hpp"
#include "smsd/rascal.hpp"
#include "smsd/smiles_parser.hpp"
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

    const auto gpuHits = smsd::gpu::batchRascalScreenAuto(query, targets, 0.5);
    require(gpuHits.size() == 1 && std::abs(gpuHits[0].score - 1.0) < 1e-6,
            "Automatic screening must use the same default bound as CPU screening");
    const auto overlapHits = smsd::gpu::batchRascalScreenAuto(
        smsd::parseSMILES("CCCO"), {smsd::parseSMILES("CCCN")}, 0.0);
    require(overlapHits.size() == 1 && std::abs(overlapHits[0].score - 0.6) < 1e-6,
            "Automatic screening scores must use the Tanimoto denominator");
    const auto cisplatin = smsd::parseSMILES("N[Pt](N)(Cl)Cl");
    const auto heavyHits = smsd::gpu::batchRascalScreenAuto(cisplatin, {cisplatin}, 0.9);
    require(heavyHits.size() == 1 && std::abs(heavyHits[0].score - 1.0) < 1e-6,
            "Heavy elements must remain in screening histograms");
    chemistry = smsd::ChemOptions();
    chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    require(std::abs(smsd::similarityUpperBound(query, target, chemistry) - 0.25) < 1e-8,
            "The explicit chemistry bound must preserve strict aromaticity");
    require(smsd::gpu::batchRascalScreenAuto(query, targets, chemistry, 0.5).empty(),
            "Automatic screening must preserve strict aromaticity");
    chemistry = smsd::ChemOptions();
    chemistry.ringMatchesRingOnly = true;
    require(smsd::screenTargets(ring, {chain}, 0.5, chemistry).empty()
            && smsd::gpu::batchRascalScreenAuto(ring, {chain}, chemistry, 0.5).empty(),
            "Explicit ring-only screening must reject chain atoms");
    chemistry = smsd::ChemOptions();
    chemistry.matchAtomType = false;
    const auto untypedHits = smsd::gpu::batchRascalScreenAuto(
        smsd::parseSMILES("C"), {smsd::parseSMILES("N")}, chemistry, 0.5);
    require(untypedHits.size() == 1 && std::abs(untypedHits[0].score - 1.0) < 1e-6,
            "Automatic untyped screening must allow different elements");
    require(smsd::gpu::batchRascalScreenAuto(smsd::MolGraph(), targets, 0.5).empty(),
            "Automatic screening must reject empty queries at positive thresholds");
#ifdef SMSD_ENABLE_METAL
    const auto packed = smsd::metal::toMetalMolecule(cisplatin.n, cisplatin.label.data());
    require(std::accumulate(std::begin(packed.labelCounts), std::end(packed.labelCounts), 0) == cisplatin.n,
            "Metal packing must preserve every atom, including platinum");
    std::cout << (smsd::metal::metalAvailable() ? "Metal device available\n" : "Metal device unavailable; CPU fallback exercised\n");
#endif
    return failures == 0 ? 0 : 1;
}

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
#include <cstring>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

// The recorded label histograms must not cost the corpus its nothrow move, which std::vector relies on.
static_assert(std::is_nothrow_move_constructible<smsd::batch::TargetCorpus>::value
              && std::is_nothrow_move_assignable<smsd::batch::TargetCorpus>::value,
              "TargetCorpus must stay nothrow movable");

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

    const auto labelled = [](std::vector<int> elements) {
        auto graph = smsd::MolGraph::Builder().atomCount(static_cast<int>(elements.size()))
            .atomicNumbers(std::vector<int>(elements.size(), 6))
            .setNeighbors(std::vector<std::vector<int>>(elements.size())).build(false);
        graph.atomicNum = std::move(elements);
        return graph;
    };
    // Labels outside the dense table, at or above its limit or below zero, take the hashed path
    // for the query or for one target, and leave the scratch table clear for the next target.
    {
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

    // The screen gives the same scores on every thread count, a corpus screens as the list of its
    // targets does, and a failing target propagates the failure of the first such target.
    {
        const std::vector<std::string> smiles{"c1ccccc1", "CCO", "CC(=O)Oc1ccccc1C(=O)O", "C1CCCCC1", "c1ccncc1",
                                              "CCN(CC)CC", "O=C(O)c1ccccc1O", "C[N+](C)(C)C", "ClC(Cl)Cl", "c1ccc2ccccc2c1"};
        std::vector<smsd::MolGraph> library;
        for (size_t i = 0; i < 3000; ++i) library.push_back(smsd::parseSMILES(smiles[(i * 7) % smiles.size()]));
        const auto aspirin = smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O");
        const smsd::ChemOptions defaults;
        const auto serial = smsd::screening_detail::similarityScores(aspirin, library, defaults, 1);
        bool same = serial.size() == library.size();
        for (size_t i = 0; same && i < library.size(); ++i)
            same = serial[i] == smsd::similarityUpperBound(aspirin, library[i], defaults);
        for (int threads : {2, 4, 0}) {
            same = same && smsd::screening_detail::similarityScores(aspirin, library, defaults, threads) == serial;
            for (double threshold : {0.0, 0.3, 0.5, 1.0})
                same = same && smsd::screenTargets(aspirin, library, threshold, defaults, threads)
                                   == smsd::screenTargets(aspirin, library, threshold, defaults);
        }
        require(same, "The similarity screen must not depend on the thread count");
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(library);
        require(corpus.screenTargets(aspirin, 0.3, defaults, 0) == smsd::screenTargets(aspirin, library, 0.3, defaults),
                "A corpus must screen as the list of its targets does");
        struct Failing {
            const std::vector<smsd::MolGraph>& graphs;
            size_t size() const { return graphs.size(); }
            const smsd::MolGraph& operator[](size_t index) const {
                if (index == 1700 || index == 2900) throw std::runtime_error(std::to_string(index));
                return graphs[index];
            }
        };
        for (int threads : {1, 2, 0}) {
            std::string thrown;
            try {
                smsd::screening_detail::similarityScores(aspirin, Failing{library}, defaults, threads);
            } catch (const std::runtime_error& error) {
                thrown = error.what();
            }
            require(thrown == "1700", "A failing target must propagate the first failure on every thread count");
        }
    }

    // A corpus scores from the label histograms it records under each label rule, bit for bit as
    // the graphs score; targets with a label outside the dense table, the empty query and a query
    // with such a label are scored on the graphs, and added targets replace the histograms.
    {
        const std::vector<std::string> smiles{"c1ccccc1", "CCO", "CC(=O)Oc1ccccc1C(=O)O", "C1CCCCC1", "c1ccncc1",
                                              "CCN(CC)CC", "O=C(O)c1ccccc1O", "C[N+](C)(C)C", "ClC(Cl)Cl", "c1ccc2ccccc2c1"};
        const int beyond = smsd::screening_detail::QueryLabelCounts::denseLimit / 4;
        std::vector<smsd::MolGraph> library;
        for (size_t i = 0; i < 3000; ++i) library.push_back(smsd::parseSMILES(smiles[(i * 7) % smiles.size()]));
        library[5] = labelled({beyond, 6, 6});
        library[1200] = labelled({-1, 8});
        const std::vector<smsd::MolGraph> queries{smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O"), smsd::MolGraph(),
                                                  smsd::parseSMILES("c1ccccc1"), labelled({beyond, 6, 8})};
        std::vector<smsd::ChemOptions> rules(5);
        rules[0].matchAtomType = false;
        rules[2].aromaticityMode = rules[4].aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
        rules[3].ringMatchesRingOnly = rules[4].ringMatchesRingOnly = true;
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(library);
        const auto exact = [&](const std::vector<smsd::MolGraph>& graphs) {
            bool same = true;
            for (const auto& chemistry : rules) {
                const smsd::screening_detail::TargetLabels labels(graphs, chemistry, 0);
                for (const auto& query : queries) {
                    const auto expected = smsd::screening_detail::similarityScores(query, graphs, chemistry, 1);
                    std::vector<double> thresholds(expected.begin(), expected.end());
                    std::sort(thresholds.begin(), thresholds.end());
                    thresholds.erase(std::unique(thresholds.begin(), thresholds.end()), thresholds.end());
                    for (size_t k = 0, distinct = thresholds.size(); k < distinct; ++k)
                        thresholds.push_back(std::nextafter(thresholds[k], 2.0));
                    for (int threads : {1, 2, 0}) {
                        const auto scores = smsd::screening_detail::similarityScores(query, graphs, labels, chemistry, threads);
                        same = same && scores.size() == expected.size()
                            && std::memcmp(scores.data(), expected.data(), scores.size() * sizeof(double)) == 0;
                        // Agreement at every score and the next double above it means equal scores.
                        for (double threshold : thresholds)
                            same = same && corpus.screenTargets(query, threshold, chemistry, threads)
                                               == smsd::screening_detail::indicesReaching(expected, threshold);
                    }
                }
            }
            return same;
        };
        require(exact(library), "Recorded label histograms must give the graph scores bit for bit");
        corpus.addTarget(smsd::parseSMILES("c1ccccc1CCN"));
        library.push_back(smsd::parseSMILES("c1ccccc1CCN"));
        require(exact(library), "Adding a target must replace the recorded label histograms");
        const std::vector<smsd::MolGraph> more{smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O"), labelled({beyond, 6})};
        corpus.addTargets(more);
        library.insert(library.end(), more.begin(), more.end());
        require(exact(library), "Adding targets must replace the recorded label histograms");

        // Concurrent first screens under each rule build every table once and agree with the graphs.
        smsd::batch::TargetCorpus fresh;
        fresh.addTargets(library);
        std::vector<int> agreed(16, 0);
        std::vector<std::thread> screens;
        for (int k = 0; k < 16; ++k)
            screens.emplace_back([&, k] {
                const auto& chemistry = rules[k % rules.size()];
                agreed[k] = fresh.screenTargets(queries[0], 0.3, chemistry, k % 3)
                            == smsd::screenTargets(queries[0], library, 0.3, chemistry);
            });
        for (auto& screen : screens) screen.join();
        require(std::all_of(agreed.begin(), agreed.end(), [](int same) { return same == 1; }),
                "Concurrent first screens must build the label histograms safely");
        const smsd::batch::TargetCorpus copy = fresh;
        require(copy.screenTargets(queries[0], 0.3, rules[1], 0) == smsd::screenTargets(queries[0], library, 0.3, rules[1]),
                "A copied corpus must screen as its targets do");
    }

    // A target whose atom count is negative, which only an edit of n gives, is scored on its graph and
    // takes no room in the histograms.
    {
        smsd::MolGraph malformed = smsd::parseSMILES("CCO");
        malformed.n = -1;
        std::vector<smsd::MolGraph> graphs(600, smsd::parseSMILES("c1ccccc1"));
        graphs[0] = graphs[300] = malformed;
        const smsd::ChemOptions chemistry;
        const auto query = smsd::parseSMILES("c1ccccc1O");
        const auto expected = smsd::screening_detail::similarityScores(query, graphs, chemistry, 1);
        const smsd::screening_detail::TargetLabels labels(graphs, chemistry, 0);
        const auto scores = smsd::screening_detail::similarityScores(query, graphs, labels, chemistry, 0);
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(graphs);
        require(!labels.held(0) && !labels.held(300)
                    && std::memcmp(scores.data(), expected.data(), scores.size() * sizeof(double)) == 0
                    && corpus.screenTargets(query, -1.0, chemistry, 0)
                           == smsd::screening_detail::indicesReaching(expected, -1.0),
                "A target with a negative atom count must be scored on its graph");
    }
#ifdef SMSD_ENABLE_METAL
    const auto packed = smsd::metal::toMetalMolecule(cisplatin.n, cisplatin.label.data());
    require(std::accumulate(std::begin(packed.labelCounts), std::end(packed.labelCounts), 0) == cisplatin.n,
            "Metal packing must preserve every atom, including platinum");
    std::cout << (smsd::metal::metalAvailable() ? "Metal device available\n" : "Metal device unavailable; CPU fallback exercised\n");
#endif
    return failures == 0 ? 0 : 1;
}

#include "search_test_helpers.hpp"
/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/smsd.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/batch.hpp"
#include "smsd/layout.hpp"

#include <atomic>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <thread>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void parallelCalls(int workerCount, const std::function<void(int)>& operation) {
    std::atomic<int> ready{0};
    std::atomic<bool> start{false};
    std::vector<std::exception_ptr> errors(workerCount);
    std::vector<std::thread> workers;
    for (int worker = 0; worker < workerCount; ++worker) {
        workers.emplace_back([&, worker] {
            ++ready;
            while (!start.load()) std::this_thread::yield();
            try { operation(worker); }
            catch (...) { errors[worker] = std::current_exception(); }
        });
    }
    while (ready.load() != workerCount) std::this_thread::yield();
    start = true;
    for (auto& worker : workers) worker.join();
    for (const auto& error : errors) if (error) std::rethrow_exception(error);
}

static void testTemplateFirstUse() {
    std::vector<std::pair<size_t, size_t>> counts(8);
    parallelCalls(8, [&](int worker) {
        counts[worker] = {smsd::detail_layout::getTemplates().size(),
                          smsd::detail_layout::getExtendedTemplates().size()};
    });
    const auto expected = std::make_pair(smsd::detail_layout::getTemplates().size(),
                                        smsd::detail_layout::getExtendedTemplates().size());
    require(expected.first > 0 && expected.second > 0, "scaffold templates must be available");
    for (const auto& count : counts)
        require(count == expected, "concurrent callers observed an incomplete template library");
}

static void testRecursiveQueryReuse() {
    const auto query = smsd::parseSMARTS("[C;$(CO)]");
    std::vector<smsd::MolGraph> targets = {
        smsd::parseSMILES(std::string(60, 'C') + "O"),
        smsd::parseSMILES("OCCO"), smsd::parseSMILES("CCCC")
    };
    for (const auto& target : targets) smsd::batch::detail::prewarmGraph(target);
    parallelCalls(6, [&](int worker) {
        const int targetIndex = worker % 3;
        for (int iteration = 0; iteration < 300; ++iteration)
            require(smsd_test::patternMatches(query, targets[targetIndex]) == (targetIndex != 2),
                    "shared recursive query reused another target's answer");
    });
}

static void checkSharedGraphs(bool initialiseDocumentedCaches) {
    for (int trial = 0; trial < 4; ++trial) {
        const auto query = smsd::parseSMILES("c1ccc(CCO)cc1");
        const auto target = smsd::parseSMILES("c1ccc(CCO)cc1O");
        if (initialiseDocumentedCaches) {
            for (const auto* graph : {&query, &target}) {
                graph->ensureCanonical();
                graph->ensureRingCounts();
                graph->ensureRingSystems();
            }
        }
        parallelCalls(6, [&](int) {
            smsd::ChemOptions options;
            options.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
            for (int iteration = 0; iteration < 10; ++iteration) {
                require(smsd_test::substructureExists(query, target, options),
                        "shared graph lost a valid substructure");
                require(!smsd_test::allSubstructureAtoms(query, target, options).empty(),
                        "shared graph lost valid enumerated mappings");
                require(target.computeRings().size() == 1, "shared graph lost its ring");
                require(target.getPharmacophoreFeatures().size() == static_cast<size_t>(target.n),
                        "shared graph has incomplete atom features");
            }
        });
    }
}

static void testColdGraphReuse() { checkSharedGraphs(false); }
static void testPartlyWarmedGraphReuse() { checkSharedGraphs(true); }

static void testParallelBatches() {
    std::vector<smsd::MolGraph> targets;
    for (int target = 0; target < 12; ++target)
        targets.push_back(smsd::parseSMILES("c1ccc(CCO)cc1O"));
    const std::vector<smsd::MolGraph> queries = {
        smsd::parseSMILES("CCO"), smsd::parseSMILES("CO"),
        smsd::parseSMILES("c1ccccc1"), smsd::parseSMILES("O")
    };
    parallelCalls(4, [&](int worker) {
        const auto matches = smsd::Substructure::find(queries[worker]).matchMany(targets, 1);
        require(matches.size() == targets.size(), "batch lost a shared target");
        for (const auto& matched : matches) require(!matched.matches.empty(), "parallel batches lost a valid match");
    });
}

static void testConcurrentCorpusUpdates() {
    smsd::batch::TargetCorpus corpus;
    corpus.addTarget(smsd::parseSMILES("CCO"));
    corpus.prewarm(1);
    const auto query = smsd::parseSMILES("CO");
    parallelCalls(4, [&](int worker) {
        if (worker == 0) {
            for (int iteration = 0; iteration < 24; ++iteration)
                corpus.addTarget(smsd::parseSMILES("CCO"));
        } else {
            for (int iteration = 0; iteration < 20; ++iteration) {
                const auto matches = smsd::Substructure::find(query).matchMany(corpus.targets(), 1);
                require(!matches.empty() && matches.size() <= 25, "invalid corpus snapshot size");
                for (const auto& matched : matches) require(!matched.matches.empty(), "corpus update corrupted a target");
            }
        }
    });
    require(corpus.size() == 25, "corpus lost an added target");
}

static void testPropertyCacheInvalidation() {
    const auto query = smsd::parseSMARTS("[C;^2]");
    auto target = smsd::parseSMILES("CC");
    require(!smsd_test::patternMatches(query, target), "ethane is not sp2 hybridised");
    target.setBondOrder(0, 1, 2);
    require(smsd_test::patternMatches(query, target), "bond mutation left stale SMARTS properties");
    const auto base = smsd::parseSMILES("CC");
    for (int order : {1, 2, 1, 2}) {
        auto copy = base;
        if (order != 1) copy.setBondOrder(0, 1, order);
        require(smsd_test::patternMatches(query, copy) == (order == 2), "copied graph reused stale SMARTS properties");
    }
}

static void testAromaticCacheInvalidation() {
    auto query = smsd::parseSMILES("c1ccccc1");
    auto target = smsd::parseSMILES("Cc1ccccc1");
    smsd::batch::detail::prewarmGraph(query);
    require(query.kekulize() && target.kekulize(), "ordinary aromatic rings should Kekulise");
    smsd::ChemOptions options;
    options.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    require(!smsd_test::allSubstructureAtoms(query, target, options).empty(),
            "aromaticity mutation left stale neighbourhood labels");
}

int main(int argc, char** argv) {
    int failures = 0;
    const std::pair<const char*, void (*)()> tests[] = {
        {"templates", testTemplateFirstUse}, {"recursive", testRecursiveQueryReuse},
        {"cold", testColdGraphReuse}, {"partly-warmed", testPartlyWarmedGraphReuse},
        {"batches", testParallelBatches}, {"corpus", testConcurrentCorpusUpdates},
        {"invalidation", testPropertyCacheInvalidation}, {"aromatic-invalidation", testAromaticCacheInvalidation}
    };
    for (const auto& test : tests) {
        if (argc > 1 && argv[1] != std::string(test.first)) continue;
        try { test.second(); std::cout << "PASS " << test.first << '\n'; }
        catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}

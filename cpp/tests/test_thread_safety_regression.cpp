/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/smsd.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/batch.hpp"
#include "smsd/layout.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>

using smsd_test::require;

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

// Every index below the lowest failing one is processed exactly once, and that failure is the one
// reported, whatever the order, the team size and the timing.
static void testDispatcherLowestFailure() {
    std::mt19937 random(8);
    for (size_t size : {size_t(1), size_t(7), size_t(300), size_t(5000)})
        for (size_t workers : {size_t(1), size_t(2), size_t(4), size_t(8)})
            for (int arrangement = 0; arrangement < 3; ++arrangement) {
                std::vector<size_t> order(size);
                std::iota(order.begin(), order.end(), size_t(0));
                if (arrangement == 1) std::reverse(order.begin(), order.end());
                if (arrangement == 2) std::shuffle(order.begin(), order.end(), random);
                for (const std::set<size_t>& failing : {std::set<size_t>{}, std::set<size_t>{size / 2},
                                                         std::set<size_t>{size - 1, size / 3, size / 2}}) {
                    std::vector<std::atomic<int>> visits(size);
                    for (auto& count : visits) count.store(0);
                    smsd::search_detail::BatchDispatcher dispatcher(size, workers, arrangement == 0 ? nullptr : &order);
                    smsd::search_detail::runTeam(workers, [&] {
                        dispatcher.work([&](size_t index) {
                            visits[index].fetch_add(1);
                            if (failing.count(index)) throw std::runtime_error(std::to_string(index));
                        });
                    });
                    const size_t lowest = failing.empty() ? size : *failing.begin();
                    std::string reported;
                    try { dispatcher.rethrowFailure(); }
                    catch (const std::runtime_error& error) { reported = error.what(); }
                    require(reported == (failing.empty() ? std::string() : std::to_string(lowest)),
                            "the dispatcher must report the lowest failing index");
                    for (size_t index = 0; index < size; ++index) {
                        require(visits[index].load() <= 1, "an index was processed twice");
                        require(index > lowest || visits[index].load() == 1, "an index below the failure was skipped");
                    }
                }
            }
}

// Short batches start their largest targets first, ties in input order, with either sort.
static void testLargestTargetsFirst() {
    struct Sized { int n; };
    std::mt19937 random(9);
    for (int largest : {5, 1000}) {
        std::vector<Sized> targets(40);
        for (auto& target : targets) target.n = static_cast<int>(random() % (largest + 1));
        std::vector<size_t> expected(targets.size());
        std::iota(expected.begin(), expected.end(), size_t(0));
        std::stable_sort(expected.begin(), expected.end(), [&](size_t first, size_t second) {
            return targets[first].n > targets[second].n;
        });
        require(smsd::search_detail::largestTargetsFirst(targets) == expected, "largest targets must start first");
    }
}

// Batches run on standard threads, short ones largest first and long ones, beyond 64 targets per
// worker, in input order; both keep the serial results and report the lowest failing target.
static void testBatchThreads() {
    struct Target { int n; size_t id; bool fail; };
    struct Query {
        std::mutex* mutex;
        std::set<std::thread::id>* threads;
        smsd::SearchResult match(const Target& target) const {
            // About 20 microseconds of work, so that a team starts after the serial start.
            const auto until = std::chrono::steady_clock::now() + std::chrono::microseconds(20);
            while (std::chrono::steady_clock::now() < until) {}
            {
                const std::lock_guard<std::mutex> lock(*mutex);
                threads->insert(std::this_thread::get_id());
            }
            if (target.fail) throw std::runtime_error(std::to_string(target.id));
            smsd::SearchResult result;
            result.nodesVisited = target.id;
            result.score = target.n;
            return result;
        }
    };
    std::mt19937 random(10);
    for (size_t size : {size_t(40), size_t(1500)})
        for (int workers : {1, 2, 4, 8}) {
            std::vector<Target> targets(size);
            for (size_t index = 0; index < size; ++index) targets[index] = {static_cast<int>(random() % 50), index, false};
            std::mutex mutex;
            std::set<std::thread::id> threads;
            const Query query{&mutex, &threads};
            const auto results = smsd::search_detail::matchMany(query, targets, workers);
            require(results.size() == size, "a batch lost a target");
            for (size_t index = 0; index < size; ++index)
                require(results[index].nodesVisited == index && results[index].score == targets[index].n,
                        "a batch result is out of place");
            // A long batch lasts long enough for its helpers to take part however threads are scheduled.
            const size_t team = std::min({size, size_t(workers), smsd::search_detail::maximumWorkers()});
            if (team == 1) require(threads.size() == 1, "a one-worker batch must stay on the calling thread");
            else if (size > 64 * team) require(threads.size() > 1, "a long batch did not use its team");
            for (size_t failing : {size / 2, size / 5, size - 1}) targets[failing].fail = true;
            std::string reported;
            try { smsd::search_detail::matchMany(query, targets, workers); }
            catch (const std::runtime_error& error) { reported = error.what(); }
            require(reported == std::to_string(size / 5), "a batch must report its lowest failing target");
        }
}

// The worker count is positive, steady and within what the platform lets this process use, and
// real searches give the same results on any number of threads, with matches built from rows.
static void testWorkerCount() {
    using namespace smsd::search_detail;
    const size_t workers = maximumWorkers();
    require(workers >= 1 && workers == maximumWorkers(), "the worker count must be positive and steady");
    const size_t available = availableProcessors(), quota = quotaProcessors();
    require(available == 0 || workers <= available, "the worker count must stay within the usable processors");
    require(quota == 0 || workers <= quota, "the worker count must stay within the CPU quota");
    std::vector<smsd::MolGraph> targets;
    const char* smiles[] = {"c1ccc(CCO)cc1O", "CC(C)Cc1ccc(cc1)C(C)C(=O)O", "c1ccc2ccccc2c1", "OCCO", "CC(=O)Oc1ccccc1C(=O)O"};
    for (int copy = 0; copy < 60; ++copy)
        for (const char* text : smiles) targets.push_back(smsd::parseSMILES(text));
    const auto query = smsd::Substructure::find(smsd::parseSMILES("c1ccccc1"));
    std::vector<smsd::SearchResult> serial;
    for (const auto& target : targets) serial.push_back(query.matchAll(target));
    struct AllMatches {
        const smsd::Substructure& query;
        smsd::SearchResult match(const smsd::MolGraph& target) const { return query.matchAll(target); }
    };
    for (int threads : {1, 2, 3, 8, 0}) {
        const auto parallel = matchMany(AllMatches{query}, targets, threads);
        require(parallel.size() == serial.size(), "a batch lost a target");
        for (size_t index = 0; index < serial.size(); ++index) {
            const auto& first = serial[index];
            const auto& second = parallel[index];
            require(first.status == second.status && first.optimal == second.optimal && first.complete == second.complete
                    && first.score == second.score && first.nodesVisited == second.nodesVisited
                    && first.matches.size() == second.matches.size(), "batch results depend on the thread count");
            for (size_t match = 0; match < first.matches.size(); ++match)
                require(first.matches[match].atomMapping == second.matches[match].atomMapping
                        && first.matches[match].bondMapping == second.matches[match].bondMapping
                        && first.matches[match].fragments == second.matches[match].fragments,
                        "batch matches depend on the thread count");
        }
    }
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
    return smsd_test::runGroups(argc, argv, {
        {"templates", testTemplateFirstUse}, {"recursive", testRecursiveQueryReuse},
        {"cold", testColdGraphReuse}, {"partly-warmed", testPartlyWarmedGraphReuse},
        {"batches", testParallelBatches}, {"corpus", testConcurrentCorpusUpdates},
        {"dispatcher", testDispatcherLowestFailure}, {"largest-first", testLargestTargetsFirst},
        {"batch-threads", testBatchThreads}, {"workers", testWorkerCount},
        {"invalidation", testPropertyCacheInvalidation}, {"aromatic-invalidation", testAromaticCacheInvalidation}});
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"

#include <iostream>
#include <set>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

using Key = std::pair<std::map<int, int>, std::vector<std::array<int, 4>>>;

static std::set<Key> keys(const smsd::SearchResult& result) {
    require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
            "Search must complete with proof");
    std::set<Key> found;
    for (const auto& match : result.matches)
        require(found.insert({match.atomMapping, match.bondMapping}).second, "No duplicate witnesses");
    return found;
}

static void refuses(const smsd::MolGraph& target) {
    const auto query = smsd::parseSMILES("CC");
    for (int operation = 0; operation < 4; ++operation) {
        bool refused = false;
        try {
            if (operation == 0) (void)smsd::Substructure::find(query).match(target);
            if (operation == 1) (void)smsd::Substructure::smarts("[#6]").match(target);
            if (operation == 2) (void)smsd::MCS::find(query).match(target);
            if (operation == 3) {
                smsd::MCSOptions options;
                options.objective = smsd::MCSObjective::MCES;
                (void)smsd::MCS::find(query).withOptions(options).match(target);
            }
        } catch (const std::invalid_argument&) {
            refused = true;
        }
        require(refused, "Incoherent topology must be clearly refused before search");
    }
}

static smsd::MolGraph sparseGraph() {
    auto graph = smsd::parseSMILES("CCC");
    graph.useDense = false;
    graph.bondPacked.clear();
    graph.sparseBondProps = {
        {smsd::MolGraph::bondKey(0, 1), {1, 0, 0}},
        {smsd::MolGraph::bondKey(1, 2), {1, 0, 0}}};
    return graph;
}

static void invalidDenseTopology() {
    auto hidden = smsd::parseSMILES("CCC");
    hidden.bondPacked[0][2] = hidden.bondPacked[2][0] = 1;
    refuses(hidden);

    auto removed = smsd::parseSMILES("CCC");
    removed.bondPacked[0][1] = removed.bondPacked[1][0] = 0;
    refuses(removed);

    auto directed = smsd::parseSMILES("CCC");
    directed.bondPacked[0][1] = 0;
    refuses(directed);

    auto missing = smsd::parseSMILES("CCC");
    missing.neighbors[0].clear();
    refuses(missing);
}

static void invalidSparseTopology() {
    auto outside = sparseGraph();
    outside.sparseBondProps[smsd::MolGraph::bondKey(0, 3)] = {1, 0, 0};
    refuses(outside);

    auto negative = sparseGraph();
    negative.sparseBondProps[-1] = {1, 0, 0};
    refuses(negative);

    auto order = sparseGraph();
    order.sparseBondProps[smsd::MolGraph::bondKey(0, 1)][0] = -1;
    refuses(order);

    auto zero = sparseGraph();
    zero.sparseBondProps[smsd::MolGraph::bondKey(0, 1)][0] = 0;
    refuses(zero);

    auto inactive = sparseGraph();
    inactive.sparseBondProps[smsd::MolGraph::bondKey(0, 2)] = {0, 0, 0};
    inactive.sparseBondProps[smsd::MolGraph::bondKey(0, 0)] = {0, 0, 0};
    require(keys(smsd::Substructure::find("CC").matchAll(inactive)).size() == 4,
            "Inactive sparse entries are not bonds");
}

static void optionalAdjacencyFallback() {
    const auto query = smsd::parseSMILES("CCC");
    const auto ordinary = smsd::parseSMILES("CCCC");
    smsd::MCSOptions options;
    options.timeoutMs = 10000;
    options.connectivity = smsd::MCSConnectivity::Fragments;
    const auto prepared = smsd::MCS::find(query).withOptions(options);
    const auto expected = keys(prepared.matchAll(ordinary, 1000));
    for (int mode = 0; mode < 4; ++mode) {
        auto stale = ordinary;
        if (mode == 0) stale.adjLong.clear();
        if (mode == 1) stale.adjLong[1][0] = 0;
        if (mode == 2) stale.adjLong[0][0] |= uint64_t(1) << 3;
        if (mode == 3) stale.adjLong[0][0] |= uint64_t(1) << 63;
        require(keys(prepared.matchAll(stale, 1000)) == expected,
                "Optional stale adjacency must fall back to actual coherent neighbours");
    }
}

static void supportedLegacyStorage() {
    auto legacy = smsd::parseSMILES("CCC");
    legacy.neighbors[1].push_back(0);
    legacy.neighbors[2].push_back(2);
    legacy.bondPacked[1][0] = 2;
    const auto result = smsd::Substructure::find("CC").matchAll(legacy);
    require(result.status == smsd::SearchStatus::Completed && result.complete,
            "Legacy listings and nonzero directional bond properties remain supported");
}

static void capturedInputsValidateBeforeSearch() {
    for (int malformed = 0; malformed < 2; ++malformed) {
        auto invalid = smsd::parseSMILES("CCC");
        if (malformed == 0) invalid.atomicNum.resize(2);
        else invalid.n = -1;
        const std::vector<smsd::MolGraph> targets{smsd::parseSMILES("CC"), invalid};
        const auto captured = smsd::batch::detail::TargetStorage::capture(targets);
        require(captured.size() == targets.size() && captured[1].n == invalid.n,
                "Owned snapshots retain graphs independently of graph-search admission");
        for (int addition = 0; addition < 3; ++addition) {
            smsd::batch::TargetCorpus corpus;
            if (addition == 0) {
                for (const auto& target : targets) corpus.addTarget(target);
            } else if (addition == 1) corpus.addTargets(targets);
            else smsd::batch::detail::CorpusAccess::addCapturedTargets(corpus, captured);
            require(corpus.size() == targets.size(), "Storage retains every original target slot");
            for (int operation = 0; operation < 6; ++operation) {
                bool refused = false;
                try {
                    if (operation == 0) (void)smsd::Substructure::find("CC").matchCorpus(corpus, false, 1);
                    if (operation == 1) (void)smsd::Substructure::smarts("[#6]").matchCorpus(corpus, false, 1);
                    if (operation == 2) (void)smsd::MCS::find("CC").matchCorpus(corpus, 2, false, 1);
                    if (operation == 3) {
                        smsd::MCSOptions options;
                        options.objective = smsd::MCSObjective::MCES;
                        (void)smsd::MCS::find("CC").withOptions(options).matchCorpus(corpus, 2, false, 1);
                    }
                    if (operation == 4) {
                        smsd::SubstructureOptions options;
                        options.chemistry.matchAtomType = false;
                        (void)smsd::Substructure::find("CC").withOptions(options).matchCorpus(corpus, false, 1);
                    }
                    if (operation == 5) {
                        smsd::MCSOptions options;
                        options.nodeLimit = 1;
                        (void)smsd::MCS::find("CC").withOptions(options).matchCorpus(corpus, 2, false, 1);
                    }
                } catch (const std::invalid_argument&) {
                    refused = true;
                }
                require(refused, "Corpus search must refuse malformed storage before its graph access or filter proof");
            }
        }
        smsd::batch::TargetCorpus optional;
        optional.addTargets(targets);
        smsd::CancellationToken cancellation;
        cancellation.cancel();
        require(!optional.substructureCandidates(smsd::parseSMILES("CC"), smsd::ChemOptions{}, 1, &cancellation),
                "A cancelled optional prefilter retains its null mask without graph admission");
        require(smsd::screen_accelerator::active() == nullptr,
                "The optional-prefilter fallback control requires no attached provider");
        require(!optional.substructureCandidates(invalid, smsd::ChemOptions{}, 1),
                "An optional prefilter without a provider retains its null mask without graph access");
    }
    smsd::batch::TargetCorpus valid;
    valid.addTarget(smsd::parseSMILES("CCC"));
    const auto found = smsd::Substructure::find("CC").matchCorpus(valid, false, 1);
    require(found.size() == 1 && found.front().thresholdReached,
            "Valid captured targets remain searchable");

    auto invalid = smsd::parseSMILES("CCC");
    invalid.neighbors[1].push_back(-1);
    smsd::batch::TargetCorpus stopped;
    stopped.addTarget(invalid);
    for (const auto& outcome : smsd::Substructure::find("CC").withTimeout(0).matchCorpus(stopped, false, 1))
        require(outcome.result && outcome.result->status == smsd::SearchStatus::TimedOut
                    && !outcome.result->optimal && !outcome.result->complete,
                "An expired no-filter search retains its stop before detailed graph admission");
    smsd::MCSOptions options;
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    for (const auto& outcome : smsd::MCS::find("CC").withOptions(options).matchCorpus(stopped, 2, false, 1))
        require(outcome.result && outcome.result->status == smsd::SearchStatus::Cancelled
                    && !outcome.result->optimal && !outcome.result->complete,
                "A cancelled no-filter search retains its stop before detailed graph admission");
}

int main() {
    try {
        invalidDenseTopology();
        invalidSparseTopology();
        optionalAdjacencyFallback();
        supportedLegacyStorage();
        capturedInputsValidateBeforeSearch();
        std::cout << "PASS coherent topology admission and optional adjacency fallback\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

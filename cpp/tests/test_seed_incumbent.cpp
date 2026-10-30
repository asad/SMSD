/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#define VertexSearch UnseededVertexSearch
#define injectiveDomainBound unseededInjectiveDomainBound
#define topKIncidentBondUpper unseededTopKIncidentBondUpper
#define childBoundAllowance unseededChildBoundAllowance
#define ChildBoundWork UnseededChildBoundWork
#include "unseeded_vertex_search.hpp"
#undef ChildBoundWork
#undef childBoundAllowance
#undef topKIncidentBondUpper
#undef injectiveDomainBound
#undef VertexSearch

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

template<class Search>
static smsd::SearchResult kernel(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                 const smsd::MCSOptions& options, bool all = false) {
    smsd::search_detail::Budget budget(options, query.n, target.n, true);
    smsd::search_detail::Collector collector(&options, all ? 10 : 1, all);
    smsd::search_detail::QueryView view(query);
    Search search(view, target, options.chemistry, &options, true, all, budget, collector,
                  -1, nullptr, true, 0, true, !all, !all);
    const auto proof = search.run();
    return collector.finish(budget, proof.first, proof.second);
}

static std::pair<size_t, size_t> objective(const smsd::SearchResult& result) {
    return result.matches.empty() ? std::make_pair(size_t(0), size_t(0))
        : std::make_pair(result.matches.front().atomMapping.size(), result.matches.front().bondMapping.size());
}

static smsd::MCSOptions fragmentOptions(uint64_t cap) {
    smsd::MCSOptions options;
    options.connectivity = smsd::MCSConnectivity::Fragments;
    options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    options.nodeLimit = cap;
    options.timeoutMs = 600000;
    return options;
}

using Before = smsd::search_detail::UnseededVertexSearch;
using After = smsd::search_detail::VertexSearch;

static size_t improved = 0;

/**
 * The seeded search is never lower, never loses a proof or a completed status, and visits no
 * more nodes than the search without the seed, of whose nodes it visits a subsequence.
 */
static void checkNeverLower(const smsd::MolGraph& query, const smsd::MolGraph& target, uint64_t cap) {
    const auto options = fragmentOptions(cap);
    const auto before = kernel<Before>(query, target, options), after = kernel<After>(query, target, options);
    assert(objective(after) >= objective(before));
    assert(!before.optimal || after.optimal);
    assert(!before.complete || after.complete);
    assert(before.status != smsd::SearchStatus::Completed || after.status == smsd::SearchStatus::Completed);
    assert(after.nodesVisited <= before.nodesVisited);
    for (const auto& match : after.matches)
        assert(smsd::MCS::validateMapping(query, target, match, options).empty());
    if (objective(after) > objective(before) || after.nodesVisited < before.nodesVisited) ++improved;
    // The same input gives the same score, status, node count and witness.
    const auto again = kernel<After>(query, target, options);
    assert(again.score == after.score && again.status == after.status && again.nodesVisited == after.nodesVisited);
    assert(again.matches.size() == after.matches.size());
    for (size_t index = 0; index < after.matches.size(); ++index) {
        assert(again.matches[index].atomMapping == after.matches[index].atomMapping);
        assert(again.matches[index].bondMapping == after.matches[index].bondMapping);
    }
}

static void checkEveryCap() {
    size_t checks = 0;
    for (const auto& texts : std::vector<std::pair<std::string, std::string>>{
             {"CC(C)CCO", "CCC(C)O"}, {"CCCCN", "CCNN"}, {"C1CC1C", "C1C2CC12"},
             {"N.C.O", "C.O.C"}, {"c1ccccc1CCN", "c1ccncc1CCO"}, {"OC(=O)CCc1ccccc1", "NC(=O)Cc1ccccc1C"}}) {
        const auto query = smsd::parseSMILES(texts.first), target = smsd::parseSMILES(texts.second);
        for (uint64_t cap = 1; cap <= 160; ++cap, ++checks) checkNeverLower(query, target, cap);
    }
    std::cout << "seed_every_cap_checks=" << checks << "\n";
}

static void checkDrugPairs() {
    size_t checks = 0;
    for (const auto& texts : std::vector<std::pair<std::string, std::string>>{
             {"Cc1cnc2c(c1-c1cc(NCc3ccc(CC(=O)O)cc3)ccc1)cccc2C(F)(F)F", "CCn1c2ccc(Cl)cc2[n+](CC)c1C=CNc1ccccc1"},
             {"COC(c1ccccc1NC(=O)Nc1ccc(C(N2CCCC2)=O)cc1)=O", "CCn1c2ccc(Cl)cc2[n+](CC)c1C=CNc1ccccc1"},
             {"COC(c1ccccc1NC(=O)Nc1ccc(C(N2CCCC2)=O)cc1)=O", "Cc1cnc2c(c1-c1cc(NCc3ccc(CC(=O)O)cc3)ccc1)cccc2C(F)(F)F"}}) {
        const auto query = smsd::parseSMILES(texts.first), target = smsd::parseSMILES(texts.second);
        for (uint64_t cap : {1u, 16u, 256u, 4096u, 65536u}) {
            checkNeverLower(query, target, cap);
            ++checks;
        }
    }
    // The seed must change something, or these checks would pass without it.
    assert(improved > 0);
    std::cout << "seed_drug_pair_checks=" << checks << " improved=" << improved << "\n";
}

/** Outside the seeded scope the two searches agree exactly. */
static void checkInactiveRoutes() {
    const auto query = smsd::parseSMILES("c1ccccc1CCN"), target = smsd::parseSMILES("c1ccncc1CCO");
    for (int guard = 0; guard < 6; ++guard) {
        auto options = fragmentOptions(48);
        if (guard == 0) options.connectivity = smsd::MCSConnectivity::Connected;
        if (guard == 1) options.nodeLimit = 0;
        if (guard == 2) options.atomWeights = std::vector<double>(static_cast<size_t>(query.n), 1.0);
        if (guard == 3) options.minFragmentAtoms = 2;
        if (guard == 4) options.maxFragments = 2;
        const bool all = guard == 5;
        const auto before = kernel<Before>(query, target, options, all), after = kernel<After>(query, target, options, all);
        assert(before.score == after.score && objective(before) == objective(after));
        assert(before.status == after.status && before.optimal == after.optimal && before.complete == after.complete);
        assert(before.nodesVisited == after.nodesVisited);
        assert(before.matches.size() == after.matches.size());
        for (size_t index = 0; index < after.matches.size(); ++index)
            assert(before.matches[index].atomMapping == after.matches[index].atomMapping);
    }
    std::cout << "seed_inactive_route_checks=6\n";
}

/**
 * A seed that reaches the upper bound still proves the optimum, at the first node; a search
 * stopped before it starts reports the stop and visits no node.
 */
static void checkBoundAndStops() {
    const auto molecule = smsd::parseSMILES("OC(=O)CCc1ccccc1N");
    const auto before = kernel<Before>(molecule, molecule, fragmentOptions(1u << 20));
    const auto after = kernel<After>(molecule, molecule, fragmentOptions(1u << 20));
    assert(after.optimal && after.complete && objective(after).first == static_cast<size_t>(molecule.n));
    assert(objective(after) == objective(before) && after.nodesVisited <= 1);
    const auto query = smsd::parseSMILES("CCO"), target = smsd::parseSMILES("CCN");
    auto stopped = fragmentOptions(32);
    stopped.timeoutMs = 0;
    const auto timed = smsd::MCS::find(query).withOptions(stopped).match(target);
    assert(timed.status == smsd::SearchStatus::TimedOut && timed.nodesVisited == 0);
    assert(!timed.optimal && !timed.complete);
    stopped.timeoutMs = 600000;
    auto token = std::make_shared<smsd::CancellationToken>();
    token->cancel();
    stopped.cancellation = token;
    const auto cancelled = smsd::MCS::find(query).withOptions(stopped).match(target);
    assert(cancelled.status == smsd::SearchStatus::Cancelled && cancelled.nodesVisited == 0);
    assert(!cancelled.optimal && !cancelled.complete);
    std::cout << "seed_bound_and_stop_checks=3\n";
}

/**
 * Short deadlines stop the search at varying points, often inside the seed; whatever the
 * outcome, the status is consistent and every witness is valid.
 */
static void checkStopsAnywhere() {
    const auto query = smsd::parseSMILES("COC(c1ccccc1NC(=O)Nc1ccc(C(N2CCCC2)=O)cc1)=O");
    const auto target = smsd::parseSMILES("Cc1cnc2c(c1-c1cc(NCc3ccc(CC(=O)O)cc3)ccc1)cccc2C(F)(F)F");
    size_t timedOut = 0;
    for (int64_t milliseconds = 1; milliseconds <= 40; ++milliseconds) {
        auto options = fragmentOptions(1u << 20);
        options.timeoutMs = milliseconds;
        const auto result = kernel<After>(query, target, options);
        assert(result.status != smsd::SearchStatus::TimedOut || (!result.optimal && !result.complete));
        for (const auto& match : result.matches)
            assert(smsd::MCS::validateMapping(query, target, match, options).empty());
        timedOut += result.status == smsd::SearchStatus::TimedOut;
    }
    std::cout << "seed_stop_anywhere_checks=40 timed_out=" << timedOut << "\n";
}

int main() {
    checkEveryCap();
    checkDrugPairs();
    checkInactiveRoutes();
    checkBoundAndStops();
    checkStopsAnywhere();
    std::cout << "seed_incumbent_controls_passed\n";
}

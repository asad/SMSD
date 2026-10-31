/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/smiles_parser.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <set>
#include <string>

namespace {

smsd::MCSOptions targetTwinOptions() {
    smsd::MCSOptions options;
    options.chemistry = smsd::ChemOptions::profile("strict");
    options.timeoutMs = 30000;
    options.nodeLimit = 100000;
    return options;
}

void singleMatchesBelongToTheCompleteOptimumSet() {
    struct Pair { const char* query; const char* target; int atoms; int bonds; bool charge; };
    const Pair pairs[] = {
        {"CCCC", "CC1CC1", 3, 2, false},
        {"CCCCC", "CCC1CC1", 4, 3, false},
        {"CCCC", "C1CCC1", 3, 2, false},
        {"CCOC", "C1COC1", 3, 2, false},
        {"CCCO", "CCOC", 3, 2, false},
        {"CC=CC", "C1=CCC1", 3, 2, false},
        {"CCCO", "C1CCN1", 3, 2, false},
        {"C[CH-]C", "[CH-]1CCC1", 3, 2, true},
    };
    for (const auto& pair : pairs) for (bool reverse : {false, true}) {
        const auto query = smsd::parseSMILES(reverse ? pair.target : pair.query);
        const auto target = smsd::parseSMILES(reverse ? pair.query : pair.target);
        auto options = targetTwinOptions();
        options.chemistry.matchFormalCharge = pair.charge;
        const auto prepared = smsd::MCS::find(query).withOptions(options);
        const auto all = prepared.matchAll(target, 10000);
        const auto single = prepared.match(target);
        assert(all.status == smsd::SearchStatus::Completed && all.optimal && all.complete);
        assert(single.status == smsd::SearchStatus::Completed && single.optimal && single.complete);
        assert(single.matches.size() == 1 && single.score == all.score && single.score == pair.atoms);
        assert(single.matches.front().atomMapping.size() == static_cast<size_t>(pair.atoms));
        assert(single.matches.front().bondMapping.size() == static_cast<size_t>(pair.bonds));
        bool retained = false;
        std::set<std::map<int, int>> mappings;
        for (const auto& match : all.matches) {
            assert(smsd::MCS::validateMapping(query, target, match, options).empty());
            assert(mappings.insert(match.atomMapping).second);
            retained = retained || match.atomMapping == single.matches.front().atomMapping;
        }
        assert(retained);
        assert(smsd::MCS::validateMapping(query, target, single.matches.front(), options).empty());
    }
}

void weightsEnumerationAndPinnedAtomsKeepTheirContracts() {
    const auto graph = smsd::parseSMILES("C1CCC1");
    auto options = targetTwinOptions();
    const auto all = smsd::MCS::find(graph).withOptions(options).matchAll(graph, 100);
    assert(all.optimal && all.complete && all.matches.size() == 8);
    options.atomWeights = {-5, 1, 1, 1};
    const auto weighted = smsd::MCS::find(graph).withOptions(options).match(graph);
    assert(weighted.optimal && weighted.complete && weighted.score == 3);
    assert(weighted.matches.front().atomMapping.size() == 3);
    assert(!weighted.matches.front().atomMapping.count(0));
    options.atomWeights.clear();
    assert(!smsd::MCS::isMappingMaximal(graph, graph, std::map<int, int>{{0, 2}}, options));
    options.excludedTargetAtoms = {0};
    const auto excluded = smsd::MCS::find(graph).withOptions(options).match(graph);
    assert(excluded.optimal && excluded.complete && excluded.score == 3);
    for (const auto& pair : excluded.matches.front().atomMapping) assert(pair.second != 0);
    assert(smsd::MCS::validateMapping(graph, graph, excluded.matches.front(), options).empty());
}

void stoppedSearchesKeepTheirStatus() {
    const auto query = smsd::parseSMILES("CCCCC");
    const auto target = smsd::parseSMILES("C1CCC1");
    auto options = targetTwinOptions();
    options.nodeLimit = 1;
    const auto capped = smsd::MCS::find(query).withOptions(options).match(target);
    assert(capped.status == smsd::SearchStatus::LimitReached && !capped.optimal && !capped.complete);
    assert(capped.nodesVisited <= 1);
    options.timeoutMs = 0;
    const auto timed = smsd::MCS::find(query).withOptions(options).match(target);
    assert(timed.status == smsd::SearchStatus::TimedOut && !timed.optimal && !timed.complete);
    options.timeoutMs = 1000;
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    const auto cancelled = smsd::MCS::find(query).withOptions(options).match(target);
    assert(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete);
}

} // namespace

int main() {
    singleMatchesBelongToTheCompleteOptimumSet();
    weightsEnumerationAndPinnedAtomsKeepTheirContracts();
    stoppedSearchesKeepTheirStatus();
    return 0;
}

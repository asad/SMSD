/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/smiles_parser.hpp"

#include <cassert>
#include <memory>

namespace {

smsd::MCSOptions ringOptions() {
    smsd::MCSOptions options;
    options.chemistry.completeRingsOnly = true;
    options.timeoutMs = 30000;
    options.nodeLimit = 100000;
    return options;
}

void wholeRingsAndTheirBranchesKeepTheMaximum() {
    struct OrdinaryPair { const char* query; const char* target; int atoms; int bonds; };
    const OrdinaryPair pairs[] = {
        {"CC1CCCCC1", "CC1CCCC1", 1, 0},
        {"CCC1CCCC1", "CCC1CCCCC1", 2, 1},
        {"CCOC1CCCCC1", "CCNC1CCCCC1", 6, 6},
        {"CC1CCC2CCCCC2C1", "CC1CCC2CCCC2C1", 7, 7},
        {"CC1CCC2(CC1)CCCC2", "CC1CCC2(CC1)CCCCC2", 7, 7},
    };
    for (const auto& pair : pairs) for (bool reverse : {false, true}) {
        const auto query = smsd::parseSMILES(reverse ? pair.target : pair.query);
        const auto target = smsd::parseSMILES(reverse ? pair.query : pair.target);
        const auto options = ringOptions();
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(!result.matches.empty());
        assert(result.matches.front().atomMapping.size() == static_cast<size_t>(pair.atoms));
        assert(result.matches.front().bondMapping.size() == static_cast<size_t>(pair.bonds));
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty());
    }
}

void ringWeightsAndEnumerationKeepTheirContracts() {
    const auto query = smsd::parseSMILES("CCC1CCCC1"), target = smsd::parseSMILES("CCC1CCCCC1");
    auto options = ringOptions();
    options.atomWeights = {0.25, 0.5, 100, 100, 100, 100, 100};
    const auto weighted = smsd::MCS::find(query).withOptions(options).match(target);
    assert(weighted.optimal && weighted.complete && weighted.score == 0.75);
    assert(weighted.matches.front().atomMapping.size() == 2);
    const auto molecule = smsd::parseSMILES("CC1CCCC1");
    options = ringOptions();
    const auto all = smsd::MCS::find(molecule).withOptions(options).matchAll(molecule, 100);
    assert(all.optimal && all.complete && all.matches.size() == 2);
    for (const auto& match : all.matches)
        assert(smsd::MCS::validateMapping(molecule, molecule, match, options).empty());
}

void ringBudgetsKeepTheirStatus() {
    const auto query = smsd::parseSMILES("CC1CCCCC1"), target = smsd::parseSMILES("CC1CCCC1");
    auto options = ringOptions();
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
    wholeRingsAndTheirBranchesKeepTheMaximum();
    ringWeightsAndEnumerationKeepTheirContracts();
    ringBudgetsKeepTheirStatus();
    return 0;
}

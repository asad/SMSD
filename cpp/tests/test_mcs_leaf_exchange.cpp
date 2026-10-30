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
#include <string>
#include <utility>

namespace {

smsd::MCSOptions strictOptions() {
    smsd::MCSOptions options;
    options.chemistry = smsd::ChemOptions::profile("strict");
    options.timeoutMs = 1000;
    options.nodeLimit = 100000;
    return options;
}

void pendantBranches() {
    struct Pair { const char* query; const char* target; int atoms; int bonds; };
    const Pair pairs[] = {
        {"CCC(C)CC", "CCCCCC", 5, 4},
        {"CC(C)(C)OCC", "CC(C)(C)NCC", 4, 3},
        {"CC(C)NCC(C)O", "CC(C)NCC(C)N", 7, 6},
        {"CC(C)(C)CC(C)(C)O", "CC(C)(C)CCC(C)(C)N", 7, 6},
        {"OC(C)(C)CC(C)(C)C", "NC(C)(C)CCC(C)(C)C", 7, 6},
    };
    for (const auto& pair : pairs) for (bool reverse : {false, true}) {
        const auto query = smsd::parseSMILES(reverse ? pair.target : pair.query);
        const auto target = smsd::parseSMILES(reverse ? pair.query : pair.target);
        const auto options = strictOptions();
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(result.matches.front().atomMapping.size() == static_cast<size_t>(pair.atoms));
        assert(result.matches.front().bondMapping.size() == static_cast<size_t>(pair.bonds));
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty());
    }
}

void isotopeWildcardsKeepAlternativePairs() {
    auto options = strictOptions();
    options.chemistry.matchIsotope = true;
    const auto query = smsd::parseSMILES("CN[13CH2]O");
    const auto target = smsd::parseSMILES("CN[12CH2]S");
    const auto result = smsd::MCS::find(query).withOptions(options).match(target);
    assert(result.optimal && result.complete);
    assert(result.matches.front().atomMapping.size() == 3);
    assert(result.matches.front().bondMapping.size() == 2);
    assert(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty());
}

void weightsAndEnumerationKeepTheirContracts() {
    const auto molecule = smsd::parseSMILES("CC(C)C");
    auto options = strictOptions();
    const auto all = smsd::MCS::find(molecule).withOptions(options).matchAll(molecule, 100);
    assert(all.optimal && all.complete && all.matches.size() == 6);
    options.atomWeights = {-5, 1, 1, 1};
    const auto weighted = smsd::MCS::find(molecule).withOptions(options).match(molecule);
    assert(weighted.optimal && weighted.complete && weighted.score == 3);
    assert(weighted.matches.front().atomMapping.size() == 3);
    assert(!weighted.matches.front().atomMapping.count(0));
}

void stoppedSearchesKeepTheirStatus() {
    const auto query = smsd::parseSMILES("CC(C)(C)CC(C)(C)O");
    const auto target = smsd::parseSMILES("CC(C)(C)CCC(C)(C)N");
    auto options = strictOptions();
    options.nodeLimit = 1;
    const auto capped = smsd::MCS::find(query).withOptions(options).match(target);
    assert(capped.status == smsd::SearchStatus::LimitReached && !capped.optimal && !capped.complete);
    assert(capped.nodesVisited <= options.nodeLimit);
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
    pendantBranches();
    isotopeWildcardsKeepAlternativePairs();
    weightsAndEnumerationKeepTheirContracts();
    stoppedSearchesKeepTheirStatus();
    return 0;
}

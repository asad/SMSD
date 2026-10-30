/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/detail/bridge_bound.hpp"

#include <cassert>
#include <memory>
#include <string>

namespace {

void branchBoundsUseTheEnabledPredicates() {
    struct Pair { const char* query; const char* target; const char* profile; int atoms; };
    const Pair pairs[] = {
        {"CCC(C)CC", "CCCCCC", "default", 5},
        {"CC(C)(C)OCC", "CC(C)(C)NCC", "strict", 4},
        {"CC=CC", "CCCC", "compat-fmcs", 4},
        {"CC=CC", "CCCC", "default", 2},
        {"CN[13CH2]O", "CN[12CH2]S", "strict", 3},
        {"C.CC", "CCC", "default", 2},
        {"N", "O", "default", 0},
    };
    for (const auto& pair : pairs) {
        const auto query = smsd::parseSMILES(pair.query), target = smsd::parseSMILES(pair.target);
        smsd::MCSOptions options;
        options.timeoutMs = 1000;
        options.chemistry = smsd::ChemOptions::profile(pair.profile);
        options.chemistry.matchIsotope = true;
        smsd::search_detail::Budget budget(options, query.n, target.n, true);
        const auto bound = smsd::search_detail::BridgeBound::create(query, target, options.chemistry, budget);
        assert(bound && bound->upper() == pair.atoms && budget.nodesVisited() == 0);
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        assert(result.optimal && result.complete);
        assert((result.matches.empty() ? 0 : result.matches.front().atomMapping.size()) == static_cast<size_t>(pair.atoms));
        if (!result.matches.empty()) for (const auto& atom : result.matches.front().atomMapping)
            assert(bound->bound(atom.first, atom.second) >= pair.atoms);
    }
}

void unsupportedOrStoppedConstructionLeavesNoTable() {
    const auto chain = smsd::parseSMILES("CCCCCC"), ring = smsd::parseSMILES("C1CCCCC1");
    smsd::MCSOptions options;
    options.timeoutMs = 1000;
    smsd::search_detail::Budget budget(options, ring.n, chain.n, true);
    assert(!smsd::search_detail::BridgeBound::create(ring, chain, options.chemistry, budget));
    const auto large = smsd::parseSMILES(std::string(513, 'C'));
    assert(!smsd::search_detail::BridgeBound::create(large, chain, options.chemistry, budget));
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    smsd::search_detail::Budget cancelled(options, chain.n, chain.n, true);
    assert(!smsd::search_detail::BridgeBound::create(chain, chain, options.chemistry, cancelled));
    assert(cancelled.status() == smsd::SearchStatus::Cancelled);
}

void largeBranchedChainsKeepTheMaximumPath() {
    const auto query = smsd::parseSMILES(std::string(64, 'C'));
    const auto target = smsd::parseSMILES(std::string(31, 'C') + "C(C)" + std::string(31, 'C'));
    for (bool reverse : {false, true}) {
        smsd::MCSOptions options;
        options.timeoutMs = 3000;
        options.nodeLimit = 1000000;
        const auto& first = reverse ? target : query;
        const auto& second = reverse ? query : target;
        const auto result = smsd::MCS::find(first).withOptions(options).match(second);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(result.matches.front().atomMapping.size() == 63 && result.matches.front().bondMapping.size() == 62);
        assert(smsd::MCS::validateMapping(first, second, result.matches.front(), options).empty());
    }
}

} // namespace

int main() {
    branchBoundsUseTheEnabledPredicates();
    unsupportedOrStoppedConstructionLeavesNoTable();
    largeBranchedChainsKeepTheMaximumPath();
    return 0;
}

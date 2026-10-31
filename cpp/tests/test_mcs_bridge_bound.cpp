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
#include <string>

namespace {

void cyclicRootsKeepTheirWholeParts() {
    const auto molecule = smsd::parseSMILES("CC1CCCCC1");
    smsd::MCSOptions options;
    options.timeoutMs = 30000;
    smsd::search_detail::Budget budget(options, molecule.n, molecule.n, true);
    const auto bound = smsd::search_detail::BridgeBound::create(molecule, molecule, options.chemistry, budget);
    assert(bound && !bound->acyclic() && bound->upper() == 7);
    assert(bound->bound(1, 1) == 7);
    assert(budget.nodesVisited() == 0);
    const auto ring = smsd::parseSMILES("C1CCCCC1");
    assert(!smsd::search_detail::BridgeBound::create(ring, molecule, options.chemistry, budget));
}

void relaxedPredicatesRemainUpperBounds() {
    struct BridgePair { const char* query; const char* target; const char* profile; bool typed; };
    const BridgePair pairs[] = {
        {"CC1CC1", "CC1CCC1", "default", true},
        {"CC=C1CCCCC1", "CCC1CCCCC1", "compat-fmcs", true},
        {"CC1CCC2CCCCC2C1", "CC1CCC2CCCC2C1", "strict", true},
        {"CC1CCC2(CC1)CCCC2", "CC1CCC2(CC1)CCCCC2", "default", true},
        {"OC1CCCCC1", "NC1CCCCC1", "default", false},
        {"CN[13CH]1CCCCC1", "CN[12CH]1CCCCC1", "default", true},
    };
    for (const auto& pair : pairs) {
        const auto query = smsd::parseSMILES(pair.query), target = smsd::parseSMILES(pair.target);
        smsd::MCSOptions options;
        options.timeoutMs = 30000;
        options.chemistry = smsd::ChemOptions::profile(pair.profile);
        options.chemistry.matchAtomType = pair.typed;
        options.chemistry.matchIsotope = true;
        smsd::search_detail::Budget budget(options, query.n, target.n, true);
        const auto bound = smsd::search_detail::BridgeBound::create(query, target, options.chemistry, budget);
        assert(bound && !bound->acyclic());
        const auto result = smsd::MCS::find(query).withOptions(options).match(target);
        assert(result.optimal && result.complete && !result.matches.empty());
        const int atoms = static_cast<int>(result.matches.front().atomMapping.size());
        assert(bound->upper() >= atoms);
        for (const auto& atom : result.matches.front().atomMapping)
            assert(bound->bound(atom.first, atom.second) >= atoms);
        assert(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty());
    }
}

void aLargeCyclicMatchCanHaveOneBondPerAtom() {
    const auto query = smsd::parseSMILES(std::string(58, 'C') + "C1CCCCC1");
    const auto target = smsd::parseSMILES(std::string(29, 'C') + "C(C)" + std::string(27, 'C') + "C1CCCCC1");
    for (bool reverse : {false, true}) {
        smsd::MCSOptions options;
        options.timeoutMs = 30000;
        options.nodeLimit = 1000000;
        const auto& first = reverse ? target : query;
        const auto& second = reverse ? query : target;
        const auto result = smsd::MCS::find(first).withOptions(options).match(second);
        assert(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete);
        assert(result.matches.front().atomMapping.size() == 63);
        assert(result.matches.front().bondMapping.size() == 63);
        assert(smsd::MCS::validateMapping(first, second, result.matches.front(), options).empty());
    }
}

} // namespace

int main() {
    cyclicRootsKeepTheirWholeParts();
    relaxedPredicatesRemainUpperBounds();
    aLargeCyclicMatchCanHaveOneBondPerAtom();
    return 0;
}

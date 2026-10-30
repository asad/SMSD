/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/query_smarts.hpp"
#include <cassert>
#include <iostream>

int main() {
    const std::vector<std::string> molecules = {"C", "N", "CC", "C=C", "C#C", "CO", "C=O",
        "[NH4+]", "[13CH4]", "[12CH4]", "C1CC1", "c1ccccc1", "c1ccncc1", "CC.O"};
    for (bool strictAromaticity : {false, true})
        for (bool anyBond : {false, true})
            for (bool charges : {false, true})
                for (bool isotopes : {false, true}) {
                    smsd::ChemOptions chemistry;
                    chemistry.aromaticityMode = strictAromaticity ? smsd::ChemOptions::AromaticityMode::STRICT
                        : smsd::ChemOptions::AromaticityMode::FLEXIBLE;
                    chemistry.matchBondOrder = anyBond ? smsd::ChemOptions::BondOrderMode::ANY
                        : smsd::ChemOptions::BondOrderMode::STRICT;
                    chemistry.matchFormalCharge = charges;
                    chemistry.matchIsotope = isotopes;
                    for (const auto& querySmiles : molecules) {
                        const auto query = smsd::parseSMILES(querySmiles);
                        const auto pattern = smsd::Substructure::smarts(smsd::toQuerySMARTS(query, chemistry));
                        auto molecular = smsd::Substructure::find(query);
                        auto options = smsd::Substructure::Options();
                        options.chemistry = chemistry;
                        molecular = molecular.withOptions(options);
                        for (const auto& targetSmiles : molecules) {
                            const auto target = smsd::parseSMILES(targetSmiles);
                            const bool expected = !molecular.match(target).matches.empty();
                            const bool actual = !pattern.match(target).matches.empty();
                            if (expected != actual) {
                                std::cerr << querySmiles << " => " << smsd::toQuerySMARTS(query, chemistry)
                                    << " on " << targetSmiles << " expected " << expected << " actual " << actual << '\n';
                                return 1;
                            }
                        }
                    }
                }
    auto query = smsd::parseSMILES("C1CC1");
    smsd::Match selected;
    selected.atomMapping = {{0, 0}, {1, 1}, {2, 2}};
    selected.bondMapping = {{0, 1, 0, 1}, {1, 2, 1, 2}};
    auto pattern = smsd::Substructure::smarts(selected.toSmarts(query));
    assert(!pattern.match(smsd::parseSMILES("CCC")).matches.empty());
    smsd::ChemOptions chemistry;
    chemistry.ringMatchesRingOnly = true;
    auto ringPattern = smsd::Substructure::smarts(smsd::toQuerySMARTS(query, chemistry));
    assert(!ringPattern.match(query).matches.empty());
    assert(ringPattern.match(smsd::parseSMILES("CCC")).matches.empty());
    chemistry.useChirality = true;
    bool rejected = false;
    try { selected.toSmarts(query, chemistry); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    const auto exported = smsd::MCS::find("CCO").matchWithSmarts("CCN");
    assert(exported.result.optimal && !exported.smarts.empty());
    assert(!smsd::Substructure::smarts(exported.smarts).match("CCCl").matches.empty());
    const auto stopped = smsd::MCS::find("CCO").withTimeout(0).matchWithSmarts("CCN");
    assert(stopped.smarts.empty() && !stopped.result.optimal);
    const int chainAtoms = 12000;
    std::vector<std::vector<int>> neighbors(chainAtoms), orders(chainAtoms);
    for (int atom = 1; atom < chainAtoms; ++atom) {
        neighbors[atom - 1].push_back(atom);
        neighbors[atom].push_back(atom - 1);
        orders[atom - 1].push_back(1);
        orders[atom].push_back(1);
    }
    const auto chain = smsd::MolGraph::Builder().atomCount(chainAtoms)
        .atomicNumbers(std::vector<int>(chainAtoms, 6)).setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(orders)).build(false);
    std::string expectedChain = "[#6]";
    for (int atom = 1; atom < chainAtoms; ++atom) expectedChain += "-,:[#6]";
    assert(smsd::toQuerySMARTS(chain) == expectedChain);
    std::cout << "Query SMARTS round trips passed\n";
}

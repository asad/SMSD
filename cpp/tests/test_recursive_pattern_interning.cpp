/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/substructure.hpp"

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

static void recursivePointers(const smsd::smarts_detail::AtomExpr& expression,
                              std::vector<const smsd::SmartsQuery*>& pointers) {
    if (expression.cachedRecursive) pointers.push_back(expression.cachedRecursive.get());
    for (const auto& child : expression.children) recursivePointers(child, pointers);
}

static void recursiveOccurrenceCounts() {
    const auto flat = smsd::parseSMARTS("[$(CC)].[$(CC)]");
    const smsd::search_detail::PatternPlan flatPlan(flat);
    assert(flatPlan.subpatternCount() == 2);
    const auto nested = smsd::parseSMARTS("[$([$([#6])])].[$([$([#6])])]");
    const smsd::search_detail::PatternPlan nestedPlan(nested);
    assert(nestedPlan.subpatternCount() == 4);
    const auto mixed = smsd::parseSMARTS("[$(CCO)]~[$(NCS)]~[$(CCO)]~[$(C(C)O)]~[$([C;$(CCO)]=O)]");
    const smsd::search_detail::PatternPlan mixedPlan(mixed);
    assert(mixedPlan.subpatternCount() == 6);
    for (const auto& atom : flat.atoms) {
        std::vector<const smsd::SmartsQuery*> payloads;
        recursivePointers(atom.expr, payloads);
        for (const auto* payload : payloads) {
            assert(flatPlan.subpatternIndex(*payload) == 0);
            assert(flatPlan.subpattern(*payload) == &flatPlan.subpatternLayout(0));
        }
    }
}

int main(int argc, char**) {
    recursiveOccurrenceCounts();
    auto first = smsd::parseSMARTS("[$(CC)].[$(CC)]");
    auto second = smsd::parseSMARTS("[$(CC)].[$(CC)]");
    std::vector<const smsd::SmartsQuery*> pointers, separate;
    for (const auto& atom : first.atoms) recursivePointers(atom.expr, pointers);
    for (const auto& atom : second.atoms) recursivePointers(atom.expr, separate);
    assert(pointers.size() == 2 && separate.size() == 2 && pointers.front() != separate.front());
    if (argc > 1 && pointers[0] != pointers[1]) {
        std::cerr << "Repeated recursive payloads are not shared within this parse\n";
        return 1;
    }
    size_t atoms = 0, nodes = 0;
    (void)smsd::smarts_detail::parsePattern("[$(CC)].[$(CC)]", 8, atoms, nodes);
    assert(atoms == 6 && nodes == 8);
    atoms = 0;
    nodes = smsd::MAX_SMARTS_EXPRESSION_NODES - 5;
    bool nodeLimit = false;
    try { (void)smsd::smarts_detail::parsePattern("[$(CC)].[$(CC)]", 8, atoms, nodes); }
    catch (const std::invalid_argument&) { nodeLimit = true; }
    assert(nodeLimit);
    const std::string large = "[$(" + std::string(6000, 'C') + "),$(" + std::string(6000, 'C') + ")]";
    bool atomLimit = false;
    try { (void)smsd::parseSMARTS(large); }
    catch (const std::invalid_argument&) { atomLimit = true; }
    assert(atomLimit);
    bool depthLimit = false;
    try { (void)smsd::parseSMARTS("[$([$([#6])])]", 1); }
    catch (const std::invalid_argument&) { depthLimit = true; }
    assert(depthLimit);
    (void)smsd::parseSMARTS("[$([$([#6])])]", 2);
    for (const std::string& target : {"CCC", "CCO", "C1CC1", "O", "CC.CC"}) {
        const auto repeated = smsd::Substructure::smarts("[$(CC),$(CC)]").matchAll(target, 32);
        const auto ordinary = smsd::Substructure::smarts("[$(CC)]").matchAll(target, 32);
        assert(repeated.status == ordinary.status && repeated.optimal == ordinary.optimal
            && repeated.complete == ordinary.complete && repeated.score == ordinary.score
            && repeated.matches.size() == ordinary.matches.size());
    }
    std::cout << "PASS recursive payload sharing, expanded limits and matching semantics\n";
}

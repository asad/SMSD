/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/substructure.hpp"

#include <memory>

int main() {
    auto pattern = smsd::parseSMARTS("[$([#6])].[$([#6])]");
    const auto ordinary = smsd::Substructure::smarts(pattern).matchAll("CC");
    smsd_test::require(ordinary.complete && ordinary.matches.size() == 2,
                      "Equivalent recursive patterns must retain both ordinary embeddings");

    auto& second = pattern.atoms[1].expr;
    smsd_test::require(second.op == smsd::smarts_detail::ExprOp::PRIM
                      && second.primType == smsd::smarts_detail::AtomPrimType::RECURSIVE,
                      "The fixture must contain a recursive atom expression");
    second.cachedRecursive = std::make_shared<smsd::SmartsQuery>(smsd::parseSMARTS("[#7]"));
    const auto mixed = smsd::Substructure::smarts(pattern).matchAll("CN");
    smsd_test::require(mixed.complete && mixed.optimal && mixed.matches.size() == 1,
                      "Distinct compiled recursive predicates must each constrain their atom");
    smsd_test::require(mixed.matches.front().atomMapping.at(0) == 0
                      && mixed.matches.front().atomMapping.at(1) == 1,
                      "The embedding must follow the carbon and nitrogen compiled predicates");
    const auto negative = smsd::Substructure::smarts(pattern).matchAll("CC");
    smsd_test::require(negative.complete && negative.optimal && negative.matches.empty(),
                      "A different recursive cache must not inherit the carbon-only domain");
    return 0;
}

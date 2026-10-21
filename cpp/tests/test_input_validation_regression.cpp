/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/batch.hpp"
#include "smsd/layout.hpp"
#include "smsd/smiles_parser.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

template<class Exception, class Operation>
static void requireError(Operation operation) {
    try { operation(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("Expected input validation error");
}

int main() {
    const auto query = smsd::parseSMILES("CCO");
    smsd::ChemOptions chemistry;
    smsd::MCSOptions options;
    for (const auto& weights : {std::vector<double>{1},
            std::vector<double>{1, std::numeric_limits<double>::quiet_NaN(), 1},
            std::vector<double>{1, std::numeric_limits<double>::infinity(), 1}}) {
        options.atomWeights = weights;
        requireError<std::invalid_argument>([&] { smsd_test::screenedCommon(query, {}, chemistry, options); });
        requireError<std::invalid_argument>([&] { smsd_test::screenedCommonSizes(query, {}, chemistry, options); });
        smsd::batch::TargetCorpus corpus;
        requireError<std::invalid_argument>([&] { smsd::MCS::find(query).withOptions(options).matchMany(corpus.targets()); });
    }
    for (int atom : {-1, query.n, 100000000}) {
        requireError<std::out_of_range>([&] { smsd_test::fragmentSmiles(query, {{atom, 0}}); });
        requireError<std::out_of_range>([&] { query.ringSystemOf(atom); });
        requireError<std::out_of_range>([&] { smsd::extractSubgraph(query, {atom}); });
    }
    if (query.ringSystemOf(0) != -1 || smsd_test::fragmentSmiles(query, {{0, 0}, {1, 1}}) != "CC"
        || smsd::extractSubgraph(query, {1, 2}).n != 2)
        throw std::runtime_error("Valid graph inputs must remain supported");

    std::vector<smsd::Point2D> coordinates{{0, 0}, {1, 0}};
    requireError<std::invalid_argument>([&] { smsd::transform::normaliseBondLength(query, coordinates); });
    coordinates.push_back({2, 0});
    smsd::transform::normaliseBondLength(query, coordinates);
    if (std::abs(coordinates[1].x - coordinates[0].x - 1.5) > 1e-8)
        throw std::runtime_error("Valid coordinates must still be scaled");

    smsd::Match selected;
    selected.atomMapping = {{0, 0}, {1, 1}};
    for (const auto& bond : {std::array<int, 4>{-1, 1, 0, 1}, {0, 100000000, 0, 1}}) {
        selected.bondMapping = {bond};
        requireError<std::out_of_range>([&] { selected.toGraph(query); });
        requireError<std::out_of_range>([&] { selected.toSmiles(query); });
    }
    for (const auto& bond : {std::array<int, 4>{0, 2, 0, 2}, {0, 0, 0, 0}, {0, 1, 1, 0}}) {
        selected.bondMapping = {bond};
        requireError<std::invalid_argument>([&] { selected.toGraph(query); });
        requireError<std::invalid_argument>([&] { selected.toSmiles(query); });
    }
    selected.bondMapping = {{0, 1, 0, 1}, {0, 1, 0, 1}};
    requireError<std::invalid_argument>([&] { selected.toGraph(query); });
    options = {};
    for (int atom : {-1, query.n, 100000000}) {
        options.excludedTargetAtoms = {atom};
        requireError<std::out_of_range>([&] { smsd::MCS::find(query).withOptions(options).match(query); });
    }
    options = {};
    options.minFragmentAtoms = 0;
    requireError<std::invalid_argument>([&] { smsd::MCS::find(query).withOptions(options).match(query); });
    options = {};
    options.maxFragments = -1;
    requireError<std::invalid_argument>([&] { smsd::MCS::find(query).withOptions(options).match(query); });
    options = {};
    options.objective = smsd::MCS::Objective::MCES;
    options.atomWeights = {1, 1, 1};
    requireError<std::invalid_argument>([&] { smsd::MCS::find(query).withOptions(options).match(query); });
    requireError<std::invalid_argument>([&] { smsd::MCS::find(query).matchAll(query, 0); });
    if (smsd::MCS::find(query).match(query).score != 3.0)
        throw std::runtime_error("Valid molecular input must remain supported");
    std::cout << "Input validation regressions passed\n";
}

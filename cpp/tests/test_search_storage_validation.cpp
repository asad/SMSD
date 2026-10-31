/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Mutation = std::pair<std::string, std::function<void(smsd::MolGraph&)>>;

static const std::vector<Mutation> mutations{
    {"negative-count", [](auto& graph) { graph.n = -1; }},
    {"larger-count", [](auto& graph) { graph.n = 4; }},
    {"negative-neighbour", [](auto& graph) { graph.neighbors[1].push_back(-1); }},
    {"large-neighbour", [](auto& graph) { graph.neighbors[1].push_back(graph.n); }},
    {"adjacency-rows", [](auto& graph) { graph.neighbors.resize(2); }},
    {"atomic-numbers", [](auto& graph) { graph.atomicNum.resize(2); }},
    {"degrees", [](auto& graph) { graph.degree.resize(2); }},
    {"aromatic-flags", [](auto& graph) { graph.aromatic.resize(2); }},
    {"ring-flags", [](auto& graph) { graph.ring.resize(2); }},
    {"dense-rows", [](auto& graph) { graph.bondPacked.resize(2); }},
    {"dense-columns", [](auto& graph) { graph.bondPacked[1].resize(2); }},
    {"tautomer-classes", [](auto& graph) { graph.tautomerClass.resize(2); }},
    {"stereo-rows", [](auto& graph) { graph.hasDbStereo = true; graph.dbStereoConf.resize(2); }},
    {"stereo-columns", [](auto& graph) { graph.hasDbStereo = true; graph.dbStereoConf[1].resize(2); }}
};

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Operation>
static void refuses(Operation operation) {
    try {
        const auto result = operation();
        if (result.status == smsd::SearchStatus::LimitReached && !result.optimal && !result.complete) return;
    } catch (const std::invalid_argument&) { return; }
      catch (const std::out_of_range&) { return; }
    throw std::runtime_error("Invalid graph storage was not clearly refused");
}

static void targetCase(const std::string& operation, const Mutation& mutation) {
    const auto query = smsd::parseSMILES("CC");
    auto target = smsd::parseSMILES("CCC");
    mutation.second(target);
    if (operation == "sub") {
        const auto prepared = smsd::Substructure::find(query).withTimeout(30000);
        refuses([&] { return prepared.match(target); });
    } else if (operation == "smarts") {
        const auto prepared = smsd::Substructure::smarts("[#6]").withTimeout(30000);
        refuses([&] { return prepared.match(target); });
    } else if (operation == "validate" || operation == "mapping") {
        bool refused = false;
        try {
            if (operation == "validate") (void)smsd::MCS::validateMapping(query, target, smsd::Match{});
            else (void)smsd::MCS::mappingMatch(query, target, {{0,0},{1,1}});
        } catch (const std::invalid_argument&) { refused = true; }
          catch (const std::out_of_range&) { refused = true; }
        require(refused, "Mapping entry points must refuse invalid graph storage");
    } else {
        smsd::MCSOptions options;
        options.timeoutMs = 30000;
        options.connectivity = smsd::MCSConnectivity::Fragments;
        options.objective = operation == "mces" ? smsd::MCSObjective::MCES : smsd::MCSObjective::MCIS;
        const auto prepared = smsd::MCS::find(query).withOptions(options);
        refuses([&] { return prepared.match(target); });
    }
}

static void queryCase(const Mutation& mutation) {
    auto query = smsd::parseSMILES("CCC");
    mutation.second(query);
    bool subRefused = false, mcsRefused = false;
    try { (void) smsd::Substructure::find(query); }
    catch (const std::invalid_argument&) { subRefused = true; }
    catch (const std::out_of_range&) { subRefused = true; }
    try { (void) smsd::MCS::find(query); }
    catch (const std::invalid_argument&) { mcsRefused = true; }
    catch (const std::out_of_range&) { mcsRefused = true; }
    require(subRefused && mcsRefused, "Prepared queries must refuse invalid graph storage");
}

static void stopped() {
    auto target = smsd::parseSMILES("CCC");
    target.neighbors[1].push_back(-1);
    const auto sub = smsd::Substructure::find("CC").withTimeout(0).match(target);
    const auto pattern = smsd::Substructure::smarts("[#6]").withTimeout(0).match(target);
    require(sub.status == smsd::SearchStatus::TimedOut && !sub.optimal && !sub.complete
        && pattern.status == smsd::SearchStatus::TimedOut && !pattern.optimal && !pattern.complete,
        "Expired searches must retain their immediate stopped result");
    smsd::MCSOptions options;
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    const auto common = smsd::MCS::find("CC").withOptions(options).match(target);
    require(common.status == smsd::SearchStatus::Cancelled && !common.optimal && !common.complete,
        "Pre-cancelled search must retain its stopped result");
}

static void supported() {
    const auto prepared = smsd::Substructure::find("CC");
    for (int mode = 0; mode < 5; ++mode) {
        auto target = smsd::parseSMILES("CCC");
        if (mode == 0) target.neighbors[1].push_back(0);
        if (mode == 1) target.neighbors[2].push_back(2);
        if (mode == 2) target.bondPacked[1][0] = 2;
        if (mode == 3) target.bondPacked[1][0] |= smsd::MolGraph::BOND_AROM_BIT;
        if (mode == 4) target.bondPacked[1][0] |= smsd::MolGraph::BOND_RING_BIT;
        const auto result = prepared.matchAll(target);
        require(result.status == smsd::SearchStatus::Completed && result.complete,
            "Supported neighbour listings/directional properties must remain accepted");
    }
    require(!smsd::MCS::validateMapping(smsd::parseSMILES("CC"), smsd::parseSMILES("CCC"),
        std::map<int,int>{{-1,0}}).empty(), "Invalid correspondences retain mapping errors");
    const auto ordinary = smsd::MCS::find("CC").match("CCC");
    require(ordinary.optimal && ordinary.complete && ordinary.score == 2,
        "Ordinary valid search must retain its optimum");
    auto optional = smsd::parseSMILES("CCC");
    optional.formalCharge.clear();
    optional.massNumber.clear();
    optional.tetraChirality.clear();
    optional.adjLong.resize(1);
    require(prepared.match(optional).status == smsd::SearchStatus::Completed,
        "Optional property and bit-adjacency fallbacks must remain supported");
}

int main(int argc, char** argv) {
    try {
        const std::string selected = argc > 1 ? argv[1] : "all";
        const std::string selectedMutation = argc > 2 ? argv[2] : "all";
        int exercised = 0;
        for (const auto& mutation : mutations) {
            if (selectedMutation != "all" && selectedMutation != mutation.first) continue;
            for (const auto& operation : {"sub", "smarts", "mcis", "mces", "query", "validate", "mapping"}) {
                if (selected != "all" && selected != operation) continue;
                if (std::string(operation) == "query") queryCase(mutation); else targetCase(operation, mutation);
                ++exercised;
            }
        }
        if (selected == "all" || selected == "stopped") { stopped(); ++exercised; }
        if (selected == "all" || selected == "supported") { supported(); ++exercised; }
        require(exercised > 0, "No matching regression case");
        std::cout << "PASS " << exercised << " storage-admission cases\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

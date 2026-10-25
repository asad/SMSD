/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/ring_finder.hpp"

#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const char* stackSmiles =
    "C12C3C4C5C6C1C17C28C32C43C54C61C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C"
    "82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17"
    "C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C4"
    "3C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C"
    "15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34"
    "C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C7"
    "1C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C1C2C3C"
    "4C6C51";

void check(bool condition, const std::string& message, int& failures) {
    if (!condition) {
        std::cerr << message << '\n';
        ++failures;
    }
}

void runCase(const std::string& name, const std::function<void()>& test, int& failures) {
    try {
        test();
    } catch (const std::exception& error) {
        check(false, name + ": " + error.what(), failures);
    }
}

smsd::MCSOptions options(smsd::MCSObjective objective, smsd::MCSConnectivity connectivity) {
    smsd::MCSOptions result;
    result.objective = objective;
    result.connectivity = connectivity;
    result.timeoutMs = 60000;
    return result;
}

/** Four-member rings linked around one larger ring. */
smsd::MolGraph fourRingChain() {
    constexpr int ringCount = 12;
    constexpr int atomCount = 1 + 3 * ringCount;
    std::vector<std::vector<int>> neighbors(atomCount), orders(atomCount);
    const auto addBond = [&](int first, int second) {
        neighbors[first].push_back(second);
        neighbors[second].push_back(first);
        orders[first].push_back(1);
        orders[second].push_back(1);
    };
    int previous = 0;
    for (int ring = 0; ring < ringCount; ++ring) {
        const int firstSide = 3 * ring + 1, secondSide = firstSide + 1;
        const int next = firstSide + 2;
        addBond(previous, firstSide);
        addBond(previous, secondSide);
        addBond(firstSide, next);
        addBond(secondSide, next);
        previous = next;
    }
    addBond(previous, 0);
    return smsd::MolGraph::Builder().atomCount(atomCount)
        .atomicNumbers(std::vector<int>(atomCount, 6)).setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(orders)).build(false);
}

void checkMaximalityUnproved(const smsd::MolGraph& query, const smsd::MolGraph& target,
                            const smsd::Match& match, const smsd::MCSOptions& settings,
                            const std::string& name, int& failures) {
    bool matchThrew = false, atomsThrew = false;
    try {
        smsd::MCS::isMappingMaximal(query, target, match, settings);
    } catch (const std::runtime_error& error) {
        matchThrew = std::string(error.what()).find("maximality") != std::string::npos;
    }
    try {
        smsd::MCS::isMappingMaximal(query, target, match.atomMapping, settings);
    } catch (const std::runtime_error& error) {
        atomsThrew = std::string(error.what()).find("maximality") != std::string::npos;
    }
    check(matchThrew, name + ": match maximality must remain unproved", failures);
    check(atomsThrew, name + ": atom mapping maximality must remain unproved", failures);
}

void checkUnproved(const smsd::MolGraph& query, const smsd::MolGraph& target,
                   const smsd::MCSOptions& settings, const smsd::SearchResult& result,
                   const std::string& name, int& failures) {
    check(result.status == smsd::SearchStatus::LimitReached,
          name + ": expected LimitReached", failures);
    check(!result.optimal, name + ": optimum must remain unproved", failures);
    check(!result.complete, name + ": search must remain incomplete", failures);
    const smsd::Match match = result.matches.empty() ? smsd::Match{} : result.matches.front();
    check(smsd::MCS::validateMapping(query, target, match, settings).empty(),
          name + ": returned mapping must remain valid", failures);
    checkMaximalityUnproved(query, target, match, settings, name, failures);
}

void checkCompleted(const smsd::MolGraph& query, const smsd::MolGraph& target,
                    const smsd::MCSOptions& settings, const smsd::SearchResult& result,
                    const std::string& name, int& failures) {
    check(result.status == smsd::SearchStatus::Completed, name + ": expected Completed", failures);
    check(result.optimal, name + ": expected proved optimum", failures);
    check(result.complete, name + ": expected complete search", failures);
    for (const auto& match : result.matches)
        check(smsd::MCS::validateMapping(query, target, match, settings).empty(),
              name + ": returned mapping must remain valid", failures);
}

} // namespace

int main() {
    int failures = 0;
    try {
        const auto stack = smsd::parseSMILES(stackSmiles);
        const auto cyclobutane = smsd::parseSMILES("C1CCC1");
        check(smsd::ringsTruncated(stack), "stack must have truncated ring perception", failures);
        check(!smsd::ringsTruncated(cyclobutane), "cyclobutane rings must be complete", failures);
        if (failures) return 1;

        for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
            for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
                for (bool reverse : {false, true}) {
                    const auto& query = reverse ? cyclobutane : stack;
                    const auto& target = reverse ? stack : cyclobutane;
                    const std::string name =
                        std::string(objective == smsd::MCSObjective::MCIS ? "MCIS" : "MCES") + "/" +
                        (connectivity == smsd::MCSConnectivity::Connected ? "Connected" : "Fragments") +
                        (reverse ? "/reverse" : "/forward");
                    runCase(name + "/complete-rings", [&] {
                        auto settings = options(objective, connectivity);
                        settings.chemistry.completeRingsOnly = true;
                        auto result = smsd::MCS::find(query).withOptions(settings).match(target);
                        checkUnproved(query, target, settings, result, name, failures);
                    }, failures);
                    runCase(name + "/unconstrained", [&] {
                        const auto settings = options(objective, connectivity);
                        auto result = smsd::MCS::find(query).withOptions(settings).match(target);
                        checkCompleted(query, target, settings, result, name, failures);
                    }, failures);
                }
            }
        }

        runCase("strict", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
            check(!settings.chemistry.completeRingsOnly, "strict control must not enable complete rings", failures);
            auto result = smsd::MCS::find(stack).withOptions(settings).match(cyclobutane);
            checkUnproved(stack, cyclobutane, settings, result, "strict", failures);
        }, failures);
        runCase("permissive", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::PERMISSIVE;
            auto result = smsd::MCS::find(stack).withOptions(settings).match(cyclobutane);
            checkCompleted(stack, cyclobutane, settings, result, "permissive", failures);
        }, failures);
        runCase("ordinary rings", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.completeRingsOnly = true;
            auto result = smsd::MCS::find(cyclobutane).withOptions(settings).match(cyclobutane);
            checkCompleted(cyclobutane, cyclobutane, settings, result, "ordinary rings", failures);
        }, failures);
        runCase("stop priority", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.completeRingsOnly = true;
            settings.timeoutMs = 0;
            auto timedOut = smsd::MCS::find(stack).withOptions(settings).match(cyclobutane);
            settings.cancellation = std::make_shared<smsd::CancellationToken>();
            settings.cancellation->cancel();
            auto cancelled = smsd::MCS::find(stack).withOptions(settings).match(cyclobutane);
            check(timedOut.status == smsd::SearchStatus::TimedOut && !timedOut.optimal && !timedOut.complete,
                  "timeout status must be retained", failures);
            check(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete,
                  "cancellation status must be retained", failures);
        }, failures);
        runCase("threaded batch", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.completeRingsOnly = true;
            auto results = smsd::MCS::find(cyclobutane).withOptions(settings)
                .matchMany(std::vector<smsd::MolGraph>{stack, cyclobutane}, 2);
            check(results.size() == 2, "batch must retain both targets", failures);
            if (results.size() != 2) return;
            checkUnproved(cyclobutane, stack, settings, results[0], "batch/truncated", failures);
            checkCompleted(cyclobutane, cyclobutane, settings, results[1], "batch/ordinary", failures);
        }, failures);
        for (bool strict : {false, true}) {
            runCase(strict ? "strict consensus" : "complete-ring consensus", [&] {
                auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
                settings.chemistry.completeRingsOnly = !strict;
                if (strict) settings.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
                const auto result = smsd::MCS::consensus({cyclobutane, stack}, settings, 1.0);
                check(result.status == smsd::SearchStatus::LimitReached && !result.optimal && !result.complete,
                      "consensus ring proof must remain incomplete", failures);
            }, failures);
        }

        const auto chain = fourRingChain();
        const auto carbon = smsd::parseSMILES("C");
        bool relevantTruncated = false;
        (void)chain.computeRelevantCyclesReportingTruncation(relevantTruncated);
        const bool hortonTruncated = smsd::ringsTruncated(chain);
        check(relevantTruncated, "ring chain must cut the relevant-cycle enumeration", failures);
        check(!hortonTruncated, "ring chain must retain complete Horton perception", failures);
        if (!relevantTruncated || hortonTruncated) return 1;

        for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
            for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
                for (bool reverse : {false, true}) {
                    const auto& query = reverse ? carbon : chain;
                    const auto& target = reverse ? chain : carbon;
                    const std::string name = std::string("relevant/")
                        + (objective == smsd::MCSObjective::MCIS ? "MCIS" : "MCES") + "/"
                        + (connectivity == smsd::MCSConnectivity::Connected ? "Connected" : "Fragments")
                        + (reverse ? "/reverse" : "/forward");
                    runCase(name, [&] {
                        auto settings = options(objective, connectivity);
                        settings.chemistry.completeRingsOnly = true;
                        const auto result = smsd::MCS::find(query).withOptions(settings).match(target);
                        checkUnproved(query, target, settings, result, name, failures);
                        settings.chemistry.completeRingsOnly = false;
                        checkCompleted(query, target, settings,
                            smsd::MCS::find(query).withOptions(settings).match(target),
                            name + "/unconstrained", failures);
                    }, failures);
                }
            }
            runCase("relevant/consensus", [&] {
                auto settings = options(objective, smsd::MCSConnectivity::Connected);
                settings.chemistry.completeRingsOnly = true;
                const auto result = smsd::MCS::consensus({chain}, settings);
                check(result.status == smsd::SearchStatus::LimitReached && !result.optimal && !result.complete,
                      "relevant-cycle consensus must remain unproved", failures);
                check(result.found && result.graph.n == chain.n,
                      "relevant-cycle consensus must retain the complete fragment", failures);
            }, failures);
        }
        runCase("relevant/strict-only", [&] {
            auto settings = options(smsd::MCSObjective::MCIS, smsd::MCSConnectivity::Connected);
            settings.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
            checkCompleted(chain, carbon, settings,
                smsd::MCS::find(chain).withOptions(settings).match(carbon), "relevant/strict-only", failures);
            const auto result = smsd::MCS::consensus({chain}, settings);
            check(result.status == smsd::SearchStatus::Completed && result.complete && !result.optimal,
                  "strict-only consensus uses the complete Horton perception", failures);
        }, failures);
    } catch (const std::exception& error) {
        check(false, std::string("fixture: ") + error.what(), failures);
    }
    return failures ? 1 : 0;
}

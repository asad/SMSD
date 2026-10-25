/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"

#include <atomic>
#include <thread>

using smsd_test::require;

static void same(const smsd::SearchResult& expected, const smsd::SearchResult& actual) {
    require(expected.status == actual.status && expected.complete == actual.complete
        && expected.optimal == actual.optimal && expected.score == actual.score
        && expected.nodesVisited == actual.nodesVisited && expected.matches.size() == actual.matches.size(),
        "captured targets changed search outcomes");
    for (size_t match = 0; match < expected.matches.size(); ++match)
        require(expected.matches[match].atomMapping == actual.matches[match].atomMapping
            && expected.matches[match].bondMapping == actual.matches[match].bondMapping
            && expected.matches[match].fragments == actual.matches[match].fragments,
            "captured targets changed original mappings");
}

static void capturedTargetsKeepValuesAndOrder() {
    auto carbon = smsd::parseSMILES("CC");
    carbon.properties["name"] = "before";
    const auto oxygen = smsd::parseSMILES("O");
    const auto expected = std::vector<smsd::MolGraph>{carbon, oxygen, carbon};
    smsd::batch::detail::MolGraphView inputs;
    inputs.push_back(carbon); inputs.push_back(oxygen); inputs.push_back(carbon);
    smsd::batch::TargetCorpus corpus;
    smsd::batch::detail::CorpusAccess::addCapturedTargets(
        corpus, smsd::batch::detail::TargetStorage::capture(inputs));
    carbon.setBondOrder(0, 1, 2);
    carbon.properties["name"] = "after";
    const auto view = smsd::batch::detail::CorpusAccess::targetView(corpus);
    require(view.size() == 3 && view[0].properties.at("name") == "before"
        && view[2].bondOrder(0, 1) == 1, "later source mutations changed captured targets");
    const auto query = smsd::Substructure::find("CC");
    const auto results = query.matchMany(view, 2);
    for (size_t index = 0; index < expected.size(); ++index) same(query.match(expected[index]), results[index]);
    auto detached = corpus.targets();
    detached[0].setBondOrder(0, 1, 2);
    require(detached[2].bondOrder(0, 1) == 1 && view[0].bondOrder(0, 1) == 1,
        "returned value snapshots must remain independently mutable");
    corpus.prewarm(1);
    require(corpus.screen(oxygen, 1.0, 1) == std::vector<int>{1}, "screening lost original target slots");
}

static void copiesAndMovesKeepViews() {
    smsd::batch::TargetCorpus corpus;
    corpus.addTargets({smsd::parseSMILES("CC"), smsd::parseSMILES("O")});
    const auto view = smsd::batch::detail::CorpusAccess::targetView(corpus);
    auto copy = corpus;
    copy.addTarget(smsd::parseSMILES("N"));
    require(copy.size() == 3 && corpus.size() == 2 && view.size() == 2, "an addition changed another snapshot");
    copy.prewarm(1);
    auto moved = std::move(copy);
    copy.addTarget(smsd::parseSMILES("F"));
    require(copy.size() == 1 && moved.size() == 3, "a moved corpus must remain reusable");
    smsd::batch::TargetCorpus assigned;
    assigned = moved;
    moved.addTarget(smsd::parseSMILES("Cl"));
    require(assigned.size() == 3 && moved.size() == 4, "assignment did not preserve independent slots");
    smsd::batch::TargetCorpus moveAssigned;
    moveAssigned = std::move(moved);
    moved.addTarget(smsd::parseSMILES("Br"));
    require(moveAssigned.size() == 4 && moved.size() == 1, "move assignment must leave a reusable corpus");
    const auto query = smsd::Substructure::find("CC");
    require(query.matchMany(view, 2)[0].matches.size() == 1, "a retained view did not survive additions and moves");
}

static void concurrentViewSurvivesAdditions() {
    smsd::batch::TargetCorpus corpus;
    corpus.addTargets({smsd::parseSMILES("c1ccc(CCO)cc1"), smsd::parseSMILES("OCCO")});
    const auto view = smsd::batch::detail::CorpusAccess::targetView(corpus);
    const auto query = smsd::Substructure::smarts("[C;$(CO)]");
    const auto first = query.match(view[0]), second = query.match(view[1]);
    std::atomic<bool> started{false};
    std::exception_ptr error;
    std::thread reader([&] {
        started = true;
        try {
            for (int iteration = 0; iteration < 20; ++iteration) {
                same(first, query.match(view[0]));
                same(second, query.match(view[1]));
            }
        } catch (...) { error = std::current_exception(); }
    });
    while (!started.load()) std::this_thread::yield();
    try { corpus.addTarget(smsd::parseSMILES("N")); }
    catch (...) { reader.join(); throw; }
    reader.join();
    if (error) std::rethrow_exception(error);
    require(view.size() == 2 && corpus.size() == 3, "a concurrent addition changed a retained view");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"capture", capturedTargetsKeepValuesAndOrder}, {"copies", copiesAndMovesKeepViews},
        {"concurrent", concurrentViewSurvivesAdditions}});
}

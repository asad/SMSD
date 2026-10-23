/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"

#include <atomic>
#include <cstdio>
#include <future>

namespace {
using smsd_test::require;

struct Fingerprints {
    int64_t targets;
    int32_t words;
    std::vector<uint64_t> bits;
    std::vector<int32_t> counts;
};

std::atomic<int> uploads{0}, calls{0}, releases{0}, minimumCalls{0};
uint64_t minimumTargets = 0;
int failure = 0;

int deviceInfo(char* text, size_t size, char*, size_t) {
    std::snprintf(text, size, "%s", "Test device");
    return 0;
}

int legacyOverlaps(const int32_t*, int32_t, int64_t, const uint64_t*, const smsd_label_count*, int32_t*,
                    char*, size_t) { return 1; }

int minimumSize(uint32_t operation, uint64_t* minimum, char*, size_t) {
    ++minimumCalls;
    *minimum = operation == 3 ? minimumTargets : UINT64_MAX;
    return 0;
}

int createFingerprints(int64_t targets, int32_t words, const uint64_t* bits, const int32_t* counts,
                       void** output, char* error, size_t size) {
    *output = nullptr;
    if (failure == 1) { std::snprintf(error, size, "%s", "allocation failed"); return 1; }
    auto value = std::make_unique<Fingerprints>();
    value->targets = targets;
    value->words = words;
    value->bits.assign(bits, bits + targets * words);
    value->counts.assign(counts, counts + targets);
    *output = value.release();
    ++uploads;
    return 0;
}

int screenFingerprints(void* handle, int32_t queries, const uint64_t* words, const int32_t* cutoffs,
                        int32_t stride, smsd_fingerprint_hit* hits, uint64_t capacity,
                        uint64_t* offsets, char* error, size_t size) {
    ++calls;
    if (failure == 2) { std::snprintf(error, size, "%s", "device lost"); return 1; }
    const auto& held = *static_cast<const Fingerprints*>(handle);
    uint64_t selected = 0;
    offsets[0] = 0;
    for (int query = 0; query < queries; ++query) {
        for (int64_t target = 0; target < held.targets; ++target) {
            int common = 0;
            for (int word = 0; word < held.words; ++word)
                common += smsd::popcount64(words[query * held.words + word] & held.bits[word * held.targets + target]);
            if (common >= cutoffs[query * stride + held.counts[target]]) {
                if (selected == capacity) return 1;
                hits[selected++] = {static_cast<int32_t>(target), common};
            }
        }
        offsets[query + 1] = selected;
    }
    if (failure == 3 && selected != 0) hits[0].common_bits = 10000;
    return 0;
}

int destroyResident(void* handle, char*, size_t) {
    delete static_cast<Fingerprints*>(handle);
    ++releases;
    return 0;
}

smsd_screen_accelerator_v1 table() {
    return {SMSD_SCREEN_ACCELERATOR_ABI, sizeof(smsd_screen_accelerator_v1), "test", "Test",
        deviceInfo, legacyOverlaps, minimumSize, createFingerprints, screenFingerprints,
        nullptr, nullptr, destroyResident};
}

void sameResults(const std::vector<smsd::SearchResult>& expected, const std::vector<smsd::SearchResult>& actual) {
    require(expected.size() == actual.size(), "one result remains for every target");
    for (size_t target = 0; target < expected.size(); ++target) {
        const auto& first = expected[target];
        const auto& second = actual[target];
        require(first.status == second.status && first.optimal == second.optimal && first.complete == second.complete
                && first.score == second.score && first.matches.size() == second.matches.size(), "result proof and score agree");
        for (size_t match = 0; match < first.matches.size(); ++match)
            require(first.matches[match].atomMapping == second.matches[match].atomMapping
                && first.matches[match].bondMapping == second.matches[match].bondMapping
                && first.matches[match].fragments == second.matches[match].fragments, "matches retain their original indices");
    }
}

void ordinaryMoleculesAndChemistryKeepMatches() {
    const auto accelerator = table();
    require(smsd::screen_accelerator::attach(&accelerator, "test").empty(), "provider attaches");
    minimumTargets = 0;
    failure = 0;
    {
        const std::vector<smsd::MolGraph> graphs{smsd::parseSMILES("CCO"), smsd::parseSMILES("CCN"),
            smsd::parseSMILES("CCCO"), smsd::parseSMILES("COC"), smsd::parseSMILES("C1CC1"),
            smsd::parseSMILES("C.O"), smsd::parseSMILES("[H]"), smsd::MolGraph{}};
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(graphs);
        for (const char* smiles : {"C", "CCO", "C.O", "C1CC1"}) {
            for (bool induced : {false, true}) {
                smsd::Substructure::Options options;
                options.induced = induced;
                const auto query = smsd::Substructure::find(smiles).withOptions(options);
                const auto survivors = query.corpusCandidates(corpus);
                require(survivors.has_value(), "typed molecule queries use the subset screen");
                sameResults(query.matchMany(graphs, 1), query.matchManyScreened(graphs, *survivors, 2));
            }
        }
        smsd::Substructure::Options strict;
        strict.chemistry = smsd::ChemOptions::profile("strict");
        strict.chemistry.useChirality = true;
        const auto ringQuery = smsd::Substructure::find("C1CC1").withOptions(strict);
        sameResults(ringQuery.matchMany(graphs), ringQuery.matchManyScreened(graphs, *ringQuery.corpusCandidates(corpus)));
        require((*smsd::Substructure::find("C").corpusCandidates(corpus))[0] == 1,
                "a carbon query remains a candidate inside a larger molecule");
    }
    smsd::screen_accelerator::detach();
}

void unsupportedQueriesAndPseudoTargetsStayOnCpu() {
    const auto accelerator = table();
    require(smsd::screen_accelerator::attach(&accelerator, "test").empty(), "provider attaches");
    {
        auto pseudo = smsd::MolGraph::Builder().atomCount(1).atomicNumbers({0}).setNeighbors({{}}).build(false);
        const std::vector<smsd::MolGraph> graphs{smsd::parseSMILES("C"), smsd::parseSMILES("O"), pseudo};
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(graphs);
        const int before = calls;
        require(!smsd::Substructure::smarts("[#6,#8]").corpusCandidates(corpus), "SMARTS predicates bypass the screen");
        require(!smsd::Substructure::find(pseudo).corpusCandidates(corpus), "pseudo query atoms bypass the screen");
        require(!smsd::Substructure::find(smsd::MolGraph{}).corpusCandidates(corpus), "empty paths admit every target");
        require(!smsd::Substructure::find("[H]").corpusCandidates(corpus), "hydrogen-only paths admit every target");
        smsd::Substructure::Options untyped;
        untyped.chemistry.matchAtomType = false;
        require(!smsd::Substructure::find("C").withOptions(untyped).corpusCandidates(corpus), "untyped queries bypass the screen");
        require(calls == before, "unsupported inputs do not call the fingerprint screen");
        const auto query = smsd::Substructure::find("C");
        const auto survivors = query.corpusCandidates(corpus);
        require(survivors && (*survivors)[2] == 1, "pseudo target atoms are always retained");
        sameResults(query.matchMany(graphs), query.matchManyScreened(graphs, *survivors));
        require(smsd::Substructure::find(pseudo).match(pseudo).matches.size() == 1
            && smsd::Substructure::find(pseudo).match(graphs[0]).matches.empty()
            && query.match(pseudo).matches.empty(), "molecule pseudo atoms retain their literal matching contract");
        minimumTargets = UINT64_MAX;
        const int bypassBefore = calls;
        require(!query.corpusCandidates(corpus) && calls == bypassBefore, "unmeasured sizes use the existing CPU route");
        minimumTargets = 0;
    }
    smsd::screen_accelerator::detach();
}

void bondAndStereoConstraintsRemainOnCpu() {
    const auto accelerator = table();
    require(smsd::screen_accelerator::attach(&accelerator, "test").empty(), "provider attaches");
    {
        auto molecule = smsd::parseSMILES("NC=O");
        auto tautomer = smsd::parseSMILES("N=CO");
        molecule.tautomerClass.assign(molecule.n, 0);
        tautomer.tautomerClass.assign(tautomer.n, 0);
        smsd::Substructure::Options options;
        options.chemistry.tautomerAware = true;
        const auto query = smsd::Substructure::find(molecule).withOptions(options);
        const std::vector<smsd::MolGraph> targets{tautomer};
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(targets);
        const auto survivors = query.corpusCandidates(corpus);
        require(survivors && (*survivors)[0], "tautomer bond relaxation remains a candidate");
        require(!query.match(tautomer).matches.empty(), "the CPU accepts the relaxed tautomer bonds");
        sameResults(query.matchMany(targets), query.matchManyScreened(targets, *survivors));
    }
    {
        const std::vector<smsd::MolGraph> targets{smsd::parseSMILES("F[C@](Cl)(Br)I"),
                                                smsd::parseSMILES("F[C@@](Cl)(Br)I")};
        smsd::Substructure::Options options;
        options.chemistry.useChirality = true;
        const auto query = smsd::Substructure::find(targets[0]).withOptions(options);
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(targets);
        const auto survivors = query.corpusCandidates(corpus);
        require(survivors && (*survivors)[0] && (*survivors)[1], "stereochemistry is refined on the CPU");
        sameResults(query.matchMany(targets), query.matchManyScreened(targets, *survivors));
    }
    smsd::screen_accelerator::detach();
}

void capturedQueriesCachesAndFailuresKeepResults() {
    const auto accelerator = table();
    require(smsd::screen_accelerator::attach(&accelerator, "test").empty(), "provider attaches");
    const int uploadedBefore = uploads;
    const int releasedBefore = releases;
    {
        const std::vector<smsd::MolGraph> graphs{smsd::parseSMILES("CCO"), smsd::parseSMILES("NN"), smsd::parseSMILES("CC")};
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(graphs);
        auto original = smsd::parseSMILES("C");
        const auto query = smsd::Substructure::find(original);
        original.atomicNum[0] = 8;
        const auto first = query.corpusCandidates(corpus);
        require(first && (*first)[0] && !(*first)[1], "screening reads the captured native query");
        std::vector<std::future<std::optional<std::vector<uint8_t>>>> pending;
        for (int index = 0; index < 3; ++index)
            pending.push_back(std::async(std::launch::async, [&] { return query.corpusCandidates(corpus); }));
        for (auto& result : pending) require(result.get() == first, "concurrent masks agree");
        require(uploads == uploadedBefore + 1, "the corpus uploads element paths once");
        auto copied = corpus;
        require(query.corpusCandidates(copied) == first && uploads == uploadedBefore + 2, "a copy owns separate device state");
        corpus.prewarm();
        require(query.corpusCandidates(corpus) == first && uploads == uploadedBefore + 3, "prewarm invalidates the path cache");
        for (int selectedFailure : {2, 3}) {
            failure = selectedFailure;
            const auto failuresBefore = smsd::screen_accelerator::failureCount();
            require(!query.corpusCandidates(corpus), "a failed screen leaves matching to the CPU");
            require(smsd::screen_accelerator::failureCount() == failuresBefore + 1, "one failed attempt is counted");
        }
        failure = 0;
        corpus.addTarget(graphs[0]);
        failure = 1;
        const auto failuresBefore = smsd::screen_accelerator::failureCount();
        require(!query.corpusCandidates(corpus), "failed allocation leaves matching to the CPU");
        require(smsd::screen_accelerator::failureCount() == failuresBefore + 1, "allocation failure is counted");
        failure = 0;
        require(query.corpusCandidates(corpus)->size() == 4 && uploads == uploadedBefore + 4,
                "added targets appear after a fresh upload");
    }
    require(releases == releasedBefore + 4, "every resident path table is released once");
    smsd::screen_accelerator::detach();
}

void screenedNegativesRespectStops() {
    const std::vector<smsd::MolGraph> graphs{smsd::parseSMILES("N"), smsd::parseSMILES("NN")};
    const std::vector<uint8_t> rejected(graphs.size(), 0);
    const auto stopped = smsd::Substructure::find("O").withTimeout(0);
    sameResults(stopped.matchMany(graphs), stopped.matchManyScreened(graphs, rejected));
    smsd::Substructure::Options options;
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    const auto cancelled = smsd::Substructure::find("O").withOptions(options);
    sameResults(cancelled.matchMany(graphs), cancelled.matchManyScreened(graphs, rejected));
    const auto completed = smsd::Substructure::find("O").matchManyScreened(graphs, rejected);
    for (const auto& result : completed)
        require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete
            && result.matches.empty() && result.score == 0 && result.nodesVisited == 0, "screen rejection is a completed proof");
}

void nodeLimitedQueriesKeepCpuRoute() {
    const auto accelerator = table();
    require(smsd::screen_accelerator::attach(&accelerator, "test").empty(), "provider attaches");
    {
        const std::vector<smsd::MolGraph> graphs{smsd::parseSMILES("C1CC1.C"), smsd::parseSMILES("CCCC")};
        smsd::batch::TargetCorpus corpus;
        corpus.addTargets(graphs);
        smsd::Substructure::Options options;
        options.nodeLimit = 1;
        const auto query = smsd::Substructure::find("CCCC").withOptions(options);
        const int minimumBefore = minimumCalls, uploadsBefore = uploads, callsBefore = calls;
        require(!query.corpusCandidates(corpus, 1), "finite node limits keep the original CPU search");
        require(minimumCalls == minimumBefore && uploads == uploadsBefore && calls == callsBefore,
                "node-limited calls do not consult or prepare an accelerator");
    }
    smsd::screen_accelerator::detach();
}
}

int main(int argc, char** argv) {
#if defined(_WIN32)
    _putenv_s("SMSD_DISABLE_GPU", "0");
#else
    unsetenv("SMSD_DISABLE_GPU");
#endif
    return smsd_test::runGroups(argc, argv, {
        {"molecules", ordinaryMoleculesAndChemistryKeepMatches},
        {"unsupported", unsupportedQueriesAndPseudoTargetsStayOnCpu},
        {"chemistry", bondAndStereoConstraintsRemainOnCpu},
        {"lifetime", capturedQueriesCachesAndFailuresKeepResults},
        {"stops", screenedNegativesRespectStops},
        {"node_limits", nodeLimitedQueriesKeepCpuRoute},
    });
}

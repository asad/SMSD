/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The screening accelerator slot: a table that computes the label overlaps on the CPU stands in
 * for a device, so that every platform checks that attached, failing, refused and disabled
 * accelerators all leave the hits of TargetCorpus screens exactly as the CPU gives them.
 */
#include "test_support.hpp"
#include "smsd/batch.hpp"
#include "smsd/screen_accelerator.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

using smsd_test::require;

int overlapCalls = 0;
int failWith = 0;  // 0 = succeed, 1 = report a failure, 2 = return an impossible overlap

int deviceInfo(char* text, size_t textSize, char*, size_t) {
    std::snprintf(text, textSize, "%s", "Test device (compute 0.0, 0 MB)");
    return 0;
}

int labelOverlaps(const int32_t* query, int32_t tableSize, int64_t targets, const uint64_t* offsets,
                  const smsd_label_count* entries, int32_t* overlaps, char* error, size_t errorSize) {
    ++overlapCalls;
    if (failWith == 1) {
        std::snprintf(error, errorSize, "%s", "device lost");
        return 1;
    }
    for (int64_t i = 0; i < targets; ++i) {
        int32_t overlap = 0;
        for (uint64_t e = offsets[i]; e < offsets[i + 1]; ++e)
            if (entries[e].label < tableSize) overlap += std::min<int32_t>(entries[e].count, query[entries[e].label]);
        overlaps[i] = failWith == 2 ? overlap + 1000 : overlap;
    }
    return 0;
}

const char* const kVersion = "test";

smsd_screen_accelerator_v1 table() {
    return {SMSD_SCREEN_ACCELERATOR_ABI, sizeof(smsd_screen_accelerator_v1), kVersion, "Test", deviceInfo,
            labelOverlaps};
}

void setGpuDisabled(bool disabled) {
#ifdef _WIN32
    _putenv_s("SMSD_DISABLE_GPU", disabled ? "1" : "0");
#else
    setenv("SMSD_DISABLE_GPU", disabled ? "1" : "0", 1);
#endif
}

/** A few hundred targets, so the scores are formed on several threads where OpenMP is present. */
smsd::batch::TargetCorpus corpus(std::vector<smsd::MolGraph>& graphs) {
    const char* smiles[] = {
        "c1ccccc1", "c1ccc(O)cc1", "CCO", "CC(=O)O", "C1CCCCC1", "c1ccncc1", "O=C(O)c1ccccc1O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "[NH4+]", "[O-]C(=O)C", "C1=CC=CC=C1C2=CC=CC=C2", "ClC(Cl)Cl",
        "OC1C(O)C(O)C(O)C(O)C1O", "c1ccc2ccccc2c1", "C#N", "[Na+].[Cl-]", "FC(F)(F)c1ccccc1", "S",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "C1CC1", "N[C@@H](C)C(=O)O", "c1cc[nH]c1", "B(O)(O)O"};
    for (const char* text : smiles) graphs.push_back(smsd::parseSMILES(text));
    smsd::batch::TargetCorpus result;
    std::vector<smsd::MolGraph> targets;
    for (int copy = 0; copy < 25; ++copy) targets.insert(targets.end(), graphs.begin(), graphs.end());
    // A target too large for a histogram count is scored on its graph, as on the CPU.
    constexpr int large = 70000;
    targets.push_back(smsd::MolGraph::Builder().atomCount(large).atomicNumbers(std::vector<int>(large, 6))
                          .setNeighbors(std::vector<std::vector<int>>(large)).build(false));
    result.addTargets(std::move(targets));
    return result;
}

std::vector<smsd::ChemOptions> chemistries() {
    smsd::ChemOptions typed;
    smsd::ChemOptions strict;
    strict.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    strict.ringMatchesRingOnly = true;
    smsd::ChemOptions untyped;
    untyped.matchAtomType = false;
    return {typed, strict, untyped};
}

/** Every threshold that matters: none, scores the corpus reaches exactly, and none reached. */
std::vector<double> thresholds(const smsd::MolGraph& query, const std::vector<smsd::MolGraph>& graphs,
                               const smsd::ChemOptions& chemistry) {
    std::vector<double> result{0.0, 1.0, 1.0000001};
    for (const auto& target : graphs) result.push_back(smsd::similarityUpperBound(query, target, chemistry));
    return result;
}

/** Screens every query under every chemistry and threshold, attached and detached, and compares the hits. */
void requireSameHits(const smsd_screen_accelerator_v1* attached, bool expectCalls) {
    std::vector<smsd::MolGraph> graphs;
    const smsd::batch::TargetCorpus targets = corpus(graphs);
    for (const auto& chemistry : chemistries()) {
        for (std::size_t q = 0; q < graphs.size(); q += 3) {
            for (double threshold : thresholds(graphs[q], graphs, chemistry)) {
                smsd::screen_accelerator::detach();
                const auto expected = targets.screenTargets(graphs[q], threshold, chemistry, 0);
                require(smsd::screen_accelerator::attach(attached, kVersion).empty(), "the table is attached");
                const int before = overlapCalls;
                const auto hits = targets.screenTargets(graphs[q], threshold, chemistry, 0);
                require(hits == expected, "the accelerated screen gives the CPU hits");
                require((overlapCalls > before) == expectCalls, "the accelerator is used only when it may be");
                smsd::screen_accelerator::detach();
            }
        }
    }
}

void attachedAcceleratorKeepsTheHits() {
    setGpuDisabled(false);
    failWith = 0;
    const auto accelerator = table();
    requireSameHits(&accelerator, true);
    require(smsd::screen_accelerator::attach(&accelerator, kVersion).empty(), "the table is attached");
    require(smsd::screen_accelerator::active() == &accelerator, "the table is active");
    require(smsd::screen_accelerator::deviceInfo() == "Test device (compute 0.0, 0 MB)", "the device is described");
    smsd::screen_accelerator::detach();
    require(smsd::screen_accelerator::active() == nullptr, "a detached table is inactive");
    require(smsd::screen_accelerator::deviceInfo().empty(), "no device without a table");
}

void failingAcceleratorFallsBackToTheCpu() {
    setGpuDisabled(false);
    const auto accelerator = table();
    for (int failure : {1, 2}) {
        failWith = failure;
        const long long before = smsd::screen_accelerator::failureCount();
        requireSameHits(&accelerator, true);
        require(smsd::screen_accelerator::failureCount() > before, "the failures are counted");
        require(!smsd::screen_accelerator::latestFailure().empty(), "the latest failure gives its reason");
    }
    failWith = 1;
    requireSameHits(&accelerator, true);
    require(smsd::screen_accelerator::latestFailure() == "device lost", "the accelerator's reason is kept");
    failWith = 0;
}

void disabledGpuKeepsTheCpu() {
    setGpuDisabled(true);
    failWith = 0;
    const auto accelerator = table();
    requireSameHits(&accelerator, false);
    require(smsd::screen_accelerator::attach(&accelerator, kVersion).empty(), "the table is attached");
    require(smsd::screen_accelerator::active() == nullptr, "SMSD_DISABLE_GPU leaves it inactive");
    smsd::screen_accelerator::detach();
    setGpuDisabled(false);
}

void unsuitableTablesAreRefused() {
    setGpuDisabled(false);
    auto wrongInterface = table();
    wrongInterface.abi_version = SMSD_SCREEN_ACCELERATOR_ABI + 1;
    auto shortTable = table();
    shortTable.struct_size = SMSD_SCREEN_ACCELERATOR_BASE_SIZE - 1;
    auto missingFunction = table();
    missingFunction.label_overlaps = nullptr;
    auto unnamed = table();
    unnamed.name = nullptr;
    auto otherVersion = table();
    otherVersion.smsd_version = "0.0.0";
    auto unknownVersion = table();
    unknownVersion.smsd_version = nullptr;
    const std::vector<const smsd_screen_accelerator_v1*> unsuitable{
        nullptr, &wrongInterface, &shortTable, &missingFunction, &unnamed, &otherVersion, &unknownVersion};
    for (const smsd_screen_accelerator_v1* refused : unsuitable) {
        require(!smsd::screen_accelerator::attach(refused, kVersion).empty(), "an unsuitable table is refused");
        require(smsd::screen_accelerator::active() == nullptr, "a refused table is not attached");
    }
    const auto accelerator = table();
    require(!smsd::screen_accelerator::attach(&accelerator, nullptr).empty(), "an unknown build version is refused");
    require(smsd::screen_accelerator::refusal(&otherVersion, kVersion).find("0.0.0") != std::string::npos,
            "the refusal names the version the table serves");
}

/** A resident provider double, retaining the exact immutable arrays supplied by the corpus. */
struct ResidentTargets {
    int64_t count = 0;
    int words = 0;
    std::vector<uint64_t> fingerprints, offsets;
    std::vector<int32_t> popcounts;
    std::vector<smsd_label_count> labels;
};

std::atomic<int> fingerprintUploads{0}, labelUploads{0}, residentCalls{0}, releases{0};
uint64_t minimumTargets = 0;
int residentFailure = 0;

int minimumSize(uint32_t, uint64_t* minimum, char*, size_t) {
    *minimum = minimumTargets;
    return 0;
}

int createFingerprints(int64_t count, int32_t words, const uint64_t* fingerprints,
                       const int32_t* popcounts, void** result, char* error, size_t errorSize) {
    if (residentFailure == 1) {
        *result = nullptr;
        std::snprintf(error, errorSize, "%s", "allocation failed");
        return 1;
    }
    auto target = std::make_unique<ResidentTargets>();
    target->count = count;
    target->words = words;
    target->fingerprints.assign(fingerprints, fingerprints + count * words);
    target->popcounts.assign(popcounts, popcounts + count);
    *result = target.release();
    ++fingerprintUploads;
    return 0;
}

int screenFingerprints(void* resident, int32_t queries, const uint64_t* words, const int32_t* cutoffs,
                       int32_t stride, smsd_fingerprint_hit* hits, uint64_t capacity, uint64_t* offsets,
                       char* error, size_t errorSize) {
    ++residentCalls;
    if (residentFailure == 2) {
        std::snprintf(error, errorSize, "%s", "device lost");
        return 1;
    }
    const auto& targets = *static_cast<ResidentTargets*>(resident);
    uint64_t selected = 0;
    offsets[0] = 0;
    for (int queryIndex = 0; queryIndex < queries; ++queryIndex) {
        for (int64_t targetIndex = 0; targetIndex < targets.count; ++targetIndex) {
            int common = 0;
            for (int wordIndex = 0; wordIndex < targets.words; ++wordIndex)
                common += smsd::popcount64(words[queryIndex * targets.words + wordIndex]
                    & targets.fingerprints[wordIndex * targets.count + targetIndex]);
            if (common < cutoffs[queryIndex * stride + targets.popcounts[targetIndex]]) continue;
            if (selected >= capacity) return 1;
            hits[selected++] = {static_cast<int32_t>(targetIndex), common};
        }
        offsets[queryIndex + 1] = selected;
    }
    if (selected && residentFailure == 3) hits[0].common_bits = 10000;
    if (selected && residentFailure == 4) hits[0].target_index = static_cast<int32_t>(targets.count);
    return 0;
}

int createLabels(int64_t count, const uint64_t* offsets, const smsd_label_count* entries,
                 void** result, char* error, size_t errorSize) {
    if (residentFailure == 1) {
        *result = nullptr;
        std::snprintf(error, errorSize, "%s", "allocation failed");
        return 1;
    }
    auto target = std::make_unique<ResidentTargets>();
    target->count = count;
    target->offsets.assign(offsets, offsets + count + 1);
    if (offsets[count]) target->labels.assign(entries, entries + offsets[count]);
    *result = target.release();
    ++labelUploads;
    return 0;
}

int residentLabelOverlaps(void* resident, int32_t queries, const int32_t* counts, int32_t tableSize,
                          int32_t* overlaps, char* error, size_t errorSize) {
    ++residentCalls;
    if (residentFailure == 2) {
        std::snprintf(error, errorSize, "%s", "device lost");
        return 1;
    }
    const auto& targets = *static_cast<ResidentTargets*>(resident);
    for (int queryIndex = 0; queryIndex < queries; ++queryIndex) {
        for (int64_t targetIndex = 0; targetIndex < targets.count; ++targetIndex) {
            int32_t overlap = 0;
            for (uint64_t entryIndex = targets.offsets[targetIndex];
                 entryIndex < targets.offsets[targetIndex + 1]; ++entryIndex) {
                const auto& entry = targets.labels[entryIndex];
                if (entry.label < tableSize)
                    overlap += std::min<int32_t>(entry.count, counts[queryIndex * tableSize + entry.label]);
            }
            overlaps[queryIndex * targets.count + targetIndex] = residentFailure == 3 ? -1 : overlap;
        }
    }
    return 0;
}

int destroyResident(void* resident, char*, size_t) {
    delete static_cast<ResidentTargets*>(resident);
    ++releases;
    return 0;
}

smsd_screen_accelerator_v1 residentTable() {
    auto result = table();
    result.minimum_targets = minimumSize;
    result.create_fingerprints = createFingerprints;
    result.fingerprint_screen = screenFingerprints;
    result.create_labels = createLabels;
    result.resident_label_overlaps = residentLabelOverlaps;
    result.destroy_resident = destroyResident;
    return result;
}

std::vector<smsd::MolGraph> smallGraphs() {
    return {smsd::MolGraph{}, smsd::parseSMILES("C"), smsd::parseSMILES("CCO"),
            smsd::parseSMILES("CC(=O)O"), smsd::parseSMILES("c1ccccc1"), smsd::parseSMILES("c1ccncc1")};
}

void residentScreensKeepBoundaryHits() {
    setGpuDisabled(false);
    minimumTargets = 0;
    residentFailure = 0;
    auto accelerator = residentTable();
    const auto graphs = smallGraphs();
    smsd::batch::TargetCorpus targets;
    targets.addTargets(graphs);
    targets.prewarm();
    for (const auto& chemistry : chemistries()) {
        for (const auto& query : graphs) {
            const auto queryFP = smsd::batch::detail::computeCircularFingerprintECFP(query, 2, 2048);
            std::vector<double> boundaries{0.0, 1.0, -std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()};
            for (const auto& graph : graphs) {
                const auto targetFP = smsd::batch::detail::computeCircularFingerprintECFP(graph, 2, 2048);
                const double score = smsd::batch::fingerprintTanimoto(queryFP, targetFP);
                boundaries.push_back(score);
                boundaries.push_back(std::nextafter(score, -std::numeric_limits<double>::infinity()));
                boundaries.push_back(std::nextafter(score, std::numeric_limits<double>::infinity()));
            }
            for (double threshold : boundaries) {
                setGpuDisabled(true);
                const auto expected = targets.screen(query, threshold);
                const auto expectedLabels = targets.screenTargets(query, threshold, chemistry);
                setGpuDisabled(false);
                require(smsd::screen_accelerator::attach(&accelerator, kVersion).empty(), "resident table attaches");
                require(targets.screen(query, threshold) == expected, "resident fingerprints keep boundary hits");
                require(targets.screenTargets(query, threshold, chemistry) == expectedLabels,
                        "resident labels keep boundary hits");
            }
        }
    }
    const auto batch = targets.screenMany(graphs, 0.5);
    const auto labelBatch = targets.screenTargetsMany(graphs, 0.5, smsd::ChemOptions{});
    require(batch.size() == graphs.size() && labelBatch.size() == graphs.size(), "batch keeps query order");
    for (std::size_t queryIndex = 0; queryIndex < graphs.size(); ++queryIndex) {
        require(batch[queryIndex] == targets.screen(graphs[queryIndex], 0.5), "fingerprint batch equals singles");
        require(labelBatch[queryIndex] == targets.screenTargets(graphs[queryIndex], 0.5, smsd::ChemOptions{}),
                "label batch equals singles");
    }
    smsd::screen_accelerator::detach();
}

void residentFailuresAndSizeBypassKeepCpuHits() {
    auto accelerator = residentTable();
    const auto graphs = smallGraphs();
    smsd::batch::TargetCorpus targets;
    targets.addTargets(graphs);
    targets.prewarm();
    setGpuDisabled(true);
    const auto expected = targets.screen(graphs[2], 0.0);
    const auto expectedLabels = targets.screenTargets(graphs[2], 0.0, smsd::ChemOptions{});
    setGpuDisabled(false);
    require(smsd::screen_accelerator::attach(&accelerator, kVersion).empty(), "resident table attaches");
    minimumTargets = std::numeric_limits<uint64_t>::max();
    const int before = residentCalls;
    const int legacyBefore = overlapCalls;
    require(targets.screen(graphs[2], 0.0) == expected, "small fingerprint corpus uses CPU");
    require(targets.screenTargets(graphs[2], 0.0, smsd::ChemOptions{}) == expectedLabels,
            "small label corpus uses CPU");
    require(residentCalls == before && overlapCalls == legacyBefore, "size bypass does not invoke a provider");
    minimumTargets = 0;
    for (int failure : {1, 2, 3, 4}) {
        residentFailure = failure;
        const auto failuresBefore = smsd::screen_accelerator::failureCount();
        require(targets.screen(graphs[2], 0.0) == expected, "failed resident fingerprint screen uses CPU hits");
        require(smsd::screen_accelerator::failureCount() == failuresBefore + 1, "one failed attempt is counted");
    }
    residentFailure = 0;
    require(targets.screen(graphs[2], 0.0) == expected, "resident screen recovers after failures");
    for (int failure : {1, 2, 3}) {
        residentFailure = failure;
        const auto failuresBefore = smsd::screen_accelerator::failureCount();
        require(targets.screenTargets(graphs[2], 0.0, smsd::ChemOptions{}) == expectedLabels,
                "failed resident label screen uses CPU hits");
        require(smsd::screen_accelerator::failureCount() == failuresBefore + 1, "one failed label attempt is counted");
    }
    residentFailure = 0;
    smsd::screen_accelerator::detach();
}

void residentStorageFollowsCorpusLifetime() {
    auto accelerator = residentTable();
    const auto graphs = smallGraphs();
    require(smsd::screen_accelerator::attach(&accelerator, kVersion).empty(), "resident table attaches");
    const int uploadedBefore = fingerprintUploads;
    const int releasedBefore = releases;
    {
        smsd::batch::TargetCorpus targets;
        targets.addTargets(graphs);
        targets.prewarm();
        const auto expected = targets.screen(graphs[2], 0.5);
        std::vector<std::future<std::vector<int>>> screens;
        for (int queryIndex = 0; queryIndex < 4; ++queryIndex)
            screens.push_back(std::async(std::launch::async, [&] { return targets.screen(graphs[2], 0.5); }));
        for (auto& screen : screens) require(screen.get() == expected, "concurrent resident reuse keeps hits");
        require(fingerprintUploads == uploadedBefore + 1, "repeated queries upload one immutable corpus");
        auto copied = targets;
        require(copied.screen(graphs[2], 0.5) == expected, "copied corpus keeps hits");
        require(fingerprintUploads == uploadedBefore + 2, "copied corpus owns distinct device storage");
        targets.prewarm();
        require(targets.screen(graphs[2], 0.5) == expected, "repeated prewarm keeps hits");
        require(fingerprintUploads == uploadedBefore + 3, "prewarm replaces resident fingerprints");
        targets.addTarget(graphs[2]);
        targets.prewarm();
        auto extended = expected;
        extended.push_back(static_cast<int>(graphs.size()));
        require(targets.screen(graphs[2], 0.5) == extended, "added targets are screened in their original order");
        require(fingerprintUploads == uploadedBefore + 4, "addition replaces resident fingerprints");
        auto moved = std::move(targets);
        require(moved.screen(graphs[2], 0.5) == extended, "moving the corpus keeps its targets and hits");
    }
    require(releases - releasedBefore == fingerprintUploads - uploadedBefore,
            "all corpus-owned allocations are released once");
    smsd::screen_accelerator::detach();
}

void oldAndPartialTablesKeepTheOriginalContract() {
    auto accelerator = residentTable();
    for (std::size_t bytes : {static_cast<std::size_t>(SMSD_SCREEN_ACCELERATOR_BASE_SIZE),
                             offsetof(smsd_screen_accelerator_v1, create_fingerprints),
                             offsetof(smsd_screen_accelerator_v1, destroy_resident)}) {
        accelerator.struct_size = static_cast<uint32_t>(bytes);
        std::unique_ptr<void, decltype(&std::free)> allocation(std::malloc(bytes), std::free);
        require(allocation != nullptr, "provider prefix is allocated");
        std::memcpy(allocation.get(), &accelerator, bytes);
        const auto* prefix = static_cast<const smsd_screen_accelerator_v1*>(allocation.get());
        require(smsd::screen_accelerator::attach(prefix, kVersion).empty(), "an existing ABI prefix is accepted");
        const auto graphs = smallGraphs();
        smsd::batch::TargetCorpus targets;
        targets.addTargets(graphs);
        const auto expected = smsd::screenTargets(graphs[2], graphs, 0.5, smsd::ChemOptions{});
        const int before = overlapCalls;
        require(targets.screenTargets(graphs[2], 0.5, smsd::ChemOptions{}) == expected,
                "a partial extension preserves the original label callback");
        require(overlapCalls == before + 1, "only the covered original callback runs");
        smsd::screen_accelerator::detach();
    }
}

void batchesCrossWorkspaceChunksAndReplaceProviders() {
    auto first = residentTable();
    auto second = residentTable();
    const auto graphs = smallGraphs();
    smsd::batch::TargetCorpus targets;
    targets.addTargets(graphs);
    targets.prewarm();
    std::vector<smsd::MolGraph> queries;
    for (std::size_t queryIndex = 0; queryIndex < 65; ++queryIndex)
        queries.push_back(graphs[queryIndex % graphs.size()]);
    require(smsd::screen_accelerator::attach(&first, kVersion).empty(), "the first provider attaches");
    const auto hits = targets.screenMany(queries, 0.5);
    const auto labels = targets.screenTargetsMany(queries, 0.5, smsd::ChemOptions{});
    for (std::size_t queryIndex = 0; queryIndex < queries.size(); ++queryIndex) {
        require(hits[queryIndex] == targets.screen(queries[queryIndex], 0.5), "fingerprint chunk keeps order");
        require(labels[queryIndex] == targets.screenTargets(queries[queryIndex], 0.5, smsd::ChemOptions{}),
                "label chunk keeps order");
    }
    const int before = releases;
    require(smsd::screen_accelerator::attach(&second, kVersion).empty(), "the replacement provider attaches");
    require(targets.screenMany(queries, 0.5) == hits, "replacement fingerprint provider keeps hits");
    require(targets.screenTargetsMany(queries, 0.5, smsd::ChemOptions{}) == labels,
            "replacement label provider keeps hits");
    require(releases == before + 2, "replacement frees each handle through its original provider");
    smsd::screen_accelerator::detach();
}

void fingerprintBatchesKeepExactThresholdsAndOrder() {
    auto accelerator = residentTable();
    residentFailure = 0;
    for (std::size_t words : {std::size_t(1), std::size_t(32), std::size_t(1025)}) {
        std::vector<std::vector<uint64_t>> targets(257, std::vector<uint64_t>(words, 0));
        for (std::size_t target = 0; target < targets.size(); ++target) {
            if (target % 5 == 0) continue;
            for (std::size_t word = 0; word < words; ++word)
                targets[target][word] = target % 5 == 1 ? UINT64_MAX
                    : (uint64_t(1) << ((target + word) % 64));
        }
        std::vector<std::vector<uint64_t>> queries;
        for (std::size_t query = 0; query < (words > 32 ? 5 : 65); ++query)
            queries.push_back(targets[query % 5]);
        smsd::screen_accelerator::ResidentCache cache;
        const auto score = [](int common, int queryBits, int targetBits) {
            return smsd::batch::tanimoto(smsd::batch::BitCounts{common, queryBits, targetBits});
        };
        for (double threshold : {-1.0, 0.0, std::nextafter(0.5, 0.0), 0.5,
                                  std::nextafter(0.5, 1.0), 1.0, std::nextafter(1.0, 2.0),
                                  std::numeric_limits<double>::quiet_NaN()}) {
            const auto actual = smsd::screen_accelerator::fingerprintScreens(
                queries, targets, threshold, cache, score, &accelerator);
            require(actual.has_value() && actual->size() == queries.size(), "prepared batch retains every query");
            for (std::size_t query = 0; query < queries.size(); ++query) {
                std::vector<int> expected;
                for (std::size_t target = 0; target < targets.size(); ++target)
                    if (smsd::batch::fingerprintTanimoto(queries[query], targets[target]) >= threshold)
                        expected.push_back(static_cast<int>(target));
                require((*actual)[query] == expected, "integer thresholds keep CPU boundary hits in target order");
            }
        }
    }
}


/** Decimal boundaries keep the CPU division and equal-score target order. */
void fingerprintBatchesKeepRoundedDecimalThresholds() {
    auto accelerator = residentTable();
    residentFailure = 0;
    const std::vector<std::vector<uint64_t>> queries{{0x1f}};
    const std::vector<std::vector<uint64_t>> targets{{0x3e1}, {0x3e1}, {0x1f}, {0}};
    smsd::screen_accelerator::ResidentCache cache;
    const auto score = [](int common, int queryBits, int targetBits) {
        return smsd::batch::tanimoto(smsd::batch::BitCounts{common, queryBits, targetBits});
    };
    const double boundary = smsd::batch::fingerprintTanimoto(queries[0], targets[0]);
    require(boundary == 0.1, "decimal fixture has one common bit out of ten union bits");
    for (double threshold : {std::nextafter(boundary, 0.0), boundary,
                              std::nextafter(boundary, 1.0)}) {
        const auto actual = smsd::screen_accelerator::fingerprintScreens(
            queries, targets, threshold, cache, score, &accelerator);
        require(actual.has_value() && actual->size() == 1, "decimal screen is retained");
        const std::vector<int> expected = threshold <= boundary
            ? std::vector<int>{0, 1, 2} : std::vector<int>{2};
        require((*actual)[0] == expected, "decimal boundary and tied target indices follow the CPU");
    }
}

/** Large counts must not overflow before forming their union. */
void labelSimilaritySupportsFullIntCounts() {
    const int maximum = std::numeric_limits<int>::max();
    require(smsd::screening_detail::overlapSimilarity(maximum, maximum, maximum) == 1.0,
            "equal large atom counts have unit similarity");
    require(smsd::screening_detail::overlapSimilarity(1, maximum, 1) == 1.0 / maximum,
            "a large atom-count sum retains its denominator");
}

} // namespace

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"attached", attachedAcceleratorKeepsTheHits},
        {"failing", failingAcceleratorFallsBackToTheCpu},
        {"disabled", disabledGpuKeepsTheCpu},
        {"refused", unsuitableTablesAreRefused},
        {"resident-boundaries", residentScreensKeepBoundaryHits},
        {"resident-failures", residentFailuresAndSizeBypassKeepCpuHits},
        {"resident-lifetime", residentStorageFollowsCorpusLifetime},
        {"old-prefix", oldAndPartialTablesKeepTheOriginalContract},
        {"resident-chunks", batchesCrossWorkspaceChunksAndReplaceProviders},
        {"fingerprint-batch-thresholds", fingerprintBatchesKeepExactThresholdsAndOrder},
        {"fingerprint-decimal-thresholds", fingerprintBatchesKeepRoundedDecimalThresholds},
        {"label-count-overflow", labelSimilaritySupportsFullIntCounts},
    });
}

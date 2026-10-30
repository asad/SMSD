/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The screening accelerator slot: a separately installed package, such as smsd-cuda, attaches
 * its table of functions (screen_accelerator_abi.h) at run time, and TargetCorpus screens then
 * have it compute the label overlaps of the targets. Nothing is attached unless a program, or the
 * smsd Python package on finding smsd-cuda, attaches a table, and SMSD_DISABLE_GPU set to a value
 * other than 0 keeps every screen on the CPU. The accelerator returns exact integer overlaps and
 * every score is formed from them here by overlapSimilarity, as on the CPU, so an accelerated
 * screen gives the CPU scores bit for bit and therefore the same hits. When the accelerator
 * fails, or returns an overlap no target can have, the screen runs on the CPU.
 */
#pragma once
#ifndef SMSD_SCREEN_ACCELERATOR_HPP
#define SMSD_SCREEN_ACCELERATOR_HPP

#include "smsd/screen_accelerator_abi.h"
#include "smsd/gpu_kernels.hpp"
#include "smsd/mol_graph.hpp"
#include "smsd/similarity_screen.hpp"

#include <atomic>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace smsd {
namespace screen_accelerator {

namespace detail {

inline std::atomic<const smsd_screen_accelerator_v1*>& slot() noexcept {
    static std::atomic<const smsd_screen_accelerator_v1*> table{nullptr};
    return table;
}

/** The accelerated screens that ran on the CPU instead, and the reason for the latest. */
struct Failures {
    std::mutex mutex;
    std::string latest;
    long long count = 0;
};

inline Failures& failures() {
    static Failures record;
    return record;
}

inline void recordFailure(const char* reason) {
    Failures& record = failures();
    const std::lock_guard<std::mutex> lock(record.mutex);
    record.latest = reason != nullptr && reason[0] != '\0' ? reason : "the accelerator failed without a reason";
    ++record.count;
}

} // namespace detail

/** Why the table cannot be attached to an SMSD build of the given version, or the empty string when it can. */
inline std::string refusal(const smsd_screen_accelerator_v1* table, const char* smsdVersion) {
    if (table == nullptr) return "no accelerator table was given";
    if (table->abi_version != SMSD_SCREEN_ACCELERATOR_ABI)
        return "the accelerator implements interface version " + std::to_string(table->abi_version)
             + ", and this build of smsd needs version " + std::to_string(SMSD_SCREEN_ACCELERATOR_ABI);
    if (table->struct_size < SMSD_SCREEN_ACCELERATOR_BASE_SIZE) return "the accelerator table is incomplete";
    if (table->name == nullptr || table->device_info == nullptr || table->label_overlaps == nullptr)
        return "the accelerator table is incomplete";
    const char* builtFor = table->smsd_version != nullptr ? table->smsd_version : "an unknown version";
    if (smsdVersion == nullptr || std::strcmp(builtFor, smsdVersion) != 0)
        return std::string("the accelerator was built for smsd ") + builtFor + ", not smsd "
             + (smsdVersion != nullptr ? smsdVersion : "of an unknown version");
    return {};
}

/**
 * Attaches the table in place of any attached before, unless refusal() gives a reason against
 * it; returns that reason, empty when the table was attached. The table must stay valid until
 * every resident handle it created has been destroyed, including after it is detached.
 */
inline std::string attach(const smsd_screen_accelerator_v1* table, const char* smsdVersion) {
    std::string reason = refusal(table, smsdVersion);
    if (reason.empty()) detail::slot().store(table, std::memory_order_release);
    return reason;
}

/** Detaches the attached table, so that screens run on the CPU. */
inline void detach() noexcept { detail::slot().store(nullptr, std::memory_order_release); }

/** The attached table, or null when none is attached or SMSD_DISABLE_GPU disables the GPU. */
inline const smsd_screen_accelerator_v1* active() noexcept {
    if (gpu_kern::detail::disabledByEnvironment()) return nullptr;
    return detail::slot().load(std::memory_order_acquire);
}

/** The number of accelerated screens that ran on the CPU instead. */
inline long long failureCount() {
    detail::Failures& record = detail::failures();
    const std::lock_guard<std::mutex> lock(record.mutex);
    return record.count;
}

/** Why the latest accelerated screen ran on the CPU instead; empty when none has. */
inline std::string latestFailure() {
    detail::Failures& record = detail::failures();
    const std::lock_guard<std::mutex> lock(record.mutex);
    return record.latest;
}

/** The active accelerator's description of its device; empty when none is active or it cannot say. */
inline std::string deviceInfo() {
    const smsd_screen_accelerator_v1* table = active();
    if (table == nullptr) return {};
    char text[256] = {};
    char error[256] = {};
    if (table->device_info(text, sizeof text, error, sizeof error) != 0) {
        error[sizeof error - 1] = '\0';
        detail::recordFailure(error);
        return {};
    }
    text[sizeof text - 1] = '\0';
    return text;
}

/** An optional callback group is read only when its complete storage is present. */
inline bool hasResidentTail(const smsd_screen_accelerator_v1* table) {
    return table != nullptr
        && table->struct_size >= offsetof(smsd_screen_accelerator_v1, destroy_resident)
                                     + sizeof(table->destroy_resident)
        && table->minimum_targets != nullptr && table->destroy_resident != nullptr;
}

inline bool hasResidentFingerprints(const smsd_screen_accelerator_v1* table) {
    return hasResidentTail(table) && table->create_fingerprints != nullptr && table->fingerprint_screen != nullptr;
}

inline bool hasResidentLabels(const smsd_screen_accelerator_v1* table) {
    return hasResidentTail(table) && table->create_labels != nullptr && table->resident_label_overlaps != nullptr;
}

/** Corpus-owned lazy device state. Copies and moves start without a device allocation. */
struct ResidentCache {
    std::mutex mutex;
    const smsd_screen_accelerator_v1* table = nullptr;
    void* resident = nullptr;
    std::vector<std::int32_t> popcounts;

    ResidentCache() = default;
    ResidentCache(const ResidentCache&) noexcept {}
    ResidentCache& operator=(const ResidentCache&) noexcept { clear(); return *this; }
    ~ResidentCache() { clear(); }

    /** The corpus lock excludes calls while its targets or fingerprints are replaced. */
    void clear() noexcept {
        if (resident != nullptr) {
            char error[256] = {};
            if (table->destroy_resident(resident, error, sizeof error) != 0) {
                error[sizeof error - 1] = '\0';
                try { detail::recordFailure(error); } catch (...) {}
            }
        }
        resident = nullptr;
        table = nullptr;
        popcounts.clear();
    }
};

/** Whether the provider has a measured crossover for this operation and target count. */
inline bool reachesMinimum(const smsd_screen_accelerator_v1* table, std::uint32_t operation,
                            std::size_t targets) {
    std::uint64_t minimum = UINT64_MAX;
    char error[256] = {};
    if (table->minimum_targets(operation, &minimum, error, sizeof error) != 0) {
        error[sizeof error - 1] = '\0';
        detail::recordFailure(error);
        return false;
    }
    return minimum != UINT64_MAX && targets >= minimum;
}

/** Bound temporary storage independently of the number of queries in a batch. */
inline std::size_t queriesPerChunk(std::size_t targets, std::size_t queries) {
    constexpr std::size_t maximumPairs = std::size_t(1) << 23;
    constexpr std::size_t maximumQueries = 64;
    return std::min(std::min(queries, maximumQueries), std::max<std::size_t>(1, maximumPairs / targets));
}

/**
 * Resident fingerprint screens. Score forms a score from (common, query count, target count)
 * through the caller's CPU expression, including while the integer threshold table is built.
 */
template<class Score>
inline std::optional<std::vector<std::vector<int>>> fingerprintScreens(
    const std::vector<std::vector<std::uint64_t>>& queries,
    const std::vector<std::vector<std::uint64_t>>& targets, double threshold,
    ResidentCache& cache, const Score& score,
    const smsd_screen_accelerator_v1* provider = nullptr) {
    const smsd_screen_accelerator_v1* table = provider != nullptr ? provider : active();
    if (!hasResidentFingerprints(table) || targets.empty() || queries.empty()
        || targets.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || (provider == nullptr && !reachesMinimum(table, 1, targets.size()))) return std::nullopt;
    const std::size_t words = targets.front().size();
    if (words == 0 || words > static_cast<std::size_t>(std::numeric_limits<int>::max() / 64)
        || targets.size() > std::numeric_limits<std::size_t>::max() / words) return std::nullopt;
    for (const auto& query : queries) if (query.size() != words) return std::nullopt;
    const int bits = static_cast<int>(words * 64);
    const std::size_t targetCount = targets.size();
    const std::lock_guard<std::mutex> lock(cache.mutex);
    if (cache.table != table) cache.clear();
    char error[256] = {};
    auto failed = [&]() {
        error[sizeof error - 1] = '\0';
        detail::recordFailure(error);
    };
    if (cache.resident == nullptr) {
        for (const auto& target : targets) if (target.size() != words) return std::nullopt;
        std::vector<std::uint64_t> wordMajor(words * targetCount);
        std::vector<std::int32_t> counts(targetCount, 0);
        for (std::size_t target = 0; target < targetCount; ++target) {
            for (std::size_t word = 0; word < words; ++word)
                wordMajor[word * targetCount + target] = targets[target][word];
            counts[target] = static_cast<std::int32_t>(smsd::popcountWords(targets[target].data(), words));
        }
        void* resident = nullptr;
        if (table->create_fingerprints(static_cast<std::int64_t>(targetCount), static_cast<std::int32_t>(words),
                                        wordMajor.data(), counts.data(), &resident, error, sizeof error) != 0) {
            failed();
            return std::nullopt;
        }
        if (resident == nullptr) {
            detail::recordFailure("the accelerator returned no resident fingerprint table");
            return std::nullopt;
        }
        cache.table = table;
        cache.resident = resident;
        cache.popcounts = std::move(counts);
    }
    std::vector<std::vector<int>> results(queries.size());
    const std::size_t chunkSize = queriesPerChunk(targetCount, queries.size());
    for (std::size_t begin = 0; begin < queries.size(); begin += chunkSize) {
        const std::size_t count = std::min(chunkSize, queries.size() - begin);
        std::vector<std::uint64_t> queryWords(count * words);
        std::vector<std::int32_t> queryCounts(count, 0);
        const std::size_t stride = static_cast<std::size_t>(bits) + 1;
        std::vector<std::int32_t> cutoffs(count * stride);
        for (std::size_t query = 0; query < count; ++query) {
            std::copy(queries[begin + query].begin(), queries[begin + query].end(), queryWords.begin() + query * words);
            queryCounts[query] = static_cast<std::int32_t>(
                smsd::popcountWords(queries[begin + query].data(), queries[begin + query].size()));
            for (int targetBits = 0; targetBits <= bits; ++targetBits) {
                int lower = 0, upper = std::min(queryCounts[query], targetBits) + 1;
                while (lower < upper) {
                    const int common = lower + (upper - lower) / 2;
                    if (score(common, queryCounts[query], targetBits) >= threshold) upper = common;
                    else lower = common + 1;
                }
                cutoffs[query * stride + targetBits] = lower;
            }
        }
        std::vector<smsd_fingerprint_hit> hits(count * targetCount);
        std::vector<std::uint64_t> offsets(count + 1, 0);
        if (table->fingerprint_screen(cache.resident, static_cast<std::int32_t>(count), queryWords.data(),
                                      cutoffs.data(), static_cast<std::int32_t>(stride), hits.data(), hits.size(),
                                      offsets.data(), error, sizeof error) != 0) {
            failed();
            return std::nullopt;
        }
        if (offsets[0] != 0 || offsets.back() > hits.size()) {
            detail::recordFailure("the accelerator returned invalid fingerprint hit offsets");
            return std::nullopt;
        }
        for (std::size_t query = 0; query < count; ++query) {
            if (offsets[query] > offsets[query + 1] || offsets[query + 1] > hits.size()) {
                detail::recordFailure("the accelerator returned invalid fingerprint hit offsets");
                return std::nullopt;
            }
            int previous = -1;
            auto& selected = results[begin + query];
            for (std::uint64_t position = offsets[query]; position < offsets[query + 1]; ++position) {
                const auto hit = hits[static_cast<std::size_t>(position)];
                if (hit.target_index <= previous || hit.target_index < 0
                    || static_cast<std::size_t>(hit.target_index) >= targetCount || hit.common_bits < 0
                    || hit.common_bits > std::min(queryCounts[query], cache.popcounts[hit.target_index])) {
                    detail::recordFailure("the accelerator returned an invalid fingerprint hit");
                    return std::nullopt;
                }
                previous = hit.target_index;
                if (!(score(hit.common_bits, queryCounts[query], cache.popcounts[hit.target_index]) >= threshold)) {
                    detail::recordFailure("the accelerator returned a fingerprint hit below the threshold");
                    return std::nullopt;
                }
                selected.push_back(hit.target_index);
            }
        }
    }
    return results;
}

/** Integer resident label overlaps, in query-major order; no score is formed on the device. */
template<class Targets, class Consume>
inline bool residentLabelOverlapsForTargets(
    const std::vector<const MolGraph*>& queries, const Targets& targets,
    const screening_detail::TargetLabels& labels, const ChemOptions& chemistry,
    ResidentCache& cache, const Consume& consume, const smsd_screen_accelerator_v1* provider = nullptr) {
    using screening_detail::QueryLabelCounts;
    using screening_detail::TargetLabels;
    static_assert(sizeof(TargetLabels::LabelCount) == sizeof(smsd_label_count)
                      && offsetof(TargetLabels::LabelCount, label) == offsetof(smsd_label_count, label)
                      && offsetof(TargetLabels::LabelCount, count) == offsetof(smsd_label_count, count),
                  "the label histogram layout must match the accelerator interface");
    const auto* table = provider != nullptr ? provider : active();
    const std::size_t targetCount = targets.size();
    if (!hasResidentLabels(table) || targetCount == 0 || queries.empty()
        || targetCount > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || labels.size() != targetCount || labels.mode() != screening_detail::labelMode(chemistry)
        || (provider == nullptr && !reachesMinimum(table, 2, targetCount))) return false;
    constexpr int tableSize = QueryLabelCounts::denseLimit;
    for (std::size_t query = 0; query < queries.size(); ++query)
        for (int atom = 0; atom < queries[query]->n; ++atom) {
            const int label = screening_detail::ubLabel(*queries[query], atom, chemistry);
            if (label < 0 || label >= tableSize) return false;
        }
    const std::lock_guard<std::mutex> lock(cache.mutex);
    if (cache.table != table) cache.clear();
    char error[256] = {};
    if (cache.resident == nullptr) {
        const TargetLabels::LabelCount* first = nullptr;
        std::vector<std::uint64_t> offsets(targetCount + 1, 0);
        for (std::size_t target = 0; target < targetCount; ++target) {
            offsets[target + 1] = offsets[target];
            if (labels.held(static_cast<int>(target)) && labels.atomCount(static_cast<int>(target)) > 0) {
                if (first == nullptr) first = labels.begin(static_cast<int>(target));
                offsets[target + 1] += static_cast<std::uint64_t>(
                    labels.end(static_cast<int>(target)) - labels.begin(static_cast<int>(target)));
            }
        }
        void* resident = nullptr;
        if (table->create_labels(static_cast<std::int64_t>(targetCount), offsets.data(),
                                  reinterpret_cast<const smsd_label_count*>(first), &resident, error, sizeof error) != 0) {
            error[sizeof error - 1] = '\0';
            detail::recordFailure(error);
            return false;
        }
        if (resident == nullptr) {
            detail::recordFailure("the accelerator returned no resident label table");
            return false;
        }
        cache.table = table;
        cache.resident = resident;
    }
    const std::size_t chunkSize = queriesPerChunk(targetCount, queries.size());
    for (std::size_t begin = 0; begin < queries.size(); begin += chunkSize) {
        const std::size_t count = std::min(chunkSize, queries.size() - begin);
        std::vector<std::int32_t> queryCounts(count * tableSize, 0);
        for (std::size_t query = 0; query < count; ++query)
            for (int atom = 0; atom < queries[begin + query]->n; ++atom)
                ++queryCounts[query * tableSize + screening_detail::ubLabel(*queries[begin + query], atom, chemistry)];
        std::vector<std::int32_t> overlaps(count * targetCount, 0);
        if (table->resident_label_overlaps(cache.resident, static_cast<std::int32_t>(count),
                                            queryCounts.data(), tableSize,
                                            overlaps.data(), error, sizeof error) != 0) {
            error[sizeof error - 1] = '\0';
            detail::recordFailure(error);
            return false;
        }
        for (std::size_t query = 0; query < count; ++query) {
            for (std::size_t target = 0; target < targetCount; ++target) {
                const int overlap = overlaps[query * targetCount + target];
                const int held = labels.held(static_cast<int>(target)) ? labels.atomCount(static_cast<int>(target)) : 0;
                if (overlap < 0 || overlap > held || overlap > std::max(0, queries[begin + query]->n)) {
                    detail::recordFailure("the accelerator returned a label overlap larger than its target");
                    return false;
                }
            }

        }
        consume(begin, count, overlaps);
    }
    return true;
}


/** Similarity-bound indices formed by the unchanged host expression. */
template<class Targets>
inline std::optional<std::vector<std::vector<int>>> residentSimilarityScreensForTargets(
    const std::vector<const MolGraph*>& queries, const Targets& targets,
    const screening_detail::TargetLabels& labels, double threshold, const ChemOptions& chemistry, int numThreads,
    ResidentCache& cache, const smsd_screen_accelerator_v1* provider = nullptr) {
    using screening_detail::QueryLabelCounts;
    const std::size_t targetCount = targets.size();
    std::vector<std::vector<int>> results(queries.size());
    auto consume = [&](std::size_t begin, std::size_t count, const std::vector<std::int32_t>& overlaps) {
        for (std::size_t query = 0; query < count; ++query) {
            const QueryLabelCounts cpuCounts(*queries[begin + query], chemistry);
            std::vector<std::uint8_t> mask(targetCount, 0);
            screening_detail::forEachTarget(static_cast<int>(targetCount), numThreads, [&](int target, std::vector<int>& scratch) {
                const double score = labels.held(target)
                    ? screening_detail::overlapSimilarity(overlaps[query * targetCount + target],
                                                           queries[begin + query]->n, labels.atomCount(target))
                    : cpuCounts.similarity(targets[target], scratch);
                mask[target] = score >= threshold;
            });
            for (std::size_t target = 0; target < targetCount; ++target)
                if (mask[target]) results[begin + query].push_back(static_cast<int>(target));
        }
    };
    if (!residentLabelOverlapsForTargets(queries, targets, labels, chemistry, cache, consume, provider))
        return std::nullopt;
    return results;
}

/**
 * The scores screening_detail::similarityScores(query, targets, labels, chemistry, numThreads)
 * gives, with the label overlaps of the targets whose histograms are held computed by the active
 * accelerator. Empty when no accelerator is active, when the screen cannot read the histograms
 * (built under another label rule or for another number of targets, or a query label outside the
 * dense table) or when the accelerator fails; the caller then screens on the CPU. Targets left to
 * their graphs are scored on the CPU, as there.
 */
template<class Targets>
inline std::optional<std::vector<double>> similarityScoresForTargets(const MolGraph& query,
                                                           const Targets& targets,
                                                           const screening_detail::TargetLabels& labels,
                                                           const ChemOptions& chemistry, int numThreads) {
    using screening_detail::QueryLabelCounts;
    using screening_detail::TargetLabels;
    static_assert(sizeof(TargetLabels::LabelCount) == sizeof(smsd_label_count)
                      && offsetof(TargetLabels::LabelCount, label) == offsetof(smsd_label_count, label)
                      && offsetof(TargetLabels::LabelCount, count) == offsetof(smsd_label_count, count),
                  "the histogram entries must have the layout the accelerator interface gives them");

    const smsd_screen_accelerator_v1* table = active();
    const std::size_t N = targets.size();
    if (table == nullptr || N == 0 || N > static_cast<std::size_t>(std::numeric_limits<int>::max())
        || labels.size() != N || labels.mode() != screening_detail::labelMode(chemistry))
        return std::nullopt;

    constexpr int tableSize = QueryLabelCounts::denseLimit;
    std::vector<std::int32_t> queryCounts(tableSize, 0);
    for (int atom = 0; atom < query.n; ++atom) {
        const int label = screening_detail::ubLabel(query, atom, chemistry);
        if (label < 0 || label >= tableSize) return std::nullopt;
        ++queryCounts[label];
    }

    // The histogram of target i is entries[offsets[i], offsets[i + 1]); a target left to its graph has none.
    const TargetLabels::LabelCount* const first = labels.begin(0);
    std::vector<std::uint64_t> offsets(N + 1);
    for (std::size_t i = 0; i < N; ++i)
        offsets[i] = static_cast<std::uint64_t>(labels.begin(static_cast<int>(i)) - first);
    offsets[N] = static_cast<std::uint64_t>(labels.end(static_cast<int>(N - 1)) - first);

    std::vector<std::int32_t> overlaps(N, 0);
    char error[256] = {};
    if (table->label_overlaps(queryCounts.data(), tableSize, static_cast<std::int64_t>(N), offsets.data(),
                              reinterpret_cast<const smsd_label_count*>(first), overlaps.data(),
                              error, sizeof error) != 0) {
        error[sizeof error - 1] = '\0';
        detail::recordFailure(error);
        return std::nullopt;
    }
    for (std::size_t i = 0; i < N; ++i) {
        const int held = labels.held(static_cast<int>(i)) ? labels.atomCount(static_cast<int>(i)) : 0;
        if (overlaps[i] < 0 || overlaps[i] > held) {
            detail::recordFailure("the accelerator returned a label overlap larger than its target");
            return std::nullopt;
        }
    }

    std::vector<double> scores(N, 0.0);
    const QueryLabelCounts cpuCounts(query, chemistry);
    screening_detail::forEachTarget(static_cast<int>(N), numThreads, [&](int i, std::vector<int>& scratch) {
        scores[i] = labels.held(i)
            ? screening_detail::overlapSimilarity(overlaps[i], query.n, labels.atomCount(i))
            : cpuCounts.similarity(targets[i], scratch);
    });
    return scores;
}

inline std::optional<std::vector<std::vector<int>>> residentSimilarityScreens(
    const std::vector<const MolGraph*>& queries, const std::vector<MolGraph>& targets,
    const screening_detail::TargetLabels& labels, double threshold, const ChemOptions& chemistry, int numThreads,
    ResidentCache& cache, const smsd_screen_accelerator_v1* provider = nullptr) {
    return ::smsd::screen_accelerator::residentSimilarityScreensForTargets<std::vector<MolGraph>>(
        queries, targets, labels, threshold, chemistry, numThreads, cache, provider);
}

inline std::optional<std::vector<double>> similarityScores(const MolGraph& query,
    const std::vector<MolGraph>& targets, const screening_detail::TargetLabels& labels,
    const ChemOptions& chemistry, int numThreads) {
    return ::smsd::screen_accelerator::similarityScoresForTargets<std::vector<MolGraph>>(
        query, targets, labels, chemistry, numThreads);
}

} // namespace screen_accelerator
} // namespace smsd

#endif // SMSD_SCREEN_ACCELERATOR_HPP

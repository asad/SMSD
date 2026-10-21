/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Batch dispatch: the work-sharing dispatcher, the thread team and matchMany(),
 * which run one query over many targets.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/search_types.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <cerrno>
#include <cstdio>
#include <memory>
#include <string>
#elif defined(_WIN32)
// Declared as the Windows headers declare it, so that this header does not include them.
extern "C" __declspec(dllimport) unsigned long __stdcall GetActiveProcessorCount(unsigned short groupNumber);
#endif

namespace smsd {
namespace search_detail {

/**
 * Hands out batch indices in contiguous chunks whose length follows the measured cost of each
 * worker's previous chunk. Expensive targets are claimed one at a time so the load stays
 * balanced; cheap targets are grouped so that workers neither contend on the shared counter
 * nor write neighbouring result slots. Chunks shrink again as the batch nears its end.
 *
 * An optional order lists the indices in the sequence they are handed out. Every index below the
 * lowest failing index is still processed, so the reported failure does not depend on the order,
 * the number of workers or their timing.
 */
class BatchDispatcher {
    using Clock = std::chrono::steady_clock;
    /** Roughly this much work per claim keeps the shared counter cold. */
    static constexpr int64_t chunkNanoseconds = 20000;
    static constexpr size_t maximumChunk = 256;
    static constexpr size_t noFailure = std::numeric_limits<size_t>::max();

    const size_t size_;
    const size_t workers_;
    const std::vector<size_t>* order_;
    // Lines of 128 bytes keep the shared counters apart on every supported processor.
    alignas(128) std::atomic<size_t> next_{0};
    alignas(128) std::atomic<size_t> firstFailure_{noFailure};
    std::mutex failureMutex_;
    std::exception_ptr failure_;

    void recordFailure(size_t index) noexcept {
        const std::lock_guard<std::mutex> lock(failureMutex_);
        if (index < firstFailure_.load(std::memory_order_relaxed)) {
            failure_ = std::current_exception();
            firstFailure_.store(index, std::memory_order_relaxed);
        }
    }
    size_t nextChunk(size_t chunk, size_t done, int64_t elapsedNanoseconds) const {
        const size_t claimed = std::min(size_, next_.load(std::memory_order_relaxed));
        const size_t tail = std::max<size_t>(1, (size_ - claimed) / (2 * workers_));
        const int64_t perIndex = std::max<int64_t>(1, elapsedNanoseconds / static_cast<int64_t>(done));
        const auto wanted = static_cast<size_t>(std::max<int64_t>(1, chunkNanoseconds / perIndex));
        return std::max<size_t>(1, std::min({wanted, chunk * 2, maximumChunk, tail}));
    }
public:
    /** order, when given, must hold each index below size once and outlive the dispatcher. */
    BatchDispatcher(size_t size, size_t workers, const std::vector<size_t>* order = nullptr)
        : size_(size), workers_(std::max<size_t>(1, workers)), order_(order && !order->empty() ? order : nullptr) {}
    BatchDispatcher(const BatchDispatcher&) = delete;
    BatchDispatcher& operator=(const BatchDispatcher&) = delete;

    /** Run task(index) for claimed indices until none remain; safe to call from every worker. */
    template<class Task>
    void work(const Task& task) noexcept {
        size_t chunk = 1;
        for (;;) {
            const size_t begin = next_.fetch_add(chunk, std::memory_order_relaxed);
            if (begin >= size_) return;
            const size_t end = std::min(size_, begin + chunk);
            const auto started = Clock::now();
            for (size_t position = begin; position < end; ++position) {
                const size_t index = order_ ? (*order_)[position] : position;
                if (index > firstFailure_.load(std::memory_order_relaxed)) {
                    // In input order every later claim is higher still, so nothing is left to report.
                    if (!order_) return;
                    continue;
                }
                try { task(index); }
                catch (...) { recordFailure(index); }
            }
            chunk = nextChunk(chunk, end - begin,
                std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - started).count());
        }
    }
    /** Rethrow the failure of the lowest failing index, if any. */
    void rethrowFailure() const {
        if (failure_) std::rethrow_exception(failure_);
    }
};

/**
 * The processors this process may use as the platform reports them, or 0 when it reports none:
 * on Linux those of its affinity mask, on Windows the active ones of every processor group.
 */
inline size_t availableProcessors() {
#if defined(__linux__) && defined(CPU_ALLOC) && defined(CPU_COUNT_S)
    struct Release { void operator()(cpu_set_t* set) const { CPU_FREE(set); } };
    // The set grows until it covers every processor the kernel knows.
    for (size_t processors = 1024; processors <= (size_t(1) << 20); processors *= 2) {
        const std::unique_ptr<cpu_set_t, Release> set(CPU_ALLOC(processors));
        if (!set) return 0;
        const size_t bytes = CPU_ALLOC_SIZE(processors);
        if (sched_getaffinity(0, bytes, set.get()) == 0) return static_cast<size_t>(CPU_COUNT_S(bytes, set.get()));
        if (errno != EINVAL) return 0;
    }
    return 0;
#elif defined(_WIN32)
    constexpr unsigned short allProcessorGroups = 0xffff;
    return static_cast<size_t>(GetActiveProcessorCount(allProcessorGroups));
#else
    return 0;
#endif
}

/**
 * The processors the cgroup v2 CPU quotas of this process allow, rounded up and the least along
 * its cgroup and the cgroups above it, or 0 when none applies. Quotas are read once per process.
 */
inline size_t quotaProcessors() {
#if defined(__linux__)
    static const size_t processors = [] {
        struct Close { void operator()(std::FILE* file) const { std::fclose(file); } };
        auto read = [](const std::string& path) {
            std::string text;
            const std::unique_ptr<std::FILE, Close> file(std::fopen(path.c_str(), "r"));
            if (!file) return text;
            char buffer[512];
            size_t length;
            while ((length = std::fread(buffer, 1, sizeof(buffer), file.get())) > 0 && text.size() < 65536)
                text.append(buffer, length);
            return text;
        };
        // The unified hierarchy is listed as "0::<path>".
        std::string path;
        const std::string groups = read("/proc/self/cgroup");
        for (size_t start = 0; start < groups.size();) {
            size_t end = groups.find('\n', start);
            if (end == std::string::npos) end = groups.size();
            if (groups.compare(start, 3, "0::") == 0) path = groups.substr(start + 3, end - start - 3);
            start = end + 1;
        }
        // A cgroup outside this namespace shows as a path up from its root; only the root is read then.
        if (path.find("..") != std::string::npos || path == "/") path.clear();
        const std::string root = "/sys/fs/cgroup";
        std::string directory = root + path;
        size_t least = 0;
        for (;;) {
            unsigned long long quota = 0, period = 0;
            const std::string limit = read(directory + "/cpu.max");
            if (std::sscanf(limit.c_str(), "%llu %llu", &quota, &period) == 2 && quota != 0 && period != 0) {
                const auto allowed = static_cast<size_t>((quota + period - 1) / period);
                if (least == 0 || allowed < least) least = allowed;
            }
            if (directory.size() <= root.size()) break;
            directory.erase(directory.rfind('/'));
        }
        return least;
    }();
    return processors;
#else
    return 0;
#endif
}

/**
 * The most workers a batch starts: one for each processor this process may use, or 64 when that
 * count is unknown, and on Linux at most what the cgroup CPU quota allows. Results do not depend on it.
 */
inline size_t maximumWorkers() {
    size_t processors = availableProcessors();
    if (processors == 0) processors = std::thread::hardware_concurrency();
    if (processors == 0) processors = 64;
    const size_t quota = quotaProcessors();
    return quota != 0 && quota < processors ? quota : processors;
}

/**
 * Run body on the calling thread and on up to workers - 1 standard threads; body must not throw.
 * No thread outlives the call, so a process that forks afterwards, as Python multiprocessing
 * does, inherits no thread team. Helpers that fail to start are left out, as the dispatcher
 * balances any team size, and every helper that started is joined before the call returns.
 */
template<class Body>
inline void runTeam(size_t workers, const Body& body) {
    if (workers <= 1) { body(); return; }
    std::vector<std::thread> helpers;
    struct Joiner {
        std::vector<std::thread>& threads;
        ~Joiner() {
            for (auto& thread : threads)
                if (thread.joinable()) thread.join();
        }
    } joiner{helpers};
    try {
        helpers.reserve(workers - 1);
        for (size_t helper = 1; helper < workers; ++helper) helpers.emplace_back([&body] { body(); });
    } catch (...) {
        // Continue with the helpers that started.
    }
    body();
}

template<class Target, class = void>
struct HasAtomCount : std::false_type {};
template<class Target>
struct HasAtomCount<Target, std::void_t<decltype(std::declval<const Target&>().n)>> : std::true_type {};

/** Indices by descending atom count, ties in input order (a stable counting sort when it is cheap). */
template<class Targets>
inline std::vector<size_t> largestTargetsFirst(const Targets& targets) {
    const size_t size = targets.size();
    std::vector<size_t> order(size);
    int largest = 0;
    for (size_t index = 0; index < size; ++index) largest = std::max(largest, static_cast<int>(targets[index].n));
    auto key = [&](size_t index) { return static_cast<size_t>(largest - std::max(0, static_cast<int>(targets[index].n))); };
    if (static_cast<size_t>(largest) <= 4 * size + 64) {
        std::vector<size_t> start(static_cast<size_t>(largest) + 2, 0);
        for (size_t index = 0; index < size; ++index) ++start[key(index) + 1];
        for (size_t bucket = 1; bucket < start.size(); ++bucket) start[bucket] += start[bucket - 1];
        for (size_t index = 0; index < size; ++index) order[start[key(index)]++] = index;
    } else {
        std::iota(order.begin(), order.end(), size_t(0));
        std::stable_sort(order.begin(), order.end(), [&](size_t first, size_t second) { return key(first) < key(second); });
    }
    return order;
}

template<class Query, class Targets>
inline std::vector<SearchResult> matchMany(const Query& query, const Targets& targets, int numThreads) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
    std::vector<SearchResult> results(targets.size());
    if (results.empty()) return results;
    const size_t requested = numThreads == 0 ? std::max(1u, std::thread::hardware_concurrency()) : static_cast<size_t>(numThreads);
    const size_t wanted = std::min(results.size(), requested);
    const size_t workers = wanted > 1 ? std::min(wanted, maximumWorkers()) : 1;
    // In a short batch one late, large target can leave the other workers idle, so larger
    // targets, which tend to cost more, start first. Long batches keep input order, which keeps
    // neighbouring results on one worker and costs nothing to set up.
    std::vector<size_t> order;
    using Target = std::decay_t<decltype(targets[0])>;
    if constexpr (HasAtomCount<Target>::value)
        if (workers > 1 && results.size() <= 64 * workers) order = largestTargetsFirst(targets);
    // A short batch of cheap targets finishes on the calling thread sooner than a team can
    // start, so targets run here in input order until about 100 microseconds of work is done.
    // Every earlier target has succeeded, so a failure here is the lowest failing index.
    size_t done = 0;
    if (workers > 1) {
        const auto started = std::chrono::steady_clock::now();
        while (done < results.size()) {
            results[done] = query.match(targets[done]);
            ++done;
            if (std::chrono::steady_clock::now() - started > std::chrono::microseconds(100)) break;
        }
        if (done == results.size()) return results;
        if (order.empty()) {
            order.resize(results.size() - done);
            std::iota(order.begin(), order.end(), done);
        } else {
            order.erase(std::remove_if(order.begin(), order.end(),
                [done](size_t index) { return index < done; }), order.end());
        }
    }
    BatchDispatcher dispatcher(order.empty() ? results.size() : order.size(), workers, &order);
    runTeam(workers, [&] {
        dispatcher.work([&](size_t index) { results[index] = query.match(targets[index]); });
    });
    dispatcher.rethrowFailure();
    return results;
}

} // namespace search_detail
} // namespace smsd

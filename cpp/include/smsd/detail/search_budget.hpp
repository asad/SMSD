/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Search time limits and the Budget that every stage of one search call shares for
 * time, cancellation and node counts.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/search_types.hpp"
#include "smsd/time_budget.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>
#include <memory>

namespace smsd {
namespace search_detail {

/**
 * The time limit of a search in milliseconds: the options' own, or, when they set none, the
 * default for a common-subgraph search of these sizes or for a substructure search.
 */
inline int64_t searchTimeLimitMillis(const SearchOptions& options, int queryAtoms, int targetAtoms, bool commonSubgraph) {
    return options.timeoutMs < 0
        ? (commonSubgraph ? std::min<int64_t>(30000, 500 + int64_t(queryAtoms) * targetAtoms * 2) : 10000)
        : options.timeoutMs;
}

/** The time from the start of a search until which the results it found in time are built: the limit and an eighth of it. */
inline int64_t resultsTimeMillis(int64_t limitMillis) {
    const int64_t grace = limitMillis > 0 ? limitMillis / 8 : 0;
    return limitMillis > std::numeric_limits<int64_t>::max() - grace ? std::numeric_limits<int64_t>::max()
                                                                     : limitMillis + grace;
}

/** Every stage and recursive query in a call shares this budget. */
class Budget {
    using Clock = std::chrono::steady_clock;
    Clock::time_point started_ = Clock::now();
    Clock::time_point deadline_;
    /** Results found in time are built until an eighth of the time limit past the deadline. */
    Clock::time_point resultsDeadline_;
    uint64_t nodeLimit_;
    // The budget shares the token, so it stays valid however long the call runs, even if the
    // caller's options drop it meanwhile, for example from another thread during consensus().
    std::shared_ptr<const CancellationToken> cancellation_;
    SearchStatus status_ = SearchStatus::Completed;
    uint64_t nodes_ = 0;
    uint32_t polls_ = 0;
public:
    /** The budget copies what it needs, so options may change or end after construction. */
    Budget(const SearchOptions& options, int queryAtoms, int targetAtoms, bool commonSubgraph)
        : nodeLimit_(options.nodeLimit), cancellation_(options.cancellation) {
        const int64_t milliseconds = searchTimeLimitMillis(options, queryAtoms, targetAtoms, commonSubgraph);
        deadline_ = detail::steadyDeadline(milliseconds, started_);
        resultsDeadline_ = detail::steadyDeadline(resultsTimeMillis(milliseconds), started_);
    }
    /** A caller that skips the clock must still read it at a bounded interval. */
    bool stopped(bool readClock = true) {
        if (status_ != SearchStatus::Completed) return true;
        if (cancellation_ && cancellation_->isCancelled()) status_ = SearchStatus::Cancelled;
        else if (readClock && Clock::now() >= deadline_) status_ = SearchStatus::TimedOut;
        return status_ != SearchStatus::Completed;
    }
    bool visit(bool readClock = true) {
        if (stopped(readClock)) return false;
        if (nodeLimit_ != 0 && nodes_ >= nodeLimit_) {
            status_ = SearchStatus::LimitReached;
            return false;
        }
        ++nodes_;
        return true;
    }
    /**
     * Whether building the results of a finished search should stop: on cancellation, or past
     * the results deadline, when a completed search is reported as timed out.
     */
    bool resultsStopped() {
        SearchStatus reason = SearchStatus::Completed;
        if (cancellation_ && cancellation_->isCancelled()) reason = SearchStatus::Cancelled;
        else if (Clock::now() >= resultsDeadline_) reason = SearchStatus::TimedOut;
        if (reason == SearchStatus::Completed) return false;
        if (status_ == SearchStatus::Completed) status_ = reason;
        return true;
    }
    /** For short node steps: cancellation is checked on each call, the clock on every sixteenth. */
    bool stoppedSampled() { return stopped((++polls_ & 15u) == 0); }
    bool visitSampled() { return visit((++polls_ & 15u) == 0); }
    SearchStatus status() const { return status_; }
    uint64_t nodesVisited() const { return nodes_; }
    /** Counts nodes that a search visited under its own work clock. */
    void countNodes(uint64_t count) { nodes_ += count; }
    /** Marks the node limit, or another limit of the search, as reached unless it already stopped. */
    void limitReached() {
        if (status_ == SearchStatus::Completed) status_ = SearchStatus::LimitReached;
    }
    int64_t elapsedMillis() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - started_).count();
    }
};

} // namespace search_detail
} // namespace smsd

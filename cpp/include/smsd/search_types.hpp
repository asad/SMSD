/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Search options, matches and results shared by Substructure and MCS.
 * Match::toGraph() and Match::toSmiles() are defined in detail/match_graph.hpp,
 * which this header includes.
 */
#pragma once
#ifndef SMSD_SEARCH_TYPES_HPP
#define SMSD_SEARCH_TYPES_HPP

#include "smsd/mol_graph.hpp"

#include <array>
#include <atomic>
#include <climits>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace smsd {

enum class SearchStatus { Completed, TimedOut, Cancelled, LimitReached };
enum class MCSObjective { MCIS, MCES };
enum class MCSConnectivity { Connected, Fragments };

/** A cancellation flag that may be shared with an active search. */
class CancellationToken {
    std::atomic<bool> cancelled_{false};
public:
    void cancel() noexcept { cancelled_.store(true, std::memory_order_release); }
    bool isCancelled() const noexcept { return cancelled_.load(std::memory_order_acquire); }
};

struct SearchOptions {
    ChemOptions chemistry;
    /** Negative selects the adaptive budget; zero requests an immediate stop. */
    int64_t timeoutMs = -1;
    /** Zero leaves node count unrestricted; the time budget still applies. */
    uint64_t nodeLimit = 0;
    std::shared_ptr<CancellationToken> cancellation;
};

struct SubstructureOptions : SearchOptions {
    bool induced = false;
};

struct MCSOptions : SearchOptions {
    MCSObjective objective = MCSObjective::MCIS;
    MCSConnectivity connectivity = MCSConnectivity::Connected;
    std::vector<double> atomWeights;
    /** Original target atom indices unavailable to this search. */
    std::vector<int> excludedTargetAtoms;
    int minFragmentAtoms = 1;
    int maxFragments = INT_MAX;
};

struct Match {
    std::map<int, int> atomMapping;
    /** Each entry is {query atom, query atom, target atom, target atom}. */
    std::vector<std::array<int, 4>> bondMapping;
    /** Connected components of selected common bonds, as query atom indices. */
    std::vector<std::vector<int>> fragments;
    /** Selected graph with original perception and hydrogen counts. */
    MolGraph toGraph(const MolGraph& query) const;
    std::string toSmiles(const MolGraph& query) const;
};

struct SearchResult {
    std::vector<Match> matches;
    SearchStatus status = SearchStatus::Completed;
    bool optimal = false;
    bool complete = false;
    /** Primary objective value; weighted scores use the supplied weight units. */
    double score = 0.0;
    uint64_t nodesVisited = 0;
    int64_t elapsedMillis = 0;
};

} // namespace smsd

#include "smsd/detail/match_graph.hpp"

#endif // SMSD_SEARCH_TYPES_HPP

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The common-bond (MCES) search. Its parts are in detail/: common_bond_classes.hpp,
 * common_bond_graphs.hpp, common_bond_conditions.hpp, common_bond_mapping_search.hpp
 * and common_bond_embedding.hpp.
 */
#pragma once

#include "smsd/bitops.hpp"
#include "smsd/detail/union_find.hpp"
#include "smsd/mol_graph.hpp"
#include "smsd/stereo.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "smsd/detail/common_bond_classes.hpp"
#include "smsd/detail/common_bond_graphs.hpp"
#include "smsd/detail/common_bond_conditions.hpp"
#include "smsd/detail/common_bond_mapping_search.hpp"
#include "smsd/detail/common_bond_embedding.hpp"

/*
 * The maximum common edge subgraph (MCES) of two molecules: the mapping from the query to the
 * target with the most common bonds and then the most atoms, or every distinct such mapping up
 * to a limit. Each mapped atom has a common bond, and every compatible bond between mapped atoms
 * is common. The common bonds form one fragment or, unless the search is connected, several
 * fragments anywhere in either molecule. Where the options also ask for a required mapping,
 * fragment sizes or counts, stereochemistry or whole rings, each mapping recorded is checked by
 * CommonBondConditions, under which it may select only some of its common bonds; the common
 * bonds still bound every selection, so the bounds hold. Every pair mapped counts towards the
 * node limit.
 *
 * CommonBondGraphs holds both molecules as int graphs with the atom and bond pairs that may be
 * mapped onto each other. The atom and bond predicates are supplied by the caller, so the
 * chemistry is exactly that of the rest of the library. Atoms, or bonds, of one molecule share
 * a kind when they match the same items of the other, so whether two items match depends on
 * their kinds alone. A class is a connected part of the relation between the kinds, so a
 * compatible pair always shares a class. Where the relation is not an equivalence, as with
 * flexible aromaticity, isotope labels left unset or ring counts compared only between ring
 * atoms, a class may hold kinds that do not all match; the items a mapping can pair within it
 * are then bounded by a maximum b-matching over the kinds that do. Every bound is a true upper
 * bound either way.
 *
 * CommonBondSearch finds the mapping. When the whole query embeds in the target, an embedding
 * is the answer. Otherwise the query atoms are ranked, fewest candidates first, greedy seeds
 * give a score to beat, and a depth-first branch and bound tries each query atom as a root,
 * highest rank first, and grows the mapping pair by pair on an explicit stack of levels. After
 * the first pair of a fragment, each pair maps a query atom next to the fragment onto a free
 * neighbour of the image of one of its mapped neighbours. A level that cannot beat the best
 * score is taken off, and the pair it tried is excluded below its later siblings, so each
 * mapping is found once. Unless the search is connected, the best score is first sought level
 * by level of common bonds, down from an upper bound, within a small budget; one search from the
 * best mapping known then finishes the work, now and then scouting the root pairs ahead for a
 * better mapping. Every loop ticks a Clock, which throws Stop when the time limit is reached or
 * the search is cancelled; the best mapping found by then is kept.
 *
 * A fragment search for one mapping without conditions runs from the molecule with fewer bond
 * incidences, and on some hard pairs the search the other way round is far shorter. Once its
 * last search has mapped ScoutLimits::race pairs, the search the other way round is raced
 * against it from the best mapping and the bound so far, which hold either way round: it maps
 * ScoutLimits::raceTrial pairs, and up to ScoutLimits::raceLonger while it holds a better
 * mapping or a lower bound than the last search had. A race that proves its mapping ends the
 * search; otherwise its best mapping and bound go to the last search, which goes on from where
 * it was. Searches that finish sooner are as before; the nodes are those of both.
 */
namespace smsd {
namespace common_bond_detail {

/** What a common-bond search found. */
struct CommonBondOutcome {
    /** Whether the search could take the graphs and the bounds it needs; when not, nothing else is set. */
    bool supported = false;
    /** Whether the search ran to its end, so that the mappings are maximum ones. */
    bool finished = false;
    /** Whether the search stopped at its node limit. */
    bool capReached = false;
    /** The query-to-target atom mapping, -1 where unmapped, or empty when nothing was found. */
    std::vector<int> mapping;
    /** The selected common bonds of the mapping, see CommonBondGraphs::commonBonds(). */
    std::vector<std::array<int, 4>> bonds;
    /** Every mapping found, the first being `mapping`, and the selected common bonds of each. */
    std::vector<std::vector<int>> mappings;
    std::vector<std::vector<std::array<int, 4>>> bondLists;
    /** The pairs mapped and the embedding atoms placed. */
    uint64_t nodes = 0;
    /** The scoutings, and the root pairs skipped as refuted, which tests read. */
    int64_t scoutings = 0, refutedRoots = 0;
};

/** What a common-bond search is asked for beyond the single best mapping of its graphs. */
struct CommonBondRequest {
    /** The most mappings held; a search that holds this many may have more of their score. */
    size_t wanted = 1;
    /** The most pairs mapped and embedding atoms placed before the search stops. */
    int64_t nodeCap = INT64_MAX;
    /** The fewest atoms of a fragment, and the most fragments. */
    int minFragmentAtoms = 1, maxFragments = INT_MAX;
    /**
     * The target that a required mapping fixes for each query atom, -1 for none, and the query
     * bonds it selects, each as its two atoms; both empty without a required mapping.
     */
    std::vector<int> requiredTarget;
    std::vector<std::array<int, 2>> requiredBonds;
    /** Whether mapped chirality, double-bond stereochemistry and whole rings are checked. */
    bool chirality = false, bondStereo = false, rings = false;
    /** Per-call relevant cycles already held by the caller, or null when none were prepared. */
    const std::vector<std::vector<int>>* queryCycles = nullptr;
    const std::vector<std::vector<int>>* targetCycles = nullptr;
};

/**
 * Searches `query` in `target` for the best MCES under the given atom and bond predicates, with
 * the target atoms in `excludedTargetAtoms` left out, as `request` asks: at most
 * `request.wanted` distinct best mappings that meet its conditions. `stopped()` is asked now and
 * then whether the time is up or the search was cancelled; the best mappings found by then are
 * the result, as they are once the search reaches its node limit. With `connected` the common
 * bonds form one fragment. The outcome is not supported when the graphs are not well formed or
 * too large, see CommonBondGraphs::supports() and CommonBondSearch::fits(). Only tests change
 * `limits`.
 */
template <class AtomMatch, class BondMatch, class BondKey, class Stopped>
inline CommonBondOutcome searchCommonBonds(const MolGraph& query, const MolGraph& target, AtomMatch atomMatches,
                                           BondMatch bondMatches, BondKey bondKey, const std::vector<int>& excludedTargetAtoms,
                                           bool connected, Stopped stopped, const CommonBondRequest& request,
                                           const ScoutLimits& limits = ScoutLimits()) {
    CommonBondOutcome outcome;
    if (!CommonBondGraphs::supports(query, target)) return outcome;
    outcome.supported = true;
    Clock clock{std::function<bool()>(std::move(stopped))};
    std::unique_ptr<CommonBondGraphs> graph;
    std::unique_ptr<CommonBondConditions> conditions;
    std::unique_ptr<CommonBondSearch> search;
    // The graphs and the search the other way round, built when first needed.
    std::unique_ptr<CommonBondGraphs> reversedGraph;
    std::unique_ptr<CommonBondSearch> reversedSearch;
    std::vector<std::vector<int>> found;
    // A mapping of `from` atoms, turned round into one of `to` atoms.
    auto turn = [](const std::vector<int>& mapping, int from, int to) {
        std::vector<int> turned(size_t(to), -1);
        for (int atom = 0; atom < from && atom < static_cast<int>(mapping.size()); ++atom)
            if (mapping[atom] >= 0) turned[size_t(mapping[atom])] = atom;
        return turned;
    };
    try {
        clock.poll();
        graph = std::make_unique<CommonBondGraphs>(query, target, atomMatches, bondMatches, bondKey, excludedTargetAtoms, clock);
        if (!CommonBondSearch::fits(*graph)) return CommonBondOutcome{};
        // Chirality and double-bond stereochemistry count only when both molecules have them, and
        // rings only when either has one.
        const bool chirality = request.chirality && CommonBondConditions::hasCentres(query)
            && CommonBondConditions::hasCentres(target);
        const bool bondStereo = request.bondStereo && CommonBondConditions::hasStereoBonds(query)
            && CommonBondConditions::hasStereoBonds(target);
        const bool rings = request.rings && (CommonBondConditions::hasRings(query, request.queryCycles)
            || CommonBondConditions::hasRings(target, request.targetCycles));
        const bool sizes = request.minFragmentAtoms > 2 || (!connected && request.maxFragments < INT_MAX);
        if (chirality || bondStereo || rings || sizes || !request.requiredTarget.empty())
            conditions = std::make_unique<CommonBondConditions>(query, target, *graph, connected, request.minFragmentAtoms,
                                                                request.maxFragments, request.requiredTarget,
                                                                request.requiredBonds, chirality, bondStereo, rings, clock,
                                                                request.queryCycles, request.targetCycles);
        search = std::make_unique<CommonBondSearch>(*graph, clock, connected, limits, conditions.get(), request.requiredTarget,
                                                    request.wanted, request.nodeCap);
        // Without a required mapping or excluded atoms, the search the other way round is of the
        // same mappings turned round, so its mappings are mappings of this search and a score
        // that no mapping beats holds for both.
        const bool symmetric = !connected && request.requiredTarget.empty() && excludedTargetAtoms.empty()
            && CommonBondGraphs::supports(target, query);
        auto reversed = [&]() -> CommonBondSearch* {
            if (!reversedGraph) {
                auto atomsReversed = [&](int targetAtom, int queryAtom) { return atomMatches(queryAtom, targetAtom); };
                auto bondsReversed = [&](int targetFirst, int targetSecond, int queryFirst, int querySecond) {
                    return bondMatches(queryFirst, querySecond, targetFirst, targetSecond);
                };
                reversedGraph = std::make_unique<CommonBondGraphs>(target, query, atomsReversed, bondsReversed, bondKey,
                                                                   excludedTargetAtoms, clock);
                if (CommonBondSearch::fits(*reversedGraph))
                    reversedSearch = std::make_unique<CommonBondSearch>(*reversedGraph, clock, connected, limits);
            }
            return reversedSearch.get();
        };
        // A fragment search grows the seed the other way round only once its levels run out of
        // budget, see CommonBondSearch::seedAlsoFrom().
        if (symmetric)
            search->seedAlsoFrom([&] {
                CommonBondSearch* other = reversed();
                const std::vector<int> seed = other ? other->greedySeed() : std::vector<int>();
                return seed.empty() ? seed : turn(seed, target.n, query.n);
            });
        // The race, see the description at the top of this file.
        if (symmetric && request.wanted == 1 && !conditions)
            search->raceAfter(limits.race, [&] {
                CommonBondSearch* other = reversed();
                if (!other) return;
                const std::vector<std::vector<int>> held = search->found();
                const int64_t best = search->bestScore(), bound = search->upperBound();
                const int64_t start = other->nodes();
                const int64_t left = request.nodeCap == INT64_MAX ? INT64_MAX : request.nodeCap - search->nodes() - start;
                const int64_t end = std::min(left, limits.raceLonger);
                other->capNodesAt(start + std::min(left, limits.raceTrial), [&, start, end, best, bound] {
                    return other->bestFound() > best || other->upperBound() < bound ? start + end : other->nodes();
                });
                std::vector<std::vector<int>> answer;
                bool proved = true;
                try {
                    answer = other->raceFrom(held.empty() ? std::vector<int>() : turn(held.front(), query.n, target.n), bound);
                } catch (const Stop&) {
                    proved = false;
                    answer = other->found();
                    // Stopped by the clock rather than its node cap, the race ends the search.
                    if (!other->capReached()) {
                        if (!answer.empty()) search->offer(turn(answer.front(), target.n, query.n), INT64_MAX);
                        throw;
                    }
                }
                const std::vector<int> mapping = answer.empty() ? std::vector<int>() : turn(answer.front(), target.n, query.n);
                search->offer(mapping, proved ? (answer.empty() ? bound : other->scoreOfMapping(answer.front())) : other->upperBound());
                other->capNodesAt(INT64_MAX);
                if (request.nodeCap != INT64_MAX) search->capNodesAt(request.nodeCap - other->nodes());
            });
        found = search->run();
        outcome.finished = true;
    } catch (const Stop&) {
        // The time ran out, the search was cancelled or it reached its node limit. What the
        // search found by then is the result.
        if (search) found = search->found();
    }
    for (auto& record : found) {
        outcome.bondLists.push_back(graph->commonBonds(record));
        record.resize(size_t(query.n));
        outcome.mappings.push_back(std::move(record));
    }
    if (!outcome.mappings.empty()) {
        outcome.mapping = outcome.mappings.front();
        outcome.bonds = outcome.bondLists.front();
    }
    if (search) {
        outcome.nodes = static_cast<uint64_t>(search->nodes());
        outcome.scoutings = search->scoutings();
        outcome.refutedRoots = search->refutedRoots();
        outcome.capReached = search->capReached();
    }
    if (reversedSearch) {
        outcome.nodes += static_cast<uint64_t>(reversedSearch->nodes());
        outcome.scoutings += reversedSearch->scoutings();
        outcome.refutedRoots += reversedSearch->refutedRoots();
    }
    return outcome;
}

/** The single best MCES, as searchCommonBonds() above with nothing else asked. */
template <class AtomMatch, class BondMatch, class BondKey, class Stopped>
inline CommonBondOutcome searchCommonBonds(const MolGraph& query, const MolGraph& target, AtomMatch atomMatches,
                                           BondMatch bondMatches, BondKey bondKey, const std::vector<int>& excludedTargetAtoms,
                                           bool connected, Stopped stopped, const ScoutLimits& limits = ScoutLimits()) {
    return searchCommonBonds(query, target, atomMatches, bondMatches, bondKey, excludedTargetAtoms, connected,
                             std::move(stopped), CommonBondRequest(), limits);
}

} // namespace common_bond_detail
} // namespace smsd

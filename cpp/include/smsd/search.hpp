/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The searches behind Substructure and MCS, and the entry points that choose one.
 * Their parts are in search_types.hpp and detail/: batch_dispatch.hpp,
 * search_budget.hpp, search_compatibility.hpp, search_results.hpp,
 * vertex_search.hpp, pattern_search.hpp and match_graph.hpp.
 */
#pragma once

#include "smsd/common_atom_search.hpp"
#include "smsd/common_bond_search.hpp"
#include "smsd/mol_graph.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/stereo.hpp"
#include "smsd/time_budget.hpp"
#include "smsd/gpu_kernels.hpp"
#include "smsd/search_types.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <functional>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <thread>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "smsd/detail/batch_dispatch.hpp"
#include "smsd/detail/search_budget.hpp"
#include "smsd/detail/search_compatibility.hpp"
#include "smsd/detail/search_results.hpp"
#include "smsd/detail/vertex_search.hpp"
#include "smsd/detail/pattern_search.hpp"
#include "smsd/detail/match_graph.hpp"

namespace smsd {

namespace search_detail {

inline std::shared_ptr<const MolGraph> ownQuery(const MolGraph& query) {
    const std::lock_guard<std::recursive_mutex> lock(query.lazyCacheMutex());
    validateSearchGraphStorage(query);
    return std::make_shared<const MolGraph>(query);
}

inline void prepare(const MolGraph* query, const MolGraph& target, const ChemOptions& chemistry) {
    if (chemistry.ringFusionMode != ChemOptions::RingFusionMode::IGNORE) {
        if (query) query->ensureRingCounts();
        target.ensureRingCounts();
    }
}

/**
 * The best MCES of one orientation by the common-bond search: one mapping, or with `all` every
 * distinct best mapping up to `limit`, each extending `required` unless it is null, with the
 * atom and bond predicates of the other searches and the shared budget for time, cancellation
 * and nodes. The matches are sorted by atom mapping and then by bond mapping. A search that
 * holds more than `limit` best mappings keeps the first `limit` it found and reports the limit
 * reached. The result is optimal only when the search ran to its end. Graphs the search cannot
 * take, as they are not well formed or too large, give no match, with the limit reached.
 */
inline SearchResult commonBondSearch(const MolGraph& query, const MolGraph& target, const MCSOptions& options,
                                     size_t limit, bool all, const Match* required, Budget& budget) {
    prepare(&query, target, options.chemistry);
    const ChemOptions& chemistry = options.chemistry;
    const bool prepareRings = chemistry.completeRingsOnly && !budget.stopped()
        && common_bond_detail::CommonBondGraphs::supports(query, target);
    RingConstraint queryRings(query, prepareRings), targetRings(target, prepareRings);
    const bool pairRings = prepareRings && !queryRings.truncated() && !targetRings.truncated();
    const bool connected = options.connectivity == MCSConnectivity::Connected || options.maxFragments == 1;
    const bool tautomers = chemistry.tautomerAware;
    common_bond_detail::CommonBondRequest request;
    request.wanted = all ? (limit == SIZE_MAX ? limit : limit + 1) : 1;
    if (options.nodeLimit != 0) {
        const uint64_t spent = budget.nodesVisited();
        request.nodeCap = spent >= options.nodeLimit ? 0
            : static_cast<int64_t>(std::min<uint64_t>(options.nodeLimit - spent, uint64_t(INT64_MAX)));
    }
    request.minFragmentAtoms = options.minFragmentAtoms;
    request.maxFragments = options.maxFragments;
    request.chirality = chemistry.useChirality;
    request.bondStereo = chemistry.useBondStereo;
    request.rings = chemistry.completeRingsOnly;
    if (prepareRings) {
        request.queryCycles = &queryRings.cycles();
        request.targetCycles = &targetRings.cycles();
    }
    if (required && !required->atomMapping.empty()) {
        request.requiredTarget.assign(size_t(query.n), -1);
        for (const auto& atom : required->atomMapping) request.requiredTarget[atom.first] = atom.second;
        for (const auto& bond : required->bondMapping) request.requiredBonds.push_back({bond[0], bond[1]});
    }
    auto outcome = common_bond_detail::searchCommonBonds(query, target,
        [&](int queryAtom, int targetAtom) {
            return (!pairRings || queryRings.ringAtom(queryAtom) == targetRings.ringAtom(targetAtom))
                && atomCompatible(query, queryAtom, target, targetAtom, chemistry);
        },
        [&](int queryFirst, int querySecond, int targetFirst, int targetSecond) {
            return (!pairRings || queryRings.ringBond(queryFirst, querySecond)
                                     == targetRings.ringBond(targetFirst, targetSecond))
                && bondCompatible(query, queryFirst, querySecond, target, targetFirst, targetSecond, chemistry);
        },
        [&](const MolGraph& graph, int first, int second) {
            const int64_t key = bondCompatibilityKey(graph, first, second, tautomers);
            const RingConstraint& rings = &graph == &query ? queryRings : targetRings;
            return pairRings ? key * 2 + (rings.ringBond(first, second) ? 1 : 0) : key;
        },
        options.excludedTargetAtoms, connected, [&budget] { return budget.stopped(); }, request);
    budget.countNodes(outcome.nodes);
    if (!outcome.supported || outcome.capReached) budget.limitReached();
    SearchResult result;
    for (size_t index = 0; index < outcome.mappings.size() && result.matches.size() < limit; ++index)
        result.matches.push_back(makeMatch(outcome.mappings[index], std::move(outcome.bondLists[index])));
    std::sort(result.matches.begin(), result.matches.end(), [](const Match& first, const Match& second) {
        return std::tie(first.atomMapping, first.bondMapping) < std::tie(second.atomMapping, second.bondMapping);
    });
    const bool proved = outcome.supported && outcome.finished && budget.status() == SearchStatus::Completed;
    const bool overflow = all && outcome.mappings.size() > limit;
    if (proved && overflow) budget.limitReached();
    result.status = budget.status();
    result.optimal = proved;
    result.complete = proved && !overflow;
    result.score = result.matches.empty() ? 0.0 : double(result.matches.front().bondMapping.size());
    result.nodesVisited = budget.nodesVisited();
    result.elapsedMillis = budget.elapsedMillis();
    return result;
}

/**
 * Whether the common-atom search takes the search: a single best MCIS in any number of fragments,
 * with the constraints the common-bond search also leaves out. A connected MCIS stays with the
 * vertex search, which is faster there.
 */
inline bool commonAtomSearchApplies(const MCSOptions& options, size_t limit, bool all, const Match* required) {
    const ChemOptions& chemistry = options.chemistry;
    return !all && limit == 1 && !required && options.objective == MCSObjective::MCIS
        && options.connectivity == MCSConnectivity::Fragments && options.maxFragments == INT_MAX
        && options.atomWeights.empty() && options.nodeLimit == 0
        && !chemistry.useChirality && !chemistry.useBondStereo && !chemistry.completeRingsOnly
        && !chemistry.tautomerAware && options.minFragmentAtoms <= 2;
}

/**
 * The single best MCIS in fragments by the common-atom search, with the atom and bond
 * predicates of the other searches and the shared budget for time and cancellation. Nothing
 * when that search cannot take the graphs. The result is optimal and complete only when the
 * search ran to its end.
 */
inline std::optional<SearchResult> commonAtomSearch(const MolGraph& query, const MolGraph& target,
                                                    const MCSOptions& options, Budget& budget) {
    prepare(&query, target, options.chemistry);
    const ChemOptions& chemistry = options.chemistry;
    auto outcome = common_atom_detail::searchCommonAtoms(query, target,
        [&](int queryAtom, int targetAtom) { return atomCompatible(query, queryAtom, target, targetAtom, chemistry); },
        [&](int queryFirst, int querySecond, int targetFirst, int targetSecond) {
            return bondCompatible(query, queryFirst, querySecond, target, targetFirst, targetSecond, chemistry);
        },
        [](const MolGraph& graph, int first, int second) { return bondCompatibilityKey(graph, first, second); },
        options.excludedTargetAtoms, options.minFragmentAtoms, [&budget] { return budget.stopped(); });
    if (!outcome.supported) return std::nullopt;
    budget.countNodes(outcome.nodes);
    SearchResult result;
    if (!outcome.mapping.empty()) result.matches.push_back(makeMatch(outcome.mapping, std::move(outcome.bonds)));
    result.status = budget.status();
    result.optimal = result.complete = outcome.finished && result.status == SearchStatus::Completed;
    result.score = result.matches.empty() ? 0.0 : double(result.matches.front().atomMapping.size());
    result.nodesVisited = budget.nodesVisited();
    result.elapsedMillis = budget.elapsedMillis();
    return result;
}

inline SearchResult commonSubgraphOriented(const MolGraph& query, const MolGraph& target,
                                          const MCSOptions& options, size_t limit, bool all,
                                          const Match* required, Budget& budget, bool allowInjectiveBound = true) {
    Collector collector(&options, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    if (options.objective == MCSObjective::MCES) return commonBondSearch(query, target, options, limit, all, required, budget);
    if (commonAtomSearchApplies(options, limit, all, required))
        if (auto result = commonAtomSearch(query, target, options, budget)) return std::move(*result);
    if (all) {
        const auto optimum = commonSubgraphOriented(query, target, options, 1, false, required, budget, false);
        for (const auto& match : optimum.matches) collector.remember(match);
        if (optimum.optimal) collector.markOptimumKnown();
        if (budget.stopped()) return collector.finish(budget, optimum.optimal, false);
    }
    prepare(&query, target, options.chemistry);
    std::pair<bool, bool> proof;
    if (!all && !required && query.n > 0 && query.n <= target.n) {
        const QueryView view(query);
        VertexSearch full(view, target, options.chemistry, &options, true, false, budget, collector, -1, nullptr, false, 2048,
                          true, false, allowInjectiveBound && !all);
        const auto fullProof = full.run();
        if (!all && collector.hasMatch()
            && collector.best().atoms == query.n
            && std::all_of(options.atomWeights.begin(), options.atomWeights.end(), [](double weight) { return weight >= 0.0; }))
            return collector.finish(budget, true, true);
        if (fullProof.second && options.atomWeights.empty()
            && collector.best().atoms == query.n - 1) {
            int64_t degreeSum = 0, minimumDegree = query.n;
            SpareVectors<int>::Lease seenLease(static_cast<size_t>(query.n), -1);
            auto& seen = seenLease.items;
            size_t work = 0;
            for (int atom = 0; atom < query.n; ++atom) {
                if (budget.stopped()) break;
                int64_t degree = 0;
                for (int neighbor : query.neighbors[atom]) {
                    if ((work++ & 255u) == 0 && budget.stopped()) break;
                    if (neighbor != atom && seen[neighbor] != atom) {
                        seen[neighbor] = atom;
                        ++degree;
                    }
                }
                degreeSum += degree;
                minimumDegree = std::min(minimumDegree, degree);
            }
            // Any n-1-atom induced graph omits all bonds incident to its omitted atom.
            if (!budget.stopped() && collector.best().bonds == degreeSum / 2 - minimumDegree)
                return collector.finish(budget, true, true);
        }
    }
    const QueryView view(query);
    VertexSearch search(view, target, options.chemistry, &options, true, all, budget, collector, -1, required,
                        true, 0, true, allowInjectiveBound && !all, allowInjectiveBound && !all);
    proof = search.run();
    return collector.finish(budget, proof.first, proof.second);
}

/**
 * Whether proof needs a complete ring source: the relevant cycles for completeRingsOnly,
 * and the ring-count perception for STRICT fusion.
 */
inline bool needsRingProof(const ChemOptions& chemistry) {
    return chemistry.completeRingsOnly || chemistry.ringFusionMode == ChemOptions::RingFusionMode::STRICT;
}

/** Whether the ring set used by the chemistry was truncated. */
inline bool ringProofTruncated(const MolGraph& graph, const ChemOptions& chemistry) {
    bool truncated = false;
    if (chemistry.completeRingsOnly) (void)recentRelevantCycles(graph, &truncated);
    return truncated || (chemistry.ringFusionMode == ChemOptions::RingFusionMode::STRICT
                         && graph.ringCountsTruncated());
}

/** Retain the search results without claiming proof from truncated ring perception. */
inline void qualifyRingProof(SearchResult& result, const MolGraph& query, const MolGraph& target,
                             const ChemOptions& chemistry) {
    if (!needsRingProof(chemistry)
        || (!result.optimal && !result.complete && result.status != SearchStatus::Completed)
        || (!ringProofTruncated(query, chemistry) && !ringProofTruncated(target, chemistry))) return;
    if (result.status == SearchStatus::Completed || result.status == SearchStatus::LimitReached)
        result.status = SearchStatus::LimitReached;
    result.optimal = result.complete = false;
}

/** Whether each atom lists every other atom once, with bonds in both directions. */
inline bool completeAtomTopology(const MolGraph& graph, Budget& budget) {
    if (graph.n < 0 || graph.neighbors.size() != static_cast<size_t>(graph.n)) return false;
    for (int atom = 0; atom < graph.n; ++atom) {
        if ((atom & 255) == 0 && budget.stopped()) return false;
        if (graph.neighbors[atom].size() != static_cast<size_t>(graph.n - 1)) return false;
    }
    std::vector<int> seen(static_cast<size_t>(graph.n), -1);
    uint64_t work = 0;
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom]) {
            if ((work++ & 255) == 0 && budget.stopped()) return false;
            if (neighbor < 0 || neighbor >= graph.n || neighbor == atom || seen[neighbor] == atom) return false;
            seen[neighbor] = atom;
            if (!graph.hasBond(atom, neighbor) || !graph.hasBond(neighbor, atom)) return false;
        }
    return true;
}

inline SearchResult commonSubgraph(const MolGraph& query, const MolGraph& target,
                                  const MCSOptions& options, size_t limit, bool all, const Match* required = nullptr,
                                  Budget* sharedBudget = nullptr) {
    validateSearchGraphShape(query);
    validateSearchGraphShape(target);
    validateOptions(options, query.n);
    for (int atom : options.excludedTargetAtoms)
        if (atom < 0 || atom >= target.n) throw std::out_of_range("Excluded target atom index out of range");
    Budget localBudget(options, query.n, target.n, true);
    Budget& budget = sharedBudget ? *sharedBudget : localBudget;
    RecursiveScope recursiveScope(budget);
    if (budget.stopped()) {
        Collector collector(&options, limit, all);
        return collector.finish(budget, false, false);
    }
    if (!validateSearchGraphStorage(query, [&budget] { return budget.stoppedSampled(); })
        || !validateSearchGraphStorage(target, [&budget] { return budget.stoppedSampled(); })) {
        Collector collector(&options, limit, all);
        return collector.finish(budget, false, false);
    }
    bool reverse = false;
    if (options.objective == MCSObjective::MCES && !required && options.excludedTargetAtoms.empty()) {
        size_t queryIncidences = 0, targetIncidences = 0;
        for (const auto& neighbors : query.neighbors) queryIncidences += neighbors.size();
        for (const auto& neighbors : target.neighbors) targetIncidences += neighbors.size();
        reverse = (options.connectivity == MCSConnectivity::Connected || options.maxFragments == 1)
                ? query.n > target.n : queryIncidences > targetIncidences;
    } else if (options.objective == MCSObjective::MCIS && !all && limit == 1 && !required
            && options.excludedTargetAtoms.empty() && options.atomWeights.empty() && options.nodeLimit == 0
            && options.connectivity == MCSConnectivity::Connected && options.maxFragments == INT_MAX
            && options.minFragmentAtoms == 1 && !options.chemistry.useChirality
            && !options.chemistry.useBondStereo && !options.chemistry.completeRingsOnly
            && !options.chemistry.tautomerAware && !options.chemistry.ringMatchesRingOnly
            && options.chemistry.ringFusionMode == ChemOptions::RingFusionMode::IGNORE
            && !hasStereoTags(query) && !hasStereoTags(target)) {
        reverse = completeAtomTopology(target, budget) && !completeAtomTopology(query, budget);
        if (budget.stopped()) {
            Collector collector(&options, limit, all);
            return collector.finish(budget, false, false);
        }
    }
    if (reverse) {
        auto result = commonSubgraphOriented(target, query, options, limit, all, nullptr, budget);
        for (auto& match : result.matches) {
            std::vector<int> mapping(query.n, -1);
            for (const auto& atom : match.atomMapping) mapping[atom.second] = atom.first;
            std::vector<std::array<int, 4>> bonds;
            bonds.reserve(match.bondMapping.size());
            for (const auto& bond : match.bondMapping) {
                if (bond[2] < bond[3]) bonds.push_back({bond[2], bond[3], bond[0], bond[1]});
                else bonds.push_back({bond[3], bond[2], bond[1], bond[0]});
            }
            match = makeMatch(mapping, std::move(bonds));
        }
        std::sort(result.matches.begin(), result.matches.end(), [](const Match& first, const Match& second) {
            return std::tie(first.atomMapping, first.bondMapping) < std::tie(second.atomMapping, second.bondMapping);
        });
        result.elapsedMillis = budget.elapsedMillis();
        qualifyRingProof(result, query, target, options.chemistry);
        return result;
    }
    auto result = commonSubgraphOriented(query, target, options, limit, all, required, budget);
    qualifyRingProof(result, query, target, options.chemistry);
    return result;
}

/**
 * Whether the target has, for each element and each minimum degree, at least as many atoms as
 * the query. An embedding maps query atoms one to one onto target atoms of their element with
 * at least as many neighbours, so a shortfall proves that none exists.
 */
inline bool atomCountsAdmitEmbedding(const MolGraph& query, const MolGraph& target,
                                    const ChemOptions& chemistry, Budget& budget) {
    if (!chemistry.matchAtomType) return true;
    // Degrees from the last bin up share it, which keeps the test necessary. Queries with more
    // elements than slots are not screened.
    constexpr int elements = 128, degrees = 8, slots = 16;
    std::array<int8_t, elements> slot;
    slot.fill(-1);
    std::array<std::array<int, 2 * degrees>, slots> counts;
    int used = 0;
    SpareVectors<int>::Lease seenLease(static_cast<size_t>(query.n), -1);
    auto& seen = seenLease.items;
    size_t work = 0;
    for (int atom = 0; atom < query.n; ++atom) {
        if ((atom & 255) == 0 && budget.stopped()) return false;
        const int element = query.atomicNum[atom];
        if (element < 0 || element >= elements) return true;
        if (slot[element] < 0) {
            if (used == slots) return true;
            counts[used].fill(0);
            slot[element] = static_cast<int8_t>(used++);
        }
        int degree = 0;
        for (int neighbor : query.neighbors[atom]) {
            if ((work++ & 255u) == 0 && budget.stopped()) return false;
            if (neighbor != atom && seen[neighbor] != atom) {
                seen[neighbor] = atom;
                ++degree;
            }
        }
        ++counts[slot[element]][std::min(degree, degrees - 1)];
    }
    for (int atom = 0; atom < target.n; ++atom) {
        if ((atom & 255) == 0 && budget.stopped()) return false;
        const int element = target.atomicNum[atom];
        if (element >= 0 && element < elements && slot[element] >= 0)
            ++counts[slot[element]][degrees + std::min(target.neighbors[atom].size(), static_cast<size_t>(degrees - 1))];
    }
    for (int index = 0; index < used; ++index) {
        int queryAtoms = 0, targetAtoms = 0;
        for (int degree = degrees - 1; degree >= 0; --degree) {
            queryAtoms += counts[index][degree];
            targetAtoms += counts[index][degrees + degree];
            if (queryAtoms > targetAtoms) return false;
        }
    }
    return true;
}

inline SearchResult patternSubstructure(const SmartsQuery& pattern, const MolGraph& target,
                                        const SubstructureOptions& options, size_t limit, bool all,
                                        Budget* sharedBudget, const PatternPlan* plan) {
    validateSearchGraphShape(target);
    Budget localBudget(options, pattern.atomCount(), target.n, false);
    Budget& budget = sharedBudget ? *sharedBudget : localBudget;
    RecursiveScope recursiveScope(budget);
    Collector collector(nullptr, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    if (!validateSearchGraphStorage(target, [&budget] { return budget.stoppedSampled(); }))
        return collector.finish(budget, false, false);
    std::optional<PatternPlan> localPlan;
    if (!plan) plan = &localPlan.emplace(pattern);
    if (pattern.atomCount() <= target.n && plan->numericRingMembership()) {
        if (!plan->elementsAvailable(target, budget)) {
            const bool stopped = budget.stopped();
            return collector.finish(budget, !stopped, !stopped);
        }
        const bool incomplete = smsd::ringsTruncated(target);
        if (budget.stopped()) return collector.finish(budget, false, false);
        if (incomplete) {
            auto result = collector.finish(budget, false, false);
            result.status = SearchStatus::LimitReached;
            return result;
        }
    }
    prepare(nullptr, target, options.chemistry);
    thread_local RecursiveMemo::Pool reusablePool;
    RecursiveMemo::Pool transientPool;
    RecursiveMemo memo(target, reusablePool.busy ? transientPool : reusablePool, plan);
    RecursiveMemoScope memoScope(memo);
    thread_local PatternWorkspace reusable;
    PatternWorkspace transient;
    PatternWorkspace& workspace = reusable.busy ? transient : reusable;
    const PatternWorkspace::Hold hold(workspace);
    PatternSearch search(pattern, *plan, target, options.chemistry, options.induced, all, limit, budget, workspace);
    const auto proof = search.run();
    auto result = search.finish(proof.first, proof.second);
    if (!search.ringIndependentProof())
        qualifyRingProof(result, target, target, options.chemistry);
    return result;
}

/**
 * A pattern plan, when given, must have been built from the same pattern. A prepared query
 * already holds the ring counts the options need.
 */
inline SearchResult substructure(const MolGraph* query, const SmartsQuery* pattern,
                                const MolGraph& target, const SubstructureOptions& options,
                                size_t limit, bool all, Budget* sharedBudget = nullptr,
                                const PatternPlan* plan = nullptr, bool queryPrepared = false,
                                bool screenedOut = false) {
    if (!query) return patternSubstructure(*pattern, target, options, limit, all, sharedBudget, plan);
    validateSearchGraphShape(*query);
    validateSearchGraphShape(target);
    const QueryView view(*query);
    Budget localBudget(options, view.size(), target.n, false);
    Budget& budget = sharedBudget ? *sharedBudget : localBudget;
    RecursiveScope recursiveScope(budget);
    Collector collector(nullptr, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    if (!validateSearchGraphStorage(*query, [&budget] { return budget.stoppedSampled(); })
        || !validateSearchGraphStorage(target, [&budget] { return budget.stoppedSampled(); }))
        return collector.finish(budget, false, false);
    // The atom counts reject a pair before any ring perception, which they do not read.
    if (screenedOut || (query && ((view.size() > target.n && target.n >= 0)
        || !atomCountsAdmitEmbedding(*query, target, options.chemistry, budget)))) {
        const bool completed = !budget.stopped();
        return collector.finish(budget, completed, completed);
    }
    prepare(queryPrepared ? nullptr : query, target, options.chemistry);
    VertexSearch search(view, target, options.chemistry, nullptr, options.induced, all, budget, collector);
    const auto proof = search.run();
    auto result = collector.finish(budget, proof.first, proof.second);
    qualifyRingProof(result, *query, target, options.chemistry);
    return result;
}

} // namespace search_detail

/**
 * The hydrogens of a nitrogen, oxygen or sulfur atom, which a tautomer shift can move: the
 * stored hydrogen count, or the valence model's count when none was supplied, as the Java
 * `hydrogenCount` reads it. Other elements carry none.
 */
inline int mobilisableProtons(const MolGraph& g, int atom) {
    const int z = g.atomicNum[atom];
    if (z != 7 && z != 8 && z != 16) return 0;
    return g.implicitHydrogenCount(atom);
}

inline bool validateTautomerConsistency(
    const MolGraph& g1, const MolGraph& g2, const std::map<int,int>& mcs) {
    if (mcs.empty()) return true;
    // Ensure tautomer classes are computed
    if (g1.tautomerClass.empty() || g2.tautomerClass.empty()) return true;

    // Group matched atoms by tautomer class and sum protons per side
    std::unordered_map<int, int> querySideProtons, targetSideProtons;
    for (auto& [qi, ti] : mcs) {
        if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
        int tc1 = g1.tautomerClass[qi];
        int tc2 = g2.tautomerClass[ti];
        if (tc1 < 0 && tc2 < 0) continue;
        // Use query-side class as group key (or target-side if query has none)
        int groupKey = (tc1 >= 0) ? tc1 : tc2 + 100000;
        querySideProtons[groupKey]  += mobilisableProtons(g1, qi);
        targetSideProtons[groupKey] += mobilisableProtons(g2, ti);
    }
    // Check each group for proton conservation
    for (auto& [grp, qProtons] : querySideProtons) {
        auto it = targetSideProtons.find(grp);
        int tProtons = (it != targetSideProtons.end()) ? it->second : 0;
        if (qProtons != tProtons) return false;
    }
    return true;
}


} // namespace smsd

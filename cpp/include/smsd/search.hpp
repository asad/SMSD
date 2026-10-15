/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#pragma once

#include "mol_graph.hpp"
#include "smarts_parser.hpp"
#include "smiles_parser.hpp"
#include "stereo.hpp"
#include "time_budget.hpp"
#include "gpu_kernels.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <functional>
#include <exception>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <thread>
#include <stdexcept>
#include <string>
#include <utility>
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

namespace search_detail {

template<class Query, class Targets>
inline std::vector<SearchResult> matchMany(const Query& query, const Targets& targets, int numThreads) {
    if (numThreads < 0) throw std::invalid_argument("Thread count must be nonnegative");
    std::vector<SearchResult> results(targets.size());
    if (results.empty()) return results;
    const size_t count = std::min(results.size(), static_cast<size_t>(numThreads == 0
        ? std::max(1u, std::thread::hardware_concurrency()) : numThreads));
    std::atomic<size_t> next{0};
    std::exception_ptr failure;
    std::mutex failureMutex;
    auto work = [&] {
        try {
            for (size_t index = next.fetch_add(1); index < targets.size(); index = next.fetch_add(1))
                results[index] = query.match(targets[index]);
        } catch (...) {
            const std::lock_guard<std::mutex> lock(failureMutex);
            if (!failure) failure = std::current_exception();
            next.store(targets.size());
        }
    };
    std::vector<std::thread> workers;
    for (size_t thread = 1; thread < count; ++thread) workers.emplace_back(work);
    work();
    for (auto& worker : workers) worker.join();
    if (failure) std::rethrow_exception(failure);
    return results;
}

/** Every stage and recursive query in a call shares this budget. */
class Budget {
    using Clock = std::chrono::steady_clock;
    Clock::time_point started_ = Clock::now();
    Clock::time_point deadline_;
    uint64_t nodeLimit_;
    std::shared_ptr<CancellationToken> cancellation_;
    SearchStatus status_ = SearchStatus::Completed;
    uint64_t nodes_ = 0;
public:
    Budget(const SearchOptions& options, int queryAtoms, int targetAtoms, bool commonSubgraph)
        : nodeLimit_(options.nodeLimit), cancellation_(options.cancellation) {
        const int64_t milliseconds = options.timeoutMs < 0
            ? (commonSubgraph ? std::min<int64_t>(30000, 500 + int64_t(queryAtoms) * targetAtoms * 2) : 10000)
            : options.timeoutMs;
        deadline_ = detail::steadyDeadline(milliseconds, started_);
    }
    bool stopped() {
        if (status_ != SearchStatus::Completed) return true;
        if (cancellation_ && cancellation_->isCancelled()) status_ = SearchStatus::Cancelled;
        else if (Clock::now() >= deadline_) status_ = SearchStatus::TimedOut;
        return status_ != SearchStatus::Completed;
    }
    bool visit() {
        if (stopped()) return false;
        if (nodeLimit_ != 0 && nodes_ >= nodeLimit_) {
            status_ = SearchStatus::LimitReached;
            return false;
        }
        ++nodes_;
        return true;
    }
    SearchStatus status() const { return status_; }
    uint64_t nodesVisited() const { return nodes_; }
    int64_t elapsedMillis() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - started_).count();
    }
};

inline thread_local Budget* recursiveBudget = nullptr;
inline bool recursiveMatch(const SmartsQuery& pattern, const MolGraph& target, int anchor);

class RecursiveScope {
    Budget* previousBudget_ = recursiveBudget;
    smarts_detail::RecursiveMatcher previousMatcher_ = smarts_detail::recursiveMatcher;
public:
    explicit RecursiveScope(Budget& budget) {
        recursiveBudget = &budget;
        smarts_detail::recursiveMatcher = recursiveMatch;
    }
    ~RecursiveScope() {
        recursiveBudget = previousBudget_;
        smarts_detail::recursiveMatcher = previousMatcher_;
    }
};

inline std::shared_ptr<const MolGraph> ownQuery(const MolGraph& query) {
    const std::lock_guard<std::recursive_mutex> lock(query.lazyCacheMutex());
    return std::make_shared<const MolGraph>(query);
}

inline bool atomCompatible(const MolGraph& gq, int qi,
                            const MolGraph& gt, int tj,
                            const ChemOptions& C) {
    // Tautomerism changes proton positions and bond orders, not elements.
    if (C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj]) return false;
    auto property = [](const auto& values, int atom) { return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : 0; };
    if (C.matchFormalCharge && property(gq.formalCharge, qi) != property(gt.formalCharge, tj)) return false;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT && gq.aromatic[qi] != gt.aromatic[tj]) return false;
    if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
    if (C.matchIsotope) {
        int qm = property(gq.massNumber, qi), tm = property(gt.massNumber, tj);
        if (qm != 0 && tm != 0 && qm != tm) return false;
    }
    if (C.useChirality) {
        const int queryTag = property(gq.tetraChirality, qi), targetTag = property(gt.tetraChirality, tj);
        const bool queryLigandsValid = gq.degree[qi] == 3 || gq.degree[qi] == 4;
        const bool targetLigandsValid = gt.degree[tj] == 3 || gt.degree[tj] == 4;
        if (queryTag != 0 && targetTag != 0 && (!queryLigandsValid || !targetLigandsValid)
            && queryTag != targetTag) return false;
    }
    if (C.ringFusionMode != ChemOptions::RingFusionMode::IGNORE && gq.ring[qi] && gt.ring[tj]) {
        if (C.ringFusionMode == ChemOptions::RingFusionMode::STRICT) {
            if (property(gq.ringCount, qi) != property(gt.ringCount, tj)) return false;
        }
    }
    return true;
}

inline bool bondCompatible(const MolGraph& g1, int qi, int qk,
                            const MolGraph& g2, int tj, int tk,
                            const ChemOptions& C) {
    int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
    if (qOrd == 0 || tOrd == 0) return false;
    // Tautomer-aware: relax bond ORDER (single/double interchangeable)
    // but still enforce aromaticity, stereo, and ring constraints.
    bool tautBondRelax = C.tautomerAware
        && !g1.tautomerClass.empty() && !g2.tautomerClass.empty()
        && g1.tautomerClass[qi] != -1 && g1.tautomerClass[qk] != -1
        && g2.tautomerClass[tj] != -1 && g2.tautomerClass[tk] != -1;
    // Strict aromaticity: check bond aromaticity directly (no ring guard)
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT
        && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk)) return false;
    if (C.ringMatchesRingOnly && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk)) return false;
    if (tautBondRelax) return true;  // bond order relaxed for tautomeric bonds
    if (C.matchBondOrder == ChemOptions::BondOrderMode::ANY) return true;
    if (qOrd == tOrd) return true;
    if (C.matchBondOrder == ChemOptions::BondOrderMode::LOOSE) return true;
    if (C.aromaticityMode == ChemOptions::AromaticityMode::FLEXIBLE) {
        bool qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
        if ((qa && ta) || (qa && (tOrd == 1 || tOrd == 2)) || (ta && (qOrd == 1 || qOrd == 2)))
            return true;
        // Resonance-aware: allow single/double interchange between atoms that
        // are both in a conjugated (aromatic) ring system.  This catches Kekulé
        // form mismatches (e.g. pyridone drawn as C=O vs C-O) without requiring
        // explicit tautomer flags.  (v6.8.0)
        if ((g1.aromatic[qi] && g1.aromatic[qk]) &&
            (g2.aromatic[tj] && g2.aromatic[tk]) &&
            (qOrd == 1 || qOrd == 2) && (tOrd == 1 || tOrd == 2))
            return true;
    }
    return false;
}

template <typename Mapping>
inline bool tetraParityCompatible(const MolGraph& gq, int qi,
                                  const MolGraph& gt, int tj,
                                  const Mapping& q2t) {
    int qs = qi < static_cast<int>(gq.tetraChirality.size()) ? gq.tetraChirality[qi] : 0;
    int ts = tj < static_cast<int>(gt.tetraChirality.size()) ? gt.tetraChirality[tj] : 0;
    if (qs == 0 || ts == 0) return true;

    std::array<int, 4> queryLigands{}, targetLigands{};
    if (!stereo_detail::tetrahedralLigands(gq.neighbors[qi], qi, queryLigands)
        || !stereo_detail::tetrahedralLigands(gt.neighbors[tj], tj, targetLigands)) return true;
    return stereo_detail::tetrahedralParityMatches(qs, queryLigands, ts, targetLigands, q2t);
}

template <typename Mapping>
inline bool mappingStereoCompatible(const MolGraph& gq, const MolGraph& gt,
                                    const ChemOptions& opts, const Mapping& q2t) {
    if (!opts.useChirality && !opts.useBondStereo) return true;
    for (int queryAtom = 0; queryAtom < gq.n; ++queryAtom) {
        const int targetAtom = q2t[queryAtom];
        if (targetAtom < 0) continue;
        if (opts.useChirality && !tetraParityCompatible(gq, queryAtom, gt, targetAtom, q2t)) return false;
        if (opts.useBondStereo) {
            for (int neighbor : gq.neighbors[queryAtom])
                if (neighbor > queryAtom && !stereo_detail::doubleBondParityMatches(gq, queryAtom, neighbor, gt, q2t)) return false;
        }
    }
    return true;
}


struct QueryView {
    const MolGraph* molecule = nullptr;
    const SmartsQuery* pattern = nullptr;
    std::unordered_map<int64_t, int> patternBondIndex;
    explicit QueryView(const MolGraph& graph) : molecule(&graph) {}
    explicit QueryView(const SmartsQuery& query) : pattern(&query) {
        for (size_t bond = 0; bond < query.bonds.size(); ++bond)
            patternBondIndex.emplace(MolGraph::bondKey(query.bonds[bond].from, query.bonds[bond].to), static_cast<int>(bond));
    }
    int size() const { return molecule ? molecule->n : pattern->atomCount(); }
    const std::vector<int>& neighbors(int atom) const {
        return molecule ? molecule->neighbors[atom] : pattern->neighbors[atom];
    }
    bool adjacent(int first, int second) const {
        return molecule ? molecule->hasBond(first, second)
            : patternBondIndex.count(MolGraph::bondKey(first, second)) != 0;
    }
    bool atomMatches(int queryAtom, const MolGraph& target, int targetAtom, const ChemOptions& chemistry) const {
        if (molecule) return atomCompatible(*molecule, queryAtom, target, targetAtom, chemistry);
        const auto truth = smarts_detail::evalAtomExprWithStereo(pattern->atoms[queryAtom].expr, target, targetAtom, [](int tag) {
            return tag == 0 ? smarts_detail::StereoTruth::TRUE_VALUE : smarts_detail::StereoTruth::UNKNOWN;
        });
        return truth != smarts_detail::StereoTruth::FALSE_VALUE;
    }
    bool bondMatches(int first, int second, const MolGraph& target,
                     int mappedFirst, int mappedSecond, const ChemOptions& chemistry) const {
        if (molecule) return bondCompatible(*molecule, first, second, target, mappedFirst, mappedSecond, chemistry);
        const auto bond = patternBondIndex.find(MolGraph::bondKey(first, second));
        return bond != patternBondIndex.end() && pattern->bondCompat(bond->second, target, mappedFirst, mappedSecond);
    }
};

inline int64_t weightedScore(double weight) {
    const double score = weight * 1000.0;
    if (!std::isfinite(score) || score <= double(INT_MIN) - 1.0 || score >= double(INT_MAX) + 1.0)
        throw std::invalid_argument("Atom weights exceed the supported millipoint score range");
    return static_cast<int64_t>(score);
}

inline void validateOptions(const MCSOptions& options, int queryAtoms) {
    if (options.minFragmentAtoms < 1 || options.maxFragments < 1)
        throw std::invalid_argument("Fragment limits must be positive");
    if (options.objective == MCSObjective::MCES && !options.atomWeights.empty())
        throw std::invalid_argument("Atom weights are supported only for MCIS");
    if (options.atomWeights.empty()) return;
    if (options.atomWeights.size() < static_cast<size_t>(queryAtoms))
        throw std::invalid_argument("Atom weights must cover every query atom");
    double positive = 0.0, negative = 0.0;
    for (int atom = 0; atom < queryAtoms; ++atom) {
        const double weight = options.atomWeights[atom];
        if (!std::isfinite(weight)) throw std::invalid_argument("Atom weights must be finite");
        if (weight > 0.0) positive += weight;
        else negative += weight;
    }
    (void)weightedScore(positive);
    (void)weightedScore(negative);
}

struct Score {
    int64_t primary = 0;
    int atoms = 0;
    int bonds = 0;
    bool operator<(const Score& other) const {
        return std::tie(primary, atoms, bonds) < std::tie(other.primary, other.atoms, other.bonds);
    }
    bool operator==(const Score& other) const {
        return primary == other.primary && atoms == other.atoms && bonds == other.bonds;
    }
};

inline Score scoreMatch(const Match& match, const MCSOptions& options) {
    Score score;
    const int atoms = static_cast<int>(match.atomMapping.size());
    const int bonds = static_cast<int>(match.bondMapping.size());
    if (options.objective == MCSObjective::MCES) {
        score.primary = bonds;
        score.atoms = atoms;
    } else {
        double weight = 0.0;
        if (!options.atomWeights.empty())
            for (const auto& atom : match.atomMapping) weight += options.atomWeights[atom.first];
        score.primary = options.atomWeights.empty() ? atoms : weightedScore(weight);
        score.atoms = atoms;
        score.bonds = bonds;
    }
    return score;
}

inline Match makeMatch(const std::vector<int>& mapping, std::vector<std::array<int, 4>> bonds) {
    Match result;
    for (size_t atom = 0; atom < mapping.size(); ++atom)
        if (mapping[atom] >= 0) result.atomMapping.emplace(static_cast<int>(atom), mapping[atom]);
    std::sort(bonds.begin(), bonds.end());
    result.bondMapping = std::move(bonds);
    std::vector<std::vector<int>> adjacent(mapping.size());
    for (const auto& bond : result.bondMapping) {
        adjacent[bond[0]].push_back(bond[1]);
        adjacent[bond[1]].push_back(bond[0]);
    }
    std::vector<bool> visited(mapping.size(), false);
    for (const auto& atom : result.atomMapping) {
        if (visited[atom.first]) continue;
        std::vector<int> component{atom.first};
        visited[atom.first] = true;
        for (size_t position = 0; position < component.size(); ++position)
            for (int neighbor : adjacent[component[position]])
                if (!visited[neighbor]) { visited[neighbor] = true; component.push_back(neighbor); }
        std::sort(component.begin(), component.end());
        result.fragments.push_back(std::move(component));
    }
    return result;
}

inline std::vector<std::array<int, 4>> mappedBonds(const QueryView& query, const MolGraph& target,
                                                const std::vector<int>& mapping,
                                                const ChemOptions& chemistry) {
    std::vector<std::array<int, 4>> bonds;
    for (int atom = 0; atom < query.size(); ++atom) {
        if (mapping[atom] < 0) continue;
        for (int neighbor : query.neighbors(atom))
            if (atom < neighbor && mapping[neighbor] >= 0
                && query.bondMatches(atom, neighbor, target, mapping[atom], mapping[neighbor], chemistry))
                bonds.push_back({atom, neighbor, mapping[atom], mapping[neighbor]});
    }
    return bonds;
}

class RingConstraint {
    const MolGraph& graph_;
    std::vector<std::vector<int>> cycles_;
public:
    RingConstraint(const MolGraph& graph, bool enabled) : graph_(graph) {
        if (enabled) cycles_ = graph.computeRelevantCycles();
    }
    bool accepts(const std::set<int>& atoms, const std::set<int64_t>& bonds) const {
        std::set<int> coveredAtoms;
        std::set<int64_t> coveredBonds;
        for (const auto& cycle : cycles_) {
            bool complete = true;
            for (size_t atom = 0; atom < cycle.size() && complete; ++atom)
                complete = atoms.count(cycle[atom]) != 0
                    && bonds.count(MolGraph::bondKey(cycle[atom], cycle[(atom + 1) % cycle.size()])) != 0;
            if (!complete) continue;
            coveredAtoms.insert(cycle.begin(), cycle.end());
            for (size_t atom = 0; atom < cycle.size(); ++atom)
                coveredBonds.insert(MolGraph::bondKey(cycle[atom], cycle[(atom + 1) % cycle.size()]));
        }
        for (int atom : atoms)
            if (graph_.ring[atom] && coveredAtoms.count(atom) == 0) return false;
        for (int64_t bond : bonds) {
            const int first = static_cast<int>(bond >> 32), second = static_cast<int>(uint32_t(bond));
            if (graph_.bondInRing(first, second) && coveredBonds.count(bond) == 0) return false;
        }
        return true;
    }
};

inline bool selectedStereoCompatible(const MolGraph& query, const MolGraph& target,
                                     const ChemOptions& chemistry, const Match& match) {
    if (!chemistry.useChirality && !chemistry.useBondStereo) return true;
    std::vector<int> mapping(query.n, -1);
    for (const auto& atom : match.atomMapping) mapping[atom.first] = atom.second;
    std::set<int64_t> selected;
    for (const auto& bond : match.bondMapping) selected.insert(MolGraph::bondKey(bond[0], bond[1]));
    struct LocalMapping {
        const std::vector<int>& mapping;
        const std::set<int64_t>& selected;
        int first, second;
        int operator[](int atom) const {
            return atom == first || atom == second || selected.count(MolGraph::bondKey(first, atom))
                || (second >= 0 && selected.count(MolGraph::bondKey(second, atom))) ? mapping[atom] : -1;
        }
    };
    if (chemistry.useChirality)
        for (const auto& atom : match.atomMapping)
            if (!tetraParityCompatible(query, atom.first, target, atom.second,
                    LocalMapping{mapping, selected, atom.first, -1})) return false;
    if (chemistry.useBondStereo)
        for (const auto& bond : match.bondMapping)
            if (!stereo_detail::doubleBondParityMatches(query, bond[0], bond[1], target,
                    LocalMapping{mapping, selected, bond[0], bond[1]})) return false;
    return true;
}

inline bool completeRings(const Match& match, const RingConstraint& query, const RingConstraint& target) {
    std::set<int> queryAtoms, targetAtoms;
    std::set<int64_t> queryBonds, targetBonds;
    for (const auto& atom : match.atomMapping) { queryAtoms.insert(atom.first); targetAtoms.insert(atom.second); }
    for (const auto& bond : match.bondMapping) {
        queryBonds.insert(MolGraph::bondKey(bond[0], bond[1]));
        targetBonds.insert(MolGraph::bondKey(bond[2], bond[3]));
    }
    return query.accepts(queryAtoms, queryBonds) && target.accepts(targetAtoms, targetBonds);
}

class Collector {
    const MCSOptions* options_;
    size_t limit_;
    bool all_;
    Score best_{};
    std::set<std::pair<std::map<int, int>, std::vector<std::array<int, 4>>>> keys_;
    std::vector<Match> matches_;
    bool truncated_ = false;
    bool optimumKnown_ = false;
public:
    Collector(const MCSOptions* options, size_t limit, bool all)
        : options_(options), limit_(limit), all_(all) {
        if (limit == 0) throw std::invalid_argument("Result limit must be positive");
    }
    bool remember(Match match) {
        const Score score = options_ ? scoreMatch(match, *options_)
            : Score{static_cast<int64_t>(match.atomMapping.size()), static_cast<int>(match.atomMapping.size()), 0};
        if (score < best_) return false;
        if (best_ < score) { best_ = score; matches_.clear(); keys_.clear(); truncated_ = false; optimumKnown_ = false; }
        auto key = std::make_pair(match.atomMapping, match.bondMapping);
        if (keys_.count(key)) return false;
        if (matches_.size() < limit_) {
            keys_.insert(std::move(key));
            matches_.push_back(std::move(match));
        } else if (all_) truncated_ = true;
        return true;
    }
    const Score& best() const { return best_; }
    bool hasMatch() const { return !matches_.empty(); }
    bool truncated() const { return truncated_; }
    void markOptimumKnown() { optimumKnown_ = true; }
    bool optimumKnown() const { return optimumKnown_; }
    bool needs(const Score& upper) const { return all_ ? !(upper < best_) : best_ < upper; }
    SearchResult finish(Budget& budget, bool proven, bool exhausted) {
        SearchResult result;
        result.matches = std::move(matches_);
        std::sort(result.matches.begin(), result.matches.end(), [](const Match& first, const Match& second) {
            return std::tie(first.atomMapping, first.bondMapping) < std::tie(second.atomMapping, second.bondMapping);
        });
        result.status = budget.status();
        if (result.status == SearchStatus::Completed && truncated_) result.status = SearchStatus::LimitReached;
        result.optimal = proven;
        result.complete = exhausted && !truncated_ && result.status == SearchStatus::Completed;
        result.score = options_ && !options_->atomWeights.empty() ? best_.primary / 1000.0 : double(best_.primary);
        result.nodesVisited = budget.nodesVisited();
        result.elapsedMillis = budget.elapsedMillis();
        return result;
    }
};

/** Atom transpositions are used only while selecting one optimal witness. */
inline std::vector<int> transpositionClasses(const MolGraph& graph, Budget& budget,
                                             const std::vector<double>* weights = nullptr) {
    for (int tag : graph.tetraChirality) if (tag != 0) return {};
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom]) if (graph.dbStereo(atom, neighbor) != 0) return {};
    auto property = [](const auto& values, int atom, int fallback = 0) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : fallback;
    };
    std::vector<int> classes(graph.n);
    std::iota(classes.begin(), classes.end(), 0);
    for (int first = 0; first < graph.n && !budget.stopped(); ++first)
        for (int second = first + 1; second < graph.n; ++second) {
            if (classes[second] != second || graph.atomicNum[first] != graph.atomicNum[second]
                || property(graph.formalCharge, first) != property(graph.formalCharge, second)
                || property(graph.massNumber, first) != property(graph.massNumber, second)
                || property(graph.hydrogenCount, first) != property(graph.hydrogenCount, second)
                || property(graph.atomClass, first) != property(graph.atomClass, second)
                || graph.ring[first] != graph.ring[second] || graph.aromatic[first] != graph.aromatic[second]
                || graph.degree[first] != graph.degree[second] || property(graph.ringCount, first) != property(graph.ringCount, second)
                || property(graph.tautomerClass, first, -1) != property(graph.tautomerClass, second, -1)
                || (weights && !weights->empty() && (*weights)[first] != (*weights)[second])) continue;
            bool equivalent = true;
            for (int neighbor : graph.neighbors[first]) {
                if (neighbor == second) continue;
                if (!graph.hasBond(second, neighbor) || graph.bondOrder(first, neighbor) != graph.bondOrder(second, neighbor)
                    || graph.bondInRing(first, neighbor) != graph.bondInRing(second, neighbor)
                    || graph.bondAromatic(first, neighbor) != graph.bondAromatic(second, neighbor)) { equivalent = false; break; }
            }
            if (equivalent) classes[second] = classes[first];
        }
    return classes;
}

class VertexSearch {
    const QueryView& query_;
    const MolGraph& target_;
    const ChemOptions& chemistry_;
    const MCSOptions* common_;
    const Match* required_;
    bool induced_;
    bool all_;
    bool allowSkipping_;
    uint64_t probeNodeLimit_;
    Budget& budget_;
    Collector& collector_;
    RingConstraint queryRings_, targetRings_;
    std::vector<std::vector<uint64_t>> domains_;
    std::vector<int> domainCounts_, mapping_;
    std::vector<int> targetTwinClass_, queryTwinClass_;
    std::vector<bool> used_;
    std::vector<bool> excluded_;
    int availableTargets_;
    int mapped_ = 0;
    int targetEdges_ = 0;
    int uniformAdjacency_ = -1;
    bool targetBitAdjacencyAvailable_ = false;
    std::vector<uint8_t> availableQuery_, availableTarget_, reachableQuery_, reachableTarget_;
    std::vector<uint64_t> availableTargetBits_, unavailableTargetBits_;
    std::vector<int> mappedQueryAtoms_, componentQueue_;

    bool domainHas(int queryAtom, int targetAtom) const {
        return (domains_[queryAtom][targetAtom >> 6] & (uint64_t(1) << (targetAtom & 63))) != 0;
    }
    bool feasible(int queryAtom, int targetAtom) const {
        if (used_[targetAtom] || !domainHas(queryAtom, targetAtom)) return false;
        if (!queryTwinClass_.empty()) for (int previous = 0; previous < queryAtom; ++previous) {
            if (queryTwinClass_[previous] != queryTwinClass_[queryAtom]) continue;
            if (mapping_[previous] == -1 || (mapping_[previous] >= 0 && mapping_[previous] >= targetAtom)) return false;
        }
        for (int previous = 0; previous < query_.size(); ++previous) {
            if (mapping_[previous] < 0) continue;
            const bool queryBond = query_.adjacent(queryAtom, previous);
            const bool targetBond = target_.hasBond(targetAtom, mapping_[previous]);
            if (queryBond) {
                if (!targetBond || !query_.bondMatches(queryAtom, previous, target_, targetAtom, mapping_[previous], chemistry_))
                    return false;
            } else if (induced_ && targetBond) return false;
        }
        return true;
    }
    /** Current adjacency and bond predicates only remove prospective targets. */
    uint64_t feasibleTargets(int queryAtom, size_t word) const {
        uint64_t candidates = domains_[queryAtom][word] & ~unavailableTargetBits_[word];
        bool hasMappedNeighbor = false;
        for (int previous : mappedQueryAtoms_) {
            if (candidates == 0) break;
            const bool queryBond = query_.adjacent(queryAtom, previous);
            if (!queryBond && !induced_) continue;
            const int mappedTarget = mapping_[previous];
            uint64_t adjacent = 0;
            if (targetBitAdjacencyAvailable_) adjacent = target_.adjLong[mappedTarget][word];
            else for (int neighbor : target_.neighbors[mappedTarget])
                if (static_cast<size_t>(neighbor >> 6) == word) adjacent |= uint64_t(1) << (neighbor & 63);
            candidates &= queryBond ? adjacent : ~adjacent;
            hasMappedNeighbor = hasMappedNeighbor || queryBond;
        }
        if (!hasMappedNeighbor) return candidates;
        uint64_t remaining = candidates;
        while (remaining != 0) {
            const int offset = smsd::ctz64(remaining);
            const int candidate = static_cast<int>(word * 64) + offset;
            remaining &= remaining - 1;
            for (int neighbor : query_.neighbors(queryAtom))
                if (mapping_[neighbor] >= 0
                    && !query_.bondMatches(queryAtom, neighbor, target_, candidate, mapping_[neighbor], chemistry_)) {
                    candidates &= ~(uint64_t(1) << offset);
                    break;
                }
        }
        return candidates;
    }
    /** A connected extension remains in its mapped root's available component. */
    int connectedBound(bool queryGraph, const std::vector<uint8_t>& available,
                       std::vector<uint8_t>& reached, int root) {
        std::fill(reached.begin(), reached.end(), uint8_t(0));
        int largest = 0;
        const int count = queryGraph ? query_.size() : target_.n;
        for (int start = root < 0 ? 0 : root; start < count; ++start) {
            if (!available[start] || reached[start]) continue;
            int head = 0, tail = 0;
            componentQueue_[tail++] = start;
            reached[start] = 1;
            while (head < tail) {
                if ((head & 255) == 0 && budget_.stopped()) return -1;
                const int atom = componentQueue_[head++];
                const auto& adjacent = queryGraph ? query_.neighbors(atom) : target_.neighbors[atom];
                for (int neighbor : adjacent)
                    if (available[neighbor] && !reached[neighbor]) {
                        reached[neighbor] = 1;
                        componentQueue_[tail++] = neighbor;
                    }
            }
            if (root >= 0) return tail;
            largest = std::max(largest, tail);
        }
        return largest;
    }
    Score feasibleUpper(const Score& conservative) {
        std::fill(availableQuery_.begin(), availableQuery_.end(), uint8_t(0));
        std::fill(availableTargetBits_.begin(), availableTargetBits_.end(), uint64_t(0));
        std::fill(unavailableTargetBits_.begin(), unavailableTargetBits_.end(), uint64_t(0));
        mappedQueryAtoms_.clear();
        for (int atom = 0; atom < query_.size(); ++atom)
            if (mapping_[atom] >= 0) {
                availableQuery_[atom] = 1;
                mappedQueryAtoms_.push_back(atom);
                const int targetAtom = mapping_[atom];
                unavailableTargetBits_[targetAtom >> 6] |= uint64_t(1) << (targetAtom & 63);
            }
        int remaining = 0;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if ((atom & 15) == 0 && budget_.stopped()) return conservative;
            if (mapping_[atom] != -2) continue;
            bool hasCandidate = false;
            for (size_t word = 0; word < domains_[atom].size(); ++word) {
                if ((word & 255) == 0 && budget_.stopped()) return conservative;
                const uint64_t candidates = feasibleTargets(atom, word);
                availableTargetBits_[word] |= candidates;
                hasCandidate = hasCandidate || candidates != 0;
            }
            if (hasCandidate) { availableQuery_[atom] = 1; ++remaining; }
        }
        int targetCount = 0;
        for (int atom = 0; atom < target_.n; ++atom) {
            availableTarget_[atom] = used_[atom]
                || (availableTargetBits_[atom >> 6] & (uint64_t(1) << (atom & 63))) != 0;
            targetCount += availableTarget_[atom] != 0;
        }
        int atoms = std::min({conservative.atoms, mapped_ + remaining, targetCount});
        const bool connected = common_->connectivity == MCSConnectivity::Connected || common_->maxFragments == 1;
        if (connected) {
            const int queryRoot = mappedQueryAtoms_.empty() ? -1 : mappedQueryAtoms_.front();
            const int queryBound = connectedBound(true, availableQuery_, reachableQuery_, queryRoot);
            const int targetBound = connectedBound(false, availableTarget_, reachableTarget_,
                                                    queryRoot < 0 ? -1 : mapping_[queryRoot]);
            if (queryBound < 0 || targetBound < 0) return conservative;
            atoms = std::min({atoms, queryBound, targetBound});
        }
        int queryBonds = 0, targetBonds = 0;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if ((atom & 255) == 0 && budget_.stopped()) return conservative;
            if (!availableQuery_[atom] || (connected && !reachableQuery_[atom])) continue;
            for (int neighbor : query_.neighbors(atom))
                if (atom < neighbor && availableQuery_[neighbor] && (!connected || reachableQuery_[neighbor])) ++queryBonds;
        }
        for (int atom = 0; atom < target_.n; ++atom) {
            if ((atom & 255) == 0 && budget_.stopped()) return conservative;
            if (!availableTarget_[atom] || (connected && !reachableTarget_[atom])) continue;
            for (int neighbor : target_.neighbors[atom])
                if (atom < neighbor && availableTarget_[neighbor] && (!connected || reachableTarget_[neighbor])) ++targetBonds;
        }
        return {common_->atomWeights.empty() ? atoms : conservative.primary, atoms,
                static_cast<int>(std::min<int64_t>({conservative.bonds, queryBonds, targetBonds,
                                                   int64_t(atoms) * (atoms - 1) / 2}))};
    }
    Score upper() {
        int remaining = 0;
        double possibleWeight = 0.0;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if (mapping_[atom] == -2 && domainCounts_[atom] != 0) ++remaining;
            if (common_ && !common_->atomWeights.empty()) {
                if (mapping_[atom] >= 0) possibleWeight += common_->atomWeights[atom];
                else if (mapping_[atom] == -2 && domainCounts_[atom] != 0 && common_->atomWeights[atom] > 0.0)
                    possibleWeight += common_->atomWeights[atom];
            }
        }
        int atoms = mapped_ + std::min(remaining, availableTargets_ - mapped_);
        if (uniformAdjacency_ >= 0) {
            std::vector<int> candidates;
            bool stopped = false;
            for (int candidate = 0; candidate < target_.n; ++candidate) {
                if ((candidate & 255) == 0 && budget_.stopped()) { stopped = true; break; }
                if (used_[candidate] || excluded_[candidate]) continue;
                bool compatibleTopology = true;
                for (int previous = 0; previous < query_.size(); ++previous)
                    if (mapping_[previous] >= 0 && target_.hasBond(candidate, mapping_[previous]) != (uniformAdjacency_ != 0)) {
                        compatibleTopology = false; break;
                    }
                if (!compatibleTopology) continue;
                for (int queryAtom = 0; queryAtom < query_.size(); ++queryAtom)
                    if (mapping_[queryAtom] == -2 && domainHas(queryAtom, candidate)) {
                        candidates.push_back(candidate); break;
                    }
            }
            std::vector<std::vector<int>> colors;
            for (int candidate : candidates) {
                if (budget_.stopped()) { stopped = true; break; }
                size_t color = 0;
                for (; color < colors.size(); ++color) {
                    bool independent = true;
                    for (int previous : colors[color])
                        if (target_.hasBond(candidate, previous) == (uniformAdjacency_ != 0)) { independent = false; break; }
                    if (independent) break;
                }
                if (color == colors.size()) colors.emplace_back();
                colors[color].push_back(candidate);
            }
            if (!stopped) atoms = std::min(atoms, mapped_ + static_cast<int>(colors.size()));
        }
        int bonds = 0;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if (mapping_[atom] == -1) continue;
            for (int neighbor : query_.neighbors(atom))
                if (neighbor > atom && mapping_[neighbor] != -1) ++bonds;
        }
        if (common_ && common_->objective == MCSObjective::MCES) return {std::min(bonds, targetEdges_), atoms, 0};
        const Score conservative{common_ && !common_->atomWeights.empty() ? weightedScore(possibleWeight) : atoms,
            atoms, static_cast<int>(std::min<int64_t>({bonds, targetEdges_, int64_t(atoms) * (atoms - 1) / 2}))};
        if (!collector_.needs(conservative)) return conservative;
        return common_ && allowSkipping_ ? feasibleUpper(conservative) : conservative;
    }
    bool record() {
        if (!common_ && mapped_ != query_.size()) return false;
        if (common_ && mapped_ == 0) return false;
        if (required_) for (const auto& atom : required_->atomMapping)
            if (mapping_[atom.first] != atom.second) return false;
        Match match = makeMatch(mapping_, mappedBonds(query_, target_, mapping_, chemistry_));
        if (common_ && common_->objective == MCSObjective::MCES) {
            if (match.bondMapping.empty()) return false;
            std::vector<int> incident(mapping_.size(), -1);
            for (const auto& bond : match.bondMapping) { incident[bond[0]] = bond[2]; incident[bond[1]] = bond[3]; }
            match = makeMatch(incident, std::move(match.bondMapping));
        }
        if (common_) {
            if (common_->connectivity == MCSConnectivity::Connected && match.fragments.size() > 1) return false;
            if (match.fragments.size() > static_cast<size_t>(common_->maxFragments)) return false;
            for (const auto& fragment : match.fragments)
                if (fragment.size() < static_cast<size_t>(common_->minFragmentAtoms)) return false;
        }
        if (query_.molecule && !selectedStereoCompatible(*query_.molecule, target_, chemistry_, match)) return false;
        if (query_.pattern && !query_.pattern->stereoCompatible(target_, mapping_)) return false;
        if (chemistry_.completeRingsOnly) {
            if (query_.molecule && !completeRings(match, queryRings_, targetRings_)) return false;
            if (query_.pattern) {
                std::set<int> atoms;
                std::set<int64_t> bonds;
                for (const auto& atom : match.atomMapping) atoms.insert(atom.second);
                for (const auto& bond : match.bondMapping) bonds.insert(MolGraph::bondKey(bond[2], bond[3]));
                if (!targetRings_.accepts(atoms, bonds)) return false;
            }
        }
        return collector_.remember(std::move(match));
    }
    struct Frame {
        int atom;
        std::vector<int> candidates;
        size_t next = 0;
        int active = -2;
        bool skipped = false;
    };
    Frame frame() const {
        int chosen = -1, bestFrontier = -1;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if (mapping_[atom] != -2) continue;
            int frontier = 0;
            for (int neighbor : query_.neighbors(atom)) frontier += mapping_[neighbor] >= 0;
            if (chosen < 0 || (frontier != 0 && bestFrontier == 0)
                || ((frontier != 0) == (bestFrontier != 0)
                    && (domainCounts_[atom] < domainCounts_[chosen]
                        || (domainCounts_[atom] == domainCounts_[chosen]
                            && query_.neighbors(atom).size() > query_.neighbors(chosen).size())))) {
                chosen = atom;
                bestFrontier = frontier;
            }
        }
        Frame result{chosen, {}};
        if (chosen < 0) return result;
        if (common_ && (common_->connectivity == MCSConnectivity::Connected || common_->maxFragments == 1)
            && mapped_ > 0 && bestFrontier == 0)
            return Frame{-1, {}};
        int neighborAnchor = -1;
        for (int neighbor : query_.neighbors(chosen))
            if (mapping_[neighbor] >= 0 && (neighborAnchor < 0
                || target_.neighbors[mapping_[neighbor]].size() < target_.neighbors[mapping_[neighborAnchor]].size()))
                neighborAnchor = neighbor;
        if (neighborAnchor >= 0) {
            for (int candidate : target_.neighbors[mapping_[neighborAnchor]])
                if (feasible(chosen, candidate)) result.candidates.push_back(candidate);
        } else {
            for (size_t word = 0; word < domains_[chosen].size(); ++word) {
                uint64_t remaining = domains_[chosen][word];
                while (remaining != 0) {
                    const int candidate = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                    remaining &= remaining - 1;
                    if (feasible(chosen, candidate)) result.candidates.push_back(candidate);
                }
            }
        }
        if (query_.molecule) {
            std::vector<int> support(target_.n, 0);
            for (int candidate : result.candidates)
                for (int neighbor : query_.neighbors(chosen)) {
                    if (mapping_[neighbor] == -1) continue;
                    for (int adjacent : target_.neighbors[candidate])
                        if (!used_[adjacent] && domainHas(neighbor, adjacent)
                            && query_.bondMatches(chosen, neighbor, target_, candidate, adjacent, chemistry_)) {
                            ++support[candidate]; break;
                        }
                }
            std::stable_sort(result.candidates.begin(), result.candidates.end(), [&](int first, int second) {
                if (support[first] != support[second]) return support[first] > support[second];
                const auto queryDegree = static_cast<int>(query_.neighbors(chosen).size());
                return std::abs(target_.degree[first] - queryDegree) < std::abs(target_.degree[second] - queryDegree);
            });
        }
        if (!all_ && !required_ && !targetTwinClass_.empty()) {
            std::set<int> represented;
            result.candidates.erase(std::remove_if(result.candidates.begin(), result.candidates.end(), [&](int candidate) {
                return !represented.insert(targetTwinClass_[candidate]).second;
            }), result.candidates.end());
        }
        return result;
    }
public:
    VertexSearch(const QueryView& query, const MolGraph& target, const ChemOptions& chemistry,
                 const MCSOptions* common, bool induced, bool all, Budget& budget, Collector& collector, int anchor = -1,
                 const Match* required = nullptr, bool allowSkipping = true, uint64_t probeNodeLimit = 0)
        : query_(query), target_(target), chemistry_(chemistry), common_(common), required_(required), induced_(induced), all_(all),
          allowSkipping_(common && allowSkipping), probeNodeLimit_(probeNodeLimit),
          budget_(budget), collector_(collector),
          queryRings_(query.molecule ? *query.molecule : target, chemistry.completeRingsOnly && query.molecule),
          targetRings_(target, chemistry.completeRingsOnly),
          domains_(query.size(), std::vector<uint64_t>((target.n + 63) >> 6, 0)),
          domainCounts_(query.size(), 0), mapping_(query.size(), -2), used_(target.n, false),
          excluded_(target.n, false), availableTargets_(target.n) {
        if (common) for (int atom : common->excludedTargetAtoms)
            if (!excluded_[atom]) { excluded_[atom] = true; --availableTargets_; }
        if (common && common->objective == MCSObjective::MCIS) {
            if (allowSkipping_) {
                availableQuery_.resize(query.size()); availableTarget_.resize(target.n);
                reachableQuery_.resize(query.size()); reachableTarget_.resize(target.n);
                availableTargetBits_.resize((target.n + 63) >> 6);
                unavailableTargetBits_.resize(availableTargetBits_.size());
                mappedQueryAtoms_.reserve(query.size());
                componentQueue_.resize(std::max(query.size(), target.n));
                targetBitAdjacencyAvailable_ = target.adjLong.size() == static_cast<size_t>(target.n)
                    && std::all_of(target.adjLong.begin(), target.adjLong.end(), [&](const auto& row) {
                        return row.size() >= availableTargetBits_.size();
                    });
            }
            bool complete = true, independent = true;
            for (int atom = 0; atom < query.size(); ++atom) {
                complete = complete && query.neighbors(atom).size() == static_cast<size_t>(query.size() - 1);
                independent = independent && query.neighbors(atom).empty();
            }
            if (complete || independent) uniformAdjacency_ = complete ? 1 : 0;
        }
        if (query.molecule && !all && !required) {
            targetTwinClass_ = transpositionClasses(target, budget);
            if (common && common->objective == MCSObjective::MCIS)
                queryTwinClass_ = transpositionClasses(*query.molecule, budget, &common->atomWeights);
        }
        for (int atom = 0; atom < target.n; ++atom)
            for (int neighbor : target.neighbors[atom]) targetEdges_ += neighbor > atom;
        std::vector<uint64_t> acceleratorDomain;
        if (query.molecule && int64_t(query.size()) * target.n > 2000 && gpu_kern::graphKernelsAvailable()) {
            std::vector<int> queryRings(query.size()), queryAromatic(query.size()), targetRings(target.n), targetAromatic(target.n);
            for (int atom = 0; atom < query.size(); ++atom) {
                queryRings[atom] = query.molecule->ring[atom]; queryAromatic[atom] = query.molecule->aromatic[atom];
            }
            for (int atom = 0; atom < target.n; ++atom) {
                targetRings[atom] = target.ring[atom]; targetAromatic[atom] = target.aromatic[atom];
            }
            const int words = (target.n + 63) >> 6;
            acceleratorDomain.assign(static_cast<size_t>(query.size()) * words, 0);
            if (!gpu_kern::domainInit(query.size(), target.n, words,
                    query.molecule->atomicNum.data(), query.molecule->formalCharge.data(), queryRings.data(), queryAromatic.data(),
                    target.atomicNum.data(), target.formalCharge.data(), targetRings.data(), targetAromatic.data(),
                    chemistry.matchAtomType && !chemistry.tautomerAware, chemistry.matchFormalCharge,
                    chemistry.ringMatchesRingOnly, chemistry.aromaticityMode == ChemOptions::AromaticityMode::STRICT,
                    acceleratorDomain.data())) acceleratorDomain.clear();
        }
        for (int atom = 0; atom < query.size() && !budget.stopped(); ++atom) {
            for (int candidate = 0; candidate < target.n; ++candidate) {
                if ((candidate & 255) == 0 && budget.stopped()) break;
                if (excluded_[candidate]) continue;
                if (!acceleratorDomain.empty() && (acceleratorDomain[static_cast<size_t>(atom) * domains_[atom].size()
                        + (candidate >> 6)] & (uint64_t(1) << (candidate & 63))) == 0) continue;
                if (required && required->atomMapping.count(atom) && required->atomMapping.at(atom) != candidate) continue;
                if (atom == 0 && anchor >= 0 && candidate != anchor) continue;
                if (!allowSkipping_ && query.neighbors(atom).size() > target.neighbors[candidate].size()) continue;
                if (query.atomMatches(atom, target, candidate, chemistry)) {
                    domains_[atom][candidate >> 6] |= uint64_t(1) << (candidate & 63);
                    ++domainCounts_[atom];
                }
            }
        }
    }
    std::pair<bool, bool> run() {
        if (budget_.stopped()) return {false, false};
        if (query_.size() == 0) { if (!common_) record(); return {true, true}; }
        if (!common_ && query_.size() > target_.n) return {true, true};
        const Score absoluteUpper = upper();
        bool optimumProven = false;
        const uint64_t startingNodes = budget_.nodesVisited();
        std::vector<Frame> pending;
        pending.push_back(frame());
        while (!pending.empty() && !budget_.stopped()) {
            if (probeNodeLimit_ && budget_.nodesVisited() - startingNodes >= probeNodeLimit_) return {optimumProven, false};
            auto& current = pending.back();
            if (current.active >= 0) {
                used_[current.active] = false;
                --mapped_;
                mapping_[current.atom] = -2;
                current.active = -2;
            } else if (current.active == -1) {
                mapping_[current.atom] = -2;
                current.active = -2;
            }
            if (current.next < current.candidates.size()) {
                if (!budget_.visit()) break;
                const int candidate = current.candidates[current.next++];
                mapping_[current.atom] = candidate;
                used_[candidate] = true;
                ++mapped_;
                current.active = candidate;
                const bool accepted = record();
                if (!allowSkipping_ && accepted && !all_ && mapped_ == query_.size()) return {true, true};
                if (common_ && collector_.hasMatch() && (collector_.best() == absoluteUpper || collector_.optimumKnown())) {
                    optimumProven = true;
                    if (!all_) return {true, true};
                    if (collector_.truncated()) return {true, false};
                }
                if (!common_ && collector_.truncated()) return {true, false};
                if (common_ && !collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.atom >= 0) pending.push_back(std::move(next));
            } else if (allowSkipping_ && !current.skipped && (!required_ || !required_->atomMapping.count(current.atom))) {
                if (!budget_.visit()) break;
                current.skipped = true;
                current.active = -1;
                mapping_[current.atom] = -1;
                if (!collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.atom >= 0) pending.push_back(std::move(next));
            } else {
                pending.pop_back();
            }
        }
        return {pending.empty() || optimumProven || collector_.optimumKnown() || (!common_ && collector_.hasMatch()), pending.empty()};
    }
};

inline bool recursiveMatch(const SmartsQuery& pattern, const MolGraph& target, int anchor) {
    if (!recursiveBudget || recursiveBudget->stopped()) return false;
    const QueryView query(pattern);
    const ChemOptions chemistry;
    Collector collector(nullptr, 1, false);
    VertexSearch search(query, target, chemistry, nullptr, false, false, *recursiveBudget, collector, anchor);
    search.run();
    return collector.hasMatch();
}

class EdgeSearch {
    const MolGraph& query_;
    const MolGraph& target_;
    const MCSOptions& options_;
    const Match* required_;
    bool all_;
    Budget& budget_;
    Collector& collector_;
    RingConstraint queryRings_, targetRings_;
    std::vector<std::pair<int, int>> queryEdges_, targetEdges_;
    std::vector<int> edgeState_, mapping_, reverse_, references_;
    std::vector<bool> requiredEdges_;
    std::vector<bool> excluded_;
    std::vector<std::array<int, 4>> selected_;
    int mapped_ = 0;

    bool feasible(int edge, int firstTarget, int secondTarget) const {
        if (excluded_[firstTarget] || excluded_[secondTarget]) return false;
        const auto [first, second] = queryEdges_[edge];
        if (required_) {
            const auto requiredFirst = required_->atomMapping.find(first), requiredSecond = required_->atomMapping.find(second);
            if ((requiredFirst != required_->atomMapping.end() && requiredFirst->second != firstTarget)
                || (requiredSecond != required_->atomMapping.end() && requiredSecond->second != secondTarget)) return false;
        }
        if ((mapping_[first] >= 0 && mapping_[first] != firstTarget)
            || (mapping_[second] >= 0 && mapping_[second] != secondTarget)) return false;
        if ((reverse_[firstTarget] >= 0 && reverse_[firstTarget] != first)
            || (reverse_[secondTarget] >= 0 && reverse_[secondTarget] != second)) return false;
        return atomCompatible(query_, first, target_, firstTarget, options_.chemistry)
            && atomCompatible(query_, second, target_, secondTarget, options_.chemistry)
            && bondCompatible(query_, first, second, target_, firstTarget, secondTarget, options_.chemistry);
    }
    Score upper() const {
        int remaining = 0;
        for (int state : edgeState_) remaining += state == 0;
        const int selected = static_cast<int>(selected_.size());
        const int bonds = selected + std::min(remaining, static_cast<int>(targetEdges_.size()) - selected);
        return {bonds, std::min({query_.n, target_.n, mapped_ + 2 * (bonds - selected)}), 0};
    }
    bool record() {
        if (selected_.empty()) return false;
        if (required_) {
            for (const auto& atom : required_->atomMapping) if (mapping_[atom.first] != atom.second) return false;
            for (size_t edge = 0; edge < requiredEdges_.size(); ++edge)
                if (requiredEdges_[edge] && edgeState_[edge] != 1) return false;
        }
        Match match = makeMatch(mapping_, selected_);
        if (options_.connectivity == MCSConnectivity::Connected && match.fragments.size() > 1) return false;
        if (match.fragments.size() > static_cast<size_t>(options_.maxFragments)) return false;
        for (const auto& fragment : match.fragments)
            if (fragment.size() < static_cast<size_t>(options_.minFragmentAtoms)) return false;
        if (!selectedStereoCompatible(query_, target_, options_.chemistry, match)) return false;
        if (options_.chemistry.completeRingsOnly && !completeRings(match, queryRings_, targetRings_)) return false;
        return collector_.remember(std::move(match));
    }
    struct Frame {
        int edge;
        size_t next = 0;
        std::vector<std::array<int, 2>> anchored;
        bool hasAnchor = false;
        bool active = false;
        bool skipped = false;
    };
    Frame frame() const {
        int chosen = -1, bestFrontier = -1;
        for (size_t edge = 0; edge < queryEdges_.size(); ++edge) {
            if (edgeState_[edge] != 0) continue;
            const auto [first, second] = queryEdges_[edge];
            const int frontier = (mapping_[first] >= 0) + (mapping_[second] >= 0);
            if (chosen < 0 || frontier > bestFrontier) { chosen = static_cast<int>(edge); bestFrontier = frontier; }
        }
        Frame result{chosen};
        if (chosen < 0) return result;
        if ((options_.connectivity == MCSConnectivity::Connected || options_.maxFragments == 1)
            && !selected_.empty() && bestFrontier == 0)
            return Frame{-1};
        const auto [first, second] = queryEdges_[chosen];
        if (mapping_[first] >= 0) {
            result.hasAnchor = true;
            for (int neighbor : target_.neighbors[mapping_[first]])
                result.anchored.push_back({mapping_[first], neighbor});
        } else if (mapping_[second] >= 0) {
            result.hasAnchor = true;
            for (int neighbor : target_.neighbors[mapping_[second]])
                result.anchored.push_back({neighbor, mapping_[second]});
        }
        return result;
    }
    bool nextCandidate(Frame& current, std::array<int, 2>& candidate) {
        const size_t count = current.hasAnchor ? current.anchored.size() : targetEdges_.size() * 2;
        while (current.next < count) {
            if ((current.next & 255) == 0 && budget_.stopped()) return false;
            const size_t position = current.next++;
            if (current.hasAnchor) candidate = current.anchored[position];
            else {
                const auto edge = targetEdges_[position / 2];
                candidate = (position & 1) == 0 ? std::array<int, 2>{edge.first, edge.second}
                    : std::array<int, 2>{edge.second, edge.first};
            }
            if (feasible(current.edge, candidate[0], candidate[1])) return true;
        }
        return false;
    }
    void assign(int queryAtom, int targetAtom) {
        if (references_[queryAtom]++ == 0) {
            mapping_[queryAtom] = targetAtom;
            reverse_[targetAtom] = queryAtom;
            ++mapped_;
        }
    }
    void unassign(int queryAtom) {
        if (--references_[queryAtom] == 0) {
            reverse_[mapping_[queryAtom]] = -1;
            mapping_[queryAtom] = -1;
            --mapped_;
        }
    }
public:
    EdgeSearch(const MolGraph& query, const MolGraph& target, const MCSOptions& options,
               bool all, Budget& budget, Collector& collector, const Match* required = nullptr)
        : query_(query), target_(target), options_(options), required_(required), all_(all), budget_(budget), collector_(collector),
          queryRings_(query, options.chemistry.completeRingsOnly), targetRings_(target, options.chemistry.completeRingsOnly),
          mapping_(query.n, -1), reverse_(target.n, -1), references_(query.n, 0), excluded_(target.n, false) {
        for (int atom : options.excludedTargetAtoms) excluded_[atom] = true;
        for (int atom = 0; atom < query.n; ++atom)
            for (int neighbor : query.neighbors[atom]) if (atom < neighbor) queryEdges_.emplace_back(atom, neighbor);
        for (int atom = 0; atom < target.n; ++atom)
            for (int neighbor : target.neighbors[atom])
                if (atom < neighbor && !excluded_[atom] && !excluded_[neighbor]) targetEdges_.emplace_back(atom, neighbor);
        edgeState_.assign(queryEdges_.size(), 0);
        requiredEdges_.assign(queryEdges_.size(), false);
        if (required) for (size_t edge = 0; edge < queryEdges_.size(); ++edge)
            for (const auto& bond : required->bondMapping)
                if (MolGraph::bondKey(queryEdges_[edge].first, queryEdges_[edge].second) == MolGraph::bondKey(bond[0], bond[1]))
                    requiredEdges_[edge] = true;
    }
    std::pair<bool, bool> run() {
        if (budget_.stopped()) return {false, false};
        if (queryEdges_.empty() || targetEdges_.empty()) return {true, true};
        Score absoluteUpper = upper();
        int queryIncident = 0, targetIncident = 0;
        for (const auto& adjacent : query_.neighbors) queryIncident += !adjacent.empty();
        for (const auto& adjacent : target_.neighbors) targetIncident += !adjacent.empty();
        absoluteUpper.atoms = std::min({absoluteUpper.atoms, queryIncident, targetIncident});
        bool optimumProven = false;
        std::vector<Frame> pending;
        pending.push_back(frame());
        while (!pending.empty() && !budget_.stopped()) {
            auto& current = pending.back();
            if (current.active) {
                const auto [first, second] = queryEdges_[current.edge];
                unassign(first);
                unassign(second);
                selected_.pop_back();
                edgeState_[current.edge] = 0;
                current.active = false;
            } else if (current.skipped && edgeState_[current.edge] == -1) {
                edgeState_[current.edge] = 0;
            }
            std::array<int, 2> candidate{};
            if (nextCandidate(current, candidate)) {
                if (!budget_.visit()) break;
                const auto [first, second] = queryEdges_[current.edge];
                assign(first, candidate[0]);
                assign(second, candidate[1]);
                selected_.push_back({first, second, candidate[0], candidate[1]});
                edgeState_[current.edge] = 1;
                current.active = true;
                record();
                if (collector_.hasMatch() && (collector_.best() == absoluteUpper || collector_.optimumKnown())) {
                    optimumProven = true;
                    if (!all_) return {true, true};
                    if (collector_.truncated()) return {true, false};
                }
                if (!collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.edge >= 0) pending.push_back(std::move(next));
            } else if (!current.skipped && !requiredEdges_[current.edge] && !budget_.stopped()) {
                if (!budget_.visit()) break;
                current.skipped = true;
                edgeState_[current.edge] = -1;
                if (!collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.edge >= 0) pending.push_back(std::move(next));
            } else pending.pop_back();
        }
        return {pending.empty() || optimumProven || collector_.optimumKnown(), pending.empty()};
    }
};

inline void prepare(const MolGraph* query, const MolGraph& target, const ChemOptions& chemistry) {
    if (chemistry.ringFusionMode != ChemOptions::RingFusionMode::IGNORE) {
        if (query) query->ensureRingCounts();
        target.ensureRingCounts();
    }
}

inline SearchResult commonSubgraph(const MolGraph& query, const MolGraph& target,
                                  const MCSOptions& options, size_t limit, bool all, const Match* required = nullptr,
                                  Budget* sharedBudget = nullptr) {
    validateOptions(options, query.n);
    for (int atom : options.excludedTargetAtoms)
        if (atom < 0 || atom >= target.n) throw std::out_of_range("Excluded target atom index out of range");
    Budget localBudget(options, query.n, target.n, true);
    Budget& budget = sharedBudget ? *sharedBudget : localBudget;
    RecursiveScope recursiveScope(budget);
    Collector collector(&options, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    if (all) {
        const auto optimum = commonSubgraph(query, target, options, 1, false, required, &budget);
        for (const auto& match : optimum.matches) collector.remember(match);
        if (optimum.optimal) collector.markOptimumKnown();
        if (budget.stopped()) return collector.finish(budget, optimum.optimal, false);
    }
    prepare(&query, target, options.chemistry);
    std::pair<bool, bool> proof;
    if (!all && !required && query.n > 0 && query.n <= target.n) {
        const QueryView view(query);
        VertexSearch full(view, target, options.chemistry, &options, options.objective == MCSObjective::MCIS,
                          false, budget, collector, -1, nullptr, false, 2048);
        full.run();
        if (!all && options.objective == MCSObjective::MCIS && collector.hasMatch()
            && collector.best().atoms == query.n
            && std::all_of(options.atomWeights.begin(), options.atomWeights.end(), [](double weight) { return weight >= 0.0; }))
            return collector.finish(budget, true, true);
    }
    if (options.objective == MCSObjective::MCES) {
        if (!all && !required && !budget.stopped()) {
            const QueryView view(query);
            VertexSearch seed(view, target, options.chemistry, &options, false, false, budget, collector,
                              -1, nullptr, true, 2048);
            seed.run();
        }
        EdgeSearch search(query, target, options, all, budget, collector, required);
        proof = search.run();
    } else {
        const QueryView view(query);
        VertexSearch search(view, target, options.chemistry, &options, true, all, budget, collector, -1, required);
        proof = search.run();
    }
    return collector.finish(budget, proof.first, proof.second);
}

inline SearchResult substructure(const MolGraph* query, const SmartsQuery* pattern,
                                const MolGraph& target, const SubstructureOptions& options,
                                size_t limit, bool all, Budget* sharedBudget = nullptr) {
    const QueryView view = query ? QueryView(*query) : QueryView(*pattern);
    Budget localBudget(options, view.size(), target.n, false);
    Budget& budget = sharedBudget ? *sharedBudget : localBudget;
    RecursiveScope recursiveScope(budget);
    Collector collector(nullptr, limit, all);
    if (budget.stopped()) return collector.finish(budget, false, false);
    prepare(query, target, options.chemistry);
    VertexSearch search(view, target, options.chemistry, nullptr, options.induced, all, budget, collector);
    const auto proof = search.run();
    return collector.finish(budget, proof.first, proof.second);
}

} // namespace search_detail
inline MolGraph Match::toGraph(const MolGraph& query) const {
    const std::lock_guard<std::recursive_mutex> lock(query.lazyCacheMutex());
    std::vector<int> originals;
    std::vector<int> reindex(query.n, -1);
    for (const auto& atom : atomMapping) {
        if (atom.first < 0 || atom.first >= query.n)
            throw std::out_of_range("Selected atom index out of range");
        reindex[atom.first] = static_cast<int>(originals.size());
        originals.push_back(atom.first);
    }
    std::set<int64_t> selected;
    for (const auto& bond : bondMapping) {
        if (bond[0] < 0 || bond[1] < 0 || bond[0] >= query.n || bond[1] >= query.n)
            throw std::out_of_range("Selected bond atom index out of range");
        if (reindex[bond[0]] < 0 || reindex[bond[1]] < 0 || !query.hasBond(bond[0], bond[1])
            || atomMapping.at(bond[0]) != bond[2] || atomMapping.at(bond[1]) != bond[3])
            throw std::invalid_argument("Selected bond is not in the atom mapping");
        if (!selected.insert(MolGraph::bondKey(bond[0], bond[1])).second)
            throw std::invalid_argument("Selected bond occurs more than once");
    }
    const int count = static_cast<int>(originals.size());
    std::vector<int> numbers(count), charges(count), masses(count), hydrogens(count), classes(count), chirality(count), ids;
    std::vector<uint8_t> rings(count), aromatic(count);
    std::vector<std::vector<int>> adjacent(count), orders(count);
    std::vector<std::vector<bool>> bondRings(count), bondAromatic(count);
    std::vector<bool> allLigands(count, true);
    if (!query.atomId.empty()) ids.resize(count);
    auto property = [](const auto& values, int atom, int fallback = 0) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : fallback;
    };
    for (int atom = 0; atom < count; ++atom) {
        const int original = originals[atom];
        numbers[atom] = query.atomicNum[original]; charges[atom] = property(query.formalCharge, original);
        masses[atom] = property(query.massNumber, original); hydrogens[atom] = property(query.hydrogenCount, original);
        classes[atom] = property(query.atomClass, original); rings[atom] = query.ring[original]; aromatic[atom] = query.aromatic[original];
        if (!ids.empty()) ids[atom] = query.atomId[original];
        for (int neighbor : query.neighbors[original]) {
            if (reindex[neighbor] < 0 || selected.count(MolGraph::bondKey(original, neighbor)) == 0) {
                allLigands[atom] = false;
                continue;
            }
            adjacent[atom].push_back(reindex[neighbor]); orders[atom].push_back(query.bondOrder(original, neighbor));
            bondRings[atom].push_back(query.bondInRing(original, neighbor));
            bondAromatic[atom].push_back(query.bondAromatic(original, neighbor));
        }
        if (allLigands[atom]) chirality[atom] = property(query.tetraChirality, original);
    }
    std::vector<std::tuple<int, int, int>> stereo;
    for (const auto& bond : bondMapping) {
        const int first = reindex[bond[0]], second = reindex[bond[1]];
        if (!allLigands[first] || !allLigands[second]) continue;
        const int configuration = query.dbStereo(bond[0], bond[1]);
        if (configuration != 0) stereo.emplace_back(first, second, configuration);
    }
    auto builder = MolGraph::Builder().atomCount(count).atomicNumbers(std::move(numbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses)).hydrogenCounts(std::move(hydrogens))
        .atomClasses(std::move(classes)).ringFlags(std::move(rings)).aromaticFlags(std::move(aromatic))
        .tetrahedralChirality(std::move(chirality)).setNeighbors(std::move(adjacent)).setBondOrders(std::move(orders))
        .bondRingFlags(std::move(bondRings)).bondAromaticFlags(std::move(bondAromatic))
        .doubleBondStereoSparse(std::move(stereo)).solvent(query.solvent_).pH(query.pH_);
    if (!ids.empty()) builder.atomIds(std::move(ids));
    MolGraph graph = builder.build(false);
    graph.hydrogenCountsKnown = query.hydrogenCountsKnown;
    query.ensureRingCounts();
    for (int atom = 0; atom < count; ++atom) {
        const int original = originals[atom];
        graph.ringCount[atom] = query.ringCount[original];
        graph.tautomerClass[atom] = property(query.tautomerClass, original, -1);
        graph.tautomerWeight[atom] = original < static_cast<int>(query.tautomerWeight.size()) ? query.tautomerWeight[original] : 1.0f;
    }
    graph.ringCountsComputed_ = true;
    return graph;
}

inline std::string Match::toSmiles(const MolGraph& g) const {
    const auto& mapping = atomMapping;
    for (const auto& atom : mapping)
        if (atom.first < 0 || atom.first >= g.n) throw std::out_of_range("Selected atom index out of range");
    std::set<int64_t> selectedBonds;
    for (const auto& bond : bondMapping) {
        if (bond[0] < 0 || bond[1] < 0 || bond[0] >= g.n || bond[1] >= g.n)
            throw std::out_of_range("Selected bond atom index out of range");
        if (mapping.count(bond[0]) == 0 || mapping.count(bond[1]) == 0
            || !g.hasBond(bond[0], bond[1]) || mapping.at(bond[0]) != bond[2] || mapping.at(bond[1]) != bond[3])
            throw std::invalid_argument("Selected bond is not in the query mapping");
        if (!selectedBonds.insert(MolGraph::bondKey(bond[0], bond[1])).second)
            throw std::invalid_argument("Selected bond occurs more than once");
    }
    if (mapping.empty()) return "";
    for (const auto& pair : mapping)
        if (pair.first < 0 || pair.first >= g.n)
            throw std::out_of_range("MCS atom index out of range");

    std::vector<int> atoms;
    atoms.reserve(mapping.size());
    for (const auto& pair : mapping) atoms.push_back(pair.first);
    const int atomCount = static_cast<int>(atoms.size());
    std::unordered_map<int, int> reindex;
    for (int atom = 0; atom < atomCount; ++atom) reindex[atoms[atom]] = atom;

    std::vector<int> atomicNumbers(atomCount), charges(atomCount), masses(atomCount), hydrogens(atomCount);
    std::vector<int> chirality(atomCount, 0), classes(atomCount, 0);
    std::vector<uint8_t> aromatic(atomCount, 0), ring(atomCount, 0);
    std::vector<std::vector<int>> neighbors(atomCount), orders(atomCount);
    std::vector<std::vector<bool>> bondRings(atomCount), bondAromatic(atomCount);
    std::vector<bool> allLigandsRetained(atomCount, true);
    auto property = [](const auto& values, int atom) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : 0;
    };
    for (int atom = 0; atom < atomCount; ++atom) {
        const int original = atoms[atom];
        atomicNumbers[atom] = g.atomicNum[original];
        charges[atom] = property(g.formalCharge, original);
        classes[atom] = property(g.atomClass, original);
        masses[atom] = property(g.massNumber, original);
        hydrogens[atom] = detail::hydrogenCountForWrite(g, original);
        for (int neighbor : g.neighbors[original]) {
            const auto retained = reindex.find(neighbor);
            if (retained == reindex.end() || selectedBonds.count(MolGraph::bondKey(original, neighbor)) == 0) {
                allLigandsRetained[atom] = false;
                continue;
            }
            neighbors[atom].push_back(retained->second);
            orders[atom].push_back(g.bondOrder(original, neighbor));
        }
        if (allLigandsRetained[atom]) chirality[atom] = property(g.tetraChirality, original);
    }

    std::vector<bool> retainedRingAtoms;
    std::vector<std::set<int>> retainedRingBonds;
    detail::detectRings(atomCount, neighbors, retainedRingAtoms, retainedRingBonds);
    for (int atom = 0; atom < atomCount; ++atom) {
        ring[atom] = retainedRingAtoms[atom];
        for (int bondIndex = 0; bondIndex < static_cast<int>(neighbors[atom].size()); ++bondIndex) {
            const int neighbor = neighbors[atom][bondIndex];
            const bool inRing = retainedRingBonds[atom].count(neighbor) != 0;
            const bool isAromatic = inRing && g.bondAromatic(atoms[atom], atoms[neighbor]);
            bondRings[atom].push_back(inRing);
            bondAromatic[atom].push_back(isAromatic);
            if (isAromatic) aromatic[atom] = 1;
            else if (orders[atom][bondIndex] == 4) orders[atom][bondIndex] = 1;
        }
        if (!allLigandsRetained[atom] || aromatic[atom] != g.aromatic[atoms[atom]]) {
            int bondOrderSum = 0;
            for (int order : orders[atom]) bondOrderSum += order == 4 ? 1 : order;
            if (aromatic[atom] && (atomicNumbers[atom] == 7 || atomicNumbers[atom] == 15))
                hydrogens[atom] += static_cast<int>(g.neighbors[atoms[atom]].size() - neighbors[atom].size());
            else
                hydrogens[atom] = atomicNumbers[atom] <= 1 ? 0
                    : detail::computeImplicitH(atomicNumbers[atom], aromatic[atom], bondOrderSum, charges[atom]);
        }
    }

    std::vector<std::tuple<int, int, int>> doubleBondStereo;
    for (int atom = 0; atom < atomCount; ++atom) {
        if (!allLigandsRetained[atom]) continue;
        for (int neighbor : neighbors[atom]) {
            if (atom >= neighbor || !allLigandsRetained[neighbor]) continue;
            const int configuration = g.dbStereo(atoms[atom], atoms[neighbor]);
            if (configuration == 0) continue;
            doubleBondStereo.emplace_back(atom, neighbor, configuration);
        }
    }
    MolGraph sub = MolGraph::Builder()
        .atomCount(atomCount).atomicNumbers(std::move(atomicNumbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses)).atomClasses(std::move(classes))
        .hydrogenCounts(std::move(hydrogens)).tetrahedralChirality(std::move(chirality))
        .ringFlags(std::move(ring)).aromaticFlags(std::move(aromatic))
        .setNeighbors(std::move(neighbors)).setBondOrders(std::move(orders))
        .bondRingFlags(std::move(bondRings)).bondAromaticFlags(std::move(bondAromatic))
        .doubleBondStereoSparse(std::move(doubleBondStereo)).build();
    if (!sub.kekulize()) {
        for (int atom = 0; atom < sub.n; ++atom) {
            if (!sub.aromatic[atom]) continue;
            sub.aromatic[atom] = 0;
            int bondOrderSum = 0;
            for (int neighbor : sub.neighbors[atom]) {
                if (sub.bondOrder(atom, neighbor) == 4) sub.setBondOrder(atom, neighbor, 1);
                sub.setBondAromatic(atom, neighbor, false);
                bondOrderSum += sub.bondOrder(atom, neighbor);
            }
            sub.hydrogenCount[atom] = detail::computeImplicitH(
                sub.atomicNum[atom], false, bondOrderSum, sub.formalCharge[atom]);
        }
    }
    sub.perceiveAromaticity();

    return toSMILES(sub);
}

namespace search_detail {

inline std::vector<std::vector<int>> componentAtoms(const MolGraph& graph) {
    std::vector<bool> visited(graph.n, false);
    std::vector<std::vector<int>> components;
    for (int atom = 0; atom < graph.n; ++atom) {
        if (visited[atom]) continue;
        visited[atom] = true;
        std::vector<int> component{atom};
        for (size_t head = 0; head < component.size(); ++head)
            for (int neighbor : graph.neighbors[component[head]])
                if (!visited[neighbor]) { visited[neighbor] = true; component.push_back(neighbor); }
        std::sort(component.begin(), component.end());
        components.push_back(std::move(component));
    }
    return components;
}

inline int countComponents(const MolGraph& graph) {
    return static_cast<int>(componentAtoms(graph).size());
}

inline std::vector<MolGraph> splitComponents(const MolGraph& graph) {
    std::vector<MolGraph> result;
    for (const auto& component : componentAtoms(graph)) {
        Match selected;
        for (int atom : component) selected.atomMapping.emplace(atom, atom);
        for (int atom : component) for (int neighbor : graph.neighbors[atom])
            if (atom < neighbor) selected.bondMapping.push_back({atom, neighbor, atom, neighbor});
        result.push_back(selected.toGraph(graph));
    }
    return result;
}

inline bool sameCanonicalGraph(const MolGraph& first, const MolGraph& second) {
    return first.n == second.n && first.getCanonicalHash() == second.getCanonicalHash()
        && first.toCanonicalSmiles() == second.toCanonicalSmiles();
}

} // namespace search_detail

inline int mobilisableProtons(const MolGraph& g, int atom) {
    int z = g.atomicNum[atom];
    // Only N, O, S can carry mobilisable protons in tautomeric systems
    if (z != 7 && z != 8 && z != 16) return 0;
    // Sum bond orders to heavy-atom neighbours.
    int bondSum = 0;
    for (int nb : g.neighbors[atom])
        bondSum += g.bondOrder(atom, nb);
    // Default valence for tautomeric heteroatoms
    int defaultVal = (z == 7) ? 3 : 2; // N=3, O=2, S=2
    const int64_t charge = g.formalCharge[atom];
    const int64_t implicitH = int64_t(defaultVal) - bondSum - (charge < 0 ? -charge : charge);
    return static_cast<int>(std::max(int64_t(0), implicitH));
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

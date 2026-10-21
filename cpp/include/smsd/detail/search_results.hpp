/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * What searches accept and report: option checks and scores, matches built from
 * mappings, ring and stereo checks of selected matches, the result collector, and
 * the atom transpositions used while choosing one optimal witness.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/search_types.hpp"
#include "smsd/detail/search_budget.hpp"
#include "smsd/detail/search_compatibility.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace smsd {
namespace search_detail {

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
        if (mapping[atom] >= 0) result.atomMapping.emplace_hint(result.atomMapping.end(), static_cast<int>(atom), mapping[atom]);
    std::sort(bonds.begin(), bonds.end());
    result.bondMapping = std::move(bonds);
    // One buffer holds both the union-find parents and the fragment numbers.
    std::vector<int> scratch(2 * mapping.size(), -1);
    int* const parent = scratch.data();
    int* const fragmentIndex = scratch.data() + mapping.size();
    auto findRoot = [&](int atom) {
        int representative = atom;
        while (parent[representative] >= 0) representative = parent[representative];
        while (atom != representative) {
            const int next = parent[atom];
            parent[atom] = representative;
            atom = next;
        }
        return representative;
    };
    for (const auto& bond : result.bondMapping) {
        int firstRoot = findRoot(bond[0]), secondRoot = findRoot(bond[1]);
        if (firstRoot == secondRoot) continue;
        if (parent[firstRoot] > parent[secondRoot]) std::swap(firstRoot, secondRoot);
        parent[firstRoot] += parent[secondRoot];
        parent[secondRoot] = firstRoot;
    }
    // Number fragments by first appearance, then fill each at its final size.
    int fragments = 0;
    for (const auto& atom : result.atomMapping) {
        int& index = fragmentIndex[findRoot(atom.first)];
        if (index < 0) index = fragments++;
    }
    result.fragments.resize(fragments);
    for (const auto& atom : result.atomMapping) {
        const int representative = findRoot(atom.first);
        auto& fragment = result.fragments[fragmentIndex[representative]];
        if (fragment.empty()) fragment.reserve(static_cast<size_t>(-parent[representative]));
        fragment.push_back(atom.first);
    }
    return result;
}

inline std::vector<std::array<int, 4>> mappedBonds(const QueryView& query, const MolGraph& target,
                                                const std::vector<int>& mapping,
                                                const ChemOptions& chemistry) {
    std::vector<std::array<int, 4>> bonds;
    size_t possible = 0;
    for (int atom = 0; atom < query.size(); ++atom)
        if (mapping[atom] >= 0)
            for (int neighbor : query.neighbors(atom)) possible += atom < neighbor && mapping[neighbor] >= 0;
    bonds.reserve(possible);
    for (int atom = 0; atom < query.size(); ++atom) {
        if (mapping[atom] < 0) continue;
        for (int neighbor : query.neighbors(atom))
            if (atom < neighbor && mapping[neighbor] >= 0
                && query.bondMatches(atom, neighbor, target, mapping[atom], mapping[neighbor], chemistry))
                bonds.push_back({atom, neighbor, mapping[atom], mapping[neighbor]});
    }
    return bonds;
}

/**
 * The bonded atom pairs of a query, lower atom first and in ascending order, and its fragments:
 * the atoms those pairs connect, numbered by lowest atom and each in ascending order. These are
 * the bonds and fragments makeMatch() gives a complete mapping that selects every query bond.
 */
struct CompleteMappingShape {
    std::vector<std::array<int, 2>> pairs;
    std::vector<std::vector<int>> fragments;
    CompleteMappingShape(int atoms, const std::vector<std::vector<int>>& neighbors) {
        size_t listed = 0;
        for (int atom = 0; atom < atoms; ++atom) listed += neighbors[atom].size();
        pairs.reserve(listed / 2);
        for (int atom = 0; atom < atoms; ++atom)
            for (int neighbor : neighbors[atom]) if (atom < neighbor) pairs.push_back({atom, neighbor});
        std::sort(pairs.begin(), pairs.end());
        std::vector<int> component(static_cast<size_t>(atoms));
        std::iota(component.begin(), component.end(), 0);
        auto root = [&](int atom) {
            while (component[atom] != atom) atom = component[atom] = component[component[atom]];
            return atom;
        };
        for (const auto& pair : pairs) {
            const int first = root(pair[0]), second = root(pair[1]);
            if (first != second) component[std::max(first, second)] = std::min(first, second);
        }
        std::vector<int> fragmentOf(static_cast<size_t>(atoms), -1);
        for (int atom = 0; atom < atoms; ++atom) {
            int& fragment = fragmentOf[root(atom)];
            if (fragment < 0) {
                fragment = static_cast<int>(fragments.size());
                fragments.emplace_back();
            }
            fragments[fragment].push_back(atom);
        }
    }
};

/**
 * The complete mappings a substructure search accepts, kept flat in the order found, and the
 * matches they give. A kept mapping must select every query bond, so that the bonds and
 * fragments of its match are those of the query's CompleteMappingShape.
 */
class CompleteMappingRows {
    std::vector<int>& values_;
    size_t width_, limit_;
    bool all_;
    size_t count_ = 0;
public:
    /** Whether, with every mapping wanted, more mappings than the limit were found. */
    bool truncated = false;
    /** Rows of `width` atoms kept in `values`, which is cleared; with `all`, at most `limit` rows. */
    CompleteMappingRows(std::vector<int>& values, size_t width, size_t limit, bool all)
        : values_(values), width_(width), limit_(limit), all_(all) { values_.clear(); }
    CompleteMappingRows(const CompleteMappingRows&) = delete;
    CompleteMappingRows& operator=(const CompleteMappingRows&) = delete;
    size_t count() const { return count_; }
    const int* row(size_t index) const { return values_.data() + index * width_; }
    /** Keeps a complete mapping or, once every mapping is wanted and the limit is kept, notes the truncation. */
    void append(const std::vector<int>& mapping) {
        if (all_ && count_ >= limit_) {
            truncated = true;
            return;
        }
        values_.insert(values_.end(), mapping.begin(), mapping.begin() + static_cast<std::ptrdiff_t>(width_));
        ++count_;
    }
    /**
     * The matches in Collector order: complete mappings share their query atoms, so map order is
     * row order. With a budget, building counts towards it: the clock is read after about 4096
     * atoms and bonds, and past the results deadline the matches built so far are returned.
     */
    std::vector<Match> build(const CompleteMappingShape& shape, Budget* budget = nullptr) const {
        std::vector<size_t> order(count_);
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&](size_t first, size_t second) {
            return std::lexicographical_compare(row(first), row(first) + width_, row(second), row(second) + width_);
        });
        const auto& pairs = shape.pairs;
        std::vector<Match> result;
        result.reserve(count_);
        size_t work = 0;
        for (size_t index : order) {
            if ((work += 1 + width_ + pairs.size()) >= 4096) {
                work = 0;
                if (budget && budget->resultsStopped()) break;
            }
            const int* source = row(index);
            Match match;
            for (size_t atom = 0; atom < width_; ++atom)
                match.atomMapping.emplace_hint(match.atomMapping.end(), static_cast<int>(atom), source[atom]);
            match.bondMapping.reserve(pairs.size());
            for (const auto& pair : pairs)
                match.bondMapping.push_back({pair[0], pair[1], source[pair[0]], source[pair[1]]});
            match.fragments = shape.fragments;
            result.push_back(std::move(match));
        }
        return result;
    }
};

using CycleList = std::vector<std::vector<int>>;

/**
 * Relevant cycles of the graphs this thread searched most recently. In a batch the shared query
 * is met once per target and each target by several searches, so the cycles are found once per
 * thread instead of each time, and no state is shared between threads.
 *
 * An entry is keyed by the graph's address, cache nonce and atom count, so it relies on every
 * change to a graph's bonds renewing the nonce through MolGraph::invalidateBondDerivedState(),
 * as the MolGraph mutators do. A new graph at a reused address has a new nonce, and a copy
 * shares its source's nonce only with the same bonds. As a cheap check of that invariant the
 * entry also keeps the number of bond ends: a graph whose bonds changed without a new nonce
 * fails an assertion, and without assertions its cycles are found again.
 */
inline std::shared_ptr<const CycleList> recentRelevantCycles(const MolGraph& graph) {
    struct Entry {
        const MolGraph* graph = nullptr;
        uint64_t nonce = 0;
        int atoms = -1;
        size_t bondEnds = 0;
        std::shared_ptr<const CycleList> cycles;
    };
    thread_local std::array<Entry, 4> recent;
    size_t bondEnds = 0;
    for (const auto& adjacent : graph.neighbors) bondEnds += adjacent.size();
    for (size_t slot = 0; slot < recent.size(); ++slot) {
        if (recent[slot].graph != &graph || recent[slot].nonce != graph.cacheNonce || recent[slot].atoms != graph.n) continue;
        assert(recent[slot].bondEnds == bondEnds && "graph bonds changed without invalidateBondDerivedState()");
        if (recent[slot].bondEnds != bondEnds) break;
        // Most recently used first, so the query outlives the targets passing through.
        std::rotate(recent.begin(), recent.begin() + slot, recent.begin() + slot + 1);
        return recent.front().cycles;
    }
    auto cycles = std::make_shared<const CycleList>(graph.computeRelevantCycles());
    std::rotate(recent.begin(), recent.end() - 1, recent.end());
    recent.front() = Entry{&graph, graph.cacheNonce, graph.n, bondEnds, cycles};
    return cycles;
}

class RingConstraint {
    const MolGraph& graph_;
    std::shared_ptr<const CycleList> cycles_;
public:
    RingConstraint(const MolGraph& graph, bool enabled) : graph_(graph) {
        if (enabled) cycles_ = recentRelevantCycles(graph);
    }
    bool accepts(const std::set<int>& atoms, const std::set<int64_t>& bonds) const {
        std::set<int> coveredAtoms;
        std::set<int64_t> coveredBonds;
        static const CycleList none;
        for (const auto& cycle : cycles_ ? *cycles_ : none) {
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

/**
 * The stereo test of a match over its selected bonds; `embedding` is the substructure rule of
 * tetraParityCompatible() and doubleBondParityMatches().
 */
inline bool selectedStereoCompatible(const MolGraph& query, const MolGraph& target,
                                     const ChemOptions& chemistry, const Match& match, bool embedding = false) {
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
                    LocalMapping{mapping, selected, atom.first, -1}, embedding)) return false;
    if (chemistry.useBondStereo)
        for (const auto& bond : match.bondMapping)
            if (!stereo_detail::doubleBondParityMatches(query, bond[0], bond[1], target,
                    LocalMapping{mapping, selected, bond[0], bond[1]}, embedding)) return false;
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

/** completeRings() for a complete mapping that selects every query bond, the bonds being `pairs`. */
inline bool completeRings(const std::vector<int>& mapping, const std::vector<std::array<int, 2>>& pairs,
                          const RingConstraint& query, const RingConstraint& target) {
    std::set<int> queryAtoms, targetAtoms;
    std::set<int64_t> queryBonds, targetBonds;
    for (size_t atom = 0; atom < mapping.size(); ++atom) {
        queryAtoms.insert(static_cast<int>(atom));
        targetAtoms.insert(mapping[atom]);
    }
    for (const auto& pair : pairs) {
        queryBonds.insert(MolGraph::bondKey(pair[0], pair[1]));
        targetBonds.insert(MolGraph::bondKey(mapping[pair[0]], mapping[pair[1]]));
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
    /** Whether the matches are already in finish() order. */
    bool ordered_ = false;
    Score scoreOf(const Match& match) const {
        return options_ ? scoreMatch(match, *options_)
            : Score{static_cast<int64_t>(match.atomMapping.size()), static_cast<int>(match.atomMapping.size()), 0};
    }
public:
    Collector(const MCSOptions* options, size_t limit, bool all)
        : options_(options), limit_(limit), all_(all) {
        if (limit == 0) throw std::invalid_argument("Match limit must be positive");
    }
    bool remember(Match match) {
        assert(!ordered_ && "a collector that adopted matches takes no more");
        const Score score = scoreOf(match);
        if (score < best_) return false;
        if (best_ < score) { best_ = score; matches_.clear(); keys_.clear(); truncated_ = false; optimumKnown_ = false; }
        if (!all_ && limit_ == 1) {
            if (!matches_.empty())
                return match.atomMapping != matches_.front().atomMapping
                    || match.bondMapping != matches_.front().bondMapping;
            matches_.push_back(std::move(match));
            return true;
        }
        auto key = std::make_pair(match.atomMapping, match.bondMapping);
        if (keys_.count(key)) return false;
        if (matches_.size() < limit_) {
            keys_.insert(std::move(key));
            matches_.push_back(std::move(match));
        } else if (all_) truncated_ = true;
        return true;
    }
    /**
     * Takes the matches of a search that kept its own: distinct, in finish() order and those
     * remember() would have kept, with whether more were found. The collector holds no match
     * before and takes none after.
     */
    void adopt(std::vector<Match> matches, bool truncated) {
        assert(matches_.empty() && keys_.empty() && !ordered_ && "only an empty collector adopts matches");
        if (!matches.empty()) {
            const Score score = scoreOf(matches.front());
            if (best_ < score) best_ = score;
        }
        matches_ = std::move(matches);
        truncated_ = truncated;
        ordered_ = true;
    }
    const Score& best() const { return best_; }
    size_t limit() const { return limit_; }
    bool hasMatch() const { return !matches_.empty(); }
    bool truncated() const { return truncated_; }
    void markOptimumKnown() { optimumKnown_ = true; }
    bool optimumKnown() const { return optimumKnown_; }
    bool needs(const Score& upper) const { return all_ ? !(upper < best_) : best_ < upper; }
    SearchResult finish(Budget& budget, bool proven, bool exhausted) {
        SearchResult result;
        result.matches = std::move(matches_);
        if (!ordered_)
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

using TranspositionKey = std::array<int, 10>;

/** Every atom property that transposable() compares, so that atoms can also be grouped by them. */
inline TranspositionKey transpositionKey(const MolGraph& graph, int atom) {
    auto property = [atom](const auto& values, int fallback = 0) {
        return atom < static_cast<int>(values.size()) ? static_cast<int>(values[atom]) : fallback;
    };
    return {graph.atomicNum[atom], property(graph.formalCharge), property(graph.massNumber), property(graph.hydrogenCount),
            property(graph.atomClass), graph.ring[atom], graph.aromatic[atom], graph.degree[atom],
            property(graph.ringCount), property(graph.tautomerClass, -1)};
}

/** Whether exchanging an earlier atom with a later one preserves every compared property and bond. */
inline bool transposable(const MolGraph& graph, int first, int second, const std::vector<double>* weights) {
    if (graph.atomicNum[first] != graph.atomicNum[second] || graph.degree[first] != graph.degree[second]
        || transpositionKey(graph, first) != transpositionKey(graph, second)
        || (weights && !weights->empty() && (*weights)[first] != (*weights)[second])) return false;
    for (int neighbor : graph.neighbors[first]) {
        if (neighbor == second) continue;
        if (!graph.hasBond(second, neighbor) || graph.bondOrder(first, neighbor) != graph.bondOrder(second, neighbor)
            || graph.bondInRing(first, neighbor) != graph.bondInRing(second, neighbor)
            || graph.bondAromatic(first, neighbor) != graph.bondAromatic(second, neighbor)) return false;
    }
    return true;
}

/** Whether a graph carries stereochemistry, which atom transpositions would not preserve. */
inline bool hasStereoTags(const MolGraph& graph) {
    for (int tag : graph.tetraChirality) if (tag != 0) return true;
    if (!graph.hasDbStereo) return false;  // dbStereo is then zero for every bond
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom]) if (graph.dbStereo(atom, neighbor) != 0) return true;
    return false;
}

/**
 * Atom transpositions are used only while selecting one optimal witness. Each atom joins the
 * class of the earliest atom transposable with it, or keeps its own.
 */
inline std::vector<int> transpositionClasses(const MolGraph& graph, Budget& budget,
                                             const std::vector<double>* weights = nullptr) {
    if (hasStereoTags(graph)) return {};
    std::vector<int> classes = SpareVectors<int>::take(static_cast<size_t>(graph.n), 0);
    std::iota(classes.begin(), classes.end(), 0);
    bool hasTranspositions = false;
    for (int first = 0; first < graph.n && !budget.stopped(); ++first)
        for (int second = first + 1; second < graph.n; ++second) {
            if (classes[second] != second || !transposable(graph, first, second, weights)) continue;
            classes[second] = classes[first];
            hasTranspositions = true;
        }
    if (!hasTranspositions) {
        SpareVectors<int>::give(classes);
        return {};
    }
    return classes;
}

} // namespace search_detail
} // namespace smsd

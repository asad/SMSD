/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The vertex search behind molecule substructure search and the MCIS: atom-by-atom
 * backtracking over candidate domains.
 *
 * Part of search.hpp.
 */
#pragma once

#include "smsd/detail/search_results.hpp"
#include "smsd/detail/bridge_bound.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace smsd {
namespace search_detail {

/** Local matching work; exhaustion retains the original child. */
struct ChildBoundWork {
    size_t remaining;
    size_t used = 0;
    bool exhausted = false;
    explicit ChildBoundWork(size_t units) : remaining(units) {}
    bool spend(size_t units = 1) {
        if (exhausted || units > remaining) { exhausted = true; return false; }
        remaining -= units;
        used += units;
        return true;
    }
};

/** A structural matching allowance; predicate calls count as units, not elapsed time. */
inline std::optional<size_t> childBoundAllowance(int queryCount, int targetCount) {
    if (queryCount < 0 || targetCount < 0) return std::nullopt;
    const size_t query = static_cast<size_t>(queryCount), target = static_cast<size_t>(targetCount);
    const size_t words = (target + 63) >> 6, maximum = std::numeric_limits<size_t>::max();
    if (query > maximum / 2 || target > (maximum - query * 2) / 3) return std::nullopt;
    const size_t reset = query * 2 + target * 3;
    if (words > maximum - reset || target > maximum - reset - words) return std::nullopt;
    const size_t structural = reset + words + target;
    if (target > maximum - words - 1 || (query != 0 && query > maximum / query)) return std::nullopt;
    const size_t perRow = words + target + 1, square = query * query;
    if (square > (maximum - structural) / perRow) return std::nullopt;
    return structural + square * perRow;
}

/**
 * Maximum injective domain cardinality, capped by an existing upper bound. Each row gives
 * possible target bits; it may include incompatible pairs, but must include every legal pair.
 * An interrupted augmentation is unknown, never a smaller completed cardinality.
 */
template<class DomainWord, class Stopped>
std::optional<int> injectiveDomainBound(int queryCount, int targetCount, int cap,
                                      DomainWord domainWord, Stopped stopped,
                                      ChildBoundWork* allowance = nullptr) {
    if (cap <= 0) return 0;
    if (stopped()) return std::nullopt;
    if (allowance) {
        if (queryCount < 0 || targetCount < 0) return std::nullopt;
        const size_t query = static_cast<size_t>(queryCount), target = static_cast<size_t>(targetCount);
        const size_t maximum = std::numeric_limits<size_t>::max();
        if (query > maximum / 2 || target > (maximum - query * 2) / 3
            || !allowance->spend(query * 2 + target * 3)) return std::nullopt;
    }
    SpareVectors<int>::Lease queryMatchLease(static_cast<size_t>(queryCount), -1);
    SpareVectors<int>::Lease targetMatchLease(static_cast<size_t>(targetCount), -1);
    SpareVectors<int>::Lease parentLease(static_cast<size_t>(targetCount), -1);
    SpareVectors<int>::Lease seenLease(static_cast<size_t>(targetCount), -1);
    SpareVectors<int>::Lease queueLease(static_cast<size_t>(queryCount), 0);
    auto& queryMatch = queryMatchLease.items;
    auto& targetMatch = targetMatchLease.items;
    auto& parent = parentLease.items;
    auto& seen = seenLease.items;
    auto& queue = queueLease.items;
    const size_t words = (static_cast<size_t>(targetCount) + 63) >> 6;
    size_t work = 0;
    auto interrupted = [&]() {
        if (allowance && !allowance->spend()) return true;
        return ((work++ & 255u) == 0) && stopped();
    };
    int cardinality = 0;
    for (int root = 0; root < queryCount; ++root) {
        if (interrupted()) return std::nullopt;
        int head = 0, tail = 1, freeTarget = -1;
        queue[0] = root;
        while (head < tail && freeTarget < 0) {
            const int queryAtom = queue[head++];
            for (size_t word = 0; word < words && freeTarget < 0; ++word) {
                if (interrupted()) return std::nullopt;
                uint64_t remaining = domainWord(queryAtom, word);
                if (allowance && allowance->exhausted) return std::nullopt;
                if (allowance && word + 1 == words && (targetCount & 63) != 0)
                    remaining &= (uint64_t(1) << (targetCount & 63)) - 1;
                while (remaining != 0) {
                    if (interrupted()) return std::nullopt;
                    const int targetAtom = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                    remaining &= remaining - 1;
                    if (seen[targetAtom] == root) continue;
                    seen[targetAtom] = root;
                    parent[targetAtom] = queryAtom;
                    if (targetMatch[targetAtom] < 0) { freeTarget = targetAtom; break; }
                    queue[tail++] = targetMatch[targetAtom];
                }
            }
        }
        if (freeTarget < 0) continue;
        while (freeTarget >= 0) {
            if (interrupted()) return std::nullopt;
            const int queryAtom = parent[freeTarget];
            const int previousTarget = queryMatch[queryAtom];
            queryMatch[queryAtom] = freeTarget;
            targetMatch[freeTarget] = queryAtom;
            freeTarget = previousTarget;
        }
        if (++cardinality == cap) {
            if (stopped()) return std::nullopt;
            return cardinality;
        }
    }
    if (stopped()) return std::nullopt;
    return cardinality;
}

/** Bond upper for at most addedAtoms live vertices, with scores 2*degreeToMapped+degreeToLive. */
template<class Stopped>
std::optional<int64_t> topKIncidentBondUpper(int64_t mappedBonds,
        const std::vector<int64_t>& liveScores, int addedAtoms, Stopped stopped) {
    if (mappedBonds < 0 || addedAtoms < 0) return std::nullopt;
    if (stopped()) return std::nullopt;
    const size_t selected = std::min(static_cast<size_t>(addedAtoms), liveScores.size());
    if (selected == 0) return mappedBonds;
    SpareVectors<int64_t>::Lease largestLease(selected, 0);
    auto& largest = largestLease.items;
    const auto greater = [](int64_t first, int64_t second) { return first > second; };
    const int64_t maximum = std::numeric_limits<int64_t>::max();
    size_t count = 0;
    int64_t sum = 0;
    for (size_t index = 0; index < liveScores.size(); ++index) {
        if ((index & 255u) == 0 && stopped()) return std::nullopt;
        const int64_t score = liveScores[index];
        if (score < 0) return std::nullopt;
        if (count < selected) {
            if (score > maximum - sum) return std::nullopt;
            sum += score;
            largest[count++] = score;
            std::push_heap(largest.begin(), largest.begin() + count, greater);
        } else if (score > largest.front()) {
            sum -= largest.front();
            if (score > maximum - sum) return std::nullopt;
            sum += score;
            std::pop_heap(largest.begin(), largest.end(), greater);
            largest.back() = score;
            std::push_heap(largest.begin(), largest.end(), greater);
        }
    }
    if (mappedBonds > maximum - sum / 2 || stopped()) return std::nullopt;
    return mappedBonds + sum / 2;
}

class VertexSearch {
    static constexpr bool certifiedSymmetryEnabled_ = true;
    const QueryView& query_;
    const MolGraph& target_;
    const ChemOptions& chemistry_;
    const MCSOptions* common_;
    const Match* required_;
    bool induced_;
    bool all_;
    bool allowSkipping_;
    const bool symmetryGraphCertified_;
    /** A main-phase cardinality bound for one capped, unweighted fragment result. */
    const bool injectiveBoundAllowed_;
    /** Element/topology transpositions and index ties for the bounded fragment policy, in both phases. */
    const bool fragmentPolicy_;
    /** Whether true leaf pairs have equivalence-only labels and no fixed future images. */
    const bool leafExchangeAllowed_;
    std::unique_ptr<BridgeBound> bridgeBound_;
    uint64_t probeNodeLimit_;
    Budget& budget_;
    Collector& collector_;
    RingConstraint queryRings_, targetRings_;
    std::vector<uint64_t> domains_;
    std::vector<int> domainCounts_, mapping_;
    std::vector<int> queryTwinClass_;
    /**
     * Of a query atom: its mapped neighbours, its neighbour count, its places in frontierAtoms_
     * and mappedAtoms_ (-1 when absent), the latest earlier atom of its transposition class (-1
     * for none) and the number of earlier atoms of the class.
     */
    struct QueryAtomState {
        int frontier = 0, degree = 0, position = -1, mappedPosition = -1, earlierTwin = -1, earlierTwins = 0;
    };
    /** Candidate support, transposition class (-1 until found) and a seen mark of a target atom. */
    struct TargetAtomState { int support = 0, twin = -1; bool seen = false; };
    std::vector<QueryAtomState> queryState_;
    std::vector<TargetAtomState> targetState_;
    /** Unmapped query atoms with a mapped neighbour, and the mapped query atoms. */
    std::vector<int> frontierAtoms_, mappedAtoms_, twinChain_;
    /** Target atoms without neighbours by transposition key, then index; filled when twins are first needed. */
    std::vector<std::pair<TranspositionKey, int>> isolatedTargets_;
    std::vector<std::vector<int>> spareCandidates_;
    bool targetTwinsWanted_ = false;
    /** The connected MCIS rule exchanges only non-terminal target atoms. */
    bool connectedTargetTwins_ = false;
    int targetTwinState_ = 0;
    bool emptyDomain_ = false;
    /** Whether the budget stopped a frame before all its candidates were tested. */
    bool frameStopped_ = false;
    /** Whether each new MCIS pair is tested against the centres and double bonds it touches. */
    bool centrePairs_ = false, doubleBondPairs_ = false;
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
    /**
     * How a molecule substructure search for every mapping keeps its complete mappings: 0 until
     * the first, then 1 as rows that become matches once the search ends, or 2 as matches in the
     * collector. A search for one mapping accepts one at most and keeps it in the collector.
     */
    int rowState_;
    std::optional<CompleteMappingShape> shape_;
    std::vector<int> rowValues_;
    std::optional<CompleteMappingRows> rows_;

    static size_t domainSize(int queryAtoms, int targetAtoms) {
        const size_t words = (static_cast<size_t>(targetAtoms) + 63) >> 6;
        if (words != 0 && static_cast<size_t>(queryAtoms) > std::numeric_limits<size_t>::max() / words)
            throw std::length_error("Domain matrix is too large");
        return static_cast<size_t>(queryAtoms) * words;
    }
    size_t domainWords() const { return (static_cast<size_t>(target_.n) + 63) >> 6; }
    size_t domainIndex(int queryAtom, size_t word) const {
        return static_cast<size_t>(queryAtom) * domainWords() + word;
    }
    bool domainHas(int queryAtom, int targetAtom) const {
        return (domains_[domainIndex(queryAtom, targetAtom >> 6)] & (uint64_t(1) << (targetAtom & 63))) != 0;
    }
    bool feasible(int queryAtom, int targetAtom) const {
        if (used_[targetAtom] || !domainHas(queryAtom, targetAtom)) return false;
        // Every earlier atom of the transposition class is visited, by its chain of earlier twins.
        for (int previous = queryState_[queryAtom].earlierTwin; previous >= 0; previous = queryState_[previous].earlierTwin)
            if (mapping_[previous] == -1 || (mapping_[previous] >= 0 && mapping_[previous] >= targetAtom)) return false;
        if (!induced_) {
            // Only mapped query neighbours constrain an extension that may add target bonds.
            for (int neighbor : query_.neighbors(queryAtom)) {
                const int mappedTarget = mapping_[neighbor];
                if (mappedTarget >= 0 && query_.adjacent(queryAtom, neighbor)
                    && (!target_.hasBond(targetAtom, mappedTarget)
                        || !query_.bondMatches(queryAtom, neighbor, target_, targetAtom, mappedTarget, chemistry_)))
                    return false;
            }
            return true;
        }
        for (int previous : mappedAtoms_) {
            const bool queryBond = query_.adjacent(queryAtom, previous);
            const bool targetBond = target_.hasBond(targetAtom, mapping_[previous]);
            if (queryBond) {
                if (!targetBond || !query_.bondMatches(queryAtom, previous, target_, targetAtom, mapping_[previous], chemistry_))
                    return false;
            } else if (induced_ && targetBond) return false;
        }
        return true;
    }
    /** Complete-ring domains also agree on topological ring membership. */
    void restrictCompleteRingDomains() {
        if (!common_ || common_->objective != MCSObjective::MCIS
            || !chemistry_.completeRingsOnly || !query_.molecule
            || queryRings_.truncated() || targetRings_.truncated()) return;
        for (int queryAtom = 0; queryAtom < query_.size(); ++queryAtom) {
            const bool queryRing = queryRings_.ringAtom(queryAtom);
            int count = 0;
            for (size_t word = 0; word < domainWords(); ++word) {
                if ((word & 63) == 0 && budget_.stopped()) return;
                uint64_t& domain = domains_[domainIndex(queryAtom, word)];
                for (uint64_t remaining = domain; remaining != 0; remaining &= remaining - 1) {
                    const int offset = smsd::ctz64(remaining);
                    const int targetAtom = static_cast<int>(word * 64) + offset;
                    if (queryRing != targetRings_.ringAtom(targetAtom))
                        domain &= ~(uint64_t(1) << offset);
                }
                count += smsd::popcount64(domain);
            }
            domainCounts_[queryAtom] = count;
        }
    }
    /** Current adjacency and bond predicates only remove prospective targets. */
    uint64_t feasibleTargets(int queryAtom, size_t word) const {
        uint64_t candidates = domains_[domainIndex(queryAtom, word)] & ~unavailableTargetBits_[word];
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
    /** A completed strict atom shortfall cannot improve the current atom/bond incumbent. */
    bool childCannotReach(int queryAtom, int targetAtom) {
        if (!injectiveBoundAllowed_ || !collector_.hasMatch()) return false;
        const int selected = targetAtom >= 0 ? 1 : 0;
        const int needed = collector_.best().atoms - mapped_ - selected;
        if (needed <= 0) return false;
        const auto units = childBoundAllowance(query_.size(), target_.n);
        if (!units) return false;
        ChildBoundWork allowance(*units);
        const size_t words = domainWords();
        if (!allowance.spend(words + static_cast<size_t>(target_.n))) return false;
        SpareVectors<uint64_t>::Lease liveLease(words, ~uint64_t(0));
        auto& live = liveLease.items;
        if (words != 0 && (target_.n & 63) != 0)
            live.back() &= (uint64_t(1) << (target_.n & 63)) - 1;
        for (int atom = 0; atom < target_.n; ++atom) {
            if ((atom & 255) == 0 && budget_.stopped()) return false;
            if (used_[atom] || excluded_[atom] || atom == targetAtom)
                live[atom >> 6] &= ~(uint64_t(1) << (atom & 63));
        }
        size_t polledWork = allowance.used;
        auto stoppedHere = [&]() {
            if (allowance.used - polledWork < 256) return false;
            polledWork = allowance.used;
            if (!budget_.stopped()) return false;
            allowance.exhausted = true;
            return true;
        };
        const auto bound = injectiveDomainBound(query_.size(), target_.n, needed,
            [&](int otherQuery, size_t word) {
                if (!allowance.spend() || mapping_[otherQuery] != -2 || otherQuery == queryAtom)
                    return uint64_t(0);
                uint64_t candidates = domains_[domainIndex(otherQuery, word)] & live[word];
                auto restrictPair = [&](int previousQuery, int previousTarget) {
                    if (!allowance.spend()) return false;
                    const bool queryBond = query_.adjacent(otherQuery, previousQuery);
                    uint64_t adjacent = 0;
                    if (targetBitAdjacencyAvailable_) {
                        if (!allowance.spend()) return false;
                        adjacent = target_.adjLong[previousTarget][word];
                    } else for (int neighbor : target_.neighbors[previousTarget]) {
                        if (!allowance.spend()) return false;
                        if (stoppedHere()) return false;
                        if (static_cast<size_t>(neighbor >> 6) == word)
                            adjacent |= uint64_t(1) << (neighbor & 63);
                    }
                    if (stoppedHere()) return false;
                    candidates &= queryBond ? adjacent : ~adjacent;
                    if (queryBond) for (uint64_t remaining = candidates; remaining != 0; remaining &= remaining - 1) {
                        if (!allowance.spend()) return false;
                        if (stoppedHere()) return false;
                        const int candidate = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                        if (!query_.bondMatches(otherQuery, previousQuery, target_, candidate, previousTarget, chemistry_))
                            candidates &= ~(uint64_t(1) << (candidate & 63));
                    }
                    return true;
                };
                for (int previous : mappedAtoms_) {
                    if (candidates == 0) break;
                    if (!restrictPair(previous, mapping_[previous])) return uint64_t(0);
                }
                if (targetAtom >= 0 && candidates != 0 && !restrictPair(queryAtom, targetAtom)) return uint64_t(0);
                return candidates;
            }, [&]() { return budget_.stopped(); }, &allowance);
        return bound && *bound < needed && !budget_.stopped();
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
    /** A selected ring atom must remain on a relevant cycle whose vertices are all available. */
    bool restrictRingAvailability(const RingConstraint& rings, std::vector<uint8_t>& available,
                                  std::vector<uint8_t>& covered) {
        std::fill(covered.begin(), covered.end(), uint8_t(0));
        size_t work = 0;
        for (const auto& cycle : rings.cycles()) {
            bool viable = true;
            for (int atom : cycle) {
                if ((work++ & 255u) == 0 && budget_.stopped()) return false;
                if (!available[atom]) { viable = false; break; }
            }
            if (viable) for (int atom : cycle) {
                if ((work++ & 255u) == 0 && budget_.stopped()) return false;
                covered[atom] = 1;
            }
        }
        for (size_t atom = 0; atom < available.size(); ++atom) {
            if ((atom & 255u) == 0 && budget_.stopped()) return false;
            if (available[atom] && rings.ringAtom(static_cast<int>(atom)) && !covered[atom]) available[atom] = 0;
        }
        return true;
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
            for (size_t word = 0; word < domainWords(); ++word) {
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
        if (chemistry_.completeRingsOnly && query_.molecule && common_->objective == MCSObjective::MCIS
            && !queryRings_.truncated() && !targetRings_.truncated()) {
            if (!restrictRingAvailability(queryRings_, availableQuery_, reachableQuery_)
                || !restrictRingAvailability(targetRings_, availableTarget_, reachableTarget_)) return conservative;
            for (int queryAtom : mappedQueryAtoms_)
                if (!availableQuery_[queryAtom] || !availableTarget_[mapping_[queryAtom]])
                    return {std::numeric_limits<int64_t>::lowest(), 0, 0};
            atoms = std::min({atoms,
                static_cast<int>(std::count(availableQuery_.begin(), availableQuery_.end(), uint8_t(1))),
                static_cast<int>(std::count(availableTarget_.begin(), availableTarget_.end(), uint8_t(1)))});
        }
        const bool connected = common_->connectivity == MCSConnectivity::Connected || common_->maxFragments == 1;
        if (connected) {
            const int queryRoot = mappedQueryAtoms_.empty() ? -1 : mappedQueryAtoms_.front();
            const int queryBound = connectedBound(true, availableQuery_, reachableQuery_, queryRoot);
            const int targetBound = connectedBound(false, availableTarget_, reachableTarget_,
                                                    queryRoot < 0 ? -1 : mapping_[queryRoot]);
            if (queryBound < 0 || targetBound < 0) return conservative;
            atoms = std::min({atoms, queryBound, targetBound});
        }
        if (injectiveBoundAllowed_ && atoms > mapped_) {
            const auto matching = injectiveDomainBound(query_.size(), target_.n, atoms - mapped_,
                [&](int queryAtom, size_t word) {
                    return mapping_[queryAtom] == -2 && availableQuery_[queryAtom]
                        ? feasibleTargets(queryAtom, word) : uint64_t(0);
                }, [&]() { return budget_.stopped(); });
            if (!matching) return conservative;
            atoms = std::min(atoms, mapped_ + *matching);
        }
        int64_t queryBonds = 0, targetBonds = 0;
        int64_t mappedQueryBonds = 0, mappedTargetBonds = 0;
        SpareVectors<int64_t>::Lease queryScoresLease(injectiveBoundAllowed_ ? static_cast<size_t>(query_.size()) : 0, 0);
        SpareVectors<int64_t>::Lease targetScoresLease(injectiveBoundAllowed_ ? static_cast<size_t>(target_.n) : 0, 0);
        auto& queryScores = queryScoresLease.items;
        auto& targetScores = targetScoresLease.items;
        size_t bondWork = 0;
        for (int atom = 0; atom < query_.size(); ++atom) {
            if ((atom & 255) == 0 && budget_.stopped()) return conservative;
            if (!availableQuery_[atom] || (connected && !reachableQuery_[atom])) continue;
            for (int neighbor : query_.neighbors(atom)) {
                if (injectiveBoundAllowed_ && (bondWork++ & 255u) == 0 && budget_.stopped()) return conservative;
                if (!availableQuery_[neighbor] || (connected && !reachableQuery_[neighbor])) continue;
                if (atom < neighbor) ++queryBonds;
                if (injectiveBoundAllowed_) {
                    if (mapping_[atom] == -2) queryScores[atom] += mapping_[neighbor] >= 0 ? 2 : 1;
                    else if (mapping_[atom] >= 0 && mapping_[neighbor] >= 0 && atom < neighbor) ++mappedQueryBonds;
                }
            }
        }
        for (int atom = 0; atom < target_.n; ++atom) {
            if ((atom & 255) == 0 && budget_.stopped()) return conservative;
            if (!availableTarget_[atom] || (connected && !reachableTarget_[atom])) continue;
            for (int neighbor : target_.neighbors[atom]) {
                if (injectiveBoundAllowed_ && (bondWork++ & 255u) == 0 && budget_.stopped()) return conservative;
                if (!availableTarget_[neighbor] || (connected && !reachableTarget_[neighbor])) continue;
                if (atom < neighbor) ++targetBonds;
                if (injectiveBoundAllowed_) {
                    if (!used_[atom]) targetScores[atom] += used_[neighbor] ? 2 : 1;
                    else if (used_[neighbor] && atom < neighbor) ++mappedTargetBonds;
                }
            }
        }
        const Score base{common_->atomWeights.empty() ? atoms : conservative.primary, atoms,
                static_cast<int>(std::min<int64_t>({conservative.bonds, queryBonds, targetBonds,
                                                   int64_t(atoms) * (atoms - 1) / 2}))};
        if (!injectiveBoundAllowed_ || atoms < mapped_) return base;
        const auto queryUpper = topKIncidentBondUpper(mappedQueryBonds, queryScores,
            atoms - mapped_, [&]() { return budget_.stopped(); });
        if (!queryUpper) return budget_.stopped() ? conservative : base;
        const auto targetUpper = topKIncidentBondUpper(mappedTargetBonds, targetScores,
            atoms - mapped_, [&]() { return budget_.stopped(); });
        if (!targetUpper) return budget_.stopped() ? conservative : base;
        if (budget_.stopped()) return conservative;
        return {base.primary, base.atoms,
                static_cast<int>(std::min<int64_t>({base.bonds, *queryUpper, *targetUpper}))};
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
                for (int previous : mappedAtoms_)
                    if (target_.hasBond(candidate, mapping_[previous]) != (uniformAdjacency_ != 0)) {
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
        if (bridgeBound_) {
            int bridgeAtoms = bridgeBound_->upper();
            for (int queryAtom : mappedAtoms_)
                bridgeAtoms = std::min(bridgeAtoms, bridgeBound_->bound(queryAtom, mapping_[queryAtom]));
            atoms = std::min(atoms, bridgeAtoms);
        }
        for (int atom = 0; atom < query_.size(); ++atom) {
            if (mapping_[atom] == -1) continue;
            for (int neighbor : query_.neighbors(atom))
                if (neighbor > atom && mapping_[neighbor] != -1) ++bonds;
        }
        if (bridgeBound_ && bridgeBound_->acyclic()) bonds = std::min(bonds, std::max(0, atoms - 1));
        const Score conservative{common_ && !common_->atomWeights.empty() ? weightedScore(possibleWeight) : atoms,
            atoms, static_cast<int>(std::min<int64_t>({bonds, targetEdges_, int64_t(atoms) * (atoms - 1) / 2}))};
        if (!collector_.needs(conservative)) return conservative;
        return common_ && allowSkipping_ ? feasibleUpper(conservative) : conservative;
    }
    /** Whether a graph stores the bond between two listed neighbours alike from either end. */
    static bool sameBothWays(const MolGraph& graph, int atom, int neighbor) {
        return !graph.useDense || graph.bondPacked[atom][neighbor] == graph.bondPacked[neighbor][atom];
    }
    /**
     * Whether rows may stand for the matches of complete mappings. Each query atom must list
     * other atoms, each once, bonded to it and listing it in return, so a complete mapping passed
     * the bond test for every query bond. Each target atom must list other atoms, each once, so
     * no complete mapping is found twice, as the collector would drop a repeat. Both graphs must
     * store each listed bond alike from either end, so the bond test reads it the same way from
     * the later atom as mappedBonds() does from the lower. The checks read each listed
     * neighbour, and for the query its list too.
     */
    bool rowsApply() const {
        const MolGraph& query = *query_.molecule;
        SpareVectors<int>::Lease seenLease(static_cast<size_t>(std::max(query.n, target_.n)), -1);
        std::vector<int>& seen = seenLease.items;
        for (int atom = 0; atom < query.n; ++atom)
            for (int neighbor : query.neighbors[atom]) {
                if (neighbor < 0 || neighbor >= query.n || neighbor == atom || seen[neighbor] == atom
                    || !query.hasBond(atom, neighbor) || !sameBothWays(query, atom, neighbor)) return false;
                seen[neighbor] = atom;
                const auto& back = query.neighbors[neighbor];
                if (std::find(back.begin(), back.end(), atom) == back.end()) return false;
            }
        std::fill(seen.begin(), seen.end(), -1);
        for (int atom = 0; atom < target_.n; ++atom)
            for (int neighbor : target_.neighbors[atom]) {
                if (neighbor < 0 || neighbor >= target_.n || neighbor == atom || seen[neighbor] == atom
                    || !sameBothWays(target_, atom, neighbor)) return false;
                seen[neighbor] = atom;
            }
        return true;
    }
    /** Whether complete mappings are kept as rows, decided at the first one. */
    bool rowsKept() {
        if (rowState_ == 0) {
            rowState_ = 2;
            if (rowsApply()) {
                shape_.emplace(query_.size(), query_.molecule->neighbors);
                rowValues_ = SpareVectors<int>::take();
                rows_.emplace(rowValues_, static_cast<size_t>(query_.size()), collector_.limit(), all_);
                rowState_ = 1;
            }
        }
        return rowState_ == 1;
    }
    /** Whether the bond test selects every query bond of the current complete mapping, as rows assume. */
    bool everyBondSelected() const {
        auto bonds = mappedBonds(query_, target_, mapping_, chemistry_);
        std::sort(bonds.begin(), bonds.end());
        std::vector<std::array<int, 4>> expected;
        for (const auto& pair : shape_->pairs) expected.push_back({pair[0], pair[1], mapping_[pair[0]], mapping_[pair[1]]});
        return bonds == expected;
    }
    /** record() for a complete mapping kept as a row, with the stereo and ring checks read from the mapping. */
    bool recordRow() {
        assert(everyBondSelected() && "a complete mapping must select every query bond");
        // A substructure search embeds the query, so its specified stereo needs specified stereo.
        if (!mappingStereoCompatible(*query_.molecule, target_, chemistry_, mapping_, !common_)) return false;
        if (chemistry_.completeRingsOnly && !completeRings(mapping_, shape_->pairs, queryRings_, targetRings_)) return false;
        rows_->append(mapping_);
        return true;
    }
    bool truncated() const { return rows_ ? rows_->truncated : collector_.truncated(); }
    bool found() const { return rows_ ? rows_->count() != 0 : collector_.hasMatch(); }
    bool record() {
        if (!common_ && mapped_ != query_.size()) return false;
        if (rowState_ != 2 && rowsKept()) return recordRow();
        if (common_ && mapped_ == 0) return false;
        if (common_ && common_->atomWeights.empty() && mapped_ < collector_.best().atoms) return false;
        if (required_) for (const auto& atom : required_->atomMapping)
            if (mapping_[atom.first] != atom.second) return false;
        Match match = makeMatch(mapping_, mappedBonds(query_, target_, mapping_, chemistry_));
        if (common_) {
            if (common_->connectivity == MCSConnectivity::Connected && match.fragments.size() > 1) return false;
            if (match.fragments.size() > static_cast<size_t>(common_->maxFragments)) return false;
            for (const auto& fragment : match.fragments)
                if (fragment.size() < static_cast<size_t>(common_->minFragmentAtoms)) return false;
        }
        if (query_.molecule && !selectedStereoCompatible(*query_.molecule, target_, chemistry_, match, !common_)) return false;
        if (query_.pattern && !query_.pattern->mappingCompatible(target_, mapping_)) return false;
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
    static bool tagged(const MolGraph& graph) {
        return std::any_of(graph.tetraChirality.begin(), graph.tetraChirality.end(), [](int tag) { return tag != 0; });
    }
    static bool configured(const MolGraph& graph) {
        if (!graph.hasDbStereo) return false;
        for (int atom = 0; atom < graph.n; ++atom)
            for (int neighbor : graph.neighbors[atom])
                if (neighbor > atom && graph.dbStereo(atom, neighbor) != 0) return true;
        return false;
    }
    /**
     * Whether pairing `atom` with `candidate` fails the stereo test of record() at a centre or
     * double bond that reads the pair: the atom itself, its mapped neighbours, and the double
     * bonds at either. Mapped pairs stay in every extension, so such a failure stays too, and the
     * candidate is passed over before its node, as in Java's VertexSearch.stereoBroken().
     */
    bool stereoBroken(int atom, int candidate) {
        mapping_[atom] = candidate;
        const bool broken = stereoBroken(atom);
        mapping_[atom] = -2;
        return broken;
    }
    bool stereoBroken(int atom) const {
        const auto& neighbors = query_.molecule->neighbors;
        if (centrePairs_) {
            if (centreBroken(atom)) return true;
            for (int neighbor : neighbors[atom])
                if (mapping_[neighbor] >= 0 && centreBroken(neighbor)) return true;
        }
        if (doubleBondPairs_)
            for (int neighbor : neighbors[atom]) {
                if (mapping_[neighbor] < 0) continue;
                if (doubleBondBroken(atom, neighbor)) return true;
                for (int other : neighbors[neighbor])
                    if (other != atom && mapping_[other] >= 0 && doubleBondBroken(neighbor, other)) return true;
            }
        return false;
    }
    /**
     * The tetrahedral test at a mapped centre. Once three ligands are mapped the permutation is
     * fixed, as a fourth can only take the remaining target ligand, and a ligand whose image is
     * not a target ligand stays so.
     */
    bool centreBroken(int centre) const {
        return !tetraParityCompatible(*query_.molecule, centre, target_, mapping_[centre], mapping_);
    }
    /**
     * The double-bond test between two mapped atoms. It reads the first mapped substituent at
     * each end; since the mapping is induced, any other gives the same relation when an end has
     * one substituent in either graph, or at most two in both, and the bond is left to record()
     * otherwise.
     */
    bool doubleBondBroken(int first, int second) const {
        const MolGraph& query = *query_.molecule;
        for (int end : {first, second}) {
            const int querySubstituents = static_cast<int>(query.neighbors[end].size()) - 1;
            const int targetSubstituents = static_cast<int>(target_.neighbors[mapping_[end]].size()) - 1;
            if (std::min(querySubstituents, targetSubstituents) > 1 && std::max(querySubstituents, targetSubstituents) > 2)
                return false;
        }
        return !stereo_detail::doubleBondParityMatches(query, first, second, target_, mapping_);
    }
    struct Frame {
        int atom;
        std::vector<int> candidates;
        size_t next = 0;
        int active = -2;
        bool skipped = false;
    };
    /** Omit a skip only when an already searched image extends every skip completion. */
    bool skipDominated(const Frame& frame) {
        if (!injectiveBoundAllowed_ || frame.active != -2 || frame.next != frame.candidates.size()
            || !common_->excludedTargetAtoms.empty() || !chemistry_.matchAtomType
            || chemistry_.matchFormalCharge || chemistry_.matchIsotope || chemistry_.useChirality
            || chemistry_.useBondStereo || chemistry_.tautomerAware || chemistry_.ringMatchesRingOnly
            || chemistry_.completeRingsOnly || chemistry_.ringFusionMode != ChemOptions::RingFusionMode::IGNORE
            || chemistry_.matchBondOrder != ChemOptions::BondOrderMode::ANY
            || chemistry_.aromaticityMode != ChemOptions::AromaticityMode::FLEXIBLE) return false;
        size_t work = 0;
        auto stopped = [&]() { return ((work++ & 255u) == 0) && budget_.stopped(); };
        for (int neighbor : query_.neighbors(frame.atom)) {
            if (stopped()) return false;
            if (mapping_[neighbor] == -2) return false;
        }
        for (int targetAtom : frame.candidates) {
            if (stopped()) return false;
            if (!feasible(frame.atom, targetAtom)) continue;
            bool closed = true;
            for (int neighbor : target_.neighbors[targetAtom]) {
                if (stopped()) return false;
                if (!used_[neighbor] && !excluded_[neighbor]) { closed = false; break; }
            }
            if (!closed) continue;
            bool unique = true;
            for (int otherQuery = 0; otherQuery < query_.size(); ++otherQuery) {
                if (stopped()) return false;
                if (otherQuery == frame.atom || mapping_[otherQuery] != -2
                    || !domainHas(otherQuery, targetAtom)) continue;
                // Ignore the earlier-twin inequalities: absence must hold in the relaxed row.
                bool available = true;
                for (int previous : mappedAtoms_) {
                    if (stopped()) return false;
                    const bool queryBond = query_.adjacent(otherQuery, previous);
                    const bool targetBond = target_.hasBond(targetAtom, mapping_[previous]);
                    if (queryBond != targetBond || (queryBond
                        && !query_.bondMatches(otherQuery, previous, target_, targetAtom, mapping_[previous], chemistry_))) {
                        available = false;
                        break;
                    }
                }
                if (available) { unique = false; break; }
            }
            if (unique && !budget_.stopped()) return true;
        }
        return false;
    }
    /** Fix a true leaf pair beside a mapped parent without losing an optimal extension. */
    bool forcePendantPair(Frame& frame) {
        if (!leafExchangeAllowed_ || queryState_[frame.atom].degree != 1) return false;
        const int queryParent = query_.neighbors(frame.atom).front();
        const int targetParent = mapping_[queryParent];
        if (targetParent < 0) return false;
        size_t tested = 0;
        for (int targetAtom : frame.candidates) {
            if ((++tested & 255u) == 0 && budget_.stopped()) return false;
            if (target_.neighbors[targetAtom].size() == 1 && target_.neighbors[targetAtom].front() == targetParent) {
                // Exchanging the partners preserves all selected edges and the bond tie.
                frame.candidates.front() = targetAtom;
                frame.candidates.resize(1);
                frame.skipped = true;
                return true;
            }
        }
        return false;
    }
    /** Whether target transpositions may thin candidates; checked when first needed. */
    bool targetTwinsAvailable() {
        if (targetTwinState_ == 0) {
            targetTwinState_ = targetTwinsWanted_ ? 2 : 1;
            if (targetTwinState_ == 2) {
                isolatedTargets_ = SpareVectors<std::pair<TranspositionKey, int>>::take();
                for (int atom = 0; atom < target_.n; ++atom)
                    if (target_.neighbors[atom].empty()) isolatedTargets_.emplace_back(transpositionKey(target_, atom, fragmentPolicy_), atom);
                std::sort(isolatedTargets_.begin(), isolatedTargets_.end());
            }
        }
        return targetTwinState_ == 2;
    }
    /**
     * The class transpositionClasses assigns a target atom: that of the earliest atom
     * transposable with it. Such an atom shares every other neighbour, so it lies within two
     * bonds or, without neighbours of its own, among the isolated atoms. An isolated atom is
     * transposable exactly with the atoms that share its key, so the earliest of those is looked
     * up instead of testing each isolated atom. -1, with nothing remembered on the way, when
     * the budget stops the search first.
     */
    int targetTwin(int atom) {
        if (connectedTargetTwins_ && target_.neighbors[atom].size() < 2)
            return targetState_[atom].twin = atom;
        twinChain_.clear();
        int work = 0;
        while (targetState_[atom].twin < 0) {
            int earliest = -1;
            auto consider = [&](int other) {
                if (other < atom && (earliest < 0 || other < earliest) && transposable(target_, other, atom, nullptr, fragmentPolicy_))
                    earliest = other;
            };
            if (!isolatedTargets_.empty()) {
                const auto key = transpositionKey(target_, atom, fragmentPolicy_);
                const auto isolated = std::lower_bound(isolatedTargets_.begin(), isolatedTargets_.end(), std::make_pair(key, INT_MIN));
                if (isolated != isolatedTargets_.end() && isolated->first == key) consider(isolated->second);
            }
            for (int neighbor : target_.neighbors[atom]) {
                consider(neighbor);
                for (int other : target_.neighbors[neighbor]) consider(other);
                // The clock is read after about 256 candidates, however many atoms a chain passes.
                work += 1 + static_cast<int>(target_.neighbors[neighbor].size());
                if (work >= 256) {
                    work = 0;
                    if (budget_.stopped()) return -1;
                }
            }
            if (++work >= 256) {
                work = 0;
                if (budget_.stopped()) return -1;
            }
            if (earliest < 0) {
                targetState_[atom].twin = atom;
                break;
            }
            twinChain_.push_back(atom);
            atom = earliest;
        }
        for (int member : twinChain_) targetState_[member].twin = targetState_[atom].twin;
        return targetState_[atom].twin;
    }
    /** Lists an atom exactly while it is unmapped and has a mapped neighbour. */
    void refreshFrontier(int atom) {
        auto& state = queryState_[atom];
        const bool wanted = mapping_[atom] == -2 && state.frontier > 0;
        if (wanted == (state.position >= 0)) return;
        if (wanted) {
            state.position = static_cast<int>(frontierAtoms_.size());
            frontierAtoms_.push_back(atom);
        } else {
            const int last = frontierAtoms_.back();
            frontierAtoms_[state.position] = last;
            queryState_[last].position = state.position;
            frontierAtoms_.pop_back();
            state.position = -1;
        }
    }
    /** Maps or unmaps a query atom, keeping the mapped atoms and the mapped-neighbour counts of its neighbours. */
    void assign(int queryAtom, int targetAtom) {
        const int delta = (targetAtom >= 0) - (mapping_[queryAtom] >= 0);
        mapping_[queryAtom] = targetAtom;
        if (delta != 0) {
            for (int neighbor : query_.neighbors(queryAtom)) {
                queryState_[neighbor].frontier += delta;
                refreshFrontier(neighbor);
            }
            auto& state = queryState_[queryAtom];
            if (delta > 0) {
                state.mappedPosition = static_cast<int>(mappedAtoms_.size());
                mappedAtoms_.push_back(queryAtom);
            } else {
                const int last = mappedAtoms_.back();
                mappedAtoms_[state.mappedPosition] = last;
                queryState_[last].mappedPosition = state.mappedPosition;
                mappedAtoms_.pop_back();
                state.mappedPosition = -1;
            }
        }
        refreshFrontier(queryAtom);
    }
    /** Keeps a finished candidate list for the next frame of this search. */
    void recycle(std::vector<int>& candidates) {
        if (candidates.capacity() == 0) return;
        candidates.clear();
        spareCandidates_.push_back(std::move(candidates));
    }
    /**
     * An empty candidate list, reused from this search or else from this thread's spares, so
     * deep or backtracking searches, and the many short searches of a batch, do not allocate
     * on each node.
     */
    std::vector<int> spareCandidates() {
        if (spareCandidates_.empty()) return SpareVectors<int>::take();
        std::vector<int> spare = std::move(spareCandidates_.back());
        spareCandidates_.pop_back();
        return spare;
    }
    Frame frame() {
        // An unmapped atom with a mapped neighbour comes first; then fewer candidates, more
        // neighbours and the lower index decide.
        auto better = [&](int atom, int chosen) {
            const int degree = queryState_[atom].degree, chosenDegree = queryState_[chosen].degree;
            return domainCounts_[atom] < domainCounts_[chosen]
                || (domainCounts_[atom] == domainCounts_[chosen]
                    && (degree > chosenDegree || (degree == chosenDegree && atom < chosen)));
        };
        int chosen = -1, bestFrontier = -1;
        for (int atom : frontierAtoms_)
            if (chosen < 0 || better(atom, chosen)) chosen = atom;
        if (chosen >= 0) bestFrontier = queryState_[chosen].frontier;
        else for (int atom = 0; atom < query_.size(); ++atom)
            if (mapping_[atom] == -2 && (chosen < 0 || better(atom, chosen))) {
                chosen = atom;
                bestFrontier = queryState_[atom].frontier;
            }
        if (chosen < 0) return Frame{-1, {}};
        if (common_ && (common_->connectivity == MCSConnectivity::Connected || common_->maxFragments == 1)
            && mapped_ > 0 && bestFrontier == 0)
            return Frame{-1, {}};
        Frame result{chosen, spareCandidates()};
        int neighborAnchor = -1;
        for (int neighbor : query_.neighbors(chosen)) {
            if (mapping_[neighbor] < 0) continue;
            const size_t targetDegree = target_.neighbors[mapping_[neighbor]].size();
            if (neighborAnchor < 0
                || targetDegree < target_.neighbors[mapping_[neighborAnchor]].size()
                || (fragmentPolicy_ && targetDegree == target_.neighbors[mapping_[neighborAnchor]].size()
                    && neighbor < neighborAnchor))
                neighborAnchor = neighbor;
        }
        result.candidates.reserve(neighborAnchor >= 0
            ? target_.neighbors[mapping_[neighborAnchor]].size() : static_cast<size_t>(domainCounts_[chosen]));
        // A test compares the candidate with each earlier twin and with each mapped atom or
        // neighbour; the clock is read after about 4096 such comparisons. Once the budget stops,
        // the frame keeps the candidates found so far and the search ends at its next node.
        const int64_t testCost = 1 + queryState_[chosen].earlierTwins
            + static_cast<int64_t>(induced_ ? mappedAtoms_.size() : query_.neighbors(chosen).size());
        int64_t work = 0;
        auto stopAfter = [&](int64_t units) {
            if (frameStopped_) return true;
            if ((work += units) < 4096) return false;
            work = 0;
            frameStopped_ = frameStopped_ || budget_.stopped();
            return frameStopped_;
        };
        if (neighborAnchor >= 0) {
            for (int candidate : target_.neighbors[mapping_[neighborAnchor]]) {
                if (stopAfter(testCost)) break;
                if (feasible(chosen, candidate)) result.candidates.push_back(candidate);
            }
        } else {
            for (size_t word = 0; word < domainWords() && !frameStopped_; ++word) {
                uint64_t remaining = domains_[domainIndex(chosen, word)];
                while (remaining != 0 && !stopAfter(testCost)) {
                    const int candidate = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                    remaining &= remaining - 1;
                    if (feasible(chosen, candidate)) result.candidates.push_back(candidate);
                }
            }
        }
        if (forcePendantPair(result)) return result;
        if (query_.molecule) {
            if (targetState_.empty())
                targetState_ = SpareVectors<TargetAtomState>::take(static_cast<size_t>(target_.n), TargetAtomState{});
            auto support = [&](int candidate) -> int& { return targetState_[candidate].support; };
            for (int candidate : result.candidates) {
                if (stopAfter(1 + static_cast<int64_t>(query_.neighbors(chosen).size()) * target_.neighbors[candidate].size()))
                    break;
                for (int neighbor : query_.neighbors(chosen)) {
                    if (mapping_[neighbor] == -1) continue;
                    for (int adjacent : target_.neighbors[candidate])
                        if (!used_[adjacent] && domainHas(neighbor, adjacent)
                            && query_.bondMatches(chosen, neighbor, target_, candidate, adjacent, chemistry_)) {
                            ++support(candidate); break;
                        }
                }
            }
            std::stable_sort(result.candidates.begin(), result.candidates.end(), [&](int first, int second) {
                if (support(first) != support(second)) return support(first) > support(second);
                const auto queryDegree = static_cast<int>(query_.neighbors(chosen).size());
                const int firstDistance = std::abs(target_.degree[first] - queryDegree);
                const int secondDistance = std::abs(target_.degree[second] - queryDegree);
                if (firstDistance != secondDistance) return firstDistance < secondDistance;
                return fragmentPolicy_ && first < second;
            });
            for (int candidate : result.candidates) support(candidate) = 0;
            if (!all_ && !required_ && result.candidates.size() > 1 && targetTwinsAvailable()) {
                if (!queryTwinClass_.empty()) {
                    // Query inequalities need the smallest feasible value of each target class.
                    auto clearMinimums = [&]() {
                        for (int candidate : result.candidates) {
                            const int twin = targetState_[candidate].twin;
                            if (twin >= 0) targetState_[twin].support = 0;
                        }
                    };
                    size_t work = 0;
                    for (int candidate : result.candidates) {
                        if ((work++ & 255u) == 0 && budget_.stopped()) {
                            clearMinimums();
                            return result;
                        }
                        const int twin = targetTwin(candidate);
                        if (twin < 0) {
                            clearMinimums();
                            return result;
                        }
                        int& minimum = targetState_[twin].support;
                        if (minimum == 0 || candidate < minimum - 1) minimum = candidate + 1;
                    }
                    if (budget_.stopped()) {
                        clearMinimums();
                        return result;
                    }
                    size_t kept = 0;
                    work = 0;
                    // Finish the linear compaction so class scratch is cleared after a stop.
                    for (int candidate : result.candidates) {
                        if ((work++ & 255u) == 0) (void)budget_.stopped();
                        int& minimum = targetState_[targetState_[candidate].twin].support;
                        if (minimum == 0) continue;
                        result.candidates[kept++] = minimum - 1;
                        minimum = 0;
                    }
                    result.candidates.resize(kept);
                    return result;
                }
                // Keep the first candidate of each class; the marks are cleared before returning.
                // Once the budget stops, the remaining candidates are kept as they are, as the
                // search ends at its next node.
                size_t kept = 0, marked = 0;
                bool stopped = false;
                for (size_t index = 0; index < result.candidates.size(); ++index) {
                    const int candidate = result.candidates[index];
                    const int twin = stopped ? -1 : targetTwin(candidate);
                    if (twin < 0) {
                        stopped = true;
                        result.candidates[kept++] = candidate;
                        continue;
                    }
                    bool& seen = targetState_[twin].seen;
                    if (seen) continue;
                    seen = true;
                    result.candidates[kept++] = candidate;
                    ++marked;
                }
                result.candidates.resize(kept);
                for (size_t index = 0; index < marked; ++index)
                    targetState_[targetState_[result.candidates[index]].twin].seen = false;
            }
        }
        return result;
    }
    /**
     * Fills the rows the general loop in the constructor would give a molecule query whose atoms
     * match only target atoms of their own element. Candidates come from per-element and
     * per-degree target bit sets, and query atoms with equal atomCompatible terms share one row
     * unless a required or anchored image singles one out. A substructure stops at its first
     * empty row, since no embedding then exists. False, with nothing filled, when an atomic
     * number lies outside the element table.
     */
    bool fillElementDomains(int anchor, const Match* required) {
        const MolGraph& query = *query_.molecule;
        constexpr int elements = 128;
        std::array<int, elements> slot;
        slot.fill(-1);
        int slots = 0, maximumDegree = 0;
        for (int atom = 0; atom < query.n; ++atom) {
            const int element = query.atomicNum[atom];
            if (element < 0 || element >= elements) return false;
            if (slot[element] < 0) slot[element] = slots++;
            maximumDegree = std::max(maximumDegree, queryState_[atom].degree);
        }
        if (allowSkipping_) maximumDegree = 0;
        const size_t words = domainWords();
        // Element rows come first. Degree row k then holds the target atoms of degree k or more;
        // the last degree row also takes every higher degree. A final row holds the shared domain.
        SpareVectors<uint64_t>::Lease bitsLease((static_cast<size_t>(slots) + maximumDegree + 2) * words, 0);
        std::vector<uint64_t>& bits = bitsLease.items;
        uint64_t* const degreeBits = bits.data() + static_cast<size_t>(slots) * words;
        uint64_t* const shared = degreeBits + static_cast<size_t>(maximumDegree + 1) * words;
        for (int atom = 0; atom < target_.n; ++atom) {
            if (excluded_[atom]) continue;
            const int element = target_.atomicNum[atom];
            const uint64_t bit = uint64_t(1) << (atom & 63);
            if (element >= 0 && element < elements && slot[element] >= 0) bits[slot[element] * words + (atom >> 6)] |= bit;
            const int degree = static_cast<int>(std::min(target_.neighbors[atom].size(), static_cast<size_t>(maximumDegree)));
            degreeBits[degree * words + (atom >> 6)] |= bit;
        }
        for (int degree = maximumDegree - 1; degree >= 0; --degree)
            for (size_t word = 0; word < words; ++word) degreeBits[degree * words + word] |= degreeBits[(degree + 1) * words + word];
        // Where candidates are tested one by one, atoms sorted by the terms atomCompatible reads
        // share a row; an element-only row costs too little to share.
        const bool elementOnly = atomCompatibleByElementOnly(chemistry_);
        SpareVectors<std::array<int, 10>>::Lease termsLease;
        SpareVectors<int>::Lease orderLease;
        std::vector<std::array<int, 10>>& terms = termsLease.items;
        std::vector<int>& order = orderLease.items;
        if (!elementOnly) {
            terms.reserve(query.n);
            for (int atom = 0; atom < query.n; ++atom) terms.push_back(atomCompatibilityTerms(query, atom));
            order.resize(query.n);
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int first, int second) {
                return std::tie(terms[first], first) < std::tie(terms[second], second);
            });
        }
        // The clock is read after 256 candidate tests or 64 words of bit operations.
        int work = 0;
        auto stopped = [&](int units) {
            if ((work += units) < 64) return false;
            work = 0;
            return budget_.stopped();
        };
        auto member = [&](size_t index) { return elementOnly ? static_cast<int>(index) : order[index]; };
        for (size_t begin = 0, end = 0; begin < static_cast<size_t>(query.n); begin = end) {
            const int representative = member(begin);
            for (end = begin + 1; !elementOnly && end < order.size()
                && terms[order[end]] == terms[representative]
                && queryState_[order[end]].degree == queryState_[representative].degree; ++end) {}
            const uint64_t* elementRow = bits.data() + static_cast<size_t>(slot[query.atomicNum[representative]]) * words;
            const uint64_t* degreeRow = degreeBits
                + static_cast<size_t>(std::min(queryState_[representative].degree, maximumDegree)) * words;
            int count = 0;
            for (size_t word = 0; word < words; ++word) {
                if (stopped(elementOnly ? 1 : 16)) return true;
                uint64_t candidates = elementRow[word] & degreeRow[word];
                if (!elementOnly) {
                    for (uint64_t remaining = candidates; remaining != 0; remaining &= remaining - 1) {
                        const int offset = smsd::ctz64(remaining);
                        if (!query_.atomMatches(representative, target_, static_cast<int>(word * 64) + offset, chemistry_))
                            candidates &= ~(uint64_t(1) << offset);
                    }
                }
                shared[word] = candidates;
                count += smsd::popcount64(candidates);
            }
            for (size_t index = begin; index < end; ++index) {
                const int atom = member(index);
                int image = -1;
                bool single = false;
                if (required) {
                    const auto mapped = required->atomMapping.find(atom);
                    if (mapped != required->atomMapping.end()) { single = true; image = mapped->second; }
                }
                if (atom == 0 && anchor >= 0) {
                    image = single && image != anchor ? -1 : anchor;
                    single = true;
                }
                uint64_t* row = domains_.data() + domainIndex(atom, 0);
                if (!single) {
                    std::copy(shared, shared + words, row);
                    domainCounts_[atom] = count;
                } else if (image >= 0 && image < target_.n && (shared[image >> 6] & (uint64_t(1) << (image & 63))) != 0) {
                    row[image >> 6] = uint64_t(1) << (image & 63);
                    domainCounts_[atom] = 1;
                }
                if (!common_ && domainCounts_[atom] == 0) {
                    emptyDomain_ = true;
                    return true;
                }
            }
        }
        return true;
    }
public:
    VertexSearch(const QueryView& query, const MolGraph& target, const ChemOptions& chemistry,
                 const MCSOptions* common, bool induced, bool all, Budget& budget, Collector& collector, int anchor = -1,
                 const Match* required = nullptr, bool allowSkipping = true, uint64_t probeNodeLimit = 0,
                 bool rowsAllowed = true, bool allowInjectiveBound = false, bool allowFragmentPolicy = false)
        : query_(query), target_(target), chemistry_(chemistry), common_(common), required_(required), induced_(induced), all_(all),
          allowSkipping_(common && allowSkipping),
          symmetryGraphCertified_(certifiedSymmetryEnabled_ && query.molecule && !all && !required && anchor < 0
              && symmetryGraphCoherent(*query.molecule, [&budget] { return budget.stoppedSampled(); })
              && symmetryGraphCoherent(target, [&budget] { return budget.stoppedSampled(); })),
          injectiveBoundAllowed_(allowInjectiveBound && symmetryGraphCertified_ && allowSkipping_
              && induced && !all && collector.limit() == 1 && !required && probeNodeLimit == 0
              && common && common->objective == MCSObjective::MCIS && common->atomWeights.empty()
              && common->connectivity == MCSConnectivity::Fragments
              && common->minFragmentAtoms == 1 && common->maxFragments == INT_MAX
              && common->nodeLimit != 0),
          fragmentPolicy_(allowFragmentPolicy && symmetryGraphCertified_ && induced && !all
              && collector.limit() == 1 && !required && anchor < 0 && common
              && common->objective == MCSObjective::MCIS && common->atomWeights.empty()
              && common->connectivity == MCSConnectivity::Fragments
              && common->minFragmentAtoms == 1 && common->maxFragments == INT_MAX
              && common->nodeLimit != 0 && common->excludedTargetAtoms.empty()
              && chemistry.matchAtomType && !chemistry.matchFormalCharge && !chemistry.matchIsotope
              && !chemistry.useChirality && !chemistry.useBondStereo && !chemistry.tautomerAware
              && !chemistry.ringMatchesRingOnly && !chemistry.completeRingsOnly
              && chemistry.ringFusionMode == ChemOptions::RingFusionMode::IGNORE
              && chemistry.matchBondOrder == ChemOptions::BondOrderMode::ANY
              && chemistry.aromaticityMode == ChemOptions::AromaticityMode::FLEXIBLE),
          leafExchangeAllowed_(symmetryGraphCertified_ && common && induced
              && common->objective == MCSObjective::MCIS && common->connectivity == MCSConnectivity::Connected
              && common->atomWeights.empty() && !chemistry.matchIsotope && !chemistry.useChirality
              && !chemistry.useBondStereo && !chemistry.tautomerAware && !chemistry.completeRingsOnly
              && !chemistry.ringMatchesRingOnly && chemistry.ringFusionMode == ChemOptions::RingFusionMode::IGNORE
              && chemistry.aromaticityMode == ChemOptions::AromaticityMode::STRICT
              && chemistry.matchBondOrder == ChemOptions::BondOrderMode::STRICT),
          probeNodeLimit_(probeNodeLimit),
          budget_(budget), collector_(collector),
          queryRings_(query.molecule ? *query.molecule : target, chemistry.completeRingsOnly && query.molecule),
          targetRings_(target, chemistry.completeRingsOnly),
          domains_(SpareVectors<uint64_t>::take(domainSize(query.size(), target.n), 0)),
          domainCounts_(SpareVectors<int>::take(static_cast<size_t>(query.size()), 0)),
          mapping_(SpareVectors<int>::take(static_cast<size_t>(query.size()), -2)),
          used_(SpareVectors<bool>::take(static_cast<size_t>(target.n), false)),
          excluded_(SpareVectors<bool>::take(static_cast<size_t>(target.n), false)), availableTargets_(target.n),
          rowState_(rowsAllowed && all && !common && !required && query.molecule ? 0 : 2) {
        queryState_ = SpareVectors<QueryAtomState>::take(static_cast<size_t>(query.size()), QueryAtomState{});
        {
            SpareVectors<int>::Lease seenLease(static_cast<size_t>(query.size()), -1);
            auto& seen = seenLease.items;
            size_t work = 0;
            for (int atom = 0; atom < query.size(); ++atom) {
                if (budget.stopped()) return;
                for (int neighbor : query.neighbors(atom)) {
                    if ((work++ & 255u) == 0 && budget.stopped()) return;
                    if (neighbor != atom && seen[neighbor] != atom) {
                        seen[neighbor] = atom;
                        ++queryState_[atom].degree;
                    }
                }
            }
        }
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
                if (targetBitAdjacencyAvailable_) {
                    for (int atom = 0; atom < target.n && targetBitAdjacencyAvailable_; ++atom) {
                        for (size_t word = 0; word < availableTargetBits_.size(); ++word) {
                            if (budget.stoppedSampled()) return;
                            uint64_t remaining = target.adjLong[atom][word];
                            while (remaining != 0) {
                                const int neighbour = static_cast<int>(word * 64) + smsd::ctz64(remaining);
                                remaining &= remaining - 1;
                                if (neighbour >= target.n || !target.hasBond(atom, neighbour)) {
                                    targetBitAdjacencyAvailable_ = false;
                                    break;
                                }
                            }
                        }
                        for (int neighbour : target.neighbors[atom]) {
                            if (budget.stoppedSampled()) return;
                            if (neighbour != atom && (target.adjLong[atom][neighbour >> 6]
                                & (uint64_t(1) << (neighbour & 63))) == 0) {
                                targetBitAdjacencyAvailable_ = false;
                                break;
                            }
                        }
                    }
                }
            }
            bool complete = true, independent = true;
            for (int atom = 0; atom < query.size(); ++atom) {
                complete = complete && queryState_[atom].degree == query.size() - 1;
                independent = independent && queryState_[atom].degree == 0;
            }
            if (complete || independent) uniformAdjacency_ = complete ? 1 : 0;
            if (induced && query.molecule) {
                centrePairs_ = chemistry.useChirality && tagged(*query.molecule) && tagged(target);
                doubleBondPairs_ = chemistry.useBondStereo && configured(*query.molecule) && configured(target);
            }
        }
        frontierAtoms_ = SpareVectors<int>::take();
        frontierAtoms_.reserve(query.size());
        mappedAtoms_ = SpareVectors<int>::take();
        mappedAtoms_.reserve(query.size());
        const bool transpositionsAllowed = symmetryGraphCertified_ && !chemistry.completeRingsOnly
            && !chemistry.ringMatchesRingOnly && chemistry.ringFusionMode == ChemOptions::RingFusionMode::IGNORE
            && (!common || (common->connectivity == MCSConnectivity::Fragments
                && common->minFragmentAtoms == 1 && common->maxFragments == INT_MAX));
        connectedTargetTwins_ = leafExchangeAllowed_ && chemistry.matchAtomType
            && common->excludedTargetAtoms.empty() && query.size() <= 512 && target.n <= 512
            && std::all_of(target.neighbors.begin(), target.neighbors.end(), [](const auto& neighbors) {
                return neighbors.size() <= 8;
            });
        if (query.molecule && !all && !required && (transpositionsAllowed || connectedTargetTwins_)
            && (fragmentPolicy_ || (!hasStereoTags(*query.molecule) && !hasStereoTags(target)))) {
            targetTwinsWanted_ = true;
            if (transpositionsAllowed && common && common->objective == MCSObjective::MCIS && common->atomWeights.empty())
                queryTwinClass_ = transpositionClasses(*query.molecule, budget, &common->atomWeights, fragmentPolicy_);
            if (!queryTwinClass_.empty()) {
                // Each atom links to the latest earlier atom of its class, which it counts with the rest.
                std::vector<int> latest(static_cast<size_t>(query.size()), -1), members(static_cast<size_t>(query.size()), 0);
                for (int atom = 0; atom < query.size(); ++atom) {
                    const int twinClass = queryTwinClass_[atom];
                    queryState_[atom].earlierTwin = latest[twinClass];
                    queryState_[atom].earlierTwins = members[twinClass]++;
                    latest[twinClass] = atom;
                }
            }
        }
        if (common) for (int atom = 0; atom < target.n; ++atom)
            for (int neighbor : target.neighbors[atom]) targetEdges_ += neighbor > atom;
        if (query.molecule && chemistry.matchAtomType && fillElementDomains(anchor, required)) {
            restrictCompleteRingDomains();
            return;
        }
        for (int atom = 0; atom < query.size() && !budget.stopped(); ++atom) {
            for (int candidate = 0; candidate < target.n; ++candidate) {
                if ((candidate & 255) == 0 && budget.stopped()) break;
                if (excluded_[candidate]) continue;
                if (required && required->atomMapping.count(atom) && required->atomMapping.at(atom) != candidate) continue;
                if (atom == 0 && anchor >= 0 && candidate != anchor) continue;
                if (!allowSkipping_ && static_cast<size_t>(queryState_[atom].degree) > target.neighbors[candidate].size()) continue;
                if (query.atomMatches(atom, target, candidate, chemistry)) {
                    domains_[domainIndex(atom, candidate >> 6)] |= uint64_t(1) << (candidate & 63);
                    ++domainCounts_[atom];
                }
            }
        }
        restrictCompleteRingDomains();
    }
    VertexSearch(const VertexSearch&) = delete;
    VertexSearch& operator=(const VertexSearch&) = delete;
    ~VertexSearch() {
        SpareVectors<uint64_t>::give(domains_);
        for (auto* spare : {&domainCounts_, &mapping_, &queryTwinClass_, &frontierAtoms_, &mappedAtoms_, &twinChain_, &rowValues_})
            SpareVectors<int>::give(*spare);
        SpareVectors<std::pair<TranspositionKey, int>>::give(isolatedTargets_);
        for (auto& spare : spareCandidates_) SpareVectors<int>::give(spare);
        spareCandidates_.clear();
        SpareVectors<std::vector<int>>::give(spareCandidates_);
        SpareVectors<QueryAtomState>::give(queryState_);
        SpareVectors<TargetAtomState>::give(targetState_);
        SpareVectors<bool>::give(used_);
        SpareVectors<bool>::give(excluded_);
    }
    /**
     * Runs the search once: whether its result is proven, then whether it ran to its end. Kept
     * rows then become the collector's matches.
     */
    std::pair<bool, bool> run() {
        const auto proof = search();
        if (rows_) collector_.adopt(rows_->build(*shape_), rows_->truncated);
        return proof;
    }
private:
    std::pair<bool, bool> search() {
        if (budget_.stopped()) return {false, false};
        if (query_.size() == 0) { if (!common_) record(); return {true, true}; }
        if (!common_ && (query_.size() > target_.n || emptyDomain_)) return {true, true};
        Score absoluteUpper = common_ ? upper() : Score{};
        if (query_.molecule && common_ && induced_ && allowSkipping_ && common_->atomWeights.empty()
            && common_->objective == MCSObjective::MCIS
            && (common_->connectivity == MCSConnectivity::Connected || common_->maxFragments == 1)
            && query_.size() >= 64 && target_.n >= 64 && collector_.needs(absoluteUpper)) {
            bridgeBound_ = BridgeBound::create(*query_.molecule, target_, chemistry_, budget_);
            if (budget_.stopped()) return {false, false};
            if (bridgeBound_) absoluteUpper = upper();
        }
        bool optimumProven = false;
        const uint64_t startingNodes = budget_.nodesVisited();
        std::vector<Frame> pending;
        // Hand the candidate lists still on the stack back to this thread however the search ends.
        struct Recycler {
            std::vector<Frame>& frames;
            ~Recycler() { for (auto& frame : frames) SpareVectors<int>::give(frame.candidates); }
        } recycler{pending};
        // Each frame on the stack holds a distinct query atom, and each candidate list is on the
        // stack or spare, so neither list grows past these sizes.
        pending.reserve(static_cast<size_t>(query_.size()) + 1);
        if (spareCandidates_.capacity() == 0) spareCandidates_ = SpareVectors<std::vector<int>>::take();
        spareCandidates_.reserve(static_cast<size_t>(query_.size()) + 1);
        pending.push_back(frame());
        // A substructure step reads the clock on every sixteenth node; between reads it does at
        // most sixteen extensions and their bounded backtracking. Cancellation is checked at each node.
        while (!pending.empty() && (!common_ || !budget_.stopped())) {
            if (probeNodeLimit_ && budget_.nodesVisited() - startingNodes >= probeNodeLimit_) return {optimumProven, false};
            auto& current = pending.back();
            if (current.active >= 0) {
                used_[current.active] = false;
                --mapped_;
                assign(current.atom, -2);
                current.active = -2;
            } else if (current.active == -1) {
                assign(current.atom, -2);
                current.active = -2;
            }
            if (current.next < current.candidates.size()) {
                if ((centrePairs_ || doubleBondPairs_) && stereoBroken(current.atom, current.candidates[current.next])) {
                    ++current.next;
                    continue;
                }
                if (childCannotReach(current.atom, current.candidates[current.next])) {
                    ++current.next;
                    continue;
                }
                if (!budget_.visit(common_ || (budget_.nodesVisited() & 15) == 0)) break;
                const int candidate = current.candidates[current.next++];
                assign(current.atom, candidate);
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
                if (!common_ && truncated()) return {true, false};
                if (common_ && !collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.atom >= 0) pending.push_back(std::move(next));
            } else if (allowSkipping_ && !current.skipped && (!required_ || !required_->atomMapping.count(current.atom))) {
                if (skipDominated(current)) {
                    current.skipped = true;
                    continue;
                }
                if (childCannotReach(current.atom, -1)) {
                    current.skipped = true;
                    continue;
                }
                if (!budget_.visit()) break;
                current.skipped = true;
                current.active = -1;
                assign(current.atom, -1);
                if (!collector_.needs(upper())) continue;
                Frame next = frame();
                if (next.atom >= 0) pending.push_back(std::move(next));
            } else {
                recycle(current.candidates);
                pending.pop_back();
            }
        }
        // A frame the budget stopped was not searched through, so the stack running empty proves nothing then.
        const bool exhausted = pending.empty() && !frameStopped_;
        return {exhausted || optimumProven || collector_.optimumKnown() || (!common_ && found()), exhausted};
    }
};

} // namespace search_detail
} // namespace smsd

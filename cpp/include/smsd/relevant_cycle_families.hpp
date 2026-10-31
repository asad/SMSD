/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Relevant cycles grouped into families, with the exact number of members of each.
 *
 * ring_finder.hpp includes this header at its end; including it first includes ring_finder.hpp.
 */
#pragma once
#ifndef SMSD_RELEVANT_CYCLE_FAMILIES_HPP
#define SMSD_RELEVANT_CYCLE_FAMILIES_HPP

#include "smsd/ring_finder.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace smsd {

/**
 * A family of relevant cycles. Each family gives its ring size, the highest-numbered atom
 * shared by its members (root), the atoms opposite it, one member and the exact number of
 * members. For an odd size, first and second are the atoms of the bond opposite root and
 * opposite is -1; for an even size, opposite is the atom opposite root and first and second
 * are its ring neighbours. first is less than second, and the prototype member starts at root
 * and reaches first before second.
 */
struct RelevantCycleFamily {
    int size = 0;
    int root = -1;
    int first = -1;
    int second = -1;
    int opposite = -1;              // -1 for odd sizes
    std::vector<int> prototype;
    std::string memberCount;        // exact, in decimal
};

/**
 * Relevant cycles grouped into families, in ascending size, then root. relevantCycleCount is
 * the exact number of members of the families listed, in decimal; truncated reports an
 * incomplete enumeration, whose families are still relevant but may not be all of them. When
 * the enumeration is complete, the relevant cycles and their number do not depend on atom
 * order; the grouping into families, and so each family's record, can.
 */
struct RelevantCycleFamilies {
    std::vector<RelevantCycleFamily> families;
    std::string relevantCycleCount = "0";
    bool truncated = false;
};

namespace detail_rf {

/** Work allowance of one family enumeration, counted as computeRelevantCycleFamiliesUncached() describes. */
inline constexpr uint64_t RELEVANT_CYCLE_FAMILY_WORK_LIMIT = uint64_t(1) << 29;

/** The most words of shortest-path atom sets that one root of a family enumeration may keep. */
inline constexpr uint64_t RELEVANT_CYCLE_FAMILY_PATH_WORDS = uint64_t(1) << 23;

/** A non-negative integer in base 2^32, lowest limb first and without leading zero limbs. */
class ExactCount {
public:
    ExactCount() = default;
    explicit ExactCount(uint32_t value) { if (value != 0) limbs_.push_back(value); }

    /** Limbs in use: the bit length divided by 32, rounded up. */
    size_t limbs() const noexcept { return limbs_.size(); }

    void add(const ExactCount& other) {
        if (limbs_.size() < other.limbs_.size()) limbs_.resize(other.limbs_.size(), 0);
        uint64_t carry = 0;
        for (size_t i = 0; i < limbs_.size(); ++i) {
            const uint64_t sum = uint64_t(limbs_[i]) + (i < other.limbs_.size() ? other.limbs_[i] : 0) + carry;
            limbs_[i] = static_cast<uint32_t>(sum);
            carry = sum >> 32;
            if (carry == 0 && i + 1 >= other.limbs_.size()) break;
        }
        if (carry != 0) limbs_.push_back(static_cast<uint32_t>(carry));
    }

    ExactCount times(const ExactCount& other) const {
        ExactCount product;
        if (limbs_.empty() || other.limbs_.empty()) return product;
        product.limbs_.assign(limbs_.size() + other.limbs_.size(), 0);
        for (size_t i = 0; i < limbs_.size(); ++i) {
            uint64_t carry = 0;
            for (size_t j = 0; j < other.limbs_.size(); ++j) {
                const uint64_t sum = uint64_t(limbs_[i]) * other.limbs_[j] + product.limbs_[i + j] + carry;
                product.limbs_[i + j] = static_cast<uint32_t>(sum);
                carry = sum >> 32;
            }
            product.limbs_[i + other.limbs_.size()] = static_cast<uint32_t>(carry);
        }
        while (!product.limbs_.empty() && product.limbs_.back() == 0) product.limbs_.pop_back();
        return product;
    }

    std::string decimal() const {
        if (limbs_.empty()) return "0";
        std::vector<uint32_t> rest(limbs_);
        std::vector<uint32_t> groups;   // base 10^9, lowest first
        while (!rest.empty()) {
            uint64_t remainder = 0;
            for (size_t i = rest.size(); i-- > 0; ) {
                const uint64_t value = (remainder << 32) | rest[i];
                rest[i] = static_cast<uint32_t>(value / 1000000000u);
                remainder = value % 1000000000u;
            }
            groups.push_back(static_cast<uint32_t>(remainder));
            while (!rest.empty() && rest.back() == 0) rest.pop_back();
        }
        std::string text = std::to_string(groups.back());
        for (size_t i = groups.size() - 1; i-- > 0; ) {
            const std::string group = std::to_string(groups[i]);
            text.append(9 - group.size(), '0').append(group);
        }
        return text;
    }

private:
    std::vector<uint32_t> limbs_;
};

/** Whether each listed neighbour of every atom lists the atom back, once. */
inline bool mutualNeighbourLists(const MolGraph& g, const std::unordered_map<int64_t, int>& edgeIndex) {
    std::vector<uint8_t> ends(edgeIndex.size(), 0);
    for (int u = 0; u < g.n; ++u)
        for (int w : g.neighbors[u]) {
            const auto it = edgeIndex.find(MolGraph::bondKey(u, w));
            if (it == edgeIndex.end() || ++ends[it->second] > 2) return false;
        }
    for (uint8_t count : ends) if (count != 2) return false;
    return true;
}

/**
 * The breadth-first levels of one root r of a family enumeration over its ring system S:
 * levels in S, and levels and shortest-path parents among the atoms of S numbered at most r,
 * the parents in the order the search reads them. An atom keeps both levels equal exactly
 * when its shortest paths among those atoms are shortest in S.
 */
struct FamilyLevels {
    std::vector<int> inSystem, belowRoot, order, reachedInSystem;
    std::vector<std::vector<int>> parents;
    uint64_t work = 0;

    void search(const MolGraph& g, const std::vector<int>& ofAtom, int root, bool withParents) {
        const size_t n = static_cast<size_t>(g.n);
        if (inSystem.size() != n) {
            inSystem.assign(n, -1);
            belowRoot.assign(n, -1);
            parents.assign(n, {});
        }
        for (int u : reachedInSystem) inSystem[u] = -1;
        for (int u : order) { belowRoot[u] = -1; parents[u].clear(); }
        const int s = ofAtom[root];
        reachedInSystem.assign(1, root);
        inSystem[root] = 0;
        for (size_t head = 0; head < reachedInSystem.size(); ++head) {
            const int u = reachedInSystem[head];
            work += 1 + g.neighbors[u].size();
            for (int w : g.neighbors[u])
                if (ofAtom[w] == s && inSystem[w] < 0) { inSystem[w] = inSystem[u] + 1; reachedInSystem.push_back(w); }
        }
        order.assign(1, root);
        belowRoot[root] = 0;
        for (size_t head = 0; head < order.size(); ++head) {
            const int u = order[head];
            work += 1 + g.neighbors[u].size();
            for (int w : g.neighbors[u]) {
                if (ofAtom[w] != s || w > root) continue;
                if (belowRoot[w] < 0) {
                    belowRoot[w] = belowRoot[u] + 1;
                    order.push_back(w);
                    if (withParents) parents[w].push_back(u);
                } else if (withParents && belowRoot[w] == belowRoot[u] + 1) {
                    parents[w].push_back(u);
                }
            }
        }
    }

    /** Both levels agree, so the atom's shortest paths below the root are shortest in its system. */
    bool shortest(int atom) const { return belowRoot[atom] >= 0 && belowRoot[atom] == inSystem[atom]; }

    /** The member whose paths to first and second follow first-read parents. */
    std::vector<int> firstCycle(int root, int first, int opposite, int second) const {
        std::vector<int> cycle;
        for (int u = first; ; u = parents[u].front()) { cycle.push_back(u); if (u == root) break; }
        std::reverse(cycle.begin(), cycle.end());
        if (opposite >= 0) cycle.push_back(opposite);
        for (int u = second; u != root; u = parents[u].front()) cycle.push_back(u);
        return cycle;
    }
};

/**
 * Relevant-cycle families of g, each system in turn and each root of a system in ascending
 * order (only the highest atom of a system that is a simple ring). A root r reads shortest
 * paths among the atoms numbered at most r whose levels agree with those of the whole system;
 * pairs of such paths that meet only at r, ending at an adjacent pair of one level or at two
 * parents of one atom, give one candidate family each, its first member and its exact count.
 * A candidate family is relevant when its first member is independent of the candidates of
 * strictly smaller size, which span all shorter cycles.
 *
 * Work, counted identically in the Java implementation: per root, the atoms and neighbours
 * that its two searches read; per shortest-path parent, the words of an atom set plus the
 * limbs of the larger count added; per pair tested, the words of an atom set; per candidate,
 * its length plus the product of the limbs of the two counts multiplied; per row added in a
 * reduction, the words of a bond set. The allowance is tested after each root and after each
 * reduction; when it is exceeded, or a root's atom sets would exceed
 * RELEVANT_CYCLE_FAMILY_PATH_WORDS, the result is truncated: the ring system in progress and
 * every later one are left out, and the families of the systems before it are kept.
 */
inline RelevantCycleFamilies computeRelevantCycleFamiliesUncached(
        const MolGraph& g, uint64_t workLimit = RELEVANT_CYCLE_FAMILY_WORK_LIMIT, uint64_t* workDone = nullptr) {
    RelevantCycleFamilies result;
    if (workDone) *workDone = 0;
    if (g.n == 0 || detail::cycleRank(g) <= 0) return result;
    const auto systems = g.cycleSystems();
    const auto edgeIndex = buildEdgeIndex(g);
    if (systems.ofAtom.empty() || !mutualNeighbourLists(g, edgeIndex)) {
        result.truncated = true;
        return result;
    }
    const auto& ofAtom = systems.ofAtom;
    const size_t systemCount = systems.atoms.size();
    std::vector<std::vector<int>> systemAtoms(systemCount);
    std::vector<int> localAtom(static_cast<size_t>(g.n), -1);
    for (int u = 0; u < g.n; ++u)
        if (ofAtom[u] >= 0) {
            localAtom[u] = static_cast<int>(systemAtoms[ofAtom[u]].size());
            systemAtoms[ofAtom[u]].push_back(u);
        }
    const size_t bondWords = (edgeIndex.size() + 63) >> 6;

    struct Candidate { int size, root, sequence, first, second, opposite; std::vector<int> cycle; ExactCount count; };
    std::vector<Candidate> kept;
    FamilyLevels levels;
    uint64_t work = 0;
    bool stopped = false;
    std::vector<uint64_t> atomSets, rows, bits;
    std::vector<ExactCount> counts;
    std::vector<int> position(static_cast<size_t>(g.n), -1);
    std::vector<int> rowByPivot(edgeIndex.size(), -1), pivotsUsed;

    for (size_t s = 0; s < systemCount && !stopped; ++s) {
        const auto& atoms = systemAtoms[s];
        const size_t setWords = (atoms.size() + 63) >> 6;
        if (static_cast<uint64_t>(atoms.size()) * setWords > RELEVANT_CYCLE_FAMILY_PATH_WORDS) {
            stopped = true;
            break;
        }
        std::vector<Candidate> found;
        const size_t firstRoot = systems.simpleRing[s] ? atoms.size() - 1 : 0;
        for (size_t at = firstRoot; at < atoms.size(); ++at) {
            const int root = atoms[at];
            levels.work = 0;
            levels.search(g, ofAtom, root, true);
            work += levels.work;
            const auto& order = levels.order;
            atomSets.assign(order.size() * setWords, 0);
            counts.assign(order.size(), ExactCount());
            for (size_t i = 0; i < order.size(); ++i) position[order[i]] = static_cast<int>(i);
            for (size_t i = 0; i < order.size(); ++i) {
                const int u = order[i];
                uint64_t* set = atomSets.data() + i * setWords;
                set[static_cast<size_t>(localAtom[u]) >> 6] |= uint64_t(1) << (localAtom[u] & 63);
                if (u == root) counts[i] = ExactCount(1);
                for (int p : levels.parents[u]) {
                    const size_t from = static_cast<size_t>(position[p]);
                    const uint64_t* parentSet = atomSets.data() + from * setWords;
                    for (size_t k = 0; k < setWords; ++k) set[k] |= parentSet[k];
                    work += setWords + std::max(counts[i].limbs(), counts[from].limbs());
                    counts[i].add(counts[from]);
                }
            }
            const size_t rootWord = static_cast<size_t>(localAtom[root]) >> 6;
            const uint64_t rootBit = uint64_t(1) << (localAtom[root] & 63);
            auto meetOnlyAtRoot = [&](int y, int z) {
                work += setWords;
                const uint64_t* first = atomSets.data() + static_cast<size_t>(position[y]) * setWords;
                const uint64_t* second = atomSets.data() + static_cast<size_t>(position[z]) * setWords;
                for (size_t k = 0; k < setWords; ++k) {
                    uint64_t shared = first[k] & second[k];
                    if (k == rootWord) shared &= ~rootBit;
                    if (shared != 0) return false;
                }
                return true;
            };
            int sequence = 0;
            auto keep = [&](int size, int y, int z, int x) {
                Candidate candidate{size, root, sequence++, y, z, x, levels.firstCycle(root, y, x, z), ExactCount()};
                const ExactCount& ofFirst = counts[static_cast<size_t>(position[y])];
                const ExactCount& ofSecond = counts[static_cast<size_t>(position[z])];
                work += candidate.cycle.size() + ofFirst.limbs() * ofSecond.limbs();
                candidate.count = ofFirst.times(ofSecond);
                found.push_back(std::move(candidate));
            };
            for (int y : order) {
                if (y == root || !levels.shortest(y)) continue;
                for (int z : g.neighbors[y]) {
                    if (z <= y || z > root || ofAtom[z] != static_cast<int>(s) || !levels.shortest(z)
                        || levels.belowRoot[z] != levels.belowRoot[y]) continue;
                    if (meetOnlyAtRoot(y, z)) keep(2 * levels.belowRoot[y] + 1, y, z, -1);
                }
            }
            for (int x : order) {
                if (x == root || !levels.shortest(x) || levels.belowRoot[x] < 2) continue;
                const auto& parents = levels.parents[x];
                for (size_t i = 0; i < parents.size(); ++i)
                    for (size_t j = i + 1; j < parents.size(); ++j) {
                        const int y = std::min(parents[i], parents[j]), z = std::max(parents[i], parents[j]);
                        if (meetOnlyAtRoot(y, z)) keep(2 * levels.belowRoot[x], y, z, x);
                    }
            }
            for (int u : order) position[u] = -1;
            if (work > workLimit) { stopped = true; break; }
        }
        if (stopped) break;

        // Relevance: each size is tested against the rows of strictly smaller sizes, then added.
        std::stable_sort(found.begin(), found.end(),
                         [](const Candidate& a, const Candidate& b) { return a.size < b.size; });
        for (int pivot : pivotsUsed) rowByPivot[pivot] = -1;
        pivotsUsed.clear();
        rows.clear();
        auto bondsOf = [&](const std::vector<int>& cycle) {
            bits.assign(bondWords, 0);
            for (size_t i = 0; i < cycle.size(); ++i) {
                const int bond = edgeIndex.at(edgeKey(cycle[i], cycle[(i + 1) % cycle.size()]));
                bits[static_cast<size_t>(bond) >> 6] ^= uint64_t(1) << (bond & 63);
            }
        };
        auto reduce = [&]() {
            for (size_t word = 0; word < bondWords; ++word)
                while (bits[word] != 0) {
                    const int pivot = static_cast<int>((word << 6) + static_cast<size_t>(ctz64(bits[word])));
                    const int row = rowByPivot[pivot];
                    if (row < 0) return pivot;
                    const uint64_t* added = rows.data() + static_cast<size_t>(row) * bondWords;
                    for (size_t k = word; k < bondWords; ++k) bits[k] ^= added[k];
                    work += bondWords;
                }
            return -1;
        };
        std::vector<char> relevant(found.size(), 0);
        for (size_t i = 0; i < found.size() && !stopped; ) {
            size_t j = i;
            while (j < found.size() && found[j].size == found[i].size) ++j;
            for (size_t k = i; k < j; ++k) {
                bondsOf(found[k].cycle);
                relevant[k] = reduce() >= 0;
                if (work > workLimit) { stopped = true; break; }
            }
            if (stopped) break;
            for (size_t k = i; k < j; ++k) {
                bondsOf(found[k].cycle);
                const int pivot = reduce();
                if (pivot >= 0) {
                    rowByPivot[pivot] = static_cast<int>(rows.size() / bondWords);
                    rows.insert(rows.end(), bits.begin(), bits.end());
                    pivotsUsed.push_back(pivot);
                }
                if (work > workLimit) { stopped = true; break; }
            }
            i = j;
        }
        // A system whose relevance test ran out is left out whole, as the Java implementation does.
        if (stopped) break;
        for (size_t k = 0; k < found.size(); ++k)
            if (relevant[k]) kept.push_back(std::move(found[k]));
    }

    std::stable_sort(kept.begin(), kept.end(), [](const Candidate& a, const Candidate& b) {
        if (a.size != b.size) return a.size < b.size;
        if (a.root != b.root) return a.root < b.root;
        return a.sequence < b.sequence;
    });
    ExactCount total;
    for (auto& candidate : kept) {
        total.add(candidate.count);
        RelevantCycleFamily family;
        family.size = candidate.size;
        family.root = candidate.root;
        family.first = candidate.first;
        family.second = candidate.second;
        family.opposite = candidate.opposite;
        family.prototype = std::move(candidate.cycle);
        family.memberCount = candidate.count.decimal();
        result.families.push_back(std::move(family));
    }
    result.relevantCycleCount = total.decimal();
    result.truncated = stopped;
    if (workDone) *workDone = work;
    return result;
}

} // namespace detail_rf

/**
 * The relevant cycles of g grouped into families, cached on the graph. Each family gives its
 * ring size, the highest-numbered atom shared by its members, the atoms opposite it, one
 * member and the exact number of members; truncated reports an incomplete enumeration.
 */
inline RelevantCycleFamilies computeRelevantCycleFamilies(const MolGraph& g) {
    if (!g.lazyCaches_.ready(MolGraph::LazyCaches::CycleFamilies)) {
        const std::lock_guard<std::recursive_mutex> lock(g.lazyCacheMutex());
        if (!g.cachedRelevantCycleFamilies_)
            g.cachedRelevantCycleFamilies_ = std::make_shared<const RelevantCycleFamilies>(
                detail_rf::computeRelevantCycleFamiliesUncached(g));
        g.lazyCaches_.publish(MolGraph::LazyCaches::CycleFamilies);
    }
    return *g.cachedRelevantCycleFamilies_;
}

/**
 * The position in families.families of the family that has ring as a member, or -1 when ring
 * is not a cycle of g or not a member of a listed family. The test is exact.
 */
inline int relevantCycleFamilyOf(const MolGraph& g, const RelevantCycleFamilies& families,
                                 const std::vector<int>& ring) {
    const size_t length = ring.size();
    if (length < 3) return -1;
    std::set<int> atoms;
    for (size_t i = 0; i < length; ++i) {
        const int u = ring[i], w = ring[(i + 1) % length];
        if (u < 0 || u >= g.n || w < 0 || w >= g.n || !atoms.insert(u).second) return -1;
        const auto& row = g.neighbors[u];
        if (std::find(row.begin(), row.end(), w) == row.end()) return -1;
    }
    const auto systems = g.cycleSystems();
    if (systems.ofAtom.empty() || systems.ofAtom[ring[0]] < 0) return -1;
    for (int u : ring) if (systems.ofAtom[u] != systems.ofAtom[ring[0]]) return -1;
    const size_t start = static_cast<size_t>(std::max_element(ring.begin(), ring.end()) - ring.begin());
    const int root = ring[start];
    detail_rf::FamilyLevels levels;
    levels.search(g, systems.ofAtom, root, false);
    // Both arcs from the root advance one level per atom, levels below the root equal to levels in the system.
    const size_t half = length / 2;
    auto at = [&](size_t offset) { return ring[(start + offset) % length]; };
    for (size_t i = 0; i <= half; ++i)
        for (int u : {at(i), at((length - i) % length)})
            if (!levels.shortest(u) || levels.belowRoot[u] != static_cast<int>(i)) return -1;
    int first, second, opposite = -1;
    if (length % 2 == 1) {
        first = at(half);
        second = at(length - half);
    } else {
        opposite = at(half);
        first = at(half - 1);
        second = at(length - half + 1);
    }
    if (first > second) std::swap(first, second);
    for (size_t i = 0; i < families.families.size(); ++i) {
        const auto& family = families.families[i];
        if (family.size == static_cast<int>(length) && family.root == root && family.first == first
            && family.second == second && family.opposite == opposite)
            return static_cast<int>(i);
    }
    return -1;
}

/**
 * Up to limit members of family, in a fixed order that starts with its prototype. family is
 * one of the families of g; a record whose atoms do not lie on shortest paths of g from its
 * root, as a family's do, gives no members.
 */
inline std::vector<std::vector<int>> relevantCycleFamilyMembers(const MolGraph& g, const RelevantCycleFamily& family,
                                                                std::size_t limit) {
    std::vector<std::vector<int>> members;
    const int root = family.root, first = family.first, second = family.second, opposite = family.opposite;
    for (int atom : {root, first, second})
        if (atom < 0 || atom >= g.n) return members;
    if (limit == 0 || opposite < -1 || opposite >= g.n || first >= second || second > root || family.size < 3)
        return members;
    const auto systems = g.cycleSystems();
    if (systems.ofAtom.empty() || systems.ofAtom[root] < 0) return members;
    detail_rf::FamilyLevels levels;
    levels.search(g, systems.ofAtom, root, true);
    if (!levels.shortest(first) || !levels.shortest(second)) return members;
    const int level = levels.belowRoot[first];
    if (opposite < 0) {
        const auto& row = g.neighbors[first];
        if (family.size != 2 * level + 1 || levels.belowRoot[second] != level
            || std::find(row.begin(), row.end(), second) == row.end()) return members;
    } else {
        const auto& parents = levels.parents[opposite];
        if (!levels.shortest(opposite) || family.size != 2 * levels.belowRoot[opposite]
            || std::find(parents.begin(), parents.end(), first) == parents.end()
            || std::find(parents.begin(), parents.end(), second) == parents.end()) return members;
    }
    // Every shortest path from an atom down its parents to the root, depth first in parent order.
    struct Frame { int atom; size_t next; };
    auto forEachPath = [&](int from, auto&& visit) {
        std::vector<Frame> stack{{from, 0}};
        std::vector<int> path{from};
        while (!stack.empty()) {
            Frame& frame = stack.back();
            if (frame.atom == root) {
                if (!visit(path)) return false;
                stack.pop_back();
                path.pop_back();
                continue;
            }
            const auto& parents = levels.parents[frame.atom];
            if (frame.next == parents.size()) { stack.pop_back(); path.pop_back(); continue; }
            const int next = parents[frame.next++];
            stack.push_back({next, 0});
            path.push_back(next);
        }
        return true;
    };
    forEachPath(first, [&](const std::vector<int>& toFirst) {
        return forEachPath(second, [&](const std::vector<int>& toSecond) {
            std::vector<int> cycle(toFirst.rbegin(), toFirst.rend());
            if (opposite >= 0) cycle.push_back(opposite);
            cycle.insert(cycle.end(), toSecond.begin(), toSecond.end() - 1);
            members.push_back(std::move(cycle));
            return members.size() < limit;
        });
    });
    return members;
}

namespace detail_rf {

/**
 * Whether families prove that cycles lists every relevant cycle of g: the enumeration is
 * complete, the cycles are distinct members of listed families, and they are as many as the
 * families' members.
 */
inline bool relevantCyclesProvedComplete(const MolGraph& g, const RelevantCycleFamilies& families,
                                         const std::vector<std::vector<int>>& cycles) {
    if (families.truncated || std::to_string(cycles.size()) != families.relevantCycleCount) return false;
    std::set<std::vector<int64_t>> distinct;
    for (const auto& cycle : cycles)
        if (!distinct.insert(canonicalCycleKey(cycle)).second || relevantCycleFamilyOf(g, families, cycle) < 0)
            return false;
    return true;
}

inline bool relevantCyclesProvedComplete(const MolGraph& g, const std::vector<std::vector<int>>& cycles) {
    return relevantCyclesProvedComplete(g, computeRelevantCycleFamiliesUncached(g), cycles);
}

} // namespace detail_rf

} // namespace smsd

#endif // SMSD_RELEVANT_CYCLE_FAMILIES_HPP

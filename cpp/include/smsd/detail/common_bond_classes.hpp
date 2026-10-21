/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Common-bond (MCES) search: the work clock and bit helpers shared by its parts, the
 * atom and bond kinds and classes of both molecules with the b-matching groups that
 * bound a mapping, and the colour refinement that labels atoms.
 *
 * Part of common_bond_search.hpp, which describes the search as a whole.
 */
#pragma once

#include "smsd/bitops.hpp"
#include "smsd/detail/union_find.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {
namespace common_bond_detail {

/** The number of radii, from 0, at which atom signatures are compared. */
constexpr int RADII = 4;
/** The most atom pairs a search takes, query atoms times target atoms, as a pair is packed in an int. */
constexpr int64_t MAX_PAIRS = int64_t(1) << 30;
/**
 * The most 64-bit words the compatible bond pairs may take, query bonds times the words of the
 * target bonds: 32 MiB, and about as much again while the bond kinds are worked out. Larger
 * graphs are declined, and their MCES search finds nothing.
 */
constexpr int64_t MAX_BOND_PAIR_WORDS = int64_t(1) << 22;
/** The odd constant that mixes the signature of an atom with those of its neighbours. */
constexpr uint64_t GOLDEN_GAMMA = 0x9e3779b97f4a7c15ULL;

/** Ends a search: the time limit was reached or the search was cancelled. */
struct Stop {};

/**
 * Counts the work of the search in units, about one for each step of an inner loop, and asks
 * whether to stop every POLL_INTERVAL units, or at once with poll(). The search ends with Stop
 * when the answer is yes, so every loop of the search ticks the clock.
 */
class Clock {
    static constexpr int64_t POLL_INTERVAL = 1024;
    std::function<bool()> stop_;
    int64_t work_ = 0;
public:
    explicit Clock(std::function<bool()> stop) : stop_(std::move(stop)) {}
    /** Counts units of work, and checks whether to stop once the interval is reached. */
    void tick(int64_t units = 1) {
        work_ += units;
        if (work_ >= POLL_INTERVAL) {
            work_ = 0;
            poll();
        }
    }
    /** Checks whether to stop now, and throws Stop if so. */
    void poll() {
        if (stop_()) throw Stop{};
    }
};

/** The number of 64-bit words that hold a number of bits. */
inline int words(int bits) { return (bits + 63) >> 6; }

/**
 * Packs a number of atoms and a number of bonds into a score, bonds << 32 | atoms, which
 * compares by bonds first. Both counts are non-negative, so scores compare as integers.
 */
inline int64_t score(int atoms, int bonds) {
    return (int64_t(bonds) << 32) | int64_t(uint32_t(atoms));
}
/** The atoms of a score. */
inline int atomsOf(int64_t value) { return static_cast<int>(static_cast<uint32_t>(value)); }
/** The bonds of a score. */
inline int bondsOf(int64_t value) { return static_cast<int>(static_cast<uint64_t>(value) >> 32); }

/** Whether bit `bit` of row `row` is set, in rows of `words` words. */
inline bool hasBit(const std::vector<uint64_t>& rows, int words, int row, int bit) {
    return ((rows[size_t(row) * size_t(words) + size_t(bit >> 6)] >> (bit & 63)) & 1U) != 0;
}
/** Sets bit `bit` of row `row`. */
inline void setBit(std::vector<uint64_t>& rows, int words, int row, int bit) {
    rows[size_t(row) * size_t(words) + size_t(bit >> 6)] |= uint64_t(1) << (bit & 63);
}
/** Clears bit `bit` of row `row`. */
inline void clearBit(std::vector<uint64_t>& rows, int words, int row, int bit) {
    rows[size_t(row) * size_t(words) + size_t(bit >> 6)] &= ~(uint64_t(1) << (bit & 63));
}

/** The first set bit of a row of `words` words, at or after bit `from`, or -1. */
inline int nextBit(const uint64_t* bits, int words, int from) {
    int wordIndex = from >> 6;
    if (wordIndex >= words) return -1;
    uint64_t word = bits[wordIndex] & (~uint64_t(0) << (from & 63));
    while (word == 0) {
        if (++wordIndex == words) return -1;
        word = bits[wordIndex];
    }
    return (wordIndex << 6) + ctz64(word);
}

/** The rank of each key among the distinct keys, from 0. */
inline std::vector<int> denseRanks(const std::vector<int64_t>& keys) {
    std::vector<int64_t> sorted(keys);
    std::sort(sorted.begin(), sorted.end());
    sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
    std::vector<int> ranks(keys.size());
    for (size_t i = 0; i < keys.size(); ++i)
        ranks[i] = static_cast<int>(std::lower_bound(sorted.begin(), sorted.end(), keys[i]) - sorted.begin());
    return ranks;
}

/** Mixes the bits of a value with a 64-bit finaliser. */
inline uint64_t mix(uint64_t value) {
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

/** A union-find over `size` items, each in a set of its own. */
inline std::vector<int> singletons(int size) {
    std::vector<int> parent(size);
    for (int i = 0; i < size; ++i) parent[i] = i;
    return parent;
}

/** The inverse of a permutation: the position of each value. */
inline std::vector<int> inverse(const std::vector<int>& permutation) {
    std::vector<int> result(permutation.size());
    for (size_t i = 0; i < permutation.size(); ++i) result[permutation[i]] = static_cast<int>(i);
    return result;
}

/**
 * Sorts the atoms of a seed frontier from `from` up to `size` by descending `rank` and merges
 * them into the atoms before them, which are in that order already. `merged` is scratch as large
 * as the frontier.
 */
inline void mergeByRank(std::vector<int>& frontier, int from, int size, const std::vector<int>& rank,
                        std::vector<int>& merged, Clock& clock) {
    if (from == size) return;
    auto higher = [&rank](int first, int second) { return rank[first] > rank[second]; };
    clock.tick(int64_t(size) + int64_t(size - from) * 16);
    std::sort(frontier.begin() + from, frontier.begin() + size, higher);
    std::merge(frontier.begin(), frontier.begin() + from, frontier.begin() + from, frontier.begin() + size, merged.begin(),
               higher);
    std::copy(merged.begin(), merged.begin() + size, frontier.begin());
}

/**
 * Fine classes of both molecules, atom kinds or fine labels, grouped by the coarse class they
 * fall in; no two of different groups match. A mapping pairs items of matching fine classes,
 * so the items it can pair within a group are at most a maximum bipartite b-matching between
 * the query and target counts of the group over the matching pairs. A group in which every
 * pair matches, as each group does when the relation is an equivalence, gives the fewer items
 * of either side. The flow is found by augmenting paths on the small graph of one group.
 */
class Groups {
public:
    /** The number of groups. */
    int count = 0;
    /** The group of each fine query class and of each fine target class, -1 for none. */
    std::vector<int> queryGroup, targetGroup;
    /** Whether every member pair of each group matches. */
    std::vector<uint8_t> complete;
    /** Whether every group is complete. */
    bool allComplete = true;

    Groups() = default;

    /** Groups the fine classes by the groups given, with the pairs that `matches` accepts. */
    template <class Matches>
    Groups(int groups, std::vector<int> queryGroups, std::vector<int> targetGroups, Matches matches)
        : count(groups), queryGroup(std::move(queryGroups)), targetGroup(std::move(targetGroups)) {
        queryMembers_ = members(queryGroup, count);
        targetMembers_ = members(targetGroup, count);
        match_.assign(count, {});
        complete.assign(count, 0);
        size_t mostQuery = 0, mostTarget = 0, mostPairs = 0;
        for (int group = 0; group < count; ++group) {
            const auto& q = queryMembers_[group];
            const auto& t = targetMembers_[group];
            std::vector<uint8_t> pairs(q.size() * t.size(), 0);
            bool every = true;
            for (size_t i = 0; i < q.size(); ++i)
                for (size_t j = 0; j < t.size(); ++j) {
                    pairs[i * t.size() + j] = matches(q[i], t[j]) ? 1 : 0;
                    every = every && pairs[i * t.size() + j];
                }
            match_[group] = std::move(pairs);
            complete[group] = every ? 1 : 0;
            allComplete = allComplete && every;
            mostQuery = std::max(mostQuery, q.size());
            mostTarget = std::max(mostTarget, t.size());
            mostPairs = std::max(mostPairs, q.size() * t.size());
        }
        supply_.assign(mostQuery, 0);
        demand_.assign(mostTarget, 0);
        flow_.assign(mostPairs, 0);
        previous_.assign(mostQuery + mostTarget, 0);
        queue_.assign(mostQuery + mostTarget, 0);
    }

    /** The most items paired off over all groups, from the counts of each fine query and target class. */
    int64_t bound(const std::vector<int>& queryCounts, const std::vector<int>& targetCounts) {
        int64_t total = 0;
        for (int group = 0; group < count; ++group) total += bound(group, queryCounts, targetCounts);
        return total;
    }

    /**
     * The most items of a group paired off, each query item with a target item of a matching
     * class, from the counts of each fine query and target class.
     */
    int bound(int group, const std::vector<int>& queryCounts, const std::vector<int>& targetCounts) {
        const auto& q = queryMembers_[group];
        const auto& t = targetMembers_[group];
        const int queryCount = static_cast<int>(q.size()), targetCount = static_cast<int>(t.size());
        int total = 0;
        if (complete[group]) {
            int queryItems = 0, targetItems = 0;
            for (int member : q) queryItems += queryCounts[member];
            for (int member : t) targetItems += targetCounts[member];
            return std::min(queryItems, targetItems);
        }
        const auto& pairs = match_[group];
        for (int i = 0; i < queryCount; ++i) supply_[i] = queryCounts[q[i]];
        for (int j = 0; j < targetCount; ++j) demand_[j] = targetCounts[t[j]];
        // A greedy flow to start, then augmenting paths, which may take back flow placed before.
        for (int i = 0; i < queryCount; ++i)
            for (int j = 0; j < targetCount; ++j) {
                const int amount = pairs[size_t(i) * targetCount + j] ? std::min(supply_[i], demand_[j]) : 0;
                flow_[size_t(i) * targetCount + j] = amount;
                supply_[i] -= amount;
                demand_[j] -= amount;
                total += amount;
            }
        // Query classes are nodes 0 to queryCount - 1, and target classes follow them.
        for (;;) {
            int tail = 0, end = -1;
            std::fill(previous_.begin(), previous_.begin() + queryCount + targetCount, -2);
            for (int i = 0; i < queryCount; ++i)
                if (supply_[i] > 0) {
                    previous_[i] = -1;
                    queue_[tail++] = i;
                }
            for (int head = 0; head < tail && end < 0; ++head) {
                const int node = queue_[head];
                if (node < queryCount) {
                    for (int j = 0; j < targetCount && end < 0; ++j) {
                        if (!pairs[size_t(node) * targetCount + j] || previous_[queryCount + j] != -2) continue;
                        previous_[queryCount + j] = node;
                        if (demand_[j] > 0) end = queryCount + j;
                        else queue_[tail++] = queryCount + j;
                    }
                } else {
                    for (int i = 0; i < queryCount; ++i)
                        if (flow_[size_t(i) * targetCount + node - queryCount] > 0 && previous_[i] == -2) {
                            previous_[i] = node;
                            queue_[tail++] = i;
                        }
                }
            }
            if (end < 0) return total;
            int amount = demand_[end - queryCount], node = end;
            for (; previous_[node] != -1; node = previous_[node])
                if (node < queryCount)
                    amount = std::min(amount, flow_[size_t(node) * targetCount + previous_[node] - queryCount]);
            amount = std::min(amount, supply_[node]);
            supply_[node] -= amount;
            demand_[end - queryCount] -= amount;
            for (node = end; previous_[node] != -1; node = previous_[node]) {
                if (node < queryCount) flow_[size_t(node) * targetCount + previous_[node] - queryCount] -= amount;
                else flow_[size_t(previous_[node]) * targetCount + node - queryCount] += amount;
            }
            total += amount;
        }
    }

private:
    /** The fine query and target classes of each group. */
    std::vector<std::vector<int>> queryMembers_, targetMembers_;
    /** Whether each query member of a group matches each target member, a row for each query member. */
    std::vector<std::vector<uint8_t>> match_;
    /** Scratch of bound(): what is left of each count, the flow on each pair, and the search for a path. */
    std::vector<int> supply_, demand_, flow_, previous_, queue_;

    static std::vector<std::vector<int>> members(const std::vector<int>& group, int count) {
        std::vector<std::vector<int>> result(count);
        for (size_t item = 0; item < group.size(); ++item)
            if (group[item] >= 0) result[group[item]].push_back(static_cast<int>(item));
        return result;
    }
};

/** A row of 64-bit words within a larger array, compared by its contents. */
struct RowView {
    const uint64_t* bits;
    int length;
};
struct RowHash {
    size_t operator()(const RowView& row) const {
        uint64_t value = 1;
        for (int i = 0; i < row.length; ++i) value = mix(value * 31U + row.bits[i]);
        return static_cast<size_t>(value);
    }
};
struct RowEqual {
    bool operator()(const RowView& first, const RowView& second) const {
        return first.length == second.length && std::equal(first.bits, first.bits + first.length, second.bits);
    }
};

/** Numbers `count` rows of `words` words so that equal rows share a number, in order of first appearance. */
inline std::vector<int> equalRows(const std::vector<uint64_t>& rows, int count, int words, Clock& clock) {
    std::unordered_map<RowView, int, RowHash, RowEqual> seen;
    seen.reserve(size_t(count));
    std::vector<int> number(count);
    for (int row = 0; row < count; ++row) {
        clock.tick(1 + words);
        const RowView view{rows.data() + size_t(row) * size_t(words), words};
        const int next = static_cast<int>(seen.size());
        number[row] = seen.emplace(view, next).first->second;
    }
    return number;
}

/** The columns of `rowCount` rows of `columnCount` bits as rows of `rowCount` bits. */
inline std::vector<uint64_t> transpose(const std::vector<uint64_t>& rows, int rowCount, int columnCount, Clock& clock) {
    const int rowWords = words(columnCount), columnWords = words(rowCount);
    std::vector<uint64_t> columns(size_t(columnCount) * size_t(columnWords), 0);
    for (int row = 0; row < rowCount; ++row) {
        clock.tick(1 + rowWords);
        const uint64_t* bits = rows.data() + size_t(row) * size_t(rowWords);
        for (int column = nextBit(bits, rowWords, 0); column >= 0; column = nextBit(bits, rowWords, column + 1))
            setBit(columns, columnWords, column, row);
    }
    return columns;
}

/**
 * The kinds of the items, atoms or bonds, of both molecules under a pair relation given as rows
 * of target items for each query item. Query items share a kind when their rows are equal and
 * target items when their columns are, so two items match exactly when their kinds do. The
 * classes are the connected parts of the relation between the kinds; a kind that matches
 * nothing has none, -1. Whether two kinds match is read from the rows, at an item of each kind,
 * so the kinds add no table of their own; the rows must outlive them.
 */
struct Kinds {
    std::vector<int> queryKind, targetKind;
    int queryKinds = 0, targetKinds = 0, classes = 0;
    std::vector<int> queryClass, targetClass;

    Kinds(const std::vector<uint64_t>& rows, int queryItems, int targetItems, Clock& clock)
        : rows_(rows), rowWords_(words(targetItems)) {
        queryKind = equalRows(rows, queryItems, rowWords_, clock);
        targetKind = equalRows(transpose(rows, queryItems, targetItems, clock), targetItems, words(queryItems), clock);
        queryKinds = countOf(queryKind);
        targetKinds = countOf(targetKind);
        queryExample_ = examples(queryKind, queryKinds);
        targetExample_ = examples(targetKind, targetKinds);
        std::vector<int> sets = singletons(queryKinds + targetKinds);
        for (int query = 0; query < queryKinds; ++query) {
            clock.tick(1 + targetKinds);
            for (int target = 0; target < targetKinds; ++target)
                if (matches(query, target)) detail::unionFindJoin(sets, query, queryKinds + target);
        }
        std::vector<int> members(size_t(queryKinds) + targetKinds, 0), setClass(size_t(queryKinds) + targetKinds, -1);
        for (int kind = 0; kind < queryKinds + targetKinds; ++kind) ++members[detail::unionFindRoot(sets, kind)];
        int next = 0;
        for (int kind = 0; kind < queryKinds + targetKinds; ++kind) setClass[kind] = members[kind] > 1 ? next++ : -1;
        classes = next;
        queryClass.assign(queryKinds, -1);
        targetClass.assign(targetKinds, -1);
        for (int kind = 0; kind < queryKinds; ++kind) queryClass[kind] = setClass[detail::unionFindRoot(sets, kind)];
        for (int kind = 0; kind < targetKinds; ++kind)
            targetClass[kind] = setClass[detail::unionFindRoot(sets, queryKinds + kind)];
    }

    /** Whether a query kind matches a target kind. */
    bool matches(int query, int target) const {
        return hasBit(rows_, rowWords_, queryExample_[query], targetExample_[target]);
    }

    /** The kinds grouped by class. */
    Groups groups() const {
        return Groups(classes, queryClass, targetClass, [this](int query, int target) { return matches(query, target); });
    }

private:
    const std::vector<uint64_t>& rows_;
    int rowWords_;
    /** The first item of each query and of each target kind. */
    std::vector<int> queryExample_, targetExample_;

    static int countOf(const std::vector<int>& kinds) {
        int result = 0;
        for (int kind : kinds) result = std::max(result, kind + 1);
        return result;
    }
    static std::vector<int> examples(const std::vector<int>& kinds, int count) {
        std::vector<int> example(count, 0);
        for (int item = static_cast<int>(kinds.size()) - 1; item >= 0; --item) example[kinds[item]] = item;
        return example;
    }
};

/**
 * Colour refinement from dense colours. Cells are split by the sorted ranks of the neighbours
 * until none splits. With individualisation, the atom of highest index in the lowest shared
 * cell is then put before its peers and the cells split again. This goes on until each atom has
 * a rank of its own, so that of atoms alike the one of lowest index ranks highest. A rank is
 * where the cell starts. A round sorts only the neighbours of atoms that changed cell, and the
 * largest part of a split cell keeps its name, so a long chain does not cost a full sort each
 * round. The atoms are kept in rank order, and each cell is a run of places in that order.
 */
class Refinement {
    const std::vector<int>& start_;
    const std::vector<int>& neighbour_;
    Clock& clock_;
    int atomCount_;
    std::vector<int> order_, place_, cell_, cellStart_, cellSize_;
    int cells_ = 0;
    /** The sorted ranks of the neighbours of each atom, at its places in the neighbour list. */
    std::vector<int> ranks_;
    /** The atoms to sort in a round, and the last round each was listed for. */
    std::vector<int> toSort_, listed_;
    int toSortCount_ = 0;
    /** The atoms that changed cell in a round. */
    std::vector<int> moved_;
    int movedCount_ = 0;
    int round_ = 1;
    std::vector<int> sorting_;
    std::vector<int64_t> byCell_;

    /** Orders atoms by the sorted ranks of their neighbours. */
    int compareNeighbours(int left, int right) const {
        const int leftDegree = start_[left + 1] - start_[left], rightDegree = start_[right + 1] - start_[right];
        if (leftDegree != rightDegree) return leftDegree < rightDegree ? -1 : 1;
        for (int leftPlace = start_[left], rightPlace = start_[right]; leftPlace < start_[left + 1]; ++leftPlace, ++rightPlace)
            if (ranks_[leftPlace] != ranks_[rightPlace]) return ranks_[leftPlace] < ranks_[rightPlace] ? -1 : 1;
        return 0;
    }

public:
    /** Starts with a cell for each colour, and every atom listed to sort. */
    Refinement(const std::vector<int>& start, const std::vector<int>& neighbour, const std::vector<int>& colour, Clock& clock)
        : start_(start), neighbour_(neighbour), clock_(clock), atomCount_(static_cast<int>(colour.size())),
          order_(atomCount_), place_(atomCount_), cell_(atomCount_), cellStart_(atomCount_, 0), cellSize_(atomCount_, 0),
          ranks_(neighbour.size()), toSort_(atomCount_), listed_(atomCount_, 1), moved_(atomCount_),
          sorting_(atomCount_), byCell_(atomCount_) {
        for (int atomColour : colour) ++cellSize_[atomColour];
        while (cells_ < atomCount_ && cellSize_[cells_] > 0) {
            cellStart_[cells_] = cells_ == 0 ? 0 : cellStart_[cells_ - 1] + cellSize_[cells_ - 1];
            ++cells_;
        }
        std::vector<int> filled(atomCount_, 0);
        for (int atom = 0; atom < atomCount_; ++atom) {
            cell_[atom] = colour[atom];
            place_[atom] = cellStart_[cell_[atom]] + filled[cell_[atom]]++;
            order_[place_[atom]] = atom;
            toSort_[atom] = atom;
        }
        toSortCount_ = atomCount_;
    }

    /** Refines the cells, and gives the rank of each atom. */
    std::vector<int> ranks(bool individualise) {
        for (int firstShared = 0;; ++round_) {
            movedCount_ = 0;
            const int shared = sortByCell();
            // Every rank is read before any cell splits.
            for (int from = 0, to; from < shared; from = to) {
                to = cellEnd(from, shared);
                readRanks(from, to);
            }
            for (int from = 0, to; from < shared; from = to) {
                to = cellEnd(from, shared);
                splitCell(from, to);
            }
            if (movedCount_ == 0) {
                while (firstShared < atomCount_ && cellSize_[cell_[order_[firstShared]]] == 1) ++firstShared;
                if (!individualise || firstShared == atomCount_) break;
                individualiseAt(firstShared);
            }
            listNeighbours();
        }
        std::vector<int> rank(atomCount_);
        for (int atom = 0; atom < atomCount_; ++atom) rank[atom] = cellStart_[cell_[atom]];
        return rank;
    }

private:
    /** Keeps the atoms to sort whose cell has other atoms, in order of cell and then of atom, and gives their number. */
    int sortByCell() {
        int shared = 0;
        for (int i = 0; i < toSortCount_; ++i) {
            const int atom = toSort_[i];
            if (cellSize_[cell_[atom]] > 1) byCell_[shared++] = (int64_t(cell_[atom]) << 32) | atom;
        }
        std::sort(byCell_.begin(), byCell_.begin() + shared);
        for (int i = 0; i < shared; ++i) toSort_[i] = static_cast<int>(byCell_[i] & 0xffffffff);
        return shared;
    }

    /** The end of the atoms to sort from place `from` that are in its cell. */
    int cellEnd(int from, int shared) const {
        int to = from + 1;
        while (to < shared && cell_[toSort_[to]] == cell_[toSort_[from]]) ++to;
        return to;
    }

    /**
     * Sorts the neighbour ranks of the atoms to sort from `from` to `to - 1`, all of one cell,
     * and of one other atom of the cell if there is one. That atom stands for the rest, which
     * have the same ranks.
     */
    void readRanks(int from, int to) {
        for (int i = from; i < to; ++i) sortNeighbourRanks(toSort_[i]);
        const int standIn = unlisted(cell_[toSort_[from]]);
        if (standIn >= 0) sortNeighbourRanks(standIn);
    }

    void sortNeighbourRanks(int atom) {
        clock_.tick(1 + start_[atom + 1] - start_[atom]);
        for (int place = start_[atom]; place < start_[atom + 1]; ++place) ranks_[place] = cellStart_[cell_[neighbour_[place]]];
        std::sort(ranks_.begin() + start_[atom], ranks_.begin() + start_[atom + 1]);
    }

    /** An atom of a cell that is not listed for this round, or -1 if there is none. */
    int unlisted(int cellName) const {
        for (int position = cellStart_[cellName]; position < cellStart_[cellName] + cellSize_[cellName]; ++position)
            if (listed_[order_[position]] != round_) return order_[position];
        return -1;
    }

    /**
     * Splits the cell of the atoms to sort from `from` to `to - 1` into runs of equal neighbour
     * ranks. The unlisted atoms of the cell form one run, that of their stand-in. Each run
     * becomes a cell, and the largest keeps the name of the cell.
     */
    void splitCell(int from, int to) {
        const int cellName = cell_[toSort_[from]], standIn = unlisted(cellName);
        const int unlistedCount = cellSize_[cellName] - (to - from);
        int count = 0;
        for (int i = from; i < to; ++i) sorting_[count++] = toSort_[i];
        if (standIn >= 0) sorting_[count++] = standIn;
        std::stable_sort(sorting_.begin(), sorting_.begin() + count,
                         [this](int left, int right) { return compareNeighbours(left, right) < 0; });
        clock_.tick(count);
        if (runEnd(0, count) == count) return;
        int largestRun = 0, largestSize = 0, standInRun = -1;
        for (int runFrom = 0, runTo; runFrom < count; runFrom = runTo) {
            runTo = runEnd(runFrom, count);
            const int size = runSize(runFrom, runTo, standIn, unlistedCount);
            if (size > largestSize) {
                largestRun = runFrom;
                largestSize = size;
            }
            for (int sortedPlace = runFrom; sortedPlace < runTo; ++sortedPlace)
                if (sorting_[sortedPlace] == standIn) standInRun = runFrom;
        }
        // The listed atoms move to their places in run order, around the run of the stand-in.
        // That run then holds the unlisted atoms.
        int renamed = -1;
        for (int runFrom = 0, runTo, position = cellStart_[cellName]; runFrom < count; runFrom = runTo) {
            runTo = runEnd(runFrom, count);
            const int size = runSize(runFrom, runTo, standIn, unlistedCount);
            const int runCell = runFrom == largestRun ? cellName : cells_++;
            cellStart_[runCell] = position;
            cellSize_[runCell] = size;
            if (runFrom == standInRun) {
                if (runCell != cellName) renamed = runCell;
            } else {
                for (int sortedPlace = runFrom; sortedPlace < runTo; ++sortedPlace) {
                    moveTo(sorting_[sortedPlace], position + sortedPlace - runFrom);
                    if (runCell != cellName) changeCell(sorting_[sortedPlace], runCell);
                }
            }
            position += size;
        }
        if (renamed >= 0)
            for (int position = cellStart_[renamed]; position < cellStart_[renamed] + cellSize_[renamed]; ++position)
                changeCell(order_[position], renamed);
    }

    int runEnd(int runFrom, int count) const {
        int runTo = runFrom + 1;
        while (runTo < count && compareNeighbours(sorting_[runFrom], sorting_[runTo]) == 0) ++runTo;
        return runTo;
    }

    /** The atoms in a run of sorted atoms; the stand-in counts for every unlisted atom. */
    int runSize(int runFrom, int runTo, int standIn, int unlistedCount) const {
        int size = 0;
        for (int sortedPlace = runFrom; sortedPlace < runTo; ++sortedPlace)
            size += sorting_[sortedPlace] == standIn ? unlistedCount : 1;
        return size;
    }

    /** Puts the atom of highest index in the shared cell that starts at place `first` in a cell of its own, before its peers. */
    void individualiseAt(int first) {
        const int cellName = cell_[order_[first]];
        int pick = first;
        clock_.tick(cellSize_[cellName]);
        for (int position = first + 1; position < first + cellSize_[cellName]; ++position)
            if (order_[position] > order_[pick]) pick = position;
        const int atom = order_[pick];
        moveTo(atom, first);
        cellStart_[cells_] = first;
        cellSize_[cells_] = 1;
        ++cellStart_[cellName];
        --cellSize_[cellName];
        changeCell(atom, cells_++);
    }

    /** Puts an atom at a place, and the atom that was there at its old place. */
    void moveTo(int atom, int newPlace) {
        const int displaced = order_[newPlace];
        order_[place_[atom]] = displaced;
        place_[displaced] = place_[atom];
        order_[newPlace] = atom;
        place_[atom] = newPlace;
    }

    void changeCell(int atom, int cellName) {
        cell_[atom] = cellName;
        moved_[movedCount_++] = atom;
    }

    /** Lists the neighbours of the atoms that changed cell, to sort in the next round. */
    void listNeighbours() {
        toSortCount_ = 0;
        for (int movedIndex = 0; movedIndex < movedCount_; ++movedIndex) {
            const int atom = moved_[movedIndex];
            clock_.tick(1 + start_[atom + 1] - start_[atom]);
            for (int place = start_[atom]; place < start_[atom + 1]; ++place)
                if (listed_[neighbour_[place]] != round_ + 1) {
                    listed_[neighbour_[place]] = round_ + 1;
                    toSort_[toSortCount_++] = neighbour_[place];
                }
        }
    }
};

/** Colour refinement of a graph, see Refinement. */
inline std::vector<int> refine(const std::vector<int>& start, const std::vector<int>& neighbour,
                               const std::vector<int>& colour, bool individualise, Clock& clock) {
    return Refinement(start, neighbour, colour, clock).ranks(individualise);
}

} // namespace common_bond_detail
} // namespace smsd

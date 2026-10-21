/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Common-bond (MCES) search: the branch and bound over atom pairs (CommonBondSearch),
 * with its greedy seeds, bounds, scouting and records.
 *
 * Part of common_bond_search.hpp, which describes the search as a whole.
 */
#pragma once

#include "smsd/detail/common_bond_conditions.hpp"
#include "smsd/detail/common_bond_graphs.hpp"
#include "smsd/detail/union_find.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace smsd {
namespace common_bond_detail {

/**
 * Whether each of the first `count` rows can have a column of its own. Each row is `words`
 * words of `rows`, and its columns are their bits. A row takes a free column, or else follows
 * an augmenting path. There is no recursion, as there may be a row for each query atom.
 */
inline bool hasDistinctTargets(const std::vector<uint64_t>& rows, int count, int words, Clock& clock) {
    const int columns = words << 6;
    std::vector<int> owner(columns, -1), visited(columns, 0), path(count), resume(count);
    std::vector<uint64_t> freeColumns(words, ~uint64_t(0));
    for (int newRow = 0; newRow < count; ++newRow) {
        clock.tick(words);
        int column = -1, step = 0;
        for (int word = 0; word < words && column < 0; ++word) {
            const uint64_t freeBits = rows[size_t(newRow) * words + word] & freeColumns[word];
            if (freeBits != 0) column = (word << 6) + ctz64(freeBits);
        }
        path[0] = newRow;
        resume[0] = 0;
        while (column < 0) {
            clock.tick();
            const int row = path[step];
            const uint64_t* bits = rows.data() + size_t(row) * words;
            column = nextBit(bits, words, resume[step]);
            while (column >= 0 && visited[column] == newRow + 1) {
                clock.tick();
                column = nextBit(bits, words, column + 1);
            }
            if (column < 0) {
                if (--step < 0) return false;
                continue;
            }
            resume[step] = column + 1;
            visited[column] = newRow + 1;
            if (owner[column] >= 0) {
                path[++step] = owner[column];
                resume[step] = 0;
                column = -1;
            }
        }
        // Each row on the path takes the column it reached. The last row takes a free column.
        freeColumns[column >> 6] &= ~(uint64_t(1) << (column & 63));
        owner[column] = path[step];
        for (int earlier = step - 1; earlier >= 0; --earlier) owner[resume[earlier] - 1] = path[earlier];
    }
    return true;
}

/** The number of common bonds that `record`, of a query of `atoms` atoms, leaves out; see CommonBondConditions::select(). */
inline int leftOut(const std::vector<int>& record, int atoms) {
    return static_cast<int>(record.size()) > atoms ? record[atoms] : 0;
}

/**
 * When the last search of a CommonBondSearch first scouts, see CommonBondSearch::scout(), and how
 * far each scouting goes. Searches use these defaults; tests lower them so that small searches
 * also scout and skip refuted root pairs.
 */
struct ScoutLimits {
    /** The pairs the last search maps before it first scouts, and the factor by which the gap to the next scouting grows. */
    int64_t after = int64_t(1) << 15, growth = 4;
    /** The most pairs mapped in one scouting, and below one root pair. */
    int64_t budget = int64_t(1) << 15, root = int64_t(1) << 10;
};

/** The MCES search on the int graphs of a CommonBondGraphs; see the description at the top of this file. */
class CommonBondSearch {
    /** The most greedy seeds grown. */
    static constexpr int SEEDS = 64;
    /** The most free targets a seed looks at next to a hub. */
    static constexpr int HUB_TARGETS = 64;
    /** The levels of likeness of two atoms, from no radius alike to all of them. */
    static constexpr int LIKENESS_LEVELS = RADII + 1;
    /** The fewest common bonds of a fragment. */
    static constexpr int MIN_BONDS = 1;
    /** The most pairs the levels of common bonds map between them, see collectByBonds(). */
    static constexpr int64_t LEVEL_BUDGET = int64_t(1) << 13;
    /** The most counts in a histogram of the bond-label bound, one for each molecule. */
    static constexpr int64_t MAX_HISTOGRAM = int64_t(1) << 24;
    /** The most atom pairs whose root bounds are kept, see rootRefuted(). */
    static constexpr int64_t MAX_ROOT_BOUNDS = int64_t(1) << 20;

    /**
     * A mapped pair of the branch and bound, one for each depth, with what it restores when
     * taken off and what it has left to try. The levels are made once and reused.
     */
    struct Level {
        int queryAtom = 0, targetAtom = 0;
        /** The common bonds and the excluded pairs before the pair was mapped, put back when it is taken off. */
        int commonBondsBefore = 0, excludedBefore = 0;
        /** The query and target bonds with no end mapped, once the pair is mapped. */
        int freeQueryBonds = 0, freeTargetBonds = 0;
        /** The level that starts the open fragment. Closing the fragment sets this to the next level. */
        int fragmentStart = 0;
        /** A score that no mapping below the level beats. */
        int64_t bound = 0;
        /** The candidates whose targets the level tries, in order, and the next to try. */
        std::vector<int> candidates;
        int nextCandidate = 0;
        /** The candidate being tried, its targets in order, and the next target to try. */
        int triedAtom = 0;
        std::vector<int> targets;
        int nextTarget = 0;
        /** Once the fragment is closed, the roots of new fragments offered are ranked below this. */
        int rootRank = 0;
        /** Whether the pair is the root of a fragment. */
        bool startsFragment = false;
    };

    /**
     * The open pairs reached from the open fragment of a state through compatible bonds. The
     * bonds that link them are counted once on each side.
     */
    struct Reach {
        /** The reached pairs as bits, rows by query atom. */
        std::vector<uint64_t> seen;
        /** The reached pairs in the order they were reached. */
        std::vector<int> pairs;
        /** The number of pairs reached from the fragment itself, which come first and give the candidates. */
        int firstLayer = 0;
        /** The stamp of each query and target bond counted as a link. */
        std::vector<int> queryBondStamp, targetBondStamp;
        int queryLinks = 0, targetLinks = 0;
        int stamp = 0;
    };

    /** The labels of one molecule in the bond-label bound, and the open atoms of a state with their counts. */
    struct Labels {
        /** The class of each atom, and the slot at each place of the neighbour lists. */
        std::vector<int> symbol, slot;
        /** The open atoms by slot and by the number of open bonds of that slot. */
        std::vector<int> histogram;
        /** The open atoms by class that may get a new bond. */
        std::vector<int> bySymbol;
        /** The open atoms by kind, and the bonds between open atoms by fine label, where these are counted. */
        std::vector<int> kindCount, fineCount, countedKinds, countedFines;
        std::vector<int> openAtoms;
        std::vector<uint8_t> isOpen;
        std::vector<uint32_t> reached;
        int openCount = 0, kinds = 0, fines = 0;
        uint32_t stamp = 0;

        void init(const Side& side, std::vector<int> slots, size_t histogramSize, int symbols, int atomKinds, int fineLabels) {
            symbol = side.atomClass;
            slot = std::move(slots);
            histogram.assign(histogramSize, 0);
            bySymbol.assign(symbols, 0);
            kindCount.assign(atomKinds, 0);
            fineCount.assign(fineLabels, 0);
            countedKinds.assign(atomKinds, 0);
            countedFines.assign(fineLabels, 0);
            openAtoms.assign(symbol.size(), 0);
            isOpen.assign(symbol.size(), 0);
            reached.assign(symbol.size(), 0);
        }
        void markOpen(int atom) {
            isOpen[atom] = 1;
            openAtoms[openCount++] = atom;
        }
        void countKind(int kind) {
            if (kindCount[kind]++ == 0) countedKinds[kinds++] = kind;
        }
        void countFine(int fine) {
            if (fineCount[fine]++ == 0) countedFines[fines++] = fine;
        }
        /** Clears the counts by kind and fine label. */
        void clearCounts() {
            while (kinds > 0) kindCount[countedKinds[--kinds]] = 0;
            while (fines > 0) fineCount[countedFines[--fines]] = 0;
        }
        /** Closes every open atom and clears the counts by kind and fine label. */
        void clearOpen() {
            for (int openIndex = 0; openIndex < openCount; ++openIndex) isOpen[openAtoms[openIndex]] = 0;
            openCount = 0;
            clearCounts();
        }
        /** Starts a new reach, clearing every mark when the stamp wraps. */
        uint32_t nextStamp() {
            if (++stamp == 0) {
                std::fill(reached.begin(), reached.end(), 0U);
                stamp = 1;
            }
            return stamp;
        }
    };

    class Embedding;

    CommonBondGraphs& graph_;
    const Side& query_;
    const Side& target_;
    Clock& clock_;
    const bool connected_;
    const int queryCount_, targetCount_, words_;
    /** The compatible pairs that are not excluded: a row for each query atom, with a bit for each target atom. */
    std::vector<uint64_t> allowed_;
    /** The mapping both ways, -1 where unmapped. */
    std::vector<int> queryToTarget_, targetToQuery_;
    /** The mapped pairs in the order they were mapped, the first depth_ of them. */
    std::vector<Level> levels_;
    int depth_ = 0, commonBonds_ = 0;
    /** The pairs excluded so far, put back as the levels that tried them go. */
    std::vector<int> excluded_;
    /** The atoms by rank, lowest first, and the rank of each atom. */
    std::vector<int> queryOrder_, queryRank_, targetOrder_, targetRank_;
    /** The cell of each atom: atoms that look alike share a cell. */
    std::vector<int> queryCell_, targetCell_;
    /** The signatures of the atoms at each radius, which tell how alike two atoms are. */
    std::vector<uint64_t> querySignature_, targetSignature_;
    /** The next atom alike to each terminal atom, see twins(). */
    std::vector<int> queryTwin_, targetTwin_;
    bool twinsKnown_ = false;
    /** The taken targets that lead the neighbours of each target atom, and the seed they were counted for. */
    std::vector<int> leadingTaken_, leadingRound_;
    /** The likest level of each query atom with any target, -1 until worked out, see likestLevel(). */
    std::vector<int> likestLevel_;
    /** The best score so far, and a score that no mapping beats. */
    int64_t bestScore_ = 0, upperBound_ = 0;
    /** The mappings held: the embeddings found, or those with the best score so far. */
    std::vector<std::vector<int>> mappings_;
    size_t maxMappings_ = 0;
    /** The seed, held at the front of the mappings from the start until the search finds it again. */
    std::vector<int> heldSeed_;
    bool holdsSeed_ = false;
    /**
     * The mapping a stopped search falls back on while the pass under way holds none: the seed,
     * or the best mapping a level of the search met below its threshold.
     */
    std::vector<int> fallback_;
    bool hasFallback_ = false;
    int64_t fallbackScore_ = 0;
    Reach reach_;
    /** The compatible pairs before any is excluded, see hasPendantPair(). */
    std::vector<uint64_t> compatible_;
    /** The labels of the query and of the target, with the open atoms of the state being bounded. */
    Labels queryLabels_, targetLabels_;
    /** The class at each slot. A label has a slot for each end, the second -1 if both ends share a class. */
    std::vector<int> slotSymbol_;
    /** The most bonds of one slot at an atom, plus one, which is the width of a histogram row. */
    int histogramWidth_ = 0;
    std::vector<int> slotCount_, slotsAt_;
    std::vector<uint64_t> unusedTargets_, reachableTargets_;
    /** The labels and the classes counted for a bound, each listed once while its flag is set. */
    std::vector<int> countedLabels_, countedSymbols_;
    std::vector<uint8_t> labelListed_, symbolListed_;
    /**
     * The union-find of the open atoms of both molecules, the target atoms after the query
     * atoms, the group of each root while its mark is the stamp, the start of each group, the
     * molecules each group has atoms of, and the atoms by group, see groupedBound().
     */
    std::vector<int> partParent_, groupMark_, groupOf_, groupStart_, groupSides_, groupItems_;
    int markStamp_ = 1;
    /** Scratch for sorting candidates and targets. */
    std::vector<int> layer_;
    std::vector<int64_t> keys_;
    /** The pairs mapped and the embedding atoms placed so far. */
    int64_t nodes_ = 0;
    /** The number of nodes at which the pass under way ends, see done(). */
    int64_t nodeLimit_ = INT64_MAX;
    /** The last bound worked out, before it was compared with the best score. */
    int64_t lastBound_ = 0;
    /**
     * The most common bonds below each root pair, by pair(), as the bound of its state showed
     * when it was last pruned, or -1. Empty when the pairs are too many to keep.
     */
    std::vector<int> rootBonds_;
    /** When the last search scouts and how far, see scout(). */
    const ScoutLimits scoutLimits_;
    /** The number of nodes at which the search under way scouts, see scout(), and the gap to the next scouting. */
    int64_t scoutAt_ = INT64_MAX, scoutGap_;
    /** The scoutings so far, and the root pairs skipped as refuted, see rootRefuted(). */
    int64_t scoutings_ = 0, refutedRoots_ = 0;
    /** The root pair the search under way is below, as the rank of its query atom and the place of its target. */
    int rootRank_ = 0, rootPlace_ = 0;
    /** What a mapping must satisfy beyond its common bonds, or null when nothing else is asked. */
    CommonBondConditions* conditions_ = nullptr;
    /** The target that a required mapping fixes for each query atom, -1 for none; empty without one. */
    std::vector<int> requiredTarget_;
    /** The most mappings the search holds, see run(). */
    size_t wanted_ = 1;
    /** The most pairs mapped and embedding atoms placed before the search stops, see countNode(). */
    int64_t nodeCap_ = INT64_MAX;
    /** Whether the search stopped at nodeCap_. */
    bool capReached_ = false;
    /** The lowest rank of a first root: with a required mapping, only that of its highest ranked atom. */
    int lowestRoot_ = 0;
    /** The selections of the mapping being recorded, and of one offered as the fallback. */
    std::vector<std::vector<int>> selections_, offered_;
    /** The mappings held, while several are wanted, see keep(). */
    std::set<std::vector<int>> held_;
    /** The number of fragments of the mapping, the open one included, see nextRoot(). */
    int fragments_ = 0;
    /** Grows a seed of the same molecules the other way round, turned round, or null; see seedAlsoFrom(). */
    std::function<std::vector<int>()> otherSeed_;

public:
    /**
     * A search of `graph`. With `connected` the common bonds form one fragment, otherwise any
     * number of fragments of at least one bond. The search scouts as `limits` say. It holds at
     * most `wanted` mappings that meet `conditions`, unless they are null, that keep the target
     * of each query atom in `requiredTarget`, -1 for none, unless it is empty, and it stops once
     * it would map more than `nodeCap` pairs.
     */
    CommonBondSearch(CommonBondGraphs& graph, Clock& clock, bool connected, const ScoutLimits& limits = ScoutLimits(),
                     CommonBondConditions* conditions = nullptr, std::vector<int> requiredTarget = {},
                     size_t wanted = 1, int64_t nodeCap = INT64_MAX)
        : graph_(graph), query_(graph.query), target_(graph.target), clock_(clock), connected_(connected),
          queryCount_(graph.query.atomCount), targetCount_(graph.target.atomCount), words_(words(targetCount_)),
          allowed_(graph.compatiblePairs()), queryToTarget_(queryCount_, -1), targetToQuery_(targetCount_, -1),
          levels_(size_t(std::min(queryCount_, targetCount_))), scoutLimits_(limits), scoutGap_(limits.after),
          conditions_(conditions), requiredTarget_(std::move(requiredTarget)), wanted_(wanted), nodeCap_(nodeCap) {
        if (!requiredTarget_.empty()) keepRequiredPairs();
        reach_.seen.assign(size_t(queryCount_) * size_t(words_), 0);
        reach_.pairs.reserve(64);
        reach_.queryBondStamp.assign(query_.bondCount, 0);
        reach_.targetBondStamp.assign(target_.bondCount, 0);
        compatible_ = allowed_;
        slotSymbol_.assign(size_t(2) * graph_.labels, 0);
        for (int label = 0; label < graph_.labels; ++label) {
            slotSymbol_[2 * label] = graph_.labelLow[label];
            slotSymbol_[2 * label + 1] = graph_.labelLow[label] == graph_.labelHigh[label] ? -1 : graph_.labelHigh[label];
        }
        histogramWidth_ = histogramWidth(graph);
        const size_t histogramSize = slotSymbol_.size() * size_t(histogramWidth_);
        queryLabels_.init(query_, slots(query_), histogramSize, graph_.atomClasses,
                          static_cast<int>(graph_.atomGroups.queryGroup.size()), static_cast<int>(graph_.labelGroups.queryGroup.size()));
        targetLabels_.init(target_, slots(target_), histogramSize, graph_.atomClasses,
                           static_cast<int>(graph_.atomGroups.targetGroup.size()), static_cast<int>(graph_.labelGroups.targetGroup.size()));
        slotCount_.assign(slotSymbol_.size(), 0);
        slotsAt_.assign(size_t(histogramWidth_), 0);
        if (int64_t(queryCount_) * targetCount_ <= MAX_ROOT_BOUNDS)
            rootBonds_.assign(size_t(queryCount_) * size_t(targetCount_), -1);
        unusedTargets_.assign(size_t(words_), 0);
        reachableTargets_.assign(size_t(words_), 0);
        leadingTaken_.assign(size_t(targetCount_), 0);
        leadingRound_.assign(size_t(targetCount_), 0);
        likestLevel_.assign(size_t(queryCount_), -1);
        countedLabels_.reserve(size_t(graph_.labels));
        countedSymbols_.reserve(size_t(graph_.atomClasses));
        labelListed_.assign(size_t(graph_.labels), 0);
        symbolListed_.assign(size_t(graph_.atomClasses), 0);
        const size_t items = size_t(queryCount_) + size_t(targetCount_);
        partParent_.assign(items, 0);
        groupMark_.assign(items, 0);
        groupOf_.assign(items, 0);
        groupStart_.assign(items + 1, 0);
        groupSides_.assign(items, 0);
        groupItems_.assign(items, 0);
    }

    /**
     * Whether the histograms of the bond-label bound of `graph` stay within MAX_HISTOGRAM
     * counts: a row of the most bonds at an atom, plus one, for each end of each label. A
     * search that needs more is declined.
     */
    static bool fits(const CommonBondGraphs& graph) {
        return int64_t(2) * graph.labels * histogramWidth(graph) <= MAX_HISTOGRAM;
    }

    /** The pairs mapped and the embedding atoms placed so far. */
    int64_t nodes() const { return nodes_; }
    /** Whether the search stopped at its node limit. */
    bool capReached() const { return capReached_; }
    /** The scoutings so far. */
    int64_t scoutings() const { return scoutings_; }
    /** The root pairs skipped so far as refuted. */
    int64_t refutedRoots() const { return refutedRoots_; }

    /**
     * A maximum mapping, as a list of one, or none when no pair is compatible or either
     * molecule has no bond. When the whole query embeds, an embedding is the result. Otherwise
     * the atoms are ranked, a greedy seed gives the score to start from, and the branch and
     * bound finds the mapping.
     */
    std::vector<std::vector<int>> run();

    /**
     * The best greedy seed, as run() grows it before the branch and bound, or empty when there is
     * none. Another search of the same molecules may start from it, see seedAlsoFrom().
     */
    std::vector<int> greedySeed() {
        clock_.poll();
        if (!anyCompatible() || std::min(query_.bondCount, target_.bondCount) < MIN_BONDS) return {};
        rank();
        fallbackScore_ = belowEveryMapping();
        seed(bestPossible());
        return hasFallback_ ? fallback_ : std::vector<int>();
    }

    /**
     * Lets the search ask `grow` for a second seed: a mapping of the same molecules grown by a
     * search the other way round and turned round, or empty. The greedy seed grows from the
     * query, so the two seeds often differ, and the one grown from the other molecule may hold
     * more common bonds. A fragment search asks only once its levels of common bonds run out of
     * budget, as the last search is then long and starts from the better seed, see
     * collectByBonds(). The seed costs no nodes.
     */
    void seedAlsoFrom(std::function<std::vector<int>()> grow) { otherSeed_ = std::move(grow); }

    /**
     * The best mappings found so far, for a search the clock stopped: the mappings held or,
     * while the pass under way holds none, the one it falls back on. Empty when nothing was
     * found. Each is complete, as a mapping is held only once it is recorded whole.
     */
    std::vector<std::vector<int>> found() const {
        if (!mappings_.empty()) return mappings_;
        if (hasFallback_) return {fallback_};
        return {};
    }

private:
    bool has(const std::vector<uint64_t>& rows, int row, int bit) const { return hasBit(rows, words_, row, bit); }
    void set(std::vector<uint64_t>& rows, int row, int bit) const { setBit(rows, words_, row, bit); }
    void clear(std::vector<uint64_t>& rows, int row, int bit) const { clearBit(rows, words_, row, bit); }

    /** A pair as one int, queryAtom * targetCount + targetAtom; pairs sort by query atom, then by target atom. */
    int pair(int queryAtom, int targetAtom) const { return queryAtom * targetCount_ + targetAtom; }
    int pairQuery(int packed) const { return packed / targetCount_; }
    int pairTarget(int packed) const { return packed % targetCount_; }

    /**
     * Keeps for each query atom of the required mapping only the pair with its target, and takes
     * these targets from every other query atom.
     */
    void keepRequiredPairs() {
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
            const int targetAtom = requiredTarget_[queryAtom];
            if (targetAtom < 0) continue;
            const bool compatible = has(allowed_, queryAtom, targetAtom);
            std::fill(allowed_.begin() + size_t(queryAtom) * words_, allowed_.begin() + size_t(queryAtom + 1) * words_, 0);
            if (compatible) set(allowed_, queryAtom, targetAtom);
            for (int other = 0; other < queryCount_; ++other)
                if (other != queryAtom) clear(allowed_, other, targetAtom);
        }
    }

    /** The width of a histogram row: the most bonds at any atom of either molecule, plus one. */
    static int histogramWidth(const CommonBondGraphs& graph) {
        return std::max(graph.query.maxDegree(), graph.target.maxDegree()) + 1;
    }

    /** The most atoms of a mapping of `bonds` common bonds: one fragment of b bonds has at most b + 1 atoms. */
    int mostAtoms(int bonds) const {
        if (bonds <= 0) return 0;
        return connected_ ? bonds + 1 : bonds + bonds / MIN_BONDS;
    }

    /** A score below that of every mapping, as every mapping has a common bond. */
    int64_t belowEveryMapping() const { return score(queryCount_ + 1, MIN_BONDS - 1); }

    bool anyCompatible() const {
        for (uint64_t word : allowed_)
            if (word != 0) return true;
        return false;
    }

    /** Whether every query atom is reached from the first by bonds. */
    bool queryConnected() const {
        if (queryCount_ == 0) return false;
        std::vector<int> queue(queryCount_);
        std::vector<uint8_t> reached(queryCount_, 0);
        int count = 1;
        reached[0] = 1;
        for (int head = 0; head < count; ++head)
            for (int place = query_.start[queue[head]]; place < query_.start[queue[head] + 1]; ++place)
                if (!reached[query_.neighbour[place]]) {
                    reached[query_.neighbour[place]] = 1;
                    queue[count++] = query_.neighbour[place];
                }
        return count == queryCount_;
    }

    /**
     * A score that no mapping beats, when the query does not fit whole. Every cap on the
     * common bonds or on the atoms holds for every mapping, so the two are capped apart: the
     * least of the scores would keep the atoms of the one with the fewest bonds, and lose a
     * lower cap on the atoms.
     */
    int64_t bestPossible() {
        int bondCap = std::min(query_.bondCount, target_.bondCount);
        int atomCap = std::min(queryCount_, targetCount_);
        // Without an embedding no mapping of a connected query has every query bond. A query in
        // parts may still have every bond.
        if (queryConnected()) bondCap = std::min(bondCap, query_.bondCount - 1);
        const int64_t label = graph_.labelBound(connected_);
        bondCap = std::min(bondCap, bondsOf(label));
        atomCap = std::min(atomCap, atomsOf(label));
        // A mapping with at most so many bonds has at most so many atoms.
        return score(std::min(atomCap, mostAtoms(bondCap)), bondCap);
    }

    /**
     * Ranks the atoms and works out their signatures. Query atoms with the fewest candidates,
     * then the most bonds, rank highest; ties are broken by colour refinement, so that the
     * order is the same for the same graph. The search tries the roots from the highest rank
     * down, and below a root maps only atoms of lower rank. Target atoms are ranked by their
     * cell only, in input order within it.
     */
    void rank() {
        querySignature_ = query_.signatures(clock_);
        targetSignature_ = target_.signatures(clock_);
        std::vector<int> queryChoices(queryCount_, 0), targetChoices(targetCount_, 0);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
            const uint64_t* row = allowed_.data() + size_t(queryAtom) * words_;
            for (int targetAtom = nextBit(row, words_, 0); targetAtom >= 0; targetAtom = nextBit(row, words_, targetAtom + 1)) {
                clock_.tick();
                ++queryChoices[queryAtom];
                ++targetChoices[targetAtom];
            }
        }
        const std::vector<int> colours = query_.invariants(queryChoices);
        queryCell_ = refine(query_.start, query_.neighbour, colours, false, clock_);
        const std::vector<int> canon = refine(query_.start, query_.neighbour, colours, true, clock_);
        // Ties among target atoms stay in input order. Breaking them one atom at a time costs a
        // round per atom on a symmetric target, and the order only has to be repeatable.
        targetCell_ = refine(target_.start, target_.neighbour, target_.invariants(targetChoices), false, clock_);
        std::vector<int> next(targetCount_, 0);
        targetRank_.assign(targetCount_, 0);
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom) next[targetCell_[targetAtom]] = targetCell_[targetAtom];
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom) targetRank_[targetAtom] = next[targetCell_[targetAtom]]++;
        targetOrder_ = inverse(targetRank_);
        std::vector<int> byRank(queryCount_);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) byRank[queryAtom] = queryAtom;
        // The most choices come first, then the fewest bonds, then the lowest canonical rank. The
        // atoms of a required mapping come last of all, so that its highest ranked atom roots
        // every mapping that holds it.
        std::stable_sort(byRank.begin(), byRank.end(), [&](int left, int right) {
            if (!requiredTarget_.empty() && (requiredTarget_[left] >= 0) != (requiredTarget_[right] >= 0))
                return requiredTarget_[left] < 0;
            if (queryChoices[left] != queryChoices[right]) return queryChoices[left] > queryChoices[right];
            if (query_.degree(left) != query_.degree(right)) return query_.degree(left) < query_.degree(right);
            return canon[left] < canon[right];
        });
        queryOrder_.assign(queryCount_, 0);
        queryRank_.assign(queryCount_, 0);
        for (int rankIndex = 0; rankIndex < queryCount_; ++rankIndex) {
            queryOrder_[rankIndex] = byRank[rankIndex];
            queryRank_[queryOrder_[rankIndex]] = rankIndex;
        }
        lowestRoot_ = requiredTarget_.empty() ? 0 : queryCount_ - 1;
    }

    /**
     * The best greedy seed, if any, as the fallback. A seed is a connected mapping grown from a
     * root pair, from each of at most SEEDS of them or until one reaches `bound`. Unless the
     * search is connected, it is a set of fragments of at least one common bond, grown one after
     * another on the free atoms. The fallback follows each better seed, so a search stopped while
     * seeding still has the best seed grown by then, or the seed being grown if that is better.
     * It needs rank() first.
     */
    void seed(int64_t bound) {
        int keptBonds = 0;
        try {
            growSeeds(bound, keptBonds);
        } catch (const Stop&) {
            // Each pair added to the seed being grown brought a common bond, so it is a fragment once it has one.
            // Its conditions are not checked once the search is to stop.
            const int64_t current = scoreOf(queryToTarget_);
            if (bondsOf(current) - keptBonds >= MIN_BONDS && current > fallbackScore_ && !conditions_) {
                fallback_ = queryToTarget_;
                hasFallback_ = true;
                fallbackScore_ = current;
            }
            throw;
        }
        std::fill(queryToTarget_.begin(), queryToTarget_.end(), -1);
        std::fill(targetToQuery_.begin(), targetToQuery_.end(), -1);
    }

    /**
     * Makes `mapping`, of score `current`, the fallback if it beats it. Under conditions its best
     * selection that meets them, if any, is the one that has to beat it.
     */
    void offerFallback(const std::vector<int>& mapping, int64_t current) {
        if (current <= fallbackScore_) return;
        if (!conditions_) {
            fallback_ = mapping;
            hasFallback_ = true;
            fallbackScore_ = current;
            return;
        }
        offered_.clear();
        const int64_t valid = conditions_->select(mapping, fallbackScore_ + 1, 1, offered_);
        if (valid >= 0) {
            fallback_ = std::move(offered_.front());
            hasFallback_ = true;
            fallbackScore_ = valid;
        }
    }

    /** The rounds of seed(), with the common bonds of the fragments kept from earlier rounds in `keptBonds`. */
    void growSeeds(int64_t bound, int& keptBonds) {
        std::vector<uint64_t> taken(size_t(queryCount_) * size_t(words_), 0);
        std::vector<int> roots(size_t(LIKENESS_LEVELS) * SEEDS, 0), rootCount(LIKENESS_LEVELS, 0);
        std::vector<int> frontier(queryCount_), merged(queryCount_), mapped(queryCount_), added(queryCount_, 0);
        int64_t bestSeed = belowEveryMapping();
        std::vector<int> bestMapping;
        for (int round = 0, keptAtoms = 0;;) {
            const int64_t before = bestSeed;
            seedRoots(roots, rootCount, taken);
            for (int level = RADII, tried = 0; level >= 0 && bestSeed != bound; --level)
                for (int rootIndex = 0; rootIndex < rootCount[level] && tried < SEEDS && bestSeed != bound; ++rootIndex, ++tried) {
                    const int64_t grown =
                        growSeed(roots[size_t(level) * SEEDS + rootIndex], frontier, merged, mapped, added, ++round);
                    const int seedAtoms = atomsOf(grown), seedBonds = bondsOf(grown);
                    if (seedBonds >= MIN_BONDS && score(keptAtoms + seedAtoms, keptBonds + seedBonds) > bestSeed) {
                        bestSeed = score(keptAtoms + seedAtoms, keptBonds + seedBonds);
                        bestMapping = queryToTarget_;
                        offerFallback(bestMapping, scoreOf(bestMapping));
                    }
                    for (int mappedIndex = 0; mappedIndex < seedAtoms; ++mappedIndex) {
                        targetToQuery_[queryToTarget_[mapped[mappedIndex]]] = -1;
                        queryToTarget_[mapped[mappedIndex]] = -1;
                    }
                }
            if (connected_ || bestSeed == before || bestSeed == bound) break;
            // The fragments found so far stay mapped while the next one grows. It shares no
            // common bond with them, since each grew until no pair next to it was left.
            keptAtoms = atomsOf(bestSeed);
            keptBonds = bondsOf(bestSeed);
            queryToTarget_ = bestMapping;
            for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
                if (queryToTarget_[queryAtom] >= 0) targetToQuery_[queryToTarget_[queryAtom]] = queryAtom;
        }
    }

    /**
     * Picks the seed roots by how alike their atoms are. The roots of each level are at
     * roots[level * SEEDS]. Alike cells grow alike seeds, so there is one root for each pair of
     * cells, marked in `taken`.
     */
    void seedRoots(std::vector<int>& roots, std::vector<int>& rootCount, std::vector<uint64_t>& taken) {
        std::fill(taken.begin(), taken.end(), 0);
        std::fill(rootCount.begin(), rootCount.end(), 0);
        const int likest = RADII;
        // The most alike roots are tried first. Once there are enough of them, the rest are not needed.
        for (int rankIndex = queryCount_ - 1; rankIndex >= 0 && rootCount[likest] < SEEDS; --rankIndex) {
            const int queryAtom = queryOrder_[rankIndex];
            for (int targetPosition = 0; targetPosition < targetCount_ && rootCount[likest] < SEEDS
                     && queryToTarget_[queryAtom] < 0; ++targetPosition) {
                const int targetAtom = targetOrder_[targetPosition];
                clock_.tick();
                if (targetToQuery_[targetAtom] >= 0 || !has(allowed_, queryAtom, targetAtom)) continue;
                const int level = likeness(queryAtom, targetAtom);
                if (rootCount[level] < SEEDS && !has(taken, queryCell_[queryAtom], targetCell_[targetAtom])) {
                    set(taken, queryCell_[queryAtom], targetCell_[targetAtom]);
                    roots[size_t(level) * SEEDS + rootCount[level]++] = pair(queryAtom, targetAtom);
                }
            }
        }
    }

    /**
     * Grows a seed from the root pair, each time by the best frontier pair. The frontier is the
     * free query atoms next to the seed, each listed once by marking it with `round` in
     * `added`, highest rank first; `merged` is scratch for it. The best pair is the same in any
     * order, and in this one bestSeedPair() soon passes over the atoms that cannot beat the pair
     * it holds. The seed is left mapped, with its query atoms in `mapped` in order. The result is
     * its score.
     */
    int64_t growSeed(int root, std::vector<int>& frontier, std::vector<int>& merged, std::vector<int>& mapped,
                     std::vector<int>& added, int round) {
        int count = 0, frontierSize = 0, seedBonds = 0;
        for (int packed = root; packed >= 0;) {
            const int queryAtom = pairQuery(packed), targetAtom = pairTarget(packed);
            queryToTarget_[queryAtom] = targetAtom;
            targetToQuery_[targetAtom] = queryAtom;
            mapped[count++] = queryAtom;
            const int listed = frontierSize;
            for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place) {
                const int queryNeighbour = query_.neighbour[place];
                if (queryToTarget_[queryNeighbour] < 0 && added[queryNeighbour] != round) {
                    added[queryNeighbour] = round;
                    frontier[frontierSize++] = queryNeighbour;
                }
            }
            mergeByRank(frontier, listed, frontierSize, queryRank_, merged, clock_);
            int64_t best = -1;
            int kept = 0;
            for (int frontierIndex = 0; frontierIndex < frontierSize; ++frontierIndex) {
                const int frontierAtom = frontier[frontierIndex];
                // The atom just mapped leaves the frontier.
                if (queryToTarget_[frontierAtom] >= 0) continue;
                frontier[kept++] = frontierAtom;
                best = bestSeedPair(frontierAtom, best, round);
            }
            frontierSize = kept;
            packed = best < 0 ? -1 : static_cast<int>(best & 0xffffffff);
            seedBonds += best < 0 ? 0 : static_cast<int>(best >> 32);
        }
        return score(count, seedBonds);
    }

    /**
     * The best pair of `frontierAtom` with a free target next to the image of one of its mapped
     * neighbours, through a compatible bond, unless the pair `held` beats it: the most alike
     * pair, then the one with the most common bonds. Ties go to the higher ranked query atom,
     * then to the lower ranked target atom. The pair is returned in the low 32 bits, with its
     * common bonds above, or -1 if there is none; `held` is read the same way. Next to a hub,
     * the first few free targets are as good for a seed as all of them. The seed of `round` only
     * takes targets, so the taken targets that lead the list of an image are passed once for the
     * seed, see leadingTaken(), and not again for each frontier atom.
     */
    int64_t bestSeedPair(int frontierAtom, int64_t held, int round) {
        int packed = held < 0 ? -1 : static_cast<int>(held & 0xffffffff);
        int links = held < 0 ? 0 : static_cast<int>(held >> 32);
        int likest = packed < 0 ? -1 : likeness(pairQuery(packed), pairTarget(packed));
        clock_.tick(1 + query_.degree(frontierAtom));
        if (packed >= 0) {
            // No pair of the atom is liker than its likest level, nor has more common bonds than the
            // atom has mapped neighbours, so its pairs are not looked at when none can beat the
            // held one. Its likest level is worked out only when it decides.
            int mappedNeighbours = 0;
            for (int queryPlace = query_.start[frontierAtom]; queryPlace < query_.start[frontierAtom + 1]; ++queryPlace)
                mappedNeighbours += queryToTarget_[query_.neighbour[queryPlace]] >= 0;
            const bool fewerBonds = mappedNeighbours < links
                || (mappedNeighbours == links && queryRank_[frontierAtom] < queryRank_[pairQuery(packed)]);
            if (fewerBonds ? likest == RADII || likestLevel(frontierAtom) <= likest
                           : likest > 0 && likestLevel(frontierAtom) < likest)
                return held;
        }
        for (int queryPlace = query_.start[frontierAtom]; queryPlace < query_.start[frontierAtom + 1]; ++queryPlace) {
            const int neighbourImage = queryToTarget_[query_.neighbour[queryPlace]];
            if (neighbourImage < 0) continue;
            for (int targetPlace = target_.start[neighbourImage] + leadingTaken(neighbourImage, round), freeCount = 0;
                 targetPlace < target_.start[neighbourImage + 1] && freeCount < HUB_TARGETS; ++targetPlace) {
                const int candidate = target_.neighbour[targetPlace];
                clock_.tick();
                if (targetToQuery_[candidate] >= 0 || !has(allowed_, frontierAtom, candidate)
                    || !graph_.bondsMatch(query_.neighbourBond[queryPlace], target_.neighbourBond[targetPlace])) continue;
                ++freeCount;
                const int level = likeness(frontierAtom, candidate);
                if (level < likest) continue;
                const int newBonds = common(frontierAtom, candidate);
                if (level > likest || newBonds > links
                    || (newBonds == links
                        && (queryRank_[frontierAtom] > queryRank_[pairQuery(packed)]
                            || (frontierAtom == pairQuery(packed) && targetRank_[candidate] < targetRank_[pairTarget(packed)])))) {
                    packed = pair(frontierAtom, candidate);
                    likest = level;
                    links = newBonds;
                }
            }
        }
        return packed < 0 ? -1 : (int64_t(links) << 32) | int64_t(packed);
    }

    /** The likest level of `queryAtom` with any target it is compatible with, worked out when first needed. */
    int likestLevel(int queryAtom) {
        if (likestLevel_[queryAtom] < 0) {
            const uint64_t* row = compatible_.data() + size_t(queryAtom) * words_;
            int level = 0;
            for (int targetAtom = common_bond_detail::nextBit(row, words_, 0); targetAtom >= 0 && level < RADII;
                 targetAtom = common_bond_detail::nextBit(row, words_, targetAtom + 1)) {
                clock_.tick();
                level = std::max(level, likeness(queryAtom, targetAtom));
            }
            likestLevel_[queryAtom] = level;
        }
        return likestLevel_[queryAtom];
    }

    /**
     * The number of taken targets that lead the list of neighbours of `targetAtom`, for the seed
     * of `round`. A seed only takes targets while it grows, so the count only rises, and it is
     * counted afresh for each seed.
     */
    int leadingTaken(int targetAtom, int round) {
        int& taken = leadingTaken_[targetAtom];
        if (leadingRound_[targetAtom] != round) {
            leadingRound_[targetAtom] = round;
            taken = 0;
        }
        while (taken < target_.degree(targetAtom) && targetToQuery_[target_.neighbour[target_.start[targetAtom] + taken]] >= 0) {
            clock_.tick();
            ++taken;
        }
        return taken;
    }

    /** The number of radii, from 0, at which two atoms look alike. */
    int likeness(int queryAtom, int targetAtom) const {
        int level = 0;
        while (level < RADII
               && querySignature_[size_t(level) * queryCount_ + queryAtom] == targetSignature_[size_t(level) * targetCount_ + targetAtom])
            ++level;
        return level;
    }

    /**
     * The mappings with the best score, at most `maxMappings`. Each is found once, as a mapping
     * is reached only from its highest ranked query atom and under the first sibling pair it
     * contains; scouting may meet a mapping before its turn, so the mappings held are kept
     * distinct. The search starts from the score `start`, which `seed` holds unless it is null.
     * It ends early once it holds enough mappings that reach `upperBound`. It needs rank()
     * first.
     */
    std::vector<std::vector<int>> collect(const std::vector<int>* seed, int64_t start, int64_t upperBound, size_t maxMappings) {
        bestScore_ = start;
        upperBound_ = upperBound;
        maxMappings_ = maxMappings;
        holdsSeed_ = seed != nullptr;
        mappings_.clear();
        held_.clear();
        if (seed) {
            heldSeed_ = *seed;
            mappings_.push_back(*seed);
            held_.insert(*seed);
        }
        findTwins();
        roots(queryCount_ - 1, targetCount_ - 1, INT64_MAX, INT64_MAX);
        return mappings_;
    }

    /**
     * Searches below each root pair in turn, from the target at place `firstPlace` of the query
     * atom of rank `firstRank`, until the pass can end or `nodeEnd` nodes are reached, and each
     * root pair for at most `rootNodes` nodes. Below a root only lower ranked query atoms are
     * mapped, so the roots of low rank have less to reach.
     */
    void roots(int firstRank, int firstPlace, int64_t nodeEnd, int64_t rootNodes) {
        const int64_t limit = nodeLimit_;
        for (int rankIndex = firstRank; rankIndex >= lowestRoot_ && !done() && nodes_ < nodeEnd && canImprove(rootBound(rankIndex));
             --rankIndex) {
            const int queryAtom = queryOrder_[rankIndex];
            for (int targetPosition = rankIndex == firstRank ? firstPlace : targetCount_ - 1;
                 targetPosition >= 0 && !done() && nodes_ < nodeEnd; --targetPosition) {
                const int targetAtom = targetOrder_[targetPosition];
                clock_.tick();
                if (!has(allowed_, queryAtom, targetAtom)) continue;
                if (rootRefuted(queryAtom, targetAtom)) {
                    ++refutedRoots_;
                    continue;
                }
                if (rootNodes == INT64_MAX) {
                    rootRank_ = rankIndex;
                    rootPlace_ = targetPosition;
                } else {
                    nodeLimit_ = std::min(limit, std::min(nodeEnd, nodes_ + rootNodes));
                }
                if (push(queryAtom, targetAtom, 0, upperBound_)) search();
                else if (!rootBonds_.empty()) rootBonds_[size_t(pair(queryAtom, targetAtom))] = bondsOf(lastBound_);
                nodeLimit_ = limit;
            }
        }
    }

    /**
     * Scouts once the search under way has mapped many pairs, to find a better mapping soon
     * when the first root pairs hold none: the state of the search is put aside, the root
     * pairs after the one it is below are searched for at most ScoutLimits::root nodes each
     * and ScoutLimits::budget in all, and the state is put back. A better mapping found raises
     * the best score, which only prunes more of what the search has left. Each gap to the next
     * scouting is ScoutLimits::growth times the one before, so scouting takes a small share of
     * a long search.
     */
    void scout() {
        scoutAt_ = INT64_MAX;
        ++scoutings_;
        // The levels of the search are put aside, and scouting gets levels of its own.
        std::vector<Level> levels(levels_.size());
        levels.swap(levels_);
        const std::vector<int> queryToTarget = queryToTarget_, targetToQuery = targetToQuery_, excluded = excluded_;
        const std::vector<uint64_t> allowed = allowed_;
        const int depth = depth_, commonBonds = commonBonds_, fragments = fragments_;
        depth_ = commonBonds_ = fragments_ = 0;
        std::fill(queryToTarget_.begin(), queryToTarget_.end(), -1);
        std::fill(targetToQuery_.begin(), targetToQuery_.end(), -1);
        excluded_.clear();
        allowed_ = compatible_;
        const int rank = rootRank_, place = rootPlace_;
        roots(place > 0 ? rank : rank - 1, place > 0 ? place - 1 : targetCount_ - 1, nodes_ + scoutLimits_.budget,
              scoutLimits_.root);
        rootRank_ = rank;
        rootPlace_ = place;
        levels_.swap(levels);
        queryToTarget_ = queryToTarget;
        targetToQuery_ = targetToQuery;
        excluded_ = excluded;
        allowed_ = allowed;
        depth_ = depth;
        commonBonds_ = commonBonds;
        fragments_ = fragments;
        scoutGap_ *= scoutLimits_.growth;
        scoutAt_ = nodes_ + scoutGap_;
    }

    /**
     * Whether a root pair can be skipped, as the bound of its state, when it was last pruned,
     * shows that nothing below it beats the best score. That bound holds for every later pass:
     * the pairs are never excluded at the top, and the upper bound only falls.
     */
    bool rootRefuted(int queryAtom, int targetAtom) const {
        if (rootBonds_.empty()) return false;
        const int bonds = rootBonds_[size_t(pair(queryAtom, targetAtom))];
        return bonds >= 0 && !canImprove(capped(std::min(queryCount_, targetCount_), bonds));
    }

    /**
     * A maximum mapping below `bound`. Unless the search is connected, the best score is sought
     * first level by level of common bonds, down from the bound to the score of the best
     * mapping known: each level asks for one mapping with at least that many bonds, and without
     * one the bound comes down. Twins can then be skipped all along, see twinExcluded(). A
     * level far above the best score prunes much and costs little, so this soon finds an
     * optimum close to the bound. Each level down repeats the work of the one above it, though,
     * so the levels together map at most LEVEL_BUDGET pairs. A level cut short refutes nothing.
     * A connected search has no levels: its bound, see groupedBound(), prunes as much near the
     * best mapping known, so the levels would mostly repeat work. One search from the best
     * mapping known then finishes the work, its threshold rising as better mappings are found.
     * When the levels run out of budget, that search starts from the seed grown the other way
     * round if that is better, see seedAlsoFrom(): the search is long, and the higher it starts
     * the more it prunes.
     */
    std::vector<std::vector<int>> collectByBonds(int64_t bound) {
        const int64_t budget = nodes_ + LEVEL_BUDGET;
        for (int triedBonds = bondsOf(bound); !connected_ && triedBonds > bondsOf(fallbackScore_); --triedBonds) {
            nodeLimit_ = budget;
            std::vector<std::vector<int>> mapping = collect(nullptr, score(0, triedBonds) - 1, bound, 1);
            nodeLimit_ = INT64_MAX;
            if (nodes_ >= budget) {
                // A mapping found in a level cut short has the bonds of the level, but perhaps not the most atoms.
                if (!mapping.empty() && scoreOf(mapping.front()) > fallbackScore_) {
                    fallback_ = mapping.front();
                    hasFallback_ = true;
                    fallbackScore_ = scoreOf(fallback_);
                }
                if (otherSeed_) {
                    const std::vector<int> other = otherSeed_();
                    if (!other.empty()) offerFallback(other, scoreOf(other));
                }
                break;
            }
            if (!mapping.empty()) {
                if (wanted_ == 1) return mapping;
                // The level holds a best mapping; one search from it finds the others of its score.
                fallback_ = mapping.front();
                hasFallback_ = true;
                fallbackScore_ = bound = scoreOf(fallback_);
                break;
            }
            // No mapping has that many bonds. With one bond fewer there are at most as many
            // atoms, or one more for each further fragment.
            bound = score(std::min(atomsOf(bound), mostAtoms(triedBonds - 1)), triedBonds - 1);
        }
        scoutAt_ = nodes_ + scoutLimits_.after;
        return collect(hasFallback_ ? &fallback_ : nullptr, fallbackScore_, bound, wanted_);
    }

    /**
     * Searches depth first from the root pair until it is taken off again. The top level tries
     * the targets of each of its candidates in turn, and then, unless the search is connected,
     * the roots of new fragments. It is taken off once nothing is left or its bound is too low.
     * The levels are an explicit stack, so large molecules do not run out of stack.
     */
    void search() {
        while (depth_ > 0) {
            clock_.tick();
            if (nodes_ >= scoutAt_) scout();
            Level& level = levels_[depth_ - 1];
            if (done() || !canImprove(level.bound)) {
                pop();
            } else if (level.nextTarget < static_cast<int>(level.targets.size())) {
                const int queryAtom = level.triedAtom, targetAtom = level.targets[level.nextTarget++];
                if (twinExcluded(queryAtom, targetAtom) || !push(queryAtom, targetAtom, level.fragmentStart, level.bound))
                    exclude(queryAtom, targetAtom);
            } else if (level.nextCandidate < static_cast<int>(level.candidates.size())) {
                // The targets tried so far are now excluded, which can lower the bound.
                if (level.nextCandidate > 0 && (level.bound = bound(level.fragmentStart, level.bound)) < 0) {
                    pop();
                } else {
                    level.triedAtom = level.candidates[level.nextCandidate++];
                    targets(level.triedAtom, level.targets);
                    level.nextTarget = 0;
                }
            } else if (!nextRoot(level)) {
                pop();
            }
        }
    }

    /**
     * Maps a pair as a new level of the fragment that starts at level `fragmentStart`, with the
     * bound of the state it extends until it has its own. False, with the pair left unmapped,
     * if the bound prunes it.
     */
    bool push(int queryAtom, int targetAtom, int fragmentStart, int64_t bound) {
        Level& level = map(queryAtom, targetAtom, fragmentStart, bound);
        record(fragmentStart);
        return isLeaf() || branch(level, fragmentStart, bound);
    }

    /**
     * Maps a pair as the level at depth_. The level records what taking the pair off restores,
     * the common bonds and the excluded pairs before it, and for the bound the query and target
     * bonds with no end mapped once it is. It has nothing to try until branch() gives it
     * candidates.
     */
    Level& map(int queryAtom, int targetAtom, int fragmentStart, int64_t bound) {
        countNode();
        Level& level = levels_[depth_];
        level.queryAtom = queryAtom;
        level.targetAtom = targetAtom;
        level.startsFragment = fragmentStart == depth_;
        if (level.startsFragment) ++fragments_;
        level.commonBondsBefore = commonBonds_;
        level.excludedBefore = static_cast<int>(excluded_.size());
        level.freeQueryBonds = depth_ > 0 ? levels_[depth_ - 1].freeQueryBonds : query_.bondCount;
        level.freeTargetBonds = depth_ > 0 ? levels_[depth_ - 1].freeTargetBonds : target_.bondCount;
        for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place)
            if (queryToTarget_[query_.neighbour[place]] < 0) --level.freeQueryBonds;
        for (int place = target_.start[targetAtom]; place < target_.start[targetAtom + 1]; ++place)
            if (targetToQuery_[target_.neighbour[place]] < 0) --level.freeTargetBonds;
        level.fragmentStart = fragmentStart;
        level.bound = bound;
        level.candidates.clear();
        level.targets.clear();
        level.nextCandidate = level.nextTarget = level.rootRank = 0;
        commonBonds_ += common(queryAtom, targetAtom);
        queryToTarget_[queryAtom] = targetAtom;
        targetToQuery_[targetAtom] = queryAtom;
        ++depth_;
        return level;
    }

    /**
     * Whether the mapping is a leaf, with nothing left to try: the search holds what it wants,
     * every query atom up to the rank of the root is mapped, or every target atom is, or the
     * mapping has the common bonds of the upper bound.
     */
    bool isLeaf() const {
        // A mapping that drops bonds for its conditions may select fewer than its common bonds,
        // and a larger mapping may then select more, up to the upper bound.
        return done() || depth_ == queryRank_[levels_[0].queryAtom] + 1 || depth_ == targetCount_
            || (commonBonds_ == bondsOf(upperBound_) && (!conditions_ || !conditions_->dropsBonds()));
    }

    /**
     * Gives the top level its own bound and the candidates whose targets it tries. Unless the
     * search is connected, a fragment with a bond also offers the roots of new fragments once
     * its candidates are done; these rank below its own root. False, with the pair taken off,
     * if the bound prunes the level.
     */
    bool branch(Level& level, int fragmentStart, int64_t bound) {
        level.bound = this->bound(fragmentStart, bound);
        if (level.bound < 0) {
            takeOff();
            return false;
        }
        candidates(level.candidates);
        clock_.tick(static_cast<int64_t>(level.candidates.size()));
        if (!connected_ && fragmentBonds(fragmentStart) >= MIN_BONDS)
            level.rootRank = queryRank_[levels_[fragmentStart].queryAtom];
        return true;
    }

    /** The common bonds of the open fragment, which starts at level `fragmentStart`. */
    int fragmentBonds(int fragmentStart) const { return commonBonds_ - levels_[fragmentStart].commonBondsBefore; }

    /** Takes the top level off, and excludes its pair for its later siblings and everything below them. */
    void pop() {
        const Level& level = levels_[depth_ - 1];
        const int queryAtom = level.queryAtom, targetAtom = level.targetAtom;
        takeOff();
        if (depth_ > 0) exclude(queryAtom, targetAtom);
    }

    /**
     * Counts a pair mapped or an embedding atom placed, and stops the search with Stop instead
     * once nodeCap_ are counted.
     */
    void countNode() {
        if (nodes_ >= nodeCap_) {
            capReached_ = true;
            throw Stop{};
        }
        ++nodes_;
    }

    /** Unmaps the pair of the top level, and allows again the pairs excluded since it was mapped. */
    void takeOff() {
        const Level& level = levels_[--depth_];
        if (level.startsFragment) --fragments_;
        while (static_cast<int>(excluded_.size()) > level.excludedBefore) {
            const int packed = excluded_.back();
            excluded_.pop_back();
            set(allowed_, pairQuery(packed), pairTarget(packed));
        }
        queryToTarget_[level.queryAtom] = -1;
        targetToQuery_[level.targetAtom] = -1;
        commonBonds_ = level.commonBondsBefore;
    }

    /** Excludes a pair, which must be allowed. */
    void exclude(int queryAtom, int targetAtom) {
        clear(allowed_, queryAtom, targetAtom);
        excluded_.push_back(pair(queryAtom, targetAtom));
    }

    /**
     * Offers the next root of a new fragment once the candidates of the level are done. This
     * closes the fragment, as every pair that extends it is excluded. A root is the next free
     * query atom ranked below the root of the closed fragment, so a mapping is reached only
     * from the highest ranked atom of each fragment. False when no root is left, when the bound
     * of those left is too low, or when the search is connected. That bound falls with their
     * rank.
     */
    bool nextRoot(Level& level) {
        if (level.rootRank == 0) return false;
        // The fragment is closed once it is set to start at the next level, where the roots
        // offered start fragments of their own. Its atoms then lie in no larger fragment, and the
        // fragments only add up, so a fragment too small or one too many leaves nothing to find below.
        if (level.fragmentStart < depth_) {
            if (hasPendantPair(level.fragmentStart)
                || (conditions_ && (depth_ - level.fragmentStart < conditions_->minFragmentAtoms()
                                    || fragments_ >= conditions_->maxFragments())))
                return false;
            level.fragmentStart = depth_;
        }
        while (--level.rootRank >= 0) {
            const int queryAtom = queryOrder_[level.rootRank];
            if (queryToTarget_[queryAtom] >= 0) continue;
            level.bound = rootBound(level.rootRank);
            if (!canImprove(level.bound)) return false;
            level.triedAtom = queryAtom;
            freeTargets(queryAtom, level.targets);
            level.nextTarget = 0;
            if (!level.targets.empty()) return true;
        }
        return false;
    }

    /**
     * Keeps the mapping if it beats the best score, or ties with it while fewer than wanted are
     * held, and otherwise keeps it as the fallback if it beats that. The open fragment, from
     * level `fragmentStart`, must have a common bond.
     */
    void record(int fragmentStart) {
        const int64_t current = score(depth_, commonBonds_);
        const bool full = mappings_.size() == maxMappings_;
        if (current < bestScore_ || (current == bestScore_ && full)) {
            // Below the threshold of a level, a mapping better than the fallback becomes the fallback.
            if (fragmentBonds(fragmentStart) >= MIN_BONDS) offerFallback(queryToTarget_, current);
            return;
        }
        // A fragment needs a common bond. The open one may still get it further down.
        if (fragmentBonds(fragmentStart) < MIN_BONDS) return;
        if (!conditions_) {
            keep(current, queryToTarget_);
            return;
        }
        // The best selections that meet the conditions score at most as much. One better than the
        // fallback serves as that while the best score is out of reach.
        selections_.clear();
        const int64_t least = std::min(fallbackScore_ + 1, full ? bestScore_ + 1 : bestScore_);
        const int64_t valid = conditions_->select(queryToTarget_, least, maxMappings_, selections_);
        if (valid < 0) return;
        if (valid < bestScore_ || (valid == bestScore_ && full)) {
            if (valid > fallbackScore_) {
                fallback_ = selections_.front();
                hasFallback_ = true;
                fallbackScore_ = valid;
            }
            return;
        }
        for (const auto& selection : selections_)
            if (!keep(valid, selection)) break;
    }

    /**
     * Holds `record`, of score `value`, which is no lower than the best: a better one replaces
     * those held. False, holding nothing, once as many as wanted are held.
     */
    bool keep(int64_t value, const std::vector<int>& record) {
        if (value > bestScore_) {
            bestScore_ = value;
            mappings_.clear();
            held_.clear();
            holdsSeed_ = false;
        } else if (holdsSeed_ && record == heldSeed_) {
            // The seed is held from the start, so it is not added again when found.
            holdsSeed_ = false;
            return true;
        }
        if (mappings_.size() == maxMappings_) return false;
        // A single mapping is found once; several may meet one met before while scouting.
        if (maxMappings_ > 1 && !held_.insert(record).second) return true;
        mappings_.push_back(record);
        return true;
    }

    /**
     * Whether the pass can end, as it holds as many mappings as wanted with a score that no
     * mapping beats, or as it has mapped as many pairs as it may.
     */
    bool done() const {
        return (mappings_.size() == maxMappings_ && bestScore_ == upperBound_) || nodes_ >= nodeLimit_;
    }

    /**
     * The bound of the state whose open fragment starts at level `fragmentStart`, or -1 if the
     * state is pruned. A mapping below adds to the open fragment only open pairs that it reaches
     * through compatible bonds, so the bound counts these. The first layer of them gives the
     * candidates.
     */
    int64_t bound(int fragmentStart, int64_t current) {
        reachStart();
        int open = 0;
        for (int i = fragmentStart; i < depth_; ++i) {
            const int queryLinks = reach_.queryLinks, targetLinks = reach_.targetLinks;
            reachLink(levels_[i].queryAtom, levels_[i].targetAtom);
            open += std::min(reach_.queryLinks - queryLinks, reach_.targetLinks - targetLinks);
        }
        reach_.firstLayer = static_cast<int>(reach_.pairs.size());
        const Level& top = levels_[depth_ - 1];
        // The atoms of the upper bound cap only a mapping with as many bonds. A mapping with
        // fewer bonds may have more atoms.
        const int most = std::min(queryCount_, targetCount_);
        // Centres that already disagree on three ligands keep disagreeing below, so a mapping
        // there leaves bonds out for them, see CommonBondConditions::persistentLoss().
        const int loss = conditions_ ? conditions_->persistentLoss(queryToTarget_) : 0;
        const int mostBonds = commonBonds_ + open + std::min(top.freeQueryBonds, top.freeTargetBonds) - loss;
        int64_t result = std::min(current, capped(most, mostBonds));
        // A new atom lies within as many bonds of the mapped atoms as a better mapping adds. One
        // that drops bonds for its conditions may keep fewer of the common bonds so far, so it may
        // add up to all its bonds.
        const int newBonds = conditions_ && conditions_->dropsBonds() ? bondsOf(result) : bondsOf(result) - commonBonds_;
        if (canImprove(result))
            result = std::min(result, edgeBound(queryRank_[levels_[fragmentStart].queryAtom], open, newBonds));
        reachFinish();
        lastBound_ = result;
        return canImprove(result) ? result : -1;
    }

    /**
     * The bound of the roots offered up to `rankLimit`. The fragment is closed, so no bond links
     * it. Centres that already disagree on three ligands lower it, as in bound().
     */
    int64_t rootBound(int rankLimit) {
        reachStart();
        const int64_t bound = edgeBound(rankLimit, 0, 0);
        const int loss = depth_ > 0 && conditions_ ? conditions_->persistentLoss(queryToTarget_) : 0;
        return loss == 0 ? bound : capped(atomsOf(bound), bondsOf(bound) - loss);
    }

    /**
     * The score of a mapping, a record (see CommonBondConditions::select()): its atoms, and the
     * compatible bonds between them that the target has between their images, less those it
     * leaves out.
     */
    int64_t scoreOf(const std::vector<int>& mapping) const {
        int atomCount = 0, bondCount = -leftOut(mapping, queryCount_);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
            if (mapping[queryAtom] >= 0) {
                ++atomCount;
                bondCount += commonBondsAt(mapping, queryAtom);
            }
        return score(atomCount, bondCount);
    }

    /**
     * The compatible bonds from a query atom to the lower mapped atoms of a mapping that the
     * target has between their images, so that each common bond counts once.
     */
    int commonBondsAt(const std::vector<int>& mapping, int queryAtom) const {
        const int targetAtom = mapping[queryAtom];
        int count = 0;
        for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place) {
            const int other = query_.neighbour[place];
            if (other >= queryAtom || mapping[other] < 0) continue;
            const int targetBond = target_.bondBetween(targetAtom, mapping[other]);
            if (targetBond >= 0 && graph_.bondsMatch(query_.neighbourBond[place], targetBond)) ++count;
        }
        return count;
    }

    /** A bound, capped by the upper bound: a mapping with the common bonds of the upper bound has no more atoms than it. */
    int64_t capped(int atomCount, int bondCount) const {
        const int mostAtomsBound = atomsOf(upperBound_), mostBonds = bondsOf(upperBound_);
        if (bondCount >= mostBonds) return score(std::min(atomCount, mostAtomsBound), mostBonds);
        return score(atomCount, bondCount);
    }

    /** Whether a state with this bound can hold a better mapping, or a tie while fewer than wanted are held. */
    bool canImprove(int64_t bound) const {
        return bound > bestScore_ || (bound == bestScore_ && mappings_.size() < maxMappings_);
    }

    /** The query atoms of the first layer of reached pairs: not terminal first, then fewer pairs, then higher rank. */
    void candidates(std::vector<int>& sorted) {
        layer_.assign(reach_.pairs.begin(), reach_.pairs.begin() + reach_.firstLayer);
        std::sort(layer_.begin(), layer_.end());
        keys_.clear();
        for (size_t from = 0, to; from < layer_.size(); from = to) {
            const int queryAtom = pairQuery(layer_[from]);
            to = from + 1;
            while (to < layer_.size() && pairQuery(layer_[to]) == queryAtom) ++to;
            keys_.push_back(candidateKey(queryAtom, static_cast<int>(to - from)));
        }
        std::sort(keys_.begin(), keys_.end());
        sorted.clear();
        for (int64_t key : keys_) sorted.push_back(candidateOf(key));
    }

    /** The sort key of a candidate: whether it is terminal, its number of pairs, then its rank in reverse. */
    int64_t candidateKey(int queryAtom, int pairs) const {
        return (query_.degree(queryAtom) == 1 ? int64_t(1) << 62 : 0) | (int64_t(pairs) << 31)
            | int64_t(queryCount_ - 1 - queryRank_[queryAtom]);
    }

    int candidateOf(int64_t key) const { return queryOrder_[queryCount_ - 1 - static_cast<int>(key & INT_MAX)]; }

    /**
     * The targets of a query atom: the free neighbours of the images of its mapped neighbours,
     * joined by compatible bonds, in pairs not excluded. The most alike come first, then by
     * descending rank. The order decides only which mappings are found first.
     */
    void targets(int queryAtom, std::vector<int>& result) {
        keys_.clear();
        for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
            const int neighbourImage = queryToTarget_[query_.neighbour[queryPlace]];
            if (neighbourImage < 0) continue;
            clock_.tick(1 + target_.degree(neighbourImage));
            for (int targetPlace = target_.start[neighbourImage]; targetPlace < target_.start[neighbourImage + 1]; ++targetPlace) {
                const int targetAtom = target_.neighbour[targetPlace];
                if (targetToQuery_[targetAtom] < 0 && has(allowed_, queryAtom, targetAtom)
                    && graph_.bondsMatch(query_.neighbourBond[queryPlace], target_.neighbourBond[targetPlace]))
                    keys_.push_back(targetKey(queryAtom, targetAtom));
            }
        }
        // A target next to two mapped neighbours is listed twice, with the same key.
        std::sort(keys_.begin(), keys_.end());
        result.clear();
        for (size_t i = 0; i < keys_.size(); ++i)
            if (i == 0 || keys_[i] != keys_[i - 1]) result.push_back(targetOf(keys_[i]));
    }

    /** The sort key of a target of a query atom: how unlike the two are, then the rank in reverse. */
    int64_t targetKey(int queryAtom, int targetAtom) const {
        return (int64_t(RADII - likeness(queryAtom, targetAtom)) << 32) | int64_t(targetCount_ - 1 - targetRank_[targetAtom]);
    }

    int targetOf(int64_t key) const { return targetOrder_[targetCount_ - 1 - static_cast<int>(key & 0xffffffff)]; }

    /** The free targets of a root, by descending rank. */
    void freeTargets(int queryAtom, std::vector<int>& result) {
        result.clear();
        clock_.tick(targetCount_);
        for (int targetPosition = targetCount_ - 1; targetPosition >= 0; --targetPosition)
            if (targetToQuery_[targetOrder_[targetPosition]] < 0 && has(allowed_, queryAtom, targetOrder_[targetPosition]))
                result.push_back(targetOrder_[targetPosition]);
    }

    /** The compatible bonds from a query atom to mapped atoms that the target has from a target atom to their images. */
    int common(int queryAtom, int targetAtom) {
        int count = 0;
        for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place) {
            const int neighbourImage = queryToTarget_[query_.neighbour[place]];
            if (neighbourImage < 0) continue;
            clock_.tick(1 + target_.degree(targetAtom));
            const int targetBond = target_.bondBetween(targetAtom, neighbourImage);
            if (targetBond >= 0 && graph_.bondsMatch(query_.neighbourBond[place], targetBond)) ++count;
        }
        return count;
    }

    /**
     * Finds the twins once, before any pair is excluded. A query that fits whole needs them only
     * while its embeddings are tried under conditions.
     */
    void findTwins() {
        if (twinsKnown_) return;
        queryTwin_ = twins(query_, true);
        targetTwin_ = twins(target_, false);
        twinsKnown_ = true;
    }

    /**
     * The twins of one molecule, each linked to the next in a ring, -1 for none. Twins are
     * terminal atoms with the same neighbour that match the same atoms by the same bonds and,
     * under conditions, carry with that neighbour no stereo tag that the conditions read, see
     * CommonBondConditions::untagged(). Swapping two takes every mapping onto one of the same
     * score.
     */
    std::vector<int> twins(const Side& side, bool onQuery) {
        std::vector<int> next(side.atomCount, -1);
        for (int parent = 0; parent < side.atomCount; ++parent) {
            if (!untagged(onQuery, parent)) continue;
            for (int place = side.start[parent]; place < side.start[parent + 1]; ++place) {
                const int first = side.neighbour[place];
                int last = first;
                if (next[first] >= 0 || side.degree(first) != 1 || !untagged(onQuery, first)) continue;
                for (int laterPlace = place + 1; laterPlace < side.start[parent + 1]; ++laterPlace) {
                    const int other = side.neighbour[laterPlace];
                    if (next[other] < 0 && side.degree(other) == 1 && untagged(onQuery, other)
                        && alike(onQuery, first, side.neighbourBond[place], other, side.neighbourBond[laterPlace])) {
                        next[last] = other;
                        last = other;
                    }
                }
                if (last != first) next[last] = first;
            }
        }
        return next;
    }

    /** Whether `atom` of one molecule carries no stereo tag that the conditions read, if there are any. */
    bool untagged(bool onQuery, int atom) const { return !conditions_ || conditions_->untagged(onQuery, atom); }

    /** Whether two atoms of one molecule, with their bonds, match alike. */
    bool alike(bool onQuery, int atom, int bond, int other, int otherBond) {
        clock_.tick(int64_t(words_) + queryCount_ + query_.bondCount + target_.bondCount);
        if (onQuery) {
            for (int word = 0; word < words_; ++word)
                if (allowed_[size_t(atom) * words_ + word] != allowed_[size_t(other) * words_ + word]) return false;
            for (int targetBond = 0; targetBond < target_.bondCount; ++targetBond)
                if (graph_.bondsMatch(bond, targetBond) != graph_.bondsMatch(otherBond, targetBond)) return false;
        } else {
            for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
                if (has(allowed_, queryAtom, atom) != has(allowed_, queryAtom, other)) return false;
            for (int queryBond = 0; queryBond < query_.bondCount; ++queryBond)
                if (graph_.bondsMatch(queryBond, bond) != graph_.bondsMatch(queryBond, otherBond)) return false;
        }
        return true;
    }

    /**
     * Whether a pair can be skipped, as a free twin of one has its pair with the other
     * excluded. Swapping the twins turns every mapping below this pair into one of the same
     * score below the excluded pair, which was tried before. So this holds only while no
     * mapping of that score could be recorded: when as many are held as wanted, and also when
     * nothing but the seed is held.
     */
    bool twinExcluded(int queryAtom, int targetAtom) const {
        if (conditions_ && !conditions_->keepsTwins()) return false;
        const bool full = mappings_.size() == maxMappings_;
        const bool onlySeed = mappings_.size() == (holdsSeed_ ? 1U : 0U);
        if (!full && !onlySeed) return false;
        for (int twin = queryTwin_[queryAtom]; twin >= 0 && twin != queryAtom; twin = queryTwin_[twin])
            if (queryToTarget_[twin] < 0 && !has(allowed_, twin, targetAtom)) return true;
        for (int twin = targetTwin_[targetAtom]; twin >= 0 && twin != targetAtom; twin = targetTwin_[twin])
            if (targetToQuery_[twin] < 0 && !has(allowed_, queryAtom, twin)) return true;
        return false;
    }

    /**
     * Whether the open fragment leaves out a pendant pair: a free terminal query atom next to a
     * mapped atom whose image has a free terminal target atom it may map onto by a compatible
     * bond. No mapping below the roots the level offers is then a maximum one, as adding that
     * pair beats it.
     */
    bool hasPendantPair(int fragmentStart) {
        if (conditions_ && !conditions_->keepsPendantPairs()) return false;
        for (int index = fragmentStart; index < depth_; ++index) {
            const int queryAtom = levels_[index].queryAtom, targetAtom = levels_[index].targetAtom;
            clock_.tick(int64_t(query_.degree(queryAtom)) * target_.degree(targetAtom));
            for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
                const int queryNeighbour = query_.neighbour[queryPlace];
                if (queryToTarget_[queryNeighbour] >= 0 || query_.degree(queryNeighbour) != 1) continue;
                for (int targetPlace = target_.start[targetAtom]; targetPlace < target_.start[targetAtom + 1]; ++targetPlace) {
                    const int targetNeighbour = target_.neighbour[targetPlace];
                    if (targetToQuery_[targetNeighbour] < 0 && target_.degree(targetNeighbour) == 1
                        && has(compatible_, queryNeighbour, targetNeighbour)
                        && graph_.bondsMatch(query_.neighbourBond[queryPlace], target_.neighbourBond[targetPlace]))
                        return true;
                }
            }
        }
        return false;
    }

    /**
     * The search for the whole query, if it may fit: the query is connected and no larger than
     * the target, and each query atom can have a target atom of its own with at least as many
     * bonds. Null if not. With `prunesTwins` it tries only embeddings that place the twins of
     * each molecule in order, see Embedding::twinsOutOfOrder().
     */
    std::unique_ptr<Embedding> embedding(bool prunesTwins);

    // ---- The reached pairs ----

    /** Starts a use with nothing reached. Every mark is cleared when the stamp would overflow. */
    void reachStart() {
        if (++reach_.stamp == INT_MAX) {
            std::fill(reach_.queryBondStamp.begin(), reach_.queryBondStamp.end(), 0);
            std::fill(reach_.targetBondStamp.begin(), reach_.targetBondStamp.end(), 0);
            reach_.stamp = 1;
        }
        reach_.pairs.clear();
        reach_.queryLinks = reach_.targetLinks = 0;
    }

    /** Reaches the open pairs next to a mapped pair. */
    void reachLink(int queryAtom, int targetAtom) {
        const int from = target_.start[targetAtom], to = target_.start[targetAtom + 1];
        const int rootRank = queryRank_[levels_[0].queryAtom];
        for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
            const int queryNeighbour = query_.neighbour[queryPlace], queryBond = query_.neighbourBond[queryPlace];
            if (queryToTarget_[queryNeighbour] >= 0 || queryRank_[queryNeighbour] > rootRank) continue;
            clock_.tick(1 + to - from);
            for (int targetPlace = from; targetPlace < to; ++targetPlace) {
                const int targetNeighbour = target_.neighbour[targetPlace], targetBond = target_.neighbourBond[targetPlace];
                if (targetToQuery_[targetNeighbour] >= 0 || !has(allowed_, queryNeighbour, targetNeighbour)
                    || !graph_.bondsMatch(queryBond, targetBond)) continue;
                if (reach_.queryBondStamp[queryBond] != reach_.stamp) {
                    reach_.queryBondStamp[queryBond] = reach_.stamp;
                    ++reach_.queryLinks;
                }
                if (reach_.targetBondStamp[targetBond] != reach_.stamp) {
                    reach_.targetBondStamp[targetBond] = reach_.stamp;
                    ++reach_.targetLinks;
                }
                if (!has(reach_.seen, queryNeighbour, targetNeighbour)) {
                    set(reach_.seen, queryNeighbour, targetNeighbour);
                    reach_.pairs.push_back(pair(queryNeighbour, targetNeighbour));
                }
            }
        }
    }


    /** Ends a use by clearing the bits of the reached pairs. */
    void reachFinish() {
        for (int packed : reach_.pairs) clear(reach_.seen, pairQuery(packed), pairTarget(packed));
    }

    // ---- The bound from the bond labels ----

    /**
     * The bound of the state, where the open atoms are those ranked up to `rankLimit` that have
     * a pair not excluded. A common bond and its image share a label. An open atom and its
     * image get no more new bonds of a label than either has open, so pairing the open atoms of
     * each label bounds its new bonds. Where a label holds fine labels that do not all match,
     * its new bonds are also at most the most bonds of matching fine labels that the open atoms
     * can pair off, and the new atoms of an atom class at most the most open atoms of matching
     * kinds. A connected extension of a mapping, of at most `newBonds` new bonds, is bounded
     * by groups of the open atoms near it instead, see groupedBound(). The `open` new bonds at
     * the mapped atoms are already counted. A fragment of b new bonds has at most b + 1 atoms,
     * and a connected extension of a mapping at most b new atoms.
     */
    int64_t edgeBound(int rankLimit, int open, int newBonds) {
        std::fill(unusedTargets_.begin(), unusedTargets_.end(), ~uint64_t(0));
        for (int i = 0; i < depth_; ++i) clear(unusedTargets_, 0, levels_[i].targetAtom);
        std::fill(reachableTargets_.begin(), reachableTargets_.end(), 0);
        if (connected_ && depth_ > 0) return groupedBound(rankLimit, open, newBonds);
        for (int rankIndex = 0; rankIndex <= rankLimit; ++rankIndex) {
            const int queryAtom = queryOrder_[rankIndex];
            if (queryToTarget_[queryAtom] < 0 && markTargets(queryAtom)) queryLabels_.markOpen(queryAtom);
        }
        for (int targetAtom = nextBit(reachableTargets_.data(), words_, 0); targetAtom >= 0;
             targetAtom = nextBit(reachableTargets_.data(), words_, targetAtom + 1))
            targetLabels_.markOpen(targetAtom);
        countOpen(queryLabels_, query_, true);
        countOpen(targetLabels_, target_, false);
        int added = open, newAtoms = 0;
        pairAll(added, newAtoms);
        queryLabels_.clearOpen();
        targetLabels_.clearOpen();
        const int mostNew = !connected_ ? added + added / MIN_BONDS : mostAtoms(added);
        return capped(depth_ + std::min(newAtoms, mostNew), commonBonds_ + added);
    }

    /** Adds the unused targets of the pairs of a free query atom to the reachable targets, and tells whether it has any. */
    bool markTargets(int queryAtom) {
        uint64_t any = 0;
        clock_.tick(words_);
        for (int word = 0; word < words_; ++word) {
            const uint64_t reachable = allowed_[size_t(queryAtom) * words_ + word] & unusedTargets_[word];
            reachableTargets_[word] |= reachable;
            any |= reachable;
        }
        return any != 0;
    }

    /**
     * The bound of a connected extension of the mapping by at most `newBonds` new bonds. A new
     * atom is joined to a mapped atom by a path of new atoms and common bonds. The first pair
     * of the path is one of the pairs reached from the mapped pairs, see reachLink(), and the
     * path goes on through open atoms and bonds that have a label, no more than `newBonds`
     * bonds in all. The query atoms so reached are open if they are ranked up to `rankLimit`
     * and have a pair not excluded, and the target atoms so reached if they are the targets of
     * such pairs of the open query atoms. The path lies in one part of the open query atoms,
     * linked by bonds that have a label, and its image in one part of the open target atoms,
     * and its first pair joins the two. So the parts fall into groups joined by the pairs
     * reached, and each group pairs off its own bonds and atoms; a group of one molecule only
     * gets nothing. Pairing off within each group pairs off no more than over all of them.
     */
    int64_t groupedBound(int rankLimit, int open, int newBonds) {
        const int queryOpen = spread(queryLabels_, query_, true, rankLimit, newBonds);
        spread(targetLabels_, target_, false, rankLimit, newBonds);
        int64_t work = static_cast<int64_t>(reach_.pairs.size());
        if (newBonds > 0)
            for (int packed : reach_.pairs) detail::unionFindJoin(partParent_, pairQuery(packed), queryCount_ + pairTarget(packed));
        // The open atoms of each group by a counting sort on the root of the group, the query atoms first.
        int groups = 0;
        const int openCount = queryOpen + targetLabels_.openCount;
        work += openCount;
        clock_.tick(work);
        for (int index = 0; index < openCount; ++index) {
            const int root = detail::unionFindRoot(partParent_, openItem(index, queryOpen));
            if (groupMark_[root] != markStamp_) {
                groupMark_[root] = markStamp_;
                groupOf_[root] = groups;
                groupStart_[++groups] = 0;
                groupSides_[groups - 1] = 0;
            }
            ++groupStart_[groupOf_[root] + 1];
            groupSides_[groupOf_[root]] |= index < queryOpen ? 1 : 2;
        }
        groupStart_[0] = 0;
        for (int group = 0; group < groups; ++group) groupStart_[group + 1] += groupStart_[group];
        for (int index = 0; index < openCount; ++index) {
            const int item = openItem(index, queryOpen);
            groupItems_[groupStart_[groupOf_[detail::unionFindRoot(partParent_, item)]]++] = item;
        }
        int added = open, newAtoms = 0;
        for (int group = 0, from = 0; group < groups; ++group) {
            const int to = groupStart_[group];
            if (groupSides_[group] == 3) {
                // The query atoms of the group come first.
                int split = from;
                while (groupItems_[split] < queryCount_) ++split;
                countAtoms<true>(queryLabels_, query_, true, groupItems_.data() + from, split - from, 0);
                countAtoms<true>(targetLabels_, target_, false, groupItems_.data() + split, to - split, queryCount_);
                pairCounted(added, newAtoms);
            }
            from = to;
        }
        nextMark();
        queryLabels_.clearOpen();
        targetLabels_.clearOpen();
        return capped(depth_ + std::min(newAtoms, added), commonBonds_ + added);
    }

    /** The union-find item of the open atom at `index`, the query atoms first and then the target atoms after queryCount_. */
    int openItem(int index, int queryOpen) const {
        return index < queryOpen ? queryLabels_.openAtoms[index] : queryCount_ + targetLabels_.openAtoms[index - queryOpen];
    }

    /** Starts new marks for the groups, clearing every mark when the stamp would overflow. */
    void nextMark() {
        if (++markStamp_ == INT_MAX) {
            std::fill(groupMark_.begin(), groupMark_.end(), 0);
            markStamp_ = 1;
        }
    }

    /**
     * Opens the atoms of one molecule in layers: first those of the pairs reached from the
     * mapped pairs, and then their neighbours through open atoms and bonds that have a label,
     * up to `newBonds` layers. Atoms linked by such bonds are joined in the union-find of the
     * groups. Each open query atom adds the targets of its pairs to the reachable targets,
     * which are the target atoms that may open. Gives the number of atoms opened.
     */
    int spread(Labels& labels, const Side& side, bool onQuery, int rankLimit, int newBonds) {
        const uint32_t stamp = labels.nextStamp();
        const int offset = onQuery ? 0 : queryCount_;
        int64_t work = 0;
        for (int i = 0; i < depth_; ++i) labels.reached[onQuery ? levels_[i].queryAtom : levels_[i].targetAtom] = stamp;
        // The first layer holds the atoms of the pairs reached from the mapped pairs.
        for (size_t i = 0; i < reach_.pairs.size() && newBonds > 0; ++i, ++work)
            reach(labels, onQuery, onQuery ? pairQuery(reach_.pairs[i]) : pairTarget(reach_.pairs[i]), rankLimit, stamp, offset);
        for (int index = 0, layerEnd = labels.openCount, distance = 1; index < labels.openCount; ++index) {
            if (index == layerEnd) {
                ++distance;
                layerEnd = labels.openCount;
            }
            const int atom = labels.openAtoms[index];
            work += 1 + side.degree(atom);
            for (int place = side.start[atom]; place < side.start[atom + 1]; ++place) {
                const int other = side.neighbour[place];
                if (labels.slot[place] < 0) continue;
                // An atom first met from the last layer would lie a bond too far.
                if (labels.reached[other] != stamp && distance < newBonds)
                    reach(labels, onQuery, other, rankLimit, stamp, offset);
                if (labels.isOpen[other]) detail::unionFindJoin(partParent_, offset + atom, offset + other);
            }
        }
        clock_.tick(work);
        return labels.openCount;
    }

    /**
     * Marks an atom met in spread() as reached, and opens it if it may be mapped: a free query
     * atom ranked up to `rankLimit` with a pair not excluded, or a reachable target atom.
     */
    void reach(Labels& labels, bool onQuery, int atom, int rankLimit, uint32_t stamp, int offset) {
        if (labels.reached[atom] == stamp) return;
        labels.reached[atom] = stamp;
        const bool mayOpen = onQuery
            ? queryToTarget_[atom] < 0 && queryRank_[atom] <= rankLimit && markTargets(atom)
            : has(reachableTargets_, 0, atom);
        if (!mayOpen) return;
        labels.markOpen(atom);
        partParent_[offset + atom] = offset + atom;
    }

    /** Counts every open atom of one molecule, see countAtoms(). */
    void countOpen(Labels& labels, const Side& side, bool onQuery) {
        countAtoms<false>(labels, side, onQuery, labels.openAtoms.data(), labels.openCount, 0);
    }

    /**
     * Counts `count` open atoms of one molecule, each less `offset`, by slot and by their
     * number of open bonds in it. One with an open bond or a link is also counted by class and,
     * where a class holds kinds that do not all match, by kind. Each bond to a higher open atom
     * of a label whose fine labels do not all match is counted by fine label. With `listed`,
     * the labels and the classes counted are listed for pairCounted().
     */
    template <bool listed>
    void countAtoms(Labels& labels, const Side& side, bool onQuery, const int* atoms, int count, int offset) {
        const Groups& labelGroups = graph_.labelGroups;
        const Groups& atomGroups = graph_.atomGroups;
        const std::vector<int>& linkStamp = onQuery ? reach_.queryBondStamp : reach_.targetBondStamp;
        int64_t work = count;
        for (int index = 0; index < count; ++index) {
            const int atom = atoms[index] - offset;
            int slots = 0;
            bool link = false;
            work += side.degree(atom);
            for (int place = side.start[atom]; place < side.start[atom + 1]; ++place) {
                const int slot = labels.slot[place], other = side.neighbour[place];
                if (slot < 0 || !labels.isOpen[other]) {
                    link = link || linkStamp[side.neighbourBond[place]] == reach_.stamp;
                    continue;
                }
                if (slotCount_[slot]++ == 0) slotsAt_[slots++] = slot;
                if (other > atom && !labelGroups.complete[slot / 2]) labels.countFine(side.fineLabel[side.neighbourBond[place]]);
            }
            for (int slotIndex = 0; slotIndex < slots; ++slotIndex) {
                const int slot = slotsAt_[slotIndex];
                ++labels.histogram[size_t(slot) * histogramWidth_ + slotCount_[slot]];
                slotCount_[slot] = 0;
                if (listed && !labelListed_[slot / 2]) {
                    labelListed_[slot / 2] = 1;
                    countedLabels_.push_back(slot / 2);
                }
            }
            if (slots > 0 || link) {
                const int symbol = labels.symbol[atom];
                if (labels.bySymbol[symbol]++ == 0 && listed && !symbolListed_[symbol]) {
                    symbolListed_[symbol] = 1;
                    countedSymbols_.push_back(symbol);
                }
                if (!atomGroups.complete[symbol]) labels.countKind(side.atomKind[atom]);
            }
        }
        clock_.tick(work);
    }

    /**
     * Adds the new bonds and atoms that the atoms counted can pair off: the bonds of each label
     * and the atoms of each class. Clears the counts.
     */
    void pairAll(int& bonds, int& atoms) {
        Groups& labelGroups = graph_.labelGroups;
        Groups& atomGroups = graph_.atomGroups;
        clock_.tick(int64_t(slotSymbol_.size()) * histogramWidth_ + int64_t(queryLabels_.bySymbol.size()));
        for (int slot = 0; slot < static_cast<int>(slotSymbol_.size()); slot += 2) {
            const int label = slot / 2;
            int pairs = slotSymbol_[slot + 1] < 0 ? slotBonds(slot) / 2 : std::min(slotBonds(slot), slotBonds(slot + 1));
            if (!labelGroups.complete[label] && pairs > 0)
                pairs = std::min(pairs, labelGroups.bound(label, queryLabels_.fineCount, targetLabels_.fineCount));
            bonds += pairs;
        }
        for (int symbol = 0; symbol < static_cast<int>(queryLabels_.bySymbol.size()); ++symbol) {
            int pairs = std::min(queryLabels_.bySymbol[symbol], targetLabels_.bySymbol[symbol]);
            if (!atomGroups.complete[symbol] && pairs > 0)
                pairs = atomGroups.bound(symbol, queryLabels_.kindCount, targetLabels_.kindCount);
            atoms += pairs;
            queryLabels_.bySymbol[symbol] = targetLabels_.bySymbol[symbol] = 0;
        }
    }

    /**
     * Adds the new bonds and atoms that the atoms listed by countAtoms() can pair off: the bonds
     * of each label listed and the atoms of each class listed. Clears the counts.
     */
    void pairCounted(int& bonds, int& atoms) {
        Groups& labelGroups = graph_.labelGroups;
        Groups& atomGroups = graph_.atomGroups;
        clock_.tick(int64_t(countedLabels_.size()) * histogramWidth_ + int64_t(countedSymbols_.size()));
        for (int label : countedLabels_) {
            const int slot = 2 * label;
            int pairs = slotSymbol_[slot + 1] < 0 ? slotBonds(slot) / 2 : std::min(slotBonds(slot), slotBonds(slot + 1));
            if (!labelGroups.complete[label] && pairs > 0)
                pairs = std::min(pairs, labelGroups.bound(label, queryLabels_.fineCount, targetLabels_.fineCount));
            bonds += pairs;
            labelListed_[label] = 0;
        }
        for (int symbol : countedSymbols_) {
            int pairs = std::min(queryLabels_.bySymbol[symbol], targetLabels_.bySymbol[symbol]);
            if (!atomGroups.complete[symbol] && pairs > 0)
                pairs = atomGroups.bound(symbol, queryLabels_.kindCount, targetLabels_.kindCount);
            atoms += pairs;
            queryLabels_.bySymbol[symbol] = targetLabels_.bySymbol[symbol] = 0;
            symbolListed_[symbol] = 0;
        }
        countedLabels_.clear();
        countedSymbols_.clear();
        queryLabels_.clearCounts();
        targetLabels_.clearCounts();
    }

    /** The most new bonds of a slot, pairing the open atoms with the most bonds of it on either side. This clears the counts of the slot. */
    int slotBonds(int slot) {
        int sum = 0;
        for (int openBonds = histogramWidth_ - 1, queryAtoms = 0, targetAtoms = 0; openBonds > 0; --openBonds) {
            const size_t index = size_t(slot) * histogramWidth_ + openBonds;
            queryAtoms += queryLabels_.histogram[index];
            targetAtoms += targetLabels_.histogram[index];
            queryLabels_.histogram[index] = 0;
            targetLabels_.histogram[index] = 0;
            sum += std::min(queryAtoms, targetAtoms);
        }
        return sum;
    }

    /** The slot at each place, or -1 for none: twice the label of its bond, plus one at the second class. */
    std::vector<int> slots(const Side& side) const {
        std::vector<int> slot(side.neighbour.size(), -1);
        for (int atom = 0; atom < side.atomCount; ++atom)
            for (int place = side.start[atom]; place < side.start[atom + 1]; ++place) {
                const int label = side.label[side.neighbourBond[place]];
                slot[place] = label < 0 ? -1 : 2 * label + (side.atomClass[atom] == slotSymbol_[2 * label] ? 0 : 1);
            }
        return slot;
    }
};

} // namespace common_bond_detail
} // namespace smsd

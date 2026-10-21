/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/bitops.hpp"
#include "smsd/common_bond_search.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

/*
 * The single best maximum common induced subgraph (MCIS) of two molecules in any number of
 * fragments: the mapping from the query to the target with the most atoms and then the most
 * bonds, in which two mapped query atoms are bonded exactly when their images are, and then by
 * compatible bonds. Each fragment has at least one or two atoms. The molecules are held as the
 * int graphs of a CommonBondGraphs, with the atom and bond predicates of the caller.
 *
 * The unmapped atoms that may still be mapped are kept in blocks: a block holds the query and
 * target atoms of one atom class with the same bonds, by bond class, to the mapped atoms. An atom
 * can be mapped only onto an atom of its own block, so the atoms a mapping can still gain are at
 * most those each block can pair off, which bounds a state. Mapping a pair splits each block by
 * the bonds to the new pair, and an atom whose block has nothing left on the other side can no
 * longer be mapped. The blocks and their bound are the label classes of C. McCreesh, P. Prosser
 * and J. Trimble, "A Partitioning Algorithm for Maximum Common Subgraph Problems", Proc. IJCAI
 * 2017, 712-719, with atom classes as vertex labels and bond classes as edge labels. A depth-first
 * branch and bound grows the mapping one fragment at a time on an
 * explicit stack of levels, one for each mapped pair. A level maps a query atom of the smallest
 * block onto each target of the block in turn, and then leaves that atom unmapped for the rest of
 * the level, so each mapping is found once. While a fragment is open, the atom is one bonded to
 * it; once none is left, the fragment is closed and the atom is one bonded to no mapped atom, the
 * root of a new fragment. A pair that some best extension of a state maps is the only one tried.
 * Greedy seeds give a score to beat. The best score is first sought level by level of common
 * atoms, down from an upper bound, while the levels grow in cost; one search from the best mapping
 * known then finishes the work. Every loop ticks a Clock, which throws Stop when the time limit is
 * reached or the search is cancelled; the best mapping found by then is kept. The Java
 * CommonAtomSearch is the same search and gives the same mappings.
 */
namespace smsd {
namespace common_atom_detail {

using common_bond_detail::Clock;
using common_bond_detail::CommonBondGraphs;
using common_bond_detail::RADII;
using common_bond_detail::Side;
using common_bond_detail::Stop;

/** A score of this search, atoms << 32 | bonds, which compares by atoms first. */
inline int64_t score(int atoms, int bonds) { return int64_t((uint64_t(uint32_t(atoms)) << 32) | uint32_t(bonds)); }
/** The atoms of a score. */
inline int atomsOf(int64_t value) { return static_cast<int>(uint64_t(value) >> 32); }
/** The bonds of a score. */
inline int bondsOf(int64_t value) { return static_cast<int>(uint32_t(uint64_t(value))); }

/** The MCIS search on the int graphs of a CommonBondGraphs; see the description at the top of this file. */
class CommonAtomSearch {
    /** The most greedy seeds grown in each round. */
    static constexpr int SEEDS = 64;
    /** The most free targets a seed looks at next to a hub. */
    static constexpr int HUB_TARGETS = 64;
    /** The levels of likeness of two atoms, from no radius alike to all of them. */
    static constexpr int LIKENESS_LEVELS = RADII + 1;
    /** The most pairs of a linked block that blockComplete() checks one by one. */
    static constexpr int64_t COMPLETE_PAIRS = 64;
    /** The least growth in cost from one level of common atoms to the next that keeps the levels going. */
    static constexpr int64_t LEVEL_GROWTH = 2;

    /**
     * A mapped pair of the branch and bound, one for each depth, or the base level with nothing
     * mapped, with what it restores when taken off and what it has left to try.
     */
    struct Level {
        /** The mapped pair, -1 at the base level. */
        int queryAtom = -1, targetAtom = -1;
        /** The level of the root of the open fragment, or one at the base level, where no fragment is open. */
        int fragmentStart = 1;
        /** The bonds between mapped atoms before the pair was mapped, put back when it is taken off. */
        int commonBondsBefore = 0;
        /** The query and target bonds with no end mapped, once the pair is mapped. */
        int freeQueryBonds = 0, freeTargetBonds = 0;
        /** The blocks of the level, from blocksFrom up to blocksTo. */
        int blocksFrom = 0, blocksTo = 0;
        /** A score that no mapping below the level beats. */
        int64_t bound = 0;
        /** A number of its own, which tells the candidates tried at this level. */
        uint64_t serial = 0;
        /** The candidate being tried, or -1, and its block. */
        int triedAtom = -1, triedBlock = 0;
        /** Whether the candidate roots a new fragment, and whether its pair is forced. */
        bool triedRoot = false, forced = false;
        /** The targets of the candidate, the first targetCount, and the place of the next to try. */
        std::vector<int> targets;
        int targetCount = 0, nextTarget = 0;
    };

    /**
     * The pool of blocks of every level on the stack, those of each level after those of its
     * parent. A block holds runs of the query and of the target atoms: the live atoms of each side
     * from its start, and as many as its count. Its links are the mapped atoms its atoms are
     * bonded to, it is open when these include an atom of the open fragment, and all its atoms have
     * one atom class. A level that leaves an atom unmapped moves it after the live atoms of its
     * block, so the runs of its parent keep the same atoms. The blocks of a level are in the order
     * of their runs on either side. What depends only on the atoms of a block, its exact pairs and
     * whether all its pairs match, is kept once known, -1 while unknown.
     */
    struct Blocks {
        std::vector<int> queryStart, queryCount, targetStart, targetCount, links, atomClass, exactPairs, complete;
        std::vector<uint8_t> open;
        int count = 0;

        void add(int queryFrom, int queryAtoms, int targetFrom, int targetAtoms, int linked, bool isOpen, int ofClass) {
            if (count == static_cast<int>(queryStart.size())) {
                const size_t capacity = std::max<size_t>(64, 2 * size_t(count));
                for (auto* column : {&queryStart, &queryCount, &targetStart, &targetCount, &links, &atomClass,
                                     &exactPairs, &complete})
                    column->resize(capacity);
                open.resize(capacity);
            }
            queryStart[count] = queryFrom;
            queryCount[count] = queryAtoms;
            targetStart[count] = targetFrom;
            targetCount[count] = targetAtoms;
            links[count] = linked;
            open[count] = isOpen ? 1 : 0;
            atomClass[count] = ofClass;
            exactPairs[count] = complete[count] = -1;
            ++count;
        }
    };

    CommonBondGraphs& graph_;
    const Side& query_;
    const Side& target_;
    Clock& clock_;
    /** The fewest atoms of a fragment, one or two. */
    const int minAtoms_;
    const int queryCount_, targetCount_, words_;
    /** The compatible atom pairs: a row for each query atom, with a bit for each target atom. */
    const std::vector<uint64_t>& compatible_;
    /** Whether every query bond matches every target bond of its class, so that blocks need no bond check. */
    bool bondClassesMatch_ = false;
    /** The mapping both ways, -1 where unmapped. */
    std::vector<int> queryToTarget_, targetToQuery_;
    /** The levels: the first holds nothing mapped, and level d the d-th pair mapped. */
    std::vector<Level> levels_;
    /** The number of mapped pairs, which is the index of the top level, and the bonds between them. */
    int depth_ = 0, commonBonds_ = 0;
    /** The serial of the last level set up, and of the level at which each query atom was last tried. */
    uint64_t serial_ = 0;
    std::vector<uint64_t> triedSerial_;
    /** The atoms of the blocks, each block a run of each array, and the place of each atom. */
    std::vector<int> queryAtoms_, targetAtoms_, queryPlace_, targetPlace_;
    Blocks blocks_;
    /** The class of the bond from the pair being mapped to each atom next to it, valid where the stamp is current. */
    std::vector<int> queryLinkClass_, targetLinkClass_;
    std::vector<uint64_t> queryLinkStamp_, targetLinkStamp_;
    uint64_t linkStamp_ = 0;
    /** The live atoms bonded to the new pair on each side, and their blocks. */
    std::vector<int> linkedQuery_, linkedQueryBlock_, linkedTarget_, linkedTargetBlock_;
    /** The atoms of live blocks by stamp, and the live block of each, for liveBonds(). */
    std::vector<uint64_t> queryLive_, targetLive_;
    std::vector<int> queryBlock_, targetBlock_;
    uint64_t liveStamp_ = 0;
    std::vector<int> queryPairCount_, targetPairCount_, touchedPairs_;
    /** The atoms of each kind in a block, for blockPairs(). */
    std::vector<int> queryKindCount_, targetKindCount_;
    /** The targets of each query atom, the query atom of each target and the targets an augmenting path passed, for exactPairs(). */
    std::array<uint64_t, 64> matchRows_{};
    std::array<int, 64> matchOwner_{};
    uint64_t matchSeen_ = 0;
    /** The forced pair found by forcedPair(), and its block. */
    int forcedBlock_ = 0, forcedQuery_ = 0, forcedTarget_ = 0;
    /** The atoms by rank, lowest first, and the rank of each atom. */
    std::vector<int> queryOrder_, queryRank_, targetOrder_, targetRank_;
    /** The cell of each atom: atoms that look alike share a cell. */
    std::vector<int> queryCell_, targetCell_;
    /**
     * The runs of target atoms of one cell in target order, from targetRunStart_[run] up to the
     * next start, and in each the first place whose atom may still be a seed root.
     */
    std::vector<int> targetRunStart_, targetRunFree_;
    /** The seedRoots() call in which each query cell had no run left to take a root from. */
    std::vector<uint64_t> queryCellClosed_;
    uint64_t rootCall_ = 0;
    /** The signatures of the atoms at each radius, radius * atomCount + atom. */
    std::vector<uint64_t> querySignature_, targetSignature_;
    /** The next atom alike to each terminal atom, and the first atom of each ring of target twins. */
    std::vector<int> queryTwin_, targetTwin_, targetTwinFirst_;
    bool twinsKnown_ = false;
    /** The stamp of each ring of twins seen and the sort keys of the targets, for targets(). */
    std::vector<uint64_t> twinSeen_;
    uint64_t twinStamp_ = 0;
    std::vector<int64_t> targetKeys_;
    /** The taken targets that lead the neighbours of each target atom, and the seed they were counted for. */
    std::vector<int> leadingTaken_, leadingRound_;
    /** The likest level of each query atom with any target, -1 until worked out, see likestLevel(). */
    std::vector<int> likestLevel_;
    /** The best score so far, and a score that no mapping beats. */
    int64_t bestScore_ = 0, upperBound_ = 0;
    /** The mapping with the best score so far, while one beats the threshold of the pass. */
    std::vector<int> best_;
    bool hasBest_ = false;
    /**
     * The mapping a stopped search falls back on while the pass under way holds none: the seed,
     * or the best mapping a level of the search met below its threshold.
     */
    std::vector<int> fallback_;
    bool hasFallback_ = false;
    int64_t fallbackScore_ = 0;
    /** The pairs mapped. */
    int64_t nodes_ = 0;

public:
    /** A search of `graph` for fragments of at least `minAtoms` atoms, one or two. */
    CommonAtomSearch(CommonBondGraphs& graph, Clock& clock, int minAtoms)
        : graph_(graph), query_(graph.query), target_(graph.target), clock_(clock), minAtoms_(std::max(1, minAtoms)),
          queryCount_(graph.query.atomCount), targetCount_(graph.target.atomCount),
          words_(common_bond_detail::words(targetCount_)), compatible_(graph.compatiblePairs()),
          queryToTarget_(queryCount_, -1), targetToQuery_(targetCount_, -1),
          levels_(size_t(std::min(queryCount_, targetCount_)) + 1), triedSerial_(queryCount_, 0),
          queryAtoms_(queryCount_, 0), targetAtoms_(targetCount_, 0), queryPlace_(queryCount_, 0),
          targetPlace_(targetCount_, 0), queryLinkClass_(queryCount_, 0), targetLinkClass_(targetCount_, 0),
          queryLinkStamp_(queryCount_, 0), targetLinkStamp_(targetCount_, 0),
          linkedQuery_(size_t(graph.query.maxDegree()), 0), linkedQueryBlock_(size_t(graph.query.maxDegree()), 0),
          linkedTarget_(size_t(graph.target.maxDegree()), 0), linkedTargetBlock_(size_t(graph.target.maxDegree()), 0),
          queryLive_(queryCount_, 0), targetLive_(targetCount_, 0), queryBlock_(queryCount_, 0),
          targetBlock_(targetCount_, 0), queryKindCount_(graph.atomGroups.queryGroup.size(), 0),
          targetKindCount_(graph.atomGroups.targetGroup.size(), 0), twinSeen_(targetCount_, 0),
          targetKeys_(targetCount_, 0), leadingTaken_(targetCount_, 0), leadingRound_(targetCount_, 0),
          likestLevel_(queryCount_, -1) {
        bondClassesMatch_ = bondClassesMatch();
    }

    /** The pairs mapped so far. */
    int64_t nodes() const { return nodes_; }

    /**
     * A maximum mapping, or empty when there is none, as no pair is compatible or, with
     * fragments of two atoms, no bond can be common. The atoms are ranked, a greedy seed gives
     * the score to start from, and the branch and bound finds the mapping.
     */
    std::vector<int> run() {
        clock_.poll();
        if (!anyCompatible()) return {};
        rank();
        const int64_t bound = bestPossible();
        fallbackScore_ = belowEveryMapping();
        seed(bound);
        return collectByAtoms(bound);
    }

    /**
     * The best mapping found so far, for a search the clock stopped: the best of the pass under
     * way or, while it holds none, the one it falls back on. Empty when nothing was found. Each is
     * complete, as a mapping is kept only once it is recorded whole.
     */
    std::vector<int> found() const {
        if (hasBest_) return best_;
        if (hasFallback_) return fallback_;
        return {};
    }

private:
    bool has(const std::vector<uint64_t>& rows, int row, int bit) const { return common_bond_detail::hasBit(rows, words_, row, bit); }
    void set(std::vector<uint64_t>& rows, int row, int bit) const { common_bond_detail::setBit(rows, words_, row, bit); }
    void clear(std::vector<uint64_t>& rows, int row, int bit) const { common_bond_detail::clearBit(rows, words_, row, bit); }

    /** A pair as one int, queryAtom * targetCount + targetAtom. */
    int pair(int queryAtom, int targetAtom) const { return queryAtom * targetCount_ + targetAtom; }
    int pairQuery(int packed) const { return packed / targetCount_; }
    int pairTarget(int packed) const { return packed % targetCount_; }

    /** A score below that of every mapping, as every mapping has an atom. */
    static int64_t belowEveryMapping() { return score(0, 0); }

    bool anyCompatible() const {
        for (uint64_t word : compatible_)
            if (word != 0) return true;
        return false;
    }

    /**
     * A maximum mapping below `bound`. The best score is sought first level by level of common
     * atoms, down from the bound to the score of the best mapping known: each level asks for one
     * mapping with at least that many atoms, and without one the bound comes down. Each level
     * down repeats the work of the one above it, which pays only while each level costs at least
     * LEVEL_GROWTH times the one above it. Once a level costs less, one search from the best
     * mapping known finishes the work, its threshold rising as better mappings are found.
     */
    std::vector<int> collectByAtoms(int64_t bound) {
        int64_t previous = 0;
        for (int triedAtoms = atomsOf(bound); triedAtoms > atomsOf(fallbackScore_); --triedAtoms) {
            const int64_t before = nodes_;
            collect(nullptr, score(triedAtoms, 0) - 1, bound);
            if (hasBest_) return best_;
            // No mapping has that many atoms.
            bound = score(triedAtoms - 1, bondsOf(bound));
            const int64_t cost = nodes_ - before;
            if (previous > 0 && cost < LEVEL_GROWTH * previous) break;
            previous = std::min<int64_t>(INT_MAX, std::max<int64_t>(1, cost));
        }
        collect(hasFallback_ ? &fallback_ : nullptr, fallbackScore_, bound);
        return hasBest_ ? best_ : std::vector<int>{};
    }

    /**
     * A score that no mapping beats. Every cap on the atoms or on the bonds holds for every
     * mapping, so the two are capped apart; every bond between mapped atoms is a common bond.
     */
    int64_t bestPossible() {
        const int64_t label = graph_.labelBound(false);
        const int atomCap = std::min(std::min(queryCount_, targetCount_), common_bond_detail::atomsOf(label));
        const int bondCap = std::min(std::min(query_.bondCount, target_.bondCount), common_bond_detail::bondsOf(label));
        // A mapping with at most so many atoms has at most so many bonds between them.
        const int64_t pairs = int64_t(atomCap) * (atomCap - 1) / 2;
        return score(atomCap, static_cast<int>(std::min<int64_t>(bondCap, pairs)));
    }

    /**
     * Whether every query bond matches every target bond of its class. The bond predicate gives
     * bonds of one kind the same pairs, so a bond of each kind stands for the others.
     */
    bool bondClassesMatch() {
        int queryKinds = 0, targetKinds = 0;
        for (int kind : query_.bondKind) queryKinds = std::max(queryKinds, kind + 1);
        for (int kind : target_.bondKind) targetKinds = std::max(targetKinds, kind + 1);
        std::vector<int> queryExample(queryKinds, -1), targetExample(targetKinds, -1);
        for (int bond = query_.bondCount - 1; bond >= 0; --bond) queryExample[query_.bondKind[bond]] = bond;
        for (int bond = target_.bondCount - 1; bond >= 0; --bond) targetExample[target_.bondKind[bond]] = bond;
        for (int queryBond : queryExample) {
            clock_.tick(1 + targetKinds);
            for (int targetBond : targetExample)
                if (queryBond >= 0 && targetBond >= 0 && query_.bondClass[queryBond] >= 0
                    && query_.bondClass[queryBond] == target_.bondClass[targetBond]
                    && !graph_.bondsMatch(queryBond, targetBond)) return false;
        }
        return true;
    }

    /**
     * Ranks the atoms and works out their signatures. Query atoms with the fewest candidates,
     * then the most bonds, rank highest; ties are broken by colour refinement, so that the order
     * is the same for the same graph. Target atoms are ranked by their cell only, in input order
     * within it. The ranks order the seeds and break ties between candidates.
     */
    void rank() {
        querySignature_ = query_.signatures(clock_);
        targetSignature_ = target_.signatures(clock_);
        std::vector<int> queryChoices(queryCount_, 0), targetChoices(targetCount_, 0);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
            const uint64_t* row = compatible_.data() + size_t(queryAtom) * words_;
            for (int targetAtom = common_bond_detail::nextBit(row, words_, 0); targetAtom >= 0;
                 targetAtom = common_bond_detail::nextBit(row, words_, targetAtom + 1)) {
                clock_.tick();
                ++queryChoices[queryAtom];
                ++targetChoices[targetAtom];
            }
        }
        const std::vector<int> colours = query_.invariants(queryChoices);
        queryCell_ = common_bond_detail::refine(query_.start, query_.neighbour, colours, false, clock_);
        const std::vector<int> canon = common_bond_detail::refine(query_.start, query_.neighbour, colours, true, clock_);
        // Ties among target atoms stay in input order. Breaking them one atom at a time costs a
        // round per atom on a symmetric target, and the order only has to be repeatable.
        targetCell_ = common_bond_detail::refine(target_.start, target_.neighbour, target_.invariants(targetChoices), false, clock_);
        std::vector<int> next(targetCount_, 0);
        targetRank_.assign(targetCount_, 0);
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom) next[targetCell_[targetAtom]] = targetCell_[targetAtom];
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom) targetRank_[targetAtom] = next[targetCell_[targetAtom]]++;
        targetOrder_ = common_bond_detail::inverse(targetRank_);
        targetRunStart_.clear();
        for (int position = 0; position < targetCount_; ++position)
            if (position == 0 || targetCell_[targetOrder_[position]] != targetCell_[targetOrder_[position - 1]])
                targetRunStart_.push_back(position);
        targetRunStart_.push_back(targetCount_);
        std::vector<int> byRank(queryCount_);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) byRank[queryAtom] = queryAtom;
        // The most choices come first, then the fewest bonds, then the lowest canonical rank.
        std::stable_sort(byRank.begin(), byRank.end(), [&](int left, int right) {
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
    }

    /**
     * The best greedy seed, if any, as the fallback: a set of fragments grown one after another,
     * each on the atoms bonded to no mapped atom on either side. A fragment is grown from a root
     * pair, from each of at most SEEDS of them or until one reaches `bound`, and the best is
     * kept. The fallback follows each better seed, so a search stopped while seeding still has
     * the best seed grown by then, or the seed being grown if that is better. It needs rank()
     * first.
     */
    void seed(int64_t bound) {
        int keptAtoms = 0;
        try {
            growSeeds(bound, keptAtoms);
        } catch (const Stop&) {
            // The fragments mapped when the clock stopped are induced, as each pair added kept them so.
            int mappedAtoms = 0;
            for (int image : queryToTarget_) mappedAtoms += image >= 0;
            if (mappedAtoms - keptAtoms >= minAtoms_ && scoreOf(queryToTarget_) > fallbackScore_) {
                fallback_ = queryToTarget_;
                hasFallback_ = true;
                fallbackScore_ = scoreOf(fallback_);
            }
            throw;
        }
        std::fill(queryToTarget_.begin(), queryToTarget_.end(), -1);
        std::fill(targetToQuery_.begin(), targetToQuery_.end(), -1);
    }

    /** The rounds of seed(), with the atoms of the fragments kept from earlier rounds in `keptAtoms`. */
    void growSeeds(int64_t bound, int& keptAtoms) {
        std::vector<uint64_t> taken(size_t(queryCount_) * size_t(words_), 0);
        std::vector<int> roots(size_t(LIKENESS_LEVELS) * SEEDS, 0), rootCount(LIKENESS_LEVELS, 0);
        targetRunFree_.assign(targetRunStart_.begin(), targetRunStart_.end() - 1);
        queryCellClosed_.assign(queryCount_, 0);
        std::vector<int> frontier(queryCount_), merged(queryCount_), mapped(queryCount_), added(queryCount_, 0);
        int64_t bestSeed = belowEveryMapping();
        for (int round = 0, keptBonds = 0;;) {
            const int64_t before = bestSeed;
            seedRoots(roots, rootCount, taken);
            for (int level = RADII, tried = 0; level >= 0 && bestSeed != bound; --level)
                for (int rootIndex = 0; rootIndex < rootCount[level] && tried < SEEDS && bestSeed != bound; ++rootIndex, ++tried) {
                    const int64_t grown =
                        growSeed(roots[size_t(level) * SEEDS + rootIndex], frontier, merged, mapped, added, ++round);
                    const int seedAtoms = atomsOf(grown), seedBonds = bondsOf(grown);
                    if (seedAtoms >= minAtoms_ && score(keptAtoms + seedAtoms, keptBonds + seedBonds) > bestSeed) {
                        bestSeed = score(keptAtoms + seedAtoms, keptBonds + seedBonds);
                        fallback_ = queryToTarget_;
                        hasFallback_ = true;
                        fallbackScore_ = scoreOf(fallback_);
                    }
                    for (int mappedIndex = 0; mappedIndex < seedAtoms; ++mappedIndex) {
                        targetToQuery_[queryToTarget_[mapped[mappedIndex]]] = -1;
                        queryToTarget_[mapped[mappedIndex]] = -1;
                    }
                }
            if (bestSeed == before || bestSeed == bound) break;
            // The fragments found so far stay mapped while the next one grows from atoms bonded to
            // none of them. Each pair added keeps the mapping induced, so the fragments stay valid.
            keptAtoms = atomsOf(bestSeed);
            keptBonds = bondsOf(bestSeed);
            queryToTarget_ = fallback_;
            for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
                if (queryToTarget_[queryAtom] >= 0) targetToQuery_[queryToTarget_[queryAtom]] = queryAtom;
        }
    }

    /**
     * Picks the seed roots by how alike their atoms are, among the free atoms bonded to no mapped
     * atom. The roots of each level are at roots[level * SEEDS]. Alike cells grow alike seeds, so
     * there is one root for each pair of cells, marked in `taken`. Each query atom takes, from
     * each run of target cells in target order, the first free compatible target bonded to no
     * mapped atom. Atoms of one cell look alike, so a run whose level is full is passed whole, and
     * a query cell with no run left is passed for the rest of the call. The kept fragments only
     * grow between calls, so a target that cannot be a root stays so and the start of each run
     * moves past it once.
     */
    void seedRoots(std::vector<int>& roots, std::vector<int>& rootCount, std::vector<uint64_t>& taken) {
        // Only the cell pairs of the last call's roots are taken.
        for (int level = 0; level < LIKENESS_LEVELS; ++level)
            for (int index = 0; index < rootCount[level]; ++index) {
                const int packed = roots[size_t(level) * SEEDS + index];
                clear(taken, queryCell_[pairQuery(packed)], targetCell_[pairTarget(packed)]);
            }
        std::fill(rootCount.begin(), rootCount.end(), 0);
        ++rootCall_;
        const int likest = RADII, runs = static_cast<int>(targetRunStart_.size()) - 1;
        // The most alike roots are tried first. Once there are enough of them, the rest are not needed.
        for (int rankIndex = queryCount_ - 1; rankIndex >= 0 && rootCount[likest] < SEEDS; --rankIndex) {
            const int queryAtom = queryOrder_[rankIndex], queryCell = queryCell_[queryAtom];
            clock_.tick();
            if (queryCellClosed_[queryCell] == rootCall_ || queryToTarget_[queryAtom] >= 0
                || mappedNeighbours(query_, queryToTarget_, queryAtom) > 0) continue;
            bool closed = true;
            for (int run = 0; run < runs && rootCount[likest] < SEEDS; ++run) {
                int& first = targetRunFree_[run];
                const int end = targetRunStart_[run + 1];
                clock_.tick();
                while (first < end && !rootTarget(targetOrder_[first])) ++first;
                if (first == end || has(taken, queryCell, targetCell_[targetOrder_[first]])) continue;
                const int level = likeness(queryAtom, targetOrder_[first]);
                if (rootCount[level] >= SEEDS) continue;
                bool rooted = false;
                for (int position = first; position < end && !rooted; ++position) {
                    const int targetAtom = targetOrder_[position];
                    clock_.tick();
                    if (!has(compatible_, queryAtom, targetAtom) || !rootTarget(targetAtom)) continue;
                    set(taken, queryCell, targetCell_[targetAtom]);
                    roots[size_t(level) * SEEDS + rootCount[level]++] = pair(queryAtom, targetAtom);
                    rooted = true;
                }
                closed = closed && rooted;
            }
            if (closed) queryCellClosed_[queryCell] = rootCall_;
        }
    }

    /** Whether a target atom is free and bonded to no mapped atom, so that it may root a seed. */
    bool rootTarget(int targetAtom) {
        return targetToQuery_[targetAtom] < 0 && mappedNeighbours(target_, targetToQuery_, targetAtom) == 0;
    }

    /** The number of mapped atoms bonded to `atom` of `side`, with `image` its mapping. */
    int mappedNeighbours(const Side& side, const std::vector<int>& image, int atom) {
        int count = 0;
        clock_.tick(1 + side.degree(atom));
        for (int place = side.start[atom]; place < side.start[atom + 1]; ++place)
            if (image[side.neighbour[place]] >= 0) ++count;
        return count;
    }

    /**
     * Grows a seed from the root pair, each time by the best frontier pair. The frontier is the
     * free query atoms next to the seed, each listed once by marking it with `round` in `added`,
     * highest rank first; `merged` is scratch for it. The best pair is the same in any order, and
     * in this one bestSeedPair() soon passes over the atoms that cannot beat the pair it holds.
     * The seed is left mapped, with its query atoms in `mapped` in order. The result is its score.
     */
    int64_t growSeed(int root, std::vector<int>& frontier, std::vector<int>& merged, std::vector<int>& mapped,
                     std::vector<int>& added, int round) {
        int count = 0, frontierSize = 0, seedBonds = 0;
        for (int64_t packed = root; packed >= 0;) {
            const int queryAtom = pairQuery(static_cast<int>(packed & 0xffffffff));
            const int targetAtom = pairTarget(static_cast<int>(packed & 0xffffffff));
            queryToTarget_[queryAtom] = targetAtom;
            targetToQuery_[targetAtom] = queryAtom;
            mapped[count++] = queryAtom;
            seedBonds += static_cast<int>(packed >> 32);
            const int listed = frontierSize;
            for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place) {
                const int queryNeighbour = query_.neighbour[place];
                if (queryToTarget_[queryNeighbour] < 0 && added[queryNeighbour] != round) {
                    added[queryNeighbour] = round;
                    frontier[frontierSize++] = queryNeighbour;
                }
            }
            common_bond_detail::mergeByRank(frontier, listed, frontierSize, queryRank_, merged, clock_);
            int64_t held = -1;
            int kept = 0;
            for (int frontierIndex = 0; frontierIndex < frontierSize; ++frontierIndex) {
                const int frontierAtom = frontier[frontierIndex];
                // The atom just mapped leaves the frontier.
                if (queryToTarget_[frontierAtom] >= 0) continue;
                frontier[kept++] = frontierAtom;
                held = bestSeedPair(frontierAtom, held, round);
            }
            frontierSize = kept;
            packed = held;
        }
        return score(count, seedBonds);
    }

    /**
     * The best pair of `frontierAtom` with a free target next to the image of one of its mapped
     * neighbours, unless the pair `held` beats it. The pair must keep the mapping induced, see
     * inducedLinks(). The most alike pair wins, then the one with the most new bonds; ties go to
     * the higher ranked query atom, then to the lower ranked target atom. The pair is returned in
     * the low 32 bits, with its new bonds above, or -1 if there is none; `held` is read the same
     * way. Next to a hub, the first few free targets are as good for a seed as all of them. The
     * seed of `round` only takes targets, so the taken targets that lead the list of an image
     * are passed once for the seed, see leadingTaken(), and not again for each frontier atom.
     */
    int64_t bestSeedPair(int frontierAtom, int64_t held, int round) {
        int packed = held < 0 ? -1 : static_cast<int>(held & 0xffffffff);
        int links = held < 0 ? 0 : static_cast<int>(held >> 32);
        int likest = packed < 0 ? -1 : likeness(pairQuery(packed), pairTarget(packed));
        clock_.tick(1 + query_.degree(frontierAtom));
        if (packed >= 0) {
            // No pair of the atom is liker than its likest level, nor has more new bonds than the
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
                if (targetToQuery_[candidate] >= 0 || !has(compatible_, frontierAtom, candidate)) continue;
                ++freeCount;
                const int level = likeness(frontierAtom, candidate);
                if (level < likest) continue;
                const int newBonds = inducedLinks(frontierAtom, candidate);
                if (newBonds < 0) continue;
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

    /**
     * The bonds that mapping `queryAtom` onto `targetAtom` adds, or -1 if the mapping would not
     * stay induced: the target atom must be bonded to the image of each mapped neighbour of the
     * query atom by a compatible bond, and to no other mapped atom. As the images of distinct
     * neighbours are distinct, the second holds when both atoms have as many mapped neighbours.
     */
    int inducedLinks(int queryAtom, int targetAtom) {
        int queryLinks = 0;
        clock_.tick(1 + int64_t(query_.degree(queryAtom)) * (1 + target_.degree(targetAtom)));
        for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
            const int neighbourImage = queryToTarget_[query_.neighbour[queryPlace]];
            if (neighbourImage < 0) continue;
            ++queryLinks;
            const int targetBond = target_.bondBetween(targetAtom, neighbourImage);
            if (targetBond < 0 || !graph_.bondsMatch(query_.neighbourBond[queryPlace], targetBond)) return -1;
        }
        return queryLinks == mappedNeighbours(target_, targetToQuery_, targetAtom) ? queryLinks : -1;
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
     * Searches for the best mapping that beats the score `start`, which `seed` holds unless it
     * is null; the result is best_, or nothing when none does and there is no seed. The search
     * ends early once it holds a mapping that reaches `upperBound`. It needs rank() first.
     */
    void collect(const std::vector<int>* seed, int64_t start, int64_t upperBound) {
        bestScore_ = start;
        upperBound_ = upperBound;
        hasBest_ = seed != nullptr;
        if (seed) best_ = *seed;
        // The twins are found once, before any search.
        if (!twinsKnown_) {
            queryTwin_ = twins(query_, true);
            targetTwin_ = twins(target_, false);
            targetTwinFirst_ = ringFirsts(targetTwin_);
            twinsKnown_ = true;
        }
        Level& base = levels_[0];
        depth_ = 0;
        commonBonds_ = 0;
        base.queryAtom = base.targetAtom = -1;
        base.fragmentStart = 1;
        base.freeQueryBonds = query_.bondCount;
        base.freeTargetBonds = target_.bondCount;
        base.serial = ++serial_;
        base.triedAtom = -1;
        base.forced = false;
        base.targetCount = base.nextTarget = 0;
        base.blocksFrom = blocks_.count = 0;
        baseBlocks();
        base.blocksTo = blocks_.count;
        base.bound = bound(base);
        search();
    }

    /**
     * Sets up the blocks of the base level: for each atom class, its query and target atoms in
     * index order. The atoms that match nothing come after every block.
     */
    void baseBlocks() {
        const int classes = graph_.atomClasses;
        std::vector<int> queryStart(size_t(classes) + 1, 0), targetStart(size_t(classes) + 1, 0);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
            if (query_.atomClass[queryAtom] >= 0) ++queryStart[query_.atomClass[queryAtom] + 1];
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom)
            if (target_.atomClass[targetAtom] >= 0) ++targetStart[target_.atomClass[targetAtom] + 1];
        for (int atomClass = 0; atomClass < classes; ++atomClass) {
            queryStart[atomClass + 1] += queryStart[atomClass];
            targetStart[atomClass + 1] += targetStart[atomClass];
        }
        std::vector<int> queryNext(queryStart.begin(), queryStart.end() - 1), targetNext(targetStart.begin(), targetStart.end() - 1);
        int queryOutside = queryStart[classes], targetOutside = targetStart[classes];
        clock_.tick(int64_t(queryCount_) + targetCount_);
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
            const int atomClass = query_.atomClass[queryAtom];
            const int place = atomClass >= 0 ? queryNext[atomClass]++ : queryOutside++;
            queryAtoms_[place] = queryAtom;
            queryPlace_[queryAtom] = place;
        }
        for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom) {
            const int atomClass = target_.atomClass[targetAtom];
            const int place = atomClass >= 0 ? targetNext[atomClass]++ : targetOutside++;
            targetAtoms_[place] = targetAtom;
            targetPlace_[targetAtom] = place;
        }
        for (int atomClass = 0; atomClass < classes; ++atomClass) {
            const int queryAtomsOfClass = queryStart[atomClass + 1] - queryStart[atomClass];
            const int targetAtomsOfClass = targetStart[atomClass + 1] - targetStart[atomClass];
            if (queryAtomsOfClass > 0 && targetAtomsOfClass > 0)
                blocks_.add(queryStart[atomClass], queryAtomsOfClass, targetStart[atomClass], targetAtomsOfClass, 0, false, atomClass);
        }
    }

    /**
     * Searches depth first from the base level until it is taken off. A level tries the targets
     * of each of its candidates in turn, see nextCandidate(). It is taken off once nothing is
     * left or its bound is too low. The levels are an explicit stack, so large molecules do not
     * run out of stack.
     */
    void search() {
        while (depth_ >= 0) {
            clock_.tick();
            Level& level = levels_[depth_];
            if (done() || !canImprove(level.bound)) {
                pop();
            } else if (level.nextTarget < level.targetCount) {
                const int targetAtom = level.targets[level.nextTarget++];
                push(level.triedAtom, targetAtom, level.triedRoot);
            } else if (!nextCandidate(level)) {
                pop();
            }
        }
    }

    /**
     * Maps `queryAtom` onto `targetAtom` as a new level, which roots a new fragment when `root`
     * is set and otherwise joins the open one. The mapping is recorded, and the level gets its
     * blocks and its bound unless it is a leaf.
     */
    void push(int queryAtom, int targetAtom, bool root) {
        ++nodes_;
        const Level& parent = levels_[depth_];
        Level& level = levels_[depth_ + 1];
        level.queryAtom = queryAtom;
        level.targetAtom = targetAtom;
        level.fragmentStart = root ? depth_ + 1 : parent.fragmentStart;
        level.commonBondsBefore = commonBonds_;
        level.freeQueryBonds = parent.freeQueryBonds;
        level.freeTargetBonds = parent.freeTargetBonds;
        clock_.tick(1 + query_.degree(queryAtom) + target_.degree(targetAtom));
        for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
            if (queryToTarget_[query_.neighbour[queryPlace]] < 0) --level.freeQueryBonds;
            else ++commonBonds_;
        }
        for (int targetPlace = target_.start[targetAtom]; targetPlace < target_.start[targetAtom + 1]; ++targetPlace)
            if (targetToQuery_[target_.neighbour[targetPlace]] < 0) --level.freeTargetBonds;
        queryToTarget_[queryAtom] = targetAtom;
        targetToQuery_[targetAtom] = queryAtom;
        ++depth_;
        level.serial = ++serial_;
        level.triedAtom = -1;
        level.forced = false;
        level.targetCount = level.nextTarget = 0;
        level.blocksFrom = blocks_.count;
        record();
        if (isLeaf()) {
            level.blocksTo = blocks_.count;
            level.bound = -1;
            return;
        }
        refine(levels_[depth_ - 1], queryAtom, targetAtom, root);
        level.blocksTo = blocks_.count;
        level.bound = bound(level);
    }

    /** Takes the top level off, putting back the blocks and the mapping of its parent; the base level ends the search. */
    void pop() {
        const Level& level = levels_[depth_];
        blocks_.count = level.blocksFrom;
        if (depth_ > 0) {
            queryToTarget_[level.queryAtom] = -1;
            targetToQuery_[level.targetAtom] = -1;
            commonBonds_ = level.commonBondsBefore;
        }
        --depth_;
    }

    /**
     * Whether the mapping is a leaf, with nothing left to add: the pass can end, or the mapping
     * has the atoms of the upper bound or every atom of either molecule.
     */
    bool isLeaf() const {
        return done() || depth_ == atomsOf(upperBound_) || depth_ == queryCount_ || depth_ == targetCount_;
    }

    /**
     * Gives `level` its next candidate and the targets to try it on, or false when nothing is
     * left. The candidate tried so far is left unmapped below the level from now on, as every
     * mapping with it was tried, which can lower the bound. A forced pair, see forcedPair(), is
     * the only one the level then tries. While the open fragment has query atoms bonded to it, the
     * next candidate is one of them, see candidate(). Once there are none the fragment is closed
     * and, unless it is too small, the next candidate is a query atom bonded to no mapped atom,
     * the root of a new fragment. A candidate with a twin tried at this level is skipped.
     */
    bool nextCandidate(Level& level) {
        for (;;) {
            clock_.tick();
            if (level.triedAtom >= 0) {
                // A forced pair holds a best extension of the state, so nothing else is needed.
                if (level.forced) return false;
                leaveUnmapped(level.triedBlock, level.triedAtom);
                triedSerial_[level.triedAtom] = level.serial;
                level.triedAtom = -1;
                level.bound = bound(level);
                if (!canImprove(level.bound)) return false;
            }
            if (forcedPair(level)) {
                level.triedAtom = forcedQuery_;
                level.triedBlock = forcedBlock_;
                level.triedRoot = false;
                level.forced = true;
                level.nextTarget = 0;
                if (level.targets.empty()) level.targets.resize(1);
                level.targets[0] = forcedTarget_;
                level.targetCount = 1;
                return true;
            }
            int64_t chosen = depth_ > 0 ? candidate(level, true) : -1;
            const bool root = chosen < 0;
            if (root) {
                if (depth_ > 0 && depth_ - level.fragmentStart + 1 < minAtoms_) return false;
                chosen = candidate(level, false);
                if (chosen < 0) return false;
            }
            const int block = static_cast<int>(chosen >> 32), atom = static_cast<int>(chosen & 0xffffffff);
            level.triedAtom = atom;
            level.triedBlock = block;
            level.triedRoot = root;
            level.nextTarget = 0;
            if (twinTried(atom, level.serial)) {
                level.targetCount = 0;
                continue;
            }
            targets(level, block, atom);
            if (level.targetCount > 0) return true;
        }
    }

    /**
     * Whether the state of `level` has a forced pair, which some best extension of the state
     * maps, so that no other pair needs to be tried: a query and a target atom of one block,
     * neither bonded to an atom of a live block, so neither to an atom that may still be mapped,
     * where every pair of the block matches, see blockComplete(). Swapping such an atom with the
     * one an extension maps onto its partner, or adding the pair when neither is mapped, keeps the
     * extension induced without losing atoms or bonds. A pair bonded to no mapped atom forms a
     * fragment of one atom, so it is forced only when such fragments are allowed. The pair is, of
     * the first block that has one, the query atom of highest rank and the first such target of
     * the run.
     */
    bool forcedPair(const Level& level) {
        for (int block = level.blocksFrom; block < level.blocksTo; ++block) {
            clock_.tick();
            const int queryAtomCount = blocks_.queryCount[block], targetAtomCount = blocks_.targetCount[block];
            if (queryAtomCount == 0 || targetAtomCount == 0 || (blocks_.links[block] == 0 && minAtoms_ > 1)
                || !blockComplete(block)) continue;
            int chosen = -1;
            for (int place = blocks_.queryStart[block]; place < blocks_.queryStart[block] + queryAtomCount; ++place) {
                const int atom = queryAtoms_[place];
                clock_.tick();
                if ((chosen < 0 || queryRank_[atom] > queryRank_[chosen])
                    && !hasLiveNeighbour(level, query_, queryPlace_, atom, true) && !twinTried(atom, level.serial))
                    chosen = atom;
            }
            if (chosen < 0) continue;
            for (int place = blocks_.targetStart[block]; place < blocks_.targetStart[block] + targetAtomCount; ++place) {
                const int atom = targetAtoms_[place];
                clock_.tick();
                if (!hasLiveNeighbour(level, target_, targetPlace_, atom, false)) {
                    forcedBlock_ = block;
                    forcedQuery_ = chosen;
                    forcedTarget_ = atom;
                    return true;
                }
            }
        }
        return false;
    }

    /** Whether `atom` of `side`, with the atoms at `place`, is bonded to an atom of a live block of `level`, one with atoms on both sides. */
    bool hasLiveNeighbour(const Level& level, const Side& side, const std::vector<int>& place, int atom, bool onQuery) {
        for (int neighbourPlace = side.start[atom]; neighbourPlace < side.start[atom + 1]; ++neighbourPlace) {
            const int block = blockOf(level, place[side.neighbour[neighbourPlace]], onQuery);
            if (block >= 0 && blocks_.queryCount[block] > 0 && blocks_.targetCount[block] > 0) return true;
        }
        return false;
    }

    /**
     * The next candidate from the open blocks of `level`, which joins the open fragment, or from
     * the others, which roots a new fragment, as its block in the high 32 bits and the atom in the
     * low ones, or -1 if there is none. The block is the one with the fewest atoms on its larger
     * side, then the earliest; the atom is the one of its query atoms with the most bonds, then of
     * highest rank.
     */
    int64_t candidate(const Level& level, bool open) {
        int chosenBlock = -1, chosenSize = INT_MAX;
        for (int block = level.blocksFrom; block < level.blocksTo; ++block) {
            const int queryAtomCount = blocks_.queryCount[block], targetAtomCount = blocks_.targetCount[block];
            if ((blocks_.open[block] != 0) != open || queryAtomCount == 0 || targetAtomCount == 0) continue;
            const int size = std::max(queryAtomCount, targetAtomCount);
            if (size < chosenSize) {
                chosenBlock = block;
                chosenSize = size;
            }
        }
        clock_.tick(level.blocksTo - level.blocksFrom);
        if (chosenBlock < 0) return -1;
        const int start = blocks_.queryStart[chosenBlock];
        int chosen = -1;
        clock_.tick(blocks_.queryCount[chosenBlock]);
        for (int place = start; place < start + blocks_.queryCount[chosenBlock]; ++place) {
            const int atom = queryAtoms_[place];
            if (chosen < 0 || query_.degree(atom) > query_.degree(chosen)
                || (query_.degree(atom) == query_.degree(chosen) && queryRank_[atom] > queryRank_[chosen]))
                chosen = atom;
        }
        return (int64_t(chosenBlock) << 32) | int64_t(chosen);
    }

    /**
     * Whether a twin of `queryAtom` was tried at the level of serial `levelSerial`. Swapping the
     * two takes every mapping below the atom onto one of the same score below the twin, which was
     * searched before. The twin is then unmapped, and the atom is left so too.
     */
    bool twinTried(int queryAtom, uint64_t levelSerial) {
        for (int twin = queryTwin_[queryAtom]; twin >= 0 && twin != queryAtom; twin = queryTwin_[twin]) {
            clock_.tick();
            if (triedSerial_[twin] == levelSerial) return true;
        }
        return false;
    }

    /**
     * Lists in `level` the targets of `queryAtom`, which is in `block`: the target atoms of the
     * block it is compatible with, by compatible bonds to the images of its mapped neighbours. The
     * most alike come first, then by descending rank. Of target twins, only the first is kept, as
     * swapping two free twins takes every mapping onto one of the same score.
     */
    void targets(Level& level, int block, int queryAtom) {
        const int start = blocks_.targetStart[block], count = blocks_.targetCount[block];
        if (static_cast<int>(level.targets.size()) < count)
            level.targets.resize(std::max<size_t>(size_t(count), 2 * level.targets.size()));
        int listed = 0;
        clock_.tick(1 + count);
        for (int place = start; place < start + count; ++place) {
            const int targetAtom = targetAtoms_[place];
            if (has(compatible_, queryAtom, targetAtom) && (bondClassesMatch_ || bondsFit(queryAtom, targetAtom)))
                targetKeys_[listed++] = targetKey(queryAtom, targetAtom);
        }
        std::sort(targetKeys_.begin(), targetKeys_.begin() + listed);
        const uint64_t stamp = ++twinStamp_;
        int kept = 0;
        for (int i = 0; i < listed; ++i) {
            const int targetAtom = targetOf(targetKeys_[i]), first = targetTwinFirst_[targetAtom];
            if (first >= 0) {
                if (twinSeen_[first] == stamp) continue;
                twinSeen_[first] = stamp;
            }
            level.targets[kept++] = targetAtom;
        }
        level.targetCount = kept;
    }

    /** Whether the bonds from `queryAtom` to its mapped neighbours match those from `targetAtom` to their images. */
    bool bondsFit(int queryAtom, int targetAtom) {
        clock_.tick(1 + int64_t(query_.degree(queryAtom)) * (1 + target_.degree(targetAtom)));
        for (int queryPlace = query_.start[queryAtom]; queryPlace < query_.start[queryAtom + 1]; ++queryPlace) {
            const int neighbourImage = queryToTarget_[query_.neighbour[queryPlace]];
            if (neighbourImage < 0) continue;
            const int targetBond = target_.bondBetween(targetAtom, neighbourImage);
            if (targetBond < 0 || !graph_.bondsMatch(query_.neighbourBond[queryPlace], targetBond)) return false;
        }
        return true;
    }

    /** The sort key of a target of a query atom: how unlike the two are, then the rank in reverse. */
    int64_t targetKey(int queryAtom, int targetAtom) const {
        return (int64_t(RADII - likeness(queryAtom, targetAtom)) << 32) | int64_t(targetCount_ - 1 - targetRank_[targetAtom]);
    }

    /** The target of a key from targetKey(). */
    int targetOf(int64_t key) const { return targetOrder_[targetCount_ - 1 - static_cast<int>(key & 0xffffffff)]; }

    /**
     * Leaves `queryAtom` of `block` unmapped below the level that holds the block: it moves to the
     * end of the live query atoms of the block, which then leave it out.
     */
    void leaveUnmapped(int block, int queryAtom) {
        blocks_.exactPairs[block] = blocks_.complete[block] = -1;
        const int last = blocks_.queryStart[block] + --blocks_.queryCount[block];
        swapQuery(queryPlace_[queryAtom], last);
    }

    /**
     * Gives the new level the blocks of its parent split by the bonds to the new pair: the query
     * atoms of a block bonded to `queryAtom` and the target atoms bonded to `targetAtom` by bonds
     * of one class form a new block, open as it is bonded to the open fragment, and the atoms
     * bonded to neither stay together. The pair leaves its blocks, and a block with nothing left
     * on one side is dropped, as none of its atoms can then be mapped. A new fragment closes the one
     * before, so the open blocks of the parent, whose query atoms were all tried, are dropped and
     * no other block stays open. Only the blocks of the pair and of its neighbours change, so the
     * others are copied as they are; the blocks of a level are in order of their runs on both
     * sides, which finds the block of an atom by bisection.
     */
    void refine(const Level& parent, int queryAtom, int targetAtom, bool root) {
        const uint64_t stamp = ++linkStamp_;
        const int queryLinked = gatherLinks(parent, query_, queryAtom, true, stamp);
        const int targetLinked = gatherLinks(parent, target_, targetAtom, false, stamp);
        const int queryPairBlock = blockOf(parent, queryPlace_[queryAtom], true);
        const int targetPairBlock = blockOf(parent, targetPlace_[targetAtom], false);
        clock_.tick(1 + parent.blocksTo - parent.blocksFrom);
        for (int block = parent.blocksFrom, queryNext = 0, targetNext = 0; block < parent.blocksTo; ++block) {
            const int queryStart = blocks_.queryStart[block], targetStart = blocks_.targetStart[block];
            int queryAtomCount = blocks_.queryCount[block], targetAtomCount = blocks_.targetCount[block];
            const int links = blocks_.links[block], atomClass = blocks_.atomClass[block];
            const bool open = blocks_.open[block] != 0;
            const int queryFirstLinked = queryNext, targetFirstLinked = targetNext;
            while (queryNext < queryLinked && linkedQueryBlock_[queryNext] == block) ++queryNext;
            while (targetNext < targetLinked && linkedTargetBlock_[targetNext] == block) ++targetNext;
            if (queryAtomCount == 0 || targetAtomCount == 0 || (root && open)) continue;
            if (block != queryPairBlock && block != targetPairBlock && queryNext == queryFirstLinked
                && targetNext == targetFirstLinked) {
                blocks_.add(queryStart, queryAtomCount, targetStart, targetAtomCount, links, open && !root, atomClass);
                continue;
            }
            if (block == queryPairBlock) swapQuery(queryPlace_[queryAtom], queryStart + --queryAtomCount);
            if (block == targetPairBlock) swapTarget(targetPlace_[targetAtom], targetStart + --targetAtomCount);
            if (queryAtomCount == 0 || targetAtomCount == 0) continue;
            // The linked atoms move to the end of the live runs, in order of bond class.
            const int queryEnd = queryStart + queryAtomCount, targetEnd = targetStart + targetAtomCount;
            int queryUnlinked = queryEnd, targetUnlinked = targetEnd;
            for (int i = queryFirstLinked; i < queryNext; ++i) swapQuery(queryPlace_[linkedQuery_[i]], --queryUnlinked);
            for (int i = targetFirstLinked; i < targetNext; ++i) swapTarget(targetPlace_[linkedTarget_[i]], --targetUnlinked);
            sortByLink(queryAtoms_, queryPlace_, queryUnlinked, queryEnd, queryLinkClass_);
            sortByLink(targetAtoms_, targetPlace_, targetUnlinked, targetEnd, targetLinkClass_);
            if (queryUnlinked > queryStart && targetUnlinked > targetStart)
                blocks_.add(queryStart, queryUnlinked - queryStart, targetStart, targetUnlinked - targetStart, links,
                            open && !root, atomClass);
            // The linked atoms of both sides are in order of bond class; equal classes form a block.
            int queryFrom = queryUnlinked, targetFrom = targetUnlinked;
            while (queryFrom < queryEnd && targetFrom < targetEnd) {
                const int queryClass = queryLinkClass_[queryAtoms_[queryFrom]];
                const int targetClass = targetLinkClass_[targetAtoms_[targetFrom]];
                int queryTo = queryFrom, targetTo = targetFrom;
                if (queryClass <= targetClass)
                    while (queryTo < queryEnd && queryLinkClass_[queryAtoms_[queryTo]] == queryClass) ++queryTo;
                if (targetClass <= queryClass)
                    while (targetTo < targetEnd && targetLinkClass_[targetAtoms_[targetTo]] == targetClass) ++targetTo;
                if (queryClass == targetClass && queryClass >= 0)
                    blocks_.add(queryFrom, queryTo - queryFrom, targetFrom, targetTo - targetFrom, links + 1, true, atomClass);
                queryFrom = queryTo;
                targetFrom = targetTo;
            }
        }
    }

    /**
     * Lists the live atoms of the blocks of `parent` that are bonded to `atom` of `side`, with
     * their blocks, in order of block, and marks each with `stamp` and the class of its bond to
     * the atom. Gives their number.
     */
    int gatherLinks(const Level& parent, const Side& side, int atom, bool onQuery, uint64_t stamp) {
        std::vector<int>& linked = onQuery ? linkedQuery_ : linkedTarget_;
        std::vector<int>& linkedBlock = onQuery ? linkedQueryBlock_ : linkedTargetBlock_;
        const std::vector<int>& place = onQuery ? queryPlace_ : targetPlace_;
        std::vector<uint64_t>& linkStampOf = onQuery ? queryLinkStamp_ : targetLinkStamp_;
        std::vector<int>& linkClass = onQuery ? queryLinkClass_ : targetLinkClass_;
        int count = 0;
        clock_.tick(1 + side.degree(atom));
        for (int neighbourPlace = side.start[atom]; neighbourPlace < side.start[atom + 1]; ++neighbourPlace) {
            const int other = side.neighbour[neighbourPlace], block = blockOf(parent, place[other], onQuery);
            if (block < 0) continue;
            linkStampOf[other] = stamp;
            linkClass[other] = side.bondClass[side.neighbourBond[neighbourPlace]];
            // An insertion by block, which keeps the order of the neighbours within a block.
            int i = count++;
            while (i > 0 && linkedBlock[i - 1] > block) {
                linked[i] = linked[i - 1];
                linkedBlock[i] = linkedBlock[i - 1];
                --i;
            }
            linked[i] = other;
            linkedBlock[i] = block;
        }
        return count;
    }

    /**
     * The block of `parent` whose live run of query atoms, or else of target atoms, holds
     * `place`, or -1 if none does. The runs are in order, so this is a bisection.
     */
    int blockOf(const Level& parent, int place, bool onQuery) {
        const std::vector<int>& start = onQuery ? blocks_.queryStart : blocks_.targetStart;
        const std::vector<int>& count = onQuery ? blocks_.queryCount : blocks_.targetCount;
        int low = parent.blocksFrom, high = parent.blocksTo - 1, found = -1;
        while (low <= high) {
            const int middle = static_cast<int>((unsigned(low) + unsigned(high)) >> 1);
            if (start[middle] <= place) {
                found = middle;
                low = middle + 1;
            } else {
                high = middle - 1;
            }
        }
        clock_.tick();
        return found >= 0 && place < start[found] + count[found] ? found : -1;
    }

    /**
     * Sorts the atoms from place `from` up to `to` of `atoms` by the class of their bond to the
     * new pair, keeping their order otherwise. There are few, so an insertion sort will do.
     */
    void sortByLink(std::vector<int>& atoms, std::vector<int>& place, int from, int to, const std::vector<int>& linkClass) {
        clock_.tick(int64_t(to - from) * (to - from));
        for (int i = from + 1; i < to; ++i) {
            const int atom = atoms[i], key = linkClass[atom];
            int j = i;
            while (j > from && linkClass[atoms[j - 1]] > key) {
                atoms[j] = atoms[j - 1];
                place[atoms[j]] = j;
                --j;
            }
            atoms[j] = atom;
            place[atom] = j;
        }
    }

    void swapQuery(int first, int second) {
        std::swap(queryAtoms_[first], queryAtoms_[second]);
        queryPlace_[queryAtoms_[first]] = first;
        queryPlace_[queryAtoms_[second]] = second;
    }

    void swapTarget(int first, int second) {
        std::swap(targetAtoms_[first], targetAtoms_[second]);
        targetPlace_[targetAtoms_[first]] = first;
        targetPlace_[targetAtoms_[second]] = second;
    }

    /**
     * Keeps the mapping if it beats the best score, and otherwise as the fallback if it beats
     * that. The open fragment must have at least minAtoms_ atoms; the closed ones had them.
     */
    void record() {
        if (depth_ - levels_[depth_].fragmentStart + 1 < minAtoms_) return;
        const int64_t current = score(depth_, commonBonds_);
        if (current > bestScore_) {
            bestScore_ = current;
            best_ = queryToTarget_;
            hasBest_ = true;
        } else if (current > fallbackScore_) {
            fallback_ = queryToTarget_;
            hasFallback_ = true;
            fallbackScore_ = current;
        }
    }

    /** Whether the pass can end, as it holds a mapping with a score that no mapping beats. */
    bool done() const { return hasBest_ && bestScore_ >= upperBound_; }

    /**
     * The bound of the state of `level`. Each block pairs off at most as many atoms as it has on
     * its smaller side, see blockPairs(). Each atom it adds brings a bond to each mapped atom its
     * block is linked to, and the bonds between the added atoms are at most the bonds with no end
     * mapped on either side, or more closely those of liveBonds(), which cost more and are counted
     * only when the bonds may decide.
     */
    int64_t bound(const Level& level) {
        int newAtoms = 0;
        int64_t newLinks = 0;
        for (int block = level.blocksFrom; block < level.blocksTo; ++block) {
            if (blocks_.queryCount[block] == 0 || blocks_.targetCount[block] == 0) continue;
            const int pairs = blockPairs(block);
            newAtoms += pairs;
            newLinks += int64_t(pairs) * blocks_.links[block];
        }
        clock_.tick(1 + level.blocksTo - level.blocksFrom);
        const int atoms = depth_ + newAtoms;
        const int64_t loose = capped(atoms, static_cast<int>(std::min<int64_t>(INT_MAX,
            commonBonds_ + newLinks + std::min(level.freeQueryBonds, level.freeTargetBonds))));
        // Bonds decide only between mappings of as many atoms as the best one.
        if (!canImprove(loose) || atomsOf(loose) > atomsOf(bestScore_)) return loose;
        return capped(atoms, static_cast<int>(std::min<int64_t>(INT_MAX, commonBonds_ + newLinks + liveBonds(level))));
    }

    /**
     * The most bonds between atoms the state of `level` may still map. A bond between two such
     * atoms has its image between atoms of the same two blocks, by a bond of the same class, so
     * for each two live blocks and bond class the fewer of these bonds on either side count.
     */
    int liveBonds(const Level& level) {
        const uint64_t stamp = ++liveStamp_;
        int live = 0;
        for (int block = level.blocksFrom; block < level.blocksTo; ++block) {
            const int queryAtomCount = blocks_.queryCount[block], targetAtomCount = blocks_.targetCount[block];
            if (queryAtomCount == 0 || targetAtomCount == 0) continue;
            clock_.tick(int64_t(queryAtomCount) + targetAtomCount);
            for (int place = blocks_.queryStart[block]; place < blocks_.queryStart[block] + queryAtomCount; ++place) {
                queryLive_[queryAtoms_[place]] = stamp;
                queryBlock_[queryAtoms_[place]] = live;
            }
            for (int place = blocks_.targetStart[block]; place < blocks_.targetStart[block] + targetAtomCount; ++place) {
                targetLive_[targetAtoms_[place]] = stamp;
                targetBlock_[targetAtoms_[place]] = live;
            }
            ++live;
        }
        const int classes = std::max(1, graph_.bondClasses);
        if (static_cast<int64_t>(queryPairCount_.size()) < int64_t(live) * classes) {
            queryPairCount_.assign(size_t(2) * live * classes, 0);
            targetPairCount_.assign(size_t(2) * live * classes, 0);
            touchedPairs_.assign(size_t(2) * live * classes, 0);
        }
        int total = 0;
        for (int block = level.blocksFrom, index = 0; block < level.blocksTo; ++block) {
            const int queryAtomCount = blocks_.queryCount[block], targetAtomCount = blocks_.targetCount[block];
            if (queryAtomCount == 0 || targetAtomCount == 0) continue;
            int touched = countBonds(query_, queryAtoms_, blocks_.queryStart[block], queryAtomCount, queryLive_, queryBlock_,
                                     stamp, index, classes, queryPairCount_, targetPairCount_, 0);
            touched = countBonds(target_, targetAtoms_, blocks_.targetStart[block], targetAtomCount, targetLive_, targetBlock_,
                                 stamp, index, classes, targetPairCount_, queryPairCount_, touched);
            for (int i = 0; i < touched; ++i) {
                const int key = touchedPairs_[i];
                total += std::min(queryPairCount_[key], targetPairCount_[key]);
                queryPairCount_[key] = targetPairCount_[key] = 0;
            }
            ++index;
        }
        return total;
    }

    /**
     * Counts in `counts`, by the live block of the other end and the bond class, the bonds of
     * `side` from the `atomCount` atoms from `start` of `atoms`, all of live block `index`, to
     * live atoms of this or a later block, each bond once. A key counted for the first time on
     * both sides is added to touchedPairs_ after the first `touched`, and the new number of keys
     * there is returned.
     */
    int countBonds(const Side& side, const std::vector<int>& atoms, int start, int atomCount,
                   const std::vector<uint64_t>& live, const std::vector<int>& blockOfAtom, uint64_t stamp, int index,
                   int classes, std::vector<int>& counts, const std::vector<int>& otherCounts, int touched) {
        int64_t work = 0;
        for (int place = start; place < start + atomCount; ++place) {
            const int atom = atoms[place];
            work += 1 + side.degree(atom);
            for (int other = side.start[atom]; other < side.start[atom + 1]; ++other) {
                const int neighbour = side.neighbour[other], bondClass = side.bondClass[side.neighbourBond[other]];
                if (live[neighbour] != stamp || bondClass < 0 || blockOfAtom[neighbour] < index
                    || (blockOfAtom[neighbour] == index && neighbour < atom)) continue;
                const int key = blockOfAtom[neighbour] * classes + bondClass;
                if (counts[key]++ == 0 && otherCounts[key] == 0) touchedPairs_[touched++] = key;
            }
        }
        clock_.tick(work);
        return touched;
    }

    /**
     * The most atoms `block` pairs off. A block linked to mapped atoms, where bonds of one class
     * need not match, gives its exact pairs, see exactPairs(). Otherwise it is the fewer of its
     * live atoms on either side or, where its class holds kinds that do not all match, the
     * b-matching of their kinds.
     */
    int blockPairs(int block) {
        const int queryStart = blocks_.queryStart[block], queryAtomCount = blocks_.queryCount[block];
        const int targetStart = blocks_.targetStart[block], targetAtomCount = blocks_.targetCount[block];
        if (!bondClassesMatch_ && blocks_.links[block] > 0 && targetAtomCount <= 64 && queryAtomCount <= 64) {
            if (blocks_.exactPairs[block] < 0) blocks_.exactPairs[block] = exactPairs(block);
            return blocks_.exactPairs[block];
        }
        if (graph_.atomGroups.complete[blocks_.atomClass[block]]) return std::min(queryAtomCount, targetAtomCount);
        clock_.tick(int64_t(queryAtomCount) + targetAtomCount);
        for (int place = queryStart; place < queryStart + queryAtomCount; ++place) ++queryKindCount_[query_.atomKind[queryAtoms_[place]]];
        for (int place = targetStart; place < targetStart + targetAtomCount; ++place) ++targetKindCount_[target_.atomKind[targetAtoms_[place]]];
        const int pairs = graph_.atomGroups.bound(blocks_.atomClass[block], queryKindCount_, targetKindCount_);
        for (int place = queryStart; place < queryStart + queryAtomCount; ++place) queryKindCount_[query_.atomKind[queryAtoms_[place]]] = 0;
        for (int place = targetStart; place < targetStart + targetAtomCount; ++place) targetKindCount_[target_.atomKind[targetAtoms_[place]]] = 0;
        return pairs;
    }

    /**
     * The most pairs of `block`, a block linked to mapped atoms with at most 64 atoms on each
     * side, that match exactly: the atoms and the bonds to the mapped atoms. Bonds of one class
     * need not match, so the pairs are a maximum bipartite matching. Such blocks are small, as
     * their atoms are bonded to a mapped atom.
     */
    int exactPairs(int block) {
        const int queryStart = blocks_.queryStart[block], queryAtomCount = blocks_.queryCount[block];
        const int targetStart = blocks_.targetStart[block], targetAtomCount = blocks_.targetCount[block];
        for (int row = 0; row < queryAtomCount; ++row) {
            const int queryAtom = queryAtoms_[queryStart + row];
            uint64_t columns = 0;
            for (int column = 0; column < targetAtomCount; ++column) {
                const int targetAtom = targetAtoms_[targetStart + column];
                if (has(compatible_, queryAtom, targetAtom) && bondsFit(queryAtom, targetAtom)) columns |= uint64_t(1) << column;
            }
            matchRows_[row] = columns;
        }
        std::fill(matchOwner_.begin(), matchOwner_.begin() + targetAtomCount, -1);
        int pairs = 0;
        for (int row = 0; row < queryAtomCount; ++row) {
            matchSeen_ = 0;
            if (augment(row)) ++pairs;
        }
        return pairs;
    }

    /**
     * Whether row `row` of matchRows_ gets a column, by an augmenting path through the columns not
     * yet seen. A path has at most one step for each row, at most 64.
     */
    bool augment(int row) {
        uint64_t columns = matchRows_[row] & ~matchSeen_;
        while (columns != 0) {
            clock_.tick();
            const int column = ctz64(columns);
            columns &= columns - 1;
            matchSeen_ |= uint64_t(1) << column;
            if (matchOwner_[column] < 0 || augment(matchOwner_[column])) {
                matchOwner_[column] = row;
                return true;
            }
        }
        return false;
    }

    /**
     * Whether every atom pair of `block` matches exactly, the atoms and the bonds to the mapped
     * atoms, which holds when bonds of one class always match and otherwise is checked for blocks
     * of at most COMPLETE_PAIRS pairs.
     */
    bool blockComplete(int block) {
        if (!graph_.atomGroups.complete[blocks_.atomClass[block]]) return false;
        if (bondClassesMatch_ || blocks_.links[block] == 0) return true;
        if (blocks_.complete[block] >= 0) return blocks_.complete[block] == 1;
        blocks_.complete[block] = 0;
        const int queryStart = blocks_.queryStart[block], queryAtomCount = blocks_.queryCount[block];
        const int targetStart = blocks_.targetStart[block], targetAtomCount = blocks_.targetCount[block];
        if (int64_t(queryAtomCount) * targetAtomCount > COMPLETE_PAIRS) return false;
        for (int queryPlace = queryStart; queryPlace < queryStart + queryAtomCount; ++queryPlace)
            for (int targetPlace = targetStart; targetPlace < targetStart + targetAtomCount; ++targetPlace)
                if (!bondsFit(queryAtoms_[queryPlace], targetAtoms_[targetPlace])) return false;
        blocks_.complete[block] = 1;
        return true;
    }

    /** The score of a mapping: its atoms, and the bonds between them. */
    int64_t scoreOf(const std::vector<int>& mapping) const {
        int atomCount = 0, bondCount = 0;
        for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
            if (mapping[queryAtom] < 0) continue;
            ++atomCount;
            for (int place = query_.start[queryAtom]; place < query_.start[queryAtom + 1]; ++place)
                if (query_.neighbour[place] < queryAtom && mapping[query_.neighbour[place]] >= 0) ++bondCount;
        }
        return score(atomCount, bondCount);
    }

    /** A bound capped by the upper bound, which caps the atoms and the bonds of every mapping apart. */
    int64_t capped(int atomCount, int bondCount) const {
        return score(std::min(atomCount, atomsOf(upperBound_)), std::min(bondCount, bondsOf(upperBound_)));
    }

    /** Whether a state with this bound can hold a better mapping. */
    bool canImprove(int64_t bound) const { return bound > bestScore_; }

    /**
     * The twins of one molecule, each linked to the next in a ring, -1 for none. Twins are
     * terminal atoms with the same neighbour that match the same atoms by the same bonds.
     * Swapping two takes every mapping onto one of the same score.
     */
    std::vector<int> twins(const Side& side, bool onQuery) {
        std::vector<int> next(side.atomCount, -1);
        for (int parent = 0; parent < side.atomCount; ++parent)
            for (int place = side.start[parent]; place < side.start[parent + 1]; ++place) {
                const int first = side.neighbour[place];
                int last = first;
                if (next[first] >= 0 || side.degree(first) != 1) continue;
                for (int laterPlace = place + 1; laterPlace < side.start[parent + 1]; ++laterPlace) {
                    const int other = side.neighbour[laterPlace];
                    if (next[other] < 0 && side.degree(other) == 1
                        && alike(onQuery, first, side.neighbourBond[place], other, side.neighbourBond[laterPlace])) {
                        next[last] = other;
                        last = other;
                    }
                }
                if (last != first) next[last] = first;
            }
        return next;
    }

    /** The first atom of the ring of twins of each atom, its lowest index, or -1 for an atom with no twin. */
    static std::vector<int> ringFirsts(const std::vector<int>& twin) {
        std::vector<int> first(twin.size(), -1);
        for (int atom = 0; atom < static_cast<int>(twin.size()); ++atom) {
            if (twin[atom] < 0 || first[atom] >= 0) continue;
            // The ring is met from its lowest atom first.
            for (int member = atom; first[member] < 0; member = twin[member]) first[member] = atom;
        }
        return first;
    }

    /** Whether two atoms of one molecule, with their bonds, match alike. */
    bool alike(bool onQuery, int atom, int bond, int other, int otherBond) {
        clock_.tick(int64_t(words_) + queryCount_ + query_.bondCount + target_.bondCount);
        if (onQuery) {
            for (int word = 0; word < words_; ++word)
                if (compatible_[size_t(atom) * words_ + word] != compatible_[size_t(other) * words_ + word]) return false;
            for (int targetBond = 0; targetBond < target_.bondCount; ++targetBond)
                if (graph_.bondsMatch(bond, targetBond) != graph_.bondsMatch(otherBond, targetBond)) return false;
        } else {
            for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom)
                if (has(compatible_, queryAtom, atom) != has(compatible_, queryAtom, other)) return false;
            for (int queryBond = 0; queryBond < query_.bondCount; ++queryBond)
                if (graph_.bondsMatch(queryBond, bond) != graph_.bondsMatch(queryBond, otherBond)) return false;
        }
        return true;
    }
};

/** What a common-atom search found. */
struct CommonAtomOutcome {
    /** Whether the search could take the graphs; when not, nothing else is set. */
    bool supported = false;
    /** Whether the search ran to its end, so that the mapping is a maximum one. */
    bool finished = false;
    /** The query-to-target atom mapping, -1 where unmapped, or empty when nothing was found. */
    std::vector<int> mapping;
    /** The bonds between mapped query atoms, each as {query atom, query atom, target atom, target atom}, the lower query atom first, in the order of the query bonds. */
    std::vector<std::array<int, 4>> bonds;
    /** The pairs mapped. */
    uint64_t nodes = 0;
};

/**
 * Searches `query` in `target` for the single best MCIS in fragments of at least
 * `minFragmentAtoms` atoms, one or two, under the given atom and bond predicates, with the
 * target atoms in `excludedTargetAtoms` left out. `stopped()` is asked now and then whether the
 * time is up or the search was cancelled; the best mapping found by then is the result. Each query
 * atom is either mapped or left unmapped, so the molecule with fewer atoms is searched as the
 * query, unless target atoms are left out; the predicates are asked in their own orientation. The
 * outcome is not supported when the graphs are not well formed or too large, see
 * CommonBondGraphs::supports().
 */
template <class AtomMatch, class BondMatch, class BondKey, class Stopped>
inline CommonAtomOutcome searchCommonAtoms(const MolGraph& query, const MolGraph& target, AtomMatch atomMatches,
                                           BondMatch bondMatches, BondKey bondKey, const std::vector<int>& excludedTargetAtoms,
                                           int minFragmentAtoms, Stopped stopped) {
    CommonAtomOutcome outcome;
    if (!CommonBondGraphs::supports(query, target)) return outcome;
    outcome.supported = true;
    const bool reversed = excludedTargetAtoms.empty() && query.n > target.n;
    Clock clock{std::function<bool()>(std::move(stopped))};
    std::unique_ptr<CommonBondGraphs> graph;
    std::unique_ptr<CommonAtomSearch> search;
    std::vector<int> found;
    try {
        clock.poll();
        if (reversed) {
            auto reversedAtoms = [&](int first, int second) { return atomMatches(second, first); };
            auto reversedBonds = [&](int firstFrom, int firstTo, int secondFrom, int secondTo) {
                return bondMatches(secondFrom, secondTo, firstFrom, firstTo);
            };
            graph = std::make_unique<CommonBondGraphs>(target, query, reversedAtoms, reversedBonds, bondKey, excludedTargetAtoms, clock);
        } else {
            graph = std::make_unique<CommonBondGraphs>(query, target, atomMatches, bondMatches, bondKey, excludedTargetAtoms, clock);
        }
        search = std::make_unique<CommonAtomSearch>(*graph, clock, minFragmentAtoms);
        found = search->run();
        outcome.finished = true;
    } catch (const Stop&) {
        // The time ran out or the search was cancelled. What the search found by then is the result.
        if (search) found = search->found();
    }
    if (!found.empty()) {
        const Side& side = reversed ? graph->target : graph->query;
        std::vector<int> mapping = found;
        if (reversed) {
            mapping.assign(size_t(side.atomCount), -1);
            for (int atom = 0; atom < static_cast<int>(found.size()); ++atom)
                if (found[atom] >= 0) mapping[found[atom]] = atom;
        }
        for (int bond = 0; bond < side.bondCount; ++bond) {
            const int first = mapping[side.bondFirst[bond]], second = mapping[side.bondSecond[bond]];
            if (first >= 0 && second >= 0) outcome.bonds.push_back({side.bondFirst[bond], side.bondSecond[bond], first, second});
        }
        outcome.mapping = std::move(mapping);
    }
    if (search) outcome.nodes = static_cast<uint64_t>(search->nodes());
    return outcome;
}

} // namespace common_atom_detail
} // namespace smsd

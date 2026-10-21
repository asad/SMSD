/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Common-bond (MCES) search: embeddings of a whole connected query, which answer
 * the search when they exist, and the entry of CommonBondSearch::run().
 *
 * Part of common_bond_search.hpp, which describes the search as a whole.
 */
#pragma once

#include "smsd/detail/common_bond_mapping_search.hpp"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace smsd {
namespace common_bond_detail {

/**
 * The embeddings of a connected query: maps of every query atom that take each query bond onto
 * a compatible target bond, with other target bonds allowed. They are found by backtracking
 * with forward checking. Placing an atom narrows the candidates of each unplaced neighbour to
 * the free neighbours of its target, and these neighbours must be able to have distinct
 * targets. The next atom placed is the one with the fewest free candidates. There is no
 * recursion, as the depth is the number of query atoms. The candidate lists live on a stack
 * that each placement trims back when it is undone.
 */
class CommonBondSearch::Embedding {
    /**
     * A placed atom, one for each depth: the place it had on the frontier, the next of its
     * candidates to try, and the candidates and counts that placing it changed as they were
     * before, so that taking it off puts them back.
     */
    struct Placement {
        int atom = 0, frontierPlace = 0, nextCandidate = 0;
        size_t stackMark = 0;
        std::vector<int> changedAtoms, oldStarts, oldLengths, oldFreeCounts;
        void enter(int placedAtom, int place) {
            atom = placedAtom;
            frontierPlace = place;
            nextCandidate = 0;
        }
        void remember(int changedAtom, int start, int length, int freeCount) {
            changedAtoms.push_back(changedAtom);
            oldStarts.push_back(start);
            oldLengths.push_back(length);
            oldFreeCounts.push_back(freeCount);
        }
    };

    CommonBondSearch& s_;
    /** The target atom each query atom is placed on, and the query atom placed on each target atom, or -1. */
    std::vector<int> placedOn_, placedAtom_;
    /** The candidate lists, and where the list of each query atom starts, -1 until a neighbour is placed, and its length. */
    std::vector<int> stack_, candidateStart_, candidateLength_;
    /** The number of candidates of each unplaced query atom that are free. */
    std::vector<int> freeCount_;
    /** The unplaced query atoms that have candidates, the first frontierSize_, in any order. */
    std::vector<int> frontier_;
    int frontierSize_ = 0;
    /** Room for the candidates of a neighbour while they are found, and for the rows of the distinct-target check. */
    std::vector<int> scratch_;
    std::vector<uint64_t> rows_;

public:
    /** Whether only embeddings with their twins in order are tried, and the first that counts ends the search. */
    const bool prunesTwins;

    /** Starts from the query atom with the fewest targets it may have, then most bonds, then lowest index. */
    Embedding(CommonBondSearch& search, const std::vector<uint64_t>& candidateRows, bool prunesTwins)
        : s_(search), placedOn_(search.queryCount_, -1), placedAtom_(search.targetCount_, -1),
          candidateStart_(search.queryCount_, -1), candidateLength_(search.queryCount_, 0),
          freeCount_(search.queryCount_, 0), frontier_(search.queryCount_, 0),
          scratch_(size_t(std::max(1, search.target_.maxDegree())), 0), prunesTwins(prunesTwins) {
        if (prunesTwins) s_.findTwins();
        const int words = s_.words_;
        int root = 0, fewest = INT_MAX;
        for (int queryAtom = 0; queryAtom < s_.queryCount_; ++queryAtom) {
            int count = 0;
            for (int i = 0; i < words; ++i) count += popcount64(candidateRows[size_t(queryAtom) * words + i]);
            if (count < fewest || (count == fewest && s_.query_.degree(queryAtom) > s_.query_.degree(root))) {
                root = queryAtom;
                fewest = count;
            }
        }
        const uint64_t* row = candidateRows.data() + size_t(root) * words;
        candidateStart_[root] = 0;
        for (int targetAtom = nextBit(row, words, 0); targetAtom >= 0; targetAtom = nextBit(row, words, targetAtom + 1))
            stack_.push_back(targetAtom);
        candidateLength_[root] = fewest;
        freeCount_[root] = fewest;
        frontier_[frontierSize_++] = root;
    }

    /**
     * The embeddings, at most `maxMappings`, held as the mappings as they are found. Under
     * conditions only an embedding that meets them with every query bond counts, as only it has
     * the score of an embedding; if none does, no mapping meets them with every query bond. When
     * it prunes twins, the first embedding that counts ends it. Called once only.
     */
    std::vector<std::vector<int>> all(size_t maxMappings) {
        s_.mappings_.clear();
        const int64_t whole = score(s_.queryCount_, s_.query_.bondCount);
        std::vector<Placement> placements(size_t(s_.queryCount_));
        placements[0].enter(frontier_[0], 0);
        for (int placed = 1; placed > 0;) {
            s_.clock_.tick();
            Placement& top = placements[placed - 1];
            if (placedOn_[top.atom] >= 0) unplace(top);
            const int targetAtom = nextTarget(top);
            if (targetAtom < 0) {
                --placed;
            } else if (place(top, targetAtom)) {
                const int next = choose();
                if (next >= 0) {
                    placements[placed++].enter(frontier_[next], next);
                } else if (placed == s_.queryCount_) {
                    if (!s_.conditions_) s_.mappings_.push_back(placedOn_);
                    else s_.conditions_->select(placedOn_, whole, maxMappings - s_.mappings_.size(), s_.mappings_);
                    if (s_.mappings_.size() == maxMappings || (prunesTwins && !s_.mappings_.empty())) return s_.mappings_;
                }
            }
        }
        return s_.mappings_;
    }

private:
    /** The next free candidate of the atom of a placement, or -1 once it has tried them all. */
    int nextTarget(Placement& placement) {
        const int start = candidateStart_[placement.atom], length = candidateLength_[placement.atom];
        while (placement.nextCandidate < length) {
            const int targetAtom = stack_[size_t(start) + placement.nextCandidate++];
            if (placedAtom_[targetAtom] < 0 && !(prunesTwins && twinsOutOfOrder(placement.atom, targetAtom)))
                return targetAtom;
        }
        return -1;
    }

    /**
     * Whether placing `queryAtom` on `targetAtom` puts twins out of order: a placed twin of the
     * query atom lies on the other side of the target atom, by index, than of the query atom, or
     * a twin of the target atom in use has its query atom on the other side of the query atom.
     * Exchanging twins of either molecule takes an embedding onto one that meets the conditions
     * alike, and the exchanges take every embedding onto one with all its twins in order: the
     * least of them, each read as the targets of the query atoms in turn.
     */
    bool twinsOutOfOrder(int queryAtom, int targetAtom) {
        for (int twin = s_.queryTwin_[queryAtom]; twin >= 0 && twin != queryAtom; twin = s_.queryTwin_[twin]) {
            s_.clock_.tick();
            if (placedOn_[twin] >= 0 && (twin < queryAtom) != (placedOn_[twin] < targetAtom)) return true;
        }
        for (int twin = s_.targetTwin_[targetAtom]; twin >= 0 && twin != targetAtom; twin = s_.targetTwin_[twin]) {
            s_.clock_.tick();
            if (placedAtom_[twin] >= 0 && (twin < targetAtom) != (placedAtom_[twin] < queryAtom)) return true;
        }
        return false;
    }

    /** Places the atom of a placement on a target atom and records the changes in the placement, to be undone on failure. */
    bool place(Placement& placement, int targetAtom) {
        const Side& query = s_.query_;
        const Side& target = s_.target_;
        s_.countNode();
        placement.stackMark = stack_.size();
        const int queryAtom = placement.atom;
        frontier_[placement.frontierPlace] = frontier_[--frontierSize_];
        placedOn_[queryAtom] = targetAtom;
        placedAtom_[targetAtom] = queryAtom;
        take(placement, targetAtom);
        const int from = target.start[targetAtom], degree = target.degree(targetAtom);
        int open = 0, narrowed = 0;
        s_.clock_.tick(1 + query.degree(queryAtom));
        for (int place = query.start[queryAtom]; place < query.start[queryAtom + 1]; ++place)
            if (placedOn_[query.neighbour[place]] < 0) ++open;
        // The candidates of each neighbour as places among the target's neighbours, for the check at the end.
        const bool checkRows = open > 1 && degree <= 64;
        if (checkRows) rows_.assign(size_t(open), 0);
        for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1]; ++queryPlace) {
            const int queryNeighbour = query.neighbour[queryPlace];
            int count = 0;
            if (placedOn_[queryNeighbour] >= 0) continue;
            // A neighbour with candidates has another placed neighbour, so it closes a ring.
            const bool closes = candidateStart_[queryNeighbour] >= 0;
            s_.clock_.tick(1 + (closes ? int64_t(degree) * (1 + query.degree(queryNeighbour)) : int64_t(degree)));
            uint64_t bits = 0;
            for (int offset = 0; offset < degree; ++offset) {
                const int candidate = target.neighbour[from + offset];
                if (placedAtom_[candidate] < 0
                    && s_.graph_.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[from + offset])
                    && mayMap(queryNeighbour, candidate)
                    && (!closes || bondsFit(queryNeighbour, candidate, queryAtom))) {
                    scratch_[count++] = candidate;
                    if (checkRows) bits |= uint64_t(1) << offset;
                }
            }
            if (count == 0) return false;
            if (!closes) frontier_[frontierSize_++] = queryNeighbour;
            placement.remember(queryNeighbour, candidateStart_[queryNeighbour], candidateLength_[queryNeighbour],
                               freeCount_[queryNeighbour]);
            candidateStart_[queryNeighbour] = static_cast<int>(stack_.size());
            candidateLength_[queryNeighbour] = count;
            stack_.insert(stack_.end(), scratch_.begin(), scratch_.begin() + count);
            freeCount_[queryNeighbour] = count;
            if (checkRows) rows_[narrowed++] = bits;
        }
        return !checkRows || hasDistinctTargets(rows_, narrowed, 1, s_.clock_);
    }

    /** Takes the atom of a placement off its target, and puts back everything that placing it changed. */
    void unplace(Placement& placement) {
        const int queryAtom = placement.atom;
        while (!placement.changedAtoms.empty()) {
            const size_t changed = placement.changedAtoms.size() - 1;
            const int changedAtom = placement.changedAtoms[changed];
            // A neighbour that had no candidates came on the frontier last.
            if (placement.oldStarts[changed] < 0) --frontierSize_;
            candidateStart_[changedAtom] = placement.oldStarts[changed];
            candidateLength_[changedAtom] = placement.oldLengths[changed];
            freeCount_[changedAtom] = placement.oldFreeCounts[changed];
            placement.changedAtoms.pop_back();
            placement.oldStarts.pop_back();
            placement.oldLengths.pop_back();
            placement.oldFreeCounts.pop_back();
        }
        stack_.resize(placement.stackMark);
        placedAtom_[placedOn_[queryAtom]] = -1;
        placedOn_[queryAtom] = -1;
        frontier_[frontierSize_++] = frontier_[placement.frontierPlace];
        frontier_[placement.frontierPlace] = queryAtom;
    }

    /**
     * The frontier place of the atom to place next, with the fewest free candidates, then the
     * most bonds, then the lowest index. It is -1 if the frontier is empty or an atom on it has
     * no free candidate left.
     */
    int choose() {
        const Side& query = s_.query_;
        s_.clock_.tick(frontierSize_);
        int pick = -1;
        for (int frontierIndex = 0; frontierIndex < frontierSize_; ++frontierIndex) {
            const int candidate = frontier_[frontierIndex], picked = pick < 0 ? -1 : frontier_[pick];
            if (freeCount_[candidate] == 0) return -1;
            if (picked < 0 || freeCount_[candidate] < freeCount_[picked]
                || (freeCount_[candidate] == freeCount_[picked]
                    && (query.degree(candidate) > query.degree(picked)
                        || (query.degree(candidate) == query.degree(picked) && candidate < picked))))
                pick = frontierIndex;
        }
        return pick;
    }

    /**
     * Takes a target atom from the free candidates of the unplaced atoms that have it,
     * remembering their counts in the placement. Such an atom may map onto the target atom
     * with compatible bonds to the targets of all its placed neighbours, and is counted from
     * the first of these.
     */
    void take(Placement& placement, int targetAtom) {
        const Side& query = s_.query_;
        const Side& target = s_.target_;
        int64_t work = 1 + target.degree(targetAtom);
        for (int targetPlace = target.start[targetAtom]; targetPlace < target.start[targetAtom + 1]; ++targetPlace) {
            const int placedNeighbour = placedAtom_[target.neighbour[targetPlace]];
            if (placedNeighbour < 0) continue;
            for (int queryPlace = query.start[placedNeighbour]; queryPlace < query.start[placedNeighbour + 1]; ++queryPlace) {
                const int queryNeighbour = query.neighbour[queryPlace];
                ++work;
                if (placedOn_[queryNeighbour] >= 0
                    || !s_.graph_.bondsMatch(query.neighbourBond[queryPlace], target.neighbourBond[targetPlace])
                    || !mayMap(queryNeighbour, targetAtom)) continue;
                work += query.degree(queryNeighbour);
                if (firstPlacedNeighbour(queryNeighbour) == placedNeighbour && bondsFit(queryNeighbour, targetAtom, placedNeighbour)) {
                    placement.remember(queryNeighbour, candidateStart_[queryNeighbour], candidateLength_[queryNeighbour],
                                       freeCount_[queryNeighbour]);
                    --freeCount_[queryNeighbour];
                }
            }
        }
        s_.clock_.tick(work);
    }

    /** Whether a query atom may map onto a target atom: they match, and the target has as many bonds. */
    bool mayMap(int queryAtom, int targetAtom) const {
        return s_.has(s_.allowed_, queryAtom, targetAtom) && s_.target_.degree(targetAtom) >= s_.query_.degree(queryAtom);
    }

    /** Whether a target atom has compatible bonds to the images of the placed neighbours of a query atom other than `skipped`. */
    bool bondsFit(int queryAtom, int targetAtom, int skipped) const {
        const Side& query = s_.query_;
        const Side& target = s_.target_;
        for (int queryPlace = query.start[queryAtom]; queryPlace < query.start[queryAtom + 1]; ++queryPlace) {
            const int queryNeighbour = query.neighbour[queryPlace], neighbourImage = placedOn_[queryNeighbour];
            if (neighbourImage < 0 || queryNeighbour == skipped) continue;
            // Search the shorter list of neighbours, as one may be a hub.
            const int targetBond = target.degree(targetAtom) <= target.degree(neighbourImage)
                ? target.bondBetween(targetAtom, neighbourImage) : target.bondBetween(neighbourImage, targetAtom);
            if (targetBond < 0 || !s_.graph_.bondsMatch(query.neighbourBond[queryPlace], targetBond)) return false;
        }
        return true;
    }

    /** The first placed neighbour of a query atom, or -1 if none is placed. */
    int firstPlacedNeighbour(int queryAtom) const {
        const Side& query = s_.query_;
        for (int i = query.start[queryAtom]; i < query.start[queryAtom + 1]; ++i)
            if (placedOn_[query.neighbour[i]] >= 0) return query.neighbour[i];
        return -1;
    }
};

inline std::unique_ptr<CommonBondSearch::Embedding> CommonBondSearch::embedding(bool prunesTwins) {
    if (queryCount_ > targetCount_ || !queryConnected()) return nullptr;
    // A row for each number of bonds holds the target atoms with at least that many, up to the
    // most a query atom has. The last row is empty.
    const int rows = std::min(query_.maxDegree(), target_.maxDegree()) + 1;
    std::vector<uint64_t> atLeast(size_t(rows + 1) * size_t(words_), 0);
    for (int targetAtom = 0; targetAtom < targetCount_; ++targetAtom)
        for (int fewestBonds = std::min(target_.degree(targetAtom), rows); fewestBonds >= 0; --fewestBonds)
            set(atLeast, fewestBonds, targetAtom);
    std::vector<uint64_t> candidateRows(size_t(queryCount_) * size_t(words_), 0);
    for (int queryAtom = 0; queryAtom < queryCount_; ++queryAtom) {
        clock_.tick(words_);
        const size_t row = size_t(std::min(query_.degree(queryAtom), rows)) * words_;
        for (int i = 0; i < words_; ++i)
            candidateRows[size_t(queryAtom) * words_ + i] = allowed_[size_t(queryAtom) * words_ + i] & atLeast[row + i];
    }
    if (!hasDistinctTargets(candidateRows, queryCount_, words_, clock_)) return nullptr;
    return std::make_unique<Embedding>(*this, candidateRows, prunesTwins);
}

inline std::vector<std::vector<int>> CommonBondSearch::run() {
    clock_.poll();
    if (!anyCompatible() || std::min(query_.bondCount, target_.bondCount) < MIN_BONDS) return {};
    // Embeddings that differ by an exchange of twins meet the conditions alike, so under
    // conditions only one of each kind is tried until one meets them. Then, if more are wanted,
    // every embedding is tried.
    std::unique_ptr<Embedding> whole = embedding(conditions_ && conditions_->keepsTwins());
    if (whole) {
        std::vector<std::vector<int>> embeddings;
        try {
            embeddings = whole->all(wanted_);
            if (whole->prunesTwins && wanted_ > 1 && !embeddings.empty()) embeddings = embedding(false)->all(wanted_);
        } catch (const Stop&) {
            // At the node limit before an embedding is found, the seeds, which take no nodes,
            // still give a mapping.
            if (capReached_ && mappings_.empty()) {
                rank();
                fallbackScore_ = belowEveryMapping();
                seed(bestPossible());
            }
            throw;
        }
        if (!embeddings.empty()) return embeddings;
    }
    rank();
    const int64_t bound = bestPossible();
    fallbackScore_ = belowEveryMapping();
    seed(bound);
    return collectByBonds(bound);
}

} // namespace common_bond_detail
} // namespace smsd

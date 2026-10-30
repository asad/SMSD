/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The ring candidates of perceiveRingSets(), read from one shortest-path tree per root.
 *
 * Part of ring_finder.hpp, which includes this header before its own code; including it
 * first includes ring_finder.hpp, which includes it again at that point.
 */
#ifndef SMSD_DETAIL_RING_CANDIDATES_HPP
#include "smsd/ring_finder.hpp"
#endif
#ifndef SMSD_DETAIL_RING_CANDIDATES_HPP
#define SMSD_DETAIL_RING_CANDIDATES_HPP

#include "smsd/bitops.hpp"
#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {
namespace detail_rf {

/** Internal switches of the ring candidate search; every one is off by default. */
struct RingSearchOptions {
    /** Where the fallback basis is not certified, finish the search and certify its basis instead. */
    bool liftRefusedBasis = false;
    /** Clear an incomplete relevant-cycle flag when the cycle families prove the list complete. */
    bool proveRelevantByFamilies = false;
    /** Tests only: a bracketed pair count has no upper bound, so only replays keep the main path. */
    bool loosenBracketsForTesting = false;
    /** Tests only: every bracketed pair count is replayed, so the decision reads the exact work. */
    bool replayAllForTesting = false;
};

/** What one ring candidate search did, for tests and diagnostics; no result reads it. */
struct RingSearchNotes {
    bool pairwise = false;              // the pairwise search itself ran
    uint64_t pairs = 0;                 // root and bond pairs of the pairwise search
    uint64_t pairsRead = 0;             // pairs read from trees, the last one included
    uint64_t nonTree = 0;               // a bond outside the root's tree
    uint64_t classA = 0;                // a tree bond whose child is reached from the level above it
    uint64_t classB1 = 0;               // ... from the child's own level
    uint64_t classB1Plus = 0;           // ... from the level below the child
    uint64_t classB2 = 0;               // ... through a descendant of the child (bracketed count)
    uint64_t simpleRing = 0;            // the one pair of a system that is a simple ring
    uint64_t baseWorkLow = 0;           // bounds on the pairwise work before the last pair
    uint64_t baseWorkHigh = 0;
    bool baseWorkHighUnbounded = false;
    bool baseWorkExact = false;
    bool lastPairSummed = false;        // stays false: the last pair's work reaches no test
    bool fellBack = false;              // the pairwise search falls back to the bounded pipeline
    uint64_t tier2Replays = 0;          // pair searches replayed to decide the fallback
    uint64_t searchWork = 0;            // work of this search, as RingCandidateSearch counts it
    uint64_t storedKeyWords = 0;        // words of the stored candidate bond sets
    uint64_t candidates = 0;            // distinct candidate cycles stored
    bool baseBasisTruncated = false;    // the basis of the pairwise search is not certified
    bool basisLifted = false;           // ... and the finished search certified its own basis
    bool relevantProvedByFamilies = false;
};

/** The fixed 64-bit hash of a candidate bond set, the same in the Java implementation. */
inline uint64_t ringCandidateKeyHash(const uint64_t* words, size_t count) noexcept {
    uint64_t hash = 0;
    for (size_t i = 0; i < count; ++i)
        hash = (((hash << 5) | (hash >> 59)) ^ words[i]) * UINT64_C(0x9E3779B97F4A7C15);
    hash ^= hash >> 33;
    hash *= UINT64_C(0xFF51AFD7ED558CCD);
    hash ^= hash >> 33;
    return hash;
}

/** Beyond this many words per root, tree-path bond sets are walked rather than kept. */
inline constexpr uint64_t RING_TREE_PATH_SET_WORDS = uint64_t(1) << 22;

/** The most candidate bond-set words that finishing a refused search may keep. */
inline constexpr uint64_t RING_LIFT_KEY_WORDS = uint64_t(1) << 24;

/**
 * The candidates of the pairwise ring search, read from one breadth-first tree per root.
 *
 * The pairwise search pairs each root with each bond of its ring system and searches the
 * system without that bond, reading neighbours in list order. Its paths and the atoms it
 * reads follow from the tree of the root over its system: a bond outside the tree leaves the
 * tree paths in place, and for a tree bond the path to its child enters the child's subtree
 * through the shortest closing bond that leaves the subtree, ties going to the earlier one in
 * the tree's depth-first order, then climbs to the child through the first listed neighbour
 * one level up. Each subtree receives that closing bond through disjoint sets over the tree.
 * The candidates, their order, their atom sequences and the work count that decides the
 * fallback are therefore those of the pairwise search. Where a pair's count is only bracketed,
 * the search of that pair is replayed when the decision needs it.
 *
 * Per root of a system with a atoms, b bonds, cycle rank r, tree depth D, w = ceil(b/64) words
 * per bond set and largest degree Delta: the tree and its numbering take O(a + b), the
 * ordering O(a + r), the covers O((a + r) alpha(a)) and the tree-path bond sets O(a w). Each
 * bond takes O(1) for its count and validity, O(w) for its bond set and O(Delta D) for the
 * climb of a valid tree bond; each stored candidate takes O(len + w), and each stored bond set
 * that shares the 64-bit hash O(w). Selection takes O(w) per row reduced. The worst case per
 * system is O(a b (w + Delta D)), against Theta(a b (a + b)) for the pairwise search.
 *
 * searchWork counts, identically in the Java implementation: per root, the atoms and
 * neighbours of the system once for the tree and once for its numbering, the closing
 * half-bonds plus 2D + 1 for their ordering, one per set-root step and one per merge of the
 * covers, and a w when the tree-path bond sets are kept; per pair, one for its count and
 * validity, and for a valid pair the atoms and neighbours its climb reads, w for its bond set
 * (len + w when it is walked), w for each stored bond set compared because its hash is equal,
 * in the order stored and up to the first equal one, and len + w when it is stored; a simple
 * ring's pair costs its exact pairwise count plus len; selection costs w per candidate it
 * reduces and w per row it adds; a replay costs its pair's search, and the tree it rebuilds for
 * the replay is not counted. Finishing a refused search builds the tree of the root it resumes
 * again, counted as any tree, and tests its bounds before each pair and before the selection.
 */
class RingCandidateSearch {
public:
    /** How a pair's search differs from the tree of its root. */
    enum class PairKind : int { NonTree, ClassA, ClassB1, ClassB1Plus, ClassB2, SimpleRing };

    /** One root and bond pair, read from the tree of the root. */
    struct PairReading {
        int edge = -1;
        PairKind kind = PairKind::NonTree;
        bool valid = false;                 // the two paths share only the root
        uint64_t searchLow = 0;             // atoms and neighbours the pair search reads, or bounds on them
        uint64_t searchHigh = 0;
        uint64_t pathTerm = 0;              // atoms of both paths
        int child = -1, parent = -1;        // a tree bond's atoms
        int coverFrom = -1, coverTo = -1;   // its closing half-bond into the child's subtree
        uint64_t climbWork = 0;             // atoms and neighbours that a valid pair's climb reads
    };

    RingCandidateSearch(const MolGraph& g, const MolGraph::CycleSystems& systems,
                        const std::unordered_map<int64_t, int>& edgeIndex,
                        const std::vector<std::pair<int, int>>& edgeList, uint64_t workLimit,
                        const RingSearchOptions& options, RingSearchNotes& notes)
        : g_(g), ofAtom_(systems.ofAtom), simpleRing_(systems.simpleRing), edgeList_(edgeList),
          workLimit_(workLimit), options_(options), notes_(notes), n_(g.n) {
        if (systems.ofAtom.empty() || edgeList.size() != edgeIndex.size()) return;
        // Each listed neighbour lists the atom back; otherwise the pairwise search runs.
        slotStart_.assign(static_cast<size_t>(n_) + 1, 0);
        for (int u = 0; u < n_; ++u)
            slotStart_[u + 1] = slotStart_[u] + static_cast<int>(g.neighbors[u].size());
        std::vector<int> slotEdge(static_cast<size_t>(slotStart_[n_]), -1);
        std::vector<uint8_t> ends(edgeIndex.size(), 0);
        for (int u = 0; u < n_; ++u) {
            const auto& row = g.neighbors[u];
            for (size_t i = 0; i < row.size(); ++i) {
                const auto it = edgeIndex.find(MolGraph::bondKey(u, row[i]));
                if (it == edgeIndex.end() || ++ends[it->second] > 2) return;
                slotEdge[slotStart_[u] + static_cast<int>(i)] = it->second;
            }
        }
        for (uint8_t count : ends) if (count != 2) return;
        totalEdges_ = edgeIndex.size();

        const size_t systemCount = systems.atoms.size();
        systemEdges_.assign(systemCount, {});
        localEdge_.assign(edgeList.size(), -1);
        for (size_t edge = 0; edge < edgeList.size(); ++edge) {
            const int s = ofAtom_[edgeList[edge].first];
            if (s < 0 || ofAtom_[edgeList[edge].second] != s) continue;
            localEdge_[edge] = static_cast<int>(systemEdges_[s].size());
            systemEdges_[s].push_back(static_cast<int>(edge));
        }
        lowestAtom_.assign(systemCount, -1);
        for (int atom = n_ - 1; atom >= 0; --atom)
            if (ofAtom_[atom] >= 0) lowestAtom_[ofAtom_[atom]] = atom;
        words_.resize(systemCount);
        systemRank_.resize(systemCount);
        systemEdgeBase_.assign(systemCount + 1, 0);
        for (size_t s = 0; s < systemCount; ++s) {
            words_[s] = (systemEdges_[s].size() + 63) >> 6;
            systemRank_[s] = static_cast<int>(systemEdges_[s].size()) - systems.atoms[s] + 1;
            systemEdgeBase_[s + 1] = systemEdgeBase_[s] + systemEdges_[s].size();
        }
        slotBond_.assign(slotEdge.size(), -1);
        cost_.resize(static_cast<size_t>(n_));
        for (int u = 0; u < n_; ++u) {
            const auto& row = g.neighbors[u];
            cost_[u] = 1 + static_cast<uint64_t>(row.size());
            for (size_t i = 0; i < row.size(); ++i) {
                const int slot = slotStart_[u] + static_cast<int>(i);
                if (ofAtom_[u] >= 0 && ofAtom_[row[i]] == ofAtom_[u]) slotBond_[slot] = localEdge_[slotEdge[slot]];
            }
        }
        // The pairs of the pairwise search, in its order, and the last of them.
        for (int v = 0; v < n_; ++v) {
            const size_t count = pairCountOf(v);
            if (count == 0) continue;
            pairCount_ += count;
            lastRoot_ = v;
            lastIndex_ = count - 1;
        }
        notes_.pairs = pairCount_;
        const size_t atoms = static_cast<size_t>(n_);
        par_.assign(atoms, UNSEEN);
        parBond_.assign(atoms, -1);
        depth_.assign(atoms, 0);
        rank_.assign(atoms, 0);
        pre_.assign(atoms, 0);
        size_.assign(atoms, 0);
        top_.assign(atoms, -1);
        subtree_.assign(atoms, 0);
        uf_.assign(atoms, 0);
        ufRank_.assign(atoms, 0);
        ufTop_.assign(atoms, 0);
        coverFrom_.assign(atoms, -1);
        coverTo_.assign(atoms, -1);
        coverBond_.assign(atoms, -1);
        b2Low_.assign(atoms, 0);
        b2High_.assign(atoms, 0);
        b2Count_.assign(atoms, 0);
        baseParent_.assign(atoms, -1);
        baseVisited_.assign(atoms, 0);
        marks_.assign(atoms, 0);
        applicable_ = true;
    }

    /** False when the graph has no ring systems or a neighbour list that is not mutual. */
    bool applicable() const noexcept { return applicable_; }

    /**
     * Reads every pair and decides whether the pairwise search falls back: it does exactly when
     * its work before the last pair exceeds the limit, since it tests the limit only before each
     * pair and its work never decreases. The sums leave out the last pair, and brackets that
     * leave the answer open are replayed.
     */
    bool decideFallback() {
        bool fallBack;
        if (!runPass(true, false)) fallBack = true;
        else if (!options_.replayAllForTesting && unbounded_ == 0 && high_ <= workLimit_) fallBack = false;
        else if (!options_.replayAllForTesting && low_ > workLimit_) fallBack = true;
        else fallBack = replay();
        notes_.baseWorkLow = low_;
        notes_.baseWorkHigh = high_;
        notes_.baseWorkHighUnbounded = unbounded_ != 0 || !passDone_;
        notes_.baseWorkExact = passDone_ && unbounded_ == 0 && low_ == high_;
        notes_.fellBack = fallBack;
        return fallBack;
    }

    /**
     * Reads the pairs left after a fallback, within the search's own work and storage bounds,
     * building the tree of the root it resumes again.
     */
    bool completeCandidates() {
        if (!passDone_) {
            treeStale_ = true;
            if (!runPass(false, true)) return false;
        }
        return withinLiftBudget();
    }

    /** Whether the search keeps within the bounds that finishing a refused search must keep. */
    bool withinLiftBudget() const noexcept {
        return notes_.searchWork <= workLimit_ && notes_.storedKeyWords <= RING_LIFT_KEY_WORDS;
    }

    /**
     * The two-phase selection of the pairwise search over the candidates found, in ascending
     * size and, within a size, in the order found: a candidate independent of strictly shorter
     * rows is relevant, and a relevant one independent of every row joins the basis until it has
     * cycleRank rows. Rows of one ring system never meet another system's bonds, so each system
     * keeps its own rows; once a system has a full set, its longer candidates are dependent.
     */
    void select(int cycleRank, std::vector<std::vector<int>>& sssr, std::vector<std::vector<int>>& relevant) {
        const size_t count = candidateSize_.size();
        int longest = 0;
        for (int size : candidateSize_) longest = std::max(longest, size);
        std::vector<size_t> start(static_cast<size_t>(longest) + 2, 0), order(count);
        for (int size : candidateSize_) ++start[static_cast<size_t>(size) + 1];
        for (size_t i = 1; i < start.size(); ++i) start[i] += start[i - 1];
        for (size_t c = 0; c < count; ++c) order[start[static_cast<size_t>(candidateSize_[c])]++] = c;

        const size_t systemCount = systemEdges_.size();
        std::vector<int> rowByPivot(systemEdgeBase_[systemCount], -1);
        std::vector<std::vector<uint64_t>> rows(systemCount);
        std::vector<std::vector<int>> rowSizes(systemCount);
        std::vector<int> longestRow(systemCount, 0);
        std::vector<uint64_t> bits;
        for (size_t c : order) {
            const int s = candidateSystem_[c], size = candidateSize_[c];
            const int filled = static_cast<int>(rowSizes[s].size());
            if (filled == systemRank_[s] && size > longestRow[s]) continue;
            const size_t words = words_[s];
            const auto key = candidateKeys_.begin() + static_cast<std::ptrdiff_t>(candidateKeyStart_[c]);
            bits.assign(key, key + static_cast<std::ptrdiff_t>(words));
            notes_.searchWork += words;
            int* pivots = rowByPivot.data() + systemEdgeBase_[s];
            // Reads the set bits upward, adding the row of each bit it may use.
            auto reduce = [&](bool shorterOnly) {
                for (size_t word = 0; word < words; ++word) {
                    uint64_t pending = bits[word];
                    while (pending != 0) {
                        const int bit = ctz64(pending);
                        const int row = pivots[(word << 6) + static_cast<size_t>(bit)];
                        if (row < 0 || (shorterOnly && rowSizes[s][row] >= size)) {
                            pending &= pending - 1;
                            continue;
                        }
                        const uint64_t* added = rows[s].data() + static_cast<size_t>(row) * words;
                        for (size_t other = word; other < words; ++other) bits[other] ^= added[other];
                        notes_.searchWork += words;
                        pending = bit == 63 ? 0 : bits[word] & (~uint64_t(0) << (bit + 1));
                    }
                }
            };
            reduce(true);
            bool nonzero = false;
            for (uint64_t word : bits) if (word != 0) { nonzero = true; break; }
            if (!nonzero) continue;
            relevant.push_back(atomsOf(c));
            reduce(false);
            int pivot = -1;
            for (size_t word = 0; word < words && pivot < 0; ++word)
                if (bits[word] != 0) pivot = static_cast<int>((word << 6) + static_cast<size_t>(ctz64(bits[word])));
            if (pivot >= 0 && static_cast<int>(sssr.size()) < cycleRank) {
                pivots[pivot] = filled;
                rows[s].insert(rows[s].end(), bits.begin(), bits.end());
                rowSizes[s].push_back(size);
                longestRow[s] = std::max(longestRow[s], size);
                sssr.push_back(atomsOf(c));
            }
        }
    }

    // ----- Inspection, for tests -----

    /** Builds the tree of root, as the pass does. */
    void buildTreeForTesting(int root) { buildTree(root, true, false); }
    /** Reads one bond of the root's system from the last tree built, walking its climb. */
    PairReading readPairForTesting(int edge) { return readPair(edge, Climb::Always); }
    /** The paths to the bond's first and second atoms of the reading just made. */
    void pathsForTesting(const PairReading& reading, std::vector<int>& first, std::vector<int>& second) const {
        pathsOf(reading, first, second);
    }
    /** The closing half-bond assigned to an atom of the last tree, or (-1, -1) for its root. */
    std::pair<int, int> coverForTesting(int atom) const { return {coverFrom_[atom], coverTo_[atom]}; }
    /** The closing half-bonds of the last tree, in its depth-first order. */
    std::vector<std::pair<int, int>> closingForTesting() const {
        std::vector<std::pair<int, int>> list;
        for (const auto& half : closing_) list.push_back({half.from, half.to});
        return list;
    }
    int depthForTesting(int atom) const { return depth_[atom]; }
    bool inSubtreeForTesting(int ancestor, int atom) const { return isAncestor(ancestor, atom); }
    const std::vector<int>& systemEdgesForTesting(int system) const { return systemEdges_[system]; }

private:
    struct HalfBond { int from, to, bond, key; };
    struct Frame { int atom; size_t next; };
    struct Slot { uint64_t hash; int system; int head; int tail; };
    enum class Climb { Never, WhenValid, Always };

    static constexpr int UNSEEN = -2;

    size_t pairCountOf(int v) const {
        const int s = ofAtom_[v];
        if (s < 0) return 0;
        if (!simpleRing_[s]) return systemEdges_[s].size();
        return v == lowestAtom_[s] && !systemEdges_[s].empty() ? 1 : 0;
    }

    bool isAncestor(int ancestor, int atom) const noexcept {
        return pre_[ancestor] <= pre_[atom] && pre_[atom] < pre_[ancestor] + size_[ancestor];
    }

    uint64_t levelSum(int64_t level) const noexcept {
        if (level < 0) return 0;
        if (level > maxDepth_) return levelPrefix_[static_cast<size_t>(maxDepth_)];
        return levelPrefix_[static_cast<size_t>(level)];
    }

    /** The tree of root over its system with its numbering, covers and tree-path bond sets. */
    void buildTree(int root, bool wantPathSets, bool charge) {
        const int s = ofAtom_[root];
        for (int u : order_) par_[u] = UNSEEN;
        order_.clear();
        root_ = root;
        rootSystem_ = s;
        par_[root] = -1;
        parBond_[root] = -1;
        depth_[root] = 0;
        order_.push_back(root);
        uint64_t atomWork = 0;
        for (size_t head = 0; head < order_.size(); ++head) {
            const int u = order_[head];
            rank_[u] = static_cast<int>(head);
            atomWork += cost_[u];
            const auto& row = g_.neighbors[u];
            const int* bonds = slotBond_.data() + slotStart_[u];
            for (size_t i = 0; i < row.size(); ++i) {
                const int w = row[i];
                if (bonds[i] < 0 || par_[w] != UNSEEN) continue;
                par_[w] = u;
                parBond_[w] = bonds[i];
                depth_[w] = depth_[u] + 1;
                order_.push_back(w);
            }
        }
        const size_t atoms = order_.size();
        prefix_.resize(atoms);
        uint64_t sum = 0;
        for (size_t i = 0; i < atoms; ++i) prefix_[i] = sum += cost_[order_[i]];
        maxDepth_ = depth_[order_.back()];
        levelPrefix_.assign(static_cast<size_t>(maxDepth_) + 1, 0);
        for (int u : order_) levelPrefix_[static_cast<size_t>(depth_[u])] += cost_[u];
        for (size_t level = 1; level < levelPrefix_.size(); ++level) levelPrefix_[level] += levelPrefix_[level - 1];
        for (int u : order_) { size_[u] = 1; subtree_[u] = cost_[u]; }
        for (size_t i = atoms; i-- > 1; ) {
            const int u = order_[i];
            size_[par_[u]] += size_[u];
            subtree_[par_[u]] += subtree_[u];
        }
        top_[root] = -1;
        for (size_t i = 1; i < atoms; ++i) {
            const int u = order_[i];
            top_[u] = par_[u] == root ? u : top_[par_[u]];
        }

        // Depth-first numbering over the tree, children in list order; every other bond of the
        // system read from an atom is a closing half-bond, listed in that order.
        int counter = 0;
        pre_[root] = counter++;
        closing_.clear();
        stack_.clear();
        stack_.push_back({root, 0});
        while (!stack_.empty()) {
            Frame& frame = stack_.back();
            const int x = frame.atom;
            const auto& row = g_.neighbors[x];
            if (frame.next == row.size()) { stack_.pop_back(); continue; }
            const size_t i = frame.next++;
            const int bond = slotBond_[static_cast<size_t>(slotStart_[x]) + i];
            const int w = row[i];
            if (bond < 0 || w == par_[x]) continue;
            if (par_[w] == x) {
                pre_[w] = counter++;
                stack_.push_back({w, 0});
            } else {
                closing_.push_back({x, w, bond, depth_[x] + depth_[w] + 1});
            }
        }

        // Closing half-bonds by path length, keeping their order within a length.
        counts_.assign(2 * static_cast<size_t>(maxDepth_) + 3, 0);
        for (const auto& half : closing_) ++counts_[static_cast<size_t>(half.key) + 1];
        for (size_t i = 1; i < counts_.size(); ++i) counts_[i] += counts_[i - 1];
        sorted_.resize(closing_.size());
        for (const auto& half : closing_) sorted_[static_cast<size_t>(counts_[static_cast<size_t>(half.key)]++)] = half;

        // Each atom's cover is the first half-bond in that order from its subtree to outside it:
        // a half-bond covers the uncovered atoms from its end up to the first ancestor of its
        // start. Covered atoms join their parent's set, whose top is its nearest uncovered atom;
        // ranks keep the sets shallow, and on equal ranks the child's root stays the root.
        for (int u : order_) {
            uf_[u] = u;
            ufRank_[u] = 0;
            ufTop_[u] = u;
            coverFrom_[u] = coverTo_[u] = coverBond_[u] = -1;
        }
        uint64_t unionWork = 0;
        auto find = [&](int u) {
            while (uf_[u] != u) { uf_[u] = uf_[uf_[u]]; u = uf_[u]; ++unionWork; }
            return u;
        };
        for (const auto& half : sorted_) {
            int u = ufTop_[find(half.to)];
            while (!isAncestor(u, half.from)) {
                coverFrom_[u] = half.from;
                coverTo_[u] = half.to;
                coverBond_[u] = half.bond;
                int keep = find(u), join = find(par_[u]);
                const int top = ufTop_[join];
                if (ufRank_[keep] < ufRank_[join]) std::swap(keep, join);
                uf_[join] = keep;
                if (ufRank_[keep] == ufRank_[join]) ++ufRank_[keep];
                ufTop_[keep] = top;
                ++unionWork;
                u = top;
            }
        }

        // The bonds of each tree path, by tree rank, when they fit.
        const size_t words = words_[s];
        pathSets_ = wantPathSets && static_cast<uint64_t>(atoms) * words <= RING_TREE_PATH_SET_WORDS;
        if (pathSets_) {
            pathSet_.assign(atoms * words, 0);
            for (size_t i = 1; i < atoms; ++i) {
                const int u = order_[i];
                const uint64_t* from = pathSet_.data() + static_cast<size_t>(rank_[par_[u]]) * words;
                uint64_t* to = pathSet_.data() + i * words;
                std::copy(from, from + words, to);
                to[static_cast<size_t>(parBond_[u]) >> 6] ^= uint64_t(1) << (parBond_[u] & 63);
            }
        }
        if (charge)
            notes_.searchWork += 2 * atomWork + closing_.size() + 2 * static_cast<uint64_t>(maxDepth_) + 1
                + unionWork + (pathSets_ ? static_cast<uint64_t>(atoms) * words : 0);
    }

    /** One bond of the root's system: its kind, validity and the count of its pair search. */
    PairReading readPair(int edge, Climb climb) {
        PairReading r;
        r.edge = edge;
        const int a = edgeList_[edge].first, b = edgeList_[edge].second;
        if (par_[a] != b && par_[b] != a) {
            // The pair search makes the tree's discoveries and stops at the later end's parent.
            r.kind = PairKind::NonTree;
            r.valid = top_[a] != top_[b];
            const int later = rank_[a] > rank_[b] ? a : b;
            r.searchLow = r.searchHigh = prefix_[static_cast<size_t>(rank_[par_[later]])];
            r.pathTerm = static_cast<uint64_t>(depth_[a]) + static_cast<uint64_t>(depth_[b]) + 2;
            return r;
        }
        const int child = par_[a] == b ? a : b, parent = child == a ? b : a;
        const int from = coverFrom_[child], to = coverTo_[child];
        r.child = child;
        r.parent = parent;
        r.coverFrom = from;
        r.coverTo = to;
        const int closing = depth_[from] + depth_[to] + 1;
        r.valid = parent == root_ || top_[from] != top_[parent];
        r.pathTerm = static_cast<uint64_t>(closing) + 1;
        climb_.assign(1, to);
        climbBond_.clear();
        if (to == child) {
            // The search reads the tree's order up to the cover's start, without the child.
            const int step = depth_[from] - depth_[child];
            r.kind = step < 0 ? PairKind::ClassA : step == 0 ? PairKind::ClassB1 : PairKind::ClassB1Plus;
            r.searchLow = r.searchHigh = prefix_[static_cast<size_t>(rank_[from])]
                - (rank_[child] < rank_[from] ? cost_[child] : 0);
            return r;
        }
        // The child's discoverer lies on level `level` of the search: every lower level is read,
        // and so is the climb, but nothing above that level and not the child itself.
        r.kind = PairKind::ClassB2;
        const int64_t level = static_cast<int64_t>(closing) - depth_[child] - 1;
        const uint64_t below = levelSum(level - 1);
        r.searchLow = below > subtree_[child] ? below - subtree_[child] : 0;
        r.searchHigh = levelSum(level) - cost_[child];
        if (climb == Climb::Never || (climb == Climb::WhenValid && !r.valid)) return r;
        uint64_t climbWork = 0;
        for (int u = to; u != child; ) {
            climbWork += cost_[u];
            const auto& row = g_.neighbors[u];
            const int* bonds = slotBond_.data() + slotStart_[u];
            int next = -1, bond = -1;
            for (size_t i = 0; i < row.size(); ++i) {
                const int y = row[i];
                if (bonds[i] < 0 || depth_[y] != depth_[u] - 1 || !isAncestor(child, y)) continue;
                next = y;
                bond = bonds[i];
                break;
            }
            climbBond_.push_back(bond);
            climb_.push_back(next);
            u = next;
        }
        if (r.valid) {
            r.climbWork = climbWork;
            r.searchLow += climbWork;
        }
        return r;
    }

    void appendTreePath(int atom, std::vector<int>& path) const {
        const size_t start = path.size();
        for (int u = atom; u != -1; u = par_[u]) path.push_back(u);
        std::reverse(path.begin() + static_cast<std::ptrdiff_t>(start), path.end());
    }

    /** The pair search's paths to the bond's first and second atoms, after reading the pair. */
    void pathsOf(const PairReading& r, std::vector<int>& first, std::vector<int>& second) const {
        first.clear();
        second.clear();
        const int a = edgeList_[r.edge].first;
        if (r.kind == PairKind::NonTree) {
            appendTreePath(a, first);
            appendTreePath(edgeList_[r.edge].second, second);
            return;
        }
        std::vector<int>& toChild = r.child == a ? first : second;
        std::vector<int>& toParent = r.child == a ? second : first;
        appendTreePath(r.coverFrom, toChild);
        toChild.insert(toChild.end(), climb_.begin(), climb_.end());
        appendTreePath(r.parent, toParent);
    }

    /** Stores a valid pair's cycle unless an earlier pair of its system gave the same bonds. */
    void storeCandidate(const PairReading& r) {
        const int s = rootSystem_;
        const size_t words = words_[s];
        const uint64_t length = r.pathTerm - 1;
        key_.assign(words, 0);
        auto flip = [&](int bond) { key_[static_cast<size_t>(bond) >> 6] ^= uint64_t(1) << (bond & 63); };
        auto addPath = [&](int atom) {
            if (pathSets_) {
                const uint64_t* set = pathSet_.data() + static_cast<size_t>(rank_[atom]) * words;
                for (size_t i = 0; i < words; ++i) key_[i] ^= set[i];
            } else {
                for (int u = atom; u != root_; u = par_[u]) flip(parBond_[u]);
            }
        };
        notes_.searchWork += pathSets_ ? words : length + words;
        if (r.kind == PairKind::NonTree) {
            addPath(edgeList_[r.edge].first);
            addPath(edgeList_[r.edge].second);
        } else {
            addPath(r.coverFrom);
            addPath(r.parent);
            flip(coverBond_[r.child]);
            for (int bond : climbBond_) flip(bond);
        }
        flip(localEdge_[r.edge]);
        uint64_t compared = 0;
        const bool stored = storeKey(s, length, compared, [&](std::vector<int>& atoms) {
            pathsOf(r, path1_, path2_);
            atoms.insert(atoms.end(), path1_.begin(), path1_.end());
            for (size_t i = path2_.size(); i-- > 1; ) atoms.push_back(path2_[i]);
        });
        notes_.searchWork += compared * words + (stored ? length + words : 0);
    }

    /**
     * Stores key_ as a candidate of system s unless the system already has it, comparing only
     * the stored bond sets with an equal hash, in the order stored; write appends the
     * candidate's atoms.
     */
    template<class Write>
    bool storeKey(int s, uint64_t length, uint64_t& compared, Write write) {
        const size_t words = words_[s];
        const uint64_t hash = ringCandidateKeyHash(key_.data(), words);
        if (table_.empty()) table_.assign(64, Slot{0, -1, -1, -1});
        auto home = [&](uint64_t value, int system) {
            return static_cast<size_t>(value ^ (static_cast<uint64_t>(system) * UINT64_C(0x9E3779B97F4A7C15)))
                & (table_.size() - 1);
        };
        size_t at = home(hash, s);
        while (table_[at].head >= 0 && !(table_[at].hash == hash && table_[at].system == s))
            at = (at + 1) & (table_.size() - 1);
        for (int c = table_[at].head; c >= 0; c = nextSameHash_[c]) {
            ++compared;
            const auto stored = candidateKeys_.begin() + static_cast<std::ptrdiff_t>(candidateKeyStart_[c]);
            if (std::equal(key_.begin(), key_.end(), stored)) return false;
        }
        const int c = static_cast<int>(candidateSize_.size());
        candidateSystem_.push_back(s);
        candidateSize_.push_back(static_cast<int>(length));
        candidateAtomStart_.push_back(candidateAtoms_.size());
        write(candidateAtoms_);
        candidateKeyStart_.push_back(candidateKeys_.size());
        candidateKeys_.insert(candidateKeys_.end(), key_.begin(), key_.end());
        notes_.storedKeyWords += words;
        ++notes_.candidates;
        nextSameHash_.push_back(-1);
        if (table_[at].head < 0) {
            table_[at] = Slot{hash, s, c, c};
            ++tableUsed_;
        } else {
            nextSameHash_[static_cast<size_t>(table_[at].tail)] = c;
            table_[at].tail = c;
        }
        if (2 * tableUsed_ > table_.size()) {
            std::vector<Slot> old(table_.size() * 2, Slot{0, -1, -1, -1});
            old.swap(table_);
            for (const auto& slot : old) {
                if (slot.head < 0) continue;
                size_t place = home(slot.hash, slot.system);
                while (table_[place].head >= 0) place = (place + 1) & (table_.size() - 1);
                table_[place] = slot;
            }
        }
        return true;
    }

    std::vector<int> atomsOf(size_t c) const {
        const auto first = candidateAtoms_.begin() + static_cast<std::ptrdiff_t>(candidateAtomStart_[c]);
        return std::vector<int>(first, first + candidateSize_[c]);
    }

    /**
     * The pairwise search of one pair, as it is written there: from start, avoiding the bond of
     * the two targets, inside their ring system; returns the atoms and neighbours it reads.
     */
    uint64_t baseSearch(int start, int firstTarget, int secondTarget,
                        std::vector<int>* firstPath, std::vector<int>* secondPath) {
        uint64_t work = 0;
        bool firstReached = start == firstTarget, secondReached = start == secondTarget;
        if (firstReached && secondReached) {
            if (firstPath) firstPath->push_back(start);
            if (secondPath) secondPath->push_back(start);
            return work;
        }
        if (++baseStamp_ == 0) { std::fill(baseVisited_.begin(), baseVisited_.end(), 0u); baseStamp_ = 1; }
        const int system = ofAtom_[firstTarget];
        baseVisited_[start] = baseStamp_;
        baseParent_[start] = -1;
        baseQueue_.assign(1, start);
        for (size_t head = 0; head < baseQueue_.size() && !(firstReached && secondReached); ++head) {
            const int u = baseQueue_[head];
            work += 1 + g_.neighbors[u].size();
            for (int next : g_.neighbors[u]) {
                if ((u == firstTarget && next == secondTarget) || (u == secondTarget && next == firstTarget)) continue;
                if (baseVisited_[next] == baseStamp_ || ofAtom_[next] != system) continue;
                baseVisited_[next] = baseStamp_;
                baseParent_[next] = u;
                if (next == firstTarget) firstReached = true;
                if (next == secondTarget) secondReached = true;
                if (firstReached && secondReached) break;
                baseQueue_.push_back(next);
            }
        }
        if (!firstReached || !firstPath) return work;
        for (int c = firstTarget; c != -1; c = baseParent_[c]) firstPath->push_back(c);
        std::reverse(firstPath->begin(), firstPath->end());
        if (!secondReached || !secondPath) return work;
        for (int c = secondTarget; c != -1; c = baseParent_[c]) secondPath->push_back(c);
        std::reverse(secondPath->begin(), secondPath->end());
        return work;
    }

    /** The one pair of a system that is a simple ring, read by the pairwise search itself. */
    uint64_t readSimplePair(int root, int edge) {
        const int a = edgeList_[edge].first, b = edgeList_[edge].second;
        path1_.clear();
        path2_.clear();
        uint64_t charge = baseSearch(root, a, b, &path1_, &path2_);
        if (path1_.empty() || path2_.empty()) return charge;
        charge += path1_.size() + path2_.size();
        if (++markStamp_ == 0) { std::fill(marks_.begin(), marks_.end(), 0u); markStamp_ = 1; }
        for (size_t i = 1; i < path1_.size(); ++i) marks_[path1_[i]] = markStamp_;
        for (size_t i = 1; i < path2_.size(); ++i) if (marks_[path2_[i]] == markStamp_) return charge;
        const uint64_t length = path1_.size() + path2_.size() - 1;
        if (length < 3) return charge;
        charge += totalEdges_;
        const int s = ofAtom_[a];
        key_.assign(words_[s], 0);
        auto flip = [&](int bond) { key_[static_cast<size_t>(bond) >> 6] ^= uint64_t(1) << (bond & 63); };
        for (const auto* path : {&path1_, &path2_})
            for (size_t i = 0; i + 1 < path->size(); ++i) {
                const auto& row = g_.neighbors[(*path)[i]];
                const auto slot = std::find(row.begin(), row.end(), (*path)[i + 1]) - row.begin();
                flip(slotBond_[static_cast<size_t>(slotStart_[(*path)[i]] + slot)]);
            }
        flip(localEdge_[edge]);
        uint64_t compared = 0;
        if (storeKey(s, length, compared, [&](std::vector<int>& atoms) {
                atoms.insert(atoms.end(), path1_.begin(), path1_.end());
                for (size_t i = path2_.size(); i-- > 1; ) atoms.push_back(path2_[i]);
            }))
            notes_.searchWork += length;
        return charge;
    }

    /**
     * Reads pairs from the cursor in the pairwise order. With decide, each pair before the last
     * adds its count or bracket to the sums, and the pass stops once the lower sum exceeds the
     * limit; with lift, it stops before a pair once the search has left its own bounds. False
     * when it stopped.
     */
    bool runPass(bool decide, bool lift) {
        for (; cursorRoot_ < n_; ++cursorRoot_, cursorIndex_ = 0) {
            const int root = cursorRoot_;
            const int s = ofAtom_[root];
            if (s < 0) continue;
            const auto& edges = systemEdges_[s];
            if (simpleRing_[s]) {
                if (root != lowestAtom_[s] || edges.empty() || cursorIndex_ > 0) continue;
                if (lift && !withinLiftBudget()) return false;
                cursorIndex_ = 1;
                const uint64_t charge = readSimplePair(root, edges.front());
                notes_.searchWork += charge;
                ++notes_.simpleRing;
                ++notes_.pairsRead;
                if (decide && !account(root, PairKind::SimpleRing, charge, charge, 0)) return false;
                continue;
            }
            if (root_ != root || treeStale_) {
                buildTree(root, true, true);
                treeStale_ = false;
            }
            while (cursorIndex_ < edges.size()) {
                if (lift && !withinLiftBudget()) return false;
                const int edge = edges[cursorIndex_++];
                const PairReading r = readPair(edge, Climb::WhenValid);
                ++notes_.pairsRead;
                ++notes_.searchWork;
                switch (r.kind) {
                    case PairKind::NonTree: ++notes_.nonTree; break;
                    case PairKind::ClassA: ++notes_.classA; break;
                    case PairKind::ClassB1: ++notes_.classB1; break;
                    case PairKind::ClassB1Plus: ++notes_.classB1Plus; break;
                    case PairKind::ClassB2: ++notes_.classB2; break;
                    case PairKind::SimpleRing: break;
                }
                if (r.valid) {
                    notes_.searchWork += r.climbWork;
                    storeCandidate(r);
                }
                const uint64_t fixed = r.pathTerm + (r.valid ? totalEdges_ : 0);
                if (decide && !account(root, r.kind, r.searchLow, r.searchHigh, fixed)) return false;
            }
        }
        passDone_ = true;
        return true;
    }

    /** Adds a pair before the last to the sums; false once the lower sum exceeds the limit. */
    bool account(int root, PairKind kind, uint64_t searchLow, uint64_t searchHigh, uint64_t fixed) {
        if (++pairsSeen_ == pairCount_) return true;   // the last pair's work reaches no test
        if (kind == PairKind::ClassB2) {
            if (options_.loosenBracketsForTesting) {
                searchLow = 0;
                searchHigh = 0;
                ++unbounded_;
            }
            b2Low_[root] += searchLow;
            b2High_[root] += searchHigh;
            ++b2Count_[root];
        }
        low_ += fixed + searchLow;
        high_ += fixed + searchHigh;
        return options_.replayAllForTesting || low_ <= workLimit_;
    }

    /**
     * Replays the bracketed pairs root by root, replacing each root's brackets by their exact
     * sum, until the sums decide the fallback; once every bracket is replayed they are exact.
     */
    bool replay() {
        treeStale_ = true;
        for (int root = 0; root < n_; ++root) {
            if (b2Count_[root] == 0) continue;
            buildTree(root, false, false);
            const auto& edges = systemEdges_[ofAtom_[root]];
            uint64_t exact = 0;
            for (size_t i = 0; i < edges.size(); ++i) {
                if (root == lastRoot_ && i == lastIndex_) continue;
                const int edge = edges[i];
                if (readPair(edge, Climb::Never).kind != PairKind::ClassB2) continue;
                const uint64_t charge = baseSearch(root, edgeList_[edge].first, edgeList_[edge].second, nullptr, nullptr);
                exact += charge;
                notes_.searchWork += charge;
                ++notes_.tier2Replays;
            }
            low_ = low_ - b2Low_[root] + exact;
            high_ = high_ - b2High_[root] + exact;
            if (options_.loosenBracketsForTesting) unbounded_ -= b2Count_[root];
            b2Count_[root] = 0;
            b2Low_[root] = b2High_[root] = 0;
            if (options_.replayAllForTesting) continue;
            if (unbounded_ == 0 && high_ <= workLimit_) return false;
            if (low_ > workLimit_) return true;
        }
        return low_ > workLimit_;
    }

    const MolGraph& g_;
    const std::vector<int>& ofAtom_;
    const std::vector<uint8_t>& simpleRing_;
    const std::vector<std::pair<int, int>>& edgeList_;
    const uint64_t workLimit_;
    const RingSearchOptions options_;
    RingSearchNotes& notes_;
    const int n_;
    bool applicable_ = false;
    uint64_t totalEdges_ = 0;

    std::vector<std::vector<int>> systemEdges_;   // global bonds of each system, in global order
    std::vector<int> localEdge_;                  // each bond's place in its system, or -1
    std::vector<int> lowestAtom_;
    std::vector<size_t> words_;                   // words of one bond set of each system
    std::vector<int> systemRank_;
    std::vector<size_t> systemEdgeBase_;
    std::vector<int> slotStart_;                  // neighbour slots of each atom
    std::vector<int> slotBond_;                   // a slot's bond in its system, or -1 outside it
    std::vector<uint64_t> cost_;                  // 1 + degree: what the pair search reads per atom

    uint64_t pairCount_ = 0;
    int lastRoot_ = -1;
    size_t lastIndex_ = 0;
    int cursorRoot_ = 0;
    size_t cursorIndex_ = 0;
    uint64_t pairsSeen_ = 0;
    bool passDone_ = false;
    bool treeStale_ = false;

    uint64_t low_ = 0, high_ = 0, unbounded_ = 0;
    std::vector<uint64_t> b2Low_, b2High_;
    std::vector<uint32_t> b2Count_;

    int root_ = -1, rootSystem_ = -1, maxDepth_ = 0;
    bool pathSets_ = false;
    std::vector<int> order_, par_, parBond_, depth_, rank_, pre_, size_, top_;
    std::vector<uint64_t> prefix_, levelPrefix_, subtree_;
    std::vector<HalfBond> closing_, sorted_;
    std::vector<int> counts_;
    std::vector<Frame> stack_;
    std::vector<int> uf_, ufRank_, ufTop_, coverFrom_, coverTo_, coverBond_;
    std::vector<uint64_t> pathSet_;
    std::vector<int> climb_, climbBond_, path1_, path2_;

    std::vector<int> candidateSystem_, candidateSize_, candidateAtoms_, nextSameHash_;
    std::vector<size_t> candidateAtomStart_, candidateKeyStart_;
    std::vector<uint64_t> candidateKeys_, key_;
    std::vector<Slot> table_;
    size_t tableUsed_ = 0;

    std::vector<int> baseParent_, baseQueue_;
    std::vector<uint32_t> baseVisited_, marks_;
    uint32_t baseStamp_ = 0, markStamp_ = 0;
};

} // namespace detail_rf
} // namespace smsd

#endif // SMSD_DETAIL_RING_CANDIDATES_HPP

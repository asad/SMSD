/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Label-histogram similarity upper bound for screening, O(V+E).
 * After Raymond, Gardiner & Willett, J. Comput.-Aided Mol. Des. 2002.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <unordered_map>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace smsd {

namespace screening_detail {
/**
 * The label two atoms must share to pair in the bound: the element, the aromatic flag under
 * strict aromaticity and the ring flag when ring atoms match ring atoms only. Tautomerism
 * changes proton positions and bond orders, not elements, so atomCompatible keeps elements apart
 * under tautomer-aware matching and the label does too. The Java MCS.similarityUpperBound uses
 * the same label.
 */
inline int ubLabel(const MolGraph& g, int i, const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    const int z = g.atomicNum[i];
    // Shifted as unsigned, so that an out-of-range negative atomic number keeps its label without undefined behaviour.
    int key = static_cast<int>(static_cast<unsigned>(z) << 2);
    if (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT) key |= g.aromatic[i] ? 2 : 0;
    if (C.ringMatchesRingOnly) key |= g.ring[i] ? 1 : 0;
    return key;
}

/** The number of label rules ubLabel applies, which labelMode numbers from zero. */
constexpr int labelModes = 5;

/** The label rule ubLabel applies under the chemistry: untyped, or one of four typed rules. */
inline int labelMode(const ChemOptions& C) {
    if (!C.matchAtomType) return 0;
    return 1 + (C.aromaticityMode == ChemOptions::AromaticityMode::STRICT ? 2 : 0) + (C.ringMatchesRingOnly ? 1 : 0);
}

/**
 * The bound for a label overlap, the sum over labels of the smaller of the two counts, between
 * graphs of n1 and n2 atoms: the overlap is capped by the smaller graph and gives the similarity
 * overlap / (n1 + n2 - overlap). Every CPU screen scores through this, so their scores agree bit for bit.
 */
inline double overlapSimilarity(int overlap, int n1, int n2) {
    if (n1 == 0 || n2 == 0) return 0.0;
    const int atomOverlap = std::min(overlap, std::min(n1, n2));
    if (atomOverlap == 0) return 0.0;
    return static_cast<double>(atomOverlap) / (n1 + n2 - atomOverlap);
}

inline int labelFrequencyUpperBound(const MolGraph& g1, const MolGraph& g2,
                                     const ChemOptions& C) {
    if (g1.n == 0 || g2.n == 0) return 0;
    std::unordered_map<int,int> freq1, freq2;
    for (int i = 0; i < g1.n; ++i) freq1[ubLabel(g1, i, C)]++;
    for (int j = 0; j < g2.n; ++j) freq2[ubLabel(g2, j, C)]++;
    int ub = 0;
    for (auto& [lab, cnt] : freq1) {
        auto it = freq2.find(lab);
        if (it != freq2.end()) ub += std::min(cnt, it->second);
    }
    return std::min(ub, std::min(g1.n, g2.n));
}

class TargetLabels;

/**
 * Query label counts held in a dense table, so one query can be bounded against many targets
 * without building hash maps per target. Gives the same value as labelFrequencyUpperBound.
 * The query graph must outlive this object.
 */
class QueryLabelCounts {
public:
    /** Labels at or above this limit, or below zero, use the hashed path. */
    static constexpr int denseLimit = 1 << 10;

    QueryLabelCounts(const MolGraph& query, const ChemOptions& chemistry)
        : query_(query), chemistry_(chemistry), counts_(denseLimit, 0) {
        for (int atom = 0; atom < query.n && dense_; ++atom) {
            const int label = ubLabel(query, atom, chemistry);
            if (label < 0 || label >= denseLimit) dense_ = false;
            else ++counts_[label];
        }
    }

    /** scratch must hold denseLimit zeros; it is left that way for the next target. */
    int upperBound(const MolGraph& target, std::vector<int>& scratch) const {
        if (query_.n == 0 || target.n == 0) return 0;
        if (!dense_) return labelFrequencyUpperBound(query_, target, chemistry_);
        int overlap = 0, atom = 0;
        for (; atom < target.n; ++atom) {
            const int label = ubLabel(target, atom, chemistry_);
            if (label < 0 || label >= denseLimit) break;
            // The first min(query, target) atoms of each label pair up.
            if (++scratch[label] <= counts_[label]) ++overlap;
        }
        const bool complete = atom == target.n;
        for (int counted = 0; counted < atom; ++counted) scratch[ubLabel(target, counted, chemistry_)] = 0;
        if (!complete) return labelFrequencyUpperBound(query_, target, chemistry_);
        return std::min(overlap, std::min(query_.n, target.n));
    }

    /** Same value as similarityUpperBound(query, target, chemistry). */
    double similarity(const MolGraph& target, std::vector<int>& scratch) const {
        if (query_.n == 0 || target.n == 0) return 0.0;
        return overlapSimilarity(upperBound(target, scratch), query_.n, target.n);
    }

    /** Whether every query label lies in the dense table. */
    bool dense() const { return dense_; }

    /**
     * The same value for a target whose histogram the table holds, read from the histogram. The
     * query must be dense and the table built under the same label rule.
     */
    double similarity(const TargetLabels& labels, int target) const;

private:
    const MolGraph& query_;
    ChemOptions chemistry_;
    std::vector<int> counts_;
    bool dense_ = true;
};

/**
 * Calls body(i, scratch) for every index i below N on up to numThreads threads (0 = the OpenMP
 * thread limit, which a larger request is capped at); scratch is the worker's table of denseLimit
 * zeros, which body must leave that way. There is at most one worker for every 256 indices, so
 * small libraries stay on the calling thread. When body throws, the exception of the first such index propagates,
 * as on one thread.
 */
template<class Body>
inline void forEachTarget(int N, int numThreads, const Body& body) {
#ifdef _OPENMP
    constexpr int minimumTargets = 256;
    const int requested = numThreads > 0 ? std::min(numThreads, omp_get_max_threads()) : omp_get_max_threads();
    const int workers = std::max(1, std::min(requested, N / minimumTargets));
    if (workers > 1) {
        // Each worker reuses one scratch table, made before the team starts, so no target allocates.
        std::vector<std::vector<int>> scratch(workers, std::vector<int>(QueryLabelCounts::denseLimit, 0));
        std::atomic<int> firstFailure{N};
        std::exception_ptr failure;
        #pragma omp parallel num_threads(workers)
        {
            std::vector<int>& table = scratch[omp_get_thread_num()];
            #pragma omp for schedule(dynamic, 128)
            for (int i = 0; i < N; ++i) {
                // Targets after a failure are skipped; every target before it is still visited.
                if (i > firstFailure.load(std::memory_order_relaxed)) continue;
                try {
                    body(i, table);
                } catch (...) {
                    #pragma omp critical(smsd_similarity_scores_failure)
                    if (i < firstFailure.load(std::memory_order_relaxed)) {
                        failure = std::current_exception();
                        firstFailure.store(i, std::memory_order_relaxed);
                    }
                }
            }
        }
        if (failure) std::rethrow_exception(failure);
        return;
    }
#else
    (void)numThreads;
#endif
    std::vector<int> scratch(QueryLabelCounts::denseLimit, 0);
    for (int i = 0; i < N; ++i) body(i, scratch);
}

/**
 * The label histogram of every target under one label rule, in flat arrays, so that a screen reads
 * each target's few distinct labels in order rather than its atom arrays from separate allocations.
 * A target with a label outside the dense table, or with an atom count that a count cannot hold
 * (negative, or above 65,535), is left to its graph. Built on up to numThreads threads, as forEachTarget runs.
 */
class TargetLabels {
public:
    struct LabelCount { std::uint16_t label, count; };

    TargetLabels(const std::vector<MolGraph>& targets, const ChemOptions& chemistry, int numThreads)
        : mode_(labelMode(chemistry)), atoms_(targets.size(), -1), offsets_(targets.size() + 1, 0) {
        const int N = static_cast<int>(targets.size());
        // Each target writes its histogram at a slot of its atom count; the slots are then packed.
        std::vector<std::size_t> slots(targets.size() + 1, 0);
        for (int i = 0; i < N; ++i) slots[i + 1] = slots[i] + (fits(targets[i]) ? targets[i].n : 0);
        std::vector<LabelCount> staged(slots[N]);
        std::vector<int> distinct(targets.size(), 0);
        forEachTarget(N, numThreads, [&](int i, std::vector<int>& scratch) {
            const MolGraph& target = targets[i];
            if (!fits(target)) return;
            LabelCount* const histogram = staged.data() + slots[i];
            int labels = 0, atom = 0;
            for (; atom < target.n; ++atom) {
                const int label = ubLabel(target, atom, chemistry);
                if (label < 0 || label >= QueryLabelCounts::denseLimit) break;
                if (scratch[label]++ == 0) histogram[labels++].label = static_cast<std::uint16_t>(label);
            }
            for (int k = 0; k < labels; ++k) {
                histogram[k].count = static_cast<std::uint16_t>(scratch[histogram[k].label]);
                scratch[histogram[k].label] = 0;
            }
            if (atom < target.n) return;
            atoms_[i] = target.n;
            distinct[i] = labels;
        });
        for (int i = 0; i < N; ++i) offsets_[i + 1] = offsets_[i] + distinct[i];
        counts_.reserve(offsets_[N]);
        for (int i = 0; i < N; ++i)
            counts_.insert(counts_.end(), staged.begin() + slots[i], staged.begin() + slots[i] + distinct[i]);
    }

    /** The label rule the histograms were built under, as labelMode numbers it. */
    int mode() const { return mode_; }
    std::size_t size() const { return atoms_.size(); }
    /** Whether the histogram of target i is held; otherwise the target is scored on its graph. */
    bool held(int i) const { return atoms_[i] >= 0; }
    int atomCount(int i) const { return atoms_[i]; }
    const LabelCount* begin(int i) const { return counts_.data() + offsets_[i]; }
    const LabelCount* end(int i) const { return counts_.data() + offsets_[i + 1]; }

private:
    static bool fits(const MolGraph& target) {
        return target.n >= 0 && target.n <= std::numeric_limits<std::uint16_t>::max();
    }

    int mode_;
    std::vector<int> atoms_;            // atom count, or -1 for a target scored on its graph
    std::vector<std::size_t> offsets_;  // target i holds counts_[offsets_[i], offsets_[i + 1])
    std::vector<LabelCount> counts_;
};

inline double QueryLabelCounts::similarity(const TargetLabels& labels, int target) const {
    int overlap = 0;
    // A label the query lacks has a count of zero there, so only shared labels add.
    for (const auto* entry = labels.begin(target); entry != labels.end(target); ++entry)
        overlap += std::min(counts_[entry->label], static_cast<int>(entry->count));
    return overlapSimilarity(overlap, query_.n, labels.atomCount(target));
}

/**
 * The similarity bound of every target to the query, in target order, on up to numThreads
 * threads as forEachTarget runs. Every score is computed as on one thread, so the scores do not
 * depend on the thread count. Targets needs size() and an operator[] whose result converts to
 * const MolGraph&; the Python binding reads SMILES strings there. When reading or scoring targets
 * throws, the exception of the first such target propagates, as on one thread.
 */
template<class Targets>
inline std::vector<double> similarityScores(const MolGraph& query, const Targets& targets,
                                            const ChemOptions& chemistry, int numThreads) {
    std::vector<double> scores(targets.size(), 0.0);
    const QueryLabelCounts queryCounts(query, chemistry);
    forEachTarget(static_cast<int>(scores.size()), numThreads, [&](int i, std::vector<int>& scratch) {
        scores[i] = queryCounts.similarity(targets[i], scratch);
    });
    return scores;
}

/**
 * The same scores, read from the label histograms of the targets where they are held. Targets left
 * to their graphs are scored there, as is every target for a query with a label outside the dense
 * table or for histograms built under another label rule or for another number of targets.
 */
inline std::vector<double> similarityScores(const MolGraph& query, const std::vector<MolGraph>& targets,
                                            const TargetLabels& labels, const ChemOptions& chemistry,
                                            int numThreads) {
    std::vector<double> scores(targets.size(), 0.0);
    const QueryLabelCounts queryCounts(query, chemistry);
    const bool histograms = queryCounts.dense() && labels.mode() == labelMode(chemistry)
                         && labels.size() == targets.size();
    forEachTarget(static_cast<int>(scores.size()), numThreads, [&](int i, std::vector<int>& scratch) {
        scores[i] = histograms && labels.held(i) ? queryCounts.similarity(labels, i)
                                                 : queryCounts.similarity(targets[i], scratch);
    });
    return scores;
}

/** The indices of the scores that reach the threshold, in order. */
inline std::vector<int> indicesReaching(const std::vector<double>& scores, double threshold) {
    std::vector<int> hits;
    for (int i = 0; i < static_cast<int>(scores.size()); ++i)
        if (scores[i] >= threshold) hits.push_back(i);
    return hits;
}

} // namespace screening_detail

inline double similarityUpperBound(const MolGraph& g1, const MolGraph& g2, const ChemOptions& chemistry) {
    return screening_detail::overlapSimilarity(screening_detail::labelFrequencyUpperBound(g1, g2, chemistry),
                                               g1.n, g2.n);
}

inline double similarityUpperBound(const MolGraph& g1, const MolGraph& g2) {
    return similarityUpperBound(g1, g2, ChemOptions());
}



/// Batch similarity-bound screening: the indices, in order, of the targets whose similarity bound
/// to the query reaches the threshold, computed on up to numThreads threads (0 = the OpenMP thread
/// limit). The hits do not depend on the thread count.
inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold, const ChemOptions& chemistry, int numThreads) {
    return screening_detail::indicesReaching(
        screening_detail::similarityScores(query, targets, chemistry, numThreads), threshold);
}

/// Batch similarity-bound screening on the calling thread.
inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold, const ChemOptions& chemistry) {
    return screenTargets(query, targets, threshold, chemistry, 1);
}

inline std::vector<int> screenTargets(const MolGraph& query,
                                       const std::vector<MolGraph>& targets,
                                       double threshold) {
    return screenTargets(query, targets, threshold, ChemOptions());
}

} // namespace smsd

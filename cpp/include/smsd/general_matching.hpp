/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Maximum-cardinality matching for undirected graphs with odd cycles.
 */
#pragma once
#ifndef SMSD_GENERAL_MATCHING_HPP
#define SMSD_GENERAL_MATCHING_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smsd { namespace detail {

// Edmonds' alternating-forest algorithm contracts odd cycles (blossoms) while
// searching for an augmenting path. The matching and parent links retain the
// original vertices, so augmentation also expands any contracted blossoms.
class GeneralMatchingSolver {
public:
    explicit GeneralMatchingSolver(const std::vector<std::vector<int>>& graph)
        : graph_(graph), n_(checkedSize(graph)), mate_(n_, -1), parent_(n_),
          base_(n_), queue_(n_), outer_(n_), blossom_(n_), ancestor_(n_) {
        for (const auto& neighbors : graph_) {
            for (int v : neighbors) {
                if (v < 0 || v >= n_)
                    throw std::out_of_range("maximumMatching: invalid vertex index");
            }
        }
    }

    std::vector<int> solve() {
        for (int root = 0; root < n_; ++root) {
            if (mate_[root] != -1) continue;
            int endpoint = augmentingPath(root);
            while (endpoint != -1) {
                const int previous = parent_[endpoint];
                const int next = previous == -1 ? -1 : mate_[previous];
                mate_[endpoint] = previous;
                if (previous != -1) mate_[previous] = endpoint;
                endpoint = next;
            }
        }
        return std::move(mate_);
    }

private:
    static int checkedSize(const std::vector<std::vector<int>>& graph) {
        if (graph.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("maximumMatching: too many vertices");
        return static_cast<int>(graph.size());
    }

    // Both vertices are outer vertices in the same alternating tree.
    int commonBase(int a, int b) {
        std::fill(ancestor_.begin(), ancestor_.end(), false);
        for (;;) {
            a = base_[a];
            ancestor_[a] = true;
            if (mate_[a] == -1) break;
            a = parent_[mate_[a]];
        }
        for (;;) {
            b = base_[b];
            if (ancestor_[b]) return b;
            b = parent_[mate_[b]];
        }
    }

    void markBranch(int vertex, int common, int child) {
        while (base_[vertex] != common) {
            blossom_[base_[vertex]] = true;
            blossom_[base_[mate_[vertex]]] = true;
            parent_[vertex] = child;
            child = mate_[vertex];
            vertex = parent_[mate_[vertex]];
        }
    }

    int augmentingPath(int root) {
        std::fill(outer_.begin(), outer_.end(), false);
        std::fill(parent_.begin(), parent_.end(), -1);
        std::iota(base_.begin(), base_.end(), 0);
        int head = 0, tail = 0;
        queue_[tail++] = root;
        outer_[root] = true;

        while (head < tail) {
            const int v = queue_[head++];
            for (int u : graph_[v]) {
                // Self loops and duplicate edges do not change the matching.
                if (base_[v] == base_[u] || mate_[v] == u) continue;
                if (u == root || (mate_[u] != -1 && parent_[mate_[u]] != -1)) {
                    const int common = commonBase(v, u);
                    std::fill(blossom_.begin(), blossom_.end(), false);
                    markBranch(v, common, u);
                    markBranch(u, common, v);
                    for (int i = 0; i < n_; ++i) {
                        if (!blossom_[base_[i]]) continue;
                        base_[i] = common;
                        if (!outer_[i]) {
                            outer_[i] = true;
                            queue_[tail++] = i;
                        }
                    }
                } else if (parent_[u] == -1) {
                    parent_[u] = v;
                    if (mate_[u] == -1) return u;
                    const int partner = mate_[u];
                    if (!outer_[partner]) {
                        outer_[partner] = true;
                        queue_[tail++] = partner;
                    }
                }
            }
        }
        return -1;
    }

    const std::vector<std::vector<int>>& graph_;
    int n_;
    std::vector<int> mate_, parent_, base_, queue_;
    std::vector<bool> outer_, blossom_, ancestor_;
};

/// Return a maximum-cardinality matching: mate[v] is v's partner or -1.
/// Input is an undirected adjacency list with symmetric edges and vertices
/// numbered [0, graph.size()). Self loops and repeated edges are allowed.
/// Uses O(V^3) time for a simple graph and O(V) auxiliary storage; handles
/// disconnected and non-bipartite graphs without external dependencies.
inline std::vector<int> maximumMatching(const std::vector<std::vector<int>>& graph) {
    return GeneralMatchingSolver(graph).solve();
}

}} // namespace smsd::detail

#endif // SMSD_GENERAL_MATCHING_HPP

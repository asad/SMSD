/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Lightweight maximum clique solver for MCS atom mapping.
 * Takes a precomputed modular product graph and finds maximum clique(s)
 * via Bron-Kerbosch with Tomita pivoting + k-core pruning.
 * McGregor bond-grow extension for completeness.
 *
 * Designed for reaction mapping: chemistry stays in Python/RDKit,
 * only the combinatorial search runs in C++. Zero MolGraph overhead.
 */
#pragma once

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include "smsd/bitops.hpp"

namespace smsd {
namespace clique {

// ---------------------------------------------------------------------------
// Input: modular product graph
// ---------------------------------------------------------------------------

/**
 * A vertex in the modular product graph represents a compatible atom pair
 * (query_atom, target_atom).  An edge between two vertices means the
 * corresponding bonds are also compatible.
 */
struct ProductVertex {
    int query_atom;   ///< 0-based atom index in query molecule
    int target_atom;  ///< 0-based atom index in target molecule
};

/**
 * Modular product graph for maximum clique = MCS.
 *
 * Build in Python from RDKit chemistry, pass to C++ for fast clique search.
 */
struct ProductGraph {
    std::vector<ProductVertex> vertices;
    std::vector<std::vector<int>> adj;  ///< adjacency list (vertex index → neighbors)
    int n = 0;                          ///< number of vertices

    void build(const std::vector<ProductVertex>& verts,
               const std::vector<std::pair<int,int>>& edges) {
        vertices = verts;
        n = static_cast<int>(vertices.size());
        adj.assign(n, {});
        for (auto& [u, v] : edges) {
            if (u >= 0 && u < n && v >= 0 && v < n) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
        for (auto& nbrs : adj) {
            std::sort(nbrs.begin(), nbrs.end());
        }
    }

    /**
     * Build product graph from atom compatibility + bond lookup tables.
     * Edge construction happens entirely in C++ — no Python loop.
     *
     * @param compat  Compatible atom pairs: (query_atom, target_atom), 0-based
     * @param bonds_a Bond lookup for molecule A: {(min_atom, max_atom): order}
     * @param bonds_b Bond lookup for molecule B: {(min_atom, max_atom): order}
     * @param bond_any  If true, any bond order matches (CompareAny mode)
     */
    void buildFromCompat(
        const std::vector<std::pair<int,int>>& compat,
        const std::map<std::pair<int,int>, int>& bonds_a,
        const std::map<std::pair<int,int>, int>& bonds_b,
        bool bond_any = true)
    {
        n = static_cast<int>(compat.size());
        vertices.resize(n);
        for (int i = 0; i < n; i++) {
            vertices[i].query_atom = compat[i].first;
            vertices[i].target_atom = compat[i].second;
        }
        adj.assign(n, {});

        for (int i = 0; i < n; i++) {
            int a1 = compat[i].first, b1 = compat[i].second;
            for (int j = i + 1; j < n; j++) {
                int a2 = compat[j].first, b2 = compat[j].second;
                if (a1 == a2 || b1 == b2) continue;

                auto ka = std::make_pair(std::min(a1,a2), std::max(a1,a2));
                auto kb = std::make_pair(std::min(b1,b2), std::max(b1,b2));
                auto it_a = bonds_a.find(ka);
                auto it_b = bonds_b.find(kb);
                bool has_a = (it_a != bonds_a.end());
                bool has_b = (it_b != bonds_b.end());

                if (has_a && has_b) {
                    // Both have bonds — check compatibility
                    if (bond_any || it_a->second == it_b->second) {
                        adj[i].push_back(j);
                        adj[j].push_back(i);
                    }
                }
                // If only one has a bond, NOT compatible (asymmetric bond)
                // If neither has a bond, no edge (non-adjacent pair)
            }
        }
        for (auto& nbrs : adj) {
            std::sort(nbrs.begin(), nbrs.end());
        }
    }
};

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

struct CliqueResult {
    /// All maximum cliques found (each as list of vertex indices)
    std::vector<std::vector<int>> cliques;
    /// Size of maximum clique
    int max_size = 0;
    /// Whether the search timed out
    bool timed_out = false;
    /// Microseconds spent
    int64_t elapsed_us = 0;
};

// ---------------------------------------------------------------------------
// Bron-Kerbosch with Tomita pivoting + k-core pruning
// ---------------------------------------------------------------------------

inline CliqueResult findMaxCliques(
    const ProductGraph& pg,
    int max_cliques = 8,
    int64_t timeout_ms = 1000,
    int incumbent_size = 0)
{
    CliqueResult result;
    if (pg.n == 0) return result;

    auto start = std::chrono::steady_clock::now();
    auto deadline = start + std::chrono::milliseconds(timeout_ms);

    // Build bitset adjacency for fast intersection
    // For graphs up to ~4000 vertices, use vector<uint64_t> bitsets
    int words = (pg.n + 63) / 64;
    std::vector<std::vector<uint64_t>> adjBits(pg.n, std::vector<uint64_t>(words, 0));
    std::vector<int> degree(pg.n, 0);

    for (int u = 0; u < pg.n; u++) {
        degree[u] = static_cast<int>(pg.adj[u].size());
        for (int v : pg.adj[u]) {
            adjBits[u][v / 64] |= (uint64_t(1) << (v % 64));
        }
    }

    // k-core pruning: iteratively remove vertices with degree < incumbent
    int best_size = incumbent_size;
    std::vector<bool> alive(pg.n, true);
    bool changed = true;
    while (changed) {
        changed = false;
        for (int u = 0; u < pg.n; u++) {
            if (alive[u] && degree[u] < best_size) {
                alive[u] = false;
                for (int v : pg.adj[u]) {
                    if (alive[v]) degree[v]--;
                }
                changed = true;
            }
        }
    }

    std::vector<int> active;
    active.reserve(pg.n);
    for (int u = 0; u < pg.n; u++) {
        if (alive[u]) active.push_back(u);
    }
    if (static_cast<int>(active.size()) < best_size) {
        return result;
    }

    // Rebuild bitsets for alive vertices only
    auto setBit = [&](std::vector<uint64_t>& bits, int v) {
        bits[v / 64] |= (uint64_t(1) << (v % 64));
    };
    auto testBit = [&](const std::vector<uint64_t>& bits, int v) -> bool {
        return (bits[v / 64] >> (v % 64)) & 1;
    };
    auto intersectCount = [&](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) -> int {
        int count = 0;
        for (int w = 0; w < words; w++) {
            count += smsd::popcount64(a[w] & b[w]);
        }
        return count;
    };
    auto intersectInto = [&](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b,
                             std::vector<uint64_t>& out) {
        for (int w = 0; w < words; w++) out[w] = a[w] & b[w];
    };
    auto isZero = [&](const std::vector<uint64_t>& bits) -> bool {
        for (int w = 0; w < words; w++) if (bits[w]) return false;
        return true;
    };
    auto bitList = [&](const std::vector<uint64_t>& bits) -> std::vector<int> {
        std::vector<int> out;
        for (int w = 0; w < words; w++) {
            uint64_t word = bits[w];
            while (word) {
                int bit = smsd::ctz64(word);
                out.push_back(w * 64 + bit);
                word &= word - 1;
            }
        }
        return out;
    };

    // BK with Tomita pivoting
    bool timeout_flag = false;
    int call_count = 0;

    std::function<void(std::vector<int>&, std::vector<uint64_t>&, std::vector<uint64_t>&)>
    bk = [&](std::vector<int>& R, std::vector<uint64_t>& P, std::vector<uint64_t>& X) {
        if (timeout_flag) return;
        if (static_cast<int>(result.cliques.size()) >= max_cliques) return;

        if (++call_count % 1024 == 0) {
            if (std::chrono::steady_clock::now() >= deadline) {
                timeout_flag = true;
                return;
            }
        }

        if (isZero(P) && isZero(X)) {
            int sz = static_cast<int>(R.size());
            if (sz > best_size) {
                best_size = sz;
                result.cliques.clear();
                result.cliques.push_back(R);
            } else if (sz == best_size && sz > 0) {
                result.cliques.push_back(R);
            }
            return;
        }

        // Pivot: vertex in P ∪ X with most neighbors in P
        std::vector<uint64_t> PuX(words);
        for (int w = 0; w < words; w++) PuX[w] = P[w] | X[w];
        auto puxList = bitList(PuX);
        if (puxList.empty()) return;

        int pivot = puxList[0];
        int pivotCount = intersectCount(adjBits[pivot], P);
        for (int v : puxList) {
            int cnt = intersectCount(adjBits[v], P);
            if (cnt > pivotCount) { pivot = v; pivotCount = cnt; }
        }

        // Candidates: P \ N(pivot)
        std::vector<uint64_t> cands(words);
        for (int w = 0; w < words; w++) cands[w] = P[w] & ~adjBits[pivot][w];
        auto candList = bitList(cands);

        for (int v : candList) {
            if (timeout_flag || static_cast<int>(result.cliques.size()) >= max_cliques) return;

            R.push_back(v);
            std::vector<uint64_t> newP(words), newX(words);
            intersectInto(P, adjBits[v], newP);
            intersectInto(X, adjBits[v], newX);
            bk(R, newP, newX);
            R.pop_back();

            // P := P \ {v}; X := X ∪ {v}
            P[v / 64] &= ~(uint64_t(1) << (v % 64));
            X[v / 64] |= (uint64_t(1) << (v % 64));
        }
    };

    // Initialize P = all active, X = empty
    std::vector<uint64_t> P(words, 0), X(words, 0);
    for (int u : active) setBit(P, u);
    std::vector<int> R;
    R.reserve(pg.n);

    bk(R, P, X);

    auto end = std::chrono::steady_clock::now();
    result.max_size = best_size;
    result.timed_out = timeout_flag;
    result.elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return result;
}

// ---------------------------------------------------------------------------
// Full MCS pipeline: greedy → seed-extend → BK → McGregor DFS
// ---------------------------------------------------------------------------

struct MCSResult {
    /// All candidate mappings: (query_atom, target_atom) pairs, 0-based
    std::vector<std::vector<std::pair<int,int>>> candidates;
    /// Best mapping size
    int best_size = 0;
    /// LFUB (label-frequency upper bound)
    int lfub = 0;
    /// Microseconds spent
    int64_t elapsed_us = 0;
};

/**
 * Complete MCS pipeline for reaction mapping.
 *
 * @param compat    Compatible atom pairs (query, target), 0-based
 * @param bonds_a   Bond lookup for mol A: {(min,max) -> order}
 * @param bonds_b   Bond lookup for mol B: {(min,max) -> order}
 * @param n_a       Number of atoms in mol A
 * @param n_b       Number of atoms in mol B
 * @param bond_any  If true, any bond order compatible
 * @param timeout_ms Total time budget
 * @param max_results Maximum candidates to return
 */
// ---------------------------------------------------------------------------
// McGregor DFS backtracking extension
// ---------------------------------------------------------------------------

/**
 * Extend a seed mapping by DFS backtracking along bond-compatible paths.
 *
 * Unlike greedy BFS (which takes the first legal extension), this tries
 * ALL candidates at each frontier atom and backtracks when stuck, finding
 * the optimal (largest) extension. Explores the full search tree within
 * the time budget.
 */
inline std::vector<std::pair<int,int>> mcgregorDFSExtend(
    const std::vector<std::pair<int,int>>& seed,
    const std::vector<std::pair<int,int>>& compat,
    const std::vector<std::vector<int>>& compat_by_query,
    const std::map<std::pair<int,int>, int>& bonds_a,
    const std::map<std::pair<int,int>, int>& bonds_b,
    const std::vector<std::vector<int>>& adj_a,
    int n_a, int n_b, bool bond_any,
    std::chrono::steady_clock::time_point deadline)
{
    if (seed.empty()) return {};

    // Initialize from seed
    std::map<int,int> current;
    std::set<int> used_q, used_t;
    for (auto& [q,t] : seed) {
        current[q] = t;
        used_q.insert(q);
        used_t.insert(t);
    }

    std::vector<std::pair<int,int>> best;
    for (auto& [q,t] : current) best.emplace_back(q,t);

    int call_count = 0;
    bool timed_out = false;

    // Get initial frontier: unmapped query atoms adjacent to mapped ones
    auto getFrontier = [&]() -> std::vector<int> {
        std::set<int> frontier_set;
        for (auto& [q,t] : current) {
            for (int nb : adj_a[q]) {
                if (!used_q.count(nb)) frontier_set.insert(nb);
            }
        }
        // Sort by most-constrained first (fewest compatible targets)
        std::vector<int> frontier(frontier_set.begin(), frontier_set.end());
        std::sort(frontier.begin(), frontier.end(), [&](int a, int b) {
            return compat_by_query[a].size() < compat_by_query[b].size();
        });
        return frontier;
    };

    // Find feasible targets for a query atom
    auto feasibleTargets = [&](int qa) -> std::vector<int> {
        std::vector<int> targets;
        for (int ci : compat_by_query[qa]) {
            int ta = compat[ci].second;
            if (used_t.count(ta)) continue;

            // Check bond compatibility with ALL mapped neighbors
            bool ok = true;
            bool has_mapped_nb = false;
            for (int nb_q : adj_a[qa]) {
                auto it = current.find(nb_q);
                if (it == current.end()) continue;
                has_mapped_nb = true;
                int nb_t = it->second;
                auto ka = std::make_pair(std::min(qa,nb_q), std::max(qa,nb_q));
                auto kb = std::make_pair(std::min(ta,nb_t), std::max(ta,nb_t));
                auto ia = bonds_a.find(ka);
                auto ib = bonds_b.find(kb);
                if (ia != bonds_a.end()) {
                    if (ib == bonds_b.end()) { ok = false; break; }
                    if (!bond_any && ia->second != ib->second) { ok = false; break; }
                }
            }
            if (ok && has_mapped_nb) targets.push_back(ta);
        }
        return targets;
    };

    // DFS backtracking over frontier atoms
    // For each frontier atom: try all feasible targets OR skip it.
    // Recurse with the updated frontier (new neighbors + remaining).
    std::set<int> skipped;

    std::function<void(std::vector<int>)> dfs = [&](std::vector<int> frontier) {
        if (timed_out) return;
        if (++call_count % 2048 == 0) {
            if (std::chrono::steady_clock::now() >= deadline) {
                timed_out = true;
                return;
            }
        }

        // Save current as best if larger
        if (static_cast<int>(current.size()) > static_cast<int>(best.size())) {
            best.clear();
            for (auto& [q,t] : current) best.emplace_back(q,t);
        }

        if (frontier.empty()) return;

        // Pick most-constrained atom from frontier
        int pick_idx = 0;
        int min_targets = INT_MAX;
        for (int i = 0; i < static_cast<int>(frontier.size()); i++) {
            int qa = frontier[i];
            if (used_q.count(qa) || skipped.count(qa)) continue;
            int count = 0;
            for (int ci : compat_by_query[qa]) {
                if (!used_t.count(compat[ci].second)) count++;
            }
            if (count < min_targets) {
                min_targets = count;
                pick_idx = i;
            }
        }

        int pick = frontier[pick_idx];
        if (used_q.count(pick) || skipped.count(pick)) return;

        auto targets = feasibleTargets(pick);

        // Remaining frontier (without pick)
        std::vector<int> remaining;
        for (int i = 0; i < static_cast<int>(frontier.size()); i++) {
            if (i != pick_idx && !used_q.count(frontier[i]) && !skipped.count(frontier[i]))
                remaining.push_back(frontier[i]);
        }

        // Branch 1: assign pick to each feasible target
        for (int ta : targets) {
            if (timed_out) return;
            current[pick] = ta;
            used_q.insert(pick);
            used_t.insert(ta);

            // Add new frontier atoms (neighbors of pick not yet mapped/skipped)
            auto next_frontier = remaining;
            for (int nb : adj_a[pick]) {
                if (!used_q.count(nb) && !skipped.count(nb)) {
                    bool already = false;
                    for (int f : next_frontier) if (f == nb) { already = true; break; }
                    if (!already) next_frontier.push_back(nb);
                }
            }

            dfs(next_frontier);

            current.erase(pick);
            used_q.erase(pick);
            used_t.erase(ta);
        }

        // Branch 2: skip pick entirely (it might block better extensions)
        skipped.insert(pick);
        dfs(remaining);
        skipped.erase(pick);
    };

    dfs(getFrontier());
    return best;
}

// Forward declaration
inline std::vector<std::vector<std::pair<int,int>>> vf2ReEmbed(
    const std::vector<std::pair<int,int>>& seed_mapping,
    const std::vector<std::pair<int,int>>& compat,
    const std::vector<std::vector<int>>& compat_by_query,
    const std::map<std::pair<int,int>, int>& bonds_a,
    const std::map<std::pair<int,int>, int>& bonds_b,
    const std::vector<std::vector<int>>& adj_a,
    const std::vector<std::vector<int>>& adj_b,
    int n_a, int n_b, bool bond_any,
    std::chrono::steady_clock::time_point deadline, int max_results);

inline MCSResult findMCSPipeline(
    const std::vector<std::pair<int,int>>& compat,
    const std::map<std::pair<int,int>, int>& bonds_a,
    const std::map<std::pair<int,int>, int>& bonds_b,
    int n_a, int n_b,
    bool bond_any = true,
    int64_t timeout_ms = 1000,
    int max_results = 8,
    const std::vector<bool>& ring_a = {},
    const std::vector<bool>& ring_b = {},
    const std::vector<bool>& arom_a = {},
    const std::vector<bool>& arom_b = {},
    int lfub_value = -1)
{
    auto start = std::chrono::steady_clock::now();
    auto deadline = start + std::chrono::milliseconds(timeout_ms);
    MCSResult result;

    if (compat.empty()) return result;
    int n = static_cast<int>(compat.size());

    // LFUB: true label-frequency upper bound passed from Python
    result.lfub = (lfub_value > 0) ? lfub_value : std::min(n_a, n_b);

    // Build adjacency index for fast neighbor lookup
    std::vector<std::vector<int>> adj_a(n_a), adj_b(n_b);
    for (auto& [pair, order] : bonds_a) {
        adj_a[pair.first].push_back(pair.second);
        adj_a[pair.second].push_back(pair.first);
    }
    for (auto& [pair, order] : bonds_b) {
        adj_b[pair.first].push_back(pair.second);
        adj_b[pair.second].push_back(pair.first);
    }

    // Index: for each query atom, which compat entries involve it
    std::vector<std::vector<int>> compat_by_query(n_a);
    for (int i = 0; i < n; i++) {
        compat_by_query[compat[i].first].push_back(i);
    }

    // --- L0.75: Greedy probe with multiple orderings for diversity ---
    auto runGreedy = [&](const std::vector<int>& order) -> std::vector<std::pair<int,int>> {
        std::set<int> used_q, used_t;
        std::map<int,int> q_to_t;
        std::vector<std::pair<int,int>> mapping;

        for (int qa : order) {
            if (used_q.count(qa)) continue;
            int best_t = -1;
            double best_score = -1.0;
            for (int ci : compat_by_query[qa]) {
                int ta = compat[ci].second;
                if (used_t.count(ta)) continue;
                double score = 0.0;
                for (int nb_q : adj_a[qa]) {
                    auto it = q_to_t.find(nb_q);
                    if (it == q_to_t.end()) continue;
                    int nb_t = it->second;
                    auto ka = std::make_pair(std::min(qa,nb_q), std::max(qa,nb_q));
                    auto kb = std::make_pair(std::min(ta,nb_t), std::max(ta,nb_t));
                    auto ia = bonds_a.find(ka);
                    auto ib = bonds_b.find(kb);
                    if (ia != bonds_a.end() && ib != bonds_b.end()) {
                        score += 4.0;
                        if (bond_any || ia->second == ib->second) score += 2.0;
                    }
                }
                if (adj_a[qa].size() == adj_b[ta].size()) score += 1.0;
                if (!ring_a.empty() && !ring_b.empty() && ring_a[qa] == ring_b[ta]) score += 1.5;
                if (!arom_a.empty() && !arom_b.empty() && arom_a[qa] == arom_b[ta]) score += 1.5;
                if (score > best_score) { best_score = score; best_t = ta; }
            }
            if (best_t >= 0) {
                mapping.emplace_back(qa, best_t);
                used_q.insert(qa); used_t.insert(best_t); q_to_t[qa] = best_t;
            }
        }
        // Largest connected component
        if (mapping.size() > 1) {
            std::set<int> mq;
            for (auto& [q,t] : mapping) mq.insert(q);
            std::set<int> vis;
            std::vector<int> bfs = {mapping[0].first};
            vis.insert(mapping[0].first);
            while (!bfs.empty()) {
                int cur = bfs.back(); bfs.pop_back();
                for (int nb : adj_a[cur]) {
                    if (mq.count(nb) && !vis.count(nb)) { vis.insert(nb); bfs.push_back(nb); }
                }
            }
            std::vector<std::pair<int,int>> conn;
            for (auto& [q,t] : mapping) if (vis.count(q)) conn.emplace_back(q,t);
            return conn;
        }
        return mapping;
    };

    std::vector<int> query_order(n_a);
    std::iota(query_order.begin(), query_order.end(), 0);
    std::sort(query_order.begin(), query_order.end(), [&](int a, int b) {
        return compat_by_query[a].size() < compat_by_query[b].size();
    });

    auto greedy = runGreedy(query_order);
    if (!greedy.empty()) {
        result.best_size = static_cast<int>(greedy.size());
        result.candidates.push_back(greedy);
    }

    // --- L1.5: Seed-and-extend from best greedy ---
    static const std::vector<std::pair<int,int>> empty_mapping;
    auto& best_greedy = result.candidates.empty()
        ? empty_mapping
        : *std::max_element(result.candidates.begin(), result.candidates.end(),
            [](const auto& a, const auto& b) { return a.size() < b.size(); });
    if (!result.candidates.empty()) {
        std::map<int,int> q_to_t;
        std::set<int> used_q, used_t;
        for (auto& [q,t] : best_greedy) { q_to_t[q]=t; used_q.insert(q); used_t.insert(t); }

        std::vector<int> frontier;
        for (auto& [q,t] : best_greedy) {
            for (int nb : adj_a[q]) {
                if (!used_q.count(nb)) frontier.push_back(nb);
            }
        }

        bool extended = true;
        while (extended && std::chrono::steady_clock::now() < deadline) {
            extended = false;
            std::vector<int> new_frontier;
            for (int qa : frontier) {
                if (used_q.count(qa)) continue;
                int best_t = -1;
                double best_score = -1.0;
                for (int ci : compat_by_query[qa]) {
                    int ta = compat[ci].second;
                    if (used_t.count(ta)) continue;
                    double score = 0.0;
                    for (int nb_q : adj_a[qa]) {
                        auto it = q_to_t.find(nb_q);
                        if (it == q_to_t.end()) continue;
                        int nb_t = it->second;
                        auto ka = std::make_pair(std::min(qa,nb_q), std::max(qa,nb_q));
                        auto kb = std::make_pair(std::min(ta,nb_t), std::max(ta,nb_t));
                        auto ia = bonds_a.find(ka);
                        auto ib = bonds_b.find(kb);
                        if (ia != bonds_a.end() && ib != bonds_b.end()) {
                            score += 3.0;
                            if (bond_any || ia->second == ib->second)
                                score += 1.0;
                        }
                    }
                    if (score > best_score) { best_score = score; best_t = ta; }
                }
                if (best_t >= 0 && best_score > 0) {
                    q_to_t[qa] = best_t;
                    used_q.insert(qa);
                    used_t.insert(best_t);
                    extended = true;
                    for (int nb : adj_a[qa]) {
                        if (!used_q.count(nb)) new_frontier.push_back(nb);
                    }
                }
            }
            frontier = new_frontier;
        }

        std::vector<std::pair<int,int>> extended_mapping;
        for (auto& [q,t] : q_to_t) extended_mapping.emplace_back(q,t);
        if (static_cast<int>(extended_mapping.size()) > result.best_size) {
            result.best_size = static_cast<int>(extended_mapping.size());
            result.candidates.push_back(extended_mapping);
        }
    } // if (!result.candidates.empty())

    // --- VF2 re-embedding: diverse embeddings ---
    // Tight per-pair deadline (100ms) for VF2 — enough for diversity,
    // fast enough to not dominate total time.
    if (!result.candidates.empty() && std::chrono::steady_clock::now() < deadline) {
        auto& best_map = *std::max_element(result.candidates.begin(), result.candidates.end(),
            [](const auto& a, const auto& b) { return a.size() < b.size(); });
        auto vf2_deadline = std::min(deadline,
            std::chrono::steady_clock::now() + std::chrono::milliseconds(100));

        {
            auto embeddings = vf2ReEmbed(
                best_map, compat, compat_by_query, bonds_a, bonds_b, adj_a, adj_b,
                n_a, n_b, bond_any, vf2_deadline, std::min(max_results, 4));

            // McGregor DFS backtracking: extend each VF2 embedding
            // by exploring ALL possible extensions with backtracking.
            // Unlike greedy BFS, this tries all candidates at each step
            // and backtracks when a choice blocks a better downstream path.
            for (auto& emb : embeddings) {
                if (std::chrono::steady_clock::now() >= deadline) break;

                auto extended = mcgregorDFSExtend(
                    emb, compat, compat_by_query, bonds_a, bonds_b,
                    adj_a, n_a, n_b, bond_any, vf2_deadline);

                if (static_cast<int>(extended.size()) > result.best_size)
                    result.best_size = static_cast<int>(extended.size());
                result.candidates.push_back(extended);
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    result.elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    return result;
}

// ---------------------------------------------------------------------------
// Lightweight VF2 re-embedding (no MolGraph dependency)
// ---------------------------------------------------------------------------

/**
 * Given a mapping (query→target atom pairs), find ALL alternative
 * embeddings of the same query subgraph in the target molecule.
 *
 * Uses VF2-style backtracking with bond-compatibility feasibility.
 * Works entirely on compat tables and bond lookups — zero MolGraph.
 */
inline std::vector<std::vector<std::pair<int,int>>> vf2ReEmbed(
    const std::vector<std::pair<int,int>>& seed_mapping,
    const std::vector<std::pair<int,int>>& compat,
    const std::vector<std::vector<int>>& compat_by_query,
    const std::map<std::pair<int,int>, int>& bonds_a,
    const std::map<std::pair<int,int>, int>& bonds_b,
    const std::vector<std::vector<int>>& adj_a,
    const std::vector<std::vector<int>>& adj_b,
    int n_a, int n_b,
    bool bond_any,
    std::chrono::steady_clock::time_point deadline,
    int max_results = 8)
{
    std::vector<std::vector<std::pair<int,int>>> results;
    if (seed_mapping.empty()) return results;

    // Extract query atom ordering (same atoms as seed, sorted by connectivity)
    std::vector<int> query_atoms;
    query_atoms.reserve(seed_mapping.size());
    for (auto& [q, t] : seed_mapping) query_atoms.push_back(q);

    // Sort query atoms: most-constrained first (fewest compatible targets)
    std::sort(query_atoms.begin(), query_atoms.end(), [&](int a, int b) {
        return compat_by_query[a].size() < compat_by_query[b].size();
    });

    int depth = static_cast<int>(query_atoms.size());

    // Build query subgraph bonds (which query-atom pairs are bonded)
    std::set<int> query_set(query_atoms.begin(), query_atoms.end());
    std::vector<std::vector<int>> query_adj(n_a);
    for (int qa : query_atoms) {
        for (int nb : adj_a[qa]) {
            if (query_set.count(nb)) query_adj[qa].push_back(nb);
        }
    }

    // Deduplicate results
    std::set<std::vector<std::pair<int,int>>> seen;
    // Add seed mapping
    {
        auto sorted_seed = seed_mapping;
        std::sort(sorted_seed.begin(), sorted_seed.end());
        seen.insert(sorted_seed);
    }

    // VF2 backtracking
    std::map<int,int> current;  // query → target
    std::set<int> used_targets;
    int call_count = 0;

    std::function<void(int)> backtrack = [&](int level) {
        if (std::chrono::steady_clock::now() >= deadline) return;
        if (static_cast<int>(results.size()) >= max_results) return;

        if (level == depth) {
            // Complete mapping found
            std::vector<std::pair<int,int>> mapping;
            for (auto& [q,t] : current) mapping.emplace_back(q, t);
            std::sort(mapping.begin(), mapping.end());
            if (!seen.count(mapping)) {
                seen.insert(mapping);
                results.push_back(mapping);
            }
            return;
        }

        int qa = query_atoms[level];

        // Find feasible targets for this query atom
        for (int ci : compat_by_query[qa]) {
            int ta = compat[ci].second;
            if (used_targets.count(ta)) continue;

            // Feasibility: check bonds to already-mapped query neighbors
            bool feasible = true;
            for (int nb_q : query_adj[qa]) {
                auto it = current.find(nb_q);
                if (it == current.end()) continue;  // not yet mapped — skip
                int nb_t = it->second;

                // Bond (qa, nb_q) must exist in A
                auto ka = std::make_pair(std::min(qa, nb_q), std::max(qa, nb_q));
                auto ia = bonds_a.find(ka);
                if (ia == bonds_a.end()) continue;  // no bond in query subgraph

                // Corresponding bond (ta, nb_t) must exist in B
                auto kb = std::make_pair(std::min(ta, nb_t), std::max(ta, nb_t));
                auto ib = bonds_b.find(kb);
                if (ib == bonds_b.end()) { feasible = false; break; }
                if (!bond_any && ia->second != ib->second) { feasible = false; break; }
            }
            if (!feasible) continue;

            // Extend
            current[qa] = ta;
            used_targets.insert(ta);
            backtrack(level + 1);
            current.erase(qa);
            used_targets.erase(ta);

            if (++call_count % 4096 == 0) {
                if (std::chrono::steady_clock::now() >= deadline) return;
            }
            if (static_cast<int>(results.size()) >= max_results) return;
        }
    };

    backtrack(0);
    return results;
}

// ---------------------------------------------------------------------------
// Substructure match: find all embeddings of query in target
// ---------------------------------------------------------------------------

/**
 * VF2-based substructure matching.
 *
 * Finds all distinct embeddings of query (n_query atoms) inside target
 * (n_target atoms). Wraps vf2ReEmbed with automatic adjacency/index
 * construction from bond maps and compat pairs.
 *
 * @param n_query     Number of atoms in query molecule
 * @param n_target    Number of atoms in target molecule
 * @param compat      Compatible atom pairs (query, target), 0-based
 * @param bonds_query Bond lookup for query: {(min,max) -> order}, 0-based
 * @param bonds_target Bond lookup for target: {(min,max) -> order}, 0-based
 * @param bond_any    If true, any bond order matches
 * @param timeout_ms  Time budget in milliseconds
 * @param max_results Maximum embeddings to return
 * @return List of embeddings, each a vector of (query_atom, target_atom) pairs
 */
inline std::vector<std::vector<std::pair<int,int>>> substructureMatch(
    int n_query, int n_target,
    const std::vector<std::pair<int,int>>& compat,
    const std::map<std::pair<int,int>, int>& bonds_query,
    const std::map<std::pair<int,int>, int>& bonds_target,
    bool bond_any = false,
    int64_t timeout_ms = 500,
    int max_results = 8)
{
    if (n_query <= 0 || n_target <= 0 || compat.empty()) return {};

    // Build compat_by_query index
    std::vector<std::vector<int>> compat_by_query(static_cast<size_t>(n_query));
    for (size_t i = 0; i < compat.size(); i++) {
        int q = compat[i].first;
        if (q >= 0 && q < n_query) {
            compat_by_query[static_cast<size_t>(q)].push_back(static_cast<int>(i));
        }
    }

    // Build adjacency lists from bond maps
    std::vector<std::vector<int>> adj_query(static_cast<size_t>(n_query));
    for (const auto& [bond, _order] : bonds_query) {
        int a = bond.first, b = bond.second;
        if (a >= 0 && a < n_query && b >= 0 && b < n_query) {
            adj_query[static_cast<size_t>(a)].push_back(b);
            adj_query[static_cast<size_t>(b)].push_back(a);
        }
    }
    std::vector<std::vector<int>> adj_target(static_cast<size_t>(n_target));
    for (const auto& [bond, _order] : bonds_target) {
        int a = bond.first, b = bond.second;
        if (a >= 0 && a < n_target && b >= 0 && b < n_target) {
            adj_target[static_cast<size_t>(a)].push_back(b);
            adj_target[static_cast<size_t>(b)].push_back(a);
        }
    }

    auto deadline = std::chrono::steady_clock::now()
        + std::chrono::milliseconds(timeout_ms);

    // Empty seed: VF2 starts from scratch
    std::vector<std::pair<int,int>> empty_seed;

    return vf2ReEmbed(
        empty_seed, compat, compat_by_query,
        bonds_query, bonds_target,
        adj_query, adj_target,
        n_query, n_target, bond_any,
        deadline, max_results);
}

/**
 * Element-based substructure match: builds compat table from element lists.
 *
 * Avoids Python O(Q×T) loop by constructing compat pairs in C++.
 * Compatible atoms = same element string.
 */
inline std::vector<std::vector<std::pair<int,int>>> substructureMatchFromElements(
    const std::vector<std::string>& elements_query,
    const std::vector<std::string>& elements_target,
    const std::map<std::pair<int,int>, int>& bonds_query,
    const std::map<std::pair<int,int>, int>& bonds_target,
    bool bond_any = false,
    int64_t timeout_ms = 500,
    int max_results = 8)
{
    const int nQ = static_cast<int>(elements_query.size());
    const int nT = static_cast<int>(elements_target.size());
    if (nQ <= 0 || nT <= 0) return {};

    // Build compat from element matching — O(Q×T) but in C++
    std::vector<std::pair<int,int>> compat;
    compat.reserve(static_cast<size_t>(nQ) * 4);  // heuristic reservation
    for (int q = 0; q < nQ; q++) {
        for (int t = 0; t < nT; t++) {
            if (elements_query[static_cast<size_t>(q)] == elements_target[static_cast<size_t>(t)]) {
                compat.emplace_back(q, t);
            }
        }
    }
    if (compat.empty()) return {};

    return substructureMatch(nQ, nT, compat, bonds_query, bonds_target,
                             bond_any, timeout_ms, max_results);
}

/**
 * Score a pair mapping by bond preservation and atom property matching.
 *
 * @param mapping       Atom pairs (query, target), 0-based
 * @param elements_a    Element symbols for molecule A
 * @param elements_b    Element symbols for molecule B
 * @param bonds_a       Bond lookup {(min,max) -> order} for A, 0-based
 * @param bonds_b       Bond lookup {(min,max) -> order} for B, 0-based
 * @param ring_a        Per-atom ring flags for A
 * @param ring_b        Per-atom ring flags for B
 * @param arom_a        Per-atom aromaticity flags for A
 * @param arom_b        Per-atom aromaticity flags for B
 * @return Composite score (higher = better mapping)
 */
inline double scorePairMapping(
    const std::vector<std::pair<int,int>>& mapping,
    const std::vector<std::string>& elements_a,
    const std::vector<std::string>& elements_b,
    const std::map<std::pair<int,int>, int>& bonds_a,
    const std::map<std::pair<int,int>, int>& bonds_b,
    const std::vector<bool>& ring_a = {},
    const std::vector<bool>& ring_b = {},
    const std::vector<bool>& arom_a = {},
    const std::vector<bool>& arom_b = {})
{
    if (mapping.empty()) return 0.0;

    double score = 0.0;
    const int nA = static_cast<int>(elements_a.size());
    const int nB = static_cast<int>(elements_b.size());

    // Build reverse lookup: target → query
    std::vector<int> targetToQuery(static_cast<size_t>(nB), -1);
    for (const auto& [q, t] : mapping) {
        if (t >= 0 && t < nB) {
            targetToQuery[static_cast<size_t>(t)] = q;
        }
    }

    // Atom-level scoring
    for (const auto& [q, t] : mapping) {
        if (q < 0 || q >= nA || t < 0 || t >= nB) continue;
        // Element match (should always be true for valid mappings)
        if (elements_a[static_cast<size_t>(q)] == elements_b[static_cast<size_t>(t)]) {
            score += 1.0;
        }
        // Ring match
        if (!ring_a.empty() && !ring_b.empty()
            && static_cast<size_t>(q) < ring_a.size()
            && static_cast<size_t>(t) < ring_b.size()
            && ring_a[static_cast<size_t>(q)] == ring_b[static_cast<size_t>(t)]) {
            score += 1.5;
        }
        // Aromatic match
        if (!arom_a.empty() && !arom_b.empty()
            && static_cast<size_t>(q) < arom_a.size()
            && static_cast<size_t>(t) < arom_b.size()
            && arom_a[static_cast<size_t>(q)] == arom_b[static_cast<size_t>(t)]) {
            score += 1.5;
        }
    }

    // Build query→target lookup for O(1) bond partner resolution
    std::vector<int> queryToTarget(static_cast<size_t>(nA), -1);
    for (const auto& [q, t] : mapping) {
        if (q >= 0 && q < nA) queryToTarget[static_cast<size_t>(q)] = t;
    }

    // Bond preservation scoring
    int bondsPreserved = 0, bondsTotal = 0;
    for (const auto& [bondKey, orderA] : bonds_a) {
        int a1 = bondKey.first, a2 = bondKey.second;
        if (a1 < 0 || a1 >= nA || a2 < 0 || a2 >= nA) continue;
        int t1 = queryToTarget[static_cast<size_t>(a1)];
        int t2 = queryToTarget[static_cast<size_t>(a2)];
        if (t1 < 0 || t2 < 0) continue;
        ++bondsTotal;
        auto key = std::make_pair(std::min(t1, t2), std::max(t1, t2));
        auto it = bonds_b.find(key);
        if (it != bonds_b.end()) {
            ++bondsPreserved;
            if (it->second == orderA) {
                score += 4.0;  // exact order match
            } else {
                score += 2.0;  // bond exists, order differs
            }
        }
    }

    // Coverage bonus
    score += static_cast<double>(mapping.size()) * 0.5;

    return score;
}

} // namespace clique
} // namespace smsd

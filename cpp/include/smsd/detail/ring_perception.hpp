/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph ring perception: ring systems from bridges, the smallest set of
 * smallest rings, relevant cycles, unique ring families, per-atom ring counts
 * and ring-system signatures.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// Ring counts per atom
// ============================================================================

inline void MolGraph::computeRingCounts(const std::vector<std::vector<int>>& rings) const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    std::fill(ringCount.begin(), ringCount.end(), 0);
    for (const auto& ring : rings) {
        for (int atom : ring) ringCount[atom]++;
    }
}

// The default computeRingCounts() counts the relevant cycles of the ring candidate search, as
// the Java implementation does; it is defined in ring_finder.hpp, which this header includes at its end.

// ============================================================================
// Ring systems from bridges
// ============================================================================

inline MolGraph::CycleSystems MolGraph::cycleSystems() const {
    CycleSystems result;
    std::vector<int> listedBy(n, -1);
    for (int atom = 0; atom < n; ++atom)
        for (int neighbor : neighbors[atom]) {
            if (neighbor == atom || listedBy[neighbor] == atom) return result;
            listedBy[neighbor] = atom;
        }
    // Iterative depth-first search: a tree bond is a bridge when no back bond from the
    // child's subtree reaches its parent or above.
    std::vector<int> order(n, -1), low(n, 0), treeParent(n, -1);
    std::vector<uint8_t> bridgeToParent(n, 0);
    struct Frame { int atom; size_t next; };
    std::vector<Frame> stack;
    int counter = 0;
    for (int root = 0; root < n; ++root) {
        if (order[root] >= 0) continue;
        order[root] = low[root] = counter++;
        stack.push_back({root, 0});
        while (!stack.empty()) {
            Frame& frame = stack.back();
            const int atom = frame.atom;
            if (frame.next < neighbors[atom].size()) {
                const int neighbor = neighbors[atom][frame.next++];
                if (neighbor == treeParent[atom]) continue;
                if (order[neighbor] >= 0) {
                    low[atom] = std::min(low[atom], order[neighbor]);
                    continue;
                }
                treeParent[neighbor] = atom;
                order[neighbor] = low[neighbor] = counter++;
                stack.push_back({neighbor, 0});
                continue;
            }
            stack.pop_back();
            const int parent = treeParent[atom];
            if (parent >= 0) {
                low[parent] = std::min(low[parent], low[atom]);
                bridgeToParent[atom] = low[atom] > order[parent];
            }
        }
    }
    auto bridge = [&](int first, int second) {
        return (treeParent[second] == first && bridgeToParent[second])
            || (treeParent[first] == second && bridgeToParent[first]);
    };
    result.ofAtom.assign(n, -1);
    std::vector<int> queue;
    for (int start = 0; start < n; ++start) {
        if (result.ofAtom[start] >= 0) continue;
        bool cyclic = false;
        for (int neighbor : neighbors[start]) if (!bridge(start, neighbor)) { cyclic = true; break; }
        if (!cyclic) continue;
        const int system = static_cast<int>(result.atoms.size());
        size_t atoms = 0, bondEnds = 0;
        queue.assign(1, start);
        result.ofAtom[start] = system;
        for (size_t head = 0; head < queue.size(); ++head) {
            const int atom = queue[head];
            ++atoms;
            for (int neighbor : neighbors[atom]) {
                if (bridge(atom, neighbor)) continue;
                ++bondEnds;
                if (result.ofAtom[neighbor] < 0) {
                    result.ofAtom[neighbor] = system;
                    queue.push_back(neighbor);
                }
            }
        }
        result.atoms.push_back(static_cast<int>(atoms));
        result.simpleRing.push_back(bondEnds / 2 == atoms ? 1 : 0);
    }
    return result;
}

// ============================================================================
// Cycle rank: the number of independent rings, |E| - |V| + components
// ============================================================================

namespace detail {

inline int bondCount(const MolGraph& g) {
    int edgeCount = 0;
    for (int i = 0; i < g.n; i++)
        edgeCount += static_cast<int>(g.neighbors[i].size());
    return edgeCount / 2;
}

inline int connectedComponentCount(const MolGraph& g) {
    std::vector<bool> visited(g.n, false);
    int components = 0;
    for (int i = 0; i < g.n; i++) {
        if (visited[i]) continue;
        components++;
        std::deque<int> bfsQ;
        bfsQ.push_back(i);
        visited[i] = true;
        while (!bfsQ.empty()) {
            int u = bfsQ.front(); bfsQ.pop_front();
            for (int v : g.neighbors[u])
                if (!visited[v]) { visited[v] = true; bfsQ.push_back(v); }
        }
    }
    return components;
}

inline int cycleRank(const MolGraph& g) {
    return bondCount(g) - g.n + connectedComponentCount(g);
}

} // namespace detail

// ============================================================================
// Smallest set of smallest rings (shortest-path cycle basis)
// ============================================================================

inline const std::vector<std::vector<int>>& MolGraph::computeRings() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (sssrComputed) return sssrRings;
    sssrComputed = true;
    if (n == 0) return sssrRings;

    const int cycleRank = detail::cycleRank(*this);
    if (cycleRank <= 0) return sssrRings;

    // For each edge (start, nb) with nb > start, the shortest path from start to nb that
    // avoids that direct edge closes a candidate ring. A candidate keeps only its size and
    // edge, and its path is traced again when the greedy selection reads it, so memory stays
    // linear in the bond count even when every candidate spans one very large ring. The
    // search buffers are shared by every edge. Such a path stays inside the ring system of
    // the edge, so the search skips other systems without changing the path it finds; a
    // bridge closes no ring, and in a system that is one simple ring the path is the ring.
    const CycleSystems systems = cycleSystems();
    const bool bySystem = !systems.ofAtom.empty();
    std::vector<int> parent(n, -2), queue;
    std::vector<uint32_t> reached(n, 0);
    queue.reserve(n);
    uint32_t stamp = 0;
    auto shortestCycle = [&](int start, int nb, std::vector<int>* path) -> size_t {
        const int system = bySystem ? systems.ofAtom[start] : -1;
        if (bySystem) {
            if (system < 0 || systems.ofAtom[nb] != system) return 0;
            if (!path && systems.simpleRing[system]) return static_cast<size_t>(systems.atoms[system]);
        }
        if (++stamp == 0) { std::fill(reached.begin(), reached.end(), 0u); stamp = 1; }
        reached[start] = stamp;
        parent[start] = -1;
        queue.clear();
        queue.push_back(start);
        bool found = false;
        for (size_t head = 0; head < queue.size() && !found; ++head) {
            const int u = queue[head];
            for (int v : neighbors[u]) {
                if (u == start && v == nb) continue;
                if (v == start && u == nb) continue;
                if (reached[v] == stamp) continue;
                if (bySystem && systems.ofAtom[v] != system) continue;
                reached[v] = stamp;
                parent[v] = u;
                if (v == nb) { found = true; break; }
                queue.push_back(v);
            }
        }
        if (!found) return 0;
        size_t size = 0;
        if (path) path->clear();
        int cur = nb;
        while (cur != start && cur != -1) {
            if (path) path->push_back(cur);
            ++size;
            cur = parent[cur];
        }
        if (cur != start) return 0;
        if (path) path->push_back(start);
        return size + 1;
    };
    struct Candidate { size_t size; int start, nb; };
    std::vector<Candidate> candidates;
    for (int start = 0; start < n; start++) {
        for (int nb : neighbors[start]) {
            if (nb <= start) continue;
            const size_t size = shortestCycle(start, nb, nullptr);
            if (size >= 3) candidates.push_back({size, start, nb});
        }
    }

    // The comparison sees only sizes, so the order equals that of sorting the paths themselves.
    // Equal sizes keep the order of their bonds, as in the Java implementation, so every
    // standard library selects the same rings.
    std::stable_sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b) {
            return a.size < b.size;
        });

    // Greedy basis: pick cycles that contribute at least one new edge
    std::vector<std::vector<int>> basis;
    std::unordered_set<int64_t> allBasis;
    std::vector<int> cycle;
    std::vector<int64_t> edges;

    for (const auto& candidate : candidates) {
        if (static_cast<int>(basis.size()) >= cycleRank) break;
        // A simple ring system yields one cycle, so once any of its bonds is in the basis
        // every candidate of that system adds no new bond.
        if (bySystem && systems.simpleRing[systems.ofAtom[candidate.start]] && !basis.empty()
            && allBasis.count(bondKey(candidate.start, candidate.nb))) continue;
        shortestCycle(candidate.start, candidate.nb, &cycle);
        edges.clear();
        for (size_t i = 0; i < cycle.size(); i++) {
            int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
            edges.push_back(bondKey(std::min(a, b), std::max(a, b)));
        }
        bool hasNewEdge = false;
        for (int64_t e : edges) {
            if (allBasis.find(e) == allBasis.end()) { hasNewEdge = true; break; }
        }
        if (hasNewEdge || basis.empty()) {
            basis.push_back(cycle);
            allBasis.insert(edges.begin(), edges.end());
        }
    }

    sssrRings = std::move(basis);
    return sssrRings;
}

// ============================================================================
// Relevant cycles (union of all minimum cycle bases) and unique ring families
// ============================================================================

inline std::vector<std::vector<std::vector<int>>> MolGraph::groupUniqueRingFamilies(
        const std::vector<std::vector<int>>& cycles,
        const std::unordered_map<int64_t, int>& edgeIndex) {
    const int edgeCount = static_cast<int>(edgeIndex.size());
    const int wordCount = (edgeCount + 63) >> 6;
    std::vector<std::vector<uint64_t>> vectors(cycles.size(), std::vector<uint64_t>(wordCount, 0));
    for (size_t cycleIndex = 0; cycleIndex < cycles.size(); ++cycleIndex) {
        const auto& cycle = cycles[cycleIndex];
        for (size_t atom = 0; atom < cycle.size(); ++atom) {
            const int edge = edgeIndex.at(bondKey(cycle[atom], cycle[(atom + 1) % cycle.size()]));
            vectors[cycleIndex][edge >> 6] ^= uint64_t(1) << (edge & 63);
        }
    }

    std::vector<int> order(cycles.size()), family(cycles.size());
    std::iota(order.begin(), order.end(), 0);
    std::iota(family.begin(), family.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int first, int second) {
        return cycles[first].size() < cycles[second].size();
    });
    std::vector<std::vector<uint64_t>> shorterBasis;
    int basisSize = 0;
    size_t begin = 0;
    while (begin < order.size()) {
        size_t end = begin + 1;
        while (end < order.size() && cycles[order[end]].size() == cycles[order[begin]].size()) ++end;
        for (size_t first = begin; first < end; ++first) {
            for (size_t second = first + 1; second < end; ++second) {
                const int firstCycle = order[first], secondCycle = order[second];
                bool sharedEdge = false;
                std::vector<uint64_t> difference(wordCount);
                for (int word = 0; word < wordCount; ++word) {
                    sharedEdge = sharedEdge || ((vectors[firstCycle][word] & vectors[secondCycle][word]) != 0);
                    difference[word] = vectors[firstCycle][word] ^ vectors[secondCycle][word];
                }
                if (sharedEdge && reduceVectorGF2(difference, shorterBasis, basisSize))
                    detail::unionFindJoin(family, firstCycle, secondCycle);
            }
        }
        for (size_t cycle = begin; cycle < end; ++cycle) {
            auto reduced = vectors[order[cycle]];
            basisSize = addToBasisGF2(reduced, shorterBasis, basisSize);
        }
        begin = end;
    }

    std::map<int, std::vector<std::vector<int>>> grouped;
    for (size_t cycle = 0; cycle < cycles.size(); ++cycle)
        grouped[detail::unionFindRoot(family, static_cast<int>(cycle))].push_back(cycles[cycle]);
    std::vector<std::vector<std::vector<int>>> result;
    for (auto& entry : grouped) result.push_back(std::move(entry.second));
    return result;
}

inline int MolGraph::findPivotBit(const std::vector<uint64_t>& vec) {
    for (size_t w = 0; w < vec.size(); w++) {
        if (vec[w] != 0) return static_cast<int>((w << 6) + smsd::ctz64(vec[w]));
    }
    return -1;
}

inline bool MolGraph::reduceVectorGF2(std::vector<uint64_t>& vec,
                            const std::vector<std::vector<uint64_t>>& basis,
                            int basisSize) {
    for (int i = 0; i < basisSize; i++) {
        int pivot = findPivotBit(basis[i]);
        if (pivot >= 0 && ((vec[pivot >> 6] >> (pivot & 63)) & 1ULL) != 0) {
            for (size_t w = 0; w < vec.size(); w++) vec[w] ^= basis[i][w];
        }
    }
    for (auto w : vec) if (w != 0) return false;
    return true;
}

inline int MolGraph::addToBasisGF2(std::vector<uint64_t>& vec,
                         std::vector<std::vector<uint64_t>>& basis,
                         int basisSize) {
    for (int i = 0; i < basisSize; i++) {
        int pivot = findPivotBit(basis[i]);
        if (pivot >= 0 && ((vec[pivot >> 6] >> (pivot & 63)) & 1ULL) != 0) {
            for (size_t w = 0; w < vec.size(); w++) vec[w] ^= basis[i][w];
        }
    }
    bool nonZero = false;
    for (auto w : vec) if (w != 0) { nonZero = true; break; }
    if (nonZero) {
        if (basisSize >= static_cast<int>(basis.size())) basis.push_back(vec);
        else basis[basisSize] = vec;
        return basisSize + 1;
    }
    return basisSize;
}

inline void MolGraph::canonicalCycleKeyInto(const std::vector<int>& cycle, std::vector<int64_t>& keys) {
    keys.clear();
    keys.reserve(cycle.size());
    for (size_t i = 0; i < cycle.size(); i++) {
        int a = cycle[i], b = cycle[(i + 1) % cycle.size()];
        int lo = std::min(a, b), hi = std::max(a, b);
        keys.push_back((static_cast<int64_t>(lo) << 32) | static_cast<int64_t>(hi));
    }
    std::sort(keys.begin(), keys.end());
}

inline std::vector<int64_t> MolGraph::canonicalCycleKey(const std::vector<int>& cycle) {
    std::vector<int64_t> keys;
    canonicalCycleKeyInto(cycle, keys);
    return keys;
}

inline size_t MolGraph::allShortestPathsCpp(int root, int target, const std::vector<int>& dist,
                                  const std::vector<std::vector<int>>& parents, std::vector<int>& current,
                                  std::vector<PathFrame>& stack, std::vector<std::vector<int>>& paths,
                                  bool& cut) {
    constexpr size_t limit = RELEVANT_CYCLE_PATHS;
    size_t count = 0;
    auto emit = [&] {
        if (count == paths.size()) paths.emplace_back();
        paths[count++].assign(current.rbegin(), current.rend());
    };
    current.assign(1, target);
    if (target == root) { emit(); return count; }
    stack.clear();
    stack.push_back({target, 0});
    while (!stack.empty() && count < limit) {
        PathFrame& frame = stack.back();
        const auto& candidates = parents[frame.node];
        while (frame.nextParent < candidates.size()
               && dist[candidates[frame.nextParent]] != dist[frame.node] - 1) ++frame.nextParent;
        if (frame.nextParent == candidates.size()) {
            stack.pop_back();
            current.pop_back();
            continue;
        }
        const int next = candidates[frame.nextParent++];
        current.push_back(next);
        if (next == root) {
            emit();
            current.pop_back();
            continue;
        }
        stack.push_back({next, 0});
    }
    // Every atom but the root has a parent one step nearer it, so a parent not yet followed
    // starts a path not yet listed.
    if (count == limit)
        for (const PathFrame& frame : stack)
            for (size_t parent = frame.nextParent; parent < parents[frame.node].size() && !cut; ++parent)
                cut = dist[parents[frame.node][parent]] == dist[frame.node] - 1;
    return count;
}

inline std::vector<std::vector<int>> MolGraph::computeRelevantCycles() const {
    bool truncated = false;
    return computeRelevantCyclesReportingTruncation(truncated);
}

inline std::vector<std::vector<int>> MolGraph::computeRelevantCyclesReportingTruncation(bool& truncated) const {
    truncated = false;
    auto& mcb = computeRings();
    if (mcb.empty()) return {};

    // Build edge index
    std::unordered_map<int64_t, int> edgeIndex;
    for (int i = 0; i < n; i++) {
        for (int j : neighbors[i]) {
            if (j > i) {
                int64_t key = bondKey(std::min(i, j), std::max(i, j));
                if (edgeIndex.find(key) == edgeIndex.end()) {
                    int idx = static_cast<int>(edgeIndex.size());
                    edgeIndex[key] = idx;
                }
            }
        }
    }
    int numEdges = static_cast<int>(edgeIndex.size());
    int longWords = (numEdges + 63) >> 6;

    // Find max cycle length in MCB
    int maxLen = 0;
    for (auto& c : mcb) if (static_cast<int>(c.size()) > maxLen) maxLen = static_cast<int>(c.size());

    // Collect all short cycles
    std::set<std::vector<int64_t>> uniqueCycleSet;
    std::vector<std::vector<int>> allCycles;

    // Add MCB cycles
    for (auto& c : mcb) {
        auto key = canonicalCycleKey(c);
        if (uniqueCycleSet.insert(key).second) allCycles.push_back(c);
    }

    // Enumerate short cycles via BFS from each vertex. The search buffers are reused for
    // every root and candidate, each root's shortest paths to an atom are listed once, and
    // a candidate is copied only when it is new, so concurrent ring searches in a batch do
    // not contend on the allocator. The reuse leaves the candidates and their order as they are.
    std::vector<int> dist(n, -1);
    std::vector<std::vector<int>> parents(n);
    std::vector<int> queue;
    queue.reserve(n);
    std::vector<std::vector<std::vector<int>>> pathsTo(n);
    std::vector<size_t> pathCount(n, 0);
    std::vector<uint8_t> pathsListed(n, 0);
    std::vector<int> pathScratch;
    std::vector<PathFrame> pathStack;
    std::vector<uint32_t> mark(n, 0);
    uint32_t stamp = 0;
    auto nextStamp = [&] {
        if (++stamp == 0) { std::fill(mark.begin(), mark.end(), 0u); stamp = 1; }
        return stamp;
    };
    uint64_t work = 0;
    std::vector<int> cycle;
    std::vector<int64_t> cycleKey;
    auto rememberCycle = [&] {
        if (static_cast<int>(cycle.size()) < 3 || static_cast<int>(cycle.size()) > maxLen) return;
        work += cycle.size();
        canonicalCycleKeyInto(cycle, cycleKey);
        if (uniqueCycleSet.find(cycleKey) != uniqueCycleSet.end()) return;
        uniqueCycleSet.insert(cycleKey);
        allCycles.push_back(cycle);
    };
    // A tree edge whose farther atom has no other parent closes no cycle: every shortest
    // path to that atom runs through the nearer one, so its pairs all overlap.
    auto onlyParent = [&](int atom, int parent) {
        for (int other : parents[atom]) if (other != parent) return false;
        return true;
    };
    // Cycles through a root lie in its ring system, and the search there reaches its atoms
    // with the same distances and parents, so the search stays inside it. A root outside
    // every system, or in a system that is one simple ring already in the basis, closes no
    // new cycle.
    const CycleSystems systems = cycleSystems();
    const bool bySystem = !systems.ofAtom.empty();
    for (int root = 0; root < n && work < RELEVANT_CYCLE_WORK_LIMIT; root++) {
        const int system = bySystem ? systems.ofAtom[root] : -1;
        if (bySystem && (system < 0 || systems.simpleRing[system])) continue;
        std::fill(dist.begin(), dist.end(), -1);
        for (auto& atomParents : parents) atomParents.clear();
        std::fill(pathsListed.begin(), pathsListed.end(), uint8_t(0));
        work += static_cast<uint64_t>(n);
        // Paths from root to atom are pathsTo[atom][0, pathCount[atom]).
        auto shortestPaths = [&](int atom) {
            if (!pathsListed[atom]) {
                pathCount[atom] = allShortestPathsCpp(root, atom, dist, parents, pathScratch, pathStack, pathsTo[atom],
                                                      truncated);
                pathsListed[atom] = 1;
                work += pathCount[atom] * static_cast<uint64_t>(dist[atom] + 1);
            }
            return pathCount[atom];
        };
        dist[root] = 0;
        queue.clear();
        queue.push_back(root);
        for (size_t head = 0; head < queue.size(); ++head) {
            int u = queue[head];
            if (dist[u] >= maxLen / 2) continue;
            work += 1 + neighbors[u].size();
            for (int v : neighbors[u]) {
                if (bySystem && systems.ofAtom[v] != system) continue;
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    parents[v].push_back(u);
                    queue.push_back(v);
                } else if (dist[v] == dist[u] + 1) {
                    parents[v].push_back(u);
                }
            }
        }

        // Odd-length cycles from edges
        for (int u = 0; u < n && work < RELEVANT_CYCLE_WORK_LIMIT; u++) {
            if (dist[u] == -1) continue;
            for (int v : neighbors[u]) {
                if (v <= u || dist[v] == -1) continue;
                int cycleLen = dist[u] + dist[v] + 1;
                if (cycleLen < 3 || cycleLen > maxLen) continue;
                if (dist[u] != dist[v] && (dist[u] > dist[v] ? onlyParent(u, v) : onlyParent(v, u))) continue;

                const size_t pathsU = shortestPaths(u);
                const size_t pathsV = shortestPaths(v);

                for (size_t pathU = 0; pathU < pathsU; ++pathU) {
                    const auto& pu = pathsTo[u][pathU];
                    const uint32_t onPath = nextStamp();
                    for (int x : pu) mark[x] = onPath;
                    work += pu.size();
                    for (size_t pathV = 0; pathV < pathsV; ++pathV) {
                        const auto& pv = pathsTo[v][pathV];
                        bool overlap = false;
                        size_t read = 1;
                        for (; read < pv.size(); read++) {
                            if (mark[pv[read]] == onPath) { overlap = true; break; }
                        }
                        work += read;
                        if (overlap) continue;

                        cycle.clear();
                        for (int x : pu) cycle.push_back(x);
                        for (int k = static_cast<int>(pv.size()) - 1; k >= 1; k--)
                            cycle.push_back(pv[k]);
                        rememberCycle();
                    }
                }
            }
        }

        // Even-length cycles from two disjoint paths to same vertex
        for (int u = 0; u < n && work < RELEVANT_CYCLE_WORK_LIMIT; u++) {
            if (dist[u] == -1 || static_cast<int>(parents[u].size()) < 2) continue;
            int cycleLen = 2 * dist[u];
            if (cycleLen < 3 || cycleLen > maxLen) continue;

            const size_t pathCountU = shortestPaths(u);
            const auto& paths = pathsTo[u];
            for (size_t pi = 0; pi < pathCountU; pi++) {
                const auto& p1 = paths[pi];
                const uint32_t interior = nextStamp();
                for (size_t k = 1; k + 1 < p1.size(); k++) mark[p1[k]] = interior;
                work += p1.size();
                for (size_t pj = pi + 1; pj < pathCountU; pj++) {
                    const auto& p2 = paths[pj];
                    bool disjoint = true;
                    size_t read = 1;
                    for (; read + 1 < p2.size(); read++) {
                        if (mark[p2[read]] == interior) { disjoint = false; break; }
                    }
                    work += read;
                    if (!disjoint) continue;

                    cycle.clear();
                    for (size_t k = 0; k < p1.size(); k++) cycle.push_back(p1[k]);
                    for (int k = static_cast<int>(p2.size()) - 2; k >= 1; k--)
                        cycle.push_back(p2[k]);
                    rememberCycle();
                }
            }
        }
    }

    // With the work allowance used up, roots or candidates may have been left unsearched.
    if (work >= RELEVANT_CYCLE_WORK_LIMIT) truncated = true;

    // Sort by length; cycles of one length keep the order in which they were found.
    std::stable_sort(allCycles.begin(), allCycles.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() < b.size();
        });

    // Build GF(2) edge-incidence vectors
    std::vector<std::vector<uint64_t>> vectors(allCycles.size(),
        std::vector<uint64_t>(longWords, 0));
    for (size_t ci = 0; ci < allCycles.size(); ci++) {
        auto& c = allCycles[ci];
        for (size_t i = 0; i < c.size(); i++) {
            int a = c[i], b = c[(i + 1) % c.size()];
            int64_t key = bondKey(std::min(a, b), std::max(a, b));
            int idx = edgeIndex.at(key);
            vectors[ci][idx >> 6] ^= 1ULL << (idx & 63);
        }
    }

    // Filter relevant cycles via Gaussian elimination on GF(2). Each basis row is stored
    // with its lowest edge as pivot, so a vector is reduced by the rows its own lowest edges
    // select, and a lowest edge without a row proves it independent.
    std::vector<std::vector<uint64_t>> rows;
    std::vector<int> rowByPivot(numEdges, -1);
    std::vector<uint64_t> reduced(longWords);
    auto reduce = [&](const std::vector<uint64_t>& cycleBits) {
        reduced = cycleBits;
        for (int word = 0; word < longWords; ++word) {
            while (reduced[word] != 0) {
                const int pivot = (word << 6) + smsd::ctz64(reduced[word]);
                const int row = rowByPivot[pivot];
                if (row < 0) return pivot;
                const auto& basisRow = rows[row];
                for (int other = word; other < longWords; ++other) reduced[other] ^= basisRow[other];
            }
        }
        return -1;
    };
    std::vector<std::vector<int>> relevant;

    size_t ci = 0;
    while (ci < allCycles.size()) {
        int len = static_cast<int>(allCycles[ci].size());
        size_t cj = ci;
        while (cj < allCycles.size() && static_cast<int>(allCycles[cj].size()) == len) cj++;

        // Test each cycle against basis of strictly shorter cycles
        for (size_t k = ci; k < cj; k++)
            if (reduce(vectors[k]) >= 0) relevant.push_back(allCycles[k]);

        // Add all cycles of this length to basis
        for (size_t k = ci; k < cj; k++) {
            const int pivot = reduce(vectors[k]);
            if (pivot < 0) continue;
            rowByPivot[pivot] = static_cast<int>(rows.size());
            rows.push_back(reduced);
        }

        ci = cj;
    }

    return relevant;
}

inline std::vector<std::vector<std::vector<int>>> MolGraph::computeURFs() {
    std::unordered_map<int64_t, int> edgeIndex;
    for (int atom = 0; atom < n; ++atom) {
        for (int neighbor : neighbors[atom]) {
            if (atom < neighbor) {
                const int index = static_cast<int>(edgeIndex.size());
                edgeIndex.emplace(bondKey(atom, neighbor), index);
            }
        }
    }
    return groupUniqueRingFamilies(computeRelevantCycles(), edgeIndex);
}

// ============================================================================
// Lazy ring counts and ring systems
// ============================================================================

inline void MolGraph::ensureRingCounts() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (ringCountsComputed_) return;
    computeRingCounts();
    ringCountsComputed_ = true;
}

inline void MolGraph::ensureRingSystems() const {
    const std::lock_guard<std::recursive_mutex> cacheLock(lazyCacheMutex());
    if (ringSysComputed_) return;
    ringSystemId_.assign(n, -1);

    // Union-find parent array (only for ring atoms)
    std::vector<int> parent(n, -1);

    // Initialise ring atoms as their own parent
    for (int i = 0; i < n; ++i) {
        if (ring[i]) parent[i] = i;
    }

    // Union ring atoms connected by edges where both endpoints are ring atoms
    for (int i = 0; i < n; ++i) {
        if (!ring[i]) continue;
        for (int nb : neighbors[i]) {
            if (nb > i && ring[nb]) detail::unionFindJoin(parent, i, nb);
        }
    }

    // Assign contiguous ring system IDs
    std::unordered_map<int, int> root2id;
    ringSysCount_ = 0;
    for (int i = 0; i < n; ++i) {
        if (!ring[i]) continue;
        int root = detail::unionFindRoot(parent, i);
        auto it = root2id.find(root);
        if (it == root2id.end()) {
            root2id[root] = ringSysCount_;
            ringSystemId_[i] = ringSysCount_++;
        } else {
            ringSystemId_[i] = it->second;
        }
    }

    // Compute signature per ring system
    ringSysSignature_.assign(ringSysCount_, 0);
    if (ringSysCount_ > 0) {
        // Collect per-system: sorted atomic numbers, bond orders, size, aromatic count
        std::vector<std::vector<int>> sysAtomicNums(ringSysCount_);
        std::vector<std::vector<int>> sysBondOrders(ringSysCount_);
        std::vector<int> sysSize(ringSysCount_, 0);
        std::vector<int> sysAromCount(ringSysCount_, 0);

        for (int i = 0; i < n; ++i) {
            int rsid = ringSystemId_[i];
            if (rsid < 0) continue;
            sysAtomicNums[rsid].push_back(atomicNum[i]);
            sysSize[rsid]++;
            if (aromatic[i]) sysAromCount[rsid]++;

            // Collect bond orders for edges within this ring system (each edge once)
            for (int nb : neighbors[i]) {
                if (nb > i && ringSystemId_[nb] == rsid) {
                    sysBondOrders[rsid].push_back(bondOrder(i, nb));
                }
            }
        }

        for (int rs = 0; rs < ringSysCount_; ++rs) {
            std::sort(sysAtomicNums[rs].begin(), sysAtomicNums[rs].end());
            std::sort(sysBondOrders[rs].begin(), sysBondOrders[rs].end());

            // FNV-1a style hash
            uint64_t h = 14695981039346656037ULL;
            auto mix = [&](uint64_t val) {
                h ^= val;
                h *= 1099511628211ULL;
            };
            mix(static_cast<uint64_t>(sysSize[rs]));
            mix(static_cast<uint64_t>(sysAromCount[rs]));
            for (int z : sysAtomicNums[rs]) mix(static_cast<uint64_t>(z));
            for (int bo : sysBondOrders[rs]) mix(static_cast<uint64_t>(bo));
            ringSysSignature_[rs] = h;
        }
    }

    ringSysComputed_ = true;
}

inline int MolGraph::ringSystemOf(int atom) const {
    if (atom < 0 || atom >= n) throw std::out_of_range("Atom index out of range");
    ensureRingSystems();
    return ringSystemId_[atom];
}

inline uint64_t MolGraph::ringSystemSig(int rsId) const {
    ensureRingSystems();
    if (rsId < 0 || rsId >= ringSysCount_) return 0;
    return ringSysSignature_[rsId];
}

} // namespace smsd

#include "smsd/ring_finder.hpp"

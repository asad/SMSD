/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/detail/search_budget.hpp"
#include "smsd/detail/search_compatibility.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <numeric>
#include <vector>

namespace smsd {
namespace search_detail {

/**
 * An atom bound for connected mappings from true bridges. Roots with only bridge
 * bonds pair their disjoint branches injectively; cyclic roots use part counts.
 * Actual predicates are used at roots and links, without raw-label equality assumptions.
 */
class BridgeBound {
    static constexpr int MAX_BRANCHES = 8;
    static constexpr int MAX_ATOMS = 512;
    static constexpr size_t MAX_PART_PAIRS = 250000;
    static constexpr size_t MAX_COUNTS_PER_GRAPH = 1000000;
    static constexpr uint64_t MAX_WORK = 50000000;
    struct Part {
        int root = -1, parent = -1, size = 0;
        std::vector<int> branches;
        std::vector<int> atomCounts;
    };
    struct Parts {
        std::vector<Part> parts;
        std::vector<int> order;
        std::vector<uint8_t> chainRoot;
        bool acyclic = true, hasBridge = false;

        bool build(const MolGraph& graph, const std::vector<int>& elements, Budget& budget) {
            if (graph.n <= 0 || graph.n > MAX_ATOMS || graph.neighbors.size() != static_cast<size_t>(graph.n))
                return false;
            std::vector<int> parent(static_cast<size_t>(graph.n), -2), queue;
            queue.reserve(static_cast<size_t>(graph.n));
            for (int atom = 0; atom < graph.n; ++atom) {
                if (budget.stopped() || graph.neighbors[atom].size() > MAX_BRANCHES) return false;
                const auto& neighbors = graph.neighbors[atom];
                for (size_t place = 0; place < neighbors.size(); ++place) {
                    const int neighbor = neighbors[place];
                    if (neighbor < 0 || neighbor >= graph.n || neighbor == atom
                        || std::find(neighbors.begin(), neighbors.begin() + place, neighbor) != neighbors.begin() + place)
                        return false;
                }
                // Certify the topology actually read by the induced search, including dense graphs.
                for (int other = 0; other < graph.n; ++other) {
                    const bool listed = std::find(neighbors.begin(), neighbors.end(), other) != neighbors.end();
                    if (graph.hasBond(atom, other) != listed
                        || graph.bondOrder(atom, other) != graph.bondOrder(other, atom)) return false;
                    if (listed && (graph.bondAromatic(atom, other) != graph.bondAromatic(other, atom)
                                   || graph.bondInRing(atom, other) != graph.bondInRing(other, atom))) return false;
                    if (listed && std::find(graph.neighbors[other].begin(), graph.neighbors[other].end(), atom)
                                      == graph.neighbors[other].end()) return false;
                }
            }
            std::vector<int> entered(static_cast<size_t>(graph.n), -1), low(static_cast<size_t>(graph.n)),
                next(static_cast<size_t>(graph.n));
            std::vector<std::vector<uint8_t>> bridge(static_cast<size_t>(graph.n));
            for (int atom = 0; atom < graph.n; ++atom) bridge[atom].resize(graph.neighbors[atom].size());
            int clock = 0;
            for (int start = 0; start < graph.n; ++start) {
                if (entered[start] >= 0) continue;
                queue.clear();
                queue.push_back(start);
                parent[start] = -1;
                entered[start] = low[start] = clock++;
                while (!queue.empty()) {
                    if (budget.stopped()) return false;
                    const int atom = queue.back();
                    if (next[atom] < static_cast<int>(graph.neighbors[atom].size())) {
                        const int neighbor = graph.neighbors[atom][next[atom]++];
                        if (neighbor == parent[atom]) continue;
                        if (entered[neighbor] < 0) {
                            parent[neighbor] = atom;
                            entered[neighbor] = low[neighbor] = clock++;
                            queue.push_back(neighbor);
                        } else low[atom] = std::min(low[atom], entered[neighbor]);
                    } else {
                        queue.pop_back();
                        const int previous = parent[atom];
                        if (previous < 0) continue;
                        low[previous] = std::min(low[previous], low[atom]);
                        if (low[atom] > entered[previous]) {
                            const auto& adjacent = graph.neighbors[previous];
                            bridge[previous][std::find(adjacent.begin(), adjacent.end(), atom) - adjacent.begin()] = 1;
                            const auto& reverse = graph.neighbors[atom];
                            bridge[atom][std::find(reverse.begin(), reverse.end(), previous) - reverse.begin()] = 1;
                            hasBridge = true;
                        }
                    }
                }
            }
            if (!hasBridge)
                for (const auto& adjacent : graph.neighbors) if (!adjacent.empty()) return false;
            chainRoot.assign(static_cast<size_t>(graph.n), 1);
            std::vector<std::vector<int>> beyond(static_cast<size_t>(graph.n));
            parts.reserve(static_cast<size_t>(graph.n) * 3);
            for (int atom = 0; atom < graph.n; ++atom) parts.push_back(Part{atom, -1, 0, {}, {}});
            for (int atom = 0; atom < graph.n; ++atom) {
                beyond[atom].assign(graph.neighbors[atom].size(), -1);
                for (size_t place = 0; place < graph.neighbors[atom].size(); ++place) {
                    if (bridge[atom][place]) {
                        beyond[atom][place] = static_cast<int>(parts.size());
                        parts.push_back(Part{graph.neighbors[atom][place], atom, 0, {}, {}});
                    } else {
                        chainRoot[atom] = 0;
                        acyclic = false;
                    }
                }
            }
            if (elements.size() > MAX_COUNTS_PER_GRAPH / parts.size()) return false;
            for (auto& part : parts) {
                if (budget.stopped()) return false;
                const auto& neighbors = graph.neighbors[part.root];
                for (size_t place = 0; place < neighbors.size(); ++place)
                    if (neighbors[place] != part.parent && beyond[part.root][place] >= 0)
                        part.branches.push_back(beyond[part.root][place]);
                part.atomCounts.resize(elements.size());
                queue.clear();
                queue.push_back(part.root);
                std::fill(parent.begin(), parent.end(), -2);
                parent[part.root] = -1;
                for (size_t head = 0; head < queue.size(); ++head) {
                    if ((head & 255u) == 0 && budget.stopped()) return false;
                    const int atom = queue[head];
                    if (!elements.empty()) ++part.atomCounts[std::lower_bound(elements.begin(), elements.end(),
                        graph.atomicNum[atom]) - elements.begin()];
                    for (int neighbor : graph.neighbors[atom])
                        if (!(atom == part.root && neighbor == part.parent) && parent[neighbor] == -2) {
                            parent[neighbor] = atom;
                            queue.push_back(neighbor);
                        }
                }
                part.size = static_cast<int>(queue.size());
            }
            order.resize(parts.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](int first, int second) {
                return parts[first].size < parts[second].size;
            });
            return true;
        }
    };
    size_t targetParts_ = 0;
    std::vector<int> pairBounds_;
    int upper_ = 0;
    bool acyclic_ = false;

public:
    /** Decline graphs without useful bridges, high degree or resource limits; publish only a finished table. */
    static std::unique_ptr<BridgeBound> create(const MolGraph& query, const MolGraph& target,
                                              const ChemOptions& chemistry, Budget& budget) {
        if (query.n <= 0 || target.n <= 0 || query.n > MAX_ATOMS || target.n > MAX_ATOMS
            || query.atomicNum.size() < static_cast<size_t>(query.n)
            || target.atomicNum.size() < static_cast<size_t>(target.n)) return nullptr;
        std::vector<int> elements;
        if (chemistry.matchAtomType) {
            elements.insert(elements.end(), query.atomicNum.begin(), query.atomicNum.begin() + query.n);
            elements.insert(elements.end(), target.atomicNum.begin(), target.atomicNum.begin() + target.n);
            std::sort(elements.begin(), elements.end());
            elements.erase(std::unique(elements.begin(), elements.end()), elements.end());
        }
        Parts queryParts, targetParts;
        if (!queryParts.build(query, elements, budget) || !targetParts.build(target, elements, budget)
            || queryParts.parts.size() > MAX_PART_PAIRS / targetParts.parts.size()) return nullptr;
        std::unique_ptr<BridgeBound> result(new BridgeBound);
        result->acyclic_ = queryParts.acyclic && targetParts.acyclic;
        result->targetParts_ = targetParts.parts.size();
        result->pairBounds_.resize(queryParts.parts.size() * targetParts.parts.size());
        std::array<int, 1 << MAX_BRANCHES> best{};
        uint64_t work = 0;
        for (int queryIndex : queryParts.order) for (int targetIndex : targetParts.order) {
            if (budget.stopped()) return nullptr;
            const Part& queryPart = queryParts.parts[queryIndex];
            const Part& targetPart = targetParts.parts[targetIndex];
            if (!atomCompatible(query, queryPart.root, target, targetPart.root, chemistry)) continue;
            if (elements.size() + 1 > MAX_WORK - work) return nullptr;
            work += elements.size() + 1;
            int bound = std::min(queryPart.size, targetPart.size);
            if (!elements.empty()) {
                bound = 0;
                for (size_t element = 0; element < elements.size(); ++element)
                    bound += std::min(queryPart.atomCounts[element], targetPart.atomCounts[element]);
            }
            if (!queryParts.chainRoot[queryPart.root] || !targetParts.chainRoot[targetPart.root]) {
                result->pairBounds_[static_cast<size_t>(queryIndex) * result->targetParts_ + targetIndex] = bound;
                if (queryIndex < query.n && targetIndex < target.n) result->upper_ = std::max(result->upper_, bound);
                continue;
            }
            const int sets = 1 << targetPart.branches.size();
            const uint64_t units = (queryPart.branches.size() + 1) * static_cast<uint64_t>(sets)
                * (targetPart.branches.size() + 1);
            if (units > MAX_WORK - work) return nullptr;
            work += units;
            std::fill(best.begin(), best.begin() + sets, -1);
            best[0] = 0;
            for (int queryBranch : queryPart.branches) {
                if (budget.stopped()) return nullptr;
                const int queryChild = queryParts.parts[queryBranch].root;
                std::array<int, MAX_BRANCHES> weights{};
                for (size_t branch = 0; branch < targetPart.branches.size(); ++branch) {
                    const int targetBranch = targetPart.branches[branch];
                    if (bondCompatible(query, queryPart.root, queryChild,
                                       target, targetPart.root, targetParts.parts[targetBranch].root, chemistry))
                        weights[branch] = result->pairBounds_[static_cast<size_t>(queryBranch) * result->targetParts_ + targetBranch];
                }
                // Supersets first: each query branch is paired at most once; skipping is retained.
                for (int used = sets - 1; used >= 0; --used) {
                    if ((used & 255) == 255 && budget.stopped()) return nullptr;
                    if (best[used] < 0) continue;
                    for (size_t branch = 0; branch < targetPart.branches.size(); ++branch)
                        if ((used & (1 << branch)) == 0 && weights[branch] != 0) {
                            const int with = used | (1 << branch);
                            best[with] = std::max(best[with], best[used] + weights[branch]);
                        }
                }
            }
            bound = std::min(bound, 1 + *std::max_element(best.begin(), best.begin() + sets));
            result->pairBounds_[static_cast<size_t>(queryIndex) * result->targetParts_ + targetIndex] = bound;
            if (queryIndex < query.n && targetIndex < target.n) result->upper_ = std::max(result->upper_, bound);
        }
        return result;
    }

    int upper() const { return upper_; }
    bool acyclic() const { return acyclic_; }
    /** The full-component bound for a root pair, rather than a directed edge side. */
    int bound(int queryAtom, int targetAtom) const {
        return pairBounds_[static_cast<size_t>(queryAtom) * targetParts_ + targetAtom];
    }
};

} // namespace search_detail
} // namespace smsd

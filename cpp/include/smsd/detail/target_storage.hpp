/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/mol_graph.hpp"

#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

namespace smsd {
namespace batch {
namespace detail {

/** Owned graph snapshots with one slot for every original target. */
class TargetStorage {
    struct Data {
        std::vector<MolGraph> graphs;
        std::vector<size_t> slots;
    };
    std::shared_ptr<Data> data_ = std::make_shared<Data>();

    void own() {
        if (!data_) { data_ = std::make_shared<Data>(); return; }
        if (data_.use_count() == 1) return;
        auto copied = std::make_shared<Data>();
        copied->slots = data_->slots;
        copied->graphs.reserve(data_->graphs.size());
        for (const auto& graph : data_->graphs) {
            const std::lock_guard<std::recursive_mutex> lock(graph.lazyCacheMutex());
            copied->graphs.push_back(graph);
        }
        data_ = std::move(copied);
    }

public:
    size_t size() const { return data_ ? data_->slots.size() : 0; }
    bool empty() const { return size() == 0; }
    const MolGraph& operator[](size_t index) const { return data_->graphs[data_->slots[index]]; }
    size_t storedSize() const { return data_ ? data_->graphs.size() : 0; }
    const MolGraph& storedGraph(size_t index) const { return data_->graphs[index]; }

    /** Capture repeated references once; no input graph remains externally owned. */
    template<class Targets>
    static TargetStorage capture(const Targets& targets) {
        static_assert(std::is_lvalue_reference<decltype(targets[size_t(0)])>::value,
                      "Captured targets must refer to persistent input graphs");
        TargetStorage result;
        result.data_->slots.reserve(targets.size());
        std::unordered_map<const MolGraph*, size_t> captured;
        for (size_t index = 0; index < targets.size(); ++index) {
            const MolGraph& graph = targets[index];
            auto found = captured.find(&graph);
            if (found == captured.end()) {
                const size_t slot = result.data_->graphs.size();
                const std::lock_guard<std::recursive_mutex> lock(graph.lazyCacheMutex());
                result.data_->graphs.push_back(graph);
                found = captured.emplace(&graph, slot).first;
            }
            result.data_->slots.push_back(found->second);
        }
        return result;
    }

    void add(MolGraph graph) {
        own();
        data_->slots.reserve(data_->slots.size() + 1);
        data_->graphs.push_back(std::move(graph));
        data_->slots.push_back(data_->graphs.size() - 1);
    }

    void add(std::vector<MolGraph> graphs) {
        own();
        data_->slots.reserve(data_->slots.size() + graphs.size());
        data_->graphs.reserve(data_->graphs.size() + graphs.size());
        for (auto& graph : graphs) {
            data_->graphs.push_back(std::move(graph));
            data_->slots.push_back(data_->graphs.size() - 1);
        }
    }

    void add(TargetStorage captured) {
        if (captured.empty()) return;
        if (empty()) { data_ = std::move(captured.data_); return; }
        own();
        captured.own();
        const size_t first = data_->graphs.size();
        data_->slots.reserve(data_->slots.size() + captured.size());
        data_->graphs.reserve(first + captured.data_->graphs.size());
        for (auto& graph : captured.data_->graphs) data_->graphs.push_back(std::move(graph));
        for (size_t slot : captured.data_->slots) data_->slots.push_back(first + slot);
    }
};

} // namespace detail
} // namespace batch
} // namespace smsd

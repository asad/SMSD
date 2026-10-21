/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/search.hpp"
#include "smsd/smiles_parser.hpp"

namespace smsd {

struct RGroupResult {
    MolGraph core;
    std::map<std::string, MolGraph> rgroups;
    SearchResult search;
};

/** An immutable molecule or SMARTS query with per-call search state. */
class Substructure {
    std::shared_ptr<const MolGraph> query_;
    std::shared_ptr<const SmartsQuery> pattern_;
    std::shared_ptr<const search_detail::PatternPlan> plan_;
    SubstructureOptions options_;
    explicit Substructure(std::shared_ptr<const MolGraph> query) : query_(std::move(query)) {}
    explicit Substructure(std::shared_ptr<const SmartsQuery> pattern)
        : pattern_(std::move(pattern)), plan_(std::make_shared<const search_detail::PatternPlan>(*pattern_)) {}
public:
    using Options = SubstructureOptions;
    static Substructure find(const MolGraph& query) { return Substructure(search_detail::ownQuery(query)); }
    static Substructure find(const std::string& smiles) { return find(parseSMILES(smiles)); }
    static Substructure smarts(const std::string& pattern) {
        return Substructure(std::make_shared<const SmartsQuery>(parseSMARTS(pattern)));
    }
    static Substructure smarts(const SmartsQuery& pattern) {
        return Substructure(std::make_shared<const SmartsQuery>(pattern));
    }
    Substructure withOptions(const Options& options) const {
        auto result = *this;
        result.options_ = options;
        // A molecule query compared by ring fusion counts its rings once here, not in every search.
        if (result.query_ && options.chemistry.ringFusionMode != ChemOptions::RingFusionMode::IGNORE)
            result.query_->ensureRingCounts();
        return result;
    }
    Substructure withTimeout(int64_t milliseconds) const {
        auto result = *this;
        result.options_.timeoutMs = milliseconds;
        return result;
    }
    SearchResult match(const MolGraph& target) const {
        return search_detail::substructure(query_.get(), pattern_.get(), target, options_, 1, false, nullptr, plan_.get(), true);
    }
    /**
     * Milliseconds from the start of a search of `target` until which the matches it found in
     * time are built: its time limit and an eighth of it, as matchAll() builds them.
     */
    int64_t resultsTimeMillis(const MolGraph& target) const {
        const int queryAtoms = query_ ? query_->n : pattern_ ? pattern_->atomCount() : 0;
        return search_detail::resultsTimeMillis(search_detail::searchTimeLimitMillis(options_, queryAtoms, target.n, false));
    }
    /** The cancellation token of the options, or null; once it is cancelled, building matches stops as well. */
    std::shared_ptr<const CancellationToken> cancellation() const { return options_.cancellation; }
    SearchResult matchAll(const MolGraph& target, size_t limit = 1000) const {
        return search_detail::substructure(query_.get(), pattern_.get(), target, options_, limit, true, nullptr, plan_.get(), true);
    }
    SearchResult match(const std::string& target) const { return match(parseSMILES(target)); }
    template<class Targets>
    std::vector<SearchResult> matchMany(const Targets& targets, int numThreads = 0) const {
        return search_detail::matchMany(*this, targets, numThreads);
    }
    SearchResult matchAll(const std::string& target, size_t limit = 1000) const {
        return matchAll(parseSMILES(target), limit);
    }
    static std::vector<RGroupResult> decomposeRGroups(const MolGraph& core, const std::vector<MolGraph>& molecules,
                                                     const Options& options = Options{}) {
        const auto query = find(core).withOptions(options);
        std::vector<RGroupResult> results;
        results.reserve(molecules.size());
        for (const auto& molecule : molecules) {
            RGroupResult result;
            result.search = query.match(molecule);
            if (!result.search.matches.empty()) {
                Match inverse;
                for (const auto& atom : result.search.matches.front().atomMapping) inverse.atomMapping.emplace(atom.second, atom.first);
                for (const auto& bond : result.search.matches.front().bondMapping)
                    inverse.bondMapping.push_back({bond[2], bond[3], bond[0], bond[1]});
                result.core = inverse.toGraph(molecule);
                std::vector<bool> visited(molecule.n, false);
                for (const auto& atom : inverse.atomMapping) visited[atom.first] = true;
                int number = 0;
                for (const auto& atom : inverse.atomMapping) for (int neighbor : molecule.neighbors[atom.first]) {
                    if (visited[neighbor]) continue;
                    visited[neighbor] = true;
                    std::vector<int> component{neighbor};
                    for (size_t head = 0; head < component.size(); ++head)
                        for (int adjacent : molecule.neighbors[component[head]])
                            if (!visited[adjacent]) { visited[adjacent] = true; component.push_back(adjacent); }
                    std::sort(component.begin(), component.end());
                    result.rgroups.emplace("R" + std::to_string(++number), extractSubgraph(molecule, component));
                }
            }
            results.push_back(std::move(result));
        }
        return results;
    }
};

} // namespace smsd

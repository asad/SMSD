/* SPDX-License-Identifier: Apache-2.0 */
#pragma once

#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/rascal.hpp"

namespace smsd_test {

inline smsd::SubstructureOptions substructureOptions(const smsd::ChemOptions& chemistry, int64_t timeoutMs = 10000) {
    smsd::SubstructureOptions options;
    options.chemistry = chemistry; options.timeoutMs = timeoutMs;
    return options;
}

inline bool substructureExists(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                const smsd::ChemOptions& chemistry = {}, int64_t timeoutMs = 10000) {
    return !smsd::Substructure::find(query).withOptions(substructureOptions(chemistry, timeoutMs)).match(target).matches.empty();
}

inline std::vector<std::pair<int, int>> substructureAtoms(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                                        const smsd::ChemOptions& chemistry = {}, int64_t timeoutMs = 10000) {
    const auto result = smsd::Substructure::find(query).withOptions(substructureOptions(chemistry, timeoutMs)).match(target);
    if (result.matches.empty()) return {};
    return {result.matches.front().atomMapping.begin(), result.matches.front().atomMapping.end()};
}

inline std::vector<std::vector<std::pair<int, int>>> allSubstructureAtoms(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                                                       const smsd::ChemOptions& chemistry = {}, int64_t timeoutMs = 10000,
                                                                       size_t limit = 1000) {
    const auto result = smsd::Substructure::find(query).withOptions(substructureOptions(chemistry, timeoutMs)).matchAll(target, limit);
    std::vector<std::vector<std::pair<int, int>>> mappings;
    for (const auto& match : result.matches) mappings.emplace_back(match.atomMapping.begin(), match.atomMapping.end());
    return mappings;
}

inline std::map<int, int> commonAtoms(const smsd::MolGraph& query, const smsd::MolGraph& target,
                                     const smsd::ChemOptions& chemistry = {}, smsd::MCSOptions options = {}) {
    options.chemistry = chemistry;
    const auto result = smsd::MCS::find(query).withOptions(options).match(target);
    return result.matches.empty() ? std::map<int, int>{} : result.matches.front().atomMapping;
}

inline bool patternMatches(const smsd::SmartsQuery& query, const smsd::MolGraph& target) {
    return !smsd::Substructure::smarts(query).match(target).matches.empty();
}

inline std::vector<std::map<int, int>> allPatternAtoms(const smsd::SmartsQuery& query, const smsd::MolGraph& target, size_t limit = 1000) {
    const auto result = smsd::Substructure::smarts(query).matchAll(target, limit);
    std::vector<std::map<int, int>> mappings;
    for (const auto& match : result.matches) mappings.push_back(match.atomMapping);
    return mappings;
}

inline bool anchoredPatternMatches(const smsd::SmartsQuery& query, const smsd::MolGraph& target, int anchor) {
    const auto result = smsd::Substructure::smarts(query).matchAll(target);
    for (const auto& match : result.matches) {
        const auto atom = match.atomMapping.find(0);
        if (atom != match.atomMapping.end() && atom->second == anchor) return true;
    }
    return false;
}

inline std::string fragmentSmiles(const smsd::MolGraph& graph, const std::map<int, int>& mapping) {
    smsd::Match selected;
    selected.atomMapping = mapping;
    for (const auto& atom : mapping) for (int neighbor : graph.neighbors.at(atom.first))
        if (atom.first < neighbor && mapping.count(neighbor))
            selected.bondMapping.push_back({atom.first, neighbor, atom.second, mapping.at(neighbor)});
    return selected.toSmiles(graph);
}

inline std::vector<bool> substructureBatch(const smsd::MolGraph& query, const std::vector<smsd::MolGraph>& targets,
                                          const smsd::ChemOptions& chemistry = {}, int threads = 0) {
    const auto results = smsd::Substructure::find(query).withOptions(substructureOptions(chemistry)).matchMany(targets, threads);
    std::vector<bool> hits;
    for (const auto& result : results) hits.push_back(!result.matches.empty());
    return hits;
}

inline std::vector<std::map<int, int>> commonBatch(const smsd::MolGraph& query, const std::vector<smsd::MolGraph>& targets,
                                                 const smsd::ChemOptions& chemistry = {}, smsd::MCSOptions options = {}, int threads = 0) {
    options.chemistry = chemistry;
    const auto results = smsd::MCS::find(query).withOptions(options).matchMany(targets, threads);
    std::vector<std::map<int, int>> mappings;
    for (const auto& result : results) mappings.push_back(result.matches.empty() ? std::map<int, int>{} : result.matches.front().atomMapping);
    return mappings;
}

inline std::vector<std::pair<int, std::map<int, int>>> screenedCommon(const smsd::MolGraph& query,
        const std::vector<smsd::MolGraph>& targets, const smsd::ChemOptions& chemistry = {}, smsd::MCSOptions options = {},
        double threshold = 0.5, int threads = 0) {
    options.chemistry = chemistry;
    const auto prepared = smsd::MCS::find(query).withOptions(options);
    const auto indices = smsd::screenTargets(query, targets, threshold, chemistry);
    std::vector<smsd::MolGraph> screened;
    for (int index : indices) screened.push_back(targets[index]);
    const auto results = prepared.matchMany(screened, threads);
    std::vector<std::pair<int, std::map<int, int>>> mappings;
    for (size_t index = 0; index < results.size(); ++index)
        if (!results[index].matches.empty()) mappings.emplace_back(indices[index], results[index].matches.front().atomMapping);
    return mappings;
}

inline std::vector<std::pair<int, int>> screenedCommonSizes(const smsd::MolGraph& query,
        const std::vector<smsd::MolGraph>& targets, const smsd::ChemOptions& chemistry = {}, smsd::MCSOptions options = {},
        double threshold = 0.5, int threads = 0) {
    const auto results = screenedCommon(query, targets, chemistry, options, threshold, threads);
    std::vector<std::pair<int, int>> sizes;
    for (const auto& result : results) sizes.emplace_back(result.first, static_cast<int>(result.second.size()));
    return sizes;
}

} // namespace smsd_test

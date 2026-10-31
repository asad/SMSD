/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/search.hpp"
#include "smsd/batch.hpp"
#include "smsd/smiles_parser.hpp"

namespace smsd {

/** An immutable maximum common subgraph query. */
class MCS {
    std::shared_ptr<const MolGraph> query_;
    MCSOptions options_;
    explicit MCS(std::shared_ptr<const MolGraph> query) : query_(std::move(query)) {}
public:
    using Options = MCSOptions;
    using Objective = MCSObjective;
    using Connectivity = MCSConnectivity;
    /** A search result and a portable SMARTS query over its best selected graph. */
    struct SmartsResult {
        SearchResult result;
        std::string smarts;
    };
    /** A progressive common graph and its supporting embeddings; never a proof of optimality. */
    struct ConsensusResult {
        Match fragment;
        int sourceIndex = -1;
        MolGraph graph;
        std::map<int, Match> mappings;
        int supportCount = 0;
        bool found = false;
        SearchStatus status = SearchStatus::Completed;
        bool optimal = false;
        bool complete = false;
        uint64_t nodesVisited = 0;
        int64_t elapsedMillis = 0;
        std::string toSmiles() const { return smiles_; }
    private:
        friend class MCS;
        std::string smiles_;
    };
    static MCS find(const MolGraph& query) { return MCS(search_detail::ownQuery(query)); }
    static MCS find(const std::string& smiles) { return find(parseSMILES(smiles)); }
    MCS withOptions(const Options& options) const {
        search_detail::validateOptions(options, query_->n);
        auto result = *this;
        result.options_ = options;
        return result;
    }
    MCS withTimeout(int64_t milliseconds) const {
        auto result = *this;
        result.options_.timeoutMs = milliseconds;
        return result;
    }
    /**
     * Milliseconds from the start of a search of `target` until which the mappings it found in
     * time are built: its time limit and an eighth of it, as Substructure::resultsTimeMillis().
     */
    int64_t resultsTimeMillis(const MolGraph& target) const {
        return search_detail::resultsTimeMillis(search_detail::searchTimeLimitMillis(options_, query_->n, target.n, true));
    }
    /** The cancellation token of the options, or null; once it is cancelled, building mappings stops as well. */
    std::shared_ptr<const CancellationToken> cancellation() const { return options_.cancellation; }
    /** Heuristic confidence for selected tautomeric relaxations, separate from the MCS score. */
    static double computeTautomerConfidence(const MolGraph& query, const MolGraph& target,
                                            const Match& match, const ChemOptions& chemistry) {
        return detail::selectedTautomerConfidence(query, target, match.atomMapping, match.bondMapping,
                                                  chemistry);
    }
    /**
     * The lexicographically smallest mapping equivalent to `mapping` under the automorphisms of
     * both graphs. Out-of-range atoms throw std::out_of_range; incomplete generators and the time
     * and orbit limits throw std::length_error and std::runtime_error, as the Java form throws
     * IllegalArgumentException and IllegalStateException with the same messages.
     */
    static std::map<int,int> canonicalizeMapping(
        const MolGraph& g1, const MolGraph& g2,
        const std::map<int,int>& mapping) {
    if (mapping.empty()) return mapping;
    auto deadline = detail::steadyDeadline(1000);
    auto expired = [&] {
        return std::chrono::steady_clock::now() >= deadline;
    };

    for (const auto& [query, target] : mapping)
        if (query < 0 || query >= g1.n || target < 0 || target >= g2.n)
            throw std::out_of_range("Mapped atom index out of range");
    if (expired())
        throw std::runtime_error("Mapping canonicalisation exceeded its time limit");
    const auto& gens1 = g1.getAutomorphismGenerators();
    const auto& gens2 = g2.getAutomorphismGenerators();
    if (expired()) throw std::runtime_error("Mapping canonicalisation exceeded its time limit");
    if (g1.automorphismGeneratorsTruncated() || g2.automorphismGeneratorsTruncated())
        throw std::length_error("Mapping canonicalisation needs complete automorphism generators");
    if (gens1.empty() && gens2.empty()) return mapping;

    using PairVec = std::vector<std::pair<int,int>>;
    size_t maximumStates = std::min<size_t>(65536, 1000000 / mapping.size());
    std::vector<PairVec> frontier{PairVec(mapping.begin(), mapping.end())};
    std::set<PairVec> seen{frontier.front()};
    auto remember = [&](PairVec candidate) {
        std::sort(candidate.begin(), candidate.end());
        if (!seen.insert(candidate).second) return;
        if (seen.size() > maximumStates)
            throw std::length_error("Mapping canonicalisation exceeded its orbit state limit");
        frontier.push_back(std::move(candidate));
    };
    for (size_t head = 0; head < frontier.size(); ++head) {
        if (expired())
            throw std::runtime_error("Mapping canonicalisation exceeded its time limit");
        // Copy because appending to the frontier may reallocate its storage.
        PairVec current = frontier[head];
        for (const auto& generator : gens1) {
            if (expired()) throw std::runtime_error("Mapping canonicalisation exceeded its time limit");
            PairVec candidate;
            candidate.reserve(current.size());
            for (const auto& [query, target] : current)
                candidate.emplace_back(generator[query], target);
            remember(std::move(candidate));
        }
        for (const auto& generator : gens2) {
            if (expired()) throw std::runtime_error("Mapping canonicalisation exceeded its time limit");
            PairVec candidate;
            candidate.reserve(current.size());
            for (const auto& [query, target] : current)
                candidate.emplace_back(query, generator[target]);
            remember(std::move(candidate));
        }
    }
    std::map<int,int> result;
    for (const auto& [query, target] : *seen.begin()) result.emplace(query, target);
    return result;
}


    SearchResult match(const MolGraph& target) const {
        return search_detail::commonSubgraph(*query_, target, options_, 1, false);
    }
    SearchResult matchAll(const MolGraph& target, size_t limit = 1000) const {
        return search_detail::commonSubgraph(*query_, target, options_, limit, true);
    }
    SearchResult match(const std::string& target) const { return match(parseSMILES(target)); }
    /** Export under this prepared query's chemistry; unsupported export rules throw. */
    SmartsResult matchWithSmarts(const MolGraph& target) const {
        auto result = match(target);
        auto pattern = result.matches.empty() ? std::string{}
            : result.matches.front().toSmarts(*query_, options_.chemistry);
        return {std::move(result), std::move(pattern)};
    }
    SmartsResult matchWithSmarts(const std::string& target) const {
        return matchWithSmarts(parseSMILES(target));
    }
    template<class Targets>
    std::vector<SearchResult> matchMany(const Targets& targets, int numThreads = 0) const {
        return search_detail::matchMany(*this, targets, numThreads);
    }
    /** Atom-count threshold decisions, retaining exact CPU results for every searched slot. */
    std::vector<CorpusSearchResult> matchCorpus(const batch::TargetCorpus& corpus, int minimumAtoms,
        bool gpuPrefilter = false, int numThreads = 0) const {
        if (minimumAtoms < 1) throw std::invalid_argument("Minimum atom count must be positive");
        if (numThreads < 0) throw std::invalid_argument("Thread count must be non-negative");
        const bool atomObjective = options_.objective == Objective::MCIS && options_.atomWeights.empty();
        const bool filter = atomObjective && options_.chemistry.matchAtomType && options_.nodeLimit == 0
            && options_.timeoutMs != 0 && !options_.cancellation && options_.excludedTargetAtoms.empty();
        auto inputs = batch::detail::CorpusAccess::atomBoundInputs(
            corpus, filter ? query_.get() : nullptr, numThreads, gpuPrefilter);
        return batch::detail::matchCorpusInputs(*this, inputs.first, inputs.second, std::nullopt,
            minimumAtoms, false, atomObjective, numThreads);
    }
    SearchResult matchAll(const std::string& target, size_t limit = 1000) const {
        return matchAll(parseSMILES(target), limit);
    }

    /** Validate the selected graph, including the objective's topology rules. */
    static std::vector<std::string> validateMapping(const MolGraph& query, const MolGraph& target,
                                                   const Match& match, const Options& options = Options{}) {
        search_detail::validateSearchGraphStorage(query);
        search_detail::validateSearchGraphStorage(target);
        search_detail::validateOptions(options, query.n);
        for (int atom : options.excludedTargetAtoms)
            if (atom < 0 || atom >= target.n) throw std::out_of_range("Excluded target atom index out of range");
        search_detail::prepare(&query, target, options.chemistry);
        std::vector<std::string> errors;
        std::vector<int> mapping(query.n, -1);
        std::set<int> usedTargets, excluded(options.excludedTargetAtoms.begin(), options.excludedTargetAtoms.end());
        for (const auto& atom : match.atomMapping) {
            if (atom.first < 0 || atom.first >= query.n || atom.second < 0 || atom.second >= target.n) {
                errors.push_back("Mapped atom index out of range");
                continue;
            }
            mapping[atom.first] = atom.second;
            if (!usedTargets.insert(atom.second).second) errors.push_back("Target atom occurs more than once");
            if (excluded.count(atom.second)) errors.push_back("Mapped target atom is excluded");
            if (!search_detail::atomCompatible(query, atom.first, target, atom.second, options.chemistry))
                errors.push_back("Mapped atoms are incompatible");
        }
        if (!errors.empty()) return errors;
        std::set<int64_t> selected;
        std::set<int> incident;
        for (const auto& bond : match.bondMapping) {
            if (bond[0] < 0 || bond[1] < 0 || bond[0] >= query.n || bond[1] >= query.n
                || mapping[bond[0]] != bond[2] || mapping[bond[1]] != bond[3]
                || mapping[bond[0]] < 0 || mapping[bond[1]] < 0) {
                errors.push_back("Selected bond is outside the atom mapping");
                continue;
            }
            if (!selected.insert(MolGraph::bondKey(bond[0], bond[1])).second)
                errors.push_back("Selected bond occurs more than once");
            incident.insert(bond[0]); incident.insert(bond[1]);
            if (!search_detail::bondCompatible(query, bond[0], bond[1], target, bond[2], bond[3], options.chemistry))
                errors.push_back("Selected bonds are incompatible");
        }
        if (!errors.empty()) return errors;
        if (options.objective == Objective::MCIS) {
            for (auto first = match.atomMapping.begin(); first != match.atomMapping.end(); ++first)
                for (auto second = std::next(first); second != match.atomMapping.end(); ++second) {
                    const bool sourceBond = query.hasBond(first->first, second->first);
                    if (sourceBond != target.hasBond(first->second, second->second))
                        errors.push_back("Mapped topology is not induced");
                    if (sourceBond != (selected.count(MolGraph::bondKey(first->first, second->first)) != 0))
                        errors.push_back("MCIS must select every induced bond");
                }
        } else if (incident.size() != match.atomMapping.size()) errors.push_back("MCES cannot retain isolated atoms");
        const Match normalized = search_detail::makeMatch(mapping, match.bondMapping);
        if (options.connectivity == Connectivity::Connected && normalized.fragments.size() > 1)
            errors.push_back("Selected graph is disconnected");
        if (normalized.fragments.size() > static_cast<size_t>(options.maxFragments)) errors.push_back("Too many fragments");
        for (const auto& fragment : normalized.fragments)
            if (fragment.size() < static_cast<size_t>(options.minFragmentAtoms)) errors.push_back("Fragment is too small");
        if (!search_detail::selectedStereoCompatible(query, target, options.chemistry, normalized))
            errors.push_back("Selected stereochemistry is incompatible");
        if (options.chemistry.completeRingsOnly) {
            search_detail::RingConstraint queryRings(query, true), targetRings(target, true);
            if (!search_detail::completeRings(normalized, queryRings, targetRings)) errors.push_back("Selected ring is incomplete");
        }
        return errors;
    }
    static Match mappingMatch(const MolGraph& query, const MolGraph& target, const std::map<int, int>& atoms,
                              const Options& options = Options{}) {
        search_detail::validateSearchGraphStorage(query);
        search_detail::validateSearchGraphStorage(target);
        std::vector<int> mapping(query.n, -1);
        for (const auto& atom : atoms) {
            if (atom.first < 0 || atom.first >= query.n || atom.second < 0 || atom.second >= target.n)
                throw std::out_of_range("Mapped atom index out of range");
            mapping[atom.first] = atom.second;
        }
        return search_detail::makeMatch(mapping, search_detail::mappedBonds(search_detail::QueryView(query), target, mapping, options.chemistry));
    }
    static std::vector<std::string> validateMapping(const MolGraph& query, const MolGraph& target,
                                                   const std::map<int, int>& mapping, const Options& options = Options{}) {
        for (const auto& atom : mapping)
            if (atom.first < 0 || atom.first >= query.n || atom.second < 0 || atom.second >= target.n)
                return {"Mapped atom index out of range"};
        return validateMapping(query, target, mappingMatch(query, target, mapping, options), options);
    }
    /** True if no valid strict extension retains the correspondence; throws when proof is unavailable. */
    static bool isMappingMaximal(const MolGraph& query, const MolGraph& target, const Match& match,
                                 const Options& options = Options{}) {
        if (!validateMapping(query, target, match, options).empty()) return false;
        auto unweighted = options;
        unweighted.atomWeights.clear();
        const auto extension = search_detail::commonSubgraph(query, target, unweighted, 1, false, &match);
        if (!extension.optimal) throw std::runtime_error("Mapping maximality search stopped before proof");
        if (extension.matches.empty()) return match.atomMapping.empty();
        const auto& best = extension.matches.front();
        return best.atomMapping.size() == match.atomMapping.size() && best.bondMapping.size() == match.bondMapping.size();
    }
    static bool isMappingMaximal(const MolGraph& query, const MolGraph& target, const std::map<int, int>& mapping,
                                 const Options& options = Options{}) {
        if (!validateMapping(query, target, mapping, options).empty()) return false;
        return isMappingMaximal(query, target, mappingMatch(query, target, mapping, options), options);
    }
    /** Sequential non-overlapping queries retain target indices and chemistry. */
    static std::vector<SearchResult> matchConstrained(const std::vector<MolGraph>& queries, const MolGraph& target,
                                                       const Options& options = Options{}) {
        Options remaining = options;
        std::set<int> excluded(options.excludedTargetAtoms.begin(), options.excludedTargetAtoms.end());
        std::vector<SearchResult> results;
        results.reserve(queries.size());
        for (const auto& query : queries) {
            remaining.excludedTargetAtoms.assign(excluded.begin(), excluded.end());
            auto result = search_detail::commonSubgraph(query, target, remaining, 1, false);
            if (!result.matches.empty())
                for (const auto& atom : result.matches.front().atomMapping) excluded.insert(atom.second);
            results.push_back(std::move(result));
        }
        return results;
    }
    /**
     * Progressive common-graph reduction under one shared budget.
     * The result always has optimal false, including empty and single-input calls.
     * complete reports whether the progressive pass finished. Neither complete nor Completed
     * proves a global optimum or the absence of a qualifying common graph.
     */
    static ConsensusResult consensus(const std::vector<MolGraph>& molecules, const Options& options = Options{},
                                      double threshold = 1.0) {
        if (!std::isfinite(threshold) || threshold <= 0.0 || threshold > 1.0)
            throw std::invalid_argument("Consensus threshold must be in (0, 1]");
        if (!options.atomWeights.empty()) throw std::invalid_argument("Consensus does not accept atom weights");
        if (!options.excludedTargetAtoms.empty()) throw std::invalid_argument("Consensus has no single target for exclusions");
        ConsensusResult result;
        int maximum = 0;
        for (const auto& molecule : molecules) maximum = std::max(maximum, molecule.n);
        search_detail::Budget budget(options, maximum, maximum, true);
        const int requiredSupport = static_cast<int>(std::ceil(threshold * molecules.size()));
        search_detail::Score best;
        SubstructureOptions embedding;
        static_cast<SearchOptions&>(embedding) = options;
        embedding.induced = options.objective == Objective::MCIS;
        auto consider = [&](const MolGraph& graph, const Match& fragment, int source) {
            if (fragment.atomMapping.empty() || !validateMapping(molecules[source], molecules[source], fragment, options).empty()) return;
            std::map<int, Match> support;
            for (size_t molecule = 0; molecule < molecules.size() && !budget.stopped(); ++molecule) {
                const auto match = search_detail::substructure(&graph, nullptr, molecules[molecule], embedding, 1, false, &budget);
                if (!match.matches.empty()) support.emplace(static_cast<int>(molecule), match.matches.front());
            }
            if (static_cast<int>(support.size()) < requiredSupport) return;
            const auto score = search_detail::scoreMatch(fragment, options);
            if (result.found && !(best < score)) return;
            best = score;
            result.fragment = fragment; result.sourceIndex = source; result.graph = graph;
            result.smiles_ = fragment.toSmiles(molecules[source]);
            result.mappings = std::move(support); result.supportCount = static_cast<int>(result.mappings.size()); result.found = true;
        };
        for (size_t seed = 0; seed < molecules.size() && !budget.stopped(); ++seed) {
            MolGraph current = *search_detail::ownQuery(molecules[seed]);
            std::vector<int> original(current.n);
            std::iota(original.begin(), original.end(), 0);
            auto sourceFragment = [&] {
                std::vector<int> mapping(molecules[seed].n, -1);
                for (int atom : original) mapping[atom] = atom;
                std::vector<std::array<int, 4>> bonds;
                for (int atom = 0; atom < current.n; ++atom)
                    for (int neighbor : current.neighbors[atom])
                        if (atom < neighbor) bonds.push_back({original[atom], original[neighbor], original[atom], original[neighbor]});
                return search_detail::makeMatch(mapping, std::move(bonds));
            };
            consider(current, sourceFragment(), static_cast<int>(seed));
            for (size_t molecule = 0; molecule < molecules.size() && !budget.stopped(); ++molecule) {
                if (molecule == seed) continue;
                auto pair = search_detail::commonSubgraph(current, molecules[molecule], options, 1, false, nullptr, &budget);
                if (pair.matches.empty()) continue;
                const auto& selected = pair.matches.front();
                std::vector<int> retained;
                for (const auto& atom : selected.atomMapping) retained.push_back(original[atom.first]);
                current = selected.toGraph(current);
                original = std::move(retained);
                consider(current, sourceFragment(), static_cast<int>(seed));
            }
        }
        result.status = budget.status();
        if (result.status == SearchStatus::Completed && search_detail::needsRingProof(options.chemistry))
            for (const auto& molecule : molecules)
                if (search_detail::ringProofTruncated(molecule, options.chemistry)) {
                    result.status = SearchStatus::LimitReached;
                    break;
                }
        result.complete = result.status == SearchStatus::Completed;
        result.nodesVisited = budget.nodesVisited(); result.elapsedMillis = budget.elapsedMillis();
        return result;
    }
};

} // namespace smsd

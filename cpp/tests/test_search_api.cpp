/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"

#include <cassert>
#include <iostream>
#include <random>
#include <set>
#include <thread>

static smsd::MolGraph graph(int atoms, unsigned edges, bool heteroatoms = false) {
    std::vector<int> numbers(atoms, 6);
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    unsigned bit = 0;
    for (int first = 0; first < atoms; ++first) {
        if (heteroatoms && first % 3 == 0) numbers[first] = 7;
        for (int second = first + 1; second < atoms; ++second, ++bit) {
            if ((edges & (1u << bit)) == 0) continue;
            neighbors[first].push_back(second); neighbors[second].push_back(first);
            orders[first].push_back(1); orders[second].push_back(1);
        }
    }
    return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(numbers)
        .setNeighbors(neighbors).setBondOrders(orders).build(false);
}

struct Expected {
    int64_t primary = 0;
    int atoms = 0, bonds = 0;
    bool operator<(const Expected& other) const {
        return std::tie(primary, atoms, bonds) < std::tie(other.primary, other.atoms, other.bonds);
    }
};

static Expected oracle(const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& options) {
    std::vector<int> mapping(query.n, -1);
    std::vector<bool> used(target.n, false);
    Expected best;
    std::function<void(int)> visit = [&](int queryAtom) {
        if (queryAtom == query.n) {
            std::set<int> atoms;
            std::vector<std::pair<int, int>> edges;
            bool valid = true;
            for (int first = 0; first < query.n; ++first) {
                if (mapping[first] < 0) continue;
                if (options.objective == smsd::MCS::Objective::MCIS) atoms.insert(first);
                for (int second = first + 1; second < query.n; ++second) {
                    if (mapping[second] < 0) continue;
                    const bool source = query.hasBond(first, second), destination = target.hasBond(mapping[first], mapping[second]);
                    if (options.objective == smsd::MCS::Objective::MCIS && source != destination) valid = false;
                    if (source && destination) {
                        edges.emplace_back(first, second);
                        atoms.insert(first); atoms.insert(second);
                    }
                }
            }
            if (!valid || atoms.empty()) return;
            std::set<int> unseen = atoms;
            int fragments = 0;
            while (!unseen.empty()) {
                std::vector<int> component{*unseen.begin()};
                unseen.erase(component.front());
                for (size_t position = 0; position < component.size(); ++position)
                    for (const auto& edge : edges) {
                        int adjacent = edge.first == component[position] ? edge.second
                            : edge.second == component[position] ? edge.first : -1;
                        if (adjacent >= 0 && unseen.erase(adjacent)) component.push_back(adjacent);
                    }
                if (component.size() < static_cast<size_t>(options.minFragmentAtoms)) return;
                ++fragments;
            }
            if ((options.connectivity == smsd::MCS::Connectivity::Connected && fragments > 1)
                || fragments > options.maxFragments) return;
            Expected score;
            score.atoms = static_cast<int>(atoms.size());
            if (options.objective == smsd::MCS::Objective::MCES) score.primary = edges.size();
            else {
                double weight = 0.0;
                for (int atom : atoms) weight += options.atomWeights.empty() ? 1.0 : options.atomWeights[atom];
                score.primary = options.atomWeights.empty() ? int64_t(atoms.size()) : static_cast<int64_t>(weight * 1000.0);
                score.bonds = static_cast<int>(edges.size());
            }
            if (best < score) best = score;
            return;
        }
        visit(queryAtom + 1);
        for (int targetAtom = 0; targetAtom < target.n; ++targetAtom) {
            if (used[targetAtom] || query.atomicNum[queryAtom] != target.atomicNum[targetAtom]) continue;
            mapping[queryAtom] = targetAtom; used[targetAtom] = true;
            visit(queryAtom + 1);
            mapping[queryAtom] = -1; used[targetAtom] = false;
        }
    };
    visit(0);
    return best;
}

static void testObjectives() {
    const auto induced = smsd::MCS::find("C1CC1").match("CCC");
    assert(induced.complete && induced.optimal && induced.score == 2.0);
    assert(induced.matches.front().bondMapping.size() == 1);
    smsd::MCS::Options options;
    options.objective = smsd::MCS::Objective::MCES;
    const auto edges = smsd::MCS::find("C1CC1").withOptions(options).match("CCC");
    assert(edges.complete && edges.optimal && edges.score == 2.0);
    assert(edges.matches.front().atomMapping.size() == 3 && edges.matches.front().bondMapping.size() == 2);
    const auto fragment = smsd::parseSMILES(edges.matches.front().toSmiles(smsd::parseSMILES("C1CC1")));
    assert(fragment.n == 3 && fragment.neighbors[0].size() + fragment.neighbors[1].size() + fragment.neighbors[2].size() == 4);
    options.connectivity = smsd::MCS::Connectivity::Fragments;
    assert(smsd::MCS::find("C.C").withOptions(options).match("CC").matches.empty());
    options.objective = smsd::MCS::Objective::MCIS;
    assert(smsd::MCS::find("C.C").withOptions(options).match("C.C").score == 2.0);
    options.atomWeights = {1, 1, 1, 10};
    assert(smsd::MCS::find("CCC.N").withOptions(options).match("CCCN").score == 12.0);
    options.objective = smsd::MCS::Objective::MCES;
    bool rejected = false;
    try { smsd::MCS::find("CCC.N").withOptions(options).match("CCCN"); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}

static void testOracle() {
    std::mt19937 random(719);
    for (int trial = 0; trial < 96; ++trial) {
        const auto query = graph(4, random() & 63u, trial % 3 == 0);
        const auto target = graph(5, random() & 1023u, trial % 3 == 0);
        for (auto objective : {smsd::MCS::Objective::MCIS, smsd::MCS::Objective::MCES}) {
            for (auto connectivity : {smsd::MCS::Connectivity::Connected, smsd::MCS::Connectivity::Fragments}) {
                smsd::MCS::Options options;
                options.objective = objective; options.connectivity = connectivity;
                if (objective == smsd::MCS::Objective::MCIS && trial % 2 == 0) options.atomWeights = {2, -3, 1, 4};
                if (trial % 4 == 0) options.minFragmentAtoms = 2;
                if (trial % 5 == 0) options.maxFragments = 1;
                const auto expected = oracle(query, target, options);
                const auto actual = smsd::MCS::find(query).withOptions(options).match(target);
                assert(actual.complete && actual.optimal && actual.status == smsd::SearchStatus::Completed);
                assert(actual.score == (options.atomWeights.empty() ? double(expected.primary) : expected.primary / 1000.0));
                if (actual.matches.empty()) assert(expected.atoms == 0);
                else {
                    assert(actual.matches.front().atomMapping.size() == static_cast<size_t>(expected.atoms));
                    if (objective == smsd::MCS::Objective::MCIS)
                        assert(actual.matches.front().bondMapping.size() == static_cast<size_t>(expected.bonds));
                }
            }
        }
    }
}

static void testResultsAndStops() {
    const auto query = smsd::Substructure::find("C");
    const auto limited = query.matchAll("CC", 1);
    assert(limited.matches.size() == 1 && limited.optimal && !limited.complete
           && limited.status == smsd::SearchStatus::LimitReached);
    const auto complete = query.matchAll("CC", 2);
    assert(complete.matches.size() == 2 && complete.complete && complete.optimal);
    assert(query.withTimeout(0).match("CC").status == smsd::SearchStatus::TimedOut);
    smsd::MCS::Options options;
    options.nodeLimit = 1;
    const auto stopped = smsd::MCS::find("CCCC").withOptions(options).match("CCCN");
    assert(stopped.status == smsd::SearchStatus::LimitReached && !stopped.optimal && !stopped.complete);
    assert(!stopped.matches.empty() && stopped.nodesVisited == 1);
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    assert(smsd::MCS::find("CCCO").withOptions(options).match("CCCN").status == smsd::SearchStatus::Cancelled);
    assert(smsd::Substructure::find(smsd::MolGraph{}).withTimeout(0).match("C").status == smsd::SearchStatus::TimedOut);
}

static void testRingsAndStereo() {
    smsd::MCS::Options options;
    options.chemistry.completeRingsOnly = true;
    const auto ring = smsd::MCS::find("c1ccccc1").withOptions(options).match("c1ccc2ccccc2c1");
    assert(ring.optimal && ring.score == 6.0);
    assert(smsd::MCS::find("CC").withOptions(options).match("C1CCCCC1").matches.empty());
    smsd::Substructure::Options sub;
    sub.chemistry.useChirality = true; sub.chemistry.useBondStereo = true;
    assert(smsd::Substructure::find("C[C@H](F)Cl").withOptions(sub).match("C[C@@H](F)Cl").matches.empty());
    assert(smsd::Substructure::find("C/C=C/C").withOptions(sub).match("C/C=C\\C").matches.empty());
    assert(smsd::Substructure::smarts("F/C=C/F").match("F/C=C\\F").matches.empty());
}

static void testConcurrentReuse() {
    const auto query = smsd::Substructure::smarts("[C;$(CO)]");
    const auto maximum = smsd::MCS::find("CCO");
    const auto good = smsd::parseSMILES("CCO"), bad = smsd::parseSMILES("CCCC");
    std::vector<std::thread> workers;
    for (int worker = 0; worker < 6; ++worker) workers.emplace_back([&] {
        for (int iteration = 0; iteration < 16; ++iteration) {
            assert(!query.match(good).matches.empty() && query.match(bad).matches.empty());
            assert(maximum.match(good).score == 3.0 && maximum.match(bad).score == 2.0);
        }
    });
    for (auto& worker : workers) worker.join();
}

static void testGraphAndWorkflows() {
    const auto source = smsd::parseSMILES("[13CH3:7]c1ccccc1");
    smsd::MCS::Options edges;
    edges.objective = smsd::MCS::Objective::MCES;
    const auto result = smsd::MCS::find(source).withOptions(edges).match("CCC");
    assert(result.optimal && !result.matches.empty());
    const auto selected = result.matches.front().toGraph(source);
    assert(selected.n == 3);
    int position = 0;
    for (const auto& atom : result.matches.front().atomMapping) {
        assert(selected.hydrogenCount[position] == source.hydrogenCount[atom.first]);
        assert(selected.aromatic[position] == source.aromatic[atom.first]);
        assert(selected.ring[position] == source.ring[atom.first]);
        assert(selected.massNumber[position] == source.massNumber[atom.first]);
        ++position;
    }
    const auto triangle = smsd::parseSMILES("C1CC1"), path = smsd::parseSMILES("CCC");
    const auto all = smsd::MCS::find(triangle).withOptions(edges).matchAll(path);
    assert(all.complete && all.matches.size() == 6);
    for (const auto& match : all.matches) {
        assert(smsd::MCS::validateMapping(triangle, path, match, edges).empty());
        assert(smsd::MCS::isMappingMaximal(triangle, path, match, edges));
    }
    auto invalid = all.matches.front();
    invalid.bondMapping.push_back({0, 2, 0, 2});
    assert(!smsd::MCS::validateMapping(triangle, path, invalid, edges).empty());
    edges.excludedTargetAtoms = {1};
    assert(smsd::MCS::find(triangle).withOptions(edges).match(path).matches.empty());
    edges.objective = smsd::MCS::Objective::MCIS;
    assert(smsd::MCS::find(triangle).withOptions(edges).match(path).score == 1.0);
    edges.excludedTargetAtoms = {3};
    bool rejected = false;
    try { smsd::MCS::find(triangle).withOptions(edges).match(path); }
    catch (const std::out_of_range&) { rejected = true; }
    assert(rejected);
    const std::vector<smsd::MolGraph> queries{smsd::parseSMILES("CC"), smsd::parseSMILES("CC")};
    const auto constrained = smsd::MCS::matchConstrained(queries, smsd::parseSMILES("CCCC"));
    assert(constrained.size() == 2 && constrained[0].matches.front().atomMapping.size() == 2
           && constrained[1].matches.front().atomMapping.size() == 2);
    std::set<int> claimed;
    for (const auto& match : constrained)
        for (const auto& atom : match.matches.front().atomMapping) assert(claimed.insert(atom.second).second);
    const std::vector<smsd::MolGraph> molecules{smsd::parseSMILES("NN"), smsd::parseSMILES("CCO"), smsd::parseSMILES("CCN")};
    const auto consensus = smsd::MCS::consensus(molecules, {}, 2.0 / 3.0);
    assert(consensus.found && consensus.supportCount >= 2 && consensus.graph.n == 2 && consensus.complete && !consensus.optimal);
    assert(!consensus.toSmiles().empty());
    for (const auto& mapping : consensus.mappings)
        assert(smsd::MCS::validateMapping(consensus.graph, molecules[mapping.first], mapping.second).empty());
    const auto decomposition = smsd::Substructure::decomposeRGroups(smsd::parseSMILES("c1ccccc1"), {smsd::parseSMILES("Cc1ccccc1")});
    assert(decomposition.size() == 1 && decomposition.front().core.n == 6 && decomposition.front().rgroups.size() == 1);
    const auto batch = smsd::Substructure::find("CC").matchMany(queries, 2);
    assert(batch.size() == 2 && !batch.front().matches.empty());
}

static void testCompleteQueryBound() {
    const std::vector<std::pair<int, int>> targetEdges{{0,1},{0,2},{0,3},{0,4},{0,6},{0,7},{0,8},
        {1,2},{1,3},{1,5},{1,6},{1,7},{1,8},{2,3},{2,4},{2,5},{2,6},{2,7},{2,8},
        {3,4},{3,5},{3,6},{3,7},{3,8},{4,6},{4,8},{5,6},{5,7},{5,8},{6,7},{6,8},{7,8}};
    std::vector<std::vector<int>> adjacent(9), complete(9);
    for (const auto& edge : targetEdges) { adjacent[edge.first].push_back(edge.second); adjacent[edge.second].push_back(edge.first); }
    for (int first = 0; first < 9; ++first) for (int second = 0; second < 9; ++second)
        if (first != second) complete[first].push_back(second);
    const auto query = smsd::MolGraph::Builder().atomCount(9).atomicNumbers(std::vector<int>(9,6)).setNeighbors(complete).build(false);
    const auto target = smsd::MolGraph::Builder().atomCount(9).atomicNumbers(std::vector<int>(9,6)).setNeighbors(adjacent).build(false);
    const auto result = smsd::MCS::find(query).withTimeout(1000).matchAll(target, 1);
    assert(result.score == 7.0 && result.optimal && !result.complete && result.status == smsd::SearchStatus::LimitReached);
    assert(result.matches.size() == 1 && result.matches.front().bondMapping.size() == 21);
}

static void testReviewChemistryWorkflows() {
    smsd::MCSOptions ringOnly; ringOnly.chemistry.ringMatchesRingOnly = true;
    assert(smsd::MCS::find("c1ccccc1CCCC").withOptions(ringOnly).match("CCCCCC").score == 4.0);
    assert(smsd::MCS::find("CCCCCC").withOptions(ringOnly).match("c1ccccc1CCCC").score == 4.0);
    smsd::SubstructureOptions fusion; fusion.chemistry.ringFusionMode = smsd::ChemOptions::RingFusionMode::STRICT;
    const auto benzeneQuery = smsd::Substructure::find("c1ccccc1").withOptions(fusion);
    assert(benzeneQuery.match("c1ccc2ccccc2c1").matches.empty());
    assert(!benzeneQuery.match("Cc1ccccc1").matches.empty());
    const auto acid = smsd::parseSMILES("OC(=O)C");
    assert(!acid.automorphismGeneratorsTruncated());
    const auto canonical = smsd::MCS::canonicalizeMapping(acid, acid, {{0, 0}});
    assert(smsd::MCS::canonicalizeMapping(acid, acid, canonical) == canonical);
    const auto mixture = smsd::parseSMILES("F/C=C/F.[NH4+].[13CH3:7]O");
    const auto components = smsd::search_detail::splitComponents(mixture);
    assert(components.size() == 3);
    assert(components[0].hydrogenCount == std::vector<int>({0, 1, 1, 0}));
    assert(components[1].hydrogenCount == std::vector<int>({4}));
    assert(components[2].massNumber[0] == 13 && components[2].atomClass[0] == 7 && components[2].hydrogenCount[0] == 3);
    smsd::SubstructureOptions stereo; stereo.chemistry.useBondStereo = true;
    assert(smsd::Substructure::find("F/C=C\\F").withOptions(stereo).match(components[0]).matches.empty());
    const std::vector<smsd::MolGraph> molecules{smsd::parseSMILES("c1ccncc1"), smsd::parseSMILES("c1ccccc1"), smsd::parseSMILES("Cc1ccccc1")};
    for (bool ringOnly : {false, true}) {
        smsd::MCSOptions options;
        options.chemistry.ringMatchesRingOnly = ringOnly;
        options.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
        const auto consensus = smsd::MCS::consensus(molecules, options);
        assert(consensus.found && consensus.graph.n == 5 && consensus.supportCount == 3);
        const std::vector<smsd::MolGraph> queries{smsd::parseSMILES("c1ccccc1"), smsd::parseSMILES("CCCC")};
        const auto constrained = smsd::MCS::matchConstrained(queries, smsd::parseSMILES("c1ccc2ccccc2c1"), options);
        assert(constrained.front().matches.front().atomMapping.size() == 6 && constrained.back().matches.empty());
    }
    smsd::MCSOptions rings; rings.chemistry.completeRingsOnly = true;
    assert(smsd::MCS::find("OC1CCCCC1N").withOptions(rings).match("O=C1NC(=O)c2ccccc12").score == 6.0);
    const auto source = smsd::parseSMILES("c1ccccc1CCCO"), target = smsd::parseSMILES("c1ccccc1NCCO");
    smsd::MCSOptions strict; strict.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    const auto maximum = smsd::MCS::find(source).withOptions(strict).match(target);
    assert(maximum.optimal && maximum.matches.front().atomMapping.size() == 6);
    assert(smsd::MCS::isMappingMaximal(source, target, maximum.matches.front(), strict));
}

int main() {
    testObjectives(); testOracle(); testResultsAndStops(); testRingsAndStereo(); testConcurrentReuse(); testGraphAndWorkflows(); testCompleteQueryBound(); testReviewChemistryWorkflows();
    std::cout << "PASS objectives, exhaustive oracle, results/stops, rings/stereo, concurrent reuse\n";
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/ring_finder.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

smsd::MCSOptions options(smsd::MCSObjective objective, smsd::MCSConnectivity connectivity) {
    smsd::MCSOptions settings;
    settings.objective = objective;
    settings.connectivity = connectivity;
    settings.chemistry.completeRingsOnly = true;
    settings.timeoutMs = 10000;
    return settings;
}

void completed(const smsd::SearchResult& result) {
    require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
            "expected a complete, proved search");
}

std::string key(const std::map<int, int>& mapping) {
    std::string result;
    for (const auto& atom : mapping)
        result += std::to_string(atom.first) + ":" + std::to_string(atom.second) + ",";
    return result;
}

/** Small single-bond carbon graphs with their explicitly listed relevant cycles. */
struct Fixture {
    int atoms;
    std::vector<std::vector<int>> cycles;
    std::vector<std::vector<bool>> edges, ringEdges;
    std::vector<bool> ringAtoms;

    Fixture(int count, const std::vector<std::array<int, 2>>& bonds,
            std::vector<std::vector<int>> rings)
        : atoms(count), cycles(std::move(rings)),
          edges(count, std::vector<bool>(count)), ringEdges(count, std::vector<bool>(count)),
          ringAtoms(count) {
        for (const auto& bond : bonds) edges[bond[0]][bond[1]] = edges[bond[1]][bond[0]] = true;
        for (const auto& cycle : cycles)
            for (size_t index = 0; index < cycle.size(); ++index) {
                const int first = cycle[index], second = cycle[(index + 1) % cycle.size()];
                ringAtoms[first] = true;
                ringEdges[first][second] = ringEdges[second][first] = true;
            }
    }

    /** Flags: absent, all topological ring flags, or alternating topological ring flags. */
    smsd::MolGraph graph(int flags) const {
        std::vector<std::vector<int>> neighbours(atoms), orders(atoms);
        std::vector<uint8_t> atomFlags(atoms);
        std::vector<std::vector<bool>> bondFlags(atoms, std::vector<bool>(atoms));
        for (int atom = 0; atom < atoms; ++atom) {
            for (int other = 0; other < atoms; ++other) if (edges[atom][other]) {
                neighbours[atom].push_back(other);
                orders[atom].push_back(1);
                bondFlags[atom][other] = flags != 0 && ringEdges[atom][other]
                    && (flags == 1 || ((atom + other) & 1) == 0);
            }
            atomFlags[atom] = flags != 0 && ringAtoms[atom] && (flags == 1 || (atom & 1) == 0);
        }
        return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(std::vector<int>(atoms, 6))
            .setNeighbors(std::move(neighbours)).setBondOrders(std::move(orders))
            .ringFlags(std::move(atomFlags)).bondRingFlags(std::move(bondFlags)).build(false);
    }

    bool complete(const std::vector<bool>& selected) const {
        std::vector<bool> coveredAtoms(atoms);
        std::vector<std::vector<bool>> coveredEdges(atoms, std::vector<bool>(atoms));
        for (const auto& cycle : cycles) {
            bool whole = true;
            for (int atom : cycle) whole = whole && selected[atom];
            if (!whole) continue;
            for (size_t index = 0; index < cycle.size(); ++index) {
                const int first = cycle[index], second = cycle[(index + 1) % cycle.size()];
                coveredAtoms[first] = true;
                coveredEdges[first][second] = coveredEdges[second][first] = true;
            }
        }
        for (int first = 0; first < atoms; ++first) if (selected[first]) {
            if (ringAtoms[first] && !coveredAtoms[first]) return false;
            for (int second = first + 1; second < atoms; ++second)
                if (selected[second] && ringEdges[first][second] && !coveredEdges[first][second])
                    return false;
        }
        return true;
    }
};

Fixture cycle(int atoms) {
    std::vector<std::array<int, 2>> edges;
    std::vector<int> ring;
    for (int atom = 0; atom < atoms; ++atom) {
        edges.push_back({atom, (atom + 1) % atoms});
        ring.push_back(atom);
    }
    return Fixture(atoms, edges, {ring});
}

std::pair<Fixture, Fixture> pair(const std::string& family) {
    if (family == "fused") return {cycle(3), Fixture(4,
        {{0, 1}, {1, 2}, {0, 2}, {1, 3}, {0, 3}}, {{0, 1, 2}, {0, 1, 3}})};
    if (family == "bridged") return {cycle(4), Fixture(5,
        {{0, 2}, {1, 2}, {0, 3}, {1, 3}, {0, 4}, {1, 4}},
        {{0, 2, 1, 3}, {0, 2, 1, 4}, {0, 3, 1, 4}})};
    if (family == "spiro") return {cycle(3), Fixture(5,
        {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {3, 4}, {0, 4}}, {{0, 1, 2}, {0, 3, 4}})};
    if (family == "mixed") return {Fixture(4,
        {{0, 1}, {1, 2}, {0, 2}, {2, 3}}, {{0, 1, 2}}), Fixture(5,
        {{0, 1}, {1, 2}, {0, 2}, {2, 3}, {3, 4}}, {{0, 1, 2}})};
    if (family == "incomplete") return {cycle(3), Fixture(3, {{0, 1}, {1, 2}}, {})};
    throw std::invalid_argument(family);
}

struct Oracle {
    int atoms = 0, bonds = 0;
    std::set<std::string> keys;
};

void enumerate(const Fixture& query, const Fixture& target, smsd::MCSConnectivity connectivity,
               int queryAtom, std::vector<int>& mapping, std::vector<bool>& used, Oracle& result) {
    if (queryAtom < query.atoms) {
        enumerate(query, target, connectivity, queryAtom + 1, mapping, used, result);
        for (int targetAtom = 0; targetAtom < target.atoms; ++targetAtom) {
            if (used[targetAtom]) continue;
            bool induced = true;
            for (int previous = 0; previous < queryAtom; ++previous)
                if (mapping[previous] >= 0
                    && query.edges[queryAtom][previous] != target.edges[targetAtom][mapping[previous]])
                    induced = false;
            if (!induced) continue;
            mapping[queryAtom] = targetAtom;
            used[targetAtom] = true;
            enumerate(query, target, connectivity, queryAtom + 1, mapping, used, result);
            used[targetAtom] = false;
            mapping[queryAtom] = -1;
        }
        return;
    }
    std::vector<bool> selected(query.atoms);
    int atoms = 0, bonds = 0, firstSelected = -1;
    std::map<int, int> pairs;
    for (int atom = 0; atom < query.atoms; ++atom) if (mapping[atom] >= 0) {
        selected[atom] = true;
        pairs.emplace(atom, mapping[atom]);
        firstSelected = atom;
        ++atoms;
        for (int other = 0; other < atom; ++other)
            if (selected[other] && query.edges[atom][other]) ++bonds;
    }
    if (atoms == 0 || !query.complete(selected) || !target.complete(used)) return;
    if (connectivity == smsd::MCSConnectivity::Connected) {
        std::vector<bool> reached(query.atoms);
        reached[firstSelected] = true;
        for (int pass = 0; pass < query.atoms; ++pass)
            for (int first = 0; first < query.atoms; ++first) if (reached[first])
                for (int second = 0; second < query.atoms; ++second)
                    if (selected[second] && query.edges[first][second]) reached[second] = true;
        for (int atom = 0; atom < query.atoms; ++atom) if (selected[atom] && !reached[atom]) return;
    }
    if (atoms < result.atoms || (atoms == result.atoms && bonds < result.bonds)) return;
    if (atoms > result.atoms || bonds > result.bonds) {
        result.atoms = atoms;
        result.bonds = bonds;
        result.keys.clear();
    }
    result.keys.insert(key(pairs));
}

Oracle oracle(const Fixture& query, const Fixture& target, smsd::MCSConnectivity connectivity) {
    Oracle result;
    std::vector<int> mapping(query.atoms, -1);
    std::vector<bool> used(target.atoms);
    enumerate(query, target, connectivity, 0, mapping, used, result);
    return result;
}

void minimal(smsd::MCSObjective objective, smsd::MCSConnectivity connectivity, bool reverse) {
    const auto built = cycle(6).graph(0);
    const auto parsed = smsd::parseSMILES("C1CCCCC1");
    for (int atom = 0; atom < built.n; ++atom) {
        require(!built.ring[atom] && parsed.ring[atom], "six-cycle fixture flags");
        for (int other : built.neighbors[atom])
            require(!built.bondInRing(atom, other), "builder bond must be unflagged");
    }
    require(!smsd::ringsTruncated(built) && !smsd::ringsTruncated(parsed), "six-cycle perception");
    const auto& query = reverse ? parsed : built;
    const auto& target = reverse ? built : parsed;
    const auto settings = options(objective, connectivity);
    const auto result = smsd::MCS::find(query).withOptions(settings).match(target);
    completed(result);
    require(result.score == 6 && !result.matches.empty(), "whole six-cycle must match");
    const auto& match = result.matches.front();
    require(match.atomMapping.size() == 6 && match.bondMapping.size() == 6, "six atoms and bonds");
    require(smsd::MCS::validateMapping(query, target, match, settings).empty(), "six-cycle validity");
}

void explicitRingOnly(smsd::MCSObjective objective, smsd::MCSConnectivity connectivity,
                      bool reverse, bool complete) {
    const auto built = cycle(6).graph(0);
    const auto parsed = smsd::parseSMILES("C1CCCCC1");
    const auto& query = reverse ? parsed : built;
    const auto& target = reverse ? built : parsed;
    auto settings = options(objective, connectivity);
    settings.chemistry.completeRingsOnly = complete;
    settings.chemistry.ringMatchesRingOnly = true;
    const auto result = smsd::MCS::find(query).withOptions(settings).match(target);
    completed(result);
    require(result.matches.empty() && result.score == 0, "explicit ring-only raw flag mismatch");
}

void exhaustive(const std::string& family, smsd::MCSConnectivity connectivity, bool reverse) {
    const auto shapes = pair(family);
    const auto& queryShape = reverse ? shapes.second : shapes.first;
    const auto& targetShape = reverse ? shapes.first : shapes.second;
    const auto expected = oracle(queryShape, targetShape, connectivity);
    const auto settings = options(smsd::MCSObjective::MCIS, connectivity);
    for (int flags = 0; flags < 3; ++flags) {
        const auto query = queryShape.graph(flags == 0 ? 1 : 0);
        const auto target = targetShape.graph(flags == 0 ? 1 : flags == 1 ? 2 : 0);
        require(!smsd::ringsTruncated(query) && !smsd::ringsTruncated(target), "small rings truncated");
        const auto prepared = smsd::MCS::find(query).withOptions(settings);
        const auto all = prepared.matchAll(target, 1000);
        completed(all);
        require(all.score == expected.atoms, "all-mapping atom optimum");
        std::set<std::string> actual;
        for (const auto& match : all.matches) {
            require(match.atomMapping.size() == static_cast<size_t>(expected.atoms)
                        && match.bondMapping.size() == static_cast<size_t>(expected.bonds),
                    "all-mapping secondary bond optimum");
            require(smsd::MCS::validateMapping(query, target, match, settings).empty(), "MCIS validity");
            require(actual.insert(key(match.atomMapping)).second, "duplicate mapping");
        }
        require(actual == expected.keys, "every optimal MCIS mapping must be retained");
        const auto one = prepared.match(target);
        completed(one);
        require(one.score == expected.atoms && (one.matches.empty() == expected.keys.empty()),
                "single-mapping optimum and presence");
        if (!one.matches.empty()) {
            require(one.matches.front().bondMapping.size() == static_cast<size_t>(expected.bonds)
                        && expected.keys.count(key(one.matches.front().atomMapping)) != 0,
                    "single mapping must belong to oracle set");
            require(smsd::MCS::validateMapping(query, target, one.matches.front(), settings).empty(),
                    "single MCIS validity");
        }
    }
}

} // namespace

int main() {
    int failures = 0;
    for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
        for (bool reverse : {false, true}) {
            for (auto objective : {smsd::MCSObjective::MCES, smsd::MCSObjective::MCIS}) {
                try { minimal(objective, connectivity, reverse); }
                catch (const std::exception& error) {
                    std::cerr << "six-cycle " << static_cast<int>(objective) << ' '
                              << static_cast<int>(connectivity) << ' ' << reverse << ": " << error.what() << '\n';
                    ++failures;
                }
            }
            for (auto objective : {smsd::MCSObjective::MCIS, smsd::MCSObjective::MCES}) {
                for (bool complete : {false, true}) {
                    try { explicitRingOnly(objective, connectivity, reverse, complete); }
                    catch (const std::exception& error) {
                        std::cerr << "explicit ring-only " << static_cast<int>(objective) << ' '
                                  << static_cast<int>(connectivity) << ' ' << reverse << ' '
                                  << complete << ": " << error.what() << '\n';
                        ++failures;
                    }
                }
            }
            for (const std::string family : {"fused", "bridged", "spiro", "mixed", "incomplete"}) {
                try { exhaustive(family, connectivity, reverse); }
                catch (const std::exception& error) {
                    std::cerr << family << ' ' << static_cast<int>(connectivity) << ' ' << reverse
                              << ": " << error.what() << '\n';
                    ++failures;
                }
            }
        }
    }
    return failures ? 1 : 0;
}

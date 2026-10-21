/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <new>
#include <set>
#include <stdexcept>

namespace {
/** While tracking, the largest single allocation made. */
bool trackAllocations = false;
size_t largestAllocation = 0;
}

void* operator new(std::size_t size) {
    if (trackAllocations) largestAllocation = std::max(largestAllocation, size);
    if (void* memory = std::malloc(size == 0 ? 1 : size)) return memory;
    throw std::bad_alloc();
}
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void* memory) noexcept { std::free(memory); }
void operator delete[](void* memory) noexcept { std::free(memory); }
void operator delete(void* memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void* memory, std::size_t) noexcept { std::free(memory); }

using AtomMapping = std::map<int, int>;
using BondMapping = std::vector<std::array<int, 4>>;
using MatchKey = std::pair<AtomMapping, BondMapping>;
using Objective = std::pair<int, int>;

using smsd_test::require;

/** Enumerate injective atom maps and every compatible selected-edge subset. */
static std::pair<Objective, std::set<MatchKey>> selectedEdgeOracle(
        const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& options) {
    Objective best{0, 0};
    std::set<MatchKey> matches;
    AtomMapping mapping;
    std::set<int> used;
    const std::set<int> excluded(options.excludedTargetAtoms.begin(), options.excludedTargetAtoms.end());
    std::function<void(int)> enumerate = [&](int queryAtom) {
        if (queryAtom < query.n) {
            enumerate(queryAtom + 1);
            for (int targetAtom = 0; targetAtom < target.n; ++targetAtom) {
                if (used.count(targetAtom) || excluded.count(targetAtom)) continue;
                if (options.chemistry.matchAtomType && query.atomicNum[queryAtom] != target.atomicNum[targetAtom]) continue;
                mapping.emplace(queryAtom, targetAtom);
                used.insert(targetAtom);
                enumerate(queryAtom + 1);
                used.erase(targetAtom);
                mapping.erase(queryAtom);
            }
            return;
        }
        BondMapping compatible;
        for (const auto& first : mapping) for (const auto& second : mapping) {
            if (first.first >= second.first || !query.hasBond(first.first, second.first)
                || !target.hasBond(first.second, second.second)) continue;
            if (options.chemistry.matchBondOrder == smsd::ChemOptions::BondOrderMode::STRICT
                && query.bondOrder(first.first, second.first) != target.bondOrder(first.second, second.second)) continue;
            compatible.push_back({first.first, second.first, first.second, second.second});
        }
        for (size_t subset = 1; subset < (size_t{1} << compatible.size()); ++subset) {
            AtomMapping incident;
            BondMapping selected;
            std::vector<std::vector<int>> neighbors(query.n);
            for (size_t edge = 0; edge < compatible.size(); ++edge) if (subset & (size_t{1} << edge)) {
                const auto& bond = compatible[edge];
                selected.push_back(bond);
                incident.emplace(bond[0], bond[2]);
                incident.emplace(bond[1], bond[3]);
                neighbors[bond[0]].push_back(bond[1]);
                neighbors[bond[1]].push_back(bond[0]);
            }
            std::set<int> visited;
            int fragments = 0;
            bool smallFragment = false;
            for (const auto& atom : incident) {
                if (!visited.insert(atom.first).second) continue;
                ++fragments;
                std::vector<int> reached{atom.first};
                for (size_t position = 0; position < reached.size(); ++position)
                    for (int neighbor : neighbors[reached[position]])
                        if (visited.insert(neighbor).second) reached.push_back(neighbor);
                smallFragment |= reached.size() < static_cast<size_t>(options.minFragmentAtoms);
            }
            if (smallFragment || fragments > options.maxFragments
                || (options.connectivity == smsd::MCSConnectivity::Connected && fragments > 1)) continue;
            const Objective score{static_cast<int>(selected.size()), static_cast<int>(incident.size())};
            if (score < best) continue;
            if (best < score) { best = score; matches.clear(); }
            matches.emplace(std::move(incident), std::move(selected));
        }
    };
    enumerate(0);
    return {best, matches};
}

/**
 * Enumerate injective atom maps and every subset of their compatible bonds, and keep the best of
 * those that MCS::validateMapping() accepts: an oracle for stereochemistry, complete rings and
 * fragment limits, under which a mapping may leave compatible bonds out.
 */
static std::pair<Objective, std::set<MatchKey>> validatedOracle(
        const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& options) {
    Objective best{0, 0};
    std::set<MatchKey> matches;
    std::vector<int> mapping(query.n, -1);
    std::vector<bool> used(target.n, false);
    std::function<void(int)> enumerate = [&](int queryAtom) {
        if (queryAtom < query.n) {
            enumerate(queryAtom + 1);
            for (int targetAtom = 0; targetAtom < target.n; ++targetAtom) {
                if (used[targetAtom]
                    || !smsd::search_detail::atomCompatible(query, queryAtom, target, targetAtom, options.chemistry)) continue;
                mapping[queryAtom] = targetAtom;
                used[targetAtom] = true;
                enumerate(queryAtom + 1);
                used[targetAtom] = false;
                mapping[queryAtom] = -1;
            }
            return;
        }
        BondMapping compatible;
        for (int first = 0; first < query.n; ++first)
            for (int second : query.neighbors[first]) {
                if (first > second || mapping[first] < 0 || mapping[second] < 0
                    || !target.hasBond(mapping[first], mapping[second])) continue;
                if (smsd::search_detail::bondCompatible(query, first, second, target, mapping[first], mapping[second],
                                                        options.chemistry))
                    compatible.push_back({first, second, mapping[first], mapping[second]});
            }
        for (size_t subset = 1; subset < (size_t{1} << compatible.size()); ++subset) {
            std::vector<int> incident(query.n, -1);
            BondMapping selected;
            for (size_t edge = 0; edge < compatible.size(); ++edge) if (subset & (size_t{1} << edge)) {
                const auto& bond = compatible[edge];
                selected.push_back(bond);
                incident[bond[0]] = bond[2];
                incident[bond[1]] = bond[3];
            }
            const auto match = smsd::search_detail::makeMatch(incident, selected);
            const Objective score{static_cast<int>(match.bondMapping.size()), static_cast<int>(match.atomMapping.size())};
            if (score < best || !smsd::MCS::validateMapping(query, target, match, options).empty()) continue;
            if (best < score) { best = score; matches.clear(); }
            matches.emplace(match.atomMapping, match.bondMapping);
        }
    };
    enumerate(0);
    return {best, matches};
}

// Every best mapping under stereochemistry, complete rings and fragment limits, including those
// that leave a compatible bond out, as the exhaustive oracle finds them.
static void validatedSelections() {
    const std::pair<const char*, const char*> pairs[] = {
        {"C[C@H](N)O", "C[C@@H](N)O"}, {"N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"}, {"F/C=C/C", "F/C=C\\C"},
        {"F/C=C/C=C/C", "F/C=C\\C=C/C"}, {"C1CC1C", "CCCC"}, {"C1CCC1", "CC1CC1"},
        {"C[C@@H]1CC1N", "C[C@H]1CC1N"}, {"OC1CC1C", "OCC1CC1"}};
    for (const auto& pair : pairs) {
        const auto query = smsd::parseSMILES(pair.first), target = smsd::parseSMILES(pair.second);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (int flags = 1; flags < 4; ++flags)
                for (bool limits : {false, true}) {
                    smsd::MCSOptions options;
                    options.objective = smsd::MCSObjective::MCES;
                    options.connectivity = connectivity;
                    options.timeoutMs = 60000;
                    options.chemistry.useChirality = options.chemistry.useBondStereo = (flags & 1) != 0;
                    options.chemistry.completeRingsOnly = (flags & 2) != 0;
                    if (limits) {
                        options.minFragmentAtoms = 3;
                        options.maxFragments = 2;
                    }
                    const auto expected = validatedOracle(query, target, options);
                    const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 100000);
                    std::set<MatchKey> actual;
                    for (const auto& match : result.matches) actual.emplace(match.atomMapping, match.bondMapping);
                    require(result.optimal && result.complete && result.status == smsd::SearchStatus::Completed,
                            "Small validated searches must prove and enumerate their answer");
                    require(result.score == expected.first.first && actual == expected.second
                            && actual.size() == result.matches.size(),
                            "Selections under stereochemistry, rings or fragment limits differ from the oracle");
                    const auto single = smsd::MCS::find(query).withOptions(options).match(target);
                    require(single.optimal && single.score == expected.first.first
                            && (single.matches.empty() || expected.second.count({single.matches.front().atomMapping,
                                                                                  single.matches.front().bondMapping})),
                            "A single best selection must be one of the oracle's");
                }
    }
}

static void labelledSelectedEdges() {
    for (const auto& pair : {std::pair<const char*, const char*>{"C1CCC1", "CCCC"},
                            {"CC(C)C", "CCCC"}, {"NC(C)O", "CC(N)C"},
                            {"CCO.CN", "CCN.CO"}, {"N=CCN", "NCC=N"},
                            {"CCO", "NNN"}, {"CO.CN", "NCCN"}}) {
        const auto query = smsd::parseSMILES(pair.first), target = smsd::parseSMILES(pair.second);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (bool atomTypes : {false, true}) for (bool bondOrders : {false, true})
                for (bool excludeLast : {false, true}) {
                    smsd::MCSOptions options;
                    options.objective = smsd::MCSObjective::MCES;
                    options.connectivity = connectivity;
                    options.timeoutMs = 10000;
                    options.chemistry.matchAtomType = atomTypes;
                    options.chemistry.matchBondOrder = bondOrders ? smsd::ChemOptions::BondOrderMode::STRICT
                        : smsd::ChemOptions::BondOrderMode::ANY;
                    if (excludeLast) options.excludedTargetAtoms = {target.n - 1};
                    const auto expected = selectedEdgeOracle(query, target, options);
                    const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 10000);
                    require(result.optimal && result.complete && result.status == smsd::SearchStatus::Completed,
                            "Small selected-edge search must prove and enumerate its answer");
                    require(result.score == expected.first.first, "Selected-edge objective differs from the independent oracle");
                    std::set<MatchKey> actual;
                    for (const auto& match : result.matches) {
                        require(match.atomMapping.size() == static_cast<size_t>(expected.first.second),
                                "Incident-atom tie differs from the independent oracle");
                        require(smsd::MCS::validateMapping(query, target, match, options).empty(), "Invalid selected-edge match");
                        require(actual.emplace(match.atomMapping, match.bondMapping).second, "Repeated selected-edge embedding");
                    }
                    require(actual == expected.second, "Raw selected-edge embeddings differ from the independent oracle");
                }
    }
}

static void fixedMappingsAndFragments() {
    const auto molecule = smsd::parseSMILES("CCO");
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.timeoutMs = 10000;
    smsd::Match extendible;
    extendible.atomMapping = {{0, 0}, {1, 1}};
    extendible.bondMapping = {{0, 1, 0, 1}};
    require(!smsd::MCS::isMappingMaximal(molecule, molecule, extendible, options), "A fixed common edge has a compatible extension");
    smsd::Match blocked;
    blocked.atomMapping = {{0, 1}, {1, 0}};
    blocked.bondMapping = {{0, 1, 1, 0}};
    require(smsd::MCS::isMappingMaximal(molecule, molecule, blocked, options), "Maximality must preserve the fixed atom orientation");
    const auto larger = smsd::parseSMILES("CCCO"), smaller = smsd::parseSMILES("CCC");
    blocked.atomMapping = {{1, 0}, {2, 1}};
    blocked.bondMapping = {{1, 2, 0, 1}};
    require(smsd::MCS::isMappingMaximal(larger, smaller, blocked, options),
            "A fixed orientation must survive unequal query and target edge counts");

    const auto query = smsd::parseSMILES("CC.CC"), target = smsd::parseSMILES("CCCC");
    options.connectivity = smsd::MCSConnectivity::Fragments;
    for (int maxFragments : {1, 2}) for (int minAtoms : {2, 3}) {
        options.maxFragments = maxFragments;
        options.minFragmentAtoms = minAtoms;
        const auto expected = selectedEdgeOracle(query, target, options);
        const auto result = smsd::MCS::find(query).withOptions(options).matchAll(target, 1000);
        std::set<MatchKey> actual;
        for (const auto& match : result.matches) actual.emplace(match.atomMapping, match.bondMapping);
        require(result.optimal && result.complete && result.score == expected.first.first && actual == expected.second,
                "Fragment limits differ from the selected-edge oracle");
    }
}

static void stopsAndResultLimits() {
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.timeoutMs = 10000;
    const auto limited = smsd::MCS::find("C1CC1").withOptions(options).matchAll("CCC", 2);
    require(limited.optimal && !limited.complete && limited.status == smsd::SearchStatus::LimitReached
            && limited.matches.size() == 2 && limited.score == 2, "A tied-result limit must retain the optimum proof");
    options.cancellation = std::make_shared<smsd::CancellationToken>();
    options.cancellation->cancel();
    const auto cancelled = smsd::MCS::find("C1CCC1").withOptions(options).match("CC(C)C");
    require(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete,
            "A cancelled selected-edge search cannot claim proof");
    options.cancellation.reset();
    // The search needs more than one pair to prove this pair; its seed still gives a mapping.
    options.nodeLimit = 1;
    const auto stopped = smsd::MCS::find("CCN(CC)CCOC(=O)c1ccc(N)cc1").withOptions(options).match("CN1CCC(CC1)OC(=O)c1ccccc1");
    require(stopped.status == smsd::SearchStatus::LimitReached && !stopped.optimal && !stopped.complete
            && !stopped.matches.empty() && stopped.nodesVisited <= 1,
            "A stopped selected-edge search cannot claim proof");
    // A node limit reached while the whole query is embedded still leaves the seed.
    const auto embedding = smsd::MCS::find("c1ccc2ccccc2c1").withOptions(options).match("c1ccccc1CC");
    require(embedding.status == smsd::SearchStatus::LimitReached && !embedding.optimal && !embedding.matches.empty()
            && smsd::MCS::validateMapping(smsd::parseSMILES("c1ccc2ccccc2c1"), smsd::parseSMILES("c1ccccc1CC"),
                                          embedding.matches.front(), options).empty(),
            "A node limit reached while embedding must keep a valid seed");
}

static void symmetricMolecularPredicates() {
    for (const auto& pair : {std::pair<const char*, const char*>{"[13CH3]C(O)N", "C(O)N"},
                            {"C[NH2+]CO", "[NH2+]CO"}, {"OC[C@H](F)Cl", "C[C@@H](F)Cl"},
                            {"OCC/C=C/F", "C/C=C\\F"}, {"OC1CCC1", "C1CCC1"},
                            {"CC(=O)CN", "C=C(O)N"}}) {
        const auto query = smsd::parseSMILES(pair.first), target = smsd::parseSMILES(pair.second);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (bool completeRings : {false, true}) {
                smsd::MCSOptions options;
                options.objective = smsd::MCSObjective::MCES;
                options.connectivity = connectivity;
                options.timeoutMs = 10000;
                options.chemistry.matchFormalCharge = true;
                options.chemistry.matchIsotope = true;
                options.chemistry.useChirality = true;
                options.chemistry.useBondStereo = true;
                options.chemistry.tautomerAware = true;
                options.chemistry.completeRingsOnly = completeRings;
                const auto forward = smsd::MCS::find(query).withOptions(options).matchAll(target, 10000);
                const auto reverse = smsd::MCS::find(target).withOptions(options).matchAll(query, 10000);
                require(forward.optimal && forward.complete && reverse.optimal && reverse.complete,
                        "Small molecular selected-edge searches must prove both orientations");
                require(forward.score == reverse.score, "Molecular MCES objective must be symmetric");
                std::set<MatchKey> expected, actual;
                for (const auto& match : forward.matches) expected.emplace(match.atomMapping, match.bondMapping);
                for (const auto& match : reverse.matches) {
                    AtomMapping atoms;
                    BondMapping bonds;
                    for (const auto& atom : match.atomMapping) atoms.emplace(atom.second, atom.first);
                    for (const auto& bond : match.bondMapping) {
                        std::array<int, 4> inverted{bond[2], bond[3], bond[0], bond[1]};
                        if (inverted[0] > inverted[1]) {
                            std::swap(inverted[0], inverted[1]);
                            std::swap(inverted[2], inverted[3]);
                        }
                        bonds.push_back(inverted);
                    }
                    std::sort(bonds.begin(), bonds.end());
                    actual.emplace(std::move(atoms), std::move(bonds));
                }
                require(actual == expected, "Molecular MCES must retain every inverted selected-edge embedding");
                for (const auto& match : forward.matches)
                    require(smsd::MCS::validateMapping(query, target, match, options).empty(), "Invalid molecular selected-edge match");
            }
    }
}


static void connectedCageAndBicyclicMappings() {
    struct Pair { const char* query; const char* target; int atoms; int bonds; };
    for (const auto& pair : {Pair{"C12C3C4C1C5C2C3C45", "C1CCC2CCCCC2C1", 8, 8},
                             Pair{"C1C2CC12", "CCCCC", 4, 3}}) {
        const auto first = smsd::parseSMILES(pair.query), second = smsd::parseSMILES(pair.target);
        for (bool reverse : {false, true}) {
            const auto& query = reverse ? second : first;
            const auto& target = reverse ? first : second;
            for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
                smsd::MCSOptions options;
                options.objective = smsd::MCSObjective::MCES;
                options.connectivity = connectivity;
                options.maxFragments = 1;
                options.timeoutMs = 10000;
                const auto result = smsd::MCS::find(query).withOptions(options).match(target);
                require(result.optimal && result.complete && result.status == smsd::SearchStatus::Completed
                        && !result.matches.empty(),
                        "Connected cage and bicyclic searches must prove both input orientations");
                require(result.score == pair.bonds && result.matches.front().atomMapping.size() == static_cast<size_t>(pair.atoms)
                        && result.matches.front().fragments.size() == 1,
                        "Connected cage or bicyclic objective differs between input orientations");
                require(smsd::MCS::validateMapping(query, target, result.matches.front(), options).empty(),
                        "Invalid cage or bicyclic selected-edge mapping");
                if (query.n <= 5 && target.n <= 5) {
                    const auto expected = selectedEdgeOracle(query, target, options);
                    const auto all = smsd::MCS::find(query).withOptions(options).matchAll(target, 10000);
                    std::set<MatchKey> actual;
                    for (const auto& match : all.matches) actual.emplace(match.atomMapping, match.bondMapping);
                    require(all.optimal && all.complete && all.score == expected.first.first && actual == expected.second,
                            "Bicyclic raw selected-edge embeddings differ from the exhaustive oracle");
                }
            }
        }
    }
}

// With scouting from the first nodes, and refuted root pairs skipped, the common-bond search still
// proves the optimum it proves with its default limits.
static void scoutingKeepsTheOptimum() {
    const std::pair<const char*, const char*> pairs[] = {
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"}, {"CC(C)Cc1ccc(C(C)C(=O)O)cc1", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"},
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"}, {"c1ccc2ccccc2c1", "c1ccccc1CC"},
        {"C1CCC2CCCCC2C1", "C12C3C4C1C5C2C3C45"}, {"NCC(=O)O.CCO", "OCC(N)=O"},
        {"OC(=O)c1ccccc1O", "Nc1ccc(cc1)S(N)(=O)=O"}, {"CCN(CC)CCOC(=O)c1ccc(N)cc1", "CN1CCC(CC1)OC(=O)c1ccccc1"}};
    smsd::common_bond_detail::ScoutLimits limits;
    limits.after = 1;
    limits.growth = 2;
    limits.budget = 32;
    limits.root = 4;
    int64_t scoutings = 0, refuted = 0;
    for (const auto& pair : pairs) {
        const auto query = smsd::parseSMILES(pair.first), target = smsd::parseSMILES(pair.second);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            for (auto order : {smsd::ChemOptions::BondOrderMode::STRICT, smsd::ChemOptions::BondOrderMode::ANY}) {
                smsd::MCSOptions options;
                options.objective = smsd::MCSObjective::MCES;
                options.connectivity = connectivity;
                options.chemistry.matchBondOrder = order;
                options.timeoutMs = 60000;
                // An unreachable node limit leaves the search and its objective unchanged.
                options.nodeLimit = UINT64_MAX;
                const auto general = smsd::MCS::find(query).withOptions(options).match(target);
                const auto& chemistry = options.chemistry;
                const auto outcome = smsd::common_bond_detail::searchCommonBonds(query, target,
                    [&](int queryAtom, int targetAtom) {
                        return smsd::search_detail::atomCompatible(query, queryAtom, target, targetAtom, chemistry);
                    },
                    [&](int queryFirst, int querySecond, int targetFirst, int targetSecond) {
                        return smsd::search_detail::bondCompatible(query, queryFirst, querySecond, target, targetFirst, targetSecond, chemistry);
                    },
                    [](const smsd::MolGraph& graph, int first, int second) {
                        return smsd::search_detail::bondCompatibilityKey(graph, first, second);
                    },
                    {}, connectivity == smsd::MCSConnectivity::Connected, [] { return false; }, limits);
                const auto atoms = std::count_if(outcome.mapping.begin(), outcome.mapping.end(), [](int image) { return image >= 0; });
                require(general.optimal && !general.matches.empty() && outcome.supported && outcome.finished,
                        "Both searches must prove their optimum");
                require(static_cast<double>(outcome.bonds.size()) == general.score
                        && static_cast<size_t>(atoms) == general.matches.front().atomMapping.size(),
                        "Scouting from the first nodes changed the proved optimum");
                scoutings += outcome.scoutings;
                refuted += outcome.refutedRoots;
            }
    }
    require(scoutings > 0 && refuted > 0, "Lowered scout limits must scout and skip refuted root pairs");
}

// The greedy seed grows from the query, and from the other molecule it may hold far more common
// bonds. Here the seed grown from the smaller molecule has 23 common bonds and the one grown from
// the larger has the optimum, 26. Once its levels of common bonds run out of budget the fragment
// search starts from the better seed, and proves the optimum within a node limit that the first
// seed alone needs over two million nodes to reach, whichever molecule is given as the query.
static void seedFromTheOtherMoleculeProvesTheChainPair() {
    const auto query = smsd::parseSMILES("CCCCCCCCCCCCCCSCC(COP(O)(OP(O)(OCC1OC(n2ccc(N)nc2=O)C(O)C1O)=O)=O)OC(=O)CCCCCCCCCCC");
    const auto target = smsd::parseSMILES("CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1");
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.connectivity = smsd::MCSConnectivity::Fragments;
    options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::ANY;
    options.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::FLEXIBLE;
    options.timeoutMs = 60000;
    options.nodeLimit = 100000;
    for (int turn = 0; turn < 2; ++turn) {
        const auto result = turn == 0 ? smsd::MCS::find(query).withOptions(options).match(target)
                                      : smsd::MCS::find(target).withOptions(options).match(query);
        require(result.optimal && !result.matches.empty(), "The chain pair must be proved within the node limit");
        require(result.matches.front().bondMapping.size() == 26 && result.matches.front().atomMapping.size() == 30,
                "The chain pair has 26 common bonds on 30 atoms");
    }
}

/**
 * The kinds of the atoms read whether two kinds match from the pair rows, so they allocate no
 * table of every pair of kinds. With as many kinds as atoms, as atoms with distinct labels have,
 * such a table would take a byte for each pair, here 144 MB, eight times the rows themselves.
 */
static void kindsAddNoPairTable() {
    namespace detail = smsd::common_bond_detail;
    const int atoms = 12000, words = detail::words(atoms);
    std::vector<uint64_t> rows(size_t(atoms) * size_t(words), 0);
    auto image = [&](int atom) { return int((int64_t(atom) * 7919) % atoms); };
    for (int atom = 0; atom < atoms; ++atom) detail::setBit(rows, words, atom, image(atom));
    detail::Clock clock([] { return false; });
    largestAllocation = 0;
    trackAllocations = true;
    const detail::Kinds kinds(rows, atoms, atoms, clock);
    trackAllocations = false;
    require(kinds.queryKinds == atoms && kinds.targetKinds == atoms && kinds.classes == atoms, "one kind and class for each atom");
    require(largestAllocation <= rows.size() * sizeof(uint64_t), "the kinds allocated a table of every pair of kinds");
    for (int atom = 0; atom < atoms; atom += 37) {
        require(kinds.matches(kinds.queryKind[atom], kinds.targetKind[image(atom)]), "a pair of the rows must match");
        require(!kinds.matches(kinds.queryKind[atom], kinds.targetKind[image(atom + 1)]), "a pair not in the rows must not match");
    }
}

int main() {
    try {
        kindsAddNoPairTable();
        labelledSelectedEdges();
        validatedSelections();
        fixedMappingsAndFragments();
        stopsAndResultLimits();
        symmetricMolecularPredicates();
        connectedCageAndBicyclicMappings();
        scoutingKeepsTheOptimum();
        seedFromTheOtherMoleculeProvesTheChainPair();
        std::cout << "PASS labelled selected-edge oracles, fixed mappings, fragment limits and stops\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}

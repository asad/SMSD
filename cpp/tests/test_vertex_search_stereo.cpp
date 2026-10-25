/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The stereo test of each new pair in the general MCIS search.
 */
#include "test_support.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <array>
#include <functional>
#include <set>
#include <string>
#include <vector>

using smsd_test::require;

namespace {

/**
 * Centres with an implicit hydrogen, four neighbours or a lone pair, partly configured
 * molecules, ring centres, and double bonds in chains, with three substituents and in rings.
 * The same molecules are in the Java tests.
 */
const std::vector<std::string> molecules{
    "C[C@H](O)CC", "C[C@@H](O)CC", "C[C@H](O)C(C)O", "C[C@@H](O)[C@H](C)O", "C[C@H](O)[C@H](C)O",
    "C[C@](F)(Cl)CO", "C[S@](=O)CC", "C[C@H]1CC[C@@H]1O", "C[C@H]1CC[C@H]1O", "C/C=C/C",
    "C/C=C\\C", "C/C(Cl)=C/C", "CC=CC", "C/C=C/[C@H](C)O", "C/C=C\\[C@@H](C)O", "O/C=C/C=C\\C",
    "C/C=C1/CC[C@H]1C"};

/** The atom mapping and the selected bonds of a match, written as the Java tests write them. */
std::string key(const smsd::Match& match) {
    std::string text = "{";
    for (const auto& atom : match.atomMapping)
        text += (text.size() > 1 ? ", " : "") + std::to_string(atom.first) + "=" + std::to_string(atom.second);
    std::vector<std::string> bonds;
    for (const auto& bond : match.bondMapping)
        bonds.push_back(std::to_string(std::min(bond[0], bond[1])) + "-" + std::to_string(std::max(bond[0], bond[1])));
    std::sort(bonds.begin(), bonds.end());
    text += "} [";
    for (size_t index = 0; index < bonds.size(); ++index) text += (index ? ", " : "") + bonds[index];
    return text + "]";
}

/**
 * Every maximum induced common subgraph that MCS::validateMapping accepts. The atom maps are
 * enumerated with only the induced-topology and bond tests, which every part of a valid MCIS
 * passes, and the stereo test is left to the validation.
 */
std::set<std::string> oracle(const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& options) {
    const auto& chemistry = options.chemistry;
    std::vector<std::vector<std::vector<int>>> bySize(static_cast<size_t>(query.n) + 1);
    std::vector<int> mapping(query.n, -1);
    std::vector<bool> used(target.n, false);
    std::function<void(int, int)> enumerate = [&](int atom, int size) {
        if (atom == query.n) {
            bySize[size].push_back(mapping);
            return;
        }
        enumerate(atom + 1, size);
        for (int image = 0; image < target.n; ++image) {
            if (used[image] || !smsd::search_detail::atomCompatible(query, atom, target, image, chemistry)) continue;
            bool induced = true;
            for (int earlier = 0; earlier < atom && induced; ++earlier) {
                if (mapping[earlier] < 0) continue;
                const bool bonded = query.hasBond(atom, earlier);
                induced = bonded == target.hasBond(image, mapping[earlier])
                    && (!bonded || smsd::search_detail::bondCompatible(query, atom, earlier, target, image,
                                                                       mapping[earlier], chemistry));
            }
            if (!induced) continue;
            mapping[atom] = image;
            used[image] = true;
            enumerate(atom + 1, size + 1);
            used[image] = false;
            mapping[atom] = -1;
        }
    };
    smsd::search_detail::prepare(&query, target, chemistry);
    enumerate(0, 0);
    for (int size = query.n; size > 0; --size) {
        std::set<std::string> keys;
        int bonds = -1;
        for (const auto& candidate : bySize[size]) {
            std::vector<std::array<int, 4>> selected;
            for (int atom = 0; atom < query.n; ++atom)
                if (candidate[atom] >= 0)
                    for (int neighbor : query.neighbors[atom])
                        if (neighbor < atom && candidate[neighbor] >= 0)
                            selected.push_back({neighbor, atom, candidate[neighbor], candidate[atom]});
            if (static_cast<int>(selected.size()) < bonds) continue;
            const auto match = smsd::search_detail::makeMatch(candidate, selected);
            if (!smsd::MCS::validateMapping(query, target, match, options).empty()) continue;
            if (static_cast<int>(selected.size()) > bonds) {
                keys.clear();
                bonds = static_cast<int>(selected.size());
            }
            keys.insert(key(match));
        }
        if (!keys.empty()) return keys;
    }
    return {};
}

smsd::MCSOptions options(bool chirality, bool bondStereo, smsd::MCSConnectivity connectivity) {
    smsd::MCSOptions result;
    result.chemistry.useChirality = chirality;
    result.chemistry.useBondStereo = bondStereo;
    result.connectivity = connectivity;
    result.timeoutMs = 60000;
    return result;
}

/** matchAll() must give every mapping of the enumeration, and match() one of them. */
std::set<std::string> requireEnumerated(const smsd::MolGraph& query, const smsd::MolGraph& target, const smsd::MCSOptions& opts,
                       const std::string& label) {
    const auto expected = oracle(query, target, opts);
    const auto prepared = smsd::MCS::find(query).withOptions(opts);
    const auto all = prepared.matchAll(target, 100000);
    require(all.status == smsd::SearchStatus::Completed && all.optimal && all.complete, label + ": matchAll stopped");
    std::set<std::string> actual;
    for (const auto& match : all.matches) require(actual.insert(key(match)).second, label + ": repeated match");
    require(actual == expected, label + ": matchAll differs from the enumeration");
    const auto one = prepared.match(target);
    require(one.optimal && one.matches.empty() == expected.empty(), label + ": match() differs");
    if (!one.matches.empty())
        require(expected.count(key(one.matches.front())) == 1, label + ": match() is not a maximum");
    return actual;
}

const char* name(smsd::MCSConnectivity connectivity) {
    return connectivity == smsd::MCSConnectivity::Connected ? " connected" : " fragments";
}

void testExhaustiveStereoOracle() {
    const std::array<std::array<bool, 2>, 3> stereo{{{true, false}, {false, true}, {true, true}}};
    size_t searches = 0;
    for (const auto& querySmiles : molecules)
        for (const auto& targetSmiles : molecules) {
            const auto query = smsd::parseSMILES(querySmiles), target = smsd::parseSMILES(targetSmiles);
            for (const auto& flags : stereo)
                for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
                    requireEnumerated(query, target, options(flags[0], flags[1], connectivity),
                                      querySmiles + " " + targetSmiles + " " + std::to_string(flags[0])
                                          + std::to_string(flags[1]) + name(connectivity));
                    ++searches;
                }
        }
    require(searches == molecules.size() * molecules.size() * 6, "every search must run");
}

/**
 * An alkene from the Builder: atom 0 doubly bonded to atom 1 in the given configuration, single
 * bonds elsewhere and no hydrogens, with the neighbours in the order given.
 */
smsd::MolGraph alkene(const std::vector<int>& elements, const std::vector<std::vector<int>>& neighbors, int configuration) {
    const int n = static_cast<int>(elements.size());
    std::vector<std::vector<int>> orders(n);
    for (int atom = 0; atom < n; ++atom)
        for (int neighbor : neighbors[atom]) orders[atom].push_back(atom <= 1 && neighbor <= 1 ? 2 : 1);
    std::vector<std::vector<int>> stereo(n, std::vector<int>(n, 0));
    stereo[0][1] = stereo[1][0] = configuration;
    return smsd::MolGraph::Builder().atomCount(n).atomicNumbers(elements).hydrogenCounts(std::vector<int>(n, 0))
        .setNeighbors(neighbors).setBondOrders(std::move(orders)).doubleBondStereo(std::move(stereo)).build(false);
}

/**
 * Double bonds with an end of two substituents in one graph and three in the other. The test
 * reads the first mapped substituent at each end, and here a later one can change its answer, so
 * these bonds are left to the test of complete mappings; testing them as each pair is made loses
 * a maximum in the first three pairs and a tied mapping in the others. The Java tests use the
 * same graphs.
 */
void testCrowdedDoubleBondEnds() {
    struct Graph {
        std::vector<int> elements;
        std::vector<std::vector<int>> neighbors;
        int configuration;
    };
    const std::vector<std::pair<Graph, Graph>> pairs{
        {{{6, 6, 6, 7, 6, 6}, {{1, 2, 3}, {0, 5}, {0}, {0, 4}, {3}, {1}}, 1},
         {{6, 6, 6, 6, 6, 6, 7, 6}, {{4, 6, 2, 1}, {0, 7}, {3, 0}, {2}, {5, 0}, {4}, {0}, {1}}, 2}},
        {{{6, 6, 7, 6, 6, 6, 6}, {{1, 2, 4, 3}, {6, 0}, {0}, {0}, {5, 0}, {4}, {1}}, 1},
         {{6, 6, 6, 6, 6, 6}, {{1, 4, 2}, {5, 0}, {0, 3}, {2}, {0}, {1}}, 1}},
        {{{6, 6, 7, 6, 6, 7, 6}, {{4, 1, 2}, {6, 0, 5}, {0, 3}, {2}, {0}, {1}, {1}}, 1},
         {{6, 6, 6, 7, 6, 7, 6, 6}, {{1, 3, 2}, {6, 0, 4, 5}, {0}, {0}, {1}, {1}, {7, 1}, {6}}, 1}},
        {{{6, 6, 6, 6, 6, 6}, {{1, 2}, {0, 4, 3}, {0}, {1}, {5, 1}, {4}}, 2},
         {{6, 6, 6, 6, 6, 6, 6}, {{1, 2}, {0, 3, 4, 5}, {0}, {1}, {1}, {6, 1}, {5}}, 1}},
        {{{6, 6, 6, 6, 6, 6, 6, 6}, {{2, 1}, {7, 0, 5, 4}, {3, 0}, {2}, {1}, {1, 6}, {5}, {1}}, 2},
         {{6, 6, 6, 6, 6, 6}, {{2, 3, 1}, {5, 0}, {0}, {4, 0}, {3}, {1}}, 2}},
        {{{6, 6, 7, 7, 6, 7, 6}, {{2, 3, 1}, {6, 5, 0, 4}, {0}, {0}, {1}, {1}, {1}}, 1},
         {{6, 6, 6, 7, 7, 6, 7}, {{1, 4, 3, 2}, {6, 5, 0}, {0}, {0}, {0}, {1}, {1}}, 1}}};
    for (size_t index = 0; index < pairs.size(); ++index) {
        const auto& [first, second] = pairs[index];
        const auto query = alkene(first.elements, first.neighbors, first.configuration);
        const auto target = alkene(second.elements, second.neighbors, second.configuration);
        for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments})
            requireEnumerated(query, target, options(false, true, connectivity), "pair " + std::to_string(index) + name(connectivity));
    }
}

/**
 * Epimers from the benchmark pairs: doxorubicin and epirubicin, two fluvastatin diastereomers
 * and two cefotaxime epimers, with chirality and double-bond stereo. A pair that breaks a centre
 * is passed over before its node, which keeps these searches small. The Java tests pin their own
 * counts of the same searches, whose results are the same.
 */
void testPairsThatBreakACentre() {
    const std::string doxorubicin = "COc1cccc2c1C(=O)c1c([O-])c3c(c([O-])c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3"
                                    "O[C@H]1C[C@H]([NH3+])[C@@H](O)[C@H](C)O1";
    const std::string epirubicin = "COc1cccc2c1C(=O)c1c([O-])c3c(c([O-])c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3"
                                   "O[C@H]1C[C@H]([NH3+])[C@H](O)[C@H](C)O1";
    const std::string fluvastatin = "CC(C)n1c(/C=C/[C@H](O)C[C@H](O)CC(=O)[O-])c(-c2ccc(F)cc2)c2ccccc21";
    const std::string diastereomer = "CC(C)n1c(/C=C/[C@@H](O)C[C@H](O)CC(=O)O)c(-c2ccc(F)cc2)c2ccccc21";
    const std::string cefotaxime = "CO/N=C(\\C(=O)N[C@H]1C(=O)N2C(C(=O)[O-])=C(COC(C)=O)CS[C@@H]12)c1csc(N)n1.[Na+]";
    const std::string epimer = "CO/N=C(\\C(=O)N[C@@H]1C(=O)N2C(C(=O)[O-])=C(COC(C)=O)CS[C@H]12)c1csc(N)n1";
    struct Case {
        std::string query, target;
        smsd::MCSConnectivity connectivity;
        size_t atoms, bonds;
        uint64_t matchNodes;
        size_t mappings;
        uint64_t allNodes;
    };
    const auto connected = smsd::MCSConnectivity::Connected, fragments = smsd::MCSConnectivity::Fragments;
    const std::vector<Case> cases{
        {doxorubicin, epirubicin, connected, 38, 42, 191, 1, 376},
        {doxorubicin, epirubicin, fragments, 38, 42, 193, 1, 378},
        {fluvastatin, diastereomer, fragments, 29, 31, 169, 4, 300},
        {cefotaxime, epimer, connected, 29, 29, 206, 1, 350}};
    for (const auto& test : cases) {
        const auto query = smsd::parseSMILES(test.query), target = smsd::parseSMILES(test.target);
        const auto prepared = smsd::MCS::find(query).withOptions(options(true, true, test.connectivity));
        const std::string label = test.query + " " + test.target;
        const auto one = prepared.match(target), all = prepared.matchAll(target, 1000);
        require(one.optimal && all.optimal && all.complete, label + ": search stopped");
        require(one.matches.front().atomMapping.size() == test.atoms, label + ": atoms");
        require(one.matches.front().bondMapping.size() == test.bonds, label + ": bonds");
        require(one.nodesVisited == test.matchNodes, label + ": match() nodes");
        require(all.matches.size() == test.mappings, label + ": mappings");
        require(all.nodesVisited == test.allNodes, label + ": matchAll() nodes");
        bool found = false;
        for (const auto& match : all.matches) found = found || key(match) == key(one.matches.front());
        require(found, label + ": match() is not among matchAll()");
    }
}

/** Small tetrahedral stars with supplied hydrogen counts and reordered ligand rows. */
smsd::MolGraph tetrahedralStar(int ligands, bool explicitHydrogen, bool permuted) {
    const int atoms = ligands + 1;
    const std::vector<int> elements{6, 9, 17, 35, explicitHydrogen ? 1 : 53};
    std::vector<int> labels(atoms), hydrogens(atoms, 0), tags(atoms, 0), ligandOrder;
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    auto image = [&](int atom) { return permuted ? (atom == 0 ? atoms - 1 : atom - 1) : atom; };
    for (int atom = 0; atom < atoms; ++atom) labels[image(atom)] = elements[atom];
    for (int ligand = 1; ligand <= ligands; ++ligand) ligandOrder.push_back(ligand);
    if (permuted) std::swap(ligandOrder[0], ligandOrder[1]);
    const int centre = image(0), configuration = explicitHydrogen ? 2 : 1;
    hydrogens[centre] = ligands == 3 ? 1 : 0;
    tags[centre] = permuted && ligands == 4 ? 3 - configuration : configuration;
    for (int ligand : ligandOrder) {
        neighbors[centre].push_back(image(ligand));
        orders[centre].push_back(1);
        neighbors[image(ligand)] = {centre};
        orders[image(ligand)] = {1};
    }
    return smsd::MolGraph::Builder().atomCount(atoms).atomicNumbers(labels).hydrogenCounts(hydrogens)
        .tetrahedralChirality(tags).setNeighbors(neighbors).setBondOrders(orders).build(false);
}

void testRelaxedLabelsAndHydrogenPermutations() {
    for (auto connectivity : {smsd::MCSConnectivity::Connected, smsd::MCSConnectivity::Fragments}) {
        for (bool stereo : {false, true}) {
            auto opts = options(stereo, stereo, connectivity);
            opts.chemistry.matchAtomType = false;
            for (int kind = 0; kind < 3; ++kind) {
                const int ligands = kind == 0 ? 3 : 4;
                const bool explicitHydrogen = kind == 2;
                const auto query = tetrahedralStar(ligands, explicitHydrogen, false);
                const auto target = tetrahedralStar(ligands, explicitHydrogen, true);
                const size_t permutations = ligands == 3 ? 6 : 24;
                const auto actual = requireEnumerated(query, target, opts, "relaxed tetrahedral permutations");
                require(actual.size() == (stereo ? permutations / 2 : permutations),
                        "Tetrahedral winding must select exactly half of the ligand permutations");
            }
            const auto query = alkene({6, 6, 9, 17, 35, 53}, {{1, 2, 3}, {0, 4, 5}, {0}, {0}, {1}, {1}}, 1);
            const auto target = alkene({6, 6, 9, 17, 35, 53}, {{3, 2, 1}, {5, 4, 0}, {0}, {0}, {1}, {1}}, 1);
            const auto actual = requireEnumerated(query, target, opts, "relaxed alkene permutations");
            require(actual.size() == (stereo ? 4U : 8U),
                    "Double-bond winding must select exactly half of the alkene permutations");
        }
        const auto implicit = tetrahedralStar(3, false, false), explicitH = tetrahedralStar(4, true, true);
        const auto opts = options(true, true, connectivity);
        requireEnumerated(implicit, explicitH, opts, "implicit to explicit hydrogen");
        requireEnumerated(explicitH, implicit, opts, "explicit to implicit hydrogen");
    }
}

void testStereoStopsRetainOnlyFeasibleWitnesses() {
    const auto query = tetrahedralStar(4, false, false), target = tetrahedralStar(4, false, true);
    for (bool all : {false, true}) {
        for (uint64_t cap : {uint64_t(1), uint64_t(2)}) {
            auto opts = options(true, true, smsd::MCSConnectivity::Connected);
            opts.chemistry.matchAtomType = false;
            opts.nodeLimit = cap;
            const auto prepared = smsd::MCS::find(query).withOptions(opts);
            const auto stopped = all ? prepared.matchAll(target, 100000) : prepared.match(target);
            require(stopped.status == smsd::SearchStatus::LimitReached && !stopped.optimal && !stopped.complete,
                    "A tiny cap must retain the stopped proof state");
            for (const auto& match : stopped.matches)
                require(smsd::MCS::validateMapping(query, target, match, opts).empty(), "Invalid capped stereo witness");
        }
        auto opts = options(true, true, smsd::MCSConnectivity::Connected);
        opts.cancellation = std::make_shared<smsd::CancellationToken>();
        opts.cancellation->cancel();
        const auto prepared = smsd::MCS::find(query).withOptions(opts);
        const auto cancelled = all ? prepared.matchAll(target, 100000) : prepared.match(target);
        require(cancelled.status == smsd::SearchStatus::Cancelled && !cancelled.optimal && !cancelled.complete
            && cancelled.matches.empty(), "Pre-cancelled stereo search must retain cancellation");
    }
}

} // namespace

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"exhaustive stereo oracle", testExhaustiveStereoOracle},
        {"crowded double-bond ends", testCrowdedDoubleBondEnds},
        {"relaxed labels and hydrogen permutations", testRelaxedLabelsAndHydrogenPermutations},
        {"stereo stops", testStereoStopsRetainOnlyFeasibleWitnesses},
        {"pairs that break a centre", testPairsThatBreakACentre},
    });
}

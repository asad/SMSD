/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/mcs.hpp"

#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

struct Graph {
    int atoms;
    std::vector<std::vector<int>> edges;
    smsd::MolGraph molecule;
};

Graph makeGraph(int atoms, const std::vector<std::pair<int, int>>& bonds,
                std::vector<int> atomicNumbers = {}) {
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    std::vector<std::vector<int>> edges(atoms, std::vector<int>(atoms));
    for (const auto& [a, b] : bonds) {
        edges[a][b] = edges[b][a] = 1;
        neighbors[a].push_back(b);
        neighbors[b].push_back(a);
        orders[a].push_back(1);
        orders[b].push_back(1);
    }
    if (atomicNumbers.empty()) atomicNumbers.assign(atoms, 6);
    auto molecule = smsd::MolGraph::Builder()
        .atomCount(atoms).atomicNumbers(atomicNumbers)
        .setNeighbors(neighbors).setBondOrders(orders).build();
    return {atoms, std::move(edges), std::move(molecule)};
}

Graph graphFromMask(int atoms, int mask) {
    std::vector<std::pair<int, int>> bonds;
    int bit = 0;
    for (int a = 0; a < atoms; ++a)
        for (int b = a + 1; b < atoms; ++b, ++bit)
            if (mask & (1 << bit)) bonds.emplace_back(a, b);
    return makeGraph(atoms, bonds);
}

bool connectedMapping(const Graph& query, const std::vector<int>& mapping) {
    int mappedCount = 0, root = -1;
    for (int atom = 0; atom < query.atoms; ++atom) {
        if (mapping[atom] < 0) continue;
        ++mappedCount;
        root = atom;
    }
    if (root < 0) return true;
    std::vector<int> seen(query.atoms), frontier{root};
    seen[root] = 1;
    for (size_t head = 0; head < frontier.size(); ++head) {
        for (int atom = 0; atom < query.atoms; ++atom) {
            if (mapping[atom] < 0 || seen[atom]
                || !query.edges[frontier[head]][atom]) continue;
            seen[atom] = 1;
            frontier.push_back(atom);
        }
    }
    return static_cast<int>(frontier.size()) == mappedCount;
}

// Independent oracle: enumerate every partial injective vertex mapping using
// only these adjacency matrices. No SMSD compatibility, bound, or search helper
// participates in the expected result or in returned-mapping validation.
int maximumAtoms(const Graph& query, const Graph& target,
                 bool induced, bool connected) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    int best = 0;
    std::function<void(int, int)> visit = [&](int atom, int count) {
        if (atom == query.atoms) {
            if (!connected || connectedMapping(query, mapping))
                best = std::max(best, count);
            return;
        }
        for (int other = 0; other < target.atoms; ++other) {
            if (used[other]) continue;
            bool valid = true;
            for (int previous = 0; previous < atom; ++previous) {
                if (mapping[previous] < 0) continue;
                bool qEdge = query.edges[atom][previous] != 0;
                bool tEdge = target.edges[other][mapping[previous]] != 0;
                if ((qEdge && !tEdge) || (induced && qEdge != tEdge)) valid = false;
            }
            if (!valid) continue;
            mapping[atom] = other;
            used[other] = 1;
            visit(atom + 1, count + 1);
            used[other] = 0;
            mapping[atom] = -1;
        }
        visit(atom + 1, count);
    };
    visit(0, 0);
    return best;
}

void validate(const Graph& query, const Graph& target,
              const std::map<int, int>& result, bool induced, bool connected) {
    std::vector<int> mapping(query.atoms, -1), used(target.atoms);
    for (const auto& [atom, other] : result) {
        assert(atom >= 0 && atom < query.atoms);
        assert(other >= 0 && other < target.atoms);
        assert(!used[other]);
        used[other] = 1;
        mapping[atom] = other;
    }
    for (const auto& [atom, other] : result) {
        for (const auto& [previous, mappedPrevious] : result) {
            bool qEdge = query.edges[atom][previous] != 0;
            bool tEdge = target.edges[other][mappedPrevious] != 0;
            assert(!qEdge || tEdge);
            assert(!induced || qEdge == tEdge);
        }
    }
    assert(!connected || connectedMapping(query, mapping));
}

smsd::MCSOptions options(bool induced, bool connected) {
    smsd::MCSOptions result;
    result.induced = induced;
    result.connectedOnly = connected;
    result.disconnectedMCS = !connected;
    result.timeoutMs = 1000;
    return result;
}

void testOptionRegressions() {
    auto fragments = makeGraph(5, {{0, 1}, {2, 3}, {3, 4}});
    auto connected = smsd::findMCS(
        fragments.molecule, fragments.molecule, {}, options(false, true));
    assert(connected.size() == 3); // Identity still obeys connectedOnly.
    validate(fragments, fragments, connected, false, true);

    auto dOpts = options(false, false);
    dOpts.connectedOnly = true; // Explicit dMCS takes precedence.
    assert(smsd::findMCS(fragments.molecule, fragments.molecule, {}, dOpts).size() == 5);
    dOpts.minFragmentSize = 3;
    dOpts.maxFragments = 1;
    auto fragment = smsd::findMCS(fragments.molecule, fragments.molecule, {}, dOpts);
    assert(fragment.size() == 3);
    validate(fragments, fragments, fragment, false, true);

    auto path = graphFromMask(3, 5), triangle = graphFromMask(3, 7);
    auto induced = smsd::findMCS(path.molecule, triangle.molecule, {}, options(true, true));
    assert(induced.size() == 2); // A cycle is not a linear-path shortcut.
    validate(path, triangle, induced, true, true);

    auto isolates = graphFromMask(3, 0), edge = graphFromMask(3, 1);
    auto forward = smsd::findMCS(isolates.molecule, edge.molecule, {}, options(false, false));
    auto reverse = smsd::findMCS(edge.molecule, isolates.molecule, {}, options(false, false));
    assert(forward.size() == 3 && reverse.size() == 2);
    validate(isolates, edge, forward, false, false);
    validate(edge, isolates, reverse, false, false);

    auto partialInduced = smsd::findMCS(
        edge.molecule, isolates.molecule, {}, options(true, false));
    assert(partialInduced.size() == 2); // Omitted neighbors do not impose degrees.
    validate(edge, isolates, partialInduced, true, false);

    auto star = graphFromMask(4, 7), edgeAndIsolates = graphFromMask(4, 1);
    auto commonEdge = smsd::findMCS(
        edgeAndIsolates.molecule, star.molecule, {}, options(false, true));
    assert(commonEdge.size() == 2); // The largest raw seed can have a worse component.
    validate(edgeAndIsolates, star, commonEdge, false, true);

    auto path4 = graphFromMask(4, 13);
    auto inducedPartial = smsd::findMCS(
        edgeAndIsolates.molecule, path4.molecule, {}, options(true, false));
    assert(inducedPartial.size() == 3); // Invalid augmented seeds must not stop exact search.
    validate(edgeAndIsolates, path4, inducedPartial, true, false);

    auto twoIsolates = graphFromMask(2, 0);
    auto unequal = smsd::findMCS(
        edge.molecule, twoIsolates.molecule, {}, options(false, false));
    assert(unequal.size() == 2); // Reverse containment is not a caller-direction proof.
    validate(edge, twoIsolates, unequal, false, false);

    // Padding bypasses the tiny exact solver and checks heuristic promotions.
    auto paddedEdge = makeGraph(9, {{0, 1}});
    auto paddedStar = makeGraph(9, {{0, 1}, {0, 2}, {0, 3}});
    auto padded = smsd::findMCS(
        paddedEdge.molecule, paddedStar.molecule, {}, options(false, true));
    assert(padded.size() == 2);
    validate(paddedEdge, paddedStar, padded, false, true);
}

void testWeightRegressions() {
    auto clique = graphFromMask(4, 63);
    auto path = makeGraph(5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}});
    auto weighted = options(false, true);
    weighted.atomWeights = {1, 1, 1, 10};
    auto result = smsd::findMCS(clique.molecule, path.molecule, {}, weighted);
    assert(result.size() == 2 && result.count(3) == 1);
    validate(clique, path, result, false, true);

    auto fragments = makeGraph(6, {{0, 1}, {1, 2}, {2, 3}, {4, 5}});
    weighted.atomWeights = {1, 1, 1, 1, 5, 5};
    auto heavy = smsd::findMCS(fragments.molecule, fragments.molecule, {}, weighted);
    assert(heavy.size() == 2 && heavy.count(4) && heavy.count(5));
    validate(fragments, fragments, heavy, false, true);

    auto bonded = makeGraph(9, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4},
        {5, 6}, {5, 7}, {5, 8}, {6, 7}, {6, 8}, {7, 8}});
    auto bondOptions = options(false, true);
    bondOptions.maximizeBonds = true;
    auto dense = smsd::findMCS(bonded.molecule, bonded.molecule, {}, bondOptions);
    assert(dense.size() == 4 && dense.count(5) && dense.count(8));
    validate(bonded, bonded, dense, false, true);

    weighted.atomWeights = {1};
    bool threw = false;
    try {
        (void)smsd::findMCS(clique.molecule, clique.molecule, {}, weighted);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw); // Identity must validate per-query options too.

    auto reject = [&](std::vector<double> weights) {
        auto invalid = options(false, true);
        invalid.atomWeights = std::move(weights);
        bool rejected = false;
        try {
            (void)smsd::findMCS(clique.molecule, clique.molecule, {}, invalid);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    };
    reject({std::numeric_limits<double>::quiet_NaN(), 0, 0, 0});
    reject({std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({-std::numeric_limits<double>::infinity(), 0, 0, 0});
    reject({std::numeric_limits<double>::max(), 0, 0, 0});
    reject({1.2e6, 1.2e6, -1.2e6, -1.2e6}); // Full score cancels, partial score overflows.
    reject({0, 0, -1.2e6, -1.2e6});

    auto atom = graphFromMask(1, 0);
    auto boundary = options(false, true);
    boundary.atomWeights = {(static_cast<double>(INT_MAX) + 0.5) / 1000.0};
    auto full = smsd::findMCS(atom.molecule, atom.molecule, {}, boundary);
    assert(smsd::detail::mcsScore(atom.molecule, full, boundary) == INT_MAX);
    boundary.atomWeights = {(static_cast<double>(INT_MIN) - 0.5) / 1000.0};
    assert(smsd::detail::mcsScore(atom.molecule, full, boundary) == INT_MIN);
    int q2t[] = {0};
    assert(smsd::detail::mcsScoreFlat(atom.molecule, q2t, 1, 1, boundary) == INT_MIN);

    // Score helpers also protect direct/internal callers that bypass API validation.
    boundary.atomWeights = {std::numeric_limits<double>::max()};
    for (bool flat : {false, true}) {
        bool rejected = false;
        try {
            if (flat) (void)smsd::detail::mcsScoreFlat(atom.molecule, q2t, 1, 1, boundary);
            else (void)smsd::detail::mcsScore(atom.molecule, full, boundary);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        assert(rejected);
    }
}

void testTautomerBound() {
    auto query = graphFromMask(4, 7), target = graphFromMask(4, 7);
    query.molecule.atomicNum.assign(4, 7);
    target.molecule.atomicNum.assign(4, 8);
    query.molecule.tautomerClass.assign(4, 0);
    target.molecule.tautomerClass.assign(4, 0);
    query.molecule.refreshAtomLabels();
    target.molecule.refreshAtomLabels();
    smsd::ChemOptions chemistry;
    chemistry.tautomerAware = true;
    // The identity witness has four compatible atoms and three compatible bonds.
    for (int atom = 0; atom < 4; ++atom)
        assert(smsd::detail::atomsCompatFast(query.molecule, atom, target.molecule, atom, chemistry));
    assert(smsd::detail::labelFrequencyUpperBound(query.molecule, target.molecule, chemistry) == 4);
    assert(smsd::detail::labelFrequencyUpperBoundDirected(query.molecule, target.molecule, chemistry) == 4);
    assert(smsd::findMCS(query.molecule, target.molecule, chemistry, options(false, true)).size() == 4);
}

void testSparseConnectedSelection() {
    auto query = makeGraph(2048, {{1500, 1501}, {1700, 1701}, {1701, 1702}});
    std::map<int, int> mapping{{1500, 1500}, {1501, 1501},
                              {1700, 1700}, {1701, 1701}, {1702, 1702}};
    auto result = smsd::detail::largestConnected(query.molecule, mapping, &query.molecule);
    assert(result.size() == 3 && result.count(1702));
    validate(query, query, result, false, true);

    auto target = makeGraph(2048, {{1500, 1501}, {1700, 1701}});
    auto common = smsd::detail::largestConnected(query.molecule, mapping, &target.molecule);
    assert(common.size() == 2 && common.count(1500));
    validate(query, target, common, false, true);
}

void testFlexibleAromaticExtension() {
    auto queryMolecule = smsd::parseSMILES("c1ccc(-c2ccccc2)cc1");
    auto targetMolecule = smsd::parseSMILES("c1ccc(Cc2ccccc2)cc1");
    auto toGraph = [](smsd::MolGraph molecule) {
        int atoms = molecule.n;
        std::vector<std::vector<int>> edges(atoms, std::vector<int>(atoms));
        for (int a = 0; a < atoms; ++a)
            for (int b = 0; b < atoms; ++b)
                edges[a][b] = molecule.bondOrder(a, b) != 0;
        return Graph{atoms, std::move(edges), std::move(molecule)};
    };
    auto query = toGraph(std::move(queryMolecule));
    auto target = toGraph(std::move(targetMolecule));
    auto checkChemistry = [&](const std::map<int, int>& mapping) {
        validate(query, target, mapping, false, true);
        for (const auto& [a, b] : mapping) {
            assert(query.molecule.atomicNum[a] == target.molecule.atomicNum[b]);
            assert(query.molecule.formalCharge[a] == target.molecule.formalCharge[b]);
            for (const auto& [k, l] : mapping) {
                int qOrder = query.molecule.bondOrder(a, k);
                if (!qOrder) continue;
                int tOrder = target.molecule.bondOrder(b, l);
                bool flexible = (query.molecule.bondAromatic(a, k) && (tOrder == 1 || tOrder == 2))
                    || (target.molecule.bondAromatic(b, l) && (qOrder == 1 || qOrder == 2));
                assert(qOrder == tOrder || flexible);
            }
        }
    };
    // A complete phenyl ring plus its attached query ring atom maps onto a
    // phenyl ring plus target methylene. FLEXIBLE atom aromaticity permits it.
    std::map<int, int> seven{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {10, 11}, {11, 12}};
    checkChemistry(seven);
    auto result = smsd::findMCS(query.molecule, target.molecule, {}, options(false, true));
    assert(result.size() >= seven.size());
    checkChemistry(result);
}

void testReverseCompleteRingFilter() {
    std::vector<int> queryElements(12, 6), targetElements(8, 6);
    queryElements[10] = queryElements[11] = 7;
    targetElements[6] = targetElements[7] = 7;
    auto query = makeGraph(12, {
        {0, 1}, {1, 2}, {2, 3}, {3, 8}, {8, 9}, {9, 0},
        {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {0, 10}, {10, 11}}, queryElements);
    auto target = makeGraph(8, {
        {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 6}, {6, 7}}, targetElements);
    smsd::ChemOptions chemistry;
    chemistry.completeRingsOnly = true;
    auto result = smsd::findMCS(query.molecule, target.molecule, chemistry, options(false, true));
    // The shared atoms of the fused query rings require all ten ring carbons.
    // The target has only six carbons, so the two pendant nitrogens are optimal.
    assert(result.size() == 2 && result.count(10) && result.count(11));
    assert(result.at(10) >= 6 && result.at(11) >= 6);
    validate(query, target, result, false, true);
}

void testExhaustiveSmallGraphs() {
    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 4; ++atoms)
        for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask)
            graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) {
        for (const auto& target : graphs) {
            for (bool induced : {false, true}) {
                for (bool connected : {false, true}) {
                    auto result = smsd::findMCS(
                        query.molecule, target.molecule, {}, options(induced, connected));
                    validate(query, target, result, induced, connected);
                    assert(static_cast<int>(result.size()) == maximumAtoms(query, target, induced, connected));
                    ++cases;
                }
            }
        }
    }
    assert(cases == 23104);
    std::cout << "Validated " << cases << " exhaustive small-graph MCS cases\n";
}

void testEnumeratedMappings() {
    auto isolates = graphFromMask(2, 0);
    auto path = graphFromMask(3, 5);
    auto result = smsd::findAllMCS(
        isolates.molecule, path.molecule, {}, options(true, false), 10);
    assert(!result.empty());
    for (const auto& mapping : result) {
        assert(mapping.size() == 2);
        validate(isolates, path, mapping, true, false);
    }

    std::vector<Graph> graphs;
    for (int atoms = 0; atoms <= 3; ++atoms)
        for (int mask = 0; mask < (1 << (atoms * (atoms - 1) / 2)); ++mask)
            graphs.push_back(graphFromMask(atoms, mask));
    size_t cases = 0;
    for (const auto& query : graphs) {
        for (const auto& target : graphs) {
            for (bool induced : {false, true}) {
                for (bool connected : {false, true}) {
                    auto mappings = smsd::findAllMCS(
                        query.molecule, target.molecule, {}, options(induced, connected), 3);
                    int expected = maximumAtoms(query, target, induced, connected);
                    assert(expected == 0 || !mappings.empty());
                    for (const auto& mapping : mappings) {
                        assert(static_cast<int>(mapping.size()) == expected);
                        validate(query, target, mapping, induced, connected);
                    }
                    ++cases;
                }
            }
        }
    }
    assert(cases == 576);
    std::cout << "Validated " << cases << " exhaustive enumeration cases\n";
}

} // namespace

int main() {
    testOptionRegressions();
    testWeightRegressions();
    testTautomerBound();
    testSparseConnectedSelection();
    testFlexibleAromaticExtension();
    testReverseCompleteRingFilter();
    testExhaustiveSmallGraphs();
    testEnumeratedMappings();
    std::cout << "MCS regressions passed\n";
}

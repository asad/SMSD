/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/query_smarts.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/substructure.hpp"
#include <exception>
#include <numeric>
#include <thread>
#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#endif
#include <set>
#include <tuple>

using smsd_test::require;
using Edge = std::tuple<int, int, int>;

static std::set<Edge> labelledEdges(const smsd::MolGraph& graph) {
    std::set<Edge> edges;
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom])
            if (atom < neighbor) {
                int first = graph.atomClass[atom], second = graph.atomClass[neighbor];
                if (second < first) std::swap(first, second);
                edges.emplace(first, second, graph.bondOrder(atom, neighbor));
            }
    return edges;
}

static void exactRoundTrip(const smsd::MolGraph& graph, const std::string& output) {
    const auto parsed = smsd::parseSMILES(output);
    require(parsed.n == graph.n, "SMILES changed the atom count");
    std::set<int> expected, actual;
    for (int atom = 0; atom < graph.n; ++atom) {
        expected.insert(graph.atomClass[atom]);
        actual.insert(parsed.atomClass[atom]);
        require(graph.atomicNum[atom] == 6 && parsed.atomicNum[atom] == 6,
                "SMILES changed an element");
    }
    require(expected == actual && expected.size() == static_cast<size_t>(graph.n),
            "SMILES changed distinct atom classes");
    require(labelledEdges(graph) == labelledEdges(parsed), "SMILES changed indexed bonds");
}

static smsd::MolGraph carbonGraph(std::vector<std::vector<int>> neighbors) {
    const int count = static_cast<int>(neighbors.size());
    std::vector<int> classes(count);
    std::iota(classes.begin(), classes.end(), 1);
    return smsd::MolGraph::Builder().atomCount(count).atomicNumbers(std::vector<int>(count, 6))
        .atomClasses(std::move(classes)).setNeighbors(std::move(neighbors)).build(false);
}

static void addBond(std::vector<std::vector<int>>& neighbors, int first, int second) {
    neighbors[first].push_back(second);
    neighbors[second].push_back(first);
}

static smsd::MolGraph pendantRings(int count) {
    std::vector<std::vector<int>> neighbors(4 * count);
    for (int atom = 1; atom < count; ++atom) addBond(neighbors, atom - 1, atom);
    for (int ring = 0; ring < count; ++ring) {
        const int first = count + 3 * ring;
        addBond(neighbors, ring, first);
        addBond(neighbors, first, first + 1);
        addBond(neighbors, first + 1, first + 2);
        addBond(neighbors, first + 2, first);
    }
    return carbonGraph(std::move(neighbors));
}

/** Isolate traversal from the separately tested canonical labeller. */
static void useIdentityWriterOrder(smsd::MolGraph& graph) {
    graph.canonicalLabel.resize(graph.n);
    std::iota(graph.canonicalLabel.begin(), graph.canonicalLabel.end(), 0);
    graph.canonicalComputed_ = true;
}

static smsd::MolGraph parallelClosures(int count) {
    std::vector<std::vector<int>> neighbors(2 * count);
    for (int atom = 1; atom < 2 * count; ++atom) addBond(neighbors, atom - 1, atom);
    for (int ring = 0; ring < count; ++ring) addBond(neighbors, ring, ring + count);
    auto graph = carbonGraph(std::move(neighbors));
    useIdentityWriterOrder(graph);
    return graph;
}

static void manyRings() {
    const auto graph = pendantRings(101);
    exactRoundTrip(graph, graph.toCanonicalSmiles());
    smsd::Match match;
    for (int atom = 0; atom < graph.n; ++atom) {
        match.atomMapping[atom] = atom;
        for (int neighbor : graph.neighbors[atom])
            if (atom < neighbor) match.bondMapping.push_back({atom, neighbor, atom, neighbor});
    }
    exactRoundTrip(graph, match.toSmiles(graph));
}

static void ringLimit() {
    const auto graph = parallelClosures(99);
    const auto text = graph.toCanonicalSmiles();
    require(text.find("%99") != std::string::npos, "99-open boundary was not exercised");
    exactRoundTrip(graph, text);
}

static void ringOverflow() {
    const auto graph = parallelClosures(100);
    try {
        graph.toCanonicalSmiles();
    } catch (const std::invalid_argument& error) {
        require(std::string(error.what()) == "SMILES output needs more than 99 simultaneously open ring closures",
                "Wrong error for excessive simultaneously open ring closures");
        return;
    }
    throw std::runtime_error("Writer accepted more than 99 simultaneously open ring closures");
}

/** The written SMARTS has the graph's shape and matches the graph it was written from. */
static void requireSmartsOf(const smsd::MolGraph& graph, const std::string& text) {
    const auto query = smsd::parseSMARTS(text);
    size_t edges = 0;
    for (int atom = 0; atom < graph.n; ++atom) edges += graph.neighbors[atom].size();
    require(query.atomCount() == graph.n && query.bonds.size() == edges / 2, "SMARTS changed the graph shape");
    require(!smsd::Substructure::smarts(text).match(graph).matches.empty(), "SMARTS does not match its source graph");
}

static void smartsRingLimit() {
    const auto graph = parallelClosures(99);
    const auto text = smsd::toQuerySMARTS(graph);
    require(text.find("%99") != std::string::npos, "99-open SMARTS boundary was not exercised");
    requireSmartsOf(graph, text);
}

static void smartsRingOverflow() {
    const auto graph = parallelClosures(100);
    try {
        (void)smsd::toQuerySMARTS(graph);
    } catch (const std::invalid_argument& error) {
        require(std::string(error.what()) == "SMARTS output needs more than 99 simultaneously open ring closures",
                "Wrong SMARTS error for excessive simultaneously open ring closures");
        return;
    }
    throw std::runtime_error("SMARTS writer accepted more than 99 simultaneously open ring closures");
}

/** More than 99 closures in all, but never more than one open, are written with reused numbers. */
static void smartsManyRings() {
    const auto graph = pendantRings(101);
    const auto text = smsd::toQuerySMARTS(graph);
    require(text.find("%") == std::string::npos, "Closed SMARTS ring numbers were not reused");
    requireSmartsOf(graph, text);
}

/** The canonical SMARTS writer shares the SMILES traversal and refuses the same graphs. */
static void canonicalSmartsRingOverflow() {
    const auto graph = parallelClosures(100);
    const std::string suffix = " output needs more than 99 simultaneously open ring closures";
    try {
        (void)smsd::toSMARTS(graph);
    } catch (const std::invalid_argument& error) {
        const std::string text = error.what();
        require(text.size() > suffix.size() && text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0,
                "Wrong canonical SMARTS error for excessive simultaneously open ring closures");
        return;
    }
    throw std::runtime_error("Canonical SMARTS writer accepted more than 99 simultaneously open ring closures");
}

static void ordinary() {
    for (const char* input : {"C1CC1.C1CC1", "c1ccc2ccccc2c1", "C12CC3CC(C1)CC(C2)C3",
            "C[C@H](N)C(=O)O", "C1[C@H](F)CCC[C@@H]1F", "C/C=C/C", "C/C=C\\C",
            "CCO.[Na+]", "[13CH4]", "[NH4+]"}) {
        const auto graph = smsd::parseSMILES(input);
        const auto output = graph.toCanonicalSmiles();
        require(smsd::parseSMILES(output).getCanonicalHash() == graph.getCanonicalHash(),
                std::string("Ordinary SMILES changed molecule: ") + input);
        std::cout << "TRACE " << input << '\t' << output << '\n';
        require(!smsd::toSMARTS(graph).empty(), "Shared SMARTS traversal returned empty");
    }
    require(smsd::parseSMILES("C1CC1.C1CC1").toCanonicalSmiles() == "C1CC1.C1CC1",
            "Ring labels were not released between components");
    const auto graph = smsd::parseSMILES(std::string(128, 'C'));
    require(smsd::parseSMILES(graph.toCanonicalSmiles()).getCanonicalHash() == graph.getCanonicalHash(),
            "Ordinary chain failed public canonical round trip");
}

struct ThreadResult {
    const smsd::MolGraph* graph;
    std::string text;
    std::exception_ptr failure;
};

static void* writeInThread(void* argument) {
    auto& result = *static_cast<ThreadResult*>(argument);
    try { result.text = result.graph->toCanonicalSmiles(); }
    catch (...) { result.failure = std::current_exception(); }
    return nullptr;
}

static void longChain() {
    const int count = 9999;
    std::vector<std::vector<int>> neighbors(count);
    for (int atom = 1; atom < count; ++atom) addBond(neighbors, atom - 1, atom);
    auto graph = carbonGraph(std::move(neighbors));
    graph.atomClass.assign(count, 0);
    useIdentityWriterOrder(graph);
    ThreadResult result{&graph, {}, {}};
#if defined(__unix__) || defined(__APPLE__)
    pthread_attr_t attributes;
    require(pthread_attr_init(&attributes) == 0, "Cannot initialise thread attributes");
    const int stackResult = pthread_attr_setstacksize(&attributes, 512 * 1024);
    if (stackResult != 0) {
        pthread_attr_destroy(&attributes);
        throw std::runtime_error("Cannot set the writer thread stack size");
    }
    pthread_t thread;
    const int createResult = pthread_create(&thread, &attributes, writeInThread, &result);
    pthread_attr_destroy(&attributes);
    require(createResult == 0, "Cannot create the writer thread");
    require(pthread_join(thread, nullptr) == 0, "Cannot join the writer thread");
#else
    std::thread thread(writeInThread, &result);
    thread.join();
#endif
    if (result.failure) std::rethrow_exception(result.failure);
    require(result.text == std::string(count, 'C'), "Long-chain writer changed atoms or branches");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {{"many-rings", manyRings}, {"ring-limit", ringLimit},
        {"ring-overflow", ringOverflow}, {"smarts-ring-limit", smartsRingLimit},
        {"smarts-ring-overflow", smartsRingOverflow}, {"smarts-many-rings", smartsManyRings},
        {"canonical-smarts-ring-overflow", canonicalSmartsRingOverflow}, {"ordinary", ordinary},
        {"long-chain", longChain}});
}

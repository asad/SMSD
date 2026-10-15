/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/smiles_parser.hpp"
#include "smsd/ring_finder.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

static void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

static int aromaticAtoms(const smsd::MolGraph& graph) {
    return std::accumulate(graph.aromatic.begin(), graph.aromatic.end(), 0);
}

static int aromaticBonds(const smsd::MolGraph& graph) {
    int count = 0;
    for (int atom = 0; atom < graph.n; ++atom)
        for (int neighbor : graph.neighbors[atom])
            if (neighbor > atom && graph.bondAromatic(atom, neighbor)) ++count;
    return count;
}

static void testChalcogenCations() {
    for (const char* text : {"c1cc[o+]cc1", "c1cc[s+]cc1", "c1cc[se+]cc1", "c1cc[te+]cc1",
            "c1ccc2cc[o+]cc2c1", "Cc1cc(C)[o+]c(C)c1", "c1ccoc1", "c1ccsc1", "c1cocn1"}) {
        auto graph = smsd::parseSMILES(text);
        const auto hash = graph.getCanonicalHash();
        require(graph.kekulize(), std::string("valid aromatic heterocycle failed Kekule assignment: ") + text);
        require(aromaticAtoms(graph) == 0 && aromaticBonds(graph) == 0, "Kekule conversion retained aromatic flags");
        graph.perceiveAromaticity();
        require(graph.getCanonicalHash() == hash, std::string("Kekule round-trip changed heterocycle: ") + text);
    }
}

static void testAtomicFailure() {
    auto graph = smsd::parseSMILES("c1ccccc1.c1cccc1");
    graph.prewarm();
    const auto before = graph;
    require(!graph.kekulize(), "odd neutral aromatic ring should fail Kekule assignment");
    require(graph.aromatic == before.aromatic && graph.bondPacked == before.bondPacked
            && graph.label == before.label && graph.canonicalLabel == before.canonicalLabel
            && graph.cachedNLF1 == before.cachedNLF1 && graph.cacheNonce == before.cacheNonce
            && graph.getCanonicalHash() == before.getCanonicalHash(),
            "failed Kekule conversion changed the graph or its cached properties");
}

static void testSaturatedCations() {
    for (const char* text : {"C1=CC=C[NH2+]C=C1", "C[N+]1(C)C=CC=CC=C1",
            "C1=CC=C[S+](C)C=C1", "C1=CC=C[O+](C)C=C1"}) {
        const auto graph = smsd::parseSMILES(text);
        require(aromaticAtoms(graph) == 0 && aromaticBonds(graph) == 0,
                std::string("saturated ring cation was marked aromatic: ") + text);
    }
    require(aromaticAtoms(smsd::parseSMILES("C1=CC=[O+]C=C1")) == 6, "pyrylium should remain aromatic");
    require(aromaticAtoms(smsd::parseSMILES("C1=CC=C[CH+]C=C1")) == 7, "tropylium should remain aromatic");
}

static void testAromaticForms() {
    struct Example { const char* aromatic; const char* kekule; int atoms, bonds; };
    for (const auto& example : {
            Example{"c1ccc2cc3cc4cc5cc6ccccc6cc5cc4cc3cc2c1", "C1=CC=C2C=C3C=C4C=C5C=C6C=CC=CC6=CC5=CC4=CC3=CC2=C1", 26, 31},
            {"c1ccc2cc3cc4cc5cc6cc7ccccc7cc6cc5cc4cc3cc2c1", "C1=CC=C2C=C3C=C4C=C5C=C6C=C7C=CC=CC7=CC6=CC5=CC4=CC3=CC2=C1", 30, 36},
            {"c1cc2cccc3c4cccc5cccc(c(c1)c23)c54", "C1=CC2=CC=CC3=C2C(=C1)C1=CC=CC2=C1C3=CC=C2", 20, 24},
            {"c1ccc2c(c1)c1cccc3c4c5ccccc5c5cccc(c2c13)c54", "C1=CC=C2C(=C1)C1=CC=CC3=C1C2=C1C=CC=C2C4=CC=CC=C4C3=C21", 26, 32},
            {"c1ccc2cc3c(cc2c1)CCCCCCCCCCCCCCC3", "C1=CC=C2C=C3CCCCCCCCCCCCCCCC3=CC2=C1", 10, 11},
            {"CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)n1c3ccccc3c3c4c(c5c6ccccc6n2c5c31)C(=O)NC4", "CN[C@@H]1C[C@H]2O[C@@](C)([C@@H]1OC)N1C3=CC=CC=C3C3=C4CNC(=O)C4=C4C5=CC=CC=C5N2C4=C31", 20, 24}}) {
        const auto aromatic = smsd::parseSMILES(example.aromatic);
        auto kekule = smsd::parseSMILES(example.kekule);
        require(aromaticAtoms(kekule) == example.atoms && aromaticBonds(kekule) == example.bonds,
                std::string("Kekule input has inconsistent aromaticity: ") + example.kekule);
        require(kekule.getCanonicalHash() == aromatic.getCanonicalHash(), "aromatic input forms have different identities");
        require(kekule.kekulize(), "fused aromatic system failed Kekule assignment");
        kekule.perceiveAromaticity();
        require(kekule.getCanonicalHash() == aromatic.getCanonicalHash(), "reperception changed fused aromatic system");
    }
}

static void testRelevantCycles() {
    struct Example { const char* text; size_t count; };
    for (const auto& example : {Example{"C1CCCCC1", 1}, {"c1ccc2ccccc2c1", 2},
            {"C12C3C4C1C5C2C3C45", 6}, {"C1CC2CC(C1)C2", 3}}) {
        const auto graph = smsd::parseSMILES(example.text);
        const auto cycles = graph.computeRelevantCycles();
        require(cycles.size() == example.count, "incorrect number of member relevant cycles");
        for (const auto& cycle : cycles) {
            require(std::set<int>(cycle.begin(), cycle.end()).size() == cycle.size(), "cycle repeats an atom");
            for (size_t atom = 0; atom < cycle.size(); ++atom)
                require(graph.bondOrder(cycle[atom], cycle[(atom + 1) % cycle.size()]) != 0,
                        "cycle contains a missing bond");
        }
    }
}

static void testUniqueRingFamilies() {
    struct Example { const char* text; size_t count; };
    for (const auto& example : {Example{"C1CCCCC1", 1}, {"c1ccc2ccccc2c1", 2},
            {"C12C3C4C1C5C2C3C45", 6}, {"c1ccccc1.c1ccccc1", 2},
            {"C1=CCCCC1.C1CCCCC1", 2}, {"C1CC2CC(C1)C2", 2}}) {
        auto graph = smsd::parseSMILES(example.text);
        const auto member = graph.computeURFs();
        const auto free = smsd::computeURFs(graph);
        require(member.size() == example.count && free.size() == example.count,
                std::string("incorrect unique ring family grouping: ") + example.text);
        if (std::string(example.text) == "C1CC2CC(C1)C2") {
            require(std::any_of(member.begin(), member.end(), [](const auto& family) { return family.size() == 2; })
                    && std::any_of(free.begin(), free.end(), [](const auto& family) { return family.size() == 2; }),
                    "rings exchangeable through a shorter ring should share a family");
        }
    }
}

int main(int argc, char** argv) {
    int failures = 0;
    const std::pair<const char*, void (*)()> tests[] = {
        {"chalcogens", testChalcogenCations}, {"atomic-kekule", testAtomicFailure},
        {"saturated-cations", testSaturatedCations}, {"aromatic-forms", testAromaticForms},
        {"cycles", testRelevantCycles}, {"families", testUniqueRingFamilies}
    };
    for (const auto& test : tests) {
        if (argc > 1 && argv[1] != std::string(test.first)) continue;
        try { test.second(); std::cout << "PASS " << test.first << '\n'; }
        catch (const std::exception& error) {
            ++failures;
            std::cerr << "FAIL " << test.first << ": " << error.what() << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}

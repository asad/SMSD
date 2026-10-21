/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/search.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

using smsd_test::require;

static void testEquivalentSmiles() {
    const std::pair<const char*, const char*> examples[] = {
        {"CC(=O)O", "CC(O)=O"}, {"C=CC", "CC=C"},
        {"[13CH3]C", "C[13CH3]"}, {"C[NH2+]CC", "CC[NH2+]C"},
        {"[C@H](F)(Cl)Br", "[C@@H](Cl)(F)Br"},
        {"C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O"},
        {"Cl/C=C/Cl", "C(=C/Cl)\\Cl"},
        {"C[C@H](F)[C@H](F)C", "C[C@@H](F)[C@@H](F)C"},
        {"C/C=C/C=C\\C", "C/C=C\\C=C\\C"},
        {"C/C=C(/F)\\C(Cl)=C/C", "FC(/C(=C\\C)Cl)=C/C"},
        {"C[C@H](F)Cl.C[C@@H](F)Cl", "C[C@@H](F)Cl.C[C@H](F)Cl"},
        {"C1[C@H](F)CCC[C@@H]1F", "F[C@H]1CCC[C@@H](F)C1"}
    };
    for (const auto& example : examples) {
        const auto first = smsd::parseSMILES(example.first);
        const auto second = smsd::parseSMILES(example.second);
        require(first.getCanonicalHash() == second.getCanonicalHash(),
                std::string("equivalent SMILES have different hashes: ") + example.first);
        require(smsd::search_detail::sameCanonicalGraph(first, second), "equivalent molecules must have the same canonical graph");
        require(smsd::toSMILES(first) == smsd::toSMILES(second),
                std::string("equivalent SMILES have different canonical strings: ") + example.first);
    }
    const auto clockwise = smsd::parseSMILES("C[C@H](N)C(=O)O");
    const auto anticlockwise = smsd::parseSMILES("C[C@@H](N)C(=O)O");
    require(clockwise.getCanonicalHash() != anticlockwise.getCanonicalHash(), "enantiomers must remain distinct");
    const auto trans = smsd::parseSMILES("Cl/C=C/Cl");
    const auto cis = smsd::parseSMILES("Cl/C=C\\Cl");
    require(trans.getCanonicalHash() != cis.getCanonicalHash(), "alkene configurations must remain distinct");
}

static smsd::MolGraph reordered(const smsd::MolGraph& graph, std::mt19937& random) {
    std::vector<int> permutation(graph.n), inverse(graph.n);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::shuffle(permutation.begin(), permutation.end(), random);
    for (int atom = 0; atom < graph.n; ++atom) inverse[permutation[atom]] = atom;
    std::vector<int> numbers(graph.n), charges(graph.n), masses(graph.n), hydrogens(graph.n), classes(graph.n);
    std::vector<uint8_t> aromatic(graph.n), rings(graph.n);
    std::vector<std::vector<int>> neighbors(graph.n), orders(graph.n);
    std::vector<std::vector<bool>> aromaticBonds(graph.n), ringBonds(graph.n);
    for (int atom = 0; atom < graph.n; ++atom) {
        const int original = permutation[atom];
        numbers[atom] = graph.atomicNum[original];
        charges[atom] = graph.formalCharge[original];
        masses[atom] = graph.massNumber[original];
        hydrogens[atom] = graph.hydrogenCount[original];
        classes[atom] = graph.atomClass[original];
        aromatic[atom] = graph.aromatic[original];
        rings[atom] = graph.ring[original];
        auto adjacent = graph.neighbors[original];
        std::shuffle(adjacent.begin(), adjacent.end(), random);
        for (int neighbor : adjacent) {
            neighbors[atom].push_back(inverse[neighbor]);
            orders[atom].push_back(graph.bondOrder(original, neighbor));
            aromaticBonds[atom].push_back(graph.bondAromatic(original, neighbor));
            ringBonds[atom].push_back(graph.bondInRing(original, neighbor));
        }
    }
    return smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(std::move(numbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses))
        .hydrogenCounts(std::move(hydrogens)).atomClasses(std::move(classes))
        .aromaticFlags(std::move(aromatic)).ringFlags(std::move(rings))
        .setNeighbors(std::move(neighbors)).setBondOrders(std::move(orders))
        .bondAromaticFlags(std::move(aromaticBonds)).bondRingFlags(std::move(ringBonds)).build(false);
}

static void testAtomAndBondOrder() {
    std::mt19937 random(1234);
    for (const char* smiles : {"CC(=O)O", "C=CC", "[13CH3]C", "C[NH2+]CC",
            "[CH2]CC", "[CH3:1][CH2:2]O", "c1ccccc1", "c1ccc2ccccc2c1", "C1CC1.C1CCCCC1"}) {
        const auto original = smsd::parseSMILES(smiles);
        const auto expectedHash = original.getCanonicalHash();
        const auto expectedSmiles = smsd::toSMILES(original);
        for (int trial = 0; trial < 12; ++trial) {
            const auto copy = reordered(original, random);
            require(copy.getCanonicalHash() == expectedHash, std::string("atom order changed hash: ") + smiles);
            require(smsd::toSMILES(copy) == expectedSmiles, std::string("atom order changed canonical SMILES: ") + smiles);
        }
    }
}

static void testVerifiedOrbits() {
    const auto propene = smsd::parseSMILES("C=CC");
    const auto& propeneOrbits = propene.getOrbit();
    require(propeneOrbits[0] != propeneOrbits[2], "different terminal carbons must have different orbits");
    const auto mixed = smsd::parseSMILES("C1CC1.C1CCCCC1");
    const auto& orbits = mixed.getOrbit();
    for (int atom = 0; atom < mixed.n; ++atom)
        require(orbits[atom] == (atom < 3 ? orbits[0] : orbits[3]), "equivalent ring atoms should share an orbit");
    require(orbits[0] != orbits[3], "different ring components must have different orbits");
    require(!mixed.automorphismGeneratorsTruncated(), "small ring symmetry search should complete");
    const auto benzene = smsd::parseSMILES("c1ccccc1");
    std::set<std::vector<int>> group;
    std::vector<int> identity(benzene.n);
    std::iota(identity.begin(), identity.end(), 0);
    group.insert(identity);
    std::vector<std::vector<int>> pending{identity};
    for (size_t position = 0; position < pending.size(); ++position) {
        for (const auto& generator : benzene.getAutomorphismGenerators()) {
            std::vector<int> product(benzene.n);
            for (int atom = 0; atom < benzene.n; ++atom) product[atom] = generator[pending[position][atom]];
            if (group.insert(product).second) pending.push_back(std::move(product));
        }
    }
    require(group.size() == 12 && !benzene.automorphismGeneratorsTruncated(), "benzene should retain its full symmetry group");
}

static void testMacrocycles() {
    for (int ringSize : {21, 22, 33}) {
        const auto graph = smsd::parseSMILES("C1" + std::string(ringSize - 1, 'C') + "1");
        require(std::count(graph.ring.begin(), graph.ring.end(), uint8_t(1)) == ringSize,
                "macrocycle atoms must be ring atoms");
        for (int atom = 0; atom < graph.n; ++atom)
            for (int neighbor : graph.neighbors[atom])
                require(graph.bondInRing(atom, neighbor), "macrocycle bonds must be ring bonds");
        const auto& memberRings = graph.computeRings();
        const auto sssr = smsd::computeSSSR(graph);
        const auto relevant = smsd::computeRelevantCycles(graph);
        require(memberRings.size() == 1 && memberRings[0].size() == static_cast<size_t>(ringSize), "member ring basis dropped macrocycle");
        require(sssr.size() == 1 && sssr[0].size() == static_cast<size_t>(ringSize), "SSSR dropped macrocycle");
        require(relevant.size() == 1 && relevant[0].size() == static_cast<size_t>(ringSize), "relevant cycles dropped macrocycle");
        require(graph.ringSystemCount() == 1, "macrocycle must form one ring system");
    }
    const auto peptide = smsd::parseSMILES("O=C1NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NC1");
    require(std::count(peptide.ring.begin(), peptide.ring.end(), uint8_t(1)) == 21,
            "cyclic heptaglycine must retain its 21 ring atoms");
    const auto fused = smsd::parseSMILES("c1ccc2c(c1)" + std::string(21, 'C') + "2");
    require(std::count(fused.ring.begin(), fused.ring.end(), uint8_t(1)) == fused.n,
            "fused benzene and macrocycle must retain every ring atom");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"equivalent-smiles", testEquivalentSmiles}, {"atom-and-bond-order", testAtomAndBondOrder},
        {"verified-orbits", testVerifiedOrbits}, {"macrocycles", testMacrocycles}});
}

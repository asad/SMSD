/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/ring_finder.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/cip.hpp"
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"

#include <algorithm>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

using smsd_test::require;

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

static const char* stereoMacrocycle = "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1";

static smsd::MolGraph withAromaticStereo(const smsd::MolGraph& graph,
                                       const std::vector<std::tuple<int, int, int>>& stereo,
                                       bool reverseAtomOrder = false) {
    auto index = [&](int atom) { return reverseAtomOrder ? graph.n - atom - 1 : atom; };
    std::vector<std::vector<int>> neighbors(graph.n), orders(graph.n);
    std::vector<std::vector<bool>> aromatic(graph.n), rings(graph.n);
    std::vector<int> numbers(graph.n), charges(graph.n), hydrogens(graph.n);
    std::vector<uint8_t> aromaticAtoms(graph.n), ringAtoms(graph.n);
    for (int original = 0; original < graph.n; ++original) {
        const int atom = index(original);
        numbers[atom] = graph.atomicNum[original]; charges[atom] = graph.formalCharge[original];
        hydrogens[atom] = graph.hydrogenCount[original];
        aromaticAtoms[atom] = graph.aromatic[original]; ringAtoms[atom] = graph.ring[original];
        for (int neighbor : graph.neighbors[original]) {
            neighbors[atom].push_back(index(neighbor));
            orders[atom].push_back(graph.bondOrder(original, neighbor));
            aromatic[atom].push_back(graph.bondAromatic(original, neighbor));
            rings[atom].push_back(graph.bondInRing(original, neighbor));
        }
    }
    std::vector<std::tuple<int, int, int>> mappedStereo;
    for (const auto& focus : stereo)
        mappedStereo.emplace_back(index(std::get<0>(focus)), index(std::get<1>(focus)), std::get<2>(focus));
    return smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(numbers)
        .formalCharges(charges).hydrogenCounts(hydrogens)
        .aromaticFlags(aromaticAtoms).ringFlags(ringAtoms).setNeighbors(neighbors)
        .setBondOrders(orders).bondAromaticFlags(aromatic).bondRingFlags(rings)
        .doubleBondStereoSparse(mappedStereo).build(false);
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

// The deprecated neighbour-label table stays under test, so its warning is silenced here alone.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
static const std::vector<std::vector<int>>& neighbourLabels(const smsd::MolGraph& graph) {
    return graph.getNLF1();
}
#if defined(_MSC_VER)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif

static void testAtomicFailure() {
    auto graph = smsd::parseSMILES("c1ccccc1.c1cccc1");
    graph.prewarm();
    neighbourLabels(graph);
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

static void testFusedResonanceForms() {
    const std::pair<const char*, const char*> examples[] = {
        {"C1=CC=C2NC=CC3=NC=CC1=C23", "C1=CC=C2NC=CC3=C2C1=CC=N3"},
        {"C1=CC=C2NC3=CC=CC=C3C3=CC=NC1=C23", "C1=CC=C2NC3=CC=CC=C3C3=C2C1=NC=C3"}
    };
    for (const auto& example : examples) {
        const auto first = smsd::parseSMILES(example.first);
        const auto second = smsd::parseSMILES(example.second);
        require(aromaticAtoms(first) == aromaticAtoms(second),
                "Equivalent fused-ring resonance forms have different aromatic atom counts");
        require(aromaticBonds(first) == aromaticBonds(second),
                "Equivalent fused-ring resonance forms have different aromatic bond counts");
        require(first.getCanonicalHash() == second.getCanonicalHash(),
                "Equivalent fused-ring resonance forms have different canonical identities");
    }
    require(aromaticAtoms(smsd::parseSMILES(examples[0].first)) == 10,
            "The donor ring must remain non-aromatic in either resonance form");
}

static void testExternalCarbonyls() {
    require(aromaticAtoms(smsd::parseSMILES("O=C1C=CC=CC1")) == 0,
            "A carbonyl must not restore conjugation to a saturated ring atom");
    require(aromaticAtoms(smsd::parseSMILES("O=C1C=CC=CC=C1")) == 7,
            "An external carbonyl must retain the zero-electron contribution of tropone");
}

static void testExternalCarbonPi() {
    for (const char* text : {
            "C=C1C=NC2=CNN=NC2=N1",
            "C1=CC(=O)C(=C2C=CN=C3NC(=N)NC=C23)C(=O)O1"}) {
        auto graph = smsd::parseSMILES(text);
        require(aromaticAtoms(graph) == 6,
                "An external carbon double bond must contribute one electron");
        const auto hash = graph.getCanonicalHash();
        require(graph.kekulize(), "An external carbon double bond prevented Kekule assignment");
        for (int atom = 0; atom < graph.n; ++atom) {
            if (graph.atomicNum[atom] != 6) continue;
            int valence = graph.hydrogenCount[atom];
            for (int neighbor : graph.neighbors[atom]) valence += graph.bondOrder(atom, neighbor);
            require(valence == 4, "Kekule assignment exceeded neutral carbon valence");
        }
        graph.perceiveAromaticity();
        require(graph.getCanonicalHash() == hash, "External carbon pi bond changed on reperception");
        require(smsd::parseSMILES(smsd::toSMILES(graph)).getCanonicalHash() == hash,
                "External carbon pi bond changed in aromatic SMILES");
    }
}

static void testAromaticDonorReperception() {
    const std::pair<const char*, int> examples[] = {
        {"C1=CC2=C3C=CC(=CC4=[O+]C(=C5C=CC(=CC1=[O+]2)O5)C=C4)O3", 20},
        {"CCC1C2=CC3=NC4=C(C(=O)OC(=O)C4=C3C)C3=NC(=CC4=NC(=CC(=N2)C1C)C(C(C)=O)=C4C)C(C)C3CCC(=O)OC", 22}
    };
    for (const auto& example : examples) {
        auto graph = smsd::parseSMILES(example.first);
        require(aromaticAtoms(graph) == example.second, "Unexpected initial oxygen ring aromaticity");
        const auto hash = graph.getCanonicalHash();
        graph.perceiveAromaticity();
        require(graph.getCanonicalHash() == hash,
                "Reperception changed aromatic donor or carbonyl donation");
        require(smsd::parseSMILES(smsd::toSMILES(graph)).getCanonicalHash() == hash,
                "Aromatic donor or carbonyl donation changed in SMILES round-trip");
    }
}

static void testAromaticStereoNormalisation() {
    const std::string block = R"mol(Benzene
  SMSD              2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
)mol";
    const auto benzene = smsd::readMolBlock(block);
    require(aromaticAtoms(benzene) == 6, "MOL benzene must be aromatic");
    for (int atom = 0; atom < benzene.n; ++atom)
        for (int neighbor : benzene.neighbors[atom])
            require(benzene.dbStereo(atom, neighbor) == 0,
                    "Aromatic bonds must not retain geometry-derived alkene stereo");
    require(benzene.getCanonicalHash() == smsd::parseSMILES("c1ccccc1").getCanonicalHash(),
            "Benzene MOL and SMILES must have the same canonical identity");

    for (int extraAtoms : {0, smsd::MolGraph::SPARSE_THRESHOLD}) {
        const std::string padding = extraAtoms == 0 ? "" : "." + std::string(extraAtoms, 'C');
        const auto graph = smsd::parseSMILES(std::string(stereoMacrocycle) + ".F/C=C/F" + padding);
        require(graph.useDense == (extraAtoms == 0), "Expected both stereo storage modes");
        require(graph.dbStereo(4, 5) != 0 && graph.dbStereo(8, 9) != 0 && graph.dbStereo(16, 17) != 0,
                "Aromatic perception discarded specified macrocycle stereo");
        require(graph.dbStereo(28, 29) != 0, "Genuine external alkene stereo must survive aromatic perception");
    }
}

static void testAromaticStereoBuilder() {
    for (int atoms : {10, smsd::MolGraph::SPARSE_THRESHOLD + 1}) {
        std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
        std::vector<std::vector<bool>> aromaticBonds(atoms);
        auto bond = [&](int first, int second, int order) {
            neighbors[first].push_back(second); orders[first].push_back(order);
            aromaticBonds[first].push_back(order == 4);
            neighbors[second].push_back(first); orders[second].push_back(order);
            aromaticBonds[second].push_back(order == 4);
        };
        for (int atom = 0; atom < 10; ++atom) bond(atom, (atom + 1) % 10, 4);
        for (int atom = 11; atom < atoms; ++atom) bond(atom - 1, atom, 1);
        std::vector<uint8_t> aromatic(atoms, 0);
        std::fill(aromatic.begin(), aromatic.begin() + 10, uint8_t(1));
        std::vector<int> hydrogens(atoms, 1);
        for (int atom = 10; atom < atoms; ++atom)
            hydrogens[atom] = 4 - static_cast<int>(neighbors[atom].size());
        for (bool denseStereo : {false, true}) {
            auto builder = smsd::MolGraph::Builder().atomCount(atoms)
                .atomicNumbers(std::vector<int>(atoms, 6)).hydrogenCounts(hydrogens)
                .aromaticFlags(aromatic).setNeighbors(neighbors).setBondOrders(orders)
                .bondAromaticFlags(aromaticBonds);
            if (denseStereo) {
                std::vector<std::vector<int>> stereo(atoms, std::vector<int>(atoms, 0));
                stereo[0][1] = stereo[1][0] = 1;
                builder.doubleBondStereo(stereo);
            } else {
                builder.doubleBondStereoSparse({{0, 1, 1}});
            }
            const auto graph = builder.build(false);
            require(graph.dbStereo(0, 1) == 1, "Builder discarded annotated aromatic focus");
        }
    }
}

static void testStereoConstrainedKekule() {
    for (const auto& example : {
            std::make_pair("c1ccccccccc1", std::make_pair(0, 9)),
            std::make_pair("c1ccc2cccc-2cc1", std::make_pair(0, 1))}) {
        auto graph = withAromaticStereo(smsd::parseSMILES(example.first),
                                       {{example.second.first, example.second.second, 1}});
        const auto hash = graph.getCanonicalHash();
        require(graph.kekulize(), "A compatible aromatic stereo focus prevented Kekule assignment");
        require(graph.bondOrder(example.second.first, example.second.second) == 2,
                "Kekule assignment moved an annotated focus to a single bond");
        require(graph.dbStereo(example.second.first, example.second.second) == 1,
                "Kekule assignment changed the stereo reference relation");
        graph.perceiveAromaticity();
        require(graph.getCanonicalHash() == hash, "Stereo-constrained Kekule round-trip changed the graph");
        auto reversed = withAromaticStereo(smsd::parseSMILES(example.first),
                                          {{example.second.first, example.second.second, 1}}, true);
        require(reversed.getCanonicalHash() == hash, "Atom reordering changed aromatic stereo identity");
        require(reversed.kekulize(), "Atom reordering prevented constrained Kekule assignment");
        reversed.perceiveAromaticity();
        require(reversed.getCanonicalHash() == hash,
                "Constrained Kekule assignment changed atom-order-invariant stereo identity");
    }
    auto conflicting = withAromaticStereo(smsd::parseSMILES("c1ccccc1"), {{0, 1, 1}, {1, 2, 2}});
    conflicting.prewarm();
    const auto before = conflicting;
    require(!conflicting.kekulize(), "Adjacent incompatible aromatic stereo foci must fail Kekule assignment");
    require(conflicting.bondPacked == before.bondPacked && conflicting.aromatic == before.aromatic
            && conflicting.dbStereoConf == before.dbStereoConf && conflicting.label == before.label
            && conflicting.cacheNonce == before.cacheNonce
            && conflicting.getCanonicalHash() == before.getCanonicalHash(),
            "Failed constrained Kekule assignment changed the graph");
}

static void testAromaticStereoDescriptors() {
    const auto graph = smsd::parseSMILES(stereoMacrocycle);
    for (const auto& focus : {std::make_pair(4, 5), std::make_pair(8, 9), std::make_pair(16, 17)})
        require(smsd::cip::assignEZ(graph, focus.first, focus.second) == smsd::cip::EZLabel::Z,
                "A specified aromatic macrocycle focus lost its known Z descriptor");
    const auto all = smsd::cip::assignAll(graph);
    for (const auto& focus : {std::make_pair(4, 5), std::make_pair(8, 9), std::make_pair(16, 17)})
        require(std::find(all.ezBonds.begin(), all.ezBonds.end(),
                          std::make_tuple(focus.first, focus.second, smsd::cip::EZLabel::Z)) != all.ezBonds.end(),
                "Bulk descriptor assignment skipped a known aromatic macrocycle Z focus");
    require(smsd::cip::assignAll(smsd::parseSMILES("c1ccccc1")).ezBonds.empty(),
            "Unannotated aromatic bonds must not acquire E/Z descriptors");
    const auto alkene = smsd::parseSMILES("F/C=C/F");
    require(smsd::cip::assignEZ(alkene, 1, 2) == smsd::cip::EZLabel::E,
            "Ordinary alkene descriptor changed");

    smsd::SubstructureOptions options;
    options.chemistry.useBondStereo = true;
    const auto self = smsd::Substructure::find(graph).withOptions(options).match(graph);
    require(!self.matches.empty(), "Annotated macrocycle must match itself");
    require(self.matches.front().toGraph(graph).getCanonicalHash() == graph.getCanonicalHash(),
            "Mapped graph extraction discarded aromatic focus stereo");
    const auto opposite = smsd::parseSMILES(
        "C1=CC2=N/C1=C/C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1");
    const auto mismatch = smsd::Substructure::find(graph).withOptions(options).match(opposite);
    require(mismatch.matches.empty() && mismatch.complete,
            "Mapped stereo must distinguish opposite macrocycle configurations");
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

// The free ring perception is the canonical one; the graph's bounded shortest-path basis
// must agree with it on the relevant cycles and on the sizes of the basis rings.
static void testRingPerceptionParity() {
    auto asSets = [](const std::vector<std::vector<int>>& rings) {
        std::set<std::set<int>> sets;
        for (const auto& ring : rings) sets.emplace(ring.begin(), ring.end());
        return sets;
    };
    auto sizes = [](const std::vector<std::vector<int>>& rings) {
        std::vector<size_t> result;
        for (const auto& ring : rings) result.push_back(ring.size());
        std::sort(result.begin(), result.end());
        return result;
    };
    for (const char* text : {"C12C3C4C1C5C4C3C25", "C1C2CC3CC1CC(C2)C3", "C1CC2CCC1C2", "c1ccc2ccccc2c1",
            "C1CCC2(CC1)CCCC2", "c1ccc2c(c1)c1ccccc1c2", "C1CC2CCC3CCCC4CCCC(C1)C2C34",
            "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C", "c1cc2cccc3c4cccc5cccc(c(c1)c23)c54",
            "O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1", "C1CCCCC1.c1ccccc1"}) {
        const auto graph = smsd::parseSMILES(text);
        require(asSets(smsd::computeRelevantCycles(graph)) == asSets(graph.computeRelevantCycles()),
                std::string("relevant cycles differ between the free and member forms: ") + text);
        require(sizes(smsd::computeSSSR(graph)) == sizes(graph.computeRings()),
                std::string("basis ring sizes differ between the free and member forms: ") + text);
    }
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"chalcogens", testChalcogenCations}, {"atomic-kekule", testAtomicFailure},
        {"ring-parity", testRingPerceptionParity},
        {"saturated-cations", testSaturatedCations}, {"aromatic-forms", testAromaticForms},
        {"fused-resonance", testFusedResonanceForms},
        {"external-carbonyls", testExternalCarbonyls},
        {"external-carbon-pi", testExternalCarbonPi},
        {"donor-reperception", testAromaticDonorReperception},
        {"aromatic-stereo", testAromaticStereoNormalisation},
        {"aromatic-stereo-builder", testAromaticStereoBuilder},
        {"stereo-kekule", testStereoConstrainedKekule},
        {"aromatic-stereo-descriptors", testAromaticStereoDescriptors},
        {"cycles", testRelevantCycles}, {"families", testUniqueRingFamilies}});
}

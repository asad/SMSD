/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using smsd_test::require;

// The deprecated three-argument form is tested here, so its warning is silenced here alone.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
static double legacyConfidence(const smsd::MolGraph& query, const smsd::MolGraph& target,
                               const std::unordered_map<int, int>& atoms) {
    return smsd::computeTautomerConfidence(query, target, atoms);
}
#if defined(_MSC_VER)
#pragma warning(pop)
#else
#pragma GCC diagnostic pop
#endif

static smsd::Match identitySelection(const smsd::MolGraph& query) {
    smsd::Match result;
    for (int atom = 0; atom < query.n; ++atom) {
        result.atomMapping.emplace(atom, atom);
        for (int neighbor : query.neighbors[atom])
            if (atom < neighbor) result.bondMapping.push_back({atom, neighbor, atom, neighbor});
    }
    return result;
}

// The atom-order test reads a benchmark corpus, whose path the build defines only when the corpus is present.
#ifdef SMSD_TAUTOMER_PERMUTATION_SMILES
/** The graph with its atoms, and the neighbours of each atom, in a random order; permutation[atom] is
 *  the original index of each new atom. */
static smsd::MolGraph reordered(const smsd::MolGraph& graph, std::mt19937& random, std::vector<int>& permutation) {
    std::vector<int> inverse(graph.n);
    permutation.resize(graph.n);
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
        hydrogens[atom] = graph.implicitHydrogenCount(original);
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

/** The lowest original atom index in each atom's class, by original atom, or -1 outside every class. */
static std::vector<int> classRepresentatives(const smsd::MolGraph& graph, const std::vector<int>& original) {
    std::map<int, int> lowest;
    for (int atom = 0; atom < graph.n; ++atom) {
        const int cls = graph.tautomerClass[atom];
        if (cls < 0) continue;
        const auto found = lowest.find(cls);
        if (found == lowest.end() || original[atom] < found->second) lowest[cls] = original[atom];
    }
    std::vector<int> representatives(graph.n, -1);
    for (int atom = 0; atom < graph.n; ++atom)
        if (graph.tautomerClass[atom] >= 0) representatives[original[atom]] = lowest[graph.tautomerClass[atom]];
    return representatives;
}

/**
 * Whether a symmetry of `graph` carries the class numbers `first` onto `second` (both by atom of
 * `graph`): an induced embedding of the graph labelled with `first` into the graph labelled with
 * `second`, each number written into the mass number of its atom.
 */
static bool symmetricNumbering(const smsd::MolGraph& graph, const std::vector<int>& first, const std::vector<int>& second) {
    auto labelled = [&](const std::vector<int>& classes) {
        std::vector<int> numbers(graph.n), charges(graph.n), masses(graph.n), hydrogens(graph.n);
        std::vector<uint8_t> aromatic(graph.n), rings(graph.n);
        std::vector<std::vector<int>> orders(graph.n);
        std::vector<std::vector<bool>> aromaticBonds(graph.n), ringBonds(graph.n);
        for (int atom = 0; atom < graph.n; ++atom) {
            numbers[atom] = graph.atomicNum[atom];
            charges[atom] = graph.formalCharge[atom];
            masses[atom] = graph.massNumber[atom] * 1000 + classes[atom] + 2;
            hydrogens[atom] = graph.implicitHydrogenCount(atom);
            aromatic[atom] = graph.aromatic[atom];
            rings[atom] = graph.ring[atom];
            for (int neighbor : graph.neighbors[atom]) {
                orders[atom].push_back(graph.bondOrder(atom, neighbor));
                aromaticBonds[atom].push_back(graph.bondAromatic(atom, neighbor));
                ringBonds[atom].push_back(graph.bondInRing(atom, neighbor));
            }
        }
        return smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(std::move(numbers))
            .formalCharges(std::move(charges)).massNumbers(std::move(masses)).hydrogenCounts(std::move(hydrogens))
            .aromaticFlags(std::move(aromatic)).ringFlags(std::move(rings)).setNeighbors(graph.neighbors)
            .setBondOrders(std::move(orders)).bondAromaticFlags(std::move(aromaticBonds))
            .bondRingFlags(std::move(ringBonds)).build(false);
    };
    smsd::SubstructureOptions options;
    options.chemistry.matchFormalCharge = options.chemistry.matchIsotope = true;
    options.chemistry.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    options.chemistry.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    options.chemistry.ringMatchesRingOnly = true;
    options.induced = true;
    options.timeoutMs = 60000;
    return !smsd::Substructure::find(labelled(first)).withOptions(options).match(labelled(second)).matches.empty();
}
#endif

/** The most specific motif sets each atom's weight and site, and each class takes one site model. */
static void motifPriorityAndClassModel() {
    using G = smsd::MolGraph;
    auto baseline = [](const char* smiles) { auto graph = smsd::parseSMILES(smiles); graph.computeTautomerClasses(); return graph; };
    const auto hydroxamic = baseline("ONC(C)=O"), hydroxamicReversed = baseline("CC(=O)NO");
    for (const auto& [graph, nitrogen, oxygen] : {std::make_tuple(&hydroxamic, 1, 4), std::make_tuple(&hydroxamicReversed, 3, 2)})
        require(graph->tautomerSite[nitrogen] == G::TAUTOMER_SITE_ACID && graph->tautomerPKa[nitrogen] == G::TPKA_HYDROXAMIC
                && graph->tautomerWeight[oxygen] == G::TW_HYDROXAMIC && graph->tautomerPKa[oxygen] == G::TPKA_HYDROXAMIC,
                "The hydroxamic acid motif must come before the amide motif in either atom order");
    const auto pyridinone = baseline("[nH]1ccccc1=O"), pyridinoneReversed = baseline("O=c1cccc[nH]1");
    for (const auto& [graph, nitrogen, oxygen] : {std::make_tuple(&pyridinone, 0, 6), std::make_tuple(&pyridinoneReversed, 6, 0)})
        for (int atom : {nitrogen, oxygen})
            require(graph->tautomerWeight[atom] == G::TW_PYRIDINONE && graph->tautomerSite[atom] == G::TAUTOMER_SITE_ACID
                    && graph->tautomerPKa[atom] == G::TPKA_PYRIDINONE,
                    "The pyridinone motif must come before the amide motif in either atom order");
    // The phenol motif reaches the carbonyl and para carbons of 2-pyridone, but its oxygen holds
    // no hydrogen, so the class keeps the N-H pKa of the pyridinone motif alone.
    for (int atom : {2, 5})
        require(pyridinone.tautomerClass[atom] == pyridinone.tautomerClass[0]
                && pyridinone.tautomerWeight[atom] == G::TW_PHENOL_QUINONE
                && pyridinone.tautomerSite[atom] == G::TAUTOMER_SITE_NONE,
                "A lactam carbonyl must not give the phenol motif an acid site");
    const auto urea = baseline("NC(N)=O");
    for (int atom : {0, 2, 3})
        require(urea.tautomerWeight[atom] == G::TW_UREA && urea.tautomerPKa[atom] == G::TPKA_UREA,
                "The urea motif must come before the amide motif");
    const auto succinimide = baseline("N1C(=O)CCC1=O");
    require(succinimide.tautomerWeight[0] == G::TW_IMIDE && succinimide.tautomerPKa[0] == G::TPKA_IMIDE
            && succinimide.tautomerWeight[1] == G::TW_IMIDE && succinimide.tautomerWeight[2] == G::TW_BARBITURIC
            && succinimide.tautomerSite[2] == G::TAUTOMER_SITE_NONE,
            "A cyclic imide must keep its imide N-H site and give its oxygens the barbituric weight");
    const auto acetylacetone = baseline("CC(=O)CC(C)=O");
    for (int atom : {2, 3, 6})
        require(acetylacetone.tautomerWeight[atom] == G::TW_DIKETONE_ENOL && acetylacetone.tautomerPKa[atom] == G::TPKA_DIKETONE_ENOL,
                "The 1,3-diketone motif must come before the keto/enol motif");
    require(acetylacetone.tautomerWeight[0] == G::TW_KETO_ENOL && acetylacetone.tautomerPKa[0] == G::TPKA_KETO_ENOL,
            "A methyl carbon outside the 1,3-diketone keeps the keto/enol motif");
    // 2-Hydroxypyridine: the phenol motif reaches the para carbon through the O-H of the
    // pyridone/hydroxypyridine motif, so that carbon carries the pKa of that motif, not the phenol pKa.
    const auto hydroxypyridine = baseline("Oc1ccccn1");
    require(hydroxypyridine.tautomerClass[0] >= 0 && hydroxypyridine.tautomerClass[4] == hydroxypyridine.tautomerClass[0]
            && hydroxypyridine.tautomerWeight[4] == G::TW_PHENOL_QUINONE
            && hydroxypyridine.tautomerSite[0] == G::TAUTOMER_SITE_ACID && hydroxypyridine.tautomerSite[4] == G::TAUTOMER_SITE_ACID
            && hydroxypyridine.tautomerPKa[0] == G::TPKA_PYRIDONE_HYDROXYPYRIDINE
            && hydroxypyridine.tautomerPKa[4] == G::TPKA_PYRIDONE_HYDROXYPYRIDINE,
            "2-Hydroxypyridine must keep the pKa of its own motif at its para carbon");
    for (double pH : {2.0, 7.4, 10.0, 12.0}) {
        const auto fractions = smsd::detail::tautomerNeutralFractions(hydroxypyridine, pH);
        const float acid = static_cast<float>(1.0 / (1.0 + std::pow(10.0, pH - G::TPKA_PYRIDONE_HYDROXYPYRIDINE)));
        for (int atom : {0, 1, 4, 6})
            require(fractions[atom] == acid, "Every site of a class must take the neutral fraction of the class");
        require(std::isnan(fractions[2]), "An atom outside every class must keep its weight at every pH");
    }
    // 4(3H)-Pyrimidinone takes its own N-H pKa in either tautomer, and the hypoxanthine N1-H too.
    for (const char* smiles : {"O=c1ccnc[nH]1", "Oc1ccncn1"}) {
        const auto pyrimidinone = baseline(smiles);
        require(pyrimidinone.tautomerSite[0] == G::TAUTOMER_SITE_ACID && pyrimidinone.tautomerPKa[0] == 8.6f
                && pyrimidinone.tautomerSite[6] == G::TAUTOMER_SITE_ACID && pyrimidinone.tautomerPKa[6] == 8.6f,
                std::string("A pyrimidin-4(3H)-one match must take pKa 8.6: ") + smiles);
    }
    require(baseline("O=c1[nH]cnc2[nH]cnc12").tautomerPKa[2] == 8.6f, "Hypoxanthine N1-H must take pKa 8.6");
    // A site needs the hydrogen it would release: no phenol site on an ether, an ester or a lactone
    // oxygen, no hydroxamic site on an N-alkoxy amide and no amide site on a tertiary amide.
    for (const char* smiles : {"COc1ccccc1", "O=c1ccc2ccccc2o1", "CON(C)C(C)=O", "CONC(C)=O", "CN(C)C=O",
                               "CN1C(=O)CCC1=O", "Cn1ccccc1=O"}) {
        const auto graph = baseline(smiles);
        for (int atom = 0; atom < graph.n; ++atom)
            require(graph.tautomerSite[atom] == G::TAUTOMER_SITE_NONE,
                    std::string("A match without the hydrogen of its site must carry no site: ") + smiles);
    }
    const auto arylEster = baseline("CC(=O)Oc1ccccc1");
    require(arylEster.tautomerClass[3] >= 0 && arylEster.tautomerSite[3] == G::TAUTOMER_SITE_NONE
            && arylEster.tautomerSite[4] == G::TAUTOMER_SITE_NONE && arylEster.tautomerSite[0] == G::TAUTOMER_SITE_ACID,
            "An aryl ester oxygen must carry no phenol site, while its acetyl keeps the keto/enol site");
    const auto phenol = baseline("Oc1ccccc1");
    require(phenol.tautomerSite[0] == G::TAUTOMER_SITE_ACID && phenol.tautomerPKa[0] == G::TPKA_PHENOL_QUINONE,
            "A phenol O-H must carry the phenol site");
    const auto weinreb = baseline("CON(C)C(C)=O");
    require(weinreb.tautomerClass[2] >= 0 && weinreb.tautomerWeight[2] == G::TW_HYDROXAMIC,
            "An N-alkoxy amide keeps the hydroxamic class and weight without its site");
    const auto acetone = baseline("CC(=O)C");
    const auto acetoneFractions = smsd::detail::tautomerNeutralFractions(acetone, 12.0);
    require(acetoneFractions[2] == static_cast<float>(1.0 / (1.0 + std::pow(10.0, 12.0 - G::TPKA_KETO_ENOL))),
            "A class with one acid site must keep the one-site formula");
    const auto amidine = baseline("CC(=N)N");
    for (float fraction : smsd::detail::tautomerNeutralFractions(amidine, 12.0))
        require(std::isnan(fraction), "A motif to which the library assigns no pKa must not be changed by pH");
}

/**
 * Random atom orders, with a fixed seed, must not change any class, weight or site: the classes
 * partition the same atoms in every order, and their numbers agree or differ only by a symmetry
 * of the graph.
 */
#ifdef SMSD_TAUTOMER_PERMUTATION_SMILES
static void atomOrderPermutation() {
    std::ifstream corpus(SMSD_TAUTOMER_PERMUTATION_SMILES);
    require(corpus.good(), "The permutation corpus must be readable");
    std::mt19937 random(20261021u);
    int molecules = 0, symmetric = 0;
    std::string line;
    while (molecules < 600 && std::getline(corpus, line)) {
        std::istringstream fields(line);
        std::string smiles;
        fields >> smiles;
        if (smiles.empty()) continue;
        smsd::MolGraph graph;
        try { graph = smsd::parseSMILES(smiles); } catch (const std::exception&) { continue; }
        ++molecules;
        const auto motifs = smsd::detail::tautomerMotifs(graph);
        std::vector<int> identity(graph.n);
        std::iota(identity.begin(), identity.end(), 0);
        const auto classes = classRepresentatives(graph, identity);
        const auto fractions = smsd::detail::tautomerNeutralFractions(graph, 12.0);
        for (int trial = 0; trial < 3; ++trial) {
            std::vector<int> permutation;
            const auto copy = reordered(graph, random, permutation);
            const auto copyMotifs = smsd::detail::tautomerMotifs(copy);
            smsd::MolGraph baseline = copy;
            baseline.computeTautomerClasses();
            for (int atom = 0; atom < copy.n; ++atom) {
                const int original = permutation[atom];
                require(copyMotifs.weight[atom] == motifs.weight[original] && copyMotifs.site[atom] == motifs.site[original]
                        && copyMotifs.pKa[atom] == motifs.pKa[original], "Atom order changed a motif weight or site: " + smiles);
                const bool inClass = baseline.tautomerClass[atom] >= 0;
                require(inClass == (copyMotifs.covered[atom] != 0), "Every atom a motif covers must be in a class: " + smiles);
                require(baseline.tautomerWeight[atom] == copyMotifs.weight[atom] && baseline.tautomerSite[atom] == copyMotifs.site[atom]
                        && baseline.tautomerPKa[atom] == copyMotifs.pKa[atom],
                        "An atom in a class must take its motif weight and site: " + smiles);
            }
            require(classRepresentatives(copy, permutation) == classes, "Atom order changed a tautomer class: " + smiles);
            std::vector<int> numbers(graph.n);
            for (int atom = 0; atom < copy.n; ++atom) numbers[permutation[atom]] = copy.tautomerClass[atom];
            if (numbers != graph.tautomerClass) {
                require(symmetricNumbering(graph, graph.tautomerClass, numbers),
                        "Atom order changed the class numbers beyond a symmetry of the graph: " + smiles);
                ++symmetric;
            }
            const auto copyFractions = smsd::detail::tautomerNeutralFractions(copy, 12.0);
            for (int atom = 0; atom < copy.n; ++atom) {
                const int original = permutation[atom];
                require(copy.tautomerWeight[atom] == graph.tautomerWeight[original]
                        && (copyFractions[atom] == fractions[original]
                            || (std::isnan(copyFractions[atom]) && std::isnan(fractions[original]))),
                        "Atom order changed a weight after the pH: " + smiles);
            }
        }
    }
    require(molecules >= 500, "The permutation test needs at least 500 molecules");
    // Classes that a symmetry exchanges take their numbers by atom order, so some orders differ.
    require(symmetric > 0, "The corpus must exercise classes that a symmetry exchanges");
    std::cout << "PASS atom order: " << molecules << " molecules in three orders each, the same classes in every order; "
              << symmetric << " orders numbered the classes differently by a symmetry\n";
}
#else
static void atomOrderPermutation() {
    std::cout << "SKIP atom order: the permutation corpus benchmarks/data/chembl_mcs_benchmark.smi is not in this "
                 "source tree\n";
}
#endif

/** The tautomer classes, by atom, of molecules each written in two atom orders; the classes must not depend on the order. */
static void classesOfReportedMolecules() {
    using Case = std::pair<const char*, std::vector<int>>;
    for (const auto& item : {Case{"ONC(C)=O", {0, 0, 0, 0, 0}}, Case{"CC(=O)NO", {0, 0, 0, 0, 0}},
                             Case{"O=C1CCC(=O)N1", {0, 0, 0, 0, 0, 0, 0}}, Case{"N1C(=O)CCC1=O", {0, 0, 0, 0, 0, 0, 0}},
                             Case{"Oc1ccccn1", {0, 0, -1, -1, 0, -1, 0}}, Case{"n1ccccc1O", {0, -1, 0, -1, -1, 0, 0}}})
        require(smsd::parseSMILES(item.first).tautomerClass == item.second,
                std::string("Tautomer classes of ") + item.first);
    // Two carbonyl classes that the mirror symmetry of 2,5-hexanedione exchanges: the one with the
    // lower atom index takes number 0.
    require(smsd::parseSMILES("CC(=O)CCC(C)=O").tautomerClass == std::vector<int>({0, -1, 0, 0, 1, -1, 1, 1}),
            "Classes exchanged by a symmetry are numbered from the lowest atom");
}

int main() {
    try {
        const auto query = smsd::parseSMILES("CC(=O)C"), target = smsd::parseSMILES("CC(O)=C");
        const auto match = identitySelection(query);
        const auto queryWeights = query.tautomerWeight, targetWeights = target.tautomerWeight;
        const auto aqueous = smsd::ChemOptions::tautomerProfile();
        auto chloroform = aqueous; chloroform.withSolvent(smsd::ChemOptions::Solvent::CHLOROFORM);
        auto acidic = aqueous; acidic.withPH(2);
        const double baseline = smsd::MCS::computeTautomerConfidence(query, target, match, aqueous);
        const double solvent = smsd::MCS::computeTautomerConfidence(query, target, match, chloroform);
        const double lowPH = smsd::MCS::computeTautomerConfidence(query, target, match, acidic);
        require(baseline > 0 && baseline < 1, "Selected keto/enol differences must affect confidence");
        require(solvent < baseline, "Solvent must affect the explicit confidence");
        require(lowPH == baseline, "A keto/enol site of pKa 20 is not ionised at pH 2, so its confidence must not move");
        require(std::abs(baseline - 0.98) < 1e-6 && std::abs(solvent - 0.73) < 1e-6,
                "Context formula differs from its baseline solvent contract");
        // Henderson-Hasselbalch: an acid site keeps its neutral form HA below its pKa, a base site
        // its neutral form B above its pKa, and the weight follows that neutral fraction.
        auto atPH = [](const char* first, const char* second, double pH) {
            const auto from = smsd::parseSMILES(first), to = smsd::parseSMILES(second);
            auto context = smsd::ChemOptions::tautomerProfile(); context.withPH(pH);
            return smsd::MCS::computeTautomerConfidence(from, to, identitySelection(from), context);
        };
        const double acid2 = atPH("ONC(C)=O", "ON=C(C)O", 2), acid74 = atPH("ONC(C)=O", "ON=C(C)O", 7.4),
                     acid12 = atPH("ONC(C)=O", "ON=C(C)O", 12);
        require(acid2 > acid74 && acid74 > acid12, "An acid site must lose its neutral form as the pH rises");
        const double hydroxamicPKa = smsd::MolGraph::TPKA_HYDROXAMIC;
        require(std::abs(acid2 - 0.85 / (1 + std::pow(10.0, 2 - hydroxamicPKa))) < 1e-6
                && std::abs(acid74 - 0.85 / (1 + std::pow(10.0, 7.4 - hydroxamicPKa))) < 1e-6
                && std::abs(acid12 - 0.1) < 1e-6,
                "A hydroxamic acid site must follow 1 / (1 + 10^(pH - pKa)) down to the 0.1 floor");
        // 4(3H)-Pyrimidinone holds the N-H site of the pyridinone motif, of pKa 8.6, alone: no
        // aromatic carbon path reaches a para carbon from its carbonyl carbon.
        const double pyrimidinone2 = atPH("O=c1ccnc[nH]1", "Oc1ccncn1", 2),
                     pyrimidinone74 = atPH("O=c1ccnc[nH]1", "Oc1ccncn1", 7.4),
                     pyrimidinone12 = atPH("O=c1ccnc[nH]1", "Oc1ccncn1", 12);
        require(pyrimidinone2 > pyrimidinone74 && pyrimidinone74 > pyrimidinone12,
                "An acid site must lose its neutral form as the pH rises");
        require(std::abs(pyrimidinone74 - 0.95 / (1 + std::pow(10.0, 7.4 - 8.6))) < 1e-6
                && std::abs(pyrimidinone12 - 0.1) < 1e-6,
                "A pyrimidin-4(3H)-one site must follow 1 / (1 + 10^(pH - 8.6))");
        // 2-Pyridone keeps the N-H pKa of 11.7 of its own motif: the phenol motif at its carbonyl
        // and para carbons has no O-H, so the weight holds up to pH 10 and falls above it.
        const double pyridone10 = atPH("[nH]1ccccc1=O", "n1ccccc1O", 10), pyridone12 = atPH("[nH]1ccccc1=O", "n1ccccc1O", 12);
        require(pyridone10 > 0.9 && pyridone12 < pyridone10, "2-Pyridone must lose its neutral form above its N-H pKa");
        for (int repeat = 0; repeat < 4; ++repeat) {
            require(smsd::MCS::computeTautomerConfidence(query, target, match, chloroform) == solvent,
                    "Repeated solvent confidence must not compound");
            require(smsd::MCS::computeTautomerConfidence(query, target, match, aqueous) == baseline,
                    "Context order must not alter the aqueous confidence");
        }
        for (auto context : {aqueous, chloroform, acidic})
            require(smsd::MCS::computeTautomerConfidence(query, query, match, context) == 1,
                    "Identical selected chemistry must have confidence one");
        auto unchanged = match;
        unchanged.bondMapping.erase(std::remove_if(unchanged.bondMapping.begin(), unchanged.bondMapping.end(),
            [&](const auto& bond) { return query.bondOrder(bond[0], bond[1]) != target.bondOrder(bond[2], bond[3]); }),
            unchanged.bondMapping.end());
        require(smsd::MCS::computeTautomerConfidence(query, target, unchanged, aqueous) == 1,
                "Unselected bond differences must not penalise confidence");
        std::vector<std::future<double>> workers;
        for (int worker = 0; worker < 8; ++worker)
            workers.push_back(std::async(std::launch::async, [&] {
                return smsd::MCS::computeTautomerConfidence(query, target, match, chloroform);
            }));
        for (auto& worker : workers) require(worker.get() == solvent, "Concurrent context reads must agree");
        require(query.tautomerWeight == queryWeights && target.tautomerWeight == targetWeights,
                "Context confidence must not mutate molecule weights");
        auto contextualQuery = query;
        contextualQuery.applySolventCorrection(smsd::ChemOptions::Solvent::CHLOROFORM);
        require(smsd::MCS::computeTautomerConfidence(contextualQuery, target, match, chloroform) == solvent,
                "Existing graph context must not be applied twice");
        std::unordered_map<int, int> legacy(match.atomMapping.begin(), match.atomMapping.end());
        require(legacyConfidence(query, target, legacy) == baseline,
                "The deprecated three-argument form must equal the four-argument form in water at pH 7.4");
        for (double value : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()}) {
            auto invalid = aqueous; invalid.pH = value;
            bool rejected = false;
            try { smsd::MCS::computeTautomerConfidence(query, target, match, invalid); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Non-finite confidence pH must be rejected");
            rejected = false;
            try { invalid.withPH(value); } catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Non-finite fluent pH must be rejected");
        }
        auto malformed = match; malformed.atomMapping[query.n] = 0;
        bool rejected = false;
        try { smsd::MCS::computeTautomerConfidence(query, target, malformed, aqueous); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "Out-of-range selected atoms must be rejected");
        legacy[query.n] = 0; rejected = false;
        try { legacyConfidence(query, target, legacy); }
        catch (const std::out_of_range&) { rejected = true; }
        require(rejected, "The three-argument form must reject out-of-range atoms");
        auto acetone = [](double pH) {
            return smsd::MolGraph::Builder().atomCount(4).atomicNumbers({6, 6, 8, 6})
                .setNeighbors({{1}, {0, 2, 3}, {1}, {1}}).setBondOrders({{1}, {1, 2, 1}, {2}, {1}})
                .pH(pH).build();
        };
        const auto neutralBuilt = acetone(7.4), acidicBuilt = acetone(2.0), clampedBuilt = acetone(-3.0);
        require(neutralBuilt.pH_ == 7.4 && acidicBuilt.pH_ == 2.0 && clampedBuilt.pH_ == 0.0,
                "Builder pH must reach the built graph, clamped to [0, 14]");
        require(neutralBuilt.tautomerWeight[2] == query.tautomerWeight[2]
                && acidicBuilt.tautomerWeight[2] == neutralBuilt.tautomerWeight[2],
                "A keto/enol site of pKa 20 must keep its builder weight between pH 2 and 7.4");
        auto hydroxamic = [](double pH) {
            return smsd::MolGraph::Builder().atomCount(5).atomicNumbers({8, 7, 6, 6, 8})
                .setNeighbors({{1}, {0, 2}, {1, 3, 4}, {2}, {2}}).setBondOrders({{1}, {1, 1}, {1, 1, 2}, {1}, {2}})
                .pH(pH).build();
        };
        const auto hydroxamicNeutral = hydroxamic(2.0), hydroxamicAlkaline = hydroxamic(12.0);
        require(hydroxamicNeutral.tautomerSite[1] == smsd::MolGraph::TAUTOMER_SITE_ACID
                && hydroxamicNeutral.tautomerPKa[1] == smsd::MolGraph::TPKA_HYDROXAMIC
                && std::abs(hydroxamicNeutral.tautomerWeight[1] - 0.85f) < 1e-5f
                && hydroxamicAlkaline.tautomerWeight[1] == 0.1f,
                "Builder pH must set the Henderson-Hasselbalch tautomer weights");
        rejected = false;
        try { smsd::MolGraph::Builder().pH(std::numeric_limits<double>::quiet_NaN()); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Non-finite builder pH must be rejected");
        // Proton conservation counts the stored hydrogens, so an aromatic N-H counts one, as in Java.
        const auto hydroxypyridine = smsd::parseSMILES("Oc1ccccn1");
        for (const char* tautomer : {"O=c1cccc[nH]1", "O=C1C=CC=CN1"}) {
            const auto identity = identitySelection(hydroxypyridine).atomMapping;
            require(smsd::validateTautomerConsistency(hydroxypyridine, smsd::parseSMILES(tautomer), identity),
                    std::string("Tautomer pair must conserve its mobile proton: ") + tautomer);
        }
        require(smsd::mobilisableProtons(smsd::parseSMILES("O=c1cccc[nH]1"), 6) == 1
                && smsd::mobilisableProtons(hydroxypyridine, 0) == 1 && smsd::mobilisableProtons(hydroxypyridine, 6) == 0,
                "Mobilisable protons must read the stored hydrogen counts");
        require(!smsd::validateTautomerConsistency(query, target, identitySelection(query).atomMapping),
                "A ketone mapped onto its enol moves a proton out of its class");
        motifPriorityAndClassModel();
        atomOrderPermutation();
        classesOfReportedMolecules();
        std::cout << "PASS selected chemistry, context, repeatability, concurrency and bounds\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}

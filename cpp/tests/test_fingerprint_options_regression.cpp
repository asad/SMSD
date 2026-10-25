/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/fingerprint_options.hpp"
#include "smsd/smiles_parser.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using smsd_test::require;
using namespace smsd::batch;

namespace {

const std::vector<std::string>& molecules() {
    static const std::vector<std::string> list = {
        "C", "CC", "CCO", "CC(=O)O", "c1ccccc1", "c1ccncc1", "Oc1ccccc1", "Cc1ccccc1", "CC(=O)Oc1ccccc1C(=O)O",
        "CC(=O)Nc1ccc(O)cc1", "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "C[C@H](N)C(=O)O",
        "C[C@@H](N)C(=O)O", "C/C=C/C(=O)O", "C/C=C\\C(=O)O", "[13CH4]", "[NH4+]", "[O-]c1ccccc1", "[Na+].CC(=O)[O-]",
        "CCO.CCO", "c1ccccc1.CCO", "C1CC1", "C1CCCCCCC1", "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3",
        "c1ccc2ccccc2c1", "C1CCC2(C1)CCCCC2", "O=C1CCCCC1", "c1ccccc1CCc1ccccc1", "NS(=O)(=O)c1ccc(Cl)cc1",
        "C#CC1(O)CCC2C3CCC4=CC(=O)CCC4(C)C3CCC21C", "OC(=O)CCC(=O)O", "CC(=O)CC(C)=O", "[H]OC", "[H]N([H])[H]",
        "BrC(Br)Br", "FC(F)(F)c1ccccc1", "CC#N", "C=N", "", "[He].[Ne]", "C1CC1.CC.[He]",
    };
    return list;
}

smsd::MolGraph parse(const std::string& smiles) {
    return smiles.empty() ? smsd::MolGraph() : smsd::parseSMILES(smiles);
}

void classicConfigurationsEqualTheFixedFunctions() {
    for (const auto& smiles : molecules()) {
        const auto graph = parse(smiles);
        for (int radius : {0, 1, 2, 3, -1})
            for (int size : {64, 1000, 2048}) {
                require(computeFingerprint(graph, FingerprintOptions::ecfp(radius, size))
                        == detail::computeCircularFingerprintECFP(graph, radius, size), "Classic ECFP bits must not change");
                require(computeFingerprintCounts(graph, FingerprintOptions::ecfp(radius, size))
                        == detail::computeCircularFingerprintECFPCounts(graph, radius, size), "Classic ECFP counts must not change");
                require(computeFingerprint(graph, FingerprintOptions::fcfp(radius, size))
                        == detail::computeCircularFingerprintFCFP(graph, radius, size), "Classic FCFP bits must not change");
                require(computeFingerprintCounts(graph, FingerprintOptions::fcfp(radius, size))
                        == detail::computeCircularFingerprintFCFPCounts(graph, radius, size), "Classic FCFP counts must not change");
            }
        for (int length : {1, 3, 7})
            for (int size : {64, 2048})
                require(computeFingerprint(graph, FingerprintOptions::path(length, size))
                        == detail::computeElementPathFingerprint(graph, length, size), "Classic path bits must not change");
        for (int size : {64, 1000, 2048}) {
            require(computeFingerprint(graph, FingerprintOptions::torsion(size))
                    == detail::computeTopologicalTorsion(graph, size), "Classic torsion bits must not change");
            require(computeFingerprintCounts(graph, FingerprintOptions::torsion(size))
                    == detail::computeTopologicalTorsionCounts(graph, size), "Classic torsion counts must not change");
        }
    }
    require(computeFingerprint(parse("CCO"), FingerprintOptions()) == detail::computeCircularFingerprintECFP(parse("CCO"), 2, 2048),
            "The default options are the classic ECFP4 configuration");
}

void overloadsAgree() {
    std::vector<smsd::MolGraph> graphs;
    for (const auto& smiles : molecules()) graphs.push_back(parse(smiles));
    auto options = FingerprintOptions::ecfp(2, 512);
    options.cipAtoms = options.ringSize = true;
    const auto batch = batchFingerprint(graphs, options, 1);
    require(batch.size() == graphs.size(), "batchFingerprint returns one fingerprint per molecule");
    for (size_t i = 0; i < graphs.size(); ++i) {
        require(batch[i] == computeFingerprint(graphs[i], options), "batchFingerprint agrees with computeFingerprint");
        if (!molecules()[i].empty())
            require(fingerprintFromSmiles(molecules()[i], options) == batch[i], "fingerprintFromSmiles agrees with computeFingerprint");
    }
}

/// SMILES parsed as they are read, as the Python batch reads its strings.
struct ParsedOnRead {
    std::vector<std::string> smiles;
    size_t size() const { return smiles.size(); }
    smsd::MolGraph operator[](size_t index) const { return parse(smiles[index]); }
};

/// The message of the exception a call throws, or the empty string.
template<class Call>
std::string failureOf(const Call& call) {
    try { call(); } catch (const std::exception& error) { return error.what(); }
    return "";
}

void batchesAgreeOnEveryThreadCount() {
    ParsedOnRead library;
    for (int copy = 0; copy < 10; ++copy)
        library.smiles.insert(library.smiles.end(), molecules().begin(), molecules().end());
    std::vector<smsd::MolGraph> graphs;
    for (const auto& smiles : library.smiles) graphs.push_back(parse(smiles));
    for (const auto& options : {FingerprintOptions::ecfp(3, 1024), FingerprintOptions::path(7, 2048), FingerprintOptions::maccs()}) {
        const auto expected = batchFingerprint(graphs, options, 1);
        for (size_t i = 0; i < graphs.size(); ++i)
            require(expected[i] == computeFingerprint(graphs[i], options), "One thread agrees with computeFingerprint");
        for (int threads : {2, 4, 0}) {
            require(batchFingerprint(graphs, options, threads) == expected, "The thread count does not change the bits");
            require(batchFingerprint(library, options, threads) == expected, "A sequence parsed on reading gives the same bits");
        }
    }
    // The first molecule that fails decides the error, however the molecules are shared out.
    ParsedOnRead failing{std::vector<std::string>(600, "CCO")};
    failing.smiles[301] = "C1CC";
    for (size_t i = 302; i < failing.smiles.size(); i += 2) failing.smiles[i] = "C(C";
    const std::string first = failureOf([] { parse("C1CC"); });
    require(!first.empty() && first != failureOf([] { parse("C(C"); }), "The two failures are told apart");
    for (int threads : {1, 2, 4, 0})
        require(failureOf([&] { batchFingerprint(failing, FingerprintOptions(), threads); }) == first,
                "The batch throws the exception of the first molecule that fails");
}

void stereoTermsSeparateIsomers() {
    const auto l = parse("C[C@H](N)C(=O)O"), d = parse("C[C@@H](N)C(=O)O");
    for (FingerprintFamily family : {FingerprintFamily::ECFP, FingerprintFamily::FCFP, FingerprintFamily::PATH, FingerprintFamily::TORSION}) {
        auto plain = FingerprintOptions::of(family), chiral = plain;
        chiral.cipAtoms = true;
        require(computeFingerprint(l, plain) == computeFingerprint(d, plain), "Enantiomers share their bits without the CIP term");
        require(computeFingerprint(l, chiral) != computeFingerprint(d, chiral), "Enantiomers differ with the CIP atom term");
    }
    auto radiusZero = FingerprintOptions::ecfp(0, 2048);
    radiusZero.cipAtoms = true;
    require(computeFingerprint(l, radiusZero) == computeFingerprint(d, radiusZero),
            "The descriptor joins after the first iteration, so radius 0 does not separate enantiomers");
    auto chiral = FingerprintOptions::ecfp(2, 4096);
    chiral.cipAtoms = true;
    const auto lCounts = computeFingerprintCounts(l, chiral), dCounts = computeFingerprintCounts(d, chiral);
    int differing = 0;
    for (size_t i = 0; i < lCounts.size(); ++i) differing += lCounts[i] != dCounts[i];
    // The centre at radius 1 and 2, and its three neighbours at radius 2: five hashes each way.
    require(differing <= 10, "Enantiomers differ only where the centre lies");

    const auto e = parse("C/C=C/C(=O)O"), z = parse("C/C=C\\C(=O)O");
    for (FingerprintFamily family : {FingerprintFamily::ECFP, FingerprintFamily::PATH, FingerprintFamily::TORSION}) {
        auto plain = FingerprintOptions::of(family), stereo = plain;
        stereo.cipBonds = true;
        require(computeFingerprint(e, plain) == computeFingerprint(z, plain), "E/Z isomers share their bits without the CIP term");
        require(computeFingerprint(e, stereo) != computeFingerprint(z, stereo), "E/Z isomers differ with the CIP bond term");
    }
}

void ringSizeClasses() {
    auto options = FingerprintOptions::ecfp(2, 2048);
    options.ringSize = true;
    const auto seven = computeFingerprint(parse("C1CCCCCC1"), options), eight = computeFingerprint(parse("C1CCCCCCC1"), options);
    const auto nine = computeFingerprint(parse("C1CCCCCCCC1"), options);
    require(seven != eight, "Ring sizes below the cap are distinct classes");
    require(eight == nine, "Rings of the cap size or more share one class");
    options.ringSizeCap = 12;
    require(computeFingerprint(parse("C1CCCCCCC1"), options) != computeFingerprint(parse("C1CCCCCCCC1"), options),
            "A higher cap separates the larger rings");
    auto plain = FingerprintOptions::ecfp(2, 2048);
    require(computeFingerprint(parse("C1CCCCCC1"), plain) == computeFingerprint(parse("C1CCCCCCC1"), plain),
            "Without the ring size term, saturated rings of different sizes share their bits");
}

void fragmentHandling() {
    auto largest = FingerprintOptions::ecfp(2, 2048);
    largest.fragments = FragmentMode::LARGEST;
    require(computeFingerprint(parse("CCO.CC(=O)O"), largest) == computeFingerprint(parse("CC(=O)O"), FingerprintOptions::ecfp()),
            "The largest fragment is the component with the most heavy atoms");
    require(computeFingerprint(parse("[Na+].CC(=O)[O-]"), largest) == computeFingerprint(parse("CC(=O)[O-]"), FingerprintOptions::ecfp()),
            "A counter-ion is dropped by the largest-fragment mode");
    auto perFragment = FingerprintOptions::ecfp(2, 2048);
    perFragment.fragments = FragmentMode::PER_FRAGMENT;
    require(computeFingerprint(parse("c1ccccc1.CCO"), perFragment) != computeFingerprint(parse("c1ccccc1.CCO"), FingerprintOptions::ecfp()),
            "Per-fragment hashing makes bits component-specific");
    require(computeFingerprint(parse("CCO"), perFragment) != computeFingerprint(parse("CCO"), FingerprintOptions::ecfp()),
            "Per-fragment hashing folds the component class even for one component");
    require(computeFingerprint(parse("CCO.CCO"), perFragment) == computeFingerprint(parse("CCO"), perFragment),
            "Identical components share one class");
}

void scaffoldHandling() {
    auto scaffold = FingerprintOptions::ecfp(2, 2048);
    scaffold.scaffold = ScaffoldMode::MURCKO;
    const auto benzene = computeFingerprint(parse("c1ccccc1"), FingerprintOptions::ecfp());
    require(computeFingerprint(parse("Cc1ccccc1"), scaffold) == benzene, "The framework of toluene is benzene, capped with hydrogen");
    require(computeFingerprint(parse("CC(=O)Nc1ccc(O)cc1"), scaffold) == benzene, "Side chains are removed down to the ring");
    require(computeFingerprint(parse("O=C1CCCCC1"), scaffold) == computeFingerprint(parse("O=C1CCCCC1"), FingerprintOptions::ecfp()),
            "An atom double bonded to a ring atom stays in the scaffold");
    require(computeFingerprint(parse("Cn1ccccc1=O"), scaffold) == computeFingerprint(parse("O=c1cccc[nH]1"), FingerprintOptions::ecfp()),
            "The scaffold of an N-methylpyridone is pyridone, as murckoScaffold gives it");
    require(computeFingerprint(parse("CC(=O)c1ccccc1"), scaffold) == benzene, "A side chain is cut at its single bond");
    require(computeFingerprint(parse("CN=C1CCCCC1"), scaffold) == computeFingerprint(parse("N=C1CCCCC1"), FingerprintOptions::ecfp()),
            "A kept exocyclic atom is capped where its own side chain was");
    require(computeFingerprint(parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C"), scaffold)
            == computeFingerprint(parse("[nH]1cnc2c1c(=O)[nH]c(=O)[nH]2"), FingerprintOptions::ecfp()),
            "The framework of caffeine is xanthine, with both carbonyls");
    {
        auto maccs = FingerprintOptions::maccs();
        maccs.scaffold = ScaffoldMode::MURCKO;
        const auto keys = smsd::bitPositions(computeFingerprint(parse("Cn1cnc2c1c(=O)n(C)c(=O)n2C"), maccs), 166);
        require(std::find(keys.begin(), keys.end(), 27) == keys.end() && std::find(keys.begin(), keys.end(), 85) == keys.end(),
                "The framework of an aromatic amide has no methylene keys");
    }
    require(computeFingerprint(parse("c1ccccc1CCc1ccccc1"), scaffold) == computeFingerprint(parse("c1ccccc1CCc1ccccc1"), FingerprintOptions::ecfp()),
            "Linkers between rings are kept");
    require(computeFingerprint(parse("CCO"), scaffold) == computeFingerprint(parse("CCO"), FingerprintOptions::ecfp()),
            "An acyclic molecule has no framework and is fingerprinted whole");
    scaffold.cipAtoms = true;
    require(computeFingerprint(parse("c1ccccc1[C@H](c1ccncc1)c1ccco1"), scaffold)
            != computeFingerprint(parse("c1ccccc1[C@@H](c1ccncc1)c1ccco1"), scaffold),
            "A linker centre that keeps its neighbours keeps its descriptor");
    require(computeFingerprint(parse("O[C@H](c1ccccc1)c1ccccc1"), scaffold) == computeFingerprint(parse("O[C@@H](c1ccccc1)c1ccccc1"), scaffold),
            "A centre that is no longer a centre after capping has no descriptor");
    auto chiral = FingerprintOptions::ecfp();
    chiral.cipAtoms = true;
    require(computeFingerprint(parse("C[C@]12CCCC[C@@H]1CCCO2"), scaffold) == computeFingerprint(parse("[C@H]12CCCC[C@@H]1CCCO2"), chiral)
            && computeFingerprint(parse("C[C@@]12CCCC[C@@H]1CCCO2"), scaffold) == computeFingerprint(parse("[C@@H]12CCCC[C@@H]1CCCO2"), chiral)
            && computeFingerprint(parse("C[C@]12CCCC[C@@H]1CCCO2"), scaffold) != computeFingerprint(parse("C[C@@]12CCCC[C@@H]1CCCO2"), scaffold),
            "A ring-fusion centre keeps its configuration when its angular methyl is capped");
    auto geometric = FingerprintOptions::ecfp(), framedGeometric = FingerprintOptions::ecfp();
    geometric.cipBonds = true;
    framedGeometric.cipBonds = true;
    framedGeometric.scaffold = ScaffoldMode::MURCKO;
    require(computeFingerprint(parse("C/C(c1ccccc1)=C/c1ccccc1"), framedGeometric) == computeFingerprint(parse("c1ccccc1/C=C\\c1ccccc1"), geometric)
            && computeFingerprint(parse("C/C(c1ccccc1)=C\\c1ccccc1"), framedGeometric) == computeFingerprint(parse("c1ccccc1/C=C/c1ccccc1"), geometric),
            "A double bond keeps its configuration when its reference substituent is capped");
    require(computeFingerprint(parse("C1C/C=C\\CC1"), geometric) == computeFingerprint(parse("C1CC=CCC1"), geometric)
            && computeFingerprint(parse("C1C/C=C/CC1"), geometric) == computeFingerprint(parse("C1CC=CCC1"), geometric),
            "A double bond in a ring of fewer than eight atoms carries no E/Z term");
    require(computeFingerprint(parse("C1CCC/C=C/CC1"), geometric) != computeFingerprint(parse("C1CCC/C=C\\CC1"), geometric),
            "Cyclooctene keeps its E/Z term");
}

void substructureKeys() {
    require(FingerprintOptions::maccsKeys().size() == 166, "There are 166 MACCS keys");
    const auto maccs = FingerprintOptions::maccs();
    require(maccs.family == FingerprintFamily::KEYS && maccs.size == 166, "The MACCS options are a key family of 166 bits");
    const auto benzene = smsd::bitPositions(computeFingerprint(parse("c1ccccc1"), maccs), 166);
    const auto has = [](const std::vector<int>& positions, int key) {
        for (int position : positions) if (position == key - 1) return true;
        return false;
    };
    require(has(benzene, 163) && has(benzene, 165) && has(benzene, 162), "Benzene sets the six-membered ring, ring and aromatic keys");
    require(!has(benzene, 164) && !has(benzene, 161), "Benzene sets neither the oxygen nor the nitrogen key");
    const auto ethanol = smsd::bitPositions(computeFingerprint(parse("CCO"), maccs), 166);
    require(has(ethanol, 164) && has(ethanol, 157) && has(ethanol, 139) && !has(ethanol, 159), "Ethanol has one oxygen, bonded to carbon as a hydroxyl");
    require(has(smsd::bitPositions(computeFingerprint(parse("CC(=O)O.O"), maccs), 166), 166), "Two fragments set key 166");
    require(has(smsd::bitPositions(computeFingerprint(parse("[13CH4]"), maccs), 166), 1), "An isotope sets key 1");
    require(has(smsd::bitPositions(computeFingerprint(parse("c1ccc2ccccc2c1"), maccs), 166), 125), "Two aromatic rings set key 125");
    require(computeFingerprint(smsd::MolGraph(), maccs) == std::vector<uint64_t>(3, 0ULL), "An empty molecule sets no key");

    const auto custom = FingerprintOptions::withKeys({FingerprintKey("c1ccccc1"), FingerprintKey("[#8]", 2), FingerprintKey("[#8]", 5)});
    const auto aspirin = parse("CC(=O)Oc1ccccc1C(=O)O");
    require(smsd::bitPositions(computeFingerprint(aspirin, custom), 3) == std::vector<int>{0, 1}, "A key is set when its minimum count is met");
    require(computeFingerprintCounts(aspirin, custom) == std::vector<int>{1, 4, 4}, "Key counts are distinct matched atom sets");
}

void orderInvariance() {
    struct Case { const char* first; const char* second; FingerprintOptions options; };
    FingerprintOptions chiral = FingerprintOptions::ecfp(), geometric = FingerprintOptions::ecfp(), perFragment = FingerprintOptions::ecfp();
    FingerprintOptions largest = FingerprintOptions::ecfp(), framework = FingerprintOptions::ecfp(), ringSize = FingerprintOptions::ecfp();
    chiral.cipAtoms = true;
    geometric.cipBonds = true;
    perFragment.fragments = FragmentMode::PER_FRAGMENT;
    largest.fragments = FragmentMode::LARGEST;
    framework.scaffold = ScaffoldMode::MURCKO;
    ringSize.ringSize = true;
    const Case cases[] = {
        {"C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", chiral},
        {"C/C=C/C(=O)O", "OC(=O)/C=C/C", geometric},
        {"c1ccccc1.CCO", "OCC.c1ccccc1", perFragment},
        {"CCO.CCN", "CCN.CCO", largest},
        {"Cc1ccccc1CC(C)c1ccccc1O", "Oc1ccccc1C(C)Cc1ccccc1C", framework},
        {"C1CC1c1ccccc1", "c1ccc(cc1)C1CC1", ringSize},
    };
    for (const auto& item : cases)
        require(computeFingerprint(parse(item.first), item.options) == computeFingerprint(parse(item.second), item.options),
                std::string("Bits must not depend on the atom order: ") + item.first);
    require(computeFingerprint(parse("C[C@@H](N)C(=O)O"), chiral) != computeFingerprint(parse("N[C@@H](C)C(=O)O"), chiral),
            "The other enantiomer still differs");
    require(computeFingerprint(parse("CCN.CCO"), largest) == computeFingerprint(parse("CCN"), largest)
            || computeFingerprint(parse("CCN.CCO"), largest) == computeFingerprint(parse("CCO"), largest),
            "A tie between unlike components of equal size picks one of them");
}

void validation() {
    const auto graph = parse("CCO");
    const auto rejects = [&](FingerprintOptions options, const char* message) {
        bool rejected = false;
        try { computeFingerprint(graph, options); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, message);
    };
    auto zero = FingerprintOptions::ecfp(2, 0);
    rejects(zero, "A zero size is rejected");
    auto radius = FingerprintOptions::ecfp(-2, 64);
    rejects(radius, "A radius below -1 is rejected");
    auto length = FingerprintOptions::path(0, 64);
    rejects(length, "A zero path length is rejected");
    auto cap = FingerprintOptions::ecfp(2, 64);
    cap.ringSizeCap = 2;
    rejects(cap, "A ring size cap below 3 is rejected");
    rejects(FingerprintOptions::of(FingerprintFamily::KEYS), "A key family without keys is rejected");
    rejects(FingerprintOptions::withKeys({FingerprintKey("C", 0)}), "A key minimum below 1 is rejected");
    bool rejected = false;
    try { computeFingerprint(graph, FingerprintOptions::withKeys({FingerprintKey("C(")})); } catch (const std::exception&) { rejected = true; }
    require(rejected, "An unparsable SMARTS key is rejected");
}

} // namespace

int main() {
    try {
        classicConfigurationsEqualTheFixedFunctions();
        overloadsAgree();
        batchesAgreeOnEveryThreadCount();
        stereoTermsSeparateIsomers();
        ringSizeClasses();
        fragmentHandling();
        scaffoldHandling();
        substructureKeys();
        orderInvariance();
        validation();
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "PASS fingerprint options regression\n";
    return 0;
}

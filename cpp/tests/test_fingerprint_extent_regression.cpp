/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/batch.hpp"
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

using smsd_test::require;

static void finiteExtent() {
    for (const auto& fixture : {std::pair<const char*, int>{"CCO", 2}, {"c1ccccc1", 3},
                                {"CC(C)(C)C", 2}, {"C1CC1.CC.[He]", 1}, {"[He].[Ne]", 0}}) {
        const auto graph = smsd::parseSMILES(fixture.first);
        for (int size : {1, 8, 17, 1000, 2048}) {
            require(smsd::batch::detail::computeCircularFingerprintECFP(graph, -1, size)
                    == smsd::batch::detail::computeCircularFingerprintECFP(graph, fixture.second, size),
                    "Whole-graph ECFP must stop at the largest finite component distance");
            require(smsd::batch::detail::computeCircularFingerprintFCFP(graph, -1, size)
                    == smsd::batch::detail::computeCircularFingerprintFCFP(graph, fixture.second, size),
                    "Whole-graph FCFP must stop at the largest finite component distance");
            require(smsd::batch::detail::computeCircularFingerprintECFPCounts(graph, -1, size)
                    == smsd::batch::detail::computeCircularFingerprintECFPCounts(graph, fixture.second, size),
                    "Whole-graph count ECFP must not add repeated full-graph rounds");
            require(smsd::batch::detail::computeCircularFingerprintFCFPCounts(graph, -1, size)
                    == smsd::batch::detail::computeCircularFingerprintFCFPCounts(graph, fixture.second, size),
                    "Whole-graph count FCFP must not add repeated full-graph rounds");
        }
    }
    const auto benzene = smsd::parseSMILES("c1ccccc1");
    const auto counts = smsd::batch::detail::computeCircularFingerprintECFPCounts(benzene, -1, 2048);
    require(std::accumulate(counts.begin(), counts.end(), 0) == 24,
            "Benzene has six radius-zero atoms and three environment expansion rounds");
    require(smsd::batch::detail::computeCircularFingerprintECFP(benzene, -1, 2048)
            != smsd::batch::detail::computeCircularFingerprintECFP(benzene, 4, 2048),
            "The redundant fourth benzene expansion must not be included");
}

static void aromaticTorsions() {
    const auto benzene = smsd::parseSMILES("c1ccccc1"), cyclohexane = smsd::parseSMILES("C1CCCCC1");
    require(smsd::batch::detail::computeTopologicalTorsion(benzene, 2048)
            != smsd::batch::detail::computeTopologicalTorsion(cyclohexane, 2048),
            "Aromatic and saturated rings must have distinct torsion types");
    require(smsd::batch::detail::computeTopologicalTorsionCounts(benzene, 2048)
            != smsd::batch::detail::computeTopologicalTorsionCounts(cyclohexane, 2048),
            "Count torsions must retain aromatic pi electrons");
    require(smsd::batch::detail::computeTopologicalTorsion(benzene, 2048)
            == smsd::batch::detail::computeTopologicalTorsion(smsd::parseSMILES("C1=CC=CC=C1"), 2048),
            "Torsions must preserve aromatic/Kekule input invariance");
}

static void validSizes() {
    const auto graph = smsd::parseSMILES("CCO");
    for (int size : {0, -1}) {
        bool rejected = false;
        try { smsd::batch::detail::computeCircularFingerprintECFP(graph, 2, size); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Fingerprint size must be positive");
    }
    for (int length : {0, -1, std::numeric_limits<int>::max()}) {
        bool rejected = false;
        try { smsd::batch::detail::computeMCSFingerprint(graph, length, 2048); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "MCS path length must be positive and representable");
    }
}

// fingerprintQuality rejects a size below 1, a position outside the size and, for packed
// words, a size beyond the words given; repeated positions and padding bits are not counted.
static void qualityArguments() {
    const std::vector<std::pair<std::vector<int>, int>> invalid = {{{}, 0}, {{}, -1}, {{-1}, 8}, {{8}, 8}};
    for (const auto& example : invalid) {
        bool rejected = false;
        try { smsd::batch::fingerprintQuality(example.first, example.second); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "fingerprint quality accepted an invalid position or size");
    }
    for (int size : {0, -1, 65}) {
        bool rejected = false;
        try { smsd::batch::fingerprintQuality(std::vector<uint64_t>(1, 0), size); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "fingerprint quality accepted a size beyond the words given");
    }
    require(smsd::batch::fingerprintQuality(std::vector<int>{7, 0, 7}, 8).setBits == 2,
            "repeated positions were counted twice");
    require(smsd::batch::fingerprintQuality(std::vector<uint64_t>{~uint64_t(0), ~uint64_t(0)}, 65).setBits == 65,
            "bits beyond the size were counted");
}

static std::vector<int> bitPositions(const std::vector<uint64_t>& fingerprint) {
    std::vector<int> result;
    for (size_t word = 0; word < fingerprint.size(); ++word)
        for (int bit = 0; bit < 64; ++bit)
            if ((fingerprint[word] & (uint64_t(1) << bit)) != 0)
                result.push_back(static_cast<int>(word * 64 + bit));
    return result;
}

static void fixedRadiusGoldens() {
    const auto graph = smsd::parseSMILES("CCO");
    require(bitPositions(smsd::batch::detail::computeCircularFingerprintECFP(graph, 0, 2048))
            == std::vector<int>{452, 1828, 1861}, "Radius-zero ECFP hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computeCircularFingerprintECFP(graph, 2, 2048))
            == std::vector<int>{273, 452, 501, 1237, 1670, 1702, 1746, 1828, 1861},
            "Fixed-radius ECFP hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computeCircularFingerprintFCFP(graph, 2, 2048))
            == std::vector<int>{357, 395, 814, 1172, 1716, 1838, 1933, 1942, 2027},
            "Fixed-radius FCFP hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computeTopologicalTorsion(smsd::parseSMILES("C1CCCCC1"), 2048))
            == std::vector<int>{1177}, "Saturated torsion hashes must match the golden bits");
}

static void pathGoldens() {
    const auto acid = smsd::parseSMILES("OC(=O)c1ccncc1");
    smsd::ChemOptions tautomeric;
    tautomeric.tautomerAware = true;
    require(bitPositions(smsd::batch::detail::computeMCSFingerprint(acid, 3, 2048))
            == std::vector<int>{51, 91, 152, 251, 287, 408, 564, 565, 593, 674, 749, 750, 805,
                                835, 837, 948, 956, 990, 1068, 1132, 1730, 1748, 1755, 1916, 1976},
            "MCS path hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computeMCSFingerprintForChemistry(acid, 3, 2048, tautomeric))
            == std::vector<int>{43, 327, 367, 445, 460, 786, 787, 819, 950, 1012, 1027, 1029, 1111,
                                1131, 1224, 1243, 1307, 1342, 1343, 1371, 1583, 1923, 1930, 2013},
            "Tautomer-aware MCS path hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computePathFingerprint(acid, 3, 1024))
            == std::vector<int>{67, 92, 167, 181, 199, 383, 522, 602, 627, 640, 663, 666, 689, 727,
                                766, 779, 827, 906, 910, 958, 990, 1019},
            "Labelled path hashes must match the golden bits");
    require(bitPositions(smsd::batch::detail::computeTopologicalTorsion(acid, 2048))
            == std::vector<int>{569, 735, 1337, 1551, 1593, 1728}, "Torsion hashes must match the golden bits");
    // A path read from either end sets the same bit, so a chain and its reversal agree.
    require(smsd::batch::detail::computeMCSFingerprint(smsd::parseSMILES("NCCO"), 7, 2048)
            == smsd::batch::detail::computeMCSFingerprint(smsd::parseSMILES("OCCN"), 7, 2048),
            "MCS path hashes must not depend on atom order");
    require(smsd::batch::detail::computePathFingerprint(smsd::parseSMILES("NCC=O"), 7, 1024)
            == smsd::batch::detail::computePathFingerprint(smsd::parseSMILES("O=CCN"), 7, 1024),
            "Labelled path hashes must not depend on atom order");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"extent", finiteExtent}, {"torsions", aromaticTorsions}, {"sizes", validSizes},
        {"quality-arguments", qualityArguments},
        {"goldens", fixedRadiusGoldens}, {"path-goldens", pathGoldens}});
}

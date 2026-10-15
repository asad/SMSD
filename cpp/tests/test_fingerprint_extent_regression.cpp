/* SPDX-License-Identifier: Apache-2.0 */
#include "smsd/batch.hpp"
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

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
            == std::vector<int>{452, 1828, 1861}, "Radius-zero ECFP hashes must stay unchanged");
    require(bitPositions(smsd::batch::detail::computeCircularFingerprintECFP(graph, 2, 2048))
            == std::vector<int>{273, 452, 501, 1237, 1670, 1702, 1746, 1828, 1861},
            "Fixed-radius ECFP hashes must stay unchanged");
    require(bitPositions(smsd::batch::detail::computeCircularFingerprintFCFP(graph, 2, 2048))
            == std::vector<int>{357, 395, 814, 1172, 1716, 1838, 1933, 1942, 2027},
            "Fixed-radius FCFP hashes must stay unchanged");
    require(bitPositions(smsd::batch::detail::computeTopologicalTorsion(smsd::parseSMILES("C1CCCCC1"), 2048))
            == std::vector<int>{1177}, "Saturated torsion hashes must stay unchanged");
}

int main(int argc, char** argv) {
    const std::string selected = argc > 1 ? argv[1] : "all";
    int failures = 0;
    for (const auto& group : {std::pair<const char*, void(*)()>{"extent", finiteExtent},
                              {"torsions", aromaticTorsions}, {"sizes", validSizes},
                              {"goldens", fixedRadiusGoldens}}) {
        if (selected != "all" && selected != group.first) continue;
        try { group.second(); std::cout << "PASS " << group.first << '\n'; }
        catch (const std::exception& error) { std::cerr << "FAIL " << group.first << ": " << error.what() << '\n'; ++failures; }
    }
    return failures == 0 ? 0 : 1;
}

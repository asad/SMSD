/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"
#include "smsd/smarts_parser.hpp"
#include <iostream>
#include <numeric>

int main() {
    int failures = 0;
    const auto require = [&](bool condition, const std::string& message) {
        if (!condition) { std::cerr << message << '\n'; ++failures; }
    };
    const auto disconnected = smsd::parseSMILES("CC.O.CC");
    const auto matched = smsd::Substructure::find(disconnected).match(disconnected);
    require(matched.status == smsd::SearchStatus::Completed && matched.matches.size() == 1
            && matched.matches.front().fragments == std::vector<std::vector<int>>{{0, 1}, {2}, {3, 4}},
            "Matched fragments must retain ordered query atoms and isolated atoms");
    for (const auto& text : {"[SiH4]", "C[SiH2]C", "[n-]1cccc1", "[2H]", "[CH3]",
                            "[nH]1cccc1", "c1ncc[nH]1", "[NH4+]", "[nH+]1ccccc1"}) {
        const auto original = smsd::parseSMILES(text);
        const auto output = smsd::toSMILES(original);
        const auto restored = smsd::parseSMILES(output);
        require(std::accumulate(original.hydrogenCount.begin(), original.hydrogenCount.end(), 0)
                == std::accumulate(restored.hydrogenCount.begin(), restored.hydrogenCount.end(), 0),
                std::string("SMILES must preserve stored H: ") + text + " -> " + output);
        smsd::SmartsWriteOptions options;
        options.includeHCount = true;
        const auto pattern = smsd::toSMARTS(original, options);
        require(smsd_test::patternMatches(smsd::parseSMARTS(pattern), original),
                std::string("H-qualified SMARTS must match its source: ") + text + " -> " + pattern);
    }
    for (const auto& text : {"c1ncc[nH]1", "c1ccc2[nH]ccc2c1", "O=c1cc[nH]c(=O)[nH]1", "Cn1cccc1", "[CH3]"}) {
        const auto original = smsd::parseSMILES(text);
        std::map<int, int> mapping;
        for (int atom = 0; atom < original.n; ++atom) mapping[atom] = atom;
        const auto output = smsd_test::fragmentSmiles(original, mapping);
        auto restored = smsd::parseSMILES(output);
        require(restored.kekulize(), std::string("Complete MCS must be kekulisable: ") + output);
        require(std::accumulate(original.hydrogenCount.begin(), original.hydrogenCount.end(), 0)
                == std::accumulate(restored.hydrogenCount.begin(), restored.hydrogenCount.end(), 0),
                std::string("Complete MCS must preserve H: ") + text + " -> " + output);
    }
    const auto phenyl = smsd_test::fragmentSmiles(smsd::parseSMILES("Cc1ccccc1"),
                                         {{1, 0}, {2, 1}, {3, 2}, {4, 3}, {5, 4}, {6, 5}});
    const auto benzene = smsd::parseSMILES(phenyl);
    require(std::accumulate(benzene.aromatic.begin(), benzene.aromatic.end(), 0) == 6
            && std::accumulate(benzene.hydrogenCount.begin(), benzene.hydrogenCount.end(), 0) == 6,
            "A retained phenyl ring must be aromatic benzene: " + phenyl);
    const auto naphthalene = smsd::parseSMILES("c1ccc2ccccc2c1");
    std::map<int, int> retainedRing;
    for (int atom : naphthalene.computeRings().front()) retainedRing[atom] = atom;
    const auto ringOutput = smsd_test::fragmentSmiles(naphthalene, retainedRing);
    const auto retainedBenzene = smsd::parseSMILES(ringOutput);
    require(std::accumulate(retainedBenzene.aromatic.begin(), retainedBenzene.aromatic.end(), 0) == 6
            && std::accumulate(retainedBenzene.hydrogenCount.begin(), retainedBenzene.hydrogenCount.end(), 0) == 6,
            "A retained fused ring must remain aromatic benzene: " + ringOutput);
    const auto partial = smsd_test::fragmentSmiles(smsd::parseSMILES("c1ccccc1"),
                                          {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}});
    const auto chain = smsd::parseSMILES(partial);
    require(std::accumulate(chain.aromatic.begin(), chain.aromatic.end(), 0) == 0,
            "An opened aromatic ring must not contain aromatic chain atoms: " + partial);
    for (int atom = 0; atom < chain.n; ++atom) {
        int valence = chain.hydrogenCount[atom];
        for (int neighbor : chain.neighbors[atom]) valence += chain.bondOrder(atom, neighbor);
        require(valence == 4, "An opened aromatic fragment must have capped carbon valences: " + partial);
    }
    const auto ethanol = smsd_test::fragmentSmiles(smsd::parseSMILES("CCCO"), {{1, 0}, {2, 1}, {3, 2}});
    const auto alcohol = smsd::parseSMILES(ethanol);
    require(std::accumulate(alcohol.hydrogenCount.begin(), alcohol.hydrogenCount.end(), 0) == 6,
            "A cut alkyl bond must be capped with H: " + ethanol);
    for (const auto& text : {"F[C@H](Cl)CCO", "C/C=C/C"}) {
        const auto original = smsd::parseSMILES(text);
        std::map<int, int> mapping;
        for (int atom = 0; atom < original.n; ++atom) mapping[atom] = atom;
        const auto output = smsd_test::fragmentSmiles(original, mapping);
        require(output.find(text[0] == 'F' ? '@' : '/') != std::string::npos
                || (text[0] != 'F' && output.find('\\') != std::string::npos),
                "Complete MCS must preserve defined stereo: " + output);
    }
    const auto inferred = smsd::MolGraph::Builder().atomCount(1).atomicNumbers({6}).setNeighbors({{}}).build();
    require(smsd::toSMILES(inferred) == "C", "Omitted Builder H counts must retain the valence fallback");
    require(smsd::toSMILES(smsd::extractSubgraph(inferred, {0})) == "C",
            "Subgraph extraction must preserve the unspecified-H fallback");
    const auto explicitZero = smsd::MolGraph::Builder().atomCount(1).atomicNumbers({6})
        .hydrogenCounts({0}).setNeighbors({{}}).build();
    require(smsd::toSMILES(explicitZero) == "[C]", "Explicit H0 must remain distinct from omitted H data");
    return failures == 0 ? 0 : 1;
}

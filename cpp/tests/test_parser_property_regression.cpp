/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/mol_reader.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using smsd_test::require;

static void testOrganicAtomTokens() {
    const std::vector<std::pair<std::string, std::string>> examples = {
        {"Cn1cccc1", "Cn1cccc1"}, {"CSc1ccccc1", "CSc1ccccc1"},
        {"Ca", "Cc1ccccc1"}, {"Na", "Nc1ccccc1"},
        {"[#6]Na", "CNc1ccccc1"}, {"Ac1ccccc1", "Cc1ccccc1"},
        {"ClCBr", "ClCBr"}, {"[Na+]", "[Na+]"},
        {"as", "c1ccsc1"}, {"[se]", "c1cc[se]c1"}
    };
    for (const auto& example : examples) {
        const auto target = smsd::parseSMILES(example.second);
        require(smsd_test::patternMatches(smsd::parseSMARTS(example.first), target),
                example.first + " should match " + example.second);
    }
}

static void testAromaticHydrogens() {
    const std::vector<std::pair<std::string, std::vector<int>>> examples = {
        {"c1ccsc1", {1, 1, 1, 0, 1}},
        {"c1ccoc1", {1, 1, 1, 0, 1}},
        {"Cn1cccc1", {3, 0, 1, 1, 1, 1}},
        {"n1ccccc1", {0, 1, 1, 1, 1, 1}},
        {"[nH]1cccc1", {1, 1, 1, 1, 1}}
    };
    for (const auto& example : examples) {
        const auto target = smsd::parseSMILES(example.first);
        require(target.hydrogenCount == example.second,
                "incorrect aromatic H count for " + example.first);
    }
    require(smsd_test::allPatternAtoms(smsd::parseSMARTS("[sH]"), smsd::parseSMILES("c1ccsc1")).empty(),
            "thiophene sulphur must not have hydrogen");
    require(smsd_test::allPatternAtoms(smsd::parseSMARTS("[nH]"),
            smsd::parseSMILES("Cn1cnc2c1c(=O)n(C)c(=O)n2C")).empty(),
            "caffeine nitrogens must not have hydrogen");
}

static void testHydrogenAndValencePrimitives() {
    const std::vector<std::pair<std::string, std::string>> examples = {
        {"[OH1]", "[H]OC"}, {"[CH3]", "[H]C([H])([H])O[H]"},
        {"[OX2H]", "[H]OC"}, {"[N;!H0]", "[H]N([H])C"},
        {"[N;X3;H2;v3]", "[H]N([H])C"},
        {"[C;v4]", "CC"}, {"[N;v3]", "CN"},
        {"[c;v4]", "c1ccccc1"}, {"[c;v4;X3]", "c1ccc2ccccc2c1"},
        {"[n;v3]", "n1ccccc1"}, {"[nH;v3]", "[nH]1cccc1"},
        {"[n;v3;H0]", "Cn1cccc1"}, {"[s;v2;H0]", "c1ccsc1"},
        {"[C;v3]", "[CH3]"}, {"[N+;v4;H4]", "[NH4+]"},
        {"[C;v4;X4;H4]", "[H]C([H])([H])[H]"}
    };
    for (const auto& example : examples) {
        const auto target = smsd::parseSMILES(example.second);
        require(smsd_test::patternMatches(smsd::parseSMARTS(example.first), target),
                example.first + " should match " + example.second);
    }
    const auto methanol = smsd::parseSMILES("[H]C([H])([H])O[H]");
    const auto molRoundTrip = smsd::readMolBlock(smsd::writeMolBlock(methanol));
    require(smsd_test::patternMatches(smsd::parseSMARTS("[CH3;v4;X4][OH1;v2;X2]"), molRoundTrip),
            "explicit-H MOL methanol should retain total H, valence and connectivity");
}

static void testLazyCanonicalProperties() {
    const auto carbon = smsd::parseSMILES("C");
    const auto oxygen = smsd::parseSMILES("O");
    require(carbon.getCanonicalHash() != oxygen.getCanonicalHash(),
            "fresh different molecules should have different hashes");
    const auto molecule = smsd::parseSMILES("CCO");
    auto labels = molecule.getCanonicalLabeling();
    std::sort(labels.begin(), labels.end());
    require(labels == std::vector<int>({0, 1, 2}),
            "fresh canonical labels should be an atom permutation");
    const auto symmetric = smsd::parseSMILES("CC");
    const auto& orbit = symmetric.getOrbit();
    require(orbit.size() == 2 && orbit[0] == orbit[1],
            "fresh equivalent atoms should share an orbit");
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"organic-atom-tokens", testOrganicAtomTokens}, {"aromatic-hydrogens", testAromaticHydrogens},
        {"hydrogen-and-valence-primitives", testHydrogenAndValencePrimitives},
        {"lazy-canonical-properties", testLazyCanonicalProperties}});
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/cip.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>

using smsd_test::require;

static void aromaticCarriers() {
    const auto first = smsd::parseSMILES("F/C=c1\\cnc2c(n1)N=NNC=2");
    const auto second = smsd::parseSMILES("F/C=c1/cnc2c(n1)N=NNC=2");
    smsd::SubstructureOptions stereo;
    stereo.chemistry.useBondStereo = true;
    require(first.getCanonicalHash() != second.getCanonicalHash(), "The two exocyclic configurations must be distinct");
    for (const auto& graph : {first, second}) {
        require(graph.dbStereo(1, 2) != 0 && graph.bondOrder(0, 1) == 1 && graph.bondOrder(2, 3) == 4,
                "The stereo double bond must have single and aromatic reference carriers");
        const auto text = smsd::toSMILES(graph);
        require(text == graph.toCanonicalSmiles(), "The member and free SMILES writers must agree");
        const auto restored = smsd::parseSMILES(text);
        require(restored.getCanonicalHash() == graph.getCanonicalHash(), "Aromatic-carrier output must retain the molecular hash");
        const auto match = smsd::Substructure::find(graph).withOptions(stereo).match(restored);
        require(match.optimal && match.complete && !match.matches.empty(), "Aromatic-carrier output must retain mapped double-bond stereo");
        smsd::MCSOptions validation;
        validation.chemistry.useBondStereo = true;
        require(smsd::MCS::validateMapping(graph, restored, match.matches.front(), validation).empty(), "Invalid mapped stereo round trip");
    }
    const auto opposite = smsd::Substructure::find(first).withOptions(stereo).match(second);
    require(opposite.optimal && opposite.complete && opposite.matches.empty(), "Opposite exocyclic stereo must not match");
    const auto equivalent = smsd::parseSMILES("F\\C=c1/cnc2c(n1)N=NNC=2");
    require(smsd::toSMILES(first) == smsd::toSMILES(equivalent), "Equivalent directional assignments must have one canonical output");
}

static void singleCarriersAtAromaticFocus() {
    // '/' and '\' are single bonds in OpenSMILES, so the configuration of an aromatic-focus double
    // bond is written on the single bonds to its substituents rather than on the ring bonds. A ring
    // double bond can be trans only in a ring of eight or more atoms, as in [10]annulene.
    const std::pair<const char*, const char*> cases[] = {
        {"F/C1=C(\\Cl)C=CC=C1", "[c]1(=[c](cccc1)/F)/Cl"},
        {"F/C1=C(\\Cl)C=CC=CC=CC=C1", "[c]1(=[c](cccccccc1)/F)/Cl"},
        {"F/C1=C(/Cl)C=CC=CC=CC=C1", "[c]1(=[c](cccccccc1)\\F)/Cl"}};
    for (const auto& item : cases) {
        const auto graph = smsd::parseSMILES(item.first);
        const auto text = smsd::toSMILES(graph);
        require(text == item.second, std::string("Directional bonds must be on single bonds: ") + text);
        require(smsd::parseSMILES(text).getCanonicalHash() == graph.getCanonicalHash(),
                "Single-carrier aromatic-focus stereo must round trip");
    }
}

static void ordinaryConjugation() {
    for (const char* input : {"F/C=C/F", "F/C=C\\F", "F/C=C/C=C/F", "F/C=C\\C=C/F",
                              "F/C=C/C(=C/F)Cl", "F/C=C/C(=C\\F)Cl"}) {
        const auto graph = smsd::parseSMILES(input);
        const auto restored = smsd::parseSMILES(smsd::toSMILES(graph));
        require(graph.getCanonicalHash() == restored.getCanonicalHash(), "Ordinary conjugated stereo must round trip");
    }
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
    auto copy = smsd::MolGraph::Builder().atomCount(graph.n).atomicNumbers(std::move(numbers))
        .formalCharges(std::move(charges)).massNumbers(std::move(masses))
        .hydrogenCounts(std::move(hydrogens)).atomClasses(std::move(classes))
        .aromaticFlags(std::move(aromatic)).ringFlags(std::move(rings))
        .setNeighbors(std::move(neighbors)).setBondOrders(std::move(orders))
        .bondAromaticFlags(std::move(aromaticBonds)).bondRingFlags(std::move(ringBonds)).build(false);
    for (int first = 0; first < graph.n; ++first) {
        for (int second : graph.neighbors[first]) {
            if (first >= second || graph.dbStereo(first, second) == 0) continue;
            int configuration = smsd::MolGraph::canonicalDoubleBondTag(graph, first, second, inverse);
            for (const auto& end : {std::make_pair(inverse[first], inverse[second]),
                                    std::make_pair(inverse[second], inverse[first])}) {
                int reference = -1, minimum = graph.n;
                for (int neighbor : copy.neighbors[end.first]) {
                    if (neighbor == end.second) continue;
                    if (reference < 0) reference = neighbor;
                    minimum = std::min(minimum, neighbor);
                }
                if (reference != minimum) configuration = 3 - configuration;
            }
            copy.dbStereoConf[inverse[first]][inverse[second]] = configuration;
            copy.dbStereoConf[inverse[second]][inverse[first]] = configuration;
        }
    }
    return copy;
}

static void aromaticFocus() {
    const char* forms[] = {
        "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1",
        "C1=CC2=N/C1=C/C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1"
    };
    const auto first = smsd::parseSMILES(forms[0]);
    const auto opposite = smsd::parseSMILES(forms[1]);
    require(first.getCanonicalHash() != opposite.getCanonicalHash(), "Opposite aromatic-focus stereo must remain distinct");
    smsd::SubstructureOptions stereo;
    stereo.chemistry.useBondStereo = true;
    require(smsd::Substructure::find(first).withOptions(stereo).match(opposite).matches.empty(),
            "Opposite aromatic-focus configurations must not match");
    std::mt19937 random(4179);
    for (const char* form : forms) {
        const auto graph = smsd::parseSMILES(form);
        const auto expected = smsd::toSMILES(graph);
        const auto restored = smsd::parseSMILES(expected);
        require(restored.getCanonicalHash() == graph.getCanonicalHash(), "Aromatic-focus stereo must round trip without extra tags");
        require(!smsd::Substructure::find(graph).withOptions(stereo).match(restored).matches.empty(),
                "Aromatic-focus output must retain mapped stereo");
        for (int trial = 0; trial < 12; ++trial) {
            const auto copy = reordered(graph, random);
            require(copy.getCanonicalHash() == graph.getCanonicalHash(), "Permutation changed aromatic-focus stereo");
            require(smsd::toSMILES(copy) == expected, "Permutation changed canonical aromatic-focus output");
        }
    }
}

static void conflictingFocus() {
    std::vector<std::vector<int>> neighbors(10), orders(10), stereo(10, std::vector<int>(10, 0));
    for (int atom = 0; atom < 10; ++atom) {
        const int next = (atom + 1) % 10;
        neighbors[atom].push_back(next);
        neighbors[next].push_back(atom);
        orders[atom].push_back(4);
        orders[next].push_back(4);
    }
    stereo[0][1] = stereo[1][0] = stereo[1][2] = stereo[2][1] = 1;
    const auto graph = smsd::MolGraph::Builder().atomCount(10).atomicNumbers(std::vector<int>(10, 6))
        .hydrogenCounts(std::vector<int>(10, 1)).aromaticFlags(std::vector<uint8_t>(10, 1))
        .ringFlags(std::vector<uint8_t>(10, 1)).setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(orders)).doubleBondStereo(std::move(stereo)).build(false);
    bool rejected = false;
    try { smsd::toSMILES(graph); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "A stereo focus cannot also be a directional reference carrier");
    for (const bool v3000 : {false, true}) {
        rejected = false;
        try {
            if (v3000) smsd::writeMolBlockV3000(graph);
            else smsd::writeMolBlock(graph);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "MOL output must reject incompatible aromatic-focus metadata");
    }
}

static void aromaticMolStereo() {
    const auto graph = smsd::parseSMILES(
        "C1=CC2=N/C1=C\\C1=N/C(=C\\C3=NC4=C(C5=N/C(=C\\2)CC5)CC(=O)C4=C3)C=C1");
    auto stereo = smsd::cip::assignAll(graph).ezBonds;
    stereo.erase(std::remove_if(stereo.begin(), stereo.end(), [&](const auto& bond) {
        return graph.atomicNum[std::get<0>(bond)] != 6 || graph.atomicNum[std::get<1>(bond)] != 6;
    }), stereo.end());
    require(stereo.size() == 3, "The macrocycle must contain three specified carbon double-bond configurations");
    for (const auto& block : {smsd::writeMolBlock(graph), smsd::writeMolBlockV3000(graph)}) {
        const auto restored = smsd::readMolBlock(block);
        require(restored.n == graph.n, "MOL output changed the macrocycle atom indices");
        for (const auto& bond : stereo) {
            const int first = std::get<0>(bond), second = std::get<1>(bond);
            require(smsd::MolGraph::canonicalDoubleBondTag(graph, first, second, {})
                    == smsd::MolGraph::canonicalDoubleBondTag(restored, first, second, {}),
                    "MOL coordinates changed a specified aromatic-focus configuration");
        }
    }
}

static std::vector<smsd::Point2D> molCoordinates(const std::string& block, int atoms) {
    std::vector<smsd::Point2D> positions(atoms);
    std::istringstream input(block);
    std::string line;
    const bool v3000 = block.find("V3000") != std::string::npos;
    bool inAtoms = false;
    int lineNumber = 0, count = 0;
    while (std::getline(input, line)) {
        if (v3000) {
            if (line == "M  V30 BEGIN ATOM") inAtoms = true;
            else if (line == "M  V30 END ATOM") inAtoms = false;
            else if (inAtoms) {
                std::istringstream fields(line);
                std::string marker, version, symbol;
                int atom;
                double x, y;
                fields >> marker >> version >> atom >> symbol >> x >> y;
                positions.at(atom - 1) = {x, y};
                ++count;
            }
        } else if (lineNumber >= 4 && lineNumber < 4 + atoms) {
            std::istringstream fields(line);
            fields >> positions[lineNumber - 4].x >> positions[lineNumber - 4].y;
            ++count;
        }
        ++lineNumber;
    }
    require(count == atoms, "MOL output has incomplete atom coordinates");
    return positions;
}

static void cyclicLigandGeometry() {
    for (const char* text : {"C/C1=C\\CCCCCC1", "C/C1=C/CCCCCC1"}) {
        const auto graph = smsd::parseSMILES(text);
        const auto hash = graph.getCanonicalHash();
        for (const auto& block : {smsd::writeMolBlock(graph), smsd::writeMolBlockV3000(graph)}) {
            const auto positions = molCoordinates(block, graph.n);
            const double axisX = positions[2].x - positions[1].x;
            const double axisY = positions[2].y - positions[1].y;
            auto side = [&](int ligand) {
                return axisX * (positions[ligand].y - positions[1].y)
                     - axisY * (positions[ligand].x - positions[1].x);
            };
            require(side(0) * side(8) < 0, "MOL double-bond ligands at one end must lie on opposite sides");
            const auto restored = smsd::readMolBlock(block);
            require(smsd::cip::assignEZ(restored, 1, 2) == smsd::cip::assignEZ(graph, 1, 2),
                    "Cyclic MOL geometry changed the specified E/Z configuration");
        }
        require(graph.getCanonicalHash() == hash, "MOL coordinate generation changed the source graph");
    }
}

static void unspecifiedAndDegenerateGeometry() {
    const auto unspecified = smsd::parseSMILES("CC1=CCCCCCC1");
    for (const auto& block : {smsd::writeMolBlock(unspecified), smsd::writeMolBlockV3000(unspecified)})
        require(smsd::cip::assignAll(smsd::readMolBlock(block)).ezBonds.empty(),
                "MOL output must not invent unspecified double-bond stereo");
    const auto graph = smsd::MolGraph::Builder().atomCount(2).atomicNumbers({6, 6})
        .hydrogenCounts({2, 2}).setNeighbors({{1}, {0}}).setBondOrders({{2}, {2}})
        .doubleBondStereo({{0, 1}, {1, 0}}).build(false);
    for (bool v3000 : {false, true}) {
        bool rejected = false;
        try {
            if (v3000) smsd::writeMolBlockV3000(graph);
            else smsd::writeMolBlock(graph);
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "A reference-less stereo focus must not produce ambiguous MOL coordinates");
    }
}

static std::string cyclicMolBlock(int atoms, int shortcut = -1, bool alternating = false) {
    std::ostringstream block;
    block << "Cycloalkene\n  SMSD              2D\n\n"
          << std::setw(3) << atoms << std::setw(3) << atoms + (shortcut >= 0 ? 1 : 0)
          << "  0  0  0  0  0  0  0  0999 V2000\n";
    for (int atom = 0; atom < atoms; ++atom) {
        const double angle = 2 * std::acos(-1.0) * atom / atoms;
        block << std::fixed << std::setprecision(4)
              << std::setw(10) << std::cos(angle) << std::setw(10) << std::sin(angle)
              << std::setw(10) << 0.0 << " C   0  0  0  0  0  0  0  0  0  0  0  0\n";
    }
    for (int atom = 0; atom < atoms; ++atom)
        block << std::setw(3) << atom + 1 << std::setw(3) << (atom + 1) % atoms + 1
              << std::setw(3) << (atom == 0 || (alternating && atom < 6 && atom % 2 == 0) ? 2 : 1) << "  0\n";
    if (shortcut >= 0) block << "  1" << std::setw(3) << shortcut + 1 << "  1  0\n";
    return block.str() + "M  END\n";
}

static void inferredRingStereo() {
    require(smsd::readMolBlock(cyclicMolBlock(6, -1, true)).dbStereo(0, 1) == 0,
            "Small aromatic ring coordinates must not create double-bond stereo");
    require(smsd::readMolBlock(cyclicMolBlock(7)).dbStereo(0, 1) == 0,
            "Nonaromatic seven-membered ring geometry must not create double-bond stereo");
    require(smsd::readMolBlock(cyclicMolBlock(8)).dbStereo(0, 1) != 0,
            "Eight-membered ring geometry must retain double-bond stereo");
    require(smsd::readMolBlock(cyclicMolBlock(8, 5, true)).dbStereo(0, 1) == 0,
            "A shorter bridged cycle must prevent inferred double-bond stereo");
    require(smsd::readMolBlock(cyclicMolBlock(10, 7)).dbStereo(0, 1) != 0,
            "A small ring away from the stereo focus must not erase its geometry");
    const auto supplied = smsd::parseSMILES("C1/C=C\\CCCC1");
    require(supplied.hasDbStereo, "Supplied small-ring stereo must remain available");
}

int main() {
    try {
        aromaticCarriers();
        ordinaryConjugation();
        singleCarriersAtAromaticFocus();
        inferredRingStereo();
        aromaticFocus();
        aromaticMolStereo();
        cyclicLigandGeometry();
        unspecifiedAndDegenerateGeometry();
        conflictingFocus();
        std::cout << "PASS aromatic reference carriers and conjugated double-bond stereo\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}

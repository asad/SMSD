/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/cip.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/smiles_parser.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>

using smsd::cip::EZLabel;
using smsd::cip::RSLabel;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static std::vector<RSLabel> labels(const smsd::MolGraph& graph) {
    auto result = smsd::cip::assignAll(graph).rsLabels;
    result.erase(std::remove(result.begin(), result.end(), RSLabel::NONE), result.end());
    std::sort(result.begin(), result.end());
    return result;
}

static std::string alanineMol(bool v3000, bool reverseBonds = false) {
    if (v3000) return
        "Alanine\n\n\n  0  0  0  0  0  0            999 V3000\n"
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 6 5 0 0 0\nM  V30 BEGIN ATOM\n"
        "M  V30 1 N 0 0 0 0\nM  V30 2 C 0 0 0 0 CFG=1\n"
        "M  V30 3 C 0 0 0 0\nM  V30 4 C 0 0 0 0\n"
        "M  V30 5 O 0 0 0 0\nM  V30 6 O 0 0 0 0\nM  V30 END ATOM\n"
        "M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\n"
        "M  V30 3 1 2 4\nM  V30 4 2 4 5\nM  V30 5 1 4 6\n"
        "M  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    std::string result = "Alanine\n\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n";
    const char* symbols[] = {"N", "C", "C", "C", "O", "O"};
    char line[160];
    for (int atom = 0; atom < 6; ++atom) {
        std::snprintf(line, sizeof(line),
            "    0.0000    0.0000    0.0000 %-3s 0  0%3d  0  0  0  0  0  0  0  0  0\n",
            symbols[atom], atom == 1 ? 1 : 0);
        result += line;
    }
    result += reverseBonds ? "  2  3  1  0\n  1  2  1  0\n" : "  1  2  1  0\n  2  3  1  0\n";
    return result + "  2  4  1  0\n  4  5  2  0\n  4  6  1  0\nM  END\n";
}

static std::string alkeneMol(bool v3000, bool either) {
    if (v3000) return std::string(
        "Butene\n\n\n  0  0  0  0  0  0            999 V3000\n"
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 4 3 0 0 0\nM  V30 BEGIN ATOM\n"
        "M  V30 1 C -1 1 0 0\nM  V30 2 C 0 0 0 0\n"
        "M  V30 3 C 1 0 0 0\nM  V30 4 C 2 -1 0 0\nM  V30 END ATOM\n"
        "M  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 2 2 3") +
        (either ? " CFG=2" : "") + "\nM  V30 3 1 3 4\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    return std::string("Butene\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "   -1.0000    1.0000    0.0000 C   0  0  0\n"
        "    0.0000    0.0000    0.0000 C   0  0  0\n"
        "    1.0000    0.0000    0.0000 C   0  0  0\n"
        "    2.0000   -1.0000    0.0000 C   0  0  0\n  1  2  1  0\n  2  3  2") +
        (either ? "  3" : "  0") + "\n  3  4  1  0\nM  END\n";
}

int main() {
    int failures = 0;
    auto run = [&](const char* name, const std::function<void()>& test) {
        try { test(); std::cout << "PASS " << name << '\n'; }
        catch (const std::exception& error) { ++failures; std::cout << "FAIL " << name << ": " << error.what() << '\n'; }
    };
    run("tetrahedral SMILES round-trip", [] {
        for (const char* text : {"C[C@@H](N)C(=O)O", "C[C@H](N)C(=O)O", "CC[C@@H](O)C",
                "[C@@H](Br)(Cl)F", "N[C@]1(C)CCCCO1", "C[C@H]1CCCCO1", "CC1.[C@H]1(N)O"}) {
            const auto graph = smsd::parseSMILES(text);
            require(labels(graph) == labels(smsd::parseSMILES(smsd::toSMILES(graph))), "SMILES changed the stereocentre");
        }
    });
    run("tetrahedral class", [] {
        require(labels(smsd::parseSMILES("N[C@TH2H](C)C(=O)O")) ==
                labels(smsd::parseSMILES("N[C@@H](C)C(=O)O")), "@TH2 must mean @@");
        require(labels(smsd::parseSMILES("N[C@TH1H](C)C(=O)O")) ==
                labels(smsd::parseSMILES("N[C@H](C)C(=O)O")), "@TH1 must mean @");
        bool rejected = false;
        try { smsd::parseSMILES("[Pt@SP1](Cl)(F)(Br)I"); } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "unsupported non-tetrahedral stereo must be rejected");
    });
    run("branch directional stereo", [] {
        for (const char* text : {"F/C=C/F", "C(\\F)=C/F", "C(=C/Cl)\\Cl"}) {
            const auto graph = smsd::parseSMILES(text);
            const auto stereo = smsd::cip::assignAll(graph).ezBonds;
            require(stereo.size() == 1 && std::get<2>(stereo[0]) == EZLabel::E, "trans branch alkene must be E");
        }
        const auto graph = smsd::parseSMILES("F/C=C(Cl)/C");
        require(std::get<2>(smsd::cip::assignAll(graph).ezBonds[0]) == EZLabel::Z, "stereo reference must follow the annotated ligand");
    });
    run("double-bond SMILES round-trip", [] {
        for (const char* text : {"C/C=C/C", "C/C=C\\C", "C(=C/Cl)\\Cl", "F/C=C(Cl)/C", "C/C=C/C=C\\C", "C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C", "C1/C=C\\CCCC1"}) {
            auto graph = smsd::parseSMILES(text);
            for (int atom = 0; atom < graph.n; ++atom) graph.atomClass[atom] = atom + 1;
            const auto restored = smsd::parseSMILES(smsd::toSMILES(graph));
            auto stereoByClass = [](const smsd::MolGraph& molecule) {
                std::map<std::pair<int, int>, EZLabel> result;
                for (const auto& bond : smsd::cip::assignAll(molecule).ezBonds) {
                    const int first = molecule.atomClass[std::get<0>(bond)];
                    const int second = molecule.atomClass[std::get<1>(bond)];
                    result[{std::min(first, second), std::max(first, second)}] = std::get<2>(bond);
                }
                return result;
            };
            require(stereoByClass(graph) == stereoByClass(restored), "SMILES changed a double-bond configuration");
        }
    });
    run("known MOL tetrahedral parity", [] {
        for (bool v3000 : {false, true}) {
            const auto graph = smsd::readMolBlock(alanineMol(v3000));
            require(smsd::cip::assignRS(graph, 1) == RSLabel::S, "MDL parity 1 alanine must be S");
        }
        require(smsd::cip::assignRS(smsd::readMolBlock(alanineMol(false, true)), 1) == RSLabel::S,
                "MDL parity must not depend on bond line order");
    });
    run("MOL wedge without atom parity", [] {
        const std::string block =
            "Alanine\n\n\n  6  5  0  0  0  0  0  0  0  0999 V2000\n"
            "    1.2990    2.2500    0.0000 N   0  0  0\n"
            "    1.2990    0.7500    0.0000 C   0  0  0\n"
            "    0.0000    0.0000    0.0000 C   0  0  0\n"
            "    2.5981    0.0000    0.0000 C   0  0  0\n"
            "    2.5981   -1.5000    0.0000 O   0  0  0\n"
            "    3.8971    0.7500    0.0000 O   0  0  0\n"
            "  1  2  1  0\n  2  3  1  6\n  2  4  1  0\n  4  5  2  0\n  4  6  1  0\nM  END\n";
        require(smsd::cip::assignRS(smsd::readMolBlock(block), 1) == RSLabel::S, "hashed wedge alanine must be S");
    });
    run("MOL tetrahedral writer", [] {
        for (const char* text : {"N[C@@H](C)C(=O)O", "[C@@H](Br)(Cl)F", "[H][C@](Br)(Cl)F", "C/C=C/[C@H](F)Cl"}) {
            const auto graph = smsd::parseSMILES(text);
            require(labels(graph) == labels(smsd::readMolBlock(smsd::writeMolBlock(graph))), "V2000 changed tetrahedral stereo");
            require(labels(graph) == labels(smsd::readMolBlock(smsd::writeMolBlockV3000(graph))), "V3000 changed tetrahedral stereo");
            std::istringstream input(smsd::writeMolBlock(graph));
            std::string block, line;
            for (int lineIndex = 0; std::getline(input, line); ++lineIndex) {
                if (lineIndex >= 4 && lineIndex < 4 + graph.n) line.replace(39, 3, "  0");
                block += line + "\n";
            }
            require(labels(graph) == labels(smsd::readMolBlock(block)), "V2000 wedge must preserve tetrahedral stereo without atom parity");
        }
    });
    run("MOL crossed double bond", [] {
        for (bool v3000 : {false, true}) {
            const auto graph = smsd::readMolBlock(alkeneMol(v3000, true));
            require(smsd::cip::assignEZ(graph, 1, 2) == EZLabel::NONE, "either must remain unspecified");
        }
    });
    run("MOL double-bond coordinates", [] {
        for (bool v3000 : {false, true}) {
            const auto graph = smsd::readMolBlock(alkeneMol(v3000, false));
            require(smsd::cip::assignEZ(graph, 1, 2) == EZLabel::E, "coordinates must encode trans butene");
        }
    });
    run("MOL double-bond writer", [] {
        for (const char* text : {"C/C=C/C", "C/C=C\\C", "C/C=C/C=C\\C", "C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C", "C1/C=C\\CCCC1"}) {
            const auto graph = smsd::parseSMILES(text);
            for (const auto& block : {smsd::writeMolBlock(graph), smsd::writeMolBlockV3000(graph)}) {
                std::istringstream input(block);
                std::string line;
                const bool v3000 = block.find("V3000") != std::string::npos;
                bool inBonds = false;
                int bondCount = 0;
                for (const auto& neighbors : graph.neighbors) bondCount += static_cast<int>(neighbors.size());
                bondCount /= 2;
                int lineIndex = 0;
                while (std::getline(input, line)) {
                    if (v3000) {
                        if (line == "M  V30 BEGIN BOND") inBonds = true;
                        else if (line == "M  V30 END BOND") inBonds = false;
                        else if (inBonds) {
                            std::istringstream fields(line);
                            std::string prefix, version;
                            int index, order;
                            fields >> prefix >> version >> index >> order;
                            if (order == 2) require(line.find("CFG=") == std::string::npos,
                                    "specified E/Z must use coordinates, not V3000 either flags");
                        }
                    } else if (lineIndex >= 4 + graph.n && lineIndex < 4 + graph.n + bondCount) {
                        if (std::stoi(line.substr(6, 3)) == 2)
                            require(std::stoi(line.substr(9, 3)) == 0, "specified E/Z must use coordinates, not V2000 stereo flags");
                    }
                    ++lineIndex;
                }
                const auto written = smsd::readMolBlock(block);
                for (const auto& bond : smsd::cip::assignAll(graph).ezBonds) {
                    const int first = std::get<0>(bond), second = std::get<1>(bond);
                    require(smsd::cip::assignEZ(written, first, second) == std::get<2>(bond), "MOL changed double-bond stereo");
                }
            }
        }
    });
    return failures == 0 ? 0 : 1;
}

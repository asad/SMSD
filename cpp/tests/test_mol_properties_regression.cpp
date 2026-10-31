/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/smarts_parser.hpp"
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <fstream>

static std::string hydrogenBlock(bool v3000, const std::string& attributes = "", const std::string& property = "") {
    if (v3000)
        return "Hydrogen\n\n\n  0  0  0  0  0  0            999 V3000\n"
            "M  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n"
            "M  V30 1 H 0 0 0 0" + attributes + "\nM  V30 END ATOM\n"
            "M  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    std::string atom = "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0";
    if (attributes == "zero-valence") atom.replace(48, 3, " 15");
    if (attributes == "zero-hydrogens") atom.replace(42, 3, "  1");
    if (attributes == "doublet") atom.replace(36, 3, "  4");
    return "Hydrogen\n\n\n  1  0  0  0  0  0            999 V2000\n" + atom + "\n" + property + "M  END\n";
}

int main() {
    for (bool v3000 : {false, true}) {
        if (smsd::readMolBlock(hydrogenBlock(v3000)).hydrogenCount != std::vector<int>{1})
            throw std::runtime_error("An isolated neutral MOL hydrogen must use valence one");
        for (const auto& attributes : {v3000 ? " VAL=-1" : "zero-valence",
                                       v3000 ? " HCOUNT=-1" : "zero-hydrogens",
                                       v3000 ? " RAD=2" : "doublet"})
            if (smsd::readMolBlock(hydrogenBlock(v3000, attributes)).hydrogenCount != std::vector<int>{0})
                throw std::runtime_error("Explicit hydrogen and radical constraints must be retained");
        for (int radical : {1, 2, 3}) {
            const auto block = v3000 ? hydrogenBlock(true, " RAD=" + std::to_string(radical))
                : hydrogenBlock(false, "", "M  RAD  1   1   " + std::to_string(radical) + "\n");
            if (smsd::readMolBlock(block).hydrogenCount != std::vector<int>{0})
                throw std::runtime_error("A radical MOL hydrogen must not gain an implicit hydrogen");
        }
        for (int charge : {-1, 1}) {
            const auto block = v3000 ? hydrogenBlock(true, " CHG=" + std::to_string(charge))
                : hydrogenBlock(false, "", "M  CHG  1   1  " + std::to_string(charge) + "\n");
            if (smsd::readMolBlock(block).hydrogenCount != std::vector<int>{0})
                throw std::runtime_error("A charged MOL hydrogen must not gain an implicit hydrogen");
        }
    }
    for (const char* text : {"CCO", "OCC(=O)O", "c1ccccc1", "[nH]1cccc1", "[NH4+]", "[H]OC", "[CH3]", "[H]", "[HH]", "[H][H]"}) {
        const auto original = smsd::parseSMILES(text);
        for (const auto& block : {smsd::writeMolBlock(original), smsd::writeMolBlockV3000(original)}) {
            const auto restored = smsd::readMolBlock(block);
            if (restored.hydrogenCount != original.hydrogenCount || restored.atomicNum != original.atomicNum
                || block.find("HCOUNT=") != std::string::npos)
                throw std::runtime_error(std::string("MOL atom properties changed: ") + text);
        }
    }
    const auto acid = smsd::readMolBlock(smsd::writeMolBlock(smsd::parseSMILES("OCC(=O)O")));
    if (!smsd_test::patternMatches(smsd::parseSMARTS("[OH1][CH2][CX3](=O)[OX2H1]"), acid))
        throw std::runtime_error("MOL hydrogen-dependent SMARTS must match glycolic acid");
    for (const char* text : {"c1ccc(cc1)c1ccccc1", "c1ccccc1c1ccncc1"}) {
        const auto graph = smsd::parseSMILES(text);
        int bridges = 0;
        for (int first = 0; first < graph.n; ++first)
            for (int second : graph.neighbors[first]) {
                if (second <= first || graph.bondInRing(first, second)) continue;
                ++bridges;
                if (graph.bondOrder(first, second) != 1 || graph.bondAromatic(first, second))
                    throw std::runtime_error("Biaryl bridge must be a single bond");
            }
        if (bridges != 1) throw std::runtime_error("Biaryl must have one inter-ring bridge");
    }
    for (bool v3000 : {false, true}) {
        std::string block = v3000 ? smsd::writeMolBlockV3000(smsd::parseSMILES("[H]O[H]"))
                                  : smsd::writeMolBlock(smsd::parseSMILES("[H]O[H]"));
        if (v3000) {
            block.replace(block.find("M  V30 1 H ") + 9, 1, "D");
            block.replace(block.find("M  V30 3 H ") + 9, 1, "T");
        } else {
            std::istringstream input(block);
            std::string line;
            block.clear();
            for (int lineIndex = 0; std::getline(input, line); ++lineIndex) {
                if (lineIndex == 4 || lineIndex == 6) line[31] = lineIndex == 4 ? 'D' : 'T';
                block += line + '\n';
            }
        }
        if (smsd::readMolBlock(block).massNumber != std::vector<int>({2, 0, 3}))
            throw std::runtime_error("D/T aliases must retain hydrogen isotopes");
    }
    std::string continuation = smsd::writeMolBlockV3000(smsd::parseSMILES("N"));
    const auto atomOffset = continuation.find("M  V30 1 N 0 0 0 0");
    continuation.replace(atomOffset, 18, "M  V30 1 N 0 0 -\nM  V30 0 0 CHG=1");
    if (smsd::readMolBlock(continuation).formalCharge[0] != 1)
        throw std::runtime_error("V3000 continuation must retain atom attributes");
    for (int bondType : {8, 9, 10}) {
        auto block = smsd::writeMolBlockV3000(smsd::parseSMILES("CO"));
        block.replace(block.find("M  V30 1 1 1 2") + 9, 1, std::to_string(bondType));
        const auto restored = smsd::readMolBlock(block);
        if (restored.bondOrder(0, 1) != (bondType == 10 ? 0 : 1))
            throw std::runtime_error("V3000 extended bond type must not be truncated");
    }
    const auto largeGraph = smsd::parseSMILES(std::string(1000, 'C'));
    const auto largeBlock = smsd::writeMolBlock(largeGraph);
    if (largeBlock.find("V3000") == std::string::npos || smsd::readMolBlock(largeBlock).n != 1000)
        throw std::runtime_error("Large MOL output must use V3000");
    if (largeBlock != smsd::writeMolBlockV3000(largeGraph))
        throw std::runtime_error("Automatic V3000 output must match the explicit writer");
    // A unique name lets concurrent runs of this test share the temporary directory.
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto sdfPath = std::filesystem::temp_directory_path()
        / ("smsd-mol-properties-regression-" + std::to_string(suffix) + ".sdf");
    {
        std::ofstream output(sdfPath);
        output << smsd::writeMolBlock(acid) << "$$$$\n\n \t\n";
    }
    const auto records = smsd::readSDF(sdfPath.string());
    std::filesystem::remove(sdfPath);
    if (records.size() != 1 || records[0].n != acid.n)
        throw std::runtime_error("SDF trailing whitespace must not add a record");
    std::cout << "MOL properties and biaryl bonds passed\n";
}

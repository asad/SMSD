/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/smiles_parser.hpp"

#include <iostream>
#include <stdexcept>

void testParserWriterParity() {
    for (const char* text : {"[13CH4]", "[NH4+]", "[C]", "[CH2]", "[CH3:7]",
            "C[C@H](N)C(=O)O", "N[C@@H](C)C(=O)O", "C/C=C/C", "C/C=C\\C",
            "C/C=C(/F)\\C(Cl)=C/C", "C1[C@H](F)CCC[C@@H]1F",
            "C[N@](CC)CCC", "c1ccc(-c2ccccc2)cc1", "[nH]1cccc1", "CCO.[Na+]"}) {
        const auto graph = smsd::parseSMILES(text);
        const auto member = graph.toCanonicalSmiles();
        if (member != smsd::toSMILES(graph))
            throw std::runtime_error(std::string("member and free writer differ: ") + text);
        if (smsd::parseSMILES(member).getCanonicalHash() != graph.getCanonicalHash())
            throw std::runtime_error(std::string("member writer changed molecule: ") + text);
    }
}

int main() {
    try { testParserWriterParity(); std::cout << "PASS parser-parity\n"; }
    catch (const std::exception& error) {
        std::cerr << "FAIL parser-parity: " << error.what() << '\n';
        return 1;
    }
    return 0;
}

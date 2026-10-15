/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/cip.hpp"
#include "smsd/smiles_parser.hpp"
#include <iostream>
#include <stdexcept>

int main() {
    struct Example { const char* smiles; int centre; smsd::cip::RSLabel expected; };
    using smsd::cip::RSLabel;
    const Example examples[] = {
        {"F[C@H](C(OC)C)C(O)CO", 1, RSLabel::S},
        {"[2H]C([2H])([C@H](O)CCF)C", 3, RSLabel::S},
        {"O[C@H](C[13CH2]C)CCCF", 1, RSLabel::R},
        {"F[C@H](C(CO)C(CBr)(CBr)CBr)C(C(O)C)C(C)C", 1, RSLabel::R},
        {"O[C@@H](c1ccccc1)C(C)(C)C", 1, RSLabel::R},
        {"O[C@@H](C1=CC=CC=C1)C(C)(C)C", 1, RSLabel::R},
        {"O=C1Nc2ccc(Cl)cc2[C@@](C#CC2CC2)(C(F)(F)F)O1", 10, RSLabel::S},
        {"C[S@](=O)c1ccccc1", 1, RSLabel::S},
        {"COc1ccc2[nH]c([S@@](=O)Cc3ncc(C)c(OC)c3C)nc2c1", 8, RSLabel::S},
        {"C/C=C/[C@@H](F)/C=C\\C", 3, RSLabel::R},
        {"C/C([C@H](F)/C(C)=C/C)=C/C", 2, RSLabel::R},
        {"C/C=C(C)/[C@@H](F)/C(C)=C\\C", 4, RSLabel::R},
        {"C/C(Cl)=C(C)/[C@@H](F)/C(C)=C(Cl)\\C", 5, RSLabel::S},
        {"N[C@@H](C)C(=O)O", 1, RSLabel::S},
        {"N[C@H](C)C(=O)O", 1, RSLabel::R}
    };
    int failures = 0;
    for (const auto& example : examples) {
        const auto graph = smsd::parseSMILES(example.smiles);
        const auto result = smsd::cip::assignRS(graph, example.centre);
        const auto all = smsd::cip::assignAll(graph);
        if (result != example.expected || all.rsLabels[example.centre] != example.expected) {
            ++failures;
            std::cout << "FAIL " << example.smiles << " expected " << smsd::cip::to_string(example.expected)
                      << " got " << smsd::cip::to_string(result) << '\n';
        } else {
            std::cout << "PASS " << example.smiles << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}

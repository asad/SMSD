/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/rdkit_adapter.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/cip.hpp"
#include <stdexcept>

int main() {
    for (const char* text : {"N[C@H](C)C(=O)O", "N[C@@H](C)C(=O)O",
                             "[C@H](F)(Cl)Br", "[C@@H](F)(Cl)Br",
                             "F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I"}) {
        const auto native = smsd::parseSMILES(text);
        const auto adapted = smsd::fromSmiles(text);
        for (int atom = 0; atom < native.n; ++atom)
            if (native.tetraChirality[atom]
                && smsd::cip::assignRS(native, atom) != smsd::cip::assignRS(adapted, atom))
                throw std::runtime_error("Adapter changed tetrahedral stereochemistry");
    }
}

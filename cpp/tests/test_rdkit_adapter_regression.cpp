/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/rdkit_adapter.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/cip.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

struct Case {
    const char* text;
    /** The same molecule as the native parser reads it once explicit hydrogens are removed. */
    const char* reference;
};

/**
 * The adapter puts the implicit ligand of a three-neighbour centre where the graph expects it,
 * so a centre whose first neighbour follows it (the first atom, or one written before its
 * neighbours) has its tag turned, and one whose first neighbour precedes it keeps it. The CIP
 * labels, and the SMILES written from the adapted graph, must agree with the native parser for
 * both kinds of centre, with and without explicit hydrogens, and in rings.
 */
void tetrahedralCentresRoundTrip() {
    for (const Case& item : {Case{"N[C@H](C)C(=O)O", "N[C@H](C)C(=O)O"}, Case{"N[C@@H](C)C(=O)O", "N[C@@H](C)C(=O)O"},
                             Case{"[C@H](F)(Cl)Br", "[C@H](F)(Cl)Br"}, Case{"[C@@H](F)(Cl)Br", "[C@@H](F)(Cl)Br"},
                             Case{"F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I"}, Case{"F[C@@](Cl)(Br)I", "F[C@@](Cl)(Br)I"},
                             Case{"[H][C@](F)(Cl)Br", "[C@H](F)(Cl)Br"}, Case{"[H][C@@](F)(Cl)Br", "[C@@H](F)(Cl)Br"},
                             Case{"[C@H]1(F)CCCO1", "[C@H]1(F)CCCO1"}, Case{"[C@@H]1(F)CCCO1", "[C@@H]1(F)CCCO1"},
                             Case{"C[C@H]1CCCCO1", "C[C@H]1CCCCO1"}, Case{"C[C@@H]1CCCCO1", "C[C@@H]1CCCCO1"},
                             Case{"O[C@H](C)[C@@H](N)C", "O[C@H](C)[C@@H](N)C"}}) {
        const auto native = smsd::parseSMILES(item.reference);
        const auto adapted = smsd::fromSmiles(item.text);
        require(native.n == adapted.n, std::string("Adapter changed the heavy atom count of ") + item.text);
        for (int atom = 0; atom < native.n; ++atom)
            require(smsd::cip::assignRS(native, atom) == smsd::cip::assignRS(adapted, atom),
                    std::string("Adapter changed the tetrahedral stereochemistry of ") + item.text);
        const auto written = smsd::parseSMILES(smsd::toSMILES(adapted));
        require(written.toCanonicalSmiles() == native.toCanonicalSmiles(),
                std::string("Adapted graph did not round trip through SMILES: ") + item.text + " gave " + written.toCanonicalSmiles());
    }
}

} // namespace

int main() {
    try {
        tetrahedralCentresRoundTrip();
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "PASS RDKit adapter tetrahedral round trip\n";
    return 0;
}

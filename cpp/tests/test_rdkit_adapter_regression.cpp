/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/rdkit_adapter.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/cip.hpp"
#include <GraphMol/Chirality.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <memory>
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

/**
 * The adapter and SMSD's own SMILES parser must build the same graph for the same molecule: the
 * same hydrogen counts (aromatic N-H, charged and isotopic atoms included), the same aromatic bond
 * orders, the same ring and aromatic flags, and the same double-bond configurations. The canonical
 * SMILES of the two graphs, which carries all of these, must also agree.
 */
void atomsBondsAndDoubleBondsAgreeWithParser() {
    for (const char* text : {"c1cc[nH]c1", "c1ccc2[nH]ccc2c1", "c1c[nH]cn1", "c1nc[nH]n1", "c1nn[nH]n1",
                             "c1cc[nH+]cc1", "C[n+]1ccccc1", "c1nnn[n-]1", "O=c1cc[nH]cc1", "Oc1ccccc1",
                             "Nc1ccccc1", "c1ccc2ccccc2c1", "[NH4+]", "[13CH4]", "CC(=O)[O-]",
                             "O=[N+]([O-])c1ccccc1", "F/C=C/F", "F/C=C\\F", "C/C=C/C(=O)O", "C/C(F)=C/C",
                             "C/C(F)=C\\C", "Cl/C(F)=C(/Br)I", "F/C=C/C=C/F", "C/C=C1\\CCOC1", "OC(=O)/C=C\\C(=O)O"}) {
        const auto native = smsd::parseSMILES(text);
        const auto adapted = smsd::fromSmiles(text);
        const std::string where = std::string(" of ") + text;
        require(native.n == adapted.n, "Adapter changed the heavy atom count" + where);
        require(adapted.hydrogenCountsKnown, "Adapter left the hydrogen counts unknown" + where);
        for (int atom = 0; atom < native.n; ++atom) {
            require(native.atomicNum[atom] == adapted.atomicNum[atom], "Adapter changed an element" + where);
            require(native.formalCharge[atom] == adapted.formalCharge[atom], "Adapter changed a charge" + where);
            require(native.massNumber[atom] == adapted.massNumber[atom], "Adapter changed an isotope" + where);
            require(native.aromatic[atom] == adapted.aromatic[atom], "Adapter changed an aromatic flag" + where);
            require(native.ring[atom] == adapted.ring[atom], "Adapter changed a ring flag" + where);
            require(native.implicitHydrogenCount(atom) == adapted.implicitHydrogenCount(atom),
                    "Adapter changed the hydrogen count of atom " + std::to_string(atom) + where);
            for (int other = atom + 1; other < native.n; ++other) {
                require(native.bondOrder(atom, other) == adapted.bondOrder(atom, other),
                        "Adapter changed a bond order" + where);
                require(native.bondAromatic(atom, other) == adapted.bondAromatic(atom, other),
                        "Adapter changed an aromatic bond flag" + where);
            }
        }
        require(native.toCanonicalSmiles() == adapted.toCanonicalSmiles(),
                "Adapted graph differs from the parsed one" + where + ": " + adapted.toCanonicalSmiles()
                    + " against " + native.toCanonicalSmiles());
    }
}

/** Whether neighbour x of a and neighbour y of b lie on the same side of the double bond a=b. */
bool sameSide(const smsd::MolGraph& g, int a, int b, int x, int y) {
    bool together = g.dbStereo(a, b) == 1;
    for (int neighbor : g.neighbors[a]) {
        if (neighbor == b) continue;
        if (neighbor != x) together = !together;
        break;
    }
    for (int neighbor : g.neighbors[b]) {
        if (neighbor == a) continue;
        if (neighbor != y) together = !together;
        break;
    }
    return together;
}

/** The adapted graph must have the canonical SMILES of the parsed reference. */
void requireSameAs(const smsd::MolGraph& adapted, const std::string& reference, const std::string& where) {
    const auto native = smsd::parseSMILES(reference);
    require(native.toCanonicalSmiles() == adapted.toCanonicalSmiles(),
            "Adapted graph differs from the parsed one" + where + ": " + adapted.toCanonicalSmiles() + " against "
                + native.toCanonicalSmiles());
}

/**
 * When hydrogens are removed, a hydrogen that is a stereo atom of a double bond is replaced by its
 * partner, and the configuration must survive. Isotopic hydrogens, an eight-membered ring, a
 * graph large enough for sparse storage, and molecules read under RDKit's newer stereo perception
 * must also agree with the parser.
 */
void stereoAtomsIsotopesRingsAndLargeGraphs() {
    // [H]/C(*)=C/F: the dummy atom and F are cis once the hydrogen is removed. After removal the
    // atoms are C 0, dummy 1, C 2 and F 3, with the double bond between 0 and 2.
    for (const char* text : {"[H]/C(*)=C/F", "[H]/C([1*])=C/F"}) {
        RDKit::SmilesParserParams keepHydrogens;
        keepHydrogens.removeHs = false;
        std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(text, keepHydrogens));
        require(mol != nullptr, std::string("RDKit did not parse ") + text);
        const auto adapted = smsd::fromRDKit(*mol);
        require(adapted.n == 4 && adapted.dbStereo(0, 2) != 0, std::string("Adapter lost the double-bond configuration of ") + text);
        require(sameSide(adapted, 0, 2, 1, 3), std::string("Adapter turned the double-bond configuration of ") + text);
    }
    const std::string longChain = "F/C=C/" + std::string(230, 'C');
    for (const std::string& text : {std::string("[2H]C([2H])([2H])O"), std::string("[2H]/C=C/F"),
                                    std::string("C1CCC/C=C/CC1"), longChain})
        requireSameAs(smsd::fromSmiles(text), text, " of " + text.substr(0, 40));
    const bool legacy = RDKit::Chirality::getUseLegacyStereoPerception();
    RDKit::Chirality::setUseLegacyStereoPerception(false);
    try {
        for (const char* text : {"F/C=C/F", "C/C(F)=C\\C", "OC(=O)/C=C\\C(=O)O", "Cl/C(F)=C(/Br)I"})
            requireSameAs(smsd::fromSmiles(text), text, std::string(" under the newer stereo perception of ") + text);
    } catch (...) {
        RDKit::Chirality::setUseLegacyStereoPerception(legacy);
        throw;
    }
    RDKit::Chirality::setUseLegacyStereoPerception(legacy);
}

/** With hydrogens kept as atoms, each heavy atom has its hydrogens as neighbours and counts none twice. */
void explicitHydrogensKept() {
    for (const char* text : {"c1cc[nH]c1", "CC(=O)O", "c1nnn[n-]1"}) {
        std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(text));
        require(mol != nullptr, std::string("RDKit did not parse ") + text);
        RDKit::MolOps::addHs(*mol);
        const auto native = smsd::parseSMILES(text);
        const auto adapted = smsd::fromRDKit(*mol, false);
        require(adapted.n == static_cast<int>(mol->getNumAtoms()), std::string("Adapter dropped hydrogen atoms of ") + text);
        for (int atom = 0; atom < native.n; ++atom) {
            require(adapted.implicitHydrogenCount(atom) == 0, std::string("Adapter counted kept hydrogens twice in ") + text);
            require(adapted.attachedHydrogenCount(atom) == native.implicitHydrogenCount(atom),
                    std::string("Adapter changed the attached hydrogens of ") + text);
        }
    }
}

/** An unspecified double bond has no configuration, and a query molecule keeps unknown hydrogen counts. */
void unspecifiedStereoAndQueries() {
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol("FC=CF"));
    require(mol != nullptr, "RDKit did not parse FC=CF");
    auto* bond = mol->getBondBetweenAtoms(1, 2);
    bond->setStereoAtoms(0, 3);
    bond->setStereo(RDKit::Bond::STEREOANY);
    require(smsd::fromRDKit(*mol).dbStereo(1, 2) == 0, "Adapter gave an unspecified double bond a configuration");
    std::unique_ptr<RDKit::RWMol> query(RDKit::SmartsToMol("[NH2]C"));
    require(query != nullptr, "RDKit did not parse the SMARTS [NH2]C");
    require(!smsd::fromRDKit(*query).hydrogenCountsKnown, "Adapter gave a query molecule known hydrogen counts");
}

} // namespace

int main() {
    try {
        tetrahedralCentresRoundTrip();
        atomsBondsAndDoubleBondsAgreeWithParser();
        stereoAtomsIsotopesRingsAndLargeGraphs();
        explicitHydrogensKept();
        unspecifiedStereoAndQueries();
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "PASS RDKit adapter tetrahedral round trip and parser parity\n";
    return 0;
}

#include "search_test_helpers.hpp"
/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/smsd.hpp"
#include "smsd/smarts_parser.hpp"
#include "smsd/smiles_parser.hpp"
#include <functional>
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void checkMolecule(const char* queryText, const char* targetText, bool expected,
                          bool chirality, bool bondStereo) {
    const auto query = smsd::parseSMILES(queryText), target = smsd::parseSMILES(targetText);
    smsd::ChemOptions chemistry;
    chemistry.useChirality = chirality;
    chemistry.useBondStereo = bondStereo;
    {

        require(smsd_test::substructureExists(query, target, chemistry) == expected, "molecular stereo existence is wrong");
        require(!smsd_test::substructureAtoms(query, target, chemistry).empty() == expected, "molecular stereo mapping is wrong");
        require(!smsd_test::allSubstructureAtoms(query, target, chemistry).empty() == expected, "molecular stereo enumeration is wrong");
    }
}

static void checkSmarts(const char* queryText, const char* targetText, bool expected) {
    const auto query = smsd::parseSMARTS(queryText);
    const auto target = smsd::parseSMILES(targetText);
    require(smsd_test::patternMatches(query, target) == expected, "SMARTS stereo existence is wrong");
    require(!smsd_test::allPatternAtoms(query, target).empty() == expected, "SMARTS stereo enumeration is wrong");
}

int main() {
    int failures = 0;
    auto run = [&](const char* name, const std::function<void()>& test) {
        try { test(); std::cout << "PASS " << name << '\n'; }
        catch (const std::exception& error) { ++failures; std::cout << "FAIL " << name << ": " << error.what() << '\n'; }
    };
    run("mapped chirality with changed CIP priorities", [] {
        checkMolecule("C[C@H](O)CC", "ClC[C@H](O)CC", true, true, false);
        checkMolecule("C[C@@H](O)CC", "ClC[C@H](O)CC", false, true, false);
        checkMolecule("N[C@@H](C)C(=O)O", "N[C@@H](CS)C(=O)O", true, true, false);
    });
    run("implicit hydrogen maps to added substituent", [] {
        checkMolecule("N[C@@H](C)C(=O)O", "N[C@@](F)(C)C(=O)O", true, true, false);
        checkMolecule("N[C@@H](C)C(=O)O", "N[C@](F)(C)C(=O)O", false, true, false);
    });
    run("mapped double-bond reference ligands", [] {
        checkMolecule("C/C=C/C", "Cl/C(C)=C/C", false, false, true);
        checkMolecule("C/C=C\\C", "Cl/C(C)=C/C", true, false, true);
        checkMolecule("C/C(Cl)=C/C", "Cl/C(C)=C\\C", true, false, true);
        checkMolecule("Cl/C(C)=C\\C", "C/C(Cl)=C/C", true, false, true);
    });
    run("linear path double-bond stereo", [] {
        checkMolecule("F/C=C/Cl", "F/C=C\\Cl", false, false, true);
        checkMolecule("F/C=C/Cl", "Cl/C=C/F", true, false, true);
    });
    run("full matchers and sparse target stereo", [] {
        const std::string padding = "." + std::string(520, 'C');
        const std::string sameChirality = "ClC[C@H](O)CC" + padding;
        checkMolecule("C[C@H](O)CC", sameChirality.c_str(), true, true, false);
        checkMolecule("C[C@@H](O)CC", sameChirality.c_str(), false, true, false);
        const std::string cis = "Cl/C(C)=C/C" + padding;
        checkMolecule("C/C=C/C", cis.c_str(), false, false, true);
        checkMolecule("C/C=C\\C", cis.c_str(), true, false, true);
    });
    run("MCS mapped stereochemistry", [] {
        const auto query = smsd::parseSMILES("C[C@H](O)CC"), target = smsd::parseSMILES("ClC[C@H](O)CC");
        smsd::ChemOptions chemistry;
        chemistry.useChirality = true;
        require(smsd_test::commonAtoms(query, target, chemistry, {}).size() == static_cast<size_t>(query.n), "MCS must preserve mapped winding despite changed CIP");
        chemistry.useChirality = false;
        chemistry.useBondStereo = true;
        const auto trans = smsd::parseSMILES("C/C=C/C"), cis = smsd::parseSMILES("Cl/C(C)=C/C");
        require(smsd_test::commonAtoms(trans, cis, chemistry, {}).size() < 4, "MCS must reject the incompatible complete alkene");
        require(smsd_test::commonAtoms(smsd::parseSMILES("C/C=C\\C"), cis, chemistry, {}).size() == 4, "MCS must retain compatible mapped alkene references");
    });
    run("SMARTS mapped tetrahedral winding", [] {
        checkSmarts("[C@@H](F)(Cl)Br", "[C@H](Cl)(F)Br", true);
        checkSmarts("[C@@H](F)(Cl)Br", "[C@@H](Cl)(F)Br", false);
        checkSmarts("C[C@H]1CCCCO1", "C[C@@H]1OCCCC1", true);
        const auto query = smsd::parseSMARTS("[C@@H](F)(Cl)Br");
        require(smsd_test::anchoredPatternMatches(query, smsd::parseSMILES("[C@H](Cl)(F)Br"), 0), "anchored stereo matching must use mapped winding");
    });
    run("SMARTS Boolean chirality", [] {
        checkSmarts("[C@@H,N](F)(Cl)Br", "[C@H](Cl)(F)Br", true);
        checkSmarts("[C@@H,N](F)(Cl)Br", "N(F)(Cl)Br", true);
        checkSmarts("[C;!@@;H1](F)(Cl)Br", "[C@H](F)(Cl)Br", true);
        checkSmarts("[C;!@@;H1](F)(Cl)Br", "[C@@H](F)(Cl)Br", false);
        checkSmarts("[C;!@@;H1](F)(Cl)Br", "[C@@H](Cl)(F)Br", true);
        checkSmarts("[C;!@@;H1](F)(Cl)Br", "[C@H](Cl)(F)Br", false);
        checkSmarts("[C@?H](F)(Cl)Br", "C(F)(Cl)Br", true);
    });
    run("SMARTS double-bond stereo", [] {
        checkSmarts("F/C=C/F", "F/C=C/F", true);
        checkSmarts("F/C=C/F", "F/C=C\\F", false);
        checkSmarts("F/C=C/F", "FC=CF", false);
        checkSmarts("C/C", "CC", true);
        checkSmarts("C/C", "C=C", false);
        checkSmarts("C/C=C/C", "Cl/C(C)=C/C", false);
        checkSmarts("C/C=C\\C", "Cl/C(C)=C/C", true);
    });
    run("SMARTS Boolean directional bonds", [] {
        checkSmarts("F/,-C=C/F", "F/C=C/F", true);
        checkSmarts("F/,-C=C/F", "F/C=C\\F", true);
        checkSmarts("F!/C=C/F", "F/C=C\\F", true);
        checkSmarts("F!/C=C/F", "F/C=C/F", false);
    });
    run("SMARTS unspecified stereo and Boolean alternatives", [] {
        checkSmarts("[C@@H](F)(Cl)Br", "C(F)(Cl)Br", false);
        checkSmarts("[C;!@@;H1](F)(Cl)Br", "C(F)(Cl)Br", true);
        checkSmarts("[C@@H,N](F)(Cl)Br", "N(F)(Cl)Br", true);
        checkSmarts("F/,-C=C/F", "FC=CF", true);
        checkSmarts("F!/C=C/F", "FC=CF", false);
        checkMolecule("C[C@H](O)CC", "CC(O)CC", true, true, false);
    });
    run("SMARTS conjugated directional carriers", [] {
        checkSmarts("C/C=C(/F)\\C(Cl)=C/C", "C\\C=C(\\F)/C(Cl)=C\\C", true);
        checkSmarts("C/C=C(/F)\\C(Cl)=C/C", "C/C=C(\\F)/C(Cl)=C/C", false);
    });
    run("recursive SMARTS stereo", [] {
        checkSmarts("[$([C@@H](F)(Cl)Br)]", "[C@H](Cl)(F)Br", true);
        checkSmarts("[$([C@@H](F)(Cl)Br)]", "[C@@H](Cl)(F)Br", false);
    });
    return failures == 0 ? 0 : 1;
}

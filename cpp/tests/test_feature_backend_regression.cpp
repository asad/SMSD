/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/batch.hpp"
#include "smsd/gpu.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

using smsd_test::require;

static void pharmacophores() {
    for (const auto& item : {std::pair<const char*, int>{"CC(=O)C", 2}, {"c1ccncc1", 3}, {"CC#N", 2}}) {
        const auto molecule = smsd::parseSMILES(item.first);
        const int features = smsd::batch::detail::classifyPharmacophore(molecule, item.second);
        require((features & 1) == 0, "An atom without attached hydrogen cannot be a donor");
        require((features & 2) != 0, "Carbonyl, pyridine and nitrile acceptors must be retained");
        require(features == molecule.getPharmacophoreFeatures()[item.second], "Cached and batch feature masks differ");
    }
    require((smsd::batch::detail::classifyPharmacophore(smsd::parseSMILES("CC#N"), 2) & 4) == 0,
            "Nitrile nitrogen is not an ionisable amine");
    require((smsd::batch::detail::classifyPharmacophore(smsd::parseSMILES("C=N"), 1) & 4) == 0,
            "An imine must not be classified as a saturated basic amine");
    for (const auto& item : {std::pair<const char*, int>{"[NH4+]", 0}, {"C[NH3+]", 1}, {"CO", 1}, {"[H]OC", 1}, {"c1cc[nH]c1", 3}})
        require((smsd::batch::detail::classifyPharmacophore(smsd::parseSMILES(item.first), item.second) & 1) != 0,
                "Attached hydrogen donors must include explicit H vertices and cations");
    require(smsd::batch::detail::computeCircularFingerprintFCFP(smsd::parseSMILES("CC(=O)C"), 0, 2048)
            != smsd::batch::detail::computeCircularFingerprintFCFP(smsd::parseSMILES("CC(O)C"), 0, 2048),
            "Ketone and hydroxyl feature fingerprints must differ");
}

static void attachedHydrogens() {
    for (const auto& pair : {std::pair<const char*, const char*>{"[NH4+]", "[NH2+]"}, {"C[NH3+]", "C[N+]"}, {"[CH2]", "C"}}) {
        const auto first = smsd::parseSMILES(pair.first), second = smsd::parseSMILES(pair.second);
        require(smsd::batch::detail::computeCircularFingerprintECFP(first, 2, 2048)
                != smsd::batch::detail::computeCircularFingerprintECFP(second, 2, 2048), "ECFP must retain stored H counts");
        require(smsd::batch::detail::computeCircularFingerprintECFPCounts(first, 2, 2048)
                != smsd::batch::detail::computeCircularFingerprintECFPCounts(second, 2, 2048), "Count ECFP must retain stored H counts");
    }
    require(smsd::detail::computeImplicitH(7, false, 0, 1) == 4, "Ammonium fallback must have four hydrogens");
    require(smsd::detail::computeImplicitH(5, false, 0, -1) == 4, "Borohydride fallback must have four hydrogens");
}

static void zeroAtomScreen() {
    const auto query = smsd::parseSMILES("CCO");
    const std::vector<smsd::MolGraph> targets{smsd::parseSMILES("CCN"), {}, smsd::parseSMILES("c1ccccc1")};
    const auto hits = smsd::gpu::batchSimilarityScreenAuto(query, targets, 0.0);
    require(hits.size() == 3 && hits[1].index == 1 && hits[1].score == 0.0, "Zero threshold must include zero-atom targets");
    const auto positive = smsd::gpu::batchSimilarityScreenAuto(query, targets, 0.01);
    for (const auto& hit : positive) require(hit.index != 1, "Positive thresholds must exclude zero-overlap targets");
}

static void backendVersions() {
    require(smsd::gpu::detail::openMPVersion(201511) == "4.5", "OpenMP 4.5 date is not its version number");
    require(smsd::gpu::detail::openMPVersion(202011) == "5.1", "OpenMP 5.1 date is not its version number");
    require(smsd::gpu::detail::openMPVersion(999999) == "999999", "Unknown OpenMP dates must remain visible");
#ifdef _OPENMP
    require(smsd::gpu::deviceInfo().find("OpenMP " + smsd::gpu::detail::openMPVersion(_OPENMP)) != std::string::npos,
            "Device information must display the specification version");
#endif
}

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"pharmacophores", pharmacophores}, {"hydrogens", attachedHydrogens},
        {"zero-atoms", zeroAtomScreen}, {"versions", backendVersions}});
}

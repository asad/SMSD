/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "test_support.hpp"
#include "smsd/mcs.hpp"

#include <cstdint>

int main() {
    const auto ordinary = smsd::parseSMILES("CC");
    smsd_test::require(ordinary.useDense, "The small-molecule fixture must use dense bond storage");
    smsd_test::require(smsd::common_bond_detail::CommonBondGraphs::supports(ordinary, ordinary),
                      "Ordinary reciprocal bond properties must remain supported");
    smsd::MCSOptions options;
    options.objective = smsd::MCSObjective::MCES;
    options.timeoutMs = 30000;
    const auto control = smsd::MCS::find(ordinary).withOptions(options).match(ordinary);
    smsd_test::require(control.optimal && control.complete && control.score == 1,
                      "Ordinary common-edge search must retain its proved result");

    for (uint8_t reverse : {uint8_t{2}, uint8_t{1 | smsd::MolGraph::BOND_AROM_BIT},
                            uint8_t{1 | smsd::MolGraph::BOND_RING_BIT}}) {
        auto asymmetric = ordinary;
        asymmetric.bondPacked[1][0] = reverse;
        smsd_test::require(!smsd::common_bond_detail::CommonBondGraphs::supports(ordinary, asymmetric)
                          && !smsd::common_bond_detail::CommonBondGraphs::supports(asymmetric, ordinary),
                          "Asymmetric bond properties must be declined in either orientation");
        const auto atoms = smsd::common_atom_detail::searchCommonAtoms(ordinary, asymmetric,
            [](int, int) { return true; }, [](int, int, int, int) { return true; },
            [](const smsd::MolGraph&, int, int) { return int64_t{0}; }, {}, 1, [] { return false; });
        smsd_test::require(!atoms.supported && !atoms.finished && atoms.mapping.empty(),
                          "The common-atom engine must share unsupported-graph admission");
        const auto result = smsd::MCS::find(ordinary).withOptions(options).match(asymmetric);
        smsd_test::require(result.status == smsd::SearchStatus::LimitReached
                          && !result.optimal && !result.complete && result.matches.empty(),
                          "Declined common-edge graphs must not be reported as a proved negative");
    }
    return 0;
}

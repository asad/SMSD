/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#define STRICT 1701
#define IGNORE 1702
#endif

constexpr int strictMacroBeforeHeaders = STRICT;
constexpr int ignoreMacroBeforeHeaders = IGNORE;

#include "smsd/chem_options.hpp"
#include "smsd/detail/graph_operations.hpp"
#include "smsd/query_smarts.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/detail/search_compatibility.hpp"
#include "smsd/search.hpp"
#include "smsd/substructure.hpp"
#include "smsd/smsd.hpp"

static_assert(STRICT == strictMacroBeforeHeaders,
              "SMSD headers changed the caller's STRICT macro");
static_assert(IGNORE == ignoreMacroBeforeHeaders,
              "SMSD headers changed the caller's IGNORE macro");
#pragma push_macro("STRICT")
#pragma push_macro("IGNORE")
#undef STRICT
#undef IGNORE

int main() {
    smsd::ChemOptions options;
    options.matchBondOrder = smsd::ChemOptions::BondOrderMode::STRICT;
    options.aromaticityMode = smsd::ChemOptions::AromaticityMode::STRICT;
    options.ringFusionMode = smsd::ChemOptions::RingFusionMode::IGNORE;
    (void)options;
    return 0;
}

#pragma pop_macro("IGNORE")
#pragma pop_macro("STRICT")

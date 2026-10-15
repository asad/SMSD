/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/batch.hpp"
#include "smsd/smiles_parser.hpp"
#include <stdexcept>

int main() {
    for (const char* name : {"strict", "STRICT", "compat-substruct", "Compat-FMCS"}) {
        const auto options = smsd::ChemOptions::profile(name);
        if (options.matchFormalCharge || options.ringMatchesRingOnly)
            throw std::runtime_error("Profile must not enable unrelated matching constraints");
    }
    bool rejected = false;
    try { smsd::ChemOptions::profile("strcit"); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) throw std::runtime_error("Unknown profile must be rejected");
    auto fingerprint = &smsd::batch::detail::computeMCSFingerprint;
    const auto words = fingerprint(smsd::parseSMILES("C"), 7, 2048);
    if (words[30] != (1ULL << 53))
        throw std::runtime_error("Default MCS fingerprint atom labels changed");
}

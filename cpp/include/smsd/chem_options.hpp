/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * ChemOptions: the chemistry rules of matching, with named profiles.
 */
#pragma once
#ifndef SMSD_CHEM_OPTIONS_HPP
#define SMSD_CHEM_OPTIONS_HPP

#include "smsd/parser_utils.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace smsd {

enum class AromaticityModel { DAYLIGHT_LIKE };

// ============================================================================
// ChemOptions -- configuration for chemical matching constraints
// ============================================================================

struct ChemOptions {
    enum class BondOrderMode   { STRICT, LOOSE, ANY };
    enum class AromaticityMode { STRICT, FLEXIBLE };
    enum class RingFusionMode  { IGNORE, PERMISSIVE, STRICT };
    enum class Solvent         { AQUEOUS, DMSO, METHANOL, CHLOROFORM, ACETONITRILE, DIETHYL_ETHER };

    // Atom matching: the defaults follow RDKit FindMCS for like-for-like comparison
    bool matchAtomType            = true;
    bool matchFormalCharge        = false;  // RDKit MCSAtomCompareParameters.MatchFormalCharge=False
    bool useChirality             = false;
    bool useBondStereo            = false;
    bool ringMatchesRingOnly      = false;  // RDKit ringMatchesRingOnly=False

    // Ring completeness
    bool completeRingsOnly        = false;

    // Isotope matching
    bool matchIsotope             = false;

    // Ring fusion
    RingFusionMode ringFusionMode = RingFusionMode::IGNORE;

    // Tautomer matching
    bool   tautomerAware = false;
    /** pH context for MCS::computeTautomerConfidence; does not change the search objective. */
    double pH            = 7.4;

    /** Set simulation pH, clamped to [0, 14]. */
    ChemOptions& withPH(double v) {
        if (!std::isfinite(v)) throw std::invalid_argument("pH must be finite");
        pH = std::max(0.0, std::min(14.0, v)); return *this;
    }

    /** Solvent context for MCS::computeTautomerConfidence. */
    Solvent solvent = Solvent::AQUEOUS;

    /** Set solvent for tautomer weight corrections (fluent). */
    ChemOptions& withSolvent(Solvent s) { solvent = s; return *this; }

    // Bond matching
    BondOrderMode  matchBondOrder  = BondOrderMode::STRICT;
    AromaticityMode aromaticityMode = AromaticityMode::FLEXIBLE;

    /**
     * @deprecated No search reads this field. Aromaticity is perceived when a graph is built, by
     * MolGraph::Builder::build() or perceiveAromaticity(), which take the model directly. It is
     * not marked with an attribute, because compilers would then warn on every copy of ChemOptions.
     */
    AromaticityModel aromaticityModel = AromaticityModel::DAYLIGHT_LIKE;

    // Fluent setters
    ChemOptions& withAromaticityMode(AromaticityMode m) { aromaticityMode = m; return *this; }
    ChemOptions& withBondOrderMode(BondOrderMode m)     { matchBondOrder = m; return *this; }
    ChemOptions& withCompleteRingsOnly(bool on)          { completeRingsOnly = on; return *this; }
    ChemOptions& withMatchIsotope(bool on)               { matchIsotope = on; return *this; }
    ChemOptions& withRingFusionMode(RingFusionMode m)    { ringFusionMode = m; return *this; }
    ChemOptions& withTautomerAware(bool on)              { tautomerAware = on; return *this; }

    /**
     * A named profile, ignoring the case of the name: "default", "strict", "pharma",
     * "compat-fmcs", "compat-substruct" or "reaction". Each starts from the defaults and
     * changes only the fields set below; an unknown name throws std::invalid_argument.
     */
    static ChemOptions profile(const std::string& profileName) {
        std::string name = profileName;
        std::transform(name.begin(), name.end(), name.begin(), parser_detail::asciiLower);
        ChemOptions opts;
        if (name == "default") {
            // The constructor defaults already follow RDKit FindMCS
        } else if (name == "strict") {
            // Strict aromaticity as well as the default strict bond orders
            opts.aromaticityMode     = AromaticityMode::STRICT;
        } else if (name == "pharma") {
            // Drug discovery screening: charge, ring-only and complete-ring matching
            opts.matchFormalCharge   = true;
            opts.ringMatchesRingOnly = true;
            opts.completeRingsOnly   = true;
        } else if (name == "compat-fmcs" || name == "compat-substruct" || name == "reaction") {
            // Loose bond orders; otherwise the defaults
            opts.matchBondOrder      = BondOrderMode::LOOSE;
        } else {
            throw std::invalid_argument("Unknown ChemOptions profile: " + profileName);
        }
        return opts;
    }

    /** Tautomer-aware matching with loose bond orders; otherwise the defaults. */
    static ChemOptions tautomerProfile() {
        ChemOptions c;
        c.tautomerAware  = true;
        c.matchBondOrder = BondOrderMode::LOOSE;
        return c;
    }

    /** @deprecated The same options as profile("compat-fmcs"). */
    [[deprecated("Use ChemOptions::profile(\"compat-fmcs\"), which returns the same options.")]]
    static ChemOptions fmcsProfile() { return profile("compat-fmcs"); }
};

} // namespace smsd

#endif // SMSD_CHEM_OPTIONS_HPP

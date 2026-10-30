/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <string>

namespace RDKit { class ROMol; }

namespace smsd {

/// Convert a molecule using the optional smsd_rdkit library.
/// Hydrogen removal changes atom indices. Each atom keeps its hydrogen count, and a double
/// bond keeps its configuration when it is given against stereo atoms. A query molecule, or one
/// whose hydrogen counts RDKit cannot compute, keeps unknown hydrogen counts. A double bond given
/// as E or Z without stereo atoms, an unspecified or atropisomeric bond, and a double bond in a
/// ring of fewer than eight atoms (whose configuration RDKit clears) carry no configuration.
MolGraph fromRDKit(const RDKit::ROMol& mol, bool removeHs = true);

/// Parse SMILES with RDKit and remove explicit hydrogens.
MolGraph fromSmiles(const std::string& smiles);

} // namespace smsd

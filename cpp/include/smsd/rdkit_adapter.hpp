/* SPDX-License-Identifier: Apache-2.0 */
#pragma once

#include "smsd/mol_graph.hpp"
#include <string>

namespace RDKit { class ROMol; }

namespace smsd {

/// Convert a molecule using the optional smsd_rdkit library.
/// Hydrogen removal changes atom indices; hydrogen counts and bond stereo
/// are not imported.
MolGraph fromRDKit(const RDKit::ROMol& mol, bool removeHs = true);

/// Parse SMILES with RDKit and remove explicit hydrogens.
MolGraph fromSmiles(const std::string& smiles);

} // namespace smsd

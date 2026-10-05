/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * RDKit adapter: converts RDKit ROMol to SMSD MolGraph.
 * This is the ONLY file that depends on RDKit.
 * The core algorithms (mol_graph.hpp, vf2pp.hpp, mcs.hpp) are RDKit-free.
 */

#ifdef SMSD_WITH_RDKIT

#include "smsd/mol_graph.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <memory>
#include <stdexcept>

namespace smsd {

/// Convert an RDKit ROMol to an SMSD MolGraph.
/// Strips explicit hydrogens by default.
MolGraph fromRDKit(const RDKit::ROMol& mol, bool removeHs = true) {
    // Optionally remove Hs
    std::unique_ptr<RDKit::RWMol> molCopy;
    const RDKit::ROMol* m = &mol;
    if (removeHs) {
        molCopy.reset(new RDKit::RWMol(mol));
        RDKit::MolOps::removeHs(*molCopy);
        m = molCopy.get();
    }

    int n = m->getNumAtoms();

    // Collect atom properties
    std::vector<int> atomicNum(n), formalCharge(n), massNumber(n, 0);
    std::vector<uint8_t> ringFlag(n, 0), aromFlag(n, 0);

    RDKit::MolOps::findSSSR(*const_cast<RDKit::ROMol*>(m));
    const auto* ri = m->getRingInfo();

    for (int i = 0; i < n; ++i) {
        const auto* atom = m->getAtomWithIdx(i);
        atomicNum[i]    = atom->getAtomicNum();
        formalCharge[i] = atom->getFormalCharge();
        massNumber[i]   = atom->getIsotope();
        aromFlag[i]     = atom->getIsAromatic() ? 1 : 0;
        ringFlag[i]     = ri->numAtomRings(i) > 0 ? 1 : 0;
    }

    // Build adjacency + per-neighbor bond properties
    std::vector<std::vector<int>> neighbors(n);
    std::vector<std::vector<int>> bondOrders(n);
    std::vector<std::vector<bool>> bondRingFlags(n);
    std::vector<std::vector<bool>> bondAromFlags(n);

    for (const auto& bond : m->bonds()) {
        int a = bond->getBeginAtomIdx();
        int b = bond->getEndAtomIdx();

        int bo;
        switch (bond->getBondType()) {
            case RDKit::Bond::SINGLE:   bo = 1; break;
            case RDKit::Bond::DOUBLE:   bo = 2; break;
            case RDKit::Bond::TRIPLE:   bo = 3; break;
            case RDKit::Bond::AROMATIC: bo = 1; break; // aromatic stored via flag
            default:                    bo = 1; break;
        }
        bool inRing = ri->numBondRings(bond->getIdx()) > 0;
        bool arom   = bond->getIsAromatic();

        // Forward edge: a → b
        neighbors[a].push_back(b);
        bondOrders[a].push_back(bo);
        bondRingFlags[a].push_back(inRing);
        bondAromFlags[a].push_back(arom);

        // Reverse edge: b → a
        neighbors[b].push_back(a);
        bondOrders[b].push_back(bo);
        bondRingFlags[b].push_back(inRing);
        bondAromFlags[b].push_back(arom);
    }

    // Tetrahedral chirality (CW = +1, CCW = -1, none = 0)
    std::vector<int> tetraChir(n, 0);
    for (int i = 0; i < n; ++i) {
        const auto* atom = m->getAtomWithIdx(i);
        auto chi = atom->getChiralTag();
        if (chi == RDKit::Atom::CHI_TETRAHEDRAL_CW)  tetraChir[i] = +1;
        if (chi == RDKit::Atom::CHI_TETRAHEDRAL_CCW) tetraChir[i] = -1;
    }

    // Build MolGraph via the Builder (handles dense/sparse bond storage,
    // label computation, adjLong, tautomers, and initDerivedFields).
    return MolGraph::Builder()
        .atomCount(n)
        .atomicNumbers(std::move(atomicNum))
        .formalCharges(std::move(formalCharge))
        .massNumbers(std::move(massNumber))
        .ringFlags(std::move(ringFlag))
        .aromaticFlags(std::move(aromFlag))
        .setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(bondOrders))
        .bondRingFlags(std::move(bondRingFlags))
        .bondAromaticFlags(std::move(bondAromFlags))
        .tetrahedralChirality(std::move(tetraChir))
        .build();
}

/// Parse SMILES string to MolGraph via RDKit.
MolGraph fromSmiles(const std::string& smiles) {
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) throw std::invalid_argument("Failed to parse SMILES: " + smiles);
    return fromRDKit(*mol);
}

} // namespace smsd

#endif // SMSD_WITH_RDKIT

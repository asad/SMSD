/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * RDKit adapter: converts RDKit ROMol to SMSD MolGraph.
 * This is the ONLY file that depends on RDKit.
 * The core algorithms (mol_graph.hpp, search.hpp, mcs.hpp) are RDKit-free.
 */

#ifdef SMSD_WITH_RDKIT

#include "smsd/rdkit_adapter.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace smsd {

/// Convert an RDKit ROMol to an SMSD MolGraph.
/// Strips explicit hydrogens by default. Each atom keeps the hydrogen count RDKit holds for it,
/// and double bonds keep the configuration given against their stereo atoms, so the graph
/// agrees with the one SMSD's SMILES parser builds for the same molecule.
MolGraph fromRDKit(const RDKit::ROMol& mol, bool removeHs) {
    // Work on a copy, so the caller's molecule is left unchanged.
    RDKit::RWMol molCopy(mol);
    // Give E and Z as cis and trans against the stereo atoms before hydrogens are removed:
    // removing a hydrogen that is a stereo atom puts its partner in its place and turns cis and
    // trans to match, but leaves E and Z as they were.
    for (auto* bond : molCopy.bonds()) {
        if (bond->getBondType() != RDKit::Bond::DOUBLE || bond->getStereoAtoms().size() != 2) continue;
        if (bond->getStereo() == RDKit::Bond::STEREOZ) bond->setStereo(RDKit::Bond::STEREOCIS);
        else if (bond->getStereo() == RDKit::Bond::STEREOE) bond->setStereo(RDKit::Bond::STEREOTRANS);
    }
    if (removeHs) RDKit::MolOps::removeHs(molCopy);
    // Hydrogen counts are taken only where they mean something and RDKit can compute them; a
    // query molecule, or one RDKit cannot cache properties for, keeps unknown counts.
    bool hydrogenCountsKnown = true;
    for (const auto* atom : molCopy.atoms()) if (atom->hasQuery()) hydrogenCountsKnown = false;
    for (const auto* bond : molCopy.bonds()) if (bond->hasQuery()) hydrogenCountsKnown = false;
    if (hydrogenCountsKnown) {
        try {
            molCopy.updatePropertyCache(false);
        } catch (const std::exception&) {
            hydrogenCountsKnown = false;
        }
    }
    const RDKit::ROMol* m = &molCopy;

    int n = m->getNumAtoms();

    // Collect atom properties
    std::vector<int> atomicNum(n), formalCharge(n), massNumber(n, 0), hydrogenCount(n, 0);
    std::vector<uint8_t> ringFlag(n, 0), aromFlag(n, 0);

    RDKit::MolOps::findSSSR(molCopy);
    const auto* ri = m->getRingInfo();

    for (int i = 0; i < n; ++i) {
        const auto* atom = m->getAtomWithIdx(i);
        atomicNum[i]    = atom->getAtomicNum();
        formalCharge[i] = atom->getFormalCharge();
        massNumber[i]   = atom->getIsotope();
        // Implicit and declared hydrogens; hydrogen atoms still in the graph are counted by it.
        if (hydrogenCountsKnown) hydrogenCount[i] = static_cast<int>(atom->getTotalNumHs(false));
        aromFlag[i]     = atom->getIsAromatic() ? 1 : 0;
        ringFlag[i]     = ri->numAtomRings(i) > 0 ? 1 : 0;
    }

    // Build adjacency + per-neighbour bond properties
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

    // The source uses a trailing virtual ligand; the graph can put it first.
    std::vector<int> tetraChir(n, 0);
    for (int i = 0; i < n; ++i) {
        const auto* atom = m->getAtomWithIdx(i);
        auto chi = atom->getChiralTag();
        if (chi == RDKit::Atom::CHI_TETRAHEDRAL_CW)  tetraChir[i] = 2;
        if (chi == RDKit::Atom::CHI_TETRAHEDRAL_CCW) tetraChir[i] = 1;
        if (tetraChir[i] && neighbors[i].size() == 3
            && !(i > 0 && neighbors[i][0] < i))
            tetraChir[i] = 3 - tetraChir[i];
    }

    // Double-bond configurations, held as SMSD's parser holds them: against the first other
    // neighbour of each end, 1 when those two lie on the same side and 2 when they do not.
    // A bond given only as E or Z without stereo atoms, or as unspecified or atropisomeric,
    // carries no configuration.
    std::vector<std::tuple<int, int, int>> doubleBondConfigurations;
    for (const auto& bond : m->bonds()) {
        if (bond->getBondType() != RDKit::Bond::DOUBLE) continue;
        const auto stereo = bond->getStereo();
        const auto& stereoAtoms = bond->getStereoAtoms();
        if ((stereo != RDKit::Bond::STEREOCIS && stereo != RDKit::Bond::STEREOTRANS) || stereoAtoms.size() != 2)
            continue;
        const int a = bond->getBeginAtomIdx(), b = bond->getEndAtomIdx();
        if (!m->getBondBetweenAtoms(a, stereoAtoms[0]) || !m->getBondBetweenAtoms(b, stereoAtoms[1])) continue;
        bool together = stereo == RDKit::Bond::STEREOCIS;
        for (int neighbor : neighbors[a]) {
            if (neighbor == b) continue;
            if (neighbor != stereoAtoms[0]) together = !together;
            break;
        }
        for (int neighbor : neighbors[b]) {
            if (neighbor == a) continue;
            if (neighbor != stereoAtoms[1]) together = !together;
            break;
        }
        doubleBondConfigurations.emplace_back(a, b, together ? 1 : 2);
    }

    // Build MolGraph via the Builder (handles dense/sparse bond storage,
    // label computation, adjLong, tautomers, and initDerivedFields).
    MolGraph::Builder builder;
    builder.atomCount(n)
        .atomicNumbers(std::move(atomicNum))
        .formalCharges(std::move(formalCharge))
        .massNumbers(std::move(massNumber))
        .ringFlags(std::move(ringFlag))
        .aromaticFlags(std::move(aromFlag))
        .setNeighbors(std::move(neighbors))
        .setBondOrders(std::move(bondOrders))
        .bondRingFlags(std::move(bondRingFlags))
        .bondAromaticFlags(std::move(bondAromFlags))
        .tetrahedralChirality(std::move(tetraChir));
    if (hydrogenCountsKnown) builder.hydrogenCounts(std::move(hydrogenCount));
    if (!doubleBondConfigurations.empty()) builder.doubleBondStereoSparse(std::move(doubleBondConfigurations));
    return builder.build();
}

/// Parse SMILES string to MolGraph via RDKit.
MolGraph fromSmiles(const std::string& smiles) {
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
    if (!mol) throw std::invalid_argument("Failed to parse SMILES: " + smiles);
    return fromRDKit(*mol);
}

} // namespace smsd

#endif // SMSD_WITH_RDKIT

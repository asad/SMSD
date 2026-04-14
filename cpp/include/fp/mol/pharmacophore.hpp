/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Pharmacophore feature classification for FCFP fingerprints.
 * Classifies atoms into: donor, acceptor, positive, negative, aromatic, hydrophobic.
 */
#pragma once
#ifndef FP_MOL_PHARMACOPHORE_HPP
#define FP_MOL_PHARMACOPHORE_HPP

#include "smsd/mol_graph.hpp"

namespace smsd {
namespace fp {
namespace mol {

/// Classify atom into pharmacophoric feature classes (FCFP invariant).
/// Returns bitmask: bit 0=donor, 1=acceptor, 2=positive, 3=negative, 4=aromatic, 5=hydrophobic.
inline int classifyPharmacophore(const MolGraph& g, int idx) {
    int z = g.atomicNum[idx];
    int charge = g.formalCharge[idx];
    bool arom = g.aromatic[idx];
    int deg = g.degree[idx];
    int features = 0;

    // H-bond donor: N-H, O-H, S-H (Rogers & Hahn 2010)
    if (z == 7 || z == 8 || z == 16) {
        int typVal = (z == 7) ? 3 : 2;
        if (deg < typVal && charge >= 0) features |= 1;
    }

    // H-bond acceptor
    bool isPyrroleTypeN = false;
    if (z == 7 && arom) {
        isPyrroleTypeN = (g.hydrogenCount[idx] > 0);
    }
    bool isAcceptorN = (z == 7 && charge <= 0 && !isPyrroleTypeN);
    bool isAcceptorS = (z == 16 && !arom);
    if (isAcceptorN || z == 8 || z == 9 || isAcceptorS) features |= 2;

    // Positive ionisable
    if (z == 7 && charge > 0) features |= 4;
    if (z == 7 && !arom && deg <= 3 && charge == 0) {
        bool isAmide = false;
        bool isAniline = false;
        for (int nb : g.neighbors[idx]) {
            if (g.atomicNum[nb] == 6) {
                if (g.aromatic[nb]) { isAniline = true; break; }
                for (int nb2 : g.neighbors[nb]) {
                    if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2) {
                        isAmide = true; break;
                    }
                }
            }
            if (isAmide) break;
        }
        if (!isAmide && !isAniline) features |= 4;
    }

    // Negative ionisable
    if (z == 8 && charge < 0) features |= 8;
    if (z == 8 && charge == 0 && deg == 1) {
        for (int nb : g.neighbors[idx]) {
            int nbZ = g.atomicNum[nb];
            if (nbZ == 6 && g.bondOrder(idx, nb) == 1) {
                for (int nb2 : g.neighbors[nb])
                    if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2)
                        features |= 8;
            }
            if (nbZ == 15 || nbZ == 16) features |= 8;
        }
    }

    // Aromatic
    if (arom) features |= 16;

    // Hydrophobic
    if (z == 6 && !arom) {
        bool hasHetero = false;
        for (int nb : g.neighbors[idx]) {
            int nbZ = g.atomicNum[nb];
            if (nbZ != 6 && nbZ != 1) { hasHetero = true; break; }
        }
        if (!hasHetero) features |= 32;
    }
    if (z == 17 || z == 35 || z == 53) features |= 32;

    return features;
}

} // namespace mol
} // namespace fp
} // namespace smsd

#endif // FP_MOL_PHARMACOPHORE_HPP

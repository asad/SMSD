/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
#pragma once
#ifndef SMSD_BOND_ENERGIES_HPP
#define SMSD_BOND_ENERGIES_HPP

#include "smsd/mol_graph.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>

namespace smsd {

// ============================================================================
// Bond dissociation energy (BDE) lookup for chemical bonds.
//
// All energies are in kJ/mol.  The table uses a compact integer key:
//   key = min(Z1,Z2) * 10000 + max(Z1,Z2) * 100 + bondOrder
// where Z1, Z2 are atomic numbers and bondOrder is 1, 2, 3, or 5 (aromatic).
//
// Data is stored as a sorted constexpr array for binary-search lookup.
// Thread-safe by construction: all data is const/constexpr.
// ============================================================================

namespace detail {

// Encode an (atomicNum1, atomicNum2, bondOrder) triple into a lookup key.
// Bond order convention: 1=single, 2=double, 3=triple, 5=aromatic.
constexpr int bdeKey(int z1, int z2, int bo) {
    int lo = (z1 < z2) ? z1 : z2;
    int hi = (z1 < z2) ? z2 : z1;
    return lo * 10000 + hi * 100 + bo;
}

struct BdeEntry {
    int    key;
    double energy;   // kJ/mol
};

// Sorted by key for binary search.
// Atomic numbers: H=1, C=6, N=7, O=8, F=9, Si=14, P=15, S=16, Cl=17, Br=35, I=53
constexpr BdeEntry BDE_TABLE[] = {
    // H-H (1-1)
    { bdeKey( 1,  1, 1),  436.0 },

    // B-H (1-5)
    { bdeKey( 1,  5, 1),  389.0 },

    // C-H (1-6)
    { bdeKey( 1,  6, 1),  413.0 },

    // N-H (1-7)
    { bdeKey( 1,  7, 1),  391.0 },

    // O-H (1-8)
    { bdeKey( 1,  8, 1),  463.0 },

    // F-H (1-9)
    { bdeKey( 1,  9, 1),  567.0 },

    // Si-H (1-14)
    { bdeKey( 1, 14, 1),  318.0 },

    // P-H (1-15)
    { bdeKey( 1, 15, 1),  322.0 },

    // S-H (1-16)
    { bdeKey( 1, 16, 1),  363.0 },

    // Cl-H (1-17)
    { bdeKey( 1, 17, 1),  431.0 },

    // Se-H (1-34)  [key=1341, between Cl-H(1171) and Br-H(1351)]
    { bdeKey( 1, 34, 1),  305.0 },

    // Br-H (1-35)
    { bdeKey( 1, 35, 1),  366.0 },

    // I-H (1-53)
    { bdeKey( 1, 53, 1),  298.0 },

    // B-C (5-6)
    { bdeKey( 5,  6, 1),  372.0 },

    // B-N (5-7)
    { bdeKey( 5,  7, 1),  444.0 },

    // B-O (5-8)
    { bdeKey( 5,  8, 1),  536.0 },

    // B-F (5-9)
    { bdeKey( 5,  9, 1),  644.0 },

    // B-Cl (5-17)
    { bdeKey( 5, 17, 1),  456.0 },

    // C-C (6-6)
    { bdeKey( 6,  6, 1),  348.0 },
    { bdeKey( 6,  6, 2),  614.0 },
    { bdeKey( 6,  6, 3),  839.0 },
    { bdeKey( 6,  6, 5),  518.0 },   // aromatic (between single and double)

    // C-N (6-7)
    { bdeKey( 6,  7, 1),  305.0 },
    { bdeKey( 6,  7, 2),  615.0 },
    { bdeKey( 6,  7, 3),  891.0 },
    { bdeKey( 6,  7, 5),  460.0 },   // aromatic C-N

    // C-O (6-8)
    { bdeKey( 6,  8, 1),  360.0 },
    { bdeKey( 6,  8, 2),  745.0 },
    { bdeKey( 6,  8, 5),  553.0 },   // aromatic C-O

    // C-F (6-9)
    { bdeKey( 6,  9, 1),  484.0 },

    // C-Si (6-14)
    { bdeKey( 6, 14, 1),  318.0 },

    // C-P (6-15)
    { bdeKey( 6, 15, 1),  264.0 },

    // C-S (6-16)
    { bdeKey( 6, 16, 1),  272.0 },
    { bdeKey( 6, 16, 2),  536.0 },  // literature range 536-573; conservative value

    // C-Cl (6-17)
    { bdeKey( 6, 17, 1),  338.0 },

    // C-Se (6-34)
    { bdeKey( 6, 34, 1),  245.0 },

    // C-Br (6-35)
    { bdeKey( 6, 35, 1),  276.0 },

    // C-I (6-53)
    { bdeKey( 6, 53, 1),  238.0 },

    // N-N (7-7)
    { bdeKey( 7,  7, 1),  163.0 },
    { bdeKey( 7,  7, 2),  418.0 },
    { bdeKey( 7,  7, 3),  941.0 },
    { bdeKey( 7,  7, 5),  291.0 },   // aromatic N-N

    // N-O (7-8)
    { bdeKey( 7,  8, 1),  201.0 },
    { bdeKey( 7,  8, 2),  607.0 },

    // N-F (7-9)
    { bdeKey( 7,  9, 1),  270.0 },

    // P-N (7-15)
    { bdeKey( 7, 15, 1),  230.0 },

    // N-Cl (7-17)
    { bdeKey( 7, 17, 1),  200.0 },

    // O-O (8-8)
    { bdeKey( 8,  8, 1),  146.0 },
    { bdeKey( 8,  8, 2),  498.0 },

    // O-F (8-9)
    { bdeKey( 8,  9, 1),  190.0 },

    // O-Si (8-14)
    { bdeKey( 8, 14, 1),  452.0 },

    // O-P (8-15)
    { bdeKey( 8, 15, 1),  335.0 },
    { bdeKey( 8, 15, 2),  544.0 },

    // O-S (8-16)
    { bdeKey( 8, 16, 1),  265.0 },
    { bdeKey( 8, 16, 2),  522.0 },

    // O-Cl (8-17)
    { bdeKey( 8, 17, 1),  203.0 },

    // Si-Si (14-14)
    { bdeKey(14, 14, 1),  226.0 },

    // P-P (15-15)
    { bdeKey(15, 15, 1),  201.0 },

    // P-S (15-16)
    { bdeKey(15, 16, 1),  230.0 },

    // S-S (16-16)
    { bdeKey(16, 16, 1),  226.0 },
    { bdeKey(16, 16, 2),  425.0 },

    // S-Cl (16-17)
    { bdeKey(16, 17, 1),  255.0 },

    // Cl-Cl (17-17)
    { bdeKey(17, 17, 1),  242.0 },

    // Se-Se (34-34)
    { bdeKey(34, 34, 1),  172.0 },

    // Br-Br (35-35)
    { bdeKey(35, 35, 1),  193.0 },

    // I-I (53-53)
    { bdeKey(53, 53, 1),  151.0 },
};

constexpr int BDE_TABLE_SIZE = sizeof(BDE_TABLE) / sizeof(BDE_TABLE[0]);

// Default energy for bonds not in the table.
constexpr double BDE_DEFAULT = 350.0;

// Binary search on the sorted constexpr table.
inline double lookupBDE(int key) {
    int lo = 0, hi = BDE_TABLE_SIZE - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (BDE_TABLE[mid].key == key) return BDE_TABLE[mid].energy;
        if (BDE_TABLE[mid].key < key) lo = mid + 1;
        else                          hi = mid - 1;
    }
    return BDE_DEFAULT;
}

} // namespace detail

// ============================================================================
// Public API
// ============================================================================

/// Effective bond order for energy lookup, mapping aromatic bonds to order 5.
inline int effectiveBondOrder(int bondOrder, bool isAromatic) {
    if (isAromatic) return 5;
    return bondOrder;
}

/// Bond dissociation energy for a given atom pair and bond order (kJ/mol).
/// Bond order: 1=single, 2=double, 3=triple, 5=aromatic.
/// Returns the default (350 kJ/mol) for unknown combinations.
inline double bondDissociationEnergy(int atomicNum1, int atomicNum2, int bondOrder) {
    int key = detail::bdeKey(atomicNum1, atomicNum2, bondOrder);
    return detail::lookupBDE(key);
}

/// Bond formation energy (negative BDE — energy released when bond forms).
inline double bondFormationEnergy(int atomicNum1, int atomicNum2, int bondOrder) {
    return -bondDissociationEnergy(atomicNum1, atomicNum2, bondOrder);
}

/// Cost of converting a bond between two atoms from one order to another.
/// Positive = net energy input, negative = net energy release.
/// If oldBondOrder == 0, this is purely formation cost.
/// If newBondOrder == 0, this is purely dissociation cost.
inline double bondChangeCost(int atomicNum1, int atomicNum2,
                             int oldBondOrder, int newBondOrder)
{
    if (oldBondOrder == newBondOrder) return 0.0;

    double breakCost = 0.0;
    double formGain  = 0.0;

    if (oldBondOrder > 0) {
        breakCost = bondDissociationEnergy(atomicNum1, atomicNum2, oldBondOrder);
    }
    if (newBondOrder > 0) {
        formGain = bondDissociationEnergy(atomicNum1, atomicNum2, newBondOrder);
    }

    // Net cost: energy to break old bond minus energy released forming new bond.
    return breakCost - formGain;
}

/// Sum of all bond dissociation energies in a molecule.
/// Each bond is counted once.
inline double totalBondEnergy(const MolGraph& mol) {
    double total = 0.0;
    for (int i = 0; i < mol.n; ++i) {
        for (int j : mol.neighbors[i]) {
            if (j <= i) continue;   // each edge once
            int bo = effectiveBondOrder(mol.bondOrder(i, j),
                                        mol.bondAromatic(i, j));
            total += bondDissociationEnergy(mol.atomicNum[i],
                                            mol.atomicNum[j], bo);
        }
    }
    return total;
}

/// Estimate the total bond-change cost for a given atom-atom mapping between
/// two molecules.  For every pair of mapped atoms (i -> mapping[i]) the
/// function compares the bond orders in g1 with those in g2.
///
/// Cost contributions:
///   - bond broken in g1 but absent in g2 (dissociation)
///   - bond absent in g1 but present in g2 (formation, counted as negative BDE)
///   - bond order change (difference of BDEs)
///
/// The mapping keys are g1 atom indices; values are g2 atom indices.
/// Unmapped atoms are ignored — only mapped pairs are evaluated.
inline double mappingBondChangeCost(const MolGraph& g1, const MolGraph& g2,
                                    const std::map<int,int>& mapping)
{
    double cost = 0.0;

    // Build reverse lookup for quick membership checks.
    std::unordered_set<int> mappedInG2;
    mappedInG2.reserve(mapping.size());
    for (auto& kv : mapping) mappedInG2.insert(kv.second);

    // Iterate over all mapped atom pairs in g1.
    for (auto it1 = mapping.begin(); it1 != mapping.end(); ++it1) {
        int a1 = it1->first;
        int b1 = it1->second;

        for (auto it2 = std::next(it1); it2 != mapping.end(); ++it2) {
            int a2 = it2->first;
            int b2 = it2->second;

            // Bond order in g1 between mapped atoms a1, a2.
            int bo1 = 0;
            if (a1 >= 0 && a1 < g1.n && a2 >= 0 && a2 < g1.n) {
                bo1 = effectiveBondOrder(g1.bondOrder(a1, a2),
                                         g1.bondAromatic(a1, a2));
            }

            // Bond order in g2 between corresponding atoms b1, b2.
            int bo2 = 0;
            if (b1 >= 0 && b1 < g2.n && b2 >= 0 && b2 < g2.n) {
                bo2 = effectiveBondOrder(g2.bondOrder(b1, b2),
                                         g2.bondAromatic(b1, b2));
            }

            if (bo1 == bo2) continue;

            // Use g1 atom types for the cost calculation (reactant perspective).
            int z1 = g1.atomicNum[a1];
            int z2 = g1.atomicNum[a2];
            cost += std::abs(bondChangeCost(z1, z2, bo1, bo2));
        }
    }

    return cost;
}

} // namespace smsd

#endif // SMSD_BOND_ENERGIES_HPP

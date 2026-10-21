/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * MolGraph construction: MolGraph::Builder::build() validates the supplied arrays
 * and fills a graph, and initDerivedFields() computes the eager derived data.
 *
 * Part of mol_graph.hpp, which includes this header; including it first includes
 * the whole of mol_graph.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

namespace smsd {

// ============================================================================
// Derived fields computed after construction
// ============================================================================

inline void MolGraph::initDerivedFields() {
    words = (n + 63) >> 6;
    adjLong.assign(n, std::vector<uint64_t>(words, 0));
    for (int i = 0; i < n; i++) {
        for (int k : neighbors[i])
            adjLong[i][k >> 6] |= uint64_t(1) << (k & 63);
    }
    computeTautomerClasses();
    applySolventCorrection(solvent_);
    adjustWeightsForPH(pH_);
    // Canonical labelling and ring counts are deferred; call
    // ensureCanonical() / ensureRingCounts() before accessing those fields.
}

// ============================================================================
// MolGraph::Builder::build(): validation and graph construction
// ============================================================================

inline MolGraph MolGraph::Builder::build(
    bool perceiveAromaticity,
    AromaticityModel model) const {
    if (n_ < 0)
        throw std::invalid_argument("atomCount must be >= 0");
    if (static_cast<int>(atomicNum_.size()) != n_)
        throw std::invalid_argument("atomicNumbers must have length == atomCount");
    if (static_cast<int>(neighbors_.size()) != n_)
        throw std::invalid_argument("neighbors must have length == atomCount");
    if (!formalCharge_.empty() && static_cast<int>(formalCharge_.size()) != n_)
        throw std::invalid_argument("formalCharges must have length == atomCount");
    if (!massNumber_.empty() && static_cast<int>(massNumber_.size()) != n_)
        throw std::invalid_argument("massNumbers must have length == atomCount");
    if (!hydrogenCount_.empty() && static_cast<int>(hydrogenCount_.size()) != n_)
        throw std::invalid_argument("hydrogenCounts must have length == atomCount");
    if (!atomClass_.empty() && static_cast<int>(atomClass_.size()) != n_)
        throw std::invalid_argument("atomClasses must have length == atomCount");
    if (!ring_.empty() && static_cast<int>(ring_.size()) != n_)
        throw std::invalid_argument("ringFlags must have length == atomCount");
    if (!aromatic_.empty() && static_cast<int>(aromatic_.size()) != n_)
        throw std::invalid_argument("aromaticFlags must have length == atomCount");
    if (!tetraChirality_.empty() && static_cast<int>(tetraChirality_.size()) != n_)
        throw std::invalid_argument("tetrahedralChirality must have length == atomCount");
    if (!bondOrders_.empty() && static_cast<int>(bondOrders_.size()) != n_)
        throw std::invalid_argument("bondOrders must have length == atomCount");
    if (!bondRings_.empty() && static_cast<int>(bondRings_.size()) != n_)
        throw std::invalid_argument("bondRingFlags must have length == atomCount");
    if (!bondAroms_.empty() && static_cast<int>(bondAroms_.size()) != n_)
        throw std::invalid_argument("bondAromaticFlags must have length == atomCount");
    if (!dbStereoConf_.empty() && static_cast<int>(dbStereoConf_.size()) != n_)
        throw std::invalid_argument("doubleBondStereo must have length == atomCount");

    for (int i = 0; i < n_; ++i) {
        if (atomicNum_[i] < 0 || atomicNum_[i] > 118)
            throw std::invalid_argument("atomicNumbers must be in [0, 118]");
        if (!massNumber_.empty() && massNumber_[i] < 0)
            throw std::invalid_argument("massNumbers must be >= 0");
        if (!hydrogenCount_.empty() && hydrogenCount_[i] < 0)
            throw std::invalid_argument("hydrogenCounts must be >= 0");
        if (!atomClass_.empty() && atomClass_[i] < 0)
            throw std::invalid_argument("atomClasses must be >= 0");
        if (!ring_.empty() && ring_[i] > 1)
            throw std::invalid_argument("ringFlags entries must be 0 or 1");
        if (!aromatic_.empty() && aromatic_[i] > 1)
            throw std::invalid_argument("aromaticFlags entries must be 0 or 1");
        if (!tetraChirality_.empty()
            && tetraChirality_[i] != 0
            && tetraChirality_[i] != 1
            && tetraChirality_[i] != 2)
            throw std::invalid_argument("tetrahedralChirality entries must be 0, 1, or 2");
    }

    std::vector<std::unordered_map<int, int>> neighborIndex(n_);
    for (int i = 0; i < n_; ++i) {
        neighborIndex[i].reserve(neighbors_[i].size());
        for (int k = 0; k < static_cast<int>(neighbors_[i].size()); ++k) {
            int nb = neighbors_[i][k];
            if (nb < 0 || nb >= n_)
                throw std::invalid_argument("neighbour index out of range");
            if (nb == i)
                throw std::invalid_argument("self-loop detected in neighbors");
            if (!neighborIndex[i].emplace(nb, k).second)
                throw std::invalid_argument("duplicate neighbour detected in adjacency list");
        }
    }

    auto classifyIntRows = [&](const std::vector<std::vector<int>>& rows,
                               const char* name) {
        std::vector<RowEncoding> enc(n_, RowEncoding::ABSENT);
        if (rows.empty()) return enc;
        for (int i = 0; i < n_; ++i) {
            const auto& row = rows[i];
            const auto deg = neighbors_[i].size();
            if (row.size() == deg) enc[i] = RowEncoding::PARALLEL;
            else if (row.size() == static_cast<size_t>(n_)) enc[i] = RowEncoding::DENSE;
            else
                throw std::invalid_argument(std::string(name)
                    + " rows must have length == neighbour count or atomCount");
        }
        return enc;
    };

    auto classifyBoolRows = [&](const std::vector<std::vector<bool>>& rows,
                                const char* name) {
        std::vector<RowEncoding> enc(n_, RowEncoding::ABSENT);
        if (rows.empty()) return enc;
        for (int i = 0; i < n_; ++i) {
            const auto& row = rows[i];
            const auto deg = neighbors_[i].size();
            if (row.size() == deg) enc[i] = RowEncoding::PARALLEL;
            else if (row.size() == static_cast<size_t>(n_)) enc[i] = RowEncoding::DENSE;
            else
                throw std::invalid_argument(std::string(name)
                    + " rows must have length == neighbour count or atomCount");
        }
        return enc;
    };

    auto bondOrderEnc = classifyIntRows(bondOrders_, "bondOrders");
    auto bondRingEnc = classifyBoolRows(bondRings_, "bondRingFlags");
    auto bondAromEnc = classifyBoolRows(bondAroms_, "bondAromaticFlags");

    auto intPropAt = [&](const std::vector<std::vector<int>>& rows,
                         const std::vector<RowEncoding>& enc,
                         int i, int k, int j, int fallback) {
        if (rows.empty()) return fallback;
        return enc[i] == RowEncoding::DENSE ? rows[i][j] : rows[i][k];
    };
    auto boolPropAt = [&](const std::vector<std::vector<bool>>& rows,
                          const std::vector<RowEncoding>& enc,
                          int i, int k, int j, bool fallback) -> bool {
        if (rows.empty()) return fallback;
        return enc[i] == RowEncoding::DENSE ? rows[i][j] : rows[i][k];
    };

    for (int i = 0; i < n_; ++i) {
        for (const auto& [nb, k] : neighborIndex[i]) {
            auto rev = neighborIndex[nb].find(i);
            if (rev == neighborIndex[nb].end())
                throw std::invalid_argument("adjacency must be symmetric");
            if (i > nb) continue;
            int rk = rev->second;
            if (!bondOrders_.empty()) {
                int forward = intPropAt(bondOrders_, bondOrderEnc, i, k, nb, 1);
                int backward = intPropAt(bondOrders_, bondOrderEnc, nb, rk, i, 1);
                if (forward <= 0 || backward <= 0)
                    throw std::invalid_argument("bondOrders for present edges must be > 0");
                if (forward != backward)
                    throw std::invalid_argument("bondOrders must be symmetric");
            }
            if (!bondRings_.empty()) {
                bool forward = boolPropAt(bondRings_, bondRingEnc, i, k, nb, false);
                bool backward = boolPropAt(bondRings_, bondRingEnc, nb, rk, i, false);
                if (forward != backward)
                    throw std::invalid_argument("bondRingFlags must be symmetric");
            }
            if (!bondAroms_.empty()) {
                bool forward = boolPropAt(bondAroms_, bondAromEnc, i, k, nb, false);
                bool backward = boolPropAt(bondAroms_, bondAromEnc, nb, rk, i, false);
                if (forward != backward)
                    throw std::invalid_argument("bondAromaticFlags must be symmetric");
            }
        }
    }

    if (!dbStereoConf_.empty()) {
        for (int i = 0; i < n_; ++i) {
            if (static_cast<int>(dbStereoConf_[i].size()) != n_)
                throw std::invalid_argument("doubleBondStereo rows must have length == atomCount");
        }
        for (int i = 0; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                int conf = dbStereoConf_[i][j];
                if (conf < 0 || conf > 2)
                    throw std::invalid_argument("doubleBondStereo entries must be 0, 1, or 2");
                if (i == j && conf != 0)
                    throw std::invalid_argument("doubleBondStereo diagonal must be zero");
                if (j <= i) continue;
                if (conf != dbStereoConf_[j][i])
                    throw std::invalid_argument("doubleBondStereo must be symmetric");
                if (conf != 0) {
                    auto it = neighborIndex[i].find(j);
                    if (it == neighborIndex[i].end())
                        throw std::invalid_argument("doubleBondStereo requires a bond");
                    int ord = intPropAt(bondOrders_, bondOrderEnc, i, it->second, j, 1);
                    if (ord != 2 && ord != 4)
                        throw std::invalid_argument("doubleBondStereo requires a double or aromatic bond");
                }
            }
        }
    }

    std::unordered_map<int64_t, int> sparseConfigurations;
    for (const auto& entry : sparseDbStereoConf_) {
        const auto [first, second, configuration] = entry;
        if (first < 0 || first >= n_ || second < 0 || second >= n_ || first == second)
            throw std::invalid_argument("doubleBondStereoSparse requires distinct atom indices in range");
        if (configuration != 1 && configuration != 2)
            throw std::invalid_argument("doubleBondStereoSparse configurations must be 1 or 2");
        const auto neighbor = neighborIndex[first].find(second);
        if (neighbor == neighborIndex[first].end())
            throw std::invalid_argument("doubleBondStereoSparse requires a bond");
        const int order = intPropAt(bondOrders_, bondOrderEnc, first, neighbor->second, second, 1);
        if (order != 2 && order != 4)
            throw std::invalid_argument("doubleBondStereoSparse requires a double or aromatic bond");
        const auto inserted = sparseConfigurations.emplace(bondKey(first, second), configuration);
        if (!inserted.second && inserted.first->second != configuration)
            throw std::invalid_argument("doubleBondStereoSparse has conflicting configurations");
    }

    if (!atomIds_.empty()) {
        if (static_cast<int>(atomIds_.size()) != n_)
            throw std::invalid_argument("atomIds must have length == atomCount");
        std::unordered_set<int> seenIds;
        seenIds.reserve(atomIds_.size());
        for (int id : atomIds_) {
            if (id <= 0)
                throw std::invalid_argument("atomIds must be positive");
            if (!seenIds.insert(id).second)
                throw std::invalid_argument("atomIds must be unique");
        }
    }

    MolGraph g;
    g.n = n_;
    g.words = (n_ + 63) >> 6;

    // Copy atom arrays
    g.atomicNum    = atomicNum_;
    g.formalCharge = formalCharge_.empty()
                     ? std::vector<int>(n_, 0) : formalCharge_;
    g.massNumber   = massNumber_.empty()
                     ? std::vector<int>(n_, 0) : massNumber_;
    g.hydrogenCount = hydrogenCount_.empty()
                     ? std::vector<int>(n_, 0) : hydrogenCount_;
    g.hydrogenCountsKnown = !hydrogenCount_.empty();
    g.atomClass    = atomClass_.empty()
                     ? std::vector<int>(n_, 0) : atomClass_;

    if (ring_.empty())     g.ring.assign(n_, uint8_t(0));     else g.ring = ring_;
    if (aromatic_.empty()) g.aromatic.assign(n_, uint8_t(0)); else g.aromatic = aromatic_;

    // Neighbours + degree
    g.neighbors = neighbors_;
    g.degree.resize(n_);
    for (int i = 0; i < n_; i++)
        g.degree[i] = static_cast<int>(g.neighbors[i].size());

    // Bond storage -- dense or sparse
    if (n_ <= SPARSE_THRESHOLD) {
        g.useDense = true;
        g.bondPacked.assign(n_, std::vector<uint8_t>(n_, 0));
        for (int i = 0; i < n_; i++) {
            const auto& nbs = g.neighbors[i];
            for (int k = 0; k < static_cast<int>(nbs.size()); k++) {
                int j = nbs[k];
                int ord = intPropAt(bondOrders_, bondOrderEnc, i, k, j, 1);
                bool inRing = boolPropAt(bondRings_, bondRingEnc, i, k, j, false);
                bool arom = boolPropAt(bondAroms_, bondAromEnc, i, k, j, false);
                g.bondPacked[i][j] = static_cast<uint8_t>((ord & 0x07)
                    | (inRing ? BOND_RING_BIT : 0)
                    | (arom ? BOND_AROM_BIT : 0));
            }
        }
    } else {
        g.useDense = false;
        for (int i = 0; i < n_; i++) {
            const auto& nbs = g.neighbors[i];
            for (int k = 0; k < static_cast<int>(nbs.size()); k++) {
                int j = nbs[k];
                if (i >= j) continue;   // store each edge once
                int ord = intPropAt(bondOrders_, bondOrderEnc, i, k, j, 1);
                int inRing = boolPropAt(bondRings_, bondRingEnc, i, k, j, false) ? 1 : 0;
                int arom = boolPropAt(bondAroms_, bondAromEnc, i, k, j, false) ? 1 : 0;
                g.sparseBondProps[bondKey(i, j)] = {ord, inRing, arom};
            }
        }
    }

    // Stereochemistry
    g.tetraChirality = tetraChirality_.empty()
                       ? std::vector<int>(n_, 0) : tetraChirality_;
    if (n_ <= SPARSE_THRESHOLD) {
        g.hasDbStereo = true;
        g.dbStereoConf.assign(n_, std::vector<int>(n_, 0));
        if (!dbStereoConf_.empty()) {
            for (int i = 0; i < n_
                 && i < static_cast<int>(dbStereoConf_.size()); i++) {
                for (int j = 0; j < n_
                     && j < static_cast<int>(dbStereoConf_[i].size()); j++) {
                    g.dbStereoConf[i][j] = dbStereoConf_[i][j];
                }
            }
        }
    } else if (!dbStereoConf_.empty()) {
        g.hasDbStereo = true;
        for (int i = 0; i < n_; ++i) {
            for (int j = i + 1; j < n_; ++j) {
                int conf = dbStereoConf_[i][j];
                if (conf != 0)
                    g.sparseDbStereo[bondKey(i, j)] = conf;
            }
        }
    }
    if (!sparseConfigurations.empty()) {
        g.hasDbStereo = true;
        if (g.useDense) {
            for (const auto& entry : sparseDbStereoConf_) {
                const auto [first, second, configuration] = entry;
                g.dbStereoConf[first][second] = g.dbStereoConf[second][first] = configuration;
            }
        } else {
            g.sparseDbStereo = std::move(sparseConfigurations);
        }
    }

    // External atom IDs (optional)
    if (!atomIds_.empty()) g.atomId = atomIds_;

    g.name = name_;
    g.programLine = programLine_;
    g.comment = comment_;
    g.properties = properties_;

    if (perceiveAromaticity) g.normalizeRingAndAromaticity(model);
    g.refreshAtomLabels();

    g.ringCount.assign(n_, 0);
    g.solvent_ = solvent_;
    g.pH_ = pH_;
    g.initDerivedFields();
    return g;
}

} // namespace smsd

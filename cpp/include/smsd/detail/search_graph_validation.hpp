/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#pragma once

#include "smsd/mol_graph.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace smsd {
namespace search_detail {

/** Checks sizes before a search uses the atom count to construct its budget or buffers. */
inline void validateSearchGraphShape(const MolGraph& graph) {
    if (graph.n < 0) throw std::invalid_argument("Search graph atom count must be non-negative");
    const size_t atoms = static_cast<size_t>(graph.n);
    if (graph.atomicNum.size() < atoms || graph.aromatic.size() < atoms || graph.ring.size() < atoms
        || graph.degree.size() < atoms || graph.neighbors.size() < atoms)
        throw std::invalid_argument("Search graph atom and adjacency arrays must cover its atom count");
    if (!graph.tautomerClass.empty() && graph.tautomerClass.size() < atoms)
        throw std::invalid_argument("Populated search graph tautomer classes must cover its atom count");
    if (graph.useDense && graph.bondPacked.size() < atoms)
        throw std::invalid_argument("Search graph dense bond rows must cover its atom count");
    if (graph.useDense && graph.hasDbStereo && graph.dbStereoConf.size() < atoms)
        throw std::invalid_argument("Search graph stereo rows must cover its atom count");
}

/** Checks array storage and reciprocal bond presence, allowing repeated and self listings. */
template<class Stop>
inline bool validateSearchGraphStorage(const MolGraph& graph, Stop stop) {
    validateSearchGraphShape(graph);
    const size_t atoms = static_cast<size_t>(graph.n);
    for (size_t atom = 0; atom < atoms; ++atom) {
        if (stop()) return false;
        if (graph.useDense && graph.bondPacked[atom].size() < atoms)
            throw std::invalid_argument("Search graph dense bond columns must cover its atom count");
        if (graph.useDense && graph.hasDbStereo && graph.dbStereoConf[atom].size() < atoms)
            throw std::invalid_argument("Search graph stereo columns must cover its atom count");
        for (int neighbour : graph.neighbors[atom]) {
            if (stop()) return false;
            if (neighbour < 0 || neighbour >= graph.n)
                throw std::out_of_range("Search graph neighbour index out of range");
        }
    }
    std::vector<int> seen(atoms, -1);
    size_t listed = 0;
    for (int atom = 0; atom < graph.n; ++atom) {
        if (stop()) return false;
        for (int neighbour : graph.neighbors[atom]) {
            if (stop()) return false;
            if (neighbour == atom || seen[neighbour] == atom) continue;
            seen[neighbour] = atom;
            if (!graph.hasBond(atom, neighbour) || !graph.hasBond(neighbour, atom))
                throw std::invalid_argument("Search graph neighbours must agree with reciprocal bond presence");
            ++listed;
        }
    }
    size_t present = 0;
    if (graph.useDense) {
        for (int atom = 0; atom < graph.n; ++atom) {
            for (int other = 0; other < graph.n; ++other) {
                if (stop()) return false;
                if (!graph.hasBond(atom, other)) continue;
                if (atom == other)
                    throw std::invalid_argument("Search graph cannot contain a self bond");
                ++present;
            }
        }
    } else {
        for (const auto& entry : graph.sparseBondProps) {
            if (stop()) return false;
            const uint64_t key = static_cast<uint64_t>(entry.first);
            const uint64_t first = key >> 32;
            const uint32_t second = static_cast<uint32_t>(key);
            if (entry.first < 0 || first >= static_cast<uint64_t>(graph.n)
                || second >= static_cast<uint32_t>(graph.n) || first > second || entry.second[0] < 0)
                throw std::invalid_argument("Search graph sparse bond entry is outside its supported range");
            if (entry.second[0] == 0) continue;
            if (first == second)
                throw std::invalid_argument("Search graph cannot contain a self bond");
            if (present > std::numeric_limits<size_t>::max() - 2)
                throw std::invalid_argument("Search graph has too many sparse bonds");
            present += 2;
        }
    }
    if (present != listed)
        throw std::invalid_argument("Search graph bond presence must agree with its neighbour lists");
    return true;
}


/** Whether graph storage supports the optional labelled-graph symmetry rules. */
template<class Stop>
inline bool symmetryGraphCoherent(const MolGraph& graph, Stop stop) {
    if (stop() || graph.n < 0) return false;
    const size_t atoms = static_cast<size_t>(graph.n);
    if (graph.atomicNum.size() < atoms || graph.aromatic.size() < atoms || graph.ring.size() < atoms
        || graph.degree.size() < atoms || graph.neighbors.size() < atoms
        || (!graph.tautomerClass.empty() && graph.tautomerClass.size() < atoms)
        || (graph.useDense && graph.bondPacked.size() < atoms)
        || (graph.useDense && graph.hasDbStereo && graph.dbStereoConf.size() < atoms)) return false;
    // Check all indices and storage before reading either direction of any bond.
    for (size_t atom = 0; atom < atoms; ++atom) {
        if (stop()) return false;
        if ((graph.useDense && graph.bondPacked[atom].size() < atoms)
            || (graph.useDense && graph.hasDbStereo && graph.dbStereoConf[atom].size() < atoms)
            || graph.degree[atom] < 0
            || static_cast<size_t>(graph.degree[atom]) != graph.neighbors[atom].size()) return false;
        for (int neighbour : graph.neighbors[atom]) {
            if (stop()) return false;
            if (neighbour < 0 || neighbour >= graph.n || neighbour == static_cast<int>(atom)) return false;
        }
    }
    std::vector<int> seen(atoms, -1);
    uint64_t listed = 0;
    for (int atom = 0; atom < graph.n; ++atom) {
        if (stop()) return false;
        for (int neighbour : graph.neighbors[atom]) {
            if (stop() || seen[neighbour] == atom) return false;
            seen[neighbour] = atom;
            if (graph.bondOrder(atom, neighbour) <= 0
                || graph.bondOrder(atom, neighbour) != graph.bondOrder(neighbour, atom)
                || graph.bondAromatic(atom, neighbour) != graph.bondAromatic(neighbour, atom)
                || graph.bondInRing(atom, neighbour) != graph.bondInRing(neighbour, atom)) return false;
            ++listed;
        }
    }
    uint64_t present = 0;
    if (graph.useDense) {
        for (int atom = 0; atom < graph.n; ++atom)
            for (int other = 0; other < graph.n; ++other) {
                if (stop()) return false;
                if (!graph.hasBond(atom, other)) continue;
                if (atom == other) return false;
                ++present;
            }
    } else {
        for (const auto& entry : graph.sparseBondProps) {
            if (stop()) return false;
            const uint64_t key = static_cast<uint64_t>(entry.first);
            const uint64_t first = key >> 32;
            const uint32_t second = static_cast<uint32_t>(key);
            if (entry.first < 0 || first >= static_cast<uint64_t>(graph.n)
                || second >= static_cast<uint32_t>(graph.n) || first > second || entry.second[0] < 0)
                return false;
            if (entry.second[0] == 0) continue;
            if (first == second) return false;
            present += 2;
        }
    }
    return !stop() && present == listed;
}

/** An owned query is checked before its snapshot is exposed to any search. */
inline void validateSearchGraphStorage(const MolGraph& graph) {
    (void)validateSearchGraphStorage(graph, [] { return false; });
}

} // namespace search_detail
} // namespace smsd

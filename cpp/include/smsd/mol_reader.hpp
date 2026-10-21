/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Header-only MDL MOL V2000 / SDF file reader and writer for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 * File paths use UTF-8 on every platform.
 */
#pragma once
#ifndef SMSD_MOL_READER_HPP
#define SMSD_MOL_READER_HPP

#include "smsd/mol_graph.hpp"
#include "smsd/layout.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// Internal implementation details for MOL reader
// ============================================================================
namespace mol_detail {

// --------------------------------------------------------------------------
// Element symbols: the table of smiles_utils.hpp. MOL files also write D and T
// for the hydrogen isotopes and never the aromatic lower-case SMILES forms.
// --------------------------------------------------------------------------
inline int elementToAtomicNum(const std::string& sym) {
    if (sym == "D" || sym == "T") return 1;
    if (!sym.empty() && sym[0] >= 'a' && sym[0] <= 'z') return -1;
    return detail::elementToAtomicNum(sym);
}

// --------------------------------------------------------------------------
// String utilities
// --------------------------------------------------------------------------
inline std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) --end;
    return s.substr(start, end - start);
}

inline int safeInt(const std::string& s, int start, int len, int defaultVal = 0) {
    if (start < 0 || start >= static_cast<int>(s.size())) return defaultVal;
    int end = std::min(start + len, static_cast<int>(s.size()));
    std::string sub = trim(s.substr(start, end - start));
    if (sub.empty()) return defaultVal;
    try {
        return std::stoi(sub);
    } catch (...) {
        return defaultVal;
    }
}

inline double safeDouble(const std::string& s, int start, int len, double defaultVal = 0.0) {
    if (start < 0 || start >= static_cast<int>(s.size())) return defaultVal;
    int end = std::min(start + len, static_cast<int>(s.size()));
    std::string sub = trim(s.substr(start, end - start));
    if (sub.empty()) return defaultVal;
    try {
        return std::stod(sub);
    } catch (...) {
        return defaultVal;
    }
}

// Split a string into lines
inline std::vector<std::string> splitLines(const std::string& s) {
    std::vector<std::string> lines;
    std::istringstream iss(s);
    std::string line;
    while (std::getline(iss, line)) {
        // Strip trailing \r for Windows line endings
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(std::move(line));
    }
    return lines;
}

inline std::vector<std::string> splitWhitespace(const std::string& s) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok) tokens.push_back(tok);
    return tokens;
}

// --------------------------------------------------------------------------
// V2000 charge mapping from the atom-block ccc field
// --------------------------------------------------------------------------
inline int v2000ChargeFromField(int ccc) {
    // V2000 spec: 0=uncharged, 1=+3, 2=+2, 3=+1, 4=doublet radical, 5=-1, 6=-2, 7=-3
    switch (ccc) {
        case 0: return 0;
        case 1: return 3;
        case 2: return 2;
        case 3: return 1;
        case 4: return 0;   // doublet radical, not a charge
        case 5: return -1;
        case 6: return -2;
        case 7: return -3;
        default: return 0;
    }
}

// --------------------------------------------------------------------------
// V2000 charge field from actual charge (for writing)
// --------------------------------------------------------------------------
inline int v2000ChargeToField(int charge) {
    switch (charge) {
        case  3: return 1;
        case  2: return 2;
        case  1: return 3;
        case  0: return 0;
        case -1: return 5;
        case -2: return 6;
        case -3: return 7;
        default: return 0;
    }
}

// --------------------------------------------------------------------------
// Parsed atom/bond structures for MOL reading
// --------------------------------------------------------------------------
struct MolAtom {
    double x = 0.0, y = 0.0, z = 0.0;
    std::string symbol;
    int atomicNum = 0;
    int formalCharge = 0;
    int massNumber = 0;
    int hCount = -1;
    int valence = -1;
    int radical = 0;
    int atomClass = 0;
    int tetraCfg = 0;
};

struct MolBond {
    int from = 0;      // 0-based
    int to = 0;        // 0-based
    int order = 1;     // 1=single, 2=double, 3=triple, 4=aromatic
    int stereo = 0;
};

/** Use ordinary Kekule bonds so implicit H does not require an HCOUNT query. */
inline MolGraph moleculeForWrite(const MolGraph& graph) {
    MolGraph encoded;
    {
        const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
        encoded = graph;
    }
    if (std::any_of(encoded.aromatic.begin(), encoded.aromatic.end(), [](uint8_t value) { return value != 0; })
        && !encoded.kekulize()) {
        for (int atom = 0; atom < graph.n; ++atom)
            for (int neighbor : graph.neighbors[atom])
                if (neighbor > atom && graph.bondOrder(atom, neighbor) == 4 && graph.dbStereo(atom, neighbor) != 0)
                    throw std::invalid_argument("Cannot preserve aromatic double-bond stereo in a Kekule form");
        const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
        return graph;
    }
    return encoded;
}

/** Preserve non-default hydrogen counts using the ordinary atom valence field. */
inline int valenceForWrite(const MolGraph& graph, int atom) {
    const int valence = detail::explicitValence(graph, atom);
    const int inferredHydrogens = graph.atomicNum[atom] == 1
        ? (graph.formalCharge[atom] == 0 ? std::max(0, 1 - valence) : 0)
        : graph.atomicNum[atom] == 0 ? 0
        : detail::computeImplicitH(graph.atomicNum[atom], graph.aromatic[atom], valence, graph.formalCharge[atom]);
    return inferredHydrogens == graph.hydrogenCount[atom] ? -1 : valence + graph.hydrogenCount[atom];
}

/** Convert between MDL's ascending ligand order and the graph's winding order. */
inline int tetrahedralParity(int parity, int centre, const std::vector<int>& neighbors,
                            const std::vector<int>& atomicNumbers) {
    if (parity != 1 && parity != 2) return 0;
    std::vector<int> mdlOrder = neighbors;
    std::sort(mdlOrder.begin(), mdlOrder.end(), [&](int first, int second) {
        const bool firstHydrogen = atomicNumbers[first] == 1;
        const bool secondHydrogen = atomicNumbers[second] == 1;
        return firstHydrogen != secondHydrogen ? !firstHydrogen : first < second;
    });
    std::vector<int> graphOrder = neighbors;
    if (neighbors.size() == 3) {
        mdlOrder.push_back(-1);
        const size_t hydrogenPosition = centre > 0 && neighbors[0] < centre ? 1 : 0;
        graphOrder.insert(graphOrder.begin() + hydrogenPosition, -1);
    }
    int inversions = 0;
    for (size_t first = 0; first < graphOrder.size(); ++first) {
        const auto firstPosition = std::find(mdlOrder.begin(), mdlOrder.end(), graphOrder[first]);
        for (size_t second = first + 1; second < graphOrder.size(); ++second) {
            const auto secondPosition = std::find(mdlOrder.begin(), mdlOrder.end(), graphOrder[second]);
            if (firstPosition > secondPosition) ++inversions;
        }
    }
    // MDL parity 1 has the opposite winding to SMILES @ in the same order.
    return inversions % 2 == 0 ? 3 - parity : parity;
}

/** Read the relation between the first explicit ligand at each double-bond end. */
inline int doubleBondGeometry(const std::vector<MolAtom>& atoms, int first, int second,
                              int firstLigand, int secondLigand) {
    const auto& firstAtom = atoms[first];
    const auto& secondAtom = atoms[second];
    const auto& firstReference = atoms[firstLigand];
    const auto& secondReference = atoms[secondLigand];
    const double axisX = secondAtom.x - firstAtom.x;
    const double axisY = secondAtom.y - firstAtom.y;
    const double axisZ = secondAtom.z - firstAtom.z;
    const double firstX = firstReference.x - firstAtom.x;
    const double firstY = firstReference.y - firstAtom.y;
    const double firstZ = firstReference.z - firstAtom.z;
    const double secondX = secondReference.x - secondAtom.x;
    const double secondY = secondReference.y - secondAtom.y;
    const double secondZ = secondReference.z - secondAtom.z;
    const double firstNormalX = axisY * firstZ - axisZ * firstY;
    const double firstNormalY = axisZ * firstX - axisX * firstZ;
    const double firstNormalZ = axisX * firstY - axisY * firstX;
    const double secondNormalX = axisY * secondZ - axisZ * secondY;
    const double secondNormalY = axisZ * secondX - axisX * secondZ;
    const double secondNormalZ = axisX * secondY - axisY * secondX;
    const double dot = firstNormalX * secondNormalX + firstNormalY * secondNormalY + firstNormalZ * secondNormalZ;
    if (std::abs(dot) < 1e-10) return 0;
    return dot > 0 ? 1 : 2;
}

/** Determine the SMILES winding from ligand vectors in the graph's order. */
inline int coordinateChirality(int centre, const std::vector<int>& neighbors,
                               const std::vector<MolAtom>& positions) {
    if (neighbors.size() != 3 && neighbors.size() != 4) return 0;
    std::vector<MolAtom> ligands;
    MolAtom hydrogen;
    for (int neighbor : neighbors) {
        MolAtom ligand;
        ligand.x = positions[neighbor].x - positions[centre].x;
        ligand.y = positions[neighbor].y - positions[centre].y;
        ligand.z = positions[neighbor].z - positions[centre].z;
        ligands.push_back(ligand);
        hydrogen.x -= ligand.x; hydrogen.y -= ligand.y; hydrogen.z -= ligand.z;
    }
    if (neighbors.size() == 3) {
        const size_t hydrogenPosition = centre > 0 && neighbors[0] < centre ? 1 : 0;
        ligands.insert(ligands.begin() + hydrogenPosition, hydrogen);
    }
    const auto& last = ligands[3];
    const double firstX = ligands[0].x - last.x, firstY = ligands[0].y - last.y, firstZ = ligands[0].z - last.z;
    const double secondX = ligands[1].x - last.x, secondY = ligands[1].y - last.y, secondZ = ligands[1].z - last.z;
    const double thirdX = ligands[2].x - last.x, thirdY = ligands[2].y - last.y, thirdZ = ligands[2].z - last.z;
    const double determinant = firstX * (secondY * thirdZ - secondZ * thirdY)
        - firstY * (secondX * thirdZ - secondZ * thirdX)
        + firstZ * (secondX * thirdY - secondY * thirdX);
    if (std::abs(determinant) < 1e-10) return 0;
    return determinant > 0 ? 1 : 2;
}

/** Allocate one directed wedge to each tetrahedral centre. */
inline std::unordered_map<int64_t, std::pair<int, int>> tetrahedralWedges(
        const MolGraph& graph, const std::vector<Point2D>& coordinates) {
    std::vector<MolAtom> positions(graph.n);
    for (int atom = 0; atom < graph.n; ++atom) {
        positions[atom].x = coordinates[atom].x;
        positions[atom].y = coordinates[atom].y;
    }
    std::vector<std::vector<std::pair<int, int>>> candidates(graph.n);
    for (int centre = 0; centre < graph.n; ++centre) {
        if (graph.tetraChirality[centre] == 0) continue;
        auto neighbors = graph.neighbors[centre];
        std::stable_sort(neighbors.begin(), neighbors.end(), [&](int first, int second) {
            if ((graph.tetraChirality[first] == 0) != (graph.tetraChirality[second] == 0))
                return graph.tetraChirality[first] == 0;
            return graph.neighbors[first].size() < graph.neighbors[second].size();
        });
        for (int neighbor : neighbors) {
            if (graph.bondOrder(centre, neighbor) != 1) continue;
            positions[neighbor].z = 1;
            const int winding = coordinateChirality(centre, graph.neighbors[centre], positions);
            positions[neighbor].z = 0;
            if (winding != 0) candidates[centre].push_back({neighbor, winding == graph.tetraChirality[centre] ? 1 : 6});
        }
    }
    std::unordered_map<int64_t, std::pair<int, int>> wedges;
    std::function<bool(int, std::vector<bool>&)> assignWedge = [&](int centre, std::vector<bool>& visited) {
        if (visited[centre]) return false;
        visited[centre] = true;
        for (const auto& candidate : candidates[centre]) {
            const int64_t key = MolGraph::bondKey(centre, candidate.first);
            const auto found = wedges.find(key);
            if (found == wedges.end() || assignWedge(found->second.first, visited)) {
                wedges[key] = {centre, candidate.second};
                return true;
            }
        }
        return false;
    };
    for (int centre = 0; centre < graph.n; ++centre) {
        if (graph.tetraChirality[centre] == 0) continue;
        std::vector<bool> visited(graph.n, false);
        if (!assignWedge(centre, visited))
            throw std::invalid_argument("Cannot represent tetrahedral stereo with the generated MOL coordinates");
    }
    return wedges;
}

enum class CoordinateFormat { FullPrecision, V2000, V3000 };

/** MOL E/Z is conveyed by coordinates, never by the bond stereo field. */
inline std::vector<Point2D> stereoCoordinates(const MolGraph& graph, CoordinateFormat format) {
    std::vector<std::pair<int, int>> stereoBonds;
    for (int atom = 0; atom < graph.n; ++atom) {
        for (int other : graph.neighbors[atom])
            if (other > atom && graph.bondOrder(atom, other) == 2 && graph.dbStereo(atom, other) != 0)
                stereoBonds.push_back({atom, other});
    }
    const bool hasTetrahedralStereo = std::any_of(graph.tetraChirality.begin(), graph.tetraChirality.end(),
                                                [](int winding) { return winding != 0; });
    if (stereoBonds.empty() && !hasTetrahedralStereo) return std::vector<Point2D>(graph.n, {0, 0});
    auto coordinates = generateCoords2D(graph);
    auto geometry = [&](int first, int second) {
        int firstLigand = -1, secondLigand = -1;
        for (int neighbor : graph.neighbors[first]) if (neighbor != second) { firstLigand = neighbor; break; }
        for (int neighbor : graph.neighbors[second]) if (neighbor != first) { secondLigand = neighbor; break; }
        if (firstLigand < 0 || secondLigand < 0) return 0;
        std::vector<MolAtom> atoms(graph.n);
        for (int atom = 0; atom < graph.n; ++atom) {
            atoms[atom].x = coordinates[atom].x;
            atoms[atom].y = coordinates[atom].y;
        }
        return doubleBondGeometry(atoms, first, second, firstLigand, secondLigand);
    };
    for (const auto& bond : stereoBonds) {
        const int first = bond.first, second = bond.second;
        if (geometry(first, second) == graph.dbStereo(first, second)) continue;
        // Reflect one component about the double-bond axis. This leaves all
        // other double-bond relations within that component unchanged.
        std::vector<bool> component(graph.n, false);
        std::vector<int> pending{second};
        component[second] = true;
        for (size_t position = 0; position < pending.size(); ++position) {
            const int atom = pending[position];
            for (int neighbor : graph.neighbors[atom]) {
                if ((atom == first && neighbor == second) || (atom == second && neighbor == first)) continue;
                if (!component[neighbor]) { component[neighbor] = true; pending.push_back(neighbor); }
            }
        }
        if (component[first]) {
            // Reflect all ligands at the endpoint together so their opposite
            // sides are retained when the cyclic double-bond relation changes.
            pending.clear();
            for (int neighbor : graph.neighbors[first])
                if (neighbor != second) pending.push_back(neighbor);
            if (pending.empty())
                throw std::invalid_argument("Double-bond stereo has no reference ligand");
        }
        const double axisX = coordinates[second].x - coordinates[first].x;
        const double axisY = coordinates[second].y - coordinates[first].y;
        const double lengthSquared = axisX * axisX + axisY * axisY;
        if (lengthSquared < 1e-10 || geometry(first, second) == 0)
            throw std::invalid_argument("Cannot represent double-bond stereo with degenerate MOL coordinates");
        for (int atom : pending) {
            const double relativeX = coordinates[atom].x - coordinates[first].x;
            const double relativeY = coordinates[atom].y - coordinates[first].y;
            const double projection = (relativeX * axisX + relativeY * axisY) / lengthSquared;
            coordinates[atom].x = coordinates[first].x + 2 * projection * axisX - relativeX;
            coordinates[atom].y = coordinates[first].y + 2 * projection * axisY - relativeY;
        }
    }
    auto ligandSide = [&](int first, int second, int ligand) {
        return (coordinates[second].x - coordinates[first].x) * (coordinates[ligand].y - coordinates[first].y)
             - (coordinates[second].y - coordinates[first].y) * (coordinates[ligand].x - coordinates[first].x);
    };
    for (const auto& bond : stereoBonds) {
        for (const auto& end : {bond, std::make_pair(bond.second, bond.first)}) {
            int reference = -1;
            for (int neighbor : graph.neighbors[end.first]) {
                if (neighbor == end.second) continue;
                if (reference < 0) { reference = neighbor; continue; }
                if (ligandSide(end.first, end.second, reference) * ligandSide(end.first, end.second, neighbor) < 0) continue;
                const double axisX = coordinates[end.second].x - coordinates[end.first].x;
                const double axisY = coordinates[end.second].y - coordinates[end.first].y;
                const double lengthSquared = axisX * axisX + axisY * axisY;
                if (lengthSquared < 1e-10)
                    throw std::invalid_argument("Cannot represent double-bond stereo with degenerate MOL coordinates");
                const double relativeX = coordinates[neighbor].x - coordinates[end.first].x;
                const double relativeY = coordinates[neighbor].y - coordinates[end.first].y;
                const double projection = (relativeX * axisX + relativeY * axisY) / lengthSquared;
                coordinates[neighbor].x = coordinates[end.first].x + 2 * projection * axisX - relativeX;
                coordinates[neighbor].y = coordinates[end.first].y + 2 * projection * axisY - relativeY;
            }
        }
    }
    if (format != CoordinateFormat::FullPrecision) {
        auto serialisedCoordinate = [&](double value) {
            if (!std::isfinite(value)) throw std::invalid_argument("MOL coordinates must be finite");
            double rounded = 0;
            if (format == CoordinateFormat::V2000) {
                char text[128];
                const int size = std::snprintf(text, sizeof(text), "%.4f", value);
                if (size < 0 || size >= static_cast<int>(sizeof(text)) || std::sscanf(text, "%lf", &rounded) != 1)
                    throw std::invalid_argument("Cannot represent a MOL coordinate");
            } else {
                std::ostringstream text;
                text << value;
                std::istringstream input(text.str());
                if (!(input >> rounded)) throw std::invalid_argument("Cannot represent a MOL coordinate");
            }
            return rounded;
        };
        for (auto& point : coordinates) {
            point.x = serialisedCoordinate(point.x);
            point.y = serialisedCoordinate(point.y);
        }
    }
    for (const auto& bond : stereoBonds) {
        if (geometry(bond.first, bond.second) != graph.dbStereo(bond.first, bond.second))
            throw std::invalid_argument("Cannot preserve double-bond stereo in MOL coordinates");
        for (const auto& end : {bond, std::make_pair(bond.second, bond.first)}) {
            int reference = -1;
            for (int neighbor : graph.neighbors[end.first]) {
                if (neighbor == end.second) continue;
                if (reference < 0) { reference = neighbor; continue; }
                if (ligandSide(end.first, end.second, reference) * ligandSide(end.first, end.second, neighbor) >= -1e-10)
                    throw std::invalid_argument("Cannot preserve opposite double-bond ligands in MOL coordinates");
            }
        }
    }
    return coordinates;
}

inline std::vector<Point2D> stereoCoordinates(const MolGraph& graph) {
    return stereoCoordinates(graph, CoordinateFormat::FullPrecision);
}

inline std::string stripV30Prefix(const std::string& line) {
    std::string trimmed = trim(line);
    const std::string prefix = "M  V30 ";
    if (trimmed.rfind(prefix, 0) == 0) return trimmed.substr(prefix.size());
    return trimmed;
}

inline std::map<std::string, std::string> parseV30Attributes(const std::vector<std::string>& tokens, size_t startIdx) {
    std::map<std::string, std::string> attrs;
    for (size_t i = startIdx; i < tokens.size(); ++i) {
        auto eq = tokens[i].find('=');
        if (eq == std::string::npos) continue;
        attrs.emplace(tokens[i].substr(0, eq), tokens[i].substr(eq + 1));
    }
    return attrs;
}

} // namespace mol_detail

inline MolGraph buildMolGraphFromParsedMol(
    const std::vector<mol_detail::MolAtom>& atoms,
    const std::vector<mol_detail::MolBond>& bonds,
    const std::string& molName,
    const std::string& programLine,
    const std::string& comment,
    std::map<std::string, std::string> properties) {
    int n = static_cast<int>(atoms.size());

    if (n == 0) {
        MolGraph g;
        g.n = 0;
        g.words = 0;
        g.name = molName;
        g.programLine = programLine;
        g.comment = comment;
        g.properties = std::move(properties);
        return g;
    }

    std::vector<std::vector<int>> neighbors(n);
    std::vector<std::vector<int>> bondOrders(n);
    std::vector<std::vector<int>> simpleAdj(n);
    std::vector<std::vector<int>> bondStereo(n);

    for (auto& b : bonds) {
        int i = b.from, j = b.to;
        if (i < 0 || i >= n || j < 0 || j >= n || i == j) continue;
        if (std::find(simpleAdj[i].begin(), simpleAdj[i].end(), j) != simpleAdj[i].end())
            continue;
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
        bondOrders[i].push_back(b.order);
        bondOrders[j].push_back(b.order);
        bondStereo[i].push_back(b.stereo);
        bondStereo[j].push_back(b.stereo);
        simpleAdj[i].push_back(j);
        simpleAdj[j].push_back(i);
    }

    std::vector<bool> atomInRing;
    std::vector<std::set<int>> bondInRingSet;
    detail::detectRings(n, simpleAdj, atomInRing, bondInRingSet);

    std::vector<std::vector<bool>> bondRingFlags(n);
    std::vector<std::vector<bool>> bondAromFlags(n);

    for (int i = 0; i < n; i++) {
        bondRingFlags[i].resize(neighbors[i].size(), false);
        bondAromFlags[i].resize(neighbors[i].size(), false);
        for (int k = 0; k < static_cast<int>(neighbors[i].size()); k++) {
            int j = neighbors[i][k];
            bondRingFlags[i][k] = bondInRingSet[i].count(j) > 0;
            bondAromFlags[i][k] = (bondOrders[i][k] == 4);
        }
    }

    std::vector<int> atomicNums(n), charges(n, 0), isotopes(n, 0), hydrogenCounts(n, 0), atomClasses(n, 0), tetraChirality(n, 0);
    std::vector<uint8_t> aromaticFlags(n, uint8_t(0)), ringFlags(n, uint8_t(0));

    for (int i = 0; i < n; i++) {
        atomicNums[i] = atoms[i].atomicNum;
        charges[i] = atoms[i].formalCharge;
        isotopes[i] = atoms[i].massNumber;
        int explicitValence = 0, explicitHydrogens = 0;
        bool aromaticAtom = false;
        for (size_t bondIndex = 0; bondIndex < neighbors[i].size(); ++bondIndex) {
            const int order = bondOrders[i][bondIndex];
            explicitValence += order == 4 ? 1 : order;
            aromaticAtom = aromaticAtom || order == 4;
            if (atoms[neighbors[i][bondIndex]].atomicNum == 1) ++explicitHydrogens;
        }
        if (atoms[i].hCount >= 0) hydrogenCounts[i] = std::max(0, atoms[i].hCount - explicitHydrogens);
        else if (atoms[i].valence >= 0) hydrogenCounts[i] = std::max(0, atoms[i].valence - explicitValence);
        else if (atomicNums[i] == 1 && charges[i] == 0 && atoms[i].radical == 0)
            hydrogenCounts[i] = std::max(0, 1 - explicitValence);
        else if (atomicNums[i] > 1)
            hydrogenCounts[i] = detail::computeImplicitH(atomicNums[i], aromaticAtom, explicitValence, charges[i]);
        atomClasses[i] = atoms[i].atomClass;
        tetraChirality[i] = atoms[i].tetraCfg;
        ringFlags[i] = atomInRing[i] ? 1 : 0;
        for (int k = 0; k < static_cast<int>(bondOrders[i].size()); k++) {
            if (bondOrders[i][k] == 4) {
                aromaticFlags[i] = 1;
                break;
            }
        }
    }

    std::vector<std::vector<const mol_detail::MolBond*>> wedgeBonds(n);
    for (const auto& bond : bonds)
        if (bond.order == 1 && (bond.stereo == 1 || bond.stereo == 6)) wedgeBonds[bond.from].push_back(&bond);
    for (int atom = 0; atom < n; ++atom) {
        tetraChirality[atom] = mol_detail::tetrahedralParity(tetraChirality[atom], atom, neighbors[atom], atomicNums);
        if (tetraChirality[atom] != 0 || wedgeBonds[atom].empty()) continue;
        auto positions = atoms;
        for (const auto* bond : wedgeBonds[atom])
            positions[bond->to].z = atoms[atom].z + (bond->stereo == 1 ? 1 : -1);
        tetraChirality[atom] = mol_detail::coordinateChirality(atom, neighbors[atom], positions);
    }

    MolGraph::Builder builder;
    builder.atomCount(n)
           .atomicNumbers(std::move(atomicNums))
           .formalCharges(std::move(charges))
           .massNumbers(std::move(isotopes))
           .hydrogenCounts(std::move(hydrogenCounts))
           .atomClasses(std::move(atomClasses))
           .ringFlags(std::move(ringFlags))
           .aromaticFlags(std::move(aromaticFlags))
           .setNeighbors(std::move(neighbors))
           .setBondOrders(std::move(bondOrders))
           .bondRingFlags(std::move(bondRingFlags))
           .bondAromaticFlags(std::move(bondAromFlags))
           .tetrahedralChirality(std::move(tetraChirality))
           .name(molName)
           .programLine(programLine)
           .comment(comment)
           .properties(std::move(properties));

    auto graph = builder.build();
    std::vector<int> ringDistance, ringQueue;
    auto hasSmallRing = [&](int first, int second) {
        if (bondInRingSet[first].count(second) == 0) return false;
        ringDistance.assign(n, -1);
        ringQueue.clear();
        ringDistance[first] = 0;
        ringQueue.push_back(first);
        for (size_t position = 0; position < ringQueue.size(); ++position) {
            const int atom = ringQueue[position];
            if (ringDistance[atom] >= 6) continue;
            for (int neighbor : graph.neighbors[atom]) {
                if ((atom == first && neighbor == second) || (atom == second && neighbor == first)) continue;
                if (ringDistance[neighbor] >= 0) continue;
                ringDistance[neighbor] = ringDistance[atom] + 1;
                if (neighbor == second) return true;
                ringQueue.push_back(neighbor);
            }
        }
        return false;
    };
    for (const auto& bond : bonds) {
        if (bond.order != 2 || bond.stereo == 3) continue;
        int firstLigand = -1, secondLigand = -1;
        for (int neighbor : graph.neighbors[bond.from]) if (neighbor != bond.to) { firstLigand = neighbor; break; }
        for (int neighbor : graph.neighbors[bond.to]) if (neighbor != bond.from) { secondLigand = neighbor; break; }
        if (firstLigand < 0 || secondLigand < 0) continue;
        // Small aromatic rings do not acquire E/Z tags from their drawing.
        if (graph.bondAromatic(bond.from, bond.to) && hasSmallRing(bond.from, bond.to)) continue;
        const int configuration = mol_detail::doubleBondGeometry(atoms, bond.from, bond.to, firstLigand, secondLigand);
        if (configuration == 0) continue;
        graph.hasDbStereo = true;
        if (graph.useDense) graph.dbStereoConf[bond.from][bond.to] = graph.dbStereoConf[bond.to][bond.from] = configuration;
        else graph.sparseDbStereo[MolGraph::bondKey(bond.from, bond.to)] = configuration;
    }

    graph.invalidateBondDerivedState();
    return graph;
}

inline MolGraph readMolBlockV3000(const std::vector<std::string>& lines) {
    if (lines.size() < 8) throw std::invalid_argument("V3000 MOL block too short");
    const std::string molName = lines.size() > 0 ? lines[0] : std::string();
    const std::string programLine = lines.size() > 1 ? lines[1] : std::string();
    const std::string comment = lines.size() > 2 ? lines[2] : std::string();

    std::vector<mol_detail::MolAtom> atoms;
    std::vector<mol_detail::MolBond> bonds;
    std::map<std::string, std::string> properties;
    std::unordered_map<int, int> atomIndexById;
    bool inAtom = false, inBond = false;
    int mEndLine = static_cast<int>(lines.size());

    for (int lineIdx = 4; lineIdx < static_cast<int>(lines.size()); ++lineIdx) {
        std::string body = mol_detail::stripV30Prefix(lines[lineIdx]);
        while (!body.empty() && body.back() == '-') {
            body.pop_back();
            if (++lineIdx >= static_cast<int>(lines.size()))
                throw std::invalid_argument("Truncated V3000 continuation line");
            body += " " + mol_detail::stripV30Prefix(lines[lineIdx]);
        }
        if (body == "BEGIN CTAB" || body == "COUNTS 0 0 0 0 0") continue;
        if (body == "END CTAB" || body == "END") { mEndLine = lineIdx; break; }
        if (body == "BEGIN ATOM") { inAtom = true; inBond = false; continue; }
        if (body == "END ATOM") { inAtom = false; continue; }
        if (body == "BEGIN BOND") { inBond = true; inAtom = false; continue; }
        if (body == "END BOND") { inBond = false; continue; }

        if (inAtom) {
            auto toks = mol_detail::splitWhitespace(body);
            if (toks.size() < 6) continue;
            mol_detail::MolAtom a;
            int atomId = std::stoi(toks[0]);
            a.symbol = toks[1];
            if (a.symbol == "R#" || (a.symbol.size() > 1 && a.symbol[0] == 'R' && std::isdigit(static_cast<unsigned char>(a.symbol[1])))) {
                a.atomicNum = 0;
                if (a.symbol != "R#") {
                    try { a.atomClass = std::stoi(a.symbol.substr(1)); } catch (...) {}
                }
            } else {
                a.atomicNum = mol_detail::elementToAtomicNum(a.symbol);
                if (a.atomicNum < 0) throw std::invalid_argument("Unknown V3000 element symbol: " + a.symbol);
                if (a.symbol == "D" || a.symbol == "T") {
                    a.massNumber = a.symbol == "D" ? 2 : 3;
                    a.symbol = "H";
                }
            }
            a.x = std::stod(toks[2]);
            a.y = std::stod(toks[3]);
            a.z = std::stod(toks[4]);
            auto attrs = mol_detail::parseV30Attributes(toks, 6);
            if (attrs.count("CHG")) a.formalCharge = std::stoi(attrs["CHG"]);
            if (attrs.count("MASS")) a.massNumber = std::stoi(attrs["MASS"]);
            if (attrs.count("RAD")) a.radical = std::stoi(attrs["RAD"]);
            if (attrs.count("HCOUNT")) {
                const int count = std::stoi(attrs["HCOUNT"]);
                if (count != 0) a.hCount = std::max(0, count);
            }
            if (attrs.count("VAL")) {
                const int valence = std::stoi(attrs["VAL"]);
                if (valence != 0) a.valence = std::max(0, valence);
            }
            if (attrs.count("CLASS")) a.atomClass = std::stoi(attrs["CLASS"]);
            if (attrs.count("CFG")) {
                int cfg = std::stoi(attrs["CFG"]);
                a.tetraCfg = (cfg == 1 ? 1 : (cfg == 2 ? 2 : 0));
            }
            atomIndexById.emplace(atomId, static_cast<int>(atoms.size()));
            atoms.push_back(std::move(a));
            continue;
        }

        if (inBond) {
            auto toks = mol_detail::splitWhitespace(body);
            if (toks.size() < 4) continue;
            mol_detail::MolBond b;
            const int bondType = std::stoi(toks[1]);
            if (bondType == 10) continue;
            b.order = bondType >= 1 && bondType <= 4 ? bondType : 1;
            int a1 = std::stoi(toks[2]);
            int a2 = std::stoi(toks[3]);
            auto it1 = atomIndexById.find(a1), it2 = atomIndexById.find(a2);
            if (it1 == atomIndexById.end() || it2 == atomIndexById.end()) {
                throw std::invalid_argument("V3000 bond references unknown atom");
            }
            b.from = it1->second;
            b.to = it2->second;
            auto attrs = mol_detail::parseV30Attributes(toks, 4);
            if (attrs.count("CFG")) {
                const int configuration = std::stoi(attrs["CFG"]);
                b.stereo = b.order == 2 ? (configuration == 2 ? 3 : 0)
                    : (configuration == 3 ? 6 : configuration);
            }
            bonds.push_back(std::move(b));
            continue;
        }
    }

    for (int lineIdx = mEndLine + 1; lineIdx < static_cast<int>(lines.size()); ++lineIdx) {
        const std::string& line = lines[lineIdx];
        if (line.size() >= 4 && line.substr(0, 4) == "$$$$") break;
        if (line.empty() || line[0] != '>') continue;
        auto lt = line.find('<');
        auto gt = line.find('>', lt == std::string::npos ? 0 : lt + 1);
        if (lt == std::string::npos || gt == std::string::npos || gt <= lt + 1) continue;
        std::string key = line.substr(lt + 1, gt - lt - 1);
        std::string value;
        bool first = true;
        int valueIdx = lineIdx + 1;
        for (; valueIdx < static_cast<int>(lines.size()); ++valueIdx) {
            const std::string& valueLine = lines[valueIdx];
            if (valueLine.empty()) break;
            if (!first) value.push_back('\n');
            value += valueLine;
            first = false;
        }
        properties[std::move(key)] = std::move(value);
        lineIdx = valueIdx;
    }

    return buildMolGraphFromParsedMol(atoms, bonds, molName, programLine, comment, std::move(properties));
}

// ============================================================================
// readMolBlock -- Read a single MOL block (V2000) from a string
// ============================================================================

inline MolGraph readMolBlock(const std::string& molBlock) {
    auto lines = mol_detail::splitLines(molBlock);

    if (lines.size() < 4) {
        throw std::invalid_argument("MOL block too short: need at least 4 lines (header + counts)");
    }

    // Lines 0-2: header (name, program/timestamp, comment)
    const std::string molName = lines.size() > 0 ? lines[0] : std::string();
    const std::string programLine = lines.size() > 1 ? lines[1] : std::string();
    const std::string comment = lines.size() > 2 ? lines[2] : std::string();
    // Line 3: counts line
    //   aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    //   aaa = atom count (cols 0-2), bbb = bond count (cols 3-5)
    const std::string& countsLine = lines[3];
    if (countsLine.size() < 6) {
        throw std::invalid_argument("Counts line too short");
    }

    int atomCount = mol_detail::safeInt(countsLine, 0, 3);
    int bondCount = mol_detail::safeInt(countsLine, 3, 3);

    if (atomCount < 0 || atomCount > 999) {
        throw std::invalid_argument("Invalid atom count: " + std::to_string(atomCount));
    }
    if (bondCount < 0 || bondCount > 999) {
        throw std::invalid_argument("Invalid bond count: " + std::to_string(bondCount));
    }

    // Check version: V2000 or V3000
    if (countsLine.size() >= 39) {
        std::string version = mol_detail::trim(countsLine.substr(33, 6));
        if (version == "V3000" || version == "v3000") {
            return readMolBlockV3000(lines);
        }
        if (!version.empty() && version != "V2000" && version != "v2000") {
            throw std::invalid_argument("Unsupported MOL version: " + version);
        }
    }

    int expectedLines = 4 + atomCount + bondCount;
    // Allow for M-lines and M END, but don't require exact count
    if (static_cast<int>(lines.size()) < 4 + atomCount + bondCount) {
        throw std::invalid_argument(
            "MOL block truncated: expected at least "
            + std::to_string(expectedLines) + " lines, got "
            + std::to_string(lines.size()));
    }

    // --- Parse atom block ---
    std::vector<mol_detail::MolAtom> atoms(atomCount);
    for (int i = 0; i < atomCount; i++) {
        int lineIdx = 4 + i;
        const std::string& line = lines[lineIdx];

        // Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // Coordinates: cols 0-9 (x), 10-19 (y), 20-29 (z)
        // Symbol: cols 31-33 (aaa)
        // Charge: cols 36-38 (ccc)

        atoms[i].x = mol_detail::safeDouble(line, 0, 10);
        atoms[i].y = mol_detail::safeDouble(line, 10, 10);
        atoms[i].z = mol_detail::safeDouble(line, 20, 10);

        // Atom symbol: columns 31-33, but could be shorter lines
        std::string sym;
        if (static_cast<int>(line.size()) > 31) {
            int symEnd = std::min(static_cast<int>(line.size()), 34);
            sym = mol_detail::trim(line.substr(31, symEnd - 31));
        }

        if (sym.empty()) {
            // Try alternative: some MOL files have the symbol right after the z coordinate
            // with just a space separator
            std::string tail;
            if (line.size() > 30) {
                tail = mol_detail::trim(line.substr(30));
            }
            if (!tail.empty()) {
                // Extract first word
                size_t sp = tail.find_first_of(" \t");
                sym = (sp != std::string::npos) ? tail.substr(0, sp) : tail;
            }
        }

        if (sym.empty()) {
            throw std::invalid_argument("Cannot parse atom symbol on line "
                                        + std::to_string(lineIdx + 1));
        }

        atoms[i].symbol = sym;
        if (sym == "R#" || (sym.size() > 1 && sym[0] == 'R' && std::isdigit(static_cast<unsigned char>(sym[1])))) {
            atoms[i].atomicNum = 0;
            if (sym != "R#") {
                try { atoms[i].atomClass = std::stoi(sym.substr(1)); } catch (...) {}
            }
        } else {
            atoms[i].atomicNum = mol_detail::elementToAtomicNum(sym);
            if (atoms[i].atomicNum < 0) {
                throw std::invalid_argument("Unknown element symbol: " + sym);
            }
            if (sym == "D" || sym == "T") {
                atoms[i].massNumber = sym == "D" ? 2 : 3;
                atoms[i].symbol = "H";
            }
        }

        // Charge from the atom-block ccc field (cols 36-38)
        int ccc = mol_detail::safeInt(line, 36, 3);
        atoms[i].formalCharge = mol_detail::v2000ChargeFromField(ccc);
        atoms[i].radical = ccc == 4 ? 2 : 0;
        int parity = mol_detail::safeInt(line, 39, 3);
        atoms[i].tetraCfg = (parity == 1 || parity == 2) ? parity : 0;
        const int hydrogenCount = mol_detail::safeInt(line, 42, 3);
        if (hydrogenCount > 0) atoms[i].hCount = hydrogenCount - 1;
        const int valence = mol_detail::safeInt(line, 48, 3);
        if (valence > 0) atoms[i].valence = valence == 15 ? 0 : valence;
    }

    // --- Parse bond block ---
    std::vector<mol_detail::MolBond> bonds(bondCount);
    for (int i = 0; i < bondCount; i++) {
        int lineIdx = 4 + atomCount + i;
        const std::string& line = lines[lineIdx];

        // Format: 111222tttsssxxxrrrccc
        // 111 = first atom (1-based), 222 = second atom (1-based), ttt = bond type
        int a1 = mol_detail::safeInt(line, 0, 3);
        int a2 = mol_detail::safeInt(line, 3, 3);
        int bt = mol_detail::safeInt(line, 6, 3, 1);

        // Validate
        if (a1 < 1 || a1 > atomCount || a2 < 1 || a2 > atomCount) {
            throw std::invalid_argument(
                "Bond references invalid atom index on line "
                + std::to_string(lineIdx + 1) + ": " + std::to_string(a1)
                + "-" + std::to_string(a2));
        }
        if (a1 == a2) {
            throw std::invalid_argument(
                "Self-loop bond on line " + std::to_string(lineIdx + 1));
        }

        bonds[i].from = a1 - 1;   // convert to 0-based
        bonds[i].to = a2 - 1;

        // Bond type: 1=single, 2=double, 3=triple, 4=aromatic, 5-8=query bonds
        if (bt >= 1 && bt <= 4) {
            bonds[i].order = bt;
        } else {
            // Treat unknown bond types as single
            bonds[i].order = 1;
        }
        int stereo = mol_detail::safeInt(line, 9, 3);
        bonds[i].stereo = bonds[i].order == 2 ? (stereo == 3 ? 3 : 0) : stereo;
    }

    // --- Parse M-lines (M CHG, M ISO, etc.) ---
    int propStart = 4 + atomCount + bondCount;
    int mEndLine = propStart;
    bool radicalProperties = false;
    std::map<std::string, std::string> properties;
    for (int lineIdx = propStart; lineIdx < static_cast<int>(lines.size()); lineIdx++) {
        const std::string& line = lines[lineIdx];

        if (line.size() >= 6 && line.substr(0, 6) == "M  END") {
            mEndLine = lineIdx;
            break;
        }
        if (line.size() >= 4 && line.substr(0, 4) == "$$$$") {
            mEndLine = lineIdx;
            break;
        }

        if (line.size() >= 1 && line[0] == 'A') {
            int atomIdx = mol_detail::safeInt(line, 3, 4) - 1;
            if (atomIdx >= 0 && atomIdx < atomCount && lineIdx + 1 < static_cast<int>(lines.size())) {
                std::string alias = mol_detail::trim(lines[++lineIdx]);
                if (alias == "R#" || (alias.size() > 1 && alias[0] == 'R'
                        && std::all_of(alias.begin() + 1, alias.end(), [](unsigned char c) { return std::isdigit(c); }))) {
                    atoms[atomIdx].atomicNum = 0;
                    if (alias != "R#") {
                        try { atoms[atomIdx].atomClass = std::stoi(alias.substr(1)); } catch (...) {}
                    }
                } else if (!alias.empty()) {
                    properties["ATOM_ALIAS_" + std::to_string(atomIdx + 1)] = alias;
                }
            }
            continue;
        }

        // M  CHGnn8 aaa vvv aaa vvv ...  (I4 fields start at column 9)
        if (!radicalProperties && line.size() >= 6
            && (line.substr(0, 6) == "M  CHG" || line.substr(0, 6) == "M  RAD")) {
            for (auto& atom : atoms) atom.radical = 0;
            radicalProperties = true;
        }
        if (line.size() >= 6 && line.substr(0, 6) == "M  CHG") {
            int count = mol_detail::safeInt(line, 6, 3);
            for (int j = 0; j < count; j++) {
                int pos = 9 + j * 8;
                int atomIdx = mol_detail::safeInt(line, pos, 4) - 1;  // 1-based to 0-based
                int charge = mol_detail::safeInt(line, pos + 4, 4);
                if (atomIdx >= 0 && atomIdx < atomCount) {
                    atoms[atomIdx].formalCharge = charge;
                }
            }
        }

        if (line.size() >= 6 && line.substr(0, 6) == "M  RAD") {
            const int count = mol_detail::safeInt(line, 6, 3);
            for (int entry = 0; entry < count; ++entry) {
                const int position = 9 + entry * 8;
                const int atomIndex = mol_detail::safeInt(line, position, 4) - 1;
                if (atomIndex >= 0 && atomIndex < atomCount)
                    atoms[atomIndex].radical = mol_detail::safeInt(line, position + 4, 4);
            }
        }

        // M  ISOnn8 aaa vvv aaa vvv ...  (I4 fields start at column 9)
        if (line.size() >= 6 && line.substr(0, 6) == "M  ISO") {
            int count = mol_detail::safeInt(line, 6, 3);
            for (int j = 0; j < count; j++) {
                int pos = 9 + j * 8;
                int atomIdx = mol_detail::safeInt(line, pos, 4) - 1;
                int isotope = mol_detail::safeInt(line, pos + 4, 4);
                if (atomIdx >= 0 && atomIdx < atomCount) {
                    atoms[atomIdx].massNumber = isotope;
                }
            }
        }

        if (line.size() >= 6 && line.substr(0, 6) == "M  RGP") {
            int count = mol_detail::safeInt(line, 6, 3);
            for (int j = 0; j < count; j++) {
                int pos = 9 + j * 8;
                int atomIdx = mol_detail::safeInt(line, pos, 4) - 1;
                int rg = mol_detail::safeInt(line, pos + 4, 4);
                if (atomIdx >= 0 && atomIdx < atomCount) {
                    atoms[atomIdx].atomicNum = 0;
                    atoms[atomIdx].atomClass = rg;
                }
            }
        }
    }

    for (int lineIdx = mEndLine + 1; lineIdx < static_cast<int>(lines.size()); ++lineIdx) {
        const std::string& line = lines[lineIdx];
        if (line.size() >= 4 && line.substr(0, 4) == "$$$$") break;
        if (line.empty() || line[0] != '>') continue;

        auto lt = line.find('<');
        auto gt = line.find('>', lt == std::string::npos ? 0 : lt + 1);
        if (lt == std::string::npos || gt == std::string::npos || gt <= lt + 1) continue;

        std::string key = line.substr(lt + 1, gt - lt - 1);
        std::string value;
        bool first = true;
        int valueIdx = lineIdx + 1;
        for (; valueIdx < static_cast<int>(lines.size()); ++valueIdx) {
            const std::string& valueLine = lines[valueIdx];
            if (valueLine.empty()) break;
            if (!first) value.push_back('\n');
            value += valueLine;
            first = false;
        }
        properties[std::move(key)] = std::move(value);
        lineIdx = valueIdx;
    }

    // --- Build MolGraph ---
    return buildMolGraphFromParsedMol(atoms, bonds, molName, programLine, comment, std::move(properties));
}

// ============================================================================
// readMolFile -- Read a MOL file from disk
// ============================================================================

inline MolGraph readMolFile(const std::string& filename) {
    std::ifstream ifs(std::filesystem::u8path(filename), std::ios::in);
    if (!ifs.is_open()) {
        throw std::invalid_argument("Cannot open MOL file: " + filename);
    }
    std::ostringstream oss;
    oss << ifs.rdbuf();
    return readMolBlock(oss.str());
}

// ============================================================================
// readSDF -- Read all molecules from an SDF file
// ============================================================================

inline std::vector<MolGraph> readSDF(const std::string& filename) {
    std::ifstream ifs(std::filesystem::u8path(filename), std::ios::in);
    if (!ifs.is_open()) {
        throw std::invalid_argument("Cannot open SDF file: " + filename);
    }

    std::vector<MolGraph> molecules;
    std::string line;
    std::ostringstream currentBlock;
    bool inBlock = false;

    while (std::getline(ifs, line)) {
        // Strip trailing \r
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.size() >= 4 && line.substr(0, 4) == "$$$$") {
            // End of molecule record
            if (inBlock) {
                std::string block = currentBlock.str();
                if (block.find_first_not_of(" \t\r\n") != std::string::npos) {
                    try {
                        molecules.push_back(readMolBlock(block));
                    } catch (const std::exception&) {
                        // Skip malformed molecules in SDF; push empty graph
                        molecules.emplace_back();
                    }
                }
            }
            currentBlock.str("");
            currentBlock.clear();
            inBlock = false;
        } else {
            if (!inBlock) inBlock = true;
            currentBlock << line << "\n";
        }
    }

    // Handle last block if file doesn't end with $$$$
    if (inBlock) {
        std::string block = currentBlock.str();
        if (block.find_first_not_of(" \t\r\n") != std::string::npos) {
            try {
                molecules.push_back(readMolBlock(block));
            } catch (const std::exception&) {
                molecules.emplace_back();
            }
        }
    }

    return molecules;
}

// ============================================================================
// writeMolBlock -- Write a MolGraph as MOL V2000 block
// ============================================================================

inline std::string writeMolBlockV3000(const MolGraph& g);

inline std::string writeMolBlock(const MolGraph& g) {
    std::ostringstream out;
    char buf[256];

    int n = g.n;
    size_t bondCount = 0;
    for (const auto& neighbors : g.neighbors) bondCount += neighbors.size();
    bondCount /= 2;
    if (n > 999 || bondCount > 999) return writeMolBlockV3000(g);
    const auto encoded = mol_detail::moleculeForWrite(g);
    const auto coordinates = mol_detail::stereoCoordinates(encoded, mol_detail::CoordinateFormat::V2000);
    const auto wedges = mol_detail::tetrahedralWedges(g, coordinates);

    // Each edge is stored once in the output.
    struct BondInfo { int from, to, order, stereo; };
    std::vector<BondInfo> bondList;
    bondList.reserve(bondCount);
    for (int i = 0; i < n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                int ord = encoded.bondOrder(i, j);
                const auto wedge = wedges.find(MolGraph::bondKey(i, j));
                if (wedge != wedges.end()) {
                    const int centre = wedge->second.first;
                    bondList.push_back({centre, centre == i ? j : i, ord, wedge->second.second});
                } else {
                    const bool unspecified = ord == 2 && g.dbStereo(i, j) == 0;
                    bondList.push_back({i, j, ord, unspecified ? 3 : 0});
                }
            }
        }
    }

    // Header lines
    out << g.name << "\n";
    out << (g.programLine.empty() ? "  SMSD    " : g.programLine) << "\n";
    out << g.comment << "\n";

    // Counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    std::snprintf(buf, sizeof(buf), "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
                  n, static_cast<int>(bondCount));
    out << buf;

    // Atom block
    for (int i = 0; i < n; i++) {
        const char* sym = (g.atomicNum[i] == 0 && i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] > 0)
            ? "R#" : detail::atomicNumToSymbol(g.atomicNum[i]);
        int ccc = mol_detail::v2000ChargeToField(g.formalCharge[i]);
        int parity = mol_detail::tetrahedralParity(g.tetraChirality[i], i, g.neighbors[i], g.atomicNum);
        const int valence = mol_detail::valenceForWrite(encoded, i);
        // Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        std::snprintf(buf, sizeof(buf),
                      "%10.4f%10.4f    0.0000 %-3s 0%3d%3d  0  0%3d  0  0  0  0  0  0\n",
                      coordinates[i].x, coordinates[i].y, sym, ccc, parity, valence < 0 ? 0 : valence == 0 ? 15 : valence);
        out << buf;
    }

    // Bond block
    for (auto& b : bondList) {
        std::snprintf(buf, sizeof(buf), "%3d%3d%3d%3d  0  0  0\n",
                      b.from + 1, b.to + 1, b.order, b.stereo);
        out << buf;
    }

    // M CHG lines carry every non-zero charge, as many MOL readers ignore the
    // atom-block ccc field and read charges from M CHG lines alone.
    {
        std::vector<std::pair<int, int>> chargeEntries;
        for (int i = 0; i < n; i++) {
            int ch = g.formalCharge[i];
            if (ch != 0) {
                chargeEntries.push_back({i + 1, ch});
            }
        }
        // Write in batches of 8 (MOL V2000 limit per M CHG line)
        for (size_t pos = 0; pos < chargeEntries.size(); pos += 8) {
            size_t batch = std::min(chargeEntries.size() - pos, size_t(8));
            std::snprintf(buf, sizeof(buf), "M  CHG%3d", static_cast<int>(batch));
            out << buf;
            for (size_t j = 0; j < batch; j++) {
                std::snprintf(buf, sizeof(buf), " %3d%4d",
                              chargeEntries[pos + j].first,
                              chargeEntries[pos + j].second);
                out << buf;
            }
            out << "\n";
        }
    }

    // M ISO lines for isotopes
    {
        std::vector<std::pair<int, int>> isoEntries;
        for (int i = 0; i < n; i++) {
            if (g.massNumber[i] != 0) {
                isoEntries.push_back({i + 1, g.massNumber[i]});
            }
        }
        for (size_t pos = 0; pos < isoEntries.size(); pos += 8) {
            size_t batch = std::min(isoEntries.size() - pos, size_t(8));
            std::snprintf(buf, sizeof(buf), "M  ISO%3d", static_cast<int>(batch));
            out << buf;
            for (size_t j = 0; j < batch; j++) {
                std::snprintf(buf, sizeof(buf), " %3d %3d",
                              isoEntries[pos + j].first,
                              isoEntries[pos + j].second);
                out << buf;
            }
            out << "\n";
        }
    }

    // M RGP lines for patent-style R-group labels on R# pseudo-atoms
    {
        std::vector<std::pair<int, int>> rgEntries;
        for (int i = 0; i < n; i++) {
            if (g.atomicNum[i] == 0 && i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] > 0) {
                rgEntries.push_back({i + 1, g.atomClass[i]});
            }
        }
        for (size_t pos = 0; pos < rgEntries.size(); pos += 8) {
            size_t batch = std::min(rgEntries.size() - pos, size_t(8));
            std::snprintf(buf, sizeof(buf), "M  RGP%3d", static_cast<int>(batch));
            out << buf;
            for (size_t j = 0; j < batch; j++) {
                std::snprintf(buf, sizeof(buf), " %3d%4d",
                              rgEntries[pos + j].first,
                              rgEntries[pos + j].second);
                out << buf;
            }
            out << "\n";
        }
    }

    out << "M  END\n";

    for (const auto& kv : g.properties) {
        out << "> <" << kv.first << ">\n";
        out << kv.second << "\n\n";
    }

    return out.str();
}

inline std::string writeMolBlockV3000(const MolGraph& g) {
    const auto encoded = mol_detail::moleculeForWrite(g);
    const auto coordinates = mol_detail::stereoCoordinates(encoded, mol_detail::CoordinateFormat::V3000);
    const auto wedges = mol_detail::tetrahedralWedges(g, coordinates);
    std::ostringstream out;
    out << g.name << "\n";
    out << (g.programLine.empty() ? "  SMSD    " : g.programLine) << "\n";
    out << g.comment << "\n";

    int bondCount = 0;
    struct BondInfo { int from, to, order, stereo; };
    std::vector<BondInfo> bondList;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                int ord = encoded.bondOrder(i, j);
                const auto wedge = wedges.find(MolGraph::bondKey(i, j));
                if (wedge != wedges.end()) {
                    const int centre = wedge->second.first;
                    bondList.push_back({centre, centre == i ? j : i, ord, wedge->second.second});
                } else {
                    const bool unspecified = ord == 2 && g.dbStereo(i, j) == 0;
                    bondList.push_back({i, j, ord, unspecified ? 3 : 0});
                }
                bondCount++;
            }
        }
    }

    out << "  0  0  0  0  0  0            999 V3000\n";
    out << "M  V30 BEGIN CTAB\n";
    out << "M  V30 COUNTS " << g.n << " " << bondCount << " 0 0 0\n";
    out << "M  V30 BEGIN ATOM\n";
    for (int i = 0; i < g.n; ++i) {
        const char* sym = (g.atomicNum[i] == 0 && i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] > 0)
            ? "R#" : detail::atomicNumToSymbol(g.atomicNum[i]);
        out << "M  V30 " << (i + 1) << " " << sym
            << " " << coordinates[i].x << " " << coordinates[i].y << " 0 0";
        if (i < static_cast<int>(g.formalCharge.size()) && g.formalCharge[i] != 0)
            out << " CHG=" << g.formalCharge[i];
        if (i < static_cast<int>(g.massNumber.size()) && g.massNumber[i] != 0)
            out << " MASS=" << g.massNumber[i];
        const int valence = mol_detail::valenceForWrite(encoded, i);
        if (valence >= 0) out << " VAL=" << (valence == 0 ? -1 : valence);
        if (i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] != 0)
            out << " CLASS=" << g.atomClass[i];
        if (i < static_cast<int>(g.tetraChirality.size()) && g.tetraChirality[i] != 0)
            out << " CFG=" << mol_detail::tetrahedralParity(g.tetraChirality[i], i, g.neighbors[i], g.atomicNum);
        out << "\n";
    }
    out << "M  V30 END ATOM\n";
    out << "M  V30 BEGIN BOND\n";
    for (size_t idx = 0; idx < bondList.size(); ++idx) {
        const auto& b = bondList[idx];
        out << "M  V30 " << (idx + 1) << " " << b.order << " " << (b.from + 1) << " " << (b.to + 1);
        if (b.stereo != 0) out << " CFG=" << (b.order == 2 ? 2 : (b.stereo == 6 ? 3 : 1));
        out << "\n";
    }
    out << "M  V30 END BOND\n";
    out << "M  V30 END CTAB\n";
    out << "M  END\n";
    for (const auto& kv : g.properties) {
        out << "> <" << kv.first << ">\n";
        out << kv.second << "\n\n";
    }
    return out.str();
}

inline std::string writeSDFRecord(const MolGraph& g) {
    std::string out = writeMolBlock(g);
    out += "$$$$\n";
    return out;
}

inline void writeSDF(const std::vector<MolGraph>& molecules, const std::string& filename) {
    std::ofstream ofs(std::filesystem::u8path(filename), std::ios::out);
    if (!ofs.is_open()) {
        throw std::invalid_argument("Cannot open output SDF file: " + filename);
    }
    for (const auto& mol : molecules) {
        ofs << writeSDFRecord(mol);
    }
}

} // namespace smsd

#endif // SMSD_MOL_READER_HPP

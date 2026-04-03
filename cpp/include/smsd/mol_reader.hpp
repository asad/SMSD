/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only MDL MOL V2000 / SDF file reader and writer for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 */
#pragma once
#ifndef SMSD_MOL_READER_HPP
#define SMSD_MOL_READER_HPP

#include "mol_graph.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
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
// Element table: symbol -> atomic number (shared with smiles_parser)
// --------------------------------------------------------------------------
inline int elementToAtomicNum(const std::string& sym) {
    static const std::unordered_map<std::string, int> table = {
        {"*",0},
        {"H",1},{"He",2},
        {"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
        {"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},
        {"K",19},{"Ca",20},
        {"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},
        {"Cu",29},{"Zn",30},{"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},
        {"Rb",37},{"Sr",38},
        {"Y",39},{"Zr",40},{"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},
        {"Ag",47},{"Cd",48},{"In",49},{"Sn",50},{"Sb",51},{"Te",52},{"I",53},{"Xe",54},
        {"Cs",55},{"Ba",56},
        {"La",57},{"Ce",58},{"Pr",59},{"Nd",60},{"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},
        {"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},{"Lu",71},
        {"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},
        {"Au",79},{"Hg",80},{"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},
        {"Fr",87},{"Ra",88},
        {"Ac",89},{"Th",90},{"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Am",95},{"Cm",96},
        {"Bk",97},{"Cf",98},{"Es",99},{"Fm",100},{"Md",101},{"No",102},{"Lr",103},
        {"Rf",104},{"Db",105},{"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},
        {"Rg",111},{"Cn",112},{"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118},
        // Common aliases
        {"D",1},{"T",1},
    };
    auto it = table.find(sym);
    return it != table.end() ? it->second : -1;
}

inline const char* atomicNumToSymbol(int z) {
    static const char* const syms[] = {
        "*","H","He",
        "Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar",
        "K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr",
        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
        "Ag","Cd","In","Sn","Sb","Te","I","Xe",
        "Cs","Ba",
        "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
        "Fr","Ra",
        "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
        "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    };
    if (z < 0 || z > 118) return "?";
    return syms[z];
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
// Ring detection (same algorithm as smiles_parser)
// --------------------------------------------------------------------------
inline void detectRings(int n,
                        const std::vector<std::vector<int>>& adj,
                        std::vector<bool>& atomInRing,
                        std::vector<std::set<int>>& bondInRing) {
    atomInRing.assign(n, false);
    bondInRing.resize(n);
    for (auto& s : bondInRing) s.clear();

    if (n == 0) return;

    std::vector<int> parent(n, -1);
    std::vector<int> depth(n, -1);
    std::vector<bool> visited(n, false);

    for (int start = 0; start < n; start++) {
        if (visited[start]) continue;

        struct Frame {
            int node;
            int parentNode;
            int neighborIdx;
        };
        std::vector<Frame> stack;
        stack.push_back({start, -1, 0});
        visited[start] = true;
        depth[start] = 0;
        parent[start] = -1;

        while (!stack.empty()) {
            auto& f = stack.back();
            if (f.neighborIdx >= static_cast<int>(adj[f.node].size())) {
                stack.pop_back();
                continue;
            }

            int nb = adj[f.node][f.neighborIdx];
            f.neighborIdx++;

            if (nb == f.parentNode) continue;

            if (visited[nb]) {
                if (depth[nb] < depth[f.node]) {
                    int cur = f.node;
                    while (cur != nb) {
                        atomInRing[cur] = true;
                        int p = parent[cur];
                        if (p >= 0) {
                            bondInRing[cur].insert(p);
                            bondInRing[p].insert(cur);
                        }
                        cur = p;
                    }
                    atomInRing[nb] = true;
                    bondInRing[f.node].insert(nb);
                    bondInRing[nb].insert(f.node);
                }
            } else {
                visited[nb] = true;
                depth[nb] = depth[f.node] + 1;
                parent[nb] = f.node;
                stack.push_back({nb, f.node, 0});
            }
        }
    }
}

// --------------------------------------------------------------------------
// V2000 charge mapping from the counts-line ccc field
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
    int hCount = 0;
    int atomClass = 0;
    int tetraCfg = 0;
};

struct MolBond {
    int from = 0;      // 0-based
    int to = 0;        // 0-based
    int order = 1;     // 1=single, 2=double, 3=triple, 4=aromatic
    int stereo = 0;
};

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
    mol_detail::detectRings(n, simpleAdj, atomInRing, bondInRingSet);

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
        hydrogenCounts[i] = atoms[i].hCount;
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

    std::vector<std::vector<int>> dbStereoMat;
    bool hasAnyStereo = false;
    for (const auto& b : bonds) {
        if (b.stereo != 0) { hasAnyStereo = true; break; }
    }
    if (hasAnyStereo && n <= MolGraph::SPARSE_THRESHOLD) {
        dbStereoMat.assign(n, std::vector<int>(n, 0));
        for (const auto& b : bonds) {
            if (b.stereo == 0) continue;
            if (b.order == 2) {
                dbStereoMat[b.from][b.to] = b.stereo;
                dbStereoMat[b.to][b.from] = b.stereo;
            }
        }
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

    if (hasAnyStereo && !dbStereoMat.empty()) builder.doubleBondStereo(std::move(dbStereoMat));
    return builder.build();
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
            }
            a.x = std::stod(toks[2]);
            a.y = std::stod(toks[3]);
            a.z = std::stod(toks[4]);
            auto attrs = mol_detail::parseV30Attributes(toks, 6);
            if (attrs.count("CHG")) a.formalCharge = std::stoi(attrs["CHG"]);
            if (attrs.count("MASS")) a.massNumber = std::stoi(attrs["MASS"]);
            if (attrs.count("HCOUNT")) a.hCount = std::max(0, std::stoi(attrs["HCOUNT"]) - 1);
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
            b.order = std::stoi(toks[1]);
            int a1 = std::stoi(toks[2]);
            int a2 = std::stoi(toks[3]);
            auto it1 = atomIndexById.find(a1), it2 = atomIndexById.find(a2);
            if (it1 == atomIndexById.end() || it2 == atomIndexById.end()) {
                throw std::invalid_argument("V3000 bond references unknown atom");
            }
            b.from = it1->second;
            b.to = it2->second;
            auto attrs = mol_detail::parseV30Attributes(toks, 4);
            if (attrs.count("CFG")) b.stereo = std::stoi(attrs["CFG"]);
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
        }

        // Charge from counts-line field (cols 36-38)
        int ccc = mol_detail::safeInt(line, 36, 3);
        atoms[i].formalCharge = mol_detail::v2000ChargeFromField(ccc);
        int parity = mol_detail::safeInt(line, 39, 3);
        atoms[i].tetraCfg = (parity == 1 || parity == 2) ? parity : 0;
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
        bonds[i].stereo = stereo;
    }

    // --- Parse M-lines (M CHG, M ISO, etc.) ---
    int propStart = 4 + atomCount + bondCount;
    int mEndLine = propStart;
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
    std::ifstream ifs(filename, std::ios::in);
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
    std::ifstream ifs(filename, std::ios::in);
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
                if (!block.empty()) {
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
        if (!block.empty()) {
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

inline std::string writeMolBlock(const MolGraph& g) {
    std::ostringstream out;
    char buf[256];

    int n = g.n;

    // Count bonds (each edge stored once)
    int bondCount = 0;
    struct BondInfo { int from, to, order, stereo; };
    std::vector<BondInfo> bondList;
    for (int i = 0; i < n; i++) {
        for (int j : g.neighbors[i]) {
            if (j > i) {
                int ord = g.bondOrder(i, j);
                if (g.bondAromatic(i, j)) ord = 4;
                int stereo = 0;
                if (!g.dbStereoConf.empty() && i < static_cast<int>(g.dbStereoConf.size())
                    && j < static_cast<int>(g.dbStereoConf[i].size())) {
                    stereo = g.dbStereoConf[i][j];
                }
                bondList.push_back({i, j, ord, stereo});
                bondCount++;
            }
        }
    }

    // Header lines
    out << g.name << "\n";
    out << (g.programLine.empty() ? "  SMSD    " : g.programLine) << "\n";
    out << g.comment << "\n";

    // Counts line: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    std::snprintf(buf, sizeof(buf), "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n",
                  n, bondCount);
    out << buf;

    // Atom block
    for (int i = 0; i < n; i++) {
        const char* sym = (g.atomicNum[i] == 0 && i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] > 0)
            ? "R#" : mol_detail::atomicNumToSymbol(g.atomicNum[i]);
        int ccc = mol_detail::v2000ChargeToField(g.formalCharge[i]);
        int parity = (i < static_cast<int>(g.tetraChirality.size())
                      && (g.tetraChirality[i] == 1 || g.tetraChirality[i] == 2))
            ? g.tetraChirality[i] : 0;
        // Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // We write 0.0 for coordinates (no 2D/3D info in MolGraph)
        std::snprintf(buf, sizeof(buf),
                      "    0.0000    0.0000    0.0000 %-3s 0%3d%3d  0  0  0  0  0  0  0  0  0\n",
                      sym, ccc, parity);
        out << buf;
    }

    // Bond block
    for (auto& b : bondList) {
        std::snprintf(buf, sizeof(buf), "%3d%3d%3d%3d  0  0  0\n",
                      b.from + 1, b.to + 1, b.order, b.stereo);
        out << buf;
    }

    // M CHG lines: always export all non-zero charges for modern compatibility.
    // Many modern mol-file parsers ignore the legacy ccc field
    // and rely exclusively on M CHG lines.
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
                int ord = g.bondOrder(i, j);
                if (g.bondAromatic(i, j)) ord = 4;
                int stereo = 0;
                if (!g.dbStereoConf.empty() && i < static_cast<int>(g.dbStereoConf.size())
                    && j < static_cast<int>(g.dbStereoConf[i].size())) {
                    stereo = g.dbStereoConf[i][j];
                }
                bondList.push_back({i, j, ord, stereo});
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
            ? "R#" : mol_detail::atomicNumToSymbol(g.atomicNum[i]);
        out << "M  V30 " << (i + 1) << " " << sym
            << " 0 0 0 0";
        if (i < static_cast<int>(g.formalCharge.size()) && g.formalCharge[i] != 0)
            out << " CHG=" << g.formalCharge[i];
        if (i < static_cast<int>(g.massNumber.size()) && g.massNumber[i] != 0)
            out << " MASS=" << g.massNumber[i];
        if (i < static_cast<int>(g.hydrogenCount.size()) && g.hydrogenCount[i] > 0)
            out << " HCOUNT=" << (g.hydrogenCount[i] + 1);
        if (i < static_cast<int>(g.atomClass.size()) && g.atomClass[i] != 0)
            out << " CLASS=" << g.atomClass[i];
        if (i < static_cast<int>(g.tetraChirality.size()) && g.tetraChirality[i] != 0)
            out << " CFG=" << (g.tetraChirality[i] == 1 ? 1 : 2);
        out << "\n";
    }
    out << "M  V30 END ATOM\n";
    out << "M  V30 BEGIN BOND\n";
    for (size_t idx = 0; idx < bondList.size(); ++idx) {
        const auto& b = bondList[idx];
        out << "M  V30 " << (idx + 1) << " " << b.order << " " << (b.from + 1) << " " << (b.to + 1);
        if (b.stereo != 0) out << " CFG=" << b.stereo;
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
    std::ofstream ofs(filename, std::ios::out);
    if (!ofs.is_open()) {
        throw std::invalid_argument("Cannot open output SDF file: " + filename);
    }
    for (const auto& mol : molecules) {
        ofs << writeSDFRecord(mol);
    }
}

} // namespace smsd

// ============================================================================
// Basic tests (compile with -DSMSD_TEST_MOL_READER)
// ============================================================================
#ifdef SMSD_TEST_MOL_READER

#include <cassert>
#include <iostream>

namespace smsd_mol_test {

inline void testBasicMolBlock() {
    // Methane: 1 carbon, no bonds
    std::string molBlock =
        "methane\n"
        "  test\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 1);
    assert(g.atomicNum[0] == 6);
    assert(g.formalCharge[0] == 0);
    std::cout << "  [PASS] testBasicMolBlock (methane)\n";
}

inline void testEthanolMolBlock() {
    // Ethanol: C-C-O
    std::string molBlock =
        "ethanol\n"
        "  test\n"
        "\n"
        "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 3);
    assert(g.atomicNum[0] == 6);
    assert(g.atomicNum[1] == 6);
    assert(g.atomicNum[2] == 8);
    assert(g.bondOrder(0, 1) == 1);
    assert(g.bondOrder(1, 2) == 1);
    assert(g.degree[1] == 2);
    std::cout << "  [PASS] testEthanolMolBlock\n";
}

inline void testBenzeneMolBlock() {
    // Benzene with aromatic bonds (type 4)
    std::string molBlock =
        "benzene\n"
        "  test\n"
        "\n"
        "  6  6  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  4  0  0  0  0\n"
        "  2  3  4  0  0  0  0\n"
        "  3  4  4  0  0  0  0\n"
        "  4  5  4  0  0  0  0\n"
        "  5  6  4  0  0  0  0\n"
        "  6  1  4  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 6);
    for (int i = 0; i < 6; i++) {
        assert(g.atomicNum[i] == 6);
        assert(g.ring[i] == true);
        assert(g.aromatic[i] == true);
        assert(g.degree[i] == 2);
    }
    assert(g.bondAromatic(0, 1) == true);
    assert(g.bondInRing(0, 1) == true);
    std::cout << "  [PASS] testBenzeneMolBlock (aromatic)\n";
}

inline void testChargesAndIsotopes() {
    // Ammonium ion with M CHG, heavy water with M ISO
    std::string molBlock =
        "charged\n"
        "  test\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  CHG  1   1   1\n"
        "M  ISO  1   2  18\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 2);
    assert(g.formalCharge[0] == 1);   // N+
    assert(g.massNumber[1] == 18);    // O-18
    std::cout << "  [PASS] testChargesAndIsotopes\n";
}

inline void testCountsLineCharge() {
    // Test charge from counts-line ccc field
    std::string molBlock =
        "ccc charge\n"
        "  test\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.formalCharge[0] == 1);  // ccc=3 means +1
    std::cout << "  [PASS] testCountsLineCharge\n";
}

inline void testDoubleBonds() {
    // Formaldehyde: C=O
    std::string molBlock =
        "formaldehyde\n"
        "  test\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  2  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.bondOrder(0, 1) == 2);
    std::cout << "  [PASS] testDoubleBonds (formaldehyde)\n";
}

inline void testTwoLetterElements() {
    // ClBr molecule
    std::string molBlock =
        "ClBr\n"
        "  test\n"
        "\n"
        "  2  1  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.atomicNum[0] == 17);  // Cl
    assert(g.atomicNum[1] == 35);  // Br
    std::cout << "  [PASS] testTwoLetterElements (Cl, Br)\n";
}

inline void testRoundTrip() {
    // Build a simple molecule, write it, read it back
    std::string molBlock =
        "acetic acid\n"
        "  test\n"
        "\n"
        "  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.0000    0.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    3.0000   -0.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  2  0  0  0  0\n"
        "  2  4  1  0  0  0  0\n"
        "M  END\n";

    auto g1 = smsd::readMolBlock(molBlock);
    std::string written = smsd::writeMolBlock(g1);
    auto g2 = smsd::readMolBlock(written);

    assert(g1.n == g2.n);
    for (int i = 0; i < g1.n; i++) {
        assert(g1.atomicNum[i] == g2.atomicNum[i]);
        assert(g1.formalCharge[i] == g2.formalCharge[i]);
        assert(g1.degree[i] == g2.degree[i]);
    }
    assert(g1.bondOrder(1, 2) == g2.bondOrder(1, 2));
    std::cout << "  [PASS] testRoundTrip (acetic acid)\n";
}

inline void testEmptyMolBlock() {
    std::string molBlock =
        "empty\n"
        "  test\n"
        "\n"
        "  0  0  0  0  0  0  0  0  0  0999 V2000\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 0);
    std::cout << "  [PASS] testEmptyMolBlock\n";
}

inline void testMalformedBlock() {
    bool threw = false;
    try {
        smsd::readMolBlock("too short");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
    std::cout << "  [PASS] testMalformedBlock\n";
}

inline void testCyclohexane() {
    // Cyclohexane: 6 C in a ring with single bonds
    std::string molBlock =
        "cyclohexane\n"
        "  test\n"
        "\n"
        "  6  6  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0  0  0  0\n"
        "  2  3  1  0  0  0  0\n"
        "  3  4  1  0  0  0  0\n"
        "  4  5  1  0  0  0  0\n"
        "  5  6  1  0  0  0  0\n"
        "  6  1  1  0  0  0  0\n"
        "M  END\n";

    auto g = smsd::readMolBlock(molBlock);
    assert(g.n == 6);
    for (int i = 0; i < 6; i++) {
        assert(g.ring[i] == true);
        assert(g.aromatic[i] == false);
        assert(g.degree[i] == 2);
    }
    assert(g.bondInRing(0, 1) == true);
    std::cout << "  [PASS] testCyclohexane (ring perception)\n";
}

inline void runAllTests() {
    std::cout << "mol_reader tests:\n";
    testBasicMolBlock();
    testEthanolMolBlock();
    testBenzeneMolBlock();
    testChargesAndIsotopes();
    testCountsLineCharge();
    testDoubleBonds();
    testTwoLetterElements();
    testRoundTrip();
    testEmptyMolBlock();
    testMalformedBlock();
    testCyclohexane();
    std::cout << "All mol_reader tests passed.\n";
}

} // namespace smsd_mol_test

#endif // SMSD_TEST_MOL_READER

#endif // SMSD_MOL_READER_HPP

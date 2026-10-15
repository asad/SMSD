/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only SMARTS pattern parser for smsd::MolGraph.
 * Zero external dependencies -- pure C++17 standard library.
 *
 * Supports:
 *   - Atom primitives: #n, A, a, D<n>, v<n>, R, r<n>, +/-, H<n>, *, elements
 *   - Extended primitives: z<n> (heteroatom neighbors), Z<n> (aliphatic
 *     heteroatom neighbors), d<n> (heavy degree), ^<n> (hybridization)
 *   - Range queries: {n-m}, {n-}, {-m} on numeric primitives
 *   - Logical operators: , (OR), & (AND high), ; (AND low), ! (NOT)
 *   - Bond primitives: - (single), = (double), # (triple), : (aromatic),
 *                       ~ (any), @ (ring bond), / (wedge), \ (dash)
 *   - Ring closures, branches ()
 *   - Recursive SMARTS: $(...)
 *   - Named predicate registry (SmartsPredicateRegistry) with 10 built-ins
 */
#pragma once
#ifndef SMSD_SMARTS_PARSER_HPP
#define SMSD_SMARTS_PARSER_HPP

#include "mol_graph.hpp"
#include "parser_utils.hpp"
#include "stereo.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <deque>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// Forward declarations
// ============================================================================
struct SmartsQuery;
SmartsQuery parseSMARTS(const std::string& smarts, int maxRecursionDepth = 20);

// ============================================================================
// Element table for SMARTS (shared inline, same as smiles_parser)
// ============================================================================
namespace smarts_detail {

using RecursiveMatcher = bool (*)(const SmartsQuery&, const MolGraph&, int);
inline thread_local RecursiveMatcher recursiveMatcher = nullptr;

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
        // Aromatic lowercase forms
        {"b",5},{"c",6},{"n",7},{"o",8},{"p",15},{"s",16},{"se",34},{"te",52},
    };
    auto it = table.find(sym);
    return it != table.end() ? it->second : -1;
}

// ============================================================================
// Atom expression tree
// ============================================================================

enum class AtomPrimType {
    TRUE_,           // matches anything (wildcard)
    ATOMIC_NUM,      // #n
    ELEMENT,         // specific element (also encodes atomicNum)
    AROMATIC,        // a
    ALIPHATIC,       // A
    DEGREE,          // D<n>
    VALENCE,         // v<n>
    IN_RING,         // R (any ring) or R<n> (ring count)
    RING_SIZE,       // r<n>
    CHARGE,          // +n / -n
    HCOUNT,          // H<n>
    RING_CONNECTIVITY, // x<n>
    MASS,            // isotope
    RECURSIVE,       // $(...)
    HETERO_NEIGHBORS,          // z<n> — count of heteroatom neighbours (not C, not H)
    ALIPHATIC_HETERO_NEIGHBORS,// Z<n> — count of aliphatic heteroatom neighbours
    TOTAL_CONNECTIVITY,        // X<n> — total connections (degree + implicit H)
    HEAVY_DEGREE,              // d<n> — non-hydrogen degree (heavy atom connections)
    HYBRIDIZATION,             // ^<n> — hybridization (0=S,1=SP,2=SP2,3=SP3,4=SP3D,5=SP3D2)
    CHIRALITY,                 // @ (CW), @@ (CCW), @? (either/unspecified)
    ATOM_CLASS,                // :n — atom map number / reaction class
};

enum class ExprOp {
    PRIM,    // leaf
    AND,
    OR,
    NOT,
};

struct AtomExpr {
    ExprOp op = ExprOp::PRIM;

    // For PRIM:
    AtomPrimType primType = AtomPrimType::TRUE_;
    int intVal = 0;                // numeric value for the primitive
    bool hasValue = false;         // whether intVal is meaningful
    // Range query support: {min-max}  (rangeQuery == true overrides intVal)
    bool rangeQuery = false;       // true when {n-m} syntax was used
    int rangeMin = 0;              // minimum (inclusive), -1 = unbounded
    int rangeMax = 0;              // maximum (inclusive), -1 = unbounded
    // For RECURSIVE:
    std::string recursiveSmarts;
    std::shared_ptr<const SmartsQuery> cachedRecursive;

    // For AND/OR: children
    std::vector<AtomExpr> children;
    // For NOT: single child in children[0]

    static AtomExpr makePrim(AtomPrimType t, int val = 0, bool hasVal = false) {
        AtomExpr e;
        e.op = ExprOp::PRIM;
        e.primType = t;
        e.intVal = val;
        e.hasValue = hasVal;
        return e;
    }

    static AtomExpr makeNot(AtomExpr child) {
        AtomExpr e;
        e.op = ExprOp::NOT;
        e.children.push_back(std::move(child));
        return e;
    }

    static AtomExpr makeAnd(std::vector<AtomExpr> ch) {
        if (ch.size() == 1) return std::move(ch[0]);
        AtomExpr e;
        e.op = ExprOp::AND;
        e.children = std::move(ch);
        return e;
    }

    static AtomExpr makeOr(std::vector<AtomExpr> ch) {
        if (ch.size() == 1) return std::move(ch[0]);
        AtomExpr e;
        e.op = ExprOp::OR;
        e.children = std::move(ch);
        return e;
    }
};

// ============================================================================
// Bond expression tree
// ============================================================================

enum class BondPrimType {
    ANY,       // ~ or implicit (default)
    SINGLE,    // -
    DOUBLE,    // =
    TRIPLE,    // #
    AROMATIC,  // :
    RING,      // @
    WEDGE,     // / (up bond / E-Z stereo)
    DASH,      // \ (down bond / E-Z stereo)
    DEFAULT_,  // implicit bond (single or aromatic in aromatic context)
};

struct BondExpr {
    ExprOp op = ExprOp::PRIM;
    BondPrimType primType = BondPrimType::DEFAULT_;
    std::vector<BondExpr> children;

    static BondExpr makePrim(BondPrimType t) {
        BondExpr e;
        e.op = ExprOp::PRIM;
        e.primType = t;
        return e;
    }

    static BondExpr makeNot(BondExpr child) {
        BondExpr e;
        e.op = ExprOp::NOT;
        e.children.push_back(std::move(child));
        return e;
    }

    static BondExpr makeAnd(std::vector<BondExpr> ch) {
        if (ch.size() == 1) return std::move(ch[0]);
        BondExpr e;
        e.op = ExprOp::AND;
        e.children = std::move(ch);
        return e;
    }

    static BondExpr makeOr(std::vector<BondExpr> ch) {
        if (ch.size() == 1) return std::move(ch[0]);
        BondExpr e;
        e.op = ExprOp::OR;
        e.children = std::move(ch);
        return e;
    }
};

// ============================================================================
// SMARTS query atom and bond
// ============================================================================

struct SmartsAtom {
    AtomExpr expr;
    bool isAromatic = false;   // hint from lowercase symbol (for default bond)
};

struct SmartsBond {
    int from = -1;
    int to = -1;
    BondExpr expr;
};

// ============================================================================
// Target-side atom property cache
// ============================================================================

inline int computeValenceRaw(const MolGraph& g, int idx) {
    int val = 0;
    for (int nb : g.neighbors[idx]) {
        int bo = g.bondOrder(idx, nb);
        val += (bo == 4 ? 1 : bo);
    }
    return val;
}

inline int computeHCountFromValence(const MolGraph& g, int idx, int explicitValence) {
    if (idx < static_cast<int>(g.hydrogenCount.size()))
        return g.hydrogenCount[idx];

    // Default valences for common elements
    static const std::unordered_map<int, std::vector<int>> defaultValences = {
        {5,  {3}},       // B
        {6,  {4}},       // C
        {7,  {3}},       // N
        {8,  {2}},       // O
        {9,  {1}},       // F
        {14, {4}},       // Si
        {15, {3, 5}},    // P
        {16, {2, 4, 6}}, // S
        {17, {1}},       // Cl
        {32, {4}},       // Ge
        {33, {3, 5}},    // As
        {34, {2, 4, 6}}, // Se
        {35, {1}},       // Br
        {52, {2, 4, 6}}, // Te
        {53, {1, 3, 5, 7}}, // I
    };

    int z = g.atomicNum[idx];
    auto it = defaultValences.find(z);
    if (it == defaultValences.end()) return 0;

    const int64_t charge = g.formalCharge[idx] < 0
        ? -int64_t(g.formalCharge[idx]) : int64_t(g.formalCharge[idx]);

    // Find the smallest default valence >= explicitValence + charge
    for (int dv : it->second) {
        const int64_t hc = dv - int64_t(explicitValence) - charge;
        if (hc >= 0) return hc;
    }
    return 0;
}

inline int computeTotalValence(const MolGraph& g, int atom, int explicitValence,
                               int implicitHydrogens) {
    int valence = explicitValence + implicitHydrogens;
    if (!g.aromatic[atom]) return valence;

    int64_t normalValence;
    switch (g.atomicNum[atom]) {
        case 5: case 7: case 15: case 33: normalValence = 3; break;
        case 6: normalValence = 4; break;
        case 8: case 16: case 34: case 52: normalValence = 2; break;
        default: return valence;
    }
    const int charge = g.formalCharge[atom];
    const int atomicNum = g.atomicNum[atom];
    if ((charge > 0 && (atomicNum == 7 || atomicNum == 15 || atomicNum == 33
                       || atomicNum == 8 || atomicNum == 16 || atomicNum == 34))
        || (charge < 0 && atomicNum == 5)) {
        normalValence += charge < 0 ? -int64_t(charge) : int64_t(charge);
    } else {
        normalValence -= charge < 0 ? -int64_t(charge) : int64_t(charge);
    }
    // A Kekulé double bond adds one, except at an occupied lone-pair donor.
    if (valence < normalValence) ++valence;
    return valence;
}

inline int computeRingMembershipRaw(const MolGraph& g, int idx) {
    g.ensureRingCounts();
    if (idx < static_cast<int>(g.ringCount.size()))
        return g.ringCount[idx];
    return g.ring[idx] ? 1 : 0;
}

inline int computeSmallestRingSizeRaw(const MolGraph& g, int idx) {
    if (!g.ring[idx]) return 0;

    int best = 999;
    for (int nb : g.neighbors[idx]) {
        if (!g.bondInRing(idx, nb)) continue;

        std::vector<int> dist(g.n, -1);
        dist[nb] = 1;
        std::deque<int> queue;
        queue.push_back(nb);
        bool found = false;

        while (!queue.empty() && !found) {
            int u = queue.front();
            queue.pop_front();
            for (int v : g.neighbors[u]) {
                if (!g.bondInRing(u, v)) continue;
                if (u == nb && v == idx) continue;
                if (v == idx) {
                    int ringSize = dist[u] + 1;
                    if (ringSize < best) best = ringSize;
                    found = true;
                    break;
                }
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    queue.push_back(v);
                }
            }
        }
    }

    return best < 999 ? best : 0;
}

inline int computeRingConnectivityRaw(const MolGraph& g, int idx) {
    int rc = 0;
    for (int nb : g.neighbors[idx]) {
        if (g.bondInRing(idx, nb)) ++rc;
    }
    return rc;
}

inline int computeHeteroNeighborsRaw(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        int z = g.atomicNum[nb];
        if (z != 6 && z != 1) ++count;
    }
    return count;
}

inline int computeAliphaticHeteroNeighborsRaw(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        int z = g.atomicNum[nb];
        if (z != 6 && z != 1 && !g.aromatic[nb]) ++count;
    }
    return count;
}

inline int computeHeavyDegreeRaw(const MolGraph& g, int idx) {
    int count = 0;
    for (int nb : g.neighbors[idx]) {
        if (g.atomicNum[nb] != 1) ++count;
    }
    return count;
}

inline int computeHybridizationRaw(const MolGraph& g, int idx) {
    if (g.aromatic[idx]) return 2; // aromatic => SP2

    int nDouble = 0;
    int nTriple = 0;
    for (int nb : g.neighbors[idx]) {
        int bo = g.bondOrder(idx, nb);
        if (bo == 2) ++nDouble;
        else if (bo == 3) ++nTriple;
    }

    if (nTriple >= 1 || nDouble >= 2) return 1; // SP
    if (nDouble == 1) return 2;                  // SP2
    return 3;                                    // SP3
}

struct SmartsTargetPropertyCache {
    const MolGraph* graph = nullptr;
    uint64_t graphNonce = 0;
    int atomCount = -1;

    std::vector<int> valence;
    std::vector<int> hcount;
    std::vector<int> ringMembership;
    std::vector<int> smallestRingSize;
    std::vector<int> ringConnectivity;
    std::vector<int> heteroNeighbors;
    std::vector<int> aliphaticHeteroNeighbors;
    std::vector<int> heavyDegree;
    std::vector<int> hybridization;

    void rebuild(const MolGraph& g) {
        graph = &g;
        graphNonce = g.cacheNonce;
        atomCount = g.n;

        valence.assign(g.n, 0);
        hcount.assign(g.n, 0);
        ringMembership.assign(g.n, 0);
        smallestRingSize.assign(g.n, 0);
        ringConnectivity.assign(g.n, 0);
        heteroNeighbors.assign(g.n, 0);
        aliphaticHeteroNeighbors.assign(g.n, 0);
        heavyDegree.assign(g.n, 0);
        hybridization.assign(g.n, 0);

        g.ensureRingCounts();

        for (int i = 0; i < g.n; ++i) {
            valence[i] = computeValenceRaw(g, i);
            ringMembership[i] = computeRingMembershipRaw(g, i);
            ringConnectivity[i] = computeRingConnectivityRaw(g, i);
            heteroNeighbors[i] = computeHeteroNeighborsRaw(g, i);
            aliphaticHeteroNeighbors[i] = computeAliphaticHeteroNeighborsRaw(g, i);
            heavyDegree[i] = computeHeavyDegreeRaw(g, i);
            hybridization[i] = computeHybridizationRaw(g, i);
        }

        for (int i = 0; i < g.n; ++i) {
            const int implicitHydrogens = computeHCountFromValence(g, i, valence[i]);
            hcount[i] = implicitHydrogens;
            for (int neighbor : g.neighbors[i]) {
                if (g.atomicNum[neighbor] == 1) ++hcount[i];
            }
            valence[i] = computeTotalValence(g, i, valence[i], implicitHydrogens);
            smallestRingSize[i] = g.ring[i] ? computeSmallestRingSizeRaw(g, i) : 0;
        }
    }
};

inline SmartsTargetPropertyCache& getSmartsTargetPropertyCache(const MolGraph& g) {
    thread_local SmartsTargetPropertyCache cache;
    if (cache.graph != &g
        || cache.graphNonce != g.cacheNonce
        || cache.atomCount != g.n)
        cache.rebuild(g);
    return cache;
}

// ============================================================================
// Cached property accessors used by SMARTS matching
// ============================================================================

inline int computeValence(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).valence[idx];
}

inline int computeHCount(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).hcount[idx];
}

inline int computeRingMembership(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).ringMembership[idx];
}

inline int computeSmallestRingSize(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).smallestRingSize[idx];
}

inline int computeRingConnectivity(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).ringConnectivity[idx];
}

inline int computeHeteroNeighbors(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).heteroNeighbors[idx];
}

inline int computeAliphaticHeteroNeighbors(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).aliphaticHeteroNeighbors[idx];
}

inline int computeHeavyDegree(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).heavyDegree[idx];
}

inline int computeHybridization(const MolGraph& g, int idx) {
    return getSmartsTargetPropertyCache(g).hybridization[idx];
}

// ============================================================================
// Expression evaluation
// ============================================================================

// Forward declarations
inline bool evalAtomExpr(const AtomExpr& expr, const MolGraph& g, int idx);
inline bool evalRecursiveSmarts(const AtomExpr& expr, const MolGraph& g, int idx);

// Helper: compare a computed value against an exact value or range
inline bool matchNumeric(const AtomExpr& expr, int computed) {
    if (expr.rangeQuery) {
        bool ok = true;
        if (expr.rangeMin >= 0) ok = ok && (computed >= expr.rangeMin);
        if (expr.rangeMax >= 0) ok = ok && (computed <= expr.rangeMax);
        return ok;
    }
    return computed == expr.intVal;
}

inline bool evalAtomPrim(const AtomExpr& expr, const MolGraph& g, int idx) {
    switch (expr.primType) {
        case AtomPrimType::TRUE_:
            return true;

        case AtomPrimType::ATOMIC_NUM:
            return g.atomicNum[idx] == expr.intVal;

        case AtomPrimType::ELEMENT:
            return g.atomicNum[idx] == expr.intVal;

        case AtomPrimType::AROMATIC:
            return g.aromatic[idx];

        case AtomPrimType::ALIPHATIC:
            return !g.aromatic[idx];

        case AtomPrimType::DEGREE:
            return matchNumeric(expr, g.degree[idx]);

        case AtomPrimType::VALENCE:
            return matchNumeric(expr, computeValence(g, idx));

        case AtomPrimType::IN_RING:
            if (!expr.hasValue && !expr.rangeQuery) return g.ring[idx];
            if (!expr.rangeQuery && expr.intVal == 0) return !g.ring[idx];
            return matchNumeric(expr, computeRingMembership(g, idx));

        case AtomPrimType::RING_SIZE:
            if (!expr.hasValue && !expr.rangeQuery) return g.ring[idx];
            return matchNumeric(expr, computeSmallestRingSize(g, idx));

        case AtomPrimType::CHARGE:
            return g.formalCharge[idx] == expr.intVal;

        case AtomPrimType::HCOUNT:
            return matchNumeric(expr, computeHCount(g, idx));

        case AtomPrimType::RING_CONNECTIVITY:
            if (!expr.hasValue && !expr.rangeQuery) return computeRingConnectivity(g, idx) > 0;
            return matchNumeric(expr, computeRingConnectivity(g, idx));

        case AtomPrimType::MASS:
            return g.massNumber[idx] == expr.intVal;

        case AtomPrimType::RECURSIVE:
            // Deferred to evalRecursiveSmarts() defined after SmartsQuery
            return evalRecursiveSmarts(expr, g, idx);

        case AtomPrimType::HETERO_NEIGHBORS:
            if (!expr.hasValue && !expr.rangeQuery) return computeHeteroNeighbors(g, idx) >= 1;
            return matchNumeric(expr, computeHeteroNeighbors(g, idx));

        case AtomPrimType::ALIPHATIC_HETERO_NEIGHBORS:
            if (!expr.hasValue && !expr.rangeQuery) return computeAliphaticHeteroNeighbors(g, idx) >= 1;
            return matchNumeric(expr, computeAliphaticHeteroNeighbors(g, idx));

        case AtomPrimType::TOTAL_CONNECTIVITY:
            return matchNumeric(expr, computeHeavyDegree(g, idx) + computeHCount(g, idx));

        case AtomPrimType::HEAVY_DEGREE:
            return matchNumeric(expr, computeHeavyDegree(g, idx));

        case AtomPrimType::HYBRIDIZATION:
            return computeHybridization(g, idx) == expr.intVal;

        case AtomPrimType::CHIRALITY:
            // 0=unspecified/@? (match any), 1=CW/@, 2=CCW/@@
            if (expr.intVal == 0) return true; // @? matches any chirality
            if (g.tetraChirality.empty() || g.tetraChirality[idx] == 0) return false;
            return g.tetraChirality[idx] == expr.intVal;

        case AtomPrimType::ATOM_CLASS:
            // Atom class/map — always matches in SMARTS matching context
            // (atom classes are used for labelling, not filtering)
            return true;
    }
    return false;
}

inline bool evalAtomExpr(const AtomExpr& expr, const MolGraph& g, int idx) {
    switch (expr.op) {
        case ExprOp::PRIM:
            return evalAtomPrim(expr, g, idx);

        case ExprOp::AND:
            for (auto& child : expr.children) {
                if (!evalAtomExpr(child, g, idx)) return false;
            }
            return true;

        case ExprOp::OR:
            for (auto& child : expr.children) {
                if (evalAtomExpr(child, g, idx)) return true;
            }
            return false;

        case ExprOp::NOT:
            return !evalAtomExpr(expr.children[0], g, idx);
    }
    return false;
}

inline bool evalBondPrim(const BondExpr& expr, const MolGraph& g, int i, int j) {
    switch (expr.primType) {
        case BondPrimType::ANY:
            return true;

        case BondPrimType::SINGLE:
            return g.bondOrder(i, j) == 1 && !g.bondAromatic(i, j);

        case BondPrimType::DOUBLE:
            return g.bondOrder(i, j) == 2;

        case BondPrimType::TRIPLE:
            return g.bondOrder(i, j) == 3;

        case BondPrimType::AROMATIC:
            return g.bondAromatic(i, j) || g.bondOrder(i, j) == 4;

        case BondPrimType::RING:
            return g.bondInRing(i, j);

        case BondPrimType::WEDGE:
        case BondPrimType::DASH:
            return g.bondOrder(i, j) == 1 && !g.bondAromatic(i, j);

        case BondPrimType::DEFAULT_: {
            // Default SMARTS bond: matches single (1) or aromatic (4)
            int ord = g.bondOrder(i, j);
            return ord == 1 || ord == 4;
        }
    }
    return false;
}

inline bool evalBondExpr(const BondExpr& expr, const MolGraph& g, int i, int j) {
    switch (expr.op) {
        case ExprOp::PRIM:
            return evalBondPrim(expr, g, i, j);

        case ExprOp::AND:
            for (auto& child : expr.children) {
                if (!evalBondExpr(child, g, i, j)) return false;
            }
            return true;

        case ExprOp::OR:
            for (auto& child : expr.children) {
                if (evalBondExpr(child, g, i, j)) return true;
            }
            return false;

        case ExprOp::NOT:
            return !evalBondExpr(expr.children[0], g, i, j);
    }
    return false;
}

/** Preserve Boolean expressions while mapped chirality is not yet known. */
enum class StereoTruth { FALSE_VALUE, TRUE_VALUE, UNKNOWN };

template <typename ChiralityEvaluator>
inline StereoTruth evalAtomExprWithStereo(const AtomExpr& expression, const MolGraph& graph,
                                          int atom, const ChiralityEvaluator& evaluateChirality) {
    if (expression.op == ExprOp::PRIM) {
        if (expression.primType == AtomPrimType::CHIRALITY) return evaluateChirality(expression.intVal);
        return evalAtomPrim(expression, graph, atom) ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
    }
    if (expression.op == ExprOp::NOT) {
        const auto result = evalAtomExprWithStereo(expression.children[0], graph, atom, evaluateChirality);
        return result == StereoTruth::UNKNOWN ? result
            : (result == StereoTruth::TRUE_VALUE ? StereoTruth::FALSE_VALUE : StereoTruth::TRUE_VALUE);
    }
    bool unknown = false;
    for (const auto& child : expression.children) {
        const auto result = evalAtomExprWithStereo(child, graph, atom, evaluateChirality);
        if (expression.op == ExprOp::AND && result == StereoTruth::FALSE_VALUE) return result;
        if (expression.op == ExprOp::OR && result == StereoTruth::TRUE_VALUE) return result;
        if (result == StereoTruth::UNKNOWN) unknown = true;
    }
    if (unknown) return StereoTruth::UNKNOWN;
    return expression.op == ExprOp::AND ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
}

inline bool hasChiralityPrimitive(const AtomExpr& expression) {
    if (expression.op == ExprOp::PRIM) return expression.primType == AtomPrimType::CHIRALITY;
    for (const auto& child : expression.children) if (hasChiralityPrimitive(child)) return true;
    return false;
}

inline bool hasDirectionalPrimitive(const BondExpr& expression) {
    if (expression.op == ExprOp::PRIM)
        return expression.primType == BondPrimType::WEDGE || expression.primType == BondPrimType::DASH;
    for (const auto& child : expression.children) if (hasDirectionalPrimitive(child)) return true;
    return false;
}

/** Evaluate a directional bond expression for one candidate slash direction. */
inline bool evalBondExprWithDirection(const BondExpr& expression, const MolGraph& graph,
                                      int first, int second, int direction) {
    if (expression.op == ExprOp::PRIM) {
        if (expression.primType == BondPrimType::WEDGE || expression.primType == BondPrimType::DASH)
            return graph.bondOrder(first, second) == 1 && !graph.bondAromatic(first, second)
                && direction == (expression.primType == BondPrimType::WEDGE ? 1 : -1);
        return evalBondPrim(expression, graph, first, second);
    }
    if (expression.op == ExprOp::NOT)
        return !evalBondExprWithDirection(expression.children[0], graph, first, second, direction);
    for (const auto& child : expression.children) {
        const bool matched = evalBondExprWithDirection(child, graph, first, second, direction);
        if (expression.op == ExprOp::AND && !matched) return false;
        if (expression.op == ExprOp::OR && matched) return true;
    }
    return expression.op == ExprOp::AND;
}

// ============================================================================
// SMARTS parser
// ============================================================================

class SmartsParser {
    const std::string& smarts_;
    int pos_ = 0;
    std::vector<std::vector<int>> neighborOrder_;
    std::vector<bool> hasIncomingBond_;

    void recordAtomOrder(int atom, int parent) {
        neighborOrder_.emplace_back();
        hasIncomingBond_.push_back(parent >= 0);
        if (parent >= 0) {
            neighborOrder_[parent].push_back(atom);
            neighborOrder_[atom].push_back(parent);
        }
    }

    char peek() const {
        if (pos_ >= static_cast<int>(smarts_.size())) return '\0';
        return smarts_[pos_];
    }

    char advance() {
        if (pos_ >= static_cast<int>(smarts_.size())) return '\0';
        return smarts_[pos_++];
    }

    bool match(char c) {
        if (peek() == c) { pos_++; return true; }
        return false;
    }

    int parseDigits() {
        return parser_detail::readDecimal(smarts_, pos_);
    }

    // Try to parse a range query {min-max}, {min-}, or {-max}.
    // Returns true if a range was parsed; fills min/max (-1 = unbounded).
    bool tryParseRange(int& rmin, int& rmax) {
        if (peek() != '{') return false;
        int saved = pos_;
        advance(); // consume '{'
        rmin = -1;
        rmax = -1;
        int d = parseDigits();
        if (d >= 0) rmin = d;
        if (peek() == '-') {
            advance(); // consume '-'
            d = parseDigits();
            if (d >= 0) rmax = d;
        } else {
            // {n} without dash means exact match: min == max
            if (rmin >= 0) rmax = rmin;
        }
        if (peek() == '}') {
            advance(); // consume '}'
            return true;
        }
        // Not a valid range, rewind
        pos_ = saved;
        return false;
    }

    // After parsing a numeric primitive, check for optional range {n-m}.
    // If range found, set expr fields accordingly. Otherwise keep exact value.
    void applyOptionalRange(AtomExpr& expr) {
        int rmin, rmax;
        if (tryParseRange(rmin, rmax)) {
            expr.rangeQuery = true;
            expr.rangeMin = rmin;
            expr.rangeMax = rmax;
            expr.hasValue = true;
        }
    }

    // Check if current position starts a two-letter element (e.g., Al, Br, Cl)
    bool isTwoLetterElement(char upper) const {
        if (pos_ + 1 >= static_cast<int>(smarts_.size())) return false;
        char next = smarts_[pos_ + 1];
        if (!std::islower(next)) return false;
        std::string twoChar{upper, next};
        return elementToAtomicNum(twoChar) >= 0;
    }

    // Parse an element symbol (uppercase or lowercase, 1 or 2 chars)
    AtomExpr parseElementSymbol() {
        char c = peek();
        if (std::isupper(c)) {
            advance();
            std::string sym(1, c);
            if (pos_ < static_cast<int>(smarts_.size()) && std::islower(peek())) {
                std::string twoChar = sym + peek();
                if (elementToAtomicNum(twoChar) >= 0) {
                    sym = twoChar;
                    advance();
                }
            }
            int z = elementToAtomicNum(sym);
            if (z < 0) throw std::invalid_argument("Unknown element in SMARTS: " + sym);
            return AtomExpr::makePrim(AtomPrimType::ELEMENT, z, true);
        }
        if (std::islower(c)) {
            advance();
            std::string sym(1, c);
            if (pos_ < static_cast<int>(smarts_.size()) && std::islower(peek())) {
                std::string twoChar = sym + peek();
                if (elementToAtomicNum(twoChar) >= 0) {
                    sym = twoChar;
                    advance();
                }
            }
            int z = elementToAtomicNum(sym);
            if (z < 0) throw std::invalid_argument("Unknown aromatic element in SMARTS: " + sym);
            // Aromatic element: AND(element, aromatic)
            std::vector<AtomExpr> parts;
            parts.push_back(AtomExpr::makePrim(AtomPrimType::ELEMENT, z, true));
            parts.push_back(AtomExpr::makePrim(AtomPrimType::AROMATIC));
            return AtomExpr::makeAnd(std::move(parts));
        }
        throw std::invalid_argument(
            std::string("Expected element symbol in SMARTS at '") + c + "'");
    }

    // Parse a single atom primitive inside brackets
    AtomExpr parsePrimitive() {
        char c = peek();

        // # followed by digits = atomic number
        if (c == '#') {
            advance();
            int num = parseDigits();
            if (num < 0) throw std::invalid_argument("Expected number after # in SMARTS");
            return AtomExpr::makePrim(AtomPrimType::ATOMIC_NUM, num, true);
        }

        // * = any atom
        if (c == '*') {
            advance();
            return AtomExpr::makePrim(AtomPrimType::TRUE_);
        }

        // A = aliphatic (but Al, Ag, As, Au, Ac, Am, Ar, At = elements)
        if (c == 'A') {
            if (isTwoLetterElement('A')) return parseElementSymbol();
            advance();
            return AtomExpr::makePrim(AtomPrimType::ALIPHATIC);
        }

        // a = aromatic atom primitive
        if (c == 'a') {
            advance();
            return AtomExpr::makePrim(AtomPrimType::AROMATIC);
        }

        // X or Q = total connectivity (degree + implicit H; but Xe = xenon)
        // Q is Ehrlich-Rarey/Abolmaali synonym for X (e.g. CQ1H2 = CX1H2)
        if (c == 'X' || c == 'Q') {
            if (c == 'X' && isTwoLetterElement('X')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;
            auto expr = AtomExpr::makePrim(AtomPrimType::TOTAL_CONNECTIVITY, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // D = degree (but Dy, Db, Ds = elements)
        if (c == 'D') {
            if (isTwoLetterElement('D')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // D without number = D1
            auto expr = AtomExpr::makePrim(AtomPrimType::DEGREE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // d = heavy degree (non-hydrogen degree) — inside brackets only
        // Note: 'd' could be deuterium in some contexts but we treat it as
        // a primitive inside bracket atoms (per OpenSMARTS convention).
        if (c == 'd') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // d without number = d1
            auto expr = AtomExpr::makePrim(AtomPrimType::HEAVY_DEGREE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // v = valence
        if (c == 'v') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;
            auto expr = AtomExpr::makePrim(AtomPrimType::VALENCE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // R = ring membership (but Rb, Re, Rf, Rg, Rh, Rn, Ru = elements)
        if (c == 'R') {
            if (isTwoLetterElement('R')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::IN_RING, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::IN_RING, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // r = ring size
        if (c == 'r') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::RING_SIZE, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::RING_SIZE, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // x = ring connectivity
        if (c == 'x') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::RING_CONNECTIVITY, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::RING_CONNECTIVITY, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // z = heteroatom neighbour count
        if (c == 'z') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                // z without number = at least 1 heteroatom neighbour
                auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // Z = aliphatic heteroatom neighbour count (but Zn, Zr = elements)
        // Only treat as primitive inside brackets when followed by digit, '{', or
        // at end/delimiter position (not followed by lowercase letter forming element).
        if (c == 'Z') {
            if (isTwoLetterElement('Z')) return parseElementSymbol();
            advance();
            int num = parseDigits();
            if (num < 0) {
                auto expr = AtomExpr::makePrim(AtomPrimType::ALIPHATIC_HETERO_NEIGHBORS, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::ALIPHATIC_HETERO_NEIGHBORS, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // ^ = hybridization (^0=S, ^1=SP, ^2=SP2, ^3=SP3, ^4=SP3D, ^5=SP3D2)
        if (c == '^') {
            advance();
            int num = parseDigits();
            if (num < 0 || num > 5) {
                throw std::invalid_argument("Expected digit 0-5 after ^ in SMARTS hybridization");
            }
            return AtomExpr::makePrim(AtomPrimType::HYBRIDIZATION, num, true);
        }

        // + = positive charge
        if (c == '+') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                num = 1;
                while (peek() == '+') { advance(); num++; }
            }
            return AtomExpr::makePrim(AtomPrimType::CHARGE, num, true);
        }

        // - = negative charge (inside brackets only, not a bond)
        if (c == '-') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                num = 1;
                while (peek() == '-') { advance(); num++; }
            }
            return AtomExpr::makePrim(AtomPrimType::CHARGE, -num, true);
        }

        // H = hydrogen count (but Hf, Hg, Ho, Hs = elements)
        if (c == 'H') {
            if (pos_ + 1 < static_cast<int>(smarts_.size())) {
                char next = smarts_[pos_ + 1];
                if (next == 'f' || next == 'g' || next == 'o' || next == 's') {
                    return parseElementSymbol();
                }
            }
            advance();
            int num = parseDigits();
            if (num < 0) num = 1;  // H without number = H1
            auto expr = AtomExpr::makePrim(AtomPrimType::HCOUNT, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // @ = chirality (CW), @@ = chirality (CCW), @? = unspecified/either
        if (c == '@') {
            advance();
            if (pos_ < static_cast<int>(smarts_.size()) && peek() == '@') {
                advance(); // @@
                return AtomExpr::makePrim(AtomPrimType::CHIRALITY, 2, true); // CCW
            }
            if (pos_ < static_cast<int>(smarts_.size()) && peek() == '?') {
                advance(); // @?
                return AtomExpr::makePrim(AtomPrimType::CHIRALITY, 0, true); // unspecified
            }
            return AtomExpr::makePrim(AtomPrimType::CHIRALITY, 1, true); // CW
        }

        // :n = atom class / map number
        if (c == ':') {
            advance();
            int num = parseDigits();
            if (num < 0) num = 0;
            return AtomExpr::makePrim(AtomPrimType::ATOM_CLASS, num, true);
        }

        // $ = recursive SMARTS
        if (c == '$') {
            advance();
            if (peek() != '(') throw std::invalid_argument("Expected '(' after '$' in recursive SMARTS");
            advance(); // consume '('
            int depth = 1;
            int start = pos_;
            while (pos_ < static_cast<int>(smarts_.size()) && depth > 0) {
                if (smarts_[pos_] == '(') depth++;
                else if (smarts_[pos_] == ')') depth--;
                if (depth > 0) pos_++;
            }
            if (depth != 0) throw std::invalid_argument("Unmatched '(' in recursive SMARTS");
            std::string subSmarts = smarts_.substr(start, pos_ - start);
            advance(); // consume ')'

            AtomExpr expr;
            expr.op = ExprOp::PRIM;
            expr.primType = AtomPrimType::RECURSIVE;
            expr.recursiveSmarts = subSmarts;
            return expr;
        }

        // Element symbol (uppercase or lowercase)
        if (std::isalpha(c)) {
            return parseElementSymbol();
        }

        // Isotope: digits at the start of a bracket atom
        if (std::isdigit(c)) {
            int num = parseDigits();
            return AtomExpr::makePrim(AtomPrimType::MASS, num, true);
        }

        throw std::invalid_argument(
            std::string("Unexpected character in SMARTS atom: '") + c + "'");
    }

    // Parse NOT-level: !prim
    AtomExpr parseNotExpr() {
        bool negated = false;
        while (peek() == '!') {
            advance();
            negated = !negated;
        }
        AtomExpr expression = parsePrimitive();
        return negated ? AtomExpr::makeNot(std::move(expression)) : std::move(expression);
    }

    // Parse AND-high: prim & prim & ... (implicit AND between consecutive primitives)
    AtomExpr parseAndHigh() {
        std::vector<AtomExpr> parts;
        parts.push_back(parseNotExpr());

        while (pos_ < static_cast<int>(smarts_.size())) {
            char c = peek();
            if (c == '&') {
                advance();
                parts.push_back(parseNotExpr());
            } else if (c == ',' || c == ';' || c == ']' || c == '\0') {
                break;
            } else {
                // Implicit AND: adjacent primitives
                parts.push_back(parseNotExpr());
            }
        }
        return AtomExpr::makeAnd(std::move(parts));
    }

    // Parse OR: ... , ...
    AtomExpr parseOr() {
        std::vector<AtomExpr> parts;
        parts.push_back(parseAndHigh());

        while (peek() == ',') {
            advance();
            parts.push_back(parseAndHigh());
        }
        return AtomExpr::makeOr(std::move(parts));
    }

    // Parse AND-low (semicolon): ... ; ...
    AtomExpr parseAndLow() {
        std::vector<AtomExpr> parts;
        parts.push_back(parseOr());

        while (peek() == ';') {
            advance();
            parts.push_back(parseOr());
        }
        return AtomExpr::makeAnd(std::move(parts));
    }

    // Parse bracket atom expression: everything between [ and ]
    AtomExpr parseBracketAtomExpr() {
        return parseAndLow();
    }

    // Bond precedence: negation, high AND, OR, low AND.
    BondExpr parseBondExpr() {
        std::vector<BondExpr> parts;
        parts.push_back(parseBondOr());
        while (peek() == ';') {
            advance();
            parts.push_back(parseBondOr());
        }
        return BondExpr::makeAnd(std::move(parts));
    }

    BondExpr parseBondOr() {
        std::vector<BondExpr> parts;
        parts.push_back(parseBondAndHigh());
        while (peek() == ',') {
            advance();
            parts.push_back(parseBondAndHigh());
        }
        return BondExpr::makeOr(std::move(parts));
    }

    BondExpr parseBondAndHigh() {
        std::vector<BondExpr> parts;
        parts.push_back(parseSingleBond());
        while (peek() == '&' || isBondChar(peek())) {
            if (peek() == '&') advance();
            parts.push_back(parseSingleBond());
        }
        return BondExpr::makeAnd(std::move(parts));
    }

    BondExpr parseSingleBond() {
        bool negated = false;
        while (peek() == '!') {
            advance();
            negated = !negated;
        }
        BondPrimType type;
        switch (peek()) {
            case '-': type = BondPrimType::SINGLE; break;
            case '=': type = BondPrimType::DOUBLE; break;
            case '#': type = BondPrimType::TRIPLE; break;
            case ':': type = BondPrimType::AROMATIC; break;
            case '~': type = BondPrimType::ANY; break;
            case '@': type = BondPrimType::RING; break;
            case '/': type = BondPrimType::WEDGE; break;
            case '\\': type = BondPrimType::DASH; break;
            default: throw std::invalid_argument("Expected bond primitive in SMARTS");
        }
        advance();
        BondExpr expression = BondExpr::makePrim(type);
        return negated ? BondExpr::makeNot(std::move(expression)) : std::move(expression);
    }

    bool isBondChar(char c) const {
        return c == '-' || c == '=' || c == '#' || c == ':' || c == '~'
            || c == '@' || c == '!' || c == '/' || c == '\\';
    }

    bool isAtomStart(char c) const {
        return c == '[' || c == '*' || std::isalpha(c);
    }

public:
    SmartsParser(const std::string& smarts) : smarts_(smarts) {}

    /** Normalise parsed tetrahedral primitives to the query adjacency order. */
    void normaliseChirality(std::vector<SmartsAtom>& atoms,
                            const std::vector<std::vector<int>>& neighbors) const {
        std::function<void(AtomExpr&)> flip = [&](AtomExpr& expression) {
            if (expression.op == ExprOp::PRIM && expression.primType == AtomPrimType::CHIRALITY
                && expression.intVal != 0) expression.intVal = 3 - expression.intVal;
            for (auto& child : expression.children) flip(child);
        };
        for (int atom = 0; atom < static_cast<int>(atoms.size()); ++atom) {
            if (!hasChiralityPrimitive(atoms[atom].expr)) continue;
            std::array<int, 4> storedOrder{};
            if (!stereo_detail::tetrahedralLigands(neighbors[atom], atom, storedOrder)) continue;
            auto parsedOrder = neighborOrder_[atom];
            if (parsedOrder.size() == 3)
                parsedOrder.insert(parsedOrder.begin() + (hasIncomingBond_[atom] ? 1 : 0), -1);
            int inversions = 0;
            for (size_t first = 0; first < parsedOrder.size(); ++first) {
                const auto firstPosition = std::find(storedOrder.begin(), storedOrder.end(), parsedOrder[first]);
                for (size_t second = first + 1; second < parsedOrder.size(); ++second) {
                    const auto secondPosition = std::find(storedOrder.begin(), storedOrder.end(), parsedOrder[second]);
                    inversions += firstPosition > secondPosition;
                }
            }
            if (inversions % 2 != 0) flip(atoms[atom].expr);
        }
    }

    // Parse the full SMARTS string into atoms and bonds
    void parse(std::vector<SmartsAtom>& atoms, std::vector<SmartsBond>& bonds) {
        if (smarts_.empty()) return;

        std::vector<int> branchStack;   // stack of parent atom indices
        int currentAtom = -1;
        struct RingClosure { int atom; BondExpr expression; int neighborSlot; };
        std::unordered_map<int, RingClosure> ringClosures;
        BondExpr pendingBond;
        bool hasPendingBond = false;

        while (pos_ < static_cast<int>(smarts_.size())) {
            char c = peek();

            if (c == '(') {
                advance();
                branchStack.push_back(currentAtom);
                hasPendingBond = false;
                continue;
            }

            if (c == ')') {
                advance();
                if (branchStack.empty()) {
                    throw std::invalid_argument("Unmatched ')' in SMARTS");
                }
                currentAtom = branchStack.back();
                branchStack.pop_back();
                hasPendingBond = false;
                continue;
            }

            // Bond symbol?
            if (isBondChar(c)) {
                pendingBond = parseBondExpr();
                hasPendingBond = true;
                continue;
            }

            // Ring closure digit or %nn
            if (std::isdigit(c) || c == '%') {
                int ringNum;
                if (c == '%') {
                    advance();
                    ringNum = parseDigits();
                    if (ringNum < 0) throw std::invalid_argument("Expected digits after % in ring closure");
                } else {
                    ringNum = advance() - '0';
                }

                if (currentAtom < 0) {
                    throw std::invalid_argument("Ring closure before any atom");
                }

                auto it = ringClosures.find(ringNum);
                if (it != ringClosures.end()) {
                    // Close the ring
                    SmartsBond bond;
                    const int openAtom = it->second.atom;
                    bond.from = hasPendingBond ? currentAtom : openAtom;
                    bond.to = hasPendingBond ? openAtom : currentAtom;
                    neighborOrder_[openAtom][it->second.neighborSlot] = currentAtom;
                    neighborOrder_[currentAtom].push_back(openAtom);
                    // Use the bond from open if it had one, or the current pending bond
                    if (hasPendingBond) {
                        bond.expr = std::move(pendingBond);
                        hasPendingBond = false;
                    } else if (it->second.expression.primType != BondPrimType::DEFAULT_
                               || it->second.expression.op != ExprOp::PRIM) {
                        bond.expr = std::move(it->second.expression);
                    } else {
                        bond.expr = BondExpr::makePrim(BondPrimType::DEFAULT_);
                    }
                    bonds.push_back(std::move(bond));
                    ringClosures.erase(it);
                } else {
                    // Open the ring
                    BondExpr ringBondExpr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    const int slot = static_cast<int>(neighborOrder_[currentAtom].size());
                    neighborOrder_[currentAtom].push_back(-1);
                    ringClosures[ringNum] = {currentAtom, std::move(ringBondExpr), slot};
                }
                continue;
            }

            // Atom: [ bracket atom ] or organic shorthand
            if (c == '[') {
                advance(); // consume [

                SmartsAtom atom;
                atom.expr = parseBracketAtomExpr();

                if (peek() != ']') {
                    throw std::invalid_argument(
                        "Expected ']' in bracket atom at pos " + std::to_string(pos_));
                }
                advance(); // consume ]

                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            // Organic shorthand atoms: B C N O P S F Cl Br I * and aromatic b c n o p s
            if (c == '*') {
                advance();
                SmartsAtom atom;
                atom.expr = AtomExpr::makePrim(AtomPrimType::TRUE_);
                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            if (std::isalpha(c)) {
                // Outside brackets, a and A are always atom wildcards.
                if (c == 'a' || c == 'A') {
                    advance();
                    SmartsAtom atom;
                    atom.isAromatic = (c == 'a');
                    atom.expr = AtomExpr::makePrim(c == 'a' ? AtomPrimType::AROMATIC
                                                           : AtomPrimType::ALIPHATIC, 0, true);
                    int newIdx = static_cast<int>(atoms.size());
                    atoms.push_back(std::move(atom));
                    recordAtomOrder(newIdx, currentAtom);
                    if (currentAtom >= 0) {
                        SmartsBond bond;
                        bond.from = currentAtom; bond.to = newIdx;
                        bond.expr = hasPendingBond ? std::move(pendingBond)
                            : BondExpr::makePrim(BondPrimType::DEFAULT_);
                        hasPendingBond = false;
                        bonds.push_back(std::move(bond));
                    }
                    currentAtom = newIdx;
                    continue;
                }
                bool aromatic = std::islower(c);
                std::string sym(1, c);
                advance();

                if ((sym == "C" && peek() == 'l') || (sym == "B" && peek() == 'r')) {
                    sym += advance();
                }

                int z = elementToAtomicNum(sym);
                if (z < 0) {
                    throw std::invalid_argument("Unknown element in SMARTS: " + sym);
                }

                SmartsAtom atom;
                atom.isAromatic = aromatic;
                if (aromatic) {
                    std::vector<AtomExpr> parts;
                    parts.push_back(AtomExpr::makePrim(AtomPrimType::ELEMENT, z, true));
                    parts.push_back(AtomExpr::makePrim(AtomPrimType::AROMATIC));
                    atom.expr = AtomExpr::makeAnd(std::move(parts));
                } else {
                    // Aliphatic organic atom: element AND aliphatic
                    std::vector<AtomExpr> parts;
                    parts.push_back(AtomExpr::makePrim(AtomPrimType::ELEMENT, z, true));
                    parts.push_back(AtomExpr::makePrim(AtomPrimType::ALIPHATIC));
                    atom.expr = AtomExpr::makeAnd(std::move(parts));
                }

                int newIdx = static_cast<int>(atoms.size());
                atoms.push_back(std::move(atom));
                recordAtomOrder(newIdx, currentAtom);

                if (currentAtom >= 0) {
                    SmartsBond bond;
                    bond.from = currentAtom;
                    bond.to = newIdx;
                    bond.expr = hasPendingBond
                        ? std::move(pendingBond)
                        : BondExpr::makePrim(BondPrimType::DEFAULT_);
                    hasPendingBond = false;
                    bonds.push_back(std::move(bond));
                }
                currentAtom = newIdx;
                continue;
            }

            // Skip whitespace
            if (std::isspace(static_cast<unsigned char>(c))) {
                advance();
                continue;
            }

            // Dot (disconnected)
            if (c == '.') {
                advance();
                currentAtom = -1;
                hasPendingBond = false;
                continue;
            }

            throw std::invalid_argument(
                std::string("Unexpected character in SMARTS: '") + c
                + "' at position " + std::to_string(pos_));
        }

        if (!branchStack.empty()) {
            throw std::invalid_argument("Unmatched '(' in SMARTS");
        }
        if (!ringClosures.empty()) {
            throw std::invalid_argument("Unclosed ring in SMARTS");
        }
    }
};

} // namespace smarts_detail

// ============================================================================
// SmartsQuery -- the parsed SMARTS pattern
// ============================================================================


struct SmartsQuery {
    std::vector<smarts_detail::SmartsAtom> atoms;
    std::vector<smarts_detail::SmartsBond> bonds;

    // Adjacency for the query graph
    std::vector<std::vector<int>> neighbors;   // neighbour atom indices
    std::vector<std::vector<int>> bondIndices; // parallel to neighbours: bond index

    int atomCount() const { return static_cast<int>(atoms.size()); }

    // Build adjacency from atoms/bonds
    void buildAdjacency();

    // Primitive compatibility helpers (kept public for backwards compatibility)
    bool atomCompat(int qi, const MolGraph& target, int ti) const;
    bool bondCompat(int bi, const MolGraph& target, int ti, int tj) const;

    /**
     * Historical/debug ordering: plain BFS from atom 0.
     *
     * The optimised matcher below uses a target-aware, anchor-preserving order;
     * this method remains intentionally simple for callers that relied on the
     * previous behaviour.
     */
    std::vector<int> matchingOrder() const;


    bool hasStereochemistry() const;
    bool stereoCompatible(const MolGraph& target, const std::vector<int>& mapping) const;

};

inline void SmartsQuery::buildAdjacency() {
    int n = atomCount();
    neighbors.assign(n, {});
    bondIndices.assign(n, {});
    for (int bi = 0; bi < static_cast<int>(bonds.size()); ++bi) {
        int a = bonds[bi].from, b = bonds[bi].to;
        neighbors[a].push_back(b);
        neighbors[b].push_back(a);
        bondIndices[a].push_back(bi);
        bondIndices[b].push_back(bi);
    }
}

inline bool SmartsQuery::hasStereochemistry() const {
    for (const auto& atom : atoms)
        if (smarts_detail::hasChiralityPrimitive(atom.expr)) return true;
    for (const auto& bond : bonds)
        if (smarts_detail::hasDirectionalPrimitive(bond.expr)) return true;
    return false;
}

inline bool SmartsQuery::stereoCompatible(const MolGraph& target, const std::vector<int>& mapping) const {
    using smarts_detail::StereoTruth;
    for (int atom = 0; atom < atomCount(); ++atom) {
        if (!smarts_detail::hasChiralityPrimitive(atoms[atom].expr)) continue;
        const int targetAtom = mapping[atom];
        std::array<int, 4> queryLigands{}, targetLigands{};
        const bool validLigands = stereo_detail::tetrahedralLigands(neighbors[atom], atom, queryLigands)
            && stereo_detail::tetrahedralLigands(target.neighbors[targetAtom], targetAtom, targetLigands);
        const int targetTag = target.tetraChirality.empty() ? 0 : target.tetraChirality[targetAtom];
        const auto matched = smarts_detail::evalAtomExprWithStereo(atoms[atom].expr, target, targetAtom, [&](int queryTag) {
            if (queryTag == 0) return StereoTruth::TRUE_VALUE;
            if (targetTag == 0) return StereoTruth::FALSE_VALUE;
            const bool compatible = validLigands
                ? stereo_detail::tetrahedralParityMatches(queryTag, queryLigands, targetTag, targetLigands, mapping)
                : queryTag == targetTag;
            return compatible ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
        });
        if (matched != StereoTruth::TRUE_VALUE) return false;
    }

    std::vector<int> variableForBond(bonds.size(), -1), allowedDirections;
    for (int bondIndex = 0; bondIndex < static_cast<int>(bonds.size()); ++bondIndex) {
        const auto& bond = bonds[bondIndex];
        if (!smarts_detail::hasDirectionalPrimitive(bond.expr)) continue;
        const int first = mapping[bond.from], second = mapping[bond.to];
        if (target.bondOrder(first, second) != 1 || target.bondAromatic(first, second)) continue;
        int allowed = 0;
        if (smarts_detail::evalBondExprWithDirection(bond.expr, target, first, second, 1)) allowed |= 1;
        if (smarts_detail::evalBondExprWithDirection(bond.expr, target, first, second, -1)) allowed |= 2;
        if (allowed == 0) return false;
        variableForBond[bondIndex] = static_cast<int>(allowedDirections.size());
        allowedDirections.push_back(allowed);
    }
    if (allowedDirections.empty()) return true;
    std::vector<std::vector<std::pair<int, int>>> constraints(allowedDirections.size());
    auto addConstraint = [&](int firstBond, int secondBond, int relation) {
        const int firstVariable = variableForBond[firstBond], secondVariable = variableForBond[secondBond];
        constraints[firstVariable].push_back({secondVariable, relation});
        constraints[secondVariable].push_back({firstVariable, relation});
    };
    auto outwardOrientation = [&](int bondIndex, int atom) { return bonds[bondIndex].from == atom ? 1 : -1; };
    for (const auto& focus : bonds) {
        const int first = focus.from, second = focus.to;
        const int mappedFirst = mapping[first], mappedSecond = mapping[second];
        if (target.bondOrder(mappedFirst, mappedSecond) != 2) continue;
        std::array<std::vector<int>, 2> carriers;
        for (int end = 0; end < 2; ++end) {
            const int atom = end == 0 ? first : second, other = end == 0 ? second : first;
            for (size_t position = 0; position < neighbors[atom].size(); ++position) {
                const int bondIndex = bondIndices[atom][position];
                if (neighbors[atom][position] != other && variableForBond[bondIndex] >= 0)
                    carriers[end].push_back(bondIndex);
            }
            for (size_t previous = 0; previous < carriers[end].size(); ++previous)
                for (size_t next = previous + 1; next < carriers[end].size(); ++next)
                    addConstraint(carriers[end][previous], carriers[end][next],
                        -outwardOrientation(carriers[end][previous], atom) * outwardOrientation(carriers[end][next], atom));
        }
        if (carriers[0].empty() || carriers[1].empty()) continue;
        const int configuration = target.dbStereo(mappedFirst, mappedSecond);
        if (configuration == 0) {
            bool firstRequired = false, secondRequired = false;
            for (int bondIndex : carriers[0]) if (allowedDirections[variableForBond[bondIndex]] != 3) firstRequired = true;
            for (int bondIndex : carriers[1]) if (allowedDirections[variableForBond[bondIndex]] != 3) secondRequired = true;
            if (firstRequired && secondRequired) return false;
            continue;
        }
        const int firstReference = stereo_detail::doubleBondReference(target, mappedFirst, mappedSecond);
        const int secondReference = stereo_detail::doubleBondReference(target, mappedSecond, mappedFirst);
        for (int firstBond : carriers[0]) {
            const auto& firstCarrier = bonds[firstBond];
            const int firstLigand = firstCarrier.from == first ? firstCarrier.to : firstCarrier.from;
            for (int secondBond : carriers[1]) {
                const auto& secondCarrier = bonds[secondBond];
                const int secondLigand = secondCarrier.from == second ? secondCarrier.to : secondCarrier.from;
                const int relation = (configuration == 1 ? 1 : -1)
                    * (mapping[firstLigand] == firstReference ? 1 : -1)
                    * (mapping[secondLigand] == secondReference ? 1 : -1)
                    * outwardOrientation(firstBond, first) * outwardOrientation(secondBond, second);
                addConstraint(firstBond, secondBond, relation);
            }
        }
    }
    std::vector<int> relativeDirections(allowedDirections.size(), 0);
    for (int start = 0; start < static_cast<int>(allowedDirections.size()); ++start) {
        if (relativeDirections[start] != 0) continue;
        relativeDirections[start] = 1;
        int allowedRoot = 3;
        std::deque<int> pending{start};
        while (!pending.empty()) {
            const int current = pending.front(); pending.pop_front();
            const int allowed = allowedDirections[current];
            allowedRoot &= relativeDirections[current] == 1 ? allowed : ((allowed & 1) << 1) | ((allowed & 2) >> 1);
            if (allowedRoot == 0) return false;
            for (const auto& constraint : constraints[current]) {
                const int required = relativeDirections[current] * constraint.second;
                if (relativeDirections[constraint.first] == 0) {
                    relativeDirections[constraint.first] = required;
                    pending.push_back(constraint.first);
                } else if (relativeDirections[constraint.first] != required) return false;
            }
        }
    }
    return true;
}

inline bool SmartsQuery::atomCompat(int qi, const MolGraph& target, int ti) const {
    return smarts_detail::evalAtomExpr(atoms[qi].expr, target, ti);
}

inline bool SmartsQuery::bondCompat(int bi, const MolGraph& target, int ti, int tj) const {
    const auto& expression = bonds[bi].expr;
    if (smarts_detail::hasDirectionalPrimitive(expression))
        return smarts_detail::evalBondExprWithDirection(expression, target, ti, tj, 1)
            || smarts_detail::evalBondExprWithDirection(expression, target, ti, tj, -1);
    return smarts_detail::evalBondExpr(expression, target, ti, tj);
}

inline std::vector<int> SmartsQuery::matchingOrder() const {
    int n = atomCount();
    if (n == 0) return {};
    std::vector<int> order;
    order.reserve(n);
    std::vector<bool> visited(n, false);
    std::deque<int> queue;
    queue.push_back(0);
    visited[0] = true;
    while (!queue.empty()) {
        int u = queue.front(); queue.pop_front();
        order.push_back(u);
        for (int v : neighbors[u]) {
            if (!visited[v]) {
                visited[v] = true;
                queue.push_back(v);
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (!visited[i]) {
            visited[i] = true;
            queue.push_back(i);
            while (!queue.empty()) {
                int u = queue.front(); queue.pop_front();
                order.push_back(u);
                for (int v : neighbors[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        queue.push_back(v);
                    }
                }
            }
        }
    }
    return order;
}

// Deferred recursive SMARTS evaluation (needs complete SmartsQuery type)
inline bool smarts_detail::evalRecursiveSmarts(const smarts_detail::AtomExpr& expr, const MolGraph& g, int idx) {
    if (expr.recursiveSmarts.empty()) return false;

    if (smarts_detail::recursiveMatcher) {
        if (expr.cachedRecursive) return smarts_detail::recursiveMatcher(*expr.cachedRecursive, g, idx);
        return smarts_detail::recursiveMatcher(parseSMARTS(expr.recursiveSmarts), g, idx);
    }
    throw std::logic_error("Recursive SMARTS evaluation requires a Substructure search context");
}

// ============================================================================
// Named predicate registry (mirrors Java Standardiser.expandNamedPredicates)
// ============================================================================

class SmartsPredicateRegistry {
    std::unordered_map<std::string, std::string> predicates_;

public:
    SmartsPredicateRegistry() = default;

    /// Register a named predicate: name -> SMARTS definition.
    void registerPredicate(const std::string& name, const std::string& smarts) {
        predicates_[name] = smarts;
    }

    /// Expand all `$name` tokens (inside brackets) with `$(definition)`.
    /// Scans for `$` followed by an identifier (alphanumeric + underscore);
    /// if the identifier is a registered predicate, replaces `$name` with
    /// `$(definition)`.  Tokens already of the form `$(...)` are left alone.
    std::string expandPredicates(const std::string& input) const {
        std::string result;
        result.reserve(input.size());
        size_t i = 0;
        while (i < input.size()) {
            if (input[i] == '$' && i + 1 < input.size() && input[i + 1] != '(') {
                // Collect the identifier after '$'
                size_t start = i + 1;
                size_t end = start;
                while (end < input.size() &&
                       (std::isalnum(static_cast<unsigned char>(input[end])) || input[end] == '_')) {
                    ++end;
                }
                std::string name = input.substr(start, end - start);
                auto it = predicates_.find(name);
                if (it != predicates_.end()) {
                    result += "$(";
                    result += it->second;
                    result += ')';
                    i = end;
                } else {
                    result += input[i];
                    ++i;
                }
            } else {
                result += input[i];
                ++i;
            }
        }
        return result;
    }

    /// Create a registry pre-loaded with the 10 built-in predicates
    /// (same set as Java Standardiser).
    static SmartsPredicateRegistry builtinRegistry() {
        SmartsPredicateRegistry reg;
        reg.registerPredicate("isAmideN",
            "[$([NX3;H2,H1;!$(NC=O)]),$([NX3;$(NC=O)])]");
        reg.registerPredicate("isCarboxylC",
            "[CX3](=O)[OX2H1,OX1-]");
        reg.registerPredicate("isSulfonamideN",
            "[NX3;$(NS(=O)=O)]");
        reg.registerPredicate("isEster",
            "[#6][CX3](=O)[OX2][#6]");
        reg.registerPredicate("isKetone",
            "[#6][CX3](=O)[#6]");
        reg.registerPredicate("isHalogen",
            "[F,Cl,Br,I]");
        reg.registerPredicate("isPositive",
            "[*+]");
        reg.registerPredicate("isNegative",
            "[*-]");
        reg.registerPredicate("isHBDonor",
            "[N!H0,O!H0,S!H0]");
        reg.registerPredicate("isHBAcceptor",
            "[N,O,S;H0;!+0]");
        return reg;
    }
};

// ============================================================================
// parseSMARTS -- top-level API
// ============================================================================

inline SmartsQuery parseSMARTS(const std::string& smarts, int maxRecursionDepth) {
    if (smarts.empty()) return SmartsQuery{};

    // Guard against unbounded recursive SMARTS nesting by tracking actual
    // nested $(...) frames rather than counting raw occurrences.
    int parenDepth = 0;
    std::vector<int> recursiveFrames;
    recursiveFrames.reserve(8);
    for (size_t i = 0; i < smarts.size(); ++i) {
        char c = smarts[i];
        if (c == '(') {
            ++parenDepth;
            if (i > 0 && smarts[i - 1] == '$') {
                recursiveFrames.push_back(parenDepth);
                if (static_cast<int>(recursiveFrames.size()) > maxRecursionDepth) {
                    throw std::invalid_argument(
                        "Recursive SMARTS nesting exceeds maximum depth of "
                        + std::to_string(maxRecursionDepth));
                }
            }
        } else if (c == ')') {
            if (!recursiveFrames.empty() && recursiveFrames.back() == parenDepth)
                recursiveFrames.pop_back();
            if (parenDepth > 0) --parenDepth;
        }
    }

    SmartsQuery query;
    smarts_detail::SmartsParser parser(smarts);
    parser.parse(query.atoms, query.bonds);
    query.buildAdjacency();
    parser.normaliseChirality(query.atoms, query.neighbors);
    std::function<void(smarts_detail::AtomExpr&)> compileRecursive = [&](smarts_detail::AtomExpr& expression) {
        if (expression.primType == smarts_detail::AtomPrimType::RECURSIVE)
            expression.cachedRecursive = std::make_shared<SmartsQuery>(parseSMARTS(expression.recursiveSmarts, maxRecursionDepth));
        for (auto& child : expression.children) compileRecursive(child);
    };
    for (auto& atom : query.atoms) compileRecursive(atom.expr);
    return query;
}

} // namespace smsd

#endif // SMSD_SMARTS_PARSER_HPP

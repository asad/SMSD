/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 * Header-only SMARTS pattern parser and matcher for smsd::MolGraph.
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
 *   - VF2-style backtracking matcher
 */
#pragma once
#ifndef SMSD_SMARTS_PARSER_HPP
#define SMSD_SMARTS_PARSER_HPP

#include "mol_graph.hpp"

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
    HETERO_NEIGHBORS,          // z<n> — count of heteroatom neighbors (not C, not H)
    ALIPHATIC_HETERO_NEIGHBORS,// Z<n> — count of aliphatic heteroatom neighbors
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
    std::string recursiveSmarts; // parsed lazily when evaluated
    mutable std::shared_ptr<SmartsQuery> cachedRecursive; // cached parse result
    mutable const MolGraph* cachedRecursiveTarget = nullptr; // last target used for recursive memo
    mutable uint64_t cachedRecursiveTargetNonce = 0;
    mutable std::vector<int8_t> cachedRecursiveState;        // -1 unknown, 0 false, 1 true per target atom

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

    int charge = std::abs(g.formalCharge[idx]);

    // Find the smallest default valence >= explicitValence + charge
    for (int dv : it->second) {
        int hc = dv - explicitValence - charge;
        if (hc >= 0) return hc;
    }
    return 0;
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
            hcount[i] = computeHCountFromValence(g, i, valence[i]);
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
            return matchNumeric(expr, g.degree[idx] + computeHCount(g, idx));

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
            // SMARTS '/' — directional bond. Matches single or double bonds.
            // When E/Z stereo is available (dbStereo != 0), checks direction;
            // otherwise matches any single/double bond (permissive fallback).
            {
                int ord = g.bondOrder(i, j);
                if (ord != 1 && ord != 2) return false;
                int stereo = g.dbStereo(i, j);
                return stereo == 0 || stereo == 1; // unspecified or same-direction
            }

        case BondPrimType::DASH:
            // SMARTS '\' — directional bond (opposite direction).
            {
                int ord = g.bondOrder(i, j);
                if (ord != 1 && ord != 2) return false;
                int stereo = g.dbStereo(i, j);
                return stereo == 0 || stereo == 2; // unspecified or opposite-direction
            }

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

// ============================================================================
// SMARTS parser
// ============================================================================

class SmartsParser {
    const std::string& smarts_;
    int pos_ = 0;

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
        int val = 0;
        bool gotDigit = false;
        while (pos_ < static_cast<int>(smarts_.size()) && std::isdigit(peek())) {
            val = val * 10 + (advance() - '0');
            gotDigit = true;
        }
        return gotDigit ? val : -1;
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

        // z = heteroatom neighbor count
        if (c == 'z') {
            advance();
            int num = parseDigits();
            if (num < 0) {
                // z without number = at least 1 heteroatom neighbor
                auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, 0, false);
                applyOptionalRange(expr);
                return expr;
            }
            auto expr = AtomExpr::makePrim(AtomPrimType::HETERO_NEIGHBORS, num, true);
            applyOptionalRange(expr);
            return expr;
        }

        // Z = aliphatic heteroatom neighbor count (but Zn, Zr = elements)
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
            expr.recursiveSmarts = subSmarts; // parsed lazily in evalRecursiveSmarts
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
        if (peek() == '!') {
            advance();
            return AtomExpr::makeNot(parseNotExpr());
        }
        return parsePrimitive();
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

    // Parse a bond expression
    BondExpr parseBondExpr() {
        std::vector<BondExpr> parts;

        while (pos_ < static_cast<int>(smarts_.size())) {
            char c = peek();
            BondExpr be;
            bool gotBond = false;

            if (c == '!') {
                advance();
                be = BondExpr::makeNot(parseSingleBond());
                gotBond = true;
            } else if (c == '-') {
                advance();
                be = BondExpr::makePrim(BondPrimType::SINGLE);
                gotBond = true;
            } else if (c == '=') {
                advance();
                be = BondExpr::makePrim(BondPrimType::DOUBLE);
                gotBond = true;
            } else if (c == '#') {
                advance();
                be = BondExpr::makePrim(BondPrimType::TRIPLE);
                gotBond = true;
            } else if (c == ':') {
                advance();
                be = BondExpr::makePrim(BondPrimType::AROMATIC);
                gotBond = true;
            } else if (c == '~') {
                advance();
                be = BondExpr::makePrim(BondPrimType::ANY);
                gotBond = true;
            } else if (c == '@') {
                advance();
                be = BondExpr::makePrim(BondPrimType::RING);
                gotBond = true;
            } else if (c == '/') {
                advance();
                be = BondExpr::makePrim(BondPrimType::WEDGE);
                gotBond = true;
            } else if (c == '\\') {
                advance();
                be = BondExpr::makePrim(BondPrimType::DASH);
                gotBond = true;
            } else {
                break;
            }

            if (gotBond) {
                // Check for logical operators between bonds
                if (peek() == ',') {
                    // OR: collect all OR'd bonds
                    std::vector<BondExpr> orParts;
                    orParts.push_back(std::move(be));
                    while (peek() == ',') {
                        advance();
                        orParts.push_back(parseSingleBond());
                    }
                    parts.push_back(BondExpr::makeOr(std::move(orParts)));
                } else if (peek() == ';') {
                    // AND-low (semicolon): e.g. -;!@ = single AND NOT ring
                    std::vector<BondExpr> andParts;
                    andParts.push_back(std::move(be));
                    while (peek() == ';') {
                        advance();
                        andParts.push_back(parseSingleBond());
                    }
                    parts.push_back(BondExpr::makeAnd(std::move(andParts)));
                } else {
                    parts.push_back(std::move(be));
                }
            }
        }

        if (parts.empty()) return BondExpr::makePrim(BondPrimType::DEFAULT_);
        if (parts.size() == 1) return std::move(parts[0]);
        return BondExpr::makeAnd(std::move(parts));
    }

    BondExpr parseSingleBond() {
        char c = peek();
        if (c == '!') {
            advance();
            return BondExpr::makeNot(parseSingleBond());
        }
        if (c == '-') { advance(); return BondExpr::makePrim(BondPrimType::SINGLE); }
        if (c == '=') { advance(); return BondExpr::makePrim(BondPrimType::DOUBLE); }
        if (c == '#') { advance(); return BondExpr::makePrim(BondPrimType::TRIPLE); }
        if (c == ':') { advance(); return BondExpr::makePrim(BondPrimType::AROMATIC); }
        if (c == '~') { advance(); return BondExpr::makePrim(BondPrimType::ANY); }
        if (c == '@') { advance(); return BondExpr::makePrim(BondPrimType::RING); }
        if (c == '/') { advance(); return BondExpr::makePrim(BondPrimType::WEDGE); }
        if (c == '\\') { advance(); return BondExpr::makePrim(BondPrimType::DASH); }
        return BondExpr::makePrim(BondPrimType::DEFAULT_);
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

    // Parse the full SMARTS string into atoms and bonds
    void parse(std::vector<SmartsAtom>& atoms, std::vector<SmartsBond>& bonds) {
        if (smarts_.empty()) return;

        std::vector<int> branchStack;   // stack of parent atom indices
        int currentAtom = -1;
        std::unordered_map<int, std::pair<int, BondExpr>> ringClosures;
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
                    bond.from = it->second.first;
                    bond.to = currentAtom;
                    // Use the bond from open if it had one, or the current pending bond
                    if (hasPendingBond) {
                        bond.expr = std::move(pendingBond);
                        hasPendingBond = false;
                    } else if (it->second.second.primType != BondPrimType::DEFAULT_
                               || it->second.second.op != ExprOp::PRIM) {
                        bond.expr = std::move(it->second.second);
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
                    ringClosures[ringNum] = {currentAtom, std::move(ringBondExpr)};
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
                // SMARTS wildcards: 'a' = any aromatic, 'A' = any aliphatic
                // 'a' is aromatic wildcard UNLESS followed by a lowercase letter
                // that forms a valid two-letter aromatic element (as, se, te)
                if (c == 'a') {
                    bool isElement = false;
                    if (pos_ + 1 < static_cast<int>(smarts_.size()) && std::islower(smarts_[pos_ + 1])) {
                        std::string twoChar = std::string(1, c) + smarts_[pos_ + 1];
                        if (elementToAtomicNum(twoChar) >= 0) isElement = true;
                    }
                    if (!isElement) {
                        advance();
                        SmartsAtom atom;
                        atom.expr = AtomExpr::makePrim(AtomPrimType::AROMATIC, 0, true);
                        int newIdx = static_cast<int>(atoms.size());
                        atoms.push_back(std::move(atom));
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
                }
                if (c == 'A') {
                    bool isElement = false;
                    if (pos_ + 1 < static_cast<int>(smarts_.size()) && std::islower(smarts_[pos_ + 1])) {
                        std::string twoChar = std::string(1, c) + smarts_[pos_ + 1];
                        if (elementToAtomicNum(twoChar) >= 0) isElement = true;
                    }
                    if (!isElement) {
                        advance();
                        SmartsAtom atom;
                        atom.expr = AtomExpr::makePrim(AtomPrimType::ALIPHATIC, 0, true);
                        int newIdx = static_cast<int>(atoms.size());
                        atoms.push_back(std::move(atom));
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
                }
                bool aromatic = std::islower(c);
                std::string sym(1, c);
                advance();

                if (pos_ < static_cast<int>(smarts_.size()) && std::islower(peek())) {
                    std::string twoChar = sym + peek();
                    if (elementToAtomicNum(twoChar) >= 0) {
                        sym = twoChar;
                        advance();
                    }
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
    std::vector<std::vector<int>> neighbors;   // neighbor atom indices
    std::vector<std::vector<int>> bondIndices; // parallel to neighbors: bond index

    int atomCount() const { return static_cast<int>(atoms.size()); }

    // Build adjacency from atoms/bonds
    void buildAdjacency();

    // Primitive compatibility helpers (kept public for backwards compatibility)
    bool atomCompat(int qi, const MolGraph& target, int ti) const;
    bool bondCompat(int bi, const MolGraph& target, int ti, int tj) const;

    /**
     * Historical/debug ordering: plain BFS from atom 0.
     *
     * The optimized matcher below uses a target-aware, anchor-preserving order;
     * this method remains intentionally simple for callers that relied on the
     * previous behavior.
     */
    std::vector<int> matchingOrder() const;

    bool matches(const MolGraph& target) const;
    std::vector<std::map<int,int>> findAll(const MolGraph& target, int maxMatches = 1000) const;

    // Match with atom 0 of the query anchored to a specific target atom
    // (used for recursive SMARTS).
    bool matchesAtAtom(const MolGraph& target, int targetAtom) const;
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

inline bool SmartsQuery::atomCompat(int qi, const MolGraph& target, int ti) const {
    return smarts_detail::evalAtomExpr(atoms[qi].expr, target, ti);
}

inline bool SmartsQuery::bondCompat(int bi, const MolGraph& target, int ti, int tj) const {
    return smarts_detail::evalBondExpr(bonds[bi].expr, target, ti, tj);
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

namespace smarts_detail {

class SmartsMatchEngine {
    const SmartsQuery& q_;
    const MolGraph& target_;
    int nQ_ = 0;
    int nT_ = 0;
    int tWords_ = 0;
    int anchorQ_ = -1;
    bool anyEmptyDomain_ = false;

    std::vector<std::vector<uint64_t>> domain_;
    std::vector<int> domainSize_;
    std::vector<int8_t> atomMemo_;   // -1 unknown, 0 false, 1 true

    std::vector<int> q2t_;
    std::vector<int> t2q_;
    std::vector<uint64_t> usedMask_;
    std::vector<std::vector<int>> candBuf_;

    std::vector<int> compId_;
    std::vector<std::vector<int>> qNbrsByTightness_;
    std::vector<std::vector<int>> qNbrBondByTightness_;
    std::vector<int> order_;

    mutable std::vector<uint64_t> availBuf_;
    mutable std::vector<uint64_t> tmpMaskBuf_;

    bool domainHas(int qi, int tj) const {
        return (domain_[qi][tj >> 6] & (uint64_t(1) << (tj & 63))) != 0;
    }

    int popcountWords(const std::vector<uint64_t>& bits) const {
        int cnt = 0;
        for (uint64_t w : bits) cnt += smsd::popcount64(w);
        return cnt;
    }

    void collectBits(const std::vector<uint64_t>& bits, std::vector<int>& out) const {
        out.clear();
        out.reserve(nT_);
        for (int w = 0; w < tWords_; ++w) {
            uint64_t word = bits[w];
            while (word != 0) {
                int bit = smsd::ctz64(word);
                int idx = (w << 6) | bit;
                if (idx < nT_) out.push_back(idx);
                word &= word - 1;
            }
        }
    }

    bool atomCompatCached(int qi, int tj) {
        const size_t memoIdx = static_cast<size_t>(qi) * static_cast<size_t>(nT_) + static_cast<size_t>(tj);
        int8_t& memo = atomMemo_[memoIdx];
        if (memo >= 0) return memo != 0;
        memo = q_.atomCompat(qi, target_, tj) ? int8_t(1) : int8_t(0);
        return memo != 0;
    }

    void buildDomains() {
        for (int qi = 0; qi < nQ_; ++qi) {
            int qDeg = static_cast<int>(q_.neighbors[qi].size());
            for (int tj = 0; tj < nT_; ++tj) {
                if (target_.degree[tj] < qDeg) continue;
                if (!atomCompatCached(qi, tj)) continue;
                domain_[qi][tj >> 6] |= uint64_t(1) << (tj & 63);
            }
            domainSize_[qi] = popcountWords(domain_[qi]);
            if (domainSize_[qi] == 0) anyEmptyDomain_ = true;
        }
    }

    void buildComponents() {
        compId_.assign(nQ_, -1);
        int cid = 0;
        for (int s = 0; s < nQ_; ++s) {
            if (compId_[s] >= 0) continue;
            std::deque<int> dq;
            dq.push_back(s);
            compId_[s] = cid;
            while (!dq.empty()) {
                int u = dq.front();
                dq.pop_front();
                for (int v : q_.neighbors[u]) {
                    if (compId_[v] < 0) {
                        compId_[v] = cid;
                        dq.push_back(v);
                    }
                }
            }
            ++cid;
        }
    }

    void buildNeighborOrder() {
        qNbrsByTightness_.assign(nQ_, {});
        qNbrBondByTightness_.assign(nQ_, {});
        for (int qi = 0; qi < nQ_; ++qi) {
            std::vector<std::pair<int,int>> items;
            items.reserve(q_.neighbors[qi].size());
            for (size_t ni = 0; ni < q_.neighbors[qi].size(); ++ni) {
                items.emplace_back(q_.neighbors[qi][ni], q_.bondIndices[qi][ni]);
            }
            std::sort(items.begin(), items.end(), [&](const auto& lhs, const auto& rhs) {
                int a = lhs.first;
                int b = rhs.first;
                if (domainSize_[a] != domainSize_[b]) return domainSize_[a] < domainSize_[b];
                if (q_.neighbors[a].size() != q_.neighbors[b].size())
                    return q_.neighbors[a].size() > q_.neighbors[b].size();
                return a < b;
            });
            auto& nbs = qNbrsByTightness_[qi];
            auto& bis = qNbrBondByTightness_[qi];
            nbs.reserve(items.size());
            bis.reserve(items.size());
            for (const auto& item : items) {
                nbs.push_back(item.first);
                bis.push_back(item.second);
            }
        }
    }

    void reachableThroughBond(const std::vector<uint64_t>& srcDomain, int queryBondIndex,
                              std::vector<uint64_t>& reachable) const {
        std::fill(reachable.begin(), reachable.end(), uint64_t(0));
        for (int w = 0; w < tWords_; ++w) {
            uint64_t bits = srcDomain[w];
            while (bits != 0) {
                int bit = smsd::ctz64(bits);
                int tj = (w << 6) | bit;
                bits &= bits - 1;
                if (tj >= nT_) continue;
                for (int tn : target_.neighbors[tj]) {
                    if (q_.bondCompat(queryBondIndex, target_, tj, tn)) {
                        reachable[tn >> 6] |= uint64_t(1) << (tn & 63);
                    }
                }
            }
        }
    }

    std::vector<int> buildOrder(int anchor, int anchorTarget = -1) const {
        if (nQ_ == 0) return {};

        std::vector<int> order;
        order.reserve(nQ_);
        std::vector<bool> picked(nQ_, false);
        std::vector<std::vector<uint64_t>> workDomain = domain_;
        if (anchor >= 0 && anchorTarget >= 0 && anchorTarget < nT_) {
            std::fill(workDomain[anchor].begin(), workDomain[anchor].end(), uint64_t(0));
            workDomain[anchor][anchorTarget >> 6] |= uint64_t(1) << (anchorTarget & 63);
        }
        std::vector<uint64_t> reachable(tWords_, 0);
        int activeComp = -1;

        auto countPickedNeighbors = [&](int qi) {
            int cnt = 0;
            for (int qk : q_.neighbors[qi]) {
                if (picked[qk]) ++cnt;
            }
            return cnt;
        };

        auto chooseBest = [&](int comp, bool frontierOnly) {
            int bestQ = -1;
            int bestSize = 0x7fffffff;
            int bestPickedNbrs = -1;
            int bestDegree = -1;
            for (int qi = 0; qi < nQ_; ++qi) {
                if (picked[qi]) continue;
                if (comp >= 0 && compId_[qi] != comp) continue;
                int pickedNbrs = countPickedNeighbors(qi);
                if (frontierOnly && pickedNbrs == 0) continue;
                int size = popcountWords(workDomain[qi]);
                int degree = static_cast<int>(q_.neighbors[qi].size());
                if (bestQ < 0 ||
                    size < bestSize ||
                    (size == bestSize && pickedNbrs > bestPickedNbrs) ||
                    (size == bestSize && pickedNbrs == bestPickedNbrs && degree > bestDegree) ||
                    (size == bestSize && pickedNbrs == bestPickedNbrs && degree == bestDegree && qi < bestQ)) {
                    bestQ = qi;
                    bestSize = size;
                    bestPickedNbrs = pickedNbrs;
                    bestDegree = degree;
                }
            }
            return bestQ;
        };

        auto forwardCheck = [&](int qi) {
            for (size_t ni = 0; ni < q_.neighbors[qi].size(); ++ni) {
                int qk = q_.neighbors[qi][ni];
                if (picked[qk]) continue;
                int bi = q_.bondIndices[qi][ni];
                reachableThroughBond(workDomain[qi], bi, reachable);
                for (int w = 0; w < tWords_; ++w) {
                    workDomain[qk][w] &= reachable[w];
                }
            }
        };

        if (anchor >= 0) {
            order.push_back(anchor);
            picked[anchor] = true;
            activeComp = compId_[anchor];
            forwardCheck(anchor);
        }

        while (static_cast<int>(order.size()) < nQ_) {
            int next = -1;
            if (activeComp >= 0) {
                next = chooseBest(activeComp, true);
                if (next < 0) next = chooseBest(activeComp, false);
            }
            if (next < 0) {
                next = chooseBest(-1, false);
                if (next >= 0) activeComp = compId_[next];
            }
            if (next < 0) break;
            order.push_back(next);
            picked[next] = true;
            forwardCheck(next);
        }

        return order;
    }

    int selectCandidates(int qi, std::vector<int>& out) const {
        bool seeded = false;
        for (size_t ni = 0; ni < q_.neighbors[qi].size(); ++ni) {
            int qk = q_.neighbors[qi][ni];
            int tk = q2t_[qk];
            if (tk < 0) continue;
            int bi = q_.bondIndices[qi][ni];
            std::fill(tmpMaskBuf_.begin(), tmpMaskBuf_.end(), uint64_t(0));
            for (int tn : target_.neighbors[tk]) {
                if (!domainHas(qi, tn)) continue;
                if (!q_.bondCompat(bi, target_, tk, tn)) continue;
                tmpMaskBuf_[tn >> 6] |= uint64_t(1) << (tn & 63);
            }
            if (!seeded) {
                availBuf_ = tmpMaskBuf_;
                seeded = true;
            } else {
                for (int w = 0; w < tWords_; ++w) availBuf_[w] &= tmpMaskBuf_[w];
            }
        }

        if (!seeded) {
            availBuf_ = domain_[qi];
        }
        for (int w = 0; w < tWords_; ++w) availBuf_[w] &= ~usedMask_[w];

        collectBits(availBuf_, out);

        const int qDeg = static_cast<int>(q_.neighbors[qi].size());
        for (size_t i = 1; i < out.size(); ++i) {
            int key = out[i];
            int keyDist = target_.degree[key] - qDeg;
            if (keyDist < 0) keyDist = -keyDist;
            size_t j = i;
            while (j > 0) {
                int prev = out[j - 1];
                int prevDist = target_.degree[prev] - qDeg;
                if (prevDist < 0) prevDist = -prevDist;
                if (prevDist <= keyDist) break;
                out[j] = prev;
                --j;
            }
            out[j] = key;
        }

        return static_cast<int>(out.size());
    }

    bool feasible(int qi, int ti) const {
        if (ti < 0 || ti >= nT_) return false;
        if (t2q_[ti] != -1) return false;
        if (!domainHas(qi, ti)) return false;
        if (target_.degree[ti] < static_cast<int>(q_.neighbors[qi].size())) return false;

        int unmappedNbrs = 0;
        int availableTargetNbrs = 0;
        for (int tn : target_.neighbors[ti]) {
            if (t2q_[tn] == -1) ++availableTargetNbrs;
        }

        for (size_t ni = 0; ni < q_.neighbors[qi].size(); ++ni) {
            int qk = q_.neighbors[qi][ni];
            int bi = q_.bondIndices[qi][ni];
            int tk = q2t_[qk];
            if (tk >= 0) {
                if (!target_.hasBond(ti, tk)) return false;
                if (!q_.bondCompat(bi, target_, ti, tk)) return false;
            } else {
                ++unmappedNbrs;
            }
        }

        if (availableTargetNbrs < unmappedNbrs) return false;

        for (size_t pi = 0; pi < qNbrsByTightness_[qi].size(); ++pi) {
            int qk = qNbrsByTightness_[qi][pi];
            int tk = q2t_[qk];
            if (tk >= 0) continue;
            int bondIndex = qNbrBondByTightness_[qi][pi];
            bool supported = false;
            for (int tn : target_.neighbors[ti]) {
                if (t2q_[tn] != -1) continue;
                if (!domainHas(qk, tn)) continue;
                if (!q_.bondCompat(bondIndex, target_, ti, tn)) continue;
                supported = true;
                break;
            }
            if (!supported) return false;
        }

        return true;
    }

    void mapPair(int qi, int ti) {
        q2t_[qi] = ti;
        t2q_[ti] = qi;
        usedMask_[ti >> 6] |= uint64_t(1) << (ti & 63);
    }

    void unmapPair(int qi, int ti) {
        q2t_[qi] = -1;
        t2q_[ti] = -1;
        usedMask_[ti >> 6] &= ~(uint64_t(1) << (ti & 63));
    }

    bool backtrackExists(int pos) {
        if (pos == nQ_) return true;
        int qi = order_[pos];
        auto& buf = candBuf_[pos];
        selectCandidates(qi, buf);
        for (int ti : buf) {
            if (!feasible(qi, ti)) continue;
            mapPair(qi, ti);
            if (backtrackExists(pos + 1)) return true;
            unmapPair(qi, ti);
        }
        return false;
    }

    void backtrackAll(int pos, std::vector<std::map<int,int>>& out, int maxMatches) {
        if (static_cast<int>(out.size()) >= maxMatches) return;
        if (pos == nQ_) {
            std::map<int,int> mapping;
            for (int qi = 0; qi < nQ_; ++qi) {
                mapping[qi] = q2t_[qi];
            }
            out.push_back(std::move(mapping));
            return;
        }
        int qi = order_[pos];
        auto& buf = candBuf_[pos];
        selectCandidates(qi, buf);
        for (int ti : buf) {
            if (!feasible(qi, ti)) continue;
            mapPair(qi, ti);
            backtrackAll(pos + 1, out, maxMatches);
            unmapPair(qi, ti);
            if (static_cast<int>(out.size()) >= maxMatches) return;
        }
    }

public:
    SmartsMatchEngine(const SmartsQuery& q, const MolGraph& target, int anchorQ = -1)
        : q_(q), target_(target), nQ_(q.atomCount()), nT_(target.n),
          tWords_((target.n + 63) >> 6), anchorQ_(anchorQ),
          domain_(q.atomCount(), std::vector<uint64_t>(((target.n + 63) >> 6), 0)),
          domainSize_(q.atomCount(), 0),
          atomMemo_(static_cast<size_t>(std::max(1, q.atomCount())) * static_cast<size_t>(std::max(1, target.n)), int8_t(-1)),
          q2t_(q.atomCount(), -1), t2q_(target.n, -1),
          usedMask_(((target.n + 63) >> 6), 0),
          candBuf_(std::max(1, q.atomCount() + 1)),
          qNbrsByTightness_(q.atomCount()), qNbrBondByTightness_(q.atomCount()),
          availBuf_(((target.n + 63) >> 6), 0),
          tmpMaskBuf_(((target.n + 63) >> 6), 0) {
        if (nQ_ == 0 || nT_ == 0) return;
        buildDomains();
        if (anyEmptyDomain_) return;
        buildComponents();
        buildNeighborOrder();
        order_ = buildOrder(anchorQ_);
        if (anchorQ_ >= 0 && !order_.empty()) assert(order_[0] == anchorQ_);
    }

    bool exists() {
        if (nQ_ == 0) return true;
        if (nT_ < nQ_ || anyEmptyDomain_ || order_.size() != static_cast<size_t>(nQ_)) return false;
        return backtrackExists(0);
    }

    bool existsAnchored(int targetAtom) {
        if (nQ_ == 0) return true;
        if (nT_ < nQ_ || anyEmptyDomain_ || targetAtom < 0 || targetAtom >= nT_) return false;
        if (anchorQ_ < 0) return false;
        if (!domainHas(anchorQ_, targetAtom)) return false;

        std::vector<int> savedOrder = order_;
        order_ = buildOrder(anchorQ_, targetAtom);
        if (order_.empty() || order_[0] != anchorQ_) {
            order_ = std::move(savedOrder);
            return false;
        }

        if (!feasible(anchorQ_, targetAtom)) {
            order_ = std::move(savedOrder);
            return false;
        }

        mapPair(anchorQ_, targetAtom);
        bool ok = backtrackExists(1);
        unmapPair(anchorQ_, targetAtom);
        order_ = std::move(savedOrder);
        return ok;
    }

    std::vector<std::map<int,int>> findAll(int maxMatches) {
        std::vector<std::map<int,int>> out;
        if (maxMatches <= 0) return out;
        if (nQ_ == 0) {
            out.push_back({});
            return out;
        }
        if (nT_ < nQ_ || anyEmptyDomain_ || order_.size() != static_cast<size_t>(nQ_)) return out;
        backtrackAll(0, out, maxMatches);
        return out;
    }
};

} // namespace smarts_detail



enum class SmartsProbeResult {
    FOUND,
    NOT_FOUND,
    UNKNOWN
};

inline SmartsProbeResult legacyBacktrackProbe(
        const SmartsQuery& query,
        const MolGraph& target,
        const std::vector<int>& order,
        int pos,
        std::vector<int>& q2t,
        std::vector<bool>& tUsed,
        std::vector<int>& mark,
        int& markStamp,
        int& nodesVisited,
        int nodeLimit) {
    if (nodesVisited++ >= nodeLimit) return SmartsProbeResult::UNKNOWN;

    int qn = query.atomCount();
    if (pos == qn) return SmartsProbeResult::FOUND;

    int qi = order[pos];
    std::vector<int> candidates;
    bool constrained = false;

    auto nextStamp = [&]() {
        if (++markStamp == std::numeric_limits<int>::max()) {
            std::fill(mark.begin(), mark.end(), 0);
            markStamp = 1;
        }
        return markStamp;
    };

    for (size_t ni = 0; ni < query.neighbors[qi].size(); ++ni) {
        int qk = query.neighbors[qi][ni];
        int tk = q2t[qk];
        if (tk < 0) continue;

        int bi = query.bondIndices[qi][ni];
        if (!constrained) {
            constrained = true;
            for (int tn : target.neighbors[tk]) {
                if (tUsed[tn]) continue;
                if (!query.atomCompat(qi, target, tn)) continue;
                if (!query.bondCompat(bi, target, tk, tn)) continue;
                candidates.push_back(tn);
            }
        } else {
            int stamp = nextStamp();
            for (int tn : target.neighbors[tk]) {
                if (!query.bondCompat(bi, target, tk, tn)) continue;
                mark[tn] = stamp;
            }
            candidates.erase(
                std::remove_if(candidates.begin(), candidates.end(),
                    [&](int c) { return mark[c] != stamp; }),
                candidates.end());
        }

        if (candidates.empty()) return SmartsProbeResult::NOT_FOUND;
    }

    if (!constrained) {
        candidates.reserve(target.n);
        for (int t = 0; t < target.n; ++t) {
            if (!tUsed[t] && query.atomCompat(qi, target, t)) {
                candidates.push_back(t);
            }
        }
    }

    for (int ti : candidates) {
        if (tUsed[ti]) continue;

        bool bondsOk = true;
        for (size_t ni = 0; ni < query.neighbors[qi].size(); ++ni) {
            int qk = query.neighbors[qi][ni];
            if (q2t[qk] >= 0) {
                int tk = q2t[qk];
                int bi = query.bondIndices[qi][ni];
                if (!target.hasBond(ti, tk) || !query.bondCompat(bi, target, ti, tk)) {
                    bondsOk = false;
                    break;
                }
            }
        }
        if (!bondsOk) continue;

        q2t[qi] = ti;
        tUsed[ti] = true;

        SmartsProbeResult sub = legacyBacktrackProbe(
            query, target, order, pos + 1, q2t, tUsed, mark, markStamp, nodesVisited, nodeLimit);

        q2t[qi] = -1;
        tUsed[ti] = false;

        if (sub == SmartsProbeResult::FOUND || sub == SmartsProbeResult::UNKNOWN) {
            return sub;
        }
    }

    return SmartsProbeResult::NOT_FOUND;
}

inline SmartsProbeResult legacyExistsProbe(
        const SmartsQuery& query,
        const MolGraph& target,
        int nodeLimit,
        bool anchored = false,
        int targetAtom = -1) {
    if (query.atomCount() == 0) return SmartsProbeResult::FOUND;
    if (target.n < query.atomCount()) return SmartsProbeResult::NOT_FOUND;

    auto order = query.matchingOrder();
    std::vector<int> q2t(query.atomCount(), -1);
    std::vector<bool> tUsed(target.n, false);
    std::vector<int> mark(std::max(1, target.n), 0);
    int markStamp = 1;
    int nodesVisited = 0;

    if (anchored) {
        if (order.empty() || order[0] != 0) return SmartsProbeResult::UNKNOWN;
        if (targetAtom < 0 || targetAtom >= target.n) return SmartsProbeResult::NOT_FOUND;
        if (!query.atomCompat(0, target, targetAtom)) return SmartsProbeResult::NOT_FOUND;
        q2t[0] = targetAtom;
        tUsed[targetAtom] = true;
        return legacyBacktrackProbe(query, target, order, 1, q2t, tUsed, mark, markStamp, nodesVisited, nodeLimit);
    }

    return legacyBacktrackProbe(query, target, order, 0, q2t, tUsed, mark, markStamp, nodesVisited, nodeLimit);
}

inline bool SmartsQuery::matches(const MolGraph& target) const {
    if (atomCount() == 0) return true;
    if (target.n < atomCount()) return false;

    constexpr int kProbeNodeLimit = 2048;
    SmartsProbeResult probe = legacyExistsProbe(*this, target, kProbeNodeLimit, false, -1);
    if (probe == SmartsProbeResult::FOUND) return true;
    if (probe == SmartsProbeResult::NOT_FOUND) return false;

    smarts_detail::SmartsMatchEngine engine(*this, target, -1);
    return engine.exists();
}

inline std::vector<std::map<int,int>> SmartsQuery::findAll(const MolGraph& target, int maxMatches) const {
    if (atomCount() == 0) return {{}};
    if (target.n < atomCount()) return {};
    smarts_detail::SmartsMatchEngine engine(*this, target, -1);
    return engine.findAll(maxMatches);
}

inline bool SmartsQuery::matchesAtAtom(const MolGraph& target, int targetAtom) const {
    if (atomCount() == 0) return true;
    if (target.n < atomCount()) return false;

    constexpr int kProbeNodeLimit = 1024;
    SmartsProbeResult probe = legacyExistsProbe(*this, target, kProbeNodeLimit, true, targetAtom);
    if (probe == SmartsProbeResult::FOUND) return true;
    if (probe == SmartsProbeResult::NOT_FOUND) return false;

    smarts_detail::SmartsMatchEngine engine(*this, target, 0);
    return engine.existsAnchored(targetAtom);
}

// Deferred recursive SMARTS evaluation (needs complete SmartsQuery type)
inline bool smarts_detail::evalRecursiveSmarts(const smarts_detail::AtomExpr& expr, const MolGraph& g, int idx) {
    if (expr.recursiveSmarts.empty()) return false;

    if (!expr.cachedRecursive) {
        expr.cachedRecursive = std::make_shared<SmartsQuery>(parseSMARTS(expr.recursiveSmarts));
    }

    if (expr.cachedRecursiveTarget != &g
        || expr.cachedRecursiveTargetNonce != g.cacheNonce
        || static_cast<int>(expr.cachedRecursiveState.size()) != g.n) {
        expr.cachedRecursiveTarget = &g;
        expr.cachedRecursiveTargetNonce = g.cacheNonce;
        expr.cachedRecursiveState.assign(g.n, int8_t(-1));
    }

    int8_t& memo = expr.cachedRecursiveState[idx];
    if (memo >= 0) return memo != 0;

    memo = expr.cachedRecursive->matchesAtAtom(g, idx) ? int8_t(1) : int8_t(0);
    return memo != 0;
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
    return query;
}

// ============================================================================
// SMARTS-based MCS: find the largest SMARTS match in a target molecule
// ============================================================================

/**
 * Find the largest substructure match of a SMARTS pattern in a target MolGraph.
 *
 * Parses the SMARTS query, enumerates all matches in the target, and returns
 * the mapping with the most matched atoms.  For simple SMARTS patterns that
 * map to connected substructures, all matches will have the same size (equal
 * to the number of query atoms), so this simply returns the first match.
 *
 * For patterns with wildcards or optional components, different embeddings
 * may cover different numbers of target atoms; this function picks the largest.
 *
 * @param smartsStr  SMARTS pattern string
 * @param target     Target molecule
 * @param maxMatches Maximum number of matches to enumerate before picking best
 * @return Mapping from SMARTS atom index to target atom index for the largest
 *         match, or empty map if no match.
 */
inline std::map<int,int> findMcsSmarts(
        const std::string& smartsStr,
        const MolGraph& target,
        int maxMatches = 1000) {
    if (smartsStr.empty() || target.n == 0) return {};

    SmartsQuery query = parseSMARTS(smartsStr);
    if (query.atomCount() == 0) return {};

    auto matches = query.findAll(target, maxMatches);
    if (matches.empty()) return {};

    // Return the largest match (most atoms mapped)
    std::map<int,int>* best = &matches[0];
    for (size_t i = 1; i < matches.size(); i++) {
        if (matches[i].size() > best->size()) {
            best = &matches[i];
        }
    }
    return *best;
}

} // namespace smsd

// ============================================================================
// Basic tests (compile with -DSMSD_TEST_SMARTS)
// ============================================================================
#ifdef SMSD_TEST_SMARTS

#include "smiles_parser.hpp"
#include <cassert>
#include <iostream>

namespace smsd_smarts_test {

inline void testWildcard() {
    auto q = smsd::parseSMARTS("*");
    auto mol = smsd::parseSMILES("C");
    assert(q.matches(mol));
    std::cout << "  [PASS] testWildcard\n";
}

inline void testElementMatch() {
    auto q = smsd::parseSMARTS("[#6]");
    auto mol = smsd::parseSMILES("C");
    assert(q.matches(mol));

    auto q2 = smsd::parseSMARTS("[#7]");
    assert(!q2.matches(mol));
    std::cout << "  [PASS] testElementMatch\n";
}

inline void testOrganicShorthand() {
    auto q = smsd::parseSMARTS("C");
    auto mol = smsd::parseSMILES("CC");
    assert(q.matches(mol));

    auto q2 = smsd::parseSMARTS("N");
    assert(!q2.matches(mol));
    std::cout << "  [PASS] testOrganicShorthand\n";
}

inline void testBondSingle() {
    auto q = smsd::parseSMARTS("C-C");
    auto mol = smsd::parseSMILES("CC");
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("C=C");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testBondSingle\n";
}

inline void testBondDouble() {
    auto q = smsd::parseSMARTS("C=O");
    auto mol = smsd::parseSMILES("CC=O");
    assert(q.matches(mol));
    std::cout << "  [PASS] testBondDouble\n";
}

inline void testBondTriple() {
    auto q = smsd::parseSMARTS("C#N");
    auto mol = smsd::parseSMILES("CC#N");
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("CC=N");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testBondTriple\n";
}

inline void testBondAny() {
    auto q = smsd::parseSMARTS("C~C");
    auto mol1 = smsd::parseSMILES("CC");
    auto mol2 = smsd::parseSMILES("C=C");
    assert(q.matches(mol1));
    assert(q.matches(mol2));
    std::cout << "  [PASS] testBondAny\n";
}

inline void testDegree() {
    // D2 = degree 2
    auto q = smsd::parseSMARTS("[D2]");
    auto mol = smsd::parseSMILES("CCC");  // middle C has degree 2
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("C");   // single C has degree 0
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testDegree\n";
}

inline void testRing() {
    auto q = smsd::parseSMARTS("[R]");
    auto mol = smsd::parseSMILES("C1CCC1");
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("CCC");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testRing\n";
}

inline void testNotRing() {
    auto q = smsd::parseSMARTS("[!R]");
    auto mol = smsd::parseSMILES("CCC");
    assert(q.matches(mol));
    std::cout << "  [PASS] testNotRing\n";
}

inline void testOrExpression() {
    // [#6,#7] = carbon or nitrogen
    auto q = smsd::parseSMARTS("[#6,#7]");
    auto molC = smsd::parseSMILES("C");
    auto molN = smsd::parseSMILES("N");
    auto molO = smsd::parseSMILES("O");
    assert(q.matches(molC));
    assert(q.matches(molN));
    assert(!q.matches(molO));
    std::cout << "  [PASS] testOrExpression\n";
}

inline void testAndExpression() {
    // [#6&D2] = carbon with degree 2
    auto q = smsd::parseSMARTS("[#6&D2]");
    auto mol = smsd::parseSMILES("CCC");
    assert(q.matches(mol));  // middle C is #6 and D2
    std::cout << "  [PASS] testAndExpression\n";
}

inline void testCharge() {
    auto q = smsd::parseSMARTS("[+1]");
    auto mol = smsd::parseSMILES("[NH4+]");
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("C");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testCharge\n";
}

inline void testTwoLetterElement() {
    auto q = smsd::parseSMARTS("[Cl]");
    auto mol = smsd::parseSMILES("CCl");
    assert(q.matches(mol));
    std::cout << "  [PASS] testTwoLetterElement\n";
}

inline void testPattern_CC() {
    auto q = smsd::parseSMARTS("CC");
    auto mol = smsd::parseSMILES("CCC");  // has C-C substructure
    auto matches = q.findAll(mol, 100);
    assert(!matches.empty());
    std::cout << "  [PASS] testPattern_CC (" << matches.size() << " matches)\n";
}

inline void testAromatic() {
    auto q = smsd::parseSMARTS("[a]");
    auto mol = smsd::parseSMILES("c1ccccc1");  // benzene
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("C1CCCCC1");  // cyclohexane
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testAromatic\n";
}

inline void testRingBond() {
    auto q = smsd::parseSMARTS("[#6]@[#6]");
    auto mol = smsd::parseSMILES("C1CCC1");  // cyclobutane: ring bond
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("CCC");  // no ring
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testRingBond\n";
}

inline void testBranches() {
    // C(=O)O = carboxylic acid
    auto q = smsd::parseSMARTS("C(=O)O");
    auto mol = smsd::parseSMILES("CC(=O)O");  // acetic acid
    assert(q.matches(mol));
    std::cout << "  [PASS] testBranches (carboxylic acid)\n";
}

inline void testNotOperator() {
    // [!#6] = not carbon
    auto q = smsd::parseSMARTS("[!#6]");
    auto molC = smsd::parseSMILES("C");
    auto molN = smsd::parseSMILES("N");
    assert(!q.matches(molC));
    assert(q.matches(molN));
    std::cout << "  [PASS] testNotOperator\n";
}

inline void testRingClosure() {
    // c1ccccc1 = aromatic 6-ring
    auto q = smsd::parseSMARTS("c1ccccc1");
    auto mol = smsd::parseSMILES("c1ccccc1");
    assert(q.matches(mol));
    std::cout << "  [PASS] testRingClosure\n";
}

inline void testEmptySmarts() {
    auto q = smsd::parseSMARTS("");
    auto mol = smsd::parseSMILES("C");
    assert(q.matches(mol));  // empty pattern matches anything
    std::cout << "  [PASS] testEmptySmarts\n";
}

inline void testInvalidSmarts() {
    bool threw = false;
    try {
        smsd::parseSMARTS("[#6");  // missing ]
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
    std::cout << "  [PASS] testInvalidSmarts\n";
}

inline void testSemicolonAnd() {
    // [#6;R] = carbon AND in ring (semicolon = low-priority AND)
    auto q = smsd::parseSMARTS("[#6;R]");
    auto mol = smsd::parseSMILES("C1CCC1");
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("CCC");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testSemicolonAnd\n";
}

inline void testMultipleMatches() {
    auto q = smsd::parseSMARTS("[#6]");
    auto mol = smsd::parseSMILES("CCCC");
    auto matches = q.findAll(mol, 100);
    assert(static_cast<int>(matches.size()) == 4);  // 4 carbons
    std::cout << "  [PASS] testMultipleMatches (" << matches.size() << " matches)\n";
}

inline void testRecursiveSmarts() {
    // [$([#6][#8])] = carbon bonded to oxygen
    auto q = smsd::parseSMARTS("[$([#6][#8])]");
    auto mol = smsd::parseSMILES("CCO");
    assert(q.matches(mol));  // C bonded to O

    auto mol2 = smsd::parseSMILES("CCC");
    assert(!q.matches(mol2));
    std::cout << "  [PASS] testRecursiveSmarts\n";
}

inline void testHeteroNeighbors() {
    // z1 in methylamine: C has 1 N neighbor
    auto q1 = smsd::parseSMARTS("[z1]");
    auto mol1 = smsd::parseSMILES("CN");
    assert(q1.matches(mol1));

    // z0 in ethane: only carbon neighbors
    auto q2 = smsd::parseSMARTS("[z0]");
    auto mol2 = smsd::parseSMILES("CC");
    assert(q2.matches(mol2));

    std::cout << "  [PASS] testHeteroNeighbors\n";
}

inline void testHeavyDegree() {
    // d2 in propane: middle C has 2 heavy neighbors
    auto q = smsd::parseSMARTS("[d2]");
    auto mol = smsd::parseSMILES("CCC");
    assert(q.matches(mol));
    std::cout << "  [PASS] testHeavyDegree\n";
}

inline void testHybridization() {
    // ^2 matches SP2 (aromatic carbon)
    auto q1 = smsd::parseSMARTS("[^2]");
    auto mol1 = smsd::parseSMILES("c1ccccc1");
    assert(q1.matches(mol1));

    // ^3 matches SP3 (methane)
    auto q2 = smsd::parseSMARTS("[^3]");
    auto mol2 = smsd::parseSMILES("C");
    assert(q2.matches(mol2));

    // ^1 matches SP (acetylene)
    auto q3 = smsd::parseSMARTS("[^1]");
    auto mol3 = smsd::parseSMILES("C#C");
    assert(q3.matches(mol3));

    std::cout << "  [PASS] testHybridization\n";
}

inline void testRangeQuery() {
    // D{2-4} matches degree 2 through 4
    auto q = smsd::parseSMARTS("[D{2-4}]");
    auto mol = smsd::parseSMILES("CC(C)(C)C"); // central C has D4
    assert(q.matches(mol));

    auto mol2 = smsd::parseSMILES("C"); // D0
    assert(!q.matches(mol2));

    // D{2-} means >= 2
    auto q2 = smsd::parseSMARTS("[D{2-}]");
    auto mol3 = smsd::parseSMILES("CCC"); // middle C has D2
    assert(q2.matches(mol3));

    std::cout << "  [PASS] testRangeQuery\n";
}

inline void runAllTests() {
    std::cout << "smarts_parser tests:\n";
    testWildcard();
    testElementMatch();
    testOrganicShorthand();
    testBondSingle();
    testBondDouble();
    testBondTriple();
    testBondAny();
    testDegree();
    testRing();
    testNotRing();
    testOrExpression();
    testAndExpression();
    testCharge();
    testTwoLetterElement();
    testPattern_CC();
    testAromatic();
    testRingBond();
    testBranches();
    testNotOperator();
    testRingClosure();
    testEmptySmarts();
    testInvalidSmarts();
    testSemicolonAnd();
    testMultipleMatches();
    testRecursiveSmarts();
    testHeteroNeighbors();
    testHeavyDegree();
    testHybridization();
    testRangeQuery();
    std::cout << "All smarts_parser tests passed.\n";
}

} // namespace smsd_smarts_test

#endif // SMSD_TEST_SMARTS

#endif // SMSD_SMARTS_PARSER_HPP

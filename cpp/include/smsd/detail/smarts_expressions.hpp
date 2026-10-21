/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMARTS expressions: the atom and bond expression trees, query atoms and bonds,
 * the hooks through which recursive SMARTS are matched, and element symbols.
 *
 * Part of smarts_parser.hpp.
 */
#pragma once

#include "smsd/mol_graph.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace smsd {

// ============================================================================
// Forward declarations
// ============================================================================
struct SmartsQuery;
SmartsQuery parseSMARTS(const std::string& smarts, int maxRecursionDepth = 20);

// ============================================================================
// Recursive matcher hooks and the element table (from smiles_utils.hpp)
// ============================================================================
namespace smarts_detail {

using RecursiveMatcher = bool (*)(const SmartsQuery&, const MolGraph&, int);
inline thread_local RecursiveMatcher recursiveMatcher = nullptr;
/** Non-zero while a recursive query parsed for one evaluation is matched; its address is not reusable. */
inline thread_local int transientRecursionDepth = 0;

inline int elementToAtomicNum(const std::string& sym) {
    return detail::elementToAtomicNum(sym);
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
    HYBRIDIZATION,             // ^<n> — hybridisation (0=S,1=SP,2=SP2,3=SP3,4=SP3D,5=SP3D2)
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
    // For ELEMENT: the required aromaticity (0 aliphatic, 1 aromatic, -1 either), tested
    // with the element as one primitive rather than as a separate AND operand.
    int aromaticity = -1;
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

    /** An element written as a symbol: aliphatic when uppercase, aromatic when lowercase. */
    static AtomExpr makeElement(int atomicNumber, bool aromatic) {
        AtomExpr e = makePrim(AtomPrimType::ELEMENT, atomicNumber, true);
        e.aromaticity = aromatic ? 1 : 0;
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
    WEDGE_OR_UNSPECIFIED, // /? (up, or no configuration)
    DASH_OR_UNSPECIFIED,  // \? (down, or no configuration)
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

} // namespace smarts_detail

} // namespace smsd

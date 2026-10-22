/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMARTS evaluation: atom and bond expressions against a target, with three-valued
 * stereo results for chirality and bond direction.
 *
 * Part of smarts_parser.hpp.
 */
#pragma once

#include "smsd/stereo.hpp"
#include "smsd/detail/smarts_properties.hpp"

namespace smsd {
namespace smarts_detail {

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
            return g.atomicNum[idx] == expr.intVal
                && (expr.aromaticity < 0 || (g.aromatic[idx] != 0) == (expr.aromaticity != 0));

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
            // 1 = @ (anticlockwise), 2 = @@ (clockwise); @? and @@? also accept no configuration
            if (g.tetraChirality.empty() || g.tetraChirality[idx] == 0) return expr.orUnspecified;
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
        case BondPrimType::WEDGE_OR_UNSPECIFIED:
        case BondPrimType::DASH_OR_UNSPECIFIED:
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

/** Whether a bond primitive is directional: /, \, /? or \?. */
inline bool isDirectionalPrimitive(BondPrimType type) {
    return type == BondPrimType::WEDGE || type == BondPrimType::DASH
        || type == BondPrimType::WEDGE_OR_UNSPECIFIED || type == BondPrimType::DASH_OR_UNSPECIFIED;
}

/** Preserve Boolean expressions while mapped chirality is not yet known. */
enum class StereoTruth { FALSE_VALUE, TRUE_VALUE, UNKNOWN };

/**
 * Whether a bond expression accepts a bond that carries no direction, as one next to a double
 * bond without a configuration does: /? and \? do, / and \ leave it open, and other primitives
 * read the bond as usual. An expression left open does not accept it, so / and \ alone, or
 * together as /,\, need a configuration, while /,- does not.
 */
inline StereoTruth evalBondExprWithoutDirection(const BondExpr& expression, const MolGraph& graph,
                                                int first, int second) {
    if (expression.op == ExprOp::PRIM) {
        if (!isDirectionalPrimitive(expression.primType))
            return evalBondPrim(expression, graph, first, second) ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
        if (graph.bondOrder(first, second) != 1 || graph.bondAromatic(first, second)) return StereoTruth::FALSE_VALUE;
        return expression.primType == BondPrimType::WEDGE_OR_UNSPECIFIED
                || expression.primType == BondPrimType::DASH_OR_UNSPECIFIED
            ? StereoTruth::TRUE_VALUE : StereoTruth::UNKNOWN;
    }
    if (expression.op == ExprOp::NOT) {
        const auto result = evalBondExprWithoutDirection(expression.children[0], graph, first, second);
        return result == StereoTruth::UNKNOWN ? result
            : (result == StereoTruth::TRUE_VALUE ? StereoTruth::FALSE_VALUE : StereoTruth::TRUE_VALUE);
    }
    bool unknown = false;
    for (const auto& child : expression.children) {
        const auto result = evalBondExprWithoutDirection(child, graph, first, second);
        if (expression.op == ExprOp::AND && result == StereoTruth::FALSE_VALUE) return result;
        if (expression.op == ExprOp::OR && result == StereoTruth::TRUE_VALUE) return result;
        if (result == StereoTruth::UNKNOWN) unknown = true;
    }
    if (unknown) return StereoTruth::UNKNOWN;
    return expression.op == ExprOp::AND ? StereoTruth::TRUE_VALUE : StereoTruth::FALSE_VALUE;
}

template <typename ChiralityEvaluator>
inline StereoTruth evalAtomExprWithStereo(const AtomExpr& expression, const MolGraph& graph,
                                          int atom, const ChiralityEvaluator& evaluateChirality) {
    if (expression.op == ExprOp::PRIM) {
        if (expression.primType == AtomPrimType::CHIRALITY) {
            // @? and @@? accept an atom without a configuration whatever the mapping.
            if (expression.orUnspecified && (graph.tetraChirality.empty() || graph.tetraChirality[atom] == 0))
                return StereoTruth::TRUE_VALUE;
            return evaluateChirality(expression);
        }
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
    if (expression.op == ExprOp::PRIM) return isDirectionalPrimitive(expression.primType);
    for (const auto& child : expression.children) if (hasDirectionalPrimitive(child)) return true;
    return false;
}

/** Evaluate a directional bond expression for one candidate slash direction, 1 up or -1 down. */
inline bool evalBondExprWithDirection(const BondExpr& expression, const MolGraph& graph,
                                      int first, int second, int direction) {
    if (expression.op == ExprOp::PRIM) {
        if (isDirectionalPrimitive(expression.primType)) {
            const bool up = expression.primType == BondPrimType::WEDGE
                || expression.primType == BondPrimType::WEDGE_OR_UNSPECIFIED;
            return graph.bondOrder(first, second) == 1 && !graph.bondAromatic(first, second)
                && direction == (up ? 1 : -1);
        }
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

} // namespace smarts_detail
} // namespace smsd

// evalRecursiveSmarts() needs the complete SmartsQuery, so it is defined there.
#include "smsd/detail/smarts_query.hpp"

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. *
 *
 * Lightweight depiction-facing API built on top of SMSD's layout engine.
 * This provides RDDepictor-style constrained coordinate generation without
 * pulling depiction policy into the renderer itself.
 */
#pragma once
#ifndef SMSD_DEPICTOR_HPP
#define SMSD_DEPICTOR_HPP

#include "smsd/layout.hpp"

#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace smsd {

using CoordMap2D = std::map<int, Point2D>;
using MatchVect2D = std::vector<std::pair<int, int>>; // (referenceAtom, moleculeAtom)

struct DepictorOptions {
    double bondLength = 1.5;
    int constrainedRefinementIters = 120;
    double constrainedStep = 0.22;
    double repulsionStrength = 0.015;
    bool kekulize = true;
};

struct ConstrainedDepictionParams {
    bool acceptFailure = false;
    bool alignOnly = false;
    bool useRingTemplates = true; // reserved for future layout backends
    int refinementIters = 120;
    double bondLength = 1.5;
    bool kekulize = true;
};

namespace detail_depictor {

inline double sqr(double v) { return v * v; }

inline double inferBondLengthFromAnchors(const MolGraph& mol,
                                         const CoordMap2D& coordMap,
                                         double fallback) {
    double sum = 0.0;
    int count = 0;
    for (const auto& anchor : coordMap) {
        int a = anchor.first;
        for (int b : mol.neighbors[a]) {
            auto it = coordMap.find(b);
            if (it == coordMap.end() || b <= a) {
                continue;
            }
            double dx = anchor.second.x - it->second.x;
            double dy = anchor.second.y - it->second.y;
            double dist = std::sqrt(dx * dx + dy * dy);
            if (dist > 1e-6) {
                sum += dist;
                count++;
            }
        }
    }
    return count > 0 ? sum / static_cast<double>(count) : fallback;
}

inline void applyRigidAlignment(std::vector<Point2D>& coords,
                                const std::vector<int>& mobileIndices,
                                const std::vector<Point2D>& reference) {
    size_t n = std::min(mobileIndices.size(), reference.size());
    if (n == 0 || coords.empty()) {
        return;
    }

    double cx1 = 0.0, cy1 = 0.0, cx2 = 0.0, cy2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        const auto& src = coords[mobileIndices[i]];
        cx1 += src.x;
        cy1 += src.y;
        cx2 += reference[i].x;
        cy2 += reference[i].y;
    }
    cx1 /= static_cast<double>(n);
    cy1 /= static_cast<double>(n);
    cx2 /= static_cast<double>(n);
    cy2 /= static_cast<double>(n);

    double theta = 0.0;
    if (n >= 2) {
        double sumCross = 0.0;
        double sumDot = 0.0;
        for (size_t i = 0; i < n; ++i) {
            const auto& src = coords[mobileIndices[i]];
            double ax = src.x - cx1;
            double ay = src.y - cy1;
            double bx = reference[i].x - cx2;
            double by = reference[i].y - cy2;
            sumDot += ax * bx + ay * by;
            sumCross += ax * by - ay * bx;
        }
        theta = std::atan2(sumCross, sumDot);
    }

    double cosT = std::cos(theta);
    double sinT = std::sin(theta);
    for (auto& point : coords) {
        double rx = point.x - cx1;
        double ry = point.y - cy1;
        point.x = cx2 + rx * cosT - ry * sinT;
        point.y = cy2 + rx * sinT + ry * cosT;
    }
}

inline void enforceAnchors(std::vector<Point2D>& coords, const CoordMap2D& coordMap) {
    for (const auto& entry : coordMap) {
        if (entry.first >= 0 && entry.first < static_cast<int>(coords.size())) {
            coords[entry.first] = entry.second;
        }
    }
}

inline void refineFreeAtoms(const MolGraph& mol,
                            std::vector<Point2D>& coords,
                            const std::vector<bool>& fixed,
                            double targetBondLength,
                            int iterations,
                            double baseStep,
                            double repulsionStrength) {
    if (iterations <= 0 || mol.n <= 1) {
        enforceAnchors(coords, CoordMap2D{});
        return;
    }

    std::vector<Point2D> delta(static_cast<size_t>(mol.n), {0.0, 0.0});
    for (int iter = 0; iter < iterations; ++iter) {
        std::fill(delta.begin(), delta.end(), Point2D{0.0, 0.0});

        for (int a = 0; a < mol.n; ++a) {
            for (int b : mol.neighbors[a]) {
                if (b <= a) {
                    continue;
                }
                double dx = coords[b].x - coords[a].x;
                double dy = coords[b].y - coords[a].y;
                double dist = std::sqrt(dx * dx + dy * dy);
                if (dist < 1e-6) {
                    dx = 1e-3 * static_cast<double>((a + b) % 3 + 1);
                    dy = 1e-3 * static_cast<double>((a + 2 * b) % 3 + 1);
                    dist = std::sqrt(dx * dx + dy * dy);
                }
                double ux = dx / dist;
                double uy = dy / dist;
                double stretch = dist - targetBondLength;
                double fx = 0.18 * stretch * ux;
                double fy = 0.18 * stretch * uy;

                if (!fixed[a]) {
                    delta[static_cast<size_t>(a)].x += fx;
                    delta[static_cast<size_t>(a)].y += fy;
                }
                if (!fixed[b]) {
                    delta[static_cast<size_t>(b)].x -= fx;
                    delta[static_cast<size_t>(b)].y -= fy;
                }
            }
        }

        for (int a = 0; a < mol.n; ++a) {
            if (fixed[a]) {
                continue;
            }
            for (int b = a + 1; b < mol.n; ++b) {
                double dx = coords[b].x - coords[a].x;
                double dy = coords[b].y - coords[a].y;
                double dist2 = dx * dx + dy * dy;
                if (dist2 < 1e-6) {
                    dist2 = 1e-6;
                    dx = 1e-3;
                    dy = 0.0;
                }
                double dist = std::sqrt(dist2);
                double cutoff = targetBondLength * 1.25;
                if (dist > cutoff) {
                    continue;
                }
                double repel = repulsionStrength / dist2;
                double fx = repel * dx / dist;
                double fy = repel * dy / dist;
                if (!fixed[a]) {
                    delta[static_cast<size_t>(a)].x -= fx;
                    delta[static_cast<size_t>(a)].y -= fy;
                }
                if (!fixed[b]) {
                    delta[static_cast<size_t>(b)].x += fx;
                    delta[static_cast<size_t>(b)].y += fy;
                }
            }
        }

        double damping = baseStep * (1.0 - 0.65 * (static_cast<double>(iter) /
                                                   static_cast<double>(iterations)));
        double maxStep = targetBondLength * 0.20;
        for (int i = 0; i < mol.n; ++i) {
            if (fixed[i]) {
                continue;
            }
            double dx = delta[static_cast<size_t>(i)].x * damping;
            double dy = delta[static_cast<size_t>(i)].y * damping;
            double stepLen = std::sqrt(dx * dx + dy * dy);
            if (stepLen > maxStep && stepLen > 1e-12) {
                double scale = maxStep / stepLen;
                dx *= scale;
                dy *= scale;
            }
            coords[static_cast<size_t>(i)].x += dx;
            coords[static_cast<size_t>(i)].y += dy;
        }
    }
}

inline void validateCoordMap(const MolGraph& mol, const CoordMap2D& coordMap) {
    for (const auto& entry : coordMap) {
        if (entry.first < 0 || entry.first >= mol.n) {
            throw std::out_of_range("coordMap atom index out of range");
        }
    }
}

} // namespace detail_depictor

inline MolGraph prepareMolForDrawing(MolGraph mol, bool kekulize = true) {
    if (kekulize) {
        mol.kekulize();
    }
    return mol;
}

inline std::vector<Point2D> compute2DCoords(const MolGraph& mol,
                                            const DepictorOptions& options = {}) {
    MolGraph prepared = prepareMolForDrawing(mol, options.kekulize);
    return generateCoords2D(prepared, options.bondLength);
}

inline std::vector<Point2D> compute2DCoords(const MolGraph& mol,
                                            const CoordMap2D& coordMap,
                                            const DepictorOptions& options = {}) {
    detail_depictor::validateCoordMap(mol, coordMap);
    if (coordMap.empty()) {
        return compute2DCoords(mol, options);
    }

    MolGraph prepared = prepareMolForDrawing(mol, options.kekulize);
    double targetBondLength = detail_depictor::inferBondLengthFromAnchors(
        prepared, coordMap, options.bondLength);
    auto coords = generateCoords2D(prepared, targetBondLength);

    std::vector<int> mobileIndices;
    std::vector<Point2D> reference;
    mobileIndices.reserve(coordMap.size());
    reference.reserve(coordMap.size());
    for (const auto& entry : coordMap) {
        mobileIndices.push_back(entry.first);
        reference.push_back(entry.second);
    }
    detail_depictor::applyRigidAlignment(coords, mobileIndices, reference);

    std::vector<bool> fixed(static_cast<size_t>(prepared.n), false);
    for (const auto& entry : coordMap) {
        fixed[static_cast<size_t>(entry.first)] = true;
    }
    detail_depictor::enforceAnchors(coords, coordMap);
    detail_depictor::refineFreeAtoms(prepared, coords, fixed, targetBondLength,
                                     options.constrainedRefinementIters,
                                     options.constrainedStep,
                                     options.repulsionStrength);
    detail_depictor::enforceAnchors(coords, coordMap);
    return coords;
}

inline std::vector<Point2D> generateDepictionMatching2DStructure(
    const MolGraph& mol,
    const MolGraph& reference,
    const MatchVect2D& refMatchVect,
    const std::vector<Point2D>& referenceCoords,
    const ConstrainedDepictionParams& params = {}) {
    if (referenceCoords.size() < static_cast<size_t>(reference.n)) {
        throw std::invalid_argument("referenceCoords must cover the reference molecule");
    }
    if (refMatchVect.empty()) {
        if (!params.acceptFailure) {
            throw std::invalid_argument("constrained depiction requires at least one match");
        }
        return compute2DCoords(mol, DepictorOptions{params.bondLength,
                                                    params.refinementIters,
                                                    0.22,
                                                    0.015,
                                                    params.kekulize});
    }

    CoordMap2D coordMap;
    std::vector<int> mobileIndices;
    std::vector<Point2D> matchedReference;
    for (const auto& match : refMatchVect) {
        int refAtom = match.first;
        int molAtom = match.second;
        if (refAtom < 0 || refAtom >= reference.n || molAtom < 0 || molAtom >= mol.n) {
            throw std::out_of_range("reference match index out of range");
        }
        coordMap[molAtom] = referenceCoords[static_cast<size_t>(refAtom)];
        mobileIndices.push_back(molAtom);
        matchedReference.push_back(referenceCoords[static_cast<size_t>(refAtom)]);
    }

    DepictorOptions depictorOptions;
    depictorOptions.bondLength = params.bondLength;
    depictorOptions.constrainedRefinementIters = params.refinementIters;
    depictorOptions.kekulize = params.kekulize;

    if (params.alignOnly) {
        auto coords = compute2DCoords(mol, depictorOptions);
        detail_depictor::applyRigidAlignment(coords, mobileIndices, matchedReference);
        return coords;
    }

    return compute2DCoords(mol, coordMap, depictorOptions);
}

inline void normalizeDepiction(const MolGraph& mol,
                               std::vector<Point2D>& coords,
                               double targetBondLength = 1.5) {
    if (coords.size() < static_cast<size_t>(mol.n)) {
        coords.resize(static_cast<size_t>(mol.n));
    }
    transform::canonicalOrientation(mol, coords);
    transform::normaliseBondLength(mol, coords, targetBondLength);
}

} // namespace smsd

#endif // SMSD_DEPICTOR_HPP

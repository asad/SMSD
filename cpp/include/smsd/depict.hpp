// SMSD — Molecular Depiction Engine
// Zero-dependency 2D layout + SVG renderer with MCS highlighting
// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2018-2026 BioInception PVT LTD
// Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
// See the NOTICE file for attribution, trademark, and algorithm IP terms.
#pragma once
#include "smsd/mol_graph.hpp"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>

namespace smsd {

// ── 2D Point ────────────────────────────────────────────────────────────────
struct Vec2 { double x, y; };
inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x+b.x, a.y+b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x-b.x, a.y-b.y}; }
inline Vec2 operator*(Vec2 a, double s) { return {a.x*s, a.y*s}; }
inline double dot(Vec2 a, Vec2 b) { return a.x*b.x + a.y*b.y; }
inline double len(Vec2 a) { return std::sqrt(dot(a,a)); }
inline Vec2 norm(Vec2 a) { double l=len(a); return l>1e-9? Vec2{a.x/l,a.y/l}: Vec2{0,0}; }
inline Vec2 perp(Vec2 a) { return {-a.y, a.x}; }

// ── Color ───────────────────────────────────────────────────────────────────
struct Color { int r, g, b; double a; };
inline std::string rgba(Color c) {
    std::ostringstream s;
    s << "rgba(" << c.r << "," << c.g << "," << c.b << "," << c.a << ")";
    return s.str();
}
inline std::string rgb(Color c) {
    char buf[8]; snprintf(buf, 8, "#%02X%02X%02X", c.r, c.g, c.b);
    return buf;
}

// ── Publication Color Palette ───────────────────────────────────────────────
namespace palette {
    constexpr Color MATCH_FILL   = {76,175,80, 0.25};    // green highlight
    constexpr Color MATCH_STROKE = {76,175,80, 1.0};     // green bonds
    constexpr Color UNMATCH_FILL = {158,158,158, 0.15};   // gray
    constexpr Color BOND         = {33,33,33, 1.0};       // near-black
    constexpr Color MAPNUM       = {21,101,192, 1.0};     // blue
    constexpr Color BACKGROUND   = {255,255,255, 1.0};    // white
    constexpr Color CARBON       = {33,33,33, 1.0};
    constexpr Color NITROGEN     = {48,63,159, 1.0};      // indigo
    constexpr Color OXYGEN       = {211,47,47, 1.0};      // red
    constexpr Color SULFUR       = {255,179,0, 1.0};      // amber
    constexpr Color PHOSPHORUS   = {255,111,0, 1.0};      // orange
    constexpr Color FLUORINE     = {0,150,136, 1.0};      // teal
    constexpr Color CHLORINE     = {0,150,136, 1.0};
    constexpr Color BROMINE      = {121,85,72, 1.0};      // brown
    constexpr Color IODINE       = {106,27,154, 1.0};     // purple
    constexpr Color OTHER        = {96,125,139, 1.0};     // blue-gray
}

inline Color atomColor(int Z) {
    switch(Z) {
        case 6: return palette::CARBON;
        case 7: return palette::NITROGEN;
        case 8: return palette::OXYGEN;
        case 16: return palette::SULFUR;
        case 15: return palette::PHOSPHORUS;
        case 9: return palette::FLUORINE;
        case 17: return palette::CHLORINE;
        case 35: return palette::BROMINE;
        case 53: return palette::IODINE;
        default: return palette::OTHER;
    }
}

inline const char* elementSymbol(int Z) {
    static const char* SYM[] = {
        "?","H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
        "Sb","Te","I","Xe"
    };
    return (Z >= 0 && Z <= 54) ? SYM[Z] : "?";
}

// ── Layout Options ──────────────────────────────────────────────────────────
struct DepictOptions {
    double bondLength   = 30.0;   // pixels
    double atomRadius   = 12.0;   // highlight circle radius
    double fontSize     = 14.0;   // atom label font size
    double mapNumSize   = 9.0;    // atom map number font size
    double padding      = 40.0;   // SVG padding
    double lineWidth    = 2.0;    // bond line width
    double matchWidth   = 3.0;    // highlighted bond width
    bool   showHydrogen = false;  // show explicit H labels
    bool   showCarbonLabels = false; // show C labels (false = standard)
    bool   showAtomIndices  = false; // show atom index numbers
    bool   showAromaticCircles = true; // inner circle for aromatic rings
    std::string fontFamily = "Helvetica, Arial, sans-serif";
};

// ── 2D Layout Engine ────────────────────────────────────────────────────────

// Find smallest set of smallest rings (SSSR) — simplified for layout
inline std::vector<std::vector<int>> findLayoutRings(const MolGraph& g) {
    std::vector<std::vector<int>> rings;
    int n = g.n;
    if (n == 0) return rings;

    // BFS-based ring detection: for each back-edge, extract ring
    std::vector<int> parent(n, -1);
    std::vector<int> depth(n, -1);
    std::vector<bool> visited(n, false);
    std::set<std::pair<int,int>> ringEdges;

    for (int start = 0; start < n; start++) {
        if (visited[start]) continue;
        std::vector<int> queue = {start};
        visited[start] = true;
        depth[start] = 0;
        for (size_t qi = 0; qi < queue.size(); qi++) {
            int u = queue[qi];
            for (int v : g.neighbors[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    depth[v] = depth[u] + 1;
                    queue.push_back(v);
                } else if (v != parent[u] && depth[v] <= depth[u]) {
                    // Back-edge found — trace ring
                    std::vector<int> ring;
                    int a = u, b = v;
                    std::vector<int> pathA, pathB;
                    while (a != -1) { pathA.push_back(a); a = parent[a]; }
                    while (b != -1) { pathB.push_back(b); b = parent[b]; }
                    // Find common ancestor
                    std::set<int> setA(pathA.begin(), pathA.end());
                    int lca = -1;
                    for (int x : pathB) { if (setA.count(x)) { lca = x; break; } }
                    if (lca == -1) continue;
                    for (int x : pathA) { ring.push_back(x); if (x == lca) break; }
                    std::vector<int> tail;
                    for (int x : pathB) { if (x == lca) break; tail.push_back(x); }
                    std::reverse(tail.begin(), tail.end());
                    for (int x : tail) ring.push_back(x);
                    if (ring.size() >= 3 && ring.size() <= 8) {
                        rings.push_back(ring);
                    }
                }
            }
        }
    }

    // Deduplicate rings by sorted atom set
    std::vector<std::vector<int>> unique;
    std::set<std::vector<int>> seen;
    for (auto& r : rings) {
        auto sorted = r;
        std::sort(sorted.begin(), sorted.end());
        if (!seen.count(sorted)) {
            seen.insert(sorted);
            unique.push_back(r);
        }
    }
    return unique;
}

// Place atoms in 2D
inline std::vector<Vec2> layout2D(const MolGraph& g, const DepictOptions& opts = {}) {
    int n = g.n;
    std::vector<Vec2> pos(n, {0, 0});
    if (n == 0) return pos;
    if (n == 1) { pos[0] = {0, 0}; return pos; }

    double BL = opts.bondLength;
    std::vector<bool> placed(n, false);

    // Step 1: Find and place rings
    auto rings = findLayoutRings(g);

    // Sort rings: largest first (fused rings placed after their parent)
    std::sort(rings.begin(), rings.end(),
        [](const auto& a, const auto& b) { return a.size() > b.size(); });

    for (auto& ring : rings) {
        int sz = (int)ring.size();
        // Check if any atom already placed (fused ring)
        int anchorIdx = -1;
        for (int i = 0; i < sz; i++) {
            if (placed[ring[i]]) { anchorIdx = i; break; }
        }

        if (anchorIdx == -1) {
            // Fresh ring — place as regular polygon centered at origin offset
            double cx = 0, cy = 0;
            // Find a good center position (away from already placed atoms)
            for (int i = 0; i < n; i++) {
                if (placed[i]) { cx = pos[i].x + BL * 2; cy = pos[i].y; break; }
            }
            double R = BL / (2.0 * std::sin(M_PI / sz));
            for (int i = 0; i < sz; i++) {
                double angle = 2.0 * M_PI * i / sz - M_PI / 2.0;
                pos[ring[i]] = {cx + R * std::cos(angle), cy + R * std::sin(angle)};
                placed[ring[i]] = true;
            }
        } else {
            // Fused ring — find shared edge and build from it
            int a1 = ring[anchorIdx];
            // Find another placed neighbor in ring
            int a2 = -1;
            for (int i = 0; i < sz; i++) {
                int j = ring[i];
                if (j != a1 && placed[j]) {
                    // Check if adjacent in ring
                    int prev = ring[(i - 1 + sz) % sz];
                    int next = ring[(i + 1) % sz];
                    if (prev == a1 || next == a1) { a2 = j; break; }
                }
            }
            if (a2 == -1) {
                // No shared edge — place as regular polygon near anchor
                double R = BL / (2.0 * std::sin(M_PI / sz));
                double cx = pos[a1].x + BL;
                double cy = pos[a1].y;
                for (int i = 0; i < sz; i++) {
                    if (!placed[ring[i]]) {
                        double angle = 2.0 * M_PI * i / sz - M_PI / 2.0;
                        pos[ring[i]] = {cx + R * std::cos(angle), cy + R * std::sin(angle)};
                        placed[ring[i]] = true;
                    }
                }
            } else {
                // Shared edge a1-a2: place remaining atoms on opposite side
                Vec2 mid = (pos[a1] + pos[a2]) * 0.5;
                Vec2 dir = norm(pos[a2] - pos[a1]);
                Vec2 outward = perp(dir);
                // Determine which side is free
                double R = BL / (2.0 * std::sin(M_PI / sz));
                Vec2 center1 = mid + outward * R * 0.8;
                Vec2 center2 = mid - outward * R * 0.8;
                // Pick side with fewer placed atoms nearby
                Vec2 center = center1;
                int unplaced = 0;
                for (int i = 0; i < sz; i++) if (!placed[ring[i]]) unplaced++;
                double startAngle = std::atan2(pos[a1].y - center.y, pos[a1].x - center.x);
                for (int k = 0; k < unplaced; k++) {
                    // Find next unplaced in ring order
                    for (int i = 0; i < sz; i++) {
                        if (!placed[ring[i]]) {
                            double angle = startAngle + 2.0 * M_PI * (k + 1) / sz;
                            pos[ring[i]] = {center.x + R * std::cos(angle), center.y + R * std::sin(angle)};
                            placed[ring[i]] = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Step 2: Place chain atoms by DFS from placed atoms
    std::vector<int> stack;
    for (int i = 0; i < n; i++) if (placed[i]) stack.push_back(i);
    // If nothing placed yet (no rings), start from atom 0
    if (stack.empty()) {
        pos[0] = {0, 0};
        placed[0] = true;
        stack.push_back(0);
    }

    while (!stack.empty()) {
        int u = stack.back(); stack.pop_back();
        int neighborCount = 0;
        for (int v : g.neighbors[u]) {
            if (placed[v]) continue;
            neighborCount++;
            // Place at angle away from existing neighbors
            double sumSin = 0, sumCos = 0;
            int usedCount = 0;
            for (int w : g.neighbors[u]) {
                if (placed[w]) {
                    double ang = std::atan2(pos[w].y - pos[u].y, pos[w].x - pos[u].x);
                    sumSin += std::sin(ang);
                    sumCos += std::cos(ang);
                    usedCount++;
                }
            }
            // Vector mean avoids circular-angle averaging trap
            double baseAngle = usedCount > 0 ? std::atan2(sumSin, sumCos) + M_PI : 0;
            // Spread multiple unplaced neighbors
            double spread = M_PI * 2.0 / 3.0; // 120 degrees
            double angle = baseAngle + (neighborCount - 1) * spread / std::max(1, (int)g.neighbors[u].size() - 1);

            pos[v] = {pos[u].x + BL * std::cos(angle), pos[u].y + BL * std::sin(angle)};
            placed[v] = true;
            stack.push_back(v);
        }
    }

    // Step 3: Force-directed refinement (resolve overlaps)
    for (int iter = 0; iter < 50; iter++) {
        bool moved = false;
        for (int i = 0; i < n; i++) {
            Vec2 force = {0, 0};
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                // Skip bonded atoms — repelling them warps ring geometry
                bool bonded = false;
                for (int nb : g.neighbors[i]) if (nb == j) { bonded = true; break; }
                if (bonded) continue;
                Vec2 d = pos[i] - pos[j];
                double dist = len(d);
                if (dist < BL * 0.8 && dist > 1e-6) {
                    // Repulsion
                    double f = (BL * 0.8 - dist) * 0.3;
                    force = force + norm(d) * f;
                    moved = true;
                }
            }
            // Don't move ring atoms much
            double dampening = g.ring[i] ? 0.1 : 0.5;
            pos[i] = pos[i] + force * dampening;
        }
        if (!moved) break;
    }

    // Step 4: Post-pass collision detection — jitter overlapping atoms apart
    const double overlapThresh = BL * 0.3;
    for (int iter = 0; iter < 20; iter++) {
        bool collisionFound = false;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                Vec2 d = pos[i] - pos[j];
                double dist = len(d);
                if (dist < overlapThresh) {
                    collisionFound = true;
                    // Apply jitter: push apart along the connecting vector,
                    // or along a deterministic diagonal if coincident
                    Vec2 dir;
                    if (dist > 1e-9) {
                        dir = norm(d);
                    } else {
                        // Deterministic jitter direction based on atom indices
                        double angle = (i * 2.654 + j * 1.337);
                        dir = {std::cos(angle), std::sin(angle)};
                    }
                    double push = (overlapThresh - dist) * 0.5 + 0.01;
                    pos[i] = pos[i] + dir * push;
                    pos[j] = pos[j] - dir * push;
                }
            }
        }
        if (!collisionFound) break;
    }

    return pos;
}

// ── SVG Renderer ────────────────────────────────────────────────────────────

struct DepictResult {
    std::string svg;
    int width, height;
};

inline DepictResult depictSVG(
    const MolGraph& g,
    const std::vector<Vec2>& pos,
    const std::map<int,int>& mapping = {},  // atom index → map number
    const std::set<int>& highlightAtoms = {},
    const std::set<std::pair<int,int>>& highlightBonds = {},
    const DepictOptions& opts = {}
) {
    // Compute bounding box
    double minX = 1e9, minY = 1e9, maxX = -1e9, maxY = -1e9;
    for (int i = 0; i < g.n; i++) {
        minX = std::min(minX, pos[i].x);
        minY = std::min(minY, pos[i].y);
        maxX = std::max(maxX, pos[i].x);
        maxY = std::max(maxY, pos[i].y);
    }
    double pad = opts.padding;
    int W = (int)(maxX - minX + 2 * pad);
    int H = (int)(maxY - minY + 2 * pad);
    if (W < 100) W = 100;
    if (H < 100) H = 100;

    auto tx = [&](double x) { return x - minX + pad; };
    auto ty = [&](double y) { return y - minY + pad; };

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\"";
    svg << " viewBox=\"0 0 " << W << " " << H << "\"";
    svg << " style=\"background:" << rgb(palette::BACKGROUND) << "\">\n";
    svg << "<defs>\n";
    svg << "  <style>text{font-family:" << opts.fontFamily << ";dominant-baseline:central;text-anchor:middle;}</style>\n";
    svg << "</defs>\n";

    // ── Highlight fills (behind everything) ─────────────────────────────────
    for (int i = 0; i < g.n; i++) {
        if (highlightAtoms.count(i)) {
            svg << "<circle cx=\"" << tx(pos[i].x) << "\" cy=\"" << ty(pos[i].y)
                << "\" r=\"" << opts.atomRadius << "\" fill=\"" << rgba(palette::MATCH_FILL) << "\" stroke=\"none\"/>\n";
        }
    }

    // ── Bonds ───────────────────────────────────────────────────────────────
    std::set<std::pair<int,int>> drawnBonds;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue; // draw each bond once
            auto key = std::make_pair(i, j);
            if (drawnBonds.count(key)) continue;
            drawnBonds.insert(key);

            double x1 = tx(pos[i].x), y1 = ty(pos[i].y);
            double x2 = tx(pos[j].x), y2 = ty(pos[j].y);

            bool isHighlighted = highlightBonds.count({i,j}) || highlightBonds.count({j,i}) ||
                                 (highlightAtoms.count(i) && highlightAtoms.count(j));
            Color bondColor = isHighlighted ? palette::MATCH_STROKE : palette::BOND;
            double width = isHighlighted ? opts.matchWidth : opts.lineWidth;

            int order = g.bondOrder(i, j);
            bool arom = g.bondAromatic(i, j);

            if (order == 1 || order == 0 || arom) {
                // Single bond (or aromatic — drawn as single + circle)
                svg << "<line x1=\"" << x1 << "\" y1=\"" << y1
                    << "\" x2=\"" << x2 << "\" y2=\"" << y2
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
            } else if (order == 2) {
                // Double bond — two parallel lines
                Vec2 d = norm(Vec2{x2-x1, y2-y1});
                Vec2 p = perp(d) * 2.5;
                svg << "<line x1=\"" << x1+p.x << "\" y1=\"" << y1+p.y
                    << "\" x2=\"" << x2+p.x << "\" y2=\"" << y2+p.y
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
                svg << "<line x1=\"" << x1-p.x << "\" y1=\"" << y1-p.y
                    << "\" x2=\"" << x2-p.x << "\" y2=\"" << y2-p.y
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
            } else if (order == 3) {
                // Triple bond — three parallel lines
                Vec2 d = norm(Vec2{x2-x1, y2-y1});
                Vec2 p = perp(d) * 3.5;
                svg << "<line x1=\"" << x1 << "\" y1=\"" << y1
                    << "\" x2=\"" << x2 << "\" y2=\"" << y2
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
                svg << "<line x1=\"" << x1+p.x << "\" y1=\"" << y1+p.y
                    << "\" x2=\"" << x2+p.x << "\" y2=\"" << y2+p.y
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width * 0.7
                    << "\" stroke-linecap=\"round\"/>\n";
                svg << "<line x1=\"" << x1-p.x << "\" y1=\"" << y1-p.y
                    << "\" x2=\"" << x2-p.x << "\" y2=\"" << y2-p.y
                    << "\" stroke=\"" << rgb(bondColor) << "\" stroke-width=\"" << width * 0.7
                    << "\" stroke-linecap=\"round\"/>\n";
            }
        }
    }

    // ── Aromatic ring circles ───────────────────────────────────────────────
    if (opts.showAromaticCircles) {
        auto rings = findLayoutRings(g);
        for (auto& ring : rings) {
            bool allAromatic = true;
            for (int a : ring) if (!g.aromatic[a]) { allAromatic = false; break; }
            if (!allAromatic) continue;
            // Compute ring center
            double cx = 0, cy = 0;
            for (int a : ring) { cx += tx(pos[a].x); cy += ty(pos[a].y); }
            cx /= ring.size(); cy /= ring.size();
            double R = 0;
            for (int a : ring) {
                double dx = tx(pos[a].x) - cx, dy = ty(pos[a].y) - cy;
                R += std::sqrt(dx*dx + dy*dy);
            }
            R = R / ring.size() * 0.55;
            svg << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"" << R
                << "\" fill=\"none\" stroke=\"" << rgb(palette::BOND)
                << "\" stroke-width=\"1\" stroke-dasharray=\"3,2\"/>\n";
        }
    }

    // ── Atom labels ─────────────────────────────────────────────────────────
    for (int i = 0; i < g.n; i++) {
        double x = tx(pos[i].x), y = ty(pos[i].y);
        int Z = g.atomicNum[i];
        bool isCarbon = (Z == 6);

        // White circle behind label to mask bonds
        if (!isCarbon || opts.showCarbonLabels) {
            svg << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << opts.fontSize * 0.5
                << "\" fill=\"" << rgb(palette::BACKGROUND) << "\" stroke=\"none\"/>\n";
        }

        // Atom symbol
        if (!isCarbon || opts.showCarbonLabels) {
            Color col = atomColor(Z);
            svg << "<text x=\"" << x << "\" y=\"" << y
                << "\" font-size=\"" << opts.fontSize << "\" fill=\"" << rgb(col)
                << "\" font-weight=\"bold\">" << elementSymbol(Z) << "</text>\n";
        }

        // Atom map number (superscript)
        if (mapping.count(i)) {
            svg << "<text x=\"" << x + opts.fontSize * 0.5 << "\" y=\"" << y - opts.fontSize * 0.5
                << "\" font-size=\"" << opts.mapNumSize << "\" fill=\"" << rgb(palette::MAPNUM)
                << "\" font-weight=\"bold\">" << mapping.at(i) << "</text>\n";
        }

        // Atom index (subscript, optional)
        if (opts.showAtomIndices) {
            svg << "<text x=\"" << x + opts.fontSize * 0.5 << "\" y=\"" << y + opts.fontSize * 0.6
                << "\" font-size=\"" << opts.mapNumSize << "\" fill=\"" << rgb(palette::UNMATCH_FILL)
                << "\">" << i << "</text>\n";
        }
    }

    svg << "</svg>\n";
    return {svg.str(), W, H};
}

// ── High-Level API ──────────────────────────────────────────────────────────

// Depict a single molecule
inline std::string depict(const MolGraph& g, const DepictOptions& opts = {}) {
    auto pos = layout2D(g, opts);
    return depictSVG(g, pos, {}, {}, {}, opts).svg;
}

// Depict molecule with MCS highlighting
inline std::string depictWithMapping(
    const MolGraph& g,
    const std::map<int,int>& atomMapping,  // this mol's atom → partner atom
    const DepictOptions& opts = {}
) {
    auto pos = layout2D(g, opts);
    std::set<int> highlighted;
    std::map<int,int> mapNums;
    int mapNum = 1;
    for (auto& [qi, tj] : atomMapping) {
        highlighted.insert(qi);
        mapNums[qi] = mapNum++;
    }
    return depictSVG(g, pos, mapNums, highlighted, {}, opts).svg;
}

// Depict a pair of molecules side-by-side with MCS highlighting
inline std::string depictPair(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& mapping,  // g1 atom → g2 atom
    const DepictOptions& opts = {}
) {
    auto pos1 = layout2D(g1, opts);
    auto pos2 = layout2D(g2, opts);

    // Offset g2 to the right of g1
    double maxX1 = -1e9;
    for (int i = 0; i < g1.n; i++) maxX1 = std::max(maxX1, pos1[i].x);
    double minX2 = 1e9;
    for (int i = 0; i < g2.n; i++) minX2 = std::min(minX2, pos2[i].x);
    double offset = maxX1 - minX2 + opts.bondLength * 3;
    for (int i = 0; i < g2.n; i++) pos2[i].x += offset;

    // Build combined molecule (virtual) for rendering
    int n = g1.n + g2.n;
    std::vector<Vec2> pos(n);
    for (int i = 0; i < g1.n; i++) pos[i] = pos1[i];
    for (int i = 0; i < g2.n; i++) pos[g1.n + i] = pos2[i];

    // Highlights
    std::set<int> highlights;
    std::map<int,int> mapNums;
    int mapNum = 1;
    for (auto& [qi, tj] : mapping) {
        highlights.insert(qi);
        highlights.insert(g1.n + tj);
        mapNums[qi] = mapNum;
        mapNums[g1.n + tj] = mapNum;
        mapNum++;
    }

    // Bounding box
    double minX = 1e9, minY = 1e9, maxX = -1e9, maxY = -1e9;
    for (int i = 0; i < n; i++) {
        minX = std::min(minX, pos[i].x); minY = std::min(minY, pos[i].y);
        maxX = std::max(maxX, pos[i].x); maxY = std::max(maxY, pos[i].y);
    }
    double pad = opts.padding;
    int W = (int)(maxX - minX + 2 * pad);
    int H = (int)(maxY - minY + 2 * pad);

    auto tx = [&](double x) { return x - minX + pad; };
    auto ty = [&](double y) { return y - minY + pad; };

    std::ostringstream svg;
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\">\n";
    svg << "<defs><style>text{font-family:" << opts.fontFamily << ";dominant-baseline:central;text-anchor:middle;}</style></defs>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";

    // Draw both molecules
    auto drawMol = [&](const MolGraph& g, int offset) {
        // Highlight fills
        for (int i = 0; i < g.n; i++) {
            int gi = offset + i;
            if (highlights.count(gi)) {
                svg << "<circle cx=\"" << tx(pos[gi].x) << "\" cy=\"" << ty(pos[gi].y)
                    << "\" r=\"" << opts.atomRadius << "\" fill=\"" << rgba(palette::MATCH_FILL) << "\"/>\n";
            }
        }
        // Bonds
        for (int i = 0; i < g.n; i++) {
            for (int j : g.neighbors[i]) {
                if (j <= i) continue;
                int gi = offset + i, gj = offset + j;
                bool hl = highlights.count(gi) && highlights.count(gj);
                Color c = hl ? palette::MATCH_STROKE : palette::BOND;
                double w = hl ? opts.matchWidth : opts.lineWidth;
                int order = g.bondOrder(i, j);

                double x1 = tx(pos[gi].x), y1 = ty(pos[gi].y);
                double x2 = tx(pos[gj].x), y2 = ty(pos[gj].y);

                if (order <= 1 || order == 4) {
                    svg << "<line x1=\"" << x1 << "\" y1=\"" << y1
                        << "\" x2=\"" << x2 << "\" y2=\"" << y2
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w << "\"/>\n";
                } else if (order == 2) {
                    Vec2 d = norm({x2-x1, y2-y1});
                    Vec2 p = perp(d) * 2.5;
                    svg << "<line x1=\"" << x1+p.x << "\" y1=\"" << y1+p.y
                        << "\" x2=\"" << x2+p.x << "\" y2=\"" << y2+p.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w << "\"/>\n";
                    svg << "<line x1=\"" << x1-p.x << "\" y1=\"" << y1-p.y
                        << "\" x2=\"" << x2-p.x << "\" y2=\"" << y2-p.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w << "\"/>\n";
                } else if (order == 3) {
                    Vec2 d = norm({x2-x1, y2-y1});
                    Vec2 p = perp(d) * 3.5;
                    svg << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w << "\"/>\n";
                    svg << "<line x1=\"" << x1+p.x << "\" y1=\"" << y1+p.y << "\" x2=\"" << x2+p.x << "\" y2=\"" << y2+p.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w*0.7 << "\"/>\n";
                    svg << "<line x1=\"" << x1-p.x << "\" y1=\"" << y1-p.y << "\" x2=\"" << x2-p.x << "\" y2=\"" << y2-p.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w*0.7 << "\"/>\n";
                }
            }
        }
        // Atom labels
        for (int i = 0; i < g.n; i++) {
            int gi = offset + i;
            double x = tx(pos[gi].x), y = ty(pos[gi].y);
            int Z = g.atomicNum[i];
            if (Z != 6 || opts.showCarbonLabels) {
                svg << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << opts.fontSize*0.5
                    << "\" fill=\"white\" stroke=\"none\"/>\n";
                svg << "<text x=\"" << x << "\" y=\"" << y << "\" font-size=\"" << opts.fontSize
                    << "\" fill=\"" << rgb(atomColor(Z)) << "\" font-weight=\"bold\">"
                    << elementSymbol(Z) << "</text>\n";
            }
            if (mapNums.count(gi)) {
                svg << "<text x=\"" << x+opts.fontSize*0.5 << "\" y=\"" << y-opts.fontSize*0.5
                    << "\" font-size=\"" << opts.mapNumSize << "\" fill=\"" << rgb(palette::MAPNUM)
                    << "\" font-weight=\"bold\">" << mapNums[gi] << "</text>\n";
            }
        }
    };

    drawMol(g1, 0);

    // Separator arrow (guard against empty molecules)
    double sepX;
    if (g1.n > 0 && g2.n > 0) {
        sepX = tx((pos[g1.n - 1].x + pos[g1.n].x) / 2);
    } else if (g1.n > 0) {
        sepX = tx(pos[g1.n - 1].x) + opts.bondLength * 1.5;
    } else if (g2.n > 0) {
        sepX = tx(pos[g1.n].x) - opts.bondLength * 1.5;
    } else {
        sepX = opts.width / 2.0;
    }
    double sepY = ty((minY + maxY) / 2);
    svg << "<text x=\"" << sepX << "\" y=\"" << sepY
        << "\" font-size=\"24\" fill=\"#666\" font-weight=\"bold\">&#x2194;</text>\n";

    drawMol(g2, g1.n);

    svg << "</svg>\n";
    return svg.str();
}

// ── Convenience: SMILES → SVG ───────────────────────────────────────────────

// Requires smiles_parser.hpp to be included before this
#ifdef SMSD_SMILES_PARSER_HPP

inline std::string depictSMILES(const std::string& smiles, const DepictOptions& opts = {}) {
    auto mol = parseSMILES(smiles);
    return depict(mol, opts);
}

inline std::string depictMCS(
    const std::string& smiles1, const std::string& smiles2,
    const std::map<int,int>& mapping,
    const DepictOptions& opts = {}
) {
    auto mol1 = parseSMILES(smiles1);
    auto mol2 = parseSMILES(smiles2);
    return depictPair(mol1, mol2, mapping, opts);
}

#endif // SMSD_SMILES_PARSER_HPP

} // namespace smsd

// SMSD Pro — Publication-Quality Molecular Depiction Engine
// Zero-dependency SVG renderer: ACS 1996, Nature, Springer journal standard
// Accepts MolGraph, SMILES, substructure/MCS mappings; auto-layout when needed
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
#include <iomanip>

namespace smsd {

// ════════════════════════════════════════════════════════════════════════════
//  2D Geometry Primitives
// ════════════════════════════════════════════════════════════════════════════
struct Vec2 { double x, y; };
inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x+b.x, a.y+b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x-b.x, a.y-b.y}; }
inline Vec2 operator*(Vec2 a, double s) { return {a.x*s, a.y*s}; }
inline double dot(Vec2 a, Vec2 b) { return a.x*b.x + a.y*b.y; }
inline double len(Vec2 a) { return std::sqrt(dot(a,a)); }
inline Vec2 norm(Vec2 a) { double l=len(a); return l>1e-9? Vec2{a.x/l,a.y/l}: Vec2{0,0}; }
inline Vec2 perp(Vec2 a) { return {-a.y, a.x}; }

// ════════════════════════════════════════════════════════════════════════════
//  Color
// ════════════════════════════════════════════════════════════════════════════
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

// ════════════════════════════════════════════════════════════════════════════
//  Jmol/CPK Element Color Palette (publication standard)
// ════════════════════════════════════════════════════════════════════════════
namespace jmol {
    constexpr Color H   = {255,255,255, 1.0};
    constexpr Color He  = {217,255,255, 1.0};
    constexpr Color Li  = {204,128,255, 1.0};
    constexpr Color Be  = {194,255,  0, 1.0};
    constexpr Color B   = {255,181,181, 1.0};
    constexpr Color C   = {144,144,144, 1.0};
    constexpr Color N   = { 48, 80,248, 1.0};
    constexpr Color O   = {255, 13, 13, 1.0};
    constexpr Color F   = {144,224, 80, 1.0};
    constexpr Color Ne  = {179,227,245, 1.0};
    constexpr Color Na  = {171, 92,242, 1.0};
    constexpr Color Mg  = {138,255,  0, 1.0};
    constexpr Color Al  = {191,166,166, 1.0};
    constexpr Color Si  = {240,200,160, 1.0};
    constexpr Color P   = {255,128,  0, 1.0};
    constexpr Color S   = {178,178,  0, 1.0};  // darkened for white bg
    constexpr Color Cl  = { 31,240, 31, 1.0};
    constexpr Color Ar  = {128,209,227, 1.0};
    constexpr Color K   = {143, 64,212, 1.0};
    constexpr Color Ca  = { 61,255,  0, 1.0};
    constexpr Color Fe  = {224,102, 51, 1.0};
    constexpr Color Cu  = {200,128, 51, 1.0};
    constexpr Color Zn  = {125,128,176, 1.0};
    constexpr Color Br  = {166, 41, 41, 1.0};
    constexpr Color Se  = {255,161,  0, 1.0};
    constexpr Color I   = {148,  0,148, 1.0};
    constexpr Color UNKNOWN = {255, 20,147, 1.0};
}

// ── Default publication palette ────────────────────────────────────────────
namespace palette {
    // Highlight colors
    constexpr Color MATCH_FILL   = { 76,175, 80, 0.25};  // translucent green
    constexpr Color MATCH_STROKE = { 46,125, 50, 1.0};   // deep green
    constexpr Color UNMATCH_FILL = {158,158,158, 0.10};   // very light gray
    // Structural colors
    constexpr Color BOND         = { 33, 33, 33, 1.0};   // near-black
    constexpr Color MAPNUM       = { 21,101,192, 1.0};   // blue
    constexpr Color BACKGROUND   = {255,255,255, 1.0};   // white
    // Highlight color set for multi-fragment / multi-match
    constexpr Color HIGHLIGHT_1  = { 76,175, 80, 0.30};  // green
    constexpr Color HIGHLIGHT_2  = { 33,150,243, 0.30};  // blue
    constexpr Color HIGHLIGHT_3  = {255,152,  0, 0.30};  // orange
    constexpr Color HIGHLIGHT_4  = {156, 39,176, 0.30};  // purple
    constexpr Color HIGHLIGHT_5  = {  0,188,212, 0.30};  // cyan
    constexpr Color HIGHLIGHT_6  = {233, 30, 99, 0.30};  // pink
}

// Element color by atomic number (Jmol standard, C darkened for skeletal)
inline Color atomColor(int Z) {
    switch (Z) {
        case  1: return jmol::H;
        case  5: return jmol::B;
        case  6: return {33, 33, 33, 1.0};  // carbon: near-black for skeletal
        case  7: return jmol::N;
        case  8: return jmol::O;
        case  9: return jmol::F;
        case 15: return jmol::P;
        case 16: return jmol::S;
        case 17: return jmol::Cl;
        case 35: return jmol::Br;
        case 34: return jmol::Se;
        case 53: return jmol::I;
        case 11: return jmol::Na;
        case 12: return jmol::Mg;
        case 20: return jmol::Ca;
        case 26: return jmol::Fe;
        case 29: return jmol::Cu;
        case 30: return jmol::Zn;
        default: return {96,125,139, 1.0};  // blue-gray
    }
}

inline const char* elementSymbol(int Z) {
    static const char* SYM[] = {
        "?","H","He","Li","Be","B","C","N","O","F","Ne",
        "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
        "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
        "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
        "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
        "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
        "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
        "Tl","Pb","Bi","Po","At","Rn"
    };
    return (Z >= 0 && Z <= 86) ? SYM[Z] : "?";
}

// ════════════════════════════════════════════════════════════════════════════
//  ACS 1996 Standard Constants
// ════════════════════════════════════════════════════════════════════════════
namespace acs {
    // All dimensions are in points (1pt = 1/72 inch); we scale to pixels.
    // At 96 DPI: 1pt = 1.333px.  We use bond length as the reference unit.
    constexpr double BOND_LENGTH_PT    = 14.4;     // standard bond length
    constexpr double LINE_WIDTH_PT     =  0.6;     // bond line width
    constexpr double BOLD_WIDTH_PT     =  2.0;     // wedge bond max width
    constexpr double MARGIN_WIDTH_PT   =  1.6;     // gap between double bond lines
    constexpr double HASH_SPACING_PT   =  2.5;     // dashed wedge stripe spacing
    constexpr double CHAIN_ANGLE       = 120.0;    // degrees
    constexpr double BOND_SPACING      =  0.18;    // fraction of bond length for double bond offset
    constexpr double FONT_SIZE_PT      = 10.0;     // atom label font size
    constexpr double SUBSCRIPT_SCALE   =  0.70;    // subscript/superscript relative size
}

// ════════════════════════════════════════════════════════════════════════════
//  Depiction Options (everything customizable)
// ════════════════════════════════════════════════════════════════════════════

enum class BondStereoType { NONE = 0, WEDGE_UP = 1, WEDGE_DOWN = 2, WAVY = 3 };
enum class AromaticStyle  { CIRCLE, KEKULE };
enum class HydrogenMode   { HETERO_ONLY, ALL, NONE, TERMINAL };

struct DepictOptions {
    // ── Geometry (scaled from ACS 1996 standard) ──
    double bondLength      = 30.0;    // pixels — the reference unit
    double lineWidth       = 0.0;     // 0 = auto from ACS ratio
    double boldWidth       = 0.0;     // 0 = auto from ACS ratio (wedge max)
    double bondSpacing     = 0.0;     // 0 = auto (18% of bondLength)
    double hashSpacing     = 0.0;     // 0 = auto from ACS ratio

    // ── Atom labels ──
    double fontSize        = 0.0;     // 0 = auto from ACS ratio
    double subscriptScale  = 0.70;    // H-count and charge subscript/superscript
    HydrogenMode hMode     = HydrogenMode::HETERO_ONLY;
    bool   showCarbonLabels = false;
    bool   showAtomIndices  = false;

    // ── Highlights ──
    double highlightRadius = 0.0;     // 0 = auto (40% of bondLength)
    double highlightOpacity = 0.30;
    double matchWidth      = 0.0;     // 0 = auto (2x lineWidth)

    // ── Aromatic rendering ──
    AromaticStyle aromaticStyle = AromaticStyle::CIRCLE;

    // ── Map numbers ──
    double mapNumSize      = 0.0;     // 0 = auto
    bool   showMapNumbers  = true;

    // ── SVG / canvas ──
    double padding         = 30.0;
    int    width           = 0;       // 0 = auto-fit
    int    height          = 0;       // 0 = auto-fit
    std::string fontFamily = "Arial, Helvetica, sans-serif";
    Color  backgroundColor = palette::BACKGROUND;

    // ── Custom colors (override defaults) ──
    Color  bondColor       = palette::BOND;
    Color  matchFill       = palette::MATCH_FILL;
    Color  matchStroke     = palette::MATCH_STROKE;
    std::map<int, Color> atomColorOverrides;   // Z -> Color
    std::map<int, Color> atomHighlightColors;  // atom index -> Color
    std::map<std::pair<int,int>, Color> bondHighlightColors; // bond -> Color

    // ── Computed values (resolved from ACS ratios) ──
    double lineWidth_() const {
        return lineWidth > 0 ? lineWidth
            : bondLength * (acs::LINE_WIDTH_PT / acs::BOND_LENGTH_PT);
    }
    double boldWidth_() const {
        return boldWidth > 0 ? boldWidth
            : bondLength * (acs::BOLD_WIDTH_PT / acs::BOND_LENGTH_PT);
    }
    double bondSpacing_() const {
        return bondSpacing > 0 ? bondSpacing : bondLength * acs::BOND_SPACING;
    }
    double hashSpacing_() const {
        return hashSpacing > 0 ? hashSpacing
            : bondLength * (acs::HASH_SPACING_PT / acs::BOND_LENGTH_PT);
    }
    double fontSize_() const {
        return fontSize > 0 ? fontSize
            : bondLength * (acs::FONT_SIZE_PT / acs::BOND_LENGTH_PT);
    }
    double mapNumSize_() const {
        return mapNumSize > 0 ? mapNumSize : fontSize_() * 0.65;
    }
    double highlightRadius_() const {
        return highlightRadius > 0 ? highlightRadius : bondLength * 0.40;
    }
    double matchWidth_() const {
        return matchWidth > 0 ? matchWidth : lineWidth_() * 2.2;
    }
};

// ════════════════════════════════════════════════════════════════════════════
//  SSSR Ring Detection (for layout and double-bond offset)
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<std::vector<int>> findLayoutRings(const MolGraph& g) {
    std::vector<std::vector<int>> rings;
    int n = g.n;
    if (n == 0) return rings;

    std::vector<int> parent(n, -1);
    std::vector<int> depth(n, -1);
    std::vector<bool> visited(n, false);

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
                    std::vector<int> ring;
                    int a = u, b = v;
                    std::vector<int> pathA, pathB;
                    while (a != -1) { pathA.push_back(a); a = parent[a]; }
                    while (b != -1) { pathB.push_back(b); b = parent[b]; }
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

    // Deduplicate by sorted atom set
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

// ════════════════════════════════════════════════════════════════════════════
//  2D Layout Engine (ring-first + chain zig-zag + force refinement)
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<Vec2> layout2D(const MolGraph& g, const DepictOptions& opts = {}) {
    int n = g.n;
    std::vector<Vec2> pos(n, {0, 0});
    if (n == 0) return pos;
    if (n == 1) return pos;

    double BL = opts.bondLength;
    std::vector<bool> placed(n, false);

    // Step 1: Place rings
    auto rings = findLayoutRings(g);
    std::sort(rings.begin(), rings.end(),
        [](const auto& a, const auto& b) { return a.size() > b.size(); });

    for (auto& ring : rings) {
        int sz = (int)ring.size();
        int anchorIdx = -1;
        for (int i = 0; i < sz; i++)
            if (placed[ring[i]]) { anchorIdx = i; break; }

        if (anchorIdx == -1) {
            double cx = 0, cy = 0;
            for (int i = 0; i < n; i++)
                if (placed[i]) { cx = pos[i].x + BL * 2.5; cy = pos[i].y; break; }
            double R = BL / (2.0 * std::sin(M_PI / sz));
            for (int i = 0; i < sz; i++) {
                double angle = 2.0 * M_PI * i / sz - M_PI / 2.0;
                pos[ring[i]] = {cx + R * std::cos(angle), cy + R * std::sin(angle)};
                placed[ring[i]] = true;
            }
        } else {
            int a1 = ring[anchorIdx];
            int a2 = -1;
            for (int i = 0; i < sz; i++) {
                int j = ring[i];
                if (j != a1 && placed[j]) {
                    int prev = ring[(i - 1 + sz) % sz];
                    int next = ring[(i + 1) % sz];
                    if (prev == a1 || next == a1) { a2 = j; break; }
                }
            }
            if (a2 == -1) {
                double R = BL / (2.0 * std::sin(M_PI / sz));
                double cx = pos[a1].x + BL;
                for (int i = 0; i < sz; i++) {
                    if (!placed[ring[i]]) {
                        double angle = 2.0 * M_PI * i / sz - M_PI / 2.0;
                        pos[ring[i]] = {cx + R * std::cos(angle), pos[a1].y + R * std::sin(angle)};
                        placed[ring[i]] = true;
                    }
                }
            } else {
                Vec2 mid = (pos[a1] + pos[a2]) * 0.5;
                Vec2 dir = norm(pos[a2] - pos[a1]);
                Vec2 outward = perp(dir);
                double R = BL / (2.0 * std::sin(M_PI / sz));
                Vec2 center1 = mid + outward * R * 0.8;
                Vec2 center2 = mid - outward * R * 0.8;
                // Pick side farther from average of already-placed ring atoms
                double score1 = 0, score2 = 0;
                int cnt = 0;
                for (int i = 0; i < sz; i++) {
                    if (placed[ring[i]] && ring[i] != a1 && ring[i] != a2) {
                        score1 += len(pos[ring[i]] - center1);
                        score2 += len(pos[ring[i]] - center2);
                        cnt++;
                    }
                }
                Vec2 center = (cnt == 0 || score1 >= score2) ? center1 : center2;
                double startAngle = std::atan2(pos[a1].y - center.y, pos[a1].x - center.x);
                int unplaced = 0;
                for (int i = 0; i < sz; i++) if (!placed[ring[i]]) unplaced++;
                int k = 0;
                for (int i = 0; i < sz; i++) {
                    if (!placed[ring[i]]) {
                        double angle = startAngle + 2.0 * M_PI * (++k) / sz;
                        pos[ring[i]] = {center.x + R * std::cos(angle), center.y + R * std::sin(angle)};
                        placed[ring[i]] = true;
                    }
                }
            }
        }
    }

    // Step 2: Place chain atoms (DFS with 120-degree zig-zag)
    std::vector<int> stack;
    for (int i = 0; i < n; i++) if (placed[i]) stack.push_back(i);
    if (stack.empty()) {
        pos[0] = {0, 0};
        placed[0] = true;
        stack.push_back(0);
    }

    while (!stack.empty()) {
        int u = stack.back(); stack.pop_back();
        int branchIdx = 0;
        for (int v : g.neighbors[u]) {
            if (placed[v]) continue;
            branchIdx++;
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
            double baseAngle = usedCount > 0 ? std::atan2(sumSin, sumCos) + M_PI : 0;
            double spread = 2.0 * M_PI / 3.0;
            double angle = baseAngle + (branchIdx - 1) * spread
                / std::max(1, (int)g.neighbors[u].size() - 1);
            pos[v] = {pos[u].x + BL * std::cos(angle), pos[u].y + BL * std::sin(angle)};
            placed[v] = true;
            stack.push_back(v);
        }
    }

    // Step 3: Force-directed refinement (non-bonded repulsion only)
    for (int iter = 0; iter < 50; iter++) {
        bool moved = false;
        for (int i = 0; i < n; i++) {
            Vec2 force = {0, 0};
            for (int j = 0; j < n; j++) {
                if (i == j) continue;
                bool bonded = false;
                for (int nb : g.neighbors[i]) if (nb == j) { bonded = true; break; }
                if (bonded) continue;
                Vec2 d = pos[i] - pos[j];
                double dist = len(d);
                if (dist < BL * 0.8 && dist > 1e-6) {
                    force = force + norm(d) * ((BL * 0.8 - dist) * 0.3);
                    moved = true;
                }
            }
            double damp = g.ring[i] ? 0.1 : 0.5;
            pos[i] = pos[i] + force * damp;
        }
        if (!moved) break;
    }

    // Step 4: Jitter coincident atoms
    const double overlapThresh = BL * 0.3;
    for (int iter = 0; iter < 20; iter++) {
        bool collision = false;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                Vec2 d = pos[i] - pos[j];
                double dist = len(d);
                if (dist < overlapThresh) {
                    collision = true;
                    Vec2 dir;
                    if (dist > 1e-9) dir = norm(d);
                    else dir = {std::cos(i * 2.654 + j * 1.337), std::sin(i * 2.654 + j * 1.337)};
                    double push = (overlapThresh - dist) * 0.5 + 0.01;
                    pos[i] = pos[i] + dir * push;
                    pos[j] = pos[j] - dir * push;
                }
            }
        }
        if (!collision) break;
    }

    return pos;
}

// ════════════════════════════════════════════════════════════════════════════
//  Label Builder (element + implicit H + charge, with text width estimate)
// ════════════════════════════════════════════════════════════════════════════
namespace detail_depict {

struct AtomLabel {
    std::string symbol;         // "N", "O", "S", etc.
    int         hCount = 0;     // implicit H to display
    int         charge = 0;     // formal charge
    bool        visible = false;
    double      textWidth = 0;  // estimated width in pixels
    bool        hLeft = false;  // place H on left side (e.g. HO-, H2N-)
};

inline AtomLabel buildLabel(const MolGraph& g, int idx, const DepictOptions& opts,
                            const std::vector<Vec2>& pos) {
    AtomLabel lbl;
    int Z = g.atomicNum[idx];
    lbl.symbol = elementSymbol(Z);

    // Formal charge
    lbl.charge = (!g.formalCharge.empty() && idx < (int)g.formalCharge.size())
        ? g.formalCharge[idx] : 0;

    // Implicit H count
    lbl.hCount = (!g.hydrogenCount.empty() && idx < (int)g.hydrogenCount.size())
        ? g.hydrogenCount[idx] : 0;

    // Visibility
    bool isCarbon = (Z == 6);
    switch (opts.hMode) {
        case HydrogenMode::ALL:
            lbl.visible = true;
            break;
        case HydrogenMode::NONE:
            lbl.visible = !isCarbon || opts.showCarbonLabels || lbl.charge != 0;
            lbl.hCount = 0;
            break;
        case HydrogenMode::TERMINAL:
            lbl.visible = !isCarbon || opts.showCarbonLabels
                || (int)g.neighbors[idx].size() <= 1 || lbl.charge != 0;
            break;
        case HydrogenMode::HETERO_ONLY:
        default:
            lbl.visible = !isCarbon || opts.showCarbonLabels || lbl.charge != 0;
            break;
    }

    if (!lbl.visible) return lbl;

    // Determine H placement: left if majority of bonds come from the right
    if (lbl.hCount > 0) {
        double rightCount = 0;
        for (int nb : g.neighbors[idx]) {
            if (pos[nb].x > pos[idx].x + 0.01) rightCount++;
        }
        lbl.hLeft = (rightCount > (int)g.neighbors[idx].size() / 2.0);
    }

    // Estimate text width (approximate: each char ~ 0.6 * fontSize)
    double charW = opts.fontSize_() * 0.6;
    double subW  = opts.fontSize_() * opts.subscriptScale * 0.5;
    lbl.textWidth = lbl.symbol.length() * charW;
    if (lbl.hCount > 0) {
        lbl.textWidth += charW;  // "H"
        if (lbl.hCount > 1) lbl.textWidth += subW;  // subscript digit
    }
    if (lbl.charge != 0) {
        lbl.textWidth += subW;  // "+" or "-"
        if (std::abs(lbl.charge) > 1) lbl.textWidth += subW;  // digit
    }

    return lbl;
}

// Compute ring centroid (for double bond offset direction)
inline Vec2 ringCentroid(const std::vector<int>& ring, const std::vector<Vec2>& pos) {
    Vec2 c = {0, 0};
    for (int a : ring) { c.x += pos[a].x; c.y += pos[a].y; }
    c.x /= ring.size();
    c.y /= ring.size();
    return c;
}

// Check if bond (i,j) is in a ring and return ring centroid
inline bool bondInRing(int i, int j, const std::vector<std::vector<int>>& rings,
                       const std::vector<Vec2>& pos, Vec2& centroid) {
    for (auto& ring : rings) {
        bool hasI = false, hasJ = false;
        for (int a : ring) {
            if (a == i) hasI = true;
            if (a == j) hasJ = true;
        }
        if (hasI && hasJ) {
            centroid = ringCentroid(ring, pos);
            return true;
        }
    }
    return false;
}

} // namespace detail_depict

// ════════════════════════════════════════════════════════════════════════════
//  SVG Renderer — Publication Quality
// ════════════════════════════════════════════════════════════════════════════

struct DepictResult {
    std::string svg;
    int width, height;
};

inline DepictResult depictSVG(
    const MolGraph& g,
    const std::vector<Vec2>& pos,
    const std::map<int,int>& mapping = {},          // atom index -> map number
    const std::set<int>& highlightAtoms = {},
    const std::set<std::pair<int,int>>& highlightBonds = {},
    const DepictOptions& opts = {}
) {
    if (g.n == 0) return {"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1\" height=\"1\"/>", 1, 1};

    // Precompute rings for double bond offset
    auto rings = findLayoutRings(g);

    // Build atom labels
    std::vector<detail_depict::AtomLabel> labels(g.n);
    for (int i = 0; i < g.n; i++)
        labels[i] = detail_depict::buildLabel(g, i, opts, pos);

    // Resolved style values
    const double LW    = opts.lineWidth_();
    const double BW    = opts.boldWidth_();
    const double BS    = opts.bondSpacing_();
    const double FS    = opts.fontSize_();
    const double HR    = opts.highlightRadius_();
    const double MW    = opts.matchWidth_();
    const double MNS   = opts.mapNumSize_();
    const double SUB   = opts.subscriptScale;

    // Bounding box
    double minX = 1e9, minY = 1e9, maxX = -1e9, maxY = -1e9;
    for (int i = 0; i < g.n; i++) {
        double extra = labels[i].visible ? labels[i].textWidth * 0.5 + 4 : 0;
        minX = std::min(minX, pos[i].x - extra);
        minY = std::min(minY, pos[i].y - FS * 0.6);
        maxX = std::max(maxX, pos[i].x + extra);
        maxY = std::max(maxY, pos[i].y + FS * 0.6);
    }
    double pad = opts.padding;
    int W = opts.width  > 0 ? opts.width  : std::max(100, (int)(maxX - minX + 2 * pad));
    int H = opts.height > 0 ? opts.height : std::max(100, (int)(maxY - minY + 2 * pad));

    // If user specifies fixed size, center the drawing
    double offsetX = pad - minX, offsetY = pad - minY;
    if (opts.width > 0) offsetX = (opts.width  - (maxX - minX)) * 0.5 - minX;
    if (opts.height > 0) offsetY = (opts.height - (maxY - minY)) * 0.5 - minY;

    auto tx = [&](double x) { return x + offsetX; };
    auto ty = [&](double y) { return y + offsetY; };

    // Get atom color (respecting overrides)
    auto getAtomColor = [&](int idx) -> Color {
        int Z = g.atomicNum[idx];
        auto it = opts.atomColorOverrides.find(Z);
        if (it != opts.atomColorOverrides.end()) return it->second;
        return atomColor(Z);
    };

    // Clip bond endpoint at label boundary
    auto clipBond = [&](int atomIdx, Vec2 from, Vec2 to) -> Vec2 {
        if (!labels[atomIdx].visible) return to;
        double halfW = labels[atomIdx].textWidth * 0.5 + 2;
        double halfH = FS * 0.55;
        Vec2 d = from - to;  // direction from endpoint toward other atom
        double dl = len(d);
        if (dl < 1e-9) return to;
        // Move endpoint outward by the label half-extent along bond direction
        double ratio = std::min(1.0, std::max(halfW, halfH) / dl);
        return to + d * ratio;
    };

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);

    // SVG header with precision rendering hints
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\"";
    svg << " viewBox=\"0 0 " << W << " " << H << "\"";
    svg << " shape-rendering=\"geometricPrecision\"";
    svg << " text-rendering=\"geometricPrecision\">\n";

    // Defs: styles and markers
    svg << "<defs>\n";
    svg << "  <style>\n";
    svg << "    text { font-family: " << opts.fontFamily << "; }\n";
    svg << "    .atom { dominant-baseline: central; text-anchor: middle; font-weight: bold; }\n";
    svg << "    .sub  { dominant-baseline: central; font-weight: bold; }\n";
    svg << "    .map  { dominant-baseline: central; text-anchor: start; font-weight: bold; }\n";
    svg << "    .idx  { dominant-baseline: central; text-anchor: start; font-weight: normal; }\n";
    svg << "  </style>\n";
    svg << "</defs>\n";

    // Background
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"" << rgb(opts.backgroundColor) << "\"/>\n";

    // ── Layer 1: Atom highlight fills ──────────────────────────────────────
    for (int i = 0; i < g.n; i++) {
        if (!highlightAtoms.count(i)) continue;
        Color hc = opts.matchFill;
        auto cit = opts.atomHighlightColors.find(i);
        if (cit != opts.atomHighlightColors.end()) hc = cit->second;
        svg << "<circle cx=\"" << tx(pos[i].x) << "\" cy=\"" << ty(pos[i].y)
            << "\" r=\"" << HR << "\" fill=\"" << rgba(hc) << "\" stroke=\"none\"/>\n";
    }

    // ── Layer 2: Bonds ────────────────────────────────────────────────────
    std::set<std::pair<int,int>> drawnBonds;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue;
            if (drawnBonds.count({i,j})) continue;
            drawnBonds.insert({i,j});

            Vec2 p1 = {tx(pos[i].x), ty(pos[i].y)};
            Vec2 p2 = {tx(pos[j].x), ty(pos[j].y)};

            // Clip at label boundaries
            Vec2 cp1 = clipBond(i, p2, p1);
            Vec2 cp2 = clipBond(j, p1, p2);

            bool isHighlighted = highlightBonds.count({i,j}) || highlightBonds.count({j,i})
                || (highlightAtoms.count(i) && highlightAtoms.count(j));

            Color bondCol = opts.bondColor;
            double width = LW;
            if (isHighlighted) {
                bondCol = opts.matchStroke;
                width = MW;
                auto bit = opts.bondHighlightColors.find({i,j});
                if (bit == opts.bondHighlightColors.end())
                    bit = opts.bondHighlightColors.find({j,i});
                if (bit != opts.bondHighlightColors.end()) bondCol = bit->second;
            }

            int order = g.bondOrder(i, j);
            bool arom = g.bondAromatic(i, j);

            // Check for stereo bond type
            BondStereoType stereo = BondStereoType::NONE;
            if (!g.tetraChirality.empty()) {
                // Wedge: if atom i is a stereo center and j is a substituent
                if (g.tetraChirality[i] != 0 && order == 1) {
                    // Assign wedge direction based on neighbor ordering
                    int nbIdx = 0;
                    for (int k = 0; k < (int)g.neighbors[i].size(); k++) {
                        if (g.neighbors[i][k] == j) { nbIdx = k; break; }
                    }
                    if (nbIdx == 0) stereo = BondStereoType::WEDGE_UP;
                    else if (nbIdx == 1) stereo = BondStereoType::WEDGE_DOWN;
                }
            }

            Vec2 bondDir = norm(cp2 - cp1);
            Vec2 bondPerp = perp(bondDir);

            if (stereo == BondStereoType::WEDGE_UP) {
                // Solid filled wedge: narrow at stereo center, wide at substituent
                Vec2 w = bondPerp * (BW * 0.5);
                svg << "<polygon points=\""
                    << cp1.x << "," << cp1.y << " "
                    << (cp2.x + w.x) << "," << (cp2.y + w.y) << " "
                    << (cp2.x - w.x) << "," << (cp2.y - w.y)
                    << "\" fill=\"" << rgb(bondCol) << "\" stroke=\"none\"/>\n";
            }
            else if (stereo == BondStereoType::WEDGE_DOWN) {
                // Dashed wedge: series of stripes widening from stereo center
                int nStripes = std::max(4, (int)(len(cp2 - cp1) / opts.hashSpacing_()));
                for (int s = 0; s < nStripes; s++) {
                    double t = (s + 0.5) / nStripes;
                    Vec2 mid = cp1 + (cp2 - cp1) * t;
                    double halfW = BW * 0.5 * t;
                    Vec2 w = bondPerp * halfW;
                    svg << "<line x1=\"" << (mid.x - w.x) << "\" y1=\"" << (mid.y - w.y)
                        << "\" x2=\"" << (mid.x + w.x) << "\" y2=\"" << (mid.y + w.y)
                        << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << LW
                        << "\" stroke-linecap=\"round\"/>\n";
                }
            }
            else if (order == 1 || order == 0 || (arom && opts.aromaticStyle == AromaticStyle::CIRCLE)) {
                // Single bond
                svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                    << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                    << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
            }
            else if (order == 2) {
                // Double bond — asymmetric toward ring interior, symmetric otherwise
                Vec2 ringC;
                bool inRing = detail_depict::bondInRing(i, j, rings, pos, ringC);

                if (inRing) {
                    // Full line on bond axis + shorter offset line toward ring center
                    svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                        << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                        << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                        << "\" stroke-linecap=\"round\"/>\n";

                    // Offset toward ring center
                    Vec2 bondMid = (p1 + p2) * 0.5;
                    Vec2 toCenter = norm(Vec2{tx(ringC.x), ty(ringC.y)} - bondMid);
                    double offsetDist = BS;
                    // Ensure offset is toward the ring (positive dot with toCenter)
                    if (dot(bondPerp, toCenter) < 0) offsetDist = -offsetDist;
                    Vec2 off = bondPerp * offsetDist;

                    // Shorten inner line by 15% on each end
                    Vec2 innerDir = bondDir;
                    double bondLen = len(cp2 - cp1);
                    double trim = bondLen * 0.15;
                    Vec2 ip1 = cp1 + innerDir * trim + off;
                    Vec2 ip2 = cp2 - innerDir * trim + off;

                    svg << "<line x1=\"" << ip1.x << "\" y1=\"" << ip1.y
                        << "\" x2=\"" << ip2.x << "\" y2=\"" << ip2.y
                        << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                        << "\" stroke-linecap=\"round\"/>\n";
                } else {
                    // Symmetric double bond
                    Vec2 off = bondPerp * (BS * 0.5);
                    svg << "<line x1=\"" << (cp1.x+off.x) << "\" y1=\"" << (cp1.y+off.y)
                        << "\" x2=\"" << (cp2.x+off.x) << "\" y2=\"" << (cp2.y+off.y)
                        << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                        << "\" stroke-linecap=\"round\"/>\n";
                    svg << "<line x1=\"" << (cp1.x-off.x) << "\" y1=\"" << (cp1.y-off.y)
                        << "\" x2=\"" << (cp2.x-off.x) << "\" y2=\"" << (cp2.y-off.y)
                        << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                        << "\" stroke-linecap=\"round\"/>\n";
                }
            }
            else if (order == 3) {
                // Triple bond: center + two offset lines
                Vec2 off = bondPerp * BS;
                svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                    << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                    << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width
                    << "\" stroke-linecap=\"round\"/>\n";
                svg << "<line x1=\"" << (cp1.x+off.x) << "\" y1=\"" << (cp1.y+off.y)
                    << "\" x2=\"" << (cp2.x+off.x) << "\" y2=\"" << (cp2.y+off.y)
                    << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width * 0.7
                    << "\" stroke-linecap=\"round\"/>\n";
                svg << "<line x1=\"" << (cp1.x-off.x) << "\" y1=\"" << (cp1.y-off.y)
                    << "\" x2=\"" << (cp2.x-off.x) << "\" y2=\"" << (cp2.y-off.y)
                    << "\" stroke=\"" << rgb(bondCol) << "\" stroke-width=\"" << width * 0.7
                    << "\" stroke-linecap=\"round\"/>\n";
            }
        }
    }

    // ── Layer 3: Aromatic ring circles ────────────────────────────────────
    if (opts.aromaticStyle == AromaticStyle::CIRCLE) {
        for (auto& ring : rings) {
            bool allArom = true;
            for (int a : ring) if (!g.aromatic[a]) { allArom = false; break; }
            if (!allArom) continue;
            Vec2 c = {0, 0};
            for (int a : ring) { c.x += tx(pos[a].x); c.y += ty(pos[a].y); }
            c.x /= ring.size(); c.y /= ring.size();
            double R = 0;
            for (int a : ring) {
                double dx = tx(pos[a].x) - c.x, dy = ty(pos[a].y) - c.y;
                R += std::sqrt(dx * dx + dy * dy);
            }
            R = R / ring.size() * 0.55;
            svg << "<circle cx=\"" << c.x << "\" cy=\"" << c.y << "\" r=\"" << R
                << "\" fill=\"none\" stroke=\"" << rgb(opts.bondColor)
                << "\" stroke-width=\"" << LW * 0.7 << "\"/>\n";
        }
    }

    // ── Layer 4: Label backgrounds (white masks) ──────────────────────────
    for (int i = 0; i < g.n; i++) {
        if (!labels[i].visible) continue;
        double x = tx(pos[i].x), y = ty(pos[i].y);
        double halfW = labels[i].textWidth * 0.5 + 2;
        double halfH = FS * 0.55;
        svg << "<rect x=\"" << (x - halfW) << "\" y=\"" << (y - halfH)
            << "\" width=\"" << (halfW * 2) << "\" height=\"" << (halfH * 2)
            << "\" fill=\"" << rgb(opts.backgroundColor)
            << "\" stroke=\"none\" rx=\"1\"/>\n";
    }

    // ── Layer 5: Atom labels with H-count and charges ─────────────────────
    for (int i = 0; i < g.n; i++) {
        if (!labels[i].visible) continue;
        double x = tx(pos[i].x), y = ty(pos[i].y);
        Color col = getAtomColor(i);
        auto& lbl = labels[i];

        // Build the label as positioned SVG text elements
        double subFS = FS * SUB;

        if (lbl.hCount > 0 && lbl.hLeft) {
            // H first: "H₂N" or "HO"
            double cursor = x - lbl.textWidth * 0.5;

            // H
            svg << "<text x=\"" << cursor << "\" y=\"" << y
                << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                << "\" fill=\"" << rgb(col) << "\">H</text>\n";
            cursor += FS * 0.6;

            // Subscript
            if (lbl.hCount > 1) {
                svg << "<text x=\"" << cursor << "\" y=\"" << (y + FS * 0.25)
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                    << "\" fill=\"" << rgb(col) << "\">" << lbl.hCount << "</text>\n";
                cursor += subFS * 0.5;
            }

            // Element symbol
            svg << "<text x=\"" << cursor << "\" y=\"" << y
                << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                << "\" fill=\"" << rgb(col) << "\">" << lbl.symbol << "</text>\n";
            cursor += lbl.symbol.length() * FS * 0.6;

            // Charge superscript
            if (lbl.charge != 0) {
                std::string chStr;
                if (lbl.charge == 1) chStr = "+";
                else if (lbl.charge == -1) chStr = "\xe2\x88\x92";  // minus sign
                else if (lbl.charge > 1) chStr = std::to_string(lbl.charge) + "+";
                else chStr = std::to_string(-lbl.charge) + "\xe2\x88\x92";
                svg << "<text x=\"" << cursor << "\" y=\"" << (y - FS * 0.3)
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                    << "\" fill=\"" << rgb(col) << "\">" << chStr << "</text>\n";
            }
        } else {
            // Element first: "NH₂" or "OH" or "S" (standard)
            double cursor = x - lbl.textWidth * 0.5;

            // Element symbol
            svg << "<text x=\"" << cursor << "\" y=\"" << y
                << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                << "\" fill=\"" << rgb(col) << "\">" << lbl.symbol << "</text>\n";
            cursor += lbl.symbol.length() * FS * 0.6;

            // H + subscript
            if (lbl.hCount > 0) {
                svg << "<text x=\"" << cursor << "\" y=\"" << y
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                    << "\" fill=\"" << rgb(col) << "\">H</text>\n";
                cursor += FS * 0.6;
                if (lbl.hCount > 1) {
                    svg << "<text x=\"" << cursor << "\" y=\"" << (y + FS * 0.25)
                        << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                        << "\" fill=\"" << rgb(col) << "\">" << lbl.hCount << "</text>\n";
                    cursor += subFS * 0.5;
                }
            }

            // Charge superscript
            if (lbl.charge != 0) {
                std::string chStr;
                if (lbl.charge == 1) chStr = "+";
                else if (lbl.charge == -1) chStr = "\xe2\x88\x92";
                else if (lbl.charge > 1) chStr = std::to_string(lbl.charge) + "+";
                else chStr = std::to_string(-lbl.charge) + "\xe2\x88\x92";
                svg << "<text x=\"" << cursor << "\" y=\"" << (y - FS * 0.3)
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                    << "\" fill=\"" << rgb(col) << "\">" << chStr << "</text>\n";
            }
        }

        // Atom map number (superscript, top-right, blue)
        if (opts.showMapNumbers && mapping.count(i)) {
            svg << "<text x=\"" << (x + lbl.textWidth * 0.5 + 2) << "\" y=\"" << (y - FS * 0.5)
                << "\" class=\"map\" font-size=\"" << MNS
                << "\" fill=\"" << rgb(palette::MAPNUM) << "\">" << mapping.at(i) << "</text>\n";
        }

        // Atom index (subscript, bottom-right, gray)
        if (opts.showAtomIndices) {
            svg << "<text x=\"" << (x + lbl.textWidth * 0.5 + 2) << "\" y=\"" << (y + FS * 0.5)
                << "\" class=\"idx\" font-size=\"" << MNS
                << "\" fill=\"rgba(120,120,120,0.6)\">" << i << "</text>\n";
        }
    }

    // Atom indices for carbon atoms (shown without label)
    if (opts.showAtomIndices) {
        for (int i = 0; i < g.n; i++) {
            if (labels[i].visible) continue;  // already drawn above
            double x = tx(pos[i].x), y = ty(pos[i].y);
            svg << "<text x=\"" << (x + 4) << "\" y=\"" << (y + FS * 0.4)
                << "\" class=\"idx\" font-size=\"" << MNS
                << "\" fill=\"rgba(120,120,120,0.5)\">" << i << "</text>\n";
        }
    }

    svg << "</svg>\n";
    return {svg.str(), W, H};
}

// ════════════════════════════════════════════════════════════════════════════
//  High-Level API
// ════════════════════════════════════════════════════════════════════════════

// Depict a single molecule → SVG string
inline std::string depict(const MolGraph& g, const DepictOptions& opts = {}) {
    auto pos = layout2D(g, opts);
    return depictSVG(g, pos, {}, {}, {}, opts).svg;
}

// Depict with MCS/substructure highlight (one molecule, matched atoms marked)
inline std::string depictWithMapping(
    const MolGraph& g,
    const std::map<int,int>& atomMapping,
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
    // Highlighted bonds: if both endpoints are matched
    std::set<std::pair<int,int>> hBonds;
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue;
            if (highlighted.count(i) && highlighted.count(j))
                hBonds.insert({i, j});
        }
    }
    return depictSVG(g, pos, mapNums, highlighted, hBonds, opts).svg;
}

// Depict a pair side-by-side with MCS/substructure highlighting
inline std::string depictPair(
    const MolGraph& g1, const MolGraph& g2,
    const std::map<int,int>& mapping,  // g1 atom → g2 atom
    const DepictOptions& opts = {}
) {
    auto pos1 = layout2D(g1, opts);
    auto pos2 = layout2D(g2, opts);

    // Offset g2 to the right
    double maxX1 = -1e9;
    for (int i = 0; i < g1.n; i++) maxX1 = std::max(maxX1, pos1[i].x);
    double minX2 = 1e9;
    for (int i = 0; i < g2.n; i++) minX2 = std::min(minX2, pos2[i].x);
    double offset = maxX1 - minX2 + opts.bondLength * 3;
    for (int i = 0; i < g2.n; i++) pos2[i].x += offset;

    // Build combined for rendering
    int n = g1.n + g2.n;
    std::vector<Vec2> pos(n);
    for (int i = 0; i < g1.n; i++) pos[i] = pos1[i];
    for (int i = 0; i < g2.n; i++) pos[g1.n + i] = pos2[i];

    // Highlights
    std::set<int> hlAtoms;
    std::map<int,int> mapNums;
    int mapNum = 1;
    for (auto& [qi, tj] : mapping) {
        hlAtoms.insert(qi);
        hlAtoms.insert(g1.n + tj);
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
    int W = opts.width  > 0 ? opts.width  : std::max(100, (int)(maxX - minX + 2 * pad));
    int H = opts.height > 0 ? opts.height : std::max(100, (int)(maxY - minY + 2 * pad));

    double offsetX = pad - minX, offsetY = pad - minY;
    if (opts.width > 0) offsetX = (opts.width - (maxX - minX)) * 0.5 - minX;
    if (opts.height > 0) offsetY = (opts.height - (maxY - minY)) * 0.5 - minY;
    auto tx = [&](double x) { return x + offsetX; };
    auto ty = [&](double y) { return y + offsetY; };

    // Resolved values
    const double LW  = opts.lineWidth_();
    const double BS  = opts.bondSpacing_();
    const double FS  = opts.fontSize_();
    const double HR  = opts.highlightRadius_();
    const double MW  = opts.matchWidth_();
    const double MNS = opts.mapNumSize_();
    const double BW  = opts.boldWidth_();
    const double SUB = opts.subscriptScale;

    auto rings1 = findLayoutRings(g1);
    auto rings2 = findLayoutRings(g2);

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\"";
    svg << " viewBox=\"0 0 " << W << " " << H << "\"";
    svg << " shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n";
    svg << "<defs><style>\n";
    svg << "  text { font-family: " << opts.fontFamily << "; }\n";
    svg << "  .atom { dominant-baseline: central; text-anchor: middle; font-weight: bold; }\n";
    svg << "  .sub  { dominant-baseline: central; font-weight: bold; }\n";
    svg << "  .map  { dominant-baseline: central; text-anchor: start; font-weight: bold; }\n";
    svg << "</style></defs>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"" << rgb(opts.backgroundColor) << "\"/>\n";

    // Lambda to draw one molecule into the combined SVG
    auto drawMol = [&](const MolGraph& g, int globalOffset,
                       const std::vector<std::vector<int>>& rings) {
        // Build labels
        std::vector<detail_depict::AtomLabel> labels(g.n);
        std::vector<Vec2> localPos(g.n);
        for (int i = 0; i < g.n; i++) localPos[i] = pos[globalOffset + i];
        for (int i = 0; i < g.n; i++)
            labels[i] = detail_depict::buildLabel(g, i, opts, localPos);

        auto clipBond = [&](int atomIdx, Vec2 from, Vec2 to) -> Vec2 {
            if (!labels[atomIdx].visible) return to;
            double halfW = labels[atomIdx].textWidth * 0.5 + 2;
            double halfH = FS * 0.55;
            Vec2 d = from - to;
            double dl = len(d);
            if (dl < 1e-9) return to;
            double ratio = std::min(1.0, std::max(halfW, halfH) / dl);
            return to + d * ratio;
        };

        auto getColor = [&](int idx) -> Color {
            auto it = opts.atomColorOverrides.find(g.atomicNum[idx]);
            if (it != opts.atomColorOverrides.end()) return it->second;
            return atomColor(g.atomicNum[idx]);
        };

        // Highlight fills
        for (int i = 0; i < g.n; i++) {
            int gi = globalOffset + i;
            if (!hlAtoms.count(gi)) continue;
            svg << "<circle cx=\"" << tx(pos[gi].x) << "\" cy=\"" << ty(pos[gi].y)
                << "\" r=\"" << HR << "\" fill=\"" << rgba(opts.matchFill) << "\" stroke=\"none\"/>\n";
        }

        // Bonds
        for (int i = 0; i < g.n; i++) {
            for (int j : g.neighbors[i]) {
                if (j <= i) continue;
                int gi = globalOffset + i, gj = globalOffset + j;
                Vec2 p1 = {tx(pos[gi].x), ty(pos[gi].y)};
                Vec2 p2 = {tx(pos[gj].x), ty(pos[gj].y)};
                Vec2 cp1 = clipBond(i, p2, p1);
                Vec2 cp2 = clipBond(j, p1, p2);

                bool hl = hlAtoms.count(gi) && hlAtoms.count(gj);
                Color c = hl ? opts.matchStroke : opts.bondColor;
                double w = hl ? MW : LW;
                int order = g.bondOrder(i, j);
                bool arom = g.bondAromatic(i, j);
                Vec2 bondDir = norm(cp2 - cp1);
                Vec2 bondPerp = perp(bondDir);

                if (order <= 1 || order == 0 || (arom && opts.aromaticStyle == AromaticStyle::CIRCLE)) {
                    svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                        << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                        << "\" stroke-linecap=\"round\"/>\n";
                } else if (order == 2) {
                    Vec2 ringC;
                    bool inRing = detail_depict::bondInRing(i, j, rings, localPos, ringC);
                    if (inRing) {
                        svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                            << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                            << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                            << "\" stroke-linecap=\"round\"/>\n";
                        Vec2 bondMid = (p1 + p2) * 0.5;
                        Vec2 toCenter = norm(Vec2{tx(ringC.x), ty(ringC.y)} - bondMid);
                        double offsetDist = BS;
                        if (dot(bondPerp, toCenter) < 0) offsetDist = -offsetDist;
                        Vec2 off = bondPerp * offsetDist;
                        double bondLen = len(cp2 - cp1);
                        double trim = bondLen * 0.15;
                        Vec2 ip1 = cp1 + bondDir * trim + off;
                        Vec2 ip2 = cp2 - bondDir * trim + off;
                        svg << "<line x1=\"" << ip1.x << "\" y1=\"" << ip1.y
                            << "\" x2=\"" << ip2.x << "\" y2=\"" << ip2.y
                            << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                            << "\" stroke-linecap=\"round\"/>\n";
                    } else {
                        Vec2 off = bondPerp * (BS * 0.5);
                        svg << "<line x1=\"" << (cp1.x+off.x) << "\" y1=\"" << (cp1.y+off.y)
                            << "\" x2=\"" << (cp2.x+off.x) << "\" y2=\"" << (cp2.y+off.y)
                            << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                            << "\" stroke-linecap=\"round\"/>\n";
                        svg << "<line x1=\"" << (cp1.x-off.x) << "\" y1=\"" << (cp1.y-off.y)
                            << "\" x2=\"" << (cp2.x-off.x) << "\" y2=\"" << (cp2.y-off.y)
                            << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                            << "\" stroke-linecap=\"round\"/>\n";
                    }
                } else if (order == 3) {
                    Vec2 off = bondPerp * BS;
                    svg << "<line x1=\"" << cp1.x << "\" y1=\"" << cp1.y
                        << "\" x2=\"" << cp2.x << "\" y2=\"" << cp2.y
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w
                        << "\" stroke-linecap=\"round\"/>\n";
                    svg << "<line x1=\"" << (cp1.x+off.x) << "\" y1=\"" << (cp1.y+off.y)
                        << "\" x2=\"" << (cp2.x+off.x) << "\" y2=\"" << (cp2.y+off.y)
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w*0.7
                        << "\" stroke-linecap=\"round\"/>\n";
                    svg << "<line x1=\"" << (cp1.x-off.x) << "\" y1=\"" << (cp1.y-off.y)
                        << "\" x2=\"" << (cp2.x-off.x) << "\" y2=\"" << (cp2.y-off.y)
                        << "\" stroke=\"" << rgb(c) << "\" stroke-width=\"" << w*0.7
                        << "\" stroke-linecap=\"round\"/>\n";
                }
            }
        }

        // Aromatic circles
        if (opts.aromaticStyle == AromaticStyle::CIRCLE) {
            for (auto& ring : rings) {
                bool allArom = true;
                for (int a : ring) if (!g.aromatic[a]) { allArom = false; break; }
                if (!allArom) continue;
                Vec2 c = {0, 0};
                for (int a : ring) { c.x += tx(pos[globalOffset+a].x); c.y += ty(pos[globalOffset+a].y); }
                c.x /= ring.size(); c.y /= ring.size();
                double R = 0;
                for (int a : ring) {
                    double dx = tx(pos[globalOffset+a].x) - c.x;
                    double dy = ty(pos[globalOffset+a].y) - c.y;
                    R += std::sqrt(dx*dx + dy*dy);
                }
                R = R / ring.size() * 0.55;
                svg << "<circle cx=\"" << c.x << "\" cy=\"" << c.y << "\" r=\"" << R
                    << "\" fill=\"none\" stroke=\"" << rgb(opts.bondColor)
                    << "\" stroke-width=\"" << LW * 0.7 << "\"/>\n";
            }
        }

        // Label backgrounds
        for (int i = 0; i < g.n; i++) {
            if (!labels[i].visible) continue;
            int gi = globalOffset + i;
            double x = tx(pos[gi].x), y = ty(pos[gi].y);
            double halfW = labels[i].textWidth * 0.5 + 2;
            double halfH = FS * 0.55;
            svg << "<rect x=\"" << (x-halfW) << "\" y=\"" << (y-halfH)
                << "\" width=\"" << (halfW*2) << "\" height=\"" << (halfH*2)
                << "\" fill=\"" << rgb(opts.backgroundColor) << "\" stroke=\"none\" rx=\"1\"/>\n";
        }

        // Atom labels
        for (int i = 0; i < g.n; i++) {
            if (!labels[i].visible) continue;
            int gi = globalOffset + i;
            double x = tx(pos[gi].x), y = ty(pos[gi].y);
            Color col = getColor(i);
            auto& lbl = labels[i];
            double subFS = FS * SUB;
            double cursor = x - lbl.textWidth * 0.5;

            if (lbl.hCount > 0 && lbl.hLeft) {
                svg << "<text x=\"" << cursor << "\" y=\"" << y
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                    << "\" fill=\"" << rgb(col) << "\">H</text>\n";
                cursor += FS * 0.6;
                if (lbl.hCount > 1) {
                    svg << "<text x=\"" << cursor << "\" y=\"" << (y + FS*0.25)
                        << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                        << "\" fill=\"" << rgb(col) << "\">" << lbl.hCount << "</text>\n";
                    cursor += subFS * 0.5;
                }
                svg << "<text x=\"" << cursor << "\" y=\"" << y
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                    << "\" fill=\"" << rgb(col) << "\">" << lbl.symbol << "</text>\n";
                cursor += lbl.symbol.length() * FS * 0.6;
            } else {
                svg << "<text x=\"" << cursor << "\" y=\"" << y
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                    << "\" fill=\"" << rgb(col) << "\">" << lbl.symbol << "</text>\n";
                cursor += lbl.symbol.length() * FS * 0.6;
                if (lbl.hCount > 0) {
                    svg << "<text x=\"" << cursor << "\" y=\"" << y
                        << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << FS
                        << "\" fill=\"" << rgb(col) << "\">H</text>\n";
                    cursor += FS * 0.6;
                    if (lbl.hCount > 1) {
                        svg << "<text x=\"" << cursor << "\" y=\"" << (y + FS*0.25)
                            << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                            << "\" fill=\"" << rgb(col) << "\">" << lbl.hCount << "</text>\n";
                        cursor += subFS * 0.5;
                    }
                }
            }

            if (lbl.charge != 0) {
                std::string chStr;
                if (lbl.charge == 1) chStr = "+";
                else if (lbl.charge == -1) chStr = "\xe2\x88\x92";
                else if (lbl.charge > 1) chStr = std::to_string(lbl.charge) + "+";
                else chStr = std::to_string(-lbl.charge) + "\xe2\x88\x92";
                svg << "<text x=\"" << cursor << "\" y=\"" << (y - FS*0.3)
                    << "\" class=\"sub\" text-anchor=\"start\" font-size=\"" << subFS
                    << "\" fill=\"" << rgb(col) << "\">" << chStr << "</text>\n";
            }

            if (opts.showMapNumbers && mapNums.count(gi)) {
                svg << "<text x=\"" << (x + lbl.textWidth*0.5 + 2) << "\" y=\"" << (y - FS*0.5)
                    << "\" class=\"map\" font-size=\"" << MNS
                    << "\" fill=\"" << rgb(palette::MAPNUM) << "\">" << mapNums[gi] << "</text>\n";
            }
        }
    };

    drawMol(g1, 0, rings1);

    // Separator: bidirectional arrow
    double sepX = (g1.n > 0 && g2.n > 0)
        ? tx((pos[g1.n - 1].x + pos[g1.n].x) / 2)
        : W / 2.0;
    double sepY = H / 2.0;
    // Draw a clean arrow line
    double arrowLen = opts.bondLength * 1.2;
    svg << "<line x1=\"" << (sepX - arrowLen*0.5) << "\" y1=\"" << sepY
        << "\" x2=\"" << (sepX + arrowLen*0.5) << "\" y2=\"" << sepY
        << "\" stroke=\"#999\" stroke-width=\"" << LW << "\" stroke-linecap=\"round\"/>\n";
    // Arrowheads
    double ah = opts.bondLength * 0.15;
    svg << "<polyline points=\""
        << (sepX + arrowLen*0.5 - ah) << "," << (sepY - ah) << " "
        << (sepX + arrowLen*0.5) << "," << sepY << " "
        << (sepX + arrowLen*0.5 - ah) << "," << (sepY + ah)
        << "\" fill=\"none\" stroke=\"#999\" stroke-width=\"" << LW
        << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n";
    svg << "<polyline points=\""
        << (sepX - arrowLen*0.5 + ah) << "," << (sepY - ah) << " "
        << (sepX - arrowLen*0.5) << "," << sepY << " "
        << (sepX - arrowLen*0.5 + ah) << "," << (sepY + ah)
        << "\" fill=\"none\" stroke=\"#999\" stroke-width=\"" << LW
        << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\"/>\n";

    drawMol(g2, g1.n, rings2);

    svg << "</svg>\n";
    return svg.str();
}

// ── Convenience: SMILES → SVG ─────────────────────────────────────────────
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
    return depictPair(parseSMILES(smiles1), parseSMILES(smiles2), mapping, opts);
}

#endif // SMSD_SMILES_PARSER_HPP

} // namespace smsd

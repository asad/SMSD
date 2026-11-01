/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * SMSD Pro — Publication-Quality Molecular Depiction Engine
 * Zero-dependency SVG renderer: ACS 1996, Nature, Springer journal standard
 * Accepts MolGraph, SMILES, substructure/MCS mappings; auto-layout when needed
 */
#pragma once
#include "smsd/mol_graph.hpp"
#include "smsd/mol_reader.hpp"
#include "smsd/layout.hpp"
#include "smsd/smiles_parser.hpp"
#include <cmath>
#include <sstream>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <iomanip>

namespace smsd {

namespace detail_depict {
inline constexpr double PI = 3.14159265358979323846;
}

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
//  Colour
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
//  Jmol/CPK Element Colour Palette (publication standard)
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
    // Highlight colours: a soft halo behind an unchanged skeleton.  Gold and mint are equally
    // light, so either can mark any match; gold is the default, and mint suits substructure hits.
    constexpr Color GOLD_HALO    = {240,192, 64, 0.40};  // brand gold
    constexpr Color MINT_HALO    = { 78,232,199, 0.40};  // brand mint
    constexpr Color MATCH_FILL   = GOLD_HALO;
    constexpr Color MATCH_STROKE = { 33, 33, 33, 1.0};   // matched bonds keep the default bond colour
    constexpr Color UNMATCH_FILL = {158,158,158, 0.10};   // very light grey
    // Structural colours
    constexpr Color BOND         = { 33, 33, 33, 1.0};   // near-black
    constexpr Color MAPNUM       = { 86,104,101, 1.0};   // muted grey-green
    constexpr Color ARROW        = {128,146,143, 1.0};   // pair arrow
    constexpr Color BACKGROUND   = {255,255,255, 1.0};   // white
    // Highlight colour set for multi-fragment / multi-match.  Every halo differs from
    // the others in hue or lightness, so the set stays apart for colour-blind readers.
    constexpr Color HIGHLIGHT_1  = GOLD_HALO;
    constexpr Color HIGHLIGHT_2  = {120,184,240, 0.50};  // sky
    constexpr Color HIGHLIGHT_3  = {232,108, 90, 0.60};  // coral
    constexpr Color HIGHLIGHT_4  = { 13, 51, 48, 0.25};  // deep teal
    constexpr Color HIGHLIGHT_5  = { 34,161,150, 0.75};  // bright teal
    constexpr Color HIGHLIGHT_6  = MINT_HALO;
}

// ── Element label colours: the customary hues, at 4.5:1 or more on white and on both halos ─
namespace ink {
    constexpr Color C     = { 33, 33, 33, 1.0};
    constexpr Color N     = { 36, 72,196, 1.0};
    constexpr Color O     = {192, 34, 28, 1.0};
    constexpr Color F     = { 42,100, 16, 1.0};
    constexpr Color Si    = {122, 90, 60, 1.0};
    constexpr Color P     = {154, 60,  0, 1.0};
    constexpr Color S     = {115, 80,  0, 1.0};
    constexpr Color Cl    = { 19,100, 52, 1.0};
    constexpr Color Se    = {146, 84,  0, 1.0};
    constexpr Color Br    = {150, 40, 34, 1.0};
    constexpr Color I     = {120, 36,150, 1.0};
    constexpr Color B     = {166, 70, 62, 1.0};
    constexpr Color Na    = {122, 72,184, 1.0};
    constexpr Color Mg    = { 74,104,  0, 1.0};
    constexpr Color Ca    = { 40,108, 26, 1.0};
    constexpr Color Fe    = {166, 66, 24, 1.0};
    constexpr Color Cu    = {138, 80, 24, 1.0};
    constexpr Color Zn    = { 92, 96,150, 1.0};
    constexpr Color OTHER = { 74, 96,108, 1.0};  // blue-grey
}

// Element label colour by atomic number
inline Color atomColor(int Z) {
    switch (Z) {
        case  1: return ink::C;
        case  5: return ink::B;
        case  6: return ink::C;
        case  7: return ink::N;
        case  8: return ink::O;
        case  9: return ink::F;
        case 14: return ink::Si;
        case 15: return ink::P;
        case 16: return ink::S;
        case 17: return ink::Cl;
        case 35: return ink::Br;
        case 34: return ink::Se;
        case 53: return ink::I;
        case 11: return ink::Na;
        case 12: return ink::Mg;
        case 20: return ink::Ca;
        case 26: return ink::Fe;
        case 29: return ink::Cu;
        case 30: return ink::Zn;
        default: return ink::OTHER;
    }
}

/** Label symbol of an element from hydrogen to radon, from the table of smiles_utils.hpp; "?" otherwise. */
inline const char* elementSymbol(int Z) {
    return (Z >= 1 && Z <= 86) ? detail::atomicNumToSymbol(Z) : "?";
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
//  Depiction Options (everything customisable)
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
    double highlightRadius = 0.0;     // 0 = auto (32% of bondLength)
    double highlightOpacity = 0.30;
    double matchWidth      = 0.0;     // 0 = auto (lineWidth)

    // ── Aromatic rendering ──
    AromaticStyle aromaticStyle = AromaticStyle::KEKULE;

    // ── Map numbers ──
    double mapNumSize      = 0.0;     // 0 = auto
    bool   showMapNumbers  = true;

    // ── SVG / canvas ──
    double padding         = 30.0;
    int    width           = 0;       // 0 = auto-fit
    int    height          = 0;       // 0 = auto-fit
    std::string fontFamily = "Arial, Helvetica, sans-serif";
    Color  backgroundColor = palette::BACKGROUND;

    // ── Custom colours (override defaults) ──
    Color  bondColor       = palette::BOND;
    Color  matchFill       = palette::MATCH_FILL;
    Color  matchStroke     = palette::MATCH_STROKE;
    std::map<int, Color> atomColorOverrides;   // Z -> Colour
    std::map<int, Color> atomHighlightColors;  // atom index -> Colour
    std::map<std::pair<int,int>, Color> bondHighlightColors; // bond -> Colour

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
        return mapNumSize > 0 ? mapNumSize : fontSize_() * 0.55;
    }
    double highlightRadius_() const {
        return highlightRadius > 0 ? highlightRadius : bondLength * 0.32;
    }
    double matchWidth_() const {
        return matchWidth > 0 ? matchWidth : lineWidth_();
    }
};

// ════════════════════════════════════════════════════════════════════════════
//  SSSR Ring Detection (for layout and double-bond offset)
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<std::vector<int>> findLayoutRings(const MolGraph& g) {
    return computeSSSR(g);
}

// ════════════════════════════════════════════════════════════════════════════
//  2D Layout Engine (ring-first + chain zig-zag + force refinement)
// ════════════════════════════════════════════════════════════════════════════
inline std::vector<Vec2> layout2D(const MolGraph& g, const DepictOptions& opts = {}) {
    const auto coordinates = generateCoords2D(g, opts.bondLength);
    std::vector<Vec2> result;
    result.reserve(coordinates.size());
    for (const auto& point : coordinates) result.push_back({point.x, point.y});
    return result;
}

// ════════════════════════════════════════════════════════════════════════════
//  Labels and drawing (element + implicit H + charge, halo, bonds, map numbers)
// ════════════════════════════════════════════════════════════════════════════
namespace detail_depict {

/** Preserve ordinary font lists without allowing CSS or SVG delimiters. */
inline std::string fontFamily(const std::string& input) {
    std::string family;
    for (unsigned char character : input) {
        if ((character >= 'A' && character <= 'Z') || (character >= 'a' && character <= 'z')
                || (character >= '0' && character <= '9') || character >= 0x80
                || character == ' ' || character == ',' || character == '-' || character == '_'
                || character == '.' || character == '\'' || character == '"') family += character;
    }
    return family.find_first_not_of(" ,.'\"_-") == std::string::npos
        ? "Arial, Helvetica, sans-serif" : family;
}

/**
 * Resolve Kekule bond orders on a local copy of an aromatic graph.  A graph without a Kekule
 * structure is returned unchanged and its aromatic rings are drawn with circles.
 */
inline const MolGraph& bondStyle(const MolGraph& graph, const DepictOptions& options, MolGraph& storage) {
    if (options.aromaticStyle != AromaticStyle::KEKULE
            || std::none_of(graph.aromatic.begin(), graph.aromatic.end(), [](uint8_t value) { return value != 0; })) return graph;
    {
        const std::lock_guard<std::recursive_mutex> cacheLock(graph.lazyCacheMutex());
        storage = graph;
    }
    return storage.kekulize() ? storage : graph;
}

/** Advance width in em of text set in Arial or Helvetica, which share their metrics. */
inline double textAdvance(const std::string& text) {
    static constexpr short upper[26] = {667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833,
                                        722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611};
    static constexpr short lower[26] = {556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833,
                                        556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500};
    double width = 0;
    for (size_t at = 0; at < text.size(); ++at) {
        const unsigned char character = text[at];
        if (character >= 'A' && character <= 'Z') width += upper[character - 'A'];
        else if (character >= 'a' && character <= 'z') width += lower[character - 'a'];
        else if (character >= '0' && character <= '9') width += 556;
        else if (character == 0xE2) { width += 584; at += 2; }  // U+2212 minus sign
        else width += 584;
    }
    return width / 1000;
}

inline std::string chargeText(int charge) {
    const std::string sign = charge > 0 ? "+" : "\xe2\x88\x92";  // minus sign
    return std::abs(charge) > 1 ? std::to_string(std::abs(charge)) + sign : sign;
}

struct TextRun {
    std::string text;
    double x, y;   // left end and vertical centre, relative to the atom
    double size;
};

struct Box { double left, top, right, bottom; };

inline bool sameColour(const Color& first, const Color& second) {
    return first.r == second.r && first.g == second.g && first.b == second.b && first.a == second.a;
}

/** Box of a text run's capitals and digits, relative to its atom. */
inline Box runBox(const TextRun& run) {
    return {run.x, run.y - run.size * 0.36, run.x + textAdvance(run.text) * run.size, run.y + run.size * 0.36};
}

/** Fraction of the way from a to b at which the segment leaves the box, or -1 when it misses it. */
inline double leaveBox(Vec2 a, Vec2 b, const Box& box) {
    const Vec2 d = b - a;
    const double p[4] = {-d.x, d.x, -d.y, d.y};
    const double q[4] = {a.x - box.left, box.right - a.x, a.y - box.top, box.bottom - a.y};
    double enter = 0, leave = 1;
    for (int side = 0; side < 4; ++side) {
        if (std::abs(p[side]) < 1e-12) {
            if (q[side] < 0) return -1;
            continue;
        }
        const double t = q[side] / p[side];
        if (p[side] < 0) enter = std::max(enter, t);
        else leave = std::min(leave, t);
        if (enter > leave) return -1;
    }
    return leave;
}

struct AtomLabel {
    std::string symbol;         // "N", "O", "S", etc.
    int         hCount = 0;     // implicit H to display
    int         charge = 0;     // formal charge
    bool        visible = false;
    std::vector<TextRun> runs;  // the symbol first, then hydrogens and charge
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

    // An isotope's mass number, which also makes a carbon visible
    const int mass = idx < (int)g.massNumber.size() ? g.massNumber[idx] : 0;

    // Visibility
    const bool plainCarbon = Z == 6 && mass <= 0;
    switch (opts.hMode) {
        case HydrogenMode::ALL:
            lbl.visible = true;
            break;
        case HydrogenMode::NONE:
            lbl.visible = !plainCarbon || opts.showCarbonLabels || lbl.charge != 0;
            lbl.hCount = 0;
            break;
        case HydrogenMode::TERMINAL:
            lbl.visible = !plainCarbon || opts.showCarbonLabels
                || (int)g.neighbors[idx].size() <= 1 || lbl.charge != 0;
            break;
        case HydrogenMode::HETERO_ONLY:
        default:
            lbl.visible = !plainCarbon || opts.showCarbonLabels || lbl.charge != 0;
            break;
    }

    if (!lbl.visible) return lbl;

    // The symbol is centred on the atom, so bonds meet the element rather than the whole label.
    const double size = opts.fontSize_(), subscript = size * opts.subscriptScale;
    const double half = textAdvance(lbl.symbol) * size * 0.5;
    lbl.runs.push_back({lbl.symbol, -half, 0, size});
    // The mass number of an isotope is a superscript before the symbol, as in 13C.
    const std::string massText = mass > 0 ? std::to_string(mass) : "";
    const double massWidth = textAdvance(massText) * subscript;
    if (mass > 0) lbl.runs.push_back({massText, -half - massWidth, -size * 0.4, subscript});
    double chargeAt = half;
    if (lbl.hCount > 0) {
        // Hydrogens take the first side free of bonds and nearby atoms: right, left, then below or above.
        std::vector<Vec2> nearby;
        for (int other = 0; other < g.n; ++other)
            if (other != idx && (g.bondOrder(idx, other) != 0 || len(pos[other] - pos[idx]) < opts.bondLength * 1.3))
                nearby.push_back(norm(pos[other] - pos[idx]));
        auto crowding = [&](Vec2 side) {
            double most = -1;
            for (Vec2 direction : nearby) most = std::max(most, dot(direction, side));
            return most;
        };
        const Vec2 side = crowding({1, 0}) < 0.55 ? Vec2{1, 0}
                        : crowding({-1, 0}) < 0.55 ? Vec2{-1, 0}
                        : crowding({0, 1}) <= crowding({0, -1}) ? Vec2{0, 1} : Vec2{0, -1};
        const std::string count = lbl.hCount > 1 ? std::to_string(lbl.hCount) : "";
        const double hydrogen = textAdvance("H") * size, digits = textAdvance(count) * subscript;
        const double x = side.x > 0 ? half : side.x < 0 ? -half - massWidth - hydrogen - digits : -hydrogen * 0.5;
        const double y = side.y * size * 0.9;
        lbl.runs.push_back({"H", x, y, size});
        if (!count.empty()) lbl.runs.push_back({count, x + hydrogen, y + size * 0.25, subscript});
        if (side.x > 0) chargeAt = half + hydrogen + digits;
    }
    if (lbl.charge != 0) lbl.runs.push_back({chargeText(lbl.charge), chargeAt, -size * 0.4, subscript});
    return lbl;
}

/** Where a bond from an atom towards another leaves the atom's label. */
inline Vec2 clipToLabel(const AtomLabel& label, Vec2 atom, Vec2 other, double gap) {
    if (!label.visible) return atom;
    double leave = 0;
    for (const auto& run : label.runs) {
        const Box box = runBox(run);
        leave = std::max(leave, leaveBox(atom, other, {atom.x + box.left - gap, atom.y + box.top - gap,
                                                       atom.x + box.right + gap, atom.y + box.bottom + gap}));
    }
    return atom + (other - atom) * std::min(leave, 0.45);
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

struct Bounds {
    double minX = 1e9, minY = 1e9, maxX = -1e9, maxY = -1e9;

    void add(double left, double top, double right, double bottom) {
        minX = std::min(minX, left); minY = std::min(minY, top);
        maxX = std::max(maxX, right); maxY = std::max(maxY, bottom);
    }
    /** An atom with its label and, when it carries a halo, the halo's radius. */
    void addAtom(Vec2 at, const AtomLabel& label, double halo) {
        add(at.x - halo, at.y - halo, at.x + halo, at.y + halo);
        if (label.visible) for (const auto& run : label.runs) {
            const Box box = runBox(run);
            add(at.x + box.left, at.y + box.top, at.x + box.right, at.y + box.bottom);
        }
    }
};

/** Canvas, label styles and background shared by single and paired drawings. */
inline void openSVG(std::ostringstream& svg, int width, int height, double viewWidth, double viewHeight,
                    const DepictOptions& opts) {
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height << "\"";
    svg << " viewBox=\"0 0 " << viewWidth << " " << viewHeight << "\" preserveAspectRatio=\"xMidYMid meet\"";
    svg << " shape-rendering=\"geometricPrecision\" text-rendering=\"geometricPrecision\">\n";
    svg << "<defs><style>\n";
    svg << "  text { font-family: " << fontFamily(opts.fontFamily) << "; }\n";
    svg << "  .sub  { dominant-baseline: central; }\n";
    svg << "  .map  { dominant-baseline: central; text-anchor: middle; }\n";
    svg << "  .idx  { dominant-baseline: central; text-anchor: start; }\n";
    svg << "</style></defs>\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"" << rgb(opts.backgroundColor) << "\"/>\n";
}

/** Relative luminance of a colour (WCAG 2). */
inline double luminance(const Color& colour) {
    auto channel = [](int value) {
        const double level = value / 255.0;
        return level <= 0.04045 ? level / 12.92 : std::pow((level + 0.055) / 1.055, 2.4);
    };
    return 0.2126 * channel(colour.r) + 0.7152 * channel(colour.g) + 0.0722 * channel(colour.b);
}

/** Contrast ratio of two colours (WCAG 2), from 1 to 21. */
inline double contrast(const Color& first, const Color& second) {
    const double a = luminance(first) + 0.05, b = luminance(second) + 0.05;
    return a > b ? a / b : b / a;
}

/**
 * Opaque colour of a hue laid over the background at a strength: up to 1 the hue is mixed with
 * the background, and from 1 to 2 it deepens towards black.
 */
inline Color shade(const Color& colour, double strength, const Color& background) {
    auto mix = [&](int over, int under) {
        return (int)std::lround(strength <= 1 ? under + (over - under) * strength : over * (2 - strength));
    };
    return {mix(colour.r, background.r), mix(colour.g, background.g), mix(colour.b, background.b), 1.0};
}

/**
 * Flat tint of a highlight colour and the hairline that edges it: the same hue, deepened until it
 * stands 1.35:1 from the tint, which for the default gold is close to the undiluted MATCH_FILL gold.
 * An opaque colour takes the strength of the match fill, and a colour of strength 0 draws nothing.
 */
struct HaloTone { Color tint, edge; bool visible = true; };
inline HaloTone haloTone(const Color& colour, const DepictOptions& opts) {
    const double strength = std::clamp(colour.a < 1.0 ? colour.a : opts.matchFill.a, 0.0, 1.0);
    const Color tint = shade(colour, strength, opts.backgroundColor);
    if (strength == 0) return {tint, tint, false};
    for (double deeper = strength; deeper < 2; deeper += 0.01) {
        const Color edge = shade(colour, deeper, opts.backgroundColor);
        if (contrast(edge, tint) >= 1.35) return {tint, edge};
    }
    return {tint, shade(colour, 2, opts.backgroundColor)};
}

inline double segmentDistance(Vec2 point, Vec2 from, Vec2 to) {
    const Vec2 along = to - from;
    const double lengthSquared = dot(along, along);
    const double t = lengthSquared > 1e-12 ? std::clamp(dot(point - from, along) / lengthSquared, 0.0, 1.0) : 0.0;
    return len(point - (from + along * t));
}

/** Marks already on a canvas, which label knockouts and map numbers are fitted round. */
struct Scene {
    struct Stroke { Vec2 from, to; double halfWidth; };
    struct Ring { Vec2 centre; double inner; bool circle; };
    struct Tint { Vec2 from, to; double radius; Color colour; };  // a band, or a disc when from == to
    std::vector<Stroke> strokes;   // bond marks and the pair arrow
    std::vector<Box> labels;       // label runs
    std::vector<Ring> rings;
    std::vector<Vec2> atoms;
    std::vector<Tint> halo;        // in drawing order
    double haloEdge = 0;
    Box view{0, 0, 0, 0};          // the view box, which map numbers keep inside

    /** The halo tint at a point, or null off the halo. */
    const Color* tintAt(Vec2 point) const {
        for (auto shape = halo.rbegin(); shape != halo.rend(); ++shape)
            if (segmentDistance(point, shape->from, shape->to) <= shape->radius) return &shape->colour;
        return nullptr;
    }
};

/**
 * Draw one molecule at canvas positions: highlight halo, bonds, ring circles and labels, adding
 * its marks to the scene.  `bonds` are the highlighted bonds; the per-atom and per-bond highlight
 * colours of the options apply only when `colourOverrides` is set.  `pixel` is the size of one
 * canvas pixel in drawing units.
 */
inline void drawMolecule(std::ostringstream& svg, Scene& scene, const MolGraph& g, const std::vector<Vec2>& at,
                         const std::vector<AtomLabel>& labels, const std::set<int>& atoms,
                         const std::set<std::pair<int,int>>& bonds, const DepictOptions& opts,
                         bool colourOverrides, double pixel) {
    const double LW  = opts.lineWidth_();
    const double BW  = opts.boldWidth_();
    const double BS  = opts.bondSpacing_();
    const double FS  = opts.fontSize_();
    const double HR  = opts.highlightRadius_();
    const double MW  = opts.matchWidth_();
    const double MNS = opts.mapNumSize_();
    const double EW  = std::max(LW * 0.6, pixel);  // halo edge, never thinner than a pixel
    const std::string background = rgb(opts.backgroundColor);
    const auto rings = findLayoutRings(g);
    auto matched = [&](int first, int second) {
        return bonds.count({first, second}) > 0 || bonds.count({second, first}) > 0;
    };
    scene.haloEdge = EW;
    scene.atoms.insert(scene.atoms.end(), at.begin(), at.end());
    for (const auto& ring : rings) {
        const Vec2 centre = ringCentroid(ring, at);
        scene.rings.push_back({centre, len(at[ring[0]] - centre) * std::cos(PI / ring.size()),
                               std::all_of(ring.begin(), ring.end(), [&](int atom) { return g.aromatic[atom] != 0; })});
    }

    // Halo: flat tints laid over the background, so overlaps keep one even tone, each shape edged
    // by a hairline in a deeper shade of its hue.  Every edge goes beneath every tint.
    std::ostringstream edges, fills;
    edges << std::fixed << std::setprecision(2);
    fills << std::fixed << std::setprecision(2);
    std::vector<std::pair<Color, HaloTone>> tones;
    auto tone = [&](const Color& colour) -> HaloTone {
        for (const auto& [key, value] : tones)
            if (sameColour(key, colour)) return value;
        tones.emplace_back(colour, haloTone(colour, opts));
        return tones.back().second;
    };
    auto bondColour = [&](int first, int second) {
        auto colourOverride = opts.bondHighlightColors.find({first, second});
        if (colourOverride == opts.bondHighlightColors.end())
            colourOverride = opts.bondHighlightColors.find({second, first});
        return colourOverrides && colourOverride != opts.bondHighlightColors.end() ? colourOverride->second : opts.matchFill;
    };
    // A small ring whose bonds all carry one halo colour is filled, so the halo stays one shape.
    for (const auto& ring : rings) {
        if (ring.size() > 8) continue;
        const Vec2 centre = ringCentroid(ring, at);
        std::vector<int> cycle = ring;
        std::sort(cycle.begin(), cycle.end(), [&](int first, int second) {
            return std::atan2(at[first].y - centre.y, at[first].x - centre.x)
                 < std::atan2(at[second].y - centre.y, at[second].x - centre.x);
        });
        const Color colour = bondColour(cycle[0], cycle[1]);
        bool filled = true;
        for (size_t k = 0; k < cycle.size() && filled; ++k) {
            const int first = cycle[k], second = cycle[(k + 1) % cycle.size()];
            filled = matched(first, second) && sameColour(bondColour(first, second), colour);
        }
        if (!filled || !tone(colour).visible) continue;
        const Color tint = tone(colour).tint;
        fills << "<polygon class=\"ring-highlight\" points=\"";
        for (int atom : cycle) fills << at[atom].x << "," << at[atom].y << " ";
        fills << "\" fill=\"" << rgb(tint) << "\"/>\n";
        scene.halo.push_back({centre, centre, len(at[ring[0]] - centre) * std::cos(PI / ring.size()), tint});
    }
    for (int first = 0; first < g.n; ++first) {
        for (int second : g.neighbors[first]) {
            if (second <= first || !matched(first, second)) continue;
            const HaloTone halo = tone(bondColour(first, second));
            if (!halo.visible) continue;
            edges << "<path d=\"M " << at[first].x << " " << at[first].y << " L " << at[second].x << " " << at[second].y
                << "\" fill=\"none\" stroke=\"" << rgb(halo.edge) << "\" stroke-width=\"" << 2 * (HR + EW)
                << "\" stroke-linecap=\"round\"/>\n";
            fills << "<path class=\"bond-highlight\" d=\"M " << at[first].x << " " << at[first].y
                << " L " << at[second].x << " " << at[second].y << "\" fill=\"none\" stroke=\"" << rgb(halo.tint)
                << "\" stroke-width=\"" << 2 * HR << "\" stroke-linecap=\"round\"/>\n";
            scene.halo.push_back({at[first], at[second], HR, halo.tint});
        }
    }
    for (int atom : atoms) {
        if (atom < 0 || atom >= g.n) continue;
        Color colour = opts.matchFill;
        const auto colourOverride = opts.atomHighlightColors.find(atom);
        if (colourOverrides && colourOverride != opts.atomHighlightColors.end()) colour = colourOverride->second;
        const HaloTone halo = tone(colour);
        if (!halo.visible) continue;
        edges << "<circle cx=\"" << at[atom].x << "\" cy=\"" << at[atom].y << "\" r=\"" << HR + EW
            << "\" fill=\"" << rgb(halo.edge) << "\"/>\n";
        fills << "<circle cx=\"" << at[atom].x << "\" cy=\"" << at[atom].y << "\" r=\"" << HR
            << "\" fill=\"" << rgb(halo.tint) << "\"/>\n";
        scene.halo.push_back({at[atom], at[atom], HR, halo.tint});
    }
    if (fills.tellp() > 0) svg << "<g class=\"highlight\">\n" << edges.str() << fills.str() << "</g>\n";

    // Bonds: matched bonds take matchStroke and matchWidth, which default to the default bond
    // colour and the line width, so a custom bondColor needs a matching matchStroke.
    std::unordered_map<int64_t, std::pair<int, int>> stereoWedges;
    if (std::any_of(g.tetraChirality.begin(), g.tetraChirality.end(), [](int winding) { return winding != 0; })) {
        std::vector<Point2D> coordinates;
        coordinates.reserve(at.size());
        // SVG's vertical axis points down; the winding calculation uses y up.
        for (const auto& point : at) coordinates.push_back({point.x, -point.y});
        stereoWedges = mol_detail::tetrahedralWedges(g, coordinates);
    }
    auto line = [&](Vec2 from, Vec2 to, const Color& colour, double width) {
        svg << "<line x1=\"" << from.x << "\" y1=\"" << from.y << "\" x2=\"" << to.x << "\" y2=\"" << to.y
            << "\" stroke=\"" << rgb(colour) << "\" stroke-width=\"" << width << "\" stroke-linecap=\"round\"/>\n";
        scene.strokes.push_back({from, to, width * 0.5});
    };
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue;
            const Vec2 p1 = at[i], p2 = at[j];
            Vec2 cp1 = clipToLabel(labels[i], p1, p2, FS * 0.08);
            Vec2 cp2 = clipToLabel(labels[j], p2, p1, FS * 0.08);
            const Color colour = matched(i, j) ? opts.matchStroke : opts.bondColor;
            const double width = matched(i, j) ? MW : LW;
            const int order = g.bondOrder(i, j);
            const bool arom = g.bondAromatic(i, j);

            BondStereoType stereo = BondStereoType::NONE;
            const auto wedge = stereoWedges.find(MolGraph::bondKey(i, j));
            if (wedge != stereoWedges.end()) {
                stereo = wedge->second.second == 1 ? BondStereoType::WEDGE_UP : BondStereoType::WEDGE_DOWN;
                if (wedge->second.first == j) std::swap(cp1, cp2);
            }
            const Vec2 bondDir = norm(cp2 - cp1);
            const Vec2 bondPerp = perp(bondDir);

            if (stereo == BondStereoType::WEDGE_UP) {
                // Solid filled wedge: narrow at stereo centre, wide at substituent
                const Vec2 w = bondPerp * (BW * 0.5);
                svg << "<polygon points=\""
                    << cp1.x << "," << cp1.y << " "
                    << (cp2.x + w.x) << "," << (cp2.y + w.y) << " "
                    << (cp2.x - w.x) << "," << (cp2.y - w.y)
                    << "\" fill=\"" << rgb(colour) << "\" stroke=\"none\"/>\n";
                scene.strokes.push_back({cp1, cp2, BW * 0.5});
            } else if (stereo == BondStereoType::WEDGE_DOWN) {
                // Dashed wedge: series of stripes widening from stereo centre
                const int nStripes = std::max(4, (int)(len(cp2 - cp1) / opts.hashSpacing_()));
                for (int s = 0; s < nStripes; s++) {
                    const double t = (s + 0.5) / nStripes;
                    const Vec2 mid = cp1 + (cp2 - cp1) * t;
                    const Vec2 w = bondPerp * (BW * 0.5 * t);
                    svg << "<line x1=\"" << (mid.x - w.x) << "\" y1=\"" << (mid.y - w.y)
                        << "\" x2=\"" << (mid.x + w.x) << "\" y2=\"" << (mid.y + w.y)
                        << "\" stroke=\"" << rgb(colour) << "\" stroke-width=\"" << LW
                        << "\" stroke-linecap=\"round\"/>\n";
                }
                scene.strokes.push_back({cp1, cp2, BW * 0.5});
            } else if (order == 2 && !arom) {
                Vec2 ringC;
                if (bondInRing(i, j, rings, at, ringC)) {
                    // Full line on the bond axis and a shorter one towards the ring centre
                    line(cp1, cp2, colour, width);
                    const Vec2 off = bondPerp * (dot(bondPerp, ringC - (p1 + p2) * 0.5) < 0 ? -BS : BS);
                    const double trim = len(cp2 - cp1) * 0.15;
                    line(cp1 + bondDir * trim + off, cp2 - bondDir * trim + off, colour, width);
                } else {
                    const Vec2 off = bondPerp * (BS * 0.5);
                    line(cp1 + off, cp2 + off, colour, width);
                    line(cp1 - off, cp2 - off, colour, width);
                }
            } else if (order == 3 && !arom) {
                const Vec2 off = bondPerp * BS;
                line(cp1, cp2, colour, width);
                line(cp1 + off, cp2 + off, colour, width * 0.7);
                line(cp1 - off, cp2 - off, colour, width * 0.7);
            } else {
                line(cp1, cp2, colour, width);
            }
        }
    }

    // Circles for aromatic rings, when asked for or when the molecule has no Kekule structure
    for (const auto& ring : rings) {
        if (!std::all_of(ring.begin(), ring.end(), [&](int atom) { return g.aromatic[atom] != 0; })) continue;
        const Vec2 centre = ringCentroid(ring, at);
        double radius = 0;
        for (int atom : ring) radius += len(at[atom] - centre);
        radius = radius / ring.size() * 0.55;
        svg << "<circle cx=\"" << centre.x << "\" cy=\"" << centre.y << "\" r=\"" << radius
            << "\" fill=\"none\" stroke=\"" << rgb(opts.bondColor)
            << "\" stroke-width=\"" << LW * 0.7 << "\"/>\n";
    }

    // Labels on the halo sit on a knockout in the colour beneath them: the tint where the ink keeps
    // 4.5:1 on it, otherwise the background.
    auto ink = [&](int atom) {
        const auto colourOverride = opts.atomColorOverrides.find(g.atomicNum[atom]);
        return colourOverride != opts.atomColorOverrides.end() ? colourOverride->second : atomColor(g.atomicNum[atom]);
    };
    auto writeRun = [&](const TextRun& run, Vec2 point, const std::string& fill, bool outline) {
        svg << "<text x=\"" << point.x + run.x << "\" y=\"" << point.y + run.y
            << "\" class=\"sub\" font-size=\"" << run.size << "\" fill=\"" << fill << "\"";
        if (outline) svg << " stroke=\"" << fill << "\"";
        svg << ">" << run.text << "</text>\n";
    };
    bool knockout = false;
    for (int atom = 0; atom < g.n; ++atom) {
        if (!labels[atom].visible) continue;
        std::vector<const Color*> under;
        for (const auto& run : labels[atom].runs) {
            const Box box = runBox(run);
            under.push_back(scene.tintAt(at[atom] + Vec2{(box.left + box.right) * 0.5, (box.top + box.bottom) * 0.5}));
            scene.labels.push_back({at[atom].x + box.left, at[atom].y + box.top, at[atom].x + box.right, at[atom].y + box.bottom});
        }
        if (std::none_of(under.begin(), under.end(), [](const Color* tint) { return tint != nullptr; })) continue;
        if (!knockout) svg << "<g class=\"knockout\" stroke-width=\"" << FS * 0.16 << "\" stroke-linejoin=\"round\">\n";
        knockout = true;
        for (size_t k = 0; k < under.size(); ++k)
            writeRun(labels[atom].runs[k], at[atom],
                     under[k] && contrast(*under[k], ink(atom)) >= 4.5 ? rgb(*under[k]) : background, true);
    }
    if (knockout) svg << "</g>\n";
    for (int atom = 0; atom < g.n; ++atom)
        if (labels[atom].visible)
            for (const auto& run : labels[atom].runs) writeRun(run, at[atom], rgb(ink(atom)), false);

    if (opts.showAtomIndices) {
        for (int atom = 0; atom < g.n; ++atom) {
            double right = 4;
            if (labels[atom].visible)
                for (const auto& run : labels[atom].runs) right = std::max(right, runBox(run).right + 2);
            svg << "<text x=\"" << (at[atom].x + right) << "\" y=\"" << (at[atom].y + FS * 0.45)
                << "\" class=\"idx\" font-size=\"" << MNS
                << "\" fill=\"rgba(120,120,120,0.6)\">" << atom << "</text>\n";
        }
    }
}

struct MapNumber {
    Vec2 atom;
    double start;      // middle of the widest gap between the atom's bonds
    std::string text;
};

inline void addMapNumbers(std::vector<MapNumber>& numbers, const MolGraph& g, const std::vector<Vec2>& at,
                          const std::map<int,int>& mapNumbers) {
    for (const auto& [atom, number] : mapNumbers) {
        if (atom < 0 || atom >= g.n) continue;
        std::vector<double> angles;
        for (int nb : g.neighbors[atom]) angles.push_back(std::atan2(at[nb].y - at[atom].y, at[nb].x - at[atom].x));
        std::sort(angles.begin(), angles.end());
        double start = -PI / 4, widest = 0;
        for (size_t k = 0; k < angles.size(); ++k) {
            const double next = k + 1 < angles.size() ? angles[k + 1] : angles[0] + 2 * PI;
            if (next - angles[k] > widest) { widest = next - angles[k]; start = angles[k] + widest * 0.5; }
        }
        numbers.push_back({at[atom], start, std::to_string(number)});
    }
}

/**
 * Draw map numbers, each at the nearest place round its atom that touches no bond, label or
 * other number.  In order of preference a place is clear of everything else; inside a ring drawn
 * without a circle; level beside a label, where it reads as part of the formula; nearer another
 * atom than its own.  Places off the halo are preferred, and the most crowded numbers choose first.
 * Every place lies inside the scene's view box.  When no place is free, the number goes in the
 * widest gap between its atom's bonds, 0.6 bond lengths out, moved inside the view box if need be.
 * `pixel` is the size of one canvas pixel in drawing units.
 */
inline void drawMapNumbers(std::ostringstream& svg, const Scene& scene, const std::vector<MapNumber>& numbers,
                           const DepictOptions& opts, double pixel) {
    const double FS = opts.fontSize_(), BL = opts.bondLength;
    // An automatic size stays legible on a small canvas, at nine pixels or more, but never
    // outgrows the atom labels.
    const double MNS = opts.mapNumSize > 0 ? opts.mapNumSize : std::min(FS, std::max(opts.mapNumSize_(), 9 * pixel));
    const double clearance = FS * 0.12, spacing = FS * 0.5;
    const std::string background = rgb(opts.backgroundColor);
    std::vector<Box> placed;
    auto overlap = [](const Box& first, const Box& second, double gap) {
        return first.left < second.right + gap && second.left < first.right + gap
            && first.top < second.bottom + gap && second.top < first.bottom + gap;
    };
    // Places lie on orbits round the atom, from 0.3 to 1.3 bond lengths out.
    constexpr int kOrbits = 20;
    auto orbitRadius = [&](int orbit) { return BL * (0.3 + orbit * 0.05); };
    // The marks within reach of one number's places; the rest cannot touch them.
    const double farthest = orbitRadius(kOrbits);
    struct Nearby {
        std::vector<const Box*> labels, numbers;
        std::vector<const Scene::Stroke*> strokes;
        std::vector<const Scene::Ring*> rings;
        std::vector<const Scene::Tint*> halo;
        std::vector<Vec2> atoms;
    };
    auto nearby = [&](const MapNumber& number, Vec2 half) {
        const double reach = farthest + std::max(half.x, half.y) + FS;
        auto boxDistance = [&](const Box& box) {
            return std::hypot(std::max({box.left - number.atom.x, 0.0, number.atom.x - box.right}),
                              std::max({box.top - number.atom.y, 0.0, number.atom.y - box.bottom}));
        };
        Nearby marks;
        for (const auto& label : scene.labels)
            if (boxDistance(label) < reach) marks.labels.push_back(&label);
        for (const auto& other : placed)
            if (boxDistance(other) < reach + FS) marks.numbers.push_back(&other);
        for (const auto& stroke : scene.strokes)
            if (segmentDistance(number.atom, stroke.from, stroke.to) < reach + stroke.halfWidth)
                marks.strokes.push_back(&stroke);
        for (const auto& ring : scene.rings)
            if (len(ring.centre - number.atom) < reach + ring.inner) marks.rings.push_back(&ring);
        for (const auto& shape : scene.halo)
            if (segmentDistance(number.atom, shape.from, shape.to) < reach + shape.radius + scene.haloEdge)
                marks.halo.push_back(&shape);
        // Another atom matters when it is nearer a place than the number's own atom.
        for (Vec2 atom : scene.atoms)
            if (len(atom - number.atom) <= 2 * farthest + clearance) marks.atoms.push_back(atom);
        return marks;
    };
    // Rank of a place in the order above, from 0, or -1 when it touches a mark.  Beside a label,
    // the rank grows with the overlap, and more on the right, where subscripts and charges go,
    // so a number forced there keeps to the corners.
    // Half a pixel in from the view box, so that rounding the printed positions keeps numbers inside.
    const Box view{scene.view.left + pixel * 0.5, scene.view.top + pixel * 0.5,
                   scene.view.right - pixel * 0.5, scene.view.bottom - pixel * 0.5};
    const bool bounded = view.right > view.left && view.bottom > view.top;
    auto rank = [&](const MapNumber& number, const Nearby& marks, const Box& box) -> double {
        if (bounded && (box.left < view.left || box.top < view.top || box.right > view.right || box.bottom > view.bottom))
            return -1;
        for (const Box* label : marks.labels)
            if (overlap(box, *label, 2 * clearance)) return -1;
        for (const Box* other : marks.numbers) {
            // Two numbers level with each other need a wider gap, or they read as one.
            const double across = std::max(other->left - box.right, box.left - other->right);
            const double down = std::max(other->top - box.bottom, box.top - other->bottom);
            if ((across < spacing && down < spacing) || (down < 0 && across < FS * 0.75)) return -1;
        }
        for (const Scene::Stroke* stroke : marks.strokes) {
            const double grow = stroke->halfWidth + clearance;
            if (leaveBox(stroke->from, stroke->to, {box.left - grow, box.top - grow, box.right + grow, box.bottom + grow}) >= 0)
                return -1;
        }
        const Vec2 centre = {(box.left + box.right) * 0.5, (box.top + box.bottom) * 0.5};
        double level = 0;
        for (const Scene::Ring* ring : marks.rings)
            if (len(centre - ring->centre) < ring->inner) {
                if (ring->circle) return -1;
                level = 1;
            }
        const double own = len(centre - number.atom);
        for (Vec2 atom : marks.atoms)
            if (len(centre - atom) <= own && len(atom - number.atom) > 1e-6) return 3;
        double beside = 0;
        for (const Box* labelBox : marks.labels) {
            const Box& label = *labelBox;
            const bool side = centre.x > label.right ? box.left < label.right + FS * 0.8
                            : centre.x < label.left && box.right > label.left - FS * 0.8;
            const double depth = std::min(box.bottom, label.bottom + FS * 0.15) - std::max(box.top, label.top - FS * 0.15);
            if (side && depth > 0)
                beside = std::max(beside, std::min(depth / (2 * FS) + (centre.x > label.right ? 0.4 : 0), 0.9));
        }
        return beside > 0 ? 2 + beside : level;
    };
    auto halfSize = [&](const MapNumber& number) { return Vec2{textAdvance(number.text) * MNS * 0.5, MNS * 0.36}; };
    // The best place for a number, given those placed so far, and its rank (4 when none is free).
    auto place = [&](const MapNumber& number, int* clear, double* bestLevel) {
        const Vec2 half = halfSize(number);
        const Nearby marks = nearby(number, half);
        auto reachesHalo = [&](Vec2 centre) {
            return std::any_of(marks.halo.begin(), marks.halo.end(), [&](const Scene::Tint* shape) {
                return segmentDistance(centre, shape->from, shape->to)
                     < shape->radius + scene.haloEdge + std::max(half.x, half.y);
            });
        };
        Vec2 best = number.atom + Vec2{std::cos(number.start), std::sin(number.start)} * (BL * 0.6);
        if (bounded) {
            best.x = std::clamp(best.x, view.left + half.x, std::max(view.left + half.x, view.right - half.x));
            best.y = std::clamp(best.y, view.top + half.y, std::max(view.top + half.y, view.bottom - half.y));
        }
        double bestScore = 1e18;
        *bestLevel = 4;
        // Places lie in 48 directions, 7.5 degrees apart, turning either way from the widest gap.
        for (int step = -23; step <= 24; ++step) {
            const double angle = number.start + step * PI / 24;
            const Vec2 direction = {std::cos(angle), std::sin(angle)};
            for (int orbit = 0; orbit <= kOrbits; ++orbit) {
                const double radius = orbitRadius(orbit);
                const Vec2 centre = number.atom + direction * radius;
                const double level = rank(number, marks,
                                          {centre.x - half.x, centre.y - half.y, centre.x + half.x, centre.y + half.y});
                if (level < 0) continue;
                if (clear && level == 0) ++*clear;
                // Further out than 0.6 bond lengths a number starts to read as a bond's, so the
                // distance counts four times over; the halo costs less, as the plate keeps it legible.
                const double score = level * BL * 10 + radius + std::max(0.0, radius - BL * 0.6) * 3
                                   + std::abs(step) * BL * 0.015
                                   + (reachesHalo(centre) ? BL * 0.15 : 0);
                if (score < bestScore) {
                    best = centre;
                    bestScore = score;
                    *bestLevel = level;
                }
            }
        }
        return best;
    };

    std::vector<int> clear(numbers.size(), 0), order(numbers.size());
    std::vector<double> levels(numbers.size());
    for (size_t k = 0; k < numbers.size(); ++k) {
        place(numbers[k], &clear[k], &levels[k]);
        order[k] = (int)k;
    }
    std::stable_sort(order.begin(), order.end(), [&](int first, int second) { return clear[first] < clear[second]; });
    std::vector<Vec2> centres(numbers.size());
    // Places the numbers in `order`; returns how many were forced beside a label or nearer another atom.
    auto layOut = [&]() {
        placed.clear();
        int forced = 0;
        for (int k : order) {
            centres[k] = place(numbers[k], nullptr, &levels[k]);
            forced += levels[k] >= 2;
            const Vec2 half = halfSize(numbers[k]);
            placed.push_back({centres[k].x - half.x, centres[k].y - half.y, centres[k].x + half.x, centres[k].y + half.y});
        }
        return forced;
    };
    // Numbers that were forced aside choose first in the next round, which often frees a place
    // for them; the round that forces the fewest is kept.
    int fewest = layOut();
    auto bestCentres = centres;
    for (int round = 0, forced = fewest; forced > 0 && round < 4; ++round) {
        std::stable_partition(order.begin(), order.end(), [&](int k) { return levels[k] >= 2; });
        forced = layOut();
        if (forced < fewest) {
            fewest = forced;
            bestCentres = centres;
        }
    }
    centres = bestCentres;
    // Each number sits on a plate in the colour beneath it, as labels do.
    for (size_t k = 0; k < numbers.size(); ++k) {
        const Color* tint = scene.tintAt(centres[k]);
        const std::string plate = tint && contrast(*tint, palette::MAPNUM) >= 4.5 ? rgb(*tint) : background;
        svg << "<text x=\"" << centres[k].x << "\" y=\"" << centres[k].y << "\" class=\"sub\" text-anchor=\"middle\" font-size=\""
            << MNS << "\" fill=\"" << plate << "\" stroke=\"" << plate << "\" stroke-width=\"" << MNS * 0.3
            << "\" stroke-linejoin=\"round\">" << numbers[k].text << "</text>\n";
        svg << "<text x=\"" << centres[k].x << "\" y=\"" << centres[k].y << "\" class=\"map\" font-size=\"" << MNS
            << "\" fill=\"" << rgb(palette::MAPNUM) << "\">" << numbers[k].text << "</text>\n";
    }
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
    const MolGraph& input,
    const std::vector<Vec2>& pos,
    const std::map<int,int>& mapping = {},          // atom index -> map number
    const std::set<int>& highlightAtoms = {},
    const std::set<std::pair<int,int>>& highlightBonds = {},
    const DepictOptions& opts = {}
) {
    MolGraph styled;
    const MolGraph& g = detail_depict::bondStyle(input, opts, styled);
    if (g.n == 0) return {"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1\" height=\"1\"/>", 1, 1};

    // Build atom labels
    std::vector<detail_depict::AtomLabel> labels(g.n);
    for (int i = 0; i < g.n; i++)
        labels[i] = detail_depict::buildLabel(g, i, opts, pos);

    // Without an explicit bond list, the bonds joining two highlighted atoms are highlighted.
    std::set<std::pair<int,int>> bonds = highlightBonds;
    if (bonds.empty())
        for (int first : highlightAtoms)
            if (first >= 0 && first < g.n)
                for (int second : g.neighbors[first])
                    if (first < second && highlightAtoms.count(second)) bonds.insert({first, second});
    std::set<int> haloAtoms = highlightAtoms;
    for (const auto& [first, second] : bonds) { haloAtoms.insert(first); haloAtoms.insert(second); }

    // Bounding box
    detail_depict::Bounds bounds;
    for (int i = 0; i < g.n; i++)
        bounds.addAtom(pos[i], labels[i], haloAtoms.count(i) ? opts.highlightRadius_() : 0);
    double minX = bounds.minX, minY = bounds.minY, maxX = bounds.maxX, maxY = bounds.maxY;
    double pad = opts.padding;
    int W = opts.width  > 0 ? opts.width  : std::max(100, (int)(maxX - minX + 2 * pad));
    int H = opts.height > 0 ? opts.height : std::max(100, (int)(maxY - minY + 2 * pad));

    // The view box keeps the canvas's aspect ratio and grows uniformly until the drawing fits.
    const double scale = std::max({1.0, (maxX - minX + 2 * pad) / W, (maxY - minY + 2 * pad) / H});
    const double viewWidth = W * scale, viewHeight = H * scale;
    const double offsetX = (viewWidth - (maxX - minX)) * 0.5 - minX;
    const double offsetY = (viewHeight - (maxY - minY)) * 0.5 - minY;
    std::vector<Vec2> at(g.n);
    for (int i = 0; i < g.n; i++) at[i] = {pos[i].x + offsetX, pos[i].y + offsetY};

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);
    detail_depict::openSVG(svg, W, H, viewWidth, viewHeight, opts);
    detail_depict::Scene scene;
    scene.view = {0, 0, viewWidth, viewHeight};
    detail_depict::drawMolecule(svg, scene, g, at, labels, highlightAtoms, bonds, opts, true, scale);
    if (opts.showMapNumbers) {
        std::vector<detail_depict::MapNumber> numbers;
        detail_depict::addMapNumbers(numbers, g, at, mapping);
        detail_depict::drawMapNumbers(svg, scene, numbers, opts, scale);
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
    const MolGraph& input1, const MolGraph& input2,
    const std::map<int,int>& mapping,  // g1 atom → g2 atom
    const DepictOptions& opts = {}
) {
    MolGraph styled1, styled2;
    const MolGraph& g1 = detail_depict::bondStyle(input1, opts, styled1);
    const MolGraph& g2 = detail_depict::bondStyle(input2, opts, styled2);
    auto pos1 = layout2D(g1, opts);
    auto pos2 = layout2D(g2, opts);

    // Matched atoms carry the same map number in both molecules.
    std::set<int> atoms1, atoms2;
    std::map<int,int> numbers1, numbers2;
    std::vector<int> partner1(g1.n, -1), partner2(g2.n, -1);
    int mapNum = 1;
    for (const auto& [qi, tj] : mapping) {
        atoms1.insert(qi);
        atoms2.insert(tj);
        numbers1[qi] = mapNum;
        numbers2[tj] = mapNum++;
        if (qi >= 0 && qi < g1.n && tj >= 0 && tj < g2.n) { partner1[qi] = tj; partner2[tj] = qi; }
    }

    // A bond is highlighted only when the mapping carries it onto a bond of the other molecule.
    auto carried = [](const MolGraph& g, const std::vector<int>& partner, const MolGraph& other) {
        std::set<std::pair<int,int>> bonds;
        for (int first = 0; first < g.n; ++first)
            for (int second : g.neighbors[first])
                if (first < second && partner[first] >= 0 && partner[second] >= 0
                        && other.bondOrder(partner[first], partner[second]) != 0) bonds.insert({first, second});
        return bonds;
    };
    const auto bonds1 = carried(g1, partner1, g2), bonds2 = carried(g2, partner2, g1);

    std::vector<detail_depict::AtomLabel> labels1(g1.n), labels2(g2.n);
    for (int i = 0; i < g1.n; i++) labels1[i] = detail_depict::buildLabel(g1, i, opts, pos1);
    for (int i = 0; i < g2.n; i++) labels2[i] = detail_depict::buildLabel(g2, i, opts, pos2);
    auto extent = [&](const std::vector<Vec2>& pos, const std::vector<detail_depict::AtomLabel>& labels,
                      const std::set<int>& atoms) {
        detail_depict::Bounds bounds;
        for (int i = 0; i < (int)pos.size(); i++)
            bounds.addAtom(pos[i], labels[i], atoms.count(i) ? opts.highlightRadius_() : 0);
        if (bounds.minX > bounds.maxX) bounds = {0, 0, 0, 0};
        return bounds;
    };
    const auto box1 = extent(pos1, labels1, atoms1), box2 = extent(pos2, labels2, atoms2);

    // Both molecules are centred on one horizontal axis, with the arrow's gap between them.
    const bool arrow = g1.n > 0 && g2.n > 0;
    const double gap = arrow ? opts.bondLength * 2.4 : 0;
    const double width1 = box1.maxX - box1.minX, width2 = box2.maxX - box2.minX;
    for (auto& p : pos1) p = {p.x - box1.minX, p.y - (box1.minY + box1.maxY) * 0.5};
    for (auto& p : pos2) p = {p.x - box2.minX + width1 + gap, p.y - (box2.minY + box2.maxY) * 0.5};
    const double halfHeight = std::max(box1.maxY - box1.minY, box2.maxY - box2.minY) * 0.5;
    const double minX = 0, maxX = width1 + gap + width2, minY = -halfHeight, maxY = halfHeight;
    double pad = opts.padding;
    int W = opts.width  > 0 ? opts.width  : std::max(100, (int)(maxX - minX + 2 * pad));
    int H = opts.height > 0 ? opts.height : std::max(100, (int)(maxY - minY + 2 * pad));

    // The view box keeps the canvas's aspect ratio and grows uniformly until the drawing fits.
    const double scale = std::max({1.0, (maxX - minX + 2 * pad) / W, (maxY - minY + 2 * pad) / H});
    const double viewWidth = W * scale, viewHeight = H * scale;
    const double offsetX = (viewWidth - (maxX - minX)) * 0.5 - minX;
    const double offsetY = (viewHeight - (maxY - minY)) * 0.5 - minY;
    std::vector<Vec2> at1(g1.n), at2(g2.n);
    for (int i = 0; i < g1.n; i++) at1[i] = {pos1[i].x + offsetX, pos1[i].y + offsetY};
    for (int i = 0; i < g2.n; i++) at2[i] = {pos2[i].x + offsetX, pos2[i].y + offsetY};

    std::ostringstream svg;
    svg << std::fixed << std::setprecision(2);
    detail_depict::openSVG(svg, W, H, viewWidth, viewHeight, opts);
    detail_depict::Scene scene;
    scene.view = {0, 0, viewWidth, viewHeight};
    detail_depict::drawMolecule(svg, scene, g1, at1, labels1, atoms1, bonds1, opts, false, scale);
    detail_depict::drawMolecule(svg, scene, g2, at2, labels2, atoms2, bonds2, opts, false, scale);

    // Double-headed arrow in the middle of the gap, with slim notched heads
    if (arrow) {
        const double y = offsetY, clearance = opts.bondLength * 0.5;
        const double left = offsetX + width1 + clearance, right = offsetX + width1 + gap - clearance;
        const double head = opts.bondLength * 0.3, barb = opts.bondLength * 0.1, notch = head * 0.7;
        const std::string colour = rgb(palette::ARROW);
        svg << "<line x1=\"" << (left + notch) << "\" y1=\"" << y << "\" x2=\"" << (right - notch) << "\" y2=\"" << y
            << "\" stroke=\"" << colour << "\" stroke-width=\"" << opts.lineWidth_() << "\"/>\n";
        for (double side : {-1.0, 1.0}) {
            const double tip = side > 0 ? right : left;
            svg << "<polygon points=\"" << tip << "," << y << " "
                << (tip - side * head) << "," << (y - barb) << " "
                << (tip - side * notch) << "," << y << " "
                << (tip - side * head) << "," << (y + barb)
                << "\" fill=\"" << colour << "\" stroke=\"none\"/>\n";
        }
        scene.strokes.push_back({{left, y}, {right, y}, barb + opts.fontSize_() * 0.3});
    }

    // The map numbers of both molecules are placed together, so the most crowded choose first.
    if (opts.showMapNumbers) {
        std::vector<detail_depict::MapNumber> numbers;
        detail_depict::addMapNumbers(numbers, g1, at1, numbers1);
        detail_depict::addMapNumbers(numbers, g2, at2, numbers2);
        detail_depict::drawMapNumbers(svg, scene, numbers, opts, scale);
    }

    svg << "</svg>\n";
    return svg.str();
}

// ── Convenience: SMILES → SVG ─────────────────────────────────────────────

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

} // namespace smsd

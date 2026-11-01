/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/depict.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

static void checkCoordinates(const smsd::MolGraph& graph,
                             const std::vector<smsd::Point2D>& coordinates, double bondLength) {
    for (int first = 0; first < graph.n; ++first)
        for (int second = first + 1; second < graph.n; ++second) {
            const double distance = std::hypot(coordinates[first].x - coordinates[second].x,
                                               coordinates[first].y - coordinates[second].y) / bondLength;
            if (distance < 0.4) throw std::runtime_error("Layout atoms overlap");
            if (graph.bondOrder(first, second) && (distance < 0.8 || distance > 1.2))
                throw std::runtime_error("Layout bond length is distorted");
        }
    if (smsd::detail_layout::countCrossings(graph, coordinates))
        throw std::runtime_error("Planar layout has crossing bonds");
}


static size_t occurrences(const std::string& text, const std::string& part) {
    size_t count = 0, offset = 0;
    while ((offset = text.find(part, offset)) != std::string::npos) { ++count; offset += part.size(); }
    return count;
}

static void checkKekuleBonds() {
    smsd::DepictOptions options;
    options.aromaticStyle = smsd::AromaticStyle::KEKULE;
    for (const char* smiles : {"c1ccccc1", "c1ccncc1", "c1ccc2ccccc2c1", "c1cc[nH]c1"}) {
        const auto graph = smsd::parseSMILES(smiles);
        const auto aromatic = graph.aromatic;
        int bonds = 0;
        for (int atom = 0; atom < graph.n; ++atom) bonds += graph.degree[atom];
        bonds /= 2;
        const auto single = smsd::depict(graph, options);
        const auto pair = smsd::depictPair(graph, graph, {}, options);
        if (occurrences(single, "<line") <= static_cast<size_t>(bonds)
                || occurrences(pair, "<line") <= static_cast<size_t>(2 * bonds))
            throw std::runtime_error("Kekule depiction must draw single and double ring bonds");
        if (graph.aromatic != aromatic) throw std::runtime_error("Depiction changed caller aromaticity");
        for (int atom = 0; atom < graph.n; ++atom) for (int neighbor : graph.neighbors[atom])
            if (graph.bondOrder(atom, neighbor) != 4) throw std::runtime_error("Depiction changed caller bond orders");
    }
}

static void checkMapNumbers() {
    const auto phenol = smsd::parseSMILES("c1ccccc1O"), benzene = smsd::parseSMILES("c1ccccc1");
    std::map<int, int> mapping;
    for (int atom = 0; atom < 7; ++atom) mapping.emplace(atom, atom);
    if (occurrences(smsd::depictWithMapping(phenol, mapping), "class=\"map\"") != 7)
        throw std::runtime_error("Every mapped atom needs its map number");
    mapping.erase(6);
    if (occurrences(smsd::depictPair(benzene, phenol, mapping), "class=\"map\"") != 12)
        throw std::runtime_error("Both mapped atoms of every pair need map numbers");
    smsd::DepictOptions options;
    options.showMapNumbers = false;
    if (occurrences(smsd::depictPair(benzene, phenol, mapping, options), "class=\"map\"") != 0)
        throw std::runtime_error("Map number option was ignored");
}

static void checkFontFamily() {
    const auto graph = smsd::parseSMILES("CCO");
    for (const char* family : {"a}</style><script>alert(1)</script><style>", "Arial; fill:red", "<>{};"}) {
        smsd::DepictOptions options; options.fontFamily = family;
        for (const auto& svg : {smsd::depict(graph, options), smsd::depictPair(graph, graph, {}, options)})
            if (svg.find("<script") != std::string::npos || svg.find("; fill:red") != std::string::npos
                    || svg.find("font-family: ;") != std::string::npos)
                throw std::runtime_error("Font family changed SVG or CSS structure");
    }
    smsd::DepictOptions options; options.fontFamily = "'Noto Sans', Helvetica, sans-serif";
    if (smsd::depict(graph, options).find(options.fontFamily) == std::string::npos)
        throw std::runtime_error("Ordinary font family lists must be preserved");
}

// Every ring of up to eight atoms keeps the interior angle of its regular polygon, 180(k-2)/k degrees.
static void checkRingsRegular(const smsd::MolGraph& graph, const std::vector<smsd::Point2D>& coordinates,
                              const std::string& smiles) {
    for (const auto& ring : smsd::computeSSSR(graph)) {
        const size_t size = ring.size();
        if (size > 8) continue;
        const double ideal = 180.0 * (size - 2) / size;
        for (size_t index = 0; index < size; ++index) {
            const auto centre = coordinates[ring[index]];
            const auto previous = coordinates[ring[(index + size - 1) % size]];
            const auto next = coordinates[ring[(index + 1) % size]];
            const double ax = previous.x - centre.x, ay = previous.y - centre.y;
            const double bx = next.x - centre.x, by = next.y - centre.y;
            const double cosine = (ax * bx + ay * by) / (std::hypot(ax, ay) * std::hypot(bx, by));
            const double angle = std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / M_PI;
            if (!std::isfinite(angle) || std::abs(angle - ideal) > 0.5)
                throw std::runtime_error("Ring is not a regular polygon in " + smiles);
        }
    }
}

// The closest approach of two non-bonded atoms, in the units of the coordinates.
static double closestNonBonded(const smsd::MolGraph& graph, const std::vector<smsd::Point2D>& coordinates) {
    double closest = std::numeric_limits<double>::infinity();
    for (int first = 0; first < graph.n; ++first)
        for (int second = first + 1; second < graph.n; ++second)
            if (!graph.bondOrder(first, second))
                closest = std::min(closest, std::hypot(coordinates[first].x - coordinates[second].x,
                                                       coordinates[first].y - coordinates[second].y));
    return closest;
}

// The smallest angle between two bonds at any atom, in degrees.
static double narrowestBondAngle(const smsd::MolGraph& graph, const std::vector<smsd::Point2D>& coordinates) {
    double narrowest = 360;
    for (int atom = 0; atom < graph.n; ++atom) {
        const auto& neighbours = graph.neighbors[atom];
        for (size_t first = 0; first < neighbours.size(); ++first)
            for (size_t second = first + 1; second < neighbours.size(); ++second) {
                const auto a = coordinates[neighbours[first]], b = coordinates[neighbours[second]], c = coordinates[atom];
                const double ax = a.x - c.x, ay = a.y - c.y, bx = b.x - c.x, by = b.y - c.y;
                const double cosine = (ax * bx + ay * by) / (std::hypot(ax, ay) * std::hypot(bx, by));
                narrowest = std::min(narrowest, std::acos(std::clamp(cosine, -1.0, 1.0)) * 180.0 / M_PI);
            }
    }
    return narrowest;
}

// Two bonds within ten degrees of a straight line leave drawn stereo unreadable (IUPAC 2006,
// ST-1.1 and ST-2): a tetrahedral centre drawn with three bonds has no two of them in a line,
// one drawn with four does not crowd them into a half-turn, and each substituent of a stereo
// double bond outside a ring stands off the bond's line, the two at one end on opposite sides,
// with the first ligands placed as the stored E/Z says.
static void checkReadableStereo(const smsd::MolGraph& graph, const std::vector<smsd::Point2D>& drawn,
                                const std::string& smiles) {
    auto direction = [&](int from, int to) { return std::atan2(drawn[to].y - drawn[from].y, drawn[to].x - drawn[from].x); };
    for (int centre = 0; centre < graph.n; ++centre) {
        const auto& neighbours = graph.neighbors[centre];
        if (graph.tetraChirality[centre] == 0 || neighbours.size() < 3) continue;
        std::vector<double> angles;
        for (int neighbour : neighbours) angles.push_back(direction(centre, neighbour));
        double straightest = 0;
        for (size_t first = 0; first < angles.size(); ++first)
            for (size_t second = first + 1; second < angles.size(); ++second)
                straightest = std::max(straightest, std::abs(std::remainder(angles[first] - angles[second], 2 * M_PI)));
        std::sort(angles.begin(), angles.end());
        double widest = angles.front() + 2 * M_PI - angles.back();
        for (size_t k = 1; k < angles.size(); ++k) widest = std::max(widest, angles[k] - angles[k - 1]);
        const double degrees = (neighbours.size() == 3 ? straightest : widest) * 180.0 / M_PI;
        if (degrees >= 170.0)
            throw std::runtime_error("A stereocentre is drawn with two bonds in a line in " + smiles);
    }
    for (int first = 0; first < graph.n; ++first)
        for (int second : graph.neighbors[first]) {
            if (second < first || graph.bondOrder(first, second) != 2 || graph.dbStereo(first, second) == 0
                    || graph.bondInRing(first, second)) continue;
            const double axisX = drawn[second].x - drawn[first].x, axisY = drawn[second].y - drawn[first].y;
            double relation = 0;
            for (const auto& [end, partner] : {std::make_pair(first, second), std::make_pair(second, first)}) {
                double previous = 0;
                for (int ligand : graph.neighbors[end]) {
                    if (ligand == partner) continue;
                    const double turn = std::remainder(direction(end, ligand) - direction(end, partner), 2 * M_PI);
                    if (std::abs(turn) * 180.0 / M_PI >= 170.0)
                        throw std::runtime_error("A stereo double bond has a substituent in line with it in " + smiles);
                    const double side = axisX * (drawn[ligand].y - drawn[end].y) - axisY * (drawn[ligand].x - drawn[end].x);
                    if (previous != 0 && previous * side > 0)
                        throw std::runtime_error("Both substituents of a double bond end lie on one side in " + smiles);
                    if (previous == 0) relation = relation == 0 ? side : relation * side;
                    previous = side;
                }
            }
            if ((relation > 0 ? 1 : 2) != graph.dbStereo(first, second))
                throw std::runtime_error("A double bond is drawn as the other isomer in " + smiles);
        }
}

static void checkReadableStereoCentres() {
    // Each drawing leaves its stereo readable.  Hyoscyamine's C3 ester must not be drawn in line
    // with C2, where it reads as the 3-epimer; nor saxagliptin's or a hydroxyethylamine's centres,
    // the thioxanthene or epothilone E/Z.
    for (const char* smiles : {
            "CN1[C@@H]2CC[C@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1",                                  // hyoscyamine
            "CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c1ccccc1)C(=O)OC",                            // cocaine
            "N#C[C@@H]1C[C@@H]2C[C@@H]2N1C(=O)[C@@H]([NH3+])C12CC3CC(CC(O)(C3)C1)C2",             // saxagliptin
            "CCC1:C:N2:C3:C(:C:C(C(=O)NC(CC4:C:C:C:C:C:4)[C@H](O)C[NH2+]CC4CCOCC4):C:C:1:3)OS(=O)(=O)CC2",
            "CC(=O)S[C@@H]1CC2=CC(=O)CC[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](CC[C@@]34CCC(=O)O4)[C@H]12",  // spironolactone
            "CCC[C@@H]1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@@]5(C)[C@H]4[C@@H](O)C[C@@]3(C)[C@@]2(O1)C(=O)CO",
            "Clc1ccc2c(c1)/C(=C/C1CN3CCC1CC3)c1ccccc1S2",
            "CN1C(=O)[C@@](c2ccc(OC(F)F)cc2)(c2cccc(/C=C/CCO)c2)N=C1N",
            "C/C(=C\\c1csc(C)n1)[C@@H]1C[C@@H]2O[C@]2(C)CCC[C@H](C)[C@H](O)[C@@H](C)C(=O)C(C)(C)[C@@H](O)CC(=O)N1"}) {
        const auto graph = smsd::parseSMILES(smiles);
        checkReadableStereo(graph, smsd::generateCoords2D(graph, 1.0), smiles);
    }
}

static void checkRegularRingCoordinates() {
    // Rings stay regular polygons and no two non-bonded atoms come closer than 0.8 bond lengths:
    // ortho substituents (aspirin) and angular methyls (testosterone) must not bend the rings
    // they clash with.
    for (const char* smiles : {"c1ccccc1", "Oc1ccccc1", "CC(=O)Nc1ccc(O)cc1", "c1ccc2ccccc2c1",
                               "c1ccc2[nH]ccc2c1", "c1ccc2cccc2cc1", "CC(=O)Oc1ccccc1C(=O)O",
                               "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O"}) {
        const auto graph = smsd::parseSMILES(smiles);
        const auto coordinates = smsd::generateCoords2D(graph, 1.0);
        checkRingsRegular(graph, coordinates, smiles);
        if (closestNonBonded(graph, coordinates) < 0.8)
            throw std::runtime_error(std::string("Non-bonded atoms clash in ") + smiles);
    }
}

static void checkNoHiddenAtoms() {
    // No atom may lie under another.  The angular methyl of a 17-ester corticosteroid must stay
    // off C1, a 16,17-acetal must not fold onto the C17 side chain, and codeine's bridges must
    // not cover each other.  The steroids and cloxazolam also keep regular rings and no crossing
    // bonds, and no bond of theirs folds to within 30 degrees of another; cloxazolam's pendant
    // aryl must not cover the fused benzo ring.
    struct Case { const char* smiles; bool regular; };
    for (const Case& molecule : {
            Case{"CCC(=O)OCC(=O)C1(OC(=O)CC)C(C)CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC21C", true},  // betamethasone dipropionate
            Case{"CCC(=O)OC1(C(=O)SCF)C(C)CC2C3CC(F)C4=CC(=O)C=CC4(C)C3(F)C(O)CC21C", true},     // fluticasone propionate
            Case{"CCCC1OC2CC3C4CCC5=CC(=O)C=CC5(C)C4C(O)CC3(C)C2(O1)C(=O)CO", true},             // budesonide
            Case{"CC1(C)OC2CC3C4CCC5=CC(=O)C=CC5(C)C4(F)C(O)CC3(C)C2(C(=O)CO)O1", true},         // triamcinolone acetonide
            Case{"O=C1CN2CCOC2(c2ccccc2Cl)c2cc(Cl)ccc2N1", true},                                // cloxazolam
            Case{"CN1CCC23c4c5ccc(OC)c4OC2C(O)C=CC3C1C5", false}}) {                           // codeine
        const auto graph = smsd::parseSMILES(molecule.smiles);
        const auto coordinates = smsd::generateCoords2D(graph, 1.0);
        if (closestNonBonded(graph, coordinates) < 0.4)
            throw std::runtime_error(std::string("An atom is hidden under another in ") + molecule.smiles);
        if (!molecule.regular) continue;
        checkRingsRegular(graph, coordinates, molecule.smiles);
        if (smsd::detail_layout::countCrossings(graph, coordinates) != 0)
            throw std::runtime_error(std::string("Bonds cross in ") + molecule.smiles);
        if (narrowestBondAngle(graph, coordinates) < 30)
            throw std::runtime_error(std::string("A bond folds onto another in ") + molecule.smiles);
    }
}

static void checkSpiroRings() {
    // Ring systems that meet at one atom keep clear of each other, with no crossing bonds and no
    // atom nearer another than half a bond length.  Two systems at a spiro centre also keep regular
    // rings: fluorescein's phthalide must turn clear of its xanthene.  Round a chelating metal,
    // where three rings cannot all be regular, the rings give way instead.
    struct Case { const char* smiles; bool regular; };
    for (const Case& molecule : {
            Case{"C1CCC2(CC1)CCCC2", true},
            Case{"CN1C(=O)C2(CC(C)(C)Cc3ccc(-c4cc(F)cc(Cl)c4)cc32)N=C1N", true},
            Case{"O=C1OC2(c3ccccc31)c1ccc(O)cc1Oc1cc(O)ccc12", true},                           // fluorescein
            Case{"O=C1OC2(c3cc(I)c([O-])c(I)c3Oc3c2cc(I)c([O-])c3I)c2c(Cl)c(Cl)c(Cl)c(Cl)c21", true},  // rose bengal
            Case{"Cl[Co+]12(Cl)NCCN1.NCCN2", true},
            Case{"[Co+3]123(NCCN1)(NCCN2)NCCN3", false},                                        // [Co(en)3]3+
            Case{"[Fe-3]123(OC(=O)C(=O)O1)(OC(=O)C(=O)O2)OC(=O)C(=O)O3", false},                // [Fe(ox)3]3-
            Case{"[Fe]123(OC(C)=CC(C)O1)(OC(C)=CC(C)O2)OC(C)=CC(C)O3", false}}) {               // [Fe(acac)3]
        const auto graph = smsd::parseSMILES(molecule.smiles);
        const auto coordinates = smsd::generateCoords2D(graph, 1.0);
        if (molecule.regular) checkRingsRegular(graph, coordinates, molecule.smiles);
        if (smsd::detail_layout::countCrossings(graph, coordinates) != 0)
            throw std::runtime_error(std::string("Bonds cross where ring systems meet in ") + molecule.smiles);
        if (closestNonBonded(graph, coordinates) < 0.5)
            throw std::runtime_error(std::string("Ring systems that meet at one atom overlap in ") + molecule.smiles);
    }
}

static void checkFoldedZChain() {
    // A long all-cis chain folds back on itself.  Each Z bond is still drawn as stored, with no
    // crossing bond and no atom nearer another than 0.4 bond lengths.
    std::string smiles = "C";
    for (int unit = 0; unit < 200; ++unit) smiles += "/C=C\\C";
    const auto graph = smsd::parseSMILES(smiles);
    const auto coordinates = smsd::generateCoords2D(graph, 1.0);
    if (smsd::detail_layout::countCrossings(graph, coordinates) != 0)
        throw std::runtime_error("Bonds cross in a folded all-cis chain");
    if (closestNonBonded(graph, coordinates) < 0.4)
        throw std::runtime_error("An atom is hidden in a folded all-cis chain");
    checkReadableStereo(graph, coordinates, "a folded all-cis chain");
}

static void checkStoredDoubleBondGeometry() {
    // Drawings show the stored isomer of every acyclic stereo double bond: the first ligands at
    // its two ends lie on one side of the bond when the stored relation is 1, and on opposite
    // sides when it is 2.  Maleic acid must not be drawn as fumaric acid, nor cis-stilbene as trans.
    for (const char* smiles : {"OC(=O)/C=C\\C(=O)O", "OC(=O)/C=C/C(=O)O", "c1ccccc1/C=C\\c1ccccc1",
                               "c1ccccc1/C=C/c1ccccc1", "C/C=C\\C", "CC/C=C\\CCO",
                               "CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O",
                               "CO/N=C(\\C(=O)NC1C(=O)N2C(C(=O)O)=C(COC(C)=O)CSC12)c1csc(N)n1"}) {
        const auto graph = smsd::parseSMILES(smiles);
        const auto drawn = smsd::layout2D(graph);
        {
            std::vector<smsd::Point2D> points;
            for (const auto& point : drawn) points.push_back({point.x, point.y});
            checkReadableStereo(graph, points, smiles);
        }
        auto firstLigand = [&](int end, int partner) {
            for (int neighbour : graph.neighbors[end])
                if (neighbour != partner) return neighbour;
            return -1;
        };
        int checked = 0;
        for (int first = 0; first < graph.n; ++first)
            for (int second : graph.neighbors[first]) {
                if (second < first || graph.bondOrder(first, second) != 2 || graph.dbStereo(first, second) == 0)
                    continue;
                const auto axis = drawn[second] - drawn[first];
                auto side = [&](int end, int ligand) {
                    const auto arm = drawn[ligand] - drawn[end];
                    return axis.x * arm.y - axis.y * arm.x;
                };
                const double relation = side(first, firstLigand(first, second)) * side(second, firstLigand(second, first));
                if ((relation > 0 ? 1 : 2) != graph.dbStereo(first, second))
                    throw std::runtime_error(std::string("A double bond is drawn as the other isomer in ") + smiles);
                ++checked;
            }
        if (checked == 0) throw std::runtime_error(std::string("No stereo double bond was read in ") + smiles);
    }
}

static void checkRingDoubleBondsInMolFiles() {
    // Drawings follow the stored E/Z of double bonds outside rings; inside a large ring they may
    // not yet, but the coordinates written to a MOL file keep it, so the file reads back as stored.
    for (const char* smiles : {"C1CCCCC/C=C/CCCC1", "C1CCCCC/C=C\\CCCC1"}) {
        const auto graph = smsd::parseSMILES(smiles);
        const auto back = smsd::readMolBlock(smsd::writeMolBlock(graph));
        if (back.n != graph.n) throw std::runtime_error(std::string("A MOL file lost atoms of ") + smiles);
        for (int first = 0; first < graph.n; ++first)
            for (int second : graph.neighbors[first])
                if (graph.bondOrder(first, second) == 2 && back.dbStereo(first, second) != graph.dbStereo(first, second))
                    throw std::runtime_error(std::string("A MOL file lost the ring E/Z of ") + smiles);
    }
}

static void checkGridOrientation() {
    // Canonical orientation puts the bonds on the 30-degree grid.  Benzene has no principal axis,
    // and a turn chosen by rounding noise would draw it at a different slant on every build.  The
    // layout rounds coordinates to 2^-20 bond lengths, so a generated hexagon's bonds may stray
    // from the grid by up to a ten-thousandth of a degree.
    const auto benzene = smsd::parseSMILES("c1ccccc1");
    std::vector<smsd::Point2D> hexagon;
    for (int corner = 0; corner < 6; ++corner)
        hexagon.push_back({std::cos(0.2 + corner * M_PI / 3), std::sin(0.2 + corner * M_PI / 3)});
    for (const auto& coordinates : {smsd::generateCoords2D(benzene, 1.0), hexagon}) {
        auto turned = coordinates;
        smsd::transform::canonicalOrientation(benzene, turned);
        for (int first = 0; first < 6; ++first) {
            const int second = (first + 1) % 6;
            const double degrees = std::atan2(turned[second].y - turned[first].y, turned[second].x - turned[first].x)
                                 * 180.0 / M_PI;
            if (std::abs(std::remainder(degrees, 30.0)) > 1e-4)
                throw std::runtime_error("Benzene's bonds must lie on the 30-degree grid");
        }
    }
}

static std::vector<double> attributeValues(const std::string& svg, const std::string& name) {
    std::vector<double> values;
    const std::string key = " " + name + "=\"";
    for (size_t offset = svg.find(key); offset != std::string::npos; offset = svg.find(key, offset + 1))
        values.push_back(std::stod(svg.substr(offset + key.size())));
    return values;
}

static void checkFixedCanvasFit() {
    // A drawing larger than a fixed canvas is scaled down uniformly: the view box
    // keeps the canvas's aspect ratio and every mark stays inside it.
    smsd::DepictOptions options;
    options.width = 240;
    options.height = 160;
    const auto chain = smsd::parseSMILES("CCCCCCCCCCCCCCCCCCCCO");
    const auto naphthol = smsd::parseSMILES("Oc1ccc2ccccc2c1");
    std::map<int, int> ends{{0, 0}, {20, 1}}, rings;
    for (int atom = 0; atom < 11; ++atom) rings.emplace(atom, atom);
    for (const auto& svg : {smsd::depictWithMapping(chain, ends, options),
                            smsd::depictPair(naphthol, naphthol, rings, options)}) {
        std::istringstream box(svg.substr(svg.find("viewBox=\"") + 9));
        double left, top, width, height;
        box >> left >> top >> width >> height;
        if (width <= options.width || std::abs(width / height - 1.5) > 1e-3)
            throw std::runtime_error("A fixed canvas must scale uniformly to fit the drawing");
        const auto radii = attributeValues(svg, "r");
        double radius = 0;
        for (double r : radii) radius = std::max(radius, r);
        for (const char* name : {"x1", "x2", "cx", "x"})
            for (double x : attributeValues(svg, name))
                if (x - radius < left - 1e-6 || x + radius > left + width + 1e-6)
                    throw std::runtime_error("A mark falls outside the fixed canvas");
        for (const char* name : {"y1", "y2", "cy", "y"})
            for (double y : attributeValues(svg, name))
                if (y - radius < top - 1e-6 || y + radius > top + height + 1e-6)
                    throw std::runtime_error("A mark falls outside the fixed canvas");
    }
}

static void checkBondHighlightBands() {
    const auto phenol = smsd::parseSMILES("c1ccccc1O");
    const auto highlighted = smsd::depictWithMapping(phenol, {{0,0},{1,1}});
    if (occurrences(highlighted, "class=\"bond-highlight\"") != 1)
        throw std::runtime_error("Matched neighbouring atoms need a continuous bond highlight");
    if (occurrences(smsd::depictWithMapping(phenol, {{0,0},{3,1}}), "class=\"bond-highlight\"") != 0)
        throw std::runtime_error("A highlight must not create a bond between non-neighbouring atoms");
    std::map<int,int> mapping;
    for (int atom = 0; atom < 6; ++atom) mapping.emplace(atom, atom);
    if (occurrences(smsd::depictPair(phenol, phenol, mapping), "class=\"bond-highlight\"") != 12)
        throw std::runtime_error("Both molecules need their matched ring bonds highlighted");
    // Hexane maps onto five of cyclohexane's six bonds; the ring closure is not matched.
    if (occurrences(smsd::depictPair(smsd::parseSMILES("CCCCCC"), smsd::parseSMILES("C1CCCCC1"), mapping),
                    "class=\"bond-highlight\"") != 10)
        throw std::runtime_error("Only bonds carried by the mapping may be highlighted");
    // A bond highlight colour sets the halo's hue; an opaque colour takes the match fill's
    // strength, so pure red at 0.40 over white is drawn flat as #FF9999.
    smsd::DepictOptions options;
    options.bondHighlightColors[{0, 1}] = {255, 0, 0, 1.0};
    const auto coloured = smsd::depictWithMapping(phenol, {{0,0},{1,1}}, options);
    const size_t band = coloured.find("class=\"bond-highlight\"");
    if (band == std::string::npos || coloured.find("stroke=\"#FF9999\"", band) == std::string::npos
            || coloured.find("rgba(") != std::string::npos || coloured.find("opacity=") != std::string::npos)
        throw std::runtime_error("A bond highlight colour must give a flat tint of its hue");
}

static std::string between(const std::string& text, const std::string& open, const std::string& close) {
    const size_t start = text.find(open);
    if (start == std::string::npos) return "";
    return text.substr(start, text.find(close, start) - start);
}

/** Stroke colours of every line, bond lines and the pair arrow's shaft. */
static std::set<std::string> lineColours(const std::string& svg) {
    std::set<std::string> colours;
    for (size_t at = svg.find("<line"); at != std::string::npos; at = svg.find("<line", at + 1)) {
        const size_t stroke = svg.find("stroke=\"", at);
        colours.insert(svg.substr(stroke + 8, 7));
    }
    return colours;
}

static void checkHighlightLayer() {
    const auto phenol = smsd::parseSMILES("c1ccccc1O");
    std::map<int, int> mapping;
    for (int atom = 0; atom < 7; ++atom) mapping.emplace(atom, atom);
    // Each molecule's halo is one layer of opaque shapes: the gold tint (MATCH_FILL, 0.40 over
    // white) and its gold edge, with the fully matched ring filled.
    for (const auto& [svg, molecules] : {std::make_pair(smsd::depictWithMapping(phenol, mapping), size_t{1}),
                                         std::make_pair(smsd::depictPair(phenol, phenol, mapping), size_t{2})}) {
        const std::string layer = between(svg, "<g class=\"highlight\">", "</g>");
        if (occurrences(svg, "<g class=\"highlight\">") != molecules || layer.empty())
            throw std::runtime_error("Each molecule's halo must be drawn as one layer");
        if (layer.find("rgba(") != std::string::npos || layer.find("opacity") != std::string::npos)
            throw std::runtime_error("The halo must be opaque, so overlapping shapes keep one tint");
        std::set<std::string> colours;
        for (const char* key : {"fill=\"#", "stroke=\"#"})
            for (size_t at = layer.find(key); at != std::string::npos; at = layer.find(key, at + 1))
                colours.insert(layer.substr(at + std::string(key).size(), 6));
        if (colours != std::set<std::string>{"F9E6B3", "F0C246"})
            throw std::runtime_error("The default halo must be the gold tint with a gold edge");
        if (occurrences(layer, "class=\"ring-highlight\"") != 1)
            throw std::runtime_error("A fully matched ring must be filled");
        // The skeleton stays near-black, matched bonds included; the pair arrow is grey-green.
        for (const auto& colour : lineColours(svg))
            if (colour != "#212121" && colour != "#80928F")
                throw std::runtime_error("Bond lines must keep the bond colour");
        if (svg.find("#4CAF50") != std::string::npos || svg.find("#2E7D32") != std::string::npos)
            throw std::runtime_error("The former green highlight must not return");
    }
}

static void checkLabelKnockouts() {
    const auto phenol = smsd::parseSMILES("c1ccccc1O");
    // The oxygen of a matched hydroxyl sits on the gold tint itself, with no white bead round it.
    const auto onHalo = smsd::depictWithMapping(phenol, {{5, 0}, {6, 1}});
    const std::string knockout = between(onHalo, "<g class=\"knockout\"", "</g>");
    if (knockout.find("fill=\"#F9E6B3\" stroke=\"#F9E6B3\">O</text>") == std::string::npos)
        throw std::runtime_error("A label on the halo must be knocked out in the tint beneath it");
    if (smsd::depictWithMapping(phenol, {{0, 0}}).find("class=\"knockout\"") != std::string::npos)
        throw std::runtime_error("A label off the halo needs no knockout");
    // On a tint too dark for the label's ink, the knockout keeps the background colour.
    smsd::DepictOptions options;
    options.atomHighlightColors[6] = {13, 51, 48, 0.9};
    const auto dark = smsd::depictSVG(phenol, smsd::layout2D(phenol, options), {}, {6}, {}, options).svg;
    if (between(dark, "<g class=\"knockout\"", "</g>").find("fill=\"#FFFFFF\" stroke=\"#FFFFFF\">O</text>") == std::string::npos)
        throw std::runtime_error("A label on a dark tint must keep a background knockout");
}

static void checkHighlightOverrides() {
    // Per-atom and per-bond colours tint the halo; the skeleton and element inks are unchanged,
    // and an element colour override still sets the label's ink.
    const auto phenol = smsd::parseSMILES("c1ccccc1O");
    smsd::DepictOptions options;
    options.atomHighlightColors[0] = {120, 184, 240, 0.5};
    options.bondHighlightColors[{0, 1}] = {255, 0, 0, 1.0};
    options.atomColorOverrides[8] = {0, 128, 0, 1.0};
    const auto svg = smsd::depictSVG(phenol, smsd::layout2D(phenol, options), {}, {0, 1}, {{0, 1}}, options).svg;
    const std::string layer = between(svg, "<g class=\"highlight\">", "</g>");
    if (layer.find("fill=\"#BCDCF8\"") == std::string::npos)
        throw std::runtime_error("An atom highlight colour must tint that atom's halo");
    if (layer.find("stroke=\"#FF9999\"") == std::string::npos)
        throw std::runtime_error("A bond highlight colour must tint that bond's halo");
    if (lineColours(svg) != std::set<std::string>{"#212121"})
        throw std::runtime_error("Highlight colours must not recolour bond lines");
    if (svg.find("fill=\"#008000\">O</text>") == std::string::npos)
        throw std::runtime_error("An element colour override must set the label colour");
}

static void checkHighlightFollowsBonds() {
    const auto benzene = smsd::parseSMILES("c1ccccc1");
    const auto positions = smsd::layout2D(benzene);
    // An explicit bond list is followed exactly; without one, bonds joining highlighted atoms are drawn.
    if (occurrences(smsd::depictSVG(benzene, positions, {}, {0, 1, 2}, {{0, 1}}).svg, "class=\"bond-highlight\"") != 1)
        throw std::runtime_error("depictSVG must highlight exactly the bonds it is given");
    if (occurrences(smsd::depictSVG(benzene, positions, {}, {0, 1, 2}, {}).svg, "class=\"bond-highlight\"") != 2)
        throw std::runtime_error("Without a bond list, bonds joining highlighted atoms are highlighted");
    // Only mapped atoms carry a halo in a pair.
    std::map<int, int> mapping{{0, 0}, {1, 1}};
    const auto pair = smsd::depictPair(smsd::parseSMILES("CCO"), smsd::parseSMILES("CCN"), mapping);
    if (occurrences(pair, "fill=\"#F9E6B3\"/>") != 4)
        throw std::runtime_error("Only the mapped atoms of a pair may carry a halo");
}

static void checkKekuleDefault() {
    // Rings are drawn as alternating double bonds by default; circles remain available, and a
    // graph without a Kekule structure falls back to them.
    if (smsd::depictSMILES("c1ccccc1").find("<circle") != std::string::npos)
        throw std::runtime_error("Aromatic rings must default to Kekule bonds");
    smsd::DepictOptions circles;
    circles.aromaticStyle = smsd::AromaticStyle::CIRCLE;
    if (smsd::depictSMILES("c1ccccc1", circles).find("<circle") == std::string::npos)
        throw std::runtime_error("The circle style must still draw ring circles");
    smsd::DepictOptions kekule;
    kekule.aromaticStyle = smsd::AromaticStyle::KEKULE;
    if (smsd::depictSMILES("c1cccc1", kekule).find("<circle") == std::string::npos)
        throw std::runtime_error("A ring without a Kekule structure must fall back to a circle");
}

/** Centres of the gold halo discs, which follow atom order in each molecule. */
static std::vector<std::pair<double, double>> haloDiscs(const std::string& svg) {
    std::vector<std::pair<double, double>> discs;
    for (size_t at = svg.find("<circle cx=\""); at != std::string::npos; at = svg.find("<circle cx=\"", at + 1))
        if (svg.find("fill=\"#F9E6B3\"", at) < svg.find("/>", at))
            discs.push_back({std::stod(svg.substr(at + 12)), std::stod(svg.substr(svg.find("cy=\"", at) + 4))});
    return discs;
}

struct SvgText { double x, y, size; std::string cls, text; bool outline, centred; };

static std::vector<SvgText> svgTexts(const std::string& svg) {
    std::vector<SvgText> texts;
    for (size_t at = svg.find("<text "); at != std::string::npos; at = svg.find("<text ", at + 1)) {
        const size_t close = svg.find('>', at), end = svg.find("</text>", close);
        const std::string tag = svg.substr(at, close - at);
        auto attribute = [&](const std::string& name) {
            const std::string key = " " + name + "=\"";
            const size_t start = tag.find(key);
            return start == std::string::npos ? std::string()
                : tag.substr(start + key.size(), tag.find('"', start + key.size()) - start - key.size());
        };
        texts.push_back({std::stod(attribute("x")), std::stod(attribute("y")), std::stod(attribute("font-size")),
                         attribute("class"), svg.substr(close + 1, end - close - 1),
                         !attribute("stroke").empty(), attribute("text-anchor") == "middle"});
    }
    return texts;
}

static void checkMapNumberPlacement() {
    // Fully mapped molecules crowd every label: each map number must stay off the labels, must not
    // sit level beside one where it reads as a subscript, must be nearer its own atom than any
    // other, and two level numbers must keep apart.  Inside a ring, where a crowded number may
    // fall back to, it cannot be read as part of a formula.
    for (const char* smiles : {"Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "C[N+](C)(C)CCOC(C)=O", "OC(=O)c1ccccc1O",
                               "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12"}) {
        const auto graph = smsd::parseSMILES(smiles);
        std::map<int, int> mapping;
        for (int atom = 0; atom < graph.n; ++atom) mapping.emplace(atom, atom);
        const auto svg = smsd::depictWithMapping(graph, mapping);
        const auto atoms = haloDiscs(svg);
        struct Box { double left, top, right, bottom; };
        std::vector<Box> labels, numbers;
        double font = 0;
        for (const auto& text : svgTexts(svg)) {
            const double width = 0.6 * text.size * text.text.size(), half = 0.36 * text.size;
            if (text.cls == "map") numbers.push_back({text.x - width / 2, text.y - half, text.x + width / 2, text.y + half});
            else if (text.cls == "sub" && !text.outline && !text.centred) {
                labels.push_back({text.x, text.y - half, text.x + width, text.y + half});
                font = std::max(font, text.size);
            }
        }
        if ((int)numbers.size() != graph.n || (int)atoms.size() != graph.n)
            throw std::runtime_error(std::string("Every mapped atom needs its halo and map number in ") + smiles);
        std::vector<std::pair<std::pair<double, double>, double>> rings;  // centre and inner radius
        for (const auto& ring : smsd::computeSSSR(graph)) {
            double x = 0, y = 0;
            for (int atom : ring) { x += atoms[atom].first / ring.size(); y += atoms[atom].second / ring.size(); }
            rings.push_back({{x, y}, std::hypot(atoms[ring[0]].first - x, atoms[ring[0]].second - y) * std::cos(M_PI / ring.size())});
        }
        auto across = [](const Box& a, const Box& b) { return std::max(a.left - b.right, b.left - a.right); };
        auto down = [](const Box& a, const Box& b) { return std::max(a.top - b.bottom, b.top - a.bottom); };
        for (size_t k = 0; k < numbers.size(); ++k) {
            const Box& number = numbers[k];
            const double cx = (number.left + number.right) / 2, cy = (number.top + number.bottom) / 2;
            const bool inRing = std::any_of(rings.begin(), rings.end(), [&](const auto& ring) {
                return std::hypot(cx - ring.first.first, cy - ring.first.second) < ring.second;
            });
            for (const auto& label : labels) {
                if (across(number, label) < 0 && down(number, label) < 0)
                    throw std::runtime_error(std::string("A map number covers a label in ") + smiles);
                if (!inRing && down(number, label) < 0 && across(number, label) < font * 0.5)
                    throw std::runtime_error(std::string("A map number reads as part of a label in ") + smiles);
            }
            for (size_t other = 0; other < numbers.size(); ++other)
                if (other != k && down(number, numbers[other]) < 0 && across(number, numbers[other]) < font * 0.5)
                    throw std::runtime_error(std::string("Two map numbers run together in ") + smiles);
            auto distance = [&](size_t atom) { return std::hypot(cx - atoms[atom].first, cy - atoms[atom].second); };
            for (size_t atom = 0; atom < atoms.size(); ++atom)
                if (atom != k && distance(atom) < distance(k))
                    throw std::runtime_error(std::string("A map number lies nearer another atom in ") + smiles);
        }
    }
}

static void checkPairWedges() {
    // Both molecules of a pair keep their wedges, and the wedge gives the drawn centre its
    // configuration: (S) for the S pair and (R) for the R pair.  A solid wedge has its points;
    // a hashed one is a row of short lines across the bond.
    const std::map<int, int> mapping{{0, 0}, {1, 1}, {2, 2}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}, {9, 10}};
    auto numbers = [](std::string text) {
        std::replace(text.begin(), text.end(), ',', ' ');
        std::istringstream stream(text);
        return std::vector<double>{std::istream_iterator<double>(stream), std::istream_iterator<double>()};
    };
    for (const bool sinister : {true, false}) {
        const auto amphetamine = smsd::parseSMILES(sinister ? "C[C@H](N)Cc1ccccc1" : "C[C@@H](N)Cc1ccccc1");
        const auto methamphetamine = smsd::parseSMILES(sinister ? "C[C@H](NC)Cc1ccccc1" : "C[C@@H](NC)Cc1ccccc1");
        const auto svg = smsd::depictPair(amphetamine, methamphetamine, mapping);
        // The second molecule has no disc for its unmapped N-methyl.
        const auto discs = haloDiscs(svg);
        std::vector<std::vector<double>> solid, lines;
        for (size_t at = svg.find("<polygon points=\""); at != std::string::npos; at = svg.find("<polygon points=\"", at + 1)) {
            const auto values = numbers(svg.substr(at + 17, svg.find('"', at + 17) - at - 17));
            if (values.size() == 6 && svg.find("fill=\"#212121\"", at) < svg.find("/>", at)) solid.push_back(values);
        }
        for (size_t at = svg.find("<line x1=\""); at != std::string::npos; at = svg.find("<line x1=\"", at + 1)) {
            std::vector<double> line;
            for (const char* key : {"x1=\"", "y1=\"", "x2=\"", "y2=\""}) line.push_back(std::stod(svg.substr(svg.find(key, at) + 4)));
            lines.push_back(line);
        }
        if (discs.size() != 20) throw std::runtime_error("Only the mapped atoms of the pair may carry a halo");
        // Centre, then the neighbours in CIP order: N, CH2Ar, CH3.  Disc indices of each molecule.
        const int order[2][4] = {{1, 2, 3, 0}, {11, 12, 13, 10}};
        for (const auto& centre : order) {
            const double x = discs[centre[0]].first, y = discs[centre[0]].second;
            auto towards = [&](double px, double py, int disc) {
                const double bx = discs[disc].first - x, by = discs[disc].second - y, length = std::hypot(bx, by);
                const double along = ((px - x) * bx + (py - y) * by) / length, off = ((px - x) * by - (py - y) * bx) / length;
                return along > 0 && along < length + 0.5 && std::abs(off) < 0.5;
            };
            int wedged = -1;
            bool hashed = false;
            for (const auto& wedge : solid)
                if (std::hypot(wedge[0] - x, wedge[1] - y) < 1e-6)
                    for (int k = 1; k < 4; ++k)
                        if (towards((wedge[2] + wedge[4]) / 2, (wedge[3] + wedge[5]) / 2, centre[k])) wedged = k;
            for (int k = 1; k < 4 && wedged < 0; ++k) {
                int stripes = 0;
                for (const auto& line : lines)
                    stripes += towards((line[0] + line[2]) / 2, (line[1] + line[3]) / 2, centre[k])
                            && !towards(line[0], line[1], centre[k]);
                if (stripes >= 3) { wedged = k; hashed = true; }
            }
            if (wedged < 0) throw std::runtime_error("Both molecules of a pair must draw their stereo wedge");
            double angles[3];
            for (int k = 0; k < 3; ++k) angles[k] = std::atan2(-(discs[centre[k + 1]].second - y), discs[centre[k + 1]].first - x);
            auto turn = [](double from, double to) { return std::fmod(to - from + 4 * M_PI, 2 * M_PI); };
            // A solid wedge puts the implicit hydrogen behind, so anticlockwise N, CH2Ar, CH3 is S;
            // a hashed wedge puts it in front and reverses the sense.
            const bool anticlockwise = turn(angles[0], angles[1]) < turn(angles[0], angles[2]);
            if ((anticlockwise != hashed) != sinister)
                throw std::runtime_error(sinister ? "The pair wedge must give (S)" : "The pair wedge must give (R)");
        }
    }
}

static void checkIsotopeLabels() {
    // An isotope keeps its mass number, set small before the symbol: 13C, 2H and 131I.
    for (const auto& [smiles, mass] : {std::make_pair("[13CH3]OC", "13"), std::make_pair("[2H]OC", "2"),
                                       std::make_pair("[131I]c1ccccc1", "131")}) {
        const auto texts = svgTexts(smsd::depictSMILES(smiles));
        const auto symbol = std::find_if(texts.begin(), texts.end(), [&](const SvgText& text) {
            return text.cls == "sub" && !text.outline && text.text == mass;
        });
        if (symbol == texts.end()) throw std::runtime_error(std::string("An isotope lost its mass number in ") + smiles);
    }
}

static void checkWedgesOutsideRings() {
    // A wedge goes on a bond outside the rings where one is free (IUPAC ST-1.2): hyoscyamine's
    // C3 takes it on the bond to the ester oxygen, testosterone's C17 on the bond to its hydroxyl.
    for (const auto& [smiles, centre, other] : {std::make_tuple("CN1[C@@H]2CC[C@H]1C[C@@H](C2)OC(=O)C(CO)c1ccccc1", 7, 9),
                                                std::make_tuple("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O", 19, 20)}) {
        const auto graph = smsd::parseSMILES(smiles);
        const auto wedges = smsd::mol_detail::tetrahedralWedges(graph, smsd::generateCoords2D(graph, 1.0));
        const auto wedge = wedges.find(smsd::MolGraph::bondKey(centre, other));
        if (wedge == wedges.end() || wedge->second.first != centre)
            throw std::runtime_error(std::string("A ring bond took the wedge from a free bond in ") + smiles);
    }
}

static void checkTransparentHighlight() {
    // A highlight colour of strength 0 draws nothing at all, not even the halo's edge.
    const auto phenol = smsd::parseSMILES("c1ccccc1O");
    smsd::DepictOptions options;
    options.atomHighlightColors[6] = {240, 192, 64, 0.0};
    const auto svg = smsd::depictSVG(phenol, smsd::layout2D(phenol, options), {}, {6}, {}, options).svg;
    if (svg.find("<g class=\"highlight\">") != std::string::npos)
        throw std::runtime_error("A transparent highlight must draw no halo");
}

/** The view box of an SVG: left, top, width and height. */
static std::vector<double> viewBox(const std::string& svg) {
    const size_t start = svg.find("viewBox=\"") + 9;
    std::istringstream box(svg.substr(start, svg.find('"', start) - start));
    std::vector<double> values{std::istream_iterator<double>(box), std::istream_iterator<double>()};
    if (values.size() != 4 || values[2] <= 0 || values[3] <= 0) throw std::runtime_error("The drawing has no view box");
    return values;
}

static void checkMapNumbersOnCanvas() {
    // Every map number stays whole on the canvas, on thumbnails too, and is never set larger
    // than the atom labels.  Cyclosporin A at 300 x 200 and taxol at 120 x 90 crowd every edge.
    for (const char* smiles : {
            "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)"
            "C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C",                        // cyclosporin A
            "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)"
            "(CO4)OC(=O)C)O)C)OC(=O)C"}) {                                                                   // taxol
        const auto graph = smsd::parseSMILES(smiles);
        std::map<int, int> mapping;
        for (int atom = 0; atom < graph.n; ++atom) mapping.emplace(atom, atom);
        for (const auto& [width, height] : {std::make_pair(300, 200), std::make_pair(120, 90)}) {
            smsd::DepictOptions options;
            options.width = width;
            options.height = height;
            const auto svg = smsd::depictWithMapping(graph, mapping, options);
            const auto box = viewBox(svg);
            double labelSize = 0;
            for (const auto& text : svgTexts(svg))
                if (text.cls == "sub" && !text.outline && !text.centred) labelSize = std::max(labelSize, text.size);
            for (const auto& text : svgTexts(svg)) {
                if (text.cls != "map") continue;
                const double halfWidth = 0.556 * text.size * text.text.size() / 2, halfHeight = 0.36 * text.size;
                if (text.x - halfWidth < box[0] || text.x + halfWidth > box[0] + box[2]
                        || text.y - halfHeight < box[1] || text.y + halfHeight > box[1] + box[3])
                    throw std::runtime_error(std::string("A map number runs off the canvas in ") + smiles);
                if (text.size > labelSize + 1e-6)
                    throw std::runtime_error(std::string("A map number outgrows the atom labels in ") + smiles);
            }
        }
    }
}

static void checkMapNumberSize() {
    // On a small fixed canvas the map numbers stay at least nine canvas pixels high, unless the
    // atom labels themselves are smaller.
    const auto first = smsd::parseSMILES("CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12");
    const auto second = smsd::parseSMILES("CCCc1nc(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(CC)CC4)ccc3OCC)nn12");
    std::map<int, int> mapping;
    for (int atom = 0; atom < first.n; ++atom) mapping.emplace(atom, atom);
    smsd::DepictOptions options;
    options.width = 600;
    options.height = 300;
    const auto svg = smsd::depictPair(first, second, mapping, options);
    const size_t start = svg.find("viewBox=\"") + 9;
    std::istringstream box(svg.substr(start, svg.find('"', start) - start));
    double x = 0, y = 0, viewWidth = 0;
    if (!(box >> x >> y >> viewWidth) || viewWidth <= 0) throw std::runtime_error("The pair has no view box");
    const double pixel = viewWidth / options.width;
    for (const auto& text : svgTexts(svg))
        if (text.cls == "map" && text.size / pixel < 9 - 1e-6)
            throw std::runtime_error("Map numbers must stay at least nine pixels high");
}

int main() {
    checkBondHighlightBands();
    checkHighlightLayer();
    checkLabelKnockouts();
    checkHighlightOverrides();
    checkHighlightFollowsBonds();
    checkKekuleDefault();
    checkMapNumberPlacement();
    checkPairWedges();
    checkRegularRingCoordinates();
    checkNoHiddenAtoms();
    checkSpiroRings();
    checkFoldedZChain();
    checkReadableStereoCentres();
    checkStoredDoubleBondGeometry();
    checkRingDoubleBondsInMolFiles();
    checkGridOrientation();
    checkTransparentHighlight();
    checkIsotopeLabels();
    checkWedgesOutsideRings();
    checkMapNumberSize();
    checkMapNumbersOnCanvas();
    checkFixedCanvasFit();
    checkKekuleBonds();
    checkMapNumbers();
    checkFontFamily();
    for (const auto& pair : {
            std::make_pair("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
            std::make_pair("F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I"),
            std::make_pair("[C@H](F)(Cl)Br", "[C@@H](F)(Cl)Br")}) {
        const auto first = smsd::depictSMILES(pair.first);
        const auto second = smsd::depictSMILES(pair.second);
        if (first == second) throw std::runtime_error("Enantiomers must have different wedge depictions");
    }
    for (const char* text : {"c1ccc2ccccc2c1", "c1ccc2[nH]ccc2c1", "Cc1ccc2ccccc2c1C",
                             "c1ccccc1CCc1ccccc1", "C1CCCCCCCCCCCCC1", "CC(=O)[O-].[Na+]",
                             "[Na+].[Cl-]", "CCCCCC.CCCCCC", "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
                             "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O"}) {
        const auto graph = smsd::parseSMILES(text);
        checkCoordinates(graph, smsd::generateCoords2D(graph), 1.5);
        std::vector<smsd::Point2D> depiction;
        for (const auto& point : smsd::layout2D(graph)) depiction.push_back({point.x, point.y});
        checkCoordinates(graph, depiction, 30.0);
        const auto scaffold = smsd::matchTemplate(graph, 1.5);
        if (!scaffold.empty()) checkCoordinates(graph, scaffold, 1.5);
    }
    std::vector<smsd::Point2D> line{{0, 0}, {std::cos(0.35), std::sin(0.35)}};
    smsd::transform::canonicalOrientation(smsd::parseSMILES("CC"), line);
    if (std::abs(line[1].y - line[0].y) > 1e-9)
        throw std::runtime_error("Canonical orientation must align the principal axis");
    std::cout << "Depiction regressions passed\n";
}

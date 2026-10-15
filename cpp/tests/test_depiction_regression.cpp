/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 */
#include "smsd/depict.hpp"
#include <iostream>
#include <stdexcept>

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

int main() {
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

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 */
#include <cmath>
#undef M_PI
#include "smsd/depict.hpp"
#include "smsd/mol_reader.hpp"

#include <chrono>
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static void requireSameGraph(const smsd::MolGraph& found, const smsd::MolGraph& expected) {
    require(found.n == expected.n && found.atomicNum == expected.atomicNum,
            "Unicode file changed the atom list");
    for (int i = 0; i < expected.n; ++i) {
        for (int j = 0; j < expected.n; ++j) {
            require(found.bondOrder(i, j) == expected.bondOrder(i, j),
                    "Unicode file changed a bond");
        }
    }
}

class TemporaryDirectory {
public:
    std::filesystem::path path;
    TemporaryDirectory() {
        const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path() / std::filesystem::u8path(
            std::string(u8"smsd-\u00e9-\u5316\u5b66-") + std::to_string(suffix));
        require(std::filesystem::create_directory(path), "temporary directory already exists");
    }
    ~TemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path, ignored);
    }
};

int main() {
    try {
        auto ring = smsd::MolGraph::Builder()
            .atomCount(6).atomicNumbers({6, 6, 6, 6, 6, 6})
            .setNeighbors({{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}})
            .setBondOrders({{1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}})
            .build();
        auto coordinates = smsd::layout2D(ring);
        require(coordinates.size() == 6, "ring layout lost an atom");
        for (int i = 0; i < ring.n; ++i) {
            const auto& point = coordinates[i];
            require(std::isfinite(point.x) && std::isfinite(point.y),
                    "ring layout produced a non-finite coordinate");
            for (int neighbor : ring.neighbors[i]) {
                require(smsd::len(point - coordinates[neighbor]) > 1.0,
                        "ring bond endpoints overlap");
            }
        }
        auto image = smsd::depictSVG(ring, coordinates);
        require(image.width > 0 && image.height > 0, "SVG dimensions are invalid");
        require(image.svg.find("<svg") != std::string::npos
                    && image.svg.find("</svg>") != std::string::npos,
                "ring depiction is not an SVG document");
        require(smsd::layout2D(smsd::MolGraph()).empty(),
                "empty molecule acquired coordinates");

        TemporaryDirectory temporary;
        auto molPath = temporary.path / std::filesystem::u8path(u8"cycle-\u00e9.mol");
        {
            std::ofstream output(molPath);
            output << smsd::writeMolBlock(ring);
            require(output.good(), "could not write the Unicode MOL fixture");
        }
        auto restored = smsd::readMolFile(molPath.u8string());
        requireSameGraph(restored, ring);

        auto sdfPath = temporary.path / std::filesystem::u8path(u8"cycles-\u5316\u5b66.sdf");
        smsd::writeSDF({ring, ring}, sdfPath.u8string());
        auto molecules = smsd::readSDF(sdfPath.u8string());
        require(molecules.size() == 2, "Unicode SDF file lost a record");
        for (const auto& molecule : molecules) {
            requireSameGraph(molecule, ring);
        }
        std::cout << "Standalone depiction and Unicode file checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

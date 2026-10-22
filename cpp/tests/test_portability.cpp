/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include <cmath>
#undef M_PI
#include "smsd/depict.hpp"
#include "smsd/mol_reader.hpp"
#include "test_support.hpp"
#include "smsd/bitops.hpp"

#include <chrono>
#include <iostream>
#include <stdexcept>

// The smsd_portable_bitops_tests build of this file forces the portable bit operations.
#if defined(SMSD_PORTABLE_BITOPS) \
    && (defined(SMSD_BITOPS_BUILTINS) || defined(SMSD_BITOPS_MSVC_X64) || defined(SMSD_BITOPS_MSVC_ARM64))
#error "SMSD_PORTABLE_BITOPS did not select the portable bit operations"
#endif

using smsd_test::require;

// popcount64 and ctz64, whichever code the target selects, against a bit-by-bit count.
static void requireBitOperations() {
    uint64_t word = 0x9E3779B97F4A7C15ULL;
    for (int round = 0; round < 4096; ++round) {
        word ^= word << 13; word ^= word >> 7; word ^= word << 17;
        for (uint64_t value : {word, word & (word >> 3), word >> (round & 63), uint64_t(1) << (round & 63)}) {
            int bits = 0, lowest = 64;
            for (int position = 63; position >= 0; --position)
                if ((value >> position) & 1) { ++bits; lowest = position; }
            require(smsd::popcount64(value) == bits && smsd::popcount64_portable(value) == bits,
                    "popcount64 miscounted a word");
            require(smsd::ctz64(value) == lowest && smsd::ctz64_portable(value) == lowest,
                    "ctz64 found the wrong lowest bit");
        }
    }
    require(smsd::popcount64(0) == 0 && smsd::ctz64(0) == 64 && smsd::ctz64_portable(0) == 64,
            "an empty word has bits");
    require(smsd::popcount64(~uint64_t(0)) == 64 && smsd::ctz64(uint64_t(1) << 63) == 63,
            "a full word lost bits");
    const std::vector<uint64_t> words = smsd::bitWords({0, 63, 64, 130}, 131);
    require(smsd::bitPositions(words) == std::vector<int>({0, 63, 64, 130}), "set bit positions changed");
    // A search that rests on the bit operations: naphthalene in anthracene.
    const auto common = smsd::MCS::find(smsd::parseSMILES("c1ccc2ccccc2c1")).match(smsd::parseSMILES("c1ccc2cc3ccccc3cc2c1"));
    require(!common.matches.empty() && common.matches[0].atomMapping.size() == 10, "naphthalene is not common to anthracene");
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

// The portable build of this file runs beside the other one, so their directories are named apart.
#if defined(SMSD_PORTABLE_BITOPS)
static const char* const BUILD_NAME = "portable-";
#else
static const char* const BUILD_NAME = "";
#endif

class TemporaryDirectory {
public:
    std::filesystem::path path;
    TemporaryDirectory() {
        const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path() / std::filesystem::u8path(
            std::string(u8"smsd-\u00e9-\u5316\u5b66-") + BUILD_NAME + std::to_string(suffix));
        require(std::filesystem::create_directory(path), "temporary directory already exists");
    }
    ~TemporaryDirectory() {
        std::error_code ignored;
        std::filesystem::remove_all(path, ignored);
    }
};

int main() {
    try {
        requireBitOperations();
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

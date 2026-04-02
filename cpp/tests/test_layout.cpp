/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Test suite for force-directed layout, stress majorisation, and template
 * matching (layout.hpp).
 */

#include "smsd/smsd.hpp"
#include "smsd/smiles_parser.hpp"
#include "smsd/layout.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// ============================================================================
// Test harness
// ============================================================================

static int g_pass = 0, g_fail = 0;

#define RUN_TEST(name) do { \
    std::cout << "  [" #name "] "; \
    try { test_##name(); g_pass++; std::cout << "PASS\n"; } \
    catch (const std::exception& e) { g_fail++; std::cout << "FAIL: " << e.what() << "\n"; } \
    catch (...) { g_fail++; std::cout << "FAIL (unknown)\n"; } \
} while(0)

#define TEST_ASSERT(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL: %s  [%s:%d]\n", msg, __FILE__, __LINE__); \
        throw std::runtime_error(msg); \
    } \
} while(0)

// ============================================================================
// Helper: create initial linear layout
// ============================================================================

static std::vector<smsd::Point2D> linearLayout(int n) {
    std::vector<smsd::Point2D> coords(n);
    for (int i = 0; i < n; i++) {
        coords[i].x = static_cast<double>(i) * 1.0;
        coords[i].y = 0.0;
    }
    return coords;
}

// ============================================================================
// Force-directed layout tests
// ============================================================================

static void test_force_directed_benzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = linearLayout(g.n);
    double stress = smsd::forceDirectedLayout(g, coords, 500, 1.5);
    TEST_ASSERT(stress >= 0.0, "Stress must be non-negative");
    // After optimisation, bonded atoms should be roughly targetBondLength apart
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue;
            double dx = coords[j].x - coords[i].x;
            double dy = coords[j].y - coords[i].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            TEST_ASSERT(dist > 0.5 && dist < 3.0,
                        "Bond length should be in reasonable range after layout");
        }
    }
}

static void test_force_directed_naphthalene() {
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    auto coords = linearLayout(g.n);
    double stress = smsd::forceDirectedLayout(g, coords);
    TEST_ASSERT(stress >= 0.0, "Stress must be non-negative for naphthalene");
    // Check no atoms are at the same position
    for (int i = 0; i < g.n; i++) {
        for (int j = i + 1; j < g.n; j++) {
            double dx = coords[j].x - coords[i].x;
            double dy = coords[j].y - coords[i].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            TEST_ASSERT(dist > 0.01, "Atoms should not overlap after layout");
        }
    }
}

static void test_force_directed_small() {
    // Ethane: 2 atoms
    auto g = smsd::parseSMILES("CC");
    auto coords = linearLayout(g.n);
    double stress = smsd::forceDirectedLayout(g, coords);
    TEST_ASSERT(stress >= 0.0, "Stress must be non-negative for ethane");
}

static void test_force_directed_single_atom() {
    auto g = smsd::parseSMILES("C");
    std::vector<smsd::Point2D> coords = {{0.0, 0.0}};
    double stress = smsd::forceDirectedLayout(g, coords);
    TEST_ASSERT(stress == 0.0, "Single atom should have zero stress");
}

static void test_force_directed_short_coords_resized() {
    auto g = smsd::parseSMILES("CCC");
    std::vector<smsd::Point2D> coords = {{0.0, 0.0}};
    double stress = smsd::forceDirectedLayout(g, coords, 50, 1.5);
    TEST_ASSERT(stress >= 0.0, "Short coordinate vector should be accepted");
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n,
                "Force-directed layout should resize short coordinate vectors");
}

// ============================================================================
// Stress majorisation tests
// ============================================================================

static void test_smacof_benzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = linearLayout(g.n);
    double stress = smsd::stressMajorisation(g, coords, 300, 1.5);
    TEST_ASSERT(stress >= 0.0, "SMACOF stress must be non-negative");
    // Bonded atoms should have reasonable distance
    for (int i = 0; i < g.n; i++) {
        for (int j : g.neighbors[i]) {
            if (j <= i) continue;
            double dx = coords[j].x - coords[i].x;
            double dy = coords[j].y - coords[i].y;
            double dist = std::sqrt(dx * dx + dy * dy);
            TEST_ASSERT(dist > 0.3, "Bonded distance should be positive after SMACOF");
        }
    }
}

static void test_smacof_naphthalene() {
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    auto coords = linearLayout(g.n);
    double stress = smsd::stressMajorisation(g, coords);
    TEST_ASSERT(stress >= 0.0, "SMACOF stress must be non-negative for naphthalene");
}

static void test_smacof_convergence() {
    // SMACOF should reduce stress relative to initial linear layout
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    auto coords = linearLayout(g.n);
    // Compute initial stress manually
    auto graphDist = smsd::detail_layout::bfsDistanceMatrix(g);
    double initialStress = 0.0;
    for (int i = 0; i < g.n; i++) {
        for (int j = i + 1; j < g.n; j++) {
            double D = graphDist[i][j] * 1.5;
            double dx = coords[i].x - coords[j].x;
            double dy = coords[i].y - coords[j].y;
            double dij = std::sqrt(dx * dx + dy * dy);
            double w = (D > 1e-10) ? 1.0 / (D * D) : 0.0;
            double delta = dij - D;
            initialStress += w * delta * delta;
        }
    }
    double finalStress = smsd::stressMajorisation(g, coords, 300, 1.5);
    TEST_ASSERT(finalStress <= initialStress + 1e-6,
                "SMACOF should not increase stress");
}

static void test_smacof_single_atom() {
    auto g = smsd::parseSMILES("C");
    std::vector<smsd::Point2D> coords = {{0.0, 0.0}};
    double stress = smsd::stressMajorisation(g, coords);
    TEST_ASSERT(stress == 0.0, "Single atom should have zero SMACOF stress");
}

static void test_smacof_short_coords_resized() {
    auto g = smsd::parseSMILES("c1ccccc1");
    std::vector<smsd::Point2D> coords = {{0.0, 0.0}, {1.0, 0.0}};
    double stress = smsd::stressMajorisation(g, coords, 100, 1.5);
    TEST_ASSERT(stress >= 0.0, "SMACOF should accept short coordinate vectors");
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n,
                "SMACOF should resize short coordinate vectors");
}

// ============================================================================
// Template matching tests
// ============================================================================

static void test_template_benzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = smsd::matchTemplate(g);
    TEST_ASSERT(!coords.empty(), "Benzene should match a template");
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n,
                "Template coords should have one point per atom");
}

static void test_template_no_match() {
    // A chain molecule should not match any template
    auto g = smsd::parseSMILES("CCCCCCCC");
    auto coords = smsd::matchTemplate(g);
    TEST_ASSERT(coords.empty(), "Linear alkane should not match any template");
}

static void test_template_scaling() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords1 = smsd::matchTemplate(g, 1.0);
    auto coords2 = smsd::matchTemplate(g, 2.0);
    TEST_ASSERT(!coords1.empty() && !coords2.empty(),
                "Both scale factors should match benzene");
    // coords2 should be ~2x the coords1
    if (!coords1.empty() && !coords2.empty()) {
        double ratio = std::abs(coords2[0].x) /
                       (std::abs(coords1[0].x) > 1e-10 ? std::abs(coords1[0].x) : 1.0);
        // Allow 10% tolerance
        TEST_ASSERT(std::abs(ratio - 2.0) < 0.3 || std::abs(coords1[0].x) < 1e-10,
                    "Template coords should scale with targetBondLength");
    }
}

// ============================================================================
// Integration test: force-directed + crossing reduction pipeline
// ============================================================================

static void test_pipeline_morphine_like() {
    // A polycyclic system: use force-directed first, then crossing reduction
    auto g = smsd::parseSMILES("c1ccc2c(c1)c1ccccc1c1ccccc21"); // fluorene-like
    auto coords = linearLayout(g.n);
    smsd::forceDirectedLayout(g, coords, 200, 1.5);
    int crossings = smsd::reduceCrossings(g, coords, 500);
    // We just verify the pipeline runs without error
    TEST_ASSERT(crossings >= 0, "Crossing count must be non-negative");
}

static void test_reduce_crossings_short_coords_resized() {
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    std::vector<smsd::Point2D> coords = {{0.0, 0.0}, {1.0, 0.0}};
    int crossings = smsd::reduceCrossings(g, coords, 100);
    TEST_ASSERT(crossings >= 0, "Crossing reduction should accept short coordinate vectors");
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n,
                "Crossing reduction should resize short coordinate vectors");
}

// ============================================================================
// main
// ============================================================================

int main() {
    std::cout << "\n=== Layout engine test suite ===\n\n";

    std::cout << "Force-directed layout:\n";
    RUN_TEST(force_directed_benzene);
    RUN_TEST(force_directed_naphthalene);
    RUN_TEST(force_directed_small);
    RUN_TEST(force_directed_single_atom);
    RUN_TEST(force_directed_short_coords_resized);

    std::cout << "\nStress majorisation (SMACOF):\n";
    RUN_TEST(smacof_benzene);
    RUN_TEST(smacof_naphthalene);
    RUN_TEST(smacof_convergence);
    RUN_TEST(smacof_single_atom);
    RUN_TEST(smacof_short_coords_resized);

    std::cout << "\nTemplate matching:\n";
    RUN_TEST(template_benzene);
    RUN_TEST(template_no_match);
    RUN_TEST(template_scaling);

    std::cout << "\nIntegration:\n";
    RUN_TEST(pipeline_morphine_like);
    RUN_TEST(reduce_crossings_short_coords_resized);

    std::cout << "\n=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n\n";
    return g_fail > 0 ? 1 : 0;
}

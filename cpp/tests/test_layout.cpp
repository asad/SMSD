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
// v6.11.0 — generateCoords2D comprehensive pipeline tests
// ============================================================================

static void test_generate2d_benzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Should generate n coordinates");
    int crossings = smsd::detail_layout::countCrossings(g, coords);
    TEST_ASSERT(crossings == 0, "Benzene layout should have zero crossings");
}

static void test_generate2d_naphthalene() {
    auto g = smsd::parseSMILES("c1ccc2ccccc2c1");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Should generate n coordinates");
    int crossings = smsd::detail_layout::countCrossings(g, coords);
    TEST_ASSERT(crossings <= 2, "Naphthalene should have minimal crossings");
}

static void test_generate2d_aspirin() {
    auto g = smsd::parseSMILES("CC(=O)Oc1ccccc1C(=O)O");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Aspirin: n coords");
    int crossings = smsd::detail_layout::countCrossings(g, coords);
    TEST_ASSERT(crossings == 0, "Aspirin should have zero crossings");
}

static void test_generate2d_caffeine() {
    auto g = smsd::parseSMILES("Cn1c(=O)c2c(ncn2C)n(C)c1=O");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Caffeine: n coords");
}

static void test_generate2d_single_atom() {
    auto g = smsd::parseSMILES("C");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(coords.size() == 1, "Single atom: 1 coord");
}

static void test_generate2d_ethanol() {
    auto g = smsd::parseSMILES("CCO");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Ethanol: n coords");
    // Chain: atoms should not overlap
    for (int i = 0; i < g.n; i++) {
        for (int j = i + 1; j < g.n; j++) {
            double dx = coords[i].x - coords[j].x;
            double dy = coords[i].y - coords[j].y;
            TEST_ASSERT(dx * dx + dy * dy > 0.01, "Atoms should not overlap");
        }
    }
}

static void test_generate2d_indole() {
    auto g = smsd::parseSMILES("c1ccc2[nH]ccc2c1");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Indole: n coords");
}

static void test_generate2d_peg_chain() {
    auto g = smsd::parseSMILES("COCCOCCOCCOCCO");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "PEG chain: n coords");
    int crossings = smsd::detail_layout::countCrossings(g, coords);
    TEST_ASSERT(crossings == 0, "PEG chain should have zero crossings");
}

static void test_generate2d_salt() {
    // Disconnected: sodium acetate
    auto g = smsd::parseSMILES("CC(=O)[O-].[Na+]");
    auto coords = smsd::generateCoords2D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Salt: n coords");
}

// ============================================================================
// v6.11.0 — 3D coordinate generation tests
// ============================================================================

static void test_generate3d_benzene() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = smsd::generateCoords3D(g);
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "3D benzene: n coords");
    // All atoms should be non-degenerate (not all at origin)
    double sumDist2 = 0.0;
    for (int i = 0; i < g.n; i++) {
        sumDist2 += coords[i].x * coords[i].x + coords[i].y * coords[i].y + coords[i].z * coords[i].z;
    }
    TEST_ASSERT(sumDist2 > 0.1, "3D coords should not be degenerate");
}

static void test_generate3d_ethane() {
    auto g = smsd::parseSMILES("CC");
    auto coords = smsd::generateCoords3D(g);
    TEST_ASSERT(coords.size() == 2, "Ethane: 2 coords");
    double dx = coords[1].x - coords[0].x;
    double dy = coords[1].y - coords[0].y;
    double dz = coords[1].z - coords[0].z;
    double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
    TEST_ASSERT(dist > 0.5, "3D bond length should be reasonable");
}

// ============================================================================
// v6.11.0 — Coordinate transform tests
// ============================================================================

static void test_translate2d() {
    std::vector<smsd::Point2D> coords = {{0, 0}, {1, 0}, {0, 1}};
    smsd::transform::translate2D(coords, 5.0, 3.0);
    TEST_ASSERT(std::abs(coords[0].x - 5.0) < 1e-10, "Translate X");
    TEST_ASSERT(std::abs(coords[0].y - 3.0) < 1e-10, "Translate Y");
    TEST_ASSERT(std::abs(coords[1].x - 6.0) < 1e-10, "Translate X2");
}

static void test_rotate2d() {
    std::vector<smsd::Point2D> coords = {{1, 0}, {-1, 0}};
    smsd::transform::rotate2D(coords, M_PI / 2.0);
    // After 90-degree rotation about centroid (0,0): (1,0) -> (0,1)
    TEST_ASSERT(std::abs(coords[0].x - 0.0) < 1e-6, "Rotate 90: x~0");
    TEST_ASSERT(std::abs(coords[0].y - 1.0) < 1e-6, "Rotate 90: y~1");
}

static void test_scale2d() {
    std::vector<smsd::Point2D> coords = {{0, 0}, {2, 0}};
    smsd::transform::scale2D(coords, 2.0);
    double dx = coords[1].x - coords[0].x;
    TEST_ASSERT(std::abs(dx - 4.0) < 1e-6, "Scale 2x: distance doubles");
}

static void test_mirror_x() {
    std::vector<smsd::Point2D> coords = {{0, 1}, {0, -1}};
    smsd::transform::mirrorX(coords);
    // Centroid y=0, mirror: (0,1) -> (0,-1)
    TEST_ASSERT(std::abs(coords[0].y - (-1.0)) < 1e-6, "Mirror X: y flipped");
}

static void test_mirror_y() {
    std::vector<smsd::Point2D> coords = {{1, 0}, {-1, 0}};
    smsd::transform::mirrorY(coords);
    TEST_ASSERT(std::abs(coords[0].x - (-1.0)) < 1e-6, "Mirror Y: x flipped");
}

static void test_center2d() {
    std::vector<smsd::Point2D> coords = {{5, 5}, {7, 5}};
    smsd::transform::center2D(coords);
    double cx = (coords[0].x + coords[1].x) / 2.0;
    double cy = (coords[0].y + coords[1].y) / 2.0;
    TEST_ASSERT(std::abs(cx) < 1e-10, "Centered X");
    TEST_ASSERT(std::abs(cy) < 1e-10, "Centered Y");
}

static void test_align2d() {
    std::vector<smsd::Point2D> a = {{0, 0}, {1, 0}, {0, 1}};
    std::vector<smsd::Point2D> ref = {{5, 5}, {6, 5}, {5, 6}};
    double rmsd = smsd::transform::align2D(a, ref);
    TEST_ASSERT(rmsd < 0.1, "Align RMSD should be near zero for congruent triangles");
}

static void test_normalise_bond_length() {
    auto g = smsd::parseSMILES("CC");
    std::vector<smsd::Point2D> coords = {{0, 0}, {3, 0}};
    smsd::transform::normaliseBondLength(g, coords, 1.5);
    double dx = coords[1].x - coords[0].x;
    TEST_ASSERT(std::abs(dx - 1.5) < 0.1, "Bond length normalised to 1.5");
}

static void test_canonical_orientation() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = smsd::generateCoords2D(g);
    smsd::transform::canonicalOrientation(g, coords);
    // Should not crash and should preserve atom count
    TEST_ASSERT(static_cast<int>(coords.size()) == g.n, "Canonical orientation preserves size");
}

static void test_bounding_box() {
    std::vector<smsd::Point2D> coords = {{-1, -2}, {3, 4}, {0, 0}};
    auto bb = smsd::transform::boundingBox2D(coords);
    TEST_ASSERT(std::abs(bb[0] - (-1.0)) < 1e-10, "BB minX");
    TEST_ASSERT(std::abs(bb[1] - (-2.0)) < 1e-10, "BB minY");
    TEST_ASSERT(std::abs(bb[2] - 3.0) < 1e-10, "BB maxX");
    TEST_ASSERT(std::abs(bb[3] - 4.0) < 1e-10, "BB maxY");
}

static void test_project_lift_roundtrip() {
    std::vector<smsd::Point3D> pts3d = {{1, 2, 3}, {4, 5, 6}};
    auto pts2d = smsd::transform::projectTo2D(pts3d);
    TEST_ASSERT(pts2d.size() == 2, "Project: correct count");
    TEST_ASSERT(std::abs(pts2d[0].x - 1.0) < 1e-10, "Project: x preserved");
    auto back = smsd::transform::liftTo3D(pts2d);
    TEST_ASSERT(std::abs(back[0].z) < 1e-10, "Lift: z=0");
}

// ============================================================================
// v6.11.0 — Overlap resolution tests
// ============================================================================

static void test_resolve_overlaps_separate() {
    std::vector<smsd::Point2D> coords = {{0, 0}, {5, 5}};
    int remaining = smsd::resolveOverlaps(coords, 0.5, 10);
    TEST_ASSERT(remaining == 0, "Non-overlapping: 0 remaining");
}

static void test_resolve_overlaps_touching() {
    std::vector<smsd::Point2D> coords = {{0, 0}, {0.1, 0.1}};
    int remaining = smsd::resolveOverlaps(coords, 0.5, 50);
    TEST_ASSERT(remaining == 0, "Overlapping pair resolved");
    double dx = coords[1].x - coords[0].x;
    double dy = coords[1].y - coords[0].y;
    double dist = std::sqrt(dx * dx + dy * dy);
    TEST_ASSERT(dist >= 0.49, "After resolution, distance >= threshold");
}

// ============================================================================
// v6.11.0 — Layout quality metric tests
// ============================================================================

static void test_layout_quality_perfect() {
    auto g = smsd::parseSMILES("c1ccccc1");
    auto coords = smsd::generateCoords2D(g);
    double q = smsd::layoutQuality(g, coords);
    // Benzene with template match should be nearly perfect
    TEST_ASSERT(q < 1.0, "Benzene quality should be very good");
}

static void test_layout_quality_bad() {
    auto g = smsd::parseSMILES("c1ccccc1");
    // All atoms at same point — terrible layout
    std::vector<smsd::Point2D> coords(6, {0, 0});
    double q = smsd::layoutQuality(g, coords);
    TEST_ASSERT(q > 10.0, "Degenerate layout should have high quality score");
}

// ============================================================================
// v6.11.0 — Extended template library tests
// ============================================================================

static void test_extended_templates_exist() {
    auto& templates = smsd::detail_layout::getExtendedTemplates();
    TEST_ASSERT(templates.size() >= 30, "Should have 30+ extended templates");
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

    std::cout << "\ngenerate_coords_2d pipeline (v6.11.0):\n";
    RUN_TEST(generate2d_benzene);
    RUN_TEST(generate2d_naphthalene);
    RUN_TEST(generate2d_aspirin);
    RUN_TEST(generate2d_caffeine);
    RUN_TEST(generate2d_single_atom);
    RUN_TEST(generate2d_ethanol);
    RUN_TEST(generate2d_indole);
    RUN_TEST(generate2d_peg_chain);
    RUN_TEST(generate2d_salt);

    std::cout << "\ngenerate_coords_3d (v6.11.0):\n";
    RUN_TEST(generate3d_benzene);
    RUN_TEST(generate3d_ethane);

    std::cout << "\nCoordinate transforms (v6.11.0):\n";
    RUN_TEST(translate2d);
    RUN_TEST(rotate2d);
    RUN_TEST(scale2d);
    RUN_TEST(mirror_x);
    RUN_TEST(mirror_y);
    RUN_TEST(center2d);
    RUN_TEST(align2d);
    RUN_TEST(normalise_bond_length);
    RUN_TEST(canonical_orientation);
    RUN_TEST(bounding_box);
    RUN_TEST(project_lift_roundtrip);

    std::cout << "\nOverlap resolution (v6.11.0):\n";
    RUN_TEST(resolve_overlaps_separate);
    RUN_TEST(resolve_overlaps_touching);

    std::cout << "\nLayout quality (v6.11.0):\n";
    RUN_TEST(layout_quality_perfect);
    RUN_TEST(layout_quality_bad);

    std::cout << "\nExtended templates (v6.11.0):\n";
    RUN_TEST(extended_templates_exist);

    std::cout << "\n=== Results: " << g_pass << " passed, " << g_fail << " failed ===\n\n";
    return g_fail > 0 ? 1 : 0;
}

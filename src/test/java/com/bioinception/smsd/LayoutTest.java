/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import org.junit.jupiter.api.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Tests for force-directed layout, stress majorisation, and template matching.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Layout Engine Tests")
public class LayoutTest extends TestBase {

  private static MolGraph graph(String smi) throws Exception {
    return new MolGraph(mol(smi));
  }

  private static MolGraph.Point2D[] linearLayout(int n) {
    MolGraph.Point2D[] coords = new MolGraph.Point2D[n];
    for (int i = 0; i < n; i++)
      coords[i] = new MolGraph.Point2D(i * 1.0, 0.0);
    return coords;
  }

  // ======================================================================
  // Force-directed layout
  // ======================================================================

  @Nested
  @DisplayName("Force-Directed Layout")
  class ForceDirected {

    @Test
    @DisplayName("Benzene: stress is non-negative")
    void testBenzeneStress() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.forceDirectedLayout(g, coords);
      assertTrue(stress >= 0.0, "Stress must be non-negative");
    }

    @Test
    @DisplayName("Benzene: atoms are separated after layout")
    void testBenzeneAtomsNotOverlapping() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      MolGraph.forceDirectedLayout(g, coords);
      for (int i = 0; i < g.atomCount(); i++) {
        for (int j = i + 1; j < g.atomCount(); j++) {
          double dx = coords[i].x - coords[j].x;
          double dy = coords[i].y - coords[j].y;
          double dist = Math.sqrt(dx * dx + dy * dy);
          assertTrue(dist > 0.01,
              "Atoms " + i + " and " + j + " should not overlap");
        }
      }
    }

    @Test
    @DisplayName("Naphthalene: runs without error")
    void testNaphthalene() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.forceDirectedLayout(g, coords);
      assertTrue(stress >= 0.0);
      assertEquals(g.atomCount(), coords.length);
    }

    @Test
    @DisplayName("Ethane: two atoms layout")
    void testEthane() throws Exception {
      MolGraph g = graph("CC");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.forceDirectedLayout(g, coords);
      assertTrue(stress >= 0.0);
    }

    @Test
    @DisplayName("Custom parameters accepted")
    void testCustomParams() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.forceDirectedLayout(g, coords, 100, 2.0);
      assertTrue(stress >= 0.0);
    }

    @Test
    @DisplayName("Short coordinate array is padded safely")
    void testShortCoordsPadded() throws Exception {
      MolGraph g = graph("CCC");
      MolGraph.Point2D[] coords = {new MolGraph.Point2D(0.0, 0.0)};
      double stress = MolGraph.forceDirectedLayout(g, coords, 50, 1.5);
      assertTrue(stress >= 0.0);
      assertNotNull(coords[0]);
    }
  }

  // ======================================================================
  // Stress majorisation (SMACOF)
  // ======================================================================

  @Nested
  @DisplayName("Stress Majorisation (SMACOF)")
  class StressMaj {

    @Test
    @DisplayName("Benzene: stress is non-negative")
    void testBenzeneStress() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.stressMajorisation(g, coords);
      assertTrue(stress >= 0.0, "SMACOF stress must be non-negative");
    }

    @Test
    @DisplayName("Naphthalene: runs without error")
    void testNaphthalene() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.stressMajorisation(g, coords);
      assertTrue(stress >= 0.0);
    }

    @Test
    @DisplayName("Atoms not collapsed after layout")
    void testAtomsNotCollapsed() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      MolGraph.stressMajorisation(g, coords);
      for (int i = 0; i < g.atomCount(); i++) {
        for (int j = i + 1; j < g.atomCount(); j++) {
          double dx = coords[i].x - coords[j].x;
          double dy = coords[i].y - coords[j].y;
          double dist = Math.sqrt(dx * dx + dy * dy);
          assertTrue(dist > 0.01,
              "Atoms " + i + " and " + j + " should not collapse");
        }
      }
    }

    @Test
    @DisplayName("Custom parameters accepted")
    void testCustomParams() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      double stress = MolGraph.stressMajorisation(g, coords, 50, 2.0);
      assertTrue(stress >= 0.0);
    }

    @Test
    @DisplayName("Short coordinate array is padded safely")
    void testShortCoordsPadded() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = {
          new MolGraph.Point2D(0.0, 0.0),
          new MolGraph.Point2D(1.0, 0.0)
      };
      double stress = MolGraph.stressMajorisation(g, coords, 50, 1.5);
      assertTrue(stress >= 0.0);
      assertNotNull(coords[0]);
    }
  }

  // ======================================================================
  // Template matching
  // ======================================================================

  @Nested
  @DisplayName("Template Matching")
  class Templates {

    @Test
    @DisplayName("Benzene matches template")
    void testBenzeneMatch() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = MolGraph.matchTemplate(g);
      assertNotNull(coords, "Benzene should match a template");
      assertEquals(6, coords.length);
    }

    @Test
    @DisplayName("Chain does not match any template")
    void testChainNoMatch() throws Exception {
      MolGraph g = graph("CCCCCCCC");
      MolGraph.Point2D[] coords = MolGraph.matchTemplate(g);
      assertNull(coords, "Linear alkane should not match any template");
    }

    @Test
    @DisplayName("Scaling changes coordinates")
    void testScaling() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] c1 = MolGraph.matchTemplate(g, 1.0);
      MolGraph.Point2D[] c2 = MolGraph.matchTemplate(g, 2.0);
      assertNotNull(c1);
      assertNotNull(c2);
      // Find a non-zero coordinate to test ratio
      for (int k = 0; k < 6; k++) {
        if (Math.abs(c1[k].x) > 0.01) {
          double ratio = Math.abs(c2[k].x) / Math.abs(c1[k].x);
          assertEquals(2.0, ratio, 0.3, "Template should scale with bond length");
          break;
        }
      }
    }

    @Test
    @DisplayName("Template coordinates are finite")
    void testCoordsFinite() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = MolGraph.matchTemplate(g);
      assertNotNull(coords);
      for (MolGraph.Point2D p : coords) {
        assertTrue(Double.isFinite(p.x), "x must be finite");
        assertTrue(Double.isFinite(p.y), "y must be finite");
      }
    }
  }

  // ======================================================================
  // Integration: layout pipeline
  // ======================================================================

  @Nested
  @DisplayName("Layout Pipeline Integration")
  class Pipeline {

    @Test
    @DisplayName("Force-directed + crossing reduction")
    void testForceThenCrossingReduction() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      MolGraph.forceDirectedLayout(g, coords, 200, 1.5);
      int crossings = MolGraph.reduceCrossings(g, coords, 500);
      assertTrue(crossings >= 0, "Crossing count must be non-negative");
    }

    @Test
    @DisplayName("SMACOF + crossing reduction")
    void testSmacofThenCrossingReduction() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      MolGraph.stressMajorisation(g, coords, 100, 1.5);
      int crossings = MolGraph.reduceCrossings(g, coords, 500);
      assertTrue(crossings >= 0, "Crossing count must be non-negative");
    }

    @Test
    @DisplayName("Template match or fallback to force-directed")
    void testTemplateOrFallback() throws Exception {
      MolGraph g = graph("c1ccccc1");
      MolGraph.Point2D[] coords = MolGraph.matchTemplate(g);
      if (coords == null) {
        coords = linearLayout(g.atomCount());
        MolGraph.forceDirectedLayout(g, coords);
      }
      assertEquals(g.atomCount(), coords.length);
    }

    @Test
    @DisplayName("Fused polycyclic: individual ring flipping reduces crossings")
    void testFusedRingFlipping() throws Exception {
      // Phenanthrene: three linearly fused rings — a good test for phase-2
      // individual ring flipping since system-level flips cannot resolve
      // crossings between the inner and outer rings.
      MolGraph g = graph("c1ccc2c(c1)cc1ccccc1c2");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      int before = MolGraph.reduceCrossings(g, coords, 0); // no iterations, just count
      int after = MolGraph.reduceCrossings(g, coords, 2000);
      assertTrue(after <= before,
          "Crossing count should not increase (before=" + before + ", after=" + after + ")");
      assertTrue(after >= 0, "Crossing count must be non-negative");
    }

    @Test
    @DisplayName("Steroid skeleton: phase-2 handles 4 fused rings")
    void testSteroidCrossingReduction() throws Exception {
      // Sterane: 4 fused 6-5-6-6 ring skeleton — bridges and fusions test
      // phase 2 individual ring flipping thoroughly.
      MolGraph g = graph("C1CCC2C(C1)CCC1C2CCC2CCCCC21");
      MolGraph.Point2D[] coords = linearLayout(g.atomCount());
      int crossings = MolGraph.reduceCrossings(g, coords, 2000);
      assertTrue(crossings >= 0, "Crossing count must be non-negative");
    }

    @Test
    @DisplayName("Crossing reduction accepts short coordinate arrays")
    void testCrossingReductionShortCoords() throws Exception {
      MolGraph g = graph("c1ccc2ccccc2c1");
      MolGraph.Point2D[] coords = {
          new MolGraph.Point2D(0.0, 0.0),
          new MolGraph.Point2D(1.0, 0.0)
      };
      int crossings = MolGraph.reduceCrossings(g, coords, 200);
      assertTrue(crossings >= 0, "Crossing count must be non-negative");
      assertNotNull(coords[0]);
    }
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;
import java.util.Random;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.api.Test;

class GraphLayoutRegressionTest {
  @ParameterizedTest
  @ValueSource(strings = {"c1ccccc1", "C1CCCCC1", "c1ccc2ccccc2c1", "c1cccc2ccccc12",
      "c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1", "c1ccccc1-c2ccccc2"})
  void mapsScaffoldCoordinatesToBondTopology(String smiles) throws Exception {
    MolGraph graph = rawGraph(smiles);
    MolGraph.Point2D[] coordinates = MolGraph.matchTemplate(graph);
    assertNotNull(coordinates);
    assertOrdinaryBondGeometry(graph, coordinates);
  }

  @Test
  void rejectsDifferentTopologyWithTheSameRingSizes() throws Exception {
    assertNull(MolGraph.matchTemplate(rawGraph("CC1CC2CCCC(C1)C2")));
  }

  @Test
  void seedsDegenerateFusedCoordinatesWithMappedTemplate() throws Exception {
    MolGraph graph = rawGraph("c1ccc2ccccc2c1");
    MolGraph.Point2D[] coordinates = new MolGraph.Point2D[graph.n];
    for (int atom = 0; atom < graph.n; atom++) coordinates[atom] = new MolGraph.Point2D(0, 0);
    MolGraph.forceDirectedLayout(graph, coordinates, 0, 1.5);
    assertOrdinaryBondGeometry(graph, coordinates);
  }

  @ParameterizedTest
  @ValueSource(strings = {"CCCCC1CCCCC1", "c1ccccc1Cc1ccccc1"})
  void writesOptimisedCoordinatesBeforeReturningCrossingCount(String smiles) throws Exception {
    MolGraph graph = rawGraph(smiles);
    for (int seed = 0; seed < 8; seed++) {
      Random random = new Random(seed);
      MolGraph.Point2D[] coordinates = new MolGraph.Point2D[graph.n];
      for (int atom = 0; atom < graph.n; atom++)
        coordinates[atom] = new MolGraph.Point2D(random.nextDouble(), random.nextDouble());
      int returned = MolGraph.reduceCrossings(graph, coordinates, 1000);
      assertEquals(crossings(graph, coordinates), returned, "returned count must describe output coordinates");
    }
  }

  @ParameterizedTest
  @ValueSource(ints = {21, 22, 33})
  void perceivesMacrocyclesWithoutAnArbitrarySizeLimit(int atoms) throws Exception {
    StringBuilder smiles = new StringBuilder("C1");
    for (int atom = 1; atom < atoms; atom++) smiles.append('C');
    smiles.append('1');
    MolGraph graph = rawGraph(smiles.toString());
    assertEquals(1, graph.computeRings().length);
    assertEquals(atoms, graph.computeRings()[0].length);
    assertEquals(1, MolGraph.computeSSSR(graph).size());
    assertEquals(atoms, MolGraph.layoutSSSR(graph).get(0).size());
  }

  private static void assertOrdinaryBondGeometry(MolGraph graph, MolGraph.Point2D[] coordinates) {
    for (int atom = 0; atom < graph.n; atom++) {
      assertTrue(Double.isFinite(coordinates[atom].x) && Double.isFinite(coordinates[atom].y));
      for (int neighbor : graph.neighbors[atom]) if (neighbor > atom)
        assertEquals(1.5, Math.hypot(coordinates[atom].x - coordinates[neighbor].x,
            coordinates[atom].y - coordinates[neighbor].y), 1e-8, "each template bond has the requested length");
      for (int other = atom + 1; other < graph.n; other++)
        assertTrue(Math.hypot(coordinates[atom].x - coordinates[other].x,
            coordinates[atom].y - coordinates[other].y) > 1e-6, "atoms must not overlap");
    }
    assertEquals(0, crossings(graph, coordinates));
  }

  private static int crossings(MolGraph graph, MolGraph.Point2D[] coordinates) {
    int count = 0;
    for (int first = 0; first < graph.n; first++)
      for (int second : graph.neighbors[first]) if (second > first)
        for (int third = first + 1; third < graph.n; third++)
          for (int fourth : graph.neighbors[third]) if (fourth > third
              && first != third && first != fourth && second != third && second != fourth) {
            double a = orientation(coordinates[first], coordinates[second], coordinates[third]);
            double b = orientation(coordinates[first], coordinates[second], coordinates[fourth]);
            double c = orientation(coordinates[third], coordinates[fourth], coordinates[first]);
            double d = orientation(coordinates[third], coordinates[fourth], coordinates[second]);
            if (a * b < -1e-12 && c * d < -1e-12) count++;
          }
    return count;
  }

  private static double orientation(MolGraph.Point2D first, MolGraph.Point2D second, MolGraph.Point2D third) {
    return (second.x - first.x) * (third.y - first.y) - (second.y - first.y) * (third.x - first.x);
  }
}

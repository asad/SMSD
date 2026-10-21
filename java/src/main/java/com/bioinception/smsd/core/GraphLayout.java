/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.core.MolGraph.Point2D;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Two-dimensional coordinates for a graph: scaffold templates, removal of bond crossings by
 * flipping ring systems and rings, force-directed relaxation and stress minimisation over graph
 * distances. Each method copies the coordinates it is given, completing missing ones on a line,
 * and writes its result back into them.
 */
final class GraphLayout {
  private GraphLayout() {}

  private static Point2D[] normaliseCoords(MolGraph g, Point2D[] coords, double spacing) {
    Point2D[] out = new Point2D[g.n];
    int copy = Math.min(coords == null ? 0 : coords.length, g.n);
    for (int i = 0; i < copy; i++) {
      Point2D p = coords[i];
      out[i] = p != null ? new Point2D(p.x, p.y) : new Point2D(i * spacing, 0.0);
    }
    for (int i = copy; i < g.n; i++) out[i] = new Point2D(i * spacing, 0.0);
    return out;
  }

  private static void copyCoords(Point2D[] src, Point2D[] dst) {
    int copy = Math.min(src.length, dst.length);
    for (int i = 0; i < copy; i++) {
      if (dst[i] == null) dst[i] = new Point2D(src[i].x, src[i].y);
      else {
        dst[i].x = src[i].x;
        dst[i].y = src[i].y;
      }
    }
  }

  private static boolean isDegenerateLayout(Point2D[] coords) {
    if (coords.length < 3) return true;
    double minX = coords[0].x, maxX = coords[0].x;
    double minY = coords[0].y, maxY = coords[0].y;
    for (Point2D p : coords) {
      minX = Math.min(minX, p.x);
      maxX = Math.max(maxX, p.x);
      minY = Math.min(minY, p.y);
      maxY = Math.max(maxY, p.y);
    }
    return (maxX - minX < 1e-6) || (maxY - minY < 1e-6);
  }

  /**
   * Reduces bond crossings by annealing: first whole ring systems, then single rings of fused
   * systems about their fusion atoms, are mirrored. Returns the crossings that remain.
   */
  static int reduceCrossings(MolGraph g, Point2D[] coords, int maxIter) {
    Point2D[] work = normaliseCoords(g, coords, 1.5);
    if (g.n < 4) return 0;

    int[][] bonds = buildBondList(g);
    List<Set<Integer>> ringSystems = RingPerception.fusedRingAtomSets(g);
    if (ringSystems.isEmpty()) {
      int result = countCrossingsImpl(bonds, work);
      copyCoords(work, coords);
      return result;
    }

    int nSys = ringSystems.size();
    int current = countCrossingsImpl(bonds, work);
    if (current == 0) {
      copyCoords(work, coords);
      return 0;
    }

    double T0 = 2.0, Tmin = 0.01, alpha = 0.995;
    Random rng = new Random(42);
    double T = T0;

    // ---- Phase 1: system-level flipping ----
    int phase1Iters = maxIter / 2;
    for (int iter = 0; iter < phase1Iters && current > 0; iter++) {
      int si = rng.nextInt(nSys);
      Set<Integer> sys = ringSystems.get(si);
      if (sys.size() < 3) continue;

      double cx = 0, cy = 0;
      for (int a : sys) { cx += work[a].x; cy += work[a].y; }
      cx /= sys.size(); cy /= sys.size();

      Map<Integer, double[]> saved = new HashMap<>();
      for (int a : sys) saved.put(a, new double[]{work[a].x, work[a].y});

      int axis = rng.nextInt(2);
      for (int a : sys) {
        if (axis == 0) work[a].x = 2.0 * cx - work[a].x;
        else           work[a].y = 2.0 * cy - work[a].y;
      }

      int newCr = countCrossingsImpl(bonds, work);
      int delta = newCr - current;
      boolean accept = (delta <= 0);
      if (!accept && T > Tmin) {
        accept = rng.nextDouble() < Math.exp(-(double) delta / T);
      }
      if (accept) { current = newCr; }
      else { for (Map.Entry<Integer, double[]> e : saved.entrySet()) {
        work[e.getKey()].x = e.getValue()[0];
        work[e.getKey()].y = e.getValue()[1];
      }}
      T *= alpha;
      if (T < Tmin) T = Tmin;
    }

    if (current == 0) {
      copyCoords(work, coords);
      return 0;
    }

    // ---- Phase 2: individual ring flipping within multi-ring systems ----
    int[][] sssr = g.computeRings();
    if (sssr.length < 2) {
      copyCoords(work, coords);
      return current;
    }

    // Build flippable rings: for each ring in a multi-ring system,
    // exclusive = atoms only in this ring, shared = atoms in ≥2 rings.
    // Shared atoms act as fixed pivots during the flip.
    List<int[]> exclusiveList = new ArrayList<>();
    List<int[]> sharedList = new ArrayList<>();

    // Group rings by system membership (which system does each SSSR ring belong to?)
    // and count per-system ring counts
    Map<Integer, List<Integer>> sysToRings = new HashMap<>();
    for (int ri = 0; ri < sssr.length; ri++) {
      for (int si = 0; si < nSys; si++) {
        Set<Integer> sys = ringSystems.get(si);
        boolean belongs = true;
        for (int a : sssr[ri]) {
          if (!sys.contains(a)) { belongs = false; break; }
        }
        if (belongs) {
          sysToRings.computeIfAbsent(si, k -> new ArrayList<>()).add(ri);
          break;
        }
      }
    }

    for (Map.Entry<Integer, List<Integer>> entry : sysToRings.entrySet()) {
      List<Integer> ringIndices = entry.getValue();
      if (ringIndices.size() < 2) continue; // single ring — handled in phase 1

      // Count how many rings each atom belongs to within this system
      Map<Integer, Integer> atomRingCount = new HashMap<>();
      for (int ri : ringIndices) {
        for (int a : sssr[ri]) atomRingCount.merge(a, 1, Integer::sum);
      }

      for (int ri : ringIndices) {
        List<Integer> excl = new ArrayList<>(), shar = new ArrayList<>();
        for (int a : sssr[ri]) {
          if (atomRingCount.getOrDefault(a, 0) >= 2) shar.add(a);
          else excl.add(a);
        }
        if (!excl.isEmpty() && !shar.isEmpty()) {
          exclusiveList.add(excl.stream().mapToInt(Integer::intValue).toArray());
          sharedList.add(shar.stream().mapToInt(Integer::intValue).toArray());
        }
      }
    }

    if (!exclusiveList.isEmpty()) {
      int nFlip = exclusiveList.size();
      T = T0; // reset temperature for phase 2
      int phase2Iters = maxIter - phase1Iters;

      for (int iter = 0; iter < phase2Iters && current > 0; iter++) {
        int fi = rng.nextInt(nFlip);
        int[] excl = exclusiveList.get(fi);
        int[] shar = sharedList.get(fi);

        // Pivot = centroid of shared (fusion) atoms
        double cx = 0, cy = 0;
        for (int a : shar) { cx += work[a].x; cy += work[a].y; }
        cx /= shar.length; cy /= shar.length;

        // Save exclusive atom coordinates
        double[][] saved = new double[excl.length][2];
        for (int i = 0; i < excl.length; i++) {
          saved[i][0] = work[excl[i]].x;
          saved[i][1] = work[excl[i]].y;
        }

        // Mirror exclusive atoms about the shared-atom centroid
        int axis = rng.nextInt(2);
        for (int a : excl) {
          if (axis == 0) work[a].x = 2.0 * cx - work[a].x;
          else           work[a].y = 2.0 * cy - work[a].y;
        }

        int newCr = countCrossingsImpl(bonds, work);
        int delta = newCr - current;
        boolean accept = (delta <= 0);
        if (!accept && T > Tmin) {
          accept = rng.nextDouble() < Math.exp(-(double) delta / T);
        }
        if (accept) { current = newCr; }
        else {
          for (int i = 0; i < excl.length; i++) {
            work[excl[i]].x = saved[i][0];
            work[excl[i]].y = saved[i][1];
          }
        }
        T *= alpha;
        if (T < Tmin) T = Tmin;
      }
    }

    copyCoords(work, coords);
    return current;
  }

  /**
   * Relaxes bond lengths towards {@code targetBondLength}, pushing unbonded atoms and crossing
   * bonds apart; returns the remaining bond-length stress.
   */
  static double forceDirectedLayout(MolGraph g, Point2D[] coords, int maxIter, double targetBondLength) {
    int n = g.n;
    if (n < 2) return 0.0;
    Point2D[] work = normaliseCoords(g, coords, targetBondLength);
    if (isDegenerateLayout(work)) {
      Point2D[] templ = matchTemplate(g, targetBondLength);
      if (templ != null && templ.length == n) work = templ;
    }

    final double kAttract  = 1.0;
    final double kRepel    = 0.5;
    final double kCrossing = 2.0;
    double step = 0.1;
    final double cooling = 0.98;
    final double convergence = 0.01;

    // Adjacency for O(1) bonded check
    boolean[][] bonded = new boolean[n][n];
    for (int i = 0; i < n; i++)
      for (int j : g.neighbors[i])
        bonded[i][j] = true;

    // Bond list
    int[][] bonds = buildBondList(g);
    int nb = bonds.length;

    for (int iter = 0; iter < maxIter; iter++) {
      double[] fx = new double[n], fy = new double[n];

      // 1. Attractive forces (bonded pairs)
      for (int[] bond : bonds) {
        int a = bond[0], b = bond[1];
        double dx = work[b].x - work[a].x;
        double dy = work[b].y - work[a].y;
        double dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 1e-10) dist = 1e-10;
        double force = kAttract * (dist - targetBondLength);
        double ux = dx / dist, uy = dy / dist;
        fx[a] += force * ux; fy[a] += force * uy;
        fx[b] -= force * ux; fy[b] -= force * uy;
      }

      // 2. Repulsive forces (non-bonded pairs)
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          if (bonded[i][j]) continue;
          double dx = work[j].x - work[i].x;
          double dy = work[j].y - work[i].y;
          double distSq = dx * dx + dy * dy;
          if (distSq < 1e-6) distSq = 1e-6;
          double force = -kRepel / distSq;
          double dist = Math.sqrt(distSq);
          double ux = dx / dist, uy = dy / dist;
          fx[i] += force * ux; fy[i] += force * uy;
          fx[j] -= force * ux; fy[j] -= force * uy;
        }
      }

      // 3. Crossing penalty forces
      for (int bi = 0; bi < nb; bi++) {
        for (int bj = bi + 1; bj < nb; bj++) {
          int a1 = bonds[bi][0], a2 = bonds[bi][1];
          int b1 = bonds[bj][0], b2 = bonds[bj][1];
          if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
          if (!segmentsCross(work[a1], work[a2], work[b1], work[b2])) continue;
          double m1x = (work[a1].x + work[a2].x) * 0.5;
          double m1y = (work[a1].y + work[a2].y) * 0.5;
          double m2x = (work[b1].x + work[b2].x) * 0.5;
          double m2y = (work[b1].y + work[b2].y) * 0.5;
          double dx = m2x - m1x, dy = m2y - m1y;
          double dist = Math.sqrt(dx * dx + dy * dy);
          if (dist < 1e-8) { dx = 0.1; dy = 0.1; dist = Math.sqrt(0.02); }
          double force = -kCrossing / (dist + 0.1);
          double ux = dx / dist, uy = dy / dist;
          fx[a1] += 0.5 * force * ux; fy[a1] += 0.5 * force * uy;
          fx[a2] += 0.5 * force * ux; fy[a2] += 0.5 * force * uy;
          fx[b1] -= 0.5 * force * ux; fy[b1] -= 0.5 * force * uy;
          fx[b2] -= 0.5 * force * ux; fy[b2] -= 0.5 * force * uy;
        }
      }

      // Update positions
      double maxDisp = 0.0;
      for (int i = 0; i < n; i++) {
        double dx = step * fx[i], dy = step * fy[i];
        double disp = Math.sqrt(dx * dx + dy * dy);
        if (disp > maxDisp) maxDisp = disp;
        work[i].x += dx;
        work[i].y += dy;
      }
      step *= cooling;
      if (maxDisp < convergence) break;
    }

    // Compute final stress
    double stress = 0.0;
    for (int[] bond : bonds) {
      double dx = work[bond[1]].x - work[bond[0]].x;
      double dy = work[bond[1]].y - work[bond[0]].y;
      double dist = Math.sqrt(dx * dx + dy * dy);
      double delta = dist - targetBondLength;
      stress += delta * delta;
    }
    copyCoords(work, coords);
    return stress;
  }

  /**
   * Places atoms so that their distances approach their graph distances times
   * {@code targetBondLength}, from the given and two perturbed starts; returns the lowest stress.
   */
  static double stressMajorisation(MolGraph g, Point2D[] coords, int maxIter, double targetBondLength) {
    int n = g.n;
    if (n < 2) return 0.0;
    Point2D[] seedCoords = normaliseCoords(g, coords, targetBondLength);
    if (isDegenerateLayout(seedCoords)) {
      Point2D[] templ = matchTemplate(g, targetBondLength);
      if (templ != null && templ.length == n) seedCoords = templ;
    }

    // BFS all-pairs shortest path
    int[][] graphDist = new int[n][n];
    for (int src = 0; src < n; src++) {
      Arrays.fill(graphDist[src], n + 1);
      graphDist[src][src] = 0;
      Deque<Integer> q = new ArrayDeque<>();
      q.add(src);
      while (!q.isEmpty()) {
        int u = q.poll();
        for (int v : g.neighbors[u]) {
          if (graphDist[src][v] > graphDist[src][u] + 1) {
            graphDist[src][v] = graphDist[src][u] + 1;
            q.add(v);
          }
        }
      }
    }

    // Target distances and weights
    double[][] D = new double[n][n];
    double[][] W = new double[n][n];
    double[] Wsum = new double[n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        D[i][j] = graphDist[i][j] * targetBondLength;
        if (i != j && D[i][j] > 1e-10) {
          W[i][j] = 1.0 / (D[i][j] * D[i][j]);
        }
        Wsum[i] += W[i][j];
      }
    }

    final double eps = 1e-6;
    Random rng = new Random(42);
    double bestStress = Double.POSITIVE_INFINITY;
    Point2D[] bestCoords = seedCoords;
    int nInit = 3;
    for (int init = 0; init < Math.max(1, nInit); init++) {
      Point2D[] runCoords = normaliseCoords(g, seedCoords, targetBondLength);
      if (init > 0) {
        for (int i = 0; i < n; i++) {
          runCoords[i].x += (rng.nextDouble() - 0.5) * targetBondLength;
          runCoords[i].y += (rng.nextDouble() - 0.5) * targetBondLength;
        }
      }

      double prevStress = 0.0;
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          double dx = runCoords[i].x - runCoords[j].x;
          double dy = runCoords[i].y - runCoords[j].y;
          double dij = Math.sqrt(dx * dx + dy * dy);
          double delta = dij - D[i][j];
          prevStress += W[i][j] * delta * delta;
        }
      }

      for (int iter = 0; iter < maxIter; iter++) {
        Point2D[] nc = new Point2D[n];
        for (int i = 0; i < n; i++) {
          double sx = 0, sy = 0;
          for (int j = 0; j < n; j++) {
            if (i == j || W[i][j] < 1e-15) continue;
            double dx = runCoords[i].x - runCoords[j].x;
            double dy = runCoords[i].y - runCoords[j].y;
            double dij = Math.sqrt(dx * dx + dy * dy);
            if (dij < 1e-10) dij = 1e-10;
            double ratio = D[i][j] / dij;
            sx += W[i][j] * (runCoords[j].x + ratio * (runCoords[i].x - runCoords[j].x));
            sy += W[i][j] * (runCoords[j].y + ratio * (runCoords[i].y - runCoords[j].y));
          }
          nc[i] = Wsum[i] > 1e-15
              ? new Point2D(sx / Wsum[i], sy / Wsum[i])
              : new Point2D(runCoords[i].x, runCoords[i].y);
        }
        runCoords = nc;

        double stress = 0.0;
        for (int i = 0; i < n; i++) {
          for (int j = i + 1; j < n; j++) {
            double dx = runCoords[i].x - runCoords[j].x;
            double dy = runCoords[i].y - runCoords[j].y;
            double dij = Math.sqrt(dx * dx + dy * dy);
            double delta = dij - D[i][j];
            stress += W[i][j] * delta * delta;
          }
        }
        if (prevStress > 1e-10 && (prevStress - stress) / prevStress < eps) {
          prevStress = stress;
          break;
        }
        prevStress = stress;
      }

      if (prevStress < bestStress) {
        bestStress = prevStress;
        bestCoords = runCoords;
      }
    }
    copyCoords(bestCoords, coords);
    return bestStress;
  }

  /**
   * Coordinates from a scaffold template with the atom count and ring sizes of the graph whose
   * bonds it matches, scaled to {@code targetBondLength}, or null.
   */
  static Point2D[] matchTemplate(MolGraph g, double targetBondLength) {
    int[][] rings = g.computeRings();
    if (rings.length == 0 && g.n != 6) return null;

    int[] sizes = new int[rings.length];
    for (int i = 0; i < rings.length; i++) sizes[i] = rings[i].length;
    Arrays.sort(sizes);

    // Try each template
    for (TemplateEntry t : SCAFFOLD_TEMPLATES) {
      if (t.atomCount != g.n) continue;
      int[] tSizes = Arrays.copyOf(t.ringSizes, t.ringSizes.length);
      Arrays.sort(tSizes);
      if (!Arrays.equals(sizes, tSizes)) continue;
      int[] mapping = templateMapping(g, t.coords);
      if (mapping == null) continue;
      Point2D[] result = new Point2D[t.atomCount];
      for (int i = 0; i < t.atomCount; i++)
        result[mapping[i]] = new Point2D(t.coords[i][0] * targetBondLength,
                                t.coords[i][1] * targetBondLength);
      return result;
    }
    return null;
  }

  /** Map the template's unit-length edges to the molecule before assigning coordinates. */
  private static int[] templateMapping(MolGraph graph, double[][] coordinates) {
    int atoms = coordinates.length;
    boolean[][] edges = new boolean[atoms][atoms];
    int[] degree = new int[atoms];
    for (int first = 0; first < atoms; first++)
      for (int second = first + 1; second < atoms; second++) {
        double distance = Math.hypot(coordinates[first][0] - coordinates[second][0],
            coordinates[first][1] - coordinates[second][1]);
        if (Math.abs(distance - 1.0) < 1e-8) {
          edges[first][second] = edges[second][first] = true;
          degree[first]++;
          degree[second]++;
        }
      }
    Integer[] order = new Integer[atoms];
    for (int atom = 0; atom < atoms; atom++) order[atom] = atom;
    Arrays.sort(order, (first, second) -> Integer.compare(degree[second], degree[first]));
    int[] mapping = new int[atoms];
    Arrays.fill(mapping, -1);
    return matchTemplateAtoms(graph, edges, degree, order, mapping, new boolean[atoms], 0) ? mapping : null;
  }

  private static boolean matchTemplateAtoms(MolGraph graph, boolean[][] edges, int[] degree,
                                            Integer[] order, int[] mapping, boolean[] used, int depth) {
    if (depth == order.length) return true;
    int templateAtom = order[depth];
    for (int graphAtom = 0; graphAtom < graph.n; graphAtom++) {
      if (used[graphAtom] || degree[templateAtom] != graph.neighbors[graphAtom].length) continue;
      boolean compatible = true;
      for (int previous = 0; previous < depth; previous++) {
        int other = order[previous];
        if (edges[templateAtom][other] != graph.hasBond(graphAtom, mapping[other])) {
          compatible = false;
          break;
        }
      }
      if (!compatible) continue;
      mapping[templateAtom] = graphAtom;
      used[graphAtom] = true;
      if (matchTemplateAtoms(graph, edges, degree, order, mapping, used, depth + 1)) return true;
      used[graphAtom] = false;
    }
    return false;
  }

  /** A scaffold drawn with unit bonds: its atom count, ring sizes and coordinates. */
  private static final class TemplateEntry {
    final int atomCount;
    final int[] ringSizes;
    final double[][] coords;

    TemplateEntry(int atomCount, int[] ringSizes, double[][] coords) {
      this.atomCount = atomCount;
      this.ringSizes = ringSizes;
      this.coords = coords;
    }
  }

  private static double[][] regularPolygonCoords(int n) {
    double r = 0.5 / Math.sin(Math.PI / n);
    double[][] coords = new double[n][2];
    for (int i = 0; i < n; i++) {
      double angle = 2.0 * Math.PI * i / n - Math.PI / 2.0;
      coords[i][0] = r * Math.cos(angle);
      coords[i][1] = r * Math.sin(angle);
    }
    return coords;
  }

  /** One drawing for each ring pattern; the first template that matches is used. */
  private static final TemplateEntry[] SCAFFOLD_TEMPLATES = {
    // Benzene, cyclohexane and piperidine.
    new TemplateEntry(6, new int[]{6}, regularPolygonCoords(6)),
    // Naphthalene and quinoline: two fused six-membered rings.
    new TemplateEntry(10, new int[]{6, 6}, buildFused66Coords()),
    // Biphenyl: two six-membered rings joined by a bond.
    new TemplateEntry(12, new int[]{6, 6}, buildBiphenylCoords()),
    // Indole and purine: a five-membered ring fused to a six-membered ring.
    new TemplateEntry(9, new int[]{5, 6}, buildFused56Coords()),
  };

  private static double[][] buildFused66Coords() {
    return buildFusedCoords(6, 1, 2);
  }

  private static double[][] buildFused56Coords() {
    return buildFusedCoords(5, 2, 3);
  }

  /** Place the new ring on the opposite side of its shared edge. */
  private static double[][] buildFusedCoords(int ringSize, int firstFusion, int secondFusion) {
    double[][] hexagon = regularPolygonCoords(6);
    double[][] coordinates = new double[6 + ringSize - 2][2];
    System.arraycopy(hexagon, 0, coordinates, 0, 6);
    double midpointX = (hexagon[firstFusion][0] + hexagon[secondFusion][0]) * 0.5;
    double midpointY = (hexagon[firstFusion][1] + hexagon[secondFusion][1]) * 0.5;
    double edgeX = hexagon[secondFusion][0] - hexagon[firstFusion][0];
    double edgeY = hexagon[secondFusion][1] - hexagon[firstFusion][1];
    double length = Math.hypot(edgeX, edgeY);
    double normalX = -edgeY / length, normalY = edgeX / length;
    if (normalX * -midpointX + normalY * -midpointY > 0) {
      normalX = -normalX;
      normalY = -normalY;
    }
    double radius = 0.5 / Math.sin(Math.PI / ringSize);
    double apothem = radius * Math.cos(Math.PI / ringSize);
    double centreX = midpointX + apothem * normalX;
    double centreY = midpointY + apothem * normalY;
    double firstX = hexagon[firstFusion][0] - centreX;
    double firstY = hexagon[firstFusion][1] - centreY;
    double secondX = hexagon[secondFusion][0] - centreX;
    double secondY = hexagon[secondFusion][1] - centreY;
    double step = Math.copySign(2.0 * Math.PI / ringSize, firstX * secondY - firstY * secondX);
    double start = Math.atan2(secondY, secondX);
    for (int vertex = 0; vertex < ringSize - 2; vertex++) {
      double angle = start + step * (vertex + 1);
      coordinates[6 + vertex][0] = centreX + radius * Math.cos(angle);
      coordinates[6 + vertex][1] = centreY + radius * Math.sin(angle);
    }
    return coordinates;
  }

  private static double[][] buildBiphenylCoords() {
    double[][] hexagon = regularPolygonCoords(6);
    double[][] coordinates = new double[12][2];
    double cosine = Math.cos(Math.PI / 6), sine = Math.sin(Math.PI / 6);
    for (int vertex = 0; vertex < 6; vertex++) {
      double x = hexagon[vertex][0] * cosine - hexagon[vertex][1] * sine;
      double y = hexagon[vertex][0] * sine + hexagon[vertex][1] * cosine;
      coordinates[vertex][0] = x;
      coordinates[vertex][1] = y;
      coordinates[6 + vertex][0] = x + 3.0;
      coordinates[6 + vertex][1] = y;
    }
    return coordinates;
  }

  private static int[][] buildBondList(MolGraph g) {
    List<int[]> bonds = new ArrayList<>();
    for (int i = 0; i < g.n; i++)
      for (int j : g.neighbors[i])
        if (j > i) bonds.add(new int[]{i, j});
    return bonds.toArray(new int[0][]);
  }

  private static boolean segmentsCross(Point2D p1, Point2D p2, Point2D p3, Point2D p4) {
    double d1x = p2.x - p1.x, d1y = p2.y - p1.y;
    double d2x = p4.x - p3.x, d2y = p4.y - p3.y;
    double denom = d1x * d2y - d1y * d2x;
    if (Math.abs(denom) < 1e-12) return false;
    double dx = p3.x - p1.x, dy = p3.y - p1.y;
    double t = (dx * d2y - dy * d2x) / denom;
    double u = (dx * d1y - dy * d1x) / denom;
    double eps = 1e-9;
    return (t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps);
  }

  private static int countCrossingsImpl(int[][] bonds, Point2D[] coords) {
    int crossings = 0, nb = bonds.length;
    for (int i = 0; i < nb; i++) {
      for (int j = i + 1; j < nb; j++) {
        int a1 = bonds[i][0], a2 = bonds[i][1];
        int b1 = bonds[j][0], b2 = bonds[j][1];
        if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
        if (segmentsCross(coords[a1], coords[a2], coords[b1], coords[b2]))
          crossings++;
      }
    }
    return crossings;
  }
}

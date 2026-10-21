/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

/**
 * Canonical labelling of a graph by partition refinement and a bounded search of its
 * individualisations, with the automorphisms met on the way, the atom orbits they prove and a
 * canonical hash. Atom, bond and stereo properties take part, so equal labellings mean equal
 * chemical graphs.
 */
final class CanonicalLabels {
  private static final int HASH_PRIME = 1000003;
  private static final int MAX_SEARCH_NODES = 50000;
  private static final int CANON_SEARCH_LIMIT = 200;
  private static final int MAX_GENERATORS = 500;

  /** The canonical labelling of a graph and the chemical symmetry found while computing it. */
  static final class Labelling {
    /** The canonical position of each atom. */
    final int[] labels;
    /** The lowest atom of the proven orbit of each atom. */
    final int[] orbits;
    /** Automorphisms that keep every atom, bond and stereo property. */
    final int[][] generators;
    /** Whether the generators may describe only part of the automorphism group. */
    final boolean generatorsTruncated;
    final long hash;

    private Labelling(int[] labels, int[] orbits, int[][] generators, boolean generatorsTruncated, long hash) {
      this.labels = labels;
      this.orbits = orbits;
      this.generators = generators;
      this.generatorsTruncated = generatorsTruncated;
      this.hash = hash;
    }
  }

  /** The partition search result, before the chemical filter of its generators. */
  private static final class CanonResult {
    final int[] canonLabel;
    final int[] orbit;
    final int[][] autGenerators;
    final boolean generatorsTruncated;

    CanonResult(int[] canonLabel, int[] orbit, int[][] autGenerators, boolean generatorsTruncated) {
      this.canonLabel = canonLabel;
      this.orbit = orbit;
      this.autGenerators = autGenerators;
      this.generatorsTruncated = generatorsTruncated;
    }
  }

  private CanonicalLabels() {}

  /** The canonical labelling of {@code graph}. */
  static Labelling of(MolGraph graph) {
    int n = graph.n;
    int[][] neighbors = graph.neighbors;
    int[] hydrogens = new int[n];
    for (int atom = 0; atom < n; atom++) hydrogens[atom] = graph.hydrogenCount(atom);
    CanonResult clResult = computeCanonicalLabelling(n, atomLabels(graph, hydrogens), graph.degree, neighbors, graph);
    int[] canonicalLabel = clResult.canonLabel;
    boolean generatorsTruncated = clResult.generatorsTruncated;
    boolean[] stereoFixed = new boolean[n];
    for (int atom = 0; atom < n; atom++) {
      if (graph.tetraChirality != null && graph.tetraChirality[atom] != 0) {
        stereoFixed[atom] = true;
        for (int ligand : neighbors[atom]) stereoFixed[ligand] = true;
      }
      for (int neighbor : neighbors[atom]) {
        if (neighbor <= atom || graph.dbStereo(atom, neighbor) == 0) continue;
        stereoFixed[atom] = stereoFixed[neighbor] = true;
        for (int ligand : neighbors[atom]) stereoFixed[ligand] = true;
        for (int ligand : neighbors[neighbor]) stereoFixed[ligand] = true;
      }
    }
    List<int[]> chemicalGenerators = new ArrayList<>();
    int[] provenOrbits = UnionFind.create(n);
    for (int[] generator : clResult.autGenerators) {
      if (!preservesChemicalProperties(graph, generator, hydrogens, stereoFixed)) {
        // Valid subgroup generators can require products of rejected topology
        // generators, so filtering does not prove completeness of the subgroup.
        generatorsTruncated = true;
        continue;
      }
      chemicalGenerators.add(generator);
      for (int atom = 0; atom < n; atom++) UnionFind.join(provenOrbits, atom, generator[atom]);
    }
    // Unproved refinement-cell equivalence cannot justify chemical orbit pruning.
    return new Labelling(canonicalLabel, buildOrbits(provenOrbits, n), chemicalGenerators.toArray(new int[0][]),
        generatorsTruncated, computeCanonicalHash(graph, canonicalLabel));
  }

  private static boolean preservesChemicalProperties(MolGraph graph, int[] generator, int[] hydrogens,
                                                     boolean[] stereoFixed) {
    for (int atom = 0; atom < graph.n; atom++) {
      int image = generator[atom];
      if (graph.label[atom] != graph.label[image] || graph.formalCharge[atom] != graph.formalCharge[image]
          || graph.massNumber[atom] != graph.massNumber[image] || hydrogens[atom] != hydrogens[image]
          || graph.ring[atom] != graph.ring[image] || graph.aromatic[atom] != graph.aromatic[image]
          || (stereoFixed[atom] && image != atom)) return false;
      for (int neighbor : graph.neighbors[atom]) {
        int mappedNeighbor = generator[neighbor];
        if (graph.bondOrder(atom, neighbor) != graph.bondOrder(image, mappedNeighbor)
            || graph.bondInRing(atom, neighbor) != graph.bondInRing(image, mappedNeighbor)
            || graph.bondAromatic(atom, neighbor) != graph.bondAromatic(image, mappedNeighbor)
            || graph.dbStereo(atom, neighbor) != graph.dbStereo(image, mappedNeighbor)) return false;
      }
    }
    return true;
  }

  /** Ranks of the atoms by element, aromaticity, ring flag, charge, mass, hydrogens and chirality. */
  private static int[] atomLabels(MolGraph graph, int[] hydrogens) {
    int n = graph.n;
    int[][] properties = new int[n][];
    Integer[] atoms = new Integer[n];
    for (int atom = 0; atom < n; atom++) {
      properties[atom] = new int[] {graph.label[atom], graph.formalCharge[atom], graph.massNumber[atom],
          hydrogens[atom], graph.tetraChirality != null && graph.tetraChirality[atom] != 0 ? 1 : 0};
      atoms[atom] = atom;
    }
    Arrays.sort(atoms, (first, second) -> compareIntArrays(properties[first], properties[second]));
    int[] labels = new int[n];
    int rank = 0;
    for (int position = 0; position < n; position++) {
      if (position > 0 && !Arrays.equals(properties[atoms[position]], properties[atoms[position - 1]])) rank++;
      labels[atoms[position]] = rank;
    }
    return labels;
  }

  /** Compares two arrays element by element, and then by length. */
  static int compareIntArrays(int[] first, int[] second) {
    for (int index = 0; index < Math.min(first.length, second.length); index++) {
      int comparison = Integer.compare(first[index], second[index]);
      if (comparison != 0) return comparison;
    }
    return Integer.compare(first.length, second.length);
  }

  private static int canonicalBondLabel(MolGraph graph, int first, int second) {
    return (graph.bondOrder(first, second) << 3) | (graph.bondInRing(first, second) ? 4 : 0)
        | (graph.bondAromatic(first, second) ? 2 : 0) | (graph.dbStereo(first, second) != 0 ? 1 : 0);
  }

  private static int canonicalTetrahedralTag(MolGraph graph, int atom, int[] ranking) {
    int tag = graph.tetraChirality != null ? graph.tetraChirality[atom] : 0;
    if (tag == 0) return 0;
    // The canonical tag is 1 for clockwise.
    tag = 3 - tag;
    int[] ligands = graph.tetrahedralLigands(atom);
    if (ligands == null) return tag;
    boolean odd = false;
    for (int first = 0; first < ligands.length; first++) {
      int firstRank = ligands[first] < 0 ? -1 : ranking[ligands[first]];
      for (int second = first + 1; second < ligands.length; second++) {
        int secondRank = ligands[second] < 0 ? -1 : ranking[ligands[second]];
        if (firstRank > secondRank) odd = !odd;
      }
    }
    return odd ? 3 - tag : tag;
  }

  private static int canonicalDoubleBondTag(MolGraph graph, int first, int second, int[] ranking) {
    int configuration = graph.dbStereo(first, second);
    if (configuration == 0) return 0;
    for (int side = 0; side < 2; side++) {
      int atom = side == 0 ? first : second, other = side == 0 ? second : first;
      int reference = -1, canonicalReference = -1;
      for (int neighbor : graph.neighbors[atom]) {
        if (neighbor == other) continue;
        if (reference < 0) reference = neighbor;
        if (canonicalReference < 0 || ranking[neighbor] < ranking[canonicalReference]) canonicalReference = neighbor;
      }
      if (reference != canonicalReference) configuration = 3 - configuration;
    }
    return configuration;
  }

  private static CanonResult computeCanonicalLabelling(int n, int[] label, int[] degree,
                                                       int[][] neighbors, MolGraph graph) {
    if (n == 0) return new CanonResult(new int[0], new int[0], new int[0][], false);
    int[] mRank = new int[n];
    System.arraycopy(label, 0, mRank, 0, n);
    int[] mNew = new int[n];
    long[] scratch = new long[16];
    for (int iter = 0; iter < 8; iter++) {
      int prevDistinct;
      { HashSet<Integer> s = new HashSet<>(); for (int v = 0; v < n; v++) s.add(mRank[v]); prevDistinct = s.size(); }
      for (int v = 0; v < n; v++) {
        int[] nbs = neighbors[v];
        int d = nbs.length;
        if (d > scratch.length) scratch = new long[d];
        for (int k = 0; k < d; k++)
          scratch[k] = ((long) mRank[nbs[k]] << 32) | canonicalBondLabel(graph, v, nbs[k]);
        Arrays.sort(scratch, 0, d);
        int h = mRank[v] * HASH_PRIME;
        for (int k = 0; k < d; k++) {
          h = h * 31 + (int) (scratch[k] >> 32);
          h = h * 31 + (int) scratch[k];
        }
        mNew[v] = h;
      }
      System.arraycopy(mNew, 0, mRank, 0, n);
      int newDistinct;
      { HashSet<Integer> s = new HashSet<>(); for (int v = 0; v < n; v++) s.add(mRank[v]); newDistinct = s.size(); }
      if (newDistinct <= prevDistinct) break;
    }
    int[] uf = UnionFind.create(n);
    Integer[] sorted = new Integer[n];
    for (int i = 0; i < n; i++) sorted[i] = i;
    final int[] mr = mRank;
    Arrays.sort(sorted, (a, b) -> {
      int c = Integer.compare(label[a], label[b]);
      if (c != 0) return c;
      c = Integer.compare(degree[a], degree[b]);
      return c != 0 ? c : Integer.compare(mr[a], mr[b]);
    });
    int[] initPerm = new int[n];
    boolean[] initCellEnd = new boolean[n];
    for (int i = 0; i < n; i++) initPerm[i] = sorted[i];
    for (int i = 0; i < n - 1; i++) {
      int vi = sorted[i], vj = sorted[i + 1];
      initCellEnd[i] = (label[vi] != label[vj] || degree[vi] != degree[vj] || mr[vi] != mr[vj]);
    }
    initCellEnd[n - 1] = true;
    if (n > CANON_SEARCH_LIMIT) {
      int cs = 0;
      for (int i = 0; i < n; i++) {
        if (initCellEnd[i]) {
          for (int j = cs + 1; j <= i; j++) UnionFind.join(uf, initPerm[cs], initPerm[j]);
          cs = i + 1;
        }
      }
      int[] canonLabel = new int[n];
      for (int pos = 0; pos < n; pos++) canonLabel[initPerm[pos]] = pos;
      return new CanonResult(canonLabel, buildOrbits(uf, n), new int[0][], true);
    }
    refinePartition(n, initPerm, initCellEnd, neighbors, label, -1, graph);
    int[] bestPerm = null;
    List<int[]> generators = new ArrayList<>();
    boolean generatorsTruncated = false;
    int nodeCount = 0;
    boolean budgetExceeded = false;
    ArrayDeque<int[]> permStack = new ArrayDeque<>();
    ArrayDeque<boolean[]> cellEndStack = new ArrayDeque<>();
    int targetCell = findSmallestNonSingleton(n, initCellEnd);
    if (targetCell < 0) {
      bestPerm = initPerm.clone();
    } else {
      int cellStart = 0;
      for (int i = 0; i < n; i++) {
        if (i > 0 && initCellEnd[i - 1]) cellStart = i;
        if (i == targetCell) break;
      }
      for (int i = targetCell; i >= cellStart; i--) {
        int[] p = initPerm.clone();
        boolean[] ce = initCellEnd.clone();
        if (i != cellStart) {
          int tmp = p[i];
          System.arraycopy(p, cellStart, p, cellStart + 1, i - cellStart);
          p[cellStart] = tmp;
        }
        ce[cellStart] = true;
        refinePartition(n, p, ce, neighbors, label, cellStart, graph);
        permStack.push(p);
        cellEndStack.push(ce);
        nodeCount++;
      }
    }
    while (!permStack.isEmpty() && !budgetExceeded) {
      int[] perm = permStack.pop();
      boolean[] cellEnd = cellEndStack.pop();
      int nc = findSmallestNonSingleton(n, cellEnd);
      if (nc < 0) {
        if (bestPerm == null || comparePerm(perm, bestPerm, n, label, neighbors, graph) < 0) {
          bestPerm = perm.clone();
        } else if (comparePerm(perm, bestPerm, n, label, neighbors, graph) == 0) {
          for (int i = 0; i < n; i++) UnionFind.join(uf, perm[i], bestPerm[i]);
          if (generators.size() < MAX_GENERATORS) {
            int[] gen = new int[n];
            for (int i = 0; i < n; i++) gen[bestPerm[i]] = perm[i];
            generators.add(gen);
          } else {
            generatorsTruncated = true;
          }
        }
        continue;
      }
      if (bestPerm != null && partialCompare(perm, bestPerm, n, cellEnd, label, neighbors) > 0) continue;
      int cellStart = 0;
      for (int i = 0; i < n; i++) {
        if (i > 0 && cellEnd[i - 1]) cellStart = i;
        if (i == nc) break;
      }
      int cellSize = nc - cellStart + 1;
      if (nodeCount + cellSize > MAX_SEARCH_NODES) { budgetExceeded = true; break; }
      for (int i = nc; i >= cellStart; i--) {
        int[] p = perm.clone();
        boolean[] ce = cellEnd.clone();
        if (i != cellStart) {
          int tmp = p[i];
          System.arraycopy(p, cellStart, p, cellStart + 1, i - cellStart);
          p[cellStart] = tmp;
        }
        ce[cellStart] = true;
        refinePartition(n, p, ce, neighbors, label, cellStart, graph);
        permStack.push(p);
        cellEndStack.push(ce);
        nodeCount++;
      }
    }
    if (bestPerm == null) bestPerm = initPerm.clone();
    int[] canonLabel = new int[n];
    for (int pos = 0; pos < n; pos++) canonLabel[bestPerm[pos]] = pos;
    int cs = 0;
    for (int i = 0; i < n; i++) {
      if (initCellEnd[i]) {
        for (int j = cs + 1; j <= i; j++) UnionFind.join(uf, initPerm[cs], initPerm[j]);
        cs = i + 1;
      }
    }
    return new CanonResult(canonLabel, buildOrbits(uf, n),
        generators.toArray(new int[0][]), generatorsTruncated || budgetExceeded);
  }

  private static int[] buildOrbits(int[] uf, int n) {
    int[] orbit = new int[n];
    for (int i = 0; i < n; i++) orbit[i] = UnionFind.root(uf, i);
    int[] orbitMin = new int[n];
    Arrays.fill(orbitMin, n);
    for (int i = 0; i < n; i++) { int r = orbit[i]; if (i < orbitMin[r]) orbitMin[r] = i; }
    for (int i = 0; i < n; i++) orbit[i] = orbitMin[orbit[i]];
    return orbit;
  }

  private static int findSmallestNonSingleton(int n, boolean[] cellEnd) {
    int bestSize = n + 1, bestEnd = -1, cellStart = 0;
    for (int i = 0; i < n; i++) {
      if (cellEnd[i]) {
        int size = i - cellStart + 1;
        if (size > 1 && size < bestSize) { bestSize = size; bestEnd = i; }
        cellStart = i + 1;
      }
    }
    return bestEnd;
  }

  private static void refinePartition(
      int n, int[] perm, boolean[] cellEnd, int[][] neighbors, int[] label, int targetPos, MolGraph graph) {
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[perm[i]] = i;
    boolean changed = true;
    int maxIter = n * 2 + 10;
    while (changed && maxIter-- > 0) {
      changed = false;
      int sStart = 0;
      for (int sEnd = 0; sEnd < n; sEnd++) {
        if (!cellEnd[sEnd]) continue;
        int cStart = 0;
        for (int cEnd = 0; cEnd < n; cEnd++) {
          if (!cellEnd[cEnd]) continue;
          int cSize = cEnd - cStart + 1;
          if (cSize <= 1) { cStart = cEnd + 1; continue; }
          int[][] count = new int[cSize][];
          for (int ci = 0; ci < cSize; ci++) {
            int atom = perm[cStart + ci], edgeCount = 0;
            int[] edges = new int[neighbors[atom].length];
            for (int neighbor : neighbors[atom]) {
              int position = inv[neighbor];
              if (position >= sStart && position <= sEnd)
                edges[edgeCount++] = canonicalBondLabel(graph, atom, neighbor);
            }
            Arrays.sort(edges, 0, edgeCount);
            count[ci] = Arrays.copyOf(edges, edgeCount);
          }
          boolean allSame = true;
          for (int ci = 1; ci < cSize; ci++) if (!Arrays.equals(count[ci], count[0])) { allSame = false; break; }
          if (allSame) { cStart = cEnd + 1; continue; }
          Integer[] idx = new Integer[cSize];
          for (int ci = 0; ci < cSize; ci++) idx[ci] = ci;
          final int fcs = cStart;
          final int[][] cnt = count;
          final int[] prm = perm;
          Arrays.sort(idx, (a, bx) -> {
            int c = compareIntArrays(cnt[a], cnt[bx]);
            if (c != 0) return c;
            c = Integer.compare(label[prm[fcs + a]], label[prm[fcs + bx]]);
            return c != 0 ? c : Integer.compare(prm[fcs + a], prm[fcs + bx]);
          });
          int[] tmp = new int[cSize];
          for (int ci = 0; ci < cSize; ci++) tmp[ci] = perm[cStart + idx[ci]];
          System.arraycopy(tmp, 0, perm, cStart, cSize);
          for (int ci = 0; ci < cSize; ci++) inv[perm[cStart + ci]] = cStart + ci;
          for (int ci = 0; ci < cSize - 1; ci++) {
            boolean wasBoundary = (cStart + ci == cEnd);
            boolean newBoundary = !Arrays.equals(cnt[idx[ci]], cnt[idx[ci + 1]]);
            cellEnd[cStart + ci] = newBoundary || wasBoundary;
            if (newBoundary && !wasBoundary) changed = true;
          }
          cellEnd[cEnd] = true;
          cStart = cEnd + 1;
        }
        sStart = sEnd + 1;
      }
    }
  }

  private static int comparePerm(int[] firstPermutation, int[] secondPermutation, int atomCount,
                                 int[] labels, int[][] neighbors, MolGraph graph) {
    int[] firstRanking = new int[atomCount], secondRanking = new int[atomCount];
    for (int position = 0; position < atomCount; position++) {
      firstRanking[firstPermutation[position]] = position;
      secondRanking[secondPermutation[position]] = position;
    }
    for (int position = 0; position < atomCount; position++) {
      int firstAtom = firstPermutation[position], secondAtom = secondPermutation[position];
      int comparison = Integer.compare(labels[firstAtom], labels[secondAtom]);
      if (comparison != 0) return comparison;
      comparison = Integer.compare(canonicalTetrahedralTag(graph, firstAtom, firstRanking),
          canonicalTetrahedralTag(graph, secondAtom, secondRanking));
      if (comparison != 0) return comparison;
      comparison = Integer.compare(neighbors[firstAtom].length, neighbors[secondAtom].length);
      if (comparison != 0) return comparison;
      long[][] signatures = new long[2][];
      for (int side = 0; side < 2; side++) {
        int atom = side == 0 ? firstAtom : secondAtom;
        int[] ranking = side == 0 ? firstRanking : secondRanking;
        signatures[side] = new long[neighbors[atom].length];
        for (int bondIndex = 0; bondIndex < neighbors[atom].length; bondIndex++) {
          int neighbor = neighbors[atom][bondIndex];
          signatures[side][bondIndex] = ((long) ranking[neighbor] << 16)
              | ((long) canonicalBondLabel(graph, atom, neighbor) << 2)
              | canonicalDoubleBondTag(graph, atom, neighbor, ranking);
        }
        Arrays.sort(signatures[side]);
      }
      for (int bondIndex = 0; bondIndex < signatures[0].length; bondIndex++) {
        comparison = Long.compare(signatures[0][bondIndex], signatures[1][bondIndex]);
        if (comparison != 0) return comparison;
      }
    }
    return 0;
  }

  private static int partialCompare(
      int[] perm, int[] bestPerm, int n, boolean[] cellEnd, int[] label, int[][] neighbors) {
    int cellStart = 0;
    for (int pos = 0; pos < n; pos++) {
      if (cellEnd[pos]) {
        if (pos == cellStart) {
          int c = Integer.compare(label[perm[pos]], label[bestPerm[pos]]);
          if (c != 0) return c;
        }
        cellStart = pos + 1;
      }
    }
    return 0;
  }

  private static long computeCanonicalHash(MolGraph g, int[] canonLabel) {
    int n = g.n;
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[canonLabel[i]] = i;
    long hash = 0;
    for (int cp = 0; cp < n; cp++) {
      int v = inv[cp];
      hash = hash * HASH_PRIME + g.label[v];
      hash = hash * HASH_PRIME + (long) g.formalCharge[v] + 512L;
      hash = hash * HASH_PRIME + g.massNumber[v];
      hash = hash * HASH_PRIME + g.hydrogenCount(v);
      hash = hash * HASH_PRIME + canonicalTetrahedralTag(g, v, canonLabel);
    }
    for (int cp = 0; cp < n; cp++) {
      int v = inv[cp];
      int[] nbs = g.neighbors[v];
      long[] edgeSigs = new long[nbs.length];
      for (int k = 0; k < nbs.length; k++) {
        int cn = canonLabel[nbs[k]];
        long edgeKey = (long) Math.min(cp, cn) * n + Math.max(cp, cn);
        long sig = edgeKey;
        sig = sig * HASH_PRIME + g.bondOrder(v, nbs[k]);
        sig = sig * 2L + (g.bondInRing(v, nbs[k]) ? 1L : 0L);
        sig = sig * 2L + (g.bondAromatic(v, nbs[k]) ? 1L : 0L);
        // include E/Z double-bond stereo in hash (symmetric: use min canonical position)
        if (cp < cn) {
          int dbConf = canonicalDoubleBondTag(g, v, nbs[k], canonLabel);
          sig = sig * HASH_PRIME + dbConf;
        }
        edgeSigs[k] = sig;
      }
      Arrays.sort(edgeSigs);
      for (long sig : edgeSigs) hash = hash * HASH_PRIME + sig;
    }
    return hash;
  }
}

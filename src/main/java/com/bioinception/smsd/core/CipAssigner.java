/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;

/**
 * CIP (Cahn-Ingold-Prelog) stereodescriptor assignment engine for MolGraph.
 *
 * Assigns R/S descriptors to tetrahedral stereocentres and E/Z descriptors
 * to double bonds using BFS digraph traversal with CIP sequence rules:
 *   Rule 1: higher atomic number gets higher priority
 *   Rule 2: higher mass number breaks ties (isotopes)
 *
 * Implements the digraph (extended connectivity tree) approach with phantom
 * (duplicate) nodes for ring closures and correct handling of multiple bonds.
 *
 * References:
 *   Prelog and Helmchen, Angew. Chem. Int. Ed. Engl. 21, 567 (1982)
 *   IUPAC Recommendations 2013 (P-92)
 *
 * @author Syed Asad Rahman
 */
public final class CipAssigner {

  /** Maximum BFS depth for the CIP digraph expansion.
   *  Must match C++ maxDepth (30) for long-chain resolution parity. */
  public static final int MAX_DEPTH = 30;

  /** Maximum number of nodes in the digraph to prevent runaway on large molecules. */
  private static final int MAX_NODES = 5000;

  private CipAssigner() {} // utility class

  // ==========================================================================
  // Default valence table (mirrors FingerprintEngine / smiles_parser.hpp)
  // ==========================================================================
  private static final int[][] DEFAULT_VALENCES = new int[54][];
  static {
    DEFAULT_VALENCES[1]  = new int[]{1};           // H
    DEFAULT_VALENCES[5]  = new int[]{3};           // B
    DEFAULT_VALENCES[6]  = new int[]{4};           // C
    DEFAULT_VALENCES[7]  = new int[]{3, 5};        // N
    DEFAULT_VALENCES[8]  = new int[]{2};           // O
    DEFAULT_VALENCES[9]  = new int[]{1};           // F
    DEFAULT_VALENCES[13] = new int[]{3};           // Al
    DEFAULT_VALENCES[14] = new int[]{4};           // Si
    DEFAULT_VALENCES[15] = new int[]{3, 5};        // P
    DEFAULT_VALENCES[16] = new int[]{2, 4, 6};     // S
    DEFAULT_VALENCES[17] = new int[]{1};           // Cl
    DEFAULT_VALENCES[33] = new int[]{3, 5};        // As
    DEFAULT_VALENCES[34] = new int[]{2, 4, 6};     // Se
    DEFAULT_VALENCES[35] = new int[]{1};           // Br
    DEFAULT_VALENCES[51] = new int[]{3, 5};        // Sb
    DEFAULT_VALENCES[52] = new int[]{2, 4, 6};     // Te
    DEFAULT_VALENCES[53] = new int[]{1, 3, 5, 7};  // I
  }

  /**
   * Compute implicit hydrogen count for an atom in the graph using the
   * OpenSMILES multi-valence model.
   */
  static int implicitH(MolGraph g, int idx) {
    int bondOrdSum = 0;
    for (int nb : g.neighbors[idx]) {
      int bo = g.bondOrder(idx, nb);
      if (bo == 4) bo = 1; // aromatic bond -> single for valence calc
      bondOrdSum += bo;
    }
    // Aromatic atoms contribute 1 pi electron to the bond order sum
    if (g.aromatic[idx]) bondOrdSum += 1;
    int z = g.atomicNum[idx];
    int[] valences = (z >= 0 && z < DEFAULT_VALENCES.length) ? DEFAULT_VALENCES[z] : null;
    if (valences == null) return 0;
    int charge = g.formalCharge[idx];
    // Charge-aware valence expansion per OpenSMILES:
    // Group 15/16 cations (N+, P+, O+, S+): valence increases by 1
    // Group 13 anions (B-): valence increases by 1
    boolean cationExpand = charge > 0 && (z == 7 || z == 15 || z == 8 || z == 16 || z == 34);
    boolean anionExpand = charge < 0 && (z == 5 || z == 13);
    for (int v : valences) {
      int target = v;
      if (cationExpand) target = v + charge;
      else if (anionExpand) target = v - charge; // charge is negative, so subtracting adds
      else target = v - Math.abs(charge);
      if (target >= bondOrdSum) return Math.max(0, target - bondOrdSum);
    }
    return 0;
  }

  // ==========================================================================
  // Digraph node for CIP priority determination
  // ==========================================================================
  private static final class DigraphNode {
    int atomIdx;       // index into MolGraph (-1 for implicit H)
    int atomicNum;
    int massNumber;
    boolean isDuplicate;
    int parent;        // index in node list (-1 for root)
    int depth;
    final List<Integer> children = new ArrayList<>(4);
  }

  /**
   * Build the CIP digraph (extended connectivity tree) rooted at the
   * specified atom. Ring closures produce duplicate (phantom) nodes.
   * Double and triple bonds produce additional phantom copies.
   */
  private static List<DigraphNode> buildDigraph(MolGraph g, int root) {
    List<DigraphNode> nodes = new ArrayList<>(128);

    DigraphNode rootNode = new DigraphNode();
    rootNode.atomIdx = root;
    rootNode.atomicNum = g.atomicNum[root];
    rootNode.massNumber = g.massNumber[root];
    rootNode.isDuplicate = false;
    rootNode.parent = -1;
    rootNode.depth = 0;
    nodes.add(rootNode);

    List<Integer> frontier = new ArrayList<>();
    frontier.add(0);

    for (int depth = 0; depth < MAX_DEPTH && !frontier.isEmpty(); ++depth) {
      List<Integer> nextFrontier = new ArrayList<>();
      for (int ni : frontier) {
        if (nodes.size() >= MAX_NODES) break;

        int atom = nodes.get(ni).atomIdx;
        int parentNodeIdx = nodes.get(ni).parent;
        int parentAtom = (parentNodeIdx >= 0) ? nodes.get(parentNodeIdx).atomIdx : -1;

        // Collect expansions: real children + phantoms for multiple bonds
        List<int[]> expansions = new ArrayList<>(); // {neighborAtom, isDuplicate}

        for (int nbr : g.neighbors[atom]) {
          if (nbr == parentAtom && !nodes.get(ni).isDuplicate) continue;

          int bo = g.bondOrder(atom, nbr);
          if (bo == 4) bo = 1; // aromatic -> single for expansion

          expansions.add(new int[]{nbr, 0}); // real child
          for (int extra = 1; extra < bo; ++extra) {
            expansions.add(new int[]{nbr, 1}); // phantom for bond order
          }
        }

        // Phantoms for bond back to parent (higher-order bonds)
        if (parentAtom >= 0 && !nodes.get(ni).isDuplicate) {
          int bo = g.bondOrder(atom, parentAtom);
          if (bo == 4) bo = 1;
          for (int extra = 1; extra < bo; ++extra) {
            expansions.add(new int[]{parentAtom, 1});
          }
        }

        for (int[] exp : expansions) {
          if (nodes.size() >= MAX_NODES) break;
          int nbrAtom = exp[0];
          boolean forceDuplicate = (exp[1] == 1);

          // Check if this atom is on the path from root to current node
          boolean onPath = false;
          if (!forceDuplicate) {
            int walk = ni;
            while (walk >= 0) {
              if (nodes.get(walk).atomIdx == nbrAtom) {
                onPath = true;
                break;
              }
              walk = nodes.get(walk).parent;
            }
          }

          boolean makeDuplicate = forceDuplicate || onPath;

          DigraphNode child = new DigraphNode();
          child.atomIdx = nbrAtom;
          child.atomicNum = g.atomicNum[nbrAtom];
          child.massNumber = makeDuplicate ? 0 : g.massNumber[nbrAtom];
          child.isDuplicate = makeDuplicate;
          child.parent = ni;
          child.depth = depth + 1;

          int childIdx = nodes.size();
          nodes.add(child);
          nodes.get(ni).children.add(childIdx);

          if (!makeDuplicate) {
            nextFrontier.add(childIdx);
          }
        }
      }
      frontier = nextFrontier;
    }

    return nodes;
  }

  /**
   * Compare two digraph sub-trees by CIP Rules 1 (atomic number) and
   * 2 (mass number) using BFS level-by-level comparison.
   *
   * @return negative if rootA has lower priority, positive if higher, 0 if equal
   */
  private static int compareBranches(List<DigraphNode> nodesA, int rootA,
                                     List<DigraphNode> nodesB, int rootB) {
    return compareBranches(nodesA, rootA, nodesB, rootB, null);
  }

  /**
   * Compare two digraph sub-trees by CIP Rules 1-2, with optional Rule 4b-c
   * (like/unlike descriptor pairing) and Rule 5 (R > S tiebreak).
   *
   * @param descriptorMap  atom index → descriptor (1=R, 2=S); null for Rules 1-2 only
   * @since 6.5.2
   */
  private static int compareBranches(List<DigraphNode> nodesA, int rootA,
                                     List<DigraphNode> nodesB, int rootB,
                                     Map<Integer, Integer> descriptorMap) {
    List<Integer> frontA = new ArrayList<>();
    frontA.add(rootA);
    List<Integer> frontB = new ArrayList<>();
    frontB.add(rootB);

    while (!frontA.isEmpty() || !frontB.isEmpty()) {
      // Collect atomic numbers at this level, sorted descending
      int[] zA = toSortedDescArray(frontA, nodesA, true);
      int[] zB = toSortedDescArray(frontB, nodesB, true);

      int len = Math.max(zA.length, zB.length);
      for (int k = 0; k < len; k++) {
        int a = k < zA.length ? zA[k] : 0;
        int b = k < zB.length ? zB[k] : 0;
        if (a != b) return a - b;
      }

      // Tie on atomic number: try mass number
      int[] mA = toSortedDescArray(frontA, nodesA, false);
      int[] mB = toSortedDescArray(frontB, nodesB, false);

      len = Math.max(mA.length, mB.length);
      for (int k = 0; k < len; k++) {
        int a = k < mA.length ? mA[k] : 0;
        int b = k < mB.length ? mB[k] : 0;
        if (a != b) return a - b;
      }

      // Expand to next BFS level
      List<Integer> nextA = new ArrayList<>(), nextB = new ArrayList<>();
      for (int i : frontA)
        for (int c : nodesA.get(i).children) nextA.add(c);
      for (int i : frontB)
        for (int c : nodesB.get(i).children) nextB.add(c);

      frontA = nextA;
      frontB = nextB;
    }

    // Rules 1-2 exhausted. Try Rule 4b-c if descriptor map is available.
    if (descriptorMap != null && !descriptorMap.isEmpty()) {
      return compareBranches_Rule4bc(nodesA, rootA, nodesB, rootB, descriptorMap);
    }
    return 0;
  }

  /**
   * CIP Rule 4b-c: compare branches by like/unlike descriptor pairing.
   * Collects R/S descriptors along each branch in BFS order and compares.
   * R(1) > S(2) > NONE(0) per Rule 5.
   * @since 6.5.2
   */
  private static int compareBranches_Rule4bc(
      List<DigraphNode> nodesA, int rootA,
      List<DigraphNode> nodesB, int rootB,
      Map<Integer, Integer> descriptorMap) {

    List<Integer> descA = collectDescriptors(nodesA, rootA, descriptorMap);
    List<Integer> descB = collectDescriptors(nodesB, rootB, descriptorMap);

    Collections.sort(descA);
    Collections.sort(descB);

    int len = Math.max(descA.size(), descB.size());
    for (int k = 0; k < len; k++) {
      int a = k < descA.size() ? descA.get(k) : 0;
      int b = k < descB.size() ? descB.get(k) : 0;
      if (a != b) {
        if (a == 0) return -1;
        if (b == 0) return 1;
        return (a < b) ? 1 : -1; // R(1) > S(2)
      }
    }
    return 0;
  }

  /** Collect CIP descriptors along a digraph branch in BFS order. */
  private static List<Integer> collectDescriptors(
      List<DigraphNode> nodes, int root,
      Map<Integer, Integer> descriptorMap) {
    List<Integer> result = new ArrayList<>();
    List<Integer> frontier = new ArrayList<>();
    frontier.add(root);
    while (!frontier.isEmpty()) {
      List<Integer> next = new ArrayList<>();
      for (int ni : frontier) {
        DigraphNode nd = nodes.get(ni);
        if (!nd.isDuplicate && nd.atomIdx >= 0) {
          Integer desc = descriptorMap.get(nd.atomIdx);
          if (desc != null && desc != 0) result.add(desc);
        }
        next.addAll(nd.children);
      }
      frontier = next;
    }
    return result;
  }

  /** Extract sorted (descending) array of atomic numbers or mass numbers from a frontier. */
  private static int[] toSortedDescArray(List<Integer> frontier, List<DigraphNode> nodes,
                                         boolean useAtomicNum) {
    int[] arr = new int[frontier.size()];
    for (int i = 0; i < arr.length; i++) {
      DigraphNode nd = nodes.get(frontier.get(i));
      arr[i] = useAtomicNum ? nd.atomicNum : nd.massNumber;
    }
    Arrays.sort(arr);
    // Reverse to get descending
    for (int i = 0, j = arr.length - 1; i < j; i++, j--) {
      int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }
    return arr;
  }

  // ==========================================================================
  // Priority computation
  // ==========================================================================

  /**
   * Compute CIP priorities for all ligands around a centre atom.
   * Returns an array of {atomIdx, rank} pairs where atomIdx is -1 for
   * implicit H and rank is 1-based (1 = lowest priority).
   */
  private static int[][] computePriorities(MolGraph g, int centre) {
    return computePriorities(g, centre, null);
  }

  /**
   * Compute CIP priorities with optional Rule 4b-c/5 descriptor map.
   * @since 6.5.2
   */
  private static int[][] computePriorities(MolGraph g, int centre,
                                            Map<Integer, Integer> descriptorMap) {
    int[] nbrs = g.neighbors[centre];
    int nExplicit = nbrs.length;

    int boSum = 0;
    for (int nb : nbrs) {
      int bo = g.bondOrder(centre, nb);
      if (bo == 4) bo = 1;
      boSum += bo;
    }
    int hCount = implicitH(g, centre);
    int total = nExplicit + hCount;

    // Build full digraph from centre (shared for all ligands)
    List<DigraphNode> digraph = buildDigraph(g, centre);

    // Find the digraph child node for each explicit neighbor
    int[][] result = new int[total][2]; // {atomIdx, rank}

    // Map each explicit neighbor to its digraph root child
    int[] ligandRoots = new int[total];
    for (int i = 0; i < nExplicit; i++) {
      result[i][0] = nbrs[i];
      ligandRoots[i] = -1;
      // Find the real (non-duplicate) child in the digraph
      for (int c : digraph.get(0).children) {
        if (digraph.get(c).atomIdx == nbrs[i] && !digraph.get(c).isDuplicate) {
          ligandRoots[i] = c;
          break;
        }
      }
      if (ligandRoots[i] < 0) {
        // Fallback: use any child matching
        for (int c : digraph.get(0).children) {
          if (digraph.get(c).atomIdx == nbrs[i]) {
            ligandRoots[i] = c;
            break;
          }
        }
      }
    }

    // Implicit H ligands: atomIdx=-1, trivial single-node digraph
    for (int i = 0; i < hCount; i++) {
      result[nExplicit + i][0] = -1;
      ligandRoots[nExplicit + i] = -2; // sentinel for H
    }

    // Sort ligand indices by CIP priority
    Integer[] order = new Integer[total];
    for (int i = 0; i < total; i++) order[i] = i;
    // Capture descriptorMap for use in lambda
    final Map<Integer, Integer> dMap = descriptorMap;

    // Comparison helper that uses full rules when descriptor map is available
    java.util.function.BiFunction<Integer, Integer, Integer> cmpLigands = (a, b) -> {
      int zA = (a < nExplicit) ? g.atomicNum[nbrs[a]] : 1;
      int zB = (b < nExplicit) ? g.atomicNum[nbrs[b]] : 1;
      if (zA != zB) return Integer.compare(zA, zB);

      int rootA = ligandRoots[a];
      int rootB = ligandRoots[b];
      if (rootA == -2 && rootB == -2) return 0;
      if (rootA == -2) return -1;
      if (rootB == -2) return 1;
      if (rootA >= 0 && rootB >= 0)
        return compareBranches(digraph, rootA, digraph, rootB, dMap);
      return 0;
    };

    Arrays.sort(order, (a, b) -> cmpLigands.apply(a, b));

    // Assign ranks
    int rank = 1;
    for (int i = 0; i < total; i++) {
      if (i > 0) {
        boolean same = (cmpLigands.apply(order[i - 1], order[i]) == 0);
        if (!same) rank = i + 1;
      }
      result[order[i]][1] = rank;
    }

    return result;
  }

  // ==========================================================================
  // R/S Assignment
  // ==========================================================================

  /**
   * Assign CIP R/S descriptors to all tetrahedral stereocentres in the molecule.
   *
   * Uses a two-pass approach per IUPAC 2013 P-92.1.4:
   *   Pass 1: Rules 1-2 (atomic number, mass number)
   *   Pass 2: Rules 4b-c (like/unlike descriptor pairing) + Rule 5 (R > S)
   *
   * Centres resolved only in pass 2 are pseudoasymmetric and receive
   * lowercase 'r' or 's'.
   *
   * @param g the molecular graph
   * @return map from atom index to 'R', 'S', 'r', or 's'; atoms that are not valid
   *         stereocentres or where priorities cannot be resolved are omitted
   * @since 6.5.2
   */
  public static Map<Integer, Character> assignRS(MolGraph g) {
    Map<Integer, Character> result = new LinkedHashMap<>();
    if (g.tetraChirality == null) return result;

    // Pass 1: Rules 1-2 only
    for (int idx = 0; idx < g.n; idx++) {
      Character c = assignRSAtom(g, idx, null);
      if (c != null) result.put(idx, c);
    }

    // Build descriptor map from pass 1 for Rule 4b-c
    Map<Integer, Integer> descriptorMap = new HashMap<>();
    for (Map.Entry<Integer, Character> e : result.entrySet()) {
      char ch = e.getValue();
      descriptorMap.put(e.getKey(), (ch == 'R' || ch == 'r') ? 1 : 2);
    }

    // Pass 2: retry NONE centres with Rules 4b-c + 5
    if (!descriptorMap.isEmpty()) {
      for (int idx = 0; idx < g.n; idx++) {
        if (g.tetraChirality[idx] == 0) continue;
        if (result.containsKey(idx)) continue; // already assigned in pass 1
        Character c = assignRSAtom(g, idx, descriptorMap);
        if (c != null) {
          // Pseudoasymmetric: resolved only by Rule 4b-c/5 → lowercase
          result.put(idx, Character.toLowerCase(c));
        }
      }
    }

    return result;
  }

  /**
   * Assign R/S to a single atom with optional descriptor map for Rule 4b-c/5.
   * @return 'R' or 'S', or null if not a valid stereocentre.
   */
  private static Character assignRSAtom(MolGraph g, int idx,
                                         Map<Integer, Integer> descriptorMap) {
    if (g.tetraChirality[idx] == 0) return null;

    int deg = g.neighbors[idx].length;
    int hc = implicitH(g, idx);
    int totalSubs = deg + hc;
    if (totalSubs != 4) return null;

    int[][] priorities = computePriorities(g, idx, descriptorMap);
    if (priorities.length != 4) return null;

    Set<Integer> ranks = new HashSet<>();
    for (int[] p : priorities) ranks.add(p[1]);
    if (ranks.size() != 4) return null;

    int[] smilesOrder;
    if (hc > 0) {
      smilesOrder = new int[4];
      smilesOrder[0] = deg;
      for (int i = 0; i < deg; i++) smilesOrder[i + 1] = i;
    } else {
      smilesOrder = new int[]{0, 1, 2, 3};
    }

    int[] cipRanks = new int[4];
    for (int k = 0; k < 4; k++) {
      cipRanks[k] = priorities[smilesOrder[k]][1];
    }

    int inversions = 0;
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++) {
        if (cipRanks[i] > cipRanks[j]) inversions++;
      }
    }
    boolean evenPerm = (inversions % 2 == 0);
    boolean smilesIsACW = (g.tetraChirality[idx] == 1);

    if (evenPerm) {
      return smilesIsACW ? 'S' : 'R';
    } else {
      return smilesIsACW ? 'R' : 'S';
    }
  }

  // ==========================================================================
  // E/Z Assignment
  // ==========================================================================

  /**
   * Assign CIP E/Z descriptors to all double bonds with stereo annotations.
   *
   * @param g the molecular graph
   * @return map from bond key (as long: bondKey(i,j) with i &lt; j) to 'E' or 'Z';
   *         bonds without stereo annotation or with unresolvable priorities are omitted
   */
  public static Map<Long, Character> assignEZ(MolGraph g) {
    Map<Long, Character> result = new LinkedHashMap<>();

    for (int i = 0; i < g.n; i++) {
      for (int j : g.neighbors[i]) {
        if (j <= i) continue;
        if (g.bondOrder(i, j) != 2) continue;

        int conf = g.dbStereo(i, j); // 1=Z, 2=E
        if (conf == 0) continue;

        // Determine higher-priority substituent on each side
        int highI = higherPrioritySub(g, i, j);
        int highJ = higherPrioritySub(g, j, i);
        if (highI == Integer.MIN_VALUE || highJ == Integer.MIN_VALUE) continue;

        // Check if the CIP high-priority atom matches the SMILES reference atom
        int refI = firstNonDBNeighbor(g, i, j);
        int refJ = firstNonDBNeighbor(g, j, i);

        int flips = 0;
        if (highI != refI && refI >= 0) flips++;
        if (highJ != refJ && refJ >= 0) flips++;

        int effectiveConf = conf;
        if (flips == 1) effectiveConf = (conf == 1) ? 2 : 1;

        long bk = bondKey(i, j);
        result.put(bk, (effectiveConf == 1) ? 'Z' : 'E');
      }
    }

    return result;
  }

  /** Find the higher-priority substituent on one end of a double bond. */
  private static int higherPrioritySub(MolGraph g, int thisEnd, int otherEnd) {
    int[][] prios = computePriorities(g, thisEnd);
    int bestRank = -1, bestAtom = Integer.MIN_VALUE;
    int secondRank = -1;
    int subCount = 0;

    for (int[] p : prios) {
      if (p[0] == otherEnd) continue;
      subCount++;
      if (p[1] > bestRank) {
        secondRank = bestRank;
        bestRank = p[1];
        bestAtom = p[0];
      } else if (p[1] > secondRank) {
        secondRank = p[1];
      }
    }

    // If two substituents have the same highest priority, E/Z is undefined
    if (subCount >= 2 && bestRank == secondRank) return Integer.MIN_VALUE;
    return bestAtom;
  }

  /** Get the first neighbor of thisEnd that is not the double-bond partner. */
  private static int firstNonDBNeighbor(MolGraph g, int thisEnd, int otherEnd) {
    for (int nb : g.neighbors[thisEnd]) {
      if (nb != otherEnd) return nb;
    }
    return -1;
  }

  /** Pack two atom indices into a canonical bond key (same as MolGraph). */
  private static long bondKey(int i, int j) {
    int lo = Math.min(i, j), hi = Math.max(i, j);
    return ((long) lo << 32) | (hi & 0xFFFFFFFFL);
  }
}

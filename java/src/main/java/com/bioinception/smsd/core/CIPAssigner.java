/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/**
 * CIP (Cahn-Ingold-Prelog) stereodescriptor assignment engine for MolGraph.
 *
 * Assigns R/S descriptors to tetrahedral stereocentres and E/Z descriptors
 * to double bonds using BFS digraph traversal with CIP sequence rules:
 *   Rule 1: higher atomic number gets higher priority
 *   Rule 2: higher mass number breaks ties (isotopes)
 *   Rule 3: Z double-bond configuration gets higher priority than E
 *
 * Implements the digraph (extended connectivity tree) approach with phantom
 * (duplicate) nodes for ring closures and multiple bonds. Aromatic ligands use
 * a Kekule representation; fractional mancude averaging and the complete
 * Rule 4/5 machinery are not implemented.
 *
 * References:
 *   Prelog and Helmchen, Angew. Chem. Int. Ed. Engl. 21, 567 (1982)
 *   IUPAC Recommendations 2013 (P-92)
 *
 * @author Syed Asad Rahman
 */
public final class CIPAssigner {

  /** Maximum BFS depth for the CIP digraph expansion.
   *  Must match C++ maxDepth (30) for long-chain resolution parity. */
  public static final int MAX_DEPTH = 30;

  /** Maximum number of nodes in the digraph to prevent runaway on large molecules. */
  private static final int MAX_NODES = 5000;

  private CIPAssigner() {} // utility class

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
    if (g.hasHydrogenCounts()) return g.hydrogenCount(idx);
    if (g.aromatic[idx]) {
      int element = g.atomicNum[idx];
      if (element == 8 || element == 16 || element == 34 || element == 52
          || ((element == 7 || element == 15) && g.neighbors[idx].length >= 3)) return 0;
    }
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

  private static boolean hasLonePairLigand(MolGraph graph, int atom, int hydrogens) {
    if (hydrogens != 0 || graph.neighbors[atom].length != 3 || graph.tetraChirality == null
        || graph.tetraChirality[atom] == 0) return false;
    int element = graph.atomicNum[atom];
    return element == 7 || element == 15 || element == 16 || element == 33
        || element == 34 || element == 51 || element == 52;
  }

  /** Prepare aromatic bond orders without changing the caller's molecule. */
  private static IAtomContainer kekuleCopy(MolGraph graph) {
    try {
      IAtomContainer copy;
      if (graph.mol != null) {
        copy = graph.mol.clone();
      } else {
        copy = SilentChemObjectBuilder.getInstance().newAtomContainer();
        for (int atom = 0; atom < graph.n; atom++) {
          IAtom value = SilentChemObjectBuilder.getInstance().newAtom();
          value.setAtomicNumber(graph.atomicNum[atom]);
          value.setFormalCharge(graph.formalCharge[atom]);
          value.setImplicitHydrogenCount(implicitH(graph, atom));
          value.setIsAromatic(graph.aromatic[atom]);
          copy.addAtom(value);
        }
        for (int atom = 0; atom < graph.n; atom++) {
          for (int neighbor : graph.neighbors[atom]) {
            if (neighbor <= atom) continue;
            int order = graph.bondOrder(atom, neighbor);
            copy.addBond(atom, neighbor, order == 2 ? IBond.Order.DOUBLE
                : order == 3 ? IBond.Order.TRIPLE : IBond.Order.SINGLE);
            copy.getBond(copy.getBondCount() - 1).setIsAromatic(order == 4);
          }
        }
      }
      Kekulization.kekulize(copy);
      return copy;
    } catch (Exception unavailable) {
      return null;
    }
  }

  /** Compare spheres in parent order, padding absent phantom children with zero. */
  private static int compareOrderedBranches(List<DigraphNode> firstNodes, int firstRoot,
                                            List<DigraphNode> secondNodes, int secondRoot,
                                            boolean compareMass) {
    List<int[]> frontier = new ArrayList<>();
    frontier.add(new int[] {firstRoot, secondRoot});
    while (!frontier.isEmpty()) {
      for (int[] pair : frontier) {
        int firstValue = pair[0] < 0 ? 0 : (compareMass
            ? firstNodes.get(pair[0]).massNumber : firstNodes.get(pair[0]).atomicNum);
        int secondValue = pair[1] < 0 ? 0 : (compareMass
            ? secondNodes.get(pair[1]).massNumber : secondNodes.get(pair[1]).atomicNum);
        if (firstValue != secondValue) return Integer.compare(firstValue, secondValue);
      }
      List<int[]> next = new ArrayList<>();
      for (int[] pair : frontier) {
        List<Integer> firstChildren = pair[0] < 0 ? Collections.emptyList() : firstNodes.get(pair[0]).children;
        List<Integer> secondChildren = pair[1] < 0 ? Collections.emptyList() : secondNodes.get(pair[1]).children;
        for (int slot = 0; slot < Math.max(firstChildren.size(), secondChildren.size()); slot++)
          next.add(new int[] {slot < firstChildren.size() ? firstChildren.get(slot) : -1,
                             slot < secondChildren.size() ? secondChildren.get(slot) : -1});
      }
      frontier = next;
    }
    return 0;
  }

  /** Rank children from the leaves while preserving all earlier atomic-number decisions. */
  private static void orderDigraph(List<DigraphNode> nodes) {
    for (int node = nodes.size() - 1; node >= 0; node--) {
      nodes.get(node).children.sort((first, second) -> {
        int atomicComparison = compareOrderedBranches(nodes, first, nodes, second, false);
        return atomicComparison != 0 ? -atomicComparison
            : -compareOrderedBranches(nodes, first, nodes, second, true);
      });
    }
  }

  /**
   * Build the CIP digraph (extended connectivity tree) rooted at the
   * specified atom. Ring closures produce duplicate (phantom) nodes.
   * Double and triple bonds produce additional phantom copies.
   */
  private static List<DigraphNode> buildDigraph(MolGraph g, int root) {
    boolean aromatic = false;
    for (boolean flag : g.aromatic) aromatic |= flag;
    IAtomContainer kekule = aromatic ? kekuleCopy(g) : null;
    if (aromatic && kekule == null) return Collections.emptyList();
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
          if (bo == 4) bo = kekule.getBond(kekule.getAtom(atom), kekule.getAtom(nbr)).getOrder().numeric();

          expansions.add(new int[]{nbr, 0}); // real child
          for (int extra = 1; extra < bo; ++extra) {
            expansions.add(new int[]{nbr, 1}); // phantom for bond order
          }
        }

        // Phantoms for bond back to parent (higher-order bonds)
        if (parentAtom >= 0 && !nodes.get(ni).isDuplicate) {
          int bo = g.bondOrder(atom, parentAtom);
          if (bo == 4) bo = kekule.getBond(kekule.getAtom(atom), kekule.getAtom(parentAtom)).getOrder().numeric();
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
        for (int hydrogen = 0; hydrogen < implicitH(g, atom) && nodes.size() < MAX_NODES; hydrogen++) {
          DigraphNode child = new DigraphNode();
          child.atomIdx = -1;
          child.atomicNum = 1;
          child.massNumber = 0;
          child.parent = ni;
          child.depth = depth + 1;
          nodes.get(ni).children.add(nodes.size());
          nodes.add(child);
        }
      }
      frontier = nextFrontier;
    }

    orderDigraph(nodes);
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
    return compareBranches(nodesA, rootA, nodesB, rootB, null, null);
  }

  /**
   * Compare two digraph sub-trees by CIP Rules 1-2, with optional Rule 3 (Z > E),
   * Rule 4b-c (like/unlike descriptor pairing) and Rule 5 (R > S tiebreak).
   *
   * @param descriptorMap  atom index → descriptor (1=R, 2=S); null for Rules 1-2 only
   * @param molGraph       molecular graph for Rule 3 E/Z lookup; null to skip Rule 3
   * @since 6.5.2
   */
  private static int compareBranches(List<DigraphNode> nodesA, int rootA,
                                     List<DigraphNode> nodesB, int rootB,
                                     Map<Integer, Integer> descriptorMap,
                                     MolGraph molGraph) {
    int atomicComparison = compareOrderedBranches(nodesA, rootA, nodesB, rootB, false);
    if (atomicComparison != 0) return atomicComparison;
    int isotopeComparison = compareOrderedBranches(nodesA, rootA, nodesB, rootB, true);
    if (isotopeComparison != 0) return isotopeComparison;
    if (molGraph != null) {
      List<Integer> frontA = new ArrayList<>(), frontB = new ArrayList<>();
      frontA.add(rootA);
      frontB.add(rootB);
      while (!frontA.isEmpty() || !frontB.isEmpty()) {
        int[] firstDescriptors = ezDescriptors(frontA, nodesA, molGraph);
        int[] secondDescriptors = ezDescriptors(frontB, nodesB, molGraph);
        for (int position = 0; position < Math.max(firstDescriptors.length, secondDescriptors.length); position++) {
          int first = position < firstDescriptors.length ? firstDescriptors[position] : 0;
          int second = position < secondDescriptors.length ? secondDescriptors[position] : 0;
          if (first != second) return Integer.compare(first, second);
        }
        List<Integer> nextA = new ArrayList<>(), nextB = new ArrayList<>();
        for (int node : frontA) nextA.addAll(nodesA.get(node).children);
        for (int node : frontB) nextB.addAll(nodesB.get(node).children);
        frontA = nextA;
        frontB = nextB;
      }
    }

    // Rules 1-3 exhausted. Try Rule 4b-c if descriptor map is available.
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

  /**
   * Collect E/Z descriptors for atoms at a BFS frontier, sorted descending.
   * Maps dbStereo values so that Z gets higher priority: Z(1)→2, E(2)→1, NONE→0.
   */
  private static int[] ezDescriptors(List<Integer> frontier, List<DigraphNode> nodes,
                                     MolGraph g) {
    int[] arr = new int[frontier.size()];
    for (int i = 0; i < arr.length; i++) {
      DigraphNode nd = nodes.get(frontier.get(i));
      int parentIdx = nd.parent;
      int ez = 0;
      if (parentIdx >= 0 && nd.atomIdx >= 0 && !nd.isDuplicate) {
        int parentAtom = nodes.get(parentIdx).atomIdx;
        if (parentAtom >= 0 && g.bondOrder(parentAtom, nd.atomIdx) == 2) {
          int conf = atomicPriorityDoubleBondConfiguration(g, parentAtom, nd.atomIdx);
          ez = (conf == 1) ? 2 : (conf == 2) ? 1 : 0;   // Z→2, E→1, NONE→0
        }
      }
      arr[i] = ez;
    }
    Arrays.sort(arr);
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
    return computePriorities(g, centre, descriptorMap, true);
  }

  private static int[][] computePriorities(MolGraph g, int centre,
                                            Map<Integer, Integer> descriptorMap, boolean includeStereo) {
    int[] nbrs = g.neighbors[centre];
    int nExplicit = nbrs.length;

    int boSum = 0;
    for (int nb : nbrs) {
      int bo = g.bondOrder(centre, nb);
      if (bo == 4) bo = 1;
      boSum += bo;
    }
    int hCount = implicitH(g, centre);
    boolean lonePair = hasLonePairLigand(g, centre, hCount);
    int total = nExplicit + hCount + (lonePair ? 1 : 0);

    // Build full digraph from centre (shared for all ligands)
    List<DigraphNode> digraph = buildDigraph(g, centre);
    if (digraph.isEmpty()) return new int[0][0];

    // Find the digraph child node for each explicit neighbour
    int[][] result = new int[total][2]; // {atomIdx, rank}

    // Map each explicit neighbour to its digraph root child
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

    for (int implicit = 0; implicit < hCount + (lonePair ? 1 : 0); implicit++) {
      boolean phantom = implicit >= hCount;
      DigraphNode node = new DigraphNode();
      node.atomIdx = phantom ? -2 : -1;
      node.atomicNum = phantom ? 0 : 1;
      node.parent = -1;
      result[nExplicit + implicit][0] = node.atomIdx;
      ligandRoots[nExplicit + implicit] = digraph.size();
      digraph.add(node);
    }

    // Sort ligand indices by CIP priority
    Integer[] order = new Integer[total];
    for (int i = 0; i < total; i++) order[i] = i;
    // Capture descriptorMap for use in lambda
    final Map<Integer, Integer> dMap = descriptorMap;

    // Comparison helper that uses full rules when descriptor map is available
    java.util.function.BiFunction<Integer, Integer, Integer> cmpLigands = (a, b) -> {
      int zA = result[a][0] >= 0 ? g.atomicNum[result[a][0]] : result[a][0] == -2 ? 0 : 1;
      int zB = result[b][0] >= 0 ? g.atomicNum[result[b][0]] : result[b][0] == -2 ? 0 : 1;
      if (zA != zB) return Integer.compare(zA, zB);

      int rootA = ligandRoots[a];
      int rootB = ligandRoots[b];
      if (rootA >= 0 && rootB >= 0)
        return compareBranches(digraph, rootA, digraph, rootB, dMap, includeStereo ? g : null);
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
    boolean lonePair = hasLonePairLigand(g, idx, hc);
    int totalSubs = deg + hc + (lonePair ? 1 : 0);
    if (totalSubs != 4) return null;

    int[][] priorities = computePriorities(g, idx, descriptorMap);
    if (priorities.length != 4) return null;

    Set<Integer> ranks = new HashSet<>();
    for (int[] p : priorities) ranks.add(p[1]);
    if (ranks.size() != 4) return null;

    int[] cipRanks = new int[4];
    boolean cdkOrder = g.tetraLigands != null && g.tetraLigands[idx] != null;
    if (cdkOrder) {
      int[] ligandOrder = g.tetraLigands[idx];
      if (ligandOrder.length != 4) return null;
      for (int k = 0; k < 4; k++) {
        boolean found = false;
        for (int[] priority : priorities) {
          int ligand = lonePair && ligandOrder[k] < 0 ? -2 : ligandOrder[k];
          if (priority[0] == ligand) {
            cipRanks[k] = priority[1];
            found = true;
            break;
          }
        }
        if (!found) return null;
      }
    } else {
      // Builder annotations retain the original implicit-H-first convention.
      for (int k = 0; k < 4; k++) {
        int ligand = hc > 0 || lonePair ? (k == 0 ? deg : k - 1) : k;
        cipRanks[k] = priorities[ligand][1];
      }
    }

    int inversions = 0;
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++) {
        if (cipRanks[i] > cipRanks[j]) inversions++;
      }
    }
    boolean evenPerm = (inversions % 2 == 0);
    boolean smilesIsACW = (g.tetraChirality[idx] == 1) != cdkOrder;

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

        // The stored relation uses the first non-double-bond neighbour at each end.
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
    return higherPrioritySub(g, thisEnd, otherEnd, true);
  }

  private static int higherPrioritySub(MolGraph g, int thisEnd, int otherEnd, boolean includeStereo) {
    int[][] prios = computePriorities(g, thisEnd, null, includeStereo);
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

  /** Rule 3 descriptors use atomic/isotopic priorities, independent of input carrier order. */
  private static int atomicPriorityDoubleBondConfiguration(MolGraph graph, int first, int second) {
    int configuration = graph.dbStereo(first, second);
    if (configuration == 0) return 0;
    int highFirst = higherPrioritySub(graph, first, second, false);
    int highSecond = higherPrioritySub(graph, second, first, false);
    if (highFirst == Integer.MIN_VALUE || highSecond == Integer.MIN_VALUE) return 0;
    if (highFirst != firstNonDBNeighbor(graph, first, second)) configuration = 3 - configuration;
    if (highSecond != firstNonDBNeighbor(graph, second, first)) configuration = 3 - configuration;
    return configuration;
  }

  /** Get the first neighbour of thisEnd that is not the double-bond partner. */
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

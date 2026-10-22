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
 *   Rule 1a: higher atomic number gets higher priority
 *   Rule 1b: a ring-closure duplicate nearer the root gets higher priority
 *   Rule 2: higher mass number breaks ties (isotopes); an unlabelled atom
 *   counts at its standard atomic weight, a duplicate at none
 *   Rule 3: Z double-bond configuration gets higher priority than E
 *   Rules 4 and 5: auxiliary descriptors, like and unlike descriptor pairs,
 *   r before s and R before S; a centre resolved by Rule 5 is pseudoasymmetric
 *
 * A first pass ranks the ligands by Rule 1a alone; ligands it leaves tied are
 * ranked over the hierarchical digraph by every rule in turn, each read
 * exhaustively before the next ({@link SequenceRules}).
 *
 * Implements the digraph (extended connectivity tree) approach with phantom
 * (duplicate) nodes for ring closures and multiple bonds. The duplicate atom of
 * a double bond in a mancude ring system takes the mean atomic number of the
 * ring neighbours it may be double bonded to (P-92.1.4.4), so labels do not
 * depend on the Kekule structure.
 *
 * References:
 *   Prelog and Helmchen, Angew. Chem. Int. Ed. Engl. 21, 567 (1982)
 *   IUPAC Recommendations 2013 (P-92)
 *   Hanson et al., JCIM 58, 1755 (2018), the auxiliary descriptor method the
 *   reference implementations follow
 *
 * @author Syed Asad Rahman
 */
public final class CIPAssigner {

  /** Maximum BFS depth for the CIP digraph expansion.
   *  Must match C++ maxDepth (30) for long-chain resolution parity. */
  public static final int MAX_DEPTH = 30;

  /** Maximum number of nodes in the digraph to prevent runaway on large molecules. */
  private static final int MAX_NODES = 5000;

  /** Scale of DigraphNode.scaledAtomicNum, so the mean of up to four atomic numbers is exact. */
  private static final int ATOMIC_NUMBER_SCALE = 12;

  private CIPAssigner() {} // utility class

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
    int scaledAtomicNum; // Rule 1 key: atomicNum * ATOMIC_NUMBER_SCALE, a mean for mancude duplicates
    boolean isRingDuplicate; // a duplicate closing a ring, rather than one of a multiple bond
    int distance;      // Rule 1b: root 1, each sphere one more; a duplicate takes its original's
    final List<Integer> children = new ArrayList<>(4);
  }

  /** Whether a tagged centre has a lone pair as its fourth ligand, see {@link MolGraph#tetrahedralLigands}. */
  private static boolean hasLonePairLigand(MolGraph graph, int atom, int hydrogens) {
    return hydrogens == 0 && graph.neighbors[atom].length == 3 && graph.tetraChirality != null
        && graph.tetraChirality[atom] != 0 && MolGraph.lonePairCentre(graph.atomicNum[atom]);
  }

  /** Prepare aromatic bond orders without changing the caller's molecule. */
  private static IAtomContainer kekuleCopy(MolGraph graph) {
    try {
      IAtomContainer copy;
      if (graph.mol != null) {
        copy = graph.mol.clone();
        // The graph's hydrogens, which differ for a ring system CDK read from ':' bonds.
        for (int atom = 0; atom < graph.n; atom++) {
          Integer supplied = copy.getAtom(atom).getImplicitHydrogenCount();
          if (supplied != null && supplied != graph.hydrogenCount(atom))
            copy.getAtom(atom).setImplicitHydrogenCount(graph.hydrogenCount(atom));
        }
      } else {
        copy = SilentChemObjectBuilder.getInstance().newAtomContainer();
        for (int atom = 0; atom < graph.n; atom++) {
          IAtom value = SilentChemObjectBuilder.getInstance().newAtom();
          value.setAtomicNumber(graph.atomicNum[atom]);
          value.setFormalCharge(graph.formalCharge[atom]);
          value.setImplicitHydrogenCount(graph.hydrogenCount(atom));
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

  /**
   * Sequence Rule 2 between two nodes Rule 1a has paired, either of which may be absent (null). Two
   * atoms without a mass number tie whatever their weights; otherwise a labelled atom counts its
   * mass number, an unlabelled atom its standard atomic weight and a duplicate nothing, as the
   * reference CIP implementations read the rule (Hanson et al., JCIM 58, 1755 (2018)). Both passes
   * use it.
   */
  private static double rule2Mass(DigraphNode node) {
    if (node.isDuplicate || node.atomicNum <= 0) return 0.0;
    return node.massNumber != 0 ? node.massNumber : atomicWeight(node.atomicNum);
  }

  private static int compareMassNumbers(DigraphNode first, DigraphNode second) {
    int firstLabel = first == null ? 0 : first.massNumber;
    int secondLabel = second == null ? 0 : second.massNumber;
    if (firstLabel == 0 && secondLabel == 0) return 0;
    double firstMass = first == null ? 0.0 : rule2Mass(first);
    double secondMass = second == null ? 0.0 : rule2Mass(second);
    return Double.compare(firstMass, secondMass);
  }

  /** Compare spheres in parent order, padding absent phantom children with zero. */
  private static int compareOrderedBranches(List<DigraphNode> firstNodes, int firstRoot,
                                            List<DigraphNode> secondNodes, int secondRoot,
                                            boolean compareMass) {
    List<int[]> frontier = new ArrayList<>();
    frontier.add(new int[] {firstRoot, secondRoot});
    while (!frontier.isEmpty()) {
      for (int[] pair : frontier) {
        if (compareMass) {
          int massComparison = compareMassNumbers(pair[0] < 0 ? null : firstNodes.get(pair[0]),
                                                  pair[1] < 0 ? null : secondNodes.get(pair[1]));
          if (massComparison != 0) return massComparison;
          continue;
        }
        int firstValue = pair[0] < 0 ? 0 : firstNodes.get(pair[0]).scaledAtomicNum;
        int secondValue = pair[1] < 0 ? 0 : secondNodes.get(pair[1]).scaledAtomicNum;
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

  /** Whether the bond is a ring bond, or an aromatic bond as Kekulization reads it. */
  private static boolean ringBond(MolGraph g, int atom, int neighbor) {
    return g.bondInRing(atom, neighbor) || (g.aromatic[atom] && g.aromatic[neighbor]
        && (g.bondOrder(atom, neighbor) == 4 || g.bondAromatic(atom, neighbor)));
  }

  /** Bond order with aromatic bonds read from the Kekule copy; 4 when that copy has none. */
  private static int kekuleOrder(MolGraph g, IAtomContainer kekule, int atom, int neighbor) {
    int order = g.bondOrder(atom, neighbor);
    if (order != 4 || kekule == null) return order;
    IBond bond = kekule.getBond(kekule.getAtom(atom), kekule.getAtom(neighbor));
    return bond == null || bond.getOrder() == null ? 4 : bond.getOrder().numeric();
  }

  /**
   * The duplicate a ring atom with exactly one double bond, partners[atom] in the Kekule bond
   * orders, carries for that bond (-1 for other atoms). A mancude ring system is what remains of
   * those atoms after repeatedly removing the ones with fewer than two others as ring neighbours.
   * There, when the partner is in the system too, the duplicate takes the mean atomic number of
   * the atom's ring neighbours in the system whenever their sum differs from the partner's atomic
   * number times their count (IUPAC 2013 P-92.1.4.4), so labels do not depend on the Kekule
   * structure even when the mean is a whole number, as for a boron and a nitrogen neighbour;
   * otherwise it is the partner. numbers[atom] is that atomic number times ATOMIC_NUMBER_SCALE.
   * The native core applies the same rule.
   */
  private static void mancudeDuplicates(MolGraph g, IAtomContainer kekule, int[] numbers, int[] partners) {
    Arrays.fill(partners, -1);
    for (int atom = 0; atom < g.n; atom++) {
      int partner = -1, doubles = 0, triples = 0;
      for (int neighbor : g.neighbors[atom]) {
        int order = kekuleOrder(g, kekule, atom, neighbor);
        if (order == 2) {
          partner = neighbor;
          doubles++;
        } else if (order == 3) {
          triples++;
        }
      }
      if ((g.ring[atom] || g.aromatic[atom]) && doubles == 1 && triples == 0) partners[atom] = partner;
    }
    boolean[] mancude = new boolean[g.n];
    int[] degree = new int[g.n];
    Deque<Integer> pending = new ArrayDeque<>();
    for (int atom = 0; atom < g.n; atom++) {
      if (partners[atom] < 0) continue;
      mancude[atom] = true;
      for (int neighbor : g.neighbors[atom])
        if (partners[neighbor] >= 0 && ringBond(g, atom, neighbor)) degree[atom]++;
      if (degree[atom] < 2) pending.push(atom);
    }
    while (!pending.isEmpty()) {
      int atom = pending.pop();
      if (!mancude[atom]) continue;
      mancude[atom] = false;
      for (int neighbor : g.neighbors[atom])
        if (mancude[neighbor] && ringBond(g, atom, neighbor) && --degree[neighbor] < 2)
          pending.push(neighbor);
    }
    for (int atom = 0; atom < g.n; atom++) {
      if (partners[atom] < 0) continue;
      int sum = 0, count = 0;
      if (mancude[atom]) {
        for (int neighbor : g.neighbors[atom]) {
          if (!mancude[neighbor] || !ringBond(g, atom, neighbor)) continue;
          sum += g.atomicNum[neighbor];
          count++;
        }
      }
      int partnerNumber = g.atomicNum[partners[atom]];
      numbers[atom] = count == 0 || !mancude[partners[atom]] || sum == partnerNumber * count
          ? partnerNumber * ATOMIC_NUMBER_SCALE : sum * ATOMIC_NUMBER_SCALE / count;
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
    int[] duplicateNumbers = new int[g.n], duplicatePartners = new int[g.n];
    mancudeDuplicates(g, kekule, duplicateNumbers, duplicatePartners);
    List<DigraphNode> nodes = new ArrayList<>(128);

    DigraphNode rootNode = new DigraphNode();
    rootNode.atomIdx = root;
    rootNode.atomicNum = g.atomicNum[root];
    rootNode.scaledAtomicNum = g.atomicNum[root] * ATOMIC_NUMBER_SCALE;
    rootNode.massNumber = g.massNumber[root];
    rootNode.isDuplicate = false;
    rootNode.parent = -1;
    rootNode.depth = 0;
    rootNode.distance = 1;
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
        List<int[]> expansions = new ArrayList<>(); // neighborAtom, isDuplicate, scaledAtomicNum

        for (int nbr : g.neighbors[atom]) {
          if (nbr == parentAtom && !nodes.get(ni).isDuplicate) continue;

          // The partner's double bond is the mancude duplicate below.
          int bo = nbr == duplicatePartners[atom] ? 1 : kekuleOrder(g, kekule, atom, nbr);
          if (bo == 4) bo = 1;
          int scaled = g.atomicNum[nbr] * ATOMIC_NUMBER_SCALE;

          expansions.add(new int[]{nbr, 0, scaled}); // real child
          for (int extra = 1; extra < bo; ++extra) {
            expansions.add(new int[]{nbr, 1, scaled}); // phantom for bond order
          }
        }

        // Phantoms for bond back to parent (higher-order bonds)
        if (parentAtom >= 0 && !nodes.get(ni).isDuplicate) {
          int bo = parentAtom == duplicatePartners[atom] ? 1 : kekuleOrder(g, kekule, atom, parentAtom);
          if (bo == 4) bo = 1;
          for (int extra = 1; extra < bo; ++extra) {
            expansions.add(new int[]{parentAtom, 1, g.atomicNum[parentAtom] * ATOMIC_NUMBER_SCALE});
          }
        }
        if (duplicatePartners[atom] >= 0)
          expansions.add(new int[]{duplicatePartners[atom], 1, duplicateNumbers[atom]});

        for (int[] exp : expansions) {
          if (nodes.size() >= MAX_NODES) break;
          int nbrAtom = exp[0];
          boolean forceDuplicate = (exp[1] == 1);

          // Check if this atom is on the path from root to current node. A duplicate
          // of an atom on the path takes that node's distance; one of an atom not
          // yet reached has none.
          int pathNode = -1;
          for (int walk = ni; walk >= 0; walk = nodes.get(walk).parent) {
            if (nodes.get(walk).atomIdx == nbrAtom) {
              pathNode = walk;
              break;
            }
          }
          boolean onPath = !forceDuplicate && pathNode >= 0;

          boolean makeDuplicate = forceDuplicate || onPath;

          DigraphNode child = new DigraphNode();
          child.atomIdx = nbrAtom;
          child.atomicNum = g.atomicNum[nbrAtom];
          child.scaledAtomicNum = exp[2];
          child.massNumber = makeDuplicate ? 0 : g.massNumber[nbrAtom];
          child.isDuplicate = makeDuplicate;
          child.isRingDuplicate = onPath;
          child.parent = ni;
          child.depth = depth + 1;
          child.distance = makeDuplicate ? (pathNode >= 0 ? nodes.get(pathNode).distance : 0) : depth + 2;

          int childIdx = nodes.size();
          nodes.add(child);
          nodes.get(ni).children.add(childIdx);

          if (!makeDuplicate) {
            nextFrontier.add(childIdx);
          }
        }
        for (int hydrogen = 0; hydrogen < g.hydrogenCount(atom) && nodes.size() < MAX_NODES; hydrogen++) {
          DigraphNode child = new DigraphNode();
          child.atomIdx = -1;
          child.atomicNum = 1;
          child.scaledAtomicNum = ATOMIC_NUMBER_SCALE;
          child.massNumber = 0;
          child.parent = ni;
          child.depth = depth + 1;
          child.distance = depth + 2;
          nodes.get(ni).children.add(nodes.size());
          nodes.add(child);
        }
      }
      frontier = nextFrontier;
    }

    orderDigraph(nodes);
    return nodes;
  }

  // ==========================================================================
  // Priority computation
  // ==========================================================================

  /**
   * The first ranking pass: the ligands of a centre by Sequence Rule 1a alone, the atomic numbers
   * of their branches compared sphere by sphere in the hierarchical order of the digraph. Every
   * later rule (1b, 2, 3, 4 and 5) is left to the hierarchical comparison of {@link SequenceRules},
   * which reads each rule exhaustively before the next.
   * Returns an array of {atomIdx, rank} pairs where atomIdx is -1 for implicit H and -2 for a lone
   * pair, and rank is 1-based (1 = lowest priority); ligands Rule 1a does not tell apart share a rank.
   */
  private static int[][] computePriorities(MolGraph g, int centre) {
    int[] nbrs = g.neighbors[centre];
    int nExplicit = nbrs.length;

    int boSum = 0;
    for (int nb : nbrs) {
      int bo = g.bondOrder(centre, nb);
      if (bo == 4) bo = 1;
      boSum += bo;
    }
    int hCount = g.hydrogenCount(centre);
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
      node.scaledAtomicNum = node.atomicNum * ATOMIC_NUMBER_SCALE;
      node.parent = -1;
      result[nExplicit + implicit][0] = node.atomIdx;
      ligandRoots[nExplicit + implicit] = digraph.size();
      digraph.add(node);
    }

    // Sort ligand indices by CIP priority
    Integer[] order = new Integer[total];
    for (int i = 0; i < total; i++) order[i] = i;

    java.util.function.BiFunction<Integer, Integer, Integer> cmpLigands = (a, b) -> {
      int zA = result[a][0] >= 0 ? g.atomicNum[result[a][0]] : result[a][0] == -2 ? 0 : 1;
      int zB = result[b][0] >= 0 ? g.atomicNum[result[b][0]] : result[b][0] == -2 ? 0 : 1;
      if (zA != zB) return Integer.compare(zA, zB);

      int rootA = ligandRoots[a];
      int rootB = ligandRoots[b];
      if (rootA >= 0 && rootB >= 0)
        return compareOrderedBranches(digraph, rootA, digraph, rootB, false);
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
   *   Pass 1: Rule 1a (atomic number), sphere by sphere
   *   Pass 2: for the centres pass 1 leaves unresolved, Rules 1a to 6 over the hierarchical
   *           digraph, with auxiliary descriptors for Rules 4 and 5 ({@link SequenceRules}),
   *           each rule read exhaustively before the next
   *
   * Centres resolved only by Rule 5 (or Rule 6) are pseudoasymmetric and
   * receive lowercase 'r' or 's'.
   *
   * @param g the molecular graph
   * @return map from atom index to 'R', 'S', 'r', or 's'; atoms that are not valid
   *         stereocentres or where priorities cannot be resolved are omitted
   */
  public static Map<Integer, Character> assignRS(MolGraph g) {
    Map<Integer, Character> result = new LinkedHashMap<>();
    if (g.tetraChirality == null) return result;

    // Pass 1: Rule 1a only
    for (int idx = 0; idx < g.n; idx++) {
      Character c = assignRSAtom(g, idx);
      if (c != null) result.put(idx, c);
    }

    // Pass 2: every rule, hierarchically, for the centres still unresolved
    for (int idx = 0; idx < g.n; idx++) {
      if (g.tetraChirality[idx] == 0 || result.containsKey(idx)) continue;
      Character c = new SequenceRules(g, idx).label();
      if (c != null) result.put(idx, c);
    }

    return result;
  }

  /**
   * Assign R/S to a single atom by Rule 1a.
   * @return 'R' or 'S', or null if not a valid stereocentre or the rules leave it unresolved.
   */
  private static Character assignRSAtom(MolGraph g, int idx) {
    if (g.tetraChirality[idx] == 0) return null;

    int deg = g.neighbors[idx].length;
    int hc = g.hydrogenCount(idx);
    boolean lonePair = hasLonePairLigand(g, idx, hc);
    int totalSubs = deg + hc + (lonePair ? 1 : 0);
    if (totalSubs != 4) return null;

    int[][] priorities = computePriorities(g, idx);
    if (priorities.length != 4) return null;

    Set<Integer> ranks = new HashSet<>();
    for (int[] p : priorities) ranks.add(p[1]);
    if (ranks.size() != 4) return null;

    // The ligand order is the one the stereo checks of the searches read.
    int[] cipRanks = new int[4];
    int[] ligandOrder = g.tetrahedralLigands(idx);
    if (ligandOrder == null || ligandOrder.length != 4) return null;
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

    int inversions = 0;
    for (int i = 0; i < 4; i++) {
      for (int j = i + 1; j < 4; j++) {
        if (cipRanks[i] > cipRanks[j]) inversions++;
      }
    }
    boolean evenPerm = (inversions % 2 == 0);
    boolean smilesIsACW = g.tetraChirality[idx] == 1;

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

        // Determine higher-priority substituent on each side by Rule 1a; substituents it leaves
        // tied go to the hierarchical comparison by Rules 1a to 3.
        int highI = higherPrioritySub(g, i, j);
        int highJ = higherPrioritySub(g, j, i);
        if (highI == TIED || highJ == TIED) {
          Character label = new SequenceRules(g, i, j).bondLabel();
          if (label != null) result.put(MolGraph.bondKey(i, j), label);
          continue;
        }
        if (highI == Integer.MIN_VALUE || highJ == Integer.MIN_VALUE) continue;

        // The stored relation uses the first non-double-bond neighbour at each end.
        int refI = g.doubleBondReference(i, j);
        int refJ = g.doubleBondReference(j, i);

        int flips = 0;
        if (highI != refI && refI >= 0) flips++;
        if (highJ != refJ && refJ >= 0) flips++;

        int effectiveConf = conf;
        if (flips == 1) effectiveConf = (conf == 1) ? 2 : 1;

        long bk = MolGraph.bondKey(i, j);
        result.put(bk, (effectiveConf == 1) ? 'Z' : 'E');
      }
    }

    return result;
  }

  /** The value of {@link #higherPrioritySub} when Rule 1a leaves the two substituents of an end tied. */
  private static final int TIED = Integer.MIN_VALUE + 1;

  /** Find the higher-priority substituent on one end of a double bond by Rule 1a. */
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

    // If two substituents have the same highest priority, Rule 1a leaves the end tied
    if (subCount >= 2 && bestRank == secondRank) return TIED;
    return bestAtom;
  }

  // ==========================================================================
  // Rules 4 and 5 with auxiliary descriptors
  // ==========================================================================

  /** Standard atomic weights by atomic number (IUPAC 2021), as the native periodic table holds them. */
  private static final double[] STANDARD_ATOMIC_WEIGHT = {
      0.0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007,
      15.999, 18.998, 20.180, 22.990, 24.305, 26.982, 28.086, 30.974,
      32.06, 35.45, 39.95, 39.098, 40.078, 44.956, 47.867, 50.942,
      51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723,
      72.630, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906,
      91.224, 92.906, 95.95, 97.0, 101.07, 102.91, 106.42, 107.87,
      112.41, 114.82, 118.71, 121.76, 127.60, 126.904, 131.29, 132.91,
      137.33, 138.91, 140.12, 140.91, 144.24, 145.0, 150.36, 151.96,
      157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97,
      178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97,
      200.59, 204.38, 207.2, 208.98, 209.0, 210.0, 222.0, 223.0,
      226.0, 227.0, 232.04, 231.04, 238.03, 237.0, 244.0, 243.0,
      247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 266.0,
      267.0, 268.0, 269.0, 270.0, 277.0, 278.0, 281.0, 282.0,
      285.0, 286.0, 289.0, 290.0, 293.0, 294.0, 294.0,
  };

  /** Standard atomic weight for Rule 2, zero outside the table. */
  static double atomicWeight(int atomicNumber) {
    return atomicNumber >= 0 && atomicNumber < STANDARD_ATOMIC_WEIGHT.length ? STANDARD_ATOMIC_WEIGHT[atomicNumber] : 0.0;
  }

  /**
   * Sequence Rules 1 to 6 (IUPAC 2013 P-92.1.3) over the hierarchical digraph of one tetrahedral
   * centre or one double bond, for the ligands Rule 1a alone leaves tied. The method is the one the
   * reference CIP implementations follow (Hanson et al., JCIM 58, 1755 (2018)):
   *
   *   - Every stereogenic unit in the digraph (a tagged tetrahedral centre or a stereogenic double
   *     bond, other than the one being labelled) receives an auxiliary descriptor at each of its nodes,
   *     found by ranking its ligands with the digraph re-rooted at that node: the edge towards the
   *     root is reversed and the subtree beyond it explored in that direction. Nodes are labelled
   *     from the deepest sphere inwards, so each sees the descriptors of the units beyond it.
   *   - The root's ligands are then ranked by Rule 1a (atomic number), Rule 1b (a ring-closure
   *     duplicate before any other node, and between two the one whose original is nearer the
   *     root), Rule 2 (the mass number of a labelled atom, against the standard atomic weight of an
   *     unlabelled one), Rule 3 (Z before E), Rule 4a (R, S, seqCis and seqTrans before r, s, E and
   *     Z before none), Rule 4b (like descriptor pairs before unlike, each branch read against its
   *     reference descriptor, the majority descriptor of the first sphere that has one), Rule 4c
   *     (r before s), Rule 5 (R before S and seqCis before seqTrans, read with both R and S as the
   *     reference) and Rule 6 (an undifferentiated reference ligand first).
   *   - Each rule is exhausted sphere by sphere before the next is consulted, and within a rule a
   *     node's branches are ordered by the rules up to and including it.
   *   - A centre whose ligands are told apart only by Rule 5 or Rule 6 is pseudoasymmetric and
   *     takes a lower-case descriptor; a double bond told apart that way at one end only is seqCis
   *     or seqTrans rather than Z or E.
   *
   * The native core applies the same rules node for node.
   */
  static final class SequenceRules {
    /** Descriptors of the units in the digraph; E and Z are constitutional, seqCis and seqTrans are not. */
    private static final int NONE = 0, UNKNOWN = 1, NS = 2, R = 3, S = 4, LOWER_R = 5, LOWER_S = 6,
        SEQ_TRANS = 7, SEQ_CIS = 8, E = 9, Z = 10;
    private static final int RULE_1A = 0, RULE_1B = 1, RULE_2 = 2, RULE_3 = 3, RULE_4A = 4, RULE_4B = 5,
        RULE_4C = 6, RULE_5 = 7, RULE_6 = 8, RULE_COUNT = 9;

    /** Recursive comparisons allowed for one centre; beyond it the centre is left unresolved. */
    private static final long COMPARISON_BUDGET = 5000000L;

    /** A rule; Rules 4b and 5 read descriptor pairs against `reference` (R or S) in their replacement form. */
    private static final class Rule {
      final int kind, reference;
      Rule(int kind, int reference) { this.kind = kind; this.reference = reference; }
    }

    /** A directed edge of the re-rooted digraph: `beg` is the node nearer the current root. */
    private static final class Edge {
      final int beg, end;
      Edge(int beg, int end) { this.beg = beg; this.end = end; }
    }

    private static final class Priority {
      final boolean unique, pseudoasymmetric;
      Priority(boolean unique, boolean pseudoasymmetric) { this.unique = unique; this.pseudoasymmetric = pseudoasymmetric; }
    }

    /** A stereogenic unit: a tagged tetrahedral centre, or a stereogenic double bond. */
    private static final class Unit {
      boolean bond;
      int first = -1, second = -1;             // the atom, or both bond atoms
      int[] ligands;                           // centre: the tag's ligand order, -1 for the implicit ligand
      int tag;                                 // centre: 1 anticlockwise, 2 clockwise
      int carrierFirst = -1, carrierSecond = -1; // bond: the reference neighbour at each end
      int configuration;                       // bond: 1 cis and 2 trans, relative to the carriers
    }

    /** Descriptor pairs of a branch, read against the reference descriptor at its head. */
    private static final class PairList {
      private final List<Integer> descriptors = new ArrayList<>();

      PairList(int reference) { add(reference); }

      void add(int descriptor) {
        int reference = pairReference(descriptor);
        if (reference != NONE) descriptors.add(reference);
      }

      /** The first pair that is like in one list and unlike in the other decides. */
      int compareTo(PairList that) {
        int count = Math.min(descriptors.size(), that.descriptors.size());
        for (int i = 1; i < count; i++) {
          boolean thisLike = descriptors.get(i).equals(descriptors.get(0));
          boolean thatLike = that.descriptors.get(i).equals(that.descriptors.get(0));
          if (thisLike && !thatLike) return 1;
          if (!thisLike && thatLike) return -1;
        }
        return 0;
      }
    }

    private final MolGraph g;
    private final int centre;
    private final int partner;           // the other atom of a double bond being labelled, -1 for a centre
    private List<DigraphNode> nodes;
    private int[] aux;                   // auxiliary descriptor by node
    private int[] up;                    // the neighbour of each node towards the current root, -1 for the root
    private final List<Unit> units = new ArrayList<>();
    private int centreUnit = -1;
    private int root = 0;                // the current root
    private int rule6Reference = -1;     // the Rule 6 reference atom, -1 for none
    private long budget = COMPARISON_BUDGET;
    private boolean exhausted = false;

    SequenceRules(MolGraph graph, int centre) {
      this.g = graph;
      this.centre = centre;
      this.partner = -1;
    }

    /** The rules for the double bond first=second, over the digraph rooted at `first`. */
    SequenceRules(MolGraph graph, int first, int second) {
      this.g = graph;
      this.centre = first;
      this.partner = second;
    }

    /**
     * The descriptor of the double bond by Rules 1a to 5, each read exhaustively before the next,
     * or null when they leave the ligands of an end tied or only Rule 5 tells the ligands of one
     * end apart (seqCis and seqTrans).
     */
    Character bondLabel() {
      if (centre < 0 || centre >= g.n || partner < 0 || partner >= g.n) return null;
      Unit unit = new Unit();
      unit.bond = true;
      unit.first = centre;
      unit.second = partner;
      unit.carrierFirst = g.doubleBondReference(centre, partner);
      unit.carrierSecond = g.doubleBondReference(partner, centre);
      unit.configuration = g.dbStereo(centre, partner);
      if (unit.configuration == 0 || unit.carrierFirst < 0 || unit.carrierSecond < 0) return null;
      nodes = buildDigraph(g, centre);
      if (nodes.isEmpty()) return null;
      aux = new int[nodes.size()];
      up = new int[nodes.size()];
      collectUnits();
      labelAuxiliaryDescriptors();
      changeRoot(0);
      List<Rule> rules = new ArrayList<>();
      for (int kind = RULE_1A; kind <= RULE_5; kind++) rules.add(new Rule(kind, NONE));
      int descriptor = labelDoubleBond(0, unit, rules);
      if (exhausted) return null;
      return descriptor == Z ? Character.valueOf('Z') : descriptor == E ? Character.valueOf('E') : null;
    }

    /** The descriptor of the centre, or null when the rules do not resolve it. */
    Character label() {
      if (centre < 0 || centre >= g.n || g.tetraChirality == null || g.tetraChirality[centre] == 0
          || g.tetrahedralLigands(centre) == null) return null;
      nodes = buildDigraph(g, centre);
      if (nodes.isEmpty()) return null;
      aux = new int[nodes.size()];
      up = new int[nodes.size()];
      collectUnits();
      if (centreUnit < 0) return null;
      labelAuxiliaryDescriptors();
      changeRoot(0);
      int descriptor = labelTetrahedral(0, units.get(centreUnit));
      if (exhausted) return null;
      switch (descriptor) {
        case R: return 'R';
        case S: return 'S';
        case LOWER_R: return 'r';
        case LOWER_S: return 's';
        default: return null;
      }
    }

    private static int threeWay(double first, double second) { return first < second ? -1 : first > second ? 1 : 0; }

    private static int pairReference(int descriptor) {
      return descriptor == R || descriptor == SEQ_CIS ? R : descriptor == S || descriptor == SEQ_TRANS ? S : NONE;
    }

    private static int rule3Order(int descriptor) { return descriptor == E ? 1 : descriptor == Z ? 2 : 0; }

    private static int rule4aOrder(int descriptor) {
      switch (descriptor) {
        case R: case S: case SEQ_CIS: case SEQ_TRANS: return 2;
        case LOWER_R: case LOWER_S: case E: case Z: return 1;
        default: return 0;
      }
    }

    private static int rule4cOrder(int descriptor) { return descriptor == LOWER_R ? 2 : descriptor == LOWER_S ? 1 : 0; }

    /** Whether the node ends its branch: a duplicate, an implicit hydrogen, or an atom with no other neighbour. */
    private boolean terminal(int node) {
      DigraphNode n = nodes.get(node);
      return n.isDuplicate || n.atomIdx < 0 || (n.parent >= 0 ? 1 : 0) + n.children.size() == 1;
    }

    /** Make `node` the root: the path from it to the original root is read in reverse. */
    private void changeRoot(int node) {
      for (int i = 0; i < nodes.size(); i++) up[i] = nodes.get(i).parent;
      up[node] = -1;
      for (int child = node, parent = nodes.get(node).parent; parent >= 0;
           child = parent, parent = nodes.get(parent).parent)
        up[parent] = child;
      root = node;
    }

    /** The edges of a node under the current root: the one towards the root first, then the rest. */
    private List<Edge> edgesOf(int node) {
      List<Edge> edges = new ArrayList<>();
      DigraphNode n = nodes.get(node);
      if (up[node] >= 0) edges.add(new Edge(up[node], node));
      if (n.parent >= 0 && n.parent != up[node]) edges.add(new Edge(node, n.parent));
      for (int child : n.children) {
        if (child == up[node]) continue;
        // The original root's multiple bonds have no duplicates in the reference digraph.
        if (node == 0 && nodes.get(child).isDuplicate && !nodes.get(child).isRingDuplicate) continue;
        edges.add(new Edge(node, child));
      }
      return edges;
    }

    private List<Edge> nonTerminalOutEdges(int node) {
      List<Edge> out = new ArrayList<>();
      for (Edge edge : edgesOf(node))
        if (edge.beg == node && !terminal(edge.end)) out.add(edge);
      return out;
    }

    /** The rules a rule orders branches with: those up to and including it; a replacement rule, itself. */
    private static List<Rule> sorterOf(Rule rule) {
      List<Rule> rules = new ArrayList<>();
      if (rule.reference != NONE) {
        rules.add(rule);
        return rules;
      }
      for (int kind = 0; kind <= rule.kind; kind++) rules.add(new Rule(kind, NONE));
      return rules;
    }

    private static List<Rule> allRules() {
      List<Rule> rules = new ArrayList<>();
      for (int kind = 0; kind < RULE_COUNT; kind++) rules.add(new Rule(kind, NONE));
      return rules;
    }

    /** All rules with `replacement` in place of its kind, read last. */
    private static List<Rule> referenceSorter(Rule replacement) {
      List<Rule> rules = new ArrayList<>();
      for (int kind = 0; kind < RULE_COUNT; kind++)
        if (kind != replacement.kind) rules.add(new Rule(kind, NONE));
      rules.add(replacement);
      return rules;
    }

    // ------------------------------------------------------------------
    // Ordering
    // ------------------------------------------------------------------

    /** Order the edges of `node` by `rules`, highest first; `top` reads every rule exhaustively. */
    private Priority prioritise(List<Rule> rules, boolean top, int node, List<Edge> edges, boolean deep) {
      boolean unique = true;
      int pseudoasymmetric = 0;
      for (int i = 0; i < edges.size(); i++) {
        for (int j = i; j > 0; j--) {
          int cmp = compareSubstituents(rules, top, node, edges.get(j - 1), edges.get(j), deep);
          if (cmp < -1 || cmp > 1) pseudoasymmetric++;
          if (cmp < 0) {
            Collections.swap(edges, j, j - 1);
          } else {
            if (cmp == 0) unique = false;
            break;
          }
        }
      }
      return new Priority(unique, pseudoasymmetric == 1);
    }

    private int compareSubstituents(List<Rule> rules, boolean top, int node, Edge first, Edge second, boolean deep) {
      boolean firstUp = first.end == node, secondUp = second.end == node;
      if (firstUp && !secondUp) return 1;
      if (!firstUp && secondUp) return -1;
      for (Rule rule : rules) {
        int cmp = top || deep ? recursiveCompare(rule, first, second) : compare(rule, first, second);
        if (cmp != 0) return cmp;
      }
      return 0;
    }

    /** Consecutive edges of a sorted list that the rules do not tell apart. */
    private List<List<Edge>> groups(List<Rule> rules, boolean top, List<Edge> sorted) {
      List<List<Edge>> result = new ArrayList<>();
      for (int i = 0; i < sorted.size(); i++) {
        if (i == 0 || compareSubstituents(rules, top, sorted.get(i - 1).beg, sorted.get(i - 1), sorted.get(i), true) != 0)
          result.add(new ArrayList<>());
        result.get(result.size() - 1).add(sorted.get(i));
      }
      return result;
    }

    /** Compare two edges by one rule at their ends only. */
    private int compare(Rule rule, Edge first, Edge second) {
      DigraphNode a = nodes.get(first.end);
      DigraphNode b = nodes.get(second.end);
      switch (rule.kind) {
        case RULE_1A:
          return threeWay(a.scaledAtomicNum, b.scaledAtomicNum);
        case RULE_1B:
          if (a.isRingDuplicate && b.isRingDuplicate) return -threeWay(a.distance, b.distance);
          if (a.isRingDuplicate) return 1;
          if (b.isRingDuplicate) return -1;
          return 0;
        case RULE_2:
          if (a.atomicNum == 0 && b.atomicNum == 0) return 0;
          if (a.atomicNum == 0 || b.atomicNum == 0) return threeWay(a.atomicNum, b.atomicNum);
          return compareMassNumbers(a, b);
        case RULE_3:
          return threeWay(rule3Order(aux[first.end]), rule3Order(aux[second.end]));
        case RULE_4A:
          return threeWay(rule4aOrder(aux[first.end]), rule4aOrder(aux[second.end]));
        case RULE_4B:
          return compareRule4b(rule, first, second);
        case RULE_4C:
          return threeWay(rule4cOrder(aux[first.end]), rule4cOrder(aux[second.end]));
        case RULE_5:
          return compareRule5(rule, first, second);
        case RULE_6:
          if (rule6Reference < 0) return 0;
          if (rule6Reference == a.atomIdx && rule6Reference != b.atomIdx) return 2;
          if (rule6Reference != a.atomIdx && rule6Reference == b.atomIdx) return -2;
          return 0;
        default:
          return 0;
      }
    }

    /** Compare two edges by one rule sphere by sphere, ordering each sphere by the rules up to it. */
    private int recursiveCompare(Rule rule, Edge first, Edge second) {
      if (exhausted || budget-- <= 0) {
        exhausted = true;
        return 0;
      }
      int cmp = compare(rule, first, second);
      if (cmp != 0) return cmp;
      List<Rule> sorter = sorterOf(rule);
      List<Edge> firstQueue = new ArrayList<>(), secondQueue = new ArrayList<>();
      firstQueue.add(first);
      secondQueue.add(second);
      for (int pos = 0; pos < firstQueue.size() && pos < secondQueue.size(); pos++) {
        first = firstQueue.get(pos);
        second = secondQueue.get(pos);
        List<Edge> firstEdges = edgesOf(first.end), secondEdges = edgesOf(second.end);
        prioritise(sorter, false, first.end, firstEdges, false);
        prioritise(sorter, false, second.end, secondEdges, false);
        int sizeDifference = threeWay(firstEdges.size(), secondEdges.size());
        for (int i = 0; i < firstEdges.size() && i < secondEdges.size(); i++) {
          if (firstEdges.get(i).end == first.end) continue; // the edges towards the root
          cmp = compare(rule, firstEdges.get(i), secondEdges.get(i));
          if (cmp != 0) return cmp;
        }
        if (sizeDifference != 0) return sizeDifference;
        prioritise(sorter, false, first.end, firstEdges, true);
        prioritise(sorter, false, second.end, secondEdges, true);
        for (int i = 0; i < firstEdges.size() && i < secondEdges.size(); i++) {
          if (firstEdges.get(i).end == first.end) continue;
          cmp = compare(rule, firstEdges.get(i), secondEdges.get(i));
          if (cmp != 0) return cmp;
          firstQueue.add(firstEdges.get(i));
          secondQueue.add(secondEdges.get(i));
        }
      }
      return 0;
    }

    // ------------------------------------------------------------------
    // Rule 4b: like descriptor pairs before unlike
    // ------------------------------------------------------------------

    /** Away from the root, a replacement rule reads the ends against its reference descriptor. */
    private int compareAgainstReference(Rule rule, Edge first, Edge second) {
      if (rule.reference == NONE) return 0;
      int a = aux[first.end], b = aux[second.end];
      if (a != NONE && b != NONE && a != NS && b != NS) {
        boolean firstLike = pairReference(rule.reference) == pairReference(a);
        boolean secondLike = pairReference(rule.reference) == pairReference(b);
        if (firstLike && !secondLike) return 1;
        if (secondLike && !firstLike) return -1;
      }
      return 0;
    }

    private int compareRule4b(Rule rule, Edge first, Edge second) {
      if (first.beg != root || second.beg != root) return compareAgainstReference(rule, first, second);
      List<Integer> firstReferences = referenceDescriptors(first.end);
      List<Integer> secondReferences = referenceDescriptors(second.end);
      if (firstReferences.isEmpty() || secondReferences.isEmpty()) return 0;
      if (firstReferences.size() == 1)
        return comparePairs(first.end, second.end, firstReferences.get(0), secondReferences.get(0));
      List<PairList> firstLists = new ArrayList<>(), secondLists = new ArrayList<>();
      for (int reference : firstReferences) {
        PairList pairs = new PairList(reference);
        fillPairs(new Rule(RULE_4B, reference), first.end, pairs);
        firstLists.add(pairs);
      }
      for (int reference : secondReferences) {
        PairList pairs = new PairList(reference);
        fillPairs(new Rule(RULE_4B, reference), second.end, pairs);
        secondLists.add(pairs);
      }
      firstLists.sort((a, b) -> b.compareTo(a));
      secondLists.sort((a, b) -> b.compareTo(a));
      for (int i = 0; i < firstLists.size() && i < secondLists.size(); i++) {
        int cmp = firstLists.get(i).compareTo(secondLists.get(i));
        if (cmp != 0) return cmp;
      }
      return 0;
    }

    /** The reference descriptor(s) of a branch: the majority of the first sphere with any, both on a tie. */
    private List<Integer> referenceDescriptors(int node) {
      List<Integer> result = new ArrayList<>();
      List<List<Integer>> level = new ArrayList<>();
      level.add(new ArrayList<>(Collections.singletonList(node)));
      while (!level.isEmpty()) {
        for (List<Integer> group : level)
          if (referenceOf(group, result)) return result;
        level = nextLevel(level);
      }
      return new ArrayList<>();
    }

    private boolean referenceOf(List<Integer> group, List<Integer> result) {
      int right = 0, left = 0;
      for (int node : group) {
        int descriptor = aux[node];
        if (descriptor == R || descriptor == SEQ_CIS) right++;
        else if (descriptor == S || descriptor == SEQ_TRANS) left++;
      }
      if (right + left == 0) return false;
      if (right > left) {
        result.add(R);
      } else if (right < left) {
        result.add(S);
      } else {
        result.add(R);
        result.add(S);
      }
      return true;
    }

    /** The next sphere, grouped by the rules up to Rule 4b so that equivalent nodes are read together. */
    private List<List<Integer>> nextLevel(List<List<Integer>> previous) {
      List<Rule> sorter = sorterOf(new Rule(RULE_4B, NONE));
      List<List<Integer>> next = new ArrayList<>();
      for (List<Integer> group : previous) {
        List<List<List<Edge>>> grouped = new ArrayList<>();
        for (int node : group) {
          List<Edge> edges = nonTerminalOutEdges(node);
          prioritise(sorter, false, node, edges, true);
          grouped.add(groups(sorter, false, edges));
        }
        int count = grouped.isEmpty() ? 0 : grouped.get(0).size();
        for (List<List<Edge>> g : grouped) count = Math.min(count, g.size());
        for (int i = 0; i < count; i++) {
          List<Integer> equivalent = new ArrayList<>();
          for (List<List<Edge>> g : grouped)
            for (Edge edge : g.get(i)) equivalent.add(edge.end);
          if (!equivalent.isEmpty()) next.add(equivalent);
        }
      }
      return next;
    }

    /** Read the descriptors of a branch outwards, each sphere ordered with the replacement rule. */
    private void fillPairs(Rule replacement, int beg, PairList pairs) {
      List<Rule> sorter = referenceSorter(replacement);
      List<Integer> queue = new ArrayList<>();
      queue.add(beg);
      for (int i = 0; i < queue.size(); i++) {
        int node = queue.get(i);
        pairs.add(aux[node]);
        List<Edge> edges = edgesOf(node);
        prioritise(sorter, false, node, edges, true);
        for (Edge edge : edges)
          if (edge.beg == node && !terminal(edge.end)) queue.add(edge.end);
      }
    }

    /** Compare two branches pair by pair, each read against its own reference descriptor. */
    private int comparePairs(int firstNode, int secondNode, int firstReference, int secondReference) {
      List<Rule> firstSorter = referenceSorter(new Rule(RULE_4B, firstReference));
      List<Rule> secondSorter = referenceSorter(new Rule(RULE_4B, secondReference));
      List<Integer> firstQueue = new ArrayList<>(), secondQueue = new ArrayList<>();
      firstQueue.add(firstNode);
      secondQueue.add(secondNode);
      for (int pos = 0; pos < firstQueue.size() && pos < secondQueue.size(); pos++) {
        int a = firstQueue.get(pos), b = secondQueue.get(pos);
        boolean firstLike = pairReference(aux[a]) == firstReference;
        boolean secondLike = pairReference(aux[b]) == secondReference;
        if (firstLike && !secondLike) return 1;
        if (!firstLike && secondLike) return -1;
        List<Edge> edges = edgesOf(a);
        prioritise(firstSorter, false, a, edges, true);
        for (Edge edge : edges)
          if (edge.beg == a && !terminal(edge.end)) firstQueue.add(edge.end);
        edges = edgesOf(b);
        prioritise(secondSorter, false, b, edges, true);
        for (Edge edge : edges)
          if (edge.beg == b && !terminal(edge.end)) secondQueue.add(edge.end);
      }
      return 0;
    }

    // ------------------------------------------------------------------
    // Rule 5: R before S, read with both references
    // ------------------------------------------------------------------

    private int compareRule5(Rule rule, Edge first, Edge second) {
      if (first.beg != root || second.beg != root) return compareAgainstReference(rule, first, second);
      PairList firstR = new PairList(R), secondR = new PairList(R), firstS = new PairList(S), secondS = new PairList(S);
      fillPairs(new Rule(RULE_5, R), first.end, firstR);
      fillPairs(new Rule(RULE_5, S), first.end, firstS);
      fillPairs(new Rule(RULE_5, R), second.end, secondR);
      fillPairs(new Rule(RULE_5, S), second.end, secondS);
      int cmpR = firstR.compareTo(secondR);
      int cmpS = firstS.compareTo(secondS);
      // Two when the references disagree, the mark of a pseudoasymmetric distinction.
      if (cmpR < 0) return cmpS < 0 ? -1 : -2;
      if (cmpR > 0) return cmpS > 0 ? 1 : 2;
      return 0;
    }

    // ------------------------------------------------------------------
    // Labelling
    // ------------------------------------------------------------------

    private void collectUnits() {
      for (int atom = 0; atom < g.n; atom++) {
        if (g.tetraChirality[atom] == 0) continue;
        int[] ligands = g.tetrahedralLigands(atom);
        if (ligands == null) continue;
        Unit unit = new Unit();
        unit.first = atom;
        unit.ligands = ligands;
        unit.tag = g.tetraChirality[atom];
        if (atom == centre) centreUnit = units.size();
        units.add(unit);
      }
      for (int atom = 0; atom < g.n; atom++) {
        for (int other : g.neighbors[atom]) {
          if (other <= atom || g.bondOrder(atom, other) != 2) continue;
          int configuration = g.dbStereo(atom, other);
          if (configuration == 0) continue;
          Unit unit = new Unit();
          unit.bond = true;
          unit.first = atom;
          unit.second = other;
          unit.carrierFirst = g.doubleBondReference(atom, other);
          unit.carrierSecond = g.doubleBondReference(other, atom);
          unit.configuration = configuration;
          if (unit.carrierFirst >= 0 && unit.carrierSecond >= 0) units.add(unit);
        }
      }
    }

    /** Label every node of every other unit, deepest first; a sphere sees only the spheres beyond it. */
    private void labelAuxiliaryDescriptors() {
      List<int[]> pending = new ArrayList<>(); // node, unit
      for (int u = 0; u < units.size(); u++) {
        Unit unit = units.get(u);
        if (!unit.bond && unit.first == centre) continue;
        if (unit.bond && partner >= 0 && ((unit.first == centre && unit.second == partner)
                                          || (unit.first == partner && unit.second == centre))) continue;
        for (int node = 0; node < nodes.size(); node++) {
          if (nodes.get(node).isDuplicate || nodes.get(node).atomIdx != unit.first) continue;
          int low = node;
          // A bond is labelled at its end nearer the root.
          int parent = nodes.get(node).parent;
          if (unit.bond && parent >= 0 && nodes.get(parent).atomIdx == unit.second) low = parent;
          pending.add(new int[] {low, u});
        }
      }
      pending.sort((a, b) -> Integer.compare(nodes.get(b[0]).depth, nodes.get(a[0]).depth));
      List<int[]> sphere = new ArrayList<>(); // node, descriptor
      int previous = Integer.MAX_VALUE;
      for (int[] entry : pending) {
        int depth = nodes.get(entry[0]).depth;
        if (depth < previous) {
          for (int[] labelled : sphere) aux[labelled[0]] = labelled[1];
          sphere.clear();
          previous = depth;
        }
        Unit unit = units.get(entry[1]);
        int descriptor;
        if (unit.bond) {
          descriptor = labelDoubleBond(entry[0], unit, allRules());
        } else {
          changeRoot(entry[0]);
          descriptor = labelTetrahedral(entry[0], unit);
        }
        boolean replaced = false;
        for (int[] labelled : sphere) {
          if (labelled[0] == entry[0]) {
            labelled[1] = descriptor;
            replaced = true;
          }
        }
        if (!replaced) sphere.add(new int[] {entry[0], descriptor});
      }
      for (int[] labelled : sphere) aux[labelled[0]] = labelled[1];
    }

    /** Whether reordering `first` into `second` is an odd permutation of the same four edges. */
    private static boolean oddPermutation(List<Edge> first, List<Edge> second) {
      List<Integer> positions = new ArrayList<>();
      for (Edge edge : second)
        for (int i = 0; i < first.size(); i++)
          if (first.get(i).beg == edge.beg && first.get(i).end == edge.end) positions.add(i);
      int inversions = 0;
      for (int i = 0; i < positions.size(); i++)
        for (int j = i + 1; j < positions.size(); j++)
          if (positions.get(i) > positions.get(j)) inversions++;
      return inversions % 2 == 1;
    }

    /** The descriptor of a tetrahedral unit at `node`, the current root. */
    private int labelTetrahedral(int node, Unit unit) {
      List<Edge> edges = edgesOf(node);
      if (edges.size() < 3) return NS;
      List<Rule> rules = allRules();
      Priority priority = prioritise(rules, true, node, edges, true);
      if (!priority.unique && edges.size() == 4) {
        List<List<Edge>> partition = groups(rules, true, edges);
        if (partition.size() == 2) {
          rule6Reference = nodes.get(edges.get(1).end).atomIdx;
          priority = prioritise(rules, true, node, edges, true);
          rule6Reference = -1;
        } else if (partition.size() == 1) {
          rule6Reference = nodes.get(edges.get(0).end).atomIdx;
          prioritise(rules, true, node, edges, true);
          List<Edge> first = new ArrayList<>(edges);
          rule6Reference = nodes.get(edges.get(1).end).atomIdx;
          priority = prioritise(rules, true, node, edges, true);
          rule6Reference = -1;
          if (oddPermutation(first, edges)) return UNKNOWN;
        }
        if (!priority.unique) return UNKNOWN;
      } else if (!priority.unique) {
        return UNKNOWN;
      }
      // The carriers in rank order: every ligand atom once, and the implicit hydrogen.
      List<Integer> carriers = new ArrayList<>();
      for (Edge edge : edges) {
        DigraphNode end = nodes.get(edge.end);
        if (end.isDuplicate && !end.isRingDuplicate) continue;
        carriers.add(end.atomIdx);
      }
      int[] ranks = new int[4];
      for (int slot = 0; slot < 4; slot++) {
        int ligand = unit.ligands[slot];
        int rank = 0; // a lone pair ranks below every carrier
        for (int i = 0; i < carriers.size(); i++) {
          if (carriers.get(i) == ligand) {
            rank = carriers.size() - i;
            break;
          }
        }
        if (ligand >= 0 && rank == 0) return UNKNOWN;
        ranks[slot] = rank;
      }
      int inversions = 0;
      for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
          if (ranks[i] > ranks[j]) inversions++;
      boolean anticlockwise = unit.tag == 1;
      boolean rectus = (inversions % 2 == 0) != anticlockwise;
      if (priority.pseudoasymmetric) return rectus ? LOWER_R : LOWER_S;
      return rectus ? R : S;
    }

    /** The descriptor of a double bond unit whose end nearer the root is `root1`, ranking its ligands by `rules`. */
    private int labelDoubleBond(int root1, Unit unit, List<Rule> rules) {
      int atom1 = nodes.get(root1).atomIdx;
      int atom2 = atom1 == unit.first ? unit.second : unit.first;
      int root2 = -1;
      for (int child : nodes.get(root1).children)
        if (!nodes.get(child).isDuplicate && nodes.get(child).atomIdx == atom2) root2 = child;
      if (root2 < 0) return UNKNOWN;
      changeRoot(root1);
      List<Edge> edges1 = withoutInternal(root1, atom2);
      changeRoot(root2);
      List<Edge> edges2 = withoutInternal(root2, atom1);
      if (edges1.isEmpty() || edges2.isEmpty()) return NS;
      int carrier1 = unit.carrierFirst, carrier2 = unit.carrierSecond;
      if (atom1 == unit.second) {
        carrier1 = unit.carrierSecond;
        carrier2 = unit.carrierFirst;
      }
      int configuration = unit.configuration;
      changeRoot(root1);
      Priority priority1 = prioritise(rules, true, root1, edges1, true);
      if (!priority1.unique) return UNKNOWN;
      if (edges1.size() > 1 && nodes.get(edges1.get(0).end).atomIdx != carrier1) configuration = 3 - configuration;
      changeRoot(root2);
      Priority priority2 = prioritise(rules, true, root2, edges2, true);
      if (!priority2.unique) return UNKNOWN;
      if (edges2.size() > 1 && nodes.get(edges2.get(0).end).atomIdx != carrier2) configuration = 3 - configuration;
      boolean pseudoasymmetric = priority1.pseudoasymmetric != priority2.pseudoasymmetric;
      if (configuration == 1) return pseudoasymmetric ? SEQ_CIS : Z;
      return pseudoasymmetric ? SEQ_TRANS : E;
    }

    /** The edges of a bond end other than those to the bond's other atom. */
    private List<Edge> withoutInternal(int node, int otherAtom) {
      List<Edge> edges = new ArrayList<>();
      for (Edge edge : edgesOf(node))
        if (nodes.get(edge.end).atomIdx != otherAtom && nodes.get(edge.beg).atomIdx != otherAtom) edges.add(edge);
      return edges;
    }
  }
}

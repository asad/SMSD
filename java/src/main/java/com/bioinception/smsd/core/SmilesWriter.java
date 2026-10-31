/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Writes the canonical SMILES of a graph: components and branches in the order of the canonical
 * labelling, ring closures reusing labels as needed, and tetrahedral and double-bond stereo written
 * relative to that order.
 */
final class SmilesWriter {
  /** The element symbols by atomic number, with "?" at 0: the one element table. */
  static final String[] ELEMENT_SYMBOL = {
      "?","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S",
      "Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
      "Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
      "Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm",
      "Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os",
      "Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa",
      "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg",
      "Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
  };

  /** Whether the element is in the SMILES organic subset, which is written without brackets. */
  private static boolean organicSubset(int atomicNumber) {
    switch (atomicNumber) {
      case 5: case 6: case 7: case 8: case 9: case 15: case 16: case 17: case 35: case 53:
        return true;
      default:
        return false;
    }
  }

  private final MolGraph graph;
  private final int n;
  private final int[][] neighbors;
  private final int[] canonicalLabel, atomicNum, formalCharge, massNumber, tetraChirality;
  private final boolean[] aromatic;

  private SmilesWriter(MolGraph graph) {
    graph.ensureCanonical();
    this.graph = graph;
    this.n = graph.n;
    this.neighbors = graph.neighbors;
    this.canonicalLabel = graph.canonicalLabel;
    this.atomicNum = graph.atomicNum;
    this.formalCharge = graph.formalCharge;
    this.massNumber = graph.massNumber;
    this.tetraChirality = graph.tetraChirality;
    this.aromatic = graph.aromatic;
  }

  /**
   * The canonical SMILES of {@code graph}, with positive map labels keyed by atom index. Labels
   * change the text only, not the canonical order.
   */
  static String write(MolGraph graph, Map<Integer, Integer> atomMaps) {
    Objects.requireNonNull(atomMaps, "atomMaps");
    for (Map.Entry<Integer, Integer> entry : atomMaps.entrySet())
      if (entry.getKey() == null || entry.getKey() < 0 || entry.getKey() >= graph.n
          || entry.getValue() == null || entry.getValue() <= 0)
        throw new IllegalArgumentException("Atom map labels require valid atom indices and positive values");
    if (graph.n == 0) return "";
    return new SmilesWriter(graph).write(atomMaps);
  }

  private String write(Map<Integer, Integer> atomMaps) {
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[canonicalLabel[i]] = i;
    int[] color = new int[n];
    @SuppressWarnings("unchecked")
    List<int[]>[] ringClosures = new List[n];
    Set<Long> ringEdges = new HashSet<>();
    int[] nextRing = {1};
    for (int ci = 0; ci < n; ci++) {
      int start = inv[ci];
      if (color[start] == 0) preScanDfs(start, -1, color, ringClosures, ringEdges, nextRing);
    }
    Map<Long, Integer> bondDirections = smilesBondDirections();
    boolean[] visited = new boolean[n];
    int[] ringLabels = new int[nextRing[0]];
    boolean[] usedRingLabels = new boolean[100];
    StringBuilder sb = new StringBuilder();
    for (int ci = 0; ci < n; ci++) {
      int start = inv[ci];
      if (visited[start]) continue;
      if (sb.length() > 0) sb.append('.');
      emitSmiles(start, -1, visited, sb, ringClosures, ringEdges, bondDirections, atomMaps,
          ringLabels, usedRingLabels);
    }
    return sb.toString();
  }

  private static final class ScanFrame {
    final int atom, parent;
    final int[] neighbors;
    int nextNeighbor;
    boolean skippedParent;

    ScanFrame(int atom, int parent, int[] neighbors) {
      this.atom = atom;
      this.parent = parent;
      this.neighbors = neighbors;
    }
  }

  private void preScanDfs(int atom, int parent, int[] color, List<int[]>[] rc,
                         Set<Long> ringEdges, int[] nextRing) {
    ArrayDeque<ScanFrame> pending = new ArrayDeque<>();
    color[atom] = 1;
    pending.push(new ScanFrame(atom, parent, sortedNeighborsByCanonical(atom)));
    while (!pending.isEmpty()) {
      ScanFrame frame = pending.peek();
      if (frame.nextNeighbor == frame.neighbors.length) {
        color[frame.atom] = 2;
        pending.pop();
        continue;
      }
      int neighbor = frame.neighbors[frame.nextNeighbor++];
      if (neighbor == frame.parent && !frame.skippedParent) {
        frame.skippedParent = true;
        continue;
      }
      if (color[neighbor] == 0) {
        color[neighbor] = 1;
        pending.push(new ScanFrame(neighbor, frame.atom, sortedNeighborsByCanonical(neighbor)));
      } else if (color[neighbor] == 1) {
        int ring = nextRing[0]++;
        if (rc[neighbor] == null) rc[neighbor] = new ArrayList<>();
        if (rc[frame.atom] == null) rc[frame.atom] = new ArrayList<>();
        rc[neighbor].add(new int[] {ring, frame.atom});
        rc[frame.atom].add(new int[] {ring, neighbor});
        ringEdges.add(MolGraph.bondKey(frame.atom, neighbor));
      }
    }
  }

  private static final class EmissionFrame {
    final int atom, childCount;
    final int[] neighbors;
    final boolean closeBranch;
    int nextNeighbor, childIndex;

    EmissionFrame(int atom, int[] neighbors, int childCount, boolean closeBranch) {
      this.atom = atom;
      this.neighbors = neighbors;
      this.childCount = childCount;
      this.closeBranch = closeBranch;
    }
  }

  private void emitSmiles(int atom, int parent, boolean[] visited, StringBuilder sb,
                          List<int[]>[] rc, Set<Long> ringEdges, Map<Long, Integer> bondDirections,
                          Map<Integer, Integer> atomMaps, int[] ringLabels, boolean[] usedRingLabels) {
    ArrayDeque<EmissionFrame> pending = new ArrayDeque<>();
    pending.push(emitAtom(atom, parent, false, visited, sb, rc, ringEdges, bondDirections,
        atomMaps, ringLabels, usedRingLabels));
    while (!pending.isEmpty()) {
      EmissionFrame frame = pending.peek();
      if (frame.nextNeighbor == frame.neighbors.length) {
        pending.pop();
        if (frame.closeBranch) sb.append(')');
        continue;
      }
      int neighbor = frame.neighbors[frame.nextNeighbor++];
      if (visited[neighbor] || ringEdges.contains(MolGraph.bondKey(frame.atom, neighbor))) continue;
      frame.childIndex++;
      boolean needParen = frame.childCount > 1 && frame.childIndex < frame.childCount;
      if (needParen) sb.append('(');
      writeBondSymbol(frame.atom, neighbor, sb, bondDirections);
      pending.push(emitAtom(neighbor, frame.atom, needParen, visited, sb, rc, ringEdges,
          bondDirections, atomMaps, ringLabels, usedRingLabels));
    }
  }

  private EmissionFrame emitAtom(int atom, int parent, boolean closeBranch, boolean[] visited,
                                StringBuilder sb, List<int[]>[] rc, Set<Long> ringEdges,
                                Map<Long, Integer> bondDirections, Map<Integer, Integer> atomMaps,
                                int[] ringLabels, boolean[] usedRingLabels) {
    visited[atom] = true;
    int[] sorted = sortedNeighborsByCanonical(atom);
    writeAtom(atom, parent, sorted, rc, ringEdges, sb, atomMaps);
    if (rc[atom] != null) {
      for (int[] entry : rc[atom]) {
        if (!visited[entry[1]]) writeBondSymbol(atom, entry[1], sb, bondDirections);
        int label = ringLabels[entry[0]];
        if (label == 0) {
          label = entry[0];
          if (label >= usedRingLabels.length || usedRingLabels[label]) {
            label = 1;
            while (label < usedRingLabels.length && usedRingLabels[label]) label++;
          }
          if (label == usedRingLabels.length)
            throw new IllegalArgumentException("SMILES output needs more than 99 simultaneously open ring closures");
          usedRingLabels[label] = true;
          ringLabels[entry[0]] = label;
        } else {
          usedRingLabels[label] = false;
          ringLabels[entry[0]] = 0;
        }
        writeRingDigit(label, sb);
      }
    }
    int childCount = 0;
    for (int neighbor : sorted)
      if (!visited[neighbor] && !ringEdges.contains(MolGraph.bondKey(atom, neighbor))) childCount++;
    return new EmissionFrame(atom, sorted, childCount, closeBranch);
  }

  private int[] sortedNeighborsByCanonical(int atom) {
    int[] sorted = neighbors[atom].clone();
    for (int a = 0; a < sorted.length - 1; a++)
      for (int b2 = a + 1; b2 < sorted.length; b2++)
        if (canonicalLabel[sorted[a]] > canonicalLabel[sorted[b2]]) {
          int tmp = sorted[a]; sorted[a] = sorted[b2]; sorted[b2] = tmp;
        }
    return sorted;
  }

  /**
   * Write an atom token, including tetrahedral stereo (@/@@) when present.
   * The token is always bracket-form when stereo is present.
   *
   * Compare the emitted ligand order with the stored winding reference.
   */
  private void writeAtom(int atom, int parent, int[] sortedChildren,
                         List<int[]>[] rc, Set<Long> ringEdges, StringBuilder sb,
                         Map<Integer, Integer> atomMaps) {
    int z = atomicNum[atom], charge = formalCharge[atom];
    boolean arom = aromatic[atom];
    int chiral = (tetraChirality != null) ? tetraChirality[atom] : 0;
    // The tag is 1 for anticlockwise (@) and 2 for clockwise (@@); below 1 is clockwise.
    if (chiral != 0) chiral = 3 - chiral;
    int mass = (massNumber != null) ? massNumber[atom] : 0;
    // Only emit @/@@ when stereo is defined.
    boolean hasStereo = chiral != 0;
    boolean hasIsotope = mass != 0;
    int hCount = graph.hydrogenCount(atom);
    Integer atomMap = atomMaps.get(atom);
    boolean aromaticStereoEnd = false;
    if (arom) for (int neighbor : neighbors[atom])
      if (graph.bondOrder(atom, neighbor) == 4 && graph.dbStereo(atom, neighbor) != 0) aromaticStereoEnd = true;
    // The bare organic-subset form is legal only without stereo, charge, isotope label and map
    // number; otherwise the atom is written in brackets.
    if (organicSubset(z) && charge == 0 && !hasStereo && !hasIsotope && !aromaticStereoEnd && atomMap == null
        && hCount == graph.defaultHydrogenCount(atom) && !(arom && hCount > 0 && (z == 7 || z == 15))) {
      String sym = ELEMENT_SYMBOL[z];
      sb.append(arom ? Character.toLowerCase(sym.charAt(0)) : sym);
    } else {
      sb.append('[');
      if (hasIsotope) sb.append(mass);
      String elemSym = (z >= 0 && z < ELEMENT_SYMBOL.length && ELEMENT_SYMBOL[z] != null) ? ELEMENT_SYMBOL[z] : "?";
      if (arom) {
        sb.append(Character.toLowerCase(elemSym.charAt(0)));
        if (elemSym.length() > 1) sb.append(elemSym.substring(1));
      } else { sb.append(elemSym); }
      if (hasStereo) {
        int[] storedLigands = graph.tetrahedralLigands(atom);
        List<Integer> traversalSeq = new ArrayList<>();
        if (parent >= 0) traversalSeq.add(parent);
        if (storedLigands != null)
          for (int ligand : storedLigands) if (ligand == -1) traversalSeq.add(-1);
        if (rc[atom] != null) {
          for (int[] entry : rc[atom]) traversalSeq.add(entry[1]);
        }
        for (int nb : sortedChildren) {
          if (!traversalSeq.contains(nb)) traversalSeq.add(nb);
        }
        boolean oddPermutation = false;
        if (storedLigands != null && traversalSeq.size() == storedLigands.length) {
          for (int first = 0; first < storedLigands.length; first++)
            for (int second = first + 1; second < storedLigands.length; second++)
              if (traversalSeq.indexOf(storedLigands[first]) > traversalSeq.indexOf(storedLigands[second]))
                oddPermutation = !oddPermutation;
        }
        sb.append((chiral == 1) ^ oddPermutation ? "@@" : "@");
      }
      // Emit implicit hydrogen count inside brackets. This is required for
      // stereo centres (e.g. [C@@H]) and is standard practice for all bracket
      // atoms — the reader cannot infer H count once brackets are in play.
      if (hCount > 0) {
        sb.append('H');
        if (hCount > 1) sb.append(hCount);
      }
      if (charge > 0) { sb.append('+'); if (charge > 1) sb.append(charge); }
      else if (charge < 0) { sb.append('-'); if (charge < -1) sb.append(-charge); }
      if (atomMap != null) sb.append(':').append(atomMap);
      sb.append(']');
    }
  }

  /** Assign consistent carrier directions across conjugated stereo elements. */
  private Map<Long, Integer> smilesBondDirections() {
    Map<Long, List<long[]>> constraints = new HashMap<>();
    for (int firstAtom = 0; firstAtom < n; firstAtom++) {
      for (int secondAtom : neighbors[firstAtom]) {
        if (firstAtom >= secondAtom || graph.dbStereo(firstAtom, secondAtom) == 0) continue;
        int firstReference = smilesStereoReference(firstAtom, secondAtom);
        int secondReference = smilesStereoReference(secondAtom, firstAtom);
        if (firstReference < 0 || secondReference < 0) continue;
        int conformation = graph.dbStereo(firstAtom, secondAtom);
        if (firstReference != graph.doubleBondReference(firstAtom, secondAtom)) conformation = 3 - conformation;
        if (secondReference != graph.doubleBondReference(secondAtom, firstAtom)) conformation = 3 - conformation;
        addSmilesDirectionConstraint(constraints, firstAtom, firstReference,
            secondAtom, secondReference, conformation);
      }
    }
    for (int atom = 0; atom < n; atom++) {
      boolean doubleBondEnd = false;
      for (int neighbor : neighbors[atom])
        if (graph.bondOrder(atom, neighbor) == 2 || (graph.bondOrder(atom, neighbor) == 4 && graph.dbStereo(atom, neighbor) != 0))
          doubleBondEnd = true;
      if (!doubleBondEnd) continue;
      int firstReference = -1;
      for (int neighbor : neighbors[atom]) {
        if (!constraints.containsKey(MolGraph.bondKey(atom, neighbor))) continue;
        if (firstReference < 0) firstReference = neighbor;
        else addSmilesDirectionConstraint(constraints, atom, firstReference, atom, neighbor, 2);
      }
    }
    List<Long> bonds = new ArrayList<>(constraints.keySet());
    bonds.sort(Comparator.comparingInt((Long key) -> Math.min(canonicalLabel[(int) (key >>> 32)], canonicalLabel[(int) (long) key]))
        .thenComparingInt(key -> Math.max(canonicalLabel[(int) (key >>> 32)], canonicalLabel[(int) (long) key])));
    Map<Long, Integer> directions = new HashMap<>();
    ArrayDeque<Long> pending = new ArrayDeque<>();
    for (long bond : bonds) {
      if (directions.containsKey(bond)) continue;
      directions.put(bond, canonicalLabel[(int) (bond >>> 32)] < canonicalLabel[(int) bond] ? 1 : -1);
      pending.add(bond);
      while (!pending.isEmpty()) {
        long current = pending.removeFirst();
        for (long[] constraint : constraints.get(current)) {
          int direction = directions.get(current) * (int) constraint[1];
          Integer previous = directions.putIfAbsent(constraint[0], direction);
          if (previous == null) pending.addLast(constraint[0]);
          else if (previous != direction) throw new IllegalArgumentException("Inconsistent double-bond stereochemistry");
        }
      }
    }
    return directions;
  }

  /**
   * The carrier of a configured double bond at one end. '/' and '\' are single bonds in
   * OpenSMILES, so it is the single, non-aromatic bond to the neighbour of lowest canonical label;
   * only an end with no such bond, such as an aromatic ring atom carrying an exocyclic double bond,
   * takes its lowest aromatic ring bond instead, as RDKit writes it, so that the configuration is
   * kept. -1 when the end has neither.
   */
  private int smilesStereoReference(int atom, int otherEnd) {
    int single = -1, any = -1;
    for (int neighbor : neighbors[atom]) {
      int order = graph.bondOrder(atom, neighbor);
      if (neighbor == otherEnd || (order != 1 && order != 4)) continue;
      if (any < 0 || canonicalLabel[neighbor] < canonicalLabel[any]) any = neighbor;
      if (order == 1 && !graph.bondAromatic(atom, neighbor)
          && (single < 0 || canonicalLabel[neighbor] < canonicalLabel[single]))
        single = neighbor;
    }
    return single >= 0 ? single : any;
  }

  /** Relate directions away from the two ends of one stereo double bond. */
  private static void addSmilesDirectionConstraint(Map<Long, List<long[]>> constraints,
      int firstAtom, int firstReference, int secondAtom, int secondReference, int conformation) {
    long firstBond = MolGraph.bondKey(firstAtom, firstReference), secondBond = MolGraph.bondKey(secondAtom, secondReference);
    int relation = (conformation == 1 ? 1 : -1)
        * (firstAtom < firstReference ? 1 : -1) * (secondAtom < secondReference ? 1 : -1);
    constraints.computeIfAbsent(firstBond, key -> new ArrayList<>()).add(new long[] {secondBond, relation});
    constraints.computeIfAbsent(secondBond, key -> new ArrayList<>()).add(new long[] {firstBond, relation});
  }

  /**
   * Write the bond token between from→to.
   * For single or aromatic carriers of a stereo double bond we emit '/' or '\'.
   * Direction is relative to the from→to traversal direction.
   */
  private void writeBondSymbol(int from, int to, StringBuilder sb, Map<Long, Integer> bondDirections) {
    int order = graph.bondOrder(from, to);
    boolean fromArom = aromatic[from], toArom = aromatic[to];
    if (order == 4 && graph.dbStereo(from, to) != 0) { sb.append('='); return; }
    if (order == 1 || order == 4) {
      Integer direction = bondDirections.get(MolGraph.bondKey(from, to));
      if (direction != null) {
        sb.append(direction * (from < to ? 1 : -1) > 0 ? '/' : '\\');
        return;
      }
    }
    if (order == 4) return;
    if (order == 2) { sb.append('='); return; }
    if (order == 3) { sb.append('#'); return; }
    if (order == 1) {
      // Single bond between TWO aromatic atoms must be emitted explicitly
      // as '-' so that biphenyl c1ccc(-c2ccccc2)cc1 is distinguished from an
      // aromatic ring fusion. Otherwise single bonds are implicit.
      if (fromArom && toArom) { sb.append('-'); return; }
    }
  }

  private void writeRingDigit(int rid, StringBuilder sb) {
    if (rid < 10) sb.append(rid);
    else { sb.append('%'); sb.append(rid); }
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Writes portable SMARTS from source graph properties and selected edges. Atom and branch order
 * follow the source indices; the text is deterministic for that numbering, rather than canonical.
 * Stereo, tautomer-aware, complete-ring and strict ring-fusion predicates cannot be exported.
 * Strict aromaticity requires source and target flags consistent with Java SMARTS preparation.
 */
public final class SmartsWriter {
  private final MolGraph graph;
  private final ChemOptions chemistry;
  private final int[] atoms;
  private final int[][] neighbors;
  private final List<Integer>[] children;
  private final List<Closure>[] closures;

  /** An edge outside the traversal tree, with a ring number only while the closure is open. */
  private static final class Closure {
    final int first, second;
    int number;
    Closure(int first, int second) {
      this.first = first;
      this.second = second;
    }
  }

  /** A portable query for every atom and bond, using the default chemistry. */
  public static String toSmarts(MolGraph graph) {
    return toSmarts(graph, new ChemOptions());
  }

  /**
   * A portable query for every atom and bond with the representable chemistry constraints.
   * Formal charge includes neutral atoms; an unlabelled isotope remains a wildcard. Ring and
   * aromaticity flags come from the source graph. Unsupported predicates throw
   * {@link IllegalArgumentException}, as do strict aromatic flags changed by SMARTS preparation.
   * Targets must use that same aromaticity context; SMARTS does not encode the induced MCIS objective.
   */
  public static String toSmarts(MolGraph graph, ChemOptions chemistry) {
    Objects.requireNonNull(graph, "graph");
    Set<Integer> atoms = new TreeSet<>();
    Set<SearchResult.Bond> bonds = new HashSet<>();
    for (int atom = 0; atom < graph.n; atom++) {
      atoms.add(atom);
      for (int neighbor : graph.neighbors[atom])
        if (atom < neighbor) bonds.add(new SearchResult.Bond(atom, neighbor));
    }
    return write(graph, atoms, bonds, chemistry);
  }

  /** The selected atoms and bonds, read in their original graph context. */
  static String write(MolGraph graph, Collection<Integer> atoms,
                      Set<SearchResult.Bond> bonds, ChemOptions chemistry) {
    return new SmartsWriter(graph, atoms, bonds, chemistry).write();
  }

  @SuppressWarnings("unchecked")
  private SmartsWriter(MolGraph graph, Collection<Integer> selectedAtoms,
                       Set<SearchResult.Bond> selectedBonds, ChemOptions chemistry) {
    this.graph = Objects.requireNonNull(graph, "graph");
    this.chemistry = ChemOptions.copyOf(Objects.requireNonNull(chemistry, "chemistry"));
    if (this.chemistry.useChirality || this.chemistry.useBondStereo || this.chemistry.tautomerAware
        || this.chemistry.completeRingsOnly
        || this.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT)
      throw new IllegalArgumentException("SMARTS export cannot represent stereo, tautomer-aware, complete-ring or strict ring-fusion matching");
    Objects.requireNonNull(this.chemistry.matchBondOrder, "matchBondOrder");
    Objects.requireNonNull(this.chemistry.aromaticityMode, "aromaticityMode");
    Objects.requireNonNull(this.chemistry.ringFusionMode, "ringFusionMode");
    atoms = selectedAtoms.stream().mapToInt(Integer::intValue).sorted().distinct().toArray();
    int[] index = new int[graph.n];
    Arrays.fill(index, -1);
    for (int atom = 0; atom < atoms.length; atom++) {
      if (atoms[atom] < 0 || atoms[atom] >= graph.n)
        throw new IllegalArgumentException("Selected atom index out of range");
      index[atoms[atom]] = atom;
    }
    List<Integer>[] adjacency = new List[atoms.length];
    children = new List[atoms.length];
    closures = new List[atoms.length];
    for (int atom = 0; atom < atoms.length; atom++) {
      adjacency[atom] = new ArrayList<>();
      children[atom] = new ArrayList<>();
      closures[atom] = new ArrayList<>();
    }
    for (SearchResult.Bond bond : selectedBonds) {
      if (bond.second() >= graph.n || index[bond.first()] < 0 || index[bond.second()] < 0
          || !graph.hasBond(bond.first(), bond.second()))
        throw new IllegalArgumentException("Selected bond must join selected source atoms");
      int first = index[bond.first()], second = index[bond.second()];
      adjacency[first].add(second);
      adjacency[second].add(first);
      bond(first, second);
    }
    neighbors = new int[atoms.length][];
    for (int atom = 0; atom < atoms.length; atom++)
      neighbors[atom] = adjacency[atom].stream().mapToInt(Integer::intValue).sorted().toArray();
    validateAromaticContext(selectedBonds);
  }

  /** Strict flags must survive the preparation used for Java SMARTS targets, in the full graph. */
  private void validateAromaticContext(Set<SearchResult.Bond> bonds) {
    if (chemistry.aromaticityMode != ChemOptions.AromaticityMode.STRICT || atoms.length == 0) return;
    IAtomContainer prepared = MatchGraph.toMolecule(graph);
    RingPerception.markRingAtomsAndBonds(prepared);
    Aromaticity.apply(Aromaticity.Model.Daylight, prepared);
    for (int atom : atoms)
      if (graph.aromatic[atom] != prepared.getAtom(atom).isAromatic())
        throw new IllegalArgumentException("Strict aromatic flags change during SMARTS target preparation");
    for (SearchResult.Bond bond : bonds)
      if (graph.bondAromatic(bond.first(), bond.second())
          != prepared.getBond(prepared.getAtom(bond.first()), prepared.getAtom(bond.second())).isAromatic())
        throw new IllegalArgumentException("Strict aromatic flags change during SMARTS target preparation");
  }

  /** Build a traversal tree without recursing on the Java stack. */
  private void plan() {
    int[] order = new int[atoms.length], parent = new int[atoms.length], cursor = new int[atoms.length];
    int[] stack = new int[atoms.length];
    Arrays.fill(order, -1);
    Arrays.fill(parent, -1);
    int next = 0;
    for (int root = 0; root < atoms.length; root++) {
      if (order[root] >= 0) continue;
      int depth = 0;
      stack[0] = root;
      order[root] = next++;
      while (depth >= 0) {
        int atom = stack[depth];
        if (cursor[atom] == neighbors[atom].length) { depth--; continue; }
        int neighbor = neighbors[atom][cursor[atom]++];
        if (neighbor == parent[atom]) continue;
        if (order[neighbor] < 0) {
          parent[neighbor] = atom;
          children[atom].add(neighbor);
          order[neighbor] = next++;
          stack[++depth] = neighbor;
        } else if (order[neighbor] < order[atom]) {
          Closure closure = new Closure(neighbor, atom);
          closures[neighbor].add(closure);
          closures[atom].add(closure);
        }
      }
    }
  }

  private String write() {
    plan();
    StringBuilder output = new StringBuilder();
    boolean[] emitted = new boolean[atoms.length], open = new boolean[100], closeBranch = new boolean[atoms.length];
    int[] stack = new int[atoms.length], nextChild = new int[atoms.length];
    for (int root = 0; root < atoms.length; root++) {
      if (emitted[root]) continue;
      if (output.length() > 0) output.append('.');
      int depth = 0;
      stack[0] = root;
      nextChild[0] = 0;
      closeBranch[0] = false;
      atom(root, output, open);
      emitted[root] = true;
      while (depth >= 0) {
        int atom = stack[depth];
        if (nextChild[depth] == children[atom].size()) {
          if (closeBranch[depth]) output.append(')');
          depth--;
          continue;
        }
        int child = children[atom].get(nextChild[depth]++);
        boolean branch = nextChild[depth] < children[atom].size();
        if (branch) output.append('(');
        output.append(bond(atom, child));
        atom(child, output, open);
        emitted[child] = true;
        stack[++depth] = child;
        nextChild[depth] = 0;
        closeBranch[depth] = branch;
      }
    }
    return output.toString();
  }

  private void atom(int atom, StringBuilder output, boolean[] open) {
    int source = atoms[atom];
    String element = chemistry.matchAtomType ? "#" + graph.atomicNum[source] : "*";
    output.append('[');
    int isotope = graph.massNumber[source];
    if (chemistry.matchIsotope && isotope > 0)
      output.append('0').append(element).append(',').append(isotope).append(element);
    else output.append(element);
    if (chemistry.aromaticityMode == ChemOptions.AromaticityMode.STRICT)
      output.append(graph.aromatic[source] ? ";a" : ";A");
    if (chemistry.matchFormalCharge) {
      int charge = graph.formalCharge[source];
      output.append(';');
      if (charge >= 0) output.append('+');
      output.append(charge);
    }
    if (chemistry.ringMatchesRingOnly) output.append(graph.ring[source] ? ";R" : ";!R");
    output.append(']');
    for (Closure closure : closures[atom]) {
      if (closure.first == atom) {
        int number = 1;
        while (number < open.length && open[number]) number++;
        if (number == open.length)
          throw new IllegalArgumentException("SMARTS output needs more than 99 simultaneously open ring closures");
        closure.number = number;
        open[number] = true;
        output.append(bond(closure.first, closure.second));
      } else open[closure.number] = false;
      if (closure.number >= 10) output.append('%');
      output.append(closure.number);
    }
  }

  private String bond(int first, int second) {
    first = atoms[first];
    second = atoms[second];
    int order = graph.bondOrder(first, second);
    boolean aromatic = graph.bondAromatic(first, second);
    if (order < 1 || order > 4 || aromatic != (order == 4))
      throw new IllegalArgumentException("Portable SMARTS requires consistent aromatic bond flags and orders");
    boolean strictOrder = chemistry.matchBondOrder == ChemOptions.BondOrderMode.STRICT;
    String result;
    if (chemistry.aromaticityMode == ChemOptions.AromaticityMode.STRICT)
      result = aromatic ? ":" : (strictOrder ? symbol(order) : "~") + ";!:";
    else if (!strictOrder) result = "~";
    else if (aromatic) result = "-,=,:";
    else {
      if (order <= 2 && graph.aromatic[first] && graph.aromatic[second])
        throw new IllegalArgumentException("Portable SMARTS cannot represent aromatic-endpoint resonance matching");
      result = symbol(order) + (order <= 2 ? ",:" : "");
    }
    if (chemistry.ringMatchesRingOnly) result += graph.bondInRing(first, second) ? ";@" : ";!@";
    return result;
  }

  private static String symbol(int order) {
    return order == 1 ? "-" : order == 2 ? "=" : order == 3 ? "#" : ":";
  }
}

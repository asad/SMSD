/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.Bond;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import org.openscience.cdk.stereo.TetrahedralChirality;

/**
 * The graph and the molecule of a match: the selected atoms and bonds of a source graph as a
 * {@link MolGraph} view for further matching, or as a CDK molecule with cut valences capped,
 * and the connected fragments of a selection. {@link SearchResult.Match} and the searches
 * delegate here.
 */
final class MatchGraph {
  private MatchGraph() {}

  /** Convert a complete graph without changing its original atom indices. */
  static IAtomContainer toMolecule(MolGraph graph) {
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int atom = 0; atom < graph.n; atom++) atoms.put(atom, atom);
    // Every bond is selected, so none is looked up.
    return materialise(graph, atoms, null);
  }

  /** Materialise exactly the selected common edges and cap cut valences. */
  static IAtomContainer toMolecule(
      MolGraph g, Map<Integer, Integer> mapping, Set<SearchResult.Bond> selected) {
    return materialise(g, mapping, Objects.requireNonNull(selected, "selected"));
  }

  /** The molecule of the selected bonds, or of every bond between mapped atoms when selected is null. */
  private static IAtomContainer materialise(
      MolGraph g, Map<Integer, Integer> mapping, Set<SearchResult.Bond> selected) {
    Objects.requireNonNull(g, "query");
    Objects.requireNonNull(mapping, "mapping");
    int[] atoms = mapping.keySet().stream().mapToInt(Integer::intValue).sorted().toArray();
    Map<Integer, Integer> reindex = new HashMap<>();
    IAtomContainer fragment = new AtomContainer();
    Map<IChemObject, IChemObject> objects = new HashMap<>();
    boolean[] allLigandsRetained = new boolean[atoms.length];
    for (int atom = 0; atom < atoms.length; atom++) {
      int original = atoms[atom];
      if (original < 0 || original >= g.n)
        throw new IndexOutOfBoundsException("MCS atom index out of range");
      reindex.put(original, atom);
      IAtom retained = new Atom(g.atomicNum[original]);
      retained.setFormalCharge(g.formalCharge[original]);
      if (g.massNumber[original] != 0)
        retained.setMassNumber(g.massNumber[original]);
      retained.setImplicitHydrogenCount(g.hydrogenCount(original));
      if (g.mol != null) {
        IAtom source = g.mol.getAtom(original);
        if (source.getPoint2d() != null)
          retained.setPoint2d(new Point2d(source.getPoint2d()));
        if (source.getPoint3d() != null)
          retained.setPoint3d(new Point3d(source.getPoint3d()));
      }
      fragment.addAtom(retained);
      if (g.mol != null)
        objects.put(g.mol.getAtom(original), fragment.getAtom(atom));
    }
    for (int atom = 0; atom < atoms.length; atom++) {
      int original = atoms[atom];
      allLigandsRetained[atom] = true;
      for (int neighbor : g.neighbors[original]) {
        Integer retainedNeighbor = reindex.get(neighbor);
        if (retainedNeighbor == null
            || (selected != null && !selected.contains(new SearchResult.Bond(original, neighbor)))) {
          allLigandsRetained[atom] = false;
          continue;
        }
        if (atom >= retainedNeighbor)
          continue;
        int order = g.bondOrder(original, neighbor);
        IBond retained = new Bond(fragment.getAtom(atom),
            fragment.getAtom(retainedNeighbor),
            order == 2       ? IBond.Order.DOUBLE
                : order == 3 ? IBond.Order.TRIPLE
                             : IBond.Order.SINGLE);
        retained.setIsAromatic(g.bondAromatic(original, neighbor));
        fragment.addBond(retained);
        if (g.mol != null)
          objects.put(g.mol.getBond(g.mol.getAtom(original), g.mol.getAtom(neighbor)),
              fragment.getBond(fragment.getBondCount() - 1));
      }
    }
    RingPerception.markRingAtomsAndBonds(fragment);
    for (IBond bond : fragment.bonds()) {
      boolean aromatic = bond.isAromatic() && bond.isInRing();
      bond.setIsAromatic(aromatic);
      if (aromatic)
        for (IAtom atom : bond.atoms()) atom.setIsAromatic(true);
    }
    for (int atom = 0; atom < atoms.length; atom++) {
      IAtom retained = fragment.getAtom(atom);
      int original = atoms[atom];
      if (allLigandsRetained[atom] && retained.isAromatic() == g.aromatic[original])
        continue;
      int hydrogens;
      if (retained.isAromatic() && (g.atomicNum[original] == 7 || g.atomicNum[original] == 15)) {
        hydrogens =
            g.hydrogenCount(original) + g.neighbors[original].length - retained.getBondCount();
      } else {
        int bondOrderSum = 0;
        for (IBond bond : retained.bonds())
          bondOrderSum += bond.isAromatic() ? 1 : bond.getOrder().numeric();
        hydrogens = MolGraph.defaultHydrogenCount(
            g.atomicNum[original], retained.isAromatic(), bondOrderSum, g.formalCharge[original]);
      }
      retained.setImplicitHydrogenCount(hydrogens);
    }
    boolean constrainedAromaticStereo = false;
    for (IBond bond : fragment.bonds()) {
      int first = fragment.indexOf(bond.getBegin()), second = fragment.indexOf(bond.getEnd());
      if (g.bondOrder(atoms[first], atoms[second]) == 4
          && g.dbStereo(atoms[first], atoms[second]) != 0
          && allLigandsRetained[first] && allLigandsRetained[second]) {
        bond.setOrder(IBond.Order.DOUBLE);
        constrainedAromaticStereo = true;
      }
    }
    try {
      Kekulization.kekulize(fragment);
    } catch (CDKException invalidAromaticFragment) {
      if (constrainedAromaticStereo)
        throw new IllegalArgumentException("Cannot retain aromatic double-bond stereochemistry",
            invalidAromaticFragment);
      for (IBond bond : fragment.bonds()) {
        if (bond.isAromatic())
          bond.setOrder(IBond.Order.SINGLE);
        bond.setIsAromatic(false);
      }
      for (IAtom atom : fragment.atoms()) {
        if (!atom.isAromatic())
          continue;
        atom.setIsAromatic(false);
        atom.setImplicitHydrogenCount(MolGraph.defaultHydrogenCount(atom.getAtomicNumber(), false,
            (int) fragment.getBondOrderSum(atom), atom.getFormalCharge()));
      }
    }
    if (g.mol != null) {
      for (IStereoElement<?, ?> stereo : g.mol.stereoElements()) {
        IStereoElement<?, ?> retained = stereo.mapStrict(objects);
        if (retained != null)
          fragment.addStereoElement(retained);
      }
    } else {
      for (int atom = 0; atom < atoms.length; atom++) {
        if (!allLigandsRetained[atom] || g.tetraChirality == null || g.tetraChirality[atoms[atom]] == 0)
          continue;
        int[] ligands = g.tetrahedralLigands(atoms[atom]);
        if (ligands == null)
          continue;
        IAtom[] retainedLigands = new IAtom[4];
        for (int ligand = 0; ligand < 4; ligand++)
          retainedLigands[ligand] =
              fragment.getAtom(ligands[ligand] < 0 ? atom : reindex.get(ligands[ligand]));
        fragment.addStereoElement(new TetrahedralChirality(
            fragment.getAtom(atom), retainedLigands,
            g.tetraChirality[atoms[atom]] == 2 ? ITetrahedralChirality.Stereo.CLOCKWISE
                                               : ITetrahedralChirality.Stereo.ANTI_CLOCKWISE));
      }
      for (IBond bond : fragment.bonds()) {
        int first = fragment.indexOf(bond.getBegin()), second = fragment.indexOf(bond.getEnd());
        int configuration = g.dbStereo(atoms[first], atoms[second]);
        if (configuration == 0 || !allLigandsRetained[first] || !allLigandsRetained[second])
          continue;
        IBond[] references = new IBond[2];
        for (int side = 0; side < 2; side++) {
          int end = side == 0 ? first : second, other = side == 0 ? second : first;
          for (int neighbor : g.neighbors[atoms[end]]) {
            if (neighbor == atoms[other])
              continue;
            references[side] =
                fragment.getBond(fragment.getAtom(end), fragment.getAtom(reindex.get(neighbor)));
            break;
          }
        }
        if (references[0] != null && references[1] != null)
          fragment.addStereoElement(
              new DoubleBondStereochemistry(bond, references,
                  configuration == 1 ? IDoubleBondStereochemistry.Conformation.TOGETHER
                                     : IDoubleBondStereochemistry.Conformation.OPPOSITE));
      }
    }
    for (IStereoElement<?, ?> stereo : fragment.stereoElements())
      if (stereo instanceof IDoubleBondStereochemistry
          && ((IDoubleBondStereochemistry) stereo).getStereoBond().getOrder() != IBond.Order.DOUBLE)
        throw new IllegalArgumentException("Double-bond stereochemistry requires a double bond");
    return fragment;
  }

  /**
   * Configures each double bond of a molecule made from the given source atoms whose ends both keep
   * a ligand but which lost a reference ligand of its source configuration: the other ligand of
   * that end, on the opposite side of the bond, becomes the reference. The molecule holds the
   * source atoms in ascending order, as {@link #toMolecule(MolGraph, Map, Set)} makes it.
   */
  static void retainDoubleBondStereo(MolGraph g, Set<Integer> kept, IAtomContainer fragment) {
    if (g.mol == null)
      return;
    Map<IAtom, IAtom> atoms = new HashMap<>();
    int index = 0;
    for (int atom : new TreeSet<>(kept)) atoms.put(g.mol.getAtom(atom), fragment.getAtom(index++));
    Set<IBond> configured = new HashSet<>();
    for (IStereoElement<?, ?> element : fragment.stereoElements())
      if (element instanceof IDoubleBondStereochemistry)
        configured.add(((IDoubleBondStereochemistry) element).getStereoBond());
    List<IStereoElement<?, ?>> added = new ArrayList<>();
    for (IStereoElement<?, ?> element : g.mol.stereoElements()) {
      if (!(element instanceof IDoubleBondStereochemistry))
        continue;
      IDoubleBondStereochemistry stereo = (IDoubleBondStereochemistry) element;
      IBond source = stereo.getStereoBond();
      IAtom begin = atoms.get(source.getBegin()), end = atoms.get(source.getEnd());
      IBond bond = begin == null || end == null ? null : fragment.getBond(begin, end);
      if (bond == null || configured.contains(bond) || bond.getOrder() != IBond.Order.DOUBLE)
        continue;
      boolean together = stereo.getStereo() == IDoubleBondStereochemistry.Conformation.TOGETHER;
      IBond[] references = new IBond[2];
      for (IBond ligand : stereo.getBonds()) {
        IAtom centre = ligand.contains(source.getBegin()) ? source.getBegin() : source.getEnd();
        IAtom retainedCentre = atoms.get(centre);
        int side = retainedCentre == bond.getBegin() ? 0 : 1;
        IAtom retainedLigand = atoms.get(ligand.getOther(centre));
        if (retainedLigand != null) {
          references[side] = fragment.getBond(retainedCentre, retainedLigand);
          continue;
        }
        int ligands = 0;
        for (IBond other : g.mol.getConnectedBondsList(centre)) {
          if (other.equals(source) || other.equals(ligand))
            continue;
          ligands++;
          IAtom retainedOther = atoms.get(other.getOther(centre));
          if (retainedOther != null)
            references[side] = fragment.getBond(retainedCentre, retainedOther);
        }
        if (ligands != 1) references[side] = null;
        together = !together;
      }
      if (references[0] != null && references[1] != null)
        added.add(new DoubleBondStereochemistry(bond, references,
            together ? IDoubleBondStereochemistry.Conformation.TOGETHER
                     : IDoubleBondStereochemistry.Conformation.OPPOSITE));
    }
    for (IStereoElement<?, ?> element : added) fragment.addStereoElement(element);
  }

  /** Preserve source chemistry while selecting atoms and common edges for further matching. */
  static MolGraph toGraph(
      MolGraph graph, Map<Integer, Integer> atoms, Set<SearchResult.Bond> bonds) {
    int[] original = atoms.keySet().stream().mapToInt(Integer::intValue).sorted().toArray();
    Map<Integer, Integer> reindex = new HashMap<>();
    for (int index = 0; index < original.length; index++) reindex.put(original[index], index);
    int size = original.length;
    int[] elements = new int[size], charges = new int[size], masses = new int[size],
          hydrogens = new int[size], tetra = new int[size];
    boolean[] rings = new boolean[size], aromatic = new boolean[size];
    int[][] adjacency = new int[size][], orders = new int[size][], stereo = new int[size][];
    boolean[][] ringBonds = new boolean[size][], aromaticBonds = new boolean[size][];
    boolean hasDoubleBondStereo = false;
    for (int index = 0; index < size; index++) {
      int source = original[index];
      elements[index] = graph.atomicNum[source];
      charges[index] = graph.formalCharge[source];
      masses[index] = graph.massNumber[source];
      hydrogens[index] = graph.hydrogenCount(source);
      rings[index] = graph.ring[source];
      aromatic[index] = graph.aromatic[source];
      List<Integer> neighbors = new ArrayList<>();
      for (int neighbor : graph.neighbors[source])
        if (reindex.containsKey(neighbor)
            && bonds.contains(new SearchResult.Bond(source, neighbor)))
          neighbors.add(neighbor);
      int count = neighbors.size();
      adjacency[index] = new int[count];
      orders[index] = new int[count];
      ringBonds[index] = new boolean[count];
      aromaticBonds[index] = new boolean[count];
      for (int edge = 0; edge < count; edge++) {
        int neighbor = neighbors.get(edge);
        adjacency[index][edge] = reindex.get(neighbor);
        orders[index][edge] = graph.bondOrder(source, neighbor);
        ringBonds[index][edge] = graph.bondInRing(source, neighbor);
        aromaticBonds[index][edge] = graph.bondAromatic(source, neighbor);
        boolean ligandsRetained = count == graph.neighbors[source].length;
        for (int ligand : graph.neighbors[neighbor])
          if (!reindex.containsKey(ligand)
              || !bonds.contains(new SearchResult.Bond(neighbor, ligand)))
            ligandsRetained = false;
        int configuration = ligandsRetained ? graph.dbStereo(source, neighbor) : 0;
        if (configuration != 0) {
          if (stereo[index] == null)
            stereo[index] = new int[size];
          stereo[index][reindex.get(neighbor)] = configuration;
          hasDoubleBondStereo = true;
        }
      }
      int[] reference = graph.tetrahedralLigands(source);
      if (graph.tetraChirality != null && graph.tetraChirality[source] != 0 && reference != null
          && count == graph.neighbors[source].length) {
        // The view reads its ligands as MolGraph.tetrahedralLigands does: an implicit ligand
        // comes second when the first neighbour precedes the atom, and first otherwise.
        List<Integer> outputOrder = new ArrayList<>(neighbors);
        if (count == 3)
          outputOrder.add(adjacency[index][0] < index ? 1 : 0, -1);
        boolean odd = false;
        for (int first = 0; first < reference.length; first++)
          for (int second = first + 1; second < reference.length; second++)
            if (outputOrder.indexOf(reference[first]) > outputOrder.indexOf(reference[second]))
              odd = !odd;
        int tag = graph.tetraChirality[source];
        tetra[index] = odd ? 3 - tag : tag;
      }
    }
    MolGraph view = new MolGraph.Builder()
                        .atomCount(size)
                        .atomicNumbers(elements)
                        .formalCharges(charges)
                        .massNumbers(masses)
                        .hydrogenCounts(hydrogens)
                        .ringFlags(rings)
                        .aromaticFlags(aromatic)
                        .neighbors(adjacency)
                        .bondOrders(orders)
                        .bondRingFlags(ringBonds)
                        .bondAromaticFlags(aromaticBonds)
                        .tetrahedralChirality(tetra)
                        .doubleBondStereo(hasDoubleBondStereo ? stereo : null)
                        .build();
    graph.ensureRingCounts();
    view.ensureRingCounts();
    for (int index = 0; index < size; index++) {
      view.ringCount[index] = graph.ringCount[original[index]];
      view.tautomerClass[index] = graph.tautomerClass[original[index]];
      view.tautomerWeight[index] = graph.tautomerWeight[original[index]];
      view.tautomerSite[index] = graph.tautomerSite[original[index]];
      view.tautomerPKa[index] = graph.tautomerPKa[original[index]];
    }
    return view;
  }

  /** Every bond of {@code graph} between mapped atoms, each mapped to the bond between the images. */
  static Map<SearchResult.Bond, SearchResult.Bond> allBonds(
      MolGraph graph, Map<Integer, Integer> atoms) {
    Map<SearchResult.Bond, SearchResult.Bond> result = new LinkedHashMap<>();
    for (int queryAtom : atoms.keySet())
      for (int neighbor : graph.neighbors[queryAtom])
        if (queryAtom < neighbor && atoms.containsKey(neighbor))
          result.put(new SearchResult.Bond(queryAtom, neighbor),
              new SearchResult.Bond(atoms.get(queryAtom), atoms.get(neighbor)));
    return result;
  }

  /** The connected components of the mapped atoms under the selected bonds, by smallest atom. */
  static List<Set<Integer>> fragments(
      Map<Integer, Integer> atoms, Collection<SearchResult.Bond> bonds) {
    Map<Integer, List<Integer>> adjacency = new TreeMap<>();
    for (int atom : atoms.keySet()) adjacency.put(atom, new ArrayList<>());
    for (SearchResult.Bond bond : bonds) {
      adjacency.get(bond.first()).add(bond.second());
      adjacency.get(bond.second()).add(bond.first());
    }
    Set<Integer> visited = new HashSet<>();
    List<Set<Integer>> result = new ArrayList<>();
    for (int root : adjacency.keySet()) {
      if (!visited.add(root))
        continue;
      Set<Integer> component = new TreeSet<>();
      ArrayDeque<Integer> queue = new ArrayDeque<>();
      queue.add(root);
      while (!queue.isEmpty()) {
        int atom = queue.remove();
        component.add(atom);
        for (int neighbor : adjacency.get(atom))
          if (visited.add(neighbor))
            queue.add(neighbor);
      }
      result.add(component);
    }
    return result;
  }
}

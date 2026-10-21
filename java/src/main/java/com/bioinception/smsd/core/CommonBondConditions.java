/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * What a mapping of {@link CommonBondSearch} must satisfy beyond its common bonds: the pairs and
 * bonds of a required mapping, the size and number of its fragments, and, where the options ask
 * for them, mapped stereochemistry and whole rings, all as {@link MCS#validateMapping} checks
 * them. The search maps atoms; a mapping here selects some of the common bonds of its atoms, and
 * its score is the number of bonds it selects, then the number of atoms. Without stereochemistry
 * and rings a mapping selects all its common bonds or none. Stereochemistry and rings are read
 * from the selected bonds alone, so a mapping may leave a common bond out to keep a centre whose
 * ligands disagree unresolved, or to keep only whole rings. The common bonds of the atoms bound
 * every selection of them, so the bounds of the search hold.
 *
 * <p>The rings are met by the selections in which each selected bond on a ring of either
 * molecule lies on a wholly selected relevant ring of that molecule. Such selections are closed
 * under union, so each set of bonds holds a largest one, found by taking off bonds on no whole
 * ring until none is left ({@link #keepWholeRings}). A stereo centre is met once at most two of
 * its ligands are selected, or its ligands agree; a stereo double bond once it is left out, or
 * one end keeps no selected ligand, or its ligands agree. Taking bonds away never breaks either,
 * and every other condition holds of a selection once it holds of a smaller one. So each best
 * selection is the largest whole-ring part of the common bonds less a least set of bonds that
 * resolves the conflicts, and {@link #select} tries these sets depth first, one conflict at a
 * time, each resolved in every least way.
 */
final class CommonBondConditions {

  /** A conflict at a tetrahedral centre, see {@link Conflict}. */
  private static final int CENTRE = 0;
  /** A conflict at a double bond. */
  private static final int DOUBLE_BOND = 1;

  /** The graphs of both molecules and their compatible pairs. */
  private final CommonBondGraphs graph;
  private final CommonBondGraphs.Side query, target;
  /** The molecules, for their stereochemistry and rings. */
  private final MolGraph queryGraph, targetGraph;
  /** Counts the work, so that a long selection can be stopped. */
  private final CommonBondGraphs.Clock clock;
  /** Whether the selected bonds form a single fragment. */
  private final boolean connected;
  /** The fewest atoms of a fragment and the most fragments. */
  private final int minFragmentAtoms, maxFragments;
  /** The target of each query atom that a required mapping fixes, -1 for none; null without one. */
  private final int[] requiredTarget;
  /** The query bonds a required mapping selects; null without one. */
  private final boolean[] requiredBond;
  /** Whether chirality, double-bond stereochemistry and whole rings are checked. */
  private final boolean chirality, bondStereo, rings;

  /** The bonds of each relevant ring of each molecule, by the bond numbers of its side. */
  private int[][] queryCycles, targetCycles;
  /** The rings through each bond and through each atom of each molecule. */
  private int[][] queryBondCycles, targetBondCycles, queryAtomCycles, targetAtomCycles;
  /** Whether each bond and atom of each molecule is on a ring, as the ring check reads it. */
  private boolean[] queryRingBond, targetRingBond, queryRingAtom, targetRingAtom;
  /** The bonds of each ring left out of the selection being checked. */
  private int[] queryMissing, targetMissing;

  /** The mapping being checked, and its number of atoms. */
  private int[] mapping;
  private int atoms;
  /** The target bond of each common bond of the mapping, by query bond. */
  private final int[] image;
  /** The query bond whose image each target bond is, while it is selected, or -1. */
  private final int[] preimage;
  /** The common bonds, the first {@link #commonCount}, in query bond order. */
  private final int[] commonBonds;
  private int commonCount;
  /** The selection under way, and the selection of a leaf once its rings are whole. */
  private final boolean[] selected, kept;
  /** Bonds still to check while rings are made whole. */
  private int[] queue;
  /** Union-find of the query atoms over the selected bonds, and the size of each set. */
  private final int[] parent, setSize;
  /** The conflicts of the selection, in the order they are resolved. */
  private final List<Conflict> conflicts = new ArrayList<>();
  /** The fewest bonds a selection needs, and the most selected so far. */
  private int leastBonds, bestBonds;
  /** The most selections wanted; ties with the best are kept until there are as many. */
  private int room;
  /** The best selections so far, as records, see {@link #record}. */
  private final List<int[]> best = new ArrayList<>();
  /** The selected bonds that the conflicts may still give up, and their number. */
  private int upperBonds;
  /** The stamp of each query bond counted by {@link #fewestToLeaveOut} or {@link #persistentLoss}, and the stamp under way. */
  private final int[] bondMark;
  private int markStamp;
  /** The query atoms with a tetrahedral tag where chirality is checked, and none otherwise. */
  private final int[] centres;

  /** A stereo centre or double bond whose selected ligands disagree. */
  private static final class Conflict {
    final int kind, atom, other, bond;
    /** The selected bonds to the ligands of a centre, or to those of each end of a double bond. */
    final int[] ligands, otherLigands;

    Conflict(int kind, int atom, int other, int bond, int[] ligands, int[] otherLigands) {
      this.kind = kind;
      this.atom = atom;
      this.other = other;
      this.bond = bond;
      this.ligands = ligands;
      this.otherLigands = otherLigands;
    }
  }

  private CommonBondConditions(MolGraph queryGraph, MolGraph targetGraph, CommonBondGraphs graph,
      MCS.Options options, boolean connected, SearchResult.Match required, boolean chirality,
      boolean bondStereo, boolean rings, CommonBondGraphs.Clock clock) {
    this.graph = graph;
    this.query = graph.query;
    this.target = graph.target;
    this.queryGraph = queryGraph;
    this.targetGraph = targetGraph;
    this.clock = clock;
    this.connected = connected;
    this.minFragmentAtoms = options.minFragmentAtoms;
    this.maxFragments = options.maxFragments;
    this.chirality = chirality;
    this.bondStereo = bondStereo;
    this.rings = rings;
    image = new int[query.bondCount];
    preimage = new int[target.bondCount];
    Arrays.fill(preimage, -1);
    commonBonds = new int[query.bondCount];
    selected = new boolean[query.bondCount];
    kept = new boolean[query.bondCount];
    bondMark = new int[query.bondCount];
    centres = chirality ? taggedAtoms(queryGraph) : new int[0];
    queue = new int[Math.max(16, query.bondCount)];
    parent = new int[query.atomCount];
    setSize = new int[query.atomCount];
    if (required != null && !required.atomMapping().isEmpty()) {
      requiredTarget = new int[query.atomCount];
      Arrays.fill(requiredTarget, -1);
      for (Map.Entry<Integer, Integer> pair : required.atomMapping().entrySet())
        requiredTarget[pair.getKey()] = pair.getValue();
      requiredBond = new boolean[query.bondCount];
      for (SearchResult.Bond bond : required.bondMapping().keySet()) {
        int queryBond = query.bondBetween(bond.first(), bond.second());
        if (queryBond >= 0)
          requiredBond[queryBond] = true;
      }
    } else {
      requiredTarget = null;
      requiredBond = null;
    }
    if (rings) {
      queryCycles = cycles(query, queryGraph);
      targetCycles = cycles(target, targetGraph);
      queryBondCycles = through(queryCycles, query.bondCount, null);
      targetBondCycles = through(targetCycles, target.bondCount, null);
      queryAtomCycles = through(queryCycles, query.atomCount, query);
      targetAtomCycles = through(targetCycles, target.atomCount, target);
      queryRingBond = ringBonds(query, queryGraph, queryBondCycles);
      targetRingBond = ringBonds(target, targetGraph, targetBondCycles);
      queryRingAtom = ringAtoms(queryGraph, queryAtomCycles);
      targetRingAtom = ringAtoms(targetGraph, targetAtomCycles);
      queryMissing = new int[queryCycles.length];
      targetMissing = new int[targetCycles.length];
    }
  }

  /**
   * The conditions of a search of {@code queryGraph} in {@code targetGraph} under
   * {@code options}, extending {@code required} unless it is null, or null when every mapping
   * with common bonds meets them. Chirality and double-bond stereochemistry count only when both
   * molecules have them, and rings only when either has one.
   */
  static CommonBondConditions of(MolGraph queryGraph, MolGraph targetGraph, CommonBondGraphs graph,
      MCS.Options options, boolean connected, SearchResult.Match required,
      CommonBondGraphs.Clock clock) {
    ChemOptions chemistry = options.chemistry;
    boolean chirality = chemistry.useChirality && hasCentres(queryGraph) && hasCentres(targetGraph);
    boolean bondStereo = chemistry.useBondStereo && hasStereoBonds(queryGraph)
        && hasStereoBonds(targetGraph);
    boolean rings = chemistry.completeRingsOnly && (hasRings(queryGraph) || hasRings(targetGraph));
    boolean sizes = options.minFragmentAtoms > 2
        || (!connected && options.maxFragments < Integer.MAX_VALUE);
    boolean seeded = required != null && !required.atomMapping().isEmpty();
    if (!chirality && !bondStereo && !rings && !sizes && !seeded)
      return null;
    return new CommonBondConditions(queryGraph, targetGraph, graph, options, connected, required,
        chirality, bondStereo, rings, clock);
  }

  /** The atoms of {@code graph} with a tetrahedral tag, in ascending order. */
  private static int[] taggedAtoms(MolGraph graph) {
    int count = 0;
    for (int tag : graph.tetraChirality) if (tag != 0) count++;
    int[] atoms = new int[count];
    count = 0;
    for (int atom = 0; atom < graph.tetraChirality.length; atom++)
      if (graph.tetraChirality[atom] != 0)
        atoms[count++] = atom;
    return atoms;
  }

  /**
   * A least number of the common bonds of {@code mapping} that the best selection of every
   * mapping extending it leaves out for its tetrahedral centres, so that the search can bound
   * the score below it. A mapped centre with at least three ligands on common bonds that
   * disagree with its image disagrees however the mapping grows, as any three ligands fix the
   * place of the fourth, so every selection leaves out all but two of these bonds. Centres that
   * share none of them leave out different bonds.
   */
  int persistentLoss(int[] mapping) {
    if (centres.length == 0)
      return 0;
    clock.tick(centres.length);
    if (++markStamp == Integer.MAX_VALUE) {
      Arrays.fill(bondMark, 0);
      markStamp = 1;
    }
    int total = 0;
    for (int atom : centres) {
      if (mapping[atom] < 0)
        continue;
      int ligands = 0;
      boolean shared = false;
      for (int ligand : queryGraph.neighbors[atom]) {
        int bond = commonBond(atom, ligand, mapping);
        if (bond >= 0) {
          ligands++;
          shared |= bondMark[bond] == markStamp;
        }
      }
      if (ligands < 3 || shared || !centreDisagrees(atom, mapping, null))
        continue;
      total += ligands - 2;
      for (int ligand : queryGraph.neighbors[atom]) {
        int bond = commonBond(atom, ligand, mapping);
        if (bond >= 0)
          bondMark[bond] = markStamp;
      }
    }
    return total;
  }

  /** The query bond between {@code atom} and {@code ligand} if it is a common bond of {@code mapping}, and -1 otherwise. */
  private int commonBond(int atom, int ligand, int[] mapping) {
    if (mapping[ligand] < 0)
      return -1;
    int bond = query.bondBetween(atom, ligand);
    if (bond < 0)
      return -1;
    int targetBond = target.bondBetween(mapping[atom], mapping[ligand]);
    return targetBond >= 0 && graph.bondsMatch(bond, targetBond) ? bond : -1;
  }

  /** Whether a mapping may select fewer than all its common bonds. */
  boolean dropsBonds() {
    return chirality || bondStereo || rings;
  }

  /**
   * Whether exchanging twin terminal atoms takes each mapping onto one that meets the conditions
   * as well, for twins that carry no stereo tag that the conditions read and whose shared atom
   * carries none either, see {@link #untagged}.
   */
  boolean keepsTwins() {
    return requiredTarget == null;
  }

  /**
   * Whether {@code atom}, of the query when {@code onQuery} and otherwise of the target, carries
   * no tetrahedral tag that the chirality check reads and no stereo double bond that the
   * double-bond check reads. Only the tags of an atom and of its bonds tie the order of its
   * neighbours to the conditions, so exchanging two terminal neighbours of such an atom, neither
   * tagged, leaves every check as it was.
   */
  boolean untagged(boolean onQuery, int atom) {
    MolGraph molecule = onQuery ? queryGraph : targetGraph;
    if (chirality && molecule.tetraChirality != null && molecule.tetraChirality[atom] != 0)
      return false;
    if (bondStereo)
      for (int neighbor : molecule.neighbors[atom])
        if (molecule.dbStereo(atom, neighbor) != 0)
          return false;
    return true;
  }

  /** Whether adding a pendant pair to a mapping that meets the conditions gives one that meets them. */
  boolean keepsPendantPairs() {
    return !chirality && !bondStereo;
  }

  /** The fewest atoms of a fragment. */
  int minFragmentAtoms() {
    return minFragmentAtoms;
  }

  /** The most fragments. */
  int maxFragments() {
    return maxFragments;
  }

  /**
   * The best selections of the common bonds of {@code mapping} that meet the conditions, if their
   * score is at least {@code least}: at most {@code room} of them are added to {@code out}, as
   * records, and their score is returned. A record is the mapping, followed, when it leaves
   * common bonds out, by their number and the query bonds left out in ascending order. Returns
   * -1, adding nothing, when no selection scores {@code least}.
   */
  long select(int[] mapping, long least, int room, List<int[]> out) {
    this.mapping = mapping;
    atoms = 0;
    for (int queryAtom = 0; queryAtom < query.atomCount; queryAtom++)
      if (mapping[queryAtom] >= 0)
        atoms++;
    if (requiredTarget != null)
      for (int queryAtom = 0; queryAtom < query.atomCount; queryAtom++)
        if (requiredTarget[queryAtom] >= 0 && mapping[queryAtom] != requiredTarget[queryAtom])
          return -1;
    commonCount = 0;
    clock.tick(query.bondCount);
    for (int queryBond = 0; queryBond < query.bondCount; queryBond++) {
      int first = mapping[query.bondFirst[queryBond]], second = mapping[query.bondSecond[queryBond]];
      int targetBond = first < 0 || second < 0 ? -1 : target.bondBetween(first, second);
      boolean common = targetBond >= 0 && graph.bondsMatch(queryBond, targetBond);
      image[queryBond] = common ? targetBond : -1;
      selected[queryBond] = common;
      if (common)
        commonBonds[commonCount++] = queryBond;
    }
    leastBonds = atoms >= CommonBondGraphs.atoms(least) ? CommonBondGraphs.bonds(least)
                                                        : CommonBondGraphs.bonds(least) + 1;
    if (commonCount < leastBonds)
      return -1;
    if (rings)
      keepWholeRings(selected);
    upperBonds = 0;
    for (int index = 0; index < commonCount; index++)
      if (selected[commonBonds[index]])
        upperBonds++;
    if (upperBonds < leastBonds)
      return -1;
    best.clear();
    bestBonds = -1;
    this.room = room;
    if (chirality || bondStereo) {
      findConflicts();
      resolve(0, 0);
    } else {
      leaf();
    }
    if (best.isEmpty())
      return -1;
    for (int index = 0; index < best.size() && index < room; index++) out.add(best.get(index));
    return CommonBondGraphs.score(atoms, bestBonds);
  }

  /** The conflicts of the selection: the centres first, by query atom, then the double bonds. */
  private void findConflicts() {
    conflicts.clear();
    if (chirality)
      for (int atom = 0; atom < query.atomCount; atom++)
        if (mapping[atom] >= 0 && centreDisagrees(atom, selected)) {
          int[] order = queryGraph.tetrahedralLigands(atom);
          conflicts.add(new Conflict(CENTRE, atom, -1, -1, selectedBonds(atom, order, -1), null));
        }
    if (bondStereo)
      for (int atom = 0; atom < query.atomCount; atom++) {
        if (mapping[atom] < 0)
          continue;
        for (int other : queryGraph.neighbors[atom]) {
          int bond = query.bondBetween(atom, other);
          if (other > atom && bond >= 0 && selected[bond] && doubleBondDisagrees(atom, other, selected))
            conflicts.add(new Conflict(DOUBLE_BOND, atom, other, bond,
                selectedBonds(atom, queryGraph.neighbors[atom], other),
                selectedBonds(other, queryGraph.neighbors[other], atom)));
        }
      }
    clock.tick(query.atomCount + conflicts.size());
  }

  /** The selected bonds from {@code atom} to the atoms of {@code ligands} other than {@code skipped}, in ascending order. */
  private int[] selectedBonds(int atom, int[] ligands, int skipped) {
    int[] bonds = new int[ligands.length];
    int count = 0;
    for (int ligand : ligands) {
      if (ligand < 0 || ligand == skipped)
        continue;
      int bond = query.bondBetween(atom, ligand);
      if (bond >= 0 && selected[bond])
        bonds[count++] = bond;
    }
    int[] result = Arrays.copyOf(bonds, count);
    Arrays.sort(result);
    return result;
  }

  /**
   * Resolves the conflicts from {@code from} on, with {@code dropped} bonds left out so far:
   * the first conflict not yet resolved is resolved in each least way, and each way is followed
   * on. Once none is left, the selection is a leaf.
   */
  private void resolve(int from, int dropped) {
    clock.tick();
    if (upperBonds - dropped - fewestToLeaveOut(from) < leastBonds || !fragmentsHold(selected))
      return;
    int open = from;
    while (open < conflicts.size() && resolved(conflicts.get(open))) open++;
    if (open == conflicts.size()) {
      leaf();
      return;
    }
    Conflict conflict = conflicts.get(open);
    if (conflict.kind == CENTRE) {
      // A centre with at most two selected ligands is unresolved.
      int selectedLigands = 0, free = 0;
      int[] choices = new int[conflict.ligands.length];
      for (int bond : conflict.ligands) {
        if (!selected[bond])
          continue;
        selectedLigands++;
        if (requiredBond == null || !requiredBond[bond])
          choices[free++] = bond;
      }
      int need = selectedLigands - 2;
      if (need <= 0 || need > free)
        return;
      leaveOut(choices, free, need, 0, open, dropped);
    } else {
      // The double bond itself, or every selected ligand of one end.
      if (selected[conflict.bond] && !required(conflict.bond)) {
        selected[conflict.bond] = false;
        resolve(open, dropped + 1);
        selected[conflict.bond] = true;
      }
      dropAll(conflict.ligands, open, dropped);
      dropAll(conflict.otherLigands, open, dropped);
    }
  }

  /**
   * A least number of further bonds that the conflicts from {@code from} on leave out: a centre
   * keeps at most two selected ligands and a double bond loses at least one bond, and conflicts
   * that share no selected bond lose different bonds.
   */
  private int fewestToLeaveOut(int from) {
    if (++markStamp == Integer.MAX_VALUE) {
      Arrays.fill(bondMark, 0);
      markStamp = 1;
    }
    int total = 0;
    for (int index = from; index < conflicts.size(); index++) {
      Conflict conflict = conflicts.get(index);
      int need;
      boolean shared = false;
      if (conflict.kind == CENTRE) {
        int selectedLigands = 0;
        for (int bond : conflict.ligands)
          if (selected[bond]) {
            selectedLigands++;
            shared |= bondMark[bond] == markStamp;
          }
        need = selectedLigands - 2;
      } else {
        need = selected[conflict.bond] && anySelected(conflict.ligands)
                && anySelected(conflict.otherLigands) ? 1 : 0;
        shared = need > 0 && (bondMark[conflict.bond] == markStamp
            || anyMarked(conflict.ligands) || anyMarked(conflict.otherLigands));
      }
      if (need <= 0 || shared)
        continue;
      total += need;
      if (conflict.kind == CENTRE) {
        mark(conflict.ligands);
      } else {
        bondMark[conflict.bond] = markStamp;
        mark(conflict.ligands);
        mark(conflict.otherLigands);
      }
    }
    return total;
  }

  /** Whether any of {@code bonds} is selected. */
  private boolean anySelected(int[] bonds) {
    for (int bond : bonds) if (selected[bond]) return true;
    return false;
  }

  /** Whether any selected bond of {@code bonds} is marked for the bound under way. */
  private boolean anyMarked(int[] bonds) {
    for (int bond : bonds) if (selected[bond] && bondMark[bond] == markStamp) return true;
    return false;
  }

  /** Marks the selected bonds of {@code bonds} for the bound under way. */
  private void mark(int[] bonds) {
    for (int bond : bonds) if (selected[bond]) bondMark[bond] = markStamp;
  }

  /** Leaves out each choice of {@code need} of the first {@code count} bonds of {@code choices}, from {@code start} on. */
  private void leaveOut(int[] choices, int count, int need, int start, int open, int dropped) {
    if (need == 0) {
      resolve(open, dropped);
      return;
    }
    for (int index = start; index <= count - need; index++) {
      selected[choices[index]] = false;
      leaveOut(choices, count, need - 1, index + 1, open, dropped + 1);
      selected[choices[index]] = true;
    }
  }

  /** Leaves out every selected bond of {@code bonds}, unless one is required or none is selected. */
  private void dropAll(int[] bonds, int open, int dropped) {
    int[] dropping = new int[bonds.length];
    int count = 0;
    for (int bond : bonds) {
      if (!selected[bond])
        continue;
      if (required(bond))
        return;
      dropping[count++] = bond;
    }
    if (count == 0)
      return;
    for (int index = 0; index < count; index++) selected[dropping[index]] = false;
    resolve(open, dropped + count);
    for (int index = 0; index < count; index++) selected[dropping[index]] = true;
  }

  /** Whether a required mapping selects query bond {@code bond}. */
  private boolean required(int bond) {
    return requiredBond != null && requiredBond[bond];
  }

  /** Whether {@code conflict} is resolved by the selection under way. */
  private boolean resolved(Conflict conflict) {
    if (conflict.kind == CENTRE)
      return !centreDisagrees(conflict.atom, selected);
    return !selected[conflict.bond] || !doubleBondDisagrees(conflict.atom, conflict.other, selected);
  }

  /**
   * Checks the selection under way once its conflicts are resolved: with its rings made whole,
   * it is kept if it meets every condition and is no worse than the best so far.
   */
  private void leaf() {
    clock.tick();
    System.arraycopy(selected, 0, kept, 0, selected.length);
    if (rings)
      keepWholeRings(kept);
    int count = 0;
    for (int index = 0; index < commonCount; index++)
      if (kept[commonBonds[index]])
        count++;
    if (count < leastBonds || !valid(kept))
      return;
    int[] record = record(kept, count);
    if (count > bestBonds) {
      best.clear();
      bestBonds = count;
    } else {
      // A tie, while fewer than wanted are held, may repeat one held.
      clock.tick(best.size());
      for (int[] other : best)
        if (Arrays.equals(other, record))
          return;
    }
    best.add(record);
    // While there is room for ties no smaller selection is of use; once full, only a larger one.
    leastBonds = best.size() < room ? count : count + 1;
  }

  /** The record of the mapping with the selection {@code selection} of {@code count} common bonds. */
  private int[] record(boolean[] selection, int count) {
    int left = commonCount - count;
    int[] record = Arrays.copyOf(mapping, query.atomCount + (left > 0 ? left + 1 : 0));
    if (left > 0) {
      record[query.atomCount] = left;
      int place = query.atomCount + 1;
      for (int index = 0; index < commonCount; index++)
        if (!selection[commonBonds[index]])
          record[place++] = commonBonds[index];
    }
    return record;
  }

  /**
   * Whether the selection meets every condition other than the bonds on rings, which it meets
   * already: the required bonds are selected, every mapped atom has a selected bond, the
   * fragments are as many and as large as asked, the mapped ring atoms lie on whole rings, and
   * the mapped stereochemistry agrees.
   */
  private boolean valid(boolean[] selection) {
    if (!fragmentsHold(selection))
      return false;
    if (rings)
      for (int atom = 0; atom < query.atomCount; atom++) {
        if (mapping[atom] < 0)
          continue;
        if (queryRingAtom[atom] && !anyWhole(queryAtomCycles[atom], queryMissing)
            || targetRingAtom[mapping[atom]]
                && !anyWhole(targetAtomCycles[mapping[atom]], targetMissing))
          return false;
      }
    if (chirality)
      for (int atom = 0; atom < query.atomCount; atom++)
        if (mapping[atom] >= 0 && centreDisagrees(atom, selection))
          return false;
    if (bondStereo)
      for (int index = 0; index < commonCount; index++) {
        int bond = commonBonds[index];
        if (selection[bond]
            && doubleBondDisagrees(query.bondFirst[bond], query.bondSecond[bond], selection))
          return false;
      }
    return true;
  }

  /**
   * Whether the selection keeps every required bond and gives every mapped atom a selected bond,
   * and its fragments are as many and as large as asked. Once one of these fails, leaving more
   * bonds out does not mend it, so a selection that fails here leads to no valid one.
   */
  private boolean fragmentsHold(boolean[] selection) {
    clock.tick(query.atomCount + commonCount);
    for (int atom = 0; atom < query.atomCount; atom++) {
      parent[atom] = atom;
      setSize[atom] = 0;
    }
    for (int index = 0; index < commonCount; index++) {
      int bond = commonBonds[index];
      if (!selection[bond]) {
        if (required(bond))
          return false;
        continue;
      }
      int first = query.bondFirst[bond], second = query.bondSecond[bond];
      setSize[first] = setSize[second] = 1;
      int firstRoot = find(first), secondRoot = find(second);
      if (firstRoot != secondRoot)
        parent[firstRoot] = secondRoot;
    }
    for (int atom = 0; atom < query.atomCount; atom++)
      if (mapping[atom] >= 0 && setSize[atom] == 0)
        return false;
    int fragments = 0;
    for (int atom = 0; atom < query.atomCount; atom++) setSize[atom] = 0;
    for (int atom = 0; atom < query.atomCount; atom++)
      if (mapping[atom] >= 0 && setSize[find(atom)]++ == 0)
        fragments++;
    if (connected && fragments != 1 || fragments > maxFragments)
      return false;
    for (int atom = 0; atom < query.atomCount; atom++)
      if (mapping[atom] >= 0 && find(atom) == atom && setSize[atom] < minFragmentAtoms)
        return false;
    return true;
  }

  /** The root of the set of {@code atom}. */
  private int find(int atom) {
    return UnionFind.root(parent, atom);
  }

  /**
   * Takes off the selected bonds of {@code selection} that lie on a ring of either molecule but
   * on no wholly selected relevant ring of it, until none is left. What is left is the largest
   * selection within it whose ring bonds all lie on whole rings, as such selections are closed
   * under union. The counts of bonds missing from each ring are left as they are for it.
   */
  private void keepWholeRings(boolean[] selection) {
    clock.tick(commonCount + queryCycles.length + targetCycles.length);
    for (int index = 0; index < commonCount; index++) {
      int bond = commonBonds[index];
      preimage[image[bond]] = selection[bond] ? bond : -1;
    }
    for (int cycle = 0; cycle < queryCycles.length; cycle++) {
      int missing = 0;
      for (int bond : queryCycles[cycle]) if (!selection[bond]) missing++;
      queryMissing[cycle] = missing;
    }
    for (int cycle = 0; cycle < targetCycles.length; cycle++) {
      int missing = 0;
      for (int bond : targetCycles[cycle]) if (preimage[bond] < 0) missing++;
      targetMissing[cycle] = missing;
    }
    int size = 0;
    for (int index = 0; index < commonCount; index++)
      if (selection[commonBonds[index]])
        size = push(size, commonBonds[index]);
    while (size > 0) {
      int bond = queue[--size];
      clock.tick();
      if (!selection[bond] || onWholeRings(bond))
        continue;
      selection[bond] = false;
      preimage[image[bond]] = -1;
      for (int cycle : queryBondCycles[bond])
        if (queryMissing[cycle]++ == 0)
          for (int other : queryCycles[cycle]) if (selection[other]) size = push(size, other);
      for (int cycle : targetBondCycles[image[bond]])
        if (targetMissing[cycle]++ == 0)
          for (int other : targetCycles[cycle])
            if (preimage[other] >= 0)
              size = push(size, preimage[other]);
    }
    for (int index = 0; index < commonCount; index++) preimage[image[commonBonds[index]]] = -1;
  }

  /** Adds {@code bond} to the queue of {@code size} bonds and gives the new size. */
  private int push(int size, int bond) {
    if (size == queue.length)
      queue = Arrays.copyOf(queue, 2 * size);
    queue[size] = bond;
    return size + 1;
  }

  /** Whether selected query bond {@code bond} and its image lie on whole rings wherever they lie on rings. */
  private boolean onWholeRings(int bond) {
    return (!queryRingBond[bond] || anyWhole(queryBondCycles[bond], queryMissing))
        && (!targetRingBond[image[bond]] || anyWhole(targetBondCycles[image[bond]], targetMissing));
  }

  /** Whether any of {@code cycles} misses no bond. */
  private static boolean anyWhole(int[] cycles, int[] missing) {
    for (int cycle : cycles) if (missing[cycle] == 0) return true;
    return false;
  }

  /**
   * Whether the tetrahedral centre at mapped query atom {@code atom} and its image disagree on
   * the ligands joined to it by selected bonds, as {@link Substructure#stereoMismatchAtom} reads
   * them.
   */
  private boolean centreDisagrees(int atom, boolean[] selection) {
    return centreDisagrees(atom, mapping, selection);
  }

  /**
   * Whether the centre at {@code atom} and its image under {@code mapping} disagree on the
   * ligands joined to it by bonds of {@code selection}, or by common bonds of the mapping when
   * {@code selection} is null.
   */
  private boolean centreDisagrees(int atom, int[] mapping, boolean[] selection) {
    int mapped = mapping[atom];
    if (queryGraph.tetraChirality[atom] == 0 || targetGraph.tetraChirality[mapped] == 0)
      return false;
    int[] queryOrder = queryGraph.tetrahedralLigands(atom);
    int[] targetOrder = targetGraph.tetrahedralLigands(mapped);
    if (queryOrder == null || targetOrder == null || queryOrder.length != 4 || targetOrder.length != 4)
      return false;
    int[] permutation = {-1, -1, -1, -1};
    boolean[] used = new boolean[4];
    int mappedCount = 0;
    for (int slot = 0; slot < 4; slot++) {
      int ligand = queryOrder[slot];
      if (ligand < 0 || mapping[ligand] < 0
          || (selection != null ? !isSelected(atom, ligand, selection) : commonBond(atom, ligand, mapping) < 0))
        continue;
      int mappedLigand = mapping[ligand];
      for (int targetSlot = 0; targetSlot < 4; targetSlot++) {
        if (targetOrder[targetSlot] == mappedLigand) {
          permutation[slot] = targetSlot;
          used[targetSlot] = true;
          mappedCount++;
          break;
        }
      }
      if (permutation[slot] < 0)
        return true;
    }
    if (mappedCount < 3)
      return false;
    for (int slot = 0; slot < 4; slot++) {
      if (permutation[slot] >= 0)
        continue;
      for (int targetSlot = 0; targetSlot < 4; targetSlot++) {
        if (!used[targetSlot]) {
          permutation[slot] = targetSlot;
          used[targetSlot] = true;
          break;
        }
      }
    }
    int inversions = 0;
    for (int first = 0; first < 4; first++)
      for (int second = first + 1; second < 4; second++)
        if (permutation[first] > permutation[second])
          inversions++;
    int queryTag = queryGraph.tetraChirality[atom], targetTag = targetGraph.tetraChirality[mapped];
    return (queryTag == targetTag) != ((inversions & 1) == 0);
  }

  /**
   * Whether the double bond between mapped query atoms {@code first} and {@code second} and its
   * image disagree on the ligands joined to its ends by selected bonds, as
   * {@link Substructure#stereoMismatchAtom} reads them.
   */
  private boolean doubleBondDisagrees(int first, int second, boolean[] selection) {
    int mappedFirst = mapping[first], mappedSecond = mapping[second];
    if (mappedFirst < 0 || mappedSecond < 0)
      return false;
    int queryConfiguration = queryGraph.dbStereo(first, second);
    int targetConfiguration = targetGraph.dbStereo(mappedFirst, mappedSecond);
    if (queryConfiguration == 0 || targetConfiguration == 0)
      return false;
    for (int end = 0; end < 2; end++) {
      int queryEnd = end == 0 ? first : second, otherQueryEnd = end == 0 ? second : first;
      int targetEnd = mapping[queryEnd], otherTargetEnd = mapping[otherQueryEnd];
      int queryReference = queryGraph.doubleBondReference(queryEnd, otherQueryEnd);
      int targetReference = targetGraph.doubleBondReference(targetEnd, otherTargetEnd);
      int queryLigand = -1;
      for (int neighbor : queryGraph.neighbors[queryEnd])
        if (neighbor != otherQueryEnd && mapping[neighbor] >= 0
            && isSelected(queryEnd, neighbor, selection)) {
          queryLigand = neighbor;
          break;
        }
      if (queryLigand < 0 || queryReference < 0 || targetReference < 0)
        return false;
      if (!targetGraph.hasBond(targetEnd, mapping[queryLigand]))
        return true;
      if ((queryLigand != queryReference) != (mapping[queryLigand] != targetReference))
        queryConfiguration = 3 - queryConfiguration;
    }
    return queryConfiguration != targetConfiguration;
  }

  /** Whether the query bond between {@code atom} and {@code other} is selected. */
  private boolean isSelected(int atom, int other, boolean[] selection) {
    int bond = query.bondBetween(atom, other);
    return bond >= 0 && selection[bond];
  }

  /** Whether some atom of {@code graph} has a tetrahedral tag. */
  private static boolean hasCentres(MolGraph graph) {
    if (graph.tetraChirality == null)
      return false;
    for (int tag : graph.tetraChirality) if (tag != 0) return true;
    return false;
  }

  /** Whether some bond of {@code graph} has a double-bond configuration. */
  private static boolean hasStereoBonds(MolGraph graph) {
    for (int atom = 0; atom < graph.n; atom++)
      for (int other : graph.neighbors[atom])
        if (other > atom && graph.dbStereo(atom, other) != 0)
          return true;
    return false;
  }

  /** Whether some atom or bond of {@code graph} is on a ring. */
  private static boolean hasRings(MolGraph graph) {
    for (int atom = 0; atom < graph.n; atom++) {
      if (graph.ring[atom])
        return true;
      for (int other : graph.neighbors[atom])
        if (graph.bondInRing(atom, other))
          return true;
    }
    return graph.computeRelevantCycles().length > 0;
  }

  /** The relevant rings of {@code graph} as bonds of {@code side}. */
  private static int[][] cycles(CommonBondGraphs.Side side, MolGraph graph) {
    int[][] rings = graph.computeRelevantCycles();
    List<int[]> result = new ArrayList<>();
    for (int[] ring : rings) {
      int[] bonds = new int[ring.length];
      boolean whole = true;
      for (int index = 0; index < ring.length && whole; index++) {
        bonds[index] = side.bondBetween(ring[index], ring[(index + 1) % ring.length]);
        whole = bonds[index] >= 0;
      }
      if (whole)
        result.add(bonds);
    }
    return result.toArray(new int[0][]);
  }

  /**
   * The rings through each item, in ascending order: through each bond when {@code side} is
   * null, otherwise through each atom of {@code side}.
   */
  private static int[][] through(int[][] cycles, int items, CommonBondGraphs.Side side) {
    List<List<Integer>> lists = new ArrayList<>();
    for (int item = 0; item < items; item++) lists.add(new ArrayList<>());
    for (int cycle = 0; cycle < cycles.length; cycle++)
      for (int bond : cycles[cycle]) {
        if (side == null) {
          lists.get(bond).add(cycle);
          continue;
        }
        // Each atom of a ring is the first or second atom of two of its bonds; list the ring once.
        for (int atom : new int[] {side.bondFirst[bond], side.bondSecond[bond]}) {
          List<Integer> list = lists.get(atom);
          if (list.isEmpty() || list.get(list.size() - 1).intValue() != cycle)
            list.add(cycle);
        }
      }
    int[][] result = new int[items][];
    for (int item = 0; item < items; item++)
      result[item] = lists.get(item).stream().mapToInt(Integer::intValue).toArray();
    return result;
  }

  /** Whether each bond of {@code side} is on a ring, as {@link MCS#validateMapping} reads it. */
  private static boolean[] ringBonds(CommonBondGraphs.Side side, MolGraph graph, int[][] bondCycles) {
    boolean[] result = new boolean[side.bondCount];
    for (int bond = 0; bond < side.bondCount; bond++)
      result[bond] = bondCycles[bond].length > 0
          || graph.bondInRing(side.bondFirst[bond], side.bondSecond[bond]);
    return result;
  }

  /** Whether each atom of {@code graph} is on a ring, as {@link MCS#validateMapping} reads it. */
  private static boolean[] ringAtoms(MolGraph graph, int[][] atomCycles) {
    boolean[] result = new boolean[graph.n];
    for (int atom = 0; atom < graph.n; atom++)
      result[atom] = graph.ring[atom] || atomCycles[atom].length > 0;
    return result;
  }
}

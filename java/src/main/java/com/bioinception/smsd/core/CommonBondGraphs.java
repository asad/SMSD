/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import java.util.Set;
import java.util.function.BooleanSupplier;

/**
 * The query and target of a {@link CommonBondSearch} as int graphs, with the atom and bond
 * pairs that may be mapped onto each other. Each molecule is a {@link Side}: its neighbour
 * lists in arrays, and a kind and a class for each atom and each bond. Atom pairs are decided
 * by the atom predicate of the substructure search, and bond pairs by
 * {@link MolGraph.ChemOps#bondsCompatible}, both under the same {@link ChemOptions}. Atoms, or
 * bonds, of one molecule share a kind when they match the same items of the other, so whether
 * two items match depends on their kinds alone. A class is a connected part of the relation
 * between the kinds, so a compatible pair always shares a class. Where the options make the
 * relation an equivalence, each class holds one kind of each molecule. Otherwise, as with
 * flexible aromaticity, isotope labels left unset or ring counts compared only between ring
 * atoms, a class may hold kinds that do not all match; the items a mapping can pair within it
 * are then bounded by a maximum b-matching over the kinds that do ({@link Groups}). Every bound
 * is a true upper bound either way. The graphs also give colour refinement for the search
 * order. A {@link Clock} counts the work of every loop and throws {@link Stop} when the search
 * is to end. The molecules are only read.
 */
final class CommonBondGraphs {

  /** The number of radii, from 0, at which atom signatures are compared. */
  static final int RADII = 4;
  /** The most atom pairs a search takes, query atoms times target atoms, as a pair is packed in an int. */
  static final long MAX_PAIRS = 1L << 30;
  /**
   * The most longs the compatible bond pairs may take, query bonds times the words of the target
   * bonds: 32 MiB, and about as much again while the bond kinds are worked out. Larger graphs are
   * declined, and their MCES search finds nothing. The native core has the same limit, so both
   * decline the same graphs.
   */
  static final long MAX_BOND_PAIR_WORDS = 1L << 22;
  /** The odd constant that mixes the signature of an atom with those of its neighbours. */
  private static final long GOLDEN_GAMMA = 0x9e3779b97f4a7c15L;

  /** The graphs of the query and of the target. */
  final Side query, target;
  /** The number of atom classes and of bond classes, each shared by both molecules. */
  final int atomClasses, bondClasses;
  /**
   * The number of bond labels shared by both molecules. The label of a bond is the classes of
   * its ends and its own class, so a common bond and its image share a label.
   */
  final int labels;
  /** The lower and the higher atom class at the ends of each label. */
  final int[] labelLow, labelHigh;
  /** The atom kinds of both molecules, grouped by atom class. */
  final Groups atomGroups;
  /**
   * The fine labels of both molecules, the kinds of the ends and of the bond, grouped by label.
   * A query and a target fine label match when their bond kinds do and their ends do either
   * way round.
   */
  final Groups labelGroups;
  /** The number of longs in a row of target bonds. */
  private final int bondWords;
  /** The compatible atom pairs: a row of target atoms for each query atom. */
  private final long[] atomPairs;
  /** The compatible bond pairs: a row of {@link #bondWords} longs for each query bond. */
  private final long[] bondPairs;

  /**
   * Prepares the search of {@code query} in {@code target} under {@code chemistry}. The atoms
   * of the target in {@code excludedTargetAtoms} match nothing. Use {@link #supports} first.
   */
  CommonBondGraphs(MolGraph query, MolGraph target, ChemOptions chemistry,
      Set<Integer> excludedTargetAtoms, Clock clock) {
    if (chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    this.query = new Side(query);
    this.target = new Side(target);
    bondWords = words(this.target.bondCount);
    int atomWords = words(this.target.atomCount);
    boolean[] excluded = new boolean[this.target.atomCount];
    for (int atom : excludedTargetAtoms)
      if (atom >= 0 && atom < excluded.length)
        excluded[atom] = true;
    clock.poll();
    boolean pairRings = chemistry.completeRingsOnly && !query.ringsTruncated() && !target.ringsTruncated();
    VertexSearch.RingIndex queryRings = pairRings ? new VertexSearch.RingIndex(query) : null;
    VertexSearch.RingIndex targetRings = pairRings ? new VertexSearch.RingIndex(target) : null;
    atomPairs = new long[this.query.atomCount * atomWords];
    for (int queryAtom = 0; queryAtom < this.query.atomCount; queryAtom++) {
      clock.tick(1 + this.target.atomCount);
      for (int targetAtom = 0; targetAtom < this.target.atomCount; targetAtom++)
        if (!excluded[targetAtom]
            && (queryRings == null
                || queryRings.ringAtom(query, queryAtom) == targetRings.ringAtom(target, targetAtom))
            && MolGraph.ChemOps.atomsCompatFast(
                query, queryAtom, target, targetAtom, chemistry))
          atomPairs[queryAtom * atomWords + (targetAtom >>> 6)] |= 1L << targetAtom;
    }
    bondPairs = bondPairs(query, target, chemistry, queryRings, targetRings, clock);
    Kinds atomKinds = new Kinds(atomPairs, this.query.atomCount, this.target.atomCount, clock);
    Kinds bondKinds = new Kinds(bondPairs, this.query.bondCount, this.target.bondCount, clock);
    atomClasses = atomKinds.classes;
    bondClasses = bondKinds.classes;
    atomGroups = atomKinds.groups();
    this.query.setKinds(atomKinds.queryKind, atomKinds.queryClass, bondKinds.queryKind,
        bondKinds.queryClass);
    this.target.setKinds(atomKinds.targetKind, atomKinds.targetClass, bondKinds.targetKind,
        bondKinds.targetClass);
    // The labels, of classes and of kinds, of every bond whose ends and itself have classes.
    long[] classKeys = new long[this.query.bondCount + this.target.bondCount];
    this.query.labelKeys(atomClasses, bondClasses, classKeys, 0);
    this.target.labelKeys(atomClasses, bondClasses, classKeys, this.query.bondCount);
    labels = shared(classKeys, this.query.bondCount);
    labelLow = new int[labels];
    labelHigh = new int[labels];
    for (int bond = 0; bond < classKeys.length; bond++) {
      Side side = bond < this.query.bondCount ? this.query : this.target;
      int sideBond = bond < this.query.bondCount ? bond : bond - this.query.bondCount;
      int label = (int) classKeys[bond];
      side.label[sideBond] = label;
      if (label >= 0) {
        int first = side.atomClass[side.bondFirst[sideBond]];
        int second = side.atomClass[side.bondSecond[sideBond]];
        labelLow[label] = Math.min(first, second);
        labelHigh[label] = Math.max(first, second);
      }
    }
    int queryFine = this.query.fineLabels(atomKinds.queryKinds, bondKinds.queryKinds);
    int targetFine = this.target.fineLabels(atomKinds.targetKinds, bondKinds.targetKinds);
    labelGroups = fineLabelGroups(queryFine, targetFine, atomKinds, bondKinds);
  }

  /**
   * Whether the search can take the two graphs: neither has a bond from an atom to itself or a
   * second bond between the same atoms, every bond is listed at both its atoms, the atom pairs
   * fit an int and the bond pairs fit {@link #MAX_BOND_PAIR_WORDS}.
   */
  static boolean supports(MolGraph query, MolGraph target) {
    if ((long) query.n * target.n > MAX_PAIRS)
      return false;
    if (!wellFormed(query) || !wellFormed(target))
      return false;
    long queryBonds = query.bondCount(), targetWords = (target.bondCount() + 63) >>> 6;
    // Divided rather than multiplied, so that no product of bond counts can overflow.
    return targetWords == 0 || queryBonds <= MAX_BOND_PAIR_WORDS / targetWords;
  }

  /** The number of longs that hold {@code bits} bits. */
  static int words(int bits) {
    return (bits + 63) >>> 6;
  }

  /**
   * Packs a number of atoms and a number of bonds into a score, {@code bonds << 32 | atoms},
   * which compares by bonds first. Both counts are non-negative, so scores compare as longs.
   */
  static long score(int atoms, int bonds) {
    return (long) bonds << 32 | atoms;
  }

  /** The atoms of a score. */
  static int atoms(long score) {
    return (int) score;
  }

  /** The bonds of a score. */
  static int bonds(long score) {
    return (int) (score >>> 32);
  }

  /** The compatible atom pairs, a row of target atoms for each query atom; the caller may change them. */
  long[] compatiblePairs() {
    return atomPairs.clone();
  }

  /** Whether query bond {@code queryBond} may be mapped onto target bond {@code targetBond}. */
  boolean bondsMatch(int queryBond, int targetBond) {
    return (bondPairs[queryBond * bondWords + (targetBond >>> 6)] & 1L << targetBond) != 0;
  }

  /**
   * The compatible bond pairs. The bond predicate reads the order, aromaticity and ring status
   * of the two bonds, whether their atoms are aromatic and, with tautomers, whether both are in
   * tautomer classes, so it is asked once for each pair of bond kinds that the molecules have,
   * on a bond of each kind. A pair must match either way round, as the rows do not say which end
   * maps onto which.
   */
  private long[] bondPairs(MolGraph queryGraph, MolGraph targetGraph, ChemOptions chemistry,
      VertexSearch.RingIndex queryRings, VertexSearch.RingIndex targetRings, Clock clock) {
    long[] keys = new long[query.bondCount + target.bondCount];
    query.bondKeys(queryGraph, keys, 0, chemistry.tautomerAware);
    target.bondKeys(targetGraph, keys, query.bondCount, chemistry.tautomerAware);
    if (queryRings != null) {
      for (int bond = 0; bond < query.bondCount; bond++)
        keys[bond] = keys[bond] * 2
            + (queryRings.ringBond(queryGraph, query.bondFirst[bond], query.bondSecond[bond]) ? 1 : 0);
      for (int bond = 0; bond < target.bondCount; bond++)
        keys[query.bondCount + bond] = keys[query.bondCount + bond] * 2
            + (targetRings.ringBond(targetGraph, target.bondFirst[bond], target.bondSecond[bond]) ? 1 : 0);
    }
    int[] kind = denseRanks(keys);
    int kinds = 0;
    for (int bondKind : kind) kinds = Math.max(kinds, bondKind + 1);
    int[] queryExample = new int[kinds], targetExample = new int[kinds];
    Arrays.fill(queryExample, -1);
    Arrays.fill(targetExample, -1);
    for (int bond = query.bondCount - 1; bond >= 0; bond--) queryExample[kind[bond]] = bond;
    for (int bond = target.bondCount - 1; bond >= 0; bond--)
      targetExample[kind[query.bondCount + bond]] = bond;
    long[] byKind = new long[kinds * bondWords];
    for (int targetBond = 0; targetBond < target.bondCount; targetBond++)
      byKind[kind[query.bondCount + targetBond] * bondWords + (targetBond >>> 6)] |=
          1L << targetBond;
    boolean[] compatible = new boolean[kinds * kinds];
    for (int queryKind = 0; queryKind < kinds; queryKind++) {
      int queryBond = queryExample[queryKind];
      for (int targetKind = 0; queryBond >= 0 && targetKind < kinds; targetKind++) {
        int targetBond = targetExample[targetKind];
        clock.tick();
        if (targetBond < 0)
          continue;
        int queryFirst = query.bondFirst[queryBond], querySecond = query.bondSecond[queryBond];
        int targetFirst = target.bondFirst[targetBond], targetSecond = target.bondSecond[targetBond];
        compatible[queryKind * kinds + targetKind] =
            (queryRings == null || queryRings.ringBond(queryGraph, queryFirst, querySecond)
                == targetRings.ringBond(targetGraph, targetFirst, targetSecond))
            && MolGraph.ChemOps.bondsCompatible(queryGraph, queryFirst, querySecond, targetGraph,
                targetFirst, targetSecond, chemistry)
            && MolGraph.ChemOps.bondsCompatible(queryGraph, queryFirst, querySecond, targetGraph,
                targetSecond, targetFirst, chemistry);
      }
    }
    long[] rows = new long[query.bondCount * bondWords];
    for (int queryBond = 0; queryBond < query.bondCount; queryBond++) {
      int queryKind = kind[queryBond];
      for (int targetKind = 0; targetKind < kinds; targetKind++) {
        clock.tick(1 + bondWords);
        if (!compatible[queryKind * kinds + targetKind])
          continue;
        for (int word = 0; word < bondWords; word++)
          rows[queryBond * bondWords + word] |= byKind[targetKind * bondWords + word];
      }
    }
    return rows;
  }

  /**
   * Replaces each key that is not negative by the place among the distinct keys of both
   * molecules that both have, the query keys first in {@code keys}, or by -1 if only one has
   * it. Returns the number of shared keys.
   */
  private static int shared(long[] keys, int queryItems) {
    int[] rank = denseRanks(keys);
    int ranks = 0;
    for (int itemRank : rank) ranks = Math.max(ranks, itemRank + 1);
    boolean[] inQuery = new boolean[ranks], inTarget = new boolean[ranks];
    for (int item = 0; item < keys.length; item++)
      if (keys[item] >= 0)
        (item < queryItems ? inQuery : inTarget)[rank[item]] = true;
    int[] place = new int[ranks];
    int count = 0;
    for (int keyRank = 0; keyRank < ranks; keyRank++)
      place[keyRank] = inQuery[keyRank] && inTarget[keyRank] ? count++ : -1;
    for (int item = 0; item < keys.length; item++)
      keys[item] = keys[item] < 0 ? -1 : place[rank[item]];
    return count;
  }

  /** The fine labels of both molecules grouped by label, with the pairs of them that match. */
  private Groups fineLabelGroups(int queryFine, int targetFine, Kinds atomKinds, Kinds bondKinds) {
    int[] queryGroup = new int[queryFine], targetGroup = new int[targetFine];
    int[][] queryEnds = new int[queryFine][], targetEnds = new int[targetFine][];
    query.describeFineLabels(queryGroup, queryEnds);
    target.describeFineLabels(targetGroup, targetEnds);
    return new Groups(labels, queryGroup, targetGroup, (queryLabel, targetLabel) -> {
      int[] q = queryEnds[queryLabel], t = targetEnds[targetLabel];
      return bondKinds.matches(q[2], t[2])
          && (atomKinds.matches(q[0], t[0]) && atomKinds.matches(q[1], t[1])
              || atomKinds.matches(q[0], t[1]) && atomKinds.matches(q[1], t[0]));
    });
  }

  /**
   * The label bound, as a {@link #score}. No mapping has more atoms than the most atoms of
   * matching kinds that can be paired off within each atom class, nor more common bonds than
   * the most bonds of matching fine labels that can be paired off within each label. A
   * {@code connected} mapping is also capped by the largest component of either molecule
   * without the bonds whose label the other lacks.
   */
  long labelBound(boolean connected) {
    long mostAtoms = atomGroups.bound(query.counts(query.atomKind, atomGroups.queryGroup.length),
        target.counts(target.atomKind, atomGroups.targetGroup.length));
    long mostBonds = labelGroups.bound(query.counts(query.fineLabel, labelGroups.queryGroup.length),
        target.counts(target.fineLabel, labelGroups.targetGroup.length));
    if (connected) {
      long queryPart = query.componentBound(), targetPart = target.componentBound();
      mostAtoms = Math.min(mostAtoms, Math.min(atoms(queryPart), atoms(targetPart)));
      mostBonds = Math.min(mostBonds, Math.min(bonds(queryPart), bonds(targetPart)));
    }
    return score((int) mostAtoms, (int) mostBonds);
  }

  /** Decides whether a fine query class and a fine target class match. */
  interface Match {
    /** Whether fine query class {@code queryClass} matches fine target class {@code targetClass}. */
    boolean matches(int queryClass, int targetClass);
  }

  /**
   * The kinds of the items, atoms or bonds, of both molecules under a pair relation given as
   * rows of target items for each query item. Query items share a kind when their rows are
   * equal and target items when their columns are, so two items match exactly when their kinds
   * do. The classes are the connected parts of the relation between the kinds; a kind that
   * matches nothing has none, -1. Whether two kinds match is read from the rows, at an item of
   * each kind, so the kinds add no table of their own.
   */
  private static final class Kinds {

    /** The kind of each query item and of each target item. */
    final int[] queryKind, targetKind;
    /** The numbers of query and target kinds, and of classes. */
    final int queryKinds, targetKinds, classes;
    /** The class of each query kind and of each target kind, or -1. */
    final int[] queryClass, targetClass;
    /** The pair relation, a row of target items for each query item, and the longs of a row. */
    private final long[] rows;
    private final int words;
    /** The first item of each query and of each target kind. */
    private final int[] queryExample, targetExample;

    /** The kinds of {@code queryItems} and {@code targetItems} items under the relation of {@code rows}. */
    Kinds(long[] rows, int queryItems, int targetItems, Clock clock) {
      this.rows = rows;
      words = words(targetItems);
      queryKind = equalRows(rows, queryItems, words, clock);
      targetKind = equalRows(transpose(rows, queryItems, targetItems, clock), targetItems,
          words(queryItems), clock);
      queryKinds = count(queryKind);
      targetKinds = count(targetKind);
      queryExample = examples(queryKind, queryKinds);
      targetExample = examples(targetKind, targetKinds);
      int[] sets = UnionFind.create(queryKinds + targetKinds);
      for (int query = 0; query < queryKinds; query++) {
        clock.tick(1 + targetKinds);
        for (int target = 0; target < targetKinds; target++)
          if (matches(query, target))
            UnionFind.join(sets, query, queryKinds + target);
      }
      int[] members = new int[queryKinds + targetKinds], setClass = new int[queryKinds + targetKinds];
      for (int kind = 0; kind < members.length; kind++) members[UnionFind.root(sets, kind)]++;
      int count = 0;
      for (int kind = 0; kind < members.length; kind++)
        setClass[kind] = members[kind] > 1 ? count++ : -1;
      classes = count;
      queryClass = new int[queryKinds];
      targetClass = new int[targetKinds];
      for (int kind = 0; kind < queryKinds; kind++) queryClass[kind] = setClass[UnionFind.root(sets, kind)];
      for (int kind = 0; kind < targetKinds; kind++)
        targetClass[kind] = setClass[UnionFind.root(sets, queryKinds + kind)];
    }

    /** Whether query kind {@code queryKind} matches target kind {@code targetKind}. */
    boolean matches(int queryKind, int targetKind) {
      int item = targetExample[targetKind];
      return (rows[queryExample[queryKind] * words + (item >>> 6)] & 1L << item) != 0;
    }

    /** The kinds grouped by class. */
    Groups groups() {
      return new Groups(classes, queryClass, targetClass, this::matches);
    }

    /** The number of kinds, one more than the highest. */
    private static int count(int[] kinds) {
      int count = 0;
      for (int kind : kinds) count = Math.max(count, kind + 1);
      return count;
    }

    /** The first item of each kind. */
    private static int[] examples(int[] kinds, int count) {
      int[] example = new int[count];
      for (int item = kinds.length - 1; item >= 0; item--) example[kinds[item]] = item;
      return example;
    }
  }

  /**
   * Fine classes of both molecules, atom kinds or fine labels, grouped by the coarse class
   * they fall in; no two of different groups match. A mapping pairs items of matching fine
   * classes, so the items it can pair within a group are at most a maximum bipartite b-matching
   * between the query and target counts of the group over the matching pairs: the most flow
   * from the query classes, each with its count, to the target classes, each with its count,
   * along those pairs. A group in which every pair matches, as each group does when the
   * relation is an equivalence, gives the fewer items of either side. The flow is found by
   * augmenting paths on the small graph of one group.
   */
  static final class Groups {

    /** The number of groups. */
    final int count;
    /** The group of each fine query class and of each fine target class, -1 for none. */
    final int[] queryGroup, targetGroup;
    /** The fine query and target classes of each group. */
    private final int[][] queryMembers, targetMembers;
    /** Whether each query member of a group matches each target member, a row for each query member. */
    private final boolean[][] match;
    /** Whether every member pair of each group matches. */
    final boolean[] complete;
    /** Whether every group is complete. */
    final boolean allComplete;
    /** Scratch of {@link #bound(int, int[], int[])}: what is left of each count, the flow on each pair, and the search for a path. */
    private final int[] supply, demand, flow, previous, queue;

    /** Groups the fine classes by the groups given, with the pairs that {@code match} accepts. */
    Groups(int count, int[] queryGroup, int[] targetGroup, Match match) {
      this.count = count;
      this.queryGroup = queryGroup;
      this.targetGroup = targetGroup;
      queryMembers = members(queryGroup, count);
      targetMembers = members(targetGroup, count);
      this.match = new boolean[count][];
      complete = new boolean[count];
      boolean all = true;
      int mostQuery = 0, mostTarget = 0, mostPairs = 0;
      for (int group = 0; group < count; group++) {
        int[] q = queryMembers[group], t = targetMembers[group];
        boolean[] pairs = new boolean[q.length * t.length];
        boolean every = true;
        for (int i = 0; i < q.length; i++) {
          for (int j = 0; j < t.length; j++) {
            pairs[i * t.length + j] = match.matches(q[i], t[j]);
            every &= pairs[i * t.length + j];
          }
        }
        this.match[group] = pairs;
        complete[group] = every;
        all &= every;
        mostQuery = Math.max(mostQuery, q.length);
        mostTarget = Math.max(mostTarget, t.length);
        mostPairs = Math.max(mostPairs, pairs.length);
      }
      allComplete = all;
      supply = new int[mostQuery];
      demand = new int[mostTarget];
      flow = new int[mostPairs];
      previous = new int[mostQuery + mostTarget];
      queue = new int[mostQuery + mostTarget];
    }

    /** The members of each group. */
    private static int[][] members(int[] group, int count) {
      int[] size = new int[count];
      for (int item : group)
        if (item >= 0)
          size[item]++;
      int[][] members = new int[count][];
      for (int g = 0; g < count; g++) members[g] = new int[size[g]];
      Arrays.fill(size, 0);
      for (int item = 0; item < group.length; item++)
        if (group[item] >= 0)
          members[group[item]][size[group[item]]++] = item;
      return members;
    }

    /** The most items paired off over all groups, from the counts of each fine query and target class. */
    long bound(int[] queryCounts, int[] targetCounts) {
      long total = 0;
      for (int group = 0; group < count; group++) total += bound(group, queryCounts, targetCounts);
      return total;
    }

    /**
     * The most items of {@code group} paired off, each query item with a target item of a
     * matching class, from the counts of each fine query and target class.
     */
    int bound(int group, int[] queryCounts, int[] targetCounts) {
      int[] q = queryMembers[group], t = targetMembers[group];
      int queryCount = q.length, targetCount = t.length, total = 0;
      if (complete[group]) {
        int queryItems = 0, targetItems = 0;
        for (int member : q) queryItems += queryCounts[member];
        for (int member : t) targetItems += targetCounts[member];
        return Math.min(queryItems, targetItems);
      }
      boolean[] pairs = match[group];
      for (int i = 0; i < queryCount; i++) supply[i] = queryCounts[q[i]];
      for (int j = 0; j < targetCount; j++) demand[j] = targetCounts[t[j]];
      // A greedy flow to start, then augmenting paths, which may take back flow placed before.
      for (int i = 0; i < queryCount; i++) {
        for (int j = 0; j < targetCount; j++) {
          int amount = pairs[i * targetCount + j] ? Math.min(supply[i], demand[j]) : 0;
          flow[i * targetCount + j] = amount;
          supply[i] -= amount;
          demand[j] -= amount;
          total += amount;
        }
      }
      // Query classes are nodes 0 to queryCount - 1, and target classes follow them.
      for (;;) {
        int tail = 0, end = -1;
        Arrays.fill(previous, 0, queryCount + targetCount, -2);
        for (int i = 0; i < queryCount; i++) {
          if (supply[i] > 0) {
            previous[i] = -1;
            queue[tail++] = i;
          }
        }
        for (int head = 0; head < tail && end < 0; head++) {
          int node = queue[head];
          if (node < queryCount) {
            for (int j = 0; j < targetCount && end < 0; j++) {
              if (!pairs[node * targetCount + j] || previous[queryCount + j] != -2)
                continue;
              previous[queryCount + j] = node;
              if (demand[j] > 0)
                end = queryCount + j;
              else
                queue[tail++] = queryCount + j;
            }
          } else {
            for (int i = 0; i < queryCount; i++) {
              if (flow[i * targetCount + node - queryCount] > 0 && previous[i] == -2) {
                previous[i] = node;
                queue[tail++] = i;
              }
            }
          }
        }
        if (end < 0)
          return total;
        int amount = demand[end - queryCount], node = end;
        for (; previous[node] != -1; node = previous[node])
          if (node < queryCount)
            amount = Math.min(amount, flow[node * targetCount + previous[node] - queryCount]);
        amount = Math.min(amount, supply[node]);
        supply[node] -= amount;
        demand[end - queryCount] -= amount;
        for (node = end; previous[node] != -1; node = previous[node]) {
          if (node < queryCount)
            flow[node * targetCount + previous[node] - queryCount] -= amount;
          else
            flow[previous[node] * targetCount + node - queryCount] += amount;
        }
        total += amount;
      }
    }
  }

  /**
   * Numbers the {@code count} rows of {@code words} longs so that equal rows share a number, in
   * order of first appearance. The first row of each number is kept in an open-addressed table
   * by the hash of its contents.
   */
  private static int[] equalRows(long[] rows, int count, int words, Clock clock) {
    int[] number = new int[count];
    int capacity = 2;
    while (capacity < 2L * count && capacity < 1 << 30)
      capacity <<= 1;
    int[] first = new int[capacity];
    Arrays.fill(first, -1);
    int mask = first.length - 1, numbers = 0;
    for (int row = 0; row < count; row++) {
      clock.tick(1 + words);
      int slot = rowHash(rows, row * words, words) & mask;
      while (first[slot] >= 0 && !sameRow(rows, first[slot] * words, row * words, words))
        slot = (slot + 1) & mask;
      if (first[slot] < 0) {
        first[slot] = row;
        number[row] = numbers++;
      } else {
        number[row] = number[first[slot]];
      }
    }
    return number;
  }

  /** A hash of the {@code words} longs of {@code rows} from {@code start}. */
  private static int rowHash(long[] rows, int start, int words) {
    long value = 1;
    for (int i = start; i < start + words; i++) value = 31 * value + rows[i];
    value ^= value >>> 32;
    value *= 0x9E3779B97F4A7C15L;
    return (int) (value >>> 32);
  }

  /** Whether the rows of {@code words} longs of {@code rows} from {@code first} and {@code second} are equal. */
  private static boolean sameRow(long[] rows, int first, int second, int words) {
    for (int i = 0; i < words; i++)
      if (rows[first + i] != rows[second + i])
        return false;
    return true;
  }

  /** The columns of {@code rowCount} rows of {@code columnCount} bits as rows of {@code rowCount} bits. */
  private static long[] transpose(long[] rows, int rowCount, int columnCount, Clock clock) {
    int words = words(columnCount), columnWords = words(rowCount);
    long[] columns = new long[columnCount * columnWords];
    for (int row = 0; row < rowCount; row++) {
      clock.tick(1 + words);
      for (int column = nextBit(rows, row * words, words, 0); column >= 0;
           column = nextBit(rows, row * words, words, column + 1))
        columns[column * columnWords + (row >>> 6)] |= 1L << row;
    }
    return columns;
  }

  /**
   * Colour refinement from dense colours. Cells are split by the sorted ranks of the
   * neighbours until none splits. With {@code individualise}, the atom of highest index in the
   * lowest shared cell is then put before its peers and the cells split again. This goes on
   * until each atom has a rank of its own, so that of atoms alike the one of lowest index ranks
   * highest. The search tries the higher ranked query atoms first and the target atoms in index
   * order, so a molecule matched with itself pairs alike atoms by index. A rank is where the
   * cell starts. A round sorts only
   * the neighbours of atoms that changed cell, and the largest part of a split cell keeps its
   * name, so a long chain does not cost a full sort each round. The clock can stop it.
   */
  static int[] refine(int[] start, int[] neighbour, int[] colour, boolean individualise,
      Clock clock) {
    return new Refinement(start, neighbour, colour, clock).ranks(individualise);
  }

  /** One colour refinement, see {@link #refine}. The atoms are kept in rank order, and each cell is a run of places in that order. */
  private static final class Refinement {

    /** The neighbour lists, as in {@link Side#start} and {@link Side#neighbour}. */
    private final int[] start, neighbour;
    /** Counts the work of each round. */
    private final Clock clock;
    /** The number of atoms. */
    private final int atomCount;
    /** The atoms in rank order, and the place of each atom in it. */
    private final int[] order, place;
    /** The cell of each atom, and the place where each cell starts and its size. */
    private final int[] cell, cellStart, cellSize;
    /** The number of cells. */
    private int cells;
    /** The sorted ranks of the neighbours of each atom, at its places in the neighbour list. */
    private final int[] ranks;
    /** The atoms to sort in a round, the first {@link #toSortCount}, and the last round each was listed for. */
    private final int[] toSort, listed;
    /** The number of atoms to sort. */
    private int toSortCount;
    /** The atoms that changed cell in a round, the first {@link #movedCount}. */
    private final int[] moved;
    /** The number of atoms that changed cell in the round. */
    private int movedCount;
    /** The round being run, from 1, which stamps {@link #listed}. */
    private int round = 1;
    /** Room to sort the atoms of one cell, and to merge them, see {@link #sortByNeighbours}. */
    private final int[] sorting, merging;
    /** Room for the atoms to sort as keys by cell, see {@link #sortByCell}. */
    private final long[] byCell;

    /** Starts with a cell for each colour, and every atom listed to sort. */
    Refinement(int[] start, int[] neighbour, int[] colour, Clock clock) {
      this.start = start;
      this.neighbour = neighbour;
      this.clock = clock;
      atomCount = colour.length;
      order = new int[atomCount];
      place = new int[atomCount];
      cell = new int[atomCount];
      cellStart = new int[atomCount];
      cellSize = new int[atomCount];
      ranks = new int[neighbour.length];
      toSort = new int[atomCount];
      listed = new int[atomCount];
      moved = new int[atomCount];
      sorting = new int[atomCount];
      merging = new int[atomCount];
      byCell = new long[atomCount];
      for (int atomColour : colour) cellSize[atomColour]++;
      while (cells < atomCount && cellSize[cells] > 0) {
        cellStart[cells] = cells == 0 ? 0 : cellStart[cells - 1] + cellSize[cells - 1];
        cells++;
      }
      int[] filled = new int[atomCount];
      for (int atom = 0; atom < atomCount; atom++) {
        cell[atom] = colour[atom];
        place[atom] = cellStart[cell[atom]] + filled[cell[atom]]++;
        order[place[atom]] = atom;
        toSort[atom] = atom;
      }
      toSortCount = atomCount;
      Arrays.fill(listed, round);
    }

    /** Refines the cells, and gives the rank of each atom. */
    int[] ranks(boolean individualise) {
      for (int firstShared = 0;; round++) {
        movedCount = 0;
        int shared = sortByCell();
        // Every rank is read before any cell splits.
        for (int from = 0, to; from < shared; from = to) {
          to = cellEnd(from, shared);
          readRanks(from, to);
        }
        for (int from = 0, to; from < shared; from = to) {
          to = cellEnd(from, shared);
          splitCell(from, to);
        }
        if (movedCount == 0) {
          while (firstShared < atomCount && cellSize[cell[order[firstShared]]] == 1)
            firstShared++;
          if (!individualise || firstShared == atomCount)
            break;
          individualise(firstShared);
        }
        listNeighbours();
      }
      int[] rank = new int[atomCount];
      for (int atom = 0; atom < atomCount; atom++) rank[atom] = cellStart[cell[atom]];
      return rank;
    }

    /** Keeps the atoms to sort whose cell has other atoms, in order of cell and then of atom, and gives their number. */
    private int sortByCell() {
      int shared = 0;
      for (int i = 0; i < toSortCount; i++) {
        int atom = toSort[i];
        if (cellSize[cell[atom]] > 1)
          byCell[shared++] = (long) cell[atom] << 32 | atom;
      }
      Arrays.sort(byCell, 0, shared);
      for (int i = 0; i < shared; i++) toSort[i] = (int) byCell[i];
      return shared;
    }

    /** The end of the atoms to sort from place {@code from} that are in its cell. */
    private int cellEnd(int from, int shared) {
      int to = from + 1;
      while (to < shared && cell[toSort[to]] == cell[toSort[from]]) to++;
      return to;
    }

    /**
     * Sorts the neighbour ranks of the atoms to sort from {@code from} to {@code to - 1}, all
     * of one cell, and of one other atom of the cell if there is one. That atom stands for the
     * rest, which have the same ranks.
     */
    private void readRanks(int from, int to) {
      for (int i = from; i < to; i++) sortNeighbourRanks(toSort[i]);
      int standIn = unlisted(cell[toSort[from]]);
      if (standIn >= 0)
        sortNeighbourRanks(standIn);
    }

    /** Puts the sorted ranks of the neighbours of {@code atom} at its places in {@link #ranks}. */
    private void sortNeighbourRanks(int atom) {
      clock.tick(1 + start[atom + 1] - start[atom]);
      for (int neighbourPlace = start[atom]; neighbourPlace < start[atom + 1]; neighbourPlace++)
        ranks[neighbourPlace] = cellStart[cell[neighbour[neighbourPlace]]];
      Arrays.sort(ranks, start[atom], start[atom + 1]);
    }

    /** An atom of cell {@code cellName} that is not listed for this round, or -1 if there is none. */
    private int unlisted(int cellName) {
      for (int position = cellStart[cellName]; position < cellStart[cellName] + cellSize[cellName];
           position++)
        if (listed[order[position]] != round)
          return order[position];
      return -1;
    }

    /**
     * Splits the cell of the atoms to sort from {@code from} to {@code to - 1} into runs of
     * equal neighbour ranks. The unlisted atoms of the cell form one run, that of their
     * stand-in. Each run becomes a cell, and the largest keeps the name of the cell.
     */
    private void splitCell(int from, int to) {
      int cellName = cell[toSort[from]], standIn = unlisted(cellName);
      int unlistedCount = cellSize[cellName] - (to - from);
      int count = 0;
      for (int i = from; i < to; i++) sorting[count++] = toSort[i];
      if (standIn >= 0)
        sorting[count++] = standIn;
      sortByNeighbours(count);
      clock.tick(count);
      if (runEnd(0, count) == count)
        return;
      int largestRun = 0, largestSize = 0, standInRun = -1;
      for (int runFrom = 0, runTo; runFrom < count; runFrom = runTo) {
        runTo = runEnd(runFrom, count);
        int size = runSize(runFrom, runTo, standIn, unlistedCount);
        if (size > largestSize) {
          largestRun = runFrom;
          largestSize = size;
        }
        for (int sortedPlace = runFrom; sortedPlace < runTo; sortedPlace++)
          if (sorting[sortedPlace] == standIn)
            standInRun = runFrom;
      }
      // The listed atoms move to their places in run order, around the run of the stand-in.
      // That run then holds the unlisted atoms.
      int renamed = -1;
      for (int runFrom = 0, runTo, position = cellStart[cellName]; runFrom < count;
           runFrom = runTo) {
        runTo = runEnd(runFrom, count);
        int size = runSize(runFrom, runTo, standIn, unlistedCount);
        int runCell = runFrom == largestRun ? cellName : cells++;
        cellStart[runCell] = position;
        cellSize[runCell] = size;
        if (runFrom == standInRun) {
          if (runCell != cellName)
            renamed = runCell;
        } else {
          for (int sortedPlace = runFrom; sortedPlace < runTo; sortedPlace++) {
            moveTo(sorting[sortedPlace], position + sortedPlace - runFrom);
            if (runCell != cellName)
              changeCell(sorting[sortedPlace], runCell);
          }
        }
        position += size;
      }
      if (renamed >= 0)
        for (int position = cellStart[renamed]; position < cellStart[renamed] + cellSize[renamed];
             position++)
          changeCell(order[position], renamed);
    }

    /** The end of the run of sorted atoms from {@code runFrom} with equal neighbour ranks. */
    private int runEnd(int runFrom, int count) {
      int runTo = runFrom + 1;
      while (runTo < count && compareNeighbours(sorting[runFrom], sorting[runTo]) == 0) runTo++;
      return runTo;
    }

    /** Orders two atoms by their numbers of neighbours, and then by the sorted ranks of their neighbours. */
    private int compareNeighbours(int left, int right) {
      int comparison =
          Integer.compare(start[left + 1] - start[left], start[right + 1] - start[right]);
      for (int leftPlace = start[left], rightPlace = start[right];
           comparison == 0 && leftPlace < start[left + 1]; leftPlace++, rightPlace++)
        comparison = Integer.compare(ranks[leftPlace], ranks[rightPlace]);
      return comparison;
    }

    /**
     * Sorts the first {@code count} atoms of {@link #sorting} by {@link #compareNeighbours},
     * keeping atoms that compare equal in their order: runs of a few atoms by insertion, then
     * merged pairwise. A stable sort has one result, so the cells split as with any other.
     */
    private void sortByNeighbours(int count) {
      final int run = 16;
      for (int from = 0; from < count; from += run) {
        int to = Math.min(count, from + run);
        for (int next = from + 1; next < to; next++) {
          int atom = sorting[next], at = next;
          while (at > from && compareNeighbours(sorting[at - 1], atom) > 0) {
            sorting[at] = sorting[at - 1];
            at--;
          }
          sorting[at] = atom;
        }
      }
      int[] source = sorting, merged = merging;
      for (int width = run; width < count; width <<= 1) {
        for (int from = 0; from < count; from += width << 1) {
          int middle = Math.min(count, from + width), to = Math.min(count, middle + width);
          int left = from, right = middle, out = from;
          while (left < middle && right < to)
            merged[out++] = compareNeighbours(source[right], source[left]) < 0 ? source[right++] : source[left++];
          while (left < middle) merged[out++] = source[left++];
          while (right < to) merged[out++] = source[right++];
        }
        int[] swap = source;
        source = merged;
        merged = swap;
      }
      if (source != sorting)
        System.arraycopy(source, 0, sorting, 0, count);
    }

    /** The atoms in the run of sorted atoms from {@code runFrom} to {@code runTo - 1}; the stand-in counts for every unlisted atom. */
    private int runSize(int runFrom, int runTo, int standIn, int unlistedCount) {
      int size = 0;
      for (int sortedPlace = runFrom; sortedPlace < runTo; sortedPlace++)
        size += sorting[sortedPlace] == standIn ? unlistedCount : 1;
      return size;
    }

    /** Puts the atom of highest index in the shared cell that starts at place {@code first} in a cell of its own, before its peers. */
    private void individualise(int first) {
      int cellName = cell[order[first]], pick = first;
      clock.tick(cellSize[cellName]);
      for (int position = first + 1; position < first + cellSize[cellName]; position++)
        if (order[position] > order[pick])
          pick = position;
      int atom = order[pick];
      moveTo(atom, first);
      cellStart[cells] = first;
      cellSize[cells] = 1;
      cellStart[cellName]++;
      cellSize[cellName]--;
      changeCell(atom, cells++);
    }

    /** Puts {@code atom} at place {@code newPlace}, and the atom that was there at its old place. */
    private void moveTo(int atom, int newPlace) {
      int displaced = order[newPlace];
      order[place[atom]] = displaced;
      place[displaced] = place[atom];
      order[newPlace] = atom;
      place[atom] = newPlace;
    }

    /** Puts {@code atom} in cell {@code cellName}, and notes that it moved. */
    private void changeCell(int atom, int cellName) {
      cell[atom] = cellName;
      moved[movedCount++] = atom;
    }

    /** Lists the neighbours of the atoms that changed cell, to sort in the next round. */
    private void listNeighbours() {
      toSortCount = 0;
      for (int movedIndex = 0; movedIndex < movedCount; movedIndex++) {
        int atom = moved[movedIndex];
        clock.tick(1 + start[atom + 1] - start[atom]);
        for (int neighbourPlace = start[atom]; neighbourPlace < start[atom + 1]; neighbourPlace++)
          if (listed[neighbour[neighbourPlace]] != round + 1) {
            listed[neighbour[neighbourPlace]] = round + 1;
            toSort[toSortCount++] = neighbour[neighbourPlace];
          }
      }
    }
  }

  /** The rank of each key among the distinct keys, from 0. */
  static int[] denseRanks(long[] keys) {
    long[] sorted = keys.clone();
    Arrays.sort(sorted);
    int size = 0;
    for (int i = 0; i < sorted.length; i++)
      if (i == 0 || sorted[i] != sorted[i - 1])
        sorted[size++] = sorted[i];
    int[] ranks = new int[keys.length];
    for (int i = 0; i < keys.length; i++) ranks[i] = Arrays.binarySearch(sorted, 0, size, keys[i]);
    return ranks;
  }

  /** Mixes the bits of {@code value} with a 64-bit finaliser. */
  private static long mix(long value) {
    value = (value ^ value >>> 30) * 0xbf58476d1ce4e5b9L;
    value = (value ^ value >>> 27) * 0x94d049bb133111ebL;
    return value ^ value >>> 31;
  }

  /**
   * Sorts the atoms of a seed frontier from {@code from} up to {@code size} by descending
   * {@code rank} and merges them into the atoms before them, which are in that order already.
   * {@code keys} and {@code merged} are scratch as large as the frontier.
   */
  static void mergeByRank(int[] frontier, int from, int size, int[] rank, long[] keys, int[] merged,
      Clock clock) {
    if (from == size)
      return;
    clock.tick((int) Math.min(Integer.MAX_VALUE, size + 16L * (size - from)));
    // A key sorts by descending rank and carries its atom.
    for (int i = from; i < size; i++)
      keys[i] = (long) (Integer.MAX_VALUE - rank[frontier[i]]) << 32 | frontier[i];
    Arrays.sort(keys, from, size);
    for (int i = from; i < size; i++)
      frontier[i] = (int) keys[i];
    int head = 0, tail = from, out = 0;
    while (head < from && tail < size)
      merged[out++] = rank[frontier[tail]] > rank[frontier[head]] ? frontier[tail++] : frontier[head++];
    while (head < from)
      merged[out++] = frontier[head++];
    while (tail < size)
      merged[out++] = frontier[tail++];
    System.arraycopy(merged, 0, frontier, 0, size);
  }

  /** The inverse of a permutation: the position of each value. */
  static int[] inverse(int[] permutation) {
    int[] inverse = new int[permutation.length];
    for (int i = 0; i < permutation.length; i++) inverse[permutation[i]] = i;
    return inverse;
  }

  /** The first set bit of a row of {@code words} longs from {@code start}, at or after bit {@code from}, or -1. */
  static int nextBit(long[] bits, int start, int words, int from) {
    int wordIndex = from >>> 6;
    if (wordIndex >= words)
      return -1;
    long word = bits[start + wordIndex] & -1L << from;
    while (word == 0) {
      if (++wordIndex == words)
        return -1;
      word = bits[start + wordIndex];
    }
    return (wordIndex << 6) + Long.numberOfTrailingZeros(word);
  }

  /**
   * Whether every neighbour list is free of the atom itself and of repeats, and every bond is
   * listed at both its atoms. The atoms that list each atom are gathered first, so that the
   * check takes time linear in the bonds even next to a hub.
   */
  private static boolean wellFormed(MolGraph graph) {
    int[] seen = new int[graph.n], listedStart = new int[graph.n + 1];
    for (int atom = 0; atom < graph.n; atom++) {
      for (int other : graph.neighbors[atom]) {
        if (other < 0 || other >= graph.n || other == atom || seen[other] == atom + 1)
          return false;
        seen[other] = atom + 1;
        listedStart[other + 1]++;
      }
    }
    // Each atom is listed by as many atoms as it lists, and then by exactly those it lists.
    for (int atom = 0; atom < graph.n; atom++) {
      if (listedStart[atom + 1] != graph.neighbors[atom].length)
        return false;
      listedStart[atom + 1] += listedStart[atom];
    }
    int[] listedBy = new int[listedStart[graph.n]], next = Arrays.copyOf(listedStart, graph.n);
    for (int atom = 0; atom < graph.n; atom++)
      for (int other : graph.neighbors[atom]) listedBy[next[other]++] = atom;
    Arrays.fill(seen, 0);
    for (int atom = 0; atom < graph.n; atom++) {
      for (int other : graph.neighbors[atom]) seen[other] = atom + 1;
      for (int place = listedStart[atom]; place < listedStart[atom + 1]; place++)
        if (seen[listedBy[place]] != atom + 1)
          return false;
    }
    return true;
  }


  /**
   * One molecule of the search as an int graph, with the atoms numbered as in the molecule and
   * the bonds numbered from the lower atom up, in the order of its neighbours.
   */
  static final class Side {

    /** The numbers of atoms and bonds. */
    final int atomCount, bondCount;
    /** The neighbours of an atom, and the bonds to them, are at places start[atom] to start[atom + 1] - 1. */
    final int[] start, neighbour, neighbourBond;
    /** The lower and the higher atom of each bond. */
    final int[] bondFirst, bondSecond;
    /** A hash of the element of each atom, for the search order and seeds. */
    final int[] symbolHash;
    /** The kind of each atom and of each bond, see {@link Kinds}. */
    int[] atomKind, bondKind;
    /** The class of each atom and of each bond, -1 for one that matches nothing. */
    int[] atomClass, bondClass;
    /** The label of each bond, see {@link CommonBondGraphs#labels}, or -1 if the other molecule lacks it. */
    final int[] label;
    /** The fine label of each bond that has a label, see {@link CommonBondGraphs#labelGroups}, or -1. */
    final int[] fineLabel;

    /** The graph of {@code graph}, which must be well formed. */
    Side(MolGraph graph) {
      atomCount = graph.n;
      start = new int[atomCount + 1];
      for (int atom = 0; atom < atomCount; atom++)
        start[atom + 1] = start[atom] + graph.neighbors[atom].length;
      neighbour = new int[start[atomCount]];
      neighbourBond = new int[start[atomCount]];
      bondCount = start[atomCount] / 2;
      bondFirst = new int[bondCount];
      bondSecond = new int[bondCount];
      int bonds = 0;
      for (int atom = 0; atom < atomCount; atom++) {
        int[] row = graph.neighbors[atom];
        for (int i = 0; i < row.length; i++) {
          neighbour[start[atom] + i] = row[i];
          if (atom < row[i]) {
            bondFirst[bonds] = atom;
            bondSecond[bonds] = row[i];
            neighbourBond[start[atom] + i] = bonds++;
          }
        }
      }
      // A place from an atom to a lower neighbour takes the bond that the neighbour numbered.
      // The bonds are gathered by their higher atom, so that this takes time linear in the
      // bonds even next to a hub.
      int[] higherStart = new int[atomCount + 1], byHigher = new int[bondCount], bondTo = new int[atomCount];
      for (int bond = 0; bond < bondCount; bond++) higherStart[bondSecond[bond] + 1]++;
      for (int atom = 0; atom < atomCount; atom++) higherStart[atom + 1] += higherStart[atom];
      int[] next = Arrays.copyOf(higherStart, atomCount);
      for (int bond = 0; bond < bondCount; bond++) byHigher[next[bondSecond[bond]]++] = bond;
      for (int atom = 0; atom < atomCount; atom++) {
        for (int i = higherStart[atom]; i < higherStart[atom + 1]; i++)
          bondTo[bondFirst[byHigher[i]]] = byHigher[i];
        for (int place = start[atom]; place < start[atom + 1]; place++)
          if (neighbour[place] < atom)
            neighbourBond[place] = bondTo[neighbour[place]];
      }
      symbolHash = new int[atomCount];
      for (int atom = 0; atom < atomCount; atom++) symbolHash[atom] = 1 + graph.atomicNum[atom];
      label = new int[bondCount];
      fineLabel = new int[bondCount];
    }

    /** Sets the kinds and classes of the atoms and bonds from those of their kinds. */
    void setKinds(int[] atomKind, int[] atomKindClass, int[] bondKind, int[] bondKindClass) {
      this.atomKind = atomKind;
      this.bondKind = bondKind;
      atomClass = new int[atomCount];
      bondClass = new int[bondCount];
      for (int atom = 0; atom < atomCount; atom++) atomClass[atom] = atomKindClass[atomKind[atom]];
      for (int bond = 0; bond < bondCount; bond++) bondClass[bond] = bondKindClass[bondKind[bond]];
    }

    /**
     * Puts the label key of each bond into {@code keys} from {@code offset}: the classes of its
     * ends, the lower first, and its own class, or -1 if any of them is missing.
     */
    void labelKeys(int atomClasses, int bondClasses, long[] keys, int offset) {
      for (int bond = 0; bond < bondCount; bond++) {
        int first = atomClass[bondFirst[bond]], second = atomClass[bondSecond[bond]];
        keys[offset + bond] = first < 0 || second < 0 || bondClass[bond] < 0 ? -1
            : ((long) Math.min(first, second) * atomClasses + Math.max(first, second)) * bondClasses
                + bondClass[bond];
      }
    }

    /**
     * Numbers the fine labels of the bonds that have a label: the kinds of the ends, the lower
     * first, and the kind of the bond. Returns their number.
     */
    int fineLabels(int atomKinds, int bondKinds) {
      long[] keys = new long[bondCount];
      for (int bond = 0; bond < bondCount; bond++) {
        int first = atomKind[bondFirst[bond]], second = atomKind[bondSecond[bond]];
        keys[bond] = label[bond] < 0 ? -1
            : ((long) Math.min(first, second) * atomKinds + Math.max(first, second)) * bondKinds
                + bondKind[bond];
      }
      long[] sorted = keys.clone();
      Arrays.sort(sorted);
      int size = 0;
      for (int i = 0; i < sorted.length; i++)
        if (sorted[i] >= 0 && (size == 0 || sorted[i] != sorted[size - 1]))
          sorted[size++] = sorted[i];
      for (int bond = 0; bond < bondCount; bond++)
        fineLabel[bond] = keys[bond] < 0 ? -1 : Arrays.binarySearch(sorted, 0, size, keys[bond]);
      return size;
    }

    /**
     * Sets, for each fine label, its label in {@code group} and in {@code ends} the kinds of its
     * ends, the lower first, and of its bond.
     */
    void describeFineLabels(int[] group, int[][] ends) {
      for (int bond = 0; bond < bondCount; bond++) {
        int fine = fineLabel[bond];
        if (fine < 0 || ends[fine] != null)
          continue;
        int first = atomKind[bondFirst[bond]], second = atomKind[bondSecond[bond]];
        group[fine] = label[bond];
        ends[fine] = new int[] {Math.min(first, second), Math.max(first, second), bondKind[bond]};
      }
    }

    /** The number of items of each fine class, from the fine class of each item, -1 for none. */
    int[] counts(int[] fine, int classes) {
      int[] counts = new int[classes];
      for (int item : fine)
        if (item >= 0)
          counts[item]++;
      return counts;
    }

    /**
     * Puts the key of each bond for the bond predicate, its order, aromaticity and ring
     * status, whether both its atoms are aromatic and, with {@code tautomers}, whether both are
     * in tautomer classes, into {@code keys} from {@code offset}. Keep it in step with
     * {@link MolGraph.ChemOps#bondsCompatible}.
     */
    void bondKeys(MolGraph graph, long[] keys, int offset, boolean tautomers) {
      for (int bond = 0; bond < bondCount; bond++) {
        int first = bondFirst[bond], second = bondSecond[bond];
        boolean tautomeric = tautomers && graph.tautomerClass != null
            && graph.tautomerClass[first] != -1 && graph.tautomerClass[second] != -1;
        keys[offset + bond] = (long) graph.bondOrder(first, second) << 4
            | (tautomeric ? 8 : 0)
            | (graph.bondAromatic(first, second) ? 4 : 0)
            | (graph.bondInRing(first, second) ? 2 : 0)
            | (graph.aromatic[first] && graph.aromatic[second] ? 1 : 0);
      }
    }

    /** The number of bonds at {@code atom}. */
    int degree(int atom) {
      return start[atom + 1] - start[atom];
    }

    /** The most bonds at any atom. */
    int maxDegree() {
      int max = 0;
      for (int atom = 0; atom < atomCount; atom++) max = Math.max(max, degree(atom));
      return max;
    }

    /**
     * The bond between atoms {@code atom} and {@code other}, or -1 if there is none. The shorter
     * list of neighbours is read.
     */
    int bondBetween(int atom, int other) {
      if (degree(atom) > degree(other)) {
        int swap = atom;
        atom = other;
        other = swap;
      }
      for (int i = start[atom]; i < start[atom + 1]; i++)
        if (neighbour[i] == other)
          return neighbourBond[i];
      return -1;
    }

    /** The dense ranks of the atoms by their number of choices, then their degree, then their element. */
    int[] invariants(int[] choices) {
      long[] keys = new long[atomCount];
      for (int atom = 0; atom < atomCount; atom++) keys[atom] = symbolHash[atom];
      int[] symbol = denseRanks(keys);
      long degrees = maxDegree() + 1;
      for (int atom = 0; atom < atomCount; atom++)
        keys[atom] = (choices[atom] * degrees + degree(atom)) * atomCount + symbol[atom];
      return denseRanks(keys);
    }

    /** Atom signatures by radius: element and degree at 0, then each mixed with the sum over its neighbours. */
    long[][] signatures(Clock clock) {
      long[][] signature = new long[RADII][atomCount];
      for (int atom = 0; atom < atomCount; atom++)
        signature[0][atom] = mix(symbolHash[atom] * 31L + degree(atom));
      for (int radius = 1; radius < RADII; radius++) {
        for (int atom = 0; atom < atomCount; atom++) {
          clock.tick(1 + degree(atom));
          long sum = 0;
          for (int i = start[atom]; i < start[atom + 1]; i++)
            sum += mix(signature[radius - 1][neighbour[i]]);
          signature[radius][atom] = mix(signature[radius - 1][atom] * GOLDEN_GAMMA + sum);
        }
      }
      return signature;
    }

    /** The most atoms and bonds of a component once the bonds whose label the other molecule lacks are removed. */
    long componentBound() {
      boolean[] seen = new boolean[atomCount];
      int[] queue = new int[atomCount];
      int mostAtoms = 0, mostBonds = 0;
      for (int root = 0; root < atomCount; root++) {
        if (seen[root] || atomClass[root] < 0)
          continue;
        int size = 1, links = 0;
        queue[0] = root;
        seen[root] = true;
        for (int head = 0; head < size; head++) {
          int atom = queue[head];
          for (int place = start[atom]; place < start[atom + 1]; place++) {
            int other = neighbour[place];
            if (label[neighbourBond[place]] < 0)
              continue;
            links++;
            if (!seen[other]) {
              seen[other] = true;
              queue[size++] = other;
            }
          }
        }
        mostAtoms = Math.max(mostAtoms, size);
        // Each bond is met from both ends.
        mostBonds = Math.max(mostBonds, links / 2);
      }
      return score(mostAtoms, mostBonds);
    }
  }

  /**
   * Counts the work of the search in units, about one for each step of an inner loop, and asks
   * whether to stop every {@link #POLL_INTERVAL} units, or at once with {@link #poll()}. The
   * search ends with {@link Stop} when the answer is yes, so every loop of the search ticks the
   * clock.
   */
  static final class Clock {

    /** The units of work between two checks. */
    private static final int POLL_INTERVAL = 1024;
    /** Whether to stop now: the time limit was reached or the search was cancelled. */
    private final BooleanSupplier stop;
    /** The units of work ticked since the last check. */
    private int work;

    /** A clock that stops the search once {@code stop} says so. */
    Clock(BooleanSupplier stop) {
      this.stop = stop;
    }

    /** Counts one unit of work. */
    void tick() {
      tick(1);
    }

    /** Counts {@code units} of work, and checks whether to stop once the interval is reached. */
    void tick(int units) {
      work += units;
      if (work >= POLL_INTERVAL) {
        work = 0;
        poll();
      }
    }

    /** Checks whether to stop now, and throws {@link Stop} if so. */
    void poll() {
      if (stop.getAsBoolean())
        throw new Stop();
    }
  }

  /** Ends a search: the time limit was reached or the search was cancelled. */
  static final class Stop extends RuntimeException {

    /** The serial version, as for any exception. */
    private static final long serialVersionUID = 1L;

    /** A signal, so no message and no stack trace. */
    Stop() {
      super(null, null, false, false);
    }
  }
}

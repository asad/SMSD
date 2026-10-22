/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.concurrent.atomic.AtomicReferenceFieldUpdater;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Matching;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;

/**
 * A molecular graph prepared for substructure and MCS searches: the atom and bond properties,
 * ring and aromaticity flags and stereo descriptors, read once from a CDK
 * {@link IAtomContainer} or supplied through a {@link Builder}, with the labels that searches
 * derive from them on first use.
 */
public final class MolGraph {
  final IAtomContainer mol;
  /** The implicit hydrogen counts, or null for a built graph that left them to the valence model. */
  private final int[] hydrogenCounts;
  final int n;
  public final int[] atomicNum, formalCharge, label, massNumber;
  int[] orbit, canonicalLabel;
  int[][] autGenerators;
  boolean autGeneratorsTruncated;
  long canonicalHash;
  private volatile boolean canonicalComputed = false;
  private volatile boolean ringCountsComputed = false;
  final boolean[] ring, aromatic;
  final int[] degree;
  public final int[] ringCount;
  final int[][] neighbors;
  long[][] adjLong;
  final int words;
  private volatile int[][] cachedNLF1, cachedNLF2, cachedNLF3;
  private volatile int[][] cachedNeighborsByDegDesc;
  private volatile int[] cachedPharmacophoreFeatures;
  /** The ECFP4 and FCFP4 fingerprints of 2048 bits, once computed, see {@link FingerprintEngine}. */
  volatile long[] cachedEcfpDefault, cachedFcfpDefault;
  private static final int SPARSE_THRESHOLD = 200;
  /** The bond properties of a graph up to the dense threshold, by both atoms. */
  final int[][] bondOrdMatrix;
  final boolean[][] bondRingMatrix, bondAromMatrix;
  /**
   * The bond properties of a larger graph: the order, ring flag and aromatic flag of the bond to
   * each neighbour, parallel to {@link #neighbors}.
   */
  final int[][] neighbourBondOrder;
  final boolean[][] neighbourBondRing, neighbourBondAromatic;
  /** The tetrahedral tags, or null when the graph has none. */
  final int[] tetraChirality;
  // CDK winding is relative to this exact ligand order; -1 is the implicit ligand.
  final int[][] tetraLigands;
  /**
   * The double-bond configurations by both atoms for a dense graph, or by neighbour as the bond
   * properties for a larger one; both null when the graph has none.
   */
  final int[][] dbStereoConf, neighbourDbStereo;
  /**
   * The tautomer class of each atom, -1 outside every class. The atoms of one motif match share a
   * class and matches that share an atom join their classes, so the partition does not depend on
   * the atom order; the numbers, from 0, depend on it only for classes that Weisfeiler-Leman
   * refinement cannot tell apart (in practice those a symmetry of the graph exchanges).
   */
  public int[]   tautomerClass;
  /**
   * Baseline per-atom relevance weight for heuristic tautomer confidence, in [0,1].
   * A value of 1.0 gives no penalty. The explicit MCS.computeTautomerConfidence
   * helper applies solvent and pH context when reading these weights.
   */
  public float[] tautomerWeight;
  /**
   * The ionisable site of each atom's tautomer motif, which the pH context of
   * MCS.computeTautomerConfidence reads: {@link #TAUTOMER_SITE_ACID} when the neutral acid HA
   * tautomerises, {@link #TAUTOMER_SITE_BASE} when the neutral base B does, and
   * {@link #TAUTOMER_SITE_NONE} when the motif has no ionisable site or no match of the motif that
   * covers the atom holds a hydrogen the site would release. No motif carries a base site.
   */
  public byte[] tautomerSite;
  /** The aqueous pKa of that site (of the conjugate acid BH+ for a base); 0 when there is none. */
  public float[] tautomerPKa;
  /**
   * Optional external atom IDs for translating between SmsdGraph 0-based
   * contiguous indices and non-contiguous IDs used by editors (e.g. BIME).
   * When null, internal indices are used directly (identity mapping).
   * When set, must have length == n. atomId[i] is the external ID of atom i.
   */
  public int[] atomId;
  public final String name;
  public final String programLine;
  public final String comment;
  public final Map<String, String> properties;
  /** The rings, once perceived, see {@link #rings()}. */
  private volatile RingPerception.Rings rings;
  /** The ring systems of the ring atoms, once found, see {@link #ringSystems()}. */
  private volatile RingPerception.AtomSystems ringSystems;
  /**
   * What SMARTS searches keep for this graph, see {@link Substructure}: null until one reads it.
   * Searches exchange it atomically, so that a kept copy they take out is read by one at a time.
   */
  private volatile Object smartsView;
  private static final AtomicReferenceFieldUpdater<MolGraph, Object> SMARTS_VIEW =
      AtomicReferenceFieldUpdater.newUpdater(MolGraph.class, Object.class, "smartsView");

  /** {@link #tautomerSite}: the motif has no ionisable site. */
  public static final byte TAUTOMER_SITE_NONE = 0;
  /** {@link #tautomerSite}: the neutral acid HA of the motif tautomerises. */
  public static final byte TAUTOMER_SITE_ACID = 1;
  /** {@link #tautomerSite}: the neutral base B of the motif tautomerises. */
  public static final byte TAUTOMER_SITE_BASE = -1;

  // ---------------------------------------------------------------------------
  // Empirical pKa-derived relevance weights of the dominant neutral tautomer; pH
  // enters only through the neutral fraction of the ionisable site below.
  // References: Sitzmann et al. JCIM 2010; Karabina et al. 2014; ChEMBL pKa db
  // ---------------------------------------------------------------------------
  public static final float TW_KETO_ENOL         = 0.98f; // pKa ~20  — keto overwhelmingly dominant
  public static final float TW_AMIDE_IMIDIC      = 0.97f; // pKa ~15-17 — amide form dominant
  public static final float TW_LACTAM_LACTIM     = 0.97f; // cyclic amide — lactam dominant
  public static final float TW_UREA              = 0.96f; // pKa ~14 — urea form dominant
  public static final float TW_PYRIDINONE        = 0.95f; // N-H pKa ~11.7 — lactam dominant
  public static final float TW_THIOAMIDE         = 0.94f; // pKa ~13 — thioamide dominant
  public static final float TW_THIONE_THIOL      = 0.94f; // pKa ~10-11 — thione dominant
  public static final float TW_ENAMINONE         = 0.85f; // vinylogous amide — amine form dominant
  public static final float TW_HYDROXAMIC        = 0.85f; // pKa ~9  — neutral hydroxamic dominant
  public static final float TW_PHENOL_QUINONE    = 0.88f; // pKa ~10 — phenol dominant
  public static final float TW_HYDROXYPYRIMIDINE = 0.90f; // lactam tautomer dominant in N-heterocycles
  public static final float TW_PURINE_NH         = 0.78f; // N9-H dominant in purines (adenine, guanine)
  public static final float TW_DIKETONE_ENOL     = 0.72f; // 1,3-diketone enolisation: pKa ~8-11, pH-sensitive, ~28% enol at pH 7.4
  public static final float TW_NITROSO_OXIME     = 0.95f; // motif T3, the imidic acid N=C-O: amide strongly dominant
  public static final float TW_CYANAMIDE         = 0.80f; // cyanamide dominant over carbodiimide
  public static final float TW_IMIDE             = 0.92f; // imide N-H dominant
  public static final float TW_AMIDINE          = 0.50f;  // symmetric — equal tautomers
  public static final float TW_GUANIDINE        = 0.50f;  // symmetric — equal tautomers
  public static final float TW_IMIDAZOLE_NH     = 0.50f;  // N1H/N3H symmetric in unsubstituted
  public static final float TW_TRIAZOLE_NH      = 0.50f;  // 1H/2H symmetric
  public static final float TW_TETRAZOLE_NH     = 0.50f;  // 1H/2H symmetric
  public static final float TW_NITRO_ACI        = 0.95f;  // nitro form strongly dominant at pH 7.4 (pKa ~3.5)

  // Tier 1b: additional transforms (T16-T30), Dhaked & Nicklaus 2024
  public static final float TW_THIOAMIDE_IMINOTHIOL = 0.95f; // thioamide form dominant (pKa ~13)
  public static final float TW_PHOSPHONATE      = 0.50f;  // symmetric P(=O)(OH) ↔ P(OH)(=O)
  /** @deprecated Use {@link #TW_NITRO_ACI} instead — unified nitro/aci-nitro weight. */
  @Deprecated public static final float TW_NITRO_ACI_NITRO  = TW_NITRO_ACI;
  public static final float TW_15_KETO_ENOL     = 0.75f;  // 1,5-shift through conjugation
  public static final float TW_FURANOSE_PYRANOSE = 0.60f;  // ring-chain: furanose/pyranose
  public static final float TW_LACTOL           = 0.65f;  // ring-chain: lactol
  public static final float TW_PHENOL_QUINONE_METHIDE = 0.92f; // phenol → para-quinone methide
  public static final float TW_TRIAZOLE_NH_SHIFT = 0.50f; // 1,2,3-triazole NH shift (symmetric)
  public static final float TW_BENZIMIDAZOLE_NH = 0.50f;  // NH-1 ↔ NH-3 shift (symmetric)
  public static final float TW_PYRIDONE_HYDROXYPYRIDINE = TW_PYRIDINONE; // unified with TW_PYRIDINONE — both are pyridone/hydroxypyridine equilibria
  public static final float TW_BARBITURIC       = 0.80f;  // tri-keto ↔ enol forms
  public static final float TW_ALLYL_SHIFT      = 0.70f;  // X-C=C-C ↔ C=C-C-X (X=OH,SH,NH)
  public static final float TW_SELENOL_SELENOKETONE = 0.94f; // C=Se ↔ C(-SeH)

  // ---------------------------------------------------------------------------
  // Aqueous pKa of the ionisable site of each motif that has one, from the ranges
  // above (midpoints where a range is given), applied in every solvent. Every site is
  // an acid, which a match of its motif carries only when the atom that would release
  // the proton bears a hydrogen. The pKa ~3.5 of the nitro/aci-nitro motif belongs to
  // the minor aci-nitro tautomer, not to the dominant nitro form from which the motif
  // arises, so that motif has no site.
  // ---------------------------------------------------------------------------
  public static final float TPKA_KETO_ENOL            = 20.0f; // ketone alpha C-H
  public static final float TPKA_AMIDE_IMIDIC         = 16.0f; // amide N-H
  public static final float TPKA_UREA                 = 14.0f; // urea N-H
  public static final float TPKA_PYRIDINONE           = 11.7f; // 2-pyridone N-H (8.6 for pyrimidin-4(3H)-one)
  public static final float TPKA_PYRIDONE_HYDROXYPYRIDINE = TPKA_PYRIDINONE; // same equilibrium
  public static final float TPKA_THIOAMIDE            = 13.0f; // thioamide N-H
  public static final float TPKA_THIONE_THIOL         = 10.5f; // thione/thiol
  public static final float TPKA_HYDROXAMIC           = 9.0f;  // hydroxamic acid O-H
  public static final float TPKA_PHENOL_QUINONE       = 10.0f; // phenol O-H
  public static final float TPKA_DIKETONE_ENOL        = 9.5f;  // 1,3-diketone C-H
  public static final float TPKA_IMIDE                = 9.5f;  // imide N-H
  public static final float TPKA_THIOAMIDE_IMINOTHIOL = 13.0f; // thioamide N-H

  /** Whether atoms {@code i} and {@code j} are bonded; every present bond has an order of at least one. */
  public boolean hasBond(int i, int j) { return (adjLong[i][j >>> 6] & 1L << (j & 63)) != 0; }

  /**
   * Whether an atom of this element with three neighbours and no hydrogen is a tetrahedral centre
   * whose lone pair is the fourth ligand: N, P, S, As, Se, Sb and Te.
   */
  static boolean lonePairCentre(int atomicNumber) {
    return atomicNumber == 7 || atomicNumber == 15 || atomicNumber == 16 || atomicNumber == 33
        || atomicNumber == 34 || atomicNumber == 51 || atomicNumber == 52;
  }

  /**
   * The ligand order of the tetrahedral tag of {@code atom}, with -1 for the implicit ligand of
   * three neighbours (a hydrogen or a lone pair), or null when the atom is not a tetrahedral
   * centre. It is a centre when it has four neighbours and no implicit hydrogen, or three
   * neighbours and either one implicit hydrogen or, for an element of {@link #lonePairCentre},
   * none. The order is the one of the CDK stereo element the graph was built from; otherwise it
   * is the neighbours in order, with the implicit ligand second when the first neighbour precedes
   * the atom, as after the atom it was written from in SMILES, and first otherwise. The native
   * core applies the same rule.
   */
  int[] tetrahedralLigands(int atom) {
    int degree = neighbors[atom].length;
    if (degree != 3 && degree != 4) return null;
    int hydrogens = hydrogenCount(atom);
    boolean centre = degree == 4 ? hydrogens == 0
        : hydrogens == 1 || hydrogens == 0 && lonePairCentre(atomicNum[atom]);
    if (!centre) return null;
    if (tetraLigands != null && tetraLigands[atom] != null) return tetraLigands[atom];
    int[] order = new int[4];
    int implicit = degree == 4 ? -1 : neighbors[atom][0] < atom ? 1 : 0;
    for (int slot = 0, position = 0; slot < 4; slot++)
      order[slot] = slot == implicit ? -1 : neighbors[atom][position++];
    return order;
  }

  public int atomCount() { return n; }

  /** The number of bonds: half the sum of the degrees. */
  long bondCount() {
    long ends = 0;
    for (int d : degree) ends += d;
    return ends / 2;
  }

  /**
   * The implicit hydrogen count of atom {@code idx}, as the native core stores it: the count the
   * CDK atom supplied when the graph was built, or else its valence or the valence model less its
   * bonds; for a graph from the {@link Builder}, the supplied count, or else the SMILES valence
   * model. Aromatic bonds are stored as order 4, so the count cannot be read off bond order sums.
   */
  public int hydrogenCount(int idx) {
    return hydrogenCounts != null ? hydrogenCounts[idx] : defaultHydrogenCount(idx);
  }

  /** The hydrogen count of a CDK atom: the supplied count, or its valence or the valence model less its bonds. */
  private static int cdkHydrogenCount(IAtomContainer mol, IAtom atom, int atomicNumber, int charge) {
    Integer supplied = atom.getImplicitHydrogenCount();
    if (supplied != null) return supplied;
    int bondOrderSum = 0;
    for (IBond bond : mol.getConnectedBondsList(atom)) {
      IBond.Order order = bond.getOrder();
      bondOrderSum += order == null || order == IBond.Order.UNSET ? 1 : order.numeric();
    }
    Integer valence = atom.getValency();
    if (valence != null) return Math.max(0, valence - bondOrderSum);
    int unpaired = mol.getConnectedSingleElectronsCount(atom);
    return atomicNumber <= 1 ? 0 : defaultHydrogenCount(atomicNumber, false, bondOrderSum + unpaired, charge);
  }

  /**
   * CDK reads a ring system written with ':' bonds between atoms written in capitals, such as
   * C1:C:C:C:C:C:1, as aromatic single bonds without a Kekulé structure and gives its atoms the
   * hydrogens of single bonds. In {@code hydrogens}, indexed as the atoms of {@code mol}, whose ring
   * flags must be set and whose atoms {@code index} numbers, those atoms take the hydrogens of
   * aromatic atoms instead, as the native SMILES reader and RDKit give them, except that a neutral
   * nitrogen or phosphorus keeps its hydrogen where the system needs it for a Kekulé structure, as
   * in pyrrole. Returns whether a count changed.
   */
  static boolean aromaticBondHydrogens(IAtomContainer mol, Map<IAtom, Integer> index, int[] hydrogens) {
    int n = mol.getAtomCount();
    int[] system = null;
    for (IBond bond : mol.bonds())
      if (bond.isAromatic() && bond.isInRing()) {
        if (system == null) system = UnionFind.create(n);
        UnionFind.join(system, index.get(bond.getBegin()), index.get(bond.getEnd()));
      }
    if (system == null) return false;
    int[] bondOrderSum = new int[n];
    boolean[] onSystem = new boolean[n], localised = new boolean[n];
    for (IBond bond : mol.bonds()) {
      int first = index.get(bond.getBegin()), second = index.get(bond.getEnd());
      IBond.Order order = bond.getOrder();
      int value = order == null || order == IBond.Order.UNSET ? 1 : order.numeric();
      bondOrderSum[first] += value;
      bondOrderSum[second] += value;
      if (!bond.isAromatic() || !bond.isInRing()) continue;
      onSystem[first] = onSystem[second] = true;
      // A double or triple aromatic bond means CDK gave the system a Kekulé structure.
      if (value > 1) localised[UnionFind.root(system, first)] = true;
    }
    // The nitrogen and phosphorus atoms whose hydrogen count changes, by system.
    Map<Integer, List<Integer>> heteroatoms = new TreeMap<>();
    boolean changed = false;
    for (int atom = 0; atom < n; atom++) {
      IAtom value = mol.getAtom(atom);
      if (!onSystem[atom] || localised[UnionFind.root(system, atom)] || !value.isAromatic()
          || value.getImplicitHydrogenCount() == null || value.getAtomicNumber() == null
          || (value.getFormalCharge() != null && value.getFormalCharge() != 0)) continue;
      int single = defaultHydrogenCount(value.getAtomicNumber(), false, bondOrderSum[atom], 0);
      int aromatic = defaultHydrogenCount(value.getAtomicNumber(), true, bondOrderSum[atom], 0);
      if (hydrogens[atom] != single || aromatic >= single) continue;
      hydrogens[atom] = aromatic;
      changed = true;
      List<Integer> list = heteroatoms.computeIfAbsent(UnionFind.root(system, atom), root -> new ArrayList<>());
      if (value.getAtomicNumber() == 7 || value.getAtomicNumber() == 15) list.add(atom);
    }
    if (!changed) return false;
    // A Kekulé structure may place a double bond on any bond between atoms of one system.
    List<List<Integer>> adjacency = new ArrayList<>(n);
    for (int atom = 0; atom < n; atom++) adjacency.add(new ArrayList<>());
    for (IBond bond : mol.bonds()) {
      int first = index.get(bond.getBegin()), second = index.get(bond.getEnd());
      if (!onSystem[first] || !onSystem[second] || UnionFind.root(system, first) != UnionFind.root(system, second))
        continue;
      adjacency.get(first).add(second);
      adjacency.get(second).add(first);
    }
    int[][] graph = new int[n][];
    for (int atom = 0; atom < n; atom++)
      graph[atom] = adjacency.get(atom).stream().mapToInt(Integer::intValue).toArray();
    long[] invariant = null;
    for (Map.Entry<Integer, List<Integer>> entry : heteroatoms.entrySet()) {
      int root = entry.getKey();
      List<Integer> list = entry.getValue();
      if (kekulisable(mol, system, root, bondOrderSum, hydrogens, graph)) continue;
      // Without a Kekulé structure, one heteroatom, else a pair, keeps its single-bond hydrogen
      // where that gives one; of several, the first by an atom invariant, so that the choice does
      // not depend on the atom order. With none, the atoms keep the aromatic reading, as natively.
      if (invariant == null) invariant = atomInvariants(mol, index);
      final long[] rank = invariant;
      list.sort((first, second) -> rank[first] != rank[second] ? Long.compare(rank[first], rank[second])
          : Integer.compare(first, second));
      boolean kekule = false;
      for (int first = 0; first < list.size() && !kekule; first++)
        kekule = keepsHydrogen(mol, system, root, bondOrderSum, hydrogens, graph, list.get(first), -1);
      for (int first = 0; first < list.size() && !kekule; first++)
        for (int second = first + 1; second < list.size() && !kekule; second++)
          kekule = keepsHydrogen(mol, system, root, bondOrderSum, hydrogens, graph, list.get(first), list.get(second));
    }
    return true;
  }

  /**
   * Atom invariants that do not depend on the atom order: element, charge, degree, hydrogens and
   * aromaticity, refined by the invariants of the neighbours until the classes stop splitting.
   */
  private static long[] atomInvariants(IAtomContainer mol, Map<IAtom, Integer> index) {
    int n = mol.getAtomCount();
    int[][] neighbours = new int[n][];
    for (int atom = 0; atom < n; atom++) {
      List<IAtom> adjacent = mol.getConnectedAtomsList(mol.getAtom(atom));
      neighbours[atom] = new int[adjacent.size()];
      for (int k = 0; k < adjacent.size(); k++) neighbours[atom][k] = index.get(adjacent.get(k));
    }
    long[] invariant = new long[n];
    for (int atom = 0; atom < n; atom++) {
      IAtom value = mol.getAtom(atom);
      Integer element = value.getAtomicNumber(), charge = value.getFormalCharge(), hydrogens = value.getImplicitHydrogenCount();
      invariant[atom] = (((element == null ? 0 : element) * 31L + (charge == null ? 0 : charge)) * 31L
          + neighbours[atom].length) * 31L + (hydrogens == null ? -1 : hydrogens) * 2L + (value.isAromatic() ? 1 : 0);
    }
    for (int classes = new HashSet<>(toList(invariant)).size(); ; ) {
      long[] refined = new long[n];
      for (int atom = 0; atom < n; atom++) {
        long[] adjacent = new long[neighbours[atom].length];
        for (int k = 0; k < adjacent.length; k++) adjacent[k] = invariant[neighbours[atom][k]];
        Arrays.sort(adjacent);
        long hash = invariant[atom];
        for (long value : adjacent) hash = hash * 1000003L + value;
        refined[atom] = hash;
      }
      int split = new HashSet<>(toList(refined)).size();
      if (split <= classes) return invariant;
      invariant = refined;
      classes = split;
    }
  }

  private static List<Long> toList(long[] values) {
    List<Long> list = new ArrayList<>(values.length);
    for (long value : values) list.add(value);
    return list;
  }

  /**
   * Gives the ring systems of {@code mol} that CDK read from ':' bonds between atoms written in
   * capitals the hydrogens {@link #aromaticBondHydrogens} gives them and a Kekulé structure where
   * one exists; the ring flags of {@code mol} must be set.
   */
  static void kekuliseAromaticBondSystems(IAtomContainer mol) {
    int n = mol.getAtomCount();
    IdentityHashMap<IAtom, Integer> index = new IdentityHashMap<>(n);
    int[] hydrogens = new int[n];
    for (int i = 0; i < n; i++) {
      IAtom atom = mol.getAtom(i);
      index.put(atom, i);
      hydrogens[i] = cdkHydrogenCount(mol, atom, atom.getAtomicNumber() == null ? 0 : atom.getAtomicNumber(),
          atom.getFormalCharge() == null ? 0 : atom.getFormalCharge());
    }
    if (!aromaticBondHydrogens(mol, index, hydrogens)) return;
    for (int i = 0; i < n; i++) {
      Integer supplied = mol.getAtom(i).getImplicitHydrogenCount();
      if (supplied != null && supplied != hydrogens[i]) mol.getAtom(i).setImplicitHydrogenCount(hydrogens[i]);
    }
    try {
      Kekulization.kekulize(mol);
    } catch (CDKException noKekuleStructure) {
      // The system keeps its single bonds, as it has no Kekulé structure with these hydrogens.
    }
  }

  /**
   * Gives {@code first} and {@code second} (unless -1) back the hydrogen of their single-bond
   * reading, one more than the aromatic one, when that gives the system a Kekulé structure.
   */
  private static boolean keepsHydrogen(IAtomContainer mol, int[] system, int root, int[] bondOrderSum,
      int[] hydrogens, int[][] graph, int first, int second) {
    hydrogens[first]++;
    if (second >= 0) hydrogens[second]++;
    if (kekulisable(mol, system, root, bondOrderSum, hydrogens, graph)) return true;
    hydrogens[first]--;
    if (second >= 0) hydrogens[second]--;
    return false;
  }

  /**
   * Whether the atoms of one aromatic ring system that need a double bond, those one bond short
   * of the valence model with their hydrogens, can each take one on a bond of the system.
   */
  private static boolean kekulisable(IAtomContainer mol, int[] system, int root, int[] bondOrderSum,
      int[] hydrogens, int[][] graph) {
    BitSet needed = new BitSet();
    for (int atom = 0; atom < graph.length; atom++) {
      if (graph[atom].length == 0 || UnionFind.root(system, atom) != root) continue;
      IAtom value = mol.getAtom(atom);
      if (value.getAtomicNumber() == null) return true;
      int charge = value.getFormalCharge() == null ? 0 : value.getFormalCharge();
      if (defaultHydrogenCount(value.getAtomicNumber(), false, bondOrderSum[atom], charge) == hydrogens[atom] + 1)
        needed.set(atom);
    }
    return needed.cardinality() % 2 == 0 && Matching.withCapacity(graph.length).perfect(graph, needed);
  }

  /** Implicit hydrogens plus explicitly represented hydrogen neighbours. */
  long attachedHydrogenCount(int atom) {
    long count = hydrogenCount(atom);
    for (int neighbor : neighbors[atom]) if (atomicNum[neighbor] == 1) count++;
    return count;
  }

  /** Hydrogen count implied by an unbracketed atom's valence. */
  int defaultHydrogenCount(int idx) {
    if (atomicNum[idx] <= 1) return 0;
    int bondOrdSum = 0;
    for (int nb : neighbors[idx]) {
      int bo = bondOrder(idx, nb);
      bondOrdSum += (bo == 4) ? 1 : bo;
    }
    return defaultHydrogenCount(atomicNum[idx], aromatic[idx], bondOrdSum, formalCharge[idx]);
  }

  /**
   * Default SMILES valences per element, in ascending order: the one Java valence table, mirroring
   * {@code defaultValences} in {@code smiles_utils.hpp}. Elements without an entry take no
   * implicit hydrogens.
   */
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
   * Hydrogen count implied by the ordinary SMILES valence model, the one Java form of
   * {@code computeImplicitH} in {@code smiles_utils.hpp}. Aromatic bonds count once in
   * {@code bondOrdSum} and an aromatic atom counts one more; a neutral aromatic atom already at
   * its lowest valence takes no hydrogens. The smallest default valence that fits, adjusted for
   * charge (cations of N, P, O, S, As and Se and anions of B and Al gain |charge|, other charged
   * atoms lose it), sets the count.
   */
  static int defaultHydrogenCount(int atomicNumber, boolean aromatic, int bondOrdSum, int charge) {
    int[] valences = atomicNumber >= 0 && atomicNumber < DEFAULT_VALENCES.length
        ? DEFAULT_VALENCES[atomicNumber] : null;
    if (valences == null) return 0;
    if (aromatic && charge == 0 && bondOrdSum >= valences[0]) return 0;
    long effective = (long) bondOrdSum + (aromatic ? 1 : 0);
    long absoluteCharge = Math.abs((long) charge);
    boolean expands = (charge > 0 && (atomicNumber == 7 || atomicNumber == 8 || atomicNumber == 15
        || atomicNumber == 16 || atomicNumber == 33 || atomicNumber == 34))
        || (charge < 0 && (atomicNumber == 5 || atomicNumber == 13));
    for (int valence : valences) {
      long hydrogens = (expands ? valence + absoluteCharge : valence - absoluteCharge) - effective;
      if (hydrogens >= 0) return (int) Math.min(hydrogens, Integer.MAX_VALUE);
    }
    return 0;
  }

  /**
   * Translate an MCS mapping from internal 0-based contiguous indices to external
   * atom IDs. If either graph has an {@code atomId} array set, the corresponding
   * keys or values in the returned map are translated. If neither graph has external
   * IDs, the input mapping is returned unchanged.
   *
   * @param mapping internal-index mapping (g1 index to g2 index)
   * @param g1      the first molecule graph (keys)
   * @param g2      the second molecule graph (values)
   * @return mapping with keys/values translated to external atom IDs where available
   */
  public static Map<Integer, Integer> translateToAtomIds(
      Map<Integer, Integer> mapping, MolGraph g1, MolGraph g2) {
    if (mapping == null || mapping.isEmpty()) return mapping;
    if (g1.atomId == null && g2.atomId == null) return mapping;
    Map<Integer, Integer> translated = new LinkedHashMap<>(mapping.size());
    for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
      int k = e.getKey(), v = e.getValue();
      if (g1.atomId != null && (k < 0 || k >= g1.atomId.length)) continue;
      if (g2.atomId != null && (v < 0 || v >= g2.atomId.length)) continue;
      int key = g1.atomId != null ? g1.atomId[k] : k;
      int val = g2.atomId != null ? g2.atomId[v] : v;
      translated.put(key, val);
    }
    return translated;
  }

  /** The orbit of each atom under the chemical automorphism group, as a copy. */
  public int[] getOrbit() { ensureCanonical(); return orbit.clone(); }

  public int[][] computeRings() { return rings().smallest; }

  /**
   * The rings, perceived on first use. Searches on several threads may share a graph, so both
   * sets are published together, once, through a volatile field: a thread sees them whole or
   * waits for the thread that perceives them.
   */
  private RingPerception.Rings rings() {
    RingPerception.Rings perceived = rings;
    if (perceived != null) return perceived;
    synchronized (this) {
      if (rings == null) rings = RingPerception.perceive(n, neighbors);
      return rings;
    }
  }

  public int numRings() { return computeRings().length; }

  /**
   * Whether the ring perception behind {@link #computeRings()} and
   * {@link #computeRelevantCycles()} stopped its exhaustive candidate search at its work bound,
   * as for a very large ring system, so that the rings come from the bounded perception and may
   * be fewer than the graph holds. The rings themselves are the same whether or not this is
   * asked. The native core reports the same through {@code smsd::ringsTruncated()}.
   *
   * @return true when the exhaustive candidate search stopped at its work bound
   */
  public boolean ringsTruncated() { return rings().truncated; }

  /**
   * The ring systems of the ring atoms, found on first use and published once, as the rings
   * are.
   */
  private RingPerception.AtomSystems ringSystems() {
    RingPerception.AtomSystems found = ringSystems;
    if (found != null) return found;
    synchronized (this) {
      if (ringSystems == null) ringSystems = RingPerception.atomSystems(this);
      return ringSystems;
    }
  }

  /** Store {@code value} as what SMARTS searches keep for this graph, and return the previous value. */
  Object exchangeSmartsView(Object value) { return SMARTS_VIEW.getAndSet(this, value); }

  /** Store {@code value} as what SMARTS searches keep for this graph if that is still {@code expected}. */
  boolean replaceSmartsView(Object expected, Object value) {
    return SMARTS_VIEW.compareAndSet(this, expected, value);
  }

  public int ringSystemOf(int atom) { return ringSystems().ofAtom[atom]; }
  public long ringSystemSig(int rsId) {
    long[] signature = ringSystems().signature;
    return rsId >= 0 && rsId < signature.length ? signature[rsId] : 0;
  }
  public int ringSystemCount() { return ringSystems().signature.length; }

  /**
   * Computes the set of relevant cycles (union of all minimum cycle bases).
   * A cycle is relevant if it cannot be expressed as the XOR (symmetric
   * difference) of strictly shorter cycles. This set is unique for a given
   * molecule -- unlike a single MCB/SSSR which may be arbitrary.
   *
   * <p>Algorithm based on Vismara (1997) and Kolodzik et al. (2012):
   * <ol>
   *   <li>Compute a single MCB via {@link #computeRings()} (Horton's algorithm)</li>
   *   <li>Enumerate all short cycles up to the maximum MCB cycle length</li>
   *   <li>A cycle is relevant iff it is not a linear combination (over GF(2))
   *       of strictly shorter cycles</li>
   * </ol>
   *
   * @return array of relevant cycles, each cycle given as an array of atom indices
   */
  public int[][] computeRelevantCycles() {
    return rings().relevant;
  }


  /**
   * Computes Unique Ring Families (URFs) as described by Kolodzik, Urbaczek
   * &amp; Rarey (2012). Two relevant cycles belong to the same URF if one can
   * be obtained from the other by a single edge exchange -- i.e. their
   * symmetric difference (XOR) has exactly two edges that form a path of
   * length 2 (sharing a common vertex).
   *
   * @return list of URFs, each URF is a list of cycles (each cycle is an array of atom indices)
   */
  public int[][][] computeURFs() {
    return RingPerception.uniqueRingFamilies(this);
  }

  /**
   * Return the Smallest Set of Smallest Rings (minimum cycle basis),
   * pre-sorted ascending by ring size.
   *
   * <p>Static convenience method that delegates to {@link #computeRings()}
   * and applies an explicit ascending-size sort to satisfy the SSSR contract.
   *
   * @param g the molecular graph
   * @return SSSR rings sorted ascending by size
   */
  public static List<List<Integer>> computeSSSR(MolGraph g) {
    return RingPerception.smallestRingsBySize(g);
  }

  /**
   * Return SSSR rings optimised for 2D coordinate generation:
   * <ul>
   *   <li>Largest ring system first</li>
   *   <li>Fused rings ordered by shared-edge adjacency within each system</li>
   *   <li>Within each system, rings sorted by size (smallest first for placement)</li>
   * </ul>
   *
   * @param g the molecular graph
   * @return layout-ordered SSSR rings
   */
  public static List<List<Integer>> layoutSSSR(MolGraph g) {
    return RingPerception.layoutOrder(g);
  }

  /**
   * 2D point for layout coordinates.
   */
  public static final class Point2D {
    public double x, y;
    public Point2D(double x, double y) { this.x = x; this.y = y; }
  }

  /**
   * Reduce bond crossings in a 2D layout via two-phase simulated annealing.
   *
   * <p><b>Phase 1 — System-level flipping:</b> Flip entire ring systems
   * (connected components of fused rings) about their centroid to find
   * the best global orientation.</p>
   *
   * <p><b>Phase 2 — Individual ring flipping:</b> Within multi-ring systems,
   * flip individual SSSR rings while keeping shared (fusion) atoms fixed
   * as pivots.  This resolves crossings that system-level moves cannot reach.</p>
   *
   * @param g        the molecular graph
   * @param coords   input/output 2D coordinates (coords[atomIdx])
   * @param maxIter  maximum total SA iterations across both phases
   * @return number of bond crossings remaining after optimisation
   */
  public static int reduceCrossings(MolGraph g, Point2D[] coords, int maxIter) {
    return GraphLayout.reduceCrossings(g, coords, maxIter);
  }

  /** Convenience overload: default maxIter = 1000. */
  public static int reduceCrossings(MolGraph g, Point2D[] coords) {
    return reduceCrossings(g, coords, 1000);
  }

  // ========================================================================
  // Force-directed layout
  // ========================================================================

  /**
   * Force-directed 2D layout minimisation.
   * Iteratively moves atoms to minimise bond-length stress, non-bonded
   * repulsion, and crossing penalties.
   *
   * @param g                the molecular graph
   * @param coords           input/output 2D coordinates
   * @param maxIter          maximum iterations (default 500)
   * @param targetBondLength desired bond length (default 1.5)
   * @return final stress energy
   */
  public static double forceDirectedLayout(MolGraph g, Point2D[] coords,
                                            int maxIter, double targetBondLength) {
    return GraphLayout.forceDirectedLayout(g, coords, maxIter, targetBondLength);
  }

  /** Convenience overload with defaults: maxIter=500, targetBondLength=1.5. */
  public static double forceDirectedLayout(MolGraph g, Point2D[] coords) {
    return forceDirectedLayout(g, coords, 500, 1.5);
  }

  // ========================================================================
  // Stress majorisation (SMACOF)
  // ========================================================================

  /**
   * Stress majorisation (SMACOF algorithm) for 2D layout.
   * Minimises: sum_{i&lt;j} w_ij * (d_ij - D_ij)^2
   * where D_ij = graph_distance(i,j) * targetBondLength.
   *
   * @param g                the molecular graph
   * @param coords           input/output 2D coordinates
   * @param maxIter          maximum iterations (default 300)
   * @param targetBondLength desired bond length (default 1.5)
   * @return final normalised stress value
   */
  public static double stressMajorisation(MolGraph g, Point2D[] coords,
                                           int maxIter, double targetBondLength) {
    return GraphLayout.stressMajorisation(g, coords, maxIter, targetBondLength);
  }

  /** Convenience overload with defaults: maxIter=300, targetBondLength=1.5. */
  public static double stressMajorisation(MolGraph g, Point2D[] coords) {
    return stressMajorisation(g, coords, 300, 1.5);
  }

  // ========================================================================
  // Template matching for common scaffolds
  // ========================================================================

  /**
   * Check if a molecule matches a known scaffold template.
   * Returns pre-computed coordinates if matched, null otherwise.
   *
   * Supported scaffolds: benzene, naphthalene, indole, purine, steroid,
   * morphinan, quinoline, biphenyl, cyclohexane, piperidine.
   *
   * @param g                the molecular graph
   * @param targetBondLength desired bond length for coordinate scaling (default 1.5)
   * @return array of Point2D coordinates if matched, null otherwise
   */
  public static Point2D[] matchTemplate(MolGraph g, double targetBondLength) {
    return GraphLayout.matchTemplate(g, targetBondLength);
  }

  /** Convenience overload: default targetBondLength = 1.5. */
  public static Point2D[] matchTemplate(MolGraph g) {
    return matchTemplate(g, 1.5);
  }

  int[][] getNLF1() {
    if (cachedNLF1 == null) cachedNLF1 = buildAllNLF(this, MolGraph::buildNLF1);
    return cachedNLF1;
  }

  /**
   * Cached neighbours sorted by descending degree for molecule query ordering.
   * Computed once per MolGraph, reused across all matcher constructions.
   */
  int[][] getNeighborsByDegDesc() {
    if (cachedNeighborsByDegDesc == null) {
      int[][] result = new int[n][];
      for (int i = 0; i < n; i++) {
        int[] nb = neighbors[i].clone();
        // Insertion sort by descending degree
        for (int a = 1; a < nb.length; a++) {
          int key = nb[a], keyDeg = degree[key]; int b = a - 1;
          while (b >= 0 && degree[nb[b]] < keyDeg) { nb[b + 1] = nb[b]; b--; }
          nb[b + 1] = key;
        }
        result[i] = nb;
      }
      cachedNeighborsByDegDesc = result;
    }
    return cachedNeighborsByDegDesc;
  }

  int[][] getNLF2() {
    if (cachedNLF2 == null) cachedNLF2 = buildAllNLF(this, MolGraph::buildNLF2);
    return cachedNLF2;
  }

  int[][] getNLF3() {
    if (cachedNLF3 == null) cachedNLF3 = buildAllNLF(this, MolGraph::buildNLF3);
    return cachedNLF3;
  }

  /**
   * Cached pharmacophore feature classification per atom for FCFP fingerprints.
   * Computed once per MolGraph, eliminating O(n×degree²) recomputation
   * on every FCFP call.
   */
  int[] getPharmacophoreFeatures() {
    if (cachedPharmacophoreFeatures == null) {
      int[] features = new int[n];
      for (int i = 0; i < n; i++)
        features[i] = FingerprintEngine.classifyPharmacophore(this, i);
      cachedPharmacophoreFeatures = features;
    }
    return cachedPharmacophoreFeatures;
  }

  /** The position of {@code j} among the neighbours of {@code i}, or -1 when they are not bonded. */
  private int slot(int i, int j) {
    int[] row = neighbors[i];
    for (int k = 0; k < row.length; k++) if (row[k] == j) return k;
    return -1;
  }

  int dbStereo(int i, int j) {
    if (dbStereoConf != null) return dbStereoConf[i][j];
    if (neighbourDbStereo == null) return 0;
    int k = slot(i, j);
    return k < 0 ? 0 : neighbourDbStereo[i][k];
  }

  /** The stored double-bond relation uses the first explicit ligand at each end. */
  int doubleBondReference(int atom, int otherEnd) {
    for (int neighbor : neighbors[atom]) if (neighbor != otherEnd) return neighbor;
    return -1;
  }

  /** Check if the molecule already has ring membership flags set (e.g. by Cycles.markRingAtomsAndBonds). */
  private static boolean hasRingFlags(IAtomContainer mol) {
    // Only check isInRing() — NOT isAromatic(). CDK's SmilesParser sets aromatic flags
    // from lowercase SMILES but does NOT set ring membership flags. We must call
    // Cycles.markRingAtomsAndBonds() even if atoms are marked aromatic.
    int n = mol.getAtomCount();
    int check = Math.min(n, 8);
    for (int i = 0; i < check; i++) {
      if (mol.getAtom(i).isInRing()) return true;
    }
    int bCheck = Math.min(mol.getBondCount(), 8);
    for (int i = 0; i < bCheck; i++) {
      if (mol.getBond(i).isInRing()) return true;
    }
    return false;
  }

  static long bondKey(int i, int j) {
    return ((long) Math.min(i, j) << 32) | Math.max(i, j);
  }

  private static int encodeBondOrder(IBond b) {
    if (b == null) return 0;
    if (b.isAromatic()) return 4;
    IBond.Order ord = b.getOrder();
    if (ord == IBond.Order.SINGLE) return 1;
    if (ord == IBond.Order.DOUBLE) return 2;
    if (ord == IBond.Order.TRIPLE) return 3;
    return 1;
  }

  int bondOrder(int i, int j) {
    if (bondOrdMatrix != null) return bondOrdMatrix[i][j];
    int k = slot(i, j);
    return k < 0 ? 0 : neighbourBondOrder[i][k];
  }

  boolean bondInRing(int i, int j) {
    if (bondRingMatrix != null) return bondRingMatrix[i][j];
    int k = slot(i, j);
    return k >= 0 && neighbourBondRing[i][k];
  }

  boolean bondAromatic(int i, int j) {
    if (bondAromMatrix != null) return bondAromMatrix[i][j];
    int k = slot(i, j);
    return k >= 0 && neighbourBondAromatic[i][k];
  }

  public MolGraph(IAtomContainer mol) {
    this.mol = mol;
    this.n = mol.getAtomCount();
    this.words = CommonBondGraphs.words(n);
    Object title = mol.getProperty("cdk:Title");
    this.name = title != null ? String.valueOf(title) : "";
    this.programLine = "";
    this.comment = "";
    this.properties = new LinkedHashMap<>();
    for (Map.Entry<Object, Object> e : mol.getProperties().entrySet()) {
      if (e.getKey() == null || e.getValue() == null) continue;
      this.properties.put(String.valueOf(e.getKey()), String.valueOf(e.getValue()));
    }

    // Ensure ring/aromaticity flags are set (SmilesParser may not call ring perception).
    // Skip the expensive CDK ring perception if the molecule already has ring flags
    // (e.g., from AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms + Aromaticity.apply).
    if (!hasRingFlags(mol)) {
      try { RingPerception.markRingAtomsAndBonds(mol); } catch (Exception ignored) {}
    }

    IdentityHashMap<IAtom, Integer> idxMap = new IdentityHashMap<>(n);
    for (int i = 0; i < n; i++) idxMap.put(mol.getAtom(i), i);

    this.atomicNum = new int[n];
    this.formalCharge = new int[n];
    this.massNumber = new int[n];
    this.label = new int[n];
    this.ring = new boolean[n];
    this.aromatic = new boolean[n];
    this.degree = new int[n];

    for (int i = 0; i < n; i++) {
      IAtom a = mol.getAtom(i);
      Integer az = a.getAtomicNumber();
      atomicNum[i] = az != null ? az : 0;
      Integer fc = a.getFormalCharge();
      formalCharge[i] = fc != null ? fc : 0;
      Integer mn = a.getMassNumber();
      massNumber[i] = mn != null ? mn : 0;
      ring[i] = a.isInRing();
      aromatic[i] = a.isAromatic();
      label[i] = (atomicNum[i] << 2) | (aromatic[i] ? 2 : 0) | (ring[i] ? 1 : 0);
    }
    this.hydrogenCounts = new int[n];
    for (int i = 0; i < n; i++)
      hydrogenCounts[i] = cdkHydrogenCount(mol, mol.getAtom(i), atomicNum[i], formalCharge[i]);
    aromaticBondHydrogens(mol, idxMap, hydrogenCounts);

    // The neighbours and bond properties come from one pass over the bonds, without the list a
    // container allocates for each atom.
    int[] degreeCount = new int[n];
    for (IBond b : mol.bonds()) {
      degreeCount[idxMap.get(b.getAtom(0))]++;
      degreeCount[idxMap.get(b.getAtom(1))]++;
    }
    this.neighbors = new int[n][];
    for (int i = 0; i < n; i++) { neighbors[i] = new int[degreeCount[i]]; degree[i] = 0; }
    boolean dense = n <= SPARSE_THRESHOLD;
    this.bondOrdMatrix = dense ? new int[n][n] : null;
    this.bondRingMatrix = dense ? new boolean[n][n] : null;
    this.bondAromMatrix = dense ? new boolean[n][n] : null;
    this.neighbourBondOrder = dense ? null : new int[n][];
    this.neighbourBondRing = dense ? null : new boolean[n][];
    this.neighbourBondAromatic = dense ? null : new boolean[n][];
    if (!dense)
      for (int i = 0; i < n; i++) {
        neighbourBondOrder[i] = new int[degreeCount[i]];
        neighbourBondRing[i] = new boolean[degreeCount[i]];
        neighbourBondAromatic[i] = new boolean[degreeCount[i]];
      }
    for (IBond b : mol.bonds()) {
      int a = idxMap.get(b.getAtom(0)), c = idxMap.get(b.getAtom(1));
      int ord = encodeBondOrder(b);
      boolean inRing = b.isInRing(), arom = b.isAromatic();
      if (dense) {
        bondOrdMatrix[a][c] = ord; bondOrdMatrix[c][a] = ord;
        bondRingMatrix[a][c] = inRing; bondRingMatrix[c][a] = inRing;
        bondAromMatrix[a][c] = arom; bondAromMatrix[c][a] = arom;
      } else {
        neighbourBondOrder[a][degree[a]] = ord; neighbourBondOrder[c][degree[c]] = ord;
        neighbourBondRing[a][degree[a]] = inRing; neighbourBondRing[c][degree[c]] = inRing;
        neighbourBondAromatic[a][degree[a]] = arom; neighbourBondAromatic[c][degree[c]] = arom;
      }
      neighbors[a][degree[a]++] = c;
      neighbors[c][degree[c]++] = a;
    }

    // Stereo is stored only when the molecule has stereo elements.
    @SuppressWarnings("unchecked")
    Iterable<IStereoElement<?, ?>> stereoElements = (Iterable) mol.stereoElements();
    boolean hasStereo = stereoElements.iterator().hasNext();
    this.tetraChirality = hasStereo ? new int[n] : null;
    this.tetraLigands = hasStereo ? new int[n][] : null;
    this.dbStereoConf = hasStereo && dense ? new int[n][n] : null;
    this.neighbourDbStereo = hasStereo && !dense ? new int[n][] : null;
    if (neighbourDbStereo != null)
      for (int i = 0; i < n; i++) neighbourDbStereo[i] = new int[degreeCount[i]];
    if (hasStereo) {
      for (IStereoElement<?, ?> se : stereoElements) {
        if (se instanceof ITetrahedralChirality) {
          ITetrahedralChirality tc = (ITetrahedralChirality) se;
          Integer idx = idxMap.get(tc.getChiralAtom());
          if (idx != null) {
            // As SMILES and the native core: 1 is anticlockwise (@) and 2 clockwise (@@).
            tetraChirality[idx] = tc.getStereo() == ITetrahedralChirality.Stereo.CLOCKWISE ? 2 : 1;
            IAtom[] ligands = tc.getLigands();
            int[] order = new int[ligands.length];
            for (int k = 0; k < ligands.length; k++)
              order[k] = ligands[k] == tc.getChiralAtom() ? -1 : idxMap.get(ligands[k]);
            tetraLigands[idx] = order;
          }
        } else if (se instanceof IDoubleBondStereochemistry) {
          IDoubleBondStereochemistry dbs = (IDoubleBondStereochemistry) se;
          IBond stereoBond = dbs.getStereoBond();
          Integer a = idxMap.get(stereoBond.getAtom(0)), c = idxMap.get(stereoBond.getAtom(1));
          if (a != null && c != null) {
            int conf = dbs.getStereo() == IDoubleBondStereochemistry.Conformation.TOGETHER ? 1 : 2;
            for (IBond reference : dbs.getBonds()) {
              if (reference.contains(mol.getAtom(a))) {
                int ligand = idxMap.get(reference.getOther(mol.getAtom(a)));
                if (ligand != doubleBondReference(a, c)) conf = 3 - conf;
              } else if (reference.contains(mol.getAtom(c))) {
                int ligand = idxMap.get(reference.getOther(mol.getAtom(c)));
                if (ligand != doubleBondReference(c, a)) conf = 3 - conf;
              }
            }
            if (dbStereoConf != null) {
              dbStereoConf[a][c] = conf;
              dbStereoConf[c][a] = conf;
            } else {
              neighbourDbStereo[a][slot(a, c)] = conf;
              neighbourDbStereo[c][slot(c, a)] = conf;
            }
          }
        }
      }
    }

    this.ringCount = new int[n];
    initDerivedFields();
  }

  MolGraph(Builder b) {
    this.mol = null;
    this.hydrogenCounts = b.hydrogenCounts != null ? b.hydrogenCounts.clone() : null;
    this.n = b.n;
    this.words = CommonBondGraphs.words(n);
    this.name = b.name != null ? b.name : "";
    this.programLine = b.programLine != null ? b.programLine : "";
    this.comment = b.comment != null ? b.comment : "";
    this.properties = b.properties != null ? new LinkedHashMap<>(b.properties) : new LinkedHashMap<>();
    this.atomicNum = b.atomicNum.clone();
    this.formalCharge = b.formalCharge != null ? b.formalCharge.clone() : new int[n];
    this.massNumber = b.massNumber != null ? b.massNumber.clone() : new int[n];
    this.ring = b.ring != null ? b.ring.clone() : new boolean[n];
    this.aromatic = b.aromatic != null ? b.aromatic.clone() : new boolean[n];
    this.label = new int[n];
    this.degree = new int[n];
    for (int i = 0; i < n; i++)
      label[i] = (atomicNum[i] << 2) | (aromatic[i] ? 2 : 0) | (ring[i] ? 1 : 0);

    this.neighbors = new int[n][];
    for (int i = 0; i < n; i++) {
      neighbors[i] = b.neighbors[i] != null ? b.neighbors[i].clone() : new int[0];
      degree[i] = neighbors[i].length;
    }

    // Each row is read from its own atom; the builder has checked that both ends agree.
    boolean dense = n <= SPARSE_THRESHOLD;
    this.bondOrdMatrix = dense ? new int[n][n] : null;
    this.bondRingMatrix = dense ? new boolean[n][n] : null;
    this.bondAromMatrix = dense ? new boolean[n][n] : null;
    this.neighbourBondOrder = dense ? null : new int[n][];
    this.neighbourBondRing = dense ? null : new boolean[n][];
    this.neighbourBondAromatic = dense ? null : new boolean[n][];
    this.dbStereoConf = dense && b.dbStereoConf != null ? new int[n][n] : null;
    this.neighbourDbStereo = !dense && b.dbStereoConf != null ? new int[n][] : null;
    for (int i = 0; i < n; i++) {
      int[] nbs = neighbors[i];
      int[] ords = b.bondOrders != null && b.bondOrders[i] != null ? b.bondOrders[i] : null;
      boolean[] rings = b.bondRings != null && b.bondRings[i] != null ? b.bondRings[i] : null;
      boolean[] aroms = b.bondAroms != null && b.bondAroms[i] != null ? b.bondAroms[i] : null;
      int[] stereo = b.dbStereoConf != null ? b.dbStereoConf[i] : null;
      boolean denseOrd = ords != null && ords.length == n;
      boolean denseRing = rings != null && rings.length == n;
      boolean denseArom = aroms != null && aroms.length == n;
      if (!dense) {
        neighbourBondOrder[i] = new int[nbs.length];
        neighbourBondRing[i] = new boolean[nbs.length];
        neighbourBondAromatic[i] = new boolean[nbs.length];
        if (neighbourDbStereo != null) neighbourDbStereo[i] = new int[nbs.length];
      } else if (dbStereoConf != null && stereo != null) {
        System.arraycopy(stereo, 0, dbStereoConf[i], 0, n);
      }
      for (int k = 0; k < nbs.length; k++) {
        int j = nbs[k];
        int order = ords != null ? (denseOrd ? ords[j] : ords[k]) : 1;
        boolean inRing = rings != null && (denseRing ? rings[j] : rings[k]);
        boolean arom = aroms != null && (denseArom ? aroms[j] : aroms[k]);
        if (dense) {
          bondOrdMatrix[i][j] = order;
          bondRingMatrix[i][j] = inRing;
          bondAromMatrix[i][j] = arom;
        } else {
          neighbourBondOrder[i][k] = order;
          neighbourBondRing[i][k] = inRing;
          neighbourBondAromatic[i][k] = arom;
          if (neighbourDbStereo != null) neighbourDbStereo[i][k] = stereo != null ? stereo[j] : 0;
        }
      }
    }

    this.tetraChirality = b.tetraChirality != null ? b.tetraChirality.clone() : null;
    this.tetraLigands = null;
    this.atomId = b.atomIds != null ? b.atomIds.clone() : null;

    this.ringCount = new int[n];
    initDerivedFields();
  }

  /** The graph of a copy of the molecule, so that later changes to the molecule do not reach it. */
  static MolGraph copyOf(IAtomContainer molecule) {
    Objects.requireNonNull(molecule, "molecule");
    try {
      return new MolGraph(molecule.clone());
    } catch (CloneNotSupportedException failure) {
      throw new IllegalArgumentException("Cannot copy molecule", failure);
    }
  }

  /** Copy matching properties without changing perception or stereo references. */
  static MolGraph copyForSearch(MolGraph source) {
    Objects.requireNonNull(source, "query");
    synchronized (source) {
      return new MolGraph(source);
    }
  }

  private MolGraph(MolGraph source) {
    try {
      mol = source.mol == null ? null : source.mol.clone();
    } catch (CloneNotSupportedException failure) {
      throw new IllegalArgumentException("Cannot copy query molecule", failure);
    }
    hydrogenCounts = source.hydrogenCounts == null ? null : source.hydrogenCounts.clone();
    n = source.n;
    words = source.words;
    atomicNum = source.atomicNum.clone();
    formalCharge = source.formalCharge.clone();
    massNumber = source.massNumber.clone();
    label = source.label.clone();
    ring = source.ring.clone();
    aromatic = source.aromatic.clone();
    degree = source.degree.clone();
    ringCount = source.ringCount.clone();
    ringCountsComputed = source.ringCountsComputed;
    neighbors = copyRows(source.neighbors);
    bondOrdMatrix = copyRows(source.bondOrdMatrix);
    bondRingMatrix = copyRows(source.bondRingMatrix);
    bondAromMatrix = copyRows(source.bondAromMatrix);
    neighbourBondOrder = copyRows(source.neighbourBondOrder);
    neighbourBondRing = copyRows(source.neighbourBondRing);
    neighbourBondAromatic = copyRows(source.neighbourBondAromatic);
    neighbourDbStereo = copyRows(source.neighbourDbStereo);
    tetraChirality = source.tetraChirality == null ? null : source.tetraChirality.clone();
    tetraLigands = copyRows(source.tetraLigands);
    dbStereoConf = copyRows(source.dbStereoConf);
    atomId = source.atomId == null ? null : source.atomId.clone();
    tautomerClass = source.tautomerClass == null ? null : source.tautomerClass.clone();
    tautomerWeight = source.tautomerWeight == null ? null : source.tautomerWeight.clone();
    tautomerSite = source.tautomerSite == null ? null : source.tautomerSite.clone();
    tautomerPKa = source.tautomerPKa == null ? null : source.tautomerPKa.clone();
    adjLong = new long[source.adjLong.length][];
    for (int row = 0; row < adjLong.length; row++) adjLong[row] = source.adjLong[row].clone();
    name = source.name;
    programLine = source.programLine;
    comment = source.comment;
    properties = new LinkedHashMap<>(source.properties);
  }

  private static int[][] copyRows(int[][] rows) {
    if (rows == null) return null;
    int[][] result = new int[rows.length][];
    for (int row = 0; row < rows.length; row++) result[row] = rows[row] == null ? null : rows[row].clone();
    return result;
  }

  private static boolean[][] copyRows(boolean[][] rows) {
    if (rows == null) return null;
    boolean[][] result = new boolean[rows.length][];
    for (int row = 0; row < rows.length; row++) result[row] = rows[row] == null ? null : rows[row].clone();
    return result;
  }

  private void initDerivedFields() {
    // Only adjLong is built eagerly (used by substructure/MCS hot paths).
    this.adjLong = new long[n][words];
    for (int i = 0; i < n; i++)
      for (int k : neighbors[i]) adjLong[i][k >>> 6] |= 1L << (k & 63);
    // Tautomer classes are O(n), cheap, and read as public fields, so they are assigned at once.
    TautomerClasses.assign(this);
    // Canonical labelling and ring counts are lazy (ensureCanonical / ensureRingCounts).
  }

  /** Lazily compute the canonical labelling, orbit partition, automorphisms and canonical hash. */
  void ensureCanonical() {
    if (canonicalComputed) return;
    synchronized (this) {
      if (canonicalComputed) return;
      CanonicalLabels.Labelling labelling = CanonicalLabels.of(this);
      this.canonicalLabel = labelling.labels;
      this.autGeneratorsTruncated = labelling.generatorsTruncated;
      this.autGenerators = labelling.generators;
      this.orbit = labelling.orbits;
      this.canonicalHash = labelling.hash;
      canonicalComputed = true;
    }
  }

  /** Lazily compute per-atom ring membership counts (needed for STRICT ring-fusion matching). */
  public void ensureRingCounts() {
    if (ringCountsComputed) return;
    synchronized (this) {
      if (ringCountsComputed) return;
      computeRingCounts();
      ringCountsComputed = true;
    }
  }

  /** Get ring counts, computing lazily if needed. */
  public int[] getRingCounts() {
    ensureRingCounts();
    return ringCount;
  }

  /**
   * The relevant cycles through each atom, found over each atom's neighbours in index order,
   * as the native core finds them. Unlike one smallest set of smallest rings they do not depend
   * on which equivalent ring a basis keeps, so cages such as cubane count every equivalent ring.
   * Of several equally short paths the candidate search keeps the first it reaches; reading the
   * neighbours in index order makes that choice follow the atom order alone, and not the order
   * in which the bonds were listed. Where a ring system has more relevant cycles of one size
   * than the search reaches, as the macrocycle of a cyclodextrin, the counts still depend on
   * the atom order.
   */
  private void computeRingCounts() {
    java.util.Arrays.fill(ringCount, 0);
    boolean indexOrder = true;
    for (int atom = 0; atom < n && indexOrder; atom++)
      for (int k = 1; k < neighbors[atom].length && indexOrder; k++)
        indexOrder = neighbors[atom][k - 1] < neighbors[atom][k];
    int[][] relevant;
    if (indexOrder) {
      relevant = computeRelevantCycles();
    } else {
      int[][] ordered = new int[n][];
      for (int atom = 0; atom < n; atom++) {
        ordered[atom] = neighbors[atom].clone();
        java.util.Arrays.sort(ordered[atom]);
      }
      relevant = RingPerception.perceive(n, ordered).relevant;
    }
    for (int[] r : relevant) {
      for (int atom : r) ringCount[atom]++;
    }
  }

  /** @deprecated Use {@link #getOrbit()}, the name the native core uses. */
  @Deprecated
  public int[] getOrbits() { return getOrbit(); }
  public int[] getCanonicalLabeling() { ensureCanonical(); return canonicalLabel.clone(); }
  public long getCanonicalHash() { ensureCanonical(); return canonicalHash; }

  /**
   * Return automorphism generators discovered during canonical labelling.
   * Each generator is a permutation array of length n where gen[i] is the
   * image of atom i. Their closure is the full chemical automorphism group only
   * when {@link #automorphismGeneratorsTruncated()} is false. Otherwise they
   * describe a proven subgroup.
   */
  public int[][] getAutomorphismGenerators() { ensureCanonical(); return autGenerators; }

  /** True when generator capture was skipped/capped or chemical-property filtering was incomplete. */
  public boolean automorphismGeneratorsTruncated() { ensureCanonical(); return autGeneratorsTruncated; }

  public String toCanonicalSmiles() {
    return toCanonicalSmiles(Collections.emptyMap());
  }

  /**
   * Write canonical SMILES with positive map labels keyed by original atom indices.
   * Labels affect the output only; canonical ranking and graph identity are unchanged.
   */
  public String toCanonicalSmiles(Map<Integer, Integer> atomMaps) {
    return SmilesWriter.write(this, atomMaps);
  }

  // ---- Builder ----

  public static final class Builder {
    private static final int ROW_ABSENT = 0;
    private static final int ROW_PARALLEL = 1;
    private static final int ROW_DENSE = 2;

    private int n;
    private int[] atomicNum, formalCharge, massNumber, hydrogenCounts;
    private boolean[] ring, aromatic;
    private int[][] neighbors, bondOrders;
    private boolean[][] bondRings, bondAroms;
    private int[] tetraChirality;
    private int[][] dbStereoConf;
    private int[] atomIds;
    private String name;
    private String programLine;
    private String comment;
    private Map<String, String> properties;

    public Builder atomCount(int n) {
      if (n < 0) throw new IllegalArgumentException("atomCount must be >= 0, got " + n);
      this.n = n;
      return this;
    }
    public Builder atomicNumbers(int[] atomicNum) { this.atomicNum = atomicNum; return this; }
    public Builder formalCharges(int[] formalCharge) { this.formalCharge = formalCharge; return this; }
    public Builder massNumbers(int[] massNumber) { this.massNumber = massNumber; return this; }
    /** Optional per-atom hydrogen counts; omitted counts follow the valence model. */
    public Builder hydrogenCounts(int[] hydrogenCounts) { this.hydrogenCounts = hydrogenCounts; return this; }
    public Builder ringFlags(boolean[] ring) { this.ring = ring; return this; }
    public Builder aromaticFlags(boolean[] aromatic) { this.aromatic = aromatic; return this; }
    public Builder neighbors(int[][] neighbors) { this.neighbors = neighbors; return this; }
    public Builder bondOrders(int[][] bondOrders) { this.bondOrders = bondOrders; return this; }
    public Builder bondRingFlags(boolean[][] bondRings) { this.bondRings = bondRings; return this; }
    public Builder bondAromaticFlags(boolean[][] bondAroms) { this.bondAroms = bondAroms; return this; }
    /**
     * Optional tetrahedral tags: 0 for none, 1 for anticlockwise (@) and 2 for clockwise (@@), as
     * in SMILES and the native core. A tag winds the neighbours in their order, with an implicit
     * hydrogen or lone pair second when the first neighbour has a lower index than the atom and
     * first otherwise.
     */
    public Builder tetrahedralChirality(int[] tetraChirality) { this.tetraChirality = tetraChirality; return this; }
    public Builder doubleBondStereo(int[][] dbStereoConf) { this.dbStereoConf = dbStereoConf; return this; }
    /** Optional external atom IDs (1-based, gaps allowed). Length must equal atomCount if provided. */
    public Builder atomIds(int[] ids) { this.atomIds = ids; return this; }
    public Builder name(String name) { this.name = name; return this; }
    public Builder programLine(String programLine) { this.programLine = programLine; return this; }
    public Builder comment(String comment) { this.comment = comment; return this; }
    public Builder properties(Map<String, String> properties) { this.properties = properties; return this; }

    /** Reject an array whose length (-1 when it is missing) is not the atom count. */
    private void requireLength(String name, int length) {
      if (length != n)
        throw new IllegalArgumentException(name + " must have length == atomCount");
    }

    /**
     * The layout of each bond property row: parallel to the neighbour list, or dense over every
     * atom. A missing row is parallel to an empty neighbour list.
     */
    private int[] rowLayouts(Object[] rows, String name) {
      int[] layout = new int[n];
      if (rows == null) return layout;
      for (int i = 0; i < n; i++) {
        Object row = rows[i];
        int length = row == null ? -1 : row instanceof int[] ? ((int[]) row).length : ((boolean[]) row).length;
        int deg = neighbors[i] != null ? neighbors[i].length : 0;
        if (length == deg || (length < 0 && deg == 0)) layout[i] = ROW_PARALLEL;
        else if (length == n) layout[i] = ROW_DENSE;
        else throw new IllegalArgumentException(name + " rows must have length == neighbour count or atomCount");
      }
      return layout;
    }

    private static int intPropAt(int[][] rows, int[] enc, int i, int k, int j, int fallback) {
      if (rows == null) return fallback;
      int[] row = rows[i];
      if (row == null) return fallback;
      return enc[i] == ROW_DENSE ? row[j] : row[k];
    }

    private static boolean boolPropAt(boolean[][] rows, int[] enc, int i, int k, int j, boolean fallback) {
      if (rows == null) return fallback;
      boolean[] row = rows[i];
      if (row == null) return fallback;
      return enc[i] == ROW_DENSE ? row[j] : row[k];
    }

    public MolGraph build() {
      requireLength("atomicNumbers", atomicNum != null ? atomicNum.length : -1);
      requireLength("neighbors", neighbors != null ? neighbors.length : -1);
      // The optional arrays, when supplied, cover every atom.
      if (formalCharge != null) requireLength("formalCharges", formalCharge.length);
      if (massNumber != null) requireLength("massNumbers", massNumber.length);
      if (hydrogenCounts != null) requireLength("hydrogenCounts", hydrogenCounts.length);
      if (ring != null) requireLength("ringFlags", ring.length);
      if (aromatic != null) requireLength("aromaticFlags", aromatic.length);
      if (tetraChirality != null) requireLength("tetrahedralChirality", tetraChirality.length);
      if (atomIds != null) requireLength("atomIds", atomIds.length);
      if (bondOrders != null) requireLength("bondOrders", bondOrders.length);
      if (bondRings != null) requireLength("bondRingFlags", bondRings.length);
      if (bondAroms != null) requireLength("bondAromaticFlags", bondAroms.length);
      if (dbStereoConf != null) requireLength("doubleBondStereo", dbStereoConf.length);

      for (int i = 0; i < n; i++) {
        if (atomicNum[i] < 0 || atomicNum[i] > 118)
          throw new IllegalArgumentException("atomicNumbers must be in [0, 118]");
        if (massNumber != null && massNumber[i] < 0)
          throw new IllegalArgumentException("massNumbers must be >= 0");
        if (hydrogenCounts != null && hydrogenCounts[i] < 0)
          throw new IllegalArgumentException("hydrogenCounts must be >= 0");
        if (tetraChirality != null && tetraChirality[i] != 0 && tetraChirality[i] != 1 && tetraChirality[i] != 2)
          throw new IllegalArgumentException("tetrahedralChirality entries must be 0, 1, or 2");
      }

      @SuppressWarnings("unchecked")
      HashMap<Integer, Integer>[] neighborIndex = new HashMap[n];
      for (int i = 0; i < n; i++) {
        neighborIndex[i] = new HashMap<>();
        if (neighbors[i] == null) continue;
        for (int k = 0; k < neighbors[i].length; k++) {
          int nb = neighbors[i][k];
          if (nb < 0 || nb >= n)
            throw new IllegalArgumentException("neighbour index " + nb + " out of range [0, " + n + ") at atom " + i);
          if (nb == i)
            throw new IllegalArgumentException("self-loop detected in neighbours at atom " + i);
          if (neighborIndex[i].put(nb, k) != null)
            throw new IllegalArgumentException("duplicate neighbour " + nb + " detected at atom " + i);
        }
      }

      int[] bondOrderEnc = rowLayouts(bondOrders, "bondOrders");
      int[] bondRingEnc = rowLayouts(bondRings, "bondRingFlags");
      int[] bondAromEnc = rowLayouts(bondAroms, "bondAromaticFlags");

      for (int i = 0; i < n; i++) {
        if (neighbors[i] == null) continue;
        for (int k = 0; k < neighbors[i].length; k++) {
          int nb = neighbors[i][k];
          Integer rk = neighborIndex[nb].get(i);
          if (rk == null)
            throw new IllegalArgumentException("adjacency must be symmetric");
          if (i > nb) continue;
          if (bondOrders != null) {
            int forward = intPropAt(bondOrders, bondOrderEnc, i, k, nb, 1);
            int backward = intPropAt(bondOrders, bondOrderEnc, nb, rk, i, 1);
            if (forward <= 0 || backward <= 0)
              throw new IllegalArgumentException("bondOrders for present edges must be > 0");
            if (forward != backward)
              throw new IllegalArgumentException("bondOrders must be symmetric");
          }
          if (bondRings != null) {
            boolean forward = boolPropAt(bondRings, bondRingEnc, i, k, nb, false);
            boolean backward = boolPropAt(bondRings, bondRingEnc, nb, rk, i, false);
            if (forward != backward)
              throw new IllegalArgumentException("bondRingFlags must be symmetric");
          }
          if (bondAroms != null) {
            boolean forward = boolPropAt(bondAroms, bondAromEnc, i, k, nb, false);
            boolean backward = boolPropAt(bondAroms, bondAromEnc, nb, rk, i, false);
            if (forward != backward)
              throw new IllegalArgumentException("bondAromaticFlags must be symmetric");
          }
        }
      }

      if (dbStereoConf != null) {
        for (int i = 0; i < n; i++) {
          if (dbStereoConf[i] != null && dbStereoConf[i].length != n)
            throw new IllegalArgumentException("doubleBondStereo rows must have length == atomCount");
        }
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            int conf = dbStereoConf[i] != null ? dbStereoConf[i][j] : 0;
            if (conf < 0 || conf > 2)
              throw new IllegalArgumentException("doubleBondStereo entries must be 0, 1, or 2");
            if (i == j && conf != 0)
              throw new IllegalArgumentException("doubleBondStereo diagonal must be zero");
            if (j <= i) continue;
            int rev = dbStereoConf[j] != null ? dbStereoConf[j][i] : 0;
            if (conf != rev)
              throw new IllegalArgumentException("doubleBondStereo must be symmetric");
            if (conf != 0) {
              Integer k = neighborIndex[i].get(j);
              if (k == null)
                throw new IllegalArgumentException("doubleBondStereo requires a bond");
              int ord = intPropAt(bondOrders, bondOrderEnc, i, k, j, 1);
              if (ord != 2 && ord != 4)
                throw new IllegalArgumentException("doubleBondStereo requires bond order 2 or 4");
            }
          }
        }
      }

      if (atomIds != null) {
        HashSet<Integer> seen = new HashSet<>();
        for (int id : atomIds) {
          if (id <= 0)
            throw new IllegalArgumentException("atomIds must be positive");
          if (!seen.add(id))
            throw new IllegalArgumentException("atomIds must be unique");
        }
      }
      return new MolGraph(this);
    }
  }

  // ---- NLF helpers (package-accessible) ----

  static final int[] EMPTY_NLF = new int[0];

  /**
   * NLF label: atomicNum + aromaticity, WITHOUT ring bit.
   * Ring compatibility is enforced separately in atom/bond checks so the
   * cached NLF tables remain reusable across option profiles. Matches C++ nlfLabel().
   */
  static int nlfLabel(MolGraph g, int idx) {
    return (g.atomicNum[idx] << 1) | (g.aromatic[idx] ? 1 : 0);
  }

  static int[] buildNLF1(MolGraph g, int idx) {
    int[] nbs = g.neighbors[idx];
    int deg = nbs.length;
    if (deg == 0) return EMPTY_NLF;
    int[] labels = new int[deg];
    for (int i = 0; i < deg; i++) labels[i] = nlfLabel(g, nbs[i]);
    return labelCounts(labels, deg);
  }

  /** Sort the first {@code count} labels and return ascending (label, frequency) pairs. */
  private static int[] labelCounts(int[] labels, int count) {
    if (count == 0) return EMPTY_NLF;
    Arrays.sort(labels, 0, count);
    int distinct = 1;
    for (int i = 1; i < count; i++) if (labels[i] != labels[i - 1]) distinct++;
    int[] arr = new int[distinct * 2];
    int k = 0;
    arr[0] = labels[0]; arr[1] = 1;
    for (int i = 1; i < count; i++) {
      if (labels[i] == labels[i - 1]) { arr[k + 1]++; }
      else { k += 2; arr[k] = labels[i]; arr[k + 1] = 1; }
    }
    return arr;
  }

  static int[] buildNLF2(MolGraph g, int idx) {
    return buildNeighborhoodNLF(g, idx, 2);
  }

  static int[] buildNLF3(MolGraph g, int idx) {
    return buildNeighborhoodNLF(g, idx, 3);
  }

  /**
   * Count atoms within the radius, excluding the root. Exact distance shells
   * are unsafe for non-induced matching: extra target edges can shorten a
   * mapped atom's distance without invalidating the query embedding.
   */
  private static int[] buildNeighborhoodNLF(MolGraph g, int idx, int radius) {
    return buildNeighborhoodNLF(g, idx, radius, new int[g.n], 1, new int[g.n]);
  }

  /**
   * The same labels with the caller's work arrays of at least {@code g.n} entries each: an atom
   * counts as seen when its {@code seen} entry equals {@code stamp}, which no entry holds before.
   */
  static int[] buildNeighborhoodNLF(MolGraph g, int idx, int radius, int[] seen, int stamp, int[] reached) {
    seen[idx] = stamp;
    reached[0] = idx;
    int levelStart = 0, end = 1;
    for (int distance = 0; distance < radius && levelStart < end; distance++) {
      int levelEnd = end;
      for (int k = levelStart; k < levelEnd; k++)
        for (int j : g.neighbors[reached[k]])
          if (seen[j] != stamp) { seen[j] = stamp; reached[end++] = j; }
      levelStart = levelEnd;
    }
    // Reuse the visit list for the labels of every reached atom except the root.
    for (int k = 1; k < end; k++) reached[k - 1] = nlfLabel(g, reached[k]);
    return labelCounts(reached, end - 1);
  }

  // ---------------------------------------------------------------------------
  // Lenient SMILES pre-sanitiser (parity with C++ ParseOptions::lenient)
  // ---------------------------------------------------------------------------

  /**
   * Best-effort sanitisation of a SMILES string before CDK parsing.
   * <ul>
   *   <li>Balances parentheses: appends missing {@code )} or strips trailing extra {@code )}</li>
   *   <li>Removes unclosed ring-digit openings (single digits and {@code %nn} notation)</li>
   * </ul>
   *
   * @param smi the raw SMILES string
   * @return a cleaned SMILES that CDK's SmilesParser is more likely to accept
   */
  public static String sanitizeSmiles(String smi) {
    if (smi == null || smi.isEmpty()) return smi;

    // --- Pass 1: balance parentheses ---
    int depth = 0;
    for (int i = 0; i < smi.length(); i++) {
      char c = smi.charAt(i);
      if (c == '[') { // skip bracket atoms entirely
        while (i < smi.length() - 1 && smi.charAt(i) != ']') i++;
      } else if (c == '(') {
        depth++;
      } else if (c == ')') {
        if (depth > 0) depth--;
        // else: will be stripped below
      }
    }

    StringBuilder sb = new StringBuilder(smi.length() + depth);
    int excess = 0; // tracks extra ')' that have no matching '('
    int openCount = 0;
    for (int i = 0; i < smi.length(); i++) {
      char c = smi.charAt(i);
      if (c == '[') {
        int start = i;
        while (i < smi.length() - 1 && smi.charAt(i) != ']') i++;
        sb.append(smi, start, i + 1);
        continue;
      }
      if (c == '(') {
        openCount++;
        sb.append(c);
      } else if (c == ')') {
        if (openCount > 0) {
          openCount--;
          sb.append(c);
        }
        // else: skip extra closing paren
      } else {
        sb.append(c);
      }
    }
    // Append missing closing parens
    for (int i = 0; i < openCount; i++) sb.append(')');

    String balanced = sb.toString();

    // --- Pass 2: remove unclosed ring digits ---
    // Track which ring digits are opened but never closed
    Map<Integer, Integer> ringFirstPos = new LinkedHashMap<>();
    Set<Integer> closedRings = new HashSet<>();
    for (int i = 0; i < balanced.length(); i++) {
      char c = balanced.charAt(i);
      if (c == '[') {
        while (i < balanced.length() - 1 && balanced.charAt(i) != ']') i++;
        continue;
      }
      int ringNum = -1;
      int consumed = 0;
      if (c == '%' && i + 2 < balanced.length()
          && Character.isDigit(balanced.charAt(i + 1))
          && Character.isDigit(balanced.charAt(i + 2))) {
        ringNum = (balanced.charAt(i + 1) - '0') * 10 + (balanced.charAt(i + 2) - '0');
        consumed = 3;
      } else if (Character.isDigit(c) && i > 0) {
        // Only treat as ring digit if preceded by an atom character
        char prev = balanced.charAt(i - 1);
        if (prev == ']' || Character.isLetter(prev) || Character.isDigit(prev)) {
          ringNum = c - '0';
          consumed = 1;
        }
      }
      if (ringNum >= 0) {
        if (ringFirstPos.containsKey(ringNum) && !closedRings.contains(ringNum)) {
          closedRings.add(ringNum);
          ringFirstPos.remove(ringNum);
        } else {
          ringFirstPos.put(ringNum, i);
        }
        if (consumed > 1) i += consumed - 1;
      }
    }

    if (ringFirstPos.isEmpty()) return balanced;

    // Rebuild string, stripping unclosed ring references
    Set<Integer> removePositions = new HashSet<>(ringFirstPos.values());
    // For %nn notation, also remove the two following digits
    Map<Integer, Integer> posToLen = new HashMap<>();
    for (Map.Entry<Integer, Integer> entry : ringFirstPos.entrySet()) {
      int pos = entry.getValue();
      if (pos < balanced.length() && balanced.charAt(pos) == '%') {
        posToLen.put(pos, 3);
      } else {
        posToLen.put(pos, 1);
      }
    }

    StringBuilder result = new StringBuilder(balanced.length());
    for (int i = 0; i < balanced.length(); i++) {
      if (posToLen.containsKey(i)) {
        i += posToLen.get(i) - 1; // skip
      } else {
        result.append(balanced.charAt(i));
      }
    }
    return result.toString();
  }

  static boolean nlfOk(int[] fq, int[] ft) {
    int fi = 0, ti = 0;
    while (fi < fq.length) {
      int fLabel = fq[fi], fFreq = fq[fi + 1];
      while (ti < ft.length && ft[ti] < fLabel) ti += 2;
      if (ti >= ft.length || ft[ti] != fLabel || ft[ti + 1] < fFreq) return false;
      fi += 2;
    }
    return true;
  }

  @FunctionalInterface
  interface NLFBuilder {
    int[] build(MolGraph g, int idx);
  }

  static int[][] buildAllNLF(MolGraph g, NLFBuilder builder) {
    int[][] nlf = new int[g.n][];
    for (int i = 0; i < g.n; i++) nlf[i] = builder.build(g, i);
    return nlf;
  }

  // ---- Atom and bond compatibility (ChemOps) ----

  /** The atom and bond matching rules that every search applies. */
  static final class ChemOps {
    /**
     * Whether two query atoms are indistinguishable to {@link #atomsCompatFast}. This reads every
     * query property that method reads under the same options, so equal atoms have equal domains.
     */
    static boolean sameAtomConstraints(MolGraph gq, int first, int second, ChemOptions C) {
      if (C.matchAtomType && gq.atomicNum[first] != gq.atomicNum[second]) return false;
      if (C.matchFormalCharge && gq.formalCharge[first] != gq.formalCharge[second]) return false;
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT && gq.aromatic[first] != gq.aromatic[second])
        return false;
      if ((C.ringMatchesRingOnly || C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE)
          && gq.ring[first] != gq.ring[second]) return false;
      if (C.matchIsotope && gq.massNumber[first] != gq.massNumber[second]) return false;
      if (C.useChirality && gq.tetraChirality != null) {
        int stereo = gq.tetraChirality[first];
        if (stereo != gq.tetraChirality[second]) return false;
        if (stereo != 0 && (gq.tetrahedralLigands(first) == null) != (gq.tetrahedralLigands(second) == null))
          return false;
      }
      return C.ringFusionMode != ChemOptions.RingFusionMode.STRICT || gq.ringCount[first] == gq.ringCount[second];
    }

    static boolean atomsCompatFast(MolGraph gq, int qi, MolGraph gt, int tj, ChemOptions C) {
      // Tautomerism changes hydrogen placement and bond order, not elemental identity.
      if (C.matchAtomType && gq.atomicNum[qi] != gt.atomicNum[tj]) return false;
      if (C.matchFormalCharge && gq.formalCharge[qi] != gt.formalCharge[tj]) return false;
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT && gq.aromatic[qi] != gt.aromatic[tj]) return false;
      if (C.ringMatchesRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
      if (C.matchIsotope) {
        int qm = gq.massNumber[qi], tm = gt.massNumber[tj];
        if (qm != 0 && tm != 0 && qm != tm) return false;
      }
      if (C.useChirality && gq.tetraChirality != null && gt.tetraChirality != null) {
        int qs = gq.tetraChirality[qi], ts = gt.tetraChirality[tj];
        if (qs != 0 && ts != 0) {
          // Preserve the comparison of explicitly supplied invalid builder annotations.
          if ((gq.tetrahedralLigands(qi) == null || gt.tetrahedralLigands(tj) == null) && qs != ts) return false;
        }
      }
      if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE && gq.ring[qi] && gt.ring[tj]) {
        if (C.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
          if (gq.ringCount[qi] != gt.ringCount[tj]) return false;
        }
      }
      return true;
    }

    static boolean bondsCompatible(MolGraph g1, int qi, int qk, MolGraph g2, int tj, int tk, ChemOptions C) {
      int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
      if (qOrd == 0 || tOrd == 0) return false;
      // Tautomer-aware: relax bond ORDER (single/double interchangeable)
      // but still enforce aromaticity, stereo, and ring constraints.
      boolean tautBondRelax = C.tautomerAware && g1.tautomerClass != null && g2.tautomerClass != null
          && g1.tautomerClass[qi] != -1 && g1.tautomerClass[qk] != -1
          && g2.tautomerClass[tj] != -1 && g2.tautomerClass[tk] != -1;
      // Strict aromaticity: aromatic bonds must match aromatic bonds (no ring guard needed)
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT
          && g1.bondAromatic(qi, qk) != g2.bondAromatic(tj, tk)) return false;
      if (C.ringMatchesRingOnly && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk)) return false;
      if (tautBondRelax) return true;  // bond order relaxed for tautomeric bonds
      if (C.matchBondOrder == ChemOptions.BondOrderMode.ANY) return true;
      if (qOrd == tOrd) return true;
      if (C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE) return true;
      if (C.aromaticityMode == ChemOptions.AromaticityMode.FLEXIBLE) {
        boolean qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
        if ((qa && ta) || (qa && (tOrd == 1 || tOrd == 2)) || (ta && (qOrd == 1 || qOrd == 2))) return true;
        // Resonance: a single and a double bond are interchangeable when both ends of each are
        // aromatic, which covers Kekule forms drawn differently (such as a pyridone drawn as C=O
        // or as C-O) without tautomer classes.
        if (g1.aromatic[qi] && g1.aromatic[qk] && g2.aromatic[tj] && g2.aromatic[tk]
            && (qOrd == 1 || qOrd == 2) && (tOrd == 1 || tOrd == 2)) return true;
      }
      return false;
    }
  }
}

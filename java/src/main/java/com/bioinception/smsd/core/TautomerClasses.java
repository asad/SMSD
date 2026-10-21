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
import java.util.function.IntConsumer;

/**
 * Tautomer classes of a graph and their relevance weights: atoms that a proton shift can join
 * share a class, and each atom of a class takes the weight of its motif, from the {@code TW_}
 * weights of {@link MolGraph}, with the ionisable site and {@code TPKA_} pKa of that motif when
 * the library assigns it one and a match of the motif that covers the atom holds the hydrogen of
 * the site.
 *
 * <p>Every motif pattern is evaluated over all of its matches. The atoms of one match share a
 * class, and matches that share an atom join their classes, so that atoms which interconvert
 * through any chain of tautomeric motifs form one class, in the manner of the mobile-hydrogen
 * groups of InChI but without requiring a mobile hydrogen. The partition is a function of the
 * graph alone, and the classes are numbered by an invariant of the graph (see
 * {@link #classNumbers}), so the partition does not depend on the order of the atoms and the
 * numbers depend on it only among classes that Weisfeiler-Leman refinement cannot tell apart. The
 * motif of an atom is the first, in a fixed priority order, whose pattern covers the atom, so the
 * weight and site do not depend on that order either. Atoms outside every class keep class -1,
 * weight 1 and no site. The native core uses the same patterns, sites, order and numbering.
 */
final class TautomerClasses {
  private TautomerClasses() {}

  /**
   * The Henderson-Hasselbalch neutral fraction of a class with an acid site of pKa
   * {@code acidPKa} when {@code hasAcid} and a base site of pKa {@code basePKa} (that of the
   * conjugate acid BH+) when {@code hasBase}: 1 / (1 + 10^(pH - pKa_a) + 10^(pKa_b - pH)), without
   * the term of a kind of site the class lacks. The relation takes the dominant site of each kind
   * and leaves out the zwitterion, so it holds for pKa_b below pKa_a. The fraction is computed in
   * double precision with {@link StrictMath#pow}, the acid term added before the base term, and
   * rounded to single precision, as the native core does, so that the two agree within the shared
   * golden tolerance of 1e-12.
   */
  static float neutralFraction(boolean hasAcid, float acidPKa, boolean hasBase, float basePKa, double pH) {
    double denominator = 1.0;
    if (hasAcid) denominator += StrictMath.pow(10.0, pH - (double) acidPKa);
    if (hasBase) denominator += StrictMath.pow(10.0, (double) basePKa - pH);
    return (float) (1.0 / denominator);
  }

  /** A tautomer weight scaled by the neutral fraction of its class, clamped to [0.1, 0.99]. */
  static float neutralFormWeight(float weight, float fraction) {
    return Math.max(0.1f, Math.min(0.99f, weight * fraction));
  }

  /**
   * The neutral fraction of each atom's class at {@code pH}, for the atoms that have an ionisable
   * site; the class takes the lowest acid pKa and the highest base pKa among the sites of its
   * atoms. Atoms without a site, or outside every class, get NaN and keep their weight at every pH.
   */
  static float[] neutralFractions(MolGraph graph, double pH) {
    int n = graph.n;
    float[] fractions = new float[n];
    Arrays.fill(fractions, Float.NaN);
    if (graph.tautomerClass == null || graph.tautomerSite == null || graph.tautomerPKa == null
        || graph.tautomerClass.length < n || graph.tautomerSite.length < n || graph.tautomerPKa.length < n)
      return fractions;
    int classes = 0;
    for (int atom = 0; atom < n; atom++) classes = Math.max(classes, graph.tautomerClass[atom] + 1);
    boolean[] hasAcid = new boolean[classes], hasBase = new boolean[classes];
    float[] acidPKa = new float[classes], basePKa = new float[classes];
    for (int atom = 0; atom < n; atom++) {
      int cls = graph.tautomerClass[atom];
      if (cls < 0) continue;
      float pKa = graph.tautomerPKa[atom];
      if (graph.tautomerSite[atom] == MolGraph.TAUTOMER_SITE_ACID) {
        if (!hasAcid[cls] || pKa < acidPKa[cls]) acidPKa[cls] = pKa;
        hasAcid[cls] = true;
      } else if (graph.tautomerSite[atom] == MolGraph.TAUTOMER_SITE_BASE) {
        if (!hasBase[cls] || pKa > basePKa[cls]) basePKa[cls] = pKa;
        hasBase[cls] = true;
      }
    }
    for (int atom = 0; atom < n; atom++) {
      int cls = graph.tautomerClass[atom];
      if (cls < 0 || graph.tautomerSite[atom] == MolGraph.TAUTOMER_SITE_NONE) continue;
      fractions[atom] = neutralFraction(hasAcid[cls], acidPKa[cls], hasBase[cls], basePKa[cls], pH);
    }
    return fractions;
  }

  /** Assigns the tautomer classes, weights and ionisable sites of {@code graph}. */
  static void assign(MolGraph graph) {
    Motifs motifs = motifs(graph);
    // Every atom a motif covers is in a class, so the weight, site and pKa of each atom are those
    // of its motif; atoms outside every class keep weight 1 and no site.
    graph.tautomerClass = classNumbers(graph, motifs);
    graph.tautomerWeight = motifs.weight;
    graph.tautomerSite = motifs.site;
    graph.tautomerPKa = motifs.pKa;
  }

  /**
   * The aqueous N-H pKa of pyrimidin-4(3H)-one, which a pyridone/hydroxypyridine or pyridinone
   * match takes in place of the 2-pyridone value when the other ring neighbour of its nitrogen is a
   * carbon bonded to a second aromatic ring nitrogen.
   */
  private static final float TPKA_PYRIMIDINONE = 8.6f;

  /** The weight, ionisable site and pKa that the motif priority gives each atom, and the classes. */
  static final class Motifs {
    /** 1 when no motif covers the atom. */
    final float[] weight;
    /** One of the {@code MolGraph.TAUTOMER_SITE_} values. */
    final byte[] site;
    /** 0 when no match of the atom's motif carries a site. */
    final float[] pKa;
    /** Whether a motif covers the atom, that is whether it has a class. */
    final boolean[] covered;
    /**
     * A union-find over the atoms in which the atoms of every motif match share a set, so the sets
     * of the covered atoms are the tautomer classes; read them with {@link UnionFind#root}.
     */
    final int[] joined;

    private Motifs(int n) {
      weight = new float[n];
      site = new byte[n];
      pKa = new float[n];
      covered = new boolean[n];
      joined = UnionFind.create(n);
      Arrays.fill(weight, 1.0f);
    }
  }

  /**
   * The motif and the class set of every atom.
   *
   * <p>Each atom takes the weight of the first motif, in the fixed priority order of
   * {@link MotifPass#run}, whose pattern covers it, and the ionisable site and pKa that the matches
   * of that motif covering it carry. Every motif claims the atoms of all of its matches before the
   * next motif is tried, so the choice does not depend on the order of the atoms or of their
   * neighbours, and the atoms that some motif covers are exactly the atoms of the classes. A match
   * of a motif with an acid site carries the site only when an atom of the match that would release
   * the proton bears a hydrogen: the N-O oxygen of a hydroxamic acid, so that an N-alkoxy amide has
   * none; the nitrogen of an imide, urea, amide, thioamide or thioamide/iminothiol; the bridging
   * carbon of a 1,3-diketone; the alpha carbon of a ketone or thione; the ring nitrogen or the
   * oxygen of a pyridone/hydroxypyridine or pyridinone; and a phenol oxygen with no second
   * heavy-atom neighbour, so that anisoles, aryl esters and lactam or lactone carbonyls have none.
   * When every such atom with a hydrogen belongs to an earlier motif, the match carries the lowest
   * acid site among them instead, so that each hydrogen keeps the pKa of its own motif; an atom
   * covered by several matches of its motif takes the lowest acid pKa among them. No motif carries
   * a base site. The order puts the most specific motif first: hydroxamic acid, imide
   * and urea before amide; 1,3-diketone, enaminone and 1,5-keto-enol before keto/enol and amide;
   * the aromatic lactam motifs (pyridone/hydroxypyridine, pyridinone, hydroxypyrimidine) before
   * amide and phenol; quinone methide before phenol; tetrazole and triazole before pyrazole;
   * benzimidazole before imidazole; guanidine before amidine. The imide motif comes before the
   * barbituric motif, whose ring nitrogen between two ring carbonyls is a cyclic imide, so that
   * the nitrogen and the carbonyl carbons keep the imide N-H site; the barbituric motif claims the
   * carbonyl oxygens. Of the motifs with the same pattern, the pyridone/hydroxypyridine motif,
   * which carries the acid site, comes before the hydroxypyrimidine motif, and the
   * furanose/pyranose motif before the lactol motif. The native core uses the same order.
   */
  static Motifs motifs(MolGraph graph) {
    MotifPass pass = new MotifPass(graph);
    pass.run();
    return pass.motifs;
  }

  /** One pass over the motifs in priority order. */
  private static final class MotifPass {
    private final MolGraph g;
    private final int n;
    private final int[][] neighbors;
    private final int[] atomicNum, formalCharge;
    private final boolean[] ring, aromatic;
    private final Motifs motifs;
    /** The motif being tried, and the motif that claimed each atom, the first whose pattern covers it. */
    private int motifIndex = -1;
    private final int[] owner;
    private float motifWeight = 1.0f;
    /** The site that the match being claimed carries, set by {@link #carry} before its atoms are claimed. */
    private byte matchSite = MolGraph.TAUTOMER_SITE_NONE;
    private float matchPKa = 0.0f;
    private boolean ownSite, earlierSite;
    private float earlierPKa;

    MotifPass(MolGraph graph) {
      g = graph;
      n = graph.n;
      neighbors = graph.neighbors;
      atomicNum = graph.atomicNum;
      formalCharge = graph.formalCharge;
      ring = graph.ring;
      aromatic = graph.aromatic;
      motifs = new Motifs(n);
      owner = new int[n];
      Arrays.fill(owner, -1);
    }

    private void claim(int atom) {
      final byte ACID = MolGraph.TAUTOMER_SITE_ACID;
      if (!motifs.covered[atom]) {
        motifs.covered[atom] = true;
        owner[atom] = motifIndex;
        motifs.weight[atom] = motifWeight;
        motifs.site[atom] = matchSite;
        motifs.pKa[atom] = matchPKa;
      } else if (owner[atom] == motifIndex && matchSite == ACID
          && (motifs.site[atom] != ACID || matchPKa < motifs.pKa[atom])) {
        // Of the matches of its motif that cover an atom, the one with the lowest acid pKa sets
        // its site, whatever the order in which they are found.
        motifs.site[atom] = ACID;
        motifs.pKa[atom] = matchPKa;
      }
    }

    /**
     * Starts the site of the next match: {@link #siteAtom} for each atom of the match that would
     * release the proton, then {@link #carry}.
     */
    private void beginSite() {
      ownSite = earlierSite = false;
      earlierPKa = 0.0f;
    }

    private void siteAtom(int atom) {
      if (g.attachedHydrogenCount(atom) == 0) return;
      if (!motifs.covered[atom] || owner[atom] == motifIndex) {
        ownSite = true;
      } else if (motifs.site[atom] == MolGraph.TAUTOMER_SITE_ACID && (!earlierSite || motifs.pKa[atom] < earlierPKa)) {
        earlierSite = true;
        earlierPKa = motifs.pKa[atom];
      }
    }

    /**
     * The match carries the acid site of pKa {@code pKa} when an atom given to {@link #siteAtom}
     * with a hydrogen belongs to this motif (or to none yet, so that the match claims it), else the
     * lowest acid site of the earlier motifs of those atoms, else no site.
     */
    private void carry(float pKa) {
      matchSite = ownSite || earlierSite ? MolGraph.TAUTOMER_SITE_ACID : MolGraph.TAUTOMER_SITE_NONE;
      matchPKa = ownSite ? pKa : earlierSite ? earlierPKa : 0.0f;
    }

    private void carries(float pKa, int... atoms) {
      beginSite();
      for (int atom : atoms) siteAtom(atom);
      carry(pKa);
    }

    /** Adds {@code atom} to the match of {@code with}. */
    private void join(int atom, int with) {
      claim(atom);
      UnionFind.join(motifs.joined, atom, with);
    }

    /** Claims the atoms of one match. */
    private void group(int... atoms) {
      claim(atoms[0]);
      for (int atom : atoms) join(atom, atoms[0]);
    }

    /** A motif without a site claims its matches with none; a motif with one carries each match first. */
    private void motif(float weight, IntConsumer match) {
      motifIndex++;
      motifWeight = weight;
      matchSite = MolGraph.TAUTOMER_SITE_NONE;
      matchPKa = 0.0f;
      for (int i = 0; i < n; i++) match.accept(i);
    }

    private int bond(int a, int b) { return g.bondOrder(a, b); }

    private boolean hasNitrogenBesides(int carbon, int except) {
      for (int k : neighbors[carbon]) if (k != except && atomicNum[k] == 7) return true;
      return false;
    }

    /** The N-O of a hydroxamic acid. */
    private boolean isHydroxyl(int nitrogen, int o) {
      return atomicNum[o] == 8 && !aromatic[o] && bond(nitrogen, o) == 1;
    }

    /** A C=O carbon, whose oxygen is not {@code except}. */
    private boolean isCarbonyl(int carbon, int except) {
      for (int k : neighbors[carbon]) if (k != except && atomicNum[k] == 8 && bond(carbon, k) == 2) return true;
      return false;
    }

    private boolean aromaticRingN(int atom) { return atomicNum[atom] == 7 && aromatic[atom] && ring[atom]; }

    /**
     * The acid pKa of a pyridone/hydroxypyridine or pyridinone match of ring nitrogen
     * {@code nitrogen} and oxygen-bearing carbon {@code carbon}: that of pyrimidin-4(3H)-one when the
     * other ring neighbour of the nitrogen is a carbon bonded to a second aromatic ring nitrogen,
     * else {@code pKa}.
     */
    private float lactamPKa(int nitrogen, int carbon, float pKa) {
      for (int m : neighbors[nitrogen]) {
        if (m == carbon || atomicNum[m] != 6 || !aromatic[m] || !ring[m]) continue;
        for (int p : neighbors[m]) if (p != nitrogen && aromaticRingN(p)) return TPKA_PYRIMIDINONE;
      }
      return pKa;
    }

    private int heavyNeighbours(int atom) {
      int count = 0;
      for (int k : neighbors[atom]) if (atomicNum[k] > 1) count++;
      return count;
    }

    /** The para carbons c of an aromatic ring carbon j (path j-a-b-c), as T4 and T24 walk them. */
    private void forParaCarbons(int i, int j, IntConsumer visit) {
      for (int a : neighbors[j]) {
        if (a == i || atomicNum[a] != 6 || !aromatic[a] || !ring[a]) continue;
        for (int b : neighbors[a]) {
          if (b == j || atomicNum[b] != 6 || !aromatic[b] || !ring[b]) continue;
          for (int c : neighbors[b])
            if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) visit.accept(c);
        }
      }
    }

    /** T1: the site is the N-H of an amide and the alpha C-H of a ketone. */
    private IntConsumer carbonyl(boolean amide, float pKa) {
      return i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2 || hasNitrogenBesides(j, i) != amide) continue;
          for (int k : neighbors[j])
            if (k != i && ((atomicNum[k] == 6 && !aromatic[k]) || atomicNum[k] == 7)) {
              beginSite();
              if (!amide || atomicNum[k] == 7) siteAtom(k);
              carry(pKa);
              group(i, k);
            }
        }
      };
    }

    /** T2: the site is the N-H of a thioamide and the non-aromatic alpha C-H of a thione. */
    private IntConsumer thiocarbonyl(boolean amide, float pKa) {
      return i -> {
        if (atomicNum[i] != 16 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2 || hasNitrogenBesides(j, i) != amide) continue;
          for (int k : neighbors[j])
            if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) {
              beginSite();
              if (amide ? atomicNum[k] == 7 : !aromatic[k]) siteAtom(k);
              carry(pKa);
              group(i, k);
            }
        }
      };
    }

    void run() {
      // T13 hydroxamic acid N(-OH)-C=O: the N, its N-O oxygens and every C=O on the N; the site is
      // the O-H of an N-O oxygen
      motif(MolGraph.TW_HYDROXAMIC, i -> {
        if (atomicNum[i] != 7) return;
        boolean hydroxyl = false, carbonyl = false;
        for (int o : neighbors[i]) hydroxyl = hydroxyl || isHydroxyl(i, o);
        for (int j : neighbors[i]) carbonyl = carbonyl || (atomicNum[j] == 6 && isCarbonyl(j, i));
        if (!hydroxyl || !carbonyl) return;
        beginSite();
        for (int o : neighbors[i]) if (isHydroxyl(i, o)) siteAtom(o);
        carry(MolGraph.TPKA_HYDROXAMIC);
        claim(i);
        for (int o : neighbors[i]) if (isHydroxyl(i, o)) join(o, i);
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 8 && bond(j, k) == 2) { join(j, i); join(k, i); }
        }
      });
      // T12 imide C(=O)-N-C(=O): the N and its carbonyl carbons; the site is the N-H
      motif(MolGraph.TW_IMIDE, i -> {
        if (atomicNum[i] != 7) return;
        int carbonyls = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && isCarbonyl(j, i)) carbonyls++;
        if (carbonyls < 2) return;
        carries(MolGraph.TPKA_IMIDE, i);
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && isCarbonyl(j, i)) join(j, i);
      });
      // T28 barbituric: ring N between ring carbonyls, with their oxygens
      motif(MolGraph.TW_BARBITURIC, i -> {
        if (atomicNum[i] != 7 || !ring[i]) return;
        int carbonyls = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 6 && ring[j] && isCarbonyl(j, i)) carbonyls++;
        if (carbonyls < 2) return;
        claim(i);
        for (int c : neighbors[i]) {
          if (atomicNum[c] != 6 || !ring[c] || !isCarbonyl(c, i)) continue;
          join(c, i);
          for (int k : neighbors[c]) if (atomicNum[k] == 8 && bond(c, k) == 2) join(k, i);
        }
      });
      // T10 urea O=C(-N)(-N): the O, its carbon and that carbon's nitrogens; the site is an N-H
      motif(MolGraph.TW_UREA, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          int nitrogens = 0;
          for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) nitrogens++;
          if (nitrogens < 2) continue;
          beginSite();
          for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) siteAtom(k);
          carry(MolGraph.TPKA_UREA);
          group(i, j);
          for (int k : neighbors[j]) if (k != i && atomicNum[k] == 7) join(k, i);
        }
      });
      // T5 1,3-diketone O=C-C-C=O: both oxygens and the bridging carbon, whose C-H is the site
      motif(MolGraph.TW_DIKETONE_ENOL, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          for (int bridge : neighbors[j]) {
            if (bridge == i || atomicNum[bridge] != 6 || aromatic[bridge]) continue;
            for (int m : neighbors[bridge]) {
              if (m == j || m == i || atomicNum[m] != 6) continue;
              boolean directlyBonded = false;
              for (int nb : neighbors[j]) if (nb == m) { directlyBonded = true; break; }
              if (directlyBonded) continue;
              for (int o2 : neighbors[m])
                if (o2 != bridge && atomicNum[o2] == 8 && !aromatic[o2] && bond(m, o2) == 2) {
                  carries(MolGraph.TPKA_DIKETONE_ENOL, bridge);
                  group(i, bridge, o2);
                }
            }
          }
        }
      });
      // T11 beta-enaminone N-C=C-C=O
      motif(MolGraph.TW_ENAMINONE, i -> {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 1) continue;
          for (int k : neighbors[j]) {
            if (k == i || atomicNum[k] != 6 || bond(j, k) != 2) continue;
            for (int m : neighbors[k]) {
              if (m == j || atomicNum[m] != 6) continue;
              for (int o : neighbors[m])
                if (o != k && atomicNum[o] == 8 && bond(m, o) == 2) group(i, j, k, m, o);
            }
          }
        }
      });
      // T21 1,5-keto-enol O=C-C=C-C
      motif(MolGraph.TW_15_KETO_ENOL, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          for (int k : neighbors[j]) {
            if (k == i || atomicNum[k] != 6 || bond(j, k) != 1) continue;
            for (int m : neighbors[k]) {
              if (m == j || atomicNum[m] != 6 || bond(k, m) != 2) continue;
              for (int p : neighbors[m])
                if (p != k && atomicNum[p] == 6 && !aromatic[p]) group(i, j, k, m, p);
            }
          }
        }
      });
      // T27 pyridone/hydroxypyridine: aromatic ring N, aromatic ring C, single-bonded O; the site
      // is the N-H or O-H
      motif(MolGraph.TW_PYRIDONE_HYDROXYPYRIDINE, i -> {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 8 && !aromatic[k] && bond(j, k) == 1) {
              carries(lactamPKa(i, j, MolGraph.TPKA_PYRIDONE_HYDROXYPYRIDINE), i, k);
              group(i, j, k);
            }
        }
      });
      // T6 pyridinone: aromatic ring N next to an aromatic C bearing an exocyclic O; the site is the
      // N-H or O-H
      motif(MolGraph.TW_PYRIDINONE, i -> {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j]) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 8 && !aromatic[k]) {
              carries(lactamPKa(i, j, MolGraph.TPKA_PYRIDINONE), i, k);
              group(i, k);
            }
        }
      });
      // T14 hydroxypyrimidine: exocyclic O on an aromatic ring C next to an aromatic ring N
      motif(MolGraph.TW_HYDROXYPYRIMIDINE, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j] || bond(i, j) != 1) continue;
          for (int k : neighbors[j])
            if (k != i && aromaticRingN(k)) group(i, j, k);
        }
      });
      // T1 keto/enol C=O with its alpha carbons, and its nitrogens for an amide
      motif(MolGraph.TW_AMIDE_IMIDIC, carbonyl(true, MolGraph.TPKA_AMIDE_IMIDIC));
      motif(MolGraph.TW_KETO_ENOL, carbonyl(false, MolGraph.TPKA_KETO_ENOL));
      // T2 thione/thiol C=S with the carbon and nitrogen neighbours of its carbon
      motif(MolGraph.TW_THIOAMIDE, thiocarbonyl(true, MolGraph.TPKA_THIOAMIDE));
      motif(MolGraph.TW_THIONE_THIOL, thiocarbonyl(false, MolGraph.TPKA_THIONE_THIOL));
      // T16 thioamide/iminothiol S=C-N; the site is the N-H
      motif(MolGraph.TW_THIOAMIDE_IMINOTHIOL, i -> {
        if (atomicNum[i] != 16 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 7) {
              carries(MolGraph.TPKA_THIOAMIDE_IMINOTHIOL, k);
              group(i, j, k);
            }
        }
      });
      // T3 nitroso/oxime N=C-O
      motif(MolGraph.TW_NITROSO_OXIME, i -> {
        if (atomicNum[i] != 7 || ring[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 8) group(i, k);
        }
      });
      // T24 phenol/quinone methide: para carbon with an exocyclic C=C; no site
      motif(MolGraph.TW_PHENOL_QUINONE_METHIDE, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j] || bond(i, j) != 1) continue;
          forParaCarbons(i, j, c -> {
            for (int d : neighbors[c])
              if (atomicNum[d] == 6 && !ring[d] && bond(c, d) == 2) group(i, j, c, d);
          });
        }
      });
      // T4 phenol/quinone: O on an aromatic ring C, its para carbon and that carbon's oxygens; the
      // site is the O-H of an oxygen with no second heavy-atom neighbour
      motif(MolGraph.TW_PHENOL_QUINONE, i -> {
        if (atomicNum[i] != 8 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
          forParaCarbons(i, j, c -> {
            beginSite();
            if (heavyNeighbours(i) == 1) siteAtom(i);
            carry(MolGraph.TPKA_PHENOL_QUINONE);
            group(i, j, c);
            for (int d : neighbors[c]) if (atomicNum[d] == 8 && !aromatic[d]) join(d, i);
          });
        }
      });
      // T15 tetrazole/triazole: aromatic ring N with two or more aromatic ring N neighbours
      motif(MolGraph.TW_TETRAZOLE_NH, i -> {
        if (!aromaticRingN(i)) return;
        int count = 0;
        for (int j : neighbors[i]) if (aromaticRingN(j)) count++;
        if (count < 2) return;
        claim(i);
        for (int j : neighbors[i]) if (aromaticRingN(j)) join(j, i);
      });
      // T25 1,2,3-triazole NH shift N-N-N
      motif(MolGraph.TW_TRIAZOLE_NH_SHIFT, i -> {
        if (!aromaticRingN(i)) return;
        int only = -1, count = 0;
        for (int j : neighbors[i]) if (aromaticRingN(j)) { only = j; count++; }
        if (count != 1) return;
        for (int k : neighbors[only])
          if (k != i && aromaticRingN(k)) group(i, only, k);
      });
      // T26 benzimidazole N-C-N through a fusion carbon
      motif(MolGraph.TW_BENZIMIDAZOLE_NH, i -> {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
          int ringNeighbours = 0;
          for (int k : neighbors[j]) if (ring[k] && aromatic[k]) ringNeighbours++;
          if (ringNeighbours < 3) continue;
          for (int k : neighbors[j])
            if (k != i && aromaticRingN(k)) group(i, j, k);
        }
      });
      // T7 pyrazole N-N, then imidazole N-C-N in a five-membered ring
      motif(MolGraph.TW_IMIDAZOLE_NH, i -> {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) if (aromaticRingN(j)) group(i, j);
      });
      motif(MolGraph.TW_IMIDAZOLE_NH, i -> {
        if (!aromaticRingN(i)) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !aromatic[j] || !ring[j]) continue;
          for (int k : neighbors[j]) {
            if (k == i || !aromaticRingN(k)) continue;
            boolean in5 = false;
            for (int[] r : g.computeRings()) {
              if (r.length != 5) continue;
              boolean hasI = false, hasK = false;
              for (int a : r) { if (a == i) hasI = true; if (a == k) hasK = true; }
              if (hasI && hasK) { in5 = true; break; }
            }
            if (in5) group(i, k);
          }
        }
      });
      // T9 guanidine: C with three or more non-aromatic N
      motif(MolGraph.TW_GUANIDINE, i -> {
        if (atomicNum[i] != 6) return;
        int count = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 7 && !aromatic[j]) count++;
        if (count < 3) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 7 && !aromatic[j]) join(j, i);
      });
      // T8 amidine N=C-N
      motif(MolGraph.TW_AMIDINE, i -> {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 7 && !aromatic[k]) group(i, j, k);
        }
      });
      // T17 phosphonate P with two or more O
      motif(MolGraph.TW_PHOSPHONATE, i -> {
        if (atomicNum[i] != 15) return;
        int count = 0;
        for (int j : neighbors[i]) if (atomicNum[j] == 8) count++;
        if (count < 2) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 8) join(j, i);
      });
      // T19 nitro/aci-nitro: non-aromatic N with two or more O and a C
      motif(MolGraph.TW_NITRO_ACI, i -> {
        if (atomicNum[i] != 7 || aromatic[i]) return;
        int oxygens = 0, carbons = 0;
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 8) oxygens++;
          else if (atomicNum[j] == 6) carbons++;
        }
        if (oxygens < 2 || carbons == 0) return;
        claim(i);
        for (int j : neighbors[i]) if (atomicNum[j] == 8 || atomicNum[j] == 6) join(j, i);
      });
      // T22 furanose/pyranose: ring O, ring C, exocyclic single-bonded O
      motif(MolGraph.TW_FURANOSE_PYRANOSE, i -> {
        if (atomicNum[i] != 8 || !ring[i] || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || !ring[j]) continue;
          for (int k : neighbors[j])
            if (k != i && atomicNum[k] == 8 && !ring[k] && bond(j, k) == 1) group(i, j, k);
        }
      });
      // T23 lactol: ring C bonded to a ring O and an exocyclic single-bonded O
      motif(MolGraph.TW_LACTOL, i -> {
        if (atomicNum[i] != 6 || !ring[i]) return;
        boolean hasRingO = false, hasExoO = false;
        for (int j : neighbors[i]) {
          hasRingO = hasRingO || (atomicNum[j] == 8 && ring[j] && !aromatic[j]);
          hasExoO = hasExoO || (atomicNum[j] == 8 && !ring[j] && !aromatic[j] && bond(i, j) == 1);
        }
        if (!hasRingO || !hasExoO) return;
        claim(i);
        for (int j : neighbors[i])
          if ((atomicNum[j] == 8 && ring[j] && !aromatic[j])
              || (atomicNum[j] == 8 && !ring[j] && !aromatic[j] && bond(i, j) == 1)) join(j, i);
      });
      // T29 allyl shift X(-H)-C=C-C, X = O, S or N
      motif(MolGraph.TW_ALLYL_SHIFT, i -> {
        if ((atomicNum[i] != 8 && atomicNum[i] != 16 && atomicNum[i] != 7) || aromatic[i] || ring[i]) return;
        int sumBO = 0;
        for (int nb : neighbors[i]) sumBO += bond(i, nb);
        int maxValence = (atomicNum[i] == 7) ? 3 : 2;
        if (!(sumBO < maxValence + Math.abs(formalCharge[i]))) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || aromatic[j] || bond(i, j) != 1) continue;
          for (int k : neighbors[j]) {
            if (k == i || atomicNum[k] != 6 || bond(j, k) != 2) continue;
            for (int m : neighbors[k])
              if (m != j && atomicNum[m] == 6 && !aromatic[m] && bond(k, m) == 1) group(i, j, k, m);
          }
        }
      });
      // T30 selenol/selenoketone C=Se with the carbon and nitrogen neighbours of its carbon
      motif(MolGraph.TW_SELENOL_SELENOKETONE, i -> {
        if (atomicNum[i] != 34 || aromatic[i]) return;
        for (int j : neighbors[i]) {
          if (atomicNum[j] != 6 || bond(i, j) != 2) continue;
          group(i, j);
          for (int k : neighbors[j]) if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) join(k, i);
        }
      });
    }
  }

  // ---------------------------------------------------------------------------
  // Class numbering, step for step as the native core numbers the classes
  // ---------------------------------------------------------------------------

  private static final long COLOUR_MASK = 0x3FFFFFFFFFFFFFFFL;

  /** The 62-bit hash of {@code value} combined into {@code hash}, shared bit for bit with the native core. */
  static long colourMix(long hash, long value) {
    long x = (hash ^ (value + 0x9E3779B97F4A7C15L)) * 0xFF51AFD7ED558CCDL;
    x ^= x >>> 33;
    x *= 0xC4CEB9FE1A85EC53L;
    x ^= x >>> 33;
    return x & COLOUR_MASK;
  }

  private static int distinctColours(long[] colours) {
    long[] sorted = colours.clone();
    Arrays.sort(sorted);
    int distinct = 0;
    for (int index = 0; index < sorted.length; index++)
      if (index == 0 || sorted[index] != sorted[index - 1]) distinct++;
    return distinct;
  }

  /** Lexicographic order of sorted colour lists, a proper prefix first; colours are below 2^62. */
  private static int compareKeys(long[] first, long[] second) {
    int common = Math.min(first.length, second.length);
    for (int index = 0; index < common; index++)
      if (first[index] != second[index]) return Long.compare(first[index], second[index]);
    return Integer.compare(first.length, second.length);
  }

  /**
   * The tautomer class of every atom: -1 for an atom that no motif covers, otherwise the number of
   * its set in {@code motifs.joined}.
   *
   * <p>The classes are numbered in increasing order of a key, the sorted colours of their atoms.
   * The colour of an atom starts from its element, charge, mass number, hydrogen count,
   * aromaticity, ring membership, degree and whether a motif covers it, and each round of
   * refinement hashes into it the bond orders and colours of its neighbours, as in
   * Weisfeiler-Leman refinement; rounds stop once the keys differ or the colours no longer split.
   * When classes still share the smallest key the refinement cannot tell them apart (in practice a
   * symmetry of the graph exchanges them): the one with the lowest atom index takes the next
   * number, its atoms are given colours of their own and the refinement resumes, so that the
   * classes tied to it are told apart consistently. The numbering is therefore the same for every
   * order of the atoms, except among classes that the refinement cannot tell apart. Colours are
   * 64-bit hashes cut to 62 bits, so that signed and unsigned comparisons agree, and every step is
   * that of the native core.
   */
  static int[] classNumbers(MolGraph g, Motifs motifs) {
    int n = g.n;
    int[] numbers = new int[n];
    Arrays.fill(numbers, -1);
    // The atoms of each class, in increasing atom order.
    int[] classOfRoot = new int[n];
    Arrays.fill(classOfRoot, -1);
    List<List<Integer>> classAtoms = new ArrayList<>();
    for (int atom = 0; atom < n; atom++) {
      if (!motifs.covered[atom]) continue;
      int root = UnionFind.root(motifs.joined, atom);
      if (classOfRoot[root] < 0) {
        classOfRoot[root] = classAtoms.size();
        classAtoms.add(new ArrayList<>());
      }
      classAtoms.get(classOfRoot[root]).add(atom);
    }
    int classes = classAtoms.size();
    if (classes == 1)
      for (int atom : classAtoms.get(0)) numbers[atom] = 0;
    if (classes <= 1) return numbers;
    int[][] members = new int[classes][];
    for (int cls = 0; cls < classes; cls++) {
      List<Integer> atoms = classAtoms.get(cls);
      members[cls] = new int[atoms.size()];
      for (int index = 0; index < members[cls].length; index++) members[cls][index] = atoms.get(index);
    }

    long[] colours = new long[n];
    for (int atom = 0; atom < n; atom++) {
      long colour = colourMix(0x5441555443L, g.atomicNum[atom]);
      colour = colourMix(colour, g.formalCharge[atom]);
      colour = colourMix(colour, g.massNumber[atom]);
      colour = colourMix(colour, g.hydrogenCount(atom));
      colour = colourMix(colour, g.aromatic[atom] ? 1 : 0);
      colour = colourMix(colour, g.ring[atom] ? 1 : 0);
      colour = colourMix(colour, g.neighbors[atom].length);
      colours[atom] = colourMix(colour, motifs.covered[atom] ? 1 : 0);
    }
    Numbering numbering = new Numbering(g, members, colours, distinctColours(colours));
    int next = 0;
    while (numbering.remaining.size() > 0) {
      while (numbering.sortRemaining() && numbering.refine()) {}
      int position = 0;
      boolean individualised = false;
      List<Integer> remaining = numbering.remaining;
      for (; position < remaining.size() && !individualised; position++) {
        int cls = remaining.get(position);
        int number = next++;
        for (int atom : members[cls]) numbers[atom] = number;
        if (position + 1 < remaining.size()
            && compareKeys(numbering.keys[cls], numbering.keys[remaining.get(position + 1)]) == 0) {
          for (int atom : members[cls])
            numbering.colours[atom] = colourMix(colourMix(0x494E444956L, number), numbering.colours[atom]);
          numbering.distinct = distinctColours(numbering.colours);
          individualised = true;
        }
      }
      remaining.subList(0, position).clear();
    }
    return numbers;
  }

  /** The colours, keys and unnumbered classes of one numbering. */
  private static final class Numbering {
    final MolGraph g;
    final int[][] members;
    long[] colours;
    int distinct;
    final long[][] keys;
    final List<Integer> remaining = new ArrayList<>();

    Numbering(MolGraph g, int[][] members, long[] colours, int distinct) {
      this.g = g;
      this.members = members;
      this.colours = colours;
      this.distinct = distinct;
      keys = new long[members.length][];
      for (int cls = 0; cls < members.length; cls++) remaining.add(cls);
    }

    /** One round of refinement; false when it splits no colour, which leaves the colours as they were. */
    boolean refine() {
      int n = g.n;
      long[] refined = new long[n];
      for (int atom = 0; atom < n; atom++) {
        int[] adjacent = g.neighbors[atom];
        long[] around = new long[adjacent.length];
        for (int index = 0; index < adjacent.length; index++)
          around[index] = colourMix(colourMix(0x424F4E44L, g.bondOrder(atom, adjacent[index])), colours[adjacent[index]]);
        Arrays.sort(around);
        long colour = colourMix(0x524F554E44L, colours[atom]);
        for (long value : around) colour = colourMix(colour, value);
        refined[atom] = colour;
      }
      int split = distinctColours(refined);
      if (split <= distinct) return false;
      colours = refined;
      distinct = split;
      return true;
    }

    /** Sorts the unnumbered classes by key and then by lowest atom; true when two share a key. */
    boolean sortRemaining() {
      for (int cls : remaining) {
        long[] key = new long[members[cls].length];
        for (int index = 0; index < key.length; index++) key[index] = colours[members[cls][index]];
        Arrays.sort(key);
        keys[cls] = key;
      }
      remaining.sort((first, second) -> {
        int order = compareKeys(keys[first], keys[second]);
        return order != 0 ? order : Integer.compare(members[first][0], members[second][0]);
      });
      boolean tied = false;
      for (int index = 1; index < remaining.size(); index++)
        tied = tied || compareKeys(keys[remaining.get(index - 1)], keys[remaining.get(index)]) == 0;
      return tied;
    }
  }
}

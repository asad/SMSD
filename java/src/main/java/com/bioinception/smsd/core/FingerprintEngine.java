/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Fingerprint computation engine for SMSD: path, MCS-aware, and circular (ECFP/FCFP) fingerprints.
 *
 * <p>All fingerprints are returned as {@code long[]} bitsets. Use the format conversion methods
 * ({@link #toBitSet}, {@link #toHex}, {@link #toBinaryString}) for interop with databases and REST APIs.
 *
 * <p>For batch workloads, pre-convert molecules to {@link MolGraph} once and use the MolGraph overloads
 * to avoid repeated IAtomContainer→MolGraph conversion overhead.
 *
 * @author Syed Asad Rahman
 */
public final class FingerprintEngine {

  /** FNV-1a 64-bit seed (offset basis). */
  static final long FNV1A_SEED = -3750763034362895579L;
  /** FNV-1a 64-bit prime. */
  static final long FNV1A_PRIME = 1099511628211L;

  private FingerprintEngine() {} // static utility

  /**
   * Does nothing: each graph keeps its own default fingerprints, which are released with the
   * graph.
   *
   * @deprecated There is no shared cache to clear; this method does nothing.
   */
  @Deprecated
  public static void clearFPCache() {}

  /** Fingerprint invariant mode for circular fingerprints. */
  public enum Mode {
    /** ECFP: atomic number, degree, implicit H count, ring, aromaticity, charge, tautomer class. */
    ECFP,
    /** FCFP: pharmacophoric feature classes (H-bond donor/acceptor, ionisable, aromatic, hydrophobic). */
    FCFP
  }

  // ---- MolGraph conversion ----

  /**
   * Standardise an {@link IAtomContainer} and convert it to a {@link MolGraph} for batch use.
   *
   * <p>When computing multiple fingerprints for the same molecule (e.g. ECFP + FCFP + counts),
   * convert once and pass the {@link MolGraph} to each fingerprint method to avoid redundant
   * IAtomContainer&rarr;MolGraph conversion overhead.
   *
   * <pre>{@code
   * MolGraph g = FingerprintEngine.toMolGraph(mol);
   * long[] ecfp = FingerprintEngine.ecfp(g, 2, 2048);
   * long[] fcfp = FingerprintEngine.fcfp(g, 2, 2048);
   * int[] ecfpC = FingerprintEngine.ecfpCounts(g, 2, 2048);
   * }</pre>
   *
   * @param mol the CDK molecule to convert
   * @return a MolGraph suitable for all fingerprint MolGraph overloads
   */
  public static MolGraph toMolGraph(IAtomContainer mol) {
    if (mol == null) throw new NullPointerException("mol");
    return new MolGraph(standardiseSafe(mol));
  }

  // ---- Path Fingerprint ----

  /**
   * Compute a path-based fingerprint with automatic standardisation.
   *
   * <p>The molecule is standardised (aromaticity, typing, implicit H) before fingerprinting
   * so that SMILES with different explicit-H representations or kekulised forms produce
   * identical fingerprints.
   *
   * @param mol        the CDK molecule to fingerprint
   * @param pathLength maximum path length (number of bonds) to enumerate
   * @param fpSize     fingerprint size in bits (should be a power of 2)
   * @return the fingerprint as a long array (bit-packed)
   */
  public static long[] pathFingerprint(IAtomContainer mol, int pathLength, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (pathLength < 0 || fpSize < 1) throw new IllegalArgumentException("Invalid fingerprint dimensions");
    mol = standardiseSafe(mol);
    int words = (fpSize - 1) / 64 + 1;
    long[] fp = new long[words];
    int n = mol.getAtomCount();
    if (n == 0) return fp;
    // Build index map, skipping explicit H atoms
    IdentityHashMap<IAtom, Integer> idxMap = new IdentityHashMap<>(n);
    int heavyCount = 0;
    int[] origToHeavy = new int[n];
    for (int i = 0; i < n; i++) {
      Integer az = mol.getAtom(i).getAtomicNumber();
      if (az != null && az == 1) { origToHeavy[i] = -1; continue; } // skip H
      origToHeavy[i] = heavyCount;
      idxMap.put(mol.getAtom(i), heavyCount++);
    }
    if (heavyCount == 0) return fp;
    int maxLength = Math.min(pathLength, heavyCount - 1);
    int[][] adj = new int[heavyCount][];
    int[] atomicNum = new int[heavyCount];
    for (int i = 0; i < n; i++) {
      if (origToHeavy[i] < 0) continue;
      int hi = origToHeavy[i];
      Integer az = mol.getAtom(i).getAtomicNumber();
      atomicNum[hi] = az != null ? az : 0;
      List<IAtom> nbAtoms = mol.getConnectedAtomsList(mol.getAtom(i));
      int cnt = 0;
      for (IAtom nb : nbAtoms) { Integer nbIdx = idxMap.get(nb); if (nbIdx != null) cnt++; }
      adj[hi] = new int[cnt];
      int k = 0;
      for (IAtom nb : nbAtoms) { Integer nbIdx = idxMap.get(nb); if (nbIdx != null) adj[hi][k++] = nbIdx; }
    }
    walkPaths(adj, null, maxLength, (path, length) -> fpHashPath(atomicNum, path, length), fp, fpSize);
    return fp;
  }

  // ---- MCS Fingerprint ----

  public static long[] mcsFingerprint(IAtomContainer mol, ChemOptions chem, int pathLength, int fpSize) {
    MolGraph g = FingerprintEngine.toMolGraph(mol);
    return mcsFingerprint(g, chem != null && chem.tautomerAware, pathLength, fpSize);
  }

  // ---- Subset & Similarity ----

  public static boolean fingerprintSubset(long[] query, long[] target) {
    int len = Math.min(query.length, target.length);
    for (int i = 0; i < len; i++) {
      if ((query[i] & ~target[i]) != 0) return false;
    }
    for (int i = len; i < query.length; i++) {
      if (query[i] != 0) return false;
    }
    return true;
  }

  /**
   * Set-bit counts of two fingerprints: bits in both, bits in the first and bits in the second.
   * Words missing from the shorter fingerprint are empty. The metrics below mirror
   * {@code smsd::batch} in {@code batch.hpp}, including their empty-fingerprint rules.
   */
  private static long[] bitCounts(long[] fp1, long[] fp2) {
    long common = 0, first = 0, second = 0;
    int shared = Math.min(fp1.length, fp2.length);
    for (int i = 0; i < shared; i++) {
      common += Long.bitCount(fp1[i] & fp2[i]);
      first += Long.bitCount(fp1[i]);
      second += Long.bitCount(fp2[i]);
    }
    for (int i = shared; i < fp1.length; i++) first += Long.bitCount(fp1[i]);
    for (int i = shared; i < fp2.length; i++) second += Long.bitCount(fp2[i]);
    return new long[] {common, first, second};
  }

  /**
   * Tanimoto coefficient: |A intersection B| / |A union B|.
   *
   * @param fp1 first fingerprint (long[] bitset)
   * @param fp2 second fingerprint (long[] bitset)
   * @return Tanimoto similarity in [0.0, 1.0]; two empty fingerprints are identical sets and
   *     score 1.0; returns 0.0 for null inputs
   */
  public static double tanimoto(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null) return 0.0;
    long[] counts = bitCounts(fp1, fp2);
    long either = counts[1] + counts[2] - counts[0];
    return either == 0 ? 1.0 : (double) counts[0] / either;
  }

  /**
   * Dice similarity coefficient: 2*|A intersection B| / (|A| + |B|).
   *
   * <p>Ranges from 0.0 (disjoint) to 1.0 (identical). Dice weights shared features
   * more heavily than Tanimoto, making it better suited for count-based fingerprints
   * and scaffold-hopping similarity searches.
   *
   * @param fp1 first fingerprint (long[] bitset)
   * @param fp2 second fingerprint (long[] bitset)
   * @return Dice similarity in [0.0, 1.0]; returns 0.0 for null inputs and when both
   *     fingerprints are empty
   */
  public static double dice(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null) return 0.0;
    long[] counts = bitCounts(fp1, fp2);
    long sum = counts[1] + counts[2];
    return sum == 0 ? 0.0 : (2.0 * counts[0]) / sum;
  }

  /**
   * Cosine similarity: |A intersection B| / sqrt(|A| * |B|).
   *
   * <p>Ranges from 0.0 (orthogonal) to 1.0 (identical). Cosine similarity normalises
   * for fingerprint density, making it useful when comparing fingerprints of
   * molecules with very different sizes.
   *
   * @param fp1 first fingerprint (long[] bitset)
   * @param fp2 second fingerprint (long[] bitset)
   * @return cosine similarity in [0.0, 1.0]; returns 0.0 for null inputs and when either
   *     fingerprint is empty
   */
  public static double cosine(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null) return 0.0;
    long[] counts = bitCounts(fp1, fp2);
    double denom = Math.sqrt((double) counts[1] * counts[2]);
    return (denom == 0.0) ? 0.0 : counts[0] / denom;
  }

  /**
   * Soergel distance: 1 - Tanimoto.
   *
   * <p>A proper metric distance in [0.0, 1.0] derived from the Tanimoto coefficient.
   * Useful for clustering and distance-based algorithms (e.g. k-nearest neighbours,
   * Taylor-Butina clustering).
   *
   * @param fp1 first fingerprint (long[] bitset)
   * @param fp2 second fingerprint (long[] bitset)
   * @return Soergel distance in [0.0, 1.0]; two empty fingerprints are at distance 0.0;
   *     returns 1.0 for null inputs
   */
  public static double soergel(long[] fp1, long[] fp2) {
    return 1.0 - tanimoto(fp1, fp2);
  }

  /**
   * Overlap coefficient: |A ∩ B| / min(|A|, |B|).
   *
   * <p>Measures the overlap between two fingerprints relative to the smaller set.
   * Returns 1.0 when the smaller fingerprint is a subset of the larger one.
   * Useful for scaffold/sub-structure containment checks.
   *
   * @param fp1 first fingerprint (long[] bitset)
   * @param fp2 second fingerprint (long[] bitset)
   * @return overlap coefficient in [0.0, 1.0]; two empty fingerprints score 1.0 and one
   *     empty fingerprint scores 0.0
   */
  public static double overlapCoefficient(long[] fp1, long[] fp2) {
    long[] counts = bitCounts(fp1, fp2);
    if (counts[1] == 0 && counts[2] == 0) return 1.0;
    long smaller = Math.min(counts[1], counts[2]);
    return smaller == 0 ? 0.0 : (double) counts[0] / smaller;
  }

  /** Density above which {@link #analyzeFpQuality(long[], int)} reports a fingerprint as saturated. */
  private static final double SATURATION_DENSITY = 0.75;

  /**
   * Analyse fingerprint quality: density, expected collision rate, saturation, as the C++
   * {@code smsd::batch::fingerprintQuality} and the Python {@code analyze_fp_quality} do.
   *
   * <p>Hashing n features uniformly into N bits sets N(1 - (1 - 1/N)^n) bits on average, so s set
   * bits suggest n = ln(1 - s/N) / ln(1 - 1/N) features, of which a share 1 - s/n, the
   * {@code expected_collision_rate}, landed on a bit that an earlier feature had set; it is 0 for
   * an empty fingerprint and 1 for a full one. Path fingerprints of a power-of-two size collide
   * more often than uniform hashing predicts, so the rate understates their collisions. A
   * fingerprint is saturated above a density of 0.75.
   *
   * @param fp   fingerprint (long[] bitset)
   * @param fpSize nominal fingerprint size in bits
   * @return ordered map with keys: set_bits, fp_size, density, expected_collision_rate, is_saturated
   */
  public static Map<String, Object> analyzeFpQuality(long[] fp, int fpSize) {
    int setBits = nominalBitCount(fp, fpSize);
    double density = (double) setBits / fpSize;
    double collisionRate = 0.0;
    if (setBits == fpSize) {
      collisionRate = 1.0;
    } else if (setBits > 0) {
      double features = Math.log1p(-density) / Math.log1p(-1.0 / fpSize);
      collisionRate = Math.min(1.0, Math.max(0.0, 1.0 - setBits / features));
    }
    Map<String, Object> result = new LinkedHashMap<>();
    result.put("set_bits", setBits);
    result.put("fp_size", fpSize);
    result.put("density", density);
    result.put("expected_collision_rate", collisionRate);
    result.put("is_saturated", density > SATURATION_DENSITY);
    return result;
  }

  // ---- Count-vector similarity metrics ----

  /**
   * Count-vector Tanimoto: sum(min(a,b)) / sum(max(a,b)).
   *
   * <p>Generalisation of the Tanimoto coefficient for integer count vectors
   * (e.g., ECFP count fingerprints). Equivalent to the binary Tanimoto when
   * all counts are 0 or 1.
   *
   * @param fp1 first count vector
   * @param fp2 second count vector (must have the same length as fp1)
   * @return count Tanimoto in [0.0, 1.0]; returns 0.0 for null/empty/mismatched inputs
   */
  public static double tanimotoCounts(int[] fp1, int[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp1.length != fp2.length) return 0.0;
    long minSum = 0, maxSum = 0;
    for (int i = 0; i < fp1.length; i++) {
      minSum += Math.min(fp1[i], fp2[i]);
      maxSum += Math.max(fp1[i], fp2[i]);
    }
    return (maxSum == 0) ? 0.0 : (double) minSum / maxSum;
  }

  /**
   * Count-vector Dice: 2 * sum(min(a,b)) / (sum(a) + sum(b)).
   *
   * @param fp1 first count vector
   * @param fp2 second count vector (must have the same length as fp1)
   * @return count Dice in [0.0, 1.0]; returns 0.0 for null/empty/mismatched inputs
   */
  public static double diceCounts(int[] fp1, int[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp1.length != fp2.length) return 0.0;
    long minSum = 0, aSum = 0, bSum = 0;
    for (int i = 0; i < fp1.length; i++) {
      minSum += Math.min(fp1[i], fp2[i]);
      aSum   += fp1[i];
      bSum   += fp2[i];
    }
    long denom = aSum + bSum;
    return (denom == 0) ? 0.0 : (2.0 * minSum) / denom;
  }

  /**
   * Count-vector Cosine: dot(a,b) / (|a| * |b|).
   *
   * @param fp1 first count vector
   * @param fp2 second count vector (must have the same length as fp1)
   * @return count cosine in [0.0, 1.0]; returns 0.0 for null/empty/mismatched inputs
   */
  public static double cosineCounts(int[] fp1, int[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp1.length != fp2.length) return 0.0;
    long dot = 0, aSqSum = 0, bSqSum = 0;
    for (int i = 0; i < fp1.length; i++) {
      dot    += (long) fp1[i] * fp2[i];
      aSqSum += (long) fp1[i] * fp1[i];
      bSqSum += (long) fp2[i] * fp2[i];
    }
    double denom = Math.sqrt((double) aSqSum * bSqSum);
    return (denom == 0.0) ? 0.0 : dot / denom;
  }

  // ---- Circular ECFP ----

  /**
   * Compute an ECFP circular fingerprint from a CDK molecule.
   *
   * @deprecated Use {@link #ecfp(IAtomContainer, int, int)}, which this method calls, or
   *     {@link #circularFingerprint(IAtomContainer, int, int, Mode)} to choose the invariants.
   */
  @Deprecated
  public static long[] circularFingerprint(IAtomContainer mol, int radius, int fpSize) {
    return ecfp(mol, radius, fpSize);
  }

  /**
   * Compute an ECFP circular fingerprint from a CDK molecule.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #ecfp(MolGraph, int, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the CDK molecule
   * @param radius maximum expansion radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   */
  public static long[] ecfp(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return ecfp(toMolGraph(mol), radius, fpSize);
  }

  public static long[] ecfp(MolGraph g, int radius, int fpSize) {
    checkCircular(g, radius, fpSize);
    // The default fingerprint is asked for repeatedly, so the graph keeps it.
    boolean defaults = radius == 2 && fpSize == 2048;
    long[] cached = defaults ? g.cachedEcfpDefault : null;
    if (cached != null) return cached.clone();
    long[] fp = new long[(fpSize + 63) / 64];
    if (g.n > 0) circular(g, ecfpAtomHashes(g), radius, fpSize, binaryStrategy(fp));
    if (defaults) g.cachedEcfpDefault = fp.clone();
    return fp;
  }

  private static void checkCircular(MolGraph g, int radius, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    if (radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
  }

  /** Record the atom hashes, then the hashes of every radius up to {@code radius}. */
  private static void circular(MolGraph g, long[] atomHash, int radius, int fpSize, BitStrategy strategy) {
    for (long h : atomHash) strategy.record(h, fpSize);
    expandRadiiImpl(g, atomHash, radius, fpSize, g.n, strategy);
  }

  /**
   * The ECFP invariant of each atom: the Rogers and Hahn properties (element, heavy degree, bond
   * order sum, mass number, charge, attached hydrogens and ring membership), aromaticity and the
   * tautomer class, hashed with FNV-1a.
   */
  private static long[] ecfpAtomHashes(MolGraph g) {
    long[] atomHash = new long[g.n];
    for (int i = 0; i < g.n; i++) {
      long h = FNV1A_SEED;
      h ^= g.atomicNum[i];       h *= FNV1A_PRIME;
      h ^= g.degree[i];          h *= FNV1A_PRIME;
      int bondOrdSum = 0;
      for (int nb : g.neighbors[i]) bondOrdSum += g.bondOrder(i, nb);
      h ^= bondOrdSum;            h *= FNV1A_PRIME;
      h ^= g.massNumber[i];       h *= FNV1A_PRIME;
      h ^= (long) g.formalCharge[i] + 4L; h *= FNV1A_PRIME;
      h ^= g.attachedHydrogenCount(i); h *= FNV1A_PRIME;
      h ^= g.ring[i] ? 1 : 0;     h *= FNV1A_PRIME;
      h ^= g.aromatic[i] ? 1 : 0; h *= FNV1A_PRIME;
      if (g.tautomerClass != null && g.tautomerClass[i] >= 0) {
        h ^= 0xCAFE00L + g.tautomerClass[i]; h *= FNV1A_PRIME;
      }
      atomHash[i] = h;
    }
    return atomHash;
  }

  /** The FCFP invariant of each atom: its pharmacophore features, degree and ring membership. */
  private static long[] fcfpAtomHashes(MolGraph g) {
    int[] pharmFeatures = g.getPharmacophoreFeatures();
    long[] atomHash = new long[g.n];
    for (int i = 0; i < g.n; i++) {
      long h = FNV1A_SEED;
      h ^= pharmFeatures[i];  h *= FNV1A_PRIME;
      h ^= g.degree[i];       h *= FNV1A_PRIME;
      h ^= g.ring[i] ? 1 : 0; h *= FNV1A_PRIME;
      atomHash[i] = h;
    }
    return atomHash;
  }

  // ---- Circular FCFP ----

  /**
   * Compute an FCFP circular fingerprint from a CDK molecule.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #fcfp(MolGraph, int, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the CDK molecule
   * @param radius maximum expansion radius (2 = FCFP4, 3 = FCFP6; -1 = unlimited)
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   */
  public static long[] fcfp(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return fcfp(toMolGraph(mol), radius, fpSize);
  }

  public static long[] fcfp(MolGraph g, int radius, int fpSize) {
    checkCircular(g, radius, fpSize);
    boolean defaults = radius == 2 && fpSize == 2048;
    long[] cached = defaults ? g.cachedFcfpDefault : null;
    if (cached != null) return cached.clone();
    long[] fp = new long[(fpSize + 63) / 64];
    if (g.n > 0) circular(g, fcfpAtomHashes(g), radius, fpSize, binaryStrategy(fp));
    if (defaults) g.cachedFcfpDefault = fp.clone();
    return fp;
  }

  // ---- Mode dispatch ----

  /**
   * Compute a circular fingerprint (ECFP or FCFP) from a CDK molecule.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #circularFingerprint(MolGraph, int, int, Mode)} to avoid repeated
   * conversion overhead.
   *
   * @param mol    the CDK molecule
   * @param radius maximum expansion radius
   * @param fpSize fingerprint size in bits
   * @param mode   ECFP or FCFP invariant mode
   * @return binary fingerprint as long[] bitset
   */
  public static long[] circularFingerprint(IAtomContainer mol, int radius, int fpSize, Mode mode) {
    if (mol == null) throw new NullPointerException("mol");
    MolGraph g = toMolGraph(mol);
    return circularFingerprint(g, radius, fpSize, mode);
  }

  public static long[] circularFingerprint(MolGraph g, int radius, int fpSize, Mode mode) {
    return (mode == Mode.FCFP) ? fcfp(g, radius, fpSize) : ecfp(g, radius, fpSize);
  }

  // ---- Count-based Circular ECFP ----

  /**
   * ECFP count fingerprint. Each element of the returned {@code int[]} array holds the
   * number of distinct substructure hashes that mapped to that bit position. This
   * preserves multiplicity information lost by binary (OR) folding.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #ecfpCounts(MolGraph, int, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the molecule
   * @param radius maximum expansion radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   */
  public static int[] ecfpCounts(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return ecfpCounts(toMolGraph(mol), radius, fpSize);
  }

  /**
   * ECFP count fingerprint from a pre-built {@link MolGraph}.
   *
   * @param g      the molecular graph
   * @param radius maximum expansion radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   */
  public static int[] ecfpCounts(MolGraph g, int radius, int fpSize) {
    checkCircular(g, radius, fpSize);
    int[] counts = new int[fpSize];
    if (g.n > 0) circular(g, ecfpAtomHashes(g), radius, fpSize, countStrategy(counts));
    return counts;
  }

  // ---- Count-based Circular FCFP ----

  /**
   * FCFP count fingerprint. Each element of the returned {@code int[]} array holds the
   * number of distinct substructure hashes that mapped to that bit position. This
   * preserves multiplicity information lost by binary (OR) folding.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #fcfpCounts(MolGraph, int, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the molecule
   * @param radius maximum expansion radius (2 = FCFP4, 3 = FCFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   */
  public static int[] fcfpCounts(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return fcfpCounts(toMolGraph(mol), radius, fpSize);
  }

  /**
   * FCFP count fingerprint from a pre-built {@link MolGraph}.
   *
   * @param g      the molecular graph
   * @param radius maximum expansion radius (2 = FCFP4, 3 = FCFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   */
  public static int[] fcfpCounts(MolGraph g, int radius, int fpSize) {
    checkCircular(g, radius, fpSize);
    int[] counts = new int[fpSize];
    if (g.n > 0) circular(g, fcfpAtomHashes(g), radius, fpSize, countStrategy(counts));
    return counts;
  }

  // ---- Count-based Mode dispatch ----

  /**
   * Count-based circular fingerprint with mode dispatch.
   *
   * @param mol    the molecule
   * @param radius maximum expansion radius
   * @param fpSize fingerprint length (number of counter bins)
   * @param mode   {@link Mode#ECFP} or {@link Mode#FCFP}
   * @return count array of length {@code fpSize}
   */
  public static int[] circularFingerprintCounts(IAtomContainer mol, int radius, int fpSize, Mode mode) {
    if (mol == null) throw new NullPointerException("mol");
    MolGraph g = toMolGraph(mol);
    return circularFingerprintCounts(g, radius, fpSize, mode);
  }

  /**
   * Count-based circular fingerprint with mode dispatch (MolGraph overload).
   *
   * @param g      the molecular graph
   * @param radius maximum expansion radius
   * @param fpSize fingerprint length (number of counter bins)
   * @param mode   {@link Mode#ECFP} or {@link Mode#FCFP}
   * @return count array of length {@code fpSize}
   */
  public static int[] circularFingerprintCounts(MolGraph g, int radius, int fpSize, Mode mode) {
    return (mode == Mode.FCFP) ? fcfpCounts(g, radius, fpSize) : ecfpCounts(g, radius, fpSize);
  }

  // ---- Topological Torsion Fingerprint (Nilakantan et al. 1987) ----

  /**
   * Compute the topological torsion atom type for atom {@code idx}:
   * encodes the atomic number, the number of heavy neighbours, the number of pi electrons and
   * ring membership as a single 32-bit integer.
   */
  private static int ttAtomType(MolGraph g, int idx) {
    int z = g.atomicNum[idx];
    // The degree counts the explicit neighbours, which are the heavy atoms of a prepared graph.
    int heavyDeg = g.degree[idx];
    // Pi electrons: count from bond orders (double bond contributes 1 pi, triple contributes 2)
    int nPi = 0;
    for (int nb : g.neighbors[idx]) {
      int bo = g.bondOrder(idx, nb);
      if (bo == 2) nPi += 1;
      else if (bo == 3) nPi += 2;
    }
    if (g.aromatic[idx]) nPi++;
    int inRing = g.ring[idx] ? 1 : 0;
    return (z << 16) | (heavyDeg << 8) | (nPi << 4) | inRing;
  }

  /**
   * Topological torsion fingerprint (binary). Enumerates all 4-atom linear paths
   * A-B-C-D, hashes the atom types using FNV-1a with canonical ordering
   * (min of forward and reverse hash), and folds into {@code fpSize} bits.
   *
   * <p>Reference: Nilakantan et al., J. Chem. Inf. Comput. Sci. 1987, 27, 82-85.
   *
   * @param g      the molecular graph
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   */
  public static long[] topologicalTorsion(MolGraph g, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    int n = g.n;
    int numWords = (fpSize + 63) / 64;
    long[] fp = new long[numWords];
    if (n < 4) return fp;

    // Precompute atom types
    int[] atomType = new int[n];
    for (int i = 0; i < n; i++) atomType[i] = ttAtomType(g, i);

    // Enumerate all 4-atom linear paths: A -> B -> C -> D
    for (int a = 0; a < n; a++) {
      for (int b : g.neighbors[a]) {
        for (int c : g.neighbors[b]) {
          if (c == a) continue;
          for (int d : g.neighbors[c]) {
            if (d == b || d == a) continue;
            // Hash forward: A-B-C-D
            long fwd = FNV1A_SEED;
            fwd ^= atomType[a]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[b]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[c]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[d]; fwd *= FNV1A_PRIME;
            // Hash reverse: D-C-B-A
            long rev = FNV1A_SEED;
            rev ^= atomType[d]; rev *= FNV1A_PRIME;
            rev ^= atomType[c]; rev *= FNV1A_PRIME;
            rev ^= atomType[b]; rev *= FNV1A_PRIME;
            rev ^= atomType[a]; rev *= FNV1A_PRIME;
            // Canonical ordering: take minimum
            long h = Long.compareUnsigned(fwd, rev) <= 0 ? fwd : rev;
            int bit = (int) (Long.remainderUnsigned(h, fpSize));
            fp[bit >> 6] |= 1L << (bit & 63);
          }
        }
      }
    }
    return fp;
  }

  /**
   * Topological torsion count fingerprint. Each element of the returned array holds
   * the number of distinct 4-atom torsion paths that hash to that bin position.
   *
   * @param g      the molecular graph
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   */
  public static int[] topologicalTorsionCounts(MolGraph g, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    int n = g.n;
    int[] counts = new int[fpSize];
    if (n < 4) return counts;

    int[] atomType = new int[n];
    for (int i = 0; i < n; i++) atomType[i] = ttAtomType(g, i);

    // Use a set to avoid double-counting the same unordered path {a,b,c,d}
    // since we enumerate both directions. Only count when a < d (canonical path direction).
    for (int a = 0; a < n; a++) {
      for (int b : g.neighbors[a]) {
        for (int c : g.neighbors[b]) {
          if (c == a) continue;
          for (int d : g.neighbors[c]) {
            if (d == b || d == a) continue;
            if (a > d) continue; // canonical: only count path once
            long fwd = FNV1A_SEED;
            fwd ^= atomType[a]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[b]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[c]; fwd *= FNV1A_PRIME;
            fwd ^= atomType[d]; fwd *= FNV1A_PRIME;
            long rev = FNV1A_SEED;
            rev ^= atomType[d]; rev *= FNV1A_PRIME;
            rev ^= atomType[c]; rev *= FNV1A_PRIME;
            rev ^= atomType[b]; rev *= FNV1A_PRIME;
            rev ^= atomType[a]; rev *= FNV1A_PRIME;
            long h = Long.compareUnsigned(fwd, rev) <= 0 ? fwd : rev;
            int idx = (int) (Long.remainderUnsigned(h, fpSize));
            counts[idx]++;
          }
        }
      }
    }
    return counts;
  }

  /**
   * Compute a topological torsion fingerprint from a CDK molecule.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #topologicalTorsion(MolGraph, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the CDK molecule
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   */
  public static long[] topologicalTorsion(IAtomContainer mol, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    return topologicalTorsion(toMolGraph(mol), fpSize);
  }

  // ---- Configurable fingerprints ----

  /**
   * The bit fingerprint of a graph under a {@link FingerprintOptions} specification. The classic
   * configuration of each family, {@link FingerprintOptions#ecfp()}, {@link FingerprintOptions#fcfp()},
   * {@link FingerprintOptions#path()} and {@link FingerprintOptions#torsion()}, gives the bits of
   * {@link #ecfp(MolGraph, int, int)}, {@link #fcfp(MolGraph, int, int)},
   * {@link #pathFingerprint(IAtomContainer, int, int)} and {@link #topologicalTorsion(MolGraph, int)};
   * {@link FingerprintOptions#maccs()} gives the 166 MACCS keys. The same options give the same bits
   * in C++ and Python.
   *
   * @param g       the molecular graph
   * @param options the fingerprint specification
   * @return binary fingerprint as long[] bitset, of {@code options.size} bits or one bit per key
   * @throws IllegalArgumentException if the options are inconsistent or a SMARTS key does not parse
   */
  public static long[] fingerprint(MolGraph g, FingerprintOptions options) {
    return OptionFingerprints.bits(g, options);
  }

  /**
   * The count fingerprint of a graph under a {@link FingerprintOptions} specification: one bin per
   * bit position of {@link #fingerprint(MolGraph, FingerprintOptions)}, or the number of distinct
   * matches of each key.
   */
  public static int[] fingerprintCounts(MolGraph g, FingerprintOptions options) {
    return OptionFingerprints.counts(g, options);
  }

  /**
   * The bit fingerprint of a CDK molecule under a {@link FingerprintOptions} specification; the
   * molecule is standardised first, as by {@link #toMolGraph(IAtomContainer)}.
   */
  public static long[] fingerprint(IAtomContainer mol, FingerprintOptions options) {
    if (mol == null) throw new NullPointerException("mol");
    return fingerprint(toMolGraph(mol), options);
  }

  /**
   * The count fingerprint of a CDK molecule under a {@link FingerprintOptions} specification; the
   * molecule is standardised first, as by {@link #toMolGraph(IAtomContainer)}.
   */
  public static int[] fingerprintCounts(IAtomContainer mol, FingerprintOptions options) {
    if (mol == null) throw new NullPointerException("mol");
    return fingerprintCounts(toMolGraph(mol), options);
  }

  // ---- Shared radius expansion ----

  /**
   * Strategy interface for recording a substructure hash during radius expansion.
   * Binary fingerprints OR a bit; count fingerprints increment an int counter.
   */
  @FunctionalInterface
  interface BitStrategy {
    /**
     * Record the hash {@code h} into the fingerprint state.
     * @param h the 64-bit substructure hash
     * @param fpSize the fingerprint size in bits (or elements)
     */
    void record(long h, int fpSize);
  }

  /** Largest finite shortest-path distance, including disconnected components. */
  private static int maxComponentExtent(MolGraph graph) {
    int extent = 0;
    int[] distances = new int[graph.n], queue = new int[graph.n];
    for (int start = 0; start < graph.n; start++) {
      Arrays.fill(distances, -1);
      distances[start] = 0;
      int next = 0, count = 1;
      queue[0] = start;
      while (next < count) {
        int atom = queue[next++];
        for (int neighbour : graph.neighbors[atom]) {
          if (distances[neighbour] >= 0) continue;
          distances[neighbour] = distances[atom] + 1;
          extent = Math.max(extent, distances[neighbour]);
          queue[count++] = neighbour;
        }
      }
    }
    return extent;
  }

  /**
   * Shared radius expansion core. Walks radii 1..maxRadius, computing neighbour-sorted
   * FNV-1a hashes, and delegates each hash to the supplied {@code strategy}.
   */
  private static void expandRadiiImpl(MolGraph g, long[] atomHash,
                                          int radius, int fpSize, int n,
                                          BitStrategy strategy) {
    int maxRadius = (radius < 0) ? maxComponentExtent(g) : radius;
    long[] prevHash = atomHash.clone();
    long[] nextHash = new long[n];
    int maxDegree = 0;
    for (int[] row : g.neighbors) maxDegree = Math.max(maxDegree, row.length);
    long[] nbHashes = new long[maxDegree];

    for (int r = 1; r <= maxRadius; r++) {
      for (int i = 0; i < n; i++) {
        int deg = g.neighbors[i].length;
        if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
        // The neighbour hashes sort as unsigned, as the native core sorts uint64_t: the sign bit
        // is flipped for the signed sort and flipped back when the hash is read.
        for (int k = 0; k < deg; k++) {
          int nb = g.neighbors[i][k];
          nbHashes[k] = (prevHash[nb] ^ (g.bondOrder(i, nb) * 31L)) ^ Long.MIN_VALUE;
        }
        Arrays.sort(nbHashes, 0, deg);

        long h = FNV1A_SEED;
        h ^= r;            h *= FNV1A_PRIME;
        h ^= prevHash[i];  h *= FNV1A_PRIME;
        for (int k = 0; k < deg; k++) { h ^= nbHashes[k] ^ Long.MIN_VALUE; h *= FNV1A_PRIME; }
        nextHash[i] = h;

        strategy.record(h, fpSize);
      }
      System.arraycopy(nextHash, 0, prevHash, 0, n);
    }
  }

  /** Binary bit-setting strategy: OR the bit into a long[] bitset. */
  private static BitStrategy binaryStrategy(long[] fp) {
    return (h, fpSize) -> {
      int bit = (int) (Long.remainderUnsigned(h, fpSize));
      long mask = 1L << (bit & 63);
      fp[bit >> 6] |= mask;
    };
  }

  /** Count-based strategy: increment an int[] counter at the hashed position. */
  private static BitStrategy countStrategy(int[] counts) {
    return (h, fpSize) -> {
      int idx = (int) (Long.remainderUnsigned(h, fpSize));
      counts[idx]++;
    };
  }

  // ---- Pharmacophore Classification (FCFP) ----

  /**
   * Classify an atom into pharmacophoric feature classes.
   * Returns bitmask: bit 0=donor, 1=acceptor, 2=positive, 3=negative, 4=aromatic, 5=hydrophobic.
   */
  static int classifyPharmacophore(MolGraph g, int idx) {
    int z = g.atomicNum[idx];
    int charge = g.formalCharge[idx];
    boolean arom = g.aromatic[idx];
    int deg = g.degree[idx];
    int features = 0;

    long attachedHydrogens = g.attachedHydrogenCount(idx);
    if ((z == 7 || z == 8 || z == 16) && attachedHydrogens > 0) features |= 1;
    boolean isPyrroleTypeN = z == 7 && arom && attachedHydrogens > 0;
    boolean isAcceptorN = (z == 7 && charge <= 0 && !isPyrroleTypeN);
    boolean isAcceptorS = (z == 16 && !arom);
    if (isAcceptorN || z == 8 || z == 9 || isAcceptorS) features |= 2;

    // Positive ionisable: basic amines only
    // Exclude: amide N (adjacent to C=O), aniline N (bonded to aromatic C)
    if (z == 7 && charge > 0) features |= 4;
    if (z == 7 && !arom && deg <= 3 && charge == 0) {
      boolean isAmide = false;
      boolean isAniline = false;
      boolean singleBonds = true;
      for (int nb : g.neighbors[idx]) {
        if (g.bondOrder(idx, nb) != 1) { singleBonds = false; break; }
        if (g.atomicNum[nb] == 6) {
          if (g.aromatic[nb]) { isAniline = true; break; } // aniline: N bonded to aromatic C
          for (int nb2 : g.neighbors[nb]) {
            if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2) {
              isAmide = true; break;
            }
          }
        }
        if (isAmide) break;
      }
      if (singleBonds && !isAmide && !isAniline) features |= 4;
    }

    // Negative ionisable: carboxylic, phosphoric, sulfonic acid
    if (z == 8 && charge < 0) features |= 8;
    if (z == 8 && charge == 0 && deg == 1) {
      for (int nb : g.neighbors[idx]) {
        int nbZ = g.atomicNum[nb];
        // Carboxylic acid: O=C with another O (require C=O, not ester C-O)
        if (nbZ == 6 && g.bondOrder(idx, nb) == 1) { // this O is single-bonded (the OH)
          for (int nb2 : g.neighbors[nb]) {
            if (nb2 != idx && g.atomicNum[nb2] == 8 && g.bondOrder(nb, nb2) == 2) features |= 8;
          }
        }
        // Phosphoric acid: O on P; Sulfonic acid: O on S
        if (nbZ == 15 || nbZ == 16) features |= 8;
      }
    }

    // Aromatic
    if (arom) features |= 16;

    // Hydrophobic: non-aromatic C with no heteroatom neighbours, or halogens (Cl, Br, I)
    if (z == 6 && !arom) {
      boolean hasHetero = false;
      for (int nb : g.neighbors[idx]) {
        int nbZ = g.atomicNum[nb];
        if (nbZ != 6 && nbZ != 1) { hasHetero = true; break; }
      }
      if (!hasHetero) features |= 32;
    }
    if (z == 17 || z == 35 || z == 53) features |= 32; // Cl, Br, I

    return features;
  }

  // ---- Format Conversions ----

  public static java.util.BitSet toBitSet(long[] fp) {
    return java.util.BitSet.valueOf(fp);
  }

  public static long[] fromBitSet(java.util.BitSet bs) {
    return bs.toLongArray();
  }

  /**
   * Convert a sparse count map {@code {bitPosition -> count}} into a dense integer
   * array of length {@code fpSize}. Useful for serialising count-based fingerprints
   * to fixed-length vectors suitable for database storage or NumPy interop.
   *
   * <p>Bit positions outside {@code [0, fpSize)} are silently ignored.
   *
   * <pre>{@code
   * Map<Integer,Integer> counts = new HashMap<>();
   * counts.put(42, 3);
   * counts.put(1000, 1);
   * int[] arr = FingerprintEngine.countsToArray(counts, 2048);
   * // arr[42] == 3, arr[1000] == 1, all others == 0
   * }</pre>
   *
   * @param counts sparse count map (bitPosition to count)
   * @param fpSize target array length
   * @return dense count array of length {@code fpSize}
   * @throws IllegalArgumentException if fpSize is not positive
   */
  public static int[] countsToArray(java.util.Map<Integer, Integer> counts, int fpSize) {
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    int[] arr = new int[fpSize];
    if (counts == null) return arr;
    for (java.util.Map.Entry<Integer, Integer> entry : counts.entrySet()) {
      int pos = entry.getKey();
      if (pos >= 0 && pos < fpSize) {
        arr[pos] = entry.getValue();
      }
    }
    return arr;
  }

  public static String toHex(long[] fp) {
    StringBuilder sb = new StringBuilder(fp.length * 16);
    for (long word : fp) {
      String hex = Long.toHexString(word);
      for (int pad = hex.length(); pad < 16; pad++) sb.append('0');
      sb.append(hex);
    }
    return sb.toString();
  }

  public static long[] fromHex(String hex) {
    int words = (hex.length() + 15) / 16;
    long[] fp = new long[words];
    for (int i = 0; i < words; i++) {
      int start = i * 16;
      int end = Math.min(start + 16, hex.length());
      fp[i] = Long.parseUnsignedLong(hex.substring(start, end), 16);
    }
    return fp;
  }

  public static String toBinaryString(long[] fp, int fpSize) {
    StringBuilder sb = new StringBuilder(fpSize);
    for (int i = 0; i < fpSize; i++) {
      sb.append((fp[i >> 6] & (1L << (i & 63))) != 0 ? '1' : '0');
    }
    return sb.toString();
  }

  // ---- Internal ----

  /**
   * Standardise the molecule (clone + atom-type + implicit H + aromaticity) for
   * fingerprint consistency. Falls back to the raw molecule on any failure.
   */
  private static IAtomContainer standardiseSafe(IAtomContainer mol) {
    try { return Standardiser.standardise(mol, Standardiser.TautomerMode.NONE); }
    catch (Exception e) { return mol; }
  }

  /** The hash of the first {@code length} atoms of a path. */
  @FunctionalInterface
  private interface PathHash {
    int of(int[] path, int length);
  }

  /**
   * Set the bit of every simple path of up to {@code maxLength} bonds over the atoms of
   * {@code adj} that {@code heavy} admits (every atom when it is null), starting from each atom.
   */
  private static void walkPaths(int[][] adj, boolean[] heavy, int maxLength, PathHash hash, long[] fp, int fpSize) {
    boolean[] visited = new boolean[adj.length];
    int[] path = new int[maxLength + 1];
    for (int start = 0; start < adj.length; start++) {
      if (heavy != null && !heavy[start]) continue;
      path[0] = start;
      visited[start] = true;
      fpSetBit(fp, hash.of(path, 1), fpSize);
      if (maxLength > 0) extendPath(adj, heavy, visited, path, 1, maxLength, hash, fp, fpSize);
      visited[start] = false;
    }
  }

  private static void extendPath(int[][] adj, boolean[] heavy, boolean[] visited, int[] path, int depth,
      int maxDepth, PathHash hash, long[] fp, int fpSize) {
    for (int nb : adj[path[depth - 1]]) {
      if (visited[nb] || (heavy != null && !heavy[nb])) continue;
      path[depth] = nb;
      visited[nb] = true;
      fpSetBit(fp, hash.of(path, depth + 1), fpSize);
      if (depth < maxDepth) extendPath(adj, heavy, visited, path, depth + 1, maxDepth, hash, fp, fpSize);
      visited[nb] = false;
    }
  }

  static int fpHashPath(int[] atomicNum, int[] path, int len) {
    int fwd = 17, rev = 17;
    for (int i = 0; i < len; i++) {
      fwd = fwd * 31 + atomicNum[path[i]];
      rev = rev * 31 + atomicNum[path[len - 1 - i]];
    }
    return Math.min(fwd & 0x7FFFFFFF, rev & 0x7FFFFFFF);
  }

  static void fpSetBit(long[] fp, int hash, int fpSize) {
    int bit = hash % fpSize;
    fp[bit >>> 6] |= 1L << (bit & 63);
  }

  public static long[] mcsFingerprint(MolGraph g, int pathLength, int fpSize) {
    return mcsFingerprint(g, false, pathLength, fpSize);
  }

  private static long[] mcsFingerprint(MolGraph g, boolean tautomerAware, int pathLength, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (pathLength < 0 || fpSize < 1) throw new IllegalArgumentException("Invalid fingerprint dimensions");
    int words = (fpSize - 1) / 64 + 1;
    long[] fp = new long[words];
    int n = g.n;
    if (n == 0) return fp;
    int[] tautClass = g.tautomerClass != null ? g.tautomerClass : new int[n];
    // Build heavy-atom mask (skip explicit H)
    boolean[] isHeavy = new boolean[n];
    int heavyCount = 0;
    for (int i = 0; i < n; i++) { isHeavy[i] = g.atomicNum[i] != 1; if (isHeavy[i]) heavyCount++; }
    if (heavyCount == 0) return fp;
    int maxLength = Math.min(pathLength, heavyCount - 1);
    walkPaths(g.neighbors, isHeavy, maxLength,
        (path, length) -> mcsFpHashPath(g, tautClass, tautomerAware, path, length), fp, fpSize);
    return fp;
  }

  static int mcsFpHashPath(MolGraph g, int[] tautClass, boolean tautAware,
      int[] path, int len) {
    int fwd = 17, rev = 17;
    for (int i = 0; i < len; i++) {
      int fi = i, ri = len - 1 - i;
      int ahF = atomHash(g, tautClass, tautAware, path[fi]);
      int ahR = atomHash(g, tautClass, tautAware, path[ri]);
      fwd = fwd * 31 + ahF;
      rev = rev * 31 + ahR;
      if (i < len - 1) {
        int bondF = bondHash(g, path[fi], path[fi + 1]);
        int bondR = bondHash(g, path[ri], path[ri - 1]);
        fwd = fwd * 31 + bondF;
        rev = rev * 31 + bondR;
      }
    }
    return Math.min(fwd & 0x7FFFFFFF, rev & 0x7FFFFFFF);
  }

  static int atomHash(MolGraph g, int[] tautClass, boolean tautAware, int atom) {
    int h = 17;
    // Tautomer-aware: atoms in tautomer classes (C/N/O) hash to same label
    int label = g.atomicNum[atom];
    if (tautAware && tautClass[atom] >= 0) label = 999; // tautomer-invariant label
    h = h * 37 + label;
    h = h * 37 + (g.ring[atom] ? 1 : 0);
    h = h * 37 + (g.aromatic[atom] ? 1 : 0);
    h = h * 37 + (tautAware ? 0 : tautClass[atom]); // collapse tautomer classes in taut mode
    h = h * 37 + g.degree[atom];
    return h;
  }

  static int bondHash(MolGraph g, int from, int to) {
    int h = 17;
    h = h * 37 + g.bondOrder(from, to);
    h = h * 37 + (g.bondInRing(from, to) ? 1 : 0);
    h = h * 37 + (g.bondAromatic(from, to) ? 1 : 0);
    return h;
  }

  /**
   * Fill rate, chi-squared uniformity over sixteen buckets and a suggested size for a 1 % false
   * positive rate.
   *
   * @deprecated Use {@link #analyzeFpQuality(long[], int)}, which the Python API shares, for the
   *     set bits, density, expected collision rate and saturation of a fingerprint.
   */
  @Deprecated
  public static Map<String, Double> analyzeFingerprintQuality(long[] fp, int fpSize) {
    int filled = nominalBitCount(fp, fpSize);
    double fillRate = (double) filled / fpSize;

    // Chi-squared test for uniform distribution (divide into 16 buckets)
    int buckets = 16;
    int[] counts = new int[buckets];
    for (int b = 0; b < fpSize; b++) {
      if ((fp[b >> 6] & (1L << (b & 63))) != 0) {
        counts[(int)((long)b * buckets / fpSize)]++;
      }
    }
    double expected = (double) filled / buckets;
    double chi2 = 0;
    for (int c : counts) chi2 += (c - expected) * (c - expected) / Math.max(1, expected);

    Map<String, Double> result = new LinkedHashMap<>();
    result.put("bitsFilled", (double) filled);
    result.put("totalBits", (double) fpSize);
    result.put("fillRate", fillRate);
    result.put("chiSquared", chi2);

    // The size at which a 1 % false positive rate is expected, counting each set bit as one
    // feature, and at least 256 bits.
    result.put("suggestedSize", Math.max(256.0, hashedFingerprintSize(filled, 0.01)));

    return result;
  }

  /**
   * Bits needed so that a feature absent from a molecule finds its bit set with probability of
   * about {@code falsePositiveRate}, when each of {@code features} features sets one uniformly
   * hashed bit: m = -n / ln(1 - p), rounded up to a multiple of 64.
   */
  private static double hashedFingerprintSize(double features, double falsePositiveRate) {
    if (!(falsePositiveRate > 0.0 && falsePositiveRate < 1.0))
      throw new IllegalArgumentException("targetFPR must lie strictly between 0 and 1");
    return Math.ceil(-features / Math.log1p(-falsePositiveRate) / 64) * 64;
  }

  /** Count only declared bits; unused padding and trailing words are excluded. */
  private static int nominalBitCount(long[] fingerprint, int size) {
    if (fingerprint == null) throw new NullPointerException("fp");
    if (size <= 0 || size > fingerprint.length * 64L)
      throw new IllegalArgumentException("fpSize is outside the fingerprint capacity");
    int count = 0, words = size / 64;
    for (int word = 0; word < words; word++) count += Long.bitCount(fingerprint[word]);
    if ((size & 63) != 0) count += Long.bitCount(fingerprint[words] & ((1L << (size & 63)) - 1));
    return count;
  }

  /**
   * A path fingerprint size for molecules of {@code avgAtomsPerMol} atoms, at which a path absent
   * from a molecule finds its bit set with probability of about {@code targetFPR}. The size allows
   * seven distinct paths per atom, more than most drug-like molecules have at the default path
   * length, each setting one hashed bit: m = -n / ln(1 - targetFPR) for n = 7 * avgAtomsPerMol,
   * rounded up to a multiple of 64. {@code librarySize} does not change the result.
   *
   * @throws IllegalArgumentException if {@code avgAtomsPerMol} is below 1, {@code targetFPR} is not
   *     strictly between 0 and 1, or the size needed exceeds the int range
   */
  public static int suggestFingerprintSize(int avgAtomsPerMol, int librarySize, double targetFPR) {
    if (avgAtomsPerMol < 1)
      throw new IllegalArgumentException("avgAtomsPerMol must be at least 1");
    double size = hashedFingerprintSize(7.0 * avgAtomsPerMol, targetFPR);
    if (size > Integer.MAX_VALUE)
      throw new IllegalArgumentException(avgAtomsPerMol + " atoms at targetFPR " + targetFPR
          + " need more bits than an int holds");
    return (int) size;
  }
}

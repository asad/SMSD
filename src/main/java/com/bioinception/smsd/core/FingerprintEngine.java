/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.Collections;
import java.util.Map;

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
 * @since 6.0
 */
public final class FingerprintEngine {

  /** FNV-1a 64-bit seed (offset basis). */
  static final long FNV1A_SEED = -3750763034362895579L;
  /** FNV-1a 64-bit prime. */
  static final long FNV1A_PRIME = 1099511628211L;

  private FingerprintEngine() {} // static utility

  // ---- Fingerprint caches (v6.5.3 perf fix) ----
  // Avoids recomputation when same MolGraph is fingerprinted 8-10x per reaction.
  // WeakHashMap: entries auto-evicted when MolGraph is GC'd.
  private static final Map<MolGraph, long[]> ecfpDefaultCache =
      Collections.synchronizedMap(new java.util.WeakHashMap<>());
  private static final Map<MolGraph, long[]> fcfpDefaultCache =
      Collections.synchronizedMap(new java.util.WeakHashMap<>());

  /** Clear fingerprint caches. */
  public static void clearFPCache() { ecfpDefaultCache.clear(); fcfpDefaultCache.clear(); }

  /** Fingerprint invariant mode for circular fingerprints. */
  public enum Mode {
    /** ECFP: atomic number, degree, implicit H count, ring, aromaticity, charge, tautomer class. */
    ECFP,
    /** FCFP: pharmacophoric feature classes (H-bond donor/acceptor, ionisable, aromatic, hydrophobic). */
    FCFP
  }

  // ---- MolGraph conversion ----

  /**
   * Pre-convert an {@link IAtomContainer} to a {@link MolGraph} for batch use.
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
   * @since 6.3.2
   */
  public static MolGraph toMolGraph(IAtomContainer mol) {
    if (mol == null) throw new NullPointerException("mol");
    return new MolGraph(mol);
  }

  // ---- Path Fingerprint ----

  /**
   * Compute a path-based fingerprint with automatic standardisation.
   *
   * <p>The molecule is standardised (aromaticity, typing, implicit H) before fingerprinting
   * so that SMILES with different explicit-H representations or kekulised forms produce
   * identical fingerprints. For pre-standardised molecules, call
   * {@link SearchEngine#pathFingerprint(IAtomContainer, int, int)} directly to skip the
   * standardisation step.
   *
   * @param mol        the CDK molecule to fingerprint
   * @param pathLength maximum path length (number of bonds) to enumerate
   * @param fpSize     fingerprint size in bits (should be a power of 2)
   * @return the fingerprint as a long array (bit-packed)
   */
  public static long[] pathFingerprint(IAtomContainer mol, int pathLength, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (pathLength <= 0) throw new IllegalArgumentException("pathLength must be > 0");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return SearchEngine.pathFingerprint(standardiseSafe(mol), pathLength, fpSize);
  }

  // ---- MCS Fingerprint ----

  public static long[] mcsFingerprint(IAtomContainer mol, ChemOptions chem, int pathLength, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (pathLength <= 0) throw new IllegalArgumentException("pathLength must be > 0");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return SearchEngine.mcsFingerprint(mol, chem, pathLength, fpSize);
  }

  // ---- Subset & Similarity ----

  public static boolean fingerprintSubset(long[] query, long[] target) {
    if (query == null) throw new NullPointerException("query");
    if (target == null) throw new NullPointerException("target");
    return SearchEngine.fingerprintSubset(query, target);
  }

  public static double tanimoto(long[] fp1, long[] fp2) {
    return SearchEngine.mcsFingerprintSimilarity(fp1, fp2);
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
   * @return Dice similarity in [0.0, 1.0]; returns 0.0 for null or empty inputs
   */
  public static double dice(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp2.length == 0) return 0.0;
    int minLen = Math.min(fp1.length, fp2.length);
    int maxLen = Math.max(fp1.length, fp2.length);
    int andBits = 0, aBits = 0, bBits = 0;
    for (int i = 0; i < minLen; i++) {
      andBits += Long.bitCount(fp1[i] & fp2[i]);
      aBits   += Long.bitCount(fp1[i]);
      bBits   += Long.bitCount(fp2[i]);
    }
    long[] longer = (fp1.length > fp2.length) ? fp1 : fp2;
    for (int i = minLen; i < maxLen; i++) {
      if (longer == fp1) aBits += Long.bitCount(longer[i]);
      else               bBits += Long.bitCount(longer[i]);
    }
    int sum = aBits + bBits;
    return (sum == 0) ? 0.0 : (2.0 * andBits) / sum;
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
   * @return cosine similarity in [0.0, 1.0]; returns 0.0 for null or empty inputs
   */
  public static double cosine(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp2.length == 0) return 0.0;
    int minLen = Math.min(fp1.length, fp2.length);
    int maxLen = Math.max(fp1.length, fp2.length);
    int andBits = 0, aBits = 0, bBits = 0;
    for (int i = 0; i < minLen; i++) {
      andBits += Long.bitCount(fp1[i] & fp2[i]);
      aBits   += Long.bitCount(fp1[i]);
      bBits   += Long.bitCount(fp2[i]);
    }
    long[] longer = (fp1.length > fp2.length) ? fp1 : fp2;
    for (int i = minLen; i < maxLen; i++) {
      if (longer == fp1) aBits += Long.bitCount(longer[i]);
      else               bBits += Long.bitCount(longer[i]);
    }
    double denom = Math.sqrt((double) aBits * bBits);
    return (denom == 0.0) ? 0.0 : andBits / denom;
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
   * @return Soergel distance in [0.0, 1.0]; returns 0.0 for null or empty inputs
   */
  public static double soergel(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null || fp1.length == 0 || fp2.length == 0) return 1.0; // max distance
    return 1.0 - tanimoto(fp1, fp2);
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

  // ---- Default valence helper (Rogers & Hahn invariant #3: implicit H count) ----

  /**
   * Returns the list of allowed valences for common organic-subset elements,
   * following the OpenSMILES / Daylight convention.  Multi-valent elements
   * (N, P, S, As, Se, Sb, Te, I) list every standard valence in ascending
   * order so that {@link #implicitH} can pick the smallest fitting one.
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
   * Computes implicit hydrogen count for an atom using the OpenSMILES
   * multi-valence model: pick the smallest allowed valence {@code v} such
   * that {@code v - |charge| >= bondOrderSum}, then
   * {@code implicitH = (v - |charge|) - bondOrderSum}.
   * Returns 0 for elements with no known default valence or when no
   * valence can accommodate the observed bond order sum.
   */
  static int implicitH(int atomicNum, int bondOrderSum, int formalCharge) {
    int[] valences = (atomicNum >= 0 && atomicNum < DEFAULT_VALENCES.length)
                     ? DEFAULT_VALENCES[atomicNum] : null;
    if (valences == null) return 0;
    int absCharge = Math.abs(formalCharge);
    for (int v : valences) {
      int capacity = v - absCharge;
      if (capacity >= bondOrderSum) return capacity - bondOrderSum;
    }
    return 0;
  }

  // ---- Circular ECFP ----

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
    return ecfp(new MolGraph(mol), radius, fpSize);
  }

  public static long[] ecfp(MolGraph g, int radius, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    if (radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
    // Cache lookup for default params (v6.5.3 perf fix: 8-10x ECFP/reaction)
    if (radius == 2 && fpSize == 2048) {
      long[] cached = ecfpDefaultCache.get(g);
      if (cached != null) return cached.clone();
    }
    int n = g.n;
    int numWords = (fpSize + 63) / 64;
    long[] fp = new long[numWords];
    if (n == 0) return fp;

    long[] atomHash = new long[n];
    for (int i = 0; i < n; i++) {
      long h = FNV1A_SEED;
      h ^= g.atomicNum[i];       h *= FNV1A_PRIME;  // R&H #1: atomic number
      h ^= g.degree[i];          h *= FNV1A_PRIME;  // R&H #2: heavy atom degree
      int bondOrdSum = 0;
      for (int nb : g.neighbors[i]) bondOrdSum += g.bondOrder(i, nb);
      h ^= bondOrdSum;            h *= FNV1A_PRIME;  // R&H #3: bond order sum (valence)
      h ^= g.massNumber[i];       h *= FNV1A_PRIME;  // R&H #4: atomic mass number
      h ^= g.formalCharge[i] + 4; h *= FNV1A_PRIME;  // R&H #5: formal charge
      int hCount = implicitH(g.atomicNum[i], bondOrdSum, g.formalCharge[i]);
      h ^= hCount;                h *= FNV1A_PRIME;  // R&H #6: attached H count
      h ^= g.ring[i] ? 1 : 0;   h *= FNV1A_PRIME;  // R&H #7: ring membership
      h ^= g.aromatic[i] ? 1 : 0; h *= FNV1A_PRIME; // Daylight extension: aromaticity
      if (g.tautomerClass != null && g.tautomerClass[i] >= 0) {
        h ^= 0xCAFE00L + g.tautomerClass[i]; h *= FNV1A_PRIME;
      }
      atomHash[i] = h;
      int bit = (int) (Long.remainderUnsigned(h, fpSize));
      fp[bit >> 6] |= 1L << (bit & 63);
    }

    long[] result = expandRadii(g, fp, atomHash, radius, fpSize, n);
    // Store in cache for default params
    if (radius == 2 && fpSize == 2048) ecfpDefaultCache.put(g, result.clone());
    return result;
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
    return fcfp(new MolGraph(mol), radius, fpSize);
  }

  public static long[] fcfp(MolGraph g, int radius, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    if (radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
    // Cache lookup for default params (v6.5.3 perf fix)
    if (radius == 2 && fpSize == 2048) {
      long[] cached = fcfpDefaultCache.get(g);
      if (cached != null) return cached.clone();
    }
    int n = g.n;
    int numWords = (fpSize + 63) / 64;
    long[] fp = new long[numWords];
    if (n == 0) return fp;

    // Use cached pharmacophore features (v6.5.3 perf fix)
    int[] pharmFeatures = g.getPharmacophoreFeatures();
    long[] atomHash = new long[n];
    for (int i = 0; i < n; i++) {
      long h = FNV1A_SEED;
      h ^= pharmFeatures[i];  h *= FNV1A_PRIME;
      h ^= g.degree[i];       h *= FNV1A_PRIME;
      h ^= g.ring[i] ? 1 : 0; h *= FNV1A_PRIME;
      atomHash[i] = h;
      int bit = (int) (Long.remainderUnsigned(h, fpSize));
      fp[bit >> 6] |= 1L << (bit & 63);
    }

    long[] result = expandRadii(g, fp, atomHash, radius, fpSize, n);
    if (radius == 2 && fpSize == 2048) fcfpDefaultCache.put(g, result.clone());
    return result;
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
    MolGraph g = new MolGraph(mol);
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
   * @since 6.0.1
   */
  public static int[] ecfpCounts(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return ecfpCounts(new MolGraph(mol), radius, fpSize);
  }

  /**
   * ECFP count fingerprint from a pre-built {@link MolGraph}.
   *
   * @param g      the molecular graph
   * @param radius maximum expansion radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   * @since 6.0.1
   */
  public static int[] ecfpCounts(MolGraph g, int radius, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    if (radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
    int n = g.n;
    int[] counts = new int[fpSize];
    if (n == 0) return counts;

    long[] atomHash = new long[n];
    for (int i = 0; i < n; i++) {
      long h = FNV1A_SEED;
      h ^= g.atomicNum[i];       h *= FNV1A_PRIME;  // R&H #1
      h ^= g.degree[i];          h *= FNV1A_PRIME;  // R&H #2
      int bondOrdSum = 0;
      for (int nb : g.neighbors[i]) bondOrdSum += g.bondOrder(i, nb);
      h ^= bondOrdSum;            h *= FNV1A_PRIME;  // R&H #3
      h ^= g.massNumber[i];       h *= FNV1A_PRIME;  // R&H #4
      h ^= g.formalCharge[i] + 4; h *= FNV1A_PRIME;  // R&H #5
      int hCount = implicitH(g.atomicNum[i], bondOrdSum, g.formalCharge[i]);
      h ^= hCount;                h *= FNV1A_PRIME;  // R&H #6
      h ^= g.ring[i] ? 1 : 0;   h *= FNV1A_PRIME;  // R&H #7
      h ^= g.aromatic[i] ? 1 : 0; h *= FNV1A_PRIME; // Daylight extension
      if (g.tautomerClass != null && g.tautomerClass[i] >= 0) {
        h ^= 0xCAFE00L + g.tautomerClass[i]; h *= FNV1A_PRIME;
      }
      atomHash[i] = h;
      int idx = (int) (Long.remainderUnsigned(h, fpSize));
      counts[idx]++;
    }

    expandRadiiImpl(g, atomHash, radius, fpSize, n, countStrategy(counts));
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
   * @since 6.0.1
   */
  public static int[] fcfpCounts(IAtomContainer mol, int radius, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    return fcfpCounts(new MolGraph(mol), radius, fpSize);
  }

  /**
   * FCFP count fingerprint from a pre-built {@link MolGraph}.
   *
   * @param g      the molecular graph
   * @param radius maximum expansion radius (2 = FCFP4, 3 = FCFP6; -1 = unlimited)
   * @param fpSize fingerprint length (number of counter bins)
   * @return count array of length {@code fpSize}
   * @since 6.0.1
   */
  public static int[] fcfpCounts(MolGraph g, int radius, int fpSize) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    if (radius < -1) throw new IllegalArgumentException("radius must be >= -1 (use -1 for unlimited)");
    int n = g.n;
    int[] counts = new int[fpSize];
    if (n == 0) return counts;

    // Use cached pharmacophore features (v6.5.3 perf fix)
    int[] pharmFeatures = g.getPharmacophoreFeatures();
    long[] atomHash = new long[n];
    for (int i = 0; i < n; i++) {
      long h = FNV1A_SEED;
      h ^= pharmFeatures[i];  h *= FNV1A_PRIME;
      h ^= g.degree[i];       h *= FNV1A_PRIME;
      h ^= g.ring[i] ? 1 : 0; h *= FNV1A_PRIME;
      atomHash[i] = h;
      int idx = (int) (Long.remainderUnsigned(h, fpSize));
      counts[idx]++;
    }

    expandRadiiImpl(g, atomHash, radius, fpSize, n, countStrategy(counts));
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
   * @since 6.0.1
   */
  public static int[] circularFingerprintCounts(IAtomContainer mol, int radius, int fpSize, Mode mode) {
    if (mol == null) throw new NullPointerException("mol");
    MolGraph g = new MolGraph(mol);
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
   * @since 6.0.1
   */
  public static int[] circularFingerprintCounts(MolGraph g, int radius, int fpSize, Mode mode) {
    return (mode == Mode.FCFP) ? fcfpCounts(g, radius, fpSize) : ecfpCounts(g, radius, fpSize);
  }

  // ---- Topological Torsion Fingerprint (Nilakantan et al. 1987) ----

  /**
   * Compute the topological torsion atom type for atom {@code idx}:
   * encodes (atomicNumber, numHeavyNeighbors, numPiElectrons, isRing)
   * as a single 32-bit integer.
   */
  private static int ttAtomType(MolGraph g, int idx) {
    int z = g.atomicNum[idx];
    // Heavy-atom degree (exclude implicit/explicit H — degree[] counts all explicit neighbors,
    // but MolGraph only stores heavy atoms so degree == numHeavyNeighbors)
    int heavyDeg = g.degree[idx];
    // Pi electrons: count from bond orders (double bond contributes 1 pi, triple contributes 2)
    int nPi = 0;
    for (int nb : g.neighbors[idx]) {
      int bo = g.bondOrder(idx, nb);
      if (bo == 2) nPi += 1;
      else if (bo == 3) nPi += 2;
    }
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
   * @since 6.0.1
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
   * @since 6.0.1
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
   * Convenience overload: compute topological torsion fingerprint from a CDK IAtomContainer.
   *
   * @param mol    the CDK molecule
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   * @since 6.0.1
   */
  /**
   * Compute a topological torsion fingerprint from a CDK molecule.
   *
   * <p>For batch workloads, pre-convert with {@link #toMolGraph(IAtomContainer)} once
   * and call {@link #topologicalTorsion(MolGraph, int)} to avoid repeated conversion overhead.
   *
   * @param mol    the CDK molecule
   * @param fpSize fingerprint size in bits
   * @return binary fingerprint as long[] bitset
   * @since 6.0.1
   */
  public static long[] topologicalTorsion(IAtomContainer mol, int fpSize) {
    if (mol == null) throw new NullPointerException("mol");
    return topologicalTorsion(new MolGraph(mol), fpSize);
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
     * @return {@code true} if this hash was previously unseen (for early termination)
     */
    boolean record(long h, int fpSize);
  }

  /**
   * Shared radius expansion core. Walks radii 1..maxRadius, computing neighbour-sorted
   * FNV-1a hashes, and delegates each hash to the supplied {@code strategy}.
   */
  private static boolean expandRadiiImpl(MolGraph g, long[] atomHash,
                                          int radius, int fpSize, int n,
                                          BitStrategy strategy) {
    int maxRadius = (radius < 0) ? n : radius;
    long[] prevHash = atomHash.clone();
    long[] nextHash = new long[n];

    for (int r = 1; r <= maxRadius; r++) {
      boolean anyNew = false;
      for (int i = 0; i < n; i++) {
        int deg = g.neighbors[i].length;
        if (deg == 0) { nextHash[i] = prevHash[i]; continue; }
        long[] nbHashes = new long[deg];
        for (int k = 0; k < deg; k++) {
          int nb = g.neighbors[i][k];
          nbHashes[k] = prevHash[nb] ^ (g.bondOrder(i, nb) * 31L);
        }
        java.util.Arrays.sort(nbHashes);

        long h = FNV1A_SEED;
        h ^= r;            h *= FNV1A_PRIME;
        h ^= prevHash[i];  h *= FNV1A_PRIME;
        for (long nh : nbHashes) { h ^= nh; h *= FNV1A_PRIME; }
        nextHash[i] = h;

        if (strategy.record(h, fpSize)) anyNew = true;
      }
      System.arraycopy(nextHash, 0, prevHash, 0, n);
      if (radius < 0 && !anyNew) break;
    }
    return true;
  }

  /** Binary bit-setting strategy: OR the bit into a long[] bitset. */
  private static BitStrategy binaryStrategy(long[] fp) {
    return (h, fpSize) -> {
      int bit = (int) (Long.remainderUnsigned(h, fpSize));
      long mask = 1L << (bit & 63);
      boolean isNew = (fp[bit >> 6] & mask) == 0;
      fp[bit >> 6] |= mask;
      return isNew;
    };
  }

  /** Count-based strategy: increment an int[] counter at the hashed position.
   *  Uses a separate HashSet to track genuinely new hashes (not just bin occupancy),
   *  preventing premature early-termination when radius=-1 (unlimited). */
  private static BitStrategy countStrategy(int[] counts) {
    java.util.Set<Long> seen = new java.util.HashSet<>();
    return (h, fpSize) -> {
      int idx = (int) (Long.remainderUnsigned(h, fpSize));
      counts[idx]++;
      return seen.add(h); // true if hash was not previously seen
    };
  }

  private static long[] expandRadii(MolGraph g, long[] fp, long[] atomHash,
                                     int radius, int fpSize, int n) {
    expandRadiiImpl(g, atomHash, radius, fpSize, n, binaryStrategy(fp));
    return fp;
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

    // H-bond donor: N-H, O-H, S-H (Rogers & Hahn 2010)
    if (z == 7 || z == 8 || z == 16) {
      int typicalValence = (z == 7) ? 3 : (z == 8) ? 2 : 2; // S typical valence 2
      if (deg < typicalValence && charge >= 0) features |= 1;
    }

    // H-bond acceptor: N (not pyrrole-type), O, F, S (not thiophene-type)
    // Pyrrole N: aromatic, has implicit H (lone pair in pi system) — NOT acceptor
    // Pyridine N: aromatic, no implicit H (lone pair in ring plane) — IS acceptor
    boolean isPyrroleTypeN = false;
    if (z == 7 && arom) {
      int bondSum = 0;
      for (int nb : g.neighbors[idx]) bondSum += g.bondOrder(idx, nb);
      int hCount = Math.max(0, 3 - bondSum - Math.abs(charge));
      isPyrroleTypeN = (hCount > 0); // has H = lone pair in pi system
    }
    boolean isAcceptorN = (z == 7 && charge <= 0 && !isPyrroleTypeN);
    boolean isAcceptorS = (z == 16 && !arom);
    if (isAcceptorN || z == 8 || z == 9 || isAcceptorS) features |= 2;

    // Positive ionisable: basic amines only
    // Exclude: amide N (adjacent to C=O), aniline N (bonded to aromatic C)
    if (z == 7 && charge > 0) features |= 4;
    if (z == 7 && !arom && deg <= 3 && charge == 0) {
      boolean isAmide = false;
      boolean isAniline = false;
      for (int nb : g.neighbors[idx]) {
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
      if (!isAmide && !isAniline) features |= 4;
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

    // Hydrophobic: non-aromatic C with no heteroatom neighbors, or halogens (Cl, Br, I)
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
   * @since 6.3.0
   */
  public static int[] countsToArray(java.util.Map<Integer, Integer> counts, int fpSize) {
    if (fpSize <= 0) throw new IllegalArgumentException("fpSize must be > 0");
    int[] arr = new int[fpSize];
    if (counts == null) return arr;
    for (var entry : counts.entrySet()) {
      int pos = entry.getKey();
      if (pos >= 0 && pos < fpSize) {
        arr[pos] = entry.getValue();
      }
    }
    return arr;
  }

  public static String toHex(long[] fp) {
    StringBuilder sb = new StringBuilder(fp.length * 16);
    for (long w : fp) sb.append(String.format("%016x", w));
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
}

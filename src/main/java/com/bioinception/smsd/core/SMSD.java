/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Main entry point for the SMSD library: substructure search, enumeration, and MCS computation.
 *
 * <p>Molecules are standardised on construction by default. For pre-standardised molecules, pass
 * {@code standardise=false} to skip preprocessing. For zero-overhead, call {@link SearchEngine}
 * directly.
 *
 * <pre>{@code
 * // Quick MCS example
 * SMSD smsd = new SMSD(mol1, mol2, new ChemOptions());
 * Map<Integer, Integer> mcs = smsd.findMCS();
 * System.out.println("MCS size: " + mcs.size());
 * }</pre>
 *
 * @author Syed Asad Rahman
 * @see ChemOptions
 * @see SearchEngine
 */
public final class SMSD {
  private static final SmilesGenerator CANONICAL_STEREO_SMILES =
      new SmilesGenerator(SmiFlavor.Canonical | SmiFlavor.Stereo);

  private final IAtomContainer query;
  private final IAtomContainer target;
  private final ChemOptions chem;
  private final String querySmarts;
  private final boolean selfMatch;
  private long substructureTimeoutMs = 8_000L;
  private long mcsTimeoutMs = 10_000L;
  /**
   * pKa-informed confidence of the most recent tautomer-aware MCS result, ∈ [0,1].
   * Updated automatically by every {@code findMCS()} call when {@code tautomerAware=true}.
   * Always {@code 1.0} for non-tautomeric searches.
   * @see #getTautomerConfidence()
   */
  private double tautomerConfidence = 1.0;

  /**
   * Create an SMSD instance from two molecule containers with automatic standardisation.
   *
   * <pre>{@code
   * IAtomContainer q = sp.parseSmiles("c1ccccc1");
   * IAtomContainer t = sp.parseSmiles("c1ccc(O)cc1");
   * SMSD smsd = new SMSD(q, t, new ChemOptions());
   * }</pre>
   *
   * @param query  the query molecule (smaller or equal in size to target for best performance)
   * @param target the target molecule to search or compare against
   * @param chem   chemical matching options; if {@code null}, defaults are used
   */
  public SMSD(IAtomContainer query, IAtomContainer target, ChemOptions chem) {
    this(query, target, chem, true);
  }

  /**
   * Create an SMSD instance with explicit control over standardisation.
   *
   * <p>Pass {@code standardise=false} for pre-processed molecules (~2-5x faster construction).
   *
   * <pre>{@code
   * // Skip standardisation for pre-processed molecules
   * SMSD smsd = new SMSD(preProcessedQuery, preProcessedTarget, opts, false);
   * }</pre>
   *
   * @param query       the query molecule
   * @param target      the target molecule
   * @param chem        chemical matching options; if {@code null}, defaults are used
   * @param standardise {@code true} to apply CDK standardisation (aromaticity, typing);
   *                    {@code false} to skip it for pre-processed molecules
   */
  public SMSD(IAtomContainer query, IAtomContainer target, ChemOptions chem, boolean standardise) {
    if (query == null) throw new NullPointerException("query");
    if (target == null) throw new NullPointerException("target");
    this.chem = chem != null ? chem : new ChemOptions();
    this.querySmarts = null;
    this.selfMatch = (query == target);
    this.query = standardise ? standardiseSafe(query) : query;
    this.target = this.selfMatch ? this.query : (standardise ? standardiseSafe(target) : target);
  }

  /**
   * Create an SMSD instance with a SMARTS query string. Target is standardised by default.
   *
   * <pre>{@code
   * SMSD smsd = new SMSD("[#6]~[#7]", targetMol, new ChemOptions());
   * boolean match = smsd.isSubstructure();
   * }</pre>
   *
   * @param querySmarts the SMARTS pattern to match against the target
   * @param target      the target molecule
   * @param chem        chemical matching options; if {@code null}, defaults are used
   */
  public SMSD(String querySmarts, IAtomContainer target, ChemOptions chem) {
    this(querySmarts, target, chem, true);
  }

  /**
   * Create an SMSD instance with a SMARTS query and explicit standardisation control.
   *
   * @param querySmarts the SMARTS pattern to match against the target
   * @param target      the target molecule
   * @param chem        chemical matching options; if {@code null}, defaults are used
   * @param standardise {@code true} to apply standardisation to the target molecule
   */
  public SMSD(String querySmarts, IAtomContainer target, ChemOptions chem, boolean standardise) {
    if (querySmarts == null || querySmarts.isEmpty())
      throw new IllegalArgumentException("querySmarts must not be null or empty");
    if (target == null) throw new NullPointerException("target");
    this.chem = chem != null ? chem : new ChemOptions();
    this.querySmarts = querySmarts;
    this.selfMatch = false;
    this.query = null;
    this.target = standardise ? standardiseSafe(target) : target;
  }

  /**
   * Convenience constructor: parse both SMILES strings and standardise.
   *
   * <pre>{@code
   * SMSD smsd = new SMSD("c1ccccc1", "c1ccc(O)cc1", new ChemOptions());
   * Map<Integer, Integer> mcs = smsd.findMCS();
   * }</pre>
   *
   * @param querySmiles  SMILES string for the query molecule
   * @param targetSmiles SMILES string for the target molecule
   * @param chem         chemical matching options; if {@code null}, defaults are used
   * @throws Exception if either SMILES string cannot be parsed
   */
  public SMSD(String querySmiles, String targetSmiles, ChemOptions chem) throws Exception {
    if (querySmiles == null || querySmiles.isEmpty())
      throw new IllegalArgumentException("querySmiles must not be null or empty");
    if (targetSmiles == null || targetSmiles.isEmpty())
      throw new IllegalArgumentException("targetSmiles must not be null or empty");
    SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
    this.chem = chem != null ? chem : new ChemOptions();
    this.querySmarts = null;
    this.selfMatch = false;
    String qSmi = this.chem.lenientSmiles ? MolGraph.sanitizeSmiles(querySmiles) : querySmiles;
    String tSmi = this.chem.lenientSmiles ? MolGraph.sanitizeSmiles(targetSmiles) : targetSmiles;
    this.query = standardiseSafe(sp.parseSmiles(qSmi));
    this.target = standardiseSafe(sp.parseSmiles(tSmi));
  }

  /**
   * Set the timeout for substructure search operations.
   *
   * <pre>{@code
   * smsd.setSubstructureTimeoutMs(5000); // 5 seconds
   * }</pre>
   *
   * @param ms timeout in milliseconds (default: 8000)
   */
  public void setSubstructureTimeoutMs(long ms) {
    if (ms <= 0) throw new IllegalArgumentException("timeout must be > 0, got " + ms);
    this.substructureTimeoutMs = ms;
  }

  /**
   * Set the timeout for MCS computation operations.
   *
   * <pre>{@code
   * smsd.setMcsTimeoutMs(30000); // 30 seconds for complex molecules
   * }</pre>
   *
   * @param ms timeout in milliseconds (default: 10000)
   */
  public void setMcsTimeoutMs(long ms) {
    if (ms <= 0) throw new IllegalArgumentException("timeout must be > 0, got " + ms);
    this.mcsTimeoutMs = ms;
  }

  // ---- Substructure ----

  /**
   * Check if the query is a substructure of the target using the default timeout.
   *
   * <pre>{@code
   * SMSD smsd = new SMSD(benzene, phenol, new ChemOptions());
   * boolean isSub = smsd.isSubstructure(); // true: benzene is in phenol
   * }</pre>
   *
   * @return {@code true} if the query is a substructure of the target
   */
  public boolean isSubstructure() {
    return isSubstructure(substructureTimeoutMs);
  }

  /**
   * Check if the query is a substructure of the target with a specified timeout.
   *
   * <pre>{@code
   * boolean found = smsd.isSubstructure(5000); // 5-second timeout
   * }</pre>
   *
   * @param timeoutMs maximum time in milliseconds before the search is abandoned
   * @return {@code true} if the query is a substructure of the target
   * @throws RuntimeException if SMARTS matching fails (for SMARTS-based queries)
   */
  public boolean isSubstructure(long timeoutMs) {
    if (querySmarts != null) {
      try { return !Standardiser.matchAll(querySmarts, target).isEmpty(); }
      catch (Exception e) { throw new RuntimeException("SMARTS matching failed", e); }
    }
    if (selfMatch) return true;
    return SearchEngine.isSubstructure(query, target, chem, timeoutMs);
  }

  /**
   * Enumerate all unique substructure mappings (query atom index to target atom index).
   *
   * <pre>{@code
   * List<Map<Integer, Integer>> maps = smsd.findAllSubstructures(10, 5000);
   * for (Map<Integer, Integer> m : maps) {
   *     System.out.println("Mapping: " + m);
   * }
   * }</pre>
   *
   * @param maxSolutions maximum number of mappings to return
   * @param timeoutMs    maximum time in milliseconds
   * @return list of atom-index mappings (query index to target index)
   * @throws UnsupportedOperationException if this instance was created with a SMARTS query
   */
  public List<Map<Integer, Integer>> findAllSubstructures(int maxSolutions, long timeoutMs) {
    requireMoleculeQuery("findAllSubstructures");
    return SearchEngine.findAllSubstructures(query, target, chem, maxSolutions, timeoutMs);
  }

  /**
   * Check substructure with telemetry (nodes visited, backtracks, pruning stats).
   *
   * <pre>{@code
   * SearchEngine.SubstructureResult result = smsd.isSubstructureWithStats(5000);
   * System.out.println("Match: " + result.exists());
   * System.out.println("Nodes visited: " + result.stats().nodesVisited());
   * }</pre>
   *
   * @param timeoutMs maximum time in milliseconds
   * @return result containing the match outcome and detailed search statistics
   * @throws UnsupportedOperationException if this instance was created with a SMARTS query
   */
  public SearchEngine.SubstructureResult isSubstructureWithStats(long timeoutMs) {
    requireMoleculeQuery("isSubstructureWithStats");
    return SearchEngine.isSubstructureWithStats(query, target, chem, timeoutMs);
  }

  /**
   * Enumerate substructure mappings with telemetry.
   *
   * <pre>{@code
   * SearchEngine.SubstructureResult result = smsd.findAllSubstructuresWithStats(100, 10000);
   * System.out.println("Found " + result.mappings().size() + " mappings");
   * System.out.println("Nodes: " + result.stats().nodesVisited());
   * }</pre>
   *
   * @param maxSolutions maximum number of mappings to enumerate
   * @param timeoutMs    maximum time in milliseconds
   * @return result containing mappings and detailed search statistics
   * @throws UnsupportedOperationException if this instance was created with a SMARTS query
   */
  public SearchEngine.SubstructureResult findAllSubstructuresWithStats(
      int maxSolutions, long timeoutMs) {
    requireMoleculeQuery("findAllSubstructuresWithStats");
    return SearchEngine.findAllSubstructuresWithStats(query, target, chem, maxSolutions, timeoutMs);
  }

  // ---- MCS ----

  /**
   * Find the Maximum Common Substructure using defaults (non-induced, connected, default timeout).
   *
   * <pre>{@code
   * SMSD smsd = new SMSD(mol1, mol2, new ChemOptions());
   * Map<Integer, Integer> mcs = smsd.findMCS();
   * System.out.println("MCS has " + mcs.size() + " atoms");
   * }</pre>
   *
   * @return mapping from query atom indices to target atom indices; empty if no common substructure
   */
  public Map<Integer, Integer> findMCS() {
    return findMCS(false, true, mcsTimeoutMs);
  }

  /**
   * Find the Maximum Common Substructure with full control over parameters.
   *
   * <pre>{@code
   * // Induced, connected MCS with 30-second timeout
   * Map<Integer, Integer> mcs = smsd.findMCS(true, true, 30000);
   * }</pre>
   *
   * @param induced       {@code true} for induced MCS (edge-preserving); {@code false} for
   *                      non-induced (allows edge mismatches)
   * @param connectedOnly {@code true} to require the MCS to be a connected subgraph
   * @param timeoutMs     maximum time in milliseconds for the computation
   * @return mapping from query atom indices to target atom indices; empty if no common substructure
   */
  public Map<Integer, Integer> findMCS(boolean induced, boolean connectedOnly, long timeoutMs) {
    requireMoleculeQuery("findMCS");
    if (timeoutMs <= 0) timeoutMs = mcsTimeoutMs;
    SearchEngine.MCSOptions opt = new SearchEngine.MCSOptions();
    opt.induced = induced;
    opt.connectedOnly = connectedOnly;
    opt.timeoutMs = timeoutMs;
    Map<Integer, Integer> result = exactIdentityMcsMapping();
    if (result.isEmpty()) result = SearchEngine.findMCS(query, target, chem, opt);
    // Compute pKa-informed tautomer confidence for tautomer-aware searches.
    tautomerConfidence = chem.tautomerAware
        ? SearchEngine.computeTautomerConfidence(query, target, result)
        : 1.0;
    return result;
  }

  private Map<Integer, Integer> exactIdentityMcsMapping() {
    if (query == null || target == null) return java.util.Collections.emptyMap();
    if (query.getAtomCount() != target.getAtomCount()) return java.util.Collections.emptyMap();

    String queryKey = canonicalStereoSmiles(query);
    String targetKey = canonicalStereoSmiles(target);
    if (queryKey == null || targetKey == null || !queryKey.equals(targetKey)) {
      return java.util.Collections.emptyMap();
    }

    MolGraph qg = SearchEngine.toMolGraph(query);
    MolGraph tg = SearchEngine.toMolGraph(target);

    Map<Integer, Integer> identity = new LinkedHashMap<>();
    for (int i = 0; i < query.getAtomCount(); i++) identity.put(i, i);
    if (!SearchEngine.validateMapping(qg, tg, identity, chem).isEmpty()) {
      return java.util.Collections.emptyMap();
    }
    return identity;
  }

  private static String canonicalStereoSmiles(IAtomContainer mol) {
    try {
      return CANONICAL_STEREO_SMILES.create(mol);
    } catch (Exception e) {
      return null;
    }
  }

  /**
   * Returns the pKa-informed tautomer confidence of the most recent {@link #findMCS()} call.
   *
   * <p>Confidence ∈ [0,1]: product of per-pair relevance weights for every cross-element
   * tautomeric atom match in the MCS.  A value of {@code 1.0} means either no tautomeric
   * flexibility was exercised, or all matched tautomeric atoms have the same element type.
   * Values below 1.0 indicate the degree to which the MCS depends on low-probability
   * tautomeric forms at the configured {@link ChemOptions#pH} (default 7.4).
   *
   * <pre>{@code
   * ChemOptions c = ChemOptions.tautomerProfile();
   * SMSD smsd = new SMSD(ketoMol, enolMol, c);
   * Map<Integer,Integer> mcs = smsd.findMCS();
   * double conf = smsd.getTautomerConfidence(); // e.g. 0.97
   * }</pre>
   *
   * @return confidence score; always {@code 1.0} when {@code tautomerAware=false}
   */
  public double getTautomerConfidence() { return tautomerConfidence; }

  /**
   * Enumerate multiple distinct MCS mappings of the maximum size.
   *
   * <p>First finds the optimal MCS size K, then collects up to {@code maxResults}
   * distinct atom-atom mappings of that same size. Useful for scaffold analysis
   * and SAR studies where alternative mappings reveal different chemical perspectives.
   *
   * <pre>{@code
   * SMSD smsd = new SMSD(mol1, mol2, new ChemOptions());
   * List<Map<Integer, Integer>> allMCS = smsd.findAllMCS(10);
   * System.out.println("Found " + allMCS.size() + " distinct MCS mappings");
   * }</pre>
   *
   * @param maxResults maximum number of distinct mappings to return (default 10)
   * @return list of distinct atom-index mappings, each of the maximum MCS size
   */
  public List<Map<Integer, Integer>> findAllMCS(int maxResults) {
    if (maxResults == 0) return java.util.Collections.emptyList();
    if (maxResults < 0) maxResults = 10;
    return findAllMCS(false, true, mcsTimeoutMs, maxResults);
  }

  /**
   * Enumerate multiple distinct MCS mappings with full control over parameters.
   *
   * <pre>{@code
   * List<Map<Integer, Integer>> allMCS = smsd.findAllMCS(true, true, 30000, 5);
   * }</pre>
   *
   * @param induced       {@code true} for induced MCS (edge-preserving)
   * @param connectedOnly {@code true} to require connected subgraph
   * @param timeoutMs     maximum time in milliseconds
   * @param maxResults    maximum number of distinct mappings to return
   * @return list of distinct atom-index mappings, each of the maximum MCS size
   */
  public List<Map<Integer, Integer>> findAllMCS(
      boolean induced, boolean connectedOnly, long timeoutMs, int maxResults) {
    requireMoleculeQuery("findAllMCS");
    SearchEngine.MCSOptions opt = new SearchEngine.MCSOptions();
    opt.induced = induced;
    opt.connectedOnly = connectedOnly;
    opt.timeoutMs = timeoutMs;
    return SearchEngine.findAllMCS(query, target, chem, opt, maxResults);
  }

  // ---- Similarity ----

  /**
   * Fast upper bound on Tanimoto atom similarity using atom-label frequency overlap.
   *
   * <p>Runs in O(V + E) time. Useful for pre-filtering molecule pairs in batch scenarios
   * before running expensive MCS computation.
   *
   * <pre>{@code
   * double ub = smsd.similarityUpperBound();
   * if (ub >= 0.5) {
   *     Map<Integer, Integer> mcs = smsd.findMCS(); // only compute MCS if promising
   * }
   * }</pre>
   *
   * @return upper bound on Tanimoto similarity in the range [0.0, 1.0]
   * @throws UnsupportedOperationException if this instance was created with a SMARTS query
   */
  public double similarityUpperBound() {
    requireMoleculeQuery("similarityUpperBound");
    return SearchEngine.similarityUpperBound(query, target, chem);
  }

  // ---- N-MCS (Multi-molecule) ----

  /**
   * Find the Maximum Common Substructure across multiple molecules with a coverage threshold.
   *
   * <pre>{@code
   * List<IAtomContainer> mols = Arrays.asList(mol1, mol2, mol3);
   * Map<Integer, Integer> nmcs = SMSD.findNMCS(mols, new ChemOptions(), 0.8, 10000);
   * System.out.println("Common core: " + nmcs.size() + " atoms");
   * }</pre>
   *
   * @param molecules list of molecules (at least 2)
   * @param chem      chemical matching options
   * @param threshold fraction of molecules that must contain the MCS (0.0 to 1.0)
   * @param timeoutMs timeout per pairwise MCS computation in milliseconds
   * @return identity mapping for the common core atom indices, or empty if no common core
   */
  public static Map<Integer, Integer> findNMCS(
      List<IAtomContainer> molecules, ChemOptions chem, double threshold, long timeoutMs) {
    if (molecules == null || molecules.isEmpty())
      throw new IllegalArgumentException("molecules must not be null or empty");
    if (threshold < 0.0 || threshold > 1.0)
      throw new IllegalArgumentException("threshold must be in [0, 1]");
    return SearchEngine.findNMCS(molecules, chem, threshold, timeoutMs);
  }

  /**
   * Find the N-molecule MCS and return the actual molecule representing the common core.
   *
   * <pre>{@code
   * IAtomContainer core = SMSD.findNMCSMolecule(mols, new ChemOptions(), 1.0, 10000);
   * if (core != null) {
   *     System.out.println("Core has " + core.getAtomCount() + " atoms");
   * }
   * }</pre>
   *
   * @param molecules list of molecules (at least 2)
   * @param chem      chemical matching options
   * @param threshold fraction of molecules that must contain the MCS (0.0 to 1.0)
   * @param timeoutMs timeout per pairwise MCS computation in milliseconds
   * @return the common core as an {@link IAtomContainer}; empty (zero atoms) if no common core
   */
  public static IAtomContainer findNMCSMolecule(
      List<IAtomContainer> molecules, ChemOptions chem, double threshold, long timeoutMs) {
    if (molecules == null || molecules.isEmpty())
      throw new IllegalArgumentException("molecules must not be null or empty");
    if (threshold < 0.0 || threshold > 1.0)
      throw new IllegalArgumentException("threshold must be in [0, 1]");
    return SearchEngine.findNMCSMolecule(molecules, chem, threshold, timeoutMs);
  }

  // ---- R-Group Decomposition ----

  /**
   * Decompose molecules into a core scaffold and R-groups.
   *
   * <pre>{@code
   * List<Map<String, IAtomContainer>> decomp =
   *     SMSD.decomposeRGroups(core, molecules, new ChemOptions(), 5000);
   * for (Map<String, IAtomContainer> d : decomp) {
   *     System.out.println("Core: " + d.get("core").getAtomCount());
   *     System.out.println("R1: " + d.getOrDefault("R1", null));
   * }
   * }</pre>
   *
   * @param core      the core scaffold molecule
   * @param molecules list of molecules to decompose
   * @param chem      chemical matching options
   * @param timeoutMs timeout for substructure matching in milliseconds
   * @return list of decomposition maps; each contains "core" and "R1", "R2", etc.
   */
  public static List<Map<String, IAtomContainer>> decomposeRGroups(
      IAtomContainer core, List<IAtomContainer> molecules, ChemOptions chem, long timeoutMs) {
    return SearchEngine.decomposeRGroups(core, molecules, chem, timeoutMs);
  }

  // ---- Fingerprints (delegated to FingerprintEngine) ----

  /**
   * Pre-convert an {@link IAtomContainer} to a {@link MolGraph} for batch fingerprinting.
   *
   * <p>When computing multiple fingerprints for the same molecule, convert once and pass
   * the {@link MolGraph} to each fingerprint method to avoid redundant conversion overhead.
   *
   * @param mol the CDK molecule to convert
   * @return a MolGraph suitable for all fingerprint MolGraph overloads
   * @since 6.3.2
   * @see FingerprintEngine#toMolGraph(IAtomContainer)
   */
  public static MolGraph toMolGraph(IAtomContainer mol) {
    return FingerprintEngine.toMolGraph(mol);
  }

  /**
   * Compute a path-based fingerprint with automatic standardisation.
   *
   * <p>The molecule is standardised before fingerprinting for consistent results
   * regardless of input representation. For pre-standardised molecules, call
   * {@link SearchEngine#pathFingerprint(IAtomContainer, int, int)} directly.
   *
   * @param mol        the CDK molecule to fingerprint
   * @param pathLength maximum path length (number of bonds) to enumerate
   * @param fpSize     fingerprint size in bits (should be a power of 2)
   * @return the fingerprint as a long array (bit-packed)
   * @see FingerprintEngine#pathFingerprint
   */
  public static long[] pathFingerprint(IAtomContainer mol, int pathLength, int fpSize) {
    return FingerprintEngine.pathFingerprint(mol, pathLength, fpSize);
  }

  /** @see FingerprintEngine#fingerprintSubset */
  public static boolean fingerprintSubset(long[] query, long[] target) {
    return FingerprintEngine.fingerprintSubset(query, target);
  }

  /** @see FingerprintEngine#mcsFingerprint */
  public static long[] mcsFingerprint(IAtomContainer mol, ChemOptions chem, int pathLength, int fpSize) {
    return FingerprintEngine.mcsFingerprint(mol, chem, pathLength, fpSize);
  }

  /** @see FingerprintEngine#tanimoto */
  public static double mcsFingerprintSimilarity(long[] fp1, long[] fp2) {
    return FingerprintEngine.tanimoto(fp1, fp2);
  }

  /** Fingerprint invariant mode. @deprecated Use {@link FingerprintEngine.Mode} instead. */
  @Deprecated
  public enum FingerprintMode { ECFP, FCFP }

  /**
   * Compute an ECFP circular fingerprint. For batch workloads, pre-convert with
   * {@link #toMolGraph(IAtomContainer)} once and use the MolGraph overload.
   * @see FingerprintEngine#circularFingerprint(IAtomContainer, int, int) */
  public static long[] circularFingerprint(IAtomContainer mol, int radius, int fpSize) {
    return FingerprintEngine.ecfp(mol, radius, fpSize);
  }

  /**
   * Compute an ECFP circular fingerprint. For batch workloads, pre-convert with
   * {@link #toMolGraph(IAtomContainer)} once and use the MolGraph overload.
   * @see FingerprintEngine#ecfp(IAtomContainer, int, int) */
  public static long[] circularFingerprintECFP(IAtomContainer mol, int radius, int fpSize) {
    return FingerprintEngine.ecfp(mol, radius, fpSize);
  }

  /**
   * Compute an FCFP circular fingerprint. For batch workloads, pre-convert with
   * {@link #toMolGraph(IAtomContainer)} once and use the MolGraph overload.
   * @see FingerprintEngine#fcfp(IAtomContainer, int, int) */
  public static long[] circularFingerprintFCFP(IAtomContainer mol, int radius, int fpSize) {
    return FingerprintEngine.fcfp(mol, radius, fpSize);
  }

  /** @see FingerprintEngine#ecfp(MolGraph, int, int) */
  public static long[] circularFingerprint(MolGraph g, int radius, int fpSize) {
    return FingerprintEngine.ecfp(g, radius, fpSize);
  }

  /** @see FingerprintEngine#ecfp(MolGraph, int, int) */
  public static long[] circularFingerprintECFP(MolGraph g, int radius, int fpSize) {
    return FingerprintEngine.ecfp(g, radius, fpSize);
  }

  /** @see FingerprintEngine#fcfp(MolGraph, int, int) */
  public static long[] circularFingerprintFCFP(MolGraph g, int radius, int fpSize) {
    return FingerprintEngine.fcfp(g, radius, fpSize);
  }

  public static long[] circularFingerprint(IAtomContainer mol, int radius, int fpSize, FingerprintMode mode) {
    return FingerprintEngine.circularFingerprint(mol, radius, fpSize,
        mode == FingerprintMode.FCFP ? FingerprintEngine.Mode.FCFP : FingerprintEngine.Mode.ECFP);
  }

  public static long[] circularFingerprint(MolGraph g, int radius, int fpSize, FingerprintMode mode) {
    return FingerprintEngine.circularFingerprint(g, radius, fpSize,
        mode == FingerprintMode.FCFP ? FingerprintEngine.Mode.FCFP : FingerprintEngine.Mode.ECFP);
  }

  // ---- Count-based circular fingerprint delegates (v6.0.1) ----

  /** @see FingerprintEngine#ecfpCounts(IAtomContainer, int, int) */
  public static int[] circularFingerprintECFPCounts(IAtomContainer mol, int radius, int fpSize) {
    return FingerprintEngine.ecfpCounts(mol, radius, fpSize);
  }

  /** @see FingerprintEngine#fcfpCounts(IAtomContainer, int, int) */
  public static int[] circularFingerprintFCFPCounts(IAtomContainer mol, int radius, int fpSize) {
    return FingerprintEngine.fcfpCounts(mol, radius, fpSize);
  }

  /** @see FingerprintEngine#ecfpCounts(MolGraph, int, int) */
  public static int[] circularFingerprintECFPCounts(MolGraph g, int radius, int fpSize) {
    return FingerprintEngine.ecfpCounts(g, radius, fpSize);
  }

  /** @see FingerprintEngine#fcfpCounts(MolGraph, int, int) */
  public static int[] circularFingerprintFCFPCounts(MolGraph g, int radius, int fpSize) {
    return FingerprintEngine.fcfpCounts(g, radius, fpSize);
  }

  /** @see FingerprintEngine#circularFingerprintCounts(IAtomContainer, int, int, FingerprintEngine.Mode) */
  public static int[] circularFingerprintCounts(IAtomContainer mol, int radius, int fpSize, FingerprintMode mode) {
    return FingerprintEngine.circularFingerprintCounts(mol, radius, fpSize,
        mode == FingerprintMode.FCFP ? FingerprintEngine.Mode.FCFP : FingerprintEngine.Mode.ECFP);
  }

  /** @see FingerprintEngine#circularFingerprintCounts(MolGraph, int, int, FingerprintEngine.Mode) */
  public static int[] circularFingerprintCounts(MolGraph g, int radius, int fpSize, FingerprintMode mode) {
    return FingerprintEngine.circularFingerprintCounts(g, radius, fpSize,
        mode == FingerprintMode.FCFP ? FingerprintEngine.Mode.FCFP : FingerprintEngine.Mode.ECFP);
  }

  /** @see FingerprintEngine#tanimoto */
  public static double fingerprintTanimoto(long[] fp1, long[] fp2) {
    return FingerprintEngine.tanimoto(fp1, fp2);
  }

  /** @see FingerprintEngine#dice */
  public static double fingerprintDice(long[] fp1, long[] fp2) {
    return FingerprintEngine.dice(fp1, fp2);
  }

  /** @see FingerprintEngine#cosine */
  public static double fingerprintCosine(long[] fp1, long[] fp2) {
    return FingerprintEngine.cosine(fp1, fp2);
  }

  /** @see FingerprintEngine#soergel */
  public static double fingerprintSoergel(long[] fp1, long[] fp2) {
    return FingerprintEngine.soergel(fp1, fp2);
  }

  /** @see FingerprintEngine#tanimotoCounts */
  public static double fingerprintTanimotoCounts(int[] fp1, int[] fp2) {
    return FingerprintEngine.tanimotoCounts(fp1, fp2);
  }

  /** @see FingerprintEngine#diceCounts */
  public static double fingerprintDiceCounts(int[] fp1, int[] fp2) {
    return FingerprintEngine.diceCounts(fp1, fp2);
  }

  /** @see FingerprintEngine#cosineCounts */
  public static double fingerprintCosineCounts(int[] fp1, int[] fp2) {
    return FingerprintEngine.cosineCounts(fp1, fp2);
  }

  /** @see FingerprintEngine#toBitSet */
  public static java.util.BitSet toBitSet(long[] fp) { return FingerprintEngine.toBitSet(fp); }
  /** @see FingerprintEngine#fromBitSet */
  public static long[] fromBitSet(java.util.BitSet bs) { return FingerprintEngine.fromBitSet(bs); }
  /** @see FingerprintEngine#toHex */
  public static String toHex(long[] fp) { return FingerprintEngine.toHex(fp); }
  /** @see FingerprintEngine#fromHex */
  public static long[] fromHex(String hex) { return FingerprintEngine.fromHex(hex); }
  /** @see FingerprintEngine#toBinaryString */
  public static String toBinaryString(long[] fp, int fpSize) { return FingerprintEngine.toBinaryString(fp, fpSize); }

  // ---- Topological torsion fingerprint delegates (v6.1) ----

  /** @see FingerprintEngine#topologicalTorsion(IAtomContainer, int) */
  public static long[] topologicalTorsion(IAtomContainer mol, int fpSize) {
    return FingerprintEngine.topologicalTorsion(mol, fpSize);
  }

  /** @see FingerprintEngine#topologicalTorsion(MolGraph, int) */
  public static long[] topologicalTorsion(MolGraph g, int fpSize) {
    return FingerprintEngine.topologicalTorsion(g, fpSize);
  }

  /** @see FingerprintEngine#topologicalTorsionCounts(MolGraph, int) */
  public static int[] topologicalTorsionCounts(MolGraph g, int fpSize) {
    return FingerprintEngine.topologicalTorsionCounts(g, fpSize);
  }

  // ---- Upper-bound delegates (v6.1) ----

  /**
   * Label-frequency upper bound (LFUB) on MCS size.
   *
   * @see SearchEngine#labelFrequencyUpperBound
   */
  public static int labelFrequencyUpperBound(MolGraph g1, MolGraph g2, ChemOptions C) {
    return SearchEngine.labelFrequencyUpperBound(g1, g2, C);
  }

  /**
   * Degree Sequence Upper Bound (DSB) on MCS size.
   *
   * @see SearchEngine#degreeSequenceUpperBound
   */
  public static int degreeSequenceUpperBound(MolGraph g1, MolGraph g2, ChemOptions C) {
    return SearchEngine.degreeSequenceUpperBound(g1, g2, C);
  }

  // ---- Scaffold MCS (Murcko) ----

  /**
   * Extract the Murcko scaffold (ring systems + linkers) from a molecule.
   *
   * <pre>{@code
   * IAtomContainer scaffold = SMSD.murckoScaffold(atorvastatin);
   * System.out.println("Scaffold atoms: " + scaffold.getAtomCount());
   * }</pre>
   *
   * @param mol the molecule to extract the scaffold from
   * @return the Murcko scaffold as an {@link IAtomContainer}
   */
  public static IAtomContainer murckoScaffold(IAtomContainer mol) {
    return SearchEngine.murckoScaffold(mol);
  }

  /**
   * Find the MCS between Murcko scaffolds of two molecules.
   *
   * <pre>{@code
   * SearchEngine.MCSOptions mcsOpts = new SearchEngine.MCSOptions();
   * mcsOpts.timeoutMs = 10000;
   * Map<Integer, Integer> scaffoldMcs =
   *     SMSD.findScaffoldMCS(mol1, mol2, new ChemOptions(), mcsOpts);
   * }</pre>
   *
   * @param m1 the first molecule
   * @param m2 the second molecule
   * @param C  chemical matching options
   * @param M  MCS options (timeout, induced, connected, etc.)
   * @return mapping from scaffold atom indices (first molecule) to scaffold atom indices (second)
   */
  public static Map<Integer, Integer> findScaffoldMCS(
      IAtomContainer m1, IAtomContainer m2, ChemOptions C, SearchEngine.MCSOptions M) {
    return SearchEngine.findScaffoldMCS(m1, m2, C, M);
  }

  // ---- Atom-Atom Mapping for Reactions ----

  


  


  // ---- MCS SMILES ----

  /**
   * Find the MCS and return it as a canonical SMILES string.
   *
   * <pre>{@code
   * SMSD smsd = new SMSD(mol1, mol2, new ChemOptions());
   * String mcsSmi = smsd.findMcsSmiles();
   * System.out.println("MCS SMILES: " + mcsSmi);
   * }</pre>
   *
   * @return canonical SMILES of the MCS, or {@code ""} if no common substructure
   */
  public String findMcsSmiles() {
    return findMcsSmiles(false, true, mcsTimeoutMs);
  }

  /**
   * Find the MCS with full parameter control and return it as a canonical SMILES string.
   *
   * @param induced       {@code true} for induced MCS (edge-preserving)
   * @param connectedOnly {@code true} to require a connected MCS subgraph
   * @param timeoutMs     maximum computation time in milliseconds
   * @return canonical SMILES of the MCS, or {@code ""} if no common substructure
   */
  public String findMcsSmiles(boolean induced, boolean connectedOnly, long timeoutMs) {
    requireMoleculeQuery("findMcsSmiles");
    SearchEngine.MCSOptions opt = new SearchEngine.MCSOptions();
    opt.induced = induced;
    opt.connectedOnly = connectedOnly;
    opt.timeoutMs = timeoutMs;
    return SearchEngine.findMcsSmiles(query, target, chem, opt);
  }

  /**
   * Extract a canonical SMILES from an existing MCS atom-index mapping.
   *
   * @param mol     the molecule from which the mapping keys originate
   * @param mapping MCS mapping (query atom index to target atom index)
   * @return canonical SMILES of the MCS subgraph, or {@code ""} if the mapping is empty
   * @see SearchEngine#mcsToSmiles(MolGraph, Map)
   */
  public static String mcsToSmiles(IAtomContainer mol, Map<Integer, Integer> mapping) {
    if (mapping == null || mapping.isEmpty()) return "";
    return SearchEngine.mcsToSmiles(new MolGraph(mol), mapping);
  }

  // ---- Accessors ----

  /**
   * Return the standardised query molecule.
   *
   * @return the query molecule (never {@code null})
   * @throws UnsupportedOperationException if this instance was created with a SMARTS query
   */
  public IAtomContainer getQuery() {
    if (querySmarts != null)
      throw new UnsupportedOperationException("getQuery is not supported for SMARTS queries.");
    return query;
  }

  /**
   * Return the standardised target molecule.
   *
   * @return the target molecule (never {@code null})
   */
  public IAtomContainer getTarget() { return target; }

  // ---- Internal ----

  private static final java.util.logging.Logger LOG =
      java.util.logging.Logger.getLogger(SMSD.class.getName());

  private static IAtomContainer standardiseSafe(IAtomContainer mol) {
    try { return Standardiser.standardise(mol, Standardiser.TautomerMode.NONE); }
    catch (Exception e) {
      LOG.warning("Standardisation failed, using raw molecule: " + e.getMessage());
      return mol;
    }
  }

  private void requireMoleculeQuery(String method) {
    if (querySmarts != null)
      throw new UnsupportedOperationException(method + " is not supported for SMARTS queries.");
  }
}

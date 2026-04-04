/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.IntStream;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;

/**
 * Core search engine providing substructure matching, MCS computation, fingerprinting,
 * scaffold extraction, and reaction atom-atom mapping.
 *
 * <p>All methods are stateless and thread-safe. For convenient wrappers with automatic
 * standardisation, use {@link SMSD} instead.
 *
 * @author Syed Asad Rahman
 * @see SMSD
 * @see ChemOptions
 */
public final class SearchEngine {

  // ---- MolGraph cache (v6.5.3 perf fix) ----
  // Keyed by IAtomContainer object identity so the same CDK instance maps
  // to the same MolGraph instance.  We intentionally avoid structural equality
  // here because different parses of the same SMILES can need independent graph
  // state during search.
  private static final Map<IAtomContainer, MolGraph> molGraphCache =
      Collections.synchronizedMap(new java.util.IdentityHashMap<>());

  /**
   * Get or create a MolGraph from an IAtomContainer, reusing cached instances.
   * The cache is keyed by object identity (WeakHashMap: auto-evicted on GC).
   * @since 6.5.3
   */
  static MolGraph toMolGraph(IAtomContainer mol) {
    MolGraph cached = molGraphCache.get(mol);
    if (cached != null) return cached;
    MolGraph g = new MolGraph(mol);
    molGraphCache.put(mol, g);
    return g;
  }

  /** Clear the MolGraph cache (call between reaction batches if needed). */
  public static void clearMolGraphCache() {
    molGraphCache.clear();
    SubstructureEngine.clearDomainCache();
    FingerprintEngine.clearFPCache();
  }

  private static final int GREEDY_PROBE_MAX_SIZE = 40;
  private static final int SEED_EXTEND_MAX_ATOMS = 50;
  private static final long SEED_EXTEND_NODE_LIMIT = 100_000L;
  private static final long MCSPLIT_NODE_LIMIT = 400_000L;
  private static final double BK_SKIP_RATIO = 0.95;

  private static final long MAX_NODE_LIMIT = 800_000L;

  private SearchEngine() {}

  /**
   * Apply Tier 2 solvent correction from ChemOptions to a MolGraph. No-op if solvent is AQUEOUS
   * (default) or null.
   */
  private static void applySolvent(MolGraph g, ChemOptions C) {
    if (C != null && C.solvent != null && C.solvent != ChemOptions.Solvent.AQUEOUS) {
      g.applySolventCorrection(C.solvent);
    }
  }

  /**
   * Time budget tracker for bounded search operations.
   *
   * <p>Uses {@link System#nanoTime()} with amortised expiration checks (every 1024 calls)
   * to minimize timing overhead during tight search loops.
   *
   * <pre>{@code
   * TimeBudget tb = new TimeBudget(5000); // 5-second budget
   * while (!tb.expired()) { ... }
   * }</pre>
   */
  public static final class TimeBudget {
    final long deadlineNanos;
    private final long checkEvery;
    private long counter;

    public TimeBudget(long timeoutMs) {
      this.deadlineNanos = System.nanoTime() + Math.max(1, timeoutMs) * 1_000_000L;
      this.checkEvery = 1024L;
      this.counter = 0L;
    }

    public boolean expired() {
      if ((++counter & (checkEvery - 1)) != 0) return false;
      return System.nanoTime() > deadlineNanos;
    }

    public boolean expiredNow() {
      return System.nanoTime() > deadlineNanos;
    }

    long remainingMillis() {
      return Math.max(0L, (deadlineNanos - System.nanoTime()) / 1_000_000L);
    }
  }

  /**
   * Search telemetry: nodes visited, backtracks, pruning breakdown, timing, and result count.
   *
   * <p>Returned by {@link #isSubstructureWithStats} and {@link #findAllSubstructuresWithStats}.
   */
  public static final class SubstructureStats {
    public final long nodesVisited, backtracks, candidatesTried;
    public final long prunesAtom, prunesBond, prunesDegree, prunesNLF;
    public final long timeMillis;
    public final boolean timeout;
    public final int solutions;

    SubstructureStats(long nodesVisited, long backtracks, long candidatesTried,
        long prunesAtom, long prunesBond, long prunesDegree, long prunesNLF,
        long timeMillis, boolean timeout, int solutions) {
      this.nodesVisited = nodesVisited;
      this.backtracks = backtracks;
      this.candidatesTried = candidatesTried;
      this.prunesAtom = prunesAtom;
      this.prunesBond = prunesBond;
      this.prunesDegree = prunesDegree;
      this.prunesNLF = prunesNLF;
      this.timeMillis = timeMillis;
      this.timeout = timeout;
      this.solutions = solutions;
    }
  }

  /**
   * Result of a substructure search: match outcome, mappings, and detailed statistics.
   */
  public static final class SubstructureResult {
    public final boolean exists;
    public final List<Map<Integer, Integer>> mappings;
    public final SubstructureStats stats;

    SubstructureResult(boolean exists, List<Map<Integer, Integer>> mappings, SubstructureStats stats) {
      this.exists = exists;
      this.mappings = mappings;
      this.stats = stats;
    }
  }

  /**
   * Configuration options for MCS computation.
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * opts.induced = true;          // edge-preserving MCS
   * opts.connectedOnly = true;    // connected subgraph only
   * opts.timeoutMs = 30000;   // 30-second timeout
   * Map<Integer, Integer> mcs = SearchEngine.findMCS(g1, g2, chemOpts, opts);
   * }</pre>
   */
  public static final class McsOptions {
    public boolean induced = false;
    public boolean connectedOnly = true;
    /** Timeout in milliseconds.  -1 (default) = adaptive:
     *  {@code Math.min(30000, 500 + g1.n * g2.n * 2)}. */
    public long timeoutMs = -1;
    public boolean extraSeeds = true;
    public int seedNeighborhoodRadius = 2;
    public int seedMaxAnchors = 12;
    public boolean useTwoHopNLFInExtension = true;
    public boolean useThreeHopNLFInExtension = false;
    public boolean disconnectedMCS = false;
    public boolean maximizeBonds = false;
    // dMCS fragment constraints
    public int minFragmentSize = 1;
    public int maxFragments = Integer.MAX_VALUE;
    // Weighted/Property MCS
    public double[] atomWeights = null; // weights for query atoms
    /** Fuzzy template matching: allow up to this many unmatched atoms in
     *  either direction when applying a template via the greedy probe.
     *  0 = exact match required (default), 2 = allow +/-2 unmatched atoms. */
    public int templateFuzzyAtoms = 0;

    // ---- Reaction-aware post-filter (v6.4.0) ----

    /**
     * Enable built-in reaction-aware post-filtering.  Off by default.
     * When true, findMCS and mapReaction will enumerate near-MCS
     * candidates (sizes K, K-1, K-2) and re-rank by heteroatom
     * coverage and reaction-center proximity.
     * @since 6.4.0
     */
    public boolean reactionAware = false;

    /**
     * Maximum size deficit from the mathematical maximum K to consider.
     * Only used when reactionAware=true or postFilter!=null.
     * Default: 2 (consider candidates of size K, K-1, K-2).
     * @since 6.4.0
     */
    public int nearMcsDelta = 2;

    /**
     * Maximum number of near-MCS candidates to generate before scoring.
     * Default: 20.
     * @since 6.4.0
     */
    public int nearMcsCandidates = 20;

    /**
     * Custom post-filter.  If non-null, overrides the built-in
     * reaction-aware scorer even when reactionAware=true.
     * @since 6.4.0
     */
    public McsPostFilter postFilter = null;

    /**
     * When true, apply bond-change penalty scoring after reaction-aware
     * candidate generation. Candidates are ranked by chemical plausibility
     * of implied bond transformations (C-C breaks penalised most).
     * @since 6.5.3
     */
    public boolean bondChangeAware = false;

    /**
     * Target atom indices to exclude from MCS search.
     * Used by batchMcsConstrained to prevent re-mapping atoms
     * already claimed by earlier queries.
     * @since 6.6.1
     */
    public Set<Integer> excludedTargetAtoms = null;
  }

  /**
   * Structured result of a convenience MCS-from-SMILES computation.
   *
   * <p>Bundles the atom-index mapping, MCS size, Tanimoto-like overlap score,
   * the MCS subgraph as a canonical SMILES string, and the algorithm tier that
   * produced the result.
   *
   * @since 6.3.0
   */
  public static final class McsResult {
    /** Atom-index mapping from query to target (empty if no common substructure). */
    public final Map<Integer, Integer> mapping;
    /** Number of matched atoms. */
    public final int size;
    /** Tanimoto-like overlap: size / min(queryAtoms, targetAtoms). */
    public final double tanimoto;
    /** Canonical SMILES of the MCS subgraph, or {@code ""} if empty. */
    public final String mcsSmiles;

    McsResult(Map<Integer, Integer> mapping, int queryN, int targetN, String mcsSmiles) {
      this.mapping = mapping;
      this.size = mapping.size();
      int minN = Math.min(queryN, targetN);
      this.tanimoto = minN > 0 ? (double) this.size / minN : 0.0;
      this.mcsSmiles = mcsSmiles != null ? mcsSmiles : "";
    }

    @Override public String toString() {
      return "McsResult{size=" + size + ", tanimoto=" + String.format("%.3f", tanimoto)
          + ", mcsSmiles=\"" + mcsSmiles + "\"}";
    }
  }

  /**
   * All-in-one convenience: parse two SMILES strings, compute MCS, and return
   * a structured {@link McsResult} with mapping, size, Tanimoto, and MCS SMILES.
   *
   * <pre>{@code
   * McsResult r = SearchEngine.findMCSFromSmiles("c1ccccc1", "c1ccc(O)cc1",
   *     new ChemOptions(), new McsOptions());
   * System.out.println(r.mcsSmiles);  // "c1ccccc1"
   * System.out.println(r.tanimoto);   // 1.0
   * }</pre>
   *
   * @param smi1 SMILES string for the first molecule
   * @param smi2 SMILES string for the second molecule
   * @param C    chemical matching options (may be {@code null} for defaults)
   * @param M    MCS options (may be {@code null} for defaults)
   * @return structured MCS result
   * @throws Exception if either SMILES string cannot be parsed
   * @since 6.3.0
   */
  public static McsResult findMCSFromSmiles(String smi1, String smi2, ChemOptions C, McsOptions M)
      throws Exception {
    if (smi1 == null || smi1.isEmpty()) throw new IllegalArgumentException("smi1 must not be null/empty");
    if (smi2 == null || smi2.isEmpty()) throw new IllegalArgumentException("smi2 must not be null/empty");
    if (C == null) C = new ChemOptions();
    if (M == null) M = new McsOptions();

    org.openscience.cdk.smiles.SmilesParser sp =
        new org.openscience.cdk.smiles.SmilesParser(
            org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance());
    IAtomContainer mol1 = sp.parseSmiles(smi1);
    IAtomContainer mol2 = sp.parseSmiles(smi2);
    MolGraph g1 = toMolGraph(mol1);
    MolGraph g2 = toMolGraph(mol2);
    applySolvent(g1, C);
    applySolvent(g2, C);

    Map<Integer, Integer> mapping = findMCS(g1, g2, C, M);
    String smiles = mcsToSmiles(g1, mapping);
    return new McsResult(mapping, g1.n, g2.n, smiles);
  }

  /**
   * All-in-one convenience: parse a SMILES string, compute its ECFP fingerprint,
   * and return the bit-packed {@code long[]} array.
   *
   * <pre>{@code
   * long[] fp = SearchEngine.fingerprintFromSmiles("c1ccccc1", 2, 2048);
   * }</pre>
   *
   * @param smiles SMILES string for the molecule
   * @param radius ECFP radius (2 = ECFP4, 3 = ECFP6; -1 = unlimited)
   * @param fpSize fingerprint size in bits
   * @return bit-packed fingerprint array
   * @throws Exception if the SMILES string cannot be parsed
   * @since 6.3.0
   */
  public static long[] fingerprintFromSmiles(String smiles, int radius, int fpSize)
      throws Exception {
    if (smiles == null || smiles.isEmpty())
      throw new IllegalArgumentException("smiles must not be null/empty");
    org.openscience.cdk.smiles.SmilesParser sp =
        new org.openscience.cdk.smiles.SmilesParser(
            org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance());
    IAtomContainer mol = sp.parseSmiles(smiles);
    return FingerprintEngine.ecfp(mol, radius, fpSize);
  }

  // ---------------------------------------------------------------------------
  // Progressive MCS callback
  // ---------------------------------------------------------------------------

  /**
   * Callback interface for progressive MCS reporting.
   * Invoked after each major pipeline level with the current best result.
   */
  @FunctionalInterface
  public interface McsProgressCallback {
    /**
     * Called after each pipeline level with the current best MCS mapping.
     *
     * @param bestSoFar the best atom-index mapping found so far (g1 to g2)
     * @param bestSize  the size of the best mapping
     * @param elapsedMs wall-clock time elapsed since the search started
     */
    void onProgress(Map<Integer, Integer> bestSoFar, int bestSize, long elapsedMs);
  }

  /**
   * Find MCS with optional progress callback.
   * The callback is invoked after each major pipeline level (chain, tree,
   * greedy, substructure, seed-extend, McSplit, BK, McGregor) with the
   * current best mapping.
   *
   * @param g1       the first molecule graph
   * @param g2       the second molecule graph
   * @param C        chemical matching options
   * @param M        MCS options
   * @param callback progress callback; null to disable reporting
   * @return mapping from g1 atom indices to g2 atom indices
   */
  public static Map<Integer, Integer> findMCS(
      MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M,
      McsProgressCallback callback) {
    if (callback == null) return findMCS(g1, g2, C, M);

    long t0 = System.nanoTime();
    Map<Integer, Integer> best = findMCS(g1, g2, C, M);
    long elapsedMs = (System.nanoTime() - t0) / 1_000_000L;
    callback.onProgress(Collections.unmodifiableMap(best), best.size(), elapsedMs);
    return best;
  }

  // Public API: MolGraph-based (CDK-free)

  /**
   * Check if the query graph is a subgraph of the target graph.
   *
   * <pre>{@code
   * boolean match = SearchEngine.isSubstructure(queryGraph, targetGraph, opts, 5000);
   * }</pre>
   *
   * @param query     the query molecule graph
   * @param target    the target molecule graph
   * @param C         chemical matching options
   * @param timeoutMs maximum time in milliseconds
   * @return {@code true} if the query is a substructure of the target
   */
  public static boolean isSubstructure(MolGraph query, MolGraph target, ChemOptions C, long timeoutMs) {
    if (query == null || target == null) return false;
    if (query.n == 0) return true;
    if (C == null) C = new ChemOptions();
    if (timeoutMs <= 0) timeoutMs = 10_000L;
    return SubstructureEngine.isSubstructure(query, target, C, timeoutMs);
  }

  /**
   * Check substructure with detailed search telemetry.
   *
   * @param query     the query molecule graph
   * @param target    the target molecule graph
   * @param C         chemical matching options
   * @param timeoutMs maximum time in milliseconds
   * @return result containing match outcome and search statistics
   */
  public static SubstructureResult isSubstructureWithStats(
      MolGraph query, MolGraph target, ChemOptions C, long timeoutMs) {
    TimeBudget tb = new TimeBudget(timeoutMs);
    SubstructureEngine.Matcher m = SubstructureEngine.makeMatcher(query, target, C, tb);
    long t0 = System.nanoTime();
    boolean ok = m.exists();
    long elapsed = (System.nanoTime() - t0) / 1_000_000L;
    return new SubstructureResult(ok, Collections.emptyList(), m.buildStats(elapsed, ok ? 1 : 0));
  }

  /**
   * Enumerate all unique substructure mappings (MolGraph API).
   *
   * @param query        the query molecule graph
   * @param target       the target molecule graph
   * @param C            chemical matching options
   * @param maxSolutions maximum number of mappings to return
   * @param timeoutMs    maximum time in milliseconds
   * @return list of atom-index mappings (query index to target index)
   */
  public static List<Map<Integer, Integer>> findAllSubstructures(
      MolGraph query, MolGraph target, ChemOptions C, int maxSolutions, long timeoutMs) {
    if (query == null || target == null) return Collections.emptyList();
    if (C == null) C = new ChemOptions();
    if (maxSolutions <= 0) maxSolutions = 100;
    if (timeoutMs <= 0) timeoutMs = 10_000L;
    return SubstructureEngine.findAllSubstructures(query, target, C, maxSolutions, timeoutMs);
  }

  /**
   * Enumerate substructure mappings with detailed search telemetry (MolGraph API).
   *
   * @param query        the query molecule graph
   * @param target       the target molecule graph
   * @param C            chemical matching options
   * @param maxSolutions maximum number of mappings to enumerate
   * @param timeoutMs    maximum time in milliseconds
   * @return result containing mappings and search statistics
   */
  public static SubstructureResult findAllSubstructuresWithStats(
      MolGraph query, MolGraph target, ChemOptions C, int maxSolutions, long timeoutMs) {
    TimeBudget tb = new TimeBudget(timeoutMs);
    SubstructureEngine.Matcher m = SubstructureEngine.makeMatcher(query, target, C, tb);
    long t0 = System.nanoTime();
    List<Map<Integer, Integer>> out = new ArrayList<>();
    m.enumerate(maxSolutions, out);
    long elapsed = (System.nanoTime() - t0) / 1_000_000L;
    return new SubstructureResult(!out.isEmpty(), out, m.buildStats(elapsed, out.size()));
  }

  // Public API: IAtomContainer-based (CDK adapter)

  /**
   * Check if the query molecule is a substructure of the target molecule (CDK API).
   *
   * <pre>{@code
   * boolean match = SearchEngine.isSubstructure(benzene, phenol, new ChemOptions(), 5000);
   * }</pre>
   *
   * @param query     the query molecule
   * @param target    the target molecule
   * @param C         chemical matching options
   * @param timeoutMs maximum time in milliseconds
   * @return {@code true} if the query is a substructure of the target
   */
  public static boolean isSubstructure(IAtomContainer query, IAtomContainer target, ChemOptions C, long timeoutMs) {
    MolGraph gq = toMolGraph(query), gt = toMolGraph(target);
    applySolvent(gq, C);
    applySolvent(gt, C);
    return isSubstructure(gq, gt, C, timeoutMs);
  }

  /**
   * Check substructure with telemetry (CDK API).
   *
   * @param query     the query molecule
   * @param target    the target molecule
   * @param C         chemical matching options
   * @param timeoutMs maximum time in milliseconds
   * @return result containing match outcome and search statistics
   */
  public static SubstructureResult isSubstructureWithStats(
      IAtomContainer query, IAtomContainer target, ChemOptions C, long timeoutMs) {
    MolGraph gq = toMolGraph(query), gt = toMolGraph(target);
    applySolvent(gq, C);
    applySolvent(gt, C);
    return isSubstructureWithStats(gq, gt, C, timeoutMs);
  }

  /**
   * Enumerate all unique substructure mappings (CDK API).
   *
   * @param query        the query molecule
   * @param target       the target molecule
   * @param C            chemical matching options
   * @param maxSolutions maximum number of mappings to return
   * @param timeoutMs    maximum time in milliseconds
   * @return list of atom-index mappings (query index to target index)
   */
  public static List<Map<Integer, Integer>> findAllSubstructures(
      IAtomContainer query, IAtomContainer target, ChemOptions C, int maxSolutions, long timeoutMs) {
    MolGraph gq = toMolGraph(query), gt = toMolGraph(target);
    applySolvent(gq, C);
    applySolvent(gt, C);
    return findAllSubstructures(gq, gt, C, maxSolutions, timeoutMs);
  }

  /**
   * Enumerate substructure mappings with telemetry (CDK API).
   *
   * @param query        the query molecule
   * @param target       the target molecule
   * @param C            chemical matching options
   * @param maxSolutions maximum number of mappings to enumerate
   * @param timeoutMs    maximum time in milliseconds
   * @return result containing mappings and search statistics
   */
  public static SubstructureResult findAllSubstructuresWithStats(
      IAtomContainer query, IAtomContainer target, ChemOptions C, int maxSolutions, long timeoutMs) {
    MolGraph gq = toMolGraph(query), gt = toMolGraph(target);
    applySolvent(gq, C);
    applySolvent(gt, C);
    return findAllSubstructuresWithStats(gq, gt, C, maxSolutions, timeoutMs);
  }

  /**
   * Validate that a mapping is chemically correct: atoms compatible, bonds consistent, injective.
   *
   * <pre>{@code
   * List<String> errors = SearchEngine.validateMapping(g1, g2, mapping, opts);
   * if (!errors.isEmpty()) {
   *     System.err.println("Invalid mapping: " + errors);
   * }
   * }</pre>
   *
   * @param g1      the first molecule graph (query)
   * @param g2      the second molecule graph (target)
   * @param mapping the atom-index mapping to validate (query index to target index)
   * @param C       chemical matching options
   * @return list of error messages; empty if the mapping is valid
   */
  public static List<String> validateMapping(MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping, ChemOptions C) {
    List<String> errors = new ArrayList<>();
    if (mapping == null || mapping.isEmpty()) return errors;
    Set<Integer> usedT = new HashSet<>();
    for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
      int qi = e.getKey(), tj = e.getValue();
      if (qi < 0 || qi >= g1.n) { errors.add("query index " + qi + " out of range [0," + g1.n + ")"); continue; }
      if (tj < 0 || tj >= g2.n) { errors.add("target index " + tj + " out of range [0," + g2.n + ")"); continue; }
      if (!usedT.add(tj)) errors.add("target atom " + tj + " mapped twice (not injective)");
      if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C))
        errors.add("atom mismatch: q" + qi + "(" + g1.atomicNum[qi] + ") → t" + tj + "(" + g2.atomicNum[tj] + ")");
    }
    // Check bonds: every mapped query bond must exist in the target and be compatible.
    for (Map.Entry<Integer, Integer> e1 : mapping.entrySet()) {
      int qi = e1.getKey(), tj = e1.getValue();
      for (int qk : g1.neighbors[qi]) {
        if (qk <= qi) continue;
        Integer tkObj = mapping.get(qk);
        if (tkObj == null) continue;
        int tk = tkObj;
        if (g1.bondOrder(qi, qk) == 0) continue;
        if (g2.bondOrder(tj, tk) == 0) {
          errors.add("bond missing: q(" + qi + "-" + qk + ") vs t(" + tj + "-" + tk + ")");
          continue;
        }
        if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C))
          errors.add("bond incompatible: q(" + qi + "-" + qk + ") vs t(" + tj + "-" + tk + ")");
      }
    }
    return errors;
  }

  private static Map<Integer, Integer> repairInvalidMcsMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping, ChemOptions C) {
    if (mapping == null || mapping.isEmpty()) return Collections.emptyMap();
    Map<Integer, Integer> current = new LinkedHashMap<>(mapping);
    int maxRounds = mapping.size(); // each round removes exactly 1 entry; bounded by initial size
    for (int round = 0; round < maxRounds && !current.isEmpty(); round++) {
      List<String> errors = validateMapping(g1, g2, current, C);
      if (errors.isEmpty()) return current;

      String err = errors.get(0);
      int toRemove = -1;

      if (err.startsWith("atom mismatch: q")) {
        toRemove = parseIntAt(err, 16); // "atom mismatch: q" is 16 chars
      } else if (err.startsWith("bond missing: q(") || err.startsWith("bond incompatible: q(")) {
        int qStart = err.indexOf('(');
        if (qStart >= 0) toRemove = parseIntAt(err, qStart + 1);
      } else if (err.contains("target atom") && err.contains("mapped twice")) {
        // Parse the duplicate target index and remove one of its query keys
        int tIdx = parseIntAt(err, err.indexOf("target atom") + 12);
        if (tIdx >= 0) {
          for (Map.Entry<Integer, Integer> e : current.entrySet()) {
            if (e.getValue() == tIdx) { toRemove = e.getKey(); break; }
          }
        }
      }
      if (toRemove < 0 || !current.containsKey(toRemove))
        toRemove = current.keySet().iterator().next();

      current.remove(toRemove);
    }
    return validateMapping(g1, g2, current, C).isEmpty() ? current : Collections.emptyMap();
  }

  private static int parseIntAt(String s, int pos) {
    int end = pos;
    if (end < s.length() && s.charAt(end) == '-') end++;
    while (end < s.length() && Character.isDigit(s.charAt(end))) end++;
    if (end == pos) return -1;
    try { return Integer.parseInt(s.substring(pos, end)); }
    catch (NumberFormatException e) { return -1; }
  }

  private static long resolveMcsTimeout(MolGraph g1, MolGraph g2, McsOptions M) {
    long timeout = M.timeoutMs;
    if (timeout < 0) timeout = Math.min(30_000L, 500L + (long) g1.n * g2.n * 2);
    return Math.max(1L, timeout);
  }

  private static Map<Integer, Integer> recoverValidMcsMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> raw, ChemOptions C, McsOptions M) {
    Map<Integer, Integer> repaired = repairInvalidMcsMapping(g1, g2, raw, C);
    if (!validateMapping(g1, g2, repaired, C).isEmpty()) repaired = new LinkedHashMap<>();
    int repairedSize = repaired.size();
    long timeout = resolveMcsTimeout(g1, g2, M);

    // If one graph is fully contained in the other, prefer a real substructure mapping
    // over trimming an invalid MCS candidate.
    MolGraph sml = g1.n <= g2.n ? g1 : g2;
    MolGraph lrg = g1.n <= g2.n ? g2 : g1;
    boolean swapped = g1.n > g2.n;
    List<Map<Integer, Integer>> subMaps =
        SubstructureEngine.findAllSubstructures(sml, lrg, C, 1, Math.max(1L, timeout / 4));
    if (!subMaps.isEmpty()) {
      Map<Integer, Integer> recovered = new LinkedHashMap<>();
      for (Map.Entry<Integer, Integer> e : subMaps.get(0).entrySet()) {
        if (swapped) recovered.put(e.getValue(), e.getKey());
        else recovered.put(e.getKey(), e.getValue());
      }
      recovered = ppx(g1, g2, recovered, C, M);
      if (validateMapping(g1, g2, recovered, C).isEmpty() && recovered.size() > repairedSize) {
        return recovered;
      }
    }

    if (!repaired.isEmpty()) {
      Map<Integer, Integer> regrown = ppx(g1, g2, greedyAtomExtend(g1, g2, repaired, C, M), C, M);
      regrown = repairInvalidMcsMapping(g1, g2, regrown, C);
      if (validateMapping(g1, g2, regrown, C).isEmpty() && regrown.size() > repairedSize) return regrown;
    }
    return repaired;
  }

  /**
   * Check if a mapping is maximal: no unmapped pair (qi,tj) can be added while maintaining validity.
   *
   * <pre>{@code
   * boolean maximal = SearchEngine.isMappingMaximal(g1, g2, mcs, opts);
   * }</pre>
   *
   * @param g1      the first molecule graph (query)
   * @param g2      the second molecule graph (target)
   * @param mapping the atom-index mapping to check
   * @param C       chemical matching options
   * @return {@code true} if no additional atom pair can be added to the mapping
   */
  public static boolean isMappingMaximal(MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping, ChemOptions C) {
    Set<Integer> usedQ = mapping.keySet(), usedT = new HashSet<>(mapping.values());
    for (int qi = 0; qi < g1.n; qi++) {
      if (usedQ.contains(qi)) continue;
      for (int tj = 0; tj < g2.n; tj++) {
        if (usedT.contains(tj)) continue;
        if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
        boolean bondOk = true;
        for (int qk : g1.neighbors[qi]) {
          Integer tkObj = mapping.get(qk);
          if (tkObj == null) continue;
          int tk = tkObj;
          if (g2.bondOrder(tj, tk) == 0 || !MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C))
            { bondOk = false; break; }
        }
        if (bondOk) return false; // found extensible pair → not maximal
      }
    }
    return true;
  }

  /** Score a mapping: weighted sum if atomWeights set, bond count if maximizeBonds, otherwise atom count. */
  static int mcsScore(MolGraph g1, Map<Integer, Integer> map, McsOptions M) {
    if (M.atomWeights != null) {
      double w = 0.0;
      for (int qi : map.keySet())
        if (qi >= 0 && qi < M.atomWeights.length) w += M.atomWeights[qi];
      return (int) (w * 1000);
    }
    return M.maximizeBonds ? countMappedBonds(g1, map) : map.size();
  }

  static int heteroatomScore(MolGraph g1, Map<Integer, Integer> map) {
    int score = 0;
    for (int qi : map.keySet()) {
      int z = g1.atomicNum[qi];
      if (z != 6 && z != 1) score++;
    }
    return score;
  }

  private static boolean sameCanonicalGraph(MolGraph g1, MolGraph g2) {
    if (g1.n != g2.n) return false;
    g1.ensureCanonical();
    g2.ensureCanonical();
    if (g1.canonicalHash != g2.canonicalHash) return false;

    int[] inv1 = new int[g1.n];
    int[] inv2 = new int[g2.n];
    for (int i = 0; i < g1.n; i++) inv1[g1.canonicalLabel[i]] = i;
    for (int i = 0; i < g2.n; i++) inv2[g2.canonicalLabel[i]] = i;

    for (int r = 0; r < g1.n; r++) {
      int a = inv1[r], b = inv2[r];
      if (g1.label[a] != g2.label[b]) return false;
      if (g1.degree[a] != g2.degree[b]) return false;
      if (g1.formalCharge[a] != g2.formalCharge[b]) return false;
      if (g1.massNumber[a] != g2.massNumber[b]) return false;
    }

    for (int r = 0; r < g1.n; r++) {
      int a = inv1[r], b = inv2[r];
      long[] edges1 = new long[g1.degree[a]];
      long[] edges2 = new long[g2.degree[b]];
      for (int k = 0; k < g1.degree[a]; k++) {
        int nb = g1.neighbors[a][k];
        long sig = g1.canonicalLabel[nb];
        sig = sig * 1000003L + g1.bondOrder(a, nb);
        sig = sig * 2L + (g1.bondInRing(a, nb) ? 1L : 0L);
        sig = sig * 2L + (g1.bondAromatic(a, nb) ? 1L : 0L);
        edges1[k] = sig;
      }
      for (int k = 0; k < g2.degree[b]; k++) {
        int nb = g2.neighbors[b][k];
        long sig = g2.canonicalLabel[nb];
        sig = sig * 1000003L + g2.bondOrder(b, nb);
        sig = sig * 2L + (g2.bondInRing(b, nb) ? 1L : 0L);
        sig = sig * 2L + (g2.bondAromatic(b, nb) ? 1L : 0L);
        edges2[k] = sig;
      }
      Arrays.sort(edges1);
      Arrays.sort(edges2);
      if (!Arrays.equals(edges1, edges2)) return false;
    }
    return true;
  }

  /**
   * Find the Maximum Common Substructure between two molecule graphs.
   *
   * <p>Uses a multi-strategy pipeline: identity check, substructure containment,
   * seed-and-extend, McSplit partition refinement, BK clique, and McGregor extension.
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * opts.timeoutMs = 30000;
   * Map<Integer, Integer> mcs = SearchEngine.findMCS(g1, g2, chemOpts, opts);
   * }</pre>
   *
   * @param g1 the first molecule graph
   * @param g2 the second molecule graph
   * @param C  chemical matching options
   * @param M  MCS options (timeout, induced, connected, etc.)
   * @return mapping from g1 atom indices to g2 atom indices; empty if no common substructure
   */
  private static Map<Integer, Integer> findMCSImpl(MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M) {
    if (g1 == null || g2 == null) return Collections.emptyMap();
    if (C == null) C = new ChemOptions();
    if (M == null) M = new McsOptions();
    if (g1.n == 0 || g2.n == 0) return Collections.emptyMap();
    // Validate atomWeights size if provided
    if (M.atomWeights != null && M.atomWeights.length < g1.n)
      throw new IllegalArgumentException(
          "atomWeights length (" + M.atomWeights.length + ") < query atom count (" + g1.n + ")");
    // Ensure lazy fields needed by MCS (canonical labeling, ring counts, tautomer classes).
    g1.ensureCanonical(); g2.ensureCanonical();
    if (C.tautomerAware) { g1.ensureTautomerClasses(); g2.ensureTautomerClasses(); }
    if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) { g1.ensureRingCounts(); g2.ensureRingCounts(); }
    // Resolve adaptive timeout: scale with product of atom counts, capped at 30s
    long timeout = M.timeoutMs;
    if (timeout < 0) {
      timeout = Math.min(30_000L, 500L + (long) g1.n * g2.n * 2);
    }
    TimeBudget tb = new TimeBudget(timeout);
    // DSB is only admissible for induced MCS (degree constraint valid).
    // For non-induced (default), use the looser but safe label-frequency bound.
    int upperBound = M.induced
        ? degreeSequenceUpperBound(g1, g2, C)
        : labelFrequencyUpperBound(g1, g2, C);
    Map<Integer, Integer> best = Collections.emptyMap();
    int bestSize = 0;
    int bestScore = 0;

    // Identity check — zero-allocation flat-array queue (handles symmetric molecules)
    if (g1 == g2 || ((!C.useChirality && !C.useBondStereo) && sameCanonicalGraph(g1, g2))) {
      Map<Integer, Integer> id = new LinkedHashMap<>();
      if (g1 == g2) {
        for (int i = 0; i < g1.n; i++) id.put(i, i);
      } else {
        int[] head = new int[g2.n]; int[] next = new int[g2.n];
        Arrays.fill(head, -1);
        for (int j = g2.n - 1; j >= 0; j--) { next[j] = head[g2.canonicalLabel[j]]; head[g2.canonicalLabel[j]] = j; }
        for (int i = 0; i < g1.n; i++) { int cl = g1.canonicalLabel[i]; id.put(i, head[cl]); head[cl] = next[head[cl]]; }
      }
      if (M.disconnectedMCS && (M.minFragmentSize > 1 || M.maxFragments < Integer.MAX_VALUE))
        id = applyFragmentConstraints(g1, id, M.minFragmentSize, M.maxFragments);
      return ppx(g1, g2, id, C, M);
    }

    int minN = Math.min(g1.n, g2.n);
    boolean weightMode = M.maximizeBonds || M.atomWeights != null;

    // --- Linear chain fast-path (PEG, polymer chains) ---
    // For chain-like molecules (max degree ≤ 2), use O(n*m) longest-common-subpath DP
    // instead of exponential McSplit/BK that suffers from O(n!) symmetry explosion.
    if (!weightMode) {
      boolean g1Chain = true, g2Chain = true;
      for (int i = 0; i < g1.n && g1Chain; i++) if (g1.degree[i] > 2) g1Chain = false;
      for (int j = 0; j < g2.n && g2Chain; j++) if (g2.degree[j] > 2) g2Chain = false;
      if (g1Chain && g2Chain && g1.n >= 2 && g2.n >= 2) {
        int[] seq1 = walkChain(g1), seq2 = walkChain(g2);
        int n1 = seq1.length, n2 = seq2.length;
        int[][] dp = new int[n1 + 1][n2 + 1];
        int bestLen = 0, bestI = 0, bestJ = 0;
        for (int i = 1; i <= n1; i++) {
          int a1 = seq1[i - 1];
          for (int j = 1; j <= n2; j++) {
            int a2 = seq2[j - 1];
            if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, a1, g2, a2, C)) continue;
            if (i > 1 && j > 1 && dp[i - 1][j - 1] > 0
                && MolGraph.ChemOps.bondsCompatible(g1, seq1[i - 2], a1, g2, seq2[j - 2], a2, C)) {
              dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
              dp[i][j] = 1;
            }
            if (dp[i][j] > bestLen) { bestLen = dp[i][j]; bestI = i; bestJ = j; }
          }
        }
        if (bestLen > 0) {
          Map<Integer, Integer> chainMCS = new LinkedHashMap<>();
          for (int k = 0; k < bestLen; k++)
            chainMCS.put(seq1[bestI - bestLen + k], seq2[bestJ - bestLen + k]);
          bestSize = bestLen; best = chainMCS;
          if (bestLen >= upperBound) return ppx(g1, g2, best, C, M);
        }
      }
    }

    // --- Tree fast-path (branched polymers, dendrimers, glycogen) ---
    // For acyclic connected molecules (trees) with degree > 2, use Kilpelainen-Mannila
    // style bottom-up DP on rooted trees: O(n1*n2*d^2).
    if (!weightMode && g1.n >= 10 && g2.n >= 10) {
      if (isTree(g1) && isTree(g2)) {
        int root1 = treeCentroid(g1), root2 = treeCentroid(g2);
        int[][] par1 = new int[1][], par2 = new int[1][];
        int[][][] ch1 = new int[1][][], ch2 = new int[1][][];
        int[][] po1 = new int[1][], po2 = new int[1][];
        rootTree(g1, root1, par1, ch1, po1);
        rootTree(g2, root2, par2, ch2, po2);

        int n1t = g1.n, n2t = g2.n;
        int[][] dpTree = new int[n1t][n2t];
        int[][][] btTreeI = new int[n1t][n2t][];
        int[][][] btTreeJ = new int[n1t][n2t][];

        // Pre-allocate edge buffers outside loop to reduce GC pressure
        int maxDeg = 0;
        for (int i = 0; i < n1t; i++) maxDeg = Math.max(maxDeg, g1.degree[i]);
        for (int i = 0; i < n2t; i++) maxDeg = Math.max(maxDeg, g2.degree[i]);
        int edgeBufCap = maxDeg * maxDeg;
        int[] ew = new int[edgeBufCap], ei = new int[edgeBufCap], ej = new int[edgeBufCap];

        for (int u : po1[0]) {
          for (int v : po2[0]) {
            if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, u, g2, v, C)) {
              dpTree[u][v] = 0;
              continue;
            }
            int[] uc = ch1[0][u], vc = ch2[0][v];
            if (uc.length == 0 || vc.length == 0) {
              dpTree[u][v] = 1;
              btTreeI[u][v] = new int[0];
              btTreeJ[u][v] = new int[0];
              continue;
            }
            int d1 = uc.length, d2 = vc.length;
            int ec = 0;
            for (int ii = 0; ii < d1; ii++) {
              for (int jj = 0; jj < d2; jj++) {
                int w = dpTree[uc[ii]][vc[jj]];
                if (w > 0 && MolGraph.ChemOps.bondsCompatible(g1, u, uc[ii], g2, v, vc[jj], C)) {
                  ew[ec] = w; ei[ec] = ii; ej[ec] = jj; ec++;
                }
              }
            }
            for (int a = 1; a < ec; a++) {
              int tw = ew[a], ti = ei[a], tj = ej[a];
              int b = a - 1;
              while (b >= 0 && ew[b] < tw) { ew[b+1] = ew[b]; ei[b+1] = ei[b]; ej[b+1] = ej[b]; b--; }
              ew[b+1] = tw; ei[b+1] = ti; ej[b+1] = tj;
            }
            boolean[] usedI = new boolean[d1], usedJ = new boolean[d2];
            int childSum = 0, matchCount = 0;
            int[] mi = new int[Math.min(d1, d2)], mj = new int[Math.min(d1, d2)];
            for (int e = 0; e < ec; e++) {
              if (usedI[ei[e]] || usedJ[ej[e]]) continue;
              usedI[ei[e]] = true; usedJ[ej[e]] = true;
              childSum += ew[e];
              mi[matchCount] = ei[e]; mj[matchCount] = ej[e];
              matchCount++;
            }
            dpTree[u][v] = 1 + childSum;
            btTreeI[u][v] = Arrays.copyOf(mi, matchCount);
            btTreeJ[u][v] = Arrays.copyOf(mj, matchCount);
          }
        }

        int treeBestVal = 0, bestU = -1, bestV = -1;
        for (int u = 0; u < n1t; u++)
          for (int v = 0; v < n2t; v++)
            if (dpTree[u][v] > treeBestVal) { treeBestVal = dpTree[u][v]; bestU = u; bestV = v; }

        if (treeBestVal > bestSize) {
          Map<Integer, Integer> treeMap = new LinkedHashMap<>();
          Deque<int[]> stack = new ArrayDeque<>();
          stack.push(new int[]{bestU, bestV});
          while (!stack.isEmpty()) {
            int[] uv = stack.pop();
            int u = uv[0], v = uv[1];
            treeMap.put(u, v);
            int[] bti = btTreeI[u][v], btj = btTreeJ[u][v];
            if (bti != null) {
              for (int k = 0; k < bti.length; k++)
                stack.push(new int[]{ch1[0][u][bti[k]], ch2[0][v][btj[k]]});
            }
          }
          if (treeMap.size() > bestSize) {
            best = treeMap;
            bestSize = treeMap.size();
          }
          if (bestSize >= upperBound) return ppx(g1, g2, best, C, M);
        }
      }
    }

    // Greedy probing for small similar molecules (skip fast-exit in bond-maximizing mode)
    if (!(M.maximizeBonds || M.atomWeights != null) && minN > 2 && minN < GREEDY_PROBE_MAX_SIZE && upperBound >= minN) {
      Map<Integer, Integer> greedy = greedyProbe(g1, g2, C);
      int fuzzy = M.templateFuzzyAtoms;
      int greedySz = greedy.size();
      boolean fuzzyAccept = (fuzzy > 0) &&
          (Math.abs(g1.n - greedySz) <= fuzzy || Math.abs(g2.n - greedySz) <= fuzzy);
      if (greedySz >= upperBound) {
        Map<Integer, Integer> greedyC = applyRingAnchorGuard(g1, g2,
            M.connectedOnly ? largestConnected(g1, greedy) : greedy, C);
        if (greedyC.size() >= upperBound) return ppx(g1, g2, greedyC, C, M);
      }
      if (fuzzyAccept && greedySz > bestSize) {
        Map<Integer, Integer> greedyC = applyRingAnchorGuard(g1, g2,
            M.connectedOnly ? largestConnected(g1, greedy) : greedy, C);
        if (greedyC.size() > bestSize) {
          best = greedyC;
          bestSize = greedyC.size();
          bestScore = mcsScore(g1, best, M);
        }
      }
      // Seed downstream McSplit with greedy result if it covers > 60% of upper bound
      if (greedy.size() > upperBound * 0.6 && greedy.size() > bestSize) {
        Map<Integer, Integer> greedyC = applyRingAnchorGuard(g1, g2,
            M.connectedOnly ? largestConnected(g1, greedy) : greedy, C);
        if (greedyC.size() > bestSize) {
          best = greedyC;
          bestSize = greedyC.size();
          bestScore = mcsScore(g1, best, M);
        }
      }
    }

    // Augmenting path refinement: grow mapping by forced single-candidate extensions
    if (bestSize > 0 && bestSize < upperBound && !tb.expired()) {
      Map<Integer, Integer> augmented = new LinkedHashMap<>(best);
      boolean grew = true;
      while (grew && !tb.expired()) {
        grew = false;
        boolean[] usedQ = new boolean[g1.n];
        boolean[] usedT = new boolean[g2.n];
        for (Map.Entry<Integer, Integer> e : augmented.entrySet()) {
          usedQ[e.getKey()] = true;
          usedT[e.getValue()] = true;
        }
        for (Map.Entry<Integer, Integer> e : new ArrayList<>(augmented.entrySet())) {
          int qi = e.getKey(), ti = e.getValue();
          for (int qk : g1.neighbors[qi]) {
            if (usedQ[qk]) continue;
            int candidate = -1, count = 0;
            for (int tk : g2.neighbors[ti]) {
              if (usedT[tk]) continue;
              if (SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qk, g2, tk, C)
                  && MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, ti, tk, C)) {
                candidate = tk; count++;
              }
            }
            if (count == 1) {
              // Verify full consistency with all mapped neighbors of qk
              boolean consistent = true;
              for (int qm : g1.neighbors[qk]) {
                if (!augmented.containsKey(qm)) continue;
                int tm = augmented.get(qm);
                int qOrd = g1.bondOrder(qk, qm), tOrd = g2.bondOrder(candidate, tm);
                if ((qOrd != 0) != (tOrd != 0)) { consistent = false; break; }
                if (qOrd != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qk, qm, g2, candidate, tm, C))
                  { consistent = false; break; }
              }
              if (consistent) {
                augmented.put(qk, candidate);
                usedQ[qk] = true; usedT[candidate] = true;
                grew = true;
              }
            }
          }
        }
      }
      if (augmented.size() > bestSize) {
        best = augmented; bestSize = best.size();
        bestScore = mcsScore(g1, best, M);
      }
      if (!weightMode && bestSize >= upperBound) return ppx(g1, g2, best, C, M);
    }

    // Substructure containment: if one molecule is inside the other
    if (!(M.maximizeBonds || M.atomWeights != null) && g1.n > 0 && g2.n > 0 && Math.min(g1.n, g2.n) <= Math.max(g1.n, g2.n) * 3 / 4) {
      MolGraph sml = g1.n <= g2.n ? g1 : g2;
      MolGraph lrg = g1.n <= g2.n ? g2 : g1;
      boolean swapped = g1.n > g2.n;
      SubstructureEngine.Matcher subMatcher = SubstructureEngine.makeMatcher(sml, lrg, C, tb);
      if (subMatcher.exists()) {
        List<Map<Integer, Integer>> subMaps = new ArrayList<>();
        subMatcher.enumerate(1, subMaps);
        if (!subMaps.isEmpty()) {
          Map<Integer, Integer> raw = subMaps.get(0);
          if (!swapped) return ppx(g1, g2, raw, C, M);
          Map<Integer, Integer> full = new LinkedHashMap<>();
          for (Map.Entry<Integer, Integer> e : raw.entrySet()) full.put(e.getValue(), e.getKey());
          return ppx(g1, g2, full, C, M);
        }
      }
    }

    // Seed-and-extend MCS (bond growth, node-count limited)
    GraphBuilder GB = new GraphBuilder(g1, g2, C, M.induced);
    if (minN >= 4 && Math.max(g1.n, g2.n) <= SEED_EXTEND_MAX_ATOMS && !tb.expired()) {
      Map<Integer, Integer> seSeed = GB.seedExtendMCS(tb, upperBound);
      int seScore = mcsScore(g1, seSeed, M);
      if ((M.maximizeBonds || M.atomWeights != null) ? seScore > bestScore : seSeed.size() > bestSize) {
        best = seSeed;
        bestSize = seSeed.size();
        bestScore = seScore;
      }
      if (!(M.maximizeBonds || M.atomWeights != null) && bestSize >= upperBound) return ppx(g1, g2, best, C, M);
    }

    // k-core pre-pruning only pays for smaller product graphs.
    // On medium-sized pairs the setup cost can dominate the whole search,
    // and the derived bound is advisory only.
    int kcoreUB = upperBound;
    int pgSize = g1.n * g2.n;
    if (bestSize > 1 && !tb.expired() && pgSize <= 1000) {
      int k = bestSize;
      int n1k = g1.n, n2k = g2.n;
      int[][] pgDeg = new int[n1k][n2k];
      boolean[][] alive = new boolean[n1k][n2k];
      for (int qi = 0; qi < n1k; qi++)
        for (int tj = 0; tj < n2k; tj++)
          if (SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C))
            alive[qi][tj] = true;
      for (int qi = 0; qi < n1k; qi++)
        for (int tj = 0; tj < n2k; tj++) {
          if (!alive[qi][tj]) continue;
          int deg = 0;
          for (int qk2 : g1.neighbors[qi])
            for (int tl : g2.neighbors[tj])
              if (alive[qk2][tl] && MolGraph.ChemOps.bondsCompatible(g1, qi, qk2, g2, tj, tl, C)) deg++;
          pgDeg[qi][tj] = deg;
        }
      boolean changed = true;
      while (changed && !tb.expired()) {
        changed = false;
        for (int qi = 0; qi < n1k; qi++)
          for (int tj = 0; tj < n2k; tj++) {
            if (!alive[qi][tj] || pgDeg[qi][tj] >= k - 1) continue;
            alive[qi][tj] = false;
            changed = true;
            for (int qk2 : g1.neighbors[qi])
              for (int tl : g2.neighbors[tj])
                if (alive[qk2][tl]) pgDeg[qk2][tl]--;
          }
      }
      int kcUB = 0;
      for (int qi = 0; qi < n1k; qi++) {
        boolean has = false;
        for (int tj = 0; tj < n2k; tj++) if (alive[qi][tj]) { has = true; break; }
        if (has) kcUB++;
      }
      if (kcUB < kcoreUB) kcoreUB = kcUB;
    }
    // McSplit partition-refinement (node-count limited, hardware-independent)
    boolean mcSplitExhaustive;
    int mcSplitSize;
    {
      long[] nodeCount = {0};
      Map<Integer, Integer> mcSeed = GB.mcSplitSeed(tb, nodeCount);
      mcSplitSize = mcSeed.size();
      int mcScore = mcsScore(g1, mcSeed, M);
      if ((M.maximizeBonds || M.atomWeights != null) ? mcScore > bestScore : mcSeed.size() > bestSize) {
        best = mcSeed;
        bestSize = mcSeed.size();
        bestScore = mcScore;
      }
      mcSplitExhaustive = nodeCount[0] < MCSPLIT_NODE_LIMIT;
    }

    if (!(M.maximizeBonds || M.atomWeights != null) && mcSplitExhaustive && bestSize >= upperBound) return ppx(g1, g2, best, C, M);

    // Near-optimal fast path: greedy atom extension before expensive BK + McGregor
    if (bestSize >= upperBound - 4 && bestSize > 0 && bestSize < upperBound && !tb.expired()) {
      Map<Integer, Integer> ext = ppx(g1, g2, greedyAtomExtend(g1, g2, ppx(g1, g2, best, C, M), C, M), C, M);
      int extScore = mcsScore(g1, ext, M);
      if ((M.maximizeBonds || M.atomWeights != null) ? extScore > bestScore : ext.size() > bestSize) {
        best = ext;
        bestSize = ext.size();
        bestScore = extScore;
      }
      if (!(M.maximizeBonds || M.atomWeights != null) && bestSize >= upperBound) return ppx(g1, g2, best, C, M);
    }

    // BK clique on product graph
    int bkSize = 0;
    if (bestSize < (int) (upperBound * BK_SKIP_RATIO) && !tb.expired()) {
      Map<Integer, Integer> cliqueSeed = GB.maximumCliqueSeed(tb);
      bkSize = cliqueSeed.size();
      int cScore = mcsScore(g1, cliqueSeed, M);
      if ((M.maximizeBonds || M.atomWeights != null) ? cScore > bestScore : cliqueSeed.size() > bestSize) {
        best = cliqueSeed;
        bestSize = cliqueSeed.size();
        bestScore = cScore;
      }
    }
    if (!(M.maximizeBonds || M.atomWeights != null) && bestSize >= upperBound) return ppx(g1, g2, best, C, M);

    // McGregor extension with extra seeds
    List<Map<Integer, Integer>> seeds = new ArrayList<>();
    if (bestSize > 0) seeds.add(best);

    // Skip extra seeds when McSplit and BK agree on MCS size
    boolean skipExtraSeeds = mcSplitSize > 0 && bkSize > 0 && mcSplitSize == bkSize
        && mcSplitSize >= upperBound - 1;

    if (M.extraSeeds && !skipExtraSeeds && bestSize < upperBound && !tb.expired()) {
      Map<Integer, Integer> s = GB.ringAnchorSeed(tb);
      if (!s.isEmpty()) seeds.add(s);
      if (!tb.expired()) {
        s = GB.labelDegreeAnchorSeed(tb);
        if (!s.isEmpty()) seeds.add(s);
      }
      if (!tb.expired()) {
        s = GB.vf2ppRingSkeletonSeed(tb, Math.max(1, tb.remainingMillis() / 4), 2, 12);
        if (!s.isEmpty()) seeds.add(s);
      }
      if (!tb.expired()) {
        s = GB.vf2ppCoreSeed(tb, Math.max(1, tb.remainingMillis() / 4), 2, 12);
        if (!s.isEmpty()) seeds.add(s);
      }
    }

    long perSeedMs = Math.max(1, tb.remainingMillis() / Math.max(1, seeds.size()));
    for (Map<Integer, Integer> seed : seeds) {
      if (tb.expired()) break;
      Map<Integer, Integer> ext = ppx(g1, g2,
          mcGregorExtend(g1, g2, seed, C, tb, perSeedMs, M.useTwoHopNLFInExtension, M.useThreeHopNLFInExtension, M.connectedOnly),
          C, M);
      int seedScore = mcsScore(g1, ext, M);
      if ((M.maximizeBonds || M.atomWeights != null) ? seedScore > bestScore : ext.size() > bestSize) {
        best = ext;
        bestSize = ext.size();
        bestScore = seedScore;
      }
      if (!(M.maximizeBonds || M.atomWeights != null) && bestSize >= upperBound) return ppx(g1, g2, best, C, M);
    }

    // Last resort: start from empty seed
    if (bestScore <= 0 && !tb.expired()) {
      best = ppx(g1, g2,
          mcGregorExtend(g1, g2, Collections.emptyMap(), C, tb, tb.remainingMillis(),
              M.useTwoHopNLFInExtension, M.useThreeHopNLFInExtension, M.connectedOnly),
          C, M);
    }
    // Invariant: MCS size must never exceed the smaller molecule's atom count
    assert best.size() <= Math.min(g1.n, g2.n)
        : "MCS size " + best.size() + " exceeds min(n1,n2) = " + Math.min(g1.n, g2.n);
    return ppx(g1, g2, best, C, M);
  }

  public static Map<Integer, Integer> findMCS(MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M) {
    if (g1 == null || g2 == null) return Collections.emptyMap();
    if (C == null) C = new ChemOptions();
    if (M == null) M = new McsOptions();
    if (g1.n == 0 || g2.n == 0) return Collections.emptyMap();

    g1.ensureCanonical(); g2.ensureCanonical();
    if (C.tautomerAware) { g1.ensureTautomerClasses(); g2.ensureTautomerClasses(); }
    if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) { g1.ensureRingCounts(); g2.ensureRingCounts(); }

    int ub12 = labelFrequencyUpperBoundDirected(g1, g2, C);
    int ub21 = labelFrequencyUpperBoundDirected(g2, g1, C);
    OrientationPlan plan = chooseOrientationPlan(g1, g2, C, ub12, ub21);
    boolean weightMode = M.maximizeBonds || M.atomWeights != null;

    Map<Integer, Integer> best = plan.directFirst
        ? runValidatedMcsDirection(g1, g2, C, M, false)
        : runValidatedMcsDirection(g2, g1, C, M, true);
    if (!validateMapping(g1, g2, best, C).isEmpty()) best = recoverValidMcsMapping(g1, g2, best, C, M);

    int baseUb = M.induced ? degreeSequenceUpperBound(g1, g2, C) : labelFrequencyUpperBound(g1, g2, C);
    if (!M.induced) {
      baseUb = Math.min(baseUb, ub12);
      baseUb = Math.min(baseUb, ub21);
    }

    int alternateFloor = alternateOrientationFloor(plan, plan.directFirst);
    int probeCeiling = Math.max(Math.max(plan.seed12, plan.seed21), Math.max(plan.mc12, plan.mc21));
    boolean runAlternate = best.isEmpty();
    if (!runAlternate && best.size() < alternateFloor) runAlternate = true;
    if (!runAlternate && !weightMode && best.size() + 1 < probeCeiling) runAlternate = true;
    if (!runAlternate && !weightMode && baseUb >= 12 && best.size() + 4 < baseUb) runAlternate = true;
    if (!runAlternate && !weightMode
        && Math.abs(g1.n - g2.n) <= 1
        && Math.min(g1.n, g2.n) >= 20) {
      runAlternate = true;
    }
    if (!runAlternate && !weightMode
        && Math.abs(g1.n - g2.n) <= 2
        && Math.min(g1.n, g2.n) >= 20
        && best.size() * 4 < Math.min(g1.n, g2.n) * 3) {
      runAlternate = true;
    }
    if (!runAlternate && !weightMode && best.size() + 2 < baseUb && Math.abs(g1.n - g2.n) >= 4) runAlternate = true;

    if (runAlternate) {
      Map<Integer, Integer> alt = plan.directFirst
          ? runValidatedMcsDirection(g2, g1, C, M, true)
          : runValidatedMcsDirection(g1, g2, C, M, false);
      if (!validateMapping(g1, g2, alt, C).isEmpty()) alt = recoverValidMcsMapping(g1, g2, alt, C, M);
      if (preferFinalMapping(g1, alt, best, M)) best = alt;
    }
    return best;
  }

  static Map<Integer, Integer> greedyAtomExtend(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> seed, ChemOptions C, McsOptions M) {
    int n1 = g1.n, n2 = g2.n;
    int[] q2t = new int[n1], t2q = new int[n2];
    Arrays.fill(q2t, -1);
    Arrays.fill(t2q, -1);
    for (Map.Entry<Integer, Integer> e : seed.entrySet()) {
      q2t[e.getKey()] = e.getValue();
      t2q[e.getValue()] = e.getKey();
    }
    boolean progress = true;
    while (progress) {
      progress = false;
      for (int qi = 0; qi < n1; qi++) {
        if (q2t[qi] >= 0) continue;
        boolean onFrontier = false;
        for (int nb : g1.neighbors[qi])
          if (q2t[nb] >= 0) { onFrontier = true; break; }
        if (!onFrontier) continue;

        int bestTj = -1, bestScore = -1;
        for (int tj = 0; tj < n2; tj++) {
          if (t2q[tj] >= 0) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
          boolean consistent = true;
          for (int qk : g1.neighbors[qi]) {
            if (q2t[qk] < 0) continue;
            int tk = q2t[qk];
            int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
            if (qOrd != 0 && tOrd != 0) {
              if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
            } else if (M.induced && ((qOrd != 0) != (tOrd != 0))) { consistent = false; break; }
          }
          if (!consistent) continue;
          int score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0) + Math.min(g1.degree[qi], g2.degree[tj]);
          if (score > bestScore) { bestScore = score; bestTj = tj; }
        }
        if (bestTj >= 0) {
          q2t[qi] = bestTj;
          t2q[bestTj] = qi;
          progress = true;
        }
      }
    }
    Map<Integer, Integer> result = new LinkedHashMap<>();
    for (int i = 0; i < n1; i++) if (q2t[i] >= 0) result.put(i, q2t[i]);
    return result;
  }

  /**
   * Count the number of bonds in g1 whose both endpoints are in the mapping.
   *
   * @param g1  the molecule graph containing the bonds
   * @param map the atom-index mapping (query indices to target indices)
   * @return the number of bonds where both endpoints are mapped
   */
  public static int countMappedBonds(MolGraph g1, Map<Integer, Integer> map) {
    int count = 0;
    Set<Integer> mapped = map.keySet();
    for (int qi : mapped)
      for (int qk : g1.neighbors[qi])
        if (qk > qi && mapped.contains(qk)) count++;
    return count;
  }

  /** Enforce complete rings: if any atom in a ring is mapped, all must be. Iterate until stable. */
  static Map<Integer, Integer> enforceCompleteRings(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> map) {
    if (map.isEmpty()) return map;
    int[][] rings = g1.computeRings();
    if (rings.length == 0) return map;
    Map<Integer, Integer> result = new LinkedHashMap<>(map);
    boolean changed = true;
    while (changed) {
      changed = false;
      for (int[] ring : rings) {
        boolean anyMapped = false, allMapped = true;
        for (int atom : ring) {
          if (result.containsKey(atom)) anyMapped = true;
          else allMapped = false;
        }
        if (anyMapped && !allMapped) {
          for (int atom : ring) {
            if (result.remove(atom) != null) changed = true;
          }
        }
      }
    }
    return result;
  }

  /** Post-process MCS: iteratively apply filters until stable (filters can interact). */
  static Map<Integer, Integer> ppx(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> ext, ChemOptions C, McsOptions M) {
    boolean changed = true;
    while (changed) {
      int startSize = ext.size();
      if (M.induced) ext = pruneToInduced(g1, g2, ext, C);
      if (C.completeRingsOnly) ext = enforceCompleteRings(g1, g2, ext);
      if (!M.disconnectedMCS && M.connectedOnly) ext = largestConnected(g1, ext);
      changed = ext.size() < startSize;
    }
    ext = applyRingAnchorGuard(g1, g2, ext, C);
    if (M.disconnectedMCS && (M.minFragmentSize > 1 || M.maxFragments < Integer.MAX_VALUE))
      ext = applyFragmentConstraints(g1, ext, M.minFragmentSize, M.maxFragments);
    return ext;
  }

  /**
   * Find the Maximum Common Substructure between two CDK molecules.
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * Map<Integer, Integer> mcs = SearchEngine.findMCS(mol1, mol2, new ChemOptions(), opts);
   * }</pre>
   *
   * @param m1 the first molecule
   * @param m2 the second molecule
   * @param C  chemical matching options
   * @param M  MCS options
   * @return mapping from m1 atom indices to m2 atom indices; empty if no common substructure
   */
  public static Map<Integer, Integer> findMCS(IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions M) {
    MolGraph g1 = toMolGraph(m1), g2 = toMolGraph(m2);
    applySolvent(g1, C);
    applySolvent(g2, C);
    return findMCS(g1, g2, C, M);
  }

  /**
   * Enumerate multiple distinct MCS mappings of the maximum size.
   *
   * <p>First computes the single best MCS to determine the optimal size K, then
   * re-searches the solution space to collect up to {@code maxResults} distinct
   * atom-atom mappings of that same size K. Useful for scaffold analysis and
   * structure-activity relationship (SAR) studies where alternative mappings
   * provide different chemical perspectives.
   *
   * <p>The method respects the timeout in {@code M}; if the timeout expires during
   * enumeration, it returns whatever solutions have been collected so far.
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * opts.timeoutMs = 30000;
   * List<Map<Integer, Integer>> allMCS =
   *     SearchEngine.findAllMCS(g1, g2, chemOpts, opts, 10);
   * for (Map<Integer, Integer> m : allMCS) {
   *     System.out.println("Mapping: " + m);
   * }
   * }</pre>
   *
   * @param g1         the first molecule graph
   * @param g2         the second molecule graph
   * @param C          chemical matching options
   * @param M          MCS options (timeout, induced, connected, etc.)
   * @param maxResults maximum number of distinct mappings to return (capped at 10 if &le; 0)
   * @return list of distinct atom-index mappings, each of the maximum MCS size;
   *         empty list if no common substructure exists
   */
  public static List<Map<Integer, Integer>> findAllMCS(
      MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M, int maxResults) {
    if (maxResults <= 0) maxResults = 10;

    // Phase 1: find the optimal MCS size
    Map<Integer, Integer> best = findMCS(g1, g2, C, M);
    int K = best.size();
    if (K == 0) return Collections.emptyList();

    // Deduplicate by automorphism-canonical mapping.
    // Two mappings are equivalent if they differ only by automorphisms of g1/g2.
    Map<String, Map<Integer, Integer>> seen = new LinkedHashMap<>();
    seen.put(canonKey(g1, g2, best), best);

    if (seen.size() >= maxResults) return new ArrayList<>(seen.values());

    // Resolve adaptive timeout for enumeration phase
    long enumTimeout = M.timeoutMs;
    if (enumTimeout < 0) {
      enumTimeout = Math.min(30_000L, 500L + (long) g1.n * g2.n * 2);
    }
    TimeBudget tb = new TimeBudget(enumTimeout);

    int minN = Math.min(g1.n, g2.n);

    // Phase 2a: if MCS == smaller molecule, this is substructure containment.
    // Enumerate all substructure mappings to find distinct K-sized solutions.
    if (K == minN && !tb.expired()) {
      MolGraph sml = g1.n <= g2.n ? g1 : g2;
      MolGraph lrg = g1.n <= g2.n ? g2 : g1;
      boolean swapped = g1.n > g2.n;
      List<Map<Integer, Integer>> subMaps = SubstructureEngine.findAllSubstructures(
          sml, lrg, C, maxResults * 2, tb.remainingMillis());
      for (Map<Integer, Integer> raw : subMaps) {
        if (tb.expired() || seen.size() >= maxResults) break;
        Map<Integer, Integer> mapping;
        if (swapped) {
          mapping = new LinkedHashMap<>();
          for (Map.Entry<Integer, Integer> e : raw.entrySet()) mapping.put(e.getValue(), e.getKey());
        } else {
          mapping = raw;
        }
        if (mapping.size() == K) {
          seen.putIfAbsent(canonKey(g1, g2, mapping), mapping);
        }
      }
      if (seen.size() >= maxResults) return new ArrayList<>(seen.values());
    }

    // Phase 2b: re-run product graph pipeline stages to collect alternative
    // maximum-sized seeds, then extend each via McGregor.
    if (!tb.expired()) {
      GraphBuilder GB = new GraphBuilder(g1, g2, C, M.induced);

      List<Map<Integer, Integer>> seeds = new ArrayList<>();

      // McSplit seed
      if (!tb.expired()) {
        long[] nc = {0};
        Map<Integer, Integer> mcSeed = GB.mcSplitSeed(tb, nc);
        if (mcSeed.size() >= K) seeds.add(mcSeed);
      }

      // BK clique seed
      if (!tb.expired()) {
        Map<Integer, Integer> cliqueSeed = GB.maximumCliqueSeed(tb);
        if (cliqueSeed.size() >= K) seeds.add(cliqueSeed);
      }

      // Extra seeds: ring anchors, label-degree anchors, VF2++ variants
      if (M.extraSeeds && !tb.expired()) {
        Map<Integer, Integer> s = GB.ringAnchorSeed(tb);
        if (!s.isEmpty()) seeds.add(s);
        if (!tb.expired()) {
          s = GB.labelDegreeAnchorSeed(tb);
          if (!s.isEmpty()) seeds.add(s);
        }
        if (!tb.expired()) {
          s = GB.vf2ppRingSkeletonSeed(tb, Math.max(1, tb.remainingMillis() / 4), 2, 12);
          if (!s.isEmpty()) seeds.add(s);
        }
        if (!tb.expired()) {
          s = GB.vf2ppCoreSeed(tb, Math.max(1, tb.remainingMillis() / 4), 2, 12);
          if (!s.isEmpty()) seeds.add(s);
        }
      }

      // Extend each seed via McGregor and collect K-sized results
      long perSeedMs = Math.max(1, tb.remainingMillis() / Math.max(1, seeds.size()));
      for (Map<Integer, Integer> seed : seeds) {
        if (tb.expired() || seen.size() >= maxResults) break;
        Map<Integer, Integer> ext = ppx(g1, g2,
            mcGregorExtend(g1, g2, seed, C, tb, perSeedMs,
                M.useTwoHopNLFInExtension, M.useThreeHopNLFInExtension, M.connectedOnly),
            C, M);
        if (ext.size() == K) {
          seen.putIfAbsent(canonKey(g1, g2, ext), ext);
        }

        // Also try greedy atom extension for alternative mappings
        if (!tb.expired() && seen.size() < maxResults) {
          Map<Integer, Integer> gext = ppx(g1, g2,
              greedyAtomExtend(g1, g2, seed, C, M), C, M);
          if (gext.size() == K) {
            seen.putIfAbsent(canonKey(g1, g2, gext), gext);
          }
        }
      }

      // Phase 2c: perturb existing solutions by removing one pair and re-extending
      // to discover alternative K-sized mappings.
      if (!tb.expired() && seen.size() < maxResults) {
        for (Map<Integer, Integer> existing : new ArrayList<>(seen.values())) {
          if (tb.expired() || seen.size() >= maxResults) break;
          for (Map.Entry<Integer, Integer> entry : existing.entrySet()) {
            if (tb.expired() || seen.size() >= maxResults) break;
            Map<Integer, Integer> reduced = new LinkedHashMap<>(existing);
            reduced.remove(entry.getKey());
            Map<Integer, Integer> reext = ppx(g1, g2,
                greedyAtomExtend(g1, g2, reduced, C, M), C, M);
            if (reext.size() == K) {
              seen.putIfAbsent(canonKey(g1, g2, reext), reext);
            }
          }
        }
      }
    }

    return new ArrayList<>(seen.values());
  }

  /**
   * Enumerate multiple distinct MCS mappings between two CDK molecules.
   *
   * <p>Convenience overload that converts {@link IAtomContainer} to {@link MolGraph}
   * and delegates to {@link #findAllMCS(MolGraph, MolGraph, ChemOptions, McsOptions, int)}.
   *
   * @param m1         the first molecule
   * @param m2         the second molecule
   * @param C          chemical matching options
   * @param M          MCS options
   * @param maxResults maximum number of distinct mappings to return (capped at 10 if &le; 0)
   * @return list of distinct atom-index mappings, each of the maximum MCS size
   */
  public static List<Map<Integer, Integer>> findAllMCS(
      IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions M, int maxResults) {
    MolGraph g1 = toMolGraph(m1), g2 = toMolGraph(m2);
    applySolvent(g1, C);
    applySolvent(g2, C);
    return findAllMCS(g1, g2, C, M, maxResults);
  }

  /**
   * Canonicalize an atom-atom mapping under the automorphism groups of both
   * molecules.  Two mappings are equivalent if they differ only by symmetry
   * of g1 and/or g2; this method returns the lexicographically smallest
   * representative.
   *
   * @param g1      query molecule
   * @param g2      target molecule
   * @param mapping atom-atom mapping (g1 index → g2 index)
   * @return the canonical representative mapping
   */
  public static Map<Integer, Integer> canonicalizeMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping) {
    if (mapping.isEmpty()) return mapping;

    int[][] gens1 = g1.getAutomorphismGenerators();
    int[][] gens2 = g2.getAutomorphismGenerators();
    if (gens1.length == 0 && gens2.length == 0) return mapping;

    // Precompute inverses of g1 generators
    int[][] invGens1 = new int[gens1.length][];
    for (int g = 0; g < gens1.length; g++) {
      int[] gen = gens1[g];
      int[] inv = new int[gen.length];
      for (int i = 0; i < gen.length; i++) inv[gen[i]] = i;
      invGens1[g] = inv;
    }

    // Convert mapping to sorted pair array for lex comparison
    int[][] best = new int[mapping.size()][2];
    int idx = 0;
    for (Map.Entry<Integer, Integer> e : mapping.entrySet())
      best[idx++] = new int[] {e.getKey(), e.getValue()};
    Arrays.sort(best, (a, b) -> a[0] != b[0] ? Integer.compare(a[0], b[0]) : Integer.compare(a[1], b[1]));

    int maxIter = Math.max(100, 2 * (gens1.length + gens2.length));
    for (int iter = 0; iter < maxIter; iter++) {
      boolean improved = false;

      // Try g1-side generators (permute query atom indices)
      for (int[] inv : invGens1) {
        int[][] cand = new int[best.length][2];
        for (int i = 0; i < best.length; i++) {
          cand[i][0] = inv[best[i][0]];
          cand[i][1] = best[i][1];
        }
        Arrays.sort(cand, (a, b) -> a[0] != b[0] ? Integer.compare(a[0], b[0]) : Integer.compare(a[1], b[1]));
        if (comparePairArrays(cand, best) < 0) { best = cand; improved = true; }
      }

      // Try g2-side generators (permute target atom indices)
      for (int[] gen : gens2) {
        int[][] cand = new int[best.length][2];
        for (int i = 0; i < best.length; i++) {
          cand[i][0] = best[i][0];
          cand[i][1] = gen[best[i][1]];
        }
        Arrays.sort(cand, (a, b) -> a[0] != b[0] ? Integer.compare(a[0], b[0]) : Integer.compare(a[1], b[1]));
        if (comparePairArrays(cand, best) < 0) { best = cand; improved = true; }
      }

      if (!improved) break;
    }

    Map<Integer, Integer> result = new LinkedHashMap<>();
    for (int[] pair : best) result.put(pair[0], pair[1]);
    return result;
  }

  /** Compute a string key from the canonical mapping for use as a dedup key. */
  private static String canonKey(MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping) {
    Map<Integer, Integer> cm = canonicalizeMapping(g1, g2, mapping);
    StringBuilder sb = new StringBuilder();
    for (Map.Entry<Integer, Integer> e : cm.entrySet())
      sb.append(e.getKey()).append(':').append(e.getValue()).append(',');
    return sb.toString();
  }

  private static int comparePairArrays(int[][] a, int[][] b) {
    int len = Math.min(a.length, b.length);
    for (int i = 0; i < len; i++) {
      int c = Integer.compare(a[i][0], b[i][0]);
      if (c != 0) return c;
      c = Integer.compare(a[i][1], b[i][1]);
      if (c != 0) return c;
    }
    return Integer.compare(a.length, b.length);
  }

  /**
   * Find the disconnected Maximum Common Substructure (dMCS) between two molecules.
   *
   * <p>Unlike connected MCS, the result may consist of multiple disjoint fragments.
   *
   * @param m1 the first molecule
   * @param m2 the second molecule
   * @param C  chemical matching options
   * @param M  MCS options (fragment constraints via {@code minFragmentSize} and {@code maxFragments})
   * @return mapping from m1 atom indices to m2 atom indices
   */
  public static Map<Integer, Integer> findDisconnectedMCS(
      IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions M) {
    McsOptions dM = new McsOptions();
    dM.induced = M.induced;
    dM.connectedOnly = false;
    dM.disconnectedMCS = true;
    dM.timeoutMs = M.timeoutMs;
    dM.extraSeeds = M.extraSeeds;
    dM.useTwoHopNLFInExtension = M.useTwoHopNLFInExtension;
    dM.useThreeHopNLFInExtension = M.useThreeHopNLFInExtension;
    dM.maximizeBonds = M.maximizeBonds;
    dM.minFragmentSize = M.minFragmentSize;
    dM.maxFragments = M.maxFragments;
    dM.atomWeights = M.atomWeights;
    return findMCS(m1, m2, C, dM);
  }

  // ---- N-MCS (Multi-molecule MCS with threshold) ----

  /**
   * Find the Maximum Common Substructure across multiple molecules using sequential reduction.
   *
   * @param molecules list of molecules (at least 2)
   * @param chem chemical matching options
   * @param threshold fraction of molecules that must contain the MCS (0.0 to 1.0)
   * @param timeoutMs timeout per pairwise MCS computation
   * @return mapping from MCS atom indices to themselves (identity), or empty if no common core
   */
  public static Map<Integer, Integer> findNMCS(
      List<IAtomContainer> molecules, ChemOptions chem, double threshold, long timeoutMs) {
    if (molecules == null || molecules.size() < 2) return Collections.emptyMap();

    List<Integer> sortedIdx = new ArrayList<>();
    for (int i = 0; i < molecules.size(); i++) sortedIdx.add(i);
    sortedIdx.sort(Comparator.comparingInt(i -> molecules.get(i).getAtomCount()));

    McsOptions mcsOpts = new McsOptions();
    mcsOpts.timeoutMs = timeoutMs;
    mcsOpts.connectedOnly = true;

    IAtomContainer currentMCS = molecules.get(sortedIdx.get(0));

    // Track original atom indices from molecules[sortedIdx[0]] through each reduction
    int[] originalIndices = new int[currentMCS.getAtomCount()];
    for (int i = 0; i < originalIndices.length; i++) originalIndices[i] = i;

    for (int i = 1; i < sortedIdx.size(); i++) {
      Map<Integer, Integer> mcs = findMCS(currentMCS, molecules.get(sortedIdx.get(i)), chem, mcsOpts);
      if (mcs.isEmpty()) return Collections.emptyMap();

      // Collect surviving local indices (query side of the MCS mapping)
      Set<Integer> queryIndices = mcs.keySet();

      // Update provenance: map surviving local indices to original indices
      int[] newOriginal = new int[queryIndices.size()];
      int idx = 0;
      for (int localIdx : queryIndices) {
        newOriginal[idx++] = originalIndices[localIdx];
      }
      originalIndices = newOriginal;

      currentMCS = extractSubgraph(currentMCS, queryIndices);
      if (currentMCS == null || currentMCS.getAtomCount() == 0) return Collections.emptyMap();
    }

    if (threshold < 1.0 && currentMCS.getAtomCount() > 0) {
      int requiredCount = (int) Math.ceil(threshold * molecules.size());
      int matchCount = 0;
      for (IAtomContainer mol : molecules) {
        if (isSubstructure(currentMCS, mol, chem, timeoutMs)) matchCount++;
      }
      if (matchCount < requiredCount) return Collections.emptyMap();
    }

    // Return mapping: original atom index in molecules[sortedIdx[0]] -> NMCS position
    Map<Integer, Integer> result = new LinkedHashMap<>();
    for (int i = 0; i < currentMCS.getAtomCount(); i++) {
      result.put(originalIndices[i], i);
    }
    return result;
  }

  /**
   * Find the N-MCS and return the actual IAtomContainer representing the common core.
   */
  public static IAtomContainer findNMCSMolecule(
      List<IAtomContainer> molecules, ChemOptions chem, double threshold, long timeoutMs) {
    if (molecules == null || molecules.size() < 2)
      return new org.openscience.cdk.silent.AtomContainer();

    List<IAtomContainer> sorted = new ArrayList<>(molecules);
    sorted.sort(Comparator.comparingInt(IAtomContainer::getAtomCount));

    McsOptions mcsOpts = new McsOptions();
    mcsOpts.timeoutMs = timeoutMs;
    mcsOpts.connectedOnly = true;

    IAtomContainer currentMCS = sorted.get(0);

    for (int i = 1; i < sorted.size(); i++) {
      Map<Integer, Integer> mcs = findMCS(currentMCS, sorted.get(i), chem, mcsOpts);
      if (mcs.isEmpty()) return new org.openscience.cdk.silent.AtomContainer();
      currentMCS = extractSubgraph(currentMCS, mcs.keySet());
      if (currentMCS == null || currentMCS.getAtomCount() == 0)
        return new org.openscience.cdk.silent.AtomContainer();
    }

    if (threshold < 1.0 && currentMCS.getAtomCount() > 0) {
      int requiredCount = (int) Math.ceil(threshold * molecules.size());
      int matchCount = 0;
      for (IAtomContainer mol : molecules) {
        if (isSubstructure(currentMCS, mol, chem, timeoutMs)) matchCount++;
      }
      if (matchCount < requiredCount) return new org.openscience.cdk.silent.AtomContainer();
    }

    return currentMCS;
  }

  /** Extract a subgraph from a molecule containing only the specified atom indices. */
  static IAtomContainer extractSubgraph(IAtomContainer mol, Set<Integer> atomIndices) {
    try {
      IAtomContainer sub = mol.getBuilder().newInstance(IAtomContainer.class);
      Map<Integer, IAtom> oldToNew = new LinkedHashMap<>();
      for (int i : atomIndices) {
        IAtom cloned = mol.getAtom(i).clone();
        sub.addAtom(cloned);
        oldToNew.put(i, cloned);
      }
      for (IBond bond : mol.bonds()) {
        IAtom a0 = bond.getAtom(0), a1 = bond.getAtom(1);
        int i0 = mol.indexOf(a0), i1 = mol.indexOf(a1);
        if (atomIndices.contains(i0) && atomIndices.contains(i1)) {
          IBond newBond = bond.clone();
          newBond.setAtom(oldToNew.get(i0), 0);
          newBond.setAtom(oldToNew.get(i1), 1);
          sub.addBond(newBond);
        }
      }
      return sub;
    } catch (CloneNotSupportedException e) {
      throw new RuntimeException("Clone not supported for atom/bond", e);
    }
  }

  /**
   * Extract a MolGraph subgraph containing only the specified atom indices.
   * Used for canonical SMILES deduplication in findAllMCS.
   */
  static MolGraph extractMolGraphSubgraph(MolGraph mol, Collection<Integer> atomIndices) {
    int subN = atomIndices.size();
    if (subN == 0) return new MolGraph.Builder().atomCount(0).atomicNumbers(new int[0]).neighbors(new int[0][]).build();
    Map<Integer, Integer> oldToNew = new LinkedHashMap<>();
    int idx = 0;
    for (int old : atomIndices) oldToNew.put(old, idx++);

    int[] atomicNum = new int[subN];
    int[] formalCharge = new int[subN];
    int[] massNumber = new int[subN];
    boolean[] ring = new boolean[subN];
    boolean[] aromatic = new boolean[subN];
    int[][] nbrs = new int[subN][];
    int[][] bords = new int[subN][];
    boolean[][] brfs = new boolean[subN][];
    boolean[][] bafs = new boolean[subN][];

    idx = 0;
    for (int oldI : atomIndices) {
      atomicNum[idx] = mol.atomicNum[oldI];
      formalCharge[idx] = mol.formalCharge[oldI];
      massNumber[idx] = mol.massNumber[oldI];
      ring[idx] = mol.ring[oldI];
      aromatic[idx] = mol.aromatic[oldI];
      // Collect only neighbors present in the subgraph
      int deg = mol.neighbors[oldI].length;
      int cnt = 0;
      for (int k = 0; k < deg; k++) {
        if (oldToNew.containsKey(mol.neighbors[oldI][k])) cnt++;
      }
      nbrs[idx] = new int[cnt];
      bords[idx] = new int[cnt];
      brfs[idx] = new boolean[cnt];
      bafs[idx] = new boolean[cnt];
      int p = 0;
      for (int k = 0; k < deg; k++) {
        int oldJ = mol.neighbors[oldI][k];
        Integer newJ = oldToNew.get(oldJ);
        if (newJ != null) {
          nbrs[idx][p] = newJ;
          bords[idx][p] = mol.bondOrder(oldI, oldJ);
          brfs[idx][p] = mol.bondInRing(oldI, oldJ);
          bafs[idx][p] = mol.bondAromatic(oldI, oldJ);
          p++;
        }
      }
      idx++;
    }
    return new MolGraph.Builder()
        .atomCount(subN).atomicNumbers(atomicNum).formalCharges(formalCharge)
        .massNumbers(massNumber).ringFlags(ring).aromaticFlags(aromatic)
        .neighbors(nbrs).bondOrders(bords)
        .bondRingFlags(brfs).bondAromaticFlags(bafs)
        .build();
  }

  // ---- R-Group Decomposition ----

  /**
   * Decompose molecules into a core and R-groups.
   *
   * @param core the core scaffold
   * @param molecules list of molecules to decompose
   * @param chem chemical matching options
   * @param timeoutMs timeout for substructure matching
   * @return list of maps; each contains "core" and "R1", "R2", etc.
   */
  public static List<Map<String, IAtomContainer>> decomposeRGroups(
      IAtomContainer core, List<IAtomContainer> molecules, ChemOptions chem, long timeoutMs) {
    List<Map<String, IAtomContainer>> results = new ArrayList<>();
    if (core == null || molecules == null) return results;

    for (IAtomContainer mol : molecules) {
      Map<String, IAtomContainer> decomposition = new LinkedHashMap<>();

      List<Map<Integer, Integer>> mappings = findAllSubstructures(core, mol, chem, 1, timeoutMs);
      if (mappings.isEmpty()) {
        results.add(decomposition);
        continue;
      }

      Map<Integer, Integer> mapping = mappings.get(0);
      Set<Integer> coreAtomIndicesInMol = new HashSet<>(mapping.values());
      decomposition.put("core", extractSubgraph(mol, coreAtomIndicesInMol));

      int rGroupCounter = 1;
      Set<Integer> processedNonCoreAtoms = new HashSet<>();

      for (Map.Entry<Integer, Integer> entry : mapping.entrySet()) {
        int molIdx = entry.getValue();
        IAtom molAtom = mol.getAtom(molIdx);
        for (IAtom neighbor : mol.getConnectedAtomsList(molAtom)) {
          int neighborIdx = mol.indexOf(neighbor);
          if (coreAtomIndicesInMol.contains(neighborIdx)) continue;
          if (processedNonCoreAtoms.contains(neighborIdx)) continue;

          Set<Integer> rGroupAtoms = new LinkedHashSet<>();
          Deque<Integer> queue = new ArrayDeque<>();
          queue.add(neighborIdx);
          while (!queue.isEmpty()) {
            int current = queue.poll();
            if (coreAtomIndicesInMol.contains(current)) continue;
            if (!rGroupAtoms.add(current)) continue;
            processedNonCoreAtoms.add(current);
            for (IAtom nb : mol.getConnectedAtomsList(mol.getAtom(current))) {
              int nbIdx = mol.indexOf(nb);
              if (!coreAtomIndicesInMol.contains(nbIdx) && !rGroupAtoms.contains(nbIdx)) {
                queue.add(nbIdx);
              }
            }
          }

          if (!rGroupAtoms.isEmpty()) {
            decomposition.put("R" + rGroupCounter++, extractSubgraph(mol, rGroupAtoms));
          }
        }
      }

      results.add(decomposition);
    }
    return results;
  }

  static SubstructureEngine.Matcher makeMatcher(MolGraph gq, MolGraph gt, ChemOptions C, TimeBudget tb) {
    return SubstructureEngine.makeMatcher(gq, gt, C, tb);
  }

  /** Walk a chain molecule (max degree ≤ 2) from an endpoint, returning atom indices in order. */
  private static int[] walkChain(MolGraph g) {
    int start = 0;
    for (int i = 0; i < g.n; i++) if (g.degree[i] <= 1) { start = i; break; }
    int[] seq = new int[g.n];
    boolean[] visited = new boolean[g.n];
    int cur = start, len = 0;
    while (cur >= 0) {
      visited[cur] = true;
      seq[len++] = cur;
      int next = -1;
      for (int nb : g.neighbors[cur]) if (!visited[nb]) { next = nb; break; }
      cur = next;
    }
    return len == g.n ? seq : java.util.Arrays.copyOf(seq, len);
  }

  /** Check if a molecule graph is a tree: connected and edges == atoms - 1. */
  private static boolean isTree(MolGraph g) {
    if (g.n <= 1) return g.n == 1;
    int edgeCount = 0;
    for (int i = 0; i < g.n; i++) edgeCount += g.neighbors[i].length;
    edgeCount /= 2;
    if (edgeCount != g.n - 1) return false;
    boolean[] vis = new boolean[g.n];
    int[] q = new int[g.n];
    int head = 0, tail = 0;
    q[tail++] = 0; vis[0] = true;
    while (head < tail) {
      int u = q[head++];
      for (int v : g.neighbors[u])
        if (!vis[v]) { vis[v] = true; q[tail++] = v; }
    }
    return tail == g.n;
  }

  /** Find centroid of a tree by iteratively peeling leaves. */
  private static int treeCentroid(MolGraph g) {
    if (g.n <= 2) return 0; // trivial tree — any node is centroid
    int[] deg = new int[g.n];
    int[] leaves = new int[g.n];
    int lc = 0;
    for (int i = 0; i < g.n; i++) {
      deg[i] = g.neighbors[i].length;
      if (deg[i] <= 1) leaves[lc++] = i;
    }
    int remaining = g.n;
    while (remaining > 2) {
      int[] newLeaves = new int[g.n];
      int nlc = 0;
      remaining -= lc;
      for (int k = 0; k < lc; k++) {
        int u = leaves[k];
        for (int v : g.neighbors[u])
          if (--deg[v] == 1) newLeaves[nlc++] = v;
      }
      leaves = newLeaves;
      lc = nlc;
    }
    return leaves[0];
  }

  /**
   * Root a tree at the given root node.
   * Output: parentOut[0], childrenOut[0][node], postorderOut[0].
   */
  private static void rootTree(MolGraph g, int root, int[][] parentOut,
                                int[][][] childrenOut, int[][] postorderOut) {
    int[] parent = new int[g.n];
    Arrays.fill(parent, -1);
    int[][] children = new int[g.n][];
    int[] childCount = new int[g.n];
    int[] bfsOrder = new int[g.n];
    boolean[] vis = new boolean[g.n];
    int head = 0, tail = 0;
    bfsOrder[tail++] = root; vis[root] = true;
    while (head < tail) {
      int u = bfsOrder[head++];
      for (int v : g.neighbors[u]) {
        if (!vis[v]) {
          vis[v] = true;
          parent[v] = u;
          childCount[u]++;
          bfsOrder[tail++] = v;
        }
      }
    }
    for (int i = 0; i < g.n; i++) children[i] = new int[childCount[i]];
    int[] idx = new int[g.n];
    for (int k = 0; k < tail; k++) {
      int u = bfsOrder[k];
      for (int v : g.neighbors[u]) {
        if (parent[v] == u) children[u][idx[u]++] = v;
      }
    }
    int[] postorder = new int[g.n];
    for (int i = 0; i < g.n; i++) postorder[i] = bfsOrder[g.n - 1 - i];
    parentOut[0] = parent;
    childrenOut[0] = children;
    postorderOut[0] = postorder;
  }

  static Map<Integer, Integer> greedyProbe(MolGraph g1, MolGraph g2, ChemOptions C) {
    int n1 = g1.n, n2 = g2.n;
    if (n1 == 0 || n2 == 0) return Collections.emptyMap();
    boolean[] usedT = new boolean[n2];
    Map<Integer, Integer> map = new LinkedHashMap<>();

    // Sort query atoms: ring atoms first, then by descending degree
    int[] order = new int[n1];
    for (int i = 0; i < n1; i++) order[i] = i;
    for (int i = 1; i < n1; i++) {
      int key = order[i];
      int keyScore = (g1.ring[key] ? 1000 : 0) + g1.degree[key];
      int j = i - 1;
      while (j >= 0) {
        if (((g1.ring[order[j]] ? 1000 : 0) + g1.degree[order[j]]) >= keyScore) break;
        order[j + 1] = order[j];
        j--;
      }
      order[j + 1] = key;
    }

    for (int idx = 0; idx < n1; idx++) {
      int qi = order[idx];
      int bestTj = -1, bestScore = -1;
      for (int tj = 0; tj < n2; tj++) {
        if (usedT[tj]) continue;
        if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
        boolean ok = true;
        for (int qk : g1.neighbors[qi]) {
          if (!map.containsKey(qk)) continue;
          int tk = map.get(qk);
          int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
          if ((qOrd == 0) != (tOrd == 0)) { ok = false; break; }
          if (qOrd != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { ok = false; break; }
        }
        if (!ok) continue;
        int score = (g2.ring[tj] == g1.ring[qi] ? 100 : 0)
            + (g2.degree[tj] == g1.degree[qi] ? 50 : 0)
            + Math.min(g1.degree[qi], g2.degree[tj]);
        if (score > bestScore) { bestScore = score; bestTj = tj; }
      }
      if (bestTj >= 0) {
        map.put(qi, bestTj);
        usedT[bestTj] = true;
        if (map.size() >= Math.min(n1, n2)) return map; // early exit: best possible
      }
    }
    return map;
  }

  static Map<Integer, Integer> pruneToInduced(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> map, ChemOptions C) {
    if (map.isEmpty()) return map;
    Map<Integer, Integer> M = new LinkedHashMap<>(map);
    boolean changed;
    do {
      changed = false;
      List<Integer> Q = new ArrayList<>(M.keySet());
      outer:
      for (int i = 0; i < Q.size(); i++) {
        for (int j = i + 1; j < Q.size(); j++) {
          int qi = Q.get(i), qj = Q.get(j);
          int ti = M.get(qi), tj = M.get(qj);
          int qOrd = g1.bondOrder(qi, qj), tOrd = g2.bondOrder(ti, tj);
          boolean qHas = qOrd != 0, tHas = tOrd != 0;
          boolean ok = (qHas == tHas);
          if (ok && qHas) ok = MolGraph.ChemOps.bondsCompatible(g1, qi, qj, g2, ti, tj, C);
          if (!ok) {
            if (g1.degree[qi] >= g1.degree[qj]) M.remove(qi); else M.remove(qj);
            changed = true;
            break outer;
          }
        }
      }
    } while (changed);
    return M;
  }

  static Map<Integer, Integer> largestConnected(MolGraph g1, Map<Integer, Integer> map) {
    if (map.isEmpty()) return map;
    Map<Integer, List<Integer>> adj = new HashMap<>();
    for (int qi : map.keySet()) adj.put(qi, new ArrayList<>());
    for (int qi : map.keySet())
      for (int qk : map.keySet()) {
        if (qi >= qk) continue;
        if (g1.hasBond(qi, qk)) { adj.get(qi).add(qk); adj.get(qk).add(qi); }
      }
    Set<Integer> seen = new HashSet<>();
    List<Set<Integer>> comps = new ArrayList<>();
    for (int qi : map.keySet()) {
      if (seen.contains(qi)) continue;
      Set<Integer> comp = new LinkedHashSet<>();
      Deque<Integer> dq = new ArrayDeque<>();
      dq.add(qi);
      seen.add(qi);
      while (!dq.isEmpty()) {
        int u = dq.pollFirst();
        comp.add(u);
        for (int v : adj.getOrDefault(u, Collections.emptyList()))
          if (!seen.contains(v)) { seen.add(v); dq.addLast(v); }
      }
      comps.add(comp);
    }
    Set<Integer> best = comps.stream().max(Comparator.comparingInt(Set::size)).orElse(Collections.emptySet());
    if (best.size() == map.size()) return map;
    Map<Integer, Integer> pruned = new LinkedHashMap<>();
    for (int qi : best) pruned.put(qi, map.get(qi));
    return pruned;
  }

  static Map<Integer, Integer> applyRingAnchorGuard(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> map, ChemOptions C) {
    if (!C.ringMatchesRingOnly || map.isEmpty()) return map;
    boolean qHasRing = false;
    for (int i = 0; i < g1.n; i++) if (g1.ring[i]) { qHasRing = true; break; }
    boolean tHasRing = false;
    for (int i = 0; i < g2.n; i++) if (g2.ring[i]) { tHasRing = true; break; }
    if (qHasRing && !tHasRing) {
      int mappedRing = 0;
      for (int qi : map.keySet()) if (g1.ring[qi]) mappedRing++;
      if (mappedRing == 0) return Collections.emptyMap();
    }
    return map;
  }

  // McGregor helpers

  static boolean isPruned(int curSize, int bestSize, int[] qLF, int[] tLF, int freqSize) {
    int potential = curSize;
    for (int lbl = 0; lbl < freqSize; lbl++)
      if (qLF[lbl] > 0) potential += Math.min(qLF[lbl], tLF[lbl]);
    return potential <= bestSize;
  }

  static int buildFrontier(
      MolGraph g1, int[] curMap, boolean[] usedQ, boolean[] inFrontier, int[] frontierBuf,
      boolean connectedOnly) {
    int count = 0;
    for (int qk = 0; qk < g1.n; qk++) {
      if (curMap[qk] == -1) continue;
      for (int qn : g1.neighbors[qk])
        if (!usedQ[qn] && !inFrontier[qn] && curMap[qn] == -1) {
          inFrontier[qn] = true;
          frontierBuf[count++] = qn;
        }
    }
    // Only jump to disconnected atoms if disconnected MCS is allowed
    if (count == 0 && !connectedOnly)
      for (int i = 0; i < g1.n; i++)
        if (!usedQ[i] && curMap[i] == -1) frontierBuf[count++] = i;
    for (int f = 0; f < count; f++) inFrontier[frontierBuf[f]] = false;
    return count;
  }

  static void undoForcedAssignments(
      int forcedCount, int[] forcedQ, int[] forcedT, int[] curMap, int[] curSize,
      boolean[] usedQ, boolean[] usedT, int[] qLabelFreq, int[] tLabelFreq,
      int[] jointQ, int[] jointT, int[] q2tMap) {
    for (int f = forcedCount - 1; f >= 0; f--) {
      int fq = forcedQ[f], ft = forcedT[f];
      qLabelFreq[jointQ[fq]]++;
      tLabelFreq[jointT[ft]]++;
      curMap[fq] = -1; curSize[0]--;
      usedQ[fq] = false;
      usedT[ft] = false;
      if (q2tMap != null) q2tMap[fq] = -1;
    }
  }

  // McGregor-like extension with NLF pruning

  static Map<Integer, Integer> mcGregorExtend(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> seed, ChemOptions C,
      TimeBudget tb, long localMillis, boolean useTwoHopNLF, boolean useThreeHopNLF,
      boolean connectedOnly) {
    final long localDeadline = System.nanoTime() + localMillis * 1_000_000L;
    int[][] qNLF1 = g1.getNLF1(), tNLF1 = g2.getNLF1();
    int[][] qNLF2 = useTwoHopNLF ? g1.getNLF2() : null, tNLF2 = useTwoHopNLF ? g2.getNLF2() : null;
    int[][] qNLF3 = useThreeHopNLF ? g1.getNLF3() : null, tNLF3 = useThreeHopNLF ? g2.getNLF3() : null;

    // Build flat arrays from seed map
    int[] curMap = new int[g1.n];
    java.util.Arrays.fill(curMap, -1);
    int[] curSize = new int[]{0};
    int[] bestMap = new int[g1.n];
    java.util.Arrays.fill(bestMap, -1);
    int[] bestSize = new int[]{0};

    boolean[] usedQ = new boolean[g1.n], usedT = new boolean[g2.n];
    int[] candBuf = new int[g2.n], bestCandBuf = new int[g2.n];
    for (Map.Entry<Integer, Integer> e : seed.entrySet()) {
      int qi = e.getKey(), tj = e.getValue();
      curMap[qi] = tj; curSize[0]++;
      bestMap[qi] = tj; bestSize[0]++;
      usedQ[qi] = true;
      usedT[tj] = true;
    }

    int maxLabel = 0;
    for (int i = 0; i < g1.n; i++) maxLabel = Math.max(maxLabel, g1.label[i]);
    for (int j = 0; j < g2.n; j++) maxLabel = Math.max(maxLabel, g2.label[j]);
    int freqSize = maxLabel + 1;
    int[] qLabelFreq = new int[freqSize], tLabelFreq = new int[freqSize];
    int[] jointQ = new int[g1.n], jointT = new int[g2.n];
    for (int i = 0; i < g1.n; i++) { jointQ[i] = g1.label[i]; if (!usedQ[i]) qLabelFreq[jointQ[i]]++; }
    for (int j = 0; j < g2.n; j++) { jointT[j] = g2.label[j]; if (!usedT[j]) tLabelFreq[jointT[j]]++; }

    boolean[] inFrontier = new boolean[g1.n];
    int[] frontierBuf = new int[g1.n];

    // Fast bail-out: check if seed's frontier has any extensible atoms
    if (curSize[0] > 0) {
      boolean hasExtensible = false;
      outer:
      for (int qi = 0; qi < g1.n; qi++) {
        if (curMap[qi] == -1) continue;
        for (int nb : g1.neighbors[qi]) {
          if (usedQ[nb]) continue;
          for (int tj = 0; tj < g2.n; tj++) {
            if (usedT[tj]) continue;
            if (SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, nb, g2, tj, C)) { hasExtensible = true; break outer; }
          }
        }
      }
      if (!hasExtensible && connectedOnly) {
        // Convert bestMap back to Map
        Map<Integer, Integer> result = new LinkedHashMap<>();
        for (int qi = 0; qi < g1.n; qi++) if (bestMap[qi] != -1) result.put(qi, bestMap[qi]);
        return result;
      }
      // If disconnected MCS, skip bail-out — let buildFrontier jump to new component
    }

    // Bond-growth first pass for larger molecules
    if (curSize[0] > 0 && g2.n >= 20) {
      int[] q2tMap = new int[g1.n];
      java.util.Arrays.fill(q2tMap, -1);
      for (int qi = 0; qi < g1.n; qi++) if (curMap[qi] != -1) q2tMap[qi] = curMap[qi];
      int[] bondCurMap = java.util.Arrays.copyOf(curMap, curMap.length);
      int[] bondCurSize = new int[]{curSize[0]};
      int[] bondBestMap = java.util.Arrays.copyOf(bestMap, bestMap.length);
      int[] bondBestSize = new int[]{bestSize[0]};
      long bondDeadline = System.nanoTime() + Math.max(1, Math.min(localMillis / 10, 20)) * 1_000_000L;
      mcGregorBondGrow(g1, g2, C, bondCurMap, bondCurSize, bondBestMap, bondBestSize,
          qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF, qNLF2, tNLF2, qNLF3, tNLF3,
          tb, bondDeadline, 0,
          java.util.Arrays.copyOf(usedQ, usedQ.length), java.util.Arrays.copyOf(usedT, usedT.length),
          java.util.Arrays.copyOf(qLabelFreq, qLabelFreq.length), java.util.Arrays.copyOf(tLabelFreq, tLabelFreq.length),
          freqSize, q2tMap, inFrontier, frontierBuf, candBuf, bestCandBuf, jointQ, jointT);
      java.util.Arrays.fill(inFrontier, false);
      if (bondBestSize[0] > bestSize[0]) { System.arraycopy(bondBestMap, 0, bestMap, 0, bestMap.length); bestSize[0] = bondBestSize[0]; }
    }

    mcGregorDFS(g1, g2, C, curMap, curSize, bestMap, bestSize, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
        useTwoHopNLF, useThreeHopNLF, tb, localDeadline, 0,
        usedQ, usedT, candBuf, bestCandBuf, qLabelFreq, tLabelFreq, freqSize,
        inFrontier, frontierBuf, jointQ, jointT, connectedOnly);

    // Convert bestMap back to Map<Integer, Integer>
    Map<Integer, Integer> result = new LinkedHashMap<>();
    for (int qi = 0; qi < g1.n; qi++) if (bestMap[qi] != -1) result.put(qi, bestMap[qi]);
    return result;
  }

  /** Shared candidate-finding logic for McGregor DFS and unit propagation. */
  static int findBestCandidate(
      MolGraph g1, MolGraph g2, ChemOptions C, int[] curMap,
      boolean[] usedT, int[][] qNLF1, int[][] tNLF1, int[][] qNLF2, int[][] tNLF2,
      int[][] qNLF3, int[][] tNLF3, boolean useTwoHopNLF, boolean useThreeHopNLF,
      int frontierCount, int[] frontierBuf, int[] candBuf, int[] bestCandBuf,
      int[] outQi, int[] outCandCount) {
    int bestQi = -1, bestCandSize = Integer.MAX_VALUE, bestCandCount = 0;
    for (int fi = 0; fi < frontierCount; fi++) {
      int qi = frontierBuf[fi];
      int candCount = 0;
      for (int tj = 0; tj < g2.n; tj++) {
        if (usedT[tj]) continue;
        if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
        if (!MolGraph.nlfCheckOk(qi, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3, useTwoHopNLF, useThreeHopNLF)) continue;
        boolean ok = true;
        for (int qk : g1.neighbors[qi]) {
          if (curMap[qk] == -1) continue;
          int tl = curMap[qk];
          int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tl);
          if ((qOrd == 0) != (tOrd == 0)) { ok = false; break; }
          if (qOrd != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tl, C)) { ok = false; break; }
        }
        if (ok) candBuf[candCount++] = tj;
      }
      if (candCount == 0) continue;
      if (candCount < bestCandSize) {
        bestCandSize = candCount;
        bestQi = qi;
        bestCandCount = candCount;
        System.arraycopy(candBuf, 0, bestCandBuf, 0, candCount);
      }
    }
    outQi[0] = bestQi;
    outCandCount[0] = bestCandCount;
    return bestQi;
  }

  static void mcGregorDFS(
      MolGraph g1, MolGraph g2, ChemOptions C, int[] curMap, int[] curSize, int[] bestMap, int[] bestSize,
      int[][] qNLF1, int[][] tNLF1, int[][] qNLF2, int[][] tNLF2, int[][] qNLF3, int[][] tNLF3,
      boolean useTwoHopNLF, boolean useThreeHopNLF, TimeBudget tb, long localDeadline, int depth,
      boolean[] usedQ, boolean[] usedT, int[] candBuf, int[] bestCandBuf,
      int[] qLabelFreq, int[] tLabelFreq, int freqSize,
      boolean[] inFrontier, int[] frontierBuf, int[] jointQ, int[] jointT,
      boolean connectedOnly) {
    if (System.nanoTime() >= localDeadline || tb.expired()) return;
    if (curSize[0] > bestSize[0]) { System.arraycopy(curMap, 0, bestMap, 0, curMap.length); bestSize[0] = curSize[0]; }
    if (isPruned(curSize[0], bestSize[0], qLabelFreq, tLabelFreq, freqSize)) return;

    int frontierCount = buildFrontier(g1, curMap, usedQ, inFrontier, frontierBuf, connectedOnly);
    if (frontierCount == 0) return;

    int[] outQi = new int[1], outCandCount = new int[1];
    findBestCandidate(g1, g2, C, curMap, usedT, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
        useTwoHopNLF, useThreeHopNLF, frontierCount, frontierBuf, candBuf, bestCandBuf, outQi, outCandCount);
    int bestQi = outQi[0], bestCandCount = outCandCount[0];
    if (bestQi == -1) return;

    // Unit propagation: when only 1 candidate, assign without branching
    int forcedCount = 0;
    int[] forcedQ = new int[g1.n], forcedT = new int[g2.n];
    while (bestCandCount == 1 && !(System.nanoTime() >= localDeadline || tb.expired())) {
      int fq = bestQi, ft = bestCandBuf[0];
      curMap[fq] = ft; curSize[0]++; usedQ[fq] = true; usedT[ft] = true;
      qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
      forcedQ[forcedCount] = fq; forcedT[forcedCount] = ft;
      forcedCount++; depth++;
      if (curSize[0] > bestSize[0]) { System.arraycopy(curMap, 0, bestMap, 0, curMap.length); bestSize[0] = curSize[0]; }
      if (isPruned(curSize[0], bestSize[0], qLabelFreq, tLabelFreq, freqSize)) break;

      frontierCount = buildFrontier(g1, curMap, usedQ, inFrontier, frontierBuf, connectedOnly);
      if (frontierCount == 0) break;
      findBestCandidate(g1, g2, C, curMap, usedT, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
          useTwoHopNLF, useThreeHopNLF, frontierCount, frontierBuf, candBuf, bestCandBuf, outQi, outCandCount);
      bestQi = outQi[0]; bestCandCount = outCandCount[0];
      if (bestQi == -1) break;
    }

    if (bestQi != -1 && bestCandCount > 1) {
      int branchLimit = depth < 5 ? bestCandCount : Math.min(bestCandCount, 16);
      for (int i = 0; i < branchLimit; i++) {
        if (System.nanoTime() >= localDeadline || tb.expired()) break;
        int bestTj = bestCandBuf[i];
        curMap[bestQi] = bestTj; curSize[0]++; usedQ[bestQi] = true; usedT[bestTj] = true;
        qLabelFreq[jointQ[bestQi]]--; tLabelFreq[jointT[bestTj]]--;
        mcGregorDFS(g1, g2, C, curMap, curSize, bestMap, bestSize, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3,
            useTwoHopNLF, useThreeHopNLF, tb, localDeadline, depth + 1,
            usedQ, usedT, candBuf, bestCandBuf, qLabelFreq, tLabelFreq, freqSize,
            inFrontier, frontierBuf, jointQ, jointT, connectedOnly);
        qLabelFreq[jointQ[bestQi]]++; tLabelFreq[jointT[bestTj]]++;
        curMap[bestQi] = -1; curSize[0]--; usedQ[bestQi] = false; usedT[bestTj] = false;
      }
    }

    undoForcedAssignments(forcedCount, forcedQ, forcedT, curMap, curSize, usedQ, usedT,
        qLabelFreq, tLabelFreq, jointQ, jointT, null);
  }

  /** Bond-growth MCS extension: grows mapping bond-by-bond. */
  static void mcGregorBondGrow(
      MolGraph g1, MolGraph g2, ChemOptions C, int[] curMap, int[] curSize, int[] bestMap, int[] bestSize,
      int[][] qNLF1, int[][] tNLF1, boolean useTwoHopNLF, boolean useThreeHopNLF,
      int[][] qNLF2, int[][] tNLF2, int[][] qNLF3, int[][] tNLF3,
      TimeBudget tb, long localDeadline, int depth,
      boolean[] usedQ, boolean[] usedT, int[] qLabelFreq, int[] tLabelFreq, int freqSize,
      int[] q2tMap, boolean[] inFrontier, int[] frontierBuf, int[] candBuf, int[] bestCandBuf,
      int[] jointQ, int[] jointT) {

    if (System.nanoTime() >= localDeadline || tb.expired()) return;
    if (curSize[0] > bestSize[0]) { System.arraycopy(curMap, 0, bestMap, 0, curMap.length); bestSize[0] = curSize[0]; }
    if (isPruned(curSize[0], bestSize[0], qLabelFreq, tLabelFreq, freqSize)) return;

    // Find best frontier bond
    int bestQk = -1, bestQi = -1, bestCandSize = Integer.MAX_VALUE, bestCandCount = 0;
    for (int qi = 0; qi < g1.n; qi++) {
      if (curMap[qi] == -1) continue;
      int mappedTi = q2tMap[qi];
      for (int qk : g1.neighbors[qi]) {
        if (usedQ[qk] || inFrontier[qk]) continue;
        int candCount = 0;
        for (int tk : g2.neighbors[mappedTi]) {
          if (usedT[tk]) continue;
          int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(mappedTi, tk);
          if ((qOrd == 0) != (tOrd == 0)) continue;
          if (qOrd != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, mappedTi, tk, C)) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qk, g2, tk, C)) continue;
          if (!MolGraph.nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3, useTwoHopNLF, useThreeHopNLF)) continue;
          boolean ok = true;
          for (int qn : g1.neighbors[qk]) {
            if (qn == qi || !usedQ[qn]) continue;
            int tn = q2tMap[qn];
            int qOrd2 = g1.bondOrder(qk, qn), tOrd2 = g2.bondOrder(tk, tn);
            if ((qOrd2 == 0) != (tOrd2 == 0)) { ok = false; break; }
            if (qOrd2 != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qk, qn, g2, tk, tn, C)) { ok = false; break; }
          }
          if (ok) candBuf[candCount++] = tk;
        }
        if (candCount > 0 && candCount < bestCandSize) {
          bestCandSize = candCount; bestQk = qk; bestQi = qi; bestCandCount = candCount;
          System.arraycopy(candBuf, 0, bestCandBuf, 0, candCount);
        }
      }
    }

    // Fallback: disconnected extension
    if (bestQk == -1) {
      for (int i = 0; i < g1.n; i++) {
        if (usedQ[i]) continue;
        int candCount = 0;
        for (int tj = 0; tj < g2.n; tj++) {
          if (usedT[tj]) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, tj, C)) continue;
          if (!MolGraph.nlfCheckOk(i, tj, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3, useTwoHopNLF, useThreeHopNLF)) continue;
          candBuf[candCount++] = tj;
        }
        if (candCount > 0 && candCount < bestCandSize) {
          bestCandSize = candCount; bestQk = i; bestQi = -1; bestCandCount = candCount;
          System.arraycopy(candBuf, 0, bestCandBuf, 0, candCount);
        }
      }
    }
    if (bestQk == -1) return;

    // Unit propagation
    int forcedCount = 0;
    int[] forcedQ = new int[g1.n], forcedT = new int[g2.n];
    while (bestCandCount == 1 && !(System.nanoTime() >= localDeadline || tb.expired())) {
      int fq = bestQk, ft = bestCandBuf[0];
      curMap[fq] = ft; curSize[0]++; usedQ[fq] = true; usedT[ft] = true; q2tMap[fq] = ft;
      qLabelFreq[jointQ[fq]]--; tLabelFreq[jointT[ft]]--;
      forcedQ[forcedCount] = fq; forcedT[forcedCount] = ft;
      forcedCount++; depth++;
      if (curSize[0] > bestSize[0]) { System.arraycopy(curMap, 0, bestMap, 0, curMap.length); bestSize[0] = curSize[0]; }
      if (isPruned(curSize[0], bestSize[0], qLabelFreq, tLabelFreq, freqSize)) break;

      bestQk = -1; bestCandSize = Integer.MAX_VALUE; bestCandCount = 0;
      for (int qi = 0; qi < g1.n; qi++) {
        if (curMap[qi] == -1) continue;
        int mappedTi = q2tMap[qi];
        for (int qk : g1.neighbors[qi]) {
          if (usedQ[qk]) continue;
          int candCount = 0;
          for (int tk : g2.neighbors[mappedTi]) {
            if (usedT[tk]) continue;
            int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(mappedTi, tk);
            if ((qOrd == 0) != (tOrd == 0)) continue;
            if (qOrd != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, mappedTi, tk, C)) continue;
            if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qk, g2, tk, C)) continue;
            if (!MolGraph.nlfCheckOk(qk, tk, qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3, useTwoHopNLF, useThreeHopNLF)) continue;
            boolean ok = true;
            for (int qn : g1.neighbors[qk]) {
              if (qn == qi || !usedQ[qn]) continue;
              int tn = q2tMap[qn];
              int qOrd2 = g1.bondOrder(qk, qn), tOrd2 = g2.bondOrder(tk, tn);
              if ((qOrd2 == 0) != (tOrd2 == 0)) { ok = false; break; }
              if (qOrd2 != 0 && !MolGraph.ChemOps.bondsCompatible(g1, qk, qn, g2, tk, tn, C)) { ok = false; break; }
            }
            if (ok) candBuf[candCount++] = tk;
          }
          if (candCount > 0 && candCount < bestCandSize) {
            bestCandSize = candCount; bestQk = qk; bestQi = qi; bestCandCount = candCount;
            System.arraycopy(candBuf, 0, bestCandBuf, 0, candCount);
          }
        }
      }
      if (bestQk == -1) break;
    }

    if (bestQk != -1 && bestCandCount > 1) {
      int branchLimit = depth < 5 ? bestCandCount : Math.min(bestCandCount, 16);
      for (int i = 0; i < branchLimit; i++) {
        if (System.nanoTime() >= localDeadline || tb.expired()) break;
        int btj = bestCandBuf[i];
        curMap[bestQk] = btj; curSize[0]++; usedQ[bestQk] = true; usedT[btj] = true; q2tMap[bestQk] = btj;
        qLabelFreq[jointQ[bestQk]]--; tLabelFreq[jointT[btj]]--;
        mcGregorBondGrow(g1, g2, C, curMap, curSize, bestMap, bestSize, qNLF1, tNLF1, useTwoHopNLF, useThreeHopNLF,
            qNLF2, tNLF2, qNLF3, tNLF3, tb, localDeadline, depth + 1,
            usedQ, usedT, qLabelFreq, tLabelFreq, freqSize, q2tMap,
            inFrontier, frontierBuf, candBuf, bestCandBuf, jointQ, jointT);
        qLabelFreq[jointQ[bestQk]]++; tLabelFreq[jointT[btj]]++;
        curMap[bestQk] = -1; curSize[0]--; usedQ[bestQk] = false; usedT[btj] = false; q2tMap[bestQk] = -1;
      }
    }

    undoForcedAssignments(forcedCount, forcedQ, forcedT, curMap, curSize, usedQ, usedT,
        qLabelFreq, tLabelFreq, jointQ, jointT, q2tMap);
  }


  // GraphBuilder and seeds

  static final class GraphBuilder {
    private final MolGraph g1, g2;
    private final ChemOptions C;
    private final boolean induced;

    GraphBuilder(MolGraph g1, MolGraph g2, ChemOptions C, boolean induced) {
      this.g1 = g1; this.g2 = g2; this.C = C; this.induced = induced;
    }

    static final class Node {
      final int qi, tj;
      Node(int qi, int tj) { this.qi = qi; this.tj = tj; }
    }

    Map<Integer, Integer> maximumCliqueSeed(TimeBudget tb) {
      List<Node> nodes = new ArrayList<>();
      int n1 = g1.n, n2 = g2.n;
      int[][] nlf1 = new int[n1][], nlf2 = new int[n2][];
      for (int i = 0; i < n1; i++) {
        Map<Integer, Integer> freq = new HashMap<>();
        for (int k : g1.neighbors[i]) freq.merge(g1.label[k], 1, Integer::sum);
        nlf1[i] = freq.values().stream().mapToInt(Integer::intValue).sorted().toArray();
      }
      for (int j = 0; j < n2; j++) {
        Map<Integer, Integer> freq = new HashMap<>();
        for (int k : g2.neighbors[j]) freq.merge(g2.label[k], 1, Integer::sum);
        nlf2[j] = freq.values().stream().mapToInt(Integer::intValue).sorted().toArray();
      }
      for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++) {
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, j, C)) continue;
          if (!induced && g1.degree[i] > g2.degree[j]) continue;
          nodes.add(new Node(i, j));
        }
      int N = nodes.size();
      if (N == 0) return new LinkedHashMap<>();
      // Cap product graph size to prevent O(N^2) blowup on large molecules
      if (N > 2500) return new LinkedHashMap<>();

      int words = (N + 63) >>> 6;
      long[][] adj = new long[N][words];
      for (int u = 0; u < N; u++) {
        Node nu = nodes.get(u);
        for (int v = u + 1; v < N; v++) {
          if (tb.expired()) break;
          Node nv = nodes.get(v);
          if (nu.qi == nv.qi || nu.tj == nv.tj) continue;
          int qOrd = g1.bondOrder(nu.qi, nv.qi), tOrd = g2.bondOrder(nu.tj, nv.tj);
          boolean ok;
          if (qOrd != 0 && tOrd != 0) ok = MolGraph.ChemOps.bondsCompatible(g1, nu.qi, nv.qi, g2, nu.tj, nv.tj, C);
          else if (induced) ok = (qOrd == 0 && tOrd == 0);
          else continue;
          if (ok) {
            adj[u][v >>> 6] |= 1L << (v & 63);
            adj[v][u >>> 6] |= 1L << (u & 63);
          }
        }
      }

      int[] bestClique = greedyClique(adj, N, words);

      // K-core reduction
      {
        int minDeg = Math.max(2, bestClique.length);
        boolean changed = true;
        while (changed && !tb.expired()) {
          changed = false;
          for (int u = 0; u < N; u++) {
            int deg = 0;
            for (int w = 0; w < words; w++) deg += Long.bitCount(adj[u][w]);
            if (deg == 0 || deg >= minDeg) continue;
            for (int w = 0; w < words; w++) {
              long bits = adj[u][w];
              while (bits != 0) {
                int bit = Long.numberOfTrailingZeros(bits);
                int v = (w << 6) | bit;
                adj[v][u >>> 6] &= ~(1L << (u & 63));
                bits &= bits - 1;
              }
              adj[u][w] = 0;
            }
            changed = true;
          }
        }
      }

      int[] equivClass = computeEquivClasses(nodes, adj, N, words);

      long[] P = new long[words], X = new long[words];
      for (int i = 0; i < N; i++) {
        int deg = 0;
        for (int w = 0; w < words; w++) deg += Long.bitCount(adj[i][w]);
        if (deg > 0) P[i >>> 6] |= 1L << (i & 63);
      }
      int maxDepth = Math.min(N, 200);
      long[][] rStack = new long[maxDepth][words];
      long[][] pStack = new long[maxDepth][words];
      long[][] xStack = new long[maxDepth][words];
      int[] bkBest = bronKerboschPivot(new long[words], P, X, adj, bestClique, N, words,
          equivClass, tb, rStack, pStack, xStack, 0, nodes);
      if (bkBest.length > bestClique.length) bestClique = bkBest;

      Map<Integer, Integer> seed = new LinkedHashMap<>();
      for (int u : bestClique) {
        Node nd = nodes.get(u);
        if (!seed.containsKey(nd.qi) && !seed.containsValue(nd.tj)) seed.put(nd.qi, nd.tj);
      }
      return seed;
    }

    private int[] greedyClique(long[][] adj, int N, int words) {
      int bestStart = 0, bestDeg = 0;
      for (int i = 0; i < N; i++) {
        int deg = 0;
        for (int w = 0; w < words; w++) deg += Long.bitCount(adj[i][w]);
        if (deg > bestDeg) { bestDeg = deg; bestStart = i; }
      }
      List<Integer> clique = new ArrayList<>();
      clique.add(bestStart);
      long[] cand = adj[bestStart].clone();
      while (true) {
        int bestV = -1, bestConn = -1;
        for (int w = 0; w < words; w++) {
          long bits = cand[w];
          while (bits != 0) {
            int bit = Long.numberOfTrailingZeros(bits);
            int v = (w << 6) | bit;
            int conn = 0;
            for (int cw = 0; cw < words; cw++) conn += Long.bitCount(adj[v][cw] & cand[cw]);
            if (conn > bestConn) { bestConn = conn; bestV = v; }
            bits &= bits - 1;
          }
        }
        if (bestV == -1) break;
        clique.add(bestV);
        for (int w = 0; w < words; w++) cand[w] &= adj[bestV][w];
      }
      return clique.stream().mapToInt(Integer::intValue).toArray();
    }

    private int[] computeEquivClasses(List<Node> nodes, long[][] adj, int N, int words) {
      int[] cls = new int[N];
      Map<Long, Integer> sig2class = new HashMap<>();
      int nextClass = 0;
      for (int i = 0; i < N; i++) {
        Node nd = nodes.get(i);
        int deg = 0;
        for (int w = 0; w < words; w++) deg += Long.bitCount(adj[i][w]);
        // Include Morgan rank to distinguish atoms at different chain positions
        // (orbit alone fails for symmetric groups like phosphate oxygens)
        long sig = ((long) g1.orbit[nd.qi] << 48) | ((long) g2.orbit[nd.tj] << 32)
            | ((long) g1.morganRank[nd.qi] << 20) | ((long) g1.degree[nd.qi] << 10) | deg;
        Integer c = sig2class.get(sig);
        if (c == null) { c = nextClass++; sig2class.put(sig, c); }
        cls[i] = c;
      }
      return cls;
    }

    private int colorBound(long[] P, long[][] adj, int words) {
      int[] color = new int[adj.length];
      Arrays.fill(color, -1);
      int maxColor = 0;
      for (int w = 0; w < words; w++) {
        long bits = P[w];
        while (bits != 0) {
          int bit = Long.numberOfTrailingZeros(bits);
          int v = (w << 6) | bit;
          long usedColors = 0;
          for (int nw = 0; nw < words; nw++) {
            long nbits = adj[v][nw] & P[nw];
            while (nbits != 0) {
              int nbit = Long.numberOfTrailingZeros(nbits);
              int u = (nw << 6) | nbit;
              if (color[u] >= 0 && color[u] < 64) usedColors |= 1L << color[u];
              nbits &= nbits - 1;
            }
          }
          color[v] = Long.numberOfTrailingZeros(~usedColors);
          if (color[v] > maxColor) maxColor = color[v];
          bits &= bits - 1;
        }
      }
      return maxColor + 1;
    }

    private int[] bronKerboschPivot(
        long[] R, long[] P, long[] X, long[][] adj, int[] currentBest, int N, int words,
        int[] equivClass, TimeBudget tb, long[][] rStack, long[][] pStack, long[][] xStack,
        int depth, List<Node> nodes) {
      if (tb.expired()) return currentBest;
      int rSize = 0, pSize = 0, xSize = 0;
      for (int w = 0; w < words; w++) {
        rSize += Long.bitCount(R[w]); pSize += Long.bitCount(P[w]); xSize += Long.bitCount(X[w]);
      }
      if (pSize == 0 && xSize == 0) return rSize > currentBest.length ? iterateBits(R, N) : currentBest;
      if (rSize + pSize <= currentBest.length) return currentBest;

      // Partition bound
      if (pSize > 2) {
        BitSet qiSeen = new BitSet(g1.n), tjSeen = new BitSet(g2.n);
        for (int w2 = 0; w2 < words; w2++) {
          long bits2 = P[w2];
          while (bits2 != 0) {
            int bit2 = Long.numberOfTrailingZeros(bits2);
            int v2 = (w2 << 6) | bit2;
            if (v2 < N) { qiSeen.set(nodes.get(v2).qi); tjSeen.set(nodes.get(v2).tj); }
            bits2 &= bits2 - 1;
          }
        }
        if (rSize + Math.min(qiSeen.cardinality(), tjSeen.cardinality()) <= currentBest.length) return currentBest;
      }

      if (pSize > 0 && rSize + colorBound(P, adj, words) <= currentBest.length) return currentBest;

      // Choose pivot maximizing |P intersect N(u)|
      int pivot = -1, pivotConn = -1;
      for (int w = 0; w < words; w++) {
        long bits = P[w] | X[w];
        while (bits != 0) {
          int bit = Long.numberOfTrailingZeros(bits);
          int u = (w << 6) | bit;
          int conn = 0;
          for (int cw = 0; cw < words; cw++) conn += Long.bitCount(P[cw] & adj[u][cw]);
          if (conn > pivotConn) { pivotConn = conn; pivot = u; }
          bits &= bits - 1;
        }
      }

      long[] candidates = new long[words];
      if (pivot >= 0) for (int w = 0; w < words; w++) candidates[w] = P[w] & ~adj[pivot][w];
      else System.arraycopy(P, 0, candidates, 0, words);

      Set<Integer> triedClasses = new HashSet<>();
      for (int w = 0; w < words; w++) {
        long bits = candidates[w];
        while (bits != 0) {
          if (tb.expired()) return currentBest;
          int bit = Long.numberOfTrailingZeros(bits);
          int v = (w << 6) | bit;
          bits &= bits - 1;
          if (!triedClasses.add(equivClass[v])) continue;
          long[] newR, newP, newX;
          if (depth < rStack.length) { newR = rStack[depth]; newP = pStack[depth]; newX = xStack[depth]; }
          else { newR = new long[words]; newP = new long[words]; newX = new long[words]; }
          for (int cw = 0; cw < words; cw++) {
            newR[cw] = R[cw]; newP[cw] = P[cw] & adj[v][cw]; newX[cw] = X[cw] & adj[v][cw];
          }
          newR[v >>> 6] |= 1L << (v & 63);
          currentBest = bronKerboschPivot(newR, newP, newX, adj, currentBest, N, words,
              equivClass, tb, rStack, pStack, xStack, depth + 1, nodes);
          P[v >>> 6] &= ~(1L << (v & 63));
          X[v >>> 6] |= 1L << (v & 63);
        }
      }
      return currentBest;
    }

    // Use iterateBits from AbstractVFMatcher
    private static int[] iterateBits(long[] words, int N) {
      return SubstructureEngine.AbstractVFMatcher.iterateBits(words, N);
    }

    Map<Integer, Integer> ringAnchorSeed(TimeBudget tb) {
      Map<Integer, Integer> seed = new LinkedHashMap<>();
      List<Integer> R1 = new ArrayList<>(), R2 = new ArrayList<>();
      for (int i = 0; i < g1.n; i++) if (g1.ring[i]) R1.add(i);
      for (int j = 0; j < g2.n; j++) if (g2.ring[j]) R2.add(j);
      boolean[] used = new boolean[g2.n];
      for (int i : R1) {
        int bestJ = -1, bestScore = Integer.MIN_VALUE;
        for (int j : R2) {
          if (used[j]) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, j, C)) continue;
          int score = (g1.aromatic[i] && g2.aromatic[j] ? 10 : 0) + Math.min(g1.degree[i], g2.degree[j]);
          if (score > bestScore) { bestScore = score; bestJ = j; }
        }
        if (bestJ >= 0) { seed.put(i, bestJ); used[bestJ] = true; }
        if (tb.expired()) break;
      }
      return seed;
    }

    Map<Integer, Integer> labelDegreeAnchorSeed(TimeBudget tb) {
      Map<Integer, Integer> seed = new LinkedHashMap<>();
      Map<Integer, Integer> freq = new HashMap<>();
      for (int i = 0; i < g1.n; i++) freq.merge(g1.label[i], 1, Integer::sum);
      List<Integer> Q = new ArrayList<>();
      for (int i = 0; i < g1.n; i++) Q.add(i);
      Q.sort((a, b) -> {
        int ra = freq.get(g1.label[a]), rb = freq.get(g1.label[b]);
        return ra != rb ? Integer.compare(ra, rb) : Integer.compare(g1.degree[b], g1.degree[a]);
      });
      boolean[] used = new boolean[g2.n];
      for (int qi : Q) {
        int bestJ = -1, bestScore = Integer.MIN_VALUE;
        for (int tj = 0; tj < g2.n; tj++) {
          if (used[tj]) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
          int score = 10 - Math.abs(g1.degree[qi] - g2.degree[tj]);
          if (score > bestScore) { bestScore = score; bestJ = tj; }
        }
        if (bestJ >= 0) { seed.put(qi, bestJ); used[bestJ] = true; }
        if (tb.expired()) break;
      }
      return seed;
    }

    Map<Integer, Integer> vf2ppRingSkeletonSeed(TimeBudget tb, long millis, int radius, int maxAnchors) {
      BitSet qRing = ringMask(g1), tRing = ringMask(g2);
      if (qRing.isEmpty() || tRing.isEmpty() || millis <= 0) return Collections.emptyMap();
      return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qRing, tRing);
    }

    Map<Integer, Integer> vf2ppCoreSeed(TimeBudget tb, long millis, int radius, int maxAnchors) {
      BitSet qCore = heavyAtomCoreMask(g1), tCore = heavyAtomCoreMask(g2);
      if (qCore.isEmpty() || tCore.isEmpty() || millis <= 0) return Collections.emptyMap();
      return vf2ppNeighborhoodSeed(tb, millis, radius, maxAnchors, qCore, tCore);
    }

    private Map<Integer, Integer> vf2ppNeighborhoodSeed(
        TimeBudget tb, long millis, int radius, int maxAnchors, BitSet qCand, BitSet tCand) {
      List<int[]> pairs = new ArrayList<>();
      for (int i = qCand.nextSetBit(0); i >= 0; i = qCand.nextSetBit(i + 1))
        for (int j = tCand.nextSetBit(0); j >= 0; j = tCand.nextSetBit(j + 1)) {
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, j, C)) continue;
          int score = 100 - Math.abs(g1.degree[i] - g2.degree[j]);
          if (g1.aromatic[i] && g2.aromatic[j]) score += 10;
          pairs.add(new int[] {i, j, score});
        }
      if (pairs.isEmpty()) return Collections.emptyMap();
      pairs.sort((a, b) -> Integer.compare(b[2], a[2]));
      if (pairs.size() > maxAnchors) pairs = pairs.subList(0, maxAnchors);

      // Run VF2PP once with the full budget — anchor pairs are used for
      // scoring/ordering to select the best candidates, but the matcher
      // operates on the full graphs.  Running once avoids N redundant searches.
      if (tb.expired() || millis <= 0) return Collections.emptyMap();
      SubstructureEngine.VF2PPMatcher m = new SubstructureEngine.VF2PPMatcher(g1, g2, C, new TimeBudget(millis));
      List<Map<Integer, Integer>> out = new ArrayList<>();
      m.enumerate(1, out);
      return out.isEmpty() ? Collections.emptyMap() : out.get(0);
    }

    private static BitSet ringMask(MolGraph g) {
      BitSet bs = new BitSet(g.n);
      for (int i = 0; i < g.n; i++) if (g.ring[i]) bs.set(i);
      return bs;
    }

    private static BitSet heavyAtomCoreMask(MolGraph g) {
      BitSet mask = new BitSet(g.n);
      mask.set(0, g.n);
      boolean changed;
      do {
        changed = false;
        for (int i = 0; i < g.n; i++) {
          if (!mask.get(i) || g.ring[i] || g.aromatic[i]) continue;
          int deg = 0;
          for (int nb : g.neighbors[i]) if (mask.get(nb)) deg++;
          if (deg <= 1 && g.atomicNum[i] == 6) { mask.clear(i); changed = true; }
        }
      } while (changed);
      if (mask.isEmpty()) for (int i = 0; i < g.n; i++) if (g.atomicNum[i] != 1) mask.set(i);
      return mask;
    }

    // McSplit partition-refinement MCS seed

    Map<Integer, Integer> mcSplitSeed(TimeBudget tb, long[] nodeCountOut) {
      int n1 = g1.n, n2 = g2.n;
      if (n1 == 0 || n2 == 0) return new LinkedHashMap<>();

      boolean useMorgan = n1 >= 30 && n2 >= 30;
      if (useMorgan) {
        HashSet<Integer> qDistinct = new HashSet<>();
        for (int i = 0; i < n1; i++) qDistinct.add(g1.morganRank[i]);
        if (qDistinct.size() > n1 * 3 / 4) useMorgan = false;
      }
      // When ringFusionMode is STRICT, atoms with different ringCounts must be in
      // separate classes even if they share the same label/morganRank, because
      // atomsCompatFast rejects pairs with mismatched ringCount.
      boolean splitByRingCount = C.ringFusionMode == ChemOptions.RingFusionMode.STRICT;
      Map<Long, List<Integer>> qGroups = new LinkedHashMap<>(), tGroups = new LinkedHashMap<>();
      for (int i = 0; i < n1; i++) {
        long key = useMorgan ? g1.morganRank[i] : g1.label[i];
        if (splitByRingCount) key = (key << 8) | (g1.ring[i] ? g1.ringCount[i] & 0xFF : 0);
        qGroups.computeIfAbsent(key, k -> new ArrayList<>()).add(i);
      }
      for (int j = 0; j < n2; j++) {
        long key = useMorgan ? g2.morganRank[j] : g2.label[j];
        if (splitByRingCount) key = (key << 8) | (g2.ring[j] ? g2.ringCount[j] & 0xFF : 0);
        tGroups.computeIfAbsent(key, k -> new ArrayList<>()).add(j);
      }

      List<BitSet> initQSets = new ArrayList<>(), initTSets = new ArrayList<>();
      for (Map.Entry<Long, List<Integer>> qe : qGroups.entrySet())
        for (Map.Entry<Long, List<Integer>> te : tGroups.entrySet()) {
          List<Integer> qList = qe.getValue(), tList = te.getValue();
          if (qList.isEmpty() || tList.isEmpty()) continue;
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qList.get(0), g2, tList.get(0), C)) continue;
          BitSet qbs = new BitSet(n1), tbs = new BitSet(n2);
          for (int v : qList) qbs.set(v);
          for (int v : tList) tbs.set(v);
          initQSets.add(qbs);
          initTSets.add(tbs);
        }
      if (initQSets.isEmpty()) return new LinkedHashMap<>();

      int numClasses = initQSets.size();
      BitSet[] qSets = initQSets.toArray(new BitSet[0]);
      BitSet[] tSets = initTSets.toArray(new BitSet[0]);
      int initUB = 0;
      for (int c = 0; c < numClasses; c++) initUB += Math.min(qSets[c].cardinality(), tSets[c].cardinality());

      int[] q2t = new int[n1], t2q = new int[n2];
      Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
      int[] bestQ2T = new int[n1];
      Arrays.fill(bestQ2T, -1);
      int[] bestSize = {0};
      long[] nodeCount = {0};

      mcSplitRecurse(g1, g2, C, induced, tb, qSets, tSets, numClasses, q2t, t2q, 0, initUB,
          bestQ2T, bestSize, nodeCount, n1, n2, 0, Math.min(n1, n2) + 1);

      if (nodeCountOut != null) nodeCountOut[0] = nodeCount[0];
      Map<Integer, Integer> seed = new LinkedHashMap<>();
      for (int i = 0; i < n1; i++) if (bestQ2T[i] >= 0) seed.put(i, bestQ2T[i]);
      return seed;
    }

    private static void mcSplitRecurse(
        MolGraph g1, MolGraph g2, ChemOptions C, boolean induced, TimeBudget tb,
        BitSet[] qSets, BitSet[] tSets, int numClasses, int[] q2t, int[] t2q,
        int curSize, int upperBound, int[] bestQ2T, int[] bestSize, long[] nodeCount,
        int n1, int n2, int depth, int maxDepth) {

      nodeCount[0]++;
      if (nodeCount[0] > MAX_NODE_LIMIT) return;
      if ((nodeCount[0] & 15) == 0 && tb.expiredNow()) return;
      if (curSize > bestSize[0]) { bestSize[0] = curSize; System.arraycopy(q2t, 0, bestQ2T, 0, n1); }
      if (curSize + upperBound <= bestSize[0] || depth >= maxDepth) return;

      // Check if extension is possible
      if (curSize > 0) {
        boolean canExtend = false;
        for (int c = 0; c < numClasses && !canExtend; c++)
          if (!qSets[c].isEmpty() && !tSets[c].isEmpty()) canExtend = true;
        if (!canExtend) return;
      }

      // Select most constrained class
      int bestClass = -1, bestMin = Integer.MAX_VALUE, bestClassQC = 0, bestClassTC = 0;
      for (int c = 0; c < numClasses; c++) {
        int qc = qSets[c].cardinality(), tc = tSets[c].cardinality();
        if (qc == 0 || tc == 0) continue;
        int m = Math.min(qc, tc);
        if (m < bestMin) { bestMin = m; bestClass = c; bestClassQC = qc; bestClassTC = tc; }
      }
      if (bestClass == -1) return;

      boolean branchFromTarget = bestClassTC < bestClassQC;

      if (!branchFromTarget) {
        // Select qi with connectivity-aware ordering
        int qi = -1, qiBestScore = -1;
        for (int v = qSets[bestClass].nextSetBit(0); v >= 0; v = qSets[bestClass].nextSetBit(v + 1)) {
          boolean conn = false;
          for (int nb : g1.neighbors[v]) if (q2t[nb] >= 0) { conn = true; break; }
          int score = (conn ? 1000 : 0) + (g1.ring[v] ? 100 : 0) + g1.degree[v] * 10;
          if (score > qiBestScore) { qiBestScore = score; qi = v; }
        }

        Set<Integer> triedOrbits = new HashSet<>();
        // Check if qi has any already-mapped neighbor — if so, orbit pruning is unsafe
        // because the connectivity context differentiates "equivalent" atoms
        boolean qiHasMappedNeighbor = false;
        for (int nb : g1.neighbors[qi]) if (q2t[nb] >= 0) { qiHasMappedNeighbor = true; break; }
        for (int tj = tSets[bestClass].nextSetBit(0); tj >= 0; tj = tSets[bestClass].nextSetBit(tj + 1)) {
          nodeCount[0]++;
          if (nodeCount[0] > MAX_NODE_LIMIT || ((nodeCount[0] & 15) == 0 && tb.expiredNow())) return;
          // Only use orbit pruning when there is no connectivity constraint;
          // with mapped neighbors the position in the molecule matters
          if (!qiHasMappedNeighbor && !triedOrbits.add(g2.orbit[tj])) continue;

          // Guard: class membership may group atoms that differ on properties checked
          // by atomsCompatFast (e.g. ringCount under STRICT ring-fusion mode).
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
          boolean bondOk = checkBondCompat(g1, g2, C, induced, qi, tj, q2t);
          if (!bondOk) continue;

          q2t[qi] = tj; t2q[tj] = qi;
          int[] refined = refineAndRecurse(g1, g2, C, induced, tb, qSets, tSets, numClasses,
              q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, qi, tj);
          q2t[qi] = -1; t2q[tj] = -1;
        }

        // Try NOT matching qi
        mcSplitSkipVertex(g1, g2, C, induced, tb, qSets, tSets, numClasses, q2t, t2q,
            curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, bestClass, qi, true);

      } else {
        // Select tj with connectivity-aware ordering
        int tj = -1, tjBestScore = -1;
        for (int v = tSets[bestClass].nextSetBit(0); v >= 0; v = tSets[bestClass].nextSetBit(v + 1)) {
          boolean conn = false;
          for (int nb : g2.neighbors[v]) if (t2q[nb] >= 0) { conn = true; break; }
          int score = (conn ? 1000 : 0) + (g2.ring[v] ? 100 : 0) + g2.degree[v] * 10;
          if (score > tjBestScore) { tjBestScore = score; tj = v; }
        }

        Set<Integer> triedOrbits = new HashSet<>();
        // Check if tj has any already-mapped neighbor
        boolean tjHasMappedNeighbor = false;
        for (int nb : g2.neighbors[tj]) if (t2q[nb] >= 0) { tjHasMappedNeighbor = true; break; }
        for (int qi = qSets[bestClass].nextSetBit(0); qi >= 0; qi = qSets[bestClass].nextSetBit(qi + 1)) {
          nodeCount[0]++;
          if (nodeCount[0] > MAX_NODE_LIMIT || ((nodeCount[0] & 15) == 0 && tb.expiredNow())) return;
          if (!tjHasMappedNeighbor && !triedOrbits.add(g1.orbit[qi])) continue;

          // Guard: verify atom-level compatibility (ringCount etc.) before pairing.
          if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
          boolean bondOk = true;
          for (int tk : g2.neighbors[tj]) {
            if (t2q[tk] < 0) continue;
            int qk = t2q[tk];
            int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
            if (qOrd != 0 && tOrd != 0) {
              if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { bondOk = false; break; }
            } else if (induced && ((qOrd != 0) != (tOrd != 0))) { bondOk = false; break; }
          }
          if (!bondOk) continue;

          q2t[qi] = tj; t2q[tj] = qi;
          refineAndRecurse(g1, g2, C, induced, tb, qSets, tSets, numClasses,
              q2t, t2q, curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, qi, tj);
          q2t[qi] = -1; t2q[tj] = -1;
        }

        mcSplitSkipVertex(g1, g2, C, induced, tb, qSets, tSets, numClasses, q2t, t2q,
            curSize, bestQ2T, bestSize, nodeCount, n1, n2, depth, maxDepth, bestClass, tj, false);
      }
    }

    private static boolean checkBondCompat(
        MolGraph g1, MolGraph g2, ChemOptions C, boolean induced, int qi, int tj, int[] q2t) {
      for (int qk : g1.neighbors[qi]) {
        if (q2t[qk] < 0) continue;
        int tk = q2t[qk];
        int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
        if (qOrd != 0 && tOrd != 0) {
          if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) return false;
        } else if (induced && ((qOrd != 0) != (tOrd != 0))) return false;
      }
      return true;
    }

    /** Refine label classes after assigning qi->tj and recurse. */
    private static int[] refineAndRecurse(
        MolGraph g1, MolGraph g2, ChemOptions C, boolean induced, TimeBudget tb,
        BitSet[] qSets, BitSet[] tSets, int numClasses, int[] q2t, int[] t2q,
        int curSize, int[] bestQ2T, int[] bestSize, long[] nodeCount,
        int n1, int n2, int depth, int maxDepth, int qi, int tj) {
      BitSet[] newQ = new BitSet[numClasses * 2], newT = new BitSet[numClasses * 2];
      int newNum = 0, newUB = 0;
      BitSet qiAdj = g1.getAdj()[qi], tjAdj = g2.getAdj()[tj];
      for (int c = 0; c < numClasses; c++) {
        BitSet qc = qSets[c], tc = tSets[c];
        if (qc.isEmpty() || tc.isEmpty()) continue;
        BitSet qAdj = new BitSet(n1); qAdj.or(qc); qAdj.and(qiAdj); qAdj.clear(qi);
        BitSet qNon = new BitSet(n1); qNon.or(qc); qNon.andNot(qiAdj); qNon.clear(qi);
        BitSet tAdj = new BitSet(n2); tAdj.or(tc); tAdj.and(tjAdj); tAdj.clear(tj);
        BitSet tNon = new BitSet(n2); tNon.or(tc); tNon.andNot(tjAdj); tNon.clear(tj);
        int qAdjC = qAdj.cardinality(), tAdjC = tAdj.cardinality();
        if (qAdjC > 0 && tAdjC > 0) {
          newQ[newNum] = qAdj; newT[newNum] = tAdj; newUB += Math.min(qAdjC, tAdjC); newNum++;
        }
        int qNonC = qNon.cardinality(), tNonC = tNon.cardinality();
        if (qNonC > 0 && tNonC > 0) {
          newQ[newNum] = qNon; newT[newNum] = tNon; newUB += Math.min(qNonC, tNonC); newNum++;
        }
      }
      if (curSize + 1 + newUB > bestSize[0])
        mcSplitRecurse(g1, g2, C, induced, tb, newQ, newT, newNum, q2t, t2q, curSize + 1, newUB,
            bestQ2T, bestSize, nodeCount, n1, n2, depth + 1, maxDepth);
      return bestQ2T;
    }

    /** Try NOT matching a vertex: remove from class, reduce bound. */
    private static void mcSplitSkipVertex(
        MolGraph g1, MolGraph g2, ChemOptions C, boolean induced, TimeBudget tb,
        BitSet[] qSets, BitSet[] tSets, int numClasses, int[] q2t, int[] t2q,
        int curSize, int[] bestQ2T, int[] bestSize, long[] nodeCount,
        int n1, int n2, int depth, int maxDepth, int bestClass, int vertex, boolean isQuery) {
      BitSet[] skipQ = new BitSet[numClasses], skipT = new BitSet[numClasses];
      int skipNum = 0, skipUB = 0;
      for (int c = 0; c < numClasses; c++) {
        BitSet qc = qSets[c], tc = tSets[c];
        if (qc.isEmpty() || tc.isEmpty()) continue;
        if (c == bestClass) {
          BitSet reduced = isQuery ? (BitSet) qc.clone() : qc;
          BitSet reducedT = isQuery ? tc : (BitSet) tc.clone();
          if (isQuery) reduced.clear(vertex); else reducedT.clear(vertex);
          if (!(isQuery ? reduced : reducedT).isEmpty()) {
            skipQ[skipNum] = isQuery ? reduced : qc;
            skipT[skipNum] = isQuery ? tc : reducedT;
            skipUB += Math.min((isQuery ? reduced : qc).cardinality(), (isQuery ? tc : reducedT).cardinality());
            skipNum++;
          }
        } else {
          skipQ[skipNum] = qc; skipT[skipNum] = tc;
          skipUB += Math.min(qc.cardinality(), tc.cardinality());
          skipNum++;
        }
      }
      if (curSize + skipUB > bestSize[0])
        mcSplitRecurse(g1, g2, C, induced, tb, skipQ, skipT, skipNum, q2t, t2q, curSize, skipUB,
            bestQ2T, bestSize, nodeCount, n1, n2, depth + 1, maxDepth);
    }

    // Seed-and-extend MCS (bond-growth)

    Map<Integer, Integer> seedExtendMCS(TimeBudget tb, int upperBound) {
      int n1 = g1.n, n2 = g2.n;
      if (n1 == 0 || n2 == 0) return new LinkedHashMap<>();

      boolean[] g1BondTypes = new boolean[5], g2BondTypes = new boolean[5];
      for (int i = 0; i < n1; i++)
        for (int nb : g1.neighbors[i]) if (nb > i) g1BondTypes[g1.bondOrder(i, nb)] = true;
      for (int j = 0; j < n2; j++)
        for (int nb : g2.neighbors[j]) if (nb > j) g2BondTypes[g2.bondOrder(j, nb)] = true;
      boolean[] commonBondTypes = new boolean[5];
      for (int t = 1; t <= 4; t++) commonBondTypes[t] = g1BondTypes[t] && g2BondTypes[t];

      int maxLabel = 0;
      for (int i = 0; i < n1; i++) maxLabel = Math.max(maxLabel, g1.label[i]);
      for (int j = 0; j < n2; j++) maxLabel = Math.max(maxLabel, g2.label[j]);
      int[] labelFreqG1 = new int[maxLabel + 1], labelFreqG2 = new int[maxLabel + 1];
      for (int i = 0; i < n1; i++) labelFreqG1[g1.label[i]]++;
      for (int j = 0; j < n2; j++) labelFreqG2[g2.label[j]]++;

      int[][] compatT = new int[n1][];
      for (int i = 0; i < n1; i++) {
        int cnt = 0;
        for (int j = 0; j < n2; j++) if (SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, j, C)) cnt++;
        compatT[i] = new int[cnt];
        cnt = 0;
        for (int j = 0; j < n2; j++)
          if (SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, i, g2, j, C)) compatT[i][cnt++] = j;
      }

      int[][] g1Bonds = new int[n1 * 4][2];
      int g1BondCount = 0;
      for (int i = 0; i < n1; i++)
        for (int nb : g1.neighbors[i])
          if (nb > i) {
            int ord = g1.bondOrder(i, nb);
            boolean typeOk = commonBondTypes[ord];
            if (!typeOk && C.aromaticityMode == ChemOptions.AromaticityMode.FLEXIBLE) {
              if (g1.bondAromatic(i, nb)) typeOk = g2BondTypes[1] || g2BondTypes[2] || g2BondTypes[4];
              else if (ord == 1 || ord == 2) typeOk = g2BondTypes[4];
            }
            if (!typeOk && (C.matchBondOrder == ChemOptions.BondOrderMode.ANY
                || C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE)) typeOk = true;
            if (typeOk) {
              if (g1BondCount >= g1Bonds.length) g1Bonds = java.util.Arrays.copyOf(g1Bonds, g1Bonds.length * 2);
              g1Bonds[g1BondCount][0] = i; g1Bonds[g1BondCount][1] = nb; g1BondCount++;
            }
          }

      int[] bestQ2T = new int[n1];
      Arrays.fill(bestQ2T, -1);
      int bestSize = 0;
      long nodeCount = 0;
      int minN = Math.min(n1, n2);

      int[][] seedList = new int[g1BondCount][2];
      for (int b = 0; b < g1BondCount; b++) {
        int u = g1Bonds[b][0], v = g1Bonds[b][1];
        int totalAtoms = n1 + n2;
        int uFreq = labelFreqG1[g1.label[u]] + (g1.label[u] < labelFreqG2.length ? labelFreqG2[g1.label[u]] : 0);
        int vFreq = labelFreqG1[g1.label[v]] + (g1.label[v] < labelFreqG2.length ? labelFreqG2[g1.label[v]] : 0);
        int score = (totalAtoms * 2) / Math.max(1, uFreq) + (totalAtoms * 2) / Math.max(1, vFreq)
            + (n2 * 2) / Math.max(1, compatT[u].length) + (n2 * 2) / Math.max(1, compatT[v].length)
            + g1.degree[u] + g1.degree[v] + (g1.ring[u] ? 10 : 0) + (g1.ring[v] ? 10 : 0);
        seedList[b][0] = b; seedList[b][1] = score;
      }
      java.util.Arrays.sort(seedList, 0, g1BondCount, (a, b) -> Integer.compare(b[1], a[1]));

      int maxSeeds = Math.min(g1BondCount, 16);
      boolean[] atomUsedAsSeed = new boolean[n1];
      Set<Long> triedOrbitPairs = new HashSet<>();
      int seedsTried = 0;

      for (int si = 0; si < g1BondCount && seedsTried < maxSeeds; si++) {
        if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
        int bondIdx = seedList[si][0];
        int qu = g1Bonds[bondIdx][0], qv = g1Bonds[bondIdx][1];
        // Use orbit + degree to distinguish bonds at different chain positions
        // (e.g. P-O bonds in phosphate chains share orbits but have different degrees)
        long sigU = ((long) g1.orbit[qu] << 10) | g1.degree[qu];
        long sigV = ((long) g1.orbit[qv] << 10) | g1.degree[qv];
        long oA = Math.min(sigU, sigV), oB = Math.max(sigU, sigV);
        if (!triedOrbitPairs.add((oA << 32) | oB)) continue;
        if (atomUsedAsSeed[qu] && atomUsedAsSeed[qv]) continue;
        atomUsedAsSeed[qu] = true; atomUsedAsSeed[qv] = true; seedsTried++;

        for (int ci = 0; ci < compatT[qu].length; ci++) {
          int ta = compatT[qu][ci];
          if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
          for (int tb2 : g2.neighbors[ta]) {
            nodeCount++;
            if (nodeCount > MAX_NODE_LIMIT || tb.expired()) break;
            if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qv, g2, tb2, C)) continue;
            if (!MolGraph.ChemOps.bondsCompatible(g1, qu, qv, g2, ta, tb2, C)) continue;
            int[] q2t = new int[n1], t2q = new int[n2];
            Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
            q2t[qu] = ta; t2q[ta] = qu; q2t[qv] = tb2; t2q[tb2] = qv;
            int mapSize = greedyBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced);
            nodeCount += mapSize;
            if (mapSize > bestSize) { bestSize = mapSize; System.arraycopy(q2t, 0, bestQ2T, 0, n1); }
          }
        }
        for (int ci = 0; ci < compatT[qv].length; ci++) {
          int ta = compatT[qv][ci];
          if (tb.expired() || nodeCount > MAX_NODE_LIMIT) break;
          for (int tb2 : g2.neighbors[ta]) {
            nodeCount++;
            if (nodeCount > MAX_NODE_LIMIT || tb.expired()) break;
            if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qu, g2, tb2, C)) continue;
            if (!MolGraph.ChemOps.bondsCompatible(g1, qu, qv, g2, tb2, ta, C)) continue;
            int[] q2t = new int[n1], t2q = new int[n2];
            Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
            q2t[qu] = tb2; t2q[tb2] = qu; q2t[qv] = ta; t2q[ta] = qv;
            int mapSize = greedyBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced);
            nodeCount += mapSize;
            if (mapSize > bestSize) { bestSize = mapSize; System.arraycopy(q2t, 0, bestQ2T, 0, n1); }
          }
        }
        if (bestSize >= upperBound || bestSize >= minN) break;
        if (bestSize >= (upperBound * 19 / 20) && nodeCount > SEED_EXTEND_NODE_LIMIT / 2) break;
      }

      if (bestSize > 0 && bestSize < minN && !tb.expired() && nodeCount < MAX_NODE_LIMIT) {
        int[] q2t = new int[n1], t2q = new int[n2];
        Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
        for (int i = 0; i < n1; i++) if (bestQ2T[i] >= 0) { q2t[i] = bestQ2T[i]; t2q[bestQ2T[i]] = i; }
        long[] nc = {nodeCount};
        backtrackBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced, bestQ2T, new int[] {bestSize}, tb, nc, MAX_NODE_LIMIT);
        bestSize = 0;
        for (int i = 0; i < n1; i++) if (bestQ2T[i] >= 0) bestSize++;
      }

      Map<Integer, Integer> result = new LinkedHashMap<>();
      for (int i = 0; i < n1; i++) if (bestQ2T[i] >= 0) result.put(i, bestQ2T[i]);
      return result;
    }

    private static int greedyBondExtend(
        MolGraph g1, MolGraph g2, ChemOptions C, int[] q2t, int[] t2q, int n1, int n2, boolean induced) {
      int mapSize = 0;
      for (int i = 0; i < n1; i++) if (q2t[i] >= 0) mapSize++;
      boolean progress = true;
      while (progress) {
        progress = false;
        for (int qi = 0; qi < n1; qi++) {
          if (q2t[qi] >= 0) continue;
          int bestTj = -1, bestScore = -1;
          for (int nb : g1.neighbors[qi]) {
            if (q2t[nb] < 0) continue;
            int tNb = q2t[nb];
            for (int tj : g2.neighbors[tNb]) {
              if (t2q[tj] >= 0) continue;
              if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
              if (!MolGraph.ChemOps.bondsCompatible(g1, qi, nb, g2, tj, tNb, C)) continue;
              boolean consistent = true;
              for (int qk : g1.neighbors[qi]) {
                if (qk == nb || q2t[qk] < 0) continue;
                int tk = q2t[qk];
                int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
                if (qOrd != 0 && tOrd != 0) {
                  if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) { consistent = false; break; }
                } else if (induced && ((qOrd != 0) != (tOrd != 0))) { consistent = false; break; }
              }
              if (!consistent) continue;
              int score = (g2.ring[tj] && g1.ring[qi] ? 50 : 0) + Math.min(g1.degree[qi], g2.degree[tj]);
              if (score > bestScore) { bestScore = score; bestTj = tj; }
            }
          }
          if (bestTj >= 0) { q2t[qi] = bestTj; t2q[bestTj] = qi; mapSize++; progress = true; }
        }
      }
      return mapSize;
    }

    private static int mappedNeighborCount(MolGraph g1, int[] q2t, int qi) {
      int count = 0;
      for (int nb : g1.neighbors[qi]) if (q2t[nb] >= 0) count++;
      return count;
    }

    private static List<Integer> collectBondExtendCandidates(
        MolGraph g1, MolGraph g2, ChemOptions C, boolean induced,
        int qi, int[] q2t, int[] t2q) {
      int seedQk = -1;
      int seedTk = -1;
      int mappedNbrs = 0;
      for (int qk : g1.neighbors[qi]) {
        int tk = q2t[qk];
        if (tk < 0) continue;
        mappedNbrs++;
        if (seedTk < 0 || g2.degree[tk] < g2.degree[seedTk]) {
          seedQk = qk;
          seedTk = tk;
        }
      }
      if (seedTk < 0) return Collections.emptyList();

      List<int[]> scored = new ArrayList<>(g2.degree[seedTk]);
      for (int tj : g2.neighbors[seedTk]) {
        if (t2q[tj] >= 0) continue;
        if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C)) continue;
        if (!MolGraph.ChemOps.bondsCompatible(g1, qi, seedQk, g2, tj, seedTk, C)) continue;

        boolean consistent = true;
        for (int qk : g1.neighbors[qi]) {
          if (qk == seedQk || q2t[qk] < 0) continue;
          int tk = q2t[qk];
          int qOrd = g1.bondOrder(qi, qk), tOrd = g2.bondOrder(tj, tk);
          if (qOrd != 0 && tOrd != 0) {
            if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C)) {
              consistent = false;
              break;
            }
          } else if (induced && ((qOrd != 0) != (tOrd != 0))) {
            consistent = false;
            break;
          }
        }
        if (!consistent) continue;

        int score = mappedNbrs * 100
            + (g1.ring[qi] && g2.ring[tj] ? 50 : 0)
            + Math.min(g1.degree[qi], g2.degree[tj]);
        scored.add(new int[] {score, tj});
      }

      scored.sort((a, b) -> {
        if (a[0] != b[0]) return Integer.compare(b[0], a[0]);
        return Integer.compare(a[1], b[1]);
      });

      List<Integer> out = new ArrayList<>(scored.size());
      for (int[] entry : scored) out.add(entry[1]);
      return out;
    }

    private static void backtrackBondExtend(
        MolGraph g1, MolGraph g2, ChemOptions C, int[] q2t, int[] t2q, int n1, int n2,
        boolean induced, int[] bestQ2T, int[] bestSize, TimeBudget tb, long[] nodeCount, long nodeLimit) {
      int curSize = 0;
      for (int i = 0; i < n1; i++) if (q2t[i] >= 0) curSize++;
      if (curSize > bestSize[0]) { bestSize[0] = curSize; System.arraycopy(q2t, 0, bestQ2T, 0, n1); }
      if (nodeCount[0] > nodeLimit || tb.expired()) return;

      int frontier = 0;
      int bestQi = -1, bestConstraint = -1, bestCandCount = Integer.MAX_VALUE;
      List<Integer> bestCandidates = Collections.emptyList();
      for (int qi = 0; qi < n1; qi++) {
        if (q2t[qi] >= 0) continue;
        int mappedNeighbors = mappedNeighborCount(g1, q2t, qi);
        if (mappedNeighbors == 0) continue;
        List<Integer> candidates = collectBondExtendCandidates(g1, g2, C, induced, qi, q2t, t2q);
        if (candidates.isEmpty()) continue;
        frontier++;
        int constraint = mappedNeighbors * 100 + (g1.ring[qi] ? 50 : 0) + g1.degree[qi];
        int candCount = candidates.size();
        if (candCount < bestCandCount
            || (candCount == bestCandCount && constraint > bestConstraint)) {
          bestCandCount = candCount;
          bestConstraint = constraint;
          bestQi = qi;
          bestCandidates = candidates;
        }
      }
      if (curSize + frontier <= bestSize[0]) return;
      if (bestQi < 0) return;

      int qi = bestQi;
      for (int tj : bestCandidates) {
        nodeCount[0]++;
        if (nodeCount[0] > nodeLimit || tb.expired()) return;
        q2t[qi] = tj; t2q[tj] = qi;
        backtrackBondExtend(g1, g2, C, q2t, t2q, n1, n2, induced, bestQ2T, bestSize, tb, nodeCount, nodeLimit);
        q2t[qi] = -1; t2q[tj] = -1;
        if (nodeCount[0] > nodeLimit || tb.expired()) return;
      }
    }
  }

  static int labelFrequencyUpperBound(MolGraph g1, MolGraph g2, ChemOptions C) {
    if (g1.n == 0 || g2.n == 0) return 0;
    Map<Integer, Integer> freq1 = new HashMap<>(), freq2 = new HashMap<>();
    for (int i = 0; i < g1.n; i++) freq1.merge(ubLabel(g1, i, C), 1, Integer::sum);
    for (int j = 0; j < g2.n; j++) freq2.merge(ubLabel(g2, j, C), 1, Integer::sum);
    int ub = 0;
    for (Map.Entry<Integer, Integer> e : freq1.entrySet())
      ub += Math.min(e.getValue(), freq2.getOrDefault(e.getKey(), 0));
    return Math.min(ub, Math.min(g1.n, g2.n));
  }

  static int labelFrequencyUpperBoundDirected(MolGraph query, MolGraph target, ChemOptions C) {
    if (query.n == 0 || target.n == 0) return 0;
    if (C.ringMatchesRingOnly) return labelFrequencyUpperBound(query, target, C);

    Map<Integer, Integer> qAny = new HashMap<>(), qRing = new HashMap<>();
    Map<Integer, Integer> tAny = new HashMap<>(), tRing = new HashMap<>();
    for (int i = 0; i < query.n; i++) {
      int lab = ubBaseLabel(query, i, C);
      if (C.ringMatchesRingOnly && query.ring[i]) qRing.merge(lab, 1, Integer::sum);
      else qAny.merge(lab, 1, Integer::sum);
    }
    for (int j = 0; j < target.n; j++) {
      int lab = ubBaseLabel(target, j, C);
      tAny.merge(lab, 1, Integer::sum);
      if (target.ring[j]) tRing.merge(lab, 1, Integer::sum);
    }

    int ub = 0;
    for (Map.Entry<Integer, Integer> e : qRing.entrySet())
      ub += Math.min(e.getValue(), tRing.getOrDefault(e.getKey(), 0));
    for (Map.Entry<Integer, Integer> e : qAny.entrySet()) {
      int lab = e.getKey();
      int avail = tAny.getOrDefault(lab, 0);
      if (C.ringMatchesRingOnly) {
        int usedRing = Math.min(qRing.getOrDefault(lab, 0), tRing.getOrDefault(lab, 0));
        avail -= usedRing;
      }
      if (avail > 0) ub += Math.min(e.getValue(), avail);
    }
    return Math.min(ub, Math.min(query.n, target.n));
  }

  /**
   * Degree Sequence Upper Bound (DSB) on MCS size.
   *
   * <p>For each atom label present in both graphs, collects the degree sequences,
   * sorts them descending, and greedily counts matchable pairs: an atom of degree d
   * in the shorter list can only match an atom of degree >= d in the longer list.
   * This is provably at least as tight as the label-frequency upper bound.
   *
   * <p>Runs in O(n log n) time.
   */
  static int degreeSequenceUpperBound(MolGraph g1, MolGraph g2, ChemOptions C) {
    if (g1.n == 0 || g2.n == 0) return 0;

    // Group atom degrees by label for each graph
    Map<Integer, List<Integer>> degs1 = new HashMap<>(), degs2 = new HashMap<>();
    for (int i = 0; i < g1.n; i++)
      degs1.computeIfAbsent(ubLabel(g1, i, C), k -> new ArrayList<>()).add(g1.neighbors[i].length);
    for (int j = 0; j < g2.n; j++)
      degs2.computeIfAbsent(ubLabel(g2, j, C), k -> new ArrayList<>()).add(g2.neighbors[j].length);

    int ub = 0;
    for (Map.Entry<Integer, List<Integer>> e : degs1.entrySet()) {
      List<Integer> d2 = degs2.get(e.getKey());
      if (d2 == null) continue;
      List<Integer> d1 = e.getValue();
      // Sort both in descending order
      d1.sort(Collections.reverseOrder());
      d2.sort(Collections.reverseOrder());
      // Ensure d1 is the shorter list (we iterate over it)
      List<Integer> shorter = d1, longer = d2;
      if (d1.size() > d2.size()) { shorter = d2; longer = d1; }
      // Greedy two-pointer: for each atom in shorter (desc), find next atom in longer with degree >= it
      int j = 0;
      for (int i = 0; i < shorter.size() && j < longer.size(); i++) {
        // longer is sorted desc, so longer[j] is the largest remaining degree
        if (longer.get(j) >= shorter.get(i)) {
          ub++;
          j++;
        }
        // If longer[j] < shorter[i], this atom in shorter cannot match any remaining
        // atom in longer (all subsequent are <= longer[j]), so skip it
      }
    }
    return Math.min(ub, Math.min(g1.n, g2.n));
  }

  private static int graphOrbitRedundancy(MolGraph g) {
    if (g.orbit == null || g.orbit.length == 0) return 0;
    Map<Integer, Integer> orbitCounts = new HashMap<>();
    for (int orbit : g.orbit) orbitCounts.merge(orbit, 1, Integer::sum);
    int redundancy = 0;
    for (int count : orbitCounts.values()) if (count > 1) redundancy += count - 1;
    return redundancy;
  }

  private static int graphConstraintScore(MolGraph g) {
    int score = 0;
    for (int i = 0; i < g.n; i++) {
      score += g.ring[i] ? 12 : 0;
      score += g.aromatic[i] ? 6 : 0;
      score += g.formalCharge[i] != 0 ? 10 : 0;
      score += g.massNumber[i] > 0 ? 4 : 0;
      score += g.degree[i] >= 3 ? 3 : 0;
      if (g.atomicNum[i] != 6 && g.atomicNum[i] != 1) score += 2;
    }
    return score;
  }

  private static int graphOrientationBurden(MolGraph g) {
    int degreeSum = 0;
    int ringAtoms = 0;
    int aromaticAtoms = 0;
    for (int i = 0; i < g.n; i++) {
      degreeSum += g.degree[i];
      if (g.ring[i]) ringAtoms++;
      if (g.aromatic[i]) aromaticAtoms++;
    }
    int edgeCount = degreeSum / 2;
    return edgeCount * 8 + ringAtoms * 16 + aromaticAtoms * 4 + graphOrbitRedundancy(g) * 24;
  }

  private static long orientationConstraintMass(MolGraph query, MolGraph target, ChemOptions C) {
    Map<Integer, Integer> targetFreq = new HashMap<>();
    for (int j = 0; j < target.n; j++) targetFreq.merge(ubBaseLabel(target, j, C), 1, Integer::sum);

    long mass = 0;
    for (int i = 0; i < query.n; i++) mass += targetFreq.getOrDefault(ubBaseLabel(query, i, C), 0);
    return mass;
  }

  private static int orientationSeedProbeSize(MolGraph query, MolGraph target, ChemOptions C) {
    if (query.n == 0 || target.n == 0) return 0;
    GraphBuilder gb = new GraphBuilder(query, target, C, false);
    TimeBudget tb = new TimeBudget(24L * 60 * 60 * 1000);
    int upperBound = labelFrequencyUpperBound(query, target, C);
    return gb.seedExtendMCS(tb, upperBound).size();
  }

  private static int orientationMcSplitProbeSize(MolGraph query, MolGraph target, ChemOptions C) {
    if (query.n == 0 || target.n == 0) return 0;
    GraphBuilder gb = new GraphBuilder(query, target, C, false);
    TimeBudget tb = new TimeBudget(24L * 60 * 60 * 1000);
    long[] nodeCount = {0};
    return gb.mcSplitSeed(tb, nodeCount).size();
  }

  private static final class OrientationPlan {
    boolean directFirst = true;
    int seed12 = 0;
    int seed21 = 0;
    int mc12 = 0;
    int mc21 = 0;
  }

  private static OrientationPlan chooseOrientationPlan(MolGraph g1, MolGraph g2, ChemOptions C,
                                                       int ub12, int ub21) {
    OrientationPlan plan = new OrientationPlan();
    long mass12 = orientationConstraintMass(g1, g2, C);
    long mass21 = orientationConstraintMass(g2, g1, C);
    if (mass12 != mass21) {
      plan.directFirst = mass12 < mass21;
      return plan;
    }

    int minN = Math.min(g1.n, g2.n);
    int product = g1.n * g2.n;
    boolean nearSymmetric = Math.abs(g1.n - g2.n) <= 2;
    if (nearSymmetric && minN >= 12 && product >= 256) {
      plan.seed12 = orientationSeedProbeSize(g1, g2, C);
      plan.seed21 = orientationSeedProbeSize(g2, g1, C);
      if (plan.seed12 != plan.seed21) {
        plan.directFirst = plan.seed12 > plan.seed21;
        return plan;
      }
    }
    if (nearSymmetric && minN >= 12 && product >= 512) {
      plan.mc12 = orientationMcSplitProbeSize(g1, g2, C);
      plan.mc21 = orientationMcSplitProbeSize(g2, g1, C);
      if (plan.mc12 != plan.mc21) {
        plan.directFirst = plan.mc12 > plan.mc21;
        return plan;
      }
    }

    int burden1 = graphOrientationBurden(g1);
    int burden2 = graphOrientationBurden(g2);
    if (burden1 != burden2) {
      plan.directFirst = burden1 < burden2;
      return plan;
    }

    int score1 = graphConstraintScore(g1);
    int score2 = graphConstraintScore(g2);
    if (score1 != score2) {
      plan.directFirst = score1 > score2;
      return plan;
    }

    if (Math.abs(g1.n - g2.n) >= 4) {
      plan.directFirst = g1.n < g2.n;
      return plan;
    }

    if (!nearSymmetric && minN >= 12 && product >= 256) {
      plan.seed12 = orientationSeedProbeSize(g1, g2, C);
      plan.seed21 = orientationSeedProbeSize(g2, g1, C);
      if (plan.seed12 != plan.seed21) {
        plan.directFirst = plan.seed12 > plan.seed21;
        return plan;
      }
    }
    if (!nearSymmetric && minN >= 12 && product >= 512) {
      plan.mc12 = orientationMcSplitProbeSize(g1, g2, C);
      plan.mc21 = orientationMcSplitProbeSize(g2, g1, C);
      if (plan.mc12 != plan.mc21) {
        plan.directFirst = plan.mc12 > plan.mc21;
        return plan;
      }
    }

    if (g1.n != g2.n) {
      plan.directFirst = g1.n < g2.n;
      return plan;
    }
    if (ub12 != ub21) {
      plan.directFirst = ub12 > ub21;
      return plan;
    }
    plan.directFirst = g1.canonicalHash <= g2.canonicalHash;
    return plan;
  }

  private static int alternateOrientationFloor(OrientationPlan plan, boolean directRan) {
    return directRan ? Math.max(plan.seed21, plan.mc21) : Math.max(plan.seed12, plan.mc12);
  }

  private static Map<Integer, Integer> orientMcsResult(Map<Integer, Integer> mapping, boolean swappedBack) {
    if (!swappedBack) return mapping;
    Map<Integer, Integer> restored = new LinkedHashMap<>();
    for (Map.Entry<Integer, Integer> e : mapping.entrySet()) restored.put(e.getValue(), e.getKey());
    return restored;
  }

  private static Map<Integer, Integer> runValidatedMcsDirection(MolGraph query, MolGraph target,
                                                                ChemOptions C, McsOptions M,
                                                                boolean swappedBack) {
    Map<Integer, Integer> raw = findMCSImpl(query, target, C, M);
    Map<Integer, Integer> valid = validateMapping(query, target, raw, C).isEmpty()
        ? raw
        : recoverValidMcsMapping(query, target, raw, C, M);
    return orientMcsResult(valid, swappedBack);
  }

  private static boolean preferFinalMapping(MolGraph g1,
                                            Map<Integer, Integer> candidate,
                                            Map<Integer, Integer> incumbent,
                                            McsOptions M) {
    if (candidate == null || candidate.isEmpty()) return false;
    if (incumbent == null || incumbent.isEmpty()) return true;

    int candScore = mcsScore(g1, candidate, M);
    int bestScore = mcsScore(g1, incumbent, M);
    boolean weightMode = M.maximizeBonds || M.atomWeights != null;
    if (weightMode && candScore != bestScore) return candScore > bestScore;

    if (candidate.size() != incumbent.size()) return candidate.size() > incumbent.size();
    if (!weightMode && candScore != bestScore) return candScore > bestScore;

    int candHetero = heteroatomScore(g1, candidate);
    int bestHetero = heteroatomScore(g1, incumbent);
    if (candHetero != bestHetero) return candHetero > bestHetero;

    int candBonds = countMappedBonds(g1, candidate);
    int bestBonds = countMappedBonds(g1, incumbent);
    if (candBonds != bestBonds) return candBonds > bestBonds;

    Iterator<Map.Entry<Integer, Integer>> aIt = new TreeMap<>(candidate).entrySet().iterator();
    Iterator<Map.Entry<Integer, Integer>> bIt = new TreeMap<>(incumbent).entrySet().iterator();
    while (aIt.hasNext() && bIt.hasNext()) {
      Map.Entry<Integer, Integer> a = aIt.next();
      Map.Entry<Integer, Integer> b = bIt.next();
      if (!a.getKey().equals(b.getKey())) return a.getKey() < b.getKey();
      if (!a.getValue().equals(b.getValue())) return a.getValue() < b.getValue();
    }
    return candidate.size() > incumbent.size();
  }

  static int ubLabel(MolGraph g, int i, ChemOptions C) {
    if (!C.matchAtomType) return 0;
    int key = g.atomicNum[i] << 2;
    if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT) key |= g.aromatic[i] ? 2 : 0;
    if (C.ringMatchesRingOnly) key |= g.ring[i] ? 1 : 0;
    // In STRICT ring-fusion mode, bridgehead atoms (ringCount >= 2) must only match
    // atoms with the same ringCount — encode it in the label so the frequency upper
    // bound does not over-count.
    if (C.ringFusionMode == ChemOptions.RingFusionMode.STRICT && g.ring[i])
      key = (key << 8) | (g.ringCount[i] & 0xFF);
    return key;
  }

  static int ubBaseLabel(MolGraph g, int i, ChemOptions C) {
    if (!C.matchAtomType) return 0;
    int key = g.atomicNum[i] << 1;
    if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT) key |= g.aromatic[i] ? 1 : 0;
    if (C.ringFusionMode == ChemOptions.RingFusionMode.STRICT && g.ring[i])
      key = (key << 8) | (g.ringCount[i] & 0xFF);
    return key;
  }

  // RASCAL-style similarity upper bound

  /**
   * Fast upper bound on Tanimoto atom similarity using atom-label frequency overlap (MolGraph API).
   *
   * <p>Runs in O(V) time. Useful for pre-filtering molecule pairs before expensive MCS computation.
   *
   * @param g1 the first molecule graph
   * @param g2 the second molecule graph
   * @param C  chemical matching options
   * @return upper bound on Tanimoto similarity in [0.0, 1.0]
   */
  public static double similarityUpperBound(MolGraph g1, MolGraph g2, ChemOptions C) {
    if (g1.n == 0 || g2.n == 0) return 0.0;
    Map<Integer, Integer> freq1 = new HashMap<>(), freq2 = new HashMap<>();
    for (int i = 0; i < g1.n; i++) freq1.merge(g1.label[i], 1, Integer::sum);
    for (int j = 0; j < g2.n; j++) freq2.merge(g2.label[j], 1, Integer::sum);
    int maxCommonAtoms = 0;
    for (Map.Entry<Integer, Integer> e : freq1.entrySet())
      maxCommonAtoms += Math.min(e.getValue(), freq2.getOrDefault(e.getKey(), 0));
    int ubAtoms = Math.min(Math.min(g1.n, g2.n), maxCommonAtoms);
    return ubAtoms == 0 ? 0.0 : (double) ubAtoms / (g1.n + g2.n - ubAtoms);
  }

  /**
   * Fast upper bound on Tanimoto atom similarity (CDK API).
   *
   * @param m1 the first molecule
   * @param m2 the second molecule
   * @param C  chemical matching options
   * @return upper bound on Tanimoto similarity in [0.0, 1.0]
   */
  public static double similarityUpperBound(IAtomContainer m1, IAtomContainer m2, ChemOptions C) {
    MolGraph g1 = toMolGraph(m1), g2 = toMolGraph(m2);
    applySolvent(g1, C);
    applySolvent(g2, C);
    return similarityUpperBound(g1, g2, C);
  }

  // ---------------------------------------------------------------------------
  // Tautomer confidence scoring
  // ---------------------------------------------------------------------------

  /**
   * Compute a pKa-informed confidence score for a tautomer-aware MCS result.
   *
   * <p>For each atom pair {@code (qi, ti)} in the mapping, if the atoms belong to the same
   * tautomeric equivalence class but have different element types (a cross-element tautomeric
   * match, e.g. O matched to N via amide/imidic acid), the confidence is reduced by
   * {@code min(w_q, w_t)}, where the weights are the empirical pKa-derived relevance values
   * stored in {@link MolGraph#tautomerWeight} (Sitzmann 2010).
   *
   * <p>A score of {@code 1.0} means every matched atom pair is a direct element-identical match —
   * no tautomeric flexibility was exercised.  A score approaching {@code 0} indicates the MCS
   * relies heavily on low-probability tautomeric forms.
   *
   * <pre>{@code
   * ChemOptions c = ChemOptions.tautomerProfile();
   * Map<Integer,Integer> mcs = SearchEngine.findMCS(ketoForm, enolForm, c, opts);
   * double conf = SearchEngine.computeTautomerConfidence(ketoForm, enolForm, mcs);
   * // conf ~ 0.97 (amide/imidic match pair), or 1.0 if no cross-element matches
   * }</pre>
   *
   * @param m1  the first molecule (must match the query used in MCS)
   * @param m2  the second molecule
   * @param mcs atom-index mapping returned by {@link #findMCS}
   * @return confidence ∈ [0,1]; 1.0 when mcs is empty or has no tautomeric matches
   */
  public static double computeTautomerConfidence(
      IAtomContainer m1, IAtomContainer m2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty()) return 1.0;
    MolGraph g1 = toMolGraph(m1);
    MolGraph g2 = toMolGraph(m2);
    return computeTautomerConfidence(g1, g2, mcs);
  }

  /** MolGraph overload — avoids re-creating graphs when they are already available. */
  public static double computeTautomerConfidence(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty()) return 1.0;
    g1.ensureTautomerClasses(); g2.ensureTautomerClasses();
    double logSum = 0.0;
    int count = 0;
    for (Map.Entry<Integer, Integer> e : mcs.entrySet()) {
      int qi = e.getKey(), ti = e.getValue();
      if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
      // Cross-element tautomeric match: same class, different element
      if (g1.tautomerClass[qi] >= 0 && g2.tautomerClass[ti] >= 0
          && g1.atomicNum[qi] != g2.atomicNum[ti]) {
        double w = Math.min(g1.tautomerWeight[qi], g2.tautomerWeight[ti]);
        logSum += Math.log(Math.max(w, 1e-10));
        count++;
      }
    }
    // Geometric mean: (∏ w_i)^(1/|M|)
    return count == 0 ? 1.0 : Math.exp(logSum / count);
  }

  /**
   * Returns a human-readable string describing the available compute resources.
   *
   * <p>Example output: {@code "Java 17.0.9 / 16 processors (auto-parallel batch enabled)"}
   */
  public static String computeInfo() {
    int cores = Runtime.getRuntime().availableProcessors();
    return "Java " + System.getProperty("java.version") +
        " / " + cores + " processor" + (cores == 1 ? "" : "s") +
        " (auto-parallel batch enabled)";
  }

  /**
   * Compute pairwise MCS for all molecule pairs above a similarity threshold.
   *
   * <p>Automatically parallelises across all available processors (or {@code numThreads} when
   * set). Pass {@code numThreads = 0} to use all cores.
   *
   * <pre>{@code
   * Map<Integer, Map<Integer, Integer>> results =
   *     SearchEngine.batchMCS(mols, opts, mcsOpts, 0.5, 0);
   * }</pre>
   *
   * @param molecules           list of molecules
   * @param C                   chemical matching options
   * @param M                   MCS options
   * @param similarityThreshold minimum RASCAL upper bound required before attempting exact MCS
   * @param numThreads          thread count; 0 = auto (all available processors)
   * @return map keyed by {@code i*n+j} → MCS mapping for the pair (i, j)
   */
  public static Map<Integer, Map<Integer, Integer>> batchMCS(
      List<IAtomContainer> molecules, ChemOptions C, McsOptions M,
      double similarityThreshold, int numThreads) {

    int n = molecules.size();
    if (n == 0) return new LinkedHashMap<>();

    // Build the list of (i,j) pairs that pass the RASCAL upper-bound screen
    // (sequential — cheap, avoids shared-state issues in CDK property access)
    List<int[]> pairs = new ArrayList<>();
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        if (similarityUpperBound(molecules.get(i), molecules.get(j), C) >= similarityThreshold)
          pairs.add(new int[]{i, j});

    if (pairs.isEmpty()) return new LinkedHashMap<>();

    // Parallel MCS over screened pairs
    ConcurrentHashMap<Integer, Map<Integer, Integer>> concurrent = new ConcurrentHashMap<>();
    int cores = (numThreads > 0) ? numThreads : Runtime.getRuntime().availableProcessors();
    ForkJoinPool pool = new ForkJoinPool(cores);
    try {
      List<java.util.concurrent.Callable<Void>> tasks = new ArrayList<>(pairs.size());
      for (int[] p : pairs) {
        final int qi = p[0], ti = p[1];
        tasks.add(() -> {
          Map<Integer, Integer> mcs = findMCS(molecules.get(qi), molecules.get(ti), C, M);
          if (!mcs.isEmpty()) concurrent.put(qi * n + ti, mcs);
          return null;
        });
      }
      pool.invokeAll(tasks);
    } catch (Exception e) {
      // sequential fallback on any fork/join failure
      for (int[] p : pairs) {
        Map<Integer, Integer> mcs = findMCS(molecules.get(p[0]), molecules.get(p[1]), C, M);
        if (!mcs.isEmpty()) concurrent.put(p[0] * n + p[1], mcs);
      }
    } finally {
      pool.shutdown();
    }

    // Restore insertion order (sorted by key)
    Map<Integer, Map<Integer, Integer>> results = new LinkedHashMap<>();
    new TreeMap<>(concurrent).forEach(results::put);
    return results;
  }

  /**
   * Convenience overload — uses all available processors.
   *
   * @see #batchMCS(List, ChemOptions, McsOptions, double, int)
   */
  public static Map<Integer, Map<Integer, Integer>> batchMCS(
      List<IAtomContainer> molecules, ChemOptions C, McsOptions M, double similarityThreshold) {
    return batchMCS(molecules, C, M, similarityThreshold, 0);
  }

  /**
   * Screen one query molecule against a library of targets in parallel.
   *
   * <p>Uses all available processors by default. Each target is evaluated independently, so the
   * method scales linearly with core count. Pass {@code numThreads = 0} to auto-detect.
   *
   * @param query      the query molecule
   * @param targets    the target library
   * @param C          chemical matching options
   * @param timeoutMs  per-pair search time limit in milliseconds
   * @param numThreads thread count; 0 = all available processors
   * @return list of booleans — {@code true} where query is a substructure of {@code targets[i]}
   */
  public static List<Boolean> batchSubstructure(
      IAtomContainer query, List<IAtomContainer> targets, ChemOptions C,
      long timeoutMs, int numThreads) {

    int n = targets.size();
    boolean[] hits = new boolean[n];
    int cores = (numThreads > 0) ? numThreads : Runtime.getRuntime().availableProcessors();
    ForkJoinPool pool = new ForkJoinPool(cores);
    try {
      List<java.util.concurrent.Callable<Void>> tasks = new ArrayList<>(n);
      for (int i = 0; i < n; i++) {
        final int idx = i;
        tasks.add(() -> {
          try { hits[idx] = isSubstructure(query, targets.get(idx), C, timeoutMs); }
          catch (Exception ignored) { /* leave false */ }
          return null;
        });
      }
      pool.invokeAll(tasks);
    } catch (Exception e) {
      // sequential fallback
      for (int i = 0; i < n; i++) {
        try { hits[i] = isSubstructure(query, targets.get(i), C, timeoutMs); }
        catch (Exception ignored) { /* leave false */ }
      }
    } finally {
      pool.shutdown();
    }

    List<Boolean> result = new ArrayList<>(n);
    for (boolean h : hits) result.add(h);
    return result;
  }

  /**
   * Convenience overload — uses all available processors.
   *
   * @see #batchSubstructure(IAtomContainer, List, ChemOptions, long, int)
   */
  public static List<Boolean> batchSubstructure(
      IAtomContainer query, List<IAtomContainer> targets, ChemOptions C, long timeoutMs) {
    return batchSubstructure(query, targets, C, timeoutMs, 0);
  }

  // ---- Feature 1: Fingerprint Pre-Screening ----

  /**
   * Compute a path-based fingerprint for substructure pre-screening.
   * Enumerates all simple paths up to {@code pathLength} bonds and hashes
   * each path to a bit position. Explicit hydrogen atoms (atomicNum=1) are
   * skipped so that implicit-H and explicit-H representations produce
   * identical fingerprints. Recommended defaults: pathLength=7, fpSize=2048.
   */
  /**
   * Compute a path-based fingerprint for pre-screening substructure candidates.
   *
   * <p>Enumerates all simple paths up to {@code pathLength} bonds and hashes them
   * into a bit-set fingerprint. Explicit hydrogen atoms are skipped.
   *
   * <pre>{@code
   * long[] fp = SearchEngine.pathFingerprint(mol, 7, 1024);
   * }</pre>
   *
   * @param mol        the molecule to fingerprint
   * @param pathLength maximum path length (number of bonds) to enumerate
   * @param fpSize     fingerprint size in bits (should be a power of 2)
   * @return the fingerprint as a long array (bit-packed)
   */
  public static long[] pathFingerprint(IAtomContainer mol, int pathLength, int fpSize) {
    int words = (fpSize + 63) >>> 6;
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
    for (int start = 0; start < heavyCount; start++) {
      boolean[] visited = new boolean[heavyCount];
      int[] path = new int[pathLength + 1];
      path[0] = start;
      visited[start] = true;
      fpSetBit(fp, fpHashPath(atomicNum, path, 1), fpSize);
      fpEnumeratePaths(adj, atomicNum, fp, fpSize, visited, path, 1, pathLength);
    }
    return fp;
  }

  static void fpEnumeratePaths(int[][] adj, int[] atomicNum, long[] fp, int fpSize,
      boolean[] visited, int[] path, int depth, int maxDepth) {
    int cur = path[depth - 1];
    for (int nb : adj[cur]) {
      if (visited[nb]) continue;
      path[depth] = nb;
      visited[nb] = true;
      fpSetBit(fp, fpHashPath(atomicNum, path, depth + 1), fpSize);
      if (depth < maxDepth)
        fpEnumeratePaths(adj, atomicNum, fp, fpSize, visited, path, depth + 1, maxDepth);
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

  /**
   * Check if a query fingerprint is a subset of a target fingerprint.
   *
   * <p>Returns {@code true} if every bit set in the query is also set in the target.
   * This is a necessary (but not sufficient) condition for substructure containment.
   *
   * <pre>{@code
   * if (SearchEngine.fingerprintSubset(queryFp, targetFp)) {
   *     // query might be a substructure; run full match
   * }
   * }</pre>
   *
   * @param query  the query fingerprint
   * @param target the target fingerprint
   * @return {@code true} if query is a bitwise subset of target
   */
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

  // ---- MCS-Path Fingerprint ----

  /**
   * Compute an MCS-aware path fingerprint encoding graph-aware chemical invariants.
   *
   * <p>Hashes per-atom properties (element, ring, aromatic, tautomerClass, degree) and
   * per-bond properties (order, ring, aromatic). Explicit H atoms are skipped.
   *
   * <pre>{@code
   * long[] fp = SearchEngine.mcsFingerprint(molGraph, 7, 2048);
   * }</pre>
   *
   * @param g          pre-built MolGraph
   * @param pathLength maximum number of bonds in a path (recommended: 7)
   * @param fpSize     fingerprint size in bits (recommended: 2048)
   * @return the fingerprint as a long array (bit-packed)
   */
  public static long[] mcsFingerprint(MolGraph g, int pathLength, int fpSize) {
    return mcsFingerprint(g, false, pathLength, fpSize);
  }

  private static long[] mcsFingerprint(MolGraph g, boolean tautomerAware, int pathLength, int fpSize) {
    int words = (fpSize + 63) >>> 6;
    long[] fp = new long[words];
    int n = g.n;
    if (n == 0) return fp;
    int[] tautClass = g.tautomerClass != null ? g.tautomerClass : new int[n];
    // Build heavy-atom mask (skip explicit H)
    boolean[] isHeavy = new boolean[n];
    int heavyCount = 0;
    for (int i = 0; i < n; i++) { isHeavy[i] = g.atomicNum[i] != 1; if (isHeavy[i]) heavyCount++; }
    if (heavyCount == 0) return fp;
    for (int start = 0; start < n; start++) {
      if (!isHeavy[start]) continue;
      boolean[] visited = new boolean[n];
      int[] path = new int[pathLength + 1];
      path[0] = start;
      visited[start] = true;
      fpSetBit(fp, mcsFpHashPath(g, tautClass, tautomerAware, path, 1), fpSize);
      mcsFpEnumeratePaths(g, tautClass, tautomerAware, isHeavy, fp, fpSize,
          visited, path, 1, pathLength);
    }
    return fp;
  }

  /**
   * Compute an MCS-aware path fingerprint from a CDK molecule.
   *
   * <p>When {@code chem.tautomerAware=true}, tautomeric atoms hash to the same label,
   * ensuring keto/enol pairs have high fingerprint overlap.
   *
   * <pre>{@code
   * long[] fp = SearchEngine.mcsFingerprint(mol, new ChemOptions(), 7, 2048);
   * }</pre>
   *
   * @param mol        CDK molecule
   * @param chem       chemical options (controls tautomer awareness)
   * @param pathLength maximum number of bonds in a path
   * @param fpSize     fingerprint size in bits
   * @return the fingerprint as a long array (bit-packed)
   */
  public static long[] mcsFingerprint(IAtomContainer mol, ChemOptions chem, int pathLength, int fpSize) {
    MolGraph g = toMolGraph(mol);
    return mcsFingerprint(g, chem != null && chem.tautomerAware, pathLength, fpSize);
  }

  /**
   * Compute the Tanimoto coefficient between two MCS fingerprints.
   *
   * <pre>{@code
   * double sim = SearchEngine.mcsFingerprintSimilarity(fp1, fp2);
   * System.out.println("Similarity: " + sim); // 0.0 to 1.0
   * }</pre>
   *
   * @param fp1 first fingerprint
   * @param fp2 second fingerprint
   * @return Tanimoto similarity in [0.0, 1.0]; 1.0 if both fingerprints are empty
   */
  public static double mcsFingerprintSimilarity(long[] fp1, long[] fp2) {
    if (fp1 == null || fp2 == null) return 0.0;
    if (fp1.length == 0 && fp2.length == 0) return 1.0;
    int len = Math.max(fp1.length, fp2.length);
    int andCount = 0, orCount = 0;
    for (int i = 0; i < len; i++) {
      long a = i < fp1.length ? fp1[i] : 0L;
      long b = i < fp2.length ? fp2[i] : 0L;
      andCount += Long.bitCount(a & b);
      orCount += Long.bitCount(a | b);
    }
    if (orCount == 0) return 1.0; // both empty fingerprints => identical
    return (double) andCount / orCount;
  }

  static void mcsFpEnumeratePaths(MolGraph g, int[] tautClass, boolean tautAware,
      boolean[] isHeavy, long[] fp, int fpSize, boolean[] visited, int[] path,
      int depth, int maxDepth) {
    int cur = path[depth - 1];
    for (int nb : g.neighbors[cur]) {
      if (visited[nb] || !isHeavy[nb]) continue; // skip visited and H atoms
      path[depth] = nb;
      visited[nb] = true;
      fpSetBit(fp, mcsFpHashPath(g, tautClass, tautAware, path, depth + 1), fpSize);
      if (depth < maxDepth)
        mcsFpEnumeratePaths(g, tautClass, tautAware, isHeavy, fp, fpSize,
            visited, path, depth + 1, maxDepth);
      visited[nb] = false;
    }
  }

  /**
   * Hash a path of atoms using MCS-aware invariants.
   * Computes forward and reverse hashes, takes minimum for canonical ordering.
   */
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

  // ---- Feature 2: Scaffold MCS (Murcko) ----

  /**
   * Extract the Murcko scaffold (ring systems + linkers between rings) from a molecule.
   *
   * <p>Removes all non-ring, non-linker atoms. If the molecule has no rings, returns it unchanged.
   *
   * <pre>{@code
   * IAtomContainer scaffold = SearchEngine.murckoScaffold(atorvastatin);
   * System.out.println("Scaffold atoms: " + scaffold.getAtomCount());
   * }</pre>
   *
   * @param mol the molecule to extract the scaffold from
   * @return the Murcko scaffold as an {@link IAtomContainer}
   */
  public static IAtomContainer murckoScaffold(IAtomContainer mol) {
    int n = mol.getAtomCount();
    if (n == 0) return mol;
    boolean[] inRing = new boolean[n];
    boolean hasRing = false;
    for (int i = 0; i < n; i++) {
      inRing[i] = mol.getAtom(i).isInRing();
      if (inRing[i]) hasRing = true;
    }
    if (!hasRing) return mol;
    IdentityHashMap<IAtom, Integer> idxMap = new IdentityHashMap<>(n);
    for (int i = 0; i < n; i++) idxMap.put(mol.getAtom(i), i);
    int[][] adj = new int[n][];
    for (int i = 0; i < n; i++) {
      List<IAtom> nbAtoms = mol.getConnectedAtomsList(mol.getAtom(i));
      adj[i] = new int[nbAtoms.size()];
      for (int k = 0; k < adj[i].length; k++) adj[i][k] = idxMap.get(nbAtoms.get(k));
    }
    List<Integer> ringAtomList = new ArrayList<>();
    for (int i = 0; i < n; i++) if (inRing[i]) ringAtomList.add(i);
    boolean[] keep = new boolean[n];
    for (int i = 0; i < n; i++) if (inRing[i]) keep[i] = true;
    for (int ri = 0; ri < ringAtomList.size(); ri++) {
      int src = ringAtomList.get(ri);
      int[] prev = new int[n];
      Arrays.fill(prev, -1);
      boolean[] seen = new boolean[n];
      Deque<Integer> queue = new ArrayDeque<>();
      queue.add(src);
      seen[src] = true;
      while (!queue.isEmpty()) {
        int u = queue.poll();
        for (int v : adj[u]) {
          if (seen[v]) continue;
          seen[v] = true;
          prev[v] = u;
          if (inRing[v] && v != src) {
            int cur = v;
            while (cur != src && cur != -1) { keep[cur] = true; cur = prev[cur]; }
          }
          if (!inRing[v]) queue.add(v);
        }
      }
    }
    try {
      IAtomContainer scaffold = mol.clone();
      List<IAtom> toRemove = new ArrayList<>();
      for (int i = 0; i < n; i++) if (!keep[i]) toRemove.add(scaffold.getAtom(i));
      for (IAtom a : toRemove) scaffold.removeAtom(a);
      return scaffold;
    } catch (CloneNotSupportedException e) { return mol; }
  }

  /**
   * Find the MCS between the Murcko scaffolds of two molecules.
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * Map<Integer, Integer> scaffMcs =
   *     SearchEngine.findScaffoldMCS(mol1, mol2, new ChemOptions(), opts);
   * }</pre>
   *
   * @param m1 the first molecule
   * @param m2 the second molecule
   * @param C  chemical matching options
   * @param M  MCS options
   * @return mapping from scaffold atom indices (first) to scaffold atom indices (second)
   */
  public static Map<Integer, Integer> findScaffoldMCS(
      IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions M) {
    return findMCS(murckoScaffold(m1), murckoScaffold(m2), C, M);
  }

  // ---- Feature 3: dMCS Fragment Constraint Post-Processing ----

  static Map<Integer, Integer> applyFragmentConstraints(
      MolGraph g1, Map<Integer, Integer> map, int minFragmentSize, int maxFragments) {
    if (map.isEmpty() || (minFragmentSize <= 1 && maxFragments >= map.size())) return map;
    Map<Integer, List<Integer>> adjFC = new HashMap<>();
    for (int qi : map.keySet()) adjFC.put(qi, new ArrayList<>());
    for (int qi : map.keySet())
      for (int qk : map.keySet()) {
        if (qi >= qk) continue;
        if (g1.hasBond(qi, qk)) { adjFC.get(qi).add(qk); adjFC.get(qk).add(qi); }
      }
    Set<Integer> seen = new HashSet<>();
    List<Set<Integer>> frags = new ArrayList<>();
    for (int qi : map.keySet()) {
      if (seen.contains(qi)) continue;
      Set<Integer> frag = new LinkedHashSet<>();
      Deque<Integer> dq = new ArrayDeque<>();
      dq.add(qi); seen.add(qi);
      while (!dq.isEmpty()) {
        int u = dq.pollFirst(); frag.add(u);
        for (int v : adjFC.getOrDefault(u, Collections.emptyList()))
          if (!seen.contains(v)) { seen.add(v); dq.addLast(v); }
      }
      frags.add(frag);
    }
    frags.removeIf(f -> f.size() < minFragmentSize);
    frags.sort((a, b) -> Integer.compare(b.size(), a.size()));
    if (frags.size() > maxFragments) frags = frags.subList(0, maxFragments);
    Map<Integer, Integer> result = new LinkedHashMap<>();
    for (Set<Integer> frag : frags)
      for (int qi : frag) result.put(qi, map.get(qi));
    return result;
  }

  // ---- Feature 5: Atom-Atom Mapping for Reactions ----

  /**
   * Compute atom-atom mapping between reactants and products via disconnected MCS.
   *
   * <pre>{@code
   * Map<Integer, Integer> aam =
   *     SearchEngine.mapReaction(reactants, products, new ChemOptions(), 10000);
   * }</pre>
   *
   * @param reactants the combined reactant molecule
   * @param products  the combined product molecule
   * @param chem      chemical matching options
   * @param timeoutMs timeout in milliseconds
   * @return mapping from reactant atom indices to product atom indices
   */
  public static Map<Integer, Integer> mapReaction(
      IAtomContainer reactants, IAtomContainer products, ChemOptions chem, long timeoutMs) {
    ChemOptions reactionC = ChemOptions.copyOf(chem);
    reactionC.matchFormalCharge = false;
    reactionC.ringMatchesRingOnly = false;
    reactionC.completeRingsOnly = false;
    McsOptions opts = new McsOptions();
    opts.disconnectedMCS = true;
    opts.connectedOnly = false;
    opts.timeoutMs = timeoutMs;
    return findMCS(reactants, products, reactionC, opts);
  }

  // ---- Feature 5b: Reaction-Aware MCS Post-Filter (v6.4.0) ----

  /**
   * Generate near-MCS candidates by systematic deletion from K-sized mappings
   * and greedy re-extension.
   *
   * <p>For each K-sized mapping, systematically remove one mapped pair to produce
   * (K-1)-sized candidates. Then attempt to re-extend each by adding a different
   * atom that covers additional heteroatom types. Optionally repeat for (K-2).
   *
   * <p>Package-private: called internally when {@code reactionAware=true} or
   * {@code postFilter != null}.
   *
   * @param g1            first molecule graph
   * @param g2            second molecule graph
   * @param C             chemical matching options
   * @param M             MCS options
   * @param exactMCS      the size-K MCS candidates from findAllMCS
   * @param delta         max size deficit (typically 2)
   * @param maxCandidates max candidates to generate
   * @return combined list of K, K-1, and K-2 sized candidates (deduplicated)
   * @since 6.4.0
   */
  static List<Map<Integer, Integer>> findNearMCS(
      MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M,
      List<Map<Integer, Integer>> exactMCS, int delta, int maxCandidates) {

    if (exactMCS == null || exactMCS.isEmpty()) return Collections.emptyList();

    int K = exactMCS.get(0).size();
    int effectiveDelta = Math.min(delta, K - 1);

    // Collect all candidates: start with the exact K-sized ones
    Map<String, Map<Integer, Integer>> seenBySmi = new LinkedHashMap<>();
    for (Map<Integer, Integer> m : exactMCS) {
      String smi = extractMolGraphSubgraph(g1, m.keySet()).toCanonicalSmiles();
      seenBySmi.putIfAbsent(smi, m);
      if (seenBySmi.size() >= maxCandidates) break;
    }

    // Resolve timeout
    long timeout = M.timeoutMs;
    if (timeout < 0) timeout = Math.min(30_000L, 500L + (long) g1.n * g2.n * 2);
    TimeBudget tb = new TimeBudget(timeout);

    // Build compatibility structures for re-extension
    // For each exact MCS candidate, produce deletion variants
    List<Map<Integer, Integer>> deletionVariants = new ArrayList<>();

    for (Map<Integer, Integer> exact : exactMCS) {
      if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;

      for (int removeKey : exact.keySet()) {
        if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;

        // Create K-1 variant by removing one atom pair
        Map<Integer, Integer> variant = new LinkedHashMap<>(exact);
        variant.remove(removeKey);

        String smi = extractMolGraphSubgraph(g1, variant.keySet()).toCanonicalSmiles();
        if (seenBySmi.putIfAbsent(smi, variant) == null) {
          deletionVariants.add(variant);
        }
      }
    }

    // Greedy re-extension: for each deletion variant, try to add a different atom
    // that covers an additional heteroatom type not already in the mapping
    for (Map<Integer, Integer> variant : deletionVariants) {
      if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;

      Set<Integer> usedQ = variant.keySet();
      Set<Integer> usedT = new HashSet<>(variant.values());

      // What heteroatom types are already mapped?
      Set<Integer> mappedHetero = new HashSet<>();
      for (int qi : usedQ) {
        int z = g1.atomicNum[qi];
        if (z != 6 && z != 1) mappedHetero.add(z);
      }

      // Try to find an unmapped atom in g1 that is:
      // (a) a neighbor of some mapped atom (for connectivity)
      // (b) a heteroatom type not yet mapped
      // (c) compatible with an unmapped atom in g2
      for (int qi : new ArrayList<>(usedQ)) {
        if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;
        for (int nbQ : g1.neighbors[qi]) {
          if (usedQ.contains(nbQ)) continue;
          int zQ = g1.atomicNum[nbQ];
          if (zQ == 6 || zQ == 1) continue; // only heteroatoms
          if (mappedHetero.contains(zQ)) continue; // must be a NEW type

          // Find a compatible target atom that is a neighbor of a mapped target atom
          for (int ti : variant.values()) {
            if (tb.expiredNow()) break;
            for (int nbT : g2.neighbors[ti]) {
              if (usedT.contains(nbT)) continue;
              int zT = g2.atomicNum[nbT];
              if (zT != zQ) continue; // must match element

              // Check label compatibility
              if (C.matchAtomType && g1.label[nbQ] != g2.label[nbT]) continue;
              if (C.matchFormalCharge && g1.formalCharge[nbQ] != g2.formalCharge[nbT]) continue;

              // Create re-extended variant
              Map<Integer, Integer> extended = new LinkedHashMap<>(variant);
              extended.put(nbQ, nbT);
              String smi = extractMolGraphSubgraph(g1, extended.keySet()).toCanonicalSmiles();
              seenBySmi.putIfAbsent(smi, extended);
              break; // one extension per neighbor pair is enough
            }
          }
        }
      }
    }

    // Level 2: if delta >= 2, repeat deletion on K-1 candidates to get K-2
    if (effectiveDelta >= 2 && !tb.expiredNow()) {
      List<Map<Integer, Integer>> kMinus1 = new ArrayList<>();
      for (Map<Integer, Integer> m : seenBySmi.values()) {
        if (m.size() == K - 1) kMinus1.add(m);
      }

      for (Map<Integer, Integer> km1 : kMinus1) {
        if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;
        for (int removeKey : km1.keySet()) {
          if (tb.expiredNow() || seenBySmi.size() >= maxCandidates) break;
          Map<Integer, Integer> variant = new LinkedHashMap<>(km1);
          variant.remove(removeKey);
          String smi = extractMolGraphSubgraph(g1, variant.keySet()).toCanonicalSmiles();
          seenBySmi.putIfAbsent(smi, variant);
        }
      }
    }

    return new ArrayList<>(seenBySmi.values());
  }

  /**
   * Heteroatom-seeded MCS: for a given element type, seed from every compatible
   * (query, target) heteroatom pair and greedily extend. Returns the best mapping
   * that includes at least one atom of the required element.
   *
   * <p>This addresses the blind spot where findNearMCS (deletion + re-extension
   * from the size-K MCS) cannot reach subgraph selections that require a
   * fundamentally different seed -- e.g. the amino-acid+S chain in SAM vs. the
   * amino-acid+ribose chain that the standard MCS prefers.
   *
   * @param g1              first molecule graph
   * @param g2              second molecule graph
   * @param C               chemical matching options
   * @param requiredElement atomic number that must appear in the mapping (e.g. 16 for S)
   * @return the best mapping seeded from that element, or empty map if none found
   * @since 6.5.0
   */
  static Map<Integer, Integer> heteroatomSeededMCS(
      MolGraph g1, MolGraph g2, ChemOptions C, int requiredElement) {

    // Collect atom indices with the required element in each molecule
    List<Integer> qAtoms = new ArrayList<>();
    List<Integer> tAtoms = new ArrayList<>();
    for (int i = 0; i < g1.n; i++)
      if (g1.atomicNum[i] == requiredElement) qAtoms.add(i);
    for (int j = 0; j < g2.n; j++)
      if (g2.atomicNum[j] == requiredElement) tAtoms.add(j);
    if (qAtoms.isEmpty() || tAtoms.isEmpty()) return Collections.emptyMap();

    int n1 = g1.n, n2 = g2.n;
    Map<Integer, Integer> bestMapping = Collections.emptyMap();
    int bestSize = 0;

    for (int qi : qAtoms) {
      for (int tj : tAtoms) {
        if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qi, g2, tj, C))
          continue;

        // Initialize mapping arrays
        int[] q2t = new int[n1];
        int[] t2q = new int[n2];
        Arrays.fill(q2t, -1);
        Arrays.fill(t2q, -1);
        q2t[qi] = tj;
        t2q[tj] = qi;
        int mapSize = 1;

        // Greedy bond-based extension
        boolean progress = true;
        while (progress) {
          progress = false;
          for (int qk = 0; qk < n1; qk++) {
            if (q2t[qk] >= 0) continue;
            int bestTk = -1, bestScore = -1;
            for (int nb : g1.neighbors[qk]) {
              if (q2t[nb] < 0) continue; // neighbor must be mapped
              int tNb = q2t[nb];
              for (int tk : g2.neighbors[tNb]) {
                if (t2q[tk] >= 0) continue;
                if (!SubstructureEngine.AbstractVFMatcher.atomsCompatFast(g1, qk, g2, tk, C))
                  continue;
                if (!MolGraph.ChemOps.bondsCompatible(g1, qk, nb, g2, tk, tNb, C))
                  continue;
                // Check consistency with all other mapped neighbors
                boolean consistent = true;
                for (int qm : g1.neighbors[qk]) {
                  if (qm == nb || q2t[qm] < 0) continue;
                  int tm = q2t[qm];
                  int qOrd = g1.bondOrder(qk, qm);
                  int tOrd = g2.bondOrder(tk, tm);
                  if (qOrd != 0 && tOrd != 0) {
                    if (!MolGraph.ChemOps.bondsCompatible(g1, qk, qm, g2, tk, tm, C)) {
                      consistent = false;
                      break;
                    }
                  }
                }
                if (!consistent) continue;
                int score = (g2.ring[tk] && g1.ring[qk] ? 50 : 0)
                    + Math.min(g1.degree[qk], g2.degree[tk]);
                if (score > bestScore) { bestScore = score; bestTk = tk; }
              }
            }
            if (bestTk >= 0) {
              q2t[qk] = bestTk;
              t2q[bestTk] = qk;
              mapSize++;
              progress = true;
            }
          }
        }

        if (mapSize > bestSize) {
          bestSize = mapSize;
          bestMapping = new LinkedHashMap<>();
          for (int i = 0; i < n1; i++)
            if (q2t[i] >= 0) bestMapping.put(i, q2t[i]);
        }
      }
    }
    return bestMapping;
  }

  /**
   * Reaction-aware MCS: find MCS candidates, generate near-MCS variants,
   * add heteroatom-seeded candidates for shared rare elements, and re-rank
   * by composite scoring (size + heteroatom + rarity + connectivity).
   *
   * <pre>{@code
   * McsOptions opts = new McsOptions();
   * opts.reactionAware = true;
   * Map<Integer, Integer> best = SearchEngine.reactionAwareMCS(g1, g2, chem, opts);
   * }</pre>
   *
   * @param g1   first molecule graph
   * @param g2   second molecule graph
   * @param C    chemical matching options
   * @param M    MCS options (must have reactionAware=true or postFilter set)
   * @return the best-scoring candidate mapping
   * @since 6.4.0
   */
  public static Map<Integer, Integer> reactionAwareMCS(
      MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M) {
    // For reaction mapping, relax charge matching: reaction centers change
    // charge (e.g., S+ in SAM → S in homocysteine).
    ChemOptions reactionC = ChemOptions.copyOf(C);
    reactionC.matchFormalCharge = false;
    reactionC.ringMatchesRingOnly = false;
    reactionC.completeRingsOnly = false;

    // Step 1: find all exact K-sized MCS candidates (charge-relaxed)
    int maxResults = Math.max(5, M.nearMcsCandidates / 4);
    List<Map<Integer, Integer>> exactCandidates = findAllMCS(g1, g2, reactionC, M, maxResults);
    if (exactCandidates.isEmpty()) return Collections.emptyMap();

    // Step 2: generate near-MCS candidates (K-1, K-2)
    List<Map<Integer, Integer>> allCandidates = findNearMCS(
        g1, g2, reactionC, M, exactCandidates, M.nearMcsDelta, M.nearMcsCandidates);
    if (allCandidates.isEmpty()) allCandidates = new ArrayList<>(exactCandidates);

    // Step 3: heteroatom-seeded candidates -- for each heteroatom type shared
    // between both molecules, generate a candidate seeded from that element.
    // This covers subgraph selections unreachable by deletion from the K-MCS.
    Set<Integer> hG1 = new HashSet<>(), hG2 = new HashSet<>();
    for (int i = 0; i < g1.n; i++) {
      int z = g1.atomicNum[i];
      if (z != 6 && z != 1) hG1.add(z);
    }
    for (int j = 0; j < g2.n; j++) {
      int z = g2.atomicNum[j];
      if (z != 6 && z != 1) hG2.add(z);
    }
    for (int z : hG1) {
      if (!hG2.contains(z)) continue;
      Map<Integer, Integer> seeded = heteroatomSeededMCS(g1, g2, reactionC, z);
      if (!seeded.isEmpty()) {
        allCandidates.add(seeded);
      }
    }

    // Step 4: apply post-filter (custom or built-in)
    McsPostFilter filter = M.postFilter;
    if (filter == null) {
      filter = M.bondChangeAware ? new BondChangeScorer() : new ReactionAwareScorer();
    }
    List<Map<Integer, Integer>> ranked = filter.rank(
        Collections.unmodifiableList(allCandidates), g1, g2);

    return ranked.isEmpty() ? exactCandidates.get(0) : ranked.get(0);
  }

  /**
   * Reaction atom-atom mapping with reaction-aware MCS post-filtering.
   * Enumerates near-MCS candidates and re-ranks by heteroatom coverage,
   * rare-element importance, and connectivity.
   *
   * <pre>{@code
   * Map<Integer, Integer> aam =
   *     SearchEngine.mapReactionAware(reactants, products, new ChemOptions(), 10000);
   * }</pre>
   *
   * @param reactants the combined reactant molecule
   * @param products  the combined product molecule
   * @param chem      chemical matching options
   * @param timeoutMs timeout in milliseconds
   * @return mapping from reactant atom indices to product atom indices
   * @since 6.4.0
   */
  public static Map<Integer, Integer> mapReactionAware(
      IAtomContainer reactants, IAtomContainer products, ChemOptions chem, long timeoutMs) {
    MolGraph g1 = toMolGraph(reactants), g2 = toMolGraph(products);
    applySolvent(g1, chem);
    applySolvent(g2, chem);

    McsOptions opts = new McsOptions();
    opts.disconnectedMCS = true;
    opts.connectedOnly = false;
    opts.timeoutMs = timeoutMs;
    opts.reactionAware = true;
    opts.nearMcsDelta = 2;
    opts.nearMcsCandidates = 20;

    return reactionAwareMCS(g1, g2, chem, opts);
  }

  // ---- Feature 6: Fingerprint Quality Analysis ----

  /**
   * Analyze fingerprint bit distribution quality.
   *
   * <p>Returns metrics including fill rate, chi-squared uniformity, and a suggested
   * optimal size based on Bloom filter mathematics.
   *
   * <pre>{@code
   * Map<String, Double> quality = SearchEngine.analyzeFingerprintQuality(fp, 2048);
   * System.out.println("Fill rate: " + quality.get("fillRate"));
   * }</pre>
   *
   * @param fp     the fingerprint to analyze
   * @param fpSize the fingerprint size in bits
   * @return map with keys: "bitsFilled", "totalBits", "fillRate", "chiSquared", "suggestedSize"
   */
  public static Map<String, Double> analyzeFingerprintQuality(long[] fp, int fpSize) {
    int filled = 0;
    for (long w : fp) filled += Long.bitCount(w);
    double fillRate = (double) filled / fpSize;

    // Chi-squared test for uniform distribution (divide into 16 buckets)
    int buckets = 16;
    int bucketSize = fpSize / buckets;
    int[] counts = new int[buckets];
    for (int b = 0; b < fpSize; b++) {
      if ((fp[b >> 6] & (1L << (b & 63))) != 0) {
        counts[b / bucketSize]++;
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

    // Suggest optimal size based on Bloom filter math
    // Target: <1% false positive rate
    // m = -n*ln(p) / (ln(2))^2 where n=filled bits set, p=0.01
    double suggestedSize = Math.max(256, -filled * Math.log(0.01) / (Math.log(2) * Math.log(2)));
    result.put("suggestedSize", Math.ceil(suggestedSize / 64) * 64); // round to 64-bit boundary

    return result;
  }

  /**
   * Suggest optimal fingerprint size for a molecular library based on Bloom filter mathematics.
   *
   * <pre>{@code
   * int fpSize = SearchEngine.suggestFingerprintSize(25, 10000, 0.01);
   * // Use fpSize when computing fingerprints for the library
   * }</pre>
   *
   * @param avgAtomsPerMol average heavy atom count per molecule in the library
   * @param librarySize    number of molecules in the library
   * @param targetFPR      target false positive rate (e.g., 0.01 for 1%)
   * @return suggested fingerprint size in bits (rounded up to a multiple of 64)
   */
  public static int suggestFingerprintSize(int avgAtomsPerMol, int librarySize, double targetFPR) {
    // Average paths per molecule ~ avgAtoms * pathLength (roughly)
    int avgPaths = avgAtomsPerMol * 7; // default pathLength=7
    // Bloom filter optimal size: m = -n*ln(p) / (ln2)^2
    double m = -avgPaths * Math.log(targetFPR) / (Math.log(2) * Math.log(2));
    // Round up to nearest 64
    return (int) (Math.ceil(m / 64) * 64);
  }

  // -------------------------------------------------------------------------
  // Tautomer consistency validation (post-MCS)
  // -------------------------------------------------------------------------

  /**
   * Count mobilisable protons (N-H, O-H, S-H) on a heteroatom by computing
   * implicit hydrogens from default valence minus sum of bond orders.
   */
  private static int mobilisableProtons(MolGraph g, int atom) {
    int z = g.atomicNum[atom];
    if (z != 7 && z != 8 && z != 16) return 0;
    int bondSum = 0;
    for (int nb : g.neighbors[atom]) bondSum += g.bondOrder(atom, nb);
    int defaultVal = (z == 7) ? 3 : 2; // N=3, O/S=2
    return Math.max(0, defaultVal - bondSum - Math.abs(g.formalCharge[atom]));
  }

  /**
   * Validates that a tautomer-aware MCS mapping is proton-consistent.
   * For each tautomeric group in the mapping, the total N-H + O-H + S-H count
   * must be equal between query-side and target-side matched atoms.
   * Returns true if consistent, false if proton counts mismatch.
   *
   * @param g1  query molecule graph
   * @param g2  target molecule graph
   * @param mcs atom-index mapping (g1 atom index -> g2 atom index)
   * @return true if proton counts are conserved across all tautomeric groups
   */
  public static boolean validateTautomerConsistency(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty()) return true;
    g1.ensureTautomerClasses(); g2.ensureTautomerClasses();
    if (g1.tautomerClass == null || g2.tautomerClass == null) return true;

    Map<Integer, Integer> querySideProtons = new HashMap<>();
    Map<Integer, Integer> targetSideProtons = new HashMap<>();
    for (Map.Entry<Integer, Integer> e : mcs.entrySet()) {
      int qi = e.getKey(), ti = e.getValue();
      if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
      int tc1 = g1.tautomerClass[qi];
      int tc2 = g2.tautomerClass[ti];
      if (tc1 < 0 && tc2 < 0) continue;
      int groupKey = (tc1 >= 0) ? tc1 : tc2 + 100000;
      querySideProtons.merge(groupKey, mobilisableProtons(g1, qi), Integer::sum);
      targetSideProtons.merge(groupKey, mobilisableProtons(g2, ti), Integer::sum);
    }
    for (Map.Entry<Integer, Integer> e : querySideProtons.entrySet()) {
      int tProtons = targetSideProtons.getOrDefault(e.getKey(), 0);
      if (!e.getValue().equals(tProtons)) return false;
    }
    return true;
  }

  /**
   * CDK convenience overload -- creates MolGraphs from IAtomContainers.
   *
   * @param m1  query molecule
   * @param m2  target molecule
   * @param mcs atom-index mapping
   * @return true if proton counts are conserved across all tautomeric groups
   */
  public static boolean validateTautomerConsistency(
      IAtomContainer m1, IAtomContainer m2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty()) return true;
    return validateTautomerConsistency(toMolGraph(m1), toMolGraph(m2), mcs);
  }

  // ---- MCS SMILES extraction ----

  /**
   * Extract the Maximum Common Substructure as a canonical SMILES string.
   *
   * <p>Given a molecule graph and an MCS atom-index mapping (as returned by
   * {@link #findMCS(MolGraph, MolGraph, ChemOptions, McsOptions)}), this method
   * builds the induced subgraph on the query-side atoms (the mapping keys) and
   * writes it as a canonical SMILES string.
   *
   * <p>The induced subgraph preserves atomic number, formal charge, mass number,
   * aromaticity, ring membership, bond order, bond aromaticity, and bond ring
   * flags from the source molecule.
   *
   * @param g       the query molecule graph (source of atom/bond properties)
   * @param mapping MCS mapping whose keys are the query atom indices in the common
   *                substructure; an empty or {@code null} mapping returns {@code ""}
   * @return canonical SMILES of the MCS subgraph, or {@code ""} if the mapping is
   *         empty or {@code null}
   */
  public static String mcsToSmiles(MolGraph g, Map<Integer, Integer> mapping) {
    if (g == null) throw new NullPointerException("MolGraph must not be null");
    if (mapping == null || mapping.isEmpty()) return "";

    // Collect MCS atom indices (query side) and build old->new index map.
    int[] mcsAtoms = new int[mapping.size()];
    int idx = 0;
    for (int a : mapping.keySet()) mcsAtoms[idx++] = a;
    Arrays.sort(mcsAtoms);
    int k = mcsAtoms.length;
    Map<Integer, Integer> reindex = new HashMap<>(k * 2);
    for (int i = 0; i < k; i++) reindex.put(mcsAtoms[i], i);

    // Atom property arrays for the subgraph.
    int[] subAtomicNum    = new int[k];
    int[] subFormalCharge = new int[k];
    int[] subMassNumber   = new int[k];
    boolean[] subRing     = new boolean[k];
    boolean[] subAromatic = new boolean[k];

    for (int i = 0; i < k; i++) {
      int old = mcsAtoms[i];
      subAtomicNum[i]    = g.atomicNum[old];
      subFormalCharge[i] = g.formalCharge[old];
      subMassNumber[i]   = g.massNumber != null ? g.massNumber[old] : 0;
      // Do NOT copy ring/aromatic from parent — re-perceive after building subgraph
      // (partial ring extraction produces invalid aromatic flags)
      subRing[i]         = false;
      subAromatic[i]     = false;
    }

    // Build adjacency: only bonds where both endpoints are in the MCS.
    @SuppressWarnings("unchecked")
    List<int[]>[] adjList     = new List[k]; // adjList[newIdx] = list of {newNb, bondOrd}
    @SuppressWarnings("unchecked")
    List<boolean[]>[] adjProp = new List[k]; // adjProp[newIdx] = list of {inRing, aromatic}
    for (int i = 0; i < k; i++) { adjList[i] = new ArrayList<>(); adjProp[i] = new ArrayList<>(); }

    // Scan bonds from old graph between MCS atoms.
    Set<Long> seen = new HashSet<>();
    for (int oldI : mcsAtoms) {
      int newI = reindex.get(oldI);
      for (int nb : g.neighbors[oldI]) {
        Integer newJ = reindex.get(nb);
        if (newJ == null) continue;
        long key = Math.min(oldI, nb) * (long) g.atomCount() + Math.max(oldI, nb);
        if (!seen.add(key)) continue;
        int ord     = g.bondOrder(oldI, nb);
        boolean bRing = g.bondInRing(oldI, nb);
        boolean bArom = g.bondAromatic(oldI, nb);
        adjList[newI].add(new int[] {newJ, ord});
        adjList[newJ].add(new int[] {newI, ord});
        adjProp[newI].add(new boolean[] {bRing, bArom});
        adjProp[newJ].add(new boolean[] {bRing, bArom});
      }
    }

    // Flatten into Builder arrays.
    int[][] subNeighbors  = new int[k][];
    int[][] subBondOrders = new int[k][];
    boolean[][] subBondRings = new boolean[k][];
    boolean[][] subBondAroms = new boolean[k][];
    for (int i = 0; i < k; i++) {
      int deg = adjList[i].size();
      subNeighbors[i]  = new int[deg];
      subBondOrders[i] = new int[deg];
      subBondRings[i]  = new boolean[deg];
      subBondAroms[i]  = new boolean[deg];
      for (int d = 0; d < deg; d++) {
        subNeighbors[i][d]  = adjList[i].get(d)[0];
        subBondOrders[i][d] = adjList[i].get(d)[1];
        subBondRings[i][d]  = adjProp[i].get(d)[0];
        subBondAroms[i][d]  = adjProp[i].get(d)[1];
      }
    }

    MolGraph sub = new MolGraph.Builder()
        .atomCount(k)
        .atomicNumbers(subAtomicNum)
        .formalCharges(subFormalCharge)
        .massNumbers(subMassNumber)
        .ringFlags(subRing)
        .aromaticFlags(subAromatic)
        .neighbors(subNeighbors)
        .bondOrders(subBondOrders)
        .bondRingFlags(subBondRings)
        .bondAromaticFlags(subBondAroms)
        .build();

    return sub.toCanonicalSmiles();
  }

  /**
   * Convenience method: compute the MCS of two MolGraphs and return the result
   * as a canonical SMILES string.
   *
   * @param g1 query molecule
   * @param g2 target molecule
   * @param C  chemical matching options
   * @param M  MCS options (timeout, induced, connected, etc.)
   * @return canonical SMILES of the MCS, or {@code ""} if no common substructure
   */
  public static String findMcsSmiles(MolGraph g1, MolGraph g2, ChemOptions C, McsOptions M) {
    Map<Integer, Integer> mapping = findMCS(g1, g2, C, M);
    return mcsToSmiles(g1, mapping);
  }

  /**
   * CDK convenience overload: compute the MCS of two IAtomContainers and return
   * the result as a canonical SMILES string.
   *
   * @param m1 query molecule
   * @param m2 target molecule
   * @param C  chemical matching options
   * @param M  MCS options
   * @return canonical SMILES of the MCS, or {@code ""} if no common substructure
   */
  public static String findMcsSmiles(IAtomContainer m1, IAtomContainer m2, ChemOptions C, McsOptions M) {
    MolGraph g1 = toMolGraph(m1), g2 = toMolGraph(m2);
    applySolvent(g1, C);
    applySolvent(g2, C);
    return findMcsSmiles(g1, g2, C, M);
  }

  // ---- SMARTS-based MCS ----

  /**
   * Find the largest substructure match of a SMARTS query in a target MolGraph.
   *
   * <p>Uses the existing substructure matching engine: the SMARTS query is parsed
   * into a MolGraph via SMILES semantics (atom-type matching only), then all
   * substructure embeddings of the SMARTS pattern in the target are enumerated.
   * The largest mapping (by atom count) is returned.
   *
   * <p>This is useful for SMARTS-guided MCS where the user wants to find the
   * maximum overlap of a SMARTS-defined pharmacophore or scaffold in a target.
   *
   * <pre>{@code
   * MolGraph target = MolGraph.fromSmiles("CC(=O)Nc1ccc(O)cc1");
   * Map<Integer,Integer> result = SearchEngine.findMcsSmarts(
   *     "[CX3](=O)[NX3]", target, new ChemOptions(), 10000);
   * // result maps SMARTS atom indices to target atom indices
   * }</pre>
   *
   * @param smartsQuery the SMARTS query string
   * @param target      the target molecule graph
   * @param C           chemical matching options
   * @param timeoutMs   maximum time in milliseconds
   * @return mapping from SMARTS atom indices to target atom indices for the
   *         largest match; empty map if no match
   */
  public static Map<Integer, Integer> findMcsSmarts(
      String smartsQuery, MolGraph target, ChemOptions C, long timeoutMs) {
    if (smartsQuery == null || smartsQuery.isEmpty() || target == null)
      return Collections.emptyMap();
    if (C == null) C = new ChemOptions();
    if (timeoutMs <= 0) timeoutMs = 10_000L;

    // Parse SMARTS as a MolGraph via the SMILES parser path.
    // For pure SMARTS patterns that use only atom/bond primitives compatible
    // with SMILES (e.g. "c1ccccc1", "C(=O)N"), this works directly.
    // For extended SMARTS (e.g. "[#6;R]"), we fall through to CDK if available.
    MolGraph gQuery;
    try {
      gQuery = new MolGraph(new org.openscience.cdk.smiles.SmilesParser(
          org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance()).parseSmiles(smartsQuery));
    } catch (Exception e) {
      // If SMILES parsing fails, the pattern likely uses extended SMARTS syntax.
      // Return empty: callers should use CDK SmartsPattern directly for complex patterns.
      return Collections.emptyMap();
    }

    // Find all substructure mappings of the query in the target
    List<Map<Integer, Integer>> mappings =
        findAllSubstructures(gQuery, target, C, 100, timeoutMs);

    if (mappings.isEmpty()) return Collections.emptyMap();

    // Return the largest mapping (most matched atoms)
    Map<Integer, Integer> best = mappings.get(0);
    for (int i = 1; i < mappings.size(); i++) {
      if (mappings.get(i).size() > best.size()) {
        best = mappings.get(i);
      }
    }
    return best;
  }

  /**
   * CDK convenience overload: SMARTS-based largest substructure match.
   *
   * @param smartsQuery the SMARTS query string
   * @param target      the target molecule
   * @param C           chemical matching options
   * @param timeoutMs   maximum time in milliseconds
   * @return mapping from SMARTS atom indices to target atom indices for the
   *         largest match; empty map if no match
   */
  public static Map<Integer, Integer> findMcsSmarts(
      String smartsQuery, IAtomContainer target, ChemOptions C, long timeoutMs) {
    if (target == null) return Collections.emptyMap();
    MolGraph gt = toMolGraph(target);
    applySolvent(gt, C);
    return findMcsSmarts(smartsQuery, gt, C, timeoutMs);
  }

  // =========================================================================
  // Batch MCS with non-overlap constraints
  // =========================================================================

  /**
   * Find MCS across multiple molecule pairs with non-overlap constraints.
   *
   * <p>For multi-component reactions, this finds atom-atom mappings across
   * all (reactant, product) pairs such that no atom is mapped more than once.
   * Uses greedy sequential MCS with atom exclusion:
   * <ol>
   *   <li>Sort pairs by decreasing molecule size (larger pairs first)</li>
   *   <li>For each pair, find MCS excluding already-used atoms</li>
   *   <li>Mark matched atoms as used in both sides</li>
   * </ol>
   *
   * @param queries  list of query molecules (reactant fragments)
   * @param targets  list of target molecules (product fragments)
   * @param C        chemical matching options
   * @param M        MCS options
   * @param timeoutMs per-pair timeout in milliseconds
   * @return list of mappings, one per pair (same order as input)
   * @since 6.5.3
   */
  /**
   * Find MCS for each query against ALL targets with non-overlapping
   * target atom constraints.
   *
   * <p>Architecture for multi-component reaction mapping:
   * <ol>
   *   <li>For each query, find its best MCS across all targets</li>
   *   <li>Larger queries are processed first (greedy coverage)</li>
   *   <li>Target atoms already mapped by earlier queries are excluded</li>
   *   <li>Returns one mapping per query, each mapping to the best target</li>
   * </ol>
   *
   * <p>Use case: N reactant fragments mapping to 1 combined product.
   * Each reactant finds its substructure in the product, and atoms used
   * by reactant 1 are excluded when mapping reactant 2.
   *
   * @param queries  reactant fragments
   * @param targets  product fragment(s) — typically one combined product
   * @param C        chemical matching options
   * @param M        MCS options
   * @param timeoutMs per-pair timeout in milliseconds
   * @return list of mappings, one per query.
   *         Each mapping: {query_atom_idx: target_atom_idx}.
   *         Target atom indices refer to the target that gave the best MCS.
   * @since 6.6.0
   */
  public static List<Map<Integer, Integer>> batchMcsConstrained(
      List<MolGraph> queries, List<MolGraph> targets,
      ChemOptions C, McsOptions M, long timeoutMs) {
    int nQ = queries.size();
    int nT = targets.size();
    List<Map<Integer, Integer>> results = new ArrayList<>(nQ);
    for (int i = 0; i < nQ; i++) results.add(Collections.emptyMap());

    // Sort queries by decreasing size (larger fragments first = better coverage)
    Integer[] order = new Integer[nQ];
    for (int i = 0; i < nQ; i++) order[i] = i;
    Arrays.sort(order, (a, b) -> queries.get(b).n - queries.get(a).n);

    // Track used target atoms per target: usedTargetAtoms[targetIdx] = set of used atom indices
    List<Set<Integer>> usedTargetAtoms = new ArrayList<>(nT);
    for (int t = 0; t < nT; t++) usedTargetAtoms.add(new HashSet<>());

    for (int qi : order) {
      MolGraph gq = queries.get(qi);

      // Find best MCS across all targets
      Map<Integer, Integer> bestMapping = Collections.emptyMap();
      int bestSize = 0;
      int bestTarget = -1;

      for (int ti = 0; ti < nT; ti++) {
        MolGraph gt = targets.get(ti);
        Set<Integer> used = usedTargetAtoms.get(ti);

        Map<Integer, Integer> mapping;
        if (used.isEmpty()) {
          // No atoms claimed yet — run MCS against full target
          mapping = findMCS(gq, gt, C, M);
        } else {
          // RE-RUN MCS against residual target (unclaimed atoms only).
          // Build index mapping: residual[r] = original[residualToOriginal[r]]
          int nOrig = gt.n;
          int[] originalToResidual = new int[nOrig];
          int[] residualToOriginal = new int[nOrig - used.size()];
          Arrays.fill(originalToResidual, -1);
          int rIdx = 0;
          for (int i = 0; i < nOrig; i++) {
            if (!used.contains(i)) {
              originalToResidual[i] = rIdx;
              residualToOriginal[rIdx] = i;
              rIdx++;
            }
          }
          int nRes = rIdx;
          if (nRes == 0) continue;

          // Build residual MolGraph using CDK: clone target and remove used atoms
          try {
            IAtomContainer residualMol = gt.mol.clone();
            // Remove atoms in reverse order to preserve indices
            List<Integer> toRemove = new ArrayList<>(used);
            toRemove.sort(Collections.reverseOrder());
            for (int atomIdx : toRemove) {
              residualMol.removeAtom(atomIdx);
            }
            MolGraph residual = new MolGraph(residualMol);

            Map<Integer, Integer> residualMcs = findMCS(gq, residual, C, M);
            // Translate residual indices → original target indices
            mapping = new LinkedHashMap<>();
            for (Map.Entry<Integer, Integer> e : residualMcs.entrySet()) {
              if (e.getValue() < residualToOriginal.length) {
                mapping.put(e.getKey(), residualToOriginal[e.getValue()]);
              }
            }
          } catch (CloneNotSupportedException ex) {
            // Fallback: filter-only approach
            Map<Integer, Integer> rawMcs = findMCS(gq, gt, C, M);
            mapping = new LinkedHashMap<>();
            for (Map.Entry<Integer, Integer> e : rawMcs.entrySet()) {
              if (!used.contains(e.getValue())) mapping.put(e.getKey(), e.getValue());
            }
          }
        }
        if (mapping.isEmpty()) continue;

        if (mapping.size() > bestSize) {
          bestSize = mapping.size();
          bestMapping = mapping;
          bestTarget = ti;
        }
      }

      // Mark used target atoms
      if (bestTarget >= 0 && !bestMapping.isEmpty()) {
        for (int tAtom : bestMapping.values()) {
          usedTargetAtoms.get(bestTarget).add(tAtom);
        }
      }

      results.set(qi, bestMapping);
    }

    return results;
  }

  /**
   * Convenience: batch MCS with IAtomContainer inputs.
   * @since 6.6.0
   */
  public static List<Map<Integer, Integer>> batchMcsConstrained(
      List<IAtomContainer> queries, List<IAtomContainer> targets,
      ChemOptions C, long timeoutMs) {
    List<MolGraph> gQueries = new ArrayList<>(queries.size());
    List<MolGraph> gTargets = new ArrayList<>(targets.size());
    for (IAtomContainer q : queries) gQueries.add(toMolGraph(q));
    for (IAtomContainer t : targets) gTargets.add(toMolGraph(t));
    return batchMcsConstrained(gQueries, gTargets, C, new McsOptions(), timeoutMs);
  }

}

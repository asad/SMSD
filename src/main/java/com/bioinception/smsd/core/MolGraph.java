/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.Random;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;

/**
 * Immutable, CDK-free molecular graph representation for fast substructure
 * and MCS computation. Can be built from a CDK {@link IAtomContainer} or
 * programmatically via {@link Builder}.
 */
public final class MolGraph {
  final IAtomContainer mol;
  final int n;
  public final int[] atomicNum, formalCharge, label, massNumber;
  int[] morganRank, orbit, canonicalLabel;
  int[][] autGenerators;
  boolean autGeneratorsTruncated;
  long canonicalHash;
  private volatile boolean canonicalComputed = false;
  private volatile boolean tautomerClassesComputed = false;
  private volatile boolean ringCountsComputed = false;
  final boolean[] ring, aromatic;
  final int[] degree;
  public final int[] ringCount;
  final int[][] neighbors;
  private volatile BitSet[] adj;          // lazy — built on first access via getAdj()
  long[][] adjLong;
  final int words;
  private volatile int[][] cachedNLF1, cachedNLF2, cachedNLF3;
  private volatile int[][] cachedNeighborsByDegDesc;
  private volatile int[] cachedPharmacophoreFeatures;
  private static final int SPARSE_THRESHOLD = 200;
  final int[][] bondOrdMatrix;
  final boolean[][] bondRingMatrix, bondAromMatrix;
  final HashMap<Long, int[]> sparseBondProps;
  final int[] tetraChirality;
  final int[][] dbStereoConf;
  public int[]   tautomerClass;
  /**
   * Per-atom pKa-informed relevance weight for tautomeric matches, in [0,1].
   * 1.0 = atom not in any tautomeric system (no confidence penalty).
   * Values below 1.0 express how likely this atom's current resonance form
   * is at physiological pH 7.4, based on empirical pKa data (Sitzmann 2010).
   * Used by SearchEngine.computeTautomerConfidence() after MCS.
   */
  public float[] tautomerWeight;
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
  ChemOptions.Solvent solvent = ChemOptions.Solvent.AQUEOUS;
  double pH = 7.4;
  int[][] sssrRings;
  int[][] rcbRings; // Relevant Cycle Basis (union of all MCBs)
  private volatile int[] ringSystemId;
  private volatile long[] ringSysSignature;
  private volatile int ringSysCount;
  private volatile boolean ringSysComputed = false;

  // ---------------------------------------------------------------------------
  // Empirical pKa-derived tautomer relevance weights at pH 7.4
  // References: Sitzmann et al. JCIM 2010; Karabina et al. 2014; ChEMBL pKa db
  // ---------------------------------------------------------------------------
  public static final float TW_KETO_ENOL         = 0.98f; // pKa ~20  — keto overwhelmingly dominant
  public static final float TW_AMIDE_IMIDIC      = 0.97f; // pKa ~15-17 — amide form dominant
  public static final float TW_LACTAM_LACTIM     = 0.97f; // cyclic amide — lactam dominant
  public static final float TW_UREA              = 0.96f; // pKa ~14 — urea form dominant
  public static final float TW_PYRIDINONE        = 0.95f; // pKa ~1  — lactam dominant
  public static final float TW_THIOAMIDE         = 0.94f; // pKa ~13 — thioamide dominant
  public static final float TW_THIONE_THIOL      = 0.94f; // pKa ~10-11 — thione dominant
  public static final float TW_ENAMINONE         = 0.85f; // vinylogous amide — amine form dominant
  public static final float TW_HYDROXAMIC        = 0.85f; // pKa ~8  — neutral hydroxamic dominant
  public static final float TW_PHENOL_QUINONE    = 0.88f; // pKa ~10 — phenol dominant
  public static final float TW_HYDROXYPYRIMIDINE = 0.90f; // lactam tautomer dominant in N-heterocycles
  public static final float TW_PURINE_NH         = 0.78f; // N9-H dominant in purines (adenine, guanine)
  public static final float TW_DIKETONE_ENOL     = 0.72f; // 1,3-diketone enolisation: pKa ~8-11, pH-sensitive, ~28% enol at pH 7.4
  public static final float TW_NITROSO_OXIME     = 0.95f; // oxime strongly dominant over nitroso (Karabina 2014)
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

  static final int HASH_PRIME = 1000003;

  public boolean hasBond(int i, int j) { return bondOrder(i, j) != 0; }

  public int atomCount() { return n; }

  /**
   * Implicit hydrogen count for atom {@code idx}.
   *
   * <p>Mirrors the C++ {@code MolGraph::hydrogenCount[]} field byte-for-byte.
   * Used by {@link FingerprintEngine#classifyPharmacophore} to distinguish
   * pyrrole-type aromatic N (has H, not an H-bond acceptor) from pyridine-type
   * aromatic N (no H, is an H-bond acceptor), and by the canonical SMILES
   * writer to emit the stereo H count inside bracket atoms.
   *
   * <p>For CDK-backed graphs this reads
   * {@code IAtom.getImplicitHydrogenCount()} directly, which is always
   * accurate regardless of how bond orders are encoded. For graphs built via
   * {@link Builder} (no CDK backing), it falls back to {@link
   * FingerprintEngine#implicitH} using a Kekulé-equivalent bond-order sum
   * (aromatic bonds counted as 1) so that pyrrole-like detection still works.
   *
   * <p>Why this method exists: SMSD internally encodes aromatic bonds as bond
   * order 4 (for bond-order matching semantics), so callers cannot compute
   * implicit H from {@code bondOrder()} summation alone — the aromatic-4
   * encoding would blow past the valence table and return 0 for every
   * aromatic atom.
   */
  public int hydrogenCount(int idx) {
    if (mol != null) {
      Integer h = mol.getAtom(idx).getImplicitHydrogenCount();
      return h != null ? h : 0;
    }
    // Builder path: no CDK backing. Reconstruct a Kekulé-equivalent bond
    // order sum by counting aromatic bonds as order 1 (single) for implicit-H
    // arithmetic. This matches the openSMILES model used by the C++ engine
    // when it populates its own hydrogenCount[] field.
    int bondOrdSum = 0;
    for (int nb : neighbors[idx]) {
      int bo = bondOrder(idx, nb);
      bondOrdSum += (bo == 4) ? 1 : bo;
    }
    return FingerprintEngine.implicitH(atomicNum[idx], bondOrdSum, formalCharge[idx]);
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

  public int[] getOrbit() { ensureCanonical(); return orbit; }

  public int[][] computeRings() {
    if (sssrRings != null) return sssrRings;
    if (n == 0) { sssrRings = new int[0][]; return sssrRings; }

    int edgeCount = 0;
    for (int i = 0; i < n; i++) edgeCount += neighbors[i].length;
    edgeCount /= 2;

    boolean[] visited = new boolean[n];
    int components = 0;
    for (int i = 0; i < n; i++) {
      if (visited[i]) continue;
      components++;
      Deque<Integer> bfsQ = new ArrayDeque<>();
      bfsQ.add(i);
      visited[i] = true;
      while (!bfsQ.isEmpty()) {
        int u = bfsQ.poll();
        for (int v : neighbors[u]) if (!visited[v]) { visited[v] = true; bfsQ.add(v); }
      }
    }
    int cycleRank = edgeCount - n + components;
    if (cycleRank <= 0) { sssrRings = new int[0][]; return sssrRings; }

    // Horton candidate generation: for every vertex v and every edge (eU, eV),
    // find shortest path from v to eU and from v to eV (excluding edge eU-eV),
    // form cycle from both paths + the edge. This is O(V*E*V) but guaranteed
    // to find all cycles needed for both SSSR and relevant cycles.

    // Build edge index
    Map<Long, Integer> edgeIdx = new HashMap<>();
    List<int[]> edgeList = new ArrayList<>();
    for (int i = 0; i < n; i++)
      for (int j : neighbors[i])
        if (j > i) {
          edgeIdx.put(bondKey(i, j), edgeList.size());
          edgeList.add(new int[]{i, j});
        }
    int totalEdges = edgeList.size();

    // Helper: BFS shortest path from start to target, avoiding edge (avoidU, avoidV)
    // Returns atom path, or empty list if unreachable
    List<boolean[]> candidateMasks = new ArrayList<>();
    List<int[]> candidateAtoms = new ArrayList<>();
    List<Integer> candidateSizes = new ArrayList<>();
    Set<String> seenMasks = new HashSet<>();

    for (int v = 0; v < n; v++) {
      for (int[] edge : edgeList) {
        int eU = edge[0], eV = edge[1];

        // Find shortest path from v to eU, avoiding edge (eU, eV)
        List<Integer> path1 = bfsPath(v, eU, eU, eV);
        if (path1.isEmpty()) continue;

        // Find shortest path from v to eV, avoiding edge (eU, eV)
        List<Integer> path2 = bfsPath(v, eV, eU, eV);
        if (path2.isEmpty()) continue;

        // Check paths are vertex-disjoint (only share vertex v)
        Set<Integer> p1Interior = new HashSet<>();
        for (int i = 1; i < path1.size(); i++) p1Interior.add(path1.get(i));
        boolean disjoint = true;
        for (int i = 1; i < path2.size(); i++) {
          if (p1Interior.contains(path2.get(i))) { disjoint = false; break; }
        }
        if (!disjoint) continue;

        // Build edge mask for this cycle
        boolean[] mask = new boolean[totalEdges];
        int edgeId = edgeIdx.get(bondKey(eU, eV));
        mask[edgeId] = true;

        boolean valid = true;
        for (int i = 0; i < path1.size() - 1 && valid; i++) {
          Integer idx = edgeIdx.get(bondKey(Math.min(path1.get(i), path1.get(i+1)),
                                            Math.max(path1.get(i), path1.get(i+1))));
          if (idx == null) valid = false; else mask[idx] = true;
        }
        for (int i = 0; i < path2.size() - 1 && valid; i++) {
          Integer idx = edgeIdx.get(bondKey(Math.min(path2.get(i), path2.get(i+1)),
                                            Math.max(path2.get(i), path2.get(i+1))));
          if (idx == null) valid = false; else mask[idx] = true;
        }
        if (!valid) continue;

        // Deduplicate by edge mask
        String maskStr = Arrays.toString(mask);
        if (!seenMasks.add(maskStr)) continue;

        int ringSize = path1.size() + path2.size() - 1;
        if (ringSize < 3 || ringSize > 20) continue;

        candidateMasks.add(mask);
        candidateSizes.add(ringSize);

        // Build atom list: path1(v→eU) + eV + reverse(path2 interior, eV→...→v excluded)
        // path1 = [v, ..., eU], path2 = [v, ..., eV]
        // Cycle: v → ... → eU → eV → ... → v (backwards through path2)
        List<Integer> atoms = new ArrayList<>(path1); // [v, ..., eU]
        // Append path2 in reverse, skipping the first element (v) since path1 starts at v
        for (int i = path2.size() - 1; i >= 1; i--) atoms.add(path2.get(i));
        candidateAtoms.add(atoms.stream().mapToInt(Integer::intValue).toArray());
      }
    }

    // Sort candidates by size (smallest first)
    Integer[] order = new Integer[candidateMasks.size()];
    for (int i = 0; i < order.length; i++) order[i] = i;
    Arrays.sort(order, Comparator.comparingInt(candidateSizes::get));

    // Interleaved SSSR + RCB pipeline (Vismara's rule, 2-phase reduction):
    // Phase 1: Check relevance (independence from STRICTLY SHORTER cycles only)
    // Phase 2: Check SSSR inclusion (independence from ALL basis cycles)
    List<int[]> sssr = new ArrayList<>();
    List<int[]> rcb = new ArrayList<>();
    List<boolean[]> rcbMasks = new ArrayList<>();
    List<boolean[]> basisMatrix = new ArrayList<>();
    List<Integer> basisSizes = new ArrayList<>();

    // O(E) leading-bit lookup for row-echelon reduction
    int[] basisByLeadingBit = new int[totalEdges];
    Arrays.fill(basisByLeadingBit, -1);

    for (int idx : order) {
      boolean[] mask = candidateMasks.get(idx);
      int cSize = candidateSizes.get(idx);

      // Phase 1: Check relevance — reduce only against STRICTLY SHORTER basis cycles
      boolean[] phase1Mask = mask.clone();
      for (int b = 0; b < totalEdges; b++) {
        if (phase1Mask[b]) {
          int basisIdx = basisByLeadingBit[b];
          if (basisIdx != -1 && basisSizes.get(basisIdx) < cSize) {
            for (int k = b; k < totalEdges; k++) {
              phase1Mask[k] ^= basisMatrix.get(basisIdx)[k];
            }
          }
        }
      }

      boolean relevant = false;
      for (boolean bit : phase1Mask) if (bit) { relevant = true; break; }

      if (relevant) {
        // Deduplicate for RCB by edge mask
        boolean dup = false;
        for (boolean[] existingMask : rcbMasks) {
          if (Arrays.equals(mask, existingMask)) { dup = true; break; }
        }
        if (!dup) {
          rcb.add(candidateAtoms.get(idx));
          rcbMasks.add(mask);
        }

        // Phase 2: Check SSSR inclusion — reduce against ALL basis cycles
        boolean[] phase2Mask = phase1Mask.clone();
        for (int b = 0; b < totalEdges; b++) {
          if (phase2Mask[b]) {
            int basisIdx = basisByLeadingBit[b];
            if (basisIdx != -1) {
              for (int k = b; k < totalEdges; k++) {
                phase2Mask[k] ^= basisMatrix.get(basisIdx)[k];
              }
            }
          }
        }

        boolean independent = false;
        int newLeadingBit = -1;
        for (int b = 0; b < totalEdges; b++) {
          if (phase2Mask[b]) { independent = true; newLeadingBit = b; break; }
        }

        if (independent && sssr.size() < cycleRank) {
          sssr.add(candidateAtoms.get(idx));
          basisMatrix.add(phase2Mask); // Add the FULLY REDUCED mask
          basisSizes.add(cSize);
          basisByLeadingBit[newLeadingBit] = basisMatrix.size() - 1;
        }
      }
    }

    sssrRings = sssr.toArray(new int[0][]);
    rcbRings = rcb.toArray(new int[0][]);
    return sssrRings;
  }

  public int numRings() { return computeRings().length; }

  // ---- Ring-system decomposition (lazy, thread-safe) ----

  private void ensureRingSystems() {
    if (ringSysComputed) return;
    synchronized (this) {
      if (ringSysComputed) return;
      int[] uf = new int[n];
      for (int i = 0; i < n; i++) uf[i] = i;
      // Unite ring atoms connected by bonds where both endpoints are ring atoms
      for (int i = 0; i < n; i++) {
        if (!ring[i]) continue;
        for (int j : neighbors[i]) {
          if (j > i && ring[j]) ufUnion(uf, i, j);
        }
      }
      // Assign contiguous IDs
      int[] idMap = new int[n];
      Arrays.fill(idMap, -1);
      int count = 0;
      int[] sysId = new int[n];
      Arrays.fill(sysId, -1);
      for (int i = 0; i < n; i++) {
        if (!ring[i]) continue;
        int root = ufFind(uf, i);
        if (idMap[root] < 0) idMap[root] = count++;
        sysId[i] = idMap[root];
      }
      // Compute FNV-1a hash per ring system
      long FNV_OFFSET = -3750763034362895579L;
      long FNV_PRIME  = 1099511628211L;
      long[] sig = new long[count];
      // Gather per-system data: size, aromCount, atomicNums, bondOrders
      int[] sysSize = new int[count];
      int[] sysArom = new int[count];
      @SuppressWarnings("unchecked")
      List<Integer>[] sysAtomicNums = new List[count];
      @SuppressWarnings("unchecked")
      List<Integer>[] sysBondOrds = new List[count];
      for (int s = 0; s < count; s++) {
        sysAtomicNums[s] = new ArrayList<>();
        sysBondOrds[s] = new ArrayList<>();
      }
      for (int i = 0; i < n; i++) {
        if (sysId[i] < 0) continue;
        int s = sysId[i];
        sysSize[s]++;
        if (aromatic[i]) sysArom[s]++;
        sysAtomicNums[s].add(atomicNum[i]);
        for (int j : neighbors[i]) {
          if (j > i && ring[j] && sysId[j] == s) {
            sysBondOrds[s].add(bondOrder(i, j));
          }
        }
      }
      for (int s = 0; s < count; s++) {
        long h = FNV_OFFSET;
        h ^= sysSize[s]; h *= FNV_PRIME;
        h ^= sysArom[s]; h *= FNV_PRIME;
        Collections.sort(sysAtomicNums[s]);
        for (int v : sysAtomicNums[s]) { h ^= v; h *= FNV_PRIME; }
        Collections.sort(sysBondOrds[s]);
        for (int v : sysBondOrds[s]) { h ^= v; h *= FNV_PRIME; }
        sig[s] = h;
      }
      ringSystemId = sysId;
      ringSysSignature = sig;
      ringSysCount = count;
      ringSysComputed = true;
    }
  }

  public int ringSystemOf(int atom) { ensureRingSystems(); return ringSystemId[atom]; }
  public long ringSystemSig(int rsId) { ensureRingSystems(); return rsId >= 0 && rsId < ringSysCount ? ringSysSignature[rsId] : 0; }
  public int ringSystemCount() { ensureRingSystems(); return ringSysCount; }

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
    computeRings(); // ensures rcbRings is populated
    if (rcbRings != null) return rcbRings;
    return new int[0][];
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
    ensureCanonical();
    int[][] rc = computeRelevantCycles();
    if (rc.length == 0) return new int[0][][];

    // Build edge index + O(1) inverse lookup
    Map<Long, Integer> edgeIndex = new HashMap<>();
    for (int i = 0; i < n; i++) {
      for (int j : neighbors[i]) {
        if (j > i) edgeIndex.computeIfAbsent(bondKey(i, j), k -> edgeIndex.size());
      }
    }
    int numEdges = edgeIndex.size();
    int longWords = (numEdges + 63) >>> 6;

    // O(1) inverse lookup: edge index → (u, v) pair
    int[][] indexToEdge = new int[numEdges][2];
    for (Map.Entry<Long, Integer> e : edgeIndex.entrySet()) {
      long k = e.getKey();
      int idx = e.getValue();
      indexToEdge[idx] = new int[] { (int)(k >>> 32), (int)(k & 0xFFFFFFFFL) };
    }

    // Build edge-incidence vectors for each relevant cycle
    long[][] vectors = new long[rc.length][longWords];
    for (int ci = 0; ci < rc.length; ci++) {
      int[] c = rc[ci];
      for (int i = 0; i < c.length; i++) {
        int a = c[i], b = c[(i + 1) % c.length];
        long key = bondKey(Math.min(a, b), Math.max(a, b));
        Integer idx = edgeIndex.get(key);
        if (idx != null) vectors[ci][idx >>> 6] ^= 1L << (idx & 63);
      }
    }

    int nrc = rc.length;
    int[] family = new int[nrc];
    for (int i = 0; i < nrc; i++) family[i] = i;

    // Generate orbit-based signatures for each cycle (lexicographically minimal)
    String[] sigs = new String[nrc];
    for (int i = 0; i < nrc; i++) sigs[i] = getCycleOrbitSignature(rc[i]);

    // Group cycles into URFs using two criteria (in order):
    //   1. Orbit-signature match: only when every atom in both cycles shares the
    //      same orbit value (vertex-transitive ring systems, e.g. cubane where
    //      all atoms have orbit=0). Without this guard, fused systems like
    //      naphthalene are incorrectly merged because both 6-rings happen to
    //      produce the same orbit sequence after rotation/reflection.
    //   2. Kolodzik single-edge exchange: symmetric difference == 2 edges
    //      sharing a common vertex.
    for (int i = 0; i < nrc; i++) {
      for (int j = i + 1; j < nrc; j++) {
        if (rc[i].length != rc[j].length) continue;

        // 1. Orbit-signature match — guarded: only apply when all atoms in
        //    both cycles belong to the same orbit (vertex-transitive).
        if (sigs[i].equals(sigs[j])) {
          int orb0 = orbit[rc[i][0]];
          boolean allSame = true;
          for (int a : rc[i]) { if (orbit[a] != orb0) { allSame = false; break; } }
          if (allSame) for (int a : rc[j]) { if (orbit[a] != orb0) { allSame = false; break; } }
          if (allSame) {
            union(family, i, j);
            continue;
          }
        }

        // 2. Kolodzik single edge exchange (handles non-symmetric but interchangeable rings)
        long[] xor = new long[longWords];
        int popcount = 0;
        for (int w = 0; w < longWords; w++) {
          xor[w] = vectors[i][w] ^ vectors[j][w];
          popcount += Long.bitCount(xor[w]);
        }
        if (popcount == 2) {
          int[] diffEdges = new int[2];
          int di = 0;
          for (int w = 0; w < longWords && di < 2; w++) {
            long bits = xor[w];
            while (bits != 0 && di < 2) {
              int bit = Long.numberOfTrailingZeros(bits);
              diffEdges[di++] = (w << 6) + bit;
              bits &= bits - 1;
            }
          }
          int[] e1 = indexToEdge[diffEdges[0]]; // O(1) instead of O(E)
          int[] e2 = indexToEdge[diffEdges[1]];
          if (e1 != null && e2 != null && sharesVertex(e1, e2)) {
            union(family, i, j);
          }
        }
      }
    }

    // Group by family
    Map<Integer, List<Integer>> groups = new HashMap<>();
    for (int i = 0; i < nrc; i++) {
      groups.computeIfAbsent(find(family, i), k -> new ArrayList<>()).add(i);
    }

    int[][][] result = new int[groups.size()][][];
    int gi = 0;
    for (List<Integer> members : groups.values()) {
      result[gi] = new int[members.size()][];
      for (int m = 0; m < members.size(); m++) {
        result[gi][m] = rc[members.get(m)];
      }
      gi++;
    }
    return result;
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
    int[][] raw = g.computeRings();
    List<List<Integer>> rings = new ArrayList<>(raw.length);
    for (int[] r : raw) {
      List<Integer> ring = new ArrayList<>(r.length);
      for (int a : r) ring.add(a);
      rings.add(ring);
    }
    rings.sort(Comparator.comparingInt(List::size));
    return rings;
  }

  /**
   * Return SSSR rings optimized for 2D coordinate generation:
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
    List<List<Integer>> rings = computeSSSR(g);
    if (rings.size() <= 1) return rings;

    int nr = rings.size();

    // Build edge sets per ring
    List<Set<Long>> edgeSets = new ArrayList<>(nr);
    for (List<Integer> r : rings) {
      Set<Long> es = new HashSet<>();
      for (int j = 0; j < r.size(); j++) {
        int a = r.get(j), b = r.get((j + 1) % r.size());
        es.add(bondKey(Math.min(a, b), Math.max(a, b)));
      }
      edgeSets.add(es);
    }

    // Build ring adjacency
    List<List<Integer>> adj = new ArrayList<>(nr);
    for (int i = 0; i < nr; i++) adj.add(new ArrayList<>());
    for (int i = 0; i < nr; i++) {
      for (int j = i + 1; j < nr; j++) {
        boolean shared = false;
        for (Long e : edgeSets.get(i)) {
          if (edgeSets.get(j).contains(e)) { shared = true; break; }
        }
        if (shared) { adj.get(i).add(j); adj.get(j).add(i); }
      }
    }

    // Connected components
    int[] comp = new int[nr];
    Arrays.fill(comp, -1);
    int nComp = 0;
    for (int i = 0; i < nr; i++) {
      if (comp[i] >= 0) continue;
      int c = nComp++;
      comp[i] = c;
      Deque<Integer> bfs = new ArrayDeque<>();
      bfs.add(i);
      while (!bfs.isEmpty()) {
        int u = bfs.poll();
        for (int v : adj.get(u)) {
          if (comp[v] < 0) { comp[v] = c; bfs.add(v); }
        }
      }
    }

    // Group by component, compute system size (unique atom count)
    List<List<Integer>> systems = new ArrayList<>(nComp);
    for (int i = 0; i < nComp; i++) systems.add(new ArrayList<>());
    for (int i = 0; i < nr; i++) systems.get(comp[i]).add(i);

    int[] sysSize = new int[nComp];
    for (int c = 0; c < nComp; c++) {
      Set<Integer> atoms = new HashSet<>();
      for (int ri : systems.get(c))
        atoms.addAll(rings.get(ri));
      sysSize[c] = atoms.size();
    }

    // Sort systems: largest first
    Integer[] sysOrder = new Integer[nComp];
    for (int i = 0; i < nComp; i++) sysOrder[i] = i;
    Arrays.sort(sysOrder, (a, b) -> Integer.compare(sysSize[b], sysSize[a]));

    List<List<Integer>> result = new ArrayList<>(nr);
    for (int ci : sysOrder) {
      List<Integer> sysRings = systems.get(ci);
      if (sysRings.size() == 1) {
        result.add(rings.get(sysRings.get(0)));
        continue;
      }

      sysRings.sort(Comparator.comparingInt(a -> rings.get(a).size()));

      boolean[] visited = new boolean[nr];
      Deque<Integer> bfs = new ArrayDeque<>();
      bfs.add(sysRings.get(0));
      visited[sysRings.get(0)] = true;
      while (!bfs.isEmpty()) {
        int u = bfs.poll();
        result.add(rings.get(u));
        List<Integer> nbrs = new ArrayList<>();
        for (int v : adj.get(u)) {
          if (!visited[v] && comp[v] == ci) nbrs.add(v);
        }
        nbrs.sort(Comparator.comparingInt(a -> rings.get(a).size()));
        for (int v : nbrs) {
          if (!visited[v]) { visited[v] = true; bfs.add(v); }
        }
      }
    }
    return result;
  }

  /**
   * 2D point for layout coordinates.
   */
  public static final class Point2D {
    public double x, y;
    public Point2D(double x, double y) { this.x = x; this.y = y; }
  }

  private static Point2D[] normalizeCoords(MolGraph g, Point2D[] coords, double spacing) {
    Point2D[] out = new Point2D[g.n];
    int copy = Math.min(coords == null ? 0 : coords.length, g.n);
    for (int i = 0; i < copy; i++) {
      Point2D p = coords[i];
      out[i] = p != null ? new Point2D(p.x, p.y) : new Point2D(i * spacing, 0.0);
    }
    for (int i = copy; i < g.n; i++) out[i] = new Point2D(i * spacing, 0.0);
    return out;
  }

  private static void copyCoords(Point2D[] src, Point2D[] dst) {
    int copy = Math.min(src.length, dst.length);
    for (int i = 0; i < copy; i++) {
      if (dst[i] == null) dst[i] = new Point2D(src[i].x, src[i].y);
      else {
        dst[i].x = src[i].x;
        dst[i].y = src[i].y;
      }
    }
  }

  private static boolean isDegenerateLayout(Point2D[] coords) {
    if (coords.length < 3) return true;
    double minX = coords[0].x, maxX = coords[0].x;
    double minY = coords[0].y, maxY = coords[0].y;
    for (Point2D p : coords) {
      minX = Math.min(minX, p.x);
      maxX = Math.max(maxX, p.x);
      minY = Math.min(minY, p.y);
      maxY = Math.max(maxY, p.y);
    }
    return (maxX - minX < 1e-6) || (maxY - minY < 1e-6);
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
   * @return number of bond crossings remaining after optimization
   * @since 6.5.1
   */
  public static int reduceCrossings(MolGraph g, Point2D[] coords, int maxIter) {
    Point2D[] work = normalizeCoords(g, coords, 1.5);
    if (g.n < 4) return 0;

    int[][] bonds = buildBondList(g);
    List<Set<Integer>> ringSystems = findRingSystemSets(g);
    if (ringSystems.isEmpty()) {
      int result = countCrossingsImpl(bonds, work);
      copyCoords(work, coords);
      return result;
    }

    int nSys = ringSystems.size();
    int current = countCrossingsImpl(bonds, work);
    if (current == 0) {
      copyCoords(work, coords);
      return 0;
    }

    double T0 = 2.0, Tmin = 0.01, alpha = 0.995;
    Random rng = new Random(42);
    double T = T0;

    // ---- Phase 1: system-level flipping ----
    int phase1Iters = maxIter / 2;
    for (int iter = 0; iter < phase1Iters && current > 0; iter++) {
      int si = rng.nextInt(nSys);
      Set<Integer> sys = ringSystems.get(si);
      if (sys.size() < 3) continue;

      double cx = 0, cy = 0;
      for (int a : sys) { cx += work[a].x; cy += work[a].y; }
      cx /= sys.size(); cy /= sys.size();

      Map<Integer, double[]> saved = new HashMap<>();
      for (int a : sys) saved.put(a, new double[]{work[a].x, work[a].y});

      int axis = rng.nextInt(2);
      for (int a : sys) {
        if (axis == 0) work[a].x = 2.0 * cx - work[a].x;
        else           work[a].y = 2.0 * cy - work[a].y;
      }

      int newCr = countCrossingsImpl(bonds, work);
      int delta = newCr - current;
      boolean accept = (delta <= 0);
      if (!accept && T > Tmin) {
        accept = rng.nextDouble() < Math.exp(-(double) delta / T);
      }
      if (accept) { current = newCr; }
      else { for (Map.Entry<Integer, double[]> e : saved.entrySet()) {
        work[e.getKey()].x = e.getValue()[0];
        work[e.getKey()].y = e.getValue()[1];
      }}
      T *= alpha;
      if (T < Tmin) T = Tmin;
    }

    if (current == 0) return 0;

    // ---- Phase 2: individual ring flipping within multi-ring systems ----
    int[][] sssr = g.computeRings();
    if (sssr.length < 2) return current; // need ≥2 rings for individual flipping

    // Build flippable rings: for each ring in a multi-ring system,
    // exclusive = atoms only in this ring, shared = atoms in ≥2 rings.
    // Shared atoms act as fixed pivots during the flip.
    List<int[]> exclusiveList = new ArrayList<>();
    List<int[]> sharedList = new ArrayList<>();

    // Group rings by system membership (which system does each SSSR ring belong to?)
    // and count per-system ring counts
    Map<Integer, List<Integer>> sysToRings = new HashMap<>();
    for (int ri = 0; ri < sssr.length; ri++) {
      for (int si = 0; si < nSys; si++) {
        Set<Integer> sys = ringSystems.get(si);
        boolean belongs = true;
        for (int a : sssr[ri]) {
          if (!sys.contains(a)) { belongs = false; break; }
        }
        if (belongs) {
          sysToRings.computeIfAbsent(si, k -> new ArrayList<>()).add(ri);
          break;
        }
      }
    }

    for (Map.Entry<Integer, List<Integer>> entry : sysToRings.entrySet()) {
      List<Integer> ringIndices = entry.getValue();
      if (ringIndices.size() < 2) continue; // single ring — handled in phase 1

      // Count how many rings each atom belongs to within this system
      Map<Integer, Integer> atomRingCount = new HashMap<>();
      for (int ri : ringIndices) {
        for (int a : sssr[ri]) atomRingCount.merge(a, 1, Integer::sum);
      }

      for (int ri : ringIndices) {
        List<Integer> excl = new ArrayList<>(), shar = new ArrayList<>();
        for (int a : sssr[ri]) {
          if (atomRingCount.getOrDefault(a, 0) >= 2) shar.add(a);
          else excl.add(a);
        }
        if (!excl.isEmpty() && !shar.isEmpty()) {
          exclusiveList.add(excl.stream().mapToInt(Integer::intValue).toArray());
          sharedList.add(shar.stream().mapToInt(Integer::intValue).toArray());
        }
      }
    }

    if (!exclusiveList.isEmpty()) {
      int nFlip = exclusiveList.size();
      T = T0; // reset temperature for phase 2
      int phase2Iters = maxIter - phase1Iters;

      for (int iter = 0; iter < phase2Iters && current > 0; iter++) {
        int fi = rng.nextInt(nFlip);
        int[] excl = exclusiveList.get(fi);
        int[] shar = sharedList.get(fi);

        // Pivot = centroid of shared (fusion) atoms
        double cx = 0, cy = 0;
        for (int a : shar) { cx += work[a].x; cy += work[a].y; }
        cx /= shar.length; cy /= shar.length;

        // Save exclusive atom coordinates
        double[][] saved = new double[excl.length][2];
        for (int i = 0; i < excl.length; i++) {
          saved[i][0] = work[excl[i]].x;
          saved[i][1] = work[excl[i]].y;
        }

        // Mirror exclusive atoms about the shared-atom centroid
        int axis = rng.nextInt(2);
        for (int a : excl) {
          if (axis == 0) work[a].x = 2.0 * cx - work[a].x;
          else           work[a].y = 2.0 * cy - work[a].y;
        }

        int newCr = countCrossingsImpl(bonds, work);
        int delta = newCr - current;
        boolean accept = (delta <= 0);
        if (!accept && T > Tmin) {
          accept = rng.nextDouble() < Math.exp(-(double) delta / T);
        }
        if (accept) { current = newCr; }
        else {
          for (int i = 0; i < excl.length; i++) {
            work[excl[i]].x = saved[i][0];
            work[excl[i]].y = saved[i][1];
          }
        }
        T *= alpha;
        if (T < Tmin) T = Tmin;
      }
    }

    copyCoords(work, coords);
    return current;
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
   * @since 6.3.3
   */
  public static double forceDirectedLayout(MolGraph g, Point2D[] coords,
                                            int maxIter, double targetBondLength) {
    int n = g.n;
    if (n < 2) return 0.0;
    Point2D[] work = normalizeCoords(g, coords, targetBondLength);
    if (isDegenerateLayout(work)) {
      Point2D[] templ = matchTemplate(g, targetBondLength);
      if (templ != null && templ.length == n) work = templ;
    }

    final double kAttract  = 1.0;
    final double kRepel    = 0.5;
    final double kCrossing = 2.0;
    double step = 0.1;
    final double cooling = 0.98;
    final double convergence = 0.01;

    // Adjacency for O(1) bonded check
    boolean[][] bonded = new boolean[n][n];
    for (int i = 0; i < n; i++)
      for (int j : g.neighbors[i])
        bonded[i][j] = true;

    // Bond list
    int[][] bonds = buildBondList(g);
    int nb = bonds.length;

    for (int iter = 0; iter < maxIter; iter++) {
      double[] fx = new double[n], fy = new double[n];

      // 1. Attractive forces (bonded pairs)
      for (int[] bond : bonds) {
        int a = bond[0], b = bond[1];
        double dx = work[b].x - work[a].x;
        double dy = work[b].y - work[a].y;
        double dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 1e-10) dist = 1e-10;
        double force = kAttract * (dist - targetBondLength);
        double ux = dx / dist, uy = dy / dist;
        fx[a] += force * ux; fy[a] += force * uy;
        fx[b] -= force * ux; fy[b] -= force * uy;
      }

      // 2. Repulsive forces (non-bonded pairs)
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          if (bonded[i][j]) continue;
          double dx = work[j].x - work[i].x;
          double dy = work[j].y - work[i].y;
          double distSq = dx * dx + dy * dy;
          if (distSq < 1e-6) distSq = 1e-6;
          double force = -kRepel / distSq;
          double dist = Math.sqrt(distSq);
          double ux = dx / dist, uy = dy / dist;
          fx[i] += force * ux; fy[i] += force * uy;
          fx[j] -= force * ux; fy[j] -= force * uy;
        }
      }

      // 3. Crossing penalty forces
      for (int bi = 0; bi < nb; bi++) {
        for (int bj = bi + 1; bj < nb; bj++) {
          int a1 = bonds[bi][0], a2 = bonds[bi][1];
          int b1 = bonds[bj][0], b2 = bonds[bj][1];
          if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
          if (!segmentsCross(work[a1], work[a2], work[b1], work[b2])) continue;
          double m1x = (work[a1].x + work[a2].x) * 0.5;
          double m1y = (work[a1].y + work[a2].y) * 0.5;
          double m2x = (work[b1].x + work[b2].x) * 0.5;
          double m2y = (work[b1].y + work[b2].y) * 0.5;
          double dx = m2x - m1x, dy = m2y - m1y;
          double dist = Math.sqrt(dx * dx + dy * dy);
          if (dist < 1e-8) { dx = 0.1; dy = 0.1; dist = Math.sqrt(0.02); }
          double force = -kCrossing / (dist + 0.1);
          double ux = dx / dist, uy = dy / dist;
          fx[a1] += 0.5 * force * ux; fy[a1] += 0.5 * force * uy;
          fx[a2] += 0.5 * force * ux; fy[a2] += 0.5 * force * uy;
          fx[b1] -= 0.5 * force * ux; fy[b1] -= 0.5 * force * uy;
          fx[b2] -= 0.5 * force * ux; fy[b2] -= 0.5 * force * uy;
        }
      }

      // Update positions
      double maxDisp = 0.0;
      for (int i = 0; i < n; i++) {
        double dx = step * fx[i], dy = step * fy[i];
        double disp = Math.sqrt(dx * dx + dy * dy);
        if (disp > maxDisp) maxDisp = disp;
        work[i].x += dx;
        work[i].y += dy;
      }
      step *= cooling;
      if (maxDisp < convergence) break;
    }

    // Compute final stress
    double stress = 0.0;
    for (int[] bond : bonds) {
      double dx = work[bond[1]].x - work[bond[0]].x;
      double dy = work[bond[1]].y - work[bond[0]].y;
      double dist = Math.sqrt(dx * dx + dy * dy);
      double delta = dist - targetBondLength;
      stress += delta * delta;
    }
    copyCoords(work, coords);
    return stress;
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
   * @since 6.4.0
   */
  public static double stressMajorisation(MolGraph g, Point2D[] coords,
                                           int maxIter, double targetBondLength) {
    int n = g.n;
    if (n < 2) return 0.0;
    Point2D[] seedCoords = normalizeCoords(g, coords, targetBondLength);
    if (isDegenerateLayout(seedCoords)) {
      Point2D[] templ = matchTemplate(g, targetBondLength);
      if (templ != null && templ.length == n) seedCoords = templ;
    }

    // BFS all-pairs shortest path
    int[][] graphDist = new int[n][n];
    for (int src = 0; src < n; src++) {
      Arrays.fill(graphDist[src], n + 1);
      graphDist[src][src] = 0;
      Deque<Integer> q = new ArrayDeque<>();
      q.add(src);
      while (!q.isEmpty()) {
        int u = q.poll();
        for (int v : g.neighbors[u]) {
          if (graphDist[src][v] > graphDist[src][u] + 1) {
            graphDist[src][v] = graphDist[src][u] + 1;
            q.add(v);
          }
        }
      }
    }

    // Target distances and weights
    double[][] D = new double[n][n];
    double[][] W = new double[n][n];
    double[] Wsum = new double[n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        D[i][j] = graphDist[i][j] * targetBondLength;
        if (i != j && D[i][j] > 1e-10) {
          W[i][j] = 1.0 / (D[i][j] * D[i][j]);
        }
        Wsum[i] += W[i][j];
      }
    }

    final double eps = 1e-6;
    Random rng = new Random(42);
    double bestStress = Double.POSITIVE_INFINITY;
    Point2D[] bestCoords = seedCoords;
    int nInit = 3;
    for (int init = 0; init < Math.max(1, nInit); init++) {
      Point2D[] runCoords = normalizeCoords(g, seedCoords, targetBondLength);
      if (init > 0) {
        for (int i = 0; i < n; i++) {
          runCoords[i].x += (rng.nextDouble() - 0.5) * targetBondLength;
          runCoords[i].y += (rng.nextDouble() - 0.5) * targetBondLength;
        }
      }

      double prevStress = 0.0;
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          double dx = runCoords[i].x - runCoords[j].x;
          double dy = runCoords[i].y - runCoords[j].y;
          double dij = Math.sqrt(dx * dx + dy * dy);
          double delta = dij - D[i][j];
          prevStress += W[i][j] * delta * delta;
        }
      }

      for (int iter = 0; iter < maxIter; iter++) {
        Point2D[] nc = new Point2D[n];
        for (int i = 0; i < n; i++) {
          double sx = 0, sy = 0;
          for (int j = 0; j < n; j++) {
            if (i == j || W[i][j] < 1e-15) continue;
            double dx = runCoords[i].x - runCoords[j].x;
            double dy = runCoords[i].y - runCoords[j].y;
            double dij = Math.sqrt(dx * dx + dy * dy);
            if (dij < 1e-10) dij = 1e-10;
            double ratio = D[i][j] / dij;
            sx += W[i][j] * (runCoords[j].x + ratio * (runCoords[i].x - runCoords[j].x));
            sy += W[i][j] * (runCoords[j].y + ratio * (runCoords[i].y - runCoords[j].y));
          }
          nc[i] = Wsum[i] > 1e-15
              ? new Point2D(sx / Wsum[i], sy / Wsum[i])
              : new Point2D(runCoords[i].x, runCoords[i].y);
        }
        runCoords = nc;

        double stress = 0.0;
        for (int i = 0; i < n; i++) {
          for (int j = i + 1; j < n; j++) {
            double dx = runCoords[i].x - runCoords[j].x;
            double dy = runCoords[i].y - runCoords[j].y;
            double dij = Math.sqrt(dx * dx + dy * dy);
            double delta = dij - D[i][j];
            stress += W[i][j] * delta * delta;
          }
        }
        if (prevStress > 1e-10 && (prevStress - stress) / prevStress < eps) {
          prevStress = stress;
          break;
        }
        prevStress = stress;
      }

      if (prevStress < bestStress) {
        bestStress = prevStress;
        bestCoords = runCoords;
      }
    }
    copyCoords(bestCoords, coords);
    return bestStress;
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
   * @since 6.3.3
   */
  public static Point2D[] matchTemplate(MolGraph g, double targetBondLength) {
    int[][] rings = g.computeRings();
    if (rings.length == 0 && g.n != 6) return null;

    int[] sizes = new int[rings.length];
    for (int i = 0; i < rings.length; i++) sizes[i] = rings[i].length;
    Arrays.sort(sizes);

    // Try each template
    for (TemplateEntry t : SCAFFOLD_TEMPLATES) {
      if (t.atomCount() != g.n) continue;
      int[] tSizes = Arrays.copyOf(t.ringSizes(), t.ringSizes().length);
      Arrays.sort(tSizes);
      if (!Arrays.equals(sizes, tSizes)) continue;
      // Match found
      Point2D[] result = new Point2D[t.atomCount()];
      for (int i = 0; i < t.atomCount(); i++)
        result[i] = new Point2D(t.coords()[i][0] * targetBondLength,
                                t.coords()[i][1] * targetBondLength);
      return result;
    }
    return null;
  }

  /** Convenience overload: default targetBondLength = 1.5. */
  public static Point2D[] matchTemplate(MolGraph g) {
    return matchTemplate(g, 1.5);
  }

  // Scaffold template storage
  private record TemplateEntry(String name, int atomCount, int[] ringSizes, double[][] coords) {}

  private static double[][] regularPolygonCoords(int n) {
    double r = 0.5 / Math.sin(Math.PI / n);
    double[][] coords = new double[n][2];
    for (int i = 0; i < n; i++) {
      double angle = 2.0 * Math.PI * i / n - Math.PI / 2.0;
      coords[i][0] = r * Math.cos(angle);
      coords[i][1] = r * Math.sin(angle);
    }
    return coords;
  }

  private static final TemplateEntry[] SCAFFOLD_TEMPLATES;
  static {
    double[][] hex = regularPolygonCoords(6);
    SCAFFOLD_TEMPLATES = new TemplateEntry[] {
      // Benzene
      new TemplateEntry("Benzene", 6, new int[]{6}, hex),
      // Cyclohexane
      new TemplateEntry("Cyclohexane", 6, new int[]{6}, hex),
      // Piperidine
      new TemplateEntry("Piperidine", 6, new int[]{6}, hex),
      // Naphthalene (10 atoms, two 6-rings fused)
      new TemplateEntry("Naphthalene", 10, new int[]{6, 6},
          buildFused66Coords()),
      // Quinoline (10 atoms, two 6-rings fused)
      new TemplateEntry("Quinoline", 10, new int[]{6, 6},
          buildFused66Coords()),
      // Biphenyl (12 atoms, two 6-rings single-bonded)
      new TemplateEntry("Biphenyl", 12, new int[]{6, 6},
          buildBiphenylCoords()),
      // Indole (9 atoms, fused 5-6)
      new TemplateEntry("Indole", 9, new int[]{5, 6},
          buildFused56Coords()),
      // Purine (9 atoms, fused 5-6)
      new TemplateEntry("Purine", 9, new int[]{5, 6},
          buildFused56Coords()),
    };
  }

  private static double[][] buildFused66Coords() {
    double[][] hex = regularPolygonCoords(6);
    double[][] coords = new double[10][2];
    System.arraycopy(hex, 0, coords, 0, 6);
    // Fuse second ring sharing edge (1,2)
    double mx = (hex[1][0] + hex[2][0]) * 0.5;
    double my = (hex[1][1] + hex[2][1]) * 0.5;
    double edx = hex[2][0] - hex[1][0], edy = hex[2][1] - hex[1][1];
    double nx = -edy, ny = edx;
    // Normalise
    double nlen = Math.sqrt(nx * nx + ny * ny);
    nx /= nlen; ny /= nlen;
    // Centre of ring A
    double cax = 0, cay = 0;
    for (double[] c : hex) { cax += c[0]; cay += c[1]; }
    cax /= 6; cay /= 6;
    if (nx * (cax - mx) + ny * (cay - my) > 0) { nx = -nx; ny = -ny; }
    double rB = 0.5 / Math.sin(Math.PI / 6.0);
    double cd = rB * Math.cos(Math.PI / 6.0);
    double cx = mx + cd * nx, cy = my + cd * ny;
    double theta = Math.atan2(edy, edx);
    for (int k = 0; k < 4; k++) {
      double angle = theta + Math.PI + (2.0 * Math.PI * (k + 1)) / 6.0;
      coords[6 + k][0] = cx + rB * Math.cos(angle);
      coords[6 + k][1] = cy + rB * Math.sin(angle);
    }
    return coords;
  }

  private static double[][] buildFused56Coords() {
    double[][] hex = regularPolygonCoords(6);
    double[][] coords = new double[9][2];
    System.arraycopy(hex, 0, coords, 0, 6);
    double mx = (hex[2][0] + hex[3][0]) * 0.5;
    double my = (hex[2][1] + hex[3][1]) * 0.5;
    double edx = hex[3][0] - hex[2][0], edy = hex[3][1] - hex[2][1];
    double nx = -edy, ny = edx;
    double nlen = Math.sqrt(nx * nx + ny * ny);
    nx /= nlen; ny /= nlen;
    double cax = 0, cay = 0;
    for (double[] c : hex) { cax += c[0]; cay += c[1]; }
    cax /= 6; cay /= 6;
    if (nx * (cax - mx) + ny * (cay - my) > 0) { nx = -nx; ny = -ny; }
    double rB = 0.5 / Math.sin(Math.PI / 5.0);
    double cd = rB * Math.cos(Math.PI / 5.0);
    double cx = mx + cd * nx, cy = my + cd * ny;
    double theta = Math.atan2(edy, edx);
    for (int k = 0; k < 3; k++) {
      double angle = theta + Math.PI + (2.0 * Math.PI * (k + 1)) / 5.0;
      coords[6 + k][0] = cx + rB * Math.cos(angle);
      coords[6 + k][1] = cy + rB * Math.sin(angle);
    }
    return coords;
  }

  private static double[][] buildBiphenylCoords() {
    double[][] hex = regularPolygonCoords(6);
    double[][] coords = new double[12][2];
    System.arraycopy(hex, 0, coords, 0, 6);
    double rB = 0.5 / Math.sin(Math.PI / 6.0);
    double offset = hex[0][0] - hex[3][0] + 1.0 + 1.5;
    for (int i = 0; i < 6; i++) {
      double angle = 2.0 * Math.PI * i / 6 - Math.PI / 2.0;
      coords[6 + i][0] = rB * Math.cos(angle) + offset;
      coords[6 + i][1] = rB * Math.sin(angle);
    }
    return coords;
  }

  private static int[][] buildBondList(MolGraph g) {
    List<int[]> bonds = new ArrayList<>();
    for (int i = 0; i < g.n; i++)
      for (int j : g.neighbors[i])
        if (j > i) bonds.add(new int[]{i, j});
    return bonds.toArray(new int[0][]);
  }

  private static boolean segmentsCross(Point2D p1, Point2D p2, Point2D p3, Point2D p4) {
    double d1x = p2.x - p1.x, d1y = p2.y - p1.y;
    double d2x = p4.x - p3.x, d2y = p4.y - p3.y;
    double denom = d1x * d2y - d1y * d2x;
    if (Math.abs(denom) < 1e-12) return false;
    double dx = p3.x - p1.x, dy = p3.y - p1.y;
    double t = (dx * d2y - dy * d2x) / denom;
    double u = (dx * d1y - dy * d1x) / denom;
    double eps = 1e-9;
    return (t > eps && t < 1.0 - eps && u > eps && u < 1.0 - eps);
  }

  private static int countCrossingsImpl(int[][] bonds, Point2D[] coords) {
    int crossings = 0, nb = bonds.length;
    for (int i = 0; i < nb; i++) {
      for (int j = i + 1; j < nb; j++) {
        int a1 = bonds[i][0], a2 = bonds[i][1];
        int b1 = bonds[j][0], b2 = bonds[j][1];
        if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2) continue;
        if (segmentsCross(coords[a1], coords[a2], coords[b1], coords[b2]))
          crossings++;
      }
    }
    return crossings;
  }

  private static List<Set<Integer>> findRingSystemSets(MolGraph g) {
    int[][] raw = g.computeRings();
    if (raw.length == 0) return Collections.emptyList();

    int nr = raw.length;
    List<Set<Long>> edgeSets = new ArrayList<>(nr);
    for (int[] r : raw) {
      Set<Long> es = new HashSet<>();
      for (int j = 0; j < r.length; j++) {
        int a = r[j], b = r[(j + 1) % r.length];
        es.add(bondKey(Math.min(a, b), Math.max(a, b)));
      }
      edgeSets.add(es);
    }

    List<List<Integer>> adj = new ArrayList<>(nr);
    for (int i = 0; i < nr; i++) adj.add(new ArrayList<>());
    for (int i = 0; i < nr; i++) {
      for (int j = i + 1; j < nr; j++) {
        for (Long e : edgeSets.get(i)) {
          if (edgeSets.get(j).contains(e)) {
            adj.get(i).add(j); adj.get(j).add(i); break;
          }
        }
      }
    }

    int[] comp = new int[nr];
    Arrays.fill(comp, -1);
    int nComp = 0;
    for (int i = 0; i < nr; i++) {
      if (comp[i] >= 0) continue;
      comp[i] = nComp;
      Deque<Integer> q = new ArrayDeque<>();
      q.add(i);
      while (!q.isEmpty()) {
        int u = q.poll();
        for (int v : adj.get(u)) {
          if (comp[v] < 0) { comp[v] = nComp; q.add(v); }
        }
      }
      nComp++;
    }

    List<Set<Integer>> systems = new ArrayList<>(nComp);
    for (int i = 0; i < nComp; i++) systems.add(new HashSet<>());
    for (int i = 0; i < nr; i++)
      for (int a : raw[i]) systems.get(comp[i]).add(a);
    return systems;
  }

  /** Reverse-lookup: given an edge index, find the (u,v) pair. */
  private static int[] edgeFromIndex(int idx, Map<Long, Integer> edgeIndex) {
    for (Map.Entry<Long, Integer> e : edgeIndex.entrySet()) {
      if (e.getValue() == idx) {
        long key = e.getKey();
        return new int[] { (int)(key >>> 32), (int)(key & 0xFFFFFFFFL) };
      }
    }
    return null;
  }

  private static boolean sharesVertex(int[] e1, int[] e2) {
    return e1[0] == e2[0] || e1[0] == e2[1] || e1[1] == e2[0] || e1[1] == e2[1];
  }

  private static int find(int[] parent, int i) {
    while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
    return i;
  }

  private static void union(int[] parent, int a, int b) {
    a = find(parent, a); b = find(parent, b);
    if (a != b) parent[a] = b;
  }

  /** Generates a lexicographically minimal orbit signature for a cycle. */
  private String getCycleOrbitSignature(int[] c) {
    int len = c.length;
    int[] bestSeq = null;
    for (int dir = 0; dir <= 1; dir++) {
      for (int start = 0; start < len; start++) {
        int[] seq = new int[len];
        for (int i = 0; i < len; i++) {
          int idx = dir == 0 ? (start + i) % len : (start - i + len) % len;
          seq[i] = orbit[c[idx]];
        }
        if (bestSeq == null || compareSeq(seq, bestSeq) < 0) bestSeq = seq;
      }
    }
    return Arrays.toString(bestSeq);
  }

  private static int compareSeq(int[] a, int[] b) {
    for (int i = 0; i < a.length; i++) {
      if (a[i] != b[i]) return Integer.compare(a[i], b[i]);
    }
    return 0;
  }

  int[][] getNLF1() {
    if (cachedNLF1 == null) cachedNLF1 = buildAllNLF(this, MolGraph::buildNLF1);
    return cachedNLF1;
  }

  /**
   * Cached neighbors sorted by descending degree for VF2++ query ordering.
   * Computed once per MolGraph, reused across all matcher constructions.
   * @since 6.5.3
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
   * @since 6.5.3
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

  int dbStereo(int i, int j) {
    return dbStereoConf != null ? dbStereoConf[i][j] : 0;
  }

  /** BFS shortest path from start to target, avoiding the edge (avoidU, avoidV). */
  private List<Integer> bfsPath(int start, int target, int avoidU, int avoidV) {
    if (start == target) {
      List<Integer> single = new ArrayList<>();
      single.add(start);
      return single;
    }
    int[] parent = new int[n];
    Arrays.fill(parent, -1);
    boolean[] vis = new boolean[n];
    vis[start] = true;
    Deque<Integer> q = new ArrayDeque<>();
    q.add(start);
    while (!q.isEmpty()) {
      int u = q.poll();
      for (int v : neighbors[u]) {
        if ((u == avoidU && v == avoidV) || (u == avoidV && v == avoidU)) continue;
        if (vis[v]) continue;
        vis[v] = true;
        parent[v] = u;
        if (v == target) {
          List<Integer> path = new ArrayList<>();
          for (int c = target; c != -1; c = parent[c]) path.add(c);
          java.util.Collections.reverse(path);
          return path;
        }
        q.add(v);
      }
    }
    return java.util.Collections.emptyList();
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
    int[] props = sparseBondProps.get(bondKey(i, j));
    return props != null ? props[0] : 0;
  }

  boolean bondInRing(int i, int j) {
    if (bondRingMatrix != null) return bondRingMatrix[i][j];
    int[] props = sparseBondProps.get(bondKey(i, j));
    return props != null && props[1] != 0;
  }

  boolean bondAromatic(int i, int j) {
    if (bondAromMatrix != null) return bondAromMatrix[i][j];
    int[] props = sparseBondProps.get(bondKey(i, j));
    return props != null && props[2] != 0;
  }

  public MolGraph(IAtomContainer mol) {
    this.mol = mol;
    this.n = mol.getAtomCount();
    this.words = (n + 63) >>> 6;
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
      try { Cycles.markRingAtomsAndBonds(mol); } catch (Exception _) {}
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

    // Build neighbors and bond properties in a single bond-iteration pass.
    // Avoids per-atom getConnectedAtomsList() calls (each allocates an ArrayList).
    int[] degreeCount = new int[n];
    for (IBond b : mol.bonds()) {
      degreeCount[idxMap.get(b.getAtom(0))]++;
      degreeCount[idxMap.get(b.getAtom(1))]++;
    }
    this.neighbors = new int[n][];
    for (int i = 0; i < n; i++) { neighbors[i] = new int[degreeCount[i]]; degree[i] = 0; }

    if (n <= SPARSE_THRESHOLD) {
      this.bondOrdMatrix = new int[n][n];
      this.bondRingMatrix = new boolean[n][n];
      this.bondAromMatrix = new boolean[n][n];
      this.sparseBondProps = null;
      for (IBond b : mol.bonds()) {
        int a = idxMap.get(b.getAtom(0)), c = idxMap.get(b.getAtom(1));
        neighbors[a][degree[a]++] = c;
        neighbors[c][degree[c]++] = a;
        int ord = encodeBondOrder(b);
        boolean inRing = b.isInRing(), arom = b.isAromatic();
        bondOrdMatrix[a][c] = ord; bondOrdMatrix[c][a] = ord;
        bondRingMatrix[a][c] = inRing; bondRingMatrix[c][a] = inRing;
        bondAromMatrix[a][c] = arom; bondAromMatrix[c][a] = arom;
      }
    } else {
      this.bondOrdMatrix = null;
      this.bondRingMatrix = null;
      this.bondAromMatrix = null;
      this.sparseBondProps = new HashMap<>(mol.getBondCount() * 3);
      for (IBond b : mol.bonds()) {
        int a = idxMap.get(b.getAtom(0)), c = idxMap.get(b.getAtom(1));
        neighbors[a][degree[a]++] = c;
        neighbors[c][degree[c]++] = a;
        sparseBondProps.put(bondKey(a, c),
            new int[] {encodeBondOrder(b), b.isInRing() ? 1 : 0, b.isAromatic() ? 1 : 0});
      }
    }

    // Stereo: only allocate if the molecule actually has stereo elements
    @SuppressWarnings("unchecked")
    Iterable<IStereoElement<?, ?>> stereoElements = (Iterable) mol.stereoElements();
    boolean hasStereo = stereoElements.iterator().hasNext();
    this.tetraChirality = hasStereo ? new int[n] : null;
    this.dbStereoConf = (hasStereo && n <= SPARSE_THRESHOLD) ? new int[n][n] : null;
    if (hasStereo) {
      for (IStereoElement<?, ?> se : stereoElements) {
        if (se instanceof ITetrahedralChirality tc) {
          Integer idx = idxMap.get(tc.getChiralAtom());
          if (idx != null)
            tetraChirality[idx] = tc.getStereo() == ITetrahedralChirality.Stereo.CLOCKWISE ? 1 : 2;
        } else if (se instanceof IDoubleBondStereochemistry dbs) {
          IBond stereoBond = dbs.getStereoBond();
          Integer a = idxMap.get(stereoBond.getAtom(0)), c = idxMap.get(stereoBond.getAtom(1));
          if (a != null && c != null && dbStereoConf != null) {
            int conf = dbs.getStereo() == IDoubleBondStereochemistry.Conformation.TOGETHER ? 1 : 2;
            dbStereoConf[a][c] = conf;
            dbStereoConf[c][a] = conf;
          }
        }
      }
    }

    this.ringCount = new int[n];
    initDerivedFields();
  }

  MolGraph(Builder b) {
    this.mol = null;
    this.n = b.n;
    this.words = (n + 63) >>> 6;
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

    if (n <= SPARSE_THRESHOLD) {
      this.bondOrdMatrix = new int[n][n];
      this.bondRingMatrix = new boolean[n][n];
      this.bondAromMatrix = new boolean[n][n];
      this.sparseBondProps = null;
      for (int i = 0; i < n; i++) {
        int[] nbs = neighbors[i];
        int[] ords = b.bondOrders != null && b.bondOrders[i] != null ? b.bondOrders[i] : null;
        boolean[] rings = b.bondRings != null && b.bondRings[i] != null ? b.bondRings[i] : null;
        boolean[] aroms = b.bondAroms != null && b.bondAroms[i] != null ? b.bondAroms[i] : null;
        boolean denseOrd = ords != null && ords.length == n;
        boolean denseRing = rings != null && rings.length == n;
        boolean denseArom = aroms != null && aroms.length == n;
        for (int k = 0; k < nbs.length; k++) {
          int j = nbs[k];
          bondOrdMatrix[i][j] = ords != null ? (denseOrd ? ords[j] : ords[k]) : 1;
          bondRingMatrix[i][j] = rings != null && (denseRing ? rings[j] : rings[k]);
          bondAromMatrix[i][j] = aroms != null && (denseArom ? aroms[j] : aroms[k]);
        }
      }
    } else {
      this.bondOrdMatrix = null;
      this.bondRingMatrix = null;
      this.bondAromMatrix = null;
      this.sparseBondProps = new HashMap<>();
      for (int i = 0; i < n; i++) {
        int[] nbs = neighbors[i];
        int[] ords = b.bondOrders != null && b.bondOrders[i] != null ? b.bondOrders[i] : null;
        boolean[] rings = b.bondRings != null && b.bondRings[i] != null ? b.bondRings[i] : null;
        boolean[] aroms = b.bondAroms != null && b.bondAroms[i] != null ? b.bondAroms[i] : null;
        boolean denseOrd = ords != null && ords.length == n;
        boolean denseRing = rings != null && rings.length == n;
        boolean denseArom = aroms != null && aroms.length == n;
        for (int k = 0; k < nbs.length; k++) {
          int j = nbs[k];
          if (i < j) {
            sparseBondProps.put(bondKey(i, j), new int[] {
                ords != null ? (denseOrd ? ords[j] : ords[k]) : 1,
                rings != null && (denseRing ? rings[j] : rings[k]) ? 1 : 0,
                aroms != null && (denseArom ? aroms[j] : aroms[k]) ? 1 : 0});
          }
        }
      }
    }

    this.tetraChirality = b.tetraChirality != null ? b.tetraChirality.clone() : new int[n];
    this.atomId = b.atomIds != null ? b.atomIds.clone() : null;
    if (n <= SPARSE_THRESHOLD) {
      this.dbStereoConf = new int[n][n];
      if (b.dbStereoConf != null)
        for (int i = 0; i < n; i++)
          if (b.dbStereoConf[i] != null)
            System.arraycopy(b.dbStereoConf[i], 0, dbStereoConf[i], 0, n);
    } else {
      this.dbStereoConf = null;
    }

    this.ringCount = new int[n];
    initDerivedFields();
  }

  private void initDerivedFields() {
    // adj (BitSet[]) is lazy — built on first getAdj() call.
    // Only adjLong is eagerly built (used by substructure/MCS hot paths).
    this.adjLong = new long[n][words];
    for (int i = 0; i < n; i++)
      for (int k : neighbors[i]) adjLong[i][k >>> 6] |= 1L << (k & 63);
    // Tautomer classes are O(n), cheap, and accessed as public fields — always compute eagerly.
    computeTautomerClasses();
    applySolventCorrection(solvent);
    adjustWeightsForPH(pH);
    tautomerClassesComputed = true;
    // Canonical labeling and ring counts are lazy (ensureCanonical / ensureRingCounts).
  }

  /** Lazily build BitSet[] adjacency (used by MCS refinement and NLF2/NLF3 builders). */
  BitSet[] getAdj() {
    if (adj == null) {
      BitSet[] a = new BitSet[n];
      for (int i = 0; i < n; i++) {
        a[i] = new BitSet(n);
        for (int k : neighbors[i]) a[i].set(k);
      }
      this.adj = a;
    }
    return adj;
  }

  /** Lazily compute Morgan ranks, canonical labeling, orbit partition, and canonical hash. */
  void ensureCanonical() {
    if (canonicalComputed) return;
    synchronized (this) {
      if (canonicalComputed) return;
      this.morganRank = computeMorganRanks(n, label, neighbors);
      CanonResult clResult = computeCanonicalLabeling(n, label, degree, neighbors);
      this.canonicalLabel = clResult.canonLabel();
      this.orbit = clResult.orbit();
      this.autGenerators = clResult.autGenerators();
      this.autGeneratorsTruncated = clResult.generatorsTruncated();
      this.canonicalHash = computeCanonicalHash(this, canonicalLabel);
      canonicalComputed = true;
    }
  }

  /** Lazily compute tautomer class assignments and pKa-based relevance weights. */
  void ensureTautomerClasses() {
    if (tautomerClassesComputed) return;
    synchronized (this) {
      if (tautomerClassesComputed) return;
      computeTautomerClasses();
      tautomerClassesComputed = true;
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

  private void computeRingCounts() {
    // Use Relevant Cycle Basis (unique/deterministic) instead of SSSR (non-deterministic).
    // This ensures ringCount is consistent across runs for symmetric molecules like cubane.
    java.util.Arrays.fill(ringCount, 0);
    int[][] relevant = computeRelevantCycles();
    for (int[] r : relevant) {
      for (int atom : r) ringCount[atom]++;
    }
  }

  /**
   * Assign a tautomer class to a set of atoms with a shared relevance weight.
   * If any atom already belongs to an existing class, merge by reusing that
   * class for all atoms in the group (prevents class fracture when
   * tautomeric motifs overlap on shared atoms).
   */
  private void tc(int[] atoms, int cls, float w, float[] wArr) {
    // Merge: if any atom already has a class, reuse it
    for (int a : atoms) if (tautomerClass[a] != -1) { cls = tautomerClass[a]; break; }
    for (int a : atoms) {
      tautomerClass[a] = cls;
      if (wArr[a] == 1.0f) wArr[a] = w;
    }
  }

  void computeTautomerClasses() {
    tautomerClass  = new int[n];
    tautomerWeight = new float[n];
    Arrays.fill(tautomerClass,  -1);
    Arrays.fill(tautomerWeight, 1.0f);
    int classId = 0;

    for (int i = 0; i < n; i++) {
      if (tautomerClass[i] != -1) continue;

      // -----------------------------------------------------------------------
      // T1: Keto/enol  C=O ↔ C-OH  (aliphatic, non-amide)
      //     Detects O of a carbonyl; groups O + adjacent C neighbours.
      //     Weight TW_KETO_ENOL: keto overwhelmingly dominant (pKa ~20).
      //     Amide carbonyls (N on C) handled separately in T3/T4 with lower weight.
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i]) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            boolean hasN = false;
            for (int k : neighbors[j]) { if (k != i && atomicNum[k] == 7) { hasN = true; break; } }
            float w = hasN ? TW_AMIDE_IMIDIC : TW_KETO_ENOL;
            int cls = -1;
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 6 && !aromatic[k]) {
                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k]==1.0f) tautomerWeight[k]=w; }
              }
            }
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 7) {
                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k]==1.0f) tautomerWeight[k]=w; }
              }
            }
            if (cls != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T2: Thione/thiol  C=S ↔ C-SH  (thione dominant, pKa ~10-11)
      //     Thioamide N-C(=S) ↔ N=C-SH gets higher weight (pKa ~13)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 16 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            boolean hasN = false;
            for (int k : neighbors[j]) { if (k != i && atomicNum[k] == 7) { hasN = true; break; } }
            float w = hasN ? TW_THIOAMIDE : TW_THIONE_THIOL;
            int cls = -1;
            for (int k : neighbors[j]) {
              if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7)) {
                if (cls == -1) { cls = classId++; tautomerClass[i] = cls; tautomerWeight[i] = w; }
                if (tautomerClass[k] == -1) { tautomerClass[k] = cls; if (tautomerWeight[k]==1.0f) tautomerWeight[k]=w; }
              }
            }
            if (cls != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T3: Nitroso/oxime  C=N-OH ↔ C(=O)-NH  (oxime slightly more stable)
      //     Requires N=C double bond to avoid matching simple amide N-C=O
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && !ring[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8) {
                int cls = classId++;
                tc(new int[]{i, k}, cls, TW_NITROSO_OXIME, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T4: Phenol/quinone  arom-C-OH ↔ para-C=O  (phenol dominant, pKa ~10)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
            for (int a : neighbors[j]) {
              if (a != i && atomicNum[a] == 6 && aromatic[a] && ring[a]) {
                for (int b : neighbors[a]) {
                  if (b != j && atomicNum[b] == 6 && aromatic[b] && ring[b]) {
                    for (int c : neighbors[b]) {
                      if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) {
                        int cls = classId++;
                        tc(new int[]{i, j, c}, cls, TW_PHENOL_QUINONE, tautomerWeight);
                        for (int d : neighbors[c]) {
                          if (atomicNum[d] == 8 && !aromatic[d] && tautomerClass[d] == -1) {
                            tautomerClass[d] = cls; tautomerWeight[d] = TW_PHENOL_QUINONE;
                          }
                        }
                        break;
                      }
                    }
                    if (tautomerClass[i] != -1) break;
                  }
                }
                if (tautomerClass[i] != -1) break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T5: 1,3-diketone enolisation  C(=O)-C-C(=O)  (pKa ~8-11, pH-sensitive)
      //     Requires exactly one carbon bridge between the two carbonyls
      //     (true 1,3-pattern: O=C-CH-C=O, not fused/directly bonded carbonyls)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            for (int bridge : neighbors[j]) {
              if (bridge != i && atomicNum[bridge] == 6 && !aromatic[bridge]) {
                for (int m : neighbors[bridge]) {
                  if (m != j && m != i && atomicNum[m] == 6) {
                    // Verify true 1,3-separation: m must NOT be a direct
                    // neighbor of j (which would mean j and m are 1,2 not 1,3).
                    boolean directlyBonded = false;
                    for (int nb : neighbors[j]) {
                      if (nb == m) { directlyBonded = true; break; }
                    }
                    if (directlyBonded) continue;
                    for (int o2 : neighbors[m]) {
                      if (o2 != bridge && atomicNum[o2] == 8 && !aromatic[o2]
                          && bondOrder(m, o2) == 2) {
                        int cls = (tautomerClass[o2] != -1) ? tautomerClass[o2] : classId++;
                        float w = TW_DIKETONE_ENOL;
                        if (tautomerClass[i]      == -1) { tautomerClass[i]      = cls; tautomerWeight[i]      = w; }
                        if (tautomerClass[bridge] == -1) { tautomerClass[bridge] = cls; tautomerWeight[bridge] = w; }
                        if (tautomerClass[o2]     == -1) { tautomerClass[o2]     = cls; tautomerWeight[o2]     = w; }
                        break;
                      }
                    }
                    if (tautomerClass[i] != -1) break;
                  }
                }
                if (tautomerClass[i] != -1) break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T6: Aromatic N adjacent to exocyclic C=O  (lactam/2-pyridinone/uracil)
      //     Lactam form dominant at pH 7.4 (pKa ~1 for 2-pyridinone)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j]) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8 && !aromatic[k] && tautomerClass[k] == -1) {
                int cls = classId++;
                tc(new int[]{i, k}, cls, TW_PYRIDINONE, tautomerWeight);
              }
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T7: Aromatic N–N or N–C–N  (imidazole, pyrazole, benzimidazole, indazole)
      //     N1H/N3H symmetric tautomers — weight 0.50
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        // Direct N–N bond (pyrazole/indazole)
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 7 && aromatic[j] && ring[j] && tautomerClass[j] == -1) {
            int cls = classId++;
            tc(new int[]{i, j}, cls, TW_IMIDAZOLE_NH, tautomerWeight);
            break;
          }
        }
        // N–C–N through aromatic C (imidazole/benzimidazole)
        // Only match when both N atoms share a 5-membered ring (skip pyrimidine)
        if (tautomerClass[i] == -1) {
          int[][] rings = computeRings();
          for (int j : neighbors[i]) {
            if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
              for (int k : neighbors[j]) {
                if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
                  // Verify both N atoms share a 5-membered ring
                  boolean in5 = false;
                  for (int[] r : rings) {
                    if (r.length != 5) continue;
                    boolean hasI = false, hasK = false;
                    for (int a : r) { if (a == i) hasI = true; if (a == k) hasK = true; }
                    if (hasI && hasK) { in5 = true; break; }
                  }
                  if (!in5) continue;
                  int cls = classId++;
                  tc(new int[]{i, k}, cls, TW_IMIDAZOLE_NH, tautomerWeight);
                  break;
                }
              }
              if (tautomerClass[i] != -1) break;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T8: Amidine  N=C-N ↔ N-C=N  (non-aromatic, symmetric, weight 0.50)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 7 && !aromatic[k]) {
                int cls = classId++;
                tc(new int[]{i, j, k}, cls, TW_AMIDINE, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T9: Guanidine  C connected to ≥3 N (symmetric, weight 0.50)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 6 && tautomerClass[i] == -1) {
        List<Integer> nNbrs = new ArrayList<>();
        for (int j : neighbors[i]) { if (atomicNum[j] == 7 && !aromatic[j]) nNbrs.add(j); }
        if (nNbrs.size() >= 3) {
          int cls = classId++;
          tautomerClass[i] = cls; tautomerWeight[i] = TW_GUANIDINE;
          for (int ni : nNbrs) { if (tautomerClass[ni]==-1) { tautomerClass[ni]=cls; tautomerWeight[ni]=TW_GUANIDINE; } }
        }
      }

      // -----------------------------------------------------------------------
      // T10: Urea  O=C(-N)(-N)  — two N on same carbonyl C (dominant, pKa ~14)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            List<Integer> nn = new ArrayList<>();
            for (int k : neighbors[j]) { if (k != i && atomicNum[k] == 7) nn.add(k); }
            if (nn.size() >= 2) {
              int cls = classId++;
              tautomerClass[i] = cls; tautomerWeight[i] = TW_UREA;
              tautomerClass[j] = cls; tautomerWeight[j] = TW_UREA;
              for (int ni : nn) { if (tautomerClass[ni]==-1) { tautomerClass[ni]=cls; tautomerWeight[ni]=TW_UREA; } }
            }
            break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T11: β-Enaminone  N-C=C-C=O  (vinylogous amide, weight 0.85)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
        outer11:
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 1) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 2) {
                for (int m : neighbors[k]) {
                  if (m != j && atomicNum[m] == 6) {
                    for (int o : neighbors[m]) {
                      if (o != k && atomicNum[o] == 8 && bondOrder(m, o) == 2) {
                        int cls = classId++;
                        tc(new int[]{i, j, k, m, o}, cls, TW_ENAMINONE, tautomerWeight);
                        break outer11;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T12: Imide  C(=O)-N-C(=O)  (N flanked by two carbonyls, pKa ~9-10)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && tautomerClass[i] == -1) {
        List<Integer> carbonyls = new ArrayList<>();
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) { carbonyls.add(j); break; }
            }
          }
        }
        if (carbonyls.size() >= 2) {
          int cls = classId++;
          tautomerClass[i] = cls; tautomerWeight[i] = TW_IMIDE;
          for (int c : carbonyls) { if (tautomerClass[c]==-1) { tautomerClass[c]=cls; tautomerWeight[c]=TW_IMIDE; } }
        }
      }

      // -----------------------------------------------------------------------
      // T13: Hydroxamic acid  N(-OH)-C=O ↔ NH-C(=O)-OH  (pKa ~8, weight 0.85)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && tautomerClass[i] == -1) {
        int oh = -1;
        for (int j : neighbors[i]) { if (atomicNum[j] == 8 && !aromatic[j] && bondOrder(i,j)==1) { oh = j; break; } }
        if (oh != -1) {
          for (int j : neighbors[i]) {
            if (atomicNum[j] == 6) {
              for (int k : neighbors[j]) {
                if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) {
                  int cls = classId++;
                  tc(new int[]{i, oh, j, k}, cls, TW_HYDROXAMIC, tautomerWeight);
                  break;
                }
              }
              if (tautomerClass[i] != -1) break;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T14: Hydroxypyrimidine/purine  exocyclic OH on arom-C with ring-N neighbour
      //      Covers adenine/guanine amino-imino, hydroxypyrimidine, cytosine
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j] && ring[j] && bondOrder(i, j) == 1) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k]) {
                int cls = classId++;
                tc(new int[]{i, j, k}, cls, TW_HYDROXYPYRIMIDINE, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T15: Tetrazole / triazole  N-H tautomers in aromatic rings with ≥2 N-N bonds
      //      (1H/2H symmetric, weight 0.50)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        List<Integer> adjN = new ArrayList<>();
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 7 && aromatic[j] && ring[j]) adjN.add(j);
        }
        if (adjN.size() >= 2) { // tetrazole (≥3N in ring); also covers 1,2,3-triazole
          for (int j : adjN) {
            if (tautomerClass[j] == -1) {
              int cls = classId++;
              tc(new int[]{i, j}, cls, TW_TETRAZOLE_NH, tautomerWeight);
              break;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T16: Thioamide/iminothiol  RC(=S)NR ↔ RC(SH)=NR  (weight 0.95)
      //      S on C with double bond to S and single bond to N
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 16 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 7) {
                int cls = classId++;
                tc(new int[]{i, j, k}, cls, TW_THIOAMIDE_IMINOTHIOL, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T17: Phosphonate ester  P(=O)(OH) ↔ P(OH)(=O)  (symmetric, weight 0.50)
      //      P with both =O and -OH
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 15 && tautomerClass[i] == -1) {
        List<Integer> oxygens = new ArrayList<>();
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 8) oxygens.add(j);
        }
        if (oxygens.size() >= 2) {
          int cls = classId++;
          tautomerClass[i] = cls; tautomerWeight[i] = TW_PHOSPHONATE;
          for (int o : oxygens) {
            if (tautomerClass[o] == -1) { tautomerClass[o] = cls; tautomerWeight[o] = TW_PHOSPHONATE; }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T19: Nitro/aci-nitro  C-N(=O)=O ↔ C=N(=O)-OH  (weight 0.95)
      //      N with two O neighbors (one double, one single/double) bonded to C
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && !aromatic[i] && tautomerClass[i] == -1) {
        List<Integer> oNbrs = new ArrayList<>();
        int cNbr = -1;
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 8) oNbrs.add(j);
          else if (atomicNum[j] == 6) cNbr = j;
        }
        if (oNbrs.size() >= 2 && cNbr != -1) {
          int cls = classId++;
          tautomerClass[i] = cls; tautomerWeight[i] = TW_NITRO_ACI_NITRO;
          tautomerClass[cNbr] = cls; tautomerWeight[cNbr] = TW_NITRO_ACI_NITRO;
          for (int o : oNbrs) {
            if (tautomerClass[o] == -1) { tautomerClass[o] = cls; tautomerWeight[o] = TW_NITRO_ACI_NITRO; }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T21: 1,5-keto-enol  C(=O)-C=C-C=C <-> C(-OH)-C=C-C=C  (weight 0.75)
      //      5-atom keto-enol shift through conjugation
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 1) {
                for (int m : neighbors[k]) {
                  if (m != j && atomicNum[m] == 6 && bondOrder(k, m) == 2) {
                    for (int p : neighbors[m]) {
                      if (p != k && atomicNum[p] == 6 && !aromatic[p]) {
                        int cls = classId++;
                        tc(new int[]{i, j, k, m, p}, cls, TW_15_KETO_ENOL, tautomerWeight);
                        break;
                      }
                    }
                    if (tautomerClass[i] != -1) break;
                  }
                }
                if (tautomerClass[i] != -1) break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T22: Ring-chain: furanose/pyranose  open-chain ↔ cyclic hemiacetal  (weight 0.60)
      //      O in ring bonded to C with OH neighbor (hemiacetal pattern)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && ring[i] && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && ring[j]) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8 && !ring[k] && bondOrder(j, k) == 1) {
                int cls = classId++;
                tc(new int[]{i, j, k}, cls, TW_FURANOSE_PYRANOSE, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T23: Ring-chain: lactol  open-chain hydroxy-aldehyde ↔ cyclic lactol  (weight 0.65)
      //      C in ring bonded to both ring-O and exocyclic OH
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 6 && ring[i] && tautomerClass[i] == -1) {
        int ringO = -1, exoOH = -1;
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 8 && ring[j] && !aromatic[j]) ringO = j;
          else if (atomicNum[j] == 8 && !ring[j] && !aromatic[j] && bondOrder(i, j) == 1) exoOH = j;
        }
        if (ringO != -1 && exoOH != -1) {
          int cls = classId++;
          tc(new int[]{ringO, i, exoOH}, cls, TW_LACTOL, tautomerWeight);
        }
      }

      // -----------------------------------------------------------------------
      // T24: Phenol/quinone methide  phenol → para-quinone methide  (weight 0.92)
      //      Aromatic C-OH where para-C has exocyclic =CH2
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 8 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j] && ring[j] && bondOrder(i, j) == 1) {
            // Walk 3 aromatic carbons to find para position
            for (int a : neighbors[j]) {
              if (a != i && atomicNum[a] == 6 && aromatic[a] && ring[a]) {
                for (int b : neighbors[a]) {
                  if (b != j && atomicNum[b] == 6 && aromatic[b] && ring[b]) {
                    for (int c : neighbors[b]) {
                      if (c != a && c != j && atomicNum[c] == 6 && aromatic[c] && ring[c]) {
                        // Check if para-C has an exocyclic C=C (methide)
                        for (int d : neighbors[c]) {
                          if (atomicNum[d] == 6 && !ring[d] && bondOrder(c, d) == 2) {
                            int cls = classId++;
                            tc(new int[]{i, j, c, d}, cls, TW_PHENOL_QUINONE_METHIDE, tautomerWeight);
                            break;
                          }
                        }
                        if (tautomerClass[i] != -1) break;
                      }
                    }
                    if (tautomerClass[i] != -1) break;
                  }
                }
                if (tautomerClass[i] != -1) break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T25: Triazole NH shift  1,2,3-triazole NH shift  (symmetric, weight 0.50)
      //      Aromatic N in ring with exactly one adjacent aromatic N (not covered by T15 tetrazole)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        List<Integer> adjNring = new ArrayList<>();
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 7 && aromatic[j] && ring[j]) adjNring.add(j);
        }
        if (adjNring.size() == 1) {
          int j = adjNring.get(0);
          // Check if the N-N bond partner also has another ring-N neighbor (triazole pattern)
          for (int k : neighbors[j]) {
            if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
              int cls = classId++;
              tc(new int[]{i, j, k}, cls, TW_TRIAZOLE_NH_SHIFT, tautomerWeight);
              break;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T26: Benzimidazole NH shift  NH-1 ↔ NH-3 shift  (symmetric, weight 0.50)
      //      Two aromatic N bridged by aromatic C in a fused ring system (distinct from T7 by ring count)
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
            // Check if this C is in a fused ring (degree >= 3 among ring atoms)
            int ringNbCount = 0;
            for (int k : neighbors[j]) { if (ring[k] && aromatic[k]) ringNbCount++; }
            if (ringNbCount >= 3) {
              for (int k : neighbors[j]) {
                if (k != i && atomicNum[k] == 7 && aromatic[k] && ring[k] && tautomerClass[k] == -1) {
                  int cls = classId++;
                  tc(new int[]{i, j, k}, cls, TW_BENZIMIDAZOLE_NH, tautomerWeight);
                  break;
                }
              }
              if (tautomerClass[i] != -1) break;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T27: Pyridone/hydroxypyridine  2-pyridinol ↔ 2-pyridone  (weight 0.85)
      //      Aromatic N in 6-membered ring with C-OH at position 2
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && aromatic[i] && ring[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && aromatic[j] && ring[j]) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8 && !aromatic[k] && bondOrder(j, k) == 1) {
                int cls = classId++;
                tc(new int[]{i, j, k}, cls, TW_PYRIDONE_HYDROXYPYRIDINE, tautomerWeight);
                break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T28: Barbituric acid  tri-keto ↔ enol forms  (weight 0.80)
      //      N in ring flanked by two carbonyl C (each C=O), forming barbiturate pattern
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 7 && ring[i] && tautomerClass[i] == -1) {
        List<Integer> carbonylC = new ArrayList<>();
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && ring[j]) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 8 && bondOrder(j, k) == 2) {
                carbonylC.add(j);
                break;
              }
            }
          }
        }
        if (carbonylC.size() >= 2) {
          int cls = classId++;
          tautomerClass[i] = cls; tautomerWeight[i] = TW_BARBITURIC;
          for (int c : carbonylC) {
            if (tautomerClass[c] == -1) { tautomerClass[c] = cls; tautomerWeight[c] = TW_BARBITURIC; }
            for (int k : neighbors[c]) {
              if (atomicNum[k] == 8 && bondOrder(c, k) == 2 && tautomerClass[k] == -1) {
                tautomerClass[k] = cls; tautomerWeight[k] = TW_BARBITURIC;
              }
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // T29: Allyl shift  X-C=C-C ↔ C=C-C-X  where X is OH, SH, NH  (weight 0.70)
      //      Validates contiguous bond path: X(-H)--C==C--C(sp3)
      //      X must carry at least one H (checked via valence/degree)
      //      and k-m must be a single bond (not another double bond system)
      // -----------------------------------------------------------------------
      if ((atomicNum[i] == 8 || atomicNum[i] == 16 || atomicNum[i] == 7)
           && !aromatic[i] && !ring[i] && tautomerClass[i] == -1) {
        // X must bear at least one H for the proton shift to occur.
        int sumBO = 0;
        for (int nb : neighbors[i]) sumBO += bondOrder(i, nb);
        int maxValence = (atomicNum[i] == 7) ? 3 : 2; // N=3, O=2, S=2
        boolean xHasH = (sumBO < maxValence + Math.abs(formalCharge[i]));
        if (!xHasH) continue; // no H on heteroatom, shift not possible
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && !aromatic[j] && bondOrder(i, j) == 1) {
            for (int k : neighbors[j]) {
              if (k != i && atomicNum[k] == 6 && bondOrder(j, k) == 2) {
                for (int m : neighbors[k]) {
                  // Require single bond k-m (sp3 terminal carbon)
                  if (m != j && atomicNum[m] == 6 && !aromatic[m] && bondOrder(k, m) == 1) {
                    int cls = classId++;
                    tc(new int[]{i, j, k, m}, cls, TW_ALLYL_SHIFT, tautomerWeight);
                    break;
                  }
                }
                if (tautomerClass[i] != -1) break;
              }
            }
            if (tautomerClass[i] != -1) break;
          }
        }
      }

      // -----------------------------------------------------------------------
      // T30: Selenol/selenoketone  C=Se ↔ C(-SeH)  (weight 0.94)
      //      Se double-bonded to C
      // -----------------------------------------------------------------------
      if (atomicNum[i] == 34 && !aromatic[i] && tautomerClass[i] == -1) {
        for (int j : neighbors[i]) {
          if (atomicNum[j] == 6 && bondOrder(i, j) == 2) {
            int cls = classId++;
            tc(new int[]{i, j}, cls, TW_SELENOL_SELENOKETONE, tautomerWeight);
            for (int k : neighbors[j]) {
              if (k != i && (atomicNum[k] == 6 || atomicNum[k] == 7) && tautomerClass[k] == -1) {
                tautomerClass[k] = cls; if (tautomerWeight[k] == 1.0f) tautomerWeight[k] = TW_SELENOL_SELENOKETONE;
              }
            }
            break;
          }
        }
      }
    }
  }

  // =========================================================================
  // Tier 2: Solvent-dependent tautomer weight corrections
  //
  // Aprotic/non-aqueous solvents shift keto-enol and related equilibria.
  // Literature pKa shifts (Bordwell 1988; Olmstead 1977):
  //   DMSO:       keto-enol +5 pKa units  -> enol more prevalent (shift -0.15)
  //               amide/lactam +3          -> imidic slightly more prevalent
  //   Methanol:   keto-enol +2             -> small shift (-0.05)
  //   Chloroform: keto-enol +8 (aprotic)   -> much more enol (-0.25)
  // =========================================================================

  /** Adjust tautomer weights for solvent effects. Call after computeTautomerClasses(). */
  void applySolventCorrection(ChemOptions.Solvent solvent) {
    if (solvent == null || solvent == ChemOptions.Solvent.AQUEOUS) return;
    if (tautomerWeight == null) return;

    float ketoShift, amideShift, generalShift;
    switch (solvent) {
      case DMSO:
        ketoShift = -0.15f;
        amideShift = -0.08f;
        generalShift = -0.05f;
        break;
      case CHLOROFORM:
        ketoShift = -0.25f;
        amideShift = -0.12f;
        generalShift = -0.08f;
        break;
      case METHANOL:
        ketoShift = -0.05f;
        amideShift = -0.02f;
        generalShift = -0.01f;
        break;
      case ACETONITRILE:
        ketoShift = -0.10f;
        amideShift = -0.05f;
        generalShift = -0.03f;
        break;
      case DIETHYL_ETHER:
        ketoShift = -0.20f;
        amideShift = -0.10f;
        generalShift = -0.06f;
        break;
      default:
        return;
    }

    for (int i = 0; i < n; i++) {
      if (tautomerClass[i] == -1) continue;
      float w = tautomerWeight[i];
      if (w >= 1.0f) continue;

      float shift;
      if (w >= 0.96f) {
        shift = (w >= 0.97f) ? amideShift : ketoShift;
        if (Math.abs(w - TW_KETO_ENOL) < 0.005f) shift = ketoShift;
      } else if (w >= 0.93f) {
        shift = amideShift;
      } else if (w >= 0.84f) {
        shift = generalShift;
      } else {
        shift = generalShift;
      }

      tautomerWeight[i] = Math.max(0.05f, Math.min(0.99f, w + shift));
    }
  }

  /**
   * Adjust tautomer weights based on pH deviation from reference (7.4).
   * Uses simplified Henderson-Hasselbalch scaling: as pH moves away from
   * physiological, weights shift toward 0.5 (equal tautomer population).
   */
  void adjustWeightsForPH(double pH) {
    if (Math.abs(pH - 7.4) < 0.01) return; // no adjustment at reference pH
    // Dampen weights toward 0.5 as pH deviates from 7.4
    float dampen = (float) (0.1 * Math.abs(pH - 7.4));
    for (int i = 0; i < n; i++) {
      if (tautomerClass[i] == -1) continue;
      float w = tautomerWeight[i];
      if (w == 1.0f) continue; // non-tautomeric atom
      // Move weight toward 0.5 but NEVER cross it (prevents preference inversion)
      float shift = (w > 0.5f) ? -dampen : dampen;
      float adjusted = w + shift;
      // Clamp: never cross 0.5, never go below 0.1 or above 0.99
      if (w > 0.5f) adjusted = Math.max(0.5f, adjusted);
      else           adjusted = Math.min(0.5f, adjusted);
      tautomerWeight[i] = Math.max(0.1f, Math.min(0.99f, adjusted));
    }
  }

  private static int[] computeMorganRanks(int n, int[] initLabel, int[][] neighbors) {
    int[] rank = new int[n];
    System.arraycopy(initLabel, 0, rank, 0, n);
    int[] rankNew = new int[n];
    int[] nbSorted = new int[16];
    for (int iter = 0; iter < 8; iter++) {
      HashSet<Integer> prevSet = new HashSet<>();
      for (int v = 0; v < n; v++) prevSet.add(rank[v]);
      int prevDistinct = prevSet.size();
      for (int v = 0; v < n; v++) {
        int[] nbs = neighbors[v];
        int deg = nbs.length;
        if (deg > nbSorted.length) nbSorted = new int[deg];
        for (int k = 0; k < deg; k++) nbSorted[k] = rank[nbs[k]];
        java.util.Arrays.sort(nbSorted, 0, deg);
        int h = rank[v] * HASH_PRIME;
        for (int k = 0; k < deg; k++) h = h * 31 + nbSorted[k];
        rankNew[v] = h;
      }
      System.arraycopy(rankNew, 0, rank, 0, n);
      HashSet<Integer> newSet = new HashSet<>();
      for (int v = 0; v < n; v++) newSet.add(rank[v]);
      if (newSet.size() <= prevDistinct) break;
    }
    return rank;
  }

  private static final int MAX_SEARCH_NODES = 50000;
  private static final int CANON_SEARCH_LIMIT = 200;
  private static final int MAX_GENERATORS = 500;

  private record CanonResult(int[] canonLabel, int[] orbit, int[][] autGenerators,
      boolean generatorsTruncated) {}

  private static CanonResult computeCanonicalLabeling(int n, int[] label, int[] degree, int[][] neighbors) {
    if (n == 0) return new CanonResult(new int[0], new int[0], new int[0][], false);
    int[] mRank = new int[n];
    System.arraycopy(label, 0, mRank, 0, n);
    int[] mNew = new int[n];
    int[] scratch = new int[16];
    for (int iter = 0; iter < 8; iter++) {
      int prevDistinct;
      { HashSet<Integer> s = new HashSet<>(); for (int v = 0; v < n; v++) s.add(mRank[v]); prevDistinct = s.size(); }
      for (int v = 0; v < n; v++) {
        int[] nbs = neighbors[v];
        int d = nbs.length;
        if (d > scratch.length) scratch = new int[d];
        for (int k = 0; k < d; k++) scratch[k] = mRank[nbs[k]];
        Arrays.sort(scratch, 0, d);
        int h = mRank[v] * HASH_PRIME;
        for (int k = 0; k < d; k++) h = h * 31 + scratch[k];
        mNew[v] = h;
      }
      System.arraycopy(mNew, 0, mRank, 0, n);
      int newDistinct;
      { HashSet<Integer> s = new HashSet<>(); for (int v = 0; v < n; v++) s.add(mRank[v]); newDistinct = s.size(); }
      if (newDistinct <= prevDistinct) break;
    }
    int[] uf = new int[n];
    for (int i = 0; i < n; i++) uf[i] = i;
    Integer[] sorted = new Integer[n];
    for (int i = 0; i < n; i++) sorted[i] = i;
    final int[] mr = mRank;
    Arrays.sort(sorted, (a, b) -> {
      int c = Integer.compare(label[a], label[b]);
      if (c != 0) return c;
      c = Integer.compare(degree[a], degree[b]);
      return c != 0 ? c : Integer.compare(mr[a], mr[b]);
    });
    int[] initPerm = new int[n];
    boolean[] initCellEnd = new boolean[n];
    for (int i = 0; i < n; i++) initPerm[i] = sorted[i];
    for (int i = 0; i < n - 1; i++) {
      int vi = sorted[i], vj = sorted[i + 1];
      initCellEnd[i] = (label[vi] != label[vj] || degree[vi] != degree[vj] || mr[vi] != mr[vj]);
    }
    initCellEnd[n - 1] = true;
    if (n > CANON_SEARCH_LIMIT) {
      int cs = 0;
      for (int i = 0; i < n; i++) {
        if (initCellEnd[i]) {
          for (int j = cs + 1; j <= i; j++) ufUnion(uf, initPerm[cs], initPerm[j]);
          cs = i + 1;
        }
      }
      int[] canonLabel = new int[n];
      for (int pos = 0; pos < n; pos++) canonLabel[initPerm[pos]] = pos;
      return new CanonResult(canonLabel, buildOrbits(uf, n), new int[0][], false);
    }
    refinePartition(n, initPerm, initCellEnd, neighbors, label, -1);
    int[] bestPerm = null;
    List<int[]> generators = new ArrayList<>();
    boolean generatorsTruncated = false;
    int nodeCount = 0;
    boolean budgetExceeded = false;
    ArrayDeque<int[]> permStack = new ArrayDeque<>();
    ArrayDeque<boolean[]> cellEndStack = new ArrayDeque<>();
    int targetCell = findSmallestNonSingleton(n, initCellEnd);
    if (targetCell < 0) {
      bestPerm = initPerm.clone();
    } else {
      int cellStart = 0;
      for (int i = 0; i < n; i++) {
        if (i > 0 && initCellEnd[i - 1]) cellStart = i;
        if (i == targetCell) break;
      }
      for (int i = targetCell; i >= cellStart; i--) {
        int[] p = initPerm.clone();
        boolean[] ce = initCellEnd.clone();
        if (i != cellStart) {
          int tmp = p[i];
          System.arraycopy(p, cellStart, p, cellStart + 1, i - cellStart);
          p[cellStart] = tmp;
        }
        ce[cellStart] = true;
        refinePartition(n, p, ce, neighbors, label, cellStart);
        permStack.push(p);
        cellEndStack.push(ce);
        nodeCount++;
      }
    }
    while (!permStack.isEmpty() && !budgetExceeded) {
      int[] perm = permStack.pop();
      boolean[] cellEnd = cellEndStack.pop();
      int nc = findSmallestNonSingleton(n, cellEnd);
      if (nc < 0) {
        if (bestPerm == null || comparePerm(perm, bestPerm, n, label, neighbors) < 0) {
          bestPerm = perm.clone();
        } else if (comparePerm(perm, bestPerm, n, label, neighbors) == 0) {
          for (int i = 0; i < n; i++) ufUnion(uf, perm[i], bestPerm[i]);
          if (generators.size() < MAX_GENERATORS) {
            int[] gen = new int[n];
            for (int i = 0; i < n; i++) gen[bestPerm[i]] = perm[i];
            generators.add(gen);
          } else {
            generatorsTruncated = true;
          }
        }
        continue;
      }
      if (bestPerm != null && partialCompare(perm, bestPerm, n, cellEnd, label, neighbors) > 0) continue;
      int cellStart = 0;
      for (int i = 0; i < n; i++) {
        if (i > 0 && cellEnd[i - 1]) cellStart = i;
        if (i == nc) break;
      }
      int cellSize = nc - cellStart + 1;
      if (nodeCount + cellSize > MAX_SEARCH_NODES) { budgetExceeded = true; break; }
      for (int i = nc; i >= cellStart; i--) {
        int[] p = perm.clone();
        boolean[] ce = cellEnd.clone();
        if (i != cellStart) {
          int tmp = p[i];
          System.arraycopy(p, cellStart, p, cellStart + 1, i - cellStart);
          p[cellStart] = tmp;
        }
        ce[cellStart] = true;
        refinePartition(n, p, ce, neighbors, label, cellStart);
        permStack.push(p);
        cellEndStack.push(ce);
        nodeCount++;
      }
    }
    if (bestPerm == null) bestPerm = initPerm.clone();
    int[] canonLabel = new int[n];
    for (int pos = 0; pos < n; pos++) canonLabel[bestPerm[pos]] = pos;
    int cs = 0;
    for (int i = 0; i < n; i++) {
      if (initCellEnd[i]) {
        for (int j = cs + 1; j <= i; j++) ufUnion(uf, initPerm[cs], initPerm[j]);
        cs = i + 1;
      }
    }
    return new CanonResult(canonLabel, buildOrbits(uf, n),
        generators.toArray(new int[0][]), generatorsTruncated);
  }

  private static int[] buildOrbits(int[] uf, int n) {
    int[] orbit = new int[n];
    for (int i = 0; i < n; i++) orbit[i] = ufFind(uf, i);
    int[] orbitMin = new int[n];
    Arrays.fill(orbitMin, n);
    for (int i = 0; i < n; i++) { int r = orbit[i]; if (i < orbitMin[r]) orbitMin[r] = i; }
    for (int i = 0; i < n; i++) orbit[i] = orbitMin[orbit[i]];
    return orbit;
  }

  private static int ufFind(int[] uf, int x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
  }

  private static void ufUnion(int[] uf, int a, int b) {
    a = ufFind(uf, a); b = ufFind(uf, b);
    if (a != b) { if (a < b) uf[b] = a; else uf[a] = b; }
  }

  private static int findSmallestNonSingleton(int n, boolean[] cellEnd) {
    int bestSize = n + 1, bestEnd = -1, cellStart = 0;
    for (int i = 0; i < n; i++) {
      if (cellEnd[i]) {
        int size = i - cellStart + 1;
        if (size > 1 && size < bestSize) { bestSize = size; bestEnd = i; }
        cellStart = i + 1;
      }
    }
    return bestEnd;
  }

  private static void refinePartition(
      int n, int[] perm, boolean[] cellEnd, int[][] neighbors, int[] label, int targetPos) {
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[perm[i]] = i;
    boolean changed = true;
    int maxIter = n * 2 + 10;
    while (changed && maxIter-- > 0) {
      changed = false;
      int sStart = 0;
      for (int sEnd = 0; sEnd < n; sEnd++) {
        if (!cellEnd[sEnd]) continue;
        int cStart = 0;
        for (int cEnd = 0; cEnd < n; cEnd++) {
          if (!cellEnd[cEnd]) continue;
          int cSize = cEnd - cStart + 1;
          if (cSize <= 1) { cStart = cEnd + 1; continue; }
          int[] count = new int[cSize];
          for (int ci = 0; ci < cSize; ci++) {
            int v = perm[cStart + ci], cnt = 0;
            for (int nb : neighbors[v]) {
              int nbPos = inv[nb];
              if (nbPos >= sStart && nbPos <= sEnd) cnt++;
            }
            count[ci] = cnt;
          }
          boolean allSame = true;
          for (int ci = 1; ci < cSize; ci++) if (count[ci] != count[0]) { allSame = false; break; }
          if (allSame) { cStart = cEnd + 1; continue; }
          Integer[] idx = new Integer[cSize];
          for (int ci = 0; ci < cSize; ci++) idx[ci] = ci;
          final int fcs = cStart;
          final int[] cnt = count, prm = perm;
          Arrays.sort(idx, (a, bx) -> {
            int c = Integer.compare(cnt[a], cnt[bx]);
            if (c != 0) return c;
            c = Integer.compare(label[prm[fcs + a]], label[prm[fcs + bx]]);
            return c != 0 ? c : Integer.compare(prm[fcs + a], prm[fcs + bx]);
          });
          int[] tmp = new int[cSize];
          for (int ci = 0; ci < cSize; ci++) tmp[ci] = perm[cStart + idx[ci]];
          System.arraycopy(tmp, 0, perm, cStart, cSize);
          for (int ci = 0; ci < cSize; ci++) inv[perm[cStart + ci]] = cStart + ci;
          for (int ci = 0; ci < cSize - 1; ci++) {
            boolean wasBoundary = (cStart + ci == cEnd);
            boolean newBoundary = (cnt[idx[ci]] != cnt[idx[ci + 1]]);
            cellEnd[cStart + ci] = newBoundary || wasBoundary;
            if (newBoundary && !wasBoundary) changed = true;
          }
          cellEnd[cEnd] = true;
          cStart = cEnd + 1;
        }
        sStart = sEnd + 1;
      }
    }
  }

  private static int comparePerm(int[] perm1, int[] perm2, int n, int[] label, int[][] neighbors) {
    int[] inv1 = new int[n], inv2 = new int[n];
    for (int i = 0; i < n; i++) { inv1[perm1[i]] = i; inv2[perm2[i]] = i; }
    for (int pos = 0; pos < n; pos++) {
      int v1 = perm1[pos], v2 = perm2[pos];
      int c = Integer.compare(label[v1], label[v2]);
      if (c != 0) return c;
      int[] nb1 = neighbors[v1], nb2 = neighbors[v2];
      c = Integer.compare(nb1.length, nb2.length);
      if (c != 0) return c;
      int[] cn1 = new int[nb1.length], cn2 = new int[nb2.length];
      for (int k = 0; k < nb1.length; k++) cn1[k] = inv1[nb1[k]];
      for (int k = 0; k < nb2.length; k++) cn2[k] = inv2[nb2[k]];
      Arrays.sort(cn1); Arrays.sort(cn2);
      for (int k = 0; k < nb1.length; k++) {
        c = Integer.compare(cn1[k], cn2[k]);
        if (c != 0) return c;
      }
    }
    return 0;
  }

  private static int partialCompare(
      int[] perm, int[] bestPerm, int n, boolean[] cellEnd, int[] label, int[][] neighbors) {
    int cellStart = 0;
    for (int pos = 0; pos < n; pos++) {
      if (cellEnd[pos]) {
        if (pos == cellStart) {
          int c = Integer.compare(label[perm[pos]], label[bestPerm[pos]]);
          if (c != 0) return c;
        }
        cellStart = pos + 1;
      }
    }
    return 0;
  }

  private static long computeCanonicalHash(MolGraph g, int[] canonLabel) {
    int n = g.n;
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[canonLabel[i]] = i;
    long hash = 0;
    for (int cp = 0; cp < n; cp++) {
      int v = inv[cp];
      hash = hash * HASH_PRIME + g.label[v];
      hash = hash * HASH_PRIME + (long) g.formalCharge[v] + 512L;
      hash = hash * HASH_PRIME + g.massNumber[v];
      // include tetrahedral chirality in hash so enantiomers get different hashes
      hash = hash * HASH_PRIME + (g.tetraChirality != null ? g.tetraChirality[v] : 0);
    }
    for (int cp = 0; cp < n; cp++) {
      int v = inv[cp];
      int[] nbs = g.neighbors[v];
      long[] edgeSigs = new long[nbs.length];
      for (int k = 0; k < nbs.length; k++) {
        int cn = canonLabel[nbs[k]];
        long edgeKey = (long) Math.min(cp, cn) * n + Math.max(cp, cn);
        long sig = edgeKey;
        sig = sig * HASH_PRIME + g.bondOrder(v, nbs[k]);
        sig = sig * 2L + (g.bondInRing(v, nbs[k]) ? 1L : 0L);
        sig = sig * 2L + (g.bondAromatic(v, nbs[k]) ? 1L : 0L);
        // include E/Z double-bond stereo in hash (symmetric: use min canonical position)
        if (cp < cn) {
          int dbConf = g.dbStereo(v, nbs[k]);
          sig = sig * HASH_PRIME + dbConf;
        }
        edgeSigs[k] = sig;
      }
      Arrays.sort(edgeSigs);
      for (long sig : edgeSigs) hash = hash * HASH_PRIME + sig;
    }
    return hash;
  }

  public int[] getOrbits() { ensureCanonical(); return orbit.clone(); }
  public int[] getCanonicalLabeling() { ensureCanonical(); return canonicalLabel.clone(); }
  public long getCanonicalHash() { ensureCanonical(); return canonicalHash; }

  /**
   * Return automorphism generators discovered during canonical labeling.
   * Each generator is a permutation array of length n where gen[i] is the
   * image of atom i.  The full automorphism group is the closure of these
   * generators.
   */
  public int[][] getAutomorphismGenerators() { ensureCanonical(); return autGenerators; }

  /** True when the generator list was capped during canonical search. */
  public boolean automorphismGeneratorsTruncated() { ensureCanonical(); return autGeneratorsTruncated; }

  // ---- Canonical SMILES generation ----

  private static final Map<Integer, String> ORGANIC_SYMBOLS = new HashMap<>();
  static {
    ORGANIC_SYMBOLS.put(5, "B"); ORGANIC_SYMBOLS.put(6, "C"); ORGANIC_SYMBOLS.put(7, "N");
    ORGANIC_SYMBOLS.put(8, "O"); ORGANIC_SYMBOLS.put(15, "P"); ORGANIC_SYMBOLS.put(16, "S");
    ORGANIC_SYMBOLS.put(9, "F"); ORGANIC_SYMBOLS.put(17, "Cl"); ORGANIC_SYMBOLS.put(35, "Br");
    ORGANIC_SYMBOLS.put(53, "I");
  }

  private static final String[] ELEMENT_SYMBOL = new String[119];
  static {
    String[] syms = {
      "?","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S",
      "Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
      "Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
      "Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm",
      "Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os",
      "Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa",
      "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg",
      "Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
    };
    System.arraycopy(syms, 0, ELEMENT_SYMBOL, 0, Math.min(syms.length, ELEMENT_SYMBOL.length));
  }

  public String toCanonicalSmiles() {
    if (n == 0) return "";
    ensureCanonical();
    int[] inv = new int[n];
    for (int i = 0; i < n; i++) inv[canonicalLabel[i]] = i;
    int[] color = new int[n];
    @SuppressWarnings("unchecked")
    List<int[]>[] ringClosures = new List[n];
    Set<Long> ringEdges = new HashSet<>();
    int[] nextRing = {1};
    for (int ci = 0; ci < n; ci++) {
      int start = inv[ci];
      if (color[start] == 0) preScanDfs(start, -1, color, ringClosures, ringEdges, nextRing);
    }
    boolean[] visited = new boolean[n];
    StringBuilder sb = new StringBuilder();
    for (int ci = 0; ci < n; ci++) {
      int start = inv[ci];
      if (visited[start]) continue;
      if (sb.length() > 0) sb.append('.');
      emitSmiles(start, -1, visited, sb, ringClosures, ringEdges);
    }
    return sb.toString();
  }

  private void preScanDfs(int atom, int parent, int[] color, List<int[]>[] rc, Set<Long> ringEdges, int[] nextRing) {
    color[atom] = 1;
    int[] sorted = sortedNeighborsByCanonical(atom);
    int skippedParent = 0;
    for (int nb : sorted) {
      if (nb == parent && skippedParent == 0) { skippedParent = 1; continue; }
      if (color[nb] == 0) {
        preScanDfs(nb, atom, color, rc, ringEdges, nextRing);
      } else if (color[nb] == 1) {
        int rid = nextRing[0]++;
        if (rc[nb] == null) rc[nb] = new ArrayList<>();
        if (rc[atom] == null) rc[atom] = new ArrayList<>();
        rc[nb].add(new int[] {rid, atom});
        rc[atom].add(new int[] {rid, nb});
        ringEdges.add(bondKey(atom, nb));
      }
    }
    color[atom] = 2;
  }

  private void emitSmiles(int atom, int parent, boolean[] visited, StringBuilder sb, List<int[]>[] rc, Set<Long> ringEdges) {
    visited[atom] = true;
    // Collect traversal neighbours in canonical order for tetrahedral parity
    int[] sorted = sortedNeighborsByCanonical(atom);
    writeAtom(atom, parent, sorted, rc, ringEdges, sb);
    if (rc[atom] != null) {
      for (int[] entry : rc[atom]) {
        if (!visited[entry[1]]) writeBondSymbol(atom, entry[1], sb);
        writeRingDigit(entry[0], sb);
      }
    }
    int childCount = 0;
    for (int nb : sorted) if (!visited[nb] && !ringEdges.contains(bondKey(atom, nb))) childCount++;
    int childIdx = 0;
    for (int nb : sorted) {
      if (visited[nb] || ringEdges.contains(bondKey(atom, nb))) continue;
      childIdx++;
      boolean needParen = childCount > 1 && childIdx < childCount;
      if (needParen) sb.append('(');
      writeBondSymbol(atom, nb, sb);
      emitSmiles(nb, atom, visited, sb, rc, ringEdges);
      if (needParen) sb.append(')');
    }
  }

  private int[] sortedNeighborsByCanonical(int atom) {
    int[] sorted = neighbors[atom].clone();
    for (int a = 0; a < sorted.length - 1; a++)
      for (int b2 = a + 1; b2 < sorted.length; b2++)
        if (canonicalLabel[sorted[a]] > canonicalLabel[sorted[b2]]) {
          int tmp = sorted[a]; sorted[a] = sorted[b2]; sorted[b2] = tmp;
        }
    return sorted;
  }

  /**
   * Write an atom token, including tetrahedral stereo (@/@@) when present.
   * The token is always bracket-form when stereo is present.
   *
   * Traversal order for parity: [parent, sorted-children...] — this is the
   * sequence of neighbours as seen when walking from parent. If parent == -1
   * (root atom), the implicit H or first neighbour fills the first position.
   */
  private void writeAtom(int atom, int parent, int[] sortedChildren,
                         List<int[]>[] rc, Set<Long> ringEdges, StringBuilder sb) {
    int z = atomicNum[atom], charge = formalCharge[atom];
    boolean arom = aromatic[atom];
    int chiral = (tetraChirality != null) ? tetraChirality[atom] : 0;
    int mass = (massNumber != null) ? massNumber[atom] : 0;
    // Only emit @/@@  when stereo is defined (1=CW, 2=CCW in our encoding)
    boolean hasStereo = chiral != 0;
    boolean hasIsotope = mass != 0;
    String sym = ORGANIC_SYMBOLS.get(z);
    // Bare organic-subset form is legal only when there is no stereo, no
    // charge, and no isotope label — otherwise we must open brackets.
    if (sym != null && charge == 0 && !hasStereo && !hasIsotope) {
      sb.append(arom ? Character.toLowerCase(sym.charAt(0)) : sym);
    } else {
      sb.append('[');
      if (hasIsotope) sb.append(mass);
      String elemSym = (z >= 0 && z < ELEMENT_SYMBOL.length && ELEMENT_SYMBOL[z] != null) ? ELEMENT_SYMBOL[z] : "?";
      if (arom) {
        sb.append(Character.toLowerCase(elemSym.charAt(0)));
        if (elemSym.length() > 1) sb.append(elemSym.substring(1));
      } else { sb.append(elemSym); }
      if (hasStereo) {
        // Build the neighbour sequence in traversal order:
        // [from (parent or implicit-H), ring-closure partners (by ring-id order), tree children]
        // then determine parity relative to canonical-label-sorted sequence.
        List<Integer> traversalSeq = new ArrayList<>();
        if (parent >= 0) traversalSeq.add(parent);
        // ring closure partners (already emitted, so they come before tree children)
        if (rc[atom] != null) {
          for (int[] entry : rc[atom]) traversalSeq.add(entry[1]);
        }
        for (int nb : sortedChildren) {
          if (!traversalSeq.contains(nb)) traversalSeq.add(nb);
        }
        // chiral==1 → CW (@@), chiral==2 → CCW (@) in our storage convention
        // emit @@ for CLOCKWISE, @ for ANTICLOCKWISE
        sb.append(chiral == 1 ? "@@" : "@");
      }
      // Emit implicit hydrogen count inside brackets. This is required for
      // stereo centres (e.g. [C@@H]) and is standard practice for all bracket
      // atoms — the reader cannot infer H count once brackets are in play.
      int hCount = computeImplicitHForWrite(atom);
      if (hCount > 0) {
        sb.append('H');
        if (hCount > 1) sb.append(hCount);
      }
      if (charge > 0) { sb.append('+'); if (charge > 1) sb.append(charge); }
      else if (charge < 0) { sb.append('-'); if (charge < -1) sb.append(-charge); }
      sb.append(']');
    }
  }

  /**
   * Compute implicit H count for an atom for use in canonical SMILES output.
   * Uses the OpenSMILES multi-valence model via
   * {@link FingerprintEngine#implicitH}. This reuses the same calculation
   * used in fingerprint invariants so canonical SMILES stays consistent
   * with the rest of the stack.
   *
   * <p>Note: Java {@link #bondOrder} returns Kekulised orders (1 or 2)
   * rather than a dedicated aromatic token like C++ uses (4), so no extra
   * "+1 for aromatic π" adjustment is required here — the Kekule sum already
   * captures the aromatic atom's contribution. We still clamp any order-4
   * tokens down to 1 defensively for molecules built directly via
   * {@code MolGraph.Builder}.
   */
  private int computeImplicitHForWrite(int atom) {
    int z = atomicNum[atom];
    int charge = formalCharge[atom];
    int boSum = 0;
    int[] nbs = neighbors[atom];
    for (int k = 0; k < nbs.length; k++) {
      int bo = bondOrder(atom, nbs[k]);
      if (bo == 4) bo = 1;
      boSum += bo;
    }
    return FingerprintEngine.implicitH(z, boSum, charge);
  }

  /**
   * Write the bond token between from→to.
   * For single bonds adjacent to a stereo double bond we emit '/' or '\'.
   * E/Z assignment: conf==1 → TOGETHER (Z, cis) → '/', conf==2 → OPPOSITE (E, trans) → '\'.
   * Direction is relative to the from→to traversal direction.
   */
  private void writeBondSymbol(int from, int to, StringBuilder sb) {
    int order = bondOrder(from, to);
    boolean fromArom = aromatic[from], toArom = aromatic[to];
    // Order 4 = aromatic bond inside an aromatic ring → implicit (no symbol).
    if (order == 4) return;
    if (order == 2) { sb.append('='); return; }
    if (order == 3) { sb.append('#'); return; }
    if (order == 1) {
      // Single bond between TWO aromatic atoms must be emitted explicitly
      // as '-' so that biphenyl c1ccc(-c2ccccc2)cc1 is distinguished from an
      // aromatic ring fusion. Otherwise single bonds are implicit.
      if (fromArom && toArom) { sb.append('-'); return; }
      // Single bond may carry E/Z direction when adjacent to a stereo double bond.
      if (dbStereoConf != null) {
        int dbConf = dbStereo(from, to);
        if (dbConf != 0) {
          // dbConf==1 → TOGETHER → '/', dbConf==2 → OPPOSITE → '\'
          sb.append(dbConf == 1 ? '/' : '\\');
        }
      }
    }
  }

  private void writeRingDigit(int rid, StringBuilder sb) {
    if (rid < 10) sb.append(rid);
    else { sb.append('%'); sb.append(rid); }
  }

  // ---- Builder ----

  public static final class Builder {
    private static final int ROW_ABSENT = 0;
    private static final int ROW_PARALLEL = 1;
    private static final int ROW_DENSE = 2;

    private int n;
    private int[] atomicNum, formalCharge, massNumber;
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
    public Builder ringFlags(boolean[] ring) { this.ring = ring; return this; }
    public Builder aromaticFlags(boolean[] aromatic) { this.aromatic = aromatic; return this; }
    public Builder neighbors(int[][] neighbors) { this.neighbors = neighbors; return this; }
    public Builder bondOrders(int[][] bondOrders) { this.bondOrders = bondOrders; return this; }
    public Builder bondRingFlags(boolean[][] bondRings) { this.bondRings = bondRings; return this; }
    public Builder bondAromaticFlags(boolean[][] bondAroms) { this.bondAroms = bondAroms; return this; }
    public Builder tetrahedralChirality(int[] tetraChirality) { this.tetraChirality = tetraChirality; return this; }
    public Builder doubleBondStereo(int[][] dbStereoConf) { this.dbStereoConf = dbStereoConf; return this; }
    /** Optional external atom IDs (1-based, gaps allowed). Length must equal atomCount if provided. */
    public Builder atomIds(int[] ids) { this.atomIds = ids; return this; }
    public Builder name(String name) { this.name = name; return this; }
    public Builder programLine(String programLine) { this.programLine = programLine; return this; }
    public Builder comment(String comment) { this.comment = comment; return this; }
    public Builder properties(Map<String, String> properties) { this.properties = properties; return this; }

    private static int classifyRow(int[] row, int deg, int atomCount, String name) {
      if (row == null) {
        if (deg == 0) return ROW_PARALLEL;
        throw new IllegalArgumentException(name + " rows must have length == neighbor count or atomCount");
      }
      if (row.length == deg) return ROW_PARALLEL;
      if (row.length == atomCount) return ROW_DENSE;
      throw new IllegalArgumentException(name + " rows must have length == neighbor count or atomCount");
    }

    private static int classifyRow(boolean[] row, int deg, int atomCount, String name) {
      if (row == null) {
        if (deg == 0) return ROW_PARALLEL;
        throw new IllegalArgumentException(name + " rows must have length == neighbor count or atomCount");
      }
      if (row.length == deg) return ROW_PARALLEL;
      if (row.length == atomCount) return ROW_DENSE;
      throw new IllegalArgumentException(name + " rows must have length == neighbor count or atomCount");
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
      if (n < 0)
        throw new IllegalArgumentException("atomCount must be >= 0");
      if (atomicNum == null || atomicNum.length != n)
        throw new IllegalArgumentException("atomicNumbers must have length == atomCount (" + n + ")");
      if (neighbors == null || neighbors.length != n)
        throw new IllegalArgumentException("neighbors must have length == atomCount (" + n + ")");
      // Validate optional arrays: if provided, size must match atomCount
      if (formalCharge != null && formalCharge.length != n)
        throw new IllegalArgumentException("formalCharges length (" + formalCharge.length + ") != atomCount (" + n + ")");
      if (massNumber != null && massNumber.length != n)
        throw new IllegalArgumentException("massNumbers length (" + massNumber.length + ") != atomCount (" + n + ")");
      if (ring != null && ring.length != n)
        throw new IllegalArgumentException("ringFlags length (" + ring.length + ") != atomCount (" + n + ")");
      if (aromatic != null && aromatic.length != n)
        throw new IllegalArgumentException("aromaticFlags length (" + aromatic.length + ") != atomCount (" + n + ")");
      if (tetraChirality != null && tetraChirality.length != n)
        throw new IllegalArgumentException("tetrahedralChirality length (" + tetraChirality.length + ") != atomCount (" + n + ")");
      if (atomIds != null && atomIds.length != n)
        throw new IllegalArgumentException("atomIds length (" + atomIds.length + ") != atomCount (" + n + ")");
      if (bondOrders != null && bondOrders.length != n)
        throw new IllegalArgumentException("bondOrders length (" + bondOrders.length + ") != atomCount (" + n + ")");
      if (bondRings != null && bondRings.length != n)
        throw new IllegalArgumentException("bondRingFlags length (" + bondRings.length + ") != atomCount (" + n + ")");
      if (bondAroms != null && bondAroms.length != n)
        throw new IllegalArgumentException("bondAromaticFlags length (" + bondAroms.length + ") != atomCount (" + n + ")");
      if (dbStereoConf != null && dbStereoConf.length != n)
        throw new IllegalArgumentException("doubleBondStereo length (" + dbStereoConf.length + ") != atomCount (" + n + ")");

      for (int i = 0; i < n; i++) {
        if (atomicNum[i] < 0 || atomicNum[i] > 118)
          throw new IllegalArgumentException("atomicNumbers must be in [0, 118]");
        if (massNumber != null && massNumber[i] < 0)
          throw new IllegalArgumentException("massNumbers must be >= 0");
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
            throw new IllegalArgumentException("neighbor index " + nb + " out of range [0, " + n + ") at atom " + i);
          if (nb == i)
            throw new IllegalArgumentException("self-loop detected in neighbors at atom " + i);
          if (neighborIndex[i].put(nb, k) != null)
            throw new IllegalArgumentException("duplicate neighbor " + nb + " detected at atom " + i);
        }
      }

      int[] bondOrderEnc = new int[n];
      int[] bondRingEnc = new int[n];
      int[] bondAromEnc = new int[n];
      if (bondOrders != null) {
        for (int i = 0; i < n; i++) {
          int deg = neighbors[i] != null ? neighbors[i].length : 0;
          bondOrderEnc[i] = classifyRow(bondOrders[i], deg, n, "bondOrders");
        }
      }
      if (bondRings != null) {
        for (int i = 0; i < n; i++) {
          int deg = neighbors[i] != null ? neighbors[i].length : 0;
          bondRingEnc[i] = classifyRow(bondRings[i], deg, n, "bondRingFlags");
        }
      }
      if (bondAroms != null) {
        for (int i = 0; i < n; i++) {
          int deg = neighbors[i] != null ? neighbors[i].length : 0;
          bondAromEnc[i] = classifyRow(bondAroms[i], deg, n, "bondAromaticFlags");
        }
      }

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
              if (ord != 2)
                throw new IllegalArgumentException("doubleBondStereo requires bond order 2");
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

  static int[] freqMapToSortedArray(Map<Integer, Integer> freq) {
    if (freq.isEmpty()) return EMPTY_NLF;
    int[] arr = new int[freq.size() * 2];
    int k = 0;
    for (Map.Entry<Integer, Integer> e : freq.entrySet()) { arr[k++] = e.getKey(); arr[k++] = e.getValue(); }
    int pairs = freq.size();
    for (int i = 1; i < pairs; i++) {
      int keyLabel = arr[i * 2], keyFreq = arr[i * 2 + 1];
      int j = i - 1;
      while (j >= 0 && arr[j * 2] > keyLabel) {
        arr[(j + 1) * 2] = arr[j * 2];
        arr[(j + 1) * 2 + 1] = arr[j * 2 + 1];
        j--;
      }
      arr[(j + 1) * 2] = keyLabel;
      arr[(j + 1) * 2 + 1] = keyFreq;
    }
    return arr;
  }

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
    // Fast path: avoid HashMap for small neighborhoods (common case: degree 1-6).
    // Collect labels, sort, then run-length encode.
    int[] labels = new int[deg];
    for (int i = 0; i < deg; i++) labels[i] = nlfLabel(g, nbs[i]);
    Arrays.sort(labels);
    // Count distinct labels
    int distinct = 1;
    for (int i = 1; i < deg; i++) if (labels[i] != labels[i - 1]) distinct++;
    int[] arr = new int[distinct * 2];
    int k = 0;
    arr[0] = labels[0]; arr[1] = 1;
    for (int i = 1; i < deg; i++) {
      if (labels[i] == labels[i - 1]) { arr[k + 1]++; }
      else { k += 2; arr[k] = labels[i]; arr[k + 1] = 1; }
    }
    return arr;
  }

  static int[] buildNLF2(MolGraph g, int idx) {
    Map<Integer, Integer> freq = new HashMap<>();
    BitSet direct = g.getAdj()[idx];
    BitSet seen = new BitSet(g.n);
    for (int nb : g.neighbors[idx])
      for (int j : g.neighbors[nb]) {
        if (j == idx || direct.get(j) || seen.get(j)) continue;
        seen.set(j);
        freq.merge(nlfLabel(g, j), 1, Integer::sum);
      }
    return freqMapToSortedArray(freq);
  }

  static int[] buildNLF3(MolGraph g, int idx) {
    Map<Integer, Integer> freq = new HashMap<>();
    BitSet level1 = g.getAdj()[idx];
    BitSet level2 = new BitSet(g.n);
    for (int v = level1.nextSetBit(0); v >= 0; v = level1.nextSetBit(v + 1))
      for (int j : g.neighbors[v]) if (j != idx && !level1.get(j)) level2.set(j);
    BitSet level3 = new BitSet(g.n);
    for (int v = level2.nextSetBit(0); v >= 0; v = level2.nextSetBit(v + 1))
      for (int j : g.neighbors[v]) if (j != idx && !level1.get(j) && !level2.get(j)) level3.set(j);
    for (int j = level3.nextSetBit(0); j >= 0; j = level3.nextSetBit(j + 1))
      freq.merge(nlfLabel(g, j), 1, Integer::sum);
    return freqMapToSortedArray(freq);
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

  static boolean nlfCheckOk(
      int qi, int tj, int[][] qNLF1, int[][] tNLF1, int[][] qNLF2, int[][] tNLF2,
      int[][] qNLF3, int[][] tNLF3, boolean useTwoHop, boolean useThreeHop) {
    if (!nlfOk(qNLF1[qi], tNLF1[tj])) return false;
    if (useTwoHop && !nlfOk(qNLF2[qi], tNLF2[tj])) return false;
    return !useThreeHop || nlfOk(qNLF3[qi], tNLF3[tj]);
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

  // ---- Bond compatibility (ChemOps) ----

  static final class ChemOps {
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
      if (C.useBondStereo) {
        int qConf = g1.dbStereo(qi, qk), tConf = g2.dbStereo(tj, tk);
        if (qConf != 0 && tConf != 0 && qConf != tConf) return false;
      }
      if (C.ringMatchesRingOnly && g1.bondInRing(qi, qk) != g2.bondInRing(tj, tk)) return false;
      if (tautBondRelax) return true;  // bond order relaxed for tautomeric bonds
      if (C.matchBondOrder == ChemOptions.BondOrderMode.ANY) return true;
      if (qOrd == tOrd) return true;
      if (C.matchBondOrder == ChemOptions.BondOrderMode.LOOSE) return true;
      if (C.aromaticityMode == ChemOptions.AromaticityMode.FLEXIBLE) {
        boolean qa = g1.bondAromatic(qi, qk), ta = g2.bondAromatic(tj, tk);
        if ((qa && ta) || (qa && (tOrd == 1 || tOrd == 2)) || (ta && (qOrd == 1 || qOrd == 2))) return true;
      }
      return false;
    }
  }
}

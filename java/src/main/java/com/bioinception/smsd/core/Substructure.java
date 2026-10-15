/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.function.BooleanSupplier;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.isomorphism.VentoFoggia;
import org.openscience.cdk.isomorphism.matchers.Expr;
import org.openscience.cdk.isomorphism.matchers.QueryAtom;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryBond;
import org.openscience.cdk.smarts.Smarts;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smarts.SmartsResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** Molecule and SMARTS substructure queries with bounded first and all-match searches. */
public final class Substructure {

  /** Options are copied when attached to a query. */
  public static final class Options {
    public ChemOptions chemistry = new ChemOptions();
    public boolean induced;
    public long timeoutMillis = -1;
    public long nodeLimit;
    public BooleanSupplier cancellation;
  }

  private final MolGraph moleculeQuery;
  private final IAtomContainer smartsQuery;
  private final Options options;

  private Substructure(MolGraph moleculeQuery, IAtomContainer smartsQuery, Options options) {
    this.moleculeQuery = moleculeQuery;
    this.smartsQuery = smartsQuery;
    this.options = copyOptions(options);
  }

  /** Prepare an independent query snapshot, retaining its atom indices. */
  public static Substructure find(MolGraph query) {
    return new Substructure(MolGraph.copyForSearch(query), null, new Options());
  }

  public static Substructure find(IAtomContainer query) {
    return new Substructure(toGraph(query), null, new Options());
  }

  /** Compile the complete SMARTS expression, including recursive and stereo predicates. */
  public static Substructure smarts(String expression) {
    Objects.requireNonNull(expression, "expression");
    IAtomContainer query = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
    SmartsResult parsed = Smarts.parseToResult(query, expression);
    if (!parsed.ok()) throw new IllegalArgumentException("Invalid SMARTS: " + parsed.getMessage());
    return new Substructure(null, query, new Options());
  }

  public Substructure withOptions(Options options) {
    return new Substructure(moleculeQuery, smartsQuery, Objects.requireNonNull(options, "options"));
  }

  public Substructure withTimeout(long timeoutMillis) {
    Options changed = copyOptions(options);
    changed.timeoutMillis = timeoutMillis;
    return withOptions(changed);
  }

  public SearchResult match(IAtomContainer target) { return run(toGraph(target), 1, false); }
  public SearchResult match(MolGraph target) { return run(target, 1, false); }
  public SearchResult matchAll(IAtomContainer target) { return matchAll(target, 1000); }
  public SearchResult matchAll(MolGraph target) { return matchAll(target, 1000); }
  public SearchResult matchAll(IAtomContainer target, int limit) { return run(toGraph(target), limit, true); }
  public SearchResult matchAll(MolGraph target, int limit) { return run(target, limit, true); }

  /** Search independent targets in input order with this immutable prepared query. */
  public List<SearchResult> matchMany(List<MolGraph> targets) {
    Objects.requireNonNull(targets, "targets");
    List<SearchResult> results = new ArrayList<>();
    for (MolGraph target : targets) results.add(match(target));
    return Collections.unmodifiableList(results);
  }

  /** A retained core, attached groups and the outcome of the substructure search. */
  public static final class Decomposition {
    private final IAtomContainer core;
    private final Map<String,IAtomContainer> rgroups;
    private final SearchResult search;
    private Decomposition(IAtomContainer core,Map<String,IAtomContainer> rgroups,SearchResult search) {
      this.core=core;this.rgroups=Collections.unmodifiableMap(new LinkedHashMap<>(rgroups));this.search=search;
    }
    public IAtomContainer core(){return core;}
    public Map<String,IAtomContainer> rgroups(){return rgroups;}
    public SearchResult search(){return search;}
  }

  /** Split each molecule while preserving no-match, timeout and cancellation outcomes. */
  public static List<Decomposition> decomposeRGroups(
      IAtomContainer core,List<IAtomContainer> molecules,ChemOptions chemistry,long timeoutMillis) {
    Objects.requireNonNull(core,"core");Objects.requireNonNull(molecules,"molecules");
    List<Decomposition> decompositions=new ArrayList<>();
    Options options = new Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    options.timeoutMillis = timeoutMillis;
    Substructure query = find(core).withOptions(options);
    for (IAtomContainer molecule : molecules) {
      MolGraph target = toGraph(molecule);
      SearchResult result = query.match(target);
      Map<String,IAtomContainer> groups=new LinkedHashMap<>();
      IAtomContainer retainedCore=null;
      if (result.found()) {
        Set<Integer> coreAtoms = new HashSet<>(result.matches().get(0).atomMapping().values());
        retainedCore=component(target,coreAtoms);
        Set<Integer> visited = new HashSet<>(coreAtoms);
        int group = 1;
        for (int coreAtom : new TreeSet<>(coreAtoms)) {
          for (int neighbor : target.neighbors[coreAtom]) {
            if (!visited.add(neighbor)) continue;
            Set<Integer> atoms = new LinkedHashSet<>();
            Deque<Integer> pending = new ArrayDeque<>();
            pending.add(neighbor);
            while (!pending.isEmpty()) {
              int atom = pending.remove();
              atoms.add(atom);
              for (int adjacent : target.neighbors[atom]) if (visited.add(adjacent)) pending.add(adjacent);
            }
            groups.put("R"+group++,component(target,atoms));
          }
        }
      }
      decompositions.add(new Decomposition(retainedCore,groups,result));
    }
    return Collections.unmodifiableList(decompositions);
  }

  private static IAtomContainer component(MolGraph graph, Set<Integer> atoms) {
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom : atoms) identity.put(atom, atom);
    return MCS.toMolecule(graph, identity, MCS.allBonds(graph, identity).keySet());
  }

  private static MolGraph toGraph(IAtomContainer molecule) {
    Objects.requireNonNull(molecule, "molecule");
    try { return new MolGraph(molecule.clone()); }
    catch (CloneNotSupportedException failure) { throw new IllegalArgumentException("Cannot copy molecule", failure); }
  }

  private static Options copyOptions(Options source) {
    Options copy = new Options();
    copy.chemistry = ChemOptions.copyOf(Objects.requireNonNull(source.chemistry, "chemistry"));
    copy.induced = source.induced;
    copy.timeoutMillis = source.timeoutMillis;
    copy.nodeLimit = source.nodeLimit;
    copy.cancellation = source.cancellation;
    if (copy.timeoutMillis < -1 || copy.nodeLimit < 0)
      throw new IllegalArgumentException("Timeout must be -1 or non-negative and node limit non-negative");
    return copy;
  }

  private SearchResult run(MolGraph target, int limit, boolean all) {
    Objects.requireNonNull(target, "target");
    if (limit < 1) throw new IllegalArgumentException("Match limit must be positive");
    int queryAtoms = moleculeQuery != null ? moleculeQuery.n : smartsQuery.getAtomCount();
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, queryAtoms, target.n), options.cancellation);
    List<Map<Integer, Integer>> mappings = new ArrayList<>();
    long nodes = 0;
    int searchLimit = all && limit < Integer.MAX_VALUE ? limit + 1 : limit;
    if (!control.stopped() && queryAtoms <= target.n) {
      if (moleculeQuery != null) {
        FrontierMatcher matcher = new FrontierMatcher(moleculeQuery, target, options.chemistry, control);
        matcher.induced = options.induced;
        matcher.nodeLimit = options.nodeLimit;
        if (all) matcher.enumerate(searchLimit, mappings);
        else if (matcher.exists()) {
          Map<Integer, Integer> mapping = new LinkedHashMap<>();
          for (int atom = 0; atom < queryAtoms; atom++) mapping.put(atom, matcher.q2t[atom]);
          mappings.add(mapping);
        }
        nodes = matcher.nodesVisited;
      } else {
        IAtomContainer molecule = MCS.toMolecule(target);
        SmartsPattern.prepare(molecule);
        SmartsCall call = new SmartsCall(control, options.nodeLimit);
        try {
          check(control);
          Pattern pattern = VentoFoggia.findSubstructure(call.query(smartsQuery));
          for (int[] embedding : pattern.matchAll(molecule)) {
            check(control);
            if (options.induced && !induced(smartsQuery, molecule, embedding)) continue;
            Map<Integer, Integer> mapping = new LinkedHashMap<>();
            for (int atom = 0; atom < embedding.length; atom++) mapping.put(atom, embedding[atom]);
            mappings.add(mapping);
            if (!all || mappings.size() >= searchLimit) break;
          }
        } catch (StoppedSearch stopped) {
          // The shared control retains the reason the search stopped.
        }
        nodes = call.nodes;
      }
    }
    boolean truncated = all && mappings.size() > limit;
    if (truncated) mappings.remove(mappings.size() - 1);
    SearchResult.Status status = truncated ? SearchResult.Status.LimitReached : control.status();
    List<SearchResult.Match> matches = new ArrayList<>();
    for (Map<Integer, Integer> mapping : mappings) matches.add(describe(mapping));
    boolean complete = status == SearchResult.Status.Completed;
    return new SearchResult(matches, status, complete || !matches.isEmpty(), complete,
        matches.isEmpty() ? 0 : queryAtoms, nodes, control.elapsedMillis());
  }

  private SearchResult.Match describe(Map<Integer, Integer> mapping) {
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    Map<Integer, Set<Integer>> adjacency = new LinkedHashMap<>();
    for (int atom : mapping.keySet()) adjacency.put(atom, new LinkedHashSet<>());
    for (int atom : mapping.keySet()) {
      Iterable<Integer> neighbors;
      if (moleculeQuery != null) {
        List<Integer> values = new ArrayList<>();
        for (int neighbor : moleculeQuery.neighbors[atom]) values.add(neighbor);
        neighbors = values;
      } else {
        List<Integer> values = new ArrayList<>();
        for (IAtom neighbor : smartsQuery.getConnectedAtomsList(smartsQuery.getAtom(atom)))
          values.add(smartsQuery.indexOf(neighbor));
        neighbors = values;
      }
      for (int neighbor : neighbors) if (neighbor > atom) {
        bonds.put(new SearchResult.Bond(atom, neighbor),
            new SearchResult.Bond(mapping.get(atom), mapping.get(neighbor)));
        adjacency.get(atom).add(neighbor);
        adjacency.get(neighbor).add(atom);
      }
    }
    List<Set<Integer>> fragments = new ArrayList<>();
    Set<Integer> visited = new HashSet<>();
    for (int atom : mapping.keySet()) if (visited.add(atom)) {
      Set<Integer> fragment = new LinkedHashSet<>();
      Deque<Integer> pending = new ArrayDeque<>();
      pending.add(atom);
      while (!pending.isEmpty()) {
        int current = pending.remove();
        fragment.add(current);
        for (int neighbor : adjacency.get(current)) if (visited.add(neighbor)) pending.add(neighbor);
      }
      fragments.add(fragment);
    }
    return new SearchResult.Match(mapping, bonds, fragments, mapping.size());
  }

  private static boolean induced(IAtomContainer query, IAtomContainer target, int[] mapping) {
    for (int first = 0; first < mapping.length; first++)
      for (int second = first + 1; second < mapping.length; second++)
        if ((query.getBond(query.getAtom(first), query.getAtom(second)) != null)
            != (target.getBond(target.getAtom(mapping[first]), target.getAtom(mapping[second])) != null)) return false;
    return true;
  }

  private static final class StoppedSearch extends RuntimeException {}

  private static void check(SearchResult.Control control) {
    if (control.stopped()) throw new StoppedSearch();
  }

  /** Per-call query predicates retain CDK recursion, stereo and component semantics. */
  private static final class SmartsCall {
    final SearchResult.Control control;
    final long nodeLimit;
    long nodes;

    SmartsCall(SearchResult.Control control, long nodeLimit) {
      this.control = control;
      this.nodeLimit = nodeLimit;
    }

    IAtomContainer query(IAtomContainer original) {
      IAtomContainer copy = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
      copy.setProperties(new HashMap<>(original.getProperties()));
      Map<IChemObject, IChemObject> objects = new IdentityHashMap<>();
      for (IAtom atom : original.atoms()) {
        check(control);
        QueryAtom source = (QueryAtom) AtomRef.deref(atom);
        QueryAtom bounded = new QueryAtom(expression(source.getExpression())) {
          @Override public boolean matches(IAtom target) {
            check(control);
            if (nodeLimit > 0 && nodes >= nodeLimit) { control.limitReached(); throw new StoppedSearch(); }
            nodes++;
            return super.matches(target);
          }
        };
        bounded.setProperties(new HashMap<>(atom.getProperties()));
        bounded.setFlags(atom.getFlags());
        copy.addAtom(bounded);
        objects.put(atom, copy.getAtom(copy.getAtomCount() - 1));
        objects.put(source, copy.getAtom(copy.getAtomCount() - 1));
      }
      for (IBond bond : original.bonds()) {
        check(control);
        QueryBond source = (QueryBond) BondRef.deref(bond);
        QueryBond bounded = new QueryBond((IAtom) objects.get(bond.getBegin()),
            (IAtom) objects.get(bond.getEnd()), expression(source.getExpression()));
        bounded.setProperties(new HashMap<>(bond.getProperties()));
        bounded.setFlags(bond.getFlags());
        copy.addBond(bounded);
        objects.put(bond, copy.getBond(copy.getBondCount() - 1));
        objects.put(source, copy.getBond(copy.getBondCount() - 1));
      }
      for (IStereoElement<?, ?> stereo : original.stereoElements()) copy.addStereoElement(stereo.mapStrict(objects));
      return copy;
    }

    Expr expression(Expr original) {
      switch (original.type()) {
        case RECURSIVE:
          return new Expr(Expr.Type.RECURSIVE, query(original.subquery()));
        case AND:
        case OR:
        case NOT:
          return new Expr(original.type(), expression(original.left()),
              original.right() == null ? null : expression(original.right()));
        default:
          return original;
      }
    }
  }

  static Matcher makeMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchResult.Control tb) {
    return new FrontierMatcher(gq, gt, C, tb);
  }

  /** Return a mapped atom whose tetrahedral or double-bond ligand relation is incompatible. */
  static int stereoMismatchAtom(MolGraph query, MolGraph target, int[] mapping, ChemOptions chemistry) {
    return stereoMismatchAtom(query,target,mapping,chemistry,null);
  }

  static int stereoMismatchAtom(MolGraph query,MolGraph target,int[] mapping,ChemOptions chemistry,
                                Set<SearchResult.Bond> selectedBonds) {
    if (chemistry.useChirality && query.tetraChirality != null && target.tetraChirality != null) {
      for (int atom = 0; atom < query.n; atom++) {
        int mapped = mapping[atom];
        if (mapped < 0 || mapped >= target.n || query.tetraChirality[atom] == 0
            || target.tetraChirality[mapped] == 0) continue;
        int[] queryOrder = query.tetrahedralLigands(atom), targetOrder = target.tetrahedralLigands(mapped);
        if (queryOrder == null || targetOrder == null || queryOrder.length != 4 || targetOrder.length != 4) continue;
        int[] permutation = new int[4];
        Arrays.fill(permutation, -1);
        boolean[] used = new boolean[4];
        int mappedCount = 0;
        for (int slot = 0; slot < 4; slot++) {
          int ligand = queryOrder[slot];
          if (ligand < 0 || mapping[ligand] < 0
              || (selectedBonds!=null&&!selectedBonds.contains(new SearchResult.Bond(atom,ligand)))) continue;
          int mappedLigand = mapping[ligand];
          for (int targetSlot = 0; targetSlot < 4; targetSlot++) {
            if (targetOrder[targetSlot] == mappedLigand) {
              permutation[slot] = targetSlot;
              used[targetSlot] = true;
              mappedCount++;
              break;
            }
          }
          if (permutation[slot] < 0) return atom;
        }
        if (mappedCount < 3) continue;
        for (int slot = 0; slot < 4; slot++) {
          if (permutation[slot] >= 0) continue;
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
            if (permutation[first] > permutation[second]) inversions++;
        // Builder tags use SMILES winding; CDK tags use its ligand order.
        int queryTag = query.mol == null ? 3 - query.tetraChirality[atom] : query.tetraChirality[atom];
        int targetTag = target.mol == null ? 3 - target.tetraChirality[mapped] : target.tetraChirality[mapped];
        if ((queryTag == targetTag) != ((inversions & 1) == 0)) return atom;
      }
    }
    if (chemistry.useBondStereo) {
      for (int atom = 0; atom < query.n; atom++) {
        if (mapping[atom] < 0) continue;
        for (int other : query.neighbors[atom])
          if (other > atom && (selectedBonds==null||selectedBonds.contains(new SearchResult.Bond(atom,other)))
              && !doubleBondParityCompatible(query,atom,other,target,mapping,selectedBonds)) return atom;
      }
    }
    return -1;
  }

  /** Compare the stored double-bond relation against whichever query ligands are mapped. */
  private static boolean doubleBondParityCompatible(MolGraph query, int first, int second,
                                                    MolGraph target, int[] mapping) {
    return doubleBondParityCompatible(query,first,second,target,mapping,null);
  }

  private static boolean doubleBondParityCompatible(MolGraph query,int first,int second,MolGraph target,
                                                    int[] mapping,Set<SearchResult.Bond> selectedBonds) {
    int mappedFirst = mapping[first], mappedSecond = mapping[second];
    if (mappedFirst < 0 || mappedSecond < 0) return true;
    int queryConfiguration = query.dbStereo(first, second);
    int targetConfiguration = target.dbStereo(mappedFirst, mappedSecond);
    if (queryConfiguration == 0 || targetConfiguration == 0) return true;
    for (int end = 0; end < 2; end++) {
      int queryEnd = end == 0 ? first : second, otherQueryEnd = end == 0 ? second : first;
      int targetEnd = mapping[queryEnd], otherTargetEnd = mapping[otherQueryEnd];
      int queryReference = query.doubleBondReference(queryEnd, otherQueryEnd);
      int targetReference = target.doubleBondReference(targetEnd, otherTargetEnd);
      int queryLigand = -1;
      for (int neighbor : query.neighbors[queryEnd])
        if (neighbor != otherQueryEnd && mapping[neighbor] >= 0
            && (selectedBonds==null||selectedBonds.contains(new SearchResult.Bond(queryEnd,neighbor)))) {
          queryLigand = neighbor; break;
        }
      if (queryLigand < 0 || queryReference < 0 || targetReference < 0) return true;
      if (!target.hasBond(targetEnd, mapping[queryLigand])) return false;
      if ((queryLigand != queryReference) != (mapping[queryLigand] != targetReference))
        queryConfiguration = 3 - queryConfiguration;
    }
    return queryConfiguration == targetConfiguration;
  }

  interface Matcher {
    boolean exists();
    void enumerate(int maxSolutions, List<Map<Integer, Integer>> out);
  }

  @SuppressWarnings("unchecked")
  abstract static class MoleculeMatcher implements Matcher {
    final MolGraph gq, gt;
    final ChemOptions C;
    final SearchResult.Control tb;
    final int Nq, Nt;
    final int[] q2t, t2q;
    final int[][] qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3;
    final boolean useNLF, useTwoHop, useThreeHop, useBitParallel, useRingOnly, useStereo;
    final boolean bitParallelSufficient;
    final long[][] domain;
    final long[] usedMask;
    final int tWords;
    final boolean singleWord;
    long nodesVisited = 0, backtracks = 0, candidatesTried = 0;
    long prunesAtom = 0, prunesBond = 0, prunesDegree = 0, prunesNLF = 0;
    boolean timedOut = false, found = false;
    boolean induced;
    long nodeLimit;
    boolean pivotEnabled = false; // Mapped-neighbour pivot: only enabled during exists()
    int[] qdeg, tdeg;
    final int[][] qNeighborsByDegDesc;
    // Pre-allocated buffers to avoid per-call allocations in selectCandidates/iterateBits
    final long[] availBuf;
    final int[] candBuf;
    // Per-depth candidate buffer stack: avoids Arrays.copyOf on every selectCandidates call
    final int[][] candBufByDepth;

    MoleculeMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchResult.Control tb) {
      this.gq = gq; this.gt = gt;
      this.qdeg = gq.degree; this.tdeg = gt.degree;
      this.C = C; this.tb = tb;
      this.Nq = gq.n; this.Nt = gt.n;
      this.tWords = (Nt + 63) >>> 6;
      this.singleWord = (tWords <= 1);
      this.q2t = new int[Nq]; this.t2q = new int[Nt];
      Arrays.fill(q2t, -1); Arrays.fill(t2q, -1);
      this.useNLF = C.matchAtomType && !C.tautomerAware;
      this.useTwoHop = useNLF && C.useTwoHopNLF && Nq > 12 && Nt > 12;
      this.useThreeHop = useNLF && C.useThreeHopNLF && Nq > 20 && Nt > 20;
      this.useBitParallel = C.useBitParallelFeasibility;
      this.useRingOnly = C.ringMatchesRingOnly; this.useStereo = C.useBondStereo;
      this.bitParallelSufficient = useBitParallel
          && C.matchBondOrder == ChemOptions.BondOrderMode.ANY && !useStereo && !useRingOnly
          && C.aromaticityMode != ChemOptions.AromaticityMode.STRICT;
      // Reuse the graph's cached neighbours sorted by degree.
      this.qNeighborsByDegDesc = gq.getNeighborsByDegDesc();
      // Use cached NLF tables from MolGraph (lazy-built, reused across calls on the same graph).
      this.qNLF1 = useNLF ? matchingNLF(gq.getNLF1(), C) : null;
      this.tNLF1 = useNLF ? matchingNLF(gt.getNLF1(), C) : null;
      this.qNLF2 = useTwoHop ? matchingNLF(gq.getNLF2(), C) : null;
      this.tNLF2 = useTwoHop ? matchingNLF(gt.getNLF2(), C) : null;
      this.qNLF3 = useThreeHop ? matchingNLF(gq.getNLF3(), C) : null;
      this.tNLF3 = useThreeHop ? matchingNLF(gt.getNLF3(), C) : null;
      // Ensure expensive lazy fields are computed only when the options actually need them.
      if (C.tautomerAware) { gq.ensureTautomerClasses(); gt.ensureTautomerClasses(); }
      if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) { gq.ensureRingCounts(); gt.ensureRingCounts(); }
      this.domain = new long[Nq][tWords];
      if (Nq > 0 && Nt > 50) {
        // Group target atoms by label before constructing the domains.
        // Pass 1: compute label keys and find max key for sizing
        int[] tLabelKeys = new int[Nt];
        int maxKey = 0;
        for (int j = 0; j < Nt; j++) {
          int key2 = C.matchAtomType ? gt.atomicNum[j] : 0;
          int ringBit = C.ringMatchesRingOnly ? (gt.ring[j] ? 1 : 0) : 0;
          key2 = (key2 << 2) | (gt.aromatic[j] ? 2 : 0) | ringBit;
          tLabelKeys[j] = key2;
          if (key2 > maxKey) maxKey = key2;
        }
        // Count entries per label key
        int[] bucketCount = new int[maxKey + 1];
        for (int j = 0; j < Nt; j++) bucketCount[tLabelKeys[j]]++;
        // Prefix-sum to get bucket start offsets
        int[] bucketStart = new int[maxKey + 2]; // +2 so bucketStart[maxKey+1] = Nt
        bucketStart[0] = 0;
        for (int k = 0; k <= maxKey; k++) bucketStart[k + 1] = bucketStart[k] + bucketCount[k];
        // Pass 2: fill flat data array (scatter using running offsets)
        int[] flatBucket = new int[Nt];
        int[] writePos = new int[maxKey + 1];
        System.arraycopy(bucketStart, 0, writePos, 0, maxKey + 1);
        for (int j = 0; j < Nt; j++) {
          int k = tLabelKeys[j];
          flatBucket[writePos[k]++] = j;
        }
        // Iterate buckets for domain init
        for (int i = 0; i < Nq; i++) {
          if (tb.stopped()) break;
          for (int k = 0; k <= maxKey; k++) {
            int bStart = bucketStart[k], bEnd = bucketStart[k + 1];
            if (bStart == bEnd) continue; // empty bucket
            if (C.matchAtomType && !C.tautomerAware
                && gq.atomicNum[i] != gt.atomicNum[flatBucket[bStart]]) { prunesAtom += (bEnd - bStart); continue; }
            for (int b = bStart; b < bEnd; b++) {
              int j = flatBucket[b];
              if (atomsCompatFast(gq, i, gt, j, C)) domain[i][j >>> 6] |= 1L << (j & 63);
              else prunesAtom++;
            }
          }
        }
      } else {
        for (int i = 0; i < Nq && !tb.stopped(); i++)
          for (int j = 0; j < Nt; j++) {
            if (atomsCompatFast(gq, i, gt, j, C)) domain[i][j >>> 6] |= 1L << (j & 63);
            else prunesAtom++;
          }
      }
      this.usedMask = new long[tWords];
      this.availBuf = new long[tWords];
      this.candBuf = new int[Nt];
      this.candBufByDepth = new int[Nq][];
    }

    /** Project cached labels onto the constraints enforced by this search. */
    static int[][] matchingNLF(int[][] tables, ChemOptions C) {
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT) return tables;
      int[][] result = new int[tables.length][];
      for (int i = 0; i < tables.length; i++) {
        int[] source = tables[i], merged = new int[source.length];
        int size = 0;
        for (int j = 0; j < source.length; j += 2) {
          int element = source[j] >>> 1;
          if (size > 0 && merged[size - 2] == element) merged[size - 1] += source[j + 1];
          else { merged[size++] = element; merged[size++] = source[j + 1]; }
        }
        result[i] = Arrays.copyOf(merged, size);
      }
      return result;
    }

    int[] candidateBuffer(int depth) {
      if (candBufByDepth[depth] == null) candBufByDepth[depth] = new int[Nt];
      return candBufByDepth[depth];
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

    static int[] iterateBits(long[] words, int maxBits) {
      int count = 0; for (long w : words) count += Long.bitCount(w);
      int[] result = new int[count]; int idx = 0;
      for (int w = 0; w < words.length; w++) {
        long bits = words[w];
        while (bits != 0) { result[idx++] = (w << 6) | Long.numberOfTrailingZeros(bits); bits &= bits - 1; }
      }
      return result;
    }

    /** Write set bits into pre-allocated buffer, return count. Zero-allocation hot path. */
    int iterateBitsInto(long[] words, int[] buf) {
      int idx = 0;
      for (int w = 0; w < words.length; w++) {
        long bits = words[w];
        while (bits != 0 && idx < buf.length) { buf[idx++] = (w << 6) | Long.numberOfTrailingZeros(bits); bits &= bits - 1; }
      }
      return idx;
    }

    int[] domainCandidates(int qi) {
      for (int w = 0; w < tWords; w++) availBuf[w] = domain[qi][w] & ~usedMask[w];
      int len = iterateBitsInto(availBuf, candBuf);
      return Arrays.copyOf(candBuf, len);
    }

    boolean greedyProbe(int[] order) {
      int matched = 0;
      for (int qi : order) {
        if (stopped()) return false;
        nodesVisited++;
        int bestTj = -1, bestDist = Integer.MAX_VALUE; boolean found2 = false;
        for (int tj : domainCandidates(qi)) {
          int dist = Math.abs(tdeg[tj] - qdeg[qi]);
          if (dist < bestDist && feasible(qi, tj)) { bestTj = tj; bestDist = dist; found2 = true; if (dist == 0) break; }
        }
        if (!found2) {
          for (int k = matched - 1; k >= 0; k--) {
            int qk = order[k], tk = q2t[qk]; q2t[qk] = -1; t2q[tk] = -1;
            usedMask[tk >>> 6] &= ~(1L << (tk & 63));
          }
          return false;
        }
        q2t[qi] = bestTj; t2q[bestTj] = qi; usedMask[bestTj >>> 6] |= 1L << (bestTj & 63); matched++;
      }
      if (stereoMismatchAtom(gq, gt, q2t, C) < 0) return true;
      for (int qi : order) {
        int tj = q2t[qi];
        q2t[qi] = -1;
        t2q[tj] = -1;
        usedMask[tj >>> 6] &= ~(1L << (tj & 63));
      }
      return false;
    }

    int[] fastisoOrder() {
      int[] order = new int[Nq]; boolean[] picked = new boolean[Nq];
      long[][] workDomain = new long[Nq][tWords];
      for (int i = 0; i < Nq; i++) System.arraycopy(domain[i], 0, workDomain[i], 0, tWords);
      long[] reachable = new long[tWords];
      for (int pos = 0; pos < Nq; pos++) {
        if (tb.stopped()) return order;
        int bestQ = -1, bestSize = Integer.MAX_VALUE, bestDeg = -1;
        for (int i = 0; i < Nq; i++) {
          if (picked[i]) continue;
          int size = 0; for (int w = 0; w < tWords; w++) size += Long.bitCount(workDomain[i][w]);
          if (size < bestSize || (size == bestSize && gq.degree[i] > bestDeg)) { bestQ = i; bestSize = size; bestDeg = gq.degree[i]; }
        }
        order[pos] = bestQ; picked[bestQ] = true;
        for (int qk : gq.neighbors[bestQ]) {
          if (picked[qk]) continue;
          // Skip neighbours whose domain is already empty.
          boolean hasAny = false;
          for (int w = 0; w < tWords && !hasAny; w++) if (workDomain[qk][w] != 0) hasAny = true;
          if (!hasAny) continue;
          Arrays.fill(reachable, 0, tWords, 0L);
          for (int w = 0; w < tWords; w++) {
            long bits = workDomain[bestQ][w];
            while (bits != 0) {
              int bit = Long.numberOfTrailingZeros(bits); int tj = (w << 6) | bit;
              if (tj < Nt) for (int rw = 0; rw < tWords; rw++) reachable[rw] |= gt.adjLong[tj][rw];
              bits &= bits - 1;
            }
          }
          for (int rw = 0; rw < tWords; rw++) workDomain[qk][rw] &= reachable[rw];
        }
      }
      return order;
    }

    void sortByDegreeProximity(int[] cands, int qi) {
      for (int i = 1; i < cands.length; i++) {
        int key = cands[i], keyDist = Math.abs(tdeg[key] - qdeg[qi]); int j = i - 1;
        while (j >= 0 && Math.abs(tdeg[cands[j]] - qdeg[qi]) > keyDist) { cands[j + 1] = cands[j]; j--; }
        cands[j + 1] = key;
      }
    }

    /** Sort first {@code len} entries of {@code buf} by degree proximity to query atom qi. */
    void sortByDegreeProximityBuf(int[] buf, int len, int qi) {
      for (int i = 1; i < len; i++) {
        int key = buf[i], keyDist = Math.abs(tdeg[key] - qdeg[qi]); int j = i - 1;
        while (j >= 0 && Math.abs(tdeg[buf[j]] - qdeg[qi]) > keyDist) { buf[j + 1] = buf[j]; j--; }
        buf[j + 1] = key;
      }
    }

    /** Mapped-neighbour pivot: prefer the candidate with most mapped-neighbour overlap. */
    void applyPivotHeuristic(int[] buf, int len, int qi) {
      if (len <= 1) return;
      // A pivot is useful with at least two mapped neighbours.
      int mappedNbrs = 0;
      for (int qk : gq.neighbors[qi]) if (q2t[qk] != -1) mappedNbrs++;
      if (mappedNbrs < 2) return;
      // Find the candidate with most adjacency overlap to mapped neighbours.
      int bestPivot = 0, bestOverlap = 0;
      for (int c = 0; c < len; c++) {
        int tj = buf[c]; int overlap = 0;
        for (int qk : gq.neighbors[qi]) {
          int tk = q2t[qk];
          if (tk != -1 && (gt.adjLong[tj][tk >>> 6] & (1L << (tk & 63))) != 0)
            overlap++;
        }
        if (overlap > bestOverlap) { bestOverlap = overlap; bestPivot = c; }
      }
      if (bestPivot != 0) { int tmp = buf[0]; buf[0] = buf[bestPivot]; buf[bestPivot] = tmp; }
    }

    boolean feasible(int qi, int tj) {
      if (gt.degree[tj] < gq.degree[qi]) { prunesDegree++; return false; }
      // Ring-only check (cheapest boolean test — run early for fast exit)
      if (useRingOnly && gq.ring[qi] != gt.ring[tj]) { prunesAtom++; return false; }
      if (useNLF && !MolGraph.nlfOk(qNLF1[qi], tNLF1[tj])) { prunesNLF++; return false; }
      if (useTwoHop && !MolGraph.nlfOk(qNLF2[qi], tNLF2[tj])) { prunesNLF++; return false; }
      if (useThreeHop && !MolGraph.nlfOk(qNLF3[qi], tNLF3[tj])) { prunesNLF++; return false; }
      if (induced) {
        for (int mappedQuery = 0; mappedQuery < Nq; mappedQuery++)
          if (q2t[mappedQuery] >= 0 && !gq.hasBond(qi, mappedQuery)
              && gt.hasBond(tj, q2t[mappedQuery])) { prunesBond++; return false; }
      }
      int[] nbSorted = qNeighborsByDegDesc[qi];
      if (useBitParallel) {
        long[] tAdj = gt.adjLong[tj];
        for (int qk : nbSorted) {
          int tk = q2t[qk];
          if (tk != -1 && (tAdj[tk >>> 6] & (1L << (tk & 63))) == 0) { prunesBond++; return false; }
        }
        if (bitParallelSufficient) return true;
      }
      for (int qk : nbSorted) {
        int tk = q2t[qk];
        if (tk != -1 && !MolGraph.ChemOps.bondsCompatible(gq, qi, qk, gt, tj, tk, C)) { prunesBond++; return false; }
      }
      return true;
    }

    int[] vf3LightOrder() {
      Map<Integer, Integer> labelFreq = new HashMap<>();
      for (int i = 0; i < Nq; i++) labelFreq.merge(gq.label[i], 1, Integer::sum);
      Integer[] order = new Integer[Nq];
      for (int i = 0; i < Nq; i++) order[i] = i;
      Arrays.sort(order, (a, b) -> {
        int freqA = labelFreq.getOrDefault(gq.label[a], 0), freqB = labelFreq.getOrDefault(gq.label[b], 0);
        if (freqA != freqB) return Integer.compare(freqA, freqB);
        if (gq.degree[a] != gq.degree[b]) return Integer.compare(gq.degree[b], gq.degree[a]);
        int domA = 0, domB = 0;
        for (int w = 0; w < tWords; w++) { domA += Long.bitCount(domain[a][w]); domB += Long.bitCount(domain[b][w]); }
        return Integer.compare(domA, domB);
      });
      int[] result = new int[Nq]; for (int i = 0; i < Nq; i++) result[i] = order[i]; return result;
    }

    /** Write candidates for query atom qi into candBufByDepth[pos] and return the count. */
    abstract int selectCandidates(int qi, int pos);
    void onMatch(int pos, int tj) {}
    void onUnmatch(int pos, int tj) {}

    boolean stopped() {
      if (nodeLimit > 0 && nodesVisited >= nodeLimit) tb.limitReached();
      if (tb.stopped()) { timedOut = true; return true; }
      return false;
    }

    /** Choose the fewest feasible targets after enforcing the current mapped neighbours. */
    int selectNext(int[] order, int position) {
      boolean hasFrontier = false;
      for (int index = position; index < order.length && !hasFrontier; index++) {
        if (tb.stopped()) { timedOut = true; return 0; }
        hasFrontier = hasMappedNeighbor(order[index]);
      }
      int bestIndex = position, bestCount = Integer.MAX_VALUE;
      for (int index = position; index < order.length; index++) {
        if (tb.stopped()) { timedOut = true; return 0; }
        int queryAtom = order[index];
        if (hasFrontier && !hasMappedNeighbor(queryAtom)) continue;
        int candidates = selectCandidates(queryAtom, position);
        int[] buffer = candBufByDepth[position];
        int compatible = 0;
        for (int candidate = 0; candidate < candidates; candidate++) {
          if (tb.stopped()) { timedOut = true; return 0; }
          int targetAtom = buffer[candidate];
          if (feasible(queryAtom, targetAtom)) buffer[compatible++] = targetAtom;
        }
        if (compatible < bestCount) {
          bestIndex = index;
          bestCount = compatible;
          System.arraycopy(buffer, 0, candBuf, 0, compatible);
          if (compatible == 0) break;
        }
      }
      int swap = order[position]; order[position] = order[bestIndex]; order[bestIndex] = swap;
      System.arraycopy(candBuf, 0, candidateBuffer(position), 0, bestCount);
      return bestCount;
    }

    boolean hasMappedNeighbor(int queryAtom) {
      for (int neighbor : gq.neighbors[queryAtom]) if (q2t[neighbor] >= 0) return true;
      return false;
    }

    void backtrack(int[] order, int pos) {
      if (stopped()) return;
      nodesVisited++; if (found) return;
      if (pos == order.length) { found = stereoMismatchAtom(gq, gt, q2t, C) < 0; return; }
      int candCount = selectNext(order, pos);
      int qi = order[pos];
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        if (tb.stopped()) { timedOut = true; return; }
        candidatesTried++;
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63); onMatch(pos, tj);
        backtrack(order, pos + 1);
        if (found || timedOut) return;
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63)); onUnmatch(pos, tj); backtracks++;
      }
    }

    void enumerateRec(int[] order, int pos, List<Map<Integer, Integer>> out, int maxSolutions) {
      if (stopped()) return;
      nodesVisited++; if (out.size() >= maxSolutions) return;
      if (pos == order.length) {
        if (stereoMismatchAtom(gq, gt, q2t, C) >= 0) return;
        Map<Integer, Integer> map = new LinkedHashMap<>(); for (int k : order) map.put(k, q2t[k]); out.add(map); return;
      }
      int candCount = selectNext(order, pos);
      int qi = order[pos];
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        if (tb.stopped()) { timedOut = true; return; }
        candidatesTried++;
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63); onMatch(pos, tj);
        enumerateRec(order, pos + 1, out, maxSolutions);
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63)); onUnmatch(pos, tj); backtracks++;
        if (out.size() >= maxSolutions || timedOut) return;
      }
    }
  }

  static final class FrontierMatcher extends MoleculeMatcher {
    FrontierMatcher(MolGraph gq, MolGraph gt, ChemOptions C, SearchResult.Control tb) {
      super(gq, gt, C, tb);
    }
    public boolean exists() {
      if (tb.stopped()) return false;
      if (Nq == 0) return true;
      pivotEnabled = true;
      int[] order = fastisoOrder();
      if (tb.stopped()) return false;
      if (Nq <= 50 && greedyProbe(order)) { found = true; pivotEnabled = false; return true; }
      backtrack(order, 0); pivotEnabled = false; return found;
    }
    public void enumerate(int maxSolutions, List<Map<Integer, Integer>> out) {
      if (tb.stopped()) return;
      if (Nq == 0) { out.add(new LinkedHashMap<>()); return; }
      enumerateRec(fastisoOrder(), 0, out, maxSolutions);
    }
    int selectCandidates(int qi, int pos) {
      int[] buf = candidateBuffer(pos);
      for (int w = 0; w < tWords; w++) availBuf[w] = domain[qi][w] & ~usedMask[w];
      for (int queryNeighbor : gq.neighbors[qi]) {
        int mappedTarget = q2t[queryNeighbor];
        if (mappedTarget < 0) continue;
        for (int word = 0; word < tWords; word++) availBuf[word] &= gt.adjLong[mappedTarget][word];
      }
      int len = iterateBitsInto(availBuf, buf);
      sortByDegreeProximityBuf(buf, len, qi);
      if (pivotEnabled) applyPivotHeuristic(buf, len, qi);
      return len;
    }

  }
}

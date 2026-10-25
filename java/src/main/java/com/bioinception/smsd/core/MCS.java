/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.function.BooleanSupplier;
import org.openscience.cdk.interfaces.IAtomContainer;

/** Exact maximum common induced subgraph and maximum common edge subgraph search. */
public final class MCS {
  public enum Objective { MCIS, MCES }
  public enum Connectivity { Connected, Fragments }

  /** Chemistry, objective and resource limits copied when preparing a query. */
  public static final class Options {
    public ChemOptions chemistry = new ChemOptions();
    public Objective objective = Objective.MCIS;
    public Connectivity connectivity = Connectivity.Connected;
    public double[] atomWeights;
    public int minFragmentAtoms = 1;
    public int maxFragments = Integer.MAX_VALUE;
    public long timeoutMillis = -1;
    public long nodeLimit;
    public BooleanSupplier cancellation;
    public Set<Integer> excludedTargetAtoms = new HashSet<>();
  }

  private final MolGraph query;
  private final Options options;

  private MCS(MolGraph query, Options options) {
    this.query = Objects.requireNonNull(query, "query");
    this.options = copy(options, query.n);
  }

  /** Prepare a snapshot; subsequent changes to the caller's graph do not affect searches. */
  public static MCS find(MolGraph query) {
    return new MCS(MolGraph.copyForSearch(query), new Options());
  }
  public static MCS find(IAtomContainer query) {
    return new MCS(MolGraph.copyOf(query), new Options());
  }
  public MCS withOptions(Options options) {
    return new MCS(query, options);
  }
  public MCS withTimeout(long millis) {
    Options changed = copy(options, query.n);
    changed.timeoutMillis = millis;
    return new MCS(query, changed);
  }
  public SearchResult match(MolGraph target) {
    return run(target, 1, false);
  }
  public SearchResult match(IAtomContainer target) {
    return match(MolGraph.copyOf(target));
  }
  public SearchResult matchAll(MolGraph target) {
    return matchAll(target, 1000);
  }
  public SearchResult matchAll(IAtomContainer target) {
    return matchAll(MolGraph.copyOf(target), 1000);
  }
  public SearchResult matchAll(IAtomContainer target, int limit) {
    return matchAll(MolGraph.copyOf(target), limit);
  }
  public SearchResult matchAll(MolGraph target, int limit) {
    return run(target, limit, true);
  }

  private SearchResult run(MolGraph target, int limit, boolean all) {
    Objects.requireNonNull(target, "target");
    if (limit < 1)
      throw new IllegalArgumentException("Match limit must be positive");
    for (int targetAtom : options.excludedTargetAtoms)
      if (targetAtom >= target.n)
        throw new IllegalArgumentException("Excluded target atom index out of range");
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, query.n, target.n, true),
        options.cancellation);
    if (control.stopped())
      return new SearchResult(
          Collections.emptyList(), control.status(), false, false, 0.0, 0, control.elapsedMillis());
    return search(query, target, options, limit, all, control);
  }

  /**
   * Searches {@code query} in {@code target} as {@link #match} and {@link #matchAll} do, under
   * {@code control}: an MCES without excluded target atoms is searched from the smaller molecule,
   * by atoms when connected and by bonds otherwise, and its mappings are turned round, as in the
   * native core.
   */
  static SearchResult search(MolGraph query, MolGraph target, Options options, int limit, boolean all,
      SearchResult.Control control) {
    if (options.objective != Objective.MCES || !options.excludedTargetAtoms.isEmpty()
        || ((options.connectivity == Connectivity.Connected || options.maxFragments == 1)
                ? query.n <= target.n : query.bondCount() <= target.bondCount()))
      return qualifyRingProof(searchOriented(query, target, options, limit, all, control),
          query, target, options.chemistry);
    SearchResult reversed = searchOriented(target, query, options, limit, all, control);
    List<SearchResult.Match> matches = new ArrayList<>();
    for (SearchResult.Match match : reversed.matches()) {
      Map<Integer, Integer> atoms = new TreeMap<>();
      match.atomMapping().forEach((targetAtom, queryAtom) -> atoms.put(queryAtom, targetAtom));
      Map<SearchResult.Bond, SearchResult.Bond> bonds = new TreeMap<>((first, second) -> {
        int order = Integer.compare(first.first(), second.first());
        return order != 0 ? order : Integer.compare(first.second(), second.second());
      });
      match.bondMapping().forEach((targetBond, queryBond) -> bonds.put(queryBond, targetBond));
      matches.add(new SearchResult.Match(
          atoms, bonds, MatchGraph.fragments(atoms, bonds.keySet()), match.score()));
    }
    matches.sort(CommonBondSearch.CANONICAL);
    return qualifyRingProof(new SearchResult(matches, reversed.status(), reversed.optimal(), reversed.complete(),
        reversed.score(), reversed.nodesVisited(), control.elapsedMillis()), query, target, options.chemistry);
  }

  /** Whether the chemistry relies on an exhaustive set of relevant cycles. */
  private static boolean needsRingProof(ChemOptions chemistry) {
    return chemistry.completeRingsOnly || chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT;
  }

  /** Retain the search results without claiming proof from truncated ring perception. */
  private static SearchResult qualifyRingProof(
      SearchResult result, MolGraph query, MolGraph target, ChemOptions chemistry) {
    if (!needsRingProof(chemistry)
        || (!result.optimal() && !result.complete() && result.status() != SearchResult.Status.Completed)
        || (!query.ringsTruncated() && !target.ringsTruncated()))
      return result;
    SearchResult.Status status = result.status();
    if (status == SearchResult.Status.Completed || status == SearchResult.Status.LimitReached)
      status = SearchResult.Status.LimitReached;
    return new SearchResult(result.matches(), status, false, false, result.score(),
        result.nodesVisited(), result.elapsedMillis());
  }

  /**
   * Search one orientation. Every MCES goes to the common-bond search, which bounds the common
   * bonds by their labels. A single best MCIS in any number of fragments, without further
   * constraints, goes to the common-atom search; any other MCIS, or a graph that search cannot
   * take, goes to the general kernel.
   */
  private static SearchResult searchOriented(
      MolGraph query, MolGraph target, Options options, int limit, boolean all, SearchResult.Control control) {
    if (options.objective == Objective.MCES)
      return CommonBondSearch.search(query, target, options, limit, all, null, control);
    if (!all && commonAtomSearchApplies(options)) {
      SearchResult result = CommonAtomSearch.search(query, target, options.chemistry,
          options.excludedTargetAtoms, options.minFragmentAtoms, control);
      if (result != null)
        return result;
    }
    return new VertexSearch(query, target, options, limit, all, control).search();
  }

  /**
   * Whether the common-atom search applies: a single best MCIS in any number of fragments,
   * without atom weights, node limit, stereochemistry, complete rings or tautomers, and with no
   * minimum fragment size above one bond. A connected MCIS stays with the general kernel, which
   * is faster there.
   */
  private static boolean commonAtomSearchApplies(Options options) {
    ChemOptions chemistry = options.chemistry;
    return options.objective == Objective.MCIS && options.connectivity == Connectivity.Fragments
        && options.maxFragments == Integer.MAX_VALUE && !hasAtomWeights(options)
        && options.nodeLimit == 0 && !chemistry.useChirality && !chemistry.useBondStereo
        && !chemistry.completeRingsOnly && !chemistry.tautomerAware
        && options.minFragmentAtoms <= 2;
  }

  /** Whether the options weight the atoms. An empty array of weights counts as none, as in C++. */
  private static boolean hasAtomWeights(Options options) {
    return options.atomWeights != null && options.atomWeights.length > 0;
  }

  /** Search independent targets in turn with this immutable prepared query. */
  public List<SearchResult> matchMany(List<MolGraph> targets) {
    return matchMany(targets, 1);
  }

  /**
   * Search independent targets on up to {@code numThreads} threads, as the native core and
   * Python do: zero uses one thread per available processor and one searches them in turn, as
   * {@link #matchMany(List)} does. Each target receives the configured per-call budget, and the
   * results keep the order of the targets.
   *
   * @throws IllegalArgumentException if {@code numThreads} is negative
   */
  public List<SearchResult> matchMany(List<MolGraph> targets, int numThreads) {
    return SearchResult.searchEach(targets, numThreads, this::match);
  }

  /**
   * Assign queries successively, never reusing a target atom selected for an earlier query. The
   * list holds one result per query, in query order, as in the native core; once the cancellation
   * token is cancelled, each remaining query returns at once with status {@code Cancelled}.
   */
  public static List<SearchResult> matchConstrained(
      List<MolGraph> queries, MolGraph target, Options options) {
    Objects.requireNonNull(queries, "queries");
    Objects.requireNonNull(target, "target");
    Set<Integer> excluded = new HashSet<>();
    List<SearchResult> results = new ArrayList<>();
    for (MolGraph query : queries) {
      Options current = copy(options, query.n);
      current.excludedTargetAtoms.addAll(excluded);
      SearchResult result = find(query).withOptions(current).match(target);
      results.add(result);
      if (result.found())
        excluded.addAll(result.matches().get(0).atomMapping().values());
    }
    return Collections.unmodifiableList(results);
  }

  /**
   * An upper bound on the atom-based Tanimoto similarity from the atom labels alone, as in the
   * native core: the number of atoms the two molecules can pair by element, by aromaticity as
   * well under strict aromaticity, and by ring membership as well when ring atoms match ring
   * atoms only, divided by the atoms of the union. Graph constraints and the other chemistry
   * options can only reduce the similarity, so they do not enter the bound. {@code null}
   * chemistry means the defaults.
   */
  public static double similarityUpperBound(
      MolGraph query, MolGraph target, ChemOptions chemistry) {
    Objects.requireNonNull(query, "query");
    Objects.requireNonNull(target, "target");
    chemistry = chemistry == null ? new ChemOptions() : chemistry;
    if (query.n == 0 || target.n == 0)
      return 0.0;
    Map<Integer, Integer> unpaired = new HashMap<>();
    for (int atom = 0; atom < query.n; atom++)
      unpaired.merge(boundLabel(query, atom, chemistry), 1, Integer::sum);
    int overlap = 0;
    for (int atom = 0; atom < target.n; atom++) {
      int label = boundLabel(target, atom, chemistry);
      Integer remaining = unpaired.get(label);
      if (remaining == null || remaining == 0)
        continue;
      unpaired.put(label, remaining - 1);
      overlap++;
    }
    overlap = Math.min(overlap, Math.min(query.n, target.n));
    return overlap == 0 ? 0.0 : (double) overlap / (query.n + target.n - overlap);
  }

  public static double similarityUpperBound(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry) {
    return similarityUpperBound(new MolGraph(query), new MolGraph(target), chemistry);
  }

  /** The label two atoms must share to pair in the bound, as the native core computes it. */
  private static int boundLabel(MolGraph graph, int atom, ChemOptions chemistry) {
    if (!chemistry.matchAtomType)
      return 0;
    int label = graph.atomicNum[atom] << 2;
    if (chemistry.aromaticityMode == ChemOptions.AromaticityMode.STRICT && graph.aromatic[atom])
      label |= 2;
    if (chemistry.ringMatchesRingOnly && graph.ring[atom])
      label |= 1;
    return label;
  }

  /** Describe an explicitly selected fragment using original query indices. */
  public static SearchResult.Match fragment(
      MolGraph query, Map<Integer, Integer> atoms, Collection<SearchResult.Bond> selectedBonds) {
    Objects.requireNonNull(query, "query");
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    Set<Integer> usedTargets = new HashSet<>();
    for (Map.Entry<Integer, Integer> pair : atoms.entrySet()) {
      if (pair.getKey() < 0 || pair.getKey() >= query.n || pair.getValue() < 0
          || !usedTargets.add(pair.getValue()))
        throw new IllegalArgumentException("Fragment atom mapping must be injective and in range");
    }
    for (SearchResult.Bond bond : selectedBonds) {
      if (!atoms.containsKey(bond.first()) || !atoms.containsKey(bond.second())
          || !query.hasBond(bond.first(), bond.second()))
        throw new IllegalArgumentException("Invalid selected fragment bond");
      bonds.put(bond, new SearchResult.Bond(atoms.get(bond.first()), atoms.get(bond.second())));
    }
    return new SearchResult.Match(atoms, bonds, MatchGraph.fragments(atoms, bonds.keySet()), atoms.size());
  }

  /**
   * Validate selected atom and bond correspondences against the declared objective, with the
   * rules and messages of the native core: the atoms first, then the selected bonds, then the
   * objective, connectivity, fragment, stereo and ring constraints.
   */
  public static List<String> validateMapping(
      MolGraph query, MolGraph target, SearchResult.Match match, Options options) {
    Options checked = copy(options, query.n);
    // Strict ring fusion compares ring counts, which are worked out only when asked for.
    if (checked.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    List<String> errors = new ArrayList<>();
    Map<Integer, Integer> atoms = match.atomMapping();
    Set<Integer> used = new HashSet<>();
    for (Map.Entry<Integer, Integer> pair : atoms.entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (queryAtom < 0 || queryAtom >= query.n || targetAtom < 0 || targetAtom >= target.n) {
        errors.add("Mapped atom index out of range");
        continue;
      }
      if (!used.add(targetAtom))
        errors.add("Target atom occurs more than once");
      if (checked.excludedTargetAtoms.contains(targetAtom))
        errors.add("Mapped target atom is excluded");
      if (!MolGraph.ChemOps.atomsCompatFast(query, queryAtom, target, targetAtom, checked.chemistry))
        errors.add("Mapped atoms are incompatible");
    }
    if (!errors.isEmpty())
      return errors;
    Set<Integer> incident = new HashSet<>();
    for (Map.Entry<SearchResult.Bond, SearchResult.Bond> pair : match.bondMapping().entrySet()) {
      SearchResult.Bond bond = pair.getKey(), mapped = pair.getValue();
      Integer first = atoms.get(bond.first()), second = atoms.get(bond.second());
      if (first == null || second == null || !new SearchResult.Bond(first, second).equals(mapped)) {
        errors.add("Selected bond is outside the atom mapping");
        continue;
      }
      incident.add(bond.first());
      incident.add(bond.second());
      // The bond rule rejects a pair that is not bonded in either molecule.
      if (!MolGraph.ChemOps.bondsCompatible(
              query, bond.first(), bond.second(), target, first, second, checked.chemistry))
        errors.add("Selected bonds are incompatible");
    }
    if (!errors.isEmpty())
      return errors;
    if (checked.objective == Objective.MCIS) {
      for (int first : atoms.keySet())
        for (int second : atoms.keySet())
          if (first < second) {
            boolean queryBond = query.hasBond(first, second);
            if (queryBond != target.hasBond(atoms.get(first), atoms.get(second)))
              errors.add("Mapped topology is not induced");
            if (queryBond != match.bondMapping().containsKey(new SearchResult.Bond(first, second)))
              errors.add("MCIS must select every induced bond");
          }
    } else if (!incident.equals(atoms.keySet()))
      errors.add("MCES cannot retain isolated atoms");
    List<Set<Integer>> components = MatchGraph.fragments(atoms, match.bondMapping().keySet());
    if (checked.connectivity == Connectivity.Connected && components.size() > 1)
      errors.add("Selected graph is disconnected");
    if (components.size() > checked.maxFragments)
      errors.add("Too many fragments");
    for (Set<Integer> component : components)
      if (component.size() < checked.minFragmentAtoms)
        errors.add("Fragment is too small");
    int[] indexed = new int[query.n];
    Arrays.fill(indexed, -1);
    atoms.forEach((queryAtom, targetAtom) -> indexed[queryAtom] = targetAtom);
    if (Substructure.stereoMismatchAtom(
            query, target, indexed, checked.chemistry, match.bondMapping().keySet())
        >= 0)
      errors.add("Selected stereochemistry is incompatible");
    if (checked.chemistry.completeRingsOnly
        && (!VertexSearch.completeRings(query, atoms.keySet(), match.bondMapping().keySet())
            || !VertexSearch.completeRings(target, new HashSet<>(atoms.values()),
                new HashSet<>(match.bondMapping().values()))))
      errors.add("Selected ring is incomplete");
    return errors;
  }

  /**
   * Prove whether any valid extension retains every supplied correspondence and selected bond.
   * An invalid match is not maximal, so it gives {@code false} rather than an exception, as in
   * the native core.
   *
   * @throws IllegalStateException if the search budget or truncated ring perception prevents proof
   */
  public static boolean isMappingMaximal(
      MolGraph query, MolGraph target, SearchResult.Match match, Options options) {
    if (!validateMapping(query, target, match, options).isEmpty())
      return false;
    Options extension = copy(options, query.n);
    extension.atomWeights = null;
    SearchResult result;
    if (extension.objective == Objective.MCES) {
      result = CommonBondSearch.search(query, target, extension, 1, false, match,
          new SearchResult.Control(
              SearchResult.Control.resolveTimeout(extension.timeoutMillis, query.n, target.n, true),
              extension.cancellation));
    } else {
      VertexSearch kernel = new VertexSearch(query, target, extension, 1, false);
      kernel.seed(match);
      result = kernel.search();
    }
    result = qualifyRingProof(result, query, target, extension.chemistry);
    if (!result.optimal())
      throw new IllegalStateException("Mapping maximality was not proved: " + result.status());
    if (!result.found())
      return match.atomCount() == 0;
    SearchResult.Match largest = result.matches().get(0);
    return largest.atomCount() == match.atomCount() && largest.bondCount() == match.bondCount();
  }

  /**
   * Prove maximality of an atom mapping, selecting its bonds as {@link #mappingMatch} does; an
   * invalid mapping gives {@code false}.
   *
   * @throws IllegalStateException if the search budget or truncated ring perception prevents proof
   */
  public static boolean isMappingMaximal(
      MolGraph query, MolGraph target, Map<Integer, Integer> mapping, Options options) {
    if (!validateMapping(query, target, mapping, options).isEmpty())
      return false;
    return isMappingMaximal(query, target, mappingMatch(query, target, mapping, options), options);
  }

  /**
   * Prove maximality of an atom mapping under chemistry options alone, as a connected induced
   * common subgraph; an invalid mapping gives {@code false}.
   *
   * @deprecated Use {@link #isMappingMaximal(MolGraph, MolGraph, Map, Options)}, which takes the
   *     objective and the connectivity, fragment and resource options as well.
   */
  @Deprecated
  public static boolean isMappingMaximal(
      MolGraph query, MolGraph target, Map<Integer, Integer> mapping, ChemOptions chemistry) {
    Options options = new Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    return isMappingMaximal(query, target, mapping, options);
  }

  /** A progressive consensus and its independently supported molecule mappings. */
  public static final class ConsensusResult {
    private static final SearchResult.Match NO_FRAGMENT =
        new SearchResult.Match(Collections.emptyMap(), Collections.emptyMap(), Collections.emptyList(), 0);
    private final MolGraph source;
    private final int sourceIndex;
    private final SearchResult.Match fragment;
    private final Map<Integer, SearchResult.Match> mappings;
    private final SearchResult.Status status;
    private final long nodesVisited, elapsedMillis;
    private ConsensusResult(MolGraph source, int sourceIndex, SearchResult.Match fragment,
        Map<Integer, SearchResult.Match> mappings, SearchResult.Status status, long nodes,
        long elapsed) {
      this.source = source;
      this.sourceIndex = sourceIndex;
      this.fragment = fragment == null ? NO_FRAGMENT : fragment;
      this.mappings = Collections.unmodifiableMap(new TreeMap<>(mappings));
      this.status = status;
      this.nodesVisited = nodes;
      this.elapsedMillis = elapsed;
    }
    /** The common graph in the atom indices of input {@link #sourceIndex()}; empty when nothing was found. */
    public SearchResult.Match fragment() {
      return fragment;
    }
    /** The input the fragment was taken from, or -1 when nothing was found. */
    public int sourceIndex() {
      return sourceIndex;
    }
    /** The common graph, compact and numbered from 0; empty when nothing was found. */
    public MolGraph graph() {
      if (!found()) return new MolGraph.Builder().atomCount(0).atomicNumbers(new int[0]).neighbors(new int[0][]).build();
      return fragment.toGraph(source);
    }
    /** The common graph as canonical SMILES, or an empty string when nothing was found. */
    public String toSmiles() {
      return found() ? fragment.toSmiles(source) : "";
    }
    public Map<Integer, SearchResult.Match> mappings() {
      return mappings;
    }
    public int supportCount() {
      return mappings.size();
    }
    public boolean found() {
      return fragment.atomCount() > 0;
    }
    public boolean optimal() {
      return false;
    }
    public boolean complete() {
      return status == SearchResult.Status.Completed;
    }
    public SearchResult.Status status() {
      return status;
    }
    public long nodesVisited() {
      return nodesVisited;
    }
    public long elapsedMillis() {
      return elapsedMillis;
    }
  }

  /**
   * Progressively reduce candidate seeds until at least ceil(threshold*n) inputs support them.
   *
   * <p>Each input seeds a candidate in turn, as in the native core: the seed itself, then the
   * common subgraph of the candidate with each other input in input order, found with the search
   * {@link #match} runs. Each candidate is embedded in every input, including its seed, to count
   * its support, and the best supported candidate by the objective is kept; a later candidate of
   * equal score does not replace it. One budget covers every search of the call. An empty list
   * gives a result that is not found, with {@code sourceIndex() == -1}.
   *
   * @throws IllegalArgumentException if the threshold lies outside (0, 1], or the options carry
   *     atom weights or excluded target atoms, which have no single target here
   */
  public static ConsensusResult consensus(
      List<MolGraph> molecules, Options options, double threshold) {
    Objects.requireNonNull(molecules, "molecules");
    Objects.requireNonNull(options, "options");
    if (!Double.isFinite(threshold) || threshold <= 0 || threshold > 1)
      throw new IllegalArgumentException("Consensus threshold must be in (0, 1]");
    if (hasAtomWeights(options))
      throw new IllegalArgumentException("Consensus does not accept atom weights");
    if (!Objects.requireNonNull(options.excludedTargetAtoms, "excludedTargetAtoms").isEmpty())
      throw new IllegalArgumentException("Consensus has no single target for exclusions");
    int largest = 0;
    for (MolGraph graph : molecules)
      largest = Math.max(largest, Objects.requireNonNull(graph, "molecule").n);
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, largest, largest, true),
        options.cancellation);
    ConsensusSearch search = new ConsensusSearch(molecules, options,
        (int) Math.ceil(threshold * molecules.size()), control);
    for (int seed = 0; seed < molecules.size() && !control.stopped(); seed++) search.reduce(seed);
    return search.result();
  }

  /** The state of one consensus call: its inputs, budget and best supported candidate so far. */
  private static final class ConsensusSearch {
    private final List<MolGraph> molecules;
    private final Options options;
    private final int required;
    private final SearchResult.Control control;
    private long nodes;
    private MolGraph bestSource;
    private int bestSourceIndex = -1;
    private SearchResult.Match best;
    private Map<Integer, SearchResult.Match> bestMappings = Collections.emptyMap();
    private long[] bestScore;

    ConsensusSearch(List<MolGraph> molecules, Options options, int required, SearchResult.Control control) {
      this.molecules = molecules;
      this.options = options;
      this.required = required;
      this.control = control;
    }

    ConsensusResult result() {
      SearchResult.Status status = control.status();
      if (status == SearchResult.Status.Completed && needsRingProof(options.chemistry))
        for (MolGraph molecule : molecules)
          if (molecule.ringsTruncated()) {
            status = SearchResult.Status.LimitReached;
            break;
          }
      return new ConsensusResult(bestSource, bestSourceIndex, best, bestMappings, status, nodes,
          control.elapsedMillis());
    }

    /** Whether the node limit of the call is spent; the control then reports the limit. */
    private boolean nodesSpent() {
      if (options.nodeLimit == 0 || nodes < options.nodeLimit) return false;
      control.limitReached();
      return true;
    }

    /** Reduce the seed by each other input in turn, considering the seed and every reduction. */
    void reduce(int seedIndex) {
      MolGraph seed = molecules.get(seedIndex);
      MolGraph current = MolGraph.copyForSearch(seed);
      int[] original = new int[seed.n];
      for (int atom = 0; atom < seed.n; atom++) original[atom] = atom;
      consider(current, sourceFragment(seed, current, original), seedIndex);
      for (int molecule = 0; molecule < molecules.size() && !control.stopped(); molecule++) {
        if (molecule == seedIndex) continue;
        if (nodesSpent()) break;
        Options pair = copy(options, current.n);
        if (options.nodeLimit != 0) pair.nodeLimit = options.nodeLimit - nodes;
        // The same search as match(), so the native core and Java reduce alike.
        SearchResult result = search(current, molecules.get(molecule), pair, 1, false, control);
        nodes += result.nodesVisited();
        if (!result.found()) continue;
        SearchResult.Match selected = result.matches().get(0);
        int[] retained = new int[selected.atomCount()];
        int index = 0;
        for (int atom : selected.atomMapping().keySet()) retained[index++] = original[atom];
        current = selected.toGraph(current);
        original = retained;
        consider(current, sourceFragment(seed, current, original), seedIndex);
      }
    }

    /** The candidate as a fragment of its seed: the kept atoms, and every bond between them. */
    private static SearchResult.Match sourceFragment(MolGraph seed, MolGraph current, int[] original) {
      Map<Integer, Integer> atoms = new TreeMap<>();
      for (int atom : original) atoms.put(atom, atom);
      List<SearchResult.Bond> bonds = new ArrayList<>();
      for (int atom = 0; atom < current.n; atom++)
        for (int neighbour : current.neighbors[atom])
          if (atom < neighbour) bonds.add(new SearchResult.Bond(original[atom], original[neighbour]));
      bonds.sort(Comparator.comparingInt(SearchResult.Bond::first).thenComparingInt(SearchResult.Bond::second));
      return fragment(seed, atoms, bonds);
    }

    /**
     * Count the inputs the candidate embeds in, and keep it when enough do and it beats the best
     * so far. A candidate the options reject, such as a disconnected one under the connected
     * objective, is not considered.
     */
    private void consider(MolGraph graph, SearchResult.Match fragment, int source) {
      if (fragment.atomCount() == 0
          || !validateMapping(molecules.get(source), molecules.get(source), fragment, options).isEmpty())
        return;
      Map<Integer, SearchResult.Match> support = new TreeMap<>();
      for (int molecule = 0; molecule < molecules.size() && !control.stopped(); molecule++) {
        if (nodesSpent()) break;
        SearchResult embedded = embed(graph, molecules.get(molecule));
        if (embedded.found()) support.put(molecule, embedded.matches().get(0));
      }
      if (support.size() < required) return;
      long[] score = score(fragment);
      if (best != null && compare(bestScore, score) >= 0) return;
      best = fragment;
      bestScore = score;
      bestSource = molecules.get(source);
      bestSourceIndex = source;
      bestMappings = support;
    }

    /** Embed the candidate in one input with the chemistry of the options, induced for an MCIS. */
    private SearchResult embed(MolGraph graph, MolGraph molecule) {
      Substructure.Options embedding = new Substructure.Options();
      embedding.chemistry = options.chemistry;
      embedding.induced = options.objective == Objective.MCIS;
      embedding.timeoutMillis = options.timeoutMillis;
      embedding.nodeLimit = options.nodeLimit == 0 ? 0 : options.nodeLimit - nodes;
      embedding.cancellation = options.cancellation;
      SearchResult result = Substructure.find(graph).withOptions(embedding).search(molecule, 1, false, control);
      nodes += result.nodesVisited();
      return result;
    }

    /** The score of a candidate as the native core ranks it: bonds then atoms for an MCES, atoms then bonds for an MCIS. */
    private long[] score(SearchResult.Match fragment) {
      return options.objective == Objective.MCES
          ? new long[] {fragment.bondCount(), fragment.atomCount(), 0}
          : new long[] {fragment.atomCount(), fragment.atomCount(), fragment.bondCount()};
    }

    private static int compare(long[] first, long[] second) {
      for (int term = 0; term < first.length; term++) {
        int order = Long.compare(first[term], second[term]);
        if (order != 0) return order;
      }
      return 0;
    }
  }

  private static Options copy(Options source, int atoms) {
    Objects.requireNonNull(source, "options");
    Options result = new Options();
    result.chemistry = ChemOptions.copyOf(Objects.requireNonNull(source.chemistry, "chemistry"));
    result.objective = Objects.requireNonNull(source.objective, "objective");
    result.connectivity = Objects.requireNonNull(source.connectivity, "connectivity");
    result.atomWeights = hasAtomWeights(source) ? source.atomWeights.clone() : null;
    if (result.atomWeights != null) {
      if (result.objective == Objective.MCES)
        throw new IllegalArgumentException("Atom weights are supported only for MCIS");
      if (result.atomWeights.length < atoms)
        throw new IllegalArgumentException("Atom weights must cover every query atom");
      for (double weight : result.atomWeights)
        if (!Double.isFinite(weight))
          throw new IllegalArgumentException("Atom weights must be finite");
    }
    SearchResult.Control.checkLimits(source.timeoutMillis, source.nodeLimit);
    if (source.minFragmentAtoms < 1 || source.maxFragments < 1)
      throw new IllegalArgumentException("Fragment limits must be positive");
    result.minFragmentAtoms = source.minFragmentAtoms;
    result.maxFragments = source.maxFragments;
    result.timeoutMillis = source.timeoutMillis;
    result.nodeLimit = source.nodeLimit;
    result.cancellation = source.cancellation;
    result.excludedTargetAtoms =
        new HashSet<>(Objects.requireNonNull(source.excludedTargetAtoms, "excludedTargetAtoms"));
    for (int targetAtom : result.excludedTargetAtoms)
      if (targetAtom < 0)
        throw new IllegalArgumentException("Excluded atom indices must be non-negative");
    return result;
  }

  /**
   * The match of an atom mapping, as the native core builds it: every query bond between mapped
   * atoms whose images the target bonds compatibly is selected, so that the mapping can be
   * checked as an induced or edge subgraph under {@code options}.
   *
   * @throws IllegalArgumentException for a mapped atom outside either graph
   */
  public static SearchResult.Match mappingMatch(
      MolGraph query, MolGraph target, Map<Integer, Integer> atoms, Options options) {
    Objects.requireNonNull(query, "query");
    Objects.requireNonNull(target, "target");
    Objects.requireNonNull(atoms, "atoms");
    ChemOptions chemistry = copy(options, query.n).chemistry;
    for (Map.Entry<Integer, Integer> pair : atoms.entrySet())
      if (pair.getKey() < 0 || pair.getKey() >= query.n || pair.getValue() < 0 || pair.getValue() >= target.n)
        throw new IllegalArgumentException("Mapped atom index out of range");
    if (chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    Map<Integer, Integer> sorted = new TreeMap<>(atoms);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    for (Map.Entry<Integer, Integer> pair : sorted.entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      for (int neighbor : query.neighbors[queryAtom]) {
        Integer image = sorted.get(neighbor);
        if (queryAtom < neighbor && image != null && target.hasBond(targetAtom, image)
            && MolGraph.ChemOps.bondsCompatible(query, queryAtom, neighbor, target, targetAtom, image, chemistry))
          bonds.put(new SearchResult.Bond(queryAtom, neighbor), new SearchResult.Bond(targetAtom, image));
      }
    }
    return new SearchResult.Match(sorted, bonds, MatchGraph.fragments(sorted, bonds.keySet()), sorted.size());
  }

  /**
   * Validate an atom mapping with the bonds {@link #mappingMatch} selects, as
   * {@link #validateMapping(MolGraph, MolGraph, SearchResult.Match, Options)} does; a mapped atom
   * outside either graph is reported rather than thrown.
   */
  public static List<String> validateMapping(
      MolGraph query, MolGraph target, Map<Integer, Integer> mapping, Options options) {
    Objects.requireNonNull(query, "query");
    Objects.requireNonNull(target, "target");
    Objects.requireNonNull(mapping, "mapping");
    for (Map.Entry<Integer, Integer> pair : mapping.entrySet())
      if (pair.getKey() < 0 || pair.getKey() >= query.n || pair.getValue() < 0 || pair.getValue() >= target.n)
        return Collections.singletonList("Mapped atom index out of range");
    return validateMapping(query, target, mappingMatch(query, target, mapping, options), options);
  }

  /**
   * Validate an atom mapping under chemistry options alone, as an induced common subgraph in any
   * number of fragments, with the messages of the native core. A {@code null} mapping is valid
   * and {@code null} chemistry means the defaults.
   *
   * @deprecated Use {@link #validateMapping(MolGraph, MolGraph, Map, Options)}, which applies the
   *     objective, connectivity, fragment, stereo and ring rules of the native core.
   */
  @Deprecated
  public static List<String> validateMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping, ChemOptions C) {
    if (mapping == null || mapping.isEmpty())
      return new ArrayList<>();
    Options options = new Options();
    options.chemistry = C == null ? new ChemOptions() : C;
    options.connectivity = Connectivity.Fragments;
    return validateMapping(g1, g2, mapping, options);
  }

  /**
   * The lexicographically smallest mapping equivalent to {@code mapping} under the automorphisms
   * of both graphs, as in the native core.
   *
   * @throws IllegalArgumentException for a mapped atom outside either graph
   * @throws IllegalStateException when either graph's automorphism generators are incomplete,
   *     or the time or orbit state limit of the canonicalisation is exceeded
   */
  public static Map<Integer, Integer> canonicalizeMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping) {
    if (mapping.isEmpty())
      return mapping;
    for (Map.Entry<Integer, Integer> pair : mapping.entrySet())
      if (pair.getKey() < 0 || pair.getKey() >= g1.n || pair.getValue() < 0 || pair.getValue() >= g2.n)
        throw new IllegalArgumentException("Mapped atom index out of range");
    int[][] queryGenerators = g1.getAutomorphismGenerators();
    int[][] targetGenerators = g2.getAutomorphismGenerators();
    if (g1.automorphismGeneratorsTruncated() || g2.automorphismGeneratorsTruncated())
      throw new IllegalStateException(
          "Mapping canonicalisation needs complete automorphism generators");
    Map<Integer, Integer> initial = new TreeMap<>(mapping);
    if (queryGenerators.length == 0 && targetGenerators.length == 0)
      return initial;

    // A local lexicographic descent can get stuck: an orbit path may first
    // require an increasing generator step before reaching its true minimum.
    SearchResult.Control budget = new SearchResult.Control(1_000, null);
    int stateLimit = Math.min(65_536, Math.max(1, 1_000_000 / mapping.size()));
    Set<Map<Integer, Integer>> seen = new HashSet<>();
    List<Map<Integer, Integer>> orbit = new ArrayList<>();
    seen.add(initial);
    orbit.add(initial);
    Map<Integer, Integer> best = initial;
    for (int index = 0; index < orbit.size(); index++) {
      if (budget.stopped())
        throw new IllegalStateException("Mapping canonicalisation exceeded its time limit");
      Map<Integer, Integer> current = orbit.get(index);
      for (int side = 0; side < 2; side++) {
        int[][] generators = side == 0 ? queryGenerators : targetGenerators;
        for (int[] generator : generators) {
          Map<Integer, Integer> next = new TreeMap<>();
          for (Map.Entry<Integer, Integer> pair : current.entrySet())
            next.put(side == 0 ? generator[pair.getKey()] : pair.getKey(),
                side == 0 ? pair.getValue() : generator[pair.getValue()]);
          if (!seen.add(next))
            continue;
          if (orbit.size() >= stateLimit)
            throw new IllegalStateException(
                "Mapping canonicalisation exceeded its orbit state limit");
          orbit.add(next);
          if (compareMappings(next, best) < 0)
            best = next;
        }
      }
    }
    return new LinkedHashMap<>(best);
  }

  private static int compareMappings(Map<Integer, Integer> left, Map<Integer, Integer> right) {
    Iterator<Map.Entry<Integer, Integer>> a = left.entrySet().iterator();
    Iterator<Map.Entry<Integer, Integer>> b = right.entrySet().iterator();
    while (a.hasNext() && b.hasNext()) {
      Map.Entry<Integer, Integer> ap = a.next();
      Map.Entry<Integer, Integer> bp = b.next();
      int keyOrder = Integer.compare(ap.getKey(), bp.getKey());
      if (keyOrder != 0)
        return keyOrder;
      int valueOrder = Integer.compare(ap.getValue(), bp.getValue());
      if (valueOrder != 0)
        return valueOrder;
    }
    return Integer.compare(left.size(), right.size());
  }

  /**
   * Tautomer confidence of an atom mapping alone.
   *
   * @deprecated Use {@link #computeTautomerConfidence(MolGraph, MolGraph, SearchResult.Match, ChemOptions)}.
   */
  @Deprecated
  public static double computeTautomerConfidence(
      IAtomContainer m1, IAtomContainer m2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty())
      return 1.0;
    MolGraph g1 = new MolGraph(m1);
    MolGraph g2 = new MolGraph(m2);
    return computeTautomerConfidence(g1, g2, mcs);
  }

  /**
   * Tautomer confidence of an atom mapping alone.
   *
   * @deprecated Use {@link #computeTautomerConfidence(MolGraph, MolGraph, SearchResult.Match, ChemOptions)},
   *     which takes the selected bonds and the solvent and pH context. This form selects every
   *     mapped bond present in both molecules and applies the tautomer profile at pH 7.4 in water,
   *     so it returns the same value as the four-argument form for that selection. Out-of-range
   *     atoms throw {@link IllegalArgumentException}.
   */
  @Deprecated
  public static double computeTautomerConfidence(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty())
      return 1.0;
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    for (Map.Entry<Integer, Integer> pair : mcs.entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (queryAtom < 0 || queryAtom >= g1.n || targetAtom < 0 || targetAtom >= g2.n)
        throw new IllegalArgumentException("Atom mapping index is outside its molecule");
    }
    for (Map.Entry<Integer, Integer> pair : new TreeMap<>(mcs).entrySet())
      for (int neighbor : g1.neighbors[pair.getKey()]) {
        Integer image = mcs.get(neighbor);
        if (neighbor > pair.getKey() && image != null && g2.hasBond(pair.getValue(), image))
          bonds.put(new SearchResult.Bond(pair.getKey(), neighbor), new SearchResult.Bond(pair.getValue(), image));
      }
    return computeTautomerConfidence(g1, g2, new SearchResult.Match(mcs, bonds, Collections.emptyList(), mcs.size()),
        ChemOptions.tautomerProfile());
  }

  /**
   * Heuristic confidence for selected tautomeric chemistry relaxations. The geometric mean
   * uses the baseline atom weights with the supplied solvent and pH context. Identical
   * selected chemistry has confidence 1. This value is not a probability or an MCS score.
   */
  public static double computeTautomerConfidence(
      MolGraph query, MolGraph target, SearchResult.Match match, ChemOptions options) {
    Objects.requireNonNull(query, "query");
    Objects.requireNonNull(target, "target");
    Objects.requireNonNull(match, "match");
    Objects.requireNonNull(options, "options");
    Objects.requireNonNull(options.solvent, "solvent");
    if (!Double.isFinite(options.pH) || options.pH < 0 || options.pH > 14)
      throw new IllegalArgumentException("pH must be finite and between 0 and 14");
    if (!options.tautomerAware)
      return 1.0;
    boolean[] relaxed = new boolean[query.n];
    for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (queryAtom < 0 || queryAtom >= query.n || targetAtom < 0 || targetAtom >= target.n)
        throw new IllegalArgumentException("Atom mapping index is outside its molecule");
      relaxed[queryAtom] = query.atomicNum[queryAtom] != target.atomicNum[targetAtom]
          || query.aromatic[queryAtom] != target.aromatic[targetAtom];
    }
    for (Map.Entry<SearchResult.Bond, SearchResult.Bond> pair : match.bondMapping().entrySet()) {
      SearchResult.Bond queryBond = pair.getKey(), targetBond = pair.getValue();
      if (queryBond.second() >= query.n || targetBond.second() >= target.n)
        throw new IllegalArgumentException("Bond mapping index is outside its molecule");
      Integer mappedFirst = match.atomMapping().get(queryBond.first());
      Integer mappedSecond = match.atomMapping().get(queryBond.second());
      if (mappedFirst == null || mappedSecond == null
          || !targetBond.equals(new SearchResult.Bond(mappedFirst, mappedSecond))
          || !query.hasBond(queryBond.first(), queryBond.second())
          || !target.hasBond(targetBond.first(), targetBond.second()))
        throw new IllegalArgumentException("Selected bonds must agree with the atom mapping and exist in both molecules");
      if (query.bondOrder(queryBond.first(), queryBond.second())
              != target.bondOrder(targetBond.first(), targetBond.second())
          || query.bondAromatic(queryBond.first(), queryBond.second())
              != target.bondAromatic(targetBond.first(), targetBond.second())) {
        relaxed[queryBond.first()] = true;
        relaxed[queryBond.second()] = true;
      }
    }
    double logSum = 0.0;
    int count = 0;
    float[] queryFractions = null, targetFractions = null;
    for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (relaxed[queryAtom] && (query.tautomerClass[queryAtom] >= 0
          || target.tautomerClass[targetAtom] >= 0)) {
        if (queryFractions == null) {
          queryFractions = TautomerClasses.neutralFractions(query, options.pH);
          targetFractions = TautomerClasses.neutralFractions(target, options.pH);
        }
        double weight = Math.min(contextualTautomerWeight(query, queryAtom, options, queryFractions),
            contextualTautomerWeight(target, targetAtom, options, targetFractions));
        logSum += Math.log(Math.max(weight, 1e-10));
        count++;
      }
    }
    return count == 0 ? 1.0 : Math.exp(logSum / count);
  }

  /**
   * The baseline weight of an atom after the solvent and pH context, in the single-precision
   * steps of the native core, so that the confidences of both languages agree within the shared
   * golden tolerance of 1e-12.
   * The solvent shift comes first; the pH then scales the weight of an atom with an ionisable
   * site by the neutral fraction of its class, from {@link TautomerClasses#neutralFractions}.
   */
  private static double contextualTautomerWeight(MolGraph graph, int atom, ChemOptions options,
      float[] neutralFractions) {
    float weight = graph.tautomerWeight[atom];
    if (graph.tautomerClass[atom] < 0 || weight >= 1.0f)
      return weight;
    if (options.solvent != ChemOptions.Solvent.AQUEOUS) {
      float ketoShift, amideShift, generalShift;
      switch (options.solvent) {
        case DMSO: ketoShift = -0.15f; amideShift = -0.08f; generalShift = -0.05f; break;
        case CHLOROFORM: ketoShift = -0.25f; amideShift = -0.12f; generalShift = -0.08f; break;
        case METHANOL: ketoShift = -0.05f; amideShift = -0.02f; generalShift = -0.01f; break;
        case ACETONITRILE: ketoShift = -0.10f; amideShift = -0.05f; generalShift = -0.03f; break;
        default: ketoShift = -0.20f; amideShift = -0.10f; generalShift = -0.06f; break;
      }
      float shift;
      if (weight >= 0.96f) {
        shift = weight >= 0.97f ? amideShift : ketoShift;
        if (Math.abs(weight - MolGraph.TW_KETO_ENOL) < 0.005f)
          shift = ketoShift;
      } else if (weight >= 0.93f) {
        shift = amideShift;
      } else {
        shift = generalShift;
      }
      weight = Math.max(0.05f, Math.min(0.99f, weight + shift));
    }
    float fraction = neutralFractions[atom];
    if (!Float.isNaN(fraction))
      weight = TautomerClasses.neutralFormWeight(weight, fraction);
    return weight;
  }

  /**
   * The Bemis-Murcko scaffold: the ring atoms, the atoms on the paths between rings and every atom
   * joined to one of these by a double or triple bond, as a molecule capped like a common
   * fragment, so that each stripped substituent is replaced by hydrogen. A stereocentre keeps its
   * configuration when all its neighbours stay, and a double bond when each end keeps a ligand.
   */
  public static IAtomContainer murckoScaffold(IAtomContainer mol) {
    MolGraph source = MolGraph.copyOf(mol);
    int n = source.n;
    boolean hasRing = false;
    for (int atom = 0; atom < n && !hasRing; atom++) hasRing = source.ring[atom];
    if (!hasRing)
      return source.mol;
    boolean[] keep = source.ring.clone();
    int[] previous = new int[n], queue = new int[n];
    boolean[] seen = new boolean[n];
    for (int start = 0; start < n; start++) {
      if (!source.ring[start])
        continue;
      // A breadth-first search over the chain atoms from each ring atom. Its path to a ring atom
      // it reaches is unique, since a second path would make the atoms on both ring atoms.
      Arrays.fill(seen, false);
      int head = 0, tail = 0;
      queue[tail++] = start;
      seen[start] = true;
      while (head < tail) {
        int atom = queue[head++];
        for (int neighbour : source.neighbors[atom]) {
          if (seen[neighbour])
            continue;
          seen[neighbour] = true;
          previous[neighbour] = atom;
          if (source.ring[neighbour])
            for (int kept = neighbour; kept != start; kept = previous[kept]) keep[kept] = true;
          else
            queue[tail++] = neighbour;
        }
      }
    }
    // An atom joined to a ring or linker atom by a double or triple bond stays with that bond, so
    // that a ring carrying an exocyclic C=O, C=N or C=S group can still be kekulised.
    boolean[] framework = keep.clone();
    for (int atom = 0; atom < n; atom++) {
      if (!framework[atom])
        continue;
      for (int neighbour : source.neighbors[atom]) {
        int order = source.bondOrder(atom, neighbour);
        if (!framework[neighbour] && (order == 2 || order == 3))
          keep[neighbour] = true;
      }
    }
    Map<Integer, Integer> retained = new TreeMap<>();
    for (int atom = 0; atom < n; atom++)
      if (keep[atom])
        retained.put(atom, atom);
    IAtomContainer scaffold =
        MatchGraph.toMolecule(source, retained, MatchGraph.allBonds(source, retained).keySet());
    MatchGraph.retainDoubleBondStereo(source, retained.keySet(), scaffold);
    return scaffold;
  }

  /** The hydrogens of a nitrogen, oxygen or sulfur atom, which a tautomer shift can move. */
  private static int mobilisableProtons(MolGraph graph, int atom) {
    int element = graph.atomicNum[atom];
    return element == 7 || element == 8 || element == 16 ? graph.hydrogenCount(atom) : 0;
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
    if (mcs == null || mcs.isEmpty())
      return true;
    if (g1.tautomerClass == null || g2.tautomerClass == null)
      return true;

    Map<Integer, Integer> querySideProtons = new HashMap<>();
    Map<Integer, Integer> targetSideProtons = new HashMap<>();
    for (Map.Entry<Integer, Integer> e : mcs.entrySet()) {
      int qi = e.getKey(), ti = e.getValue();
      if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n)
        continue;
      int tc1 = g1.tautomerClass[qi];
      int tc2 = g2.tautomerClass[ti];
      if (tc1 < 0 && tc2 < 0)
        continue;
      int groupKey = (tc1 >= 0) ? tc1 : tc2 + 100000;
      querySideProtons.merge(groupKey, mobilisableProtons(g1, qi), Integer::sum);
      targetSideProtons.merge(groupKey, mobilisableProtons(g2, ti), Integer::sum);
    }
    for (Map.Entry<Integer, Integer> e : querySideProtons.entrySet()) {
      int tProtons = targetSideProtons.getOrDefault(e.getKey(), 0);
      if (!e.getValue().equals(tProtons))
        return false;
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
    if (mcs == null || mcs.isEmpty())
      return true;
    return validateTautomerConsistency(new MolGraph(m1), new MolGraph(m2), mcs);
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.lang.ref.SoftReference;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BooleanSupplier;
import org.openscience.cdk.AtomRef;
import org.openscience.cdk.BondRef;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.graph.Cycles;
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
  /** The compiled searches of the SMARTS query, shared by the queries that differ only in options. */
  private final SmartsSearches smartsSearches;
  /** Pairs of atomic number and the number of SMARTS atoms that require it. */
  private final int[] requiredElements;
  /** The neighbours of each SMARTS atom in bond order, as CDK lists connected atoms. */
  private final int[][] smartsNeighbors;
  private final Options options;
  private final QueryTables queryTables;
  /** The query bonds as matches report them, by first atom and then in neighbour order. */
  private final SearchResult.QueryBonds queryBonds;
  /** The work counted for building one match: one, and one per query atom and listed query bond. */
  private final int matchWork;
  /** The connected fragment of each query atom, numbered by first appearance in atom order. */
  private final int[] fragmentOf;
  /** The fragments as unmodifiable sorted sets, in fragment order. */
  private final List<Set<Integer>> fragments;

  private Substructure(MolGraph moleculeQuery, IAtomContainer smartsQuery, SmartsSearches smartsSearches,
                       Options options) {
    this.moleculeQuery = moleculeQuery;
    this.smartsQuery = smartsQuery;
    this.smartsSearches = smartsSearches;
    this.requiredElements = smartsQuery == null ? null : requiredElements(smartsQuery, smartsSearches.numericRingCounts);
    this.smartsNeighbors = smartsQuery == null ? null : neighbors(smartsQuery);
    this.options = copyOptions(options);
    // Query-only search data is prepared once here rather than in every call.
    this.queryTables = moleculeQuery == null ? null : new QueryTables(moleculeQuery, this.options.chemistry);
    int[][] neighbors = moleculeQuery != null ? moleculeQuery.neighbors : smartsNeighbors;
    int queryAtoms = neighbors.length, listed = 0;
    // Every embedding retains the complete query topology.
    Set<SearchResult.Bond> bonds = new LinkedHashSet<>();
    for (int atom = 0; atom < queryAtoms; atom++)
      for (int neighbor : neighbors[atom]) if (neighbor > atom) {
        bonds.add(new SearchResult.Bond(atom, neighbor));
        listed++;
      }
    this.queryBonds = new SearchResult.QueryBonds(new ArrayList<>(bonds));
    this.matchWork = 1 + queryAtoms + listed;
    List<Set<Integer>> sets = new ArrayList<>();
    int[] component = new int[queryAtoms], pending = new int[queryAtoms];
    Arrays.fill(component, -1);
    for (int atom = 0; atom < queryAtoms; atom++) if (component[atom] < 0) {
      Set<Integer> fragment = new TreeSet<>();
      int first = 0, end = 0;
      component[atom] = sets.size();
      pending[end++] = atom;
      while (first < end) {
        int current = pending[first++];
        fragment.add(current);
        for (int neighbor : neighbors[current]) if (component[neighbor] < 0) {
          component[neighbor] = sets.size();
          pending[end++] = neighbor;
        }
      }
      sets.add(Collections.unmodifiableSet(fragment));
    }
    this.fragmentOf = component;
    this.fragments = Collections.unmodifiableList(sets);
  }

  private static int[] requiredElements(IAtomContainer query, boolean numericRingCounts) {
    Map<Integer, Integer> counts = new TreeMap<>();
    int[] remaining = numericRingCounts ? new int[] {MAX_SMARTS_EXPRESSION_NODES} : null;
    for (IAtom atom : query.atoms()) {
      IAtom source = AtomRef.deref(atom);
      if (!(source instanceof QueryAtom)) continue;
      Expr expression = ((QueryAtom) source).getExpression();
      int element = requiredElement(expression);
      if (element < 0 && numericRingCounts) element = ringProofElement(expression, remaining, 0);
      if (element >= 0) counts.merge(element, 1, Integer::sum);
    }
    int[] pairs = new int[counts.size() * 2];
    int index = 0;
    for (Map.Entry<Integer, Integer> entry : counts.entrySet()) {
      pairs[index++] = entry.getKey();
      pairs[index++] = entry.getValue();
    }
    return pairs;
  }

  /** The atomic number every matching atom has, or -1 when the expression does not fix one. */
  private static int requiredElement(Expr expression) {
    if (expression == null) return -1;
    switch (expression.type()) {
      case ELEMENT:
      case ALIPHATIC_ELEMENT:
      case AROMATIC_ELEMENT:
        return expression.value();
      case AND: {
        int left = requiredElement(expression.left()), right = requiredElement(expression.right());
        if (left >= 0 && right >= 0 && left != right) return -1;
        return Math.max(left, right);
      }
      case OR: {
        int left = requiredElement(expression.left());
        return left >= 0 && left == requiredElement(expression.right()) ? left : -1;
      }
      default:
        return -1;
    }
  }

  /** A necessary anchor element for the numeric-ring absence screen; no recursive match is evaluated. */
  private static int ringProofElement(Expr expression, int[] remaining, int depth) {
    if (expression == null || depth >= 64 || remaining[0] == 0) return -1;
    remaining[0]--;
    switch (expression.type()) {
      case ELEMENT:
      case ALIPHATIC_ELEMENT:
      case AROMATIC_ELEMENT:
        return expression.value() >= 0 && expression.value() < 256 ? expression.value() : -1;
      case RECURSIVE: {
        IAtomContainer nested = expression.subquery();
        if (nested == null || nested.getAtomCount() == 0) return -1;
        IAtom anchor = AtomRef.deref(nested.getAtom(0));
        return anchor instanceof QueryAtom
            ? ringProofElement(((QueryAtom) anchor).getExpression(), remaining, depth + 1) : -1;
      }
      case AND: {
        int left = ringProofElement(expression.left(), remaining, depth + 1);
        int right = ringProofElement(expression.right(), remaining, depth + 1);
        if (left >= 0 && right >= 0 && left != right) return -1;
        return Math.max(left, right);
      }
      case OR: {
        int left = ringProofElement(expression.left(), remaining, depth + 1);
        int right = ringProofElement(expression.right(), remaining, depth + 1);
        return left >= 0 && left == right ? left : -1;
      }
      default:
        return -1;
    }
  }

  /** Neighbour lists from one pass over the bonds; a container's own lookups scan every bond. */
  private static int[][] neighbors(IAtomContainer query) {
    int atoms = query.getAtomCount();
    Map<IAtom, Integer> index = new IdentityHashMap<>();
    for (int atom = 0; atom < atoms; atom++) index.put(AtomRef.deref(query.getAtom(atom)), atom);
    int[][] ends = new int[query.getBondCount()][];
    int[] degree = new int[atoms];
    int bond = 0;
    for (IBond edge : query.bonds()) {
      int first = index.get(AtomRef.deref(edge.getBegin())), second = index.get(AtomRef.deref(edge.getEnd()));
      ends[bond++] = new int[] {first, second};
      degree[first]++;
      if (second != first) degree[second]++;
    }
    int[][] neighbors = new int[atoms][];
    for (int atom = 0; atom < atoms; atom++) neighbors[atom] = new int[degree[atom]];
    Arrays.fill(degree, 0);
    for (int[] end : ends) {
      neighbors[end[0]][degree[end[0]]++] = end[1];
      if (end[1] != end[0]) neighbors[end[1]][degree[end[1]]++] = end[0];
    }
    return neighbors;
  }

  /** An injective embedding needs at least as many target atoms of each required element. */
  private boolean elementsAvailable(MolGraph target) {
    for (int index = 0; index < requiredElements.length; index += 2) {
      int element = requiredElements[index], available = 0;
      for (int atom = 0; atom < target.n && available < requiredElements[index + 1]; atom++)
        if (target.atomicNum[atom] == element) available++;
      if (available < requiredElements[index + 1]) return false;
    }
    return true;
  }

  /** Prepare an independent query snapshot, retaining its atom indices. */
  public static Substructure find(MolGraph query) {
    return new Substructure(MolGraph.copyForSearch(query), null, null, new Options());
  }

  public static Substructure find(IAtomContainer query) {
    return new Substructure(MolGraph.copyOf(query), null, null, new Options());
  }

  /** Longest SMARTS accepted, in characters; the same limit as in the native core. */
  public static final int MAX_SMARTS_LENGTH = 100_000;
  /** Most atoms in a SMARTS, counting those of its recursive SMARTS; the native SMILES and SMARTS atom limit. */
  public static final int MAX_SMARTS_ATOMS = 10_000;
  /**
   * Most atom and bond expression nodes in a SMARTS, counting those of its recursive SMARTS:
   * each primitive test and each logical operator joining two operands counts once.
   */
  public static final int MAX_SMARTS_EXPRESSION_NODES = 100_000;
  /** Most nodes in one atom or bond expression of a SMARTS; a recursive SMARTS counts one node there. */
  public static final int MAX_SMARTS_NODES_PER_EXPRESSION = 10_000;
  /** Deepest nesting of recursive SMARTS; the native parser's default. */
  public static final int MAX_SMARTS_RECURSION_DEPTH = 20;

  /**
   * Compile the complete SMARTS expression, including recursive and stereo predicates.
   *
   * @throws IllegalArgumentException for invalid SMARTS, or SMARTS beyond {@link #MAX_SMARTS_LENGTH},
   *     {@link #MAX_SMARTS_ATOMS}, {@link #MAX_SMARTS_EXPRESSION_NODES},
   *     {@link #MAX_SMARTS_NODES_PER_EXPRESSION} or {@link #MAX_SMARTS_RECURSION_DEPTH}, counted
   *     as the native parser counts them
   */
  public static Substructure smarts(String expression) {
    Objects.requireNonNull(expression, "expression");
    IAtomContainer query = SmartsInput.parse(expression);
    return new Substructure(null, query, new SmartsSearches(query), new Options());
  }

  public Substructure withOptions(Options options) {
    return new Substructure(moleculeQuery, smartsQuery, smartsSearches,
        Objects.requireNonNull(options, "options"));
  }

  public Substructure withTimeout(long timeoutMillis) {
    Options changed = copyOptions(options);
    changed.timeoutMillis = timeoutMillis;
    return withOptions(changed);
  }

  public SearchResult match(IAtomContainer target) { return run(MolGraph.copyOf(target), 1, false); }
  public SearchResult match(MolGraph target) { return run(target, 1, false); }
  public SearchResult matchAll(IAtomContainer target) { return matchAll(target, 1000); }
  public SearchResult matchAll(MolGraph target) { return matchAll(target, 1000); }
  public SearchResult matchAll(IAtomContainer target, int limit) { return run(MolGraph.copyOf(target), limit, true); }
  public SearchResult matchAll(MolGraph target, int limit) { return run(target, limit, true); }

  /** Search independent targets in input order with this immutable prepared query. */
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

  /** A retained core, attached groups and the outcome of the substructure search. */
  public static final class Decomposition {
    private final IAtomContainer core;
    private final Map<String, IAtomContainer> rgroups;
    private final SearchResult search;
    private Decomposition(IAtomContainer core, Map<String, IAtomContainer> rgroups, SearchResult search) {
      this.core = core;
      this.rgroups = Collections.unmodifiableMap(new LinkedHashMap<>(rgroups));
      this.search = search;
    }
    public IAtomContainer core() { return core; }
    public Map<String, IAtomContainer> rgroups() { return rgroups; }
    public SearchResult search() { return search; }
  }

  /**
   * Split each molecule around the first embedding of {@code core}, preserving no-match, timeout
   * and cancellation outcomes in {@code search()}. {@code core()} is the matched target atoms
   * with the images of the query bonds only, as in the native core, so a chain core matched in a
   * ring gives a chain; it is {@code null} without a match. The attached groups are the connected
   * components of the remaining atoms, numbered from the lowest core atom outwards.
   */
  public static List<Decomposition> decomposeRGroups(
      IAtomContainer core, List<IAtomContainer> molecules, Options options) {
    Objects.requireNonNull(core, "core");
    Objects.requireNonNull(molecules, "molecules");
    List<Decomposition> decompositions = new ArrayList<>();
    Substructure query = find(core).withOptions(Objects.requireNonNull(options, "options"));
    for (IAtomContainer molecule : molecules) {
      MolGraph target = MolGraph.copyOf(molecule);
      SearchResult result = query.match(target);
      Map<String, IAtomContainer> groups = new LinkedHashMap<>();
      IAtomContainer retainedCore = null;
      if (result.found()) {
        SearchResult.Match match = result.matches().get(0);
        Map<Integer, Integer> coreAtoms = new TreeMap<>();
        for (int atom : match.atomMapping().values()) coreAtoms.put(atom, atom);
        retainedCore = MatchGraph.toMolecule(target, coreAtoms, new HashSet<>(match.bondMapping().values()));
        Set<Integer> visited = new HashSet<>(coreAtoms.keySet());
        int group = 1;
        for (int coreAtom : coreAtoms.keySet()) {
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
            groups.put("R" + group++, component(target, atoms));
          }
        }
      }
      decompositions.add(new Decomposition(retainedCore, groups, result));
    }
    return Collections.unmodifiableList(decompositions);
  }

  /**
   * As {@link #decomposeRGroups(IAtomContainer, List, Options)} with the given chemistry, or the
   * defaults for {@code null}, and time limit.
   */
  public static List<Decomposition> decomposeRGroups(
      IAtomContainer core, List<IAtomContainer> molecules, ChemOptions chemistry, long timeoutMillis) {
    Options options = new Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    options.timeoutMillis = timeoutMillis;
    return decomposeRGroups(core, molecules, options);
  }

  private static IAtomContainer component(MolGraph graph, Set<Integer> atoms) {
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom : atoms) identity.put(atom, atom);
    return MatchGraph.toMolecule(graph, identity, MatchGraph.allBonds(graph, identity).keySet());
  }

  private static Options copyOptions(Options source) {
    Options copy = new Options();
    copy.chemistry = ChemOptions.copyOf(Objects.requireNonNull(source.chemistry, "chemistry"));
    copy.induced = source.induced;
    copy.timeoutMillis = source.timeoutMillis;
    copy.nodeLimit = source.nodeLimit;
    copy.cancellation = source.cancellation;
    SearchResult.Control.checkLimits(copy.timeoutMillis, copy.nodeLimit);
    return copy;
  }

  private SearchResult run(MolGraph target, int limit, boolean all) {
    Objects.requireNonNull(target, "target");
    if (limit < 1) throw new IllegalArgumentException("Match limit must be positive");
    int queryAtoms = moleculeQuery != null ? moleculeQuery.n : smartsQuery.getAtomCount();
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, queryAtoms, target.n, false),
        options.cancellation);
    return search(target, limit, all, control);
  }

  /**
   * The search of {@link #match} and {@link #matchAll} under a control the caller owns, so that
   * one budget can cover several searches, as it does in {@link MCS#consensus}; the node limit of
   * the options bounds this search alone.
   */
  SearchResult search(MolGraph target, int limit, boolean all, SearchResult.Control control) {
    int queryAtoms = moleculeQuery != null ? moleculeQuery.n : smartsQuery.getAtomCount();
    List<Map<Integer, Integer>> mappings = new ArrayList<>();
    long nodes = 0;
    boolean ringDependent = false;
    int searchLimit = all && limit < Integer.MAX_VALUE ? limit + 1 : limit;
    if (!control.stopped() && queryAtoms <= target.n) {
      if (moleculeQuery != null) {
        // Pairs the element and degree counts rule out end without a search node.
        if (queryTables.admits(target)) return qualifyRingProof(searchMolecule(target, limit, all, control), target);
      } else if (elementsAvailable(target)) {
        if (smartsSearches.numericRingCounts) {
          boolean incomplete = target.ringsTruncated();
          if (control.stopped())
            return new SearchResult(Collections.<SearchResult.Match>emptyList(), control.status(), false, false, 0, 0, control.elapsedMillis());
          if (incomplete)
            return new SearchResult(Collections.<SearchResult.Match>emptyList(), SearchResult.Status.LimitReached,
                false, false, 0, 0, control.elapsedMillis());
        }
        ringDependent = queryAtoms > 0;
        SmartsTarget view = SmartsTarget.lend(target);
        SmartsCall call = smartsSearches.lend();
        call.start(control, options.nodeLimit, view,
            options.chemistry.completeRingsOnly || options.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT
                ? queryAtoms : 0);
        InducedCheck induced = options.induced ? new InducedCheck(smartsNeighbors, target, control) : null;
        CompleteRingCheck completeRings = null;
        try {
          check(control);
          for (int[] embedding : call.pattern.matchAll(view.molecule)) {
            check(control);
            if (induced != null && !induced.holds(embedding)) continue;
            if (options.chemistry.completeRingsOnly) {
              if (completeRings == null)
                completeRings = new CompleteRingCheck(null, smartsNeighbors, target, control);
              if (!completeRings.holds(embedding)) continue;
            }
            mappings.add(new Embedding(embedding));
            if (!all || mappings.size() >= searchLimit) break;
          }
        } catch (CommonBondGraphs.Stop stopped) {
          // The shared control retains the reason the search stopped.
        }
        nodes = call.nodes;
        if (call.emptyDomain) ringDependent = false;
        // Returned only after an ordinary end or a stop, so that neither holds the state of a failure.
        smartsSearches.giveBack(call.finish());
        view.giveBack();
      }
    }
    boolean truncated = all && mappings.size() > limit;
    if (truncated) mappings.remove(mappings.size() - 1);
    List<int[]> embeddings = new ArrayList<>(mappings.size());
    for (Map<Integer, Integer> mapping : mappings) {
      // Each mapping holds every query atom, in ascending order, so fragments keep their order.
      int[] targets = new int[mapping.size()];
      for (Map.Entry<Integer, Integer> atom : mapping.entrySet()) targets[atom.getKey()] = atom.getValue();
      embeddings.add(targets);
    }
    SearchResult result = finish(describe(embeddings, null, control), mappings.size(), truncated, queryAtoms, nodes, control);
    return ringDependent ? qualifyRingProof(result, target) : result;
  }

  /** Complete-ring coverage reused within one search of full query embeddings. */
  static final class CompleteRingCheck {
    private final int[][] neighbors;
    private final MolGraph target;
    private final SearchResult.Control control;
    private final boolean queryComplete;
    private final VertexSearch.RingIndex targetRings;
    private final int[] targetToQuery;
    private int stamp;

    CompleteRingCheck(MolGraph query, int[][] neighbors, MolGraph target,
                      SearchResult.Control control) {
      this.neighbors = neighbors;
      this.target = target;
      this.control = control;
      boolean complete = true;
      if (query != null) {
        Set<Integer> atoms = new HashSet<>();
        Set<SearchResult.Bond> bonds = new HashSet<>();
        for (int atom = 0; atom < neighbors.length; atom++) {
          atoms.add(atom);
          for (int neighbor : neighbors[atom]) if (neighbor > atom)
            bonds.add(new SearchResult.Bond(atom, neighbor));
        }
        complete = VertexSearch.completeRings(query, atoms, bonds);
      }
      queryComplete = complete;
      targetRings = complete && !control.stopped() ? new VertexSearch.RingIndex(target) : null;
      targetToQuery = new int[target.n];
      Arrays.fill(targetToQuery, -1);
    }

    boolean holds(int[] mapping) {
      if (control.stopped() || !queryComplete) return false;
      if (stamp == Integer.MAX_VALUE) {
        targetRings.clearStamps();
        stamp = 0;
      }
      ++stamp;
      for (int atom = 0; atom < mapping.length; atom++) targetToQuery[mapping[atom]] = atom;
      try {
        for (int ring = 0; ring < targetRings.cycles.length; ring++) {
          if (control.stopped()) return false;
          int[] cycle = targetRings.cycles[ring];
          boolean whole = true;
          for (int index = 0; index < cycle.length && whole; index++) {
            int first = targetToQuery[cycle[index]];
            int second = targetToQuery[cycle[(index + 1) % cycle.length]];
            whole = first >= 0 && second >= 0 && bonded(first, second);
          }
          if (whole) targetRings.cover(ring, stamp);
        }
        for (int atom = 0; atom < mapping.length; atom++) {
          if (control.stopped() || !targetRings.atomCovered(target, mapping[atom], stamp))
            return false;
          for (int neighbor : neighbors[atom]) if (neighbor > atom
              && !targetRings.bondCovered(target, mapping[atom], mapping[neighbor], stamp))
            return false;
        }
        return !control.stopped();
      } finally {
        for (int targetAtom : mapping) targetToQuery[targetAtom] = -1;
      }
    }

    /** Membership of the query bond set, whose bonds are listed by their lower atom. */
    private boolean bonded(int first, int second) {
      int low = Math.min(first, second), high = Math.max(first, second);
      for (int neighbor : neighbors[low]) if (neighbor == high) return true;
      return false;
    }
  }

  /** Retain results without claiming proof from truncated ring perception. */
  private SearchResult qualifyRingProof(SearchResult result, MolGraph target) {
    ChemOptions chemistry = options.chemistry;
    if ((!chemistry.completeRingsOnly && chemistry.ringFusionMode != ChemOptions.RingFusionMode.STRICT)
        || (!result.optimal() && !result.complete() && result.status() != SearchResult.Status.Completed)
        || ((moleculeQuery == null || !moleculeQuery.ringProofIncomplete(chemistry)) && !target.ringProofIncomplete(chemistry)))
      return result;
    SearchResult.Status status = result.status();
    if (status == SearchResult.Status.Completed || status == SearchResult.Status.LimitReached)
      status = SearchResult.Status.LimitReached;
    return new SearchResult(result.matches(), status, false, false, result.score(),
        result.nodesVisited(), result.elapsedMillis());
  }

  /** The molecule search of a target that the atom counts admit. */
  private SearchResult searchMolecule(MolGraph target, int limit, boolean all, SearchResult.Control control) {
    SubgraphSearch matcher = new SubgraphSearch(moleculeQuery, queryTables, target, options.chemistry, control);
    matcher.induced = options.induced;
    matcher.nodeLimit = options.nodeLimit;
    List<int[]> embeddings = new ArrayList<>();
    // With several fragments, each embedding lists them in the order the search first reached them.
    List<int[]> fragmentOrders = all && fragments.size() > 1 ? new ArrayList<>() : null;
    if (all) {
      matcher.enumerate(limit < Integer.MAX_VALUE ? limit + 1 : limit, embeddings, fragmentOf, fragmentOrders);
    } else if (matcher.exists()) {
      embeddings.add(matcher.q2t.clone());
    }
    boolean truncated = all && embeddings.size() > limit;
    if (truncated) {
      embeddings.remove(embeddings.size() - 1);
      if (fragmentOrders != null) fragmentOrders.remove(fragmentOrders.size() - 1);
    }
    return finish(describe(embeddings, fragmentOrders, control), embeddings.size(), truncated, moleculeQuery.n,
        matcher.nodesVisited, control);
  }

  /** The result of a search that kept {@code kept} embeddings, and found one more when truncated. */
  private static SearchResult finish(List<SearchResult.Match> matches, int kept, boolean truncated, int queryAtoms,
                                     long nodes, SearchResult.Control control) {
    // Matches left unbuilt past the results deadline make the result a timed-out one.
    SearchResult.Status status = truncated && matches.size() == kept
        ? SearchResult.Status.LimitReached : control.status();
    boolean complete = status == SearchResult.Status.Completed;
    return new SearchResult(matches, status, complete || !matches.isEmpty(), complete,
        matches.isEmpty() ? 0 : queryAtoms, nodes, control.elapsedMillis());
  }

  /**
   * Matches of complete embeddings, built while the control allows: the clock is read after about
   * 4096 atoms and bonds, and past the results deadline the matches built so far are returned.
   * The fragments of each match follow its entry of {@code fragmentOrders}, or fragment order
   * when that is null.
   */
  private List<SearchResult.Match> describe(List<int[]> embeddings, List<int[]> fragmentOrders,
                                            SearchResult.Control control) {
    List<SearchResult.Match> matches = new ArrayList<>(embeddings.size());
    long work = 0;
    for (int index = 0; index < embeddings.size(); index++) {
      if ((work += matchWork) >= 4096) {
        work = 0;
        if (control.resultsStopped()) break;
      }
      List<Set<Integer>> ordered = fragments;
      if (fragmentOrders != null) {
        List<Set<Integer>> listed = new ArrayList<>(fragments.size());
        for (int fragment : fragmentOrders.get(index)) listed.add(fragments.get(fragment));
        ordered = Collections.unmodifiableList(listed);
      }
      int[] targets = embeddings.get(index);
      matches.add(SearchResult.Match.of(targets, queryBonds, ordered, targets.length));
    }
    return matches;
  }

  /**
   * Whether an embedding is induced: the target atoms bonded to each image are, among the
   * images, exactly the images of its pattern neighbours. Each pattern atom reads its own
   * neighbours and those of its image, and the control is polled as the reads add up.
   */
  private static final class InducedCheck {
    private final int[][] queryNeighbors;
    private final MolGraph target;
    private final SearchResult.Control control;
    /** The pattern atom mapped to each target atom, or -1; reset after every test. */
    private final int[] owner;
    /** Holds atom + 1 at the pattern neighbours of the atom being tested; only that atom writes the value. */
    private final int[] marks;
    /** The distinct pattern neighbours of each atom other than itself. */
    private final int[] distinct;

    InducedCheck(int[][] queryNeighbors, MolGraph target, SearchResult.Control control) {
      this.queryNeighbors = queryNeighbors;
      this.target = target;
      this.control = control;
      owner = new int[target.n];
      Arrays.fill(owner, -1);
      marks = new int[queryNeighbors.length];
      distinct = new int[queryNeighbors.length];
      for (int atom = 0; atom < queryNeighbors.length; atom++)
        for (int neighbor : queryNeighbors[atom])
          if (neighbor != atom && marks[neighbor] != atom + 1) {
            marks[neighbor] = atom + 1;
            distinct[atom]++;
          }
    }

    boolean holds(int[] mapping) {
      for (int atom = 0; atom < mapping.length; atom++) owner[mapping[atom]] = atom;
      try {
        long work = 0;
        for (int atom = 0; atom < mapping.length; atom++) {
          for (int neighbor : queryNeighbors[atom]) marks[neighbor] = atom + 1;
          int[] adjacent = target.neighbors[mapping[atom]];
          if ((work += 1 + queryNeighbors[atom].length + adjacent.length) >= 4096) {
            work = 0;
            check(control);
          }
          int mapped = 0;
          for (int neighbor : adjacent) {
            int other = owner[neighbor];
            if (other < 0) continue;
            if (marks[other] != atom + 1) return false;
            mapped++;
          }
          if (mapped != distinct[atom]) return false;
        }
        return true;
      } finally {
        for (int image : mapping) owner[image] = -1;
      }
    }
  }

  private static void check(SearchResult.Control control) {
    if (control.stopped()) throw new CommonBondGraphs.Stop();
  }

  /**
   * The compiled searches of one SMARTS query, built once rather than in every call. A call takes
   * one, or compiles another when every one is in use, and returns it when it ends, so that each
   * is used by one call at a time.
   */
  private static final class SmartsSearches {
    private final IAtomContainer query;
    private final boolean numericRingCounts;
    /** A free search, which calls made in turn take and return without allocating. */
    private final AtomicReference<SmartsCall> free;
    /** Further free searches, left by calls that ran at the same time. */
    private final Queue<SmartsCall> spare = new ConcurrentLinkedQueue<>();

    SmartsSearches(IAtomContainer query) {
      this.query = query;
      numericRingCounts = numericRingCounts(query);
      free = new AtomicReference<>(new SmartsCall(query));
    }

    /** Walk all boolean branches and recursive payloads before evaluating a numeric ring count. */
    private static boolean numericRingCounts(IAtomContainer query) {
      Deque<Expr> pending = new ArrayDeque<>();
      Set<Expr> seen = Collections.newSetFromMap(new IdentityHashMap<Expr, Boolean>());
      for (IAtom atom : query.atoms()) pending.push(((QueryAtom) AtomRef.deref(atom)).getExpression());
      while (!pending.isEmpty()) {
        Expr expression = pending.pop();
        if (!seen.add(expression)) continue;
        if (expression.type() == Expr.Type.RING_COUNT && expression.value() > 0) return true;
        switch (expression.type()) {
          case AND:
          case OR:
            pending.push(expression.right());
            pending.push(expression.left());
            break;
          case NOT:
            pending.push(expression.left());
            break;
          case RECURSIVE:
            for (IAtom atom : expression.subquery().atoms())
              pending.push(((QueryAtom) AtomRef.deref(atom)).getExpression());
            break;
          default:
            break;
        }
      }
      return false;
    }

    SmartsCall lend() {
      SmartsCall call = free.getAndSet(null);
      if (call == null) call = spare.poll();
      return call != null ? call : new SmartsCall(query);
    }

    void giveBack(SmartsCall call) {
      if (!free.compareAndSet(null, call)) spare.offer(call);
    }
  }

  /**
   * A target as SMARTS searches read it: a CDK copy with ring flags and the aromaticity model that
   * SmartsPattern applies, and ring data of its atoms worked out on first use. Matching writes to
   * the copy, as CDK caches its adjacency there, so one search at a time reads it. A graph keeps
   * its copy, held softly, from the second SMARTS search that reads it.
   */
  private static final class SmartsTarget {
    /** What a graph holds once a SMARTS search has read it, while it keeps no copy or lends it. */
    private static final Object READ = new Object();
    final MolGraph graph;
    final IAtomContainer molecule;
    /** The reference the graph keeps to this copy, or null when it keeps none. */
    private SoftReference<SmartsTarget> kept;
    /** The smallest ring size of each atom, 0 outside rings and -1 until worked out. */
    private int[] smallestRing;
    /** The distances and the queue of the ring searches of {@link #smallestRing}, made on first use. */
    private int[] distance, queue;
    /** The number of rings through each atom, once counted. */
    private int[] ringCounts;

    private SmartsTarget(MolGraph graph) {
      this.graph = graph;
      molecule = MatchGraph.toMolecule(graph);
      // As SmartsPattern.prepare(): ring flags, marked here without one call per atom of a
      // path, then the aromaticity model that SmartsPattern applies.
      RingPerception.markRingAtomsAndBonds(molecule);
      Aromaticity.apply(Aromaticity.Model.Daylight, molecule);
    }

    /** The copy the graph keeps, when no other search reads it, or else a new copy. */
    static SmartsTarget lend(MolGraph graph) {
      Object state = graph.exchangeSmartsView(READ);
      SmartsTarget view = state instanceof SoftReference ? (SmartsTarget) ((SoftReference<?>) state).get() : null;
      if (view != null) return view;
      view = new SmartsTarget(graph);
      // Most graphs are read once, so a copy is kept from the second read on.
      if (state != null) view.kept = new SoftReference<>(view);
      return view;
    }

    /** Keep the copy for later searches, unless another search has returned one first. */
    void giveBack() {
      if (kept != null) graph.replaceSmartsView(READ, kept);
    }

    /**
     * The size of the smallest ring through {@code atom}, 0 when it is in none: the shortest path
     * back to it over ring bonds from each of its ring neighbours, as the native core reads SMARTS
     * r. CDK's own test can take a ring nearby for one through the atom.
     */
    int smallestRing(int atom) {
      if (smallestRing == null) {
        smallestRing = new int[graph.n];
        Arrays.fill(smallestRing, -1);
      }
      if (smallestRing[atom] >= 0) return smallestRing[atom];
      int best = 0;
      if (graph.ring[atom]) {
        if (distance == null) {
          distance = new int[graph.n];
          queue = new int[graph.n];
        }
        for (int neighbor : graph.neighbors[atom]) {
          if (!graph.bondInRing(atom, neighbor)) continue;
          Arrays.fill(distance, -1);
          int head = 0, tail = 0;
          distance[neighbor] = 1;
          queue[tail++] = neighbor;
          search:
          while (head < tail) {
            int current = queue[head++];
            for (int next : graph.neighbors[current]) {
              if (!graph.bondInRing(current, next) || current == neighbor && next == atom) continue;
              if (next == atom) {
                if (best == 0 || distance[current] + 1 < best) best = distance[current] + 1;
                break search;
              }
              if (distance[next] < 0) {
                distance[next] = distance[current] + 1;
                queue[tail++] = next;
              }
            }
          }
        }
      }
      return smallestRing[atom] = best;
    }

    /**
     * Numeric SMARTS ring membership in the deterministic minimum basis shared with the
     * native core, counted once per target. Search admission requires the basis to be complete.
     */
    int ringCount(int atom) {
      if (ringCounts == null) {
        int[] counts = new int[molecule.getAtomCount()];
        for (int[] path : graph.numericRingBasis())
          for (int ringAtom : path) counts[ringAtom]++;
        ringCounts = counts;
      }
      return ringCounts[atom];
    }
  }

  /** The atom mapping of an embedding, in query atom order, read from it without copying. */
  private static final class Embedding extends AbstractMap<Integer, Integer> {
    private final int[] images;

    Embedding(int[] images) {
      this.images = images;
    }

    @Override public int size() {
      return images.length;
    }

    @Override public boolean containsKey(Object key) {
      return key instanceof Integer && (Integer) key >= 0 && (Integer) key < images.length;
    }

    @Override public Integer get(Object key) {
      if (!containsKey(key)) return null;
      return images[(Integer) key];
    }

    @Override public Set<Entry<Integer, Integer>> entrySet() {
      return new AbstractSet<Entry<Integer, Integer>>() {
        @Override public int size() {
          return images.length;
        }

        @Override public Iterator<Entry<Integer, Integer>> iterator() {
          return new Iterator<Entry<Integer, Integer>>() {
            private int atom;

            @Override public boolean hasNext() {
              return atom < images.length;
            }

            @Override public Entry<Integer, Integer> next() {
              if (atom >= images.length) throw new NoSuchElementException();
              int query = atom++;
              return new SimpleImmutableEntry<>(query, images[query]);
            }
          };
        }
      };
    }
  }

  /**
   * A compiled SMARTS search: a copy of the query whose atom tests count nodes and poll the
   * control of the current call, and its pattern. The query predicates retain CDK recursion,
   * stereo and component semantics. One call at a time uses it, see {@link SmartsSearches}.
   */
  private static final class SmartsCall {
    /** The per-target-atom memo tables of one call use at most sixteen MiB. */
    private static final int RECURSIVE_MEMO_BYTES = 1 << 24;
    private static final int RECURSIVE_MEMO_TABLES = 1024;
    final Pattern pattern;
    private final IAtomContainer rootQuery;
    /** Predicate rejections already observed by this call, used only for ring-proof qualification. */
    private BitSet[] rejectedTargets;
    private int[] rejectionCounts;
    private boolean[] nonemptyDomains;
    private boolean emptyDomain;
    /** The control, node limit, counts and target of the current call. */
    SearchResult.Control control;
    long nodeLimit;
    long nodes;
    int polls;
    SmartsTarget view;
    IdentityHashMap<Expr, byte[]> recursiveResults;
    int recursiveBytes;

    SmartsCall(IAtomContainer query) {
      rootQuery = query;
      pattern = VentoFoggia.findSubstructure(query(query));
    }

    void start(SearchResult.Control control, long nodeLimit, SmartsTarget view, int trackedAtoms) {
      this.control = control;
      this.nodeLimit = nodeLimit;
      this.view = view;
      nodes = 0;
      polls = 0;
      recursiveBytes = 0;
      rejectedTargets = trackedAtoms > 0 ? new BitSet[trackedAtoms] : null;
      rejectionCounts = trackedAtoms > 0 ? new int[trackedAtoms] : null;
      nonemptyDomains = trackedAtoms > 0 ? new boolean[trackedAtoms] : null;
      emptyDomain = false;
    }

    /** End the current call, releasing its control and target. */
    SmartsCall finish() {
      control = null;
      view = null;
      recursiveResults = null;
      rejectedTargets = null;
      rejectionCounts = null;
      nonemptyDomains = null;
      return this;
    }

    /** Record an existing predicate result without evaluating another atom or charging a node. */
    private void observedDomain(int queryAtom, int targetAtom, boolean matches) {
      if (rejectedTargets == null || queryAtom < 0 || nonemptyDomains[queryAtom] || emptyDomain) return;
      if (matches) {
        nonemptyDomains[queryAtom] = true;
        rejectedTargets[queryAtom] = null;
        return;
      }
      BitSet rejected = rejectedTargets[queryAtom];
      if (rejected == null) rejected = rejectedTargets[queryAtom] = new BitSet(view.graph.n);
      if (!rejected.get(targetAtom)) {
        rejected.set(targetAtom);
        if (++rejectionCounts[queryAtom] == view.graph.n) emptyDomain = true;
      }
    }

    /** Whether the expression, outside its recursive SMARTS, has a test of the type. */
    static boolean reads(Expr expression, Expr.Type type) {
      switch (expression.type()) {
        case AND:
        case OR:
          return reads(expression.left(), type) || reads(expression.right(), type);
        case NOT:
          return reads(expression.left(), type);
        default:
          return expression.type() == type;
      }
    }

    /**
     * The expression at {@code atom}, with the smallest ring size and the ring count read from the
     * target's caches. With {@code plainNegation} a negation holds where its operand fails;
     * otherwise, as in CDK, it also holds where the operand is a stereo test, or a disjunction led
     * by one, which is checked once the embedding is complete.
     */
    boolean evaluate(Expr expression, IAtom atom, boolean plainNegation) {
      switch (expression.type()) {
        case AND:
          return evaluate(expression.left(), atom, plainNegation)
              && evaluate(expression.right(), atom, plainNegation);
        case OR:
          return evaluate(expression.left(), atom, plainNegation)
              || evaluate(expression.right(), atom, plainNegation);
        case NOT:
          return !evaluate(expression.left(), atom, plainNegation)
              || !plainNegation && stereoTest(expression.left());
        case RING_SMALLEST:
          return view.smallestRing(atom.getIndex()) == expression.value();
        case RING_COUNT:
          return atom.isInRing() && view.ringCount(atom.getIndex()) == expression.value();
        case RECURSIVE:
          return recursive(expression, atom);
        default:
          return expression.matches(atom);
      }
    }

    /** A recursive predicate depends on this target and anchor, independently of the outer mapping. */
    private boolean recursive(Expr expression, IAtom atom) {
      check(control);
      byte[] results = recursiveResults == null ? null : recursiveResults.get(expression);
      if (results == null) {
        int atoms = view.molecule.getAtomCount();
        if (atoms > RECURSIVE_MEMO_BYTES - recursiveBytes
            || recursiveResults != null && recursiveResults.size() >= RECURSIVE_MEMO_TABLES)
          return expression.matches(atom);
        if (recursiveResults == null) recursiveResults = new IdentityHashMap<>();
        results = new byte[atoms];
        recursiveResults.put(expression, results);
        recursiveBytes += atoms;
      }
      int anchor = atom.getIndex();
      if (results[anchor] != 0) return results[anchor] == 2;
      boolean found = expression.matches(atom);
      // A stopped nested search throws before a partial answer can enter the memo.
      results[anchor] = (byte) (found ? 2 : 1);
      return found;
    }

    private static boolean stereoTest(Expr expression) {
      return expression.type() == Expr.Type.STEREOCHEMISTRY
          || expression.type() == Expr.Type.OR && expression.left().type() == Expr.Type.STEREOCHEMISTRY;
    }

    IAtomContainer query(IAtomContainer original) {
      IAtomContainer copy = new QueryAtomContainer(SilentChemObjectBuilder.getInstance());
      copy.setProperties(new HashMap<>(original.getProperties()));
      Map<IChemObject, IChemObject> objects = new IdentityHashMap<>();
      for (IAtom atom : original.atoms()) {
        QueryAtom source = (QueryAtom) AtomRef.deref(atom);
        Expr expression = expression(source.getExpression());
        // Ring sizes and ring counts are read from the caches of the target. An atom that reads
        // ring sizes negates plainly, and the others as CDK does.
        boolean ringSizes = reads(expression, Expr.Type.RING_SMALLEST);
        boolean cached = ringSizes || reads(expression, Expr.Type.RING_COUNT)
            || reads(expression, Expr.Type.RECURSIVE);
        final int queryAtom = original == rootQuery ? copy.getAtomCount() : -1;
        QueryAtom bounded = new QueryAtom(expression) {
          @Override public boolean matches(IAtom target) {
            // Atom tests are short, so the clock and cancellation are read on every sixteenth.
            if ((++polls & 15) == 0) check(control);
            if (nodeLimit > 0 && nodes >= nodeLimit) { control.limitReached(); throw new CommonBondGraphs.Stop(); }
            nodes++;
            boolean matches = cached ? evaluate(getExpression(), target, ringSizes) : super.matches(target);
            observedDomain(queryAtom, target.getIndex(), matches);
            return matches;
          }
        };
        bounded.setProperties(new HashMap<>(atom.getProperties()));
        bounded.setFlags(atom.getFlags());
        copy.addAtom(bounded);
        objects.put(atom, copy.getAtom(copy.getAtomCount() - 1));
        objects.put(source, copy.getAtom(copy.getAtomCount() - 1));
      }
      for (IBond bond : original.bonds()) {
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
          // A negation has one operand; a retained right operand is not part of it.
          return new Expr(original.type(), expression(original.left()),
              original.type() == Expr.Type.NOT || original.right() == null ? null : expression(original.right()));
        default:
          return original;
      }
    }
  }

  /**
   * Return a mapped atom whose tetrahedral or double-bond ligand relation is incompatible. An
   * unspecified query centre or bond matches any target; an unspecified target centre or bond
   * matches a specified query one as well, as in a common-subgraph search.
   */
  static int stereoMismatchAtom(MolGraph query, MolGraph target, int[] mapping, ChemOptions chemistry) {
    return stereoMismatchAtom(query, target, mapping, chemistry, null, false);
  }

  /** The same test over the selected bonds alone, or over every bond when {@code selectedBonds} is null. */
  static int stereoMismatchAtom(MolGraph query, MolGraph target, int[] mapping, ChemOptions chemistry,
                                Set<SearchResult.Bond> selectedBonds) {
    return stereoMismatchAtom(query, target, mapping, chemistry, selectedBonds, false);
  }

  /**
   * The test of an embedding, as a substructure search makes it: a specified query centre or
   * double-bond configuration needs a specified, matching one in the target, as a SMARTS
   * chirality primitive or pair of directional bonds does. The native core applies the same rule.
   */
  static int stereoMismatchAtom(MolGraph query, MolGraph target, int[] mapping, ChemOptions chemistry,
                                boolean embedding) {
    return stereoMismatchAtom(query, target, mapping, chemistry, null, embedding);
  }

  private static int stereoMismatchAtom(MolGraph query, MolGraph target, int[] mapping, ChemOptions chemistry,
                                        Set<SearchResult.Bond> selectedBonds, boolean embedding) {
    if (chemistry.useChirality && query.tetraChirality != null) {
      for (int atom = 0; atom < query.n; atom++) {
        int mapped = mapping[atom];
        if (mapped < 0 || mapped >= target.n || query.tetraChirality[atom] == 0) continue;
        int targetTag = target.tetraChirality == null ? 0 : target.tetraChirality[mapped];
        if (targetTag == 0) {
          if (embedding) return atom;
          continue;
        }
        int[] queryOrder = query.tetrahedralLigands(atom), targetOrder = target.tetrahedralLigands(mapped);
        if (queryOrder == null || targetOrder == null || queryOrder.length != 4 || targetOrder.length != 4) continue;
        int[] permutation = new int[4];
        Arrays.fill(permutation, -1);
        boolean[] used = new boolean[4];
        int mappedCount = 0;
        for (int slot = 0; slot < 4; slot++) {
          int ligand = queryOrder[slot];
          if (ligand < 0 || mapping[ligand] < 0
              || (selectedBonds != null && !selectedBonds.contains(new SearchResult.Bond(atom, ligand)))) continue;
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
        int queryTag = query.tetraChirality[atom];
        if ((queryTag == targetTag) != ((inversions & 1) == 0)) return atom;
      }
    }
    if (chemistry.useBondStereo) {
      for (int atom = 0; atom < query.n; atom++) {
        if (mapping[atom] < 0) continue;
        for (int other : query.neighbors[atom])
          if (other > atom && (selectedBonds == null || selectedBonds.contains(new SearchResult.Bond(atom, other)))
              && !doubleBondParityCompatible(query, atom, other, target, mapping, selectedBonds, embedding)) return atom;
      }
    }
    return -1;
  }

  /**
   * Compare the stored double-bond relation against whichever query ligands are mapped. With
   * {@code embedding}, a configured query bond needs a configured target bond.
   */
  private static boolean doubleBondParityCompatible(MolGraph query, int first, int second, MolGraph target,
                                                    int[] mapping, Set<SearchResult.Bond> selectedBonds,
                                                    boolean embedding) {
    int mappedFirst = mapping[first], mappedSecond = mapping[second];
    if (mappedFirst < 0 || mappedSecond < 0) return true;
    int queryConfiguration = query.dbStereo(first, second);
    if (queryConfiguration == 0) return true;
    int targetConfiguration = target.dbStereo(mappedFirst, mappedSecond);
    if (targetConfiguration == 0) return !embedding;
    for (int end = 0; end < 2; end++) {
      int queryEnd = end == 0 ? first : second, otherQueryEnd = end == 0 ? second : first;
      int targetEnd = mapping[queryEnd], otherTargetEnd = mapping[otherQueryEnd];
      int queryReference = query.doubleBondReference(queryEnd, otherQueryEnd);
      int targetReference = target.doubleBondReference(targetEnd, otherTargetEnd);
      int queryLigand = -1;
      for (int neighbor : query.neighbors[queryEnd])
        if (neighbor != otherQueryEnd && mapping[neighbor] >= 0
            && (selectedBonds == null || selectedBonds.contains(new SearchResult.Bond(queryEnd, neighbor)))) {
          queryLigand = neighbor;
          break;
        }
      if (queryLigand < 0 || queryReference < 0 || targetReference < 0) return true;
      if (!target.hasBond(targetEnd, mapping[queryLigand])) return false;
      if ((queryLigand != queryReference) != (mapping[queryLigand] != targetReference))
        queryConfiguration = 3 - queryConfiguration;
    }
    return queryConfiguration == targetConfiguration;
  }

  /**
   * Query-side search data for one prepared query and chemistry. It is immutable once built, so
   * concurrent searches with the same prepared query share it safely.
   */
  static final class QueryTables {
    final int[][] neighborsByDegree, labels1, labels2, labels3;
    /**
     * The lowest-indexed atom that also has the same degree, ring membership and neighbourhood
     * labels, so that both pass the mapping-independent tests against the same targets.
     */
    final int[] rootClass;
    /**
     * The lowest-indexed atom of each class of query atoms with equal atom constraints, and the
     * class of each query atom. Atoms of one class have equal domains.
     */
    final int[] representatives, classIndex;
    /**
     * Constraint class by atomic number when the element is the only atom constraint, so a
     * single pass over the target builds every domain; otherwise null.
     */
    final int[] elementClass;
    /**
     * The slot of each atomic number in the atom count screen, or null when the query is not
     * screened: without the element test, with an atomic number outside 0 to 127, or with more
     * than 16 elements, the limits of the native core.
     */
    private final byte[] elementSlot;
    /** The query atoms of each slot by degree, degrees of 7 and more counted at 7. */
    private final int[] degreeCounts;

    @SuppressWarnings("deprecation")
    QueryTables(MolGraph query, ChemOptions chemistry) {
      byte[] slots = null;
      int[] counts = null;
      int used = 0;
      if (chemistry.matchAtomType) {
        slots = new byte[128];
        Arrays.fill(slots, (byte) -1);
        counts = new int[16 * 8];
        for (int atom = 0; atom < query.n; atom++) {
          int element = query.atomicNum[atom];
          if (element < 0 || element >= slots.length || (slots[element] < 0 && used == 16)) {
            slots = null;
            break;
          }
          if (slots[element] < 0) slots[element] = (byte) used++;
          counts[slots[element] * 8 + Math.min(query.degree[atom], 7)]++;
        }
      }
      elementSlot = slots;
      degreeCounts = slots == null ? null : Arrays.copyOf(counts, used * 8);
      boolean labels = chemistry.matchAtomType && !chemistry.tautomerAware;
      neighborsByDegree = query.getNeighborsByDegDesc();
      labels1 = labels ? matchingNLF(query.getNLF1(), chemistry) : null;
      labels2 = labels && chemistry.useTwoHopNLF && query.n > 12
          ? matchingNLF(query.getNLF2(), chemistry) : null;
      labels3 = labels && chemistry.useThreeHopNLF && query.n > 20
          ? matchingNLF(query.getNLF3(), chemistry) : null;
      // Ring counts take part in the atom constraints compared below.
      if (chemistry.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) query.ensureRingCounts();
      classIndex = new int[query.n];
      rootClass = new int[query.n];
      int[] classAtoms = new int[query.n], rootRepresentatives = new int[query.n];
      int classes = 0, rootClasses = 0, maxElement = 0;
      for (int atom = 0; atom < query.n; atom++) {
        int match = classes;
        for (int index = 0; index < classes; index++)
          if (MolGraph.ChemOps.sameAtomConstraints(query, classAtoms[index], atom, chemistry)) {
            match = index;
            break;
          }
        if (match == classes) classAtoms[classes++] = atom;
        classIndex[atom] = match;
        maxElement = Math.max(maxElement, query.atomicNum[atom]);
        int rootMatch = atom;
        for (int index = 0; index < rootClasses; index++)
          if (sameRootTests(query, rootRepresentatives[index], atom)) {
            rootMatch = rootRepresentatives[index];
            break;
          }
        if (rootMatch == atom) rootRepresentatives[rootClasses++] = atom;
        rootClass[atom] = rootMatch;
      }
      representatives = Arrays.copyOf(classAtoms, classes);
      boolean elementOnly = chemistry.matchAtomType && !chemistry.matchFormalCharge
          && chemistry.aromaticityMode != ChemOptions.AromaticityMode.STRICT && !chemistry.ringMatchesRingOnly
          && !chemistry.matchIsotope && !chemistry.useChirality
          && chemistry.ringFusionMode == ChemOptions.RingFusionMode.IGNORE;
      for (int atom = 0; atom < query.n && elementOnly; atom++) elementOnly = query.atomicNum[atom] >= 0;
      if (elementOnly) {
        // Here atomsCompatFast compares only atomic numbers, and each class is one element.
        elementClass = new int[maxElement + 1];
        Arrays.fill(elementClass, -1);
        for (int index = 0; index < classes; index++) elementClass[query.atomicNum[representatives[index]]] = index;
      } else {
        elementClass = null;
      }
    }

    /**
     * Whether the target has, for each element and each minimum degree, at least as many atoms
     * as the query. Under the element test an embedding maps query atoms one to one onto target
     * atoms of their element and of at least their degree, as atomsCompatFast and rootFeasible
     * require, so a shortfall proves that none exists. Keep the three in step.
     */
    boolean admits(MolGraph target) {
      if (elementSlot == null) return true;
      int[] available = new int[degreeCounts.length];
      for (int atom = 0; atom < target.n; atom++) {
        int element = target.atomicNum[atom];
        int slot = element >= 0 && element < elementSlot.length ? elementSlot[element] : -1;
        if (slot >= 0) available[slot * 8 + Math.min(target.degree[atom], 7)]++;
      }
      // Atoms of a degree also serve every lower one, so counts accumulate from the highest.
      for (int slot = 0; slot < degreeCounts.length; slot += 8) {
        int queryAtoms = 0, targetAtoms = 0;
        for (int degree = 7; degree >= 0; degree--) {
          queryAtoms += degreeCounts[slot + degree];
          targetAtoms += available[slot + degree];
          if (queryAtoms > targetAtoms) return false;
        }
      }
      return true;
    }

    private boolean sameRootTests(MolGraph query, int first, int second) {
      return classIndex[first] == classIndex[second]
          && query.degree[first] == query.degree[second] && query.ring[first] == query.ring[second]
          && (labels1 == null || Arrays.equals(labels1[first], labels1[second]))
          && (labels2 == null || Arrays.equals(labels2[first], labels2[second]))
          && (labels3 == null || Arrays.equals(labels3[first], labels3[second]));
    }

    private static int[][] matchingNLF(int[][] tables, ChemOptions chemistry) {
      int[][] result = new int[tables.length][];
      for (int atom = 0; atom < tables.length; atom++) result[atom] = SubgraphSearch.matchingNLF(tables[atom], chemistry);
      return result;
    }
  }

  static final class SubgraphSearch {
    final MolGraph gq, gt;
    final ChemOptions C;
    final SearchResult.Control tb;
    final int Nq, Nt;
    final int[] q2t, t2q;
    final int[][] qNLF1, tNLF1, qNLF2, tNLF2, qNLF3, tNLF3;
    final boolean useNLF, useTwoHop, useThreeHop, useBitParallel, useRingOnly, useStereo;
    final boolean bitParallelSufficient;
    final long[][] domain;
    /** Lazily built domain rows restricted by the tests that do not depend on the mapping. */
    final long[][] rootRows;
    final int[] rootClass;
    /** A query atom with no compatible target, or -1. */
    final int emptyDomainAtom;
    final long[] usedMask;
    final int tWords;
    final boolean singleWord;
    long nodesVisited = 0;
    boolean timedOut = false, found = false;
    boolean induced;
    long nodeLimit;
    int[] qdeg, tdeg;
    final int[][] qNeighborsByDegDesc;
    // Pre-allocated buffers to avoid per-call allocations in availableCandidates/iterateBitsInto
    final long[] availBuf;
    final int[] candBuf;
    // Per-depth candidate buffer stack: avoids Arrays.copyOf on every availableCandidates call
    final int[][] candBufByDepth;
    /** The work arrays of target neighbourhood labels: atoms seen by the current build hold its stamp. */
    private int[] labelSeen, labelReached;
    private int labelStamp;
    /** The fragment of each query atom while enumerating, and the number of fragments. */
    private int[] fragmentOf;
    private int fragmentCount;
    /** Ring coverage belongs to this invocation, never to the prepared query. */
    private CompleteRingCheck completeRings;

    private boolean completeRingEmbedding() {
      if (tb.stopped()) return false;
      if (completeRings == null) completeRings = new CompleteRingCheck(gq, gq.neighbors, gt, tb);
      return completeRings.holds(q2t);
    }

    @SuppressWarnings("deprecation")
    SubgraphSearch(MolGraph gq, QueryTables query, MolGraph gt, ChemOptions C, SearchResult.Control tb) {
      this.gq = gq; this.gt = gt;
      this.qdeg = gq.degree; this.tdeg = gt.degree;
      this.C = C; this.tb = tb;
      this.Nq = gq.n; this.Nt = gt.n;
      this.tWords = CommonBondGraphs.words(Nt);
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
      this.qNeighborsByDegDesc = query.neighborsByDegree;
      this.qNLF1 = useNLF ? query.labels1 : null;
      this.qNLF2 = useTwoHop ? query.labels2 : null;
      this.qNLF3 = useThreeHop ? query.labels3 : null;
      // Target labels are built per atom when a feasibility test first needs them.
      this.tNLF1 = useNLF ? new int[Nt][] : null;
      this.tNLF2 = useTwoHop ? new int[Nt][] : null;
      this.tNLF3 = useThreeHop ? new int[Nt][] : null;
      // Ring counts are worked out only when the options need them.
      if (C.ringFusionMode != ChemOptions.RingFusionMode.IGNORE) { gq.ensureRingCounts(); gt.ensureRingCounts(); }
      // Query atoms with equal atom constraints have equal domains, so each distinct row is built
      // once and shared; domain rows are never modified after construction.
      this.domain = new long[Nq][];
      int classes = query.representatives.length;
      long[][] classRows = query.elementClass != null ? elementRows(query.elementClass, classes) : new long[classes][];
      int empty = -1;
      for (int k = 0; k < classes; k++) {
        if (classRows[k] == null) classRows[k] = tb.stopped() ? new long[tWords] : domainRow(query.representatives[k]);
        if (empty < 0 && bitCount(classRows[k]) == 0) empty = query.representatives[k];
      }
      for (int i = 0; i < Nq; i++) domain[i] = classRows[query.classIndex[i]];
      this.emptyDomainAtom = empty;
      this.usedMask = new long[tWords];
      this.availBuf = new long[tWords];
      this.candBuf = new int[Nt];
      this.candBufByDepth = new int[Nq][];
      this.rootRows = new long[Nq][];
      this.rootClass = query.rootClass;
    }

    /** Domains of every element class in one pass over the target. */
    long[][] elementRows(int[] elementClass, int classes) {
      long[][] rows = new long[classes][tWords];
      for (int tj = 0; tj < Nt; tj++) {
        int element = gt.atomicNum[tj];
        int k = element >= 0 && element < elementClass.length ? elementClass[element] : -1;
        if (k >= 0) rows[k][tj >>> 6] |= 1L << (tj & 63);
      }
      return rows;
    }

    /** Target atoms compatible with one query atom before any mapping is made. */
    long[] domainRow(int qi) {
      long[] row = new long[tWords];
      // The element test repeats the first test of atomsCompatFast to skip most calls.
      boolean element = C.matchAtomType;
      int atomicNumber = gq.atomicNum[qi];
      for (int tj = 0; tj < Nt; tj++)
        if ((!element || gt.atomicNum[tj] == atomicNumber) && MolGraph.ChemOps.atomsCompatFast(gq, qi, gt, tj, C))
          row[tj >>> 6] |= 1L << (tj & 63);
      return row;
    }

    /** Project cached labels onto the constraints enforced by this search. */
    static int[] matchingNLF(int[] source, ChemOptions C) {
      if (C.aromaticityMode == ChemOptions.AromaticityMode.STRICT) return source;
      int[] merged = new int[source.length];
      int size = 0;
      for (int j = 0; j < source.length; j += 2) {
        int element = source[j] >>> 1;
        if (size > 0 && merged[size - 2] == element) merged[size - 1] += source[j + 1];
        else { merged[size++] = element; merged[size++] = source[j + 1]; }
      }
      return Arrays.copyOf(merged, size);
    }

    int[] targetNLF1(int tj) {
      int[] labels = tNLF1[tj];
      if (labels == null) tNLF1[tj] = labels = matchingNLF(MolGraph.buildNLF1(gt, tj), C);
      return labels;
    }

    int[] targetNLF2(int tj) {
      int[] labels = tNLF2[tj];
      if (labels == null) tNLF2[tj] = labels = matchingNLF(neighbourhoodLabels(tj, 2), C);
      return labels;
    }

    int[] targetNLF3(int tj) {
      int[] labels = tNLF3[tj];
      if (labels == null) tNLF3[tj] = labels = matchingNLF(neighbourhoodLabels(tj, 3), C);
      return labels;
    }

    /** The labels within {@code radius} bonds of a target atom, built in this search's work arrays. */
    int[] neighbourhoodLabels(int tj, int radius) {
      if (labelSeen == null) {
        labelSeen = new int[Nt];
        labelReached = new int[Nt];
      }
      if (++labelStamp == 0) {
        Arrays.fill(labelSeen, 0);
        labelStamp = 1;
      }
      return MolGraph.buildNeighborhoodNLF(gt, tj, radius, labelSeen, labelStamp, labelReached);
    }

    int[] candidateBuffer(int depth) {
      if (candBufByDepth[depth] == null) candBufByDepth[depth] = new int[Nt];
      return candBufByDepth[depth];
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

    /** Write the unused targets of qi that pass the tests independent of the mapping; return their count. */
    int domainCandidates(int qi, int[] buf) {
      long[] row = rootRow(qi);
      for (int w = 0; w < tWords; w++) availBuf[w] = row[w] & ~usedMask[w];
      return iterateBitsInto(availBuf, buf);
    }

    /** The first feasible unused target of qi closest in degree, or -1. */
    int greedyTarget(int qi) {
      int bestTj = -1, bestDist = Integer.MAX_VALUE;
      for (int candidate = 0, candidates = domainCandidates(qi, candBuf); candidate < candidates; candidate++) {
        int tj = candBuf[candidate];
        int dist = Math.abs(tdeg[tj] - qdeg[qi]);
        if (dist < bestDist && mappedFeasible(qi, tj)) { bestTj = tj; bestDist = dist; if (dist == 0) break; }
      }
      return bestTj;
    }

    boolean greedyProbe(int[] order) {
      int matched = 0;
      for (int qi : order) {
        if (stopped()) return false;
        nodesVisited++;
        int bestTj = greedyTarget(qi);
        if (bestTj < 0) {
          for (int k = matched - 1; k >= 0; k--) {
            int qk = order[k], tk = q2t[qk]; q2t[qk] = -1; t2q[tk] = -1;
            usedMask[tk >>> 6] &= ~(1L << (tk & 63));
          }
          return false;
        }
        q2t[qi] = bestTj; t2q[bestTj] = qi; usedMask[bestTj >>> 6] |= 1L << (bestTj & 63); matched++;
      }
      if (stereoMismatchAtom(gq, gt, q2t, C, true) < 0
          && (!C.completeRingsOnly || completeRingEmbedding())) return true;
      for (int qi : order) {
        int tj = q2t[qi];
        q2t[qi] = -1;
        t2q[tj] = -1;
        usedMask[tj >>> 6] &= ~(1L << (tj & 63));
      }
      return false;
    }

    /**
     * Order query atoms by fewest remaining targets, then highest degree, then lowest index,
     * narrowing each unpicked neighbour of a picked atom to targets adjacent to its candidates.
     */
    int[] orderQueryAtoms() {
      if (emptyDomainAtom >= 0 && !tb.stopped()) {
        // Every search stops at a first atom without targets, before reading the rest of the order.
        int[] order = new int[Nq];
        order[0] = emptyDomainAtom;
        for (int i = 0, next = 1; i < Nq; i++) if (i != emptyDomainAtom) order[next++] = i;
        return order;
      }
      int[] order = new int[Nq], size = new int[Nq];
      boolean[] picked = new boolean[Nq], narrowed = new boolean[Nq];
      // Rows are shared with the domain until they are first narrowed.
      long[][] workDomain = new long[Nq][];
      for (int i = 0; i < Nq; i++) { workDomain[i] = domain[i]; size[i] = bitCount(domain[i]); }
      long[] reachable = new long[tWords];
      for (int pos = 0; pos < Nq; pos++) {
        if (tb.stopped()) return order;
        int bestQ = nextQueryAtom(picked, size);
        order[pos] = bestQ; picked[bestQ] = true;
        if (pos == 0 && !hasRootCandidate(bestQ)) {
          // Every search stops at this first atom without reading the rest of the order.
          for (int i = 0, next = 1; i < Nq; i++) if (!picked[i]) order[next++] = i;
          return order;
        }
        narrowNeighbours(bestQ, picked, size, workDomain, narrowed, reachable);
      }
      return order;
    }

    /** Narrow the unpicked neighbours of a picked atom to targets adjacent to its candidates. */
    void narrowNeighbours(int bestQ, boolean[] picked, int[] size, long[][] workDomain, boolean[] narrowed,
                          long[] reachable) {
      boolean reachableReady = false;
      for (int qk : gq.neighbors[bestQ]) {
        // Skip picked neighbours and those whose domain is already empty.
        if (picked[qk] || size[qk] == 0) continue;
        if (!reachableReady) { reachableFrom(workDomain[bestQ], reachable); reachableReady = true; }
        if (!narrowed[qk]) { workDomain[qk] = workDomain[qk].clone(); narrowed[qk] = true; }
        long[] row = workDomain[qk];
        for (int rw = 0; rw < tWords; rw++) row[rw] &= reachable[rw];
        size[qk] = bitCount(row);
      }
    }

    /** Whether any domain target of a query atom is feasible while nothing is mapped. */
    boolean hasRootCandidate(int qi) {
      for (long word : rootRow(qi)) if (word != 0) return true;
      return false;
    }

    /** The unpicked query atom with the fewest targets, then the highest degree, then the lowest index. */
    int nextQueryAtom(boolean[] picked, int[] size) {
      int bestQ = -1, bestSize = Integer.MAX_VALUE, bestDeg = -1;
      for (int i = 0; i < Nq; i++) {
        if (picked[i]) continue;
        if (size[i] < bestSize || (size[i] == bestSize && qdeg[i] > bestDeg)) { bestQ = i; bestSize = size[i]; bestDeg = qdeg[i]; }
      }
      return bestQ;
    }

    /** Write the union of the target neighbourhoods of every atom in {@code row}. */
    void reachableFrom(long[] row, long[] reachable) {
      Arrays.fill(reachable, 0L);
      for (int w = 0; w < tWords; w++) {
        long bits = row[w];
        while (bits != 0) {
          int tj = (w << 6) | Long.numberOfTrailingZeros(bits);
          if (tj < Nt) for (int rw = 0; rw < tWords; rw++) reachable[rw] |= gt.adjLong[tj][rw];
          bits &= bits - 1;
        }
      }
    }

    static int bitCount(long[] row) {
      int count = 0;
      for (long word : row) count += Long.bitCount(word);
      return count;
    }

    /** Sort first {@code len} entries of {@code buf} by degree proximity to query atom qi. */
    void sortByDegreeProximityBuf(int[] buf, int len, int qi) {
      for (int i = 1; i < len; i++) {
        int key = buf[i], keyDist = Math.abs(tdeg[key] - qdeg[qi]); int j = i - 1;
        while (j >= 0 && Math.abs(tdeg[buf[j]] - qdeg[qi]) > keyDist) { buf[j + 1] = buf[j]; j--; }
        buf[j + 1] = key;
      }
    }

    /**
     * Domain targets of qi that pass the degree, ring and neighbourhood label tests. These tests
     * do not depend on the mapping, so each row is built once, when the search first needs it.
     */
    long[] rootRow(int qi) {
      long[] row = rootRows[qi];
      if (row != null) return row;
      // Atoms with equal mapping-independent tests share one read-only row.
      if (rootClass[qi] != qi) return rootRows[qi] = rootRow(rootClass[qi]);
      row = new long[tWords];
      long[] candidates = domain[qi];
      for (int w = 0; w < tWords; w++) {
        long bits = candidates[w];
        while (bits != 0) {
          int tj = (w << 6) | Long.numberOfTrailingZeros(bits);
          if (rootFeasible(qi, tj)) row[w] |= 1L << (tj & 63);
          bits &= bits - 1;
        }
      }
      return rootRows[qi] = row;
    }

    boolean rootFeasible(int qi, int tj) {
      if (gt.degree[tj] < gq.degree[qi]) return false;
      // The ring test is the cheapest, so it runs first.
      if (useRingOnly && gq.ring[qi] != gt.ring[tj]) return false;
      if (useNLF && !MolGraph.nlfOk(qNLF1[qi], targetNLF1(tj))) return false;
      if (useTwoHop && !MolGraph.nlfOk(qNLF2[qi], targetNLF2(tj))) return false;
      if (useThreeHop && !MolGraph.nlfOk(qNLF3[qi], targetNLF3(tj))) return false;
      return true;
    }

    /** The feasibility tests that depend on the atoms already mapped. */
    boolean mappedFeasible(int qi, int tj) {
      if (induced) {
        for (int mappedQuery = 0; mappedQuery < Nq; mappedQuery++)
          if (q2t[mappedQuery] >= 0 && !gq.hasBond(qi, mappedQuery)
              && gt.hasBond(tj, q2t[mappedQuery])) return false;
      }
      int[] nbSorted = qNeighborsByDegDesc[qi];
      if (useBitParallel) {
        long[] tAdj = gt.adjLong[tj];
        for (int qk : nbSorted) {
          int tk = q2t[qk];
          if (tk != -1 && (tAdj[tk >>> 6] & (1L << (tk & 63))) == 0) return false;
        }
        if (bitParallelSufficient) return true;
      }
      for (int qk : nbSorted) {
        int tk = q2t[qk];
        if (tk != -1 && !MolGraph.ChemOps.bondsCompatible(gq, qi, qk, gt, tj, tk, C)) return false;
      }
      return true;
    }

    boolean stopped() {
      if (nodeLimit > 0 && nodesVisited >= nodeLimit) tb.limitReached();
      if (tb.stopped()) { timedOut = true; return true; }
      return false;
    }

    /**
     * Choose the fewest feasible targets after enforcing the current mapped neighbours, the first
     * in order among equals, and keep them in this depth's buffer ordered by degree proximity.
     */
    int selectNext(int[] order, int position) {
      if (tb.stopped()) { timedOut = true; return 0; }
      // Atoms are mapped in order positions, so nothing is mapped at the first position.
      if (position == 0) return selectRoot(order);
      boolean hasFrontier = false;
      for (int index = position; index < order.length && !hasFrontier; index++)
        hasFrontier = hasMappedNeighbor(order[index]);
      int[] buffer = candidateBuffer(position);
      int bestIndex = position, bestCount = Integer.MAX_VALUE;
      for (int index = position; index < order.length; index++) {
        if (tb.stoppedSampled()) { timedOut = true; return 0; }
        int queryAtom = order[index];
        if (hasFrontier && !hasMappedNeighbor(queryAtom)) continue;
        int compatible = feasibleCandidates(queryAtom, buffer, bestCount);
        if (compatible < 0) return 0;
        if (compatible < bestCount) {
          bestIndex = index;
          bestCount = compatible;
          System.arraycopy(buffer, 0, candBuf, 0, compatible);
          if (compatible == 0) break;
        }
      }
      int swap = order[position]; order[position] = order[bestIndex]; order[bestIndex] = swap;
      System.arraycopy(candBuf, 0, buffer, 0, bestCount);
      // Filtering before the stable sort keeps the order of sorting every candidate first.
      sortByDegreeProximityBuf(buffer, bestCount, order[position]);
      return bestCount;
    }

    /**
     * The first atom: nothing is mapped, so its candidates are the targets of its root row. Atoms
     * of one root class share that row and their degree, so after the first of a class in order
     * none can have fewer targets, and only the first can be chosen.
     */
    int selectRoot(int[] order) {
      boolean[] counted = new boolean[Nq];
      int bestIndex = 0, bestCount = Integer.MAX_VALUE;
      for (int index = 0; index < order.length; index++) {
        if (tb.stoppedSampled()) { timedOut = true; return 0; }
        int queryAtom = order[index];
        if (counted[rootClass[queryAtom]]) continue;
        counted[rootClass[queryAtom]] = true;
        int count = bitCount(rootRow(queryAtom));
        if (count < bestCount) {
          bestIndex = index;
          bestCount = count;
          if (count == 0) break;
        }
      }
      int swap = order[0]; order[0] = order[bestIndex]; order[bestIndex] = swap;
      int[] buffer = candidateBuffer(0);
      int candidates = availableCandidates(rootRow(order[0]), order[0], buffer);
      sortByDegreeProximityBuf(buffer, candidates, order[0]);
      return candidates;
    }

    /**
     * Write the feasible candidates of one query atom to {@code buffer} in target order, up to
     * {@code bound} of them, since an atom with as many cannot have the fewest; -1 when stopped.
     */
    int feasibleCandidates(int queryAtom, int[] buffer, int bound) {
      // The candidate row already applies the tests that do not depend on the mapping.
      int candidates = availableCandidates(rootRow(queryAtom), queryAtom, buffer);
      int compatible = 0;
      for (int candidate = 0; candidate < candidates && compatible < bound; candidate++) {
        if (tb.stoppedSampled()) { timedOut = true; return -1; }
        int targetAtom = buffer[candidate];
        if (mappedFeasible(queryAtom, targetAtom)) buffer[compatible++] = targetAtom;
      }
      return compatible;
    }

    boolean hasMappedNeighbor(int queryAtom) {
      for (int neighbor : gq.neighbors[queryAtom]) if (q2t[neighbor] >= 0) return true;
      return false;
    }

    void backtrack(int[] order, int pos) {
      if (stopped()) return;
      nodesVisited++; if (found) return;
      if (pos == order.length) {
        found = stereoMismatchAtom(gq, gt, q2t, C, true) < 0
            && (!C.completeRingsOnly || completeRingEmbedding());
        return;
      }
      int candCount = selectNext(order, pos);
      int qi = order[pos];
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        if (tb.stoppedSampled()) { timedOut = true; return; }
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63);
        backtrack(order, pos + 1);
        if (found || timedOut) return;
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63));
      }
    }

    void enumerateRec(int[] order, int pos, List<int[]> out, List<int[]> fragmentOrders, int maxSolutions) {
      if (stopped()) return;
      nodesVisited++; if (out.size() >= maxSolutions) return;
      if (pos == order.length) {
        if (stereoMismatchAtom(gq, gt, q2t, C, true) >= 0) return;
        if (C.completeRingsOnly && !completeRingEmbedding()) return;
        out.add(q2t.clone());
        if (fragmentOrders != null) fragmentOrders.add(fragmentOrder(order));
        return;
      }
      int candCount = selectNext(order, pos);
      int qi = order[pos];
      for (int k = 0; k < candCount; k++) {
        int tj = candBufByDepth[pos][k];
        if (tb.stoppedSampled()) { timedOut = true; return; }
        q2t[qi] = tj; t2q[tj] = qi; usedMask[tj >>> 6] |= 1L << (tj & 63);
        enumerateRec(order, pos + 1, out, fragmentOrders, maxSolutions);
        q2t[qi] = -1; t2q[tj] = -1; usedMask[tj >>> 6] &= ~(1L << (tj & 63));
        if (out.size() >= maxSolutions || timedOut) return;
      }
    }

    /** The query fragments in the order their atoms first appear in the mapping order. */
    int[] fragmentOrder(int[] order) {
      int[] listed = new int[fragmentCount];
      boolean[] seen = new boolean[fragmentCount];
      for (int index = 0, count = 0; index < order.length && count < fragmentCount; index++) {
        int fragment = fragmentOf[order[index]];
        if (!seen[fragment]) {
          seen[fragment] = true;
          listed[count++] = fragment;
        }
      }
      return listed;
    }

    public boolean exists() {
      if (tb.stopped()) return false;
      if (Nq == 0) return true;
      int[] order = orderQueryAtoms();
      if (tb.stopped()) return false;
      if (Nq <= 50 && greedyProbe(order)) { found = true; return true; }
      backtrack(order, 0); return found;
    }

    /**
     * Collect up to {@code maxSolutions} complete embeddings, each a copy of q2t. When
     * {@code fragmentOrders} is given, it receives for each embedding the fragments, numbered by
     * {@code fragmentOf}, in the order the search reached their first atoms.
     */
    void enumerate(int maxSolutions, List<int[]> out, int[] fragmentOf, List<int[]> fragmentOrders) {
      if (tb.stopped()) return;
      if (Nq == 0) { out.add(new int[0]); return; }
      this.fragmentOf = fragmentOf;
      for (int fragment : fragmentOf) fragmentCount = Math.max(fragmentCount, fragment + 1);
      enumerateRec(orderQueryAtoms(), 0, out, fragmentOrders, maxSolutions);
    }
    /** Write unused targets of {@code row} adjacent to every mapped neighbour of qi, in index order. */
    int availableCandidates(long[] row, int qi, int[] buf) {
      for (int w = 0; w < tWords; w++) availBuf[w] = row[w] & ~usedMask[w];
      for (int queryNeighbor : gq.neighbors[qi]) {
        int mappedTarget = q2t[queryNeighbor];
        if (mappedTarget < 0) continue;
        for (int word = 0; word < tWords; word++) availBuf[word] &= gt.adjLong[mappedTarget][word];
      }
      return iterateBitsInto(availBuf, buf);
    }

  }
}

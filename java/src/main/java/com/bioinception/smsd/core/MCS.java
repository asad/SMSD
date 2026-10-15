/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.function.BooleanSupplier;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.ITetrahedralChirality;

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
    return new MCS(copyMolecule(query), new Options());
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
    return match(copyMolecule(target));
  }
  public SearchResult matchAll(MolGraph target) {
    return matchAll(target, 1000);
  }
  public SearchResult matchAll(IAtomContainer target) {
    return matchAll(copyMolecule(target), 1000);
  }
  public SearchResult matchAll(IAtomContainer target, int limit) {
    return matchAll(copyMolecule(target), limit);
  }
  public SearchResult matchAll(MolGraph target, int limit) {
    return run(target, limit, true);
  }

  private static MolGraph copyMolecule(IAtomContainer molecule) {
    Objects.requireNonNull(molecule, "molecule");
    try {
      return new MolGraph(molecule.clone());
    } catch (CloneNotSupportedException failure) {
      throw new IllegalArgumentException("Cannot copy molecule", failure);
    }
  }

  private SearchResult run(MolGraph target, int limit, boolean all) {
    Objects.requireNonNull(target, "target");
    if (limit < 1)
      throw new IllegalArgumentException("limit must be positive");
    for (int targetAtom : options.excludedTargetAtoms)
      if (targetAtom >= target.n)
        throw new IllegalArgumentException("Excluded target atom index out of range");
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, query.n, target.n),
        options.cancellation);
    if (control.stopped())
      return new SearchResult(
          Collections.emptyList(), control.status(), false, false, 0.0, 0, control.elapsedMillis());
    return new Kernel(query, target, options, limit, all, control).search();
  }

  /** Search independent targets with this immutable prepared query. */
  public List<SearchResult> matchMany(List<MolGraph> targets) {
    Objects.requireNonNull(targets, "targets");
    List<SearchResult> results = new ArrayList<>();
    for (MolGraph target : targets) results.add(match(target));
    return Collections.unmodifiableList(results);
  }

  /** Assign queries successively without reusing previously selected target atoms. */
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
      if (result.status() == SearchResult.Status.Cancelled)
        break;
    }
    return Collections.unmodifiableList(results);
  }

  /** Chemistry-only atom similarity bound; graph constraints can only reduce it. */
  public static double similarityUpperBound(
      MolGraph query, MolGraph target, ChemOptions chemistry) {
    Objects.requireNonNull(query, "query");
    Objects.requireNonNull(target, "target");
    chemistry = chemistry == null ? new ChemOptions() : chemistry;
    if (chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
      query.ensureRingCounts();
      target.ensureRingCounts();
    }
    int[] assigned = new int[target.n];
    Arrays.fill(assigned, -1);
    int count = 0;
    for (int atom = 0; atom < query.n; atom++)
      if (augmentAtom(query, target, chemistry, atom, assigned, new boolean[target.n]))
        count++;
    int union = query.n + target.n - count;
    return union == 0 ? 0.0 : (double) count / union;
  }

  public static double similarityUpperBound(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry) {
    return similarityUpperBound(new MolGraph(query), new MolGraph(target), chemistry);
  }

  private static boolean augmentAtom(MolGraph query, MolGraph target, ChemOptions chemistry,
      int queryAtom, int[] assigned, boolean[] visited) {
    for (int targetAtom = 0; targetAtom < target.n; targetAtom++) {
      if (visited[targetAtom]
          || !Substructure.MoleculeMatcher.atomsCompatFast(
              query, queryAtom, target, targetAtom, chemistry))
        continue;
      visited[targetAtom] = true;
      if (assigned[targetAtom] < 0
          || augmentAtom(query, target, chemistry, assigned[targetAtom], assigned, visited)) {
        assigned[targetAtom] = queryAtom;
        return true;
      }
    }
    return false;
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
    return new SearchResult.Match(atoms, bonds, fragments(atoms, bonds.keySet()), atoms.size());
  }

  /** Validate selected atom and bond correspondences against the declared objective. */
  public static List<String> validateMapping(
      MolGraph query, MolGraph target, SearchResult.Match match, Options options) {
    Options checked = copy(options, query.n);
    List<String> errors = new ArrayList<>();
    Map<Integer, Integer> atoms = match.atomMapping();
    Set<Integer> used = new HashSet<>();
    for (Map.Entry<Integer, Integer> pair : atoms.entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (queryAtom < 0 || queryAtom >= query.n || targetAtom < 0 || targetAtom >= target.n) {
        errors.add("Atom index out of range");
        continue;
      }
      if (!used.add(targetAtom))
        errors.add("Atom mapping is not injective");
      if (checked.excludedTargetAtoms.contains(targetAtom))
        errors.add("Excluded target atom is mapped");
      if (!Substructure.MoleculeMatcher.atomsCompatFast(
              query, queryAtom, target, targetAtom, checked.chemistry))
        errors.add("Atom mismatch");
    }
    if (!errors.isEmpty())
      return errors;
    Set<Integer> incident = new HashSet<>();
    for (Map.Entry<SearchResult.Bond, SearchResult.Bond> pair : match.bondMapping().entrySet()) {
      SearchResult.Bond bond = pair.getKey(), mapped = pair.getValue();
      if (!atoms.containsKey(bond.first()) || !atoms.containsKey(bond.second())) {
        errors.add("Selected bond endpoint is unmapped");
        continue;
      }
      SearchResult.Bond expected =
          new SearchResult.Bond(atoms.get(bond.first()), atoms.get(bond.second()));
      if (!expected.equals(mapped) || !query.hasBond(bond.first(), bond.second())
          || !target.hasBond(mapped.first(), mapped.second())
          || !MolGraph.ChemOps.bondsCompatible(query, bond.first(), bond.second(), target,
              atoms.get(bond.first()), atoms.get(bond.second()), checked.chemistry))
        errors.add("Selected bond mismatch");
      incident.add(bond.first());
      incident.add(bond.second());
    }
    if (checked.objective == Objective.MCIS) {
      for (int first : atoms.keySet())
        for (int second : atoms.keySet())
          if (first < second) {
            boolean queryBond = query.hasBond(first, second),
                    targetBond = target.hasBond(atoms.get(first), atoms.get(second));
            if (queryBond != targetBond)
              errors.add("Induced bond mismatch");
            if (queryBond && !match.bondMapping().containsKey(new SearchResult.Bond(first, second)))
              errors.add("Induced bond is not selected");
          }
    } else if (!incident.equals(atoms.keySet()))
      errors.add("MCES cannot contain isolated atoms");
    List<Set<Integer>> components = fragments(atoms, match.bondMapping().keySet());
    if (checked.connectivity == Connectivity.Connected && components.size() > 1)
      errors.add("Common graph is disconnected");
    if (components.size() > checked.maxFragments)
      errors.add("Too many common fragments");
    for (Set<Integer> component : components)
      if (component.size() < checked.minFragmentAtoms)
        errors.add("Common fragment is too small");
    if (checked.chemistry.completeRingsOnly
        && (!completeRings(query, atoms.keySet(), match.bondMapping().keySet())
            || !completeRings(target, new HashSet<>(atoms.values()),
                new HashSet<>(match.bondMapping().values()))))
      errors.add("Incomplete retained ring");
    int[] indexed = new int[query.n];
    Arrays.fill(indexed, -1);
    atoms.forEach((queryAtom, targetAtom) -> indexed[queryAtom] = targetAtom);
    if (Substructure.stereoMismatchAtom(
            query, target, indexed, checked.chemistry, match.bondMapping().keySet())
        >= 0)
      errors.add("Mapped stereochemistry mismatch");
    return errors;
  }

  /** Prove whether any valid extension retains every supplied correspondence and selected bond. */
  public static boolean isMappingMaximal(
      MolGraph query, MolGraph target, SearchResult.Match match, Options options) {
    if (!validateMapping(query, target, match, options).isEmpty())
      throw new IllegalArgumentException("Invalid common graph mapping");
    Options extension = copy(options, query.n);
    extension.atomWeights = null;
    Kernel kernel = new Kernel(query, target, extension, 1, false);
    kernel.seed(match);
    SearchResult result = kernel.search();
    if (!result.optimal())
      throw new IllegalStateException("Mapping maximality was not proved: " + result.status());
    if (!result.found())
      return true;
    SearchResult.Match largest = result.matches().get(0);
    return largest.atomCount() == match.atomCount() && largest.bondCount() == match.bondCount();
  }

  public static boolean isMappingMaximal(
      MolGraph query, MolGraph target, Map<Integer, Integer> mapping, ChemOptions chemistry) {
    Options options = new Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    return isMappingMaximal(
        query, target, fragment(query, mapping, allBonds(query, mapping).keySet()), options);
  }

  /** A progressive consensus and its independently supported molecule mappings. */
  public static final class ConsensusResult {
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
      this.fragment = fragment;
      this.mappings = Collections.unmodifiableMap(new TreeMap<>(mappings));
      this.status = status;
      this.nodesVisited = nodes;
      this.elapsedMillis = elapsed;
    }
    public SearchResult.Match fragment() {
      return fragment;
    }
    public int sourceIndex() {
      return sourceIndex;
    }
    public MolGraph graph() {
      return fragment.toGraph(source);
    }
    public String toSmiles() {
      return fragment.toSmiles(source);
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

  /** Progressively reduce candidate seeds until at least ceil(threshold*n) inputs support them. */
  public static ConsensusResult consensus(
      List<MolGraph> molecules, Options options, double threshold) {
    Objects.requireNonNull(molecules, "molecules");
    Objects.requireNonNull(options, "options");
    if (molecules.isEmpty() || !Double.isFinite(threshold) || threshold <= 0 || threshold > 1)
      throw new IllegalArgumentException("Consensus requires molecules and a threshold in (0,1]");
    if (options.atomWeights != null)
      throw new IllegalArgumentException("Consensus has no single query atom-weight domain");
    int largest = 0;
    for (MolGraph graph : molecules)
      largest = Math.max(largest, Objects.requireNonNull(graph, "molecule").n);
    SearchResult.Control control = new SearchResult.Control(
        SearchResult.Control.resolveTimeout(options.timeoutMillis, largest, largest),
        options.cancellation);
    int required = (int) Math.ceil(threshold * molecules.size());
    MolGraph bestSource = molecules.get(0);
    int bestSourceIndex = 0;
    SearchResult.Match best = fragment(bestSource, Collections.emptyMap(), Collections.emptyList());
    Map<Integer, SearchResult.Match> bestMappings = Collections.emptyMap();
    long nodes = 0;
    for (int seedIndex = 0; seedIndex < molecules.size() && !control.stopped(); seedIndex++) {
      MolGraph source = molecules.get(seedIndex);
      Map<Integer, Integer> initial = new TreeMap<>();
      for (int atom = 0; atom < source.n; atom++) initial.put(atom, atom);
      Map<SearchResult.Bond, SearchResult.Bond> initialBonds = allBonds(source, initial);
      if (options.objective == Objective.MCES) {
        Set<Integer> incident = new HashSet<>();
        for (SearchResult.Bond bond : initialBonds.keySet()) {
          incident.add(bond.first());
          incident.add(bond.second());
        }
        initial.keySet().retainAll(incident);
      }
      SearchResult.Match current = fragment(source, initial, initialBonds.keySet());
      while (current.atomCount() > 0 && !control.stopped()) {
        MolGraph graph = current.toGraph(source);
        int[] originals =
            current.atomMapping().keySet().stream().mapToInt(Integer::intValue).toArray();
        Map<Integer, SearchResult.Match> supported = new TreeMap<>();
        SearchResult.Match reduction = null;
        for (int targetIndex = 0; targetIndex < molecules.size() && !control.stopped();
            targetIndex++) {
          Options pair = copy(options, graph.n);
          if (options.nodeLimit != 0) {
            if (nodes >= options.nodeLimit) {
              control.limitReached();
              break;
            }
            pair.nodeLimit = options.nodeLimit - nodes;
          }
          SearchResult result =
              new Kernel(graph, molecules.get(targetIndex), pair, 1, false, control).search();
          nodes += result.nodesVisited();
          if (!result.found())
            continue;
          SearchResult.Match match = result.matches().get(0);
          if (match.atomCount() == graph.n && match.bondCount() == current.bondCount())
            supported.put(targetIndex, match);
          else if (reduction == null || better(match, reduction, options.objective))
            reduction = match;
        }
        if (supported.size() >= required) {
          if (better(current, best, options.objective)) {
            best = current;
            bestSource = source;
            bestSourceIndex = seedIndex;
            bestMappings = supported;
          }
          break;
        }
        if (reduction == null || control.stopped())
          break;
        Map<Integer, Integer> reducedAtoms = new TreeMap<>();
        for (int atom : reduction.atomMapping().keySet())
          reducedAtoms.put(originals[atom], originals[atom]);
        List<SearchResult.Bond> reducedBonds = new ArrayList<>();
        for (SearchResult.Bond bond : reduction.bondMapping().keySet())
          reducedBonds.add(
              new SearchResult.Bond(originals[bond.first()], originals[bond.second()]));
        current = fragment(source, reducedAtoms, reducedBonds);
      }
    }
    return new ConsensusResult(bestSource, bestSourceIndex, best, bestMappings, control.status(),
        nodes, control.elapsedMillis());
  }

  private static boolean better(
      SearchResult.Match left, SearchResult.Match right, Objective objective) {
    int primary = objective == Objective.MCES
        ? Integer.compare(left.bondCount(), right.bondCount())
        : Integer.compare(left.atomCount(), right.atomCount());
    return primary > 0
        || (primary == 0
            && (objective == Objective.MCES ? left.atomCount() > right.atomCount()
                                            : left.bondCount() > right.bondCount()));
  }

  private static Options copy(Options source, int atoms) {
    Objects.requireNonNull(source, "options");
    Options result = new Options();
    result.chemistry = ChemOptions.copyOf(Objects.requireNonNull(source.chemistry, "chemistry"));
    result.objective = Objects.requireNonNull(source.objective, "objective");
    result.connectivity = Objects.requireNonNull(source.connectivity, "connectivity");
    result.atomWeights = source.atomWeights == null ? null : source.atomWeights.clone();
    if (result.atomWeights != null) {
      if (result.objective == Objective.MCES)
        throw new IllegalArgumentException("Atom weights require the MCIS objective");
      if (result.atomWeights.length < atoms)
        throw new IllegalArgumentException("Weights must cover every query atom");
      for (double weight : result.atomWeights)
        if (!Double.isFinite(weight))
          throw new IllegalArgumentException("Atom weights must be finite");
    }
    if (source.timeoutMillis < -1 || source.nodeLimit < 0 || source.minFragmentAtoms < 1
        || source.maxFragments < 1)
      throw new IllegalArgumentException("Invalid search limits");
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

  /** Convert a complete graph without changing its original atom indices. */
  static IAtomContainer toMolecule(MolGraph graph) {
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int atom = 0; atom < graph.n; atom++) atoms.put(atom, atom);
    return toMolecule(graph, atoms, allBonds(graph, atoms).keySet());
  }

  /** Materialise exactly the selected common edges and cap cut valences. */
  static IAtomContainer toMolecule(
      MolGraph g, Map<Integer, Integer> mapping, Set<SearchResult.Bond> selected) {
    Objects.requireNonNull(g, "query");
    Objects.requireNonNull(mapping, "mapping");
    int[] atoms = mapping.keySet().stream().mapToInt(Integer::intValue).sorted().toArray();
    Map<Integer, Integer> reindex = new HashMap<>();
    IAtomContainer fragment = new org.openscience.cdk.silent.AtomContainer();
    Map<org.openscience.cdk.interfaces.IChemObject, org.openscience.cdk.interfaces.IChemObject>
        objects = new IdentityHashMap<>();
    boolean[] allLigandsRetained = new boolean[atoms.length];
    for (int atom = 0; atom < atoms.length; atom++) {
      int original = atoms[atom];
      if (original < 0 || original >= g.n)
        throw new IndexOutOfBoundsException("MCS atom index out of range");
      reindex.put(original, atom);
      IAtom retained = new org.openscience.cdk.silent.Atom(g.atomicNum[original]);
      retained.setFormalCharge(g.formalCharge[original]);
      if (g.massNumber[original] != 0)
        retained.setMassNumber(g.massNumber[original]);
      retained.setImplicitHydrogenCount(g.hydrogenCount(original));
      fragment.addAtom(retained);
      if (g.mol != null)
        objects.put(g.mol.getAtom(original), fragment.getAtom(atom));
    }
    for (int atom = 0; atom < atoms.length; atom++) {
      int original = atoms[atom];
      allLigandsRetained[atom] = true;
      for (int neighbor : g.neighbors[original]) {
        Integer retainedNeighbor = reindex.get(neighbor);
        if (retainedNeighbor == null
            || !selected.contains(new SearchResult.Bond(original, neighbor))) {
          allLigandsRetained[atom] = false;
          continue;
        }
        if (atom >= retainedNeighbor)
          continue;
        int order = g.bondOrder(original, neighbor);
        IBond retained = new org.openscience.cdk.silent.Bond(fragment.getAtom(atom),
            fragment.getAtom(retainedNeighbor),
            order == 2       ? IBond.Order.DOUBLE
                : order == 3 ? IBond.Order.TRIPLE
                             : IBond.Order.SINGLE);
        retained.setIsAromatic(g.bondAromatic(original, neighbor));
        fragment.addBond(retained);
        if (g.mol != null)
          objects.put(g.mol.getBond(g.mol.getAtom(original), g.mol.getAtom(neighbor)),
              fragment.getBond(fragment.getBondCount() - 1));
      }
    }
    org.openscience.cdk.graph.Cycles.markRingAtomsAndBonds(fragment);
    for (IBond bond : fragment.bonds()) {
      boolean aromatic = bond.isAromatic() && bond.isInRing();
      bond.setIsAromatic(aromatic);
      if (aromatic)
        for (IAtom atom : bond.atoms()) atom.setIsAromatic(true);
    }
    for (int atom = 0; atom < atoms.length; atom++) {
      IAtom retained = fragment.getAtom(atom);
      int original = atoms[atom];
      if (allLigandsRetained[atom] && retained.isAromatic() == g.aromatic[original])
        continue;
      int hydrogens;
      if (retained.isAromatic() && (g.atomicNum[original] == 7 || g.atomicNum[original] == 15)) {
        hydrogens =
            g.hydrogenCount(original) + g.neighbors[original].length - retained.getBondCount();
      } else {
        int bondOrderSum = 0;
        for (IBond bond : retained.bonds())
          bondOrderSum += bond.isAromatic() ? 1 : bond.getOrder().numeric();
        hydrogens = MolGraph.defaultHydrogenCount(
            g.atomicNum[original], retained.isAromatic(), bondOrderSum, g.formalCharge[original]);
      }
      retained.setImplicitHydrogenCount(hydrogens);
    }
    try {
      org.openscience.cdk.aromaticity.Kekulization.kekulize(fragment);
    } catch (org.openscience.cdk.exception.CDKException invalidAromaticFragment) {
      for (IBond bond : fragment.bonds()) {
        if (bond.isAromatic())
          bond.setOrder(IBond.Order.SINGLE);
        bond.setIsAromatic(false);
      }
      for (IAtom atom : fragment.atoms()) {
        if (!atom.isAromatic())
          continue;
        atom.setIsAromatic(false);
        atom.setImplicitHydrogenCount(MolGraph.defaultHydrogenCount(atom.getAtomicNumber(), false,
            (int) fragment.getBondOrderSum(atom), atom.getFormalCharge()));
      }
    }
    if (g.mol != null) {
      for (IStereoElement<?, ?> stereo : g.mol.stereoElements()) {
        IStereoElement<?, ?> retained = stereo.mapStrict(objects);
        if (retained != null)
          fragment.addStereoElement(retained);
      }
    } else {
      for (int atom = 0; atom < atoms.length; atom++) {
        if (!allLigandsRetained[atom] || g.tetraChirality[atoms[atom]] == 0)
          continue;
        int[] ligands = g.tetrahedralLigands(atoms[atom]);
        if (ligands == null)
          continue;
        IAtom[] retainedLigands = new IAtom[4];
        for (int ligand = 0; ligand < 4; ligand++)
          retainedLigands[ligand] =
              fragment.getAtom(ligands[ligand] < 0 ? atom : reindex.get(ligands[ligand]));
        fragment.addStereoElement(new org.openscience.cdk.stereo.TetrahedralChirality(
            fragment.getAtom(atom), retainedLigands,
            g.tetraChirality[atoms[atom]] == 2 ? ITetrahedralChirality.Stereo.CLOCKWISE
                                               : ITetrahedralChirality.Stereo.ANTI_CLOCKWISE));
      }
      for (IBond bond : fragment.bonds()) {
        int first = fragment.indexOf(bond.getBegin()), second = fragment.indexOf(bond.getEnd());
        int configuration = g.dbStereo(atoms[first], atoms[second]);
        if (configuration == 0 || !allLigandsRetained[first] || !allLigandsRetained[second])
          continue;
        IBond[] references = new IBond[2];
        for (int side = 0; side < 2; side++) {
          int end = side == 0 ? first : second, other = side == 0 ? second : first;
          for (int neighbor : g.neighbors[atoms[end]]) {
            if (neighbor == atoms[other])
              continue;
            references[side] =
                fragment.getBond(fragment.getAtom(end), fragment.getAtom(reindex.get(neighbor)));
            break;
          }
        }
        if (references[0] != null && references[1] != null)
          fragment.addStereoElement(
              new org.openscience.cdk.stereo.DoubleBondStereochemistry(bond, references,
                  configuration == 1 ? IDoubleBondStereochemistry.Conformation.TOGETHER
                                     : IDoubleBondStereochemistry.Conformation.OPPOSITE));
      }
    }
    return fragment;
  }

  /** Preserve source chemistry while selecting atoms and common edges for further matching. */
  static MolGraph toGraph(
      MolGraph graph, Map<Integer, Integer> atoms, Set<SearchResult.Bond> bonds) {
    int[] original = atoms.keySet().stream().mapToInt(Integer::intValue).sorted().toArray();
    Map<Integer, Integer> reindex = new HashMap<>();
    for (int index = 0; index < original.length; index++) reindex.put(original[index], index);
    int size = original.length;
    int[] elements = new int[size], charges = new int[size], masses = new int[size],
          hydrogens = new int[size], tetra = new int[size];
    boolean[] rings = new boolean[size], aromatic = new boolean[size];
    int[][] adjacency = new int[size][], orders = new int[size][], stereo = new int[size][];
    boolean[][] ringBonds = new boolean[size][], aromaticBonds = new boolean[size][];
    boolean hasDoubleBondStereo = false;
    for (int index = 0; index < size; index++) {
      int source = original[index];
      elements[index] = graph.atomicNum[source];
      charges[index] = graph.formalCharge[source];
      masses[index] = graph.massNumber[source];
      hydrogens[index] = graph.hydrogenCount(source);
      rings[index] = graph.ring[source];
      aromatic[index] = graph.aromatic[source];
      List<Integer> neighbors = new ArrayList<>();
      for (int neighbor : graph.neighbors[source])
        if (reindex.containsKey(neighbor)
            && bonds.contains(new SearchResult.Bond(source, neighbor)))
          neighbors.add(neighbor);
      int count = neighbors.size();
      adjacency[index] = new int[count];
      orders[index] = new int[count];
      ringBonds[index] = new boolean[count];
      aromaticBonds[index] = new boolean[count];
      for (int edge = 0; edge < count; edge++) {
        int neighbor = neighbors.get(edge);
        adjacency[index][edge] = reindex.get(neighbor);
        orders[index][edge] = graph.bondOrder(source, neighbor);
        ringBonds[index][edge] = graph.bondInRing(source, neighbor);
        aromaticBonds[index][edge] = graph.bondAromatic(source, neighbor);
        boolean ligandsRetained = count == graph.neighbors[source].length;
        for (int ligand : graph.neighbors[neighbor])
          if (!reindex.containsKey(ligand)
              || !bonds.contains(new SearchResult.Bond(neighbor, ligand)))
            ligandsRetained = false;
        int configuration = ligandsRetained ? graph.dbStereo(source, neighbor) : 0;
        if (configuration != 0) {
          if (stereo[index] == null)
            stereo[index] = new int[size];
          stereo[index][reindex.get(neighbor)] = configuration;
          hasDoubleBondStereo = true;
        }
      }
      int[] reference = graph.tetrahedralLigands(source);
      if (graph.tetraChirality != null && graph.tetraChirality[source] != 0 && reference != null
          && count == graph.neighbors[source].length) {
        List<Integer> outputOrder = new ArrayList<>();
        if (count == 3)
          outputOrder.add(-1);
        outputOrder.addAll(neighbors);
        boolean odd = false;
        for (int first = 0; first < reference.length; first++)
          for (int second = first + 1; second < reference.length; second++)
            if (outputOrder.indexOf(reference[first]) > outputOrder.indexOf(reference[second]))
              odd = !odd;
        int tag =
            graph.mol == null ? graph.tetraChirality[source] : 3 - graph.tetraChirality[source];
        tetra[index] = odd ? 3 - tag : tag;
      }
    }
    MolGraph view = new MolGraph.Builder()
                        .atomCount(size)
                        .atomicNumbers(elements)
                        .formalCharges(charges)
                        .massNumbers(masses)
                        .hydrogenCounts(hydrogens)
                        .ringFlags(rings)
                        .aromaticFlags(aromatic)
                        .neighbors(adjacency)
                        .bondOrders(orders)
                        .bondRingFlags(ringBonds)
                        .bondAromaticFlags(aromaticBonds)
                        .tetrahedralChirality(tetra)
                        .doubleBondStereo(hasDoubleBondStereo ? stereo : null)
                        .build();
    graph.ensureRingCounts();
    view.ensureRingCounts();
    graph.ensureTautomerClasses();
    view.ensureTautomerClasses();
    for (int index = 0; index < size; index++) {
      view.ringCount[index] = graph.ringCount[original[index]];
      view.tautomerClass[index] = graph.tautomerClass[original[index]];
      view.tautomerWeight[index] = graph.tautomerWeight[original[index]];
    }
    return view;
  }

  static Map<SearchResult.Bond, SearchResult.Bond> allBonds(
      MolGraph graph, Map<Integer, Integer> atoms) {
    Map<SearchResult.Bond, SearchResult.Bond> result = new LinkedHashMap<>();
    for (int queryAtom : atoms.keySet())
      for (int neighbor : graph.neighbors[queryAtom])
        if (queryAtom < neighbor && atoms.containsKey(neighbor))
          result.put(new SearchResult.Bond(queryAtom, neighbor),
              new SearchResult.Bond(atoms.get(queryAtom), atoms.get(neighbor)));
    return result;
  }

  static List<Set<Integer>> fragments(
      Map<Integer, Integer> atoms, Collection<SearchResult.Bond> bonds) {
    Map<Integer, List<Integer>> adjacency = new TreeMap<>();
    for (int atom : atoms.keySet()) adjacency.put(atom, new ArrayList<>());
    for (SearchResult.Bond bond : bonds) {
      adjacency.get(bond.first()).add(bond.second());
      adjacency.get(bond.second()).add(bond.first());
    }
    Set<Integer> visited = new HashSet<>();
    List<Set<Integer>> result = new ArrayList<>();
    for (int root : adjacency.keySet()) {
      if (!visited.add(root))
        continue;
      Set<Integer> component = new TreeSet<>();
      ArrayDeque<Integer> queue = new ArrayDeque<>();
      queue.add(root);
      while (!queue.isEmpty()) {
        int atom = queue.remove();
        component.add(atom);
        for (int neighbor : adjacency.get(atom))
          if (visited.add(neighbor))
            queue.add(neighbor);
      }
      result.add(component);
    }
    return result;
  }

  /** Per-invocation state; prepared queries do not retain mutable search data. */
  private static final class Kernel {
    private final MolGraph query, target;
    private final Options options;
    private final SearchResult.Control control;
    private final boolean all;
    private final int limit;
    private final int[] queryToTarget, targetToQuery;
    private final int[][] domains;
    private final long[][] domainBits;
    private final long[] usedTargetBits, availableTargetBits;
    private final boolean[] availableQuery, availableTarget, reachableQuery, reachableTarget;
    private final int[] queryQueue, targetQueue;
    private final int[] mappedQueryAtoms;
    private int mappedCount;
    private final double[] weights;
    private final List<SearchResult.Bond> queryEdges, targetEdges;
    private final boolean[] selectedQueryEdges, selectedTargetEdges, excludedEdges;
    private final boolean[] excludedAtoms;
    private final List<SearchResult.Match> best = new ArrayList<>();
    private final Set<String> bestKeys = new HashSet<>();
    private double bestWeight;
    private long nodes;
    private int bestAtoms, bestBonds, initialAtoms, initialBonds;
    private boolean overflow, upperBoundReached;

    Kernel(MolGraph query, MolGraph target, Options options, int limit, boolean all) {
      this(query, target, options, limit, all, null);
    }

    Kernel(MolGraph query, MolGraph target, Options options, int limit, boolean all,
        SearchResult.Control sharedControl) {
      this.query = query;
      this.target = target;
      this.options = options;
      this.limit = limit;
      this.all = all;
      control = sharedControl != null
          ? sharedControl
          : new SearchResult.Control(
                SearchResult.Control.resolveTimeout(options.timeoutMillis, query.n, target.n),
                options.cancellation);
      queryToTarget = new int[query.n];
      Arrays.fill(queryToTarget, -1);
      targetToQuery = new int[target.n];
      Arrays.fill(targetToQuery, -1);
      excludedAtoms = new boolean[query.n];
      domainBits = new long[query.n][target.words];
      usedTargetBits = new long[target.words];
      availableTargetBits = new long[target.words];
      availableQuery = new boolean[query.n];
      availableTarget = new boolean[target.n];
      reachableQuery = new boolean[query.n];
      reachableTarget = new boolean[target.n];
      queryQueue = new int[query.n];
      targetQueue = new int[target.n];
      mappedQueryAtoms = new int[query.n];
      weights = new double[query.n];
      for (int atom = 0; atom < query.n; atom++)
        weights[atom] = options.atomWeights == null ? 1.0 : options.atomWeights[atom];
      if (control.stopped()) {
        domains = new int[query.n][];
        queryEdges = Collections.emptyList();
        targetEdges = Collections.emptyList();
        selectedQueryEdges = new boolean[0];
        selectedTargetEdges = new boolean[0];
        excludedEdges = new boolean[0];
        return;
      }
      if (options.chemistry.ringFusionMode == ChemOptions.RingFusionMode.STRICT) {
        query.ensureRingCounts();
        target.ensureRingCounts();
      }
      if (options.chemistry.tautomerAware) {
        query.ensureTautomerClasses();
        target.ensureTautomerClasses();
      }
      domains = new int[query.n][];
      Arrays.fill(domains, new int[0]);
      for (int queryAtom = 0; queryAtom < query.n && !control.stopped(); queryAtom++) {
        int[] candidates = new int[target.n];
        int count = 0;
        for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
          if (!options.excludedTargetAtoms.contains(targetAtom)
              && Substructure.MoleculeMatcher.atomsCompatFast(
                  query, queryAtom, target, targetAtom, options.chemistry))
            candidates[count++] = targetAtom;
        for (int index = 0; index < count; index++) {
          int targetAtom = candidates[index];
          domainBits[queryAtom][targetAtom >>> 6] |= 1L << (targetAtom & 63);
        }
        domains[queryAtom] = Arrays.copyOf(candidates, count);
      }
      queryEdges = edges(query);
      targetEdges = edges(target);
      selectedQueryEdges = new boolean[queryEdges.size()];
      excludedEdges = new boolean[queryEdges.size()];
      selectedTargetEdges = new boolean[targetEdges.size()];
    }

    SearchResult search() {
      if (!control.stopped()) {
        if (options.objective == Objective.MCIS)
          atoms(initialAtoms, initialBonds);
        else
          edges(0, initialAtoms, initialBonds);
      }
      boolean proved = control.status() == SearchResult.Status.Completed;
      if (proved && overflow && all)
        control.limitReached();
      double score = best.isEmpty() ? 0.0 : best.get(0).score();
      return new SearchResult(best, control.status(), proved, proved && !(overflow && all), score,
          nodes, control.elapsedMillis());
    }

    void seed(SearchResult.Match match) {
      match.atomMapping().forEach((queryAtom, targetAtom) -> {
        queryToTarget[queryAtom] = targetAtom;
        targetToQuery[targetAtom] = queryAtom;
        mappedQueryAtoms[mappedCount++] = queryAtom;
        usedTargetBits[targetAtom >>> 6] |= 1L << (targetAtom & 63);
      });
      initialAtoms = match.atomCount();
      initialBonds = match.bondCount();
      for (SearchResult.Bond bond : match.bondMapping().keySet())
        selectedQueryEdges[queryEdges.indexOf(bond)] = true;
      for (SearchResult.Bond bond : match.bondMapping().values())
        selectedTargetEdges[targetEdges.indexOf(bond)] = true;
    }

    private boolean visit() {
      if (upperBoundReached || control.stopped())
        return false;
      if (options.nodeLimit != 0 && nodes >= options.nodeLimit) {
        control.limitReached();
        return false;
      }
      nodes++;
      return true;
    }

    private static List<SearchResult.Bond> edges(MolGraph graph) {
      List<SearchResult.Bond> result = new ArrayList<>();
      for (int atom = 0; atom < graph.n; atom++)
        for (int neighbor : graph.neighbors[atom])
          if (atom < neighbor)
            result.add(new SearchResult.Bond(atom, neighbor));
      return result;
    }

    private long feasibleTargets(int queryAtom, int word) {
      long candidates = domainBits[queryAtom][word] & ~usedTargetBits[word];
      boolean hasMappedNeighbor = false;
      for (int index = 0; index < mappedCount && candidates != 0; index++) {
        int mappedQuery = mappedQueryAtoms[index];
        long adjacent = target.adjLong[queryToTarget[mappedQuery]][word];
        if (query.bondOrder(queryAtom, mappedQuery) != 0) {
          candidates &= adjacent;
          hasMappedNeighbor = true;
        } else {
          candidates &= ~adjacent;
        }
      }
      if (!hasMappedNeighbor)
        return candidates;
      long remaining = candidates;
      while (remaining != 0) {
        int offset = Long.numberOfTrailingZeros(remaining);
        int targetAtom = (word << 6) + offset;
        remaining &= remaining - 1;
        for (int neighbor : query.neighbors[queryAtom]) {
          int mappedTarget = queryToTarget[neighbor];
          if (mappedTarget >= 0 && !MolGraph.ChemOps.bondsCompatible(
                  query, queryAtom, neighbor, target, targetAtom, mappedTarget, options.chemistry)) {
            candidates &= ~(1L << offset);
            break;
          }
        }
      }
      return candidates;
    }

    private void atoms(int mappedAtoms, int mappedBonds) {
      if (!visit())
        return;
      int remaining = 0;
      Arrays.fill(availableTargetBits, 0L);
      for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
        availableQuery[queryAtom] = queryToTarget[queryAtom] >= 0;

      int chosen = -1, smallestDomain = Integer.MAX_VALUE, mostNeighbors = -1;
      for (int queryAtom = 0; queryAtom < query.n; queryAtom++) {
        if (excludedAtoms[queryAtom] || queryToTarget[queryAtom] >= 0)
          continue;
        int count = 0;
        for (int word = 0; word < target.words; word++) {
          long candidates = feasibleTargets(queryAtom, word);
          availableTargetBits[word] |= candidates;
          count += Long.bitCount(candidates);
        }
        if (count > 0) {
          remaining++;
          availableQuery[queryAtom] = true;
        }
        int mappedNeighbors = 0;
        for (int neighbor : query.neighbors[queryAtom])
          if (queryToTarget[neighbor] >= 0)
            mappedNeighbors++;
        if (options.connectivity == Connectivity.Connected && mappedAtoms > 0
            && mappedNeighbors == 0)
          continue;
        if (count < smallestDomain
            || (count == smallestDomain && mappedNeighbors > mostNeighbors)) {
          chosen = queryAtom;
          smallestDomain = count;
          mostNeighbors = mappedNeighbors;
        }
      }
      int targetCapacity = 0;
      for (long candidates : availableTargetBits)
        targetCapacity += Long.bitCount(candidates);
      int atomBound = mappedAtoms + Math.min(remaining, targetCapacity);
      for (int targetAtom = 0; targetAtom < target.n; targetAtom++)
        availableTarget[targetAtom] = targetToQuery[targetAtom] >= 0
            || (availableTargetBits[targetAtom >>> 6] & (1L << (targetAtom & 63))) != 0;
      // Feasible domains only shrink; every connected extension stays in these components.
      if (options.connectivity == Connectivity.Connected) {
        int queryRoot = mappedCount == 0 ? -1 : mappedQueryAtoms[0];
        int targetRoot = mappedCount == 0 ? -1 : queryToTarget[queryRoot];
        atomBound = Math.min(atomBound,
            connectedBound(query, availableQuery, reachableQuery, queryQueue, queryRoot));
        atomBound = Math.min(atomBound,
            connectedBound(target, availableTarget, reachableTarget, targetQueue, targetRoot));
      }
      int queryBondBound = 0, targetBondBound = 0;
      for (int edge = 0; edge < queryEdges.size(); edge++) {
        SearchResult.Bond bond = queryEdges.get(edge);
        if (availableQuery[bond.first()] && availableQuery[bond.second()]
            && (options.connectivity != Connectivity.Connected
                || (reachableQuery[bond.first()] && reachableQuery[bond.second()])))
          queryBondBound++;
      }
      for (int edge = 0; edge < targetEdges.size(); edge++) {
        SearchResult.Bond bond = targetEdges.get(edge);
        if (availableTarget[bond.first()] && availableTarget[bond.second()]
            && (options.connectivity != Connectivity.Connected
                || (reachableTarget[bond.first()] && reachableTarget[bond.second()])))
          targetBondBound++;
      }
      int bondBound = (int) Math.min(Math.min(queryBondBound, targetBondBound),
          (long) atomBound * (atomBound - 1) / 2);
      double weightBound = options.atomWeights == null ? 0.0 : weightBound();
      if (!Double.isNaN(weightBound)) {
        int boundOrder = compare(atomBound, bondBound, weightBound);
        if (boundOrder < 0 || (!all && boundOrder == 0 && !best.isEmpty()))
          return;
      }
      consider(mappedAtoms, mappedBonds,
          options.atomWeights == null ? mappedAtoms : orderedWeight());
      if (chosen < 0)
        return;
      int addedBonds = 0;
      for (int neighbor : query.neighbors[chosen])
        if (queryToTarget[neighbor] >= 0)
          addedBonds++;
      for (int word = 0; word < target.words; word++) {
        long candidates = feasibleTargets(chosen, word);
        while (candidates != 0) {
          int targetAtom = (word << 6) + Long.numberOfTrailingZeros(candidates);
          candidates &= candidates - 1;
          queryToTarget[chosen] = targetAtom;
          targetToQuery[targetAtom] = chosen;
          usedTargetBits[word] |= 1L << (targetAtom & 63);
          mappedQueryAtoms[mappedCount++] = chosen;
          atoms(mappedAtoms + 1, mappedBonds + addedBonds);
          mappedCount--;
          usedTargetBits[word] &= ~(1L << (targetAtom & 63));
          targetToQuery[targetAtom] = -1;
          queryToTarget[chosen] = -1;
          if (upperBoundReached || control.stopped())
            return;
        }
      }
      excludedAtoms[chosen] = true;
      atoms(mappedAtoms, mappedBonds);
      excludedAtoms[chosen] = false;
    }

    private static int connectedBound(
        MolGraph graph, boolean[] available, boolean[] reached, int[] queue, int root) {
      Arrays.fill(reached, false);
      int largest = 0;
      for (int start = root < 0 ? 0 : root; start < graph.n; start++) {
        if (!available[start] || reached[start])
          continue;
        int head = 0, tail = 0;
        queue[tail++] = start;
        reached[start] = true;
        while (head < tail) {
          int atom = queue[head++];
          for (int neighbor : graph.neighbors[atom])
            if (available[neighbor] && !reached[neighbor]) {
              reached[neighbor] = true;
              queue[tail++] = neighbor;
            }
        }
        if (root >= 0)
          return tail;
        largest = Math.max(largest, tail);
      }
      return largest;
    }

    private double orderedWeight() {
      double result = 0.0;
      for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
        if (queryToTarget[queryAtom] >= 0)
          result += weights[queryAtom];
      return result;
    }

    private double weightBound() {
      double result = 0.0;
      for (int queryAtom = 0; queryAtom < query.n; queryAtom++) {
        double contribution = queryToTarget[queryAtom] >= 0 ? weights[queryAtom]
            : excludedAtoms[queryAtom]                      ? 0.0
                                                            : Math.max(0.0, weights[queryAtom]);
        result = Math.nextUp(result + contribution);
      }
      return result;
    }

    private boolean endpointCompatible(int queryAtom, int targetAtom) {
      if (queryToTarget[queryAtom] >= 0)
        return queryToTarget[queryAtom] == targetAtom;
      if (targetToQuery[targetAtom] >= 0)
        return false;
      return Arrays.binarySearch(domains[queryAtom], targetAtom) >= 0;
    }

    private void edges(int depth, int mappedAtoms, int mappedBonds) {
      if (!visit())
        return;
      int remaining = 0, chosen = -1, smallestDomain = Integer.MAX_VALUE;
      for (int edgeIndex = 0; edgeIndex < queryEdges.size(); edgeIndex++) {
        if (selectedQueryEdges[edgeIndex] || excludedEdges[edgeIndex])
          continue;
        remaining++;
        SearchResult.Bond bond = queryEdges.get(edgeIndex);
        if (options.connectivity == Connectivity.Connected && mappedBonds > 0
            && queryToTarget[bond.first()] < 0 && queryToTarget[bond.second()] < 0)
          continue;
        int candidates = 0;
        for (int targetIndex = 0; targetIndex < targetEdges.size(); targetIndex++) {
          if (selectedTargetEdges[targetIndex])
            continue;
          SearchResult.Bond targetBond = targetEdges.get(targetIndex);
          for (int orientation = 0; orientation < 2; orientation++) {
            int first = orientation == 0 ? targetBond.first() : targetBond.second();
            int second = orientation == 0 ? targetBond.second() : targetBond.first();
            if (endpointCompatible(bond.first(), first) && endpointCompatible(bond.second(), second)
                && MolGraph.ChemOps.bondsCompatible(
                    query, bond.first(), bond.second(), target, first, second, options.chemistry))
              candidates++;
          }
        }
        if (candidates < smallestDomain) {
          chosen = edgeIndex;
          smallestDomain = candidates;
        }
      }
      int bondBound = mappedBonds + Math.min(remaining, targetEdges.size() - mappedBonds);
      if (bondBound < bestBonds
          || (bondBound == bestBonds && Math.min(query.n, target.n) < bestAtoms))
        return;
      consider(mappedAtoms, mappedBonds, 0L);
      if (chosen < 0)
        return;
      SearchResult.Bond queryBond = queryEdges.get(chosen);
      int queryFirst = queryBond.first(), querySecond = queryBond.second();
      for (int targetIndex = 0; targetIndex < targetEdges.size(); targetIndex++) {
        if (selectedTargetEdges[targetIndex])
          continue;
        SearchResult.Bond targetBond = targetEdges.get(targetIndex);
        for (int orientation = 0; orientation < 2; orientation++) {
          int targetFirst = orientation == 0 ? targetBond.first() : targetBond.second();
          int targetSecond = orientation == 0 ? targetBond.second() : targetBond.first();
          if (!endpointCompatible(queryFirst, targetFirst)
              || !endpointCompatible(querySecond, targetSecond)
              || !MolGraph.ChemOps.bondsCompatible(query, queryFirst, querySecond, target,
                  targetFirst, targetSecond, options.chemistry))
            continue;
          boolean firstNew = queryToTarget[queryFirst] < 0,
                  secondNew = queryToTarget[querySecond] < 0;
          queryToTarget[queryFirst] = targetFirst;
          targetToQuery[targetFirst] = queryFirst;
          queryToTarget[querySecond] = targetSecond;
          targetToQuery[targetSecond] = querySecond;
          selectedQueryEdges[chosen] = true;
          selectedTargetEdges[targetIndex] = true;
          edges(depth + 1, mappedAtoms + (firstNew ? 1 : 0) + (secondNew ? 1 : 0), mappedBonds + 1);
          selectedQueryEdges[chosen] = false;
          selectedTargetEdges[targetIndex] = false;
          if (firstNew) {
            queryToTarget[queryFirst] = -1;
            targetToQuery[targetFirst] = -1;
          }
          if (secondNew) {
            queryToTarget[querySecond] = -1;
            targetToQuery[targetSecond] = -1;
          }
          if (upperBoundReached || control.stopped())
            return;
        }
      }
      excludedEdges[chosen] = true;
      edges(depth + 1, mappedAtoms, mappedBonds);
      excludedEdges[chosen] = false;
    }

    private int compare(int atoms, int bonds, double weight) {
      if (options.objective == Objective.MCES) {
        int order = Integer.compare(bonds, bestBonds);
        return order != 0 ? order : Integer.compare(atoms, bestAtoms);
      }
      if (options.atomWeights != null) {
        int order = Double.compare(weight, bestWeight);
        if (order != 0)
          return order;
      }
      int order = Integer.compare(atoms, bestAtoms);
      return order != 0 ? order : Integer.compare(bonds, bestBonds);
    }

    private void consider(int atomCount, int bondCount, double weight) {
      int objectiveOrder = compare(atomCount, bondCount, weight);
      if (atomCount == 0 || objectiveOrder < 0 || (!all && objectiveOrder == 0 && !best.isEmpty()))
        return;
      Map<Integer, Integer> atoms = new TreeMap<>();
      for (int queryAtom = 0; queryAtom < query.n; queryAtom++)
        if (queryToTarget[queryAtom] >= 0)
          atoms.put(queryAtom, queryToTarget[queryAtom]);
      Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
      for (int index = 0; index < queryEdges.size(); index++) {
        SearchResult.Bond bond = queryEdges.get(index);
        if (options.objective == Objective.MCIS
                ? atoms.containsKey(bond.first()) && atoms.containsKey(bond.second())
                : selectedQueryEdges[index])
          bonds.put(bond,
              new SearchResult.Bond(queryToTarget[bond.first()], queryToTarget[bond.second()]));
      }
      if ((options.chemistry.useChirality || options.chemistry.useBondStereo)
          && Substructure.stereoMismatchAtom(query, target, queryToTarget, options.chemistry,
                 options.objective == Objective.MCES ? bonds.keySet() : null)
              >= 0)
        return;
      List<Set<Integer>> components = fragments(atoms, bonds.keySet());
      if (options.connectivity == Connectivity.Connected && components.size() != 1)
        return;
      if (components.size() > options.maxFragments)
        return;
      for (Set<Integer> component : components)
        if (component.size() < options.minFragmentAtoms)
          return;
      if (options.chemistry.completeRingsOnly
          && (!completeRings(query, atoms.keySet(), bonds.keySet())
              || !completeRings(
                  target, new HashSet<>(atoms.values()), new HashSet<>(bonds.values()))))
        return;
      int comparison = compare(atomCount, bondCount, weight);
      if (comparison > 0) {
        best.clear();
        bestKeys.clear();
        overflow = false;
        bestAtoms = atomCount;
        bestBonds = bondCount;
        bestWeight = weight;
      }
      String key = atoms.toString() + bonds.keySet();
      if (!bestKeys.add(key))
        return;
      if (best.size() >= limit) {
        if (all)
          overflow = true;
        return;
      }
      double score = options.objective == Objective.MCES ? bondCount
          : options.atomWeights == null                  ? atomCount
                                                         : weight;
      best.add(new SearchResult.Match(atoms, bonds, components, score));
      if (!all && options.atomWeights == null && atomCount == Math.min(query.n, target.n)
          && bondCount == Math.min(queryEdges.size(), targetEdges.size()))
        upperBoundReached = true;
    }
  }

  private static boolean completeRings(
      MolGraph graph, Set<Integer> selectedAtoms, Set<SearchResult.Bond> selectedBonds) {
    Set<Integer> coveredAtoms = new HashSet<>();
    Set<SearchResult.Bond> coveredBonds = new HashSet<>();
    Set<Integer> ringAtoms = new HashSet<>();
    Set<SearchResult.Bond> ringBonds = new HashSet<>();
    for (int[] ring : graph.computeRelevantCycles()) {
      for (int index = 0; index < ring.length; index++) {
        ringAtoms.add(ring[index]);
        ringBonds.add(new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]));
      }
      boolean complete = true;
      for (int index = 0; index < ring.length; index++)
        if (!selectedBonds.contains(
                new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]))) {
          complete = false;
          break;
        }
      if (complete)
        for (int index = 0; index < ring.length; index++) {
          coveredAtoms.add(ring[index]);
          coveredBonds.add(new SearchResult.Bond(ring[index], ring[(index + 1) % ring.length]));
        }
    }
    for (int atom : selectedAtoms)
      if ((graph.ring[atom] || ringAtoms.contains(atom)) && !coveredAtoms.contains(atom))
        return false;
    for (SearchResult.Bond bond : selectedBonds)
      if ((graph.bondInRing(bond.first(), bond.second()) || ringBonds.contains(bond))
          && !coveredBonds.contains(bond))
        return false;
    return true;
  }

  public static List<String> validateMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping, ChemOptions C) {
    List<String> errors = new ArrayList<>();
    if (mapping == null || mapping.isEmpty())
      return errors;
    if (C == null)
      C = new ChemOptions();
    Set<Integer> usedT = new HashSet<>();
    for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
      int qi = e.getKey(), tj = e.getValue();
      if (qi < 0 || qi >= g1.n) {
        errors.add("query index " + qi + " out of range [0," + g1.n + ")");
        continue;
      }
      if (tj < 0 || tj >= g2.n) {
        errors.add("target index " + tj + " out of range [0," + g2.n + ")");
        continue;
      }
      if (!usedT.add(tj))
        errors.add("target atom " + tj + " mapped twice (not injective)");
      if (!Substructure.MoleculeMatcher.atomsCompatFast(g1, qi, g2, tj, C))
        errors.add("atom mismatch: q" + qi + "(" + g1.atomicNum[qi] + ") → t" + tj + "("
            + g2.atomicNum[tj] + ")");
    }
    // Check bonds: every mapped query bond must exist in the target and be compatible.
    for (Map.Entry<Integer, Integer> e1 : mapping.entrySet()) {
      int qi = e1.getKey(), tj = e1.getValue();
      if (qi < 0 || qi >= g1.n || tj < 0 || tj >= g2.n)
        continue;
      for (int qk : g1.neighbors[qi]) {
        if (qk <= qi)
          continue;
        Integer tkObj = mapping.get(qk);
        if (tkObj == null)
          continue;
        int tk = tkObj;
        if (tk < 0 || tk >= g2.n)
          continue;
        if (g1.bondOrder(qi, qk) == 0)
          continue;
        if (g2.bondOrder(tj, tk) == 0) {
          errors.add("bond missing: q(" + qi + "-" + qk + ") vs t(" + tj + "-" + tk + ")");
          continue;
        }
        if (!MolGraph.ChemOps.bondsCompatible(g1, qi, qk, g2, tj, tk, C))
          errors.add("bond incompatible: q(" + qi + "-" + qk + ") vs t(" + tj + "-" + tk + ")");
      }
    }
    if (C.useChirality || C.useBondStereo) {
      int[] indexed = new int[g1.n];
      Arrays.fill(indexed, -1);
      for (Map.Entry<Integer, Integer> pair : mapping.entrySet())
        if (pair.getKey() >= 0 && pair.getKey() < g1.n && pair.getValue() >= 0
            && pair.getValue() < g2.n)
          indexed[pair.getKey()] = pair.getValue();
      int mismatch = Substructure.stereoMismatchAtom(g1, g2, indexed, C);
      if (mismatch >= 0)
        errors.add("atom mismatch: q" + mismatch + " (mapped stereochemistry)");
    }
    return errors;
  }

  public static Map<Integer, Integer> canonicalizeMapping(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mapping) {
    if (mapping.isEmpty())
      return mapping;
    int[][] queryGenerators = g1.getAutomorphismGenerators();
    int[][] targetGenerators = g2.getAutomorphismGenerators();
    if (g1.automorphismGeneratorsTruncated() || g2.automorphismGeneratorsTruncated())
      throw new IllegalStateException(
          "Cannot canonicalize mapping with incomplete automorphism generators");
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
        throw new IllegalStateException("Mapping canonicalization exceeded its time budget");
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
                "Mapping canonicalization exceeded its orbit state limit");
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

  public static double computeTautomerConfidence(
      IAtomContainer m1, IAtomContainer m2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty())
      return 1.0;
    MolGraph g1 = new MolGraph(m1);
    MolGraph g2 = new MolGraph(m2);
    return computeTautomerConfidence(g1, g2, mcs);
  }

  public static double computeTautomerConfidence(
      MolGraph g1, MolGraph g2, Map<Integer, Integer> mcs) {
    if (mcs == null || mcs.isEmpty())
      return 1.0;
    g1.ensureTautomerClasses();
    g2.ensureTautomerClasses();
    double logSum = 0.0;
    int count = 0;
    for (Map.Entry<Integer, Integer> e : mcs.entrySet()) {
      int qi = e.getKey(), ti = e.getValue();
      if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n)
        continue;
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
    query.ensureTautomerClasses();
    target.ensureTautomerClasses();
    double logSum = 0.0;
    int count = 0;
    for (Map.Entry<Integer, Integer> pair : match.atomMapping().entrySet()) {
      int queryAtom = pair.getKey(), targetAtom = pair.getValue();
      if (relaxed[queryAtom] && (query.tautomerClass[queryAtom] >= 0
          || target.tautomerClass[targetAtom] >= 0)) {
        double weight = Math.min(contextualTautomerWeight(query, queryAtom, options),
            contextualTautomerWeight(target, targetAtom, options));
        logSum += Math.log(Math.max(weight, 1e-10));
        count++;
      }
    }
    return count == 0 ? 1.0 : Math.exp(logSum / count);
  }

  private static double contextualTautomerWeight(MolGraph graph, int atom, ChemOptions options) {
    double weight = graph.tautomerWeight[atom];
    if (weight >= 1.0)
      return weight;
    double ketoShift, amideShift, generalShift;
    switch (options.solvent) {
      case DMSO: ketoShift = -0.15; amideShift = -0.08; generalShift = -0.05; break;
      case CHLOROFORM: ketoShift = -0.25; amideShift = -0.12; generalShift = -0.08; break;
      case METHANOL: ketoShift = -0.05; amideShift = -0.02; generalShift = -0.01; break;
      case ACETONITRILE: ketoShift = -0.10; amideShift = -0.05; generalShift = -0.03; break;
      case DIETHYL_ETHER: ketoShift = -0.20; amideShift = -0.10; generalShift = -0.06; break;
      default: ketoShift = amideShift = generalShift = 0.0;
    }
    if (options.solvent != ChemOptions.Solvent.AQUEOUS) {
      double shift = weight >= 0.96f
          ? (Math.abs(weight - MolGraph.TW_KETO_ENOL) < 0.005 || weight < 0.97f
                  ? ketoShift : amideShift)
          : weight >= 0.93f ? amideShift : generalShift;
      weight = Math.max(0.05, Math.min(0.99, weight + shift));
    }
    double dampening = 0.1 * Math.abs(options.pH - 7.4);
    if (dampening >= 0.001) {
      weight = weight > 0.5 ? Math.max(0.5, weight - dampening)
                            : Math.min(0.5, weight + dampening);
      weight = Math.max(0.1, Math.min(0.99, weight));
    }
    return weight;
  }

  public static IAtomContainer murckoScaffold(IAtomContainer mol) {
    MolGraph source=copyMolecule(mol);
    mol=source.mol;
    int n = mol.getAtomCount();
    if (n == 0)
      return mol;
    boolean[] inRing = new boolean[n];
    boolean hasRing = false;
    for (int i = 0; i < n; i++) {
      inRing[i] = mol.getAtom(i).isInRing();
      if (inRing[i])
        hasRing = true;
    }
    if (!hasRing)
      return mol;
    IdentityHashMap<IAtom, Integer> idxMap = new IdentityHashMap<>(n);
    for (int i = 0; i < n; i++) idxMap.put(mol.getAtom(i), i);
    int[][] adj = new int[n][];
    for (int i = 0; i < n; i++) {
      List<IAtom> nbAtoms = mol.getConnectedAtomsList(mol.getAtom(i));
      adj[i] = new int[nbAtoms.size()];
      for (int k = 0; k < adj[i].length; k++) adj[i][k] = idxMap.get(nbAtoms.get(k));
    }
    List<Integer> ringAtomList = new ArrayList<>();
    for (int i = 0; i < n; i++)
      if (inRing[i])
        ringAtomList.add(i);
    boolean[] keep = new boolean[n];
    for (int i = 0; i < n; i++)
      if (inRing[i])
        keep[i] = true;
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
          if (seen[v])
            continue;
          seen[v] = true;
          prev[v] = u;
          if (inRing[v] && v != src) {
            int cur = v;
            while (cur != src && cur != -1) {
              keep[cur] = true;
              cur = prev[cur];
            }
          }
          if (!inRing[v])
            queue.add(v);
        }
      }
    }
    Map<Integer,Integer> retained=new TreeMap<>();
    for(int atom=0;atom<n;atom++)if(keep[atom])retained.put(atom,atom);
    return toMolecule(source,retained,allBonds(source,retained).keySet());
  }

  private static int mobilisableProtons(MolGraph graph,int atom) {
    int element=graph.atomicNum[atom];
    return element==7||element==8||element==16?graph.hydrogenCount(atom):0;
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
    g1.ensureTautomerClasses();
    g2.ensureTautomerClasses();
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

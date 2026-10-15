/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.*;
import java.util.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/** Immutable test fixtures and text helpers available on Java 8. */
public final class TestSupport {
  private TestSupport() {}

  public static Map<Integer, Integer> mapping(int... pairs) {
    if ((pairs.length & 1) != 0)
      throw new IllegalArgumentException("Incomplete mapping pair");
    Map<Integer, Integer> values = new LinkedHashMap<>();
    for (int i = 0; i < pairs.length; i += 2) {
      if (values.put(pairs[i], pairs[i + 1]) != null)
        throw new IllegalArgumentException("Duplicate mapping key: " + pairs[i]);
    }
    return Collections.unmodifiableMap(values);
  }

  @SafeVarargs
  public static <T> List<T> list(T... elements) {
    List<T> values = new ArrayList<>(elements.length);
    for (T element : elements) values.add(Objects.requireNonNull(element));
    return Collections.unmodifiableList(values);
  }

  @SafeVarargs
  public static <T> Set<T> set(T... elements) {
    Set<T> values = new LinkedHashSet<>();
    for (T element : elements) {
      if (!values.add(Objects.requireNonNull(element)))
        throw new IllegalArgumentException("Duplicate set element: " + element);
    }
    return Collections.unmodifiableSet(values);
  }

  public static String repeat(String text, int count) {
    if (count < 0)
      throw new IllegalArgumentException("Negative repeat count");
    StringBuilder result = new StringBuilder();
    for (int i = 0; i < count; i++) result.append(text);
    return result.toString();
  }

  public static boolean isBlank(String text) {
    return text.codePoints().allMatch(Character::isWhitespace);
  }

  /** Use the public search facade while older molecular fixtures retain concise assertions. */
  public static Map<Integer, Integer> firstMapping(SearchResult result) {
    return result.found() ? result.matches().get(0).atomMapping() : Collections.emptyMap();
  }
  public static List<Map<Integer, Integer>> mappings(SearchResult result) {
    List<Map<Integer, Integer>> resultMaps = new ArrayList<>();
    for (SearchResult.Match match : result.matches()) resultMaps.add(match.atomMapping());
    return resultMaps;
  }
  private static int[] neighbors(MolGraph graph, int atom) {
    List<Integer> result = new ArrayList<>();
    for (int other = 0; other < graph.atomCount(); other++)
      if (graph.hasBond(atom, other))
        result.add(other);
    return result.stream().mapToInt(Integer::intValue).toArray();
  }
  public static MolGraph graph(IAtomContainer molecule) {
    return new MolGraph(molecule);
  }
  public static void clearCaches() {
    FingerprintEngine.clearFPCache();
  }
  private static MCS.Options configured(ChemOptions chemistry, MCS.Options source) {
    MCS.Options options = new MCS.Options();
    if (source != null) {
      options.objective = source.objective;
      options.connectivity = source.connectivity;
      options.atomWeights = source.atomWeights == null ? null : source.atomWeights.clone();
      options.timeoutMillis = source.timeoutMillis;
      options.nodeLimit = source.nodeLimit;
      options.cancellation = source.cancellation;
      options.minFragmentAtoms = source.minFragmentAtoms;
      options.maxFragments = source.maxFragments;
      options.excludedTargetAtoms = new HashSet<>(source.excludedTargetAtoms);
    }
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    return options;
  }
  public static Map<Integer, Integer> mcsMapping(
      MolGraph query, MolGraph target, ChemOptions chemistry, MCS.Options options) {
    return firstMapping(MCS.find(query).withOptions(configured(chemistry, options)).match(target));
  }
  public static Map<Integer, Integer> mcsMapping(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, MCS.Options options) {
    return mcsMapping(graph(query), graph(target), chemistry, options);
  }
  public static Map<Integer, Integer> mcsMapping(MolGraph query, MolGraph target) {
    return mcsMapping(query, target, null, null);
  }
  public static Map<Integer, Integer> mcsMapping(IAtomContainer query, IAtomContainer target) {
    return mcsMapping(query, target, null, null);
  }
  public static List<Map<Integer, Integer>> mcsMappings(
      MolGraph query, MolGraph target, ChemOptions chemistry, MCS.Options options, int limit) {
    return mappings(MCS.find(query)
            .withOptions(configured(chemistry, options))
            .matchAll(target, Math.max(1, limit)));
  }
  public static List<Map<Integer, Integer>> mcsMappings(IAtomContainer query, IAtomContainer target,
      ChemOptions chemistry, MCS.Options options, int limit) {
    return mcsMappings(graph(query), graph(target), chemistry, options, limit);
  }
  public static SearchResult substructureResult(
      MolGraph query, MolGraph target, ChemOptions chemistry, long timeout) {
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    options.timeoutMillis = timeout;
    return Substructure.find(query).withOptions(options).match(target);
  }
  public static SearchResult substructureResult(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, long timeout) {
    return substructureResult(graph(query), graph(target), chemistry, timeout);
  }
  public static boolean contains(
      MolGraph query, MolGraph target, ChemOptions chemistry, long timeout) {
    return substructureResult(query, target, chemistry, timeout).found();
  }
  public static boolean contains(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, long timeout) {
    return contains(graph(query), graph(target), chemistry, timeout);
  }
  public static SearchResult substructureResults(
      MolGraph query, MolGraph target, ChemOptions chemistry, int limit, long timeout) {
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry == null ? new ChemOptions() : chemistry;
    options.timeoutMillis = timeout;
    return Substructure.find(query).withOptions(options).matchAll(
        target, limit <= 0 ? 1000 : limit);
  }
  public static SearchResult substructureResults(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, int limit, long timeout) {
    return substructureResults(graph(query), graph(target), chemistry, limit, timeout);
  }
  public static List<Map<Integer, Integer>> substructureMappings(
      MolGraph query, MolGraph target, ChemOptions chemistry, int limit, long timeout) {
    return mappings(substructureResults(query, target, chemistry, limit, timeout));
  }
  public static List<Map<Integer, Integer>> substructureMappings(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, int limit, long timeout) {
    return substructureMappings(graph(query), graph(target), chemistry, limit, timeout);
  }
  public static Map<Integer, Integer> substructureMapping(
      MolGraph query, MolGraph target, ChemOptions chemistry, long timeout) {
    return firstMapping(substructureResult(query, target, chemistry, timeout));
  }
  public static Map<Integer, Integer> substructureMapping(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, long timeout) {
    return firstMapping(substructureResult(query, target, chemistry, timeout));
  }
  public static Map<Integer, Integer> substructureMapping(MolGraph query, MolGraph target) {
    return substructureMapping(query, target, null, -1);
  }
  public static Map<Integer, Integer> substructureMapping(
      IAtomContainer query, IAtomContainer target) {
    return substructureMapping(query, target, null, -1);
  }
  public static int mappedBondCount(MolGraph query, Map<Integer, Integer> mapping) {
    int count = 0;
    for (int atom : mapping.keySet())
      for (int neighbor : neighbors(query, atom))
        if (atom < neighbor && mapping.containsKey(neighbor))
          count++;
    return count;
  }
  public static int atomBound(MolGraph query, MolGraph target, ChemOptions chemistry) {
    double similarity = MCS.similarityUpperBound(query, target, chemistry);
    return (int) Math.round(
        similarity * (query.atomCount() + target.atomCount()) / (1 + similarity));
  }
  public static Map<Integer, Integer> largestConnected(
      MolGraph query, Map<Integer, Integer> mapping) {
    Map<Integer, Integer> best = Collections.emptyMap();
    Set<Integer> visited = new HashSet<>();
    for (int root : mapping.keySet()) {
      if (!visited.add(root))
        continue;
      Map<Integer, Integer> component = new TreeMap<>();
      ArrayDeque<Integer> queue = new ArrayDeque<>();
      queue.add(root);
      while (!queue.isEmpty()) {
        int atom = queue.remove();
        component.put(atom, mapping.get(atom));
        for (int neighbor : neighbors(query, atom))
          if (mapping.containsKey(neighbor) && visited.add(neighbor))
            queue.add(neighbor);
      }
      if (component.size() > best.size())
        best = component;
    }
    return best;
  }
  public static Map<Integer, Integer> fragmentMapping(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, MCS.Options options) {
    options = configured(chemistry, options);
    options.connectivity = MCS.Connectivity.Fragments;
    return mcsMapping(query, target, chemistry, options);
  }
  public static String fragmentSmiles(MolGraph query, Map<Integer, Integer> atoms) {
    List<SearchResult.Bond> selected = new ArrayList<>();
    for (int atom : atoms.keySet())
      for (int neighbor : neighbors(query, atom))
        if (atom < neighbor && atoms.containsKey(neighbor))
          selected.add(new SearchResult.Bond(atom, neighbor));
    return MCS.fragment(query, atoms, selected).toSmiles(query);
  }
  public static String fragmentSmiles(IAtomContainer query, Map<Integer, Integer> atoms) {
    return fragmentSmiles(graph(query), atoms);
  }
  public static String mcsSmiles(
      MolGraph query, MolGraph target, ChemOptions chemistry, MCS.Options options) {
    SearchResult result = MCS.find(query).withOptions(configured(chemistry, options)).match(target);
    return result.found() ? result.matches().get(0).toSmiles(query) : "";
  }
  public static String mcsSmiles(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, MCS.Options options) {
    return mcsSmiles(graph(query), graph(target), chemistry, options);
  }
  public static long[] fingerprintFromSmiles(String smiles, int radius, int size) throws Exception {
    return FingerprintEngine.ecfp(
        new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles), radius, size);
  }

  /** A standardised pair fixture; all searching is delegated to the public prepared queries. */
  private static IAtomContainer standardised(IAtomContainer molecule) {
    try {
      return Standardiser.standardise(molecule, Standardiser.TautomerMode.NONE);
    } catch (Exception failure) {
      return molecule;
    }
  }
  public static final class MoleculePair {
    private final IAtomContainer query, target;
    private final String smarts;
    private final ChemOptions chemistry;
    private long substructureTimeout = -1, mcsTimeout = -1;
    private Map<Integer, Integer> last = Collections.emptyMap();
    public MoleculePair(IAtomContainer query, IAtomContainer target, ChemOptions chemistry) {
      this(query, target, chemistry, true);
    }
    public MoleculePair(
        IAtomContainer query, IAtomContainer target, ChemOptions chemistry, boolean standardise) {
      Objects.requireNonNull(query, "query");
      Objects.requireNonNull(target, "target");
      this.query = standardise ? standardised(query) : query;
      this.target = standardise ? standardised(target) : target;
      this.chemistry = chemistry == null ? new ChemOptions() : chemistry;
      this.smarts = null;
    }
    public MoleculePair(String smarts, IAtomContainer target, ChemOptions chemistry) {
      this(smarts, target, chemistry, true);
    }
    public MoleculePair(
        String smarts, IAtomContainer target, ChemOptions chemistry, boolean standardise) {
      Objects.requireNonNull(smarts, "query");
      Objects.requireNonNull(target, "target");
      this.query = null;
      this.target = standardise ? standardised(target) : target;
      this.chemistry = chemistry == null ? new ChemOptions() : chemistry;
      this.smarts = smarts;
    }
    public MoleculePair(String query, String target, ChemOptions chemistry) throws Exception {
      this(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(query),
          new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(target), chemistry);
    }
    private Substructure prepared(long timeout) {
      Substructure.Options options = new Substructure.Options();
      options.chemistry = chemistry;
      options.timeoutMillis = timeout;
      return (smarts == null ? Substructure.find(query) : Substructure.smarts(smarts))
          .withOptions(options);
    }
    public SearchResult substructureResult(long timeout) {
      return prepared(timeout).match(target);
    }
    public SearchResult substructureResults(int limit, long timeout) {
      return prepared(timeout).matchAll(target, limit);
    }
    public boolean contains() {
      return contains(substructureTimeout);
    }
    public boolean contains(long timeout) {
      return prepared(timeout).match(target).found();
    }
    public List<Map<Integer, Integer>> substructureMappings(int limit, long timeout) {
      return mappings(prepared(timeout).matchAll(target, Math.max(1, limit)));
    }
    public Map<Integer, Integer> mcsMapping() {
      return mcsMapping(true, mcsTimeout);
    }
    public Map<Integer, Integer> mcsMapping(boolean connected, long timeout) {
      if (smarts != null)
        throw new IllegalStateException("MCS requires a molecule query");
      MCS.Options options = new MCS.Options();
      options.chemistry = chemistry;
      options.timeoutMillis = timeout;
      options.connectivity = connected ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
      last = firstMapping(MCS.find(query).withOptions(options).match(target));
      return last;
    }
    public List<Map<Integer, Integer>> mcsMappings(int limit) {
      return mcsMappings(limit, true, mcsTimeout);
    }
    public List<Map<Integer, Integer>> mcsMappings(int limit, boolean connected, long timeout) {
      MCS.Options options = new MCS.Options();
      options.chemistry = chemistry;
      options.timeoutMillis = timeout;
      options.connectivity = connected ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
      return mappings(MCS.find(query).withOptions(options).matchAll(target, Math.max(1, limit)));
    }
    public String mcsSmiles() {
      return TestSupport.fragmentSmiles(query, mcsMapping());
    }
    public String mcsSmiles(boolean connected, long timeout) {
      return TestSupport.fragmentSmiles(query, mcsMapping(connected, timeout));
    }
    public double similarityUpperBound() {
      return MCS.similarityUpperBound(query, target, chemistry);
    }
    public double getTautomerConfidence() {
      return MCS.computeTautomerConfidence(query, target, last);
    }
    public void setSubstructureTimeoutMs(long timeout) {
      substructureTimeout = timeout;
    }
    public void setMCSTimeoutMs(long timeout) {
      mcsTimeout = timeout;
    }
    public IAtomContainer query() {
      if (smarts != null)
        throw new UnsupportedOperationException("SMARTS has no molecule query");
      return query;
    }
    public IAtomContainer target() {
      return target;
    }
  }

  public static Map<Integer, Map<Integer, Integer>> pairMappings(List<IAtomContainer> molecules,
      ChemOptions chemistry, MCS.Options options, double threshold) {
    return pairMappings(molecules, chemistry, options, threshold, 0);
  }
  public static Map<Integer, Map<Integer, Integer>> pairMappings(List<IAtomContainer> molecules,
      ChemOptions chemistry, MCS.Options options, double threshold, int threads) {
    Map<Integer, Map<Integer, Integer>> result = new LinkedHashMap<>();
    int size = molecules.size();
    for (int first = 0; first < size; first++)
      for (int second = first + 1; second < size; second++)
        if (MCS.similarityUpperBound(molecules.get(first), molecules.get(second), chemistry)
            >= threshold)
          result.put(first * size + second,
              mcsMapping(molecules.get(first), molecules.get(second), chemistry, options));
    return result;
  }
  public static List<Map<Integer, Integer>> constrainedMappings(List<MolGraph> queries,
      List<MolGraph> targets, ChemOptions chemistry, MCS.Options source, int timeout) {
    List<Map<Integer, Integer>> result = new ArrayList<>();
    List<Set<Integer>> excluded = new ArrayList<>();
    for (MolGraph target : targets) excluded.add(new HashSet<>());
    for (MolGraph query : queries) {
      SearchResult.Match best = null;
      int bestTarget = -1;
      for (int index = 0; index < targets.size(); index++) {
        MCS.Options options = configured(chemistry, source);
        options.timeoutMillis = timeout;
        options.excludedTargetAtoms.addAll(excluded.get(index));
        SearchResult search = MCS.find(query).withOptions(options).match(targets.get(index));
        if (!search.found())
          continue;
        SearchResult.Match current = search.matches().get(0);
        if (best == null || Double.compare(current.score(), best.score()) > 0
            || (Double.compare(current.score(), best.score()) == 0
                && current.atomCount() > best.atomCount())) {
          best = current;
          bestTarget = index;
        }
      }
      result.add(best == null ? Collections.emptyMap() : best.atomMapping());
      if (bestTarget >= 0)
        excluded.get(bestTarget).addAll(best.atomMapping().values());
    }
    return result;
  }
  public static List<Map<String, IAtomContainer>> decompositions(
      IAtomContainer core, List<IAtomContainer> molecules, ChemOptions chemistry, long timeout) {
    List<Map<String, IAtomContainer>> results = new ArrayList<>();
    for (Substructure.Decomposition decomposition :
        Substructure.decomposeRGroups(core, molecules, chemistry, timeout)) {
      Map<String, IAtomContainer> components = new LinkedHashMap<>(decomposition.rgroups());
      if (decomposition.core() != null)
        components.put("core", decomposition.core());
      results.add(components);
    }
    return results;
  }
  public static Map<Integer, Integer> consensusMapping(
      List<IAtomContainer> molecules, ChemOptions chemistry, double threshold, long timeout) {
    MCS.ConsensusResult result = consensus(molecules, chemistry, threshold, timeout);
    return result.fragment().atomMapping();
  }
  public static IAtomContainer consensusMolecule(
      List<IAtomContainer> molecules, ChemOptions chemistry, double threshold, long timeout) {
    MCS.ConsensusResult result = consensus(molecules, chemistry, threshold, timeout);
    return result.fragment().toMolecule(graph(molecules.get(result.sourceIndex())));
  }
  private static MCS.ConsensusResult consensus(
      List<IAtomContainer> molecules, ChemOptions chemistry, double threshold, long timeout) {
    List<MolGraph> graphs = new ArrayList<>();
    for (IAtomContainer molecule : molecules) graphs.add(graph(molecule));
    MCS.Options options = new MCS.Options();
    options.chemistry = chemistry;
    options.timeoutMillis = timeout;
    return MCS.consensus(graphs, options, threshold);
  }
  public static Map<Integer, Integer> scaffoldMapping(
      IAtomContainer query, IAtomContainer target, ChemOptions chemistry, MCS.Options options) {
    return mcsMapping(MCS.murckoScaffold(query), MCS.murckoScaffold(target), chemistry, options);
  }
}

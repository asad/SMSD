/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.lang.ref.WeakReference;
import java.lang.reflect.Field;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class ObjectiveChemistryRegressionTest {
  private static MolGraph parse(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }

  private static MolGraph graph(int n, int[][] edges) {
    int[] elements = new int[n];
    Arrays.fill(elements, 6);
    List<List<Integer>> adjacency = new ArrayList<>();
    for (int i = 0; i < n; i++) adjacency.add(new ArrayList<>());
    for (int[] edge : edges) {
      adjacency.get(edge[0]).add(edge[1]);
      adjacency.get(edge[1]).add(edge[0]);
    }
    int[][] neighbors = adjacency.stream().map(row -> row.stream().mapToInt(Integer::intValue).toArray())
        .toArray(int[][]::new);
    return new MolGraph.Builder().atomCount(n).atomicNumbers(elements).neighbors(neighbors).build();
  }

  private static SearchEngine.MCSOptions weighted(double... weights) {
    var options = new SearchEngine.MCSOptions();
    options.atomWeights = weights;
    options.timeoutMs = 2_000;
    return options;
  }

  private static double score(Map<Integer, Integer> mapping, double[] weights) {
    return mapping.keySet().stream().mapToDouble(i -> weights[i]).sum();
  }

  @Test
  void tautomerMatchingPreservesElements() throws Exception {
    var chemistry = ChemOptions.tautomerProfile();
    MolGraph ketone = parse("CC(=O)C");
    for (String other : List.of("CC(=S)C", "CC(=O)N")) {
      MolGraph target = parse(other);
      Map<Integer, Integer> mapping = SearchEngine.findMCS(ketone, target, chemistry, new SearchEngine.MCSOptions());
      assertTrue(mapping.size() <= 3, "Only three atoms of each element can be shared with " + other);
      for (var pair : mapping.entrySet())
        assertEquals(ketone.atomicNum[pair.getKey()], target.atomicNum[pair.getValue()]);
      assertFalse(SearchEngine.isSubstructure(ketone, target, chemistry, 1_000));
    }
    assertEquals(4, SearchEngine.findMCS(ketone, parse("CC(O)=C"), chemistry,
        new SearchEngine.MCSOptions()).size());
  }

  @Test
  void tetrahedralMatchingUsesConfigurationRatherThanTraversalTag() throws Exception {
    MolGraph query = parse("N[C@@H](C)C(=O)O");
    MolGraph equivalent = parse("C[C@H](N)C(=O)O");
    MolGraph opposite = parse("C[C@@H](N)C(=O)O");
    var chemistry = new ChemOptions();
    chemistry.useChirality = true;
    Map<Integer, Integer> witness = Map.of(0, 2, 1, 1, 2, 0, 3, 3, 4, 4, 5, 5);
    assertTrue(SearchEngine.validateMapping(query, equivalent, witness, chemistry).isEmpty(),
        "The swapped N and methyl traversal and reversed @ tag describe the same configuration");
    assertFalse(SearchEngine.validateMapping(query, opposite, witness, chemistry).isEmpty());
    assertTrue(SearchEngine.isSubstructure(query, equivalent, chemistry, 1_000));
    assertFalse(SearchEngine.isSubstructure(query, opposite, chemistry, 1_000));
    assertEquals(6, SearchEngine.findMCS(query, equivalent, chemistry, new SearchEngine.MCSOptions()).size());
    assertEquals(3, SearchEngine.findMCS(query, opposite, chemistry, new SearchEngine.MCSOptions()).size());
  }

  @Test
  void stereoLigandOrderSurvivesBondStorageReordering() throws Exception {
    var molecule = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("N[C@@H](C)C(=O)O");
    var reordered = molecule.clone();
    var bonds = new org.openscience.cdk.interfaces.IBond[reordered.getBondCount()];
    for (int i = 0; i < bonds.length; i++) bonds[i] = reordered.getBond(bonds.length - 1 - i);
    reordered.setBonds(bonds);
    MolGraph query = new MolGraph(molecule), target = new MolGraph(reordered);
    assertEquals(Map.of(1, 'S'), CIPAssigner.assignRS(query));
    assertEquals(CIPAssigner.assignRS(query), CIPAssigner.assignRS(target));
    var chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertTrue(SearchEngine.validateMapping(query, target,
        Map.of(0,0,1,1,2,2,3,3,4,4,5,5), chemistry).isEmpty());
    assertTrue(SearchEngine.isSubstructure(query, target, chemistry, 1_000));
  }

  @Test
  void tetrahedralPermutationIsCheckedAfterNeighborsAreMapped() throws Exception {
    MolGraph molecule = parse("F[C@](Cl)(Br)I");
    var chemistry = new ChemOptions();
    chemistry.useChirality = true;
    chemistry.matchAtomType = false;
    assertFalse(SearchEngine.validateMapping(molecule, molecule,
        Map.of(0,2,1,1,2,0,3,3,4,4), chemistry).isEmpty(), "Swapping two ligands reverses winding");
    assertTrue(SearchEngine.validateMapping(molecule, molecule,
        Map.of(0,2,1,1,2,3,3,0,4,4), chemistry).isEmpty(), "A three-ligand rotation preserves winding");
    for (var engine : ChemOptions.MatcherEngine.values()) {
      chemistry.matcherEngine = engine;
      var mappings = SearchEngine.findAllSubstructures(molecule, molecule, chemistry, 100, 1_000);
      assertEquals(12, mappings.size(), "Half of the 4! ligand permutations preserve stereo for " + engine);
      for (var mapping : mappings)
        assertTrue(SearchEngine.validateMapping(molecule, molecule, mapping, chemistry).isEmpty());
    }
  }

  @Test
  void unspecifiedStereoRetainsItsExistingWildcardBehavior() throws Exception {
    MolGraph specified = parse("N[C@@H](C)C(=O)O"), unspecified = parse("NC(C)C(=O)O");
    var chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertTrue(SearchEngine.isSubstructure(specified, unspecified, chemistry, 1_000));
    assertTrue(SearchEngine.isSubstructure(unspecified, specified, chemistry, 1_000));
  }

  @Test
  void builderWindingAccountsForItsNeighborPermutation() {
    int[] elements = {6, 9, 17, 35, 53};
    MolGraph query = new MolGraph.Builder().atomCount(5).atomicNumbers(elements)
        .neighbors(new int[][] {{1,2,3,4},{0},{0},{0},{0}})
        .tetrahedralChirality(new int[] {1,0,0,0,0}).build();
    MolGraph target = new MolGraph.Builder().atomCount(5).atomicNumbers(elements)
        .neighbors(new int[][] {{2,1,3,4},{0},{0},{0},{0}})
        .tetrahedralChirality(new int[] {2,0,0,0,0}).build();
    var chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertEquals(CIPAssigner.assignRS(query), CIPAssigner.assignRS(target));
    assertTrue(SearchEngine.isSubstructure(query, target, chemistry, 1_000));
  }

  @Test
  void negativeBridgeDoesNotMakeAnIdentityMappingOptimal() throws Exception {
    MolGraph chain = parse("CCC");
    var options = weighted(1, -5, 1);
    Map<Integer, Integer> mapping = SearchEngine.findMCS(chain, chain, new ChemOptions(), options);
    assertEquals(1.0, score(mapping, options.atomWeights));
    assertEquals(1, mapping.size(), "A connected singleton scores better than the negative bridge");
    options.connectedOnly = false;
    assertEquals(2.0, score(SearchEngine.findMCS(chain, chain, new ChemOptions(), options), options.atomWeights));
    options.atomWeights = new double[] {-1, -2, -1};
    assertTrue(SearchEngine.findMCS(chain, chain, new ChemOptions(), options).isEmpty(),
        "The empty subgraph has score zero");
  }

  @Test
  void positiveWeightsCanPreferFewerMappedAtoms() {
    MolGraph query = graph(5, new int[][] {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});
    MolGraph target = graph(3, new int[][] {{0,1},{0,2},{1,2}});
    var options = weighted(1, 1, 1, 1, 100);
    options.connectedOnly = false;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(query, target, new ChemOptions(), options);
    assertEquals(102.0, score(mapping, options.atomWeights),
        "The heavy isolate and a clique edge embed together non-induced in a triangle");
    assertTrue(mapping.containsKey(4));
    assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void fragmentLimitsUseTheSelectedObjective() throws Exception {
    MolGraph components = parse("CCCC.CC");
    var options = weighted(1, 1, 1, 1, 5, 5);
    options.disconnectedMCS = true;
    options.maxFragments = 1;
    Map<Integer, Integer> mapping = SearchEngine.findMCS(components, components, new ChemOptions(), options);
    assertEquals(10.0, score(mapping, options.atomWeights));
    assertEquals(Set.of(4, 5), mapping.keySet());

    MolGraph bondComponents = graph(9, new int[][] {
        {0,1},{1,2},{2,3},{3,4},{5,6},{5,7},{5,8},{6,7},{6,8},{7,8}});
    var bondOptions = new SearchEngine.MCSOptions();
    bondOptions.disconnectedMCS = true;
    bondOptions.maxFragments = 1;
    bondOptions.maximizeBonds = true;
    assertEquals(6, SearchEngine.countMappedBonds(bondComponents,
        SearchEngine.findMCS(bondComponents, bondComponents, new ChemOptions(), bondOptions)));
  }

  @Test
  void smallWeightsKeepTheirPrecision() throws Exception {
    MolGraph components = parse("C.C");
    var options = weighted(0.0001, 0.0002);
    assertEquals(Set.of(1), SearchEngine.findMCS(components, components, new ChemOptions(), options).keySet());
  }

  @Test
  void cachedGraphDoesNotStronglyRetainTheWeakKey() throws Exception {
    SearchEngine.clearMolGraphCache();
    var molecule = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CCO");
    MolGraph graph = SearchEngine.toMolGraph(molecule);
    assertSame(graph, SearchEngine.toMolGraph(molecule));
    Field field = SearchEngine.class.getDeclaredField("molGraphCache");
    field.setAccessible(true);
    Map<?, ?> cache = (Map<?, ?>) field.get(null);
    Object cached = cache.get(molecule);
    assertInstanceOf(WeakReference.class, cached,
        "A strong cached MolGraph retains its IAtomContainer key and prevents weak-key eviction");
    assertSame(graph, ((WeakReference<?>) cached).get());
    SearchEngine.clearMolGraphCache();
  }

  @Test
  void canonicalizationFindsTheGlobalGeneratorOrbitMinimum() {
    MolGraph cycle = graph(4, new int[][] {{0,1},{1,2},{2,3},{3,0}});
    Map<Integer, Integer> mapping = Map.of(0, 2, 1, 3, 2, 0);
    Map<Integer, Integer> expected = orbitMinimum(cycle, cycle, mapping);
    assertEquals(expected, SearchEngine.canonicalizeMapping(cycle, cycle, mapping));
    for (int[] generator : cycle.getAutomorphismGenerators()) {
      Map<Integer, Integer> transformed = new TreeMap<>();
      mapping.forEach((q, t) -> transformed.put(q, generator[t]));
      assertEquals(expected, SearchEngine.canonicalizeMapping(cycle, cycle, transformed));
    }
  }

  @Test
  void objectiveModesAgreeWithAnIndependentThreeVertexOracle() {
    MolGraph[] graphs = new MolGraph[8];
    for (int mask = 0; mask < graphs.length; mask++) {
      List<int[]> edges = new ArrayList<>();
      int bit = 0;
      for (int a = 0; a < 3; a++)
        for (int b = a + 1; b < 3; b++)
          if ((mask & (1 << bit++)) != 0) edges.add(new int[] {a,b});
      graphs[mask] = graph(3, edges.toArray(int[][]::new));
    }
    for (int q = 0; q < graphs.length; q++) {
      for (int t = 0; t < graphs.length; t++) {
        for (boolean connected : new boolean[] {true, false}) {
          for (boolean induced : new boolean[] {true, false}) {
            for (double[] weights : new double[][] {{1,-5,2},{0.0001,0.0003,-0.0002},{1,2,3},null}) {
              var options = new SearchEngine.MCSOptions();
              options.connectedOnly = connected;
              options.induced = induced;
              options.atomWeights = weights;
              options.maximizeBonds = weights == null;
              options.timeoutMs = 1_000;
              Map<Integer, Integer> actual = SearchEngine.findMCS(graphs[q], graphs[t], new ChemOptions(), options);
              double expected = exhaustiveObjective(graphs[q], graphs[t], weights, connected, induced);
              double actualScore = weights == null ? SearchEngine.countMappedBonds(graphs[q], actual) : score(actual, weights);
              assertEquals(expected, actualScore, 1e-12,
                  "masks=" + q + "," + t + " connected=" + connected + " induced=" + induced + " weights=" + Arrays.toString(weights));
              assertTrue(SearchEngine.validateMapping(graphs[q], graphs[t], actual, new ChemOptions()).isEmpty());
            }
          }
        }
      }
    }
  }

  @Test
  void incompleteSymmetryGeneratorsAreReportedExplicitly() {
    MolGraph cycle = graph(4, new int[][] {{0,1},{1,2},{2,3},{3,0}});
    cycle.ensureCanonical();
    cycle.autGeneratorsTruncated = true;
    assertThrows(IllegalStateException.class,
        () -> SearchEngine.canonicalizeMapping(cycle, cycle, Map.of(0,0,1,1)));
  }

  @Test
  void skippedCanonicalSearchReportsIncompleteGenerators() {
    MolGraph large = graph(201, new int[][] {});
    assertTrue(large.automorphismGeneratorsTruncated(), "Large-graph refinement does not enumerate its automorphism group");
    assertThrows(IllegalStateException.class,
        () -> SearchEngine.canonicalizeMapping(large, large, Map.of(200,200)));
  }

  @Test
  void topologySymmetriesCannotExchangeDifferentChemicalProperties() throws Exception {
    MolGraph isotope = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6,6})
        .neighbors(new int[][] {{1},{0}}).massNumbers(new int[] {12,13}).build();
    MolGraph charge = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6,6})
        .neighbors(new int[][] {{},{}}).formalCharges(new int[] {1,-1}).build();
    MolGraph alternatingBonds = parse("C1=CC=C1");
    for (MolGraph molecule : List.of(isotope,charge,alternatingBonds)) {
      for (int[] generator : molecule.getAutomorphismGenerators()) {
        for (int atom = 0; atom < molecule.n; atom++) {
          assertEquals(molecule.formalCharge[atom], molecule.formalCharge[generator[atom]]);
          assertEquals(molecule.massNumber[atom], molecule.massNumber[generator[atom]]);
          assertEquals(molecule.hydrogenCount(atom), molecule.hydrogenCount(generator[atom]));
          for (int neighbor : molecule.neighbors[atom])
            assertEquals(molecule.bondOrder(atom,neighbor), molecule.bondOrder(generator[atom],generator[neighbor]));
        }
      }
      assertTrue(molecule.automorphismGeneratorsTruncated(),
          "Discarded topology generators cannot prove the full chemical automorphism group");
    }
    assertNotEquals(isotope.getOrbits()[0], isotope.getOrbits()[1]);
    assertNotEquals(charge.getOrbits()[0], charge.getOrbits()[1]);
  }

  @Test
  void annotatedCentersCannotBeMergedByOppositeStereoComponentExchange() throws Exception {
    MolGraph molecule = parse("F[C@](Cl)(Br)I.F[C@@](Cl)(Br)I");
    assertTrue(molecule.automorphismGeneratorsTruncated());
    for (int[] generator : molecule.getAutomorphismGenerators()) {
      assertEquals(1, generator[1]);
      assertEquals(6, generator[6]);
    }
    assertThrows(IllegalStateException.class,
        () -> SearchEngine.canonicalizeMapping(molecule,molecule,Map.of(6,6)));
  }

  @Test
  void safeSymmetryAwayFromAnAnnotatedCenterRemainsAvailable() throws Exception {
    MolGraph molecule = parse("CC(C)[C@H](F)Cl");
    assertFalse(molecule.automorphismGeneratorsTruncated());
    assertEquals(Map.of(0,0), SearchEngine.canonicalizeMapping(molecule,molecule,Map.of(2,2)));
  }

  @Test
  void nonfiniteWeightsAreRejected() throws Exception {
    MolGraph query = parse("CC");
    for (double invalid : new double[] {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY}) {
      var options = weighted(1, invalid);
      assertThrows(IllegalArgumentException.class,
          () -> SearchEngine.findMCS(query, query, new ChemOptions(), options));
    }
  }

  @Test
  void targetExclusionsApplyToIdentityAndDoNotLeakIntoLaterDomains() throws Exception {
    MolGraph molecule = parse("CCC");
    var chemistry = new ChemOptions();
    var options = weighted(100, 1, 1);
    options.excludedTargetAtoms = Set.of(0);
    Map<Integer, Integer> mapping = SearchEngine.findMCS(molecule, molecule, chemistry, options);
    assertFalse(mapping.values().contains(0));
    assertEquals(101.0, score(mapping, options.atomWeights),
        "Excluded target0 does not prevent mapping the high-weight query0 to target1");
    assertTrue(mapping.containsKey(0));
    options.atomWeights = null;
    options.induced = true;
    for (var engine : ChemOptions.MatcherEngine.values()) {
      chemistry.matcherEngine = engine;
      assertEquals(2, SearchEngine.findMCS(molecule, molecule, chemistry, options).size());
      assertEquals(3, SearchEngine.findMCS(molecule, molecule, chemistry, new SearchEngine.MCSOptions()).size());
    }
    assertNull(chemistry.mcsExcludedTargetAtoms, "Caller chemical options must remain unchanged");
    assertEquals(Set.of(0), options.excludedTargetAtoms);
  }

  @Test
  void excludedIndicesAreValidatedBeforeSearch() throws Exception {
    MolGraph molecule = parse("CC");
    var options = new SearchEngine.MCSOptions();
    for (int invalid : new int[] {-1, 2}) {
      options.excludedTargetAtoms = Set.of(invalid);
      assertThrows(IllegalArgumentException.class,
          () -> SearchEngine.findMCS(molecule, molecule, new ChemOptions(), options));
    }
  }

  @Test
  void constrainedBatchSupportsBuilderTargetsAndRetainsOriginalIndices() {
    MolGraph query = graph(2, new int[][] {{0,1}});
    MolGraph target = graph(4, new int[][] {{0,1},{1,2},{2,3}});
    var options = new SearchEngine.MCSOptions();
    var mappings = SearchEngine.batchMCSConstrained(List.of(query,query), List.of(target),
        new ChemOptions(), options, 1_000);
    assertEquals(2, mappings.size());
    Set<Integer> used = new HashSet<>();
    for (var mapping : mappings) {
      assertEquals(2, mapping.size());
      assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
      for (int atom : mapping.values()) assertTrue(used.add(atom), "Target atoms cannot be reused");
    }
    assertEquals(Set.of(0,1,2,3), used);
    assertNull(options.excludedTargetAtoms);
  }

  @Test
  void constrainedBatchRanksTargetsByObjective() {
    MolGraph query = new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6,6,7})
        .neighbors(new int[][] {{},{},{}}).build();
    MolGraph carbons = graph(2, new int[][] {{0,1}});
    MolGraph nitrogen = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {7})
        .neighbors(new int[][] {{}}).build();
    var options = weighted(1,1,100);
    options.connectedOnly = false;
    var mappings = SearchEngine.batchMCSConstrained(List.of(query), List.of(carbons,nitrogen),
        new ChemOptions(), options, 1_000);
    assertEquals(Map.of(2,0), mappings.get(0), "The single nitrogen scores100 versus two carbons scoring2");
  }

  @Test
  void constrainedBatchUsesItsPerPairTimeoutWithoutMutatingOptions() throws Exception {
    MolGraph query = parse("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl");
    MolGraph target = parse("CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6");
    query.ensureCanonical();
    target.ensureCanonical();
    var options = new SearchEngine.MCSOptions();
    options.timeoutMs = 10_000;
    long started = System.nanoTime();
    var mappings = SearchEngine.batchMCSConstrained(List.of(query), List.of(target),
        new ChemOptions(), options, 50);
    long elapsedMillis = (System.nanoTime() - started) / 1_000_000;
    assertTrue(elapsedMillis < 500, "The batch argument must override the10s option; elapsed=" + elapsedMillis);
    assertTrue(SearchEngine.validateMapping(query, target, mappings.get(0), new ChemOptions()).isEmpty());
    assertEquals(10_000, options.timeoutMs);
  }

  @Test
  void multipleMCSMappingsRetainTheIncumbentObjective() throws Exception {
    MolGraph query = parse("CCCC"), target = parse("CC");
    var options = weighted(100,1,1,1);
    var mappings = SearchEngine.findAllMCS(query, target, new ChemOptions(), options, 20);
    assertFalse(mappings.isEmpty());
    for (var mapping : mappings) {
      assertEquals(2, mapping.size());
      assertEquals(101.0, score(mapping, options.atomWeights));
      assertTrue(SearchEngine.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
    }
  }

  private static double exhaustiveObjective(MolGraph query, MolGraph target, double[] weights,
                                            boolean connected, boolean induced) {
    double best = 0.0;
    // Each base-4 digit is either an unmapped query atom or a target index.
    for (int code = 0; code < 64; code++) {
      int[] mapping = new int[3];
      int remaining = code;
      boolean valid = true;
      int size = 0;
      for (int q = 0; q < 3; q++) {
        mapping[q] = remaining % 4 - 1;
        remaining /= 4;
        if (mapping[q] >= 0) {
          size++;
          for (int prior = 0; prior < q; prior++) if (mapping[prior] == mapping[q]) valid = false;
        }
      }
      for (int a = 0; a < 3; a++) {
        if (mapping[a] < 0) continue;
        for (int b = a + 1; b < 3; b++) {
          if (mapping[b] < 0) continue;
          boolean queryEdge = query.hasBond(a,b), targetEdge = target.hasBond(mapping[a],mapping[b]);
          if (queryEdge && !targetEdge || induced && queryEdge != targetEdge) valid = false;
        }
      }
      if (!valid) continue;
      if (connected && size > 1) {
        Set<Integer> reached = new HashSet<>();
        for (int q = 0; q < 3; q++) if (mapping[q] >= 0) { reached.add(q); break; }
        boolean changed;
        do {
          changed = false;
          for (int a = 0; a < 3; a++) {
            if (mapping[a] < 0 || reached.contains(a)) continue;
            for (int b : reached.toArray(Integer[]::new))
              if (query.hasBond(a,b)) { reached.add(a); changed = true; break; }
          }
        } while (changed);
        if (reached.size() != size) continue;
      }
      double value = 0.0;
      for (int a = 0; a < 3; a++) {
        if (mapping[a] < 0) continue;
        if (weights != null) value += weights[a];
        else for (int b = a + 1; b < 3; b++) if (mapping[b] >= 0 && query.hasBond(a,b)) value++;
      }
      best = Math.max(best, value);
    }
    return best;
  }

  private static Map<Integer, Integer> orbitMinimum(MolGraph query, MolGraph target,
                                                    Map<Integer, Integer> mapping) {
    Comparator<Map<Integer, Integer>> order = Comparator.comparing(Object::toString);
    Set<Map<Integer, Integer>> visited = new HashSet<>();
    List<Map<Integer, Integer>> work = new ArrayList<>();
    work.add(new TreeMap<>(mapping));
    visited.add(work.get(0));
    for (int i = 0; i < work.size(); i++) {
      Map<Integer, Integer> current = work.get(i);
      for (int[] generator : query.getAutomorphismGenerators()) {
        Map<Integer, Integer> next = new TreeMap<>();
        current.forEach((q, t) -> next.put(generator[q], t));
        if (visited.add(next)) work.add(next);
      }
      for (int[] generator : target.getAutomorphismGenerators()) {
        Map<Integer, Integer> next = new TreeMap<>();
        current.forEach((q, t) -> next.put(q, generator[t]));
        if (visited.add(next)) work.add(next);
      }
    }
    return work.stream().min(order).orElseThrow();
  }
}

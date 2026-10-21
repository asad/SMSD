/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.core.MCS;

import com.bioinception.smsd.TestSupport;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static com.bioinception.smsd.TestSupport.list;
import static com.bioinception.smsd.TestSupport.mapping;
import static com.bioinception.smsd.TestSupport.set;

import java.lang.ref.WeakReference;
import java.lang.reflect.Field;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class ObjectiveChemistryRegressionTest {
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

  private static MCS.Options weighted(double... weights) {
    MCS.Options options = new MCS.Options();
    options.atomWeights = weights;
    options.timeoutMillis = 2_000;
    return options;
  }

  private static double score(Map<Integer, Integer> mapping, double[] weights) {
    return mapping.keySet().stream().mapToDouble(i -> weights[i]).sum();
  }

  @Test
  void tautomerMatchingPreservesElements() throws Exception {
    ChemOptions chemistry = ChemOptions.tautomerProfile();
    MolGraph ketone = rawGraph("CC(=O)C");
    for (String other : list("CC(=S)C", "CC(=O)N")) {
      MolGraph target = rawGraph(other);
      Map<Integer, Integer> mapping = TestSupport.mcsMapping(ketone, target, chemistry, new MCS.Options());
      assertTrue(mapping.size() <= 3, "Only three atoms of each element can be shared with " + other);
      for (Map.Entry<Integer, Integer> pair : mapping.entrySet())
        assertEquals(ketone.atomicNum[pair.getKey()], target.atomicNum[pair.getValue()]);
      assertFalse(TestSupport.contains(ketone, target, chemistry, 1_000));
    }
    assertEquals(4, TestSupport.mcsMapping(ketone, rawGraph("CC(O)=C"), chemistry,
        new MCS.Options()).size());
  }

  @Test
  void tetrahedralMatchingUsesConfigurationRatherThanTraversalTag() throws Exception {
    MolGraph query = rawGraph("N[C@@H](C)C(=O)O");
    MolGraph equivalent = rawGraph("C[C@H](N)C(=O)O");
    MolGraph opposite = rawGraph("C[C@@H](N)C(=O)O");
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    Map<Integer, Integer> witness = mapping(0, 2, 1, 1, 2, 0, 3, 3, 4, 4, 5, 5);
    assertTrue(MCS.validateMapping(query, equivalent, witness, chemistry).isEmpty(),
        "The swapped N and methyl traversal and reversed @ tag describe the same configuration");
    assertFalse(MCS.validateMapping(query, opposite, witness, chemistry).isEmpty());
    assertTrue(TestSupport.contains(query, equivalent, chemistry, 1_000));
    assertFalse(TestSupport.contains(query, opposite, chemistry, 1_000));
    assertEquals(6, TestSupport.mcsMapping(query, equivalent, chemistry, new MCS.Options()).size());
    Map<Integer, Integer> partial = TestSupport.mcsMapping(query, opposite, chemistry, new MCS.Options());
    assertEquals(5, partial.size());
    assertTrue(MCS.validateMapping(query, opposite, partial, chemistry).isEmpty());
  }

  @Test
  void stereoLigandOrderSurvivesBondStorageReordering() throws Exception {
    org.openscience.cdk.interfaces.IAtomContainer molecule = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("N[C@@H](C)C(=O)O");
    org.openscience.cdk.interfaces.IAtomContainer reordered = molecule.clone();
    org.openscience.cdk.interfaces.IBond[] bonds = new org.openscience.cdk.interfaces.IBond[reordered.getBondCount()];
    for (int i = 0; i < bonds.length; i++) bonds[i] = reordered.getBond(bonds.length - 1 - i);
    reordered.setBonds(bonds);
    MolGraph query = new MolGraph(molecule), target = new MolGraph(reordered);
    assertEquals(Collections.singletonMap(1, 'S'), CIPAssigner.assignRS(query));
    assertEquals(CIPAssigner.assignRS(query), CIPAssigner.assignRS(target));
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertTrue(MCS.validateMapping(query, target,
        mapping(0,0,1,1,2,2,3,3,4,4,5,5), chemistry).isEmpty());
    assertTrue(TestSupport.contains(query, target, chemistry, 1_000));
  }

  @Test
  void tetrahedralPermutationIsCheckedAfterNeighborsAreMapped() throws Exception {
    MolGraph molecule = rawGraph("F[C@](Cl)(Br)I");
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    chemistry.matchAtomType = false;
    assertFalse(MCS.validateMapping(molecule, molecule,
        mapping(0,2,1,1,2,0,3,3,4,4), chemistry).isEmpty(), "Swapping two ligands reverses winding");
    assertTrue(MCS.validateMapping(molecule, molecule,
        mapping(0,2,1,1,2,3,3,0,4,4), chemistry).isEmpty(), "A three-ligand rotation preserves winding");
    for (String engine : new String[] {"molecule"}) {
      List<Map<Integer, Integer>> mappings = TestSupport.substructureMappings(molecule, molecule, chemistry, 100, 1_000);
      assertEquals(12, mappings.size(), "Half of the 4! ligand permutations preserve stereo for " + engine);
      for (Map<Integer, Integer> mapping : mappings)
        assertTrue(MCS.validateMapping(molecule, molecule, mapping, chemistry).isEmpty());
    }
  }

  @Test
  void specifiedQueryCentresNeedSpecifiedTargetCentres() throws Exception {
    // An embedding follows the SMARTS rule; a common-subgraph search keeps the unspecified wildcard.
    MolGraph specified = rawGraph("N[C@@H](C)C(=O)O"), unspecified = rawGraph("NC(C)C(=O)O");
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertFalse(TestSupport.contains(specified, unspecified, chemistry, 1_000));
    assertTrue(TestSupport.contains(unspecified, specified, chemistry, 1_000));
    MCS.Options options = new MCS.Options();
    options.chemistry = chemistry;
    assertEquals(6, MCS.find(specified).withOptions(options).match(unspecified).matches().get(0).atomCount());
    chemistry.useChirality = false;
    assertTrue(TestSupport.contains(specified, unspecified, chemistry, 1_000));
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
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    assertEquals(CIPAssigner.assignRS(query), CIPAssigner.assignRS(target));
    assertTrue(TestSupport.contains(query, target, chemistry, 1_000));
  }

  @Test
  void negativeBridgeDoesNotMakeAnIdentityMappingOptimal() throws Exception {
    MolGraph chain = rawGraph("CCC");
    MCS.Options options = weighted(1, -5, 1);
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(chain, chain, new ChemOptions(), options);
    assertEquals(1.0, score(mapping, options.atomWeights));
    assertEquals(1, mapping.size(), "A connected singleton scores better than the negative bridge");
    options.connectivity = (false) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
    assertEquals(2.0, score(TestSupport.mcsMapping(chain, chain, new ChemOptions(), options), options.atomWeights));
    options.atomWeights = new double[] {-1, -2, -1};
    assertTrue(TestSupport.mcsMapping(chain, chain, new ChemOptions(), options).isEmpty(),
        "The empty subgraph has score zero");
  }

  @Test
  void positiveWeightsCanPreferFewerMappedAtoms() {
    MolGraph query = graph(5, new int[][] {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});
    MolGraph target = graph(3, new int[][] {{0,1},{0,2},{1,2}});
    MCS.Options options = weighted(1, 1, 1, 1, 100);
    options.connectivity = (false) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, new ChemOptions(), options);
    assertEquals(100.0, score(mapping, options.atomWeights),
        "The heavy isolate cannot share an induced image with a triangle edge");
    assertTrue(mapping.containsKey(4));
    assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
  }

  @Test
  void fragmentLimitsUseTheSelectedObjective() throws Exception {
    MolGraph components = rawGraph("CCCC.CC");
    MCS.Options options = weighted(1, 1, 1, 1, 5, 5);
    options.connectivity = (true) ? MCS.Connectivity.Fragments : MCS.Connectivity.Connected;
    options.maxFragments = 1;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(components, components, new ChemOptions(), options);
    assertEquals(10.0, score(mapping, options.atomWeights));
    assertEquals(set(4, 5), mapping.keySet());

    MolGraph bondComponents = graph(9, new int[][] {
        {0,1},{1,2},{2,3},{3,4},{5,6},{5,7},{5,8},{6,7},{6,8},{7,8}});
    MCS.Options bondOptions = new MCS.Options();
    bondOptions.connectivity = (true) ? MCS.Connectivity.Fragments : MCS.Connectivity.Connected;
    bondOptions.maxFragments = 1;
    bondOptions.objective = (true) ? MCS.Objective.MCES : MCS.Objective.MCIS;
    assertEquals(6, TestSupport.mappedBondCount(bondComponents,
        TestSupport.mcsMapping(bondComponents, bondComponents, new ChemOptions(), bondOptions)));
  }

  @Test
  void singleFragmentEdgeSearchKeepsTheDenserComponent() {
    MolGraph components = graph(9, new int[][] {
        {0,1},{1,2},{2,3},{3,4},{5,6},{5,7},{5,8},{6,7},{6,8},{7,8}});
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    options.maxFragments = 1;
    options.nodeLimit = 10_000;
    options.timeoutMillis = 60_000;
    SearchResult result = MCS.find(components).withOptions(options).match(components);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    SearchResult.Match match = result.matches().get(0);
    assertEquals(6, match.bondCount());
    assertEquals(set(5, 6, 7, 8), match.atomMapping().keySet());
    assertEquals(1, match.fragments().size());
    assertTrue(MCS.validateMapping(components, components, match, options).isEmpty());
  }

  @Test
  void smallWeightsKeepTheirPrecision() throws Exception {
    MolGraph components = rawGraph("C.C");
    MCS.Options options = weighted(0.0001, 0.0002);
    assertEquals(set(1), TestSupport.mcsMapping(components, components, new ChemOptions(), options).keySet());
  }

  @Test
  void preparedGraphSupportsConcurrentReuseWithoutAGlobalCache() throws Exception {
    MolGraph query=rawGraph("CCO");MCS prepared=MCS.find(query);
    java.util.concurrent.ExecutorService executor=java.util.concurrent.Executors.newFixedThreadPool(2);
    try {
      java.util.concurrent.Future<SearchResult> first=executor.submit(()->prepared.match(rawGraph("CCCO")));
      java.util.concurrent.Future<SearchResult> second=executor.submit(()->prepared.match(rawGraph("OCCC")));
      assertEquals(3,first.get().matches().get(0).atomCount());
      assertEquals(3,second.get().matches().get(0).atomCount());
      assertEquals(3,query.n);
      assertTrue(prepared.match(query).optimal());
    } finally {executor.shutdownNow();}
  }

  @Test
  void canonicalizationFindsTheGlobalGeneratorOrbitMinimum() {
    MolGraph cycle = graph(4, new int[][] {{0,1},{1,2},{2,3},{3,0}});
    Map<Integer, Integer> mapping = mapping(0, 2, 1, 3, 2, 0);
    Map<Integer, Integer> expected = orbitMinimum(cycle, cycle, mapping);
    assertEquals(expected, MCS.canonicalizeMapping(cycle, cycle, mapping));
    for (int[] generator : cycle.getAutomorphismGenerators()) {
      Map<Integer, Integer> transformed = new TreeMap<>();
      mapping.forEach((q, t) -> transformed.put(q, generator[t]));
      assertEquals(expected, MCS.canonicalizeMapping(cycle, cycle, transformed));
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
      graphs[mask] = graph(3, edges.toArray(new int[edges.size()][]));
    }
    for (int q = 0; q < graphs.length; q++) {
      for (int t = 0; t < graphs.length; t++) {
        for (boolean connected : new boolean[] {true, false}) {
          for (boolean induced : new boolean[] {true}) {
            for (double[] weights : new double[][] {{1,-5,2},{0.0001,0.0003,-0.0002},{1,2,3}}) {
              MCS.Options options = new MCS.Options();
              options.connectivity = (connected) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
              options.atomWeights = weights;
              options.objective = (weights == null) ? MCS.Objective.MCES : MCS.Objective.MCIS;
              options.timeoutMillis = 1_000;
              Map<Integer, Integer> actual = TestSupport.mcsMapping(graphs[q], graphs[t], new ChemOptions(), options);
              double expected = exhaustiveObjective(graphs[q], graphs[t], weights, connected, induced);
              double actualScore = weights == null ? TestSupport.mappedBondCount(graphs[q], actual) : score(actual, weights);
              assertEquals(expected, actualScore, 1e-12,
                  "masks=" + q + "," + t + " connected=" + connected + " induced=" + induced + " weights=" + Arrays.toString(weights));
              assertTrue(MCS.validateMapping(graphs[q], graphs[t], actual, new ChemOptions()).isEmpty());
            }
          }
        }
      }
    }
  }

  @Test
  void cancellationAndLargeOppositeWeightsRetainValidDeterministicScores() throws Exception {
    MolGraph chain=rawGraph("CCC(C)CCC"),ring=rawGraph("C1CCCCCC1");
    java.util.concurrent.atomic.AtomicInteger checks=new java.util.concurrent.atomic.AtomicInteger();
    MCS.Options limited=new MCS.Options();limited.cancellation=()->checks.incrementAndGet()>50;
    SearchResult stopped=MCS.find(chain).withOptions(limited).match(ring);
    assertEquals(SearchResult.Status.Cancelled,stopped.status());assertFalse(stopped.optimal());assertTrue(stopped.found());
    assertTrue(MCS.validateMapping(chain,ring,stopped.matches().get(0),limited).isEmpty());
    for(int queryMask=0;queryMask<8;queryMask++)for(int targetMask=0;targetMask<8;targetMask++) {
      MolGraph query=threeVertexGraph(queryMask),target=threeVertexGraph(targetMask);
      for(double[] weights:new double[][]{{1e16,-1e16,3},{-1e16,1e16,3},{1e-300,3e-300,-2e-300}}) {
        MCS.Options options=weighted(weights);options.connectivity=MCS.Connectivity.Fragments;
        SearchResult result=MCS.find(query).withOptions(options).match(target);
        assertTrue(result.optimal());
        assertEquals(exhaustiveObjective(query,target,weights,false,true),result.score(),0.0);
        if(result.found())assertEquals(score(result.matches().get(0).atomMapping(),weights),result.score(),0.0);
      }
    }
  }

  private static MolGraph threeVertexGraph(int mask) {
    List<int[]> edges=new ArrayList<>();int bit=0;
    for(int first=0;first<3;first++)for(int second=first+1;second<3;second++)
      if((mask&(1<<bit++))!=0)edges.add(new int[]{first,second});
    return graph(3,edges.toArray(new int[edges.size()][]));
  }

  @Test
  void incompleteSymmetryGeneratorsAreReportedExplicitly() {
    MolGraph cycle = graph(4, new int[][] {{0,1},{1,2},{2,3},{3,0}});
    cycle.ensureCanonical();
    cycle.autGeneratorsTruncated = true;
    assertThrows(IllegalStateException.class,
        () -> MCS.canonicalizeMapping(cycle, cycle, mapping(0,0,1,1)));
  }

  @Test
  void skippedCanonicalSearchReportsIncompleteGenerators() {
    MolGraph large = graph(201, new int[][] {});
    assertTrue(large.automorphismGeneratorsTruncated(), "Large-graph refinement does not enumerate its automorphism group");
    assertThrows(IllegalStateException.class,
        () -> MCS.canonicalizeMapping(large, large, mapping(200,200)));
  }

  @Test
  void molecularSymmetriesPreserveChemicalProperties() throws Exception {
    MolGraph isotope = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6,6})
        .neighbors(new int[][] {{1},{0}}).massNumbers(new int[] {12,13}).build();
    MolGraph charge = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6,6})
        .neighbors(new int[][] {{},{}}).formalCharges(new int[] {1,-1}).build();
    MolGraph alternatingBonds = rawGraph("C1=CC=C1");
    for (MolGraph molecule : list(isotope,charge,alternatingBonds)) {
      for (int[] generator : molecule.getAutomorphismGenerators()) {
        for (int atom = 0; atom < molecule.n; atom++) {
          assertEquals(molecule.formalCharge[atom], molecule.formalCharge[generator[atom]]);
          assertEquals(molecule.massNumber[atom], molecule.massNumber[generator[atom]]);
          assertEquals(molecule.hydrogenCount(atom), molecule.hydrogenCount(generator[atom]));
          for (int neighbor : molecule.neighbors[atom])
            assertEquals(molecule.bondOrder(atom,neighbor), molecule.bondOrder(generator[atom],generator[neighbor]));
        }
      }
      assertFalse(molecule.automorphismGeneratorsTruncated());
    }
    assertNotEquals(isotope.getOrbits()[0], isotope.getOrbits()[1]);
    assertNotEquals(charge.getOrbits()[0], charge.getOrbits()[1]);
    assertNotEquals(MCS.canonicalizeMapping(isotope,isotope,mapping(0,0)),
        MCS.canonicalizeMapping(isotope,isotope,mapping(1,1)));
    assertEquals(MCS.canonicalizeMapping(alternatingBonds,alternatingBonds,mapping(0,0)),
        MCS.canonicalizeMapping(alternatingBonds,alternatingBonds,mapping(2,2)));
  }

  @Test
  void annotatedCentersCannotBeMergedByOppositeStereoComponentExchange() throws Exception {
    MolGraph molecule = rawGraph("F[C@](Cl)(Br)I.F[C@@](Cl)(Br)I");
    assertFalse(molecule.automorphismGeneratorsTruncated());
    for (int[] generator : molecule.getAutomorphismGenerators()) {
      assertEquals(1, generator[1]);
      assertEquals(6, generator[6]);
    }
    assertEquals(mapping(6,6), MCS.canonicalizeMapping(molecule,molecule,mapping(6,6)));
  }

  @Test
  void safeSymmetryAwayFromAnAnnotatedCenterRemainsAvailable() throws Exception {
    MolGraph molecule = rawGraph("CC(C)[C@H](F)Cl");
    assertFalse(molecule.automorphismGeneratorsTruncated());
    assertEquals(mapping(0,0), MCS.canonicalizeMapping(molecule,molecule,mapping(2,2)));
  }

  @Test
  void nonfiniteWeightsAreRejected() throws Exception {
    MolGraph query = rawGraph("CC");
    for (double invalid : new double[] {Double.NaN, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY}) {
      MCS.Options options = weighted(1, invalid);
      assertThrows(IllegalArgumentException.class,
          () -> TestSupport.mcsMapping(query, query, new ChemOptions(), options));
    }
  }

  @Test
  void targetExclusionsApplyToIdentityAndDoNotLeakIntoLaterDomains() throws Exception {
    MolGraph molecule = rawGraph("CCC");
    ChemOptions chemistry = new ChemOptions();
    MCS.Options options = weighted(100, 1, 1);
    options.excludedTargetAtoms = set(0);
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(molecule, molecule, chemistry, options);
    assertFalse(mapping.values().contains(0));
    assertEquals(101.0, score(mapping, options.atomWeights),
        "Excluded target0 does not prevent mapping the high-weight query0 to target1");
    assertTrue(mapping.containsKey(0));
    options.atomWeights = null;
    for (String engine : new String[] {"molecule"}) {
      assertEquals(2, TestSupport.mcsMapping(molecule, molecule, chemistry, options).size());
      assertEquals(3, TestSupport.mcsMapping(molecule, molecule, chemistry, new MCS.Options()).size());
    }
    Substructure.Options substructure=new Substructure.Options();substructure.chemistry=chemistry;
    assertEquals(3,Substructure.find(molecule).withOptions(substructure).match(molecule).matches().get(0).atomCount());
    assertEquals(set(0), options.excludedTargetAtoms);
  }

  @Test
  void excludedIndicesAreValidatedBeforeSearch() throws Exception {
    MolGraph molecule = rawGraph("CC");
    MCS.Options options = new MCS.Options();
    for (int invalid : new int[] {-1, 2}) {
      options.excludedTargetAtoms = set(invalid);
      assertThrows(IllegalArgumentException.class,
          () -> TestSupport.mcsMapping(molecule, molecule, new ChemOptions(), options));
    }
  }

  @Test
  void nullTargetExclusionsAreRejectedBeforeSearch() throws Exception {
    MolGraph molecule = rawGraph("CC");
    MCS.Options options = new MCS.Options();
    options.excludedTargetAtoms = new HashSet<>();
    options.excludedTargetAtoms.add(null);
    assertThrows(NullPointerException.class,
        () -> TestSupport.mcsMapping(molecule, molecule, new ChemOptions(), options));
    assertEquals(Collections.singleton(null), options.excludedTargetAtoms);
    assertEquals(2, TestSupport.mcsMapping(molecule, molecule, new ChemOptions(),
        new MCS.Options()).size());
  }

  @Test
  void constrainedBatchSupportsBuilderTargetsAndRetainsOriginalIndices() {
    MolGraph query = graph(2, new int[][] {{0,1}});
    MolGraph target = graph(4, new int[][] {{0,1},{1,2},{2,3}});
    MCS.Options options = new MCS.Options();
    List<Map<Integer, Integer>> mappings = TestSupport.constrainedMappings(list(query,query), list(target),
        new ChemOptions(), options, 1_000);
    assertEquals(2, mappings.size());
    Set<Integer> used = new HashSet<>();
    for (Map<Integer, Integer> mapping : mappings) {
      assertEquals(2, mapping.size());
      assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
      for (int atom : mapping.values()) assertTrue(used.add(atom), "Target atoms cannot be reused");
    }
    assertEquals(set(0,1,2,3), used);
    assertTrue(options.excludedTargetAtoms.isEmpty());
  }

  @Test
  void constrainedBatchRanksTargetsByObjective() {
    MolGraph query = new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6,6,7})
        .neighbors(new int[][] {{},{},{}}).build();
    MolGraph carbons = graph(2, new int[][] {{0,1}});
    MolGraph nitrogen = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {7})
        .neighbors(new int[][] {{}}).build();
    MCS.Options options = weighted(1,1,100);
    options.connectivity = (false) ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
    List<Map<Integer, Integer>> mappings = TestSupport.constrainedMappings(list(query), list(carbons,nitrogen),
        new ChemOptions(), options, 1_000);
    assertEquals(mapping(2,0), mappings.get(0), "The single nitrogen scores100 versus two carbons scoring2");
  }

  @Test
  void constrainedBatchUsesItsPerPairTimeoutWithoutMutatingOptions() throws Exception {
    MolGraph query = rawGraph("c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl");
    MolGraph target = rawGraph("CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6");
    query.ensureCanonical();
    target.ensureCanonical();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 10_000;
    long started = System.nanoTime();
    List<Map<Integer, Integer>> mappings = TestSupport.constrainedMappings(list(query), list(target),
        new ChemOptions(), options, 50);
    long elapsedMillis = (System.nanoTime() - started) / 1_000_000;
    assertTrue(elapsedMillis < 500, "The batch argument must override the10s option; elapsed=" + elapsedMillis);
    assertTrue(MCS.validateMapping(query, target, mappings.get(0), new ChemOptions()).isEmpty());
    assertEquals(10_000, options.timeoutMillis);
  }

  @Test
  void multipleMCSMappingsRetainTheIncumbentObjective() throws Exception {
    MolGraph query = rawGraph("CCCC"), target = rawGraph("CC");
    MCS.Options options = weighted(100,1,1,1);
    List<Map<Integer, Integer>> mappings = TestSupport.mcsMappings(query, target, new ChemOptions(), options, 20);
    assertFalse(mappings.isEmpty());
    for (Map<Integer, Integer> mapping : mappings) {
      assertEquals(2, mapping.size());
      assertEquals(101.0, score(mapping, options.atomWeights));
      assertTrue(MCS.validateMapping(query, target, mapping, new ChemOptions()).isEmpty());
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
            for (int b : reached.toArray(new Integer[reached.size()]))
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
    return work.stream().min(order).orElseThrow(NoSuchElementException::new);
  }
}

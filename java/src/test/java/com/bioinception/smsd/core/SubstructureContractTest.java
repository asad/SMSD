/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smarts.SmartsPattern;
import org.openscience.cdk.smiles.SmilesParser;

class SubstructureContractTest {
  private static IAtomContainer molecule(String smiles) throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles);
  }

  @ParameterizedTest
  @CsvSource(value = {
      "[CX3](=O)[NX3]|CC(=O)Nc1ccc(O)cc1|true",
      "[#6;R]|CCO|false",
      "[#6;R]|C1CCCCC1|true",
      "[$(C(=O)N)]|CC(=O)N|true",
      "[C;!$(C=O)]|C=O|false",
      "[$([C;$(C(=O)N)])]|CC(=O)N|true",
      "[N+;H3]|C[NH3+]|true",
      "[N+;H3]|CN|false",
      "[13C]|[13CH4]|true",
      "[13C]|C|false",
      "C-,=C|C=C|true",
      "C-,=C|CCC|true",
      "C-,=C|C#C|false",
      "N[C@@H](C)C(=O)O|N[C@@H](C)C(=O)O|true",
      "N[C@@H](C)C(=O)O|N[C@H](C)C(=O)O|false",
      "F/C=C/Cl|F/C=C/Cl|true",
      "F/C=C/Cl|F/C=C\\Cl|false",
      "[$(N[C@@H](C)C(=O)O)]|N[C@@H](C)C(=O)O|true",
      "[$(N[C@@H](C)C(=O)O)]|N[C@H](C)C(=O)O|false",
      "(C).(O)|CO|false",
      "(C).(O)|C.O|true"
  }, delimiter = '|')
  void preservesCompleteSmartsPredicates(String query, String target, boolean expected) throws Exception {
    IAtomContainer parsed = molecule(target);
    assertEquals(expected, SmartsPattern.create(query).matches(parsed.clone()));
    SearchResult result = Substructure.smarts(query).withTimeout(1000).match(parsed);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertEquals(expected, result.found());
    assertTrue(result.complete());
    assertTrue(result.optimal());
  }

  @Test
  void firstMatchDoesNotEnumerateDisconnectedEmbeddings() throws Exception {
    Substructure.smarts("C").match(molecule("CC"));
    String chain = String.join("", Collections.nCopies(40, "C"));
    Substructure query = Substructure.smarts("C.C.C.C").withTimeout(100);
    IAtomContainer target = molecule(chain);
    SearchResult result = assertTimeout(Duration.ofSeconds(1), () -> query.match(target));
    assertTrue(result.found());
    assertEquals(1, result.matches().size());
    assertEquals(4, result.matches().get(0).atomCount());
    assertEquals(4, result.matches().get(0).fragments().size());
  }

  @Test
  void distinguishesExhaustionFromTheResultLimit() throws Exception {
    Substructure query = Substructure.smarts("O");
    SearchResult one = query.matchAll(molecule("CO"), 1);
    assertEquals(SearchResult.Status.Completed, one.status());
    assertTrue(one.complete());
    SearchResult limited = query.matchAll(molecule("COCO"), 1);
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertEquals(1, limited.matches().size());
    assertTrue(limited.optimal());
    assertFalse(limited.complete());
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void sharesTimeoutAndNodeLimitsAcrossBothQueryTypes(boolean smarts) throws Exception {
    MolGraph target = new MolGraph(molecule("CCCCCCCCCCCC"));
    Substructure query = smarts ? Substructure.smarts("CCCC") : Substructure.find(new MolGraph(molecule("CCCC")));
    SearchResult expired = query.withTimeout(0).match(target);
    assertEquals(SearchResult.Status.TimedOut, expired.status());
    assertFalse(expired.complete());
    assertFalse(expired.optimal());
    Substructure.Options options = new Substructure.Options();
    options.nodeLimit = 1;
    SearchResult limited = query.withOptions(options).match(target);
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertTrue(limited.nodesVisited() <= 1);
    assertFalse(limited.complete());
  }

  @Test
  void propagatesCancellationInsideRecursiveQueries() throws Exception {
    AtomicInteger checks = new AtomicInteger();
    Substructure.Options options = new Substructure.Options();
    options.cancellation = () -> checks.incrementAndGet() > 200;
    options.timeoutMillis = 1000;
    SearchResult result = Substructure.smarts("[$(C1CCCCCCCCCCCCCCCCCCC1)]").withOptions(options)
        .match(molecule(String.join("", Collections.nCopies(100, "C"))));
    assertEquals(SearchResult.Status.Cancelled, result.status());
    assertTrue(result.nodesVisited() > 1);
    assertFalse(result.found());
    assertFalse(result.complete());
  }

  @Test
  void reportsSelectedBondsAndTheirFragments() throws Exception {
    SearchResult result = Substructure.find(molecule("CC.O")).match(molecule("CCO"));
    assertTrue(result.found());
    SearchResult.Match match = result.matches().get(0);
    assertEquals(3, match.atomCount());
    assertEquals(1, match.bondCount());
    assertEquals(2, match.fragments().size());
    assertEquals(new HashSet<>(Arrays.asList(0, 1, 2)), match.atomMapping().keySet());
    assertThrows(UnsupportedOperationException.class, () -> match.atomMapping().put(0, 0));
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void optionallyRequiresAnInducedEmbedding(boolean smarts) throws Exception {
    Substructure query = smarts ? Substructure.smarts("CCC") : Substructure.find(molecule("CCC"));
    IAtomContainer triangle = molecule("C1CC1");
    assertTrue(query.match(triangle).found());
    Substructure.Options options = new Substructure.Options();
    options.induced = true;
    assertFalse(query.withOptions(options).match(triangle).found());
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void moleculeEmbeddingsAgreeWithAllThreeAtomPermutations(boolean induced) {
    for (int queryMask = 0; queryMask < 8; queryMask++) {
      MolGraph query = threeAtomGraph(queryMask);
      Substructure.Options options = new Substructure.Options();
      options.induced = induced;
      Substructure prepared = Substructure.find(query).withOptions(options);
      for (int targetMask = 0; targetMask < 8; targetMask++) {
        MolGraph target = threeAtomGraph(targetMask);
        Set<List<Integer>> expected = new HashSet<>();
        for (int first = 0; first < 3; first++)
          for (int second = 0; second < 3; second++)
            for (int third = 0; third < 3; third++) {
              if (first == second || first == third || second == third) continue;
              int[] mapping = {first, second, third};
              boolean valid = true;
              for (int atom = 0; atom < 3; atom++)
                for (int neighbor = atom + 1; neighbor < 3; neighbor++) {
                  boolean queryBond = query.hasBond(atom, neighbor);
                  boolean targetBond = target.hasBond(mapping[atom], mapping[neighbor]);
                  if (queryBond && !targetBond || induced && queryBond != targetBond) valid = false;
                }
              if (valid) expected.add(Arrays.asList(first, second, third));
            }
        SearchResult result = prepared.matchAll(target);
        Set<List<Integer>> actual = new HashSet<>();
        for (SearchResult.Match match : result.matches())
          actual.add(Arrays.asList(match.atomMapping().get(0), match.atomMapping().get(1),
              match.atomMapping().get(2)));
        assertEquals(SearchResult.Status.Completed, result.status());
        assertEquals(expected, actual, "query=" + queryMask + ", target=" + targetMask);
      }
    }
  }

  private static MolGraph threeAtomGraph(int mask) {
    List<List<Integer>> rows = Arrays.asList(new ArrayList<>(), new ArrayList<>(), new ArrayList<>());
    int bit = 0;
    for (int atom = 0; atom < 3; atom++)
      for (int neighbor = atom + 1; neighbor < 3; neighbor++)
        if ((mask & (1 << bit++)) != 0) {
          rows.get(atom).add(neighbor);
          rows.get(neighbor).add(atom);
        }
    int[][] neighbors = rows.stream().map(row -> row.stream().mapToInt(Integer::intValue).toArray())
        .toArray(int[][]::new);
    return new MolGraph.Builder().atomCount(3).atomicNumbers(new int[] {6, 6, 6})
        .neighbors(neighbors).build();
  }

  @Test
  void queryOptionsAreAnIndependentSnapshot() throws Exception {
    Substructure.Options options = new Substructure.Options();
    IAtomContainer target = molecule("[NH4+]");
    Substructure query = Substructure.find(molecule("N")).withOptions(options);
    options.chemistry.matchFormalCharge = true;
    assertTrue(query.match(target).found());
    assertFalse(query.withOptions(options).match(target).found());
  }

  @Test
  void compiledQueriesCanBeReusedConcurrently() throws Exception {
    Substructure query = Substructure.smarts("[$(C(=O)[O,N])]");
    ExecutorService workers = Executors.newFixedThreadPool(4);
    try {
      List<Future<Boolean>> results = new ArrayList<>();
      for (int call = 0; call < 12; call++) {
        final boolean expected = call % 2 == 0;
        results.add(workers.submit(() -> query.match(molecule(expected ? "CC(=O)N" : "CCO")).found() == expected));
      }
      for (Future<Boolean> result : results) assertTrue(result.get(2, TimeUnit.SECONDS));
    } finally { workers.shutdownNow(); }
  }

  @Test
  void emptyQueryIsACompleteEmbedding() throws Exception {
    MolGraph empty = new MolGraph.Builder().atomCount(0).atomicNumbers(new int[0]).neighbors(new int[0][]).build();
    SearchResult result = Substructure.find(empty).match(molecule("C"));
    assertTrue(result.found());
    assertTrue(result.complete());
    assertTrue(result.matches().get(0).atomMapping().isEmpty());
  }

  @Test
  void decomposesAttachedGroupsWithoutLosingAtomProperties() throws Exception {
    List<Substructure.Decomposition> groups = Substructure.decomposeRGroups(molecule("c1ccccc1"),
        Arrays.asList(molecule("[13CH3]c1ccccc1"), molecule("Nc1ccccc1"), molecule("CCCCCC")),
        new ChemOptions(), 1000);
    assertEquals(6, groups.get(0).core().getAtomCount());
    assertEquals(13, groups.get(0).rgroups().get("R1").getAtom(0).getMassNumber());
    assertEquals(7, groups.get(1).rgroups().get("R1").getAtom(0).getAtomicNumber());
    assertNull(groups.get(2).core());assertTrue(groups.get(2).rgroups().isEmpty());
    assertFalse(groups.get(2).search().found());assertEquals(SearchResult.Status.Completed,groups.get(2).search().status());
  }

  @Test
  void decompositionDoesNotTreatTimeoutAsAnAbsentCore() throws Exception {
    Substructure.Decomposition result=Substructure.decomposeRGroups(molecule("CC"),
        Collections.singletonList(molecule("CCC")),new ChemOptions(),0).get(0);
    assertNull(result.core());assertTrue(result.rgroups().isEmpty());
    assertEquals(SearchResult.Status.TimedOut,result.search().status());assertFalse(result.search().complete());
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.Duration;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;

class SubstructureContractTest {
  @ParameterizedTest
  @CsvSource(value = {
      "[CX3](=O)[NX3]|CC(=O)Nc1ccc(O)cc1|true",
      "[#6;R]|CCO|false",
      "[#6;R]|C1CCCCC1|true",
      "[$(C(=O)N)]|CC(=O)N|true",
      "[C;!$(C=O)]|C=O|false",
      "[O;!$([OH])]|CC(=O)Nc1ccc(O)cc1|true",
      "[O;!$([OH])]|CO|false",
      "[#7,#8,#16;-0,-1;!$([o,s,nX3]);!$([Nv5,Pv4,Pv5,Sv4,Sv6])]|CC(=O)Nc1ccc(O)cc1|true",
      "[Cl]C|CCO|false",
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
    IAtomContainer parsed = rawMolecule(target);
    assertEquals(expected, SmartsPattern.create(query).matches(parsed.clone()));
    SearchResult result = Substructure.smarts(query).withTimeout(1000).match(parsed);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertEquals(expected, result.found());
    assertTrue(result.complete());
    assertTrue(result.optimal());
  }

  @Test
  void firstMatchDoesNotEnumerateDisconnectedEmbeddings() throws Exception {
    Substructure.smarts("C").match(rawMolecule("CC"));
    String chain = String.join("", Collections.nCopies(40, "C"));
    Substructure query = Substructure.smarts("C.C.C.C").withTimeout(100);
    IAtomContainer target = rawMolecule(chain);
    SearchResult result = assertTimeout(Duration.ofSeconds(1), () -> query.match(target));
    assertTrue(result.found());
    assertEquals(1, result.matches().size());
    assertEquals(4, result.matches().get(0).atomCount());
    assertEquals(4, result.matches().get(0).fragments().size());
  }

  @Test
  void distinguishesExhaustionFromTheResultLimit() throws Exception {
    Substructure query = Substructure.smarts("O");
    SearchResult one = query.matchAll(rawMolecule("CO"), 1);
    assertEquals(SearchResult.Status.Completed, one.status());
    assertTrue(one.complete());
    SearchResult limited = query.matchAll(rawMolecule("COCO"), 1);
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertEquals(1, limited.matches().size());
    assertTrue(limited.optimal());
    assertFalse(limited.complete());
  }

  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void sharesTimeoutAndNodeLimitsAcrossBothQueryTypes(boolean smarts) throws Exception {
    MolGraph target = new MolGraph(rawMolecule("CCCCCCCCCCCC"));
    Substructure query = smarts ? Substructure.smarts("CCCC") : Substructure.find(new MolGraph(rawMolecule("CCCC")));
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
        .match(rawMolecule(String.join("", Collections.nCopies(100, "C"))));
    assertEquals(SearchResult.Status.Cancelled, result.status());
    assertTrue(result.nodesVisited() > 1);
    assertFalse(result.found());
    assertFalse(result.complete());
  }

  @Test
  void reportsSelectedBondsAndTheirFragments() throws Exception {
    SearchResult result = Substructure.find(rawMolecule("CC.O")).match(rawMolecule("CCO"));
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
    Substructure query = smarts ? Substructure.smarts("CCC") : Substructure.find(rawMolecule("CCC"));
    IAtomContainer triangle = rawMolecule("C1CC1");
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

  /** Induced SMARTS embeddings of small graphs equal those of an exhaustive check. */
  @ParameterizedTest
  @ValueSource(booleans = {false, true})
  void smartsEmbeddingsAgreeWithAnExhaustiveCheck(boolean induced) {
    Substructure.Options options = new Substructure.Options();
    options.induced = induced;
    for (int queryAtoms = 1; queryAtoms <= 4; queryAtoms++)
      for (int queryMask = 0; queryMask < 1 << pairs(queryAtoms); queryMask++) {
        MolGraph query = graph(queryAtoms, queryMask);
        Substructure prepared = Substructure.smarts(pattern(query)).withOptions(options);
        for (int targetAtoms = queryAtoms; targetAtoms <= 4; targetAtoms++)
          for (int targetMask = 0; targetMask < 1 << pairs(targetAtoms); targetMask++) {
            MolGraph target = graph(targetAtoms, targetMask);
            SearchResult result = prepared.matchAll(target);
            Set<List<Integer>> actual = new HashSet<>();
            for (SearchResult.Match match : result.matches())
              actual.add(new ArrayList<>(match.atomMapping().values()));
            assertEquals(SearchResult.Status.Completed, result.status());
            assertEquals(embeddings(query, target, induced), actual,
                "query=" + queryAtoms + "/" + queryMask + ", target=" + targetAtoms + "/" + targetMask);
          }
      }
  }

  /** A SMARTS pattern that keeps the graph's atom order: dot-separated atoms bonded by ring closures. */
  private static String pattern(MolGraph graph) {
    List<StringBuilder> closures = new ArrayList<>();
    for (int atom = 0; atom < graph.n; atom++) closures.add(new StringBuilder());
    int label = 10;
    for (int atom = 0; atom < graph.n; atom++)
      for (int neighbor : graph.neighbors[atom])
        if (atom < neighbor) {
          closures.get(atom).append('%').append(label);
          closures.get(neighbor).append('%').append(label++);
        }
    StringBuilder smarts = new StringBuilder();
    for (int atom = 0; atom < graph.n; atom++) smarts.append(atom == 0 ? "" : ".").append("[#6]").append(closures.get(atom));
    return smarts.toString();
  }

  private static Set<List<Integer>> embeddings(MolGraph query, MolGraph target, boolean induced) {
    Set<List<Integer>> found = new HashSet<>();
    int[] mapping = new int[query.n];
    boolean[] used = new boolean[target.n];
    extend(query, target, induced, 0, mapping, used, found);
    return found;
  }

  private static void extend(MolGraph query, MolGraph target, boolean induced, int atom, int[] mapping,
      boolean[] used, Set<List<Integer>> found) {
    if (atom == query.n) {
      List<Integer> images = new ArrayList<>();
      for (int image : mapping) images.add(image);
      found.add(images);
      return;
    }
    for (int image = 0; image < target.n; image++) {
      if (used[image]) continue;
      boolean valid = true;
      for (int previous = 0; previous < atom && valid; previous++) {
        boolean queryBond = query.hasBond(atom, previous), targetBond = target.hasBond(image, mapping[previous]);
        valid = (!queryBond || targetBond) && (!induced || queryBond == targetBond);
      }
      if (!valid) continue;
      mapping[atom] = image;
      used[image] = true;
      extend(query, target, induced, atom + 1, mapping, used, found);
      used[image] = false;
    }
  }

  private static int pairs(int atoms) {
    return atoms * (atoms - 1) / 2;
  }

  private static MolGraph graph(int atoms, int mask) {
    List<List<Integer>> rows = new ArrayList<>();
    for (int atom = 0; atom < atoms; atom++) rows.add(new ArrayList<>());
    int bit = 0;
    for (int atom = 0; atom < atoms; atom++)
      for (int neighbor = atom + 1; neighbor < atoms; neighbor++)
        if ((mask & (1 << bit++)) != 0) {
          rows.get(atom).add(neighbor);
          rows.get(neighbor).add(atom);
        }
    int[][] neighbors = rows.stream().map(row -> row.stream().mapToInt(Integer::intValue).toArray())
        .toArray(int[][]::new);
    int[] elements = new int[atoms];
    Arrays.fill(elements, 6);
    return new MolGraph.Builder().atomCount(atoms).atomicNumbers(elements).neighbors(neighbors).build();
  }

  /** The induced test reads each atom's neighbours, so a large star SMARTS embeds well within its limit. */
  @Test
  void inducedStarSmartsEmbedsInLinearTime() throws Exception {
    int leaves = 3000;
    String star = "C" + String.join("", Collections.nCopies(leaves - 1, "(C)")) + "C";
    Substructure.Options options = new Substructure.Options();
    options.induced = true;
    options.timeoutMillis = 30_000;
    Substructure query = Substructure.smarts(star).withOptions(options);
    IAtomContainer target = rawMolecule(star);
    // Comparing every pair of atoms through the container would take about 30 s here.
    SearchResult result = assertTimeout(Duration.ofSeconds(5), () -> query.match(target));
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.found());
    assertEquals(leaves + 1, result.matches().get(0).atomCount());
    assertEquals(leaves, result.matches().get(0).bondCount());
    // A bond between two leaves is not part of an induced star.
    assertFalse(Substructure.smarts("C(C)(C)C").withOptions(options).match(rawMolecule("CC1CC1")).found());
    assertTrue(Substructure.smarts("C(C)(C)C").match(rawMolecule("CC1CC1")).found());
  }

  /** Every benchmark pattern with a negation compiles where CDK compiles it and finds what CDK finds. */
  @Test
  void negatedBenchmarkPatternsMatchAsCdkDoes() throws Exception {
    Path data = Paths.get("benchmarks", "data", "ehrlich_rarey_smarts.txt");
    assumeTrue(Files.exists(data), "the benchmark patterns are not available");
    String[] targets = {"CC(=O)Nc1ccc(O)cc1", "OC(=O)c1ccccc1OC(C)=O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "C[NH3+].[O-]C(=O)CC(F)(F)F", "c1ccc2[nH]ccc2c1", "NS(=O)(=O)c1ccc(Cl)cc1", "C=CC#N",
        "CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1", "O=C1NC(=O)C(N1)=O", "CC(C)(C)OC(=O)N[C@@H](Cc1ccccc1)C(=O)O"};
    int compared = 0;
    for (String line : Files.readAllLines(data, StandardCharsets.UTF_8)) {
      if (line.isEmpty() || line.startsWith("#")) continue;
      String pattern = line.split("\t")[0];
      if (!pattern.contains("!")) continue;
      SmartsPattern reference;
      try {
        reference = SmartsPattern.create(pattern);
      } catch (IllegalArgumentException unsupported) {
        continue;
      }
      Substructure query = Substructure.smarts(pattern).withTimeout(10_000);
      for (String smiles : targets) {
        IAtomContainer target = rawMolecule(smiles);
        SearchResult result = query.match(target);
        assertEquals(SearchResult.Status.Completed, result.status(), pattern + " in " + smiles);
        assertEquals(reference.matches(target.clone()), result.found(), pattern + " in " + smiles);
        compared++;
      }
    }
    assertTrue(compared > 3000, "compared " + compared);
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
    IAtomContainer target = rawMolecule("[NH4+]");
    Substructure query = Substructure.find(rawMolecule("N")).withOptions(options);
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
        results.add(workers.submit(() -> query.match(rawMolecule(expected ? "CC(=O)N" : "CCO")).found() == expected));
      }
      for (Future<Boolean> result : results) assertTrue(result.get(2, TimeUnit.SECONDS));
    } finally { workers.shutdownNow(); }
  }

  @Test
  void emptyQueryIsACompleteEmbedding() throws Exception {
    MolGraph empty = new MolGraph.Builder().atomCount(0).atomicNumbers(new int[0]).neighbors(new int[0][]).build();
    SearchResult result = Substructure.find(empty).match(rawMolecule("C"));
    assertTrue(result.found());
    assertTrue(result.complete());
    assertTrue(result.matches().get(0).atomMapping().isEmpty());
  }

  @Test
  void decomposesAttachedGroupsWithoutLosingAtomProperties() throws Exception {
    List<Substructure.Decomposition> groups = Substructure.decomposeRGroups(rawMolecule("c1ccccc1"),
        Arrays.asList(rawMolecule("[13CH3]c1ccccc1"), rawMolecule("Nc1ccccc1"), rawMolecule("CCCCCC")),
        new ChemOptions(), 1000);
    assertEquals(6, groups.get(0).core().getAtomCount());
    assertEquals(13, groups.get(0).rgroups().get("R1").getAtom(0).getMassNumber());
    assertEquals(7, groups.get(1).rgroups().get("R1").getAtom(0).getAtomicNumber());
    assertNull(groups.get(2).core());assertTrue(groups.get(2).rgroups().isEmpty());
    assertFalse(groups.get(2).search().found());assertEquals(SearchResult.Status.Completed,groups.get(2).search().status());
  }

  @Test
  void decompositionCoreKeepsOnlyTheImagesOfTheQueryBonds() throws Exception {
    // A chain core matched in a ring is a chain: the ring-closing target bond is not a core bond.
    Substructure.Decomposition chain = Substructure.decomposeRGroups(rawMolecule("CCCCCC"),
        Collections.singletonList(rawMolecule("C1CCCCC1C")), new Substructure.Options()).get(0);
    assertTrue(chain.search().found());
    assertEquals(6, chain.core().getAtomCount());
    assertEquals(5, chain.core().getBondCount());
    assertEquals(Collections.singleton("R1"), chain.rgroups().keySet());
    assertEquals(1, chain.rgroups().get("R1").getAtomCount());
    // The options form carries the induced flag, under which a four-atom chain does not embed in cyclobutane.
    Substructure.Options induced = new Substructure.Options();
    induced.induced = true;
    Substructure.Decomposition none = Substructure.decomposeRGroups(rawMolecule("CCCC"),
        Collections.singletonList(rawMolecule("C1CCC1")), induced).get(0);
    assertFalse(none.search().found());
    assertNull(none.core());
    Substructure.Decomposition ring = Substructure.decomposeRGroups(rawMolecule("CCCC"),
        Collections.singletonList(rawMolecule("C1CCC1")), new Substructure.Options()).get(0);
    assertEquals(3, ring.core().getBondCount());
    assertTrue(ring.rgroups().isEmpty());
    assertEquals(6, Substructure.decomposeRGroups(rawMolecule("c1ccccc1"),
        Collections.singletonList(rawMolecule("Cc1ccccc1")), new ChemOptions(), 1000).get(0).core().getBondCount());
  }

  @Test
  void decompositionDoesNotTreatTimeoutAsAnAbsentCore() throws Exception {
    Substructure.Decomposition result=Substructure.decomposeRGroups(rawMolecule("CC"),
        Collections.singletonList(rawMolecule("CCC")),new ChemOptions(),0).get(0);
    assertNull(result.core());assertTrue(result.rgroups().isEmpty());
    assertEquals(SearchResult.Status.TimedOut,result.search().status());assertFalse(result.search().complete());
  }
}

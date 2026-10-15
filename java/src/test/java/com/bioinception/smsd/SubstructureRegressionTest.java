/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.Substructure;
import com.bioinception.smsd.core.SearchResult;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

/** Regression coverage for pruning, cache isolation, and complete enumeration. */
class SubstructureRegressionTest {
  private static final long TIMEOUT_MS = 5_000;

  private static MolGraph trisaccharide() throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(
        "OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1OC1OC(CO)C(O)C(O)C1O"));
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void largeConnectedQueryMatchesItself(int pruning) throws Exception {
    assertTrue(Substructure.find(trisaccharide()).withOptions(searchOptions(options(pruning))).withTimeout(1_000).match(trisaccharide()).found());
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void largeConnectedQueryEnumeratesACompleteSelfMatch(int pruning) throws Exception {
    MolGraph query = trisaccharide();
    MolGraph target = trisaccharide();
    ChemOptions chemistry = options(pruning);
    List<Map<Integer, Integer>> mappings = Substructure.find(query).withOptions(searchOptions(chemistry)).withTimeout(1_000).matchAll(target, 1).matches().stream().map(SearchResult.Match::atomMapping).collect(java.util.stream.Collectors.toList());
    assertEquals(1, mappings.size());
    assertEquals(query.atomCount(), mappings.get(0).size());
    assertTrue(validationErrors(query,target,mappings.get(0),chemistry).isEmpty());
  }

  private static ChemOptions options(int pruning) {
    ChemOptions options = new ChemOptions();
    options.useTwoHopNLF = pruning > 0;
    options.useThreeHopNLF = pruning > 1;
    return options;
  }

  private static MolGraph graph(int[] elements, int[][] neighbors) {
    return new MolGraph.Builder().atomCount(elements.length)
        .atomicNumbers(elements).neighbors(neighbors).build();
  }

  private static MolGraph isolatedCarbons(int count) {
    int[] elements = new int[count];
    Arrays.fill(elements, 6);
    return graph(elements, new int[count][]);
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void mutableOptionsDoNotReuseIncompatibleDomains(int pruning) {
    MolGraph query = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {{}}).formalCharges(new int[] {1})
        .massNumbers(new int[] {13}).tetrahedralChirality(new int[] {1}).build();
    MolGraph target = new MolGraph.Builder().atomCount(1).atomicNumbers(new int[] {6})
        .neighbors(new int[][] {{}}).formalCharges(new int[] {0})
        .massNumbers(new int[] {14}).tetrahedralChirality(new int[] {2})
        .ringFlags(new boolean[] {true}).aromaticFlags(new boolean[] {true}).build();
    ChemOptions options = options(pruning);

    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.matchFormalCharge = true;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.matchFormalCharge = false;
    options.matchIsotope = true;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.matchIsotope = false;
    options.useChirality = true;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.useChirality = false;
    options.ringMatchesRingOnly = true;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.ringMatchesRingOnly = false;
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void relaxedAtomTypesAlsoRelaxNeighborLabels(int pruning) {
    MolGraph query = graph(new int[] {6, 8}, new int[][] {{1}, {0}});
    MolGraph target = graph(new int[] {6, 7}, new int[][] {{1}, {0}});
    ChemOptions options = options(pruning);
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.matchAtomType = false;
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    assertEquals(2, Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).matchAll(target, 10).matches().stream().map(SearchResult.Match::atomMapping).collect(java.util.stream.Collectors.toList()).size());
    options.matchAtomType = true;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void flexibleAromaticityAlsoRelaxesNeighborLabels(int pruning) {
    int[][] neighbors = {{1}, {0}};
    MolGraph query = graph(new int[] {6, 6}, neighbors);
    MolGraph target = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(neighbors).aromaticFlags(new boolean[] {true, true}).build();
    ChemOptions options = options(pruning);
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void extraTargetEdgesCannotInvalidateMultiHopPruning(int pruning) {
    // In a clique all mapped path neighbors are distance 1. Exact distance
    // shells at hops 2 and 3 therefore cannot be used as rejection tests.
    int count = 22;
    int[] elements = new int[count];
    Arrays.fill(elements, 6);
    int[][] path = new int[count][], clique = new int[count][];
    for (int i = 0; i < count; i++) {
      path[i] = i == 0 ? new int[] {1} : i == count - 1 ? new int[] {i - 1}
          : new int[] {i - 1, i + 1};
      clique[i] = new int[count - 1];
      for (int j = 0, k = 0; j < count; j++) if (j != i) clique[i][k++] = j;
    }
    MolGraph query = graph(elements, path), target = graph(elements, clique);
    ChemOptions options = options(pruning);
    options.useThreeHopNLF = true;
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    List<Map<Integer, Integer>> maps = Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).matchAll(target, 1).matches().stream().map(SearchResult.Match::atomMapping).collect(java.util.stream.Collectors.toList());
    assertEquals(1, maps.size());
    assertEquals(count, maps.get(0).size());
    assertTrue(validationErrors(query,target,maps.get(0),options).isEmpty());
    options.useTwoHopNLF = false;
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void enumerationIncludesTargetsPast4096(int pruning) {
    List<Map<Integer, Integer>> maps = Substructure.find(isolatedCarbons(1)).withOptions(searchOptions(options(pruning))).withTimeout(TIMEOUT_MS).matchAll(isolatedCarbons(4_200), 5_000).matches().stream().map(SearchResult.Match::atomMapping).collect(java.util.stream.Collectors.toList());
    assertEquals(4_200, maps.size());
    assertEquals(4_199, maps.get(maps.size() - 1).get(0));
  }

  @ParameterizedTest
  @ValueSource(ints = {0, 1, 2})
  void anyBondOrderStillChecksStrictAromaticity(int pruning) {
    MolGraph query = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}}).aromaticFlags(new boolean[] {true, true})
        .bondAromaticFlags(new boolean[][] {{true}, {true}}).build();
    MolGraph target = new MolGraph.Builder().atomCount(2).atomicNumbers(new int[] {6, 6})
        .neighbors(new int[][] {{1}, {0}}).aromaticFlags(new boolean[] {true, true})
        .bondAromaticFlags(new boolean[][] {{false}, {false}}).build();
    ChemOptions options = options(pruning);
    options.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    assertFalse(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).match(target).found());
    assertTrue(Substructure.find(query).withOptions(searchOptions(options)).withTimeout(TIMEOUT_MS).matchAll(target, 10).matches().stream().map(SearchResult.Match::atomMapping).collect(java.util.stream.Collectors.toList()).isEmpty());
  }
  private static Substructure.Options searchOptions(ChemOptions chemistry) {
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry;
    return options;
  }

  private static List<String> validationErrors(MolGraph query, MolGraph target,
                                               Map<Integer, Integer> mapping, ChemOptions chemistry) {
    List<String> errors = new java.util.ArrayList<>();
    if (new java.util.HashSet<>(mapping.values()).size() != mapping.size()) errors.add("mapping is not injective");
    for (Map.Entry<Integer, Integer> atom : mapping.entrySet()) {
      if (chemistry.matchAtomType && query.atomicNum[atom.getKey()] != target.atomicNum[atom.getValue()])
        errors.add("atom type differs");
      for (int neighbor : mapping.keySet())
        if (query.hasBond(atom.getKey(), neighbor) && !target.hasBond(atom.getValue(), mapping.get(neighbor)))
          errors.add("query bond absent");
    }
    return errors;
  }
}

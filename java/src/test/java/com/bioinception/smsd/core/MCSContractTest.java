/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class MCSContractTest {
  private static MolGraph parse(String smiles) throws Exception {
    return new MolGraph(
        new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }
  @ParameterizedTest
  @CsvSource({"C1CC1,CCC", "CCC,C1CC1"})
  void objectivesDistinguishInducedVerticesFromSelectedEdges(String first, String second)
      throws Exception {
    MolGraph query = parse(first), target = parse(second);
    SearchResult induced = MCS.find(query).match(target);
    assertEquals(2, induced.matches().get(0).atomCount());
    assertEquals(1, induced.matches().get(0).bondCount());
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    SearchResult edges = MCS.find(query).withOptions(options).match(target);
    SearchResult.Match match = edges.matches().get(0);
    assertTrue(edges.optimal());
    assertTrue(edges.complete());
    assertEquals(3, match.atomCount());
    assertEquals(2, match.bondCount());
    assertEquals(2, match.toMolecule(query).getBondCount());
    assertEquals(3, new MolGraph(match.toMolecule(query)).atomCount());
    assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    assertTrue(MCS.isMappingMaximal(query, target, match, options));
  }
  @Test
  void mcesDoesNotPadWithIsolatedAtoms() throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    SearchResult result = MCS.find(parse("CC.O")).withOptions(options).match(parse("CC.O"));
    assertEquals(2, result.matches().get(0).atomCount());
    assertEquals(1, result.matches().get(0).bondCount());
    assertFalse(MCS.find(parse("C.O")).withOptions(options).match(parse("C.O")).found());
    options.atomWeights = new double[] {1, 1, 1};
    assertThrows(
        IllegalArgumentException.class, () -> MCS.find(parse("CC.O")).withOptions(options));
  }
  @Test
  void selectedEdgesDefineFragmentsAndOutput() throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    SearchResult.Match match =
        MCS.find(parse("CC.CC")).withOptions(options).match(parse("CCCC")).matches().get(0);
    assertEquals(2, match.fragments().size());
    assertEquals(2, match.bondCount());
    assertEquals(2, match.toMolecule(parse("CC.CC")).getBondCount());
    assertTrue(match.toSmiles(parse("CC.CC")).contains("."));
  }
  @Test
  void ringsAreCompleteOnBothSidesAndFusionAtomsMaySurvive() throws Exception {
    MCS.Options options = new MCS.Options();
    options.chemistry.completeRingsOnly = true;
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph benzene = parse("c1ccccc1"), naphthalene = parse("c1ccc2ccccc2c1");
      SearchResult result = MCS.find(reverse ? naphthalene : benzene)
                                .withOptions(options)
                                .match(reverse ? benzene : naphthalene);
      assertEquals(6, result.matches().get(0).atomCount());
      assertTrue(result.optimal());
      SearchResult empty = MCS.find(reverse ? parse("CCCCCC") : benzene)
                               .withOptions(options)
                               .match(reverse ? benzene : parse("CCCCCC"));
      assertFalse(empty.found());
      assertTrue(empty.optimal());
    }
  }
  @Test
  void ringConstraintsDoNotEraseAnAcyclicOptimum() throws Exception {
    MCS.Options options = new MCS.Options();
    options.chemistry.ringMatchesRingOnly = true;
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph query = parse("c1ccccc1CCCC"), target = parse("CCCCCC");
      assertEquals(4,
          MCS.find(reverse ? target : query)
              .withOptions(options)
              .match(reverse ? query : target)
              .matches()
              .get(0)
              .atomCount());
    }
  }
  @Test
  void fragmentBoundsArePartOfOptimisation() throws Exception {
    MCS.Options options = new MCS.Options();
    options.connectivity = MCS.Connectivity.Fragments;
    options.minFragmentAtoms = 3;
    options.maxFragments = 1;
    SearchResult result =
        MCS.find(parse("CCCC.OCCO")).withOptions(options).match(parse("CCCCC.OCCCO"));
    assertEquals(4, result.matches().get(0).atomCount());
    assertEquals(1, result.matches().get(0).fragments().size());
    assertFalse(MCS.find(parse("CCO")).withOptions(options).match(parse("C")).found());
  }
  @Test
  void resultLimitAndBudgetStopsHaveDifferentProofFacts() throws Exception {
    SearchResult cap = MCS.find(parse("CC")).matchAll(parse("CCC"), 1);
    assertEquals(SearchResult.Status.LimitReached, cap.status());
    assertTrue(cap.optimal());
    assertFalse(cap.complete());
    assertEquals(1, cap.matches().size());
    SearchResult deadline = MCS.find(parse("CC")).withTimeout(0).match(parse("CCC"));
    assertEquals(SearchResult.Status.TimedOut, deadline.status());
    assertFalse(deadline.optimal());
    assertFalse(deadline.complete());
    assertEquals(0, deadline.nodesVisited());
    MCS.Options options = new MCS.Options();
    options.cancellation = () -> true;
    assertEquals(SearchResult.Status.Cancelled,
        MCS.find(parse("CC")).withOptions(options).match(parse("CCC")).status());
    options.cancellation = null;
    options.nodeLimit = 3;
    SearchResult limited = MCS.find(parse("CCO")).withOptions(options).match(parse("CCCO"));
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertFalse(limited.optimal());
    assertTrue(limited.found());
    assertTrue(MCS.validateMapping(parse("CCO"), parse("CCCO"), limited.matches().get(0), options)
            .isEmpty());
  }
  @Test
  void preparedOptionsAreCopiedAndOriginalTargetIndicesSurvive() throws Exception {
    MolGraph query = parse("CC"), target = parse("CCCC");
    MCS.Options options = new MCS.Options();
    options.excludedTargetAtoms.add(0);
    MCS prepared = MCS.find(query).withOptions(options);
    options.excludedTargetAtoms.clear();
    assertFalse(prepared.match(target).matches().get(0).atomMapping().containsValue(0));
    assertEquals(
        2, MCS.matchConstrained(Arrays.asList(query, query), target, new MCS.Options()).size());
    List<SearchResult> constrained =
        MCS.matchConstrained(Arrays.asList(query, query), target, new MCS.Options());
    Set<Integer> used = new HashSet<>();
    for (SearchResult result : constrained)
      for (int atom : result.matches().get(0).atomMapping().values()) assertTrue(used.add(atom));
    assertEquals(4, used.size());
  }
  @Test
  void consensusHonoursSupportThresholdAndReportsHeuristicOptimality() throws Exception {
    List<MolGraph> molecules =
        Arrays.asList(parse("N"), parse("c1ccccc1"), parse("Cc1ccccc1"), parse("Oc1ccccc1"));
    MCS.ConsensusResult result = MCS.consensus(molecules, new MCS.Options(), 0.75);
    assertTrue(result.found());
    assertEquals(6, result.fragment().atomCount());
    assertTrue(result.supportCount() >= 3);
    assertFalse(result.optimal());
    assertEquals(6, result.graph().atomCount());
    assertFalse(result.mappings().containsKey(0));
    assertFalse(MCS.consensus(molecules, new MCS.Options(), 1.0).found());
    MCS.ConsensusResult cut =
        MCS.consensus(Arrays.asList(parse("CCO"), parse("CCN")), new MCS.Options(), 1.0);
    assertTrue(cut.found());
    assertEquals("CC", cut.toSmiles());
  }
  @Test
  void sourceGraphViewsPreserveChemistryAndDeterminedStereo() throws Exception {
    MolGraph query = parse("N[C@@H](C)C(=O)O");
    SearchResult.Match match = MCS.find(query).match(query).matches().get(0);
    MolGraph view = match.toGraph(query);
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = true;
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry;
    assertTrue(Substructure.find(query).withOptions(options).match(view).found());
    assertTrue(Substructure.smarts("N[C@@H](C)C(=O)O").match(view).found());
    assertFalse(Substructure.smarts("N[C@H](C)C(=O)O").match(view).found());
    assertEquals(query.getCanonicalHash(), view.getCanonicalHash());
    assertEquals(
        CIPAssigner.assignRS(query), CIPAssigner.assignRS(parse(view.toCanonicalSmiles())));
    assertEquals(
        CIPAssigner.assignRS(query), CIPAssigner.assignRS(new MolGraph(MCS.toMolecule(view))));
    MolGraph aromatic = parse("c1ccncc1");
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int atom = 0; atom < aromatic.n; atom++)
      if (aromatic.atomicNum[atom] == 6)
        atoms.put(atom, atom);
    MolGraph partial =
        MCS.fragment(aromatic, atoms, MCS.allBonds(aromatic, atoms).keySet()).toGraph(aromatic);
    for (int atom = 0; atom < partial.n; atom++) {
      assertTrue(partial.ring[atom]);
      assertTrue(partial.aromatic[atom]);
    }
  }
  @Test
  void originalTargetExclusionsWorkAcrossLargeAtomIndices() throws Exception {
    MolGraph query = parse("CCCC"), target = parse(String.join("", Collections.nCopies(80, "C")));
    MCS.Options options = new MCS.Options();
    for (int atom = 0; atom < 66; atom++) options.excludedTargetAtoms.add(atom);
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal());
    assertEquals(4, result.matches().get(0).atomCount());
    assertTrue(result.matches().get(0).atomMapping().values().stream().allMatch(atom -> atom >= 66));
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  @Test
  void smallGraphsAgreeWithIndependentVertexEnumeration() {
    for (int queryMask = 0; queryMask < 8; queryMask++)
      for (int targetMask = 0; targetMask < 8; targetMask++)
        for (MCS.Objective objective : MCS.Objective.values())
          for (MCS.Connectivity connectivity : MCS.Connectivity.values()) {
            MolGraph query = graph(queryMask), target = graph(targetMask);
            MCS.Options options = new MCS.Options();
            options.objective = objective;
            options.connectivity = connectivity;
            int[] expected = oracle(query, target, objective, connectivity, new int[] {-1, -1, -1},
                new boolean[3], 0, new int[] {0, 0});
            SearchResult result = MCS.find(query).withOptions(options).match(target);
            int atoms = result.found() ? result.matches().get(0).atomCount() : 0,
                bonds = result.found() ? result.matches().get(0).bondCount() : 0;
            assertArrayEquals(expected, new int[] {atoms, bonds},
                "masks=" + queryMask + "," + targetMask + " " + objective + " " + connectivity);
            assertTrue(result.optimal());
            if (result.found())
              assertTrue(
                  MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
          }
  }
  @Test
  void mcesStereoUsesSelectedCarriersRatherThanEveryMappedEndpoint() throws Exception {
    MolGraph query = parse("CC/C=C/CC"), target = parse("CC/C=C\\CC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    options.chemistry.useBondStereo = true;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    SearchResult.Match match = result.matches().get(0);
    assertTrue(result.optimal());
    assertEquals(4, match.bondCount());
    assertEquals(6, match.atomCount());
    assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    assertEquals(4, match.toMolecule(query).getBondCount());
  }
  @Test
  void largerSmallGraphsAgreeWithIndependentEnumeration() {
    Random random = new Random(5197);
    for (int sample = 0; sample < 80; sample++) {
      int size = sample < 40 ? 4 : 5, possibilities = 1 << (size * (size - 1) / 2);
      MolGraph query = graph(size, random.nextInt(possibilities)),
               target = graph(size, random.nextInt(possibilities));
      for (MCS.Objective objective : MCS.Objective.values())
        for (MCS.Connectivity connectivity : MCS.Connectivity.values()) {
          MCS.Options options = new MCS.Options();
          options.objective = objective;
          options.connectivity = connectivity;
          int[] indexed = new int[size];
          Arrays.fill(indexed, -1);
          int[] expected = oracle(query, target, objective, connectivity, indexed,
              new boolean[size], 0, new int[] {0, 0});
          SearchResult result = MCS.find(query).withOptions(options).match(target);
          int atoms = result.found() ? result.matches().get(0).atomCount() : 0,
              bonds = result.found() ? result.matches().get(0).bondCount() : 0;
          assertArrayEquals(expected, new int[] {atoms, bonds},
              "sample=" + sample + " " + objective + " " + connectivity);
          assertTrue(result.optimal());
          if (result.found())
            assertTrue(
                MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
        }
    }
  }
  private static MolGraph graph(int mask) {
    return graph(3, mask);
  }
  private static MolGraph graph(int size, int mask) {
    List<List<Integer>> rows = new ArrayList<>();
    for (int atom = 0; atom < size; atom++) rows.add(new ArrayList<>());
    int bit = 0;
    for (int first = 0; first < size; first++)
      for (int second = first + 1; second < size; second++)
        if ((mask & (1 << bit++)) != 0) {
          rows.get(first).add(second);
          rows.get(second).add(first);
        }
    int[][] neighbors = rows.stream()
                            .map(row -> row.stream().mapToInt(Integer::intValue).toArray())
                            .toArray(int[][] ::new);
    int[] elements = new int[size];
    Arrays.fill(elements, 6);
    return new MolGraph.Builder()
        .atomCount(size)
        .atomicNumbers(elements)
        .neighbors(neighbors)
        .build();
  }
  private static int[] oracle(MolGraph query, MolGraph target, MCS.Objective objective,
      MCS.Connectivity connectivity, int[] mapping, boolean[] used, int depth, int[] best) {
    if (depth < query.n) {
      oracle(query, target, objective, connectivity, mapping, used, depth + 1, best);
      for (int atom = 0; atom < target.n; atom++)
        if (!used[atom]) {
          mapping[depth] = atom;
          used[atom] = true;
          oracle(query, target, objective, connectivity, mapping, used, depth + 1, best);
          used[atom] = false;
          mapping[depth] = -1;
        }
      return best;
    }
    Map<Integer, Integer> atoms = new TreeMap<>();
    List<SearchResult.Bond> bonds = new ArrayList<>();
    for (int first = 0; first < query.n; first++)
      if (mapping[first] >= 0) {
        atoms.put(first, mapping[first]);
        for (int second = first + 1; second < query.n; second++)
          if (mapping[second] >= 0) {
            boolean queryBond = query.hasBond(first, second),
                    targetBond = target.hasBond(mapping[first], mapping[second]);
            if (objective == MCS.Objective.MCIS && queryBond != targetBond)
              return best;
            if (queryBond && targetBond)
              bonds.add(new SearchResult.Bond(first, second));
          }
      }
    if (objective == MCS.Objective.MCES) {
      Set<Integer> incident = new HashSet<>();
      for (SearchResult.Bond bond : bonds) {
        incident.add(bond.first());
        incident.add(bond.second());
      }
      atoms.keySet().retainAll(incident);
    }
    List<Set<Integer>> fragments = MCS.fragments(atoms, bonds);
    if (connectivity == MCS.Connectivity.Connected && fragments.size() > 1)
      return best;
    int atomCount = atoms.size(), bondCount = bonds.size();
    int primary = objective == MCS.Objective.MCIS ? atomCount : bondCount,
        oldPrimary = objective == MCS.Objective.MCIS ? best[0] : best[1];
    int secondary = objective == MCS.Objective.MCIS ? bondCount : atomCount,
        oldSecondary = objective == MCS.Objective.MCIS ? best[1] : best[0];
    if (primary > oldPrimary || (primary == oldPrimary && secondary > oldSecondary)) {
      best[0] = atomCount;
      best[1] = bondCount;
    }
    return best;
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

class MCSContractTest {
  @Test
  void validateMappingReportsASelectedBondOutsideTheAtomMapping() throws Exception {
    MolGraph query = rawGraph("CCO"), target = rawGraph("CCO");
    Map<Integer, Integer> atoms = new TreeMap<>();
    atoms.put(0, 0);
    atoms.put(1, 1);
    Map<SearchResult.Bond, SearchResult.Bond> bonds = new LinkedHashMap<>();
    bonds.put(new SearchResult.Bond(0, 1), new SearchResult.Bond(0, 1));
    bonds.put(new SearchResult.Bond(1, 2), new SearchResult.Bond(1, 2));
    // A match built from its parts, as one read from JSON is, may select a bond at an unmapped atom.
    SearchResult.Match match = new SearchResult.Match(atoms, bonds, Collections.emptyList(), 2);
    assertEquals(Collections.singletonList("Selected bond is outside the atom mapping"),
        MCS.validateMapping(query, target, match, new MCS.Options()));
  }

  @ParameterizedTest
  @CsvSource({"C1CC1,CCC", "CCC,C1CC1"})
  void objectivesDistinguishInducedVerticesFromSelectedEdges(String first, String second)
      throws Exception {
    MolGraph query = rawGraph(first), target = rawGraph(second);
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
    SearchResult result = MCS.find(rawGraph("CC.O")).withOptions(options).match(rawGraph("CC.O"));
    assertEquals(2, result.matches().get(0).atomCount());
    assertEquals(1, result.matches().get(0).bondCount());
    assertFalse(MCS.find(rawGraph("C.O")).withOptions(options).match(rawGraph("C.O")).found());
    options.atomWeights = new double[] {1, 1, 1};
    assertThrows(
        IllegalArgumentException.class, () -> MCS.find(rawGraph("CC.O")).withOptions(options));
  }
  @Test
  void selectedEdgesDefineFragmentsAndOutput() throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    SearchResult.Match match =
        MCS.find(rawGraph("CC.CC")).withOptions(options).match(rawGraph("CCCC")).matches().get(0);
    assertEquals(2, match.fragments().size());
    assertEquals(2, match.bondCount());
    assertEquals(2, match.toMolecule(rawGraph("CC.CC")).getBondCount());
    assertTrue(match.toSmiles(rawGraph("CC.CC")).contains("."));
  }
  @ParameterizedTest
  @CsvSource({"C1CC1,CCC,Connected,6", "CCC,CCCC,Connected,4",
      "CC.CC,CCCC,Fragments,8", "C1CCC1,C1CCC1,Connected,8", "CC=O,O=CCC,Connected,1",
      "C.C,CC,Connected,0", "CC,C.C,Fragments,0", "CC.C.C,CC.C.C,Fragments,2",
      "CCC.O,CC.O,Fragments,4"})
  void mcesEnumerationRetainsAllOptimalCorrespondences(
      String querySmiles, String targetSmiles, MCS.Connectivity connectivity, int expected)
      throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = connectivity;
    SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 1000);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(expected, result.matches().size());
    Set<String> correspondences = new HashSet<>();
    for (SearchResult.Match match : result.matches()) {
      assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
      assertEquals(result.score(), match.score());
      assertTrue(correspondences.add(match.atomMapping().toString() + match.bondMapping()));
    }
  }
  @ParameterizedTest
  @CsvSource({"1,LimitReached,false", "5,LimitReached,false", "6,Completed,true"})
  void mcesEnumerationDistinguishesAnExactResultCountFromTruncation(
      int limit, SearchResult.Status status, boolean complete) throws Exception {
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    SearchResult result = MCS.find(rawGraph("C1CC1")).withOptions(options).matchAll(rawGraph("CCC"), limit);
    assertEquals(status, result.status());
    assertTrue(result.optimal());
    assertEquals(complete, result.complete());
    assertEquals(limit, result.matches().size());
    for (SearchResult.Match match : result.matches())
      assertEquals(2, match.bondCount());
  }
  @Test
  void ringsAreCompleteOnBothSidesAndFusionAtomsMaySurvive() throws Exception {
    MCS.Options options = new MCS.Options();
    options.chemistry.completeRingsOnly = true;
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph benzene = rawGraph("c1ccccc1"), naphthalene = rawGraph("c1ccc2ccccc2c1");
      SearchResult result = MCS.find(reverse ? naphthalene : benzene)
                                .withOptions(options)
                                .match(reverse ? benzene : naphthalene);
      assertEquals(6, result.matches().get(0).atomCount());
      assertTrue(result.optimal());
      SearchResult empty = MCS.find(reverse ? rawGraph("CCCCCC") : benzene)
                               .withOptions(options)
                               .match(reverse ? benzene : rawGraph("CCCCCC"));
      assertFalse(empty.found());
      assertTrue(empty.optimal());
    }
  }
  @Test
  void ringConstraintsDoNotEraseAnAcyclicOptimum() throws Exception {
    MCS.Options options = new MCS.Options();
    options.chemistry.ringMatchesRingOnly = true;
    for (boolean reverse : new boolean[] {false, true}) {
      MolGraph query = rawGraph("c1ccccc1CCCC"), target = rawGraph("CCCCCC");
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
        MCS.find(rawGraph("CCCC.OCCO")).withOptions(options).match(rawGraph("CCCCC.OCCCO"));
    assertEquals(4, result.matches().get(0).atomCount());
    assertEquals(1, result.matches().get(0).fragments().size());
    assertFalse(MCS.find(rawGraph("CCO")).withOptions(options).match(rawGraph("C")).found());
  }
  @Test
  void resultLimitAndBudgetStopsHaveDifferentProofFacts() throws Exception {
    SearchResult cap = MCS.find(rawGraph("CC")).matchAll(rawGraph("CCC"), 1);
    assertEquals(SearchResult.Status.LimitReached, cap.status());
    assertTrue(cap.optimal());
    assertFalse(cap.complete());
    assertEquals(1, cap.matches().size());
    SearchResult deadline = MCS.find(rawGraph("CC")).withTimeout(0).match(rawGraph("CCC"));
    assertEquals(SearchResult.Status.TimedOut, deadline.status());
    assertFalse(deadline.optimal());
    assertFalse(deadline.complete());
    assertEquals(0, deadline.nodesVisited());
    MCS.Options options = new MCS.Options();
    options.cancellation = () -> true;
    assertEquals(SearchResult.Status.Cancelled,
        MCS.find(rawGraph("CC")).withOptions(options).match(rawGraph("CCC")).status());
    options.cancellation = null;
    options.nodeLimit = 3;
    SearchResult limited = MCS.find(rawGraph("CCO")).withOptions(options).match(rawGraph("CCCO"));
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertFalse(limited.optimal());
    assertTrue(limited.found());
    assertTrue(MCS.validateMapping(rawGraph("CCO"), rawGraph("CCCO"), limited.matches().get(0), options)
            .isEmpty());
  }
  @Test
  void preparedOptionsAreCopiedAndOriginalTargetIndicesSurvive() throws Exception {
    MolGraph query = rawGraph("CC"), target = rawGraph("CCCC");
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
  void constrainedSearchReturnsOneResultPerQueryAfterCancellation() throws Exception {
    MolGraph query = rawGraph("CC"), target = rawGraph("CCCC");
    MCS.Options options = new MCS.Options();
    options.cancellation = () -> true;
    List<SearchResult> results =
        MCS.matchConstrained(Arrays.asList(query, query, query), target, options);
    assertEquals(3, results.size());
    for (SearchResult result : results) {
      assertEquals(SearchResult.Status.Cancelled, result.status());
      assertFalse(result.found());
    }
  }

  @Test
  @SuppressWarnings("deprecation")
  void mappingFormsFollowTheNativeCore() throws Exception {
    MolGraph query = rawGraph("CCO"), target = rawGraph("CCN");
    Map<Integer, Integer> incompatible = new TreeMap<>();
    incompatible.put(0, 0);
    incompatible.put(1, 1);
    incompatible.put(2, 2);
    MCS.Options options = new MCS.Options();
    // An invalid mapping is reported, and is not maximal, rather than thrown.
    assertEquals(Collections.singletonList("Mapped atoms are incompatible"),
        MCS.validateMapping(query, target, incompatible, options));
    assertFalse(MCS.isMappingMaximal(query, target, incompatible, options));
    assertFalse(MCS.isMappingMaximal(query, target, incompatible, new ChemOptions()));
    assertFalse(MCS.isMappingMaximal(query, target, MCS.mappingMatch(query, target, incompatible, options), options));
    Map<Integer, Integer> backbone = new TreeMap<>(incompatible);
    backbone.remove(2);
    SearchResult.Match match = MCS.mappingMatch(query, target, backbone, options);
    assertEquals(1, match.bondCount());
    assertTrue(MCS.validateMapping(query, target, backbone, options).isEmpty());
    assertTrue(MCS.isMappingMaximal(query, target, backbone, options));
    assertTrue(MCS.isMappingMaximal(query, target, backbone, new ChemOptions()));
    // The map form selects only the bonds the target has, so a missing bond is an induced-topology error.
    MolGraph chain = rawGraph("CCC"), apart = rawGraph("C.CC");
    Map<Integer, Integer> identity = new TreeMap<>();
    for (int atom = 0; atom < 3; atom++) identity.put(atom, atom);
    assertEquals(Arrays.asList("Mapped topology is not induced", "MCIS must select every induced bond",
            "Selected graph is disconnected"),
        MCS.validateMapping(chain, apart, identity, options));
    // The chemistry-only form checks an induced common subgraph in any number of fragments.
    assertEquals(Arrays.asList("Mapped topology is not induced", "MCIS must select every induced bond"),
        MCS.validateMapping(chain, apart, identity, new ChemOptions()));
    assertEquals(Collections.singletonList("Mapped atom index out of range"),
        MCS.validateMapping(chain, apart, Collections.singletonMap(0, 3), options));
    assertThrows(IllegalArgumentException.class, () -> MCS.mappingMatch(chain, apart, Collections.singletonMap(0, 3), options));
    assertThrows(IllegalArgumentException.class, () -> MCS.canonicalizeMapping(chain, apart, Collections.singletonMap(0, 3)));
    assertThrows(IllegalArgumentException.class, () -> MCS.canonicalizeMapping(chain, apart, Collections.singletonMap(-1, 0)));
    assertEquals(identity, MCS.canonicalizeMapping(chain, chain, identity));
  }

  @Test
  void consensusHonoursSupportThresholdAndReportsHeuristicOptimality() throws Exception {
    List<MolGraph> molecules =
        Arrays.asList(rawGraph("N"), rawGraph("c1ccccc1"), rawGraph("Cc1ccccc1"), rawGraph("Oc1ccccc1"));
    MCS.ConsensusResult result = MCS.consensus(molecules, new MCS.Options(), 0.75);
    assertTrue(result.found());
    assertEquals(6, result.fragment().atomCount());
    assertTrue(result.supportCount() >= 3);
    assertFalse(result.optimal());
    assertEquals(6, result.graph().atomCount());
    assertFalse(result.mappings().containsKey(0));
    MCS.ConsensusResult unsupported = MCS.consensus(molecules, new MCS.Options(), 1.0);
    assertFalse(unsupported.found());
    assertEquals(-1, unsupported.sourceIndex());
    assertEquals(0, unsupported.fragment().atomCount());
    assertEquals(0, unsupported.graph().atomCount());
    assertEquals("", unsupported.toSmiles());
    assertTrue(unsupported.mappings().isEmpty());
    MCS.ConsensusResult cut =
        MCS.consensus(Arrays.asList(rawGraph("CCO"), rawGraph("CCN")), new MCS.Options(), 1.0);
    assertTrue(cut.found());
    assertEquals("CC", cut.toSmiles());
  }

  @Test
  void consensusArgumentChecksFollowTheNativeCore() throws Exception {
    // An empty list is not found; exclusions and atom weights have no single target here.
    MCS.ConsensusResult none = MCS.consensus(new ArrayList<>(), new MCS.Options(), 1.0);
    assertFalse(none.found());
    assertEquals(-1, none.sourceIndex());
    assertEquals(0, none.supportCount());
    assertEquals(SearchResult.Status.Completed, none.status());
    List<MolGraph> molecules = Arrays.asList(rawGraph("CCO"), rawGraph("CCN"));
    MCS.Options excluded = new MCS.Options();
    excluded.excludedTargetAtoms.add(0);
    assertThrows(IllegalArgumentException.class, () -> MCS.consensus(molecules, excluded, 1.0));
    MCS.Options weighted = new MCS.Options();
    weighted.atomWeights = new double[] {1, 1, 1};
    assertThrows(IllegalArgumentException.class, () -> MCS.consensus(molecules, weighted, 1.0));
    for (double threshold : new double[] {0.0, 1.5, Double.NaN})
      assertThrows(IllegalArgumentException.class, () -> MCS.consensus(molecules, new MCS.Options(), threshold));
  }
  @Test
  void matchManyOnThreadsKeepsTheSequentialResultsInOrder() throws Exception {
    List<MolGraph> targets = new ArrayList<>();
    for (String smiles : new String[] {"CCO", "CCN", "c1ccccc1O", "CC(=O)O", "C", "OCCO", "CCCCO", "N"})
      targets.add(rawGraph(smiles));
    MCS mcs = MCS.find(rawGraph("CCO"));
    Substructure substructure = Substructure.find(rawGraph("CO"));
    for (int threads : new int[] {0, 1, 3, 16}) {
      List<SearchResult> common = mcs.matchMany(targets, threads), embedded = substructure.matchMany(targets, threads);
      List<SearchResult> commonInTurn = mcs.matchMany(targets), embeddedInTurn = substructure.matchMany(targets);
      assertEquals(targets.size(), common.size());
      for (int index = 0; index < targets.size(); index++) {
        assertEquals(commonInTurn.get(index).matches(), common.get(index).matches(), threads + " threads");
        assertEquals(embeddedInTurn.get(index).found(), embedded.get(index).found(), threads + " threads");
        assertEquals(embeddedInTurn.get(index).matches(), embedded.get(index).matches(), threads + " threads");
      }
    }
    assertTrue(mcs.matchMany(new ArrayList<>(), 4).isEmpty());
    assertThrows(IllegalArgumentException.class, () -> mcs.matchMany(targets, -1));
    assertThrows(IllegalArgumentException.class, () -> substructure.matchMany(targets, -1));
    assertThrows(NullPointerException.class, () -> substructure.matchMany(Arrays.asList(targets.get(0), null), 2));
  }

  @Test
  void matchManyOnThreadsReadsTheRingsOfColdGraphsWhole() throws Exception {
    // The threads share one prepared query, and each target several times, whose rings are
    // perceived on first use. A thread that saw them half set would drop rings from a
    // complete-ring MCES, so every result must equal that of a search in turn.
    String[] molecules = {"c1ccc2c(c1)ccc1ccccc12", "C1CC2CCC3CCCC4CCC(C1)C2C34", "C1CC2CC1C1CCCCC21",
        "c1ccc2[nH]ccc2c1", "O=C1CCC2CC(O)CCC2C1", "C1CCC2CCCCC2C1", "c1ccccc1CC",
        "OC1CCC2C(C1)CCC1C2CCC2CCCC12", "C1CCC(CC1)C1CCCCC1", "CCCCCC"};
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 60_000;
    String query = "C1CC2C(C1)CCC1C2CCC2CCCCC12";
    String[] patterns = {"[R2]", "[r5]", "[x3;R]", "[R]@[R]@[R]", "[!R]-[R]", "[r6;!r5]~[R2]"};
    for (int round = 0; round < 8; round++) {
      // Each cold target appears three times, so that the threads share its rings as well.
      List<MolGraph> repeated = repeatedCold(molecules), separate = new ArrayList<>();
      for (int copy = 0; copy < 3; copy++)
        for (String smiles : molecules) separate.add(rawGraph(smiles));
      List<SearchResult> threaded = MCS.find(rawGraph(query)).withOptions(options).matchMany(repeated, 8);
      List<SearchResult> inTurn = MCS.find(rawGraph(query)).withOptions(options).matchMany(separate);
      for (int index = 0; index < repeated.size(); index++) {
        assertEquals(SearchResult.Status.Completed, threaded.get(index).status(), "round " + round);
        assertEquals(inTurn.get(index).score(), threaded.get(index).score(), "round " + round);
        assertEquals(inTurn.get(index).matches(), threaded.get(index).matches(), "round " + round);
      }
      for (String pattern : patterns) {
        List<MolGraph> cold = repeatedCold(molecules);
        List<SearchResult> embedded = Substructure.smarts(pattern).matchMany(cold, 8);
        List<SearchResult> embeddedInTurn = Substructure.smarts(pattern).matchMany(separate);
        for (int index = 0; index < cold.size(); index++) {
          assertEquals(embeddedInTurn.get(index).found(), embedded.get(index).found(), pattern);
          assertEquals(embeddedInTurn.get(index).matches(), embedded.get(index).matches(), pattern);
        }
      }
    }
  }

  /** A fresh graph of each molecule, the list holding each three times. */
  private static List<MolGraph> repeatedCold(String[] molecules) throws Exception {
    List<MolGraph> once = new ArrayList<>(), repeated = new ArrayList<>();
    for (String smiles : molecules) once.add(rawGraph(smiles));
    for (int copy = 0; copy < 3; copy++) repeated.addAll(once);
    return repeated;
  }

  @Test
  void consensusReducesWithTheSearchOfMatch() throws Exception {
    // Each reduction is the search match() runs, from the smaller molecule for MCES and by the
    // common-atom search for a fragment MCIS, so a pair reduces to the pairwise optimum.
    MolGraph large = rawGraph("CC(=O)Oc1ccccc1C(=O)O"), small = rawGraph("OC(=O)c1ccccc1");
    for (MCS.Objective objective : MCS.Objective.values()) {
      MCS.Options options = new MCS.Options();
      options.objective = objective;
      options.connectivity = MCS.Connectivity.Fragments;
      SearchResult pair = MCS.find(large).withOptions(options).match(small);
      MCS.ConsensusResult result = MCS.consensus(Arrays.asList(large, small), options, 1.0);
      assertTrue(result.found(), objective.toString());
      assertEquals(pair.matches().get(0).atomCount(), result.fragment().atomCount(), objective.toString());
      assertEquals(pair.matches().get(0).bondCount(), result.fragment().bondCount(), objective.toString());
    }
  }

  @Test
  void sourceGraphViewsPreserveChemistryAndDeterminedStereo() throws Exception {
    MolGraph query = rawGraph("N[C@@H](C)C(=O)O");
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
        CIPAssigner.assignRS(query), CIPAssigner.assignRS(rawGraph(view.toCanonicalSmiles())));
    assertEquals(
        CIPAssigner.assignRS(query), CIPAssigner.assignRS(new MolGraph(MatchGraph.toMolecule(view))));
    MolGraph aromatic = rawGraph("c1ccncc1");
    Map<Integer, Integer> atoms = new TreeMap<>();
    for (int atom = 0; atom < aromatic.n; atom++)
      if (aromatic.atomicNum[atom] == 6)
        atoms.put(atom, atom);
    MolGraph partial =
        MCS.fragment(aromatic, atoms, MatchGraph.allBonds(aromatic, atoms).keySet()).toGraph(aromatic);
    for (int atom = 0; atom < partial.n; atom++) {
      assertTrue(partial.ring[atom]);
      assertTrue(partial.aromatic[atom]);
    }
  }
  @Test
  void originalTargetExclusionsWorkAcrossLargeAtomIndices() throws Exception {
    MolGraph query = rawGraph("CCCC"), target = rawGraph(String.join("", Collections.nCopies(80, "C")));
    MCS.Options options = new MCS.Options();
    for (int atom = 0; atom < 66; atom++) options.excludedTargetAtoms.add(atom);
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal());
    assertEquals(4, result.matches().get(0).atomCount());
    assertTrue(result.matches().get(0).atomMapping().values().stream().allMatch(atom -> atom >= 66));
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  @ParameterizedTest
  @CsvSource({"Connected", "Fragments"})
  void mcesOriginalIndicesAndExclusionsSurviveLeadingIsolatedAtoms(
      MCS.Connectivity connectivity) throws Exception {
    MolGraph query = rawGraph(String.join("", Collections.nCopies(66, "C.")) + "CCC");
    MolGraph target = rawGraph(String.join("", Collections.nCopies(66, "C.")) + "CCC.CCC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = connectivity;
    SearchResult unexcluded = MCS.find(query).withOptions(options).matchAll(target, 10);
    assertEquals(4, unexcluded.matches().size());
    for (int atom = 66; atom < 69; atom++) options.excludedTargetAtoms.add(atom);
    SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 10);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(2, result.matches().size());
    for (SearchResult.Match match : result.matches()) {
      assertEquals(3, match.atomCount());
      assertEquals(2, match.bondCount());
      assertTrue(match.atomMapping().keySet().stream().allMatch(atom -> atom >= 66));
      assertTrue(match.atomMapping().values().stream().allMatch(atom -> atom >= 69));
      assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    }
  }

  @ParameterizedTest
  @CsvSource({"Connected", "Fragments"})
  void mcesCandidateBondsUseOriginalEdgeIndicesAcrossWords(
      MCS.Connectivity connectivity) throws Exception {
    MolGraph query = rawGraph("CCCC");
    MolGraph target = rawGraph(String.join("", Collections.nCopies(66, "NN.")) + "CCC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = connectivity;
    SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 10);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(4, result.matches().size());
    for (SearchResult.Match match : result.matches()) {
      assertEquals(3, match.atomCount());
      assertEquals(2, match.bondCount());
      assertTrue(match.atomMapping().values().stream().allMatch(atom -> atom >= 132));
      assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    }
  }

  @ParameterizedTest
  @CsvSource({"Connected,192,7,7", "Fragments,32,9,10"})
  void mcesProvesAspirinParacetamolOptimumWithinNodeBudget(
      MCS.Connectivity connectivity, int nodeLimit, int bonds, int atoms) throws Exception {
    MolGraph query = rawGraph("CC(=O)Oc1ccccc1C(=O)O"), target = rawGraph("CC(=O)Nc1ccc(O)cc1");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = connectivity;
    options.nodeLimit = nodeLimit;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(bonds, result.score());
    assertEquals(bonds, result.matches().get(0).bondCount());
    assertEquals(atoms, result.matches().get(0).atomCount());
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  @ParameterizedTest
  @CsvSource({"CCC[C@H](F)Cl,C[C@H](F)Cl,tetra,3,4,1",
      "CCC[C@H](F)Cl,C[C@@H](F)Cl,tetra,2,3,3",
      "CCC/C=C/CC,C/C=C\\C,double,2,4,56",
      "c1ccc2ccccc2c1,c1ccccc1,rings,6,6,24",
      "c1ccccc1CC,CCC,rings,1,2,4", "NNO,CC,elements,1,2,4"})
  void mcesOrientationPreservesSelectedChemistry(String querySmiles, String targetSmiles,
      String mode, int bonds, int atoms, int mappings) throws Exception {
    MolGraph query = rawGraph(querySmiles), target = rawGraph(targetSmiles);
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    options.chemistry.useChirality = mode.equals("tetra");
    options.chemistry.useBondStereo = mode.equals("double");
    options.chemistry.completeRingsOnly = mode.equals("rings");
    options.chemistry.matchAtomType = !mode.equals("elements");
    SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 1000);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(bonds, result.score());
    assertEquals(mappings, result.matches().size());
    for (SearchResult.Match match : result.matches()) {
      assertEquals(bonds, match.bondCount());
      assertEquals(atoms, match.atomCount());
      assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    }
  }

  @Test
  void mcesOrientationRetainsOriginalIndicesExclusionsAndStops() throws Exception {
    MolGraph query = rawGraph("NN.CCCC"), target = rawGraph("N.CCC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 4);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(4, result.matches().size());
    for (SearchResult.Match match : result.matches()) {
      assertTrue(match.atomMapping().keySet().stream().allMatch(atom -> atom >= 2));
      assertTrue(match.atomMapping().values().stream().allMatch(atom -> atom >= 1));
      assertEquals(2, match.toMolecule(query).getBondCount());
      assertTrue(MCS.validateMapping(query, target, match, options).isEmpty());
    }
    SearchResult capped = MCS.find(query).withOptions(options).matchAll(target, 1);
    assertEquals(SearchResult.Status.LimitReached, capped.status());
    assertTrue(capped.optimal());
    assertFalse(capped.complete());
    assertEquals(1, capped.matches().size());
    options.excludedTargetAtoms.add(2);
    SearchResult excluded = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Completed, excluded.status());
    assertTrue(excluded.optimal());
    assertFalse(excluded.found());
    options.excludedTargetAtoms.clear();
    options.timeoutMillis = 0;
    SearchResult timedOut = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.TimedOut, timedOut.status());
    assertFalse(timedOut.optimal());
    assertEquals(0, timedOut.nodesVisited());
    options.timeoutMillis = -1;
    options.nodeLimit = 1;
    // The seeds already prove the small pair; a larger reversed pair needs more than one node.
    MolGraph larger = rawGraph("NN.CCCCC1CCCC1"), smaller = rawGraph("N.CCC1CCC1");
    SearchResult limited = MCS.find(larger).withOptions(options).match(smaller);
    assertEquals(SearchResult.Status.LimitReached, limited.status());
    assertFalse(limited.optimal());
    assertEquals(1, limited.nodesVisited());
    assertTrue(MCS.validateMapping(larger, smaller, limited.matches().get(0), options).isEmpty());
    // A node limit reached while the whole query is embedded still leaves the seed.
    MolGraph naphthalene = rawGraph("c1ccc2ccccc2c1"), ethylbenzene = rawGraph("c1ccccc1CC");
    SearchResult embedding = MCS.find(naphthalene).withOptions(options).match(ethylbenzene);
    assertEquals(SearchResult.Status.LimitReached, embedding.status());
    assertFalse(embedding.optimal());
    assertTrue(embedding.found());
    assertTrue(
        MCS.validateMapping(naphthalene, ethylbenzene, embedding.matches().get(0), options).isEmpty());
    options.cancellation = () -> true;
    SearchResult cancelled = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
    assertFalse(cancelled.optimal());
    assertEquals(0, cancelled.nodesVisited());
  }

  /**
   * The greedy seed grows from the query, and from the other molecule it may hold far more
   * common bonds. Here the seed grown from the smaller molecule has 23 common bonds and the one
   * grown from the larger has the optimum, 26. Once its levels of common bonds run out of budget
   * the fragment search starts from the better seed, and proves the optimum within a node limit
   * that the first seed alone needs over two million nodes to reach, whichever molecule is given
   * as the query.
   */
  @Test
  void seedFromTheOtherMoleculeProvesTheChainPair() throws Exception {
    MolGraph query =
        rawGraph("CCCCCCCCCCCCCCSCC(COP(O)(OP(O)(OCC1OC(n2ccc(N)nc2=O)C(O)C1O)=O)=O)OC(=O)CCCCCCCCCCC");
    MolGraph target = rawGraph("CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    options.chemistry.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    options.chemistry.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    options.timeoutMillis = 60_000;
    options.nodeLimit = 100_000;
    for (int turn = 0; turn < 2; turn++) {
      SearchResult result = turn == 0 ? MCS.find(query).withOptions(options).match(target)
                                      : MCS.find(target).withOptions(options).match(query);
      String label = "turn " + turn + " nodes " + result.nodesVisited();
      assertTrue(result.optimal() && !result.matches().isEmpty(), label);
      assertEquals(26, result.matches().get(0).bondCount(), label);
      assertEquals(30, result.matches().get(0).atomCount(), label);
    }
  }

  @Test
  void reversedMcesQueriesCanBeReusedConcurrently() throws Exception {
    MolGraph query = rawGraph("NN.CCCC"), target = rawGraph("N.CCC");
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = MCS.Connectivity.Fragments;
    MCS prepared = MCS.find(query).withOptions(options);
    Set<SearchResult.Match> expected = new HashSet<>(prepared.matchAll(target, 4).matches());
    java.util.concurrent.ExecutorService executor =
        java.util.concurrent.Executors.newFixedThreadPool(2);
    try {
      java.util.concurrent.Future<SearchResult> first =
          executor.submit(() -> prepared.matchAll(target, 4));
      java.util.concurrent.Future<SearchResult> second =
          executor.submit(() -> prepared.matchAll(target, 4));
      for (SearchResult result : Arrays.asList(first.get(), second.get())) {
        assertEquals(SearchResult.Status.Completed, result.status());
        assertTrue(result.optimal());
        assertTrue(result.complete());
        assertEquals(expected, new HashSet<>(result.matches()));
      }
    } finally {
      executor.shutdownNow();
    }
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
    MolGraph query = rawGraph("CC/C=C/CC"), target = rawGraph("CC/C=C\\CC");
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
    List<Set<Integer>> fragments = MatchGraph.fragments(atoms, bonds);
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

  @ParameterizedTest
  @CsvSource({"false,Connected", "true,Connected", "false,Fragments", "true,Fragments"})
  void connectedMCESRetainsOptimumInBothInputOrders(
      boolean reverse, MCS.Connectivity connectivity) throws Exception {
    MolGraph cage = rawGraph("C12C3C4C1C5C2C3C45"), bicyclic = rawGraph("C1CCC2CCCCC2C1");
    MolGraph query = reverse ? bicyclic : cage, target = reverse ? cage : bicyclic;
    MCS.Options options = new MCS.Options();
    options.objective = MCS.Objective.MCES;
    options.connectivity = connectivity;
    if (connectivity == MCS.Connectivity.Fragments)
      options.maxFragments = 1;
    options.timeoutMillis = 500;
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertEquals(SearchResult.Status.Completed, result.status());
    assertTrue(result.optimal());
    assertTrue(result.complete());
    assertEquals(8, result.matches().get(0).bondCount());
    assertEquals(8, result.matches().get(0).atomCount());
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  /** The best selections an exhaustive enumeration finds, see {@link #validatedOracle}. */
  private static final class Oracle {
    int bonds, atoms;
    final Set<String> keys = new TreeSet<>();
  }

  /**
   * Enumerates injective atom maps and every subset of their compatible bonds, and keeps the best
   * of those that {@link MCS#validateMapping} accepts: an oracle for stereochemistry, complete
   * rings and fragment limits, under which a mapping may leave compatible bonds out.
   */
  private static Oracle validatedOracle(MolGraph query, MolGraph target, MCS.Options options) {
    Oracle oracle = new Oracle();
    int[] mapping = new int[query.n];
    Arrays.fill(mapping, -1);
    enumerate(query, target, options, 0, mapping, new boolean[target.n], oracle);
    return oracle;
  }

  private static void enumerate(MolGraph query, MolGraph target, MCS.Options options, int queryAtom,
      int[] mapping, boolean[] used, Oracle oracle) {
    if (queryAtom < query.n) {
      enumerate(query, target, options, queryAtom + 1, mapping, used, oracle);
      for (int targetAtom = 0; targetAtom < target.n; targetAtom++) {
        if (used[targetAtom] || !MolGraph.ChemOps.atomsCompatFast(
                query, queryAtom, target, targetAtom, options.chemistry))
          continue;
        mapping[queryAtom] = targetAtom;
        used[targetAtom] = true;
        enumerate(query, target, options, queryAtom + 1, mapping, used, oracle);
        used[targetAtom] = false;
        mapping[queryAtom] = -1;
      }
      return;
    }
    List<int[]> compatible = new ArrayList<>();
    for (int first = 0; first < query.n; first++)
      for (int second : query.neighbors[first])
        if (first < second && mapping[first] >= 0 && mapping[second] >= 0
            && target.hasBond(mapping[first], mapping[second])
            && MolGraph.ChemOps.bondsCompatible(query, first, second, target, mapping[first],
                mapping[second], options.chemistry))
          compatible.add(new int[] {first, second});
    for (int subset = 1; subset < 1 << compatible.size(); subset++) {
      Map<Integer, Integer> atoms = new TreeMap<>();
      List<SearchResult.Bond> bonds = new ArrayList<>();
      for (int edge = 0; edge < compatible.size(); edge++)
        if ((subset & 1 << edge) != 0) {
          int[] bond = compatible.get(edge);
          atoms.put(bond[0], mapping[bond[0]]);
          atoms.put(bond[1], mapping[bond[1]]);
          bonds.add(new SearchResult.Bond(bond[0], bond[1]));
        }
      if (bonds.size() < oracle.bonds || bonds.size() == oracle.bonds && atoms.size() < oracle.atoms)
        continue;
      SearchResult.Match match = MCS.fragment(query, atoms, bonds);
      if (!MCS.validateMapping(query, target, match, options).isEmpty())
        continue;
      if (bonds.size() > oracle.bonds || atoms.size() > oracle.atoms) {
        oracle.bonds = bonds.size();
        oracle.atoms = atoms.size();
        oracle.keys.clear();
      }
      oracle.keys.add(key(match));
    }
  }

  /** The atoms and the selected bonds of {@code match}, in a form that compares as text. */
  private static String key(SearchResult.Match match) {
    List<String> bonds = new ArrayList<>();
    for (SearchResult.Bond bond : match.bondMapping().keySet()) bonds.add(bond.toString());
    Collections.sort(bonds);
    return match.atomMapping() + " " + bonds;
  }

  @Test
  void selectionsUnderStereoRingsAndFragmentLimitsMatchAnExhaustiveOracle() throws Exception {
    // Every best mapping, including those that leave a compatible bond out.
    String[][] pairs = {{"C[C@H](N)O", "C[C@@H](N)O"}, {"N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"},
        {"F/C=C/C", "F/C=C\\C"}, {"F/C=C/C=C/C", "F/C=C\\C=C/C"}, {"C1CC1C", "CCCC"},
        {"C1CCC1", "CC1CC1"}, {"C[C@@H]1CC1N", "C[C@H]1CC1N"}, {"OC1CC1C", "OCC1CC1"}};
    for (String[] pair : pairs) {
      MolGraph query = rawGraph(pair[0]), target = rawGraph(pair[1]);
      for (MCS.Connectivity connectivity : MCS.Connectivity.values())
        for (int flags = 1; flags < 4; flags++)
          for (boolean limits : new boolean[] {false, true}) {
            MCS.Options options = new MCS.Options();
            options.objective = MCS.Objective.MCES;
            options.connectivity = connectivity;
            options.timeoutMillis = 60_000;
            options.chemistry.useChirality = options.chemistry.useBondStereo = (flags & 1) != 0;
            options.chemistry.completeRingsOnly = (flags & 2) != 0;
            if (limits) {
              options.minFragmentAtoms = 3;
              options.maxFragments = 2;
            }
            String label = pair[0] + " " + pair[1] + " " + connectivity + " " + flags + " " + limits;
            Oracle expected = validatedOracle(query, target, options);
            SearchResult result = MCS.find(query).withOptions(options).matchAll(target, 100_000);
            assertTrue(result.optimal() && result.complete(), label);
            Set<String> actual = new TreeSet<>();
            for (SearchResult.Match match : result.matches()) actual.add(key(match));
            assertEquals(expected.bonds, (int) result.score(), label);
            assertEquals(expected.keys, actual, label);
            assertEquals(actual.size(), result.matches().size(), label);
            SearchResult single = MCS.find(query).withOptions(options).match(target);
            assertTrue(single.optimal(), label);
            if (single.found())
              assertTrue(expected.keys.contains(key(single.matches().get(0))), label);
          }
    }
  }
}

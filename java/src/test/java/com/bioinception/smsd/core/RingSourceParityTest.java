/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawGraph;
import static org.junit.jupiter.api.Assertions.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.params.provider.CsvSource;

class RingSourceParityTest {
  // An 80-layer hexagonal prism stack: a synthetic cage whose ring system exceeds the
  // ring-perception work limit, so it has no exact ring basis.
  private static final String STACK_SMILES = "C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21";
  private static MolGraph stack, relevantOnly;

  @BeforeAll
  static void graphs() throws Exception {
    stack = rawGraph(STACK_SMILES);
    relevantOnly = fourRingChain();
  }

  /** Four-member rings linked around one larger ring, shared with the native fixture. */
  private static MolGraph fourRingChain() {
    int ringCount = 12, atomCount = 1 + 3 * ringCount;
    List<List<Integer>> adjacency = new ArrayList<>();
    for (int atom = 0; atom < atomCount; atom++) adjacency.add(new ArrayList<Integer>());
    int previous = 0;
    for (int ring = 0; ring < ringCount; ring++) {
      int firstSide = 3 * ring + 1, secondSide = firstSide + 1, next = firstSide + 2;
      addBond(adjacency, previous, firstSide); addBond(adjacency, previous, secondSide);
      addBond(adjacency, firstSide, next); addBond(adjacency, secondSide, next);
      previous = next;
    }
    addBond(adjacency, previous, 0);
    int[][] neighbors = new int[atomCount][], orders = new int[atomCount][];
    int[] elements = new int[atomCount];
    Arrays.fill(elements, 6);
    for (int atom = 0; atom < atomCount; atom++) {
      neighbors[atom] = new int[adjacency.get(atom).size()];
      orders[atom] = new int[neighbors[atom].length];
      Arrays.fill(orders[atom], 1);
      for (int neighbor = 0; neighbor < neighbors[atom].length; neighbor++)
        neighbors[atom][neighbor] = adjacency.get(atom).get(neighbor);
    }
    return new MolGraph.Builder().atomCount(atomCount).atomicNumbers(elements)
        .neighbors(neighbors).bondOrders(orders).build();
  }

  private static void addBond(List<List<Integer>> adjacency, int first, int second) {
    adjacency.get(first).add(second);
    adjacency.get(second).add(first);
  }

  /** The existing topology beside independent atoms, for necessary-element count checks. */
  private static MolGraph withIsolatedAtoms(MolGraph graph, int[] elements) {
    int atomCount = graph.n + elements.length;
    int[] atomicNumbers = Arrays.copyOf(graph.atomicNum, atomCount);
    System.arraycopy(elements, 0, atomicNumbers, graph.n, elements.length);
    int[][] neighbors = new int[atomCount][], orders = new int[atomCount][];
    for (int atom = 0; atom < atomCount; atom++) {
      neighbors[atom] = atom < graph.n ? graph.neighbors[atom].clone() : new int[0];
      orders[atom] = new int[neighbors[atom].length];
      Arrays.fill(orders[atom], 1);
    }
    return new MolGraph.Builder().atomCount(atomCount).atomicNumbers(atomicNumbers)
        .neighbors(neighbors).bondOrders(orders).build();
  }

  private static void row(String id, String value) {
    System.out.println("RING_SOURCE_PARITY\t" + id + "\t" + value);
  }

  private static void status(String id, SearchResult result, boolean incomplete) {
    System.out.println("RING_SOURCE_ACTUAL\t" + id + "\t" + result.status()
        + "\t" + result.optimal() + "\t" + result.complete());
    System.out.println("RING_SOURCE_MATCH_ACTUAL\t" + id + "\t" + result.found()
        + "\t" + result.matches().size() + "\t" + result.nodesVisited());
    System.out.flush();
    assertEquals(incomplete ? SearchResult.Status.LimitReached : SearchResult.Status.Completed, result.status(), id);
    assertEquals(!incomplete, result.optimal(), id + "/optimal");
    assertEquals(!incomplete, result.complete(), id + "/complete");
    row(id, result.status().toString());
  }

  @Test
  void stoppedSearchKeepsItsStatusBeforeRingAdmission() {
    Substructure.Options options = new Substructure.Options();
    options.chemistry.completeRingsOnly = true;
    options.timeoutMillis = 0;
    SearchResult timedOut = Substructure.smarts("[#7;R2]").withOptions(options).match(stack);
    assertEquals(SearchResult.Status.TimedOut, timedOut.status());
    assertFalse(timedOut.optimal());
    assertFalse(timedOut.complete());
    options.timeoutMillis = 60_000;
    options.cancellation = () -> true;
    SearchResult cancelled = Substructure.smarts("[#7;R2]").withOptions(options).match(stack);
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
    assertFalse(cancelled.optimal());
    assertFalse(cancelled.complete());
  }

  @ParameterizedTest
  @CsvSource(value = {
      "[$([#7]-[#6]);R2]|true", "[$([$([#7])]);R2]|true",
      "[$([#7]),$([#7]-[#6]);R2]|true", "[$([#7]),$([#6]);R2]|false",
      "[!$([#7]);R2]|false", "[$([#6]-[#7]);R2]|false",
      "[$([#6]-[#6]);R2]|false", "[$([#7]);$([#6]);R2]|false"
  }, delimiter = '|')
  void recursiveAnchorProofKeepsBooleanAlternativesConservative(String pattern, boolean absent) {
    Substructure.Options options = new Substructure.Options();
    options.timeoutMillis = 60_000;
    options.chemistry.completeRingsOnly = true;
    SearchResult result = Substructure.smarts(pattern).withOptions(options).match(stack);
    status("stack/recursive-anchor/" + pattern, result, !absent);
    assertFalse(result.found());
    assertEquals(0, result.nodesVisited());
  }

  @Test
  void recursiveRequirementsCountRootAtomsRatherThanOverlappingPayloads() {
    Substructure.Options options = new Substructure.Options();
    options.timeoutMillis = 60_000;
    options.chemistry.completeRingsOnly = true;
    MolGraph oneNitrogen = withIsolatedAtoms(stack, new int[] {7});
    SearchResult cardinality = Substructure.smarts("[$([#7]-[#6]);R2].[$([#7]-[#6])]")
        .withOptions(options).match(oneNitrogen);
    status("stack/recursive-anchor/root-cardinality", cardinality, false);
    assertFalse(cardinality.found());
    assertEquals(0, cardinality.nodesVisited());
    SearchResult sharedAnchor = Substructure.smarts("[$([#7]-[#6]);$([#7]-[#6]);R2]")
        .withOptions(options).match(oneNitrogen);
    status("stack/recursive-anchor/shared-anchor", sharedAnchor, true);
    assertFalse(sharedAnchor.found());
    assertEquals(0, sharedAnchor.nodesVisited());
  }

  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  void ringSourceParityTable(boolean basisIncomplete) throws Exception {
    MolGraph graph = basisIncomplete ? stack : relevantOnly;
    String name = basisIncomplete ? "stack" : "relevant-only";
    assertEquals(basisIncomplete, graph.ringsTruncated());
    assertTrue(graph.relevantCyclesTruncated());
    int[][] orderedNeighbors = new int[graph.n][];
    for (int atom = 0; atom < graph.n; atom++) {
      orderedNeighbors[atom] = graph.neighbors[atom].clone();
      Arrays.sort(orderedNeighbors[atom]);
    }
    MolGraph countSource = new MolGraph.Builder().atomCount(graph.n)
        .atomicNumbers(graph.atomicNum.clone()).neighbors(orderedNeighbors).build();
    boolean countsIncomplete = countSource.ringsTruncated();
    assertEquals(countsIncomplete, graph.ringCountsIncomplete(), name + "/count-source");
    row(name + "/count-source", countsIncomplete ? "Incomplete" : "Complete");
    MolGraph carbon = rawGraph("C");
    for (boolean strict : new boolean[] {false, true}) {
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = 60_000;
      options.chemistry.completeRingsOnly = !strict;
      if (strict) options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      String id = name + (strict ? "/strict" : "/complete");
      boolean incomplete = strict ? countsIncomplete : graph.relevantCyclesTruncated();
      status(id + "/find", MCS.find(graph).withOptions(options).match(carbon), incomplete);
      // With one carbon beside it the consensus fragment stays writable: the cage alone needs
      // more simultaneously open ring closures than SMILES output allows.
      MCS.ConsensusResult consensus = MCS.consensus(java.util.Arrays.asList(graph, carbon), options, 1.0);
      assertEquals(incomplete ? SearchResult.Status.LimitReached : SearchResult.Status.Completed, consensus.status(), id);
      assertFalse(consensus.optimal(), id);
      assertEquals(!incomplete, consensus.complete(), id);
      row(id + "/consensus", consensus.status().toString());
    }
    Substructure.Options options = new Substructure.Options();
    options.timeoutMillis = 60_000;
    options.chemistry.completeRingsOnly = true;
    String[] patterns = {"[#7;R2]", "[#6;D99]", "[R2]", "[#6;R2]", "[#7,#6;R2]",
        "[!#7;R2]", "[$([#7]);R2]", "[#7,#7;R2]", "[#7;#6;R2]"};
    for (int item = 0; item < patterns.length; item++) {
      String pattern = patterns[item];
      boolean independentAbsence = item == 0 || item == 1 || item == 6 || item == 7;
      SearchResult result = Substructure.smarts(pattern).withOptions(options).match(graph);
      boolean needsRelevantProof = pattern.equals("[R2]") || pattern.equals("[#6;R2]")
          || pattern.equals("[#7,#6;R2]") || pattern.equals("[!#7;R2]");
      boolean incomplete = !independentAbsence
          && (basisIncomplete || (needsRelevantProof && graph.relevantCyclesTruncated()));
      status(name + "/smarts/" + pattern, result, incomplete);
      if (independentAbsence || item == 6 || item == 8) assertFalse(result.found());
      if (item == 0 || item == 6 || item == 7) assertEquals(0, result.nodesVisited());
    }
    MolGraph oneNitrogen = withIsolatedAtoms(graph, new int[] {7});
    SearchResult cardinality = Substructure.smarts("[#7;R2].[#7]").withOptions(options).match(oneNitrogen);
    status(name + "/element-cardinality", cardinality, false);
    assertFalse(cardinality.found());
    assertEquals(0, cardinality.nodesVisited());
    if (basisIncomplete) {
      int[] firstSixteen = new int[16];
      StringBuilder pattern = new StringBuilder();
      for (int element = 1; element <= 17; element++) {
        if (element <= 16) firstSixteen[element - 1] = element;
        if (pattern.length() > 0) pattern.append('.');
        pattern.append("[#").append(element).append(element == 17 ? ";R2]" : "]");
      }
      MolGraph firstSixteenPresent = withIsolatedAtoms(graph, firstSixteen);
      SearchResult seventeenth = Substructure.smarts(pattern.toString()).withOptions(options).match(firstSixteenPresent);
      status(name + "/seventeenth-element", seventeenth, false);
      assertFalse(seventeenth.found());
      assertEquals(0, seventeenth.nodesVisited());
    }
    for (int term = 0; term < 3; term++) {
      FingerprintOptions fingerprint = FingerprintOptions.ecfp();
      if (term == 0) fingerprint.ringSize = true;
      if (term == 1) fingerprint.cipBonds = true;
      if (term == 2) fingerprint = FingerprintOptions.withKeys(Collections.singletonList(
          FingerprintOptions.Key.special(FingerprintOptions.Key.Kind.AROMATIC_RINGS, 1)));
      final FingerprintOptions selected = fingerprint;
      if (basisIncomplete) {
        RingBasisUnavailableException error = assertThrows(RingBasisUnavailableException.class,
            () -> FingerprintEngine.fingerprint(graph, selected));
        assertEquals(RingBasisUnavailableException.MESSAGE, error.getMessage());
      } else {
        long[] words = assertDoesNotThrow(() -> FingerprintEngine.fingerprint(graph, selected));
        StringBuilder values = new StringBuilder();
        for (long word : words) values.append(Long.toUnsignedString(word)).append(',');
        row(name + "/fingerprint/" + term + "/words", values.toString());
      }
      row(name + "/fingerprint/" + term, basisIncomplete ? "Refused" : "Completed");
    }
    for (int minimum : new int[] {1, 2}) {
      FingerprintOptions keys = FingerprintOptions.withKeys(Collections.singletonList(
          new FingerprintOptions.Key("[R2]", minimum)));
      if (basisIncomplete) {
        RingBasisUnavailableException error = assertThrows(RingBasisUnavailableException.class,
            () -> FingerprintEngine.fingerprint(graph, keys));
        assertEquals(RingBasisUnavailableException.MESSAGE, error.getMessage());
      } else assertDoesNotThrow(() -> FingerprintEngine.fingerprint(graph, keys));
      row(name + "/key/" + minimum, basisIncomplete ? "Refused" : "Completed");
      if (basisIncomplete) {
        RingBasisUnavailableException error = assertThrows(RingBasisUnavailableException.class,
            () -> FingerprintEngine.fingerprintCounts(graph, keys));
        assertEquals(RingBasisUnavailableException.MESSAGE, error.getMessage());
      } else assertDoesNotThrow(() -> FingerprintEngine.fingerprintCounts(graph, keys));
      row(name + "/key-count/" + minimum, basisIncomplete ? "Refused" : "Completed");
    }
  }
}

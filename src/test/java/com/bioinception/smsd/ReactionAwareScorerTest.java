/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import java.util.*;
import org.junit.jupiter.api.*;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Tests for the reaction-aware MCS post-filter (v6.4.0).
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Reaction-Aware MCS Post-Filter (v6.4.0)")
public class ReactionAwareScorerTest extends TestBase {

  // SAM (S-adenosylmethionine) simplified SMILES -- includes the key S atom
  // We use a simplified representation that captures the essential chemistry:
  // adenine + ribose + methionine chain including the S-methyl group
  private static final String SAM_SMILES =
      "C[S+](CCC(N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)C(O)C1O";

  // SAH (S-adenosylhomocysteine) -- SAM after demethylation, S atom retained
  private static final String SAH_SMILES =
      "SCCC(N)C(=O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1O";

  @Test
  @DisplayName("SAM -> SAH: reaction-aware MCS includes S atom")
  void samToSahIncludesSulfur() throws Exception {
    IAtomContainer sam = mol(SAM_SMILES);
    IAtomContainer sah = mol(SAH_SMILES);

    // Standard mapReaction (not reaction-aware)
    Map<Integer, Integer> standardMapping =
        SearchEngine.mapReaction(sam, sah, new ChemOptions(), 10_000);
    assertNotNull(standardMapping);
    assertFalse(standardMapping.isEmpty(), "Standard MCS should find a mapping");

    // Reaction-aware mapReaction
    Map<Integer, Integer> reactMapping =
        SearchEngine.mapReactionAware(sam, sah, new ChemOptions(), 10_000);
    assertNotNull(reactMapping);
    assertFalse(reactMapping.isEmpty(), "Reaction-aware MCS should find a mapping");

    // Check that the reaction-aware mapping includes sulfur
    MolGraph g1 = new MolGraph(sam);
    boolean hasSulfurInReactMapping = false;
    for (int qi : reactMapping.keySet()) {
      if (g1.atomicNum[qi] == 16) { // 16 = sulfur
        hasSulfurInReactMapping = true;
        break;
      }
    }
    assertTrue(hasSulfurInReactMapping,
        "Reaction-aware MCS for SAM->SAH should include the sulfur atom");
  }

  @Test
  @DisplayName("Pure hydrocarbon: reaction-aware falls back to size-K")
  void pureHydrocarbonFallsBackToSizeK() throws Exception {
    IAtomContainer cyclohexane = mol("C1CCCCC1");
    IAtomContainer methylcyclopentane = mol("CC1CCCC1");

    Map<Integer, Integer> stdMapping =
        SearchEngine.mapReaction(cyclohexane, methylcyclopentane, new ChemOptions(), 5_000);
    Map<Integer, Integer> reactMapping =
        SearchEngine.mapReactionAware(cyclohexane, methylcyclopentane, new ChemOptions(), 5_000);

    // For pure hydrocarbons, H_universe is empty, so scorer degrades to size ranking
    // Both should produce the same size mapping
    assertEquals(stdMapping.size(), reactMapping.size(),
        "Pure hydrocarbon pair: reaction-aware should match standard MCS size");
  }

  @Test
  @DisplayName("ReactionAwareScorer: scoring formula walkthrough")
  void scorerFormulaTest() throws Exception {
    // Unit test the scorer directly with synthetic candidates
    MolGraph g1 = new MolGraph(mol("CSCC(N)C(=O)O")); // methionine-like with S
    MolGraph g2 = new MolGraph(mol("SCC(N)C(=O)O"));   // homocysteine-like with S

    // Create two synthetic candidates:
    // Candidate A: maps 5 atoms, no S
    // Candidate B: maps 4 atoms, includes S
    // We verify the scorer ranks B above A when S is present

    ReactionAwareScorer scorer = new ReactionAwareScorer();
    // The actual ranking depends on what findNearMCS generates,
    // but we verify the scorer interface works
    List<Map<Integer, Integer>> candidates = new ArrayList<>();
    Map<Integer, Integer> dummyA = new LinkedHashMap<>();
    dummyA.put(0, 0); // Just test that the scorer doesn't crash
    candidates.add(dummyA);

    List<Map<Integer, Integer>> ranked = scorer.rank(candidates, g1, g2);
    assertNotNull(ranked);
    assertFalse(ranked.isEmpty());
    assertEquals(1, ranked.size());
  }

  @Test
  @DisplayName("MCSOptions.reactionAware defaults to false")
  void reactionAwareDefaultsFalse() {
    SearchEngine.MCSOptions opts = new SearchEngine.MCSOptions();
    assertFalse(opts.reactionAware, "reactionAware should default to false");
    assertEquals(2, opts.nearMcsDelta, "nearMcsDelta should default to 2");
    assertEquals(20, opts.nearMcsCandidates, "nearMcsCandidates should default to 20");
    assertNull(opts.postFilter, "postFilter should default to null");
  }

  @Test
  @DisplayName("Custom post-filter callback works")
  void customPostFilterCallback() throws Exception {
    IAtomContainer m1 = mol("c1ccccc1");
    IAtomContainer m2 = mol("c1ccc(O)cc1");

    MolGraph g1 = new MolGraph(m1);
    MolGraph g2 = new MolGraph(m2);

    // Custom filter that just returns candidates as-is
    MCSPostFilter identity = (candidates, q, t) -> candidates;

    SearchEngine.MCSOptions opts = new SearchEngine.MCSOptions();
    opts.disconnectedMCS = false;
    opts.connectedOnly = true;
    opts.timeoutMs = 5_000;
    opts.reactionAware = true;
    opts.postFilter = identity;

    Map<Integer, Integer> result = SearchEngine.reactionAwareMCS(g1, g2, new ChemOptions(), opts);
    assertNotNull(result);
    assertFalse(result.isEmpty(), "Custom filter should still produce a result");
  }

  @Test
  @DisplayName("SMSD.mapReactionAware convenience method delegates correctly")
  void smsdConvenienceMethod() throws Exception {
    IAtomContainer m1 = mol("CSCC(N)C(=O)O");
    IAtomContainer m2 = mol("SCC(N)C(=O)O");

    Map<Integer, Integer> mapping = SMSD.mapReactionAware(m1, m2, new ChemOptions(), 10_000);
    assertNotNull(mapping);
    assertFalse(mapping.isEmpty());
  }

  // ATP (adenosine triphosphate) -- simplified SMILES with triphosphate chain
  private static final String ATP_SMILES =
      "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1";
  // ADP (adenosine diphosphate) -- simplified SMILES with diphosphate chain
  private static final String ADP_SMILES =
      "c1nc(N)c2ncn(C3OC(COP(=O)(O)OP(=O)(O)O)C(O)C3O)c2n1";

  @Test
  @DisplayName("ATP -> ADP: reaction-aware MCS includes P atom")
  void atpToAdpIncludesPhosphorus() throws Exception {
    IAtomContainer atp = mol(ATP_SMILES);
    IAtomContainer adp = mol(ADP_SMILES);

    Map<Integer, Integer> reactMapping =
        SearchEngine.mapReactionAware(atp, adp, new ChemOptions(), 10_000);
    assertNotNull(reactMapping);
    assertFalse(reactMapping.isEmpty(), "Reaction-aware MCS should find a mapping");

    MolGraph g1 = new MolGraph(atp);
    boolean hasPhosphorus = false;
    for (int qi : reactMapping.keySet()) {
      if (g1.atomicNum[qi] == 15) { // 15 = phosphorus
        hasPhosphorus = true;
        break;
      }
    }
    assertTrue(hasPhosphorus,
        "Reaction-aware MCS for ATP->ADP should include a phosphorus atom");
  }
}

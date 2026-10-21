/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.core.MCS;

import com.bioinception.smsd.TestSupport;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import java.util.Arrays;
import java.util.Map;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.interfaces.IAtomContainer;

class SearchChemistryRegressionTest {
  private static ChemOptions strictFusion() {
    ChemOptions chemistry = new ChemOptions();
    chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
    return chemistry;
  }

  @Test
  void strictFusionBoundAllowsRingToChainMatches() throws Exception {
    MolGraph query = new MolGraph(rawMolecule("C1CCCCC1"));
    MolGraph target = new MolGraph(rawMolecule("CCCCCC"));
    ChemOptions chemistry = strictFusion();
    assertTrue(MCS.similarityUpperBound(query, target, chemistry) >= 5.0 / 7.0);
    chemistry.ringMatchesRingOnly = true;
    assertEquals(0.0, MCS.similarityUpperBound(query, target, chemistry));
  }

  @Test
  void strictFusionScreenRetainsAValidRingToChainPair() throws Exception {
    IAtomContainer query = rawMolecule("C1CCCCC1");
    IAtomContainer target = rawMolecule("CCCCCC");
    Map<Integer, Map<Integer, Integer>> mappings = TestSupport.pairMappings(
        Arrays.asList(query, target), strictFusion(), new MCS.Options(), 0.5, 1);
    assertEquals(1, mappings.size());
    assertEquals(5, mappings.get(1).size());
  }

  @ParameterizedTest
  @CsvSource({
      "OC1CCCCC1O, NCCCCCCN, 5",
      "C1CCC(CC1)CCC1CCCCC1, CCCCCCCCCCCC, 12"
  })
  void strictFusionSearchKeepsValidRingToChainMappings(
      String querySmiles, String targetSmiles, int expectedSize) throws Exception {
    MolGraph query = new MolGraph(rawMolecule(querySmiles));
    MolGraph target = new MolGraph(rawMolecule(targetSmiles));
    ChemOptions chemistry = strictFusion();
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 2_000;
    Map<Integer, Integer> mapping = TestSupport.mcsMapping(query, target, chemistry, options);
    assertEquals(expectedSize, mapping.size());
    assertTrue(MCS.validateMapping(query, target, mapping, chemistry).isEmpty());
  }

  @ParameterizedTest
  @CsvSource({
      "C1CCCCC1, C1CCC2CCCCC2C1, false, 6",
      "CCC, C1CC1, true, 2",
      "CCC, C1CC1, false, 3"
  })
  void bondGrowthPreservesMappedBondsAndTheRequestedInducedPolicy(
      String querySmiles, String targetSmiles, boolean induced, int expectedSize) throws Exception {
    MolGraph query = new MolGraph(rawMolecule(querySmiles));
    MolGraph target = new MolGraph(rawMolecule(targetSmiles));
    query.ensureCanonical();
    target.ensureCanonical();
    ChemOptions chemistry = new ChemOptions();
    MCS.Options options=new MCS.Options();options.chemistry=chemistry;
    options.objective=induced?MCS.Objective.MCIS:MCS.Objective.MCES;
    SearchResult result=MCS.find(query).withOptions(options).match(target);
    Map<Integer,Integer> mapping=result.matches().get(0).atomMapping();
    assertTrue(result.optimal());
    assertEquals(expectedSize, mapping.size());
    if (induced) {
      for (int queryAtom : mapping.keySet())
        for (int otherQueryAtom : mapping.keySet())
          assertEquals(query.hasBond(queryAtom, otherQueryAtom),
              target.hasBond(mapping.get(queryAtom), mapping.get(otherQueryAtom)));
    }
  }

  private static final String[] PREPARED_SEARCH_MOLECULES = {
      "C[C@H](N)C(=O)O", "[13CH3]C(=O)[O-].[Na+]", "c1ccc2c(c1)[nH]c1ccccc12",
      "C/C=C/C(=O)Nc1ccc(Cl)cc1", "C1CC2CCC1C2", "O=C1NC(=O)C(=O)N1", "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
      "CC(C)(C)OC(=O)N1CCC(CC1)c1nc2ccccc2s1"};

  private static ChemOptions[] preparedSearchChemistries() {
    ChemOptions any = new ChemOptions();
    any.matchBondOrder = ChemOptions.BondOrderMode.ANY;
    ChemOptions strict = new ChemOptions();
    strict.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
    strict.matchFormalCharge = strict.matchIsotope = strict.ringMatchesRingOnly = true;
    ChemOptions stereo = new ChemOptions();
    stereo.useChirality = stereo.useBondStereo = true;
    ChemOptions fusion = strictFusion();
    fusion.useThreeHopNLF = true;
    ChemOptions untyped = new ChemOptions();
    untyped.matchAtomType = false;
    ChemOptions tautomer = new ChemOptions();
    tautomer.tautomerAware = true;
    return new ChemOptions[] {new ChemOptions(), any, strict, stereo, fusion, untyped, tautomer};
  }

  @Test
  void sharedSearchRowsEqualPerAtomCompatibility() throws Exception {
    for (ChemOptions chemistry : preparedSearchChemistries())
      for (String querySmiles : PREPARED_SEARCH_MOLECULES)
        for (String targetSmiles : PREPARED_SEARCH_MOLECULES) {
          MolGraph query = new MolGraph(rawMolecule(querySmiles)), target = new MolGraph(rawMolecule(targetSmiles));
          Substructure.SubgraphSearch search = new Substructure.SubgraphSearch(query,
              new Substructure.QueryTables(query, chemistry), target, chemistry, new SearchResult.Control(60_000, null));
          for (int queryAtom = 0; queryAtom < query.atomCount(); queryAtom++) {
            long[] compatible = new long[search.tWords], rootCompatible = new long[search.tWords];
            for (int targetAtom = 0; targetAtom < target.atomCount(); targetAtom++) {
              if (!MolGraph.ChemOps.atomsCompatFast(query, queryAtom, target, targetAtom, chemistry)) continue;
              compatible[targetAtom >>> 6] |= 1L << (targetAtom & 63);
              if (search.rootFeasible(queryAtom, targetAtom)) rootCompatible[targetAtom >>> 6] |= 1L << (targetAtom & 63);
            }
            String context = querySmiles + " in " + targetSmiles + " atom " + queryAtom;
            assertArrayEquals(compatible, search.domain[queryAtom], context);
            assertArrayEquals(rootCompatible, search.rootRow(queryAtom), context);
          }
        }
  }

  @Test
  void neighbourhoodLabelsCountEveryAtomWithinTheRadius() throws Exception {
    for (String smiles : PREPARED_SEARCH_MOLECULES) {
      MolGraph graph = new MolGraph(rawMolecule(smiles));
      for (int atom = 0; atom < graph.atomCount(); atom++)
        for (int radius = 1; radius <= 3; radius++) {
          int[] distance = new int[graph.atomCount()];
          Arrays.fill(distance, -1);
          distance[atom] = 0;
          java.util.ArrayDeque<Integer> pending = new java.util.ArrayDeque<>(Arrays.asList(atom));
          java.util.TreeMap<Integer, Integer> counts = new java.util.TreeMap<>();
          while (!pending.isEmpty()) {
            int current = pending.remove();
            if (distance[current] == radius) continue;
            for (int neighbor : graph.neighbors[current])
              if (distance[neighbor] < 0) {
                distance[neighbor] = distance[current] + 1;
                counts.merge(MolGraph.nlfLabel(graph, neighbor), 1, Integer::sum);
                pending.add(neighbor);
              }
          }
          int[] expected = new int[counts.size() * 2];
          int index = 0;
          for (Map.Entry<Integer, Integer> entry : counts.entrySet()) {
            expected[index++] = entry.getKey();
            expected[index++] = entry.getValue();
          }
          int[] actual = radius == 1 ? MolGraph.buildNLF1(graph, atom)
              : radius == 2 ? MolGraph.buildNLF2(graph, atom) : MolGraph.buildNLF3(graph, atom);
          assertArrayEquals(expected, actual, smiles + " atom " + atom + " radius " + radius);
        }
    }
  }
}

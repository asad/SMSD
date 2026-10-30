/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

class MCSForestBoundTest {
  private static MolGraph molecule(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }

  @ParameterizedTest
  @CsvSource({
      "CCC(C)CC, CCCCCC, default, 5",
      "CC(C)(C)OCC, CC(C)(C)NCC, strict, 4",
      "CC=CC, CCCC, compat-fmcs, 4",
      "CC=CC, CCCC, default, 2",
      "CN[13CH2]O, CN[12CH2]S, strict, 3",
      "C.CC, CCC, default, 2",
      "N, O, default, 0"
  })
  void branchBoundsUseTheEnabledPredicates(String first, String second, String profile, int atoms) throws Exception {
    MolGraph query = molecule(first), target = molecule(second);
    MCS.Options options = new MCS.Options();
    options.chemistry = ChemOptions.profile(profile);
    options.chemistry.matchIsotope = true;
    options.timeoutMillis = 1000;
    SearchResult.Control control = new SearchResult.Control(1000, null);
    BridgeBound bound = BridgeBound.create(query, target, options.chemistry, control);
    assertNotNull(bound);
    assertEquals(atoms, bound.upper());
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal() && result.complete());
    assertEquals(atoms, result.found() ? result.matches().get(0).atomCount() : 0);
    if (result.found())
      result.matches().get(0).atomMapping().forEach((queryAtom, targetAtom) ->
          assertTrue(bound.bound(queryAtom, targetAtom) >= atoms));
  }

  @Test
  void unsupportedOrStoppedConstructionLeavesNoTable() throws Exception {
    MolGraph chain = molecule("CCCCCC");
    SearchResult.Control control = new SearchResult.Control(1000, null);
    assertNull(BridgeBound.create(molecule("C1CCCCC1"), chain, new ChemOptions(), control));
    assertNull(BridgeBound.create(molecule(repeat('C', 513)), chain, new ChemOptions(), control));
    SearchResult.Control cancelled = new SearchResult.Control(1000, () -> true);
    assertNull(BridgeBound.create(chain, chain, new ChemOptions(), cancelled));
    assertEquals(SearchResult.Status.Cancelled, cancelled.status());
  }

  @Test
  void largeBranchedChainsKeepTheMaximumPath() throws Exception {
    MolGraph query = molecule(repeat('C', 64));
    MolGraph target = molecule(repeat('C', 31) + "C(C)" + repeat('C', 31));
    for (boolean reverse : new boolean[] {false, true}) {
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = 3000;
      options.nodeLimit = 1000000;
      MolGraph first = reverse ? target : query, second = reverse ? query : target;
      SearchResult result = MCS.find(first).withOptions(options).match(second);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal() && result.complete());
      assertEquals(63, result.matches().get(0).atomCount());
      assertEquals(62, result.matches().get(0).bondCount());
      assertTrue(MCS.validateMapping(first, second, result.matches().get(0), options).isEmpty());
    }
  }

  private static String repeat(char value, int count) {
    char[] result = new char[count];
    java.util.Arrays.fill(result, value);
    return new String(result);
  }
}

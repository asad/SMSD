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

class MCSBridgeBoundTest {
  private static MolGraph molecule(String smiles) throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }

  @Test
  void cyclicRootsKeepTheirWholeParts() throws Exception {
    MolGraph molecule = molecule("CC1CCCCC1");
    SearchResult.Control control = new SearchResult.Control(1000, null);
    BridgeBound bound = BridgeBound.create(molecule, molecule, new ChemOptions(), control);
    assertNotNull(bound);
    assertFalse(bound.acyclic());
    assertEquals(7, bound.upper());
    assertEquals(7, bound.bound(1, 1));
    assertNull(BridgeBound.create(molecule("C1CCCCC1"), molecule, new ChemOptions(), control));
  }

  @ParameterizedTest
  @CsvSource({
      "CC1CC1, CC1CCC1, default, true",
      "CC=C1CCCCC1, CCC1CCCCC1, compat-fmcs, true",
      "CC1CCC2CCCCC2C1, CC1CCC2CCCC2C1, strict, true",
      "CC1CCC2(CC1)CCCC2, CC1CCC2(CC1)CCCCC2, default, true",
      "OC1CCCCC1, NC1CCCCC1, default, false",
      "CN[13CH]1CCCCC1, CN[12CH]1CCCCC1, default, true"
  })
  void relaxedPredicatesRemainUpperBounds(String first, String second, String profile, boolean typed) throws Exception {
    MolGraph query = molecule(first), target = molecule(second);
    MCS.Options options = new MCS.Options();
    options.timeoutMillis = 1000;
    options.chemistry = ChemOptions.profile(profile);
    options.chemistry.matchAtomType = typed;
    options.chemistry.matchIsotope = true;
    BridgeBound bound = BridgeBound.create(query, target, options.chemistry, new SearchResult.Control(1000, null));
    assertNotNull(bound);
    assertFalse(bound.acyclic());
    SearchResult result = MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal() && result.complete() && result.found());
    int atoms = result.matches().get(0).atomCount();
    assertTrue(bound.upper() >= atoms);
    result.matches().get(0).atomMapping().forEach((queryAtom, targetAtom) ->
        assertTrue(bound.bound(queryAtom, targetAtom) >= atoms));
    assertTrue(MCS.validateMapping(query, target, result.matches().get(0), options).isEmpty());
  }

  @Test
  void aLargeCyclicMatchCanHaveOneBondPerAtom() throws Exception {
    MolGraph query = molecule(repeat('C', 58) + "C1CCCCC1");
    MolGraph target = molecule(repeat('C', 29) + "C(C)" + repeat('C', 27) + "C1CCCCC1");
    for (boolean reverse : new boolean[] {false, true}) {
      MCS.Options options = new MCS.Options();
      options.timeoutMillis = 3000;
      options.nodeLimit = 1000000;
      MolGraph first = reverse ? target : query, second = reverse ? query : target;
      SearchResult result = MCS.find(first).withOptions(options).match(second);
      assertEquals(SearchResult.Status.Completed, result.status());
      assertTrue(result.optimal() && result.complete());
      assertEquals(63, result.matches().get(0).atomCount());
      assertEquals(63, result.matches().get(0).bondCount());
      assertTrue(MCS.validateMapping(first, second, result.matches().get(0), options).isEmpty());
    }
  }

  private static String repeat(char value, int count) {
    char[] result = new char[count];
    java.util.Arrays.fill(result, value);
    return new String(result);
  }
}

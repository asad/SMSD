/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.InputStream;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

@DisplayName("Aromaticity parity")
public class AromaticityParityTest extends TestBase {

  private static final ObjectMapper MAPPER = new ObjectMapper();

  record AromaticityCase(
      String name,
      String kekule,
      int aromatic_atom_count,
      int aromatic_bond_count,
      int ring_atom_count) {}

  static Stream<Arguments> parityCorpus() throws Exception {
    return loadCorpus().stream().map(entry -> Arguments.arguments(entry.name(), entry));
  }

  @ParameterizedTest(name = "{0}")
  @MethodSource("parityCorpus")
  void cdkDaylightParityCorpus(String ignoredName, AromaticityCase entry) throws Exception {
    IAtomContainer mol = mol(entry.kekule());

    assertEquals(entry.aromatic_atom_count(), aromaticAtomCount(mol));
    assertEquals(entry.aromatic_bond_count(), aromaticBondCount(mol));
    assertEquals(entry.ring_atom_count(), ringAtomCount(mol));
  }

  private static List<AromaticityCase> loadCorpus() throws Exception {
    try (InputStream in =
        AromaticityParityTest.class
            .getClassLoader()
            .getResourceAsStream("com/bioinception/smsd/aromaticity_parity.json")) {
      assertNotNull(in, "Missing aromaticity parity corpus resource");
      return MAPPER.readValue(in, new TypeReference<List<AromaticityCase>>() {});
    }
  }

  private static int aromaticAtomCount(IAtomContainer mol) {
    int count = 0;
    for (IAtom atom : mol.atoms()) {
      if (atom.isAromatic()) count++;
    }
    return count;
  }

  private static int aromaticBondCount(IAtomContainer mol) {
    int count = 0;
    for (IBond bond : mol.bonds()) {
      if (bond.isAromatic()) count++;
    }
    return count;
  }

  private static int ringAtomCount(IAtomContainer mol) {
    int count = 0;
    for (IAtom atom : mol.atoms()) {
      if (atom.isInRing()) count++;
    }
    return count;
  }
}

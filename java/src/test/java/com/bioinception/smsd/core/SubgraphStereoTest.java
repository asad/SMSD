/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.InputStream;
import java.util.*;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Stereo kept by the fragments of a selection and by R-group decomposition, from the fixture also
 * read by the C++ and Python tests. Java fragments replace lost neighbours by hydrogen, so they are
 * compared with the capped SMILES of each case, which carries the same stereo.
 */
class SubgraphStereoTest {
  private static JsonNode fixture;

  @BeforeAll
  static void load() throws Exception {
    try (InputStream in = SubgraphStereoTest.class.getClassLoader()
        .getResourceAsStream("com/bioinception/smsd/subgraph_stereo.json")) {
      assertNotNull(in, "Missing shared subgraph stereo resource");
      fixture = new ObjectMapper().readTree(in);
    }
  }

  /** SMSD's own canonical SMILES, which writes tetrahedral and double-bond marks. */
  private static String absolute(IAtomContainer molecule) {
    return new MolGraph(molecule).toCanonicalSmiles();
  }

  private static boolean hasStereo(IAtomContainer molecule) {
    return molecule.stereoElements().iterator().hasNext();
  }

  @Test
  void fragmentsKeepOnlyTheStereoTheirAtomsDefine() throws Exception {
    for (JsonNode item : fixture.get("subgraphs")) {
      String name = item.get("name").asText();
      MolGraph graph = MolGraph.copyOf(rawMolecule(item.get("smiles").asText()));
      Map<Integer, Integer> atoms = new TreeMap<>();
      for (JsonNode atom : item.get("atoms")) atoms.put(atom.asInt(), atom.asInt());
      IAtomContainer fragment = MatchGraph.toMolecule(graph, atoms, MatchGraph.allBonds(graph, atoms).keySet());
      assertEquals(absolute(rawMolecule(item.get("capped").asText())), absolute(fragment), name);
      assertEquals(item.get("stereo").asBoolean(), hasStereo(fragment), name);
    }
  }

  @Test
  void rGroupsKeepOnlyTheStereoTheirAtomsDefine() throws Exception {
    for (JsonNode item : fixture.get("rgroups")) {
      String name = item.get("name").asText();
      Substructure.Decomposition decomposition = Substructure.decomposeRGroups(
          rawMolecule(item.get("core").asText()),
          Collections.singletonList(rawMolecule(item.get("smiles").asText())), new Substructure.Options()).get(0);
      assertTrue(decomposition.search().found(), name);
      List<String> groups = new ArrayList<>(), expected = new ArrayList<>();
      boolean stereo = false;
      for (IAtomContainer group : decomposition.rgroups().values()) {
        groups.add(absolute(group));
        stereo |= hasStereo(group);
      }
      for (JsonNode smiles : item.get("capped")) expected.add(absolute(rawMolecule(smiles.asText())));
      Collections.sort(groups);
      Collections.sort(expected);
      assertEquals(expected, groups, name);
      assertEquals(item.get("stereo").asBoolean(), stereo, name);
    }
  }
}

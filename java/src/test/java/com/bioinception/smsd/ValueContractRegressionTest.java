/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd;

import static com.bioinception.smsd.TestSupport.list;
import static com.bioinception.smsd.TestSupport.mapping;
import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.cli.SMSDcli.MolIO.Query;
import com.bioinception.smsd.core.SearchEngine.MCSProfiledResult;
import com.bioinception.smsd.core.SearchEngine.MCSResult;
import com.bioinception.smsd.core.SearchEngine.MCSStageTimers;
import com.bioinception.smsd.core.SearchEngine.SubstructureResult;
import com.bioinception.smsd.core.SearchEngine.SubstructureStats;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.node.ObjectNode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/** Stable constructor, accessor, JSON and value contracts across supported Java versions. */
class ValueContractRegressionTest {
  private static final ObjectMapper JSON = new ObjectMapper();
  private static final String STATS_JSON = "{\"nodesVisited\":1,\"backtracks\":2,"
      + "\"candidatesTried\":3,\"prunesAtom\":4,\"prunesBond\":5,\"prunesDegree\":6,"
      + "\"prunesNLF\":7,\"timeMillis\":8,\"timeout\":true,\"solutions\":9}";
  private static final String STATS_TEXT = "SubstructureStats[nodesVisited=1, backtracks=2, "
      + "candidatesTried=3, prunesAtom=4, prunesBond=5, prunesDegree=6, prunesNLF=7, "
      + "timeMillis=8, timeout=true, solutions=9]";
  private static final String RESULT_JSON = "{\"mapping\":{\"0\":2,\"1\":3},\"size\":2,"
      + "\"overlapCoefficient\":0.5,\"mcsSmiles\":\"CC\"}";
  private static final String RESULT_TEXT = "MCSResult[mapping={0=2, 1=3}, size=2, "
      + "overlapCoefficient=0.5, mcsSmiles=CC]";
  private static final String TIMERS_JSON = "{\"orientationUs\":1,\"seedsUs\":2,\"mcSplitUs\":3,"
      + "\"bkUs\":4,\"mcGregorUs\":5,\"repairUs\":6,\"totalUs\":7,\"bestAfterGreedy\":8,"
      + "\"bestAfterSeed\":9,\"bestAfterBK\":10,\"bestAfterMcGregor\":11}";
  private static final String TIMERS_TEXT = "MCSStageTimers[orientationUs=1, seedsUs=2, "
      + "mcSplitUs=3, bkUs=4, mcGregorUs=5, repairUs=6, totalUs=7, bestAfterGreedy=8, "
      + "bestAfterSeed=9, bestAfterBK=10, bestAfterMcGregor=11]";

  private static SubstructureStats stats() {
    return new SubstructureStats(1, 2, 3, 4, 5, 6, 7, 8, true, 9);
  }

  private static MCSStageTimers timers() {
    return new MCSStageTimers(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11);
  }

  static Stream<Arguments> values() {
    Map<Integer, Integer> atoms = mapping(0, 2, 1, 3);
    MCSResult result = new MCSResult(atoms, 2, 0.5, "CC");
    return Stream.of(
        Arguments.of(stats(), STATS_JSON, 820997278, STATS_TEXT),
        Arguments.of(new SubstructureResult(true, list(atoms), stats()),
            "{\"exists\":true,\"mappings\":[{\"0\":2,\"1\":3}],\"stats\":" + STATS_JSON + "}",
            822181354, "SubstructureResult[exists=true, mappings=[{0=2, 1=3}], stats="
                + STATS_TEXT + "]"),
        Arguments.of(result, RESULT_JSON, -1138630306, RESULT_TEXT),
        Arguments.of(timers(), TIMERS_JSON, -320062458, TIMERS_TEXT),
        Arguments.of(new MCSProfiledResult(result, timers()),
            "{\"result\":" + RESULT_JSON + ",\"timers\":" + TIMERS_JSON + "}", -1257863576,
            "MCSProfiledResult[result=" + RESULT_TEXT + ", timers=" + TIMERS_TEXT + "]"),
        Arguments.of(new Query(true, "C*N", null),
            "{\"isSmarts\":true,\"text\":\"C*N\",\"container\":null}", 3221768,
            "Query[isSmarts=true, text=C*N, container=null]"),
        Arguments.of(new SubstructureResult(false, null, null),
            "{\"exists\":false,\"mappings\":null,\"stats\":null}", 1188757,
            "SubstructureResult[exists=false, mappings=null, stats=null]"),
        Arguments.of(new MCSResult(null, 0, Double.NaN, null),
            "{\"mapping\":null,\"size\":0,\"overlapCoefficient\":\"NaN\",\"mcsSmiles\":null}",
            2131230720, "MCSResult[mapping=null, size=0, overlapCoefficient=NaN, mcsSmiles=null]"),
        Arguments.of(new MCSProfiledResult(null, null), "{\"result\":null,\"timers\":null}", 0,
            "MCSProfiledResult[result=null, timers=null]"),
        Arguments.of(new Query(false, null, null),
            "{\"isSmarts\":false,\"text\":null,\"container\":null}", 1188757,
            "Query[isSmarts=false, text=null, container=null]"),
        Arguments.of(new MCSResult(atoms, 2, -0.0, "CC"),
            "{\"mapping\":{\"0\":2,\"1\":3},\"size\":2,\"overlapCoefficient\":-0.0,\"mcsSmiles\":\"CC\"}",
            -2147360418, "MCSResult[mapping={0=2, 1=3}, size=2, overlapCoefficient=-0.0, mcsSmiles=CC]"),
        Arguments.of(new MCSResult(atoms, 2, 0.0, "CC"),
            "{\"mapping\":{\"0\":2,\"1\":3},\"size\":2,\"overlapCoefficient\":0.0,\"mcsSmiles\":\"CC\"}",
            123230, "MCSResult[mapping={0=2, 1=3}, size=2, overlapCoefficient=0.0, mcsSmiles=CC]"));
  }

  @ParameterizedTest(name = "{0}")
  @MethodSource("values")
  void valueAndJsonContracts(Object value, String expectedJson, int expectedHash,
                             String expectedText) throws Exception {
    JsonNode expected = JSON.readTree(expectedJson);
    assertEquals(expected, JSON.readTree(JSON.writeValueAsString(value)));
    Object restored = JSON.readValue(expectedJson, value.getClass());
    assertEquals(value, restored);
    assertEquals(restored, value);
    assertEquals(expectedHash, value.hashCode());
    assertEquals(expectedHash, restored.hashCode());
    assertEquals(expectedText, value.toString());
    assertEquals(expectedText, restored.toString());
    assertNotEquals(value, null);
    assertNotEquals(value, expectedText);
  }

  @Test
  void canonicalConstructorsRetainComponentAccessors() {
    SubstructureStats statistics = stats();
    assertArrayEquals(new long[] {1, 2, 3, 4, 5, 6, 7, 8}, new long[] {
        statistics.nodesVisited(), statistics.backtracks(), statistics.candidatesTried(),
        statistics.prunesAtom(), statistics.prunesBond(), statistics.prunesDegree(),
        statistics.prunesNLF(), statistics.timeMillis()});
    assertTrue(statistics.timeout());
    assertEquals(9, statistics.solutions());
    MCSStageTimers stages = timers();
    assertArrayEquals(new long[] {1, 2, 3, 4, 5, 6, 7}, new long[] {
        stages.orientationUs(), stages.seedsUs(), stages.mcSplitUs(), stages.bkUs(),
        stages.mcGregorUs(), stages.repairUs(), stages.totalUs()});
    assertArrayEquals(new int[] {8, 9, 10, 11}, new int[] {stages.bestAfterGreedy(),
        stages.bestAfterSeed(), stages.bestAfterBK(), stages.bestAfterMcGregor()});
    Map<Integer, Integer> atoms = mapping(0, 2, 1, 3);
    MCSResult result = new MCSResult(atoms, 2, 0.5, "CC");
    assertSame(atoms, result.mapping());
    assertEquals(2, result.size());
    assertEquals(0.5, result.overlapCoefficient());
    assertEquals("CC", result.mcsSmiles());
    List<Map<Integer, Integer>> maps = list(atoms);
    SubstructureResult substructure = new SubstructureResult(true, maps, statistics);
    assertTrue(substructure.exists());
    assertSame(maps, substructure.mappings());
    assertSame(statistics, substructure.stats());
    MCSProfiledResult profile = new MCSProfiledResult(result, stages);
    assertSame(result, profile.result());
    assertSame(stages, profile.timers());
    Query query = new Query(true, "C*N", null);
    assertTrue(query.isSmarts());
    assertEquals("C*N", query.text());
    assertNull(query.container());
  }

  @Test
  void constructorCollectionsRetainReferenceSemantics() {
    Map<Integer, Integer> atoms = new HashMap<>();
    List<Map<Integer, Integer>> maps = new ArrayList<>();
    MCSResult result = new MCSResult(atoms, 0, 0.0, "");
    SubstructureResult substructure = new SubstructureResult(false, maps, null);
    atoms.put(1, 2);
    maps.add(atoms);
    assertSame(atoms, result.mapping());
    assertEquals(mapping(1, 2), result.mapping());
    assertSame(maps, substructure.mappings());
    assertSame(atoms, substructure.mappings().get(0));
  }

  @Test
  void signedZeroIsDistinctAndNanPayloadsCompareEqual() {
    MCSResult positive = new MCSResult(null, 0, 0.0, null);
    MCSResult negative = new MCSResult(null, 0, -0.0, null);
    assertNotEquals(positive, negative);
    assertNotEquals(positive.hashCode(), negative.hashCode());
    MCSResult canonical = new MCSResult(null, 0, Double.NaN, null);
    MCSResult payload = new MCSResult(null, 0, Double.longBitsToDouble(0x7ff8000000000001L), null);
    assertEquals(canonical, payload);
    assertEquals(canonical.hashCode(), payload.hashCode());
  }

  @Test
  void referenceAndResultComponentsParticipateInEquality() {
    Map<Integer, Integer> atoms = mapping(0, 2, 1, 3);
    MCSResult result = new MCSResult(atoms, 2, 0.5, "CC");
    assertNotEquals(result, new MCSResult(mapping(0, 3, 1, 2), 2, 0.5, "CC"));
    assertNotEquals(result, new MCSResult(atoms, 1, 0.5, "CC"));
    assertNotEquals(result, new MCSResult(atoms, 2, 0.25, "CC"));
    assertNotEquals(result, new MCSResult(atoms, 2, 0.5, "CO"));
    SubstructureResult substructure = new SubstructureResult(true, list(atoms), stats());
    assertNotEquals(substructure, new SubstructureResult(false, list(atoms), stats()));
    assertNotEquals(substructure, new SubstructureResult(true, list(mapping(0, 1)), stats()));
    assertNotEquals(substructure, new SubstructureResult(true, list(atoms), null));
    MCSProfiledResult profile = new MCSProfiledResult(result, timers());
    assertNotEquals(profile, new MCSProfiledResult(null, timers()));
    assertNotEquals(profile, new MCSProfiledResult(result, null));
    Query query = new Query(true, "C*N", null);
    assertNotEquals(query, new Query(false, "C*N", null));
    assertNotEquals(query, new Query(true, "C*O", null));
    org.openscience.cdk.interfaces.IAtomContainer container =
        org.openscience.cdk.silent.SilentChemObjectBuilder.getInstance().newAtomContainer();
    Query molecule = new Query(true, "C*N", container);
    assertSame(container, molecule.container());
    assertEquals(molecule, new Query(true, "C*N", container));
    assertNotEquals(query, molecule);
  }

  @Test
  void everyNumericTelemetryComponentParticipatesInEquality() throws Exception {
    for (Object value : new Object[] {stats(), timers()}) {
      ObjectNode source = (ObjectNode) JSON.valueToTree(value);
      Iterator<Map.Entry<String, JsonNode>> fields = source.fields();
      while (fields.hasNext()) {
        Map.Entry<String, JsonNode> field = fields.next();
        ObjectNode changed = source.deepCopy();
        if (field.getValue().isBoolean()) changed.put(field.getKey(), !field.getValue().asBoolean());
        else changed.put(field.getKey(), field.getValue().asLong() + 1);
        Object alternative = JSON.treeToValue(changed, value.getClass());
        assertNotEquals(value, alternative, field.getKey());
      }
    }
  }
}

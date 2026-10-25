/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.StringJoiner;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/**
 * The built-in named SMARTS predicates: the atoms each one matches, and definitions identical to
 * those of the native core. The same cases are in cpp/tests/test_parsers.cpp.
 */
class NamedPredicateRegressionTest {

  /** The expansion of [$name] for each built-in predicate, as the native core writes it. */
  private static Map<String, String> expansions() {
    Map<String, String> expected = new LinkedHashMap<>();
    expected.put("isAmideN", "[$([#7X3][#6]=O)]");
    expected.put("isCarboxylC", "[$([CX3](=O)[OX2H1,OX1-])]");
    expected.put("isSulfonamideN", "[$([NX3;$(NS(=O)=O)])]");
    expected.put("isEster", "[$([#6X3;H1,$([#6][#6])](=O)[#8X2][#6;!$([#6]=O)])]");
    expected.put("isKetone", "[$([#6X3](=O)([#6])[#6])]");
    expected.put("isHalogen", "[$([F,Cl,Br,I])]");
    expected.put("isPositive", "[$([+,+2,+3,+4,+5,+6,+7,+8,+9,+10,+11,+12,+13,+14,+15])]");
    expected.put("isNegative", "[$([-,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15])]");
    expected.put("isHBDonor", "[$([#7,#8,#16;!H0])]");
    expected.put("isHBAcceptor", "[$([#7&!X3&!v5,#7&A&!v5&!$(NC=O),#8,#16&!v4&!v6;+0,-,-2,-3])]");
    return expected;
  }

  @Test
  void builtInPredicatesExpandToTheirDefinitions() {
    Standardiser.PredicateRegistry registry = new Standardiser.PredicateRegistry();
    for (Map.Entry<String, String> entry : expansions().entrySet())
      assertEquals(entry.getValue(),
          Standardiser.expandNamedPredicates("[$" + entry.getKey() + "]", registry), entry.getKey());
  }

  @Test
  void builtInPredicatesAreDefinedAsInTheNativeCore() throws Exception {
    Path header = Paths.get("cpp", "include", "smsd", "smarts_parser.hpp");
    // In a source checkout, where java/pom.xml sits beside cpp/, the header must be present.
    assumeTrue(Files.exists(Paths.get("java", "pom.xml")), "not run from a source checkout");
    assertTrue(Files.exists(header), header + " is missing from the source checkout");
    String source = new String(Files.readAllBytes(header), StandardCharsets.UTF_8);
    Matcher definition =
        Pattern.compile("registerPredicate\\(\"(\\w+)\",\\s*\"([^\"]*)\"\\)").matcher(source);
    Map<String, String> nativeDefinitions = new LinkedHashMap<>();
    while (definition.find()) nativeDefinitions.put(definition.group(1), definition.group(2));
    Standardiser.PredicateRegistry registry = new Standardiser.PredicateRegistry();
    Map<String, String> javaDefinitions = new LinkedHashMap<>();
    for (String name : expansions().keySet()) javaDefinitions.put(name, registry.get(name));
    assertEquals(javaDefinitions, nativeDefinitions);
  }

  /**
   * Names inside recursive SMARTS expand as in the native core, and a name ends at the first
   * character that is not an ASCII letter, digit or underscore; other text is kept.
   */
  @Test
  void namesInsideRecursiveSmartsExpand() throws Exception {
    Standardiser.PredicateRegistry registry = new Standardiser.PredicateRegistry();
    String query = Standardiser.expandNamedPredicates("[$(C[$isHalogen])]", registry);
    assertEquals("[$(C[$([F,Cl,Br,I])])]", query);
    assertTrue(Substructure.smarts(query).match(rawMolecule("CCl")).found());
    assertFalse(Substructure.smarts(query).match(rawMolecule("CCO")).found());
    assertEquals("[$([#6])]", Standardiser.expandNamedPredicates("[$([#6])]", registry));
    assertEquals("[$isUnknown]", Standardiser.expandNamedPredicates("[$isUnknown]", registry));
    assertEquals("[$([F,Cl,Br,I])\u00c4]",
        Standardiser.expandNamedPredicates("[$isHalogen\u00c4]", registry));
    registry.put("is\u00c4", "[#8]");
    assertEquals("[$is\u00c4]", Standardiser.expandNamedPredicates("[$is\u00c4]", registry));
  }

  /** Each predicate with a target and the target atoms it matches, none when the list is empty. */
  @ParameterizedTest(name = "{0} {1}")
  @CsvSource(value = {
      "isHBAcceptor|CC(=O)C|2",
      "isHBAcceptor|c1ccncc1|3",
      "isHBAcceptor|CO|1",
      "isHBAcceptor|CC[O-]|2",
      "isHBAcceptor|Nc1ccccc1|0",
      "isHBAcceptor|CC#N|2",
      "isHBAcceptor|CC(=O)N|2",
      "isHBAcceptor|c1c[nH]cn1|4",
      "isHBAcceptor|CSC|1",
      "isHBAcceptor|C[S-]|1",
      "isHBAcceptor|CS(=O)C|2",
      "isHBAcceptor|CS(=O)(=O)C|2 3",
      "isHBAcceptor|C[N+](C)(C)C|",
      "isHBAcceptor|c1cc[nH+]cc1|",
      "isHBAcceptor|Cn1cccc1|",
      "isHBAcceptor|CCl|",
      "isHBAcceptor|CN(=O)=O|2 3",
      "isHBAcceptor|CN(C)(C)=O|4",
      "isHBAcceptor|CC(=S)N|2 3",
      "isAmideN|CC(=O)N|3",
      "isAmideN|NC(=O)N|0 3",
      "isAmideN|O=C1CCCN1|5",
      "isAmideN|CN|",
      "isAmideN|CNC|",
      "isAmideN|Nc1ccccc1|",
      "isAmideN|CS(=O)(=O)N|",
      "isAmideN|NCC(C)=O|",
      "isAmideN|O=c1cc[nH]c(=O)[nH]1|4 7",
      "isAmideN|O=c1cccc[nH]1|6",
      "isAmideN|CC(=O)n1cccc1|3",
      "isAmideN|O=c1cc[nH]cc1|",
      "isPositive|[NH4+]|0",
      "isPositive|[Fe+2]|0",
      "isPositive|[Fe+3]|0",
      "isPositive|[Fe+10]|0",
      "isPositive|[NH3+]CC(=O)[O-]|0",
      "isPositive|CCO|",
      "isPositive|C[O-]|",
      "isNegative|C[O-]|1",
      "isNegative|[O-2]|0",
      "isNegative|[S-2]|0",
      "isNegative|[Fe-15]|0",
      "isNegative|[NH3+]CC(=O)[O-]|4",
      "isNegative|CCO|",
      "isNegative|[NH4+]|",
      "isKetone|CC(=O)C|1",
      "isKetone|CC(=O)c1ccccc1|1",
      "isKetone|O=c1ccoc2ccccc12|1",
      "isKetone|O=C1C=COC2=CC=CC=C12|1",
      "isKetone|CC=O|",
      "isKetone|CC(=O)O|",
      "isKetone|CC(=O)OC|",
      "isEster|CC(=O)OC|1",
      "isEster|O=COC|1",
      "isEster|O=C1CCCO1|1",
      "isEster|O=c1ccc2ccccc2o1|1",
      "isEster|O=C1C=CC2=CC=CC=C2O1|1",
      "isEster|CC(=O)O|",
      "isEster|COC(=O)OC|",
      "isEster|COC(=O)N|",
      "isEster|CC(=O)SC|",
      "isEster|CC(=O)OC(C)=O|",
      "isHBDonor|c1cc[nH]c1|3",
      "isHBDonor|CO|1",
      "isHBDonor|C[NH3+]|1",
      "isHBDonor|CS|1",
      "isHBDonor|COC|",
      "isHBDonor|c1ccncc1|",
      "isCarboxylC|CC(=O)O|1",
      "isCarboxylC|CC(=O)[O-]|1",
      "isCarboxylC|CNC(=O)O|2",
      "isCarboxylC|CC(=O)OC|",
      "isSulfonamideN|CS(=O)(=O)NC|4",
      "isSulfonamideN|CNC|",
      "isHalogen|FC(Cl)(Br)I|0 2 3 4",
      "isHalogen|CCO|"
  }, delimiter = '|')
  void builtInPredicatesMatchTheAtomsTheyName(String name, String target, String atoms)
      throws Exception {
    String query = Standardiser.expandNamedPredicates(
        "[$" + name + "]", new Standardiser.PredicateRegistry());
    TreeSet<Integer> matched = new TreeSet<>();
    for (SearchResult.Match match : Substructure.smarts(query).matchAll(rawMolecule(target)).matches())
      matched.add(match.atomMapping().get(0));
    StringJoiner found = new StringJoiner(" ");
    for (int atom : matched) found.add(String.valueOf(atom));
    assertEquals(atoms == null ? "" : atoms, found.toString());
  }
}

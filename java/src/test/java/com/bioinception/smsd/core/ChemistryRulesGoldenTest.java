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
import org.openscience.cdk.silent.SilentChemObjectBuilder;

/** Shared chemistry-rule fixture, also read by the C++ and Python tests. */
class ChemistryRulesGoldenTest {
  private static JsonNode rules;

  @BeforeAll
  static void load() throws Exception {
    try (InputStream in = ChemistryRulesGoldenTest.class.getClassLoader()
        .getResourceAsStream("com/bioinception/smsd/chemistry_rules.json")) {
      assertNotNull(in, "Missing shared chemistry rules resource");
      rules = new ObjectMapper().readTree(in);
    }
  }

  /** The molecule of a fixture, where an empty SMILES is an empty molecule. */
  private static IAtomContainer molecule(String smiles) throws Exception {
    if (smiles.isEmpty()) return SilentChemObjectBuilder.getInstance().newAtomContainer();
    return rawMolecule(smiles);
  }

  private static int[] ints(JsonNode values) {
    int[] result = new int[values.size()];
    for (int i = 0; i < result.length; i++) result[i] = values.get(i).asInt();
    return result;
  }

  private static List<Integer> list(JsonNode values) {
    List<Integer> result = new ArrayList<>();
    for (JsonNode value : values) result.add(value.asInt());
    return result;
  }

  private static long[] words(JsonNode positions, int size) {
    long[] words = new long[(size - 1) / 64 + 1];
    for (int position : ints(positions)) words[position >>> 6] |= 1L << (position & 63);
    return words;
  }

  private static List<Integer> positions(long[] words, int size) {
    List<Integer> result = new ArrayList<>();
    for (int bit = 0; bit < Math.min(size, words.length * 64); bit++)
      if ((words[bit >>> 6] & 1L << (bit & 63)) != 0) result.add(bit);
    return result;
  }

  @Test
  void similarityMetricsShareTheirEmptyFingerprintRules() {
    for (JsonNode item : rules.get("similarity").get("cases")) {
      String name = item.get("name").asText();
      int size = item.get("fp_size").asInt();
      long[] first = words(item.get("first"), size), second = words(item.get("second"), size);
      List<long[][]> forms = new ArrayList<>();
      forms.add(new long[][] {first, second});
      // An empty fingerprint is the same set whether it has zero words or none.
      if (item.get("first").size() == 0) forms.add(new long[][] {new long[0], second});
      if (item.get("second").size() == 0) forms.add(new long[][] {first, new long[0]});
      for (long[][] form : forms) {
        assertEquals(item.get("tanimoto").asDouble(), FingerprintEngine.tanimoto(form[0], form[1]), 1e-12, name);
        assertEquals(item.get("overlap").asDouble(), FingerprintEngine.overlapCoefficient(form[0], form[1]), 1e-12, name);
        assertEquals(item.get("dice").asDouble(), FingerprintEngine.dice(form[0], form[1]), 1e-12, name);
        assertEquals(item.get("cosine").asDouble(), FingerprintEngine.cosine(form[0], form[1]), 1e-12, name);
        assertEquals(item.get("soergel").asDouble(), FingerprintEngine.soergel(form[0], form[1]), 1e-12, name);
      }
    }
  }

  @Test
  void fingerprintQualitySharesItsDefinitions() {
    JsonNode section = rules.get("fingerprint_quality");
    double tolerance = section.get("tolerance").asDouble();
    for (JsonNode item : section.get("cases")) {
      String name = item.get("name").asText();
      int size = item.get("fp_size").asInt(), setBits = item.get("set_bits").asInt();
      // Positions 0 to set_bits - 1, with every padding bit beyond the size set as well.
      long[] words = new long[(size - 1) / 64 + 1];
      for (int bit = 0; bit < setBits; bit++) words[bit >>> 6] |= 1L << (bit & 63);
      if (size % 64 != 0) words[words.length - 1] |= -1L << (size % 64);
      Map<String, Object> quality = FingerprintEngine.analyzeFpQuality(words, size);
      assertEquals(Arrays.asList("set_bits", "fp_size", "density", "expected_collision_rate", "is_saturated"),
          new ArrayList<>(quality.keySet()), name);
      assertEquals(setBits, quality.get("set_bits"), name);
      assertEquals(size, quality.get("fp_size"), name);
      assertEquals(item.get("density").asDouble(), (Double) quality.get("density"), tolerance, name);
      assertEquals(item.get("expected_collision_rate").asDouble(),
          (Double) quality.get("expected_collision_rate"), tolerance, name);
      assertEquals(item.get("is_saturated").asBoolean(), quality.get("is_saturated"), name);
    }
  }

  @Test
  void circularSimilarityOfEmptyMolecules() throws Exception {
    for (JsonNode item : rules.get("circular_similarity").get("cases")) {
      long[] query = FingerprintEngine.ecfp(molecule(item.get("query").asText()), 2, 2048);
      for (int index = 0; index < item.get("targets").size(); index++) {
        String target = item.get("targets").get(index).asText();
        long[] fingerprint = FingerprintEngine.ecfp(molecule(target), 2, 2048);
        assertEquals(item.get("tanimoto").get(index).asDouble(), FingerprintEngine.tanimoto(query, fingerprint), 1e-12,
            item.get("query").asText() + " to " + target);
      }
    }
  }

  /** Every field of a profile as the fixture names it. */
  private static Map<String, String> profileFields(ChemOptions options) {
    Map<String, String> fields = new TreeMap<>();
    fields.put("match_atom_type", String.valueOf(options.matchAtomType));
    fields.put("match_formal_charge", String.valueOf(options.matchFormalCharge));
    fields.put("use_chirality", String.valueOf(options.useChirality));
    fields.put("use_bond_stereo", String.valueOf(options.useBondStereo));
    fields.put("ring_matches_ring_only", String.valueOf(options.ringMatchesRingOnly));
    fields.put("complete_rings_only", String.valueOf(options.completeRingsOnly));
    fields.put("match_isotope", String.valueOf(options.matchIsotope));
    fields.put("ring_fusion_mode", options.ringFusionMode.name());
    fields.put("tautomer_aware", String.valueOf(options.tautomerAware));
    fields.put("pH", String.valueOf(options.pH));
    fields.put("solvent", options.solvent.name());
    fields.put("match_bond_order", options.matchBondOrder.name());
    fields.put("aromaticity_mode", options.aromaticityMode.name());
    return fields;
  }

  @Test
  void chemistryProfilesChangeOnlyTheirListedFields() {
    JsonNode section = rules.get("chemistry_profiles");
    for (JsonNode item : section.get("cases")) {
      String name = item.get("name").asText();
      Map<String, String> expected = new TreeMap<>();
      section.get("defaults").fields().forEachRemaining(field -> expected.put(field.getKey(), field.getValue().asText()));
      item.get("changes").fields().forEachRemaining(field -> expected.put(field.getKey(), field.getValue().asText()));
      ChemOptions options = item.get("profile").asText().equals("tautomer_profile")
          ? ChemOptions.tautomerProfile() : ChemOptions.profile(item.get("profile").asText());
      assertEquals(expected, profileFields(options), name);
    }
  }

  /** The chemistry of a fixture: the defaults with the listed fields changed. */
  private static ChemOptions chemistry(JsonNode fields) {
    ChemOptions chemistry = new ChemOptions();
    if (fields.has("match_atom_type")) chemistry.matchAtomType = fields.get("match_atom_type").asBoolean();
    if (fields.has("match_formal_charge")) chemistry.matchFormalCharge = fields.get("match_formal_charge").asBoolean();
    if (fields.has("ring_matches_ring_only")) chemistry.ringMatchesRingOnly = fields.get("ring_matches_ring_only").asBoolean();
    if (fields.has("tautomer_aware")) chemistry.tautomerAware = fields.get("tautomer_aware").asBoolean();
    if (fields.has("aromaticity_mode"))
      chemistry.aromaticityMode = ChemOptions.AromaticityMode.valueOf(fields.get("aromaticity_mode").asText());
    return chemistry;
  }

  @Test
  void similarityUpperBoundUsesTheSharedAtomLabels() throws Exception {
    JsonNode section = rules.get("similarity_upper_bound");
    double tolerance = section.get("tolerance").asDouble();
    for (JsonNode item : section.get("cases")) {
      String name = item.get("name").asText();
      ChemOptions chemistry = chemistry(item.get("chemistry"));
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      assertEquals(item.get("bound").asDouble(), MCS.similarityUpperBound(query, target, chemistry), tolerance, name);
      assertEquals(item.get("bound").asDouble(), MCS.similarityUpperBound(target, query, chemistry), tolerance,
          name + " reversed");
    }
  }

  @Test
  void consensusReducesEachSeedInInputOrder() throws Exception {
    for (JsonNode item : rules.get("consensus").get("cases")) {
      String name = item.get("name").asText();
      List<MolGraph> molecules = new ArrayList<>();
      for (JsonNode smiles : item.get("molecules")) molecules.add(new MolGraph(molecule(smiles.asText())));
      MCS.Options options = new MCS.Options();
      if (item.has("objective")) options.objective = MCS.Objective.valueOf(item.get("objective").asText());
      MCS.ConsensusResult result = MCS.consensus(molecules, options, item.get("threshold").asDouble());
      assertEquals(item.get("found").asBoolean(), result.found(), name);
      assertTrue(result.complete(), name);
      if (!result.found()) continue;
      assertEquals(item.get("atoms").asInt(), result.fragment().atomCount(), name);
      assertEquals(item.get("bonds").asInt(), result.fragment().bondCount(), name);
      assertEquals(item.get("support").asInt(), result.supportCount(), name);
      assertEquals(item.get("source_index").asInt(), result.sourceIndex(), name);
      assertEquals(list(item.get("mappings")), new ArrayList<>(result.mappings().keySet()), name);
      assertEquals(item.get("atoms").asInt(), result.graph().n, name);
    }
  }

  @Test
  void implicitHydrogensFollowTheSharedValenceTable() {
    for (JsonNode item : rules.get("implicit_hydrogens").get("cases")) {
      assertEquals(item.get("hydrogens").asInt(),
          MolGraph.defaultHydrogenCount(item.get("atomic_number").asInt(), item.get("aromatic").asBoolean(),
              item.get("bond_order_sum").asInt(), item.get("charge").asInt()),
          item.get("note").asText());
    }
  }

  private static List<Integer> matchedAtoms(String smarts, MolGraph target) {
    List<Integer> atoms = new ArrayList<>();
    for (SearchResult.Match match : Substructure.smarts(smarts).matchAll(target).matches())
      atoms.add(match.atomMapping().get(0));
    Collections.sort(atoms);
    return atoms;
  }

  /** The same graph built without hydrogen counts, so every count follows the valence model. */
  private static MolGraph withoutHydrogenCounts(MolGraph graph) {
    int[][] orders = new int[graph.n][];
    boolean[][] aromaticBonds = new boolean[graph.n][];
    for (int atom = 0; atom < graph.n; atom++) {
      orders[atom] = new int[graph.neighbors[atom].length];
      aromaticBonds[atom] = new boolean[graph.neighbors[atom].length];
      for (int k = 0; k < graph.neighbors[atom].length; k++) {
        orders[atom][k] = graph.bondOrder(atom, graph.neighbors[atom][k]);
        aromaticBonds[atom][k] = graph.bondAromatic(atom, graph.neighbors[atom][k]);
      }
    }
    return new MolGraph.Builder().atomCount(graph.n).atomicNumbers(graph.atomicNum.clone())
        .formalCharges(graph.formalCharge.clone()).ringFlags(graph.ring.clone()).aromaticFlags(graph.aromatic.clone())
        .neighbors(graph.neighbors).bondOrders(orders).bondAromaticFlags(aromaticBonds).build();
  }

  @Test
  void smartsHydrogenConnectivityAndValence() throws Exception {
    for (JsonNode item : rules.get("smarts").get("cases"))
      assertEquals(list(item.get("atoms")),
          matchedAtoms(item.get("smarts").asText(), new MolGraph(molecule(item.get("smiles").asText()))),
          item.get("name").asText() + " " + item.get("smarts").asText());
  }

  @Test
  void smartsOnGraphsWithoutHydrogenCounts() throws Exception {
    for (JsonNode item : rules.get("smarts_without_hydrogen_counts").get("cases")) {
      MolGraph graph = withoutHydrogenCounts(new MolGraph(molecule(item.get("smiles").asText())));
      assertEquals(list(item.get("atoms")), matchedAtoms(item.get("smarts").asText(), graph),
          item.get("name").asText() + " " + item.get("smarts").asText());
    }
  }

  @Test
  @SuppressWarnings("deprecation")
  void tautomerConfidenceFormsAgree() throws Exception {
    double tolerance = rules.get("tautomer_confidence").get("tolerance").asDouble();
    Map<String, Double> confidences = new HashMap<>();
    for (JsonNode item : rules.get("tautomer_confidence").get("cases")) {
      String name = item.get("name").asText();
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      Map<Integer, Integer> atoms = new LinkedHashMap<>();
      for (JsonNode pair : item.get("atoms")) atoms.put(pair.get(0).asInt(), pair.get(1).asInt());
      ChemOptions chemistry = ChemOptions.tautomerProfile()
          .withSolvent(ChemOptions.Solvent.valueOf(item.get("solvent").asText()))
          .withPH(item.get("pH").asDouble());
      if (item.has("error")) {
        assertThrows(IllegalArgumentException.class, () -> MCS.computeTautomerConfidence(query, target, atoms), name);
        SearchResult.Match unchecked = new SearchResult.Match(atoms, Collections.emptyMap(), Collections.emptyList(), atoms.size());
        assertThrows(IllegalArgumentException.class,
            () -> MCS.computeTautomerConfidence(query, target, unchecked, chemistry), name);
        continue;
      }
      List<SearchResult.Bond> bonds = new ArrayList<>();
      for (Map.Entry<Integer, Integer> pair : atoms.entrySet())
        for (int neighbor : query.neighbors[pair.getKey()]) {
          Integer image = atoms.get(neighbor);
          if (neighbor > pair.getKey() && image != null && target.hasBond(pair.getValue(), image))
            bonds.add(new SearchResult.Bond(pair.getKey(), neighbor));
        }
      SearchResult.Match match = MCS.fragment(query, atoms, bonds);
      double confidence = MCS.computeTautomerConfidence(query, target, match, chemistry);
      confidences.put(name, confidence);
      assertEquals(item.get("confidence").asDouble(), confidence, tolerance, name);
      if (item.has("three_argument"))
        assertEquals(item.get("three_argument").asDouble(), MCS.computeTautomerConfidence(query, target, atoms),
            tolerance, name + " (three-argument form)");
    }
    // A case with same_as differs from the named case only in the order of the atoms.
    for (JsonNode item : rules.get("tautomer_confidence").get("cases")) {
      if (!item.has("same_as")) continue;
      String name = item.get("name").asText(), reference = item.get("same_as").asText();
      assertTrue(confidences.containsKey(reference), name + " names no case");
      assertEquals(confidences.get(reference), confidences.get(name), tolerance, name + " against " + reference);
    }
  }

  private static int[][] rows(JsonNode rows) {
    int[][] result = new int[rows.size()][];
    for (int row = 0; row < result.length; row++) result[row] = ints(rows.get(row));
    return result;
  }

  @Test
  void tetrahedralCentresShareTheirLigandRule() throws Exception {
    JsonNode section = rules.get("tetrahedral_centres");
    for (JsonNode item : section.get("graphs")) {
      int[] elements = ints(item.get("atomic_numbers"));
      int atom = item.get("atom").asInt();
      int[] tags = new int[elements.length];
      tags[atom] = 1;
      MolGraph.Builder builder = new MolGraph.Builder().atomCount(elements.length).atomicNumbers(elements)
          .neighbors(rows(item.get("neighbors"))).bondOrders(rows(item.get("orders"))).tetrahedralChirality(tags);
      if (item.has("hydrogens")) builder.hydrogenCounts(ints(item.get("hydrogens")));
      int[] ligands = builder.build().tetrahedralLigands(atom);
      String name = item.get("name").asText();
      if (item.get("ligands").isNull()) assertNull(ligands, name);
      else assertArrayEquals(ints(item.get("ligands")), ligands, name);
    }
    Substructure.Options options = new Substructure.Options();
    options.chemistry.useChirality = true;
    for (JsonNode item : section.get("matches")) {
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      assertEquals(item.get("found").asBoolean(),
          Substructure.find(query).withOptions(options).match(target).found(), item.get("name").asText());
    }
    for (JsonNode item : section.get("labels")) {
      String smiles = item.get("smiles").asText();
      assertEquals(expectedLabels(item), cipLabels(smiles), smiles);
    }
  }

  /** The R/S labels of a molecule by atom, as the fixture writes them. */
  private static Map<String, String> cipLabels(String smiles) throws Exception {
    Map<String, String> found = new TreeMap<>();
    for (Map.Entry<Integer, Character> label : CIPAssigner.assignRS(new MolGraph(molecule(smiles))).entrySet())
      if (label.getValue() == 'R' || label.getValue() == 'S')
        found.put(String.valueOf(label.getKey()), String.valueOf(label.getValue()));
    return found;
  }

  /** All four descriptors of a molecule by atom, the pseudoasymmetric r and s included. */
  private static Map<String, String> allCipLabels(String smiles) throws Exception {
    Map<String, String> found = new TreeMap<>();
    for (Map.Entry<Integer, Character> label : CIPAssigner.assignRS(new MolGraph(molecule(smiles))).entrySet())
      found.put(String.valueOf(label.getKey()), String.valueOf(label.getValue()));
    return found;
  }

  private static Map<String, String> expectedLabels(JsonNode item) {
    Map<String, String> expected = new TreeMap<>();
    item.get("labels").fields().forEachRemaining(label -> expected.put(label.getKey(), label.getValue().asText()));
    return expected;
  }

  @Test
  void mancudeCipLabelsDoNotDependOnTheKekuleStructure() throws Exception {
    for (JsonNode item : rules.get("cip_mancude").get("cases"))
      assertEquals(expectedLabels(item), cipLabels(item.get("smiles").asText()), item.get("name").asText());
  }

  @Test
  void auxiliaryDescriptorsFollowSequenceRules4And5() throws Exception {
    for (JsonNode item : rules.get("cip_auxiliary_descriptors").get("cases"))
      assertEquals(expectedLabels(item), allCipLabels(item.get("smiles").asText()), item.get("name").asText());
  }

  @Test
  void aromaticityIsPerceivedFromTheStructure() throws Exception {
    for (JsonNode item : rules.get("aromatic_input").get("cases")) {
      MolGraph graph = new MolGraph(Standardiser.standardise(molecule(item.get("smiles").asText()),
          Standardiser.TautomerMode.NONE));
      int aromatic = 0;
      for (int atom = 0; atom < graph.n; atom++) if (graph.aromatic[atom]) aromatic++;
      assertEquals(item.get("aromatic_atoms").asInt(), aromatic, item.get("name").asText());
    }
  }

  @Test
  void firstPassRanksByRule1aAndLeavesTheRestToTheHierarchicalRules() throws Exception {
    for (JsonNode item : rules.get("cip_hierarchical_ranking").get("cases")) {
      String name = item.get("name").asText();
      assertEquals(expectedLabels(item), allCipLabels(item.get("smiles").asText()), name);
      MolGraph graph = new MolGraph(molecule(item.get("smiles").asText()));
      Map<Long, Character> found = CIPAssigner.assignEZ(graph);
      Map<String, String> bonds = new TreeMap<>(), expected = new TreeMap<>();
      for (int first = 0; first < graph.n; first++)
        for (int second : graph.neighbors[first])
          if (second > first && found.containsKey(MolGraph.bondKey(first, second)))
            bonds.put(first + "-" + second, String.valueOf(found.get(MolGraph.bondKey(first, second))));
      item.get("bonds").fields().forEachRemaining(bond -> expected.put(bond.getKey(), bond.getValue().asText()));
      assertEquals(expected, bonds, name + " E/Z");
    }
  }

  @Test
  void smartsFollowsDaylightSemantics() throws Exception {
    for (JsonNode item : rules.get("smarts_daylight").get("cases")) {
      MolGraph target = new MolGraph(molecule(item.get("smiles").asText()));
      assertEquals(item.get("embeddings").asInt(),
          Substructure.smarts(item.get("smarts").asText()).matchAll(target).matches().size(),
          item.get("name").asText() + " " + item.get("smarts").asText());
    }
    for (JsonNode item : rules.get("smarts_daylight").get("rejected"))
      assertThrows(IllegalArgumentException.class, () -> Substructure.smarts(item.get("smarts").asText()),
          item.get("name").asText() + " " + item.get("smarts").asText());
  }

  /** The options of a fixture case: the classic configuration of its family, then its keyword arguments. */
  private static FingerprintOptions fingerprintOptionsOf(JsonNode spec) {
    FingerprintOptions options = spec.has("maccs") && spec.get("maccs").asBoolean()
        ? FingerprintOptions.maccs()
        : FingerprintOptions.of(FingerprintOptions.Family.valueOf(spec.get("family").asText().toUpperCase(Locale.ROOT)));
    Iterator<Map.Entry<String, JsonNode>> fields = spec.fields();
    while (fields.hasNext()) {
      Map.Entry<String, JsonNode> field = fields.next();
      String key = field.getKey();
      JsonNode value = field.getValue();
      switch (key) {
        case "family": case "maccs": break;
        case "radius": options.radius = value.asInt(); break;
        case "path_length": options.pathLength = value.asInt(); break;
        case "size": options.size = value.asInt(); break;
        case "ring_size_cap": options.ringSizeCap = value.asInt(); break;
        case "element": options.element = value.asBoolean(); break;
        case "degree": options.degree = value.asBoolean(); break;
        case "valence": options.valence = value.asBoolean(); break;
        case "formal_charge": options.formalCharge = value.asBoolean(); break;
        case "isotope": options.isotope = value.asBoolean(); break;
        case "hydrogen_count": options.hydrogenCount = value.asBoolean(); break;
        case "tautomer_class": options.tautomerClass = value.asBoolean(); break;
        case "aromaticity": options.aromaticity = value.asBoolean(); break;
        case "bond_order": options.bondOrder = value.asBoolean(); break;
        case "aromatic_bond": options.aromaticBond = value.asBoolean(); break;
        case "ring_bond": options.ringBond = value.asBoolean(); break;
        case "ring_membership": options.ringMembership = value.asBoolean(); break;
        case "ring_size": options.ringSize = value.asBoolean(); break;
        case "cip_atoms": options.cipAtoms = value.asBoolean(); break;
        case "cip_bonds": options.cipBonds = value.asBoolean(); break;
        case "fragments":
          options.fragments = FingerprintOptions.Fragments.valueOf(value.asText().toUpperCase(Locale.ROOT));
          break;
        case "scaffold":
          options.scaffold = FingerprintOptions.Scaffold.valueOf(value.asText().toUpperCase(Locale.ROOT));
          break;
        case "keys":
          options.keys = new ArrayList<>();
          for (JsonNode item : value) options.keys.add(new FingerprintOptions.Key(item.get(0).asText(), item.get(1).asInt()));
          options.size = options.keys.size();
          break;
        default:
          throw new IllegalArgumentException("Unknown fingerprint option " + key);
      }
    }
    return options;
  }

  private static List<List<Integer>> nonZero(int[] counts) {
    List<List<Integer>> pairs = new ArrayList<>();
    for (int i = 0; i < counts.length; i++) if (counts[i] != 0) pairs.add(Arrays.asList(i, counts[i]));
    return pairs;
  }

  private static List<List<Integer>> pairs(JsonNode items) {
    List<List<Integer>> pairs = new ArrayList<>();
    for (JsonNode item : items) pairs.add(Arrays.asList(item.get(0).asInt(), item.get(1).asInt()));
    return pairs;
  }

  @Test
  void fingerprintOptionsBitsAreSharedAcrossLanguages() throws Exception {
    for (JsonNode item : rules.get("fingerprint_options").get("cases")) {
      String name = item.get("name").asText();
      IAtomContainer container = molecule(item.get("smiles").asText());
      MolGraph graph = FingerprintEngine.toMolGraph(container);
      FingerprintOptions options = fingerprintOptionsOf(item.get("options"));
      int size = options.family == FingerprintOptions.Family.KEYS ? options.keys.size() : options.size;
      long[] bits = FingerprintEngine.fingerprint(graph, options);
      assertEquals(list(item.get("bits")), positions(bits, size), name + " bits");
      int[] counts = FingerprintEngine.fingerprintCounts(graph, options);
      assertEquals(pairs(item.get("counts")), nonZero(counts), name + " counts");
      if (!item.has("classic") || !item.get("classic").asBoolean()) continue;
      // The classic configuration of a family gives the bits of its fixed method.
      switch (options.family) {
        case ECFP:
          assertArrayEquals(FingerprintEngine.ecfp(graph, options.radius, size), bits, name + " equals ecfp");
          assertArrayEquals(FingerprintEngine.ecfpCounts(graph, options.radius, size), counts, name + " equals ecfpCounts");
          break;
        case FCFP:
          assertArrayEquals(FingerprintEngine.fcfp(graph, options.radius, size), bits, name + " equals fcfp");
          assertArrayEquals(FingerprintEngine.fcfpCounts(graph, options.radius, size), counts, name + " equals fcfpCounts");
          break;
        case PATH:
          assertArrayEquals(FingerprintEngine.pathFingerprint(container, options.pathLength, size), bits,
              name + " equals pathFingerprint");
          break;
        case TORSION:
          assertArrayEquals(FingerprintEngine.topologicalTorsion(graph, size), bits, name + " equals topologicalTorsion");
          assertArrayEquals(FingerprintEngine.topologicalTorsionCounts(graph, size), counts,
              name + " equals topologicalTorsionCounts");
          break;
        default:
          break;
      }
    }
  }

  @Test
  void explicitAromaticBondsTakeTheSharedHydrogenCounts() throws Exception {
    for (JsonNode item : rules.get("aromatic_bond_hydrogens").get("cases")) {
      MolGraph graph = new MolGraph(molecule(item.get("smiles").asText()));
      List<Integer> hydrogens = new ArrayList<>();
      for (int atom = 0; atom < graph.n; atom++) hydrogens.add(graph.hydrogenCount(atom));
      assertEquals(list(item.get("hydrogens")), hydrogens, item.get("smiles").asText());
    }
  }

  private static ChemOptions stereoChemistry(JsonNode item) {
    ChemOptions chemistry = new ChemOptions();
    chemistry.useChirality = item.get("chirality").asBoolean();
    chemistry.useBondStereo = item.get("bond_stereo").asBoolean();
    return chemistry;
  }

  @Test
  void substructureEmbeddingsNeedSpecifiedStereo() throws Exception {
    JsonNode section = rules.get("stereo_embeddings");
    for (JsonNode item : section.get("cases")) {
      Substructure.Options options = new Substructure.Options();
      options.chemistry = stereoChemistry(item);
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      Substructure prepared = Substructure.find(query).withOptions(options);
      assertEquals(item.get("found").asBoolean(), prepared.match(target).found(), item.get("name").asText());
      assertEquals(item.get("found").asBoolean(), prepared.matchAll(target).found(), item.get("name").asText());
    }
    for (JsonNode item : section.get("common_subgraphs")) {
      MCS.Options options = new MCS.Options();
      options.chemistry = stereoChemistry(item);
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertEquals(item.get("atoms").asInt(), result.found() ? result.matches().get(0).atomCount() : 0,
          item.get("name").asText());
    }
  }

  @Test
  void tautomerConsistencyCountsStoredHydrogens() throws Exception {
    for (JsonNode item : rules.get("tautomer_consistency").get("cases")) {
      Map<Integer, Integer> mapping = new TreeMap<>();
      for (JsonNode pair : item.get("atoms")) mapping.put(pair.get(0).asInt(), pair.get(1).asInt());
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      assertEquals(item.get("consistent").asBoolean(), MCS.validateTautomerConsistency(query, target, mapping),
          item.get("name").asText());
    }
  }

  @Test
  void tautomerClassesDoNotDependOnTheAtomOrder() throws Exception {
    for (JsonNode item : rules.get("tautomer_classes").get("cases")) {
      MolGraph graph = new MolGraph(molecule(item.get("smiles").asText()));
      assertArrayEquals(ints(item.get("classes")), graph.tautomerClass, item.get("name").asText());
    }
  }

  @Test
  void murckoScaffoldsCapStrippedAtomsWithHydrogen() throws Exception {
    for (JsonNode item : rules.get("murcko_scaffolds").get("cases")) {
      IAtomContainer scaffold = MCS.murckoScaffold(molecule(item.get("smiles").asText()));
      List<Integer> hydrogens = new ArrayList<>();
      for (org.openscience.cdk.interfaces.IAtom atom : scaffold.atoms()) hydrogens.add(atom.getImplicitHydrogenCount());
      assertEquals(list(item.get("hydrogens")), hydrogens, item.get("smiles").asText());
      assertEquals(item.get("scaffold").asText(), new MolGraph(scaffold).toCanonicalSmiles(), item.get("smiles").asText());
    }
  }

  @Test
  void consensusFollowsTheNativeReductionOrder() throws Exception {
    for (JsonNode item : rules.get("consensus_native_order").get("cases")) {
      String name = item.get("name").asText();
      List<MolGraph> molecules = new ArrayList<>();
      for (JsonNode smiles : item.get("molecules")) molecules.add(new MolGraph(molecule(smiles.asText())));
      MCS.Options options = new MCS.Options();
      options.objective = MCS.Objective.valueOf(item.get("objective").asText());
      options.connectivity = MCS.Connectivity.valueOf(item.get("connectivity").asText());
      MCS.ConsensusResult result = MCS.consensus(molecules, options, item.get("threshold").asDouble());
      assertEquals(item.get("found").asBoolean(), result.found(), name);
      assertEquals(item.get("source_index").asInt(), result.sourceIndex(), name);
      assertEquals(list(item.get("atoms")), new ArrayList<>(result.fragment().atomMapping().keySet()), name);
      assertEquals(item.get("bonds").asInt(), result.fragment().bondCount(), name);
      assertEquals(list(item.get("support")), new ArrayList<>(result.mappings().keySet()), name);
      assertEquals(list(item.get("support")).size(), result.supportCount(), name);
      assertEquals(item.get("smiles").asText(), result.toSmiles(), name);
    }
  }

  @Test
  void strictRingFusionComparesRelevantCycleCounts() throws Exception {
    JsonNode section = rules.get("strict_ring_fusion");
    for (JsonNode item : section.get("ring_counts")) {
      MolGraph graph = new MolGraph(molecule(item.get("smiles").asText()));
      assertArrayEquals(ints(item.get("counts")), graph.getRingCounts(), item.get("name").asText());
      assertEquals(item.get("truncated").asBoolean(), graph.ringsTruncated(), item.get("name").asText());
    }
    for (JsonNode item : section.get("common_bonds")) {
      String name = item.get("name").asText();
      MCS.Options options = new MCS.Options();
      options.objective = MCS.Objective.MCES;
      options.chemistry.ringFusionMode = ChemOptions.RingFusionMode.STRICT;
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      for (int turn = 0; turn < 2; turn++) {
        SearchResult result = turn == 0 ? MCS.find(query).withOptions(options).match(target)
                                        : MCS.find(target).withOptions(options).match(query);
        assertTrue(result.optimal(), name);
        assertEquals(item.get("bonds").asInt(), result.matches().get(0).bondCount(), name);
        assertEquals(item.get("atoms").asInt(), result.matches().get(0).atomCount(), name);
      }
    }
  }

  @Test
  void doubleBondLigandsAreJoinedToTheirEnd() throws Exception {
    for (JsonNode item : rules.get("double_bond_ligands").get("cases")) {
      String name = item.get("name").asText();
      MCS.Options options = new MCS.Options();
      options.objective = MCS.Objective.MCES;
      options.connectivity = "connected".equals(item.get("connectivity").asText())
          ? MCS.Connectivity.Connected : MCS.Connectivity.Fragments;
      options.chemistry.useBondStereo = true;
      MolGraph query = new MolGraph(molecule(item.get("query").asText()));
      MolGraph target = new MolGraph(molecule(item.get("target").asText()));
      SearchResult result = MCS.find(query).withOptions(options).match(target);
      assertTrue(result.optimal() && result.found(), name);
      assertEquals(item.get("bonds").asInt(), result.matches().get(0).bondCount(), name);
      assertEquals(item.get("atoms").asInt(), result.matches().get(0).atomCount(), name);
      assertEquals(Collections.emptyList(), MCS.validateMapping(query, target, result.matches().get(0), options), name);
    }
  }

  @Test
  void elementPathAndMcsFingerprintBitsAreUnchanged() throws Exception {
    for (JsonNode item : rules.get("path_fingerprints").get("cases")) {
      String name = item.get("name").asText();
      IAtomContainer container = molecule(item.get("smiles").asText());
      int length = item.get("path_length").asInt(), size = item.get("path_fp_size").asInt();
      assertEquals(list(item.get("path_fingerprint")),
          positions(FingerprintEngine.pathFingerprint(container, length, size), size), name + " path bits");
      ChemOptions chemistry = new ChemOptions();
      chemistry.tautomerAware = item.get("tautomer_aware").asBoolean();
      assertEquals(list(item.get("mcs_fingerprint")),
          positions(FingerprintEngine.mcsFingerprint(container, chemistry, length, size), size), name + " MCS bits");
    }
  }
}

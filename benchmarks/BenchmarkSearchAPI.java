/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Standardiser;
import com.bioinception.smsd.core.Substructure;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/** Shared-corpus diagnostics for the public prepared search APIs, on Java 8 or later. */
public final class BenchmarkSearchAPI {
  private static final ObjectMapper JSON = new ObjectMapper();
  private static final SmilesParser PARSER = new SmilesParser(SilentChemObjectBuilder.getInstance());
  private static final Set<String> PHASES = new HashSet<>(
      Arrays.asList("prepare", "prepared", "end_to_end", "batch"));
  private static volatile Object retained;

  private BenchmarkSearchAPI() {}

  private static final class Case {
    final String id, operation, query, target;
    final List<String> phases;
    final int limit, batchSize;
    final ChemOptions chemistry;
    final MCS.Options maximum;
    final Substructure.Options substructure;

    Case(String[] columns) {
      id = columns[0];
      operation = columns[1];
      if (!Arrays.asList("mcs", "sub", "smarts").contains(operation))
        throw new IllegalArgumentException("Unknown operation " + operation);
      chemistry = chemistry(columns[4]);
      if (!columns[5].equals("0") && !columns[5].equals("1"))
        throw new IllegalArgumentException("induced must be 0 or 1");
      long timeout = Long.parseLong(columns[6]);
      limit = Integer.parseInt(columns[7]);
      batchSize = Integer.parseInt(columns[8]);
      if (timeout < 0 || limit < 0 || batchSize < 1)
        throw new IllegalArgumentException("Invalid timeout, limit or batch size");
      query = columns[9];
      target = columns[10];
      maximum = new MCS.Options();
      maximum.chemistry = chemistry;
      maximum.objective = MCS.Objective.valueOf(columns[2]);
      maximum.connectivity = MCS.Connectivity.valueOf(columns[3]);
      maximum.timeoutMillis = timeout;
      maximum.nodeLimit = 0;
      maximum.minFragmentAtoms = 1;
      maximum.maxFragments = Integer.MAX_VALUE;
      if (!columns[11].isEmpty()) {
        String[] terms = columns[11].split(",");
        maximum.atomWeights = new double[terms.length];
        for (int atom = 0; atom < terms.length; atom++) {
          double weight = Double.parseDouble(terms[atom]);
          if (!Double.isFinite(weight)) throw new IllegalArgumentException("Non-finite weight");
          maximum.atomWeights[atom] = weight;
        }
      }
      substructure = new Substructure.Options();
      substructure.chemistry = chemistry;
      substructure.induced = columns[5].equals("1");
      substructure.timeoutMillis = timeout;
      substructure.nodeLimit = 0;
      phases = Arrays.asList(columns[12].split(","));
      if (phases.isEmpty() || !PHASES.containsAll(phases)
          || new HashSet<>(phases).size() != phases.size())
        throw new IllegalArgumentException("Invalid phases " + columns[12]);
      if (phases.contains("batch") && limit != 0)
        throw new IllegalArgumentException("Batch phase requires limit 0");
    }
  }

  private static ChemOptions chemistry(String profile) {
    ChemOptions options = new ChemOptions();
    options.matchAtomType = true;
    options.matchFormalCharge = false;
    options.matchIsotope = false;
    options.useChirality = false;
    options.useBondStereo = false;
    options.ringMatchesRingOnly = false;
    options.completeRingsOnly = false;
    options.ringFusionMode = ChemOptions.RingFusionMode.IGNORE;
    options.tautomerAware = false;
    options.matchBondOrder = ChemOptions.BondOrderMode.STRICT;
    options.aromaticityMode = ChemOptions.AromaticityMode.FLEXIBLE;
    switch (profile) {
      case "default": break;
      case "rings": options.completeRingsOnly = true;
      case "strict":
        options.aromaticityMode = ChemOptions.AromaticityMode.STRICT;
        options.matchFormalCharge = true;
        options.matchIsotope = true;
        options.ringMatchesRingOnly = true;
        break;
      case "any": options.matchBondOrder = ChemOptions.BondOrderMode.ANY; break;
      case "tautomer":
        options.tautomerAware = true;
        options.matchBondOrder = ChemOptions.BondOrderMode.LOOSE;
        break;
      case "stereo": options.useChirality = options.useBondStereo = true; break;
      default: throw new IllegalArgumentException("Unknown profile " + profile);
    }
    return options;
  }

  private static final class Molecule {
    final IAtomContainer container;
    final MolGraph graph;
    Molecule(String smiles) throws Exception {
      container = Standardiser.standardise(
          PARSER.parseSmiles(smiles), Standardiser.TautomerMode.NONE);
      graph = new MolGraph(container);
    }
  }

  private static final class Inputs {
    final Molecule query, target;
    final List<MolGraph> batch;
    Inputs(Case row) throws Exception {
      query = row.operation.equals("smarts") ? null : new Molecule(row.query);
      target = new Molecule(row.target);
      batch = Collections.nCopies(row.batchSize, target.graph);
    }
  }

  private static final class Prepared {
    final MCS maximum;
    final Substructure substructure;
    final Case row;
    Prepared(Case row, Inputs inputs) {
      this.row = row;
      maximum = row.operation.equals("mcs")
          ? MCS.find(inputs.query.graph).withOptions(row.maximum) : null;
      substructure = maximum != null ? null
          : (row.operation.equals("smarts") ? Substructure.smarts(row.query)
                                              : Substructure.find(inputs.query.graph))
                .withOptions(row.substructure);
    }
    SearchResult match(MolGraph target) {
      return maximum != null
          ? (row.limit == 0 ? maximum.match(target) : maximum.matchAll(target, row.limit))
          : (row.limit == 0 ? substructure.match(target)
                            : substructure.matchAll(target, row.limit));
    }
    List<SearchResult> batch(List<MolGraph> targets) {
      return maximum != null ? maximum.matchMany(targets) : substructure.matchMany(targets);
    }
  }

  private static final class Measurement {
    final long elapsed;
    final List<SearchResult> results;
    final Object prepared;
    Measurement(long elapsed, List<SearchResult> results, Object prepared) {
      this.elapsed = elapsed;
      this.results = results;
      this.prepared = prepared;
    }
  }

  private static Measurement measure(Case row, String phase, Inputs inputs, Prepared query)
      throws Exception {
    long started = System.nanoTime();
    if (phase.equals("prepare")) {
      Prepared prepared = new Prepared(row, inputs);
      long elapsed = System.nanoTime() - started;
      return new Measurement(elapsed, Collections.emptyList(), prepared);
    }
    if (phase.equals("batch")) {
      List<SearchResult> results = query.batch(inputs.batch);
      long elapsed = System.nanoTime() - started;
      return new Measurement(elapsed, results, query);
    }
    if (phase.equals("end_to_end")) {
      Inputs parsed = new Inputs(row);
      Prepared prepared = new Prepared(row, parsed);
      SearchResult result = prepared.match(parsed.target.graph);
      long elapsed = System.nanoTime() - started;
      return new Measurement(elapsed, Collections.singletonList(result), prepared);
    }
    SearchResult result = query.match(inputs.target.graph);
    long elapsed = System.nanoTime() - started;
    return new Measurement(elapsed, Collections.singletonList(result), query);
  }

  private static Map<String, Object> graph(Molecule molecule) {
    if (molecule == null) return null;
    List<Map<String, Object>> atoms = new ArrayList<>();
    for (int index = 0; index < molecule.container.getAtomCount(); index++) {
      IAtom atom = molecule.container.getAtom(index);
      Map<String, Object> item = new LinkedHashMap<>();
      item.put("element", atom.getAtomicNumber());
      item.put("charge", atom.getFormalCharge() == null ? 0 : atom.getFormalCharge());
      item.put("isotope", atom.getMassNumber() == null ? 0 : atom.getMassNumber());
      item.put("hydrogens", molecule.graph.hydrogenCount(index));
      item.put("aromatic", atom.isAromatic());
      item.put("ring", atom.isInRing());
      atoms.add(item);
    }
    List<List<Object>> bonds = new ArrayList<>();
    for (IBond bond : molecule.container.bonds()) {
      int first = molecule.container.indexOf(bond.getBegin());
      int second = molecule.container.indexOf(bond.getEnd());
      int order = bond.isAromatic() ? 4
          : bond.getOrder() == null || bond.getOrder() == IBond.Order.UNSET ? 1
                                                                         : bond.getOrder().numeric();
      bonds.add(Arrays.<Object>asList(Math.min(first, second), Math.max(first, second), order,
          bond.isAromatic(), bond.isInRing()));
    }
    bonds.sort((left, right) -> {
      int first = Integer.compare((Integer) left.get(0), (Integer) right.get(0));
      return first != 0 ? first : Integer.compare((Integer) left.get(1), (Integer) right.get(1));
    });
    Map<String, Object> result = new LinkedHashMap<>();
    result.put("atoms", atoms);
    result.put("bonds", bonds);
    return result;
  }

  private static Map<String, Object> result(SearchResult search) {
    Map<String, Object> result = new LinkedHashMap<>();
    result.put("status", search.status().toString());
    result.put("optimal", search.optimal());
    result.put("complete", search.complete());
    result.put("score", search.score());
    result.put("nodes_visited", search.nodesVisited());
    List<Map<String, Object>> matches = new ArrayList<>();
    for (SearchResult.Match match : search.matches()) {
      Map<String, Object> item = new LinkedHashMap<>();
      List<List<Integer>> atoms = new ArrayList<>(), bonds = new ArrayList<>();
      for (Map.Entry<Integer, Integer> atom : match.atomMapping().entrySet())
        atoms.add(Arrays.asList(atom.getKey(), atom.getValue()));
      for (Map.Entry<SearchResult.Bond, SearchResult.Bond> bond : match.bondMapping().entrySet())
        bonds.add(Arrays.asList(bond.getKey().first(), bond.getKey().second(),
            match.atomMapping().get(bond.getKey().first()),
            match.atomMapping().get(bond.getKey().second())));
      item.put("atoms", atoms);
      item.put("bonds", bonds);
      item.put("fragments", match.fragments());
      matches.add(item);
    }
    result.put("matches", matches);
    return result;
  }

  private static void validate(Case row, Inputs inputs, List<SearchResult> results) {
    for (SearchResult search : results)
      for (SearchResult.Match match : search.matches()) {
        if (new HashSet<>(match.atomMapping().values()).size() != match.atomCount())
          throw new IllegalStateException("Non-injective result");
        for (int targetAtom : match.atomMapping().values())
          if (targetAtom < 0 || targetAtom >= inputs.target.graph.atomCount())
            throw new IllegalStateException("Target index outside molecule");
        if (inputs.query == null) continue;
        List<String> errors = row.operation.equals("mcs")
            ? MCS.validateMapping(inputs.query.graph, inputs.target.graph, match, row.maximum)
            : MCS.validateMapping(inputs.query.graph, inputs.target.graph,
                  match.atomMapping(), row.chemistry);
        if (!errors.isEmpty()) throw new IllegalStateException("Invalid mapping: " + errors);
      }
  }

  private static List<Case> load(Path input) throws Exception {
    List<Case> cases = new ArrayList<>();
    String header = "case_id\toperation\tobjective\tconnectivity\tprofile\tinduced\ttimeout_ms\tlimit\tbatch_size\tquery\ttarget\tweights\tphases";
    boolean first = true;
    Set<String> identifiers = new HashSet<>();
    for (String line : Files.readAllLines(input, StandardCharsets.UTF_8)) {
      if (line.trim().isEmpty() || line.startsWith("#")) continue;
      if (first) {
        if (!line.equals(header)) throw new IllegalArgumentException("Unexpected TSV header");
        first = false;
        continue;
      }
      String[] columns = line.split("\t", -1);
      if (columns.length != 13) throw new IllegalArgumentException("Expected 13 columns: " + line);
      try {
        Case row = new Case(columns);
        if (!identifiers.add(row.id)) throw new IllegalArgumentException("Duplicate case ID");
        cases.add(row);
      } catch (RuntimeException failure) {
        throw new IllegalArgumentException("Case " + columns[0] + ": " + failure.getMessage(), failure);
      }
    }
    if (first) throw new IllegalArgumentException("Missing TSV header");
    return cases;
  }

  public static void main(String[] arguments) throws Exception {
    Map<String, String> flags = new LinkedHashMap<>();
    Set<String> accepted = new HashSet<>(
        Arrays.asList("--cases", "--output", "--warmup", "--repeats", "--phase", "--case"));
    for (int argument = 0; argument < arguments.length; argument += 2) {
      if (argument + 1 >= arguments.length || !accepted.contains(arguments[argument])
          || flags.put(arguments[argument], arguments[argument + 1]) != null)
        throw new IllegalArgumentException("Expected --cases FILE --output FILE --warmup N --repeats N [--phase PHASE] [--case ID]");
    }
    if (!flags.keySet().containsAll(Arrays.asList("--cases", "--output", "--warmup", "--repeats")))
      throw new IllegalArgumentException("Missing required arguments");
    int warmup = Integer.parseInt(flags.get("--warmup"));
    int repeats = Integer.parseInt(flags.get("--repeats"));
    if (warmup < 0 || repeats < 1) throw new IllegalArgumentException("Invalid warmup/repeat count");
    if (flags.containsKey("--phase") && !PHASES.contains(flags.get("--phase")))
      throw new IllegalArgumentException("Unknown phase " + flags.get("--phase"));
    Path output = Paths.get(flags.get("--output"));
    if (output.getParent() != null) Files.createDirectories(output.getParent());
    int selected = 0;
    try (BufferedWriter writer = Files.newBufferedWriter(output, StandardCharsets.UTF_8)) {
      for (Case row : load(Paths.get(flags.get("--cases")))) {
        if (flags.containsKey("--case") && !row.id.equals(flags.get("--case"))) continue;
        List<String> phases = new ArrayList<>();
        for (String phase : row.phases)
          if (!flags.containsKey("--phase") || phase.equals(flags.get("--phase"))) phases.add(phase);
        if (phases.isEmpty()) continue;
        selected++;
        try {
          Inputs inputs = new Inputs(row);
          String graphPhase = phases.contains("prepared") ? "prepared" : phases.get(0);
          for (String phase : phases) {
            Prepared query = phase.equals("prepare") || phase.equals("end_to_end")
                ? null : new Prepared(row, inputs);
            for (int iteration = 0; iteration < warmup; iteration++) {
              Measurement ignored = measure(row, phase, inputs, query);
              retained = ignored;
            }
            for (int sample = 0; sample < repeats; sample++) {
              Measurement measured = measure(row, phase, inputs, query);
              retained = measured;
              validate(row, inputs, measured.results);
              Map<String, Object> line = new LinkedHashMap<>();
              line.put("case_id", row.id);
              line.put("language", "java");
              line.put("phase", phase);
              line.put("sample", sample);
              line.put("elapsed_ns", measured.elapsed);
              List<Map<String, Object>> results = new ArrayList<>();
              for (SearchResult result : measured.results) results.add(result(result));
              line.put("results", results);
              line.put("query_graph", sample == 0 && phase.equals(graphPhase) ? graph(inputs.query) : null);
              line.put("target_graph", sample == 0 && phase.equals(graphPhase) ? graph(inputs.target) : null);
              writer.write(JSON.writeValueAsString(line));
              writer.newLine();
            }
            writer.flush();
          }
        } catch (Exception failure) {
          throw new IllegalStateException("Case " + row.id + ": " + failure.getMessage(), failure);
        }
      }
    }
    if (selected == 0) throw new IllegalArgumentException("No selected case/phase");
  }
}

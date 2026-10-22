/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.MCS;
import com.bioinception.smsd.core.MolGraph;
import com.bioinception.smsd.core.Substructure;
import com.bioinception.smsd.core.SearchResult;
import com.bioinception.smsd.core.Standardiser;
import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.databind.ObjectMapper;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.stream.IntStream;
import java.io.InputStream;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IDoubleBondStereochemistry;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.stereo.DoubleBondStereochemistry;
import uk.ac.ebi.beam.Bond;
import uk.ac.ebi.beam.Edge;
import uk.ac.ebi.beam.Graph;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

/**
 * Picocli-based CLI for substructure and MCS search with JSON output, SDF batch processing,
 * Tanimoto similarity, and MCS fragment export.
 */
@Command(
    name = "smsd",
    mixinStandardHelpOptions = true,
    exitCodeOnExecutionException = 2,
    footer = "Exit status: 0 = match, 1 = no match, 2 = invalid input, error or incomplete search.",
    versionProvider = SMSDcli.ManifestVersion.class,
    description =
        "SMSD — Substructure & MCS search engine by BioInception PVT LTD"
            + " (com.bioinceptionlabs).%nSupports substructure matching, MCS computation, SMARTS"
            + " queries, SDF batch processing,%nTanimoto similarity, and JSON/SMILES/MDL mapping"
            + " export.")
public class SMSDcli implements Callable<Integer> {

  /** The {@code --version} text, with the release version read from the JAR manifest. */
  public static final class ManifestVersion implements CommandLine.IVersionProvider {
    @Override
    public String[] getVersion() {
      Package cli = SMSDcli.class.getPackage();
      return new String[] {versionLine(cli == null ? null : cli.getImplementationVersion())};
    }
  }

  /** The version line for a manifest Implementation-Version, which is absent outside a JAR. */
  static String versionLine(String implementationVersion) {
    String version = implementationVersion == null || implementationVersion.trim().isEmpty()
        ? "(development build)" : implementationVersion.trim();
    return "SMSD Pro " + version + " — com.bioinceptionlabs";
  }

  @Option(names = {"--Q"}, required = true,
      description = "Query type: SMI, SIG (SMARTS), MOL, ML2, PDB, CML.")
  String qType;

  @Option(names = {"--q"}, required = true,
      description = "Query SMILES/SMARTS string or file path (by type).")
  String q;

  @Option(names = {"--T"}, required = true,
      description = "Target type: SMI, MOL, ML2, PDB, CML, SDF (batch).")
  String tType;

  @Option(names = {"--t"}, required = true,
      description = "Target SMILES string or file path (by type).")
  String t;

  @Option(names = {"--mode"}, description = "Mode: sub (substructure) or mcs. Default: sub",
      defaultValue = "sub")
  String mode;

  @Option(names = {"--objective"}, defaultValue = "MCIS", description = "MCS objective: MCIS or MCES.")
  MCS.Objective objective = MCS.Objective.MCIS;

  @Option(names = {"--fragments"}, description = "Allow disconnected common fragments in MCS mode.")
  boolean fragments;

  @Option(names = {"--induced"}, description = "Require an induced substructure embedding.")
  boolean induced;

  @Option(names = {"-m", "--mappings"},
      description = "Enumerate substructure atom mappings.")
  boolean allMappings;

  @Option(names = {"--timeout"}, description = "Time limit in ms (default: 10000).")
  long timeoutMs = 10_000L;

  @Option(names = {"-a", "--add-hydrogens"},
      description = "Add explicit hydrogens before matching.")
  boolean addHydrogens;

  @Option(names = {"-r", "--remove-hydrogens"},
      description = "Remove explicit hydrogens before matching.")
  boolean removeHydrogens;

  @Option(names = {"--json"}, description = "Write JSON to file or '-' for stdout.")
  String jsonOut;

  @Option(names = {"--json-pretty"}, description = "Pretty-print JSON.")
  boolean jsonPretty;

  @Option(names = {"--map-out"},
      description = "Write mapping export to file or '-' for stdout.")
  String mapOut;

  @Option(names = {"--map-format"},
      description = "Mapping format: json|smiles|smarts|mdl. Default: json",
      defaultValue = "json")
  String mapFormat;

  @Option(names = {"--map-pretty"}, description = "Pretty-print JSON mapping export.")
  boolean mapPretty;

  @Option(names = {"-o", "--output"},
      description = "Write MCS fragment as SMILES to file or '-' for stdout (MCS mode only).")
  String mcsFragmentOut;

  @Option(names = {"--benchmark"},
      description = "Run a quick speed test: substructure + MCS on the given query/target pair.")
  boolean benchmark;

  @Option(names = {"--threads"},
      description = "Worker threads for SDF batch mode. 0 = auto (all available processors). Default: 0.")
  int numThreads = 0;

  private static final ObjectMapper MAPPER = new ObjectMapper();

  public static void main(String[] args) {
    System.exit(new CommandLine(new SMSDcli()).execute(args));
  }

  @Override
  public Integer call() throws Exception {
    if (!"sub".equalsIgnoreCase(mode) && !"mcs".equalsIgnoreCase(mode))
      throw new CommandLine.ParameterException(new CommandLine(this), "Unknown mode: " + mode);
    if ("mcs".equalsIgnoreCase(mode) && "SIG".equalsIgnoreCase(qType))
      throw new CommandLine.ParameterException(new CommandLine(this), "SMARTS queries require substructure mode");
    if (benchmark) return runBenchmark();
    if ("SDF".equalsIgnoreCase(tType)) return runBatch();
    return runSingle();
  }

  private int runSingle() throws Exception {
    MolIO.Query query = MolIO.loadQuery(qType, q);
    IAtomContainer tm = applyHydrogenOptions(MolIO.loadTarget(tType, t));
    tm = prepare(tm);
    IAtomContainer qm = query.isSmarts() ? null : prepare(applyHydrogenOptions(query.container()));
    Substructure substructure = createSubstructure(query, qm, new ChemOptions());

    Map<String, Object> result = new LinkedHashMap<>();
    result.put("query", q);
    result.put("target", t);
    result.put("mode", mode);
    result.put("timeout_ms", timeoutMs);

    return "sub".equalsIgnoreCase(mode)
        ? runSubstructure(substructure, qm, tm, result)
        : runMCS(qm, tm, result);
  }

  private int runSubstructure(
      Substructure query, IAtomContainer qm, IAtomContainer tm, Map<String, Object> result)
      throws Exception {
    if (!allMappings) {
      SearchResult search = query.withTimeout(timeoutMs).match(tm);
      result.put("exists", search.found());
      recordOutcome(result, search);
      writeJson(result);
      return exitCode(search);
    }
    SearchResult search = query.withTimeout(timeoutMs).matchAll(tm, 10_000);
    recordOutcome(result, search);
    List<Map<Integer, Integer>> maps = mappings(search);
    List<Map<String, Object>> items = new ArrayList<>();
    int idx = 0;
    for (Map<Integer, Integer> m : maps) {
      Map<String, Object> row = new LinkedHashMap<>();
      row.put("index", idx++);
      row.put("pairs", toPairs(m));
      row.put("query_atoms", new ArrayList<>(m.keySet()));
      row.put("target_atoms", new ArrayList<>(m.values()));
      items.add(row);
    }
    result.put("mappings", items);
    if (qm != null && tm != null) {
      result.put("tanimoto", tanimoto(qm, tm, maps.isEmpty() ? 0 : maps.get(0).size()));
    }
    if (mapOut != null) writeMappingExport(qm, tm, maps);
    writeJson(result);
    return exitCode(search);
  }

  private int runMCS(IAtomContainer qm, IAtomContainer tm, Map<String, Object> result)
      throws Exception {
    if (qm == null) throw new IllegalArgumentException("SMARTS queries require substructure mode");
    SearchResult search = MCS.find(qm).withOptions(mcsOptions(new ChemOptions())).withTimeout(timeoutMs).match(tm);
    Map<Integer, Integer> m = firstMapping(search);
    recordOutcome(result, search);
    result.put("objective", objective.toString());
    result.put("mcs_size", m.size());
    result.put("pairs", toPairs(m));
    result.put("query_atoms", new ArrayList<>(m.keySet()));
    result.put("target_atoms", new ArrayList<>(m.values()));

    String mcsSmi = search.found() ? search.matches().get(0).toSmiles(qm) : null;
    if (qm != null && tm != null) {
      result.put("tanimoto", mcsSimilarity(qm, tm, search));
      result.put("similarity_upper_bound", similarityUpperBound(qm, tm));
      if (mcsSmi != null) result.put("mcs_smiles", mcsSmi);
    }
    if (mapOut != null) writeMCSMappingExport(qm, tm, search);
    writeJson(result);

    if (mcsFragmentOut != null && mcsSmi != null) {
      if ("-".equals(mcsFragmentOut)) {
        System.out.println(mcsSmi);
      } else {
        java.nio.file.Files.write(java.nio.file.Paths.get(mcsFragmentOut),
            (mcsSmi + "\n").getBytes(StandardCharsets.UTF_8));
      }
    }
    return exitCode(search);
  }

  private static int exitCode(SearchResult search) {
    return search.complete() ? (search.found() ? 0 : 1) : 2;
  }

  @SuppressWarnings("unchecked")
  private int runBatch() throws Exception {
    MolIO.Query query = MolIO.loadQuery(qType, q);
    List<IAtomContainer> targets = MolIO.loadTargetsSDF(t);
    ChemOptions chem = new ChemOptions();
    int n = targets.size();

    // Apply hydrogen options to the query ONCE here, outside the parallel block.
    // CDK IAtomContainer is not thread-safe for concurrent modification; calling
    // applyHydrogenOptions(query.container()) from N threads simultaneously
    // would corrupt it. Snapshot the result into a final local.
    final IAtomContainer queryMol = query.isSmarts() ? null : prepare(applyHydrogenOptions(query.container()));
    final Substructure substructure = createSubstructure(query, queryMol, chem).withTimeout(timeoutMs);
    final MCS maximum = queryMol == null ? null : MCS.find(queryMol).withOptions(mcsOptions(chem)).withTimeout(timeoutMs);

    // Resolve thread count: 0 = all available processors
    int cores = (numThreads > 0) ? numThreads : Runtime.getRuntime().availableProcessors();

    // Pre-allocate result slots so parallel writes stay at distinct indices
    Map<String, Object>[] results = new Map[n];
    int[] exitCodes = new int[n];
    Arrays.fill(exitCodes, 2);

    ForkJoinPool pool = new ForkJoinPool(cores);
    try {
      List<Callable<Void>> tasks = new ArrayList<>(n);
      for (int i = 0; i < n; i++) {
        final int idx = i;
        tasks.add(() -> {
          results[idx] = batchRow(idx, targets.get(idx), substructure, maximum, queryMol, exitCodes);
          return null;
        });
      }
      pool.invokeAll(tasks);
    } catch (Exception e) {
      // Every row left unwritten by the pool is searched in turn, so the output is complete.
      if (e instanceof InterruptedException) Thread.currentThread().interrupt();
      for (int i = 0; i < n; i++)
        if (results[i] == null) results[i] = batchRow(i, targets.get(i), substructure, maximum, queryMol, exitCodes);
    } finally {
      pool.shutdown();
    }

    Map<String, Object> output = new LinkedHashMap<>();
    output.put("query", q);
    output.put("target_file", t);
    output.put("target_count", n);
    output.put("threads_used", cores);
    output.put("mode", mode);
    output.put("results", Arrays.asList(results));
    writeJson(output);
    boolean found = false;
    for (int code : exitCodes) {
      if (code == 2) return 2;
      found |= code == 0;
    }
    return found ? 0 : 1;
  }

  /** The JSON row of one batch target: its search outcome, or the error that stopped it. */
  private Map<String, Object> batchRow(int index, IAtomContainer target, Substructure substructure, MCS maximum,
      IAtomContainer queryMol, int[] exitCodes) {
    Map<String, Object> row = new LinkedHashMap<>();
    row.put("target_index", index);
    try {
      IAtomContainer tm = prepare(applyHydrogenOptions(target));
      if ("sub".equalsIgnoreCase(mode)) {
        SearchResult search = substructure.match(tm);
        exitCodes[index] = exitCode(search);
        recordOutcome(row, search);
        row.put("exists", search.found());
      } else {
        if (maximum == null) throw new IllegalArgumentException("SMARTS queries require substructure mode");
        SearchResult search = maximum.match(tm);
        exitCodes[index] = exitCode(search);
        recordOutcome(row, search);
        Map<Integer, Integer> m = firstMapping(search);
        row.put("mcs_size", m.size());
        row.put("pairs", toPairs(m));
        if (queryMol != null && tm != null)
          row.put("tanimoto", mcsSimilarity(queryMol, tm, search));
      }
    } catch (Exception e) {
      row.put("error", e.getMessage());
    }
    return row;
  }

  private int runBenchmark() throws Exception {
    MolIO.Query query = MolIO.loadQuery(qType, q);
    IAtomContainer target = prepare(applyHydrogenOptions(MolIO.loadTarget(tType, t)));
    IAtomContainer queryMolecule = query.isSmarts() ? null : prepare(applyHydrogenOptions(query.container()));
    Substructure substructure = createSubstructure(query, queryMolecule, new ChemOptions()).withTimeout(timeoutMs);
    MCS maximum = queryMolecule == null ? null : MCS.find(queryMolecule)
        .withOptions(mcsOptions(new ChemOptions())).withTimeout(timeoutMs);
    int warmupRuns = 3, timedRuns = 10;
    PrintWriter writer = new PrintWriter(System.out, true);
    writer.println("SMSD Benchmark (com.bioinceptionlabs)");
    writer.printf("Query: %s%nTarget: %s%nWarmup: %d runs, Timed: %d runs%n", q, t, warmupRuns, timedRuns);
    for (int run = 0; run < warmupRuns; run++) {
      substructure.match(target);
      if (maximum != null) maximum.match(target);
    }
    long substructureNanos = 0;
    SearchResult result = null;
    boolean foundMatch = false, incomplete = false;
    for (int run = 0; run < timedRuns; run++) {
      long started = System.nanoTime();
      result = substructure.match(target);
      foundMatch |= result.found();
      incomplete |= !result.complete();
      substructureNanos += System.nanoTime() - started;
    }
    writer.printf("Substructure: %s, %s, average %.3f ms%n", result.found() ? "FOUND" : "NOT FOUND",
        result.status(), substructureNanos / (double) timedRuns / 1_000_000);
    if (maximum != null) {
      long mcsNanos = 0;
      for (int run = 0; run < timedRuns; run++) {
        long started = System.nanoTime();
        result = maximum.match(target);
        foundMatch |= result.found();
        incomplete |= !result.complete();
        mcsNanos += System.nanoTime() - started;
      }
      writer.printf("MCS: %s, score %.3f, %s, average %.3f ms%n", objective, result.score(),
          result.status(), mcsNanos / (double) timedRuns / 1_000_000);
    }
    return incomplete ? 2 : foundMatch ? 0 : 1;
  }

  private Substructure createSubstructure(MolIO.Query query, IAtomContainer molecule, ChemOptions chemistry) {
    Substructure.Options options = new Substructure.Options();
    options.chemistry = chemistry;
    options.induced = induced;
    return (query.isSmarts() ? Substructure.smarts(query.text()) : Substructure.find(molecule)).withOptions(options);
  }

  private MCS.Options mcsOptions(ChemOptions chemistry) {
    MCS.Options options = new MCS.Options();
    options.chemistry = chemistry;
    options.objective = objective;
    options.connectivity = fragments ? MCS.Connectivity.Fragments : MCS.Connectivity.Connected;
    return options;
  }

  private static IAtomContainer prepare(IAtomContainer molecule) {
    try { return Standardiser.standardise(molecule, Standardiser.TautomerMode.NONE); }
    catch (Exception failure) { return molecule; }
  }

  private static Map<Integer, Integer> firstMapping(SearchResult result) {
    return result.found() ? result.matches().get(0).atomMapping() : java.util.Collections.emptyMap();
  }

  private static List<Map<Integer, Integer>> mappings(SearchResult result) {
    List<Map<Integer, Integer>> mappings = new ArrayList<>();
    for (SearchResult.Match match : result.matches()) mappings.add(match.atomMapping());
    return mappings;
  }

  private static void recordOutcome(Map<String, Object> output, SearchResult result) {
    output.put("status", result.status().toString());
    output.put("optimal", result.optimal());
    output.put("complete", result.complete());
    output.put("score", result.score());
    output.put("nodes_visited", result.nodesVisited());
    output.put("elapsed_ms", result.elapsedMillis());
    if (result.found()) {
      SearchResult.Match match = result.matches().get(0);
      List<Map<String, Object>> bonds = new ArrayList<>();
      for (Map.Entry<SearchResult.Bond, SearchResult.Bond> bond : match.bondMapping().entrySet()) {
        Map<String, Object> pair = new LinkedHashMap<>();
        pair.put("query", Arrays.asList(bond.getKey().first(), bond.getKey().second()));
        pair.put("target", Arrays.asList(bond.getValue().first(), bond.getValue().second()));
        bonds.add(pair);
      }
      output.put("bonds", bonds);
      output.put("fragments", match.fragments());
    }
  }

  private double mcsSimilarity(IAtomContainer query, IAtomContainer target, SearchResult result) {
    if (objective == MCS.Objective.MCIS) return tanimoto(query, target, firstMapping(result).size());
    int shared = result.found() ? result.matches().get(0).bondCount() : 0;
    int total = query.getBondCount() + target.getBondCount() - shared;
    return total > 0 ? (double) shared / total : 0;
  }

  /** The bond-count bound for MCES, or the atom-label bound of the search's chemistry for MCIS. */
  private double similarityUpperBound(IAtomContainer query, IAtomContainer target) {
    if (objective == MCS.Objective.MCES) {
      int shared = Math.min(query.getBondCount(), target.getBondCount());
      int total = query.getBondCount() + target.getBondCount() - shared;
      return total > 0 ? (double) shared / total : 0;
    }
    return MCS.similarityUpperBound(query, target, new ChemOptions());
  }

  private IAtomContainer applyHydrogenOptions(IAtomContainer mol) throws Exception {
    if (mol == null) return null;
    if (addHydrogens) AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
    if (removeHydrogens) mol = AtomContainerManipulator.removeHydrogens(mol);
    return mol;
  }

  private static double tanimoto(IAtomContainer query, IAtomContainer target, int commonAtoms) {
    int denom = query.getAtomCount() + target.getAtomCount() - commonAtoms;
    return denom > 0 ? (double) commonAtoms / denom : 0.0;
  }

  static IAtomContainer buildSubgraph(IAtomContainer target, Set<Integer> atomIdxs)
      throws CloneNotSupportedException {
    IAtomContainer copy = target.clone();
    Set<IAtom> selected = new HashSet<>();
    for (int atom : atomIdxs) selected.add(copy.getAtom(atom));
    return AtomContainerManipulator.extractSubstructure(copy, selected);
  }

  private void writeMappingExport(
      IAtomContainer qm, IAtomContainer tm, List<Map<Integer, Integer>> maps) throws Exception {
    boolean isStdout = "-".equals(mapOut);
    OutputStream os = isStdout ? System.out : new FileOutputStream(mapOut);
    try {
      OutputUtil.writeMappings(qm, tm, maps, parseOutType(mapFormat), os, mapPretty);
    } finally {
      if (!isStdout) os.close();
    }
  }

  private void writeMCSMappingExport(
      IAtomContainer query, IAtomContainer target, SearchResult search) throws Exception {
    if (parseOutType(mapFormat) != OutputUtil.OutType.MOL) {
      writeMappingExport(query, target, mappings(search));
      return;
    }
    boolean isStdout = "-".equals(mapOut);
    OutputStream output = isStdout ? System.out : new FileOutputStream(mapOut);
    try {
      if (search.found())
        OutputUtil.writeMolecule(search.matches().get(0).toMolecule(query), OutputUtil.OutType.MOL, output);
    } finally {
      if (!isStdout) output.close();
    }
  }

  private static List<List<Integer>> toPairs(Map<Integer, Integer> m) {
    List<List<Integer>> out = new ArrayList<>(m.size());
    for (Map.Entry<Integer, Integer> e : m.entrySet())
      out.add(Arrays.asList(e.getKey(), e.getValue()));
    return out;
  }

  private void writeJson(Map<String, Object> data) throws Exception {
    if (jsonOut == null) return;
    String payload = jsonPretty
        ? MAPPER.writerWithDefaultPrettyPrinter().writeValueAsString(data)
        : MAPPER.writeValueAsString(data);
    if ("-".equals(jsonOut)) {
      PrintWriter pw = new PrintWriter(new OutputStreamWriter(System.out, StandardCharsets.UTF_8));
      pw.println(payload);
      pw.flush();
    } else {
      java.nio.file.Files.write(java.nio.file.Paths.get(jsonOut), payload.getBytes(StandardCharsets.UTF_8));
    }
  }

  private static OutputUtil.OutType parseOutType(String s) {
    switch ((s == null ? "json" : s).trim().toLowerCase(Locale.ROOT)) {
      case "json": return OutputUtil.OutType.JSON;
      case "smiles":
      case "smi": return OutputUtil.OutType.SMI;
      case "smarts": return OutputUtil.OutType.SMARTS;
      case "mdl":
      case "mol":
      case "sdf": return OutputUtil.OutType.MOL;
      default: throw new IllegalArgumentException("Unknown --map-format: " + s);
    }
  }

  // ========================================================================
  // Static inner class: MolIO
  // ========================================================================

  /** Molecule I/O utilities for loading molecules from SMILES, SMARTS, and file formats. */
  public static final class MolIO {

    private MolIO() {}

    private static final ThreadLocal<SmilesParser> SP = ThreadLocal.withInitial(
        () -> new SmilesParser(SilentChemObjectBuilder.getInstance()));

    /** Result of loading a query: either a SMARTS string or a parsed molecule. */
    public static final class Query {
      private final boolean isSmarts;
      private final String text;
      private final IAtomContainer container;

      @JsonCreator
      public Query(@JsonProperty("isSmarts") boolean isSmarts,
          @JsonProperty("text") String text, @JsonProperty("container") IAtomContainer container) {
        this.isSmarts = isSmarts;
        this.text = text;
        this.container = container;
      }

      @JsonProperty("isSmarts")
      public boolean isSmarts() { return isSmarts; }

      @JsonProperty("text")
      public String text() { return text; }

      @JsonProperty("container")
      public IAtomContainer container() { return container; }

      @Override
      public final boolean equals(Object other) {
        if (this == other) return true;
        if (!(other instanceof Query)) return false;
        Query query = (Query) other;
        return isSmarts == query.isSmarts && Objects.equals(text, query.text)
            && Objects.equals(container, query.container);
      }

      @Override
      public final int hashCode() {
        return 31 * (31 * Boolean.hashCode(isSmarts) + Objects.hashCode(text)) + Objects.hashCode(container);
      }

      @Override
      public final String toString() {
        return "Query[isSmarts=" + isSmarts + ", text=" + text + ", container=" + container + "]";
      }
    }

    /**
     * Load a query molecule or SMARTS pattern.
     * File input must be non-empty; CML and PDB input must contain exactly one
     * molecule/model. Use SDF target input for batch searches.
     *
     * @param type input type: SMI (SMILES), SIG (SMARTS), MOL, ML2, PDB, or CML
     * @param value SMILES/SMARTS string or file path
     */
    public static Query loadQuery(String type, String value) throws Exception {
      Objects.requireNonNull(type, "Query type must not be null");
      Objects.requireNonNull(value, "Query value must not be null");
      switch (type.toUpperCase(Locale.ROOT)) {
        case "SIG": return new Query(true, value, null);
        case "SMI": return new Query(false, value, readSmiles(value));
        default: return new Query(false, value, readFile(type, value));
      }
    }

    /**
     * Load a single target molecule.
     * CML and PDB input must contain exactly one non-empty molecule/model.
     * Files with multiple molecules are rejected; SDF supports batch targets.
     *
     * @param type input type: SMI (SMILES), MOL, ML2, PDB, or CML
     * @param value SMILES string or file path
     */
    public static IAtomContainer loadTarget(String type, String value) throws Exception {
      Objects.requireNonNull(type, "Target type must not be null");
      Objects.requireNonNull(value, "Target value must not be null");
      switch (type.toUpperCase(Locale.ROOT)) {
        case "SMI": return readSmiles(value);
        case "SDF": throw new IllegalArgumentException("Use loadTargetsSDF for multi-molecule SDF.");
        default: return readFile(type, value);
      }
    }

    /**
     * The SMILES with each atom written in capitals on a ring bond written ':' written aromatic, as
     * the native reader reads it, so that CDK gives it the hydrogens and Kekulé structure of an
     * aromatic atom while an atom in brackets keeps the hydrogens written for it. An atom with a
     * double or triple ring bond as well keeps its symbol. The text is unchanged when CDK cannot
     * read it, or its aromatic form.
     */
    static String aromaticBondAtoms(String value) {
      if (value.indexOf(':') < 0) return value;
      IAtomContainer raw;
      try {
        raw = SP.get().parseSmiles(value);
        new MolGraph(raw);  // marks the ring atoms and bonds of raw
      } catch (CDKException | RuntimeException unreadable) {
        return value;
      }
      boolean[] aromaticRing = new boolean[raw.getAtomCount()], multipleRing = new boolean[raw.getAtomCount()];
      for (IBond bond : raw.bonds()) {
        if (!bond.isInRing()) continue;
        boolean[] marks = bond.isAromatic() ? aromaticRing
            : bond.getOrder() != null && bond.getOrder().numeric() > 1 ? multipleRing : null;
        if (marks != null) marks[raw.indexOf(bond.getBegin())] = marks[raw.indexOf(bond.getEnd())] = true;
      }
      StringBuilder text = new StringBuilder(value);
      boolean changed = false;
      int atom = 0;
      for (int at = 0; at < text.length() && !Character.isWhitespace(text.charAt(at)); at++) {
        char c = text.charAt(at);
        boolean rewrite = atom < aromaticRing.length && aromaticRing[atom] && !multipleRing[atom];
        if (c == '[') {
          int close = text.indexOf("]", at);
          if (close < 0) return value;
          int symbol = at + 1;
          while (symbol < close && Character.isDigit(text.charAt(symbol))) symbol++;
          int end = symbol < close && Character.isUpperCase(text.charAt(symbol)) ? symbol + 1 : symbol;
          if (end < close && Character.isLowerCase(text.charAt(end))) end++;
          if (rewrite && Arrays.asList("B", "C", "N", "O", "P", "S", "Se", "As").contains(text.substring(symbol, end))) {
            text.setCharAt(symbol, Character.toLowerCase(text.charAt(symbol)));
            changed = true;
          }
          atom++;
          at = close;
        } else if (c == '%') {
          at = at + 1 < text.length() && text.charAt(at + 1) == '(' ? text.indexOf(")", at) : at + 2;
          if (at < 0) return value;
        } else if (Character.isLetter(c) || c == '*') {
          boolean halogen = at + 1 < text.length() && (c == 'C' && text.charAt(at + 1) == 'l' || c == 'B' && text.charAt(at + 1) == 'r');
          if (rewrite && !halogen && "BCNOPS".indexOf(c) >= 0) {
            text.setCharAt(at, Character.toLowerCase(c));
            changed = true;
          }
          if (halogen) at++;
          atom++;
        }
      }
      if (!changed || atom != raw.getAtomCount()) return value;
      try {
        SP.get().parseSmiles(text.toString());
        return text.toString();
      } catch (CDKException noKekuleStructure) {
        return value;
      }
    }

    /** Retain aromatic directional carriers omitted by the CDK conversion. */
    private static IAtomContainer readSmiles(String value) throws Exception {
      value = aromaticBondAtoms(value);
      IAtomContainer molecule;
      Graph graph = null;
      try {
        molecule = SP.get().parseSmiles(value);
      } catch (InvalidSmilesException rejected) {
        if (value.indexOf('/') < 0 && value.indexOf('\\') < 0) throw rejected;
        try {
          graph = Graph.parse(value, false, new HashSet<String>());
        } catch (IOException invalid) {
          throw rejected;
        }
        if (!hasAromaticStereoFocus(graph)) throw rejected;
        // The CDK localiser honours explicit aromatic double bonds and their carriers.
        SmilesParser unlocalised = new SmilesParser(SilentChemObjectBuilder.getInstance());
        unlocalised.kekulise(false);
        try {
          molecule = unlocalised.parseSmiles(value);
          Kekulization.kekulize(molecule);
        } catch (CDKException invalid) {
          throw rejected;
        }
      }
      if (value.indexOf('/') < 0 && value.indexOf('\\') < 0) return molecule;
      boolean hasAromaticBond = false;
      for (IBond bond : molecule.bonds())
        if (bond.isAromatic()) {
          hasAromaticBond = true;
          break;
        }
      if (!hasAromaticBond) return molecule;
      if (graph == null)
        graph = Graph.parse(value, false, new HashSet<String>()).kekule();
      Set<IBond> existing = new HashSet<>();
      for (IStereoElement<?, ?> element : molecule.stereoElements())
        if (element instanceof IDoubleBondStereochemistry)
          existing.add(((IDoubleBondStereochemistry) element).getStereoBond());
      for (Edge edge : graph.edges()) {
        if (edge.bond() != Bond.DOUBLE) continue;
        int first = edge.either(), second = edge.other(first);
        IBond focus = molecule.getBond(molecule.getAtom(first), molecule.getAtom(second));
        if (existing.contains(focus)) continue;
        Edge firstCarrier = smilesCarrier(graph, first, second);
        Edge secondCarrier = smilesCarrier(graph, second, first);
        if (firstCarrier == null || secondCarrier == null) continue;
        Bond firstDirection = firstCarrier.bond(first);
        Bond secondDirection = secondCarrier.bond(second);
        if (!aromaticSmilesCarrier(graph, firstCarrier)
            && !aromaticSmilesCarrier(graph, secondCarrier)) continue;
        IBond firstReference = molecule.getBond(
            molecule.getAtom(first), molecule.getAtom(firstCarrier.other(first)));
        IBond secondReference = molecule.getBond(
            molecule.getAtom(second), molecule.getAtom(secondCarrier.other(second)));
        molecule.addStereoElement(new DoubleBondStereochemistry(
            focus, new IBond[] {firstReference, secondReference},
            smilesDirection(firstDirection) == smilesDirection(secondDirection)
                ? IDoubleBondStereochemistry.Conformation.TOGETHER
                : IDoubleBondStereochemistry.Conformation.OPPOSITE));
      }
      return molecule;
    }

    private static boolean hasAromaticStereoFocus(Graph graph) {
      for (Edge edge : graph.edges()) {
        if (edge.bond() != Bond.DOUBLE) continue;
        int first = edge.either(), second = edge.other(first);
        if (graph.atom(first).aromatic() && graph.atom(second).aromatic()
            && smilesCarrier(graph, first, second) != null
            && smilesCarrier(graph, second, first) != null) return true;
      }
      return false;
    }

    private static boolean aromaticSmilesCarrier(
        Graph graph, Edge carrier) {
      int first = carrier.either();
      return carrier.bond() == Bond.UP_AROMATIC
          || carrier.bond() == Bond.DOWN_AROMATIC
          || (graph.atom(first).aromatic() && graph.atom(carrier.other(first)).aromatic());
    }

    private static Edge smilesCarrier(
        Graph graph, int atom, int otherEnd) {
      Edge carrier = null;
      int direction = 0;
      for (Edge edge : graph.edges(atom)) {
        if (edge.other(atom) == otherEnd) continue;
        int current = smilesDirection(edge.bond(atom));
        if (current == 0) continue;
        if (carrier == null) {
          carrier = edge;
          direction = current;
        } else if (current == direction) return null;
      }
      return carrier;
    }

    private static int smilesDirection(Bond bond) {
      switch (bond) {
        case UP:
        case UP_AROMATIC: return 1;
        case DOWN:
        case DOWN_AROMATIC: return -1;
        default: return 0;
      }
    }

    /** Maximum number of molecules to load from a single SDF file. */
    private static final int SDF_BATCH_LIMIT = 100_000;

    /** Load all molecules from an SDF file for batch processing. */
    public static List<IAtomContainer> loadTargetsSDF(String path) throws Exception {
      Objects.requireNonNull(path, "SDF path must not be null");
      List<IAtomContainer> out = new ArrayList<>();
      try (FileInputStream fis = new FileInputStream(path);
          IteratingSDFReader it =
              new IteratingSDFReader(fis, SilentChemObjectBuilder.getInstance())) {
        while (it.hasNext()) {
          if (out.size() >= SDF_BATCH_LIMIT)
            throw new IllegalStateException("SDF batch limit exceeded (" + SDF_BATCH_LIMIT + " molecules)");
          IAtomContainer molecule = it.next();
          correctMdlRadicalHydrogens(molecule);
          out.add(molecule);
        }
      }
      return out;
    }

    /** Apply explicit MOL valence to freshly read radical atoms. */
    private static void correctMdlRadicalHydrogens(IAtomContainer molecule) {
      for (IAtom atom : molecule.atoms()) {
        Integer hydrogens = atom.getImplicitHydrogenCount();
        Integer valence = atom.getValency();
        int unpaired = molecule.getConnectedSingleElectronsCount(atom);
        if (hydrogens == null || hydrogens != 0 || valence == null || unpaired == 0) continue;
        int bondValence = 0;
        for (IBond bond : molecule.getConnectedBondsList(atom)) {
          IBond.Order order = bond.getOrder();
          if (order == null || order == IBond.Order.UNSET) {
            bondValence = -1;
            break;
          }
          bondValence += order.numeric();
        }
        if (bondValence >= 0 && valence > bondValence && valence < (long) bondValence + unpaired)
          atom.setImplicitHydrogenCount(valence - bondValence);
      }
    }

    private static IAtomContainer readFile(String type, String path) throws Exception {
      Function<InputStream, ISimpleChemObjectReader> factory;
      switch (type.toUpperCase(Locale.ROOT)) {
        case "MOL": factory = MDLV2000Reader::new; break;
        case "ML2": factory = Mol2Reader::new; break;
        case "PDB": factory = PDBReader::new; break;
        case "CML": factory = CMLReader::new; break;
        default: throw new IllegalArgumentException("Unsupported file type: " + type);
      }
      try (FileInputStream fis = new FileInputStream(path);
          ISimpleChemObjectReader reader = factory.apply(fis)) {
        IAtomContainer molecule;
        if ("CML".equalsIgnoreCase(type) || "PDB".equalsIgnoreCase(type)) {
          IChemFile file = reader.read(
              DefaultChemObjectBuilder.getInstance().newInstance(IChemFile.class));
          List<IAtomContainer> molecules = file == null
              ? Collections.emptyList() : ChemFileManipulator.getAllAtomContainers(file);
          if (molecules.size() != 1) {
            throw new IOException(type + " input must contain exactly one molecule/model; found "
                + molecules.size() + ". Use SDF for batch targets.");
          }
          molecule = molecules.get(0);
        } else {
          molecule = reader.read(
              DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
        }
        if (molecule == null || molecule.getAtomCount() == 0)
          throw new IOException(type + " input does not contain a non-empty molecule.");
        if ("MOL".equalsIgnoreCase(type)) correctMdlRadicalHydrogens(molecule);
        return molecule;
      }
    }
  }

  // ========================================================================
  // Static inner class: OutputUtil
  // ========================================================================

  /** Output helpers for mapping export in JSON, SMILES, SMARTS, and MDL V2000 formats. */
  public static final class OutputUtil {

    private OutputUtil() {}

    public enum OutType { JSON, SMI, SMARTS, MOL }

    /** Write a single molecule as SMILES or MDL V2000 MOL. */
    public static void writeMolecule(IAtomContainer mol, OutType type, OutputStream os)
        throws Exception {
      os=nonClosing(os);
      if (type == OutType.MOL) {
        IAtomContainer output = mol.clone();
        boolean has2D = true, has3D = true;
        for (IAtom atom : output.atoms()) {
          has2D &= atom.getPoint2d() != null;
          has3D &= atom.getPoint3d() != null;
        }
        boolean generated = !has2D && !has3D;
        if (generated)
          new org.openscience.cdk.layout.StructureDiagramGenerator().generateCoordinates(output);
        if (!has3D && (has2D || generated)) prepareMolStereoCoordinates(output, generated);
        try (MDLV2000Writer w = new MDLV2000Writer(os)) {
          if (!has3D)
            w.getSetting(MDLV2000Writer.OptForceWriteAs2DCoordinates).setSetting("true");
          w.write(output);
        }
      } else {
        try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
          w.write(new MolGraph(mol).toCanonicalSmiles());
          w.write(System.lineSeparator());
          w.flush();
        }
      }
    }

    /** Correct generated MOL geometry and verify every supplied double-bond relation. */
    private static void prepareMolStereoCoordinates(IAtomContainer molecule, boolean adjust) {
      List<IDoubleBondStereochemistry> elements = new ArrayList<>();
      for (IStereoElement<?, ?> element : molecule.stereoElements())
        if (element instanceof IDoubleBondStereochemistry)
          elements.add((IDoubleBondStereochemistry) element);
      for (IDoubleBondStereochemistry stereo : elements) {
        if (molStereoGeometry(stereo) == stereo.getStereo()) continue;
        if (!adjust)
          throw new IllegalArgumentException("MOL coordinates conflict with double-bond stereochemistry");
        IAtom first = stereo.getStereoBond().getBegin(), second = stereo.getStereoBond().getEnd();
        Set<IAtom> component = new HashSet<>();
        List<IAtom> pending = new ArrayList<>();
        component.add(second);
        pending.add(second);
        for (int position = 0; position < pending.size(); position++) {
          IAtom atom = pending.get(position);
          for (IBond bond : atom.bonds()) {
            if (bond.equals(stereo.getStereoBond())) continue;
            IAtom neighbor = bond.getOther(atom);
            if (component.add(neighbor)) pending.add(neighbor);
          }
        }
        if (component.contains(first)) {
          pending.clear();
          pending.add(molStereoReferences(stereo)[0]);
        }
        javax.vecmath.Point2d firstPoint = first.getPoint2d(), secondPoint = second.getPoint2d();
        double axisX = secondPoint.x - firstPoint.x, axisY = secondPoint.y - firstPoint.y;
        double lengthSquared = axisX * axisX + axisY * axisY;
        if (lengthSquared < 1e-10 || molStereoGeometry(stereo) == null)
          throw new IllegalArgumentException("Cannot represent double-bond stereo with degenerate MOL coordinates");
        for (IAtom atom : pending) {
          javax.vecmath.Point2d point = atom.getPoint2d();
          double relativeX = point.x - firstPoint.x, relativeY = point.y - firstPoint.y;
          double projection = (relativeX * axisX + relativeY * axisY) / lengthSquared;
          point.x = firstPoint.x + 2 * projection * axisX - relativeX;
          point.y = firstPoint.y + 2 * projection * axisY - relativeY;
        }
      }
      for (IDoubleBondStereochemistry stereo : elements)
        prepareMolStereoLigands(stereo, adjust);
      for (IAtom atom : molecule.atoms()) {
        javax.vecmath.Point2d point = atom.getPoint2d();
        point.x = MolCoordinateFormat.round(point.x);
        point.y = MolCoordinateFormat.round(point.y);
      }
      for (IDoubleBondStereochemistry stereo : elements) {
        if (molStereoGeometry(stereo) != stereo.getStereo())
          throw new IllegalArgumentException("Cannot preserve double-bond stereo in MOL coordinates");
        prepareMolStereoLigands(stereo, false);
        // SDG marks uncorrected cyclic bonds as unknown; verified explicit tags remain defined.
        if (stereo.getStereoBond().getStereo() == IBond.Stereo.E_OR_Z)
          stereo.getStereoBond().setStereo(IBond.Stereo.NONE);
      }
    }

    /** Use the writer's float conversion and four-decimal formatter before validation. */
    private static final class MolCoordinateFormat extends MDLV2000Writer {
      static double round(double coordinate) {
        return Double.parseDouble(formatMDLFloat((float) coordinate).trim());
      }
    }

    private static void prepareMolStereoLigands(
        IDoubleBondStereochemistry stereo, boolean adjust) {
      IBond focus = stereo.getStereoBond();
      IAtom[] references = molStereoReferences(stereo);
      javax.vecmath.Point2d first = focus.getBegin().getPoint2d(), second = focus.getEnd().getPoint2d();
      double axisX = second.x - first.x, axisY = second.y - first.y;
      double lengthSquared = axisX * axisX + axisY * axisY;
      for (int side = 0; side < 2; side++) {
        IAtom end = side == 0 ? focus.getBegin() : focus.getEnd();
        javax.vecmath.Point2d reference = references[side].getPoint2d();
        double referenceSide = axisX * (reference.y - first.y) - axisY * (reference.x - first.x);
        for (IBond carrier : end.bonds()) {
          if (carrier.equals(focus) || carrier.getOther(end).equals(references[side])) continue;
          javax.vecmath.Point2d point = carrier.getOther(end).getPoint2d();
          double ligandSide = axisX * (point.y - first.y) - axisY * (point.x - first.x);
          if (!Double.isFinite(referenceSide) || !Double.isFinite(ligandSide)
              || Math.abs(referenceSide) < 1e-8 || Math.abs(ligandSide) < 1e-8)
            throw new IllegalArgumentException("Cannot represent double-bond stereo with degenerate MOL coordinates");
          if ((referenceSide > 0) != (ligandSide > 0)) continue;
          if (!adjust)
            throw new IllegalArgumentException("MOL coordinates have inconsistent double-bond ligands");
          double relativeX = point.x - first.x, relativeY = point.y - first.y;
          double projection = (relativeX * axisX + relativeY * axisY) / lengthSquared;
          point.x = first.x + 2 * projection * axisX - relativeX;
          point.y = first.y + 2 * projection * axisY - relativeY;
        }
      }
    }

    private static IAtom[] molStereoReferences(IDoubleBondStereochemistry stereo) {
      IBond focus = stereo.getStereoBond();
      IAtom[] references = new IAtom[2];
      for (IBond carrier : stereo.getBonds()) {
        if (carrier.contains(focus.getBegin())) references[0] = carrier.getOther(focus.getBegin());
        if (carrier.contains(focus.getEnd())) references[1] = carrier.getOther(focus.getEnd());
      }
      if (focus.getOrder() != IBond.Order.DOUBLE || references[0] == null || references[1] == null)
        throw new IllegalArgumentException("Double-bond stereochemistry requires a double bond and two reference ligands");
      return references;
    }

    private static IDoubleBondStereochemistry.Conformation molStereoGeometry(
        IDoubleBondStereochemistry stereo) {
      IAtom[] references = molStereoReferences(stereo);
      javax.vecmath.Point2d first = stereo.getStereoBond().getBegin().getPoint2d();
      javax.vecmath.Point2d second = stereo.getStereoBond().getEnd().getPoint2d();
      javax.vecmath.Point2d firstReference = references[0].getPoint2d();
      javax.vecmath.Point2d secondReference = references[1].getPoint2d();
      double axisX = second.x - first.x, axisY = second.y - first.y;
      double firstSide = axisX * (firstReference.y - first.y) - axisY * (firstReference.x - first.x);
      double secondSide = axisX * (secondReference.y - second.y) - axisY * (secondReference.x - second.x);
      if (!Double.isFinite(firstSide) || !Double.isFinite(secondSide)
          || Math.abs(firstSide) < 1e-8 || Math.abs(secondSide) < 1e-8) return null;
      return (firstSide > 0) == (secondSide > 0)
          ? IDoubleBondStereochemistry.Conformation.TOGETHER
          : IDoubleBondStereochemistry.Conformation.OPPOSITE;
    }

    /** Write substructure mappings in the requested format. */
    public static void writeMappings(
        IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> mappings,
        OutType type, OutputStream os, boolean prettyJson) throws Exception {
      writeMappings(query, target, mappings, type, os, prettyJson, Double.NaN);
    }

    /** Write substructure mappings, optionally including similarity upper bound in JSON output. */
    public static void writeMappings(
        IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> mappings,
        OutType type, OutputStream os, boolean prettyJson, double similarityUpperBound)
        throws Exception {
      os=nonClosing(os);
      Objects.requireNonNull(target, "target");
      Objects.requireNonNull(mappings, "mappings");
      switch (type) {
        case JSON:
          writeMappingsJSON(query, target, mappings, os, prettyJson, similarityUpperBound);
          break;
        case SMI:
        case SMARTS:
          writeMappingsSMILES(query, target, mappings, os);
          break;
        case MOL:
          writeMappingsMOL(target, mappings, os);
          break;
        default: throw new IllegalArgumentException("Unsupported: " + type);
      }
    }

    private static OutputStream nonClosing(OutputStream output) {
      return new FilterOutputStream(Objects.requireNonNull(output,"output")) {
        @Override public void close()throws IOException {flush();}
      };
    }

    private static void writeMappingsJSON(
        IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> mappings,
        OutputStream os, boolean pretty, double similarityUpperBound) throws Exception {
      Map<String, Object> root = new LinkedHashMap<>();
      root.put("query_size", query != null ? query.getAtomCount() : null);
      root.put("target_size", target.getAtomCount());
      if (!Double.isNaN(similarityUpperBound)) {
        root.put("similarity_upper_bound", similarityUpperBound);
      }
      List<Map<String, Object>> rows = new ArrayList<>();
      int idx = 0;
      for (Map<Integer, Integer> m : mappings) {
        Map<String, Object> row = new LinkedHashMap<>();
        row.put("index", idx++);
        row.put("pairs", toPairs(m));
        row.put("query_atoms", new ArrayList<>(m.keySet()));
        row.put("target_atoms", new ArrayList<>(m.values()));
        rows.add(row);
      }
      root.put("mappings", rows);
      String json = pretty
          ? MAPPER.writerWithDefaultPrettyPrinter().writeValueAsString(root)
          : MAPPER.writeValueAsString(root);
      try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
        w.write(json);
        w.flush();
      }
    }

    private static void writeMappingsSMILES(
        IAtomContainer query, IAtomContainer target, List<Map<Integer, Integer>> mappings,
        OutputStream os) throws Exception {
      MolGraph queryGraph = query == null ? null : new MolGraph(query);
      MolGraph targetGraph = new MolGraph(target);
      Map<Integer, Integer> queryLabels = atomMapLabels(query);
      Map<Integer, Integer> targetLabels = atomMapLabels(target);
      try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
        for (Map<Integer, Integer> mapping : mappings) {
          Map<Integer, Integer> mappedQuery = new HashMap<>(queryLabels);
          Map<Integer, Integer> mappedTarget = new HashMap<>(targetLabels);
          int label = 1;
          for (Map.Entry<Integer, Integer> pair : mapping.entrySet()) {
            if (queryGraph != null) mappedQuery.put(pair.getKey(), label);
            mappedTarget.put(pair.getValue(), label++);
          }
          String querySmiles = queryGraph == null ? "" : queryGraph.toCanonicalSmiles(mappedQuery);
          String targetSmiles = targetGraph.toCanonicalSmiles(mappedTarget);
          w.write(querySmiles.isEmpty() ? targetSmiles : querySmiles + "\t" + targetSmiles);
          w.write(System.lineSeparator());
        }
        w.flush();
      }
    }

    private static void writeMappingsMOL(
        IAtomContainer target, List<Map<Integer, Integer>> mappings, OutputStream os)
        throws Exception {
      OutputStream nonClosing = os;
      for (Map<Integer, Integer> m : mappings) {
        writeMolecule(buildSubgraph(target, new HashSet<>(m.values())), OutType.MOL, nonClosing);
        OutputStreamWriter ow = new OutputStreamWriter(nonClosing, StandardCharsets.UTF_8);
        ow.write(System.lineSeparator());
        ow.flush();
      }
    }

    private static Map<Integer, Integer> atomMapLabels(IAtomContainer molecule) {
      Map<Integer, Integer> labels = new HashMap<>();
      if (molecule != null)
        for (int atom = 0; atom < molecule.getAtomCount(); atom++) {
          int label = molecule.getAtom(atom).getMapIdx();
          if (label > 0) labels.put(atom, label);
        }
      return labels;
    }

  }
}

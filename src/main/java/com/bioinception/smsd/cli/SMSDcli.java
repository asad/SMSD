/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms. */
package com.bioinception.smsd.cli;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import com.bioinception.smsd.core.Standardiser;
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
import java.util.IdentityHashMap;
import java.util.LinkedHashMap;
import java.util.List;
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
import org.openscience.cdk.io.CMLReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
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
    version = "SMSD Pro 6.12.3 — com.bioinceptionlabs",
    description =
        "SMSD — Substructure & MCS search engine by BioInception PVT LTD"
            + " (com.bioinceptionlabs).%nSupports substructure matching, MCS computation, SMARTS"
            + " queries, SDF batch processing,%nTanimoto similarity, and JSON/SMILES/MDL mapping"
            + " export.")
public class SMSDcli implements Callable<Integer> {

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

  @Option(names = {"-m", "--mappings"},
      description = "Enumerate all unique substructure mappings.")
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
  private static final SmilesGenerator SMIGEN = new SmilesGenerator(SmiFlavor.Unique);

  public static void main(String[] args) {
    System.exit(new CommandLine(new SMSDcli()).execute(args));
  }

  @Override
  public Integer call() throws Exception {
    if (benchmark) return runBenchmark();
    if ("SDF".equalsIgnoreCase(tType)) return runBatch();
    return runSingle();
  }

  private int runSingle() throws Exception {
    MolIO.Query query = MolIO.loadQuery(qType, q);
    IAtomContainer tm = applyHydrogenOptions(MolIO.loadTarget(tType, t));
    SMSD smsd = createSMSD(query, tm, new ChemOptions());
    IAtomContainer qm = query.isSmarts() ? null : applyHydrogenOptions(query.container());

    Map<String, Object> result = new LinkedHashMap<>();
    result.put("query", q);
    result.put("target", t);
    result.put("mode", mode);
    result.put("timeout_ms", timeoutMs);

    return "sub".equalsIgnoreCase(mode)
        ? runSubstructure(smsd, qm, tm, result)
        : runMCS(smsd, qm, tm, result);
  }

  private int runSubstructure(
      SMSD smsd, IAtomContainer qm, IAtomContainer tm, Map<String, Object> result)
      throws Exception {
    if (!allMappings) {
      boolean ok = smsd.isSubstructure(timeoutMs);
      result.put("exists", ok);
      writeJson(result);
      return ok ? 0 : 1;
    }
    List<Map<Integer, Integer>> maps = smsd.findAllSubstructures(10_000, timeoutMs);
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
    return 0;
  }

  private int runMCS(SMSD smsd, IAtomContainer qm, IAtomContainer tm, Map<String, Object> result)
      throws Exception {
    Map<Integer, Integer> m = smsd.findMCS(false, true, timeoutMs);
    result.put("mcs_size", m.size());
    result.put("pairs", toPairs(m));
    result.put("query_atoms", new ArrayList<>(m.keySet()));
    result.put("target_atoms", new ArrayList<>(m.values()));

    String mcsSmi = (qm != null && tm != null) ? mcsFragmentSmiles(tm, m) : null;
    if (qm != null && tm != null) {
      result.put("tanimoto", tanimoto(qm, tm, m.size()));
      result.put("similarity_upper_bound", smsd.similarityUpperBound());
      if (mcsSmi != null) result.put("mcs_smiles", mcsSmi);
    }
    writeJson(result);

    if (mcsFragmentOut != null && mcsSmi != null) {
      if ("-".equals(mcsFragmentOut)) {
        System.out.println(mcsSmi);
      } else {
        java.nio.file.Files.writeString(java.nio.file.Path.of(mcsFragmentOut), mcsSmi + "\n");
      }
    }
    return 0;
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
    final IAtomContainer queryMol = query.isSmarts() ? null : applyHydrogenOptions(query.container());

    // Resolve thread count: 0 = all available processors
    int cores = (numThreads > 0) ? numThreads : Runtime.getRuntime().availableProcessors();

    // Pre-allocate result slots so parallel writes stay at distinct indices
    Map<String, Object>[] results = new Map[n];

    ForkJoinPool pool = new ForkJoinPool(cores);
    try {
      List<java.util.concurrent.Callable<Void>> tasks = new java.util.ArrayList<>(n);
      for (int i = 0; i < n; i++) {
        final int idx = i;
        tasks.add(() -> {
          Map<String, Object> row = new LinkedHashMap<>();
          row.put("target_index", idx);
          try {
            IAtomContainer tm = applyHydrogenOptions(targets.get(idx));
            try { tm = Standardiser.standardise(tm, Standardiser.TautomerMode.NONE); }
            catch (Exception _) { }
            SMSD smsd = query.isSmarts()
                ? new SMSD(query.text(), tm, chem)
                : new SMSD(queryMol, tm, chem);
            if ("sub".equalsIgnoreCase(mode)) {
              row.put("exists", smsd.isSubstructure(timeoutMs));
            } else {
              Map<Integer, Integer> m = smsd.findMCS(false, true, timeoutMs);
              row.put("mcs_size", m.size());
              row.put("pairs", toPairs(m));
              if (queryMol != null && tm != null)
                row.put("tanimoto", tanimoto(queryMol, tm, m.size()));
            }
          } catch (Exception e) {
            row.put("error", e.getMessage());
          }
          results[idx] = row;
          return null;
        });
      }
      pool.invokeAll(tasks);
    } catch (Exception e) {
      // Sequential fallback — guarantees output even if fork/join fails
      for (int i = 0; i < n; i++) {
        if (results[i] != null) continue;
        Map<String, Object> row = new LinkedHashMap<>();
        row.put("target_index", i);
        try {
          IAtomContainer tm = applyHydrogenOptions(targets.get(i));
          try { tm = Standardiser.standardise(tm, Standardiser.TautomerMode.NONE); }
          catch (Exception _) { }
          SMSD smsd = query.isSmarts()
              ? new SMSD(query.text(), tm, chem)
              : new SMSD(queryMol, tm, chem);
          if ("sub".equalsIgnoreCase(mode)) {
            row.put("exists", smsd.isSubstructure(timeoutMs));
          } else {
            Map<Integer, Integer> m = smsd.findMCS(false, true, timeoutMs);
            row.put("mcs_size", m.size());
            row.put("pairs", toPairs(m));
            if (queryMol != null && tm != null)
              row.put("tanimoto", tanimoto(queryMol, tm, m.size()));
          }
        } catch (Exception ex) {
          row.put("error", ex.getMessage());
        }
        results[i] = row;
      }
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
    return 0;
  }

  private int runBenchmark() throws Exception {
    MolIO.Query query = MolIO.loadQuery(qType, q);
    IAtomContainer tm = applyHydrogenOptions(MolIO.loadTarget(tType, t));
    ChemOptions chem = new ChemOptions();
    int warmupRuns = 3, timedRuns = 10;
    PrintWriter pw = new PrintWriter(System.out, true);
    pw.println("SMSD Benchmark (com.bioinceptionlabs)");
    pw.println("======================================");
    pw.printf("Query:  %s%nTarget: %s%nWarmup: %d runs, Timed: %d runs%n%n",
        q, t, warmupRuns, timedRuns);

    for (int i = 0; i < warmupRuns; i++) {
      SMSD smsd = createSMSD(query, tm, chem);
      smsd.isSubstructure(timeoutMs);
      if (!query.isSmarts()) smsd.findMCS(false, true, timeoutMs);
    }

    long subTotal = 0;
    boolean subResult = false;
    for (int i = 0; i < timedRuns; i++) {
      SMSD smsd = createSMSD(query, tm, chem);
      long t0 = System.nanoTime();
      subResult = smsd.isSubstructure(timeoutMs);
      subTotal += System.nanoTime() - t0;
    }
    pw.printf("Substructure: %s  avg %.3f ms%n",
        subResult ? "FOUND" : "NOT FOUND", (subTotal / (double) timedRuns) / 1_000_000.0);

    if (!query.isSmarts()) {
      long mcsTotal = 0;
      int mcsSize = 0;
      for (int i = 0; i < timedRuns; i++) {
        SMSD smsd = createSMSD(query, tm, chem);
        long t0 = System.nanoTime();
        Map<Integer, Integer> m = smsd.findMCS(false, true, timeoutMs);
        mcsTotal += System.nanoTime() - t0;
        mcsSize = m.size();
      }
      pw.printf("MCS:          size %d  avg %.3f ms%n",
          mcsSize, (mcsTotal / (double) timedRuns) / 1_000_000.0);
      pw.printf("Sim UB:       %.4f%n", createSMSD(query, tm, chem).similarityUpperBound());
    }
    pw.println("======================================");
    return 0;
  }

  // ---- helpers ----

  private SMSD createSMSD(MolIO.Query query, IAtomContainer tm, ChemOptions chem) throws Exception {
    return query.isSmarts()
        ? new SMSD(query.text(), tm, chem)
        : new SMSD(applyHydrogenOptions(query.container()), tm, chem);
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

  static IAtomContainer buildSubgraph(IAtomContainer target, Set<Integer> atomIdxs) {
    IAtomContainer sub = target.getBuilder().newInstance(IAtomContainer.class);
    Map<Integer, IAtom> idx2new = new HashMap<>();
    for (int i : atomIdxs) {
      IAtom a = target.getAtom(i);
      IAtom b = sub.newAtom(a.getAtomicNumber() != null ? a.getAtomicNumber() : 0);
      b.setImplicitHydrogenCount(a.getImplicitHydrogenCount());
      b.setIsAromatic(a.isAromatic());
      b.setPoint2d(a.getPoint2d());
      b.setPoint3d(a.getPoint3d());
      b.setCharge(a.getCharge());
      b.setFormalCharge(a.getFormalCharge());
      idx2new.put(i, b);
    }
    for (IBond bo : target.bonds()) {
      int ai = target.indexOf(bo.getAtom(0));
      int bi = target.indexOf(bo.getAtom(1));
      if (idx2new.containsKey(ai) && idx2new.containsKey(bi))
        sub.newBond(idx2new.get(ai), idx2new.get(bi), bo.getOrder());
    }
    return sub;
  }

  private static String mcsFragmentSmiles(IAtomContainer target, Map<Integer, Integer> mapping) {
    if (mapping == null || mapping.isEmpty() || target == null) return null;
    try {
      return SMIGEN.create(buildSubgraph(target, new HashSet<>(mapping.values())));
    } catch (Exception e) {
      return null;
    }
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
      PrintWriter pw = new PrintWriter(System.out);
      pw.println(payload);
      pw.flush();
    } else {
      java.nio.file.Files.writeString(java.nio.file.Path.of(jsonOut), payload);
    }
  }

  private static OutputUtil.OutType parseOutType(String s) {
    return switch ((s == null ? "json" : s).trim().toLowerCase()) {
      case "json"                -> OutputUtil.OutType.JSON;
      case "smiles", "smi"      -> OutputUtil.OutType.SMI;
      case "smarts"             -> OutputUtil.OutType.SMARTS;
      case "mdl", "mol", "sdf"  -> OutputUtil.OutType.MOL;
      default -> throw new IllegalArgumentException("Unknown --map-format: " + s);
    };
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
    public record Query(boolean isSmarts, String text, IAtomContainer container) {}

    /**
     * Load a query molecule or SMARTS pattern.
     *
     * @param type input type: SMI (SMILES), SIG (SMARTS), MOL, ML2, PDB, or CML
     * @param value SMILES/SMARTS string or file path
     */
    public static Query loadQuery(String type, String value) throws Exception {
      Objects.requireNonNull(type, "Query type must not be null");
      Objects.requireNonNull(value, "Query value must not be null");
      return switch (type.toUpperCase()) {
        case "SIG" -> new Query(true, value, null);
        case "SMI" -> new Query(false, value, SP.get().parseSmiles(value));
        default    -> new Query(false, value, readFile(type, value));
      };
    }

    /**
     * Load a single target molecule.
     *
     * @param type input type: SMI (SMILES), MOL, ML2, PDB, or CML
     * @param value SMILES string or file path
     */
    public static IAtomContainer loadTarget(String type, String value) throws Exception {
      Objects.requireNonNull(type, "Target type must not be null");
      Objects.requireNonNull(value, "Target value must not be null");
      return switch (type.toUpperCase()) {
        case "SMI" -> SP.get().parseSmiles(value);
        case "SDF" -> throw new IllegalArgumentException("Use loadTargetsSDF for multi-molecule SDF.");
        default    -> readFile(type, value);
      };
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
          out.add(it.next());
        }
      }
      return out;
    }

    private static IAtomContainer readFile(String type, String path) throws Exception {
      Function<InputStream, ISimpleChemObjectReader> factory = switch (type.toUpperCase()) {
        case "MOL" -> MDLV2000Reader::new;
        case "ML2" -> Mol2Reader::new;
        case "PDB" -> PDBReader::new;
        case "CML" -> CMLReader::new;
        default -> throw new IllegalArgumentException("Unsupported file type: " + type);
      };
      try (FileInputStream fis = new FileInputStream(path);
          ISimpleChemObjectReader reader = factory.apply(fis)) {
        return reader.read(
            DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
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
      if (type == OutType.MOL) {
        try (MDLV2000Writer w = new MDLV2000Writer(os)) { w.write(mol); }
      } else {
        try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
          w.write(new SmilesGenerator(SmiFlavor.Unique).create(mol));
          w.write(System.lineSeparator());
          w.flush();
        }
      }
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
      Objects.requireNonNull(target, "target");
      Objects.requireNonNull(mappings, "mappings");
      switch (type) {
        case JSON -> writeMappingsJSON(query, target, mappings, os, prettyJson, similarityUpperBound);
        case SMI, SMARTS -> writeMappingsSMILES(query, target, mappings, os);
        case MOL -> writeMappingsMOL(target, mappings, os);
        default -> throw new IllegalArgumentException("Unsupported: " + type);
      }
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
      SmilesGenerator sg = new SmilesGenerator(
          SmiFlavor.Unique | SmiFlavor.UseAromaticSymbols | SmiFlavor.AtomAtomMap);
      Map<IAtom, Object> qBackup = snapshotAtomMaps(query);
      Map<IAtom, Object> tBackup = snapshotAtomMaps(target);
      try (OutputStreamWriter w = new OutputStreamWriter(os, StandardCharsets.UTF_8)) {
        for (Map<Integer, Integer> m : mappings) {
          assignAtomMaps(query, target, m);
          try {
            String qs = query != null ? sg.create(query) : "";
            String ts = sg.create(target);
            w.write(qs.isEmpty() ? ts : qs + "\t" + ts);
            w.write(System.lineSeparator());
          } finally {
            restoreAtomMaps(query, qBackup);
            restoreAtomMaps(target, tBackup);
          }
        }
        w.flush();
      }
    }

    private static void writeMappingsMOL(
        IAtomContainer target, List<Map<Integer, Integer>> mappings, OutputStream os)
        throws Exception {
      OutputStream nonClosing = new FilterOutputStream(os) {
        @Override public void close() throws IOException { flush(); }
      };
      for (Map<Integer, Integer> m : mappings) {
        try (MDLV2000Writer w = new MDLV2000Writer(nonClosing)) {
          w.write(buildSubgraph(target, new HashSet<>(m.values())));
        }
        OutputStreamWriter ow = new OutputStreamWriter(nonClosing, StandardCharsets.UTF_8);
        ow.write(System.lineSeparator());
        ow.flush();
      }
    }

    private static void assignAtomMaps(
        IAtomContainer query, IAtomContainer target, Map<Integer, Integer> mapping) {
      int n = 1;
      for (Map.Entry<Integer, Integer> e : mapping.entrySet()) {
        if (query != null) query.getAtom(e.getKey()).setProperty(CDKConstants.ATOM_ATOM_MAPPING, n);
        target.getAtom(e.getValue()).setProperty(CDKConstants.ATOM_ATOM_MAPPING, n);
        n++;
      }
    }

    private static Map<IAtom, Object> snapshotAtomMaps(IAtomContainer mol) {
      if (mol == null) return Collections.emptyMap();
      Map<IAtom, Object> snap = new IdentityHashMap<>();
      for (IAtom a : mol.atoms()) snap.put(a, a.getProperty(CDKConstants.ATOM_ATOM_MAPPING));
      return snap;
    }

    private static void restoreAtomMaps(IAtomContainer mol, Map<IAtom, Object> snapshot) {
      if (mol == null) return;
      snapshot.forEach((a, v) -> a.setProperty(CDKConstants.ATOM_ATOM_MAPPING, v));
    }
  }
}

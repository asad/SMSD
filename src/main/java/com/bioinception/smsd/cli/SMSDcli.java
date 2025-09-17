/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd.cli;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import picocli.CommandLine;
import picocli.CommandLine.Command;
import picocli.CommandLine.Option;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * SMSDcli
 * -------
 * Picocli-based CLI entrypoint. Supports substructure and MCS modes with JSON output.
 * Adds a global --timeout option (ms) that safely bounds runtime.
 */
@Command(name = "smsd",
        mixinStandardHelpOptions = true,
        version = "3.0.0",
        description = "SMSD (CDK-compat) — substructure & MCS with advanced SMARTS and JSON.")
public class SMSDcli implements Callable<Integer> {

    // ---- shared flags
    @Option(names={"--Q"}, required = true, description = "Query type (SMI, MOL, SDF). For now: SMI only.")
    String qType;

    @Option(names={"--q"}, required = true, description = "Query string or path (by type).")
    String q;

    @Option(names={"--T"}, required = true, description = "Target type (SMI, MOL, SDF). For now: SMI only.")
    String tType;

    @Option(names={"--t"}, required = true, description = "Target string or path (by type).")
    String t;

    @Option(names={"-m", "--mappings"}, description = "Return all unique substructure mappings (JSON).")
    boolean allMappings;

    @Option(names={"--json"}, description = "Write JSON to file path or '-' for stdout.")
    String jsonOut;

    @Option(names={"--json-pretty"}, description = "Pretty-print JSON.")
    boolean jsonPretty;

    @Option(names={"--timeout"}, description = "Time limit in milliseconds (applies to the chosen mode). Default: 10000")
    long timeoutMs = 10_000L;

    @Option(names={"--mode"}, description = "Mode: sub (substructure), mcs (MCS). Default: sub", defaultValue = "sub")
    String mode;

    public static void main(String[] args) {
        int code = new CommandLine(new SMSDcli()).execute(args);
        System.exit(code);
    }

    @Override
    public Integer call() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer qm = sp.parseSmiles(q);
        IAtomContainer tm = sp.parseSmiles(t);

        ChemOptions chem = new ChemOptions();
        SMSD smsd = new SMSD(qm, tm, chem);

        Map<String,Object> result = new LinkedHashMap<>();
        result.put("query", q);
        result.put("target", t);
        result.put("mode", mode);
        result.put("timeout_ms", timeoutMs);

        if ("sub".equalsIgnoreCase(mode)) {
            if (allMappings) {
                List<Map<Integer,Integer>> maps = smsd.findAllSubstructures(10_000, timeoutMs);
                List<Map<String,Object>> items = new ArrayList<>();
                int idx = 0;
                for (Map<Integer,Integer> m : maps) {
                    Map<String,Object> row = new LinkedHashMap<>();
                    row.put("index", idx++);
                    row.put("pairs", toPairs(m));
                    row.put("query_atoms", new ArrayList<>(m.keySet()));
                    row.put("target_atoms", new ArrayList<>(m.values()));
                    row.put("target_bonds", Collections.emptyList()); // keep shape
                    items.add(row);
                }
                result.put("mappings", items);
                writeJson(result);
                return 0;
            } else {
                boolean ok = smsd.isSubstructure(timeoutMs);
                result.put("exists", ok);
                writeJson(result);
                return ok ? 0 : 1;
            }
        } else {
            // mcs
            Map<Integer,Integer> m = smsd.findMCS(false, true, timeoutMs);
            result.put("mcs_type", "MCCS");
            result.put("pairs", toPairs(m));
            result.put("query_atoms", new ArrayList<>(m.keySet()));
            result.put("target_atoms", new ArrayList<>(m.values()));
            writeJson(result);
            return 0;
        }
    }

    private static List<List<Integer>> toPairs(Map<Integer,Integer> m) {
        List<List<Integer>> out = new ArrayList<>();
        for (Map.Entry<Integer,Integer> e : m.entrySet()) {
            out.add(Arrays.asList(e.getKey(), e.getValue()));
        }
        return out;
    }

    private void writeJson(Map<String,Object> data) throws Exception {
        if (jsonOut == null) return;
        String payload = (jsonPretty ? prettyJson(data) : compactJson(data));
        if ("-".equals(jsonOut)) {
            PrintWriter pw = new PrintWriter(System.out);
            pw.println(payload);
            pw.flush();
        } else {
            java.nio.file.Files.write(java.nio.file.Path.of(jsonOut), payload.getBytes(java.nio.charset.StandardCharsets.UTF_8));
        }
    }

    private static String compactJson(Object o) throws Exception {
        return new com.fasterxml.jackson.databind.ObjectMapper().writeValueAsString(o);
    }
    private static String prettyJson(Object o) throws Exception {
        return new com.fasterxml.jackson.databind.ObjectMapper().writerWithDefaultPrettyPrinter().writeValueAsString(o);
    }
}

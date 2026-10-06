/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import com.bioinception.smsd.core.MolGraph;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.*;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.condition.EnabledIfSystemProperty;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * External benchmark tests using repository-maintained corpora:
 * <ul>
 *   <li>Tautobase — 468 tautomer pairs (Chodera/Wahl-Sander)</li>
 *   <li>Derived MoleculeNet random pairs — 1000 checked-in rows</li>
 *   <li>Legacy derived neighbor pairs — 1000 rows, variable similarity and repeated source IDs</li>
 *   <li>Stress pairs — 12 adversarial graph-theory hard cases</li>
 *   <li>Ehrlich-Rarey SMARTS — 1400 substructure patterns</li>
 * </ul>
 *
 * Run with: mvn test -Dtest=ExternalBenchmarkTest -Dbenchmark=true
 *
 * @author Syed Asad Rahman
 */
@DisplayName("External Benchmark Tests")
@Timeout(value = 30, unit = TimeUnit.MINUTES, threadMode = Timeout.ThreadMode.SEPARATE_THREAD)
public class ExternalBenchmarkTest extends TestBase {

    private static final String DATA_DIR = "benchmarks/data";
    private static final long MCS_TIMEOUT_MS = BenchmarkRunSettings.TIMEOUT_MS;
    private static final long SUB_TIMEOUT_MS = BenchmarkRunSettings.TIMEOUT_MS;


    /** Each search has a bounded budget; witness validation and checkpoints are outside its timer. */
    private static Map<Integer, Integer> benchmarkMCS(IAtomContainer query, IAtomContainer target,
            ChemOptions chemistry, String suite, String row) {
        Map<Integer, Integer> result = Collections.emptyMap();
        for (int trial = 0; trial < BenchmarkRunSettings.WARMUP + BenchmarkRunSettings.ROUNDS; trial++) {
            BenchmarkRunSettings.checkInterrupted();
            long start = System.nanoTime();
            try {
                result = new SMSD(query, target, chemistry).findMCS(false, false, MCS_TIMEOUT_MS);
            } catch (Exception exception) {
                BenchmarkRunSettings.checkpoint(suite, row, "ERROR", 0, System.nanoTime() - start,
                    exception.toString());
                throw new AssertionError("Benchmark search failed for " + row, exception);
            }
            long elapsed = System.nanoTime() - start;
            BenchmarkRunSettings.validate(query, target, result, chemistry, row);
            if (trial >= BenchmarkRunSettings.WARMUP)
                BenchmarkRunSettings.checkpoint(suite, row + "/" + trial, "OK", result.size(), elapsed,
                    "induced=false connected=false");
        }
        return result;
    }

    // ======================================================================
    // Helper: load TSV pairs (skip comment lines starting with #)
    // ======================================================================
    private static List<String[]> loadTsvPairs(String filename) throws IOException {
        Path p = Paths.get(DATA_DIR, filename);
        if (!Files.exists(p)) {
            System.err.println("SKIP: " + p + " not found");
            return Collections.emptyList();
        }
        return Files.readAllLines(p, StandardCharsets.UTF_8).stream()
                .filter(line -> !line.startsWith("#") && !line.isBlank())
                .map(line -> line.split("\t"))
                .filter(parts -> parts.length >= 2)
                .collect(Collectors.toList());
    }

    // ======================================================================
    // 1. Tautobase — Tautomer-Aware MCS (468 pairs)
    // ======================================================================
    @Nested
    @DisplayName("Tautobase Tautomer Benchmark (468 pairs)")
    class TautobaseBenchmark {

        @Test
        @DisplayName("Tautomer-aware MCS: measure full recovery and partial gains")
        @EnabledIfSystemProperty(named = "benchmark", matches = "true")
        void tautobaseFullRecovery() throws Exception {
            Path p = Paths.get(DATA_DIR, "chodera_tautobase_subset.txt");
            if (!Files.exists(p)) {
                System.err.println("SKIP: " + p + " not found");
                return;
            }

            List<String> lines = Files.readAllLines(p, StandardCharsets.UTF_8);
            int total = 0, fullMatch = 0, partialMatch = 0, failed = 0;
            int totalGain = 0;

            ChemOptions tautOpts = new ChemOptions();
            tautOpts.tautomerAware = true;
            ChemOptions strictOpts = new ChemOptions();

            for (String line : lines) {
                BenchmarkRunSettings.checkInterrupted();
                if (line.startsWith("name") || line.isBlank()) continue;
                String[] parts = line.split(",\\s*");
                if (parts.length < 3) continue;

                String smi1 = parts[1].trim();
                String smi2 = parts[2].trim();

                try {
                    IAtomContainer m1 = mol(smi1);
                    IAtomContainer m2 = mol(smi2);
                    if (m1 == null || m2 == null) continue;

                    int maxAtoms = Math.min(m1.getAtomCount(), m2.getAtomCount());

                    // Tautomer-aware MCS
                    SMSD smsdTaut = new SMSD(m1, m2, tautOpts);
                    Map<Integer,Integer> tautMap = benchmarkMCS(m1, m2, tautOpts, "tautobase-tautomer", parts[0]);
                    int tautMCS = tautMap.size();

                    // Strict MCS (for comparison)
                    SMSD smsdStrict = new SMSD(m1, m2, strictOpts);
                    Map<Integer,Integer> strictMap = benchmarkMCS(m1, m2, strictOpts, "tautobase-strict", parts[0]);
                    int strictMCS = strictMap.size();

                    total++;
                    int gain = tautMCS - strictMCS;
                    if (gain > 0) totalGain += gain;

                    if (tautMCS >= maxAtoms) {
                        fullMatch++;
                    } else if (tautMCS > strictMCS) {
                        partialMatch++;
                    } else {
                        failed++;
                    }
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("tautobase-tautomer", parts[0], "PARSE_ERROR", 0, 0, e.toString());
                }
            }

            System.out.printf("Tautobase results: %d pairs, %d full match (%.1f%%), " +
                            "%d partial gain, %d no gain, total +%d atoms recovered%n",
                    total, fullMatch, 100.0 * fullMatch / Math.max(total, 1),
                    partialMatch, failed, totalGain);

            // Coverage is checked; partial mappings need not show tautomer gain.
            assertTrue(total > 100, "Should parse at least 100 Tautobase pairs");
        }

        @Test
        @DisplayName("Tautobase budget crossings and valid mappings")
        @EnabledIfSystemProperty(named = "benchmark", matches = "true")
        void tautobaseNoTimeout() throws Exception {
            Path p = Paths.get(DATA_DIR, "chodera_tautobase_subset.txt");
            if (!Files.exists(p)) return;

            List<String> lines = Files.readAllLines(p, StandardCharsets.UTF_8);
            int timeouts = 0, total = 0;
            long totalTimeUs = 0;

            for (String line : lines) {
                BenchmarkRunSettings.checkInterrupted();
                if (line.startsWith("name") || line.isBlank()) continue;
                String[] parts = line.split(",\\s*");
                if (parts.length < 3) continue;

                try {
                    IAtomContainer m1 = mol(parts[1].trim());
                    IAtomContainer m2 = mol(parts[2].trim());
                    if (m1 == null || m2 == null) continue;

                    long t0 = System.nanoTime();
                    SMSD smsd = new SMSD(m1, m2, new ChemOptions());
                    benchmarkMCS(m1, m2, new ChemOptions(), "tautobase-timing", parts[0]);
                    long elapsed = (System.nanoTime() - t0) / 1000;
                    totalTimeUs += elapsed;
                    total++;

                    if (elapsed > MCS_TIMEOUT_MS * 1000) timeouts++;
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("external-parse", Arrays.toString(parts), "PARSE_ERROR", 0, 0, e.toString());
                }
            }

            System.out.printf("Tautobase timing: %d pairs, avg %.1f us, %d elapsed budget crossings%n",
                    total, (double) totalTimeUs / Math.max(total, 1), timeouts);
            assertTrue(total > 100, "Should parse at least 100 Tautobase pairs");
        }
    }

    // ======================================================================
    // 2. Dalke-style MCS Benchmark (Random + Nearest-Neighbor)
    // ======================================================================
    @Nested
    @DisplayName("Dalke-style MCS Benchmark (2000 pairs)")
    class DalkeBenchmark {

        @Test
        @DisplayName("Random pairs: measure validated MCS and LFUB attainment")
        @EnabledIfSystemProperty(named = "benchmark", matches = "true")
        void randomPairsLfubRate() throws Exception {
            List<String[]> pairs = loadTsvPairs("dalke_random_pairs.tsv");
            if (pairs.isEmpty()) return;

            int total = 0, lfubHit = 0, timeouts = 0;
            long totalTimeUs = 0;

            for (String[] parts : pairs) {
                BenchmarkRunSettings.checkInterrupted();
                try {
                    IAtomContainer m1 = mol(parts[0]);
                    IAtomContainer m2 = mol(parts[1]);
                    if (m1 == null || m2 == null) continue;

                    ChemOptions opts = new ChemOptions();
                    long t0 = System.nanoTime();
                    SMSD smsd = new SMSD(m1, m2, opts);
                    Map<Integer, Integer> mcsMap = benchmarkMCS(m1, m2, opts, "random", Integer.toString(total));
                    long elapsed = (System.nanoTime() - t0) / 1000;
                    totalTimeUs += elapsed;
                    total++;

                    int mcsSize = mcsMap.size();
                    // Compute LFUB via static method
                    MolGraph g1 = new MolGraph(m1);
                    MolGraph g2 = new MolGraph(m2);
                    int lfub = SMSD.labelFrequencyUpperBound(g1, g2, opts);
                    if (mcsSize == lfub) lfubHit++;
                    if (elapsed > MCS_TIMEOUT_MS * 1000) timeouts++;
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("external-parse", Arrays.toString(parts), "PARSE_ERROR", 0, 0, e.toString());
                }
            }

            double lfubRate = 100.0 * lfubHit / Math.max(total, 1);
            System.out.printf("Dalke random: %d pairs, LFUB hit %.1f%%, avg %.1f us, %d elapsed budget crossings%n",
                    total, lfubRate, (double) totalTimeUs / Math.max(total, 1), timeouts);

            assertTrue(total > 500, "Should parse at least 500 random pairs");
        }

        @Test
        @DisplayName("Neighbor pairs: measure validated MCS sizes")
        @EnabledIfSystemProperty(named = "benchmark", matches = "true")
        void nnPairsMCSQuality() throws Exception {
            List<String[]> pairs = loadTsvPairs("dalke_nn_pairs.tsv");
            if (pairs.isEmpty()) return;

            int total = 0, mcsGe5 = 0, timeouts = 0;
            long totalTimeUs = 0;

            for (String[] parts : pairs) {
                BenchmarkRunSettings.checkInterrupted();
                try {
                    IAtomContainer m1 = mol(parts[0]);
                    IAtomContainer m2 = mol(parts[1]);
                    if (m1 == null || m2 == null) continue;

                    long t0 = System.nanoTime();
                    SMSD smsd = new SMSD(m1, m2, new ChemOptions());
                    Map<Integer, Integer> nnMap = benchmarkMCS(m1, m2, new ChemOptions(), "neighbor", Integer.toString(total));
                    long elapsed = (System.nanoTime() - t0) / 1000;
                    totalTimeUs += elapsed;
                    total++;

                    if (nnMap.size() >= 5) mcsGe5++;
                    if (elapsed > MCS_TIMEOUT_MS * 1000) timeouts++;
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("external-parse", Arrays.toString(parts), "PARSE_ERROR", 0, 0, e.toString());
                }
            }

            double qualRate = 100.0 * mcsGe5 / Math.max(total, 1);
            System.out.printf("Dalke NN: %d pairs, MCS>=5 %.1f%%, avg %.1f us, %d elapsed budget crossings%n",
                    total, qualRate, (double) totalTimeUs / Math.max(total, 1), timeouts);

            assertTrue(total > 500, "Should parse at least 500 NN pairs");
        }
    }

    // ======================================================================
    // 3. Stress Pairs — Adversarial Hard Cases (12 pairs)
    // ======================================================================
    @Nested
    @DisplayName("Stress Test Hard Cases (12 pairs)")
    class StressTestBenchmark {

        @Test
        @DisplayName("All stress pairs execute with bounded budgets")
        void stressPairsNoTimeout() throws Exception {
            List<String[]> pairs = loadTsvPairs("stress_pairs.tsv");
            if (pairs.isEmpty()) return;

            int total = 0, timeouts = 0;

            for (String[] parts : pairs) {
                BenchmarkRunSettings.checkInterrupted();
                String name = parts.length > 2 ? parts[2] : "pair_" + total;
                try {
                    IAtomContainer m1 = mol(parts[0]);
                    IAtomContainer m2 = mol(parts[1]);
                    if (m1 == null || m2 == null) continue;

                    long t0 = System.nanoTime();
                    SMSD smsd = new SMSD(m1, m2, new ChemOptions());
                    Map<Integer, Integer> stressMap = benchmarkMCS(m1, m2, new ChemOptions(), "stress", name);
                    long elapsedMs = (System.nanoTime() - t0) / 1_000_000;
                    int mcs = stressMap.size();
                    total++;

                    System.out.printf("  %-35s MCS=%2d  time=%,dms%n", name, mcs, elapsedMs);

                    if (elapsedMs > MCS_TIMEOUT_MS) {
                        timeouts++;
                        System.out.printf("  ** ELAPSED BUDGET CROSSING on %s%n", name);
                    }
                } catch (Exception e) {
                    System.out.printf("  %-35s PARSE ERROR: %s%n", name, e.getMessage());
                    BenchmarkRunSettings.checkpoint("stress", name, "PARSE_ERROR", 0, 0, e.toString());
                }
            }

            System.out.printf("Stress test: %d pairs, %d elapsed budget crossings%n", total, timeouts);
            assertTrue(total > 0, "Should parse at least 1 stress pair");

        }

        @Test
        @DisplayName("Self-match pairs return exact atom count")
        void selfMatchExact() throws Exception {
            List<String[]> pairs = loadTsvPairs("stress_pairs.tsv");
            if (pairs.isEmpty()) return;

            for (String[] parts : pairs) {
                BenchmarkRunSettings.checkInterrupted();
                if (parts.length < 3) continue;
                String name = parts[2];
                if (!name.contains("self")) continue;

                try {
                    IAtomContainer m1 = mol(parts[0]);
                    IAtomContainer m2 = mol(parts[1]);
                    if (m1 == null || m2 == null) continue;

                    SMSD smsd = new SMSD(m1, m2, new ChemOptions());
                    Map<Integer, Integer> selfMap = benchmarkMCS(m1, m2, new ChemOptions(), "stress-self", name);
                    int mcs = selfMap.size();
                    int expected = m1.getAtomCount();

                    assertEquals(expected, mcs,
                            name + ": self-match should return all " + expected + " atoms");
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("stress-self", name, "PARSE_ERROR", 0, 0, e.toString());
                }
            }
        }
    }

    // ======================================================================
    // 4. Ehrlich-Rarey SMARTS — Substructure Search (1400 patterns)
    // ======================================================================
    @Nested
    @DisplayName("Ehrlich-Rarey Substructure Benchmark (1400 SMARTS)")
    class EhrlichRareyBenchmark {

        @Test
        @DisplayName("All SMARTS patterns: coverage and observed latency")
        @EnabledIfSystemProperty(named = "benchmark", matches = "true")
        void smartsQueryPerformance() throws Exception {
            Path p = Paths.get(DATA_DIR, "ehrlich_rarey_smarts.txt");
            if (!Files.exists(p)) {
                System.err.println("SKIP: " + p + " not found");
                return;
            }

            // Load SMARTS patterns
            List<String> smarts = Files.readAllLines(p, StandardCharsets.UTF_8).stream()
                    .filter(line -> !line.startsWith("#") && !line.isBlank())
                    .map(line -> line.split("\t")[0].trim())
                    .filter(s -> !s.isEmpty())
                    .collect(Collectors.toList());

            // Target molecule: a moderately complex drug (ibuprofen)
            IAtomContainer target = mol("CC(C)Cc1ccc(CC(C)C(O)=O)cc1");

            int total = 0, matched = 0;
            long totalTimeNs = 0;
            List<Long> times = new ArrayList<>();

            for (String sma : smarts) {
                BenchmarkRunSettings.checkInterrupted();
                try {
                    ChemOptions smOpts = new ChemOptions();
                    long t0 = System.nanoTime();
                    boolean hit = false;
                    for (int trial = 0; trial < BenchmarkRunSettings.WARMUP + BenchmarkRunSettings.ROUNDS; trial++) {
                        BenchmarkRunSettings.checkInterrupted();
                        long trialStart = System.nanoTime();
                        hit = new SMSD(sma, target, smOpts).isSubstructure(SUB_TIMEOUT_MS);
                        if (trial >= BenchmarkRunSettings.WARMUP)
                            BenchmarkRunSettings.checkpoint("smarts", sma + "/" + trial, "OK", 0,
                                System.nanoTime() - trialStart, "matched=" + hit);
                    }
                    long elapsed = System.nanoTime() - t0;
                    totalTimeNs += elapsed;
                    times.add(elapsed / 1000); // microseconds
                    total++;
                    if (hit) matched++;
                } catch (Exception e) {
                    BenchmarkRunSettings.checkpoint("smarts", sma, "UNSUPPORTED", 0, 0, e.toString());
                }
            }

            Collections.sort(times);
            long medianUs = times.isEmpty() ? 0 : times.get(times.size() / 2);

            System.out.printf("Ehrlich-Rarey: %d patterns, %d matched, " +
                            "median %d us, total %.1f ms%n",
                    total, matched, medianUs, totalTimeNs / 1e6);

            assertTrue(total > 500, "Should parse at least 500 SMARTS patterns");
        }
    }
}

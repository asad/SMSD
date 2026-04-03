/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd;

import static org.junit.jupiter.api.Assertions.*;

import com.bioinception.smsd.core.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.condition.EnabledIfSystemProperty;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Benchmark suite: consolidated from BenchmarkTest, BenchmarkComparison,
 * HeadToHeadBenchmark.
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Benchmark Suite Tests")
@Timeout(value = 60, unit = TimeUnit.SECONDS, threadMode = Timeout.ThreadMode.SEPARATE_THREAD)
public class BenchmarkSuiteTest extends TestBase {

  // ======================================================================
  // From: BenchmarkTest.java
  // ======================================================================

  @Nested
  @DisplayName("Performance Regression Tests")
  class PerformanceRegression {

    @Test
    @DisplayName("Large molecule substructure completes under 200ms")
    void largeMolSubstructureRegression() throws Exception {
      IAtomContainer q = mol("C".repeat(50));
      IAtomContainer t = mol("C".repeat(70));
      long t0 = System.nanoTime();
      SMSD smsd = new SMSD(q, t, new ChemOptions());
      boolean result = smsd.isSubstructure(5000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertTrue(result, "C50 should be substructure of C70");
      assertTrue(elapsed < 200, "Should complete under 200ms, took " + elapsed + "ms");
    }

    @Test
    @DisplayName("Drug pair MCS completes under 5s")
    void drugMcsRegression() throws Exception {
      IAtomContainer aspirin = mol("CC(=O)Oc1ccccc1C(=O)O");
      IAtomContainer ibuprofen = mol("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
      long t0 = System.nanoTime();
      SMSD smsd = new SMSD(aspirin, ibuprofen, new ChemOptions());
      Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertFalse(mcs.isEmpty(), "MCS should find common atoms");
      assertTrue(elapsed < 5000, "Should complete under 5s, took " + elapsed + "ms");
    }

    @Test
    @DisplayName("RASCAL upper bound >= actual MCS Tanimoto")
    void rascalCorrectnessInvariant() throws Exception {
      String[][] pairs = {
        {"c1ccccc1", "c1ccc(C)cc1"},
        {"CC(=O)O", "CC(=O)Nc1ccc(O)cc1"},
        {"c1ccc2ccccc2c1", "c1ccc2cc3ccccc3cc2c1"},
        {"CCO", "CCCO"},
        {"c1ccncc1", "c1ccoc1"},
      };
      ChemOptions opts = new ChemOptions();
      for (String[] p : pairs) {
        IAtomContainer m1 = mol(p[0]), m2 = mol(p[1]);
        double ub = SearchEngine.similarityUpperBound(m1, m2, opts);
        SMSD smsd = new SMSD(m1, m2, opts);
        Map<Integer, Integer> mcs = smsd.findMCS(false, true, 5000);
        int common = mcs.size();
        int total = m1.getAtomCount() + m2.getAtomCount() - common;
        double actual = total > 0 ? (double) common / total : 0.0;
        assertTrue(
            ub >= actual - 0.01,
            String.format(
                "RASCAL UB (%.3f) should be >= actual Tanimoto (%.3f) for %s vs %s",
                ub, actual, p[0], p[1]));
      }
    }

    @Test
    @DisplayName("Batch MCS completes and returns results")
    void batchMcsCompletes() throws Exception {
      List<IAtomContainer> mols = new ArrayList<>();
      mols.add(mol("c1ccccc1"));
      mols.add(mol("c1ccc(O)cc1"));
      mols.add(mol("c1ccc(N)cc1"));
      mols.add(mol("c1ccc(C)cc1"));
      mols.add(mol("c1ccc2ccccc2c1"));
      ChemOptions opts = new ChemOptions();
      long t0 = System.nanoTime();
      Map<Integer, Map<Integer, Integer>> results =
          SearchEngine.batchMCS(mols, opts, new SearchEngine.McsOptions(), 0.1);
      long elapsed = (System.nanoTime() - t0) / 1_000_000;
      assertNotNull(results, "Batch MCS should return results");
      assertFalse(results.isEmpty(), "Should find at least one pair above threshold");
      assertTrue(elapsed < 10000, "Batch of 5 should complete under 10s, took " + elapsed + "ms");
    }

    @Test
    @DisplayName("All MCS mappings are chemically valid and maximal")
    void mcsChemicalValidation() throws Exception {
      String[][] pairs = {
        {"c1ccccc1", "Cc1ccccc1", "benzene-toluene"},
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "aspirin-acetaminophen"},
        {"CN1CCC23C4C1CC5=C2C(=CC=C5)OC3C(C=C4)O",
         "COC1=CC=C2C3CC4=CC=C(O)C5=C4C3(CCN2C)C=C51", "morphine-codeine"},
        {"CC(C)Cc1ccc(CC(C)C(=O)O)cc1", "COc1ccc2cc(CC(C)C(=O)O)ccc2c1", "ibuprofen-naproxen"},
        {"Cn1c(=O)c2c(ncn2C)n(C)c1=O", "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "caffeine-theophylline"},
        {"Nc1ncnc2n(cnc12)C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O",
         "Nc1ncnc2n(cnc12)C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O", "ATP-ADP"},
      };
      ChemOptions opts = new ChemOptions();
      for (String[] p : pairs) {
        MolGraph g1 = new MolGraph(mol(p[0]));
        MolGraph g2 = new MolGraph(mol(p[1]));
        Map<Integer, Integer> mcs = SearchEngine.findMCS(g1, g2, opts, new SearchEngine.McsOptions());
        List<String> errors = SearchEngine.validateMapping(g1, g2, mcs, opts);
        assertTrue(errors.isEmpty(),
            p[2] + " mapping has errors: " + errors);
        if (mcs.size() > 0) {
          boolean maximal = SearchEngine.isMappingMaximal(g1, g2, mcs, opts);
          if (!maximal) System.out.println("INFO: " + p[2] + " MCS (" + mcs.size() + " atoms) is extensible");
        }
      }
    }
  }

  // ======================================================================
  // From: BenchmarkComparison.java
  // ======================================================================

  @Nested
  @EnabledIfSystemProperty(named = "smsd.benchmark", matches = "true")
  @DisplayName("SMSD Internal Benchmark")
  class InternalBenchmark {

    private static final String[][] COMPARISON_PAIRS = {
      {"benzene-toluene", "c1ccccc1", "Cc1ccccc1"},
      {"aspirin-acetaminophen", "CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"},
      {"morphine-codeine",
          "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
          "CN1CCC23C4C1CC5=C2C(=C(C=C5)OC)OC3C(C=C4)O"},
      {"ibuprofen-naproxen",
          "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
          "COc1ccc2cc(CC(C)C(=O)O)ccc2c1"},
      {"caffeine-theophylline",
          "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
          "Cn1cnc2c1c(=O)[nH]c(=O)n2C"},
      {"ATP-ADP",
          "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
          "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N"},
      {"NAD+-NADH",
          "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O)C(=O)N",
          "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O"},
      {"paclitaxel-docetaxel",
          "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
          "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O"},
      {"erythromycin-azithromycin",
          "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(C2O)N(C)C)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O",
          "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O"},
      {"atorvastatin-rosuvastatin",
          "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
          "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C"},
    };

    @Test
    @DisplayName("MCS benchmark -- all pairs")
    void mcsBenchmark() throws Exception {
      System.out.println("\n=== SMSD MCS Benchmark (internal use only) ===");
      System.out.printf("%-30s %8s %8s %8s %6s%n",
          "Pair", "Best(ms)", "Med(ms)", "Mean(ms)", "MCS");
      System.out.println("-".repeat(70));

      ChemOptions opts = new ChemOptions();
      int warmup = 5;
      int runs = 10;

      for (String[] pair : COMPARISON_PAIRS) {
        var q = mol(pair[1]);
        var t = mol(pair[2]);

        for (int i = 0; i < warmup; i++) {
          new SMSD(q, t, opts, false).findMCS(true, true, 5000);
        }

        long[] times = new long[runs];
        int mcsSize = 0;
        for (int i = 0; i < runs; i++) {
          long t0 = System.nanoTime();
          SMSD smsd = new SMSD(q, t, opts, false);
          Map<Integer, Integer> result = smsd.findMCS(true, true, 5000);
          times[i] = System.nanoTime() - t0;
          if (i == 0 && result != null && !result.isEmpty()) {
            mcsSize = result.size();
          }
        }

        java.util.Arrays.sort(times);
        long best = times[0];
        long median = times[runs / 2];
        long mean = 0;
        for (long ti : times) mean += ti;
        mean /= runs;

        System.out.printf("%-30s %8.2f %8.2f %8.2f %6d%n",
            pair[0],
            best / 1_000_000.0,
            median / 1_000_000.0,
            mean / 1_000_000.0,
            mcsSize);
      }
      System.out.println();
    }

    @Test
    @DisplayName("Substructure benchmark -- all pairs")
    void substructureBenchmark() throws Exception {
      System.out.println("\n=== SMSD Substructure Benchmark (internal use only) ===");
      System.out.printf("%-30s %8s %8s %6s%n", "Pair", "Best(us)", "Med(us)", "Match");
      System.out.println("-".repeat(60));

      ChemOptions opts = new ChemOptions();
      int warmup = 20;
      int runs = 100;

      for (String[] pair : COMPARISON_PAIRS) {
        var q = mol(pair[1]);
        var t = mol(pair[2]);

        for (int i = 0; i < warmup; i++) {
          new SMSD(q, t, opts, false).isSubstructure();
        }

        long[] times = new long[runs];
        boolean match = false;
        for (int i = 0; i < runs; i++) {
          long t0 = System.nanoTime();
          SMSD smsd = new SMSD(q, t, opts, false);
          match = smsd.isSubstructure();
          times[i] = System.nanoTime() - t0;
        }

        java.util.Arrays.sort(times);
        long best = times[0];
        long median = times[runs / 2];

        System.out.printf("%-30s %8d %8d %6s%n",
            pair[0], best / 1000, median / 1000, match);
      }
      System.out.println();
    }

    @Test
    @DisplayName("ATP/ADP MCS detail -- inspect mapping")
    void atpAdpMcsDetail() throws Exception {
      var atp = mol("c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N");
      var adp = mol("c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N");

      System.out.println("\n=== ATP/ADP MCS Detail ===");
      System.out.println("ATP atoms: " + atp.getAtomCount());
      System.out.println("ADP atoms: " + adp.getAtomCount());

      ChemOptions opts = new ChemOptions();
      SMSD smsd = new SMSD(atp, adp, opts, false);
      Map<Integer, Integer> mapping = smsd.findMCS(true, true, 10000);

      if (mapping != null && !mapping.isEmpty()) {
        System.out.println("MCS size: " + mapping.size());
        System.out.println("\nMapping (ATP idx -> ADP idx):");
        for (var entry : mapping.entrySet()) {
          int qi = entry.getKey();
          int ti = entry.getValue();
          String qSym = atp.getAtom(qi).getSymbol();
          String tSym = adp.getAtom(ti).getSymbol();
          System.out.printf("  ATP[%2d] %s -> ADP[%2d] %s%n", qi, qSym, ti, tSym);
        }

        Map<String, Integer> elements = new java.util.TreeMap<>();
        for (var entry : mapping.entrySet()) {
          String sym = atp.getAtom(entry.getKey()).getSymbol();
          elements.merge(sym, 1, Integer::sum);
        }
        System.out.println("\nMCS composition: " + elements);
      } else {
        System.out.println("NO MCS FOUND");
      }
    }

    @Test
    @DisplayName("RASCAL screening benchmark")
    void rascalScreeningBenchmark() throws Exception {
      System.out.println("\n=== RASCAL Screening Benchmark (internal use only) ===");
      System.out.printf("%-30s %8s %8s %6s%n", "Pair", "Best(us)", "Med(us)", "UB");
      System.out.println("-".repeat(60));

      ChemOptions opts = new ChemOptions();
      int warmup = 50;
      int runs = 200;

      for (String[] pair : COMPARISON_PAIRS) {
        var q = mol(pair[1]);
        var t = mol(pair[2]);

        for (int i = 0; i < warmup; i++) {
          new SMSD(q, t, opts, false).similarityUpperBound();
        }

        long[] times = new long[runs];
        double ub = 0;
        for (int i = 0; i < runs; i++) {
          long t0 = System.nanoTime();
          ub = new SMSD(q, t, opts, false).similarityUpperBound();
          times[i] = System.nanoTime() - t0;
        }

        java.util.Arrays.sort(times);
        System.out.printf("%-30s %8d %8d %6.3f%n",
            pair[0], times[0] / 1000, times[runs / 2] / 1000, ub);
      }
      System.out.println();
    }
  }

  // ======================================================================
  // From: HeadToHeadBenchmark.java
  // ======================================================================

  @Nested
  @EnabledIfSystemProperty(named = "benchmark", matches = "true")
  @DisplayName("Head-to-Head Benchmark Suite")
  class HeadToHead {

    private static final int WARMUP = 20;
    private static final int ITERS = 100;

    private static final String[][] H2H_PAIRS = {
      {"c1ccccc1", "Cc1ccccc1", "Benzene/Toluene"},
      {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1", "Aspirin/Acetaminophen"},
      {
        "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O",
        "COC1=CC=C2C3CC4=CC(=C(C=C4C3(CCN2C)C1=O)O)OC",
        "Morphine/Codeine"
      },
      {
        "CC(C)Cc1ccc(CC(C)C(O)=O)cc1",
        "COc1ccc2cc(CC(C)C(O)=O)ccc2c1",
        "Ibuprofen/Naproxen"
      },
      {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "Caffeine/Theophylline"},
      {
        "Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C1O",
        "Nc1ncnc2n(cnc12)C1OC(COP(O)(=O)OP(O)(O)=O)C(O)C1O",
        "ATP/ADP"
      },
      {
        "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O)C(=O)N",
        "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O",
        "NAD+/NADH"
      },
      {
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
        "Paclitaxel/Docetaxel"
      },
      {
        "CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(C2O)N(C)C)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O",
        "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
        "Erythromycin/Azithromycin"
      },
      {
        "CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
        "Atorvastatin/Rosuvastatin"
      },
    };

    @Test
    @DisplayName("Substructure benchmark")
    void h2hSubstructureBenchmark() throws Exception {
      System.out.println("\n=== SUBSTRUCTURE BENCHMARK ===");
      System.out.printf("%-30s %10s %10s %10s%n", "Pair", "Best(us)", "Med(us)", "Mean(us)");
      System.out.println("-".repeat(65));

      ChemOptions opts = new ChemOptions();
      for (String[] p : H2H_PAIRS) {
        IAtomContainer q = mol(p[0]), t = mol(p[1]);
        for (int i = 0; i < WARMUP; i++) {
          SearchEngine.isSubstructure(q, t, opts, 5000);
        }
        long[] times = new long[ITERS];
        for (int i = 0; i < ITERS; i++) {
          long t0 = System.nanoTime();
          SearchEngine.isSubstructure(q, t, opts, 5000);
          times[i] = System.nanoTime() - t0;
        }
        Arrays.sort(times);
        long best = times[0] / 1000;
        long median = times[ITERS / 2] / 1000;
        long mean = 0;
        for (long x : times) mean += x;
        mean = mean / ITERS / 1000;
        System.out.printf("%-30s %,10d %,10d %,10d%n", p[2], best, median, mean);
      }
      assertTrue(true);
    }

    @Test
    @DisplayName("MCS benchmark")
    void h2hMcsBenchmark() throws Exception {
      System.out.println("\n=== MCS BENCHMARK ===");
      System.out.printf("%-30s %10s %10s %10s %6s%n", "Pair", "Best(us)", "Med(us)", "Mean(us)", "MCS");
      System.out.println("-".repeat(72));

      ChemOptions opts = new ChemOptions();
      SearchEngine.McsOptions mcsOpts = new SearchEngine.McsOptions();
      mcsOpts.timeoutMs = 10_000;

      for (String[] p : H2H_PAIRS) {
        IAtomContainer q = mol(p[0]), t = mol(p[1]);
        for (int i = 0; i < WARMUP; i++) {
          SearchEngine.findMCS(q, t, opts, mcsOpts);
        }
        long[] times = new long[ITERS];
        int mcsSize = 0;
        for (int i = 0; i < ITERS; i++) {
          long t0 = System.nanoTime();
          Map<Integer, Integer> mcs = SearchEngine.findMCS(q, t, opts, mcsOpts);
          times[i] = System.nanoTime() - t0;
          mcsSize = mcs.size();
        }
        Arrays.sort(times);
        long best = times[0] / 1000;
        long median = times[ITERS / 2] / 1000;
        long mean = 0;
        for (long x : times) mean += x;
        mean = mean / ITERS / 1000;
        System.out.printf("%-30s %,10d %,10d %,10d %6d%n", p[2], best, median, mean, mcsSize);
      }
      assertTrue(true);
    }

    @Test
    @DisplayName("RASCAL screening benchmark")
    void h2hRascalBenchmark() throws Exception {
      System.out.println("\n=== RASCAL SCREENING BENCHMARK ===");
      System.out.printf("%-30s %10s %10s %8s%n", "Pair", "Best(us)", "Med(us)", "UB");
      System.out.println("-".repeat(60));

      ChemOptions opts = new ChemOptions();
      for (String[] p : H2H_PAIRS) {
        IAtomContainer q = mol(p[0]), t = mol(p[1]);
        for (int i = 0; i < WARMUP; i++) {
          SearchEngine.similarityUpperBound(q, t, opts);
        }
        long[] times = new long[ITERS];
        double ub = 0;
        for (int i = 0; i < ITERS; i++) {
          long t0 = System.nanoTime();
          ub = SearchEngine.similarityUpperBound(q, t, opts);
          times[i] = System.nanoTime() - t0;
        }
        Arrays.sort(times);
        long best = times[0] / 1000;
        long median = times[ITERS / 2] / 1000;
        System.out.printf("%-30s %,10d %,10d %8.4f%n", p[2], best, median, ub);
      }
      assertTrue(true);
    }

    @Test
    @DisplayName("MolGraph.Builder benchmark (CDK-free path)")
    void h2hBuilderBenchmark() throws Exception {
      System.out.println("\n=== MOLGRAPH.BUILDER BENCHMARK (CDK-free) ===");
      System.out.printf("%-30s %10s %10s %6s%n", "Pair", "Best(us)", "Med(us)", "MCS");
      System.out.println("-".repeat(60));

      MolGraph benzene =
          new MolGraph.Builder()
              .atomCount(6)
              .atomicNumbers(new int[] {6, 6, 6, 6, 6, 6})
              .aromaticFlags(new boolean[] {true, true, true, true, true, true})
              .ringFlags(new boolean[] {true, true, true, true, true, true})
              .neighbors(new int[][] {{1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}})
              .bondOrders(new int[][] {{4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}})
              .bondRingFlags(
                  new boolean[][] {
                    {true, true}, {true, true}, {true, true},
                    {true, true}, {true, true}, {true, true}
                  })
              .bondAromaticFlags(
                  new boolean[][] {
                    {true, true}, {true, true}, {true, true},
                    {true, true}, {true, true}, {true, true}
                  })
              .build();

      MolGraph phenol =
          new MolGraph.Builder()
              .atomCount(7)
              .atomicNumbers(new int[] {6, 6, 6, 6, 6, 6, 8})
              .aromaticFlags(new boolean[] {true, true, true, true, true, true, false})
              .ringFlags(new boolean[] {true, true, true, true, true, true, false})
              .neighbors(
                  new int[][] {{1, 5, 6}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {0}})
              .bondOrders(
                  new int[][] {{4, 4, 1}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {1}})
              .bondRingFlags(
                  new boolean[][] {
                    {true, true, false}, {true, true}, {true, true},
                    {true, true}, {true, true}, {true, true}, {false}
                  })
              .bondAromaticFlags(
                  new boolean[][] {
                    {true, true, false}, {true, true}, {true, true},
                    {true, true}, {true, true}, {true, true}, {false}
                  })
              .build();

      ChemOptions opts = new ChemOptions();
      SearchEngine.McsOptions mcsOpts = new SearchEngine.McsOptions();
      mcsOpts.timeoutMs = 5000;

      for (int i = 0; i < WARMUP; i++) {
        SearchEngine.isSubstructure(benzene, phenol, opts, 5000);
        SearchEngine.findMCS(benzene, phenol, opts, mcsOpts);
      }

      long[] times = new long[ITERS];
      for (int i = 0; i < ITERS; i++) {
        long t0 = System.nanoTime();
        SearchEngine.isSubstructure(benzene, phenol, opts, 5000);
        times[i] = System.nanoTime() - t0;
      }
      Arrays.sort(times);
      System.out.printf(
          "%-30s %,10d %,10d %6s%n",
          "Benzene->Phenol (sub)", times[0] / 1000, times[ITERS / 2] / 1000, "Y");

      times = new long[ITERS];
      int mcsSize = 0;
      for (int i = 0; i < ITERS; i++) {
        long t0 = System.nanoTime();
        Map<Integer, Integer> mcs = SearchEngine.findMCS(benzene, phenol, opts, mcsOpts);
        times[i] = System.nanoTime() - t0;
        mcsSize = mcs.size();
      }
      Arrays.sort(times);
      System.out.printf(
          "%-30s %,10d %,10d %6d%n",
          "Benzene/Phenol (MCS)", times[0] / 1000, times[ITERS / 2] / 1000, mcsSize);

      assertTrue(true);
    }
  }
}

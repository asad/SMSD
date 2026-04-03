/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * SMSD Pro 6.4.0 — ZINC20 Tautomer Over-Matching Benchmark
 * =========================================================
 * Paper Item 15 (Future Work): quantify how tautomer-aware MCS relaxation
 * causes over-matching (matching atoms that violate proton conservation).
 *
 * Since ZINC20 cannot be downloaded in CI, this benchmark uses the
 * tautomer section (lines 801-900) of diverse_molecules.txt as a proxy
 * for drug-like ZINC20 tautomeric pairs.
 *
 * For each consecutive pair of tautomers:
 *   1. Run MCS with ChemOptions.tautomerProfile() (tautomer-aware)
 *   2. Run MCS with default ChemOptions()          (strict)
 *   3. Call SearchEngine.validateTautomerConsistency() on the tautomer result
 *   4. Compute a TautConf score = avg(tautomerWeight) for mapped tautomeric atoms
 *
 * Compile and run:
 *   cd <project-root>
 *   mvn package -DskipTests
 *   javac -cp target/smsd-6.4.0.jar:target/dependency/* \
 *         benchmarks/benchmark_tautomer_zinc.java -d benchmarks/
 *   java  -cp target/smsd-6.4.0.jar:target/dependency/*:benchmarks/ \
 *         benchmark_tautomer_zinc benchmarks/diverse_molecules.txt
 */

import com.bioinception.smsd.core.*;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class benchmark_tautomer_zinc {

    // --- Configuration ---
    static final int TIMEOUT_MS = 10_000;
    static final int TAUT_SECTION_START = 801;  // 1-based molecule index
    static final int TAUT_SECTION_END   = 900;

    record Molecule(String smiles, String name, IAtomContainer mol) {}

    record PairResult(
        String nameA, String nameB, String smiA, String smiB,
        int tautMcsSize, int defaultMcsSize, int overMatchDelta,
        boolean protonConsistent, double tautConfScore,
        double tautTimeMs, double defaultTimeMs
    ) {}

    private static final SmilesParser SP =
        new SmilesParser(SilentChemObjectBuilder.getInstance());

    // ======================================================================
    // Load molecules — filter to tautomer section
    // ======================================================================
    static List<Molecule> loadTautomerMolecules(Path path) throws Exception {
        List<Molecule> mols = new ArrayList<>();
        int lineNum = 0;
        int molIdx  = 0;
        for (String line : Files.readAllLines(path)) {
            lineNum++;
            line = line.trim();
            if (line.isEmpty() || line.startsWith("#")) continue;
            molIdx++;
            if (molIdx < TAUT_SECTION_START || molIdx > TAUT_SECTION_END) continue;

            String[] parts = line.split("\t", 2);
            String smi  = parts[0].trim();
            String name = parts.length > 1 ? parts[1].trim() : "mol_" + molIdx;
            try {
                IAtomContainer mol = Standardiser.standardise(
                    SP.parseSmiles(smi), Standardiser.TautomerMode.NONE);
                mols.add(new Molecule(smi, name, mol));
            } catch (Exception e) {
                System.err.printf("  SKIP: %s [%s] -> %s%n", name, smi, e.getMessage());
            }
        }
        return mols;
    }

    // ======================================================================
    // Compute TautConf score for a mapping
    // ======================================================================
    static double computeTautConfScore(MolGraph g1, MolGraph g2, Map<Integer,Integer> mcs) {
        g1.ensureTautomerClasses();
        g2.ensureTautomerClasses();
        if (g1.tautomerClass == null || g2.tautomerClass == null) return 1.0;
        if (mcs == null || mcs.isEmpty()) return 1.0;

        double weightSum = 0.0;
        int tautAtomCount = 0;
        for (Map.Entry<Integer,Integer> e : mcs.entrySet()) {
            int qi = e.getKey(), ti = e.getValue();
            if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
            boolean isTaut = (g1.tautomerClass[qi] >= 0) || (g2.tautomerClass[ti] >= 0);
            if (isTaut) {
                double w1 = (g1.tautomerWeight != null && qi < g1.tautomerWeight.length)
                    ? g1.tautomerWeight[qi] : 1.0;
                double w2 = (g2.tautomerWeight != null && ti < g2.tautomerWeight.length)
                    ? g2.tautomerWeight[ti] : 1.0;
                weightSum += (w1 + w2) / 2.0;
                tautAtomCount++;
            }
        }
        return tautAtomCount > 0 ? weightSum / tautAtomCount : 1.0;
    }

    // ======================================================================
    // Run one pair: tautomer-aware vs default
    // ======================================================================
    static PairResult benchmarkPair(Molecule a, Molecule b) {
        ChemOptions tautOpts    = ChemOptions.tautomerProfile();
        ChemOptions defaultOpts = new ChemOptions();
        SearchEngine.McsOptions mcsOpts = new SearchEngine.McsOptions();
        mcsOpts.timeoutMillis = TIMEOUT_MS;

        // Tautomer-aware MCS
        Map<Integer,Integer> tautMcs = null;
        long t0 = System.nanoTime();
        try {
            tautMcs = SearchEngine.findMCS(a.mol(), b.mol(), tautOpts, mcsOpts);
        } catch (Exception e) {
            System.err.printf("  TAUT MCS error: %s vs %s -> %s%n",
                a.name(), b.name(), e.getMessage());
        }
        double tautTimeMs = (System.nanoTime() - t0) / 1_000_000.0;

        // Default MCS
        Map<Integer,Integer> defMcs = null;
        t0 = System.nanoTime();
        try {
            defMcs = SearchEngine.findMCS(a.mol(), b.mol(), defaultOpts, mcsOpts);
        } catch (Exception e) {
            System.err.printf("  DEF MCS error: %s vs %s -> %s%n",
                a.name(), b.name(), e.getMessage());
        }
        double defTimeMs = (System.nanoTime() - t0) / 1_000_000.0;

        int tautSize = tautMcs != null ? tautMcs.size() : 0;
        int defSize  = defMcs  != null ? defMcs.size()  : 0;
        int overMatchDelta = tautSize - defSize;

        // Validate proton consistency
        boolean consistent = true;
        if (tautMcs != null && !tautMcs.isEmpty()) {
            consistent = SearchEngine.validateTautomerConsistency(a.mol(), b.mol(), tautMcs);
        }

        // TautConf score
        double tautConf = 1.0;
        if (tautMcs != null && !tautMcs.isEmpty()) {
            tautConf = computeTautConfScore(new MolGraph(a.mol()), new MolGraph(b.mol()), tautMcs);
        }

        return new PairResult(a.name(), b.name(), a.smiles(), b.smiles(),
            tautSize, defSize, overMatchDelta,
            consistent, tautConf, tautTimeMs, defTimeMs);
    }

    // ======================================================================
    // Main
    // ======================================================================
    public static void main(String[] args) throws Exception {
        Path molPath = args.length > 0 ? Path.of(args[0]) : Path.of("benchmarks/diverse_molecules.txt");

        System.err.printf("Loading tautomer molecules from %s (lines %d-%d) ...%n",
            molPath, TAUT_SECTION_START, TAUT_SECTION_END);
        List<Molecule> mols = loadTautomerMolecules(molPath);
        System.err.printf("  %d valid tautomer molecules loaded%n", mols.size());

        if (mols.size() < 2) {
            System.err.println("ERROR: Need at least 2 molecules. Aborting.");
            System.exit(1);
        }

        // Form consecutive pairs (keto/enol, amide/iminol, etc.)
        List<int[]> pairIndices = new ArrayList<>();
        for (int i = 0; i + 1 < mols.size(); i += 2) {
            pairIndices.add(new int[]{i, i + 1});
        }
        System.err.printf("  %d tautomer pairs formed%n", pairIndices.size());

        // Warmup JIT
        System.err.println("JVM warmup ...");
        if (mols.size() >= 2) {
            ChemOptions warmOpts = ChemOptions.tautomerProfile();
            SearchEngine.McsOptions warmMcs = new SearchEngine.McsOptions();
            warmMcs.timeoutMillis = 2000;
            for (int w = 0; w < 5; w++) {
                try {
                    SearchEngine.findMCS(mols.get(0).mol(), mols.get(1).mol(), warmOpts, warmMcs);
                } catch (Exception ignored) {}
            }
        }

        // Run benchmark
        System.err.println("Running tautomer over-matching benchmark ...");
        List<PairResult> results = new ArrayList<>();
        for (int i = 0; i < pairIndices.size(); i++) {
            int[] pi = pairIndices.get(i);
            PairResult pr = benchmarkPair(mols.get(pi[0]), mols.get(pi[1]));
            results.add(pr);
            System.err.printf("  [%2d/%d] %-40s taut=%2d def=%2d delta=%+d consistent=%s tautConf=%.3f%n",
                i + 1, pairIndices.size(),
                pr.nameA() + " / " + pr.nameB(),
                pr.tautMcsSize(), pr.defaultMcsSize(), pr.overMatchDelta(),
                pr.protonConsistent() ? "PASS" : "FAIL", pr.tautConfScore());
        }

        // ================================================================
        // Report
        // ================================================================
        System.out.println();
        System.out.println("=".repeat(78));
        System.out.println("SMSD Pro 6.4.0 — ZINC20 Tautomer Over-Matching Benchmark Results");
        System.out.println("=".repeat(78));
        System.out.printf("Pairs tested:              %d%n", results.size());

        // Header
        System.out.printf("%n%-42s %5s %5s %6s %6s %8s%n",
            "Pair", "Taut", "Def", "Delta", "Valid", "TautConf");
        System.out.println("-".repeat(78));

        int overMatchCount = 0;
        int failCount      = 0;
        double totalTautConf = 0.0;
        double totalTautTime = 0.0, totalDefTime = 0.0;

        for (PairResult pr : results) {
            String pair = pr.nameA() + " / " + pr.nameB();
            if (pair.length() > 42) pair = pair.substring(0, 39) + "...";
            System.out.printf("%-42s %5d %5d %+6d %6s %8.3f%n",
                pair, pr.tautMcsSize(), pr.defaultMcsSize(), pr.overMatchDelta(),
                pr.protonConsistent() ? "PASS" : "FAIL", pr.tautConfScore());

            if (pr.overMatchDelta() > 0) overMatchCount++;
            if (!pr.protonConsistent()) failCount++;
            totalTautConf += pr.tautConfScore();
            totalTautTime += pr.tautTimeMs();
            totalDefTime  += pr.defaultTimeMs();
        }

        System.out.println("-".repeat(78));
        int n = results.size();
        double avgTautConf = n > 0 ? totalTautConf / n : 0.0;
        double passRate    = n > 0 ? 100.0 * (n - failCount) / n : 0.0;
        double overMatchPct= n > 0 ? 100.0 * overMatchCount / n : 0.0;

        System.out.printf("%nSummary Statistics:%n");
        System.out.printf("  Total pairs:                     %d%n", n);
        System.out.printf("  Pairs with over-matching:        %d (%.1f%%)%n", overMatchCount, overMatchPct);
        System.out.printf("  Proton consistency PASS rate:    %.1f%% (%d/%d)%n", passRate, n - failCount, n);
        System.out.printf("  Proton consistency FAIL count:   %d%n", failCount);
        System.out.printf("  Average TautConf score:          %.4f%n", avgTautConf);
        System.out.printf("  Total tautomer MCS time:         %.1f ms%n", totalTautTime);
        System.out.printf("  Total default MCS time:          %.1f ms%n", totalDefTime);
        System.out.printf("  Avg tautomer MCS time per pair:  %.2f ms%n", n > 0 ? totalTautTime / n : 0);
        System.out.printf("  Avg default MCS time per pair:   %.2f ms%n", n > 0 ? totalDefTime / n : 0);

        // Over-matching details
        List<PairResult> overMatched = results.stream()
            .filter(r -> r.overMatchDelta() > 0)
            .sorted(Comparator.comparingInt(PairResult::overMatchDelta).reversed())
            .collect(Collectors.toList());

        if (!overMatched.isEmpty()) {
            System.out.printf("%nOver-matched pairs (tautomer MCS > default MCS):%n");
            for (PairResult pr : overMatched) {
                System.out.printf("  %s / %s: delta=%+d (taut=%d, def=%d) consistent=%s tautConf=%.3f%n",
                    pr.nameA(), pr.nameB(), pr.overMatchDelta(),
                    pr.tautMcsSize(), pr.defaultMcsSize(),
                    pr.protonConsistent() ? "PASS" : "FAIL", pr.tautConfScore());
            }
        }

        // False positives: over-matched AND proton-inconsistent
        List<PairResult> falsePositives = results.stream()
            .filter(r -> r.overMatchDelta() > 0 && !r.protonConsistent())
            .collect(Collectors.toList());
        System.out.printf("%nFalse positives (over-matched + proton violation): %d%n", falsePositives.size());
        for (PairResult pr : falsePositives) {
            System.out.printf("  %s / %s: delta=%+d tautConf=%.3f%n",
                pr.nameA(), pr.nameB(), pr.overMatchDelta(), pr.tautConfScore());
        }

        System.out.println("=".repeat(78));
    }
}

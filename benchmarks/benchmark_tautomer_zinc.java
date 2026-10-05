/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * SMSD tautomer feature diagnostics on the curated molecule pool
 * =========================================================
 * Compare returned sizes and the engine's proton-consistency diagnostics.
 * Validator failures are not independent proof of chemical false positives.
 *
 * This benchmark uses the labelled tautomer section of diverse_molecules.txt.
 * The historical filename does not identify an original ZINC20 corpus.
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
 *   javac -cp target/smsd-7.2.0-jar-with-dependencies.jar \
 *         benchmarks/benchmark_tautomer_zinc.java -d build/local-benchmarks/java
 *   java  -cp target/smsd-7.2.0-jar-with-dependencies.jar:build/local-benchmarks/java \
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
    static final int TIMEOUT_MS = Integer.getInteger("smsd.benchmark.timeoutMs", 10_000);

    record Molecule(String smiles, String name, IAtomContainer mol) {}

    record PairResult(
        String nameA, String nameB, String smiA, String smiB,
        int tautMCSSize, int defaultMCSSize, int overMatchDelta,
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
        boolean inTautomerSection = false;
        for (String line : Files.readAllLines(path)) {
            lineNum++;
            line = line.trim();
            if (line.startsWith("# SECTION ")) {
                inTautomerSection = line.startsWith("# SECTION 7:");
                continue;
            }
            if (line.isEmpty() || line.startsWith("#")) continue;
            molIdx++;
            if (!inTautomerSection) continue;

            String[] parts = line.split("\t", 2);
            String smi  = parts[0].trim();
            String name = parts.length > 1 ? parts[1].trim() : "mol_" + molIdx;
            try {
                IAtomContainer mol = Standardiser.standardise(
                    SP.parseSmiles(smi), Standardiser.TautomerMode.NONE);
                mols.add(new Molecule(smi, name, mol));
            } catch (Exception e) {
                System.err.printf("  SKIP: %s [%s] -> %s%n", name, smi, e.getMessage());
                mols.add(new Molecule(smi, name, null));
            }
        }
        return mols;
    }

    // ======================================================================
    // Compute TautConf score for a mapping
    // ======================================================================
    static double computeTautConfScore(MolGraph g1, MolGraph g2, Map<Integer,Integer> mcs) {
        if (mcs == null || mcs.isEmpty()) return 1.0;
        // The public confidence API initializes the exported tautomer annotations.
        // Keep this benchmark's arithmetic mean below rather than its geometric mean.
        SearchEngine.computeTautomerConfidence(g1, g2, mcs);
        if (g1.tautomerClass == null || g2.tautomerClass == null) return 1.0;

        double weightSum = 0.0;
        int tautAtomCount = 0;
        for (Map.Entry<Integer,Integer> e : mcs.entrySet()) {
            int qi = e.getKey(), ti = e.getValue();
            if (qi < 0 || qi >= g1.atomCount() || ti < 0 || ti >= g2.atomCount()) continue;
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
        SearchEngine.MCSOptions mcsOpts = new SearchEngine.MCSOptions();
        mcsOpts.timeoutMs = TIMEOUT_MS;

        // Tautomer-aware MCS
        Map<Integer,Integer> tautMCS = null;
        long t0 = System.nanoTime();
        try {
            tautMCS = SearchEngine.findMCS(a.mol(), b.mol(), tautOpts, mcsOpts);
        } catch (Exception e) {
            System.err.printf("  TAUT MCS error: %s vs %s -> %s%n",
                a.name(), b.name(), e.getMessage());
        }
        double tautTimeMs = (System.nanoTime() - t0) / 1_000_000.0;

        // Default MCS
        Map<Integer,Integer> defMCS = null;
        t0 = System.nanoTime();
        try {
            defMCS = SearchEngine.findMCS(a.mol(), b.mol(), defaultOpts, mcsOpts);
        } catch (Exception e) {
            System.err.printf("  DEF MCS error: %s vs %s -> %s%n",
                a.name(), b.name(), e.getMessage());
        }
        double defTimeMs = (System.nanoTime() - t0) / 1_000_000.0;

        int tautSize = tautMCS != null ? tautMCS.size() : 0;
        int defSize  = defMCS  != null ? defMCS.size()  : 0;
        int overMatchDelta = tautSize - defSize;

        // Validate proton consistency
        boolean consistent = true;
        if (tautMCS != null && !tautMCS.isEmpty()) {
            consistent = SearchEngine.validateTautomerConsistency(a.mol(), b.mol(), tautMCS);
        }

        // TautConf score
        double tautConf = 1.0;
        if (tautMCS != null && !tautMCS.isEmpty()) {
            tautConf = computeTautConfScore(new MolGraph(a.mol()), new MolGraph(b.mol()), tautMCS);
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

        System.err.printf("Loading labelled tautomer section from %s ...%n", molPath);
        List<Molecule> mols = loadTautomerMolecules(molPath);
        System.err.printf("  %d tautomer records loaded%n", mols.size());

        if (mols.size() < 2) {
            System.err.println("ERROR: Need at least 2 molecules. Aborting.");
            System.exit(1);
        }

        // Form consecutive pairs (keto/enol, amide/iminol, etc.)
        List<int[]> pairIndices = new ArrayList<>();
        for (int i = 0; i + 1 < mols.size(); i += 2) {
            if (mols.get(i).mol() != null && mols.get(i+1).mol() != null) pairIndices.add(new int[]{i, i + 1});
            else System.err.printf("  SKIP PAIR: %s / %s%n", mols.get(i).name(), mols.get(i+1).name());
        }
        System.err.printf("  %d tautomer pairs formed%n", pairIndices.size());

        // Warmup JIT
        System.err.println("JVM warmup ...");
        if (!pairIndices.isEmpty()) {
            ChemOptions warmOpts = ChemOptions.tautomerProfile();
            SearchEngine.MCSOptions warmMCS = new SearchEngine.MCSOptions();
            warmMCS.timeoutMs = 2000;
            for (int w = 0; w < 5; w++) {
                try {
                    int[] first = pairIndices.get(0);
                    SearchEngine.findMCS(mols.get(first[0]).mol(), mols.get(first[1]).mol(), warmOpts, warmMCS);
                } catch (Exception ignored) {}
            }
        }

        // Run benchmark
        System.err.println("Running tautomer feature diagnostics ...");
        List<PairResult> results = new ArrayList<>();
        for (int i = 0; i < pairIndices.size(); i++) {
            int[] pi = pairIndices.get(i);
            PairResult pr = benchmarkPair(mols.get(pi[0]), mols.get(pi[1]));
            results.add(pr);
            System.err.printf("  [%2d/%d] %-40s taut=%2d def=%2d delta=%+d consistent=%s tautConf=%.3f%n",
                i + 1, pairIndices.size(),
                pr.nameA() + " / " + pr.nameB(),
                pr.tautMCSSize(), pr.defaultMCSSize(), pr.overMatchDelta(),
                pr.protonConsistent() ? "PASS" : "FAIL", pr.tautConfScore());
        }

        // ================================================================
        // Report
        // ================================================================
        System.out.println();
        System.out.println("=".repeat(78));
        System.out.println("SMSD tautomer feature diagnostics on the curated molecule pool");
        System.out.println("=".repeat(78));
        System.out.printf("Pairs tested:              %d%n", results.size());

        // Header
        System.out.printf("%n%-42s %5s %5s %6s %6s %8s%n",
            "Pair", "Taut", "Def", "Delta", "Proton", "TautConf");
        System.out.println("-".repeat(78));

        int overMatchCount = 0;
        int failCount      = 0;
        double totalTautConf = 0.0;
        double totalTautTime = 0.0, totalDefTime = 0.0;

        for (PairResult pr : results) {
            String pair = pr.nameA() + " / " + pr.nameB();
            if (pair.length() > 42) pair = pair.substring(0, 39) + "...";
            System.out.printf("%-42s %5d %5d %+6d %6s %8.3f%n",
                pair, pr.tautMCSSize(), pr.defaultMCSSize(), pr.overMatchDelta(),
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
        System.out.printf("  Pairs with positive atom delta:        %d (%.1f%%)%n", overMatchCount, overMatchPct);
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
            System.out.printf("%nPositive atom deltas (tautomer MCS > default MCS):%n");
            for (PairResult pr : overMatched) {
                System.out.printf("  %s / %s: delta=%+d (taut=%d, def=%d) consistent=%s tautConf=%.3f%n",
                    pr.nameA(), pr.nameB(), pr.overMatchDelta(),
                    pr.tautMCSSize(), pr.defaultMCSSize(),
                    pr.protonConsistent() ? "PASS" : "FAIL", pr.tautConfScore());
            }
        }

        // Positive size deltas with the engine's proton-consistency failure.
        List<PairResult> falsePositives = results.stream()
            .filter(r -> r.overMatchDelta() > 0 && !r.protonConsistent())
            .collect(Collectors.toList());
        System.out.printf("%nPositive atom deltas with proton-consistency failures: %d%n", falsePositives.size());
        for (PairResult pr : falsePositives) {
            System.out.printf("  %s / %s: delta=%+d tautConf=%.3f%n",
                pr.nameA(), pr.nameB(), pr.overMatchDelta(), pr.tautConfScore());
        }

        System.out.println("=".repeat(78));
    }
}

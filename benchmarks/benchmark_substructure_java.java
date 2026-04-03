/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * Substructure Search Benchmark: SMSD vs CDK (Java-to-Java)
 * ==========================================================
 *
 * Fair like-for-like comparison of substructure search performance.
 * Both SMSD and CDK run in the same JVM process with identical molecules.
 *
 * Compile and run:
 *   cd <project-root>
 *   mvn package -DskipTests
 *   javac -cp target/smsd-*.jar:target/dependency/* \
 *         benchmarks/benchmark_substructure_java.java -d benchmarks/
 *   java -cp target/smsd-*.jar:target/dependency/*:benchmarks/ \
 *         benchmark_substructure_java
 *
 * Output: results_substructure.tsv
 */

import com.bioinception.smsd.core.*;
import com.bioinception.smsd.core.MolGraph;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.isomorphism.Pattern;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.*;
import java.nio.file.*;
import java.util.*;

public class benchmark_substructure_java {

    // --- Configuration ---
    static final int WARMUP = 20;
    static final int ITERS  = 100;
    static final long TIMEOUT_MS = 5_000;
    static final String DEFAULT_PAIRS_FILE = "benchmarks/substructure_pairs.tsv";

    static final SmilesParser SP =
        new SmilesParser(SilentChemObjectBuilder.getInstance());

    static class TestPair {
        final String querySmi, targetSmi, label;
        TestPair(String q, String t, String l) { querySmi = q; targetSmi = t; label = l; }
    }

    static List<TestPair> loadPairs(Path path) throws IOException {
        List<TestPair> pairs = new ArrayList<>();
        for (String line : Files.readAllLines(path)) {
            line = line.trim();
            if (line.isEmpty() || line.startsWith("#")) continue;
            String[] cols = line.split("\t", 3);
            if (cols.length < 3) {
                System.err.println("  SKIP (bad format): " + line);
                continue;
            }
            pairs.add(new TestPair(cols[0].trim(), cols[1].trim(), cols[2].trim()));
        }
        return pairs;
    }

    static IAtomContainer prepareMol(IAtomContainer mol, Aromaticity arom)
            throws Exception {
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
        CDKHydrogenAdder.getInstance(
            SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(mol);
        arom.apply(mol);
        return mol;
    }

    // =====================================================================
    public static void main(String[] args) throws Exception {
        String pairsFile = args.length > 0 ? args[0] : DEFAULT_PAIRS_FILE;
        List<TestPair> pairs = loadPairs(Path.of(pairsFile));
        System.out.printf("Substructure Benchmark: SMSD vs CDK (Java-to-Java)%n");
        System.out.printf("Loaded %d pairs from %s%n", pairs.size(), pairsFile);
        System.out.println("WARMUP=" + WARMUP + " ITERS=" + ITERS);
        System.out.println("================================================\n");

        Aromaticity arom = new Aromaticity(
            ElectronDonation.daylight(), Cycles.or(Cycles.all(), Cycles.relevant()));

        // Pre-parse all molecules
        IAtomContainer[][] mols = new IAtomContainer[pairs.size()][2];
        for (int i = 0; i < pairs.size(); i++) {
            TestPair p = pairs.get(i);
            try {
                mols[i][0] = prepareMol(SP.parseSmiles(p.querySmi), arom);
                mols[i][1] = prepareMol(SP.parseSmiles(p.targetSmi), arom);
            } catch (Exception e) {
                System.err.printf("  SKIP parse error [%s]: %s%n", p.label, e.getMessage());
                mols[i] = null;
            }
        }

        // TSV output
        PrintWriter tsv = new PrintWriter(
            new FileWriter("benchmarks/results_substructure.tsv"));
        tsv.println("Pair\tSMSD_us\tSMSD_cached_us\tCDK_us\tSpeedup\tCached_Speedup\tSMSD_match\tCDK_match");

        ChemOptions opts = new ChemOptions();

        for (int i = 0; i < pairs.size(); i++) {
            if (mols[i] == null) continue; // skipped due to parse error
            String label = pairs.get(i).label;
            IAtomContainer query  = mols[i][0];
            IAtomContainer target = mols[i][1];

            // Pre-build MolGraph for cached path
            MolGraph qGraph = new MolGraph(query);
            MolGraph tGraph = new MolGraph(target);

            // --- SMSD (IAtomContainer API, rebuilds MolGraph each call) ---
            for (int w = 0; w < WARMUP; w++)
                SearchEngine.isSubstructure(query, target, opts, TIMEOUT_MS);
            long[] smsdTimes = new long[ITERS];
            boolean smsdMatch = false;
            for (int it = 0; it < ITERS; it++) {
                long t0 = System.nanoTime();
                smsdMatch = SearchEngine.isSubstructure(query, target, opts, TIMEOUT_MS);
                smsdTimes[it] = System.nanoTime() - t0;
            }
            Arrays.sort(smsdTimes);
            double smsdMedianUs = smsdTimes[ITERS / 2] / 1000.0;

            // --- SMSD (MolGraph API, cached graphs — pure search time) ---
            for (int w = 0; w < WARMUP; w++)
                SearchEngine.isSubstructure(qGraph, tGraph, opts, TIMEOUT_MS);
            long[] cachedTimes = new long[ITERS];
            for (int it = 0; it < ITERS; it++) {
                long t0 = System.nanoTime();
                SearchEngine.isSubstructure(qGraph, tGraph, opts, TIMEOUT_MS);
                cachedTimes[it] = System.nanoTime() - t0;
            }
            Arrays.sort(cachedTimes);
            double cachedMedianUs = cachedTimes[ITERS / 2] / 1000.0;

            // --- CDK DfPattern ---
            for (int w = 0; w < WARMUP; w++) {
                Pattern pattern = DfPattern.findSubstructure(query);
                pattern.matches(target);
            }
            long[] cdkTimes = new long[ITERS];
            boolean cdkMatch = false;
            for (int it = 0; it < ITERS; it++) {
                long t0 = System.nanoTime();
                Pattern pattern = DfPattern.findSubstructure(query);
                cdkMatch = pattern.matches(target);
                cdkTimes[it] = System.nanoTime() - t0;
            }
            Arrays.sort(cdkTimes);
            double cdkMedianUs = cdkTimes[ITERS / 2] / 1000.0;

            double speedup = cdkMedianUs / smsdMedianUs;
            double cachedSpeedup = cdkMedianUs / cachedMedianUs;

            System.out.printf("%-40s  SMSD: %7.1f  cached: %7.1f  CDK: %7.1f us  %.1fx / %.1fx  [%s vs %s]%n",
                label, smsdMedianUs, cachedMedianUs, cdkMedianUs, speedup, cachedSpeedup,
                smsdMatch ? "Y" : "N", cdkMatch ? "Y" : "N");

            tsv.printf("%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\t%s%n",
                label, smsdMedianUs, cachedMedianUs, cdkMedianUs, speedup, cachedSpeedup,
                smsdMatch ? "Y" : "N", cdkMatch ? "Y" : "N");
        }

        tsv.close();
        System.out.println("\nResults written to benchmarks/results_substructure.tsv");
    }
}

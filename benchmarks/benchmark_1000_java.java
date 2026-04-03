/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * SMSD 1000-Molecule Benchmark — Java Side
 * ==========================================
 *
 * Reads diverse_molecules.txt, generates 500 random + 500 systematic pairs,
 * and benchmarks SMSD MCS on each pair (5 rounds, 10s timeout).
 *
 * Compile and run:
 *   cd <project-root>
 *   mvn package -DskipTests
 *   javac -cp target/smsd-6.0.0.jar:target/dependency/* benchmarks/benchmark_1000_java.java -d benchmarks/
 *   java -cp target/smsd-6.0.0.jar:target/dependency/*:benchmarks/ benchmark_1000_java benchmarks/diverse_molecules.txt
 *
 * Output:
 *   benchmark_smsd_results.tsv — per-pair timing and MCS sizes
 *   benchmark_smsd_summary.txt — aggregate statistics
 */

import com.bioinception.smsd.core.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.stream.*;

public class benchmark_1000_java {

    // --- Configuration ---
    static final int ROUNDS = 5;
    static final int TIMEOUT_MS = 10_000;
    static final int N_RANDOM_PAIRS = 500;
    static final int N_SYSTEMATIC_PAIRS = 500;
    static final long SEED = 42L;

    // --- Data structures ---
    record Molecule(String smiles, String name, IAtomContainer mol) {}
    record Pair(int idxA, int idxB, String nameA, String nameB, String smiA, String smiB, String pairType) {}
    record BenchResult(double medianTimeMs, double meanTimeMs, int mcsSize, int completed, int rounds) {}

    private static final SmilesParser SP = new SmilesParser(SilentChemObjectBuilder.getInstance());

    // ======================================================================
    // Molecule loading
    // ======================================================================
    static List<Molecule> loadMolecules(Path path) throws Exception {
        List<Molecule> mols = new ArrayList<>();
        int invalidCount = 0;
        for (String line : Files.readAllLines(path)) {
            line = line.trim();
            if (line.isEmpty() || line.startsWith("#")) continue;
            String[] parts = line.split("\t", 2);
            String smi = parts[0].trim();
            String name = parts.length > 1 ? parts[1].trim() : "mol_" + mols.size();
            try {
                IAtomContainer mol = Standardiser.standardise(
                    SP.parseSmiles(smi), Standardiser.TautomerMode.NONE);
                mols.add(new Molecule(smi, name, mol));
            } catch (Exception e) {
                invalidCount++;
                System.err.printf("  INVALID: %s -> %s (%s)%n", name, smi, e.getMessage());
            }
        }
        if (invalidCount > 0) {
            System.err.printf("  %d invalid SMILES skipped out of %d%n",
                invalidCount, mols.size() + invalidCount);
        }
        return mols;
    }

    // ======================================================================
    // Pair generation (mirrors Python logic)
    // ======================================================================
    static List<Pair> generatePairs(List<Molecule> mols) {
        Random rng = new Random(SEED);
        int n = mols.size();
        List<Pair> pairs = new ArrayList<>();
        Set<Long> seen = new HashSet<>();

        // --- Random pairs ---
        int attempts = 0;
        while (pairs.size() < N_RANDOM_PAIRS && attempts < N_RANDOM_PAIRS * 20) {
            int i = rng.nextInt(n);
            int j = rng.nextInt(n);
            if (i == j) { attempts++; continue; }
            long key = Math.min(i, j) * 100_000L + Math.max(i, j);
            if (seen.contains(key)) { attempts++; continue; }
            seen.add(key);
            pairs.add(new Pair(i, j, mols.get(i).name(), mols.get(j).name(),
                               mols.get(i).smiles(), mols.get(j).smiles(), "random"));
        }

        // --- Systematic: size-matched pairs ---
        Integer[] sizeOrder = new Integer[n];
        for (int i = 0; i < n; i++) sizeOrder[i] = i;
        Arrays.sort(sizeOrder, Comparator.comparingInt(x -> mols.get(x).smiles().length()));

        int sysCount = 0;
        for (int k = 0; k < sizeOrder.length - 1 && sysCount < N_SYSTEMATIC_PAIRS; k += 2) {
            int i = sizeOrder[k], j = sizeOrder[k + 1];
            long key = Math.min(i, j) * 100_000L + Math.max(i, j);
            if (seen.contains(key)) continue;
            seen.add(key);
            pairs.add(new Pair(i, j, mols.get(i).name(), mols.get(j).name(),
                               mols.get(i).smiles(), mols.get(j).smiles(), "size_matched"));
            sysCount++;
        }

        return pairs;
    }

    // ======================================================================
    // SMSD MCS benchmark
    // ======================================================================
    static BenchResult benchmarkPair(Molecule a, Molecule b) {
        ChemOptions opts = new ChemOptions();
        SearchEngine.McsOptions mcsOpts = new SearchEngine.McsOptions();
        mcsOpts.timeoutMillis = TIMEOUT_MS;

        double[] times = new double[ROUNDS];
        int[] mcsSizes = new int[ROUNDS];
        int completed = 0;

        // Warmup (2 rounds)
        for (int w = 0; w < 2; w++) {
            try {
                SearchEngine.findMCS(a.mol(), b.mol(), opts, mcsOpts);
            } catch (Exception e) {
                // ignore warmup errors
            }
        }

        // Timed rounds
        for (int r = 0; r < ROUNDS; r++) {
            long t0 = System.nanoTime();
            try {
                Map<Integer, Integer> mcs = SearchEngine.findMCS(a.mol(), b.mol(), opts, mcsOpts);
                long elapsed = System.nanoTime() - t0;
                times[r] = elapsed / 1_000_000.0; // ms
                mcsSizes[r] = mcs != null ? mcs.size() : 0;
                completed++;
            } catch (Exception e) {
                long elapsed = System.nanoTime() - t0;
                times[r] = elapsed / 1_000_000.0;
                mcsSizes[r] = -1;
            }
        }

        Arrays.sort(times);
        double medianTime = times[ROUNDS / 2];
        double meanTime = 0;
        for (double t : times) meanTime += t;
        meanTime /= ROUNDS;

        int maxMcs = Arrays.stream(mcsSizes).max().orElse(-1);

        return new BenchResult(medianTime, meanTime, maxMcs, completed, ROUNDS);
    }

    // ======================================================================
    // Output
    // ======================================================================
    static void writeResults(List<Pair> pairs, Map<String, BenchResult> results, Path outPath) throws IOException {
        try (PrintWriter pw = new PrintWriter(Files.newBufferedWriter(outPath))) {
            pw.println("pair_id\tpair_type\tname_a\tname_b\tsmi_a\tsmi_b\t" +
                       "smsd_median_ms\tsmsd_mean_ms\tsmsd_mcs_size\tsmsd_completed");
            for (int i = 0; i < pairs.size(); i++) {
                Pair p = pairs.get(i);
                String key = p.nameA() + "__vs__" + p.nameB();
                BenchResult br = results.getOrDefault(key, new BenchResult(-1, -1, -1, 0, ROUNDS));
                pw.printf("%d\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%d\t%d/%d%n",
                    i + 1, p.pairType(), p.nameA(), p.nameB(), p.smiA(), p.smiB(),
                    br.medianTimeMs(), br.meanTimeMs(), br.mcsSize(), br.completed(), br.rounds());
            }
        }
    }

    static void writeSummary(List<Pair> pairs, Map<String, BenchResult> results,
                             int totalMols, Path outPath) throws IOException {
        double[] allMedians = results.values().stream()
            .mapToDouble(BenchResult::medianTimeMs)
            .filter(t -> t >= 0)
            .toArray();
        int[] allMcs = results.values().stream()
            .mapToInt(BenchResult::mcsSize)
            .filter(s -> s >= 0)
            .toArray();
        long completedAll = results.values().stream()
            .filter(r -> r.completed() == r.rounds())
            .count();

        Arrays.sort(allMedians);
        Arrays.sort(allMcs);

        StringBuilder sb = new StringBuilder();
        sb.append("=".repeat(70)).append("\n");
        sb.append("SMSD 1000-Molecule Benchmark — Java Results\n");
        sb.append("=".repeat(70)).append("\n");
        sb.append(String.format("Total molecules loaded: %d%n", totalMols));
        sb.append(String.format("Total pairs tested:     %d%n", pairs.size()));
        sb.append(String.format("Rounds per pair:        %d%n", ROUNDS));
        sb.append(String.format("Timeout per pair:       %d ms%n", TIMEOUT_MS));
        sb.append("\n");

        if (allMedians.length > 0) {
            double medOfMed = allMedians[allMedians.length / 2];
            double meanOfMed = Arrays.stream(allMedians).average().orElse(-1);
            double minTime = allMedians[0];
            double maxTime = allMedians[allMedians.length - 1];

            sb.append("--- SMSD Timing ---\n");
            sb.append(String.format("  Median of medians:  %.3f ms%n", medOfMed));
            sb.append(String.format("  Mean of medians:    %.3f ms%n", meanOfMed));
            sb.append(String.format("  Fastest pair:       %.3f ms%n", minTime));
            sb.append(String.format("  Slowest pair:       %.3f ms%n", maxTime));
            sb.append(String.format("  Completion rate:    %d/%d (%.1f%%)%n",
                completedAll, results.size(), 100.0 * completedAll / results.size()));
        }
        if (allMcs.length > 0) {
            double medMcs = allMcs[allMcs.length / 2];
            double meanMcs = Arrays.stream(allMcs).average().orElse(-1);
            sb.append(String.format("  Median MCS size:    %.1f%n", medMcs));
            sb.append(String.format("  Mean MCS size:      %.1f%n", meanMcs));
        }
        sb.append("\n");

        // Per-section breakdown
        String[][] sections = {
            {"tiny", "1", "50"}, {"aromatic", "51", "150"},
            {"druglike", "151", "350"}, {"complex", "351", "550"},
            {"large", "551", "700"}, {"very_large", "701", "800"},
            {"tautomer", "801", "900"}, {"symmetric", "901", "950"},
            {"edge", "951", "1050"}
        };

        sb.append("--- Per-Section Breakdown (median time ms) ---\n");
        sb.append(String.format("%-15s %12s %12s %8s%n", "Section", "Median(ms)", "Mean(ms)", "Pairs"));
        sb.append("-".repeat(50)).append("\n");

        for (String[] sec : sections) {
            String secName = sec[0];
            int lo = Integer.parseInt(sec[1]);
            int hi = Integer.parseInt(sec[2]);

            List<Double> secTimes = new ArrayList<>();
            int secPairCount = 0;
            for (Pair p : pairs) {
                boolean inSection = (p.idxA() + 1 >= lo && p.idxA() + 1 <= hi)
                                 || (p.idxB() + 1 >= lo && p.idxB() + 1 <= hi);
                if (inSection) {
                    String key = p.nameA() + "__vs__" + p.nameB();
                    BenchResult br = results.get(key);
                    if (br != null && br.medianTimeMs() >= 0) {
                        secTimes.add(br.medianTimeMs());
                    }
                    secPairCount++;
                }
            }

            if (!secTimes.isEmpty()) {
                Collections.sort(secTimes);
                double secMedian = secTimes.get(secTimes.size() / 2);
                double secMean = secTimes.stream().mapToDouble(Double::doubleValue).average().orElse(-1);
                sb.append(String.format("%-15s %12.3f %12.3f %8d%n", secName, secMedian, secMean, secPairCount));
            } else {
                sb.append(String.format("%-15s %12s %12s %8d%n", secName, "N/A", "N/A", secPairCount));
            }
        }

        // Slowest 10 pairs
        sb.append("\n--- Slowest 10 Pairs ---\n");
        List<Map.Entry<String, BenchResult>> sorted = results.entrySet().stream()
            .sorted((a, b) -> Double.compare(b.getValue().medianTimeMs(), a.getValue().medianTimeMs()))
            .limit(10)
            .collect(Collectors.toList());
        for (Map.Entry<String, BenchResult> e : sorted) {
            sb.append(String.format("  %-50s  %.3f ms  MCS=%d%n",
                e.getKey(), e.getValue().medianTimeMs(), e.getValue().mcsSize()));
        }

        String text = sb.toString();
        Files.writeString(outPath, text);
        System.out.println(text);
    }

    // ======================================================================
    // Main
    // ======================================================================
    public static void main(String[] args) throws Exception {
        Path molPath = args.length > 0 ? Path.of(args[0]) : Path.of("benchmarks/diverse_molecules.txt");
        Path outTsv = molPath.resolveSibling("benchmark_smsd_results.tsv");
        Path outSummary = molPath.resolveSibling("benchmark_smsd_summary.txt");

        System.err.printf("Loading molecules from %s ...%n", molPath);
        List<Molecule> mols = loadMolecules(molPath);
        System.err.printf("  %d valid molecules loaded%n", mols.size());

        if (mols.size() < 10) {
            System.err.println("ERROR: Too few valid molecules. Aborting.");
            System.exit(1);
        }

        System.err.printf("Generating %d random + %d systematic pairs ...%n",
            N_RANDOM_PAIRS, N_SYSTEMATIC_PAIRS);
        List<Pair> pairs = generatePairs(mols);
        System.err.printf("  Generated %d pairs%n", pairs.size());

        // JVM warmup — parse a few molecules to get JIT going
        System.err.println("JVM warmup ...");
        if (mols.size() >= 2) {
            ChemOptions warmOpts = new ChemOptions();
            SearchEngine.McsOptions warmMcs = new SearchEngine.McsOptions();
            warmMcs.timeoutMillis = 2000;
            for (int w = 0; w < 10; w++) {
                try {
                    SearchEngine.findMCS(mols.get(0).mol(), mols.get(1).mol(), warmOpts, warmMcs);
                } catch (Exception e) { /* ignore */ }
            }
        }

        // Run benchmark
        System.err.printf("Running SMSD MCS benchmark (%d rounds, %d ms timeout) ...%n",
            ROUNDS, TIMEOUT_MS);
        Map<String, BenchResult> results = new LinkedHashMap<>();
        int total = pairs.size();
        for (int i = 0; i < total; i++) {
            Pair p = pairs.get(i);
            String key = p.nameA() + "__vs__" + p.nameB();
            Molecule molA = mols.get(p.idxA());
            Molecule molB = mols.get(p.idxB());

            BenchResult br = benchmarkPair(molA, molB);
            results.put(key, br);

            if ((i + 1) % 50 == 0 || i == total - 1) {
                System.err.printf("  %d/%d pairs done (latest: %.3f ms, MCS=%d)%n",
                    i + 1, total, br.medianTimeMs(), br.mcsSize());
            }
        }

        // Output
        System.err.printf("Writing results to %s ...%n", outTsv);
        writeResults(pairs, results, outTsv);

        System.err.printf("Writing summary to %s ...%n", outSummary);
        writeSummary(pairs, results, mols.size(), outSummary);

        System.err.println("Done.");
    }
}

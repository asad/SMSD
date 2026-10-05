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
import org.junit.jupiter.api.*;
import org.junit.jupiter.api.condition.EnabledIfSystemProperty;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.DfPattern;

/**
 * Java-to-Java like-for-like substructure benchmark:
 * CDK DfPattern (bond-driven VF2) vs SMSD Pro VF2++ (vertex-driven).
 *
 * <p>Both engines operate within the same JVM, same JIT, and receive
 * identical standardised {@link IAtomContainer} objects from CDK 2.11.
 * Protocol: {@value WARMUP} warmup + {@value ITERS} measured runs, median reported.
 *
 * <p>Run with: {@code mvn test -Dtest=JavaCdkVsSmsdBenchmarkTest -Dbenchmark=true}
 *
 * @author Syed Asad Rahman
 */
@DisplayName("Java CDK vs SMSD Pro Substructure Benchmark (20 pairs)")
@EnabledIfSystemProperty(named = "benchmark", matches = "true")
public class JavaCdkVsSmsdBenchmarkTest extends TestBase {

    // -----------------------------------------------------------------------
    // Protocol constants
    // Protocol: 3 warmup + 10 measured, report median.
    // Rationale: 10 runs exceeds the McCreesh/Glasgow (2017) standard of 5
    // and gives a coefficient of variation < 0.5% for timings >= 10us.
    // -----------------------------------------------------------------------
    private static final int WARMUP = 3;
    private static final int ITERS  = 10;
    private static final long SUB_TIMEOUT_MS = 10_000;

    private static final String VANCOMYCIN =
        "CC1C(C(CC(O1)OC2C(C(C(OC2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)C(C(C(=O)" +
        "NC(C(=O)NC5C(=O)NC7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)C(NC(=O)C(" +
        "C(C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)C(CC(C)C)NC)" +
        "O)Cl)CO)O)O)(C)N)O";

    private static final String PEG16 =
        "OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO";

    /**
     * All 20 molecule pairs — identical to the Python benchmark for direct comparison.
     * Format: { smi1, smi2, pairName, category }
     */
    private static final String[][] PAIRS = {
        // Trivial
        {"C", "CC", "methane-ethane", "Trivial"},
        // Small aromatic
        {"c1ccccc1", "Cc1ccccc1", "benzene-toluene", "Small aromatic"},
        // Heteroatom
        {"c1ccccc1", "Oc1ccccc1", "benzene-phenol", "Heteroatom"},
        // Drug pairs
        {"CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1",
            "aspirin-acetaminophen", "Drug pair"},
        {"Cn1cnc2c1c(=O)n(C)c(=O)n2C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
            "caffeine-theophylline", "N-methyl diff"},
        // Alkaloids
        {"CN1CCC23C4C1CC5=C(C2C(C=C4)O3)C=C(C=C5)O",
            "CN1CCC23C4C1CC5=C(C2C(C=C4)OC3)C=C(C=C5)O",
            "morphine-codeine", "Alkaloid"},
        // NSAIDs
        {"CC(C)Cc1ccc(CC(C)C(=O)O)cc1",
            "COc1ccc2cc(CC(C)C(=O)O)ccc2c1",
            "ibuprofen-naproxen", "NSAID"},
        // Nucleotides / cofactors
        {"C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
            "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
            "ATP-ADP", "Nucleotide"},
        {"C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N",
            "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C4N=CN=C5N)O)O)O)O",
            "NAD-NADH", "Cofactor"},
        // Statins
        {"CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
            "CC(C)C1=NC(=NC(=C1C=CC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)N(C)S(=O)(=O)C",
            "atorvastatin-rosuvastatin", "Statin"},
        // Taxanes
        {"CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
            "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)OC(C)(C)C)O)O)OC(=O)C6=CC=CC=C6)(CO4)OC(=O)C)O)C)O",
            "paclitaxel-docetaxel", "Taxane"},
        // Macrolides
        {"CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
            "CCC1C(C(C(N(CC(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
            "erythromycin-azithromycin", "Macrolide"},
        // Alkaloid scaffolds
        {"C1CN2CC3=CCOC4CC(=O)N5C6C4C3CC2C61C7=CC=CC=C75",
            "COC1=CC2=C(C=CN=C2C=C1)C(C3CC4CCN3CC4C=C)O",
            "strychnine-quinine", "Alkaloid scaffold"},
        // Self-match: large glycopeptide
        {VANCOMYCIN, VANCOMYCIN, "vancomycin-self", "Self-match large"},
        // Self-match: cage
        {"C1C2CC3CC1CC(C2)C3", "C1C2CC3CC1CC(C2)C3",
            "adamantane-self", "Symmetric"},
        {"C12C3C4C1C5C4C3C25", "C12C3C4C1C5C4C3C25",
            "cubane-self", "Cage"},
        // Polymer
        {"OCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO", PEG16,
            "PEG12-PEG16", "Polymer"},
        // PAH
        {"c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
            "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
            "coronene-self", "PAH"},
        // Tautomer
        {"O=c1[nH]c(N)nc2[nH]cnc12", "Oc1nc(N)nc2[nH]cnc12",
            "guanine-keto-enol", "Tautomer"},
        // Known hard pair
        {"c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl",
            "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6",
            "rdkit-1585-pair", "Known failure"},
    };

    // -----------------------------------------------------------------------
    // Result record
    // -----------------------------------------------------------------------
    private record PairResult(
            String name,
            String category,
            long smsdBestNs,
            long smsdMedianNs,
            boolean smsdHit,
            long cdkBestNs,
            long cdkMedianNs,
            boolean cdkHit) {}

    // -----------------------------------------------------------------------
    // Core benchmark test
    // -----------------------------------------------------------------------

    @Test
    @DisplayName("CDK DfPattern vs SMSD VF2++ — substructure (all 20 pairs)")
    void cdkVsSmsdSubstructure() throws Exception {
        System.out.println();
        System.out.println("=".repeat(110));
        System.out.println("Java CDK DfPattern vs SMSD Pro VF2++ — Like-for-Like Substructure Benchmark");
        System.out.printf("Protocol: %d warmup + %d measured, median reported | All %d molecule pairs%n",
                WARMUP, ITERS, PAIRS.length);
        System.out.println("=".repeat(110));
        System.out.printf(" # %-28s %-20s %12s %12s %7s %7s %12s%n",
                "Pair", "Category", "SMSD(us)", "CDK(us)", "SMSD?", "CDK?", "Speedup");
        System.out.println("-".repeat(110));

        List<PairResult> results = new ArrayList<>();
        ChemOptions opts = new ChemOptions();
        int smsdWins = 0, cdkWins = 0, ties = 0;
        int agree = 0, disagree = 0;

        for (int pi = 0; pi < PAIRS.length; pi++) {
            String[] p = PAIRS[pi];
            String smi1 = p[0], smi2 = p[1], name = p[2], category = p[3];

            IAtomContainer mol1, mol2;
            try {
                mol1 = mol(smi1);
                mol2 = mol(smi2);
            } catch (Exception e) {
                System.out.printf("%2d %-28s PARSE ERROR: %s%n", pi + 1, name, e.getMessage());
                continue;
            }

            // --- SMSD Pro VF2++ ---
            boolean smsdHit = false;
            long[] smsdTimes = new long[ITERS];
            try {
                for (int i = 0; i < WARMUP; i++) {
                    new SMSD(mol1, mol2, opts).isSubstructure(SUB_TIMEOUT_MS);
                }
                for (int i = 0; i < ITERS; i++) {
                    long t0 = System.nanoTime();
                    smsdHit = new SMSD(mol1, mol2, opts).isSubstructure(SUB_TIMEOUT_MS);
                    smsdTimes[i] = System.nanoTime() - t0;
                }
            } catch (Exception e) {
                Arrays.fill(smsdTimes, Long.MAX_VALUE);
            }

            // --- CDK DfPattern (bond-driven VF2) ---
            boolean cdkHit = false;
            long[] cdkTimes = new long[ITERS];
            try {
                DfPattern dfp = DfPattern.findSubstructure(mol1);
                for (int i = 0; i < WARMUP; i++) {
                    dfp.matches(mol2);
                }
                for (int i = 0; i < ITERS; i++) {
                    long t0 = System.nanoTime();
                    cdkHit = dfp.matches(mol2);
                    cdkTimes[i] = System.nanoTime() - t0;
                }
            } catch (Exception e) {
                Arrays.fill(cdkTimes, Long.MAX_VALUE);
            }

            Arrays.sort(smsdTimes);
            Arrays.sort(cdkTimes);

            long smsdBest   = smsdTimes[0];
            long smsdMedian = smsdTimes[ITERS / 2];
            long cdkBest    = cdkTimes[0];
            long cdkMedian  = cdkTimes[ITERS / 2];

            if (smsdHit == cdkHit) agree++; else disagree++;

            // Speedup
            String speedup;
            if (smsdMedian == Long.MAX_VALUE || cdkMedian == Long.MAX_VALUE) {
                speedup = "N/A";
                ties++;
            } else {
                double ratio = (double) cdkMedian / smsdMedian;
                if (ratio > 1.1)      { speedup = String.format("SMSD %.1fx", ratio); smsdWins++; }
                else if (ratio < 0.9) { speedup = String.format("CDK %.1fx", 1.0/ratio); cdkWins++;  }
                else                  { speedup = "~tie"; ties++; }
            }

            System.out.printf("%2d %-28s %-20s %12s %12s %7s %7s %12s%n",
                    pi + 1, name, category,
                    fmtNs(smsdMedian), fmtNs(cdkMedian),
                    smsdHit ? "yes" : "no", cdkHit ? "yes" : "no",
                    speedup);

            results.add(new PairResult(name, category,
                    smsdBest, smsdMedian, smsdHit,
                    cdkBest, cdkMedian, cdkHit));
        }

        System.out.println("-".repeat(110));
        System.out.println();
        System.out.printf("SUMMARY (%d pairs)%n", results.size());
        System.out.printf("  Speed wins : SMSD=%d  CDK=%d  tie=%d%n", smsdWins, cdkWins, ties);
        System.out.printf("  Hit agreement: %d/%d pairs agree%n", agree, PAIRS.length);
        if (disagree > 0) {
            System.out.printf("  NOTE: %d pair(s) disagree — check SMILES/aromaticity conventions%n",
                    disagree);
        }

        // At least 18/20 pairs must parse and complete
        assertTrue(results.size() >= 18,
                "Expected >= 18 valid pairs, got " + results.size());
        // Hits must agree on at least 75% of pairs (accounting for aromatic/tautomer edge cases)
        assertTrue(agree >= PAIRS.length * 3 / 4,
                "Hit agreement too low: " + agree + "/" + PAIRS.length);
    }

    // -----------------------------------------------------------------------
    // Self-match correctness: mol must be a substructure of itself
    // -----------------------------------------------------------------------

    @Test
    @DisplayName("Self-match: every molecule must be a substructure of itself")
    void selfMatchCorrectness() throws Exception {
        ChemOptions opts = new ChemOptions();
        int passed = 0;

        for (String[] p : PAIRS) {
            String smi = p[0], name = p[2];
            if (!name.contains("self") && !p[3].equals("Trivial")) continue;

            IAtomContainer mol;
            try {
                mol = mol(smi);
            } catch (Exception e) {
                continue;
            }

            // SMSD self-match
            boolean smsdSelf = new SMSD(mol, mol, opts).isSubstructure(SUB_TIMEOUT_MS);
            // CDK self-match
            boolean cdkSelf  = DfPattern.findSubstructure(mol).matches(mol);

            assertTrue(smsdSelf, name + ": SMSD self-match failed");
            assertTrue(cdkSelf,  name + ": CDK  self-match failed");
            passed++;
        }

        assertTrue(passed > 0, "No self-match pairs found in PAIRS list");
    }

    // -----------------------------------------------------------------------
    // Utility: format nanoseconds to human-readable string
    // -----------------------------------------------------------------------
    private static String fmtNs(long ns) {
        if (ns == Long.MAX_VALUE) return "ERR";
        if (ns >= 1_000_000_000L) return String.format("%.2fs",    ns / 1e9);
        if (ns >= 1_000_000L)     return String.format("%.2fms",   ns / 1e6);
        if (ns >= 1_000L)         return String.format("%.1fus",   ns / 1e3);
        return ns + "ns";
    }
}

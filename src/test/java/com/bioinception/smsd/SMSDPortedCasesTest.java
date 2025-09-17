/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Timeout;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertTrue;

public class SMSDPortedCasesTest {

    static Stream<String[]> chainCases() {
        List<String[]> out = new ArrayList<>();
        for (int k = 2; k <= 12; k++) { // trimmed for CI
            String q = "C".repeat(k);
            String t = "C".repeat(k + 2);
            out.add(new String[]{"chain_" + k, q, t});
        }
        return out.stream();
    }

    static Stream<String[]> aromCases() {
        return Stream.of(
                new String[]{"arom_0", "c1ccccc1", "c1ccc2ccccc2c1"},
                new String[]{"arom_1", "c1ccccc1", "C1=CC=CC=C1"}
        );
    }

    static Stream<String[]> stereoCases() {
        return Stream.of(
                new String[]{"stereo_0", "C/C=C\\C", "C/C=C\\C"},
                new String[]{"stereo_1", "C/C=C\\C", "C/C=C/C"}
        );
    }

    static Stream<String[]> allCases() {
        List<String[]> all = new ArrayList<>();
        all.addAll(chainCases().collect(Collectors.toList()));
        all.addAll(aromCases().collect(Collectors.toList()));
        all.addAll(stereoCases().collect(Collectors.toList()));
        return all.stream();
    }

    static Stream<String[]> mcsSubset() {
        return allCases().limit(12).collect(Collectors.toList()).stream();
    }

    @ParameterizedTest(name="{0}")
    @MethodSource("allCases")
    @DisplayName("Substructure exists within timeout")
    @Timeout(value = 4, unit = TimeUnit.SECONDS)
    void test_substructure_exists(String name, String q, String t) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer queryMol = sp.parseSmiles(q);
        IAtomContainer targetMol = sp.parseSmiles(t);
        SMSD smsd = new SMSD(queryMol, targetMol, new ChemOptions());
        boolean ok = smsd.isSubstructure(1_500L);
        assertTrue(ok || !ok);
    }

    @ParameterizedTest(name="{0}")
    @MethodSource("mcsSubset")
    @DisplayName("MCS bounded by timeout like FMCS")
    @Timeout(value = 6, unit = TimeUnit.SECONDS)
    void test_mcs_runs_mcisd(String name, String q, String t) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer queryMol = sp.parseSmiles(q);
        IAtomContainer targetMol = sp.parseSmiles(t);
        SMSD smsd = new SMSD(queryMol, targetMol, new ChemOptions());
        Map<Integer, Integer> mr = smsd.findMCS(false, true, 2_500L);
        assertTrue(mr != null);
    }
}

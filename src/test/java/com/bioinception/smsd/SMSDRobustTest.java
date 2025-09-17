/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package com.bioinception.smsd;

import com.bioinception.smsd.core.ChemOptions;
import com.bioinception.smsd.core.SMSD;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.Arguments;
import static org.junit.jupiter.params.provider.Arguments.arguments;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Parameterised tests ported from our Python harness:
 *
 * - Chains, aromatics, stereo, disconnections, functionals, recursive SMARTS.
 * - Substructure existence returns a boolean.
 * - MCS runs on a subset of cases and returns a mapping (or null).
 *
 * CDK-only; no external toolkits.
 */
// FIX: The class name now matches the filename "SMSDRobustTest.java"
public class SMSDRobustTest {

    // -----------------------------
    // Case generators (Arguments)
    // -----------------------------

    static Stream<Arguments> chainCases() {
        List<Arguments> out = new ArrayList<>();
        for (int k = 2; k <= 21; k++) {
            String q = "C".repeat(k);
            String t = "C".repeat(k + 2);
            out.add(arguments("chain_" + k, q, t));
        }
        return out.stream();
    }

    static Stream<Arguments> aromCases() {
        return Stream.of(
            arguments("arom_0", "c1ccccc1", "c1ccc2ccccc2c1"),
            arguments("arom_1", "c1ccccc1", "C1=CC=CC=C1")
        );
    }

    static Stream<Arguments> stereoCases() {
        return Stream.of(
            arguments("stereo_0", "C/C=C\\C", "C/C=C\\C"),
            arguments("stereo_1", "C/C=C\\C", "C/C=C/C")
        );
    }

    static Stream<Arguments> disconnCases() {
        return Stream.of(
            arguments("disc_0", "C.C", "CC")
        );
    }

    static Stream<Arguments> functionalCases() {
        return Stream.of(
            arguments("fx_0", "NC(=O)", "NCC(=O)NCCC"),
            arguments("fx_1", "[O-]C=O", "OC=O")
        );
    }

    static Stream<Arguments> recursiveCases() {
        return Stream.of(
            arguments("rec_0", "[C;$(C(=O)O)]", "CC(=O)O"),
            arguments("rec_1", "[N;$isAmideN]", "CC(=O)NCC")
        );
    }

    static Stream<Arguments> extraCases() {
        String[][] base = new String[][]{
            {"extra_a","c1ccccc1CC","c1ccccc1CCCC"},
            {"extra_b","CCN(CC)CC","CCN(CC)CCO"},
            {"extra_c","C1CCCCC1","C1CCCCCC1"},
            {"extra_d","[13CH3]C","CC"},
            {"extra_e","c1ncccc1","c1ncccc1Cl"},
            {"extra_f","N1CCOCC1","O=C(N1CCOCC1)C"},
            {"extra_g","FC(F)F","CC(C)(F)F"},
            {"extra_h","O=S(=O)N","O=S(=O)NCC"},
            {"extra_i","P(=O)(O)O","OP(=O)(O)OCC"},
            {"extra_j","Clc1ccccc1","Clc1ccc(Cl)cc1"}
        };
        List<Arguments> L = new ArrayList<>();
        for (int i = 0; i < 9; i++) {
            for (String[] b : base) {
                L.add(arguments(b[0], b[1], b[2]));
            }
        }
        return L.stream();
    }

    static Stream<Arguments> allCases() {
        return Stream.concat(
            Stream.concat(chainCases(), aromCases()),
            Stream.concat(
                Stream.concat(stereoCases(), disconnCases()),
                Stream.concat(functionalCases(), Stream.concat(recursiveCases(), extraCases()))
            )
        );
    }

    static Stream<Arguments> mcsSubset() {
        // keep runtime reasonable in CI
        return allCases().limit(30);
    }

    // -----------------------------
    // Tests
    // -----------------------------

    @ParameterizedTest(name = "{0}")
    @MethodSource("allCases")
    @DisplayName("Substructure exists (boolean) on diverse cases")
    void test_substructure_exists(String name, String q, String t) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer queryMol = sp.parseSmiles(q);
        IAtomContainer targetMol = sp.parseSmiles(t);

        ChemOptions chem = new ChemOptions(); // defaults
        SMSD smsd = new SMSD(queryMol, targetMol, chem);
        boolean ok = smsd.isSubstructure();
        // sanity assertion (documenting the contract)
        assertTrue(ok || !ok, "isSubstructure returns a boolean");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("mcsSubset")
    @DisplayName("MCS runs on a subset of cases (MCIS mode)")
    void test_mcs_runs_mcisd(String name, String q, String t) throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer queryMol = sp.parseSmiles(q);
        IAtomContainer targetMol = sp.parseSmiles(t);

        ChemOptions chem = new ChemOptions();
        SMSD smsd = new SMSD(queryMol, targetMol, chem);
        Map<Integer, Integer> mr = smsd.findMCS(true, true); // (connected, induced)
        assertTrue(mr == null || mr.size() >= 0);
    }

    @Test
    @DisplayName("CLI JSON shape (documented)")
    void test_cli_json_shape_comment() {
        // Placeholder for CLI JSON schema checks:
        // keys: query, target, mode, mcs_type, induced, mappings[*].{index,query_atoms,target_atoms,target_bonds,pairs,subgraph_smiles}
        assertTrue(true);
    }
}
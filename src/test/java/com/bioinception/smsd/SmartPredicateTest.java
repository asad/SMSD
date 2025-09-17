/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */

package com.bioinception.smsd;

import com.bioinception.smsd.core.*;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;

public class SmartPredicateTest {

    @Test
    public void testNamedPredicateExpansionAndMatch() throws Exception {
        String sig = "[C;$isKetone]";
        PredicateRegistry reg = new PredicateRegistry();
        String expanded = SmartPredicateProcessor.expandNamedPredicates(sig, reg);
        assertTrue(expanded.contains("$("));
        IAtomContainer t = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CC(=O)C");
        t = Standardiser.standardise(t, Standardiser.TautomerMode.NONE);
        java.util.List<int[]> hits = SmartsUtil.matchAll(expanded, t);
        assertFalse(hits.isEmpty());
    }
}

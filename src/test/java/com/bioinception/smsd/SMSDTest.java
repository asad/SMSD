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

import java.util.Map;

public class SMSDTest {

    @Test
    public void testSubstructureExists() throws Exception {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer q = sp.parseSmiles("CCO");
        IAtomContainer t = sp.parseSmiles("CCNCCO");
        SMSD smsd = new SMSD(Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
                             Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
                             ChemOptions.profile("compat-substruct"));
        assertTrue(smsd.isSubstructure());
    }

    @Test
    public void testSubstructureNotExists() throws Exception {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer q = sp.parseSmiles("c1ccccc1");
        IAtomContainer t = sp.parseSmiles("CCCCCC");
        SMSD smsd = new SMSD(Standardiser.standardise(q, Standardiser.TautomerMode.NONE),
                             Standardiser.standardise(t, Standardiser.TautomerMode.NONE),
                             ChemOptions.profile("strict"));
        assertFalse(smsd.isSubstructure());
    }

    @Test
    public void testMCSBenzeneNaphthalene() throws Exception {
        SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer a = sp.parseSmiles("c1ccccc1");
        IAtomContainer b = sp.parseSmiles("c1ccc2ccccc2c1");
        SMSD smsd = new SMSD(Standardiser.standardise(a, Standardiser.TautomerMode.NONE),
                             Standardiser.standardise(b, Standardiser.TautomerMode.NONE),
                             ChemOptions.profile("compat-fmcs"));
        Map<Integer,Integer> m = smsd.findMCS(true, true);
        assertNotNull(m);
        assertEquals(6, m.size());
    }
}

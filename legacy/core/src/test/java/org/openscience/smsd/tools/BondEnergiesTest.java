/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.tools;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond.Order;

/**
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * test-smsd
 */
public class BondEnergiesTest {

    @Test
    public void testGetInstance() throws Exception {
        BondEnergies energies = BondEnergies.getInstance();
        Assert.assertNotNull(energies);
    }

    /**
     * Test of getEnergies method, of class BondEnergies.
     */
    @Test
    public void testGetEnergies() {
        //////System.out.println("getEnergies");
        IAtom sourceAtom = new Atom("C");
        IAtom targetAtom = new Atom("C");
        Order bondOrder = Order.SINGLE;
        BondEnergies instance = new BondEnergies();
        Integer expResult = 346;
        Integer result = instance.getEnergies(sourceAtom, targetAtom, bondOrder);
        Assert.assertEquals(expResult, result);
    }
}

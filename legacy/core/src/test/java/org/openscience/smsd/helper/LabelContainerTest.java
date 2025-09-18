/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.helper;

import org.junit.Assert;
import org.junit.Test;

/**
 * Unit testing for the {@link LabelContainer} class.
 *
 * Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 * test-smsd
 */
public class LabelContainerTest {

    @Test
    public void testGetInstance() {
        Assert.assertNotNull(LabelContainer.getInstance());
    }

    /**
     * Test of addLabel method, of class LabelContainer.
     */
    @Test
    public void testAddLabel() {
        //////System.out.println("addLabel");
        String label = "R3";
        LabelContainer instance = new LabelContainer();
        instance.addLabel(label);
        Assert.assertEquals(3, instance.getSize());
        Integer expectedValue = 2;
        Assert.assertEquals(expectedValue, instance.getLabelID("R3"));
    }

    /**
     * Test of getLabelID method, of class LabelContainer.
     */
    @Test
    public void testGetLabelID() {
        //////System.out.println("getLabelID");
        String label = "R3";
        LabelContainer instance = new LabelContainer();
        instance.addLabel(label);
        Integer expectedValue = 2;
        Assert.assertEquals(expectedValue, instance.getLabelID("R3"));
    }

    /**
     * Test of getLabel method, of class LabelContainer.
     */
    @Test
    public void testGetLabel() {
        //////System.out.println("getLabel");
        String label = "R3";
        LabelContainer instance = new LabelContainer();
        instance.addLabel(label);
        Integer index = 2;
        String result = instance.getLabel(index);
        Assert.assertEquals(label, result);
    }

    /**
     * Test of getSize method, of class LabelContainer.
     */
    @Test
    public void testGetSize() {
        //////System.out.println("getSize");
        String label = "R3";
        LabelContainer instance = new LabelContainer();
        instance.addLabel(label);
        int expectedValue = 3;
        int result = instance.getSize();
        Assert.assertEquals(expectedValue, result);
    }
}

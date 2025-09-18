/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.smsd.algorithm.rgraph;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Assert;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

/**
 *  test-smsd
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *  java1.8+
 */
public class CDKRMapHandlerTest {

    public CDKRMapHandlerTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getSource method, of class CDKRMapHandler.
     */
    @Test
    public void testGetSource() {
        ////////System.out.println("getSource");
        IAtomContainer expResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        CDKRMapHandler rMapHandler = new CDKRMapHandler();
        rMapHandler.setSource(expResult);
        IAtomContainer result = rMapHandler.getSource();
        Assert.assertEquals(expResult, result);
    }

    /**
     * Test of setSource method, of class CDKRMapHandler.
     */
    @Test
    public void testSetSource() {
        ////////System.out.println("setSource");
        IAtomContainer expResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        CDKRMapHandler rMapHandler = new CDKRMapHandler();
        rMapHandler.setSource(expResult);
        IAtomContainer result = rMapHandler.getSource();
        Assert.assertEquals(expResult, result);
    }

    /**
     * Test of getTarget method, of class CDKRMapHandler.
     */
    @Test
    public void testGetTarget() {
        ////////System.out.println("getTarget");
        IAtomContainer expResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        CDKRMapHandler rMapHandler = new CDKRMapHandler();
        rMapHandler.setTarget(expResult);
        IAtomContainer result = rMapHandler.getTarget();
        Assert.assertEquals(expResult, result);
    }

    /**
     * Test of setTarget method, of class CDKRMapHandler.
     */
    @Test
    public void testSetTarget() {
        ////////System.out.println("setTarget");
        IAtomContainer expResult = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        CDKRMapHandler rMapHandler = new CDKRMapHandler();
        rMapHandler.setTarget(expResult);
        IAtomContainer result = rMapHandler.getTarget();
        Assert.assertEquals(expResult, result);
    }

    /**
     * Test of calculateOverlapsAndReduce method, of class CDKRMapHandler.
     * @throws java.lang.Exception
     */
    @Test
    public void testCalculateOverlapsAndReduce() throws Exception {
        ////////System.out.println("calculateOverlapsAndReduce");
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        SmilesParser sp = new SmilesParser(builder);
        IAtomContainer Molecule1 = sp.parseSmiles("O1C=CC=C1");
        IAtomContainer Molecule2 = sp.parseSmiles("C1CCCC1");
        CDKRMapHandler instance = new CDKRMapHandler();
        instance.calculateOverlapsAndReduce(Molecule1, Molecule2, true, false,false);
        Assert.assertNotNull(instance.getMappings().size());
    }

    /**
     * Test of calculateOverlapsAndReduceExactMatch method, of class CDKRMapHandler.
     * @throws java.lang.Exception
     */
    @Test
    public void testCalculateOverlapsAndReduceExactMatch() throws Exception {
        ////////System.out.println("calculateOverlapsAndReduceExactMatch");
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        SmilesParser sp = new SmilesParser(builder);
        IAtomContainer Molecule1 = sp.parseSmiles("O1C=CC=C1");
        IAtomContainer Molecule2 = sp.parseSmiles("O1C=CC=C1");
        CDKRMapHandler instance = new CDKRMapHandler();
        instance.calculateOverlapsAndReduceExactMatch(Molecule1, Molecule2, true, false,false);
        // TODO review the generated test code and remove the default call to fail.
        Assert.assertNotNull(instance.getMappings());
    }

    /**
     * Test of getMappings method, of class CDKRMapHandler.
     * @throws org.openscience.cdk.exception.InvalidSmilesException
     */
    @Test
    public void testGetMappings() throws InvalidSmilesException, CDKException {
        ////////System.out.println("getMappings");
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        SmilesParser sp = new SmilesParser(builder);
        IAtomContainer Molecule1 = sp.parseSmiles("O1C=CC=C1");
        IAtomContainer Molecule2 = sp.parseSmiles("O1C=CC=C1");
        CDKRMapHandler instance = new CDKRMapHandler();
        instance.calculateOverlapsAndReduceExactMatch(Molecule1, Molecule2, true, false,false);
        List<Map<Integer, Integer>> result = instance.getMappings();
        Assert.assertEquals(2, result.size());
    }

    /**
     * Test of setMappings method, of class CDKRMapHandler.
     */
    @Test
    public void testSetMappings() {
        ////////System.out.println("setMappings");
        Map<Integer, Integer> map = new TreeMap<>();
        map.put(0, 0);
        map.put(1, 1);

        List<Map<Integer, Integer>> mappings = new ArrayList<>();
        mappings.add(map);
        CDKRMapHandler instance = new CDKRMapHandler();
        instance.setMappings(mappings);
        Assert.assertNotNull(instance.getMappings());
    }

    /**
     * Test of isTimeoutFlag method, of class CDKRMapHandler.
     */
    @Test
    public void testIsTimeoutFlag() {
        ////////System.out.println("isTimeoutFlag");
        CDKRMapHandler instance = new CDKRMapHandler();
        boolean expResult = true;
        instance.setTimeout(true);
        boolean result = instance.isTimeout();
        Assert.assertEquals(expResult, result);
    }

    /**
     * Test of setTimeoutFlag method, of class CDKRMapHandler.
     */
    @Test
    public void testSetTimeoutFlag() {
        ////////System.out.println("setTimeoutFlag");
        boolean timeoutFlag = false;
        CDKRMapHandler instance = new CDKRMapHandler();
        instance.setTimeout(timeoutFlag);
        Assert.assertNotSame(true, instance.isTimeout());
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.vflib;

import org.openscience.cdk.DefaultChemObjectBuilder;
import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.interfaces.AbstractSubGraphTest;
import org.openscience.smsd.tools.TimeManager;

/**
 *
 * @author Asad
 */
public class VF2SubstructureTest extends AbstractSubGraphTest{
    private static final long serialVersionUID = 1L;
    
    public VF2SubstructureTest() {
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
     * Test of setTimeManager method, of class VF2Substructure.
     */
    @Test
    public void testSetTimeManager() {
        System.out.println("setTimeManager");
        TimeManager aTimeManager = null;
        VF2Substructure instance = new VF2Substructure();
        instance.setTimeManager(aTimeManager);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isSubgraph method, of class VF2Substructure.
     */
    @Test
    public void testIsSubgraph() {
        System.out.println("isSubgraph");
        boolean shouldMatchBonds = false;
        boolean shouldMatchRings = false;
        VF2Substructure instance = new VF2Substructure();
        boolean expResult = false;
        boolean result = instance.isSubgraph(shouldMatchBonds, shouldMatchRings);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Substructure.
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() {
        System.out.println("set");
        IAtomContainer source = null;
        IAtomContainer target = null;
        VF2Substructure instance = new VF2Substructure();
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Substructure.
     */
    @Test
    public void testSet_IQueryAtomContainer_IAtomContainer() {
        System.out.println("set");
        IQueryAtomContainer source = null;
        IAtomContainer target = null;
        VF2Substructure instance = new VF2Substructure();
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllMapping method, of class VF2Substructure.
     */
    @Test
    public void testGetAllMapping() {
        System.out.println("getAllMapping");
        VF2Substructure instance = new VF2Substructure();
        List expResult = null;
        List result = instance.getAllMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstMapping method, of class VF2Substructure.
     */
    @Test
    public void testGetFirstMapping() {
        System.out.println("getFirstMapping");
        VF2Substructure instance = new VF2Substructure();
        Map expResult = null;
        Map result = instance.getFirstMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllAtomMapping method, of class VF2Substructure.
     */
    @Test
    public void testGetAllAtomMapping() {
        System.out.println("getAllAtomMapping");
        VF2Substructure instance = new VF2Substructure();
        List expResult = null;
        List result = instance.getAllAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2Substructure.
     */
    @Test
    public void testGetFirstAtomMapping() {
        System.out.println("getFirstAtomMapping");
        VF2Substructure instance = new VF2Substructure();
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.getFirstAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isBondMatchFlag method, of class VF2Substructure.
     */
    @Test
    public void testIsBondMatchFlag() {
        System.out.println("isBondMatchFlag");
        VF2Substructure instance = new VF2Substructure();
        boolean expResult = false;
        boolean result = instance.isBondMatchFlag();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setBondMatchFlag method, of class VF2Substructure.
     */
    @Test
    public void testSetBondMatchFlag() {
        System.out.println("setBondMatchFlag");
        boolean shouldMatchBonds = false;
        VF2Substructure instance = new VF2Substructure();
        instance.setBondMatchFlag(shouldMatchBonds);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isTimeOut method, of class VF2Substructure.
     */
    @Test
    public void testIsTimeOut() {
        System.out.println("isTimeOut");
        VF2Substructure instance = new VF2Substructure();
        boolean expResult = false;
        boolean result = instance.isTimeOut();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isMatchRings method, of class VF2Substructure.
     */
    @Test
    public void testIsMatchRings() {
        System.out.println("isMatchRings");
        VF2Substructure instance = new VF2Substructure();
        boolean expResult = false;
        boolean result = instance.isMatchRings();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of setMatchRings method, of class VF2Substructure.
     */
    @Test
    public void testSetMatchRings() {
        System.out.println("setMatchRings");
        boolean shouldMatchRings = false;
        VF2Substructure instance = new VF2Substructure();
        instance.setMatchRings(shouldMatchRings);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
        /**
     * Test linker to linker match
     * @throws InvalidSmilesException
     */
    @Test
    public void testLinkersSystemMatch() throws InvalidSmilesException {
        System.out.println("testLinkersSystemMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(5, smsd1.getFirstMapping().size());
    }

    /**
     * Test linker to ring match
     * @throws InvalidSmilesException
     */
    @Test
    public void testLinkersVsRingsMatch() throws InvalidSmilesException {
        System.out.println("testLinkersVsRingsMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");


        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, true);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(4, smsd1.getFirstMapping().size());
    }
}

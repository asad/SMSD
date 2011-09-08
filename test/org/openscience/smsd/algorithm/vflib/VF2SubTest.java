/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.vflib;

import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.tools.TimeManager;

/**
 *
 * @author Asad
 */
public class VF2SubTest {
    
    public VF2SubTest() {
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
     * Test of setTimeManager method, of class VF2Sub.
     */
    @Test
    public void testSetTimeManager() {
        System.out.println("setTimeManager");
        TimeManager aTimeManager = null;
        VF2Sub instance = null;
        instance.setTimeManager(aTimeManager);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isSubgraph method, of class VF2Sub.
     */
    @Test
    public void testIsSubgraph() {
        System.out.println("isSubgraph");
        VF2Sub instance = null;
        boolean expResult = false;
        boolean result = instance.isSubgraph();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Sub.
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() {
        System.out.println("set");
        IAtomContainer source = null;
        IAtomContainer target = null;
        VF2Sub instance = null;
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of set method, of class VF2Sub.
     */
    @Test
    public void testSet_IQueryAtomContainer_IAtomContainer() {
        System.out.println("set");
        IQueryAtomContainer source = null;
        IAtomContainer target = null;
        VF2Sub instance = null;
        instance.set(source, target);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllMapping method, of class VF2Sub.
     */
    @Test
    public void testGetAllMapping() {
        System.out.println("getAllMapping");
        VF2Sub instance = null;
        List expResult = null;
        List result = instance.getAllMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstMapping method, of class VF2Sub.
     */
    @Test
    public void testGetFirstMapping() {
        System.out.println("getFirstMapping");
        VF2Sub instance = null;
        Map expResult = null;
        Map result = instance.getFirstMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllAtomMapping method, of class VF2Sub.
     */
    @Test
    public void testGetAllAtomMapping() {
        System.out.println("getAllAtomMapping");
        VF2Sub instance = null;
        List expResult = null;
        List result = instance.getAllAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2Sub.
     */
    @Test
    public void testGetFirstAtomMapping() {
        System.out.println("getFirstAtomMapping");
        VF2Sub instance = null;
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.getFirstAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isBondMatchFlag method, of class VF2Sub.
     */
    @Test
    public void testIsBondMatchFlag() {
        System.out.println("isBondMatchFlag");
        VF2Sub instance = null;
        boolean expResult = false;
        boolean result = instance.isBondMatchFlag();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isTimeOut method, of class VF2Sub.
     */
    @Test
    public void testIsTimeOut() {
        System.out.println("isTimeOut");
        VF2Sub instance = null;
        boolean expResult = false;
        boolean result = instance.isTimeOut();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isMatchRings method, of class VF2Sub.
     */
    @Test
    public void testIsMatchRings() {
        System.out.println("isMatchRings");
        VF2Sub instance = null;
        boolean expResult = false;
        boolean result = instance.isMatchRings();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.interfaces;

import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.smsd.AtomAtomMapping;

/**
 *
 * @author Asad
 */
public class IResultsTest {
    
    public IResultsTest() {
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
     * Test of getAllAtomMapping method, of class IResults.
     */
    @Test
    public void testGetAllAtomMapping() {
        System.out.println("getAllAtomMapping");
        IResults instance = new IResultsImpl();
        List expResult = null;
        List result = instance.getAllAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAllMapping method, of class IResults.
     */
    @Test
    public void testGetAllMapping() {
        System.out.println("getAllMapping");
        IResults instance = new IResultsImpl();
        List expResult = null;
        List result = instance.getAllMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstAtomMapping method, of class IResults.
     */
    @Test
    public void testGetFirstAtomMapping() {
        System.out.println("getFirstAtomMapping");
        IResults instance = new IResultsImpl();
        AtomAtomMapping expResult = null;
        AtomAtomMapping result = instance.getFirstAtomMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getFirstMapping method, of class IResults.
     */
    @Test
    public void testGetFirstMapping() {
        System.out.println("getFirstMapping");
        IResults instance = new IResultsImpl();
        Map expResult = null;
        Map result = instance.getFirstMapping();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    public class IResultsImpl implements IResults {

        public List<AtomAtomMapping> getAllAtomMapping() {
            return null;
        }

        public List<Map<Integer, Integer>> getAllMapping() {
            return null;
        }

        public AtomAtomMapping getFirstAtomMapping() {
            return null;
        }

        public Map<Integer, Integer> getFirstMapping() {
            return null;
        }
    }
}

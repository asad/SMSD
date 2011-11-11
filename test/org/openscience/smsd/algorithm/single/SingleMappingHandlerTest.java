/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.openscience.smsd.algorithm.single;

import java.io.IOException;
import java.io.InputStream;
import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *
 * @author Asad
 */
public class SingleMappingHandlerTest {

    public SingleMappingHandlerTest() {
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
     * Test of set method, of class SingleMappingHandler.
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() {
        System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class SingleMappingHandler.
     * @throws Exception 
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IMolecule source = new Molecule();
        source.addAtom(atomSource);
        IMolecule target = new Molecule();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class SingleMappingHandler.
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() {
        System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IMolecule source = new Molecule();
        source.addAtom(atomSource);
        IMolecule target = new Molecule();
        target.addAtom(atomTarget);

        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of searchMCS method, of class SingleMappingHandler.
     */
    @Test
    public void testSearchMCS() {
        System.out.println("searchMCS");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getAllMapping());
        Assert.assertEquals(1, instance.getAllMapping().size());
    }

    /**
     * Test of getAllMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetAllMapping() {
        System.out.println("getAllMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getAllMapping());
    }

    /**
     * Test of getFirstMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetFirstMapping() {
        System.out.println("getFirstMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getFirstMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetAllAtomMapping() {
        System.out.println("getAllAtomMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getAllAtomMapping());
    }

    /**
     * Test of getFirstAtomMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetFirstAtomMapping() {
        System.out.println("getFirstAtomMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }
}

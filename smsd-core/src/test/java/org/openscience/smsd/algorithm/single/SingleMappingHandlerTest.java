
/* Copyright (C) 2009-2014 Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.single;

import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 * @cdk.require java1.6+
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
        ////////System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class SingleMappingHandler.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        ////////System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class SingleMappingHandler.
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() {
        ////////System.out.println("set");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);

        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of searchMCS method, of class SingleMappingHandler.
     */
    @Test
    public void testSearchMCS() {
        ////////System.out.println("searchMCS");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getAllAtomMapping());
        Assert.assertEquals(1, instance.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testgetAllAtomMapping() {
        ////////System.out.println("getAllAtomMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getAllAtomMapping());
    }

    /**
     * Test of getFirstMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetFirstMapping() {
        ////////System.out.println("getFirstMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetAllAtomMapping() {
        ////////System.out.println("getAllAtomMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getAllAtomMapping());
    }

    /**
     * Test of getFirstAtomMapping method, of class SingleMappingHandler.
     */
    @Test
    public void testGetFirstAtomMapping() {
        ////////System.out.println("getFirstAtomMapping");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler(source, target, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }
}

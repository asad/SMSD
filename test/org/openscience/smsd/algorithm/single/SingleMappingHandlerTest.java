/* Copyright (C) 2010  Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work,
 * which includes - but is not limited to - adding the above copyright notice to
 * the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received rAtomCount copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.algorithm.single;

import java.io.IOException;
import java.io.InputStream;
import junit.framework.Assert;
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
import org.openscience.smsd.interfaces.AbstractMCSAlgorithmTest;
import org.openscience.smsd.interfaces.Algorithm;


/**
 * Unit testing for the {@link SingleMappingHandler} class.
 * 
 * @author     egonw
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 */
public class SingleMappingHandlerTest extends AbstractMCSAlgorithmTest {

    @BeforeClass
    public static void setMCSAlgorithm() {
        AbstractMCSAlgorithmTest.setMCSAlgorithm(new SingleMappingHandler());
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class SingleMappingHandler.
     * @throws CDKException
     * @throws IOException 
     */
    @Test
    public void testSet_String_String() throws CDKException, IOException {
        System.out.println("set");
        String molfile = "data/mdl/decalin.mol";
        String queryfile = "data/mdl/decalin.mol";
        Molecule query = new Molecule();
        Molecule target = new Molecule();

        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(molfile);
        MDLV2000Reader reader = new MDLV2000Reader(ins, Mode.STRICT);
        reader.read(query);
        ins = this.getClass().getClassLoader().getResourceAsStream(queryfile);
        reader = new MDLV2000Reader(ins, Mode.STRICT);
        reader.read(target);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.DEFAULT, true, false);
        smsd1.setChemFilters(true, true, true);
        double score = 1.0;
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.0001);
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

        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of searchMCS method, of class SingleMappingHandler.
     */
    @Test
    @Override
    public void testSearchMCS() {
        System.out.println("searchMCS");
        IAtom atomSource = new Atom("R");
        IAtom atomTarget = new Atom("R");
        IAtomContainer source = new AtomContainer();
        source.addAtom(atomSource);
        IAtomContainer target = new AtomContainer();
        target.addAtom(atomTarget);
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
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
        SingleMappingHandler instance = new SingleMappingHandler();
        instance.set(source, target);
        instance.searchMCS(true, false);
        Assert.assertNotNull(instance.getFirstAtomMapping());
    }
}

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
package org.openscience.smsd.algorithm.vflib;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.InputStream;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader.Mode;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.interfaces.AbstractMCSAlgorithmTest;


/**
 * Unit testing for the {@link VF2lib} class.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd
 * @cdk.require java1.6+
 */
public class VFlibMCSHandlerTest extends AbstractMCSAlgorithmTest {

    public VFlibMCSHandlerTest() {
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

    @BeforeClass
    public static void setMCSAlgorithm() {
        AbstractMCSAlgorithmTest.setMCSAlgorithm(
                new VF2MCSPlusHandler());
    }

    /**
     * Test of searchMCS method, of class VF2lib.
     */
    @Test
    @Override
    public void testSearchMCS() {
        System.out.println("searchMCS");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = null;
        try {
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(VFlibMCSHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        IAtomContainer query = null;
        try {
            query = sp.parseSmiles("Nc1ccccc1");
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(VFlibMCSHandlerTest.class.getName()).log(Level.SEVERE, null, ex);
        }

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws InvalidSmilesException {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());

    }

    /**
     * Test of set method, of class VF2lib.
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IMolecule query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2lib.
     * @throws CDKException
     */
    @Test
    public void testSet_String_String() throws CDKException {
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

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);

        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws InvalidSmilesException {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        VF2lib instance = new VF2lib();
        instance.set(query, target);
        instance.searchMCS(true);
        assertNotNull(instance.getFirstMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllMapping method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllMapping() throws InvalidSmilesException {
        System.out.println("getAllMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(4, smsd1.getAllMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());

        Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstMapping method, of class VF2lib.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstMapping() throws InvalidSmilesException {
        System.out.println("getFirstMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2lib smsd1 = new VF2lib();
        smsd1.set(query, target);
        smsd1.searchMCS(true);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(7, smsd1.getFirstMapping().size());
    }
}

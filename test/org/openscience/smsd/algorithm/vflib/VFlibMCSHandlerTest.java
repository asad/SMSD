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
 * Unit testing for the {@link VF2MCSHandler} class.
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
                new VF2MCSHandler());
    }

    /**
     * Test of searchMCS method, of class VF2SubStructure.
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

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws InvalidSmilesException {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

    }

    /**
     * Test of set method, of class VF2SubStructure.
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IMolecule query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2SubStructure.
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

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);

        assertNotNull(smsd1.getFirstMapping());
    }

    /**
     * Test of set method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws InvalidSmilesException {
        System.out.println("set");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        VF2MCSHandler instance = new VF2MCSHandler();
        instance.set(query, target);
        instance.searchMCS(true, false);
        assertNotNull(instance.getFirstMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllMapping method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllMapping() throws InvalidSmilesException {
        System.out.println("getAllMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(4, smsd1.getAllMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstMapping method, of class VF2SubStructure.
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstMapping() throws InvalidSmilesException {
        System.out.println("getFirstMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(7, smsd1.getFirstMapping().size());
    }

    /**
     * Test of ring to ring match
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSystemMatch() throws InvalidSmilesException {
        System.out.println("testRingSystemMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)N1CCc2c(C1)cccc2)Br");
        IAtomContainer target = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)NCCOc1ccccc1)Br");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, true);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(24, smsd1.getFirstMapping().size());
        
        
        VF2MCSHandler smsd2 = new VF2MCSHandler();
        smsd2.set(query, target);
        smsd2.searchMCS(true, false);
        assertNotNull(smsd2.getFirstMapping());

        assertEquals(27, smsd2.getFirstMapping().size());
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

    /**
     * Test ring size match
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSizeMatch() throws InvalidSmilesException {
        System.out.println("testRingSizeMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C1=CC2=C3C(C=CC4=CC=CC(C=C2)=C34)=C1");
        IAtomContainer query = sp.parseSmiles("C1\\C=C/C=C/C=C\\C2=CC=CC(=C2)\\C=C/1");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, false);
        assertNotNull(smsd1.getFirstMapping());

        assertEquals(6, smsd1.getFirstMapping().size());
    }
    
    /**
     * Bug report by John Gerlits <jgerlits@utah.gov>
     * Cl should not match
     * Test ring size match
     * @throws InvalidSmilesException
     */
    @Test
    public void testVFMCSClMappingBugReportByJohn() throws InvalidSmilesException {
        System.out.println("testVFMCSClMappingBug");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3ccc(c4c3cccc4)Cl");
        IAtomContainer query = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3cccc4c3cccc4Cl");

        VF2MCSHandler smsd1 = new VF2MCSHandler();
        smsd1.set(query, target);
        smsd1.searchMCS(true, true);
        assertNotNull(smsd1.getFirstMapping());
        Assert.assertEquals(27, query.getAtomCount());
        Assert.assertEquals(27, target.getAtomCount());
        assertEquals(26, smsd1.getFirstMapping().size());
    }
}

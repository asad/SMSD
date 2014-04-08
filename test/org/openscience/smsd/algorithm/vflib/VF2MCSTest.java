
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
package org.openscience.smsd.algorithm.vflib;

import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.helper.MoleculeInitializer;
import org.openscience.smsd.interfaces.Algorithm;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * @cdk.module test-smsd @cdk.require java1.6+
 */
public class VF2MCSTest {

    public VF2MCSTest() {
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
     * Test of searchMCS method, of class VF2SubStructure.
     */
    @Test
    public void testSearchMCS() {
        //////System.out.println("1");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = null;
        try {
            target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(VF2MCSTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        IAtomContainer query = null;
        try {
            query = sp.parseSmiles("Nc1ccccc1");
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(VF2MCSTest.class.getName()).log(Level.SEVERE, null, ex);
        }

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws InvalidSmilesException {
        //////System.out.println("2");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());

    }

    /**
     * Test of set method, of class VF2SubStructure.
     *
     * @throws Exception
     */
    @Test
    public void testSet_IMolecule_IMolecule() throws Exception {
        //////System.out.println("3");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
    }

    /**
     * Test of set method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_AtomContainer_AtomContainer() throws InvalidSmilesException {
        //////System.out.println("4");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        Isomorphism instance = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(instance.getFirstAtomMapping());
    }

    /**
     * Test of getAllAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        //////System.out.println("5");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetAllAtomMapping() throws InvalidSmilesException {
        //////System.out.println("6");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        //////System.out.println("7");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetFirstAtomMapping() throws InvalidSmilesException {
        //////System.out.println("8");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of ring to ring match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSystemMatch() throws InvalidSmilesException, CDKException {
        //////System.out.println("9");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)N1CCc2c(C1)cccc2)Br");
        IAtomContainer target = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)NCCOc1ccccc1)Br");

        /*
         Set Ring perception ON
         */
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(24, smsd1.getFirstAtomMapping().getCount());

        VF2MCS smsd2 = new VF2MCS(query, target, true, false, false);
        assertNotNull(smsd2.getFirstAtomMapping());
        assertEquals(27, smsd2.getFirstAtomMapping().getCount());
    }

    /**
     * Test linker to linker match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testLinkersSystemMatch() throws InvalidSmilesException {
        //////System.out.println("10");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(5, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test linker to ring match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testLinkersVsRingsMatch() throws InvalidSmilesException, CDKException {
        //////System.out.println("11");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");

        /*
         Set Ring perception ON
         */
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(4, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test ring size match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSizeMatch() throws InvalidSmilesException {
        //////System.out.println("12");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C1=CC2=C3C(C=CC4=CC=CC(C=C2)=C34)=C1");
        IAtomContainer query = sp.parseSmiles("C1\\C=C/C=C/C=C\\C2=CC=CC(=C2)\\C=C/1");
        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        assertEquals(15, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Bug report by John Gerlits <jgerlits@utah.gov> Cl should not match Test
     * ring size match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testVFMCSClMappingBugReportByJohn() throws InvalidSmilesException {
        //////System.out.println("13");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3ccc(c4c3cccc4)Cl");
        IAtomContainer query = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3cccc4c3cccc4Cl");

        Isomorphism smsd1 = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd1.getFirstAtomMapping());
        Assert.assertEquals(27, query.getAtomCount());
        Assert.assertEquals(27, target.getAtomCount());
        assertEquals(26, smsd1.getFirstAtomMapping().getCount());
    }
}

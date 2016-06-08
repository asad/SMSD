
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
     * Test of set method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testSet_IAtomContainer_IAtomContainer() throws InvalidSmilesException {
        //System.out.println("1");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());

    }

    /**
     * Test of getAllAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException {
        //System.out.println("2");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(7, smsd.getFirstAtomMapping().getCount());
        assertEquals(4, smsd.getAllAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetAllAtomMapping() throws InvalidSmilesException {
        //System.out.println("3");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(4, smsd.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException {
        //System.out.println("4");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(7, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getFirstAtomMapping method, of class VF2SubStructure.
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testgetFirstAtomMapping() throws InvalidSmilesException {
        //System.out.println("5");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(7, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * Test of ring to ring match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSystemMatch() throws InvalidSmilesException, CDKException {
        //System.out.println("6");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)N1CCc2c(C1)cccc2)Br");
        IAtomContainer target = sp.parseSmiles("c1(ccc2c(c1)c(c([nH]2)C(=O)N)S(=O)(=O)N1CCOC(C1)C(=O)NCCOc1ccccc1)Br");

        /*
         Set Ring perception ON
         */
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(24, smsd.getFirstAtomMapping().getCount());

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
        //System.out.println("7");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(5, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * Test linker to ring match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testLinkersVsRingsMatch() throws InvalidSmilesException, CDKException {
        //System.out.println("8");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer query = sp.parseSmiles("OP(O)(=O)S1=CC=CC=C1");
        IAtomContainer target = sp.parseSmiles("OP(O)(S)=O");

        /*
         Set Ring perception ON
         */
        MoleculeInitializer.initializeMolecule(query);
        MoleculeInitializer.initializeMolecule(target);

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(4, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * Test ring size match
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testRingSizeMatch() throws InvalidSmilesException {
        //System.out.println("9");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C1=CC2=C3C(C=CC4=CC=CC(C=C2)=C34)=C1");
        IAtomContainer query = sp.parseSmiles("C1\\C=C/C=C/C=C\\C2=CC=CC(=C2)\\C=C/1");
        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, false);
        assertEquals(15, query.getAtomCount());
        assertEquals(16, target.getAtomCount());
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(6, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * Bug report by John Gerlits <jgerlits@utah.gov> Cl should not match Test
     * ring size match
     *
     * @BUG in the MCSPlus
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testVFMCSClMappingBugReportByJohn() throws InvalidSmilesException {
        //System.out.println("10");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3ccc(c4c3cccc4)Cl");
        IAtomContainer query = sp.parseSmiles("CCCCCn1c2c(cccc2)c(c1)C(=O)c3cccc4c3cccc4Cl");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, true, true, true);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(27, query.getAtomCount());
        assertEquals(27, target.getAtomCount());
        assertEquals(26, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * ComplexCase
     *
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testComplexCaseR03165() throws InvalidSmilesException {
        //System.out.println("11");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("OCC1=C(CC(O)=O)C(CCC(O)=O)=C(CC2=C(CC(O)=O)C(CCC(O)=O)=C(CC3=C(CC(O)=O)C(CCC(O)=O)=C(CC4=C(CC(O)=O)C(CCC(O)=O)=CN4)N3)N2)N1");
        IAtomContainer query = sp.parseSmiles("OC(=O)CCC1=C2CC3=C(CCC(O)=O)C(CC(O)=O)=C(CC4=C(CC(O)=O)C(CCC(O)=O)=C(CC5=C(CC(O)=O)C(CCC(O)=O)=C(CC(N2)=C1CC(O)=O)N5)N4)N3");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, false, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(60, query.getAtomCount());
        assertEquals(61, target.getAtomCount());
        assertEquals(55, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * ComplexCase
     *
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testComplexCaseR09087() throws InvalidSmilesException {
        //System.out.println("12");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("OP(O)(=O)O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, false, false, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(40, query.getAtomCount());
        assertEquals(36, target.getAtomCount());
        assertEquals(36, smsd.getFirstAtomMapping().getCount());
    }

    /**
     * ComplexCase
     *
     *
     * @throws InvalidSmilesException
     */
    @Test
    public void testMappingCoverage() throws InvalidSmilesException {
        //System.out.println("12");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC(=O)C=C");
        IAtomContainer target = sp.parseSmiles("CC1CC(CC=C1)C(C)=O");

        Isomorphism smsd = new Isomorphism(query, target, Algorithm.VFLibMCS, false, false, false);
        smsd.setChemFilters(true, true, false);
        assertNotNull(smsd.getFirstAtomMapping());
        assertEquals(5, query.getAtomCount());
        assertEquals(10, target.getAtomCount());
        assertEquals(5, smsd.getFirstAtomMapping().getCount());
    }
}

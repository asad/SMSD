/* 
 * Copyright (C) 2009-2014 Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import static org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator.createAnyAtomAnyBondContainer;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Unit testing for the {@link Substructure} class.
 *
 * @author Syed Asad Rahman
 * @author egonw
 * @cdk.module test-smsd
 */
public class SubstructureTest {

    /**
     * Tests if the Substructure can be instantiated without throwing
     * exceptions.
     *
     * @throws CDKException
     */
    @Test
    public void testSubStructureSearchAlgorithms() throws CDKException {
        Assert.assertNotNull(
                new Substructure(new AtomContainer(), new AtomContainer(),
                        true, false, false, true));
        Assert.assertNotNull(
                new Substructure(new QueryAtomContainer(DefaultChemObjectBuilder.getInstance()), new AtomContainer(), true));
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testInit_3args_1() throws InvalidSmilesException, CDKException {
//        //////System.out.println("init");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQuery());
        assertNotNull(smsd1.getTarget());
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testInit_3args_2() throws InvalidSmilesException, CDKException {
//        //////System.out.println("init");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQuery());
        assertNotNull(smsd1.getTarget());
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     *
     * @throws Exception
     */
    @Test
    public void testInit_3args_3() throws Exception {
//        //////System.out.println("init");
//        String sourceMolFileName = "";
//        String targetMolFileName = "";
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQuery());
        assertNotNull(smsd1.getTarget());

    }

    /**
     * Test of setChemFilters method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testSetChemFilters() throws InvalidSmilesException, CDKException {

//        //////System.out.println("setChemFilters");
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer mol = builder.newInstance(IAtomContainer.class);
        IAtom a1 = builder.newInstance(IAtom.class, "C");
        a1.setFormalCharge(0);
        mol.addAtom(a1);
        IAtom a2 = builder.newInstance(IAtom.class, "C");
        a2.setFormalCharge(0);
        mol.addAtom(a2);
        IAtom a3 = builder.newInstance(IAtom.class, "C");
        a3.setFormalCharge(0);
        mol.addAtom(a3);
        IAtom a4 = builder.newInstance(IAtom.class, "O");
        a4.setFormalCharge(0);
        mol.addAtom(a4);
        IAtom a5 = builder.newInstance(IAtom.class, "C");
        a5.setFormalCharge(0);
        mol.addAtom(a5);
        IAtom a6 = builder.newInstance(IAtom.class, "C");
        a6.setFormalCharge(0);
        mol.addAtom(a6);
        IAtom a7 = builder.newInstance(IAtom.class, "C");
        a7.setFormalCharge(0);
        mol.addAtom(a7);
        IAtom a8 = builder.newInstance(IAtom.class, "C");
        a8.setFormalCharge(0);
        mol.addAtom(a8);
        IBond b1 = builder.newInstance(IBond.class, a2, a1, IBond.Order.SINGLE);
        mol.addBond(b1);
        IBond b2 = builder.newInstance(IBond.class, a3, a2, IBond.Order.SINGLE);
        mol.addBond(b2);
        IBond b3 = builder.newInstance(IBond.class, a4, a3, IBond.Order.SINGLE);
        mol.addBond(b3);
        IBond b4 = builder.newInstance(IBond.class, a5, a4, IBond.Order.SINGLE);
        mol.addBond(b4);
        IBond b5 = builder.newInstance(IBond.class, a6, a5, IBond.Order.SINGLE);
        mol.addBond(b5);
        IBond b6 = builder.newInstance(IBond.class, a7, a6, IBond.Order.SINGLE);
        mol.addBond(b6);
        IBond b7 = builder.newInstance(IBond.class, a8, a6, IBond.Order.DOUBLE);
        mol.addBond(b7);
        IAtomContainer query = mol;

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, true, false, true, false);
        Assert.assertTrue(smsd1.isSubgraph());
        smsd1.setChemFilters(true, true, true);
        assertEquals(1, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getStereoScore method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetStereoScore() throws InvalidSmilesException, CDKException {
//        //////System.out.println("getStereoScore");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");
        IAtomContainer target = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        Assert.assertTrue(smsd1.isSubgraph());
        smsd1.setChemFilters(true, true, true);

        Integer score = 111;//as per new score
        assertEquals(score, smsd1.getStereoScore(0));
    }

    /**
     * Test of getFirstMapping method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetFirstMapping() throws InvalidSmilesException, CDKException {
//        //////System.out.println("getFirstMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        smsd1.setChemFilters(true, true, true);
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getAllAtomMapping method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testgetAllAtomMapping() throws InvalidSmilesException, CDKException {
//        //////System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("ONc1ccccc1");
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        Substructure smsd1 = new Substructure(query, target, true, true, true, true);
        smsd1.setChemFilters(false, true, true);
        assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class
     * SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException, CDKException {
//        //////System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false, false, false, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertTrue(smsd1.isSubgraph());
        assertEquals(7, smsd1.getFirstAtomMapping().getCount());
    }

    /**
     * Test of getAllAtomMapping method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException, CDKException {
//        //////System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//       query="Nc1ccccc1"
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IAtomContainer mol = builder.newInstance(IAtomContainer.class);
        IAtom a1 = builder.newInstance(IAtom.class, "N");
        a1.setFormalCharge(0);
        mol.addAtom(a1);
        IAtom a2 = builder.newInstance(IAtom.class, "C");
        a2.setFormalCharge(0);
        mol.addAtom(a2);
        IAtom a3 = builder.newInstance(IAtom.class, "C");
        a3.setFormalCharge(0);
        mol.addAtom(a3);
        IAtom a4 = builder.newInstance(IAtom.class, "C");
        a4.setFormalCharge(0);
        mol.addAtom(a4);
        IAtom a5 = builder.newInstance(IAtom.class, "C");
        a5.setFormalCharge(0);
        mol.addAtom(a5);
        IAtom a6 = builder.newInstance(IAtom.class, "C");
        a6.setFormalCharge(0);
        mol.addAtom(a6);
        IAtom a7 = builder.newInstance(IAtom.class, "C");
        a7.setFormalCharge(0);
        mol.addAtom(a7);
        IBond b1 = builder.newInstance(IBond.class, a2, a1, IBond.Order.SINGLE);
        mol.addBond(b1);
        IBond b2 = builder.newInstance(IBond.class, a3, a2, IBond.Order.SINGLE);
        mol.addBond(b2);
        IBond b3 = builder.newInstance(IBond.class, a4, a3, IBond.Order.SINGLE);
        mol.addBond(b3);
        IBond b4 = builder.newInstance(IBond.class, a5, a4, IBond.Order.SINGLE);
        mol.addBond(b4);
        IBond b5 = builder.newInstance(IBond.class, a6, a5, IBond.Order.SINGLE);
        mol.addBond(b5);
        IBond b6 = builder.newInstance(IBond.class, a7, a6, IBond.Order.SINGLE);
        mol.addBond(b6);
        IBond b7 = builder.newInstance(IBond.class, a7, a2, IBond.Order.SINGLE);
        mol.addBond(b7);

        IAtomContainer query = mol;

//      target = "C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C";
        builder = DefaultChemObjectBuilder.getInstance();
        mol = builder.newInstance(IAtomContainer.class);
        a1 = builder.newInstance(IAtom.class, "C");
        a1.setFormalCharge(0);
        mol.addAtom(a1);
        a2 = builder.newInstance(IAtom.class, "C");
        a2.setFormalCharge(0);
        mol.addAtom(a2);
        a3 = builder.newInstance(IAtom.class, "C");
        a3.setFormalCharge(0);
        mol.addAtom(a3);
        a4 = builder.newInstance(IAtom.class, "N");
        a4.setFormalCharge(0);
        mol.addAtom(a4);
        a5 = builder.newInstance(IAtom.class, "C");
        a5.setFormalCharge(0);
        mol.addAtom(a5);
        a6 = builder.newInstance(IAtom.class, "C");
        a6.setFormalCharge(0);
        mol.addAtom(a6);
        a7 = builder.newInstance(IAtom.class, "C");
        a7.setFormalCharge(0);
        mol.addAtom(a7);
        IAtom a8 = builder.newInstance(IAtom.class, "C");
        a8.setFormalCharge(0);
        mol.addAtom(a8);
        IAtom a9 = builder.newInstance(IAtom.class, "C");
        a9.setFormalCharge(0);
        mol.addAtom(a9);
        IAtom a10 = builder.newInstance(IAtom.class, "C");
        a10.setFormalCharge(0);
        mol.addAtom(a10);
        IAtom a11 = builder.newInstance(IAtom.class, "N");
        a11.setFormalCharge(0);
        mol.addAtom(a11);
        IAtom a12 = builder.newInstance(IAtom.class, "O");
        a12.setFormalCharge(0);
        mol.addAtom(a12);
        IAtom a13 = builder.newInstance(IAtom.class, "C");
        a13.setFormalCharge(0);
        mol.addAtom(a13);
        IAtom a14 = builder.newInstance(IAtom.class, "C");
        a14.setFormalCharge(0);
        mol.addAtom(a14);
        IAtom a15 = builder.newInstance(IAtom.class, "C");
        a15.setFormalCharge(0);
        mol.addAtom(a15);
        IAtom a16 = builder.newInstance(IAtom.class, "C");
        a16.setFormalCharge(0);
        mol.addAtom(a16);
        IAtom a17 = builder.newInstance(IAtom.class, "C");
        a17.setFormalCharge(0);
        mol.addAtom(a17);
        IAtom a18 = builder.newInstance(IAtom.class, "C");
        a18.setFormalCharge(0);
        mol.addAtom(a18);
        IAtom a19 = builder.newInstance(IAtom.class, "C");
        a19.setFormalCharge(0);
        mol.addAtom(a19);
        IAtom a20 = builder.newInstance(IAtom.class, "C");
        a20.setFormalCharge(0);
        mol.addAtom(a20);
        b1 = builder.newInstance(IBond.class, a2, a1, IBond.Order.SINGLE);
        mol.addBond(b1);
        b2 = builder.newInstance(IBond.class, a3, a2, IBond.Order.DOUBLE);
        mol.addBond(b2);
        b3 = builder.newInstance(IBond.class, a4, a3, IBond.Order.SINGLE);
        mol.addBond(b3);
        b4 = builder.newInstance(IBond.class, a5, a4, IBond.Order.SINGLE);
        mol.addBond(b4);
        b5 = builder.newInstance(IBond.class, a6, a5, IBond.Order.SINGLE);
        mol.addBond(b5);
        b6 = builder.newInstance(IBond.class, a7, a6, IBond.Order.SINGLE);
        mol.addBond(b6);
        b7 = builder.newInstance(IBond.class, a8, a7, IBond.Order.SINGLE);
        mol.addBond(b7);
        IBond b8 = builder.newInstance(IBond.class, a9, a8, IBond.Order.SINGLE);
        mol.addBond(b8);
        IBond b9 = builder.newInstance(IBond.class, a10, a9, IBond.Order.SINGLE);
        mol.addBond(b9);
        IBond b10 = builder.newInstance(IBond.class, a10, a5, IBond.Order.SINGLE);
        mol.addBond(b10);
        IBond b11 = builder.newInstance(IBond.class, a11, a9, IBond.Order.SINGLE);
        mol.addBond(b11);
        IBond b12 = builder.newInstance(IBond.class, a12, a11, IBond.Order.SINGLE);
        mol.addBond(b12);
        IBond b13 = builder.newInstance(IBond.class, a13, a11, IBond.Order.SINGLE);
        mol.addBond(b13);
        IBond b14 = builder.newInstance(IBond.class, a14, a13, IBond.Order.DOUBLE);
        mol.addBond(b14);
        IBond b15 = builder.newInstance(IBond.class, a15, a14, IBond.Order.SINGLE);
        mol.addBond(b15);
        IBond b16 = builder.newInstance(IBond.class, a16, a15, IBond.Order.SINGLE);
        mol.addBond(b16);
        IBond b17 = builder.newInstance(IBond.class, a17, a16, IBond.Order.DOUBLE);
        mol.addBond(b17);
        IBond b18 = builder.newInstance(IBond.class, a18, a17, IBond.Order.SINGLE);
        mol.addBond(b18);
        IBond b19 = builder.newInstance(IBond.class, a19, a18, IBond.Order.DOUBLE);
        mol.addBond(b19);
        IBond b20 = builder.newInstance(IBond.class, a20, a19, IBond.Order.SINGLE);
        mol.addBond(b20);

        IAtomContainer target = mol;

        Substructure smsd1 = new Substructure(query, target, true, false, true, true);
        Assert.assertTrue(smsd1.isSubgraph());
        smsd1.setChemFilters(true, true, true);

        assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getQuery method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetReactantMolecule() throws InvalidSmilesException, CDKException {
        //////System.out.println("getQuery");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true, false, true, false);
        smsd1.setChemFilters(true, true, true);
        assertEquals(7, smsd1.getQuery().getAtomCount());
    }

    /**
     * Test of getTarget method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetProductMolecule() throws InvalidSmilesException, CDKException {
        //////System.out.println("getTarget");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true, false, true, false);
        smsd1.setChemFilters(true, true, true);

        assertEquals(20, smsd1.getTarget().getAtomCount());
    }

    /**
     * Test of getTanimotoSimilarity method, of class
     * SubStructureSearchAlgorithms.
     *
     * @throws Exception
     */
    @Test
    public void testGetTanimotoSimilarity() throws Exception {
        //////System.out.println("getTanimotoSimilarity");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true, false, true, false);
        Assert.assertTrue(smsd1.isSubgraph());
        smsd1.setChemFilters(true, true, true);

        double score = 0.35;
        assertEquals(score, smsd1.getTanimotoSimilarity(), 0);
    }

    /**
     * Test of isStereoMisMatch method, of class SubStructureSearchAlgorithms.
     *
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testIsStereoMisMatch() throws InvalidSmilesException, CDKException {
        //////System.out.println("isStereoMisMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false, false, true, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertTrue(smsd1.isSubgraph());
        assertEquals(false, smsd1.isStereoMisMatch());
    }

    /**
     * Test of getEuclideanDistance method, of class
     * SubStructureSearchAlgorithms.
     *
     * @throws Exception
     */
    @Test
    public void testGetEuclideanDistance() throws Exception {
        //////System.out.println("getEuclideanDistance");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        double score = 3.605;
        Substructure smsd1 = new Substructure(query, target, true, false, true, false);
        smsd1.setChemFilters(true, true, true);
        //////System.out.println("smsd1.isSubgraph() " + smsd1.isSubgraph());
        Assert.assertTrue(smsd1.isSubgraph());
        assertEquals(score, smsd1.getEuclideanDistance(), 0.005);
    }

    @Test
    public void testImpossibleQuery() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C");
        Substructure smsd = new Substructure(query, target, true, false, true, false);
        boolean foundMatches = smsd.isSubgraph();
        Assert.assertFalse(foundMatches);
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond Sensitive.
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondSensitive() throws Exception {
        //////System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 0.0;
        Substructure smsd1 = new Substructure(molecule1, molecule2, true, false, true, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertFalse(smsd1.isSubgraph());
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond InSensitive.
     *
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondInSensitive() throws Exception {
        //////System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]"
                + "2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");
        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]"
                + "2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 1.0;
        Substructure smsd1 = new Substructure(molecule1, molecule2, false, false, true, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertTrue(smsd1.isSubgraph());
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }

    @Test
    public void testQueryAtomContainerDefault() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Substructure smsd = new Substructure(query, target, true, false, true, false);

        boolean foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Substructure(queryContainer, target, false);
        foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);
    }

    @Test
    public void testQueryAtomContainerSubstructure() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Substructure smsd = new Substructure(query, target, true, false, true, false);

        boolean foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd1 = new Substructure(queryContainer, target, false, true, true, false);
        Assert.assertTrue(smsd1.isSubgraph());
    }

    @Test
    public void testMatchAtomCount() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();

        Substructure smsd = new Substructure(query, target, true, false, false, true);
        Assert.assertTrue(smsd.isSubgraph());
//        for (AtomAtomMapping m : smsd.getAllAtomMapping()) {
//            //////System.out.println(m.getMappingsByIndex());
//        }
        Assert.assertEquals(18, smsd.getAllAtomMapping().size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd2 = new Substructure(queryContainer, target, true);
        Assert.assertTrue(smsd2.isSubgraph());
        Assert.assertEquals(18, smsd2.getAllAtomMapping().size());
    }

    @Test
    public void testIsSubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("[#8]P([#8])(=O)[#8]-[#6@H]-1-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H]-1-[#8]P([#8])([#8])=O");
        IAtomContainer target = sp.parseSmiles("[#8]P([#8])(=O)[#8]-[#6@@H]-1-[#6@H](-[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@H](-[#8]P([#8])(=O)[#8]P([#8])([#8])=O)-[#6@@H](-[#8]P([#8])([#8])=O)-[#6@@H]-1-[#8]P([#8])([#8])=O");
        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();

        Substructure smsd = new Substructure(query, target, true, false, true, true);
        Assert.assertTrue(smsd.isSubgraph());
//        for (AtomAtomMapping m : smsd.getAllAtomMapping()) {
//            //////System.out.println(m.getMappingsByIndex());
//        }
        Assert.assertEquals(26, smsd.getAllAtomMapping().size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd2 = new Substructure(queryContainer, target, false);
        Assert.assertTrue(smsd2.isSubgraph());
    }

    @Test
    public void testRingSubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String s = "";
        String t = "";

//      s="C[C@@]12CC3=C(CCC([O-])=O)[C@](C)(CC([O-])=O)C(\\C=C4/[NH2+][C@@](C)([C@@H]5N=C(\\C=C(/[NH2+]1)C(CCC([O-])=O)=C2CC([O-])=O)[C@](C)(CCC([O-])=O)[C@H]5CC([O-])=O)[C@@](C)(CC([O-])=O)[C@@H]4CCC([O-])=O)=N3");
//      t="C[C@@]12CC3=C(CCC([O-])=O)[C@](C)(CC([O-])=O)C(\\C=C4/N[C@@](C)(C5=C(CC([O-])=O)[C@@](C)(CCC([O-])=O)C(CC(=N1)C(CCC([O-])=O)=C2CC([O-])=O)=N5)[C@@](C)(CC([O-])=O)[C@@H]4CCC([O-])=O)=N3");
////        s="O=C([O-])CCC1C=2N=C(C(=C3[NH2+]C(C)(C4N=C(C(=C5N=C(C2)C(C)(C)C5CCC(=O)[O-])C)C(C)(CCC(=O)[O-])C4CC(=O)[O-])C(C)(CC(=O)[O-])C3CCC(=O)[O-])C)C1(C)CC(=O)[O-]");
////        t="O=C([O-])CCC1C=2N=C(C(=C3[NH2+]C(C)(C4N=C(C(=C5N=C(C2)C(C)(C)C5CCC(=O)[O-])C)C(C)(CCC(=O)[O-])C4CC(=O)[O-])C(C)(CC(=O)N)C3CCC(=O)[O-])C)C1(C)CC(=O)N");
//        s = "O=C([O-])CC1=C(C2=[N+]3C1(C)CC4=C(CCC(=O)[O-])C(C=5C=C6N7C(C)(C8[N+](=C(C2)C(C)(CCC(=O)[O-])C8CC(=O)[O-])[Co-2]73[N+]45)C(C)(CC(=O)[O-])C6CCC(=O)[O-])(C)CC(=O)[O-])CCC(=O)[O-]";
//        t = "O=C([O-])CC1=C(C2=[N+]3C1(C)CC4=C(CCC(=O)[O-])C(C=5C=C6N7C(C8=C(CC(=O)[O-])C(C(=[N+]8[Co-2]73[N+]45)C2)(C)CCC(=O)[O-])(C)C(C)(CC(=O)[O-])C6CCC(=O)[O-])(C)CC(=O)[O-])CCC(=O)[O-]";
//        
        s = "CC1=C2CC3=C(C)C(C=C)=C(CC4=C(C)C(C=C)=C(CC5=C(C)C(CCC([O-])=O)=C(CC(N2)=C1CCC([O-])=O)N5)N4)N3";
        t = "CC1=C2CC3=C(C)C(CCC([O-])=O)=C(CC4=C(CCC([O-])=O)C(C)=C(CC5=C(CCC([O-])=O)C(C)=C(CC(N2)=C1CCC([O-])=O)N5)N4)N3";

        IAtomContainer query = sp.parseSmiles(s);
        IAtomContainer target = sp.parseSmiles(t);

        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();

        Substructure smsd = new Substructure(query, target, false, true, false, true);
        Assert.assertTrue(smsd.isSubgraph());
        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd2 = new Substructure(queryContainer, target, false);
        Assert.assertFalse(smsd2.isSubgraph());
    }

    @Test
    public void testQuerySubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String s = "";
        String t = "";

        s = "CC1=C2CC3=C(C)C(C=C)=C(CC4=C(C)C(C=C)=C(CC5=C(C)C(CCC([O-])=O)=C(CC(N2)=C1CCC([O-])=O)N5)N4)N3";
        t = "CC1=C2CC3=C(C)C(CCC([O-])=O)=C(CC4=C(CCC([O-])=O)C(C)=C(CC5=C(CCC([O-])=O)C(C)=C(CC(N2)=C1CCC([O-])=O)N5)N4)N3";

        IAtomContainer query = sp.parseSmiles(s);
        IAtomContainer target = sp.parseSmiles(t);

        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }

        QueryAtomContainer createAnyAtomAnyBondContainer1 = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(query, false);

        for (i = 0; i < query.getAtomCount(); i++) {
            createAnyAtomAnyBondContainer1.getAtom(i).setID(query.getAtom(i).getID());
        }

        QueryAtomContainer createAnyAtomAnyBondContainer2 = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(target, false);
        for (i = 0; i < target.getAtomCount(); i++) {
            createAnyAtomAnyBondContainer2.getAtom(i).setID(target.getAtom(i).getID());
        }

//        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd = new Substructure(createAnyAtomAnyBondContainer1, target, true);
        ////System.out.println("Subgraph " + smsd.isSubgraph());
        Assert.assertTrue(smsd.isSubgraph());
    }

    @Test
    public void testisSubgraph() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        String s = "";
        String t = "";

        s = "OC1CCC2(C3=C(CCC2C1C)C4(C)CCC(C(C)CCC(=C)C(C)C)C4(C)CC3)C";
        t = "OC1CCC23CC43CCC5(C)C(CCC5(C)C4CCC2C1C)C(C)CCC(=C)C(C)C";

        IAtomContainer query = sp.parseSmiles(s);
        IAtomContainer target = sp.parseSmiles(t);

        int i = 1;
        for (IAtom a : query.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }
        //////System.out.println();
        i = 1;
        for (IAtom a : target.atoms()) {
            a.setID(String.valueOf(i));
            //System.out.print(i + ",");
            i++;
        }

        QueryAtomContainer createAnyAtomAnyBondContainer1 = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(query, false);

        for (i = 0; i < query.getAtomCount(); i++) {
            createAnyAtomAnyBondContainer1.getAtom(i).setID(query.getAtom(i).getID());
        }

        QueryAtomContainer createAnyAtomAnyBondContainer2 = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(target, false);
        for (i = 0; i < target.getAtomCount(); i++) {
            createAnyAtomAnyBondContainer2.getAtom(i).setID(target.getAtom(i).getID());
        }

//        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        Substructure smsd = new Substructure(createAnyAtomAnyBondContainer1, target, true);
        ////System.out.println("Subgraph " + smsd.isSubgraph());
        Assert.assertTrue(smsd.isSubgraph());
    }

    @Test
    public void testSolutionCount() throws Exception {
        ////System.out.println("33");
        String g1 = "NC1CCCCC1";
        String g2 = "NC1CCC(N)CC1";
        String g3 = "CNC1CCC(N)CC1";
        String g4 = "CNC1CCC(CC1)NC";

//        String g1 = "OC1CCCCC1";
//        String g2 = "OC1CCC(O)CC1";
//        String g3 = "COC1CCC(O)CC1";
//        String g4 = "COC1CCC(CC1)OC";
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer ac1 = smilesParser.parseSmiles(g1);
        IAtomContainer ac2 = smilesParser.parseSmiles(g2);
        IAtomContainer ac3 = smilesParser.parseSmiles(g3);
        IAtomContainer ac4 = smilesParser.parseSmiles(g4);

        Substructure overlap = new Substructure(ac1, ac2, true, false, false, true);
        overlap.setChemFilters(false, false, true);
        junit.framework.Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Substructure(ac1, ac3, true, false, false, true);
        overlap.setChemFilters(false, false, true);
        junit.framework.Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Substructure(ac1, ac4, true, false, false, true);
        overlap.setChemFilters(false, false, true);
        junit.framework.Assert.assertEquals(4, overlap.getAllAtomMapping().size());
        overlap = new Substructure(ac1, ac1, true, false, false, true);
        overlap.setChemFilters(false, false, true);
        junit.framework.Assert.assertEquals(2, overlap.getAllAtomMapping().size());

//        SmilesGenerator aromatic = SmilesGenerator.unique().aromatic();
//        ////System.out.println("SMILES Q :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnQuery()));
//        ////System.out.println("SMILES T :" + aromatic.create(overlap.getFirstAtomMapping().getMapCommonFragmentOnTarget()));
//        ////System.out.println("SMILES Common:" + overlap.getFirstAtomMapping().getCommonFragmentAsSMILES());
//
    }
}

/* 
 * Copyright (C) 2009-2011 Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.cdk.smsd;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Unit testing for the {@link Substructure} class.
 * @author     Syed Asad Rahman
 * @author     egonw
 * @cdk.module test-smsd
 */
public class SubstructureTest {

    /**
     * Tests if the CDKMCS can be instantiated without throwing exceptions.
     */
    @Test
    public void testSubStructureSearchAlgorithms() {
        Assert.assertNotNull(
                new Substructure(new AtomContainer(), new AtomContainer(),
                true));
        Assert.assertNotNull(
                new Substructure(new AtomContainer(), new AtomContainer(),
                false));
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testInit_3args_1() throws InvalidSmilesException, CDKException {
        System.out.println("init");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IMolecule target = sp.parseSmiles("C\\C=C/OCC=C");
        IMolecule query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQueryContainer());
        assertNotNull(smsd1.getTargetContainer());
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testInit_3args_2() throws InvalidSmilesException, CDKException {
        System.out.println("init");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQueryContainer());
        assertNotNull(smsd1.getTargetContainer());
    }

    /**
     * Test of init method, of class SubStructureSearchAlgorithms.
     * @throws Exception 
     */
    @Test
    public void testInit_3args_3() throws Exception {
        System.out.println("init");
//        String sourceMolFileName = "";
//        String targetMolFileName = "";
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/OCC=C");
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.setChemFilters(true, false, false);
        assertNotNull(smsd1.getQueryContainer());
        assertNotNull(smsd1.getTargetContainer());

    }

    /**
     * Test of setChemFilters method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testSetChemFilters() throws InvalidSmilesException, CDKException {

        System.out.println("setChemFilters");
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IMolecule mol = builder.newInstance(IMolecule.class);
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

        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.findSubgraphs();
        smsd1.setChemFilters(true, true, true);

        assertEquals(1, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getStereoScore method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetStereoScore() throws InvalidSmilesException, CDKException {
        System.out.println("getStereoScore");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CCCOCC(C)=C");
        IAtomContainer target = sp.parseSmiles("CCCOCC(C)=C");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.findSubgraphs();
        smsd1.setChemFilters(true, true, true);

        Integer score = 111;//as per new score
        assertEquals(score, smsd1.getStereoScore(0));
    }

    /**
     * Test of getFirstMapping method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetFirstMapping() throws InvalidSmilesException, CDKException {
        System.out.println("getFirstMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.findSubgraphs();
        smsd1.setChemFilters(true, true, true);
        assertEquals(7, smsd1.getFirstAtomMapping().size());
    }

    /**
     * Test of getAllMapping method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetAllMapping() throws InvalidSmilesException, CDKException {
        System.out.println("getAllMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");


        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.setChemFilters(true, true, true);
        smsd1.findSubgraphs();
        assertEquals(4, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getFirstAtomMapping method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetFirstAtomMapping() throws InvalidSmilesException, CDKException {
        System.out.println("getFirstAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.setChemFilters(true, true, true);
        smsd1.findSubgraphs();
        assertEquals(7, smsd1.getFirstAtomMapping().size());
    }

    /**
     * Test of getAllAtomMapping method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetAllAtomMapping() throws InvalidSmilesException, CDKException {
        System.out.println("getAllAtomMapping");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//       query="Nc1ccccc1"
        IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
        IMolecule mol = builder.newInstance(IMolecule.class);
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
        mol = builder.newInstance(IMolecule.class);
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

        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.findSubgraphs();
        smsd1.setChemFilters(true, true, true);

        assertEquals(2, smsd1.getAllAtomMapping().size());
    }

    /**
     * Test of getQueryContainer method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetReactantMolecule() throws InvalidSmilesException, CDKException {
        System.out.println("getQueryContainer");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.setChemFilters(true, true, true);
        assertEquals(7, smsd1.getQueryContainer().getAtomCount());
    }

    /**
     * Test of getTargetContainer method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testGetProductMolecule() throws InvalidSmilesException, CDKException {
        System.out.println("getTargetContainer");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.setChemFilters(true, true, true);

        assertEquals(20, smsd1.getTargetContainer().getAtomCount());
    }

    /**
     * Test of getTanimotoSimilarity method, of class SubStructureSearchAlgorithms.
     * @throws Exception
     */
    @Test
    public void testGetTanimotoSimilarity() throws Exception {
        System.out.println("getTanimotoSimilarity");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.findSubgraph();
        smsd1.setChemFilters(true, true, true);

        double score = 0.35;
        assertEquals(score, smsd1.getTanimotoSimilarity(), 0);
    }

    /**
     * Test of isStereoMisMatch method, of class SubStructureSearchAlgorithms.
     * @throws InvalidSmilesException
     * @throws CDKException
     */
    @Test
    public void testIsStereoMisMatch() throws InvalidSmilesException, CDKException {
        System.out.println("isStereoMisMatch");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        Substructure smsd1 = new Substructure(query, target, false);
        smsd1.setChemFilters(true, true, true);
        smsd1.findSubgraphs();
        assertEquals(false, smsd1.isStereoMisMatch());
    }

    /**
     * Test of getEuclideanDistance method, of class SubStructureSearchAlgorithms.
     * @throws Exception
     */
    @Test
    public void testGetEuclideanDistance() throws Exception {
        System.out.println("getEuclideanDistance");
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C");
        IAtomContainer query = sp.parseSmiles("Nc1ccccc1");

        double score = 3.605;
        Substructure smsd1 = new Substructure(query, target, true);
        smsd1.setChemFilters(true, true, true);
        smsd1.findSubgraph();
        assertEquals(score, smsd1.getEuclideanDistance(), 0.005);
    }

    @Test
    public void testImpossibleQuery() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C");
        Substructure smsd = new Substructure(query, target, true);
        boolean foundMatches = smsd.findSubgraph();
        Assert.assertFalse(foundMatches);
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond Sensitive.
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondSensitive() throws Exception {
        System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 0.0;
        Substructure smsd1 = new Substructure(molecule1, molecule2, true);
        smsd1.setChemFilters(true, true, true);
        Assert.assertFalse(smsd1.findSubgraph());
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }

    /**
     * Test Tanimoto NAD+ & NADH for Bond InSensitive.
     * @throws Exception
     */
    @Test
    public void testNADPlusNADHBondInSensitive() throws Exception {
        System.out.println("getTanimoto for NAD+ and NADH");
        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles("NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        IAtomContainer molecule2 = smilesParser.parseSmiles("NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O");

        double score = 1.0;
        Substructure smsd1 = new Substructure(molecule1, molecule2, false);
        smsd1.setChemFilters(true, true, true);
        Assert.assertTrue(smsd1.findSubgraph());
        Assert.assertEquals(score, smsd1.getTanimotoSimilarity(), 0.001);
    }
}

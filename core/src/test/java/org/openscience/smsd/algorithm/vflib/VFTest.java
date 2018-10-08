/* 
 * Copyright (C) 2009-2014 Syed Asad Rahman <s9asad@gmail.com>
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

import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import static org.junit.Assert.assertEquals;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.algorithm.vflib.vf2.sub.Pattern;
import org.openscience.smsd.algorithm.vflib.vf2.sub.VF;

/**
 * Unit testing for the {@link VF} class.
 *
 * @author Syed Asad Rahman <s9asad@gmail.com>
 * @author egonw test-smsd
 */
public class VFTest {

    public VFTest() {
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

    @Test
    public void testQueryAtomContainerDefault() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");
        Pattern findIdentical = VF.findSubstructure(query, true, true, true);

        boolean foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        findIdentical = VF.findSubstructure(queryContainer);
        foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);
    }

    @Test
    public void testMatchDefault() throws CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");
        IAtomContainer target = sp.parseSmiles("OP(O)(=O)O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");
        Pattern findIdentical = VF.findSubstructure(query, false, false, false);

        boolean foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        findIdentical = VF.findSubstructure(queryContainer);
        foundMatches = findIdentical.matches(target);
        Assert.assertTrue(foundMatches);
    }

    @Test
    /*
     * Imp test case
     */
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

        Pattern findIdentical = VF.findSubstructure(query, true, true, true);
        List<Map<IAtom, IAtom>> matchAll = findIdentical.matchAll(target);
        Assert.assertTrue(!matchAll.isEmpty());
        Assert.assertEquals(18, matchAll.size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        findIdentical = VF.findSubstructure(queryContainer);
        matchAll = findIdentical.matchAll(target);
        Assert.assertEquals(18, matchAll.size());
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
//        Pattern pattern = Ullmann.findSubstructure(query);
        Pattern pattern = VF.findSubstructure(query, true, true, true);
        List<Map<IAtom, IAtom>> matchAll = pattern.matchAll(target);
        Assert.assertEquals(768, matchAll.size());

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        pattern = VF.findSubstructure(queryContainer);
        Assert.assertTrue(pattern.matches(target));
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
        IAtomContainer query = sp.parseSmiles("OP(O)(=O)O[C@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");
        IAtomContainer target = sp.parseSmiles("OP(O)(=O)O[C@@H]1[C@H](OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@H](OP(O)(=O)OP(O)(O)=O)[C@@H](OP(O)(O)=O)[C@@H]1OP(O)(O)=O");

        Pattern find = VF.findSubstructure(query, false, false, false);
        List<Map<IAtom, IAtom>> matchAll = find.matchAll(target);
        Assert.assertTrue(!matchAll.isEmpty());
        assertEquals(36, query.getAtomCount());
        assertEquals(40, target.getAtomCount());
        assertEquals(559872, matchAll.size());
        assertEquals(36, matchAll.iterator().next().size());
    }

}

/* 
 * Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.mcsplus.MCSPlusHandler;
import org.openscience.smsd.algorithm.rgraph.CDKMCSHandler;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.smsd.algorithm.vflib.VF2MCSPlusHandler;
import org.openscience.smsd.global.TimeOut;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.interfaces.ITimeOut;

/**
 *  <p>This class implements the Isomorphism- a multipurpose structure comparison tool.
 *  It allows users to, i) find the maximal common substructure(s) (MCS);
 *  ii) perform the mapping of a substructure in another structure, and;
 *  iii) map two isomorphic structures.</p>
 *
 *  <p>It also comes with various published algorithms. The user is free to
 *  choose his favorite algorithm to perform MCS or substructure search.
 *  For example:</p>
 * <OL>
 * <lI>0: Default,
 * <lI>1: Default_1,
 * <lI>2: MCSPlus,
 * <lI>3: VFLibMCS,
 * <lI>4: CDKMCS
 * </OL>
 *  <p>It also has a set of robust chemical filters (i.e. bond energy, fragment
 *  count, stereo & bond match) to sort the reported MCS solutions in a chemically
 *  relevant manner. Each comparison can be made with or without using the bond
 *  sensitive mode and with implicit or explicit hydrogens.</p>
 *
 *  <p>If you are using <font color="#FF0000">Isomorphism, please cite Rahman <i>et.al. 2009</i></font>
 *  {@cdk.cite SMSD2009}. The Isomorphism algorithm is described in this paper.
 *  </p>
 *
 * <p>An example for <b>MCS search</b>:</p>
 *  <font color="#003366">
 *  <pre>
 *  SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 *  // Benzene
 *  IAtomContainer query = sp.parseSmiles("C1=CC=CC=C1");
 *  // Napthalene
 *  IAtomContainer target = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
 *  //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
 *  //Bond Sensitive is set true
 * 
 *  // set molecules and/or remove hydrogens
 *  Isomorphism comparison = new Isomorphism(query, target, Algorithm.VFLibMCS, true);
 *  // set chemical filter true
 *  comparison.setChemFilters(true, true, true);
 *
 *  //Get similarity score
 *   System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 *  // Print the mapping between molecules
 *   System.out.println(" Mappings: ");
 *   for (Map.Entry<IAtom, IAtom> mapping : comparison.getMappings().entrySet()) {
 *      IAtom eAtom = query.getKey();
 *      IAtom pAtom = target.getValue();
 *      System.out.println(eAtom.getSymbol() + " " + pAtom.getSymbol());
 *   }
 *   System.out.println("");
 *  
 *
 *  </pre>
 *  </font>
 *
 * @cdk.require java1.5+
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
@TestClass("org.openscience.cdk.smsd.factory.SubStructureSearchAlgorithmsTest")
public final class Isomorphism extends BaseMapping implements ITimeOut, Serializable {

    static final long serialVersionUID = 10278639972837495L;
    private Algorithm algorithmType;
    private final ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Isomorphism.class);
    private double bondSensitiveCDKMCSTimeOut = 1.00;//mins
    private double bondInSensitiveCDKMCSTimeOut = 1.00;//mins
    private double bondSensitiveMCSPlusTimeOut = 1.00;//mins
    private double bondInSensitiveMCSPlusTimeOut = 2.00;//mins
    private double bondSensitiveVFTimeOut = 2.00;//mins
    private double bondInSensitiveVFTimeOut = 2.00;//mins
    private boolean matchBonds;
    private boolean matchRings;

    /**
     * Initialize query and target molecules.
     *
     * Note: Here its assumed that hydrogens are implicit
     * and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector 
     * before initializing calling this method.
     * 
     * @param query query molecule
     * @param target target molecule
     * This is the algorithm factory and entry port for all the MCS algorithm in the Isomorphism
     * supported algorithm {@link org.openscience.cdk.smsd.interfaces.Algorithm} types:
     * <OL>
     * <lI>0: Default,
     * <lI>1: Default_1,
     * <lI>2: MCSPlus,
     * <lI>3: VFLibMCS,
     * <lI>4: CDKMCS
     * </OL>
     * @param algorithmType {@link org.openscience.cdk.smsd.interfaces.Algorithm}
     */
    @TestMethod("testIsomorphismTest")
    public Isomorphism(
            IQueryAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType) {
        this.algorithmType = algorithmType;
        this.mol1 = query;
        this.mol2 = target;
        this.mcsList = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        mcsBuilder(mol1, mol2);
    }

    /**
     * Initialize query and target molecules.
     *
     * Note: Here its assumed that hydrogens are implicit
     * and user has called these two methods
     * percieveAtomTypesAndConfigureAtoms and CDKAromicityDetector 
     * before initializing calling this method.
     * 
     * @param query query mol
     * @param target target mol
     * This is the algorithm factory and entry port for all the MCS algorithm in the Isomorphism
     * supported algorithm {@link org.openscience.cdk.smsd.interfaces.Algorithm} types:
     * <OL>
     * <lI>0: Default,
     * <lI>1: Default_1,
     * <lI>2: MCSPlus,
     * <lI>3: VFLibMCS,
     * <lI>4: CDKMCS
     * </OL>
     * @param algorithmType {@link org.openscience.cdk.smsd.interfaces.Algorithm}
     * @param bondTypeFlag Match bond types (i.e. double to double etc)
     * @param matchRings Match ring atoms and ring size 
     */
    @TestMethod("testIsomorphismTest")
    public Isomorphism(
            IAtomContainer query,
            IAtomContainer target,
            Algorithm algorithmType,
            boolean bondTypeFlag,
            boolean matchRings) {
        this.algorithmType = algorithmType;
        this.mol1 = query;
        this.mol2 = target;
        if (matchRings) {
            try {
                initializeMolecule(mol1);
                initializeMolecule(mol2);
            } catch (CDKException ex) {
                Logger.error(Level.SEVERE, null, ex);
            }
        }
        this.mcsList = Collections.synchronizedList(new ArrayList<AtomAtomMapping>());
        this.matchBonds = bondTypeFlag;
        this.matchRings = matchRings;
        setTime(bondTypeFlag);
        mcsBuilder(mol1, mol2);
    }

    private synchronized void mcsBuilder(IAtomContainer mol1, IAtomContainer mol2) {

        int rBondCount = mol1.getBondCount();
        int pBondCount = mol2.getBondCount();

        int rAtomCount = mol1.getAtomCount();
        int pAtomCount = mol2.getAtomCount();

        if ((rBondCount == 0 && rAtomCount > 0) || (pBondCount == 0 && pAtomCount > 0)) {
            singleMapping();
        } else {
            chooseAlgorithm();
        }
    }

    /**
     * Returns bond maps between sourceAtomCount and targetAtomCount molecules based on the atoms
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mappings mappings between sourceAtomCount and targetAtomCount molecule atoms
     * @return bond maps between sourceAtomCount and targetAtomCount molecules based on the atoms
     */
    public synchronized List<Map<IBond, IBond>> makeBondMapsOfAtomMaps(IAtomContainer ac1,
            IAtomContainer ac2, List<AtomAtomMapping> mappings) {
        List<Map<IBond, IBond>> bondMaps = Collections.synchronizedList(new ArrayList<Map<IBond, IBond>>());
        for (AtomAtomMapping mapping : mappings) {
            bondMaps.add(makeBondMapOfAtomMap(ac1, ac2, mapping));
        }
        return bondMaps;
    }

    /**
     *
     * Returns bond map between sourceAtomCount and targetAtomCount molecules based on the atoms
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mapping mappings between sourceAtomCount and targetAtomCount molecule atoms
     * @return bond map between sourceAtomCount and targetAtomCount molecules based on the atoms
     */
    private synchronized Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2,
            AtomAtomMapping mapping) {

        Map<IBond, IBond> bondbondMappingMap = Collections.synchronizedMap(new HashMap<IBond, IBond>());

        for (Map.Entry<IAtom, IAtom> map1 : mapping.getMappings().entrySet()) {
            for (Map.Entry<IAtom, IAtom> map2 : mapping.getMappings().entrySet()) {
                if (map1.getKey() != map2.getKey()) {
                    IBond bond1 = ac1.getBond(map1.getKey(), map2.getKey());
                    IBond bond2 = ac2.getBond(map1.getValue(), map2.getValue());
                    if (bond1 != null && bond2 != null && !bondbondMappingMap.containsKey(bond1)) {
                        bondbondMappingMap.put(bond1, bond2);
                    }
                }
            }
        }
//        System.out.println("Mol Map size:" + bondbondMappingMap.size());
        return bondbondMappingMap;
    }

    private synchronized void chooseAlgorithm() {

        switch (algorithmType) {
            case CDKMCS:
                cdkMCSAlgorithm();
                break;
            case DEFAULT:
                defaultMCSAlgorithm();
                break;
            case MCSPlus:
                mcsPlusAlgorithm();
                break;
            case VFLibMCS:
                vfLibMCSAlgorithm();
                break;
        }
    }

    private synchronized void cdkMCSAlgorithm() {
//        double time = System.currentTimeMillis();
        CDKMCSHandler mcs = null;
        mcs = new CDKMCSHandler();
        mcs.set(mol1, mol2);
        mcs.searchMCS(isMatchBonds(), isMatchRings());
        clearMaps();
        mcsList.addAll(mcs.getAllAtomMapping());
//        System.out.println("\nCDKMCS used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));
    }

    private synchronized void mcsPlusAlgorithm() {
//        double time = System.currentTimeMillis();
        MCSPlusHandler mcs = null;
        mcs = new MCSPlusHandler();
        mcs.set(mol1, mol2);
        mcs.searchMCS(isMatchBonds(), isMatchRings());
        clearMaps();
        mcsList.addAll(mcs.getAllAtomMapping());
//        System.out.println("\nMCSPlus used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));
    }

    private synchronized void vfLibMCSAlgorithm() {
//        double time = System.currentTimeMillis();
        VF2MCSPlusHandler mcs = null;
        mcs = new VF2MCSPlusHandler();
        mcs.set(mol1, mol2);
        mcs.searchMCS(isMatchBonds(), isMatchRings());
        clearMaps();
        mcsList.addAll(mcs.getAllAtomMapping());
//        System.out.println("\nVFLib used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));
    }

    private synchronized void singleMapping() {
        SingleMappingHandler mcs = null;
        mcs = new SingleMappingHandler();
        mcs.set(mol1, mol2);
        mcs.searchMCS(isMatchBonds(), isMatchRings());
        clearMaps();
        mcsList.addAll(mcs.getAllAtomMapping());
    }

    private synchronized void defaultMCSAlgorithm() {
        try {
            cdkMCSAlgorithm();
            if (getMappingCount() == 0 || isTimeOut()) {
//                System.out.println("\nCDKMCS hit by timeout\n");
//                double time = System.currentTimeMillis();
                vfLibMCSAlgorithm();
//                System.out.println("\nVF Lib used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private synchronized void setTime(boolean bondTypeFlag) {
        if (bondTypeFlag) {
            TimeOut tmo = TimeOut.getInstance();
            tmo.setCDKMCSTimeOut(getBondSensitiveCDKMCSTimeOut());
            tmo.setMCSPlusTimeout(getBondSensitiveMCSPlusTimeOut());
            tmo.setVFTimeout(getBondSensitiveVFTimeOut());
        } else {
            TimeOut tmo = TimeOut.getInstance();
            tmo.setCDKMCSTimeOut(getBondInSensitiveCDKMCSTimeOut());
            tmo.setMCSPlusTimeout(getBondInSensitiveMCSPlusTimeOut());
            tmo.setVFTimeout(getBondInSensitiveVFTimeOut());
        }
    }

    public synchronized boolean isTimeOut() {
        return TimeOut.getInstance().isTimeOutFlag();
    }

    public synchronized void resetTimeOut() {
        TimeOut.getInstance().setTimeOutFlag(false);
    }

    private synchronized void clearMaps() {
        this.mcsList.clear();
    }

    /** 
     *
     * @return true if query is a subgraph of the target
     */
    @TestMethod("testIsSubgraph")
    public synchronized boolean isSubgraph() {

        float mappingSize = 0;
        if (getMappingCount() > 0) {
            mappingSize = getAllAtomMapping().iterator().next().getCount();
        } else {
            return false;
        }
        int sourceAtomCount = mol1.getAtomCount();
        int targetAtomCount = mol2.getAtomCount();

        if (mappingSize == sourceAtomCount && mappingSize <= targetAtomCount) {
            if (mappingSize == 1) {
                return true;
            } else if (!getAllBondMaps().isEmpty()
                    && getAllBondMaps().iterator().next().size() == mol1.getBondCount()) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return the matchBonds
     */
    public synchronized boolean isMatchBonds() {
        return matchBonds;
    }

    /**
     * @return the allBondMCS
     */
    public synchronized List<Map<IBond, IBond>> getAllBondMaps() {
        if (!mcsList.isEmpty()) {
            return makeBondMapsOfAtomMaps(mol1, mol2, mcsList);
        }
        return new ArrayList<Map<IBond, IBond>>();
    }

    @Override
    public synchronized double getBondSensitiveCDKMCSTimeOut() {
        return this.bondSensitiveCDKMCSTimeOut;
    }

    @Override
    public synchronized void setBondSensitiveCDKMCSTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveCDKMCSTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public synchronized double getBondInSensitiveCDKMCSTimeOut() {
        return bondInSensitiveCDKMCSTimeOut;
    }

    @Override
    public synchronized void setBondInSensitiveCDKMCSTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveCDKMCSTimeOut = bondInSensitiveTimeOut;
    }

    @Override
    public synchronized double getBondSensitiveMCSPlusTimeOut() {
        return this.bondSensitiveMCSPlusTimeOut;
    }

    @Override
    public synchronized void setBondSensitiveMCSPlusTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveMCSPlusTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public synchronized double getBondInSensitiveMCSPlusTimeOut() {
        return this.bondInSensitiveMCSPlusTimeOut;
    }

    @Override
    public synchronized void setBondInSensitiveMCSPlusTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveMCSPlusTimeOut = bondInSensitiveTimeOut;
    }

    @Override
    public synchronized double getBondSensitiveVFTimeOut() {
        return this.bondSensitiveVFTimeOut;
    }

    @Override
    public synchronized void setBondSensitiveVFTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveVFTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public synchronized double getBondInSensitiveVFTimeOut() {
        return this.bondInSensitiveVFTimeOut;
    }

    @Override
    public synchronized void setBondInSensitiveVFTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveVFTimeOut = bondInSensitiveTimeOut;
    }

    /**
     * @return the matchRings
     */
    public boolean isMatchRings() {
        return matchRings;
    }

    /**
     * @param matchRings the matchRings to set
     */
    public void setMatchRings(boolean matchRings) {
        this.matchRings = matchRings;
    }
}

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

import java.io.IOException;
import java.io.Serializable;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
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
import org.openscience.smsd.algorithm.vflib.VFlibMCSHandler;
import org.openscience.smsd.filters.ChemicalFilters;
import org.openscience.smsd.global.TimeOut;
import org.openscience.smsd.interfaces.AbstractMCS;
import org.openscience.smsd.interfaces.Algorithm;

/**
 *  <p>This class implements the Isomorphism- a multipurpose structure comparison tool.
 *  It allows users to, i) find the maximal common substructure(s) (MCS);
 *  ii) perform the mapping of a substructure in another structure, and;
 *  iii) map two isomorphic structures.</p>
 *
 *  <p>It also comes with various published algorithms. The user is free to
 *  choose his favorite algorithm to perform MCS or substructure search.
 *  For example 0: Isomorphism algorithm, 1: MCSPlus, 2: VFLibMCS, 3: CDKMCS, 4:
 *  Substructure</p>
 *
 *  <p>It also has a set of robust chemical filters (i.e. bond energy, fragment
 *  count, stereo & bond match) to sort the reported MCS solutions in a chemically
 *  relevant manner. Each comparison can be made with or without using the bond
 *  sensitive mode and with implicit or explicit hydrogens.</p>
 *
 *  <p>If you are using <font color="#FF0000">Isomorphism, please cite Rahman <i>et.al. 2009</i></font>
 *  {@cdk.cite SMSD2009}. The Isomorphism algorithm is described in this paper.
 *  </p>
 *
 *
 * <p>An example for <b>Substructure search</b>:</p>
 *  <font color="#003366">
 *  <pre>
 *  SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 *  // Benzene
 *  IAtomContainer A1 = sp.parseSmiles("C1=CC=CC=C1");
 *  // Napthalene
 *  IAtomContainer A2 = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
 *  //Turbo mode search
 *  //Bond Sensitive is set true
 *  Isomorphism comparison = new Isomorphism(Algorithm.Substructure, true);
 *  // set molecules, remove hydrogens, clean and configure molecule
 *  comparison.init(A1, A2, true, true);
 *  // set chemical filter true
 *  comparison.setChemFilters(false, false, false);
 *  if (comparison.isSubgraph()) {
 *  //Get similarity score
 *   System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 *   System.out.println("A1 is a subgraph of A2:  " + comparison.isSubgraph());
 *  //Get Modified GraphAtomContainer
 *   IAtomContainer Mol1 = comparison.getQueryMolecule();
 *   IAtomContainer Mol2 = comparison.getTargetMolecule();
 *  // Print the mapping between molecules
 *   System.out.println(" Mappings: ");
 *   for (Map.Entry <Integer, Integer> mapping : comparison.getFirstMapping().entrySet()) {
 *      System.out.println((mapping.getKey() + 1) + " " + (mapping.getValue() + 1));
 *
 *      IAtom eAtom = Mol1.getAtom(mapping.getKey());
 *      IAtom pAtom = Mol2.getAtom(mapping.getValue());
 *      System.out.println(eAtom.getSymbol() + " " + pAtom.getSymbol());
 *   }
 *   System.out.println("");
 *  }
 *
 *  </pre>
 *  </font>
 *
 * <p>An example for <b>MCS search</b>:</p>
 *  <font color="#003366">
 *  <pre>
 *  SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
 *  // Benzene
 *  IAtomContainer A1 = sp.parseSmiles("C1=CC=CC=C1");
 *  // Napthalene
 *  IAtomContainer A2 = sp.parseSmiles("C1=CC2=C(C=C1)C=CC=C2");
 *  //{ 0: Default Isomorphism Algorithm, 1: MCSPlus Algorithm, 2: VFLibMCS Algorithm, 3: CDKMCS Algorithm}
 *  //Bond Sensitive is set true
 *  Isomorphism comparison = new Isomorphism(Algorithm.DEFAULT, true);
 *  // set molecules, remove hydrogens, clean and configure molecule
 *  comparison.init(A1, A2, true, true);
 *  // set chemical filter true
 *  comparison.setChemFilters(true, true, true);
 *
 *  //Get similarity score
 *  System.out.println("Tanimoto coefficient:  " + comparison.getTanimotoSimilarity());
 *  System.out.println("A1 is a subgraph of A2:  " + comparison.isSubgraph());
 *  //Get Modified GraphAtomContainer
 *  IAtomContainer Mol1 = comparison.getQueryMolecule();
 *  IAtomContainer Mol2 = comparison.getTargetMolecule();
 *  // Print the mapping between molecules
 *  System.out.println(" Mappings: ");
 *  for (Map.Entry <Integer, Integer> mapping : comparison.getFirstMapping().entrySet()) {
 *      System.out.println((mapping.getKey() + 1) + " " + (mapping.getValue() + 1));
 *
 *      IAtom eAtom = Mol1.getAtom(mapping.getKey());
 *      IAtom pAtom = Mol2.getAtom(mapping.getValue());
 *      System.out.println(eAtom.getSymbol() + " " + pAtom.getSymbol());
 *  }
 *  System.out.println("");
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
public final class Isomorphism extends AbstractMCS implements IAtomAtomMapping, Serializable {

    static final long serialVersionUID = 10278639972837495L;
    private List<Map<Integer, Integer>> allMCS = null;
    private Map<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private List<Map<IBond, IBond>> allBondMCS = null;
    private Map<IBond, IBond> firstBondMCS = null;
    private IAtomContainer target = null;
    private IAtomContainer query = null;
    private List<Double> stereoScore = null;
    private List<Integer> fragmentSize = null;
    private List<Double> bEnergies = null;
    private Algorithm algorithmType;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Isomorphism.class);
    private double bondSensitiveCDKMCSTimeOut = 1.00;//mins
    private double bondInSensitiveCDKMCSTimeOut = 2.00;//mins
    private double bondSensitiveMCSPlusTimeOut = 1.00;//mins
    private double bondInSensitiveMCSPlusTimeOut = 2.00;//mins
    private double bondSensitiveVFTimeOut = 2.00;//mins
    private double bondInSensitiveVFTimeOut = 5.00;//mins
    private boolean matchBonds = false;

    /**
     * This is the algorithm factory and entry port for all the MCS algorithm in the Isomorphism
     * supported algorithm {@link org.openscience.cdk.smsd.interfaces.Algorithm} types:
     * <OL>
     * <lI>0: Default,
     * <lI>1: MCSPlus,
     * <lI>2: VFLibMCS,
     * <lI>3: CDKMCS,
     * <lI>4: Substructure
     * </OL>
     * @param algorithmType {@link org.openscience.cdk.smsd.interfaces.Algorithm}
     * @param bondTypeFlag
     */
    @TestMethod("testSubStructureSearchAlgorithms")
    public Isomorphism(Algorithm algorithmType, boolean bondTypeFlag) {
        this.algorithmType = algorithmType;
        firstSolution = new TreeMap<Integer, Integer>();
        allMCS = new ArrayList<Map<Integer, Integer>>();
        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        firstAtomMCS = new HashMap<IAtom, IAtom>();
        allBondMCS = new ArrayList<Map<IBond, IBond>>();
        firstBondMCS = new HashMap<IBond, IBond>();
        setTime(bondTypeFlag);
        setMatchBonds(bondTypeFlag);
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

        if (!allAtomMCS.isEmpty() && !firstAtomMCS.isEmpty() && firstAtomMCS.size() > 1) {
            setAllBondMaps(makeBondMapsOfAtomMaps(mol1, mol2, allAtomMCS));
            if (getAllBondMaps().iterator().hasNext()) {
                setFirstBondMap(getAllBondMaps().iterator().next());
            }
        }
    }

    private synchronized void mcsBuilder(IQueryAtomContainer mol1, IAtomContainer mol2) {

        int rBondCount = mol1.getBondCount();
        int pBondCount = mol2.getBondCount();

        int rAtomCount = mol1.getAtomCount();
        int pAtomCount = mol2.getAtomCount();

        if ((rBondCount == 0 && rAtomCount > 0) || (pBondCount == 0 && pAtomCount > 0)) {
            singleMapping();
        } else {
            chooseAlgorithm();
        }

        if (!allAtomMCS.isEmpty() && !firstAtomMCS.isEmpty() && firstAtomMCS.size() > 1) {
            setAllBondMaps(makeBondMapsOfAtomMaps(mol1, mol2, allAtomMCS));
            if (getAllBondMaps().iterator().hasNext()) {
                setFirstBondMap(getAllBondMaps().iterator().next());
            }
        }
    }

    /**
     * Returns bond maps between sourceAtomCount and targetAtomCount molecules based on the atoms
     * @param ac1 sourceAtomCount molecule
     * @param ac2 targetAtomCount molecule
     * @param mappings mappings between sourceAtomCount and targetAtomCount molecule atoms
     * @return bond maps between sourceAtomCount and targetAtomCount molecules based on the atoms
     */
    public static List<Map<IBond, IBond>> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, List<Map<IAtom, IAtom>> mappings) {
        List<Map<IBond, IBond>> bondMaps = new ArrayList<Map<IBond, IBond>>();
        for (Map<IAtom, IAtom> mapping : mappings) {
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
    public static Map<IBond, IBond> makeBondMapOfAtomMap(IAtomContainer ac1, IAtomContainer ac2, Map<IAtom, IAtom> mapping) {
        Map<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (Map.Entry<IAtom, IAtom> mapS : mapping.entrySet()) {
            IAtom indexI = mapS.getKey();
            IAtom indexJ = mapS.getValue();

            for (Map.Entry<IAtom, IAtom> mapD : mapping.entrySet()) {
                IAtom indexIPlus = mapD.getKey();
                IAtom indexJPlus = mapD.getValue();

                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {
                    IBond ac1Bond = ac1.getBond(indexI, indexIPlus);
                    if (ac1Bond != null) {
                        IBond ac2Bond = ac2.getBond(indexJ, indexJPlus);
                        if (ac2Bond != null) {
                            maps.put(ac1Bond, ac2Bond);
                        }
                    }
                }
            }
        }

//        System.out.println("bond Map size:" + maps.size());

        return maps;

    }

    private void chooseAlgorithm() {

        switch (algorithmType) {
            case CDKMCS:
                cdkMCSAlgorithm();
                break;
            case DEFAULT:
                defaultMCSAlgorithm();
                break;
            case DEFAULT_1:
                defaultMCSAlgorithm_1();
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
        CDKMCSHandler mcs = null;
        mcs = new CDKMCSHandler();
        mcs.set(query, target);
        mcs.searchMCS(isMatchBonds());

        clearMaps();

        firstSolution.putAll(mcs.getFirstMapping());
        allMCS.addAll(mcs.getAllMapping());

        firstAtomMCS.putAll(mcs.getFirstAtomMapping());
        allAtomMCS.addAll(mcs.getAllAtomMapping());

    }

    private synchronized void mcsPlusAlgorithm() {
//        double time = System.currentTimeMillis();
        MCSPlusHandler mcs = null;
        mcs = new MCSPlusHandler();

        mcs.set(query, target);
        mcs.searchMCS(isMatchBonds());

        clearMaps();

        firstSolution.putAll(mcs.getFirstMapping());
        allMCS.addAll(mcs.getAllMapping());

        firstAtomMCS.putAll(mcs.getFirstAtomMapping());
        allAtomMCS.addAll(mcs.getAllAtomMapping());
//        System.out.println("\nMCSPlus used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));

    }

    private void vfLibMCS() {
        VFlibMCSHandler mcs = null;
        mcs = new VFlibMCSHandler();
        mcs.set(query, target);
        mcs.searchMCS(isMatchBonds());

        clearMaps();
        firstSolution.putAll(mcs.getFirstMapping());
        allMCS.addAll(mcs.getAllMapping());

        firstAtomMCS.putAll(mcs.getFirstAtomMapping());
        allAtomMCS.addAll(mcs.getAllAtomMapping());
    }

    private void singleMapping() {
        SingleMappingHandler mcs = null;

        mcs = new SingleMappingHandler();
        mcs.set(query, target);
        mcs.searchMCS(isMatchBonds());

        clearMaps();
        firstSolution.putAll(mcs.getFirstMapping());
        allMCS.addAll(mcs.getAllMapping());

        firstAtomMCS.putAll(mcs.getFirstAtomMapping());
        allAtomMCS.addAll(mcs.getAllAtomMapping());

    }

    private int getHCount(IAtomContainer molecule) {
        int count = 0;
        for (IAtom atom : molecule.atoms()) {
            if (atom.getSymbol().equalsIgnoreCase("H")) {
                ++count;
            }
        }
        return count;
    }

    private void defaultMCSAlgorithm() {
        try {
            if (isMatchBonds()) {
                cdkMCSAlgorithm();
                if (getFirstMapping() == null || isTimeOut()) {
                    vfLibMCS();
                }
            } else {
                try {
                    mcsPlusAlgorithm();
                } catch (Exception e) {
                    cdkMCSAlgorithm();
                }
                if (getFirstMapping() == null || isTimeOut()) {
                    vfLibMCS();
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void defaultMCSAlgorithm_1() {
        try {
            cdkMCSAlgorithm();
            if (getFirstMapping() == null || isTimeOut()) {
//                System.out.println("\nCDKMCS hit by timeout\n");
//                double time = System.currentTimeMillis();
                vfLibMCS();
//                System.out.println("\nVF Lib used\n" + ((System.currentTimeMillis() - time) / (60 * 1000)));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void vfLibMCSAlgorithm() {
        vfLibMCS();
    }

    private void setTime(boolean bondTypeFlag) {
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

    public boolean isTimeOut() {
        return TimeOut.getInstance().isTimeOutFlag();
    }

    public void resetTimeOut() {
        TimeOut.getInstance().setTimeOutFlag(false);
    }

    private void clearMaps() {
        this.firstSolution.clear();
        this.allMCS.clear();
        this.allAtomMCS.clear();
        this.firstAtomMCS.clear();
    }

    /**
     *
     * @param query 
     *
     */
    @Override
    public void init(IQueryAtomContainer query, IAtomContainer target) throws CDKException {
        this.query = query;
        this.target = target;
        mcsBuilder(query, target);
    }

    @Override
    public void init(IAtomContainer query, IAtomContainer target) throws CDKException {
        this.query = query;
        this.target = target;
        mcsBuilder(query, target);
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testSetChemFilters")
    public void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (firstAtomMCS != null) {
            ChemicalFilters chemFilter = new ChemicalFilters(allMCS, allAtomMCS, firstSolution, firstAtomMCS, getQueryMolecule(), getTargetMolecule());

            if (energyFilter) {
                try {
                    chemFilter.sortResultsByEnergies();
                    this.bEnergies = chemFilter.getSortedEnergy();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }

            if (stereoFilter && firstAtomMCS.size() > 1) {
                try {
                    chemFilter.sortResultsByStereoAndBondMatch();
                    this.stereoScore = chemFilter.getStereoMatches();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }

            if (fragmentFilter) {
                chemFilter.sortResultsByFragments();
                this.fragmentSize = chemFilter.getSortedFragment();
            }
        }
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetFragmentSize")
    public synchronized Integer getFragmentSize(int Key) {
        return (fragmentSize != null && !fragmentSize.isEmpty())
                ? fragmentSize.get(Key) : null;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetStereoScore")
    public synchronized Integer getStereoScore(int Key) {
        return (stereoScore != null && !stereoScore.isEmpty()) ? stereoScore.get(Key).intValue() : null;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetEnergyScore")
    public synchronized Double getEnergyScore(int Key) {
        return (bEnergies != null && !bEnergies.isEmpty()) ? bEnergies.get(Key) : null;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetFirstMapping")
    public synchronized Map<Integer, Integer> getFirstMapping() {
        return firstSolution.isEmpty() ? null : firstSolution;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetAllMapping")
    public synchronized List<Map<Integer, Integer>> getAllMapping() {
        return allMCS.isEmpty() ? null : allMCS;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetFirstAtomMapping")
    public synchronized Map<IAtom, IAtom> getFirstAtomMapping() {
        return firstAtomMCS.isEmpty() ? null : firstAtomMCS;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetAllAtomMapping")
    public synchronized List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS.isEmpty() ? null : allAtomMCS;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetreactantMolecule")
    public IAtomContainer getQueryMolecule() {
        return query;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetproductMolecule")
    public IAtomContainer getTargetMolecule() {
        return target;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetTanimotoSimilarity")
    public double getTanimotoSimilarity() throws IOException {
        double tanimoto = getTanimotoAtomSimilarity() + getTanimotoBondSimilarity();
        if (tanimoto > 0 && getQueryMolecule().getBondCount() > 0
                && getTargetMolecule().getBondCount() > 0) {
            tanimoto /= 2;
        }
        return tanimoto;
    }

    public double getTanimotoAtomSimilarity() throws IOException {
        int decimalPlaces = 4;
        int rAtomCount = 0;
        int pAtomCount = 0;
        double tanimotoAtom = 0.0;

        if (getFirstMapping() != null && !getFirstMapping().isEmpty()) {

            rAtomCount = getQueryMolecule().getAtomCount();
            pAtomCount = getTargetMolecule().getAtomCount();

            double matchCount = getFirstMapping().size();
            tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
            BigDecimal tan = new BigDecimal(tanimotoAtom);
            tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
            tanimotoAtom = tan.doubleValue();
        }
        return tanimotoAtom;
    }

    public double getTanimotoBondSimilarity() throws IOException {
        int decimalPlaces = 4;
        int rBondCount = 0;
        int pBondCount = 0;
        double tanimotoAtom = 0.0;

        if (getFirstBondMap() != null && !getFirstBondMap().isEmpty()) {
            rBondCount = getQueryMolecule().getBondCount();
            pBondCount = getTargetMolecule().getBondCount();

            double matchCount = getFirstBondMap().size();
            tanimotoAtom = (matchCount) / (rBondCount + pBondCount - matchCount);
            BigDecimal tan = new BigDecimal(tanimotoAtom);
            tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
            tanimotoAtom = tan.doubleValue();
        }
        return tanimotoAtom;
    }

    /** {@inheritDoc}
     *
     */
    @Override
    @TestMethod("testIsStereoMisMatch")
    public boolean isStereoMisMatch() {
        boolean flag = false;
        IAtomContainer reactant = getQueryMolecule();
        IAtomContainer product = getTargetMolecule();
        int Score = 0;

        for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.entrySet()) {
            IAtom indexI = mappingI.getKey();
            IAtom indexJ = mappingI.getValue();
            for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.entrySet()) {

                IAtom indexIPlus = mappingJ.getKey();
                IAtom indexJPlus = mappingJ.getValue();
                if (!indexI.equals(indexIPlus) && !indexJ.equals(indexJPlus)) {

                    IAtom sourceAtom1 = indexI;
                    IAtom sourceAtom2 = indexIPlus;

                    IBond rBond = reactant.getBond(sourceAtom1, sourceAtom2);

                    IAtom targetAtom1 = indexJ;
                    IAtom targetAtom2 = indexJPlus;
                    IBond pBond = product.getBond(targetAtom1, targetAtom2);

                    if ((rBond != null && pBond != null) && (rBond.getStereo() != pBond.getStereo())) {
                        Score++;
                    }
                }
            }
        }
        if (Score > 0) {
            flag = true;
        }
        return flag;
    }

    /** {@inheritDoc}
     *
     */
    @Override
    @TestMethod("testIsSubgraph")
    public boolean isSubgraph() {

        IAtomContainer reactant = getQueryMolecule();
        IAtomContainer product = getTargetMolecule();

        float mappingSize = 0;
        if (firstSolution != null && !firstSolution.isEmpty()) {
            mappingSize = firstSolution.size();
        } else {
            return false;
        }
        int sourceAtomCount = reactant.getAtomCount();
        int targetAtomCount = product.getAtomCount();

        if (mappingSize == sourceAtomCount && mappingSize <= targetAtomCount) {
            if (!getFirstBondMap().isEmpty()
                    && getFirstBondMap().size() == reactant.getBondCount()) {
                return true;
            } else if (mappingSize == 1) {
                return true;
            }
        }
        return false;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetEuclideanDistance")
    public double getEuclideanDistance() throws IOException {
        int decimalPlaces = 4;
        double sourceAtomCount = 0;
        double targetAtomCount = 0;
        double euclidean = -1;


        if (getFirstMapping() != null || !getFirstMapping().isEmpty()) {

            sourceAtomCount = getQueryMolecule().getAtomCount();
            targetAtomCount = getTargetMolecule().getAtomCount();

            double common = getFirstMapping().size();
            euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
            BigDecimal dist = new BigDecimal(euclidean);
            dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
            euclidean = dist.doubleValue();
        }
        return euclidean;
    }

    /**
     * @return the matchBonds
     */
    public boolean isMatchBonds() {
        return matchBonds;
    }

    /**
     * @param matchBonds the matchBonds to set
     */
    public void setMatchBonds(boolean matchBonds) {
        this.matchBonds = matchBonds;
    }

    /**
     * @return the allBondMCS
     */
    public List<Map<IBond, IBond>> getAllBondMaps() {
        return allBondMCS;
    }

    /**
     * @param allBondMCS the allBondMCS to set
     */
    private void setAllBondMaps(List<Map<IBond, IBond>> allBondMCS) {
        this.allBondMCS = allBondMCS;
    }

    /**
     * @return the firstBondMCS
     */
    public Map<IBond, IBond> getFirstBondMap() {
        return firstBondMCS;
    }

    /**
     * @param firstBondMCS the firstBondMCS to set
     */
    private void setFirstBondMap(Map<IBond, IBond> firstBondMCS) {
        this.firstBondMCS = firstBondMCS;
    }

    @Override
    public double getBondSensitiveCDKMCSTimeOut() {
        return this.bondSensitiveCDKMCSTimeOut;
    }

    @Override
    public void setBondSensitiveCDKMCSTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveCDKMCSTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public double getBondInSensitiveCDKMCSTimeOut() {
        return bondInSensitiveCDKMCSTimeOut;
    }

    @Override
    public void setBondInSensitiveCDKMCSTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveCDKMCSTimeOut = bondInSensitiveTimeOut;
    }

    @Override
    public double getBondSensitiveMCSPlusTimeOut() {
        return this.bondSensitiveMCSPlusTimeOut;
    }

    @Override
    public void setBondSensitiveMCSPlusTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveMCSPlusTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public double getBondInSensitiveMCSPlusTimeOut() {
        return this.bondInSensitiveMCSPlusTimeOut;
    }

    @Override
    public void setBondInSensitiveMCSPlusTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveMCSPlusTimeOut = bondInSensitiveTimeOut;
    }

    @Override
    public double getBondSensitiveVFTimeOut() {
        return this.bondSensitiveVFTimeOut;
    }

    @Override
    public void setBondSensitiveVFTimeOut(double bondSensitiveTimeOut) {
        this.bondSensitiveVFTimeOut = bondSensitiveTimeOut;
    }

    @Override
    public double getBondInSensitiveVFTimeOut() {
        return this.bondInSensitiveVFTimeOut;
    }

    @Override
    public void setBondInSensitiveVFTimeOut(double bondInSensitiveTimeOut) {
        this.bondInSensitiveVFTimeOut = bondInSensitiveTimeOut;
    }
}

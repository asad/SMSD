/* Copyright (C) 2009-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.algorithm.single.SingleMappingHandler;
import org.openscience.smsd.algorithm.vflib.substructure.AtomMapping;
import org.openscience.smsd.algorithm.vflib.substructure.VF2;
import org.openscience.smsd.filters.ChemicalFilters;
import org.openscience.smsd.global.TimeOut;

/**
 * This is an ultra fast method to report if query
 * is a substructure for target molecule. If this case is true
 * then it returns only all mapping.
 *
 * This is much faster than {@link
 * org.openscience.cdk.smsd.algorithm.vflib.VFlibHandler} class
 * as it only reports first match and backtracks.
 *
 * This class should only be used to report if a query
 * graph is a substructure of the target graph.
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Substructure implements IAtomAtomMapping {

    private List<Double> stereoScore = null;
    private List<Integer> fragmentSize = null;
    private List<Double> bEnergies = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private Map<Integer, Integer> firstIndexMCS = null;
    private List<Map<Integer, Integer>> allIndexMCS = null;
    private IAtomContainer mol1 = null;
    private IAtomContainer mol2 = null;
    private List<AtomMapping> vfLibSolutions = null;
    private int vfMappingSize = -1;
    private boolean bond_Match_Flag = false;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(Substructure.class);

    /**
     * Constructor for VF Substructure Algorithm 
     */
    public Substructure() {
        allAtomMCS = new ArrayList<Map<IAtom, IAtom>>();
        firstAtomMCS = new HashMap<IAtom, IAtom>();
        firstIndexMCS = new TreeMap<Integer, Integer>();
        allIndexMCS = new ArrayList<Map<Integer, Integer>>();

        TimeOut tmo = TimeOut.getInstance();
        tmo.setCDKMCSTimeOut(0.15);
    }

    private void setFirstMappings() {
        if (!allAtomMCS.isEmpty()) {
            firstAtomMCS.putAll(allAtomMCS.get(0));
            firstIndexMCS.putAll(allIndexMCS.get(0));
        }
    }

    /** {@inheritDoc}
     *
     * Set the VFLib MCS software
     *
     * @param reactant
     * @param product
     */
    @Override
    public void init(IAtomContainer reactant, IAtomContainer product) {
        this.mol1 = reactant;
        this.mol2 = product;
    }

    /** {@inheritDoc}
     *
     * Set the VFLib MCS software
     *
     * @param reactant
     * @param product
     */
    @Override
    public void init(IQueryAtomContainer reactant, IAtomContainer product) {
        this.mol1 = reactant;
        this.mol2 = product;
    }

    private boolean hasMap(Map<Integer, Integer> map, List<Map<Integer, Integer>> mapGlobal) {
        for (Map<Integer, Integer> test : mapGlobal) {
            if (test.equals(map)) {
                return true;
            }
        }
        return false;
    }

    /** {@inheritDoc}
     *
     * @return 
     */
    @Override
    public List<Map<IAtom, IAtom>> getAllAtomMapping() {
        return allAtomMCS;
    }

    /** {@inheritDoc}
     * @return 
     */
    @Override
    public List<Map<Integer, Integer>> getAllMapping() {
        return allIndexMCS;
    }

    /** {@inheritDoc}
     * @return 
     */
    @Override
    public Map<IAtom, IAtom> getFirstAtomMapping() {
        return firstAtomMCS;
    }

    /** {@inheritDoc}
     * @return 
     */
    @Override
    public Map<Integer, Integer> getFirstMapping() {
        return firstIndexMCS;
    }

    private void setVFMappings() {
        int counter = 0;
        for (AtomMapping solution : vfLibSolutions) {
            Map<IAtom, IAtom> atomatomMapping = new HashMap<IAtom, IAtom>();
            Map<Integer, Integer> indexindexMapping = new TreeMap<Integer, Integer>();
            if (solution.getSize() > vfMappingSize) {
                this.vfMappingSize = solution.getSize();
                counter = 0;
            }
            for (Map.Entry<IAtom, IAtom> mapping : solution.getAtomMapping().entrySet()) {
                IAtom qAtom = null;
                IAtom tAtom = null;

                qAtom = mapping.getKey();
                tAtom = mapping.getValue();

                Integer qIndex = Integer.valueOf(getReactantMol().getAtomNumber(qAtom));
                Integer tIndex = Integer.valueOf(getProductMol().getAtomNumber(tAtom));
                if (qIndex != -1 && tIndex != -1) {
                    atomatomMapping.put(qAtom, tAtom);
                    indexindexMapping.put(qIndex, tIndex);
                } else {
                    try {
                        throw new CDKException("Atom index pointing to NULL");
                    } catch (CDKException ex) {
                        Logger.error(Level.SEVERE, null, ex);
                    }
                }
            }
            if (!atomatomMapping.isEmpty() && !hasMap(indexindexMapping, allIndexMCS)
                    && indexindexMapping.size() == vfMappingSize) {
                allAtomMCS.add(counter, atomatomMapping);
                allIndexMCS.add(counter, indexindexMapping);
                counter++;
            }
        }
    }

    public boolean isSubgraph(boolean shouldMatchBonds) throws CDKException {

        setBondMatchFlag(shouldMatchBonds);
        if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
            return false;
        } else {
            if (getReactantMol() != null && getProductMol() != null) {

                VF2 mapper = new VF2();
                if (!mapper.isomorphism(getReactantMol(), getProductMol(), shouldMatchBonds).isEmpty()) {
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }

    /**
     * 
     * @param shouldMatchBonds
     * @return
     */
    public boolean findSubgraphs(boolean shouldMatchBonds) {

//        System.out.println("Mol1 Size: -> " + getReactantMol().getAtomCount());
//        System.out.println("Mol2 Size. -> " + getProductMol().getAtomCount());
        setBondMatchFlag(shouldMatchBonds);


        if (getReactantMol().getAtomCount() == 1 || getProductMol().getAtomCount() == 1) {
            singleMapping(shouldMatchBonds);
        } else {
            if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
                return false;
            } else {

                vfLibSolutions = new ArrayList<AtomMapping>();
                if (getProductMol() != null) {
                    VF2 mapper = new VF2();
                    List<AtomMapping> atomMappings = mapper.isomorphisms(getReactantMol(), getProductMol(), shouldMatchBonds);
                    if (!atomMappings.isEmpty()) {
                        vfLibSolutions.addAll(atomMappings);
                    } else {
                        return false;
                    }
                }
                setVFMappings();
            }
            if (!allAtomMCS.isEmpty()) {
                setFirstMappings();
            }
        }
        return (!allIndexMCS.isEmpty() && allIndexMCS.iterator().next().size() == getReactantMol().getAtomCount()) ? true : false;
    }

    /**
     * 
     * @param shouldMatchBonds
     * @return
     */
    public boolean findSubgraph(boolean shouldMatchBonds) {

        setBondMatchFlag(shouldMatchBonds);
        if (getReactantMol().getAtomCount() == 1 || getProductMol().getAtomCount() == 1) {
            singleMapping(shouldMatchBonds);
        } else {
            if (getReactantMol().getAtomCount() > getProductMol().getAtomCount()) {
                return false;
            } else {
                vfLibSolutions = new ArrayList<AtomMapping>();
                if (getProductMol() != null) {
                    VF2 mapper = new VF2();
                    AtomMapping atomMapping = mapper.isomorphism(getReactantMol(), getProductMol(), shouldMatchBonds);
                    if (!atomMapping.isEmpty()) {
                        vfLibSolutions.add(atomMapping);
                    } else {
                        return false;
                    }
                }
                setVFMappings();
            }
            if (!allAtomMCS.isEmpty()) {
                setFirstMappings();
            }
        }
        return (!allIndexMCS.isEmpty() && allIndexMCS.iterator().next().size() == getReactantMol().getAtomCount()) ? true : false;
    }

    private void singleMapping(boolean shouldMatchBonds) {
        SingleMappingHandler mcs = null;

        mcs = new SingleMappingHandler();
        mcs.set(getReactantMol(), getProductMol());
        mcs.searchMCS(shouldMatchBonds);

        firstIndexMCS.putAll(mcs.getFirstMapping());
        allIndexMCS.addAll(mcs.getAllMapping());

        firstAtomMCS.putAll(mcs.getFirstAtomMapping());
        allAtomMCS.addAll(mcs.getAllAtomMapping());

    }

    /**
     * @return the shouldMatchBonds
     */
    public boolean isBondMatchFlag() {
        return bond_Match_Flag;
    }

    /**
     * @param shouldMatchBonds the shouldMatchBonds to set
     */
    public void setBondMatchFlag(boolean shouldMatchBonds) {
        this.bond_Match_Flag = shouldMatchBonds;
    }

    private IAtomContainer getReactantMol() {
        return mol1;
    }

    private IAtomContainer getProductMol() {
        return mol2;
    }

    @Override
    public void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (firstAtomMCS != null) {
            ChemicalFilters chemFilter = new ChemicalFilters(allIndexMCS, allAtomMCS,
                    getFirstMapping(), getFirstAtomMapping(), getReactantMol(), getProductMol());

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

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        return (fragmentSize != null && !fragmentSize.isEmpty())
                ? fragmentSize.get(Key) : null;
    }

    @Override
    public synchronized Integer getStereoScore(int Key) {
        return (stereoScore != null && !stereoScore.isEmpty()) ? stereoScore.get(Key).intValue() : null;
    }

    @Override
    public synchronized Double getEnergyScore(int Key) {
        return (bEnergies != null && !bEnergies.isEmpty()) ? bEnergies.get(Key) : null;
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

            sourceAtomCount = this.mol1.getAtomCount();
            targetAtomCount = this.mol2.getAtomCount();

            double common = getFirstMapping().size();
            euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
            BigDecimal dist = new BigDecimal(euclidean);
            dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
            euclidean = dist.doubleValue();
        }
        return euclidean;
    }

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetTanimotoSimilarity")
    public double getTanimotoSimilarity() throws IOException {
        int decimalPlaces = 4;
        int rAtomCount = 0;
        int pAtomCount = 0;
        double tanimotoAtom = 0.0;

        if (getFirstMapping() != null && !getFirstMapping().isEmpty()) {

            rAtomCount = this.mol1.getAtomCount();
            pAtomCount = this.mol2.getAtomCount();

            double matchCount = getFirstMapping().size();
            tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
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
        IAtomContainer reactant = this.mol1;
        IAtomContainer product = this.mol2;
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

    @Override
    public IAtomContainer getTargetMolecule() {
        return this.mol1;
    }

    @Override
    public IAtomContainer getQueryMolecule() {
        return this.mol2;
    }
}

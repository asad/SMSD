/* Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
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

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.smsd.filters.ChemicalFilters;
import org.openscience.smsd.interfaces.IAtomMapping;

/**
 *
 * @cdk.require java1.5+
 *
 * @cdk.module smsd @cdk.githash
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
public class BaseMapping implements IAtomMapping {

    protected List<AtomAtomMapping> mcsList;
    protected boolean matchBonds;
    protected boolean matchRings;
    protected boolean subgraph;
    protected IAtomContainer mol1;
    protected IAtomContainer mol2;
    private List<Double> stereoScoreList;
    private List<Integer> fragmentSizeList;
    private List<Double> bondEnergiesList;
    private final static ILoggingTool Logger =
            LoggingToolFactory.createLoggingTool(BaseMapping.class);

    @Override
    public synchronized void setChemFilters(boolean stereoFilter, boolean fragmentFilter, boolean energyFilter) {

        if (getMappingCount() > 0) {
            ChemicalFilters chemFilter = new ChemicalFilters(mcsList, mol1, mol2);

            if (energyFilter) {
                try {
                    chemFilter.sortResultsByEnergies();
                    this.bondEnergiesList = chemFilter.getSortedEnergy();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }

            if (stereoFilter) {
                try {
                    chemFilter.sortResultsByStereoAndBondMatch();
                    this.stereoScoreList = chemFilter.getStereoMatches();
                } catch (CDKException ex) {
                    Logger.error(Level.SEVERE, null, ex);
                }
            }

            if (fragmentFilter) {
                chemFilter.sortResultsByFragments();
                this.fragmentSizeList = chemFilter.getSortedFragment();
            }
        }
    }

    @Override
    public synchronized Integer getFragmentSize(int Key) {
        return (fragmentSizeList != null && !fragmentSizeList.isEmpty())
                ? fragmentSizeList.get(Key) : null;
    }

    @Override
    public synchronized Integer getStereoScore(int Key) {
        return (stereoScoreList != null && !stereoScoreList.isEmpty()) ? stereoScoreList.get(Key).intValue() : null;
    }

    @Override
    public synchronized Double getEnergyScore(int Key) {
        return (bondEnergiesList != null && !bondEnergiesList.isEmpty()) ? bondEnergiesList.get(Key) : null;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    @TestMethod("testGetTanimotoSimilarity")
    public synchronized double getTanimotoSimilarity() {
        int decimalPlaces = 4;
        double rAtomCount;
        double pAtomCount;
        double tanimotoAtom = 0.0;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = mcsList.iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                rAtomCount = (double) this.mcsList.iterator().next().getQuery().getAtomCount();
                pAtomCount = (double) this.mcsList.iterator().next().getTarget().getAtomCount();

                double matchCount = (double) firstAtomMCS.getCount();
                tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
                BigDecimal tan = new BigDecimal(tanimotoAtom);
                tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                tanimotoAtom = tan.doubleValue();
            }
        }
        return tanimotoAtom;
    }

    /**
     * {@inheritDoc}
     *
     */
    @Override
    @TestMethod("testIsStereoMisMatch")
    public synchronized boolean isStereoMisMatch() {
        boolean flag = false;
        IAtomContainer reactant = this.mol1;
        IAtomContainer product = this.mol2;
        int stereoMisMatchScore = 0;
        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = mcsList.iterator().next();
            for (Map.Entry<IAtom, IAtom> mappingI : firstAtomMCS.getMappings().entrySet()) {
                IAtom indexI = mappingI.getKey();
                IAtom indexJ = mappingI.getValue();
                for (Map.Entry<IAtom, IAtom> mappingJ : firstAtomMCS.getMappings().entrySet()) {

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
                            stereoMisMatchScore++;
                        }
                    }
                }
            }
        }
        if (stereoMisMatchScore > 0) {
            flag = true;
        }
        return flag;
    }

    @Override
    public synchronized int getMappingCount() {
        return this.mcsList.isEmpty() ? 0 : this.mcsList.size();
    }

    /**
     * {@inheritDoc}
     */
    @TestMethod("testGetEuclideanDistance")
    @Override
    public synchronized double getEuclideanDistance() {
        int decimalPlaces = 4;
        double sourceAtomCount;
        double targetAtomCount;
        double euclidean = -1.;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = mcsList.iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                sourceAtomCount = (double) this.mcsList.iterator().next().getQuery().getAtomCount();
                targetAtomCount = (double) this.mcsList.iterator().next().getTarget().getAtomCount();

                double common = (double) firstAtomMCS.getCount();
                euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
                BigDecimal dist = new BigDecimal(euclidean);
                dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                euclidean = dist.doubleValue();
            }
        }
        return euclidean;
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(new ArrayList<AtomAtomMapping>(mcsList));
    }

    /**
     * {@inheritDoc}
     *
     * @return
     */
    @Override
    public synchronized AtomAtomMapping getFirstAtomMapping() {
        return mcsList.isEmpty() ? new AtomAtomMapping(mol1, mol2)
                : mcsList.iterator().next();
    }

    @Override
    public synchronized IAtomContainer getQueryContainer() {
        return this.mol1;
    }

    @Override
    public synchronized IAtomContainer getTargetContainer() {
        return this.mol2;
    }

    /**
     * Returns true if bond are to be matched.
     *
     * @return true if bond are to be matched
     */
    protected boolean isMatchBonds() {
        return matchBonds;
    }

    /**
     * Returns true if rings are to be matched.
     *
     * @return true if rings are to be matched
     */
    protected boolean isMatchRings() {
        return matchRings;
    }

    /**
     * Returns true if Query is a subgraph of the Target.
     *
     * @return true if Query is a subgraph of the Target
     */
    public synchronized boolean isSubgraph() {
        return this.subgraph;
    }
}

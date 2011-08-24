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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.filters.ChemicalFilters;
import org.openscience.smsd.ring.HanserRingFinder;

/**
 *
 * @cdk.require java1.5+
 *
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
public class BaseMapping implements IAtomMapping {

    protected List<AtomAtomMapping> mcsList;
    protected IAtomContainer mol1;
    protected IAtomContainer mol2;
    protected IQueryAtomContainer queryMol;
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

    /** {@inheritDoc}
     */
    @Override
    @TestMethod("testGetTanimotoSimilarity")
    public synchronized double getTanimotoSimilarity() {
        int decimalPlaces = 4;
        int rAtomCount = 0;
        int pAtomCount = 0;
        double tanimotoAtom = 0.0;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = mcsList.iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                rAtomCount = this.mol1.getAtomCount();
                pAtomCount = this.mol2.getAtomCount();

                double matchCount = firstAtomMCS.getCount();
                tanimotoAtom = (matchCount) / (rAtomCount + pAtomCount - matchCount);
                BigDecimal tan = new BigDecimal(tanimotoAtom);
                tan = tan.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                tanimotoAtom = tan.doubleValue();
            }
        }
        return tanimotoAtom;
    }

    /** {@inheritDoc}
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

    /** {@inheritDoc}
     */
    @TestMethod("testGetEuclideanDistance")
    @Override
    public synchronized double getEuclideanDistance() {
        int decimalPlaces = 4;
        double sourceAtomCount = 0;
        double targetAtomCount = 0;
        double euclidean = -1.;

        if (getMappingCount() > 0) {
            AtomAtomMapping firstAtomMCS = mcsList.iterator().next();

            if (!firstAtomMCS.isEmpty()) {

                sourceAtomCount = this.mol1.getAtomCount();
                targetAtomCount = this.mol2.getAtomCount();

                double common = firstAtomMCS.getCount();
                euclidean = Math.sqrt(sourceAtomCount + targetAtomCount - 2 * common);
                BigDecimal dist = new BigDecimal(euclidean);
                dist = dist.setScale(decimalPlaces, BigDecimal.ROUND_HALF_UP);
                euclidean = dist.doubleValue();
            }
        }
        return euclidean;
    }

    /** {@inheritDoc}
     *
     * @return 
     */
    @Override
    public synchronized List<AtomAtomMapping> getAllAtomMapping() {
        return Collections.unmodifiableList(new ArrayList<AtomAtomMapping>(mcsList));
    }

    /** {@inheritDoc}
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
     * Prepare the target molecule for analysis. 
     * <p/> We perform ring perception and aromaticity detection and set up
     * the appropriate properties. Right now, this function is called each time 
     * we need to do a query and this is
     * inefficient.
     *
     * @param atomContainer Atom container where rings are to be marked
     * @throws CDKException if there is a problem in ring perception or aromaticity detection, 
     * which is usually related to a timeout in the ring finding code.
     */
    protected void initializeMolecule(IAtomContainer atomContainer) throws CDKException {
        // Code copied from
        // org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
        Map<String, Integer> valencesTable = new HashMap<String, Integer>();
        valencesTable.put("H", 1);
        valencesTable.put("Li", 1);
        valencesTable.put("Be", 2);
        valencesTable.put("B", 3);
        valencesTable.put("C", 4);
        valencesTable.put("N", 5);
        valencesTable.put("O", 6);
        valencesTable.put("F", 7);
        valencesTable.put("Na", 1);
        valencesTable.put("Mg", 2);
        valencesTable.put("Al", 3);
        valencesTable.put("Si", 4);
        valencesTable.put("P", 5);
        valencesTable.put("S", 6);
        valencesTable.put("Cl", 7);
        valencesTable.put("K", 1);
        valencesTable.put("Ca", 2);
        valencesTable.put("Ga", 3);
        valencesTable.put("Ge", 4);
        valencesTable.put("As", 5);
        valencesTable.put("Se", 6);
        valencesTable.put("Br", 7);
        valencesTable.put("Rb", 1);
        valencesTable.put("Sr", 2);
        valencesTable.put("In", 3);
        valencesTable.put("Sn", 4);
        valencesTable.put("Sb", 5);
        valencesTable.put("Te", 6);
        valencesTable.put("I", 7);
        valencesTable.put("Cs", 1);
        valencesTable.put("Ba", 2);
        valencesTable.put("Tl", 3);
        valencesTable.put("Pb", 4);
        valencesTable.put("Bi", 5);
        valencesTable.put("Po", 6);
        valencesTable.put("At", 7);
        valencesTable.put("Fr", 1);
        valencesTable.put("Ra", 2);
        valencesTable.put("Cu", 2);
        valencesTable.put("Mn", 2);
        valencesTable.put("Co", 2);

        // do all ring perception
//        AllRingsFinder arf = new AllRingsFinder();
//        IRingSet allRings;
//        try {
//            allRings = arf.findAllRings(atomContainer);
//        } catch (CDKException e) {
//            Logger.error(Level.SEVERE, null, e);
//            throw new CDKException(e.toString(), e);
//        }

        // sets Hanser Ring information
        IRingSet allRings = HanserRingFinder.getRingSet(atomContainer);

        // sets SSSR information
        SSSRFinder finder = new SSSRFinder(atomContainer);
        IRingSet sssr = finder.findEssentialRings();

        for (IAtom atom : atomContainer.atoms()) {

            // add a property to each ring atom that will be an array of
            // Integers, indicating what size ring the given atom belongs to
            // Add SSSR ring counts
            if (allRings.contains(atom)) { // it's in a ring
                atom.setFlag(CDKConstants.ISINRING, true);
                // lets find which ring sets it is a part of
                List<Integer> ringsizes = new ArrayList<Integer>();
                IRingSet currentRings = allRings.getRings(atom);
                int min = 0;
                for (int i = 0; i < currentRings.getAtomContainerCount(); i++) {
                    int size = currentRings.getAtomContainer(i).getAtomCount();
                    if (min > size) {
                        min = size;
                    }
                    ringsizes.add(size);
                }
                atom.setProperty(CDKConstants.RING_SIZES, ringsizes);
                atom.setProperty(CDKConstants.SMALLEST_RINGS, sssr.getRings(atom));
            } else {
                atom.setFlag(CDKConstants.ISINRING, false);
            }

            // determine how many rings bonds each atom is a part of
            int hCount;
            if (atom.getImplicitHydrogenCount() == CDKConstants.UNSET) {
                hCount = 0;
            } else {
                hCount = atom.getImplicitHydrogenCount();
            }

            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);
            int total = hCount + connectedAtoms.size();
            for (IAtom connectedAtom : connectedAtoms) {
                if (connectedAtom.getSymbol().equals("H")) {
                    hCount++;
                }
            }
            atom.setProperty(CDKConstants.TOTAL_CONNECTIONS, total);
            atom.setProperty(CDKConstants.TOTAL_H_COUNT, hCount);

            if (valencesTable.get(atom.getSymbol()) != null) {
                int formalCharge = atom.getFormalCharge() == CDKConstants.UNSET ? 0 : atom.getFormalCharge();
                atom.setValency(valencesTable.get(atom.getSymbol()) - formalCharge);
            }
        }

        for (IBond bond : atomContainer.bonds()) {
            if (allRings.getRings(bond).getAtomContainerCount() > 0) {
                bond.setFlag(CDKConstants.ISINRING, true);
            }
        }

        for (IAtom atom : atomContainer.atoms()) {
            List<IAtom> connectedAtoms = atomContainer.getConnectedAtomsList(atom);

            int counter = 0;
            IAtom any;
            for (IAtom connectedAtom : connectedAtoms) {
                any = connectedAtom;
                if (any.getFlag(CDKConstants.ISINRING)) {
                    counter++;
                }
            }
            atom.setProperty(CDKConstants.RING_CONNECTIONS, counter);
        }

        // check for atomaticity
        try {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atomContainer);
            CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        } catch (CDKException e) {
            Logger.error(Level.SEVERE, null, e);
            throw new CDKException(e.toString(), e);
        }
    }
}

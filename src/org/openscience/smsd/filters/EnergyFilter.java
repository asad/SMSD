/* Copyright (C) 2009-2010  Syed Asad Rahman <asad@ebi.ac.uk>
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
package org.openscience.smsd.filters;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.smsd.tools.BondEnergies;

/**
 * Filter based on energies.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * @cdk.module smsd
 */
public class EnergyFilter extends BaseFilter implements IChemicalFilter<Double> {

//    public static final Double MAX_ENERGY = Double.MAX_VALUE;
    public static final Double MAX_ENERGY = 99999999.99;
    private List<Double> bEnergies = null;

    public EnergyFilter(IAtomContainer rMol, IAtomContainer pMol) {
        super(rMol, pMol);
        bEnergies = new ArrayList<Double>();
    }

    @Override
    public Double sortResults(
            Map<Integer, Map<Integer, Integer>> allEnergyMCS,
            Map<Integer, Map<IAtom, IAtom>> allAtomEnergyMCS,
            Map<Integer, Double> energySelectionMap) throws CDKException {

        for (Integer Key : allEnergyMCS.keySet()) {
            Map<Integer, Integer> mcsAtom = allEnergyMCS.get(Key);
            Double Energies = getMappedMoleculeEnergies(mcsAtom);
            energySelectionMap.put(Key, Energies);
        }

        energySelectionMap = sortMapByValueInAscendingOrder(energySelectionMap);

        double lowestEnergyScore = MAX_ENERGY;
        for (Integer key : energySelectionMap.keySet()) {
            lowestEnergyScore = energySelectionMap.get(key);
            break;
        }
        return lowestEnergyScore;
    }

    @Override
    public List<Double> getScores() {
        return Collections.unmodifiableList(bEnergies);
    }

    @Override
    public void clearScores() {
        bEnergies.clear();
    }

    @Override
    public void addScore(int counter, Double value) {
        bEnergies.add(counter, value);
    }

    @Override
    public void fillMap(Map<Integer, Double> energySelectionMap) {
        int Index = 0;
        for (Double score : bEnergies) {
            energySelectionMap.put(Index, score);
            Index++;
        }
    }

    private synchronized Double getMappedMoleculeEnergies(Map<Integer, Integer> MCSAtomSolution) throws CDKException {

//      System.out.println("\nSort By Energies");
        double totalBondEnergy = -9999.0;

        IAtomContainer educt = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, rMol);
        IAtomContainer product = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, pMol);

        for (IAtom eAtom : educt.atoms()) {
            eAtom.setFlag(0, false);
        }

        for (IAtom pAtom : product.atoms()) {
            pAtom.setFlag(0, false);
        }

        if (MCSAtomSolution != null) {
            for (Map.Entry<Integer, Integer> map : MCSAtomSolution.entrySet()) {
                int eNum = map.getKey();
                int pNum = map.getValue();

                IAtom eAtom = educt.getAtom(eNum);
                IAtom pAtom = product.getAtom(pNum);

                eAtom.setFlag(0, true);
                pAtom.setFlag(0, true);
            }
        }

        if (MCSAtomSolution != null) {
            totalBondEnergy = getEnergy(educt, product);
        }
        return totalBondEnergy;
    }

    private double getEnergy(IAtomContainer Educt, IAtomContainer product) throws CDKException {
        Double eEnergy = 0.0;
        BondEnergies bondEnergy = BondEnergies.getInstance();
        for (int i = 0; i
                < Educt.getBondCount(); i++) {
            IBond bond = Educt.getBond(i);
            eEnergy += getBondEnergy(bond, bondEnergy);
        }
        Double pEnergy = 0.0;
        for (int j = 0; j
                < product.getBondCount(); j++) {
            IBond bond = product.getBond(j);
            pEnergy += getBondEnergy(bond, bondEnergy);
        }
        return (eEnergy + pEnergy);
    }

    private double getBondEnergy(IBond bond, BondEnergies bondEnergy) {
        double energy = 0.0;
        if ((bond.getAtom(0).getFlag(0) == true && bond.getAtom(1).getFlag(0) == false)
                || (bond.getAtom(0).getFlag(0) == false && bond.getAtom(1).getFlag(0) == true)) {
            Integer val = bondEnergy.getEnergies(bond.getAtom(0), bond.getAtom(1), bond.getOrder());
            if (val != null) {
                energy = val;
            }
        }
        return energy;
    }
}

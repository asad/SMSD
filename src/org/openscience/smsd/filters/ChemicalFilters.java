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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * A set of filters applied to the results.
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * @cdk.module smsd
 *
 */
@TestClass("org.openscience.cdk.smsd.filters.ChemicalFiltersTest")
public class ChemicalFilters {

    private List<Map<Integer, Integer>> allMCS = null;
    private Map<Integer, Integer> firstSolution = null;
    private List<Map<IAtom, IAtom>> allAtomMCS = null;
    private Map<IAtom, IAtom> firstAtomMCS = null;
    private IAtomContainer rMol = null;
    private IAtomContainer pMol = null;
    private IChemicalFilter<Double> energyFilter;
    private IChemicalFilter<Integer> fragmentFilter;
    private IChemicalFilter<Double> stereoFilter;

    public ChemicalFilters(List<Map<Integer, Integer>> allMCS,
            List<Map<IAtom, IAtom>> allAtomMCS,
            Map<Integer, Integer> firstSolution,
            Map<IAtom, IAtom> firstAtomMCS,
            IAtomContainer sourceMol,
            IAtomContainer targetMol) {
        this.allAtomMCS = allAtomMCS;
        this.allMCS = allMCS;
        this.firstAtomMCS = firstAtomMCS;
        this.firstSolution = firstSolution;
        this.rMol = sourceMol;
        this.pMol = targetMol;

        energyFilter = new EnergyFilter(rMol, pMol);
        fragmentFilter = new FragmentFilter(rMol, pMol);
        stereoFilter = new StereoFilter(rMol, pMol);
    }

    private void clear(Map<Integer, Map<Integer, Integer>> sortedAllMCS,
            Map<Integer, Map<IAtom, IAtom>> sortedAllAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Integer> fragmentScoreMap,
            Map<Integer, Double> energySelectionMap) {
        sortedAllMCS.clear();
        sortedAllAtomMCS.clear();
        stereoScoreMap.clear();
        fragmentScoreMap.clear();
        energySelectionMap.clear();
    }

    /**
     * Sort MCS solution by bond breaking energy.
     *
     * @throws CDKException
     */
    @TestMethod("testSortResultsByEnergies")
    public synchronized void sortResultsByEnergies() throws CDKException {
        Map<Integer, Map<Integer, Integer>> allEnergyMCS = new TreeMap<Integer, Map<Integer, Integer>>();
        Map<Integer, Map<IAtom, IAtom>> allEnergyAtomMCS = new TreeMap<Integer, Map<IAtom, IAtom>>();

        Map<Integer, Double> stereoScoreMap = new TreeMap<Integer, Double>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Double> energySelectionMap = new TreeMap<Integer, Double>();

        initializeMaps(allEnergyMCS, allEnergyAtomMCS, stereoScoreMap, fragmentScoreMap, energySelectionMap);

        double lowestEnergyScore = energyFilter.sortResults(
                allEnergyMCS, allEnergyAtomMCS, energySelectionMap);
        clear();

        int counter = 0;
        for (Map.Entry<Integer, Double> map : energySelectionMap.entrySet()) {
            if (lowestEnergyScore == map.getValue().doubleValue()) {
                addSolution(counter, map.getKey(),
                        allEnergyAtomMCS,
                        allEnergyMCS,
                        stereoScoreMap,
                        energySelectionMap,
                        fragmentScoreMap);
                counter++;
            }
        }

        if (lowestEnergyScore != EnergyFilter.MAX_ENERGY) {
            firstSolution.putAll(allMCS.get(0));
            firstAtomMCS.putAll(allAtomMCS.get(0));
            clear(allEnergyMCS, allEnergyAtomMCS, stereoScoreMap, fragmentScoreMap, energySelectionMap);
        }
    }

    /**
     * Sort solution by ascending order of the fragment count.
     */
    @TestMethod("testSortResultsByFragments")
    public synchronized void sortResultsByFragments() {
        Map<Integer, Map<Integer, Integer>> allFragmentMCS = new TreeMap<Integer, Map<Integer, Integer>>();
        Map<Integer, Map<IAtom, IAtom>> allFragmentAtomMCS = new TreeMap<Integer, Map<IAtom, IAtom>>();

        Map<Integer, Double> stereoScoreMap = new TreeMap<Integer, Double>();
        Map<Integer, Double> energyScoreMap = new TreeMap<Integer, Double>();
        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();

        initializeMaps(allFragmentMCS,
                allFragmentAtomMCS,
                stereoScoreMap,
                fragmentScoreMap,
                energyScoreMap);

        try {
            int minFragmentScore = fragmentFilter.sortResults(
                    allFragmentMCS, allFragmentAtomMCS, fragmentScoreMap);

            boolean flag = false;
            if (minFragmentScore < 9999) {
                flag = true;
                clear();
            }
            int counter = 0;
            for (Map.Entry<Integer, Integer> map : fragmentScoreMap.entrySet()) {
                if (minFragmentScore == map.getValue()) {
                    addSolution(counter, map.getKey(),
                            allFragmentAtomMCS,
                            allFragmentMCS,
                            stereoScoreMap,
                            energyScoreMap,
                            fragmentScoreMap);
                    counter++;
                }
            }

            if (flag) {
                firstSolution.putAll(allMCS.get(0));
                firstAtomMCS.putAll(allAtomMCS.get(0));
                clear(allFragmentMCS, allFragmentAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);
            }
        } catch (CDKException c) {
            // actually, never thrown, but in the interface
        }

    }

    /**
     * Sort MCS solution by stereo and bond type matches.
     * @throws CDKException
     */
    @TestMethod("testSortResultsByStereoAndBondMatch")
    public synchronized void sortResultsByStereoAndBondMatch() throws CDKException {
        Map<Integer, Map<Integer, Integer>> allStereoMCS = new TreeMap<Integer, Map<Integer, Integer>>();
        Map<Integer, Map<IAtom, IAtom>> allStereoAtomMCS = new HashMap<Integer, Map<IAtom, IAtom>>();

        Map<Integer, Integer> fragmentScoreMap = new TreeMap<Integer, Integer>();
        Map<Integer, Double> energyScoreMap = new TreeMap<Integer, Double>();
        Map<Integer, Double> stereoScoreMap = new HashMap<Integer, Double>();

        initializeMaps(allStereoMCS,
                allStereoAtomMCS,
                stereoScoreMap,
                fragmentScoreMap,
                energyScoreMap);

        double highestStereoScore = stereoFilter.sortResults(
                allStereoMCS, allStereoAtomMCS, stereoScoreMap);

        if (highestStereoScore != 0) {
            boolean flag = false;

            //Higher Score is mapped preferred over lower


            double secondhigestStereoScore = highestStereoScore;
            for (Integer key : stereoScoreMap.keySet()) {
                if (secondhigestStereoScore < highestStereoScore
                        && stereoScoreMap.get(key) > secondhigestStereoScore) {
                    secondhigestStereoScore = stereoScoreMap.get(key);
                } else if (secondhigestStereoScore == highestStereoScore
                        && stereoScoreMap.get(key) < secondhigestStereoScore) {
                    secondhigestStereoScore = stereoScoreMap.get(key);
                }
            }

            if (!stereoScoreMap.isEmpty()) {
                flag = true;
                clear();
            }

            /*Put back the sorted solutions*/

            int counter = 0;
            for (Integer I : stereoScoreMap.keySet()) {
                if (highestStereoScore == stereoScoreMap.get(I)) {
                    addSolution(counter, I,
                            allStereoAtomMCS,
                            allStereoMCS,
                            stereoScoreMap,
                            energyScoreMap,
                            fragmentScoreMap);
                    counter++;

                }
            }
            if (flag) {
                firstSolution.putAll(allMCS.get(0));
                firstAtomMCS.putAll(allAtomMCS.get(0));
                clear(allStereoMCS, allStereoAtomMCS, stereoScoreMap, fragmentScoreMap, energyScoreMap);
            }
        }
    }

    /**
     * Return sorted energy in ascending order.
     * @return sorted bond breaking energy
     */
    @TestMethod("testGetSortedEnergy")
    public List<Double> getSortedEnergy() {
        return energyFilter.getScores();
    }

    /**
     * Return sorted fragment in ascending order of the size.
     * @return sorted fragment count
     */
    @TestMethod("testGetSortedFragment")
    public List<Integer> getSortedFragment() {
        return fragmentFilter.getScores();
    }

    /**
     * Return Stereo matches in descending order.
     * @return sorted stereo matches
     */
    @TestMethod("testGetStereoMatches")
    public List<Double> getStereoMatches() {
        return stereoFilter.getScores();
    }

    private void initializeMaps(
            Map<Integer, Map<Integer, Integer>> sortedAllMCS,
            Map<Integer, Map<IAtom, IAtom>> sortedAllAtomMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Integer> fragmentScoreMap,
            Map<Integer, Double> energySelectionMap) {

        Integer Index = 0;
        for (Map<IAtom, IAtom> atomsMCS : allAtomMCS) {
            sortedAllAtomMCS.put(Index, atomsMCS);
            fragmentScoreMap.put(Index, 0);
            energySelectionMap.put(Index, 0.0);
            stereoScoreMap.put(Index, 0.0);
            Index++;
        }

        Index = 0;
        for (Map<Integer, Integer> MCS : allMCS) {
            sortedAllMCS.put(Index, MCS);
            Index++;
        }

        energyFilter.fillMap(energySelectionMap);
        fragmentFilter.fillMap(fragmentScoreMap);
        stereoFilter.fillMap(stereoScoreMap);

    }

    private void addSolution(int counter, int key,
            Map<Integer, Map<IAtom, IAtom>> allFragmentAtomMCS,
            Map<Integer, Map<Integer, Integer>> allFragmentMCS,
            Map<Integer, Double> stereoScoreMap,
            Map<Integer, Double> energyScoreMap,
            Map<Integer, Integer> fragmentScoreMap) {

        allAtomMCS.add(counter, allFragmentAtomMCS.get(key));
        allMCS.add(counter, allFragmentMCS.get(key));
        stereoFilter.addScore(counter, stereoScoreMap.get(key));
        fragmentFilter.addScore(counter, fragmentScoreMap.get(key));
        energyFilter.addScore(counter, energyScoreMap.get(key));

    }

    private void clear() {
        firstSolution.clear();
        allMCS.clear();
        allAtomMCS.clear();
        firstAtomMCS.clear();
        energyFilter.clearScores();
        fragmentFilter.clearScores();
        stereoFilter.clearScores();
    }
}

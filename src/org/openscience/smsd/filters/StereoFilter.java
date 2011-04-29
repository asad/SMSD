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
package org.openscience.smsd.filters;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.isomorphism.matchers.IQueryAtom;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

/**
 * Filter on stereo and bond matches.
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 * @cdk.module smsd
 */
public class StereoFilter extends BaseFilter implements IChemicalFilter<Double> {

    private List<Double> stereoScore = null;

    public StereoFilter(IAtomContainer rMol, IAtomContainer pMol) {
        super(rMol, pMol);
        stereoScore = new ArrayList<Double>();
    }

    @Override
    public Double sortResults(
            Map<Integer, Map<Integer, Integer>> allStereoMCS,
            Map<Integer, Map<IAtom, IAtom>> allStereoAtomMCS,
            Map<Integer, Double> stereoScoreMap) throws CDKException {

        getStereoBondChargeMatch(stereoScoreMap, allStereoMCS, allStereoAtomMCS);

        stereoScoreMap = sortMapByValueInDescendingOrder(stereoScoreMap);
        double highestStereoScore =
                stereoScoreMap.isEmpty() ? 0
                : stereoScoreMap.values().iterator().next();
        return highestStereoScore;
    }

    @Override
    public List<Double> getScores() {
        return Collections.unmodifiableList(stereoScore);
    }

    @Override
    public void clearScores() {
        stereoScore.clear();
    }

    @Override
    public void addScore(int counter, Double score) {
        stereoScore.add(counter, score);
    }

    @Override
    public void fillMap(Map<Integer, Double> stereoScoreMap) {
        int Index = 0;
        for (Double score : stereoScore) {
            stereoScoreMap.put(Index, score);
            Index++;
        }
    }

    private boolean getStereoBondChargeMatch(Map<Integer, Double> stereoScoreMap,
            Map<Integer, Map<Integer, Integer>> allStereoMCS,
            Map<Integer, Map<IAtom, IAtom>> allStereoAtomMCS) throws CDKException {

        boolean stereoMatchFlag = false;
        IAtomContainer reactant = rMol;
        IAtomContainer product = pMol;
        CDKHueckelAromaticityDetector.detectAromaticity(reactant);
        CDKHueckelAromaticityDetector.detectAromaticity(product);

        for (Integer Key : allStereoMCS.keySet()) {
            try {
                double score = 0.0;
                //            System.out.println("\nStart score " + score);
                Map<Integer, Integer> atomsMCS = allStereoMCS.get(Key);
                Map<IAtom, IAtom> atomMapMCS = allStereoAtomMCS.get(Key);
                double atomScore = getAtomScore(score, atomMapMCS, reactant, product);
                Map<IBond, IBond> bondMaps = makeBondMapsOfAtomMaps(rMol, pMol, atomsMCS);
                double ringScore = 0.0;
                if (rMol.getBondCount() > 1
                        && pMol.getBondCount() > 1) {
                    List<Object> subgraphRList = getMappedFragment(rMol, atomMapMCS.keySet());

                    double rscore = getRingMatchScore(subgraphRList);
                    List<Object> subgraphPList = getMappedFragment(pMol, atomMapMCS.values());
                    double pscore = getRingMatchScore(subgraphPList);
                    ringScore = rscore + pscore;
                }
                double bondScore = getBondScore(score, bondMaps);

                score = atomScore + ringScore + bondScore;
                if (!stereoMatchFlag) {
                    stereoMatchFlag = true;
                }
                stereoScoreMap.put(Key, score);
            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(StereoFilter.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return stereoMatchFlag;
    }

    private Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2,
            Map<Integer, Integer> mappings) {

        Map<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (IAtom atoms : ac1.atoms()) {

            int ac1AtomNumber = ac1.getAtomNumber(atoms);

            if (mappings.containsKey(ac1AtomNumber)) {

                int ac2AtomNumber = mappings.get(ac1AtomNumber);

                List<IAtom> connectedAtoms = ac1.getConnectedAtomsList(atoms);

                for (IAtom cAtoms : connectedAtoms) {
                    int ac1ConnectedAtomNumber = ac1.getAtomNumber(cAtoms);

                    if (mappings.containsKey(ac1ConnectedAtomNumber)) {
                        {
                            int ac2ConnectedAtomNumber = mappings.get(ac1ConnectedAtomNumber);

                            IBond ac1Bond = ac1.getBond(atoms, cAtoms);
                            IBond ac2Bond = ac2.getBond(ac2.getAtom(ac2AtomNumber),
                                    ac2.getAtom(ac2ConnectedAtomNumber));

                            if (ac2Bond == null) {
                                ac2Bond = ac2.getBond(ac2.getAtom(ac2ConnectedAtomNumber), ac2.getAtom(ac2AtomNumber));
                            }

                            if (ac1Bond != null && ac2Bond != null) {
                                maps.put(ac1Bond, ac2Bond);
                            }
                        }
                    }
                }
            }
        }
//        System.out.println("Mol Map size:" + maps.size());
        return maps;

    }

    private double getAtomScore(double score, Map<IAtom, IAtom> atomMapMCS, IAtomContainer reactant,
            IAtomContainer product) {
        for (Map.Entry<IAtom, IAtom> mappings : atomMapMCS.entrySet()) {
            IAtom rAtom = mappings.getKey();
            IAtom pAtom = mappings.getValue();

            int rHCount = 0;
            int pHCount = 0;
            double rBO = reactant.getBondOrderSum(rAtom);
            double pBO = product.getBondOrderSum(pAtom);

            if (rAtom.getImplicitHydrogenCount() != null) {
                rHCount = rAtom.getImplicitHydrogenCount();
            }
            if (pAtom.getImplicitHydrogenCount() != null) {
                pHCount = pAtom.getImplicitHydrogenCount();
            }

            int HScore = Math.abs(rHCount - pHCount);
            double BOScore = Math.abs(rBO - pBO);

            if (rHCount != pHCount) {
                score -= HScore;
            } else {
                score += HScore;
            }

            if (rBO != pBO) {
                score -= BOScore;
            } else {
                score += BOScore;
            }

            if (rAtom.getFormalCharge() == pAtom.getFormalCharge()) {
                score += 5.0;
            }
        }
        return score;
    }

    private double getBondScore(double score, Map<IBond, IBond> bondMaps) {
        for (Map.Entry<IBond, IBond> matchedBonds : bondMaps.entrySet()) {

            IBond RBond = matchedBonds.getKey();
            IBond PBond = matchedBonds.getValue();

            score += getBondTypeMatches(RBond, PBond);
        }
        return score;
    }

    private double getBondTypeMatches(IBond queryBond, IBond targetBond) {
        double score = 0;

        if (targetBond instanceof IQueryBond && queryBond instanceof IBond) {
            IQueryBond bond = (IQueryBond) targetBond;
            IQueryAtom atom1 = (IQueryAtom) (targetBond.getAtom(0));
            IQueryAtom atom2 = (IQueryAtom) (targetBond.getAtom(1));
            if (bond.matches(queryBond)) {
                // ok, bonds match
                if (atom1.matches(queryBond.getAtom(0)) && atom2.matches(queryBond.getAtom(1))
                        || atom1.matches(queryBond.getAtom(1)) && atom2.matches(queryBond.getAtom(0))) {
                    // ok, atoms match in either order
                    score += 4;
                }
            } else {
                score -= 4;
            }
        } else if (queryBond instanceof IQueryBond && targetBond instanceof IBond) {
            IQueryBond bond = (IQueryBond) queryBond;
            IQueryAtom atom1 = (IQueryAtom) (queryBond.getAtom(0));
            IQueryAtom atom2 = (IQueryAtom) (queryBond.getAtom(1));
            if (bond.matches(targetBond)) {
                // ok, bonds match
                if (atom1.matches(targetBond.getAtom(0)) && atom2.matches(targetBond.getAtom(1))
                        || atom1.matches(targetBond.getAtom(1)) && atom2.matches(targetBond.getAtom(0))) {
                    // ok, atoms match in either order
                    score += 4;
                }
            } else {
                score -= 4;
            }
        } else {

            int reactantBondType = convertBondOrder(queryBond);
            int productBondType = convertBondOrder(targetBond);
            int rStereo = convertBondStereo(queryBond);
            int pStereo = convertBondStereo(targetBond);
            if ((queryBond.getFlag(CDKConstants.ISAROMATIC) == targetBond.getFlag(CDKConstants.ISAROMATIC))
                    && (reactantBondType == productBondType)) {
                score += 8;
            } else if (queryBond.getFlag(CDKConstants.ISAROMATIC) && targetBond.getFlag(CDKConstants.ISAROMATIC)) {
                score += 4;
            }

            if (reactantBondType == productBondType) {
                score += productBondType;
            } else {
                score -= 4 * Math.abs(reactantBondType - productBondType);
            }

            if (rStereo != 4 || pStereo != 4 || rStereo != 3 || pStereo != 3) {
                if (rStereo == pStereo) {
                    score += 1;
                } else {
                    score -= 1;
                }
            }

        }
        return score;
    }

    /**
     * Get stereo value as integer
     * @param bond
     * @return
     */
    public static int convertBondStereo(IBond bond) {
        int value = 0;
        switch (bond.getStereo()) {
            case UP:
                value = 1;
                break;
            case UP_INVERTED:
                value = 1;
                break;
            case DOWN:
                value = 6;
                break;
            case DOWN_INVERTED:
                value = 6;
                break;
            case UP_OR_DOWN:
                value = 4;
                break;
            case UP_OR_DOWN_INVERTED:
                value = 4;
                break;
            case E_OR_Z:
                value = 3;
                break;
            default:
                value = 0;
        }
        return value;
    }

    /**
     *Get bond order value as integer
     * @param bond
     * @return
     */
    public static int convertBondOrder(IBond bond) {
        int value = 0;
        switch (bond.getOrder()) {
            case QUADRUPLE:
                value = 4;
                break;
            case TRIPLE:
                value = 3;
                break;
            case DOUBLE:
                value = 2;
                break;
            case SINGLE:
                value = 1;
                break;
            default:
                value = 0;
        }
        return value;
    }

    private double getRingMatchScore(List<Object> list) {
        double lScore = 0;
        List<IAtom> listMap = (List<IAtom>) list.get(0);
        IAtomContainer ac = (IAtomContainer) list.get(1);
//        HanserRingFinder ringFinder = new HanserRingFinder();
        IRingSet sssr = null;
        try {
            SSSRFinder finder = new SSSRFinder(ac);
            sssr = finder.findEssentialRings();
            RingSetManipulator.sort(sssr);
//            System.out.println("Ring length " + sssr.getAtomContainerCount());
            lScore = getRingMatch(sssr, listMap);
        } catch (Exception ex) {
            Logger.getLogger(StereoFilter.class.getName()).log(Level.SEVERE, null, ex);
        }
        return lScore;
    }

    private double getRingMatch(IRingSet rings, List<IAtom> atoms) {
        double score = 0.0;
        for (IAtom a : atoms) {
            for (IAtomContainer ring : rings.atomContainers()) {
                if (ring.contains(a)) {
                    score += 10;
                } else {
                    score -= 10;
                }
            }
        }
        return score;
    }

    private List<Object> getMappedFragment(IAtomContainer molecule, Collection<IAtom> atomsMCS) throws CloneNotSupportedException {
        IAtomContainer subgraphContainer = molecule.getBuilder().newInstance(IAtomContainer.class, molecule);
        List<IAtom> list = new ArrayList<IAtom>(atomsMCS.size());
        for (IAtom atom : atomsMCS) {
            int post = molecule.getAtomNumber(atom);
//            System.out.println("Atom to be removed " + post);
            list.add(subgraphContainer.getAtom(post));
        }

        List<IAtom> rlist = new ArrayList<IAtom>();
        for (IAtom atoms : subgraphContainer.atoms()) {
            if (!list.contains(atoms)) {
                rlist.add(atoms);
            }
        }

        for (IAtom atoms : rlist) {
            subgraphContainer.removeAtomAndConnectedElectronContainers(atoms);
        }
        List<Object> l = new ArrayList<Object>();
        l.add(list);
        l.add(subgraphContainer);
        return l;
    }
}

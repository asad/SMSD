/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.gui;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;
import java.util.List;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 *  java1.8+
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 */
public class SubGraphMapper {

    /**
     *
     * @param mainMolecule
     * @param matchedCore
     * @return container
     */
    public static IAtomContainer getNeedle(IAtomContainer mainMolecule, List<Integer> matchedCore) {
        IAtomContainer needle = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);

        assert matchedCore != null;
//        System.out.println("Number of matched atoms = " + matchedCore.size());
        for (int i = 0; i < matchedCore.size(); i++) {
            for (int j = 0; j < matchedCore.size(); j++) {

                int iIndex = matchedCore.get(i);
                int jIndex = matchedCore.get(j);

                if (i != j && iIndex != jIndex) {
                    IBond bond = mainMolecule.getBond(mainMolecule.getAtom(iIndex), mainMolecule.getAtom(jIndex));
                    if (bond != null) {
                        needle.addBond(bond);
                    }
                }
            }
        }

        return needle;
    }

    public static List<List<Integer>> getNeedleList(List<Integer> matchedCore) {
        List<List<Integer>> matchedList = new ArrayList<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (int i = 0; i < matchedCore.size(); i += 2) {

            rList.add(matchedCore.get(i));
            pList.add(matchedCore.get(i + 1));

        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

    public static List<List<Integer>> getNeedleList(Map<Integer, Integer> matchedCore) {
        List<List<Integer>> matchedList = new ArrayList<List<Integer>>();
        List<Integer> rList = new LinkedList<Integer>();
        List<Integer> pList = new LinkedList<Integer>();

        for (Map.Entry<Integer, Integer> M : matchedCore.entrySet()) {
            rList.add(M.getKey());
            pList.add(M.getValue());
        }

        matchedList.add(0, rList);
        matchedList.add(1, pList);
        return matchedList;
    }

    public static Map<IBond, IBond> makeBondMapsOfAtomMaps(IAtomContainer ac1, IAtomContainer ac2, TreeMap<Integer, Integer> mappings) {

        HashMap<IBond, IBond> maps = new HashMap<IBond, IBond>();

        for (Map.Entry<Integer, Integer> mapS : mappings.entrySet()) {
            int indexI = mapS.getKey();
            int indexJ = mapS.getValue();

            for (Map.Entry<Integer, Integer> mapD : mappings.entrySet()) {

                int indexIPlus = mapD.getKey();
                int indexJPlus = mapD.getValue();

                if (indexI != indexIPlus && indexJ != indexJPlus) {

                    IBond ac1Bond = ac1.getBond(ac1.getAtom(indexI), ac1.getAtom(indexIPlus));

                    if (ac1Bond != null) {

                        IBond ac2Bond = ac2.getBond(ac2.getAtom(indexJ), ac2.getAtom(indexJPlus));

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

    /**
     * This function will returns a subgraph extracted from the source or target
     * molecules
     *
     * @param container source/target container
     * @param matches source/target mapping
     * @return mapped subgraph/substructure
     */
    public IAtomContainer getMatchedSubgraph(IAtomContainer container, Collection<Integer> matches) {
        IAtomContainer needle = container.getBuilder().newInstance(IAtomContainer.class, container);
        List<IAtom> atomListToBeRemoved = new ArrayList<>(matches.size());
        matches.stream().map((index) -> needle.getAtom(index)).forEach((atomToBeRemoved) -> {
            atomListToBeRemoved.add(atomToBeRemoved);
        });

        atomListToBeRemoved.stream().forEach((removeAtom) -> {
            needle.removeAtom(removeAtom);
        });
        atomListToBeRemoved.clear();
        return needle;
    }
}

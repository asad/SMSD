/* Copyright (C) 2009-2015  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: Syed Asad Rahman <asad@ebi.ac.uk>
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
package gui;

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
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK. e-mail: asad@ebi.ac.uk
 */
public class SubGraphMapper {

    /**
     *
     * @param mainMolecule
     * @param matchedCore
     * @return
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
     * This function will returns a subgraph extracted from the source or target molecules
     *
     * @param container source/target container
     * @param matches source/target mapping
     * @return mapped subgraph/substructure
     */
    public IAtomContainer getMatchedSubgraph(IAtomContainer container, Collection<Integer> matches) {
        IAtomContainer needle = container.getBuilder().newInstance(IAtomContainer.class, container);
        List<IAtom> atomListToBeRemoved = new ArrayList<IAtom>(matches.size());
        for (Integer index : matches) {
            IAtom atomToBeRemoved = needle.getAtom(index.intValue());
            atomListToBeRemoved.add(atomToBeRemoved);
        }

        for (IAtom removeAtom : atomListToBeRemoved) {
            needle.removeAtomAndConnectedElectronContainers(removeAtom);
        }
        atomListToBeRemoved.clear();
        return needle;
    }
}

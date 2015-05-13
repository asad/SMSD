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

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;

public class MCS {

    public static IAtomContainer getMcsAsNewContainer(IAtomContainer mol1, IAtomContainer mol2) throws CDKException, CloneNotSupportedException, IOException {
        Isomorphism mcs = new Isomorphism(mol1, mol2, org.openscience.smsd.interfaces.Algorithm.DEFAULT, true, false, true);
        mcs.setChemFilters(true, true, true);

        System.out.println("MCS Size: " + mcs.getTanimotoSimilarity());
        System.out.println("MCS First Map: " + mcs.getFirstAtomMapping().toString());
        System.out.println("MCS First size: " + mcs.getFirstAtomMapping().getCount());

        mol1 = mcs.getQuery();
        mol2 = mcs.getTarget();

        IAtomContainer mcsmolecule = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class, mol1);
        Map<Integer, Integer> indexMapping = getIndexMapping(mcs.getFirstAtomMapping());

        List<IAtom> atomsToBeRemoved = new ArrayList<>();
        for (IAtom atom : mcsmolecule.atoms()) {
            int index = mcsmolecule.getAtomNumber(atom);
//            System.out.println("index: " +index);
            if (!indexMapping.containsKey(index)) {
                atomsToBeRemoved.add(atom);
            }
        }

        for (IAtom atom : atomsToBeRemoved) {
            mcsmolecule.removeAtomAndConnectedElectronContainers(atom);
        }

        StringWriter stringWriter = new StringWriter();
        try (SMILESWriter smilesWriter = new SMILESWriter(stringWriter)) {
            smilesWriter.write(mcsmolecule);
            smilesWriter.close();
        }
        System.out.println("MCS SMILES: " + stringWriter.toString());

        return mcsmolecule.clone();
    }

    public static void main(String[] args) throws CDKException, CloneNotSupportedException, IOException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Men_12
        IAtomContainer A1 = sp.parseSmiles("c1cccc(COC(=O)NC(CC(C)C)C(=O)NC(CCc2ccccc2)C(=O)COC)c1");
        // Tri_06
        IAtomContainer A2 = sp.parseSmiles("c1cccc(COC(=O)NC(CC(C)C)C(=O)NCC#N)c1");

        getMcsAsNewContainer(A2, A1);

        getMcsAsNewContainer(A1, A2);
    }

    private static Map<Integer, Integer> getIndexMapping(AtomAtomMapping aam) {
        Map<IAtom, IAtom> mappings = aam.getMappingsByAtoms();
        Map<Integer, Integer> mapping = new TreeMap<>();
        for (IAtom keys : mappings.keySet()) {
            mapping.put(aam.getQueryIndex(keys), aam.getTargetIndex(mappings.get(keys)));
        }
        return mapping;
    }
}

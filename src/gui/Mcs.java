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
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;

public class Mcs {

    public static IAtomContainer getMcsAsNewContainer(IAtomContainer mol1, IAtomContainer mol2) throws CDKException, CloneNotSupportedException, IOException {
        Isomorphism mcs = new Isomorphism(mol1, mol2, org.openscience.smsd.interfaces.Algorithm.DEFAULT, true);
        mcs.setChemFilters(true, true, true);

        System.out.println("MCS Size: " + mcs.getTanimotoSimilarity());
        System.out.println("MCS First Map: " + mcs.getFirstAtomMapping().toString());
        System.out.println("MCS First size: " + mcs.getFirstAtomMapping().getCount());

        mol1 = mcs.getQueryContainer();
        mol2 = mcs.getTargetContainer();

        IMolecule mcsmolecule = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, mol1);
        Map<Integer, Integer> indexMapping = getIndexMapping(mcs.getFirstAtomMapping());

        List<IAtom> atomsToBeRemoved = new ArrayList<IAtom>();
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
        SMILESWriter smilesWriter = new SMILESWriter(stringWriter);
        smilesWriter.write(mcsmolecule);
        smilesWriter.close();

        System.out.println("MCS SMILES: " + stringWriter.toString());

        return (IAtomContainer) mcsmolecule.clone();
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
        Map<IAtom, IAtom> mappings = aam.getMappings();
        Map<Integer, Integer> mapping = new TreeMap<Integer, Integer>();
        for (IAtom keys : mappings.keySet()) {
            mapping.put(aam.getQueryIndex(keys), aam.getTargetIndex(mappings.get(keys)));
        }
        return mapping;
    }
}

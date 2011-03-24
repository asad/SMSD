package gui;

import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.Isomorphism;

public class Mcs {

    public static IAtomContainer getMcsAsNewContainer(IAtomContainer mol1, IAtomContainer mol2) throws CDKException, CloneNotSupportedException, IOException {
        Isomorphism mcs = new Isomorphism(org.openscience.smsd.interfaces.Algorithm.DEFAULT, true);
        mcs.init(mol1, mol2);
        mcs.setChemFilters(true, true, true);

        System.out.println("MCS Size: " + mcs.getTanimotoAtomSimilarity());
        System.out.println("MCS First Map: " + mcs.getFirstMapping());
        System.out.println("MCS First size: " + mcs.getFirstAtomMapping().size());

        mol1 = mcs.getReactantMolecule();
        mol2 = mcs.getProductMolecule();

        IMolecule mcsmolecule = DefaultChemObjectBuilder.getInstance().newInstance(IMolecule.class, mol1);

        List<IAtom> atomsToBeRemoved = new ArrayList<IAtom>();
        for (IAtom atom : mcsmolecule.atoms()) {
            int index = mcsmolecule.getAtomNumber(atom);
//            System.out.println("index: " +index);
            if (!mcs.getFirstMapping().containsKey(index)) {
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

    @Test
    public static void main(String[] args) throws CDKException, CloneNotSupportedException, IOException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        // Men_12
        IAtomContainer A1 = sp.parseSmiles("c1cccc(COC(=O)NC(CC(C)C)C(=O)NC(CCc2ccccc2)C(=O)COC)c1");
        // Tri_06
        IAtomContainer A2 = sp.parseSmiles("c1cccc(COC(=O)NC(CC(C)C)C(=O)NCC#N)c1");

        getMcsAsNewContainer(A2, A1);

        getMcsAsNewContainer(A1, A2);
    }
}

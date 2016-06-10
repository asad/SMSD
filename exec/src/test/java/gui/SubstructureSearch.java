package gui;

import java.util.Map;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;


/*
 * This example file covers the substructure search part of the SMSD. A result
 * is obtained if query is a subgraph of template structure. Query molecule is
 * always the first molecule and template is the second molecule.
 *
 * SubstructureSearch.java @modified 4 April, 2011, 11.06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK contact asad@ebi.ac.uk
 *
 */
public class SubstructureSearch {

    /**
     *
     * @param args
     */
    public static void main(String[] args) {

        try {

            String query = "CC";
            String target = "C1CCC12CCCC2";
            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

            IAtomContainer mol1 = sp.parseSmiles(query);
            IAtomContainer mol2 = sp.parseSmiles(target);

            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);

            mol1 = ExtAtomContainerManipulator.removeHydrogens(mol1);
            mol2 = ExtAtomContainerManipulator.removeHydrogens(mol2);

//	Calling the main algorithm to perform MCS cearch
            ExtAtomContainerManipulator.aromatizeDayLight(mol1);
            ExtAtomContainerManipulator.aromatizeDayLight(mol2);

            boolean bondSensitive = true;
            boolean ringMatcher = true;
            boolean allMatch = false;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            /*
             * Test if mol1 is a substructure of mol2
             */
            Substructure comparison = new Substructure(mol1, mol2, bondSensitive, ringMatcher, true, allMatch);
            ////System.out.println("Is a Subgraph " + comparison.isSubgraph());
            if (comparison.isSubgraph()) {
                comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

                AtomAtomMapping mapping = comparison.getFirstAtomMapping();
                for (Map.Entry<IAtom, IAtom> aam : mapping.getMappingsByAtoms().entrySet()) {
                    //Get the mapped atom number in Query Molecule
                    int queryMappingNumber = mapping.getQueryIndex(aam.getKey());
                    //Get the mapped atom number in Target Molecule
                    int targetMappingNumber = mapping.getTargetIndex(aam.getValue());

                    //Get the mapped atom in Query Molecule
                    IAtom queryAtom = aam.getKey();
                    //Get the mapped atom in Target Molecule
                    IAtom targetAtom = aam.getValue();
                    //Print mapped atoms
                    System.out.println(
                            queryAtom.getSymbol()
                            + "(" + queryMappingNumber + "), "
                            + targetAtom.getSymbol()
                            + "(" + targetMappingNumber + ")");

                }
                ////System.out.println("");
                ////System.out.println("");

                ////System.out.println("Stereo Match: " + comparison.getStereoScore(0));
                ////System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                ////System.out.println("Fragment Size: " + comparison.getFragmentSize(0));
                ////System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                ////System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                ////System.out.println("");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        ////System.out.println("");
        ////System.out.println("");
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     */
    private static void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        ////System.out.println("Molecule 1");
        for (int i = 0; i
                < Molecule1.getAtomCount(); i++) {
            //System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        ////System.out.println();
        ////System.out.println("Molecule 2");
        for (int i = 0; i
                < Molecule2.getAtomCount(); i++) {
            //System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }

        ////System.out.println();
    }
}

package example;

import java.util.Map;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.tools.GraphMolecule;


/* This example file covers the substructure search part of the SMSD.
 * A result is obtained if query is a subgraph of template structure.
 * Query molecule is always the first molecule and template is the second molecule.
 *
 * SubstructureSearch.java
 * @modified 4 April, 2011, 11.06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * contact asad@ebi.ac.uk
 *
 */
public class SubstructureSearch {

    /**
     * 
     * @param args
     * @throws CDKException  
     */
    public static void main(String[] args) throws CDKException {




        String query = "CC";
        String target = "C1CCC12CCCC2";
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IAtomContainer mol1 = sp.parseSmiles(query);
        IAtomContainer mol2 = sp.parseSmiles(target);

        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);

        mol1 = AtomContainerManipulator.removeHydrogens(mol1);
        mol2 = AtomContainerManipulator.removeHydrogens(mol2);

//	Calling the main algorithm to perform MCS cearch

        CDKHueckelAromaticityDetector.detectAromaticity(mol1);
        CDKHueckelAromaticityDetector.detectAromaticity(mol2);

        mol1 = new GraphMolecule(mol1);
        mol2 = new GraphMolecule(mol2);

        boolean bondSensitive = true;
        boolean stereoMatch = true;
        boolean fragmentMinimization = true;
        boolean energyMinimization = true;

        Substructure comparison = new Substructure(mol1, mol2, bondSensitive, false, true);

        if (comparison.isSubgraph()) {
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);
            int count_final_sol = 0;
            System.out.println("Output of the final Mappings: ");

            if (!comparison.getAllAtomMapping().isEmpty()) {
                for (AtomAtomMapping aams : comparison.getAllAtomMapping()) {
                    int final_solution_size = aams.getCount();
                    System.out.println("Final mapping Nr. " + (count_final_sol + 1) + " Size:" + final_solution_size);

                    for (Map.Entry<IAtom, IAtom> mapping : aams.getMappings().entrySet()) {


                        //Get the mapped atom in Query Molecule
                        IAtom queryAtom = mapping.getKey();
                        //Get the mapped atom in Target Molecule
                        IAtom targetAtom = mapping.getValue();


                        //Get the mapped atom number in Query Molecule
                        int queryMappingNumber = aams.getQueryIndex(queryAtom);
                        //Get the mapped atom number in Target Molecule
                        int targetMappingNumber = aams.getTargetIndex(targetAtom);
                        //Print mapped atom numbers
                        System.out.println(queryMappingNumber + " "
                                + (targetMappingNumber));
                        //Print mapped atoms
                        System.out.println(queryAtom.getSymbol() + " "
                                + targetAtom.getSymbol());
                    }
                    System.out.println("");

                    System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol));
                    System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol));
                    System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                    System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                    count_final_sol++;

                }


                System.out.println("");
            }
        }
        System.out.println("");
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     */
    private static void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        System.out.println("Molecule 1");

        for (int i = 0; i
                < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        System.out.println();
        System.out.println("Molecule 2");
        for (int i = 0; i
                < Molecule2.getAtomCount(); i++) {

            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }

        System.out.println();

    }
}
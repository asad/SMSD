package example;

import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;


/*
 * This example file covers the MCS search part of the SMSD.
 * A result is obtained if query and template structure has atleast one common MCS.
 * Query molecule is always the first molecule and template is the second molecule.
 * Use the isSubgraph() method to find if query is a subgraph/substructure of the template mol.
 *
 *
 *
 * MCSSearch.java
 * @modified 5 November, 2011, 11.06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * contact asad@ebi.ac.uk
 *
 */
public class RingTest {

    /**
     * Creates a new instance of SMSD
     */
    public RingTest() {
        boolean first_MCS = false;
        try {
            String query = "C1=CC=C2C(=C1)C=CC=N2";//Quinoline 
            String target = "C1=CC=C2C(=C1)C=CN2";//Indole

            SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

            IAtomContainer mol1 = sp.parseSmiles(query);
            IAtomContainer mol2 = sp.parseSmiles(target);

            /*
             * perceive Atom Types
             */
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol1);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);

            /*
             * Remove Hydrogens
             */
            mol1 = ExtAtomContainerManipulator.removeHydrogens(mol1);
            mol2 = ExtAtomContainerManipulator.removeHydrogens(mol2);

            /*
             * Detect aromatic compounds and rings
             */
            ExtAtomContainerManipulator.aromatizeDayLight(mol1);
            ExtAtomContainerManipulator.aromatizeDayLight(mol2);

            /*
             * Calling the main algorithm to perform MCS search
             */
            boolean bondSensitive = true;
            boolean ringMatch = true;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            Isomorphism comparison = new Isomorphism(mol1, mol2, Algorithm.DEFAULT, bondSensitive, ringMatch, true);
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

            //Print all MCS solutions if first_MCS is false
            if (first_MCS == false) {
                int count_final_sol = 0;
                ////System.out.println("Output of the final Mappings: ");
                try {
                    if (!comparison.getAllAtomMapping().isEmpty()) {
                        for (AtomAtomMapping aams : comparison.getAllAtomMapping()) {
                            int final_solution_size = aams.getCount();
                            ////System.out.println("Final mapping Nr. " + (count_final_sol + 1) + " Size:" + final_solution_size);

                            for (Map.Entry<IAtom, IAtom> mapping : aams.getMappingsByAtoms().entrySet()) {

                                //Get the mapped atom in Query AtomContainer
                                IAtom queryAtom = mapping.getKey();
                                //Get the mapped atom in Target AtomContainer
                                IAtom targetAtom = mapping.getValue();

                                //Get the mapped atom number in Query AtomContainer
                                int queryMappingNumber = aams.getQueryIndex(queryAtom);
                                //Get the mapped atom number in Target AtomContainer
                                int targetMappingNumber = aams.getTargetIndex(targetAtom);
                                System.out.println(queryAtom.getSymbol() + "(" + queryMappingNumber + "), "
                                        + targetAtom.getSymbol() + "(" + targetMappingNumber + ") ");
                            }
                            ////System.out.println("");

                            ////System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol));
                            ////System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                            ////System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol));
                            ////System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                            ////System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                            count_final_sol++;

                        }

                        ////System.out.println("");
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }

            //Output only first solution
            if (first_MCS) {
                try {
                    AtomAtomMapping firstAtomMapping = comparison.getFirstAtomMapping();
                    for (Map.Entry<IAtom, IAtom> mapping : firstAtomMapping.getMappingsByAtoms().entrySet()) {

                        //Get the mapped atom in Query AtomContainer
                        IAtom queryAtom = mapping.getKey();
                        //Get the mapped atom in Target AtomContainer
                        IAtom targetAtom = mapping.getValue();
                        int queryMappingNumber = firstAtomMapping.getQueryIndex(queryAtom);
                        //Get the mapped atom number in Target AtomContainer
                        int targetMappingNumber = firstAtomMapping.getTargetIndex(targetAtom);
//                                Print mapped atom numbers
                        System.out.println(queryAtom.getSymbol() + "(" + queryMappingNumber + "), "
                                + targetAtom.getSymbol() + "(" + targetMappingNumber + ") ");
                    }
                    ////System.out.println("");

                    ////System.out.println("");
                    ////System.out.println("Stereo Match: " + comparison.getStereoScore(0));
                    ////System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    ////System.out.println("Fragment Size: " + comparison.getFragmentSize(0));
                    ////System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                    ////System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                    ////System.out.println("");
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                ////System.out.println("");
            }
            ////System.out.println("");

        } catch (Exception ex) {
            Logger.getLogger(RingTest.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * @param args the command line arguments
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {

        new RingTest();
    }
}

package example;

import java.io.File;
import java.io.FileInputStream;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
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
public class MCSSearch {

    /**
     * Creates a new instance of SMSD
     */
    public MCSSearch() {
    }

    /**
     * @param args the command line arguments
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        try {

            String mol1 = "Data/ATP.mol";
            String mol2 = "Data/ADP.mol";

            boolean first_MCS = true;

            boolean exists = (new File(mol1)).exists();
            if (!exists) {
                //System.err.println("Error: The Assigned File Path is not Correct " + mol1);
                //System.exit(1);
            }

            exists = (new File(mol2)).exists();
            if (!exists) {
                //System.err.println("Error: The Assigned File Path is not Correct " + mol2);
                //System.exit(1);
            }

            MDLV2000Reader molQuery = new MDLV2000Reader(new FileInputStream(mol1));
            IAtomContainer query = molQuery.read(new AtomContainer());

            MDLV2000Reader molTarget = new MDLV2000Reader(new FileInputStream(mol2));
            IAtomContainer target = molTarget.read(new AtomContainer());

            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(query);
            ExtAtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);

            query = ExtAtomContainerManipulator.removeHydrogens(query);
            target = ExtAtomContainerManipulator.removeHydrogens(target);

//	Calling the main algorithm to perform MCS cearch
            ExtAtomContainerManipulator.aromatizeMolecule(query);
            ExtAtomContainerManipulator.aromatizeMolecule(target);

            query = new AtomContainer(query);
            target = new AtomContainer(target);

            boolean bondSensitive = true;
            boolean ringMatch = false;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            Isomorphism comparison = new Isomorphism(query, target, Algorithm.DEFAULT, bondSensitive, ringMatch, false);
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
//                                Print mapped atom numbers
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
                        //Print mapped atom numbers
                        //Get the mapped atom number in Query AtomContainer
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
            Logger.getLogger(MCSSearch.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

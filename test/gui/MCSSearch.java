package gui;

import java.io.File;
import java.io.FileInputStream;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.AtomContainer;
import org.openscience.smsd.tools.Utility;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;


/*
 * This example file covers the MCS search part of the SMSD. A result is
 * obtained if query and template structure has atleast one common MCS. Query
 * molecule is always the first molecule and template is the second molecule.
 * Use the isSubgraph() method to find if query is a subgraph/substructure of
 * the template mol.
 *
 *
 *
 * MCSSearch.java @modified 27 June, 2012, 11.06 AM
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK contact asad@ebi.ac.uk
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

            boolean first_MCS = false;

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

            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(query);
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(target);

            query = AtomContainerManipulator.removeHydrogens(query);
            target = AtomContainerManipulator.removeHydrogens(target);

//	Calling the main algorithm to perform MCS cearch
            Utility.aromatizeDayLight(query);
            Utility.aromatizeDayLight(target);

            boolean bondSensitive = true;
            boolean ringmatch = false;
            boolean stereoMatch = true;
            boolean fragmentMinimization = true;
            boolean energyMinimization = true;

            Isomorphism comparison = new Isomorphism(query, target, Algorithm.DEFAULT, bondSensitive, ringmatch, true);
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

            //Print all MCS solutions if first_MCS is false
            if (!first_MCS) {
                int count_final_sol = 0;
                //System.out.println("Output of the final Mappings: ");
                try {
                    if (comparison.getAllAtomMapping() != null) {
                        for (AtomAtomMapping mapping : comparison.getAllAtomMapping()) {
                            int final_solution_size = mapping.getCount();
                            //System.out.println("Final mapping Nr. " + (count_final_sol + 1) + " Size:" + final_solution_size);

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
                                //System.out.println(
//                                        queryAtom.getSymbol()
//                                        + "(" + queryMappingNumber + "), "
//                                        + targetAtom.getSymbol()
//                                        + "(" + targetMappingNumber + ")");
                            }
                            //System.out.println("");

                            //System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol));
                            //System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                            //System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol));
                            //System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                            //System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                            count_final_sol++;

                        }

                        //System.out.println("");
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }

            //Output only first solution
            if (first_MCS) {
                try {
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
                    //System.out.println("");

                    //System.out.println("");
                    //System.out.println("Stereo Match: " + comparison.getStereoScore(0));
                    //System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    //System.out.println("Fragment Size: " + comparison.getFragmentSize(0));
                    //System.out.println("Tanimoto Similarity Score: " + comparison.getTanimotoSimilarity());
                    //System.out.println("Tanimoto Euclidean Distance: " + comparison.getEuclideanDistance());
                    //System.out.println("");
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                //System.out.println("");
            }
            //System.out.println("");

            //System.out.println("");
        } catch (Exception ex) {
            Logger.getLogger(MCSSearch.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     *
     * @param Molecule1
     * @param Molecule2
     */
    private static void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        //System.out.println("Molecule 1");
        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            //System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        //System.out.println();
        //System.out.println("Molecule 2");
        for (int i = 0; i < Molecule2.getAtomCount(); i++) {

            //System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }
        //System.out.println();

    }
}

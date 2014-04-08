package example;

import java.util.Map;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.smsd.tools.Utility;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Substructure;


/* This example file covers the substructure search part of the SMSD.
 * A result is obtained if query is a subgraph of template structure.
 * Query molecule is always the first molecule and template is the second molecule.
 *
 * SubstructureSearch.java
 * @modified 5 November, 2011, 11.06 AM
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
        Utility.aromatizeDayLight(mol1);
        Utility.aromatizeDayLight(mol2);

        mol1 = new AtomContainer(mol1);
        mol2 = new AtomContainer(mol2);

        boolean bondSensitive = true;
        boolean stereoMatch = true;
        boolean fragmentMinimization = true;
        boolean energyMinimization = true;

        Substructure comparison = new Substructure(mol1, mol2, bondSensitive, false, true,false);

        if (comparison.isSubgraph()) {
            comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);
            int count_final_sol = 0;
            ////System.out.println("Output of the final Mappings: ");

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
                        //Print mapped atom numbers
                        System.out.println(queryMappingNumber + " "
                                + (targetMappingNumber));
                        //Print mapped atoms
                        System.out.println(queryAtom.getSymbol() + " "
                                + targetAtom.getSymbol());
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
        }
        ////System.out.println("");
    }
}

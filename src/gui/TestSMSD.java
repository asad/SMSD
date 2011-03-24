/*
 * TestIsomorphism.java
 *
 * Created on January 28, 2007, 2:06 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */
package gui;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.IChemObjectReader;
import org.openscience.cdk.io.MDLReader;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.isomorphism.matchers.IQueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainer;
import org.openscience.cdk.isomorphism.matchers.QueryAtomContainerCreator;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.Substructure;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.GraphMolecule;

/**
 *
 * @author Syed Asad Rahman, EMBL-EBI, Cambridge, UK
 * @contact asad@ebi.ac.uk
 */
public class TestSMSD {

    static private ImageGenerator imageGenerator = null;

    /** Creates a new instance of TestIsomorphism */
    public TestSMSD() {
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     * @throws org.openscience.cdk.exception.CDKException
     * @throws Exception 
     */
    public static void main(String[] args) throws IOException, CDKException, Exception {

        MCSTest();
        //        testIsomorphismAdpAtpSubgraph();
        //        testAnyAtomAnyBondCase();
        //        testS();
        //        testQueryAtomContainerDefault();
        //        testCOAExampleMols();
        //          testCoverage();

    }

    /**
     * 
     * @param Molecule1
     * @param Molecule2
     */
    private static void printMolecules(IAtomContainer Molecule1, IAtomContainer Molecule2) {

        System.out.println("GraphMolecule 1");

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol() + " ");
        }

        System.out.println();
        System.out.println("GraphMolecule 2");
        for (int i = 0; i < Molecule2.getAtomCount(); i++) {

            System.out.print(Molecule2.getAtom(i).getSymbol() + " ");
        }
        System.out.println();

    }

    private static void printMolecules(IAtomContainer Molecule1) {

        System.out.println("GraphMolecule 1: " + Molecule1.getAtomCount());

        for (int i = 0; i < Molecule1.getAtomCount(); i++) {

            System.out.print(Molecule1.getAtom(i).getSymbol());
//            System.out.print(Molecule1.getAtom(i).getSymbol() + "[" + Molecule1.getAtom(i).getProperty(CANONICAL_LABEL) + "]");

        }

        System.out.println();


    }

    public static void MCSTest() throws CDKException, Exception {

//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/MCSTest/Benchmark/1.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/MCSTest/Benchmark/30.mol");

        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ADP.mol");
        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ATP.mol");



//        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
//        sdg.setMolecule(new GraphMolecule(query));
//        sdg.generateCoordinates();
//        query = sdg.getMolecule();
//
//        sdg = new StructureDiagramGenerator();
//        sdg.setMolecule(new GraphMolecule(target));
//        sdg.generateCoordinates();
//        target = sdg.getMolecule();


        Isomorphism smsd = new Isomorphism(Algorithm.VFLibMCS, true);
        smsd.init(query, target);
        smsd.setChemFilters(false, false, false);

        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();

        System.out.println("Mol1 Size. " + query.getAtomCount());
        System.out.println("Mol2 Size. " + target.getAtomCount());
        System.out.println("matched Size. " + smsd.getTanimotoSimilarity());
        System.out.println("Matches: " + smsd.getFirstAtomMapping().size());

//        printMapping(smsd, query, target);
        generateImage("VFLibMCS", query, target, smsd);


        smsd = new Isomorphism(Algorithm.CDKMCS, true);
        smsd.init(query, target);
        smsd.setChemFilters(false, false, false);

        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();

        System.out.println("Mol1 Size. " + query.getAtomCount());
        System.out.println("Mol2 Size. " + target.getAtomCount());
        System.out.println("matched Size. " + smsd.getTanimotoSimilarity());
        System.out.println("Matches: " + smsd.getFirstAtomMapping().size());

//        printMapping(smsd, query, target);
        generateImage("CDKMCS", query, target, smsd);

    }

    public static void testIsomorphismAdpAtpSubgraph() throws Exception {

        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());


//        String target = "CC(C)(COP([O-])([O-])(O*)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N)C(O)C(=O)NCCC(=O)NCCS";
//        String query = "CC(C)(COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C1N=CN=C2N)C(O)C(=O)NCCC(=O)NCCS";
////////
////        String query = "Nc1ccccc1";
////        String target = "C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C";
////        //Phosphate matches
//        String query = "[H]OC1C(COP([O-])(=O)OP([O-])(*)=O)OC(C1O)n1cnc2c(N)ncnc12";
//        String target = "[H]C(O)OC1C(COP([O-])(=O)OP([O-])(*)=O)OC(C1O)n1cnc2c(N)ncnc12";
//
////        String query= "[O-][CH]=[O]";//"CC(=O)C=O";//"O=C(CN1C=NC=CC1=O)NCC2=CC=CC=C2";
////        String target= "[O-][CH]=[O]";//"C[C@@H](O)C-O";//"CC(C)(C)N1C2=C(C=N1)C(=O)N(C=N2)CC(=O)NCC3=CC=CC=C3Cl";
//
//        //String query = "[H][O][C@]([H])([C]([H])=[O])[C]([H])([H])[H]";
//        String query = "[H][O][C]([H])([H])[C](=[O])[C]([H])([H])[O][P]([O-])([O-])=[O]";
//        String target = "[H][O][C@]1([H])[C@]([H])([O][P]([O-])([O-])=[O])[O][C@@]([H])([C]([H])([H])[H])[C@@]([H])([O][H])[C@@]1([H])[O][H]";
//
//        String query="[H]C([H])(C([O-])=O)C([H])([H])C(=O)SCCN";
//        String target ="[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)SCC";
//
//        String target = "[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)S[*]";
//        String query = "[H]C([H])(C([O-])=O)C([H])([H])C(=O)S[*]";
////
////        String query = //"Nc1ncnc2n(cnc12)[C@H]1O[C@@H](COP([O-])(=O)OP([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1O";
////                "[H]OP(O)([O-])([O-])OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12";
//////                "Nc1ncnc2n(cnc12)[C@H]1O[C@@H](COP([O-])(=O)OP(O)([O-])=O)[C@H](OP([O-])([O-])=O)[C@@H]1O";
////        String target = "[H]OP(O)([O-])([O-])OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12";
//
//
//        String target = "[H][O][C+]([O][P]([O-])(=[O])[O][C]([H])([H])[C]1([H])[O][C]([H])([n]2[c]([H])[n][c]3[c]([n][c]"
//                + "([H])[n][c]23)[N]([H])[H])[C]([H])([O][H])[C]1([H])[O][H])[c]1[c]([H])[c]([H])[c]([H])[n+]([c]1[H])"
//                + "[C]1([H])[O][C]([H])([C]([H])([H])[O][P]([O-])(=[O])[O][P]([O-])(=[O])[O][C]([H])([H])[C]2([H])[O][C]"
//                + "([H])([n]3[c]([H])[n][c]4[c]([n][c]([H])[n][c]34)[N]([H])[H])[C@]([H])([O][H])[C@]2([H])[O][H])[C@@]"
//                + "([H])([O][H])[C@@]1([H])[O][H]";
//        String query = "[H][O][C@@]1([H])[C]([H])([O][C]([H])([C]([H])([H])[O][P]([O-])(=[O])[O][P]([O-])(=[O])[O][C]"
//                + "([H])([H])[C]2([H])[O][C]([H])([n+]3[c]([H])[c]([H])[c]([H])[c]([c]3[H])[C](=[O])[N]([H])[H])[C@]"
//                + "([H])([O][H])[C@]2([H])[O][H])[C@@]1([H])[O][H])[n]1[c]([H])[n][c]2[c]([n][c]([H])[n][c]12)[N]([H])[H]";
//////
//
////        //This is from MACiE reaction 192.02
//        String query = "[H][N](*)[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]1([H])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](*)=[O])[C]([H])([H])[S][S]*)[C]([H])([H])[S][S][C]1([H])[H])[C]([H])([H])[S-]";
//        String target = "[H][N](*)[C](=[O])[C]1([H])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](=[O])[C]([H])([*])[N]([H])[C](=[O])[C]([H])([N]([H])[C](*)=[O])[C]([H])([H])[S][S]*)[C]([H])([H])[S-])[C]([H])([H])[S][S][C]1([H])[H]";

//        /*Ring matches*/
//        String target = "[H][O][C@@]([H])([C]([H])([H])[H])[C@]([H])([O][H])[C@]1([H])[N]([H])[C]2([O][H])[C]([O-])=[N][C](=[N][C]2=[N][C]1([H])[H])[N]([H])[H]";
//        String query = "[H][O][C@@]([H])([C]([H])([H])[H])[C@]([H])([O][H])[C@]1([H])[N+]([H])=[c]2[c]([O-])[n][c]([n][c]2=[N][C]1([H])[H])[N]([H])[H]";
////

//        String target = "C\\C=C/Nc1cccc(c1)N(O)\\C=C\\C\\C=C\\C=C/C";
//        String query = "Nc1ccccc1";

//        //Subgraph match
//
        String query = "CC";
        String target = "C1CCC12CCCC2";
//        //Single MApping test
////        String target = "O";
////        String query = "Nc1ncnc2n(cnc12)C1OC(COP([O-])=O)C(O)C1O.Nc1ncnc2n(cnc12)C1OC(COP([O-])(=O)OP([O-])(=O)OCC2OC([C@H](O)[C@@H]2O)[n+]2cccc([CH+]O)c2)[C@@H](O)[C@H]1O";
//
//
////        String query = "C1CC1";
//        String target = "CC(C)C";
//
////        String query = "CC";
////        String query="C1CC1";//"C1CC1(CC1CC1)";//Gillian
//        String query="C1CCCCC1";//Hexane
//         String query = "C1=CC=CC=C1";//BENZENE
////        String target = "C1CCC12CCCC2";
//        String target="C1CCCCC1";//Hexane//
//        String target = "C1=CC=CC=C1";//BENZENE
////        String target = "C1CC1(CC1CC1)";//"C1CC1";//gillian

//        String query = "CCCc1nn(C)c2c1nc(nc2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1";//cdk bug expected 9 got 8
//        String target = "NC1=NC2=C(N=CN2[C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)C(=O)N1";//cdk bug expected 9 got 8
////
//        String query = "O=C(CN1C=NC=CC1=O)NCC2=CC=CC=C2";//rguha query
//        String target = "CC(C)(C)N1C2=C(C=N1)C(=O)N(C=N2)CC(=O)NCC3=CC=CC=C3Cl";

//        //Sameer hassan
//        String query = "(c1cc(c[n+](c1)[C@H]2[C@@H]([C@@H]([C@H](O2)CO[P@@](=O)([O-])O[P@@](=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)n4cnc5c4ncnc5N)O)O)O)O)C(=O)N)";//sameer query
//        String target = "(C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C4N=CN=C5N)O)O)O)O)";



        IAtomContainer mol1 = sp.parseSmiles(query);
        IAtomContainer mol2 = sp.parseSmiles(target);

//        CanonicalLabeler canLabler = new CanonicalLabeler();
//        canLabler.canonLabel(molecule);

//
//        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(mol1);

        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(mol1));
        sdg.generateCoordinates();
        mol1 = sdg.getMolecule();

        sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(mol2));
        sdg.generateCoordinates();
        mol2 = sdg.getMolecule();


        boolean bondSensitive = true;
        boolean removeHydrogen = true;
        boolean stereoMatch = false;//true;
        boolean fragmentMinimization = false;//true;
        boolean energyMinimization = false;

//        Isomorphism comparison = new Isomorphism(Algorithm.CDKMCS, bondSensitive);
//        comparison.init(mol1, mol2, removeHydrogen, true);
//        comparison.setChemFilters(stereoMatch, fragmentMinimization, energyMinimization);

        Substructure ss = new Substructure();
        ss.init(mol1, mol2);

//      Get modified Query and Target Molecules as Mappings will correspond to these molecules
//        Assert.assertEquals(true, comparison.isSubgraph());
//        Assert.assertEquals(7, comparison.getFirstAtomMapping().size());
////        Assert.assertEquals(2, comparison.getAllMapping().size());

        Assert.assertEquals(true, ss.isSubgraph(bondSensitive));

//        mol1 = comparison.getReactantMolecule();
//        mol2 = comparison.getProductMolecule();

        System.out.println("Mol1 Size. " + mol1.getAtomCount());
        System.out.println("Mol2 Size. " + mol2.getAtomCount());

//        Double score = 610.0;
//        Assert.assertEquals(score, comparison.getEnergyScore(0));
//        Assert.assertNotNull(comparison.getFirstMapping());
//        Assert.assertEquals(2, comparison.getAllAtomMapping().size());

//        generateImage("TestCase", mol1, mol2, comparison);
//        printMapping(comparison, mol1, mol2);

    }

    private static void printMapping(Isomorphism comparison, IAtomContainer query, IAtomContainer target) {

        int count_final_sol = 0;
        System.out.println("Output of the final Mappings: ");
        try {
            if (comparison.getAllMapping() != null) {

                for (Map<Integer, Integer> final_solution : comparison.getAllMapping()) {
                    int final_solution_size = final_solution.size();
                    System.out.println("Final mapping Nr. " + ++count_final_sol + " Size:" + final_solution_size);

                    for (Map.Entry<Integer, Integer> mapping : final_solution.entrySet()) {
                        System.out.println((mapping.getKey() + 1) + " " + (mapping.getValue() + 1));

                        IAtom eAtom = query.getAtom(mapping.getKey());
                        IAtom pAtom = target.getAtom(mapping.getValue());

                        System.out.println(eAtom.getSymbol() + " "
                                + pAtom.getSymbol());
                    }
                    System.out.println("");

                    System.out.println("Stereo Match: " + comparison.getStereoScore(count_final_sol - 1));
                    System.out.println("Stereo different: " + comparison.isStereoMisMatch());
                    System.out.println("Fragment Size: " + comparison.getFragmentSize(count_final_sol - 1));
                }

                System.out.println("");
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private static void generateImage(String outPutFileName, IAtomContainer query, IAtomContainer target, Isomorphism smsd) throws Exception {

        imageGenerator = new ImageGenerator();

        ////set the format right for the Tanimoto score (only two digits printed)
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(2);
        nf.setMinimumFractionDigits(2);
        System.out.println("Output of the final Mappings: ");
        int counter = 1;
        for (Map<Integer, Integer> mapping : smsd.getAllMapping()) {
            String tanimoto = nf.format(smsd.getTanimotoSimilarity());
            String stereo = "NA";
            if (smsd.getStereoScore(counter - 1) != null) {
                stereo = nf.format(smsd.getStereoScore(counter - 1));
            }
            String label = "Scores [" + "Tanimoto: " + tanimoto + ", Stereo: " + stereo + "]";
            imageGenerator.addImages(query, target, label, mapping);
            counter++;
        }
        String filePNG = System.getProperty("user.dir") + File.separator + outPutFileName;
        imageGenerator.createImage(filePNG, "Query", "Target");

    }

    @Test
    public static void testS() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("[H]C([H])(C([O-])=O)C([H])([H])C(=O)SCCN");
        IAtomContainer queryac = sp.parseSmiles("[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)SCC");
        QueryAtomContainer query = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(queryac, false);

        Isomorphism comparison = new Isomorphism(Algorithm.VFLibMCS, true);
        comparison.init(query, target);
        comparison.setChemFilters(true, true, true);

        System.out.println("[H]C([H])(C([O-])=O)C([H])([H])C(=O)SCCN should be a subgraph of "
                + "[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)SCC " + comparison.isSubgraph());
        System.out.println("[H]C([H])(C([O-])=O)C([H])([H])C(=O)SCCN should be a isomorph of "
                + "[H]C([H])([H])[C@@]([H])(C([O-])=O)C(=O)SCC " + comparison.getTanimotoSimilarity());
    }

    @Test
    public static void testAnyAtomAnyBondCase() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer target = sp.parseSmiles("O1C=CC=C1");
        IAtomContainer queryac = sp.parseSmiles("C1CCCC1");
        QueryAtomContainer query = QueryAtomContainerCreator.createAnyAtomAnyBondContainer(queryac, false);

        Isomorphism comparison = new Isomorphism(Algorithm.VFLibMCS, true);
        comparison.init(query, target);
        comparison.setChemFilters(true, true, true);

        System.out.println("C1CCCC1 should be a subgraph of O1C=CC=C1 " + comparison.isSubgraph());
        System.out.println("C1CCCC1 should be a isomorph of O1C=CC=C1 " + comparison.getTanimotoSimilarity());
    }

    @Test
    public void testQueryTool() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles("CC(=O)OC(=O)C");
        SMARTSQueryTool querytool = new SMARTSQueryTool("O=CO");

        boolean status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        int nmatch = querytool.countMatches();
        Assert.assertEquals(2, nmatch);

        List<Integer> map1 = new ArrayList<Integer>();
        map1.add(1);
        map1.add(2);
        map1.add(3);

        List<Integer> map2 = new ArrayList<Integer>();
        map2.add(3);
        map2.add(4);
        map2.add(5);

        List<List<Integer>> mappings = querytool.getMatchingAtoms();
        List<Integer> ret1 = mappings.get(0);
//        sort(ret1);
        for (int i = 0; i < 3; i++) {
            Assert.assertEquals(map1.get(i), ret1.get(i));
        }

        List<Integer> ret2 = mappings.get(1);
//        sort(ret2);
        for (int i = 0; i < 3; i++) {
            Assert.assertEquals(map2.get(i), ret2.get(i));
        }
    }

    @Test
    public void testQueryToolSingleAtomCase() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles("C1CCC12CCCC2");
        SMARTSQueryTool querytool = new SMARTSQueryTool("C");

        boolean status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        int nmatch = querytool.countMatches();
        Assert.assertEquals(8, nmatch);
    }

    @Test
    public void testQueryToolResetSmarts() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles("C1CCC12CCCC2");
        SMARTSQueryTool querytool = new SMARTSQueryTool("C");

        boolean status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        int nmatch = querytool.countMatches();
        Assert.assertEquals(8, nmatch);

        querytool.setSmarts("CC");
        status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        nmatch = querytool.countMatches();
        Assert.assertEquals(18, nmatch);

        List<List<Integer>> umatch = querytool.getUniqueMatchingAtoms();
        Assert.assertEquals(9, umatch.size());
    }

    @Test
    public void testUniqueQueries() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles("c1ccccc1CCCNCCCc1ccccc1");
        CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        SMARTSQueryTool querytool = new SMARTSQueryTool("c1ccccc1");

        boolean status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        int nmatch = querytool.countMatches();
        Assert.assertEquals(24, nmatch);

        List<List<Integer>> umatch = querytool.getUniqueMatchingAtoms();
        Assert.assertEquals(2, umatch.size());
    }

    @Test
    public void testQuery() throws Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer atomContainer = sp.parseSmiles("c12cc(CCN)ccc1c(COC)ccc2");
        CDKHueckelAromaticityDetector.detectAromaticity(atomContainer);
        SMARTSQueryTool querytool = new SMARTSQueryTool("c12ccccc1cccc2");

        boolean status = querytool.matches(atomContainer);
        Assert.assertTrue(status);

        int nmatch = querytool.countMatches();
        Assert.assertEquals(4, nmatch);

        List<List<Integer>> umatch = querytool.getUniqueMatchingAtoms();
        Assert.assertEquals(1, umatch.size());
    }

    /**
     * Note that we don't test the generated SMILES against the
     * molecule obtained from the factory since the factory derived
     * molecule does not have an explicit hydrogen, which it really should
     * have.
     *
     * @cdk.bug 1985811
     */
    @Test
    public void testIndoleAgainstItself() throws Exception {

        IMolecule indole = MoleculeFactory.makeIndole();

        SmilesGenerator generator = new SmilesGenerator();
        generator.setUseAromaticityFlag(true);
        String indoleSmiles = generator.createSMILES(indole);

        SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        indole = smilesParser.parseSmiles(indoleSmiles);

        SMARTSQueryTool querytool = new SMARTSQueryTool(indoleSmiles);
        Assert.assertTrue(querytool.matches(indole));
    }

    /**
     * @cdk.bug 2149621
     */
    @Test
    public void testMethane() throws Exception {
        IMolecule methane =
                NoNotificationChemObjectBuilder.getInstance().newInstance(IMolecule.class);
        IAtom carbon = methane.getBuilder().newInstance(IAtom.class, Elements.CARBON);
        methane.addAtom(carbon);

        SMARTSQueryTool sqt = new SMARTSQueryTool("CC");
        boolean matches = sqt.matches(methane);
        Assert.assertFalse(matches);

    }

    public static void testQueryAtomContainerDefault() throws CDKException, Exception {
        Isomorphism smsd = new Isomorphism(Algorithm.MCSPlus, true);
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer query = sp.parseSmiles("CC");
        IAtomContainer target = sp.parseSmiles("C1CCC12CCCC2");

//        smsd.init(query, target, false);
//        boolean foundMatches = smsd.isSubgraph();
//        Assert.assertTrue(foundMatches);

        IQueryAtomContainer queryContainer = QueryAtomContainerCreator.createSymbolAndBondOrderQueryContainer(query);
        smsd = new Isomorphism(Algorithm.DEFAULT, true);
        smsd.init(queryContainer, target);
        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();
        generateImage("testQuery ", query, target, smsd);
        boolean foundMatches = smsd.isSubgraph();
        Assert.assertTrue(foundMatches);
    }

    public static void testCOAExampleMols() throws CDKException, Exception {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = sp.parseSmiles("[H][C]([H])([C]([O-])=[O])[C]([H])([H])[C](=[O])[S][*]");
//        IAtomContainer target = sp.parseSmiles("[H][C]([H])([H])[C@@]([H])([C]([O-])=[O])[C](=[O])[S][*]");


//
//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/2-methylhexane.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/hexane.mol");

//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ADP.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/NAD_1a4z.mol");


        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ADP.mol");
        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/NAD_1a4z.mol");


        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(query));
        sdg.generateCoordinates();
        query = sdg.getMolecule();

        sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(target));
        sdg.generateCoordinates();
        target = sdg.getMolecule();


        Isomorphism smsd = new Isomorphism(Algorithm.VFLibMCS, false);
        smsd.init(query, target);
        smsd.setChemFilters(true, true, true);

        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();

        System.out.println("Mol1 Size. " + query.getAtomCount());
        System.out.println("Mol2 Size. " + target.getAtomCount());
        System.out.println("matched Size. " + smsd.getTanimotoSimilarity());
        System.out.println("Matches: " + smsd.getFirstAtomMapping().size());

//        printMapping(smsd, query, target);
        generateImage("testCOA", query, target, smsd);

//        boolean foundMatches = smsd.isSubgraph();
//        Assert.assertTrue(foundMatches);
    }

    public static void testCoverage() throws CDKException, Exception {
//        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
//        IAtomContainer query = sp.parseSmiles("[H][C]([H])([C]([O-])=[O])[C]([H])([H])[C](=[O])[S][*]");
//        IAtomContainer target = sp.parseSmiles("[H][C]([H])([H])[C@@]([H])([C]([O-])=[O])[C](=[O])[S][*]");
//

//
//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/2-methylhexane.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/hexane.mol");

        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ADP.mol");
        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/NAD_1a4z.mol");

//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/M00002.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/M00004.mol");



        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(query));
        sdg.generateCoordinates();
        query = sdg.getMolecule();

        sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(target));
        sdg.generateCoordinates();
        target = sdg.getMolecule();


        Isomorphism smsd = new Isomorphism(Algorithm.VFLibMCS, false);
        smsd.init(query, target);
        smsd.setChemFilters(true, true, true);

        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();

        System.out.println("Mol1 Size. " + query.getAtomCount());
        System.out.println("Mol2 Size. " + target.getAtomCount());
        System.out.println("matched Size. " + smsd.getTanimotoSimilarity());
        System.out.println("Matches: " + smsd.getFirstAtomMapping().size());

//        printMapping(smsd, query, target);
        generateImage("testQueryTargetVFLibMCS", query, target, smsd);


//        IAtomContainer query = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/ADP.mol");
//        IAtomContainer target = readMol("/Users/Asad/Software/GITROOT/SMSD-CMD/Data/NAD_1a4z.mol");


        query = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/M00002.mol");
        target = readMol("/Users/Asad/Software/GITROOT/SMSD-GUI/Data/M00004.mol");



        sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(query));
        sdg.generateCoordinates();
        query = sdg.getMolecule();

        sdg = new StructureDiagramGenerator();
        sdg.setMolecule(new GraphMolecule(target));
        sdg.generateCoordinates();
        target = sdg.getMolecule();


        smsd = new Isomorphism(Algorithm.MCSPlus, false);
        smsd.init(query, target);
        smsd.setChemFilters(true, true, true);

        query = smsd.getReactantMolecule();
        target = smsd.getProductMolecule();

        System.out.println("Mol1 Size. " + query.getAtomCount());
        System.out.println("Mol2 Size. " + target.getAtomCount());
        System.out.println("matched Size. " + smsd.getTanimotoSimilarity());
        System.out.println("Matches: " + smsd.getFirstAtomMapping().size());

//        printMapping(smsd, query, target);
        generateImage("testTargetQueryMCSPLUS", query, target, smsd);
    }

    /**
     * 
     * @param molFile
     * @return
     */
    public static IMolecule readMol(String molFile) {
        IMolecule mol = null;
        StructureDiagramGenerator sdg = new StructureDiagramGenerator();
        try {
            MDLReader reader = new MDLReader(new FileReader(molFile), IChemObjectReader.Mode.RELAXED);
            IMolecule molObject = (IMolecule) reader.read(new GraphMolecule());
            reader.close();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molObject);
            CDKHueckelAromaticityDetector.detectAromaticity(molObject);
            mol=new GraphMolecule(molObject);
            writeMol(mol);
            mol = sdg.getMolecule();
            setID(mol);
        } catch (Exception ex) {
            Logger.getLogger(TestSMSD.class.getName()).log(Level.SEVERE, null, ex);
        }
        return mol;
    }

    private static void writeMol(IMolecule molecule) throws IOException, CDKException {
        MDLV2000Writer writer = new MDLV2000Writer(new FileWriter(new File("output.mol")));
        writer.write((GraphMolecule) molecule);
        writer.close();
    }

    private static void setID(IMolecule mol) {
        int index = 1;
        for (IAtom atom : mol.atoms()) {
            atom.setID(String.valueOf(index++));
        }
    }
}

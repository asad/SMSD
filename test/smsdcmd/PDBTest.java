package smsdcmd;

import org.junit.Test;

public class PDBTest {
    
    public static String MOL_DATA_FOLDER = "test/data/mdl";
    
    public static String PDB_DATA_FOLDER = "test/data/pdb";
    
    @Test
    public void bondOrderTest() {
        ArgumentHandler arguments = new ArgumentHandler();
        arguments.setQueryType("PDB");
        arguments.setTargetType("MOL");
        arguments.setQueryFilepath(PDB_DATA_FOLDER + "/NAD.pdb");
        arguments.setTargetFilepath(MOL_DATA_FOLDER + "/ADP.mol");
        arguments.setMatchBondType(true);
        arguments.setOutputSubgraph(true);
        arguments.setOutputFilepath("--");
        arguments.setOutputFiletype("SMI");
        arguments.setImage(true);
        
        SMSDcmd.run(arguments);
    }

}

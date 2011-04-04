package smsdcmd;

import cmd.ArgumentHandler;
import cmd.SMSDcmd;
import java.io.IOException;
import java.io.StringWriter;

import junit.framework.Assert;

import org.junit.Test;

import cmd.ArgumentHandler;
import cmd.SMSDcmd;

public class PDBTest {
    
    public static String MOL_DATA_FOLDER = "test/data/mdl";
    
    public static String PDB_DATA_FOLDER = "test/data/pdb";
    
    @Test
    public void bondOrderTest() throws IOException {
        ArgumentHandler arguments = new ArgumentHandler();
        arguments.setQueryType("MOL");
        arguments.setTargetType("PDB");
        arguments.setQueryFilepath(MOL_DATA_FOLDER + "/ADP.mol");
        arguments.setTargetFilepath(PDB_DATA_FOLDER + "/NAD.pdb");
        arguments.setMatchBondType(true);
        arguments.setOutputSubgraph(true);
        arguments.setOutputFilepath("--");
        arguments.setOutputFiletype("SMI");
        arguments.setImage(false);
        
        StringWriter writer = new StringWriter();
        arguments.setOutputWriter(writer);
        SMSDcmd.run(arguments);
        writer.close();
        String actual = writer.toString().trim();
        String expected = "NC3=NC=NC1=C3(N=CN1C2OC(COP(=O)(O)OP(=O)(O)O)C(O)C2(O))";
        Assert.assertEquals(expected, actual);
    }

}

package reader.mdl;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Asad
 */
public class SDFReaderTest {

    public static String SDF_DATA_FOLDER = "test/data/mdl";

    @Before
    public void setUp() {
    }

    @Test
    public void testSDFReader() throws FileNotFoundException {
        FileReader in = new FileReader(SDF_DATA_FOLDER + File.separator + "trupti.sdf");
        IteratingMDLReader iteratingMDLReader = new IteratingMDLReader(
                in, NoNotificationChemObjectBuilder.getInstance());
        int i = 0;
        while (iteratingMDLReader.hasNext()) {
            IAtomContainer ac = (IAtomContainer) iteratingMDLReader.next();
            boolean flag = ConnectivityChecker.isConnected(ac);
            if (!flag) {
                System.out.println("error:file not connectted " + ac.getProperty(CDKConstants.TITLE));
            }
            i++;
        }
        System.out.println("total files read: " + i);
    }
}

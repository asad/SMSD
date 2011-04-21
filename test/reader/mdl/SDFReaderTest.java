package reader.mdl;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.libio.md.MDMolecule;
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
    public void testSDFReader() throws FileNotFoundException, CDKException, IOException {
        FileReader in = new FileReader(SDF_DATA_FOLDER + File.separator + "trupti.sdf");
        IteratingMDLReader iteratingMDLReader = new IteratingMDLReader(
                in, NoNotificationChemObjectBuilder.getInstance());
        int i = 0;
        IMoleculeSet connectedSet = NoNotificationChemObjectBuilder.getInstance().newInstance(IMoleculeSet.class);
        IMoleculeSet disconnectedSet = NoNotificationChemObjectBuilder.getInstance().newInstance(IMoleculeSet.class);

        while (iteratingMDLReader.hasNext()) {
            IAtomContainer ac = (IAtomContainer) iteratingMDLReader.next();
            boolean flag = ConnectivityChecker.isConnected(ac);
            IMolecule mol = new MDMolecule(ac);
            mol.setID((String) ac.getProperty(CDKConstants.TITLE));
            mol.setProperty(CDKConstants.TITLE, mol.getID());
            if (!flag) {
                disconnectedSet.addMolecule(mol);
                System.out.println("error:file not connectted " + ac.getProperty(CDKConstants.TITLE));
            } else {
                connectedSet.addMolecule(mol);
            }
            i++;
        }

        SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(SDF_DATA_FOLDER + File.separator + "connected.sdf"));
        sdfWriter.write(connectedSet);
        sdfWriter.close();
        sdfWriter = new SDFWriter(new FileOutputStream(SDF_DATA_FOLDER + File.separator + "disconnected.sdf"));
        sdfWriter.write(disconnectedSet);
        sdfWriter.close();
        System.out.println("total files read: " + i);
    }
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package gui.helper;

import java.io.File;
import javax.swing.filechooser.*;

/**
 *
 * @author sar
 */
public class molFileFilter extends FileFilter {
    //Accept all directories and all gif, jpg, tiff, or png files.

    @Override
    public boolean accept(File f) {
        if (f.isDirectory()) {
            return true;
        }

        String extension = FileFilterUtility.getExtension(f);
        if (extension != null) {
            if (extension.equals(FileFilterUtility.mol)) {
                return true;
            } else if (extension.equals(FileFilterUtility.sdf)) {
                return true;
            } else if (extension.equals(FileFilterUtility.cml)) {
                return true;
            } else {
                return false;
            }
        }

        return false;
    }

    //The description of this filter
    @Override
    public String getDescription() {
        return ".mol,.sdf,.cml";
    }
}
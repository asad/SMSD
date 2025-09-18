/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.gui.helper;

import java.io.File;
import javax.swing.filechooser.*;

/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
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
            switch (extension) {
                case FileFilterUtility.mol:
                    return true;
                case FileFilterUtility.sdf:
                    return true;
                case FileFilterUtility.cml:
                    return true;
                case FileFilterUtility.pdb:
                    return true;
                default:
                    return false;
            }
        }

        return false;
    }

    //The description of this filter
    @Override
    public String getDescription() {
        return ".mol,.sdf,.cml,.pdb";
    }
}

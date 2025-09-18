/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package uk.ac.ebi.smsd.gui.helper;

import java.io.File;
import javax.swing.ImageIcon;

/**
 *
 *  java1.8+
 *
 * 
 * 
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 *
 */
public class FileFilterUtility {

//    public final static String jpeg = "jpeg";
//    public final static String jpg = "jpg";
//    public final static String gif = "gif";
//    public final static String tiff = "tiff";
//    public final static String tif = "tif";
//    public final static String png = "png";
    public final static String mol = "mol";
    public final static String cml = "cml";
    public final static String sdf = "sdf";
    public final static String mol2 = "ml2";
    public final static String pdb = "pdb";

    /*
     * Get the extension of a file.
     */
    public static String getExtension(File f) {
        String ext = null;
        String s = f.getName();
        int i = s.lastIndexOf('.');

        if (i > 0 && i < s.length() - 1) {
            ext = s.substring(i + 1).toLowerCase();
        }
        return ext;
    }

    /**
     * Returns an ImageIcon, or null if the path was invalid.
     *
     * @param path
     * @return ImageIcon
     */
    protected static ImageIcon createImageIcon(String path) {
        java.net.URL imgURL = FileFilterUtility.class.getResource(path);
        if (imgURL != null) {
            return new ImageIcon(imgURL);
        } else {
            System.err.println("Couldn't find file: " + path);
            return null;
        }
    }
}

/* Copyright (C) 2009-2018  Syed Asad Rahman <s9asad@gmail.com>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package uk.ac.ebi.smsd.gui.helper;

import java.io.File;
import java.util.List;
import java.util.ArrayList;

import javax.swing.JFileChooser;

/**
 * An export filter for JCP file formats
 *
 *  jchempaint
 * @author Egon Willighagen
 * 2003-04-01
 */
public class FileExportFilter extends javax.swing.filechooser.FileFilter {

    // only those extensions are given here that are *not* on JCPFileFilter
    public final static String bmp = "bmp";
    public final static String png = "png";
    public final static String jpg = "jpg";
    protected List<String> types;

    public FileExportFilter(String type) {
        super();
        types = new ArrayList<String>();
        types.add(type);
    }

    /**
     * Adds the JCPFileFilter to the JFileChooser object.
     *
     * @param chooser
     */
    public static void addChoosableFileFilters(JFileChooser chooser) {
        chooser.addChoosableFileFilter(new FileExportFilter(FileExportFilter.bmp));
        chooser.addChoosableFileFilter(new FileExportFilter(FileExportFilter.jpg));
        chooser.addChoosableFileFilter(new FileExportFilter(FileExportFilter.png));
    }

    /**
     * The description of this filter.
     */
    public String getDescription() {
        String type = (String) types.get(0);
        String result = "Unknown";
        switch (type) {
            case png:
                result = "PNG";
                break;
            case bmp:
                result = "BMP";
                break;
            case jpg:
                result = "JPEG";
                break;
        }
        return result;
    }

    // Accept all directories and all gif, jpg, or tiff files.
    @Override
    public boolean accept(File f) {
        boolean accepted = false;
        if (f.isDirectory()) {
            accepted = true;
        }

        String extension = getExtension(f);
        if (extension != null) {
            if (types.contains(extension)) {
                accepted = true;
            }
        }
        return accepted;
    }

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

    public String getType() {
        return (String) types.get(0);
    }

    public void setType(String type) {
        types.add(type);
    }
}

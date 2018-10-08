/* Copyright (C) 2009-2018  Syed Asad Rahman <asad@ebi.ac.uk>
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

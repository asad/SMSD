package gui.helper;

import java.io.File;
import java.util.List;
import java.util.ArrayList;

import javax.swing.JFileChooser;

/**
 * An export filter for JCP file formats
 *
 * @cdk.module jchempaint
 * @author  Egon Willighagen
 * @cdk.created 2003-04-01
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
        if (type.equals(png)) {
            result = "PNG";
        } else if (type.equals(bmp)) {
            result = "BMP";
        } else if (type.equals(jpg)) {
            result = "JPEG";
        }
        return result;
    }

    // Accept all directories and all gif, jpg, or tiff files.
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

package gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;
import javax.vecmath.Vector2d;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.reactionblast.graphics.direct.DirectMoleculeDrawer;
import org.openscience.reactionblast.graphics.direct.Params;
import org.openscience.reactionblast.graphics.direct.layout.SingleMoleculeLayout;
import org.openscience.reactionblast.graphics.direct.layout.ZoomToFitGridLayout;

public class ImageGenerator {

    private class QueryTargetPair {

        public final IAtomContainer query;
        public final IAtomContainer target;
        public final IAtomContainer querySubgraph;
        public final IAtomContainer targetSubgraph;
        public final String label;

        public QueryTargetPair(
                IAtomContainer query, IAtomContainer target,
                IAtomContainer querySubgraph, IAtomContainer targetSubgraph,
                String label) {
            this.query = query;
            this.target = target;
            this.querySubgraph = querySubgraph;
            this.targetSubgraph = targetSubgraph;
            this.label = label;
        }
    }
    public final static int SUB_IMAGE_WIDTH = 300;
    public final static int SUB_IMAGE_HEIGHT = 300;
    private List<QueryTargetPair> queryTargetPairs;
    private Params params;

    public ImageGenerator() {
        queryTargetPairs = new ArrayList<QueryTargetPair>();
        params = new Params();
    }

    public void addImages(
            IAtomContainer query,
            IAtomContainer target,
            String label,
            Map<Integer, Integer> maxac) throws IOException, Exception {

        SingleMoleculeLayout msl = new SingleMoleculeLayout(params);
        msl.layout(query, new Vector2d(0.0, 0.0));
        msl.layout(target, new Vector2d(0.0, 0.0));

        IAtomContainer cloneOfQuery = (IAtomContainer) (new org.openscience.cdk.Molecule(query).clone());
        IAtomContainer cloneOfTarget = (IAtomContainer) (new org.openscience.cdk.Molecule(target).clone());

        IAtomContainer querySubgraph = query.getBuilder().newInstance(IAtomContainer.class, cloneOfQuery);
        IAtomContainer targetSubgraph = target.getBuilder().newInstance(IAtomContainer.class, cloneOfTarget);
        List<IAtom> n1 = new ArrayList<IAtom>(query.getAtomCount());
        List<IAtom> n2 = new ArrayList<IAtom>(target.getAtomCount());

        for (Map.Entry<Integer, Integer> aMaps : maxac.entrySet()) {
            IAtom qAtom = cloneOfQuery.getAtom(aMaps.getKey());
            IAtom tAtom = cloneOfTarget.getAtom(aMaps.getValue());
            qAtom.setID(aMaps.getKey().toString());
            tAtom.setID(aMaps.getValue().toString());
            n1.add(qAtom);
            n2.add(tAtom);
        }

        for (IAtom atom : cloneOfQuery.atoms()) {
            if (!n1.contains(atom)) {
                querySubgraph.removeAtomAndConnectedElectronContainers(atom);
            }
        }

        for (IAtom atom : cloneOfTarget.atoms()) {
            if (!n2.contains(atom)) {
                targetSubgraph.removeAtomAndConnectedElectronContainers(atom);
            }
        }

        queryTargetPairs.add(
                new QueryTargetPair(
                cloneOfQuery, cloneOfTarget, querySubgraph, targetSubgraph, label));
    }

    public void createImage(String outImageFileName, String qName, String tName) {

        // layout, and set the highlight subgraphs
        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
        IAtomContainerSet leftHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        IAtomContainerSet rightHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        for (QueryTargetPair pair : queryTargetPairs) {
            moleculeDrawer.addHighlights(pair.querySubgraph);
            moleculeDrawer.addHighlights(pair.targetSubgraph);
            leftHandMoleculeSet.addAtomContainer(pair.query);
            rightHandMoleculeSet.addAtomContainer(pair.target);
        }

        // calculate the total dimensions of the final image
        int width = SUB_IMAGE_WIDTH * 2;
        int height = SUB_IMAGE_HEIGHT * queryTargetPairs.size();

        // make the image, and draw the molecules on it
        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        List<IAtomContainer> mols = new ArrayList<IAtomContainer>();
        for (QueryTargetPair pair : queryTargetPairs) {
            mols.add(pair.query);
            mols.add(pair.target);
        }
        ZoomToFitGridLayout layoutDrawer = new ZoomToFitGridLayout(moleculeDrawer, queryTargetPairs.size(), 2);
        layoutDrawer.layout(mols, new Dimension(SUB_IMAGE_WIDTH, SUB_IMAGE_HEIGHT), g);

        float labelX = SUB_IMAGE_WIDTH / 2;
        float labelY = 15;
        g.setColor(Color.BLACK);
        for (QueryTargetPair pair : queryTargetPairs) {
            g.drawString(pair.label, labelX, labelY);
            labelY += SUB_IMAGE_HEIGHT;
        }
        g.dispose();

        try {
            ImageIO.write((RenderedImage) image, "PNG",
                    new File(outImageFileName + ".png"));
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }

    public RenderedImage createImage() {

        // layout, and set the highlight subgraphs
        DirectMoleculeDrawer moleculeDrawer = new DirectMoleculeDrawer();
        IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
        IAtomContainerSet leftHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        IAtomContainerSet rightHandMoleculeSet = builder.newInstance(IAtomContainerSet.class);
        for (QueryTargetPair pair : queryTargetPairs) {
            moleculeDrawer.addHighlights(pair.querySubgraph);
            moleculeDrawer.addHighlights(pair.targetSubgraph);
            leftHandMoleculeSet.addAtomContainer(pair.query);
            rightHandMoleculeSet.addAtomContainer(pair.target);
        }

        // calculate the total dimensions of the final image
        int width = SUB_IMAGE_WIDTH * 2;
        int height = SUB_IMAGE_HEIGHT * queryTargetPairs.size();

        // make the image, and draw the molecules on it
        Image image = moleculeDrawer.makeBlankImage(width, height);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        List<IAtomContainer> mols = new ArrayList<IAtomContainer>();
        for (QueryTargetPair pair : queryTargetPairs) {
            mols.add(pair.query);
            mols.add(pair.target);
        }
        ZoomToFitGridLayout layoutDrawer = new ZoomToFitGridLayout(moleculeDrawer, queryTargetPairs.size(), 2);
        layoutDrawer.layout(mols, new Dimension(SUB_IMAGE_WIDTH, SUB_IMAGE_HEIGHT), g);

        float labelX = SUB_IMAGE_WIDTH / 2;
        float labelY = 15;
        g.setColor(Color.BLACK);
        for (QueryTargetPair pair : queryTargetPairs) {
            g.drawString(pair.label, labelX, labelY);
            labelY += SUB_IMAGE_HEIGHT;
        }
        g.dispose();

        return (RenderedImage) image;

    }
}

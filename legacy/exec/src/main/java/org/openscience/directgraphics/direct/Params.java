/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct;

import java.awt.Color;

public class Params {

    public BondStrokeCap bondStrokeCap = BondStrokeCap.BUTT;
    public BondStrokeJoin bondStrokeJoin = BondStrokeJoin.MITRE;
    public XAlign leftRightAlignment = XAlign.CENTER;
    public YAlign topBottomAlignment = YAlign.CENTER;
    public int bondLength = 30;
    public int borderX = 20;
    public int borderY = 20;
    public int plusGap = 20;
    public int arrowLength = 30;
    public int arrowGap = 10;
    public int arrowHeadLength = 10;
    public boolean drawBounds = false;
    public boolean drawCarbons = false;
    public boolean drawExplicitHydrogens = true;
    public boolean drawImplicitHydrogens = true;
    public boolean drawTerminalCarbons = true;
    public int atomSymbolFontSize = 10;
    public int plusFontSize = 14;
    public boolean drawMappings = true;
    public int subgraphBoxXBorder = 1;
    public int subgraphBoxYBorder = 2;
    public double doubleBondGap = 2.0;
    public int subscriptHeight = 2;
    public int subscriptTextSize = 9;
    public boolean drawAromaticCircles = true;
    public double ringProportion = 0.75;
    public float bondStrokeWidth = 1.1f;
    public double offsetBondDistanceProportion = 0.75;
    public int filledWedgeWidth = 6;
    public double wiggleLineWidth = 4.0;
    public boolean drawAtomID = false;
    public int atomIDFontSize = 7;
    public double labelYGap = 10.0;
    public int moleculeLabelFontSize = 7;
    public int leftToRightMoleculeLabelFontSize = 9;
    public int topToBottomMoleculeLabelFontSize = 8;
    public boolean drawMoleculeID = true;
    public boolean drawLonePairs = true;
    public double electronRadius = 1.0;
    public double bondMarkLength = 6.0;
    public boolean drawSubgraphBoxes = true;
    public double doubleMarkGap = 1.0;
    public int lonePairSeparation = 4;
    public boolean drawHighlights = true;
    public double highlightRadius = 8.0;
    public Color highlightColor = Color.BLUE;
    public boolean highlightsAbove = true;
    public boolean highlightsBelow = false;
    public float highlightAlpha = 0.15f;
    public float highlightBondStroke = 4.0f;
    public boolean drawSubgraphMappingLines = false;
    public boolean colorSubgraphBoxes = true;
    public boolean drawReactionID = false;
    public boolean layoutLeftToRight = true;
    public boolean highlightSubgraphs = false;
    public boolean drawBondStereoChanges = true;
    public double arrowHeadAngle = 45.0;
    public double circularHighlightBorder = 5.0;
    public boolean useCircularHighlight = false;
    public double circularHighlightMinRadius = 10.0;
    public boolean circularHighlightIsConcentric = true;
    public boolean circularHighlightTransparentFilled = false;
    public boolean useAntialias = true;
    public double tripleBondGap = 2.5;
    public boolean drawRS = false;
    public int chiralSymbolFontSize = 9;
    public float dashedWedgeStroke = 1.0f;
    public double dashedGapFactor = 0.1;
    public double dashedWidthFactor = 0.2;
    public double dashedWedgeWidth = 6.0;
    public int arrowHeadIndent = 5;
    public int arrowBodyWidth = 5;
    public boolean drawFatArrow = false;
    public boolean drawArrowFilled = false;
    public ArrowType arrowType = ArrowType.FORWARD;
    public boolean alignMolecules = false;
    public MoleculeAlignMethod moleculeAlignMethod = MoleculeAlignMethod.MAX_AXIS;
    public boolean circularHighlightShowAtoms = true;
    public boolean drawBondFormedCleavedMarks = true;
    public boolean drawBondOrderChangedMarks = true;
    public boolean drawLabelPanel = false;
    public String labelPanelFont = "ROMAN";
    public int labelPanelFontSize = 14;
    public boolean shouldCrop = true;
    public double labelPanelHeight = 20.0;
    public double labelGap = 10.0;

    public static enum MoleculeAlignMethod {
        MAX_AXIS,
        MIN_AREA;

        private MoleculeAlignMethod() {
        }
    }

    public static enum ArrowType {
        FORWARD,
        BACKWARD,
        BIDIRECTIONAL;

        private ArrowType() {
        }
    }

    public static enum YAlign {
        TOP,
        CENTER,
        BOTTOM;

        private YAlign() {
        }
    }

    public static enum XAlign {
        LEFT,
        CENTER,
        RIGHT;

        private XAlign() {
        }
    }

    public static enum BondStrokeJoin {
        BEVEL,
        MITRE,
        ROUND;

        private BondStrokeJoin() {
        }
    }

    public static enum BondStrokeCap {
        BUTT,
        ROUND,
        SQUARE;

        private BondStrokeCap() {
        }
    }

}

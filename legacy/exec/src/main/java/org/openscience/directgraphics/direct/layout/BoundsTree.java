/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
package org.openscience.directgraphics.direct.layout;

import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class BoundsTree
        implements Iterable<Rectangle2D> {

    private Rectangle2D root = null;
    private String rootLabel;
    private Map<String, Rectangle2D> childMap;

    public BoundsTree(String rootLabel) {
        this.rootLabel = rootLabel;
        this.childMap = new HashMap<>();
    }

    public BoundsTree(String rootLabel, String firstLabel, Rectangle2D firstBox) {
        this(rootLabel);
        this.add(firstLabel, firstBox);
    }

    public /* varargs */ BoundsTree(String rootLabel, BoundsTree... boundsTrees) {
        this(rootLabel);
        for (BoundsTree tree : boundsTrees) {
            this.add(rootLabel, tree);
        }
    }

    public BoundsTree getSubtree(String prefix) {
        BoundsTree subtree = new BoundsTree(this.rootLabel);
        for (String label : this.childMap.keySet()) {
            if (!label.startsWith(prefix)) {
                continue;
            }
            subtree.add(label, this.childMap.get(label));
        }
        return subtree;
    }

    public Rectangle2D getRoot() {
        return this.root;
    }

    public void add(String label, Rectangle2D bounds) {
        boolean isEmpty = bounds.getCenterX() == 0.0 && bounds.getCenterY() == 0.0 && bounds.getWidth() == 0.0 && bounds.getHeight() == 0.0;
        this.childMap.put(label, bounds);
        if (this.root == null && !isEmpty) {
            this.root = new Rectangle2D.Double(bounds.getMinX(), bounds.getMinY(), bounds.getWidth(), bounds.getHeight());
            this.childMap.put(this.rootLabel, this.root);
        } else if (!isEmpty) {
            this.root.add(bounds);
        }
        if (this.root != null) {
            // empty if block
        }
    }

    public void add(String label, Point2D point) {
        Rectangle2D.Double bounds = new Rectangle2D.Double(point.getX(), point.getY(), 0.0, 0.0);
        this.childMap.put(label, bounds);
        if (this.root == null) {
            this.root = new Rectangle2D.Double(bounds.getMinX(), bounds.getMinY(), bounds.getWidth(), bounds.getHeight());
        } else {
            this.root.add(point);
        }
    }

    public Rectangle2D getBounds(List<String> labels) {
        Rectangle2D totalBounds = null;
        for (String label : labels) {
            Rectangle2D bounds = this.get(label);
            if (bounds == null) {
                continue;
            }
            if (totalBounds == null) {
                totalBounds = new Rectangle2D.Double(bounds.getMinX(), bounds.getMinY(), bounds.getWidth(), bounds.getHeight());
                continue;
            }
            totalBounds.add(bounds);
        }
        if (totalBounds == null) {
            return new Rectangle2D.Double(0.0, 0.0, 100.0, 100.0);
        }
        return totalBounds;
    }

    public void setRoot(Rectangle2D root) {
        this.root = root;
    }

    public void add(String prefix, BoundsTree tree) {
        for (String label : tree.getBoundLabels()) {
            this.add(prefix + "_" + label, tree.get(label));
        }
    }

    public List<String> getBoundLabels() {
        return new ArrayList<>(this.childMap.keySet());
    }

    public void shift(double dx, double dy) {
        for (String key : this.childMap.keySet()) {
            Rectangle2D bounds = this.childMap.get(key);
            bounds.setRect(bounds.getMinX() + dx, bounds.getMinY() + dy, bounds.getWidth(), bounds.getHeight());
        }
    }

    public Rectangle2D get(String label) {
        return this.childMap.get(label);
    }

    public double getWidth() {
        if (this.root == null) {
            return 0.0;
        }
        return this.root.getWidth();
    }

    public double getHeight() {
        if (this.root == null) {
            return 0.0;
        }
        return this.root.getHeight();
    }

    @Override
    public Iterator<Rectangle2D> iterator() {
        return this.childMap.values().iterator();
    }

    public BoundsTree transform(AffineTransform transform) {
        BoundsTree transformedTree = new BoundsTree(this.rootLabel);
        for (String key : this.childMap.keySet()) {
            Rectangle2D shape = this.childMap.get(key);
            transformedTree.add(key, transform.createTransformedShape(shape).getBounds2D());
        }
        return transformedTree;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (String key : this.childMap.keySet()) {
            Rectangle2D rect = this.get(key);
            sb.append(key).append("=").append(String.format("[(%2.0f, %2.0f), (%2.0f, %2.0f)]", rect.getMinX(), rect.getMinY(), rect.getMaxX(), rect.getMaxY()));
            sb.append("\n");
        }
        return sb.toString();
    }
}

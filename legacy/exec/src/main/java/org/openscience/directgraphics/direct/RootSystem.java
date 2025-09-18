/*
 * SPDX-License-Identifier: Apache-2.0
 * © 2025 BioInception PVT LTD.
 */
/*
 * Decompiled with CFR 0_114.
 * 
 * Could not load the following classes:
 *  org.openscience.cdk.interfaces.IAtom
 *  org.openscience.cdk.interfaces.IBond
 */
package org.openscience.directgraphics.direct;

import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;

public class RootSystem {

    private final List<IAtom> roots = new ArrayList<>();
    private final List<IAtom> leaves = new ArrayList<>();

    public void addRoot(IAtom root) {
        if (this.roots.contains(root)) {
            return;
        }
        this.roots.add(root);
    }

    public void addRootsFromBond(IBond bond) {
        this.addRoot(bond.getAtom(0));
        this.addRoot(bond.getAtom(1));
    }

    public void addLeaf(IAtom leaf) {
        if (this.leaves.contains(leaf)) {
            return;
        }
        this.leaves.add(leaf);
    }

    public List<IAtom> getRoots() {
        return this.roots;
    }

    public List<IAtom> getLeaves() {
        return this.leaves;
    }

    public RootSystem merge(RootSystem otherRootSystem) {
        RootSystem merged = new RootSystem();
        merged.roots.addAll(this.roots);
        merged.roots.addAll(otherRootSystem.roots);
        merged.leaves.addAll(this.leaves);
        merged.leaves.addAll(otherRootSystem.leaves);
        return merged;
    }

    private void printAtomList(List<IAtom> atoms, StringBuilder sb) {
        sb.append("{");
        for (int index = 0; index < atoms.size(); ++index) {
            IAtom root = atoms.get(index);
            sb.append(root.getID());
            if (index >= atoms.size() - 1) {
                continue;
            }
            sb.append(",");
        }
        sb.append("}");
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Roots ");
        this.printAtomList(this.roots, sb);
        sb.append(" Leaves ");
        this.printAtomList(this.leaves, sb);
        return sb.toString();
    }
}

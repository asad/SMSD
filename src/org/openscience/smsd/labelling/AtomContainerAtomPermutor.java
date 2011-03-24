/* Copyright (C) 2010  Gilleain Torrance <gilleain.torrance@gmail.com>
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

package org.openscience.smsd.labelling;

import java.util.Iterator;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 * @cdk.module smsd
 * @cdk.githash
 */

public class AtomContainerAtomPermutor extends Permutor 
    implements Iterator<IAtomContainer>{

    private IAtomContainer original;
    
    public AtomContainerAtomPermutor(IAtomContainer atomContainer) {
        super(atomContainer.getAtomCount());
        original = atomContainer;
    }
    
    @Override
    public IAtomContainer next() {
        int[] p = super.getNextPermutation();
        return AtomContainerAtomPermutor.permute(p, original);
    }
    
    public static IAtomContainer permute(int[] p, IAtomContainer atomContainer) {
        boolean useA = false;
        if (useA) {
            return permuteA(p, atomContainer);
        } else {
            return permuteB(p, atomContainer);
        }
    }

    private static IAtomContainer permuteA(int[] p, IAtomContainer atomContainer) {
        IAtomContainer permutedContainer = null;
        try {
            permutedContainer = 
                atomContainer.getBuilder().newInstance(IAtomContainer.class);
            for (int i = 0; i < p.length; i++) {
                IAtom atom = atomContainer.getAtom(p[i]);
                permutedContainer.addAtom((IAtom) atom.clone());
            }
            for (IBond bond : atomContainer.bonds()) {
                IBond clonedBond = (IBond) bond.clone();
                clonedBond.setAtoms(new IAtom[clonedBond.getAtomCount()]);
                int i = 0;
                for (IAtom atom : bond.atoms()) {
                    int index = atomContainer.getAtomNumber(atom);
                    IAtom permutedAtom = permutedContainer.getAtom(p[index]);
                    clonedBond.setAtom(permutedAtom, i++);
                }
                permutedContainer.addBond(clonedBond);
            }

        } catch (CloneNotSupportedException cne) {
            //?
            System.out.println(cne);
        }
    
        return permutedContainer;
    }
    
    private static IAtomContainer permuteB(int[] p, IAtomContainer atomContainer) {
        IAtomContainer permutedContainer = null;
        try {
            permutedContainer = (IAtomContainer) atomContainer.clone();
            int n = atomContainer.getAtomCount();
            IAtom[] permutedAtoms = new IAtom[n];
            for (int originalIndex = 0; originalIndex < n; originalIndex++) {
                // get the newly cloned atom 
                IAtom atom = permutedContainer.getAtom(originalIndex);
                
                // permute the index
                int newIndex = p[originalIndex];
                
                // put the atom in the new place
                permutedAtoms[newIndex] = atom;
            }
            permutedContainer.setAtoms(permutedAtoms);
        } catch (CloneNotSupportedException cne) {
            //?
            System.out.println(cne);
        }
        return permutedContainer;
    }
   
    @Override
    public void remove() {
        // can just increase rank....
    }

}

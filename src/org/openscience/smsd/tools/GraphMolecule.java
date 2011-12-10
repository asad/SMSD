/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 1997-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
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
package org.openscience.smsd.tools;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.Serializable;

/**
 * Represents the concept of a chemical molecule, an object composed of 
 * atoms connected by bonds.
 *
 * @cdk.module  data
 * @cdk.githash
 *
 * @author      steinbeck
 * @cdk.created 2000-10-02
 *
 * @cdk.keyword molecule
 */
public class GraphMolecule extends GraphAtomContainer implements Serializable, IAtomContainer, Cloneable {

    /**
     * Determines if a de-serialized object is compatible with this class.
     *
     * This value must only be changed if and only if the new version
     * of this class is incompatible with the old version. See Sun docs
     * for <a href=http://java.sun.com/products/jdk/1.1/docs/guide
     * /serialization/spec/version.doc.html>details</a>.
     */
    private static final long serialVersionUID = 6451193093484831136L;

    /**
     *  Creates an GraphMolecule without Atoms and Bonds.
     */
    public GraphMolecule() {
        super();
    }

    /**
     * Constructs a GraphMolecule with
     * a shallow copy of the atoms and bonds of an GraphAtomContainer.
     *
     * @param   container  An GraphMolecule to copy the atoms and bonds from
     */
    public GraphMolecule(IAtomContainer container) {
        super(container);
    }

    /**
     * Returns a one line string representation of this Atom.
     * Methods is conform RFC #9.
     *
     * @return  The string representation of this Atom
     */
    @Override
    public synchronized String toString() {
        StringBuilder description = new StringBuilder();
        description.append("GraphMolecule(");
        description.append(hashCode());
        if (getID() != null) {
            description.append(", ID=").append(getID());
        }
        description.append(", ").append(super.toString());
        description.append(')');
        return description.toString();
    }

    @Override
    public synchronized Object clone() throws CloneNotSupportedException {
        return super.clone();
    }
}

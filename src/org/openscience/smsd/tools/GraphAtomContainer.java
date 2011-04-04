/*  $RCSfile$
 *  $Author$
 *  $Date$
 *  $Revision$
 *
 *  Copyright (C) 1997-2007  Christoph Steinbeck
 *                      2011 Syed Asad Rahman<asad@ebi.ac.uk>
 *  Contact: cdk-devel@lists.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.smsd.tools;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.SingleElectron;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectChangeEvent;
import org.openscience.cdk.interfaces.IChemObjectListener;
import org.openscience.cdk.interfaces.IElectronContainer;
import org.openscience.cdk.interfaces.ILonePair;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.isomorphism.matchers.IQueryBond;

/**
 *  Base class for all chemical objects that maintain a list of Atoms and
 *  ElectronContainers. <p>
 *
 *  Looping over all Bonds in the GraphAtomContainer is typically done like: <pre>
 * Iterator iter = atomContainer.bonds();
 * while (iter.hasNext()) {
 *   IBond aBond = (IBond) iter.next();
 * }
 *
 *  </pre>
 *
 * @cdk.module data
 * @cdk.githash
 *
 * @author steinbeck, Syed Asad Rahman<asad@ebi.ac.uk>
 * @cdk.created 2000-10-02
 */
public class GraphAtomContainer extends ChemObject
        implements IAtomContainer, IChemObjectListener, Serializable, Cloneable {

    /**
     * Determines if a de-serialized object is compatible with this class.
     *
     * This value must only be changed if and only if the new version
     * of this class is incompatible with the old version. See Sun docs
     * for <a href=http://java.sun.com/products/jdk/1.1/docs/guide
     * /serialization/spec/version.doc.html>details</a>.
     */
    private static final long serialVersionUID = 5678100348445919254L;
    /**
     *  Internal array of atoms.
     */
    protected List<IAtom> atoms;
    /**
     *  Internal array of bonds.
     */
    protected List<IBond> bonds;
    /**
     *  Internal array of lone pairs.
     */
    protected List<ILonePair> lonePairs;
    /**
     *  Internal array of single electrons.
     */
    protected List<ISingleElectron> singleElectrons;
    /**
     * Internal list of atom parities.
     */
    protected List<IStereoElement> stereoElements;
    /**
     * Stores Adjacency Map with List of neighbors
     */
    private Map<IAtom, List<IAtom>> atomAdjacencyMap;
    /**
     * Stores Adjacency Map with List of neighbors
     */
    private Map<IAtom, List<IBond>> bondAdjacencyMap;

    /**
     * Constructs an GraphAtomContainer with a copy of the atoms and electronContainers
     * of another GraphAtomContainer (A shallow copy, i.e., with the same objects as in
     * the original GraphAtomContainer).
     *
     * @param  container  An GraphAtomContainer to copy the atoms and electronContainers from
     */
    public GraphAtomContainer(IAtomContainer container) {
        this.atoms = new ArrayList<IAtom>();
        this.bonds = new ArrayList<IBond>();
        this.lonePairs = new ArrayList<ILonePair>();
        this.singleElectrons = new ArrayList<ISingleElectron>();
        this.atomAdjacencyMap = new HashMap<IAtom, List<IAtom>>();
        this.bondAdjacencyMap = new HashMap<IAtom, List<IBond>>();

        this.stereoElements = new ArrayList<IStereoElement>();

        for (int f = 0; f < container.getAtomCount(); f++) {
            atoms.add(f, container.getAtom(f));
            container.getAtom(f).addListener(this);
            atomAdjacencyMap.put(container.getAtom(f), new ArrayList<IAtom>());
            bondAdjacencyMap.put(container.getAtom(f), new ArrayList<IBond>());
        }
        for (int f = 0; f < container.getBondCount(); f++) {
            bonds.add(f, container.getBond(f));
            container.getBond(f).addListener(this);
        }
        for (int f = 0; f < container.getLonePairCount(); f++) {
            lonePairs.add(f, container.getLonePair(f));
            container.getLonePair(f).addListener(this);
        }
        for (int f = 0; f < container.getSingleElectronCount(); f++) {
            singleElectrons.add(f, container.getSingleElectron(f));
            container.getSingleElectron(f).addListener(this);
        }
        /*sort atoms by atom symbol*/
        Collections.sort(atoms, new Mycompare());
        /*generate adjacency map*/
        generateAdjacencyMap();

        container.setFlags(this.getFlags());
        container.setID(this.getID());
    }

    /**
     *  Constructs an empty GraphAtomContainer that will contain a certain number of
     *  atoms and electronContainers. It will set the starting array lengths to the
     *  defined values, but will not create any Atom or ElectronContainer's.
     *
     */
    public GraphAtomContainer() {
        this.atoms = new ArrayList<IAtom>();
        this.bonds = new ArrayList<IBond>();
        this.lonePairs = new ArrayList<ILonePair>();
        this.singleElectrons = new ArrayList<ISingleElectron>();
        this.stereoElements = new ArrayList<IStereoElement>();
        this.atomAdjacencyMap = new HashMap<IAtom, List<IAtom>>();
        this.bondAdjacencyMap = new HashMap<IAtom, List<IBond>>();
    }

    /** {@inheritDoc} */
    @Override
    public void addStereoElement(IStereoElement element) {
        stereoElements.add(stereoElements.size(), element);
    }

    /** {@inheritDoc} */
    @Override
    public Iterable<IStereoElement> stereoElements() {
        return new Iterable<IStereoElement>() {

            @Override
            public Iterator<IStereoElement> iterator() {
                return stereoElements.iterator();
            }
        };
    }

    /**
     *  Sets the array of atoms of this GraphAtomContainer.
     *
     *@param  atoms  The array of atoms to be assigned to this GraphAtomContainer
     *@see           #getAtom
     */
    @Override
    public void setAtoms(IAtom[] atoms) {
        int counter = this.atoms.size();
        for (int i = 0; i < atoms.length; i++) {
            IAtom atom = atoms[i];
            this.atoms.add(counter++, atom);
            atom.addListener(this);
            atomAdjacencyMap.put(atom, new ArrayList<IAtom>());
            bondAdjacencyMap.put(atom, new ArrayList<IBond>());
        }
        Collections.sort(this.atoms, new Mycompare());
        notifyChanged();
    }

    /**
     * Sets the array of bonds of this GraphAtomContainer.
     *
     * @param  bonds  The array of bonds to be assigned to
     *                             this GraphAtomContainer
     * @see  #getBond
     */
    @Override
    public void setBonds(IBond[] bonds) {
        int count = this.bonds.size();
        for (IBond bond : bonds) {
            this.bonds.add(count++, bond);
            bond.addListener(this);
        }
        generateAdjacencyMap();
    }

    /**
     *  Sets the atom at position <code>number</code> in [0,..].
     *
     *@param  number  The position of the atom to be set.
     *@param  atom    The atom to be stored at position <code>number</code>
     *@see            #getAtom(int)
     */
    @Override
    public void setAtom(int number, IAtom atom) {
        atom.addListener(this);
        atoms.set(number, atom);
        Collections.sort(atoms, new Mycompare());
        atomAdjacencyMap.put(atom, new ArrayList<IAtom>());
        bondAdjacencyMap.put(atom, new ArrayList<IBond>());
        notifyChanged();
    }

    /**
     *  Get the atom at position <code>number</code> in [0,..].
     *
     *@param  number  The position of the atom to be retrieved.
     *@return         The atomAt value
     * @see #setAtom(int, org.openscience.cdk.interfaces.IAtom)
     * @see #setAtoms(org.openscience.cdk.interfaces.IAtom[])
     *
     */
    @Override
    public IAtom getAtom(int number) {
        return atoms.get(number);
    }

    /**
     *  Get the bond at position <code>number</code> in [0,..].
     *
     *@param  number  The position of the bond to be retrieved.
     *@return         The bondAt value
     */
    @Override
    public IBond getBond(int number) {
        return bonds.get(number);
    }

    /**
     *  Get the lone pair at position <code>number</code> in [0,..].
     *
     *@param  number  The position of the LonePair to be retrieved.
     *@return         The lone pair number
     */
    @Override
    public ILonePair getLonePair(int number) {
        return lonePairs.get(number);
    }

    /**
     *  Get the single electron at position <code>number</code> in [0,..].
     *
     *@param  number  The position of the SingleElectron to be retrieved.
     *@return         The single electron number
     */
    @Override
    public ISingleElectron getSingleElectron(int number) {
        return singleElectrons.get(number);
    }

    /**
     *  Returns an Iterable for looping over all atoms in this container.
     *
     *@return    An Iterable with the atoms in this container
     */
    @Override
    public Iterable<IAtom> atoms() {
        return new Iterable<IAtom>() {

            @Override
            public Iterator<IAtom> iterator() {
                return new AtomIterator();
            }
        };
    }

    /**
     * The inner AtomIterator class.
     *
     */
    private class AtomIterator implements Iterator<IAtom> {

        private int pointer = 0;

        @Override
        public boolean hasNext() {
            return pointer < getAtomCount();
        }

        @Override
        public IAtom next() {
            return atoms.get(pointer++);
        }

        @Override
        public void remove() {
            removeAtom(--pointer);
        }
    }

    /**
     *  Returns an Iterable for looping over all bonds in this container.
     *
     *@return    An Iterable with the bonds in this container
     */
    @Override
    public Iterable<IBond> bonds() {
        return new Iterable<IBond>() {

            @Override
            public Iterator<IBond> iterator() {
                return new BondIterator();
            }
        };
    }

    /**
     * The inner BondIterator class.
     *
     */
    private class BondIterator implements Iterator<IBond> {

        private int pointer = 0;

        @Override
        public boolean hasNext() {
            return pointer < getBondCount();
        }

        @Override
        public IBond next() {
            return bonds.get(pointer++);
        }

        @Override
        public void remove() {
            removeBond(--pointer);
        }
    }

    /**
     *  Returns an Iterable for looping over all lone pairs in this container.
     *
     *@return    An Iterable with the lone pairs in this container
     */
    @Override
    public Iterable<ILonePair> lonePairs() {
        return new Iterable<ILonePair>() {

            @Override
            public Iterator<ILonePair> iterator() {
                return new LonePairIterator();
            }
        };
    }

    /**
     * The inner LonePairIterator class.
     *
     */
    private class LonePairIterator implements Iterator<ILonePair> {

        private int pointer = 0;

        @Override
        public boolean hasNext() {
            return pointer < getLonePairCount();
        }

        @Override
        public ILonePair next() {
            return lonePairs.get(pointer++);
        }

        @Override
        public void remove() {
            removeLonePair(--pointer);
        }
    }

    /**
     *  Returns an Iterable for looping over all single electrons in this container.
     *
     *@return    An Iterable with the single electrons in this container
     */
    @Override
    public Iterable<ISingleElectron> singleElectrons() {
        return new Iterable<ISingleElectron>() {

            @Override
            public Iterator<ISingleElectron> iterator() {
                return new SingleElectronIterator();
            }
        };
    }

    /**
     * The inner SingleElectronIterator class.
     *
     */
    private class SingleElectronIterator implements Iterator<ISingleElectron> {

        private int pointer = 0;

        @Override
        public boolean hasNext() {
            return pointer < getSingleElectronCount();
        }

        @Override
        public ISingleElectron next() {
            return singleElectrons.get(pointer++);
        }

        @Override
        public void remove() {
            removeSingleElectron(--pointer);
        }
    }

    /**
     *  Returns an Iterable for looping over all electron containers in this container.
     *
     *@return    An Iterable with the electron containers in this container
     */
    @Override
    public Iterable<IElectronContainer> electronContainers() {
        return new Iterable<IElectronContainer>() {

            @Override
            public Iterator<IElectronContainer> iterator() {
                return new ElectronContainerIterator();
            }
        };
    }

    /**
     * The inner ElectronContainerIterator class.
     *
     */
    private class ElectronContainerIterator implements Iterator<IElectronContainer> {

        private int pointer = 0;

        @Override
        public boolean hasNext() {
            return pointer < (getBondCount() + getLonePairCount() + getSingleElectronCount());
        }

        @Override
        public IElectronContainer next() {
            if (pointer < getBondCount()) {
                return bonds.get(pointer++);
            } else if (pointer < getBondCount() + getLonePairCount()) {
                return lonePairs.get((pointer++) - getBondCount());
            } else if (pointer < getBondCount() + getLonePairCount() + getSingleElectronCount()) {
                return singleElectrons.get((pointer++) - getBondCount() - getLonePairCount());
            }
            return null;
        }

        @Override
        public void remove() {
            if (pointer <= getBondCount()) {
                removeBond(--pointer);
            } else if (pointer <= getBondCount() + getLonePairCount()) {
                removeLonePair((--pointer) - getBondCount());
            } else if (pointer <= getBondCount() + getLonePairCount() + getSingleElectronCount()) {
                removeSingleElectron((--pointer) - getBondCount() - getLonePairCount());
            }
        }
    }

    /**
     *  Returns the atom at position 0 in the container.
     *
     *@return    The atom at position 0 .
     */
    @Override
    public IAtom getFirstAtom() {
        return atoms.get(0);
    }

    /**
     *  Returns the atom at the last position in the container.
     *
     *@return    The atom at the last position
     */
    @Override
    public IAtom getLastAtom() {
        return getAtomCount() > 0 ? (Atom) atoms.get(getAtomCount() - 1) : null;
    }

    /**
     *  Returns the position of a given atom in the atoms array. It returns -1 if
     *  the atom does not exist.
     *
     *@param  atom  The atom to be sought
     *@return       The Position of the atom in the atoms array in [0,..].
     */
    @Override
    public int getAtomNumber(IAtom atom) {
        for (int f = 0; f < atoms.size(); f++) {
            if (atoms.get(f).equals(atom)) {
                return f;
            }
        }
        return -1;
    }

    /**
     *  Returns the position of the bond between two given atoms in the
     *  electronContainers array. It returns -1 if the bond does not exist.
     *
     *@param  atom1  The first atom
     *@param  atom2  The second atom
     *@return        The Position of the bond between a1 and a2 in the
     *               electronContainers array.
     */
    @Override
    public int getBondNumber(IAtom atom1, IAtom atom2) {
        return (getBondNumber(getBond(atom1, atom2)));
    }

    /**
     *  Returns the position of a given bond in the electronContainers array. It
     *  returns -1 if the bond does not exist.
     *
     *@param  bond  The bond to be sought
     *@return       The Position of the bond in the electronContainers array in [0,..].
     */
    @Override
    public int getBondNumber(IBond bond) {
        for (int f = 0; f < getBondCount(); f++) {
            if (bonds.get(f).equals(bond)) {
                return f;
            }
        }
        return -1;
    }

    /**
     *  Returns the position of a given lone pair in the lone pair array. 
     *  It returns -1 if the lone pair does not exist.
     *
     *@param  lonePair  The lone pair to be sought
     *@return       The Position of the lone pair in the array..
     */
    @Override
    public int getLonePairNumber(ILonePair lonePair) {
        for (int f = 0; f < getLonePairCount(); f++) {
            if (lonePairs.get(f).equals(lonePair)) {
                return f;
            }
        }
        return -1;
    }

    /**
     *  Returns the position of a given single electron in the single electron array. 
     *  It returns -1 if the single electron does not exist.
     *
     *@param  singleElectron  The single electron to be sought
     *@return       The Position of the single electron in the array.
     */
    @Override
    public int getSingleElectronNumber(ISingleElectron singleElectron) {
        for (int f = 0; f < getSingleElectronCount(); f++) {
            if (singleElectrons.get(f).equals(singleElectron)) {
                return f;
            }
        }
        return -1;
    }

    /**
     *  Returns the ElectronContainer at position <code>number</code> in the
     *  container.
     *
     * @param  number  The position of the ElectronContainer to be returned.
     * @return         The ElectronContainer at position <code>number</code>.
     */
    @Override
    public IElectronContainer getElectronContainer(int number) {
        if (number < this.getBondCount()) {
            return bonds.get(number);
        }
        number -= this.getBondCount();
        if (number < this.getLonePairCount()) {
            return lonePairs.get(number);
        }
        number -= this.getLonePairCount();
        if (number < this.getSingleElectronCount()) {
            return singleElectrons.get(number);
        }
        return null;
    }

    /**
     * Returns the bond that connects the two given atoms.
     *
     * @param  atom1  The first atom
     * @param  atom2  The second atom
     * @return        The bond that connects the two atoms
     */
    @Override
    public IBond getBond(IAtom atom1, IAtom atom2) {
        Collection<IBond> common = new HashSet<IBond>(this.bondAdjacencyMap.get(atom1));
        common.retainAll(this.bondAdjacencyMap.get(atom2));
        return common.isEmpty() ? null : common.iterator().next();
    }

    /**
     *  Returns the number of Atoms in this Container.
     *
     *@return    The number of Atoms in this Container
     */
    @Override
    public int getAtomCount() {
        return this.atoms.size();
    }

    /**
     *  Returns the number of Bonds in this Container.
     *
     *@return    The number of Bonds in this Container
     */
    @Override
    public int getBondCount() {
        return this.bonds.size();
    }

    /**
     *  Returns the number of LonePairs in this Container.
     *
     *@return    The number of LonePairs in this Container
     */
    @Override
    public int getLonePairCount() {
        return this.lonePairs.size();
    }

    /**
     *  Returns the number of the single electrons in this container,
     *
     *@return       The number of SingleElectron objects of this GraphAtomContainer
     */
    @Override
    public int getSingleElectronCount() {
        return this.singleElectrons.size();
    }

    /**
     * Returns the number of ElectronContainers in this Container.
     *
     * @return    The number of ElectronContainers in this Container
     */
    @Override
    public int getElectronContainerCount() {
        return this.getBondCount() + this.getLonePairCount() + this.getSingleElectronCount();
    }

    /**
     *  Returns an ArrayList of all atoms connected to the given atom.
     *
     *@param  atom  The atom the bond partners are searched of.
     *@return       The ArrayList with the connected atoms
     */
    @Override
    public List<IAtom> getConnectedAtomsList(IAtom atom) {
        return this.atomAdjacencyMap.get(atom);
    }

    /**
     *  Returns an ArrayList of all Bonds connected to the given atom.
     *
     *@param  atom  The atom the connected bonds are searched of
     *@return       The ArrayList with connected atoms
     */
    @Override
    public List<IBond> getConnectedBondsList(IAtom atom) {
        return this.bondAdjacencyMap.get(atom);
    }

    /**
     * Returns the array of lone pairs connected to an atom.
     *
     * @param atom The atom for which to get lone pairs
     * @return The array of LonePairs of this GraphAtomContainer
     * @see #getElectronContainer
     * @see #electronContainers()
     * @see #getBond
     */
    @Override
    public List<ILonePair> getConnectedLonePairsList(IAtom atom) {
        List<ILonePair> lps = new ArrayList<ILonePair>();
        for (int i = 0; i < getLonePairCount(); i++) {
            if (lonePairs.get(i).contains(atom)) {
                lps.add(lonePairs.get(i));
            }
        }
        return lps;
    }

    /**
     *  Returns an array of all SingleElectron connected to the given atom.
     *
     *@param  atom  The atom on which the single electron is located
     *@return       The array of SingleElectron of this GraphAtomContainer
     */
    @Override
    public List<ISingleElectron> getConnectedSingleElectronsList(IAtom atom) {
        List<ISingleElectron> lps = new ArrayList<ISingleElectron>();
        for (int i = 0; i < getSingleElectronCount(); i++) {
            if (singleElectrons.get(i).contains(atom)) {
                lps.add(singleElectrons.get(i));
            }
        }
        return lps;
    }

    /**
     *  Returns an ArrayList of all electronContainers connected to the given atom.
     *
     *@param  atom  The atom the connected electronContainers are searched of
     *@return       The ArrayList with the  connected atoms
     */
    @Override
    public List<IElectronContainer> getConnectedElectronContainersList(IAtom atom) {
        List<IElectronContainer> lps = new ArrayList<IElectronContainer>();
        for (int i = 0; i < getBondCount(); i++) {
            if (bonds.get(i).contains(atom)) {
                lps.add(bonds.get(i));
            }
        }
        for (int i = 0; i < getLonePairCount(); i++) {
            if (lonePairs.get(i).contains(atom)) {
                lps.add(lonePairs.get(i));
            }
        }
        for (int i = 0; i < getSingleElectronCount(); i++) {
            if (singleElectrons.get(i).contains(atom)) {
                lps.add(singleElectrons.get(i));
            }
        }
        return lps;
    }

    /**
     *  Returns the number of atoms connected to the given atom.
     *
     *@param  atom  The atom the number of bond partners are searched of.
     *@return       The the size of connected atoms
     */
    @Override
    public int getConnectedAtomsCount(IAtom atom) {
        return this.atomAdjacencyMap.get(atom).size();
    }

    /**
     *  Returns the number of Bonds for a given Atom.
     *
     *@param  atom  The atom
     *@return       The number of Bonds for this atom
     */
    @Override
    public int getConnectedBondsCount(IAtom atom) {
        return this.bondAdjacencyMap.get(atom).size();
    }

    /**
     *  Returns the number of connected atoms (degree) to the given atom.
     *
     *@param  atomNumber  The atomnumber the degree is searched for
     *@return             The number of connected atoms (degree)
     */
    @Override
    public int getConnectedBondsCount(int atomNumber) {
        return this.bondAdjacencyMap.get(atoms.get(atomNumber)).size();
    }

    /**
     *  Returns the number of LonePairs for a given Atom.
     *
     *@param  atom  The atom
     *@return       The number of LonePairs for this atom
     */
    @Override
    public int getConnectedLonePairsCount(IAtom atom) {
        int count = 0;
        for (int i = 0; i < getLonePairCount(); i++) {
            if (lonePairs.get(i).contains(atom)) {
                ++count;
            }
        }
        return count;
    }

    /**
     *  Returns the sum of the SingleElectron for a given Atom.
     *
     *@param  atom  The atom on which the single electron is located
     *@return       The array of SingleElectron of this GraphAtomContainer
     */
    @Override
    public int getConnectedSingleElectronsCount(IAtom atom) {
        int count = 0;
        for (int i = 0; i < getSingleElectronCount(); i++) {
            if (singleElectrons.get(i).contains(atom)) {
                ++count;
            }
        }
        return count;
    }

    /**
     * Returns the sum of the bond orders for a given Atom.
     *
     * @param  atom  The atom
     * @return       The number of bond orders for this atom
     * 
     * @deprecated   Replaced by <code>AtomContainerManipulator#getBondOrderSum(IAtomContainer, IAtom)</code>
     */
    @Override
    public double getBondOrderSum(IAtom atom) {
        double count = 0;
        for (int i = 0; i < getBondCount(); i++) {
            if (bonds.get(i).contains(atom)) {
                if (bonds.get(i).getOrder().equals(IBond.Order.SINGLE)) {
                    count += 1;
                } else if (bonds.get(i).getOrder().equals(IBond.Order.DOUBLE)) {
                    count += 2;
                } else if (bonds.get(i).getOrder().equals(IBond.Order.TRIPLE)) {
                    count += 3;
                } else if (bonds.get(i).getOrder().equals(IBond.Order.QUADRUPLE)) {
                    count += 4;
                }
            }
        }
        return count;
    }

    /**
     * Returns the maximum bond order that this atom currently has in the context
     * of this GraphAtomContainer.
     *
     * @param  atom  The atom
     * @return       The maximum bond order that this atom currently has
     */
    @Override
    public Order getMaximumBondOrder(IAtom atom) {
        IBond.Order max = IBond.Order.SINGLE;
        List<IBond> l = this.bondAdjacencyMap.get(atom);
        for (IBond bond : l) {
            if (bond.getOrder().ordinal() > max.ordinal()) {
                max = bond.getOrder();
            }
        }
        return max;
    }

    /**
     *  Returns the minimum bond order that this atom currently has in the context
     *  of this GraphAtomContainer.
     *
     *@param  atom  The atom
     *@return       The minimum bond order that this atom currently has
     */
    @Override
    public Order getMinimumBondOrder(IAtom atom) {
        IBond.Order min = IBond.Order.QUADRUPLE;
        List<IBond> l = this.bondAdjacencyMap.get(atom);
        for (IBond bond : l) {
            if (bond.getOrder().ordinal() < min.ordinal()) {
                min = bond.getOrder();
            }
        }
        return min;
    }

    /**
     *  Adds all atoms and electronContainers of a given atomcontainer to this
     *  container.
     *
     *@param  atomContainer  The atomcontainer to be added
     */
    @Override
    public void add(IAtomContainer atomContainer) {
        for (int f = 0; f < atomContainer.getAtomCount(); f++) {
            if (!contains(atomContainer.getAtom(f))) {
                addAtom(atomContainer.getAtom(f));
            }
        }
        for (int f = 0; f < atomContainer.getBondCount(); f++) {
            if (!contains(atomContainer.getBond(f))) {
                addBond(atomContainer.getBond(f));
            }
        }
        for (int f = 0; f < atomContainer.getLonePairCount(); f++) {
            if (!contains(atomContainer.getLonePair(f))) {
                addLonePair(atomContainer.getLonePair(f));
            }
        }
        for (int f = 0; f < atomContainer.getSingleElectronCount(); f++) {
            if (!contains(atomContainer.getSingleElectron(f))) {
                addSingleElectron(atomContainer.getSingleElectron(f));
            }
        }
        notifyChanged();
    }

    /**
     *  Adds an atom to this container.
     *
     *@param  atom  The atom to be added to this container
     */
    @Override
    public void addAtom(IAtom atom) {
        if (contains(atom)) {
            return;
        }
        atom.addListener(this);
        atoms.add(atoms.size(), atom);
        Collections.sort(atoms, new Mycompare());
        atomAdjacencyMap.put(atom, new ArrayList<IAtom>());
        bondAdjacencyMap.put(atom, new ArrayList<IBond>());
        notifyChanged();
    }

    /**
     *  Adds a Bond to this GraphAtomContainer.
     *
     *@param  bond  The bond to added to this container
     */
    @Override
    public void addBond(IBond bond) {
        if (contains(bond)) {
            return;
        }
        bond.addListener(this);
        bonds.add(bonds.size(), bond);
        updateAdjacencyMap(bond);
        notifyChanged();
    }

    /**
     *  Adds a lone pair to this GraphAtomContainer.
     *
     *@param  lonePair  The LonePair to added to this container
     */
    @Override
    public void addLonePair(ILonePair lonePair) {
        lonePair.addListener(this);
        lonePairs.add(lonePairs.size(), lonePair);
        notifyChanged();
    }

    /**
     *  Adds a single electron to this GraphAtomContainer.
     *
     *@param  singleElectron  The SingleElectron to added to this container
     */
    @Override
    public void addSingleElectron(ISingleElectron singleElectron) {
        singleElectron.addListener(this);
        singleElectrons.add(singleElectrons.size(), singleElectron);
        notifyChanged();
    }

    /**
     *  Adds a ElectronContainer to this GraphAtomContainer.
     *
     *@param  electronContainer  The ElectronContainer to added to this container
     */
    @Override
    public void addElectronContainer(IElectronContainer electronContainer) {
        if (electronContainer instanceof IBond) {
            this.addBond((IBond) electronContainer);
        }
        if (electronContainer instanceof ILonePair) {
            this.addLonePair((ILonePair) electronContainer);
        }
        if (electronContainer instanceof ISingleElectron) {
            this.addSingleElectron((ISingleElectron) electronContainer);
        }
    }

    /**
     *  Removes all atoms and electronContainers of a given atomcontainer from this
     *  container.
     *
     *@param  atomContainer  The atomcontainer to be removed
     */
    @Override
    public void remove(IAtomContainer atomContainer) {
        for (int f = 0; f < atomContainer.getAtomCount(); f++) {
            removeAtom(atomContainer.getAtom(f));
        }
        for (int f = 0; f < atomContainer.getBondCount(); f++) {
            removeBond(atomContainer.getBond(f));
        }
        for (int f = 0; f < atomContainer.getLonePairCount(); f++) {
            removeLonePair(atomContainer.getLonePair(f));
        }
        for (int f = 0; f < atomContainer.getSingleElectronCount(); f++) {
            removeSingleElectron(atomContainer.getSingleElectron(f));
        }
    }

    /**
     *  Removes the atom at the given position from the GraphAtomContainer. Note that
     *  the electronContainers are unaffected: you also have to take care of
     *  removing all electronContainers to this atom from the container manually.
     *
     *@param  position  The position of the atom to be removed.
     */
    @Override
    public void removeAtom(int position) {
        atoms.get(position).removeListener(this);
        atoms.remove(position);
        generateAdjacencyMap();
        notifyChanged();
    }

    /**
     *  Removes the given atom from the GraphAtomContainer. Note that the
     *  electronContainers are unaffected: you also have to take care of removing
     *  all electronContainers to this atom from the container.
     *
     *@param  atom  The atom to be removed
     */
    @Override
    public void removeAtom(IAtom atom) {
        int position = getAtomNumber(atom);
        if (position != -1) {
            removeAtom(position);
        }
    }

    /**
     *  Removes the bond at the given position from the GraphAtomContainer.
     *
     *@param  position  The position of the bond to be removed.
     * @return removed bond 
     */
    @Override
    public IBond removeBond(int position) {
        IBond bond = bonds.get(position);
        bond.removeListener(this);
        bonds.remove(bond);
        generateAdjacencyMap();
        notifyChanged();
        return bond;
    }

    /**
     * Removes the bond that connects the two given atoms.
     *
     * @param  atom1  The first atom
     * @param  atom2  The second atom
     * @return        The bond that connects the two atoms
     */
    @Override
    public IBond removeBond(IAtom atom1, IAtom atom2) {
        int pos = getBondNumber(atom1, atom2);
        IBond bond = null;
        if (pos != -1) {
            bond = bonds.get(pos);
            removeBond(pos);
        }
        return bond;
    }

    /**
     * Removes the bond from this container.
     *
     * @param  bond   The bond to be removed.
     */
    @Override
    public void removeBond(IBond bond) {
        int pos = getBondNumber(bond);
        if (pos != -1) {
            removeBond(pos);
        }
    }

    /**
     *  Removes the lone pair at the given position from the GraphAtomContainer.
     *
     *@param  position  The position of the LonePair to be removed.
     */
    @Override
    public ILonePair removeLonePair(int position) {
        ILonePair lp = lonePairs.get(position);
        lp.removeListener(this);
        lonePairs.remove(lp);
        notifyChanged();
        return lp;
    }

    /**
     *  Removes the lone pair from the GraphAtomContainer.
     *
     *@param  lonePair  The LonePair to be removed.
     */
    @Override
    public void removeLonePair(ILonePair lonePair) {
        int pos = getLonePairNumber(lonePair);
        if (pos != -1) {
            removeLonePair(pos);
        }
    }

    /**
     *  Removes the single electron at the given position from the GraphAtomContainer.
     *
     *@param  position  The position of the SingleElectron to be removed.
     */
    @Override
    public ISingleElectron removeSingleElectron(int position) {
        ISingleElectron se = singleElectrons.get(position);
        se.removeListener(this);
        singleElectrons.remove(se);
        notifyChanged();
        return se;
    }

    /**
     *  Removes the single electron from the GraphAtomContainer.
     *
     *@param  singleElectron  The SingleElectron to be removed.
     */
    @Override
    public void removeSingleElectron(ISingleElectron singleElectron) {
        int pos = getSingleElectronNumber(singleElectron);
        if (pos != -1) {
            removeSingleElectron(pos);
        }
    }

    /**
     * Removes the bond at the given position from this container.
     *
     * @param  number  The position of the bond in the electronContainers array
     * @return           Bond that was removed
     */
    @Override
    public IElectronContainer removeElectronContainer(int number) {
        if (number < this.getBondCount()) {
            return removeBond(number);
        }
        number -= this.getBondCount();
        if (number < this.getLonePairCount()) {
            return removeLonePair(number);
        }
        number -= this.getLonePairCount();
        if (number < this.getSingleElectronCount()) {
            return removeSingleElectron(number);
        }
        return null;
    }

    /**
     * Removes this ElectronContainer from this container.
     *
     * @param electronContainer The electronContainer to be removed
     */
    @Override
    public void removeElectronContainer(IElectronContainer electronContainer) {
        if (electronContainer instanceof IBond) {
            removeBond((IBond) electronContainer);
        } else if (electronContainer instanceof ILonePair) {
            removeLonePair((ILonePair) electronContainer);
        } else if (electronContainer instanceof ISingleElectron) {
            removeSingleElectron((ISingleElectron) electronContainer);
        }
    }

    /**
     *  Removes the given atom and all connected electronContainers from the
     *  GraphAtomContainer.
     *
     *@param  atom  The atom to be removed
     */
    @Override
    public void removeAtomAndConnectedElectronContainers(IAtom atom) {
        int position = getAtomNumber(atom);
        if (position != -1) {
            for (int i = 0; i < getBondCount(); i++) {
                if (bonds.get(i).contains(atom)) {
                    removeBond(i);
                    --i;
                }
            }
            for (int i = 0; i < getLonePairCount(); i++) {
                if (lonePairs.get(i).contains(atom)) {
                    removeLonePair(i);
                    --i;
                }
            }
            for (int i = 0; i < getSingleElectronCount(); i++) {
                if (singleElectrons.get(i).contains(atom)) {
                    removeSingleElectron(i);
                    --i;
                }
            }
            removeAtom(position);
        }
        notifyChanged();
    }

    /**
     * Removes all atoms and bond from this container.
     */
    @Override
    public void removeAllElements() {
        removeAllElectronContainers();
        for (int f = 0; f < getAtomCount(); f++) {
            getAtom(f).removeListener(this);
        }
        atoms.clear();
        atoms = new ArrayList<IAtom>();
        clearAdjacencyMap();
        this.atomAdjacencyMap = new HashMap<IAtom, List<IAtom>>();
        this.bondAdjacencyMap = new HashMap<IAtom, List<IBond>>();
        notifyChanged();
    }

    /**
     *  Removes electronContainers from this container.
     */
    @Override
    public void removeAllElectronContainers() {
        removeAllBonds();
        for (int f = 0; f < getLonePairCount(); f++) {
            getLonePair(f).removeListener(this);
        }
        for (int f = 0; f < getSingleElectronCount(); f++) {
            getSingleElectron(f).removeListener(this);
        }
        lonePairs.clear();
        singleElectrons.clear();
        lonePairs = new ArrayList<ILonePair>();
        singleElectrons = new ArrayList<ISingleElectron>();
        notifyChanged();
    }

    /**
     *  Removes all Bonds from this container.
     */
    @Override
    public void removeAllBonds() {
        for (int f = 0; f < getBondCount(); f++) {
            getBond(f).removeListener(this);
        }
        bonds.clear();
        bonds = new ArrayList<IBond>();

        this.atomAdjacencyMap = new HashMap<IAtom, List<IAtom>>();
        this.bondAdjacencyMap = new HashMap<IAtom, List<IBond>>();

        generateAdjacencyMap();
        notifyChanged();
    }

    /**
     *  Adds a bond to this container.
     *
     *@param  atom1   Id of the first atom of the Bond in [0,..]
     *@param  atom2   Id of the second atom of the Bond in [0,..]
     *@param  order   Bondorder
     *@param  stereo  Stereochemical orientation
     */
    @Override
    public void addBond(int atom1, int atom2, IBond.Order order,
            IBond.Stereo stereo) {
        IBond bond = getBuilder().newInstance(IBond.class, getAtom(atom1), getAtom(atom2), order, stereo);
        if (contains(bond)) {
            return;
        }
        addBond(bond);
        /* no notifyChanged() here because addBond(bond) does 
        it already */
    }

    /**
     *  Adds a bond to this container.
     *
     *@param  atom1  Id of the first atom of the Bond in [0,..]
     *@param  atom2  Id of the second atom of the Bond in [0,..]
     *@param  order  Bondorder
     */
    @Override
    public void addBond(int atom1, int atom2, IBond.Order order) {
        IBond bond = getBuilder().newInstance(IBond.class, getAtom(atom1), getAtom(atom2), order);
        addBond(bond);
        /* no notifyChanged() here because addBond(bond) does 
        it already */
    }

    /**
     *  Adds a LonePair to this Atom.
     *
     *@param  atomIndex  The atom number to which the LonePair is added in [0,..]
     */
    @Override
    public void addLonePair(int atomIndex) {
        ILonePair lonePair = getBuilder().newInstance(ILonePair.class, atoms.get(atomIndex));
        lonePair.addListener(this);
        addLonePair(lonePair);
        /* no notifyChanged() here because addElectronContainer() does 
        it already */
    }

    /**
     *  Adds a LonePair to this Atom.
     *
     *@param  atomIndex  The atom number to which the LonePair is added in [0,..]
     */
    @Override
    public void addSingleElectron(int atomIndex) {
        ISingleElectron singleElectron = getBuilder().newInstance(ISingleElectron.class, atoms.get(atomIndex));
        singleElectron.addListener(this);
        addSingleElectron(singleElectron);
        /* no notifyChanged() here because addSingleElectron() does 
        it already */
    }

    /**
     *  True, if the GraphAtomContainer contains the given atom object.
     *
     *@param  atom  the atom this GraphAtomContainer is searched for
     *@return       true if the GraphAtomContainer contains the given atom object
     */
    @Override
    public boolean contains(IAtom atom) {
        for (int i = 0; i < getAtomCount(); i++) {
            if (atom.equals(atoms.get(i))) {
                return true;
            }
        }
        return false;
    }

    /**
     *  True, if the GraphAtomContainer contains the given bond object.
     *
     *@param  bond  the bond this GraphAtomContainer is searched for
     *@return       true if the GraphAtomContainer contains the given bond object
     */
    @Override
    public boolean contains(IBond bond) {
        for (int i = 0; i < getBondCount(); i++) {
            if (bond.equals(bonds.get(i))) {
                return true;
            }
        }
        return false;
    }

    /**
     *  True, if the GraphAtomContainer contains the given LonePair object.
     *
     *@param  lonePair  the LonePair this GraphAtomContainer is searched for
     *@return           true if the GraphAtomContainer contains the given LonePair object
     */
    @Override
    public boolean contains(ILonePair lonePair) {
        for (int i = 0; i < getLonePairCount(); i++) {
            if (lonePair.equals(lonePairs.get(i))) {
                return true;
            }
        }
        return false;
    }

    /**
     *  True, if the GraphAtomContainer contains the given SingleElectron object.
     *
     *@param  singleElectron  the LonePair this GraphAtomContainer is searched for
     *@return           true if the GraphAtomContainer contains the given LonePair object
     */
    @Override
    public boolean contains(ISingleElectron singleElectron) {
        for (int i = 0; i < getSingleElectronCount(); i++) {
            if (singleElectron.equals(singleElectrons.get(i))) {
                return true;
            }
        }
        return false;
    }

    /**
     *  True, if the GraphAtomContainer contains the given ElectronContainer object.
     *
     *@param  electronContainer ElectronContainer that is searched for
     *@return                   true if the GraphAtomContainer contains the given bond object
     */
    @Override
    public boolean contains(IElectronContainer electronContainer) {
        if (electronContainer instanceof IBond) {
            return contains((IBond) electronContainer);
        }
        if (electronContainer instanceof ILonePair) {
            return contains((ILonePair) electronContainer);
        }
        if (electronContainer instanceof ISingleElectron) {
            return contains((SingleElectron) electronContainer);
        }
        return false;
    }

    /**
     *  Returns a one line string representation of this Container. This method is
     *  conform RFC #9.
     *
     *@return    The string representation of this Container
     */
    @Override
    public String toString() {
        StringBuilder stringContent = new StringBuilder(64);
        stringContent.append("GraphAtomContainer(");
        stringContent.append(this.hashCode());
        if (getAtomCount() > 0) {
            stringContent.append(", #A:").append(getAtomCount());
            for (int i = 0; i < getAtomCount(); i++) {
                stringContent.append(", ").append(getAtom(i).toString());
            }
        }
        if (getBondCount() > 0) {
            stringContent.append(", #B:").append(getBondCount());
            for (int i = 0; i < getBondCount(); i++) {
                stringContent.append(", ").append(getBond(i).toString());
            }
        }
        if (getLonePairCount() > 0) {
            stringContent.append(", #LP:").append(getLonePairCount());
            for (int i = 0; i < getLonePairCount(); i++) {
                stringContent.append(", ").append(getLonePair(i).toString());
            }
        }
        if (getSingleElectronCount() > 0) {
            stringContent.append(", #SE:").append(getSingleElectronCount());
            for (int i = 0; i < getSingleElectronCount(); i++) {
                stringContent.append(", ").append(getSingleElectron(i).toString());
            }
        }
        if (stereoElements.size() > 0) {
            stringContent.append(", ST:[#").append(stereoElements.size());
            for (IStereoElement elements : stereoElements) {
                stringContent.append(", ").append(elements.toString());
            }
            stringContent.append(']');
        }
        stringContent.append(')');
        return stringContent.toString();
    }

    /**
     * Clones this GraphAtomContainer object and its content.
     *
     * @return    The cloned object
     * @see       #shallowCopy
     */
    @Override
    public Object clone() throws CloneNotSupportedException {
//
        /**
         * Stores Adjacency Map with List of neighbors
         */
        Map<IAtom, List<IAtom>> cloneAtomAdjacencyMap = new HashMap<IAtom, List<IAtom>>();
        /**
         * Stores Adjacency Map with List of neighbors
         */
        Map<IAtom, List<IBond>> cloneBondAdjacencyMap = new HashMap<IAtom, List<IBond>>();

        List<IBond> clonedBondsT = new ArrayList<IBond>();

        IAtomContainer clone = new GraphAtomContainer();

//        System.out.println("super Atom count " + this.getAtomCount());
        Map<IAtom, IAtom> mapAtomCloneMap = new HashMap<IAtom, IAtom>();

        for (int f = 0; f < this.atoms.size(); f++) {
            IAtom cloneAtm = new Atom(this.atoms.get(f));
            clone.addAtom(cloneAtm);
            cloneAtomAdjacencyMap.put(cloneAtm, new ArrayList<IAtom>());
            cloneBondAdjacencyMap.put(cloneAtm, new ArrayList<IBond>());
            mapAtomCloneMap.put(this.atoms.get(f), cloneAtm);
        }

        for (int f = 0; f < this.bonds.size(); f++) {
            IBond bond = this.bonds.get(f);
            IAtom[] newAtoms = new IAtom[bond.getAtomCount()];
            for (int j = 0; j < bond.getAtomCount(); ++j) {
                newAtoms[j] = mapAtomCloneMap.get(bond.getAtom(j));
            }
            Order bondOrder = bond.getOrder();
            IBond cloneBond = null;
            if (bondOrder != null) {
                cloneBond = new Bond(newAtoms[0], newAtoms[1], bondOrder);
            } else {
                IQueryBond qBond = (IQueryBond) bond;
                if (qBond.getOrder() == null) {
                    cloneBond = new Bond(newAtoms[0], newAtoms[1], Order.SINGLE);
                } else {
                    cloneBond = new Bond(newAtoms[0], newAtoms[1], qBond.getOrder());
                }
            }
            cloneBond.setAtoms(newAtoms);
            clone.addBond(cloneBond);
            clonedBondsT.add(cloneBond);
        }

        ILonePair lp;
        ILonePair newLp;
        for (int i = 0; i < this.lonePairs.size(); ++i) {
            lp = this.lonePairs.get(i);
            newLp = (ILonePair) lp.clone();
            if (lp.getAtom() != null) {
                IAtom atom = lp.getAtom();
                newLp.setAtom(mapAtomCloneMap.get(atom));
            }
            clone.addLonePair(newLp);
        }
        ISingleElectron se;
        ISingleElectron newSe;
        for (int i = 0; i < this.singleElectrons.size(); ++i) {
            se = this.singleElectrons.get(i);
            newSe = (ISingleElectron) se.clone();
            if (se.getAtom() != null) {
                IAtom atom = se.getAtom();
                newSe.setAtom(mapAtomCloneMap.get(atom));
            }
            clone.addSingleElectron(newSe);
        }

        generateAdjacencyMap(cloneAtomAdjacencyMap, cloneBondAdjacencyMap, clonedBondsT);
        clonedBondsT.clear();
        clone.setFlags(this.getFlags());
        clone.setID(this.getID());
        clone.notifyChanged();
        return clone;

    }

    /**
     *  Called by objects to which this object has
     *  registered as a listener.
     *
     *@param  event  A change event pointing to the source of the change
     */
    @Override
    public void stateChanged(IChemObjectChangeEvent event) {
        notifyChanged(event);
    }

    /**
     * Sorts the List by symbols of the elements
     */
    class Mycompare implements Comparator<IAtom> {

        @Override
        public int compare(IAtom o1, IAtom o2) {
            if (o1.getSymbol().equals("H") || o2.getSymbol().equals("H")) {
                return 1;
            }
            return (o1.getSymbol()).compareTo(o2.getSymbol());
        }
    }

    private void generateAdjacencyMap(Map<IAtom, List<IAtom>> newAtomAdjacencyMap, Map<IAtom, List<IBond>> newBondAdjacencyMap, List<IBond> newBonds) {
        for (IBond bond : newBonds) {
            IAtom atom0 = bond.getAtom(0);
            IAtom atom1 = bond.getAtom(1);

            if (newAtomAdjacencyMap.containsKey(atom0)) {
                List<IAtom> l = newAtomAdjacencyMap.get(atom0);
                if (!l.contains(atom1)) {
                    l.add(atom1);
                    newAtomAdjacencyMap.put(atom0, l);
                }
                List<IBond> l1 = newBondAdjacencyMap.get(atom0);
                if (!l1.contains(bond)) {
                    l1.add(bond);
                    newBondAdjacencyMap.put(atom0, l1);
                }
            }

            if (newAtomAdjacencyMap.containsKey(atom1)) {
                List<IAtom> l = newAtomAdjacencyMap.get(atom1);
                if (!l.contains(atom0)) {
                    l.add(atom0);
                    newAtomAdjacencyMap.put(atom1, l);
                }
                List<IBond> l1 = newBondAdjacencyMap.get(atom1);
                if (!l1.contains(bond)) {
                    l1.add(bond);
                    newBondAdjacencyMap.put(atom1, l1);
                }
            }
        }
    }

    private void generateAdjacencyMap() {
        generateAdjacencyMap(this.atomAdjacencyMap, this.bondAdjacencyMap, this.bonds);
    }

    private void clearAdjacencyMap() {
        this.atomAdjacencyMap.clear();
        this.bondAdjacencyMap.clear();
    }

    private void updateAdjacencyMap(IBond bond) {
        IAtom atom0 = bond.getAtom(0);
        IAtom atom1 = bond.getAtom(1);

        if (this.atomAdjacencyMap.containsKey(atom0)) {
            List<IAtom> l = this.atomAdjacencyMap.get(atom0);
            if (!l.contains(atom1)) {
                l.add(atom1);
                this.atomAdjacencyMap.put(atom0, l);
            }
            List<IBond> l1 = this.bondAdjacencyMap.get(atom0);
            if (!l1.contains(bond)) {
                l1.add(bond);
                this.bondAdjacencyMap.put(atom0, l1);
            }
        }

        if (this.atomAdjacencyMap.containsKey(atom1)) {
            List<IAtom> l = this.atomAdjacencyMap.get(atom1);
            if (!l.contains(atom0)) {
                l.add(atom0);
                this.atomAdjacencyMap.put(atom1, l);
            }
            List<IBond> l1 = this.bondAdjacencyMap.get(atom1);
            if (!l1.contains(bond)) {
                l1.add(bond);
                this.bondAdjacencyMap.put(atom1, l1);
            }
        }
    }
}

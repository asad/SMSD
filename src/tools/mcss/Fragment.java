
/* 
 * Copyright (C) 2009-2013 Syed Asad Rahman <asad@ebi.ac.uk>
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
package tools.mcss;

import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.ShortestPathFingerprinter;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 *
 */
public class Fragment implements Comparable<Fragment> {

    public IAtomContainer getContainer() {
        return container;
    }

    public void setContainer(IAtomContainer container) {
        this.container = container;
    }
    private BitSet fingerprint;
    private long fingerprintAsLong;
    private IAtomContainer container;

    public Fragment(IAtomContainer container) throws CDKException {
        if (container == null) {
            throw new CDKException("NULL container not supported");
        }
        ShortestPathFingerprinter spf = new ShortestPathFingerprinter(1024);
        this.container = container;
        this.fingerprint = spf.getBitFingerprint(container).asBitSet();
        this.fingerprintAsLong = convert(this.fingerprint);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Fragment other = (Fragment) obj;
        if (this.fingerprint != other.fingerprint && (this.fingerprint == null || !this.fingerprint.equals(other.fingerprint))) {
            return false;
        }
        if (this.fingerprintAsLong != other.fingerprintAsLong) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 47 * hash + (this.fingerprint != null ? this.fingerprint.hashCode() : 0);
        hash = 47 * hash + (int) (this.fingerprintAsLong ^ (this.fingerprintAsLong >>> 32));
        return hash;
    }

    @Override
    public int compareTo(Fragment t) {

        if (this.fingerprintAsLong == t.fingerprintAsLong) {
            return 0;
        } else if (this.fingerprintAsLong > t.fingerprintAsLong) {
            return 1;
        } else {
            return -1;
        }
    }

    private long convert(BitSet bits) {
        long value = 0L;
        if (bits == null || bits.isEmpty()) {
            return value;
        }
        for (int i = 0; i < bits.length(); ++i) {
            value += bits.get(i) ? (1L << i) : 0L;
        }
        return value;
    }
}

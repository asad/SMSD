
package cmd.pdb;

import java.io.IOException;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.AtomTypeFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.protein.data.PDBAtom;
import org.openscience.cdk.tools.manipulator.BondManipulator;

public class LigandHelper {

    private static AtomTypeFactory atomTypeFactory;
    
    public enum BondOrderMethod { HET_DICT, SATURATE };

    public static void addMissingBondOrders(IAtomContainer ligand) {
        addMissingBondOrders(ligand, BondOrderMethod.HET_DICT);
    }
    
    public static void addMissingBondOrders(IAtomContainer ligand, BondOrderMethod method) {
        if (method == BondOrderMethod.HET_DICT) {
            useBondDictionary(ligand);
        } else {
            saturate(ligand);
        }
    }
    
    private static void useBondDictionary(IAtomContainer ligand) {
        try {
            BondTypeFactory bondFactory = BondTypeFactory.getInstance();
//            String resName = ligand.getID();  // TODO : this should be true...
            String resName = null;
            for (IBond bond : ligand.bonds()) {
                IAtom a0 = bond.getAtom(0);
                String idA;
                if (a0 instanceof PDBAtom) {
                    idA = ((PDBAtom) a0).getName();
                } else {
                    idA = a0.getID();
                }
                
                IAtom a1 = bond.getAtom(1);
                String idB;
                if (a1 instanceof PDBAtom) {
                    idB = ((PDBAtom) a1).getName();
                } else {
                    idB = a1.getID();
                }

                // a hack... atom IDs are 'resName.atomName'
                if (resName == null) {
                    if (a0 instanceof PDBAtom) {
                        resName = ((PDBAtom) a0).getResName();
                    } else {
                        resName = idA.split("\\.")[0];
                    }
                }
                
                if (!(a0 instanceof PDBAtom)) {
                    idA = idA.split("\\.")[1];
                    idB = idB.split("\\.")[1];
                }
                
                IBond.Order order = bondFactory.getBondOrder(resName, idA, idB); 
//                System.out.println("Searching for " + resName + " " + idA + " " + idB + " found " + order);
                bond.setOrder(order);
            }
        } catch (IOException ioe) {
            // TODO
            ioe.printStackTrace();
        }
        
    }

    private static AtomTypeFactory getAtomTypeFactory(IChemObjectBuilder builder) {
        if (atomTypeFactory == null) {
            atomTypeFactory = 
                AtomTypeFactory.getInstance(
//                   "org/openscience/cdk/dict/data/cdk-atom-types.owl", builder);
                     "org/openscience/cdk/config/data/structgen_atomtypes.xml", builder);
        }
        return atomTypeFactory;
    }

    private static void saturate(IAtomContainer atomContainer) {
        AtomTypeFactory typeFactory = 
            getAtomTypeFactory(DefaultChemObjectBuilder.getInstance());
        
        // set some lookup information
        int numberOfAtoms = atomContainer.getAtomCount();
        int[] degreeLookup = new int[numberOfAtoms];
        int[] hCountLookup = new int[numberOfAtoms];
        double[] bosLookup = new double[numberOfAtoms];
        int maxDegree = 0;
        for (int atomIndex = 0; atomIndex < atomContainer.getAtomCount(); atomIndex++) {
            IAtom atom = atomContainer.getAtom(atomIndex);
            int degree = atomContainer.getConnectedAtomsCount(atom);
            degreeLookup[atomIndex] = degree;
            hCountLookup[atomIndex] =
                atom.getImplicitHydrogenCount() == CDKConstants.UNSET ? 0
                    : atom.getImplicitHydrogenCount();
            bosLookup[atomIndex] = atomContainer.getBondOrderSum(atom);
            if (degree > maxDegree) maxDegree = degree;
        }
        
        for (int degree = 1; degree <= maxDegree; degree++) {
            for (int atomIndex = 0; atomIndex < numberOfAtoms; atomIndex++) {
                if (degreeLookup[atomIndex] != degree) continue;
                IAtom atom = atomContainer.getAtom(atomIndex);
                
                IAtomType[] atomTypes1 = typeFactory.getAtomTypes(atom.getSymbol());
                if (atomTypes1.length == 0) continue;

                int hCount = hCountLookup[atomIndex];
                double bosAtom = bosLookup[atomIndex]; 
                double bosType1 = 
                    atomTypes1[0].getBondOrderSum() == CDKConstants.UNSET ? 0.0
                        : atomTypes1[0].getBondOrderSum();
                // check aromatic atoms
                boolean isAromatic = atom.getFlag(CDKConstants.ISAROMATIC); 
                if (isAromatic && bosAtom < bosType1 - hCount) {
                    for (IAtom partner : atomContainer.getConnectedAtomsList(atom)) {
                        int partnerIndex = atomContainer.getAtomNumber(partner);
                        IAtomType[] atomTypes2 = 
                            typeFactory.getAtomTypes(partner.getSymbol());
                        if (atomTypes2.length == 0)
                            return;

                        int partnerHCount = hCountLookup[partnerIndex]; 
                        IBond bond = atomContainer.getBond(partner, atom);
                        double bosPartner = bosLookup[partnerIndex];
                        boolean bondIsAromatic = bond.getFlag(CDKConstants.ISAROMATIC);
                        double bosType2 = 
                            (atomTypes2[0].getBondOrderSum() == CDKConstants.UNSET)? 
                                    0.0 : atomTypes2[0].getBondOrderSum();
                        if (bondIsAromatic && bosPartner < bosType2 - partnerHCount) {
                            bond = atomContainer.getBond(atom, partner);
                            BondManipulator.increaseBondOrder(bond);
                            bosLookup[atomIndex]++;
                            bosLookup[partnerIndex]++;
                            break;
                        }
                    }
                }
                
                // check non-aromatic atoms
                bosAtom = bosLookup[atomIndex];
                if (bosAtom < bosType1 - hCount) {
                    for (IAtom partner : atomContainer.getConnectedAtomsList(atom)) {
                        int partnerIndex = atomContainer.getAtomNumber(partner);
                        IAtomType[] atomTypes2 = 
                            typeFactory.getAtomTypes(partner.getSymbol());
                        if (atomTypes2.length == 0)
                            return;

                        double bos2 = 
                            (atomTypes2[0].getBondOrderSum() == CDKConstants.UNSET)? 
                                    0.0 : atomTypes2[0].getBondOrderSum();
                        int hc2 = hCountLookup[partnerIndex];
                        double acbos2 = bosLookup[partnerIndex];

                        if (acbos2 < bos2 - hc2) {
                            IBond bond = atomContainer.getBond(atom, partner);
                            BondManipulator.increaseBondOrder(bond);
                            bosLookup[atomIndex]++;
                            bosLookup[partnerIndex]++;
                            break;
                        }
                    }
                }
            }
        }
    }
    
}

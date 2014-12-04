/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Test;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.smsd.AtomAtomMapping;
import org.openscience.smsd.Isomorphism;
import org.openscience.smsd.interfaces.Algorithm;
import org.openscience.smsd.tools.ExtAtomContainerManipulator;

/**
 *
 * @author Asad <asad@ebi.ac.uk>
 */
public class UnionTest {

    @Test
    public void unionMolecules() throws IOException, CDKException {
        SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer mol1 = sp.parseSmiles("OOC1=CC=CC=C1");
        IAtomContainer mol2 = sp.parseSmiles("c1ccc(cc1)c2ccccc2");
        int i = 0;
        for (IAtom atom1 : mol1.atoms()) {
            atom1.setID(String.valueOf((i++)));
        }
        int j = 0;
        for (IAtom atom2 : mol2.atoms()) {
            atom2.setID(String.valueOf((j++)));
        }

        ExtAtomContainerManipulator.aromatizeMolecule(mol1);
        ExtAtomContainerManipulator.aromatizeMolecule(mol2);

        Isomorphism isomorphism = new Isomorphism(mol1, mol2, Algorithm.DEFAULT, true, false, false);
        isomorphism.setChemFilters(false, false, false);

        int combinations = 1;

        List<String> acSet = new ArrayList<String>();

        if (isomorphism.getFirstAtomMapping() != null) {

            for (AtomAtomMapping mapping : isomorphism.getAllAtomMapping()) {

                IAtomContainer union = new AtomContainer();

                for (IAtom atom : mol1.atoms()) {
                    union.addAtom(atom);
                }

                for (IBond bond : mol1.bonds()) {
                    union.addBond(bond);
                }

                for (IBond bond : mol2.bonds()) {
                    IAtom a1 = bond.getAtom(0);
                    IAtom a2 = bond.getAtom(1);

                    if (!mapping.getMappingsByAtoms().containsValue(a1)
                            && !mapping.getMappingsByAtoms().containsValue(a2)) {
                        if (!union.contains(a1)) {
                            union.addAtom(a1);
                        }
                        if (!union.contains(a2)) {
                            union.addAtom(a2);
                        }
                        union.addBond(bond);
                    } else if (mapping.getMappingsByAtoms().containsValue(a1)
                            && !mapping.getMappingsByAtoms().containsValue(a2)) {
                        if (!union.contains(a2)) {
                            union.addAtom(a2);
                        }
                        union.addBond(new Bond(a2, getKey(a1, mapping.getMappingsByAtoms()), bond.getOrder(), bond.getStereo()));
                    } else if (!mapping.getMappingsByAtoms().containsValue(a1)
                            && mapping.getMappingsByAtoms().containsValue(a2)) {
                        if (!union.contains(a1)) {
                            union.addAtom(a1);
                        }
                        union.addBond(new Bond(a1, getKey(a2, mapping.getMappingsByAtoms()), bond.getOrder(), bond.getStereo()));
                    }
                }
                /*check if this combination is chemically valid*/
                if (isChemicallyValid(union)) {
                    String molSMILES = getSMILES(union).toString();
                    if (!acSet.contains(molSMILES)) {
                        acSet.add(molSMILES);
                    }
                }

            }
        }

//        for (String container : acSet) {
        ////System.out.println("\n-------------" + " Combination " + combinations++ + "--------------------");
        ////System.out.println("Query SMILES " + getSMILES(mol1).toString() + ", count " + mol1.getAtomCount());
        ////System.out.println("Target SMILES " + getSMILES(mol2).toString() + ", count " + mol2.getAtomCount());
        ////System.out.println("Union SMILES " + container + ", count " + sp.parseSmiles(container).getAtomCount());
//        }
    }

    public IAtom getKey(IAtom a1, Map<IAtom, IAtom> map) {
        for (Map.Entry<IAtom, IAtom> v : map.entrySet()) {
            if (v.getValue() == a1) {
                return v.getKey();
            }
        }
        return null;
    }

    /**
     *
     * @param molecule
     * @return
     * @throws CDKException
     */
    public String getSMILES(IAtomContainer molecule) throws CDKException {

        String smiles = "";
        if (molecule.getAtomCount() == 0) {
            return smiles;
        }

        SmilesGenerator sg = new SmilesGenerator().aromatic();
        try {
            addImplicitHydrogens(molecule);
        } catch (Exception ex) {
            Logger.getLogger(UnionTest.class.getName()).log(Level.SEVERE, null, ex);
        }
        smiles = sg.create(molecule);
        return smiles;
    }

    private boolean isChemicallyValid(IAtomContainer union) throws CDKException {
        for (IAtom atom : union.atoms()) {
            if ((union.getConnectedBondsCount(atom) + atom.getFormalCharge())
                    > atom.getFormalNeighbourCount()) {
                return false;
            }
        }
        return true;
    }

    /**
     * Convenience method that perceives atom types (CDK scheme) and adds
     * implicit hydrogens accordingly. It does not create 2D or 3D coordinates
     * for the new hydrogens.
     *
     * @param container to which implicit hydrogens are added.
     */
    protected void addImplicitHydrogens(IAtomContainer container) throws Exception {
        CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(container.getBuilder());
        int atomCount = container.getAtomCount();
        String[] originalAtomTypeNames = new String[atomCount];
        for (int i = 0; i < atomCount; i++) {
            IAtom atom = container.getAtom(i);
            IAtomType type = matcher.findMatchingAtomType(container, atom);
            originalAtomTypeNames[i] = atom.getAtomTypeName();
            atom.setAtomTypeName(type.getAtomTypeName());
        }
        CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(container.getBuilder());
        hAdder.addImplicitHydrogens(container);
        // reset to the original atom types
        for (int i = 0; i < atomCount; i++) {
            IAtom atom = container.getAtom(i);
            atom.setAtomTypeName(originalAtomTypeNames[i]);
        }
    }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.lang.ref.WeakReference;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

class ChemicalUtilityRegressionTest {
  @ParameterizedTest @ValueSource(strings={"C[CH2]","[OH]","[C]","Cl[Ru](Cl)(Cl)Cl","c1ccccc1[Ru]"})
  void standardisationPreservesSuppliedHydrogenCounts(String smiles)throws Exception {
    IAtomContainer original=rawMolecule(smiles),standard=Standardiser.standardise(original,Standardiser.TautomerMode.NONE);
    for(int atom=0;atom<original.getAtomCount();atom++)
      assertEquals(original.getAtom(atom).getImplicitHydrogenCount(),standard.getAtom(atom).getImplicitHydrogenCount(),"atom="+atom+" "+smiles);
  }
  @ParameterizedTest @CsvSource({"C,0,4","c1ccccc1,0,1","c1cc[nH]c1,3,1","c1ccncc1,3,0","[NH4+],0,4","Cl[Ru](Cl)(Cl)Cl,1,0"})
  void nullHydrogensUseRealBondValence(String smiles,int atom,int expected)throws Exception {
    IAtomContainer original=rawMolecule(smiles);original.getAtom(atom).setImplicitHydrogenCount(null);
    assertEquals(expected,new MolGraph(original).hydrogenCount(atom));
  }
  @Test void programmaticPyrroleHasTheSameFunctionalFingerprint()throws Exception {
    IAtomContainer original=rawMolecule("c1cc[nH]c1"),unset=original.clone();
    for(org.openscience.cdk.interfaces.IAtom atom:unset.atoms())atom.setImplicitHydrogenCount(null);
    assertArrayEquals(FingerprintEngine.fcfp(new MolGraph(original),0,2048),FingerprintEngine.fcfp(new MolGraph(unset),0,2048));
  }
  @Test void aromaticAndKekuleTautomersConserveTheSameMobileProtons()throws Exception {
    MolGraph query=new MolGraph(rawMolecule("Oc1ccccn1"));
    MCS.Options options=new MCS.Options();options.chemistry=ChemOptions.tautomerProfile();
    for(String text:Arrays.asList("O=c1cccc[nH]1","O=C1C=CC=CN1")) {
      MolGraph target=new MolGraph(rawMolecule(text));SearchResult result=MCS.find(query).withOptions(options).match(target);
      assertEquals(7,result.matches().get(0).atomCount());
      assertTrue(MCS.validateTautomerConsistency(query,target,result.matches().get(0).atomMapping()),text);
    }
  }
  @ParameterizedTest @ValueSource(strings={"c1ccccc1","Cc1ccccc1","Cc1ccccc1C","O=C(O)c1ccccc1"})
  void substitutedBenzeneScaffoldsAreClosedShell(String text)throws Exception {
    IAtomContainer input=rawMolecule(text);Cycles.markRingAtomsAndBonds(input);
    IAtomContainer scaffold=MCS.murckoScaffold(input);
    assertEquals(6,scaffold.getAtomCount());int hydrogens=0;
    for(org.openscience.cdk.interfaces.IAtom atom:scaffold.atoms())hydrogens+=atom.getImplicitHydrogenCount();
    assertEquals(6,hydrogens,new MolGraph(scaffold).toCanonicalSmiles());
    assertEquals(new MolGraph(rawMolecule("c1ccccc1")).getCanonicalHash(),new MolGraph(scaffold).getCanonicalHash());
  }
  @Test void editedCdkContainersAreReadByFreshInvocations()throws Exception {
    IAtomContainer query=rawMolecule("CCO"),target=rawMolecule("CCN");Substructure prepared=Substructure.find(query);
    assertFalse(prepared.match(target).found());
    target.getAtom(2).setAtomicNumber(8);target.getAtom(2).setSymbol("O");
    assertTrue(prepared.match(target).found());assertEquals(3,MCS.find(query).match(target).matches().get(0).atomCount());
  }
  private static WeakReference<MolGraph> releasedTarget()throws Exception {
    MolGraph query=new MolGraph(rawMolecule("CC")),target=new MolGraph(rawMolecule("CCC"));
    assertTrue(Substructure.find(query).match(target).found());return new WeakReference<>(target);
  }
  @Test void releasedPreparedQueriesDoNotPinTheirTargets()throws Exception {
    WeakReference<MolGraph> target=releasedTarget();
    for(int attempt=0;attempt<20&&target.get()!=null;attempt++){System.gc();Thread.sleep(10);}
    assertNull(target.get(),"Completed searches must release targets after the prepared query is discarded");
  }
}

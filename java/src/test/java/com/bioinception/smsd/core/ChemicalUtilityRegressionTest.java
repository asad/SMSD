/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.core;

import com.bioinception.smsd.cli.SMSDcli;
import com.bioinception.smsd.cli.SMSDcli.OutputUtil;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.lang.ref.WeakReference;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import picocli.CommandLine;
import static org.junit.jupiter.api.Assertions.*;

class ChemicalUtilityRegressionTest {
  private static IAtomContainer molecule(String smiles)throws Exception {
    return new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles);
  }
  @ParameterizedTest @ValueSource(strings={"C[CH2]","[OH]","[C]","Cl[Ru](Cl)(Cl)Cl","c1ccccc1[Ru]"})
  void standardisationPreservesSuppliedHydrogenCounts(String smiles)throws Exception {
    IAtomContainer original=molecule(smiles),standard=Standardiser.standardise(original,Standardiser.TautomerMode.NONE);
    for(int atom=0;atom<original.getAtomCount();atom++)
      assertEquals(original.getAtom(atom).getImplicitHydrogenCount(),standard.getAtom(atom).getImplicitHydrogenCount(),"atom="+atom+" "+smiles);
  }
  @ParameterizedTest @CsvSource({"C,0,4","c1ccccc1,0,1","c1cc[nH]c1,3,1","c1ccncc1,3,0","[NH4+],0,4","Cl[Ru](Cl)(Cl)Cl,1,0"})
  void nullHydrogensUseRealBondValence(String smiles,int atom,int expected)throws Exception {
    IAtomContainer original=molecule(smiles);original.getAtom(atom).setImplicitHydrogenCount(null);
    assertEquals(expected,new MolGraph(original).hydrogenCount(atom));
  }
  @Test void programmaticPyrroleHasTheSameFunctionalFingerprint()throws Exception {
    IAtomContainer original=molecule("c1cc[nH]c1"),unset=original.clone();
    for(org.openscience.cdk.interfaces.IAtom atom:unset.atoms())atom.setImplicitHydrogenCount(null);
    assertArrayEquals(FingerprintEngine.fcfp(new MolGraph(original),0,2048),FingerprintEngine.fcfp(new MolGraph(unset),0,2048));
  }
  @Test void aromaticAndKekuleTautomersConserveTheSameMobileProtons()throws Exception {
    MolGraph query=new MolGraph(molecule("Oc1ccccn1"));
    MCS.Options options=new MCS.Options();options.chemistry=ChemOptions.tautomerProfile();
    for(String text:Arrays.asList("O=c1cccc[nH]1","O=C1C=CC=CN1")) {
      MolGraph target=new MolGraph(molecule(text));SearchResult result=MCS.find(query).withOptions(options).match(target);
      assertEquals(7,result.matches().get(0).atomCount());
      assertTrue(MCS.validateTautomerConsistency(query,target,result.matches().get(0).atomMapping()),text);
    }
  }
  @ParameterizedTest @ValueSource(strings={"c1ccccc1","Cc1ccccc1","Cc1ccccc1C","O=C(O)c1ccccc1"})
  void substitutedBenzeneScaffoldsAreClosedShell(String text)throws Exception {
    IAtomContainer input=molecule(text);Cycles.markRingAtomsAndBonds(input);
    IAtomContainer scaffold=MCS.murckoScaffold(input);
    assertEquals(6,scaffold.getAtomCount());int hydrogens=0;
    for(org.openscience.cdk.interfaces.IAtom atom:scaffold.atoms())hydrogens+=atom.getImplicitHydrogenCount();
    assertEquals(6,hydrogens,new MolGraph(scaffold).toCanonicalSmiles());
    assertEquals(new MolGraph(molecule("c1ccccc1")).getCanonicalHash(),new MolGraph(scaffold).getCanonicalHash());
  }
  private static final class OwnedStream extends ByteArrayOutputStream {
    boolean closed;
    @Override public void close(){closed=true;}
  }
  @ParameterizedTest @ValueSource(strings={"JSON","SMI","SMARTS","MOL"})
  void mappingWritersLeaveCallerOwnedStreamsOpen(String type)throws Exception {
    IAtomContainer query=molecule("CC"),target=molecule("CCC");OwnedStream output=new OwnedStream();
    Map<Integer,Integer> mapping=new TreeMap<>();mapping.put(0,0);mapping.put(1,1);
    OutputUtil.writeMappings(query,target,Collections.singletonList(mapping),OutputUtil.OutType.valueOf(type),output,false);
    assertFalse(output.closed,type);assertTrue(output.size()>0);
  }
  @ParameterizedTest @ValueSource(strings={"SMI","MOL"})
  void moleculeWritersLeaveCallerOwnedStreamsOpen(String type)throws Exception {
    OwnedStream output=new OwnedStream();OutputUtil.writeMolecule(molecule("CC"),OutputUtil.OutType.valueOf(type),output);
    assertFalse(output.closed,type);assertTrue(output.size()>0);
  }
  @Test void editedCdkContainersAreReadByFreshInvocations()throws Exception {
    IAtomContainer query=molecule("CCO"),target=molecule("CCN");Substructure prepared=Substructure.find(query);
    assertFalse(prepared.match(target).found());
    target.getAtom(2).setAtomicNumber(8);target.getAtom(2).setSymbol("O");
    assertTrue(prepared.match(target).found());assertEquals(3,MCS.find(query).match(target).matches().get(0).atomCount());
  }
  @Test void cliUsesTheCappedCommonFragment()throws Exception {
    PrintStream previous=System.out;ByteArrayOutputStream output=new ByteArrayOutputStream();
    try {
      System.setOut(new PrintStream(output,true,"UTF-8"));
      assertEquals(0,new CommandLine(new SMSDcli()).execute("--Q","SMI","--q","CCO","--T","SMI","--t","CCCO","--mode","mcs","--json","-"));
    } finally {System.setOut(previous);}
    String smiles=new com.fasterxml.jackson.databind.ObjectMapper().readTree(output.toByteArray()).get("mcs_smiles").asText();
    assertEquals(new MolGraph(molecule("CCO")).getCanonicalHash(),new MolGraph(molecule(smiles)).getCanonicalHash());
  }
  private static WeakReference<MolGraph> releasedTarget()throws Exception {
    MolGraph query=new MolGraph(molecule("CC")),target=new MolGraph(molecule("CCC"));
    assertTrue(Substructure.find(query).match(target).found());return new WeakReference<>(target);
  }
  @Test void releasedPreparedQueriesDoNotPinTheirTargets()throws Exception {
    WeakReference<MolGraph> target=releasedTarget();
    for(int attempt=0;attempt<20&&target.get()!=null;attempt++){System.gc();Thread.sleep(10);}
    assertNull(target.get(),"Completed searches must release targets after the prepared query is discarded");
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.cli;

import com.bioinception.smsd.cli.SMSDcli.OutputUtil;
import com.bioinception.smsd.core.MolGraph;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import picocli.CommandLine;
import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

/** CLI writers leave caller-owned streams open, and MCS output reports the capped fragment. */
class CLIOutputRegressionTest {
  private static final class OwnedStream extends ByteArrayOutputStream {
    boolean closed;
    @Override public void close(){closed=true;}
  }
  @ParameterizedTest @ValueSource(strings={"JSON","SMI","SMARTS","MOL"})
  void mappingWritersLeaveCallerOwnedStreamsOpen(String type)throws Exception {
    IAtomContainer query=rawMolecule("CC"),target=rawMolecule("CCC");OwnedStream output=new OwnedStream();
    Map<Integer,Integer> mapping=new TreeMap<>();mapping.put(0,0);mapping.put(1,1);
    OutputUtil.writeMappings(query,target,Collections.singletonList(mapping),OutputUtil.OutType.valueOf(type),output,false);
    assertFalse(output.closed,type);assertTrue(output.size()>0);
  }
  @ParameterizedTest @ValueSource(strings={"SMI","MOL"})
  void moleculeWritersLeaveCallerOwnedStreamsOpen(String type)throws Exception {
    OwnedStream output=new OwnedStream();OutputUtil.writeMolecule(rawMolecule("CC"),OutputUtil.OutType.valueOf(type),output);
    assertFalse(output.closed,type);assertTrue(output.size()>0);
  }
  @Test void cliUsesTheCappedCommonFragment()throws Exception {
    PrintStream previous=System.out;ByteArrayOutputStream output=new ByteArrayOutputStream();
    try {
      System.setOut(new PrintStream(output,true,"UTF-8"));
      assertEquals(0,new CommandLine(new SMSDcli()).execute("--Q","SMI","--q","CCO","--T","SMI","--t","CCCO","--mode","mcs","--json","-"));
    } finally {System.setOut(previous);}
    String smiles=new com.fasterxml.jackson.databind.ObjectMapper().readTree(output.toByteArray()).get("mcs_smiles").asText();
    assertEquals(new MolGraph(rawMolecule("CCO")).getCanonicalHash(),new MolGraph(rawMolecule(smiles)).getCanonicalHash());
  }
}

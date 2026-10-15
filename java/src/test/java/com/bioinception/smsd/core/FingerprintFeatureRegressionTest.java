/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.core;

import java.util.Arrays;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class FingerprintFeatureRegressionTest {
  private static MolGraph graph(String smiles)throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }
  @ParameterizedTest @CsvSource({"C=O,1,2","C#N,1,2","C=N,1,3","N,0,7","[NH4+],0,5","CN(C)C,1,6",
      "c1ccncc1,3,18","c1cc[nH]c1,3,17","[O-],0,10","CO,1,3","[H]O[H],1,3","[H]N([H])[H],1,7"})
  void functionalBitsUseActualAttachedHydrogensAndBondOrders(String smiles,int atom,int expected)throws Exception {
    MolGraph molecule=graph(smiles);
    assertEquals(expected,FingerprintEngine.classifyPharmacophore(molecule,atom),smiles);
    assertEquals(expected,molecule.getPharmacophoreFeatures()[atom],smiles);
  }
  @ParameterizedTest @CsvSource({"[CH3],1650","C,1307","[C],1991","[NH4+],1931","[SiH2],1573","[SiH4],1571"})
  void initialCircularInvariantUsesTheSuppliedHydrogenCount(String smiles,int bit)throws Exception {
    MolGraph molecule=graph(smiles);
    int[] expected=new int[2048];expected[bit]=1;
    assertArrayEquals(expected,FingerprintEngine.ecfpCounts(molecule,0,2048),smiles);
    assertEquals(bitSet(bit),FingerprintEngine.toBitSet(FingerprintEngine.ecfp(molecule,0,2048)),smiles);
  }
  @Test void explicitHydrogenNeighboursContributeToTheAttachedCount()throws Exception {
    MolGraph ammonia=graph("[H]N([H])[H]");int[] expected=new int[2048];
    expected[514]=3;expected[1945]=1;
    assertArrayEquals(expected,FingerprintEngine.ecfpCounts(ammonia,0,2048));
    assertEquals(7,FingerprintEngine.classifyPharmacophore(ammonia,1));
  }
  @Test void storedRadicalHydrogensRemainDistinctFromTheOrdinaryValence()throws Exception {
    assertFalse(Arrays.equals(FingerprintEngine.ecfp(graph("[CH3]"),2,2048),FingerprintEngine.ecfp(graph("C"),2,2048)));
  }
  @ParameterizedTest @CsvSource(delimiter='|',value={
      "CC(=O)C|32,0,2,32|44,446,1548|44,416,417,446,549,585,659,1548,1658|484,1596,1838|484,1009,1057,1134,1343,1375,1596,1780,1838",
      "CC(O)C|32,0,3,32|310,452,1828|310,325,452,600,1080,1763,1828,1900,1987|395,1596,1838|84,210,389,395,694,1343,1596,1780,1838",
      "c1ccncc1|16,16,16,18,16,16|528,790|44,417,454,528,790,995,1310,1641,2027|610,1116|35,217,610,852,1116,1523,1561,1916,1926",
      "CC#N|32,0,2|473,562,1828|473,562,698,745,794,867,988,1697,1828|357,484,1838|11,357,484,1141,1172,1838,1908,1952",
      "C=N|0,3|448,1032|448,726,755,758,947,1032|395,1166|395,555,719,1166,1167,1257",
      "[NH4+]|5|1931|1931|1632|1632",
      "C[NH3+]|0,5|1616,1828|509,774,1616,1750,1828,1877|393,1166|266,368,393,598,1166,1838",
      "CO|0,3|452,1828|452,1241,1305,1558,1828,1942|395,1166|30,395,700,1166,1238,1750",
      "[H]OC|0,3,0|514,1828,1974|67,167,289,514,1035,1270,1433,1828,1974|16,1166|16,1101,1166,1645,1803,1994",
      "c1cc[nH]c1|16,16,16,17,16|528,1749|44,59,528,746,1363,1376,1749,1987|610,961|15,465,610,961,1468,1779,1900,1926",
      "[CH2]|32|1309|1309|599|599",
      "C|32|1307|1307|599|599"})
  void circularAndFunctionalBitsMatchTheSharedDefinitions(String smiles,String features,
      String ecfp0,String ecfp2,String fcfp0,String fcfp2)throws Exception {
    MolGraph molecule=graph(smiles);
    assertArrayEquals(numbers(features),molecule.getPharmacophoreFeatures(),smiles);
    assertArrayEquals(numbers(ecfp0),FingerprintEngine.toBitSet(FingerprintEngine.ecfp(molecule,0,2048)).stream().toArray(),smiles);
    assertArrayEquals(numbers(ecfp2),FingerprintEngine.toBitSet(FingerprintEngine.ecfp(molecule,2,2048)).stream().toArray(),smiles);
    assertArrayEquals(numbers(fcfp0),FingerprintEngine.toBitSet(FingerprintEngine.fcfp(molecule,0,2048)).stream().toArray(),smiles);
    assertArrayEquals(numbers(fcfp2),FingerprintEngine.toBitSet(FingerprintEngine.fcfp(molecule,2,2048)).stream().toArray(),smiles);
  }
  private static java.util.BitSet bitSet(int bit){java.util.BitSet result=new java.util.BitSet();result.set(bit);return result;}
  private static int[] numbers(String text) {
    return Arrays.stream(text.split(",")).mapToInt(Integer::parseInt).toArray();
  }

}

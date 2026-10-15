/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * See the NOTICE file for attribution.
 */
package com.bioinception.smsd.core;

import java.util.*;
import java.util.concurrent.*;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import static org.junit.jupiter.api.Assertions.*;

class TautomerConfidenceRegressionTest {
  private static MolGraph graph(String smiles)throws Exception {
    return new MolGraph(new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles(smiles));
  }
  private static SearchResult.Match match(MolGraph query,MolGraph target) {
    MCS.Options options=new MCS.Options();options.chemistry=ChemOptions.tautomerProfile();
    SearchResult result=MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal());assertEquals(4,result.matches().get(0).atomCount());
    return result.matches().get(0);
  }
  @Test void selectedKetoEnolRelaxationsUseReadTimeContext()throws Exception {
    MolGraph query=graph("CC(=O)C"),target=graph("CC(O)=C");SearchResult.Match mapping=match(query,target);
    float[] queryWeights=query.tautomerWeight.clone(),targetWeights=target.tautomerWeight.clone();
    ChemOptions aqueous=ChemOptions.tautomerProfile();
    ChemOptions solvent=ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.CHLOROFORM);
    ChemOptions acidic=ChemOptions.tautomerProfile().withPH(2);
    double baseline=MCS.computeTautomerConfidence(query,target,mapping,aqueous);
    double solventConfidence=MCS.computeTautomerConfidence(query,target,mapping,solvent);
    double acidicConfidence=MCS.computeTautomerConfidence(query,target,mapping,acidic);
    assertTrue(baseline>0&&baseline<1);
    assertTrue(solventConfidence<baseline);assertTrue(acidicConfidence<baseline);
    for(int repeat=0;repeat<4;repeat++) {
      assertEquals(solventConfidence,MCS.computeTautomerConfidence(query,target,mapping,solvent));
      assertEquals(baseline,MCS.computeTautomerConfidence(query,target,mapping,aqueous));
    }
    assertArrayEquals(queryWeights,query.tautomerWeight);assertArrayEquals(targetWeights,target.tautomerWeight);
    assertEquals(1.0,MCS.computeTautomerConfidence(query,target,mapping.atomMapping()),"The three-argument utility retains its element-only contract");
  }
  @Test void identicalSelectedChemistryHasNoPenaltyInAnyContext()throws Exception {
    MolGraph molecule=graph("CC(=O)C");SearchResult.Match identity=match(molecule,molecule);
    for(ChemOptions.Solvent solvent:ChemOptions.Solvent.values())
      assertEquals(1.0,MCS.computeTautomerConfidence(molecule,molecule,identity,ChemOptions.tautomerProfile().withSolvent(solvent).withPH(0)));
  }
  @Test void unselectedBondOrderDifferencesDoNotAffectConfidence()throws Exception {
    MolGraph query=graph("CC(=O)C"),target=graph("CC(O)=C");SearchResult.Match mapping=match(query,target);
    Map<SearchResult.Bond,SearchResult.Bond> sameBonds=new LinkedHashMap<>();
    for(Map.Entry<SearchResult.Bond,SearchResult.Bond> pair:mapping.bondMapping().entrySet()) {
      SearchResult.Bond first=pair.getKey(),second=pair.getValue();
      if(query.bondOrder(first.first(),first.second())==target.bondOrder(second.first(),second.second()))sameBonds.put(first,second);
    }
    SearchResult.Match selected=new SearchResult.Match(mapping.atomMapping(),sameBonds,mapping.fragments(),mapping.score());
    assertEquals(1.0,MCS.computeTautomerConfidence(query,target,selected,ChemOptions.tautomerProfile()));
  }
  @Test void concurrentContextReadsLeaveSearchAndWeightsUnchanged()throws Exception {
    MolGraph query=graph("CC(=O)C"),target=graph("CC(O)=C");SearchResult.Match mapping=match(query,target);
    float[] queryWeights=query.tautomerWeight.clone(),targetWeights=target.tautomerWeight.clone();
    ChemOptions options=ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.DMSO).withPH(6);
    double expected=MCS.computeTautomerConfidence(query,target,mapping,options);
    ExecutorService workers=Executors.newFixedThreadPool(4);
    try {
      List<Callable<Double>> calls=new ArrayList<>();
      for(int repeat=0;repeat<8;repeat++)calls.add(()->MCS.computeTautomerConfidence(query,target,mapping,options));
      for(Future<Double> value:workers.invokeAll(calls))assertEquals(expected,value.get());
    } finally {workers.shutdownNow();}
    MCS.Options searchOptions=new MCS.Options();searchOptions.chemistry=options;
    SearchResult contextual=MCS.find(query).withOptions(searchOptions).match(target);
    assertEquals(mapping.atomMapping(),contextual.matches().get(0).atomMapping());
    assertEquals(mapping.bondMapping(),contextual.matches().get(0).bondMapping());
    assertArrayEquals(queryWeights,query.tautomerWeight);assertArrayEquals(targetWeights,target.tautomerWeight);
  }
  @Test void malformedSelectedBondsAreRejected()throws Exception {
    MolGraph query=graph("CCO"),target=graph("CCO");
    Map<Integer,Integer> atoms=new LinkedHashMap<>();
    for(int atom=0;atom<3;atom++)atoms.put(atom,atom);
    for(Map<SearchResult.Bond,SearchResult.Bond> bonds:Arrays.asList(
        Collections.singletonMap(new SearchResult.Bond(0,1),new SearchResult.Bond(1,2)),
        Collections.singletonMap(new SearchResult.Bond(0,2),new SearchResult.Bond(0,2)),
        Collections.singletonMap(new SearchResult.Bond(0,1),new SearchResult.Bond(0,3)))) {
      SearchResult.Match malformed=new SearchResult.Match(atoms,bonds,Collections.emptyList(),3);
      assertThrows(IllegalArgumentException.class,
          ()->MCS.computeTautomerConfidence(query,target,malformed,ChemOptions.tautomerProfile()));
    }
  }
  @ParameterizedTest @ValueSource(doubles={Double.NaN,Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY})
  void nonFinitePHIsRejected(double pH)throws Exception {
    assertThrows(IllegalArgumentException.class,()->new ChemOptions().withPH(pH));
    MolGraph molecule=graph("CC(=O)C");SearchResult.Match mapping=match(molecule,molecule);
    ChemOptions options=ChemOptions.tautomerProfile();options.pH=pH;
    assertThrows(IllegalArgumentException.class,()->MCS.computeTautomerConfidence(molecule,molecule,mapping,options));
  }
}

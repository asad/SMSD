/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
package com.bioinception.smsd.core;

import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.*;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;
import static com.bioinception.smsd.TestSupport.rawGraph;
import static com.bioinception.smsd.TestSupport.rawMolecule;
import static org.junit.jupiter.api.Assertions.*;

class TautomerConfidenceRegressionTest {
  private static SearchResult.Match match(MolGraph query,MolGraph target) {
    MCS.Options options=new MCS.Options();options.chemistry=ChemOptions.tautomerProfile();
    SearchResult result=MCS.find(query).withOptions(options).match(target);
    assertTrue(result.optimal());assertEquals(4,result.matches().get(0).atomCount());
    return result.matches().get(0);
  }
  @Test void selectedKetoEnolRelaxationsUseReadTimeContext()throws Exception {
    MolGraph query=rawGraph("CC(=O)C"),target=rawGraph("CC(O)=C");SearchResult.Match mapping=match(query,target);
    float[] queryWeights=query.tautomerWeight.clone(),targetWeights=target.tautomerWeight.clone();
    ChemOptions aqueous=ChemOptions.tautomerProfile();
    ChemOptions solvent=ChemOptions.tautomerProfile().withSolvent(ChemOptions.Solvent.CHLOROFORM);
    ChemOptions acidic=ChemOptions.tautomerProfile().withPH(2);
    double baseline=MCS.computeTautomerConfidence(query,target,mapping,aqueous);
    double solventConfidence=MCS.computeTautomerConfidence(query,target,mapping,solvent);
    double acidicConfidence=MCS.computeTautomerConfidence(query,target,mapping,acidic);
    assertTrue(baseline>0&&baseline<1);
    assertTrue(solventConfidence<baseline);
    assertEquals(baseline,acidicConfidence,"A keto/enol site of pKa 20 is not ionised at pH 2, so its confidence does not move");
    for(int repeat=0;repeat<4;repeat++) {
      assertEquals(solventConfidence,MCS.computeTautomerConfidence(query,target,mapping,solvent));
      assertEquals(baseline,MCS.computeTautomerConfidence(query,target,mapping,aqueous));
    }
    assertArrayEquals(queryWeights,query.tautomerWeight);assertArrayEquals(targetWeights,target.tautomerWeight);
    assertEquals(baseline,MCS.computeTautomerConfidence(query,target,mapping.atomMapping()),"The deprecated three-argument form equals the four-argument form in water at pH 7.4");
  }
  private static double identityConfidence(String first,String second,double pH)throws Exception {
    MolGraph query=rawGraph(first),target=rawGraph(second);
    Map<Integer,Integer> atoms=new LinkedHashMap<>();
    Map<SearchResult.Bond,SearchResult.Bond> bonds=new LinkedHashMap<>();
    for(int atom=0;atom<query.n;atom++) {
      atoms.put(atom,atom);
      for(int neighbor:query.neighbors[atom])
        if(atom<neighbor&&target.hasBond(atom,neighbor))bonds.put(new SearchResult.Bond(atom,neighbor),new SearchResult.Bond(atom,neighbor));
    }
    SearchResult.Match identity=new SearchResult.Match(atoms,bonds,Collections.emptyList(),atoms.size());
    return MCS.computeTautomerConfidence(query,target,identity,ChemOptions.tautomerProfile().withPH(pH));
  }
  @Test void ionisableSitesFollowHendersonHasselbalch()throws Exception {
    double acid2=identityConfidence("ONC(C)=O","ON=C(C)O",2),acid74=identityConfidence("ONC(C)=O","ON=C(C)O",7.4);
    double acid12=identityConfidence("ONC(C)=O","ON=C(C)O",12);
    assertTrue(acid2>acid74&&acid74>acid12,"An acid site loses its neutral form as the pH rises");
    double hydroxamicPKa=MolGraph.TPKA_HYDROXAMIC;
    assertEquals(0.85/(1+Math.pow(10,2-hydroxamicPKa)),acid2,1e-6);
    assertEquals(0.85/(1+Math.pow(10,7.4-hydroxamicPKa)),acid74,1e-6);
    assertEquals(0.1,acid12,1e-6,"The weight stops at its 0.1 floor");
    // 4(3H)-Pyrimidinone holds the N-H site of the pyridinone motif, of pKa 8.6, alone: no aromatic
    // carbon path reaches a para carbon from its carbonyl carbon.
    double pyrimidinone2=identityConfidence("O=c1ccnc[nH]1","Oc1ccncn1",2);
    double pyrimidinone74=identityConfidence("O=c1ccnc[nH]1","Oc1ccncn1",7.4);
    double pyrimidinone12=identityConfidence("O=c1ccnc[nH]1","Oc1ccncn1",12);
    assertTrue(pyrimidinone2>pyrimidinone74&&pyrimidinone74>pyrimidinone12,"An acid site loses its neutral form as the pH rises");
    assertEquals(0.95/(1+Math.pow(10,7.4-8.6)),pyrimidinone74,1e-6);
    assertEquals(0.1,pyrimidinone12,1e-6);
    // 2-Pyridone keeps the N-H pKa of 11.7 of its own motif: the phenol motif at its carbonyl and
    // para carbons has no O-H, so the weight holds up to pH 10 and falls above it.
    double pyridone10=identityConfidence("[nH]1ccccc1=O","n1ccccc1O",10);
    assertTrue(pyridone10>0.9&&identityConfidence("[nH]1ccccc1=O","n1ccccc1O",12)<pyridone10,
        "2-Pyridone loses its neutral form above its N-H pKa");
    MolGraph hydroxamic=rawGraph("ONC(C)=O");
    assertEquals(MolGraph.TAUTOMER_SITE_ACID,hydroxamic.tautomerSite[1]);
    assertEquals(MolGraph.TPKA_HYDROXAMIC,hydroxamic.tautomerPKa[1]);
    assertEquals(MolGraph.TW_HYDROXAMIC,hydroxamic.tautomerWeight[1],"Graphs keep the baseline weight; the pH applies when it is read");
    MolGraph nitro=rawGraph("C[N+](=O)[O-]");
    for(int atom=0;atom<nitro.n;atom++)assertEquals(MolGraph.TAUTOMER_SITE_NONE,nitro.tautomerSite[atom]);
  }
  @Test void mostSpecificMotifSetsWeightAndSiteInEitherAtomOrder()throws Exception {
    for(int[] atoms:new int[][]{{0,1,4},{1,3,2}}) {
      MolGraph graph=rawGraph(atoms[0]==0?"ONC(C)=O":"CC(=O)NO");
      assertEquals(MolGraph.TAUTOMER_SITE_ACID,graph.tautomerSite[atoms[1]]);
      assertEquals(MolGraph.TPKA_HYDROXAMIC,graph.tautomerPKa[atoms[1]],"The hydroxamic acid motif comes before the amide motif");
      assertEquals(MolGraph.TW_HYDROXAMIC,graph.tautomerWeight[atoms[2]]);
      assertEquals(MolGraph.TPKA_HYDROXAMIC,graph.tautomerPKa[atoms[2]]);
    }
    for(String smiles:new String[]{"[nH]1ccccc1=O","O=c1cccc[nH]1"}) {
      MolGraph graph=rawGraph(smiles);
      for(int atom:new int[]{0,6}) {
        assertEquals(MolGraph.TW_PYRIDINONE,graph.tautomerWeight[atom],smiles);
        assertEquals(MolGraph.TAUTOMER_SITE_ACID,graph.tautomerSite[atom],smiles);
        assertEquals(MolGraph.TPKA_PYRIDINONE,graph.tautomerPKa[atom],"The pyridinone motif comes before the amide motif: "+smiles);
      }
    }
    // The phenol motif reaches the carbonyl and para carbons of 2-pyridone, but its oxygen holds no
    // hydrogen, so the class keeps the N-H pKa of the pyridinone motif alone.
    MolGraph pyridone=rawGraph("[nH]1ccccc1=O");
    for(int atom:new int[]{2,5}) {
      assertEquals(pyridone.tautomerClass[0],pyridone.tautomerClass[atom]);
      assertEquals(MolGraph.TW_PHENOL_QUINONE,pyridone.tautomerWeight[atom]);
      assertEquals(MolGraph.TAUTOMER_SITE_NONE,pyridone.tautomerSite[atom],"A lactam carbonyl gives the phenol motif no acid site");
    }
    MolGraph urea=rawGraph("NC(N)=O");
    for(int atom:new int[]{0,2,3}) {
      assertEquals(MolGraph.TW_UREA,urea.tautomerWeight[atom]);
      assertEquals(MolGraph.TPKA_UREA,urea.tautomerPKa[atom],"The urea motif comes before the amide motif");
    }
    MolGraph succinimide=rawGraph("N1C(=O)CCC1=O");
    assertEquals(MolGraph.TPKA_IMIDE,succinimide.tautomerPKa[0],"A cyclic imide keeps its imide N-H site");
    assertEquals(MolGraph.TW_IMIDE,succinimide.tautomerWeight[1]);
    assertEquals(MolGraph.TW_BARBITURIC,succinimide.tautomerWeight[2]);
    assertEquals(MolGraph.TAUTOMER_SITE_NONE,succinimide.tautomerSite[2]);
    MolGraph acetylacetone=rawGraph("CC(=O)CC(C)=O");
    for(int atom:new int[]{2,3,6})
      assertEquals(MolGraph.TPKA_DIKETONE_ENOL,acetylacetone.tautomerPKa[atom],"The 1,3-diketone motif comes before the keto/enol motif");
    assertEquals(MolGraph.TPKA_KETO_ENOL,acetylacetone.tautomerPKa[0]);
  }
  @Test void eachClassTakesOneSiteModel()throws Exception {
    // 2-Hydroxypyridine: the phenol motif reaches the para carbon through the O-H of the
    // pyridone/hydroxypyridine motif, so that carbon carries the pKa of that motif, not the phenol pKa.
    MolGraph hydroxypyridine=rawGraph("Oc1ccccn1");
    assertEquals(hydroxypyridine.tautomerClass[0],hydroxypyridine.tautomerClass[4]);
    assertEquals(MolGraph.TW_PHENOL_QUINONE,hydroxypyridine.tautomerWeight[4]);
    for(int atom:new int[]{0,4}) {
      assertEquals(MolGraph.TAUTOMER_SITE_ACID,hydroxypyridine.tautomerSite[atom]);
      assertEquals(MolGraph.TPKA_PYRIDONE_HYDROXYPYRIDINE,hydroxypyridine.tautomerPKa[atom]);
    }
    for(double pH:new double[]{2,7.4,10,12}) {
      float[] fractions=TautomerClasses.neutralFractions(hydroxypyridine,pH);
      float acid=(float)(1.0/(1.0+StrictMath.pow(10.0,pH-MolGraph.TPKA_PYRIDONE_HYDROXYPYRIDINE)));
      for(int atom:new int[]{0,1,4,6})assertEquals(acid,fractions[atom],"pH "+pH);
      assertTrue(Float.isNaN(fractions[2]),"An atom outside every class keeps its weight");
    }
    // 4(3H)-Pyrimidinone takes its own N-H pKa in either tautomer, and the hypoxanthine N1-H too.
    for(String smiles:new String[]{"O=c1ccnc[nH]1","Oc1ccncn1"}) {
      MolGraph pyrimidinone=rawGraph(smiles);
      for(int atom:new int[]{0,6}) {
        assertEquals(MolGraph.TAUTOMER_SITE_ACID,pyrimidinone.tautomerSite[atom],smiles);
        assertEquals(8.6f,pyrimidinone.tautomerPKa[atom],smiles);
      }
    }
    assertEquals(8.6f,rawGraph("O=c1[nH]cnc2[nH]cnc12").tautomerPKa[2],"Hypoxanthine N1-H");
    float[] acetone=TautomerClasses.neutralFractions(rawGraph("CC(=O)C"),12);
    assertEquals((float)(1.0/(1.0+StrictMath.pow(10.0,12-MolGraph.TPKA_KETO_ENOL))),acetone[2],"One acid site keeps the one-site formula");
    for(float fraction:TautomerClasses.neutralFractions(rawGraph("CC(=N)N"),12))
      assertTrue(Float.isNaN(fraction),"A motif to which the library assigns no pKa is not changed by pH");
  }
  @Test void aSiteNeedsTheHydrogenItWouldRelease()throws Exception {
    // No phenol site on an ether, an ester or a lactone oxygen, no hydroxamic site on an N-alkoxy
    // amide and no amide site on a tertiary amide.
    for(String smiles:new String[]{"COc1ccccc1","O=c1ccc2ccccc2o1","CON(C)C(C)=O","CONC(C)=O","CN(C)C=O",
        "CN1C(=O)CCC1=O","Cn1ccccc1=O"}) {
      MolGraph graph=rawGraph(smiles);
      for(int atom=0;atom<graph.n;atom++)assertEquals(MolGraph.TAUTOMER_SITE_NONE,graph.tautomerSite[atom],smiles);
    }
    MolGraph arylEster=rawGraph("CC(=O)Oc1ccccc1");
    assertTrue(arylEster.tautomerClass[3]>=0);
    assertEquals(MolGraph.TAUTOMER_SITE_NONE,arylEster.tautomerSite[3],"An aryl ester oxygen carries no phenol site");
    assertEquals(MolGraph.TAUTOMER_SITE_NONE,arylEster.tautomerSite[4]);
    assertEquals(MolGraph.TAUTOMER_SITE_ACID,arylEster.tautomerSite[0],"The acetyl keeps the keto/enol site");
    MolGraph phenol=rawGraph("Oc1ccccc1");
    assertEquals(MolGraph.TAUTOMER_SITE_ACID,phenol.tautomerSite[0]);
    assertEquals(MolGraph.TPKA_PHENOL_QUINONE,phenol.tautomerPKa[0]);
    MolGraph weinreb=rawGraph("CON(C)C(C)=O");
    assertTrue(weinreb.tautomerClass[2]>=0);
    assertEquals(MolGraph.TW_HYDROXAMIC,weinreb.tautomerWeight[2],"An N-alkoxy amide keeps the hydroxamic class and weight");
  }
  /** The molecule with its atoms and bonds in a random order; permutation[atom] is each new atom's original index. */
  private static IAtomContainer reordered(IAtomContainer molecule,Random random,int[] permutation) {
    List<Integer> order=new ArrayList<>();
    for(int atom=0;atom<molecule.getAtomCount();atom++)order.add(atom);
    Collections.shuffle(order,random);
    IAtomContainer copy=SilentChemObjectBuilder.getInstance().newAtomContainer();
    for(int atom=0;atom<order.size();atom++) {permutation[atom]=order.get(atom);copy.addAtom(molecule.getAtom(order.get(atom)));}
    List<IBond> bonds=new ArrayList<>();
    for(IBond bond:molecule.bonds())bonds.add(bond);
    Collections.shuffle(bonds,random);
    for(IBond bond:bonds)copy.addBond(bond);
    return copy;
  }
  /** The lowest original atom index in each atom's class, by original atom, or -1 outside every class. */
  private static int[] classRepresentatives(MolGraph graph,int[] original) {
    Map<Integer,Integer> lowest=new HashMap<>();
    for(int atom=0;atom<graph.n;atom++)
      if(graph.tautomerClass[atom]>=0)lowest.merge(graph.tautomerClass[atom],original[atom],Math::min);
    int[] representatives=new int[graph.n];
    Arrays.fill(representatives,-1);
    for(int atom=0;atom<graph.n;atom++)
      if(graph.tautomerClass[atom]>=0)representatives[original[atom]]=lowest.get(graph.tautomerClass[atom]);
    return representatives;
  }
  /**
   * Whether a symmetry of the graph carries the class numbers {@code first} onto {@code second}
   * (both by atom): an induced embedding of the graph labelled with {@code first} into the graph
   * labelled with {@code second}, each number written into the mass number of its atom.
   */
  private static boolean symmetricNumbering(MolGraph graph,int[] first,int[] second) {
    SearchResult embedding=Substructure.find(classLabelled(graph,first)).withOptions(symmetryOptions())
        .match(classLabelled(graph,second));
    return !embedding.matches().isEmpty();
  }
  private static Substructure.Options symmetryOptions() {
    Substructure.Options options=new Substructure.Options();
    options.chemistry.matchFormalCharge=true;options.chemistry.matchIsotope=true;
    options.chemistry.aromaticityMode=ChemOptions.AromaticityMode.STRICT;
    options.chemistry.matchBondOrder=ChemOptions.BondOrderMode.STRICT;
    options.chemistry.ringMatchesRingOnly=true;
    options.induced=true;options.timeoutMillis=60_000;
    return options;
  }
  private static MolGraph classLabelled(MolGraph graph,int[] classes) {
    int n=graph.n;
    int[] masses=new int[n],hydrogens=new int[n];
    int[][] orders=new int[n][];
    boolean[][] aromaticBonds=new boolean[n][],ringBonds=new boolean[n][];
    for(int atom=0;atom<n;atom++) {
      masses[atom]=graph.massNumber[atom]*1000+classes[atom]+2;
      hydrogens[atom]=graph.hydrogenCount(atom);
      int[] adjacent=graph.neighbors[atom];
      orders[atom]=new int[adjacent.length];aromaticBonds[atom]=new boolean[adjacent.length];ringBonds[atom]=new boolean[adjacent.length];
      for(int k=0;k<adjacent.length;k++) {
        orders[atom][k]=graph.bondOrder(atom,adjacent[k]);
        aromaticBonds[atom][k]=graph.bondAromatic(atom,adjacent[k]);
        ringBonds[atom][k]=graph.bondInRing(atom,adjacent[k]);
      }
    }
    return new MolGraph.Builder().atomCount(n).atomicNumbers(graph.atomicNum).formalCharges(graph.formalCharge)
        .massNumbers(masses).hydrogenCounts(hydrogens).aromaticFlags(graph.aromatic).ringFlags(graph.ring)
        .neighbors(graph.neighbors).bondOrders(orders).bondAromaticFlags(aromaticBonds).bondRingFlags(ringBonds).build();
  }
  /**
   * Random atom orders, with a fixed seed, change no class, weight or site: the classes partition
   * the same atoms in every order, and their numbers agree or differ only by a symmetry of the graph.
   */
  @Test void atomOrderDoesNotChangeClassesWeightsOrSites()throws Exception {
    Random random=new Random(20261021L);
    int molecules=0,symmetric=0;
    for(String line:Files.readAllLines(Paths.get("benchmarks/data/chembl_mcs_benchmark.smi"),StandardCharsets.UTF_8)) {
      if(molecules>=600)break;
      String smiles=line.trim().split("\\s+")[0];
      if(smiles.isEmpty())continue;
      IAtomContainer molecule;MolGraph graph;
      try {molecule=rawMolecule(smiles);graph=new MolGraph(molecule);} catch(Exception unparsable) {continue;}
      molecules++;
      TautomerClasses.Motifs motifs=TautomerClasses.motifs(graph);
      int[] identity=new int[graph.n];
      for(int atom=0;atom<graph.n;atom++)identity[atom]=atom;
      int[] classes=classRepresentatives(graph,identity);
      float[] fractions=TautomerClasses.neutralFractions(graph,12);
      for(int trial=0;trial<3;trial++) {
        int[] permutation=new int[graph.n];
        MolGraph copy=new MolGraph(reordered(molecule,random,permutation));
        TautomerClasses.Motifs copyMotifs=TautomerClasses.motifs(copy);
        for(int atom=0;atom<copy.n;atom++) {
          int original=permutation[atom];
          assertEquals(motifs.weight[original],copyMotifs.weight[atom],"Atom order changed a motif weight: "+smiles);
          assertEquals(motifs.site[original],copyMotifs.site[atom],"Atom order changed a motif site: "+smiles);
          assertEquals(motifs.pKa[original],copyMotifs.pKa[atom],"Atom order changed a motif pKa: "+smiles);
          assertEquals(copyMotifs.covered[atom],copy.tautomerClass[atom]>=0,"Every atom a motif covers is in a class: "+smiles);
          assertEquals(copyMotifs.weight[atom],copy.tautomerWeight[atom],"An atom in a class takes its motif weight: "+smiles);
          assertEquals(copyMotifs.site[atom],copy.tautomerSite[atom],smiles);
          assertEquals(copyMotifs.pKa[atom],copy.tautomerPKa[atom],smiles);
        }
        assertArrayEquals(classes,classRepresentatives(copy,permutation),"Atom order changed a tautomer class: "+smiles);
        int[] numbers=new int[graph.n];
        for(int atom=0;atom<copy.n;atom++)numbers[permutation[atom]]=copy.tautomerClass[atom];
        if(!Arrays.equals(numbers,graph.tautomerClass)) {
          assertTrue(symmetricNumbering(graph,graph.tautomerClass,numbers),
              "Atom order changed the class numbers beyond a symmetry of the graph: "+smiles);
          symmetric++;
        }
        float[] copyFractions=TautomerClasses.neutralFractions(copy,12);
        for(int atom=0;atom<copy.n;atom++)
          assertEquals(fractions[permutation[atom]],copyFractions[atom],"Atom order changed a neutral fraction: "+smiles);
      }
    }
    assertTrue(molecules>=500,"The permutation test needs at least 500 molecules, found "+molecules);
    // Classes that a symmetry exchanges take their numbers by atom order, so some orders differ.
    assertTrue(symmetric>0,"The corpus exercises classes that a symmetry exchanges");
  }
  /** The tautomer classes, by atom, of molecules sensitive to atom order, each in two atom orders. */
  @Test void classesOfReportedMoleculesDoNotDependOnTheAtomOrder()throws Exception {
    Object[][] cases={
        {"ONC(C)=O",new int[]{0,0,0,0,0}},{"CC(=O)NO",new int[]{0,0,0,0,0}},
        {"O=C1CCC(=O)N1",new int[]{0,0,0,0,0,0,0}},{"N1C(=O)CCC1=O",new int[]{0,0,0,0,0,0,0}},
        {"Oc1ccccn1",new int[]{0,0,-1,-1,0,-1,0}},{"n1ccccc1O",new int[]{0,-1,0,-1,-1,0,0}},
        // Two carbonyl classes that the mirror symmetry exchanges: the lower atoms take number 0.
        {"CC(=O)CCC(C)=O",new int[]{0,-1,0,0,1,-1,1,1}}};
    for(Object[] item:cases)
      assertArrayEquals((int[])item[1],rawGraph((String)item[0]).tautomerClass,(String)item[0]);
  }
  @Test void identicalSelectedChemistryHasNoPenaltyInAnyContext()throws Exception {
    MolGraph molecule=rawGraph("CC(=O)C");SearchResult.Match identity=match(molecule,molecule);
    for(ChemOptions.Solvent solvent:ChemOptions.Solvent.values())
      assertEquals(1.0,MCS.computeTautomerConfidence(molecule,molecule,identity,ChemOptions.tautomerProfile().withSolvent(solvent).withPH(0)));
  }
  @Test void unselectedBondOrderDifferencesDoNotAffectConfidence()throws Exception {
    MolGraph query=rawGraph("CC(=O)C"),target=rawGraph("CC(O)=C");SearchResult.Match mapping=match(query,target);
    Map<SearchResult.Bond,SearchResult.Bond> sameBonds=new LinkedHashMap<>();
    for(Map.Entry<SearchResult.Bond,SearchResult.Bond> pair:mapping.bondMapping().entrySet()) {
      SearchResult.Bond first=pair.getKey(),second=pair.getValue();
      if(query.bondOrder(first.first(),first.second())==target.bondOrder(second.first(),second.second()))sameBonds.put(first,second);
    }
    SearchResult.Match selected=new SearchResult.Match(mapping.atomMapping(),sameBonds,mapping.fragments(),mapping.score());
    assertEquals(1.0,MCS.computeTautomerConfidence(query,target,selected,ChemOptions.tautomerProfile()));
  }
  @Test void concurrentContextReadsLeaveSearchAndWeightsUnchanged()throws Exception {
    MolGraph query=rawGraph("CC(=O)C"),target=rawGraph("CC(O)=C");SearchResult.Match mapping=match(query,target);
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
    MolGraph query=rawGraph("CCO"),target=rawGraph("CCO");
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
    MolGraph molecule=rawGraph("CC(=O)C");SearchResult.Match mapping=match(molecule,molecule);
    ChemOptions options=ChemOptions.tautomerProfile();options.pH=pH;
    assertThrows(IllegalArgumentException.class,()->MCS.computeTautomerConfidence(molecule,molecule,mapping,options));
  }
}

/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <iostream>
#include <set>
#include <stdexcept>
static void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
static void check(const smsd::SearchResult& result, size_t expected, int atoms, int bonds, double score) {
 require(result.status == smsd::SearchStatus::Completed && result.optimal && result.complete,
     "Enumeration must complete with proof");
 require(result.matches.size() == expected && result.score == score, "Every optimum embedding must remain");
 std::set<std::map<int,int>> maps;
 for (const auto& match : result.matches) {
  require(match.atomMapping.size() == static_cast<size_t>(atoms) && match.bondMapping.size() == static_cast<size_t>(bonds),
      "Each selected atom and bond must appear exactly once");
  require(maps.insert(match.atomMapping).second, "Embeddings must be distinct");
 }
}
int main() {
 try {
  const auto target=smsd::parseSMILES("CCC");
  for (bool both : {false,true}) {
   auto query=smsd::parseSMILES("CC"); query.neighbors[0].push_back(1);
   if (both) query.neighbors[1].push_back(0);
   for (bool induced : {false,true}) for (bool atomType : {false,true}) {
    smsd::SubstructureOptions options; options.timeoutMs=1000; options.induced=induced; options.chemistry.matchAtomType=atomType;
    auto result=smsd::Substructure::find(query).withOptions(options).matchAll(target,100);
    check(result,4,2,1,2);
    for (const auto& match:result.matches) require(smsd::MCS::validateMapping(smsd::parseSMILES("CC"),target,match).empty(),"Sub witnesses must validate against the same bond graph");
   }
   for (auto connectivity : {smsd::MCSConnectivity::Connected,smsd::MCSConnectivity::Fragments}) {
    smsd::MCSOptions options; options.timeoutMs=1000; options.connectivity=connectivity;
    auto prepared=smsd::MCS::find(query).withOptions(options);
    auto result=prepared.matchAll(target,100); check(result,4,2,1,2);
    for (const auto& match:result.matches) require(smsd::MCS::validateMapping(smsd::parseSMILES("CC"),target,match,options).empty(),"MCIS witnesses must validate against the same bond graph");
    auto capped=prepared.matchAll(target,2);require(capped.matches.size()==2 && !capped.complete && capped.status==smsd::SearchStatus::LimitReached,"Enumeration cap must remain visible");
    options.atomWeights={1,2}; check(smsd::MCS::find(query).withOptions(options).matchAll(target,100),4,2,1,3);
   }
   const auto mapped=smsd::MCS::mappingMatch(query,target,{{0,0},{1,1}});
   require(mapped.bondMapping.size()==1,"mappingMatch must select one correspondence per bond");
  }
  auto query=smsd::parseSMILES("CCC");query.neighbors[0].push_back(1);query.neighbors[0].push_back(1);
  smsd::MCSOptions options;options.timeoutMs=1000;options.connectivity=smsd::MCSConnectivity::Fragments;
  check(smsd::MCS::find(query).withOptions(options).matchAll(smsd::parseSMILES("CC"),100),4,2,1,2);
  auto sameTerms=smsd::parseSMILES("C(C)(C)C.C(C)C");
  sameTerms.neighbors[4].push_back(5); sameTerms.neighbors[5].push_back(4);
  sameTerms.degree[4]=3; sameTerms.degree[5]=2;
  const auto plain=smsd::parseSMILES("C(C)(C)C.C(C)C");
  smsd::SubstructureOptions charged; charged.timeoutMs=1000;charged.chemistry.matchFormalCharge=true;
  check(smsd::Substructure::find(sameTerms).withOptions(charged).matchAll(plain,100),12,7,5,7);
  std::cout<<"PASS duplicate query bonds and complete embeddings\n";return 0;
 } catch (const std::exception& error) {std::cerr<<error.what()<<'\n';return 1;}
}

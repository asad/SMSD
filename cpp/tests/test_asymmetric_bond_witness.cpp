/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/mcs.hpp"
#include "smsd/substructure.hpp"
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Witness = std::pair<std::map<int,int>,std::vector<std::array<int,4>>>;
struct Expected { int atoms=0,bonds=0;std::set<Witness> witnesses; };
struct Fixture { std::string name;smsd::MolGraph query,target;bool strictAromatic=false,ringOnly=false; };
static int failures=0,checks=0;
static void require(bool condition,const std::string& context) {++checks;if(!condition){++failures;std::cerr<<context<<'\n';}}
static Witness witness(const smsd::Match& match) {auto bonds=match.bondMapping;std::sort(bonds.begin(),bonds.end());return {match.atomMapping,bonds};}
static std::string text(const Witness& key) {std::ostringstream out;for(const auto& atom:key.first)out<<atom.first<<'>'<<atom.second<<',';out<<'/';for(const auto& bond:key.second)out<<bond[0]<<'-'<<bond[1]<<'>'<<bond[2]<<'-'<<bond[3]<<',';return out.str();}

// Independent exhaustive injection oracle for these small element/order/ring/aromatic fixtures.
// Its bond orientation is the public selected-bond convention: lower query index first.
static Expected oracle(const Fixture& fixture,bool substructure,bool induced,bool connected) {
    const auto& query=fixture.query;const auto& target=fixture.target;
    Expected best;std::vector<int> mapping(query.n,-1);std::vector<bool> used(target.n,false);
    auto compatible=[&](int queryFirst,int querySecond,int targetFirst,int targetSecond) {
        if(query.bondOrder(queryFirst,querySecond)!=target.bondOrder(targetFirst,targetSecond))return false;
        if(fixture.strictAromatic && query.bondAromatic(queryFirst,querySecond)!=target.bondAromatic(targetFirst,targetSecond))return false;
        return !fixture.ringOnly || query.bondInRing(queryFirst,querySecond)==target.bondInRing(targetFirst,targetSecond);
    };
    std::function<void(int)> visit=[&](int queryAtom) {
        if(queryAtom==query.n) {
            Witness key;
            for(int atom=0;atom<query.n;++atom)if(mapping[atom]>=0)key.first[atom]=mapping[atom];
            if(key.first.empty() || substructure && key.first.size()!=static_cast<size_t>(query.n))return;
            for(const auto& first:key.first)for(const auto& second:key.first)if(first.first<second.first) {
                const bool source=query.hasBond(first.first,second.first),image=target.hasBond(first.second,second.second);
                if(source && (!image || !compatible(first.first,second.first,first.second,second.second)))return;
                if(induced && source!=image)return;
                if(source)key.second.push_back({first.first,second.first,first.second,second.second});
            }
            if(connected) {
                std::set<int> reached{key.first.begin()->first};bool changed=true;
                while(changed){changed=false;for(const auto& bond:key.second)if(reached.count(bond[0])!=reached.count(bond[1])){reached.insert(bond[0]);reached.insert(bond[1]);changed=true;}}
                if(reached.size()!=key.first.size())return;
            }
            const int atoms=static_cast<int>(key.first.size()),bonds=static_cast<int>(key.second.size());
            if(std::make_pair(atoms,bonds)>std::make_pair(best.atoms,best.bonds)){best.atoms=atoms;best.bonds=bonds;best.witnesses.clear();}
            if(atoms==best.atoms && bonds==best.bonds)best.witnesses.insert(std::move(key));
            return;
        }
        if(!substructure)visit(queryAtom+1);
        for(int targetAtom=0;targetAtom<target.n;++targetAtom) {
            if(used[targetAtom] || query.atomicNum[queryAtom]!=target.atomicNum[targetAtom])continue;
            if(fixture.ringOnly && query.ring[queryAtom]!=target.ring[targetAtom])continue;
            if(fixture.strictAromatic && query.aromatic[queryAtom]!=target.aromatic[targetAtom])continue;
            used[targetAtom]=true;mapping[queryAtom]=targetAtom;visit(queryAtom+1);mapping[queryAtom]=-1;used[targetAtom]=false;
        }
    };
    visit(0);return best;
}

static smsd::MolGraph permute(const smsd::MolGraph& original,const std::vector<int>& order) {
    auto graph=original;std::vector<int> inverse(order.size());
    for(size_t atom=0;atom<order.size();++atom)inverse[order[atom]]=static_cast<int>(atom);
    for(int atom=0;atom<original.n;++atom) {
        const int old=order[atom];
        graph.atomicNum[atom]=original.atomicNum[old];graph.formalCharge[atom]=original.formalCharge[old];
        graph.massNumber[atom]=original.massNumber[old];graph.hydrogenCount[atom]=original.hydrogenCount[old];
        graph.atomClass[atom]=original.atomClass[old];graph.ring[atom]=original.ring[old];graph.aromatic[atom]=original.aromatic[old];
        graph.neighbors[atom].clear();for(int neighbor:original.neighbors[old])graph.neighbors[atom].push_back(inverse[neighbor]);
        graph.degree[atom]=static_cast<int>(graph.neighbors[atom].size());
        std::fill(graph.adjLong[atom].begin(),graph.adjLong[atom].end(),uint64_t(0));
        for(int neighbor:graph.neighbors[atom])graph.adjLong[atom][neighbor>>6]|=uint64_t(1)<<(neighbor&63);
        for(int other=0;other<original.n;++other)graph.bondPacked[atom][other]=original.bondPacked[old][order[other]];
    }
    graph.sssrRings.clear();graph.sssrComputed=false;graph.ringCount.clear();graph.ringCountsComputed_=false;
    graph.cachedNeighborsByDegDesc.clear();graph.invalidateBondDerivedState();graph.refreshAtomLabels();return graph;
}

static void check(const std::string& label,const Fixture& fixture,const smsd::SearchResult& result,
                  const Expected& expected,bool all,bool trace) {
    require(result.status==smsd::SearchStatus::Completed && result.optimal && result.complete,label+": no complete proof");
    require(result.score==expected.atoms,label+": wrong atom score");
    std::set<Witness> found;
    for(const auto& match:result.matches) {
        const auto key=witness(match);
        require(found.insert(key).second,label+": repeated witness");
        require(expected.witnesses.count(key)!=0,label+": unexpected atom/bond witness "+text(key));
    }
    if(all)require(found==expected.witnesses,label+": exact witness set differs");
    else require(found.empty()==expected.witnesses.empty() && found.size()<=1,label+": wrong single result cardinality");
    if(trace){std::cout<<"ORDINARY|"<<label<<'|'<<int(result.status)<<'|'<<result.optimal<<'|'<<result.complete<<'|'<<result.score<<'|'<<result.nodesVisited;for(const auto& match:result.matches)std::cout<<'|'<<text(witness(match));std::cout<<'\n';}
}

static void molecular(Fixture fixture) {
    std::vector<int> permutation(fixture.query.n);std::iota(permutation.begin(),permutation.end(),0);int index=0;
    do {
        auto current=fixture;current.query=permute(fixture.query,permutation);
        smsd::ChemOptions chemistry;
        if(fixture.strictAromatic)chemistry.aromaticityMode=smsd::ChemOptions::AromaticityMode::STRICT;
        chemistry.ringMatchesRingOnly=fixture.ringOnly;
        for(bool induced:{false,true}) {
            smsd::SubstructureOptions options;options.timeoutMs=1000;options.induced=induced;options.chemistry=chemistry;
            const auto prepared=smsd::Substructure::find(current.query).withOptions(options);
            const auto expected=oracle(current,true,induced,false);
            const std::string label=fixture.name+"/perm"+std::to_string(index)+"/Sub"+(induced?"Induced":"");
            check(label+"/single",current,prepared.match(current.target),expected,false,fixture.name=="ordinary");
            check(label+"/all",current,prepared.matchAll(current.target,1000),expected,true,fixture.name=="ordinary");
        }
        for(bool connected:{false,true}) {
            smsd::MCSOptions options;options.timeoutMs=1000;options.chemistry=chemistry;options.connectivity=connected?smsd::MCSConnectivity::Connected:smsd::MCSConnectivity::Fragments;
            const auto prepared=smsd::MCS::find(current.query).withOptions(options);
            const auto expected=oracle(current,false,true,connected);
            const std::string label=fixture.name+"/perm"+std::to_string(index)+"/MCIS"+(connected?"Connected":"FragmentsDefault");
            check(label+"/single",current,prepared.match(current.target),expected,false,fixture.name=="ordinary");
            check(label+"/all",current,prepared.matchAll(current.target,1000),expected,true,fixture.name=="ordinary");
        }
        ++index;
    }while(std::next_permutation(permutation.begin(),permutation.end()));
}

static void smarts() {
    auto target=smsd::parseSMILES("CCC");target.bondPacked[1][0]=2;
    const Fixture fixture{"smartsexplicit",smsd::parseSMILES("CC"),target};
    for(bool induced:{false,true}) {
        smsd::SubstructureOptions options;options.timeoutMs=1000;options.induced=induced;
        const auto prepared=smsd::Substructure::smarts("[#6]-[#6]").withOptions(options);
        const auto expected=oracle(fixture,true,induced,false);
        check("SMARTS-explicit/single",fixture,prepared.match(target),expected,false,false);
        check("SMARTS-explicit/all",fixture,prepared.matchAll(target,1000),expected,true,false);
    }
    const auto recursive=smsd::Substructure::smarts("[$([#6]-[#6])]").withTimeout(1000).matchAll(target,8);
    std::set<int> atoms;for(const auto& match:recursive.matches)if(match.atomMapping.size()==1)atoms.insert(match.atomMapping.at(0));
    require(recursive.status==smsd::SearchStatus::Completed && recursive.optimal && recursive.complete
            && atoms==std::set<int>({0,1,2}) && recursive.matches.size()==3,"Recursive single-bond SMARTS must find every carbon with a canonical outgoing single bond");
}

int main(){
    const auto ethane=smsd::parseSMILES("CC"),propane=smsd::parseSMILES("CCC");
    std::vector<Fixture> fixtures{{"ordinary",ethane,propane}};
    for(int property=0;property<3;++property)for(bool querySide:{false,true}) {
        Fixture fixture{(querySide?"query_":"target_")+std::string(property==0?"order":property==1?"aromatic":"ring"),ethane,propane};
        auto& graph=querySide?fixture.query:fixture.target;
        if(property==0)graph.bondPacked[1][0]=2;
        if(property==1){graph.bondPacked[1][0]|=smsd::MolGraph::BOND_AROM_BIT;fixture.strictAromatic=true;}
        if(property==2){graph.bondPacked[1][0]|=smsd::MolGraph::BOND_RING_BIT;fixture.ringOnly=true;}
        fixtures.push_back(fixture);
    }
    Fixture partial{"partial_mcis",smsd::parseSMILES("CCN"),propane};partial.target.bondPacked[1][0]=2;fixtures.push_back(partial);
    Fixture degree{"different_degree",smsd::parseSMILES("CCO"),smsd::parseSMILES("CCCO")};degree.target.bondPacked[2][1]=2;fixtures.push_back(degree);
    for(const auto& fixture:fixtures)molecular(fixture);
    smarts();
    std::cout<<"checks="<<checks<<" failures="<<failures<<'\n';return failures==0?0:1;
}

// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2018-2026 BioInception PVT LTD
// Real SMSD native MCS benchmark. Parsing/warmup/validation are outside timers.
// c++ -std=c++17 -O3 -Icpp/include benchmarks/benchmark_cpp.cpp -o benchmark_cpp
// ./benchmark_cpp [pairs.tsv] [timeout_ms=10000] [repeats=3] [warmup=1]
// TSV fields: SMILES1, SMILES2, optional name. Comments start with '#'.
// SMSD mapping API does not expose cancellation; report budget crossings only.
#include "smsd/smsd.hpp"
#include "smsd/smiles_parser.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Pair { std::string first, second, name; };

static bool validMapping(const smsd::MolGraph& a, const smsd::MolGraph& b,
                         const std::map<int,int>& mapping) {
    std::set<int> targets;
    for (auto [q,t] : mapping) {
        if (q < 0 || q >= a.n || t < 0 || t >= b.n || !targets.insert(t).second
            || a.atomicNum[q] != b.atomicNum[t]) return false;
        for (int n : a.neighbors[q]) {
            auto other = mapping.find(n);
            if (other != mapping.end() && !b.hasBond(t,other->second)) return false;
        }
    }
    if (mapping.empty()) return true;
    std::set<int> visited;
    std::vector<int> pending{mapping.begin()->first};
    while (!pending.empty()) {
        int q = pending.back(); pending.pop_back();
        if (!visited.insert(q).second) continue;
        for (int n : a.neighbors[q]) if (mapping.count(n)) pending.push_back(n);
    }
    return visited.size() == mapping.size();
}

int main(int argc, char** argv) {
    try {
        if (argc > 5) throw std::invalid_argument("usage: benchmark_cpp [pairs.tsv] [timeout_ms] [repeats] [warmup]");
        int timeout = argc > 2 ? std::stoi(argv[2]) : 10000;
        int repeats = argc > 3 ? std::stoi(argv[3]) : 3;
        int warmup = argc > 4 ? std::stoi(argv[4]) : 1;
        if (timeout < 1 || repeats < 1 || warmup < 0) throw std::invalid_argument("invalid measurement settings");
        std::vector<Pair> pairs{{"C","CC","methane-ethane"},
            {"c1ccccc1","Cc1ccccc1","benzene-toluene"},
            {"CC(=O)O","CC(=O)OC","acid-ester"},
            {"C1CCCCC1","C1CCC(O)CC1","cyclohexane-cyclohexanol"},
            {"C1C2CC3CC1CC(C2)C3","C1C2CC3CC1CC(C2)C3","adamantane-self"},
            {"OCCOCCOCCO","OCCOCCOCCOCCOCCO","PEG-short-long"}};
        if (argc > 1) {
            pairs.clear(); std::ifstream input(argv[1]);
            if (!input) throw std::runtime_error("cannot read pair file");
            std::string line;
            while (std::getline(input,line)) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream row(line); Pair pair;
                std::getline(row,pair.first,'\t'); std::getline(row,pair.second,'\t'); std::getline(row,pair.name,'\t');
                if (pair.second.empty()) throw std::runtime_error("invalid TSV pair");
                pairs.push_back(std::move(pair));
            }
        }
        smsd::ChemOptions chem; chem.matchBondOrder=smsd::ChemOptions::BondOrderMode::ANY;
        chem.aromaticityMode=smsd::ChemOptions::AromaticityMode::FLEXIBLE;
        smsd::MCSOptions options; options.timeoutMs=timeout; options.maxStage=5;
        options.connectedOnly=true; options.maximizeBonds=false; options.induced=false;
        std::cout << "# engine=SMSD native; policy=bond-any; objective=atoms; connected=true; timeout_status=unknown\n"
                  << "pair\ttrial\telapsed_us\tatoms\tvalid\tbudget_reached\n";
        for (const auto& pair : pairs) {
            auto a=smsd::parseSMILES(pair.first), b=smsd::parseSMILES(pair.second);
            if (a.n == 0 || b.n == 0) throw std::runtime_error("empty input graph: "+pair.name);
            for (int i=0;i<warmup;++i) smsd::findMCS(a,b,chem,options);
            for (int trial=0;trial<repeats;++trial) {
                auto start=std::chrono::steady_clock::now();
                auto mapping=smsd::findMCS(a,b,chem,options);
                double us=std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-start).count();
                std::cout << pair.name << '\t' << trial << '\t' << us << '\t' << mapping.size()
                          << '\t' << validMapping(a,b,mapping) << '\t' << (us >= timeout*1000.0) << std::endl;
            }
        }
        return 0;
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}

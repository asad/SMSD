/* SPDX-License-Identifier: Apache-2.0 */
// Reproduce the connectivity-filter comparison locally (no hosted runner):
//   mkdir -p /tmp/smsd-mcs-components-baseline
//   git archive 49d7303 cpp/include | tar -x -C /tmp/smsd-mcs-components-baseline
//   clang++ -std=c++17 -O3 -I/tmp/smsd-mcs-components-baseline/cpp/include \
//       benchmarks/benchmark_mcs_components.cpp -o /tmp/mcs-components-before
//   clang++ -std=c++17 -O3 -Icpp/include \
//       benchmarks/benchmark_mcs_components.cpp -o /tmp/mcs-components-after
//   /tmp/mcs-components-before
//   /tmp/mcs-components-after
// Run each executable three times and compare median milliseconds per case.
// Baseline cpp/include/smsd/mcs.hpp SHA-256:
//   5d293389c224d2c2e1a315872cb53dbe83953a615cc9e73654d0ce5ea4e62344
// This measures largestConnected only; graph setup is outside the timed region.
// It does not establish an end-to-end MCS speedup.
#include "smsd/mcs.hpp"
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <map>
#include <vector>

smsd::MolGraph path(int atoms, bool fragments) {
    std::vector<std::vector<int>> neighbors(atoms), orders(atoms);
    for (int atom = 1; atom < atoms; ++atom) {
        if (fragments && atom % 20 == 0) continue;
        neighbors[atom-1].push_back(atom); orders[atom-1].push_back(1);
        neighbors[atom].push_back(atom-1); orders[atom].push_back(1);
    }
    return smsd::MolGraph::Builder().atomCount(atoms)
        .atomicNumbers(std::vector<int>(atoms,6))
        .setNeighbors(neighbors).setBondOrders(orders).build();
}
int main() {
    volatile size_t sink=0;
    for (auto spec : {std::array<int,3>{100,100,20000}, {1000,1000,2000}, {1000,10,20000}, {100000,10,1000}}) {
        for (bool fragmented : {false,true}) {
            auto g=path(spec[0],fragmented);
            std::map<int,int> m;
            for(int atom=0;atom<spec[1];++atom)m.emplace(atom,atom);
            auto start=std::chrono::steady_clock::now();
            for(int i=0;i<spec[2];++i)sink+=smsd::detail::largestConnected(g,m,&g).size();
            auto ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            std::cout<<"atoms="<<spec[0]<<" mapped="<<spec[1]<<" fragmented="<<fragmented<<" iterations="<<spec[2]<<" ms="<<ms<<"\n";
        }
    }
    assert(sink == 4860000);
    std::cout<<"sink="<<sink<<"\n";
}

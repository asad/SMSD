/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 *
 * SMSD Pro 6.0.0 — ZINC20 Tautomer Over-Matching Benchmark (C++)
 * ===============================================================
 * Paper Item 15 (Future Work): quantify how tautomer-aware MCS relaxation
 * causes over-matching (matching atoms that violate proton conservation).
 *
 * Uses the tautomer section (lines 801-900) of diverse_molecules.txt,
 * forming consecutive keto/enol pairs.  For each pair:
 *   1. Run MCS with ChemOptions::tautomerProfile() (tautomer-aware)
 *   2. Run MCS with default ChemOptions()          (strict)
 *   3. Call validateTautomerConsistency() on the tautomer result
 *   4. Compute TautConf = avg(tautomerWeight) for mapped tautomeric atoms
 *
 * Compilation:
 *   clang++ -std=c++17 -O2 -I cpp/include benchmarks/benchmark_tautomer.cpp \
 *           -o /tmp/benchmark_tautomer && /tmp/benchmark_tautomer
 */

#include "smsd/smsd.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace smsd;

// ============================================================================
// Molecule record
// ============================================================================
struct Mol {
    std::string smiles;
    std::string name;
    MolGraph    graph;
};

// ============================================================================
// Load tautomer molecules (section 7: mol indices 801-900, 1-based)
// ============================================================================
static std::vector<Mol> loadTautomerMolecules(const std::string& path) {
    std::vector<Mol> mols;
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        fprintf(stderr, "ERROR: cannot open %s\n", path.c_str());
        return mols;
    }

    constexpr int TAUT_START = 801;
    constexpr int TAUT_END   = 900;

    std::string line;
    int molIdx = 0;
    ParseOptions popts;
    popts.lenient = true;

    while (std::getline(ifs, line)) {
        // Trim
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' '))
            line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        molIdx++;
        if (molIdx < TAUT_START || molIdx > TAUT_END) continue;

        // Split on tab
        std::string smi, name;
        auto tab = line.find('\t');
        if (tab != std::string::npos) {
            smi  = line.substr(0, tab);
            name = line.substr(tab + 1);
        } else {
            smi  = line;
            name = "mol_" + std::to_string(molIdx);
        }

        try {
            MolGraph g = parseSMILES(smi, popts);
            mols.push_back({smi, name, std::move(g)});
        } catch (const std::exception& e) {
            fprintf(stderr, "  SKIP: %s [%s] -> %s\n", name.c_str(), smi.c_str(), e.what());
        }
    }
    return mols;
}

// ============================================================================
// TautConf score: average tautomer weight for mapped tautomeric atoms
// ============================================================================
static double computeTautConfScore(const MolGraph& g1, const MolGraph& g2,
                                   const std::map<int,int>& mcs) {
    if (mcs.empty()) return 1.0;
    if (g1.tautomerClass.empty() || g2.tautomerClass.empty()) return 1.0;

    double weightSum = 0.0;
    int tautCount = 0;
    for (auto& [qi, ti] : mcs) {
        if (qi < 0 || qi >= g1.n || ti < 0 || ti >= g2.n) continue;
        bool isTaut = (g1.tautomerClass[qi] >= 0) || (g2.tautomerClass[ti] >= 0);
        if (isTaut) {
            double w1 = (qi < (int)g1.tautomerWeight.size()) ? g1.tautomerWeight[qi] : 1.0;
            double w2 = (ti < (int)g2.tautomerWeight.size()) ? g2.tautomerWeight[ti] : 1.0;
            weightSum += (w1 + w2) / 2.0;
            tautCount++;
        }
    }
    return tautCount > 0 ? weightSum / tautCount : 1.0;
}

// ============================================================================
// Pair result
// ============================================================================
struct PairResult {
    std::string nameA, nameB;
    int tautMcsSize;
    int defaultMcsSize;
    int overMatchDelta;
    bool protonConsistent;
    double tautConfScore;
    double tautTimeMs;
    double defaultTimeMs;
};

// ============================================================================
// Benchmark one pair
// ============================================================================
static PairResult benchmarkPair(Mol& a, Mol& b) {
    ChemOptions tautOpts = ChemOptions::tautomerProfile();
    ChemOptions defOpts;
    McsOptions mcsOpts;
    mcsOpts.timeoutMs = 10000;

    // Tautomer-aware MCS
    auto t0 = std::chrono::high_resolution_clock::now();
    std::map<int,int> tautMcs;
    try {
        tautMcs = findMCS(a.graph, b.graph, tautOpts, mcsOpts);
    } catch (...) {}
    auto t1 = std::chrono::high_resolution_clock::now();
    double tautTimeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // Default MCS
    t0 = std::chrono::high_resolution_clock::now();
    std::map<int,int> defMcs;
    try {
        defMcs = findMCS(a.graph, b.graph, defOpts, mcsOpts);
    } catch (...) {}
    t1 = std::chrono::high_resolution_clock::now();
    double defTimeMs = std::chrono::duration<double, std::milli>(t1 - t0).count();

    int tautSize = (int)tautMcs.size();
    int defSize  = (int)defMcs.size();

    bool consistent = validateTautomerConsistency(a.graph, b.graph, tautMcs);
    double tautConf = computeTautConfScore(a.graph, b.graph, tautMcs);

    return {a.name, b.name, tautSize, defSize,
            tautSize - defSize, consistent, tautConf,
            tautTimeMs, defTimeMs};
}

// ============================================================================
// Main
// ============================================================================
int main(int argc, char** argv) {
    std::string molPath = "benchmarks/diverse_molecules.txt";
    if (argc > 1) molPath = argv[1];

    fprintf(stderr, "Loading tautomer molecules from %s ...\n", molPath.c_str());
    auto mols = loadTautomerMolecules(molPath);
    fprintf(stderr, "  %zu valid tautomer molecules loaded\n", mols.size());

    if (mols.size() < 2) {
        fprintf(stderr, "ERROR: Need at least 2 molecules. Aborting.\n");
        return 1;
    }

    // Form consecutive pairs
    std::vector<std::pair<int,int>> pairs;
    for (int i = 0; i + 1 < (int)mols.size(); i += 2) {
        pairs.push_back({i, i + 1});
    }
    fprintf(stderr, "  %zu tautomer pairs formed\n", pairs.size());
    fprintf(stderr, "Running tautomer over-matching benchmark ...\n");

    std::vector<PairResult> results;
    for (int i = 0; i < (int)pairs.size(); i++) {
        auto [ai, bi] = pairs[i];
        PairResult pr = benchmarkPair(mols[ai], mols[bi]);
        results.push_back(pr);
        fprintf(stderr, "  [%2d/%zu] %-40s taut=%2d def=%2d delta=%+d %s tautConf=%.3f\n",
            i + 1, pairs.size(),
            (pr.nameA + " / " + pr.nameB).substr(0, 40).c_str(),
            pr.tautMcsSize, pr.defaultMcsSize, pr.overMatchDelta,
            pr.protonConsistent ? "PASS" : "FAIL", pr.tautConfScore);
    }

    // ================================================================
    // Report
    // ================================================================
    printf("\n");
    printf("==============================================================================\n");
    printf("SMSD Pro 6.0.0 — ZINC20 Tautomer Over-Matching Benchmark Results (C++)\n");
    printf("==============================================================================\n");
    printf("Pairs tested:              %zu\n", results.size());

    printf("\n%-42s %5s %5s %6s %6s %8s\n",
        "Pair", "Taut", "Def", "Delta", "Valid", "TautConf");
    printf("------------------------------------------------------------------------------\n");

    int overMatchCount = 0, failCount = 0;
    double totalTautConf = 0.0, totalTautTime = 0.0, totalDefTime = 0.0;

    for (auto& pr : results) {
        std::string pair = pr.nameA + " / " + pr.nameB;
        if (pair.size() > 42) pair = pair.substr(0, 39) + "...";
        printf("%-42s %5d %5d %+6d %6s %8.3f\n",
            pair.c_str(), pr.tautMcsSize, pr.defaultMcsSize, pr.overMatchDelta,
            pr.protonConsistent ? "PASS" : "FAIL", pr.tautConfScore);

        if (pr.overMatchDelta > 0) overMatchCount++;
        if (!pr.protonConsistent) failCount++;
        totalTautConf += pr.tautConfScore;
        totalTautTime += pr.tautTimeMs;
        totalDefTime  += pr.defaultTimeMs;
    }

    printf("------------------------------------------------------------------------------\n");
    int n = (int)results.size();
    double avgTautConf = n > 0 ? totalTautConf / n : 0.0;
    double passRate    = n > 0 ? 100.0 * (n - failCount) / n : 0.0;
    double overMatchPct= n > 0 ? 100.0 * overMatchCount / n : 0.0;

    printf("\nSummary Statistics:\n");
    printf("  Total pairs:                     %d\n", n);
    printf("  Pairs with over-matching:        %d (%.1f%%)\n", overMatchCount, overMatchPct);
    printf("  Proton consistency PASS rate:    %.1f%% (%d/%d)\n", passRate, n - failCount, n);
    printf("  Proton consistency FAIL count:   %d\n", failCount);
    printf("  Average TautConf score:          %.4f\n", avgTautConf);
    printf("  Total tautomer MCS time:         %.1f ms\n", totalTautTime);
    printf("  Total default MCS time:          %.1f ms\n", totalDefTime);
    printf("  Avg tautomer MCS time per pair:  %.2f ms\n", n > 0 ? totalTautTime / n : 0.0);
    printf("  Avg default MCS time per pair:   %.2f ms\n", n > 0 ? totalDefTime / n : 0.0);

    // Over-matched pairs detail
    std::vector<PairResult*> overMatched;
    for (auto& pr : results) {
        if (pr.overMatchDelta > 0) overMatched.push_back(&pr);
    }
    std::sort(overMatched.begin(), overMatched.end(),
        [](const PairResult* a, const PairResult* b) {
            return a->overMatchDelta > b->overMatchDelta;
        });

    if (!overMatched.empty()) {
        printf("\nOver-matched pairs (tautomer MCS > default MCS):\n");
        for (auto* pr : overMatched) {
            printf("  %s / %s: delta=%+d (taut=%d, def=%d) %s tautConf=%.3f\n",
                pr->nameA.c_str(), pr->nameB.c_str(), pr->overMatchDelta,
                pr->tautMcsSize, pr->defaultMcsSize,
                pr->protonConsistent ? "PASS" : "FAIL", pr->tautConfScore);
        }
    }

    // False positives
    int fpCount = 0;
    printf("\nFalse positives (over-matched + proton violation):\n");
    for (auto& pr : results) {
        if (pr.overMatchDelta > 0 && !pr.protonConsistent) {
            printf("  %s / %s: delta=%+d tautConf=%.3f\n",
                pr.nameA.c_str(), pr.nameB.c_str(),
                pr.overMatchDelta, pr.tautConfScore);
            fpCount++;
        }
    }
    if (fpCount == 0) printf("  (none)\n");

    printf("==============================================================================\n");

    return 0;
}

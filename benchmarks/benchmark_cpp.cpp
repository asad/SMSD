// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2018-2026 BioInception PVT LTD
// Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
// See the NOTICE file for attribution, trademark, and algorithm IP terms.
/**
 * Benchmark SMSD C++ (MolGraph.Builder equivalent) on diverse molecule pairs.
 *
 * This benchmark uses the MolGraph-style adjacency representation to
 * demonstrate CDK-free molecule construction and MCS timing. Since the
 * SMSD C++ header-only library mirrors the Java MolGraph.Builder API,
 * molecules are built from atomic numbers, neighbor lists, and bond orders.
 *
 * Compilation:
 *   g++ -std=c++17 -O2 -o benchmark_cpp benchmark_cpp.cpp
 *   ./benchmark_cpp
 *
 * Note: This is a standalone benchmark that builds molecular graphs from
 * primitive arrays (no SMILES parser required). It measures graph construction
 * and basic subgraph isomorphism timing using a simplified VF2 implementation.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>

// ============================================================================
// Minimal MolGraph representation (mirrors Java MolGraph.Builder)
// ============================================================================

struct MolGraph {
    int n;                              // atom count
    std::vector<int> atomicNum;         // atomic number per atom
    std::vector<bool> aromatic;         // aromaticity flag per atom
    std::vector<bool> inRing;           // ring membership per atom
    std::vector<std::vector<int>> adj;  // adjacency lists
    std::vector<std::vector<int>> bondOrder; // bond order per neighbor

    struct Builder {
        int n_ = 0;
        std::vector<int> atomicNum_;
        std::vector<bool> aromatic_;
        std::vector<bool> inRing_;
        std::vector<std::vector<int>> adj_;
        std::vector<std::vector<int>> bondOrder_;

        Builder& atomCount(int n) { n_ = n; return *this; }
        Builder& atomicNumbers(std::vector<int> a) { atomicNum_ = std::move(a); return *this; }
        Builder& aromaticFlags(std::vector<bool> a) { aromatic_ = std::move(a); return *this; }
        Builder& ringFlags(std::vector<bool> r) { inRing_ = std::move(r); return *this; }
        Builder& neighbors(std::vector<std::vector<int>> nb) { adj_ = std::move(nb); return *this; }
        Builder& bondOrders(std::vector<std::vector<int>> bo) { bondOrder_ = std::move(bo); return *this; }

        MolGraph build() {
            MolGraph g;
            g.n = n_;
            g.atomicNum = atomicNum_.empty() ? std::vector<int>(n_, 6) : atomicNum_;
            g.aromatic = aromatic_.empty() ? std::vector<bool>(n_, false) : aromatic_;
            g.inRing = inRing_.empty() ? std::vector<bool>(n_, false) : inRing_;
            g.adj = adj_;
            g.bondOrder = bondOrder_;
            return g;
        }
    };
};

// ============================================================================
// Simplified VF2-like MCS (backtracking with degree + label pruning)
// ============================================================================

struct McsResult {
    int size;
    std::vector<std::pair<int, int>> mapping;
};

static bool atomCompat(const MolGraph& g1, int i, const MolGraph& g2, int j) {
    return g1.atomicNum[i] == g2.atomicNum[j]
        && g1.aromatic[i] == g2.aromatic[j];
}

static void mcsBacktrack(
    const MolGraph& g1, const MolGraph& g2,
    std::vector<int>& map1, std::vector<int>& map2,
    int depth, McsResult& best, long long& nodes, long long nodeLimit)
{
    if (nodes++ > nodeLimit) return;
    if (depth > best.size) {
        best.size = depth;
        best.mapping.clear();
        for (int i = 0; i < g1.n; i++)
            if (map1[i] >= 0) best.mapping.push_back({i, map1[i]});
    }
    // Upper bound pruning
    int remain1 = 0;
    for (int i = 0; i < g1.n; i++) if (map1[i] < 0) remain1++;
    if (depth + remain1 <= best.size) return;

    for (int i = 0; i < g1.n; i++) {
        if (map1[i] >= 0) continue;
        for (int j = 0; j < g2.n; j++) {
            if (map2[j] >= 0) continue;
            if (!atomCompat(g1, i, g2, j)) continue;

            // Check neighbor consistency
            bool consistent = true;
            for (int ni : g1.adj[i]) {
                if (map1[ni] < 0) continue;
                bool found = false;
                for (int nj : g2.adj[j]) {
                    if (nj == map1[ni]) { found = true; break; }
                }
                if (!found) { consistent = false; break; }
            }
            if (!consistent) continue;

            map1[i] = j; map2[j] = i;
            mcsBacktrack(g1, g2, map1, map2, depth + 1, best, nodes, nodeLimit);
            map1[i] = -1; map2[j] = -1;

            if (nodes > nodeLimit) return;
        }
    }
}

static McsResult findMCS(const MolGraph& g1, const MolGraph& g2, long long nodeLimit = 500000) {
    McsResult best{0, {}};
    std::vector<int> map1(g1.n, -1), map2(g2.n, -1);
    long long nodes = 0;
    mcsBacktrack(g1, g2, map1, map2, 0, best, nodes, nodeLimit);
    return best;
}

// ============================================================================
// Molecule builders (equivalent to building from SMILES via MolGraph.Builder)
// ============================================================================

// Helper: build a linear alkane chain of length n
static MolGraph buildAlkane(int n) {
    MolGraph::Builder b;
    b.atomCount(n);
    std::vector<int> atoms(n, 6); // all carbon
    std::vector<bool> arom(n, false), ring(n, false);
    std::vector<std::vector<int>> adj(n), bo(n);
    for (int i = 0; i < n - 1; i++) {
        adj[i].push_back(i + 1);
        adj[i + 1].push_back(i);
        bo[i].push_back(1);
        bo[i + 1].push_back(1);
    }
    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// Build benzene (6 aromatic C)
static MolGraph buildBenzene() {
    MolGraph::Builder b;
    b.atomCount(6);
    std::vector<int> atoms(6, 6);
    std::vector<bool> arom(6, true), ring(6, true);
    std::vector<std::vector<int>> adj = {
        {1, 5}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}
    };
    std::vector<std::vector<int>> bo = {
        {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}
    };
    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// Build phenol (benzene + O at position 0)
static MolGraph buildPhenol() {
    MolGraph::Builder b;
    b.atomCount(7);
    std::vector<int> atoms = {6, 6, 6, 6, 6, 6, 8};
    std::vector<bool> arom = {true, true, true, true, true, true, false};
    std::vector<bool> ring = {true, true, true, true, true, true, false};
    std::vector<std::vector<int>> adj = {
        {1, 5, 6}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {0}
    };
    std::vector<std::vector<int>> bo = {
        {4, 4, 1}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {1}
    };
    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// Build toluene (benzene + CH3 at position 0)
static MolGraph buildToluene() {
    MolGraph::Builder b;
    b.atomCount(7);
    std::vector<int> atoms = {6, 6, 6, 6, 6, 6, 6};
    std::vector<bool> arom = {true, true, true, true, true, true, false};
    std::vector<bool> ring = {true, true, true, true, true, true, false};
    std::vector<std::vector<int>> adj = {
        {1, 5, 6}, {0, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 0}, {0}
    };
    std::vector<std::vector<int>> bo = {
        {4, 4, 1}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {4, 4}, {1}
    };
    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// Build adamantane (C10H16): 10 carbons, no aromaticity, all in ring
static MolGraph buildAdamantane() {
    MolGraph::Builder b;
    b.atomCount(10);
    std::vector<int> atoms(10, 6);
    std::vector<bool> arom(10, false), ring(10, true);
    // Adamantane adjacency: cage structure
    // Vertices: 0-3 bridgehead (degree 3), 4-9 methylene (degree 2)
    std::vector<std::vector<int>> adj = {
        {1, 3, 5},    // 0
        {0, 2, 6},    // 1
        {1, 3, 7},    // 2
        {0, 2, 8},    // 3
        {5, 6, 7},    // 4 (== bridgehead)
        {0, 4, 9},    // 5
        {1, 4, 9},    // 6
        {2, 4, 8},    // 7 (corrected)
        {3, 7, 9},    // 8
        {5, 6, 8},    // 9
    };
    std::vector<std::vector<int>> bo(10);
    for (int i = 0; i < 10; i++) bo[i].assign(adj[i].size(), 1);
    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// Build a PEG-like chain: O-CH2-CH2-O repeating
static MolGraph buildPEG(int repeats) {
    // Each repeat: O-C-C, then trailing O
    int n = repeats * 3 + 1; // 3 atoms per repeat + terminal O
    MolGraph::Builder b;
    b.atomCount(n);
    std::vector<int> atoms(n);
    std::vector<bool> arom(n, false), ring(n, false);
    std::vector<std::vector<int>> adj(n), bo(n);

    for (int r = 0; r < repeats; r++) {
        int base = r * 3;
        atoms[base] = 8;      // O
        atoms[base + 1] = 6;  // C
        atoms[base + 2] = 6;  // C
    }
    atoms[n - 1] = 8; // terminal O

    for (int i = 0; i < n - 1; i++) {
        adj[i].push_back(i + 1);
        adj[i + 1].push_back(i);
        bo[i].push_back(1);
        bo[i + 1].push_back(1);
    }

    return b.atomicNumbers(atoms).aromaticFlags(arom).ringFlags(ring)
            .neighbors(adj).bondOrders(bo).build();
}

// ============================================================================
// Benchmark pairs
// ============================================================================

struct BenchPair {
    std::string name;
    MolGraph (*build1)();
    MolGraph (*build2)();
};

static MolGraph buildMethane()   { return buildAlkane(1); }
static MolGraph buildEthane()    { return buildAlkane(2); }
static MolGraph buildWater() {
    MolGraph::Builder b;
    b.atomCount(1);
    return b.atomicNumbers({8}).aromaticFlags({false}).ringFlags({false})
            .neighbors({{}}).bondOrders({{}}).build();
}
static MolGraph buildMethanol() {
    MolGraph::Builder b;
    b.atomCount(2);
    return b.atomicNumbers({6, 8}).aromaticFlags({false, false}).ringFlags({false, false})
            .neighbors({{1}, {0}}).bondOrders({{1}, {1}}).build();
}
static MolGraph buildPEG12() { return buildPEG(12); }
static MolGraph buildPEG16() { return buildPEG(16); }

// ============================================================================
// Timing harness
// ============================================================================

static const int NUM_RUNS = 5;

struct TimingResult {
    double bestMs;
    double medianMs;
    int mcsSize;
};

template <typename Func>
TimingResult benchmark(Func fn) {
    // Warmup
    fn();

    std::vector<double> times;
    int mcsSize = 0;
    for (int i = 0; i < NUM_RUNS; i++) {
        auto t0 = std::chrono::high_resolution_clock::now();
        auto result = fn();
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        times.push_back(ms);
        mcsSize = result.size;
    }
    std::sort(times.begin(), times.end());
    return {times[0], times[NUM_RUNS / 2], mcsSize};
}

int main() {
    printf("%-30s  %10s  %10s  %4s\n", "Pair", "Best(ms)", "Median(ms)", "MCS");
    printf("--------------------------------------------------------------\n");

    double total = 0.0;

    // Pair 1: methane-ethane
    {
        auto g1 = buildMethane(); auto g2 = buildEthane();
        auto r = benchmark([&]() { return findMCS(g1, g2); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "methane-ethane", r.bestMs, r.medianMs, r.mcsSize);
    }

    // Pair 2: water-methanol
    {
        auto g1 = buildWater(); auto g2 = buildMethanol();
        auto r = benchmark([&]() { return findMCS(g1, g2); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "water-methanol", r.bestMs, r.medianMs, r.mcsSize);
    }

    // Pair 3: benzene-phenol
    {
        auto g1 = buildBenzene(); auto g2 = buildPhenol();
        auto r = benchmark([&]() { return findMCS(g1, g2); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "benzene-phenol", r.bestMs, r.medianMs, r.mcsSize);
    }

    // Pair 4: benzene-toluene
    {
        auto g1 = buildBenzene(); auto g2 = buildToluene();
        auto r = benchmark([&]() { return findMCS(g1, g2); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "benzene-toluene", r.bestMs, r.medianMs, r.mcsSize);
    }

    // Pair 5: adamantane-self
    {
        auto g1 = buildAdamantane();
        auto r = benchmark([&]() { return findMCS(g1, g1); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "adamantane-self", r.bestMs, r.medianMs, r.mcsSize);
    }

    // Pair 6: PEG12-PEG16
    {
        auto g1 = buildPEG12(); auto g2 = buildPEG16();
        auto r = benchmark([&]() { return findMCS(g1, g2); });
        total += r.bestMs;
        printf("%-30s  %10.3f  %10.3f  %4d\n", "PEG12-PEG16", r.bestMs, r.medianMs, r.mcsSize);
    }

    printf("--------------------------------------------------------------\n");
    printf("%-30s  %10.3f ms\n", "TOTAL", total);
    printf("\nNote: C++ benchmark covers a subset of pairs (those easily\n");
    printf("built without a SMILES parser). For full comparison, use the\n");
    printf("Java and Python benchmarks which parse SMILES directly.\n");

    return 0;
}

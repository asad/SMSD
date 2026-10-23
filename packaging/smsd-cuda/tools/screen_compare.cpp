/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Compares the CUDA screen of the C++ library (smsd::gpu::batchSimilarityScreenAuto, built with
 * SMSD_BUILD_CUDA) with the CPU screen over the molecules of SMILES files repeated to each size:
 * the hits must be the same and each GPU score must equal the CPU score bit for bit. Writes one
 * TSV row per screen with the median times.
 *
 *   screen_compare --sizes 10000,100000,1000000 --repeats 3 --out cpp.tsv a.smi b.txt
 */
#include "smsd/gpu.hpp"
#include "smsd/similarity_screen.hpp"
#include "smsd/smiles_parser.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<smsd::MolGraph> readMolecules(const std::vector<std::string>& paths, int& rejected) {
    std::vector<smsd::MolGraph> molecules;
    for (const auto& path : paths) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot read " + path);
        std::string line;
        while (std::getline(input, line)) {
            std::istringstream fields(line);
            std::string smiles;
            if (!(fields >> smiles) || smiles[0] == '#') continue;
            try { molecules.push_back(smsd::parseSMILES(smiles)); } catch (const std::exception&) { ++rejected; }
        }
    }
    return molecules;
}

template<class F>
double medianMs(int repeats, F&& run) {
    std::vector<double> times;
    for (int r = 0; r < repeats; ++r) {
        const auto start = std::chrono::steady_clock::now();
        run();
        times.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
    }
    std::sort(times.begin(), times.end());
    return times[times.size() / 2];
}

std::vector<int> splitSizes(const std::string& text) {
    std::vector<int> sizes;
    std::istringstream input(text);
    std::string value;
    while (std::getline(input, value, ',')) sizes.push_back(std::stoi(value));
    return sizes;
}

} // namespace

int main(int argc, char** argv) {
    std::string sizesText = "10000,100000,1000000", out = "screen_compare_cpp.tsv";
    int repeats = 3, queryCount = 4;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--sizes" && i + 1 < argc) sizesText = argv[++i];
        else if (arg == "--repeats" && i + 1 < argc) repeats = std::max(1, std::atoi(argv[++i]));
        else if (arg == "--queries" && i + 1 < argc) queryCount = std::max(1, std::atoi(argv[++i]));
        else if (arg == "--out" && i + 1 < argc) out = argv[++i];
        else files.push_back(arg);
    }
    if (files.empty()) { std::cerr << "usage: screen_compare [--sizes a,b] [--repeats n] [--out file] smiles...\n"; return 2; }
    if (!smsd::gpu::isAvailable()) { std::cerr << "no CUDA device: " << smsd::gpu::deviceInfo() << '\n'; return 2; }
    std::cout << smsd::gpu::deviceInfo() << '\n';

    int rejected = 0;
    const std::vector<smsd::MolGraph> pool = readMolecules(files, rejected);
    if (pool.empty()) { std::cerr << "no molecules were read\n"; return 2; }
    std::cout << pool.size() << " molecules, " << rejected << " lines rejected\n";
    std::vector<const smsd::MolGraph*> queries;
    for (int q = 0; q < queryCount; ++q) queries.push_back(&pool[(static_cast<size_t>(q) * pool.size()) / queryCount]);
    const smsd::ChemOptions chemistry;
    const double thresholds[] = {0.0, 0.5, 0.8, 1.0};

    std::ofstream table(out);
    table << "size\tquery\tthreshold\tcpu_hits\tgpu_hits\tequal_hits\tequal_scores\tcpu_ms\tgpu_ms\tcpu_over_gpu\n";
    long mismatches = 0;
    for (const int size : splitSizes(sizesText)) {
        std::vector<smsd::MolGraph> targets;
        targets.reserve(size);
        for (int i = 0; i < size; ++i) targets.push_back(pool[static_cast<size_t>(i) % pool.size()]);
        for (size_t q = 0; q < queries.size(); ++q) {
            for (const double threshold : thresholds) {
                std::vector<double> scores;
                std::vector<smsd::gpu::ScreenHit> gpu;
                const double cpuMs = medianMs(repeats, [&] {
                    scores = smsd::screening_detail::similarityScores(*queries[q], targets, chemistry, 0);
                });
                const std::vector<int> cpu = smsd::screening_detail::indicesReaching(scores, threshold);
                const double gpuMs = medianMs(repeats, [&] {
                    gpu = smsd::gpu::batchSimilarityScreenAuto(*queries[q], targets, chemistry, threshold);
                });
                bool sameHits = gpu.size() == cpu.size(), sameScores = sameHits;
                for (size_t k = 0; sameHits && k < gpu.size(); ++k) {
                    sameHits = gpu[k].index == cpu[k];
                    sameScores = sameScores && sameHits
                        && std::memcmp(&gpu[k].score, &scores[cpu[k]], sizeof(double)) == 0;
                }
                mismatches += !(sameHits && sameScores);
                table << size << '\t' << q << '\t' << threshold << '\t' << cpu.size() << '\t' << gpu.size() << '\t'
                      << sameHits << '\t' << sameScores << '\t' << cpuMs << '\t' << gpuMs << '\t'
                      << (gpuMs > 0 ? cpuMs / gpuMs : 0.0) << '\n';
            }
        }
        table.flush();
        std::cout << "size " << size << " done, mismatches so far " << mismatches << '\n';
    }
    std::cout << "screen_compare: " << mismatches << " mismatches\n";
    return mismatches == 0 ? 0 : 1;
}

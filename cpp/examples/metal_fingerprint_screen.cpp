/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Packed ECFP4 fingerprint screening on a Metal GPU, with the CPU screen as the fallback.
 *
 * On macOS, build it from the source tree with Metal enabled:
 *
 *   cmake -S cpp -B build/metal-example -DSMSD_BUILD_METAL=ON -DSMSD_BUILD_EXAMPLES=ON -DSMSD_BUILD_TESTS=OFF
 *   cmake --build build/metal-example --target smsd_metal_fingerprint_screen
 *   build/metal-example/smsd_metal_fingerprint_screen
 *
 * Linking smsd_metal defines SMSD_ENABLE_METAL. Without it, on any platform, the same file
 * builds against the headers alone and runs the CPU screen:
 *
 *   c++ -std=c++17 -I cpp/include cpp/examples/metal_fingerprint_screen.cpp -o metal_fingerprint_screen
 *
 * Both paths give the same hits, in query order and, for each query, in target order.
 * Set SMSD_DISABLE_GPU=1 to keep the screen on the CPU in a Metal build.
 */
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "smsd/smsd.hpp"
#include "smsd/screen_accelerator.hpp"
#ifdef SMSD_ENABLE_METAL
#include "smsd/metal_gpu.hpp"
#endif

using Fingerprints = std::vector<std::vector<std::uint64_t>>;
using Hits = std::vector<std::vector<int>>;

/** The CPU screen: for each query, the targets whose Tanimoto similarity reaches the threshold. */
static Hits cpuScreen(const Fingerprints& queries, const Fingerprints& targets, double threshold) {
    Hits hits(queries.size());
    for (std::size_t query = 0; query < queries.size(); ++query)
        for (std::size_t target = 0; target < targets.size(); ++target)
            if (smsd::batch::fingerprintTanimoto(queries[query], targets[target]) >= threshold)
                hits[query].push_back(static_cast<int>(target));
    return hits;
}

int main() {
    try {
        const std::vector<std::string> querySmiles{
            "CC(=O)Oc1ccccc1C(=O)O",   // aspirin
            "Oc1ccccc1",               // phenol
        };
        const std::vector<std::string> targetSmiles{
            "OC(=O)c1ccccc1O",         // salicylic acid
            "CC(=O)Oc1ccccc1",         // phenyl acetate
            "c1ccccc1",                // benzene
            "Oc1ccc(O)cc1",            // hydroquinone
            "CC(=O)Nc1ccc(O)cc1",      // paracetamol
        };
        const double threshold = 0.3;

        // Every query and target uses the same settings, so all fingerprints have 32 words.
        const auto options = smsd::batch::FingerprintOptions::ecfp(2, 2048);
        Fingerprints queries, targets;
        for (const auto& smiles : querySmiles) queries.push_back(smsd::batch::fingerprintFromSmiles(smiles, options));
        for (const auto& smiles : targetSmiles) targets.push_back(smsd::batch::fingerprintFromSmiles(smiles, options));

        std::optional<Hits> hits;
        std::string path = "CPU";
        std::string reason = "this build does not link smsd_metal";

#ifdef SMSD_ENABLE_METAL
        // The same Tanimoto expression as the CPU screen, formed from the bit counts.
        const auto score = [](int common, int queryBits, int targetBits) {
            return smsd::batch::tanimoto(smsd::batch::BitCounts{common, queryBits, targetBits});
        };
        // Keeps the target fingerprints on the GPU for later screens of the same targets.
        smsd::screen_accelerator::ResidentCache cache;
        const long long failuresBefore = smsd::screen_accelerator::failureCount();
        hits = smsd::screen_accelerator::fingerprintScreens(
            queries, targets, threshold, cache, score, smsd::metal::fingerprintAccelerator());
        if (hits) {
            path = smsd::metal::metalDeviceInfo();   // for example "Metal GPU: Apple M2 Pro"
            reason.clear();
        } else if (smsd::screen_accelerator::failureCount() > failuresBefore) {
            reason = smsd::screen_accelerator::latestFailure();
        } else {
            reason = "the screen was declined";
        }
#endif

        // std::nullopt means the GPU did not screen; the caller runs its CPU screen instead.
        const Hits cpu = cpuScreen(queries, targets, threshold);
        if (!hits) hits = cpu;

        std::cout << "Screen ran on: " << path;
        if (!reason.empty()) std::cout << " (" << reason << ")";
        std::cout << "\n";
        std::cout << "Tanimoto threshold: " << threshold << "\n";
        for (std::size_t query = 0; query < queries.size(); ++query) {
            std::cout << "Query " << query << "  " << querySmiles[query] << "\n";
            if ((*hits)[query].empty()) std::cout << "  no hits\n";
            for (int hit : (*hits)[query]) {
                const auto target = static_cast<std::size_t>(hit);
                std::cout << "  target " << target << "  Tanimoto " << std::fixed << std::setprecision(3)
                          << smsd::batch::fingerprintTanimoto(queries[query], targets[target])
                          << std::defaultfloat << "  " << targetSmiles[target] << "\n";
            }
        }

        if (*hits != cpu) {
            std::cerr << "The screen gave different hits from the CPU screen\n";
            return 1;
        }
        std::cout << "The hits are identical to those of the CPU screen.\n";
        return 0;
    } catch (const std::exception& failure) {
        std::cerr << failure.what() << '\n';
        return 1;
    }
}

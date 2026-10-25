/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#define SMSD_ENABLE_METAL
#include "test_support.hpp"
#include "smsd/gpu.hpp"
#include "../metal/metal_batch.mm"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>

int main() {
    try {
        if (!smsd::metal::metalAvailable()) return 77;
        @autoreleasepool {
            auto& context = MetalContext::instance();
            auto& probe = DeviceProbe::instance();
            smsd_test::require(context.valid && probe.found, "Physical Metal pipeline must compile");
            smsd::metal::MetalMolecule query{}, target{};
            query.n = 5; query.labelCounts[0] = 1; query.labelCounts[1] = 4;
            target.n = 6; target.labelCounts[0] = 1; target.labelCounts[2] = 5;
            const std::vector<smsd::metal::MetalMolecule> targets{target, target, query, {}};
            id<MTLBuffer> queryBuffer = [probe.device newBufferWithBytes:&query length:sizeof query options:MTLResourceStorageModeShared];
            id<MTLBuffer> targetBuffer = [probe.device newBufferWithBytes:targets.data() length:targets.size() * sizeof target options:MTLResourceStorageModeShared];
            id<MTLBuffer> resultBuffer = [probe.device newBufferWithLength:targets.size() * sizeof(int) options:MTLResourceStorageModeShared];
            smsd_test::require(queryBuffer && targetBuffer && resultBuffer, "Physical buffers must be created");
            id<MTLCommandBuffer> command = [context.queue commandBuffer];
            id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
            smsd_test::require(command && encoder, "Physical command must be created");
            [encoder setComputePipelineState:context.pipeline];
            [encoder setBuffer:targetBuffer offset:0 atIndex:0];
            [encoder setBuffer:queryBuffer offset:0 atIndex:1];
            [encoder setBuffer:resultBuffer offset:0 atIndex:2];
            const uint32_t count = static_cast<uint32_t>(targets.size());
            [encoder setBytes:&count length:sizeof count atIndex:3];
            const NSUInteger width = std::min<NSUInteger>(count, context.pipeline.maxTotalThreadsPerThreadgroup);
            smsd_test::require(width != 0, "Physical thread group must be nonempty");
            [encoder dispatchThreadgroups:MTLSizeMake((count + width - 1) / width, 1, 1)
                    threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
            [encoder endEncoding];
            [command commit];
            [command waitUntilCompleted];
            smsd_test::require(command.status == MTLCommandBufferStatusCompleted, "Physical command must complete without fallback");
            const int expected[]{1, 1, 5, 0};
            smsd_test::require(resultBuffer.contents && std::memcmp(resultBuffer.contents, expected, sizeof expected) == 0,
                               "Physical shader must return exact integer overlaps");

            const auto floatHits = smsd::metal::batchSimilarityScreen(query, {target}, 0.1f);
            smsd_test::require(floatHits.size() == 1 && floatHits[0].score == 1.0f / 10.0f,
                               "Direct Metal float threshold semantics must stay unchanged");
            smsd::metal::MetalMolecule largeQuery{}, largeTarget{};
            largeQuery.n = (1 << 24) + 1; largeQuery.labelCounts[0] = largeQuery.n;
            largeTarget.n = (1 << 24) + 3; largeTarget.labelCounts[0] = largeTarget.n;
            const auto largeHits = smsd::metal::batchSimilarityScreen(largeQuery, {largeTarget}, 0.0f);
            smsd_test::require(largeHits.size() == 1 && largeHits[0].score ==
                static_cast<float>(largeQuery.n) / static_cast<float>(largeTarget.n),
                "Direct Metal scores must retain float operand division");

            const auto graph = [](const std::vector<int>& elements) {
                return smsd::MolGraph::Builder().atomCount(static_cast<int>(elements.size()))
                    .atomicNumbers(elements).setNeighbors(std::vector<std::vector<int>>(elements.size())).build(false);
            };
            const auto graphQuery = graph({6, 7, 7, 7, 7});
            const auto graphTarget = graph({6, 8, 8, 8, 8, 8});
            const std::vector<smsd::MolGraph> graphTargets{graphTarget, graphTarget, graphQuery, graph({})};
            const auto scores = smsd::screening_detail::similarityScores(graphQuery, graphTargets, smsd::ChemOptions(), 1);
            for (double threshold : {std::nextafter(scores[0], 0.0), scores[0], std::nextafter(scores[0], 1.0),
                    0.0, -0.0, 1.0, -std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
                const auto indices = smsd::screening_detail::indicesReaching(scores, threshold);
                const auto hits = smsd::gpu::batchSimilarityScreenAuto(graphQuery, graphTargets, threshold);
                smsd_test::require(hits.size() == indices.size(), "Metal threshold membership must follow host doubles");
                for (size_t hit = 0; hit < hits.size(); ++hit)
                    smsd_test::require(hits[hit].index == indices[hit] &&
                        std::memcmp(&hits[hit].score, &scores[indices[hit]], sizeof(double)) == 0,
                        "Metal ordered hit scores must equal CPU scores bit for bit");
            }
            std::cout << "Physical Metal integer overlaps and host thresholds passed\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}

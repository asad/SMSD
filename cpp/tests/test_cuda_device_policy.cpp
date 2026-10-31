/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The CUDA device rules on every platform, with the device properties and environment values
 * given directly: the compute capability floor, SMSD_CUDA_MIN_TARGETS parsing, the Tesla T4
 * default, and SMSD_DISABLE_GPU overriding both.
 */
#include "test_support.hpp"
#include "smsd/cuda_device_policy.hpp"

#include <cstdint>
#include <string>

namespace {

using smsd_test::require;
namespace policy = smsd::cuda_policy;

constexpr std::uint64_t kNone = policy::kNoAutomaticUse;

void capabilityFloorIsSevenPointFive() {
    require(!policy::supportedCapability(6, 1), "compute 6.1 is not supported");
    require(!policy::supportedCapability(7, 0), "compute 7.0 is not supported");
    require(!policy::supportedCapability(7, 2), "compute 7.2 is not supported");
    require(policy::supportedCapability(7, 5), "compute 7.5 is supported");
    require(policy::supportedCapability(8, 0), "compute 8.0 is supported");
    require(policy::supportedCapability(12, 0), "a later major version is supported");

    char reason[256] = {};
    require(!policy::explainCapability("Tesla P100-PCIE-16GB", 6, 0, reason, sizeof reason),
            "an older device is refused");
    require(std::string(reason) == "the CUDA device Tesla P100-PCIE-16GB has compute capability 6.0; "
                                   "CUDA screening needs 7.5 or later", "the reason names the device and floor");
    char untouched[8] = "same";
    require(policy::explainCapability("Tesla T4", 7, 5, untouched, sizeof untouched)
                && std::string(untouched) == "same", "a supported device leaves the reason untouched");
    char shortReason[12] = {};
    require(!policy::explainCapability(nullptr, 5, 2, shortReason, sizeof shortReason)
                && std::string(shortReason) == "the CUDA de", "a short buffer is truncated and terminated");
    require(!policy::explainCapability("", 3, 5, nullptr, 0), "a null buffer is accepted");
}

void minimumTargetsParsing() {
    require(policy::requestedMinimumTargets(nullptr) == 0, "unset asks for nothing");
    require(policy::requestedMinimumTargets("") == 0, "empty asks for nothing");
    require(policy::requestedMinimumTargets("0") == 0, "zero asks for nothing");
    require(policy::requestedMinimumTargets("000") == 0, "zeros ask for nothing");
    require(policy::requestedMinimumTargets("1") == 1, "one is accepted");
    require(policy::requestedMinimumTargets("5000") == 5000, "a positive integer is accepted");
    require(policy::requestedMinimumTargets("0100") == 100, "leading zeros keep the decimal value");
    for (const char* invalid : {"-1", "+5", " 5", "5 ", "5k", "1e5", "1.0", "0x10", "abc", "5\n"})
        require(policy::requestedMinimumTargets(invalid) == 0,
                std::string("a non-decimal value is ignored: ") + invalid);
    require(policy::requestedMinimumTargets("18446744073709551615") == kNone, "the 64-bit maximum is parsed");
    require(policy::requestedMinimumTargets("18446744073709551616") == 0, "64-bit overflow is ignored");
    require(policy::requestedMinimumTargets("99999999999999999999999") == 0, "a long overflow is ignored");
}

void disableRuleMatchesTheLibrary() {
    require(!policy::disabled(nullptr), "unset leaves the GPU enabled");
    require(!policy::disabled(""), "empty leaves the GPU enabled");
    require(!policy::disabled("0"), "0 leaves the GPU enabled");
    for (const char* value : {"1", "yes", "true", "false", "no", "00", " "})
        require(policy::disabled(value), std::string("any other value disables the GPU: ") + value);
}

void automaticMinimumIsUnchangedWhenUnset() {
    require(policy::automaticMinimumTargets(1, "Tesla T4", 7, 5, nullptr, nullptr) == 100000,
            "a Tesla T4 keeps its crossover");
    require(policy::automaticMinimumTargets(1, "NVIDIA A100-SXM4-40GB", 8, 0, nullptr, nullptr) == kNone,
            "another device is not used automatically");
    require(policy::automaticMinimumTargets(1, "NVIDIA GeForce RTX 2080", 7, 5, nullptr, nullptr) == kNone,
            "another compute 7.5 device is not used automatically");
    require(policy::automaticMinimumTargets(1, nullptr, 7, 5, nullptr, nullptr) == kNone,
            "an unnamed device is not used automatically");
    for (std::uint32_t operation : {2u, 3u})
        require(policy::automaticMinimumTargets(operation, "Tesla T4", 7, 5, nullptr, nullptr) == kNone,
                "label and substructure screens stay on the CPU");
    for (const char* invalid : {"", "0", "-5", "5k", "18446744073709551616"}) {
        require(policy::automaticMinimumTargets(1, "Tesla T4", 7, 5, nullptr, invalid) == 100000,
                std::string("an ignored value keeps the Tesla T4 crossover: ") + invalid);
        require(policy::automaticMinimumTargets(1, "NVIDIA A100-SXM4-40GB", 8, 0, nullptr, invalid) == kNone,
                std::string("an ignored value does not enable another device: ") + invalid);
    }
}

void minimumTargetsExtendsSupportedDevices() {
    require(policy::automaticMinimumTargets(1, "NVIDIA A100-SXM4-40GB", 8, 0, nullptr, "5000") == 5000,
            "a supported device is used from the requested count");
    require(policy::automaticMinimumTargets(1, "NVIDIA GeForce RTX 2080", 7, 5, "0", "1") == 1,
            "SMSD_DISABLE_GPU=0 leaves the request in force");
    require(policy::automaticMinimumTargets(1, "Tesla T4", 7, 5, nullptr, "250000") == 250000,
            "the request replaces the Tesla T4 crossover");
    require(policy::automaticMinimumTargets(1, "Tesla T4", 7, 5, nullptr, "20") == 20,
            "the request may lower the Tesla T4 crossover");
    require(policy::automaticMinimumTargets(1, "Tesla V100-SXM2-16GB", 7, 0, nullptr, "1") == kNone,
            "an older device is never used");
    require(policy::automaticMinimumTargets(1, "Tesla P4", 6, 1, nullptr, "1") == kNone,
            "an older device is never used");
    for (std::uint32_t operation : {2u, 3u})
        require(policy::automaticMinimumTargets(operation, "NVIDIA A100-SXM4-40GB", 8, 0, nullptr, "1") == kNone,
                "the request covers fingerprint screening only");
}

void disableOverridesMinimumTargets() {
    for (const char* value : {"1", "yes", "false"}) {
        require(policy::automaticMinimumTargets(1, "NVIDIA A100-SXM4-40GB", 8, 0, value, "1") == kNone,
                std::string("SMSD_DISABLE_GPU overrides SMSD_CUDA_MIN_TARGETS: ") + value);
        require(policy::automaticMinimumTargets(1, "Tesla T4", 7, 5, value, nullptr) == kNone,
                std::string("SMSD_DISABLE_GPU overrides the Tesla T4 crossover: ") + value);
    }
}

} // namespace

int main(int argc, char** argv) {
    return smsd_test::runGroups(argc, argv, {
        {"capability-floor", capabilityFloorIsSevenPointFive},
        {"minimum-targets-parsing", minimumTargetsParsing},
        {"disable-rule", disableRuleMatchesTheLibrary},
        {"unset-is-unchanged", automaticMinimumIsUnchangedWhenUnset},
        {"minimum-targets", minimumTargetsExtendsSupportedDevices},
        {"disable-overrides", disableOverridesMinimumTargets},
    });
}

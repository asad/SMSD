/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * Which CUDA devices the screens use, and from how many targets automatic screening uses one.
 * Plain C++ with no CUDA types, so the CUDA sources apply these rules and every platform tests
 * them. A device below compute capability 7.5 is not used, and its screens run on the CPU.
 * Automatic fingerprint screening uses a Tesla T4 from 100,000 targets. SMSD_CUDA_MIN_TARGETS
 * set to a positive decimal integer N extends it to any supported device from N targets; any
 * other value is ignored. SMSD_DISABLE_GPU set to a value other than 0 overrides both.
 */
#pragma once
#ifndef SMSD_CUDA_DEVICE_POLICY_HPP
#define SMSD_CUDA_DEVICE_POLICY_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>

namespace smsd {
namespace cuda_policy {

/** The oldest compute capability the CUDA screens are built for. */
constexpr int kMinimumMajor = 7;
constexpr int kMinimumMinor = 5;

/** The automatic fingerprint-screening crossover measured on a Tesla T4. */
constexpr std::uint64_t kTeslaT4MinimumTargets = 100000;

/** The minimum that means automatic screening does not use the device. */
constexpr std::uint64_t kNoAutomaticUse = std::numeric_limits<std::uint64_t>::max();

/** Whether a device of this compute capability can run the screens. */
constexpr bool supportedCapability(int major, int minor) noexcept {
    return major > kMinimumMajor || (major == kMinimumMajor && minor >= kMinimumMinor);
}

/**
 * True when a device of this compute capability can run the screens. Otherwise false, with
 * the reason written to reason (truncated to size, always terminated when size is not zero).
 */
inline bool explainCapability(const char* name, int major, int minor, char* reason,
                              std::size_t size) noexcept {
    if (supportedCapability(major, minor)) return true;
    if (reason != nullptr && size != 0)
        std::snprintf(reason, size, "the CUDA device %s has compute capability %d.%d; CUDA screening "
                      "needs %d.%d or later", name != nullptr && name[0] != '\0' ? name : "in use",
                      major, minor, kMinimumMajor, kMinimumMinor);
    return false;
}

/** Whether the two texts are equal; a null text equals nothing. */
constexpr bool sameText(const char* first, const char* second) noexcept {
    if (first == nullptr || second == nullptr) return false;
    for (; *first != '\0' && *first == *second; ++first, ++second) {}
    return *first == *second;
}

/** Whether this value of SMSD_DISABLE_GPU (null when unset) keeps every screen on the CPU. */
constexpr bool disabled(const char* value) noexcept {
    return value != nullptr && value[0] != '\0' && !sameText(value, "0");
}

/**
 * The target count this value of SMSD_CUDA_MIN_TARGETS (null when unset) asks for, or 0 when
 * it asks for none: unset, empty, zero, anything but decimal digits (signs and spaces
 * included), or a number too large for 64 bits.
 */
constexpr std::uint64_t requestedMinimumTargets(const char* value) noexcept {
    if (value == nullptr || value[0] == '\0') return 0;
    std::uint64_t result = 0;
    for (const char* digit = value; *digit != '\0'; ++digit) {
        if (*digit < '0' || *digit > '9') return 0;
        const std::uint64_t next = static_cast<std::uint64_t>(*digit - '0');
        if (result > (kNoAutomaticUse - next) / 10) return 0;
        result = result * 10 + next;
    }
    return result;
}

/**
 * The fewest targets from which automatic screening of the operation (1 fingerprint, 2 label,
 * 3 substructure screening) uses this device, or kNoAutomaticUse. disableSetting and
 * minimumSetting are the values of SMSD_DISABLE_GPU and SMSD_CUDA_MIN_TARGETS, null when unset.
 */
constexpr std::uint64_t automaticMinimumTargets(std::uint32_t operation, const char* name, int major,
                                                int minor, const char* disableSetting,
                                                const char* minimumSetting) noexcept {
    if (operation != 1 || disabled(disableSetting) || !supportedCapability(major, minor))
        return kNoAutomaticUse;
    const std::uint64_t requested = requestedMinimumTargets(minimumSetting);
    if (requested != 0) return requested;
    if (major == 7 && minor == 5 && sameText(name, "Tesla T4"))
        return kTeslaT4MinimumTargets;
    return kNoAutomaticUse;
}

} // namespace cuda_policy
} // namespace smsd

#endif // SMSD_CUDA_DEVICE_POLICY_HPP

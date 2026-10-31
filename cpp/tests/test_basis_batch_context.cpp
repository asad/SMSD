/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 */
#include "smsd/fingerprint_options.hpp"

#include <bitset>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

namespace {
const char* basisCause =
    "Exact ring basis unavailable: ring perception reached its work limit before the minimum cycle basis was complete";
// A 32-layer hexagonal prism stack, whose exact basis the ring perception completes.
const char* stack32Smiles = "C12C3C4C5C6C1C17C28C32C43C54C61C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C17C28C32C43C64C51C15C46C34C23C82C71C1C2C3C4C6C51";
// An 80-layer stack: a synthetic cage whose ring system exceeds the ring-perception work limit.
const char* cageSmiles = "C12C3C4C5C6C1C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C17C68C56C45C34C21C12C43C54C65C86C71C1C6C5C4C3C21";

void require(bool condition, const char* message) {
    if (!condition) throw std::logic_error(message);
}

struct Failure {
    std::string type;
    std::string message;
    std::string cause;
    bool runtime;
};

template<class Call> Failure failureOf(const Call& call) {
    try {
        call();
    } catch (const std::exception& error) {
        Failure failure{typeid(error).name(), error.what(), {},
                        dynamic_cast<const std::runtime_error*>(&error) != nullptr};
        try {
            std::rethrow_if_nested(error);
        } catch (const std::exception& cause) {
            failure.cause = cause.what();
        }
        return failure;
    }
    throw std::logic_error("The fixture must raise an exception");
}

/** Readable graph sequence with a later distinguishable native parsing failure. */
struct Library {
    const smsd::MolGraph& ordinary;
    const smsd::MolGraph& incomplete;
    bool earlierInvalid = false;
    size_t size() const { return 32; }
    const smsd::MolGraph& operator[](size_t inputIndex) const {
        if (earlierInvalid && inputIndex == 3) {
            (void)smsd::parseSMILES("C(C");
            throw std::logic_error("The earlier malformed fixture must fail to parse");
        }
        if (inputIndex == 17) {
            (void)smsd::parseSMILES("C1CC");
            throw std::logic_error("The later malformed fixture must fail to parse");
        }
        return inputIndex == 9 ? incomplete : ordinary;
    }
};

template<class Molecules> void batch(const Molecules& molecules,
                                    const smsd::batch::FingerprintOptions& options,
                                    int threads, bool counts) {
    if (counts) {
        // The existing native primitive used by the Python count-batch binding.
        (void)smsd::batch::detail::fingerprintEach(molecules, threads,
            [&options](const smsd::MolGraph& molecule) {
                return smsd::batch::computeFingerprintCounts(molecule, options);
            });
    } else {
        (void)smsd::batch::batchFingerprint(molecules, options, threads);
    }
}

void basisContext(const smsd::MolGraph& ordinary, const smsd::MolGraph& incomplete) {
    auto options = smsd::batch::FingerprintOptions::ecfp();
    options.ringSize = true;
    (void)smsd::batch::computeFingerprint(ordinary, options);
    const auto earlier = failureOf([] { (void)smsd::parseSMILES("C(C"); });
    const auto later = failureOf([] { (void)smsd::parseSMILES("C1CC"); });
    require(earlier.message != later.message, "The parse failures must be distinguishable");
    for (bool counts : {false, true}) {
        const auto single = failureOf([&] {
            if (counts) (void)smsd::batch::computeFingerprintCounts(incomplete, options);
            else (void)smsd::batch::computeFingerprint(incomplete, options);
        });
        require(single.runtime && single.message == basisCause && single.cause.empty(),
                "The single-call runtime_error cause stays unchanged");
        for (int threads : {1, 2, 4, 0}) {
            const auto actual = failureOf([&] { batch(Library{ordinary, incomplete}, options, threads, counts); });
            std::cout << "basis counts=" << counts << " threads=" << threads
                      << " message=" << actual.message << " cause=" << actual.cause << std::endl;
            require(actual.runtime && actual.message == single.message + " (input index 9, zero-based)",
                    "The first unavailable-basis input must be identified");
            require(actual.cause == single.message, "The original single-call cause is retained");
            const auto invalid = failureOf([&] { batch(Library{ordinary, incomplete, true}, options, threads, counts); });
            require(invalid.type == earlier.type && invalid.message == earlier.message && invalid.cause.empty(),
                    "An earlier ordinary error keeps its exact type and message");
        }
    }
}

void otherRuntimeErrors(const smsd::MolGraph& ordinary, const smsd::MolGraph& incomplete) {
    const std::vector<smsd::MolGraph> molecules(32, ordinary);
    for (int threads : {1, 2, 4, 0}) {
        const auto actual = failureOf([&] {
            (void)smsd::batch::detail::fingerprintEach(molecules, threads,
                [](const smsd::MolGraph&) -> std::vector<int> {
                    throw std::runtime_error("ordinary callback failure");
                });
        });
        require(actual.message == "ordinary callback failure" && actual.cause.empty(),
                "An unrelated runtime_error must not acquire basis context");
    }
    // A key that counts ring memberships is refused with the ring-basis cause, and a batch names its input.
    const auto keys = smsd::batch::FingerprintOptions::withKeys({smsd::batch::FingerprintKey("[R2]")});
    const auto single = failureOf([&] { (void)smsd::batch::computeFingerprint(incomplete, keys); });
    require(single.runtime && single.message == basisCause && single.cause.empty(),
            "The ring-counting key names the ring-basis cause");
    for (bool counts : {false, true}) for (int threads : {1, 2, 4, 0}) {
        const auto actual = failureOf([&] { batch(Library{ordinary, incomplete}, keys, threads, counts); });
        require(actual.message == single.message + " (input index 9, zero-based)" && actual.cause == single.message,
                "A ring-caused key refusal identifies its input");
    }
}

/** The embedding limit settles a key whose minimum is met and refuses a count it left open. */
void embeddingLimit(const smsd::MolGraph& stack32) {
    const auto key = smsd::batch::FingerprintKey("C~C~C~C~C~C~C~C~C~C~C~C");
    const auto presence = smsd::batch::computeFingerprint(stack32, smsd::batch::FingerprintOptions::withKeys({key}));
    int bits = 0;
    for (auto word : presence) bits += static_cast<int>(std::bitset<64>(word).count());
    require(bits == 1, "A key whose minimum is met is present although the search stopped");
    const char* unsettled = "a substructure key matched more than 100000 embeddings before its count was settled";
    const auto counts = smsd::batch::FingerprintOptions::withKeys({key});
    const auto single = failureOf([&] { (void)smsd::batch::computeFingerprintCounts(stack32, counts); });
    require(single.runtime && single.message == unsettled && single.cause.empty(),
            "An unsettled count keeps the 8.0.3 message");
    const auto ordinary = smsd::parseSMILES("CCO");
    for (int threads : {1, 2, 4, 0}) {
        const auto actual = failureOf([&] {
            (void)smsd::batch::detail::fingerprintEach(Library{ordinary, stack32}, threads,
                [&counts](const smsd::MolGraph& molecule) {
                    return smsd::batch::computeFingerprintCounts(molecule, counts);
                });
        });
        require(actual.message == std::string(unsettled) + " (input index 9, zero-based)" && actual.cause == unsettled,
                "An unsettled key count identifies its input");
    }
}
}

int main() {
    const auto ordinary = smsd::parseSMILES("CCO");
    const auto incomplete = smsd::parseSMILES(cageSmiles);
    require(smsd::ringsTruncated(incomplete), "The cage genuinely lacks a complete basis");
    require(!smsd::ringsTruncated(ordinary), "The ordinary prefix has a complete basis");
    const auto stack32 = smsd::parseSMILES(stack32Smiles);
    require(!smsd::ringsTruncated(stack32), "The 32-layer stack has a complete exact basis");
    const auto& basis = smsd::detail_rf::completeBasis(stack32);
    size_t weight = 0;
    for (const auto& cycle : basis) weight += cycle.size();
    require(basis.size() == 187 && weight == 750, "The 32-layer stack basis has full rank and minimum weight");
    (void)smsd::batch::computeFingerprint(stack32, smsd::batch::FingerprintOptions::maccs());
    basisContext(ordinary, incomplete);
    otherRuntimeErrors(ordinary, incomplete);
    embeddingLimit(stack32);
#ifdef _OPENMP
    std::cout << "basis_batch_context_controls_passed openmp=1" << std::endl;
#else
    std::cout << "basis_batch_context_controls_passed openmp=0" << std::endl;
#endif
}

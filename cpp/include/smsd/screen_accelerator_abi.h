/*
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) 2018-2026 BioInception PVT LTD
 * Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
 * See the NOTICE file for attribution, trademark, and algorithm IP terms.
 *
 * The C interface between SMSD and a separately installed screening accelerator, such as the
 * smsd-cuda package. It holds plain C types only, so the accelerator may be built with another
 * compiler or toolkit than SMSD itself. SMSD forms every score from the integer overlaps the
 * accelerator returns, so an accelerated screen gives exactly the hits of the CPU screen.
 */
#ifndef SMSD_SCREEN_ACCELERATOR_ABI_H
#define SMSD_SCREEN_ACCELERATOR_ABI_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The interface version this header describes. */
#define SMSD_SCREEN_ACCELERATOR_ABI 1

/** The name of the Python capsule that carries a pointer to the accelerator's table. */
#define SMSD_SCREEN_ACCELERATOR_CAPSULE "smsd.screen_accelerator"

/** One entry of a target's label histogram: an atom label and the number of atoms carrying it. */
typedef struct smsd_label_count {
    uint16_t label;
    uint16_t count;
} smsd_label_count;

/** A selected target and the integer fingerprint intersection used to score it on the host. */
typedef struct smsd_fingerprint_hit {
    int32_t target_index;
    int32_t common_bits;
} smsd_fingerprint_hit;

/**
 * The functions an accelerator provides. Each returns 0 on success; on failure it returns a
 * non-zero value and writes a NUL-terminated reason of at most error_size bytes into error.
 * SMSD accepts a table whose abi_version is SMSD_SCREEN_ACCELERATOR_ABI, whose struct_size covers
 * the required prefix and whose smsd_version is the version of the SMSD build it is attached to.
 * An optional callback group is used only when struct_size covers every member of that group.
 */
typedef struct smsd_screen_accelerator_v1 {
    uint32_t abi_version;
    uint32_t struct_size;
    /** The SMSD version the accelerator was built for, such as "8.0.0". */
    const char* smsd_version;
    /** A short name for the accelerator, such as "CUDA". */
    const char* name;
    /** Writes a description of the device, such as "Tesla T4 (compute 7.5, 14930 MB)". */
    int (*device_info)(char* text, size_t text_size, char* error, size_t error_size);
    /**
     * For each of the num_targets targets, whose histogram is entries[offsets[i], offsets[i + 1]),
     * writes overlaps[i], the sum over its entries of the smaller of the entry's count and
     * query_counts[label]. Every label is below table_size, the length of query_counts.
     */
    int (*label_overlaps)(const int32_t* query_counts, int32_t table_size, int64_t num_targets,
                          const uint64_t* offsets, const smsd_label_count* entries, int32_t* overlaps,
                          char* error, size_t error_size);
    /** Optional tail: operation 1 is fingerprint screening, 2 is label screening, 3 is substructure screening.
     * UINT64_MAX means automatic use has no measured crossover for this device. */
    int (*minimum_targets)(uint32_t operation, uint64_t* minimum, char* error, size_t error_size);
    /** Uploads immutable target fingerprints in word-major order and their popcounts once. */
    int (*create_fingerprints)(int64_t num_targets, int32_t num_words, const uint64_t* word_major,
                               const int32_t* target_popcounts, void** resident, char* error, size_t error_size);
    /** Query words and cutoff rows are query-major. Output hits are grouped by query, then
     * strictly increasing target index; hit_offsets has num_queries + 1 entries. */
    int (*fingerprint_screen)(void* resident, int32_t num_queries, const uint64_t* query_words,
                              const int32_t* minimum_common, int32_t cutoff_stride,
                              smsd_fingerprint_hit* hits, uint64_t hit_capacity, uint64_t* hit_offsets,
                              char* error, size_t error_size);
    /** Uploads immutable target label offsets and entries once. */
    int (*create_labels)(int64_t num_targets, const uint64_t* offsets, const smsd_label_count* entries,
                         void** resident, char* error, size_t error_size);
    /** Query counts and output integer overlaps are query-major. */
    int (*resident_label_overlaps)(void* resident, int32_t num_queries, const int32_t* query_counts,
                                   int32_t table_size, int32_t* overlaps, char* error, size_t error_size);
    /** Consumes a resident handle, including when cleanup reports an error. The table and its
     * library must remain loaded until every handle it created has been destroyed. */
    int (*destroy_resident)(void* resident, char* error, size_t error_size);
} smsd_screen_accelerator_v1;

/** The required prefix; providers compiled before the optional tail remain compatible. */
#define SMSD_SCREEN_ACCELERATOR_BASE_SIZE offsetof(smsd_screen_accelerator_v1, minimum_targets)

#ifdef __cplusplus
}
#endif

#endif /* SMSD_SCREEN_ACCELERATOR_ABI_H */

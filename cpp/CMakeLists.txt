cmake_minimum_required(VERSION 3.18)
project(smsd VERSION 8.0.3 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Options ----
option(SMSD_BUILD_TESTS   "Build unit tests" ON)
option(SMSD_BUILD_PYTHON  "Build Python bindings via pybind11" OFF)
option(SMSD_BUILD_OPENMP  "Enable OpenMP for batch processing" ON)
option(SMSD_WITH_RDKIT    "Build RDKit adapter for SMILES/MOL parsing" OFF)
option(SMSD_HEADER_CHECK_WERROR "Treat warnings in the header self-containment check as errors" OFF)
# SMSD_BUILD_CUDA:  ON=always build, OFF=never build, AUTO=build if nvcc found (default)
set(SMSD_BUILD_CUDA  "AUTO" CACHE STRING "CUDA batch screening: ON / OFF / AUTO (default)")
# SMSD_BUILD_METAL: ON=always build, OFF=never build, AUTO=build on macOS if Metal found (default)
set(SMSD_BUILD_METAL "AUTO" CACHE STRING "Metal/MPS batch screening (macOS): ON / OFF / AUTO (default)")
foreach(accelerator IN ITEMS CUDA METAL)
    set_property(CACHE SMSD_BUILD_${accelerator} PROPERTY STRINGS AUTO ON OFF)
    string(TOUPPER "${SMSD_BUILD_${accelerator}}" accelerator_mode)
    if(NOT accelerator_mode STREQUAL "AUTO")
        if(SMSD_BUILD_${accelerator})
            set(accelerator_mode "ON")
        else()
            set(accelerator_mode "OFF")
        endif()
    endif()
    set(_SMSD_${accelerator}_MODE "${accelerator_mode}")
endforeach()

# ---- Header-only library (no RDKit dependency) ----
add_library(smsd INTERFACE)
find_package(Threads REQUIRED)
target_link_libraries(smsd INTERFACE Threads::Threads)
target_compile_features(smsd INTERFACE cxx_std_17)
target_compile_options(smsd INTERFACE $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/utf-8>)
target_include_directories(smsd INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

# ---- Optional: OpenMP for batch processing ----
set(SMSD_PACKAGE_REQUIRES_OPENMP FALSE)
if(SMSD_BUILD_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "SMSD: OpenMP found (version ${OpenMP_CXX_VERSION})")
        # Propagate OpenMP flags via the INTERFACE library so that any target
        # linking smsd automatically picks up OpenMP support for batch.hpp.
        target_link_libraries(smsd INTERFACE OpenMP::OpenMP_CXX)
        set(SMSD_PACKAGE_REQUIRES_OPENMP TRUE)
    else()
        message(STATUS "SMSD: OpenMP not found -- batch.hpp will use sequential fallback")
    endif()
endif()

# ---- Optional: RDKit adapter ----
set(SMSD_PACKAGE_REQUIRES_RDKIT FALSE)
if(SMSD_WITH_RDKIT)
    find_package(RDKit REQUIRED)
    add_library(smsd_rdkit STATIC src/rdkit_adapter.cpp)
    target_compile_definitions(smsd_rdkit PRIVATE SMSD_WITH_RDKIT)
    target_compile_features(smsd_rdkit PUBLIC cxx_std_20)
    set_target_properties(smsd_rdkit PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(smsd_rdkit PUBLIC smsd RDKit::SmilesParse RDKit::GraphMol
        RDKit::RingDecomposerLib RDKit::Depictor)
    set(SMSD_PACKAGE_REQUIRES_RDKIT TRUE)
endif()

# ---- Optional: CUDA batch screening ----
# AUTO: probe for nvcc / CUDAToolkit and enable if found
# ON:   require CUDA (fails configuration if not found)
# OFF:  skip entirely
if(_SMSD_CUDA_MODE STREQUAL "AUTO")
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
        include(CheckLanguage)
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
            set(_SMSD_CUDA_ENABLED TRUE)
            message(STATUS "SMSD: CUDA ${CUDAToolkit_VERSION} found — GPU batch screening enabled (auto)")
        else()
            set(_SMSD_CUDA_ENABLED FALSE)
            message(STATUS "SMSD: CUDA compiler unavailable — GPU batch screening disabled")
        endif()
    else()
        set(_SMSD_CUDA_ENABLED FALSE)
        message(STATUS "SMSD: CUDA not found — GPU batch screening disabled (CPU/OpenMP only)")
    endif()
elseif(_SMSD_CUDA_MODE STREQUAL "ON")
    find_package(CUDAToolkit REQUIRED)
    set(_SMSD_CUDA_ENABLED TRUE)
    message(STATUS "SMSD: CUDA ${CUDAToolkit_VERSION} found — GPU batch screening enabled (forced)")
else()
    set(_SMSD_CUDA_ENABLED FALSE)
    message(STATUS "SMSD: CUDA disabled by SMSD_BUILD_CUDA=OFF")
endif()

if(_SMSD_CUDA_ENABLED)
    # Keep compilation and final device links on the same architecture set.
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        if(NOT "$ENV{CUDAARCHS}" STREQUAL "")
            set(CMAKE_CUDA_ARCHITECTURES "$ENV{CUDAARCHS}" CACHE STRING "CUDA architectures")
        else()
            set(CMAKE_CUDA_ARCHITECTURES "75;80;86;89;90" CACHE STRING "CUDA architectures")
        endif()
    endif()
    enable_language(CUDA)
    add_library(smsd_cuda STATIC cuda/batch_screen.cu cuda/graph_kernels.cu)
    target_link_libraries(smsd_cuda PUBLIC smsd CUDA::cudart)
    set_target_properties(smsd_cuda PROPERTIES
        CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
        CUDA_SEPARABLE_COMPILATION ON)
    # Expose SMSD_ENABLE_CUDA so gpu.hpp picks up the right branch
    target_compile_definitions(smsd_cuda PUBLIC SMSD_ENABLE_CUDA)
endif()

# ---- Optional: Metal/MPS batch screening (macOS only) ----
# AUTO: enable on macOS when Metal.framework is present (always true on macOS 10.13+)
# ON:   require Metal — fails configuration on non-macOS
# OFF:  skip entirely
set(_SMSD_METAL_ENABLED FALSE)
if(_SMSD_METAL_MODE STREQUAL "ON" AND NOT APPLE)
    message(FATAL_ERROR "SMSD: SMSD_BUILD_METAL=ON but Metal is only available on macOS/Apple platforms.")
endif()

if(APPLE AND NOT _SMSD_METAL_MODE STREQUAL "OFF")
    find_library(METAL_LIBRARY    Metal)
    find_library(FOUNDATION_LIBRARY Foundation)
    if(METAL_LIBRARY AND FOUNDATION_LIBRARY)
        set(_SMSD_METAL_ENABLED TRUE)
        message(STATUS "SMSD: Metal.framework found — GPU batch screening enabled on macOS")
    else()
        if(_SMSD_METAL_MODE STREQUAL "ON")
            message(FATAL_ERROR "SMSD: SMSD_BUILD_METAL=ON but Metal.framework was not found.")
        else()
            message(STATUS "SMSD: Metal.framework not found — GPU batch screening disabled (CPU/OpenMP only)")
        endif()
    endif()
elseif(NOT APPLE)
    message(STATUS "SMSD: Metal not applicable (non-Apple platform)")
else()
    message(STATUS "SMSD: Metal disabled by SMSD_BUILD_METAL=OFF")
endif()

if(_SMSD_METAL_ENABLED)
    # Objective-C++ requires the OBJCXX language to be enabled.
    enable_language(OBJCXX)
    add_library(smsd_metal STATIC metal/metal_batch.mm metal/metal_graph_kernels.mm metal/metal_fingerprint.mm)
    # ARC (Automatic Reference Counting) for the Objective-C objects.
    target_compile_options(smsd_metal PRIVATE -fobjc-arc -fobjc-arc-exceptions)
    target_link_libraries(smsd_metal PUBLIC smsd ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
    # Expose SMSD_ENABLE_METAL so gpu.hpp picks up the Metal branch.
    target_compile_definitions(smsd_metal PUBLIC SMSD_ENABLE_METAL)
    target_compile_definitions(smsd_metal PRIVATE SMSD_METAL_SMSD_VERSION="${PROJECT_VERSION}")
endif()

# ---- Optional: Python bindings ----
if(SMSD_BUILD_PYTHON)
    # Discover one interpreter and its extension-module headers together.
    # Development.Module avoids requiring an embedding libpython on Unix.
    set(PYBIND11_FINDPYTHON ON)
    find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(_smsd bindings/pybind11/smsd_bindings.cpp)
    target_link_libraries(_smsd PRIVATE smsd)
    # The extension reports the project version, which the package checks on import.
    target_compile_definitions(_smsd PRIVATE SMSD_VERSION_INFO="${PROJECT_VERSION}")
    if(SMSD_WITH_RDKIT)
        target_link_libraries(_smsd PRIVATE smsd_rdkit)
    endif()
    if(_SMSD_CUDA_ENABLED)
        target_link_libraries(_smsd PRIVATE smsd_cuda)
    endif()
    if(_SMSD_METAL_ENABLED)
        target_link_libraries(_smsd PRIVATE smsd_metal)
    endif()
    # Install the Python extension module (required by scikit-build-core)
    install(TARGETS _smsd LIBRARY DESTINATION smsd RUNTIME DESTINATION smsd)
endif()

# ---- Tests ----
if(SMSD_BUILD_TESTS)
    enable_testing()

    # Builds one test program against the header-only library, keeps assert()
    # active in optimised builds and registers it with CTest under its own name.
    function(smsd_add_test name)
        add_executable(${name} ${ARGN})
        target_link_libraries(${name} PRIVATE smsd)
        # Large test translation units exceed MSVC's default object-file section count.
        target_compile_options(${name} PRIVATE
            $<$<CXX_COMPILER_ID:MSVC>:/UNDEBUG /bigobj>
            $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-UNDEBUG>)
        add_test(NAME ${name} COMMAND ${name})
    endfunction()

    # Links the optional RDKit, CUDA and Metal back ends that this build enables.
    function(smsd_link_backends target)
        if(SMSD_WITH_RDKIT)
            target_link_libraries(${target} PRIVATE smsd_rdkit)
        endif()
        if(_SMSD_CUDA_ENABLED)
            target_link_libraries(${target} PRIVATE smsd_cuda)
        endif()
        if(_SMSD_METAL_ENABLED)
            target_link_libraries(${target} PRIVATE smsd_metal)
        endif()
    endfunction()

    # Every public header must compile on its own, so a program can include
    # just the part of the library it uses. One translation unit is generated
    # per header under include/smsd (detail/ excluded), holding nothing but
    # that #include, and compiled as strict C++17 with the common warnings
    # enabled, as errors when SMSD_HEADER_CHECK_WERROR is on. The objects are
    # built by the CTest entry, not by the default build, and a header added
    # later is picked up at the next configure.
    file(GLOB _smsd_public_headers CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/include/smsd/*.hpp")
    set(_smsd_header_units "")
    foreach(header IN LISTS _smsd_public_headers)
        get_filename_component(header_name "${header}" NAME)
        get_filename_component(header_stem "${header}" NAME_WE)
        set(header_unit "${CMAKE_CURRENT_BINARY_DIR}/header-check/${header_stem}.cpp")
        file(CONFIGURE OUTPUT "${header_unit}"
            CONTENT "#include \"smsd/${header_name}\"\n" @ONLY)
        list(APPEND _smsd_header_units "${header_unit}")
    endforeach()
    add_library(smsd_header_check OBJECT EXCLUDE_FROM_ALL ${_smsd_header_units})
    target_link_libraries(smsd_header_check PRIVATE smsd)
    set_target_properties(smsd_header_check PROPERTIES CXX_EXTENSIONS OFF)
    if(MSVC)
        target_compile_options(smsd_header_check PRIVATE /W4 /permissive-
            $<$<BOOL:${SMSD_HEADER_CHECK_WERROR}>:/WX>)
    else()
        target_compile_options(smsd_header_check PRIVATE -Wall -Wextra
            $<$<BOOL:${SMSD_HEADER_CHECK_WERROR}>:-Werror>)
    endif()
    add_test(NAME smsd_header_self_containment_tests
        COMMAND "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}"
            --target smsd_header_check --config $<CONFIG>)

    # Search API tests
    smsd_add_test(smsd_search_api_tests tests/test_search_api.cpp)
    smsd_add_test(smsd_smarts_basics_tests tests/test_smarts_basics.cpp)
    smsd_add_test(smsd_feature_backend_regression_tests tests/test_feature_backend_regression.cpp)
    smsd_add_test(smsd_device_workspace_regression_tests tests/test_device_workspace_regression.cpp)
    smsd_add_test(smsd_tautomer_confidence_tests tests/test_tautomer_confidence.cpp)
    # The atom-order test reads a benchmark corpus, which the Python source
    # distribution leaves out; without it the test reports a skip.
    set(_smsd_permutation_corpus "${CMAKE_CURRENT_SOURCE_DIR}/../benchmarks/data/chembl_mcs_benchmark.smi")
    if(EXISTS "${_smsd_permutation_corpus}")
        target_compile_definitions(smsd_tautomer_confidence_tests PRIVATE
            SMSD_TAUTOMER_PERMUTATION_SMILES="${_smsd_permutation_corpus}")
    endif()
    smsd_add_test(smsd_fingerprint_extent_regression_tests tests/test_fingerprint_extent_regression.cpp)
    smsd_add_test(smsd_fingerprint_options_regression_tests tests/test_fingerprint_options_regression.cpp)
    smsd_add_test(smsd_search_bound_regression_tests tests/test_search_bound_regression.cpp)
    smsd_add_test(smsd_mces_contract_regression_tests tests/test_mces_contract_regression.cpp)
    # Chemistry rules shared with the Java and Python tests.
    smsd_add_test(smsd_chemistry_rules_golden_tests tests/test_chemistry_rules_golden.cpp)
    target_compile_definitions(smsd_chemistry_rules_golden_tests PRIVATE
        SMSD_CHEMISTRY_RULES_JSON="${CMAKE_CURRENT_SOURCE_DIR}/../java/src/test/resources/com/bioinception/smsd/chemistry_rules.json")
    add_test(NAME smsd_cmake_option_regression_tests COMMAND "${CMAKE_COMMAND}"
        "-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}" "-DBINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}/option-tests"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cmake_options.cmake")
    if(APPLE)
        enable_language(OBJCXX)
        find_library(SMSD_TEST_FOUNDATION Foundation REQUIRED)
        find_library(SMSD_TEST_METAL Metal REQUIRED)
        add_executable(smsd_metal_failure_regression_tests tests/test_metal_failure_regression.mm)
        target_link_libraries(smsd_metal_failure_regression_tests PRIVATE smsd ${SMSD_TEST_FOUNDATION} ${SMSD_TEST_METAL})
        target_compile_options(smsd_metal_failure_regression_tests PRIVATE -fobjc-arc)
        foreach(failure IN ITEMS pipeline buffer command encoder status disabled)
            add_test(NAME smsd_metal_failure_${failure}_tests COMMAND smsd_metal_failure_regression_tests ${failure})
        endforeach()
        add_executable(smsd_metal_fingerprint_regression_tests tests/test_metal_fingerprint_regression.mm)
        target_link_libraries(smsd_metal_fingerprint_regression_tests PRIVATE smsd ${SMSD_TEST_FOUNDATION} ${SMSD_TEST_METAL})
        target_compile_options(smsd_metal_fingerprint_regression_tests PRIVATE -fobjc-arc -fobjc-arc-exceptions)
        target_compile_definitions(smsd_metal_fingerprint_regression_tests PRIVATE SMSD_METAL_SMSD_VERSION="${PROJECT_VERSION}")
        foreach(mode IN ITEMS contract concurrent device queue pipeline function buffer contents exception command encoder status disabled)
            add_test(NAME smsd_metal_fingerprint_${mode}_tests COMMAND smsd_metal_fingerprint_regression_tests ${mode})
        endforeach()
    endif()

    # Core test suite (substructure, MCS, ring perception and chemistry)
    smsd_add_test(smsd_core_tests tests/test_core.cpp)
    smsd_link_backends(smsd_core_tests)

    if(SMSD_WITH_RDKIT)
        smsd_add_test(smsd_rdkit_adapter_regression_tests tests/test_rdkit_adapter_regression.cpp)
        target_link_libraries(smsd_rdkit_adapter_regression_tests PRIVATE smsd_rdkit)
    endif()

    # Parser test suite (SMARTS, hydrogen handling and the comprehensive SMILES
    # cases, which are compiled from test_smiles_comprehensive.cpp as a second
    # translation unit)
    smsd_add_test(smsd_parser_tests tests/test_parsers.cpp tests/test_smiles_comprehensive.cpp)
    smsd_add_test(smsd_parser_property_regression_tests tests/test_parser_property_regression.cpp)
    smsd_add_test(smsd_parser_robustness_regression_tests tests/test_parser_robustness_regression.cpp)
    smsd_add_test(smsd_member_smiles_regression_tests tests/test_member_smiles_regression.cpp)
    add_test(NAME smsd_member_smiles_parity_tests COMMAND smsd_member_smiles_regression_tests parity)
    smsd_add_test(smsd_thread_safety_regression_tests tests/test_thread_safety_regression.cpp)
    smsd_add_test(smsd_corpus_storage_regression_tests tests/test_corpus_storage_regression.cpp)
    smsd_add_test(smsd_graph_identity_regression_tests tests/test_graph_identity_regression.cpp)
    smsd_add_test(smsd_aromatic_ring_regression_tests tests/test_aromatic_ring_regression.cpp)
    smsd_add_test(smsd_aromatic_stereo_writer_tests tests/test_aromatic_stereo_writer.cpp)
    smsd_add_test(smsd_stereo_io_regression_tests tests/test_stereo_io_regression.cpp)
    smsd_add_test(smsd_stereo_matching_regression_tests tests/test_stereo_matching_regression.cpp)
    smsd_add_test(smsd_input_validation_regression_tests tests/test_input_validation_regression.cpp)
    smsd_add_test(smsd_mol_properties_regression_tests tests/test_mol_properties_regression.cpp)
    smsd_add_test(smsd_depiction_regression_tests tests/test_depiction_regression.cpp)

    # Batch and GPU test suite (OpenMP, fingerprints, CUDA and Metal)
    smsd_add_test(smsd_batch_gpu_tests tests/test_batch_gpu.cpp)
    smsd_link_backends(smsd_batch_gpu_tests)

    smsd_add_test(smsd_cip_tests tests/test_cip.cpp)
    smsd_add_test(smsd_layout_tests tests/test_layout.cpp)
    smsd_add_test(smsd_portability_tests tests/test_portability.cpp)
    # The same checks on the portable bit operations that targets without an intrinsic use.
    smsd_add_test(smsd_portable_bitops_tests tests/test_portability.cpp)
    target_compile_definitions(smsd_portable_bitops_tests PRIVATE SMSD_PORTABLE_BITOPS)
    smsd_add_test(smsd_general_matching_tests tests/test_general_matching.cpp)
    smsd_add_test(smsd_substructure_regression_tests tests/test_substructure_regression.cpp)
    smsd_add_test(smsd_gpu_domain_regression_tests tests/test_gpu_domain_regression.cpp)
    smsd_add_test(smsd_search_enumeration_tests tests/test_search_enumeration.cpp)
    smsd_add_test(smsd_assignment_regression_tests tests/test_assignment_regression.cpp)
    smsd_add_test(smsd_mcs_regression_tests tests/test_mcs_regression.cpp)
    smsd_add_test(smsd_vertex_search_stereo_tests tests/test_vertex_search_stereo.cpp)
    smsd_add_test(smsd_mcs_ring_proof_status_tests tests/test_mcs_ring_proof_status.cpp)
    smsd_add_test(smsd_mcs_complete_ring_domain_tests tests/test_mcs_complete_ring_domains.cpp)

    if(_SMSD_CUDA_ENABLED)
        smsd_add_test(smsd_cuda_host_link_tests tests/test_cuda_host_link.cpp)
        target_link_libraries(smsd_cuda_host_link_tests PRIVATE smsd_cuda)
    endif()

    if(_SMSD_METAL_ENABLED OR _SMSD_CUDA_ENABLED)
        smsd_add_test(smsd_native_gpu_domain_tests tests/test_gpu_domain_regression.cpp)
        target_compile_definitions(smsd_native_gpu_domain_tests PRIVATE SMSD_TEST_NATIVE_GPU)
        if(_SMSD_METAL_ENABLED)
            target_link_libraries(smsd_native_gpu_domain_tests PRIVATE smsd_metal)
        else()
            target_link_libraries(smsd_native_gpu_domain_tests PRIVATE smsd_cuda)
        endif()
        set_tests_properties(smsd_native_gpu_domain_tests PROPERTIES SKIP_RETURN_CODE 77)
    endif()

    smsd_add_test(smsd_screening_regression_tests tests/test_screening_regression.cpp)
    smsd_add_test(smsd_screen_accelerator_tests tests/test_screen_accelerator.cpp)
    smsd_add_test(smsd_legacy_screening_threshold_tests tests/test_legacy_screening_thresholds.cpp)
    if(_SMSD_CUDA_ENABLED)
        smsd_add_test(smsd_cuda_screening_threshold_tests tests/test_legacy_screening_thresholds.cpp)
        target_compile_definitions(smsd_cuda_screening_threshold_tests PRIVATE SMSD_TEST_NATIVE_GPU)
        target_link_libraries(smsd_cuda_screening_threshold_tests PRIVATE smsd_cuda)
        set_tests_properties(smsd_cuda_screening_threshold_tests PROPERTIES SKIP_RETURN_CODE 77)
    endif()
    smsd_add_test(smsd_corpus_substructure_screen_tests tests/test_corpus_substructure_screen.cpp)
    smsd_add_test(smsd_fragment_output_regression_tests tests/test_fragment_output_regression.cpp)
    smsd_add_test(smsd_option_fingerprint_regression_tests tests/test_option_fingerprint_regression.cpp)
    if(_SMSD_METAL_ENABLED)
        smsd_add_test(smsd_metal_resident_device_tests tests/test_metal_resident_device.cpp)
        target_link_libraries(smsd_metal_resident_device_tests PRIVATE smsd_metal)
        set_tests_properties(smsd_metal_resident_device_tests PROPERTIES SKIP_RETURN_CODE 77)
        add_executable(smsd_metal_threshold_regression_tests tests/test_metal_threshold_regression.mm)
        target_link_libraries(smsd_metal_threshold_regression_tests PRIVATE smsd ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
        target_compile_options(smsd_metal_threshold_regression_tests PRIVATE -fobjc-arc)
        add_test(NAME smsd_metal_threshold_regression_tests COMMAND smsd_metal_threshold_regression_tests)
        set_tests_properties(smsd_metal_threshold_regression_tests PROPERTIES SKIP_RETURN_CODE 77)
        smsd_add_test(smsd_metal_screening_regression_tests tests/test_screening_regression.cpp)
        target_link_libraries(smsd_metal_screening_regression_tests PRIVATE smsd_metal)
    endif()
endif()

# ---- Install ----
set(SMSD_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/smsd")
configure_package_config_file(
    cmake/smsdConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/smsdConfig.cmake"
    INSTALL_DESTINATION "${SMSD_INSTALL_CMAKEDIR}")
set(SMSD_PACKAGE_ARCH_OPTION "")
if(NOT SMSD_WITH_RDKIT)
    set(SMSD_PACKAGE_ARCH_OPTION ARCH_INDEPENDENT)
endif()
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/smsdConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
    ${SMSD_PACKAGE_ARCH_OPTION})

install(DIRECTORY include/smsd DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS smsd EXPORT smsdTargets)
if(SMSD_WITH_RDKIT)
    install(TARGETS smsd_rdkit EXPORT smsdTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
endif()
install(EXPORT smsdTargets NAMESPACE smsd:: DESTINATION "${SMSD_INSTALL_CMAKEDIR}")
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/smsdConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/smsdConfigVersion.cmake"
    DESTINATION "${SMSD_INSTALL_CMAKEDIR}")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE"
    "${CMAKE_CURRENT_SOURCE_DIR}/../NOTICE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/smsd")

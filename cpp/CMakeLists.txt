cmake_minimum_required(VERSION 3.18)
project(smsd VERSION 6.10.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Options ----
option(SMSD_BUILD_TESTS   "Build unit tests" ON)
option(SMSD_BUILD_PYTHON  "Build Python bindings via pybind11" OFF)
option(SMSD_BUILD_OPENMP  "Enable OpenMP for batch processing" ON)
option(SMSD_WITH_RDKIT    "Build RDKit adapter for SMILES/MOL parsing" OFF)
# SMSD_BUILD_CUDA:  ON=always build, OFF=never build, AUTO=build if nvcc found (default)
set(SMSD_BUILD_CUDA  "AUTO" CACHE STRING "CUDA batch screening: ON / OFF / AUTO (default)")
# SMSD_BUILD_METAL: ON=always build, OFF=never build, AUTO=build on macOS if Metal found (default)
set(SMSD_BUILD_METAL "AUTO" CACHE STRING "Metal/MPS batch screening (macOS): ON / OFF / AUTO (default)")

# ---- Header-only library (no RDKit dependency) ----
add_library(smsd INTERFACE)
target_include_directories(smsd INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# ---- Optional: OpenMP for batch processing ----
if(SMSD_BUILD_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "SMSD: OpenMP found (version ${OpenMP_CXX_VERSION})")
        # Propagate OpenMP flags via the INTERFACE library so that any target
        # linking smsd automatically picks up OpenMP support for batch.hpp.
        target_link_libraries(smsd INTERFACE OpenMP::OpenMP_CXX)
    else()
        message(STATUS "SMSD: OpenMP not found -- batch.hpp will use sequential fallback")
    endif()
endif()

# ---- Optional: RDKit adapter ----
if(SMSD_WITH_RDKIT)
    find_package(RDKit REQUIRED)
    add_library(smsd_rdkit STATIC src/rdkit_adapter.cpp)
    target_link_libraries(smsd_rdkit PUBLIC smsd RDKit::SmilesParse RDKit::GraphMol
        RDKit::RingDecomposerLib RDKit::Depictor)
    target_include_directories(smsd_rdkit PUBLIC ${RDKit_INCLUDE_DIRS})
endif()

# ---- Optional: CUDA batch screening ----
# AUTO: probe for nvcc / CUDAToolkit and enable if found
# ON:   require CUDA (fails configuration if not found)
# OFF:  skip entirely
if(SMSD_BUILD_CUDA STREQUAL "AUTO")
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
        set(_SMSD_CUDA_ENABLED TRUE)
        message(STATUS "SMSD: CUDA ${CUDAToolkit_VERSION} found — GPU batch screening enabled (auto)")
    else()
        set(_SMSD_CUDA_ENABLED FALSE)
        message(STATUS "SMSD: CUDA not found — GPU batch screening disabled (CPU/OpenMP only)")
    endif()
elseif(SMSD_BUILD_CUDA)
    find_package(CUDAToolkit REQUIRED)
    set(_SMSD_CUDA_ENABLED TRUE)
    message(STATUS "SMSD: CUDA ${CUDAToolkit_VERSION} found — GPU batch screening enabled (forced)")
else()
    set(_SMSD_CUDA_ENABLED FALSE)
    message(STATUS "SMSD: CUDA disabled by SMSD_BUILD_CUDA=OFF")
endif()

if(_SMSD_CUDA_ENABLED)
    enable_language(CUDA)
    add_library(smsd_cuda STATIC cuda/batch_screen.cu cuda/graph_kernels.cu)
    target_link_libraries(smsd_cuda PUBLIC smsd CUDA::cudart)
    set_target_properties(smsd_cuda PROPERTIES
        CUDA_ARCHITECTURES "70;80;86;89;90"   # Volta, Ampere, Ada, Hopper
        CUDA_SEPARABLE_COMPILATION ON)
    target_compile_features(smsd_cuda PUBLIC cxx_std_17)
    # Expose SMSD_ENABLE_CUDA so gpu.hpp picks up the right branch
    target_compile_definitions(smsd_cuda PUBLIC SMSD_ENABLE_CUDA)
    target_include_directories(smsd_cuda PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
endif()

# ---- Optional: Metal/MPS batch screening (macOS only) ----
# AUTO: enable on macOS when Metal.framework is present (always true on macOS 10.13+)
# ON:   require Metal — fails configuration on non-macOS
# OFF:  skip entirely
set(_SMSD_METAL_ENABLED FALSE)
if(SMSD_BUILD_METAL STREQUAL "ON" AND NOT APPLE)
    message(FATAL_ERROR "SMSD: SMSD_BUILD_METAL=ON but Metal is only available on macOS/Apple platforms.")
endif()

if(APPLE AND NOT (SMSD_BUILD_METAL STREQUAL "OFF"))
    find_library(METAL_LIBRARY    Metal)
    find_library(FOUNDATION_LIBRARY Foundation)
    if(METAL_LIBRARY AND FOUNDATION_LIBRARY)
        set(_SMSD_METAL_ENABLED TRUE)
        message(STATUS "SMSD: Metal.framework found — GPU batch screening enabled on macOS")
    else()
        if(SMSD_BUILD_METAL STREQUAL "ON")
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
    add_library(smsd_metal STATIC metal/metal_batch.mm metal/metal_graph_kernels.mm)
    # ARC (Automatic Reference Counting) for the Objective-C objects.
    target_compile_options(smsd_metal PRIVATE -fobjc-arc)
    target_link_libraries(smsd_metal PUBLIC smsd ${METAL_LIBRARY} ${FOUNDATION_LIBRARY})
    # Expose SMSD_ENABLE_METAL so gpu.hpp picks up the Metal branch.
    target_compile_definitions(smsd_metal PUBLIC SMSD_ENABLE_METAL)
    target_include_directories(smsd_metal PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    # Mark .mm as Objective-C++ (CMake may infer this from the extension,
    # but being explicit avoids edge cases with older CMake versions).
    set_source_files_properties(metal/metal_batch.mm PROPERTIES
        LANGUAGE OBJCXX)
endif()

# ---- Optional: Python bindings ----
if(SMSD_BUILD_PYTHON)
    find_package(pybind11 REQUIRED)
    pybind11_add_module(_smsd bindings/pybind11/smsd_bindings.cpp)
    target_link_libraries(_smsd PRIVATE smsd)
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

    # Core test suite (substructure, MCS, ring perception, chemistry, v6.4.1 features)
    add_executable(smsd_core_tests tests/test_core.cpp)
    target_link_libraries(smsd_core_tests PRIVATE smsd)
    if(SMSD_WITH_RDKIT)
        target_link_libraries(smsd_core_tests PRIVATE smsd_rdkit)
    endif()
    if(_SMSD_CUDA_ENABLED)
        target_link_libraries(smsd_core_tests PRIVATE smsd_cuda)
    endif()
    if(_SMSD_METAL_ENABLED)
        target_link_libraries(smsd_core_tests PRIVATE smsd_metal)
    endif()
    add_test(NAME smsd_core_tests COMMAND smsd_core_tests)

    # Parser test suite (SMARTS, hydrogen handling, comprehensive SMILES)
    # The comprehensive SMILES tests (441) are compiled from test_smiles_comprehensive.cpp
    # as a second translation unit, gated by SMSD_TEST_SMILES_COMPREHENSIVE.
    add_executable(smsd_parser_tests
        tests/test_parsers.cpp
        tests/test_smiles_comprehensive.cpp)
    target_compile_definitions(smsd_parser_tests PRIVATE
        SMSD_TEST_SMILES_COMPREHENSIVE SMSD_PARSERS_SUITE)
    target_link_libraries(smsd_parser_tests PRIVATE smsd)
    add_test(NAME smsd_parser_tests COMMAND smsd_parser_tests)

    # Batch & GPU test suite (OpenMP, fingerprint, CUDA/Metal)
    add_executable(smsd_batch_gpu_tests tests/test_batch_gpu.cpp)
    target_link_libraries(smsd_batch_gpu_tests PRIVATE smsd)
    if(SMSD_WITH_RDKIT)
        target_link_libraries(smsd_batch_gpu_tests PRIVATE smsd_rdkit)
    endif()
    if(_SMSD_CUDA_ENABLED)
        target_link_libraries(smsd_batch_gpu_tests PRIVATE smsd_cuda)
    endif()
    if(_SMSD_METAL_ENABLED)
        target_link_libraries(smsd_batch_gpu_tests PRIVATE smsd_metal)
    endif()
    add_test(NAME smsd_batch_gpu_tests COMMAND smsd_batch_gpu_tests)

    # CIP stereo descriptor test suite
    add_executable(smsd_cip_tests tests/test_cip.cpp)
    target_link_libraries(smsd_cip_tests PRIVATE smsd)
    add_test(NAME smsd_cip_tests COMMAND smsd_cip_tests)

    # Layout test suite
    add_executable(smsd_layout_tests tests/test_layout.cpp)
    target_link_libraries(smsd_layout_tests PRIVATE smsd)
    add_test(NAME smsd_layout_tests COMMAND smsd_layout_tests)
endif()

# ---- Install ----
install(DIRECTORY include/smsd DESTINATION include)
install(TARGETS smsd EXPORT smsdTargets)
install(EXPORT smsdTargets NAMESPACE smsd:: DESTINATION lib/cmake/smsd)

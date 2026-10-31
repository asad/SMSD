# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# The smsd-cuda package: one shared library holding the CUDA screening accelerator, which the
# smsd package loads at run time. It is built only where the CUDA toolkit is installed and never
# by the smsd build itself.
cmake_minimum_required(VERSION 3.18)
project(smsd_cuda LANGUAGES CXX)
include(GNUInstallDirs)

# The SMSD sources the package is built from. Their version is the smsd version it serves, and
# the package version must equal it.
set(SMSD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../cpp" CACHE PATH "The SMSD C++ sources (the cpp folder)")
if(NOT EXISTS "${SMSD_SOURCE_DIR}/include/smsd/screen_accelerator_abi.h")
    message(FATAL_ERROR "smsd-cuda: the SMSD sources were not found at ${SMSD_SOURCE_DIR}. "
        "Build it from an SMSD checkout, or set SMSD_SOURCE_DIR to its cpp folder.")
endif()
file(STRINGS "${SMSD_SOURCE_DIR}/CMakeLists.txt" _smsd_project REGEX "^project\\(smsd VERSION [0-9]+\\.[0-9]+\\.[0-9]+")
string(REGEX REPLACE "^project\\(smsd VERSION ([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" SMSD_VERSION "${_smsd_project}")
if(NOT SMSD_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "smsd-cuda: the smsd version could not be read from ${SMSD_SOURCE_DIR}/CMakeLists.txt.")
endif()
if(DEFINED SKBUILD_PROJECT_VERSION AND NOT SKBUILD_PROJECT_VERSION STREQUAL SMSD_VERSION)
    message(FATAL_ERROR "smsd-cuda ${SKBUILD_PROJECT_VERSION} must be built from the smsd "
        "${SKBUILD_PROJECT_VERSION} sources, and these are smsd ${SMSD_VERSION}.")
endif()

include(CheckLanguage)
check_language(CUDA)
if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "smsd-cuda needs the NVIDIA CUDA toolkit 12, and no CUDA compiler was found. "
        "Install the toolkit and put nvcc on PATH, or set CUDACXX or CMAKE_CUDA_COMPILER to it. "
        "The smsd package itself needs no CUDA; this package is built only where CUDA is installed.")
endif()
enable_language(CUDA)
find_package(CUDAToolkit 12.0 REQUIRED)
if(NOT CUDAToolkit_VERSION_MAJOR EQUAL 12)
    message(FATAL_ERROR "smsd-cuda depends on nvidia-cuda-runtime-cu12 and needs the CUDA 12 toolkit; "
        "found CUDA ${CUDAToolkit_VERSION}.")
endif()
# The package requires a runtime at least this new, so it is built with a toolkit no newer.
set(SMSD_CUDA_RUNTIME_MINIMUM "12.8" CACHE STRING "The oldest nvidia-cuda-runtime-cu12 the package accepts")
if("${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR}" VERSION_GREATER "${SMSD_CUDA_RUNTIME_MINIMUM}")
    message(WARNING "smsd-cuda is built with CUDA ${CUDAToolkit_VERSION} but accepts runtimes from "
        "${SMSD_CUDA_RUNTIME_MINIMUM}; build with CUDA ${SMSD_CUDA_RUNTIME_MINIMUM} or raise the bound in pyproject.toml.")
endif()

# CUDA 12.8 adds Blackwell device code for data-centre (10.0) and RTX (12.0) GPUs. Each has its
# own code and PTX, so later GPUs of either family can compile the kernels when they load.
if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "12.8")
    set(_smsd_cuda_architectures "75-real;80-real;86-real;89-real;90-real;100;120")
else()
    set(_smsd_cuda_architectures "75-real;80-real;86-real;89-real;90")
endif()
set(SMSD_CUDA_ARCHITECTURES "${_smsd_cuda_architectures}" CACHE STRING "The CUDA architectures to build for")
option(SMSD_CUDA_BINARY_MMA_PROTOTYPE "Build the opt-in SM80-SM90 binary MMA prototype" OFF)
option(SMSD_CUDA_BUILD_LAYOUT_TESTS "Build the CPU binary MMA layout oracle" OFF)
message(STATUS "smsd-cuda ${SMSD_VERSION}: CUDA ${CUDAToolkit_VERSION}, architectures ${SMSD_CUDA_ARCHITECTURES}")

add_library(smsd_cuda_screen SHARED src/screen_accelerator.cu)
target_include_directories(smsd_cuda_screen PRIVATE "${SMSD_SOURCE_DIR}/include")
target_compile_definitions(smsd_cuda_screen PRIVATE SMSD_CUDA_SMSD_VERSION="${SMSD_VERSION}")
target_compile_features(smsd_cuda_screen PRIVATE cuda_std_17)
if(SMSD_CUDA_BINARY_MMA_PROTOTYPE)
    target_compile_definitions(smsd_cuda_screen PRIVATE SMSD_CUDA_BINARY_MMA_PROTOTYPE=1)
endif()
if(SMSD_CUDA_BUILD_LAYOUT_TESTS)
    enable_testing()
    add_executable(smsd_cuda_binary_mma_layout_tests tests/test_binary_mma_layout.cpp)
    target_compile_features(smsd_cuda_binary_mma_layout_tests PRIVATE cxx_std_17)
    add_test(NAME smsd_cuda_binary_mma_layout_tests COMMAND smsd_cuda_binary_mma_layout_tests)
endif()
target_link_libraries(smsd_cuda_screen PRIVATE CUDA::cudart)
set_target_properties(smsd_cuda_screen PROPERTIES
    CUDA_ARCHITECTURES "${SMSD_CUDA_ARCHITECTURES}"
    CUDA_RUNTIME_LIBRARY Shared
    CXX_VISIBILITY_PRESET hidden
    CUDA_VISIBILITY_PRESET hidden
    # The CUDA runtime from nvidia-cuda-runtime-cu12, installed beside smsd_cuda in site-packages.
    INSTALL_RPATH "$ORIGIN/../nvidia/cuda_runtime/lib")
install(TARGETS smsd_cuda_screen LIBRARY DESTINATION smsd_cuda RUNTIME DESTINATION smsd_cuda)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
    "${CMAKE_CURRENT_SOURCE_DIR}/NOTICE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/smsd-cuda")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/licenses/cub/LICENSE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/smsd-cuda/licenses/cub")
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/licenses/cccl/LICENSE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/smsd-cuda/licenses/cccl")

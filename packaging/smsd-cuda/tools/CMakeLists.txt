# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# screen_compare: the C++ library built with SMSD_BUILD_CUDA=ON, comparing its CUDA screen with
# the CPU screen. Built only where the CUDA toolkit is installed.
cmake_minimum_required(VERSION 3.18)
project(smsd_cuda_tools LANGUAGES CXX)

include(CheckLanguage)
check_language(CUDA)
if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "screen_compare needs the NVIDIA CUDA toolkit, and no CUDA compiler was found. "
        "Install the toolkit and put nvcc on PATH, or set CUDACXX or CMAKE_CUDA_COMPILER to it.")
endif()
# The CUDA objects of the library are linked into the program here, so CUDA is enabled at this level.
enable_language(CUDA)

set(SMSD_BUILD_CUDA ON CACHE STRING "CUDA batch screening" FORCE)
set(SMSD_BUILD_METAL OFF CACHE STRING "Metal batch screening" FORCE)
set(SMSD_BUILD_TESTS OFF CACHE BOOL "Build unit tests" FORCE)
set(SMSD_BUILD_PYTHON OFF CACHE BOOL "Build Python bindings" FORCE)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../../../cpp" smsd)

add_executable(screen_compare screen_compare.cpp)
target_link_libraries(screen_compare PRIVATE smsd smsd_cuda)

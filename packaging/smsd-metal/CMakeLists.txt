# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
#
# The smsd-metal package: one shared library holding the Apple Metal screening accelerator, which
# the smsd package loads at run time. It is built on macOS for Apple silicon and never by the smsd
# build itself.
cmake_minimum_required(VERSION 3.18)
if(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "The architectures to build for")
endif()
project(smsd_metal LANGUAGES CXX OBJCXX)
include(GNUInstallDirs)

# The SMSD sources the package is built from. Their version is the smsd version it serves, and
# the package version must equal it.
set(SMSD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../cpp" CACHE PATH "The SMSD C++ sources (the cpp folder)")
if(NOT EXISTS "${SMSD_SOURCE_DIR}/include/smsd/screen_accelerator_abi.h"
   OR NOT EXISTS "${SMSD_SOURCE_DIR}/metal/metal_fingerprint.mm")
    message(FATAL_ERROR "smsd-metal: the SMSD sources were not found at ${SMSD_SOURCE_DIR}. "
        "Build it from an SMSD checkout, or set SMSD_SOURCE_DIR to its cpp folder.")
endif()
file(STRINGS "${SMSD_SOURCE_DIR}/CMakeLists.txt" _smsd_project REGEX "^project\\(smsd VERSION [0-9]+\\.[0-9]+\\.[0-9]+")
string(REGEX REPLACE "^project\\(smsd VERSION ([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" SMSD_VERSION "${_smsd_project}")
if(NOT SMSD_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "smsd-metal: the smsd version could not be read from ${SMSD_SOURCE_DIR}/CMakeLists.txt.")
endif()
if(DEFINED SKBUILD_PROJECT_VERSION AND NOT SKBUILD_PROJECT_VERSION STREQUAL SMSD_VERSION)
    message(FATAL_ERROR "smsd-metal ${SKBUILD_PROJECT_VERSION} must be built from the smsd "
        "${SKBUILD_PROJECT_VERSION} sources, and these are smsd ${SMSD_VERSION}.")
endif()

if(NOT APPLE)
    message(FATAL_ERROR "smsd-metal needs Apple Metal and is built on macOS only. "
        "The smsd package itself needs no Metal.")
endif()
if(NOT CMAKE_OSX_ARCHITECTURES STREQUAL "arm64")
    message(FATAL_ERROR "smsd-metal is built for Apple silicon (arm64) only, "
        "and CMAKE_OSX_ARCHITECTURES is ${CMAKE_OSX_ARCHITECTURES}.")
endif()
find_library(SMSD_METAL_FRAMEWORK Metal REQUIRED)
find_library(SMSD_FOUNDATION_FRAMEWORK Foundation REQUIRED)
message(STATUS "smsd-metal ${SMSD_VERSION}: ${CMAKE_OSX_ARCHITECTURES}, macOS ${CMAKE_OSX_DEPLOYMENT_TARGET}")

# The resident fingerprint screen and the device probe are the SMSD Metal sources; the package
# adds the resident label kernel and the table smsd attaches.
add_library(smsd_metal_screen SHARED
    src/screen_accelerator.mm
    "${SMSD_SOURCE_DIR}/metal/metal_batch.mm"
    "${SMSD_SOURCE_DIR}/metal/metal_fingerprint.mm")
target_include_directories(smsd_metal_screen PRIVATE "${SMSD_SOURCE_DIR}/include")
target_compile_definitions(smsd_metal_screen PRIVATE SMSD_METAL_SMSD_VERSION="${SMSD_VERSION}")
target_compile_features(smsd_metal_screen PRIVATE cxx_std_17)
target_compile_options(smsd_metal_screen PRIVATE -fobjc-arc -fobjc-arc-exceptions)
target_link_libraries(smsd_metal_screen PRIVATE ${SMSD_METAL_FRAMEWORK} ${SMSD_FOUNDATION_FRAMEWORK})
set_target_properties(smsd_metal_screen PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    OBJCXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON)
install(TARGETS smsd_metal_screen LIBRARY DESTINATION smsd_metal RUNTIME DESTINATION smsd_metal)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE"
    "${CMAKE_CURRENT_SOURCE_DIR}/NOTICE"
    DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/smsd-metal")

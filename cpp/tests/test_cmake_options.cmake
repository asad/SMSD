# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2018-2026 BioInception PVT LTD
# Algorithm Copyright (c) 2009-2026 Syed Asad Rahman
# See the NOTICE file for attribution, trademark, and algorithm IP terms.
foreach(value IN ITEMS OFF off 0 FALSE false NO no)
    execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${BINARY_DIR}/false-${value}"
        -DSMSD_BUILD_TESTS=OFF -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_OPENMP=OFF
        "-DSMSD_BUILD_CUDA=${value}" "-DSMSD_BUILD_METAL=${value}"
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
    if(NOT result EQUAL 0 OR NOT output MATCHES "CUDA disabled" OR NOT output MATCHES "Metal disabled|Metal not applicable")
        message(FATAL_ERROR "False accelerator option '${value}' was not honoured: ${output}${error}")
    endif()
endforeach()
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${BINARY_DIR}/lower-auto"
    -DSMSD_BUILD_TESTS=OFF -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_OPENMP=OFF
    -DSMSD_BUILD_CUDA=auto -DSMSD_BUILD_METAL=OFF
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "Lowercase CUDA auto did not configure: ${output}${error}")
endif()
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${SOURCE_DIR}" -B "${BINARY_DIR}/non-apple-true"
    -DCMAKE_SYSTEM_NAME=Linux -DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY
    -DSMSD_BUILD_TESTS=OFF -DSMSD_BUILD_PYTHON=OFF -DSMSD_BUILD_OPENMP=OFF
    -DSMSD_BUILD_CUDA=OFF -DSMSD_BUILD_METAL=TRUE
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(result EQUAL 0 OR NOT error MATCHES "Metal is only available")
    message(FATAL_ERROR "True Metal option on a non-Apple platform was not rejected: ${output}${error}")
endif()

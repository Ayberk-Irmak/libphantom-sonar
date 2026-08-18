# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required(VERSION 3.20)

project(phantom-sonar
        VERSION 0.14.0
        DESCRIPTION "Real-time ocean acoustics, ray tracing and sonar pulse analysis"
        LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(PHANTOM_BUILD_TESTS    "Build the unit test binary"                ON)
option(PHANTOM_BUILD_EXAMPLES "Build the example programs"                ON)
option(PHANTOM_BUILD_BENCH    "Build the hot-path benchmark"              ON)
option(PHANTOM_BUILD_TOOLS    "Build the Bellhop comparison tool"         ON)
option(PHANTOM_REAL_FLOAT     "Use float instead of double for Real"      OFF)
option(PHANTOM_NO_EXCEPTIONS "Build with -fno-exceptions -fno-rtti"          OFF)

# Reject near-miss option names rather than ignoring them.
#
# CMake accepts any -D silently, so a mistyped -DPHANTOM_USE_FLOAT=ON configures
# a perfectly ordinary double build and every test passes. That is worse than a
# failure: it produces a "verified at float precision" claim from a run that
# never used float. This happened during v0.10 development.
foreach(_typo PHANTOM_USE_FLOAT PHANTOM_FLOAT PHANTOM_SINGLE_PRECISION
              PHANTOM_USE_DOUBLE PHANTOM_TESTS PHANTOM_EXAMPLES
              PHANTOM_BENCH PHANTOM_SANITIZERS PHANTOM_WARNINGS_AS_ERRORS)
    if(DEFINED ${_typo})
        message(FATAL_ERROR
            "Unknown option ${_typo}. CMake would have ignored it and configured "
            "a default build. Did you mean one of: PHANTOM_REAL_FLOAT, "
            "PHANTOM_BUILD_TESTS, PHANTOM_BUILD_EXAMPLES, PHANTOM_BUILD_BENCH, "
            "PHANTOM_SANITIZE, PHANTOM_WERROR?")
    endif()
endforeach()
option(PHANTOM_WERROR         "Treat warnings as errors"                  OFF)
option(PHANTOM_SANITIZE       "Build with ASan + UBSan (tests only)"      OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Warning set. Deliberately strict: this library is meant to run in a control
# loop, so implicit conversions and shadowing are bugs waiting for sea trials.
# ---------------------------------------------------------------------------
if(MSVC)
    set(PHANTOM_WARNINGS /W4 /permissive-)
    if(PHANTOM_WERROR)
        list(APPEND PHANTOM_WARNINGS /WX)
    endif()
else()
    set(PHANTOM_WARNINGS
        -Wall -Wextra -Wpedantic
        -Wshadow -Wconversion -Wsign-conversion
        -Wold-style-cast -Wcast-qual -Wdouble-promotion
        -Wnull-dereference -Wformat=2)
    if(PHANTOM_WERROR)
        list(APPEND PHANTOM_WARNINGS -Werror)
    endif()
endif()

# ---------------------------------------------------------------------------
# The library
# ---------------------------------------------------------------------------
add_library(phantom_sonar
    src/ray_tracer.cpp
    src/channel.cpp
    src/coverage.cpp
    src/fft.cpp
    src/waveform.cpp
    src/matched_filter.cpp
    src/ping_analyzer.cpp
    src/echo_synth.cpp
    src/eigenray.cpp
    src/boundary.cpp
    src/reverberation.cpp
    src/array.cpp
    src/beamformer.cpp
    src/tracker.cpp
    src/imm.cpp
    src/ctrv.cpp
    src/comm.cpp
    src/c_api.cpp)
add_library(phantom::sonar ALIAS phantom_sonar)

target_include_directories(phantom_sonar
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)

target_compile_features(phantom_sonar PUBLIC cxx_std_20)
target_compile_options(phantom_sonar PRIVATE ${PHANTOM_WARNINGS})

if(PHANTOM_REAL_FLOAT)
    target_compile_definitions(phantom_sonar PUBLIC PHANTOM_REAL_FLOAT)
endif()

set_target_properties(phantom_sonar PROPERTIES
    OUTPUT_NAME phantom-sonar
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

# Freestanding / embedded build mode.
#
# The library never throws and never uses RTTI, so this should be a no-op -- and
# that is exactly why it is worth having as a build, because if it ever stops
# being true the failure is a link error here rather than an unwinder pulled
# into a firmware image that has no business carrying one.
#
# It is applied to the library only. The tests use no exceptions either, but the
# C++ standard library headers they include are not all exception-free, so the
# constraint that matters is the one on what ships.
if(PHANTOM_NO_EXCEPTIONS AND NOT MSVC)
    target_compile_options(phantom_sonar PRIVATE -fno-exceptions -fno-rtti)
    target_compile_definitions(phantom_sonar PUBLIC PHANTOM_NO_EXCEPTIONS=1)
endif()

# NOTE: -ffast-math is deliberately NOT enabled anywhere. It would let the
# compiler reassociate the Snell invariant and break bit-reproducibility across
# targets, which is exactly what a deterministic real-time library must not do.

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(PHANTOM_BUILD_TESTS)
    enable_testing()
    add_executable(phantom_tests
        tests/test_main.cpp
        tests/test_sound_speed.cpp
        tests/test_profile.cpp
        tests/test_ray_tracer.cpp
        tests/test_channel.cpp
        tests/test_fft.cpp
        tests/test_waveform.cpp
        tests/test_matched_filter.cpp
        tests/test_ping_analyzer.cpp
        tests/test_echo_synth.cpp
        tests/test_eigenray.cpp
        tests/test_boundary.cpp
        tests/test_reverberation.cpp
        tests/test_array.cpp
        tests/test_beamformer.cpp
        tests/test_tracker.cpp
        tests/test_real_profiles.cpp
        tests/test_imm.cpp
        tests/test_ctrv.cpp
        tests/test_comm.cpp)
    target_link_libraries(phantom_tests PRIVATE phantom_sonar)
    target_include_directories(phantom_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_compile_options(phantom_tests PRIVATE ${PHANTOM_WARNINGS})

    if(PHANTOM_SANITIZE AND NOT MSVC)
        target_compile_options(phantom_tests PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(phantom_tests PRIVATE -fsanitize=address,undefined)
    endif()

    add_test(NAME phantom_unit_tests COMMAND phantom_tests)

    # The C ABI, compiled by a C COMPILER in C11 mode -- not by C++ in C mode.
    # A header that only ever meets a C++ compiler quietly accumulates C++-isms
    # (a default argument, a `bool`, a scoped enum) and nobody finds out until a
    # real C caller arrives. This target is what finds out.
    enable_language(C)
    add_executable(phantom_c_tests tests/c/test_c_api.c)
    set_target_properties(phantom_c_tests PROPERTIES
        C_STANDARD 11
        C_STANDARD_REQUIRED ON
        C_EXTENSIONS OFF)
    target_link_libraries(phantom_c_tests PRIVATE phantom_sonar)
    if(NOT MSVC)
        target_compile_options(phantom_c_tests PRIVATE
            -Wall -Wextra -Wpedantic -Wstrict-prototypes -Wmissing-prototypes)
        if(PHANTOM_WERROR)
            target_compile_options(phantom_c_tests PRIVATE -Werror)
        endif()
        target_link_libraries(phantom_c_tests PRIVATE m)
    endif()
    if(PHANTOM_SANITIZE AND NOT MSVC)
        target_compile_options(phantom_c_tests PRIVATE -fsanitize=address,undefined)
        target_link_options(phantom_c_tests PRIVATE -fsanitize=address,undefined)
    endif()
    add_test(NAME phantom_c_abi_tests COMMAND phantom_c_tests)
endif()

# ---------------------------------------------------------------------------
# Examples and benchmark
# ---------------------------------------------------------------------------
if(PHANTOM_BUILD_EXAMPLES)
    add_executable(munk_simulation examples/munk_simulation.cpp)
    target_link_libraries(munk_simulation PRIVATE phantom_sonar)
    target_compile_options(munk_simulation PRIVATE ${PHANTOM_WARNINGS})

    add_executable(ping_intercept examples/ping_intercept.cpp)
    target_link_libraries(ping_intercept PRIVATE phantom_sonar)
    target_compile_options(ping_intercept PRIVATE ${PHANTOM_WARNINGS})

    add_executable(countermeasure_loop examples/countermeasure_loop.cpp)
    target_link_libraries(countermeasure_loop PRIVATE phantom_sonar)
    target_compile_options(countermeasure_loop PRIVATE ${PHANTOM_WARNINGS})
endif()

if(PHANTOM_BUILD_BENCH)
    add_executable(phantom_bench bench/bench_hotpath.cpp)
    target_link_libraries(phantom_bench PRIVATE phantom_sonar)
    target_compile_options(phantom_bench PRIVATE ${PHANTOM_WARNINGS})
endif()

if(PHANTOM_BUILD_TOOLS)
    # Reads a Bellhop .env file and traces it with this library, so both codes
    # consume the same input. See tools/bellhop_compare/README.md.
    add_executable(phantom_trace tools/bellhop_compare/phantom_trace.cpp)
    target_link_libraries(phantom_trace PRIVATE phantom_sonar)
    target_compile_options(phantom_trace PRIVATE ${PHANTOM_WARNINGS})
endif()

# ---------------------------------------------------------------------------
# Static analysis (optional targets, skipped silently when tools are absent)
# ---------------------------------------------------------------------------
find_program(PHANTOM_CLANG_TIDY NAMES clang-tidy)
if(PHANTOM_CLANG_TIDY)
    add_custom_target(tidy
        COMMAND ${PHANTOM_CLANG_TIDY}
                -p ${CMAKE_BINARY_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}/src/ray_tracer.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/channel.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/src/coverage.cpp
        COMMENT "Running clang-tidy")
endif()

find_program(PHANTOM_CPPCHECK NAMES cppcheck)
if(PHANTOM_CPPCHECK)
    add_custom_target(cppcheck
        COMMAND ${PHANTOM_CPPCHECK}
                --enable=warning,style,performance,portability
                --std=c++20 --error-exitcode=1 --inline-suppr --quiet
                -I ${CMAKE_CURRENT_SOURCE_DIR}/include
                ${CMAKE_CURRENT_SOURCE_DIR}/src
        COMMENT "Running cppcheck")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Install / package config
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS phantom_sonar
        EXPORT phantom-sonar-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/phantom DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT phantom-sonar-targets
        FILE phantom-sonar-targets.cmake
        NAMESPACE phantom::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phantom-sonar)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/phantom-sonar-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/phantom-sonar-config-version.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/phantom-sonar)

message(STATUS "phantom-sonar ${PROJECT_VERSION}  Real=${PHANTOM_REAL_FLOAT}  build=${CMAKE_BUILD_TYPE}")

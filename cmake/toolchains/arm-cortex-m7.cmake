# SPDX-License-Identifier: Apache-2.0
#
# Cortex-M7 (with double-precision FPU) cross-compilation.
#
#   cmake -S . -B build-m7 \
#         -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/arm-cortex-m7.cmake \
#         -DPHANTOM_BUILD_TESTS=OFF -DPHANTOM_BUILD_EXAMPLES=OFF -DPHANTOM_BUILD_BENCH=OFF
#
# A note on target selection: M0 and M3 have no FPU at all, so this library
# would run there only through soft-float emulation -- viable but slow enough
# that the real-time claims stop meaning anything. M4F has a single-precision
# FPU, so pair it with -DPHANTOM_REAL_FLOAT=ON and read docs/validation.md on
# what that costs in accuracy. M7 with -mfpu=fpv5-d16 is the target that runs
# the default double-precision build honestly.

set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(PHANTOM_ARM_PREFIX arm-none-eabi- CACHE STRING "cross toolchain prefix")

set(CMAKE_C_COMPILER   ${PHANTOM_ARM_PREFIX}gcc)
set(CMAKE_CXX_COMPILER ${PHANTOM_ARM_PREFIX}g++)
set(CMAKE_AR           ${PHANTOM_ARM_PREFIX}ar)
set(CMAKE_OBJCOPY      ${PHANTOM_ARM_PREFIX}objcopy)
set(CMAKE_SIZE         ${PHANTOM_ARM_PREFIX}size)

# No hosted runtime to link against when probing the compiler.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(PHANTOM_M7_FLAGS "-mcpu=cortex-m7 -mthumb -mfpu=fpv5-d16 -mfloat-abi=hard")

# -ffunction-sections/-fdata-sections let the linker drop unused engines.
# NOT -ffast-math: it would let the compiler reassociate the Snell invariant and
# break bit-reproducibility against the host build.
set(CMAKE_C_FLAGS_INIT   "${PHANTOM_M7_FLAGS} -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS_INIT "${PHANTOM_M7_FLAGS} -ffunction-sections -fdata-sections -fno-exceptions -fno-rtti")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

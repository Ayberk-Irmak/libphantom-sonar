# SPDX-License-Identifier: Apache-2.0
#
# Cortex-M7 with a single-precision FPU, the class of part an AUV payload
# computer actually carries.
#
#   cmake -S . -B build-m7 -DCMAKE_TOOLCHAIN_FILE=cmake/cortex-m7.cmake \
#         -DCMAKE_BUILD_TYPE=Release -DPHANTOM_REAL_FLOAT=ON \
#         -DPHANTOM_NO_EXCEPTIONS=ON
#
# PHANTOM_REAL_FLOAT is not optional advice here. The M7's FPU is single
# precision; double arithmetic is emulated in software and costs roughly two
# orders of magnitude, which turns a real-time budget into a non-real-time one
# without any error message. docs/validation.md quantifies what float costs in
# accuracy so the trade can be made deliberately.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER   arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_AR           arm-none-eabi-ar)
set(CMAKE_RANLIB       arm-none-eabi-ranlib)

# There is no runtime to link against when probing the compiler, so ask CMake
# for a static library instead of a test executable.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

set(PHANTOM_M7_FLAGS "-mcpu=cortex-m7 -mthumb -mfpu=fpv5-sp-d16 -mfloat-abi=hard -ffreestanding")
set(CMAKE_C_FLAGS_INIT   "${PHANTOM_M7_FLAGS}")
set(CMAKE_CXX_FLAGS_INIT "${PHANTOM_M7_FLAGS}")

# -ffunction-sections and -fdata-sections let the final link drop everything the
# firmware does not call, which is most of a library this size.
set(CMAKE_C_FLAGS_RELEASE_INIT   "-Os -ffunction-sections -fdata-sections")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Os -ffunction-sections -fdata-sections")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

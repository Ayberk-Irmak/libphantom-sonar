#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# Downloads and builds Bellhop (Acoustics Toolbox) for cross-validation.
#
# LICENSING: the Acoustics Toolbox is GPL-3.0. libphantom-sonar is Apache-2.0
# and does NOT link against it, include any of its code, or redistribute it.
# This script fetches it on your machine and `compare.py` runs bellhop.exe as a
# separate process, exchanging text files. No GPL obligations attach to
# libphantom-sonar as a result. Nothing from the toolbox is vendored into this
# repository.
#
#   ./tools/bellhop_compare/setup_bellhop.sh [install_dir]
#
# Default install dir: build/acoustics-toolbox
set -euo pipefail

AT_VERSION="at_2026_7"
AT_URL="http://oalib.hlsresearch.com/AcousticsToolbox/${AT_VERSION}.zip"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
INSTALL_DIR="${1:-${REPO_ROOT}/build/acoustics-toolbox}"

echo "==> Acoustics Toolbox setup"
echo "    version : ${AT_VERSION}"
echo "    target  : ${INSTALL_DIR}"

if ! command -v gfortran >/dev/null 2>&1; then
    cat >&2 <<'EOF'
error: gfortran not found. Install it first:

    Debian/Ubuntu/Kali : sudo apt-get install gfortran
    Fedora             : sudo dnf install gcc-gfortran
    macOS (Homebrew)   : brew install gcc
EOF
    exit 1
fi
echo "    fortran : $(gfortran --version | head -1)"

mkdir -p "${INSTALL_DIR}"
cd "${INSTALL_DIR}"

if [ ! -f "${AT_VERSION}.zip" ]; then
    echo "==> Downloading (~40 MB)"
    curl -fSL -o "${AT_VERSION}.zip" "${AT_URL}"
else
    echo "==> Archive already present, skipping download"
fi

if [ ! -d at ]; then
    echo "==> Extracting"
    unzip -q -o "${AT_VERSION}.zip"
fi

cd at

# Two deviations from the shipped Makefile, both needed on Linux:
#
#   -Wa,-q     is a macOS assembler flag and is rejected by GNU as.
#   serial     the build MUST NOT be parallelised. Fortran module files
#              (.mod) are compile-time dependencies that the Makefile does not
#              declare, so `make -j` races and fails with "Cannot open module
#              file". This costs a few minutes and is not optional.
echo "==> Building (serial; this takes a few minutes)"
FFLAGS_LINUX="-march=native -std=gnu -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -I../misc -I../tslib"

make clean >/dev/null 2>&1 || true
if ! make FC=gfortran FFLAGS="${FFLAGS_LINUX}" > build.log 2>&1; then
    echo "error: build failed; last 30 lines of build.log:" >&2
    tail -30 build.log >&2
    exit 1
fi

BELLHOP="${INSTALL_DIR}/at/Bellhop/bellhop.exe"
if [ ! -x "${BELLHOP}" ]; then
    echo "error: bellhop.exe was not produced; see ${INSTALL_DIR}/at/build.log" >&2
    exit 1
fi

cat <<EOF

==> Done.

    bellhop.exe : ${BELLHOP}

Run the cross-validation with:

    cmake -S . -B build && cmake --build build
    python3 tools/bellhop_compare/compare.py --bellhop "${BELLHOP}"

EOF

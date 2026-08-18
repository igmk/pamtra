#!/usr/bin/env bash
# Builds a private, static, single-threaded OpenBLAS and installs it (with a
# .pc file) into $PREFIX. Point meson at it with:
#
#   PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH" pip install .
#
# meson.build's `dependency('openblas', required: false)` call already tries
# pkg-config first, so this needs no meson.build changes -- it just has to
# win that pkg-config lookup ahead of any system/conda/Homebrew openblas.
#
# Why this exists: a normal dynamic OpenBLAS, if bundled into a PyPI wheel,
# would collide at runtime with numpy/scipy's own bundled OpenBLAS copy in
# the same process (duplicate global symbols, shared thread-pool state).
# Building it USE_THREAD=0 (PAMTRA only calls LAPACK for small per-particle
# T-matrix solves, not large GEMMs, so single-threaded is plenty fast) and
# statically with hidden visibility means its symbols are resolved at link
# time into pyPamtraLib's .so but never added to its exported dynamic symbol
# table -- so nothing else in the process can ever see or collide with this
# copy, regardless of what numpy/scipy bundle.
set -euo pipefail

VERSION="${OPENBLAS_VERSION:-0.3.34}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-"$ROOT_DIR/build-deps/openblas"}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

echo "Building static sequential OpenBLAS $VERSION -> $PREFIX"

curl -sL "https://github.com/OpenMathLib/OpenBLAS/releases/download/v${VERSION}/OpenBLAS-${VERSION}.tar.gz" \
  -o "$WORK_DIR/OpenBLAS.tar.gz"
tar -xzf "$WORK_DIR/OpenBLAS.tar.gz" -C "$WORK_DIR"

# CCOMMON_OPT/FCOMMON_OPT must be environment variables, not `make VAR=val`
# command-line args: Makefile.system appends its own required flags to them
# (e.g. -DMAX_PARALLEL_NUMBER=...) via `+=`, and GNU Make silently drops
# makefile `+=` appends to a variable that was set on the command line
# (only an `override` directive in the makefile could change it back, and
# Makefile.system doesn't have one) -- so a command-line CCOMMON_OPT locks
# out those required flags instead of just adding ours on top.
export CCOMMON_OPT="-fvisibility=hidden"
export FCOMMON_OPT="-fvisibility=hidden"

# BUILD_BFLOAT16=0: bfloat16 kernels are on by default on arm64, PAMTRA never
# needs them (only double-precision BLAS/LAPACK), and they fail to build
# under this toolchain (BGEMM_P/BGEMM_Q undeclared in gemm.c). This one is
# passed on the command line deliberately, to override Makefile.system's
# own arm64 default of BUILD_BFLOAT16=1.
make -C "$WORK_DIR/OpenBLAS-${VERSION}" -j"$JOBS" \
  USE_THREAD=0 \
  USE_LOCKING=1 \
  NO_SHARED=1 \
  BUILD_BFLOAT16=0

make -C "$WORK_DIR/OpenBLAS-${VERSION}" install PREFIX="$PREFIX"

# `make install` unconditionally creates libopenblas.{dylib,so}(.0) symlinks
# for the shared library even though NO_SHARED=1 means one was never built,
# leaving them dangling. A linker resolving `-lopenblas` could hit one of
# these before falling back to the .a, so remove them to force static
# linking deterministically.
find "$PREFIX/lib" -maxdepth 1 \( -name '*.dylib' -o -name '*.so*' \) ! -exec test -e {} \; -delete

echo "Done. Static lib + pkg-config file:"
find "$PREFIX" -name 'libopenblas*' -o -name 'openblas.pc'

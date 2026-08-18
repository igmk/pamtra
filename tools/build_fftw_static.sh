#!/usr/bin/env bash
# Builds a private, static FFTW and installs it (with a .pc file) into
# $PREFIX. Point meson at it with:
#
#   PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH" pip install .
#
# meson.build's `fftw3dep = dependency('fftw3')` call already tries
# pkg-config first, so this needs no meson.build changes -- it just has to
# win that pkg-config lookup ahead of any system/conda/Homebrew fftw3.
#
# Why this exists: for an eventual PyPI wheel, PAMTRA's own dependencies
# need to be self-contained. Unlike OpenBLAS, FFTW has no known symbol-
# collision risk with other commonly-bundled packages, so it doesn't need
# OpenBLAS's static+hidden-visibility treatment for correctness -- it's
# built static here purely to keep the wheel's runtime dependency surface
# minimal and avoid giving auditwheel/delocate one more shared object to
# chase down and bundle.
set -euo pipefail

VERSION="${FFTW_VERSION:-3.3.10}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-"$ROOT_DIR/build-deps/fftw"}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

echo "Building static FFTW $VERSION -> $PREFIX"

curl -sL "https://www.fftw.org/fftw-${VERSION}.tar.gz" -o "$WORK_DIR/fftw.tar.gz"
tar -xzf "$WORK_DIR/fftw.tar.gz" -C "$WORK_DIR"

cd "$WORK_DIR/fftw-${VERSION}"

# --enable-static is FFTW's default, spelled out here for clarity.
# --disable-shared: skip building the .so/.dylib entirely, so -lfftw3
# unambiguously resolves to the static .a (mirrors the dangling-symlink
# lesson from build_openblas_static.sh: without this, `make install` would
# still leave a shared-library placeholder around to trip up later).
# --with-pic: the static objects need to be relocatable to link into
# pyPamtraLib's own shared object.
# Fortran wrappers are intentionally left enabled (the default) --
# src/convolution.f90 includes FFTW's fftw3.f header directly, and that
# header's installation isn't worth risking over the tiny build-time
# saving from --disable-fortran.
./configure --prefix="$PREFIX" \
  --enable-static \
  --disable-shared \
  --with-pic

make -j"$JOBS"
make install

echo "Done. Static lib + pkg-config file:"
find "$PREFIX" -name 'libfftw3*' -o -name 'fftw3.pc'

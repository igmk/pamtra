#!/usr/bin/env bash
# Builds a private, dynamic HDF5 -> netCDF-C -> netCDF-Fortran stack into one
# shared $PREFIX. Point meson at it with:
#
#   PKG_CONFIG_PATH="$PREFIX/lib/pkgconfig:$PKG_CONFIG_PATH" pip install .
#
# meson.build's `dependency('netcdf')`/`dependency('netcdf-fortran')` calls
# already try pkg-config first, so this needs no meson.build changes.
#
# Unlike tools/build_openblas_static.sh and tools/build_fftw_static.sh, this
# stack is built as ordinary DYNAMIC libraries, not static. OpenBLAS needed
# static+hidden-visibility specifically to avoid colliding with another
# package's bundled copy in the same process (see that script's header);
# HDF5/netCDF have no such collision risk (no other commonly-bundled PyPI
# package embeds them), so the standard `auditwheel repair`/`delocate-wheel`
# path -- which detects a wheel's dynamic library dependencies and bundles
# them automatically -- is the right, lower-risk tool here. It's also what
# netCDF4's and h5py's own PyPI wheels already do. That repair step is a
# separate, later part of the cibuildwheel pipeline, not this script.
#
# zlib is intentionally NOT built here: HDF5 links against whatever zlib
# headers/lib configure finds on the build machine (manylinux images and
# macOS both ship one), and it's stable/ABI-compatible enough across
# versions that building our own adds little. auditwheel/delocate bundle
# (or, per manylinux policy, may skip as an allowed baseline lib) whichever
# zlib the build actually picked up.
#
# NOTE for local (non-wheel) testing after running this script: these are
# dynamic libraries installed to a non-standard prefix, so unlike the static
# OpenBLAS/FFTW builds, `import pyPamtra` needs that prefix on the loader
# search path at runtime:
#   macOS:  DYLD_LIBRARY_PATH="$PREFIX/lib" python3 -c "import pyPamtra"
#   Linux:  LD_LIBRARY_PATH="$PREFIX/lib" python3 -c "import pyPamtra"
# The eventual wheel doesn't need this -- auditwheel/delocate rewrite the
# built .so to load its own bundled copies via a relative rpath instead.
set -euo pipefail

HDF5_VERSION="${HDF5_VERSION:-1.14.6}"
NETCDF_C_VERSION="${NETCDF_C_VERSION:-4.10.1}"
NETCDF_FORTRAN_VERSION="${NETCDF_FORTRAN_VERSION:-4.6.4}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-"$ROOT_DIR/build-deps/netcdf"}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT

# Each stage's configure step compiles and runs small test programs linked
# against the previous stage's just-built shared library, before it's
# installed anywhere the OS would normally find it -- so both the loader
# path and the compiler/linker flags need to point at $PREFIX from the
# start, not just at the end.
export CPPFLAGS="-I$PREFIX/include"
export LDFLAGS="-L$PREFIX/lib"
export LD_LIBRARY_PATH="$PREFIX/lib:${LD_LIBRARY_PATH:-}"
export DYLD_LIBRARY_PATH="$PREFIX/lib:${DYLD_LIBRARY_PATH:-}"

echo "=== Building HDF5 $HDF5_VERSION -> $PREFIX ==="
curl -sL "https://github.com/HDFGroup/hdf5/releases/download/hdf5_${HDF5_VERSION}/hdf5-${HDF5_VERSION}.tar.gz" \
  -o "$WORK_DIR/hdf5.tar.gz"
tar -xzf "$WORK_DIR/hdf5.tar.gz" -C "$WORK_DIR"
(
  cd "$WORK_DIR/hdf5-${HDF5_VERSION}"
  ./configure --prefix="$PREFIX" \
    --enable-shared \
    --disable-static \
    --enable-fortran=no \
    --enable-cxx=no \
    --enable-build-mode=production
  make -j"$JOBS"
  make install
)

echo "=== Building netCDF-C $NETCDF_C_VERSION -> $PREFIX ==="
curl -sL "https://github.com/Unidata/netcdf-c/archive/refs/tags/v${NETCDF_C_VERSION}.tar.gz" \
  -o "$WORK_DIR/netcdf-c.tar.gz"
tar -xzf "$WORK_DIR/netcdf-c.tar.gz" -C "$WORK_DIR"
(
  cd "$WORK_DIR/netcdf-c-${NETCDF_C_VERSION}"
  # --disable-dap/--disable-byterange: PAMTRA never does remote/OPeNDAP
  # netCDF access (confirmed via grep across src/ and python/pyPamtra/),
  # so drop the libcurl dependency entirely rather than bundle it unused.
  # --disable-nczarr: PAMTRA only reads/writes classic netCDF files.
  ./configure --prefix="$PREFIX" \
    --enable-shared \
    --disable-static \
    --disable-dap \
    --disable-byterange \
    --disable-nczarr
  make -j"$JOBS"
  make install
)

echo "=== Building netCDF-Fortran $NETCDF_FORTRAN_VERSION -> $PREFIX ==="
curl -sL "https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v${NETCDF_FORTRAN_VERSION}.tar.gz" \
  -o "$WORK_DIR/netcdf-fortran.tar.gz"
tar -xzf "$WORK_DIR/netcdf-fortran.tar.gz" -C "$WORK_DIR"
(
  cd "$WORK_DIR/netcdf-fortran-${NETCDF_FORTRAN_VERSION}"
  ./configure --prefix="$PREFIX" \
    --enable-shared \
    --disable-static
  make -j"$JOBS"
  make install
)

echo "Done. Shared libs + pkg-config files:"
find "$PREFIX/lib" -maxdepth 1 -name '*hdf5*' -o -name '*netcdf*'
find "$PREFIX/lib/pkgconfig" -name '*.pc'

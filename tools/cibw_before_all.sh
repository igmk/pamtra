#!/usr/bin/env bash
# cibuildwheel `before-all` hook: builds the two bundled dependencies wheel
# builds actually need (static OpenBLAS, static FFTW) into a fixed location
# OUTSIDE the checked-out source tree. That matters: a prefix under the
# repo (the build scripts' own default) can trip meson.build's
# netcdf-fortran includedir handling on absolute paths that resolve inside
# the source tree -- see the meson.build commit that fixed this for local
# testing. Building outside the tree sidesteps the whole question rather
# than relying on that one fix being complete.
#
# tools/build_netcdf_stack.sh is NOT called here: wheel builds pass
# -Dbuild_cli=false (see [tool.cibuildwheel] in pyproject.toml), which
# drops the standalone pamtra CLI executable -- the only thing that needs
# netCDF-C/-Fortran/HDF5 at all (pyPamtraLib gets its NetCDF I/O purely
# through the Python netCDF4 package). That script still exists for local/
# dev use building the full CLI-included package.
#
# $PAMTRA_DEPS_PREFIX (set here, read by [tool.cibuildwheel.environment]'s
# PKG_CONFIG_PATH in pyproject.toml) must stay in sync with that setting.
set -euo pipefail

PAMTRA_DEPS_PREFIX="${PAMTRA_DEPS_PREFIX:-/tmp/pamtra-deps}"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [ "$(uname)" = "Darwin" ]; then
  # GitHub-hosted macOS runners ship gfortran already, but cibuildwheel's
  # manylinux-equivalent (a fresh container) does not -- brew install is
  # a harmless no-op if it's already present.
  brew install gcc
else
  # manylinux images are minimal containers with no Fortran compiler and
  # no zlib headers (only the runtime .so) by default.
  (yum install -y gcc-gfortran zlib-devel) || (dnf install -y gcc-gfortran zlib-devel)
fi

PREFIX="$PAMTRA_DEPS_PREFIX/openblas" "$ROOT_DIR/tools/build_openblas_static.sh"
PREFIX="$PAMTRA_DEPS_PREFIX/fftw" "$ROOT_DIR/tools/build_fftw_static.sh"

#!/usr/bin/env bash
# cibuildwheel `before-all` hook: builds all three bundled dependencies
# (static OpenBLAS, static FFTW, dynamic HDF5/netCDF-C/netCDF-Fortran) into
# a fixed location OUTSIDE the checked-out source tree. That matters: a
# prefix under the repo (the three build scripts' own default) can trip
# meson.build's netcdf-fortran includedir handling on absolute paths that
# resolve inside the source tree -- see the meson.build commit that fixed
# this for local testing. Building outside the tree sidesteps the whole
# question rather than relying on that one fix being complete.
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
  # manylinux images are minimal containers with no Fortran compiler, no
  # zlib headers (only the runtime .so), and no libxml2-config by default
  # (netCDF-C's configure wants the latter even with --disable-dap set --
  # some other feature, not just DAP, depends on it).
  (yum install -y gcc-gfortran zlib-devel libxml2-devel) || (dnf install -y gcc-gfortran zlib-devel libxml2-devel)
fi

PREFIX="$PAMTRA_DEPS_PREFIX/openblas" "$ROOT_DIR/tools/build_openblas_static.sh"
PREFIX="$PAMTRA_DEPS_PREFIX/fftw" "$ROOT_DIR/tools/build_fftw_static.sh"
PREFIX="$PAMTRA_DEPS_PREFIX/netcdf" "$ROOT_DIR/tools/build_netcdf_stack.sh"

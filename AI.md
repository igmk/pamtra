# AI.md

This file provides guidance to AI coding assistants (Claude Code, Codex, Gemini CLI, Cursor,
etc.) when working with code in this repository.

## What this is

PAMTRA (Passive and Active Microwave TRAnsfer) is a Fortran 90 radiative transfer model with a
Python wrapper (`pyPamtra`). It simulates radiometer and radar measurements of the cloudy
atmosphere. The core physics/RT engine lives in `src/` (Fortran), and is exposed to Python via
`f2py`-generated bindings (`pyPamtraLib`). There is also a standalone Fortran CLI binary (`pamtra`)
that reads namelist + profile files directly, for users who don't need Python.

## Build systems

**Meson (`meson.build`)** is the canonical build for normal use, driving both the pip-installable
`pyPamtra` package (`pyproject.toml` uses `mesonpy` as the build backend) and the standalone
`pamtra` Fortran CLI binary — both are declared as build targets in the same `meson.build` and
produced by the same `pip install .`.

**`Makefile`** is kept alongside it, but scoped down to one purpose: manual/HPC deployments (e.g.
DKRZ Levante via `install_levante_readmefirst.sh`) where hand-tuned linker flags against
cluster-specific module paths are easier to express as Makefile variables than to wire through
meson/pkg-config. Don't treat it as a second general-purpose build path — for everything else
(local dev, CI, packaging) use meson. When adding or removing a Fortran source file, update
**both**: the `OBJECTS` list in `Makefile` (order matters — dependencies must be listed before
dependents) and the `common_fsources`/`fsources`/`pamtra_fsources`/`f2py_sources` lists in
`meson.build`.

**Gotcha**: the Makefile compiles directly into `src/` (its `OBJDIR`/`SRCDIR` are both `src/` —
see the comment on that in `Makefile` — unlike meson, which builds in an isolated build dir). If
`make`/`make pamtra` has ever been run in a checkout, gfortran's subsequent meson build can pick up
those leftover `.mod`/`.o` files instead of building fresh ones (gfortran implicitly also searches
the directory of the source file it's compiling for module files), and fail with something like
`Cannot read module file '../src/foo.mod' ... created by a different version of GNU Fortran`. Fix:
`make clean` (or `rm -f src/*.o src/*.mod`) before building with meson/pip.

### Common commands

```bash
# Installs pyPamtra *and* puts the pamtra CLI binary on PATH (in the env's bin/)
pip install .

# Editable/dev install -- rebuilds the Python extension on change, but meson-python's
# editable installs don't install non-Python targets, so the pamtra CLI binary is NOT
# placed on PATH this way. Build/run it straight out of the build dir instead:
pip install --no-build-isolation -e . -Cbuild-dir=build
meson setup build --reconfigure && meson compile -C build   # -> build/pamtra
# (the pixi env wraps this as `pixi run build-cli`)

# Full install onto a cluster (DKRZ Levante) via the legacy Makefile, for
# reference on module/env setup and cluster-specific linker flags
bash install_levante_readmefirst.sh
```

### Testing

```bash
pip install ".[test]"
pytest tests/
```

`tests/` covers the pure-Python wrapper logic (descriptor file, meteoSI, namelist round-trips),
a couple of already-f2py-exposed Fortran routines called directly through the compiled
`pyPamtraLib` extension (`test_fortran_functions.py`), and golden-output regression tests
(`test_regression.py`) that run a small end-to-end scenario and compare against stored reference
arrays in `tests/golden/`. The Fortran core is otherwise built around shared `vars_*`/`settings`
module state rather than pure functions, so it isn't a realistic target for per-subroutine unit
tests — golden-output regression is the main safety net for it. To intentionally update the
reference values after a real physics change: `python tests/generate_golden_data.py`.

GitHub Actions CI (`.github/workflows/ci.yml`) runs the pip install + pytest above on Linux and
macOS, `make pamtra` on Linux (so the retained Makefile doesn't silently bit-rot), and a Sphinx
docs build, on every push/PR. Beyond that, correctness is checked via the
example scripts/notebooks in `examples/` (e.g. `examples/run_all_examples.py`,
`examples/pamtra_vs_pyPamtra.py`) and by comparing the standalone Fortran binary output against
the Python wrapper output.

### Required runtime environment

`pyPamtra` uses the `PAMTRA_DATADIR` environment variable to point at a directory containing
scattering databases and surface reflection catalogues. If it's set (`PAMTRA_DATADIR=""` counts
as set, and opts out of everything below -- several scattering methods and the built-in surface
emissivity defaults need no external data), it's used as-is. If it's unset entirely, `core.py`
auto-downloads and caches the data via [pamtra_data.py](python/pamtra_data.py) (`pooch`-based;
top-level, not inside the `pyPamtra` package, precisely so it's importable before `PAMTRA_DATADIR`
exists) and sets the env var itself -- `import pyPamtra` just works on a fresh machine, at the
cost of a one-time download on first import. If that download fails (no network, etc.), import
raises `RuntimeError` with next steps. The same fetch is available standalone via the
`pamtra-fetch-data` console script (see doc/source/installation.rst). The standalone Fortran CLI
binary has no such auto-fetch -- it always needs `PAMTRA_DATADIR` set manually.

External library dependencies for the Fortran build: LAPACK/BLAS (or OpenBLAS), FFTW3, NetCDF
(Fortran bindings), and a Fortran 90 compiler (gfortran assumed by both build systems).

`tools/build_openblas_static.sh` builds a private, static, single-threaded OpenBLAS with hidden
symbol visibility and installs its pkg-config file; point `meson.build`'s pkg-config-based
`dependency('openblas', ...)` lookup at it with `PKG_CONFIG_PATH=<prefix>/lib/pkgconfig pip install
.` (no meson.build changes needed for discovery). This exists for eventual PyPI wheel builds: a
normal dynamic OpenBLAS bundled into a wheel would collide at runtime with numpy/scipy's own
bundled copy in the same process. See [RELEASING.md](RELEASING.md) for the "why not PyPI" context
this is a building block for.

`tools/build_fftw_static.sh` and `tools/build_netcdf_stack.sh` are the same idea for the other two
dependencies -- static for FFTW (no collision risk like OpenBLAS, just kept static to give
auditwheel/delocate one less shared object to chase), dynamic for HDF5/netCDF-C/netCDF-Fortran
(large, complex chain where the standard auditwheel/delocate bundling path is lower-risk than
statically linking it ourselves; also what netCDF4's and h5py's own PyPI wheels do). Only the
standalone `pamtra` CLI executable needs netCDF at all -- `pyPamtraLib` (what `import pyPamtra`
loads) needs no direct netCDF link, since its NetCDF I/O goes through the pure-Python `netCDF4`
package instead. **Local testing gotcha**: unlike the static builds, these are ordinary dynamic
libraries installed to a non-standard prefix, so running anything against them locally (not
through a repaired wheel) needs `LD_LIBRARY_PATH`/`DYLD_LIBRARY_PATH` set to `<prefix>/lib` --
without it, the loader can silently resolve to a same-SONAME system copy instead (e.g. Debian's
own `libnetcdff.so.7` package) and produce corrupted output rather than an error. The eventual
wheel doesn't have this problem: `auditwheel`/`delocate` rewrite the built library to load its own
bundled copy via a relative rpath.

## Architecture

### Fortran core (`src/`)

~107 Fortran source files, compiled in dependency order (see `common_fsources` in `meson.build` for
the canonical ordering). Rough grouping:

- **Settings & I/O plumbing**: `settings.f90` (namelist parsing, global run settings),
  `parse_options.f90`/`getopt.f90` (CLI args), `descriptor_file.f90` (hydrometeor descriptor file
  parsing), `vars_atmosphere.f90`/`vars_output.f90`/`vars_rt.f90`/`vars_index.f90` (shared state
  passed between routines instead of function args — this is old-style Fortran global state, be
  careful about ordering/reentrancy when calling from Python), `report_module.f90` (error/warning
  reporting used throughout).
- **Gas absorption**: `gasabs_module.f90`, `rosen98_gasabs.f90`, `mpm93.f90`.
- **Hydrometeor microphysics / PSDs**: `drop_size_dist.f90`, `make_dist*.f90`, `calc_moment.f90`,
  `make_mass_size.f90`, `make_soft_spheroid.f90`, `vars_hydroFullSpec.f90` (explicit/"full
  spectrum" particle size distributions as an alternative to parameterized PSDs).
- **Scattering**: `mie_spheres.f90`, `tmatrix*.f90`/`tmatrix_lpq.f`/`tmatrix_amplq.lp.f` (T-matrix,
  largely translated legacy Fortran 77), `rayleigh_gans.f90`, `scatProperties.f90`,
  `dda_db_liu.f90`/`dda_db_hong.f90`/`hongdb.f90` (precomputed DDA scattering databases),
  `scatdb.c` (C module for one of the scattering DBs), `refractive_index.f90`, `eps_ice.f90`,
  `eps_water.f90`, `eps_mix.f90`.
- **Surface emissivity/reflection**: `sfc_optics.f90`, `land_sfc_optics.f90`,
  `ocean_sfc_optics.f90`, `tessem2.f90`, `telsem2.f90`, `fastemx.f90`,
  `large_scale_correction_module.f90`/`small_scale_correction_module.f90`,
  `foam_utility_module.f90`.
- **Radiative transfer solvers**: `rt4.f90`/`radtran4.f90`/`radintg4.f90`/`radscat4.f90` (the RT4
  discrete-ordinate solver), `run_rt.f90` (orchestrates a full RT call per profile/frequency).
- **Radar simulator**: `radar_simulator.f90`, `radar_spectrum.f90`,
  `radar_spectral_broadening.f90`, `radar_moments.f90`, `radar_hildebrand_sekhon.f90`,
  `dia2vel.f90`, `rescale_spectra.f90`.
- **Entry points**: `pamtra.f90` (standalone CLI `program pamtra` — reads namelist + profile file,
  loops over profiles/frequencies, writes NetCDF via `write_nc_results.f90`), `pyPamtraLib.f90`
  (the subroutine(s) f2py wraps to expose the same RT engine to Python, driven per-profile instead
  of via files).

The standalone binary and the Python path share almost all of the physics code but have separate
entry points and separate "fill settings" logic (`settings_read` from a namelist file for the CLI
vs `settings_fill_default` + direct attribute assignment from Python — see
`python/pyPamtra/libWrapper.py`). When changing shared global state in `vars_*` modules, check
both entry points.

### Python wrapper (`python/pyPamtra/`)

- `core.py` — the main `pyPamtra` class. Holds run settings (`self.set`), namelist settings
  (`self.nmlSet`), profile data, and dimensions; the class is the primary user-facing API.
- `libWrapper.py` — thin glue calling into the compiled `pyPamtraLib` f2py extension
  (`PamtraFortranWrapper`/`parallelPamtraFortranWrapper`); pushes Python-side settings into the
  Fortran module-level globals before invoking the RT engine, and includes a `multiprocessing`
  parallel wrapper for looping over profiles.
- `descriptorFile.py` — `pamDescriptorFile`: parses/builds the hydrometeor descriptor file format
  (one row per hydrometeor species, describing PSD shape, mass-size and area-size relations,
  scattering method, etc. — see `descriptorfiles/*.txt` for examples/format). This is a required
  input alongside the atmospheric profile for any RT run.
- `importer.py` — readers for various external model/profile formats (large file; grep for the
  specific format function rather than reading it all).
- `meteoSI.py` — meteorological unit conversions/thermodynamics helpers.
- `fortranNamelist/` — generic Fortran namelist reader/writer used to parse `*.nml` files
  independent of f2py.
- `plot.py`, `tools.py` — plotting helpers and misc utilities (e.g. `sftp2Cluster`).

### Descriptor files

Hydrometeor properties (ice, snow, rain, etc.) are configured via descriptor files
(`descriptorfiles/*.txt`), not hardcoded — each row defines one hydrometeor species' PSD
parameterization, mass-size/area-size power laws, and which scattering method/database to use.
Both the Fortran CLI (`descriptor_file.f90`) and Python (`descriptorFile.py`) parse this same text
format; keep them in sync if the format changes.

### Examples as documentation

`examples/` contains runnable scripts and notebooks that double as integration tests/usage docs
(e.g. `pamtra_vs_pyPamtra.py` cross-checks the standalone binary against the Python wrapper on the
same input). When changing the public Python API or the descriptor file format, check these for
breakage.

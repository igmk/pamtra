# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

PAMTRA (Passive and Active Microwave TRAnsfer) is a Fortran 90 radiative transfer model with a
Python wrapper (`pyPamtra`). It simulates radiometer and radar measurements of the cloudy
atmosphere. The core physics/RT engine lives in `src/` (Fortran), and is exposed to Python via
`f2py`-generated bindings (`pyPamtraLib`). There is also a standalone Fortran CLI binary (`pamtra`)
that reads namelist + profile files directly, for users who don't need Python.

## Build systems

There are **two parallel build systems** — know which one you're touching:

1. **Meson (`meson.build`)** — the modern, canonical build, used for the pip-installable
   `pyPamtra` package (`pyproject.toml` uses `mesonpy` as the build backend). This is what
   `pip install .` uses.
2. **`Makefile`** — the legacy/manual build, still used for building the standalone `pamtra`
   Fortran binary and for quick local iteration without going through meson-python.

When adding or removing a Fortran source file, **update both**: the `OBJECTS` list in `Makefile`
(order matters — dependencies must be listed before dependents) and the `fsources`/`f2py_sources`
lists in `meson.build`.

### Common commands

```bash
# Modern build (installs pyPamtra into the current Python env)
pip install .

# Editable/dev install (rebuilds Fortran extension on change, useful while iterating)
pip install --no-build-isolation -e . -Cbuild-dir=build

# Legacy Makefile build: Fortran binary + Python extension
make            # builds both `bin/pamtra` and the pyPamtraLib.so (equivalent to `make pamtra py`)
make pamtra     # standalone Fortran CLI binary only -> bin/pamtra
make py         # python extension only -> python/pyPamtra/pyPamtraLib*.so
make pamtraDebug   # debug build with -g -fbacktrace -fbounds-check (for gdb/valgrind)
make clean      # remove build artifacts

# Full install onto a cluster (DKRZ Levante), for reference on module/env setup
bash install_levante_readmefirst.sh
```

There is no local test suite (no `pytest`/`unittest` files in the repo) and no CI workflow
configured. Correctness is currently checked via the example scripts/notebooks in `examples/`
(e.g. `examples/run_all_examples.py`, `examples/pamtra_vs_pyPamtra.py`) and by comparing the
standalone Fortran binary output against the Python wrapper output.

### Required runtime environment

`pyPamtra` requires the `PAMTRA_DATADIR` environment variable to point at a directory containing
scattering databases and surface reflection catalogues (downloaded separately, see `readme.md`).
Importing `pyPamtra` without it set raises `RuntimeError` at import time
([python/pyPamtra/core.py](python/pyPamtra/core.py)).

External library dependencies for the Fortran build: LAPACK/BLAS (or OpenBLAS), FFTW3, NetCDF
(Fortran bindings), and a Fortran 90 compiler (gfortran assumed by both build systems).

## Architecture

### Fortran core (`src/`)

~107 Fortran source files, compiled in dependency order (see `OBJECTS` in `Makefile` for the
canonical ordering). Rough grouping:

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

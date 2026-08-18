project( 'pyPamtra', 'c',
  # Note that the git commit hash cannot be added dynamically here
  version: '1.0.3',
  license: 'GPL v3',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=debugoptimized',
    # TODO: the below -Wno flags are all needed to silence warnings in
    # f2py-generated code. This should be fixed in f2py itself.
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',  # silence "conversion from REAL(8) to INTEGER(4)"
    #'fortran_std=legacy',
  ],
)

# Adding at project level causes many spurious -lgfortran flags.
add_languages('fortran', native: false)

fc = meson.get_compiler('fortran')
cc = meson.get_compiler('c')

# This argument is called -Wno-unused-but-set-variable by GCC, however Clang
# doesn't recognize that.
if cc.has_argument('-Wno-unused-but-set-variable')
  add_global_arguments('-Wno-unused-but-set-variable', language : 'c')
endif

# For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

# Add more link arguments
add_project_link_arguments('-lquadmath', language: ['c', 'fortran'])
add_project_link_arguments('-lz', language: ['c', 'fortran'])

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation()
py3_dep = py3.dependency()
message(py3.full_path())
message(py3.get_install_dir())

incdir_numpy = run_command(py3,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py3,
  ['-c', 'import os; os.chdir(".."); from numpy import f2py; print(f2py.get_include())'],
  check : true
).stdout().strip()

inc_dirs = include_directories(incdir_numpy, incdir_f2py)

# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in a meson.build file::
#
#   py3.extension_module('_name',
#     'source_fname',
#     numpy_nodepr_api)
#
c_flags = ['-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION']

# Platform detection to set more flags for Windows systems
is_windows = host_machine.system() == 'windows'
is_mac = host_machine.system() == 'darwin'
is_mingw = is_windows and cc.get_id() == 'gcc'

cython_c_args = []
if is_windows
  # For mingw-w64, link statically against the UCRT.
  gcc_link_args = ['-lucrt', '-static']
  if is_mingw
    add_global_link_arguments(gcc_link_args, language: ['c'])

    # Force gcc to float64 long doubles for compatibility with MSVC
    # builds, for C only.
    add_global_arguments('-mlong-double-64', language: 'c')

    # Make fprintf("%zd") work (see https://github.com/rgommers/scipy/issues/118)
    add_global_arguments('-D__USE_MINGW_ANSI_STDIO=1', language: ['c'])

    # Manual add of MS_WIN64 macro when not using MSVC.
    # https://bugs.python.org/issue28267
    bitness = run_command('_gcc_build_bitness.py').stdout().strip()
    if bitness == '64'
      add_global_arguments('-DMS_WIN64', language: ['c', 'fortran'])
    endif

    # Silence warnings emitted by PyOS_snprintf for (%zd), see
    # https://github.com/rgommers/scipy/issues/118.
    # Use as c_args for extensions containing Cython code
    cython_c_args += ['-Wno-format-extra-args', '-Wno-format']

    # Windows may need the directory for Python.h added.  This location
    # does not have the same name for all installs.  This is the one
    # for Windows on GitHub actions
    incdir_py = run_command(py3,
      ['-c', 'import os; os.chdir(".."); import sys; print(os.path.join(sys.prefix, "include"))'],
      check : true).stdout().strip()
    inc_dirs = include_directories(incdir_numpy, incdir_f2py, incdir_py)
  endif

  if meson.get_compiler('fortran').get_id() == 'gcc'
    add_global_link_arguments(gcc_link_args, language: ['fortran'])

    # Flag needed to work around BLAS and LAPACK Gfortran dependence on
    # undocumented C feature when passing single character string
    # arguments.
    # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90329
    # https://github.com/wch/r-source/blob/838f9d5a7be08f2a8c08e47bcd28756f5d0aac90/src/gnuwin32/MkRules.rules#L121
    add_global_arguments('-fno-optimize-sibling-calls',
      language: ['fortran'])
  endif
elif is_mac
  c_flags += ['-FALIGN_FUNCTIONS=8']
endif

# Check the python headers
cc.check_header('Python.h', dependencies: [py3_dep], required: true)

# Unlike distutils, meson doesn't yet include some of the f2py stuff
fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  c_args: c_flags,
  dependencies: py3_dep,
  include_directories: [incdir_numpy, incdir_f2py])

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [incdir_numpy, incdir_f2py])

fftw3dep = dependency('fftw3')
#blasdep = dependency('blas')
#lapackdep = dependency('lapack')

fs = import('fs')
openblasdep = dependency('openblas', required: false)

if not openblasdep.found()
  # A conda/pixi environment's own libopenblas package doesn't ship a
  # pkg-config file, but it must be preferred over any system copy: e.g.
  # linking against Homebrew's openblas from inside a conda/pixi build
  # pulls in a second, different OpenMP runtime that aborts at dlopen
  # time when mixed with the env's own libraries. sys.prefix is
  # authoritative for whatever Python environment is actually building.
  py_prefix = run_command(py3, ['-c', 'import sys; print(sys.prefix)'], check: true).stdout().strip()
  foreach libname : ['libopenblas.dylib', 'libopenblas.so', 'libopenblas.0.dylib', 'libopenblas.so.0']
    candidate = py_prefix / 'lib' / libname
    if not openblasdep.found() and fs.exists(candidate)
      # conda-forge's libopenblas has an @rpath-relative install name (unlike
      # e.g. Homebrew's, which uses an absolute one), so a raw path here
      # isn't enough -- dyld also needs a matching -rpath at link time, or
      # the binary fails to even start with "Library not loaded: @rpath/...".
      # This can't rely on some *other* dependency's pkg-config-derived
      # rpath happening to cover the same directory: on a system where
      # netcdf/fftw resolve to Homebrew but openblas falls back to this
      # env, their lib dirs differ.
      openblasdep = declare_dependency(link_args: [candidate, '-Wl,-rpath,' + py_prefix / 'lib'])
    endif
  endforeach
endif

if not openblasdep.found() and is_mac
  # openblas is keg-only in Homebrew (macOS ships BLAS/LAPACK via the
  # Accelerate framework instead), so its .pc file is not symlinked into
  # the default pkg-config search path and dependency('openblas') fails
  # unless the user manually exports PKG_CONFIG_PATH. Only reached when
  # there's no env-local openblas above (plain `pip install .` with
  # Homebrew-installed dependencies, not a conda/pixi build).
  brew_prefix = run_command('brew', ['--prefix', 'openblas'], check: false)
  if brew_prefix.returncode() == 0
    openblas_root = brew_prefix.stdout().strip()
    openblasdep = declare_dependency(
      include_directories: include_directories(openblas_root / 'include'),
      link_args: ['-L' + openblas_root / 'lib', '-lopenblas', '-Wl,-rpath,' + openblas_root / 'lib'],
    )
  endif
endif
if not openblasdep.found()
  error('openblas not found. On macOS, install it with `brew install openblas`; ' +
        'on other platforms, make sure its pkg-config file is on PKG_CONFIG_PATH.')
endif

netcdfdep = dependency('netcdf')

# Best-effort early warning for the same class of problem as the
# netcdf-fortran fallback below: if pkg-config resolved `netcdf` to a
# copy outside the active Python environment even though that
# environment has its own netcdf.pc, pkg-config itself is most likely
# not installed in the environment, so PATH fell back to a system/
# Homebrew pkg-config that found a different, possibly ABI-incompatible
# netcdf (e.g. linked against a different libcurl/libhdf5 than the rest
# of the environment). That mismatch otherwise only surfaces much later,
# as a confusing "undefined reference" error at link time.
netcdf_env_check_prefix = run_command(py3, ['-c', 'import sys; print(sys.prefix)'], check: true).stdout().strip()
netcdf_own_pc = netcdf_env_check_prefix / 'lib' / 'pkgconfig' / 'netcdf.pc'
netcdf_libdir = netcdfdep.get_variable(pkgconfig: 'libdir', default_value: '')
if fs.exists(netcdf_own_pc) and netcdf_libdir != '' and not netcdf_libdir.startswith(netcdf_env_check_prefix)
  warning('''pkg-config resolved `netcdf` to @0@, outside the active Python environment (@1@), even though that environment has its own copy at @2@. This usually means pkg-config itself is not installed in the environment, so a system/Homebrew pkg-config was used instead. This can cause confusing link errors later (e.g. undefined curl symbols) if the two netcdf/hdf5 builds are ABI-incompatible. Fix: install pkg-config into the active environment, e.g. `conda install pkg-config`.'''.format(netcdf_libdir, netcdf_env_check_prefix, netcdf_own_pc))
endif

netcdffdep = dependency('netcdf-fortran', required: false)

if netcdffdep.found()
  # pkg-config drops "-I/usr/include" from --cflags when that's the
  # resolved includedir, since it's a default search path for the C
  # compiler -- but gfortran, unlike cc, does NOT search /usr/include
  # for .mod files unless told to explicitly. On Debian/Ubuntu (where
  # netcdf-fortran.pc's includedir is exactly /usr/include, e.g. via
  # apt's libnetcdff-dev), that silently drops the one flag gfortran
  # needs to find netcdf.mod, breaking the `pamtra` executable target
  # below with "Fatal Error: Cannot open module file 'netcdf.mod'"
  # even though pkg-config resolution "succeeded". Re-add it by hand,
  # using get_variable() (unaffected by the --cflags stripping) as the
  # source of truth for where the .mod files actually live.
  netcdff_includedir = netcdffdep.get_variable(pkgconfig: 'includedir', default_value: '')
  if netcdff_includedir != ''
    netcdffdep = declare_dependency(
      dependencies: netcdffdep,
      include_directories: include_directories(netcdff_includedir),
    )
  endif
endif

if not netcdffdep.found()
  # Same class of problem as the openblas fallback above: the active
  # Python environment may have netcdf-fortran installed (e.g. via
  # `conda install netcdf-fortran`) without its .pc file ever being found
  # -- either because that environment has no pkg-config package of its
  # own, so PATH falls back to a system/Homebrew pkg-config that doesn't
  # know about the env's lib/pkgconfig at all, or PKG_CONFIG_PATH just
  # isn't exported. sys.prefix is authoritative for whatever Python
  # environment is actually building.
  py_prefix = run_command(py3, ['-c', 'import sys; print(sys.prefix)'], check: true).stdout().strip()
  foreach libname : ['libnetcdff.dylib', 'libnetcdff.so']
    candidate = py_prefix / 'lib' / libname
    if not netcdffdep.found() and fs.exists(candidate)
      # See the matching comment on the openblas fallback above: this needs
      # its own -rpath too, since conda-forge's libnetcdff also uses an
      # @rpath-relative install name.
      netcdffdep = declare_dependency(
        include_directories: include_directories(py_prefix / 'include'),
        link_args: [candidate, '-Wl,-rpath,' + py_prefix / 'lib'],
      )
    endif
  endforeach
endif

if not netcdffdep.found() and is_mac
  brew_prefix = run_command('brew', ['--prefix', 'netcdf-fortran'], check: false)
  if brew_prefix.returncode() == 0
    netcdff_root = brew_prefix.stdout().strip()
    netcdffdep = declare_dependency(
      include_directories: include_directories(netcdff_root / 'include'),
      link_args: ['-L' + netcdff_root / 'lib', '-lnetcdff', '-Wl,-rpath,' + netcdff_root / 'lib'],
    )
  endif
endif

if not netcdffdep.found()
  error('netcdf-fortran not found. Install it (e.g. `conda install -c conda-forge netcdf-fortran` ' +
        'or `brew install netcdf-fortran`), and either make sure pkg-config itself is installed ' +
        'in the active environment, or export PKG_CONFIG_PATH to include its lib/pkgconfig dir.')
endif

# List sources
csources = ['src/scatdb.c']
# Shared by the pyPamtraLib python extension and the standalone `pamtra`
# CLI binary -- their own entry points/glue (pyPamtraLib.f90 vs.
# pamtra.f90 + write_nc_results.f90) are added separately below.
common_fsources = ['src/kinds.f90',
            'src/mod_fastem4_coef.f90', # seems unused
            'src/zlib_stuff.f90',
            'src/mod_io_strings.f90',
            'src/rho_air.f90',
#            'src/versionNumber.auto.f90', # this needs to be created by a script
            'src/getopt.f90',#1
            'src/random.f90',
            'src/e_sat_gg_water.f90',
            'src/equcom.f90',   # TODO check used by
            'src/equare.f90',   # TODO check used by
            'src/vars_index.f90',
            'src/rt_utilities.f90',
            'src/nan.f90',
            'src/constants.f90', #2
            'src/conversions.f90',
            'src/dsort.f90',
            'src/report_module.f90',
            'src/smooth_savitzky_golay.f90',
            'src/radar_hildebrand_sekhon.f90',
            'src/viscosity_air.f90',
            'src/interpolation.f90',
            'src/gasabs_module.f90',
            'src/rosen98_gasabs.f90',
            'src/eps_ice.f90',
            'src/eps_water.f90',
            'src/ref_ice.f90',
            'src/ref_water.f90',
            'src/eps_mix.f90',
            'src/refractive_index.f90', # seems unused
            'src/mie_scat_utilities.f90',
            'src/scat_utilities.f90',
            'src/mpm93.f90',
            'src/mod_mwatlas_nt_bin.f90',
            'src/surface.f90',
            'src/tessem2.f90',
            'src/liu.f90',
            'src/fresnel.f90',
            'src/hyperbolic_step.f90',
            'src/slope_variance.f90',
            'src/reflection_correction_module.f90',
            'src/large_scale_correction_module.f90',
            'src/small_scale_correction_module.f90',
            'src/foam_utility_module.f90',
            'src/azimuth_emissivity_module.f90',
            'src/descriptor_file.f90',
            'src/settings.f90',
            'src/land_emis_ssmi.f90', # uses equcom equare
            #'src/land_emis.f90', # dropped by older versions already
            'src/dia2vel.f90',
            'src/mie_spheres.f90',
            'src/dda_db_liu.f90',
            'src/hongdb.f90',
            'src/dda_db_hong.f90',
            'src/telsem2.f90',
            'src/fastemx.f90',
            'src/convolution.f90',
            'src/radmat.f90',
            'src/parse_options.f90',
            'src/vars_atmosphere.f90',
            'src/vars_rt.f90',
            'src/vars_hydroFullSpec.f90',
            'src/vars_output.f90',
            'src/collect_output.f90',
            'src/save_active.f90',
            'src/get_gasabs.f90',
            'src/sfc_matrices.f90',
            'src/land_sfc_optics.f90',
            'src/ocean_sfc_optics.f90',
            'src/sfc_optics.f90',
            'src/drop_size_dist.f90',
            'src/make_dist.f90',
            'src/make_dist_param.f90',
            'src/make_mass_size.f90',
            'src/calc_moment.f90',
            'src/make_soft_spheroid.f90',
            'src/check_print.f90',
            'src/tmatrix_lpq.f',
            'src/tmatrix_amplq.lp.f',
            'src/tmatrix.f90',
            'src/rayleigh_gans.f90',
            'src/run_rt.f90',
            'src/scatProperties.f90', # this gives compiling errors # -fallow-argument-mismatch? 
            'src/hydrometeor_extinction.f90',
            'src/rescale_spectra.f90',
            'src/radar_moments.f90',
            'src/radar_spectrum.f90',
            'src/radar_spectral_broadening.f90',
            'src/radar_simulator.f90',
            'src/radar_hildebrand_sekhon.f90',
            'src/get_scat_mat.f90',
            'src/radscat4.f90',
            'src/radtran4.f90',
            'src/radintg4.f90',
            'src/scatcnv.f90',
            'src/rt4.f90',
            'src/deallocate_everything.f90', # is this really necessary? only pyPamtra stuff
            ]

fsources = common_fsources + ['src/pyPamtraLib.f90']
pamtra_fsources = common_fsources + ['src/write_nc_results.f90']

f2py_sources = ['src/report_module.f90',
                'src/vars_index.f90',
                'src/viscosity_air.f90',
                'src/convolution.f90',
                'src/deallocate_everything.f90',
                'src/vars_output.f90',
                'src/vars_atmosphere.f90',
                'src/settings.f90',
                'src/descriptor_file.f90',
                'src/vars_hydroFullSpec.f90',
                'src/radar_moments.f90',
                'src/eps_water.f90',
                'src/radar_hildebrand_sekhon.f90',
                'src/dia2vel.f90',
                'src/pyPamtraLib.f90'
                ]

versionNumberAuto = custom_target(
  'versionNumber.auto.f90',
  output : 'versionNumber.auto.f90',
  input : '_version_creator.py',
  command : [py3, '@INPUT@'])



# Make a custom target, output should not contain a path segment.
# F2PY controls output names, so it may change at some point.
fortran_pamtra_source = custom_target(
   'interface_module',
   input : [f2py_sources],
   output : ['pyPamtraLib-f2pywrappers2.f90','pyPamtraLibmodule.c'],
   command : [py3, '-m', 'numpy.f2py', '@INPUT@', '-m', 'pyPamtraLib', '--lower', '--f2cmap', '../src/.f2py_f2cmap'])


extra_fflags_list = ['-fallow-argument-mismatch', '-fallow-invalid-boz'] # TODO check if these are really needed for modern gfortran
extra_fflags = ['-cpp']
foreach ff : extra_fflags_list
  if fc.has_argument(ff) # older gfortran do not have certain flags
    extra_fflags += ff
  endif
endforeach

# Restrict pyPamtraLib's exported dynamic symbol table to just the one
# symbol CPython's import machinery actually needs (dlopen + dlsym for
# PyInit_pyPamtraLib). Nothing else in the process ever looks up symbols in
# this .so by name, so this is safe -- and it matters for OpenBLAS: if
# another package in the same process (e.g. numpy/scipy) bundles its own
# copy of OpenBLAS, two same-named global symbols (dgemm_, thread-pool
# state, ...) loaded into one process can collide. `-fvisibility=hidden`
# passed to OpenBLAS's own build (see tools/build_openblas_static.sh)
# handles most of it, but not all: several of its ARM64 kernels are
# hand-written assembly with their own `.globl` directives, which ignore
# C-compiler visibility flags entirely and still end up global in
# libopenblas.a. A link-time whitelist here closes that gap regardless of
# where a leak comes from.
pyPamtraLib_link_args = []
pyPamtraLib_link_depends = []
if is_mac
  pyPamtraLib_link_args = ['-Wl,-exported_symbol,_PyInit_pyPamtraLib']
elif not is_windows
  pyPamtraLib_map = meson.current_source_dir() / 'src/pyPamtraLib.map'
  pyPamtraLib_link_args = ['-Wl,--version-script=' + pyPamtraLib_map]
  pyPamtraLib_link_depends = [pyPamtraLib_map]
endif
# Windows is intentionally left alone: PAMTRA has no Windows build to test
# this against (see RELEASING.md), and MSVC/MinGW handle symbol export via
# .def files / __declspec, not GNU-ld version scripts or ld64 flags.

# Declare the fortran extension module
py3.extension_module('pyPamtraLib', # extension module should have the same name ftmatrix as the target of f2py to be linked... at least when build with meson
  [fsources, csources, fortran_pamtra_source, fortranobject_c, versionNumberAuto],
  c_args: c_flags,
  fortran_args: extra_fflags,
  include_directories: inc_dirs,
  link_with: fortranobject_lib,
  link_args: pyPamtraLib_link_args,
  link_depends: pyPamtraLib_link_depends,
  dependencies : [py3_dep, fortranobject_dep, fftw3dep, openblasdep, netcdfdep],
  subdir: 'pyPamtra/',
  install : true)

# Standalone Fortran CLI binary (reads namelist + profile files directly,
# writes NetCDF output) -- shares almost all physics code with pyPamtraLib
# above, but has its own entry point and NetCDF-writing routine.
executable('pamtra',
  ['src/pamtra.f90', pamtra_fsources, csources, versionNumberAuto],
  c_args: c_flags,
  fortran_args: extra_fflags,
  dependencies : [fftw3dep, openblasdep, netcdfdep, netcdffdep],
  # Link with the Fortran driver (gfortran), not the project-default C
  # driver: the T-matrix code uses quad precision, and only gfortran's
  # own link spec pulls in the libgcc routines (__multf3 etc.) needed for
  # that -- cc doesn't add them, so a plain executable() here fails to
  # link (unlike the pyPamtraLib bundle, which hides the same gap behind
  # -undefined dynamic_lookup).
  link_language: 'fortran',
  install : true)

# Declare and copy python sources
py3.install_sources([
  'python/pyPamtra/__init__.py',
  'python/pyPamtra/core.py',
  'python/pyPamtra/descriptorFile.py',
  'python/pyPamtra/importer.py',
  'python/pyPamtra/libWrapper.py',
  'python/pyPamtra/meteoSI.py',
  'python/pyPamtra/pamtra.py',
  'python/pyPamtra/plot.py',
  'python/pyPamtra/tools.py',
  ],
  pure: false,
  subdir: 'pyPamtra/')

py3.install_sources([
  'python/pyPamtra/fortranNamelist/__init__.py',
  'python/pyPamtra/fortranNamelist/namelist.py',
  'python/pyPamtra/fortranNamelist/namelist_settings.py',
  ],
  pure: false,
  subdir: 'pyPamtra/fortranNamelist/')

# Top-level (not under pyPamtra/) so it can be imported/run before
# PAMTRA_DATADIR is set -- pyPamtra/__init__.py raises at import time if
# it isn't, so a helper that sets it up can't live inside that package.
py3.install_sources(['python/pamtra_data.py'], pure: true)

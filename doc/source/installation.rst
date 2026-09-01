..  _installation:


Installation
============

PAMTRA is installed with pip, which compiles the Fortran core and builds both
the ``pyPamtra`` Python extension and the standalone ``pamtra`` command line
binary (:ref:`pamtra`) in one step via `meson-python
<https://mesonbuild.com/meson-python/>`_. This has replaced the old
``make`` / ``make pyinstall`` workflow.

.. note::
   A regular ``pip install .`` puts the ``pamtra`` binary in the same
   ``bin/`` directory as ``python``/``pip`` themselves, so it's already on
   ``PATH`` whenever that environment is active -- no extra step needed. An
   *editable* install (``pip install -e .``, see below) does **not** install
   it, since meson-python's editable-install support only covers the Python
   extension; build/run it straight out of the build directory instead (see
   :ref:`pamtra`).

.. warning::
   If this checkout was ever built with the legacy ``Makefile`` (``make`` /
   ``make pamtra``, still used for HPC deployments -- see below), run
   ``make clean`` first. The Makefile compiles directly into ``src/``, and a
   subsequent ``pip install .`` can pick up those leftover ``.mod`` files
   instead of building fresh ones, failing with something like
   ``Cannot read module file '../src/foo.mod' ... created by a different
   version of GNU Fortran``.


pip install (recommended, prebuilt wheels)
*********************************************

For **Linux (x86_64)** and **macOS (Apple Silicon / arm64) only**::

    pip install pamtra

This installs a self-contained wheel with FFTW and OpenBLAS already bundled
in -- no system libraries, compiler, or conda/pixi environment needed. It
does **not** include the standalone ``pamtra`` CLI binary (:ref:`pamtra`),
which needs netCDF-Fortran (not bundled into the wheel); use one of the
from-source installs below if you need it.

Not available for **Windows** or **Intel macOS** (``osx-64``) -- no wheels
are built for either platform. Use conda-forge below (covers Intel macOS,
and includes the CLI binary) or WSL2 below (Windows) instead.

.. warning::
   If you previously used the legacy ``make pyinstall`` workflow (from
   before PAMTRA switched to ``pip``/meson-python), it copied the Python
   package to ``~/lib/python/pyPamtra`` -- a path with no meaning to pip
   itself, only useful if manually added to ``PYTHONPATH`` (a common setup
   on shared/HPC systems without permission to install into the environment
   directly). If that's still on your ``PYTHONPATH``, it will **shadow** a
   freshly ``pip install``ed PAMTRA: Python resolves ``import pyPamtra`` to
   whichever copy comes first on ``sys.path``, and a stale
   ``PYTHONPATH``-injected directory typically wins over the environment's
   own site-packages. This silently reintroduces old bugs/missing features
   with no error to explain why. Fix: remove ``~/lib/python/pyPamtra`` (or
   the whole ``~/lib/python`` entry from ``PYTHONPATH``, if nothing else
   still needs it) by hand -- ``pip uninstall`` has no knowledge of it,
   since ``make pyinstall`` never went through pip in the first place.


conda-forge (prebuilt, includes the CLI binary)
*************************************************

For **Linux (x86_64)** and **macOS (both Intel and Apple Silicon)**::

    conda install -c conda-forge pamtra

or, with `mamba <https://mamba.readthedocs.io/>`_::

    mamba install -c conda-forge pamtra

This is a self-contained install like the pip wheel above -- no system
libraries, compiler, or manually-managed environment needed -- but with one
important difference: it **does** include the standalone ``pamtra`` CLI
binary (:ref:`pamtra`). The pip wheel omits it (see the warning above)
specifically because bundling netCDF-Fortran into a wheel is impractical;
the conda-forge package doesn't have that constraint, since conda already
manages netCDF-Fortran (and every other dependency) as a regular package,
so the recipe just declares it as a dependency and builds the CLI normally
alongside the ``pyPamtra`` Python extension. In short: pip gives you
``pyPamtra`` only; conda-forge gives you ``pyPamtra`` *and* the ``pamtra``
binary.

Not available for **Windows** -- the conda-forge recipe currently skips it,
since PAMTRA's Fortran code has never been built against MSVC. Use WSL2
below instead.

Get the code
*************

The version control system git (http://git-scm.com/) is used to keep track
of the code. Get a copy of the model with::

    git clone https://github.com/igmk/pamtra.git
    cd pamtra


conda-forge / pixi, building from source (dev/editable installs)
**********************************************************************

If you need an editable/dev install (see below), or want to build on a
platform the prebuilt conda-forge package doesn't cover, its own
dependencies (openblas, fftw, netcdf, a matching C/Fortran compiler pair)
are all available as standalone `conda-forge <https://conda-forge.org/>`_
packages too, including working macOS (both Intel and Apple Silicon)
builds. This avoids needing a system package manager (apt/brew) at all, and
is the same install path used by PAMTRA's CI.

With `conda <https://docs.conda.io/>`_ or `mamba <https://mamba.readthedocs.io/>`_::

    conda create -n pamtra -c conda-forge python numpy scipy netcdf4 matplotlib \
        meson meson-python cython pkg-config fftw libopenblas libnetcdf \
        c-compiler fortran-compiler
    conda activate pamtra
    pip install --no-deps --no-build-isolation .

Or, if you use `pixi <https://pixi.sh/>`_ (also conda-forge-based), the
repository already has a ``pixi.toml`` with this dependency set defined ---
just run::

    pixi install
    pixi run install
    pixi run test    # optional, runs the test suite

``pixi run install`` and the manual ``pip install`` above both build with
``--no-build-isolation``, so the compiler/library versions actually pinned
in your conda/pixi environment are used instead of a fresh isolated build
environment.


Linux (Ubuntu), apt
********************

Install the system libraries needed to compile PAMTRA::

    sudo apt install git gfortran libopenblas-dev libfftw3-dev libnetcdff-dev

Create and activate a virtual environment::

    sudo apt install python3-venv
    python3 -m venv pamtraenv
    source pamtraenv/bin/activate

Install the Python build and runtime dependencies::

    pip install numpy scipy matplotlib netcdf4 xarray meson numexpr cython

Then install PAMTRA itself::

    pip install .

.. warning::
    On some Linux systems, OpenBLAS is not thread-safe when run with
    multiple parallel jobs. If you see hangs or crashes, set::

        export OPENBLAS_NUM_THREADS=1

    before starting python.


macOS, Homebrew
*****************

Install the required libraries with `Homebrew <https://brew.sh>`_::

    brew install openblas pkgconf netcdf fftw

Then install PAMTRA, pointing the C compiler at the Homebrew ``gcc`` that
matches your ``gfortran`` (adjust the version number, e.g. ``gcc-14``, to
whatever ``brew install gcc`` provides on your system)::

    env CC=gcc-14 pip install .

.. note::
   ``openblas`` is keg-only in Homebrew (macOS ships BLAS/LAPACK via the
   Accelerate framework instead), so its ``.pc`` file is not on the default
   ``pkg-config`` search path. The build automatically falls back to
   ``brew --prefix openblas`` to locate it, so you do **not** need to
   manually export ``PKG_CONFIG_PATH`` for openblas.


Windows, WSL2
**************

On Windows, install `WSL2 <https://learn.microsoft.com/windows/wsl/install>`_
with an Ubuntu distribution, then follow the Linux instructions above
verbatim inside the WSL2 Ubuntu shell -- there is no separate native Windows
build.


DKRZ Levante HPC
*****************

::

    module load git
    spack load /fwvsvi # python3.9.9
    python -m venv pamtraenv
    source pamtraenv/bin/activate
    pip install numpy scipy matplotlib netcdf4 cython xarray meson

    git clone https://github.com/igmk/pamtra.git
    cd pamtra

    spack load /bcn7mbu # gcc 11.2
    spack load /tpmfvwu # openblas 0.3.18 gcc 11.2
    spack load /fnfhvr6 # fftw 3.10.10
    spack load /jn6xcuy # netcdf-fortran 4.6.1 gcc 11.2
    pip install .

The exact spack hashes may change over time; use ``spack find`` to look up
the current ones if a ``spack load`` fails. See also
``install_levante_readmefirst.sh`` in the repository root, which automates
an equivalent module/env setup using the legacy ``Makefile`` build --
kept around specifically for HPC deployments like this one, where
hand-tuned linker flags against cluster module paths are simpler to
express as Makefile variables than through meson/pkg-config.

For Jupyter support::

    pip install ipykernel
    python -m ipykernel install --user --name=pamtra-kernel --display-name="pamtra kernel"


Editable / development install
*******************************

While developing PAMTRA, an editable install avoids a full reinstall after
every Fortran change::

    pip install --no-build-isolation -e . -Cbuild-dir=build


Download data
*************

This data includes the land surface emissivity maps and some scattering
databases. Many features (e.g. Mie-sphere scattering, the built-in surface
emissivity defaults) work without it.

For ``pyPamtra``, nothing to do here: if ``PAMTRA_DATADIR`` isn't set at all,
``import pyPamtra`` downloads and caches the data automatically (via
`pooch <https://www.fatiando.org/pooch/>`_, with a checksum check), the first
time only. If you'd rather not have the first import trigger a ~250 MB
download, fetch it ahead of time the same way::

  export PAMTRA_DATADIR=$(pamtra-fetch-data)

and add that line to your shell startup file. To explicitly skip the data
entirely (rather than let it auto-download), set ``PAMTRA_DATADIR=""``
before importing.

For the standalone ``pamtra`` binary (:ref:`pamtra`, no pip/Python dependency
at runtime, so no auto-download either), download and unpack the data
manually::

  wget -q -O data.tar.bz2 https://github.com/igmk/pamtra/releases/download/data-v1/pamtra_data.tar.bz2
  tar xjf data.tar.bz2
  rm data.tar.bz2
  echo 'export PAMTRA_DATADIR="wherever/it/is/"' >> ~/.bashrc
  source ~/.bashrc


Start PAMTRA
************

You can start using pyPamtra in python with ::

  import pyPamtra


Build documentation
********************

The documentation is built using Sphinx. Install the build dependency
with pip::

    pip install sphinx

Then build it using the Makefile in the ``doc`` directory::

  cd doc
  make html

..  _installation:


Installation
============

PAMTRA is installed with pip, which compiles the Fortran core and builds the
``pyPamtra`` Python extension in one step via `meson-python
<https://mesonbuild.com/meson-python/>`_. This has replaced the old
``make`` / ``make pyinstall`` workflow.

.. note::
   The standalone Fortran command line binary (:ref:`pamtra`) is still built
   with the legacy ``Makefile`` (``make pamtra``), since it links directly
   against the NetCDF Fortran interface via ``nf-config``. The pip install
   only builds the Python package (:ref:`pyPamtra`).


Get the code
*************

The version control system git (http://git-scm.com/) is used to keep track
of the code. Get a copy of the model with::

    git clone https://github.com/igmk/pamtra.git
    cd pamtra


Linux (Ubuntu)
**************

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


macOS
*****

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
the equivalent setup for the legacy ``Makefile`` build.

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

Although PAMTRA can be used without additional data by simply setting the
data path to an empty directory::

  export PAMTRA_DATADIR=""

it is recommended that you download the data. The data includes the land
surface emissivity maps and some scattering databases. They can be found on
the servers of University of Cologne

  https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS

Download and unpack the data::

  wget -q -O data.tar.bz2 https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS/download
  tar xjf data.tar.bz2
  rm data.tar.bz2

and set the ``$PAMTRA_DATADIR`` variable, e.g. by adding it to your shell
startup file::

  echo 'export PAMTRA_DATADIR="wherever/it/is/"' >> ~/.bashrc
  source ~/.bashrc


Start PAMTRA
************

You can start using pyPamtra in python with ::

  import pyPamtra


Build documentation
********************

The documentation is built using Sphinx. Install the build dependencies
with pip::

    pip install sphinx numpydoc sphinx_rtd_theme

Then build it using the Makefile in the ``doc`` directory::

  cd doc
  make html

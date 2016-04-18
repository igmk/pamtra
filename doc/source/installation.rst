..  _installation:


Installation
============


Dependencies on Ubuntu
**********************

On Ubuntu, install the following packages in advance in order to compile the Fortran Pamtra version::

    sudo apt-get install gcc gfortran libnetcdff5 netcdf-bin libnetcdf-dev liblapack-dev libzip-dev git gitk

The model is tested with gcc version 4.8.2.

To get the Python version, the following packages are additionally required::

    sudo apt-get install python-pandas python-numpy python-matplotlib python-scipy python-netcdf python-pip

The Netcdf4 python package is not in the repositories, install with::

    sudo pip install netcdf4


Dependencies on Mac OS X
************************

On Mac OS X, it is recommended to use brew (http://brew.sh) to install gfortran (via gcc) and netcdf ::

    brew install gcc
    brew install netcdf --enable-fortran

For the Python version, it is recommended not to use OS X's default python version,
but to install an independent one, e.g. with brew or conda
(https://www.continuum.io/downloads). In addition, the following packages are required::

    pip install pandas numpy scipy matplotlib netcdf4

Please note that netcdf4 must be installed using pip even if you use the conda
package manager. The reason is that conda brings its own netcdf library, but without
the fortran libraries which are required by the fortran part of PAMTRA.

Get model from git repository
*****************************
The version control system git (http://git-scm.com/) is used to keep track of the code. Get a copy of the model with::

    git clone --recursive https://github.com/igmk/pamtra

The very basics of git can be found here https://try.github.io/levels/1/challenges/1 .
"--recursive" is required because gut submodules are used.

Build Pamtra
*******************
Simply type ::

  make

to build :ref:`pamtra` and :ref:`pyPamtra`. You can build them also separately with ::

  make pamtra

and ::

  make py

Usually superuser permission are required to install python routines. To avoid
that a local python library folder is used in ~/lib/python/ and this path has to
be added to the $PYTHONPATH variable of your shell (assuming you are using Ubuntu
and bash) ::

  echo 'export PYTHONPATH=$PYTHONPATH:$HOME/lib/python' >> ~/.bashrc

For Mac OS X, do ::

    echo 'export PYTHONPATH=$PYTHONPATH:$HOME/lib/python' >> ~/.bash_profile


Then, the python routines can be installed with ::

  make pyinstall

To start using pyPamtra, you have to open a new bash session or source the ~/.bashrc ::

  source ~/.bashrc

You can start using pyPamtra in python with ::

  import pyPamtra

Build documentation
*******************

Several package have to be installed to be able to build the documentation. The documentation is build using sphinx ::

    sudo apt-get install python-sphinx

In addition, the numpydoc is required ::

    sudo apt-get install python-numpydoc

If not available try ::

    sudo easy_install numpydoc

In addition, the sphinx-fortran-extension is required which can be found in the tools folder of Pamtra::

    cd tools/sphinx-fortran-extension
    sudo python setup.py install

if you do not have root permissions you can also use instead of the last line::

    python setup.py install --user

Eventually, you can build the documentation by using the Makefile in the pamtra main directory with ::

  make htmldoc

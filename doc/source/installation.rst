..  _installation:


Installation
============


Dependencies
************

On Ubuntu, install the following packages in advance in order to compile the Fortran Pamtra version::

    sudo apt-get install gcc gfortran libnetcdff5 netcdf-bin libnetcdf-dev liblapack-dev libzip-dev git gitk

The model is tested with gcc version 4.8.2.

To get the Python version, the following packages are additionally required::

    sudo apt-get install python-pandas python-numpy python-matplotlib python-scipy python-netcdf python-pip

The Netcdf4 python package is not in the repositories, install with::

    sudo pip install netcdf4

Get model from git repository
*****************************
The version control system git (http://git-scm.com/) is used to keep track of the code. Get a copy of the model with::

    git clone ssh://your_user_name@gop/srv/git_repos/pamtra.git/

The very basics of git can be found here https://try.github.io/levels/1/challenges/1 

Build Pamtra
*******************
Simply type ::

  make

to build :ref:`pamtra` and :ref:`pyPamtra`. You can build them also separately with ::

  make pamtra

and ::

  make py

Usually superuser permission are required to install python routines. To avoid that, a local python library folder has to be created ::
  
  mkdir -p ~/lib/python/

and this path has to be added to the $PYTHONPATH variable of your shell (assuming you are using Ubuntu and bash) ::

  echo 'export PYTHONPATH=$PYTHONPATH:$HOME/lib/python' >> ~/.bashrc

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

In addition, the sphinx-fortran-extension is required ::

    cd /tmp
    git clone https://github.com/maahn/sphinx-fortran-extension
    cd /tmp/sphinx-fortran-extension
    sudo python setup.py install
    
if you do not have root permissions you can also use instead of the last line::

    python setup.py install --user

Eventually, you can build the documentation by using the Makefile in the pamtra main directory with ::

  make htmldoc

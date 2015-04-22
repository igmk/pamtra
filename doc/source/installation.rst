
Installation
============


Dependencies
************

On Ubuntu, install the following packages in advance in order to copmpile the Fortran Pamtra version::

    sudo apt-get install gcc gfortran libnetcdff5 netcdf-bin libnetcdf-dev liblapack-dev git gitk

The model is tested with gcc version 4.8.2.

To get the Python version, the following packages are additionally required::

    sudo apt-get install python-pandas python-numpy python-matplotlib python-scipy python-netcdf python-pip

The Netcdf4 python package is not in the repositiries, install with::

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

to build pamtra and pyPamtra. You can build them also sperately with ::

  make pamtra

and ::

  make py

The python routines have to be installed with::

  make pyinstall

.. note:: Add python path etc

Build documentation
*******************

Several package have to be installed to be able to build the documentation. The documentation is build using sphinx::

    sudo apt-get install python-sphinx

In addition, the numpydoc is required::

    sudo apt-get install python-numpydoc

If not available try::

    sudo easy_install numpydoc

In addition, the sphinx-fortran-extension is required::

    cd /tmp
    git clone https://github.com/paulromano/sphinx-fortran-extension
    cd /tmp/sphinx-fortran-extension
    sudo python setup.py install

Eventually, you can build the documentation with::

  make htmldoc

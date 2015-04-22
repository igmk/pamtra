
Installation
============


Dependencies
************

liblapack, netcdf fortran headers

python: f2py, numpy, matplotlib, pandas(?), paramiko (optional), netcdf4

Build documentation
*******************

Several package have to be installed to be able to build the documentation. The documentation is build using sphinx::

    sudo apt-get install python-sphinx

In addition, the numpydoc is required::

    sudo apt-get python-numpydoc

If not available try::

    sudo easy_install numpydoc

In addition, the sphinx-fortran-extension is required::

    cd /tmp
    git clone https://github.com/paulromano/sphinx-fortran-extension
    cd /tmp/sphinx-fortran-extension
    sudo python setup.py install

..  _installation:


Installation
============


Virtual Machine
***************

For testing, a virtual machine is provided. You have to install the following software in order to 
run the virtual machine with PAMTRA on your computer:

* Virtualbox https://www.virtualbox.org/
* Vagrant https://www.vagrantup.com/

Change the directory to tools/v in the terminal (on linux and 
mac os `cd tools/virtual_machine`) and type 

    vagrant up

This command will configure the virtual machine based on the information provided
in `Vagrantfile`. It will download a lot of data, so don't do that with a bad 
internet connection and make sure you have some GB space on your hard disk.

When the machine is running, you can connect to it with 

    vagrant ssh

in order to have terminal access. When you are done, do 

    vagrant halt

from the host computer to turn the virtual machine off. You can start with 

    vagrant up

again. To delete the virtual machine and its containing data, do 

    vagrant destroy

from the host. Note that the `shared` folder in the virtual machine's home directory links to 
the `virtual_machine` of the host system. I.e. you can use that directory for file exchange and for editing files with you favorite editor. 

When the virtual machine is running, you should be able to access a jupyter notebook
at http://localhost:8880 from your host machine.  The password is `jupyter`. All examples in shared/examples
should work and produce a plot. Jupyter notebook is a python environment, but can be 
used with other languages, too. More information at http://jupyter.org. The default python version is 3.6, but 2.7 is installed which needs to be used with PAMTRA. From the command line, you can switch with 

    conda activate py27

to python 2.7, and go back to 3.6 with 

    conda activate base

If you miss any of your favorite python packages, do

    conda install missing_module_name

Make sure that the correct python version is selected with `conda activate `. 

The `update.sh` script in the shared folder can be used to update the PAMTRA model
and the data provided on the FTP server if they get updated.


Dependencies on Ubuntu
**********************

On fresh Ubuntu 16.04 installation, the following packages need to be installed to get PAMTRA from github repository and to compile and install PAMTRA::

    sudo apt install git gfortran libnetcdf-dev libnetcdff-dev liblapack-dev libfftw3-dev python-dev python-numpy 

The model is tested with gcc version 4.8.2.

Not all numpy versions work with PAMTRA, 1.11.3 and 1.12.1 are confirmed to work. 1.13.3 and 1.14.3 do not work. You can check the numpy version in python ::

    import numpy
    print numpy.__version__

To ensure to have the correct version installed, do::

    sudo pip install numpy==1.12.1

if you do not have root permissions you can also use instead of the last line::

    pip install --user numpy==1.12.1

Although not required for comppilation and installation, to use PAMTRA, some additional python packages need to be installed on your system or python environemnt.::

    sudo apt install python-matplotlib python-pandas python-scipy python-netcdf

More recent Ubuntu versions have `python-netcdf4` instead of `python-netcdf`. For older versions use pip ::

    sudo apt install python-pip
    sudo pip install netcdf4==1.3.1 # newer version have problems with `import netCDF3`

If the system's numpy version is to recent, use pip ::

    pip install --user netcdf4

.. warning::
    As of June 2018, do NOT use conda for Ubuntu because the provided libgfortran 
    library does not work with PAMTRA.


Dependencies on Mac OS X
************************

On Mac OS X, it is recommended to use brew (http://brew.sh) to install gfortran (via gcc) and netcdf ::

    brew install gcc
    brew install fftw
    brew install netcdf --enable-fortran

For the Python version, it is recommended not to use OS X's default python version,
but to install an independent one, e.g. with brew or conda
(https://www.continuum.io/downloads). Note that pyPamtra does not support Python3 yet.
In addition, the following packages are required::

    pip install pandas numpy==1.12.1 scipy matplotlib netcdf4

Please note that netcdf4 must be installed using pip even if you use the conda
package manager. The reason is that conda brings its own netcdf library, but without
the fortran libraries which are required by the fortran part of PAMTRA. Similar to Ubuntu teh most recent numpy versions do not work with PAMTRA. 

Installation on Microsoft Windows 10 with windows subsystem for linux
**********************************************************************
To install windows subsystem for linux follow the instructions on ::

https://docs.microsoft.com/de-de/windows/wsl/install-win10

Install ubuntu 16.04 from the Microsoft Store. After configuration, you need to install additional packages within the ubuntu linux system ::

  sudo apt update
  sudo apt install git make

Afterwards, follow the instructions for Ubuntu.


Get model from git repository
*****************************
The version control system git (http://git-scm.com/) is used to keep track of the code. Get a copy of the model with::

    git clone --recursive https://github.com/igmk/pamtra

The very basics of git can be found here https://try.github.io/levels/1/challenges/1 .
"--recursive" is required because git submodules are used.


Build PAMTRA
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

Start PAMTRA
******************
Although PAMTRA can be used without additional data by simply setting the data path to an empty directory, 

  echo 'export PAMTRA_DATADIR=""'' >> ~/.bashrc

it is recommended that you download the data. The data includes the land surface emissivity maps and some scattering databases. They can be found on the servers of University of Cologne

  https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS

Unpack the data and set the $PAMTRA_DATADIR variables ::

  echo 'export PAMTRA_DATADIR="wherever/it/is/"'' >> ~/.bashrc

To start using pyPamtra, you have to open a new bash session or source the ~/.bashrc ::

  source ~/.bashrc

You can start using pyPamtra in python with ::

  import pyPamtra

Build documentation
*******************

Several package have to be installed to be able to build the documentation. The documentation is build using sphinx ::

    sudo apt install python-sphinx

In addition, the numpydoc is required ::

    sudo apt install python-numpydoc

If not available try ::

    sudo easy_install numpydoc

In addition, the sphinx-fortran-extension is required which can be found in the tools folder of PAMTRA::

    cd tools/sphinx-fortran-extension
    sudo python setup.py install

if you do not have root permissions you can also use instead of the last line::

    python setup.py install --user

Eventually, you can build the documentation by using the Makefile in the PAMTRA main directory with ::

  make htmldoc

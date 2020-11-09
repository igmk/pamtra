..  _installation:


Installation
============

.. note::
   Note that PAMTRA finally supports Python 2.7 *and* Python 3.6 - 3.8 as of August 2020.


Quickstart: Virtual Machine
***************************

This is the recommended way to test the model. A virtual machine is provided 
which makes the installation 100% automatic. You have to install the following
 software in order to run the virtual machine with PAMTRA on your computer:

* Virtualbox https://www.virtualbox.org/
* Vagrant https://www.vagrantup.com/

When using Mac OS or Windows, download the installation routines frotm the
websites, on Linux use ::

    sudo apt-get install virtualbox vagrant

Change the directory to ´tools/virtual_machine in the terminal (on linux and 
mac os `cd tools/virtual_machine`) and type ::

    vagrant up

This command will configure the virtual machine based on the information provided
in `Vagrantfile`. It will download a lot of data, so don't do that with a bad 
internet connection and make sure you have some GB space on your hard disk.

When the machine is running, you can connect to it with ::

    vagrant ssh

in order to have terminal access. When you are done, do ::

    vagrant halt

from the host computer to turn the virtual machine off. You can start with ::

    vagrant up

again. To delete the virtual machine and its containing data, do ::

    vagrant destroy

from the host. Note that the `shared` folder in the virtual machine's home directory links to 
the `virtual_machine` of the host system. I.e. you can use that directory for file exchange and for editing files with you favorite editor. 

When the virtual machine is running, you should be able to access a jupyter notebook
at http://localhost:8880 from your host machine.  The password is `jupyter`. All examples in shared/examples
should work and produce a plot. Jupyter notebook is a python environment, but can be 
used with other languages, too. More information at http://jupyter.org. The default python version is 3.8, but 2.7 is installed as well, PAMTRA is installed for both environments. From the command line, you can switch with ::

    conda activate py27

to python 2.7, and go back to 3.8 with ::

    conda activate base

If you miss any of your favorite python packages, do ::

    pip install missing_module_name

Make sure that the correct python version is selected with `conda activate `. When using conda instaed of pip, make sure that libgfortran is not installed by conda. For some reason, PAMTRA does not like conda's libgfortran on Linux.

The `update.sh` script in the shared folder can be used to update the PAMTRA model
and the data provided on the FTP server if they get updated.




Install Dependencies
********************


Ubuntu 20.04
------------

On a fresh Ubuntu 20.04 installation, the following packages need to be installed to get PAMTRA from the github repository and to compile and install PAMTRA::

    sudo apt install git gfortran libnetcdf-dev libnetcdff-dev liblapack-dev libfftw3-dev python3-dev python3-numpy 

Replace python3-X with python-X when using Python 2.7 instead of 3.

Although not required for compilation and installation, to use PAMTRA, some additional python packages need to be installed on your system or python environemnt.::

    sudo apt install python3-matplotlib python3-pandas python3-scipy python3-netcdf4

.. warning::
    As of August 2020, do NOT use conda for Ubuntu because the provided libgfortran 
    library does not work with PAMTRA.

Now, follow -  :ref:`Get model from git repository`


Mac OS X
--------

On Mac OS X, it is recommended to use brew (http://brew.sh) to install gfortran (via gcc) and netcdf ::

    brew install gcc
    brew install fftw
    brew install netcdf

For the Python version, it is recommended not to use OS X's default python version,
but to install an independent one, e.g. with brew or conda
(https://www.continuum.io/downloads). 
In addition, the following packages are required::

    pip install pandas numpy scipy matplotlib netcdf4

or ::

    conda install pandas numpy scipy matplotlib netcdf4

Now, follow -  :ref:`Get model from git repository`


Microsoft Windows 10 with windows subsystem for linux
-----------------------------------------------------
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

Download data
*************
Although PAMTRA can be used without additional data by simply setting the data path to an empty directory, :: 

  echo 'export PAMTRA_DATADIR=""' >> ~/.bashrc

it is recommended that you download the data. The data includes the land surface emissivity maps and some scattering databases. They can be found on the servers of University of Cologne

  https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS

Download and unpack the data::

  wget -q -O data.tar.bz2 https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS/download
  tar xjf data.tar.bz2
  rm data.tar.bz2

and set the $PAMTRA_DATADIR variables ::

  echo 'export PAMTRA_DATADIR="wherever/it/is/"' >> ~/.bashrc

To start using pyPamtra, you have to open a new bash session or source the ~/.bashrc ::

  source ~/.bashrc


Start PAMTRA
************

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

Eventually, you can build the documentation by using the Makefile in the PAMTRA main directory with ::

  make htmldoc

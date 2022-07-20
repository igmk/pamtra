#!/bin/bash

# Hi,
# this script performs a series of commands that have successfully installed pamtra on Levante (DKRZ HLRE-4 supercomputing machine)

# It also modifies your user .bashrc and .bash_profile in order to always make the required modules available for using pamtra already in the login nodes (not recommended) or on the DKRZ jupyterhub
# if you do not want to modify your profile, just comment out those lines

# TOCHECK: those modules will be needed to be loaded by the compute nodes. 
# If you submit a batch script that runs pamtra remember to set the required environment variables and load the modules and packages

# So far there is only the default compiling option available.
# This uses gcc as default compiler, and openblas for the implementation of blas and lapack (beware of performance issues with openblas which is worked around by setting environment variable OPENBLAS_NUM_THREAD=1)

# to install pamtra on your $HOME just run this script with
# bash install_levante_readmefirst.sh

# enjoy

#################################################################################
# First, load the required modules.
# This block might need some changes if the available libraries vary
spack load gcc@11.2.0%gcc@11.2.0
module load python3/2022.01-gcc-11.2.0
spack load /fnfhvr6 # one of the 2 copies of fftw@3.3.10%gcc available on Levante
spack load openblas%gcc
module load netcdf-fortran/4.5.3-gcc-11.2.0
#################################################################################

# Create a PAMTRA_DATADIR, change the PAMTRA_DATA_FOLDER variable if you want it to point somewhere else
THISDIR=$(pwd)
echo $THISDIR
PAMTRA_DATA_FOLDER=$THISDIR/pamtra_data/
mkdir -p $PAMTRA_DATA_FOLDER

# If you do not need pamtra_data on scattering and emissivity, or the example data, you can comment out the relevant lines from this block
wget -q -O data.tar.bz2 https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS/download
wget -q -O example_data.tar.bz2 https://uni-koeln.sciebo.de/s/28700CuFssmin8q/download
tar -xjvf example_data.tar.bz2 -C $PAMTRA_DATA_FOLDER
tar -xjvf data.tar.bz2 -C $PAMTRA_DATA_FOLDER
rm example_data.tar.bz2
rm data.tar.bz2


# This is needed
echo "export PYTHONPATH=$PYTHONPATH:$HOME/lib/python" >> ${HOME}/.bashrc
echo "export PAMTRA_DATADIR=${PAMTRA_DATA_FOLDER}" >> ${HOME}/.bashrc
echo "export OPENBLAS_NUM_THREADS=1" >> ${HOME}/.bashrc # this is important only for the openblas option
echo "source .bashrc" >> ${HOME}/.bash_profile

# Now, let's modify the Makefile for openblas
sed 's/-llapack//g' Makefile > Makefile.levante
sed -i 's/-lblas/-lopenblas/g' Makefile.levante

# and install
make -f Makefile.levante
make pyinstall -f Makefile.levante

# Make pamtra readily available, even on the login node
source $HOME/.bashrc

# Bye
echo "If no error message was printed, you should be good to go! Thank you for using pamtra"
echo "Remember to load the python module if you want to use pamtra in a batch script!"



# #!/bin/bash

mkdir -p ~/software
mkdir -p /vagrant/data
mkdir -p ~/data

ln -s /vagrant ~/shared
chmod +x ~/shared/update.sh


# # Sytem tools!
# Ubuntu 20.04
# sudo apt-get update -y 
# sudo apt-get install -y -qq tree git gcc gfortran liblapack-dev libfftw3-dev  
# sudo apt-get install -y -qq libhdf5-dev make libnetcdff7 libnetcdf-dev 
# sudo apt-get install -y -qq netcdf-bin gsl-bin libgsl0-dev libgsl23 libgslcblas0 lftp

# ubtuntu 16.04
sudo apt-get update -y 
sudo apt-get install -y -qq tree git gcc gfortran liblapack-dev libfftw3-dev  
sudo apt-get install -y -qq libnetcdff5 libnetcdf-dev netcdf-bin gsl-bin 
sudo apt-get install -y -qq libgsl0-dev libgsl0ldbl lftp
sudo apt-get install -y -qq libhdf5-dev

# guest additions update
# cd /opt 
# sudo wget -c -q http://download.virtualbox.org/virtualbox/5.2.12/VBoxGuestAdditions_5.2.12.iso \
#                        -O VBoxGuestAdditions_5.2.12.iso
# sudo mount VBoxGuestAdditions_5.2.12.iso -o loop /mnt
# cd /mnt
# sudo sh VBoxLinuxAdditions.run --nox11 -- --force
# cd /opt
# sudo rm *.iso
# sudo /etc/init.d/vboxadd setup
# sudo chkconfig --add vboxadd
# sudo chkconfig vboxadd on
# cd ~

wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda

echo ". /home/vagrant/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
echo "conda activate base">> ~/.bashrc
. /home/vagrant/miniconda/etc/profile.d/conda.sh
conda activate base

echo "export PYTHONPATH=$PYTHONPATH:$HOME/lib/python" >> ~/.bashrc
export PYTHONPATH=$PYTHONPATH:$HOME/lib/python

# download PAMTRA data
cd ~/data
wget -q -O data.tar.bz2 https://uni-koeln.sciebo.de/s/As5fqDdPCOx4JbS/download
tar xjf data.tar.bz2
rm -rf data.tar.bz2
cd 

echo "export PAMTRA_DATADIR=$HOME/data" >> ~/.bashrc
export PAMTRA_DATADIR=$HOME/data

# # python 3.x
conda activate base
conda config --add channels conda-forge 
# Pamtra1 hates apparently conda's libgfortran... so use pip!
pip install  --no-cache-dir -q cython numpy scipy xarray dask numba jupyter matplotlib ipython pytest netcdf4 
# pip install -q arm_pyart
conda install jupyter 
#launch jupyter on startup
mkdir -p /home/vagrant/.jupyter
(crontab -l ; echo "@reboot cd /home/vagrant; source ~/.bashrc;  /home/vagrant/miniconda/bin/jupyter notebook  --ip 0.0.0.0 --no-browser >> jupyter.log")| crontab -
echo '{
  "NotebookApp": {
    "password": "sha1:b269d96dff97:82c6b6c3ba8438a203a2d635b0c6d53b38f51ae9"
  }
}' >> /home/vagrant/.jupyter/jupyter_notebook_config.json

pip install  --no-cache-dir -q jupyter_contrib_nbextensions
jupyter contrib nbextension install --user

nohup /home/vagrant/miniconda/bin/jupyter notebook  --ip 0.0.0.0 --no-browser >> jupyter.log 2>&1 &


#now python 2.7
conda create -y --name py27 python=2.7
conda activate py27
conda config --add channels conda-forge 
# Pamtra1 hates apparently conda's libgfortran... so use pip!
# Python 2.7 is not maintained any more, so let's specify some versions, just to be sure
pip install  --no-cache-dir -q numpy==1.13.3 scipy==1.2.3 llvmlite==0.31.0 numba==0.47.0 ipykernel==4.10.1 matplotlib==2.2.5 ipython==5.10.0 xarray==0.11.3  netcdf4==1.3.1
python -m ipykernel install --user


# clean up
conda clean -a -y


# pamtra 
cd ~/software
git clone --recursive https://github.com/igmk/pamtra
cd ~/software/pamtra


#for python 3.x
cd ~/software/pamtra
conda activate
make clean
make pyinstall
cd ~/software
python -c "import matplotlib; matplotlib.use('Agg'); exec(open('pamtra/examples/pyPamTest_ssrg_hogan2014.py').read());plt.savefig('test_pamtra_3.png')"


# for python 2.7
cd ~/software/pamtra
conda activate py27
make clean
make pyinstall
cd ~/software
python -c "import matplotlib; matplotlib.use('Agg'); execfile ('pamtra/examples/pyPamTest_ssrg_hogan2014.py');plt.savefig('test_pamtra_27.png')"




# # pamtra2
# conda activate base
# cd ~/software
# git clone --recursive https://github.com/maahn/pamtra2
# cd ~/software/pamtra2
# rm -rf build
# python setup.py install
# pytest -q

# #libradtran
# conda activate py27
# cd ~/software
# wget -q http://www.libradtran.org/download/libRadtran-2.0.2.tar.gz
# tar -xzf libradtran-2.0.2.tar.gz
# cd ~/software/libradtran-2.0.2
# #to make up for the hard link in the tar archive
# # cp src_py/uvspec_lex.l src/uvspec_lex.l
# export F77=gfortran
# make clean
# ./configure
# make

#libradtran
# conda activate py27
# cd ~/software
# wget -q http://www.libradtran.org/download/history/libRadtran-1.7.tar.gz
# # sometime it is compressed, sometimes not... strange! So try both
# tar -xf libRadtran-1.7.tar.gz
# tar -xzf libRadtran-1.7.tar.gz
# # cp libRadtran-1.7/src/uvspec_lex.l libRadtran-1.7/python/uvspec_lex.l
# cd ~/software/libRadtran-1.7
# #apply patch to allow more heights
# patch -p1 < /vagrant/libradtran.patch 
# # export F77=gfortran-4.7
# ./configure
# make    
# make check
# sudo make install
# echo "export LIBRADTRAN_DATA_FILES=/usr/local/share/libRadtran/data/"  >> ~/.bashrc
# export LIBRADTRAN_DATA_FILES=/usr/local/share/libRadtran/data/

# echo 'DOWNLOAD DATA... THIS TAKES VERY LONG!'
# echo 'check the arm-summer-2018/group3/virtual_machine/data folder to see the progress!'
# #download DATA
# cd ~/shared/data
# lftp -v -c 'open ftp://ftp.cdc.noaa.gov/Public/mmaahn/arm-summer-school/ && mirror && exit'




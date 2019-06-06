# #!/bin/bash
set -e

echo '############################################################'
echo updating pamtra
echo '############################################################'

# pamtra
source activate py27
cd ~/software/pamtra
git pull
make clean
make pyinstall

# echo '############################################################'
# echo pamtra2
# echo '############################################################'

# # pamtra2
# source activate base
# cd ~/software/pamtra2
# git pull
# rm -rf build
# python setup.py install



# #!/bin/bash
set -e

echo '############################################################'
echo updating pamtra for Python 2.7
echo '############################################################'

# pamtra
source activate py27
cd ~/software/pamtra
git pull
make clean
make pyinstall

# #!/bin/bash
set -e

echo '############################################################'
echo updating pamtra for Python 3
echo '############################################################'

# pamtra
source activate
cd ~/software/pamtra
git pull
make clean
make pyinstall
#!/bin/bash -l
#PBS -q default
#PBS -l nodes=1:ppn=8
#PBS -l mem=100mb
#PBS -l walltime=05:00:00
#PBS -A AG-Crewell
#PBS -t 0-3



module load gnu/4.4.4
module load python 

export NCORES=`cat $PBS_NODEFILE | wc -l`
cd $PBS_O_WORKDIR

echo "python bin/ppserver.py -t 30 -s 'pyPamtra' -w $NCORES -a"

python bin/ppserver.py -t 30 -s 'pyPamtra' -w $NCORES -a 


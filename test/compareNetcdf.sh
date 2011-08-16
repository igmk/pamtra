#!/usr/bin/env bash
#wrapper for diff to compare netcdf pamtra output with reference output.

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

status=0

for fname in $(ls tmp)
   do
   fname=$(basename tmp/$fname .txt)
   python dumpNetcdf.py tmp/$fname > tmp/$fname.txt
   status=$[$?+$status]
   diff tmp/$fname.txt referenceOutput/$1/$fname.txt
   status=$[$?+$status]
done


exit $status
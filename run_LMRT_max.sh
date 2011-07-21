#!/bin/bash

cd dirname $0

outdir=$2
# file=$1
for file in ls $1*gz
do
	tmpFile=${file%.*}

	echo "decompress $file > $tmpFile..."
	zcat $file > $tmpFile
	echo "decompress $file > $tmpFile ... done"

		for frequency in 10.65 # 21.0 36.5
		do
			echo "starting ./pamtra $(basename $tmpFile) $frequency..."
			./pamtra $(basename $tmpFile) $frequency &
			echo "starting ./pamtra $(basename $tmpFile) $frequency...done"
		wait

		for g in $outdir/*nc
			do
			echo "compress profile $g..."
			gzip $g &
		done
		wait

		done

	echo "deleting $tmpFile..."
	rm $tmpFile
	echo "deleting $tmpFile... done"
done
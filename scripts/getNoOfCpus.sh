#!/bin/bash
saveCPUs=0.8

noP=$(grep processor /proc/cpuinfo | wc|awk '{print $1}')
load=$(uptime| awk -F, '{print $(NF-1)}')
useP=$(echo $noP\-$saveCPUs\-\($load\%$noP\) |bc | cut -d '.' -f1)

if [[ $useP =~ ^-?[0-9]+$ ]]
then
	:
else
	useP=0
fi
if  [ $useP -gt 0 ]
then
	echo $useP
fi
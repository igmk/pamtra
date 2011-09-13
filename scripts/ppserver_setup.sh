#!/bin/bash

hosts="albe   euros   refoli   notos   roumet   helm   irifi trombe   habagat   embat caspar"

timeout=300

#how many cpus shall be saved?
saveCPUs=1

totalWorkers=0

OkHosts=()

if [ $1 == "create" ]
then

	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" ${host} "echo 2>&1" 
		then
			noP=`ssh -o "BatchMode yes" $host grep processor /proc/cpuinfo | wc|awk '{print $1}'`
			load=`ssh -o "BatchMode yes" $host uptime| awk -F, '{print $(NF)}'`
			useP=`echo $noP\-\($load\%$noP\) |bc | cut -d '.' -f1`
			if [[ $useP =~ ^-?[0-9]+$ ]]
			then
				:
			else
				useP=0
			fi
			if  [ $useP -gt $saveCPUs ]
			then
				useP=`echo $useP -$saveCPUs|bc`
				echo "$host: nohup ppserver.py -t $timeout -w $useP &"
				ssh -o "BatchMode yes" $host ppserver.py -t $timeout -w $useP &
				if [ $? -eq 0 ]
				then
					OkHosts="\"$host\",$OkHosts" 
					totalWorkers=`echo $totalWorkers \+ $useP| bc`
				fi
			else
				echo $host: not enough CPUs available: $useP
			fi
		else
			echo "No access to ${host}"
		fi
	done
	echo "pyHosts=($OkHosts)"
	echo "in total $totalWorkers nodes"
elif [ $1 == "kill" ]


then

	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" ${host} "echo 2>&1" 
		then

			echo "ssh -o 'BatchMode yes' $host killall ppserver.py"
			ssh -o "BatchMode yes" $host killall ppserver.py &
		else
			echo "No access to ${host}"
		fi
	done
	
fi

#!/bin/bash

hosts="albe   euros   refoli   notos   roumet   helm   irifi trombe   habagat   embat caspar"


hosts="roumet habagat"
hosts=`ls /net`

timeout=300

#how many cpus shall be saved?
saveCPUs=100

totalWorkers=0
OkHosts=()

if [ "$1" == "create" ]
then

	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" ${host} "echo 2>&1" 
		then
			:
		else
			echo "No access to ${host}"
			continue
		fi
					
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

	done
	echo "pyHosts=($OkHosts)"
	echo "in total $totalWorkers nodes"
	exit

elif [ "$1" == "test" ]
then

	for host in $hosts
	do
		errors=0
	
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" ${host} "echo > /dev/null 2>&1" 
		then
			echo "${host}: Connected"
		else
			echo "${host}: No access"
			continue
		fi
		
		for module in "numpy" "pp" "pyPamtra" "pyPamtraLib"
		do
			if ssh -o "BatchMode yes" -o "ConnectTimeout 2" $host "python -c 'import $module;'> /dev/null 2>&1"
			then
				echo "${host} : Found $module"
			else
				echo "${host} : No $module"
				errors=$[$errors+1]
				continue
			fi
		done
		if [[ "$errors" -eq 0 ]]
		then
			OkHosts="$host $OkHosts" 
		fi
	done
	echo "hosts=\"$OkHosts\""
	exit

elif [ "$1" == "kill" ]
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
	exit
fi

echo "Options:"
echo "create"
echo "test"
echo "kill"
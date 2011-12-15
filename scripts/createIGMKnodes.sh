#!/bin/bash

hosts=`ls /net`
hosts="zapod viracao trombe tivano tempest seistan schirokko roumet refoli rebat ramier qeb orkan notos nemere molan meltemi kaus karif irifi habagat gargal garbi gales frog forano1 euros embat chinook aure albe"
ppTimeout=1800

#how many cpus shall be saved?
saveCPUs=0.98

user=hatpro

totalWorkers=0
OkHosts=()

if [ "$1" == "create" ]
then

	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" 
		then
			:
		else
			echo "No access to ${host}"
			continue
		fi
		echo "$host: checking no of available cpus"
		noP=`ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host grep processor /proc/cpuinfo | wc|awk '{print $1}'`
		load=`ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host uptime| awk -F, '{print $(NF-1)}'`
		useP=`echo $noP\-$saveCPUs\-\($load\%$noP\) |bc | cut -d '.' -f1`
		if [[ $useP =~ ^-?[0-9]+$ ]]
		then
			:
		else
			useP=0
		fi
		if  [ $useP -gt 0 ]
		then
			useP=`echo $useP |bc`
			echo "$host: nice -n 18 ppserver.py -t $ppTimeout -w $useP &"
# 			export PYTHONPATH=:/home/$user/python/lib:/home/$user/python/libs-local/$host/ && <- now in bashrc of hatpro
			ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "/usr/bin/nice -n 18  ~/python/bin/ppserver.py -t $ppTimeout -w $useP -s 'pyPamtra' &"&
			if [ $? -eq 0 ] 
			then
				OkHosts="\"$host\",$OkHosts" 
				totalWorkers=`echo $totalWorkers \+ $useP| bc`
				echo "$host: OK"
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
	
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo > /dev/null 2>&1" 
		then
			echo "${host}: Connected"
		else
			echo "${host}: No access"
			continue
		fi
		
		for module in "numpy" "pp" "pyPamtra" "pyPamtraLib" "pyPamtraLibWrapper"
		do
			if ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "python -c 'import $module;'> /dev/null 2>&1"
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
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" 
		then

			echo "ssh -o 'BatchMode yes' $user@$host killall ppserver.py"
			ssh -o "BatchMode yes" $user@$host killall ppserver.py &
		else
			echo "No access to ${host}"
		fi
	done
	exit

elif [ "$1" == "killpython" ]
then

	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" 
		then

			echo "ssh -o 'BatchMode yes' $user@$host killall python"
			ssh -o "BatchMode yes" $user@$host killall python &
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
echo "killpython"
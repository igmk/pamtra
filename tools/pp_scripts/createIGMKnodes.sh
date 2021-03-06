#!/bin/bash

source nodesConfig.sh

hosts=$igmkNodes

user="$igmkUser"

port="$igmkPort"
saveCPUs="$igmkSaveCPUs"


ppTimeout=10000
totalWorkers=0
OkHosts=()

killpython ()
{
	host=$1
	if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" 
	then

		echo "ssh -o 'BatchMode yes' $user@$host killall python"
		ssh -o "BatchMode yes" $user@$host killall python &
	else
		echo "No access to ${host}"
	fi
}

killpp ()
{
	host=$1
	if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" 
	then

		echo "ssh -o 'BatchMode yes' $user@$host killall ppserver.py"
		ssh -o "BatchMode yes" $user@$host killall ppserver.py &
	else
		echo "No access to ${host}"
	fi
}

create_parallel () {
	local  __okhosts=$1
	local  __noworker=$2
	if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" $user@${host} "echo 2>&1" > /dev/null 2>&1
	then
		:
	else
# 		echo "No access to ${host}"
		continue
	fi
# 	echo "$host: checking no of available cpus"
	local noP=`ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host grep processor /proc/cpuinfo | wc|awk '{print $1}' 2>/dev/null`
	local load=`ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host uptime| awk -F, '{print $(NF-1)}' 2>/dev/null`
	local useP=`echo $noP\-$saveCPUs\-\($load\%$noP\) |bc | cut -d '.' -f1`
	if [[ $useP =~ ^-?[0-9]+$ ]]
	then
		:
	else
		useP=0
	fi
	if  [ $useP -gt 0 ]
	then
		if ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "python -c 'import pyPamtraLibWrapper;'> /dev/null 2>&1"
		then
# 			echo "${host} : Found pyPamtraLibWrapper"
			local useP=`echo $useP |bc`
# 			echo "$host: /usr/bin/nice -n 18  ~/python/bin/ppserver.py -r -t $ppTimeout -w $useP -s 'pyPamtra' &"
# 				export PYTHONPATH=:/home/$user/python/lib:/home/$user/python/libs-local/$host/ && <- now in bashrc of hatpro
			ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "/usr/bin/nice -n 18  ~/python/bin/ppserver.py -r -t $ppTimeout -w $useP -s 'pyPamtra' -p $port 2>/dev/null &"&
			if [ $? -eq 0 ] 
			then
# 				OkHosts="\"$host\",$OkHosts" 
# 				totalWorkers=`echo $totalWorkers \+ $useP| bc`
				echo -n "\"$host:$port\"," 
			fi
		else
# 			echo "${host} : No pyPamtraLibWrapper"
			:
		fi
	else
# 		echo $host: not enough CPUs available: $useP
		:
	fi
}

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
			if ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "python -c 'import pyPamtraLibWrapper;'> /dev/null 2>&1"
			then
				echo "${host} : Found pyPamtraLibWrapper"
				useP=`echo $useP |bc`
				echo "$host: /usr/bin/nice -n 18  ~/python/bin/ppserver.py -r -t $ppTimeout -w $useP -s 'pyPamtra' &"
# 				export PYTHONPATH=:/home/$user/python/lib:/home/$user/python/libs-local/$host/ && <- now in bashrc of hatpro
				ssh -o "BatchMode yes" -o "ConnectTimeout 2" $user@$host "/usr/bin/nice -n 18  ~/python/bin/ppserver.py -r -t $ppTimeout -w $useP -s 'pyPamtra' -p $port &"&
				if [ $? -eq 0 ] 
				then
					OkHosts="\"$host:$port\",$OkHosts" 
					totalWorkers=`echo $totalWorkers \+ $useP| bc`
					echo "$host: OK"
				fi
			else
				echo "${host} : No pyPamtraLibWrapper"
			fi
		else
			echo $host: not enough CPUs available: $useP
		fi

	done
	echo "pyHosts=($OkHosts)"
	echo "in total $totalWorkers nodes"
	exit

elif [ "$1" == "create_parallel" ]
then
	echo -n "pyHosts=("
	for host in $hosts
	do
		create_parallel $host &
	done
	wait
	echo ")"
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
		killpp $host &
	done
	wait
	exit

elif [ "$1" == "killpython" ]
then

	for host in $hosts
	do
		killpython $host &
	done
	wait
	exit
fi


echo "Options:"
echo "create          - create ppserver"
echo "create_parallel - less output, but runs paralell"
echo "test            - test pyPamtra installation"
echo "kill            - kill ppserver"
echo "killpython      - kill python"

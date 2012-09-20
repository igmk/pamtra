#!/bin/bash

hosts="albe euros refoli notos roumet   helm   irifi trombe   habagat   embat caspar"

#  hosts=`ls /net`
hosts="roumet habagat"
hosts="karif"


# hosts="orkan ramier rebat" 
# 
# hosts="roumet"
# hosts="lomar karif orkan ramier rebat"
# hosts="ramier"
# hosts="albe"

# sudo apt-get install python-dev
# zapod viracao trombe tivano terral tempest taifun seistan schirokko refoli ramier qeb notos molan kaus irifi habagat garbi gales frog forano1 euros embat ecir chinook caspar albe

timeout=300

user=hatpro
branch="rt4"
date=`date "+%Y%m%d"`


OkHosts=()
testFailHosts=()
NoNumpyHosts=()
NotCompileHosts=()
#install pp with:
# python setup.py install --install-lib /home/hatpro/python/lib --install-scripts /home/hatpro/python/bin

if [ "$1" == "prepare" ]
then
	mkdir -p /tmp/pyPamtra
	cd /tmp/pyPamtra
	git clone /home/mmaahn/projects/pamtra -b $branch
	cd /tmp/pyPamtra/pamtra/src/
	make precompile
	tar -czf /tmp/pyPamtra.tar.gz /tmp/pyPamtra
	rm -rf /tmp/pyPamtra
	scp /tmp/pyPamtra.tar.gz $user@localhost:/home/hatpro/python/pyPamtra.tar.gz
	exit
fi


if [ "$1" == "install" ]
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
		if ssh -o "BatchMode yes" -o "ConnectTimeout 2" $host "python -c 'import numpy;'> /dev/null 2>&1"
		then
			:
		else
			echo "############################"
			echo "No numpy found: $host"
			echo "############################"
			NoNumpyHosts="$host $NoNumpyHosts"
			continue
		fi
		ssh -o "BatchMode yes" $user@$host "mkdir -p /tmp/pyPamtraCopy &&tar -xzf /home/hatpro/python/pyPamtra.tar.gz -C /tmp/pyPamtraCopy/"
		#make py
		ssh -o "BatchMode yes" $user@$host "cd /tmp/pyPamtraCopy/tmp/pyPamtra/pamtra/src && make pyNode"
		if [ $? -gt 0 ]
		then
			echo "############################"
			echo "Could not compile: $host"
			echo "############################"
			ssh -o "BatchMode yes" $user@$host rm -rf /tmp/pyPamtraCopy
			NotCompileHosts="$host $NotCompileHosts"
			continue
		fi
		#make pyinstall
		ssh -o "BatchMode yes" $user@$host "mkdir -p ~/python/libs-local/$host/ && mv /tmp/pyPamtraCopy/tmp/pyPamtra/pamtra/py/*so ~/python/libs-local/$host/ && mv /tmp/pyPamtraCopy/tmp/pyPamtra/pamtra/py/*py ~/python/lib/ && mv /tmp/pyPamtraCopy/tmp/pyPamtra/pamtra/py/libs/*/*.py ~/python/lib/"
		#make pytest
		ssh -o "BatchMode yes" $user@$host "export PYTHONPATH=:/home/$user/python/lib:/home/$user/python/libs-local/$host/ && cd /tmp/pyPamtraCopy/tmp/pyPamtra/pamtra/src && make pytest123"
		if [ $? -gt 0 ]
		then
			echo "############################"
			echo "make pytest Fehler bei $host"
			echo "############################"
			ssh -o "BatchMode yes" $user@$host rm -rf /tmp/pyPamtraCopy
			ssh -o "BatchMode yes" $user@$host rm -rvf ~/python/libs-local/$host/
			testFailHosts="$host $testFailHosts"
			continue
		fi
		echo "$host OK"
		OkHosts="$host $OkHosts"
		ssh -o "BatchMode yes" $user@$host rm -rf /tmp/pyPamtraCopy

	done
	echo "hosts=$OkHosts"
	echo "test failed: $testFailHosts"
	echo "not compiled: $NotCompileHosts"
	echo "No Numpy: $NoNumpyHosts"
	exit
fi

if [ "$1" == "installHACK" ]
then
	hosts="albe"
	for host in $hosts
	do
		if ssh -q -q -o "BatchMode=yes" -o "ConnectTimeout 2" ${host} "echo 2>&1" 
		then
			:
		else
			echo "No access to ${host}"
			continue
		fi
		ssh -o "BatchMode yes" $user@$host "mkdir -p ~/python/libs-local/$host/ && cp ~/python/libs-local/karif/* ~/python/libs-local/$host/ "
		if [ $? -gt 0 ]
		then
			echo "############################"
			echo "Could not copy $host"
			echo "############################"
			ssh -o "BatchMode yes" $user@$host rm -rf /tmp/pyPamtraCopy
			NotCompileHosts="$host $NotCompileHosts"
			continue
		fi
		echo "$host OK"
		OkHosts="$host $OkHosts"
	done
	echo "hosts=$OkHosts"
	echo "not copied: $NotCompileHosts"
	exit
fi

if [ "$1" == "remove" ]
then
	ssh -o "BatchMode yes" $user@karif "rm -rvf ~/python/libs-local"
	exit
fi

echo "Options:"
echo "install"
echo "prepare"
echo "remove (machine dependend code only)"

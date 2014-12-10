
servers="aure austre aziab barat bayamos belat  brise cers chinook diablo ecir euros farou fuga  garbi gargal ghibli habagat  helm junta kalmen karif kaus leveche lips lodas magellan melan meltemi metel molan nemere notos orkan ostro purga qeb ramier rebat refoli sarma seistan sharki vento viracao zapod roumet marin"

for server in $servers; 
  do 
  ncpu=$(ssh $server "sh -c 'grep -c processor /proc/cpuinfo'")
  echo $server $ncpu
  for i in $(seq 1 $ncpu)
    do
    ssh -n -f $server "sh -c 'cd /home/mmaahn/projects/pamtra/tools/picklePam; nohup nice -19 ./picklePam1.py /net/marin/mmaahn/pamPickle -m 602 -s2> $server.log 2>&1 &'"
  done
done
date

#for server in $servers; do echo $server && ssh -n -f $server "sh -c 'killall -9 python'"; done; date
#for server in $servers; do echo $server && ssh $server "sh -c 'uptime'"; done; date



#for server in $servers; do echo $server && ssh -n -f $server "sh -c 'cd /home/mmaahn/projects/pamtra/py/; nohup nice -19 ./picklePam.py /net/marin/mmaahn/pamPickle -n -1 -m 300 > $server.log 2>&1 &'"; done; date

#8 cores
#servers="austre habagat nemere roumet junta helm magellan aziab qeb ramier rebat sarma sharki"

#for server in $servers; do echo $server && ssh -n -f $server "sh -c 'cd /home/mmaahn/projects/pamtra/tools/picklePam; nohup nice -19 ./picklePam1.py /net/marin/mmaahn/pamPickle -m 302 -s 2 > $server.log 2>&1 &'"; done; date

#for server in $servers; do echo $server && ssh -n -f $server "sh -c 'grep -c processor /proc/cpuinfo'"; done; date


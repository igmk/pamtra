#!/usr/bin/env python

from __future__ import division
import pyPamtraLibWrapper 
import pyPamtraLib
import cPickle as pickle
import numpy as np
import time
import datetime
import hashlib
import sys
import glob
import os
import multiprocessing
import argparse
import sys, signal, time


parser = argparse.ArgumentParser(description='run pam\'s pickles')
parser.add_argument("pickel_path", help="path to pam pickles")
#parser.add_argument("-v", "--verbosity", action="count", default=0)
parser.add_argument("-m", "--max_age", type=int, default=300, help="maximum idle runtime of service in s (default 300)")
parser.add_argument("-t", "--timeout", type=int, default=60, help="maximum runtime per job in s (default 60)")
parser.add_argument("-n", "--ncpu", type=int, default=0, help="number of cpus, 0 = auto, negative = do not use n cpus. (default 0)")
parser.add_argument("-s", "--sleep", type=int, default=10, help="sleep time before checking directory again (default 10)")
args = parser.parse_args()


print args.pickel_path

maxAge = args.max_age
timeOut = args.timeout
nCpu = args.ncpu
sleepTime = args.ncpu
resetStartTime = True
host = os.uname()[1]
keepRunning = multiprocessing.Value("i",1) #variable to handle termination

if nCpu <= 0: nCpu = multiprocessing.cpu_count() + nCpu

#for server in $(ls /net); do echo $server && ssh -n -f $server "sh -c 'cd /home/mmaahn/projects/pamtra/py/; nohup ./picklePam.py path >> $server.log 2>&1 &'"; done
def processData(fname):
    global keepRunning
    if keepRunning.value: 
      try:
        os.rename(fname,fname+"."+host)
      except OSError:
        print "cannot rename %s probably already taken"%fname
        return
      #try:
      with open(fname+"."+host, 'r') as f:
        inputAr = pickle.load(f)
      #except IOError:
        #print "cannot open %s"%fname
        #continue
      print "running ", fname, inputAr[0]
      try: result = pyPamtraLibWrapper.parallelPamtraFortranWrapper(*inputAr,returnModule=False)
      except:
        print "FAILED", fname, inputAr[0]
        result = [None, host]
      print "done ", fname, inputAr[0]
      fname2 = fname.replace(".job",".result")
      with open(fname2+".tmp", 'w') as f:
        pickle.dump(result, f)
      os.rename(fname2+".tmp",fname2)
      try: os.remove(fname+"."+host)
      except OSError: pass #whoe removed teh file??
    return

startTime = time.time()


#allFiles   = []

def cleanUp(*args,**kwargs):
  global keepRunning
  keepRunning.value = False
  print "I WAS KILLED. Please be patient and wait for the last job to exit:", args
  return


for sig in [signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGQUIT]:
    signal.signal(sig, cleanUp)


    
while ((time.time() - startTime) < maxAge) and keepRunning.value:
  pathIn = "%s/*.job"%args.pickel_path
  nFiles = len(glob.glob(pathIn))
  if nFiles > 0: 
    print "creating pool"
    pool = multiprocessing.Pool(processes=nCpu,maxtasksperchild=100)
    jobs = []

    #if nFiles >100: nFiles=100
    try:
      for fname in sorted(glob.glob(pathIn))[:100]: 
        jobs.append([pool.apply_async(processData,(fname,)),fname])
        #allFiles.append(fname)
      #pool.map(processData,[]*nCpu*10)
      for jj, (job, fname) in enumerate(jobs):
        try:
          job.get(timeout=timeOut) 
          if resetStartTime: startTime = time.time() # reset startTime
        except multiprocessing.TimeoutError:
          print "TIMEOUT", jj
          fname2 = fname.replace(".job",".result")
          os.remove(fname+"."+host)
          with open(fname2+".tmp", 'w') as f:
            pickle.dump([None,host], f)
          os.rename(fname2+".tmp",fname2)  
        #allFiles.remove(fname)
      pool.terminate()
      pool.join()
      del pool, jobs
    except KeyboardInterrupt:
      pool.terminate()
      pool.join()
      del pool, jobs
      print "TERMINATED: KeyboardInterrupt"  
      #cleanUp()
      break
  print "waiting..."
  if keepRunning.value: time.sleep(sleepTime)
  continue      
print "EXIT, time is over", datetime.datetime.now()

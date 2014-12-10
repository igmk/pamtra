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
import argparse
import sys, signal, time
from functools import wraps
import errno

parser = argparse.ArgumentParser(description='run pam\'s pickles')
parser.add_argument("pickel_path", help="path to pam pickles")
#parser.add_argument("-v", "--verbosity", action="count", default=0)
parser.add_argument("-m", "--max_age", type=int, default=300, help="maximum idle runtime of service in s (default 300)")
parser.add_argument("-t", "--timeout", type=int, default=60,  help="maximum runtime per job in s (default 60)")
parser.add_argument("-s", "--sleep",   type=int, default=10,  help="sleep time before checking directory again (default 10)")
args = parser.parse_args()


print args.pickel_path

maxAge = args.max_age
timeOutSec = args.timeout
sleepTime = args.sleep
resetStartTime = True
host = os.uname()[1]
keepRunning = True #variable to handle termination



def cleanUp(*args,**kwargs):
  global keepRunning
  keepRunning = False
  print "I WAS KILLED. Please be patient and wait for the last job to exit:", args
  return



for sig in [signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGQUIT]:
    signal.signal(sig, cleanUp)

class TimeoutError(Exception):
    pass

def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            print("TIMEOUT")
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wraps(func)(wrapper)
    return decorator

#for server in $(ls /net); do echo $server && ssh -n -f $server "sh -c 'cd /home/mmaahn/projects/pamtra/py/; nohup ./picklePam.py path >> $server.log 2>&1 &'"; done
@timeout(timeOutSec)
def processData(fname):
    global keepRunning
    if keepRunning: 
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
      os.remove(fname+"."+host)
    return

startTime = time.time()


#allFiles   = []


    
    



    
while ((time.time() - startTime) < maxAge) and keepRunning:
  pathIn = "%s/*.job"%args.pickel_path
  nFiles = len(glob.glob(pathIn))
  if nFiles > 0: 
    #if nFiles >100: nFiles=100
    try:
      for fname in sorted(glob.glob(pathIn))[:100]: 
        if resetStartTime: startTime = time.time() # reset startTime
        processData(fname)
    except KeyboardInterrupt:
      print "TERMINATED: KeyboardInterrupt"  
      cleanUp()
      break
  if keepRunning: 
    print "waiting..."
    time.sleep(sleepTime)
  continue      
  sys.stdout.flush()

print "EXIT, time is over", datetime.datetime.now()

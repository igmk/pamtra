import pp
import numpy
import sys
##sys.path.remove('/usr/share/pyshared')

#def get_val():
    #return numpy.random.rand()

#job_server = pp.Server(ppservers= ("*",))
#print "Starting pp with", job_server.get_ncpus(), "workers"

#jobs = []
#for i in xrange(500):
    #jobs.append(job_server.submit(get_val,
                                  #(),
                                  #(),
                                  #("numpy",)))
#for i,job in enumerate(jobs):
    #print '%d: %f' % (i,job())
    
    
    
    
ttt = time.time()
    
import sys
import pyPamtra
import numpy as np

inputFile="../test/referenceProfile/testProfiles.dat"

t = pyPamtra.pyPamtra()


t.readProfile("/home/mmaahn/projects/pamtra/profiles/rt_comp.dat")

t.set["data_path"]='/home/mech/models/pamtra/data/'
t.set["verbose"]=0

t.runParallelPamtra([24,30,],pp_servers=("*",),pp_local_workers=4,pp_deltaF=1,pp_deltaX=1,pp_deltaY = 0)
#35,45,55,65,75,85,90,100,120
#t.runParallelPamtra([24,30,35,45,55,65,75,85,90,100,120],pp_servers=(),pp_local_workers=4,pp_deltaF=1,pp_deltaX=1,pp_deltaY = 0)

#t.runPamtra([24,30,35,45,55,65,75,85,90,100,120])

print time.time() - ttt, " s"
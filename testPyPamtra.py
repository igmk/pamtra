import pp
import numpy
import sys
#sys.path.remove('/usr/share/pyshared')

def get_val():
    return numpy.random.rand()

job_server = pp.Server(ppservers= ("*",))
print "Starting pp with", job_server.get_ncpus(), "workers"

jobs = []
for i in xrange(500):
    jobs.append(job_server.submit(get_val,
                                  (),
                                  (),
                                  ("numpy",)))
for i,job in enumerate(jobs):
    print '%d: %f' % (i,job())
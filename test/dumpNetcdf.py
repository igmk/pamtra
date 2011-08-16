import scipy.io.netcdf as nc
import numpy as np
import sys
#dump all netcdf values to stdout.

fnameIn =sys.argv[1]

N = nc.NetCDFFile(fnameIn)
for var in N.variables:
   print np.around(N.variables[var][:],4)

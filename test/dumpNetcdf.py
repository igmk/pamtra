#try: 
	#import scipy.io.netcdf as nc
	#ncOpen = nc.NetCDFFile
#except: 
import netCDF4 as nc
ncOpen = nc.Dataset

import numpy as np
import sys
#dump all netcdf values to stdout.

fnameIn =sys.argv[1]

N = ncOpen(fnameIn)

vars = N.variables.keys()
vars.sort()

for var in vars:
	print var
	try: print np.around(np.array(N.variables[var][:],dtype=float),decimals=4)
	except: print np.array(N.variables[var][:],dtype=str)
N.close()
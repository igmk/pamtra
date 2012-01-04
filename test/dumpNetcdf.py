# -*- coding: utf-8 -*-

try:
	import netCDF4 as nc
	pyNc = True
except:
	try: 
		import Scientific.IO.NetCDF as nc
		pyNc = False
	except:
		import netCDF3 as nc
		pyNc = True
		
	
	
if pyNc: ncOpen = nc.Dataset
else: ncOpen = nc.NetCDFFile

import numpy as np
import sys
#dump all netcdf values to stdout.

fnameIn =sys.argv[1]

N = ncOpen(fnameIn,"r")

vars = N.variables.keys()
vars.sort()

for var in vars:
	print var
	try: print np.around(np.array(N.variables[var][:],dtype=float),decimals=4)
	except: print np.array(N.variables[var][:],dtype=str)
N.close()
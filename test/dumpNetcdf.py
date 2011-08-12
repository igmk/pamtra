import scipy.io.netcdf as nc
import sys
#dump all netcdf values to stdout.

fnameIn =sys.argv[1]

N = nc.NetCDFFile(fnameIn)
for var in N.variables:
   print N.variables[var][:]

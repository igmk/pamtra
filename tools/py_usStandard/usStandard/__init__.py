# -*- coding: utf-8 -*-

# M.Maahn, IGMK 2012

import numpy as np
try:
  from collections.abc import Iterable
except ImportError:
  from collections import Iterable
from .usStandardAtmosphere import atmosphere

   

def usStandard(height):
  """
  Compute the US standard atmosphere                             

  Input  : height      [m]      height as single value or list

  Output : density [kg/m^3]   air density
	  pressure [pa]      air pressure
	  temperature [K]    air temperature
  source: http://www.pdas.com/programs/atmos.f90
  """                                                          

  #check whether height is float or array:
  if isinstance(height, Iterable):
    density = np.ones_like(height,dtype=float)
    pressure = np.ones_like(height,dtype=float)
    temperature = np.ones_like(height,dtype=float)    
    for ii,hh in enumerate(height):
      #Fortran processes only first entry of vector, so make sure it is not a vector
      assert not isinstance(hh, Iterable)
      density[ii], pressure[ii], temperature[ii] = atmosphere(hh/1000.)
  else:
    density, pressure, temperature = atmosphere(height/1000.)
    
  #results of Fortran programm are normed to standard conditions:
  density = density * 1.2250
  pressure = pressure * 101325
  temperature = temperature * 288.15

  return density, pressure, temperature
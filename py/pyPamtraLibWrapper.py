import pyPamtraLib
import os
import logging


logging.basicConfig(filename='/tmp/pyPamtraLibWrapper.log',level=logging.WARNING) #change INFO to DEBUG if needed

#this wrapper is needed because pp cannot work with fortran modules directly. returns results from pamtra AND name of the host (for pp statistics)
def PamtraFortranWrapper(*args):
    
  host = os.uname()[1]
  
  logging.debug('Starting pyPamtraLib... '+host)
  result = pyPamtraLib.pypamtralib(*args)
  logging.debug('Done pyPamtraLib... '+host)

  return result, host


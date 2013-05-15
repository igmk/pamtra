import pyPamtraLib
import os
import logging
import collections
import numpy
import random
import string

logging.basicConfig(filename='/tmp/pyPamtraLibWrapper.log',level=logging.WARNING) #change WARNING to INFO or DEBUG if needed

def PamtraFortranWrapper(nmlSettings,nmlDefaultSettings,nmlFile,*pamtraArgs):
  """
  this wrapper is needed because pp cannot work with fortran modules directly. returns results from pamtra AND name of the host (for pp statistics)
  In addition, it takes care of the nml file generation
  IN
  nmlSettings		nml File Settings for Pamtra 
  nmlDefaultSettings	default nml File Settings for Pamtra 
  nmlFile		nml Filename 
  *pamtraArgs		list of pamtra arguments
  OUT
  result		list of pamtra results
  host			hostname
  """
  host = os.uname()[1]

  if nmlFile == "TMPFILE": 
  
    """
    Write OrderedDict nmlSettings to nmlFile. Type is taken from nmlDefaultSettings
    """
    #create random name
    nmlFileName = nmlSettings["inoutput_mode"]["tmp_path"] +"/pyPamtra_namelist_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nml.tmp"
    logging.debug("opening: "+ nmlFileName)
    f = open(nmlFileName,"w")
    for keygr in nmlSettings.keys():
      if keygr not in nmlDefaultSettings.keys(): 
	logging.warning("Warning can not parse setting: "+str(keygr))
	continue
      f.write("&%s\n\r"%keygr)
      for key in nmlSettings[keygr].keys():
	logging.debug("write: "+ str(keygr) + ": " +str(key))
	if key not in nmlDefaultSettings[keygr].keys(): 
	  logging.warning("Warning can not parse setting: "+str(key))
	  continue
	if type(nmlDefaultSettings[keygr][key])==bool:
	  value = str(nmlSettings[keygr][key]).lower()
	  f.write("%s=.%s.\n\r"%(key,value,))
	elif type(nmlDefaultSettings[keygr][key]) in [int,numpy.int32,numpy.int64]:
	  value = int(nmlSettings[keygr][key])
	  f.write("%s=%i\n\r"%(key,value,))  
	elif type(nmlDefaultSettings[keygr][key]) in [float,numpy.float32,numpy.float64]:
	  value = numpy.float64(nmlSettings[keygr][key])
	  f.write("%s=%f\n\r"%(key,value,))  
	elif type(nmlDefaultSettings[keygr][key]) in [str]:
	  value = str(nmlSettings[keygr][key])
	  f.write("%s='%s'\n\r"%(key,value,))
	else:
	  logging.warning("cannot determine type of nml key "+ key)
      f.write("/\n\r")
    f.close()
  else:
    logging.debug("taking non-temporary nml file: "+nmlFile)
    nmlFileName = nmlFile
    
  logging.debug('Starting pyPamtraLib... '+host)
  logging.debug(nmlFileName)
  try: 
    result = pyPamtraLib.pypamtralib(nmlFileName,*pamtraArgs)
  except:
    logging.warning('PyPamtraLib crashed '+host)
  finally:
    if nmlFile == "TMP": 
      logging.debug("removing "+nmlFileName)
      try: 
	os.remove(nmlFileName)
      except:
        logging.warning("removing "+nmlFileName+" failed!")
  logging.debug('Done pyPamtraLib... '+host)
  return result, host

# -*- coding: utf-8 -*-
from __future__ import division, print_function

import os
#import logging
#import collections
import numpy# as np
#import random
#import string
import copy
import sys

try:
    from .pyPamtraLib import *
except ImportError:
    print('PAMTRA FORTRAN LIBRARY NOT AVAILABLE!')



is3 = sys.version_info[0] == 3

#logging.basicConfig(filename='/tmp/pyPamtraLibWrapper.log',level=logging.WARNING) #change WARNING to INFO or DEBUG if needed

def PamtraFortranWrapper(
  sets,
  nmlSets,
  descriptorFile,
  descriptorFile4D,
  descriptorFileFS,
  profile,
  returnModule=True
  ):
  from . import pyPamtraLib

  report_module.verbose = sets["verbose"]

  #make sure the shape of the profiles is the same!
  for key in list(profile.keys()):
    if type(profile[key]) == numpy.ndarray:
      assert profile[key].shape[0] == profile["lat"].shape[0]
      assert profile[key].shape[1] == profile["lat"].shape[1]
  #todo: make shape check for all variables! seg faults can easily be created here!

  #Pamtra can handle not more than maxlay.
  assert profile["max_nlyrs"] <= settings.maxlay

  #be sure everything is cleaned up before we start
  error = deallocate_everything.do_deallocate_everything()
  if error > 0: raise RuntimeError("Error in deallocate everything")

  settings.in_python = True
  settings.nfrq = len(sets["freqs"])
  settings.freqs[:len(sets["freqs"])] = sets["freqs"]

  settings.settings_fill_default()

  #temporary fixes:
  settings.input_file = "test_mie.dat".ljust(300)

  if 'salinity' in nmlSets:
    raise DeprecationWarning("nmlSets['salinity'] is deprecated. Use 2D profile 'sfc_salinity' instead.")
  if 'ground_type' in nmlSets:
    raise DeprecationWarning("nmlSets['ground_type'] is deprecated. Use 2D profile 'sfc_refl' instead.")

  #loop through settings

  for key in list(nmlSets.keys()):
    isList = getattr(settings, key.lower()).size > 1
    if type(nmlSets[key]) == str:
      if sets["pyVerbose"] > 3: print(("settings."+key.lower() +"[:] = '" + str(nmlSets[key])+"'"))
      setattr(settings, key.lower(), nmlSets[key].ljust(lenFortStrAr(getattr(settings,key.lower()))  ))
    else:
      if isList:
        if sets["pyVerbose"] > 3: print(("settings."+key.lower() +"[:] = numpy.array(nmlSets['"+key+"']).tolist()"))
        getattr(settings, key.lower())[0:numpy.asarray(nmlSets[key]).size] = numpy.asarray(nmlSets[key]).flatten()
      else:
        if sets["pyVerbose"] > 3: print(("settings."+key.lower() +" = " + str(nmlSets[key])))
        setattr(settings, key.lower(), nmlSets[key])

  #see whether it worked:
  if sets["pyVerbose"] > 3:
    print("Fortran view on settings variables")
    settings.print_settings()


  #deal with the descriptor_file
  descriptor_file.n_hydro = len(descriptorFile)

  #allocation of string array does not work via f2py. Thus we allocate the arrays in Fortran:
  error = descriptor_file.allocate_descriptor_file()
  #(looks like we can reallocate some variables later to include 4D data)
  if error > 0: raise RuntimeError("Error in allocate_descriptor_file")

  descFileCharLength=15

  for name in descriptorFile.dtype.names:
    #1D data
    if name in ["moment_in","liq_ice"]:#,
      if sets["pyVerbose"] > 3: print(("descriptor_file."+name +"_arr = descriptorFile['"+name+"'].tolist()"))
      setattr(descriptor_file, name +"_arr", descriptorFile[name].tolist())
    #1d Strings, these are ugly...
    elif name in ["hydro_name","dist_name","scat_name","vel_size_mod"]:
      if sets["pyVerbose"] > 3: print("setattr(descriptor_file, name+"+_str+", ','.join(descriptorFile[name])))")
      thisStr = ','.join(descriptorFile[name])
      maxLen = lenFortStrAr(getattr(descriptor_file, name+"_str"))
      assert len(thisStr) <= maxLen
      setattr(descriptor_file, name+"_str", thisStr.ljust(maxLen))
    #potential 4D data
    else:
      if sets["pyVerbose"] > 3: print(("descriptor_file."+name +"_arr = [[[descriptorFile['"+name+"'].tolist()]]]"))
      setattr(descriptor_file, name +"_arr", [[[descriptorFile[name].tolist()]]])
  for name4d in list(descriptorFile4D.keys()):
    assert descriptorFile4D[name4d].shape[0] == profile["lat"].shape[0]
    assert descriptorFile4D[name4d].shape[1] == profile["lat"].shape[1]
    if sets["pyVerbose"] > 3: print(("descriptor_file."+name4d +"_arr = descriptorFile4D['"+name4d+"'].tolist()"))
    setattr(descriptor_file, name4d +"_arr", descriptorFile4D[name4d].tolist())


  descriptor_file.process_descriptor_file()

  #see whether it worked:
  if sets["pyVerbose"] > 3:
    print("Fortran view on descriptor_file variables")
    descriptor_file.printdescriptorvars()


  vars_atmosphere.atmo_ngridx = profile["ngridx"]
  vars_atmosphere.atmo_ngridy = profile["ngridy"]
  vars_atmosphere.atmo_max_nlyrs = profile["max_nlyrs"]

  error = vars_atmosphere.allocate_atmosphere_vars()
  if error > 0: raise RuntimeError("Error in allocate_vars_atmosphere")

  if sets["pyVerbose"] > 8:
    for key in list(profile.keys()):
      if key not in ["noutlevels"]:
        try: 
          print(key, getattr(vars_atmosphere, "atmo_"+key), profile[key])
        except AttributeError:
          print(key, getattr(vars_atmosphere, key), profile[key])


  #return  dict(),pyPamtraLib
  #deal with the atmospheric input_file
  for key in list(profile.keys()):

    assert type(profile[key]) != numpy.ma.core.MaskedArray

    if key in ["ngridx","ngridy","max_nlyrs"]:
      continue
    elif key in ["noutlevels"]:
      if sets["pyVerbose"] > 3: print(("settings."+key +" = profile['"+key+"']"))
      setattr(settings, key, profile[key])
    elif key in ["sfc_type","sfc_model","sfc_salinity","sfc_slf","sfc_sif"]:
      if sets["pyVerbose"] > 3: print(("vars_atmosphere."+key +" = profile['"+key+"']"))
      setattr(vars_atmosphere, key, profile[key])
    elif key in ["sfc_refl"]:
      if sets["pyVerbose"] > 3: print(("vars_atmosphere."+key +" = profile['"+key+"'])"))
      setattr(vars_atmosphere, key, profile[key])
    elif type(profile[key]) in [int, float, str]:
      if sets["pyVerbose"] > 3: print(("vars_atmosphere.atmo_"+key +" = profile['"+key+"'].tolist()"))
      setattr(vars_atmosphere, "atmo_"+key, profile[key])
    elif type(profile[key]) == numpy.ndarray:
      if sets["pyVerbose"] > 3: print(("vars_atmosphere.atmo_"+key +" = profile['"+key+"'].tolist()"))
      setattr(vars_atmosphere, "atmo_"+key, profile[key].tolist())
    else:
      raise TypeError("do not understand type of "+ key+": " + str(type(profile[key])))
    #vars_atmosphere.atmo_max_nlyr

  if sets["pyVerbose"] > 8:
    for key in list(profile.keys()):
      if key not in ["noutlevels"]:
        print(key, getattr(vars_atmosphere, "atmo_"+key, profile[key]))
  #see whether it worked:
  if sets["pyVerbose"] > 3:
    print("Fortran view on vars_atmosphere variables")
    vars_atmosphere.print_vars_atmosphere()


  error = vars_atmosphere.fillmissing_atmosphere_vars()
  if error > 0: raise RuntimeError("Error in fillmissing_atmosphere_vars")

  #see whether it worked:
  if sets["pyVerbose"] > 3:
    print("Fortran view on vars_atmosphere variables")
    vars_atmosphere.print_vars_atmosphere()


  if nmlSets["hydro_fullspec"]:
    assert descriptorFileFS
    #assert nmlSets["save_psd"] = False #because there is nothing to save!

    error = vars_hydrofullspec.allocate_hydrofs_vars(descriptorFileFS["d_ds"].shape[-1])
    if error > 0: raise RuntimeError("Error in allocate_hydrofs_vars")

    for key in list(descriptorFileFS.keys()):
      assert descriptorFileFS[key].shape[0] == profile["lat"].shape[0]
      assert descriptorFileFS[key].shape[1] == profile["lat"].shape[1]
      if sets["pyVerbose"] > 3: print(("vars_hydrofullspec.hydrofs_"+key +" = descriptorFileFS['"+key+"'].tolist()"))
      setattr(vars_hydrofullspec, "hydrofs_"+key, descriptorFileFS[key].tolist())

    if sets["pyVerbose"] > 3:
      print("Fortran view on hydro_fullspec variables")
      vars_hydrofullspec.print_hydrofs_vars()

  ##now, finally rund the model
  pamError = pypamtralib.run_pamtra()
  #if error > 0: raise RuntimeError("Error in run_pamtra")

  ##process the results!
  results = dict()
  for key in ["tb","Ze","emissivity","Att_hydro","Att_atmo","radar_hgt","radar_moments","radar_edges","radar_slopes","radar_quality","radar_snr", "radar_spectra","radar_vel","psd_d","psd_deltad","psd_n","psd_mass","psd_area","psd_bscat","kextatmo","scatter_matrix","extinct_matrix","emis_vector","angles_deg"]:
    if sets["pyVerbose"] > 3: print(("allocTest = vars_output.out_"+key.lower()+" is None"))
    allocTest = getattr(vars_output, "out_"+key.lower()) is None
    if not allocTest:
      if sets["pyVerbose"] > 3: print(("results['"+key+"'] = copy.deepcopy(vars_output.out_"+key.lower()+")"))
      results[key] = copy.deepcopy(getattr(vars_output, "out_"+key.lower()))
    else:
      if sets["pyVerbose"] > 3: print("filling key", key)
      if key in ["radar_quality"]: results[key] = -9999
      else: results[key] = -9999.


  results["pamtraVersion"] = str(pypamtralib.gitversion.astype('U40')).strip()
  results["pamtraHash"] = str(pypamtralib.githash.astype('U40')).strip()  

  if sets["pyVerbose"] > 2: "processed results"

  if returnModule:
    return results,pamError, pyPamtraLib
  else:
    del pyPamtraLib
    return results,pamError


def lenFortStrAr(arr):
  '''
  get string length of a fortran string array

  This change is probably related to numpy 1.12 -> 1.16
  '''
  try:
    return len(arr)
  except TypeError:
    return arr.dtype.itemsize




# def _str_py2f(array,length=None):
#   # the byte order of fortran and numpy string arrays is different, this here works sometimes...
#   #if len(array.shape) > 2: raise NotImplemented("Can only handle 1D lists of strings")
#   if length is None: length = array.shape[1]
#   return numpy.lib.stride_tricks.as_strided(array,strides=(length,1))

# def _strList2charArray(strList,charLength=None,arrayLength=None):
#   #makes from list strList an aray of type "s1" that Fortran can handle it.
#   if arrayLength:
#     dim1 = arrayLength
#   else:
#     dim1 = len(strList)
#   if charLength:
#     dim2 = charLength
#   else:
#     dim2 = numpy.max(list(map(len,strList)))
#   charArray = numpy.zeros((dim1,dim2),dtype="S1")
#   for ss,string in enumerate(strList):
#     charArray[ss,:len(string)] = list(string)
#     charArray[ss,len(string):]= " "
#   return charArray

# def setFortranStrList(fortranList,pythonList,charLength=None):
#   #this routine takes care of all the oddities if yo transfer a list of strings from python to fortran
#   if len(fortranList.shape) > 2: raise NotImplemented("Can only handle 1D lists of strings and fortranList must be allocated")
#   if not charLength:
#     charLength = fortranList.shape[1]
#   #pythonList = _strList2charArray(pythonList,charLength=charLength)
#   for pp,pythonStr in enumerate(pythonList):
#     _str_py2f(fortranList)[pp][0:len(pythonStr)] = list(pythonStr)
#     #we have to fill teh rest of the variable with spaces, otherwise it contains only random!
#     _str_py2f(fortranList)[pp][len(pythonStr):] = " "
#   return

def parallelPamtraFortranWrapper(indices, *args, **kwargs):
  if args[0]["pyVerbose"] > 1: print('starting', __name__, 'parent process:', os.getppid(), 'process id:', os.getpid())
  results, pamError = PamtraFortranWrapper(*args, **kwargs)
  host = os.uname()[1]
  return indices, results, pamError, host
  #return indices, dict()


#def PamtraFortranWrapper_OLD(nmlSets,nmlDefaultSettings,nmlFile,*pamtraArgs):
  #"""
  #this wrapper is needed because pp cannot work with fortran modules directly. returns results from pamtra AND name of the host (for pp statistics)
  #In addition, it takes care of the nml file generation
  #IN
  #nmlSets		nml File Settings for Pamtra
  #nmlDefaultSettings	default nml File Settings for Pamtra
  #nmlFile		nml Filename
  #*pamtraArgs		list of pamtra arguments
  #OUT
  #result		list of pamtra results
  #host			hostname
  #"""
  #host = os.uname()[1]

  #if nmlFile == "TMPFILE":

    #"""
    #Write OrderedDict nmlSets to nmlFile. Type is taken from nmlDefaultSettings
    #"""
    ##create random name
    #logging.debug("opening: "+ nmlFileName)
    #f = open(nmlFileName,"w")
    #for keygr in nmlSets.keys():
      #if keygr not in nmlDefaultSettings.keys():
	#logging.warning("Warning can not parse setting: "+str(keygr))
	#continue
      #f.write("&%s\n\r"%keygr)
      #for key in nmlSets[keygr].keys():
	#logging.debug("write: "+ str(keygr) + ": " +str(key))
	#if key not in nmlDefaultSettings[keygr].keys():
	  #logging.warning("Warning can not parse setting: "+str(key))
	  #continue
	#if type(nmlDefaultSettings[keygr][key])==bool:
	  #value = str(nmlSets[keygr][key]).lower()
	  #f.write("%s=.%s.\n\r"%(key,value,))
	#elif type(nmlDefaultSettings[keygr][key]) in [int,numpy.int32,numpy.int64]:
	  #value = int(nmlSets[keygr][key])
	  #f.write("%s=%i\n\r"%(key,value,))
	#elif type(nmlDefaultSettings[keygr][key]) in [float,numpy.float32,numpy.float64]:
	  #value = numpy.float64(nmlSets[keygr][key])
	  #f.write("%s=%f\n\r"%(key,value,))
	#elif type(nmlDefaultSettings[keygr][key]) in [str]:
	  #value = str(nmlSets[keygr][key])
	  #f.write("%s='%s'\n\r"%(key,value,))
	#else:
	  #logging.warning("cannot determine type of nml key "+ key)
      #f.write("/\n\r")
    #f.close()
  #else:
    #logging.debug("taking non-temporary nml file: "+nmlFile)
    #nmlFileName = nmlFile

  #logging.debug('Starting .. '+host)
  #logging.debug(nmlFileName)
  #try:
    #result = pypamtralib(nmlFileName,*pamtraArgs)
  #except:
    #logging.warning('PyPamtraLib crashed '+host)
  #finally:
    #if nmlFile == "TMP":
      #logging.debug("removing "+nmlFileName)
      #try:
	#os.remove(nmlFileName)
      #except:
        #logging.warning("removing "+nmlFileName+" failed!")
  #logging.debug('Done .. '+host)
  #return result, host

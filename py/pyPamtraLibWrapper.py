import pyPamtraLib
import os
#import logging
#import collections
import numpy as np
#import random
#import string

#logging.basicConfig(filename='/tmp/pyPamtraLibWrapper.log',level=logging.WARNING) #change WARNING to INFO or DEBUG if needed

def PamtraFortranWrapper(
  settings,
  nmlSettings,
  descriptorFile,
  profile,
  returnModule=True
  ):
    
  pyPamtraLib.report_module.verbose = settings["verbose"]
   
       
  #be sure everything is cleaned up before we start
  error = pyPamtraLib.deallocate_everything.do_deallocate_everything()
  if error > 0: raise RuntimeError("Error in deallocate everything")
  
  pyPamtraLib.settings.settings_fill_default()    
  
  pyPamtraLib.settings.in_python = True
  pyPamtraLib.settings.nfrq = len(settings["freqs"])
  pyPamtraLib.settings.freqs[:len(settings["freqs"])] = settings["freqs"]
  
  #temporary fixes:
  pyPamtraLib.settings.input_file[:] = "test_mie.dat"

  #loop through settings
  for key in nmlSettings.keys():
    exec("foo = pyPamtraLib.settings."+key)
    if type(nmlSettings[key]) == str:
      if settings["pyVerbose"] > 3: print("pyPamtraLib.settings."+key +"[:] = '" + str(nmlSettings[key])+"'") 
      exec("pyPamtraLib.settings."+key +"[:] = '" + nmlSettings[key]+"'")
    else:
      if settings["pyVerbose"] > 3: print("pyPamtraLib.settings."+key +" = " + str(nmlSettings[key]))
      exec("pyPamtraLib.settings."+key +" = " + str(nmlSettings[key]))
    
  #see whether it worked:
  if settings["pyVerbose"] > 3:
    print "Fortran view on settings variables"
    pyPamtraLib.settings.print_settings()
    

  #deal with the descriptor_file
  pyPamtraLib.descriptor_file.n_hydro = len(descriptorFile.data)
  
  #allocation of string array does not work via f2py. Thus we allocate the arrays in Fortran:
  error = pyPamtraLib.descriptor_file.allocate_descriptor_file()    
  #(looks like we can reallocate some variables later to include 4D data)
  if error > 0: raise RuntimeError("Error in allocate_descriptor_file")

  descFileCharLength=15
  
  for name in descriptorFile.data.dtype.names:
    #1D data
    if name in ["moment_in","liq_ice"]:#,
      if settings["pyVerbose"] > 3: print("pyPamtraLib.descriptor_file."+name +"_arr = descriptorFile.data['"+name+"'].tolist()")
      exec("pyPamtraLib.descriptor_file."+name +"_arr = descriptorFile.data['"+name+"'].tolist()")
    #1d Strings, these are ugly...
    elif name in ["hydro_name","dist_name","scat_name","vel_size_mod"]:
      if settings["pyVerbose"] > 3: print("setFortranStrList(pyPamtraLib.descriptor_file."+name+"_arr,descriptorFile.data['"+name+"'])")
      exec("setFortranStrList(pyPamtraLib.descriptor_file."+name+"_arr,descriptorFile.data['"+name+"'])")
    #potential 4D data
    else:
      if settings["pyVerbose"] > 3: print("pyPamtraLib.descriptor_file."+name +"_arr = [[[descriptorFile.data['"+name+"'].tolist()]]]")
      exec("pyPamtraLib.descriptor_file."+name +"_arr = [[[descriptorFile.data['"+name+"'].tolist()]]]")
  for name4d in descriptorFile.data4D.keys():
    if settings["pyVerbose"] > 3: print("pyPamtraLib.descriptor_file."+name4d +"_arr = descriptorFile.data4D['"+name4d+"'].tolist()")
    exec("pyPamtraLib.descriptor_file."+name4d +"_arr = descriptorFile.data4D['"+name4d+"'].tolist()")
  
  #see whether it worked:
  if settings["pyVerbose"] > 3:
    print "Fortran view on descriptor_file variables"
    pyPamtraLib.descriptor_file.printdescriptorvars()
    
    
  pyPamtraLib.vars_atmosphere.atmo_ngridx = profile["ngridx"]  
  pyPamtraLib.vars_atmosphere.atmo_ngridy = profile["ngridy"]  
  pyPamtraLib.vars_atmosphere.atmo_max_nlyrs = profile["max_nlyrs"]  
  
  error = pyPamtraLib.vars_atmosphere.allocate_atmosphere_vars()  
  if error > 0: raise RuntimeError("Error in allocate_vars_atmosphere")
  
  #return  dict(),pyPamtraLib    
  #deal with the atmospheric input_file
  for key in profile.keys():
    if key in ["ngridx","ngridy","max_nlyrs"]:
      continue
    
    elif type(profile[key]) in [int, float, str]:
      if settings["pyVerbose"] > 3: print("pyPamtraLib.vars_atmosphere.atmo_"+key +" = profile['"+key+"'].tolist()")
      exec("pyPamtraLib.vars_atmosphere.atmo_"+key +" = profile['"+key+"']")
    elif type(profile[key]) == np.ndarray:
      if settings["pyVerbose"] > 3: print("pyPamtraLib.vars_atmosphere.atmo_"+key +" = profile['"+key+"'].tolist()")
      exec("pyPamtraLib.vars_atmosphere.atmo_"+key +" = profile['"+key+"'].tolist()")
    #pyPamtraLib.vars_atmosphere.atmo_max_nlyr
  
  pyPamtraLib.vars_atmosphere.fillmissing_atmosphere_vars()  
  
  #see whether it worked:
  if settings["pyVerbose"] > 3:
    print "Fortran view on vars_atmosphere variables"
    pyPamtraLib.vars_atmosphere.print_vars_atmosphere()
    
    
  ##now, finally rund the model  
  error = pyPamtraLib.pypamtralib.run_pamtra() 
  #if error > 0: raise RuntimeError("Error in run_pamtra")
  
  ##process the results!
  results = dict()
  for key in ["tb","Ze","Att_hydro","Att_atmo","radar_hgt","radar_moments","radar_slopes","radar_quality","radar_snr", "radar_spectra","radar_vel","psd_d_bound","psd_f","psd_mass","psd_area","angles_deg"]:
    if settings["pyVerbose"] > 3: print("allocTest = pyPamtraLib.vars_output.out_"+key.lower()+" == None")
    exec("allocTest = pyPamtraLib.vars_output.out_"+key.lower()+" == None")
    if not allocTest:
      if settings["pyVerbose"] > 3: print("results['"+key+"'] = pyPamtraLib.vars_output.out_"+key.lower())
      exec("results['"+key+"'] = pyPamtraLib.vars_output.out_"+key.lower())
    else:
      if settings["pyVerbose"] > 3: print "filling key", key
      if key in ["radar_quality"]: results[key] = -9999
      else: results[key] = -9999.
    
    
  results["pamtraVersion"] = "".join(list(pyPamtraLib.pypamtralib.gitversion)).strip()      
  results["pamtraHash"] = "".join(list(pyPamtraLib.pypamtralib.githash)).strip()  
  
  if settings["pyVerbose"] > 2: "processed results"
        
  if returnModule: return results, pyPamtraLib
  else: return results

def _str_py2f(array,length=None):
  # the byte order of fortran and numpy string arrays is different, this here works sometimes...
  #if len(array.shape) > 2: raise NotImplemented("Can only handle 1D lists of strings")
  if length==None: length = array.shape[1]
  return np.lib.stride_tricks.as_strided(array,strides=(length,1))
  
def _strList2charArray(strList,charLength=None,arrayLength=None):
  #makes from list strList an aray of type "s1" that Fortran can handle it.
  if arrayLength:
    dim1 = arrayLength
  else:  
    dim1 = len(strList)
  if charLength:
    dim2 = charLength
  else:
    dim2 = np.max(map(len,strList))
  charArray = np.zeros((dim1,dim2),dtype="S1")
  for ss,string in enumerate(strList):
    charArray[ss,:len(string)] = list(string)
    charArray[ss,len(string):]= " "
  return charArray
  
def setFortranStrList(fortranList,pythonList,charLength=None):
  #this routine takes care of all the oddities if yo transfer a list of strings from python to fortran
  if len(fortranList.shape) > 2: raise NotImplemented("Can only handle 1D lists of strings and fortranList must be allocated")
  if not charLength:
    charLength = fortranList.shape[1]
  #pythonList = _strList2charArray(pythonList,charLength=charLength)
  for pp,pythonStr in enumerate(pythonList):
    _str_py2f(fortranList)[pp][0:len(pythonStr)] = list(pythonStr) 
    #we have to fill teh rest of the variable with spaces, otherwise it contains only random!
    _str_py2f(fortranList)[pp][len(pythonStr):] = " "
  return
  
def parallelPamtraFortranWrapper(indices, *args, **kwargs):
  print 'starting', __name__, 'parent process:', os.getppid(), 'process id:', os.getpid()
  results = PamtraFortranWrapper(*args, **kwargs)
  return indices, results
  #return indices, dict()
  
  
#def PamtraFortranWrapper_OLD(nmlSettings,nmlDefaultSettings,nmlFile,*pamtraArgs):
  #"""
  #this wrapper is needed because pp cannot work with fortran modules directly. returns results from pamtra AND name of the host (for pp statistics)
  #In addition, it takes care of the nml file generation
  #IN
  #nmlSettings		nml File Settings for Pamtra 
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
    #Write OrderedDict nmlSettings to nmlFile. Type is taken from nmlDefaultSettings
    #"""
    ##create random name
    #nmlFileName = nmlSettings["inoutput_mode"]["tmp_path"] +"/pyPamtra_namelist_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nml.tmp"
    #logging.debug("opening: "+ nmlFileName)
    #f = open(nmlFileName,"w")
    #for keygr in nmlSettings.keys():
      #if keygr not in nmlDefaultSettings.keys(): 
	#logging.warning("Warning can not parse setting: "+str(keygr))
	#continue
      #f.write("&%s\n\r"%keygr)
      #for key in nmlSettings[keygr].keys():
	#logging.debug("write: "+ str(keygr) + ": " +str(key))
	#if key not in nmlDefaultSettings[keygr].keys(): 
	  #logging.warning("Warning can not parse setting: "+str(key))
	  #continue
	#if type(nmlDefaultSettings[keygr][key])==bool:
	  #value = str(nmlSettings[keygr][key]).lower()
	  #f.write("%s=.%s.\n\r"%(key,value,))
	#elif type(nmlDefaultSettings[keygr][key]) in [int,numpy.int32,numpy.int64]:
	  #value = int(nmlSettings[keygr][key])
	  #f.write("%s=%i\n\r"%(key,value,))  
	#elif type(nmlDefaultSettings[keygr][key]) in [float,numpy.float32,numpy.float64]:
	  #value = numpy.float64(nmlSettings[keygr][key])
	  #f.write("%s=%f\n\r"%(key,value,))  
	#elif type(nmlDefaultSettings[keygr][key]) in [str]:
	  #value = str(nmlSettings[keygr][key])
	  #f.write("%s='%s'\n\r"%(key,value,))
	#else:
	  #logging.warning("cannot determine type of nml key "+ key)
      #f.write("/\n\r")
    #f.close()
  #else:
    #logging.debug("taking non-temporary nml file: "+nmlFile)
    #nmlFileName = nmlFile
    
  #logging.debug('Starting pyPamtraLib... '+host)
  #logging.debug(nmlFileName)
  #try: 
    #result = pyPamtraLib.pypamtralib(nmlFileName,*pamtraArgs)
  #except:
    #logging.warning('PyPamtraLib crashed '+host)
  #finally:
    #if nmlFile == "TMP": 
      #logging.debug("removing "+nmlFileName)
      #try: 
	#os.remove(nmlFileName)
      #except:
        #logging.warning("removing "+nmlFileName+" failed!")
  #logging.debug('Done pyPamtraLib... '+host)
  #return result, host

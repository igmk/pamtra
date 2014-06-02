# -*- coding: utf-8 -*-

#todo: function to calculate PIA?
from __future__ import division
import numpy as np
import datetime
import random
import os
import sys
import string
import glob
import warnings

try: 
  import numexpr as ne
  neAvail = True
except:
  warnings.warn("numexpr not available", Warning)
  neAvail = False


import pyPamtra
import meteoSI

#import csv
#import pickle
#import time,calendar,datetime
#import warnings
#import sys
#import os
from copy import deepcopy

missingNumber =-9999.

def readWrfDataset(fname,kind):
  import netCDF4
  
  '''
  import wrf Dataset with fname of kind
  '''
  
  if kind not in ["gce"]:
    raise TypeError("unknown wrf data type")
  
  variables = ['XLONG', 'XLAT', 'U10', 'V10', 'TSK', 'PSFC', 'T', 'P', 'PB', 'QVAPOR', 'QCLOUD', 'QICE', 'QRAIN', 'QGRAUP', 'QSNOW', 'LANDMASK', 'PH', 'PHB']

  data = dict()

  ncFile = netCDF4.Dataset(fname,"r")

  timeStr = "".join(ncFile.variables["Times"][:].tolist()[0])
  data["Times"] = netCDF4.date2num(datetime.datetime.strptime(timeStr,"%Y-%m-%d_%H:%M:%S"),"seconds since 1970-01-01 00:00:00")
  for var in variables:    
    data[var] = ncFile.variables[var][0,:].T
  ncFile.close()

  data["PH"] = (data["PH"]+data["PHB"])/9.81 #now "real m
  del data["PHB"]
  
  data["P"] = data["P"]+data["PB"] # now in Pa
  del data["PB"]

  Cp = 7.*meteoSI.Rair/2. 
  RdCp =  meteoSI.Rair/Cp
  p0 = 100000
  data["T"] = data["T"] + 300 # now in K
  data["T"] = data["T"]*(data["P"]/p0)**RdCp # potential temp to temp
  
  #add surface data
  data["TSK"] = data["TSK"].reshape(np.append(list(np.shape(data["TSK"])),1))
  data["PSFC"] = data["PSFC"].reshape(np.append(list(np.shape(data["PSFC"])),1))

  data["T"] = np.concatenate((data["TSK"],data["T"]),axis=-1)
  data["P"] = np.concatenate((data["PSFC"],data["P"]),axis=-1)
  
  #translate variables
  varPairs = [["Times","timestamp"],["XLAT","lat"],["XLONG","lon"],["LANDMASK","lfrac"],["U10","wind10u"],["V10","wind10v"],["PH","hgt_lev"],["P","press_lev"],["T","temp_lev"],["QVAPOR","q"],["QCLOUD","cwc_q"],["QICE","iwc_q"],["QRAIN","rwc_q"],["QSNOW","swc_q"],["QGRAUP","gwc_q"]]
  
  pamData = dict()
  for wrfVar,pamVar in varPairs:
    pamData[pamVar] = data[wrfVar]
  
  pam = pyPamtra.pyPamtra()
  pam.createProfile(**pamData)
  del data
  
  return pam

  
def readCosmoDe1MomDataset(fnames,kind,forecastIndex = 1,colIndex=0,tmpDir="/tmp/",fnameInTar="",concatenateAxis=1,debug=False,verbosity=0,df_kind="default",constantFields=None):
  '''
  import wrf Dataset with fname of kind 
  
  fnames = str , fileNames, wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file
  kind = kind of Cosmo file, right now only collum netcdf files are implemented.
  forecastIndex = 1 #take the forecast being between 3 and 5.75 hours old.
  colIndex: which collum should be taken? list allowed!
  fnameInTar = if nc.gz file in tar file, name of nc.gz file (wildcards allowed!)
  concatenateAxis: axis to concatenate the grids
  debug: stop and load debugger on exception
  '''
  import netCDF4
  
  if kind == "gop_collumn":

    variables1Dx = ["time1h"]
    variables1Dy = ["fr_land","latitude","longitude",]
    variables2D = [ "hfl","hhl"]
    variables3D = ["u_10m","v_10m","t_2m","surface_air_pressure","t_s"]
    variables4D = ["temperature","p","qv","qc","qi","qi","qr","qs","qg"]
    
    nHydro = 5
    
    data = dict()
    
    files = np.sort(glob.glob(fnames))
    if len(files) == 0: raise RuntimeError( "no files found")
    files.sort()
    
    ffOK = 0 #successfull runs
    
    for ff, fname in enumerate(files):
      if verbosity>0: print fname
      try:
        if fname.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fname.split(".")[-1]=="tar":
            if verbosity>3:print "tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz"
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fname)
            if verbosity>2:print "created ", gzFile
          else:
            gzFile = fname
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print "created ", tmpFile
          ncFile = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print "opend ", tmpFile    
        else:  
          ncFile = netCDF4.Dataset(fname,"r")
          if verbosity>1:print "opend ", fname

        #import pdb;pdb.set_trace()
        dataSingle = dict()  
        for var in variables1Dx:  
          dataSingle[var] = ncFile.variables[var][:]
        for var in variables1Dy:  
          dataSingle[var] = ncFile.variables[var][[colIndex]]
        for var in variables2D:  
          dataSingle[var] = ncFile.variables[var][[colIndex],::-1] #reverse height order
        for var in variables3D:  
          dataSingle[var] = ncFile.variables[var][[colIndex],forecastIndex,:]
        for var in variables4D:
          dataSingle[var] = np.swapaxes(ncFile.variables[var][[colIndex],:,forecastIndex,:],1,2)[...,::-1]#reverse height order  
          
        ncFile.close()
        if verbosity>1:print "closed nc"
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
          os.system("rm -f "+tmpFile+"*")
        shape2D = (np.shape(dataSingle["latitude"])[0],np.shape(dataSingle["time1h"])[0],)
        shape3Dplus = (np.shape(dataSingle["latitude"])[0],np.shape(dataSingle["time1h"])[0],np.shape(dataSingle["temperature"])[2]+1)
        shape3D = (np.shape(dataSingle["latitude"])[0],np.shape(dataSingle["time1h"])[0],np.shape(dataSingle["temperature"])[2])
        
        time1h = np.zeros(shape2D)
        time1h[:] = dataSingle["time1h"]
        latitude = np.zeros(shape2D)
        latitude.T[:] = dataSingle["latitude"]
        longitude = np.zeros(shape2D)
        longitude.T[:] = dataSingle["longitude"]
        fr_land =  np.zeros(shape2D)
        fr_land.T[:] = dataSingle["fr_land"]
        hfl =  np.zeros(shape3D)
        hhl =  np.zeros(shape3Dplus)
        for t in xrange(shape2D[1]):
          hfl[:,t,:] = dataSingle["hfl"]
          hhl[:,t,:] = dataSingle["hhl"]
          
        dataSingle["time1h"] = time1h
        dataSingle["latitude"] = latitude
        dataSingle["longitude"] = longitude
        dataSingle["fr_land"] = fr_land
        dataSingle["hhl"] = hhl
        dataSingle["hfl"] = hfl
        
        del time1h, longitude, latitude, fr_land, hhl, hfl
        
        
        if ffOK == 0: #if the first file is broken, checking for ff==0 would fail!
          data = deepcopy(dataSingle)
        else:
          for key in data.keys():
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)
          
        ffOK += 1  
        
      #except IOError:
      except Exception as inst:
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
          os.system("rm -f "+tmpFile+"*")
        print "ERROR:", fname      
        print type(inst)     # the exception instance
        print inst.args      # arguments stored in .args
        print inst
        if debug: import pdb;pdb.set_trace()



    #shapes have changed!
    shape2D = np.shape(data["t_2m"])
    shape3Dplus = np.shape(data["hhl"])
    shape3D = np.shape(data["hfl"])
      
    #data["p_hl"] = np.zeros(shape3Dplus)
    #data["t_hl"] = np.zeros(shape3Dplus)
    
    #data["t_hl"][...,0] = data["t_2m"]
    #data["t_hl"][...,1:-1] = data["temperature"][...,0:-1]+0.5*np.diff(data["temperature"],axis=-1)
    #data["t_hl"][...,-1] = data["temperature"][...,-1] + 0.25*(data["temperature"][...,-1]-data["temperature"][...,-2])
    
    #data["p_hl"][...,0] = data["surface_air_pressure"]

    #press1 = deepcopy(data["p"])
    #press1[press1==missingNumber]=1
    
    #p0 = press1[...,0:-1]
    #p1 = press1[...,1:]  
    #dz = np.diff(data["hfl"],axis=-1)
    #if neAvail: xp = ne.evaluate("-1.*log(p1/p0)/dz")
    #else: xp = -1.*np.log(p1/p0)/dz
        
    #xp[xp==0] = 9999
      
    #if neAvail: data["p_hl"][...,1:-1] = ne.evaluate("-1.*p0/xp*(exp(-xp*dz)-1.)/dz")
    #else: data["p_hl"][...,1:-1] = -1.*p0/xp*(np.exp(-xp*dz)-1.)/dz
    
    #data["p_hl"][...,-1] = -1.*p1[...,-1]/xp[...,-1]*(np.exp(-xp[...,-1]*dz[...,-1])-1.)/dz[...,-1]  #check!
    #print "TO DO: check calculation of p half level!"

    #we can also fill the arrays partly!
    data["temp_lev"] = np.zeros(shape3Dplus) + np.nan
    data["temp_lev"][...,0] = data["t_2m"]
    data["press_lev"] = np.zeros(shape3Dplus) + np.nan
    data["press_lev"][...,0] = data["surface_air_pressure"]
    
    data["relhum"] = meteoSI.q2rh(data["qv"],data["temperature"],data["p"])
    
    data["hydro_q"] = np.zeros(data["qc"].shape + (nHydro,)) + np.nan
    data["hydro_q"][...,0] = data["qc"]
    data["hydro_q"][...,1] = data["qi"]
    data["hydro_q"][...,2] = data["qr"]
    data["hydro_q"][...,3] = data["qs"]
    data["hydro_q"][...,4] = data["qg"]
    
    #translate variables
    varPairs = [["time1h","timestamp"],["latitude","lat"],["longitude","lon"],["fr_land","lfrac"],["u_10m","wind10u"],["v_10m","wind10v"],["hhl","hgt_lev"],["p","press"],["temperature","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["t_s","groundtemp"],["temp_lev","temp_lev"],["press_lev","press_lev"]]
  
  elif kind == "gop_fields_SynSatMic":

    forecastIndex = 0
    variables3D = ["U_10M","V_10M","T_G","surface_air_pressure"]
    variables4D = ["T","P","QV","QC","QI","QI","QR","QS","QG"]
   
    nHydro = 5
    
    conFields = ncToDict(constantFields)
    data = dict()
    
    files = np.sort(glob.glob(fnames))
    if len(files) == 0: raise RuntimeError( "no files found")
    files.sort()
    
    ffOK = 0 #successfull runs
    
    for ff, fname in enumerate(files):
      if verbosity>0: print fname
      try:
        if fname.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fname.split(".")[-1]=="tar":
            if verbosity>3:print "tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz"
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fname)
            if verbosity>2:print "created ", gzFile
          else:
            gzFile = fname
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print "created ", tmpFile
          ncFile = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print "opend ", tmpFile    
        else:  
          ncFile = netCDF4.Dataset(fname,"r")
          if verbosity>1:print "opend ", fname

        #import pdb;pdb.set_trace()
        dataSingle = dict()  
        for var in variables3D:  
          dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,1)
        for var in variables4D:
          dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,2)[...,::-1]#reverse height order  
        timestamp =   ncFile.variables["time"][:]
        ncFile.close()
        if verbosity>1:print "closed nc"
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
          os.system("rm -f "+tmpFile+"*")
        shape3Dplus = tuple(np.array(dataSingle["T"].shape) + np.array([0,0,1]))
        shape3D = dataSingle["T"].shape
        shape2D = shape3D[:2]
        
        dataSingle["timestamp"] = np.zeros(shape2D)
        dataSingle["timestamp"][:] = timestamp
        
        for key in ["lon","lat"]:
          dataSingle[key]  = np.zeros(shape2D)
          dataSingle[key][:] = np.swapaxes(conFields[key],0,1)
        for key in ["FR_LAND","HSURF"]:
          dataSingle[key]  = np.zeros(shape2D)
          dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)
        for key in ["HHL"]:
          dataSingle[key]  = np.zeros(shape3Dplus)
          dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1]
            
    
        if ffOK == 0: #if the first file is broken, checking for ff==0 would fail!
          data = deepcopy(dataSingle)
        else:
          for key in data.keys():
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)
          
        ffOK += 1  
        

      #except IOError:
      except Exception as inst:
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
          os.system("rm -f "+tmpFile+"*")
        print "ERROR:", fname      
        print type(inst)     # the exception instance
        print inst.args      # arguments stored in .args
        print inst
        if debug: import pdb;pdb.set_trace()



    #shapes may have changed!
    shape3Dplus = tuple(np.array(data["T"].shape) + np.array([0,0,1]))
    shape3D = data["T"].shape
    shape2D = shape3D[:2]
    
    #we can also fill the arrays partly!
    data["press_lev"] = np.zeros(shape3Dplus) + np.nan
    data["press_lev"][...,0] = data["surface_air_pressure"]
    
    data["relhum"] = meteoSI.q2rh(data["QV"],data["T"],data["P"])
    
    data["hydro_q"] = np.zeros(data["QC"].shape + (nHydro,)) + np.nan
    data["hydro_q"][...,0] = data["QC"]
    data["hydro_q"][...,1] = data["QI"]
    data["hydro_q"][...,2] = data["QR"]
    data["hydro_q"][...,3] = data["QS"]
    data["hydro_q"][...,4] = data["QG"]
    
    varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],["FR_LAND","lfrac"],["U_10M","wind10u"],["V_10M","wind10v"],["HHL","hgt_lev"],["P","press"],["T","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["T_G","groundtemp"],["press_lev","press_lev"]]    
    
  else:
        raise TypeError("unknown Cosmo DE data type "+kind)

  pamData = dict()
  for cosmoVar,pamVar in varPairs:
    pamData[pamVar] = data[cosmoVar]
    
  pam = pyPamtra.pyPamtra()
  pam.set["pyVerbose"]= verbosity
  
  pam = descriptorFile_cosmo_1mom(pam,kind=df_kind)  
  pam.createProfile(**pamData)
  del data


  return pam
  
def createUsStandardProfile(**kwargs):
  '''
  Function to create clear sky US Standard Atmosphere.
  
  hgt_lev is teh only mandatory variables
  humidity will be set to zero if not provided, all other variables are guessed by "createProfile"
  
  values provided in kwargs will be passed to "createProfile", however, press_lev and temp_lev will overwritte us staandard if provided 
  
  '''    
  
  pamData = _createUsStandardProfile(**kwargs)
  
  pam = pyPamtra.pyPamtra()

  pam.createProfile(**pamData)
  del kwargs
  
  return pam
  
def _createUsStandardProfile(**kwargs):
  '''
  HELPER
  Function to create clear sky US Standard Atmosphere.
  
  hgt_lev is teh only mandatory variables
  humidity will be set to zero if not provided, all other variables are guessed by "createProfile"
  
  values provided in kwargs will be passed to "createProfile", however, press_lev and temp_lev will overwritte us staandard if provided 
  
  '''    
  
  import usStandard #see in tools
  
  assert "hgt_lev" in kwargs.keys() #hgt_lev is mandatory
  
  pamData = dict()
  
  density = np.zeros_like(kwargs["hgt_lev"])
  pamData["press_lev"] = np.zeros_like(kwargs["hgt_lev"])
  pamData["temp_lev"] = np.zeros_like(kwargs["hgt_lev"])
  
  if len(np.shape(kwargs["hgt_lev"]))==1:
    density[:], pamData["press_lev"][:], pamData["temp_lev"][:]  =  usStandard.usStandard(kwargs["hgt_lev"])
  elif  len(np.shape(kwargs["hgt_lev"]))==2:
    for xx in range(np.shape(kwargs["hgt_lev"])[0]):
      density[xx], pamData["press_lev"][xx], pamData["temp_lev"][xx]  =  usStandard.usStandard(kwargs["hgt_lev"][xx])
  elif  len(np.shape(kwargs["hgt_lev"]))==3:
    for xx in range(np.shape(kwargs["hgt_lev"])[0]):
      for yy in range(np.shape(kwargs["hgt_lev"])[1]):
        density[xx,yy], pamData["press_lev"][xx,yy], pamData["temp_lev"][xx,yy]  =  usStandard.usStandard(kwargs["hgt_lev"][xx,yy])
  else: raise IOError("hgt_lev has wrong number of dimensions")
  
  for kk in kwargs.keys():
        pamData[kk] = np.array(kwargs[kk])
  
  if ("relhum_lev" not in kwargs.keys()) and ("q" not in kwargs.keys()):
    pamData["relhum_lev"] = np.zeros_like(kwargs["hgt_lev"])
  
  return pamData
  
  
def descriptorFile_cosmo_1mom(pamObject, kind="default"):
  if kind == "default":
    pamObject.df.addHydrometeor(('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere', 'khvorostyanov01_drops', -99.0))
    pamObject.df.addHydrometeor(('iwc_q', -99.0, -1, 917.0, 130.0, 3.0, 0.684, 2.0, 3, 1, 'mono_cosmo_ice', -99.0, -99.0, -99.0, -99.0, -99.0, -99.0, 'mie-sphere', 'heymsfield10_particles', -99.0))
    pamObject.df.addHydrometeor(('rwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 8000000.0, -99.0, 0.00012, 0.006, 'mie-sphere', 'khvorostyanov01_drops', -99.0))
    pamObject.df.addHydrometeor(('swc_q', -99.0, -1, 200.0, 0.038, 2.0, 0.3971, 1.88, 3, 50, 'exp_field_t', -99.0, -99.0, -99.0, -99.0, 5.1e-11, 0.02, 'mie-sphere', 'heymsfield10_particles', -99.0))
    pamObject.df.addHydrometeor(('gwc_q', -99.0, -1, 400.0, 169.6, 3.1, -99.0, -99.0, 3, 50, 'exp', -99.0, -99.0, 4000000.0, -99.0, 1e-10, 0.01, 'mie-sphere', 'khvorostyanov01_spheres', -99.0))
  else:
    raise ValueError("Do not know kind "+ kind)
  return pamObject
    
    
    
#helper function    
def ncToDict(ncFilePath,keys='all',joinDimension='time',offsetKeys={},ncLib='netCDF4',tmpDir="/tmp/",skipFiles=[]):
  '''
  load keys of netcdf file into dictionary
  by using wildcards in ncFiles, multiple files are possible, they are joint using joinDimension, offsets of e.g. time vector can be corrected by offsetKeys with value to be corrected as key as correction key in value
  gzip compressed netcdf with extension .gz are possible
  '''
  if ncLib == 'netCDF4':
    import netCDF4 as nc
    openNc = nc.Dataset
  elif ncLib == 'netCDF3':
    import netCDF3 as nc
    openNc = nc.Dataset
  elif ncLib == 'Scientific.IO.NetCDF':
    import Scientific.IO.NetCDF as nc
    openNc = nc.NetCDFFile
  else:
    raise ImportError('ncLib must be netCDF4, netCDF3 or Scientific.IO.NetCDF')
  joinDimensionNumber = dict()
  joinedData = dict()
  if type(ncFilePath) == list:
    ncFiles = ncFilePath
  else:
    ncFiles = glob.glob(ncFilePath)
  ncFiles.sort()
  if len(ncFiles) == 0:
    raise IOError('No files found: ' + str(ncFilePath))

  for ncFile in deepcopy(ncFiles):
    if ncFile.split("/")[-1] in skipFiles or ncFile in skipFiles:
      ncFiles.remove(ncFile)
      print "skipping:", ncFile
    
  for nn, ncFile in enumerate(ncFiles):
    tmpFile=False
    if ncFile.split(".")[-1]=="gz":
      tmpFile = True
      gzFile = deepcopy(ncFile)
      ncFile = tmpDir+"/maxLibs_netcdf_"+''.join(randomLib.choice(string.ascii_uppercase + string.digits) for x in range(5))+".tmp.nc"
      print 'uncompressing', nn+1,'of',len(ncFiles), gzFile, ncFile
      returnValue = os.system("zcat "+gzFile+">"+ncFile)
      assert returnValue == 0
    else:
      print 'opening', nn+1,'of',len(ncFiles), ncFile
    try: ncData = openNc(ncFile,'r')
    except: raise RuntimeError("Could not open file: '" + ncFile+"'")
    if nn == 0:
      if keys == 'all':
        keys = ncData.variables.keys()
      #make sure the join dimension is actually present!
      assert joinDimension in ncData.dimensions.keys()
      #get the axis to join the arrays
      for key in keys:
        joinDimensionNumber[key] = -9999
        for dd, dim in enumerate(ncData.variables[key].dimensions):
          if dim == joinDimension:
            joinDimensionNumber[key] = dd
    #get and join the data
    for key in keys:
      #special sausage for nc files with time offset
      if key in offsetKeys.keys():
        data = ncData.variables[key][:] + ncData.variables[offsetKeys[key]].getValue()
      else:
        if ncData.variables[key].shape == ():
          data = ncData.variables[key].getValue()
        else:
          data = ncData.variables[key][:]
      if nn == 0:
        joinedData[key] = data
      elif joinDimensionNumber[key] != -9999:
        joinedData[key] = np.ma.concatenate((joinedData[key],data),axis=joinDimensionNumber[key])
      else:
        #skip since data is already in joinedData
        continue
    ncData.close()
    if (tmpFile and ncFile.split(".")[-2] == "tmp"):
      #os.system("rm -f "+ncFile)
      os.remove(ncFile)
  return joinedData

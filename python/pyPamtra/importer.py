# -*- coding: utf-8 -*-

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


from .core import pyPamtra
from .meteoSI import Rair, q2rh


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

  Cp = 7.*Rair/2.
  RdCp =  Rair/Cp
  p0 = 100000
  data["T"] = data["T"] + 300 # now in K
  data["T"] = data["T"]*(data["P"]/p0)**RdCp # potential temp to temp

  #add surface data
  data["TSK"] = data["TSK"].reshape(np.append(list(np.shape(data["TSK"])),1))
  data["PSFC"] = data["PSFC"].reshape(np.append(list(np.shape(data["PSFC"])),1))

  data["T"] = np.concatenate((data["TSK"],data["T"]),axis=-1)
  data["P"] = np.concatenate((data["PSFC"],data["P"]),axis=-1)


  #translate variables
  varPairs = [["Times","timestamp"],["XLAT","lat"],["XLONG","lon"],["U10","wind10u"],["V10","wind10v"],["PH","hgt_lev"],["P","press_lev"],["T","temp_lev"],["QVAPOR","q"],["QCLOUD","cwc_q"],["QICE","iwc_q"],["QRAIN","rwc_q"],["QSNOW","swc_q"],["QGRAUP","gwc_q"]]

  pamData = dict()

  # surface properties
  pamData['sfc_type'] = np.zeros(data['PSFC'].shape)
  pamData['sfc_type'] = np.around(data['LANDMASK'])
  pamData['sfc_model'] = np.zeros(data['PSFC'].shape)
  pamData['sfc_refl'] = np.chararray(data['PSFC'].shape)
  pamData['sfc_refl'][:] = 'F'

  for wrfVar,pamVar in varPairs:
    pamData[pamVar] = data[wrfVar]

  pam = pyPamtra()
  pam.createProfile(**pamData)
  del data

  return pam


def readCosmoDe1MomDataset(fnames,kind,descriptorFile,forecastIndex = 1,colIndex=0,tmpDir="/tmp/",fnameInTar="",concatenateAxis=1,debug=False,verbosity=0,df_kind="default",constantFields=None,maxLevel=0,subGrid=None):
  '''
  convert cosmo 1-moment dataset with fname of kind to pamtra object

  Parameters
  ----------
  fname : str
    file name, wildCards allowed! can be either nc, nc and gz, or nc.gz files of name fnameInTar within tar file
  kind : str
    kind of Cosmo file, right now only collum and synsat netcdf files are implemented.
  descriptorFile : str or pyPamtra descriptorFile object
    Pamtra descriptor file
  forecastIndex : int, optional
    e.g. 1 means take the forecast being between 3 and 5.75 hours old (default 1).
  colIndex: list or int, optional
    which collum should be taken? list allowed! (default 0)
  tmpDir : str, optional
    temporary directory (default "/tmp/")
  fnameInTar : str, optional
    if nc.gz file in tar file, name of nc.gz file (wildcards allowed!) (default "")
  concatenateAxis : in,t optional
    axis to concatenate the grids (default 1)
  debug: bool optional
    stop and load debugger on exception (default False)
  verbosity : int, optional
    verbosity level (default 0)
  constantFields : str, optional
    fiel with constant fields (default "")s
  subGrid : int array, optional
    array with indices [lon_start,lon_end,lat_start,lat_end] ((1,1) in model corresponds to (0,0) in python!) (default None)

  Returns
  --------

  pam : pyPamtra Object

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

        if maxLevel  == 0: maxLevel = ncFile.variables["hfl"].shape[1]

        #import pdb;pdb.set_trace()
        dataSingle = dict()
        for var in variables1Dx:
          dataSingle[var] = ncFile.variables[var][:]
        for var in variables1Dy:
          dataSingle[var] = ncFile.variables[var][[colIndex]]
        for var in [ "hfl"]:
          dataSingle[var] = ncFile.variables[var][[colIndex],::-1][...,:maxLevel] #reverse height order and cut heights
        for var in [ "hhl"]:
          dataSingle[var] = ncFile.variables[var][[colIndex],::-1][...,:maxLevel+1] #reverse height order and cut heights
        for var in variables3D:
          dataSingle[var] = ncFile.variables[var][[colIndex],forecastIndex,:]
        for var in variables4D:
          dataSingle[var] = np.swapaxes(ncFile.variables[var][[colIndex],:,forecastIndex,:],1,2)[...,::-1][...,:maxLevel]#reverse height order

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

    data["relhum"] = q2rh(data["qv"],data["temperature"],data["p"]) * 100.

    data["hydro_q"] = np.zeros(data["qc"].shape + (nHydro,)) + np.nan
    data["hydro_q"][...,0] = data["qc"]
    data["hydro_q"][...,1] = data["qi"]
    data["hydro_q"][...,2] = data["qr"]
    data["hydro_q"][...,3] = data["qs"]
    data["hydro_q"][...,4] = data["qg"]

    #translate variables
    varPairs = [["time1h","timestamp"],["latitude","lat"],["longitude","lon"],["u_10m","wind10u"],["v_10m","wind10v"],["hhl","hgt_lev"],["p","press"],["temperature","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["t_s","groundtemp"],["temp_lev","temp_lev"],["press_lev","press_lev"]]

  elif kind == "gop_fields_SynSatMic":
    assert constantFields
    forecastIndex = 0
    variables3D = ["U_10M","V_10M","T_G","surface_air_pressure"]
    variables4D = ["T","P","QV","QC","QI","QR","QS","QG"]

    nHydro = 5

    conFields = ncToDict(constantFields)
    data = dict()
    if maxLevel  == 0: maxLevel = conFields["HHL"].shape[1] - 1
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
          if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,1)
          else:
	    dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
        for var in variables4D:
          if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel]#reverse height order
          else:
	    dataSingle[var] = np.swapaxes(ncFile.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]
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
          if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)
          else:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
        for key in ["FR_LAND","HSURF"]:
	  dataSingle[key]  = np.zeros(shape2D)
          if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)
          else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
        for key in ["HHL"]:
	  dataSingle[key]  = np.zeros(shape3Dplus)
          if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1]
          else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]

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

    data["relhum"] = q2rh(data["QV"],data["T"],data["P"]) * 100.

    data["hydro_q"] = np.zeros(data["QC"].shape + (nHydro,)) + np.nan
    data["hydro_q"][...,0] = data["QC"]
    data["hydro_q"][...,1] = data["QI"]
    data["hydro_q"][...,2] = data["QR"]
    data["hydro_q"][...,3] = data["QS"]
    data["hydro_q"][...,4] = data["QG"]

    varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],["U_10M","wind10u"],["V_10M","wind10v"],["HHL","hgt_lev"],["P","press"],["T","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["T_G","groundtemp"],["press_lev","press_lev"]]

  else:
        raise TypeError("unknown Cosmo DE data type "+kind)

  pamData = dict()
  for cosmoVar,pamVar in varPairs:
    pamData[pamVar] = data[cosmoVar]
  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['FR_LAND'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity
  if type(descriptorFile) == str:
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  pam.createProfile(**pamData)
  del data


  return pam

def readCosmoDe2MomDataset(fnamesA,descriptorFile,fnamesN=None,kind='new',forecastIndex = 1,tmpDir="/tmp/",fnameInTar="",debug=False,verbosity=0,df_kind="default",constantFields=None,maxLevel=0,subGrid=None):
  '''
  import cosmo 2-moment dataset

  fnamesA = str , fileNames of atmospheric variables, wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file
  descriptorFile = Pamtra descriptor file
  fnamesN = str , fileNames of number concentrations, wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file (only required for the old style of NARVAL cases)
  kind = str, can be either new (default, new NARVAL calculations) or old (old NARVAL calculations)
  forecastIndex = 1 #take the forecast being between 3 and 5.75 hours old.
  fnameInTar = if nc.gz file in tar file, name of nc.gz file (wildcards allowed!)
  debug: stop and load debugger on exception
  subGrid: array with indices [lon_start,lon_end,lat_start,lat_end] ((1,1) in model corresponds to (0,0) in python!)
  '''
  import netCDF4

  assert constantFields
  forecastIndex = 0
  nHydro = 6
  ffOK = 0 #successfull runs

  conFields = ncToDict(constantFields)
  data = dict()
  if maxLevel  == 0: maxLevel = conFields["HHL"].shape[1] - 1
  filesA = np.sort(glob.glob(fnamesA))
  if len(filesA) == 0: raise RuntimeError( "no files found")
  filesA.sort()

  if kind == 'new':
    variables3D = ["T_G","PS","U_10M","V_10M"]
    variables4D = ["T","P","QV","QC","QI","QR","QS","QG","QH","QNCLOUD","QNICE","QNRAIN","QNSNOW","QNGRAUPEL","QNHAIL"]


    for ff, fnameA in enumerate(filesA):
      if verbosity>0: print fnameA
      try:
	if fnameA.split(".")[-1]!="nc":
	  tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
	  if fnameA.split(".")[-1]=="tar":
	    if verbosity>3:print "tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz"
	    os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
	    gzFile = tmpFile+".gz"
	    if os.stat(gzFile).st_size == 0:
	      os.system("rm -f "+tmpFile+"*")
	      raise IOError("fnameInTar not found in "+fnameA)
	    if verbosity>2:print "created ", gzFile
	  else:
	    gzFile = fnameA
	  os.system("zcat "+gzFile+">"+tmpFile)
	  if verbosity>1:print "created ", tmpFile
	  ncFileA = netCDF4.Dataset(tmpFile,"r")
	  if verbosity>1:print "opend ", tmpFile
	else:
	  ncFileA = netCDF4.Dataset(fnameA,"r")
	  if verbosity>1:print "opend ", fnameA

	#import pdb;pdb.set_trace()
	dataSingle = dict()
	for var in variables3D:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,1)
	  else:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for var in variables4D:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel]#reverse height order
	  else:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]
	timestamp = ncFileA.variables["time"][:]
	ncFileA.close()
	if verbosity>1:print "closed nc"
	if fnameA.split(".")[-1]!="nc":
	  if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
	  os.system("rm -f "+tmpFile+"*")
	shape3Dplus = tuple(np.array(dataSingle["T"].shape) + np.array([0,0,1]))
	shape3D = dataSingle["T"].shape
	shape2D = shape3D[:2]


	dataSingle["timestamp"] = np.zeros(shape2D)
	dataSingle["timestamp"][:] = timestamp

	for key in ["lon","lat"]:
	  dataSingle[key]  = np.zeros(shape2D)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for key in ["FR_LAND","HSURF"]:
	  dataSingle[key]  = np.zeros(shape2D)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for key in ["HHL"]:
	  dataSingle[key]  = np.zeros(shape3Dplus)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1]
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]

	if ffOK == 0: #if the first file is broken, checking for ff==0 would fail!
	  data = deepcopy(dataSingle)
	else:
	  for key in data.keys():
	    data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

	ffOK += 1


      #except IOError:
      except Exception as inst:
	if fnameA.split(".")[-1]!="nc":
	  if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
	  os.system("rm -f "+tmpFile+"*")
	print "ERROR:", fnameA
	print type(inst)     # the exception instance
	print inst.args      # arguments stored in .args
	print inst
	if debug: import pdb;pdb.set_trace()

    data["hydro_n"] = np.zeros(data["QNCLOUD"].shape + (nHydro,)) + np.nan
    data["hydro_n"][...,0] = data["QNCLOUD"]
    data["hydro_n"][...,1] = data["QNICE"]
    data["hydro_n"][...,2] = data["QNRAIN"]
    data["hydro_n"][...,3] = data["QNSNOW"]
    data["hydro_n"][...,4] = data["QNGRAUPEL"]
    data["hydro_n"][...,5] = data["QNHAIL"]

  elif kind == 'old':

    filesN = np.sort(glob.glob(fnamesN))
    if len(filesN) == 0: raise RuntimeError( "no files found")
    filesN.sort()

    if len(filesA) != len(filesN): raise RuntimeError("number of atmospheric files does not equal number of concentration files!")

    variables3D = ["T_G","PS"]
    variables3Dplus1 = ["U_10M","V_10M"]
    variables4D = ["T","P","QV","QC","QI","QR","QS","QG","QH"]
    variables4DN = ['QNC','QNI','QNR','QNS','QNG','QNH']


    for ff, fnameA in enumerate(filesA):
      fnameN = filesN[ff]
      if verbosity>0: print fnameA, fnameN
      try:
	if fnameA.split(".")[-1]!="nc":
	  tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
	  if fnameA.split(".")[-1]=="tar":
	    if verbosity>3:print "tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz"
	    os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
	    gzFile = tmpFile+".gz"
	    if os.stat(gzFile).st_size == 0:
	      os.system("rm -f "+tmpFile+"*")
	      raise IOError("fnameInTar not found in "+fnameA)
	    if verbosity>2:print "created ", gzFile
	  else:
	    gzFile = fnameA
	  os.system("zcat "+gzFile+">"+tmpFile)
	  if verbosity>1:print "created ", tmpFile
	  ncFileA = netCDF4.Dataset(tmpFile,"r")
	  if verbosity>1:print "opend ", tmpFile
	else:
	  ncFileA = netCDF4.Dataset(fnameA,"r")
	  if verbosity>1:print "opend ", fnameA

	# number density files
	if fnameN.split(".")[-1]!="nc":
	  tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
	  if fnameN.split(".")[-1]=="tar":
	    if verbosity>3:print "tar -O --wildcards -x "+fnameInTar+" -f "+fnameN+">"+tmpFile+".gz"
	    os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameN+">"+tmpFile+".gz")
	    gzFile = tmpFile+".gz"
	    if os.stat(gzFile).st_size == 0:
	      os.system("rm -f "+tmpFile+"*")
	      raise IOError("fnameInTar not found in "+fnameN)
	    if verbosity>2:print "created ", gzFile
	  else:
	    gzFile = fnameN
	  os.system("zcat "+gzFile+">"+tmpFile)
	  if verbosity>1:print "created ", tmpFile
	  ncFileN = netCDF4.Dataset(tmpFile,"r")
	  if verbosity>1:print "opend ", tmpFile
	else:
	  ncFileN = netCDF4.Dataset(fnameN,"r")
	  if verbosity>1:print "opend ", fnameN

	#import pdb;pdb.set_trace()
	dataSingle = dict()
	for var in variables3D:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,1)
	  else:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for var in variables3Dplus1:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex,0],0,1)
	  else:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex,0],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for var in variables4D:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel]#reverse height order
	  else:
	    dataSingle[var] = np.swapaxes(ncFileA.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]
	timestamp = ncFileA.variables["time"][:]
	ncFileA.close()
	if verbosity>1:print "closed nc"
	if fnameA.split(".")[-1]!="nc":
	  if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
	  os.system("rm -f "+tmpFile+"*")
	for var in variables4DN:
	  if subGrid == None:
	    dataSingle[var] = np.swapaxes(ncFileN.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel]#reverse height order
	  else:
	    dataSingle[var] = np.swapaxes(ncFileN.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]
	ncFileN.close()
	if verbosity>1:print "closed nc"
	if fnameN.split(".")[-1]!="nc":
	  if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
	  os.system("rm -f "+tmpFile+"*")
	shape3Dplus = tuple(np.array(dataSingle["T"].shape) + np.array([0,0,1]))
	shape3D = dataSingle["T"].shape
	shape2D = shape3D[:2]


	dataSingle["timestamp"] = np.zeros(shape2D)
	dataSingle["timestamp"][:] = timestamp

	for key in ["lon","lat"]:
	  dataSingle[key]  = np.zeros(shape2D)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for key in ["FR_LAND","HSURF"]:
	  dataSingle[key]  = np.zeros(shape2D)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)[subGrid[0]:subGrid[1],subGrid[2]:subGrid[3]]
	for key in ["HHL"]:
	  dataSingle[key]  = np.zeros(shape3Dplus)
	  if subGrid == None:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1]
	  else:
	    dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]

	if ffOK == 0: #if the first file is broken, checking for ff==0 would fail!
	  data = deepcopy(dataSingle)
	else:
	  for key in data.keys():
	    data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

	ffOK += 1


      #except IOError:
      except Exception as inst:
	if fnameA.split(".")[-1]!="nc":
	  if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
	  os.system("rm -f "+tmpFile+"*")
	print "ERROR:", fnameA
	print type(inst)     # the exception instance
	print inst.args      # arguments stored in .args
	print inst
	if debug: import pdb;pdb.set_trace()

    data["hydro_n"] = np.zeros(data["QNC"].shape + (nHydro,)) + np.nan
    data["hydro_n"][...,0] = data["QNC"]
    data["hydro_n"][...,1] = data["QNI"]
    data["hydro_n"][...,2] = data["QNR"]
    data["hydro_n"][...,3] = data["QNS"]
    data["hydro_n"][...,4] = data["QNG"]
    data["hydro_n"][...,5] = data["QNH"]

  #shapes may have changed!
  shape3Dplus = tuple(np.array(data["T"].shape) + np.array([0,0,1]))
  shape3D = data["T"].shape
  shape2D = shape3D[:2]

  #we can also fill the arrays partly!
  data["press_lev"] = np.zeros(shape3Dplus) + np.nan
  data["press_lev"][...,0] = data["PS"]

  data["relhum"] = q2rh(data["QV"],data["T"],data["P"]) * 100.

  data["hydro_q"] = np.zeros(data["QC"].shape + (nHydro,)) + np.nan
  data["hydro_q"][...,0] = data["QC"]
  data["hydro_q"][...,1] = data["QI"]
  data["hydro_q"][...,2] = data["QR"]
  data["hydro_q"][...,3] = data["QS"]
  data["hydro_q"][...,4] = data["QG"]
  data["hydro_q"][...,5] = data["QH"]

  varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],["U_10M","wind10u"],["V_10M","wind10v"],["HHL","hgt_lev"],["P","press"],["T","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["hydro_n","hydro_n"],["T_G","groundtemp"],["press_lev","press_lev"]]

  pamData = dict()
  for cosmoVar,pamVar in varPairs:
    pamData[pamVar] = data[cosmoVar]

  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['FR_LAND'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity
  if type(descriptorFile) == str:
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  pam.createProfile(**pamData)

  return pam

def readCosmoDe2MomDatasetOnFlightTrack(fnameA,descriptorFile,tmpDir="/tmp/",debug=False,verbosity=0,df_kind="default",maxLevel=0):
  '''
  import cosmo 2-moment dataset extracted on a HALO flight track^

  fnamesA = str , fileNames of atmospheric variables, wildCards allowed! can be either nc file, nc,gz file or nc.gz file
  descriptorFile = Pamtra descriptor file
  debug: stop and load debugger on exception
  '''
  import netCDF4

  nHydro = 6

#  filesA = np.sort(glob.glob(fnamesA))
#  if len(filesA) == 0: raise RuntimeError( "no files found")
#  filesA.sort()

  variables1D = ["T_G","PS","U_10M","V_10M","FR_LAND"]
  variables2D = ["T","P","QV","QC","QI","QR","QS","QG","QH","QNCLOUD","QNICE","QNRAIN","QNSNOW","QNGRAUPEL","QNHAIL"]

  if verbosity>0: print fnameA
  try:
    if fnameA.split(".")[-1]!="nc":
      tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
      gzFile = fnameA
      os.system("zcat "+gzFile+">"+tmpFile)
      if verbosity>1:print "created ", tmpFile
      data = ncToDict(tmpFile)
      if verbosity>1:print "read ", tmpFile
    else:
      data = ncToDict(fnameA)
      if verbosity>1:print "read ", fnameA

    if debug: import pdb;pdb.set_trace()
    #dataSingle = dict()
    #for var in variables1D:
      #data[var] = np.swapaxes(ncFileA.variables[var][0],0,1)
    #for var in variables2D:
      #data[var] = np.swapaxes(ncFileA.variables[var][0],0,2)[...,::-1][...,:maxLevel]#reverse height order
    if verbosity>1:print "closed nc"
    if fnameA.split(".")[-1]!="nc":
      if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      os.system("rm -f "+tmpFile+"*")

    #for key in ["cosmoTime","cosmoLon","cosmoLat"]:
      #dataSingle[key]  = np.zeros(shape2D)
      #dataSingle[key][:] = np.swapaxes(conFields[key],0,1)
    #for key in ["FR_LAND","HSURF"]:
      #dataSingle[key]  = np.zeros(shape2D)
      #dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,1)
    #for key in ["HHL"]:
      #dataSingle[key]  = np.zeros(shape3Dplus)
      #dataSingle[key][:] = np.swapaxes(conFields[key][forecastIndex],0,2)[...,::-1][...,:maxLevel+1]
    #data = deepcopy(dataSingle)
  except Exception as inst:
    if fnameA.split(".")[-1]!="nc":
      if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      os.system("rm -f "+tmpFile+"*")
    print "ERROR:", fnameA
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()

  shape3D = (data["T"].shape[0], 1, data["T"].shape[1])
  shape3Dplus = (shape3D[0], shape3D[1], shape3D[2] + 1)
  shape2D = shape3D[:2]

  if maxLevel  == 0: maxLevel = data["HHL"].shape[1] - 1

  pamData = dict()

  pamData["hydro_n"] = np.zeros(data["QNCLOUD"].shape + (nHydro,)) + np.nan
  pamData["hydro_n"][...,0] = data["QNCLOUD"][:,::-1]
  pamData["hydro_n"][...,1] = data["QNICE"][:,::-1]
  pamData["hydro_n"][...,2] = data["QNRAIN"][:,::-1]
  pamData["hydro_n"][...,3] = data["QNSNOW"][:,::-1]
  pamData["hydro_n"][...,4] = data["QNGRAUPEL"][:,::-1]
  pamData["hydro_n"][...,5] = data["QNHAIL"][:,::-1]

  pamData["hydro_q"] = np.zeros(data["QC"].shape + (nHydro,)) + np.nan
  pamData["hydro_q"][...,0] = data["QC"][:,::-1]
  pamData["hydro_q"][...,1] = data["QI"][:,::-1]
  pamData["hydro_q"][...,2] = data["QR"][:,::-1]
  pamData["hydro_q"][...,3] = data["QS"][:,::-1]
  pamData["hydro_q"][...,4] = data["QG"][:,::-1]
  pamData["hydro_q"][...,5] = data["QH"][:,::-1]

  #we can also fill the arrays partly!
  pamData["press_lev"] = np.zeros(shape3Dplus) + np.nan
  pamData["press_lev"][:,0,0] = data["PS"]

  pamData["relhum"] = np.expand_dims((q2rh(data["QV"],data["T"],data["P"]) * 100.)[:,::-1],axis=1)

#  import pdb;pdb.set_trace()
  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['FR_LAND'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  varPairs = [["cosmoTime","timestamp"],["cosmoLat","lat"],["cosmoLon","lon"],["sfc_type","sfc_type"],["sfc_model","sfc_model"],["U_10M","wind10u"],["V_10M","wind10v"],["T_G","groundtemp"]]

  for cosmoVar,pamVar in varPairs:
    pamData[pamVar] = np.expand_dims(data[cosmoVar],axis=1)

  varPairs = [["HHL","hgt_lev"],["P","press"],["T","temp"]]

  for cosmoVar,pamVar in varPairs:
    pamData[pamVar] = np.expand_dims(data[cosmoVar][:,::-1],axis=1)

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity
  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  pam.createProfile(**pamData)

  return pam

def readCosmoReAn2km(constantFields,fname,descriptorFile,forecastIndex = 1,tmpDir="/tmp/",fnameInTar="",debug=False,verbosity=0,df_kind="default",maxLevel=51,subGrid=None):

  import pygrib
  import os
  import time
  import datetime

  variables2DC = ['rlat','rlon','fr_land','soiltyp']
  variables3DC = ['hhl']

  variables2D = ['t_g','sp','10u','10v']
  variables3D = ['t','pp','q']
  variables4D = ['qc','qi','qr','qs','qg']

  nhydro = len(variables4D)

  data = dict()

  try:
    grbsC = pygrib.open(constantFields)

    # determine dimensions
    nLon = grbsC.select(shortName='rlon')[0]['Ni']
    nLat = grbsC.select(shortName='rlon')[0]['Nj']
    nLev = len(grbsC(shortName='hhl'))

    shape2D = (nLat,nLon)
    shape3D = (nLat,nLon,nLev-1)
    shape3Dplus = (nLat,nLon,nLev)
    shape4D = (nLat,nLon,nLev-1,nhydro)

    grbsC.seek(0)

    data['hhl'] = np.zeros(shape3Dplus) + np.nan

    for var in variables3DC:
      if verbosity>1: print var
      selected_grbs = grbsC(shortName=var)
      for i in range(nLev-maxLevel,nLev):
	data[var][...,nLev-1-i] = selected_grbs[i].values

    selected_grbs = grbsC(shortName=variables2DC)
    for var in selected_grbs:
      if verbosity>1: print var.shortName
      data[var.shortName] = var.values

    grbsC.close()

    grbs = pygrib.open(fname)

    selected_grbs = grbs(shortName=variables2D)
    for var in selected_grbs:
      if verbosity>1: print var.shortName
      data[var.shortName] = var.values

    for var in variables3D:
      if verbosity>1: print var
      data[var] = np.zeros(shape3D) + np.nan
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
	data[selected_grbs[i].shortName][...,nLev-2-i] = selected_grbs[i].values

    data['hydro_q'] = np.zeros(shape4D) + np.nan

    for j,var in enumerate(variables4D):
      if verbosity>1: print var
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
	data['hydro_q'][...,nLev-2-i,j] = selected_grbs[i].values

    grbs.close()
  except Exception as inst:
    if fname.split(".")[-1]!="nc":
      if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      os.system("rm -f "+tmpFile+"*")
    print "ERROR:", fname
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()

  def calc_p0(height):

    # parameters taken from COSMO documentation
    psl = 100000 #Pa
    tsl = 288.15 #K
    beta = 42    #K
    g = 9.80665  #m s-2
    rd = 287.05  #J kg-1 K-1

    # calculate reference pressure (equation 4 in documentation for beta != 0)
    p0 = psl * np.exp(-tsl/beta*(1-np.sqrt(1-2*beta*g*height/(rd*tsl**2))))

    return p0

  # some conversions and filling of needed variables

  data['hfl'] = (data['hhl'][...,1:]+data['hhl'][...,:-1])/2.
  pref = calc_p0(data['hfl'])
  data['press'] = pref + data['pp']*100.
  data['relhum'] = q2rh(data['q']/(1+data['q']),data['t'],data['press'])*100.

  data['timestamp'] = np.zeros(shape2D)
  data['timestamp'][:,:] = time.mktime(datetime.datetime.strptime(os.path.basename(fname), "laf%Y%m%d%H%M0000").timetuple())

  pam = pyPamtra()
  pam.df.readFile(descriptorFile)
  varPairs = [["timestamp","timestamp"],["rlat","lat"],["rlon","lon"],["10u","wind10u"],["10v","wind10v"],
	      ["press","press"],["t","temp"],["relhum","relhum"],["t_g","groundtemp"],['hhl','hgt_lev'],['hydro_q','hydro_q']]

  pamData = dict()

  for dataVar,pamVar in varPairs:
    pamData[pamVar] = data[dataVar]
  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['fr_land'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  del data

  pam.createProfile(**pamData)

  return pam

def readCosmoReAn6km(constantFields,fname,descriptorFile,forecastIndex = 1,tmpDir="/tmp/",fnameInTar="",debug=False,verbosity=0,df_kind="default",maxLevel=41,subGrid=None):

  import pygrib
  import os
  import time
  import datetime

  variables2DC = ['RLAT','RLON','lsm']
  variables3DC = ['HHL']

  variables2D = ['T_G','sp','10u','10v']
  variables3D = ['t','PP','q']
  variables4D = ['QC','QI','QR','QS']

  nhydro = len(variables4D)

  data = dict()

  try:
    # determine dimensions
    constantFile = constantFields+'cosmo_cordex_lon'
    grbsC = pygrib.open(constantFile)
    nLon = grbsC.select(shortName='RLON')[0]['Ni']
    nLat = grbsC.select(shortName='RLON')[0]['Nj']
    data['rlon'] = grbsC.select(shortName='RLON')[0].values
    grbsC.close()

    constantFile = constantFields+'cosmo_cordex_lat'
    grbsC = pygrib.open(constantFile)
    data['rlat'] = grbsC.select(shortName='RLAT')[0].values
    grbsC.close()

    constantFile = constantFields+'fr_land_cordex'
    grbsC = pygrib.open(constantFile)
    data['fr_land'] = grbsC.select(shortName='lsm')[0].values
    grbsC.close()

    constantFile = constantFields+'cosmo_cordex_HHL'
    grbsC = pygrib.open(constantFile)
    selected_grbs = grbsC(shortName='HHL')
    nLev = len(selected_grbs)
    shape2D = (nLat,nLon)
    shape3D = (nLat,nLon,nLev-1)
    shape3Dplus = (nLat,nLon,nLev)
    shape4D = (nLat,nLon,nLev-1,nhydro)

    data['hhl'] = np.zeros(shape3Dplus) + np.nan
    for i in range(nLev-maxLevel,nLev):
	data['hhl'][...,nLev-1-i] = selected_grbs[i].values
    grbsC.close()

    grbs = pygrib.open(fname)


    for var in variables3D:
      if verbosity>1: print var
      data[var] = np.zeros(shape3D) + np.nan
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
	data[selected_grbs[i].shortName][...,nLev-2-i] = selected_grbs[i].values

    data['hydro_q'] = np.zeros(shape4D) + np.nan

    for j,var in enumerate(variables4D):
      if verbosity>1: print var
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
	data['hydro_q'][...,nLev-2-i,j] = selected_grbs[i].values


    grbs.close()

    grbs = pygrib.open(fname+'00')

    selected_grbs = grbs(shortName=variables2D)
    for var in selected_grbs:
      if verbosity>1: print var.shortName
      data[var.shortName] = var.values


    grbs.close()

  except Exception as inst:
    #if fname.split(".")[-1]!="nc":
      #if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      #os.system("rm -f "+tmpFile+"*")
    print "ERROR:", fname
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()

  def calc_p0(height):

    # parameters taken from COSMO documentation
    psl = 100000 #Pa
    tsl = 288.15 #K
    beta = 42    #K
    g = 9.80665  #m s-2
    rd = 287.05  #J kg-1 K-1

    # calculate reference pressure (equation 4 in documentation for beta != 0)
    p0 = psl * np.exp(-tsl/beta*(1-np.sqrt(1-2*beta*g*height/(rd*tsl**2))))

    return p0

  # some conversions and filling of needed variables

  data['hfl'] = (data['hhl'][...,1:]+data['hhl'][...,:-1])/2.
  pref = calc_p0(data['hfl'])
  data['press'] = pref + data['PP']*100.
  data['relhum'] = q2rh(data['q']/(1+data['q']),data['t'],data['press'])*100.

  data['timestamp'] = np.zeros(shape2D)
  data['timestamp'][:,:] = time.mktime(datetime.datetime.strptime(os.path.basename(fname), "laf%Y%m%d%H%M").timetuple())

  pam = pyPamtra()
  pam.df.readFile(descriptorFile)
  varPairs = [["timestamp","timestamp"],["rlat","lat"],["rlon","lon"],["10u","wind10u"],["10v","wind10v"],
	      ["press","press"],["t","temp"],["relhum","relhum"],["T_G","groundtemp"],['hhl','hgt_lev'],['hydro_q','hydro_q']]

  pamData = dict()


  for dataVar,pamVar in varPairs:
    pamData[pamVar] = data[dataVar]

  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['fr_land'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  del data

  pam.createProfile(**pamData)

  return pam

__ICDN_regridding_remarks = """
Remarks
=======
This importer is adjusted to read regridded ICON-LEM simulations as they where run by Daniel Klocke within the HErZ-NARVALII framework.
The Output consists of a *_fg_* and a *_cloud_* netcdf file.
From the *_fg_* file, which contains most atmospheric variables, we need data on grid 2.
From the *_cloud_* file, which contains qg, we need data on grid 1.
Also there is the extpar_narval_nestTropAtl_R1250m.nc file wich contains time constant parameter.

Regridding
----------
Use cdo to regrid the unstructured ICON data on a latlon grid.
Nearest neighbor interpolation is used for regridding.
The interpolation is defined in a Grid file.
Example (filename "test_grid"):

```
gridtype   = lonlat
gridsize   = 8000
xname      = lon
xlongname  = longitude
xunits     = degrees_east
yname      = lat
ylongname  = latitude
yunits     = degree_north
xsize      = 100
ysize      = 80
xfirst     = -59.0
xinc       = 0.05
yfirst     = 10.0
yinc       = 0.05
```

Convert unstructured data to latlon in one cdo command:
    cdo remapnn,test_grid -selgrid,2 -setgrid,narval_nestTropAtl_R1250m.nc /dei4_NARVALII_2016081700_fg_DOM02_ML_0016.nc dei4_NARVALII_2016081700_fg_DOM02_ML_0016_test_grid.nc

Here, "narval_nestTropAtl_R1250m.nc" is the description of the unstructured data.
Current location:
/data/hamp/models/ICON/HErZ-NARVALII/GRIDS/narval_nestTropAtl_R1250m.nc

The corresponding command for _cloud_ data is:
    cdo remapnn,test_grid -selgrid,1 -setgrid,narval_nestTropAtl_R1250m.nc /dei4_NARVALII_2016081700_cloud_DOM02_ML_0016.nc dei4_NARVALII_2016081700_cloud_DOM02_ML_0016_test_grid.nc

Optimized regridding
--------------------
To speedup the regridding of multiple ICON output one can save the interpolation weights in an intermediate netcdf file.
All ICON output must be given on the same unstructured grid.

Generate the weights:
    cdo gennn,test_grid -selgrid,2 -setgrid,narval_nestTropAtl_R1250m.nc dei4_NARVALII_2016081700_fg_DOM02_ML_0016.nc test_grid_weights.nc

The same weights can be used for _fg_ and _cloud_ (unconfirmed).
Apply weights:
    cdo remap,test_grid,test_grid_weights.nc -selgrid,2 -setgrid,narval_nestTropAtl_R1250m.nc dei4_NARVALII_2016081700_fg_DOM02_ML_0016.nc dei4_NARVALII_2016081700_fg_DOM02_ML_0016_test_grid.nc
    cdo remap,test_grid,test_grid_weights.nc -selgrid,1 -setgrid,narval_nestTropAtl_R1250m.nc dei4_NARVALII_2016081700_cloud_DOM02_ML_0016.nc dei4_NARVALII_2016081700_cloud_DOM02_ML_0016_test_grid.nc
    cdo remapnn,test_grid,test_grid_weights.nc -setgrid,narval_nestTropAtl_R1250m extpar_narval_nestTropAtl_R1250m.nc extpar_narval_nestTropAtl_R1250m_test_grid.nc
"""

def readIconLem1MomDataset(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0):
  '''
  import ICON LEM 1-moment dataset
  It is Icon with cosmo physics

  fname_fg = str , fileName of atmospheric variables ("_fg_" file, no wildCards allowed! must be nc file.
              A corresponding "_cloud_" has do exist!
  descriptorFile = Pamtra descriptor file
  debug: stop and load debugger on exception
  constantFields: str, file name of nc containing constant fields (ie. FR_LAND)
  '''
  import netCDF4

  assert constantFields
  forecastIndex = 0 # time step in forecast
  nHydro = 5

  data = dict()

  variables2D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
  variables3D = ["t_g","pres_sfc"]
  variables4D_10m = ["u_10m","v_10m"]
  variables4D = ["temp","pres","qv","qc","qi","qr","qs"]
  variables4D_cloud = ["qg"]

  if verbosity>0: print fname_fg

  if not fname_fg.endswith('.nc'):
    raise IOError("fname_fg has to be .nc.", fname_fg)
  if not '_fg_' in os.path.basename(fname_fg):
    raise IOError("fname_fg has to contain '_fg_'.", fname_fg)

  dataSingle = dict()
  try:
    ncFile_const = netCDF4.Dataset(constantFields, "r")
    if verbosity > 1: print "opened ", constantFields

    for var in variables2D_const:
      # nc dimensions: lat, lon; target dimensions: lon, lat
      assert ncFile_const.variables[var].dimensions == ('lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_const.variables[var],0,1)

    ncFile_const.close()
    if verbosity > 1: print "closed const nc"

    ncFile_fg = netCDF4.Dataset(fname_fg, "r")
    if verbosity > 1: print "opened ", fname_fg

    fname_cloud = '_cloud_'.join(fname_fg.rsplit('_fg_',1)) # right-replace _fg_ with _cloud
    ncFile_cloud = netCDF4.Dataset(fname_cloud, "r")
    if verbosity > 1: print "opened ", fname_cloud

    if maxLevel == 0:
      assert ncFile_fg.variables["z_ifc"].dimensions == ('lat', 'height_3', 'lon')
      assert ncFile_fg.dimensions['height'].size + 1 == ncFile_fg.dimensions['height_3'].size
      maxLevel = ncFile_fg.variables["z_ifc"].shape[1] - 1

    for var in variables3D:
      # nc dimensions: time, lat, lon; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_fg.variables[var][forecastIndex],0,1)

    for var in variables4D_10m:
      # nc dimensions: time, lat, height_5 lon; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'height_5', 'lon')
      assert ncFile_fg.dimensions['height_5'].size == 1
      dataSingle[var] = np.swapaxes(ncFile_fg.variables[var][forecastIndex, :, 0, :],0,1)

    for var in variables4D:
      # nc dimensions: time, lat, height, lon; target dimensions: lon, lat, level
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'height', 'lon')
      dataSingle[var] = np.transpose(ncFile_fg.variables[var][forecastIndex], (2, 0, 1))[...,::-1][...,:maxLevel]#reverse height order

    for var in variables4D_cloud:
      # nc dimensions: time, lat, height, lon; target dimensions: lon, lat, level
      assert ncFile_cloud.variables[var].dimensions == ('time', 'lat', 'height', 'lon')
      dataSingle[var] = np.transpose(ncFile_cloud.variables[var][forecastIndex], (2, 0, 1))[...,::-1][...,:maxLevel]#reverse height order

    shape3D = dataSingle["temp"].shape
    shape3Dplus = (shape3D[0], shape3D[1], shape3D[2] + 1)
    shape2D = shape3D[:2]

    date_times = netCDF4.num2date(ncFile_fg.variables["time"][:], ncFile_fg.variables["time"].units) # convert from any given reference time
    timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp

    dataSingle["timestamp"] = np.zeros(shape2D)
    dataSingle["timestamp"][:] = timestamp

    # 1D variables, right angled lat-lon grid
    #lat, lon = np.meshgrid(ncFile_fg.variables['lat'][:], ncFile_fg.variables['lon'][:])
    #assert lat.shape == shape2D
    #dataSingle['lat'] = lat
    #dataSingle['lon'] = lon # use exact lon_2 and lat_2 from constantFields.(the coordinates found by nearest neighbor interpolation)

    for key in ["z_ifc"]:
      # nc dimensions: lat, height, lon; target dimensions: lon, lat, level
      assert ncFile_fg.variables[key].dimensions == ('lat', 'height_3', 'lon')
      dataSingle[key] = np.transpose(ncFile_fg[key],(2, 0, 1))[...,::-1][...,:maxLevel+1]#reverse height order
      assert dataSingle[key].shape == shape3Dplus

    ncFile_fg.close()
    if verbosity > 1: print "closed fg nc"
    ncFile_cloud.close()
    if verbosity > 1: print "closed cloud nc"

    data = dataSingle

  #except IOError:
  except Exception as inst:
    print "ERROR:", fname_fg
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()
    raise


  #shapes may have changed!
  shape3Dplus = tuple(np.array(data["temp"].shape) + np.array([0,0,1]))
  shape3D = data["temp"].shape
  shape2D = shape3D[:2]

  #we can also fill the arrays partly!
  data["press_lev"] = np.zeros(shape3Dplus) + np.nan
  data["press_lev"][...,0] = data["pres_sfc"]

  data["relhum"] = q2rh(data["qv"],data["temp"],data["pres"]) * 100.

  data["hydro_q"] = np.zeros(data["qc"].shape + (nHydro,)) + np.nan
  data["hydro_q"][...,0] = data["qc"]
  data["hydro_q"][...,1] = data["qi"]
  data["hydro_q"][...,2] = data["qr"]
  data["hydro_q"][...,3] = data["qs"]
  data["hydro_q"][...,4] = data["qg"]

  varPairs = [["timestamp","timestamp"],["lat_2","lat"],["lon_2","lon"],
    ["u_10m","wind10u"],["v_10m","wind10v"],
    ["z_ifc","hgt_lev"],["pres","press"],["temp","temp"],["relhum","relhum"],["hydro_q","hydro_q"],
    ["t_g","groundtemp"],["press_lev","press_lev"]]

  pamData = dict()
  for iconVar, pamVar in varPairs:
    pamData[pamVar] = data[iconVar]

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity
  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  # surface properties
  pamData['sfc_type'] = np.around(data['FR_LAND'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  pam.createProfile(**pamData)

  return pam

# Add regridding remarks
readIconLem1MomDataset.__doc__ += __ICDN_regridding_remarks


def readIconLem2MomDataset(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0):
  '''
  import ICON LEM 2-moment dataset
  It is Icon with cosmo physics

  fname_fg = str , fileName of atmospheric variables ("_fg_" file, no wildCards allowed! must be nc file.
  descriptorFile = Pamtra descriptor file
  debug: stop and load debugger on exception
  constantFields: str, file name of nc containing constant fields (ie. FR_LAND)
  '''
  import netCDF4

  assert constantFields
  forecastIndex = 0 # time step in forecast
  nHydro = 6

  data = dict()

  variables2D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
  variables3D = ["t_g","pres_sfc"]
  variables4D_10m = ["u_10m","v_10m"]
  variables4D = ["temp","pres","qv","qc","qi","qr","qs","qg","qh","qnc","qni","qnr","qns","qng","qnh"]

  if verbosity>0: print fname_fg

  if not fname_fg.endswith('.nc'):
    raise IOError("fname_fg has to be .nc.", fname_fg)
  if not '_fg_' in os.path.basename(fname_fg):
    raise IOError("fname_fg has to contain '_fg_'.", fname_fg)

  dataSingle = dict()
  try:
    ncFile_const = netCDF4.Dataset(constantFields, "r")
    if verbosity > 1: print "opened ", constantFields

    for var in variables2D_const:
      # nc dimensions: lat, lon; target dimensions: lon, lat
      assert ncFile_const.variables[var].dimensions == ('lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_const.variables[var],0,1)

    ncFile_const.close()
    if verbosity > 1: print "closed const nc"

    ncFile_fg = netCDF4.Dataset(fname_fg, "r")
    if verbosity > 1: print "opened ", fname_fg


    if maxLevel == 0:
      assert ncFile_fg.variables["z_ifc"].dimensions == ('lat', 'height_3', 'lon')
      assert ncFile_fg.dimensions['height'].size + 1 == ncFile_fg.dimensions['height_3'].size
      maxLevel = ncFile_fg.variables["z_ifc"].shape[1] - 1

    for var in variables3D:
      # nc dimensions: time, lat, lon; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_fg.variables[var][forecastIndex],0,1)

    for var in variables4D_10m:
      # nc dimensions: time, lat, height_5 lon; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'height_5', 'lon')
      assert ncFile_fg.dimensions['height_5'].size == 1
      dataSingle[var] = np.swapaxes(ncFile_fg.variables[var][forecastIndex, :, 0, :],0,1)

    for var in variables4D:
      # nc dimensions: time, lat, height, lon; target dimensions: lon, lat, level
      assert ncFile_fg.variables[var].dimensions == ('time', 'lat', 'height', 'lon')
      dataSingle[var] = np.transpose(ncFile_fg.variables[var][forecastIndex], (2, 0, 1))[...,::-1][...,:maxLevel]#reverse height order

    shape3D = dataSingle["temp"].shape
    shape3Dplus = (shape3D[0], shape3D[1], shape3D[2] + 1)
    shape2D = shape3D[:2]

    date_times = netCDF4.num2date(ncFile_fg.variables["time"][:], ncFile_fg.variables["time"].units) # convert from any given reference time
    timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp

    dataSingle["timestamp"] = np.zeros(shape2D)
    dataSingle["timestamp"][:] = timestamp

    # 1D variables, right angled lat-lon grid
    #lat, lon = np.meshgrid(ncFile_fg.variables['lat'][:], ncFile_fg.variables['lon'][:])
    #assert lat.shape == shape2D
    #dataSingle['lat'] = lat
    #dataSingle['lon'] = lon # use exact lon_2 and lat_2 from constantFields.(the coordinates found by nearest neighbor interpolation)

    for key in ["z_ifc"]:
      # nc dimensions: lat, height, lon; target dimensions: lon, lat, level
      assert ncFile_fg.variables[key].dimensions == ('lat', 'height_3', 'lon')
      dataSingle[key] = np.transpose(ncFile_fg[key],(2, 0, 1))[...,::-1][...,:maxLevel+1]#reverse height order
      assert dataSingle[key].shape == shape3Dplus

    ncFile_fg.close()
    if verbosity > 1: print "closed fg nc"

    data = dataSingle

  #except IOError:
  except Exception as inst:
    print "ERROR:", fname_fg
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()
    raise


  #shapes may have changed!
  shape3Dplus = tuple(np.array(data["temp"].shape) + np.array([0,0,1]))
  shape3D = data["temp"].shape
  shape2D = shape3D[:2]

  #we can also fill the arrays partly!
  data["press_lev"] = np.zeros(shape3Dplus) + np.nan
  data["press_lev"][...,0] = data["pres_sfc"]

  data["relhum"] = q2rh(data["qv"],data["temp"],data["pres"]) * 100.

  data["hydro_q"] = np.zeros(data["qc"].shape + (nHydro,)) + np.nan
  data["hydro_q"][...,0] = data["qc"]
  data["hydro_q"][...,1] = data["qi"]
  data["hydro_q"][...,2] = data["qr"]
  data["hydro_q"][...,3] = data["qs"]
  data["hydro_q"][...,4] = data["qg"]
  data["hydro_q"][...,5] = data["qh"]

  pamData["hydro_n"] = np.zeros(data["qnc"].shape + (nHydro,)) + np.nan
  pamData["hydro_n"][...,0] = data["qnc"]
  pamData["hydro_n"][...,1] = data["qni"]
  pamData["hydro_n"][...,2] = data["qnr"]
  pamData["hydro_n"][...,3] = data["qns"]
  pamData["hydro_n"][...,4] = data["qng"]
  pamData["hydro_n"][...,5] = data["qnh"]

  varPairs = [["timestamp","timestamp"],["lat_2","lat"],["lon_2","lon"],
    ["u_10m","wind10u"],["v_10m","wind10v"],
    ["z_ifc","hgt_lev"],["pres","press"],["temp","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["hydro_n","hydro_n"],
    ["t_g","groundtemp"],["press_lev","press_lev"]]

  pamData = dict()
  for iconVar,pamVar in varPairs:
    pamData[pamVar] = data[iconVar]

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity
  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  # surface properties
  pamData['sfc_type'] = np.around(data['FR_LAND'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.full(shape2D, 'F', '|S1')

  pam.createProfile(**pamData)

  return pam

# Add regridding remarks
readIconLem2MomDataset.__doc__ += __ICDN_regridding_remarks


def readMesoNH(fnameBase,fnameExt,dataDir=".",debug=False,verbosity=0,dimX=160,dimY=160,dimZ=25,subGrid=None):

  variables = ['ALTITUDE','CLOUD','GEO','GRAUPEL','ICE','PRESSURE','RAIN','RHO','SEAPRES',
	      'SEATEMP','SNOW','SURFPRES','TEMP','TEMP2M','TEMPGRD','VAPOR','VAPOR2M','WINDLVLKB','ZSBIS']
  def vapor2rh(rv,t,p):
    rh = q2rh(rv/(1+rv),t,p)

    return rh

  nhydro = 5
  shape2D = (dimx,dimy,)
  shape3D = (dimx,dimy,dimz,)
  shape3Dplus = (dimx,dimy,dimz+1,)
  shape4D = (dimx,dimy,dimz,nhydro,)
  shape3DH = (dimx,dimy,nhydro+1,)
  shape3DM = (dimz,dimx,dimy,)

  try:
    data_mesonh = dict()
    for name in variables:
      fname = dataDir+'/'+fnameBase+name+fnameExt
      data_mesonh[name] = genfromtxt(fname)

    data = dict()
    data['hgt_lev'] = np.zeros(shape3Dplus)
    # GEO
    data['lon'] = data_mesonh['GEO'][:,0].reshape(shape2D)
    data['lat'] = data_mesonh['GEO'][:,1].reshape(shape2D)
    # time
    data['timestamp'] = np.zeros(shape2D)
    data['timestamp'][:,:] = 1000922400
    # 2D variables
    #data[''] = data_mesonh['SEAPRES'].reshape(shape2D)
    data['lfrac'] = np.zeros((shape2D))
    seatemp_tmp = data_mesonh['SEATEMP'].reshape(shape2D)
    data['lfrac'][np.where(seatemp_tmp > 998.)] = 1.
    data['groundtemp'] = data_mesonh['TEMPGRD'].reshape(shape2D)
    #data[''] = data_mesonh['TEMP2M'].reshape(shape2D)
    data['wind10u'] = (data_mesonh['WINDLVLKB']/np.sqrt(2.)).reshape(shape2D)
    data['wind10v'] = data['wind10u']
    data['hgt_lev'][...,0] = np.swapaxes(data_mesonh['ZSBIS'].reshape(shape2D),0,1)
    # 3D variables
    data['press'] = np.rollaxis(data_mesonh['PRESSURE'].reshape(shape3DM),0,3)
    data['temp'] = np.rollaxis(data_mesonh['TEMP'].reshape(shape3DM),0,3)
    data['hgt'] = np.rollaxis(data_mesonh['ALTITUDE'].reshape(shape3DM),0,3)
    data['hgt_lev'][:,:,1:] = (data['hgt'][:,:,1:]+data['hgt'][:,:,:-1])*0.5
    data['hgt_lev'][:,:,dimz] = data['hgt'][:,:,-1]+0.25*(data['hgt'][:,:,-1]-data['hgt'][:,:,-2])
    data['relhum'] = np.rollaxis(vapor2rh(data_mesonh['VAPOR'],data_mesonh['TEMP'],data_mesonh['PRESSURE']).reshape(shape3DM),0,3)*100.
    data['hydro_q'] = np.zeros(shape4D) + np.nan
    #data[''] = data_mesonh['RHO'].reshape(shape3DM)
    data['hydro_q'][...,0] = np.rollaxis(data_mesonh['CLOUD'].reshape(shape3DM),0,3)
    data['hydro_q'][...,1] = np.rollaxis(data_mesonh['ICE'].reshape(shape3DM),0,3)
    data['hydro_q'][...,2] = np.rollaxis(data_mesonh['RAIN'].reshape(shape3DM),0,3)
    data['hydro_q'][...,3] = np.rollaxis(data_mesonh['SNOW'].reshape(shape3DM),0,3)
    data['hydro_q'][...,4] = np.rollaxis(data_mesonh['GRAUPEL'].reshape(shape3DM),0,3)

    del data_mesonh

  except Exception as inst:
    if fname.split(".")[-1]!="nc":
      if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      os.system("rm -f "+tmpFile+"*")
    print "ERROR:", fname
    print type(inst)     # the exception instance
    print inst.args      # arguments stored in .args
    print inst
    if debug: import pdb;pdb.set_trace()

  #data['hydro_wp'] = np.zeros(shape3DH) + np.nan

  pam = pyPamtra()
  pam.df.readFile("../descriptorfiles/descriptor_file_meso-nh.txt")
  #pam.df.readFile("../descriptorfiles/descriptor_file_COSMO_1mom.txt")
  varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],["wind10u","wind10u"],["wind10v","wind10v"],
	      ["press","press"],["temp","temp"],["relhum","relhum"],["groundtemp","groundtemp"],['hgt_lev','hgt_lev'],['hydro_q','hydro_q']]

  pamData = dict()
  for dataVar,pamVar in varPairs:
    pamData[pamVar] = data[dataVar]
  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(data['lfrac'])
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  del data

  pam.createProfile(**pamData)

  return pam

def createUsStandardProfile(pam=pyPamtra(),**kwargs):
  '''
  Function to create clear sky US Standard Atmosphere.

  hgt_lev is teh only mandatory variables
  humidity will be set to zero if not provided, all other variables are guessed by "createProfile"

  values provided in kwargs will be passed to "createProfile", however, press_lev and temp_lev will overwritte us staandard if provided

  '''

  pamData = _createUsStandardProfile(**kwargs)

  pam.createProfile(**pamData)
  del kwargs

  return pam

def _createUsStandardProfile(**kwargs):
  '''
  HELPER
  Function to create clear sky US Standard Atmosphere.

  hgt_lev is the only mandatory variables
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

  noFiles = len(ncFiles)
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
      if noFiles > 1: assert joinDimension in ncData.dimensions.keys()
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

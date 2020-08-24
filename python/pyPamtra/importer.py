# -*- coding: utf-8 -*-

"""
Collection of importer modules for the different input types.

They all return a pyPamtra object that holds the profile information in the .p object.

In addition ncToDict is provide to easily read necdf files into a dictionary.
"""

from __future__ import division, print_function
import numpy as np
import datetime
import random
import os
import sys
import string
import glob
import warnings
import pdb

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
  """
  
  read WRF data set

  Parameters
  ----------
  fname : {str}
    path and name of file  atmospheric variables
  kind : {str}
    kind of data file (up to now only gce is implemented)
  
  Returns
  -------
  pam : pyPamtra Object

  Raises
  ------
  TypeError
  """

  import netCDF4

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
  """
  read COMSO one moment data set

  Parameters
  ----------
  fnames : {str}
    file name, wildCards allowed! can be either nc, nc and gz, or nc.gz files of name fnameInTar within tar file
  kind : {str}
    kind of COSMO file, right now only collum and synsat netcdf files are implemented
  descriptorFile : {str}
    path and name of descriptor file
  forecastIndex : {int}, optional
    index of forecast to be read (the default is 1)
  colIndex : {int or list of int}, optional
    which collum should be taken? list allowed! (the default is 0)
  tmpDir : {str}, optional
    [description] (the default is "/tmp/")
  fnameInTar : {str}, optional
    [description] (the default is "")
  concatenateAxis : {int}, optional
    concatenation along which axis (the default is 1)
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  df_kind : {str}, optional
    kind of data set. This defines the name of variables and dimensions. (the default is "default")
  constantFields : {str}, optional
    path and name of file with constant fields (the default is None)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)
  subGrid : {[int,int,int,int]}, optional
    array with indices [lon_start,lon_end,lat_start,lat_end] ((1,1) in model corresponds to (0,0) in python!) (the default is None)
 
  Returns
  -------
  pam : pyPamtra Object

  Raises
  ------
  IOError
  TypeError
  """

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
      if verbosity>0: print(fname)
      try:
        if fname.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fname.split(".")[-1]=="tar":
            if verbosity>3:print("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fname)
            if verbosity>2:print("created ", gzFile)
          else:
            gzFile = fname
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print("created ", tmpFile)
          ncFile = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print("opend ", tmpFile)
        else:
          ncFile = netCDF4.Dataset(fname,"r")
          if verbosity>1:print("opend ", fname)

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
        if verbosity>1:print("closed nc")
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
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
        for t in range(shape2D[1]):
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
          for key in list(data.keys()):
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

        ffOK += 1

      #except IOError:
      except Exception as inst:
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
          os.system("rm -f "+tmpFile+"*")
        print("ERROR:", fname)
        print(type(inst))     # the exception instance
        print(inst.args)      # arguments stored in .args
        print(inst)
        if debug: import pdb;pdb.set_trace()



    #shapes have changed!
    shape2D = np.shape(data["t_2m"])
    shape3Dplus = np.shape(data["hhl"])
    shape3D = np.shape(data["hfl"])

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
      if verbosity>0: print(fname)
      try:
        if fname.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fname.split(".")[-1]=="tar":
            if verbosity>3:print("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fname+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fname)
            if verbosity>2:print("created ", gzFile)
          else:
            gzFile = fname
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print("created ", tmpFile)
          ncFile = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print("opend ", tmpFile)
        else:
          ncFile = netCDF4.Dataset(fname,"r")
          if verbosity>1:print("opend ", fname)

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
        if verbosity>1:print("closed nc")
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
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
          for key in list(data.keys()):
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

        ffOK += 1


      #except IOError:
      except Exception as inst:
        if fname.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
          os.system("rm -f "+tmpFile+"*")
        print("ERROR:", fname)
        print(type(inst))     # the exception instance
        print(inst.args)      # arguments stored in .args
        print(inst)
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
  """
  import COSMO 2-moment dataset
  
  Parameters
  ----------
  fnamesA : {str}
    path and names of files with atmospheric variables. wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file
  descriptorFile : {str}
    path and name of descriptor file
  fnamesN : {str}, optional
    path and names of files with total numbers.  wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file (the default is None)
  kind : {str}, optional
    can be either new (default, new NARVAL calculations) or old (old NARVAL calculations) (the default is 'new')
  forecastIndex : {int}, optional
    index of forecast to be read (the default is 1)
  tmpDir : {str}, optional
    [description] (the default is "/tmp/")
  fnameInTar : {str}, optional
    [description] (the default is "")
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  df_kind : {str}, optional
    kind of data set. This defines the name of variables and dimensions. (the default is "default")
  constantFields : {str}, optional
    path and name of file with constant fields (the default is None)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)
  subGrid : {[int,int,int,int]}, optional
    array with indices [lon_start,lon_end,lat_start,lat_end] ((1,1) in model corresponds to (0,0) in python!) (the default is None)
  
  Returns
  -------
  pam : pyPamtra Object

  Raises
  ------
  IOError
    on file error

  """

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
      if verbosity>0: print(fnameA)
      try:
        if fnameA.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fnameA.split(".")[-1]=="tar":
            if verbosity>3:print("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fnameA)
            if verbosity>2:print("created ", gzFile)
          else:
            gzFile = fnameA
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print("created ", tmpFile)
          ncFileA = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print("opend ", tmpFile)
        else:
          ncFileA = netCDF4.Dataset(fnameA,"r")
          if verbosity>1:print("opend ", fnameA)

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
        if verbosity>1:print("closed nc")
        if fnameA.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
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
          for key in list(data.keys()):
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

        ffOK += 1


      #except IOError:
      except Exception as inst:
        if fnameA.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
          os.system("rm -f "+tmpFile+"*")
        print("ERROR:", fnameA)
        print(type(inst))     # the exception instance
        print(inst.args)      # arguments stored in .args
        print(inst)
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
      if verbosity>0: print(fnameA, fnameN)
      try:
        if fnameA.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fnameA.split(".")[-1]=="tar":
            if verbosity>3:print("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameA+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fnameA)
            if verbosity>2:print("created ", gzFile)
          else:
            gzFile = fnameA
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print("created ", tmpFile)
          ncFileA = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print("opend ", tmpFile)
        else:
          ncFileA = netCDF4.Dataset(fnameA,"r")
          if verbosity>1:print("opend ", fnameA)

        # number density files
        if fnameN.split(".")[-1]!="nc":
          tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
          if fnameN.split(".")[-1]=="tar":
            if verbosity>3:print("tar -O --wildcards -x "+fnameInTar+" -f "+fnameN+">"+tmpFile+".gz")
            os.system("tar -O --wildcards -x "+fnameInTar+" -f "+fnameN+">"+tmpFile+".gz")
            gzFile = tmpFile+".gz"
            if os.stat(gzFile).st_size == 0:
              os.system("rm -f "+tmpFile+"*")
              raise IOError("fnameInTar not found in "+fnameN)
            if verbosity>2:print("created ", gzFile)
          else:
            gzFile = fnameN
          os.system("zcat "+gzFile+">"+tmpFile)
          if verbosity>1:print("created ", tmpFile)
          ncFileN = netCDF4.Dataset(tmpFile,"r")
          if verbosity>1:print("opend ", tmpFile)
        else:
          ncFileN = netCDF4.Dataset(fnameN,"r")
          if verbosity>1:print("opend ", fnameN)

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
        if verbosity>1:print("closed nc")
        if fnameA.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
          os.system("rm -f "+tmpFile+"*")
        for var in variables4DN:
          if subGrid == None:
            dataSingle[var] = np.swapaxes(ncFileN.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel]#reverse height order
          else:
            dataSingle[var] = np.swapaxes(ncFileN.variables[var][forecastIndex],0,2)[...,::-1][...,:maxLevel][subGrid[0]:subGrid[1],subGrid[2]:subGrid[3],:]
        ncFileN.close()
        if verbosity>1:print("closed nc")
        if fnameN.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
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
          for key in list(data.keys()):
            data[key] = np.ma.concatenate((data[key],dataSingle[key],),axis=concatenateAxis)

        ffOK += 1


      #except IOError:
      except Exception as inst:
        if fnameA.split(".")[-1]!="nc":
          if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
          os.system("rm -f "+tmpFile+"*")
        print("ERROR:", fnameA)
        print(type(inst))     # the exception instance
        print(inst.args)      # arguments stored in .args
        print(inst)
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

def readCosmoDe2MomDatasetOnFlightTrack(fnameA,descriptorFile,tmpDir="/tmp/",debug=False,verbosity=0,maxLevel=0):
  """  
  import COSMO 2-moment dataset extracted on a HALO flight track
  
  Parameters
  ----------
  fnameA : {str}
    path and names of file with atmospheric variables, wildCards allowed! can be either nc file, nc,gz file or nc.gz file of name fnameInTar within tar file
  descriptorFile : {str}
    path and name of descriptor file
  tmpDir : {str}, optional
    directory to unpack temporary files (the default is "/tmp/")
  debug : {bool}, optional
    stop and load debugger on exception (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)

  Returns
  --------
  pam : pyPamtra Object

  """
  import netCDF4

  nHydro = 6

  variables1D = ["T_G","PS","U_10M","V_10M","FR_LAND"]
  variables2D = ["T","P","QV","QC","QI","QR","QS","QG","QH","QNCLOUD","QNICE","QNRAIN","QNSNOW","QNGRAUPEL","QNHAIL"]

  if verbosity>0: print(fnameA)
  try:
    if fnameA.split(".")[-1]!="nc":
      tmpFile = tmpDir+"/pyPamtraImport_netcdf_"+''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(5))+".nc"
      gzFile = fnameA
      os.system("zcat "+gzFile+">"+tmpFile)
      if verbosity>1:print("created ", tmpFile)
      data = ncToDict(tmpFile)
      if verbosity>1:print("read ", tmpFile)
    else:
      data = ncToDict(fnameA)
      if verbosity>1:print("read ", fnameA)

    if debug: import pdb;pdb.set_trace()
    if verbosity>1:print("closed nc")
    if fnameA.split(".")[-1]!="nc":
      if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
      os.system("rm -f "+tmpFile+"*")

  except Exception as inst:
    if fnameA.split(".")[-1]!="nc":
      if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
      os.system("rm -f "+tmpFile+"*")
    print("ERROR:", fnameA)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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

def readCosmoReAn2km(constantFields,fname,descriptorFile,forecastIndex = 1,tmpDir="/tmp/",fnameInTar="",debug=False,verbosity=0,maxLevel=51,subGrid=None):
  """
  reads COSMO 2 km Re-Analyses
   
  Parameters
  ----------
  constantFields : {str}
    path and file name of constant file
  fname : {str}
    path and file name of data file
  descriptorFile : {str}
    path and name of descriptor file
  forecastIndex : {int}, optional
    index of forecast to be read (the default is 1)
  tmpDir : {str}, optional
    temporary directory to unpack data (the default is "/tmp/", which [default_description])
  fnameInTar : {str}, optional
    file name within tar file (the default is "")
  debug : {bool}, optional
    switch on debugging mode(the default is False)
  verbosity : {int}, optional
    [description] (the default is 0)
  maxLevel : {int}, optional
    index of maximum level (the default is 51)
  subGrid : {[int,int,int,int]}, optional
    read subgrid of whole field (the default is None)
  
  Returns
  -------
  pyPamtra object
  """

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
      if verbosity>1: print(var)
      selected_grbs = grbsC(shortName=var)
      for i in range(nLev-maxLevel,nLev):
        data[var][...,nLev-1-i] = selected_grbs[i].values

    selected_grbs = grbsC(shortName=variables2DC)
    for var in selected_grbs:
      if verbosity>1: print(var.shortName)
      data[var.shortName] = var.values

    grbsC.close()

    grbs = pygrib.open(fname)

    selected_grbs = grbs(shortName=variables2D)
    for var in selected_grbs:
      if verbosity>1: print(var.shortName)
      data[var.shortName] = var.values

    for var in variables3D:
      if verbosity>1: print(var)
      data[var] = np.zeros(shape3D) + np.nan
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
        data[selected_grbs[i].shortName][...,nLev-2-i] = selected_grbs[i].values

    data['hydro_q'] = np.zeros(shape4D) + np.nan

    for j,var in enumerate(variables4D):
      if verbosity>1: print(var)
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
        data['hydro_q'][...,nLev-2-i,j] = selected_grbs[i].values

    grbs.close()
  except Exception as inst:
    if fname.split(".")[-1]!="nc":
      if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
      os.system("rm -f "+tmpFile+"*")
    print("ERROR:", fname)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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
      if verbosity>1: print(var)
      data[var] = np.zeros(shape3D) + np.nan
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
        data[selected_grbs[i].shortName][...,nLev-2-i] = selected_grbs[i].values

    data['hydro_q'] = np.zeros(shape4D) + np.nan

    for j,var in enumerate(variables4D):
      if verbosity>1: print(var)
      selected_grbs = grbs(shortName=var)
      for i in range(nLev-maxLevel,nLev-1):
        data['hydro_q'][...,nLev-2-i,j] = selected_grbs[i].values


    grbs.close()

    grbs = pygrib.open(fname+'00')

    selected_grbs = grbs(shortName=variables2D)
    for var in selected_grbs:
      if verbosity>1: print(var.shortName)
      data[var.shortName] = var.values


    grbs.close()

  except Exception as inst:
    #if fname.split(".")[-1]!="nc":
      #if verbosity>1:print "removing ", glob.glob(tmpFile+"*")
      #os.system("rm -f "+tmpFile+"*")
    print("ERROR:", fname)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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
  """
  import ICON LEM 1-moment dataset

  It is Icon with cosmo physics
  
  Parameters
  ----------
  fname_fg : {str}
    path and name of file with atmospheric variables.
  descriptorFile : {str}
    path and name of descriptor file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  constantFields : {str}, optional
    path and name of file with constant fields (the default is None)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)
  
  Returns
  -------
  pam : pyPamtra Object
  
  Raises
  ------
  IOError
  """
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

  if verbosity>0: print(fname_fg)

  if not fname_fg.endswith('.nc'):
    raise IOError("fname_fg has to be .nc.", fname_fg)
  if not '_fg_' in os.path.basename(fname_fg):
    raise IOError("fname_fg has to contain '_fg_'.", fname_fg)

  dataSingle = dict()
  try:
    ncFile_const = netCDF4.Dataset(constantFields, "r")
    if verbosity > 1: print("opened ", constantFields)

    for var in variables2D_const:
      # nc dimensions: lat, lon; target dimensions: lon, lat
      assert ncFile_const.variables[var].dimensions == ('lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_const.variables[var],0,1)

    ncFile_const.close()
    if verbosity > 1: print("closed const nc")

    ncFile_fg = netCDF4.Dataset(fname_fg, "r")
    if verbosity > 1: print("opened ", fname_fg)

    fname_cloud = '_cloud_'.join(fname_fg.rsplit('_fg_',1)) # right-replace _fg_ with _cloud
    ncFile_cloud = netCDF4.Dataset(fname_cloud, "r")
    if verbosity > 1: print("opened ", fname_cloud)

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
    if verbosity > 1: print("closed fg nc")
    ncFile_cloud.close()
    if verbosity > 1: print("closed cloud nc")

    data = dataSingle

  #except IOError:
  except Exception as inst:
    print("ERROR:", fname_fg)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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
  assert np.all(np.logical_or(0<=pamData['sfc_type'], pamData['sfc_type'] <=1))
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][pamData['sfc_type'] == 0] = 'F' # ocean
  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'S' # land

  pam.createProfile(**pamData)

  return pam

# Add regridding remarks
readIconLem1MomDataset.__doc__ += __ICDN_regridding_remarks


def readIconLem2MomDataset(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0):
  """
  import ICON LEM 2-moment dataset

  It is ICON with cosmo physics
  
  Parameters
  ----------
  fname_fg : {str}
    path and name of file with atmospheric variables.
  descriptorFile : {str}
    path and name of descriptor file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  constantFields : {str}, optional
    path and name of file with constant fields (the default is None)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)
  
  Returns
  -------
  pam : pyPamtra Object
  
  Raises
  ------
  IOError
  """
  import netCDF4

  assert constantFields
  forecastIndex = 0 # time step in forecast
  nHydro = 6

  data = dict()

  variables2D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
  variables3D = ["t_g","pres_sfc"]
  variables4D_10m = ["u_10m","v_10m"]
  variables4D = ["temp","pres","qv","qc","qi","qr","qs","qg","qh","qnc","qni","qnr","qns","qng","qnh"]

  if verbosity>0: print(fname_fg)

  if not fname_fg.endswith('.nc'):
    raise IOError("fname_fg has to be .nc.", fname_fg)
  if not '_fg_' in os.path.basename(fname_fg):
    raise IOError("fname_fg has to contain '_fg_'.", fname_fg)

  dataSingle = dict()
  try:
    ncFile_const = netCDF4.Dataset(constantFields, "r")
    if verbosity > 1: print("opened ", constantFields)

    for var in variables2D_const:
      # nc dimensions: lat, lon; target dimensions: lon, lat
      assert ncFile_const.variables[var].dimensions == ('lat', 'lon')
      dataSingle[var] = np.swapaxes(ncFile_const.variables[var],0,1)

    ncFile_const.close()
    if verbosity > 1: print("closed const nc")

    ncFile_fg = netCDF4.Dataset(fname_fg, "r")
    if verbosity > 1: print("opened ", fname_fg)


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
    if verbosity > 1: print("closed fg nc")

    data = dataSingle

  #except IOError:
  except Exception as inst:
    print("ERROR:", fname_fg)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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
  assert np.all(np.logical_or(0<=pamData['sfc_type'], pamData['sfc_type'] <=1))
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][pamData['sfc_type'] == 0] = 'F' # ocean
  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'S' # land

  pam.createProfile(**pamData)

  return pam

# Add regridding remarks
readIconLem2MomDataset.__doc__ += __ICDN_regridding_remarks

def readHIRHAM(dataFile,additionalFile,topoFile,descriptorFile,grid=[0,200,0,218],timestep=0,debug=False,verbosity=0):
  """  
  read  HIRHAM output into PAMTRA

  it reads the output produce by Ines Haberstedt at the AWI in Potsdam
  
  Parameters
  ----------
  dataFile : {str}
    path and name of data file
  additionalFile : {str}
    path and name of additional file
  topoFile : {str}
    path and name of file with topography information
  descriptorFile : {str}
    path and name of descriptor file
  grid : {[int,int,int,int]}, optional
    subgrid to read (the default is [0,200,0,218] the full grid)
  timestep : {int}, optional
    whcih if the 8 available time steps should be read (the default is 0, max is 7)
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  
  Returns
  -------
  pam : pyPamtra Object
  """
  import netCDF4

  data = dict() #  this will hold the data before we create the pamtra object data

  a = grid[0]     # 0 # 70
  b = grid[1]     #150 # 82
  # y-axis 0:218
  c = grid[2]     #80 # 50
  d = grid[3]    #218 # 70
  # time
  t = timestep

  fo = netCDF4.Dataset(dataFile, 'r')

  data['lon'] = fo['XLONG'][a:b,c:d]
  data['lat'] = fo['XLAT'][a:b,c:d]
  ts = (datetime.datetime.strptime(str(fo['DateTime'][t]),'%Y%m%d%H')- datetime.datetime(1970,1,1)).total_seconds()
  data['hgt'] = np.moveaxis(fo['GHT'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1] # geopotnetial height [m] used as hgt
  data['temp'] = np.moveaxis(fo['TT'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1] # temperature [K]
  data['relhum'] = np.moveaxis(fo['RH'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1] # relative humidity [%]
  Frain_ls = np.moveaxis(fo['LSRAIN'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]     # Flux large scale cloud rain in kg m^-2 s^ 
  Frain_c = np.moveaxis(fo['CCRAIN'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]      # Flux convective cloud rain in kg m^-2 s^
  Fsnow_ls = np.moveaxis(fo['LSSNOW'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]     # Flux large scale cloud snow in kg m^-2 s^
  Fsnow_c = np.moveaxis(fo['CCSNOW'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]      # Flux convective cloud snow in kg m^-2 s^
  qvapor = np.moveaxis(fo['QVAPOR'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1] #[:,:,c:d,g:h]         # Water vapor mixing ratio in kg kg-1
  qcloud = np.moveaxis(fo['QCLOUD'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]              # Cloud water mixing ratio in kg kg-1
  qice = np.moveaxis(fo['QICE'][t,...],[0,1,2],[2,0,1])[a:b,c:d,::-1]                  # Ice mixing ratio in kg kg-1
  # LANDMASK = (fo['LANDKASK'][t,...])[a:b,c:d]            # 0 for water 1 for land (more detailed variable in /data/mod/hirham/5/3hr/topo_lsm.nc)
  data['groundtemp']= (fo['SKINTEMP'][t,...])[a:b,c:d]             # Surface skin  temperature (time, lat,lon)
  data['wind10u'] = (fo['U10M'][t,...])[a:b,c:d] # 10 m windspeed zonal component [m/s]
  data['wind10v'] = (fo['V10M'][t,...])[a:b,c:d] # 10 m windspeed meridional component [m/s]
  a_mid = (fo['hyam'][:])[::-1]  # A coefficients at layer midpoint [Pa]
  b_mid = (fo['hybm'][:])[::-1]  # B coefficients at layer midpoints [1]

  fo.close()

  # loading additional files needed for PAMTRA simulations
  # surface pressure
  fo2 = netCDF4.Dataset(additionalFile,'r')
  ps = (fo2['ps'][0,...])[a:b,c:d]
  seaicefraction = (fo2['seaice'][0,...])[a:b,c:d] # [0,1] 0 = ocean, 1 = seaice
  sealandfraction = (fo2['slf'][0,...])[a:b,c:d] # [0,1] 0 = sea/water, 1 = land; compared to LANDMASk it is a continuous value
  fo2.close()

  # surf_geopotential
  fo_inv = netCDF4.Dataset(topoFile,'r')
  surf_geopotential = fo_inv['fis'][a:b,c:d]#[c:d,g:h] # surface geopotential (orography) in m^2/s^2
  fo_inv.close()

  # construct all shapes needed for PAMTRA
  shape2D = ps.shape
  shape3D = data['hgt'].shape
  shape3Dplus = (shape3D[0],shape3D[1],shape3D[2]+1)
  shape4D = (shape3D[0],shape3D[1],shape3D[2],4)

  # let's set the time for each profile
  data['timestamp'] = np.zeros(shape2D)
  data['timestamp'][:] = ts

  # Calculating pressure variable from a and b hybrid coefficients:
  data['press'] = np.zeros(shape3D)

  for i in range(shape3D[2]): 
      if i == 0:
          data['press'][:,:,i] = a_mid[i] + b_mid[i]*ps
      else:
          # data['press'][:,:,i] = a_mid[i-1]*0.5 + (b_mid[i-1]*0.5)*ps
          data['press'][:,:,i] = a_mid[i] + (b_mid[i])*ps

  Rspec = 8.31432e3 # [Nm/kmol K]
  rho_dry =  data['press']/(Rspec * data['temp']) # IGL: P = rho*Rspec*T

  # calculate height of levels
  data['hgt_lev'] = np.zeros(shape3Dplus)
  data['hgt_lev'][...,0] = surf_geopotential/9.81

  for i in range(shape3D[2]-1):
    data['hgt_lev'][...,i+1] = (data['hgt'][...,i+1] + data['hgt'][...,i])*0.5

  data['hgt_lev'][...,-1] = data['hgt'][...,-1] + (data['hgt'][...,-1] - data['hgt'][...,-2])*0.5

  ################################################################################################################
  #               Calculate mixing ratios from fluxes from HIRHAM5. Scheme is like in HIRLAM
  ################################################################################################################

  # Constants:
  R = 287.05 # [J/kg K] specific gas constant of dry air
  # snow:
  Cpr = 1. # assumtion that we have snow/ice/rain in the whole grid cell
  rho = 1000. # density of water [kg/m^3] 
  a11 = 3.29
  b10 = 0.16

  #rain:
  a10 = 90.8
  n0r = 8.e6 # [m^-4] intercept parameter
  n0s = 3.e6 # [m^-4] intercept parameter

  # ice - same constants as for snow

  # Mass mixing ratio of snow within the fraction Cpr of the grid cell covered with snow is obtained from the snow fall rate:

  # rho - air density --> calculate from IGL: rho_dry = P/R*T = 1.225 kg/m^3 (P -sea level and T at standard atmosp pressure)
  rho0 = 1.3           # rho0 - reference density of air: 1.3 kg/m^3
  rho_i = 500.         # rhoi - density of cloud ice: 500 kg/m^3
  rho_s = 100.         # rhos - bulk density of snow: 100 kg/m^3
  rho_r = 1000.        # rhow - density of water: 1000 kg/m^3

  # rs - mass mixing ratio of snow
  # rr - mass mixing ratio of rain
  # ri - mas mixing ratio of falling ice
  # Fsnow - snow flux 
  # Frain - rain flux
  # Fitop - grid cell mean sedimantation flux
  # Cpr - fraction of the grid cell covered with snow

  ## 1.) snow:
  qsnow = ((Fsnow_ls + np.abs(Fsnow_c))/(Cpr * a11))**(1/(1 + b10))

  # rs = ((Fsnow_ls)/Cpr * a11)**(1./(1. + b10))

  ##2.) rain;
  ## vr - mass-weighted fall velocity of rain drops parameterized according to Kessler (1969)

  qrain = ( (Frain_ls + np.abs(Frain_c)) / ( Cpr*a10 *  n0r**(-1./8.)  * np.sqrt(rho0/rho_dry) ) )**(8./9.)

  #-------------------------------------------------------------------------------------------------

  data['hydro_q'] = np.zeros(shape4D)

  data['hydro_q'][...,0] = qcloud
  data['hydro_q'][...,1] = qice
  data['hydro_q'][...,2] = qrain
  data['hydro_q'][...,3] = qsnow

  pamData = dict()  # Creating pamtra dictionary where all variables will be added

  varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],["wind10u","wind10u"],["wind10v","wind10v"],["hgt","hgt"],
    ["hgt_lev","hgt_lev"],["press","press"],["temp","temp"],["relhum","relhum"],["hydro_q","hydro_q"],["groundtemp","groundtemp"]]

  for hirhamVar,pamVar in varPairs:
      pamData[pamVar] = data[hirhamVar]

  # surface properties
  pamData['sfc_type'] = np.zeros(shape2D)
  pamData['sfc_type'] = np.around(sealandfraction)
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'

  pamData['sfc_refl'][(pamData['sfc_type'] == 1)] = 'S'
  pamData['sfc_model'][(pamData['sfc_type'] == 1)] = 0.
  # pamData['sfc_type'][(pamData['sfc_type'] == 0.) & (pamData['groundtemp'] < 270.)] = 2.
  pamData['sfc_slf'] = sealandfraction
  pamData['sfc_sif'] = seaicefraction

  pam = pyPamtra()
  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)
  pam.createProfile(**pamData)

  return pam

def readECMWF(fname,constantFile,descriptorFile,landseamask,debug=False,verbosity=0,grid=[0,-1,0,-1]):
  """
  read ECMWF IFS data

  reads the data produced for the NAWDEX campaign by Heini Wernli at ETH Zuerich
  
  Q          specific humidity                   kg/kg
  LWC    specific liquid water content    kg/kg
  IWC    specific ice water content        kg/kg
  RWC    specific rain water content      kg/kg
  SWC     specific snow water content   kg/kg
  T            temperature                           K    -
  U            u wind component                m/s
  V            v wind component                m/s
  PS          surface pressure                 Pa
  SLP        mean sea level pressure      Pa
  P            pressure                              Pa
  SKT        skin temperature                  K
  u_wind_10m  10 meter wind u component  m/s
  v_wind_10m  10 meter wind v component  m/s
  lsm   land-sea-mask            0-1   (0=sea, 1=land)
  
  Parameters
  ----------
  fname : {str}
    path and name of data file
  constantFile : {str}
    path and name of constant fields
  descriptorFile : {str}
    path and name of descriptor file
  landseamask : {str}
    path and name of landseamask file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  grid : {[int,int,int,int]}, optional
    read only subgrid of whole model field
  
  Returns
  -------
  pam : pyPamtra Object
  """
  import netCDF4

  nHydro = 4

  variables1D = ["SKT","SLP","U_10M","V_10M","lsm"]
  variables2D = ["T","P","Q","LWC","IWC","RWC","SWC"]

  if verbosity>0: print(fname)
  if debug: import pdb;pdb.set_trace()

  dataTmp = ncToDict(fname)
  dataC = ncToDict(constantFile,keys=['aklev','bklev','aklay','bklay'])
  pamData = dict()

  a = grid[0]  
  b = grid[1]  
  c = grid[2]  
  d = grid[3]  


  data = dict()
  data['T'] = dataTmp['T'][0,:,c:d,a:b] + 273.15 # -> K

  shape3D = (data['T'].shape[2],data['T'].shape[1],data['T'].shape[0])
  shape3Dplus = (shape3D[0],shape3D[1],shape3D[2]+1)
  shape2D = shape3D[:2]

  dataLSM = ncToDict(landseamask)
  data['LSM'] = dataLSM['LSM'][0,0,c:d,a:b]

  # bring everything to desired units
  data['PS'] = dataTmp['PS'][0,0,c:d,a:b] * 100. # -> Pa

  # calculate pressure on half-levels with SLP as level 0
  pamData['press_lev'] = np.empty(shape3Dplus)
  pamData['press'] = np.empty(shape3D)
  pamData['hgt_lev'] = np.empty(shape3Dplus)

  pamData['press_lev'][:,:,0] = np.swapaxes(data['PS'][...],0,1)

  pamData['temp'] = np.empty(shape3D)
  pamData['temp'] = np.swapaxes(data['T'][:,:,:],0,2)
  pamData['temp_lev'] = np.empty(shape3Dplus)
  pamData['temp_lev'][...,1:-1] = (pamData['temp'][...,1:] + pamData['temp'][...,0:-1])*0.5
  pamData['temp_lev'][...,-1] = pamData['temp_lev'][...,-2]
  pamData['temp_lev'][...,0] = np.swapaxes(dataTmp['SKT'][0,0,c:d,a:b],0,1)

  for i in range(shape2D[0]):
    for j in range(shape2D[1]):
      pamData['press_lev'][i,j,1::] = data['PS'][j,i]*dataC['bklev']+dataC['aklev']*100.
      pamData['press'][i,j,:] = data['PS'][j,i]*dataC['bklay']+dataC['aklay']*100.
      # pamData['hgt_lev'][i,j,:] = (10**(np.log10(pamData['press_lev'][i,j,:]/data['PS'][0,0,j,i])/5.2558797)-1.)/-6.8755856e-6
      pamData['hgt_lev'][i,j,:] = (((data['PS'][j,i]/pamData['press_lev'][i,j,:])**(1./5.257)-1.)* pamData['temp_lev'][i,j,:])/0.0065# (10**(np.log10(pamData['press_lev'][i,j,:]/data['PS'][0,0,j,i])/5.2558797)-1.)/-6.8755856e-6

  pamData['relhum'] = np.empty(shape3D)
  pamData["relhum"][:,:,:] = (q2rh(np.swapaxes(dataTmp["Q"][0,:,c:d,a:b],0,2),pamData["temp"][:,:,:],pamData["press"][:,:,:]) * 100.)#[:,-2::-1]

  pamData["hydro_q"] = np.zeros(shape3D + (nHydro,)) + np.nan
  pamData["hydro_q"][:,:,:,0] = np.swapaxes(dataTmp["LWC"][0,:,c:d,a:b],0,2)
  pamData["hydro_q"][:,:,:,1] = np.swapaxes(dataTmp["IWC"][0,:,c:d,a:b],0,2)
  pamData["hydro_q"][:,:,:,2] = np.swapaxes(dataTmp["RWC"][0,:,c:d,a:b],0,2)
  pamData["hydro_q"][:,:,:,3] = np.swapaxes(dataTmp["SWC"][0,:,c:d,a:b],0,2)

  varPairs = [["U10","wind10u"],["V10","wind10v"],["SKT","groundtemp"]]

  for ecmwfVar,pamVar in varPairs:
    pamData[pamVar] = np.swapaxes(dataTmp[ecmwfVar][0,0,c:d,a:b],0,1)


  # surface properties
  pamData['sfc_type'] = np.around(np.swapaxes(data['LSM'],0,1))
  pamData['sfc_model'] = np.zeros(shape2D)
  pamData['sfc_refl'] = np.chararray(shape2D)
  pamData['sfc_refl'][:] = 'F'
  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'S'

  pam = pyPamtra()
  pam.set["pyVerbose"]= verbosity

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  pam.createProfile(**pamData)

  return pam

def readMesoNH(fnameBase,fnameExt,dataDir=".",debug=False,verbosity=0,dimX=160,dimY=160,dimZ=25,subGrid=None):
  """
  read MesoNH output as produced by Jean-Pierre Chaboureau for the SIMGEO study

  Parameters
  ----------
  fnameBase : {str}
    base name of file name
  fnameExt : {str}
    extension of file name
  dataDir : {str}, optional
    path to data directory (the default is ".")
  debug : {bool}, optional
    stop and load debugger on exception (the default is False) (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  dimX : {int}, optional
    dimension in x (the default is 160)
  dimY : {int}, optional
    dimension in y (the default is 160)
  dimZ : {int}, optional
    dimension in z (the default is 25)
  subGrid : {[int,int,int,int]}, optional
    only read subgrid (the default is None)
  
  Returns
  -------
  pam : pyPamtra Object
  """
  variables = ['ALTITUDE','CLOUD','GEO','GRAUPEL','ICE','PRESSURE','RAIN','RHO','SEAPRES',
              'SEATEMP','SNOW','SURFPRES','TEMP','TEMP2M','TEMPGRD','VAPOR','VAPOR2M','WINDLVLKB','ZSBIS']
  def vapor2rh(rv,t,p):
    rh = q2rh(rv/(1+rv),t,p)

    return rh

  if debug: import pdb;pdb.set_trace()
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
      if verbosity>1:print("removing ", glob.glob(tmpFile+"*"))
      os.system("rm -f "+tmpFile+"*")
    print("ERROR:", fname)
    print(type(inst))     # the exception instance
    print(inst.args)      # arguments stored in .args
    print(inst)
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

def readIcon1momMeteogram(fname, descriptorFile, debug=False, verbosity=0, timeidx=None, hydro_content=[1.,1.,1.,1.,1.]):
  '''
  import ICON LEM 2-moment dataset output cross section over a site (Meteogram)
  
  WARNING: This importer has been designed upon Icon meteogram output generated by Dr. Vera Schemann
  Might not work on other files.
  
  fname = filename of the Icon output
  descriptorFile = pyPamtra descriptorFile object or filename of a proper formatted descriptorFile
  debug = flag causing stop and load of debugger upon raised exception
  verbosity = pyPamtra.pyVerbose verbosity level
  timeidx = indexes of timestep to be computed. Used to slice just a few profiles out of the .nc file for rapid testing
  
  Pamtra assumes that the apart from height, the first given dimension is latitude, here coordinates do not change across the dimension, but time does, thus generating a meteogram
  
  03/03/2018 - Davide Ori - dori@uni-koeln.de
  '''
  
  import netCDF4
  
  ICON_file = netCDF4.Dataset(fname, mode='r')
  vals = ICON_file.variables

  ## THIS PART IS TROUBLESOME, BUT PEOPLE MIGHT NEED STATION DETAILS ##
#  station = ICON_file.station.split()
#  lon = float([s for s in station if '_lon' in s][0].split('=')[-1])
#  lat = float([s for s in station if '_lat' in s][0].split('=')[-1])
#  height = float([s for s in station if '_hsurf' in s][0].split('=')[-1])
#  name = [s for s in station if '_name' in s][0].split('=')[-1]
#  frland = float([s for s in station if '_frland' in s][0].split('=')[-1])
#  fc = float([s for s in station if '_fc' in s][0].split('=')[-1])
#  soiltype = int([s for s in station if '_soiltype' in s][0].split('=')[-1])
#  tile_frac=[',  u'1.]',
#  tile_luclass=[4]
  
  pamData = dict() # empty dictionary to store pamtra Data

  hgt_key = 'height'
  if hgt_key in vals:
    hgt_key = 'height'
  elif 'heights' in vals:
    hgt_key = 'heights'
  else:
    raise AttributeError('ICON file does not have a valid height label (I know only height, heights, height_2 and heights_2)')
  
  Nh = len(vals[hgt_key])-1
  Nt = len(vals['time'])
  
  #if vals.has_key('QG'): # It might be either 3 or 2 ice classes
  nhydros = 5
  #else:
  #  nhydros = 4
  if timeidx is None:
    timeidx = timeidx = np.arange(0,Nt)
  shapeSFC = (len(timeidx),)
    
  date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
  timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
  pamData['timestamp'] = timestamp  
  pamData['hgt_lev'] = np.tile(np.flip(vals[hgt_key],0),(len(timeidx),1)) # heights at which fields are defined

  pamData['press']  = np.flip(vals['P'][timeidx],1)    # pressure 
  pamData['temp']   = np.flip(vals['T'][timeidx],1)    # temperature
  wind_u = np.flip(vals['U'][timeidx],1)               # zonal wind speed
  wind_v = np.flip(vals['V'][timeidx],1)               # meridional wind speed
  pamData['wind_uv'] = np.hypot(wind_u, wind_v)
  wind_w = np.flip(vals['W'][timeidx],1)               # vertical wind speed
  pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
  pamData['relhum'] = np.flip(vals['REL_HUM'][timeidx],1)

  # Read hydrometeors content
  hydro_cmpl = np.zeros((len(timeidx),Nh,nhydros))
  hydro_cmpl[...,0] = hydro_content[0]*np.flip(vals['QC'][timeidx],1)   # specific cloud water content
  hydro_cmpl[...,1] = hydro_content[1]*np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
  hydro_cmpl[...,2] = hydro_content[2]*np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
  hydro_cmpl[...,3] = hydro_content[3]*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
  if 'QG' in vals:
    hydro_cmpl[...,4] = hydro_content[4]*np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
  else:
    hydro_cmpl[...,4] = 0.0*np.flip(vals['QC'][timeidx],1)   # graupel set to 0 (avoid double descriptorFile)

  pamData["hydro_q"] = hydro_cmpl
  
  # surface properties
  pamData['wind10u'] = vals['U10M'][timeidx]
  pamData['wind10v'] = vals['V10M'][timeidx]
  pamData['groundtemp'] = vals['T_S'][timeidx]

# surface properties
  pamData['sfc_type']  = np.ones(pamData['groundtemp'].shape)
  pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
  pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
  pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land

  pam = pyPamtra()
  pam.set['pyVerbose'] = verbosity

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else: # is an iterable containing arrays
    for df in descriptorFile:
      pam.df.addHydrometeor(df)
  
  pam.createProfile(**pamData)
  
  return pam


def readIcon2momMeteogram(fname, descriptorFile, debug=False, verbosity=0, timeidx=None, hydro_content=[1.,1.,1.,1.,1.,1.]):
  '''
  import ICON LEM 2-moment dataset output cross section over a site (Meteogram)
  
  WARNING: This importer has been designed upon Icon meteogram output generated by Dr. Vera Schemann
  Might not work on other files.
  
  fname = filename of the Icon output
  descriptorFile = pyPamtra descriptorFile object or filename of a proper formatted descriptorFile
  debug = flag causing stop and load of debugger upon raised exception
  verbosity = pyPamtra.pyVerbose verbosity level
  timeidx = indexes of timestep to be computed. Used to slice just a few profiles out of the .nc file for rapid testing
  
  Pamtra assumes that the apart from height, the first given dimension is latitude, here coordinates do not change across the dimension, but time does, thus generating a meteogram
  
  03/03/2018 - Davide Ori - dori@uni-koeln.de
  '''
  
  import netCDF4
  
  ICON_file = netCDF4.Dataset(fname, mode='r')
  vals = ICON_file.variables

  ## THIS PART IS TROUBLESOME, BUT PEOPLE MIGHT NEED STATION DETAILS ##
#  station = ICON_file.station.split()
#  lon = float([s for s in station if '_lon' in s][0].split('=')[-1])
#  lat = float([s for s in station if '_lat' in s][0].split('=')[-1])
#  height = float([s for s in station if '_hsurf' in s][0].split('=')[-1])
#  name = [s for s in station if '_name' in s][0].split('=')[-1]
#  frland = float([s for s in station if '_frland' in s][0].split('=')[-1])
#  fc = float([s for s in station if '_fc' in s][0].split('=')[-1])
#  soiltype = int([s for s in station if '_soiltype' in s][0].split('=')[-1])
#  tile_frac=[',  u'1.]',
#  tile_luclass=[4]
  
  pamData = dict() # empty dictionary to store pamtra Data

  hgt_key = 'height'
  if hgt_key in vals:
    hgt_key = 'height'
  elif 'heights' in vals:
    hgt_key = 'heights'
  else:
    raise AttributeError('ICON file does not have a valid height label (I know only height, heights, height_2 and heights_2)')
  
  Nh = len(vals[hgt_key])-1
  Nt = len(vals['time'])
  nhydros = 6
  if timeidx is None:
    timeidx = timeidx = np.arange(0,Nt)
  shapeSFC = (len(timeidx),)
    
  date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
  timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
  pamData['timestamp'] = timestamp  
  pamData['hgt_lev'] = np.tile(np.flip(vals[hgt_key],0),(len(timeidx),1)) # heights at which fields are defined

  pamData['press']  = np.flip(vals['P'][timeidx],1)    # pressure 
  pamData['temp']   = np.flip(vals['T'][timeidx],1)    # temperature
  wind_u = np.flip(vals['U'][timeidx],1)               # zonal wind speed
  wind_v = np.flip(vals['V'][timeidx],1)               # meridional wind speed
  pamData['wind_uv'] = np.hypot(wind_u, wind_v)
  wind_w = np.flip(vals['W'][timeidx],1)               # vertical wind speed
  pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
  pamData['relhum'] = np.flip(vals['REL_HUM'][timeidx],1)

  # Read hydrometeors content
  hydro_cmpl = np.zeros((len(timeidx),Nh,nhydros))
  hydro_cmpl[...,0] = hydro_content[0]*np.flip(vals['QC'][timeidx],1)   # specific cloud water content
  hydro_cmpl[...,1] = hydro_content[1]*np.flip(vals['QI'][timeidx],1)   # specific cloud ice content
  hydro_cmpl[...,2] = hydro_content[2]*np.flip(vals['QR'][timeidx],1)   # rain mixing ratio
  hydro_cmpl[...,3] = hydro_content[3]*np.flip(vals['QS'][timeidx],1)   # snow mixing ratio
  hydro_cmpl[...,4] = hydro_content[4]*np.flip(vals['QG'][timeidx],1)   # graupel mixing ratio
  hydro_cmpl[...,5] = hydro_content[5]*np.flip(vals['QH'][timeidx],1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

  # Read hydrometeors number concentration
  hydro_num_cmpl = np.zeros((len(timeidx),Nh,nhydros))
  hydro_num_cmpl[...,0] = hydro_content[0]*np.flip(vals['QNC'][timeidx],1)  # number concentration of cloud water
  hydro_num_cmpl[...,1] = hydro_content[1]*np.flip(vals['QNI'][timeidx],1)  # number concentration ice
  hydro_num_cmpl[...,2] = hydro_content[2]*np.flip(vals['QNR'][timeidx],1)  # number concentration droplets
  hydro_num_cmpl[...,3] = hydro_content[3]*np.flip(vals['QNS'][timeidx],1)  # number concentration snow
  hydro_num_cmpl[...,4] = hydro_content[4]*np.flip(vals['QNG'][timeidx],1)  # number concentration graupel
  hydro_num_cmpl[...,5] = hydro_content[5]*np.flip(vals['QNH'][timeidx],1)  # number concentration hail 

  pamData["hydro_q"] = hydro_cmpl
  pamData["hydro_n"] = hydro_num_cmpl
  
  # surface properties
  pamData['wind10u'] = vals['U10M'][timeidx]
  pamData['wind10v'] = vals['V10M'][timeidx]
  pamData['groundtemp'] = vals['T_S'][timeidx]

# surface properties
  pamData['sfc_type']  = np.ones(pamData['groundtemp'].shape)
  pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
  pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
  pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land

  pam = pyPamtra()
  pam.set['pyVerbose'] = verbosity

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else: # is an iterable containing arrays
    for df in descriptorFile:
      pam.df.addHydrometeor(df)
  
  pam.createProfile(**pamData)
  
  return pam

def readIcon2momOnFlightTrack(fname, descriptorFile, debug=False, verbosity=0):
  '''
  import ICON LEM 2-moment dataset output along a flight track
  
  WARNING: This importer has been designed upon Icon output generated by Dr. Vera Schemann
  Might not work on other files.
  
  fname = filename of the Icon output
  descriptorFile = pyPamtra descriptorFile object or filename of a proper formatted descriptorFile
  debug = flag causing stop and load of debugger upon raised exception
  verbosity = pyPamtra.pyVerbose verbosity level

  Pamtra assumes that the apart from height, the first given dimension is latitude, here coordinates do not change across the dimension, but time does, thus generating a meteogram
  
  16/05/2018 - Mario Mech - mario.mech@uni-koeln.de based on readIcon2momMeteogram by Davide Ori
  '''
  
  import netCDF4
  
  ICON_file = netCDF4.Dataset(fname, mode='r')
  vals = ICON_file.variables

  pamData = dict() # empty dictionary to store pamtra Data
  
  Nh = len(vals['height_2d'])
  Nx = len(vals['lat'])
  nhydros = 6
  # if timeidx is None:
  latidx = np.arange(0,Nx)
  # shapeSFC = (len(timeidx),)
    
  # date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
  # timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
  # pamData['timestamp'] = timestamp
  
#  variables1D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
#  variables2D = ["t_g","pres_sfc"]
#  variables3D_10m = ["u_10m","v_10m"]
#  variables3D = ["temp","pres","qv","qc","qi","qr","qs","qg","qh","qnc","qni","qnr","qns","qng","qnh"]
  
  pamData['hgt'] = np.swapaxes(np.flip(vals['height_2d'],0),0,1) # heights at which fields are defined

  pamData['press'] = np.swapaxes(np.flip(vals['P'],0),0,1)    # pressure 
  pamData['temp'] = np.swapaxes(np.flip(vals['T'],0),0,1)    # temperature
  #pamData['wind_u']  = np.flip(vals['U'][timeidx],1)    # zonal wind speed
  #pamData['wind_v']  = np.flip(vals['V'][timeidx],1)    # meridional wind speed
  # wind_w = np.flip(vals['W'][timeidx],1)    # vertical wind speed
  # pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
  pamData['relhum'] = np.swapaxes(np.flip(vals['REL_HUM'],0),0,1)

  # Read hydrometeors content
  hydro_cmpl = np.zeros((len(latidx),Nh,nhydros))
  hydro_cmpl[...,0] = np.swapaxes(np.flip(vals['QC'],0),0,1)   # specific cloud water content
  hydro_cmpl[...,1] = np.swapaxes(np.flip(vals['QI'],0),0,1)   # specific cloud ice content
  hydro_cmpl[...,2] = np.swapaxes(np.flip(vals['QR'],0),0,1)   # rain mixing ratio
  hydro_cmpl[...,3] = np.swapaxes(np.flip(vals['QS'],0),0,1)   # snow mixing ratio
  hydro_cmpl[...,4] = np.swapaxes(np.flip(vals['QG'],0),0,1)   # graupel mixing ratio
  hydro_cmpl[...,5] = np.swapaxes(np.flip(vals['QH'],0),0,1)   # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

  # Read hydrometeors number concentration
  hydro_num_cmpl = np.zeros((len(latidx),Nh,nhydros))
  hydro_num_cmpl[...,0] = np.swapaxes(np.flip(vals['QNC'],0),0,1)  # number concentration of cloud water
  hydro_num_cmpl[...,1] = np.swapaxes(np.flip(vals['QNI'],0),0,1)  # number concentration ice
  hydro_num_cmpl[...,2] = np.swapaxes(np.flip(vals['QNR'],0),0,1)  # number concentration droplets
  hydro_num_cmpl[...,3] = np.swapaxes(np.flip(vals['QNS'],0),0,1)  # number concentration snow
  hydro_num_cmpl[...,4] = np.swapaxes(np.flip(vals['QNG'],0),0,1)  # number concentration graupel
  hydro_num_cmpl[...,5] = np.swapaxes(np.flip(vals['QNH'],0),0,1)  # number concentration hail 

  pamData["hydro_q"] = hydro_cmpl
  pamData["hydro_n"] = hydro_num_cmpl
  
  #pamData["lfrac"] = np.array([1]) ####
  #pamData["groundtemp"] = pamData['temp'][149]
  #pamData['phalf'] = vals['PHALF'][:]# pressure on the half levels

  pam = pyPamtra()
  pam.set['pyVerbose'] = verbosity
  
  # surface properties
  pamData['lat'] = vals['lat'][:]
  pamData['lon'] = vals['lon'][:]
#  pamData['lon'] = lon
  # pdb.set_trace()  
  pamData['wind10u'] = vals['u10'][:]
  pamData['wind10v'] = vals['v10'][:]
  pamData['groundtemp'] = vals['TS'][:]
#  pamData['sfc_type'] = np.around(frland*np.ones(shapeSFC))
#  assert np.all(np.logical_or(0<=pamData['sfc_type'], pamData['sfc_type'] <=1))
#  pamData['sfc_model'] = np.zeros(shapeSFC)
#  pamData['sfc_refl'] = np.chararray(shapeSFC)
#  pamData['sfc_refl'][pamData['sfc_type'] == 0] = 'F' # ocean
#  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'L' # land

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)
  
  pam.createProfile(**pamData)
  
  return pam
  
  

def createUsStandardProfile(pam=pyPamtra(),**kwargs):
  '''
  Function to create clear sky US Standard Atmosphere.

  hgt_lev is the only mandatory variables
  humidity will be set to zero if not provided, all other variables are guessed by "createProfile"

  values provided in kwargs will be passed to "createProfile", however, press_lev and temp_lev will overwrite us standard if provided

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

  values provided in kwargs will be passed to "createProfile", however, press_lev and temp_lev will overwrite us standard if provided

  '''

  import usStandard #see in tools

  assert "hgt_lev" in list(kwargs.keys()) #hgt_lev is mandatory

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

  for kk in list(kwargs.keys()):
        pamData[kk] = np.array(kwargs[kk])

  if ("relhum_lev" not in list(kwargs.keys())) and ("q" not in list(kwargs.keys())):
    pamData["relhum_lev"] = np.zeros_like(kwargs["hgt_lev"])

  return pamData

#helper function
def ncToDict(ncFilePath,keys='all',joinDimension='time',offsetKeys={},ncLib='netCDF4',tmpDir="/tmp/",skipFiles=[]):
  '''
  Load keys of netcdf file into dictionary that is returned.
  By using wildcards in ncFilePath, multiple files are possible. They are joint using joinDimension. 
  Offsets of e.g.,, time vector can be corrected by offsetKeys with value to be corrected as key as correction key in value.
  gzip compressed netcdf with extension .gz are possible.
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
      print("skipping:", ncFile)

  for nn, ncFile in enumerate(ncFiles):
    tmpFile=False
    if ncFile.split(".")[-1]=="gz":
      tmpFile = True
      gzFile = deepcopy(ncFile)
      ncFile = tmpDir+"/maxLibs_netcdf_"+''.join(randomLib.choice(string.ascii_uppercase + string.digits) for x in range(5))+".tmp.nc"
      print('uncompressing', nn+1,'of',len(ncFiles), gzFile, ncFile)
      returnValue = os.system("zcat "+gzFile+">"+ncFile)
      assert returnValue == 0
    else:
      print('opening', nn+1,'of',len(ncFiles), ncFile)
    try: ncData = openNc(ncFile,'r')
    except: raise RuntimeError("Could not open file: '" + ncFile+"'")
    if nn == 0:
      if keys == 'all':
        keys = list(ncData.variables.keys())
      #make sure the join dimension is actually present!
      if noFiles > 1: assert joinDimension in list(ncData.dimensions.keys())
      #get the axis to join the arrays
      for key in keys:
        joinDimensionNumber[key] = -9999
        for dd, dim in enumerate(ncData.variables[key].dimensions):
          if dim == joinDimension:
            joinDimensionNumber[key] = dd
    #get and join the data
    for key in keys:
      #special sausage for nc files with time offset
      if key in list(offsetKeys.keys()):
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

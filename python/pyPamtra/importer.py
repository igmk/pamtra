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
    stop and load debugger on exception (the default is False)
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
    stop and load debugger on exception (the default is False)
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
This importer is adjusted to read regridded ICON-SRM simulations as they where run by Daniel Klocke within the HErZ-NARVALII framework.
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

def readIconNwp1MomDataset(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0,diagnostic=False):
  """
  import ICON SRM 1-moment dataset

  It is Icon with cosmo physics

  Parameters
  ----------
  fname_fg : {str}
    path and name of file with atmospheric variables.
  descriptorFile : {str}
    path and name of descriptor file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False)
  verbosity : {int}, optional
    verbosity of the module (the default is 0) (the default is 0)
  constantFields : {str}, optional
    path and name of file with constant fields (the default is None)
  maxLevel : {int}, optional
    Number of maximum levels to consider. (the default is 0)
  diagnostic : {bool}
    Switch to use diagnostic qc and qi instead of prognostic ones is set to True (the default is False)

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

  if diagnostic:
    variables4D.remove('qc')
    variables4D.remove('qi')
    variables4D_cloud.append('tot_qc_dia')
    variables4D_cloud.append('tot_qi_dia')

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

  data["hydro_q"] = np.zeros(data["qr"].shape + (nHydro,)) + np.nan
  if diagnostic:
    data["hydro_q"][...,0] = data["tot_qc_dia"]
    data["hydro_q"][...,1] = data["tot_qi_dia"]
  else:
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
readIconNwp1MomDataset.__doc__ += __ICDN_regridding_remarks


def readIconNwp2MomDataset(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0):
  """
  import ICON SRM 2-moment dataset

  It is ICON with cosmo physics

  Parameters
  ----------
  fname_fg : {str}
    path and name of file with atmospheric variables.
  descriptorFile : {str}
    path and name of descriptor file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False)
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
readIconNwp2MomDataset.__doc__ += __ICDN_regridding_remarks

def readIconNwp1MomDataset_cells(fname_fg,descriptorFile,debug=False,verbosity=0,constantFields=None,maxLevel=0):
  """
  import ICON SRM 1-moment dataset where certain cell were sellected

  Such data gan be generated using $(cdo -selgridcell)
  It is ICON with cosmo physics

  Use PAMTRAs internal "lat" index for cell. "lon"-dimension is 1.

  Parameters
  ----------
  fname_fg : {str}
    path and name of file with atmospheric variables.
    fname_fg has to include "_fg_".
    A corresponding "_cloud_" has do exist.
    No wildCards allowed. Must be nc file.
  descriptorFile : {str}
    path and name of descriptor file
  debug : {bool}, optional
    stop and load debugger on exception (the default is False)
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

  variables2D_const = ["FR_LAND", 'lon', 'lat'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
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
      # nc dimensions: cell; target dimensions: lon, lat
      assert ncFile_const.variables[var].dimensions == ('cell', )
      dataSingle[var] = ncFile_const.variables[var][:][np.newaxis, :]

    ncFile_const.close()
    if verbosity > 1: print("closed const nc")

    ncFile_fg = netCDF4.Dataset(fname_fg, "r")
    if verbosity > 1: print("opened ", fname_fg)

    fname_cloud = '_cloud_'.join(fname_fg.rsplit('_fg_',1)) # right-replace _fg_ with _cloud
    ncFile_cloud = netCDF4.Dataset(fname_cloud, "r")
    if verbosity > 1: print("opened ", fname_cloud)

    if maxLevel == 0:
      assert ncFile_fg.variables["z_ifc"].dimensions == ('height_3', 'cell')
      assert ncFile_fg.dimensions['height'].size + 1 == ncFile_fg.dimensions['height_3'].size
      maxLevel = ncFile_fg.variables["z_ifc"].shape[0] - 1

    for var in variables3D:
      # nc dimensions: time, cell; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'cell')
      dataSingle[var] = ncFile_fg.variables[var][:][forecastIndex, np.newaxis, :]

    for var in variables4D_10m:
      # nc dimensions: time, height_5, cell; target dimensions: lon, lat
      assert ncFile_fg.variables[var].dimensions == ('time', 'height_5', 'cell')
      assert ncFile_fg.dimensions['height_5'].size == 1
      dataSingle[var] = ncFile_fg.variables[var][:][forecastIndex, np.newaxis, 0, :]

    for var in variables4D:
      # nc dimensions: time, height, cell; target dimensions: lon, lat, level
      assert ncFile_fg.variables[var].dimensions == ('time', 'height', 'cell')
      dataSingle[var] = np.transpose(ncFile_fg.variables[var][forecastIndex], (1, 0))[np.newaxis, :, ::-1][...,:maxLevel]#reverse height order

    for var in variables4D_cloud:
      # nc dimensions: time, height, cell; target dimensions: lon, lat, level
      assert ncFile_cloud.variables[var].dimensions == ('time', 'height', 'cell')
      dataSingle[var] = np.transpose(ncFile_cloud.variables[var][forecastIndex], (1, 0))[np.newaxis, :, ::-1][...,:maxLevel]#reverse height order

    shape3D = dataSingle["temp"].shape
    assert shape3D[0] == 1, 'Lon dimension should not be used'
    shape3Dplus = (shape3D[0], shape3D[1], shape3D[2] + 1)
    shape2D = shape3D[:2]

    date_times = netCDF4.num2date(ncFile_fg.variables["time"][:], ncFile_fg.variables["time"].units) # convert from any given reference time
    timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp

    dataSingle["timestamp"] = np.zeros(shape2D)
    dataSingle["timestamp"][:] = timestamp

    for key in ["z_ifc"]:
      # nc dimensions: height, cell; target dimensions: lon, lat, level
      assert ncFile_fg.variables[key].dimensions == ('height_3', 'cell')
      dataSingle[key] = np.transpose(ncFile_fg[key], (1, 0))[np.newaxis, : ,::-1][...,:maxLevel+1]#reverse height order
      assert dataSingle[key].shape == shape3Dplus

    ncFile_fg.close()
    if verbosity > 1: print("closed fg nc")
    ncFile_cloud.close()
    if verbosity > 1: print("closed cloud nc")

    data = dataSingle

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

  varPairs = [["timestamp","timestamp"],["lat","lat"],["lon","lon"],
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
  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'L' # land

  pam.createProfile(**pamData)

  return pam

# Add regridding remarks
readIconNwp1MomDataset_cells.__doc__ += __ICDN_regridding_remarks

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
    stop and load debugger on exception (the default is False)
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
    stop and load debugger on exception (the default is False)
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
    stop and load debugger on exception (the default is False)
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
  import ICON LEM 1-moment dataset output cross section over a site (Meteogram)

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
  pamData['wind_w'] = -0.5*(wind_w[:,:-1]+wind_w[:,1:]) #sign change due to conventions: for model upward is positive; in pamtra downward velocity is positive 
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
  pamData['wind_w'] = -0.5*(wind_w[:,:-1]+wind_w[:,1:]) #sign change due to conventions: for model upward is positive; in pamtra downward velocity is positive 

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

def readIcon2momOnFlightTrack(fname, descriptorFile, time=0, kind='processed', constant_file=None, debug=False, verbosity=0):
  '''
  import ICON LEM 2-moment dataset output along a flight track
  
  WARNING: This importer has been originally designed upon Icon output processed by Dr. Vera Schemann
  It has been later on extended to read extracted flight tracks without further processing.
  Might not work on other files.

  fname = filename of the Icon output
  descriptorFile = pyPamtra descriptorFile object or filename of a proper formatted descriptorFile
  processed = flag to specify whether the files have been processed or not. (default = True)
  debug = flag causing stop and load of debugger upon raised exception
  verbosity = pyPamtra.pyVerbose verbosity level

  Pamtra assumes that the apart from height, the first given dimension is latitude, here coordinates do not change across the dimension, but time does, thus generating a meteogram

  16/05/2018 - Mario Mech - mario.mech@uni-koeln.de based on readIcon2momMeteogram by Davide Ori
  '''

  import netCDF4

  ICON_file = netCDF4.Dataset(fname, mode='r')
  vals = ICON_file.variables

  pamData = dict() # empty dictionary to store pamtra Data

  if kind == 'processed':
  
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

  elif kind == 'exact':


    ICONC_file = netCDF4.Dataset(constant_file, mode='r')
    cvals = ICONC_file.variables

    time = 0
    Nh = len(cvals['height'])
    Nx = len(cvals['clon'])
    nhydros = 6
      
    # date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
    # timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
    # pamData['timestamp'] = timestamp
    
  #  variables1D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
    variables2D = ['t_g','fr_land','fr_seaice']
    variables3D_10m = ['u_10m','v_10m']
    variables3D = ['hgt','temp','pres','qv']
    variables3Dplus = ['hgt_lev']
    variables3Dqx = ['qc','qi','qr','qs','qg','qh']
    variables3Dqnx = ['qnc','qni','qnr','qns','qng','qnh']
    
    # pamData['hgt'] = np.swapaxes(np.flip(vals['height_2d'],0),0,1) # heights at which fields are defined
    # pamData['press'] = np.swapaxes(np.flip(vals['P'],0),0,1)    # pressure 
    # pamData['temp'] = np.swapaxes(np.flip(vals['T'],0),0,1)    # temperature
    pamData['hgt'] = np.swapaxes(cvals['z_mc'][::-1,:],0,1)
    pamData['hgt_lev'] = np.swapaxes(cvals['z_ifc'][::-1,:],0,1)
    pamData['press'] = np.swapaxes(vals['pres'][time,::-1,:],0,1)
    pamData['temp'] = np.swapaxes(vals['temp'][time,::-1,:],0,1)
    pamData['relhum'] = np.swapaxes(q2rh(vals['qv'][time,::-1,:],vals['temp'][time,::-1,:],vals['pres'][time,::-1,:]) * 100.,0,1)
    #pamData['wind_u']  = np.flip(vals['U'][timeidx],1)    # zonal wind speed
    #pamData['wind_v']  = np.flip(vals['V'][timeidx],1)    # meridional wind speed
    # wind_w = np.flip(vals['W'][timeidx],1)    # vertical wind speed
    # pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
    # pamData['relhum'] = np.swapaxes(np.flip(vals['REL_HUM'],0),0,1)

    # Read hydrometeors content
    pamData['hydro_q'] = np.zeros((Nx,Nh,nhydros))

    for i,qi in enumerate(variables3Dqx):
      pamData['hydro_q'][...,i] = np.swapaxes(vals[qi][time,::-1,...],0,1)

    # pamData['hydro_q'][...,0] = np.swapaxes(vals['qc'][time,...],0,1)  # specific cloud water content
    # pamData['hydro_q'][...,1] = np.swapaxes(vals['qi'][time,...],0,1)  # specific cloud ice content
    # pamData['hydro_q'][...,2] = np.swapaxes(vals['qr'][time,...],0,1) # rain mixing ratio
    # pamData['hydro_q'][...,3] = np.swapaxes(vals['qs'][time,...],0,1) # snow mixing ratio
    # pamData['hydro_q'][...,4] = np.swapaxes(vals['qg'][time,...],0,1)  # graupel mixing ratio
    # pamData['hydro_q'][...,5] = np.swapaxes(vals['qh'][time,...],0,1) # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

    # Read hydrometeors number concentration
    pamData['hydro_n'] = np.zeros((Nx,Nh,nhydros))

    for i,ni in enumerate(variables3Dqnx):
      pamData['hydro_n'][...,i] = np.swapaxes(vals[ni][time,::-1,...],0,1)
    # pamData['hydro_n'][...,0] = np.swapaxes(vals['qnc'][time,...],0,1)  # number concentration of cloud water
    # pamData['hydro_n'][...,1] = np.swapaxes(vals['qni'][time,...],0,1)  # number concentration ice
    # pamData['hydro_n'][...,2] = np.swapaxes(vals['qnr'][time,...],0,1)  # number concentration droplets
    # pamData['hydro_n'][...,3] = np.swapaxes(vals['qns'][time,...],0,1)  # number concentration snow
    # pamData['hydro_n'][...,4] = np.swapaxes(vals['qng'][time,...],0,1)  # number concentration graupel
    # pamData['hydro_n'][...,5] = np.swapaxes(vals['qnh'][time,...],0,1)  # number concentration hail 
    
    # surface properties
    pamData['lat'] = np.rad2deg(cvals['clat'][:])
    pamData['lon'] = np.rad2deg(cvals['clon'][:])
    pamData['wind10u'] = vals['u_10m'][time,0,:]
    pamData['wind10v'] = vals['v_10m'][time,0,:]
    pamData['groundtemp'] = vals['t_g'][time,:]

    pamData['sfc_slf'] = cvals['fr_land'][:]
    pamData['sfc_sif'] = vals['fr_seaice'][0,...]

  # surface properties
    pamData['sfc_type']  = np.around(pamData['sfc_slf'])
    pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
    pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
    pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land
    pamData['sfc_type'][(pamData['sfc_type'] == 0) & (pamData['sfc_sif'] > 0)] = 2

    ICONC_file.close()

  elif kind == 'original':

    Nh = len(vals['height'])
    Nx = len(vals['clon'])
    Ny = len(vals['clat'])
    nhydros = 6
      
    # date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
    # timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
    # pamData['timestamp'] = timestamp
    
    variables2D = ['t_g','fr_land','fr_seaice']
    variables3D_10m = ['u_10m','v_10m']
    variables3D = ['hgt','temp','pres','qv']
    variables3Dplus = ['hgt_lev']
    variables3Dqx = ['qc','qi','qr','qs','qg','qh']
    variables3Dqnx = ['qnc','qni','qnr','qns','qng','qnh']
    
    pamData['hgt'] = np.swapaxes(vals['z_mc'][::-1,:],0,1)
    pamData['hgt_lev'] = np.swapaxes(vals['z_ifc'][::-1,:],0,1)
    pamData['press'] = np.swapaxes(vals['pres'][time,::-1,:],0,1)
    pamData['temp'] = np.swapaxes(vals['temp'][time,::-1,:],0,1)
    pamData['relhum'] = np.swapaxes(q2rh(vals['qv'][time,::-1,:],vals['temp'][time,::-1,:],vals['pres'][time,::-1,:]) * 100.,0,1)

    # Read hydrometeors content
    pamData['hydro_q'] = np.zeros((Nx,Nh,nhydros))

    for i,qi in enumerate(variables3Dqx):
      pamData['hydro_q'][...,i] = np.swapaxes(vals[qi][time,::-1,...],0,1)

    # Read hydrometeors number concentration
    pamData['hydro_n'] = np.zeros((Nx,Nh,nhydros))

    for i,ni in enumerate(variables3Dqnx):
      pamData['hydro_n'][...,i] = np.swapaxes(vals[ni][time,::-1,...],0,1)
    
    # surface properties
    pamData['lat'] = np.rad2deg(vals['clat'][:])
    pamData['lon'] = np.rad2deg(vals['clon'][:])
    pamData['wind10u'] = vals['u_10m'][time,0,:]
    pamData['wind10v'] = vals['v_10m'][time,0,:]
    pamData['groundtemp'] = vals['t_g'][time,:]

    pamData['sfc_slf'] = vals['fr_land']
    pamData['sfc_sif'] = vals['fr_seaice'][0,...]

  # surface properties
    pamData['sfc_type']  = np.around(pamData['sfc_slf'])
    pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
    pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
    pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land
    pamData['sfc_type'][(pamData['sfc_type'] == 0) & (pamData['sfc_sif'] > 0)] = 2

  else:

    time = 0
    Nh = len(vals['height'])
    Nx = len(vals['clat'])
    nhydros = 6
      
    # date_times = netCDF4.num2date(vals["time"][timeidx], vals["time"].units) # datetime autoconversion
    # timestamp = netCDF4.date2num(date_times, "seconds since 1970-01-01 00:00:00") # to unix epoch timestamp (python uses nanoseconds int64 internally)
    # pamData['timestamp'] = timestamp
    
  #  variables1D_const = ["FR_LAND", 'lon_2', 'lat_2'] # 'topography_c' would be the ICON equivalent to the COSMO 'HSURF'
    variables2D = ['t_g','fr_land','fr_seaice']
    variables3D_10m = ['u_10m','v_10m']
    variables3D = ['hgt','temp','pres','qv']
    variables3Dplus = ['hgt_lev']
    variables3Dqx = ['qc','qi','qr','qs','qg','qh']
    variables3Dqnx = ['qnc','qni','qnr','qns','qng','qnh']
    
    # pamData['hgt'] = np.swapaxes(np.flip(vals['height_2d'],0),0,1) # heights at which fields are defined
    # pamData['press'] = np.swapaxes(np.flip(vals['P'],0),0,1)    # pressure 
    # pamData['temp'] = np.swapaxes(np.flip(vals['T'],0),0,1)    # temperature
    pamData['hgt'] = np.swapaxes(vals['z_mc'][::-1,:],0,1)
    pamData['hgt_lev'] = np.swapaxes(vals['z_ifc'][::-1,:],0,1)
    pamData['press'] = np.swapaxes(vals['pres'][time,::-1,:],0,1)
    pamData['temp'] = np.swapaxes(vals['temp'][time,::-1,:],0,1)
    pamData['relhum'] = np.swapaxes(q2rh(vals['qv'][time,::-1,:],vals['temp'][time,::-1,:],vals['pres'][time,::-1,:]) * 100.,0,1)
    #pamData['wind_u']  = np.flip(vals['U'][timeidx],1)    # zonal wind speed
    #pamData['wind_v']  = np.flip(vals['V'][timeidx],1)    # meridional wind speed
    # wind_w = np.flip(vals['W'][timeidx],1)    # vertical wind speed
    # pamData['wind_w'] = 0.5*(wind_w[:,:-1]+wind_w[:,1:])
    # pamData['relhum'] = np.swapaxes(np.flip(vals['REL_HUM'],0),0,1)

    # Read hydrometeors content
    pamData['hydro_q'] = np.zeros((Nx,Nh,nhydros))

    for i,qi in enumerate(variables3Dqx):
      pamData['hydro_q'][...,i] = np.swapaxes(vals[qi][time,::-1,...],0,1)

    # pamData['hydro_q'][...,0] = np.swapaxes(vals['qc'][time,...],0,1)  # specific cloud water content
    # pamData['hydro_q'][...,1] = np.swapaxes(vals['qi'][time,...],0,1)  # specific cloud ice content
    # pamData['hydro_q'][...,2] = np.swapaxes(vals['qr'][time,...],0,1) # rain mixing ratio
    # pamData['hydro_q'][...,3] = np.swapaxes(vals['qs'][time,...],0,1) # snow mixing ratio
    # pamData['hydro_q'][...,4] = np.swapaxes(vals['qg'][time,...],0,1)  # graupel mixing ratio
    # pamData['hydro_q'][...,5] = np.swapaxes(vals['qh'][time,...],0,1) # graupel mixing ratio # TODO report probably error encoding long name ...  should be hail mixing ratio

    # Read hydrometeors number concentration
    pamData['hydro_n'] = np.zeros((Nx,Nh,nhydros))

    for i,ni in enumerate(variables3Dqnx):
      pamData['hydro_n'][...,i] = np.swapaxes(vals[ni][time,::-1,...],0,1)
    # pamData['hydro_n'][...,0] = np.swapaxes(vals['qnc'][time,...],0,1)  # number concentration of cloud water
    # pamData['hydro_n'][...,1] = np.swapaxes(vals['qni'][time,...],0,1)  # number concentration ice
    # pamData['hydro_n'][...,2] = np.swapaxes(vals['qnr'][time,...],0,1)  # number concentration droplets
    # pamData['hydro_n'][...,3] = np.swapaxes(vals['qns'][time,...],0,1)  # number concentration snow
    # pamData['hydro_n'][...,4] = np.swapaxes(vals['qng'][time,...],0,1)  # number concentration graupel
    # pamData['hydro_n'][...,5] = np.swapaxes(vals['qnh'][time,...],0,1)  # number concentration hail 
    
    # surface properties
    pamData['lat'] = np.rad2deg(vals['clat'][:])
    pamData['lon'] = np.rad2deg(vals['clon'][:])
    pamData['wind10u'] = vals['u_10m'][time,0,:]
    pamData['wind10v'] = vals['v_10m'][time,0,:]
    pamData['groundtemp'] = vals['t_g'][time,:]

    pamData['sfc_slf'] = vals['fr_land']
    pamData['sfc_sif'] = vals['fr_seaice'][0,...]

  # surface properties
    pamData['sfc_type']  = np.around(pamData['sfc_slf'])
    pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
    pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
    pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land
    pamData['sfc_type'][(pamData['sfc_type'] == 0) & (pamData['sfc_sif'] > 0)] = 2

  pam = pyPamtra()
  pam.set['pyVerbose'] = verbosity

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else:
    for df in descriptorFile:
      pam.df.addHydrometeor(df)

  pam.createProfile(**pamData)
  
  ICON_file.close()

  return pam
  
def readICON1momNWP(fname, descriptorFile, debug=False, verbosity=0,grid=[0,463,0,6653]):
  '''
  import ICON LEM 1-moment dataset output cross section over a site (Meteogram)
  
  WARNING: This importer has been designed upon ICON 1 mom NWP output generated by Dr. Helene Bresson.
  Will most likely not work on other files.
  
  fname = filename of the ICON output
  descriptorFile = pyPamtra descriptorFile object or filename of a proper formatted descriptorFile
  debug = flag causing stop and load of debugger upon raised exception
  verbosity = pyPamtra.pyVerbose verbosity level
  
  22/03/2020 - Mario Mech - mario.mech@uni-koeln.de
  '''

  dataTmp = ncToDict(fname)  
  
  # variables2D = ['fr_land','t_g','u_10m','v_10m']
  # variables3D = ['z_mc','temp','pres','qv','qc','qi','qr','qs','qg']
  # variables3Dplus = ['z_ifc']

  variables = ['fr_land','fr_seaice','t_g','u_10m','v_10m','temp','pres','qv','qc','qi','qr','qs','qg','z_mc','z_ifc']

  data = dict()
  nhydro = 5

  a = grid[0]  
  b = grid[1]  
  c = grid[2]  
  d = grid[3]
  nlat = b-a
  nlon = d-c

  data['lat'] = dataTmp['lat'][a:b]
  data['lon'] = dataTmp['lon'][c:d]
  for key in variables:
    data[key] = dataTmp[key][...,a:b,c:d]

  data['time'] = dataTmp['time']

  del dataTmp
  # construct all shapes needed for PAMTRA
  shape2D = data['fr_land'].shape
  shape3D = (data['z_mc'].shape[1],data['z_mc'].shape[2],data['z_mc'].shape[0])
  shape3Dplus = (shape3D[0],shape3D[1],shape3D[2]+1)
  shape4D = (shape3D[0],shape3D[1],shape3D[2],nhydro)

  pamData = dict() # empty dictionary to store pamtra Data

  pamData['timestamp'] = np.zeros(shape2D) + data['time'][0]
  pamData['lat'] = np.reshape(np.repeat(data['lat'],nlon),(nlat,nlon))
  pamData['lon'] = np.reshape(np.repeat(data['lon'],nlat),(nlon,nlat)).T

  pamData['hgt_lev'] = np.moveaxis(data['z_ifc'][::-1,:,:],0,-1) # heights at which fields are defined
  pamData['press']  = np.moveaxis(data['pres'][0,::-1,:,:],0,-1)  # pressure 
  pamData['temp']   = np.moveaxis(data['temp'][0,::-1,:,:],0,-1)    # temperature
  pamData['relhum'] = np.moveaxis(q2rh(data['qv'][0,::-1,:,:],data['temp'][0,::-1,:,:],data['pres'][0,::-1,:,:]) * 100.,0,-1)

  # Read hydrometeors
  pamData["hydro_q"] = np.zeros(shape4D)

  pamData["hydro_q"][...,0] = np.moveaxis(data['qc'][0,::-1,:,:],0,-1)
  pamData["hydro_q"][...,1] = np.moveaxis(data['qi'][0,::-1,:,:],0,-1)
  pamData["hydro_q"][...,2] = np.moveaxis(data['qr'][0,::-1,:,:],0,-1)
  pamData["hydro_q"][...,3] = np.moveaxis(data['qs'][0,::-1,:,:],0,-1)
  pamData["hydro_q"][...,4] = np.moveaxis(data['qg'][0,::-1,:,:],0,-1)
  
  # surface properties
  pamData['wind10u'] = data['u_10m'][0,0,...]             # zonal wind speed
  pamData['wind10v'] = data['v_10m'][0,0,...]               # meridional wind speed
  pamData['groundtemp'] = data['t_g'][0,...]
  pamData['sfc_slf'] = data['fr_land']
  pamData['sfc_sif'] = data['fr_seaice'][0,...]

# surface properties
  pamData['sfc_type']  = np.around(pamData['sfc_slf'])
  pamData['sfc_model'] = np.zeros(pamData['groundtemp'].shape)
  pamData['sfc_refl']  = np.chararray(pamData['groundtemp'].shape)
  pamData['sfc_refl'][:] = 'S' # land  'F' # ocean 'L' lambertian, land
  pamData['sfc_type'][(pamData['sfc_type'] == 0) & (pamData['sfc_sif'] > 0)] = 2

  pam = pyPamtra()
  pam.set['pyVerbose'] = verbosity

  if isinstance(descriptorFile, str):
    pam.df.readFile(descriptorFile)
  else: # is an iterable containing arrays
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

def read_AVAPS_dropsonde(file_raw_sondes, sst=np.nan, opt_T=np.nan, opt_P=np.nan,
  opt_RH=np.nan, opt_Z=np.nan, verbose=False):

  """
    PAMTRA importer for raw AVAPS dropsonde data (preferred ending: 'PQC.nc') in clear sky
    conditions over oceans.

    Parameters
    ----------
    file_raw_sondes : str
      File of raw dropsonde data (preferred ending: 'PQC.nc').
    sst : float, optional
      Sea surface temperature (in K). If this parameter is not set, the lowest dropsonde
      temperature measurement is used (eventually extrapolated).
    opt_T : float, optional
      Optional temperature (e.g. measured by BAHAMAS) (in K) at flight altitude.
    opt_P: float, optional
      Optional pressure (e.g. measured by BAHAMAS) (in Pa) at flight altitude.
    opt_RH: float, optional
      Optional relative humidity (e.g. measured by BAHAMAS) (in %) at flight altitude.
    opt_Z: float, optional
      Optional flight altitude during sonde launch (e.g. taken from BAHAMAS) (in m).
    verbose : bool, optional
      If True some extra information is printed. May be extensive if many dropsondes are imported!
  """

  # FUNCTIONS

  def fill_gaps(old_var, verbose=False):
    # old variable gets linearly interpolated for each sonde launch. The function is ignoring nan 
    # values at the surface and above the launch altitude.

    new_var = deepcopy(old_var)
    # create flag variable indicating if an entry of old_var has been changed: if = 0: not interpol.
    interp_flag = np.zeros(old_var.shape)

    # identify regions of nan values in the middle of the drop. Extrapolation will be handled in another 
    # function.
    # identify the highest non-nan entry so we can cut the values above that highest entry:
    # identify the lowest non-nan entry for similar reasons:
    non_nan_idx = [idx for idx, x in np.ndenumerate(old_var) if (not np.isnan(x))]
    non_nan_idx = np.where(~np.isnan(old_var))[0]
    limits = np.array([non_nan_idx[0], non_nan_idx[-1]])

    temp_var = deepcopy(old_var)
    temp_var = temp_var[limits[0]:limits[1]+1]		# will be the variable where the gaps are filled

    interp_flag_temp = np.zeros(temp_var.shape)

    # identify mid-drop-nan-values: need values after and before the nan:
    nan_idx = np.argwhere(np.isnan(temp_var))
    interp_flag_temp[nan_idx] = 1

    if nan_idx.size == 0:
      return new_var, interp_flag

    else: # correct nan values: find the hole size via subtraction of subsequent indices
      hole_size = np.zeros((len(nan_idx)+1,)).astype(int)
      k = 0		# index to address a hole ('hole number')

      for m in range(0, len(temp_var)-1):

        if not np.isnan(temp_var[m+1] - temp_var[m]):
          hole_size[k] = 0

        elif np.isnan(temp_var[m+1] - temp_var[m]): # k only incremented if an END of a hole has been identified:
          if len(nan_idx) == 1: 	# handled seperately in case that merely one nan value exists in temp_var
            hole_size[k] = 1
            break

          else:
            if (not np.isnan(temp_var[m+1])) & (np.isnan(temp_var[m])): # END of a hole
              k = k + 1
              continue
            hole_size[k] = hole_size[k] + 1		# k won't be incremented until m finds another non-nan value

        else:
          raise RuntimeError("Something unexpected happened when trying to find nan values in the middle " +
            "of the dropsonde launch. Please contact 'a.walbroel@uni-koeln.de'. \n")

      # holes have been identified: edit the FIRST hole (editing depends on the size of the hole...)
      c = 0 		# dummy variable needed for the right jumps in hole_size and nan_idx. c is used to address
                      # nan_idx and therefore new_var...

      # meanwhile 'd' just runs through the array hole_size:
      for d in range(0, len(hole_size)):
        for L in range(0, hole_size[d]):		# range(0, 1): L = 0
          temp_var[nan_idx[c] + L] = temp_var[nan_idx[c] - 1] + (L + 1)*(temp_var[int(nan_idx[c] +
            hole_size[d])] - temp_var[nan_idx[c]-1]) / (hole_size[d] + 1)

        c = c + int(hole_size[d])
        if c > len(hole_size)-1:
          break

    # overwrite the possibly holey section:
    new_var[limits[0]:limits[1]+1] = temp_var
    # update interp_flag
    interp_flag[limits[0]:limits[1]+1] = interp_flag_temp

    return new_var, interp_flag

  def std_extrapol_BAH(old_dict, ill_keys, bah_dict, old_ipflag_dict=dict(), verbose=False):
  # Will extrapolate some atmospheric variables to the ceiling of the dropsonde; old_ipflag will be updated.
  # Needs the old variable, the interpolation flag (should've been generated by fill_gaps()), the key and 
  # height levels as INPUT

    new_dict = old_dict
    n_alt = len(new_dict['Z'])

    new_ipflag_dict = old_ipflag_dict

    # To get the obs_height: Floor BAHAMAS altitude to the next lowest 100 m.
    drop_alt = np.floor(np.asarray(bah_dict['Z'])/100)*100
    obs_height = np.max(np.unique(drop_alt))		# omit repeated values; this value is used for the top of the extrapolation

    # BAHAMAS (or other optional measurement platform at flight altitude) temperature, pressure, relative humidity at launch time:
    bah_T = bah_dict['T']
    bah_P = bah_dict['P']
    bah_RH = bah_dict['RH']
    bah_alt = bah_dict['Z']

    ceiling = obs_height	# last entry of altitude
    if ceiling > 15000:
      raise ValueError("Dropsonde launch altitude appears to be > 15000 m. Extrapolation is aborted because the " +
        "tropopause may intervene.\n")
      return new_dict, new_ipflag_dict

    # Any value above obs_height will be deleted: So if e.g. Z has got values above obs_height, delete them:
    # Find the first index that overshoots obs_height:
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      overshoot = np.argwhere(new_dict['Z'] >= obs_height)
    if len(overshoot) > 0:
      overshoot = overshoot[0][0] + 1

      for key in new_dict.keys():
        if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
          continue

        if new_dict[key].ndim > 0:		# otherwise: error when using len()

          if len(new_dict[key]) == n_alt:
            new_dict[key] = new_dict[key][:overshoot]	# limit variable to obs_height
            if key in ill_keys or key == 'Z':
              new_ipflag_dict[key] = new_ipflag_dict[key][:overshoot]

    # At the end of 'Z' there may still be nans -> so we don't know to which altitude meteorological variables belong to in this region:
    # Therefore: delete it and replace by extrapolation:
    n_alt = len(new_dict['Z'])
    last_nonnan_alt = np.argwhere(~np.isnan(new_dict['Z']))[-1][0]
    if ceiling - new_dict['Z'][last_nonnan_alt] > 1000:
      if verbose:
        print("WARNING: highest GPS altitude measurement is at least 1000 m below the aircraft. Extrapolation may be erroneous.\n")

    for key in new_dict.keys():
      if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
        continue

      if new_dict[key].ndim > 0:		# otherwise: error when using len()

        if len(new_dict[key]) == n_alt:
          new_dict[key] = new_dict[key][:last_nonnan_alt+1]	# limit variable to obs_height
          if key in ill_keys or key == 'Z':
            new_ipflag_dict[key] = new_ipflag_dict[key][:last_nonnan_alt+1]

    # Extend the old height grid up to the ceiling if the distance is greater than 10 meters:
    alt = new_dict['Z']
    n_alt = len(alt)
    alt = np.append(alt, np.arange(alt[np.argwhere(~np.isnan(alt))[-1]]+10, ceiling+11, 10))
    n_alt_new = len(alt)

    # Update the altitude variable in the dictionary: & update ipflag for gpsalt:
    new_dict['Z'] = alt
    new_ipflag_dict['Z'] = np.append(new_ipflag_dict['Z'], np.ones((n_alt_new - n_alt,)))

    launch_time = datetime.datetime.utcfromtimestamp(new_dict['launch_time']).strftime("%Y-%m-%d %H:%M:%S") # for printing

    for key in ill_keys:

      new_var = new_dict[key]
      # Must be expanded to the new height grid (up to the new ceiling)
      new_var = np.append(new_var, np.full((n_alt_new - n_alt,), np.nan), axis=0) # append nans at the top of the profile

      if not new_ipflag_dict: # in case fill_gaps(...) hasn't been called before this one, it's assumed that nothing has been interpolated yet.
        new_ipflag_dict[key] = np.zeros(new_var.shape)

      else: # new_ipflag also has to be extended to the new hgt grid:
        new_ipflag_dict[key] = np.append(new_ipflag_dict[key], np.zeros((n_alt_new - n_alt,)), axis=0)


      if key == 'T': 
        # If BAHAMAS Temperature measurement is available use it as target in case only the top 15 % of measurements
        # are missing. Otherwise:
		# Temperature: If dropsondes with measurements (ipflag = 0) from that day exist, estimate their average T gradient.
		# If the extrapolated dropsonde temperature then deviates from BAH T by more than 5 K, use the ICAO std atmosphere
		# as T gradient:
        # Standard atmosphere (shifted accordingly to avoid a jump between the last known value and the extrapolation).
        # ICAO standard atmosphere taken from:
        # https://www.dwd.de/DE/service/lexikon/begriffe/S/Standardatmosphaere_pdf.pdf?__blob=publicationFile&v=3
        ICAO_standard_T = 288.15 - 0.0065*alt

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]


        if alt[idx] < 0.6*ceiling:
          if verbose:
            print("Insufficient amount of measurements for temperature extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no temperature measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var			# then just overwrite the dictionary entry with the nonedited (but extended) variable
          continue

        if alt[idx] < 0.85*ceiling: # then use BAHAMAS temperature as extrapolation target:
          new_var[idx+1:] = (new_var[idx] + (bah_T - new_var[idx]) / (bah_alt - alt[idx]) * (alt[idx+1:] - alt[idx]))

        else:
          # Or use mean T gradient of highest 20 measurements and continue with this gradient:
          # Compute mean T gradient of highest 20 measurements:
          mean_T_grad = np.mean(np.asarray([(new_var[idx-19:idx+1] - new_var[idx-20:idx]) / (alt[idx-19:idx+1] - alt[idx-20:idx])]))
          extra_T = 288.15 + mean_T_grad*alt

          new_var[idx+1:] = extra_T[idx+1:] - (extra_T[idx] - new_var[idx])

          if np.abs(new_var[-1] - bah_T) > 5:    # then the deviation from BAHAMAS T is too great and we use the ICAO std atmosphere
            new_var[idx+1:] = ICAO_standard_T[idx+1:] - (ICAO_standard_T[idx] - new_var[idx])


        new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'P':
        # Pressure: use hydrostatic eq. with scale height H = R <T> / g0, R = 287 J kg^-1 K^-1, g0 = 9.81 m s^-2, using 
        # the vertican mean temperature <T>:
        # p(z) = p0 exp( -z / H) (Holton, p.21);

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]


        if alt[idx] < ceiling/3:
          if verbose:
            print("Insufficient amount of measurements for pressure extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no pressure measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        # MAKE SURE THAT mean TEMPERATURE CAPTURES THE ACTUAL MEAN TEMPERATURE!!
        if np.count_nonzero(~np.isnan(new_dict['T'][:])) / float(n_alt_new) <= 0.75:  # in this case you may expect that a mean temperature would yield
                                                                                      # a bad representation of the true scale height. 0.75 was chosen arbitrarily.
          T_icao_0 = 12*np.cos(4*np.pi*np.nanmean(new_dict['lat'][:])/360) + 288.15	# strongly simplified meridional surface temperature structure
          H = 287 * np.mean(288.15 - 0.0065*alt) / 9.81		# using the ICAO standard atmosphere to compute the mean temperature

          if verbose:
            print("Warning: Because insufficient non-nan temperature values were given for launch " +
              launch_time + ", '" + str(np.mean(288.15 - 0.0065*alt)) +
              " K' was assumed to be the mean temperature for hydrostatic pressure calculation. \n")

        else:
          H = 287 * np.nanmean(new_dict['T'][:]) / 9.81		# scale height

          # Find index of lowest non nan pressure measurement:
          l_idx = np.argwhere(~np.isnan(new_var[:]))[0]
          p_ref = new_var[l_idx]		# in Pa
          alt_ref = alt[l_idx]
          p_hydrostat = p_ref * np.exp(-(alt - alt_ref) / H)		# in Pa

          new_var[idx+1:] = p_hydrostat[idx+1:] - (p_hydrostat[idx] - new_var[idx])
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'u_wind' or key == 'v_wind':
        # Wind: idea: Fill nan values with the mean wind gradient of the highest 20 (non-nan)measurents. 
        # It will only be extrapolated if the the last non-nan entry is higher than 0.80*ceiling:
        # Other idea: just keep the same wind value

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]


        if alt[idx] < 0.8*ceiling:
          if verbose:
            print("Insufficient amount of measurements for wind extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no wind measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        else:

          # Alternative: just use the latest value for higher altitudes:
          new_var[idx+1:] = new_var[idx]
          new_var[idx+1:] = new_var[idx+1:]
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'RH':
        # Relative humidity (RH): Idea: Linearly interpolate to the BAHAMAS value

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]


        if alt[idx] < 0.65*ceiling:
          if verbose:
            print("Insufficient amount of measurements for relative humidity extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no rel. hum. measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        else:
          new_var[idx+1:] = (new_var[idx] + (bah_RH - new_var[idx]) / (bah_alt - alt[idx]) * (alt[idx+1:] - alt[idx]))
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag
          new_var[np.argwhere(new_var[:] < 0)] = 0.0

      new_dict[key] = new_var

    return new_dict, new_ipflag_dict, obs_height

  def std_extrapol(old_dict, ill_keys, old_ipflag_dict=dict(), verbose=False):
  # Will extrapolate some atmospheric variables to the ceiling of the dropsonde; old_ipflag will be updated.
  # Needs the old variable, the interpolation flag (should've been generated by fill_gaps()), 
  # the key and height levels as INPUT

    new_dict = old_dict
    n_alt = len(new_dict['Z'])

    new_ipflag_dict = old_ipflag_dict

    # To get the obs_height: find highest ... (?)
    if np.isnan(new_dict['reference_alt']):
      # Select the highest non nan index of T or P. 
      highest_nonnan_Z = np.argwhere(~np.isnan(new_dict['Z']))[-1]
      obs_height = np.floor(new_dict['Z'][highest_nonnan_Z[0]]/100)*100
    else:	# use the reference_alt
      obs_height = (np.floor(new_dict['reference_alt']/100)*100)[0]
		

    ceiling = obs_height	# last entry of altitude
    if ceiling > 15000:
      raise ValueError("Dropsonde launch altitude appears to be > 15000 m. Extrapolation is aborted because the " +
        "tropopause may intervene.\n")
      return new_dict, new_ipflag_dict

    # Any value above obs_height will be deleted: So if e.g. Z has got values above obs_height, delete them:
    # Find the first index that overshoots obs_height:
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      overshoot = np.argwhere(new_dict['Z'] >= obs_height)
    if len(overshoot) > 0:
      overshoot = overshoot[0][0] + 1

      for key in new_dict.keys():
        if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
          continue

        if new_dict[key].ndim > 0:		# otherwise: error when using len()

          if len(new_dict[key]) == n_alt:
            new_dict[key] = new_dict[key][:overshoot]	# limit variable to obs_height
            if key in ill_keys or key == 'Z':
              new_ipflag_dict[key] = new_ipflag_dict[key][:overshoot]

    # At the end of 'Z' there may still be nans -> so we don't know to which altitude meteorological variables belong to in this region:
    # Therefore: delete it and replace by extrapolation:
    n_alt = len(new_dict['Z'])
    last_nonnan_alt = np.argwhere(~np.isnan(new_dict['Z']))[-1][0]
    if ceiling - new_dict['Z'][last_nonnan_alt] > 1000:
      if verbose:
        print("WARNING: highest GPS altitude measurement is at least 1000 m below the aircraft. Extrapolation may be erroneous.\n")

    for key in new_dict.keys():
      if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
        continue

      if new_dict[key].ndim > 0:		# otherwise: error when using len()

        if len(new_dict[key]) == n_alt:
          new_dict[key] = new_dict[key][:last_nonnan_alt+1]	# limit variable to obs_height
          if key in ill_keys or key == 'Z':
            new_ipflag_dict[key] = new_ipflag_dict[key][:last_nonnan_alt+1]


    # Extend the old height grid up to the ceiling if the distance is greater than 10 meters:
    alt = new_dict['Z']
    n_alt = len(alt)
    alt = np.append(alt, np.arange(alt[np.argwhere(~np.isnan(alt))[-1]]+10, ceiling+11, 10))
    n_alt_new = len(alt)

    # Update the altitude variable in the dictionary: & update ipflag for gpsalt:
    new_dict['Z'] = alt
    new_ipflag_dict['Z'] = np.append(new_ipflag_dict['Z'], np.ones((n_alt_new - n_alt,)))


    launch_time = datetime.datetime.utcfromtimestamp(new_dict['launch_time']).strftime("%Y-%m-%d %H:%M:%S") # for printing

    for key in ill_keys:

      new_var = new_dict[key]
      # Must be expanded to the new height grid (up to the new ceiling)
      new_var = np.append(new_var, np.full((n_alt_new - n_alt,), np.nan), axis=0) # append nans at the top of the profile

      if not new_ipflag_dict: # in case fill_gaps(...) wasn't called before this one, it's assumed that nothing has been interpolated yet.
        new_ipflag_dict[key] = np.zeros(new_var.shape)

      else: # new_ipflag also has to be extended to the new hgt grid:
        new_ipflag_dict[key] = np.append(new_ipflag_dict[key], np.zeros((n_alt_new - n_alt,)), axis=0)


      if key == 'T':
        # Temperature: If dropsondes with measurements (ipflag = 0) from that day exist, estimate their average T gradient.
        # In case more than 15 % of T measurements must be extrapolated, the ICAO standard atmosphere is used as T gradient.
        # ICAO standard atmosphere taken from:
        # https://www.dwd.de/DE/service/lexikon/begriffe/S/Standardatmosphaere_pdf.pdf?__blob=publicationFile&v=3
        ICAO_standard_T = 288.15 - 0.0065*alt

        # find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]

        if alt[idx] < 0.6*ceiling:
          if verbose:
            print("Insufficient amount of measurements for temperature extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no temperature measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var			# then just overwrite the dictionary entry with the nonedited (but extended) variable
          continue

        if alt[idx] < 0.85*ceiling: # then use standard atmosphere (ICAO):
          new_var[idx+1:] = ICAO_standard_T[idx+1:] + (new_var[idx] - ICAO_standard_T[idx])

        else:
          # Use mean T gradient of highest 20 measurements and continue with this gradient:
          # Compute mean T gradient of highest 20 measurements:
          mean_T_grad = np.mean(np.asarray([(new_var[idx-19:idx+1] - new_var[idx-20:idx]) / (alt[idx-19:idx+1] - alt[idx-20:idx])]))
          T_continued = 288.15 + mean_T_grad*alt

          new_var[idx+1:] = T_continued[idx+1:] - (T_continued[idx] - new_var[idx])


        new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'P':
        # Pressure: use hydrostatic eq. with scale height H = R <T> / g0, R = 287 J kg^-1 K^-1, g0 = 9.81 m s^-2, using the vertican mean temperature <T>:
        # p(z) = p0 exp( -z / H) (Holton, p.21);

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]

        if alt[idx] < ceiling/3:
          if verbose:
            print("Insufficient amount of measurements for pressure extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no pressure measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        # MAKE SURE THAT mean TEMPERATURE CAPTURES THE ACTUAL MEAN TEMPERATURE!!
        if np.count_nonzero(~np.isnan(new_dict['T'][:])) / float(n_alt_new) <= 0.75: 		# in this case you may expect that a mean temperature would yield
																								# a bad representation of the true scale height. 0.75 was chosen arbitrarily.
          T_icao_0 = 12*np.cos(4*np.pi*np.nanmean(new_dict['lat'][:])/360) + 288.15	# strongly simplified meridional surface temperature structure
          H = 287 * np.mean(288.15 - 0.0065*alt) / 9.81		# using the ICAO standard atmosphere to compute the mean temperature

          if verbose:
            print("Warning: Because insufficient non-nan temperature values were given for launch " +
              launch_time + ", '" + str(np.mean(288.15 - 0.0065*alt)) +
              " K' was assumed to be the mean temperature for hydrostatic pressure calculation. Can possibly be " +
              "avoided if the temperature is extrapolated before the pressure.\n")

        else:
          H = 287 * np.nanmean(new_dict['T'][:]) / 9.81		# scale height

          # Find index of lowest non nan pressure measurement:
          l_idx = np.argwhere(~np.isnan(new_var[:]))[0]
          p_ref = new_var[l_idx]		# in Pa
          alt_ref = alt[l_idx]
          p_hydrostat = p_ref * np.exp(-(alt - alt_ref) / H)		# in Pa

          new_var[idx+1:] = p_hydrostat[idx+1:] - (p_hydrostat[idx] - new_var[idx])
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'u_wind' or key == 'v_wind':
        # Wind: idea: Fill nan values with the mean wind gradient of the highest 20 (non-nan)measurents. It will only be extrapolated if the the last non-nan entry
        # is higher than 0.80*ceiling:
        # Other idea: Just keep the same wind value

        # Find highest non nan value if it lies below the ceiling:
        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]

        if alt[idx] < 0.8*ceiling:
          if verbose:
            print("Insufficient amount of measurements for wind extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no wind measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        else:

          # Alternative: Just use the latest value for higher altitudes:
          new_var[idx+1:] = new_var[idx]
          new_var[idx+1:] = new_var[idx+1:]
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag


      elif key == 'RH':
        # Relative humidity (RH): Idea: fill nan values with the mean RH of the highest 10 measurements but only 
        # if the highest measurement exceeds or is equal to 0.95*ceiling!
        # Other idea: Use RH approx 0 for greater altitudes.
        noise_strength = 0/2		# percent

        idx = np.argwhere(~np.isnan(new_var)).flatten()[-1]

        if alt[idx] < 0.65*ceiling:
          if verbose:
            print("Insufficient amount of measurements for relative humidity extrapolation at the top of the dropsonde grid (" + launch_time +
              "). There are no rel. hum. measurements above " + str(alt[idx]) + " m.\n")
          new_dict[key] = new_var
          continue

        elif alt[idx] >= 0.95*ceiling:
          new_var[idx+1:] = np.mean(new_var[idx-9:idx+1])

        else:
          new_var[idx+1:] = 1.5
          new_ipflag_dict[key][idx+1:] = 1		# setting the interpol flag
          new_var[np.argwhere(new_var[:] < 0)] = 0.0

      new_dict[key] = new_var

    return new_dict, new_ipflag_dict, obs_height

  def regridding(new_dict, obs_height, ill_keys, resolution=10):
    '''
      Regridding variables specified in ill_keys to a uniform grid 
      (from the surface up to obs_height) with a user-defined 
      resolution (in meters, default=10).
    '''

    new_alt = np.arange(0, obs_height+1, 10)

    for key in ill_keys:
      new_dict[key] = np.interp(new_alt, new_dict['Z'], new_dict[key])

    new_dict['Z'] = new_alt

    return new_dict

  def repair_surface(old_dict, ill_keys, old_ipflag_dict=dict(), verbose=False):
  # Filling nan values at the surface if the gap to the surface isn't too large (e.g. measurements below 200 m
  # must exist (roughly 15-20 seconds before splash).
    new_dict = old_dict
    alt = old_dict['Z']
    n_alt = len(alt)
    new_ipflag_dict = old_ipflag_dict
    launch_time = datetime.datetime.utcfromtimestamp(new_dict['launch_time']).strftime("%Y-%m-%d %H:%M:%S")

    lim = 200		# if there are no measurements below this altitude then the extrapolation at the surface won't be performed

    if ill_keys == ['Z']:
      threshold_list = [ill_keys, [200], ['m']]
    else:
      threshold_list = [ill_keys, [5.0, 4000.0, 50.0, 0.1, 0.1, 5.0, 5.0, 1.0],
        ['K', 'hPa', '%', 'deg', 'deg', 'm/s', 'm/s', 'm/s']]	# used to check if surface value deviates siginificantly from lowest measurement

    for key in ill_keys:

      new_var = new_dict[key]

      if not new_ipflag_dict: # in case fill_gaps(...) wasn't called before this one, it's assumed that nothing has been interpolated yet.
        new_ipflag_dict[key] = np.zeros(new_var.shape)


      # find the first non-nan entry
      idx = np.argwhere(~np.isnan(new_var[:]))[0][0]

      if alt[idx] < lim:
        sfc_gap = np.arange(0,idx)

        if len(sfc_gap) == 0:
          continue
        else:

          # create mean gradient of the variable of 10 measurements above the lowest measurement, or, if grid is too coarse, take 200-400 m average:
          if alt[idx+10] > 2*lim: # default: if alt[idx+10] > 400: # take lowest measurement to 400 m mean gradient
            # find index closest to 400 m:
            idx2 = np.argmin(np.abs(alt - 2*lim))

          else: # take mean grad. of 10 measurem. above lowest measurement:
            idx2 = idx+10

          mean_grad = np.mean([new_var[j+1] - new_var[j] for j in range(idx,idx2)])	# theoretically, should never be nan because fill_gaps
																					# should've fixed the holes between the first and last measurement
          for j in sfc_gap:
            new_var[idx-j-1] = new_var[idx] - mean_grad*(j+1)


          # check if sfc value not too far off the lowest measurement:
          if key == 'RH':
            if np.any(new_var[sfc_gap] < 0):
              new_var[sfc_gap] = 0
              if verbose:
                print("Caution, '" + key + "' surface repair resulted in negative values. Manually setting the missing values at the ground to 0 for launch "
                  + launch_time + ".\n")

            elif np.any(new_var[sfc_gap] > 100):
              new_var[sfc_gap] = 100
              if verbose:
                print("Caution, '" + key + "' surface repair resulted in >100 %. Manually setting the missing values at the ground to 100 for launch "
                  + launch_time + ".\n")

          threshold = threshold_list[1][threshold_list[0].index(key)]
          si_unit = threshold_list[2][threshold_list[0].index(key)]
          if np.abs(new_var[0] - new_var[idx]) > threshold:
            if verbose:
              print("Caution, '" + key + "' surface value deviates more than " + str(threshold) + " " + si_unit + " from the lowest measurement (launch "
                + launch_time + ").\n")


          new_ipflag_dict[key][sfc_gap] = 1

      else:
        if verbose:
         print("No measurements below " + str(lim) + " m. Extrapolation of '" + key + "', launch " + launch_time +
           " would eventually lead to wrong assumptions at the surface. Therefore aborted.\n")
        continue

    return new_dict, new_ipflag_dict

  def mark_outliers(sonde_dict, ill_keys): # mark outliers: outliers defined when exceeding certain thresholds

    new_dict = sonde_dict

    # Thresholds are defined by change of meteorol. variable with altitude: e.g. delta p / delta z
    thresholds = [0.065, 40, 2.5, 1, 1]	# in [K/m, Pa/m, %/m, ms^-1/m, ms^-1/m]
    dz = new_dict['Z'][1:] - new_dict['Z'][:-1]	# delta z

    for key in ill_keys:
      if key == 'lat' or key == 'lon':
        continue

      met_threshold = thresholds[ill_keys.index(key)]		# threshold for key

      d_met = new_dict[key][1:] - new_dict[key][:-1]		# change of meteorological variable 'key'

      with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        exceed_idx = np.argwhere(np.abs(d_met / dz) >= met_threshold)
      new_dict[key][exceed_idx] = float('nan')

    return new_dict

  def readDropsondeNCraw(filename, verbose=False):
    """
      Dropsonde netCDF file to dictionary.
      Variables with a designated unit that is not an SI unit will be converted to SI units. 
      Rel. humidity will be in % though.
      Time will be given in seconds since 1970-01-01 00:00:00 UTC.

      Parameters
      ----------
      filename : str
        Filename of the raw dropsonde files (suffix "PQC.nc"!)
      verbose : bool
        Additional comments printed if True.
    """

    import netCDF4 as nc

    file_nc = nc.Dataset(filename)

    dropsonde_dict = dict()
    dropsonde_dict['fillValues'] = dict()
    for nc_keys in file_nc.variables.keys():
      nc_var = file_nc.variables[nc_keys]
      dropsonde_dict[nc_keys] = np.asarray(nc_var)

      if hasattr(nc_var, 'missing_value'):
        dropsonde_dict['fillValues'][nc_keys] = nc_var.missing_value

        # type of the nc_var:
        ncvar_type = type(dropsonde_dict[nc_keys][0])

        # find where the current variable has missing values and set them to nan:
        missing_idx = np.argwhere(dropsonde_dict[nc_keys] == dropsonde_dict['fillValues'][nc_keys])

        if ((ncvar_type == np.float32) or (ncvar_type == np.float64)):
          dropsonde_dict[nc_keys][missing_idx] = float('nan')

      # converting units: time stamps will be handled seperately.
      if hasattr(nc_var, 'units'):
        if nc_var.units == 'degC':
          dropsonde_dict[nc_keys] = dropsonde_dict[nc_keys] + 273.15
          if verbose: print("From degC to K: " + str(nc_keys))
        elif nc_var.units == 'hPa':
          dropsonde_dict[nc_keys] = dropsonde_dict[nc_keys]*100
          if verbose: print("From hPa to Pa: " + str(nc_keys))
        elif nc_var.units == 'gram/kg':
          dropsonde_dict[nc_keys] = dropsonde_dict[nc_keys]/1000
          if verbose: print("From g/kg to kg/kg: " + str(nc_keys))

    time_base = datetime.datetime.strptime(file_nc.variables['launch_time'].units[14:-4], "%Y-%m-%d %H:%M:%S") # time base given in the units attribute
    dropsonde_dict['launch_time'] = (time_base - datetime.datetime(2020,1,1)).total_seconds() + dropsonde_dict['launch_time']
    dropsonde_dict['reference_time'] = (time_base - datetime.datetime(2020,1,1)).total_seconds() + dropsonde_dict['reference_time']
    dropsonde_dict['time'] = (time_base - datetime.datetime(2020,1,1)).total_seconds() + dropsonde_dict['time']
    if verbose: print("\n")

    # Convert to internal convention: Temperature = T, Pressure = P, relative humidity = RH, altitude = Z
    dropsonde_dict['T'] = dropsonde_dict['tdry']
    dropsonde_dict['P'] = dropsonde_dict['pres']
    dropsonde_dict['RH'] = dropsonde_dict['rh']
    dropsonde_dict['Z'] = dropsonde_dict['gpsalt']

    return dropsonde_dict



  # Dropsonde quality control: if there are small gaps of measurements between launch altitude and surface, they will be filled via interpolation.
  # Additionally, the sondes will be extrapolated to a certain altitude.

  # BAHAMAS data (optional):
  if np.all([~np.isnan(opt_T), ~np.isnan(opt_P), ~np.isnan(opt_RH), ~np.isnan(opt_Z)]):
    opt_dict = {'T': opt_T, 'P': opt_P, 'RH': opt_RH, 'Z': opt_Z}
  elif np.all([np.isnan(opt_T), np.isnan(opt_P), np.isnan(opt_RH), np.isnan(opt_Z)]):
    opt_dict = dict()
  else:
    raise ValueError("Optional measurements at flight altitude must include all four variables: Temperature (in K), pressure (in Pa), relative " +
      "humidity (in %), flight altitude (in m). If not all measurements can be provided, leave optional data input empty.")

  failed_sondes = []				# lists the filenames of failed sondes (nearly no measurements)
  stuck_sondes = []					# lists the filenames of stuck sondes (sondes that don't show values above 1500 m)

  sonde_dict = readDropsondeNCraw(file_raw_sondes, verbose)
  launch_date = datetime.datetime.utcfromtimestamp(sonde_dict['launch_time']).strftime("%Y-%m-%d")
  dropsonde_date = (datetime.datetime.strptime(launch_date, "%Y-%m-%d")).strftime("%Y%m%d")	# date displayed in the filename
    # --> comfy way to find the right BAHAMAS data for std_extrapol
  if verbose:
    print("########## Launch time: " + datetime.datetime.utcfromtimestamp(sonde_dict['launch_time']).strftime("%Y-%m-%d %H:%M:%S") + " ##########\n")
    print("Input: ", file_raw_sondes)

  # Add another condition that checks if e.g. nearly no measurements exist at all (for T, P and RH):
  if np.any([np.count_nonzero(~np.isnan(sonde_dict['T'])) < 0.1*len(sonde_dict['T']),
    np.count_nonzero(~np.isnan(sonde_dict['P'])) < 0.1*len(sonde_dict['P']),
    np.count_nonzero(~np.isnan(sonde_dict['RH'])) < 0.1*len(sonde_dict['RH']),
    np.count_nonzero(~np.isnan(sonde_dict['lat'])) < 0.05*len(sonde_dict['lat']),
    np.count_nonzero(~np.isnan(sonde_dict['lon'])) < 0.05*len(sonde_dict['lon']),
    np.count_nonzero(~np.isnan(sonde_dict['u_wind'])) < 0.05*len(sonde_dict['u_wind']),
    np.count_nonzero(~np.isnan(sonde_dict['v_wind'])) < 0.05*len(sonde_dict['v_wind'])]):
    failed_sondes.append(file_raw_sondes)
    raise ValueError("One PAMTRA-critical variable measurement failed. Skipping this dropsonde. \n")

  # Add yet another condition that checks if the sonde got stuck mid air:
  if not np.any(sonde_dict['Z'][~np.isnan(sonde_dict['Z'])] < 1500):	# then I assume that the whole launch was doomed
    raise ValueError("'gpsalt' doesn't seem to include any values below 1500 m. That is insufficient for extrapolation at the surface. \n")
    stuck_sondes.append(file_raw_sondes)

  # Subsequent variables will be cured from holey nan value disease...:
  # pressure, temperature, relhum, wind (u & v & w), lat, lon.
  # Other variables, of which you can expect a linear interpolated over gaps to be applicable, may be added.
  ill_keys = ['T', 'P', 'RH', 'lat', 'lon', 'u_wind', 'v_wind']

  sonde_ipflag = dict()		# will contain the interpolation flags for interpolated nan values in the middle of the drop

  # It's known that the altitude axis (gpsalt) does have some broken values (e.g. sudden jumps or exceeding the aircraft altitude).
  # Before filliing the gaps, the altitude axis must be fixed:
  sonde_dict['Z'], sonde_ipflag['Z'] = fill_gaps(sonde_dict['Z'], verbose=verbose)
  for key in ill_keys:
    sonde_dict[key], sonde_ipflag[key] = fill_gaps(sonde_dict[key], verbose=verbose)	# altitude must be passed to check for dimensions of the to-be-cured variable...

  sonde_dict['ipflag'] = sonde_ipflag

  # The raw dropsonde files show an altitude variable with increases from [0 to -1] in general: but probably due to gps tracking,
  # the altitude decreases at the "top": this must be filtered out:
  # Find highest index where altitude[idx+1] - altitude[idx] > 0:
  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    altitude_stop = np.argwhere(sonde_dict['Z'][1:] - sonde_dict['Z'][:-1] > 0)[-1][0] + 2
    # +2 because it's used as indexing [... : altitude_stop] => +1 and because the array size had been reduced by 1 during argwhere(...)

  # If the lowest non nan value of the altitude coordinate is < 0: cut the rest off:
  # Find lowest non nan value of altitude:
  lowest = np.argwhere(~np.isnan(sonde_dict['Z'][:]))[0][0]

  ndata = len(sonde_dict['Z'])

  # Then cut each variable at this altitude index:
  if sonde_dict['Z'][lowest] < 0:
    for key in sonde_dict.keys():
      if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
        continue

      if sonde_dict[key].ndim > 0:		# otherwise: error when using len()

        if len(sonde_dict[key]) == ndata:
          sonde_dict[key] = sonde_dict[key][lowest:altitude_stop]
          if key in ill_keys or key == 'Z':
            sonde_dict['ipflag'][key] = sonde_dict['ipflag'][key][lowest:altitude_stop]

  else:
    for key in sonde_dict.keys():
      if key in ['trajectory', 'fillValues', 'ipflag']:	# skip these ones ... it s not interesting anyway
        continue

      if sonde_dict[key].ndim > 0:		# otherwise: error when using len()

        if len(sonde_dict[key]) == ndata:
          sonde_dict[key] = sonde_dict[key][:altitude_stop]
          if key in ill_keys or key == 'Z':
            sonde_dict['ipflag'][key] = sonde_dict['ipflag'][key][:altitude_stop]

  # The altitude index may be a bit broken... needs to be fixed. mark them as nan and let it run through fill_gaps again:
  dz = sonde_dict['Z'][1:] - sonde_dict['Z'][:-1]		# dz[i] = z[i+1] - z[i]

  for k in range(len(dz)):
    if (dz[k] < 0) or (dz[k] > 15):		# this filters too big jumps in the altitude coordinate
      sonde_dict['Z'][k+1] = float('nan')

  sonde_dict['Z'], sonde_ipflag['Z'] = fill_gaps(sonde_dict['Z'], verbose=verbose)

  # Perform surface repair for altitude coordinate:
  sonde_dict, sonde_dict['ipflag'] = repair_surface(sonde_dict, ['Z'], sonde_dict['ipflag'], verbose=verbose)
  # for some reasons, 'alt' is sort of unused but an assigned key in the dictionary. So ... we can also just set it to gpsalt:
  sonde_dict['alt'] = sonde_dict['Z']

  # now we still need to handle the nan values at the surface: If there are no non-nan values in the lowest 5 % of the variable
  # -> don't interpolate because the assumption would eventually lead to senseless surface values:
  sonde_dict, sonde_dict['ipflag'] = repair_surface(sonde_dict, ill_keys, sonde_dict['ipflag'], verbose=verbose)

  # Extrapolating the ill_keys to the ceiling of the dropsondes (e.g. below aircraft altitude):
  # CAUTION: it is expected that the dropsondes start BELOW THE TROPOPAUSE!
  # It is advisable to do the TEMPERATURE EXTRAPOLATION FIRST because it improves the pressure extrapolation.
  if opt_dict:
    sonde_dict, sonde_dict['ipflag'], obs_height = std_extrapol_BAH(sonde_dict, ill_keys, opt_dict, sonde_dict['ipflag'], verbose=verbose)

  else:
    sonde_dict, sonde_dict['ipflag'], obs_height = std_extrapol(sonde_dict, ill_keys, sonde_dict['ipflag'], verbose=verbose)

  # Regridding to a uniform vertical grid with a user-specified resolution:
  sonde_dict = regridding(sonde_dict, obs_height, ill_keys, 10)

  # Find outliers and mark them (as nan): afterwards fill them again
  sonde_dict = mark_outliers(sonde_dict, ['T', 'P', 'RH', 'u_wind', 'v_wind'])
  for key in ill_keys:
    sonde_dict[key], sonde_ipflag[key] = fill_gaps(sonde_dict[key], verbose=verbose)


  # Check for remaining NaNs in meteorological variables:
  # RH, T and P (and the height coordinate Z):
  if (np.any([np.isnan(sonde_dict['T']), np.isnan(sonde_dict['P']), np.isnan(sonde_dict['RH']),
    np.isnan(sonde_dict['Z'])])):
    if verbose:
      print('     NaN-counts: %d, %d, %d, %d' % (
        np.isnan(sonde_dict['T']).sum(),
        np.isnan(sonde_dict['P']).sum(),
        np.isnan(sonde_dict['RH']).sum(),
        np.isnan(sonde_dict['Z']).sum()))
    raise ValueError("Nans still remained in T, P, RH or Z. Aborting execution of dropsonde %s. Please contact 'a.walbroel@uni-koeln.de'."
      % file_raw_sondes)

  # Check if surface winds are not nan: (index [1] is used because [0] is at 0 m above sea level)
  sfc_wind_available = ~np.isnan(sonde_dict['u_wind'][1] + sonde_dict['v_wind'][1])
  if not sfc_wind_available:
    if verbose:
      print("WARNING, NaN in u_wind[1] or v_wind[1]. Surface winds won't be included in PAMTRA for %s !"%(
        file_raw_sondes))
      print('     NaN-counts (u,v): %d, %d' % (
        np.isnan(sonde_dict['u_wind'][1]),
        np.isnan(sonde_dict['v_wind'][1])))


  # Create pamtra object; change settings:
  pam = pyPamtra()

  pam.nmlSet['passive'] = True						# passive simulation
  pam.nmlSet['active'] = False						# False: no radar simulation

  # Create dictionary for PAMTRA data:
  pamData = dict()

  Nh = len(sonde_dict['Z'])		# number of height levels
  Nx = 1						# number of sondes (is always 1 because one sonde is considered at a time)
  nhydros = 0					# number of hydrometeors (here: clear sky assumed)
  shape2d = [Nx, Nx]
  shape3d = shape2d + [Nh]

  pamData['hgt_lev'] = np.broadcast_to(sonde_dict['Z'], shape3d)	  # height
  pamData['press_lev'] = np.broadcast_to(sonde_dict['P'][:], shape3d)	  # pressure
  pamData['temp_lev'] = np.broadcast_to(sonde_dict['T'][:], shape3d)	  # temperature
  pamData['relhum_lev'] = np.broadcast_to(sonde_dict['RH'][:], shape3d)  # relative humidity

  if sfc_wind_available:
    pamData['wind10u'] = np.broadcast_to(sonde_dict['u_wind'][1], shape2d)
    pamData['wind10v'] = np.broadcast_to(sonde_dict['v_wind'][1], shape2d)

  if not np.any([np.isnan(sonde_dict['reference_lon']), np.isnan(sonde_dict['reference_lat'])]):
    pamData['lon'] = np.broadcast_to(sonde_dict['reference_lon'], shape2d)
    pamData['lat'] = np.broadcast_to(sonde_dict['reference_lat'], shape2d)

  else:		# use highest non nan lat/lon values:
    pamData['lon'] = np.broadcast_to(sonde_dict['lon'][~np.isnan(sonde_dict['lon'])][-1], shape2d)
    pamData['lat'] = np.broadcast_to(sonde_dict['lat'][~np.isnan(sonde_dict['lat'])][-1], shape2d)

  pamData['timestamp'] = np.broadcast_to(sonde_dict['launch_time'], shape2d)
  if np.isnan(sst):
    pamData['groundtemp'] = np.broadcast_to(sonde_dict['T'][0], shape2d)
  else:
    pamData['groundtemp'] = np.broadcast_to(sst, shape2d)

  # Obseravtion height:
  obs_height = np.asarray(obs_height).flatten()
  pamData['obs_height'] = np.broadcast_to(obs_height, shape2d + [len(obs_height), ])

  # Surface type & reflectivity:
  pamData['sfc_type'] = np.zeros(shape2d)			# 0: ocean, 1: land
  pamData['sfc_refl'] = np.chararray(shape2d)
  pamData['sfc_refl'][:] = 'F'
  pamData['sfc_refl'][pamData['sfc_type'] == 1] = 'S'

  # 4d variables: hydrometeors:
  shape4d = [Nx, Nx, Nh-1, 1]			# potentially 5 hydrometeor classes with this setting
  pamData['hydro_q'] = np.zeros(shape4d)
  pamData['hydro_q'][...,0] = 0 # CLOUD

  # Descriptorfile must be included. otherwise, pam.p.nhydro would be 0 which is not permitted.
  descriptorFile = np.array([
    #['hydro_name' 'as_ratio' 'liq_ice' 'rho_ms' 'a_ms' 'b_ms' 'alpha_as' 'beta_as' 'moment_in' 'nbin' 'dist_name' 'p_1' 'p_2' 'p_3' 'p_4' 'd_1' 'd_2' 'scat_name' 'vel_size_mod' 'canting']
    ('cwc_q', -99.0, 1, -99.0, -99.0, -99.0, -99.0, -99.0, 3, 1, 'mono', -99.0, -99.0, -99.0, -99.0, 2e-05, -99.0, 'mie-sphere', 'khvorostyanov01_drops', -99.0)],
    dtype=[('hydro_name', 'S15'), ('as_ratio', '<f8'), ('liq_ice', '<i8'), ('rho_ms', '<f8'), ('a_ms', '<f8'), ('b_ms', '<f8'), ('alpha_as', '<f8'), ('beta_as', '<f8'),
    ('moment_in', '<i8'), ('nbin', '<i8'), ('dist_name', 'S15'), ('p_1', '<f8'), ('p_2', '<f8'), ('p_3', '<f8'),
    ('p_4', '<f8'), ('d_1', '<f8'), ('d_2', '<f8'), ('scat_name', 'S15'), ('vel_size_mod', 'S30'), ('canting', '<f8')]
    )
  for hyd in descriptorFile: pam.df.addHydrometeor(hyd)

  pam.createProfile(**pamData)

  return pam

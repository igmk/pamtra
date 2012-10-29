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

	
def readCosmoDe1MomDataset(fnames,kind,forecastIndex = 1,colIndex=0,tmpDir="/tmp/",fnameInTar="",concatenateAxis=1,debug=False,verbosity=0):
	import netCDF4
	
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
	
	if kind not in ["collum"]:
		raise TypeError("unknown Cosmo DE data type")

	variables1Dx = ["time1h"]
	variables1Dy = ["fr_land","latitude","longitude",]
	variables2D = [ "hfl","hhl"]
	variables3D = ["u_10m","v_10m","t_2m","surface_air_pressure"]
	variables4D = ["temperature","p","qv","qc","qi","qi","qr","qs","qg"]
	
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
	  
	data["p_hl"] = np.zeros(shape3Dplus)
	data["t_hl"] = np.zeros(shape3Dplus)
	
	data["t_hl"][...,0] = data["t_2m"]
	data["t_hl"][...,1:-1] = data["temperature"][...,0:-1]+0.5*np.diff(data["temperature"],axis=-1)
	data["t_hl"][...,-1] = data["temperature"][...,-1] + 0.25*(data["temperature"][...,-1]-data["temperature"][...,-2])
	
	data["p_hl"][...,0] = data["surface_air_pressure"]

	press1 = deepcopy(data["p"])
	press1[press1==missingNumber]=1
	
	p0 = press1[...,0:-1]
	p1 = press1[...,1:]	
	dz = np.diff(data["hfl"],axis=-1)
	if neAvail: xp = ne.evaluate("-1.*log(p1/p0)/dz")
	else: xp = -1.*np.log(p1/p0)/dz
			
	xp[xp==0] = 9999
		
	if neAvail: data["p_hl"][...,1:-1] = ne.evaluate("-1.*p0/xp*(exp(-xp*dz)-1.)/dz")
	else: data["p_hl"][...,1:-1] = -1.*p0/xp*(np.exp(-xp*dz)-1.)/dz
	
	data["p_hl"][...,-1] = -1.*p1[...,-1]/xp[...,-1]*(np.exp(-xp[...,-1]*dz[...,-1])-1.)/dz[...,-1]  #check!
	print "TO DO: check calculation of p half level!"
	del data["p"], data["temperature"], p0, p1, dz, xp, press1

	#translate variables
	varPairs = [["time1h","timestamp"],["latitude","lat"],["longitude","lon"],["fr_land","lfrac"],["u_10m","wind10u"],["v_10m","wind10v"],["hhl","hgt_lev"],["p_hl","press_lev"],["t_hl","temp_lev"],["qv","q"],["qc","cwc_q"],["qi","iwc_q"],["qr","rwc_q"],["qs","swc_q"],["qg","gwc_q"]]
	
	pamData = dict()
	for cosmoVar,pamVar in varPairs:
		pamData[pamVar] = data[cosmoVar]

	pam = pyPamtra.pyPamtra()
	pam.set["pyVerbose"]= verbosity
	pam.createProfile(**pamData)
	del data
	
	return pam
	
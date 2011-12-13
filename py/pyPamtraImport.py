# -*- coding: utf-8 -*-

#todo: function to calculate PIA?
from __future__ import division
import numpy as np
import datetime

import pyPamtra
import meteoSI

#import csv
#import pickle
#import time,calendar,datetime
#import warnings
#import sys
#import os
#from copy import deepcopy

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
